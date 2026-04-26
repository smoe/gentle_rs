//! Virtual protein-gel model and rendering primitives.
//!
//! This module keeps the protein molecular-weight gel path separate from the
//! DNA pool-gel renderer so workflow examples can render SDS-PAGE-like figures
//! from first-class protein derivation reports without overloading bp-based
//! migration semantics.

use std::sync::LazyLock;

use svg::Document;
use svg::node::element::{Circle, Line, Rectangle, Text};

const SVG_WIDTH: f32 = 1900.0;
const SVG_HEIGHT: f32 = 840.0;
const GEL_LEFT: f32 = 90.0;
const GEL_RIGHT: f32 = SVG_WIDTH - 520.0;
const GEL_TOP: f32 = 90.0;
const GEL_BOTTOM: f32 = SVG_HEIGHT - 110.0;
const DETAIL_PANEL_LEFT: f32 = GEL_RIGHT + 72.0;
const DETAIL_PANEL_TOP: f32 = GEL_TOP + 8.0;
const DETAIL_PANEL_WIDTH: f32 = SVG_WIDTH - DETAIL_PANEL_LEFT - 36.0;

#[derive(Clone, Debug)]
pub struct ProteinGelSample {
    pub name: String,
    pub detail: Option<String>,
    pub molecular_weight_kda: f32,
}

#[derive(Clone, Debug)]
pub struct ProteinGelBand {
    pub kda: f32,
    pub intensity: f32,
    pub label: String,
}

#[derive(Clone, Debug)]
pub struct ProteinGelLane {
    pub name: String,
    pub detail: Option<String>,
    pub is_ladder: bool,
    pub bands: Vec<ProteinGelBand>,
}

#[derive(Clone, Debug)]
pub struct ProteinGelLayout {
    pub lanes: Vec<ProteinGelLane>,
    pub selected_ladders: Vec<String>,
    pub sample_count: usize,
    pub protein_count: usize,
    pub range_min_kda: f32,
    pub range_max_kda: f32,
    pub notes: Vec<String>,
}

#[derive(Clone, Copy, Debug)]
struct ProteinLadderBand {
    kda: f32,
    relative_strength: f32,
}

#[derive(Clone, Debug)]
struct ProteinLadder {
    name: &'static str,
    bands: &'static [ProteinLadderBand],
}

#[derive(Clone, Debug)]
struct ProteinLadderCatalog {
    ladders: Vec<ProteinLadder>,
}

static PROTEIN_LADDERS: LazyLock<ProteinLadderCatalog> =
    LazyLock::new(ProteinLadderCatalog::default);

impl ProteinLadderCatalog {
    fn default() -> Self {
        Self {
            ladders: vec![
                ProteinLadder {
                    name: "Protein Ladder 10-250 kDa",
                    bands: &[
                        ProteinLadderBand {
                            kda: 10.0,
                            relative_strength: 0.40,
                        },
                        ProteinLadderBand {
                            kda: 15.0,
                            relative_strength: 0.55,
                        },
                        ProteinLadderBand {
                            kda: 25.0,
                            relative_strength: 0.70,
                        },
                        ProteinLadderBand {
                            kda: 37.0,
                            relative_strength: 0.82,
                        },
                        ProteinLadderBand {
                            kda: 50.0,
                            relative_strength: 0.95,
                        },
                        ProteinLadderBand {
                            kda: 75.0,
                            relative_strength: 1.00,
                        },
                        ProteinLadderBand {
                            kda: 100.0,
                            relative_strength: 0.95,
                        },
                        ProteinLadderBand {
                            kda: 150.0,
                            relative_strength: 0.82,
                        },
                        ProteinLadderBand {
                            kda: 250.0,
                            relative_strength: 0.62,
                        },
                    ],
                },
                ProteinLadder {
                    name: "Protein Ladder 10-100 kDa",
                    bands: &[
                        ProteinLadderBand {
                            kda: 10.0,
                            relative_strength: 0.46,
                        },
                        ProteinLadderBand {
                            kda: 15.0,
                            relative_strength: 0.62,
                        },
                        ProteinLadderBand {
                            kda: 20.0,
                            relative_strength: 0.72,
                        },
                        ProteinLadderBand {
                            kda: 25.0,
                            relative_strength: 0.82,
                        },
                        ProteinLadderBand {
                            kda: 37.0,
                            relative_strength: 0.94,
                        },
                        ProteinLadderBand {
                            kda: 50.0,
                            relative_strength: 1.00,
                        },
                        ProteinLadderBand {
                            kda: 75.0,
                            relative_strength: 0.92,
                        },
                        ProteinLadderBand {
                            kda: 100.0,
                            relative_strength: 0.74,
                        },
                    ],
                },
                ProteinLadder {
                    name: "Protein Ladder 25-200 kDa",
                    bands: &[
                        ProteinLadderBand {
                            kda: 25.0,
                            relative_strength: 0.60,
                        },
                        ProteinLadderBand {
                            kda: 37.0,
                            relative_strength: 0.76,
                        },
                        ProteinLadderBand {
                            kda: 50.0,
                            relative_strength: 0.90,
                        },
                        ProteinLadderBand {
                            kda: 75.0,
                            relative_strength: 1.00,
                        },
                        ProteinLadderBand {
                            kda: 100.0,
                            relative_strength: 0.94,
                        },
                        ProteinLadderBand {
                            kda: 150.0,
                            relative_strength: 0.80,
                        },
                        ProteinLadderBand {
                            kda: 200.0,
                            relative_strength: 0.66,
                        },
                    ],
                },
            ],
        }
    }

    fn names_sorted(&self) -> Vec<String> {
        let mut names = self
            .ladders
            .iter()
            .map(|l| l.name.to_string())
            .collect::<Vec<_>>();
        names.sort_unstable();
        names
    }

    fn get(&self, name: &str) -> Option<&ProteinLadder> {
        self.ladders
            .iter()
            .find(|ladder| ladder.name.eq_ignore_ascii_case(name.trim()))
    }

    fn range_uncovered(lo: f32, hi: f32, min_kda: f32, max_kda: f32) -> f32 {
        let uncovered_low = if lo > min_kda { lo - min_kda } else { 0.0 };
        let uncovered_high = if hi < max_kda { max_kda - hi } else { 0.0 };
        uncovered_low + uncovered_high
    }

    fn choose_for_range(&self, min_kda: f32, max_kda: f32, max_ladders: usize) -> Vec<String> {
        if self.ladders.is_empty() || max_ladders == 0 || min_kda <= 0.0 || max_kda <= 0.0 {
            return vec![];
        }
        let (min_kda, max_kda) = if min_kda <= max_kda {
            (min_kda, max_kda)
        } else {
            (max_kda, min_kda)
        };
        let mut ladders = self
            .ladders
            .iter()
            .filter_map(|ladder| {
                let lo = ladder.bands.iter().map(|band| band.kda).reduce(f32::min)?;
                let hi = ladder.bands.iter().map(|band| band.kda).reduce(f32::max)?;
                Some((ladder.name.to_string(), lo, hi))
            })
            .collect::<Vec<_>>();
        ladders.sort_by(|a, b| a.0.cmp(&b.0));

        let mut best_primary: Option<(String, f32, f32, f32)> = None;
        for (name, lo, hi) in &ladders {
            let uncovered = Self::range_uncovered(*lo, *hi, min_kda, max_kda);
            let span_delta = (*lo - min_kda).abs() + (*hi - max_kda).abs();
            let score = uncovered * 1_000.0 + span_delta;
            match &best_primary {
                None => best_primary = Some((name.clone(), *lo, *hi, score)),
                Some((_n, _l, _h, best_score)) if score < *best_score => {
                    best_primary = Some((name.clone(), *lo, *hi, score))
                }
                _ => {}
            }
        }

        let Some((primary_name, primary_lo, primary_hi, primary_score)) = best_primary else {
            return vec![];
        };
        if max_ladders == 1 {
            return vec![primary_name];
        }
        if Self::range_uncovered(primary_lo, primary_hi, min_kda, max_kda) == 0.0 {
            return vec![primary_name];
        }

        let mut best_secondary: Option<(String, f32)> = None;
        for (name, lo, hi) in ladders.iter().filter(|(name, _, _)| *name != primary_name) {
            let combined_lo = primary_lo.min(*lo);
            let combined_hi = primary_hi.max(*hi);
            let uncovered = Self::range_uncovered(combined_lo, combined_hi, min_kda, max_kda);
            let span_delta = (combined_lo - min_kda).abs() + (combined_hi - max_kda).abs();
            let score = uncovered * 1_000.0 + span_delta + primary_score / 10.0;
            match &best_secondary {
                None => best_secondary = Some((name.clone(), score)),
                Some((_n, best_score)) if score < *best_score => {
                    best_secondary = Some((name.clone(), score))
                }
                _ => {}
            }
        }

        if let Some((secondary_name, _)) = best_secondary {
            vec![primary_name, secondary_name]
        } else {
            vec![primary_name]
        }
    }
}

fn normalize_ladder_name(name: &str) -> String {
    name.trim().to_ascii_lowercase()
}

fn ladder_lane_names_for_display(selected_ladders: &[String]) -> Vec<String> {
    match selected_ladders {
        [] => vec![],
        [name] => vec![name.clone(), name.clone()],
        [left, right] => vec![left.clone(), right.clone()],
        more => {
            let mut names = vec![more[0].clone()];
            names.extend(more[1..more.len() - 1].iter().cloned());
            names.push(more[more.len() - 1].clone());
            names
        }
    }
}

fn format_kda_label(kda: f32) -> String {
    if (kda - kda.round()).abs() < 0.05 {
        format!("{} kDa", kda.round() as usize)
    } else {
        format!("{kda:.1} kDa")
    }
}

fn protein_mass_y_for_range(
    range_min_kda: f32,
    range_max_kda: f32,
    kda: f32,
    top: f32,
    bottom: f32,
) -> f32 {
    let min_kda = range_min_kda.max(0.1) as f64;
    let max_kda = range_max_kda.max(range_min_kda + 0.1) as f64;
    let kda = kda.clamp(range_min_kda.max(0.1), range_max_kda.max(0.2)) as f64;
    let log_min = min_kda.log10();
    let log_max = max_kda.log10();
    let denom = (log_max - log_min).max(1e-6);
    let f = ((log_max - kda.log10()) / denom) as f32;
    top + f.clamp(0.0, 1.0) * (bottom - top)
}

fn protein_mass_y(layout: &ProteinGelLayout, kda: f32, top: f32, bottom: f32) -> f32 {
    protein_mass_y_for_range(layout.range_min_kda, layout.range_max_kda, kda, top, bottom)
}

fn resolve_ladder_names(requested: &[String], min_kda: f32, max_kda: f32) -> Vec<String> {
    let names = PROTEIN_LADDERS.names_sorted();
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

    PROTEIN_LADDERS.choose_for_range(min_kda, max_kda, 2)
}

fn protein_ladder_band(kda: f32, intensity: f32) -> ProteinGelBand {
    ProteinGelBand {
        kda,
        intensity,
        label: format_kda_label(kda),
    }
}

pub fn build_protein_gel_layout(
    samples: &[ProteinGelSample],
    requested_ladders: &[String],
    notes: Vec<String>,
) -> Result<ProteinGelLayout, String> {
    if samples.is_empty() {
        return Err("Protein gel needs at least one sample lane".to_string());
    }

    let mut normalized_samples: Vec<ProteinGelSample> = vec![];
    let mut all_band_kdas: Vec<f32> = vec![];
    for (sample_idx, sample) in samples.iter().enumerate() {
        let name = sample.name.trim().to_string();
        if name.is_empty() {
            return Err(format!(
                "Protein sample lane {} has an empty name",
                sample_idx + 1
            ));
        }
        let kda = sample.molecular_weight_kda;
        if !kda.is_finite() || kda <= 0.0 {
            return Err(format!(
                "Protein sample lane '{}' has an invalid molecular weight ({kda})",
                name
            ));
        }
        all_band_kdas.push(kda);
        normalized_samples.push(ProteinGelSample {
            name,
            detail: sample.detail.clone(),
            molecular_weight_kda: kda,
        });
    }

    let pool_min = all_band_kdas
        .iter()
        .copied()
        .reduce(f32::min)
        .unwrap_or(1.0);
    let pool_max = all_band_kdas
        .iter()
        .copied()
        .reduce(f32::max)
        .unwrap_or(pool_min);
    let selected_ladders = resolve_ladder_names(requested_ladders, pool_min, pool_max);
    if selected_ladders.is_empty() {
        return Err("No protein ladders available for protein-gel rendering".to_string());
    }

    let mut lanes: Vec<ProteinGelLane> = vec![];
    let mut all_band_values: Vec<f32> = all_band_kdas.clone();

    let display_ladders = ladder_lane_names_for_display(&selected_ladders);
    let left_ladders = display_ladders
        .iter()
        .take(display_ladders.len().div_ceil(2))
        .cloned()
        .collect::<Vec<_>>();
    let right_ladders = display_ladders
        .iter()
        .skip(display_ladders.len().div_ceil(2))
        .cloned()
        .collect::<Vec<_>>();

    for ladder_name in &left_ladders {
        let Some(ladder) = PROTEIN_LADDERS.get(ladder_name) else {
            continue;
        };
        let bands = ladder
            .bands
            .iter()
            .map(|band| {
                all_band_values.push(band.kda);
                protein_ladder_band(
                    band.kda,
                    (0.20 + 0.80 * band.relative_strength).clamp(0.20, 1.0),
                )
            })
            .collect::<Vec<_>>();
        lanes.push(ProteinGelLane {
            name: ladder.name.to_string(),
            detail: Some("protein ladder".to_string()),
            is_ladder: true,
            bands,
        });
    }
    if lanes.is_empty() {
        return Err("None of the selected protein ladders could be resolved".to_string());
    }

    for sample in &normalized_samples {
        lanes.push(ProteinGelLane {
            name: sample.name.clone(),
            detail: sample.detail.clone(),
            is_ladder: false,
            bands: vec![ProteinGelBand {
                kda: sample.molecular_weight_kda,
                intensity: 0.88,
                label: format_kda_label(sample.molecular_weight_kda),
            }],
        });
    }

    for ladder_name in &right_ladders {
        let Some(ladder) = PROTEIN_LADDERS.get(ladder_name) else {
            continue;
        };
        let bands = ladder
            .bands
            .iter()
            .map(|band| {
                all_band_values.push(band.kda);
                protein_ladder_band(
                    band.kda,
                    (0.20 + 0.80 * band.relative_strength).clamp(0.20, 1.0),
                )
            })
            .collect::<Vec<_>>();
        lanes.push(ProteinGelLane {
            name: ladder.name.to_string(),
            detail: Some("protein ladder".to_string()),
            is_ladder: true,
            bands,
        });
    }

    let min_band = all_band_values
        .iter()
        .copied()
        .reduce(f32::min)
        .unwrap_or(pool_min);
    let max_band = all_band_values
        .iter()
        .copied()
        .reduce(f32::max)
        .unwrap_or(pool_max);
    let range_min_kda = ((min_band as f64) * 0.72).floor().max(1.0) as f32;
    let mut range_max_kda = ((max_band as f64) * 1.28).ceil().max(2.0) as f32;
    if range_max_kda <= range_min_kda {
        range_max_kda = range_min_kda + 1.0;
    }

    Ok(ProteinGelLayout {
        lanes,
        selected_ladders,
        sample_count: normalized_samples.len(),
        protein_count: normalized_samples.len(),
        range_min_kda,
        range_max_kda,
        notes,
    })
}

pub fn export_protein_gel_svg(layout: &ProteinGelLayout) -> String {
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
                .set("fill", "#101317"),
        )
        .add(
            Rectangle::new()
                .set("x", DETAIL_PANEL_LEFT - 18.0)
                .set("y", GEL_TOP)
                .set("width", DETAIL_PANEL_WIDTH + 18.0)
                .set("height", gel_height)
                .set("rx", 10)
                .set("ry", 10)
                .set("fill", "#f3f4f6"),
        );

    let mut tick_kdas: Vec<f32> = vec![];
    for lane in layout.lanes.iter().filter(|l| l.is_ladder) {
        for band in &lane.bands {
            if !tick_kdas
                .iter()
                .any(|existing| (*existing - band.kda).abs() < 0.05)
            {
                tick_kdas.push(band.kda);
            }
        }
    }
    if tick_kdas.is_empty() {
        tick_kdas.push(layout.range_min_kda);
        tick_kdas.push(layout.range_max_kda);
    }
    tick_kdas.sort_by(|a, b| a.total_cmp(b));
    let mut accepted_ticks: Vec<f32> = vec![];
    let mut last_y: Option<f32> = None;
    for kda in tick_kdas.iter().rev() {
        let y = protein_mass_y(layout, *kda, GEL_TOP, GEL_BOTTOM);
        if last_y.map(|v| (v - y).abs() >= 16.0).unwrap_or(true) {
            accepted_ticks.push(*kda);
            last_y = Some(y);
        }
        if accepted_ticks.len() >= 20 {
            break;
        }
    }
    for kda in accepted_ticks {
        let y = protein_mass_y(layout, kda, GEL_TOP, GEL_BOTTOM);
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
                Text::new(format_kda_label(kda))
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
                    .set("font-size", 11)
                    .set("fill", "#0f172a"),
            );

        if let Some(detail) = lane.detail.as_deref() {
            doc = doc.add(
                Text::new(detail.to_string())
                    .set("x", x)
                    .set("y", GEL_BOTTOM + 40.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#334155"),
            );
        }

        for band in &lane.bands {
            let y = protein_mass_y(layout, band.kda, GEL_TOP + 14.0, GEL_BOTTOM - 14.0);
            let width = if lane.is_ladder {
                34.0 + 20.0 * band.intensity
            } else {
                42.0
            };
            let height = if lane.is_ladder {
                3.0 + 3.0 * band.intensity
            } else {
                8.0
            };
            let fill = if lane.is_ladder { "#e5e7eb" } else { "#d97706" };
            doc = doc.add(
                Rectangle::new()
                    .set("x", x - width * 0.5)
                    .set("y", y - height * 0.5)
                    .set("width", width)
                    .set("height", height)
                    .set("rx", 2)
                    .set("ry", 2)
                    .set("fill", fill)
                    .set("opacity", (0.45 + 0.55 * band.intensity).clamp(0.35, 1.0)),
            );
            if !lane.is_ladder {
                doc = doc.add(
                    Text::new(band.label.clone())
                        .set("x", x + 40.0)
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
        "Protein Gel Preview ({} sample lane(s), {} protein(s)) | ladder: {}",
        layout.sample_count, layout.protein_count, ladder_caption
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
                "range {}..{} | lanes: {} | notes: {}",
                format_kda_label(layout.range_min_kda),
                format_kda_label(layout.range_max_kda),
                lane_count,
                layout.notes.len()
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
                .set("stroke", "#d97706")
                .set("stroke-width", 1.5)
                .set("opacity", 0.25),
        );
    }

    let mut header_y = DETAIL_PANEL_TOP + 24.0;
    doc = doc
        .add(
            Text::new("Selection notes")
                .set("x", DETAIL_PANEL_LEFT)
                .set("y", header_y)
                .set("font-family", "monospace")
                .set("font-size", 12)
                .set("font-weight", 700)
                .set("fill", "#0f172a"),
        )
        .add(
            Text::new("Protein products are derived from the selected transcript features")
                .set("x", DETAIL_PANEL_LEFT + 4.0)
                .set("y", header_y + 16.0)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#475569"),
        );
    header_y += 36.0;

    for (idx, line) in layout.notes.iter().enumerate() {
        doc = doc.add(
            Text::new(line.clone())
                .set("x", DETAIL_PANEL_LEFT + 4.0)
                .set("y", header_y + idx as f32 * 14.0)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#111827"),
        );
    }
    let lane_details_start = header_y + layout.notes.len() as f32 * 14.0 + 22.0;
    doc = doc.add(
        Text::new("Lane details")
            .set("x", DETAIL_PANEL_LEFT)
            .set("y", lane_details_start)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("font-weight", 700)
            .set("fill", "#0f172a"),
    );
    let mut lane_detail_y = lane_details_start + 16.0;
    for lane in layout.lanes.iter().filter(|lane| !lane.is_ladder) {
        if lane_detail_y > GEL_BOTTOM - 18.0 {
            break;
        }
        let detail = lane
            .detail
            .as_deref()
            .filter(|text| !text.trim().is_empty())
            .unwrap_or("");
        let line = if detail.is_empty() {
            format!("{} | {}", lane.name, lane.bands[0].label)
        } else {
            format!("{} | {} | {}", lane.name, detail, lane.bands[0].label)
        };
        doc = doc.add(
            Text::new(line)
                .set("x", DETAIL_PANEL_LEFT + 4.0)
                .set("y", lane_detail_y)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#111827"),
        );
        lane_detail_y += 14.0;
    }

    doc.to_string()
}

#[derive(Clone, Debug)]
pub struct Protein2dGelSpot {
    pub name: String,
    pub detail: Option<String>,
    pub molecular_weight_kda: f32,
    pub isoelectric_point: f32,
}

#[derive(Clone, Debug)]
pub struct Protein2dGelLayout {
    pub spots: Vec<Protein2dGelSpot>,
    pub selected_ladders: Vec<String>,
    pub sample_count: usize,
    pub protein_count: usize,
    pub range_min_pi: f32,
    pub range_max_pi: f32,
    pub range_min_kda: f32,
    pub range_max_kda: f32,
    pub notes: Vec<String>,
}

fn format_pi_label(pi: f32) -> String {
    if (pi - pi.round()).abs() < 0.05 {
        format!("{} pI", pi.round() as usize)
    } else {
        format!("{pi:.1} pI")
    }
}

fn protein_pi_x(layout: &Protein2dGelLayout, pi: f32, left: f32, right: f32) -> f32 {
    let min_pi = layout.range_min_pi.min(layout.range_max_pi);
    let max_pi = layout.range_max_pi.max(layout.range_min_pi);
    let pi = pi.clamp(min_pi, max_pi);
    let span = (max_pi - min_pi).max(0.1);
    left + ((pi - min_pi) / span) * (right - left)
}

fn protein_pi_tick_step(layout: &Protein2dGelLayout) -> f32 {
    let span = (layout.range_max_pi - layout.range_min_pi).abs();
    if span <= 6.0 { 0.5 } else { 1.0 }
}

fn protein_pi_tick_values(layout: &Protein2dGelLayout) -> Vec<f32> {
    let step = protein_pi_tick_step(layout);
    let start = (layout.range_min_pi / step).ceil() * step;
    let mut ticks = vec![];
    let mut tick = start;
    let mut guard = 0usize;
    while tick <= layout.range_max_pi + 0.000_1 && guard < 40 {
        ticks.push((tick * 2.0).round() / 2.0);
        tick += step;
        guard += 1;
    }
    if ticks.is_empty() {
        ticks.push(layout.range_min_pi);
        if layout.range_max_pi > layout.range_min_pi {
            ticks.push(layout.range_max_pi);
        }
    }
    ticks
}

fn protein_spot_fill_color(layout: &Protein2dGelLayout, pi: f32) -> String {
    let min_pi = layout.range_min_pi.min(layout.range_max_pi);
    let max_pi = layout.range_max_pi.max(layout.range_min_pi);
    let span = (max_pi - min_pi).max(0.1);
    let t = ((pi - min_pi) / span).clamp(0.0, 1.0);
    let start = [37.0, 99.0, 235.0];
    let end = [249.0, 115.0, 22.0];
    let rgb = [
        start[0] + (end[0] - start[0]) * t,
        start[1] + (end[1] - start[1]) * t,
        start[2] + (end[2] - start[2]) * t,
    ];
    format!(
        "#{:02x}{:02x}{:02x}",
        rgb[0].round() as u8,
        rgb[1].round() as u8,
        rgb[2].round() as u8
    )
}

pub fn build_protein_2d_gel_layout(
    samples: &[Protein2dGelSpot],
    requested_ladders: &[String],
    notes: Vec<String>,
) -> Result<Protein2dGelLayout, String> {
    if samples.is_empty() {
        return Err("Protein 2D gel needs at least one sample spot".to_string());
    }

    let mut normalized_samples: Vec<Protein2dGelSpot> = vec![];
    let mut all_band_kdas: Vec<f32> = vec![];
    let mut all_pi_values: Vec<f32> = vec![];
    for (sample_idx, sample) in samples.iter().enumerate() {
        let name = sample.name.trim().to_string();
        if name.is_empty() {
            return Err(format!(
                "Protein 2D sample spot {} has an empty name",
                sample_idx + 1
            ));
        }
        let kda = sample.molecular_weight_kda;
        if !kda.is_finite() || kda <= 0.0 {
            return Err(format!(
                "Protein 2D sample spot '{}' has an invalid molecular weight ({kda})",
                name
            ));
        }
        let pi = sample.isoelectric_point;
        if !pi.is_finite() {
            return Err(format!(
                "Protein 2D sample spot '{}' has an invalid pI ({pi})",
                name
            ));
        }
        let pi = pi.clamp(0.0, 14.0);
        all_band_kdas.push(kda);
        all_pi_values.push(pi);
        normalized_samples.push(Protein2dGelSpot {
            name,
            detail: sample.detail.clone(),
            molecular_weight_kda: kda,
            isoelectric_point: pi,
        });
    }

    let pool_min = all_band_kdas
        .iter()
        .copied()
        .reduce(f32::min)
        .unwrap_or(1.0);
    let pool_max = all_band_kdas
        .iter()
        .copied()
        .reduce(f32::max)
        .unwrap_or(pool_min);
    let selected_ladders = resolve_ladder_names(requested_ladders, pool_min, pool_max);
    if selected_ladders.is_empty() {
        return Err("No protein ladders available for protein 2D-gel rendering".to_string());
    }

    let min_band = all_band_kdas
        .iter()
        .copied()
        .reduce(f32::min)
        .unwrap_or(pool_min);
    let max_band = all_band_kdas
        .iter()
        .copied()
        .reduce(f32::max)
        .unwrap_or(pool_max);
    let range_min_kda = ((min_band as f64) * 0.72).floor().max(1.0) as f32;
    let mut range_max_kda = ((max_band as f64) * 1.28).ceil().max(2.0) as f32;
    if range_max_kda <= range_min_kda {
        range_max_kda = range_min_kda + 1.0;
    }

    let min_pi = all_pi_values
        .iter()
        .copied()
        .reduce(f32::min)
        .unwrap_or(0.0);
    let max_pi = all_pi_values
        .iter()
        .copied()
        .reduce(f32::max)
        .unwrap_or(min_pi);
    let range_min_pi = (min_pi - 1.0).floor().clamp(0.0, 13.0);
    let mut range_max_pi = (max_pi + 1.0).ceil().clamp(range_min_pi + 2.0, 14.0);
    if range_max_pi <= range_min_pi {
        range_max_pi = range_min_pi + 2.0;
    }

    Ok(Protein2dGelLayout {
        spots: normalized_samples,
        selected_ladders,
        sample_count: samples.len(),
        protein_count: samples.len(),
        range_min_pi,
        range_max_pi,
        range_min_kda,
        range_max_kda,
        notes,
    })
}

pub fn export_protein_2d_gel_svg(layout: &Protein2dGelLayout) -> String {
    let plot_left = GEL_LEFT;
    let plot_right = GEL_RIGHT;
    let plot_top = GEL_TOP;
    let plot_bottom = GEL_BOTTOM;
    let plot_width = plot_right - plot_left;
    let plot_height = plot_bottom - plot_top;
    let pi_ticks = protein_pi_tick_values(layout);
    let mut kda_ticks: Vec<f32> = vec![];
    for ladder_name in &layout.selected_ladders {
        let Some(ladder) = PROTEIN_LADDERS.get(ladder_name) else {
            continue;
        };
        for band in ladder.bands {
            if !kda_ticks
                .iter()
                .any(|existing| (*existing - band.kda).abs() < 0.05)
            {
                kda_ticks.push(band.kda);
            }
        }
    }
    if kda_ticks.is_empty() {
        kda_ticks.push(layout.range_min_kda);
        kda_ticks.push(layout.range_max_kda);
    }
    kda_ticks.sort_by(|a, b| a.total_cmp(b));
    let mut accepted_kda_ticks: Vec<f32> = vec![];
    let mut last_y: Option<f32> = None;
    for kda in kda_ticks.iter().rev() {
        let y = protein_mass_y_for_range(
            layout.range_min_kda,
            layout.range_max_kda,
            *kda,
            plot_top,
            plot_bottom,
        );
        if last_y.map(|v| (v - y).abs() >= 16.0).unwrap_or(true) {
            accepted_kda_ticks.push(*kda);
            last_y = Some(y);
        }
        if accepted_kda_ticks.len() >= 20 {
            break;
        }
    }

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
                .set("x", plot_left)
                .set("y", plot_top)
                .set("width", plot_width)
                .set("height", plot_height)
                .set("rx", 10)
                .set("ry", 10)
                .set("fill", "#101317"),
        )
        .add(
            Rectangle::new()
                .set("x", DETAIL_PANEL_LEFT - 18.0)
                .set("y", plot_top)
                .set("width", DETAIL_PANEL_WIDTH + 18.0)
                .set("height", plot_height)
                .set("rx", 10)
                .set("ry", 10)
                .set("fill", "#f3f4f6"),
        );

    for pi in &pi_ticks {
        let x = protein_pi_x(layout, *pi, plot_left, plot_right);
        doc = doc
            .add(
                Line::new()
                    .set("x1", x)
                    .set("y1", plot_top)
                    .set("x2", x)
                    .set("y2", plot_bottom)
                    .set("stroke", "#2d3238")
                    .set("stroke-width", 1),
            )
            .add(
                Text::new(format_pi_label(*pi))
                    .set("x", x)
                    .set("y", plot_bottom + 18.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 11)
                    .set("fill", "#334155"),
            );
    }

    for kda in &accepted_kda_ticks {
        let y = protein_mass_y_for_range(
            layout.range_min_kda,
            layout.range_max_kda,
            *kda,
            plot_top,
            plot_bottom,
        );
        doc = doc
            .add(
                Line::new()
                    .set("x1", plot_left)
                    .set("y1", y)
                    .set("x2", plot_right)
                    .set("y2", y)
                    .set("stroke", "#2d3238")
                    .set("stroke-width", 1),
            )
            .add(
                Text::new(format_kda_label(*kda))
                    .set("x", plot_right + 12.0)
                    .set("y", y + 4.0)
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("fill", "#374151"),
            );
    }

    for (idx, spot) in layout.spots.iter().enumerate() {
        let x = protein_pi_x(
            layout,
            spot.isoelectric_point,
            plot_left + 22.0,
            plot_right - 22.0,
        );
        let y = protein_mass_y_for_range(
            layout.range_min_kda,
            layout.range_max_kda,
            spot.molecular_weight_kda,
            plot_top + 18.0,
            plot_bottom - 18.0,
        );
        let fill = protein_spot_fill_color(layout, spot.isoelectric_point);
        doc = doc
            .add(
                Circle::new()
                    .set("cx", x)
                    .set("cy", y)
                    .set("r", 13.0)
                    .set("fill", fill)
                    .set("stroke", "#f8fafc")
                    .set("stroke-width", 2.0)
                    .set("opacity", 0.94),
            )
            .add(
                Text::new(format!("{}", idx + 1))
                    .set("x", x)
                    .set("y", y + 4.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 10)
                    .set("font-weight", 700)
                    .set("fill", "#ffffff"),
            );
    }

    let ladder_caption = if layout.selected_ladders.is_empty() {
        "auto ladder".to_string()
    } else {
        layout.selected_ladders.join(" + ")
    };
    let title = format!(
        "Protein 2D Gel Preview ({} sample spot(s)) | axis refs: {}",
        layout.sample_count, ladder_caption
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
                "x=pI | y=log(kDa) | pI {}..{} | kDa {}..{} | notes: {}",
                format_pi_label(layout.range_min_pi),
                format_pi_label(layout.range_max_pi),
                format_kda_label(layout.range_min_kda),
                format_kda_label(layout.range_max_kda),
                layout.notes.len()
            ))
            .set("x", GEL_LEFT)
            .set("y", 62.0)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#334155"),
        );

    let mut header_y = DETAIL_PANEL_TOP + 24.0;
    doc = doc
        .add(
            Text::new("Selection notes")
                .set("x", DETAIL_PANEL_LEFT)
                .set("y", header_y)
                .set("font-family", "monospace")
                .set("font-size", 12)
                .set("font-weight", 700)
                .set("fill", "#0f172a"),
        )
        .add(
            Text::new("Protein spots are derived from the selected transcript features")
                .set("x", DETAIL_PANEL_LEFT + 4.0)
                .set("y", header_y + 16.0)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#475569"),
        );
    header_y += 36.0;

    for (idx, line) in layout.notes.iter().enumerate() {
        doc = doc.add(
            Text::new(line.clone())
                .set("x", DETAIL_PANEL_LEFT + 4.0)
                .set("y", header_y + idx as f32 * 14.0)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#111827"),
        );
    }
    let spot_details_start = header_y + layout.notes.len() as f32 * 14.0 + 22.0;
    doc = doc.add(
        Text::new("Spot details")
            .set("x", DETAIL_PANEL_LEFT)
            .set("y", spot_details_start)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("font-weight", 700)
            .set("fill", "#0f172a"),
    );
    let mut spot_detail_y = spot_details_start + 16.0;
    for (idx, spot) in layout.spots.iter().enumerate() {
        if spot_detail_y > plot_bottom - 18.0 {
            break;
        }
        let detail = spot
            .detail
            .as_deref()
            .filter(|text| !text.trim().is_empty())
            .unwrap_or("");
        let line = if detail.is_empty() {
            format!(
                "{:>2}. {} | {} | {}",
                idx + 1,
                spot.name,
                format_pi_label(spot.isoelectric_point),
                format_kda_label(spot.molecular_weight_kda)
            )
        } else {
            format!(
                "{:>2}. {} | {} | {} | {}",
                idx + 1,
                spot.name,
                detail,
                format_pi_label(spot.isoelectric_point),
                format_kda_label(spot.molecular_weight_kda)
            )
        };
        doc = doc.add(
            Text::new(line)
                .set("x", DETAIL_PANEL_LEFT + 4.0)
                .set("y", spot_detail_y)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#111827"),
        );
        spot_detail_y += 14.0;
    }

    doc.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_gel_layout_renders_selection_notes_and_labels() {
        let layout = build_protein_gel_layout(
            &[
                ProteinGelSample {
                    name: "NM_005427.4".to_string(),
                    detail: Some("636 aa".to_string()),
                    molecular_weight_kda: 68.7,
                },
                ProteinGelSample {
                    name: "NM_001204184.2".to_string(),
                    detail: Some("617 aa".to_string()),
                    molecular_weight_kda: 66.9,
                },
            ],
            &[],
            vec![
                "Source: test_files/tp73.ncbi.gb".to_string(),
                "Selection: gene=TP73 / transcript_id /^NM_|^XM_/".to_string(),
            ],
        )
        .expect("protein gel layout");
        let svg = export_protein_gel_svg(&layout);
        assert!(svg.contains("Protein Gel Preview"));
        assert!(svg.contains("Selection notes"));
        assert!(svg.contains("Lane details"));
        assert!(svg.contains("NM_005427.4"));
        assert!(svg.contains("kDa"));
    }

    #[test]
    fn protein_2d_gel_layout_renders_pis_and_spots() {
        let layout = build_protein_2d_gel_layout(
            &[
                Protein2dGelSpot {
                    name: "NM_005427.4".to_string(),
                    detail: Some("636 aa".to_string()),
                    molecular_weight_kda: 68.7,
                    isoelectric_point: 6.3,
                },
                Protein2dGelSpot {
                    name: "NM_001204184.2".to_string(),
                    detail: Some("617 aa".to_string()),
                    molecular_weight_kda: 66.9,
                    isoelectric_point: 6.7,
                },
            ],
            &[],
            vec![
                "Source: test_files/tp73.ncbi.gb".to_string(),
                "Selection: gene=TP73 / transcript_id /^NM_/".to_string(),
            ],
        )
        .expect("protein 2D gel layout");
        let svg = export_protein_2d_gel_svg(&layout);
        assert!(svg.contains("Protein 2D Gel Preview"));
        assert!(svg.contains("Selection notes"));
        assert!(svg.contains("Spot details"));
        assert!(svg.contains("NM_005427.4"));
        assert!(svg.contains("pI"));
        assert!(svg.contains("kDa"));
    }
}
