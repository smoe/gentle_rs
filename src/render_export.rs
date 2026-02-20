use crate::{
    dna_sequence::DNAsequence, engine::DisplaySettings, feature_location::feature_is_reverse,
    restriction_enzyme::RestrictionEnzymeKey,
};
use gb_io::seq::Feature;
use std::collections::{HashMap, HashSet};
use svg::node::element::path::Data;
use svg::node::element::{Circle, Line, Path, Rectangle, Text};
use svg::Document;

const W: f32 = 1200.0;
const H: f32 = 700.0;
const RE_LABEL_BASE_OFFSET: f32 = 76.0;
const RE_LABEL_ROW_HEIGHT: f32 = 12.0;
const FEATURE_SIDE_MARGIN: f32 = 28.0;
const FEATURE_LANE_GAP: f32 = 16.0;
const FEATURE_BLOCK_HEIGHT: f32 = 10.0;
const REGULATORY_SIDE_MARGIN: f32 = 6.0;
const REGULATORY_LANE_GAP: f32 = 8.0;
const REGULATORY_BLOCK_HEIGHT: f32 = 7.0;
const REGULATORY_GROUP_GAP: f32 = 8.0;

#[derive(Clone, Debug)]
struct FeatureVm {
    from: usize,
    to: usize,
    label: String,
    color: &'static str,
    is_reverse: bool,
    is_pointy: bool,
    is_regulatory: bool,
}

fn feature_name(feature: &Feature) -> String {
    for k in [
        "label",
        "name",
        "standard_name",
        "gene",
        "protein_id",
        "product",
        "region_name",
        "bound_moiety",
    ] {
        if let Some(s) = feature.qualifier_values(k.into()).next() {
            let label_text = s.to_string();
            if !label_text.trim().is_empty() {
                return label_text;
            }
        }
    }
    if is_regulatory_feature(feature) {
        return String::new();
    }
    match feature.location.find_bounds() {
        Ok((from, to)) => format!("{from}..{to}"),
        Err(_) => String::new(),
    }
}

fn feature_color(feature: &Feature) -> &'static str {
    if is_vcf_track_feature(feature) {
        let class = vcf_variant_class(feature)
            .unwrap_or_else(|| "OTHER".to_string())
            .to_ascii_uppercase();
        return match class.as_str() {
            "SNP" => "#e17f0f",
            "INS" => "#238c64",
            "DEL" => "#b42d2d",
            "SV" => "#5a5a1e",
            _ => "#5a5a5a",
        };
    }
    match feature.kind.to_string().to_ascii_uppercase().as_str() {
        "CDS" => "#cc1f1f",
        "GENE" => "#1f4fcc",
        "MRNA" => "#b4640a",
        "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND" => "#238023",
        _ => "#6e6e6e",
    }
}

fn feature_pointy(feature: &Feature) -> bool {
    matches!(
        feature.kind.to_string().to_ascii_uppercase().as_str(),
        "CDS" | "GENE"
    )
}

fn is_tfbs_feature(feature: &Feature) -> bool {
    matches!(
        feature.kind.to_string().to_ascii_uppercase().as_str(),
        "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND"
    )
}

fn is_track_feature(feature: &Feature) -> bool {
    feature
        .qualifier_values("gentle_generated".into())
        .any(|value| {
            matches!(
                value.trim().to_ascii_lowercase().as_str(),
                "genome_bed_track" | "genome_bigwig_track" | "genome_vcf_track" | "blast_hit_track"
            )
        })
}

fn is_vcf_track_feature(feature: &Feature) -> bool {
    feature
        .qualifier_values("gentle_generated".into())
        .any(|value| value.trim().eq_ignore_ascii_case("genome_vcf_track"))
}

fn is_regulatory_feature(feature: &Feature) -> bool {
    let kind = feature.kind.to_string().to_ascii_uppercase();
    is_tfbs_feature(feature) || kind.contains("REGULATORY") || is_track_feature(feature)
}

fn feature_qualifier_f64(feature: &Feature, key: &str) -> Option<f64> {
    feature
        .qualifier_values(key.into())
        .next()
        .and_then(|v| v.trim().parse::<f64>().ok())
}

fn feature_qualifier_text(feature: &Feature, key: &str) -> Option<String> {
    feature
        .qualifier_values(key.into())
        .map(str::trim)
        .find(|value| !value.is_empty())
        .map(ToOwned::to_owned)
}

fn vcf_variant_class(feature: &Feature) -> Option<String> {
    if !is_vcf_track_feature(feature) {
        return None;
    }
    feature_qualifier_text(feature, "vcf_variant_class").or_else(|| {
        let reference = feature_qualifier_text(feature, "vcf_ref")?;
        let alt = feature_qualifier_text(feature, "vcf_alt")?;
        let ref_len = reference.trim().len().max(1);
        let alt_trimmed = alt.trim();
        if alt_trimmed.starts_with('<')
            || alt_trimmed.ends_with('>')
            || alt_trimmed.contains('[')
            || alt_trimmed.contains(']')
            || alt_trimmed == "*"
        {
            return Some("SV".to_string());
        }
        let alt_len = alt_trimmed.len().max(1);
        if ref_len == 1 && alt_len == 1 {
            Some("SNP".to_string())
        } else if alt_len > ref_len {
            Some("INS".to_string())
        } else if alt_len < ref_len {
            Some("DEL".to_string())
        } else {
            Some("OTHER".to_string())
        }
    })
}

fn tfbs_feature_passes_display_filter(feature: &Feature, display: &DisplaySettings) -> bool {
    if !is_tfbs_feature(feature) {
        return true;
    }
    if display.tfbs_display_use_llr_bits {
        let Some(value) = feature_qualifier_f64(feature, "llr_bits") else {
            return false;
        };
        if value < display.tfbs_display_min_llr_bits {
            return false;
        }
    }
    if display.tfbs_display_use_llr_quantile {
        let Some(value) = feature_qualifier_f64(feature, "llr_quantile") else {
            return false;
        };
        if value < display.tfbs_display_min_llr_quantile {
            return false;
        }
    }
    if display.tfbs_display_use_true_log_odds_bits {
        let value = feature_qualifier_f64(feature, "true_log_odds_bits")
            .or_else(|| feature_qualifier_f64(feature, "log_odds_ratio_bits"));
        let Some(value) = value else {
            return false;
        };
        if value < display.tfbs_display_min_true_log_odds_bits {
            return false;
        }
    }
    if display.tfbs_display_use_true_log_odds_quantile {
        let value = feature_qualifier_f64(feature, "true_log_odds_quantile")
            .or_else(|| feature_qualifier_f64(feature, "log_odds_ratio_quantile"));
        let Some(value) = value else {
            return false;
        };
        if value < display.tfbs_display_min_true_log_odds_quantile {
            return false;
        }
    }
    true
}

fn vcf_feature_passes_display_filter(feature: &Feature, display: &DisplaySettings) -> bool {
    if !is_vcf_track_feature(feature) {
        return true;
    }
    let class = vcf_variant_class(feature)
        .unwrap_or_else(|| "OTHER".to_string())
        .to_ascii_uppercase();
    let class_ok = match class.as_str() {
        "SNP" => display.vcf_display_show_snp,
        "INS" => display.vcf_display_show_ins,
        "DEL" => display.vcf_display_show_del,
        "SV" => display.vcf_display_show_sv,
        _ => display.vcf_display_show_other,
    };
    if !class_ok {
        return false;
    }
    if display.vcf_display_pass_only {
        let filter = feature_qualifier_text(feature, "vcf_filter")
            .unwrap_or_default()
            .trim()
            .to_ascii_uppercase();
        if filter != "PASS" {
            return false;
        }
    }
    if display.vcf_display_use_min_qual || display.vcf_display_use_max_qual {
        let qual =
            feature_qualifier_f64(feature, "vcf_qual").or_else(|| feature_qualifier_f64(feature, "score"));
        let Some(qual) = qual else {
            return false;
        };
        if display.vcf_display_use_min_qual && qual < display.vcf_display_min_qual {
            return false;
        }
        if display.vcf_display_use_max_qual && qual > display.vcf_display_max_qual {
            return false;
        }
    }
    if !display.vcf_display_required_info_keys.is_empty() {
        let info = feature_qualifier_text(feature, "vcf_info").unwrap_or_default();
        let info_keys = info
            .split(';')
            .map(str::trim)
            .filter(|entry| !entry.is_empty())
            .map(|entry| {
                entry
                    .split_once('=')
                    .map(|(key, _)| key.trim().to_ascii_uppercase())
                    .unwrap_or_else(|| entry.to_ascii_uppercase())
            })
            .collect::<HashSet<_>>();
        if display
            .vcf_display_required_info_keys
            .iter()
            .map(|key| key.trim().to_ascii_uppercase())
            .any(|key| !info_keys.contains(&key))
        {
            return false;
        }
    }
    true
}

fn collect_features(dna: &DNAsequence, display: &DisplaySettings) -> Vec<FeatureVm> {
    let mut ret = Vec::new();
    for feature in dna.features() {
        if feature.kind.to_string().to_ascii_uppercase() == "SOURCE" {
            continue;
        }
        let kind = feature.kind.to_string().to_ascii_uppercase();
        if kind == "CDS" && !display.show_cds_features {
            continue;
        }
        if kind == "GENE" && !display.show_gene_features {
            continue;
        }
        if kind == "MRNA" && !display.show_mrna_features {
            continue;
        }
        if is_tfbs_feature(feature) {
            if !display.show_tfbs {
                continue;
            }
            if !tfbs_feature_passes_display_filter(feature, display) {
                continue;
            }
        }
        if is_vcf_track_feature(feature) && !vcf_feature_passes_display_filter(feature, display) {
            continue;
        }
        let Ok((from, to)) = feature.location.find_bounds() else {
            continue;
        };
        if from < 0 || to < 0 {
            continue;
        }
        ret.push(FeatureVm {
            from: from as usize,
            to: to as usize,
            label: feature_name(feature),
            color: feature_color(feature),
            is_reverse: feature_is_reverse(feature),
            is_pointy: feature_pointy(feature),
            is_regulatory: is_regulatory_feature(feature),
        });
    }
    ret.sort_by(|a, b| {
        a.from
            .cmp(&b.from)
            .then(a.to.cmp(&b.to))
            .then(a.label.cmp(&b.label))
    });
    ret
}

fn bp_to_x(bp: usize, len: usize, left: f32, right: f32) -> f32 {
    if len == 0 {
        return left;
    }
    let frac = bp as f32 / len as f32;
    left + (right - left) * frac
}

fn lane_allocate(lanes: &mut Vec<f32>, start: f32, end: f32, padding: f32) -> usize {
    for (idx, lane_end) in lanes.iter_mut().enumerate() {
        if start >= *lane_end + padding {
            *lane_end = end;
            return idx;
        }
    }
    lanes.push(end);
    lanes.len() - 1
}

fn estimate_text_width(label: &str) -> f32 {
    (label.chars().count().max(1) as f32) * 6.5
}

fn re_color(key: &RestrictionEnzymeKey) -> &'static str {
    match key.number_of_cuts() {
        1 => "#8b0000",
        2 => "#00008b",
        3 => "#006400",
        _ => "#000000",
    }
}

pub fn export_linear_svg(dna: &DNAsequence, display: &DisplaySettings) -> String {
    let len = dna.len();
    let left = 60.0;
    let right = W - 60.0;
    let baseline = H * 0.5;

    let mut doc = Document::new()
        .set("viewBox", (0, 0, W, H))
        .set("width", W)
        .set("height", H)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", W)
                .set("height", H)
                .set("fill", "#ffffff"),
        );
    let mut labels: Vec<Text> = vec![];

    if display.show_gc_contents {
        for region in dna.gc_content().regions() {
            let x1 = bp_to_x(region.from(), len, left, right);
            let x2 = bp_to_x(region.to(), len, left, right).max(x1 + 1.0);
            let g = (region.gc() * 255.0).round() as u8;
            let r = 255u8.saturating_sub(g);
            let color = format!("#{:02x}{:02x}00", r, g);
            doc = doc.add(
                Rectangle::new()
                    .set("x", x1)
                    .set("y", baseline - 8.0)
                    .set("width", x2 - x1)
                    .set("height", 6.0)
                    .set("fill", color),
            );
        }
    }

    if display.show_methylation_sites {
        for site in dna.methylation_sites().sites() {
            let x = bp_to_x(*site, len, left, right);
            doc = doc.add(
                Line::new()
                    .set("x1", x)
                    .set("y1", baseline - 12.0)
                    .set("x2", x)
                    .set("y2", baseline - 2.0)
                    .set("stroke", "#8b0000")
                    .set("stroke-width", 1),
            );
        }
    }

    doc = doc.add(
        Line::new()
            .set("x1", left)
            .set("y1", baseline)
            .set("x2", right)
            .set("y2", baseline)
            .set("stroke", "#000000")
            .set("stroke-width", 2),
    );

    if display.show_open_reading_frames {
        for orf in dna.open_reading_frames() {
            let start = (orf.from().min(orf.to())).max(0) as usize;
            let end = (orf.from().max(orf.to())).max(0) as usize;
            let x1 = bp_to_x(start, len, left, right);
            let x2 = bp_to_x(end, len, left, right).max(x1 + 1.0);
            let level = (orf.frame().unsigned_abs()) as f32;
            let y = if orf.is_reverse() {
                baseline + 10.0 + 7.0 * level
            } else {
                baseline - 10.0 - 7.0 * level
            };
            let color = match orf.frame() {
                -1 => "#ff7777",
                -2 => "#77cc77",
                -3 => "#7777ff",
                1 => "#8b0000",
                2 => "#006400",
                3 => "#00008b",
                _ => "#444444",
            };
            doc = doc.add(
                Rectangle::new()
                    .set("x", x1)
                    .set("y", y - 2.5)
                    .set("width", x2 - x1)
                    .set("height", 5.0)
                    .set("fill", color),
            );
            let label = format!("ORF {}", orf.frame());
            if (x2 - x1) >= estimate_text_width(&label) + 6.0 {
                labels.push(
                    Text::new(label)
                        .set("x", (x1 + x2) / 2.0)
                        .set("y", y + 3.0)
                        .set("text-anchor", "middle")
                        .set("font-family", "monospace")
                        .set("font-size", 8)
                        .set("fill", "#111111"),
                );
            }
        }
    }

    if display.show_features {
        let regulatory_tracks_near_baseline = display.regulatory_tracks_near_baseline;
        let features = collect_features(dna, display);
        let mut lane_top_by_idx: Vec<usize> = vec![0; features.len()];
        let mut lane_bottom_by_idx: Vec<usize> = vec![0; features.len()];
        let mut lane_regulatory_top_by_idx: Vec<usize> = vec![0; features.len()];
        let mut top_lane_ends: Vec<f32> = vec![];
        let mut bottom_lane_ends: Vec<f32> = vec![];
        let mut regulatory_top_lane_ends: Vec<f32> = vec![];

        let mut top_order: Vec<usize> = features
            .iter()
            .enumerate()
            .filter_map(|(idx, f)| {
                if !f.is_reverse && !f.is_regulatory {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect();
        top_order.sort_by(|a, b| {
            let fa = &features[*a];
            let fb = &features[*b];
            let span_a = fa.to.saturating_sub(fa.from).saturating_add(1);
            let span_b = fb.to.saturating_sub(fb.from).saturating_add(1);
            span_b.cmp(&span_a).then(fa.from.cmp(&fb.from))
        });
        for idx in top_order {
            let f = &features[idx];
            let x1 = bp_to_x(f.from, len, left, right);
            let x2 = bp_to_x(f.to, len, left, right).max(x1 + 1.0);
            lane_top_by_idx[idx] = lane_allocate(&mut top_lane_ends, x1, x2, 4.0);
        }

        let mut bottom_order: Vec<usize> = features
            .iter()
            .enumerate()
            .filter_map(|(idx, f)| {
                if f.is_reverse && !f.is_regulatory {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect();
        bottom_order.sort_by(|a, b| {
            let fa = &features[*a];
            let fb = &features[*b];
            let span_a = fa.to.saturating_sub(fa.from).saturating_add(1);
            let span_b = fb.to.saturating_sub(fb.from).saturating_add(1);
            span_b.cmp(&span_a).then(fa.from.cmp(&fb.from))
        });
        for idx in bottom_order {
            let f = &features[idx];
            let x1 = bp_to_x(f.from, len, left, right);
            let x2 = bp_to_x(f.to, len, left, right).max(x1 + 1.0);
            lane_bottom_by_idx[idx] = lane_allocate(&mut bottom_lane_ends, x1, x2, 4.0);
        }

        let mut regulatory_top_order: Vec<usize> = features
            .iter()
            .enumerate()
            .filter_map(|(idx, f)| if f.is_regulatory { Some(idx) } else { None })
            .collect();
        regulatory_top_order.sort_by(|a, b| {
            let fa = &features[*a];
            let fb = &features[*b];
            let span_a = fa.to.saturating_sub(fa.from).saturating_add(1);
            let span_b = fb.to.saturating_sub(fb.from).saturating_add(1);
            span_b.cmp(&span_a).then(fa.from.cmp(&fb.from))
        });
        for idx in regulatory_top_order {
            let f = &features[idx];
            if regulatory_tracks_near_baseline {
                lane_regulatory_top_by_idx[idx] = 0;
            } else {
                let x1 = bp_to_x(f.from, len, left, right);
                let x2 = bp_to_x(f.to, len, left, right).max(x1 + 1.0);
                lane_regulatory_top_by_idx[idx] =
                    lane_allocate(&mut regulatory_top_lane_ends, x1, x2, 0.0);
            }
        }

        let top_regular_extent = if top_lane_ends.is_empty() {
            0.0
        } else {
            FEATURE_SIDE_MARGIN
                + (top_lane_ends.len().saturating_sub(1) as f32) * FEATURE_LANE_GAP
                + FEATURE_BLOCK_HEIGHT * 0.5
        };
        let regulatory_group_gap = if !regulatory_tracks_near_baseline
            && !top_lane_ends.is_empty()
            && !regulatory_top_lane_ends.is_empty()
        {
            REGULATORY_GROUP_GAP
        } else {
            0.0
        };

        for (idx, f) in features.iter().enumerate() {
            let x1 = bp_to_x(f.from, len, left, right);
            let x2 = bp_to_x(f.to, len, left, right).max(x1 + 1.0);
            let (y, block_height) = if f.is_regulatory {
                let y = if regulatory_tracks_near_baseline {
                    baseline - REGULATORY_SIDE_MARGIN
                } else {
                    let lane = lane_regulatory_top_by_idx[idx];
                    baseline
                        - top_regular_extent
                        - regulatory_group_gap
                        - REGULATORY_SIDE_MARGIN
                        - lane as f32 * REGULATORY_LANE_GAP
                };
                (y, REGULATORY_BLOCK_HEIGHT)
            } else if f.is_reverse {
                let lane = lane_bottom_by_idx[idx];
                let y = baseline + FEATURE_SIDE_MARGIN + lane as f32 * FEATURE_LANE_GAP;
                (y, FEATURE_BLOCK_HEIGHT)
            } else {
                let lane = lane_top_by_idx[idx];
                let y = baseline - FEATURE_SIDE_MARGIN - lane as f32 * FEATURE_LANE_GAP;
                (y, FEATURE_BLOCK_HEIGHT)
            };
            let half_height = block_height * 0.5;
            doc = doc.add(
                Rectangle::new()
                    .set("x", x1)
                    .set("y", y - half_height)
                    .set("width", x2 - x1)
                    .set("height", block_height)
                    .set("fill", f.color),
            );

            if f.is_pointy {
                let arrow_dx = 6.0;
                let data = if f.is_reverse {
                    Data::new()
                        .move_to((x1, y - half_height))
                        .line_to((x1 - arrow_dx, y))
                        .line_to((x1, y + half_height))
                        .close()
                } else {
                    Data::new()
                        .move_to((x2, y - half_height))
                        .line_to((x2 + arrow_dx, y))
                        .line_to((x2, y + half_height))
                        .close()
                };
                doc = doc.add(Path::new().set("d", data).set("fill", f.color));
            }

            if f.label.trim().is_empty() {
                continue;
            }
            let (text_y, anchor) = if f.is_regulatory {
                (y - half_height - 2.0, "start")
            } else if f.is_reverse {
                (y + 16.0, "start")
            } else {
                (y - 10.0, "start")
            };
            labels.push(
                Text::new(f.label.clone())
                    .set("x", x1)
                    .set("y", text_y)
                    .set("text-anchor", anchor)
                    .set("font-family", "monospace")
                    .set("font-size", 10)
                    .set("fill", "#111111"),
            );
        }
    }

    if display.show_restriction_enzymes {
        let mut keys: Vec<RestrictionEnzymeKey> =
            dna.restriction_enzyme_groups().keys().cloned().collect();
        keys.sort();
        let mut top_label_lanes: Vec<f32> = vec![];
        let mut bottom_label_lanes: Vec<f32> = vec![];

        for (i, key) in keys.iter().enumerate() {
            let x = bp_to_x(key.pos().max(0) as usize, len, left, right);
            let color = re_color(key);
            doc = doc.add(
                Line::new()
                    .set("x1", x)
                    .set("y1", baseline - 8.0)
                    .set("x2", x)
                    .set("y2", baseline + 8.0)
                    .set("stroke", color)
                    .set("stroke-width", 1),
            );
            if let Some(names) = dna.restriction_enzyme_groups().get(key) {
                let label = names.join(",");
                let label_width = estimate_text_width(&label);
                let label_left = x - label_width / 2.0;
                let label_right = x + label_width / 2.0;
                let place_top = i % 2 == 0;
                let lane = if place_top {
                    lane_allocate(&mut top_label_lanes, label_left, label_right, 6.0)
                } else {
                    lane_allocate(&mut bottom_label_lanes, label_left, label_right, 6.0)
                };
                let y = if place_top {
                    baseline - RE_LABEL_BASE_OFFSET - RE_LABEL_ROW_HEIGHT * lane as f32
                } else {
                    baseline + RE_LABEL_BASE_OFFSET + RE_LABEL_ROW_HEIGHT * lane as f32
                };
                labels.push(
                    Text::new(label)
                        .set("x", x)
                        .set("y", y)
                        .set("text-anchor", "middle")
                        .set("font-family", "monospace")
                        .set("font-size", 9)
                        .set("fill", color),
                );
            }
        }
    }

    labels.push(
        Text::new(
            dna.name()
                .clone()
                .unwrap_or_else(|| "<no name>".to_string()),
        )
        .set("x", 12)
        .set("y", 24)
        .set("font-family", "monospace")
        .set("font-size", 16)
        .set("fill", "#111111"),
    );
    labels.push(
        Text::new(format!("{} bp", len))
            .set("x", W - 12.0)
            .set("y", 24)
            .set("text-anchor", "end")
            .set("font-family", "monospace")
            .set("font-size", 14)
            .set("fill", "#444444"),
    );
    for label in labels {
        doc = doc.add(label);
    }

    doc.to_string()
}

fn pos2xy(pos: usize, len: usize, cx: f32, cy: f32, r: f32) -> (f32, f32) {
    if len == 0 {
        return (cx, cy);
    }
    let angle =
        2.0 * std::f32::consts::PI * (pos as f32 / len as f32) - std::f32::consts::FRAC_PI_2;
    (cx + r * angle.cos(), cy + r * angle.sin())
}

fn circular_arc_path(
    from: usize,
    to: usize,
    len: usize,
    cx: f32,
    cy: f32,
    r: f32,
) -> Option<String> {
    if len == 0 {
        return None;
    }
    let mut delta = if to >= from {
        to - from
    } else {
        (len - from) + to
    };
    if delta == 0 {
        return None;
    }
    if delta >= len {
        delta %= len;
        if delta == 0 {
            return None;
        }
    }
    let large_arc = if delta as f32 > (len as f32 / 2.0) {
        1
    } else {
        0
    };
    let sweep = 1;
    let (x1, y1) = pos2xy(from, len, cx, cy, r);
    let (x2, y2) = pos2xy(to, len, cx, cy, r);
    Some(format!(
        "M {x1:.3} {y1:.3} A {r:.3} {r:.3} 0 {large_arc} {sweep} {x2:.3} {y2:.3}"
    ))
}

pub fn export_circular_svg(dna: &DNAsequence, display: &DisplaySettings) -> String {
    let len = dna.len();
    let cx = W * 0.56;
    let cy = H * 0.52;
    let r = H.min(W) * 0.28;

    let mut doc = Document::new()
        .set("viewBox", (0, 0, W, H))
        .set("width", W)
        .set("height", H)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", W)
                .set("height", H)
                .set("fill", "#ffffff"),
        );

    doc = doc.add(
        Text::new(
            dna.name()
                .clone()
                .unwrap_or_else(|| "<no name>".to_string()),
        )
        .set("x", 20)
        .set("y", 34)
        .set("font-family", "monospace")
        .set("font-size", 16)
        .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new(format!("{} bp", len))
            .set("x", 20)
            .set("y", 54)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#444444"),
    );

    doc = doc.add(
        Circle::new()
            .set("cx", cx)
            .set("cy", cy)
            .set("r", r)
            .set("fill", "none")
            .set("stroke", "#000000")
            .set("stroke-width", 2),
    );

    if display.show_gc_contents {
        for region in dna.gc_content().regions() {
            let (x1, y1) = pos2xy(region.from(), len, cx, cy, r * 0.72);
            let (x2, y2) = pos2xy(region.to(), len, cx, cy, r * 0.72);
            let g = (region.gc() * 255.0).round() as u8;
            let rr = 255u8.saturating_sub(g);
            let color = format!("#{:02x}{:02x}00", rr, g);
            doc = doc.add(
                Line::new()
                    .set("x1", x1)
                    .set("y1", y1)
                    .set("x2", x2)
                    .set("y2", y2)
                    .set("stroke", color)
                    .set("stroke-width", 6),
            );
        }
    }

    if display.show_methylation_sites {
        for site in dna.methylation_sites().sites() {
            let (x1, y1) = pos2xy(*site, len, cx, cy, r);
            let (x2, y2) = pos2xy(*site, len, cx, cy, r * 0.9);
            doc = doc.add(
                Line::new()
                    .set("x1", x1)
                    .set("y1", y1)
                    .set("x2", x2)
                    .set("y2", y2)
                    .set("stroke", "#8b0000")
                    .set("stroke-width", 1),
            );
        }
    }

    if display.show_features {
        for f in collect_features(dna, display) {
            let mid = (f.from + f.to) / 2;
            let span = f.to.saturating_sub(f.from).saturating_add(1);
            let length_fraction = if len == 0 {
                0.0
            } else {
                (span as f32 / len as f32).clamp(0.0, 1.0)
            };
            let offset = (0.14 * (1.0 - 0.6 * length_fraction)).clamp(0.05, 0.14);
            let band = if f.is_reverse {
                1.0 - offset
            } else {
                1.0 + offset
            };
            if let Some(path_d) = circular_arc_path(f.from, f.to, len, cx, cy, r * band) {
                doc = doc.add(
                    Path::new()
                        .set("d", path_d)
                        .set("fill", "none")
                        .set("stroke", f.color)
                        .set("stroke-width", 5),
                );
            }
            let label_band = if f.is_reverse {
                1.0 - (offset + 0.08)
            } else {
                1.0 + (offset + 0.08)
            };
            let (lx, ly) = pos2xy(mid, len, cx, cy, r * label_band);
            doc = doc.add(
                Text::new(f.label)
                    .set("x", lx)
                    .set("y", ly)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 10)
                    .set("fill", "#111111"),
            );
        }
    }

    if display.show_open_reading_frames {
        for orf in dna.open_reading_frames() {
            let start = (orf.from().min(orf.to())).max(0) as usize;
            let end = (orf.from().max(orf.to())).max(0) as usize;
            let frame = orf.frame().unsigned_abs() as f32;
            let rr = if orf.is_reverse() {
                r * (0.92 - 0.03 * frame)
            } else {
                r * (1.08 + 0.03 * frame)
            };
            let color = match orf.frame() {
                -1 => "#ff7777",
                -2 => "#77cc77",
                -3 => "#7777ff",
                1 => "#8b0000",
                2 => "#006400",
                3 => "#00008b",
                _ => "#444444",
            };
            if let Some(path_d) = circular_arc_path(start, end, len, cx, cy, rr) {
                doc = doc.add(
                    Path::new()
                        .set("d", path_d)
                        .set("fill", "none")
                        .set("stroke", color)
                        .set("stroke-width", 2),
                );
            }
        }
    }

    if display.show_restriction_enzymes {
        let mut keys: Vec<RestrictionEnzymeKey> =
            dna.restriction_enzyme_groups().keys().cloned().collect();
        keys.sort();
        let mut last_right_y = -1e9f32;
        let mut last_left_y = 1e9f32;
        for key in &keys {
            let pos = key.pos().max(0) as usize;
            let font_color = re_color(key);
            let raw_label = match dna.restriction_enzyme_groups().get(key) {
                Some(names) => names.join(", "),
                None => continue,
            };
            let right_side = pos < len / 2;
            let label = if right_side {
                format!("{pos} {raw_label}")
            } else {
                format!("{raw_label} {pos}")
            };

            let (x1, y1) = pos2xy(pos, len, cx, cy, r);
            let (x2, y2) = pos2xy(pos, len, cx, cy, r * 1.15);
            let (mut x3, _) = pos2xy(pos, len, cx, cy, r * 1.25);
            let mut y3 = y2;
            if right_side {
                y3 = y3.max(last_right_y + 11.0).min(H - 10.0);
                last_right_y = y3;
            } else {
                y3 = y3.min(last_left_y - 11.0).max(10.0);
                last_left_y = y3;
            }
            let (mut x4, _) = pos2xy(pos, len, cx, cy, r * 1.28);
            x3 = x3.max(0.0).min(W);
            y3 = y3.max(0.0).min(H);
            x4 = x4.max(0.0).min(W);
            let y4 = y3;

            doc = doc.add(
                Line::new()
                    .set("x1", x1)
                    .set("y1", y1)
                    .set("x2", x2)
                    .set("y2", y2)
                    .set("stroke", "#808080")
                    .set("stroke-width", 1),
            );
            doc = doc.add(
                Line::new()
                    .set("x1", x2)
                    .set("y1", y2)
                    .set("x2", x3)
                    .set("y2", y3)
                    .set("stroke", "#808080")
                    .set("stroke-width", 1),
            );
            doc = doc.add(
                Text::new(label)
                    .set("x", x4)
                    .set("y", y4)
                    .set("text-anchor", if right_side { "start" } else { "end" })
                    .set("dominant-baseline", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", font_color),
            );
        }
    }

    doc.to_string()
}

pub fn export_svg_pair(
    dna: &DNAsequence,
    display: &DisplaySettings,
) -> HashMap<&'static str, String> {
    let mut out = HashMap::new();
    out.insert("linear", export_linear_svg(dna, display));
    out.insert("circular", export_circular_svg(dna, display));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::DisplaySettings;
    use gb_io::{seq::Location, FeatureKind};
    #[cfg(feature = "snapshot-tests")]
    use std::fs;

    fn push_tfbs_feature(
        dna: &mut DNAsequence,
        label: &str,
        start: usize,
        end: usize,
        llr_bits: f64,
    ) {
        dna.features_mut().push(gb_io::seq::Feature {
            kind: FeatureKind::from("TFBS"),
            location: Location::simple_range(start as i64, end as i64),
            qualifiers: vec![
                ("label".into(), Some(label.to_string())),
                ("llr_bits".into(), Some(format!("{llr_bits:.3}"))),
                ("llr_quantile".into(), Some("1.000".to_string())),
                ("true_log_odds_bits".into(), Some("1.000".to_string())),
                ("true_log_odds_quantile".into(), Some("1.000".to_string())),
            ],
        });
    }

    fn push_vcf_track_feature(
        dna: &mut DNAsequence,
        label: &str,
        start: usize,
        end: usize,
        class: &str,
        filter: &str,
        qual: f64,
        info: &str,
    ) {
        dna.features_mut().push(gb_io::seq::Feature {
            kind: FeatureKind::from("track"),
            location: Location::simple_range(start as i64, end as i64),
            qualifiers: vec![
                ("label".into(), Some(label.to_string())),
                ("gentle_generated".into(), Some("genome_vcf_track".to_string())),
                ("vcf_variant_class".into(), Some(class.to_string())),
                ("vcf_filter".into(), Some(filter.to_string())),
                ("vcf_qual".into(), Some(format!("{qual:.3}"))),
                ("vcf_info".into(), Some(info.to_string())),
            ],
        });
    }

    #[test]
    #[cfg(feature = "snapshot-tests")]
    fn snapshot_linear_svg() {
        let mut dna = DNAsequence::from_sequence(&"ATGC".repeat(80)).unwrap();
        dna.update_computed_features();
        let svg = export_linear_svg(&dna, &DisplaySettings::default());
        let expected = include_str!("../tests/snapshots/linear/minimal.svg");
        assert_eq!(svg, expected);
    }

    #[test]
    #[cfg(feature = "snapshot-tests")]
    fn snapshot_circular_svg() {
        let mut dna = DNAsequence::from_sequence(&"ATGC".repeat(80)).unwrap();
        dna.set_circular(true);
        dna.update_computed_features();
        let svg = export_circular_svg(&dna, &DisplaySettings::default());
        let expected = include_str!("../tests/snapshots/circular/minimal.svg");
        assert_eq!(svg, expected);
    }

    #[test]
    fn tfbs_display_filter_applies_to_svg_export() {
        let mut dna_linear = DNAsequence::from_sequence(&"ATGC".repeat(80)).unwrap();
        dna_linear.update_computed_features();
        push_tfbs_feature(&mut dna_linear, "TFBS low", 10, 20, -2.0);
        push_tfbs_feature(&mut dna_linear, "TFBS high", 30, 40, 2.0);

        let mut display = DisplaySettings::default();
        display.show_tfbs = true;
        display.tfbs_display_use_llr_bits = true;
        display.tfbs_display_min_llr_bits = 0.0;
        display.tfbs_display_use_llr_quantile = false;
        display.tfbs_display_use_true_log_odds_bits = false;
        display.tfbs_display_use_true_log_odds_quantile = false;

        let linear_svg = export_linear_svg(&dna_linear, &display);
        assert!(linear_svg.contains("TFBS high"));
        assert!(!linear_svg.contains("TFBS low"));

        let mut dna_circular = dna_linear.clone();
        dna_circular.set_circular(true);
        let circular_svg = export_circular_svg(&dna_circular, &display);
        assert!(circular_svg.contains("TFBS high"));
        assert!(!circular_svg.contains("TFBS low"));
    }

    #[test]
    fn vcf_display_filter_applies_to_svg_export() {
        let mut dna_linear = DNAsequence::from_sequence(&"ATGC".repeat(80)).unwrap();
        dna_linear.update_computed_features();
        push_vcf_track_feature(
            &mut dna_linear,
            "VCF keep",
            10,
            20,
            "SNP",
            "PASS",
            40.0,
            "AC=1;AN=2",
        );
        push_vcf_track_feature(
            &mut dna_linear,
            "VCF drop",
            30,
            40,
            "DEL",
            "q10",
            10.0,
            "AN=2",
        );

        let mut display = DisplaySettings::default();
        display.vcf_display_show_snp = true;
        display.vcf_display_show_ins = false;
        display.vcf_display_show_del = false;
        display.vcf_display_show_sv = false;
        display.vcf_display_show_other = false;
        display.vcf_display_pass_only = true;
        display.vcf_display_use_min_qual = true;
        display.vcf_display_min_qual = 20.0;
        display.vcf_display_required_info_keys = vec!["AC".to_string()];

        let linear_svg = export_linear_svg(&dna_linear, &display);
        assert!(linear_svg.contains("VCF keep"));
        assert!(!linear_svg.contains("VCF drop"));
        assert!(linear_svg.contains("#e17f0f"));

        let mut dna_circular = dna_linear.clone();
        dna_circular.set_circular(true);
        let circular_svg = export_circular_svg(&dna_circular, &display);
        assert!(circular_svg.contains("VCF keep"));
        assert!(!circular_svg.contains("VCF drop"));
    }

    #[test]
    #[cfg(feature = "snapshot-tests")]
    #[ignore]
    fn write_snapshots() {
        let mut dna_linear = DNAsequence::from_sequence(&"ATGC".repeat(80)).unwrap();
        dna_linear.update_computed_features();
        let linear = export_linear_svg(&dna_linear, &DisplaySettings::default());
        fs::write("tests/snapshots/linear/minimal.svg", linear).unwrap();

        let mut dna_circular = DNAsequence::from_sequence(&"ATGC".repeat(80)).unwrap();
        dna_circular.set_circular(true);
        dna_circular.update_computed_features();
        let circular = export_circular_svg(&dna_circular, &DisplaySettings::default());
        fs::write("tests/snapshots/circular/minimal.svg", circular).unwrap();
    }
}
