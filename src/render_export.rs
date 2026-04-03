//! Shared export surfaces (SVG and snapshot pathways).

use crate::{
    dna_display::DnaDisplay,
    dna_sequence::DNAsequence,
    engine::DisplaySettings,
    feature_location::{collect_location_ranges_usize, feature_is_reverse},
    gc_contents::GcContents,
    restriction_enzyme::RestrictionEnzymeKey,
};
use gb_io::seq::Feature;
use std::collections::{HashMap, HashSet};
use svg::Document;
use svg::node::element::path::Data;
use svg::node::element::{Circle, Line, Path, Rectangle, Text};

const W: f32 = 1200.0;
const H: f32 = 700.0;
const RE_LABEL_BASE_OFFSET: f32 = 76.0;
const RE_LABEL_ROW_HEIGHT: f32 = 12.0;
const FEATURE_SIDE_MARGIN: f32 = 22.0;
const FEATURE_LANE_GAP: f32 = 13.0;
const FEATURE_BLOCK_HEIGHT: f32 = 10.0;
const REGULATORY_SIDE_MARGIN: f32 = 5.0;
const REGULATORY_LANE_GAP: f32 = 6.0;
const REGULATORY_BLOCK_HEIGHT: f32 = 6.0;
const REGULATORY_GROUP_GAP: f32 = 5.0;
const FEATURE_LANE_PADDING: f32 = 5.0;
const REGULATORY_LANE_PADDING: f32 = 2.0;
const LINEAR_SVG_HEADER_HEIGHT: f32 = 44.0;
const LINEAR_SVG_BOTTOM_PADDING: f32 = 28.0;
const LINEAR_SVG_MIN_HEIGHT: f32 = 180.0;
const SVG_TEXT_ASCENT: f32 = 10.0;
const SVG_TEXT_DESCENT: f32 = 4.0;
const VARIATION_MARKER_STROKE_WIDTH: f32 = 2.0;
const VARIATION_MARKER_OVERSHOOT_PX: f32 = 5.0;
const VARIATION_MARKER_RADIUS: f32 = 2.5;
const LINEAR_TSS_TICK_HALF_HEIGHT: f32 = 8.0;
const LINEAR_TSS_ARROW_LENGTH: f32 = 14.0;
const LINEAR_TSS_ARROW_HALF_HEIGHT: f32 = 4.5;
const REGULATORY_LABEL_DEDUP_DISTANCE_PX: f32 = 18.0;
const FEATURE_LABEL_DEDUP_DISTANCE_PX: f32 = 26.0;
const CIRCULAR_FEATURE_STROKE_WIDTH: f32 = 6.0;
const CIRCULAR_VARIATION_MARKER_STROKE_WIDTH: f32 = 2.5;
const CIRCULAR_VARIATION_MARKER_RADIUS: f32 = 4.0;
const CIRCULAR_VARIATION_MARKER_LABEL_OFFSET: f32 = 22.0;
const CIRCULAR_TSS_STEM_LENGTH: f32 = 10.0;
const CIRCULAR_TSS_ARROW_LENGTH: f32 = 12.0;
const CIRCULAR_TSS_ARROW_HALF_WIDTH: f32 = 4.0;
const CIRCULAR_TITLE_FONT_SIZE: usize = 18;
const CIRCULAR_SUBTITLE_FONT_SIZE: usize = 14;
const CIRCULAR_FEATURE_LABEL_FONT_SIZE: usize = 12;

#[derive(Clone, Debug)]
struct FeatureVm {
    from: usize,
    to: usize,
    label: String,
    color: &'static str,
    is_gene: bool,
    is_promoter: bool,
    is_reverse: bool,
    is_pointy: bool,
    is_regulatory: bool,
    is_variation: bool,
    has_transcription_direction: bool,
    is_fallback_label: bool,
}

#[derive(Clone, Debug)]
struct LinearSvgFeatureLayout {
    features: Vec<FeatureVm>,
    lane_top_by_idx: Vec<usize>,
    lane_bottom_by_idx: Vec<usize>,
    lane_regulatory_top_by_idx: Vec<usize>,
    top_regular_extent: f32,
    regulatory_group_gap: f32,
    top_extent: f32,
    bottom_extent: f32,
}

#[derive(Clone, Copy, Debug)]
struct LinearExportViewport {
    start_bp: usize,
    end_bp_exclusive: usize,
    span_bp: usize,
}

fn feature_name(feature: &Feature) -> (String, bool) {
    let kind = feature.kind.to_string().to_ascii_uppercase();
    if is_regulatory_feature(feature) {
        if let Some(reg_class) = feature_qualifier_text(feature, "regulatory_class") {
            let note = feature_qualifier_text(feature, "note");
            return (
                if let Some(note) = note {
                    format!("{reg_class}: {note}")
                } else {
                    reg_class
                },
                false,
            );
        }
    }
    if is_tfbs_feature(feature) {
        for key in ["bound_moiety", "standard_name", "name", "tf_id", "label"] {
            if let Some(value) = feature_qualifier_text(feature, key) {
                return (value, false);
            }
        }
    }
    if matches!(kind.as_str(), "GENE" | "MRNA" | "CDS") {
        for key in ["gene", "gene_name", "standard_name", "name", "locus_tag"] {
            if let Some(value) = feature_qualifier_text(feature, key) {
                return (value, false);
            }
        }
    }
    for k in [
        "label",
        "name",
        "standard_name",
        "gene",
        "gene_name",
        "protein_id",
        "product",
        "region_name",
        "bound_moiety",
    ] {
        if let Some(s) = feature.qualifier_values(k).next() {
            let label_text = s.to_string();
            if !label_text.trim().is_empty() {
                if matches!(kind.as_str(), "GENE" | "MRNA" | "CDS")
                    && looks_like_accession_label(&label_text)
                {
                    return (String::new(), false);
                }
                return (label_text, false);
            }
        }
    }
    (
        match feature.location.find_bounds() {
            Ok((from, to)) => format!("{from}..{to}"),
            Err(_) => String::new(),
        },
        true,
    )
}

fn feature_has_visible_export_label(feature: &FeatureVm) -> bool {
    !feature.label.trim().is_empty() && !feature.is_fallback_label
}

fn should_skip_nearby_repeated_label(
    placed_labels: &[(String, f32)],
    label: &str,
    x: f32,
    distance_px: f32,
) -> bool {
    placed_labels.iter().any(|(placed_label, placed_x)| {
        placed_label == label && (placed_x - x).abs() <= distance_px
    })
}

fn looks_like_accession_label(label: &str) -> bool {
    let trimmed = label.trim();
    if trimmed.is_empty() {
        return false;
    }
    let upper = trimmed.to_ascii_uppercase();
    [
        "ENST", "ENSG", "ENSP", "ENSE", "NM_", "NR_", "XM_", "XR_", "NP_", "XP_", "YP_", "WP_",
        "NC_", "NG_", "NT_", "NW_", "AP_",
    ]
    .iter()
    .any(|prefix| upper.starts_with(prefix))
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
    if is_regulatory_feature(feature) {
        let regulatory_class = feature_qualifier_text(feature, "regulatory_class")
            .unwrap_or_default()
            .to_ascii_lowercase();
        if regulatory_class.contains("silencer") || regulatory_class.contains("repressor") {
            return "#be3232";
        }
        if regulatory_class.contains("enhancer") || regulatory_class.contains("activator") {
            return "#2d9641";
        }
    }
    match feature.kind.to_string().to_ascii_uppercase().as_str() {
        "CDS" => "#cc1f1f",
        "GENE" => "#1f4fcc",
        "MRNA" => "#b4640a",
        "PROMOTER" => "#43aaa1",
        "VARIATION" => "#e17f0f",
        "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND" => "#238023",
        _ => "#6e6e6e",
    }
}

fn feature_pointy(feature: &Feature) -> bool {
    matches!(
        feature.kind.to_string().to_ascii_uppercase().as_str(),
        "CDS" | "GENE" | "MRNA" | "PROMOTER"
    )
}

fn has_transcription_direction(feature: &Feature) -> bool {
    matches!(
        feature.kind.to_string().to_ascii_uppercase().as_str(),
        "CDS" | "GENE" | "MRNA" | "PROMOTER"
    )
}

fn suppress_transcription_direction_marker(feature: &Feature) -> bool {
    for key in ["gentle_export_tss_marker", "gentle_hide_tss_marker"] {
        for value in feature.qualifier_values(key) {
            let normalized = value.trim().to_ascii_lowercase();
            if matches!(
                normalized.as_str(),
                "none" | "hide" | "hidden" | "off" | "false"
            ) {
                return true;
            }
        }
    }
    false
}

fn is_tfbs_feature(feature: &Feature) -> bool {
    matches!(
        feature.kind.to_string().to_ascii_uppercase().as_str(),
        "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND"
    )
}

fn is_variation_feature(feature: &Feature) -> bool {
    feature.kind.to_string().to_ascii_uppercase() == "VARIATION"
}

fn is_promoter_feature(feature: &Feature) -> bool {
    feature.kind.to_string().to_ascii_uppercase() == "PROMOTER"
}

fn is_track_feature(feature: &Feature) -> bool {
    feature.qualifier_values("gentle_generated").any(|value| {
        matches!(
            value.trim().to_ascii_lowercase().as_str(),
            "genome_bed_track" | "genome_bigwig_track" | "genome_vcf_track" | "blast_hit_track"
        )
    })
}

fn is_vcf_track_feature(feature: &Feature) -> bool {
    feature
        .qualifier_values("gentle_generated")
        .any(|value| value.trim().eq_ignore_ascii_case("genome_vcf_track"))
}

fn has_regulatory_hint(feature: &Feature) -> bool {
    for key in [
        "regulatory_class",
        "regulation",
        "function",
        "note",
        "label",
    ] {
        for value in feature.qualifier_values(key) {
            let lower = value.to_ascii_lowercase();
            if lower.contains("regulatory")
                || lower.contains("enhancer")
                || lower.contains("promoter")
                || lower.contains("silencer")
                || lower.contains("insulator")
                || lower.contains("atac")
                || lower.contains("chip")
            {
                return true;
            }
        }
    }
    false
}

fn is_regulatory_feature(feature: &Feature) -> bool {
    let kind = feature.kind.to_string().to_ascii_uppercase();
    is_tfbs_feature(feature)
        || kind.contains("REGULATORY")
        || (kind == "MISC_FEATURE" && has_regulatory_hint(feature))
        || is_track_feature(feature)
}

fn feature_qualifier_f64(feature: &Feature, key: &str) -> Option<f64> {
    feature
        .qualifier_values(key)
        .next()
        .and_then(|v| v.trim().parse::<f64>().ok())
}

fn feature_qualifier_text(feature: &Feature, key: &str) -> Option<String> {
    feature
        .qualifier_values(key)
        .map(|value| value.split_whitespace().collect::<Vec<_>>().join(" "))
        .map(|value| value.trim().to_string())
        .find(|value| !value.is_empty())
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
        let qual = feature_qualifier_f64(feature, "vcf_qual")
            .or_else(|| feature_qualifier_f64(feature, "score"));
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

fn feature_max_view_span_bp(feature: &Feature, regulatory_max_view_span_bp: usize) -> usize {
    let kind = feature.kind.to_string().to_ascii_uppercase();
    if kind == "SOURCE" {
        return 0;
    }
    if is_regulatory_feature(feature) {
        return regulatory_max_view_span_bp;
    }
    match kind.as_str() {
        "CDS" | "GENE" | "MRNA" => 2_000_000,
        "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND" => 200_000,
        "TRACK" => 250_000,
        "RESTRICTION_SITE" => 100_000,
        _ => 500_000,
    }
}

fn normalize_linear_export_viewport(
    dna: &DNAsequence,
    display: &DisplaySettings,
) -> LinearExportViewport {
    let seq_len = dna.len();
    if seq_len == 0 || dna.is_circular() {
        return LinearExportViewport {
            start_bp: 0,
            end_bp_exclusive: seq_len,
            span_bp: seq_len.max(1),
        };
    }
    let requested_span = display.linear_view_span_bp;
    let span_bp = if requested_span == 0 || requested_span > seq_len {
        seq_len
    } else {
        requested_span
    };
    let max_start = seq_len.saturating_sub(span_bp);
    let start_bp = display.linear_view_start_bp.min(max_start);
    let end_bp_exclusive = start_bp.saturating_add(span_bp).min(seq_len);
    LinearExportViewport {
        start_bp,
        end_bp_exclusive,
        span_bp: end_bp_exclusive.saturating_sub(start_bp).max(1),
    }
}

fn absolute_bp_to_view_x(
    absolute_bp: usize,
    viewport: LinearExportViewport,
    left: f32,
    right: f32,
) -> f32 {
    let local_bp = absolute_bp.saturating_sub(viewport.start_bp);
    bp_to_x(local_bp, viewport.span_bp, left, right)
}

fn clip_linear_bounds_to_viewport(
    from: usize,
    to: usize,
    sequence_length: usize,
    viewport: LinearExportViewport,
) -> Option<(usize, usize)> {
    if sequence_length == 0 || viewport.end_bp_exclusive <= viewport.start_bp {
        return None;
    }
    let mut from = from;
    let mut to = to;
    if to < from {
        std::mem::swap(&mut from, &mut to);
    }
    if from >= sequence_length {
        return None;
    }
    let to = to.min(sequence_length.saturating_sub(1));
    if to < from {
        return None;
    }
    let clip_from = from.max(viewport.start_bp);
    let clip_to = to.min(viewport.end_bp_exclusive.saturating_sub(1));
    (clip_to >= clip_from).then_some((clip_from, clip_to))
}

fn feature_bounds_in_viewport(
    feature: &Feature,
    sequence_length: usize,
    viewport: LinearExportViewport,
) -> Option<(usize, usize)> {
    let mut ranges = Vec::new();
    collect_location_ranges_usize(&feature.location, &mut ranges);
    if ranges.is_empty() {
        if let Ok((raw_from, raw_to)) = feature.location.find_bounds() {
            if raw_from >= 0 && raw_to >= 0 {
                ranges.push((raw_from as usize, raw_to as usize));
            }
        }
    }
    let mut clipped_start: Option<usize> = None;
    let mut clipped_end: usize = 0;
    for (from, to) in ranges {
        let Some((clip_from, clip_to)) =
            clip_linear_bounds_to_viewport(from, to, sequence_length, viewport)
        else {
            continue;
        };
        clipped_start = Some(clipped_start.map_or(clip_from, |value| value.min(clip_from)));
        clipped_end = clipped_end.max(clip_to);
    }
    clipped_start.map(|start| (start, clipped_end))
}

fn orf_bounds_in_viewport(
    from: usize,
    to: usize,
    sequence_length: usize,
    viewport: LinearExportViewport,
) -> Option<(usize, usize)> {
    clip_linear_bounds_to_viewport(from, to, sequence_length, viewport)
}

fn collect_features(
    dna: &DNAsequence,
    display: &DisplaySettings,
    view_span_bp: usize,
    viewport: LinearExportViewport,
) -> Vec<FeatureVm> {
    let mut ret = Vec::new();
    let sequence_length = dna.len();
    for feature in dna.features() {
        if feature.kind.to_string().to_ascii_uppercase() == "SOURCE" {
            continue;
        }
        let max_view_span =
            feature_max_view_span_bp(feature, display.regulatory_feature_max_view_span_bp);
        if max_view_span > 0 && view_span_bp > max_view_span {
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
        let Some((from, to)) = feature_bounds_in_viewport(feature, sequence_length, viewport)
        else {
            continue;
        };
        let is_variation = is_variation_feature(feature);
        let (label, is_fallback_label) = feature_name(feature);
        ret.push(FeatureVm {
            from,
            to,
            label,
            color: feature_color(feature),
            is_gene: kind == "GENE",
            is_promoter: is_promoter_feature(feature),
            is_reverse: feature_is_reverse(feature),
            is_pointy: feature_pointy(feature),
            is_regulatory: is_regulatory_feature(feature),
            is_variation,
            has_transcription_direction: has_transcription_direction(feature)
                && !suppress_transcription_direction_marker(feature),
            is_fallback_label,
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

fn linear_transcription_direction_glyph(
    feature: &FeatureVm,
    viewport: LinearExportViewport,
    left: f32,
    right: f32,
    y: f32,
) -> (Path, Path) {
    let tss_pos = if feature.is_promoter {
        if feature.is_reverse {
            feature.from
        } else {
            feature.to
        }
    } else {
        if feature.is_reverse {
            feature.to
        } else {
            feature.from
        }
    };
    let x = absolute_bp_to_view_x(tss_pos, viewport, left, right);
    let stem_y = if feature.is_reverse {
        y + LINEAR_TSS_TICK_HALF_HEIGHT
    } else {
        y - LINEAR_TSS_TICK_HALF_HEIGHT
    };
    let hook_y = if feature.is_reverse {
        stem_y + LINEAR_TSS_ARROW_HALF_HEIGHT * 1.3
    } else {
        stem_y - LINEAR_TSS_ARROW_HALF_HEIGHT * 1.3
    };
    let elbow_x = if feature.is_reverse {
        x - LINEAR_TSS_ARROW_LENGTH * 0.22
    } else {
        x + LINEAR_TSS_ARROW_LENGTH * 0.22
    };
    let tip_x = if feature.is_reverse {
        x - LINEAR_TSS_ARROW_LENGTH
    } else {
        x + LINEAR_TSS_ARROW_LENGTH
    };
    let shaft = Data::new()
        .move_to((x, y))
        .line_to((x, stem_y))
        .quadratic_curve_to((x, hook_y, elbow_x, hook_y))
        .line_to((tip_x, hook_y));
    let tick = Path::new()
        .set("d", shaft)
        .set("fill", "none")
        .set("stroke", feature.color)
        .set("stroke-width", 2)
        .set("stroke-linecap", "round")
        .set("stroke-linejoin", "round")
        .set("data-gentle-role", "linear-transcription-start-tick");
    let (head_tip_x, base_x) = if feature.is_reverse {
        (
            tip_x - LINEAR_TSS_ARROW_LENGTH * 0.14,
            tip_x + LINEAR_TSS_ARROW_LENGTH * 0.18,
        )
    } else {
        (
            tip_x + LINEAR_TSS_ARROW_LENGTH * 0.14,
            tip_x - LINEAR_TSS_ARROW_LENGTH * 0.18,
        )
    };
    let tri = Data::new()
        .move_to((head_tip_x, hook_y))
        .line_to((base_x, hook_y - LINEAR_TSS_ARROW_HALF_HEIGHT))
        .line_to((base_x, hook_y + LINEAR_TSS_ARROW_HALF_HEIGHT))
        .close();
    let arrow = Path::new()
        .set("d", tri)
        .set("fill", feature.color)
        .set("stroke", "none")
        .set("data-gentle-role", "linear-transcription-start-arrow");
    (tick, arrow)
}

fn estimate_text_width(label: &str) -> f32 {
    (label.chars().count().max(1) as f32) * 6.5
}

fn truncate_label_to_chars(label: &str, max_chars: usize) -> String {
    if max_chars == 0 {
        return String::new();
    }
    let total = label.chars().count();
    if total <= max_chars {
        return label.to_string();
    }
    if max_chars == 1 {
        return "…".to_string();
    }
    let mut out = label
        .chars()
        .take(max_chars.saturating_sub(1))
        .collect::<String>();
    out.push('…');
    out
}

fn re_color(key: &RestrictionEnzymeKey) -> &'static str {
    match key.number_of_cuts() {
        1 => "#8b0000",
        2 => "#00008b",
        3 => "#006400",
        _ => "#000000",
    }
}

fn compute_linear_svg_feature_layout(
    dna: &DNAsequence,
    display: &DisplaySettings,
    viewport: LinearExportViewport,
    left: f32,
    right: f32,
) -> LinearSvgFeatureLayout {
    let regulatory_tracks_near_baseline = display.regulatory_tracks_near_baseline;
    let features = collect_features(dna, display, viewport.span_bp, viewport);
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
            if !f.is_reverse && !f.is_regulatory && !f.is_variation {
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
        fa.from.cmp(&fb.from).then(span_b.cmp(&span_a))
    });
    for idx in top_order {
        let f = &features[idx];
        let x1 = absolute_bp_to_view_x(f.from, viewport, left, right);
        let x2 = absolute_bp_to_view_x(f.to, viewport, left, right).max(x1 + 1.0);
        lane_top_by_idx[idx] = lane_allocate(&mut top_lane_ends, x1, x2, FEATURE_LANE_PADDING);
    }

    let mut bottom_order: Vec<usize> = features
        .iter()
        .enumerate()
        .filter_map(|(idx, f)| {
            if f.is_reverse && !f.is_regulatory && !f.is_variation {
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
        fa.from.cmp(&fb.from).then(span_b.cmp(&span_a))
    });
    for idx in bottom_order {
        let f = &features[idx];
        let x1 = absolute_bp_to_view_x(f.from, viewport, left, right);
        let x2 = absolute_bp_to_view_x(f.to, viewport, left, right).max(x1 + 1.0);
        lane_bottom_by_idx[idx] =
            lane_allocate(&mut bottom_lane_ends, x1, x2, FEATURE_LANE_PADDING);
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
        fa.from.cmp(&fb.from).then(span_b.cmp(&span_a))
    });
    for idx in regulatory_top_order {
        let f = &features[idx];
        if regulatory_tracks_near_baseline {
            lane_regulatory_top_by_idx[idx] = 0;
        } else {
            let x1 = absolute_bp_to_view_x(f.from, viewport, left, right);
            let x2 = absolute_bp_to_view_x(f.to, viewport, left, right).max(x1 + 1.0);
            lane_regulatory_top_by_idx[idx] = lane_allocate(
                &mut regulatory_top_lane_ends,
                x1,
                x2,
                REGULATORY_LANE_PADDING,
            );
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
    let top_regular_margin = if regulatory_tracks_near_baseline {
        FEATURE_SIDE_MARGIN.max(
            REGULATORY_SIDE_MARGIN
                + REGULATORY_BLOCK_HEIGHT * 0.5
                + FEATURE_BLOCK_HEIGHT * 0.5
                + 1.0,
        )
    } else {
        FEATURE_SIDE_MARGIN
    };

    let mut top_extent = 0.0f32;
    let mut bottom_extent = 0.0f32;
    for (idx, f) in features.iter().enumerate() {
        if f.is_variation {
            if !f.label.trim().is_empty() {
                top_extent = top_extent.max(12.0 + SVG_TEXT_ASCENT);
            }
        } else if f.is_regulatory {
            let lane = lane_regulatory_top_by_idx[idx] as f32;
            let center_offset = if regulatory_tracks_near_baseline {
                REGULATORY_SIDE_MARGIN
            } else {
                top_regular_extent
                    + regulatory_group_gap
                    + REGULATORY_SIDE_MARGIN
                    + lane * REGULATORY_LANE_GAP
            };
            top_extent = top_extent.max(center_offset + REGULATORY_BLOCK_HEIGHT * 0.5);
            if !f.label.trim().is_empty() {
                top_extent = top_extent
                    .max(center_offset + REGULATORY_BLOCK_HEIGHT * 0.5 + 2.0 + SVG_TEXT_ASCENT);
            }
        } else if f.is_reverse {
            let lane = lane_bottom_by_idx[idx] as f32;
            let center_offset = FEATURE_SIDE_MARGIN + lane * FEATURE_LANE_GAP;
            bottom_extent = bottom_extent.max(center_offset + FEATURE_BLOCK_HEIGHT * 0.5);
            if !f.label.trim().is_empty() && !f.is_gene {
                bottom_extent = bottom_extent.max(center_offset + 16.0 + SVG_TEXT_DESCENT);
            }
        } else {
            let lane = lane_top_by_idx[idx] as f32;
            let center_offset = top_regular_margin + lane * FEATURE_LANE_GAP;
            top_extent = top_extent.max(center_offset + FEATURE_BLOCK_HEIGHT * 0.5);
            if !f.label.trim().is_empty() && !f.is_gene {
                top_extent = top_extent.max(center_offset + 10.0 + SVG_TEXT_ASCENT);
            }
        }
    }

    LinearSvgFeatureLayout {
        features,
        lane_top_by_idx,
        lane_bottom_by_idx,
        lane_regulatory_top_by_idx,
        top_regular_extent,
        regulatory_group_gap,
        top_extent,
        bottom_extent,
    }
}

fn compute_linear_svg_re_label_extents(
    dna: &DNAsequence,
    display: &DisplaySettings,
    viewport: LinearExportViewport,
    left: f32,
    right: f32,
) -> (f32, f32) {
    if !display.show_restriction_enzymes {
        return (0.0, 0.0);
    }
    let mut keys: Vec<RestrictionEnzymeKey> = dna
        .restriction_enzyme_groups()
        .iter()
        .filter_map(|(key, names)| {
            DnaDisplay::restriction_group_matches_mode(
                display.restriction_enzyme_display_mode,
                &display.preferred_restriction_enzymes,
                key,
                names,
            )
            .then_some(key.clone())
        })
        .collect();
    keys.sort();
    let mut top_label_lanes: Vec<f32> = vec![];
    let mut bottom_label_lanes: Vec<f32> = vec![];
    for (i, key) in keys.iter().enumerate() {
        let Some(pos) = usize::try_from(key.pos().max(0)).ok() else {
            continue;
        };
        let Some(mate_pos) = usize::try_from(key.mate_pos().max(0)).ok() else {
            continue;
        };
        let in_view = (pos >= viewport.start_bp && pos < viewport.end_bp_exclusive)
            || (mate_pos >= viewport.start_bp && mate_pos < viewport.end_bp_exclusive);
        if !in_view {
            continue;
        }
        let x = absolute_bp_to_view_x(pos, viewport, left, right);
        if let Some(names) = dna.restriction_enzyme_groups().get(key) {
            let label = names.join(",");
            let label_width = estimate_text_width(&label);
            let label_left = x - label_width / 2.0;
            let label_right = x + label_width / 2.0;
            if i % 2 == 0 {
                let _ = lane_allocate(&mut top_label_lanes, label_left, label_right, 6.0);
            } else {
                let _ = lane_allocate(&mut bottom_label_lanes, label_left, label_right, 6.0);
            }
        }
    }
    let top_extent = if top_label_lanes.is_empty() {
        0.0
    } else {
        RE_LABEL_BASE_OFFSET
            + (top_label_lanes.len().saturating_sub(1) as f32) * RE_LABEL_ROW_HEIGHT
            + SVG_TEXT_ASCENT
    };
    let bottom_extent = if bottom_label_lanes.is_empty() {
        0.0
    } else {
        RE_LABEL_BASE_OFFSET
            + (bottom_label_lanes.len().saturating_sub(1) as f32) * RE_LABEL_ROW_HEIGHT
            + SVG_TEXT_DESCENT
    };
    (top_extent, bottom_extent)
}

pub fn export_linear_svg(dna: &DNAsequence, display: &DisplaySettings) -> String {
    let len = dna.len();
    let left = 60.0;
    let right = W - 60.0;
    let viewport = normalize_linear_export_viewport(dna, display);
    let feature_layout = display
        .show_features
        .then(|| compute_linear_svg_feature_layout(dna, display, viewport, left, right));
    let (re_top_extent, re_bottom_extent) =
        compute_linear_svg_re_label_extents(dna, display, viewport, left, right);
    let mut top_extent = 0.0f32;
    let mut bottom_extent = 0.0f32;
    if display.show_gc_contents {
        top_extent = top_extent.max(8.0);
    }
    if display.show_methylation_sites {
        top_extent = top_extent.max(12.0);
    }
    if display.show_open_reading_frames {
        let mut top_orf_extent = 0.0f32;
        let mut bottom_orf_extent = 0.0f32;
        for orf in dna.open_reading_frames() {
            let level = (orf.frame().unsigned_abs()) as f32;
            let extent = 12.5 + 7.0 * level;
            if orf.is_reverse() {
                bottom_orf_extent = bottom_orf_extent.max(extent);
            } else {
                top_orf_extent = top_orf_extent.max(extent);
            }
        }
        top_extent = top_extent.max(top_orf_extent);
        bottom_extent = bottom_extent.max(bottom_orf_extent);
    }
    if let Some(layout) = &feature_layout {
        top_extent = top_extent.max(layout.top_extent);
        bottom_extent = bottom_extent.max(layout.bottom_extent);
    }
    top_extent = top_extent.max(re_top_extent);
    bottom_extent = bottom_extent.max(re_bottom_extent);

    let baseline = LINEAR_SVG_HEADER_HEIGHT + top_extent;
    let canvas_height = (baseline + bottom_extent + LINEAR_SVG_BOTTOM_PADDING)
        .max(LINEAR_SVG_MIN_HEIGHT)
        .ceil();

    let mut doc = Document::new()
        .set("viewBox", (0, 0, W, canvas_height))
        .set("width", W)
        .set("height", canvas_height)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", W)
                .set("height", canvas_height)
                .set("fill", "#ffffff"),
        );
    let mut labels: Vec<Text> = vec![];

    if display.show_gc_contents {
        let gc_contents = GcContents::new_from_sequence_with_bin_size(
            dna.forward_bytes(),
            display.gc_content_bin_size_bp,
        );
        for region in gc_contents.regions() {
            let Some((from, to)) =
                clip_linear_bounds_to_viewport(region.from(), region.to(), len, viewport)
            else {
                continue;
            };
            let x1 = absolute_bp_to_view_x(from, viewport, left, right);
            let x2 = absolute_bp_to_view_x(to, viewport, left, right).max(x1 + 1.0);
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
            if *site < viewport.start_bp || *site >= viewport.end_bp_exclusive {
                continue;
            }
            let x = absolute_bp_to_view_x(*site, viewport, left, right);
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
            let Some((view_start, view_end)) = orf_bounds_in_viewport(start, end, len, viewport)
            else {
                continue;
            };
            let x1 = absolute_bp_to_view_x(view_start, viewport, left, right);
            let x2 = absolute_bp_to_view_x(view_end, viewport, left, right).max(x1 + 1.0);
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
        if let Some(layout) = &feature_layout {
            let regulatory_tracks_near_baseline = display.regulatory_tracks_near_baseline;
            let features = &layout.features;
            let lane_top_by_idx = &layout.lane_top_by_idx;
            let lane_bottom_by_idx = &layout.lane_bottom_by_idx;
            let lane_regulatory_top_by_idx = &layout.lane_regulatory_top_by_idx;
            let top_regular_extent = layout.top_regular_extent;
            let regulatory_group_gap = layout.regulatory_group_gap;
            let mut placed_regulatory_labels: Vec<(String, f32)> = Vec::new();
            let mut placed_top_feature_labels: Vec<(String, f32)> = Vec::new();
            let mut placed_bottom_feature_labels: Vec<(String, f32)> = Vec::new();

            for (idx, f) in features.iter().enumerate() {
                let x1 = absolute_bp_to_view_x(f.from, viewport, left, right);
                let x2 = absolute_bp_to_view_x(f.to, viewport, left, right).max(x1 + 1.0);
                if f.is_variation {
                    let x = (x1 + x2) * 0.5;
                    doc = doc.add(
                        Line::new()
                            .set("x1", x)
                            .set("y1", baseline - VARIATION_MARKER_OVERSHOOT_PX)
                            .set("x2", x)
                            .set("y2", baseline + VARIATION_MARKER_OVERSHOOT_PX)
                            .set("stroke", f.color)
                            .set("stroke-width", VARIATION_MARKER_STROKE_WIDTH),
                    );
                    doc = doc.add(
                        Circle::new()
                            .set("cx", x)
                            .set("cy", baseline)
                            .set("r", VARIATION_MARKER_RADIUS)
                            .set("fill", "#ffffff")
                            .set("stroke", f.color)
                            .set("stroke-width", 1.5),
                    );
                    if feature_has_visible_export_label(f) {
                        labels.push(
                            Text::new(f.label.clone())
                                .set("x", x)
                                .set("y", baseline - 10.0)
                                .set("text-anchor", "middle")
                                .set("font-family", "monospace")
                                .set("font-size", 9)
                                .set("fill", "#111111"),
                        );
                    }
                    continue;
                }
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
                    let required_margin = if regulatory_tracks_near_baseline {
                        REGULATORY_SIDE_MARGIN
                            + REGULATORY_BLOCK_HEIGHT * 0.5
                            + FEATURE_BLOCK_HEIGHT * 0.5
                            + 1.0
                    } else {
                        FEATURE_SIDE_MARGIN
                    };
                    let y = baseline
                        - required_margin.max(FEATURE_SIDE_MARGIN)
                        - lane as f32 * FEATURE_LANE_GAP;
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

                if f.has_transcription_direction {
                    let (tick, arrow) =
                        linear_transcription_direction_glyph(f, viewport, left, right, y);
                    doc = doc.add(tick);
                    doc = doc.add(arrow);
                }

                if !feature_has_visible_export_label(f) {
                    continue;
                }
                if f.is_gene {
                    let inline_chars = (((x2 - x1) - 4.0).max(0.0) / 6.5).floor() as usize;
                    let inline_label = truncate_label_to_chars(f.label.trim(), inline_chars.max(1));
                    labels.push(
                        Text::new(inline_label)
                            .set("x", (x1 + x2) * 0.5)
                            .set("y", y + 3.0)
                            .set("text-anchor", "middle")
                            .set("font-family", "monospace")
                            .set("font-size", 10)
                            .set("fill", "#f8f8f8"),
                    );
                    continue;
                }
                let (text_y, anchor) = if f.is_regulatory {
                    (y - half_height - 2.0, "start")
                } else if f.is_reverse {
                    (y + 16.0, "start")
                } else {
                    (y - 10.0, "start")
                };
                if f.is_regulatory
                    && should_skip_nearby_repeated_label(
                        &placed_regulatory_labels,
                        &f.label,
                        x1,
                        REGULATORY_LABEL_DEDUP_DISTANCE_PX,
                    )
                {
                    continue;
                }
                if !f.is_regulatory
                    && !f.is_gene
                    && should_skip_nearby_repeated_label(
                        if f.is_reverse {
                            &placed_bottom_feature_labels
                        } else {
                            &placed_top_feature_labels
                        },
                        &f.label,
                        x1,
                        FEATURE_LABEL_DEDUP_DISTANCE_PX,
                    )
                {
                    continue;
                }
                labels.push(
                    Text::new(f.label.clone())
                        .set("x", x1)
                        .set("y", text_y)
                        .set("text-anchor", anchor)
                        .set("font-family", "monospace")
                        .set("font-size", 10)
                        .set("fill", "#111111"),
                );
                if f.is_regulatory {
                    placed_regulatory_labels.push((f.label.clone(), x1));
                } else if !f.is_gene {
                    if f.is_reverse {
                        placed_bottom_feature_labels.push((f.label.clone(), x1));
                    } else {
                        placed_top_feature_labels.push((f.label.clone(), x1));
                    }
                }
            }
        }
    }

    if display.show_restriction_enzymes {
        let mut keys: Vec<RestrictionEnzymeKey> = dna
            .restriction_enzyme_groups()
            .iter()
            .filter_map(|(key, names)| {
                DnaDisplay::restriction_group_matches_mode(
                    display.restriction_enzyme_display_mode,
                    &display.preferred_restriction_enzymes,
                    key,
                    names,
                )
                .then_some(key.clone())
            })
            .collect();
        keys.sort();
        let mut top_label_lanes: Vec<f32> = vec![];
        let mut bottom_label_lanes: Vec<f32> = vec![];

        for (i, key) in keys.iter().enumerate() {
            let Some(pos) = usize::try_from(key.pos().max(0)).ok() else {
                continue;
            };
            let Some(mate_pos) = usize::try_from(key.mate_pos().max(0)).ok() else {
                continue;
            };
            let in_view = (pos >= viewport.start_bp && pos < viewport.end_bp_exclusive)
                || (mate_pos >= viewport.start_bp && mate_pos < viewport.end_bp_exclusive);
            if !in_view {
                continue;
            }
            let x = absolute_bp_to_view_x(pos, viewport, left, right);
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
        Text::new(
            if viewport.start_bp == 0 && viewport.end_bp_exclusive >= len {
                format!("{} bp", len)
            } else {
                format!(
                    "{}..{} ({} bp view of {} bp)",
                    viewport.start_bp.saturating_add(1),
                    viewport.end_bp_exclusive,
                    viewport.span_bp,
                    len
                )
            },
        )
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

fn circular_feature_band_radius(base_radius: f32, feature: &FeatureVm, seq_len: usize) -> f32 {
    let span = feature.to.saturating_sub(feature.from).saturating_add(1);
    let length_fraction = if seq_len == 0 {
        0.0
    } else {
        (span as f32 / seq_len as f32).clamp(0.0, 1.0)
    };
    let offset = (0.14 * (1.0 - 0.6 * length_fraction)).clamp(0.05, 0.14);
    let band = if feature.is_reverse {
        1.0 - offset
    } else {
        1.0 + offset
    };
    base_radius * band
}

fn circular_tangent_unit_vector(
    pos: usize,
    is_reverse: bool,
    len: usize,
    cx: f32,
    cy: f32,
    r: f32,
) -> (f32, f32) {
    if len == 0 {
        return (1.0, 0.0);
    }
    let next_pos = if is_reverse {
        pos.checked_sub(1).unwrap_or(len.saturating_sub(1))
    } else {
        (pos + 1) % len
    };
    let (x1, y1) = pos2xy(pos, len, cx, cy, r);
    let (x2, y2) = pos2xy(next_pos, len, cx, cy, r);
    let dx = x2 - x1;
    let dy = y2 - y1;
    let norm = (dx * dx + dy * dy).sqrt();
    if norm <= f32::EPSILON {
        (1.0, 0.0)
    } else {
        (dx / norm, dy / norm)
    }
}

fn circular_radial_unit_vector(pos: usize, len: usize, cx: f32, cy: f32, r: f32) -> (f32, f32) {
    let (x, y) = pos2xy(pos, len, cx, cy, r);
    let dx = x - cx;
    let dy = y - cy;
    let norm = (dx * dx + dy * dy).sqrt();
    if norm <= f32::EPSILON {
        (0.0, -1.0)
    } else {
        (dx / norm, dy / norm)
    }
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
    let r = H.min(W) * 0.33;

    let mut doc = Document::new()
        .set("viewBox", (0, 0, W, H))
        .set("width", W)
        .set("height", H);

    doc = doc.add(
        Text::new(
            dna.name()
                .clone()
                .unwrap_or_else(|| "<no name>".to_string()),
        )
        .set("x", 20)
        .set("y", 34)
        .set("font-family", "monospace")
        .set("font-size", CIRCULAR_TITLE_FONT_SIZE)
        .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new(format!("{} bp", len))
            .set("x", 20)
            .set("y", 54)
            .set("font-family", "monospace")
            .set("font-size", CIRCULAR_SUBTITLE_FONT_SIZE)
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
        let gc_contents = GcContents::new_from_sequence_with_bin_size(
            dna.forward_bytes(),
            display.gc_content_bin_size_bp,
        );
        for region in gc_contents.regions() {
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
        for f in collect_features(
            dna,
            display,
            len,
            normalize_linear_export_viewport(dna, display),
        ) {
            let band_radius = circular_feature_band_radius(r, &f, len);
            let mid = (f.from + f.to) / 2;
            if f.is_variation {
                let marker_pos = mid.min(len.saturating_sub(1));
                let (inner_x, inner_y) = pos2xy(marker_pos, len, cx, cy, r * 0.94);
                let (outer_x, outer_y) = pos2xy(marker_pos, len, cx, cy, r * 1.08);
                let (label_x, label_y) = pos2xy(
                    marker_pos,
                    len,
                    cx,
                    cy,
                    r * 1.08 + CIRCULAR_VARIATION_MARKER_LABEL_OFFSET,
                );
                doc = doc.add(
                    Line::new()
                        .set("x1", inner_x)
                        .set("y1", inner_y)
                        .set("x2", outer_x)
                        .set("y2", outer_y)
                        .set("stroke", f.color)
                        .set("stroke-width", CIRCULAR_VARIATION_MARKER_STROKE_WIDTH)
                        .set("data-gentle-role", "variation-marker-line"),
                );
                doc = doc.add(
                    Circle::new()
                        .set("cx", outer_x)
                        .set("cy", outer_y)
                        .set("r", CIRCULAR_VARIATION_MARKER_RADIUS)
                        .set("fill", "#ffffff")
                        .set("stroke", f.color)
                        .set("stroke-width", 2)
                        .set("data-gentle-role", "variation-marker-dot"),
                );
                if feature_has_visible_export_label(&f) {
                    doc = doc.add(
                        Text::new(f.label)
                            .set("x", label_x)
                            .set("y", label_y)
                            .set("text-anchor", "middle")
                            .set("font-family", "monospace")
                            .set("font-size", CIRCULAR_FEATURE_LABEL_FONT_SIZE)
                            .set("fill", "#111111"),
                    );
                }
                continue;
            }
            if let Some(path_d) = circular_arc_path(f.from, f.to, len, cx, cy, band_radius) {
                doc = doc.add(
                    Path::new()
                        .set("d", path_d)
                        .set("fill", "none")
                        .set("stroke", f.color)
                        .set("stroke-width", CIRCULAR_FEATURE_STROKE_WIDTH),
                );
            }
            if f.has_transcription_direction {
                let tss_pos = if f.is_promoter {
                    if f.is_reverse { f.from } else { f.to }
                } else if f.is_reverse {
                    f.to
                } else {
                    f.from
                }
                .min(len.saturating_sub(1));
                let (tick_dx, tick_dy) = circular_radial_unit_vector(tss_pos, len, cx, cy, r);
                let (tick_x, tick_y) = pos2xy(tss_pos, len, cx, cy, band_radius);
                let (dir_x, dir_y) =
                    circular_tangent_unit_vector(tss_pos, f.is_reverse, len, cx, cy, band_radius);
                let stem_x = tick_x + tick_dx * CIRCULAR_TSS_STEM_LENGTH;
                let stem_y = tick_y + tick_dy * CIRCULAR_TSS_STEM_LENGTH;
                let shaft_tip_x = stem_x + dir_x * (CIRCULAR_TSS_ARROW_LENGTH * 0.55);
                let shaft_tip_y = stem_y + dir_y * (CIRCULAR_TSS_ARROW_LENGTH * 0.55);
                let shaft = Data::new()
                    .move_to((tick_x, tick_y))
                    .line_to((stem_x, stem_y))
                    .line_to((shaft_tip_x, shaft_tip_y));
                doc = doc.add(
                    Path::new()
                        .set("d", shaft)
                        .set("fill", "none")
                        .set("stroke", f.color)
                        .set("stroke-width", 2)
                        .set("stroke-linecap", "round")
                        .set("stroke-linejoin", "round")
                        .set("data-gentle-role", "transcription-start-arrow-shaft"),
                );
                let tip_x = stem_x + dir_x * (CIRCULAR_TSS_ARROW_LENGTH * 0.7);
                let tip_y = stem_y + dir_y * (CIRCULAR_TSS_ARROW_LENGTH * 0.7);
                let back_x = stem_x + dir_x * (CIRCULAR_TSS_ARROW_LENGTH * 0.1);
                let back_y = stem_y + dir_y * (CIRCULAR_TSS_ARROW_LENGTH * 0.1);
                let tri = Data::new()
                    .move_to((tip_x, tip_y))
                    .line_to((
                        back_x + tick_dx * CIRCULAR_TSS_ARROW_HALF_WIDTH,
                        back_y + tick_dy * CIRCULAR_TSS_ARROW_HALF_WIDTH,
                    ))
                    .line_to((
                        back_x - tick_dx * CIRCULAR_TSS_ARROW_HALF_WIDTH,
                        back_y - tick_dy * CIRCULAR_TSS_ARROW_HALF_WIDTH,
                    ))
                    .close();
                doc = doc.add(
                    Path::new()
                        .set("d", tri)
                        .set("fill", f.color)
                        .set("stroke", "none")
                        .set("data-gentle-role", "transcription-start-arrow"),
                );
            }
            let label_radius = if f.is_reverse {
                band_radius - 22.0
            } else {
                band_radius + 22.0
            };
            if feature_has_visible_export_label(&f) {
                let (lx, ly) = pos2xy(mid, len, cx, cy, label_radius.max(0.0));
                doc = doc.add(
                    Text::new(f.label)
                        .set("x", lx)
                        .set("y", ly)
                        .set("text-anchor", "middle")
                        .set("font-family", "monospace")
                        .set("font-size", CIRCULAR_FEATURE_LABEL_FONT_SIZE)
                        .set("fill", "#111111"),
                );
            }
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
        let mut keys: Vec<RestrictionEnzymeKey> = dna
            .restriction_enzyme_groups()
            .iter()
            .filter_map(|(key, names)| {
                DnaDisplay::restriction_group_matches_mode(
                    display.restriction_enzyme_display_mode,
                    &display.preferred_restriction_enzymes,
                    key,
                    names,
                )
                .then_some(key.clone())
            })
            .collect();
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
    use crate::{
        engine::{DisplaySettings, RestrictionEnzymeDisplayMode},
        enzymes::active_restriction_enzymes,
    };
    use gb_io::seq::Location;
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
            kind: "TFBS".into(),
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

    fn push_variation_feature(dna: &mut DNAsequence, label: &str, start: usize, end: usize) {
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "variation".into(),
            location: Location::simple_range(start as i64, end as i64),
            qualifiers: vec![("label".into(), Some(label.to_string()))],
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
            kind: "track".into(),
            location: Location::simple_range(start as i64, end as i64),
            qualifiers: vec![
                ("label".into(), Some(label.to_string())),
                (
                    "gentle_generated".into(),
                    Some("genome_vcf_track".to_string()),
                ),
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
    fn restriction_display_mode_applies_to_svg_export() {
        let mut dna = DNAsequence::from_sequence("GAATTCAAAAGAATTCAAAAGGATCCAAA").unwrap();
        *dna.restriction_enzymes_mut() = active_restriction_enzymes();
        dna.set_max_restriction_enzyme_sites(None);
        dna.update_computed_features();

        let mut display = DisplaySettings::default();
        display.show_restriction_enzymes = true;
        display.restriction_enzyme_display_mode = RestrictionEnzymeDisplayMode::PreferredOnly;
        display.preferred_restriction_enzymes = vec!["EcoRI".to_string()];

        let preferred_svg = export_linear_svg(&dna, &display);
        assert!(preferred_svg.contains("EcoRI"));
        assert!(!preferred_svg.contains("BamHI"));

        display.restriction_enzyme_display_mode = RestrictionEnzymeDisplayMode::AllInView;
        let all_svg = export_linear_svg(&dna, &display);
        assert!(all_svg.contains("EcoRI"));
        assert!(all_svg.contains("BamHI"));
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
    fn linear_svg_export_sizes_height_to_content() {
        let dna = DNAsequence::from_sequence("ATGCATGCATGCATGC").expect("sequence");
        let svg = export_linear_svg(&dna, &DisplaySettings::default());
        let height_marker = "height=\"";
        let height_start =
            svg.find(height_marker).expect("svg height attribute") + height_marker.len();
        let height_rest = &svg[height_start..];
        let height_end = height_rest.find('"').expect("svg height terminator");
        let height_value = height_rest[..height_end]
            .parse::<f32>()
            .expect("svg height parses");
        assert!(
            height_value < H,
            "content-sized linear export should be shorter than fixed 700px canvas, got {height_value}"
        );
        assert!(
            height_value >= LINEAR_SVG_MIN_HEIGHT,
            "content-sized linear export should respect minimum height, got {height_value}"
        );
    }

    #[test]
    fn linear_svg_export_honors_active_viewport_and_clips_features() {
        let mut dna = DNAsequence::from_sequence(&"ATGC".repeat(250)).expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: Location::simple_range(390, 420),
            qualifiers: vec![("label".into(), Some("visible_tx".to_string()))],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: Location::simple_range(10, 20),
            qualifiers: vec![("label".into(), Some("hidden_tx".to_string()))],
        });

        let mut display = DisplaySettings::default();
        display.show_gene_features = false;
        display.show_restriction_enzymes = false;
        display.show_gc_contents = false;
        display.linear_view_start_bp = 400;
        display.linear_view_span_bp = 100;

        let svg = export_linear_svg(&dna, &display);
        assert!(svg.contains("401..500 (100 bp view of 1000 bp)"));
        assert!(svg.contains("visible_tx"));
        assert!(!svg.contains("hidden_tx"));
    }

    #[test]
    fn linear_svg_export_renders_variation_as_baseline_marker() {
        let mut dna = DNAsequence::from_sequence(&"ATGC".repeat(80)).expect("sequence");
        push_variation_feature(&mut dna, "rs9923231", 120, 121);

        let svg = export_linear_svg(&dna, &DisplaySettings::default());
        assert!(svg.contains("rs9923231"));
        assert!(svg.contains("<circle"));
        assert!(svg.contains("stroke=\"#e17f0f\""));
        assert!(!svg.contains("<rect fill=\"#e17f0f\""));
    }

    #[test]
    fn linear_svg_export_marks_transcription_start_and_direction() {
        let mut dna = DNAsequence::from_sequence(&"ATGC".repeat(200)).expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "CDS".into(),
            location: Location::simple_range(20, 80),
            qualifiers: vec![("label".into(), Some("luc2".to_string()))],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: Location::Complement(Box::new(Location::simple_range(120, 180))),
            qualifiers: vec![("label".into(), Some("VKORC1".to_string()))],
        });

        let svg = export_linear_svg(&dna, &DisplaySettings::default());
        assert!(
            svg.matches("data-gentle-role=\"linear-transcription-start-tick\"")
                .count()
                >= 2
        );
        assert!(
            svg.matches("data-gentle-role=\"linear-transcription-start-arrow\"")
                .count()
                >= 2
        );
    }

    #[test]
    fn linear_svg_export_hides_fallback_coordinate_labels_and_compacts_nearby_tfbs_names() {
        let mut dna = DNAsequence::from_sequence(&"ATGC".repeat(120)).expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: Location::simple_range(30, 90),
            qualifiers: vec![],
        });
        push_tfbs_feature(&mut dna, "HNF4A", 40, 45, 2.0);
        push_tfbs_feature(&mut dna, "HNF4A", 42, 47, 2.0);

        let mut display = DisplaySettings::default();
        display.show_tfbs = true;
        let svg = export_linear_svg(&dna, &display);
        assert!(!svg.contains("30..90"));
        assert_eq!(svg.matches("HNF4A").count(), 1);
    }

    #[test]
    fn linear_svg_export_prefers_gene_symbols_over_accession_like_transcript_labels() {
        let mut dna = DNAsequence::from_sequence(&"ATGC".repeat(160)).expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: Location::Complement(Box::new(Location::simple_range(40, 120))),
            qualifiers: vec![
                ("label".into(), Some("ENST00000498155".to_string())),
                ("transcript_id".into(), Some("ENST00000498155".to_string())),
                ("gene".into(), Some("VKORC1".to_string())),
            ],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: Location::Complement(Box::new(Location::simple_range(44, 118))),
            qualifiers: vec![
                ("label".into(), Some("ENST00000420057".to_string())),
                ("transcript_id".into(), Some("ENST00000420057".to_string())),
                ("gene".into(), Some("VKORC1".to_string())),
            ],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: Location::simple_range(10, 70),
            qualifiers: vec![
                ("label".into(), Some("ENST00000624508".to_string())),
                ("transcript_id".into(), Some("ENST00000624508".to_string())),
            ],
        });

        let svg = export_linear_svg(&dna, &DisplaySettings::default());
        assert_eq!(svg.matches("VKORC1").count(), 1);
        assert!(!svg.contains("ENST00000498155"));
        assert!(!svg.contains("ENST00000420057"));
        assert!(!svg.contains("ENST00000624508"));
    }

    #[test]
    fn mrna_features_are_rendered_with_directional_pointed_bars() {
        let mrna = gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: Location::simple_range(20, 80),
            qualifiers: vec![],
        };
        assert!(feature_pointy(&mrna));
    }

    #[test]
    fn promoter_features_are_rendered_with_directional_pointed_bars() {
        let promoter = gb_io::seq::Feature {
            kind: "promoter".into(),
            location: Location::simple_range(20, 80),
            qualifiers: vec![],
        };
        assert!(feature_pointy(&promoter));
        assert!(has_transcription_direction(&promoter));
    }

    #[test]
    fn promoter_can_hide_export_tss_marker_via_qualifier() {
        let mut dna = DNAsequence::from_sequence(&"ATGC".repeat(200)).expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "promoter".into(),
            location: Location::simple_range(20, 80),
            qualifiers: vec![
                ("label".into(), Some("insert".to_string())),
                ("gentle_export_tss_marker".into(), Some("none".to_string())),
            ],
        });
        dna.set_circular(true);

        let svg = export_circular_svg(&dna, &DisplaySettings::default());
        assert!(!svg.contains("data-gentle-role=\"transcription-start-arrow-shaft\""));
        assert!(!svg.contains("data-gentle-role=\"transcription-start-arrow\""));
        assert!(svg.contains("insert"));
    }

    #[test]
    fn circular_svg_export_uses_transparent_background_and_marks_variation_and_tss() {
        let mut dna = DNAsequence::from_sequence(&"ATGC".repeat(200)).expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "CDS".into(),
            location: Location::simple_range(20, 80),
            qualifiers: vec![("label".into(), Some("luc2".to_string()))],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "CDS".into(),
            location: Location::Complement(Box::new(Location::simple_range(120, 180))),
            qualifiers: vec![("label".into(), Some("bla".to_string()))],
        });
        push_variation_feature(&mut dna, "rs9923231", 130, 131);
        dna.set_circular(true);

        let svg = export_circular_svg(&dna, &DisplaySettings::default());
        assert!(!svg.contains("<rect fill=\"#ffffff\""));
        assert!(svg.contains("data-gentle-role=\"variation-marker-line\""));
        assert!(svg.contains("data-gentle-role=\"variation-marker-dot\""));
        assert!(svg.contains("rs9923231"));
        assert!(
            svg.matches("data-gentle-role=\"transcription-start-arrow\"")
                .count()
                >= 2
        );
        assert!(
            svg.matches("data-gentle-role=\"transcription-start-arrow-shaft\"")
                .count()
                >= 2
        );
    }

    #[test]
    fn linear_svg_export_prefers_tfbs_bound_moiety_label() {
        let mut dna = DNAsequence::from_sequence(&"ATGC".repeat(80)).expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "TFBS".into(),
            location: Location::simple_range(40, 49),
            qualifiers: vec![
                ("label".into(), Some("TFBS MA0079.5".to_string())),
                ("tf_id".into(), Some("MA0079.5".to_string())),
                ("bound_moiety".into(), Some("SP1".to_string())),
                ("llr_bits".into(), Some("12.5".to_string())),
                ("llr_quantile".into(), Some("0.999".to_string())),
                ("true_log_odds_quantile".into(), Some("0.999".to_string())),
            ],
        });

        let mut display = DisplaySettings::default();
        display.show_tfbs = true;
        let svg = export_linear_svg(&dna, &display);
        assert!(svg.contains("SP1"));
        assert!(!svg.contains("TFBS MA0079.5"));
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
