//! Shared SVG renderer for continuous TF motif score-track plots.
//!
//! This keeps the Promoter design score-track figure available outside egui so
//! GUI, CLI, shell, and future adapters can all export the same visual rather
//! than only the raw per-position arrays.

use crate::engine::{TfbsScoreTrackReport, TfbsScoreTrackValueKind};

const SVG_WIDTH: f64 = 1180.0;
const SVG_MARGIN_LEFT: f64 = 32.0;
const SVG_MARGIN_RIGHT: f64 = 28.0;
const SVG_HEADER_HEIGHT: f64 = 82.0;
const SVG_FOOTER_HEIGHT: f64 = 42.0;
const SVG_ROW_HEIGHT: f64 = 96.0;
const SVG_ROW_GAP: f64 = 14.0;
const SVG_OVERLAY_ROW_HEIGHT: f64 = 34.0;
const SVG_OVERLAY_ROW_GAP: f64 = 8.0;
const SVG_LABEL_WIDTH: f64 = 320.0;
const SVG_LOGO_WIDTH: f64 = 84.0;
const SVG_LOGO_HEIGHT: f64 = 24.0;
const SVG_TRACK_PADDING_X: f64 = 10.0;
const SVG_TRACK_PADDING_Y: f64 = 8.0;
const SVG_MAX_POINTS_PER_POLYLINE: usize = 1400;
const SVG_TSS_MARKER_ROW_HEIGHT: f64 = 18.0;
const SVG_TSS_MARKER_MAX_LANES: usize = 3;
const CORR_SVG_WIDTH: f64 = 1360.0;
const CORR_MARGIN_LEFT: f64 = 28.0;
const CORR_MARGIN_RIGHT: f64 = 28.0;
const CORR_HEADER_HEIGHT: f64 = 112.0;
const CORR_FOOTER_HEIGHT: f64 = 54.0;
const CORR_PANEL_GAP: f64 = 34.0;
const CORR_CELL_SIZE: f64 = 42.0;
const CORR_LABEL_GAP: f64 = 8.0;

#[derive(Clone)]
struct SvgTssMarkerLabel {
    label: String,
    position_0based: usize,
    is_reverse: bool,
}

fn escape_svg_text(raw: &str) -> String {
    raw.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn sampled_point_indices(point_count: usize, max_points: usize) -> Vec<usize> {
    if point_count <= 1 || point_count <= max_points {
        return (0..point_count).collect();
    }
    let stride = point_count.div_ceil(max_points).max(1);
    let mut indices = (0..point_count).step_by(stride).collect::<Vec<_>>();
    let last_idx = point_count - 1;
    if indices.last().copied() != Some(last_idx) {
        indices.push(last_idx);
    }
    indices
}

fn format_bp_ticks(start_0based: usize, end_0based_exclusive: usize) -> (String, String, String) {
    let mid = start_0based + (end_0based_exclusive.saturating_sub(start_0based) / 2);
    (
        start_0based.to_string(),
        mid.to_string(),
        end_0based_exclusive.to_string(),
    )
}

fn format_positive_fraction_as_percent(fraction: f64) -> String {
    format!("{:.1}%", fraction.max(0.0) * 100.0)
}

fn format_tf_track_label(tf_id: &str, tf_name: Option<&str>) -> String {
    let trimmed_name = tf_name.map(str::trim).unwrap_or_default();
    if trimmed_name.is_empty() {
        return tf_id.to_string();
    }
    if trimmed_name.eq_ignore_ascii_case(tf_id) {
        trimmed_name.to_string()
    } else {
        format!("{trimmed_name} ({tf_id})")
    }
}

fn tfbs_logo_base_color(base: char) -> &'static str {
    match base {
        'A' => "#228b22",
        'C' => "#1e90ff",
        'G' => "#ff8c00",
        'T' => "#dc143c",
        _ => "#94a3b8",
    }
}

fn render_tf_track_logo_svg(
    columns: &[crate::engine::JasparExpertColumn],
    x: f64,
    y: f64,
    width: f64,
    height: f64,
) -> String {
    let mut svg = String::new();
    svg.push_str(&format!(
        "<g data-gentle-role=\"tfbs-score-track-logo\"><rect x=\"{x:.1}\" y=\"{y:.1}\" width=\"{width:.1}\" height=\"{height:.1}\" rx=\"4\" fill=\"#ffffff\" stroke=\"#cbd5e1\" stroke-width=\"1\"/></g>\n"
    ));
    if columns.is_empty() {
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">logo</text>\n",
            x + width * 0.5,
            y + height * 0.5,
        ));
        return svg;
    }
    let margin_x = 4.0_f64;
    let margin_y = 3.0_f64;
    let left = x + margin_x;
    let baseline = y + height - margin_y;
    let usable_width = (width - margin_x * 2.0).max(8.0);
    let usable_height = (height - margin_y * 2.0).max(8.0);
    let column_width = (usable_width / columns.len() as f64).max(3.0);
    let bits_to_px = usable_height / 2.0;
    svg.push_str(&format!(
        "<line x1=\"{left:.1}\" y1=\"{baseline:.1}\" x2=\"{:.1}\" y2=\"{baseline:.1}\" stroke=\"#cbd5e1\" stroke-width=\"0.8\"/>\n",
        x + width - margin_x,
    ));
    svg.push_str(&format!(
        "<line x1=\"{left:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#e2e8f0\" stroke-width=\"0.8\" stroke-dasharray=\"2 2\"/>\n",
        baseline - usable_height,
        x + width - margin_x,
        baseline - usable_height,
    ));
    for column in columns {
        let x_left = left + (column.position_1based as f64 - 1.0) * column_width;
        let x_center = x_left + column_width * 0.5;
        let information_height =
            (column.information_content_bits * bits_to_px).clamp(0.0, usable_height);
        if information_height > 0.75 {
            let alpha =
                ((column.information_content_bits / 2.0).clamp(0.24, 0.62) * 100.0).round() / 100.0;
            svg.push_str(&format!(
                "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" rx=\"1.5\" fill=\"#cbd5e1\" fill-opacity=\"{alpha:.2}\" stroke=\"#cbd5e1\" stroke-opacity=\"{:.2}\" stroke-width=\"0.4\"/>\n",
                x_left + column_width * 0.08,
                baseline - information_height,
                (column_width * 0.84).max(1.5),
                information_height,
                (alpha * 0.8).clamp(0.12, 0.55),
            ));
        }
        let mut rows = vec![
            ('A', column.a_logo_bits),
            ('C', column.c_logo_bits),
            ('G', column.g_logo_bits),
            ('T', column.t_logo_bits),
        ];
        rows.sort_by(|left, right| left.1.total_cmp(&right.1));
        let mut used_height = 0.0_f64;
        for (base, bits) in rows {
            if bits <= 0.0001 {
                continue;
            }
            let letter_height = (bits * bits_to_px).clamp(0.0, usable_height - used_height);
            if letter_height <= 0.75 {
                continue;
            }
            svg.push_str(&format!(
                "<text x=\"{x_center:.2}\" y=\"{:.2}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"{:.2}\" font-weight=\"700\" fill=\"{}\" data-gentle-role=\"tfbs-score-track-logo-letter\">{}</text>\n",
                baseline - used_height,
                letter_height,
                tfbs_logo_base_color(base),
                base,
            ));
            used_height += letter_height * 0.86;
        }
    }
    svg
}

fn format_track_normalization_summary(
    score_kind: TfbsScoreTrackValueKind,
    normalization: &crate::engine::TfbsScoreTrackNormalizationReference,
) -> String {
    match score_kind {
        TfbsScoreTrackValueKind::LlrBackgroundQuantile
        | TfbsScoreTrackValueKind::LlrBackgroundTailLog10
        | TfbsScoreTrackValueKind::TrueLogOddsBackgroundQuantile
        | TfbsScoreTrackValueKind::TrueLogOddsBackgroundTailLog10 => format!(
            "theory max {:.2} | peak q {:.6} | -log10 tail {:.2}",
            normalization.theoretical_max_score,
            normalization.observed_peak_modeled_quantile,
            normalization.observed_peak_modeled_tail_log10,
        ),
        _ => format!(
            "p99 {:.2} | Δp99 {:+.2} | bg+ {}",
            normalization.p99_score,
            normalization.observed_peak_delta_from_p99,
            format_positive_fraction_as_percent(normalization.positive_fraction)
        ),
    }
}

fn format_report_normalization_note(report: &TfbsScoreTrackReport) -> Option<String> {
    let normalization = report
        .tracks
        .iter()
        .find_map(|track| track.normalization_reference.as_ref())?;
    Some(format!(
        "normalization={} {}bp deterministic background | chance_model={}",
        normalization.background_model,
        normalization.random_sequence_length_bp,
        normalization.chance_model
    ))
}

fn format_overlay_track_meta(track: &crate::engine::TfbsScoreTrackOverlayTrack) -> String {
    let mut parts = vec![format!("{} interval(s)", track.interval_count)];
    if let Some(max_score) = track.max_score {
        parts.push(format!("max score {:.2}", max_score));
    }
    match track.source_file_name.as_deref() {
        Some(file_name) if file_name != track.track_name => {
            parts.push(format!("file {}", file_name));
        }
        _ => {}
    }
    format!("{} | {}", track.source_kind, parts.join(" | "))
}

fn overlay_track_interval_fill_opacity(
    interval: &crate::engine::TfbsScoreTrackOverlayInterval,
    track: &crate::engine::TfbsScoreTrackOverlayTrack,
) -> f64 {
    match (interval.score, track.max_score) {
        (Some(score), Some(max_score)) if max_score.is_finite() && max_score > 0.0 => {
            ((score / max_score).clamp(0.0, 1.0) * 0.5 + 0.25).clamp(0.25, 0.85)
        }
        _ => 0.6,
    }
}

fn correlation_fill_color(value: f64) -> &'static str {
    if value >= 0.85 {
        "#0f766e"
    } else if value >= 0.65 {
        "#14b8a6"
    } else if value >= 0.35 {
        "#5eead4"
    } else if value >= 0.10 {
        "#ccfbf1"
    } else if value > -0.10 {
        "#fff7ed"
    } else if value > -0.35 {
        "#e0e7ff"
    } else if value > -0.65 {
        "#a5b4fc"
    } else if value > -0.85 {
        "#818cf8"
    } else {
        "#5b21b6"
    }
}

fn correlation_text_color(value: f64) -> &'static str {
    if value.abs() >= 0.45 {
        "#fffbeb"
    } else {
        "#111827"
    }
}

fn correlation_value(
    row: &crate::engine::TfbsScoreTrackCorrelationRow,
    metric: crate::engine::TfbsScoreTrackCorrelationMetric,
    smoothed: bool,
) -> f64 {
    match (metric, smoothed) {
        (crate::engine::TfbsScoreTrackCorrelationMetric::Pearson, true) => row.smoothed_pearson,
        (crate::engine::TfbsScoreTrackCorrelationMetric::Pearson, false) => row.raw_pearson,
        (crate::engine::TfbsScoreTrackCorrelationMetric::Spearman, true) => row.smoothed_spearman,
        (crate::engine::TfbsScoreTrackCorrelationMetric::Spearman, false) => row.raw_spearman,
    }
}

fn cross_strand_correlation_value(
    cell: &crate::engine::TfbsScoreTrackCrossStrandCorrelationCell,
    metric: crate::engine::TfbsScoreTrackCorrelationMetric,
    smoothed: bool,
) -> f64 {
    match (metric, smoothed) {
        (crate::engine::TfbsScoreTrackCorrelationMetric::Pearson, true) => cell.smoothed_pearson,
        (crate::engine::TfbsScoreTrackCorrelationMetric::Pearson, false) => cell.raw_pearson,
        (crate::engine::TfbsScoreTrackCorrelationMetric::Spearman, true) => cell.smoothed_spearman,
        (crate::engine::TfbsScoreTrackCorrelationMetric::Spearman, false) => cell.raw_spearman,
    }
}

fn directional_correlation_value(
    summary: &crate::engine::TfbsScoreTrackDirectionalSummary,
    metric: crate::engine::TfbsScoreTrackCorrelationMetric,
    smoothed: bool,
) -> f64 {
    match (metric, smoothed) {
        (crate::engine::TfbsScoreTrackCorrelationMetric::Pearson, true) => summary.smoothed_pearson,
        (crate::engine::TfbsScoreTrackCorrelationMetric::Pearson, false) => summary.raw_pearson,
        (crate::engine::TfbsScoreTrackCorrelationMetric::Spearman, true) => {
            summary.smoothed_spearman
        }
        (crate::engine::TfbsScoreTrackCorrelationMetric::Spearman, false) => summary.raw_spearman,
    }
}

fn cross_strand_lookup_cell<'a>(
    row: &'a crate::engine::TfbsScoreTrackCrossStrandCorrelationRow,
    row_is_left: bool,
    row_strand: crate::engine::TfbsScoreTrackStrandComponent,
    col_strand: crate::engine::TfbsScoreTrackStrandComponent,
) -> Option<&'a crate::engine::TfbsScoreTrackCrossStrandCorrelationCell> {
    let (left_strand, right_strand) = if row_is_left {
        (row_strand, col_strand)
    } else {
        (col_strand, row_strand)
    };
    row.cells
        .iter()
        .find(|cell| cell.left_strand == left_strand && cell.right_strand == right_strand)
}

fn global_score_bounds(report: &TfbsScoreTrackReport) -> (f64, f64) {
    let mut min_score = 0.0_f64;
    let mut max_score = report.global_max_score.max(0.0);
    for track in &report.tracks {
        for score in track
            .forward_scores
            .iter()
            .chain(track.reverse_scores.iter())
        {
            if !score.is_finite() {
                continue;
            }
            min_score = min_score.min(*score);
            max_score = max_score.max(*score);
        }
    }
    if (max_score - min_score).abs() < f64::EPSILON {
        if max_score <= 0.0 {
            max_score = 1.0;
        } else {
            min_score = 0.0;
        }
    }
    (min_score, max_score)
}

fn score_to_y(score: f64, min_score: f64, max_score: f64, top: f64, bottom: f64) -> f64 {
    if (max_score - min_score).abs() < f64::EPSILON {
        return bottom;
    }
    let fraction = ((score - min_score) / (max_score - min_score)).clamp(0.0, 1.0);
    bottom - fraction * (bottom - top)
}

fn polyline_points(
    scores: &[f64],
    left: f64,
    right: f64,
    top: f64,
    bottom: f64,
    min_score: f64,
    max_score: f64,
) -> String {
    if scores.is_empty() {
        return String::new();
    }
    let indices = sampled_point_indices(scores.len(), SVG_MAX_POINTS_PER_POLYLINE);
    let denom = scores.len().saturating_sub(1).max(1) as f64;
    indices
        .into_iter()
        .map(|idx| {
            let x = left + (idx as f64 / denom) * (right - left);
            let y = score_to_y(scores[idx], min_score, max_score, top, bottom);
            format!("{x:.2},{y:.2}")
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn summarize_tss_label(raw: &str) -> String {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return "TSS".to_string();
    }
    let char_count = trimmed.chars().count();
    if char_count <= 16 {
        return trimmed.to_string();
    }
    let shortened = trimmed.chars().take(13).collect::<String>();
    format!("{shortened}...")
}

fn summarize_svg_tss_markers(report: &TfbsScoreTrackReport) -> Vec<SvgTssMarkerLabel> {
    let mut markers: Vec<SvgTssMarkerLabel> = vec![];
    for marker in &report.tss_markers {
        if let Some(existing) = markers.iter_mut().find(|existing| {
            existing.position_0based == marker.position_0based
                && existing.is_reverse == marker.is_reverse
        }) {
            if existing.label != marker.label && !existing.label.contains("(+") {
                existing.label = format!("{} (+1)", existing.label);
            }
            continue;
        }
        markers.push(SvgTssMarkerLabel {
            label: marker.label.clone(),
            position_0based: marker.position_0based,
            is_reverse: marker.is_reverse,
        });
    }
    markers.sort_by(|left, right| {
        left.position_0based
            .cmp(&right.position_0based)
            .then(left.label.cmp(&right.label))
    });
    markers
}

fn tss_marker_x(
    marker_position_0based: usize,
    view_start_0based: usize,
    view_end_0based_exclusive: usize,
    plot_left: f64,
    plot_right: f64,
) -> f64 {
    let span = view_end_0based_exclusive
        .saturating_sub(view_start_0based)
        .max(1) as f64;
    let clamped_pos = marker_position_0based.clamp(view_start_0based, view_end_0based_exclusive);
    plot_left
        + ((clamped_pos.saturating_sub(view_start_0based) as f64) / span) * (plot_right - plot_left)
}

fn layout_tss_marker_lanes(
    markers: &[SvgTssMarkerLabel],
    view_start_0based: usize,
    view_end_0based_exclusive: usize,
    plot_left: f64,
    plot_right: f64,
) -> usize {
    if markers.is_empty() {
        return 0;
    }
    let mut lane_right_edges = vec![f64::NEG_INFINITY; SVG_TSS_MARKER_MAX_LANES];
    for marker in markers {
        let x = tss_marker_x(
            marker.position_0based,
            view_start_0based,
            view_end_0based_exclusive,
            plot_left,
            plot_right,
        );
        let label = summarize_tss_label(&marker.label);
        let width = (label.chars().count() as f64) * 6.0 + 18.0;
        let left = x.min(plot_right - 4.0);
        let right = (left + width).min(plot_right);
        let mut assigned = false;
        for lane_idx in 0..SVG_TSS_MARKER_MAX_LANES {
            if left > lane_right_edges[lane_idx] + 10.0 {
                lane_right_edges[lane_idx] = right;
                assigned = true;
                break;
            }
        }
        if !assigned {
            lane_right_edges[SVG_TSS_MARKER_MAX_LANES - 1] =
                lane_right_edges[SVG_TSS_MARKER_MAX_LANES - 1].max(right);
        }
    }
    lane_right_edges
        .iter()
        .rposition(|edge| edge.is_finite())
        .map(|idx| idx + 1)
        .unwrap_or(0)
}

fn render_tss_marker_annotations(
    markers: &[SvgTssMarkerLabel],
    view_start_0based: usize,
    view_end_0based_exclusive: usize,
    plot_left: f64,
    plot_right: f64,
    marker_top: f64,
    line_top: f64,
    line_bottom: f64,
) -> String {
    if markers.is_empty() {
        return String::new();
    }
    let mut lane_right_edges = vec![f64::NEG_INFINITY; SVG_TSS_MARKER_MAX_LANES];
    let mut svg = String::new();
    for marker in markers {
        let x = tss_marker_x(
            marker.position_0based,
            view_start_0based,
            view_end_0based_exclusive,
            plot_left,
            plot_right,
        );
        let label = summarize_tss_label(&marker.label);
        let label_width = (label.chars().count() as f64) * 6.0 + 18.0;
        let mut lane_idx = SVG_TSS_MARKER_MAX_LANES - 1;
        let left = x.min(plot_right - 4.0);
        let right = (left + label_width).min(plot_right);
        for candidate in 0..SVG_TSS_MARKER_MAX_LANES {
            if left > lane_right_edges[candidate] + 10.0 {
                lane_idx = candidate;
                lane_right_edges[candidate] = right;
                break;
            }
        }
        let label_y = marker_top + lane_idx as f64 * SVG_TSS_MARKER_ROW_HEIGHT + 10.0;
        let stem_y = line_top - 8.0;
        let hook_y = label_y + 2.0;
        let elbow_x = if marker.is_reverse { x - 4.0 } else { x + 4.0 };
        let tip_x = if marker.is_reverse {
            x - 14.0
        } else {
            x + 14.0
        };
        let color = if marker.is_reverse {
            "#b45309"
        } else {
            "#0e7490"
        };
        svg.push_str(&format!(
            "<line x1=\"{x:.2}\" y1=\"{line_top:.2}\" x2=\"{x:.2}\" y2=\"{line_bottom:.2}\" stroke=\"{color}\" stroke-width=\"1.6\" stroke-dasharray=\"6 4\" opacity=\"0.95\" data-gentle-role=\"tfbs-score-track-tss-line\"/>\n"
        ));
        svg.push_str(&format!(
            "<path d=\"M{x:.2},{line_top:.2} L{x:.2},{stem_y:.2} Q{x:.2},{hook_y:.2} {elbow_x:.2},{hook_y:.2} L{tip_x:.2},{hook_y:.2}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"2.4\" stroke-linecap=\"round\" stroke-linejoin=\"round\" data-gentle-role=\"tfbs-score-track-tss-arrow\"/>\n"
        ));
        let (head_tip_x, base_x) = if marker.is_reverse {
            (tip_x - 2.2, tip_x + 3.4)
        } else {
            (tip_x + 2.2, tip_x - 3.4)
        };
        svg.push_str(&format!(
            "<path d=\"M{head_tip_x:.2},{hook_y:.2} L{base_x:.2},{:.2} L{base_x:.2},{:.2} z\" fill=\"{color}\" stroke=\"none\" data-gentle-role=\"tfbs-score-track-tss-arrow\"/>\n",
            hook_y - 3.4,
            hook_y + 3.4
        ));
        let label_x = if marker.is_reverse {
            tip_x - label_width - 6.0
        } else {
            tip_x + 6.0
        }
        .clamp(plot_left, (plot_right - label_width).max(plot_left));
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"16\" rx=\"3\" fill=\"#fffdf8\" fill-opacity=\"0.92\" stroke=\"{}\" stroke-width=\"0.8\" data-gentle-role=\"tfbs-score-track-tss-label-bg\"/>\n",
            label_x - 4.0,
            label_y - 11.0,
            label_width + 8.0,
            color
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"{}\">{}</text>\n",
            label_x,
            label_y,
            color,
            escape_svg_text(&format!("TSS {}", label))
        ));
    }
    svg
}

/// Render one stacked TFBS score-track figure as SVG text.
pub fn render_tfbs_score_tracks_svg(report: &TfbsScoreTrackReport) -> String {
    let summarized_tss_markers = summarize_svg_tss_markers(report);
    let tss_marker_lanes = layout_tss_marker_lanes(
        &summarized_tss_markers,
        report.view_start_0based,
        report.view_end_0based_exclusive,
        SVG_MARGIN_LEFT + SVG_LABEL_WIDTH,
        SVG_WIDTH - SVG_MARGIN_RIGHT,
    );
    let header_height = SVG_HEADER_HEIGHT + (tss_marker_lanes as f64) * SVG_TSS_MARKER_ROW_HEIGHT;
    let track_count = report.tracks.len().max(1);
    let overlay_count = report.overlay_tracks.len();
    let motif_block_height =
        track_count as f64 * SVG_ROW_HEIGHT + track_count.saturating_sub(1) as f64 * SVG_ROW_GAP;
    let overlay_block_height = if overlay_count == 0 {
        0.0
    } else {
        SVG_ROW_GAP
            + overlay_count as f64 * SVG_OVERLAY_ROW_HEIGHT
            + overlay_count.saturating_sub(1) as f64 * SVG_OVERLAY_ROW_GAP
    };
    let svg_height = header_height + SVG_FOOTER_HEIGHT + motif_block_height + overlay_block_height;
    let content_left = SVG_MARGIN_LEFT;
    let content_right = SVG_WIDTH - SVG_MARGIN_RIGHT;
    let plot_left = content_left + SVG_LABEL_WIDTH;
    let plot_right = content_right;
    let (min_score, max_score) = global_score_bounds(report);
    let zero_visible = min_score < 0.0 && max_score > 0.0;
    let (bp_start_text, bp_mid_text, bp_end_text) =
        format_bp_ticks(report.view_start_0based, report.view_end_0based_exclusive);
    let title = "Continuous TF motif score tracks";
    let target_label = if report.target_label.trim().is_empty() {
        report.seq_id.as_str()
    } else {
        report.target_label.trim()
    };
    let subtitle = format!(
        "target={} | span={}..{} | motifs={} | score={}{}",
        target_label,
        report.view_start_0based,
        report.view_end_0based_exclusive,
        report.tracks.len(),
        report.score_kind.as_str(),
        if report.clip_negative && report.score_kind.supports_negative_values() {
            " | positive-only"
        } else {
            ""
        }
    );
    let legend = if overlay_count > 0 {
        "forward strand = teal | reverse strand = amber | imported BED intervals = violet"
    } else {
        "forward strand = teal | reverse strand = amber"
    };
    let score_range_label = format!("{min_score:.2} .. {max_score:.2}");
    let normalization_note = format_report_normalization_note(report);

    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{SVG_WIDTH:.0}\" height=\"{svg_height:.0}\" viewBox=\"0 0 {SVG_WIDTH:.0} {svg_height:.0}\">\n"
    ));
    svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#fffdf8\"/>\n");
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"28\" font-family=\"monospace\" font-size=\"16\" fill=\"#111827\">{}</text>\n",
        content_left,
        escape_svg_text(title)
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"45\" font-family=\"monospace\" font-size=\"11\" fill=\"#475569\">{}</text>\n",
        content_left,
        escape_svg_text(&subtitle)
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"60\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">{}</text>\n",
        content_left,
        escape_svg_text(legend)
    ));
    if let Some(note) = normalization_note {
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"73\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">{}</text>\n",
            content_left,
            escape_svg_text(&note)
        ));
    }
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"60\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">score range {}</text>\n",
        content_right,
        escape_svg_text(&score_range_label)
    ));

    for (row_idx, track) in report.tracks.iter().enumerate() {
        let row_top = header_height + row_idx as f64 * (SVG_ROW_HEIGHT + SVG_ROW_GAP);
        let row_bottom = row_top + SVG_ROW_HEIGHT;
        let plot_top = row_top + SVG_TRACK_PADDING_Y;
        let plot_bottom = row_bottom - 20.0;
        let label = format_tf_track_label(&track.tf_id, track.tf_name.as_deref());
        let max_position = track
            .max_position_0based
            .map(|value| format!(" @ {value}"))
            .unwrap_or_default();
        let title_x = content_left + 12.0;
        let logo_x = content_left + 12.0;
        let logo_y = row_top + 30.0;
        let meta_x = logo_x + SVG_LOGO_WIDTH + 12.0;
        let meta = format!(
            "{} bp motif | {} windows | max {:.2}{}",
            track.motif_length_bp, track.scored_window_count, track.max_score, max_position
        );
        let normalization_meta = track.normalization_reference.as_ref().map(|normalization| {
            format_track_normalization_summary(report.score_kind, normalization)
        });
        let row_fill = if row_idx % 2 == 0 {
            "#fffaf0"
        } else {
            "#fffefb"
        };
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"6\" fill=\"{}\" stroke=\"#e2e8f0\" stroke-width=\"1\"/>\n",
            content_left,
            row_top,
            content_right - content_left,
            SVG_ROW_HEIGHT,
            row_fill
        ));
        svg.push_str(&render_tf_track_logo_svg(
            &track.motif_logo_columns,
            logo_x,
            logo_y,
            SVG_LOGO_WIDTH,
            SVG_LOGO_HEIGHT,
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"13\" fill=\"#0f172a\">{}</text>\n",
            title_x,
            row_top + 24.0,
            escape_svg_text(&label)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">{}</text>\n",
            meta_x,
            row_top + 46.0,
            escape_svg_text(&meta)
        ));
        if let Some(normalization_meta) = normalization_meta {
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">{}</text>\n",
                meta_x,
                row_top + 62.0,
                escape_svg_text(&normalization_meta)
            ));
        }

        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"4\" fill=\"#ffffff\" stroke=\"#cbd5e1\" stroke-width=\"1\"/>\n",
            plot_left,
            plot_top,
            plot_right - plot_left,
            plot_bottom - plot_top
        ));

        for fraction in [0.25_f64, 0.5, 0.75] {
            let y = plot_bottom - fraction * (plot_bottom - plot_top);
            svg.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#f1f5f9\" stroke-width=\"1\"/>\n",
                plot_left, y, plot_right, y
            ));
        }

        let baseline_y = if zero_visible {
            score_to_y(0.0, min_score, max_score, plot_top, plot_bottom)
        } else {
            plot_bottom
        };
        if zero_visible {
            svg.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#94a3b8\" stroke-width=\"1\" stroke-dasharray=\"4 3\"/>\n",
                plot_left, baseline_y, plot_right, baseline_y
            ));
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"9\" fill=\"#64748b\">0.0</text>\n",
                plot_right - 6.0,
                baseline_y - 4.0,
            ));
        } else {
            svg.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#94a3b8\" stroke-width=\"1\"/>\n",
                plot_left, baseline_y, plot_right, baseline_y
            ));
        }

        if track.scored_window_count == 0 {
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"11\" fill=\"#64748b\">motif longer than selected range</text>\n",
                (plot_left + plot_right) * 0.5,
                (plot_top + plot_bottom) * 0.5
            ));
            continue;
        }

        let forward_points = polyline_points(
            &track.forward_scores,
            plot_left + SVG_TRACK_PADDING_X,
            plot_right - SVG_TRACK_PADDING_X,
            plot_top + 4.0,
            plot_bottom - 4.0,
            min_score,
            max_score,
        );
        let reverse_points = polyline_points(
            &track.reverse_scores,
            plot_left + SVG_TRACK_PADDING_X,
            plot_right - SVG_TRACK_PADDING_X,
            plot_top + 4.0,
            plot_bottom - 4.0,
            min_score,
            max_score,
        );
        if !forward_points.is_empty() {
            svg.push_str(&format!(
                "<polyline fill=\"none\" stroke=\"#0e7490\" stroke-width=\"1.8\" stroke-linejoin=\"round\" stroke-linecap=\"round\" points=\"{}\"/>\n",
                forward_points
            ));
        }
        if !reverse_points.is_empty() {
            svg.push_str(&format!(
                "<polyline fill=\"none\" stroke=\"#b45309\" stroke-width=\"1.6\" stroke-linejoin=\"round\" stroke-linecap=\"round\" points=\"{}\"/>\n",
                reverse_points
            ));
        }
    }

    if !report.overlay_tracks.is_empty() {
        let overlay_block_top = header_height + motif_block_height + SVG_ROW_GAP;
        for (overlay_idx, overlay_track) in report.overlay_tracks.iter().enumerate() {
            let row_top = overlay_block_top
                + overlay_idx as f64 * (SVG_OVERLAY_ROW_HEIGHT + SVG_OVERLAY_ROW_GAP);
            let row_bottom = row_top + SVG_OVERLAY_ROW_HEIGHT;
            let plot_top = row_top + 6.0;
            let plot_bottom = row_bottom - 6.0;
            let row_fill = if overlay_idx % 2 == 0 {
                "#faf5ff"
            } else {
                "#fdfbff"
            };
            svg.push_str(&format!(
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"6\" fill=\"{}\" stroke=\"#e9d5ff\" stroke-width=\"1\" data-gentle-role=\"tfbs-score-track-overlay-row\"/>\n",
                content_left,
                row_top,
                content_right - content_left,
                SVG_OVERLAY_ROW_HEIGHT,
                row_fill
            ));
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"11\" fill=\"#581c87\">{}</text>\n",
                content_left + 12.0,
                row_top + 15.0,
                escape_svg_text(&overlay_track.display_label)
            ));
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"9\" fill=\"#7c3aed\">{}</text>\n",
                content_left + 12.0,
                row_top + 27.0,
                escape_svg_text(&format_overlay_track_meta(overlay_track))
            ));
            svg.push_str(&format!(
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"4\" fill=\"#ffffff\" stroke=\"#ddd6fe\" stroke-width=\"1\" data-gentle-role=\"tfbs-score-track-overlay-plot\"/>\n",
                plot_left,
                plot_top,
                plot_right - plot_left,
                plot_bottom - plot_top
            ));
            let lane_top = plot_top + 5.0;
            let lane_bottom = plot_bottom - 5.0;
            let lane_height = (lane_bottom - lane_top).max(2.0);
            for interval in &overlay_track.intervals {
                let start_fraction = if report.view_end_0based_exclusive > report.view_start_0based
                {
                    (interval
                        .start_0based
                        .saturating_sub(report.view_start_0based) as f64)
                        / (report.view_end_0based_exclusive - report.view_start_0based) as f64
                } else {
                    0.0
                };
                let end_fraction = if report.view_end_0based_exclusive > report.view_start_0based {
                    (interval
                        .end_0based_exclusive
                        .saturating_sub(report.view_start_0based) as f64)
                        / (report.view_end_0based_exclusive - report.view_start_0based) as f64
                } else {
                    1.0
                };
                let x1 = plot_left
                    + SVG_TRACK_PADDING_X
                    + (plot_right - plot_left - SVG_TRACK_PADDING_X * 2.0)
                        * start_fraction.clamp(0.0, 1.0);
                let x2 = plot_left
                    + SVG_TRACK_PADDING_X
                    + (plot_right - plot_left - SVG_TRACK_PADDING_X * 2.0)
                        * end_fraction.clamp(0.0, 1.0);
                let width = (x2 - x1).max(2.0);
                let opacity = overlay_track_interval_fill_opacity(interval, overlay_track);
                svg.push_str(&format!(
                    "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" rx=\"2\" fill=\"#7c3aed\" fill-opacity=\"{:.2}\" stroke=\"#5b21b6\" stroke-width=\"0.6\" data-gentle-role=\"tfbs-score-track-overlay-interval\"/>\n",
                    x1,
                    lane_top,
                    width,
                    lane_height,
                    opacity
                ));
            }
        }
    }

    if tss_marker_lanes > 0 {
        let overlay_bottom = if overlay_count > 0 {
            header_height + motif_block_height + overlay_block_height - 10.0
        } else {
            header_height + motif_block_height - 10.0
        };
        svg.push_str(&render_tss_marker_annotations(
            &summarized_tss_markers,
            report.view_start_0based,
            report.view_end_0based_exclusive,
            plot_left + SVG_TRACK_PADDING_X,
            plot_right - SVG_TRACK_PADDING_X,
            SVG_HEADER_HEIGHT - 2.0,
            header_height - 6.0,
            overlay_bottom,
        ));
    }

    let axis_y = svg_height - SVG_FOOTER_HEIGHT + 8.0;
    svg.push_str(&format!(
        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#475569\" stroke-width=\"1\"/>\n",
        plot_left, axis_y, plot_right, axis_y
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">{}</text>\n",
        plot_left,
        axis_y + 16.0,
        escape_svg_text(&bp_start_text)
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">{}</text>\n",
        (plot_left + plot_right) * 0.5,
        axis_y + 16.0,
        escape_svg_text(&bp_mid_text)
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">{}</text>\n",
        plot_right,
        axis_y + 16.0,
        escape_svg_text(&bp_end_text)
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">base-pair position in selected span</text>\n",
        content_left,
        svg_height - 8.0
    ));
    svg.push_str("</svg>\n");
    svg
}

/// Render one TFBS-track correlation summary as SVG text.
pub fn render_tfbs_score_track_correlation_svg(
    report: &TfbsScoreTrackReport,
    metric: crate::engine::TfbsScoreTrackCorrelationMetric,
    signal_source: crate::engine::TfbsScoreTrackCorrelationSignalSource,
) -> String {
    let cross_strand_correlation = report
        .cross_strand_correlation_summary
        .as_ref()
        .filter(|_| {
            signal_source == crate::engine::TfbsScoreTrackCorrelationSignalSource::MaxStrands
        });
    let correlation = report
        .correlation_summaries
        .iter()
        .find(|summary| summary.signal_source == signal_source)
        .or_else(|| {
            report.correlation_summary.as_ref().filter(|summary| {
                summary.signal_source == signal_source
                    || signal_source
                        == crate::engine::TfbsScoreTrackCorrelationSignalSource::MaxStrands
            })
        });
    if correlation.is_none() && cross_strand_correlation.is_none() {
        return format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"960\" height=\"220\" viewBox=\"0 0 960 220\">\n\
<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#fffdf8\"/>\n\
<text x=\"28\" y=\"34\" font-family=\"monospace\" font-size=\"16\" fill=\"#111827\">TFBS track correlation</text>\n\
<text x=\"28\" y=\"58\" font-family=\"monospace\" font-size=\"11\" fill=\"#475569\">seq={} | score={} | span={}..{}</text>\n\
<text x=\"28\" y=\"108\" font-family=\"monospace\" font-size=\"12\" fill=\"#64748b\">No pairwise correlation summary is available for signal={} in this report.</text>\n\
</svg>\n",
            escape_svg_text(&report.seq_id),
            escape_svg_text(report.score_kind.as_str()),
            report.view_start_0based,
            report.view_end_0based_exclusive,
            escape_svg_text(signal_source.as_str()),
        );
    }
    if report.tracks.is_empty() {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"960\" height=\"220\" viewBox=\"0 0 960 220\">\n\
<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#fffdf8\"/>\n\
<text x=\"28\" y=\"34\" font-family=\"monospace\" font-size=\"16\" fill=\"#111827\">TFBS track correlation</text>\n\
<text x=\"28\" y=\"108\" font-family=\"monospace\" font-size=\"12\" fill=\"#64748b\">No tracks were available to correlate.</text>\n\
</svg>\n"
            .to_string();
    }

    if let Some(cross_strand_correlation) = cross_strand_correlation {
        let components = report
            .tracks
            .iter()
            .enumerate()
            .flat_map(|(track_idx, track)| {
                let label = format_tf_track_label(&track.tf_id, track.tf_name.as_deref());
                [
                    (
                        track_idx,
                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                        format!("{label} F"),
                    ),
                    (
                        track_idx,
                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                        format!("{label} R"),
                    ),
                ]
            })
            .collect::<Vec<_>>();
        let matrix_size = components.len();
        let label_width = components
            .iter()
            .map(|(_, _, label)| label.chars().count() as f64 * 6.8)
            .fold(160.0_f64, f64::max)
            .min(280.0);
        let panel_size = matrix_size as f64 * CORR_CELL_SIZE;
        let panel_width = label_width + CORR_LABEL_GAP + panel_size;
        let content_width =
            CORR_MARGIN_LEFT + panel_width * 2.0 + CORR_PANEL_GAP + CORR_MARGIN_RIGHT;
        let svg_width = CORR_SVG_WIDTH.max(content_width);

        let mut ranked_cells = vec![];
        for row in &cross_strand_correlation.rows {
            for cell in &row.cells {
                ranked_cells.push((row, cell));
            }
        }
        ranked_cells.sort_by(|(left_row, left_cell), (right_row, right_cell)| {
            cross_strand_correlation_value(right_cell, metric, true)
                .abs()
                .total_cmp(&cross_strand_correlation_value(left_cell, metric, true).abs())
                .then(
                    cross_strand_correlation_value(right_cell, metric, false)
                        .abs()
                        .total_cmp(&cross_strand_correlation_value(left_cell, metric, false).abs()),
                )
                .then(left_row.left_tf_id.cmp(&right_row.left_tf_id))
                .then(left_row.right_tf_id.cmp(&right_row.right_tf_id))
                .then(
                    left_cell
                        .left_strand
                        .as_str()
                        .cmp(right_cell.left_strand.as_str()),
                )
                .then(
                    left_cell
                        .right_strand
                        .as_str()
                        .cmp(right_cell.right_strand.as_str()),
                )
        });
        let list_rows = ranked_cells.len().min(8);
        let list_height = 34.0 + list_rows as f64 * 18.0;
        let svg_height = CORR_HEADER_HEIGHT + panel_size + list_height + CORR_FOOTER_HEIGHT;

        let left_panel_left = CORR_MARGIN_LEFT;
        let right_panel_left = left_panel_left + panel_width + CORR_PANEL_GAP;
        let matrix_top = CORR_HEADER_HEIGHT;
        let matrix_left = |panel_left: f64| panel_left + label_width + CORR_LABEL_GAP;

        let mut cross_row_lookup = std::collections::HashMap::<
            (usize, usize),
            &crate::engine::TfbsScoreTrackCrossStrandCorrelationRow,
        >::new();
        for row in &cross_strand_correlation.rows {
            let left_idx = report
                .tracks
                .iter()
                .position(|track| track.tf_id == row.left_tf_id);
            let right_idx = report
                .tracks
                .iter()
                .position(|track| track.tf_id == row.right_tf_id);
            if let (Some(left_idx), Some(right_idx)) = (left_idx, right_idx) {
                cross_row_lookup.insert((left_idx.min(right_idx), left_idx.max(right_idx)), row);
            }
        }

        let mut svg = String::new();
        svg.push_str(&format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{svg_width:.0}\" height=\"{svg_height:.0}\" viewBox=\"0 0 {svg_width:.0} {svg_height:.0}\">\n"
        ));
        svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#fffdf8\"/>\n");
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"34\" font-family=\"monospace\" font-size=\"16\" fill=\"#111827\">TFBS track correlation</text>\n",
            CORR_MARGIN_LEFT
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"54\" font-family=\"monospace\" font-size=\"11\" fill=\"#475569\">seq={} | span={}..{} | score={} | pair_count={}</text>\n",
            CORR_MARGIN_LEFT,
            escape_svg_text(&report.seq_id),
            report.view_start_0based,
            report.view_end_0based_exclusive,
            escape_svg_text(report.score_kind.as_str()),
            cross_strand_correlation.pair_count
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"71\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">Smoothed {} is the main neighborhood-synchrony view; raw {} stays alongside it as the strict per-position check.</text>\n",
            CORR_MARGIN_LEFT,
            escape_svg_text(metric.as_str()),
            escape_svg_text(metric.as_str())
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"86\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">rows and columns list each motif twice in F then R order | smoothing={} {} bp</text>\n",
            CORR_MARGIN_LEFT,
            escape_svg_text(&cross_strand_correlation.smoothing_method),
            cross_strand_correlation.smoothing_window_bp
        ));

        for (panel_left, panel_title, is_smoothed) in [
            (left_panel_left, metric.display_label(true), true),
            (right_panel_left, metric.display_label(false), false),
        ] {
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"106\" font-family=\"monospace\" font-size=\"12\" fill=\"#0f172a\">{}</text>\n",
                panel_left,
                escape_svg_text(panel_title)
            ));
            let grid_left = matrix_left(panel_left);
            svg.push_str(&format!(
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"6\" fill=\"#ffffff\" stroke=\"#cbd5e1\" stroke-width=\"1\"/>\n",
                grid_left,
                matrix_top,
                panel_size,
                panel_size
            ));
            for tf_group_idx in 0..report.tracks.len() {
                let block_x = grid_left + (tf_group_idx * 2) as f64 * CORR_CELL_SIZE;
                let block_y = matrix_top + (tf_group_idx * 2) as f64 * CORR_CELL_SIZE;
                let block_size = CORR_CELL_SIZE * 2.0;
                svg.push_str(&format!(
                    "<rect x=\"{block_x:.1}\" y=\"{block_y:.1}\" width=\"{block_size:.1}\" height=\"{block_size:.1}\" fill=\"#fef3c7\" fill-opacity=\"0.18\" stroke=\"#92400e\" stroke-width=\"1.4\" stroke-dasharray=\"3 2\"/>\n"
                ));
            }
            for (idx, (_, _, label)) in components.iter().enumerate() {
                let y = matrix_top + idx as f64 * CORR_CELL_SIZE + CORR_CELL_SIZE * 0.65;
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#334155\">{}</text>\n",
                    grid_left - 6.0,
                    y,
                    escape_svg_text(label)
                ));
                let column_x = grid_left + idx as f64 * CORR_CELL_SIZE + CORR_CELL_SIZE * 0.35;
                let label_y = matrix_top - 8.0;
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" transform=\"rotate(-45 {:.1} {:.1})\" font-family=\"monospace\" font-size=\"9\" fill=\"#334155\">{}</text>\n",
                    column_x,
                    label_y,
                    column_x,
                    label_y,
                    escape_svg_text(label)
                ));
            }

            for row_idx in 0..matrix_size {
                for col_idx in 0..matrix_size {
                    let x = grid_left + col_idx as f64 * CORR_CELL_SIZE;
                    let y = matrix_top + row_idx as f64 * CORR_CELL_SIZE;
                    let (row_track_idx, row_strand, _) = &components[row_idx];
                    let (col_track_idx, col_strand, _) = &components[col_idx];
                    let value = if row_track_idx == col_track_idx {
                        if row_strand == col_strand {
                            1.0
                        } else {
                            report.tracks[*row_track_idx]
                                .directional_summary
                                .as_ref()
                                .map(|summary| {
                                    directional_correlation_value(summary, metric, is_smoothed)
                                })
                                .unwrap_or(0.0)
                        }
                    } else {
                        let key = (
                            (*row_track_idx).min(*col_track_idx),
                            (*row_track_idx).max(*col_track_idx),
                        );
                        let row_is_left = row_track_idx < col_track_idx;
                        cross_row_lookup
                            .get(&key)
                            .and_then(|row| {
                                cross_strand_lookup_cell(row, row_is_left, *row_strand, *col_strand)
                            })
                            .map(|cell| cross_strand_correlation_value(cell, metric, is_smoothed))
                            .unwrap_or(0.0)
                    };
                    svg.push_str(&format!(
                        "<rect x=\"{x:.1}\" y=\"{y:.1}\" width=\"{CORR_CELL_SIZE:.1}\" height=\"{CORR_CELL_SIZE:.1}\" fill=\"{}\" stroke=\"#e2e8f0\" stroke-width=\"0.8\"/>\n",
                        correlation_fill_color(value)
                    ));
                    svg.push_str(&format!(
                        "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"8.2\" fill=\"{}\">{:+.2}</text>\n",
                        x + CORR_CELL_SIZE * 0.5,
                        y + CORR_CELL_SIZE * 0.60,
                        correlation_text_color(value),
                        value
                    ));
                    svg.push_str(&format!(
                        "<rect x=\"{x:.1}\" y=\"{y:.1}\" width=\"{CORR_CELL_SIZE:.1}\" height=\"{CORR_CELL_SIZE:.1}\" fill=\"none\" stroke=\"#cbd5e1\" stroke-width=\"1\"/>\n"
                    ));
                }
            }
            for separator_idx in (2..matrix_size).step_by(2) {
                let separator_x = grid_left + separator_idx as f64 * CORR_CELL_SIZE;
                let separator_y = matrix_top + separator_idx as f64 * CORR_CELL_SIZE;
                svg.push_str(&format!(
                    "<line x1=\"{separator_x:.1}\" y1=\"{matrix_top:.1}\" x2=\"{separator_x:.1}\" y2=\"{:.1}\" stroke=\"#475569\" stroke-width=\"2.2\" opacity=\"0.95\"/>\n",
                    matrix_top + panel_size
                ));
                svg.push_str(&format!(
                    "<line x1=\"{grid_left:.1}\" y1=\"{separator_y:.1}\" x2=\"{:.1}\" y2=\"{separator_y:.1}\" stroke=\"#475569\" stroke-width=\"2.2\" opacity=\"0.95\"/>\n",
                    grid_left + panel_size
                ));
            }
        }

        let legend_y = matrix_top + panel_size + 20.0;
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">color scale: strong negative -> neutral -> strong positive | rows/columns are ordered F then R for every motif, so each TF-pair block reads F-F / F-R / R-F / R-R</text>\n",
            CORR_MARGIN_LEFT,
            legend_y
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#92400e\">amber dashed 2x2 diagonal blocks mark one motif's own F/R self-comparison; the off-diagonal cells inside those blocks are the mutual F↔R direction scores</text>\n",
            CORR_MARGIN_LEFT,
            legend_y + 14.0
        ));
        for (idx, (value, label)) in [
            (-1.0, "-1.0"),
            (-0.5, "-0.5"),
            (0.0, "0"),
            (0.5, "+0.5"),
            (1.0, "+1.0"),
        ]
        .iter()
        .enumerate()
        {
            let x = CORR_MARGIN_LEFT + 430.0 + idx as f64 * 66.0;
            svg.push_str(&format!(
                "<rect x=\"{x:.1}\" y=\"{:.1}\" width=\"18\" height=\"12\" fill=\"{}\" stroke=\"#cbd5e1\" stroke-width=\"1\"/>\n",
                legend_y - 10.0,
                correlation_fill_color(*value)
            ));
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"9\" fill=\"#475569\">{}</text>\n",
                x + 24.0,
                legend_y,
                label
            ));
        }

        let list_top = legend_y + 36.0;
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"12\" fill=\"#0f172a\">Top synchronized strand pairs</text>\n",
            CORR_MARGIN_LEFT,
            list_top
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">ranked by |smoothed {}|, then |raw {}|; offset uses strand-specific primary peak start positions</text>\n",
            CORR_MARGIN_LEFT,
            list_top + 14.0,
            escape_svg_text(metric.as_str()),
            escape_svg_text(metric.as_str())
        ));
        let mut line_y = list_top + 32.0;
        for (row, cell) in ranked_cells.into_iter().take(list_rows) {
            let left_label = format_tf_track_label(&row.left_tf_id, row.left_tf_name.as_deref());
            let right_label = format_tf_track_label(&row.right_tf_id, row.right_tf_name.as_deref());
            let offset = cell
                .signed_primary_peak_offset_bp
                .map(|value| format!("{value:+} bp"))
                .unwrap_or_else(|| "n/a".to_string());
            let summary = format!(
                "{} {} vs {} {} | smoothed {:+.3} | raw {:+.3} | offset {}",
                left_label,
                cell.left_strand.short_label(),
                right_label,
                cell.right_strand.short_label(),
                cross_strand_correlation_value(cell, metric, true),
                cross_strand_correlation_value(cell, metric, false),
                offset
            );
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#334155\">{}</text>\n",
                CORR_MARGIN_LEFT,
                line_y,
                escape_svg_text(&summary)
            ));
            line_y += 18.0;
        }
        svg.push_str("</svg>\n");
        return svg;
    }

    let correlation = correlation.expect("checked above");

    let labels = report
        .tracks
        .iter()
        .map(|track| format_tf_track_label(&track.tf_id, track.tf_name.as_deref()))
        .collect::<Vec<_>>();
    let matrix_size = labels.len();
    let label_width = labels
        .iter()
        .map(|label| label.chars().count() as f64 * 6.8)
        .fold(160.0_f64, f64::max)
        .min(280.0);
    let panel_size = matrix_size as f64 * CORR_CELL_SIZE;
    let panel_width = label_width + CORR_LABEL_GAP + panel_size;
    let content_width = CORR_MARGIN_LEFT + panel_width * 2.0 + CORR_PANEL_GAP + CORR_MARGIN_RIGHT;
    let svg_width = CORR_SVG_WIDTH.max(content_width);
    let list_rows = correlation.rows.len().min(8);
    let list_height = 34.0 + list_rows as f64 * 18.0;
    let svg_height = CORR_HEADER_HEIGHT + panel_size + list_height + CORR_FOOTER_HEIGHT;

    let left_panel_left = CORR_MARGIN_LEFT;
    let right_panel_left = left_panel_left + panel_width + CORR_PANEL_GAP;
    let matrix_top = CORR_HEADER_HEIGHT;
    let matrix_left = |panel_left: f64| panel_left + label_width + CORR_LABEL_GAP;

    let mut row_lookup = std::collections::HashMap::<
        (usize, usize),
        &crate::engine::TfbsScoreTrackCorrelationRow,
    >::new();
    for row in &correlation.rows {
        let left_idx = report
            .tracks
            .iter()
            .position(|track| track.tf_id == row.left_tf_id);
        let right_idx = report
            .tracks
            .iter()
            .position(|track| track.tf_id == row.right_tf_id);
        if let (Some(left_idx), Some(right_idx)) = (left_idx, right_idx) {
            row_lookup.insert((left_idx.min(right_idx), left_idx.max(right_idx)), row);
        }
    }

    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{svg_width:.0}\" height=\"{svg_height:.0}\" viewBox=\"0 0 {svg_width:.0} {svg_height:.0}\">\n"
    ));
    svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#fffdf8\"/>\n");
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"34\" font-family=\"monospace\" font-size=\"16\" fill=\"#111827\">TFBS track correlation</text>\n",
        CORR_MARGIN_LEFT
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"54\" font-family=\"monospace\" font-size=\"11\" fill=\"#475569\">seq={} | span={}..{} | score={} | pair_count={}</text>\n",
        CORR_MARGIN_LEFT,
        escape_svg_text(&report.seq_id),
        report.view_start_0based,
        report.view_end_0based_exclusive,
        escape_svg_text(report.score_kind.as_str()),
        correlation.pair_count
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"71\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">Smoothed {} is the main neighborhood-synchrony view; raw {} stays alongside it as the strict per-position check.</text>\n",
        CORR_MARGIN_LEFT,
        escape_svg_text(metric.as_str()),
        escape_svg_text(metric.as_str())
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"86\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">signal={} | smoothing={} {} bp | diagonal fixed at 1.000</text>\n",
        CORR_MARGIN_LEFT,
        escape_svg_text(correlation.signal_source.summary_label()),
        escape_svg_text(&correlation.smoothing_method),
        correlation.smoothing_window_bp
    ));

    for (panel_left, panel_title, is_smoothed) in [
        (left_panel_left, metric.display_label(true), true),
        (right_panel_left, metric.display_label(false), false),
    ] {
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"106\" font-family=\"monospace\" font-size=\"12\" fill=\"#0f172a\">{}</text>\n",
            panel_left,
            escape_svg_text(panel_title)
        ));
        let grid_left = matrix_left(panel_left);
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"6\" fill=\"#ffffff\" stroke=\"#cbd5e1\" stroke-width=\"1\"/>\n",
            grid_left,
            matrix_top,
            panel_size,
            panel_size
        ));
        for (idx, label) in labels.iter().enumerate() {
            let y = matrix_top + idx as f64 * CORR_CELL_SIZE + CORR_CELL_SIZE * 0.65;
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#334155\">{}</text>\n",
                grid_left - 6.0,
                y,
                escape_svg_text(label)
            ));
            let column_x = grid_left + idx as f64 * CORR_CELL_SIZE + CORR_CELL_SIZE * 0.35;
            let label_y = matrix_top - 8.0;
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" transform=\"rotate(-45 {:.1} {:.1})\" font-family=\"monospace\" font-size=\"9\" fill=\"#334155\">{}</text>\n",
                column_x,
                label_y,
                column_x,
                label_y,
                escape_svg_text(label)
            ));
        }

        for row_idx in 0..matrix_size {
            for col_idx in 0..matrix_size {
                let value = if row_idx == col_idx {
                    1.0
                } else {
                    let lookup = row_lookup
                        .get(&(row_idx.min(col_idx), row_idx.max(col_idx)))
                        .copied();
                    match lookup {
                        Some(row) => correlation_value(row, metric, is_smoothed),
                        None => 0.0,
                    }
                };
                let x = grid_left + col_idx as f64 * CORR_CELL_SIZE;
                let y = matrix_top + row_idx as f64 * CORR_CELL_SIZE;
                svg.push_str(&format!(
                    "<rect x=\"{x:.1}\" y=\"{y:.1}\" width=\"{CORR_CELL_SIZE:.1}\" height=\"{CORR_CELL_SIZE:.1}\" fill=\"{}\" stroke=\"#e2e8f0\" stroke-width=\"1\"/>\n",
                    correlation_fill_color(value)
                ));
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"9\" fill=\"{}\">{:+.2}</text>\n",
                    x + CORR_CELL_SIZE * 0.5,
                    y + CORR_CELL_SIZE * 0.60,
                    correlation_text_color(value),
                    value
                ));
            }
        }
    }

    let legend_y = matrix_top + panel_size + 20.0;
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">color scale: strong negative -> neutral -> strong positive</text>\n",
        CORR_MARGIN_LEFT,
        legend_y
    ));
    for (idx, (value, label)) in [
        (-1.0, "-1.0"),
        (-0.5, "-0.5"),
        (0.0, "0"),
        (0.5, "+0.5"),
        (1.0, "+1.0"),
    ]
    .iter()
    .enumerate()
    {
        let x = CORR_MARGIN_LEFT + 220.0 + idx as f64 * 66.0;
        svg.push_str(&format!(
            "<rect x=\"{x:.1}\" y=\"{:.1}\" width=\"18\" height=\"12\" fill=\"{}\" stroke=\"#cbd5e1\" stroke-width=\"1\"/>\n",
            legend_y - 10.0,
            correlation_fill_color(*value)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"9\" fill=\"#475569\">{}</text>\n",
            x + 24.0,
            legend_y,
            label
        ));
    }

    let list_top = legend_y + 22.0;
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"12\" fill=\"#0f172a\">Top synchronized pairs</text>\n",
        CORR_MARGIN_LEFT,
        list_top
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">ranked by |smoothed {}|, then |raw {}|; offset uses primary peak start positions</text>\n",
        CORR_MARGIN_LEFT,
        list_top + 14.0,
        escape_svg_text(metric.as_str()),
        escape_svg_text(metric.as_str())
    ));
    let mut ranked_rows = correlation.rows.iter().collect::<Vec<_>>();
    ranked_rows.sort_by(|left, right| {
        correlation_value(right, metric, true)
            .abs()
            .total_cmp(&correlation_value(left, metric, true).abs())
            .then(
                correlation_value(right, metric, false)
                    .abs()
                    .total_cmp(&correlation_value(left, metric, false).abs()),
            )
            .then(left.left_tf_id.cmp(&right.left_tf_id))
            .then(left.right_tf_id.cmp(&right.right_tf_id))
    });
    let mut line_y = list_top + 32.0;
    for row in ranked_rows.into_iter().take(list_rows) {
        let left_label = format_tf_track_label(&row.left_tf_id, row.left_tf_name.as_deref());
        let right_label = format_tf_track_label(&row.right_tf_id, row.right_tf_name.as_deref());
        let offset = row
            .signed_primary_peak_offset_bp
            .map(|value| format!("{value:+} bp"))
            .unwrap_or_else(|| "n/a".to_string());
        let summary = format!(
            "{} vs {} | smoothed {:+.3} | raw {:+.3} | offset {}",
            left_label,
            right_label,
            correlation_value(row, metric, true),
            correlation_value(row, metric, false),
            offset
        );
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#334155\">{}</text>\n",
            CORR_MARGIN_LEFT,
            line_y,
            escape_svg_text(&summary)
        ));
        line_y += 18.0;
    }
    svg.push_str("</svg>\n");
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{
        JasparExpertColumn, TfbsScoreTrackNormalizationReference, TfbsScoreTrackOverlayInterval,
        TfbsScoreTrackOverlayTrack, TfbsScoreTrackReport, TfbsScoreTrackRow,
        TfbsScoreTrackValueKind,
    };

    fn sample_logo_columns(len: usize) -> Vec<JasparExpertColumn> {
        (0..len)
            .map(|idx| JasparExpertColumn {
                position_1based: idx + 1,
                total_count: 10.0,
                dominant_base: if idx % 2 == 0 { "G" } else { "C" }.to_string(),
                a_count: 1.0,
                c_count: if idx % 2 == 0 { 2.0 } else { 6.0 },
                g_count: if idx % 2 == 0 { 6.0 } else { 2.0 },
                t_count: 1.0,
                a_fraction: 0.1,
                c_fraction: if idx % 2 == 0 { 0.2 } else { 0.6 },
                g_fraction: if idx % 2 == 0 { 0.6 } else { 0.2 },
                t_fraction: 0.1,
                information_content_bits: 1.4,
                a_logo_bits: 0.14,
                c_logo_bits: if idx % 2 == 0 { 0.28 } else { 0.84 },
                g_logo_bits: if idx % 2 == 0 { 0.84 } else { 0.28 },
                t_logo_bits: 0.14,
            })
            .collect()
    }

    #[test]
    fn render_tfbs_score_tracks_svg_contains_track_labels_and_axes() {
        let report = TfbsScoreTrackReport {
            schema: "gentle.tfbs_score_tracks.v1".to_string(),
            target_kind: "seq_id".to_string(),
            target_label: "tp73_upstream".to_string(),
            seq_id: "tp73_upstream".to_string(),
            source_sequence_length_bp: 4000,
            sequence_length_bp: 4000,
            scan_topology: crate::engine::InlineSequenceTopology::Linear,
            generated_at_unix_ms: 0,
            op_id: None,
            run_id: None,
            score_kind: TfbsScoreTrackValueKind::LlrBits,
            view_start_0based: 0,
            view_end_0based_exclusive: 40,
            clip_negative: true,
            motifs_requested: vec!["TP73".to_string(), "SP1".to_string()],
            global_max_score: 8.5,
            tss_markers: vec![crate::engine::TfbsScoreTrackTssMarker {
                feature_id: 7,
                feature_kind: "mRNA".to_string(),
                label: "NM_TP73".to_string(),
                position_0based: 12,
                is_reverse: false,
            }],
            overlay_tracks: vec![TfbsScoreTrackOverlayTrack {
                source_kind: "BED".to_string(),
                track_name: "cutrun_tp73".to_string(),
                display_label: "cutrun_tp73 (tp73_cutrun.bed)".to_string(),
                source_file_name: Some("tp73_cutrun.bed".to_string()),
                interval_count: 2,
                max_score: Some(42.0),
                intervals: vec![
                    TfbsScoreTrackOverlayInterval {
                        start_0based: 8,
                        end_0based_exclusive: 14,
                        label: Some("peak_a".to_string()),
                        score: Some(21.0),
                        strand: None,
                    },
                    TfbsScoreTrackOverlayInterval {
                        start_0based: 22,
                        end_0based_exclusive: 28,
                        label: Some("peak_b".to_string()),
                        score: Some(42.0),
                        strand: None,
                    },
                ],
            }],
            correlation_summary: None,
            correlation_summaries: vec![],
            cross_strand_correlation_summary: None,
            tracks: vec![
                TfbsScoreTrackRow {
                    tf_id: "MA0828.2".to_string(),
                    tf_name: Some("p73".to_string()),
                    motif_length_bp: 10,
                    motif_logo_columns: sample_logo_columns(10),
                    track_start_0based: 0,
                    scored_window_count: 4,
                    max_score: 8.5,
                    max_position_0based: Some(12),
                    normalization_reference: Some(TfbsScoreTrackNormalizationReference {
                        background_model: "uniform_random_dna".to_string(),
                        chance_model: "quantized_iid_uniform_window_dp".to_string(),
                        random_sequence_length_bp: 10_000,
                        random_seed: 7,
                        sample_count: 100,
                        mean_score: 1.2,
                        stddev_score: 0.8,
                        p95_score: 2.4,
                        p99_score: 3.1,
                        positive_fraction: 0.045,
                        observed_peak_empirical_quantile: 1.0,
                        observed_peak_modeled_quantile: 0.999999,
                        observed_peak_modeled_tail_probability: 1e-6,
                        observed_peak_modeled_tail_log10: 6.0,
                        observed_peak_delta_from_p95: 6.1,
                        observed_peak_delta_from_p99: 5.4,
                        theoretical_min_score: -4.2,
                        theoretical_max_score: 12.3,
                    }),
                    directional_summary: None,
                    top_peaks: vec![],
                    forward_scores: vec![0.0, 3.0, 8.5, 2.0],
                    reverse_scores: vec![0.0, 1.0, 2.0, 0.5],
                },
                TfbsScoreTrackRow {
                    tf_id: "MA0079.3".to_string(),
                    tf_name: Some("SP1".to_string()),
                    motif_length_bp: 9,
                    motif_logo_columns: sample_logo_columns(9),
                    track_start_0based: 0,
                    scored_window_count: 4,
                    max_score: 4.0,
                    max_position_0based: Some(8),
                    normalization_reference: Some(TfbsScoreTrackNormalizationReference {
                        background_model: "uniform_random_dna".to_string(),
                        chance_model: "quantized_iid_uniform_window_dp".to_string(),
                        random_sequence_length_bp: 10_000,
                        random_seed: 7,
                        sample_count: 100,
                        mean_score: 0.9,
                        stddev_score: 0.6,
                        p95_score: 1.8,
                        p99_score: 2.2,
                        positive_fraction: 0.022,
                        observed_peak_empirical_quantile: 0.997,
                        observed_peak_modeled_quantile: 0.9985,
                        observed_peak_modeled_tail_probability: 0.0015,
                        observed_peak_modeled_tail_log10: 2.8239,
                        observed_peak_delta_from_p95: 2.2,
                        observed_peak_delta_from_p99: 1.8,
                        theoretical_min_score: -3.1,
                        theoretical_max_score: 8.4,
                    }),
                    directional_summary: None,
                    top_peaks: vec![],
                    forward_scores: vec![0.5, 4.0, 1.0, 0.0],
                    reverse_scores: vec![0.0, 0.8, 2.0, 0.3],
                },
            ],
        };

        let svg = render_tfbs_score_tracks_svg(&report);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("Continuous TF motif score tracks"));
        assert!(svg.contains("tp73_upstream"));
        assert!(svg.contains("p73 (MA0828.2)"));
        assert!(svg.contains("SP1 (MA0079.3)"));
        assert!(svg.contains("data-gentle-role=\"tfbs-score-track-logo\""));
        assert!(svg.contains("data-gentle-role=\"tfbs-score-track-logo-letter\""));
        assert!(svg.contains("base-pair position in selected span"));
        assert!(svg.contains("#0e7490"));
        assert!(svg.contains("#b45309"));
        assert!(svg.contains("imported BED intervals = violet"));
        assert!(svg.contains("data-gentle-role=\"tfbs-score-track-overlay-row\""));
        assert!(svg.contains("data-gentle-role=\"tfbs-score-track-overlay-interval\""));
        assert!(svg.contains("cutrun_tp73 (tp73_cutrun.bed)"));
        assert!(svg.contains("data-gentle-role=\"tfbs-score-track-tss-line\""));
        assert!(svg.contains("data-gentle-role=\"tfbs-score-track-tss-arrow\""));
        assert!(svg.contains("data-gentle-role=\"tfbs-score-track-tss-label-bg\""));
        assert!(svg.contains("TSS NM_TP73"));
        assert!(svg.contains(
            "normalization=uniform_random_dna 10000bp deterministic background | chance_model=quantized_iid_uniform_window_dp"
        ));
        assert!(svg.contains("Δp99 +5.40"));
        assert!(svg.contains("bg+ 4.5%"));
    }

    #[test]
    fn render_tfbs_score_tracks_svg_supports_negative_ranges() {
        let report = TfbsScoreTrackReport {
            schema: "gentle.tfbs_score_tracks.v1".to_string(),
            target_kind: "seq_id".to_string(),
            target_label: "tp73_context".to_string(),
            seq_id: "tp73_context".to_string(),
            source_sequence_length_bp: 80,
            sequence_length_bp: 80,
            scan_topology: crate::engine::InlineSequenceTopology::Linear,
            generated_at_unix_ms: 0,
            op_id: None,
            run_id: None,
            score_kind: TfbsScoreTrackValueKind::LlrBits,
            view_start_0based: 10,
            view_end_0based_exclusive: 30,
            clip_negative: false,
            motifs_requested: vec!["REST".to_string()],
            global_max_score: 4.0,
            tss_markers: vec![],
            overlay_tracks: vec![],
            correlation_summary: None,
            correlation_summaries: vec![],
            cross_strand_correlation_summary: None,
            tracks: vec![TfbsScoreTrackRow {
                tf_id: "REST".to_string(),
                tf_name: Some("REST".to_string()),
                motif_length_bp: 8,
                motif_logo_columns: sample_logo_columns(8),
                track_start_0based: 10,
                scored_window_count: 4,
                max_score: 4.0,
                max_position_0based: Some(14),
                normalization_reference: Some(TfbsScoreTrackNormalizationReference {
                    background_model: "uniform_random_dna".to_string(),
                    chance_model: "quantized_iid_uniform_window_dp".to_string(),
                    random_sequence_length_bp: 10_000,
                    random_seed: 7,
                    sample_count: 100,
                    mean_score: -0.5,
                    stddev_score: 1.1,
                    p95_score: 2.0,
                    p99_score: 3.0,
                    positive_fraction: 0.019,
                    observed_peak_empirical_quantile: 0.991,
                    observed_peak_modeled_quantile: 0.992,
                    observed_peak_modeled_tail_probability: 0.008,
                    observed_peak_modeled_tail_log10: 2.0969,
                    observed_peak_delta_from_p95: 2.0,
                    observed_peak_delta_from_p99: 1.0,
                    theoretical_min_score: -5.0,
                    theoretical_max_score: 6.0,
                }),
                directional_summary: None,
                top_peaks: vec![],
                forward_scores: vec![-2.0, -1.0, 2.0, 4.0],
                reverse_scores: vec![-1.5, 0.0, 1.5, 2.5],
            }],
        };

        let svg = render_tfbs_score_tracks_svg(&report);
        assert!(svg.contains("score range -2.00 .. 4.00"));
        assert!(svg.contains("positive-only") == false);
        assert!(svg.contains("stroke-dasharray=\"4 3\""));
        assert!(svg.contains(">0.0</text>"));
    }

    #[test]
    fn render_tfbs_score_track_correlation_svg_contains_dual_heatmaps() {
        let report = TfbsScoreTrackReport {
            schema: "gentle.tfbs_score_tracks.v1".to_string(),
            target_kind: "seq_id".to_string(),
            target_label: "tert_promoter".to_string(),
            seq_id: "tert_promoter".to_string(),
            source_sequence_length_bp: 1201,
            sequence_length_bp: 1201,
            scan_topology: crate::engine::InlineSequenceTopology::Linear,
            generated_at_unix_ms: 0,
            op_id: None,
            run_id: None,
            score_kind: TfbsScoreTrackValueKind::LlrBackgroundTailLog10,
            view_start_0based: 0,
            view_end_0based_exclusive: 1201,
            clip_negative: true,
            motifs_requested: vec!["SP1".to_string(), "MYC".to_string(), "PATZ1".to_string()],
            global_max_score: 4.2,
            tss_markers: vec![],
            overlay_tracks: vec![],
            correlation_summary: Some(crate::engine::TfbsScoreTrackCorrelationSummary {
                signal_source: crate::engine::TfbsScoreTrackCorrelationSignalSource::MaxStrands,
                smoothing_method: "centered_boxcar".to_string(),
                smoothing_window_bp: 25,
                pair_count: 3,
                rows: vec![
                    crate::engine::TfbsScoreTrackCorrelationRow {
                        left_tf_id: "MA0079.5".to_string(),
                        left_tf_name: Some("SP1".to_string()),
                        right_tf_id: "MA0147.4".to_string(),
                        right_tf_name: Some("MYC".to_string()),
                        overlap_window_count: 100,
                        raw_pearson: 0.42,
                        smoothed_pearson: 0.81,
                        raw_spearman: 0.55,
                        smoothed_spearman: 0.88,
                        signed_primary_peak_offset_bp: Some(80),
                    },
                    crate::engine::TfbsScoreTrackCorrelationRow {
                        left_tf_id: "MA0079.5".to_string(),
                        left_tf_name: Some("SP1".to_string()),
                        right_tf_id: "MA1961.2".to_string(),
                        right_tf_name: Some("PATZ1".to_string()),
                        overlap_window_count: 100,
                        raw_pearson: 0.38,
                        smoothed_pearson: 0.96,
                        raw_spearman: 0.48,
                        smoothed_spearman: 0.97,
                        signed_primary_peak_offset_bp: Some(-152),
                    },
                    crate::engine::TfbsScoreTrackCorrelationRow {
                        left_tf_id: "MA0147.4".to_string(),
                        left_tf_name: Some("MYC".to_string()),
                        right_tf_id: "MA1961.2".to_string(),
                        right_tf_name: Some("PATZ1".to_string()),
                        overlap_window_count: 100,
                        raw_pearson: 0.21,
                        smoothed_pearson: 0.44,
                        raw_spearman: 0.27,
                        smoothed_spearman: 0.52,
                        signed_primary_peak_offset_bp: Some(-72),
                    },
                ],
            }),
            correlation_summaries: vec![crate::engine::TfbsScoreTrackCorrelationSummary {
                signal_source: crate::engine::TfbsScoreTrackCorrelationSignalSource::MaxStrands,
                smoothing_method: "centered_boxcar".to_string(),
                smoothing_window_bp: 25,
                pair_count: 3,
                rows: vec![
                    crate::engine::TfbsScoreTrackCorrelationRow {
                        left_tf_id: "MA0079.5".to_string(),
                        left_tf_name: Some("SP1".to_string()),
                        right_tf_id: "MA0147.4".to_string(),
                        right_tf_name: Some("MYC".to_string()),
                        overlap_window_count: 100,
                        raw_pearson: 0.42,
                        smoothed_pearson: 0.81,
                        raw_spearman: 0.55,
                        smoothed_spearman: 0.88,
                        signed_primary_peak_offset_bp: Some(80),
                    },
                    crate::engine::TfbsScoreTrackCorrelationRow {
                        left_tf_id: "MA0079.5".to_string(),
                        left_tf_name: Some("SP1".to_string()),
                        right_tf_id: "MA1961.2".to_string(),
                        right_tf_name: Some("PATZ1".to_string()),
                        overlap_window_count: 100,
                        raw_pearson: 0.38,
                        smoothed_pearson: 0.96,
                        raw_spearman: 0.48,
                        smoothed_spearman: 0.97,
                        signed_primary_peak_offset_bp: Some(-152),
                    },
                    crate::engine::TfbsScoreTrackCorrelationRow {
                        left_tf_id: "MA0147.4".to_string(),
                        left_tf_name: Some("MYC".to_string()),
                        right_tf_id: "MA1961.2".to_string(),
                        right_tf_name: Some("PATZ1".to_string()),
                        overlap_window_count: 100,
                        raw_pearson: 0.21,
                        smoothed_pearson: 0.44,
                        raw_spearman: 0.27,
                        smoothed_spearman: 0.52,
                        signed_primary_peak_offset_bp: Some(-72),
                    },
                ],
            }],
            cross_strand_correlation_summary: Some(
                crate::engine::TfbsScoreTrackCrossStrandCorrelationSummary {
                    smoothing_method: "centered_boxcar".to_string(),
                    smoothing_window_bp: 25,
                    pair_count: 3,
                    rows: vec![
                        crate::engine::TfbsScoreTrackCrossStrandCorrelationRow {
                            left_tf_id: "MA0079.5".to_string(),
                            left_tf_name: Some("SP1".to_string()),
                            right_tf_id: "MA0147.4".to_string(),
                            right_tf_name: Some("MYC".to_string()),
                            cells: vec![
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.41,
                                    smoothed_pearson: 0.82,
                                    raw_spearman: 0.54,
                                    smoothed_spearman: 0.87,
                                    signed_primary_peak_offset_bp: Some(80),
                                },
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.22,
                                    smoothed_pearson: 0.61,
                                    raw_spearman: 0.26,
                                    smoothed_spearman: 0.59,
                                    signed_primary_peak_offset_bp: Some(319),
                                },
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.19,
                                    smoothed_pearson: 0.52,
                                    raw_spearman: 0.07,
                                    smoothed_spearman: 0.51,
                                    signed_primary_peak_offset_bp: Some(239),
                                },
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.38,
                                    smoothed_pearson: 0.77,
                                    raw_spearman: 0.44,
                                    smoothed_spearman: 0.79,
                                    signed_primary_peak_offset_bp: Some(0),
                                },
                            ],
                        },
                        crate::engine::TfbsScoreTrackCrossStrandCorrelationRow {
                            left_tf_id: "MA0079.5".to_string(),
                            left_tf_name: Some("SP1".to_string()),
                            right_tf_id: "MA1961.2".to_string(),
                            right_tf_name: Some("PATZ1".to_string()),
                            cells: vec![
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.38,
                                    smoothed_pearson: 0.96,
                                    raw_spearman: 0.48,
                                    smoothed_spearman: 0.97,
                                    signed_primary_peak_offset_bp: Some(-152),
                                },
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.24,
                                    smoothed_pearson: 0.73,
                                    raw_spearman: 0.24,
                                    smoothed_spearman: 0.73,
                                    signed_primary_peak_offset_bp: Some(-1),
                                },
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.39,
                                    smoothed_pearson: 0.87,
                                    raw_spearman: 0.39,
                                    smoothed_spearman: 0.87,
                                    signed_primary_peak_offset_bp: Some(-72),
                                },
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.34,
                                    smoothed_pearson: 0.90,
                                    raw_spearman: 0.34,
                                    smoothed_spearman: 0.90,
                                    signed_primary_peak_offset_bp: Some(-1),
                                },
                            ],
                        },
                        crate::engine::TfbsScoreTrackCrossStrandCorrelationRow {
                            left_tf_id: "MA0147.4".to_string(),
                            left_tf_name: Some("MYC".to_string()),
                            right_tf_id: "MA1961.2".to_string(),
                            right_tf_name: Some("PATZ1".to_string()),
                            cells: vec![
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.18,
                                    smoothed_pearson: 0.44,
                                    raw_spearman: 0.21,
                                    smoothed_spearman: 0.46,
                                    signed_primary_peak_offset_bp: Some(87),
                                },
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.21,
                                    smoothed_pearson: 0.62,
                                    raw_spearman: 0.20,
                                    smoothed_spearman: 0.62,
                                    signed_primary_peak_offset_bp: Some(318),
                                },
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Forward,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.11,
                                    smoothed_pearson: 0.48,
                                    raw_spearman: 0.10,
                                    smoothed_spearman: 0.48,
                                    signed_primary_peak_offset_bp: Some(87),
                                },
                                crate::engine::TfbsScoreTrackCrossStrandCorrelationCell {
                                    left_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    right_strand:
                                        crate::engine::TfbsScoreTrackStrandComponent::Reverse,
                                    overlap_window_count: 100,
                                    raw_pearson: 0.18,
                                    smoothed_pearson: 0.41,
                                    raw_spearman: 0.19,
                                    smoothed_spearman: 0.42,
                                    signed_primary_peak_offset_bp: Some(-72),
                                },
                            ],
                        },
                    ],
                },
            ),
            tracks: vec![
                TfbsScoreTrackRow {
                    tf_id: "MA0079.5".to_string(),
                    tf_name: Some("SP1".to_string()),
                    motif_length_bp: 10,
                    motif_logo_columns: sample_logo_columns(10),
                    track_start_0based: 0,
                    scored_window_count: 8,
                    max_score: 4.0,
                    max_position_0based: Some(220),
                    normalization_reference: None,
                    directional_summary: None,
                    top_peaks: vec![],
                    forward_scores: vec![0.0; 8],
                    reverse_scores: vec![0.0; 8],
                },
                TfbsScoreTrackRow {
                    tf_id: "MA0147.4".to_string(),
                    tf_name: Some("MYC".to_string()),
                    motif_length_bp: 10,
                    motif_logo_columns: sample_logo_columns(10),
                    track_start_0based: 0,
                    scored_window_count: 8,
                    max_score: 4.0,
                    max_position_0based: Some(300),
                    normalization_reference: None,
                    directional_summary: None,
                    top_peaks: vec![],
                    forward_scores: vec![0.0; 8],
                    reverse_scores: vec![0.0; 8],
                },
                TfbsScoreTrackRow {
                    tf_id: "MA1961.2".to_string(),
                    tf_name: Some("PATZ1".to_string()),
                    motif_length_bp: 10,
                    motif_logo_columns: sample_logo_columns(10),
                    track_start_0based: 0,
                    scored_window_count: 8,
                    max_score: 4.0,
                    max_position_0based: Some(148),
                    normalization_reference: None,
                    directional_summary: None,
                    top_peaks: vec![],
                    forward_scores: vec![0.0; 8],
                    reverse_scores: vec![0.0; 8],
                },
            ],
        };

        let svg = render_tfbs_score_track_correlation_svg(
            &report,
            crate::engine::TfbsScoreTrackCorrelationMetric::Spearman,
            crate::engine::TfbsScoreTrackCorrelationSignalSource::MaxStrands,
        );
        assert!(svg.contains("TFBS track correlation"));
        assert!(svg.contains("Smoothed Spearman rho"));
        assert!(svg.contains("Raw Spearman rho"));
        assert!(svg.contains("SP1 (MA0079.5)"));
        assert!(svg.contains("MYC (MA0147.4)"));
        assert!(svg.contains("PATZ1 (MA1961.2)"));
        assert!(svg.contains("Top synchronized strand pairs"));
        assert!(svg.contains("rows and columns list each motif twice in F then R order"));
        assert!(svg.contains("SP1 (MA0079.5) F"));
        assert!(svg.contains("SP1 (MA0079.5) R"));
        assert!(svg.contains("PATZ1 (MA1961.2) F"));
        assert!(svg.contains("SP1 (MA0079.5) F vs PATZ1 (MA1961.2) F"));
        assert!(svg.contains("offset +80 bp"));
        assert!(svg.contains("offset -152 bp"));
    }
}
