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
const SVG_ROW_HEIGHT: f64 = 88.0;
const SVG_ROW_GAP: f64 = 14.0;
const SVG_LABEL_WIDTH: f64 = 260.0;
const SVG_TRACK_PADDING_X: f64 = 10.0;
const SVG_TRACK_PADDING_Y: f64 = 8.0;
const SVG_MAX_POINTS_PER_POLYLINE: usize = 1400;
const SVG_TSS_MARKER_ROW_HEIGHT: f64 = 18.0;
const SVG_TSS_MARKER_MAX_LANES: usize = 3;

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
            "<line x1=\"{x:.2}\" y1=\"{line_top:.2}\" x2=\"{x:.2}\" y2=\"{line_bottom:.2}\" stroke=\"{color}\" stroke-width=\"1.25\" stroke-dasharray=\"6 4\" opacity=\"0.85\" data-gentle-role=\"tfbs-score-track-tss-line\"/>\n"
        ));
        svg.push_str(&format!(
            "<path d=\"M{x:.2},{line_top:.2} L{x:.2},{stem_y:.2} Q{x:.2},{hook_y:.2} {elbow_x:.2},{hook_y:.2} L{tip_x:.2},{hook_y:.2}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"2\" stroke-linecap=\"round\" stroke-linejoin=\"round\" data-gentle-role=\"tfbs-score-track-tss-arrow\"/>\n"
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
    let svg_height = header_height
        + SVG_FOOTER_HEIGHT
        + track_count as f64 * SVG_ROW_HEIGHT
        + track_count.saturating_sub(1) as f64 * SVG_ROW_GAP;
    let content_left = SVG_MARGIN_LEFT;
    let content_right = SVG_WIDTH - SVG_MARGIN_RIGHT;
    let plot_left = content_left + SVG_LABEL_WIDTH;
    let plot_right = content_right;
    let (min_score, max_score) = global_score_bounds(report);
    let zero_visible = min_score < 0.0 && max_score > 0.0;
    let (bp_start_text, bp_mid_text, bp_end_text) =
        format_bp_ticks(report.view_start_0based, report.view_end_0based_exclusive);
    let title = "Continuous TF motif score tracks";
    let subtitle = format!(
        "seq={} | span={}..{} | motifs={} | score={}{}",
        report.seq_id,
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
    let legend = "forward strand = teal | reverse strand = amber";
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
        let meta = format!(
            "{} windows | max {:.2}{}",
            track.scored_window_count, track.max_score, max_position
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
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"13\" fill=\"#0f172a\">{}</text>\n",
            content_left + 12.0,
            row_top + 24.0,
            escape_svg_text(&label)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">{}</text>\n",
            content_left + 12.0,
            row_top + 40.0,
            escape_svg_text(&meta)
        ));
        if let Some(normalization_meta) = normalization_meta {
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">{}</text>\n",
                content_left + 12.0,
                row_top + 56.0,
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

    if tss_marker_lanes > 0 {
        let stack_bottom = header_height
            + report.tracks.len().max(1) as f64 * (SVG_ROW_HEIGHT + SVG_ROW_GAP)
            - SVG_ROW_GAP
            - 10.0;
        svg.push_str(&render_tss_marker_annotations(
            &summarized_tss_markers,
            report.view_start_0based,
            report.view_end_0based_exclusive,
            plot_left + SVG_TRACK_PADDING_X,
            plot_right - SVG_TRACK_PADDING_X,
            SVG_HEADER_HEIGHT - 2.0,
            header_height - 6.0,
            stack_bottom,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{
        TfbsScoreTrackNormalizationReference, TfbsScoreTrackReport, TfbsScoreTrackRow,
        TfbsScoreTrackValueKind,
    };

    #[test]
    fn render_tfbs_score_tracks_svg_contains_track_labels_and_axes() {
        let report = TfbsScoreTrackReport {
            schema: "gentle.tfbs_score_tracks.v1".to_string(),
            seq_id: "tp73_upstream".to_string(),
            sequence_length_bp: 4000,
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
            correlation_summary: None,
            tracks: vec![
                TfbsScoreTrackRow {
                    tf_id: "MA0828.2".to_string(),
                    tf_name: Some("p73".to_string()),
                    motif_length_bp: 10,
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
                    top_peaks: vec![],
                    forward_scores: vec![0.0, 3.0, 8.5, 2.0],
                    reverse_scores: vec![0.0, 1.0, 2.0, 0.5],
                },
                TfbsScoreTrackRow {
                    tf_id: "MA0079.3".to_string(),
                    tf_name: Some("SP1".to_string()),
                    motif_length_bp: 9,
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
        assert!(svg.contains("base-pair position in selected span"));
        assert!(svg.contains("#0e7490"));
        assert!(svg.contains("#b45309"));
        assert!(svg.contains("data-gentle-role=\"tfbs-score-track-tss-line\""));
        assert!(svg.contains("data-gentle-role=\"tfbs-score-track-tss-arrow\""));
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
            seq_id: "tp73_context".to_string(),
            sequence_length_bp: 80,
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
            correlation_summary: None,
            tracks: vec![TfbsScoreTrackRow {
                tf_id: "REST".to_string(),
                tf_name: Some("REST".to_string()),
                motif_length_bp: 8,
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
}
