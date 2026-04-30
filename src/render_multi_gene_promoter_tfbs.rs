//! Shared SVG renderer for multi-gene promoter TFBS comparisons.
//!
//! This renderer keeps promoter-set comparisons adapter-neutral: the same
//! small-multiples promoter-aligned figure can be exported from CLI, shell,
//! ClawBio, and future frontends without reimplementing TF-track plotting.

use crate::engine::MultiGenePromoterTfbsReport;

const SVG_WIDTH: f64 = 1320.0;
const SVG_MARGIN_LEFT: f64 = 36.0;
const SVG_MARGIN_RIGHT: f64 = 28.0;
const SVG_HEADER_HEIGHT: f64 = 72.0;
const SVG_FOOTER_HEIGHT: f64 = 26.0;
const GENE_HEADER_HEIGHT: f64 = 34.0;
const GENE_FOOTER_HEIGHT: f64 = 22.0;
const GENE_PANEL_GAP: f64 = 18.0;
const TRACK_ROW_HEIGHT: f64 = 56.0;
const TRACK_ROW_GAP: f64 = 8.0;
const LABEL_WIDTH: f64 = 256.0;
const PLOT_PADDING_X: f64 = 10.0;
const PLOT_PADDING_Y: f64 = 7.0;
const MAX_POLYLINE_POINTS: usize = 1400;

fn escape_svg_text(raw: &str) -> String {
    raw.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn format_tf_track_label(tf_id: &str, tf_name: Option<&str>) -> String {
    let trimmed_name = tf_name.map(str::trim).unwrap_or_default();
    if trimmed_name.is_empty() || trimmed_name.eq_ignore_ascii_case(tf_id) {
        tf_id.to_string()
    } else {
        format!("{trimmed_name} ({tf_id})")
    }
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

fn score_to_y(score: f64, min_score: f64, max_score: f64, top: f64, bottom: f64) -> f64 {
    if (max_score - min_score).abs() < f64::EPSILON {
        return bottom;
    }
    let fraction = ((score - min_score) / (max_score - min_score)).clamp(0.0, 1.0);
    bottom - fraction * (bottom - top)
}

fn sequence_position_to_x(
    position_0based: usize,
    sequence_length_bp: usize,
    left: f64,
    right: f64,
) -> f64 {
    if sequence_length_bp <= 1 {
        return left;
    }
    let width = (right - left).max(1.0);
    let fraction = position_0based.min(sequence_length_bp.saturating_sub(1)) as f64
        / sequence_length_bp.saturating_sub(1) as f64;
    left + fraction * width
}

fn polyline_points(
    values: &[f64],
    sequence_length_bp: usize,
    left: f64,
    right: f64,
    top: f64,
    bottom: f64,
    min_score: f64,
    max_score: f64,
) -> String {
    if values.is_empty() {
        return String::new();
    }
    let indices = sampled_point_indices(values.len(), MAX_POLYLINE_POINTS);
    indices
        .into_iter()
        .map(|idx| {
            let x = sequence_position_to_x(idx, sequence_length_bp, left, right);
            let y = score_to_y(values[idx], min_score, max_score, top, bottom);
            format!("{x:.2},{y:.2}")
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn global_score_bounds(report: &MultiGenePromoterTfbsReport) -> (f64, f64) {
    let mut min_score = if report.clip_negative {
        0.0
    } else {
        f64::INFINITY
    };
    let mut max_score = 0.0_f64;
    for gene in &report.genes {
        for track in &gene.tfbs_score_tracks.tracks {
            for value in track
                .forward_scores
                .iter()
                .chain(track.reverse_scores.iter())
            {
                min_score = min_score.min(*value);
                max_score = max_score.max(*value);
            }
        }
    }
    if !min_score.is_finite() {
        min_score = 0.0;
    }
    if max_score <= min_score {
        max_score = min_score + 1.0;
    }
    (min_score, max_score)
}

/// Render a promoter-aligned multi-gene TFBS small-multiples figure as SVG.
pub fn render_multi_gene_promoter_tfbs_svg(report: &MultiGenePromoterTfbsReport) -> String {
    let track_rows = report
        .genes
        .iter()
        .map(|gene| gene.tfbs_score_tracks.tracks.len().max(1))
        .sum::<usize>();
    let track_block_height = track_rows as f64 * TRACK_ROW_HEIGHT
        + track_rows.saturating_sub(report.genes.len()) as f64 * TRACK_ROW_GAP;
    let panel_header_height = report.genes.len() as f64 * GENE_HEADER_HEIGHT;
    let panel_footer_height = report.genes.len() as f64 * GENE_FOOTER_HEIGHT;
    let panel_gap_height = report.genes.len().saturating_sub(1) as f64 * GENE_PANEL_GAP;
    let svg_height = SVG_HEADER_HEIGHT
        + SVG_FOOTER_HEIGHT
        + panel_header_height
        + panel_footer_height
        + track_block_height
        + panel_gap_height;
    let content_left = SVG_MARGIN_LEFT;
    let content_right = SVG_WIDTH - SVG_MARGIN_RIGHT;
    let plot_left = content_left + LABEL_WIDTH;
    let plot_right = content_right;
    let (min_score, max_score) = global_score_bounds(report);
    let zero_visible = min_score < 0.0 && max_score > 0.0;

    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{SVG_WIDTH:.0}\" height=\"{svg_height:.0}\" viewBox=\"0 0 {SVG_WIDTH:.0} {svg_height:.0}\">\n"
    ));
    svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#fffdf8\"/>\n");
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"28\" font-family=\"monospace\" font-size=\"16\" fill=\"#111827\">{}</text>\n",
        content_left,
        escape_svg_text("Multi-gene promoter TF motif score tracks")
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"45\" font-family=\"monospace\" font-size=\"11\" fill=\"#475569\">{}</text>\n",
        content_left,
        escape_svg_text(&format!(
            "genome={} | genes={} | motifs={} | promoter windows are transcription-aligned | score={}{}",
            report.genome_id,
            report.returned_gene_count,
            report.motifs_requested.len(),
            report.score_kind.as_str(),
            if report.clip_negative && report.score_kind.supports_negative_values() {
                " | positive-only"
            } else {
                ""
            }
        ))
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"60\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">{}</text>\n",
        content_left,
        escape_svg_text("forward strand = teal | reverse strand = amber | x-axis is relative to transcript TSS")
    ));

    let mut cursor_y = SVG_HEADER_HEIGHT;
    for gene in &report.genes {
        let track_count = gene.tfbs_score_tracks.tracks.len().max(1);
        let panel_height = GENE_HEADER_HEIGHT
            + GENE_FOOTER_HEIGHT
            + track_count as f64 * TRACK_ROW_HEIGHT
            + track_count.saturating_sub(1) as f64 * TRACK_ROW_GAP;
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"8\" fill=\"#fffaf0\" stroke=\"#e2e8f0\" stroke-width=\"1\"/>\n",
            content_left,
            cursor_y,
            content_right - content_left,
            panel_height
        ));
        let gene_title = format!(
            "{} | {} | {}:{}-{} | strand={} | TSS={}",
            gene.display_label,
            gene.transcript_id,
            gene.chromosome,
            gene.promoter_start_1based,
            gene.promoter_end_1based,
            gene.strand,
            gene.tss_1based
        );
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"12\" fill=\"#0f172a\">{}</text>\n",
            content_left + 12.0,
            cursor_y + 18.0,
            escape_svg_text(&gene_title)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">{}</text>\n",
            content_left + 12.0,
            cursor_y + 31.0,
            escape_svg_text(&format!(
                "query='{}' | match={} | promoter_length={} bp | orientation={}",
                gene.gene_query,
                if gene.used_fuzzy_gene_match { "fuzzy" } else { "exact" },
                gene.promoter_length_bp,
                gene.sequence_orientation
            ))
        ));

        let tss_position_0based = gene
            .tfbs_score_tracks
            .tss_markers
            .first()
            .map(|marker| marker.position_0based)
            .unwrap_or(0);
        for (track_idx, track) in gene.tfbs_score_tracks.tracks.iter().enumerate() {
            let row_top = cursor_y
                + GENE_HEADER_HEIGHT
                + track_idx as f64 * (TRACK_ROW_HEIGHT + TRACK_ROW_GAP);
            let row_bottom = row_top + TRACK_ROW_HEIGHT;
            let plot_top = row_top + PLOT_PADDING_Y;
            let plot_bottom = row_bottom - 16.0;
            let plot_content_left = plot_left + PLOT_PADDING_X;
            let plot_content_right = plot_right - PLOT_PADDING_X;
            let row_fill = if track_idx % 2 == 0 {
                "#fffdf9"
            } else {
                "#fffefc"
            };
            svg.push_str(&format!(
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"5\" fill=\"{}\" stroke=\"#f1f5f9\" stroke-width=\"1\"/>\n",
                content_left + 6.0,
                row_top,
                content_right - content_left - 12.0,
                TRACK_ROW_HEIGHT,
                row_fill
            ));
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"11\" fill=\"#0f172a\">{}</text>\n",
                content_left + 12.0,
                row_top + 18.0,
                escape_svg_text(&format_tf_track_label(&track.tf_id, track.tf_name.as_deref()))
            ));
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"9\" fill=\"#64748b\">{}</text>\n",
                content_left + 12.0,
                row_top + 33.0,
                escape_svg_text(&format!(
                    "max {:.2} | windows {} | motif {} bp",
                    track.max_score, track.scored_window_count, track.motif_length_bp
                ))
            ));
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
            svg.push_str(&format!(
                "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#94a3b8\" stroke-width=\"1\"{} />\n",
                plot_left,
                baseline_y,
                plot_right,
                baseline_y,
                if zero_visible {
                    " stroke-dasharray=\"4 3\""
                } else {
                    ""
                }
            ));
            if track.scored_window_count == 0 {
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"11\" fill=\"#64748b\">motif longer than selected promoter</text>\n",
                    (plot_left + plot_right) * 0.5,
                    (plot_top + plot_bottom) * 0.5
                ));
                continue;
            }

            let forward_points = polyline_points(
                &track.forward_scores,
                gene.promoter_length_bp,
                plot_content_left,
                plot_content_right,
                plot_top + 3.0,
                plot_bottom - 3.0,
                min_score,
                max_score,
            );
            let reverse_points = polyline_points(
                &track.reverse_scores,
                gene.promoter_length_bp,
                plot_content_left,
                plot_content_right,
                plot_top + 3.0,
                plot_bottom - 3.0,
                min_score,
                max_score,
            );
            if !forward_points.is_empty() {
                svg.push_str(&format!(
                    "<polyline fill=\"none\" stroke=\"#0e7490\" stroke-width=\"1.7\" stroke-linejoin=\"round\" stroke-linecap=\"round\" points=\"{}\"/>\n",
                    forward_points
                ));
            }
            if !reverse_points.is_empty() {
                svg.push_str(&format!(
                    "<polyline fill=\"none\" stroke=\"#b45309\" stroke-width=\"1.5\" stroke-linejoin=\"round\" stroke-linecap=\"round\" points=\"{}\"/>\n",
                    reverse_points
                ));
            }

            if gene.promoter_length_bp > 1 && tss_position_0based < gene.promoter_length_bp {
                let tss_x = sequence_position_to_x(
                    tss_position_0based,
                    gene.promoter_length_bp,
                    plot_content_left,
                    plot_content_right,
                );
                svg.push_str(&format!(
                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"#7c3aed\" stroke-width=\"1\" stroke-dasharray=\"3 3\"/>\n",
                    tss_x, plot_top, tss_x, plot_bottom
                ));
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"9\" fill=\"#7c3aed\">TSS</text>\n",
                    tss_x,
                    plot_top - 2.0
                ));
            }

            if track_idx + 1 == gene.tfbs_score_tracks.tracks.len() {
                let relative_start = -(tss_position_0based as i64);
                let relative_end =
                    gene.promoter_length_bp
                        .saturating_sub(tss_position_0based + 1) as i64;
                let tss_x = sequence_position_to_x(
                    tss_position_0based,
                    gene.promoter_length_bp,
                    plot_content_left,
                    plot_content_right,
                );
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"9\" fill=\"#64748b\">{}</text>\n",
                    plot_content_left,
                    row_bottom - 2.0,
                    escape_svg_text(&format!("{relative_start:+} bp"))
                ));
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"9\" fill=\"#64748b\">0 bp</text>\n",
                    tss_x,
                    row_bottom - 2.0
                ));
                svg.push_str(&format!(
                    "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"9\" fill=\"#64748b\">{}</text>\n",
                    plot_content_right,
                    row_bottom - 2.0,
                    escape_svg_text(&format!("{relative_end:+} bp"))
                ));
            }
        }

        cursor_y += panel_height + GENE_PANEL_GAP;
    }

    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#64748b\">{}</text>\n",
        content_left,
        svg_height - 10.0,
        escape_svg_text(&format!(
            "summary rows={} | warnings={}",
            report.summary_rows.len(),
            report.warnings.len()
        ))
    ));
    svg.push_str("</svg>\n");
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{
        InlineSequenceTopology, MultiGenePromoterTfbsGeneReport, MultiGenePromoterTfbsSummaryRow,
        TfbsScoreTrackReport, TfbsScoreTrackRow, TfbsScoreTrackTssMarker, TfbsScoreTrackValueKind,
    };

    #[test]
    fn render_multi_gene_promoter_tfbs_svg_contains_gene_labels_and_tss() {
        let report = MultiGenePromoterTfbsReport {
            schema: "gentle.multi_gene_promoter_tfbs.v1".to_string(),
            generated_at_unix_ms: 0,
            op_id: None,
            run_id: None,
            genome_id: "ToyGenome".to_string(),
            upstream_bp: 100,
            downstream_bp: 20,
            score_kind: TfbsScoreTrackValueKind::LlrBackgroundTailLog10,
            clip_negative: true,
            motifs_requested: vec!["SP1".to_string()],
            gene_queries_requested: vec![],
            returned_gene_count: 1,
            genes: vec![MultiGenePromoterTfbsGeneReport {
                gene_query: "TERT".to_string(),
                occurrence: 1,
                transcript_id_requested: None,
                display_label: "TERT".to_string(),
                gene_id: Some("ENSG_TERT".to_string()),
                gene_name: Some("TERT".to_string()),
                transcript_id: "ENST_TERT".to_string(),
                chromosome: "5".to_string(),
                strand: "+".to_string(),
                promoter_start_1based: 1000,
                promoter_end_1based: 1120,
                promoter_length_bp: 121,
                tss_1based: 1100,
                sequence_orientation: "transcription_aligned".to_string(),
                used_fuzzy_gene_match: false,
                tfbs_score_tracks: TfbsScoreTrackReport {
                    schema: "gentle.tfbs_score_tracks.v1".to_string(),
                    target_kind: "inline_sequence".to_string(),
                    target_label: "TERT".to_string(),
                    seq_id: "TERT".to_string(),
                    source_sequence_length_bp: 121,
                    sequence_length_bp: 121,
                    scan_topology: InlineSequenceTopology::Linear,
                    generated_at_unix_ms: 0,
                    op_id: None,
                    run_id: None,
                    score_kind: TfbsScoreTrackValueKind::LlrBackgroundTailLog10,
                    view_start_0based: 0,
                    view_end_0based_exclusive: 121,
                    clip_negative: true,
                    motifs_requested: vec!["SP1".to_string()],
                    global_max_score: 4.0,
                    tss_markers: vec![TfbsScoreTrackTssMarker {
                        feature_id: usize::MAX,
                        feature_kind: "genome_promoter_slice".to_string(),
                        label: "ENST_TERT".to_string(),
                        position_0based: 100,
                        is_reverse: false,
                    }],
                    overlay_tracks: vec![],
                    correlation_summary: None,
                    correlation_summaries: vec![],
                    cross_strand_correlation_summary: None,
                    tracks: vec![TfbsScoreTrackRow {
                        tf_id: "SP1".to_string(),
                        tf_name: Some("SP1".to_string()),
                        motif_length_bp: 10,
                        motif_logo_columns: vec![],
                        track_start_0based: 0,
                        scored_window_count: 4,
                        max_score: 4.0,
                        max_position_0based: Some(2),
                        normalization_reference: None,
                        directional_summary: None,
                        top_peaks: vec![],
                        forward_scores: vec![0.0, 2.0, 4.0, 1.0],
                        reverse_scores: vec![0.0, 1.0, 0.5, 0.0],
                    }],
                },
            }],
            summary_rows: vec![MultiGenePromoterTfbsSummaryRow {
                gene_label: "TERT".to_string(),
                gene_query: "TERT".to_string(),
                transcript_id: "ENST_TERT".to_string(),
                tf_id: "SP1".to_string(),
                tf_name: Some("SP1".to_string()),
                max_score: 4.0,
                peak_position_0based: Some(2),
                peak_position_promoter_relative_bp: Some(-98),
                peak_genomic_position_1based: Some(1002),
                positive_fraction: 0.75,
            }],
            warnings: vec![],
        };
        let svg = render_multi_gene_promoter_tfbs_svg(&report);
        let plot_left = SVG_MARGIN_LEFT + LABEL_WIDTH;
        let plot_right = SVG_WIDTH - SVG_MARGIN_RIGHT;
        let expected_zero_x = sequence_position_to_x(
            100,
            121,
            plot_left + PLOT_PADDING_X,
            plot_right - PLOT_PADDING_X,
        );
        let row_bottom = SVG_HEADER_HEIGHT + GENE_HEADER_HEIGHT + TRACK_ROW_HEIGHT;
        assert!(svg.contains("Multi-gene promoter TF motif score tracks"));
        assert!(svg.contains("TERT"));
        assert!(svg.contains("TSS"));
        assert!(svg.contains("SP1"));
        assert!(svg.contains(&format!(
            "<text x=\"{expected_zero_x:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"9\" fill=\"#64748b\">0 bp</text>",
            row_bottom - 2.0
        )));
    }
}
