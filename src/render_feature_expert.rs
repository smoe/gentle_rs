//! Feature expert-view SVG renderer.

use crate::feature_expert::{
    FeatureExpertView, RestrictionSiteExpertView, SplicingExpertView, TfbsExpertView,
};
use std::collections::HashMap;
use svg::Document;
use svg::node::element::path::Data;
use svg::node::element::{Circle, Line, Path, Rectangle, Text};

const W: f32 = 1200.0;
const H: f32 = 700.0;

fn nuc_color(base_idx: usize) -> &'static str {
    match base_idx {
        0 => "#2b8a3e", // A
        1 => "#1d4ed8", // C
        2 => "#d97706", // G
        3 => "#b91c1c", // T
        _ => "#555555",
    }
}

fn wrap_text(input: &str, max_chars: usize) -> Vec<String> {
    let mut out = Vec::new();
    let mut current = String::new();
    for token in input.split_whitespace() {
        let next_len = if current.is_empty() {
            token.len()
        } else {
            current.len() + 1 + token.len()
        };
        if !current.is_empty() && next_len > max_chars {
            out.push(current);
            current = token.to_string();
        } else {
            if !current.is_empty() {
                current.push(' ');
            }
            current.push_str(token);
        }
    }
    if !current.is_empty() {
        out.push(current);
    }
    if out.is_empty() {
        out.push(String::new());
    }
    out
}

fn match_base_to_idx(base: Option<char>) -> Option<usize> {
    match base.map(|c| c.to_ascii_uppercase()) {
        Some('A') => Some(0),
        Some('C') => Some(1),
        Some('G') => Some(2),
        Some('T') => Some(3),
        _ => None,
    }
}

fn render_tfbs(view: &TfbsExpertView) -> String {
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

    let left = 84.0;
    let right = W - 54.0;
    let chart_top = 120.0;
    let chart_bottom = H - 190.0;
    let chart_height = chart_bottom - chart_top;
    let baseline = chart_bottom;
    let n = view.columns.len().max(1);
    let available = (right - left).max(1.0);
    let step = (available / n as f32).max(8.0);
    let bar_width = (step * 0.72).clamp(5.0, 22.0);
    let scale = chart_height / 2.0;

    doc = doc.add(
        Text::new(format!(
            "TFBS expert: {} ({}) on {}:{}..{} ({})",
            view.tf_name.clone().unwrap_or_else(|| view.tf_id.clone()),
            view.tf_id,
            view.seq_id,
            view.start_1based,
            view.end_1based,
            view.strand
        ))
        .set("x", left)
        .set("y", 36)
        .set("font-family", "monospace")
        .set("font-size", 18)
        .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new(format!(
            "matched={} | motif_length={} | llr_bits={:.4} | true_log_odds_bits={:.4}",
            view.matched_sequence,
            view.motif_length,
            view.llr_total_bits.unwrap_or_default(),
            view.true_log_odds_total_bits.unwrap_or_default()
        ))
        .set("x", left)
        .set("y", 58)
        .set("font-family", "monospace")
        .set("font-size", 13)
        .set("fill", "#444444"),
    );

    doc = doc.add(
        Line::new()
            .set("x1", left)
            .set("y1", baseline)
            .set("x2", right)
            .set("y2", baseline)
            .set("stroke", "#333333")
            .set("stroke-width", 1),
    );
    for tick in [0.0_f32, 1.0, 2.0] {
        let y = baseline - tick * scale;
        doc = doc.add(
            Line::new()
                .set("x1", left - 6.0)
                .set("y1", y)
                .set("x2", right)
                .set("y2", y)
                .set("stroke", if tick == 0.0 { "#333333" } else { "#e5e7eb" })
                .set("stroke-width", 1),
        );
        doc = doc.add(
            Text::new(format!("{tick:.0}"))
                .set("x", left - 10.0)
                .set("y", y + 4.0)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#555555"),
        );
    }
    doc = doc.add(
        Text::new("information content (bits)")
            .set("x", left - 48.0)
            .set("y", chart_top - 10.0)
            .set("font-family", "monospace")
            .set("font-size", 11)
            .set("fill", "#555555"),
    );

    let mut match_points: Vec<(f32, f32)> = Vec::new();

    for (idx, col) in view.columns.iter().enumerate() {
        let x_center = left + step * (idx as f32 + 0.5);
        let x = x_center - bar_width / 2.0;
        let ic = col.information_content_bits.clamp(0.0, 2.0) as f32;
        let full_height = ic * scale;
        let mut y_cursor = baseline;
        let match_idx = match_base_to_idx(col.match_base);
        let mut match_y: Option<f32> = None;

        for base_idx in 0..4 {
            let frac = col.frequencies[base_idx].clamp(0.0, 1.0) as f32;
            let segment_h = (full_height * frac).max(0.0);
            if segment_h <= 0.0 {
                continue;
            }
            let y_top = y_cursor - segment_h;
            doc = doc.add(
                Rectangle::new()
                    .set("x", x)
                    .set("y", y_top)
                    .set("width", bar_width)
                    .set("height", segment_h)
                    .set("fill", nuc_color(base_idx)),
            );
            if match_idx == Some(base_idx) {
                match_y = Some(y_top + segment_h * 0.5);
            }
            y_cursor = y_top;
        }

        if let Some(y) = match_y {
            match_points.push((x_center, y));
        }

        if idx % 2 == 0 || n <= 24 {
            doc = doc.add(
                Text::new(format!("{}", col.index_1based))
                    .set("x", x_center)
                    .set("y", baseline + 16.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 10)
                    .set("fill", "#4b5563"),
            );
        }

        if let Some(base) = col.match_base {
            doc = doc.add(
                Text::new(base.to_ascii_uppercase().to_string())
                    .set("x", x_center)
                    .set("y", baseline + 32.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 11)
                    .set("fill", "#111111"),
            );
        }
    }

    if match_points.len() >= 2 {
        let mut data = Data::new().move_to((match_points[0].0, match_points[0].1));
        for (x, y) in match_points.iter().skip(1) {
            data = data.line_to((*x, *y));
        }
        doc = doc.add(
            Path::new()
                .set("d", data)
                .set("fill", "none")
                .set("stroke", "#111111")
                .set("stroke-width", 2.0),
        );
        for (x, y) in match_points {
            doc = doc.add(
                Circle::new()
                    .set("cx", x)
                    .set("cy", y)
                    .set("r", 2.3)
                    .set("fill", "#111111"),
            );
        }
    }

    let legend_items = [("A", 0usize), ("C", 1usize), ("G", 2usize), ("T", 3usize)];
    let mut lx = left;
    let ly = H - 132.0;
    for (label, idx) in legend_items {
        doc = doc.add(
            Rectangle::new()
                .set("x", lx)
                .set("y", ly - 9.0)
                .set("width", 12)
                .set("height", 12)
                .set("fill", nuc_color(idx)),
        );
        doc = doc.add(
            Text::new(label)
                .set("x", lx + 16.0)
                .set("y", ly + 1.0)
                .set("font-family", "monospace")
                .set("font-size", 12)
                .set("fill", "#111111"),
        );
        lx += 52.0;
    }
    doc = doc.add(
        Text::new("black line = matched sequence path across motif columns")
            .set("x", left + 230.0)
            .set("y", ly + 1.0)
            .set("font-family", "monospace")
            .set("font-size", 11)
            .set("fill", "#4b5563"),
    );

    for (line_idx, line) in wrap_text(&view.instruction, 125).into_iter().enumerate() {
        doc = doc.add(
            Text::new(line)
                .set("x", left)
                .set("y", H - 86.0 + line_idx as f32 * 14.0)
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#374151"),
        );
    }

    doc.to_string()
}

fn render_restriction(view: &RestrictionSiteExpertView) -> String {
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

    let title = format!(
        "Restriction expert: {} @ {}:{}..{}",
        view.selected_enzyme
            .clone()
            .unwrap_or_else(|| view.enzyme_names.join(",")),
        view.seq_id,
        view.recognition_start_1based,
        view.recognition_end_1based
    );
    doc = doc.add(
        Text::new(title)
            .set("x", 90)
            .set("y", 42)
            .set("font-family", "monospace")
            .set("font-size", 18)
            .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new(format!(
            "cut_pos={} | cuts_for_enzyme={} | recognition_iupac={}",
            view.cut_pos_1based,
            view.number_of_cuts_for_enzyme,
            view.recognition_iupac
                .clone()
                .unwrap_or_else(|| "-".to_string())
        ))
        .set("x", 90)
        .set("y", 66)
        .set("font-family", "monospace")
        .set("font-size", 12)
        .set("fill", "#4b5563"),
    );

    let top = if view.site_sequence.is_empty() {
        view.recognition_iupac.clone().unwrap_or_default()
    } else {
        view.site_sequence.clone()
    };
    let bottom = if view.site_sequence_complement.is_empty() {
        String::new()
    } else {
        view.site_sequence_complement.clone()
    };

    let n = top.chars().count().max(1);
    let step = ((W - 220.0) / n as f32).clamp(16.0, 28.0);
    let total = step * n as f32;
    let start_x = (W - total) * 0.5;
    let top_y = 260.0;
    let bottom_y = 320.0;
    let rail_left = start_x - 18.0;
    let rail_right = start_x + total + 18.0;

    doc = doc.add(
        Line::new()
            .set("x1", rail_left)
            .set("y1", top_y - 10.0)
            .set("x2", rail_right)
            .set("y2", top_y - 10.0)
            .set("stroke", "#111111")
            .set("stroke-width", 2),
    );
    doc = doc.add(
        Line::new()
            .set("x1", rail_left)
            .set("y1", bottom_y + 10.0)
            .set("x2", rail_right)
            .set("y2", bottom_y + 10.0)
            .set("stroke", "#111111")
            .set("stroke-width", 2),
    );
    doc = doc.add(
        Text::new("5'")
            .set("x", rail_left - 32.0)
            .set("y", top_y - 6.0)
            .set("font-family", "monospace")
            .set("font-size", 14)
            .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new("3'")
            .set("x", rail_right + 10.0)
            .set("y", top_y - 6.0)
            .set("font-family", "monospace")
            .set("font-size", 14)
            .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new("3'")
            .set("x", rail_left - 32.0)
            .set("y", bottom_y + 14.0)
            .set("font-family", "monospace")
            .set("font-size", 14)
            .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new("5'")
            .set("x", rail_right + 10.0)
            .set("y", bottom_y + 14.0)
            .set("font-family", "monospace")
            .set("font-size", 14)
            .set("fill", "#111111"),
    );

    for (idx, ch) in top.chars().enumerate() {
        let x = start_x + step * (idx as f32 + 0.5);
        doc = doc.add(
            Text::new(ch.to_string())
                .set("x", x)
                .set("y", top_y)
                .set("text-anchor", "middle")
                .set("font-family", "monospace")
                .set("font-size", 22)
                .set("fill", "#0f172a"),
        );
    }
    for (idx, ch) in bottom.chars().enumerate() {
        let x = start_x + step * (idx as f32 + 0.5);
        doc = doc.add(
            Text::new(ch.to_string())
                .set("x", x)
                .set("y", bottom_y)
                .set("text-anchor", "middle")
                .set("font-family", "monospace")
                .set("font-size", 22)
                .set("fill", "#334155"),
        );
    }

    let cut_x = start_x + step * view.cut_index_0based as f32;
    doc = doc.add(
        Line::new()
            .set("x1", cut_x)
            .set("y1", top_y - 42.0)
            .set("x2", cut_x)
            .set("y2", bottom_y + 32.0)
            .set("stroke", "#b91c1c")
            .set("stroke-width", 2.5),
    );
    doc = doc.add(
        Text::new("cut")
            .set("x", cut_x + 4.0)
            .set("y", top_y - 46.0)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#b91c1c"),
    );

    if let Some(overlap) = view.overlap_bp {
        doc = doc.add(
            Text::new(format!("overhang_bp={overlap}"))
                .set("x", 90)
                .set("y", 402)
                .set("font-family", "monospace")
                .set("font-size", 12)
                .set("fill", "#4b5563"),
        );
    }

    for (line_idx, line) in wrap_text(&view.instruction, 125).into_iter().enumerate() {
        doc = doc.add(
            Text::new(line)
                .set("x", 90)
                .set("y", H - 84.0 + line_idx as f32 * 14.0)
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#374151"),
        );
    }

    doc.to_string()
}

fn render_splicing(view: &SplicingExpertView) -> String {
    let lane_count = view.transcripts.len().max(1);
    let exon_count = view.unique_exons.len().max(1);
    let lane_height = 30.0_f32;
    let lane_gap = 14.0_f32;
    let chart_top = 126.0_f32;
    let chart_height =
        lane_count as f32 * lane_height + (lane_count.saturating_sub(1) as f32) * lane_gap;
    let matrix_top = chart_top + chart_height + 78.0;
    let matrix_row_h = 17.0_f32;
    let matrix_h = lane_count as f32 * matrix_row_h + 28.0;
    let footer_top = matrix_top + matrix_h + 36.0;
    let dyn_h = (footer_top + 126.0).max(H);
    let left = 240.0_f32;
    let right = W - 56.0_f32;
    let map_width = (right - left).max(1.0);
    let region_start = view.region_start_1based.max(1);
    let region_end = view.region_end_1based.max(region_start);
    let region_span = (region_end - region_start + 1).max(1) as f32;
    let x_for = |pos_1based: usize| -> f32 {
        let clamped = pos_1based.clamp(region_start, region_end);
        let rel = (clamped.saturating_sub(region_start)) as f32 / region_span;
        left + rel * map_width
    };

    let exon_meta: HashMap<(usize, usize), (usize, bool)> = view
        .unique_exons
        .iter()
        .map(|exon| {
            (
                (exon.start_1based, exon.end_1based),
                (exon.support_transcript_count, exon.constitutive),
            )
        })
        .collect();

    let mut doc = Document::new()
        .set("viewBox", (0, 0, W, dyn_h))
        .set("width", W)
        .set("height", dyn_h)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", W)
                .set("height", dyn_h)
                .set("fill", "#ffffff"),
        );

    doc = doc
        .add(
            Text::new(format!(
                "Splicing expert: {} ({}) on {}:{}..{} ({})",
                view.group_label,
                view.transcript_count,
                view.seq_id,
                view.region_start_1based,
                view.region_end_1based,
                view.strand
            ))
            .set("x", 88)
            .set("y", 38)
            .set("font-family", "monospace")
            .set("font-size", 18)
            .set("fill", "#111111"),
        )
        .add(
            Text::new(format!(
                "target_feature={} | unique_exons={} | donor/acceptor markers shown",
                view.target_feature_id, view.unique_exon_count
            ))
            .set("x", 88)
            .set("y", 60)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#4b5563"),
        )
        .add(
            Line::new()
                .set("x1", left)
                .set("y1", chart_top - 16.0)
                .set("x2", right)
                .set("y2", chart_top - 16.0)
                .set("stroke", "#6b7280")
                .set("stroke-width", 1.0),
        )
        .add(
            Text::new(format!("{} bp", region_start))
                .set("x", left)
                .set("y", chart_top - 21.0)
                .set("text-anchor", "start")
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#4b5563"),
        )
        .add(
            Text::new(format!("{} bp", region_end))
                .set("x", right)
                .set("y", chart_top - 21.0)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#4b5563"),
        );

    for junction in &view.junctions {
        let x1 = x_for(junction.donor_1based);
        let x2 = x_for(junction.acceptor_1based);
        if (x2 - x1).abs() < 2.0 {
            continue;
        }
        let mid_x = (x1 + x2) * 0.5;
        let span = (x2 - x1).abs().max(1.0);
        let lift = (22.0 + span * 0.04).clamp(24.0, 84.0);
        let apex_y = chart_top - lift;
        let stroke_w = (1.0 + junction.support_transcript_count as f32 * 0.35).clamp(1.2, 4.4);
        let path = Data::new()
            .move_to((x1, chart_top - 2.0))
            .quadratic_curve_to((mid_x, apex_y, x2, chart_top - 2.0));
        doc = doc.add(
            Path::new()
                .set("d", path)
                .set("fill", "none")
                .set("stroke", "#64748b")
                .set("stroke-width", stroke_w),
        );
    }

    for (lane_idx, lane) in view.transcripts.iter().enumerate() {
        let y = chart_top + lane_idx as f32 * (lane_height + lane_gap) + lane_height * 0.5;
        if lane.has_target_feature {
            doc = doc.add(
                Rectangle::new()
                    .set("x", left - 4.0)
                    .set("y", y - lane_height * 0.5 - 3.0)
                    .set("width", map_width + 8.0)
                    .set("height", lane_height + 6.0)
                    .set("fill", "#fef9c3")
                    .set("fill-opacity", 0.35),
            );
        }
        doc = doc.add(
            Text::new(format!(
                "{} ({})",
                lane.transcript_id,
                if lane.strand.trim().is_empty() {
                    "?"
                } else {
                    lane.strand.as_str()
                }
            ))
            .set("x", left - 10.0)
            .set("y", y + 4.0)
            .set("text-anchor", "end")
            .set("font-family", "monospace")
            .set("font-size", 11)
            .set(
                "fill",
                if lane.has_target_feature {
                    "#854d0e"
                } else {
                    "#111827"
                },
            ),
        );

        for intron in &lane.introns {
            let x1 = x_for(intron.start_1based);
            let x2 = x_for(intron.end_1based);
            if x2 <= x1 {
                continue;
            }
            doc = doc.add(
                Line::new()
                    .set("x1", x1)
                    .set("y1", y)
                    .set("x2", x2)
                    .set("y2", y)
                    .set("stroke", "#6b7280")
                    .set("stroke-width", 1.2),
            );
        }

        for exon in &lane.exons {
            let x1 = x_for(exon.start_1based);
            let x2 = x_for(exon.end_1based);
            if x2 <= x1 {
                continue;
            }
            let (support, constitutive) = exon_meta
                .get(&(exon.start_1based, exon.end_1based))
                .copied()
                .unwrap_or((1, false));
            let fill = if constitutive { "#2563eb" } else { "#b45309" };
            let opacity = if constitutive {
                0.92
            } else if support <= 1 {
                0.65
            } else {
                0.78
            };
            doc = doc
                .add(
                    Rectangle::new()
                        .set("x", x1)
                        .set("y", y - 9.5)
                        .set("width", (x2 - x1).max(1.0))
                        .set("height", 19.0)
                        .set("fill", fill)
                        .set("fill-opacity", opacity)
                        .set("stroke", "#1f2937")
                        .set("stroke-width", 0.6),
                )
                .add(
                    Text::new(format!("{support}"))
                        .set("x", x1 + 2.0)
                        .set("y", y - 11.0)
                        .set("font-family", "monospace")
                        .set("font-size", 8)
                        .set("fill", "#334155"),
                );
        }
    }

    for marker in &view.boundaries {
        let lane_idx = view
            .transcripts
            .iter()
            .position(|lane| lane.transcript_feature_id == marker.transcript_feature_id);
        let Some(lane_idx) = lane_idx else {
            continue;
        };
        let y = chart_top + lane_idx as f32 * (lane_height + lane_gap) + lane_height * 0.5;
        let x = x_for(marker.position_1based);
        let (tick_top, tick_bottom, color, dy) = if marker.side.eq_ignore_ascii_case("donor") {
            (y - 18.0, y - 8.0, "#be123c", -20.0)
        } else {
            (y + 8.0, y + 18.0, "#0f766e", 28.0)
        };
        doc = doc
            .add(
                Line::new()
                    .set("x1", x)
                    .set("y1", tick_top)
                    .set("x2", x)
                    .set("y2", tick_bottom)
                    .set("stroke", color)
                    .set("stroke-width", 1.4),
            )
            .add(
                Text::new(format!(
                    "{}:{}{}",
                    marker.side,
                    marker.motif_2bp,
                    if marker.canonical { "" } else { "*" }
                ))
                .set("x", x + 2.0)
                .set("y", y + dy)
                .set("font-family", "monospace")
                .set("font-size", 8)
                .set("fill", color),
            );
    }

    let matrix_label_x = 88.0_f32;
    let matrix_left = 250.0_f32;
    let matrix_cell_w = ((right - matrix_left - 6.0) / exon_count as f32).clamp(7.0, 16.0);
    let matrix_cell_h = 12.0_f32;
    doc = doc.add(
        Text::new("Transcript vs exon matrix")
            .set("x", matrix_label_x)
            .set("y", matrix_top - 14.0)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#111827"),
    );
    for (col_idx, _) in view.unique_exons.iter().enumerate() {
        let x = matrix_left + col_idx as f32 * matrix_cell_w + matrix_cell_w * 0.5;
        doc = doc.add(
            Text::new(format!("{}", col_idx + 1))
                .set("x", x)
                .set("y", matrix_top - 2.0)
                .set("text-anchor", "middle")
                .set("font-family", "monospace")
                .set("font-size", 8)
                .set("fill", "#475569"),
        );
    }
    for (row_idx, row) in view.matrix_rows.iter().enumerate() {
        let y = matrix_top + row_idx as f32 * matrix_row_h;
        doc = doc.add(
            Text::new(row.transcript_id.clone())
                .set("x", matrix_label_x)
                .set("y", y + 10.0)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#1f2937"),
        );
        for (col_idx, present) in row.exon_presence.iter().copied().enumerate() {
            let x = matrix_left + col_idx as f32 * matrix_cell_w;
            doc = doc.add(
                Rectangle::new()
                    .set("x", x)
                    .set("y", y)
                    .set("width", matrix_cell_w - 1.0)
                    .set("height", matrix_cell_h)
                    .set("fill", if present { "#334155" } else { "#f1f5f9" })
                    .set("stroke", "#cbd5e1")
                    .set("stroke-width", 0.5),
            );
        }
    }

    let mut event_y = footer_top;
    doc = doc.add(
        Text::new("Event summary")
            .set("x", 88.0)
            .set("y", event_y)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#111827"),
    );
    event_y += 14.0;
    for event in &view.events {
        doc = doc.add(
            Text::new(format!("{}: {}", event.event_type, event.count))
                .set("x", 88.0)
                .set("y", event_y)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#334155"),
        );
        event_y += 12.0;
        for detail in event.details.iter().take(2) {
            doc = doc.add(
                Text::new(format!("  - {detail}"))
                    .set("x", 88.0)
                    .set("y", event_y)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#475569"),
            );
            event_y += 11.0;
        }
    }

    for (line_idx, line) in wrap_text(&view.instruction, 140).into_iter().enumerate() {
        doc = doc.add(
            Text::new(line)
                .set("x", 88)
                .set("y", dyn_h - 64.0 + line_idx as f32 * 14.0)
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#374151"),
        );
    }

    doc.to_string()
}

pub fn render_feature_expert_svg(view: &FeatureExpertView) -> String {
    match view {
        FeatureExpertView::Tfbs(tfbs) => render_tfbs(tfbs),
        FeatureExpertView::RestrictionSite(re) => render_restriction(re),
        FeatureExpertView::Splicing(splicing) => render_splicing(splicing),
    }
}
