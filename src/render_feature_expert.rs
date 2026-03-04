//! Feature expert-view SVG renderer.

use crate::feature_expert::{
    FeatureExpertView, IsoformArchitectureExpertView, RestrictionSiteExpertView,
    SplicingExpertView, TfbsExpertView,
};
use std::collections::{BTreeMap, BTreeSet, HashMap};
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

fn render_isoform_architecture(view: &IsoformArchitectureExpertView) -> String {
    let lane_count = view.transcript_lanes.len().max(1);
    let lane_h = 24.0_f32;
    let lane_gap = 10.0_f32;
    let top_header = 114.0_f32;
    let exon_chart_h =
        lane_count as f32 * lane_h + (lane_count.saturating_sub(1) as f32) * lane_gap;
    let protein_chart_top = top_header + exon_chart_h + 92.0;
    let protein_chart_h =
        lane_count as f32 * lane_h + (lane_count.saturating_sub(1) as f32) * lane_gap;
    let footer_top = protein_chart_top + protein_chart_h + 64.0;
    let dyn_h = (footer_top + 120.0).max(H + 140.0);

    let label_x = 44.0_f32;
    let left = 340.0_f32;
    let right = W - 44.0_f32;
    let width = (right - left).max(1.0);

    let region_start = view.region_start_1based.max(1);
    let region_end = view.region_end_1based.max(region_start);
    let region_span = (region_end - region_start + 1).max(1) as f32;
    let dominant_strand_is_reverse = {
        let mut plus = 0usize;
        let mut minus = 0usize;
        for lane in view.transcript_lanes.iter().filter(|lane| lane.mapped) {
            match lane.strand.trim() {
                "+" => plus += 1,
                "-" => minus += 1,
                _ => {}
            }
        }
        if plus == 0 && minus == 0 {
            for lane in &view.transcript_lanes {
                match lane.strand.trim() {
                    "+" => plus += 1,
                    "-" => minus += 1,
                    _ => {}
                }
            }
        }
        minus > plus
    };
    let left_axis_bp = if dominant_strand_is_reverse {
        region_end
    } else {
        region_start
    };
    let right_axis_bp = if dominant_strand_is_reverse {
        region_start
    } else {
        region_end
    };
    let x_for_genomic = |pos_1based: usize| -> f32 {
        let clamped = pos_1based.clamp(region_start, region_end);
        let rel = if dominant_strand_is_reverse {
            (region_end.saturating_sub(clamped)) as f32 / region_span
        } else {
            (clamped.saturating_sub(region_start)) as f32 / region_span
        };
        left + rel * width
    };

    let aa_max = view
        .protein_lanes
        .iter()
        .flat_map(|lane| lane.domains.iter().map(|domain| domain.end_aa))
        .chain(
            view.protein_lanes
                .iter()
                .filter_map(|lane| lane.expected_length_aa),
        )
        .chain(
            view.protein_lanes
                .iter()
                .filter_map(|lane| lane.reference_end_aa),
        )
        .max()
        .unwrap_or(1)
        .max(1);
    let aa_span = aa_max as f32;
    let x_for_aa = |aa_1based: usize| -> f32 {
        let clamped = aa_1based.clamp(1, aa_max);
        left + ((clamped.saturating_sub(1)) as f32 / aa_span) * width
    };
    let x_for_aa_f = |aa_1based: f32| -> f32 {
        let clamped = aa_1based.clamp(1.0, aa_max as f32);
        left + ((clamped - 1.0) / aa_span) * width
    };
    let top_geometry_kind = if view
        .transcript_geometry_mode
        .trim()
        .eq_ignore_ascii_case("cds")
    {
        "cds"
    } else {
        "exon"
    };
    let top_geometry_title = if top_geometry_kind == "cds" {
        "A) transcript exons (faint) + CDS (solid) architecture (coordinate-true)"
    } else {
        "A) transcript / exon architecture (coordinate-true)"
    };

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
                "Isoform architecture: {} panel '{}' on {}",
                view.gene_symbol, view.panel_id, view.seq_id
            ))
            .set("x", label_x)
            .set("y", 36)
            .set("font-family", "monospace")
            .set("font-size", 18)
            .set("fill", "#111111"),
        )
        .add(
            Text::new(format!(
                "genomic span {}..{} | isoforms={} | protein max={} aa",
                view.region_start_1based,
                view.region_end_1based,
                view.transcript_lanes.len(),
                aa_max
            ))
            .set("x", label_x)
            .set("y", 58)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#4b5563"),
        )
        .add(
            Text::new(format!(
                "display orientation: transcript 5'->3' left-to-right (dominant strand {}) | top-panel geometry: {}",
                if dominant_strand_is_reverse { "-" } else { "+" },
                top_geometry_kind
            ))
            .set("x", label_x)
            .set("y", 74)
            .set("font-family", "monospace")
            .set("font-size", 11)
            .set("fill", "#4b5563"),
        )
        .add(
            Text::new(top_geometry_title)
            .set("x", label_x)
            .set("y", 88)
            .set("font-family", "monospace")
            .set("font-size", 13)
            .set("fill", "#111827"),
        )
        .add(
            Line::new()
                .set("x1", left)
                .set("y1", top_header - 16.0)
                .set("x2", right)
                .set("y2", top_header - 16.0)
                .set("stroke", "#6b7280")
                .set("stroke-width", 1.0),
        )
        .add(
            Text::new(format!("{} bp", left_axis_bp))
                .set("x", left)
                .set("y", top_header - 20.0)
                .set("text-anchor", "start")
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#4b5563"),
        )
        .add(
            Text::new(format!("{} bp", right_axis_bp))
                .set("x", right)
                .set("y", top_header - 20.0)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#4b5563"),
        );

    for (idx, lane) in view.transcript_lanes.iter().enumerate() {
        let y = top_header + idx as f32 * (lane_h + lane_gap) + lane_h * 0.5;
        let label = lane
            .transcript_id
            .as_deref()
            .map(|tx| format!("{} ({tx})", lane.label))
            .unwrap_or_else(|| lane.label.clone());
        doc = doc.add(
            Text::new(label)
                .set("x", left - 12.0)
                .set("y", y + 3.5)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", if lane.mapped { "#111827" } else { "#b45309" }),
        );
        for intron in &lane.introns {
            let xa = x_for_genomic(intron.start_1based);
            let xb = x_for_genomic(intron.end_1based);
            let (x1, x2) = if xa <= xb { (xa, xb) } else { (xb, xa) };
            if (x2 - x1).abs() < 0.5 {
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
        if top_geometry_kind == "cds" {
            for exon in &lane.transcript_exons {
                let xa = x_for_genomic(exon.start_1based);
                let xb = x_for_genomic(exon.end_1based);
                let (x1, x2) = if xa <= xb { (xa, xb) } else { (xb, xa) };
                if (x2 - x1).abs() < 0.5 {
                    continue;
                }
                doc = doc.add(
                    Rectangle::new()
                        .set("x", x1)
                        .set("y", y - 7.5)
                        .set("width", (x2 - x1).max(1.0))
                        .set("height", 15.0)
                        .set("fill", "#93c5fd")
                        .set("fill-opacity", 0.25)
                        .set("stroke", "#60a5fa")
                        .set("stroke-width", 0.35),
                );
            }
        }
        for exon in &lane.exons {
            let xa = x_for_genomic(exon.start_1based);
            let xb = x_for_genomic(exon.end_1based);
            let (x1, x2) = if xa <= xb { (xa, xb) } else { (xb, xa) };
            if (x2 - x1).abs() < 0.5 {
                continue;
            }
            doc = doc.add(
                Rectangle::new()
                    .set("x", x1)
                    .set("y", y - 7.5)
                    .set("width", (x2 - x1).max(1.0))
                    .set("height", 15.0)
                    .set("fill", if lane.mapped { "#2563eb" } else { "#f59e0b" })
                    .set("fill-opacity", if lane.mapped { 0.85 } else { 0.45 })
                    .set("stroke", "#1f2937")
                    .set("stroke-width", 0.5),
            );
        }
        if let Some(tag) = lane.transactivation_class.as_deref() {
            doc = doc.add(
                Text::new(format!("TA={tag}"))
                    .set("x", right + 6.0)
                    .set("y", y + 3.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#374151"),
            );
        }
    }

    let protein_axis_y = protein_chart_top - 14.0;
    let genome_rail_y = top_header + exon_chart_h + 24.0;
    let mut boundary_set: BTreeSet<usize> = BTreeSet::new();
    for lane in &view.transcript_lanes {
        let boundary_exons = if top_geometry_kind == "cds" && !lane.transcript_exons.is_empty() {
            &lane.transcript_exons
        } else {
            &lane.exons
        };
        for exon in boundary_exons {
            if exon.end_1based >= exon.start_1based {
                boundary_set.insert(exon.start_1based);
                boundary_set.insert(exon.end_1based);
            }
        }
    }
    let boundaries = boundary_set.into_iter().collect::<Vec<_>>();
    if boundaries.len() >= 2 {
        let mut ribbon_bins: BTreeMap<(i32, i32, i32, i32), (f32, f32, f32, f32, usize)> =
            BTreeMap::new();
        let quantize_coord = |value: f32| -> i32 { (value * 1000.0).round() as i32 };
        let rail_x_a = x_for_genomic(*boundaries.first().unwrap_or(&region_start));
        let rail_x_b = x_for_genomic(*boundaries.last().unwrap_or(&region_end));
        let rail_left = rail_x_a.min(rail_x_b);
        let rail_right = rail_x_a.max(rail_x_b);
        doc = doc
            .add(
                Text::new("genome boundary rail")
                    .set("x", label_x)
                    .set("y", genome_rail_y + 4.0)
                    .set("font-family", "monospace")
                    .set("font-size", 10)
                    .set("fill", "#4b5563"),
            )
            .add(
                Line::new()
                    .set("x1", rail_left)
                    .set("y1", genome_rail_y)
                    .set("x2", rail_right)
                    .set("y2", genome_rail_y)
                    .set("stroke", "#6b7280")
                    .set("stroke-width", 1.0),
            );
        for boundary in &boundaries {
            let x = x_for_genomic(*boundary);
            doc = doc.add(
                Line::new()
                    .set("x1", x)
                    .set("y1", genome_rail_y - 3.0)
                    .set("x2", x)
                    .set("y2", genome_rail_y + 3.0)
                    .set("stroke", "#6b7280")
                    .set("stroke-width", 0.8)
                    .set("stroke-opacity", 0.8),
            );
        }
        for flank in boundaries.windows(2) {
            let flank_start = flank[0];
            let flank_end_exclusive = flank[1];
            if flank_end_exclusive <= flank_start {
                continue;
            }
            for lane in &view.transcript_lanes {
                if lane.cds_to_protein_segments.is_empty() {
                    continue;
                }
                for segment in &lane.cds_to_protein_segments {
                    if segment.aa_end < segment.aa_start {
                        continue;
                    }
                    let seg_start = segment.genomic_start_1based.min(segment.genomic_end_1based);
                    let seg_end = segment.genomic_start_1based.max(segment.genomic_end_1based);
                    // Avoid boundary-only tail overlaps: when flank starts exactly at
                    // a segment's inclusive end, this would otherwise create a tiny
                    // duplicate ribbon that visually looks like a second edge mapping.
                    if seg_end <= flank_start || seg_start >= flank_end_exclusive {
                        continue;
                    }
                    let overlap_start = flank_start.max(seg_start);
                    let overlap_end_exclusive = flank_end_exclusive.min(seg_end.saturating_add(1));
                    if overlap_end_exclusive <= overlap_start {
                        continue;
                    }

                    let seg_nt = seg_end.saturating_sub(seg_start).saturating_add(1).max(1) as f32;
                    let aa_span = segment
                        .aa_end
                        .saturating_sub(segment.aa_start)
                        .saturating_add(1)
                        .max(1) as f32;
                    let start_offset_nt = overlap_start.saturating_sub(seg_start) as f32;
                    let end_offset_nt = overlap_end_exclusive.saturating_sub(seg_start) as f32;
                    let aa_left_f = segment.aa_start as f32 + (start_offset_nt / seg_nt) * aa_span;
                    let aa_right_f = segment.aa_start as f32 + (end_offset_nt / seg_nt) * aa_span;
                    let aa_x_a = x_for_aa_f(aa_left_f);
                    let aa_x_b = x_for_aa_f(aa_right_f);
                    let (aa_left, aa_right) = if aa_x_a <= aa_x_b {
                        (aa_x_a, aa_x_b)
                    } else {
                        (aa_x_b, aa_x_a)
                    };
                    let g_x_a = x_for_genomic(overlap_start);
                    let g_x_b = x_for_genomic(overlap_end_exclusive);
                    let (g_left, g_right) = if g_x_a <= g_x_b {
                        (g_x_a, g_x_b)
                    } else {
                        (g_x_b, g_x_a)
                    };
                    if (g_right - g_left).abs() < 0.4 && (aa_right - aa_left).abs() < 0.4 {
                        continue;
                    }
                    let key = (
                        quantize_coord(g_left),
                        quantize_coord(g_right),
                        quantize_coord(aa_left),
                        quantize_coord(aa_right),
                    );
                    let entry = ribbon_bins
                        .entry(key)
                        .or_insert((g_left, g_right, aa_left, aa_right, 0usize));
                    entry.4 += 1;
                }
            }
        }
        for (_key, (g_left, g_right, aa_left, aa_right, support_count)) in ribbon_bins {
            let opacity = (0.14 + support_count.saturating_sub(1) as f32 * 0.04).min(0.42);
            let ribbon = Data::new()
                .move_to((g_left, genome_rail_y + 0.5))
                .line_to((g_right, genome_rail_y + 0.5))
                .line_to((aa_right, protein_axis_y + 0.5))
                .line_to((aa_left, protein_axis_y + 0.5))
                .close();
            let mut path = Path::new()
                .set("d", ribbon)
                .set("fill", "#64748b")
                .set("fill-opacity", opacity);
            if support_count > 1 {
                path = path
                    .set("stroke", "#475569")
                    .set("stroke-opacity", 0.18)
                    .set("stroke-width", 0.3);
            }
            doc = doc.add(path);
        }
    }

    doc = doc
        .add(
            Text::new("B) protein domain architecture")
                .set("x", label_x)
                .set("y", protein_chart_top - 32.0)
                .set("font-family", "monospace")
                .set("font-size", 13)
                .set("fill", "#111827"),
        )
        .add(
            Line::new()
                .set("x1", left)
                .set("y1", protein_chart_top - 14.0)
                .set("x2", right)
                .set("y2", protein_chart_top - 14.0)
                .set("stroke", "#6b7280")
                .set("stroke-width", 1.0),
        )
        .add(
            Text::new("1 aa")
                .set("x", left)
                .set("y", protein_chart_top - 18.0)
                .set("text-anchor", "start")
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#4b5563"),
        )
        .add(
            Text::new(format!("{aa_max} aa"))
                .set("x", right)
                .set("y", protein_chart_top - 18.0)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#4b5563"),
        );

    for (idx, lane) in view.protein_lanes.iter().enumerate() {
        let y = protein_chart_top + idx as f32 * (lane_h + lane_gap) + lane_h * 0.5;
        doc = doc.add(
            Text::new(lane.label.clone())
                .set("x", left - 12.0)
                .set("y", y + 3.5)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#111827"),
        );
        let max_domain_end = lane
            .domains
            .iter()
            .map(|domain| domain.end_aa.max(domain.start_aa))
            .max()
            .unwrap_or(1)
            .max(1);
        let lane_start = lane.reference_start_aa.unwrap_or(1).clamp(1, aa_max);
        let inferred_end = lane
            .expected_length_aa
            .unwrap_or(max_domain_end)
            .max(max_domain_end);
        let lane_end = lane
            .reference_end_aa
            .unwrap_or(inferred_end)
            .clamp(lane_start, aa_max);
        let rail_left = x_for_aa(lane_start);
        let rail_right = x_for_aa(lane_end);
        doc = doc.add(
            Line::new()
                .set("x1", rail_left)
                .set("y1", y)
                .set("x2", rail_right)
                .set("y2", y)
                .set("stroke", "#94a3b8")
                .set("stroke-width", 2.0),
        );
        for domain in &lane.domains {
            let domain_start = domain.start_aa.max(lane_start);
            let domain_end = domain.end_aa.max(domain.start_aa).min(lane_end);
            if domain_end < domain_start {
                continue;
            }
            let x1 = x_for_aa(domain_start);
            let x2 = x_for_aa(domain_end);
            let fill = domain.color_hex.as_deref().unwrap_or("#7c3aed");
            doc = doc
                .add(
                    Rectangle::new()
                        .set("x", x1)
                        .set("y", y - 7.0)
                        .set("width", (x2 - x1).max(1.0))
                        .set("height", 14.0)
                        .set("fill", fill)
                        .set("fill-opacity", 0.82)
                        .set("stroke", "#1f2937")
                        .set("stroke-width", 0.5),
                )
                .add(
                    Text::new(domain.name.clone())
                        .set("x", x1 + 2.0)
                        .set("y", y - 9.0)
                        .set("font-family", "monospace")
                        .set("font-size", 8)
                        .set("fill", "#334155"),
                );
        }
        if let Some(tag) = lane.transactivation_class.as_deref() {
            doc = doc.add(
                Text::new(format!("TA={tag}"))
                    .set("x", right + 6.0)
                    .set("y", y + 3.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#374151"),
            );
        }
    }

    let mut y = footer_top;
    if let Some(source) = view
        .panel_source
        .as_deref()
        .filter(|value| !value.trim().is_empty())
    {
        doc = doc.add(
            Text::new(format!("panel source: {source}"))
                .set("x", label_x)
                .set("y", y)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#334155"),
        );
        y += 14.0;
    }
    for warning in view.warnings.iter().take(6) {
        doc = doc.add(
            Text::new(format!("warning: {warning}"))
                .set("x", label_x)
                .set("y", y)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#b45309"),
        );
        y += 12.0;
    }
    for (line_idx, line) in wrap_text(&view.instruction, 140).into_iter().enumerate() {
        doc = doc.add(
            Text::new(line)
                .set("x", label_x)
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
        FeatureExpertView::IsoformArchitecture(isoform) => render_isoform_architecture(isoform),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature_expert::{
        IsoformArchitectureCdsAaSegment, IsoformArchitectureProteinDomain,
        IsoformArchitectureProteinLane, IsoformArchitectureTranscriptLane, SplicingRange,
    };

    fn isoform_test_view(transcript_strand: &str) -> IsoformArchitectureExpertView {
        IsoformArchitectureExpertView {
            seq_id: "tp53".to_string(),
            panel_id: "panel".to_string(),
            gene_symbol: "TP53".to_string(),
            transcript_geometry_mode: "exon".to_string(),
            panel_source: Some("test".to_string()),
            region_start_1based: 101,
            region_end_1based: 200,
            instruction: "test".to_string(),
            transcript_lanes: vec![IsoformArchitectureTranscriptLane {
                isoform_id: "i1".to_string(),
                label: "iso1".to_string(),
                transcript_id: Some("tx1".to_string()),
                transcript_feature_id: Some(1),
                strand: transcript_strand.to_string(),
                transcript_exons: vec![
                    SplicingRange {
                        start_1based: 105,
                        end_1based: 125,
                    },
                    SplicingRange {
                        start_1based: 145,
                        end_1based: 165,
                    },
                ],
                exons: vec![
                    SplicingRange {
                        start_1based: 110,
                        end_1based: 120,
                    },
                    SplicingRange {
                        start_1based: 150,
                        end_1based: 160,
                    },
                ],
                introns: vec![SplicingRange {
                    start_1based: 120,
                    end_1based: 150,
                }],
                mapped: true,
                transactivation_class: None,
                cds_to_protein_segments: vec![IsoformArchitectureCdsAaSegment {
                    genomic_start_1based: 110,
                    genomic_end_1based: 120,
                    aa_start: 10,
                    aa_end: 20,
                }],
                note: None,
            }],
            protein_lanes: vec![IsoformArchitectureProteinLane {
                isoform_id: "i1".to_string(),
                label: "iso1".to_string(),
                transcript_id: Some("tx1".to_string()),
                expected_length_aa: Some(100),
                reference_start_aa: None,
                reference_end_aa: None,
                domains: vec![IsoformArchitectureProteinDomain {
                    name: "dbd".to_string(),
                    start_aa: 10,
                    end_aa: 60,
                    color_hex: Some("#ff0000".to_string()),
                }],
                transactivation_class: None,
            }],
            warnings: vec![],
        }
    }

    fn extract_exon_x_positions(svg: &str) -> Vec<f32> {
        let mut xs = Vec::new();
        for part in svg.split("<rect").skip(1) {
            if !part.contains("fill=\"#2563eb\"") || !part.contains("height=\"15\"") {
                continue;
            }
            let Some(attr_start) = part.find("x=\"") else {
                continue;
            };
            let rest = &part[attr_start + 3..];
            let Some(attr_end) = rest.find('"') else {
                continue;
            };
            if let Ok(value) = rest[..attr_end].parse::<f32>() {
                xs.push(value);
            }
        }
        xs
    }

    fn connector_ribbon_count(svg: &str) -> usize {
        svg.matches("fill=\"#64748b\"").count()
    }

    #[test]
    fn isoform_renderer_keeps_forward_strand_exon_order_left_to_right() {
        let svg = render_isoform_architecture(&isoform_test_view("+"));
        let exon_x = extract_exon_x_positions(&svg);
        assert_eq!(exon_x.len(), 2);
        assert!(exon_x[0] < exon_x[1]);
        assert!(svg.contains("dominant strand +"));
        assert!(svg.contains("101 bp"));
        assert!(svg.contains("200 bp"));
    }

    #[test]
    fn isoform_renderer_flips_reverse_strand_exons_but_keeps_them_visible() {
        let svg = render_isoform_architecture(&isoform_test_view("-"));
        let exon_x = extract_exon_x_positions(&svg);
        assert_eq!(exon_x.len(), 2);
        assert!(exon_x[0] > exon_x[1]);
        assert!(svg.contains("dominant strand -"));
        assert!(svg.contains("200 bp"));
        assert!(svg.contains("101 bp"));
    }

    #[test]
    fn isoform_renderer_draws_cds_to_protein_connector_guides() {
        let view = isoform_test_view("+");
        let lane_count = view.transcript_lanes.len().max(1) as f32;
        let lane_h = 24.0_f32;
        let lane_gap = 10.0_f32;
        let top_header = 114.0_f32;
        let exon_chart_h = lane_count * lane_h + (lane_count - 1.0) * lane_gap;
        let protein_chart_top = top_header + exon_chart_h + 92.0;
        let protein_axis_y = protein_chart_top - 14.0 + 0.5;
        let svg = render_isoform_architecture(&view);
        assert!(svg.contains("genome boundary rail"));
        assert!(svg.contains("fill=\"#64748b\""));
        assert!(svg.contains("fill-opacity=\"0.14\""));
        assert!(svg.contains(&format!(",{protein_axis_y}")));
        assert!(!svg.contains(",233.5"));
    }

    #[test]
    fn isoform_renderer_merges_duplicate_connector_ribbons() {
        let single = isoform_test_view("+");
        let single_svg = render_isoform_architecture(&single);
        let single_count = connector_ribbon_count(&single_svg);

        let mut duplicated = isoform_test_view("+");
        let mut duplicated_lane = duplicated.transcript_lanes[0].clone();
        duplicated_lane.isoform_id = "i2".to_string();
        duplicated_lane.label = "iso2".to_string();
        duplicated_lane.transcript_id = Some("tx2".to_string());
        duplicated.transcript_lanes.push(duplicated_lane);

        let duplicated_svg = render_isoform_architecture(&duplicated);
        let duplicated_count = connector_ribbon_count(&duplicated_svg);
        assert_eq!(duplicated_count, single_count);
        assert!(duplicated_svg.contains("fill-opacity=\"0.18\""));
        assert!(duplicated_svg.contains("stroke=\"#475569\""));
    }
}
