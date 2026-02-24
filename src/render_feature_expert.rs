use crate::feature_expert::{FeatureExpertView, RestrictionSiteExpertView, TfbsExpertView};
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

pub fn render_feature_expert_svg(view: &FeatureExpertView) -> String {
    match view {
        FeatureExpertView::Tfbs(tfbs) => render_tfbs(tfbs),
        FeatureExpertView::RestrictionSite(re) => render_restriction(re),
    }
}
