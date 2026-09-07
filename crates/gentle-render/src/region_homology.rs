//! Deterministic SVG rendering for query-referenced genomic-region homology.

use gentle_protocol::{
    GenomicRegionHomologyLocusClass, GenomicRegionHomologyScreenReport,
    GenomicRegionHomologySupportClass,
};
use svg::Document;
use svg::node::element::{Group, Line, Rectangle, Text};

const WIDTH: f32 = 1400.0;
const LEFT: f32 = 230.0;
const RIGHT: f32 = 36.0;
const TILE_BP: usize = 100;
const CHAR_WIDTH: f32 = 9.2;
const ROW_PITCH: f32 = 17.0;
const MAX_RENDERED_ROWS: usize = 100;

fn support_color(class: GenomicRegionHomologySupportClass) -> &'static str {
    match class {
        GenomicRegionHomologySupportClass::ExpectedOrtholog => "#15803d",
        GenomicRegionHomologySupportClass::CrossSpeciesUnassigned => "#0e7490",
        GenomicRegionHomologySupportClass::SameGenomeNonself => "#b45309",
    }
}

fn row_color(class: GenomicRegionHomologyLocusClass) -> &'static str {
    match class {
        GenomicRegionHomologyLocusClass::Query => "#111827",
        GenomicRegionHomologyLocusClass::ExpectedOrtholog => "#166534",
        GenomicRegionHomologyLocusClass::CrossSpeciesUnassigned => "#155e75",
        GenomicRegionHomologyLocusClass::SameGenomeSelf => "#4b5563",
        GenomicRegionHomologyLocusClass::SameGenomeNonself => "#92400e",
    }
}

fn substitution_color(base: char) -> &'static str {
    match base.to_ascii_uppercase() {
        'A' => "#15803d",
        'C' => "#2563eb",
        'G' => "#d97706",
        'T' | 'U' => "#dc2626",
        '-' => "#6b7280",
        _ => "#7c3aed",
    }
}

fn short_label(value: &str, max_chars: usize) -> String {
    if value.chars().count() <= max_chars {
        return value.to_string();
    }
    let mut out = value
        .chars()
        .take(max_chars.saturating_sub(1))
        .collect::<String>();
    out.push('…');
    out
}

/// Render a query-first, insertion-free alignment plus whole-region support.
pub fn render_genomic_region_homology_svg(report: &GenomicRegionHomologyScreenReport) -> String {
    let query_len = report.query.sequence.len().max(1);
    let rendered_rows = report
        .alignment_rows
        .iter()
        .filter(|row| row.locus_class != GenomicRegionHomologyLocusClass::Query)
        .take(MAX_RENDERED_ROWS)
        .collect::<Vec<_>>();
    let tile_count = query_len.div_ceil(TILE_BP).max(1);
    let rows_per_tile = rendered_rows.len() + 1;
    let overview_top = 92.0_f32;
    let overview_height = 58.0_f32;
    let alignment_top = overview_top + overview_height + 52.0;
    let tile_height = 42.0 + rows_per_tile as f32 * ROW_PITCH + 28.0;
    let footer_top = alignment_top + tile_count as f32 * tile_height + 18.0;
    let footer_lines = 4 + report.non_claims.len().min(6);
    let height = footer_top + footer_lines as f32 * 16.0 + 36.0;
    let track_width = WIDTH - LEFT - RIGHT;
    let x_for = |position: usize| LEFT + position as f32 / query_len as f32 * track_width;

    let mut document = Document::new()
        .set("viewBox", (0, 0, WIDTH, height))
        .set("width", WIDTH)
        .set("height", height)
        .set("data-gentle-schema", report.schema.as_str())
        .set("data-report-sha256", report.content_sha256.as_str())
        .set("data-query-length", query_len)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", WIDTH)
                .set("height", height)
                .set("fill", "#ffffff"),
        )
        .add(
            Text::new(format!(
                "Genomic-region homology: {} / {}",
                report.query.set_id, report.query.region.region_id
            ))
            .set("x", 36)
            .set("y", 38)
            .set("font-family", "sans-serif")
            .set("font-size", 22)
            .set("font-weight", 600)
            .set("fill", "#111827"),
        )
        .add(
            Text::new(format!(
                "{} bp | {} targets | {} retained loci | query-referenced; target insertions omitted",
                query_len,
                report.targets.len(),
                report.loci.len()
            ))
            .set("x", 36)
            .set("y", 64)
            .set("font-family", "sans-serif")
            .set("font-size", 12)
            .set("fill", "#4b5563"),
        );

    document = document
        .add(
            Text::new("Whole-region exact support")
                .set("x", 36)
                .set("y", overview_top + 13.0)
                .set("font-family", "sans-serif")
                .set("font-size", 12)
                .set("font-weight", 600)
                .set("fill", "#374151"),
        )
        .add(
            Rectangle::new()
                .set("x", LEFT)
                .set("y", overview_top)
                .set("width", track_width)
                .set("height", overview_height)
                .set("fill", "#f8fafc")
                .set("stroke", "#cbd5e1"),
        );
    for (lane_index, class) in [
        GenomicRegionHomologySupportClass::ExpectedOrtholog,
        GenomicRegionHomologySupportClass::CrossSpeciesUnassigned,
        GenomicRegionHomologySupportClass::SameGenomeNonself,
    ]
    .into_iter()
    .enumerate()
    {
        let y = overview_top + 8.0 + lane_index as f32 * 16.0;
        document = document.add(
            Text::new(class.as_str())
                .set("x", 36)
                .set("y", y + 9.0)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", support_color(class)),
        );
        for block in report
            .conserved_blocks
            .iter()
            .filter(|block| block.support_class == class)
        {
            let x1 = x_for(block.query_start_0based);
            let x2 = x_for(block.query_end_0based_exclusive);
            document = document.add(
                Rectangle::new()
                    .set("x", x1)
                    .set("y", y)
                    .set("width", (x2 - x1).max(1.0))
                    .set("height", 10)
                    .set("fill", support_color(class))
                    .set("fill-opacity", 0.78)
                    .set("data-block-id", block.block_id.as_str())
                    .set("data-support-class", class.as_str())
                    .set("data-query-start", block.query_start_0based)
                    .set("data-query-end", block.query_end_0based_exclusive),
            );
        }
    }

    for tile_index in 0..tile_count {
        let start = tile_index * TILE_BP;
        let end = ((tile_index + 1) * TILE_BP).min(query_len);
        let tile_top = alignment_top + tile_index as f32 * tile_height;
        document = document
            .add(
                Text::new(format!("query {}..{}", start + 1, end))
                    .set("x", 36)
                    .set("y", tile_top)
                    .set("font-family", "monospace")
                    .set("font-size", 11)
                    .set("fill", "#374151"),
            )
            .add(
                Line::new()
                    .set("x1", LEFT)
                    .set("x2", LEFT + (end - start) as f32 * CHAR_WIDTH)
                    .set("y1", tile_top + 7.0)
                    .set("y2", tile_top + 7.0)
                    .set("stroke", "#d1d5db"),
            );
        let query_y = tile_top + 27.0;
        document = document
            .add(
                Text::new("QUERY")
                    .set("x", 36)
                    .set("y", query_y)
                    .set("font-family", "monospace")
                    .set("font-size", 11)
                    .set("font-weight", 700)
                    .set("fill", "#111827"),
            )
            .add(
                Text::new(report.query.sequence[start..end].to_string())
                    .set("x", LEFT)
                    .set("y", query_y)
                    .set("xml:space", "preserve")
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("fill", "#111827")
                    .set("data-row-class", "query"),
            );
        for (row_index, row) in rendered_rows.iter().enumerate() {
            let y = query_y + (row_index + 1) as f32 * ROW_PITCH;
            let label = format!(
                "{} | {}",
                short_label(&row.target_genome_id, 16),
                short_label(&row.subject_id, 22)
            );
            let mut group = Group::new()
                .set("data-row-id", row.row_id.as_str())
                .set("data-row-class", row.locus_class.as_str())
                .set("data-omitted-insertions", row.omitted_insertion_ids.len())
                .add(
                    Text::new(label)
                        .set("x", 36)
                        .set("y", y)
                        .set("font-family", "monospace")
                        .set("font-size", 10)
                        .set("fill", row_color(row.locus_class)),
                );
            let projection = row.query_projection.as_bytes();
            for position in start..end {
                let symbol = projection.get(position).copied().unwrap_or(b' ') as char;
                if symbol == ' ' {
                    continue;
                }
                let fill = if symbol == '.' {
                    "#9ca3af"
                } else {
                    substitution_color(symbol)
                };
                group = group.add(
                    Text::new(symbol.to_string())
                        .set("x", LEFT + (position - start) as f32 * CHAR_WIDTH)
                        .set("y", y)
                        .set("font-family", "monospace")
                        .set("font-size", 12)
                        .set("fill", fill),
                );
            }
            document = document.add(group);
        }
    }

    let target_row_count = report
        .alignment_rows
        .iter()
        .filter(|row| row.locus_class != GenomicRegionHomologyLocusClass::Query)
        .count();
    let hidden = target_row_count.saturating_sub(rendered_rows.len());
    let summary = format!(
        "Rows {}{} | omitted insertions {} | projection conflicts {} | same-genome non-self loci {} ({:.1}% query coverage)",
        rendered_rows.len(),
        if hidden > 0 {
            format!(" (+{hidden} retained in JSON)")
        } else {
            String::new()
        },
        report.omitted_insertions.len(),
        report
            .alignment_rows
            .iter()
            .map(|row| row.conflicts.len())
            .sum::<usize>(),
        report.same_genome_nonself_locus_count,
        report.same_genome_nonself_query_coverage_percent,
    );
    document = document.add(
        Text::new(summary)
            .set("x", 36)
            .set("y", footer_top)
            .set("font-family", "sans-serif")
            .set("font-size", 11)
            .set("fill", "#374151"),
    );
    document = document.add(
        Text::new("Legend: . exact match | coloured base substitution | - target deletion | blank no accepted HSP")
            .set("x", 36)
            .set("y", footer_top + 18.0)
            .set("font-family", "sans-serif")
            .set("font-size", 11)
            .set("fill", "#4b5563"),
    );
    for (index, statement) in report.non_claims.iter().take(6).enumerate() {
        document = document.add(
            Text::new(format!("Non-claim: {statement}"))
                .set("x", 36)
                .set("y", footer_top + 42.0 + index as f32 * 16.0)
                .set("font-family", "sans-serif")
                .set("font-size", 10)
                .set("fill", "#6b7280"),
        );
    }
    document.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle_protocol::{GenomicRegionHomologyAlignmentRow, GenomicRegionHomologyQueryBinding};

    #[test]
    fn homology_svg_renders_query_once_and_keeps_structural_row_markers() {
        let report = GenomicRegionHomologyScreenReport {
            schema: gentle_protocol::GENOMIC_REGION_HOMOLOGY_SCREEN_SCHEMA.to_string(),
            content_sha256: "sha256:synthetic".to_string(),
            query: GenomicRegionHomologyQueryBinding {
                set_id: "fixture".to_string(),
                sequence: "AACCGGTT".to_string(),
                ..Default::default()
            },
            alignment_rows: vec![
                GenomicRegionHomologyAlignmentRow {
                    row_id: "query".to_string(),
                    locus_class: GenomicRegionHomologyLocusClass::Query,
                    query_projection: "AACCGGTT".to_string(),
                    ..Default::default()
                },
                GenomicRegionHomologyAlignmentRow {
                    row_id: "ortholog".to_string(),
                    target_genome_id: "target".to_string(),
                    subject_id: "chr1".to_string(),
                    locus_class: GenomicRegionHomologyLocusClass::ExpectedOrtholog,
                    query_projection: "....A...".to_string(),
                    omitted_insertion_ids: vec!["insertion_1".to_string()],
                    ..Default::default()
                },
            ],
            ..Default::default()
        };

        let svg = render_genomic_region_homology_svg(&report);
        assert_eq!(svg.matches("data-row-class=\"query\"").count(), 1);
        assert!(svg.contains("data-row-id=\"ortholog\""));
        assert!(svg.contains("data-row-class=\"expected_ortholog\""));
        assert!(svg.contains("data-omitted-insertions=\"1\""));
    }
}
