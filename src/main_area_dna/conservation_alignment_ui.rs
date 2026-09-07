//! Viewport-sized layout of query-referenced homology alignments.

use eframe::egui;
use gentle_protocol::GenomicRegionHomologyScreenReport;

/// Returns the number of instantiated rows and the first visible query tile.
pub(super) fn render_alignment(
    ui: &mut egui::Ui,
    report: &GenomicRegionHomologyScreenReport,
    jump_to_query_base: Option<usize>,
) -> (usize, Option<usize>) {
    const TILE_BP: usize = 100;
    let query_len = report.query.sequence.len();
    let rows_per_tile = report.alignment_rows.len().saturating_add(1);
    let row_count = query_len.div_ceil(TILE_BP).saturating_mul(rows_per_tile);
    let height = ui.text_style_height(&egui::TextStyle::Monospace).max(20.0);
    let mut scroll = egui::ScrollArea::both()
        .id_salt(("conservation_alignment", &report.content_sha256))
        .max_height(360.0)
        .auto_shrink([false, false]);
    if let Some(base) = jump_to_query_base {
        let tile = base.min(query_len.saturating_sub(1)) / TILE_BP;
        scroll = scroll.vertical_scroll_offset(
            (tile * rows_per_tile) as f32 * (height + ui.spacing().item_spacing.y),
        );
    }
    let mut rendered = 0;
    let mut first_tile = None;
    let output = scroll.show_rows(ui, height, row_count, |ui, visible| {
        for index in visible {
            rendered += 1;
            let start = (index / rows_per_tile) * TILE_BP;
            first_tile.get_or_insert(start);
            let end = (start + TILE_BP).min(query_len);
            let row = index % rows_per_tile;
            if row == 0 {
                ui.add_sized(
                    [1000.0, height],
                    egui::Label::new(
                        egui::RichText::new(format!("Query {}..{}", start + 1, end)).monospace(),
                    )
                    .truncate(),
                );
            } else if let Some(row) = report.alignment_rows.get(row - 1) {
                ui.horizontal(|ui| {
                    let label = if row.locus_class
                        == gentle_protocol::GenomicRegionHomologyLocusClass::Query
                    {
                        "QUERY".to_string()
                    } else {
                        format!("{}:{}", row.target_genome_id, row.subject_id)
                    };
                    ui.add_sized(
                        [230.0, height],
                        egui::Label::new(egui::RichText::new(&label).monospace()).truncate(),
                    )
                    .on_hover_text(label);
                    let sequence = row
                        .query_projection
                        .get(start..end)
                        .unwrap_or("[invalid projection]");
                    ui.add_sized(
                        [800.0, height],
                        egui::Label::new(egui::RichText::new(sequence).monospace()).extend(),
                    );
                });
            }
        }
    });
    if jump_to_query_base.is_some() {
        ui.scroll_to_rect(output.inner_rect, Some(egui::Align::Center));
    }
    (rendered, first_tile)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conservation_alignment_instantiates_only_visible_rows_and_jumps_to_block() {
        let report = GenomicRegionHomologyScreenReport {
            content_sha256: "synthetic".into(),
            query: gentle_protocol::GenomicRegionHomologyQueryBinding {
                sequence: "A".repeat(100_000),
                ..Default::default()
            },
            alignment_rows: (0..80)
                .map(|_| gentle_protocol::GenomicRegionHomologyAlignmentRow {
                    query_projection: ".".repeat(100_000),
                    ..Default::default()
                })
                .collect(),
            ..Default::default()
        };
        let ctx = egui::Context::default();
        let mut stats = (0, None);
        for _ in 0..2 {
            let _ = ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(1100.0, 480.0),
                    )),
                    ..Default::default()
                },
                |ui| {
                    stats = render_alignment(ui, &report, Some(90_000));
                },
            );
        }
        assert!(stats.0 > 0 && stats.0 < 50, "{stats:?}");
        assert!(
            stats
                .1
                .is_some_and(|start| (89_900..=90_000).contains(&start)),
            "{stats:?}"
        );
    }
}
