//! Viewport-sized layout of query-referenced homology alignments.

use eframe::egui;
use gentle_protocol::GenomicRegionHomologyScreenReport;

/// Render a compact ordered-block matrix. Rows are transcript-derived promoter
/// windows and the horizontal axis is the selected query region.
pub(super) fn render_promoter_similarity_matrix(
    ui: &mut egui::Ui,
    report: &GenomicRegionHomologyScreenReport,
) -> usize {
    let Some(matrix) = report.promoter_similarity_matrix.as_ref() else {
        return 0;
    };
    let query_len = report.query.sequence.len().max(1);
    let row_height = 22.0;
    let label_width = 230.0;
    let track_width = 760.0;
    let mut rendered = 0;
    egui::ScrollArea::vertical()
        .id_salt(("promoter_similarity_matrix", &report.content_sha256))
        .max_height(360.0)
        .show_rows(ui, row_height, matrix.rows.len(), |ui, visible| {
            for index in visible {
                let row = &matrix.rows[index];
                rendered += 1;
                ui.horizontal(|ui| {
                    let gene = row
                        .gene_names
                        .first()
                        .or_else(|| row.gene_ids.first())
                        .map(String::as_str)
                        .unwrap_or("unlabelled gene");
                    let label = format!(
                        "{} | {} tx | {:.0}%",
                        gene,
                        row.transcript_ids.len(),
                        row.query_coverage_percent
                    );
                    ui.add_sized(
                        [label_width, row_height],
                        egui::Label::new(egui::RichText::new(&label).monospace()).truncate(),
                    )
                    .on_hover_text(format!(
                        "{}:{}-{} {}\nGenes: {}\nTranscripts: {}",
                        row.chromosome,
                        row.promoter_start_0based.saturating_add(1),
                        row.promoter_end_0based_exclusive,
                        row.strand.human_value(),
                        row.gene_names.join(", "),
                        row.transcript_ids.join(", ")
                    ));
                    let (response, painter) = ui.allocate_painter(
                        egui::vec2(track_width, row_height - 3.0),
                        egui::Sense::hover(),
                    );
                    painter.rect_filled(response.rect, 0.0, egui::Color32::from_gray(245));
                    for block in &row.blocks {
                        let x1 = response.rect.left()
                            + block.query_start_0based.min(query_len) as f32 / query_len as f32
                                * response.rect.width();
                        let x2 = response.rect.left()
                            + block.query_end_0based_exclusive.min(query_len) as f32
                                / query_len as f32
                                * response.rect.width();
                        let rect = egui::Rect::from_min_max(
                            egui::pos2(x1, response.rect.top()),
                            egui::pos2(x2.max(x1 + 3.0), response.rect.bottom()),
                        );
                        let alpha = (45.0
                            + 210.0 * block.identity_percent.clamp(0.0, 100.0) as f32 / 100.0)
                            .round() as u8;
                        painter.rect_filled(
                            rect,
                            0.0,
                            egui::Color32::from_rgba_unmultiplied(37, 99, 235, alpha),
                        );
                        painter.rect_stroke(
                            rect,
                            0.0,
                            egui::Stroke::new(
                                if block.order_break_before { 2.0 } else { 0.5 },
                                if block.order_break_before {
                                    egui::Color32::from_rgb(220, 38, 38)
                                } else {
                                    egui::Color32::from_rgb(29, 78, 216)
                                },
                            ),
                            egui::StrokeKind::Inside,
                        );
                        painter.text(
                            rect.left_top() + egui::vec2(3.0, 2.0),
                            egui::Align2::LEFT_TOP,
                            block.target_order.to_string(),
                            egui::FontId::monospace(10.0),
                            egui::Color32::WHITE,
                        );
                    }
                });
            }
        });
    rendered
}

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
            let mut output = ctx.run_ui(
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
            output.textures_delta.clear();
        }
        assert!(stats.0 > 0 && stats.0 < 50, "{stats:?}");
        assert!(
            stats
                .1
                .is_some_and(|start| (89_900..=90_000).contains(&start)),
            "{stats:?}"
        );
    }

    #[test]
    fn promoter_matrix_virtualizes_rows_and_accepts_order_breaks() {
        let report = GenomicRegionHomologyScreenReport {
            content_sha256: "synthetic-matrix".into(),
            query: gentle_protocol::GenomicRegionHomologyQueryBinding {
                sequence: "A".repeat(1_000),
                ..Default::default()
            },
            promoter_similarity_matrix: Some(gentle_protocol::PromoterSimilarityMatrix {
                rows: (0..1_000)
                    .map(|index| gentle_protocol::PromoterSimilarityMatrixRow {
                        row_id: format!("row-{index}"),
                        gene_names: vec![format!("GENE{index}")],
                        transcript_ids: vec![format!("TX{index}")],
                        blocks: vec![gentle_protocol::PromoterSimilarityBlock {
                            query_start_0based: 100,
                            query_end_0based_exclusive: 300,
                            target_order: 2,
                            identity_percent: 90.0,
                            order_break_before: true,
                            ..Default::default()
                        }],
                        ..Default::default()
                    })
                    .collect(),
                ..Default::default()
            }),
            ..Default::default()
        };
        let ctx = egui::Context::default();
        let mut rendered = 0;
        for _ in 0..2 {
            let mut output = ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(1_100.0, 480.0),
                    )),
                    ..Default::default()
                },
                |ui| rendered = render_promoter_similarity_matrix(ui, &report),
            );
            output.textures_delta.clear();
        }
        assert!(rendered > 0 && rendered < 40, "rendered {rendered}");
    }
}
