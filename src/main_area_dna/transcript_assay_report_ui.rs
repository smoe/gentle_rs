//! Bounded navigation and engine-owned coverage explanations for saved assay panels.

use crate::engine::{GentleEngine, TranscriptAssayPanelReport};
use eframe::egui;
use std::ops::Range;

fn page_range(total: usize, page_size: usize, page: usize) -> Range<usize> {
    let last_page = total.saturating_sub(1) / page_size;
    let start = page.min(last_page) * page_size;
    start..start.saturating_add(page_size).min(total)
}

/// Page only the presentation, never the report or its biological coverage set.
pub(super) fn page(
    ui: &mut egui::Ui,
    report_id: &str,
    section: &str,
    total: usize,
    page_size: usize,
) -> Range<usize> {
    let id = ui.make_persistent_id(("transcript_assay_report_page", report_id, section));
    let pages = total.div_ceil(page_size).max(1);
    let mut current = ui.data(|data| data.get_temp::<usize>(id).unwrap_or(0));
    current = current.min(pages - 1);
    ui.push_id(id, |ui| {
        ui.horizontal_wrapped(|ui| {
            ui.label(section);
            if pages > 1 {
                if ui
                    .add_enabled(current > 0, egui::Button::new("Previous"))
                    .clicked()
                {
                    current -= 1;
                }
                let mut displayed = current + 1;
                ui.add(
                    egui::DragValue::new(&mut displayed)
                        .range(1..=pages)
                        .prefix("Page "),
                );
                current = displayed - 1;
                ui.label(format!("of {pages}"));
                if ui
                    .add_enabled(current + 1 < pages, egui::Button::new("Next"))
                    .clicked()
                {
                    current += 1;
                }
            }
            let range = page_range(total, page_size, current);
            let first = if total == 0 { 0 } else { range.start + 1 };
            ui.label(format!("{first}-{} of {total}", range.end));
        });
    });
    ui.data_mut(|data| data.insert_temp(id, current));
    page_range(total, page_size, current)
}

/// Use the same record/digest coverage accounting as the experimental handoff.
pub(super) fn coverage(ui: &mut egui::Ui, report: &TranscriptAssayPanelReport) {
    let summary = GentleEngine::experimental_assay_coverage_summary(report);
    ui.group(|ui| {
        ui.strong("Coverage scope: what this panel establishes");
        ui.label(format!(
            "Requested universe: {} | annotation release: {}",
            summary.coverage_universe_kind,
            summary.annotation_release.as_deref().unwrap_or("not recorded")
        ));
        for line in &summary.summary_lines {
            ui.label(line);
        }
        for (label, ids) in [
            ("Unresolved requested targets", &report.coverage_resolution.unresolved_target_ids),
            ("Ambiguous requested targets", &report.coverage_resolution.ambiguous_target_ids),
        ] {
            if !ids.is_empty() {
                ui.label(format!("{label}: {}", ids.join(", ")));
            }
        }
        if !report.unresolved_group_pairs.is_empty() || !report.unresolved_coverage_target_pairs.is_empty() {
            ui.label(format!(
                "Recorded unseparated cDNA-class pairs: {}; unseparated coverage-target pairs: {}.",
                report.unresolved_group_pairs.len(), report.unresolved_coverage_target_pairs.len()
            ));
        }
        ui.small("Coverage is not isoform discrimination. Identical mature cDNAs cannot be separated by sequence-based primers. A completed design objective is not whole-reference specificity, experimental validation or order approval.");
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(
        ctx: &egui::Context,
        events: Vec<egui::Event>,
        render: impl FnOnce(&mut egui::Ui),
    ) -> egui::FullOutput {
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1600.0, 12000.0),
            )),
            events,
            ..Default::default()
        });
        crate::egui_compat::show_central_panel_for_test_context(
            ctx,
            egui::CentralPanel::default(),
            render,
        );
        crate::egui_compat::end_test_pass(ctx)
    }

    fn labels(output: &egui::FullOutput) -> Vec<(String, egui::Rect)> {
        fn collect(shape: &egui::epaint::Shape, labels: &mut Vec<(String, egui::Rect)>) {
            match shape {
                egui::epaint::Shape::Text(text) => labels.push((
                    text.galley.job.text.clone(),
                    egui::Rect::from_min_size(text.pos, text.galley.size()),
                )),
                egui::epaint::Shape::Vec(shapes) => {
                    for shape in shapes {
                        collect(shape, labels);
                    }
                }
                _ => {}
            }
        }
        let mut result = Vec::new();
        for shape in &output.shapes {
            collect(&shape.shape, &mut result);
        }
        result
    }

    fn rendered_text(ctx: &egui::Context, render: impl FnOnce(&mut egui::Ui)) -> String {
        labels(&frame(ctx, vec![], render))
            .into_iter()
            .map(|(text, _)| text)
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn pages_cover_all_transcripts_assays_and_bands_without_duplicates() {
        for (total, size) in [
            (95_usize, 40_usize),
            (25, 24),
            (161, 80),
            (80, 80),
            (1, 40),
            (0, 40),
        ] {
            let rows = (0..total.div_ceil(size))
                .flat_map(|page| page_range(total, size, page))
                .collect::<Vec<_>>();
            assert_eq!(rows, (0..total).collect::<Vec<_>>());
            assert!(page_range(total, size, usize::MAX).end <= total);
        }
    }

    #[test]
    fn page_state_is_report_and_section_scoped_and_clamps_after_shrink() {
        let ctx = egui::Context::default();
        let text = rendered_text(&ctx, |ui| {
            let id =
                ui.make_persistent_id(("transcript_assay_report_page", "panel_a", "Transcripts"));
            ui.data_mut(|data| data.insert_temp(id, 2_usize));
            assert_eq!(page(ui, "panel_a", "Transcripts", 95, 40), 80..95);
            assert_eq!(page(ui, "panel_b", "Transcripts", 95, 40), 0..40);
            assert_eq!(page(ui, "panel_a", "Assays", 25, 24), 0..24);
        });
        assert!(text.contains("81-95 of 95"), "{text}");
        rendered_text(&ctx, |ui| {
            assert_eq!(page(ui, "panel_a", "Transcripts", 3, 40), 0..3);
            assert_eq!(page(ui, "panel_b", "Transcripts", 0, 40), 0..0);
        });
    }

    #[test]
    fn coverage_renders_engine_accounting_without_claiming_all_annotation_was_assessed() {
        let mut report: TranscriptAssayPanelReport = serde_json::from_str(include_str!("../../docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_sybr_juc_panel.report.json")).unwrap();
        report.coverage_resolution.annotated_transcript_count = report.transcript_rows.len() + 1;
        report.coverage_resolution.annotated_equivalence_group_count =
            Some(report.equivalence_groups.len() + 1);
        report.coverage_resolution.excluded_annotated_transcript_ids = vec!["outside_scope".into()];
        report.coverage_universe.kind =
            crate::engine::TranscriptAssayCoverageUniverseKind::ExplicitTranscripts;
        report.coverage_resolution.unresolved_target_ids = vec!["missing_annotation".into()];
        let before = serde_json::to_value(&report).unwrap();
        let summary = GentleEngine::experimental_assay_coverage_summary(&report);
        let text = rendered_text(&egui::Context::default(), |ui| coverage(ui, &report));
        for line in summary.summary_lines {
            assert!(text.contains(&line), "{text}");
        }
        assert!(text.contains("outside_scope"));
        assert!(text.contains("unassessed, not uncovered"));
        assert!(text.contains("Unresolved requested targets: missing_annotation"));
        assert!(text.contains("Coverage is not isoform discrimination"));
        assert_eq!(serde_json::to_value(&report).unwrap(), before);
    }

    #[test]
    fn panel_matrix_navigation_reaches_rows_assays_and_bands_beyond_old_limits() {
        use crate::engine::{
            TranscriptAssayBandSizeRow, TranscriptAssayPanelAssay,
            TranscriptAssayPanelTranscriptRow,
        };
        // Synthetic presentation-only records, not designed or validated assays.
        let report = TranscriptAssayPanelReport {
            report_id: "large_ui_panel".into(),
            assay_kind: crate::engine::TranscriptAssayKind::EndpointRtPcr,
            transcript_rows: (0..95)
                .map(|i| TranscriptAssayPanelTranscriptRow {
                    transcript_id: format!("matrix_transcript_{i}"),
                    transcript_feature_id: i,
                    ..Default::default()
                })
                .collect(),
            selected_assays: (0..25)
                .map(|i| TranscriptAssayPanelAssay {
                    assay_id: format!("assay_{i}"),
                    rank: i + 1,
                    ..Default::default()
                })
                .collect(),
            band_size_matrix: (0..81)
                .map(|i| TranscriptAssayBandSizeRow {
                    transcript_id: format!("band_transcript_{i}"),
                    ..Default::default()
                })
                .collect(),
            ..Default::default()
        };
        let before = serde_json::to_value(&report).unwrap();
        let ctx = egui::Context::default();
        let render = |ui: &mut egui::Ui| {
            super::super::MainAreaDna::render_transcript_assay_panel_report(ui, &report)
        };
        frame(&ctx, vec![], render);
        let mut output = frame(&ctx, vec![], render);
        for section in ["Transcripts", "Transcripts", "Assays", "Band rows"] {
            let labels = labels(&output);
            let section_rect = labels.iter().find(|(text, _)| text == section).unwrap().1;
            let next = labels
                .iter()
                .find(|(text, rect)| {
                    text == "Next" && (rect.center().y - section_rect.center().y).abs() < 5.0
                })
                .unwrap()
                .1
                .center();
            let events = vec![
                egui::Event::PointerMoved(next),
                egui::Event::PointerButton {
                    pos: next,
                    button: egui::PointerButton::Primary,
                    pressed: true,
                    modifiers: egui::Modifiers::NONE,
                },
                egui::Event::PointerButton {
                    pos: next,
                    button: egui::PointerButton::Primary,
                    pressed: false,
                    modifiers: egui::Modifiers::NONE,
                },
            ];
            frame(&ctx, events, render);
            output = frame(&ctx, vec![], render);
        }
        let labels = labels(&output)
            .into_iter()
            .map(|(text, _)| text)
            .collect::<Vec<_>>();
        for expected in [
            "81-95 of 95",
            "25-25 of 25",
            "81-81 of 81",
            "matrix_transcript_94",
            "A25",
            "band_transcript_80",
            "Not assessed",
        ] {
            assert!(
                labels.iter().any(|text| text == expected),
                "Missing {expected}: {labels:?}"
            );
        }
        assert_eq!(serde_json::to_value(&report).unwrap(), before);
    }
}
