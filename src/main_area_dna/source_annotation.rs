//! Read-only source comparison shared by the DNA map and Splicing Structure tab.

use super::locus_inspector::{AnnotationFilter, membership_color, membership_label};
use super::*;
use gentle_engine::transcript_presentation::projection::LocalAnnotationRow;
use gentle_protocol::transcript_presentation::TranscriptInterval;

impl MainAreaDna {
    pub(super) fn render_annotation_filters(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            for (mode, label) in [
                (AnnotationFilter::All, "All sources"),
                (AnnotationFilter::Ensembl, "Ensembl chains"),
                (AnnotationFilter::RefSeq, "RefSeq chains"),
                (AnnotationFilter::Shared, "Shared exon chains"),
                (AnnotationFilter::SourceOnly, "Source-only chains"),
            ] {
                ui.selectable_value(&mut self.splicing_locus_source_filter, mode, label);
            }
        });
    }

    pub(super) fn selected_annotation_chain(
        &self,
        feature_id: Option<usize>,
    ) -> Option<(Vec<TranscriptInterval>, i8)> {
        let dna = self.dna.try_read().ok()?;
        let feature = dna.features().get(feature_id?)?;
        if !GentleEngine::is_splicing_transcript_feature(feature) {
            return None;
        }
        let mut ranges = Vec::new();
        collect_location_ranges_usize(&feature.location, &mut ranges);
        Some((
            ranges
                .into_iter()
                .map(|(start, end)| TranscriptInterval {
                    start_1based: start as u64 + 1,
                    end_1based: end as u64,
                })
                .collect(),
            if feature_is_reverse(feature) { -1 } else { 1 },
        ))
    }

    /// Navigation affects the viewport only, never PCR selection or feature identities.
    pub(super) fn focus_source_annotation_row(
        &mut self,
        report: &GeneLocusEvidenceDisplayReport,
        row: &LocalAnnotationRow,
    ) -> Result<(), String> {
        self.verify_splicing_locus_binding(report)?;
        let start = row
            .exons
            .iter()
            .map(|e| e.start_1based)
            .min()
            .ok_or("This chain is outside the loaded locus")?;
        let end = row.exons.iter().map(|e| e.end_1based).max().unwrap();
        self.set_linear_viewport(start as usize - 1, (end - start + 1) as usize);
        Ok(())
    }

    pub(super) fn render_source_annotation_comparison(
        &mut self,
        ui: &mut egui::Ui,
        surface: &str,
        feature_id: Option<usize>,
    ) {
        let has_annotation = self
            .splicing_locus_report
            .as_ref()
            .is_some_and(|r| r.transcript_presentation.is_some());
        ui.push_id(("source_annotation", surface, self.panel_scope_key()), |ui| {
            egui::CollapsingHeader::new("Source annotation comparison (read-only)")
                .default_open(has_annotation)
                .show(ui, |ui| {
                    if ui.small_button("Load source comparison report...").on_hover_text("Load a hash-bound gene-locus report. Build one in Splicing Expert > Locus figure using annotation sources.").clicked() {
                        self.load_splicing_locus_report_json_dialog();
                    }
                    if !self.splicing_locus_status.is_empty() { ui.small(&self.splicing_locus_status); }
                    let Some(report) = self.splicing_locus_report.clone().filter(|r| r.transcript_presentation.is_some()) else {
                        ui.small("Supply Ensembl/RefSeq annotations in Splicing Expert > Locus figure, or load its report here. A GenBank file alone does not establish cross-source agreement.");
                        return;
                    };
                    if self.is_circular() {
                        ui.small("Source comparison requires an anchored linear sequence.");
                        return;
                    }
                    if let Err(error) = self.cached_locus_binding(&report) {
                        ui.colored_label(egui::Color32::DARK_RED, format!("Comparison withheld: {error}"));
                        return;
                    }
                    let prepared = self.splicing_locus_presentation.clone();
                    if let Some(error) = &prepared.annotation_error {
                        ui.colored_label(egui::Color32::DARK_RED, error);
                        return;
                    }
                    let Some(comparison) = &prepared.annotation else { return; };
                    self.render_annotation_filters(ui);
                    ui.small("Blue: Ensembl. Orange: RefSeq. Green: exact full exon chain in both. Thin: exon; thick: CDS; tick: annotated start. CDS/phase alternatives stay separate.");
                    ui.small("Comparison records are NOT added to primer-design targets or loaded features. Source-only is not biological absence; an unsupplied provider remains unassessed.");
                    let selected_chain = self.selected_annotation_chain(feature_id);
                    let (start, span, len) = self.current_linear_viewport();
                    let end = start.saturating_add(span).min(len);
                    ui.horizontal_wrapped(|ui| {
                        ui.small(format!("DNA viewport: {}..{} bp. Click to highlight; double-click to fit a chain.", start + 1, end));
                        if ui.small_button("Fit comparison locus").clicked() {
                            if let Err(error) = self.verify_splicing_locus_binding(&report) {
                                self.splicing_locus_status = error;
                            } else {
                                let start = report.locus_local_start_1based.max(1);
                                self.set_linear_viewport(start - 1, report.locus_local_end_1based.saturating_sub(start) + 1);
                            }
                        }
                        if ui.small_button("Clear comparison highlight").clicked() {
                            self.source_annotation_selected_structure = None;
                        }
                    });
                    let rows: Vec<_> = comparison.rows.iter().filter(|r| self.splicing_locus_source_filter.accepts(Some(r.comparison.membership))).collect();
                    let visible = rows.iter().filter(|r| r.exons.iter().any(|e| e.end_1based > start as u64 && e.start_1based <= end as u64)).count();
                    ui.small(format!("{} structure rows match the filter; {visible} have exons in this viewport. Bold border also marks a full-chain match to the selected project transcript.", rows.len()));
                    let max_height = (ui.available_height() * 0.35).clamp(52.0, 180.0);
                    egui::ScrollArea::vertical().id_salt("rows").max_height(max_height).auto_shrink([false, true])
                        .show_rows(ui, 46.0, rows.len(), |ui, range| {
                            for row in &rows[range] {
                                let selected = self.source_annotation_selected_structure.as_ref() == Some(&row.structure_id);
                                let matches_project = selected_chain.as_ref().is_some_and(|(exons, strand)| row.matches_local_chain(exons, *strand));
                                let color = membership_color(row.comparison.membership);
                                let ids = row.records.iter().map(|r| r.structure.transcript_id.as_str()).collect::<Vec<_>>().join(" / ");
                                let (rect, response) = ui.allocate_exact_size(egui::vec2(ui.available_width(), 46.0), egui::Sense::click());
                                let painter = ui.painter_at(rect);
                                if selected || matches_project {
                                    painter.rect_stroke(rect.shrink(1.0), 2.0, egui::Stroke::new(2.0, color), egui::StrokeKind::Inside);
                                }
                                painter.text(rect.left_top() + egui::vec2(4.0, 2.0), egui::Align2::LEFT_TOP, format!("{} | {ids}", membership_label(row.comparison.membership)), egui::FontId::proportional(11.0), color);
                                let axis = rect.shrink2(egui::vec2(5.0, 0.0));
                                let x = |position: u64| axis.left() + axis.width() * ((position as f64 - (start + 1) as f64) / span.max(1) as f64).clamp(0.0, 1.0) as f32;
                                let y = rect.bottom() - 12.0;
                                if let (Some(lo), Some(hi)) = (row.exons.iter().map(|e| e.start_1based).min(), row.exons.iter().map(|e| e.end_1based).max()) {
                                    if hi > start as u64 && lo <= end as u64 {
                                        painter.line_segment([egui::pos2(x(lo), y), egui::pos2(x(hi + 1), y)], egui::Stroke::new(1.0, color));
                                    }
                                }
                                for (i, h) in row.exons.iter().map(|i| (i, 5.0)).chain(row.cds.iter().flatten().map(|c| (&c.interval, 11.0))) {
                                    if i.end_1based > start as u64 && i.start_1based <= end as u64 {
                                        painter.rect_filled(egui::Rect::from_min_max(egui::pos2(x(i.start_1based), y - h / 2.0), egui::pos2(x(i.end_1based + 1).max(x(i.start_1based) + 1.0), y + h / 2.0)), 1.0, color);
                                    }
                                }
                                if let Some(pos) = row.annotated_start_1based.filter(|p| *p > start as u64 && *p <= end as u64) {
                                    painter.line_segment([egui::pos2(x(pos), y - 8.0), egui::pos2(x(pos), y + 8.0)], egui::Stroke::new(1.5, color));
                                }
                                #[cfg(feature = "gui-test-support")]
                                crate::gui_test_support::register_response(&response, "annotation.comparison.row", if surface == "dna_map" { "window.dna_viewer" } else { "window.splicing_expert" }, Some(&crate::gui_test_support::pseudonymous_subject_scope(&[self.seq_id.as_deref().unwrap_or("unnamed"), &row.structure_id])), crate::gui_test_support::GuiTestWidgetKind::Row, selected);
                                if response.clicked() { self.source_annotation_selected_structure = Some(row.structure_id.clone()); }
                                if response.double_clicked() {
                                    if let Err(error) = self.focus_source_annotation_row(&report, row) { self.splicing_locus_status = error; }
                                }
                                response.on_hover_ui(|ui| {
                                    ui.label(&ids);
                                    ui.label(format!("Genomic strand {}; local strand {}. Full-chain ID: {}", row.genomic_strand, row.local_strand, row.comparison.exon_chain_id));
                                    let spans = |intervals: &[TranscriptInterval]| intervals.iter().map(|i| format!("{}..{}", i.start_1based, i.end_1based)).collect::<Vec<_>>().join(", ");
                                    ui.label(format!("Local exons: {}", spans(&row.exons)));
                                    if let Some(cds) = &row.cds {
                                        ui.label(format!("Local CDS: {}", cds.iter().map(|c| format!("{}..{} (source phase {})", c.interval.start_1based, c.interval.end_1based, c.phase.map(|p| p.to_string()).unwrap_or_else(|| "unassessed".into()))).collect::<Vec<_>>().join(", ")));
                                        ui.small("Phases refer to original source boundaries, not a new frame at a clipped edge.");
                                    } else { ui.label("CDS unassessed"); }
                                    for record in &row.records {
                                        ui.label(format!("{} | designations {:?} | source exon IDs {:?}", record.structure.transcript_id, record.structure.designations, record.structure.exons.iter().map(|e| &e.source_exon_id).collect::<Vec<_>>()));
                                        if let Some(source) = comparison.sources.iter().find(|s| s.source_id == record.structure.source_id) {
                                            ui.label(format!("{:?} {} | {} | {}", source.provider, source.release, source.accession, source.source_id));
                                            ui.monospace(format!("Annotation SHA-256 {}\nLocus SHA-256 {}", source.annotation_sha256, source.locus_sequence_sha256));
                                        }
                                    }
                                    if !row.complete_chain_in_locus { ui.label("Clipped at the locus boundary; not eligible for exact loaded-transcript matching."); }
                                });
                            }
                        });
                });
        });
    }
}
