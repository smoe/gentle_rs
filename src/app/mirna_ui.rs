//! Graphical microRNA target-site scan specialist window.

use super::*;
use crate::mirna::{
    MirnaRegionClass, MirnaSeedClass, MirnaSeedMotif, MirnaTargetHit, MirnaTargetScanReport,
};

#[derive(Clone)]
pub(super) struct MirnaTargetScanPanelState {
    pub(super) show_panel: bool,
    mirna: String,
    mature_sequence: String,
    target: String,
    transcript_filter: String,
    include_3utr: bool,
    include_exons: bool,
    include_introns: bool,
    include_boundaries: bool,
    include_8mer: bool,
    include_7mer_m8: bool,
    include_7mer_a1: bool,
    include_6mer: bool,
    boundary_flank_bp: String,
    species_note: String,
    evidence_note: String,
    comparison_sequences: String,
    selected_hit_index: usize,
    status: String,
    report: Option<MirnaTargetScanReport>,
    error: Option<String>,
}

impl Default for MirnaTargetScanPanelState {
    fn default() -> Self {
        Self {
            show_panel: false,
            mirna: "hsa-miR-96-5p".to_string(),
            mature_sequence: String::new(),
            target: "TP73".to_string(),
            transcript_filter: String::new(),
            include_3utr: true,
            include_exons: true,
            include_introns: true,
            include_boundaries: true,
            include_8mer: true,
            include_7mer_m8: true,
            include_7mer_a1: true,
            include_6mer: true,
            boundary_flank_bp: "25".to_string(),
            species_note: String::new(),
            evidence_note: "rat Tp73 PMID 37099528 orthologous experimental context".to_string(),
            comparison_sequences: String::new(),
            selected_hit_index: 0,
            status: String::new(),
            report: None,
            error: None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ComparisonSequence {
    label: String,
    sequence: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct MirnaComparisonHit {
    label: String,
    seed_class: MirnaSeedClass,
    motif: String,
    start_0based: usize,
}

impl GENtleApp {
    pub(super) fn open_mirna_target_scan_dialog(&mut self) {
        let was_open = self.mirna_panel.show_panel;
        self.mirna_panel.show_panel = true;
        if let Some((seq_id, _)) = self.active_dna_window_context()
            && self.mirna_panel.target.trim().is_empty()
        {
            self.mirna_panel.target = seq_id;
        }
        self.mark_window_open_or_focus(Self::mirna_target_scan_viewport_id(), was_open);
    }

    pub(super) fn render_mirna_target_scan_dialog(&mut self, ctx: &egui::Context) {
        if !self.mirna_panel.show_panel {
            return;
        }
        let mut open = self.mirna_panel.show_panel;
        let viewport_id = Self::mirna_target_scan_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "microRNA Target Scan",
            Self::hosted_mirna_target_scan_window_id(),
            viewport_id,
            Vec2::new(1120.0, 780.0),
            Vec2::new(720.0, 500.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            self.render_mirna_target_scan_contents(ui);
        });
        self.clear_viewport_foreground_request_after_render(viewport_id);
        self.finalize_viewport_open_probe(viewport_id, "microRNA Target Scan");
        self.mirna_panel.show_panel = open && self.mirna_panel.show_panel;
    }

    fn render_mirna_target_scan_contents(&mut self, ui: &mut Ui) {
        if self.render_specialist_window_nav_with_close(
            ui,
            Some(("Close", "Close microRNA Target Scan")),
        ) {
            self.mirna_panel.show_panel = false;
        }
        ui.heading("microRNA Target Scan");
        self.render_mirna_target_scan_inputs(ui);
        ui.separator();
        if let Some(error) = &self.mirna_panel.error {
            ui.colored_label(egui::Color32::RED, error);
        }
        if !self.mirna_panel.status.is_empty() {
            ui.small(self.mirna_panel.status.as_str());
        }
        let report = self.mirna_panel.report.clone();
        if let Some(report) = report {
            self.render_mirna_target_scan_report(ui, &report);
        }
    }

    fn render_mirna_target_scan_inputs(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.label("microRNA");
            ui.text_edit_singleline(&mut self.mirna_panel.mirna);
            ui.label("mature sequence");
            ui.text_edit_singleline(&mut self.mirna_panel.mature_sequence);
        });
        ui.horizontal(|ui| {
            ui.label("target gene or sequence id");
            ui.text_edit_singleline(&mut self.mirna_panel.target);
            if ui.button("Use active sequence").clicked()
                && let Some((seq_id, _)) = self.active_dna_window_context()
            {
                self.mirna_panel.target = seq_id;
            }
            ui.label("transcript filter");
            ui.text_edit_singleline(&mut self.mirna_panel.transcript_filter);
        });
        ui.horizontal_wrapped(|ui| {
            ui.checkbox(&mut self.mirna_panel.include_3utr, "3'UTR");
            ui.checkbox(&mut self.mirna_panel.include_exons, "exons/CDS");
            ui.checkbox(&mut self.mirna_panel.include_introns, "introns");
            ui.checkbox(
                &mut self.mirna_panel.include_boundaries,
                "splice boundaries",
            );
            ui.separator();
            ui.checkbox(&mut self.mirna_panel.include_8mer, "8mer");
            ui.checkbox(&mut self.mirna_panel.include_7mer_m8, "7mer-m8");
            ui.checkbox(&mut self.mirna_panel.include_7mer_a1, "7mer-A1");
            ui.checkbox(&mut self.mirna_panel.include_6mer, "6mer");
            ui.label("boundary flank bp");
            ui.add(
                egui::TextEdit::singleline(&mut self.mirna_panel.boundary_flank_bp)
                    .desired_width(48.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("species note");
            ui.text_edit_singleline(&mut self.mirna_panel.species_note);
        });
        ui.horizontal(|ui| {
            ui.label("evidence note");
            ui.text_edit_singleline(&mut self.mirna_panel.evidence_note);
        });
        egui::CollapsingHeader::new("Ortholog or candidate target snippets")
            .default_open(false)
            .show(ui, |ui| {
                ui.small("Paste FASTA or one raw DNA/RNA sequence per line. These snippets are scanned only for side-by-side conservation/context inspection.");
                ui.add(
                    egui::TextEdit::multiline(&mut self.mirna_panel.comparison_sequences)
                        .desired_rows(5)
                        .code_editor(),
                );
            });
        ui.horizontal(|ui| {
            if ui.button("Scan").clicked() {
                self.run_mirna_target_scan_from_panel();
            }
            if ui.button("Clear result").clicked() {
                self.mirna_panel.report = None;
                self.mirna_panel.error = None;
                self.mirna_panel.status.clear();
            }
        });
    }

    fn run_mirna_target_scan_from_panel(&mut self) {
        let Ok(boundary_flank_bp) = self.mirna_panel.boundary_flank_bp.trim().parse::<usize>()
        else {
            self.mirna_panel.error =
                Some("Boundary flank must be a non-negative integer.".to_string());
            return;
        };
        let mut regions = vec![];
        if self.mirna_panel.include_3utr {
            regions.push(MirnaRegionClass::ThreePrimeUtr);
        }
        if self.mirna_panel.include_exons {
            regions.push(MirnaRegionClass::CodingExon);
            regions.push(MirnaRegionClass::NoncodingExon);
        }
        if self.mirna_panel.include_introns {
            regions.push(MirnaRegionClass::Intron);
        }
        if self.mirna_panel.include_boundaries {
            regions.push(MirnaRegionClass::ExonIntronBoundary);
            regions.push(MirnaRegionClass::IntronExonBoundary);
        }
        let mut seed_classes = vec![];
        if self.mirna_panel.include_8mer {
            seed_classes.push(MirnaSeedClass::EightMer);
        }
        if self.mirna_panel.include_7mer_m8 {
            seed_classes.push(MirnaSeedClass::SevenMerM8);
        }
        if self.mirna_panel.include_7mer_a1 {
            seed_classes.push(MirnaSeedClass::SevenMerA1);
        }
        if self.mirna_panel.include_6mer {
            seed_classes.push(MirnaSeedClass::SixMer);
        }
        if regions.is_empty() || seed_classes.is_empty() {
            self.mirna_panel.error =
                Some("Select at least one region class and one seed class.".to_string());
            return;
        }
        let evidence_notes = self
            .mirna_panel
            .evidence_note
            .trim()
            .is_empty()
            .then(Vec::new)
            .unwrap_or_else(|| vec![self.mirna_panel.evidence_note.trim().to_string()]);
        let command = ShellCommand::MirnaScanTarget {
            mirna: self.mirna_panel.mirna.trim().to_string(),
            target: self.mirna_panel.target.trim().to_string(),
            mature_sequence: nonempty_text(&self.mirna_panel.mature_sequence),
            transcript_filter: nonempty_text(&self.mirna_panel.transcript_filter),
            regions,
            seed_classes,
            boundary_flank_bp,
            species_note: nonempty_text(&self.mirna_panel.species_note),
            evidence_notes,
        };
        let result = self
            .engine
            .write()
            .map_err(|_| "Could not lock GENtle engine for microRNA scan.".to_string())
            .and_then(|mut engine| {
                execute_shell_command_with_options(
                    &mut engine,
                    &command,
                    &ShellExecutionOptions::default(),
                )
            })
            .and_then(|run| {
                serde_json::from_value::<MirnaTargetScanReport>(run.output)
                    .map_err(|error| format!("Could not parse microRNA scan report: {error}"))
            });
        match result {
            Ok(report) => {
                let hit_count = mirna_flat_hits(&report).len();
                self.mirna_panel.selected_hit_index = 0;
                self.mirna_panel.status = format!("Scan complete: {hit_count} candidate hit(s).");
                self.mirna_panel.report = Some(report);
                self.mirna_panel.error = None;
            }
            Err(error) => {
                self.mirna_panel.error = Some(error);
                self.mirna_panel.status.clear();
            }
        }
    }

    fn render_mirna_target_scan_report(&mut self, ui: &mut Ui, report: &MirnaTargetScanReport) {
        ui.horizontal(|ui| {
            ui.label(format!(
                "{} / {}",
                report.query_mirna_id, report.target_gene_or_sequence_id
            ));
            ui.separator();
            ui.small(format!("schema {}", report.schema));
        });
        for warning in &report.warnings {
            ui.colored_label(egui::Color32::from_rgb(150, 96, 0), warning);
        }
        let hits = mirna_flat_hits(report);
        if hits.is_empty() {
            ui.label("No exact seed candidates were found in the selected region classes.");
            self.render_mirna_comparison_panel(ui, report);
            return;
        }
        let selected = self
            .mirna_panel
            .selected_hit_index
            .min(hits.len().saturating_sub(1));
        self.mirna_panel.selected_hit_index = selected;
        ui.columns(2, |columns| {
            columns[0].vertical(|ui| {
                ui.label("Candidate hits");
                egui::ScrollArea::vertical()
                    .max_height(220.0)
                    .show(ui, |ui| {
                        for (idx, hit) in hits.iter().enumerate() {
                            let label = format!(
                                "{} {} {}..{} {}",
                                hit.region_class.as_str(),
                                hit.seed_class.as_str(),
                                hit.genomic_start_1based,
                                hit.genomic_end_1based,
                                hit.transcript_id.as_deref().unwrap_or("-")
                            );
                            if ui.selectable_label(idx == selected, label).clicked() {
                                self.mirna_panel.selected_hit_index = idx;
                            }
                        }
                    });
                self.render_mirna_summary_counts(ui, report);
            });
            columns[1].vertical(|ui| {
                self.render_mirna_comparison_panel(ui, report);
            });
        });
        ui.separator();
        let hit = hits[self
            .mirna_panel
            .selected_hit_index
            .min(hits.len().saturating_sub(1))];
        let motif = report
            .seed_motif_table
            .iter()
            .find(|motif| motif.seed_class == hit.seed_class);
        if let Some(motif) = motif {
            render_mirna_pairing_graphic(ui, report, hit, motif);
        }
        ui.separator();
        ui.label("Splicing and transcript-context interpretation");
        ui.small(mirna_region_discussion(hit.region_class));
        ui.small("This is a deterministic sequence-evidence view. Conservation or region class can prioritize follow-up, but neither proves repression, altered splicing, nor direct AGO occupancy.");
    }

    fn render_mirna_summary_counts(&self, ui: &mut Ui, report: &MirnaTargetScanReport) {
        egui::CollapsingHeader::new("Summary counts")
            .default_open(false)
            .show(ui, |ui| {
                egui::Grid::new("mirna_summary_counts_grid")
                    .striped(true)
                    .show(ui, |ui| {
                        ui.strong("region");
                        ui.strong("seed");
                        ui.strong("count");
                        ui.end_row();
                        for count in &report.summary_counts {
                            ui.label(count.region_class.as_str());
                            ui.label(count.seed_class.as_str());
                            ui.label(count.count.to_string());
                            ui.end_row();
                        }
                    });
            });
    }

    fn render_mirna_comparison_panel(&self, ui: &mut Ui, report: &MirnaTargetScanReport) {
        ui.label("Ortholog / candidate snippets");
        let comparison_hits = mirna_comparison_hits(
            &self.mirna_panel.comparison_sequences,
            &report.seed_motif_table,
        );
        if comparison_hits.is_empty() {
            ui.small(
                "Paste ortholog or candidate target snippets above to compare motif conservation.",
            );
            return;
        }
        egui::ScrollArea::vertical()
            .max_height(220.0)
            .show(ui, |ui| {
                for hit in comparison_hits.iter().take(80) {
                    ui.horizontal(|ui| {
                        ui.colored_label(seed_class_color(hit.seed_class), hit.seed_class.as_str());
                        ui.monospace(format!("{}:{} {}", hit.label, hit.start_0based, hit.motif));
                    });
                }
                if comparison_hits.len() > 80 {
                    ui.small(format!(
                        "{} additional comparison hit(s) hidden.",
                        comparison_hits.len() - 80
                    ));
                }
            });
    }
}

fn render_mirna_pairing_graphic(
    ui: &mut Ui,
    report: &MirnaTargetScanReport,
    hit: &MirnaTargetHit,
    motif: &MirnaSeedMotif,
) {
    let desired = Vec2::new(ui.available_width().max(420.0), 260.0);
    let (rect, _) = ui.allocate_exact_size(desired, egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 6.0, egui::Color32::from_rgb(248, 249, 250));
    painter.rect_stroke(
        rect,
        6.0,
        egui::Stroke::new(1.0, egui::Color32::from_gray(190)),
        egui::StrokeKind::Outside,
    );
    let mature = report.mature_sequence.replace('T', "U");
    let bases: Vec<char> = mature.chars().collect();
    let (seed_start_1based, seed_end_1based) = seed_bounds(motif);
    let cell_w = ((rect.width() - 40.0) / bases.len().max(1) as f32).clamp(14.0, 28.0);
    let x0 = rect.left() + 20.0;
    let mirna_y = rect.top() + 44.0;
    painter.text(
        Pos2::new(x0, rect.top() + 14.0),
        egui::Align2::LEFT_CENTER,
        "microRNA 5'->3'",
        egui::FontId::monospace(13.0),
        egui::Color32::from_gray(50),
    );
    for (idx, base) in bases.iter().enumerate() {
        let pos_1based = idx + 1;
        let x = x0 + idx as f32 * cell_w;
        let cell = egui::Rect::from_min_size(Pos2::new(x, mirna_y), Vec2::new(cell_w - 2.0, 26.0));
        let in_seed = pos_1based >= seed_start_1based && pos_1based <= seed_end_1based;
        painter.rect_filled(
            cell,
            3.0,
            if in_seed {
                egui::Color32::from_rgb(255, 231, 166)
            } else {
                egui::Color32::WHITE
            },
        );
        painter.rect_stroke(
            cell,
            3.0,
            egui::Stroke::new(1.0, egui::Color32::from_gray(170)),
            egui::StrokeKind::Outside,
        );
        painter.text(
            cell.center(),
            egui::Align2::CENTER_CENTER,
            base,
            egui::FontId::monospace(15.0),
            egui::Color32::BLACK,
        );
    }
    let target_y = mirna_y + 92.0;
    let motif_chars: Vec<char> = motif.target_motif.chars().collect();
    painter.text(
        Pos2::new(x0, target_y - 28.0),
        egui::Align2::LEFT_CENTER,
        format!(
            "target {} {}..{} ({})",
            hit.region_class.as_str(),
            hit.genomic_start_1based,
            hit.genomic_end_1based,
            hit.strand
        ),
        egui::FontId::monospace(13.0),
        egui::Color32::from_gray(50),
    );
    let target_x0 = x0 + seed_start_1based.saturating_sub(1) as f32 * cell_w;
    let a1_anchor_index = mirna_a1_anchor_target_index(motif);
    for (idx, base) in motif_chars.iter().enumerate() {
        let x = target_x0 + idx as f32 * cell_w;
        let cell = egui::Rect::from_min_size(Pos2::new(x, target_y), Vec2::new(cell_w - 2.0, 26.0));
        let is_a1_anchor = Some(idx) == a1_anchor_index;
        painter.rect_filled(
            cell,
            3.0,
            if is_a1_anchor {
                egui::Color32::from_rgb(95, 101, 114)
            } else {
                seed_class_color(hit.seed_class)
            },
        );
        painter.text(
            cell.center(),
            egui::Align2::CENTER_CENTER,
            base,
            egui::FontId::monospace(15.0),
            egui::Color32::WHITE,
        );
        if is_a1_anchor {
            painter.text(
                Pos2::new(cell.center().x, cell.bottom() + 12.0),
                egui::Align2::CENTER_CENTER,
                "A1",
                egui::FontId::monospace(10.0),
                egui::Color32::from_gray(80),
            );
        }
    }
    for (mirna_col, target_col) in mirna_seed_pair_columns(motif) {
        let top_x = x0 + mirna_col as f32 * cell_w + cell_w * 0.5;
        let bottom_x = target_x0 + target_col as f32 * cell_w + cell_w * 0.5;
        painter.line_segment(
            [
                Pos2::new(top_x, mirna_y + 28.0),
                Pos2::new(bottom_x, target_y),
            ],
            egui::Stroke::new(1.0, egui::Color32::from_gray(80)),
        );
    }
    if let Some(target_col) = a1_anchor_index {
        let top_x = x0 + cell_w * 0.5;
        let bottom_x = target_x0 + target_col as f32 * cell_w + cell_w * 0.5;
        painter.line_segment(
            [
                Pos2::new(top_x, mirna_y + 28.0),
                Pos2::new(bottom_x, target_y),
            ],
            egui::Stroke::new(1.0, egui::Color32::from_gray(150)),
        );
    }
    let context_text = format!("context: {}", hit.region_context_sequence);
    painter.text(
        Pos2::new(x0, target_y + 48.0),
        egui::Align2::LEFT_CENTER,
        context_text,
        egui::FontId::monospace(12.0),
        egui::Color32::from_gray(40),
    );
    let bar = egui::Rect::from_min_size(
        Pos2::new(x0, rect.bottom() - 36.0),
        Vec2::new(rect.width() - 40.0, 10.0),
    );
    painter.rect_filled(bar, 4.0, egui::Color32::from_gray(220));
    let marker_x = bar.left() + bar.width() * 0.5;
    painter.circle_filled(
        Pos2::new(marker_x, bar.center().y),
        6.0,
        region_class_color(hit.region_class),
    );
    painter.text(
        Pos2::new(x0, rect.bottom() - 14.0),
        egui::Align2::LEFT_CENTER,
        format!(
            "local 0-based [{}..{}), matched {}",
            hit.local_start_0based, hit.local_end_0based_exclusive, hit.matched_sequence
        ),
        egui::FontId::monospace(12.0),
        egui::Color32::from_gray(50),
    );
}

fn nonempty_text(value: &str) -> Option<String> {
    let trimmed = value.trim();
    (!trimmed.is_empty()).then(|| trimmed.to_string())
}

fn mirna_flat_hits(report: &MirnaTargetScanReport) -> Vec<&MirnaTargetHit> {
    report
        .grouped_hits
        .iter()
        .flat_map(|group| group.hits.iter())
        .collect()
}

fn seed_bounds(motif: &MirnaSeedMotif) -> (usize, usize) {
    let Some((start, end)) = motif.mirna_positions_1based.split_once("..") else {
        return (2, 8);
    };
    (
        start.trim().parse::<usize>().unwrap_or(2),
        end.trim().parse::<usize>().unwrap_or(8),
    )
}

fn mirna_a1_anchor_target_index(motif: &MirnaSeedMotif) -> Option<usize> {
    matches!(
        motif.seed_class,
        MirnaSeedClass::EightMer | MirnaSeedClass::SevenMerA1
    )
    .then(|| motif.target_motif.len().saturating_sub(1))
}

pub(super) fn mirna_seed_pair_columns(motif: &MirnaSeedMotif) -> Vec<(usize, usize)> {
    let (seed_start_1based, seed_end_1based) = seed_bounds(motif);
    let seed_len = seed_end_1based.saturating_sub(seed_start_1based) + 1;
    let paired_target_len = motif
        .target_motif
        .len()
        .saturating_sub(usize::from(mirna_a1_anchor_target_index(motif).is_some()));
    let pair_count = seed_len.min(paired_target_len);
    (0..pair_count)
        .map(|idx| {
            (
                seed_start_1based.saturating_sub(1) + idx,
                paired_target_len.saturating_sub(1 + idx),
            )
        })
        .collect()
}

fn seed_class_color(seed_class: MirnaSeedClass) -> egui::Color32 {
    match seed_class {
        MirnaSeedClass::EightMer => egui::Color32::from_rgb(120, 58, 180),
        MirnaSeedClass::SevenMerM8 => egui::Color32::from_rgb(32, 116, 188),
        MirnaSeedClass::SevenMerA1 => egui::Color32::from_rgb(23, 133, 108),
        MirnaSeedClass::SixMer => egui::Color32::from_rgb(132, 120, 42),
    }
}

fn region_class_color(region_class: MirnaRegionClass) -> egui::Color32 {
    match region_class {
        MirnaRegionClass::ThreePrimeUtr => egui::Color32::from_rgb(55, 143, 74),
        MirnaRegionClass::CodingExon => egui::Color32::from_rgb(43, 102, 190),
        MirnaRegionClass::NoncodingExon => egui::Color32::from_rgb(39, 142, 154),
        MirnaRegionClass::Intron => egui::Color32::from_rgb(128, 91, 168),
        MirnaRegionClass::ExonIntronBoundary | MirnaRegionClass::IntronExonBoundary => {
            egui::Color32::from_rgb(190, 103, 32)
        }
        MirnaRegionClass::WholeTranscript | MirnaRegionClass::WholeGene => {
            egui::Color32::from_gray(90)
        }
    }
}

pub(super) fn mirna_region_discussion(region_class: MirnaRegionClass) -> &'static str {
    match region_class {
        MirnaRegionClass::ThreePrimeUtr => {
            "3'UTR candidates match the canonical model, but transcript-specific UTR usage and alternative polyadenylation determine whether the site is present."
        }
        MirnaRegionClass::CodingExon => {
            "CDS candidates can affect translation or mRNA stability and may also overlap exonic splicing enhancers/silencers; inspect codon phase and transcript usage before treating this as a simple repression site."
        }
        MirnaRegionClass::NoncodingExon => {
            "Noncoding exon candidates may behave like UTR sites and can also intersect exon definition signals, so splice-isoform-specific inclusion matters."
        }
        MirnaRegionClass::Intron => {
            "Intronic candidates are pre-mRNA-context evidence; consider nuclear AGO/RBP competition, intron retention, and nearby branch/acceptor/donor features rather than only mature 3'UTR repression."
        }
        MirnaRegionClass::ExonIntronBoundary => {
            "Exon-to-intron boundary candidates sit near splice donor neighborhoods and may affect exon definition or splice-site accessibility if the site is present in the relevant pre-mRNA."
        }
        MirnaRegionClass::IntronExonBoundary => {
            "Intron-to-exon boundary candidates sit near splice acceptor neighborhoods and should be interpreted together with exon phase, branchpoint distance, and transcript-specific exon inclusion."
        }
        MirnaRegionClass::WholeTranscript | MirnaRegionClass::WholeGene => {
            "Whole-record scans are discovery views; assign each hit back to transcript structure before making a mechanistic interpretation."
        }
    }
}

fn parse_comparison_sequences(raw: &str) -> Vec<ComparisonSequence> {
    let mut records = vec![];
    let mut label: Option<String> = None;
    let mut seq = String::new();
    for line in raw.lines().map(str::trim).filter(|line| !line.is_empty()) {
        if let Some(rest) = line.strip_prefix('>') {
            if !seq.is_empty() {
                records.push(ComparisonSequence {
                    label: label
                        .clone()
                        .unwrap_or_else(|| format!("candidate_{}", records.len() + 1)),
                    sequence: normalize_candidate_sequence(&seq),
                });
                seq.clear();
            }
            label = Some(rest.trim().to_string());
        } else if label.is_none() && !seq.is_empty() {
            records.push(ComparisonSequence {
                label: format!("candidate_{}", records.len() + 1),
                sequence: normalize_candidate_sequence(&seq),
            });
            seq = line.to_string();
        } else {
            seq.push_str(line);
        }
    }
    if !seq.is_empty() {
        records.push(ComparisonSequence {
            label: label.unwrap_or_else(|| format!("candidate_{}", records.len() + 1)),
            sequence: normalize_candidate_sequence(&seq),
        });
    }
    records
        .into_iter()
        .filter(|record| !record.sequence.is_empty())
        .collect()
}

fn normalize_candidate_sequence(raw: &str) -> String {
    raw.chars()
        .filter_map(|c| match c.to_ascii_uppercase() {
            'A' | 'C' | 'G' | 'T' => Some(c.to_ascii_uppercase()),
            'U' => Some('T'),
            _ => None,
        })
        .collect()
}

pub(super) fn mirna_comparison_hits(
    raw: &str,
    motifs: &[MirnaSeedMotif],
) -> Vec<MirnaComparisonHit> {
    let mut hits = vec![];
    for record in parse_comparison_sequences(raw) {
        for motif in motifs {
            let mut offset = 0usize;
            while let Some(pos) = record.sequence[offset..].find(&motif.target_motif) {
                let start = offset + pos;
                hits.push(MirnaComparisonHit {
                    label: record.label.clone(),
                    seed_class: motif.seed_class,
                    motif: motif.target_motif.clone(),
                    start_0based: start,
                });
                offset = start + 1;
            }
        }
    }
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    fn motif_with_positions(
        seed_class: MirnaSeedClass,
        target_motif: &str,
        positions: &str,
    ) -> MirnaSeedMotif {
        MirnaSeedMotif {
            seed_class,
            mirna_positions_1based: positions.to_string(),
            mirna_seed_5p: "TTGGCAC".to_string(),
            target_motif: target_motif.to_string(),
        }
    }

    fn motif(seed_class: MirnaSeedClass, target_motif: &str) -> MirnaSeedMotif {
        motif_with_positions(seed_class, target_motif, "2..8")
    }

    #[test]
    fn mirna_region_discussion_keeps_splicing_context_visible() {
        assert!(mirna_region_discussion(MirnaRegionClass::CodingExon).contains("splicing"));
        assert!(mirna_region_discussion(MirnaRegionClass::IntronExonBoundary).contains("splice"));
        assert!(mirna_region_discussion(MirnaRegionClass::Intron).contains("pre-mRNA"));
    }

    #[test]
    fn comparison_hits_scan_fasta_and_raw_candidate_sequences() {
        let hits = mirna_comparison_hits(
            ">rat_tp73\nAACCTGCCAAATTT\n>candidate_2\nGGGGTGCCAAAG\n",
            &[
                motif(MirnaSeedClass::SevenMerA1, "TGCCAAA"),
                motif(MirnaSeedClass::EightMer, "GTGCCAAA"),
            ],
        );
        assert!(
            hits.iter()
                .any(|hit| hit.label == "rat_tp73" && hit.start_0based == 4)
        );
        assert!(
            hits.iter()
                .any(|hit| hit.label == "candidate_2"
                    && hit.seed_class == MirnaSeedClass::EightMer)
        );
    }

    #[test]
    fn seed_pair_columns_keep_a1_anchor_out_of_watson_crick_seed_pairs() {
        assert_eq!(
            mirna_seed_pair_columns(&motif(MirnaSeedClass::EightMer, "GTGCCAAA")),
            vec![(1, 6), (2, 5), (3, 4), (4, 3), (5, 2), (6, 1), (7, 0)]
        );
        assert_eq!(
            mirna_seed_pair_columns(&motif_with_positions(
                MirnaSeedClass::SevenMerA1,
                "TGCCAAA",
                "2..7"
            )),
            vec![(1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
        );
        assert_eq!(
            mirna_seed_pair_columns(&motif(MirnaSeedClass::SevenMerM8, "GTGCCAA")),
            vec![(1, 6), (2, 5), (3, 4), (4, 3), (5, 2), (6, 1), (7, 0)]
        );
        assert_eq!(
            mirna_seed_pair_columns(&motif_with_positions(
                MirnaSeedClass::SixMer,
                "TGCCAA",
                "2..7"
            )),
            vec![(1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
        );
    }
}
