//! Read-only, snapshot-bound review of retained specificity receipts in PCR Designer.
//! Tests construct synthetic presentation records inline; they are not search or assay evidence.

use super::{MainAreaDna, transcript_assay_report_ui::page};
use crate::engine::{
    Operation, TranscriptAssayPanelReport, TranscriptAssayPanelSpecificityMatrix as Matrix,
    TranscriptAssaySpecificityMatrixCell as Cell,
    TranscriptAssaySpecificityMatrixSelection as Receipt,
};
use crate::genomes::BlastDatabaseIndexKind;
use eframe::egui;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::sync::{Arc, Mutex, Weak, mpsc};
use std::time::Duration;

type Reply = Result<Arc<Matrix>, String>;
type Inbox = Mutex<Option<mpsc::Receiver<Reply>>>;

fn hosted_review_id() -> egui::Id {
    egui::Id::unique("pcr_designer_specificity_inbox")
}

#[derive(Debug)]
struct Binding {
    project: (u64, u64),
    report: Arc<TranscriptAssayPanelReport>,
}

#[derive(Debug)]
struct Task {
    binding: Binding,
    receiver: Arc<Inbox>,
}

#[derive(Debug, Default)]
pub(super) struct Review {
    task: Option<Task>,
    accepted: Option<(Binding, Arc<Matrix>)>,
    status: String,
}

impl Clone for Review {
    // A cloned DNA view must not consume another view's one-shot reply.
    fn clone(&self) -> Self {
        Self::default()
    }
}

impl MainAreaDna {
    // PCR Designer borrows its DNA view. Closing just the specialist must also
    // disconnect its pending reply, without locking the DNA window or the project.
    pub(crate) fn disconnect_hosted_specificity_review(ctx: &egui::Context) {
        ctx.data_mut(|data| {
            let id = hosted_review_id();
            if let Some(weak) = data.remove_temp::<Weak<Inbox>>(id)
                && let Some(inbox) = weak.upgrade()
            {
                match inbox.try_lock() {
                    Ok(mut slot) => {
                        *slot = None;
                    }
                    Err(std::sync::TryLockError::WouldBlock) => {
                        data.insert_temp(id, weak);
                    }
                    Err(std::sync::TryLockError::Poisoned(e)) => {
                        *e.into_inner() = None;
                    }
                }
            }
        });
    }

    pub(super) fn invalidate_transcript_assay_specificity_review(&mut self) {
        let review = &mut self.transcript_assay_specificity_review;
        let had_review = review.task.is_some() || review.accepted.is_some();
        *review = Review::default();
        if had_review {
            review.status = "Specificity review dropped: displayed panel changed".into();
        }
    }

    fn start_transcript_assay_specificity_review(&mut self, ctx: &egui::Context) {
        let review = &mut self.transcript_assay_specificity_review;
        if review.task.is_some() {
            return;
        }
        let Some(report) = self.cached_transcript_assay_panel_report.clone() else {
            return;
        };
        let Some(engine) = self.engine.clone() else {
            review.status = "No project engine".into();
            return;
        };
        let project = match engine.try_read() {
            Ok(g) => (g.instance_id(), g.structural_revision()),
            Err(_) => {
                review.status = "Project busy; retry shortly".into();
                return;
            }
        };
        let operation = Operation::InspectTranscriptAssayPanelSpecificityMatrix {
            panel_report_id: report.report_id.clone(),
        };
        let (tx, rx) = mpsc::channel();
        review.accepted = None;
        match std::thread::Builder::new()
            .name("gentle-panel-specificity-review".into())
            .stack_size(16 * 1024 * 1024)
            .spawn(move || {
                let reply =
                    crate::background_engine::execute_read_only_operation_on_engine_snapshot(
                        &engine, operation,
                    )
                    .map_err(|e| e.message)
                    .and_then(|result| {
                        result
                            .transcript_assay_specificity_matrix
                            .map(Arc::from)
                            .ok_or_else(|| "Inspection returned no specificity matrix".to_string())
                    });
                let _ = tx.send(reply);
            }) {
            Ok(_) => {
                let receiver = Arc::new(Mutex::new(Some(rx)));
                if self.rendering_pcr_specialist {
                    Self::disconnect_hosted_specificity_review(ctx);
                    ctx.data_mut(|d| d.insert_temp(hosted_review_id(), Arc::downgrade(&receiver)));
                }
                review.task = Some(Task {
                    binding: Binding { project, report },
                    receiver,
                });
                review.status = "Reading retained specificity receipts...".into();
            }
            Err(e) => review.status = format!("Could not start specificity review: {e}"),
        }
    }

    pub(super) fn poll_transcript_assay_specificity_review(&mut self, ctx: &egui::Context) {
        let review = &mut self.transcript_assay_specificity_review;
        let Some(task) = &review.task else { return };
        // Check before receiving: a completed reply stays in the channel during contention.
        let project = match &self.engine {
            Some(engine) => match engine.try_read() {
                Ok(g) => Some((g.instance_id(), g.structural_revision())),
                Err(_) => {
                    ctx.request_repaint_after(Duration::from_millis(100));
                    return;
                }
            },
            None => None,
        };
        let reply = match task.receiver.try_lock() {
            Ok(slot) => match slot.as_ref().map(mpsc::Receiver::try_recv) {
                None => Err("Specificity review dropped: PCR Designer closed".into()),
                Some(Ok(reply)) => reply,
                Some(Err(mpsc::TryRecvError::Empty)) => {
                    ctx.request_repaint_after(Duration::from_millis(100));
                    return;
                }
                Some(Err(mpsc::TryRecvError::Disconnected)) => {
                    Err("Specificity review worker disconnected".into())
                }
            },
            Err(std::sync::TryLockError::WouldBlock) => {
                ctx.request_repaint_after(Duration::from_millis(100));
                return;
            }
            Err(_) => Err("Specificity review channel unavailable".into()),
        };
        let binding = review.task.take().unwrap().binding;
        review.status = if project != Some(binding.project) {
            "Specificity review dropped: project changed".into()
        } else if !self
            .cached_transcript_assay_panel_report
            .as_ref()
            .is_some_and(|r| Arc::ptr_eq(r, &binding.report))
        {
            "Specificity review dropped: displayed panel changed".into()
        } else {
            match reply {
                Ok(matrix) if matrix.panel_report_id == binding.report.report_id => {
                    review.accepted = Some((binding, matrix));
                    "Retained specificity matrix".into()
                }
                Ok(_) => "Specificity review dropped: matrix report ID does not match the displayed panel".into(),
                Err(message) => message,
            }
        };
        ctx.request_repaint();
    }

    pub(super) fn specificity_review_button(&mut self, ui: &mut egui::Ui) {
        let enabled = self.transcript_assay_specificity_review.task.is_none()
            && self
                .cached_transcript_assay_panel_report
                .as_ref()
                .is_some_and(|r| !r.report_id.trim().is_empty());
        if ui
            .add_enabled(enabled, egui::Button::new("Review specificity"))
            .clicked()
        {
            self.start_transcript_assay_specificity_review(ui.ctx());
            ui.ctx().request_repaint();
        }
    }

    pub(super) fn show_transcript_assay_specificity_review(&mut self, ui: &mut egui::Ui) {
        let review = &mut self.transcript_assay_specificity_review;
        if let Some((binding, _)) = &review.accepted {
            let stale = match &self.engine {
                Some(engine) => match engine.try_read() {
                    Ok(g) => binding.project != (g.instance_id(), g.structural_revision()),
                    Err(_) => {
                        ui.ctx().request_repaint_after(Duration::from_millis(100));
                        false
                    }
                },
                None => true,
            };
            if stale {
                review.accepted = None;
                review.status = "Project changed since review; review again".into();
            } else if !self
                .cached_transcript_assay_panel_report
                .as_ref()
                .is_some_and(|r| Arc::ptr_eq(r, &binding.report))
            {
                review.accepted = None;
                review.status = "Specificity review dropped: displayed panel changed".into();
            }
        }
        if !review.status.is_empty() {
            ui.label(&review.status);
        }
        if let Some((_, matrix)) = &review.accepted {
            render_transcript_assay_specificity_matrix(ui, matrix);
        }
    }
}

const RECEIPTS_PER_PAGE: usize = 6;
const ROWS_PER_PAGE: usize = 24;

fn field(ui: &mut egui::Ui, name: &str, value: impl std::fmt::Display) {
    ui.label(format!("{name}: {value}"));
}

fn render_transcript_assay_specificity_matrix(ui: &mut egui::Ui, matrix: &Matrix) {
    field(ui, "Panel report", &matrix.panel_report_id);
    field(
        ui,
        "Current panel design digest",
        &matrix.current_panel_design_digest,
    );
    for line in &matrix.nonclaims {
        ui.small(line);
    }

    for (kind, name, unassessed) in [
        (
            BlastDatabaseIndexKind::GenomicDna,
            "Genomic",
            &matrix.unassessed_genomic_assay_ids,
        ),
        (
            BlastDatabaseIndexKind::TranscriptomeCdna,
            "Transcriptome",
            &matrix.unassessed_transcriptome_assay_ids,
        ),
    ] {
        ui.push_id((&matrix.panel_report_id, name), |ui| {
            render_dimension(ui, matrix, kind, name, unassessed);
        });
    }
    ui.strong("Same-oligo citations, not panel evidence");
    let range = page(
        ui,
        &matrix.panel_report_id,
        "Same-oligo citations",
        matrix.unlinked_pair_summaries.len(),
        ROWS_PER_PAGE,
    );
    for citation in &matrix.unlinked_pair_summaries[range] {
        field(ui, "Summary", &citation.summary_id);
        field(ui, "Pair binding", &citation.pair_binding_sha256);
        field(
            ui,
            "Same-oligo current assays",
            citation.same_oligo_current_assay_ids.join(", "),
        );
        ui.label(&citation.reason);
    }
}

fn render_dimension(
    ui: &mut egui::Ui,
    matrix: &Matrix,
    kind: BlastDatabaseIndexKind,
    name: &str,
    unassessed: &[String],
) {
    let id = &matrix.panel_report_id;
    let receipts: Vec<_> = matrix
        .selections
        .iter()
        .filter(|s| s.reference.index_kind == kind)
        .collect();
    ui.separator();
    ui.strong(format!("{name} specificity"));
    let gaps = if unassessed.len() > ROWS_PER_PAGE {
        page(
            ui,
            id,
            &format!("{name} unassessed assays"),
            unassessed.len(),
            ROWS_PER_PAGE,
        )
    } else {
        0..unassessed.len()
    };
    ui.label(format!(
        "Unassessed {} assays: {}",
        name.to_ascii_lowercase(),
        if unassessed.is_empty() {
            "none listed".into()
        } else {
            unassessed[gaps].join(", ")
        }
    ));
    if receipts.is_empty() {
        ui.label(format!("{name}: no retained receipts"));
    } else {
        ui.small("Receipt order is not chronological.");
    }
    let cols = page(
        ui,
        id,
        &format!("{name} receipts"),
        receipts.len(),
        RECEIPTS_PER_PAGE,
    );
    let current = page(
        ui,
        id,
        &format!("{name} assays"),
        matrix.assays.len(),
        ROWS_PER_PAGE,
    );
    let current_ids: HashSet<_> = matrix.assays.iter().map(|a| a.assay_id.as_str()).collect();
    let history: Vec<_> = receipts
        .iter()
        .flat_map(|r| &r.cells)
        .map(|c| c.assay_id.as_str())
        .filter(|assay| !current_ids.contains(assay))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let old = page(
        ui,
        id,
        &format!("{name} retained-only assays"),
        history.len(),
        ROWS_PER_PAGE,
    );
    let visible = &receipts[cols];
    let lookup: Vec<HashMap<_, _>> = visible
        .iter()
        .map(|s| s.cells.iter().map(|c| (c.assay_id.as_str(), c)).collect())
        .collect();
    egui::ScrollArea::horizontal()
        .id_salt("matrix-scroll")
        .show(ui, |ui| {
            egui::Grid::new("receipt-grid")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("Assay ID");
                    for receipt in visible {
                        ui.vertical(|ui| receipt_header(ui, receipt));
                    }
                    ui.end_row();
                    for (heading, ids) in [
                        (
                            "Current saved design",
                            matrix.assays[current]
                                .iter()
                                .map(|a| a.assay_id.as_str())
                                .collect::<Vec<_>>(),
                        ),
                        (
                            "Only in retained receipts (not in the saved design)",
                            history[old].to_vec(),
                        ),
                    ] {
                        ui.strong(heading);
                        ui.end_row();
                        for assay in ids {
                            ui.label(assay);
                            for column in &lookup {
                                ui.vertical(|ui| render_cell(ui, column.get(assay).copied()));
                            }
                            ui.end_row();
                        }
                    }
                });
        });
    for receipt in visible {
        ui.push_id(&receipt.selection_id, |ui| {
            ui.collapsing(format!("Receipt details: {}", receipt.selection_id), |ui| {
                field(ui, "Reference", serde_json::json!(receipt.reference));
                field(ui, "Policy", serde_json::json!(receipt.policy));
                field(ui, "Selection ID", &receipt.selection_id);
                field(ui, "Acceptance ID", &receipt.acceptance_id);
                field(ui, "Panel design digest", &receipt.panel_design_digest);
                for diagnostic in &receipt.diagnostics {
                    ui.label(diagnostic);
                }
                let range = page(ui, id, "Receipt cells", receipt.cells.len(), ROWS_PER_PAGE);
                for cell in &receipt.cells[range] {
                    ui.strong(&cell.assay_id);
                    field(
                        ui,
                        "Report ID",
                        cell.report_id.as_deref().unwrap_or("not recorded"),
                    );
                    field(
                        ui,
                        "Report content SHA-256",
                        cell.report_content_sha256
                            .as_deref()
                            .unwrap_or("not recorded"),
                    );
                    for diagnostic in &cell.diagnostics {
                        ui.label(diagnostic);
                    }
                }
            });
        });
    }
}

fn yes_no(value: bool) -> &'static str {
    if value { "yes" } else { "no" }
}

fn receipt_header(ui: &mut egui::Ui, receipt: &Receipt) {
    let r = &receipt.reference;
    let mut chars = r.content_fingerprint.chars();
    let mut brief: String = chars.by_ref().take(18).collect();
    if chars.next().is_some() {
        brief.push_str("...");
    }
    ui.set_max_width(220.);
    ui.strong(&r.genome_id);
    ui.label(&receipt.selection_id);
    ui.label(r.index_kind.as_str());
    field(
        ui,
        "Assembly",
        r.assembly.as_deref().unwrap_or("not recorded"),
    );
    field(
        ui,
        "Release",
        r.release.as_deref().unwrap_or("not recorded"),
    );
    ui.label(format!("{}: {brief}", r.fingerprint_algorithm))
        .on_hover_text(&r.content_fingerprint);
    ui.label(format!(
        "Active: {} | Superseded: {}",
        yes_no(receipt.is_active),
        yes_no(receipt.is_superseded)
    ));
    field(
        ui,
        "Replaces",
        receipt.replaces_selection_id.as_deref().unwrap_or("none"),
    );
    field(
        ui,
        "Matches current panel design",
        yes_no(receipt.matches_current_panel_design),
    );
    field(
        ui,
        "Recorded acceptance",
        receipt.recorded_acceptance_status.as_str(),
    );
}

fn render_cell(ui: &mut egui::Ui, cell: Option<&Cell>) {
    let Some(cell) = cell else {
        ui.weak("not recorded in this receipt");
        return;
    };
    if cell
        .diagnostics
        .iter()
        .any(|d| d == "assay_not_in_retained_selection")
    {
        ui.weak("not in receipt: unassessed");
        return;
    }
    if !cell.retained_binding_verified {
        ui.strong("Unverified binding");
    }
    if !cell.matches_current_assay {
        ui.strong("Not a current assay match");
    }
    let verdict = serde_json::json!(cell.verdict);
    field(
        ui,
        "Verdict",
        verdict.as_str().expect("serialized verdict is a string"),
    );
    field(
        ui,
        "Recorded status",
        cell.recorded_status.as_deref().unwrap_or("not recorded"),
    );
    field(
        ui,
        "Retained binding verified",
        yes_no(cell.retained_binding_verified),
    );
    field(
        ui,
        "Matches current assay",
        yes_no(cell.matches_current_assay),
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence::DNAsequence;
    use crate::engine::{
        GentleEngine, PRIMER_DESIGN_REPORTS_METADATA_KEY, PrimerDesignPairRecord,
        PrimerDesignPrimerRecord, PrimerSpecificityMultiVerdict,
        PrimerSpecificityReferenceIdentity, TranscriptAssayPanelAssay,
        TranscriptAssayPanelSpecificityAcceptanceStatus, TranscriptAssaySpecificityMatrixAssay,
        TranscriptAssaySpecificityMatrixUnlinkedSummary,
    };
    use std::sync::RwLock;
    use std::time::Instant;

    const RAW_ID: &str = "P '\u{03b2}";
    const PANEL: &str = "P___";

    fn matrix(assays: usize, per_kind: usize, old: usize) -> Matrix {
        let assays: Vec<_> = (0..assays)
            .map(|i| TranscriptAssaySpecificityMatrixAssay {
                assay_id: format!("assay-{i:02}"),
                assay_rank: i + 1,
                pair_binding_sha256: format!("pair-{i}"),
            })
            .collect();
        let mut selections = Vec::new();
        for kind in [
            BlastDatabaseIndexKind::GenomicDna,
            BlastDatabaseIndexKind::TranscriptomeCdna,
        ] {
            for i in 0..per_kind {
                let cells = assays
                    .iter()
                    .map(|a| a.assay_id.clone())
                    .chain((0..old).rev().map(|i| format!("old-{i:02}")))
                    .map(|id| Cell {
                        matches_current_assay: i != 0 && !id.starts_with("old-"),
                        assay_id: id,
                        report_id: Some("child-report".into()),
                        report_content_sha256: Some("synthetic-child-hash".into()),
                        recorded_status: Some("pass".into()),
                        verdict: if i == 0 {
                            PrimerSpecificityMultiVerdict::Incomplete
                        } else {
                            PrimerSpecificityMultiVerdict::Pass
                        },
                        retained_binding_verified: i != 0,
                        diagnostics: if i == 0 {
                            vec!["synthetic_cell_diagnostic".into()]
                        } else {
                            vec![]
                        },
                    })
                    .collect();
                selections.push(Receipt {
                    selection_id: format!("{}-{i:02}", kind.as_str()),
                    panel_design_digest: format!("design-{i}"),
                    reference: PrimerSpecificityReferenceIdentity {
                        genome_id: "same-genome".into(),
                        index_kind: kind,
                        assembly: Some("assembly-A".into()),
                        release: Some("release-B".into()),
                        fingerprint_algorithm: "synthetic-algorithm".into(),
                        content_fingerprint: format!("digest-{i:02}-01234567890123456789"),
                        subject_annotation_fingerprint: Some("annotation-fingerprint".into()),
                        subject_annotation_fingerprint_algorithm: Some(
                            "annotation-algorithm".into(),
                        ),
                    },
                    policy: Default::default(),
                    acceptance_id: format!("accept-{i}"),
                    recorded_acceptance_status:
                        TranscriptAssayPanelSpecificityAcceptanceStatus::Pass,
                    is_active: i == 0,
                    is_superseded: i == 0,
                    replaces_selection_id: Some("previous-selection".into()),
                    matches_current_panel_design: i != 0,
                    cells,
                    diagnostics: vec!["synthetic_receipt_diagnostic".into()],
                });
            }
        }
        Matrix {
            schema: "gentle.transcript_assay_panel_specificity_matrix.v1".into(),
            panel_report_id: PANEL.into(),
            current_panel_design_digest: "current-digest".into(),
            assays,
            selections,
            unassessed_genomic_assay_ids: vec!["genomic-gap".into()],
            unassessed_transcriptome_assay_ids: vec!["transcriptome-gap".into()],
            unlinked_pair_summaries: vec![citation(0)],
            nonclaims: vec!["Synthetic limitation (verbatim).".into()],
        }
    }

    fn citation(i: usize) -> TranscriptAssaySpecificityMatrixUnlinkedSummary {
        TranscriptAssaySpecificityMatrixUnlinkedSummary {
            summary_id: format!("citation-{i}"),
            pair_binding_sha256: "citation-only-binding".into(),
            same_oligo_current_assay_ids: vec!["assay-00".into()],
            reason: "Same oligos only; no verified panel association.".into(),
        }
    }

    type Labels = Vec<(String, egui::Rect)>;
    fn frame(
        ctx: &egui::Context,
        events: Vec<egui::Event>,
        render: impl FnOnce(&mut egui::Ui),
    ) -> Labels {
        fn collect(shape: &egui::epaint::Shape, labels: &mut Labels) {
            match shape {
                egui::epaint::Shape::Text(t) => labels.push((
                    t.galley.job.text.clone(),
                    egui::Rect::from_min_size(t.pos, t.galley.size()),
                )),
                egui::epaint::Shape::Vec(shapes) => {
                    for s in shapes {
                        collect(s, labels);
                    }
                }
                _ => {}
            }
        }
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(3200., 24000.),
            )),
            events,
            ..Default::default()
        });
        crate::egui_compat::show_central_panel_for_test_context(
            ctx,
            egui::CentralPanel::default(),
            render,
        );
        let mut labels = vec![];
        for s in &crate::egui_compat::end_test_pass(ctx).shapes {
            collect(&s.shape, &mut labels);
        }
        labels
    }

    fn painted(ctx: &egui::Context, matrix: &Matrix) -> Labels {
        frame(ctx, vec![], |ui| {
            render_transcript_assay_specificity_matrix(ui, matrix)
        })
    }
    fn text(labels: &Labels) -> String {
        labels
            .iter()
            .map(|(s, _)| s.as_str())
            .collect::<Vec<_>>()
            .join("\n")
    }
    fn click_at(pos: egui::Pos2) -> Vec<egui::Event> {
        let mut events = vec![egui::Event::PointerMoved(pos)];
        for pressed in [true, false] {
            events.push(egui::Event::PointerButton {
                pos,
                button: egui::PointerButton::Primary,
                pressed,
                modifiers: egui::Modifiers::NONE,
            });
        }
        events
    }
    fn next(ctx: &egui::Context, matrix: &Matrix, section: &str) {
        let labels = painted(ctx, matrix);
        let row = labels
            .iter()
            .find(|(s, _)| s == section)
            .unwrap()
            .1
            .center()
            .y;
        let pos = labels
            .iter()
            .find(|(s, r)| s == "Next" && (r.center().y - row).abs() < 5.)
            .unwrap()
            .1
            .center();
        frame(ctx, click_at(pos), |ui| {
            render_transcript_assay_specificity_matrix(ui, matrix)
        });
        painted(ctx, matrix);
    }

    #[test]
    fn transcript_assay_specificity_review_renders_fields_history_and_citations_without_mutation() {
        let ctx = egui::Context::default();
        ctx.memory_mut(|m| m.set_everything_is_visible(true));
        let mut m = matrix(2, 2, 2);
        m.selections[0].cells[0]
            .diagnostics
            .push("assay_not_in_retained_selection".into());
        m.selections[1].cells.remove(1);
        // Only the history-only IDs are sorted; the other axes follow the engine.
        m.assays.swap(0, 1);
        m.selections.swap(0, 1);
        m.selections.swap(2, 3);
        let before = serde_json::to_vec(&m).unwrap();
        painted(&ctx, &m);
        let all = text(&painted(&ctx, &m));
        assert_eq!(serde_json::to_vec(&m).unwrap(), before);
        let (genomic, rest) = all.split_once("Transcriptome specificity").unwrap();
        let (transcriptome, citations) = rest
            .split_once("Same-oligo citations, not panel evidence")
            .unwrap();
        for (dimension, kind, gap) in [
            (genomic, "genomic_dna", "genomic-gap"),
            (transcriptome, "transcriptome_cdna", "transcriptome-gap"),
        ] {
            for expected in [
                kind,
                gap,
                "same-genome",
                "Assembly: assembly-A",
                "Release: release-B",
                "Receipt order is not chronological",
                "synthetic-algorithm: digest-00-01234567...",
                "Active: yes | Superseded: yes",
                "Replaces: previous-selection",
                "Matches current panel design: no",
                "Matches current panel design: yes",
                "Recorded acceptance: pass",
                "Verdict: incomplete",
                "Verdict: pass",
                "Recorded status: pass",
                "Retained binding verified: no",
                "Retained binding verified: yes",
                "Matches current assay: no",
                "Matches current assay: yes",
                "Unverified binding",
                "Not a current assay match",
                "Selection ID:",
                "Acceptance ID: accept-0",
                "Panel design digest: design-0",
                "Report ID: child-report",
                "Report content SHA-256: synthetic-child-hash",
                "synthetic_receipt_diagnostic",
                "synthetic_cell_diagnostic",
            ] {
                assert!(
                    dimension.contains(expected),
                    "missing {expected}: {dimension}"
                );
            }
            let (current, history) = dimension
                .split_once("Current saved design")
                .unwrap()
                .1
                .split_once("Only in retained receipts (not in the saved design)")
                .unwrap();
            assert!(current.find("assay-01").unwrap() < current.find("assay-00").unwrap());
            assert!(history.find("old-00").unwrap() < history.find("old-01").unwrap());
            assert!(!dimension.contains("citation-only-binding"));
            assert!(
                dimension.find(&format!("{kind}-01")).unwrap()
                    < dimension.find(&format!("{kind}-00")).unwrap()
            );
        }
        assert!(!genomic.contains("transcriptome-gap") && !transcriptome.contains("genomic-gap"));
        assert!(
            genomic.contains("not recorded in this receipt")
                && genomic.contains("not in receipt: unassessed")
        );
        assert!(all.contains(&format!(
            "Reference: {}",
            serde_json::json!(m.selections[0].reference)
        )));
        assert!(all.contains(&format!(
            "Policy: {}",
            serde_json::json!(m.selections[0].policy)
        )));
        for r in &m.selections {
            assert!(all.contains(&r.reference.content_fingerprint));
        }
        assert!(
            all.contains("Current panel design digest: current-digest")
                && all.contains(&format!("Panel report: {PANEL}"))
        );
        assert!(all.contains(&m.nonclaims[0]));
        assert!(
            citations.contains("citation-0")
                && citations.contains("citation-only-binding")
                && citations.contains(&m.unlinked_pair_summaries[0].reason)
        );
        for forbidden in [
            "panel pass",
            "order-ready",
            "0 off-targets",
            "probe evidence",
        ] {
            assert!(!all.contains(forbidden));
        }
    }

    #[test]
    fn transcript_assay_specificity_review_empty_receipts_and_assays_still_show_evidence_gaps() {
        let ctx = egui::Context::default();
        for n in [2, 0] {
            let m = matrix(n, 0, 0);
            painted(&ctx, &m);
            let all = text(&painted(&ctx, &m));
            for expected in [
                "Genomic: no retained receipts",
                "Transcriptome: no retained receipts",
                "genomic-gap",
                "transcriptome-gap",
                "0-0 of 0",
                &m.nonclaims[0],
            ] {
                assert!(all.contains(expected), "missing {expected}");
            }
        }
    }

    #[test]
    fn transcript_assay_specificity_review_pages_both_axes_independently_to_the_last_cell() {
        let ctx = egui::Context::default();
        let m = matrix(61, 13, 7);
        let before = serde_json::to_vec(&m).unwrap();
        painted(&ctx, &m);
        for kind in ["Genomic", "Transcriptome"] {
            for section in ["receipts", "assays"] {
                for _ in 0..2 {
                    next(&ctx, &m, &format!("{kind} {section}"));
                }
            }
            let all = text(&painted(&ctx, &m));
            assert!(all.contains("49-61 of 61") && all.contains("13-13 of 13"));
            if kind == "Genomic" {
                assert!(all.contains("1-24 of 61") && all.contains("1-6 of 13"));
            }
        }
        let all = text(&painted(&ctx, &m));
        for line in [
            "genomic_dna-12",
            "transcriptome_cdna-12",
            "assay-60",
            "old-06",
            "1-7 of 7",
        ] {
            assert!(all.lines().any(|s| s == line), "missing {line}");
        }
        assert_eq!(serde_json::to_vec(&m).unwrap(), before);
    }

    #[test]
    fn transcript_assay_specificity_review_pages_retained_only_rows_and_citations() {
        let ctx = egui::Context::default();
        let mut m = matrix(0, 1, 25);
        m.unlinked_pair_summaries = (0..25).map(citation).collect();
        painted(&ctx, &m);
        next(&ctx, &m, "Genomic retained-only assays");
        next(&ctx, &m, "Same-oligo citations");
        let all = text(&painted(&ctx, &m));
        assert_eq!(all.matches("25-25 of 25").count(), 2);
        assert!(all.contains("old-24") && all.contains("Summary: citation-24"));
    }

    type Shared = Arc<RwLock<GentleEngine>>;
    fn area() -> (MainAreaDna, Shared) {
        let mut engine = GentleEngine::new();
        let primer = |sequence: &str| PrimerDesignPrimerRecord {
            sequence: sequence.into(),
            length_bp: sequence.len(),
            anneal_length_bp: sequence.len(),
            ..Default::default()
        };
        let report = TranscriptAssayPanelReport {
            report_id: PANEL.into(),
            selected_assays: vec![TranscriptAssayPanelAssay {
                assay_id: "saved-assay".into(),
                rank: 1,
                primer_pair: PrimerDesignPairRecord {
                    forward: primer("ACGTACGT"),
                    reverse: primer("TGCATGCA"),
                    ..Default::default()
                },
                ..Default::default()
            }],
            ..Default::default()
        };
        engine.state_mut().metadata.insert(
            PRIMER_DESIGN_REPORTS_METADATA_KEY.into(),
            serde_json::json!({"transcript_assay_panels": {(PANEL): report}}),
        );
        // Both this fixture and a later reopened project start at revision zero.
        let shared = Arc::new(RwLock::new(GentleEngine::from_state(
            engine.state().clone(),
        )));
        let mut area = MainAreaDna::new(
            DNAsequence::from_sequence("ACGT").unwrap(),
            Some("synthetic".into()),
            Some(shared.clone()),
        );
        area.show_transcript_assay_panel_report(RAW_ID);
        assert_eq!(
            area.cached_transcript_assay_panel_report
                .as_ref()
                .unwrap()
                .report_id,
            PANEL
        );
        (area, shared)
    }
    fn snapshot(shared: &Shared) -> (serde_json::Value, usize, u64) {
        let g = shared.read().unwrap();
        (
            serde_json::to_value(g.state()).unwrap(),
            g.journal_len(),
            g.mutation_revision(),
        )
    }
    fn inject(area: &mut MainAreaDna, shared: &Shared) -> mpsc::Sender<Reply> {
        let g = shared.read().unwrap();
        let (tx, rx) = mpsc::channel();
        area.transcript_assay_specificity_review.task = Some(Task {
            binding: Binding {
                project: (g.instance_id(), g.structural_revision()),
                report: area.cached_transcript_assay_panel_report.clone().unwrap(),
            },
            receiver: Arc::new(Mutex::new(Some(rx))),
        });
        tx
    }
    fn finish(area: &mut MainAreaDna, ctx: &egui::Context) {
        let deadline = Instant::now() + Duration::from_secs(10);
        while area.transcript_assay_specificity_review.task.is_some() && Instant::now() < deadline {
            assert!(area.primer_design_task.is_none());
            area.poll_primer_design_task(ctx);
            std::thread::sleep(Duration::from_millis(5));
        }
        assert!(
            area.transcript_assay_specificity_review.task.is_none(),
            "{}",
            area.transcript_assay_specificity_review.status
        );
    }

    #[test]
    fn transcript_assay_specificity_review_worker_uses_displayed_canonical_id_and_preserves_state()
    {
        let (mut area, shared) = area();
        let before = snapshot(&shared);
        let displayed = area.cached_transcript_assay_panel_report.clone().unwrap();
        area.transcript_assay_panel_ui.report_id = "different unsaved text".into();
        let expected = shared
            .read()
            .unwrap()
            .inspect_transcript_assay_panel_specificity_matrix(RAW_ID)
            .unwrap();
        let ctx = egui::Context::default();
        area.start_transcript_assay_specificity_review(&ctx);
        assert!(area.transcript_assay_specificity_review.task.is_some());
        assert!(area.primer_design_task.is_none());
        {
            let _held = shared.write().unwrap();
            for _ in 0..3 {
                area.poll_primer_design_task(&ctx);
                frame(&ctx, vec![], |ui| {
                    area.specificity_review_button(ui);
                    area.show_transcript_assay_specificity_review(ui);
                });
                assert!(area.transcript_assay_specificity_review.task.is_some());
            }
        }
        finish(&mut area, &ctx);
        let (_, got) = area
            .transcript_assay_specificity_review
            .accepted
            .as_ref()
            .unwrap();
        assert_eq!(got.panel_report_id, displayed.report_id);
        assert_eq!(got.unassessed_genomic_assay_ids, ["saved-assay"]);
        assert_eq!(
            serde_json::to_value(got.as_ref()).unwrap(),
            serde_json::to_value(expected).unwrap()
        );
        assert!(Arc::ptr_eq(
            &displayed,
            area.cached_transcript_assay_panel_report.as_ref().unwrap()
        ));
        assert!(area.primer_design_task.is_none());
        assert_eq!(snapshot(&shared), before);
    }

    #[test]
    fn transcript_assay_specificity_review_button_starts_only_one_read() {
        let (mut area, shared) = area();
        let before = snapshot(&shared);
        let ctx = egui::Context::default();
        frame(&ctx, vec![], |ui| area.specificity_review_button(ui));
        let labels = frame(&ctx, vec![], |ui| area.specificity_review_button(ui));
        let pos = labels
            .iter()
            .find(|(s, _)| s == "Review specificity")
            .unwrap()
            .1
            .center();
        frame(&ctx, click_at(pos), |ui| area.specificity_review_button(ui));
        let inbox = area
            .transcript_assay_specificity_review
            .task
            .as_ref()
            .unwrap()
            .receiver
            .clone();
        frame(&ctx, click_at(pos), |ui| area.specificity_review_button(ui));
        assert!(Arc::ptr_eq(
            &inbox,
            &area
                .transcript_assay_specificity_review
                .task
                .as_ref()
                .unwrap()
                .receiver
        ));
        finish(&mut area, &ctx);
        assert!(area.transcript_assay_specificity_review.accepted.is_some());
        assert_eq!(snapshot(&shared), before);
    }

    #[test]
    fn transcript_assay_specificity_review_busy_admission_and_missing_panel_are_explicit() {
        let (mut area, shared) = area();
        let ctx = egui::Context::default();
        {
            let _held = shared.write().unwrap();
            area.start_transcript_assay_specificity_review(&ctx);
            assert!(area.transcript_assay_specificity_review.task.is_none());
            assert_eq!(
                area.transcript_assay_specificity_review.status,
                "Project busy; retry shortly"
            );
        }
        shared
            .write()
            .unwrap()
            .state_mut()
            .metadata
            .remove(PRIMER_DESIGN_REPORTS_METADATA_KEY);
        let expected = shared
            .read()
            .unwrap()
            .inspect_transcript_assay_panel_specificity_matrix(PANEL)
            .unwrap_err()
            .message;
        area.start_transcript_assay_specificity_review(&ctx);
        finish(&mut area, &ctx);
        assert!(area.transcript_assay_specificity_review.accepted.is_none());
        assert_eq!(area.transcript_assay_specificity_review.status, expected);
    }

    #[test]
    fn transcript_assay_specificity_review_rejects_changed_arc_revision_owner_or_canonical_id() {
        let ctx = egui::Context::default();
        for case in ["displayed panel", "revision", "owner", "report ID"] {
            let (mut area, shared) = area();
            let mut m = matrix(0, 0, 0);
            let tx = inject(&mut area, &shared);
            match case {
                "displayed panel" => {
                    let same_id =
                        (**area.cached_transcript_assay_panel_report.as_ref().unwrap()).clone();
                    area.cached_transcript_assay_panel_report = Some(Arc::new(same_id));
                }
                "revision" => {
                    let mut g = shared.write().unwrap();
                    g.state_mut()
                        .metadata
                        .get_mut(PRIMER_DESIGN_REPORTS_METADATA_KEY)
                        .unwrap()["transcript_assay_panels"][PANEL]["source_seq_id"] =
                        "redesigned".into();
                }
                "owner" => {
                    let mut g = shared.write().unwrap();
                    let replacement = GentleEngine::from_state(g.state().clone());
                    assert_eq!(g.structural_revision(), replacement.structural_revision());
                    assert_ne!(g.instance_id(), replacement.instance_id());
                    *g = replacement;
                }
                _ => m.panel_report_id = "wrong-panel".into(),
            }
            tx.send(Ok(Arc::new(m))).unwrap();
            area.poll_transcript_assay_specificity_review(&ctx);
            let r = &area.transcript_assay_specificity_review;
            assert!(r.accepted.is_none() && r.task.is_none());
            assert!(
                r.status.contains("dropped")
                    && r.status.contains(match case {
                        "owner" | "revision" => "project changed",
                        other => other,
                    }),
                "{}",
                r.status
            );
        }
    }

    #[test]
    fn transcript_assay_specificity_review_queued_reply_survives_contention_then_invalidates_on_paint()
     {
        let (mut area, shared) = area();
        let ctx = egui::Context::default();
        let tx = inject(&mut area, &shared);
        tx.send(Ok(Arc::new(matrix(0, 0, 0)))).unwrap();
        {
            let _held = shared.write().unwrap();
            area.poll_transcript_assay_specificity_review(&ctx);
            assert!(area.transcript_assay_specificity_review.task.is_some());
            assert!(area.transcript_assay_specificity_review.accepted.is_none());
        }
        area.poll_transcript_assay_specificity_review(&ctx);
        {
            let mut held = shared.write().unwrap();
            held.state_mut()
                .metadata
                .insert("changed".into(), true.into());
            frame(&ctx, vec![], |ui| {
                area.show_transcript_assay_specificity_review(ui)
            });
            assert!(area.transcript_assay_specificity_review.accepted.is_some());
        }
        let labels = frame(&ctx, vec![], |ui| {
            area.show_transcript_assay_specificity_review(ui)
        });
        assert!(area.transcript_assay_specificity_review.accepted.is_none());
        assert!(text(&labels).contains("Project changed since review; review again"));
    }

    #[test]
    fn transcript_assay_specificity_review_reload_and_closed_views_disconnect_late_replies() {
        let (mut area, shared) = area();
        let before = snapshot(&shared);
        let tx = inject(&mut area, &shared);
        area.show_transcript_assay_panel_report(RAW_ID);
        assert!(tx.send(Ok(Arc::new(matrix(0, 0, 0)))).is_err());
        let ctx = egui::Context::default();
        area.set_rendering_pcr_specialist(true);
        area.start_transcript_assay_specificity_review(&ctx);
        MainAreaDna::disconnect_hosted_specificity_review(&ctx);
        area.poll_transcript_assay_specificity_review(&ctx);
        assert!(area.transcript_assay_specificity_review.accepted.is_none());
        assert!(
            area.transcript_assay_specificity_review
                .status
                .contains("closed")
        );
        let tx = inject(&mut area, &shared);
        drop(area);
        assert!(tx.send(Ok(Arc::new(matrix(0, 0, 0)))).is_err());
        assert_eq!(snapshot(&shared), before);
    }
}
