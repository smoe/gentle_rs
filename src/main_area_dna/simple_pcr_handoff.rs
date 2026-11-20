//! Explicit, session-only copy of a validated saved pair into the simple PCR form.
//!
//! Preparation and overwrite confirmation never apply operations or persist UI
//! metadata. PCR execution remains the separate existing PCR button.

use super::*;
use crate::engine::SimplePcrPairPreparation;

#[derive(Default)]
pub(super) struct SimplePcrHandoffUiState {
    pending: Option<PendingSimplePcrHandoff>,
    prepared: Option<SimplePcrPairPreparation>,
    pub(super) focus_requested: bool,
    pub(super) owner_viewport: Option<egui::ViewportId>,
}

struct PendingSimplePcrHandoff {
    report: Arc<PrimerDesignReport>,
    prepared: SimplePcrPairPreparation,
    engine_instance: u64,
    previous_forward: String,
    previous_reverse: String,
}

impl MainAreaDna {
    /// Returns true when copied, false when explicit overwrite approval is needed.
    pub(super) fn request_simple_pcr_handoff(
        &mut self,
        report: Arc<PrimerDesignReport>,
        pair_rank: usize,
    ) -> Result<bool, String> {
        self.simple_pcr_handoff.pending = None;
        if !self.dna_presentation_mode.allows_engine_shell_panels() {
            return Err("Choose Region, Gene or cDNA view to access the PCR form".to_string());
        }
        let template = self.seq_id.as_deref().ok_or("No active template")?;
        let engine = self.engine.as_ref().ok_or("No engine attached")?;
        let (prepared, engine_instance) = {
            let engine = engine.read().map_err(|_| "Engine lock poisoned")?;
            (
                engine
                    .prepare_simple_pcr_pair(&report, pair_rank, template)
                    .map_err(|error| error.message)?,
                engine.instance_id(),
            )
        };
        let pending = PendingSimplePcrHandoff {
            report,
            prepared,
            engine_instance,
            previous_forward: self.pcr_forward.clone(),
            previous_reverse: self.pcr_reverse.clone(),
        };
        let needs_confirmation = (!self.pcr_forward.is_empty() || !self.pcr_reverse.is_empty())
            && (self.pcr_forward != pending.prepared.forward_primer
                || self.pcr_reverse != pending.prepared.reverse_primer);
        if needs_confirmation {
            self.op_status =
                "Review and confirm replacement of the existing PCR primers".to_string();
            self.simple_pcr_handoff.pending = Some(pending);
            return Ok(false);
        }
        self.apply_simple_pcr_handoff(pending)?;
        Ok(true)
    }

    pub(super) fn confirm_simple_pcr_handoff(&mut self) -> Result<(), String> {
        let pending = self
            .simple_pcr_handoff
            .pending
            .take()
            .ok_or("No pending PCR primer replacement")?;
        self.apply_simple_pcr_handoff(pending)
    }

    pub(super) fn cancel_simple_pcr_handoff(&mut self) {
        self.simple_pcr_handoff.pending = None;
        self.op_status = "PCR primer replacement cancelled; input fields unchanged".to_string();
    }

    pub(super) fn focus_prepared_simple_pcr_form(&self, ctx: &egui::Context) {
        if let Some(viewport) = self.simple_pcr_handoff.owner_viewport {
            ctx.send_viewport_cmd_to(viewport, egui::ViewportCommand::Visible(true));
            ctx.send_viewport_cmd_to(viewport, egui::ViewportCommand::Focus);
            ctx.request_repaint_of(viewport);
        }
        ctx.request_repaint_of(egui::ViewportId::ROOT);
    }

    fn apply_simple_pcr_handoff(&mut self, pending: PendingSimplePcrHandoff) -> Result<(), String> {
        if !self.dna_presentation_mode.allows_engine_shell_panels() {
            return Err("Sequence view changed; request the handoff again".to_string());
        }
        if self.pcr_forward != pending.previous_forward
            || self.pcr_reverse != pending.previous_reverse
        {
            return Err("PCR fields changed; request the handoff again".to_string());
        }
        let template = self.seq_id.as_deref().ok_or("No active template")?;
        let engine = self.engine.as_ref().ok_or("No engine attached")?;
        {
            let engine = engine.read().map_err(|_| "Engine lock poisoned")?;
            if engine.instance_id() != pending.engine_instance {
                return Err("Project changed; request the handoff again".to_string());
            }
            let current = engine
                .prepare_simple_pcr_pair(&pending.report, pending.prepared.pair_rank, template)
                .map_err(|error| error.message)?;
            if current != pending.prepared {
                return Err("PCR source changed; request the handoff again".to_string());
            }
        }
        self.pcr_forward = pending.prepared.forward_primer.clone();
        self.pcr_reverse = pending.prepared.reverse_primer.clone();
        self.show_engine_ops = true;
        self.simple_pcr_handoff.focus_requested = true;
        self.op_status = format!(
            "Prepared simple PCR from report '{}' pair {} on '{}' (backend {} -> {}). Review the form and click PCR separately; no specificity/order approval.",
            pending.prepared.report_id,
            pending.prepared.pair_rank,
            pending.prepared.template,
            pending.prepared.backend_requested,
            pending.prepared.backend_used,
        );
        self.simple_pcr_handoff.prepared = Some(pending.prepared);
        self.op_error_popup = None;
        Ok(())
    }

    pub(super) fn render_simple_pcr_handoff_confirmation(&mut self, ui: &mut egui::Ui) {
        let Some(pending) = self.simple_pcr_handoff.pending.as_ref() else {
            return;
        };
        let mut confirm = false;
        let mut cancel = false;
        ui.group(|ui| {
            ui.strong("Replace existing simple-PCR primers?");
            ui.label(format!(
                "Report '{}', pair {}, template '{}'",
                pending.prepared.report_id, pending.prepared.pair_rank, pending.prepared.template
            ));
            ui.monospace(format!(
                "Forward: {} -> {}",
                pending.previous_forward, pending.prepared.forward_primer
            ));
            ui.monospace(format!(
                "Reverse: {} -> {}",
                pending.previous_reverse, pending.prepared.reverse_primer
            ));
            ui.small("This only fills the form. It does not run PCR or approve the assay.");
            ui.horizontal(|ui| {
                confirm = ui.button("Replace PCR primers").clicked();
                cancel = ui.button("Cancel primer replacement").clicked();
            });
        });
        if cancel {
            self.cancel_simple_pcr_handoff();
        } else if confirm {
            match self.confirm_simple_pcr_handoff() {
                Ok(()) => self.focus_prepared_simple_pcr_form(ui.ctx()),
                Err(message) => {
                    self.op_status = message.clone();
                    self.op_error_popup = Some(message);
                }
            }
        }
    }

    pub(super) fn render_simple_pcr_handoff_source(&self, ui: &mut egui::Ui) {
        let Some(source) = self.simple_pcr_handoff.prepared.as_ref() else {
            return;
        };
        if self.pcr_forward != source.forward_primer || self.pcr_reverse != source.reverse_primer {
            ui.small("PCR primers were edited after preparation; the saved-pair binding no longer applies.");
            return;
        }
        ui.small(format!(
            "Prepared from '{}' pair {}. Binding checked at preparation, not renewed specificity approval.",
            source.report_id, source.pair_rank,
        )).on_hover_text(format!("Report: {}\nTemplate at preparation: {}",
            source.report_sha256, source.template_sha256));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::simple_pcr_handoff::tests::{fixture, install};

    fn area() -> (MainAreaDna, Arc<PrimerDesignReport>) {
        let (engine, report) = fixture();
        let dna = engine.state().sequences[&report.template].clone();
        let engine = Arc::new(RwLock::new(engine));
        let area = MainAreaDna::new(dna, Some(report.template.clone()), Some(engine));
        (area, Arc::new(report))
    }

    fn snapshot(area: &MainAreaDna) -> serde_json::Value {
        let engine = area.engine.as_ref().unwrap().read().unwrap();
        serde_json::json!({"state": engine.state(), "history": engine.history_summary()})
    }

    #[test]
    fn simple_pcr_handoff_gui_prepares_only_basic_fields_without_project_mutation() {
        let (mut area, report) = area();
        area.pcr_adv_forward = "advanced unchanged".to_string();
        area.pcr_adv_reverse = "advanced reverse unchanged".to_string();
        area.pcr_adv_anneal_len = "19".to_string();
        area.pcr_adv_max_mismatch = "2".to_string();
        area.pcr_adv_3prime_exact = "5".to_string();
        area.pcr_mut_position = "11".to_string();
        area.pcr_mut_ref = "C".to_string();
        area.pcr_mut_alt = "T".to_string();
        area.pcr_unique = true;
        let other = (
            area.pcr_adv_forward.clone(),
            area.pcr_adv_reverse.clone(),
            area.pcr_adv_anneal_len.clone(),
            area.pcr_adv_max_mismatch.clone(),
            area.pcr_adv_3prime_exact.clone(),
            area.pcr_mut_position.clone(),
            area.pcr_mut_ref.clone(),
            area.pcr_mut_alt.clone(),
            area.pcr_unique,
        );
        let before = snapshot(&area);
        assert!(area.request_simple_pcr_handoff(report, 7).unwrap());
        assert_eq!(area.pcr_forward, "AACCG");
        assert_eq!(area.pcr_reverse, "TAGTGC");
        assert!(area.show_engine_ops && area.simple_pcr_handoff.focus_requested);
        assert!(area.op_status.contains("pair 7"));
        assert_eq!(
            other,
            (
                area.pcr_adv_forward.clone(),
                area.pcr_adv_reverse.clone(),
                area.pcr_adv_anneal_len.clone(),
                area.pcr_adv_max_mismatch.clone(),
                area.pcr_adv_3prime_exact.clone(),
                area.pcr_mut_position.clone(),
                area.pcr_mut_ref.clone(),
                area.pcr_mut_alt.clone(),
                area.pcr_unique
            )
        );
        assert_eq!(before, snapshot(&area));
    }

    #[test]
    fn simple_pcr_handoff_gui_overwrite_requires_confirmation_and_cancel_preserves_inputs() {
        let (mut area, report) = area();
        area.pcr_forward = "old forward".to_string();
        area.pcr_reverse = "old reverse".to_string();
        let before = snapshot(&area);
        assert!(
            !area
                .request_simple_pcr_handoff(Arc::clone(&report), 7)
                .unwrap()
        );
        assert_eq!(area.pcr_forward, "old forward");
        area.cancel_simple_pcr_handoff();
        assert_eq!(area.pcr_reverse, "old reverse");
        assert!(area.confirm_simple_pcr_handoff().is_err());
        assert!(!area.request_simple_pcr_handoff(report, 7).unwrap());
        area.confirm_simple_pcr_handoff().unwrap();
        assert_eq!(area.pcr_reverse, "TAGTGC");
        assert_eq!(before, snapshot(&area));
    }

    #[test]
    fn simple_pcr_handoff_gui_confirmation_refuses_edited_destination() {
        let (mut area, report) = area();
        area.pcr_forward = "old".to_string();
        assert!(!area.request_simple_pcr_handoff(report, 7).unwrap());
        area.pcr_forward = "new manual edit".to_string();
        assert!(area.confirm_simple_pcr_handoff().is_err());
        assert_eq!(area.pcr_forward, "new manual edit");
        assert!(area.pcr_reverse.is_empty());
        assert!(area.simple_pcr_handoff.pending.is_none());
    }

    #[test]
    fn simple_pcr_handoff_gui_confirmation_refuses_changed_template_or_report() {
        for change_report in [false, true] {
            let (mut area, report) = area();
            area.pcr_forward = "old".to_string();
            assert!(
                !area
                    .request_simple_pcr_handoff(Arc::clone(&report), 7)
                    .unwrap()
            );
            {
                let mut engine = area.engine.as_ref().unwrap().write().unwrap();
                if change_report {
                    let mut changed = report.as_ref().clone();
                    changed.backend.used = "changed".to_string();
                    install(&mut engine, &changed);
                } else {
                    // Change only the interior, retaining both valid primer sites.
                    engine.state_mut().sequences.insert(
                        report.template.clone(),
                        DNAsequence::from_sequence("AACCGGCTTTGCATGCACTA").unwrap(),
                    );
                }
            }
            let before = snapshot(&area);
            assert!(area.confirm_simple_pcr_handoff().is_err());
            assert_eq!(area.pcr_forward, "old");
            assert!(area.pcr_reverse.is_empty());
            assert_eq!(before, snapshot(&area));
        }
    }

    #[test]
    fn simple_pcr_handoff_gui_refuses_stale_preview_bad_rank_tails_and_wrong_template() {
        let (mut area, report) = area();
        for rank in [0, 1, usize::MAX] {
            assert!(
                area.request_simple_pcr_handoff(Arc::clone(&report), rank)
                    .is_err()
            );
        }
        area.seq_id = Some("other".to_string());
        assert!(
            area.request_simple_pcr_handoff(Arc::clone(&report), 7)
                .is_err()
        );
        area.seq_id = Some(report.template.clone());
        area.dna_presentation_mode = DnaPresentationMode::Genome;
        assert!(
            area.request_simple_pcr_handoff(Arc::clone(&report), 7)
                .is_err()
        );
        area.dna_presentation_mode = DnaPresentationMode::Region;
        let mut changed = report.as_ref().clone();
        changed.pairs[0].forward.non_annealing_5prime_tail_bp = 1;
        install(
            &mut area.engine.as_ref().unwrap().write().unwrap(),
            &changed,
        );
        assert!(area.request_simple_pcr_handoff(report, 7).is_err());
        assert!(
            area.request_simple_pcr_handoff(Arc::new(changed), 7)
                .is_err()
        );
        assert!(area.pcr_forward.is_empty() && area.pcr_reverse.is_empty());
    }

    #[test]
    fn simple_pcr_handoff_gui_confirmation_cannot_cross_project_instances() {
        let (mut area, report) = area();
        area.pcr_forward = "old".to_string();
        assert!(!area.request_simple_pcr_handoff(report, 7).unwrap());
        area.engine = Some(Arc::new(RwLock::new(fixture().0)));
        let before = snapshot(&area);
        assert!(area.confirm_simple_pcr_handoff().is_err());
        assert_eq!(area.pcr_forward, "old");
        assert!(area.pcr_reverse.is_empty());
        assert_eq!(before, snapshot(&area));
    }

    #[test]
    fn simple_pcr_handoff_gui_confirmation_cannot_bypass_read_only_sequence_view() {
        let (mut area, report) = area();
        area.pcr_forward = "old".to_string();
        assert!(!area.request_simple_pcr_handoff(report, 7).unwrap());
        area.dna_presentation_mode = DnaPresentationMode::Chromosomal;
        let before = snapshot(&area);
        assert!(area.confirm_simple_pcr_handoff().is_err());
        assert_eq!(area.pcr_forward, "old");
        assert!(area.pcr_reverse.is_empty());
        assert_eq!(before, snapshot(&area));
    }

    #[test]
    fn simple_pcr_handoff_gui_pending_state_is_not_persisted_in_project_ui_metadata() {
        let (mut area, report) = area();
        area.pcr_forward = "old".to_string();
        assert!(!area.request_simple_pcr_handoff(report, 7).unwrap());
        let state = area.current_engine_ops_state();
        assert!(
            !serde_json::to_string(&state)
                .unwrap()
                .contains("simple_pcr_handoff")
        );
        area.apply_engine_ops_state(state);
        assert!(area.simple_pcr_handoff.pending.is_none());
    }

    #[test]
    fn simple_pcr_handoff_gui_idempotent_preparation_does_not_prompt_for_identical_fields() {
        let (mut area, report) = area();
        assert!(
            area.request_simple_pcr_handoff(Arc::clone(&report), 7)
                .unwrap()
        );
        let before = snapshot(&area);
        assert!(area.request_simple_pcr_handoff(report, 7).unwrap());
        assert!(area.simple_pcr_handoff.pending.is_none());
        assert_eq!(before, snapshot(&area));
    }

    #[test]
    fn simple_pcr_handoff_gui_report_preview_paints_without_copying_or_running_pcr() {
        let (mut area, report) = area();
        area.primer_design_ui.report_id = report.report_id.clone();
        area.primer_design_ui.specificity_pair_rank_1based = "7".to_string();
        let before = snapshot(&area);
        let ctx = egui::Context::default();
        let output = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                area.render_primer_design_report_preview(ui);
            });
        });
        assert!(area.pcr_forward.is_empty() && area.pcr_reverse.is_empty());
        assert!(area.simple_pcr_handoff.pending.is_none());
        assert!(output.platform_output.commands.is_empty());
        assert_eq!(before, snapshot(&area));
    }
}
