//! Selection-formula and ROI-formula helpers for `MainAreaDna`.
//!
//! This submodule keeps the sequence-window formula entry path in one place:
//! shared numeric parsers, feature-relative ROI resolution, and the small UI
//! controls that apply those formulas back into the window state.

use super::*;

impl MainAreaDna {
    pub(super) fn parse_optional_usize_text(
        raw: &str,
        field_name: &str,
    ) -> Result<Option<usize>, String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            Ok(None)
        } else {
            trimmed
                .parse::<usize>()
                .map(Some)
                .map_err(|_| format!("Invalid {field_name}: expected an integer"))
        }
    }

    pub(super) fn parse_positive_usize_text(raw: &str, field_name: &str) -> Result<usize, String> {
        let trimmed = raw.trim();
        let parsed = trimmed
            .parse::<usize>()
            .map_err(|_| format!("Invalid {field_name}: expected an integer"))?;
        if parsed == 0 {
            return Err(format!("Invalid {field_name}: expected >= 1"));
        }
        Ok(parsed)
    }

    pub(super) fn parse_optional_f64_text(
        raw: &str,
        field_name: &str,
    ) -> Result<Option<f64>, String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            Ok(None)
        } else {
            trimmed
                .parse::<f64>()
                .map(Some)
                .map_err(|_| format!("Invalid {field_name}: expected a number"))
        }
    }

    pub(super) fn resolve_roi_range_inputs_0based(
        &self,
        roi_start_raw: &str,
        roi_end_raw: &str,
        field_prefix: &str,
    ) -> Result<(usize, usize), String> {
        let dna = self
            .dna
            .read()
            .map_err(|_| format!("Could not read DNA while parsing {field_prefix}"))?;
        resolve_formula_roi_range_inputs_0based_on_sequence(
            &dna,
            roi_start_raw,
            roi_end_raw,
            field_prefix,
        )
    }

    pub(super) fn resolve_selection_formula_range_0based(
        &self,
        raw: &str,
    ) -> Result<(usize, usize), String> {
        let dna = self.dna.read().map_err(|_| {
            "Cannot apply selection formula: could not read active sequence".to_string()
        })?;
        resolve_selection_formula_range_0based_on_sequence(&dna, raw)
    }

    pub(super) fn apply_selection_formula(&mut self) {
        let (start, end_exclusive) =
            match self.resolve_selection_formula_range_0based(&self.selection_formula_text) {
                Ok(value) => value,
                Err(err) => {
                    self.op_status = err;
                    return;
                }
            };
        match self.set_selection_range_0based(start, end_exclusive) {
            Ok(()) => {
                self.op_status = format!(
                    "Selected region from formula: {}..{} (0-based, end-exclusive)",
                    start, end_exclusive
                );
            }
            Err(err) => {
                self.op_status = err;
            }
        }
    }

    pub(super) fn render_selection_formula_inline_controls(
        &mut self,
        ui: &mut egui::Ui,
        desired_width: f32,
    ) {
        ui.vertical(|ui| {
            ui.label("Selection formula")
                .on_hover_text("Set map/text selection from feature-relative formula range");
            ui.horizontal(|ui| {
                let field_width = (ui.available_width() - 96.0).clamp(72.0, desired_width);
                let response = ui
                    .add(
                        egui::TextEdit::singleline(&mut self.selection_formula_text)
                            .desired_width(field_width)
                            .hint_text("=CDS.start+10 .. CDS.end-500"),
                    )
                    .on_hover_text(
                        "Excel-like range formula. Examples: `=CDS.start+10 .. CDS.end-500`, `=gene[label=TP73].start to gene[label=TP73].end`",
                    );
                if response.changed() {
                    self.save_engine_ops_state();
                }
                let apply_clicked = ui
                    .button("Apply Sel")
                    .on_hover_text(
                        "Resolve selection formula and apply it as current map/text selection",
                    )
                    .clicked();
                let enter_pressed = ui.input(|i| i.key_pressed(egui::Key::Enter));
                if apply_clicked || (enter_pressed && response.lost_focus()) {
                    self.apply_selection_formula();
                }
            });
        });
    }

    pub(super) fn resolve_and_apply_primer_roi_fields(&mut self) {
        match self.resolve_roi_range_inputs_0based(
            &self.primer_design_ui.roi_start_0based,
            &self.primer_design_ui.roi_end_0based,
            "primer_design",
        ) {
            Ok((start, end_exclusive)) => {
                self.primer_design_ui.roi_start_0based = start.to_string();
                self.primer_design_ui.roi_end_0based = end_exclusive.to_string();
                self.op_status = format!(
                    "Resolved primer ROI formula to {}..{} (0-based, end-exclusive)",
                    start, end_exclusive
                );
            }
            Err(err) => self.op_status = err,
        }
    }

    pub(super) fn resolve_and_apply_qpcr_roi_fields(&mut self) {
        match self.resolve_roi_range_inputs_0based(
            &self.qpcr_design_ui.roi_start_0based,
            &self.qpcr_design_ui.roi_end_0based,
            "qpcr_design",
        ) {
            Ok((start, end_exclusive)) => {
                self.qpcr_design_ui.roi_start_0based = start.to_string();
                self.qpcr_design_ui.roi_end_0based = end_exclusive.to_string();
                self.op_status = format!(
                    "Resolved qPCR ROI formula to {}..{} (0-based, end-exclusive)",
                    start, end_exclusive
                );
            }
            Err(err) => self.op_status = err,
        }
    }

    pub(super) fn parse_required_usize_text(raw: &str, field_name: &str) -> Result<usize, String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(format!("Invalid {field_name}: expected an integer"));
        }
        trimmed
            .parse::<usize>()
            .map_err(|_| format!("Invalid {field_name}: expected an integer"))
    }
}
