//! Thin GUI adapter for previewing and applying exact feature-boundary edits.

use super::*;
use gentle_protocol::{
    FeatureLocationEditReport, FeatureLocationEditRequest, RelatedFeatureBoundaryReason,
};

#[derive(Clone, Debug, Default)]
pub(super) struct FeatureLocationEditorUiState {
    pub(super) selected_feature_index: Option<usize>,
    pub(super) start_1based: String,
    pub(super) end_1based_inclusive: String,
    pub(super) preview: Option<FeatureLocationEditReport>,
    pub(super) status: String,
}

impl MainAreaDna {
    pub(super) fn feature_supports_location_edit(&self, feature_index: usize) -> bool {
        self.dna
            .read()
            .ok()
            .and_then(|dna| dna.features().get(feature_index).cloned())
            .is_some_and(|feature| {
                crate::feature_location::exact_feature_location_snapshot(&feature).is_ok()
            })
    }

    fn feature_location_editor_options(&self) -> Vec<(usize, String, bool)> {
        let Ok(dna) = self.dna.read() else {
            return vec![];
        };
        dna.features()
            .iter()
            .enumerate()
            .map(|(index, feature)| {
                let label = feature
                    .qualifier_values("label")
                    .next()
                    .or_else(|| feature.qualifier_values("gene").next())
                    .or_else(|| feature.qualifier_values("transcript_id").next())
                    .map(str::to_string)
                    .unwrap_or_else(|| feature.location.to_string());
                (
                    index,
                    format!("{index}: {} ({label})", feature.kind),
                    crate::feature_location::exact_feature_location_snapshot(feature).is_ok(),
                )
            })
            .collect()
    }

    fn seed_feature_location_editor(&mut self, feature_index: usize) {
        let snapshot = self
            .dna
            .read()
            .ok()
            .and_then(|dna| dna.features().get(feature_index).cloned())
            .and_then(|feature| {
                crate::feature_location::exact_feature_location_snapshot(&feature).ok()
            });
        self.feature_location_editor_ui.selected_feature_index = Some(feature_index);
        self.feature_location_editor_ui.preview = None;
        match snapshot {
            Some(snapshot) => {
                self.feature_location_editor_ui.start_1based = snapshot.start_1based.to_string();
                self.feature_location_editor_ui.end_1based_inclusive =
                    snapshot.end_1based_inclusive.to_string();
                self.feature_location_editor_ui.status.clear();
            }
            None => {
                self.feature_location_editor_ui.start_1based.clear();
                self.feature_location_editor_ui.end_1based_inclusive.clear();
                self.feature_location_editor_ui.status =
                    "This feature has a compound or fuzzy location and cannot be edited here."
                        .to_string();
            }
        }
    }

    pub fn focus_feature_location_editor(&mut self, feature_index: Option<usize>) {
        let options = self.feature_location_editor_options();
        let selected = feature_index
            .filter(|index| options.iter().any(|(candidate, _, _)| candidate == index))
            .or(self.feature_location_editor_ui.selected_feature_index)
            .or_else(|| {
                self.get_selected_feature_id()
                    .filter(|index| options.iter().any(|(candidate, _, _)| candidate == index))
            })
            .or_else(|| {
                options
                    .iter()
                    .find_map(|(index, _, supported)| supported.then_some(*index))
            });
        if let Some(selected) = selected {
            self.seed_feature_location_editor(selected);
        } else {
            self.feature_location_editor_ui.status =
                "No editable exact feature is available on this sequence.".to_string();
        }
        self.show_feature_location_editor = true;
    }

    pub fn close_feature_location_editor(&mut self) -> bool {
        let was_open = self.show_feature_location_editor;
        self.show_feature_location_editor = false;
        was_open
    }

    pub fn feature_location_editor_is_open(&self) -> bool {
        self.show_feature_location_editor
    }

    fn feature_location_editor_request(
        &self,
        expected_feature_fingerprint_sha256: Option<String>,
    ) -> Result<FeatureLocationEditRequest, String> {
        let seq_id = self
            .seq_id
            .clone()
            .ok_or_else(|| "No active sequence is attached to this window".to_string())?;
        let feature_index = self
            .feature_location_editor_ui
            .selected_feature_index
            .ok_or_else(|| "Select a feature first".to_string())?;
        let start_1based = self
            .feature_location_editor_ui
            .start_1based
            .trim()
            .parse::<i64>()
            .map_err(|_| "Start must be a positive 1-based coordinate".to_string())?;
        let end_1based_inclusive = self
            .feature_location_editor_ui
            .end_1based_inclusive
            .trim()
            .parse::<i64>()
            .map_err(|_| "End must be a positive 1-based inclusive coordinate".to_string())?;
        if start_1based < 1 || end_1based_inclusive < 1 {
            return Err("Start and end must both be at least 1".to_string());
        }
        let new_start_0based = start_1based
            .checked_sub(1)
            .ok_or_else(|| "Start must be at least 1".to_string())?;
        if end_1based_inclusive < start_1based {
            return Err("End must be greater than or equal to start".to_string());
        }
        Ok(FeatureLocationEditRequest {
            seq_id,
            feature_index,
            new_start_0based,
            new_end_0based_exclusive: end_1based_inclusive,
            expected_feature_fingerprint_sha256,
        })
    }

    pub(super) fn run_feature_location_preview(&mut self) {
        let request = match self.feature_location_editor_request(None) {
            Ok(request) => request,
            Err(err) => {
                self.feature_location_editor_ui.status = err;
                self.feature_location_editor_ui.preview = None;
                return;
            }
        };
        let Some(result) =
            self.apply_operation_with_feedback_and_result(Operation::PreviewFeatureLocationEdit {
                request,
            })
        else {
            self.feature_location_editor_ui.status = self.op_status.clone();
            self.feature_location_editor_ui.preview = None;
            return;
        };
        self.feature_location_editor_ui.preview = result
            .feature_location_edit_report
            .map(|report| report.as_ref().clone());
        self.feature_location_editor_ui.status =
            "Preview is current. Review related annotations before applying.".to_string();
    }

    fn run_feature_location_apply(&mut self) {
        let Some(preview) = self.feature_location_editor_ui.preview.clone() else {
            self.feature_location_editor_ui.status =
                "Run Preview before applying the edit.".to_string();
            return;
        };
        let request = match self.feature_location_editor_request(Some(
            preview.before_feature_fingerprint_sha256.clone(),
        )) {
            Ok(request) => request,
            Err(err) => {
                self.feature_location_editor_ui.status = err;
                return;
            }
        };
        if request.feature_index != preview.feature_index
            || i64::try_from(preview.after.start_0based).ok() != Some(request.new_start_0based)
            || i64::try_from(preview.after.end_0based_exclusive).ok()
                != Some(request.new_end_0based_exclusive)
        {
            self.feature_location_editor_ui.status =
                "Coordinates changed after preview; preview again before applying.".to_string();
            self.feature_location_editor_ui.preview = None;
            return;
        }
        let Some(result) = self
            .apply_operation_with_feedback_and_result(Operation::EditFeatureLocation { request })
        else {
            self.feature_location_editor_ui.status = self.op_status.clone();
            return;
        };
        let feature_index = preview.feature_index;
        let status = result
            .messages
            .first()
            .cloned()
            .unwrap_or_else(|| "Feature location updated.".to_string());
        self.seed_feature_location_editor(feature_index);
        self.feature_location_editor_ui.status = status;
    }

    fn related_reason_label(reason: RelatedFeatureBoundaryReason) -> &'static str {
        match reason {
            RelatedFeatureBoundaryReason::SharesOldStartBoundary => "shares old start",
            RelatedFeatureBoundaryReason::SharesOldEndBoundary => "shares old end",
        }
    }

    pub(super) fn render_feature_location_editor(&mut self, ctx: &egui::Context) {
        if !self.show_feature_location_editor {
            return;
        }
        let options = self.feature_location_editor_options();
        let mut open = self.show_feature_location_editor;
        let spec = crate::egui_compat::HostedWindowSpec::new(
            "Feature Location Editor",
            egui::Id::new(("feature_location_editor", self.panel_scope_key())),
            egui::vec2(620.0, 520.0),
            egui::vec2(460.0, 340.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            ui.label(
                "Edit an exact Range or Complement(Range). Compound and fuzzy GenBank locations require a dedicated curation workflow.",
            );
            ui.separator();
            let selected_label = self
                .feature_location_editor_ui
                .selected_feature_index
                .and_then(|selected| {
                    options
                        .iter()
                        .find(|(index, _, _)| *index == selected)
                        .map(|(_, label, _)| label.clone())
                })
                .unwrap_or_else(|| "Select feature".to_string());
            let mut next_selected = None;
            egui::ComboBox::from_id_salt((
                "feature_location_editor_feature",
                self.panel_scope_key(),
            ))
            .selected_text(selected_label)
            .width(430.0)
            .show_ui(ui, |ui| {
                for (index, label, supported) in &options {
                    let response = ui.add_enabled(
                        *supported,
                        egui::Button::new(label).selected(
                            self.feature_location_editor_ui.selected_feature_index == Some(*index),
                        ),
                    );
                    if response.clicked() {
                        next_selected = Some(*index);
                    }
                }
            });
            if let Some(index) = next_selected {
                self.seed_feature_location_editor(index);
            }
            ui.horizontal(|ui| {
                ui.label("Start (1-based)");
                if ui
                    .text_edit_singleline(&mut self.feature_location_editor_ui.start_1based)
                    .changed()
                {
                    self.feature_location_editor_ui.preview = None;
                }
                ui.label("End (1-based, inclusive)");
                if ui
                    .text_edit_singleline(&mut self.feature_location_editor_ui.end_1based_inclusive)
                    .changed()
                {
                    self.feature_location_editor_ui.preview = None;
                }
            });
            ui.horizontal(|ui| {
                if ui.button("Preview").clicked() {
                    self.run_feature_location_preview();
                }
                let apply = ui.add_enabled(
                    self.feature_location_editor_ui.preview.is_some(),
                    egui::Button::new("Apply"),
                );
                if apply
                    .on_hover_text("Apply is enabled only for the exact coordinates just previewed")
                    .clicked()
                {
                    self.run_feature_location_apply();
                }
            });
            if !self.feature_location_editor_ui.status.is_empty() {
                ui.label(&self.feature_location_editor_ui.status);
            }
            if let Some(preview) = self.feature_location_editor_ui.preview.as_ref() {
                ui.separator();
                ui.label(egui::RichText::new("Preview").strong());
                ui.monospace(format!(
                    "{} {}..{} -> {}..{} ({:?})",
                    preview.feature_kind,
                    preview.before.start_1based,
                    preview.before.end_1based_inclusive,
                    preview.after.start_1based,
                    preview.after.end_1based_inclusive,
                    preview.after.strand
                ));
                ui.small(format!(
                    "5'={}  3'={}  fingerprint={}",
                    preview.after.five_prime_position_1based,
                    preview.after.three_prime_position_1based,
                    preview.before_feature_fingerprint_sha256
                ));
                ui.separator();
                ui.label(egui::RichText::new("Related annotations to review").strong());
                if preview.related_features.is_empty() {
                    ui.small("No other feature shares either old boundary.");
                } else {
                    egui::ScrollArea::vertical()
                        .max_height(180.0)
                        .show(ui, |ui| {
                            for related in &preview.related_features {
                                let reasons = related
                                    .match_reasons
                                    .iter()
                                    .copied()
                                    .map(Self::related_reason_label)
                                    .collect::<Vec<_>>()
                                    .join(", ");
                                ui.monospace(format!(
                                    "{} {} {} [{}]; not modified",
                                    related.feature_index,
                                    related.feature_kind,
                                    related.location_display,
                                    reasons
                                ));
                            }
                        });
                }
            }
        });
        self.show_feature_location_editor = open;
    }
}
