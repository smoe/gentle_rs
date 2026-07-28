//! Thin GUI adapter for previewing and applying exact feature-boundary edits.

use super::*;
use gentle_protocol::{
    FeatureLocationCompoundWarning, FeatureLocationEditReport, FeatureLocationEditRequest,
    FeatureLocationEditStrand, FeatureLocationIntervalBoundaryRole, FeatureRecordCreateRequest,
    FeatureRecordCurationOutcome, FeatureRecordCurationReport, FeatureRecordCurationRequest,
    FeatureRecordDeleteRequest, FeatureRecordQualifier, FeatureRecordReviewEvidence,
    RelatedFeatureBoundaryReason,
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) enum FeatureEditorMode {
    #[default]
    Location,
    Create,
    Delete,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct FeatureRecordQualifierUiRow {
    pub(super) key: String,
    pub(super) value: String,
    pub(super) has_value: bool,
}

#[derive(Clone, Debug)]
pub(super) struct FeatureRecordCurationUiState {
    pub(super) mode: FeatureEditorMode,
    pub(super) create_feature_kind: String,
    pub(super) create_start_1based: String,
    pub(super) create_end_1based_inclusive: String,
    pub(super) create_strand: FeatureLocationEditStrand,
    pub(super) create_qualifiers: Vec<FeatureRecordQualifierUiRow>,
    pub(super) delete_feature_index: Option<usize>,
    pub(super) preview_request: Option<FeatureRecordCurationRequest>,
    pub(super) preview: Option<FeatureRecordCurationReport>,
    pub(super) status: String,
}

impl Default for FeatureRecordCurationUiState {
    fn default() -> Self {
        Self {
            mode: FeatureEditorMode::Location,
            create_feature_kind: "misc_feature".to_string(),
            create_start_1based: "1".to_string(),
            create_end_1based_inclusive: "1".to_string(),
            create_strand: FeatureLocationEditStrand::Forward,
            create_qualifiers: Vec::new(),
            delete_feature_index: None,
            preview_request: None,
            preview: None,
            status: String::new(),
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct FeatureLocationEditorSegmentOption {
    pub(super) segment_index: usize,
    pub(super) label: String,
    pub(super) start_1based: usize,
    pub(super) end_1based_inclusive: usize,
}

#[derive(Clone, Debug, Default)]
pub(super) struct FeatureLocationEditorUiState {
    pub(super) selected_feature_index: Option<usize>,
    pub(super) selected_segment_index: Option<usize>,
    pub(super) segment_options: Vec<FeatureLocationEditorSegmentOption>,
    pub(super) start_1based: String,
    pub(super) end_1based_inclusive: String,
    pub(super) preview: Option<FeatureLocationEditReport>,
    pub(super) status: String,
    pub(super) record: FeatureRecordCurationUiState,
}

impl MainAreaDna {
    pub(super) fn feature_location_edit_unavailable_reason(
        &self,
        feature_index: usize,
    ) -> Option<String> {
        let Ok(dna) = self.dna.read() else {
            return Some("Sequence state is temporarily unavailable.".to_string());
        };
        let Some(feature) = dna.features().get(feature_index) else {
            return Some("The selected feature is no longer available.".to_string());
        };
        crate::feature_location::editable_feature_location(feature, dna.len(), dna.is_circular())
            .err()
            .map(|error| error.message)
    }

    fn feature_location_editor_options(&self) -> Vec<(usize, String, Option<String>)> {
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
                    crate::feature_location::editable_feature_location(
                        feature,
                        dna.len(),
                        dna.is_circular(),
                    )
                    .err()
                    .map(|error| error.message),
                )
            })
            .collect()
    }

    fn feature_record_editor_options(&self) -> Vec<(usize, String)> {
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
                (index, format!("{index}: {} ({label})", feature.kind))
            })
            .collect()
    }

    pub(super) fn invalidate_feature_record_preview(&mut self) {
        self.feature_location_editor_ui.record.preview_request = None;
        self.feature_location_editor_ui.record.preview = None;
    }

    pub fn focus_feature_record_create_editor(
        &mut self,
        selected_range_0based: Option<(usize, usize)>,
    ) {
        let sequence_len = self.dna.read().map(|dna| dna.len()).unwrap_or_default();
        let range = selected_range_0based
            .filter(|(start, end)| start < end && *end <= sequence_len)
            .or_else(|| (sequence_len > 0).then_some((0, 1)));
        if let Some((start, end)) = range {
            self.feature_location_editor_ui.record.create_start_1based =
                start.saturating_add(1).to_string();
            self.feature_location_editor_ui
                .record
                .create_end_1based_inclusive = end.to_string();
            self.feature_location_editor_ui.record.status.clear();
        } else {
            self.feature_location_editor_ui.record.status =
                "The active sequence has no bases available for an exact feature range."
                    .to_string();
        }
        self.feature_location_editor_ui.record.mode = FeatureEditorMode::Create;
        self.invalidate_feature_record_preview();
        self.show_feature_location_editor = true;
    }

    pub fn focus_feature_record_delete_editor(&mut self, feature_index: Option<usize>) {
        let options = self.feature_record_editor_options();
        self.feature_location_editor_ui.record.delete_feature_index = feature_index
            .filter(|index| options.iter().any(|(candidate, _)| candidate == index))
            .or(self.feature_location_editor_ui.record.delete_feature_index)
            .filter(|index| options.iter().any(|(candidate, _)| candidate == index))
            .or_else(|| {
                self.get_selected_feature_id()
                    .filter(|index| options.iter().any(|(candidate, _)| candidate == index))
            })
            .or_else(|| options.first().map(|(index, _)| *index));
        self.feature_location_editor_ui.record.mode = FeatureEditorMode::Delete;
        self.feature_location_editor_ui.record.status =
            if self
                .feature_location_editor_ui
                .record
                .delete_feature_index
                .is_some()
            {
                String::new()
            } else {
                "No feature record is available to delete.".to_string()
            };
        self.invalidate_feature_record_preview();
        self.show_feature_location_editor = true;
    }

    fn seed_feature_location_editor(&mut self, feature_index: usize) {
        let editability = self.dna.read().ok().and_then(|dna| {
            let feature = dna.features().get(feature_index)?;
            Some(crate::feature_location::editable_feature_location(
                feature,
                dna.len(),
                dna.is_circular(),
            ))
        });
        self.feature_location_editor_ui.selected_feature_index = Some(feature_index);
        self.feature_location_editor_ui.selected_segment_index = None;
        self.feature_location_editor_ui.segment_options.clear();
        self.feature_location_editor_ui.preview = None;
        match editability {
            Some(Ok(crate::feature_location::EditableFeatureLocation::Simple(snapshot))) => {
                self.feature_location_editor_ui.start_1based = snapshot.start_1based.to_string();
                self.feature_location_editor_ui.end_1based_inclusive =
                    snapshot.end_1based_inclusive.to_string();
                self.feature_location_editor_ui.status.clear();
            }
            Some(Ok(crate::feature_location::EditableFeatureLocation::FlatCompound(view))) => {
                self.feature_location_editor_ui.segment_options = view
                    .segments
                    .iter()
                    .map(|segment| FeatureLocationEditorSegmentOption {
                        segment_index: segment.segment_index,
                        label: format!(
                            "Stored {} / biological {}: {}..{}",
                            segment.segment_index + 1,
                            segment.biological_segment_number,
                            segment.snapshot.start_1based,
                            segment.snapshot.end_1based_inclusive
                        ),
                        start_1based: segment.snapshot.start_1based,
                        end_1based_inclusive: segment.snapshot.end_1based_inclusive,
                    })
                    .collect();
                let first = self.feature_location_editor_ui.segment_options[0].clone();
                self.feature_location_editor_ui.selected_segment_index = Some(first.segment_index);
                self.feature_location_editor_ui.start_1based = first.start_1based.to_string();
                self.feature_location_editor_ui.end_1based_inclusive =
                    first.end_1based_inclusive.to_string();
                self.feature_location_editor_ui.status.clear();
            }
            Some(Err(error)) => {
                self.feature_location_editor_ui.start_1based.clear();
                self.feature_location_editor_ui.end_1based_inclusive.clear();
                self.feature_location_editor_ui.status = error.message;
            }
            None => {
                self.feature_location_editor_ui.start_1based.clear();
                self.feature_location_editor_ui.end_1based_inclusive.clear();
                self.feature_location_editor_ui.status =
                    "The selected feature is no longer available.".to_string();
            }
        }
    }

    pub(super) fn select_feature_location_editor_segment(&mut self, segment_index: usize) {
        let Some(option) = self
            .feature_location_editor_ui
            .segment_options
            .iter()
            .find(|option| option.segment_index == segment_index)
            .cloned()
        else {
            self.feature_location_editor_ui.status =
                format!("Segment index {segment_index} is not available.");
            return;
        };
        self.feature_location_editor_ui.selected_segment_index = Some(segment_index);
        self.feature_location_editor_ui.start_1based = option.start_1based.to_string();
        self.feature_location_editor_ui.end_1based_inclusive =
            option.end_1based_inclusive.to_string();
        self.feature_location_editor_ui.preview = None;
        self.feature_location_editor_ui.status.clear();
    }

    pub fn focus_feature_location_editor(&mut self, feature_index: Option<usize>) {
        self.feature_location_editor_ui.record.mode = FeatureEditorMode::Location;
        let options = self.feature_location_editor_options();
        let selected = feature_index
            .filter(|index| options.iter().any(|(candidate, _, _)| candidate == index))
            .or(self.feature_location_editor_ui.selected_feature_index)
            .or_else(|| {
                self.get_selected_feature_id()
                    .filter(|index| options.iter().any(|(candidate, _, _)| candidate == index))
            })
            .or_else(|| {
                options.iter().find_map(|(index, _, unavailable_reason)| {
                    unavailable_reason.is_none().then_some(*index)
                })
            });
        if let Some(selected) = selected {
            self.seed_feature_location_editor(selected);
        } else {
            self.feature_location_editor_ui.status =
                "No editable exact simple or flat compound feature is available.".to_string();
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
            segment_index: self.feature_location_editor_ui.selected_segment_index,
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
            || preview
                .compound_context
                .as_ref()
                .map(|context| context.segment_index)
                != request.segment_index
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

    fn feature_record_create_request(
        &self,
        expected_annotation_state_fingerprint_sha256: Option<String>,
    ) -> Result<FeatureRecordCurationRequest, String> {
        let seq_id = self
            .seq_id
            .clone()
            .ok_or_else(|| "No active sequence is attached to this window".to_string())?;
        let start_1based = self
            .feature_location_editor_ui
            .record
            .create_start_1based
            .trim()
            .parse::<i64>()
            .map_err(|_| "Start must be a positive 1-based coordinate".to_string())?;
        let end_1based_inclusive = self
            .feature_location_editor_ui
            .record
            .create_end_1based_inclusive
            .trim()
            .parse::<i64>()
            .map_err(|_| "End must be a positive 1-based inclusive coordinate".to_string())?;
        if start_1based < 1 || end_1based_inclusive < 1 {
            return Err("Start and end must both be at least 1".to_string());
        }
        if end_1based_inclusive < start_1based {
            return Err("End must be greater than or equal to start".to_string());
        }
        let feature_kind = self
            .feature_location_editor_ui
            .record
            .create_feature_kind
            .trim();
        if feature_kind.is_empty() {
            return Err("Feature kind must not be empty".to_string());
        }
        let qualifiers = self
            .feature_location_editor_ui
            .record
            .create_qualifiers
            .iter()
            .map(|row| FeatureRecordQualifier {
                key: row.key.clone(),
                value: row.has_value.then(|| row.value.clone()),
            })
            .collect();
        Ok(FeatureRecordCurationRequest::Create(
            FeatureRecordCreateRequest {
                seq_id,
                feature_kind: feature_kind.to_string(),
                start_0based: start_1based - 1,
                end_0based_exclusive: end_1based_inclusive,
                strand: self.feature_location_editor_ui.record.create_strand,
                qualifiers,
                expected_annotation_state_fingerprint_sha256,
            },
        ))
    }

    fn feature_record_delete_request(
        &self,
        expected_feature_fingerprint_sha256: Option<String>,
        expected_annotation_state_fingerprint_sha256: Option<String>,
    ) -> Result<FeatureRecordCurationRequest, String> {
        let seq_id = self
            .seq_id
            .clone()
            .ok_or_else(|| "No active sequence is attached to this window".to_string())?;
        let feature_index = self
            .feature_location_editor_ui
            .record
            .delete_feature_index
            .ok_or_else(|| "Select a feature record first".to_string())?;
        Ok(FeatureRecordCurationRequest::Delete(
            FeatureRecordDeleteRequest {
                seq_id,
                feature_index,
                expected_feature_fingerprint_sha256,
                expected_annotation_state_fingerprint_sha256,
            },
        ))
    }

    pub(super) fn run_feature_record_preview(&mut self) {
        let request = match self.feature_location_editor_ui.record.mode {
            FeatureEditorMode::Create => self.feature_record_create_request(None),
            FeatureEditorMode::Delete => self.feature_record_delete_request(None, None),
            FeatureEditorMode::Location => {
                self.feature_location_editor_ui.record.status =
                    "Choose Create or Delete first.".to_string();
                return;
            }
        };
        let request = match request {
            Ok(request) => request,
            Err(error) => {
                self.feature_location_editor_ui.record.status = error;
                self.invalidate_feature_record_preview();
                return;
            }
        };
        let Some(result) = self.apply_operation_with_feedback_and_result(
            Operation::PreviewFeatureRecordCuration {
                request: request.clone(),
            },
        ) else {
            self.feature_location_editor_ui.record.status = self.op_status.clone();
            self.invalidate_feature_record_preview();
            return;
        };
        self.feature_location_editor_ui.record.preview_request = Some(request);
        self.feature_location_editor_ui.record.preview = result
            .feature_record_curation_report
            .map(|report| report.as_ref().clone());
        self.feature_location_editor_ui.record.status =
            "Preview is current. Review the complete record and related annotations."
                .to_string();
    }

    pub(super) fn run_feature_record_apply(&mut self) {
        let Some(preview) = self.feature_location_editor_ui.record.preview.clone() else {
            self.feature_location_editor_ui.record.status =
                "Run Preview before applying feature-record curation.".to_string();
            return;
        };
        let Some(preview_request) = self
            .feature_location_editor_ui
            .record
            .preview_request
            .clone()
        else {
            self.feature_location_editor_ui.record.status =
                "Preview request is unavailable; preview again.".to_string();
            return;
        };
        let current_request = match self.feature_location_editor_ui.record.mode {
            FeatureEditorMode::Create => self.feature_record_create_request(None),
            FeatureEditorMode::Delete => self.feature_record_delete_request(None, None),
            FeatureEditorMode::Location => {
                self.feature_location_editor_ui.record.status =
                    "Choose Create or Delete first.".to_string();
                return;
            }
        };
        let current_request = match current_request {
            Ok(request) => request,
            Err(error) => {
                self.feature_location_editor_ui.record.status = error;
                return;
            }
        };
        if current_request != preview_request {
            self.feature_location_editor_ui.record.status =
                "Feature record fields changed after preview; preview again before applying."
                    .to_string();
            self.invalidate_feature_record_preview();
            return;
        }
        let request = match current_request {
            FeatureRecordCurationRequest::Create(mut request) => {
                request.expected_annotation_state_fingerprint_sha256 = Some(
                    preview
                        .before_annotation_state_fingerprint_sha256
                        .clone(),
                );
                FeatureRecordCurationRequest::Create(request)
            }
            FeatureRecordCurationRequest::Delete(mut request) => {
                let FeatureRecordCurationOutcome::Delete {
                    deleted_feature, ..
                } = &preview.outcome
                else {
                    self.feature_location_editor_ui.record.status =
                        "Preview outcome no longer matches Delete mode.".to_string();
                    self.invalidate_feature_record_preview();
                    return;
                };
                request.expected_feature_fingerprint_sha256 =
                    Some(deleted_feature.feature_fingerprint_sha256.clone());
                request.expected_annotation_state_fingerprint_sha256 = Some(
                    preview
                        .before_annotation_state_fingerprint_sha256
                        .clone(),
                );
                FeatureRecordCurationRequest::Delete(request)
            }
        };
        let Some(result) = self.apply_operation_with_feedback_and_result(
            Operation::ApplyFeatureRecordCuration { request },
        ) else {
            self.feature_location_editor_ui.record.status = self.op_status.clone();
            return;
        };
        let outcome = result
            .feature_record_curation_report
            .as_ref()
            .map(|report| report.outcome.clone());
        let status = result
            .messages
            .first()
            .cloned()
            .unwrap_or_else(|| "Feature-record curation applied.".to_string());
        self.invalidate_feature_record_preview();
        match outcome {
            Some(FeatureRecordCurationOutcome::Create {
                created_feature_index: Some(feature_index),
                ..
            }) => {
                self.focus_feature(feature_index);
            }
            Some(FeatureRecordCurationOutcome::Delete { .. }) => {
                self.feature_location_editor_ui.record.delete_feature_index =
                    self.feature_record_editor_options().first().map(|(index, _)| *index);
            }
            _ => {}
        }
        self.feature_location_editor_ui.record.status = status;
    }

    fn related_reason_label(reason: RelatedFeatureBoundaryReason) -> &'static str {
        match reason {
            RelatedFeatureBoundaryReason::SharesOldStartBoundary => "shares old start",
            RelatedFeatureBoundaryReason::SharesOldEndBoundary => "shares old end",
        }
    }

    fn boundary_role_label(role: FeatureLocationIntervalBoundaryRole) -> &'static str {
        match role {
            FeatureLocationIntervalBoundaryRole::Start0Based => "start",
            FeatureLocationIntervalBoundaryRole::End0BasedExclusive => "end",
        }
    }

    fn compound_warning_label(warning: &FeatureLocationCompoundWarning) -> String {
        match warning {
            FeatureLocationCompoundWarning::OverlappingSegments {
                edited_segment_index,
                overlapping_segment_index,
            } => format!(
                "Segment {} overlaps segment {}.",
                edited_segment_index + 1,
                overlapping_segment_index + 1
            ),
            FeatureLocationCompoundWarning::EstablishedDirectionBroken {
                first_segment_index,
                second_segment_index,
                established_direction,
            } => format!(
                "Stored {:?} direction is broken between segments {} and {}.",
                established_direction,
                first_segment_index + 1,
                second_segment_index + 1
            ),
            FeatureLocationCompoundWarning::CdsCodingLengthDeltaNotDivisibleByThree {
                signed_length_delta_nt,
            } => format!(
                "CDS length changes by {signed_length_delta_nt} nt; /codon_start is unchanged."
            ),
        }
    }

    fn render_feature_record_preview(
        ui: &mut egui::Ui,
        preview: &FeatureRecordCurationReport,
    ) {
        ui.separator();
        ui.label(egui::RichText::new("Preview").strong());
        match &preview.outcome {
            FeatureRecordCurationOutcome::Create {
                proposed_feature,
                created_feature_index,
            } => {
                ui.monospace(format!(
                    "Create {} at {}",
                    proposed_feature.feature_kind, proposed_feature.location_display
                ));
                ui.small(match created_feature_index {
                    Some(index) => format!("Applied feature index: {index}"),
                    None => "Preview does not reserve a feature index.".to_string(),
                });
                for qualifier in &proposed_feature.qualifiers {
                    ui.monospace(match qualifier.value.as_deref() {
                        Some(value) => format!("/{}={value}", qualifier.key),
                        None => format!("/{}", qualifier.key),
                    });
                }
            }
            FeatureRecordCurationOutcome::Delete {
                deleted_feature_index,
                deleted_feature,
                shifted_feature_count,
            } => {
                ui.monospace(format!(
                    "Delete {}: {} at {}",
                    deleted_feature_index,
                    deleted_feature.feature_kind,
                    deleted_feature.location_display
                ));
                ui.small(format!(
                    "{} subsequent feature index{} will shift down by one.",
                    shifted_feature_count,
                    if *shifted_feature_count == 1 {
                        ""
                    } else {
                        "es"
                    }
                ));
                for qualifier in &deleted_feature.qualifiers {
                    ui.monospace(match qualifier.value.as_deref() {
                        Some(value) => format!("/{}={value}", qualifier.key),
                        None => format!("/{}", qualifier.key),
                    });
                }
            }
        }
        ui.small(format!(
            "annotation state: {} -> {}",
            preview.before_annotation_state_fingerprint_sha256,
            preview.after_annotation_state_fingerprint_sha256
        ));
        ui.separator();
        ui.label(egui::RichText::new("Annotations to review").strong());
        if preview.review_candidates.is_empty() {
            ui.small("No location overlap or recognized shared identifier was found.");
        } else {
            egui::ScrollArea::vertical()
                .max_height(160.0)
                .show(ui, |ui| {
                    for candidate in &preview.review_candidates {
                        let evidence = candidate
                            .evidence
                            .iter()
                            .map(|evidence| match evidence {
                                FeatureRecordReviewEvidence::OverlappingLocation => {
                                    "location overlap".to_string()
                                }
                                FeatureRecordReviewEvidence::SharedIdentifier {
                                    qualifier_key,
                                    qualifier_value,
                                    ..
                                } => {
                                    format!("shared /{qualifier_key}={qualifier_value}")
                                }
                            })
                            .collect::<Vec<_>>()
                            .join(", ");
                        ui.monospace(format!(
                            "{} {} {} [{}]; not modified",
                            candidate.feature_index,
                            candidate.feature_kind,
                            candidate.location_display,
                            evidence
                        ));
                    }
                });
        }
        ui.small(
            "Overlap and shared identifiers are informational review evidence, not dependency claims.",
        );
    }

    fn render_feature_record_create_editor(&mut self, ui: &mut egui::Ui) {
        ui.label("Append one exact simple GenBank feature after preview.");
        let mut changed = false;
        ui.horizontal(|ui| {
            ui.label("Kind");
            changed |= ui
                .text_edit_singleline(
                    &mut self
                        .feature_location_editor_ui
                        .record
                        .create_feature_kind,
                )
                .changed();
            ui.label("Strand");
            let before = self.feature_location_editor_ui.record.create_strand;
            egui::ComboBox::from_id_salt((
                "feature_record_create_strand",
                self.panel_scope_key(),
            ))
            .selected_text(match before {
                FeatureLocationEditStrand::Forward => "Forward",
                FeatureLocationEditStrand::Reverse => "Reverse",
            })
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut self.feature_location_editor_ui.record.create_strand,
                    FeatureLocationEditStrand::Forward,
                    "Forward",
                );
                ui.selectable_value(
                    &mut self.feature_location_editor_ui.record.create_strand,
                    FeatureLocationEditStrand::Reverse,
                    "Reverse",
                );
            });
            changed |= before != self.feature_location_editor_ui.record.create_strand;
        });
        ui.horizontal(|ui| {
            ui.label("Start (1-based)");
            changed |= ui
                .text_edit_singleline(
                    &mut self
                        .feature_location_editor_ui
                        .record
                        .create_start_1based,
                )
                .changed();
            ui.label("End (1-based, inclusive)");
            changed |= ui
                .text_edit_singleline(
                    &mut self
                        .feature_location_editor_ui
                        .record
                        .create_end_1based_inclusive,
                )
                .changed();
        });
        ui.label(egui::RichText::new("Ordered qualifiers").strong());
        let mut remove_row = None;
        for (index, row) in self
            .feature_location_editor_ui
            .record
            .create_qualifiers
            .iter_mut()
            .enumerate()
        {
            ui.horizontal(|ui| {
                changed |= ui
                    .add(egui::TextEdit::singleline(&mut row.key).hint_text("key"))
                    .changed();
                changed |= ui.checkbox(&mut row.has_value, "Value").changed();
                if row.has_value {
                    changed |= ui
                        .add(egui::TextEdit::singleline(&mut row.value).hint_text("value"))
                        .changed();
                }
                if ui.small_button("x").on_hover_text("Remove qualifier").clicked() {
                    remove_row = Some(index);
                }
            });
        }
        if let Some(index) = remove_row {
            self.feature_location_editor_ui
                .record
                .create_qualifiers
                .remove(index);
            changed = true;
        }
        if ui
            .small_button("+")
            .on_hover_text("Add qualifier")
            .clicked()
        {
            self.feature_location_editor_ui
                .record
                .create_qualifiers
                .push(FeatureRecordQualifierUiRow::default());
            changed = true;
        }
        if changed {
            self.invalidate_feature_record_preview();
        }
        ui.horizontal(|ui| {
            if ui.button("Preview").clicked() {
                self.run_feature_record_preview();
            }
            let apply = ui.add_enabled(
                self.feature_location_editor_ui.record.preview.is_some(),
                egui::Button::new("Create"),
            );
            if apply
                .on_hover_text("Create is enabled only for the exact record just previewed")
                .clicked()
            {
                self.run_feature_record_apply();
            }
        });
        if !self.feature_location_editor_ui.record.status.is_empty() {
            ui.label(&self.feature_location_editor_ui.record.status);
        }
        if let Some(preview) = self.feature_location_editor_ui.record.preview.as_ref() {
            Self::render_feature_record_preview(ui, preview);
        }
    }

    fn render_feature_record_delete_editor(
        &mut self,
        ui: &mut egui::Ui,
        options: &[(usize, String)],
    ) {
        ui.label("Delete one complete feature record of any existing location shape.");
        let selected_label = self
            .feature_location_editor_ui
            .record
            .delete_feature_index
            .and_then(|selected| {
                options
                    .iter()
                    .find(|(index, _)| *index == selected)
                    .map(|(_, label)| label.clone())
            })
            .unwrap_or_else(|| "Select feature".to_string());
        let mut next_selected = None;
        egui::ComboBox::from_id_salt((
            "feature_record_delete_feature",
            self.panel_scope_key(),
        ))
        .selected_text(selected_label)
        .width(430.0)
        .show_ui(ui, |ui| {
            for (index, label) in options {
                if ui
                    .selectable_label(
                        self.feature_location_editor_ui
                            .record
                            .delete_feature_index
                            == Some(*index),
                        label,
                    )
                    .clicked()
                {
                    next_selected = Some(*index);
                }
            }
        });
        if let Some(index) = next_selected {
            self.feature_location_editor_ui.record.delete_feature_index = Some(index);
            self.invalidate_feature_record_preview();
            self.feature_location_editor_ui.record.status.clear();
        }
        ui.horizontal(|ui| {
            if ui.button("Preview").clicked() {
                self.run_feature_record_preview();
            }
            let apply = ui.add_enabled(
                self.feature_location_editor_ui.record.preview.is_some(),
                egui::Button::new("Delete"),
            );
            if apply
                .on_hover_text("Delete is enabled only for the exact record just previewed")
                .clicked()
            {
                self.run_feature_record_apply();
            }
        });
        if !self.feature_location_editor_ui.record.status.is_empty() {
            ui.label(&self.feature_location_editor_ui.record.status);
        }
        if let Some(preview) = self.feature_location_editor_ui.record.preview.as_ref() {
            Self::render_feature_record_preview(ui, preview);
        }
    }

    pub(super) fn render_feature_location_editor(&mut self, ctx: &egui::Context) {
        if !self.show_feature_location_editor {
            return;
        }
        let options = self.feature_location_editor_options();
        let record_options = self.feature_record_editor_options();
        let mut open = self.show_feature_location_editor;
        let spec = crate::egui_compat::HostedWindowSpec::new(
            "Feature Editor",
            egui::Id::new(("feature_location_editor", self.panel_scope_key())),
            egui::vec2(680.0, 560.0),
            egui::vec2(460.0, 340.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            let previous_mode = self.feature_location_editor_ui.record.mode;
            ui.horizontal(|ui| {
                ui.selectable_value(
                    &mut self.feature_location_editor_ui.record.mode,
                    FeatureEditorMode::Location,
                    "Location",
                );
                ui.selectable_value(
                    &mut self.feature_location_editor_ui.record.mode,
                    FeatureEditorMode::Create,
                    "Create",
                );
                ui.selectable_value(
                    &mut self.feature_location_editor_ui.record.mode,
                    FeatureEditorMode::Delete,
                    "Delete",
                );
            });
            if previous_mode != self.feature_location_editor_ui.record.mode {
                self.feature_location_editor_ui.preview = None;
                self.invalidate_feature_record_preview();
                self.feature_location_editor_ui.record.status.clear();
                if self.feature_location_editor_ui.record.mode == FeatureEditorMode::Delete
                    && self
                        .feature_location_editor_ui
                        .record
                        .delete_feature_index
                        .is_none()
                {
                    self.feature_location_editor_ui.record.delete_feature_index =
                        record_options.first().map(|(index, _)| *index);
                }
            }
            ui.separator();
            match self.feature_location_editor_ui.record.mode {
                FeatureEditorMode::Create => {
                    self.render_feature_record_create_editor(ui);
                    return;
                }
                FeatureEditorMode::Delete => {
                    self.render_feature_record_delete_editor(ui, &record_options);
                    return;
                }
                FeatureEditorMode::Location => {}
            }
            ui.label("Edit one exact simple range or one segment in a flat Join/Order annotation.");
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
                for (index, label, unavailable_reason) in &options {
                    let response = ui.add_enabled(
                        unavailable_reason.is_none(),
                        egui::Button::new(label).selected(
                            self.feature_location_editor_ui.selected_feature_index == Some(*index),
                        ),
                    );
                    let response = if let Some(reason) = unavailable_reason {
                        response.on_disabled_hover_text(reason)
                    } else {
                        response
                    };
                    if response.clicked() {
                        next_selected = Some(*index);
                    }
                }
            });
            if let Some(index) = next_selected {
                self.seed_feature_location_editor(index);
            }
            if !self.feature_location_editor_ui.segment_options.is_empty() {
                let selected_segment_label = self
                    .feature_location_editor_ui
                    .selected_segment_index
                    .and_then(|selected| {
                        self.feature_location_editor_ui
                            .segment_options
                            .iter()
                            .find(|option| option.segment_index == selected)
                    })
                    .map(|option| option.label.clone())
                    .unwrap_or_else(|| "Select segment".to_string());
                let mut next_segment = None;
                egui::ComboBox::from_id_salt((
                    "feature_location_editor_segment",
                    self.panel_scope_key(),
                ))
                .selected_text(selected_segment_label)
                .width(430.0)
                .show_ui(ui, |ui| {
                    for option in &self.feature_location_editor_ui.segment_options {
                        if ui
                            .selectable_label(
                                self.feature_location_editor_ui.selected_segment_index
                                    == Some(option.segment_index),
                                &option.label,
                            )
                            .clicked()
                        {
                            next_segment = Some(option.segment_index);
                        }
                    }
                });
                if let Some(segment_index) = next_segment {
                    self.select_feature_location_editor_segment(segment_index);
                }
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
                if let Some(context) = preview.compound_context.as_ref() {
                    ui.small(format!(
                        "{:?}; stored segment {} of {}; biological segment {}",
                        context.container_kind,
                        context.segment_index + 1,
                        context.total_segments,
                        context.biological_segment_number
                    ));
                }
                for warning in &preview.compound_validation_warnings {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 120, 0),
                        Self::compound_warning_label(warning),
                    );
                }
                ui.separator();
                ui.label(egui::RichText::new("Related annotations to review").strong());
                if preview.related_features.is_empty()
                    && preview.related_segment_boundaries.is_empty()
                {
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
                            for related in &preview.related_segment_boundaries {
                                ui.monospace(format!(
                                    "{} {}{}: edited {} = related {} at boundary {}; not modified",
                                    related.related_feature_index,
                                    related.related_feature_kind,
                                    related
                                        .related_segment_index
                                        .map(|index| format!(" segment {}", index + 1))
                                        .unwrap_or_default(),
                                    Self::boundary_role_label(related.edited_boundary),
                                    Self::boundary_role_label(related.related_boundary),
                                    related.boundary_coordinate_0based
                                ));
                            }
                        });
                }
            }
        });
        self.show_feature_location_editor = open;
    }
}
