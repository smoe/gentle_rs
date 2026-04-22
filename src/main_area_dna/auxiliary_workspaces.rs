//! Dotplot, splicing, and auxiliary workspace support for `MainAreaDna`.
//!
//! This submodule groups the heavyweight secondary workspaces and their shared
//! dotplot support so the main sequence-window file can shrink toward clearer
//! GUI slices ahead of a later crate split.
//!
//! Look here for:
//! - dotplot UI state and export helpers
//! - secondary workspace rendering shared by splicing/RNA-read/dotplot panes
//! - window-sizing/default helpers for extracted auxiliary tools

use super::*;
use crate::engine::DotplotOverlayQuerySpec;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct DotplotOpsUiState {
    pub(super) half_window_bp: String,
    pub(super) word_size: String,
    pub(super) step_bp: String,
    pub(super) max_mismatches: String,
    pub(super) tile_bp: String,
    pub(super) mode: DotplotMode,
    pub(super) overlay_enabled: bool,
    pub(super) overlay_transcript_feature_ids: Vec<usize>,
    pub(super) overlay_x_axis_mode: DotplotOverlayXAxisMode,
    pub(super) overlay_anchor_exon: Option<DotplotOverlayAnchorExonRef>,
    pub(super) dotplot_id: String,
    pub(super) display_density_threshold: f32,
    pub(super) display_intensity_gain: f32,
    pub(super) reference_seq_id: String,
    pub(super) reference_span_start_0based: String,
    pub(super) reference_span_end_0based: String,
    pub(super) show_flexibility_track: bool,
    pub(super) flex_model: FlexibilityModel,
    pub(super) flex_bin_bp: String,
    pub(super) flex_smoothing_bp: String,
    pub(super) flex_track_id: String,
}

impl Default for DotplotOpsUiState {
    fn default() -> Self {
        Self {
            half_window_bp: String::new(),
            word_size: "7".to_string(),
            step_bp: "1".to_string(),
            max_mismatches: "0".to_string(),
            tile_bp: String::new(),
            mode: DotplotMode::SelfReverseComplement,
            overlay_enabled: false,
            overlay_transcript_feature_ids: vec![],
            overlay_x_axis_mode: DotplotOverlayXAxisMode::PercentLength,
            overlay_anchor_exon: None,
            dotplot_id: "dotplot_primary".to_string(),
            display_density_threshold: 0.0,
            display_intensity_gain: 1.0,
            reference_seq_id: String::new(),
            reference_span_start_0based: String::new(),
            reference_span_end_0based: String::new(),
            show_flexibility_track: true,
            flex_model: FlexibilityModel::AtRichness,
            flex_bin_bp: "25".to_string(),
            flex_smoothing_bp: "75".to_string(),
            flex_track_id: "flex_primary".to_string(),
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct DotplotComputeDiagnostics {
    pub(super) seq_id: String,
    pub(super) reference_seq_id: Option<String>,
    pub(super) mode: DotplotMode,
    pub(super) query_span_start_0based: usize,
    pub(super) query_span_end_0based: usize,
    pub(super) query_span_bp: usize,
    pub(super) reference_span_start_0based: usize,
    pub(super) reference_span_end_0based: usize,
    pub(super) reference_span_bp: usize,
    pub(super) word_size: usize,
    pub(super) step_bp: usize,
    pub(super) max_mismatches: usize,
    pub(super) tile_bp: Option<usize>,
    pub(super) query_windows: usize,
    pub(super) reference_windows: usize,
    pub(super) estimated_pair_evaluations: usize,
}

#[derive(Clone, Debug)]
pub(super) struct DotplotOverlayChoice {
    pub(super) transcript_feature_id: usize,
    pub(super) transcript_id: String,
    pub(super) label: String,
    pub(super) mode: DotplotMode,
    pub(super) estimated_length_bp: usize,
    pub(super) derived_seq_id: String,
    pub(super) color_rgb: [u8; 3],
}

impl MainAreaDna {
    const DOTPLOT_AUTO_COMPUTE_DEBOUNCE_MS: u64 = 300;

    pub(super) fn dotplot_window_identity_seq_id(&self) -> Option<String> {
        self.seq_id
            .clone()
            .or_else(|| self.current_dotplot_query_seq_id())
    }

    pub(super) fn dotplot_overlay_series_color(index: usize) -> [u8; 3] {
        const PALETTE: [[u8; 3]; 8] = [
            [29, 78, 216],
            [220, 38, 38],
            [5, 150, 105],
            [217, 119, 6],
            [124, 58, 237],
            [190, 24, 93],
            [8, 145, 178],
            [71, 85, 105],
        ];
        PALETTE[index % PALETTE.len()]
    }

    pub(super) fn normalize_dotplot_transcript_token(raw: &str) -> String {
        let mut out = String::new();
        for c in raw.chars() {
            if c.is_ascii_alphanumeric() {
                out.push(c.to_ascii_lowercase());
            } else if matches!(c, '_' | '-' | '.') {
                out.push('_');
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "feature".to_string()
        } else {
            trimmed.to_string()
        }
    }

    pub(super) fn schedule_dotplot_auto_compute(&mut self) {
        self.dotplot_auto_compute_dirty = true;
        self.dotplot_auto_compute_due_at =
            Some(Instant::now() + Duration::from_millis(Self::DOTPLOT_AUTO_COMPUTE_DEBOUNCE_MS));
    }

    pub(super) fn clear_dotplot_auto_compute(&mut self) {
        self.dotplot_auto_compute_dirty = false;
        self.dotplot_auto_compute_due_at = None;
    }

    pub(super) fn dotplot_overlay_reference_view(&self) -> Option<Arc<SplicingExpertView>> {
        if self.using_dotplot_query_override()
            || !Self::dotplot_mode_requires_reference(self.dotplot_ui.mode)
        {
            return None;
        }
        let reference_seq_id = self.dotplot_ui.reference_seq_id.trim();
        self.current_dotplot_reference_splicing_view()
            .filter(|view| view.seq_id == reference_seq_id)
    }

    pub(super) fn dotplot_overlay_choices(&self) -> Vec<DotplotOverlayChoice> {
        let Some(view) = self.dotplot_overlay_reference_view() else {
            return vec![];
        };
        view.transcripts
            .iter()
            .enumerate()
            .map(|(index, lane)| {
                let label = if lane.label.trim().is_empty() {
                    lane.transcript_id.trim().to_string()
                } else {
                    lane.label.trim().to_string()
                };
                let estimated_length_bp = lane
                    .exons
                    .iter()
                    .map(|exon| {
                        exon.end_1based
                            .saturating_add(1)
                            .saturating_sub(exon.start_1based)
                    })
                    .sum::<usize>();
                let transcript_token =
                    Self::normalize_dotplot_transcript_token(&lane.transcript_id);
                DotplotOverlayChoice {
                    transcript_feature_id: lane.transcript_feature_id,
                    transcript_id: lane.transcript_id.clone(),
                    label,
                    mode: if lane.strand.trim() == "-" {
                        DotplotMode::PairReverseComplement
                    } else {
                        DotplotMode::PairForward
                    },
                    estimated_length_bp,
                    derived_seq_id: format!(
                        "{}__mrna__f{}__{}",
                        view.seq_id,
                        lane.transcript_feature_id + 1,
                        transcript_token
                    ),
                    color_rgb: Self::dotplot_overlay_series_color(index),
                }
            })
            .collect()
    }

    pub(super) fn sync_dotplot_overlay_selection(&mut self) -> bool {
        let choices = self.dotplot_overlay_choices();
        let allowed = choices
            .iter()
            .map(|choice| choice.transcript_feature_id)
            .collect::<BTreeSet<_>>();
        let before = self.dotplot_ui.overlay_transcript_feature_ids.clone();
        self.dotplot_ui
            .overlay_transcript_feature_ids
            .retain(|feature_id| allowed.contains(feature_id));
        if self.dotplot_ui.overlay_enabled
            && !choices.is_empty()
            && self.dotplot_ui.overlay_transcript_feature_ids.is_empty()
        {
            self.dotplot_ui.overlay_transcript_feature_ids = choices
                .iter()
                .map(|choice| choice.transcript_feature_id)
                .collect();
        }
        self.dotplot_ui
            .overlay_transcript_feature_ids
            .sort_unstable();
        self.dotplot_ui.overlay_transcript_feature_ids.dedup();
        before != self.dotplot_ui.overlay_transcript_feature_ids
    }

    pub(super) fn dotplot_overlay_anchor_choices(
        &self,
    ) -> Vec<(DotplotOverlayAnchorExonRef, usize)> {
        let Some(view) = self.dotplot_overlay_reference_view() else {
            return vec![];
        };
        let selected_feature_ids = self
            .dotplot_ui
            .overlay_transcript_feature_ids
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        if selected_feature_ids.len() < 2 {
            return vec![];
        }
        view.unique_exons
            .iter()
            .enumerate()
            .filter_map(|(idx, exon)| {
                let support_count = view
                    .matrix_rows
                    .iter()
                    .filter(|row| {
                        selected_feature_ids.contains(&row.transcript_feature_id)
                            && row.exon_presence.get(idx).copied().unwrap_or(false)
                    })
                    .count();
                if support_count < 2 {
                    None
                } else {
                    Some((
                        DotplotOverlayAnchorExonRef {
                            start_1based: exon.start_1based,
                            end_1based: exon.end_1based,
                        },
                        support_count,
                    ))
                }
            })
            .collect()
    }

    pub(super) fn sync_dotplot_overlay_anchor_selection(&mut self) -> bool {
        let choices = self.dotplot_overlay_anchor_choices();
        let before = self.dotplot_ui.overlay_anchor_exon.clone();
        if let Some(selected) = self.dotplot_ui.overlay_anchor_exon.as_ref()
            && !choices.iter().any(|(exon, _)| exon == selected)
        {
            self.dotplot_ui.overlay_anchor_exon = None;
        }
        before != self.dotplot_ui.overlay_anchor_exon
    }

    pub(super) fn current_dotplot_query_seq_id(&self) -> Option<String> {
        let override_id = self.dotplot_query_override_seq_id.trim();
        if !override_id.is_empty() {
            Some(override_id.to_string())
        } else {
            self.seq_id.clone()
        }
    }

    pub(super) fn using_dotplot_query_override(&self) -> bool {
        !self.dotplot_query_override_seq_id.trim().is_empty()
    }

    pub(super) fn current_dotplot_query_label(&self) -> String {
        if self.using_dotplot_query_override() {
            let query_id = self.dotplot_query_override_seq_id.trim();
            let source = self.dotplot_query_override_source_label.trim();
            if source.is_empty() {
                query_id.to_string()
            } else {
                format!("{query_id} ({source})")
            }
        } else {
            self.seq_id
                .clone()
                .unwrap_or_else(|| "<no-seq-id>".to_string())
        }
    }

    pub(super) fn current_dotplot_query_length(&self) -> usize {
        if let Some(query_seq_id) = self.current_dotplot_query_seq_id() {
            if self.seq_id.as_deref() == Some(query_seq_id.as_str()) {
                return self.dna.read().ok().map(|dna| dna.len()).unwrap_or(0);
            }
            return self
                .engine
                .as_ref()
                .and_then(|engine| {
                    engine.read().ok().and_then(|guard| {
                        guard
                            .state()
                            .sequences
                            .get(&query_seq_id)
                            .map(|dna| dna.len())
                    })
                })
                .unwrap_or(0);
        }
        0
    }

    pub(super) fn current_dotplot_owner_seq_id(&self) -> Option<String> {
        if self.dotplot_ui.overlay_enabled
            && let Some(reference_seq_id) =
                Some(self.dotplot_ui.reference_seq_id.trim()).filter(|value| !value.is_empty())
        {
            return Some(reference_seq_id.to_string());
        }
        self.current_dotplot_query_seq_id()
    }

    pub(super) fn dotplot_view_is_overlay(view: &DotplotView) -> bool {
        view.series_count > 1 || view.query_series.len() > 1
    }

    pub(super) fn dotplot_view_total_point_count(view: &DotplotView) -> usize {
        if Self::dotplot_view_is_overlay(view) {
            view.query_series
                .iter()
                .map(|series| series.point_count)
                .sum()
        } else {
            view.point_count
        }
    }

    pub(super) fn dotplot_overlay_max_query_span(view: &DotplotView) -> usize {
        view.query_series
            .iter()
            .map(|series| {
                series
                    .span_end_0based
                    .saturating_sub(series.span_start_0based)
                    .max(1)
            })
            .max()
            .unwrap_or(1)
    }

    pub(super) fn clear_dotplot_query_override(&mut self) {
        self.dotplot_query_override_seq_id.clear();
        self.dotplot_query_override_source_label.clear();
        self.invalidate_dotplot_cache();
    }

    pub(super) fn dotplot_selection_sync_enabled_for_view(&self, view: &DotplotView) -> bool {
        if Self::dotplot_view_is_overlay(view) {
            return false;
        }
        self.seq_id.as_deref() == Some(view.seq_id.as_str())
    }

    pub(super) fn dotplot_ids_for_active_sequence(&self) -> Vec<String> {
        let Some(seq_id) = self.current_dotplot_owner_seq_id() else {
            return vec![];
        };
        let Some(engine) = self.engine.as_ref() else {
            return vec![];
        };
        let Ok(guard) = engine.try_read() else {
            return vec![];
        };
        guard
            .list_dotplot_views(Some(seq_id.as_str()))
            .into_iter()
            .map(|row| row.dotplot_id)
            .collect()
    }

    pub(super) fn flexibility_track_ids_for_active_sequence(&self) -> Vec<String> {
        let Some(seq_id) = self.current_dotplot_owner_seq_id() else {
            return vec![];
        };
        let Some(engine) = self.engine.as_ref() else {
            return vec![];
        };
        let Ok(guard) = engine.try_read() else {
            return vec![];
        };
        guard
            .list_flexibility_tracks(Some(seq_id.as_str()))
            .into_iter()
            .map(|row| row.track_id)
            .collect()
    }

    pub(super) fn dotplot_reference_sequence_ids(&self) -> Vec<String> {
        let Some(engine) = self.engine.as_ref() else {
            return vec![];
        };
        let Ok(guard) = engine.try_read() else {
            return vec![];
        };
        let mut ids: Vec<String> = guard.state().sequences.keys().cloned().collect();
        ids.sort_by_key(|value| value.to_ascii_lowercase());
        ids
    }

    pub(super) fn current_dotplot_reference_splicing_view(
        &self,
    ) -> Option<Arc<SplicingExpertView>> {
        self.rna_read_mapping_window_view
            .clone()
            .or_else(|| self.splicing_expert_window_view.clone())
    }

    pub(super) fn dotplot_transcript_reference_choices(&self) -> Vec<(usize, String)> {
        let Some(view) = self.current_dotplot_reference_splicing_view() else {
            return vec![];
        };
        view.transcripts
            .iter()
            .map(|lane| {
                let label = if lane.label.trim().is_empty() {
                    lane.transcript_id.trim()
                } else {
                    lane.label.trim()
                };
                (
                    lane.transcript_feature_id,
                    format!(
                        "n-{} {} ({})",
                        lane.transcript_feature_id + 1,
                        label,
                        lane.strand.trim()
                    ),
                )
            })
            .collect()
    }

    pub(super) fn set_dotplot_reference_to_sequence(
        &mut self,
        reference_seq_id: &str,
    ) -> Result<(), String> {
        let reference_seq_id = reference_seq_id.trim();
        if reference_seq_id.is_empty() {
            return Err("Reference sequence id cannot be empty".to_string());
        }
        let reference_len = if self.seq_id.as_deref() == Some(reference_seq_id) {
            self.dna.read().ok().map(|dna| dna.len()).unwrap_or(0)
        } else {
            let Some(engine) = self.engine.as_ref() else {
                return Err("No engine attached while resolving dotplot reference".to_string());
            };
            let guard = engine.read().map_err(|_| {
                "Engine lock poisoned while resolving dotplot reference".to_string()
            })?;
            guard
                .state()
                .sequences
                .get(reference_seq_id)
                .map(|dna| dna.len())
                .ok_or_else(|| format!("Reference sequence '{}' not found", reference_seq_id))?
        };
        if reference_len == 0 {
            return Err(format!(
                "Reference sequence '{}' is empty and cannot be used for dotplot comparison",
                reference_seq_id
            ));
        }
        self.dotplot_ui.reference_seq_id = reference_seq_id.to_string();
        if let Some(view) = self
            .current_dotplot_reference_splicing_view()
            .filter(|view| view.seq_id == reference_seq_id)
        {
            self.dotplot_ui.reference_span_start_0based =
                view.region_start_1based.saturating_sub(1).to_string();
            self.dotplot_ui.reference_span_end_0based = view.region_end_1based.to_string();
        } else {
            self.dotplot_ui.reference_span_start_0based = "0".to_string();
            self.dotplot_ui.reference_span_end_0based = reference_len.to_string();
        }
        self.invalidate_dotplot_cache();
        Ok(())
    }

    pub(super) fn maybe_normalize_dotplot_reference_from_current_input(&mut self) -> bool {
        let reference_seq_id = self.dotplot_ui.reference_seq_id.trim().to_string();
        if reference_seq_id.is_empty() {
            return false;
        }
        self.set_dotplot_reference_to_sequence(&reference_seq_id)
            .is_ok()
    }

    pub(super) fn use_current_query_as_dotplot_reference(&mut self) -> Result<(), String> {
        let Some(query_seq_id) = self.current_dotplot_query_seq_id() else {
            return Err("No active dotplot query is available".to_string());
        };
        self.dotplot_ui.mode = DotplotMode::PairForward;
        self.set_dotplot_reference_to_sequence(&query_seq_id)?;
        self.op_status = format!(
            "Using current query '{}' as the dotplot reference in pair-forward mode.",
            query_seq_id
        );
        Ok(())
    }

    pub(super) fn use_annotated_transcript_as_dotplot_reference(
        &mut self,
        transcript_feature_id: usize,
    ) -> Result<(), String> {
        let Some(view) = self.current_dotplot_reference_splicing_view() else {
            return Err(
                "Annotated transcript references are only available when a splicing or RNA-read mapping locus is active."
                    .to_string(),
            );
        };
        let lane = view
            .transcripts
            .iter()
            .find(|lane| lane.transcript_feature_id == transcript_feature_id)
            .ok_or_else(|| {
                format!(
                    "Annotated transcript n-{} is not available in the current locus",
                    transcript_feature_id + 1
                )
            })?;
        let lane_label = if lane.label.trim().is_empty() {
            lane.transcript_id.trim()
        } else {
            lane.label.trim()
        };
        let derived_seq_id = self
            .derive_transcript_sequence_for_dotplot_reference(
                view.seq_id.clone(),
                transcript_feature_id,
                &format!(
                    "annotated transcript n-{} '{}'",
                    transcript_feature_id + 1,
                    lane_label
                ),
            )
            .ok_or_else(|| {
                format!(
                    "Could not derive annotated transcript n-{} for dotplot reference use",
                    transcript_feature_id + 1
                )
            })?;
        self.dotplot_ui.mode = DotplotMode::PairForward;
        self.set_dotplot_reference_to_sequence(&derived_seq_id)?;
        self.op_status = format!(
            "{}\nUsing '{}' as the current dotplot reference in pair-forward mode.",
            self.op_status, derived_seq_id
        );
        Ok(())
    }

    pub(super) fn render_dotplot_reference_quick_actions(
        &mut self,
        ui: &mut egui::Ui,
        id_namespace: &str,
        save_state: &mut bool,
    ) {
        if ui
            .small_button("Use query as ref")
            .on_hover_text(
                "Set the current query sequence as the reference too and reset the reference span to the full query length. This is useful for self-vs-self density checks and should give a clean diagonal in pair-forward mode.",
            )
            .clicked()
        {
            match self.use_current_query_as_dotplot_reference() {
                Ok(()) => *save_state = true,
                Err(message) => self.op_status = message,
            }
        }
        if let Some(view) = self.current_dotplot_reference_splicing_view() {
            let locus_seq_id = view.seq_id.clone();
            let locus_label = view.group_label.clone();
            if ui
                .small_button("Use locus DNA")
                .on_hover_text(
                    "Restore the current locus DNA sequence as the pairwise reference and reset ref_start/ref_end to the current splicing/RNA-mapping ROI.",
                )
                .clicked()
            {
                self.dotplot_ui.mode = DotplotMode::PairForward;
                match self.set_dotplot_reference_to_sequence(&locus_seq_id) {
                    Ok(()) => {
                        self.op_status = format!(
                            "Using locus DNA '{}' for {} as the dotplot reference in pair-forward mode.",
                            locus_seq_id, locus_label
                        );
                        *save_state = true;
                    }
                    Err(message) => self.op_status = message,
                }
            }
        }
        let transcript_choices = self.dotplot_transcript_reference_choices();
        if !transcript_choices.is_empty() {
            let mut selected_feature_id: Option<usize> = None;
            let transcript_combo = egui::ComboBox::from_id_salt(format!(
                "dotplot_transcript_reference_{id_namespace}_{}",
                self.dotplot_window_identity_seq_id()
                    .unwrap_or_else(|| "<none>".to_string())
            ))
            .selected_text("Annotated mRNA ref...")
            .show_ui(ui, |ui| {
                for (feature_id, label) in &transcript_choices {
                    if ui.selectable_label(false, label).clicked() {
                        selected_feature_id = Some(*feature_id);
                    }
                }
            });
            transcript_combo.response.on_hover_text(
                "Derive one annotated transcript/mRNA sequence for the current locus and use it as the dotplot reference in pair-forward mode.",
            );
            if let Some(feature_id) = selected_feature_id {
                match self.use_annotated_transcript_as_dotplot_reference(feature_id) {
                    Ok(()) => *save_state = true,
                    Err(message) => self.op_status = message,
                }
            }
        }
    }

    pub(super) fn comparable_pair_forward_point_count(&self, view: &DotplotView) -> Option<usize> {
        if !matches!(view.mode, DotplotMode::PairReverseComplement) {
            return None;
        }
        let Some(engine) = self.engine.as_ref() else {
            return None;
        };
        let Ok(guard) = engine.try_read() else {
            return None;
        };
        guard
            .list_dotplot_views(Some(view.seq_id.as_str()))
            .into_iter()
            .filter(|row| row.mode == DotplotMode::PairForward)
            .filter(|row| row.reference_seq_id == view.reference_seq_id)
            .filter(|row| row.span_start_0based == view.span_start_0based)
            .filter(|row| row.span_end_0based == view.span_end_0based)
            .filter(|row| row.reference_span_start_0based == view.reference_span_start_0based)
            .filter(|row| row.reference_span_end_0based == view.reference_span_end_0based)
            .filter(|row| row.word_size == view.word_size)
            .filter(|row| row.step_bp == view.step_bp)
            .filter(|row| row.max_mismatches == view.max_mismatches)
            .max_by_key(|row| (row.generated_at_unix_ms, row.dotplot_id.clone()))
            .map(|row| row.point_count)
    }

    pub(super) fn dotplot_mode_requires_reference(mode: DotplotMode) -> bool {
        matches!(
            mode,
            DotplotMode::PairForward | DotplotMode::PairReverseComplement
        )
    }

    pub(super) fn add_small_uint_text_edit(
        ui: &mut egui::Ui,
        value: &mut String,
        expected_digits: usize,
    ) -> egui::Response {
        let expected_digits = expected_digits.clamp(1, 8);
        let width = 20.0 + expected_digits as f32 * 8.0;
        ui.add_sized(
            [width, ui.spacing().interact_size.y],
            egui::TextEdit::singleline(value),
        )
    }

    pub(super) fn dotplot_window_count(span_bp: usize, word_size: usize, step_bp: usize) -> usize {
        if word_size == 0 || step_bp == 0 || span_bp < word_size {
            0
        } else {
            span_bp
                .saturating_sub(word_size)
                .checked_div(step_bp)
                .unwrap_or(0)
                .saturating_add(1)
        }
    }

    pub(super) fn dotplot_reference_hit_envelope(
        view: &DotplotView,
        padding_bp: usize,
    ) -> Option<(usize, usize)> {
        if view.reference_span_end_0based <= view.reference_span_start_0based {
            return None;
        }
        let mut min_hit = usize::MAX;
        let mut max_hit = 0usize;
        if Self::dotplot_view_is_overlay(view) {
            for series in &view.query_series {
                for point in &series.points {
                    min_hit = min_hit.min(point.y_0based);
                    max_hit = max_hit.max(point.y_0based);
                }
            }
        } else {
            for point in &view.points {
                min_hit = min_hit.min(point.y_0based);
                max_hit = max_hit.max(point.y_0based);
            }
        }
        if min_hit == usize::MAX {
            return None;
        }
        let ref_start = view.reference_span_start_0based;
        let ref_end = view.reference_span_end_0based;
        let padded_start = min_hit.saturating_sub(padding_bp).max(ref_start);
        let padded_end = max_hit
            .saturating_add(1)
            .saturating_add(padding_bp)
            .min(ref_end);
        if padded_end <= padded_start {
            return None;
        }
        Some((padded_start, padded_end))
    }

    pub(super) fn fit_dotplot_reference_span_to_loaded_hits(&mut self) -> bool {
        self.ensure_dotplot_cache_current();
        let Some(view) = self.dotplot_cached_view.clone() else {
            self.op_status = "No loaded dotplot payload to fit.".to_string();
            return false;
        };
        if !Self::dotplot_mode_requires_reference(view.mode) {
            self.op_status =
                "Reference span fitting is only available in pair dotplot modes.".to_string();
            return false;
        }
        let Some((fit_start, fit_end)) =
            Self::dotplot_reference_hit_envelope(&view, DOTPLOT_HIT_ENVELOPE_PADDING_BP)
        else {
            self.op_status = format!("Dotplot '{}' has no hit envelope to fit.", view.dotplot_id);
            return false;
        };
        self.dotplot_ui.reference_span_start_0based = fit_start.to_string();
        self.dotplot_ui.reference_span_end_0based = fit_end.to_string();
        self.dotplot_last_compute_status = format!(
            "Fitted reference span to hit envelope: [{}..{}] (padding={} bp). Recompute to refresh payload.",
            fit_start.saturating_add(1),
            fit_end,
            DOTPLOT_HIT_ENVELOPE_PADDING_BP
        );
        true
    }

    pub(super) fn auto_fit_reference_span_for_view_if_requested(
        auto_fit_requested: bool,
        view: &DotplotView,
        padding_bp: usize,
    ) -> Option<(usize, usize)> {
        if !auto_fit_requested || !Self::dotplot_mode_requires_reference(view.mode) {
            return None;
        }
        let (fit_start, fit_end) = Self::dotplot_reference_hit_envelope(view, padding_bp)?;
        if fit_start == view.reference_span_start_0based
            && fit_end == view.reference_span_end_0based
        {
            return None;
        }
        Some((fit_start, fit_end))
    }

    pub(super) fn dotplot_quantile(sorted_values: &[f32], q: f32) -> f32 {
        if sorted_values.is_empty() {
            return 0.0;
        }
        let q = q.clamp(0.0, 1.0);
        let idx = ((sorted_values.len().saturating_sub(1)) as f32 * q).round() as usize;
        sorted_values
            .get(idx.min(sorted_values.len().saturating_sub(1)))
            .copied()
            .unwrap_or(0.0)
    }

    pub(super) fn recommend_dotplot_display_from_view(
        view: &DotplotView,
        overlay_x_axis_mode: DotplotOverlayXAxisMode,
    ) -> Option<(f32, f32, String)> {
        let overlay_mode = Self::dotplot_view_is_overlay(view);
        let total_points = Self::dotplot_view_total_point_count(view);
        if total_points == 0 {
            return None;
        }
        let query_span = if overlay_mode {
            let max_query_span = Self::dotplot_overlay_max_query_span(view);
            let average_query_span = view
                .query_series
                .iter()
                .map(|series| {
                    series
                        .span_end_0based
                        .saturating_sub(series.span_start_0based)
                        .max(1)
                })
                .sum::<usize>()
                .checked_div(view.query_series.len().max(1))
                .unwrap_or(1)
                .max(1);
            overlay_x_axis_mode.plot_query_span_bp(max_query_span, average_query_span)
        } else {
            view.span_end_0based
                .saturating_sub(view.span_start_0based)
                .max(1)
        };
        let reference_span = view
            .reference_span_end_0based
            .saturating_sub(view.reference_span_start_0based)
            .max(1);
        let cols: usize = 192;
        let aspect = reference_span as f32 / query_span as f32;
        let rows = ((cols as f32 * aspect).round() as usize).clamp(48, 384);
        let total_cells = cols.saturating_mul(rows).max(1);

        let query_span_max = query_span.saturating_sub(1).max(1);
        let reference_span_max = reference_span.saturating_sub(1).max(1);
        let mut cells: HashMap<(usize, usize), usize> = HashMap::new();
        if overlay_mode {
            let max_query_span = Self::dotplot_overlay_max_query_span(view);
            for series in &view.query_series {
                for point in &series.points {
                    let y_local = point
                        .y_0based
                        .saturating_sub(view.reference_span_start_0based)
                        .min(reference_span_max);
                    let x_frac = overlay_x_axis_mode.point_fraction(
                        point.x_0based,
                        series.span_start_0based,
                        series.span_end_0based,
                        max_query_span,
                    );
                    let y_frac = (y_local as f32 / reference_span_max as f32).clamp(0.0, 1.0);
                    let x_cell = ((x_frac * (cols.saturating_sub(1)) as f32).round() as usize)
                        .min(cols.saturating_sub(1));
                    let y_cell = ((y_frac * (rows.saturating_sub(1)) as f32).round() as usize)
                        .min(rows.saturating_sub(1));
                    let entry = cells.entry((x_cell, y_cell)).or_insert(0);
                    *entry = entry.saturating_add(1);
                }
            }
        } else {
            for point in &view.points {
                let x_local = point
                    .x_0based
                    .saturating_sub(view.span_start_0based)
                    .min(query_span_max);
                let y_local = point
                    .y_0based
                    .saturating_sub(view.reference_span_start_0based)
                    .min(reference_span_max);
                let x_frac = (x_local as f32 / query_span_max as f32).clamp(0.0, 1.0);
                let y_frac = (y_local as f32 / reference_span_max as f32).clamp(0.0, 1.0);
                let x_cell = ((x_frac * (cols.saturating_sub(1)) as f32).round() as usize)
                    .min(cols.saturating_sub(1));
                let y_cell = ((y_frac * (rows.saturating_sub(1)) as f32).round() as usize)
                    .min(rows.saturating_sub(1));
                let entry = cells.entry((x_cell, y_cell)).or_insert(0);
                *entry = entry.saturating_add(1);
            }
        }
        if cells.is_empty() {
            return None;
        }
        let mut counts: Vec<usize> = cells.values().copied().collect();
        counts.sort_unstable();
        let max_count = (*counts.last().unwrap_or(&1)).max(1) as f32;
        let mut densities: Vec<f32> = counts
            .into_iter()
            .map(|count| (count as f32 / max_count).clamp(0.0, 1.0))
            .collect();
        densities.sort_by(|left, right| left.partial_cmp(right).unwrap_or(Ordering::Equal));
        let occupancy = cells.len() as f32 / total_cells as f32;

        let q_target = if occupancy > 0.35 {
            0.70
        } else if occupancy > 0.20 {
            0.55
        } else if occupancy > 0.08 {
            0.35
        } else {
            0.15
        };
        let mut threshold = if total_points >= 800 {
            Self::dotplot_quantile(&densities, q_target)
        } else {
            0.0
        };
        if occupancy < 0.02 {
            threshold = threshold.min(0.05);
        }
        threshold = threshold.clamp(0.0, 0.90);
        let q90 = Self::dotplot_quantile(&densities, 0.90);
        let q90_post = if q90 <= threshold {
            0.05
        } else {
            ((q90 - threshold) / (1.0 - threshold)).clamp(0.05, 1.0)
        };
        let mut gain = (0.85 / q90_post).clamp(0.25, 8.0);
        if total_points < 300 {
            gain = gain.max(1.5);
        }

        let summary = format!(
            "Auto contrast: threshold={:.3}, gain={:.3}, occupancy={:.2}%, q90={:.3}, points={}{}{}, cells={}/{}{}",
            threshold,
            gain,
            occupancy * 100.0,
            q90,
            total_points,
            if overlay_mode { ", series=" } else { "" },
            if overlay_mode {
                view.query_series.len().to_string()
            } else {
                String::new()
            },
            cells.len(),
            total_cells,
            if overlay_mode {
                format!(", overlay_x={}", overlay_x_axis_mode.as_str())
            } else {
                String::new()
            }
        );
        Some((threshold, gain, summary))
    }

    pub(super) fn apply_auto_dotplot_display_contrast(&mut self) {
        self.ensure_dotplot_cache_current();
        let Some(view) = self.dotplot_cached_view.as_ref() else {
            self.dotplot_last_compute_status =
                "Auto contrast skipped: no dotplot payload loaded".to_string();
            return;
        };
        let Some((threshold, gain, summary)) =
            Self::recommend_dotplot_display_from_view(view, self.dotplot_ui.overlay_x_axis_mode)
        else {
            self.dotplot_last_compute_status = format!(
                "Auto contrast skipped: payload '{}' has no usable point density",
                view.dotplot_id
            );
            return;
        };
        self.dotplot_ui.display_density_threshold = threshold;
        self.dotplot_ui.display_intensity_gain = gain;
        self.dotplot_last_compute_status = summary;
        self.save_engine_ops_state();
    }

    pub(super) fn resolve_dotplot_span_bounds_for_status(
        seq_len: usize,
        start_opt: Option<usize>,
        end_opt: Option<usize>,
        label: &str,
    ) -> Result<(usize, usize), String> {
        if seq_len == 0 {
            return Err(format!("Selected {label} sequence is empty"));
        }
        let start = start_opt.unwrap_or(0);
        let end = end_opt.unwrap_or(seq_len);
        if start >= seq_len {
            return Err(format!(
                "{label} span start ({start}) must be within sequence length ({seq_len})"
            ));
        }
        if end > seq_len {
            return Err(format!(
                "{label} span end ({end}) must be <= sequence length ({seq_len})"
            ));
        }
        if end <= start {
            return Err(format!(
                "{label} span end ({end}) must be > span start ({start})"
            ));
        }
        Ok((start, end))
    }

    pub(super) fn default_dotplot_half_window_bp(&self) -> Option<usize> {
        let query_len = self.current_dotplot_query_length();
        if query_len == 0 {
            return None;
        }
        let mut largest_len = query_len;
        if Self::dotplot_mode_requires_reference(self.dotplot_ui.mode) {
            let reference_seq_id = self.dotplot_ui.reference_seq_id.trim();
            if !reference_seq_id.is_empty()
                && let Some(engine) = self.engine.as_ref()
                && let Ok(guard) = engine.read()
                && let Some(reference_dna) = guard.state().sequences.get(reference_seq_id)
            {
                largest_len = largest_len.max(reference_dna.len());
            }
        }
        Some(largest_len.max(1))
    }

    pub(super) fn normalize_dotplot_half_window_default_if_needed(&mut self) -> bool {
        let raw = self.dotplot_ui.half_window_bp.trim();
        if !raw.is_empty() {
            return false;
        }
        let Some(default_half_window_bp) = self.default_dotplot_half_window_bp() else {
            return false;
        };
        let normalized = default_half_window_bp.to_string();
        if self.dotplot_ui.half_window_bp == normalized {
            return false;
        }
        self.dotplot_ui.half_window_bp = normalized;
        true
    }

    pub(super) fn resolve_dotplot_half_window_bp(&self, field_name: &str) -> Result<usize, String> {
        if self.dotplot_ui.half_window_bp.trim().is_empty() {
            return self.default_dotplot_half_window_bp().ok_or_else(|| {
                format!("Could not resolve {field_name}; active sequence is empty")
            });
        }
        Self::parse_positive_usize_text(&self.dotplot_ui.half_window_bp, field_name)
    }

    pub(super) fn build_dotplot_compute_diagnostics(
        &self,
    ) -> Result<DotplotComputeDiagnostics, String> {
        let Some(seq_id) = self.current_dotplot_query_seq_id() else {
            return Err("No active sequence selected for dotplot diagnostics".to_string());
        };
        let half_window_bp = self.resolve_dotplot_half_window_bp("dotplot half_window_bp")?;
        let word_size =
            Self::parse_positive_usize_text(&self.dotplot_ui.word_size, "dotplot word")?;
        let step_bp = Self::parse_positive_usize_text(&self.dotplot_ui.step_bp, "dotplot step")?;
        let max_mismatches =
            Self::parse_optional_usize_text(&self.dotplot_ui.max_mismatches, "dotplot mismatches")?
                .unwrap_or(0);
        let tile_bp = Self::parse_optional_usize_text(&self.dotplot_ui.tile_bp, "dotplot tile_bp")?
            .filter(|value| *value > 0);
        let (query_span_start_0based, query_span_end_0based) = self
            .default_dotplot_span_for_view(half_window_bp)
            .ok_or_else(|| "Active sequence is empty; query span unavailable".to_string())?;
        let query_span_bp = query_span_end_0based.saturating_sub(query_span_start_0based);
        let mode = self.dotplot_ui.mode;

        let (
            reference_seq_id,
            reference_span_start_0based,
            reference_span_end_0based,
            reference_span_bp,
        ) = if Self::dotplot_mode_requires_reference(mode) {
            let ref_id = self.dotplot_ui.reference_seq_id.trim();
            if ref_id.is_empty() {
                return Err("Pair mode requires reference_seq_id".to_string());
            }
            let Some(engine) = self.engine.as_ref() else {
                return Err("No engine attached while resolving reference span".to_string());
            };
            let guard = engine
                .read()
                .map_err(|_| "Engine lock poisoned while resolving reference span".to_string())?;
            let ref_len = guard
                .state()
                .sequences
                .get(ref_id)
                .map(|dna| dna.len())
                .ok_or_else(|| format!("Reference sequence '{ref_id}' not found"))?;
            let ref_start_opt = Self::parse_optional_usize_text(
                &self.dotplot_ui.reference_span_start_0based,
                "dotplot ref_start",
            )?;
            let ref_end_opt = Self::parse_optional_usize_text(
                &self.dotplot_ui.reference_span_end_0based,
                "dotplot ref_end",
            )?;
            let (ref_start, ref_end) = Self::resolve_dotplot_span_bounds_for_status(
                ref_len,
                ref_start_opt,
                ref_end_opt,
                "Reference",
            )?;
            (
                Some(ref_id.to_string()),
                ref_start,
                ref_end,
                ref_end.saturating_sub(ref_start),
            )
        } else {
            (
                None,
                query_span_start_0based,
                query_span_end_0based,
                query_span_bp,
            )
        };

        let query_windows = Self::dotplot_window_count(query_span_bp, word_size, step_bp);
        let reference_windows = Self::dotplot_window_count(reference_span_bp, word_size, step_bp);
        let estimated_pair_evaluations = query_windows.saturating_mul(reference_windows);
        Ok(DotplotComputeDiagnostics {
            seq_id,
            reference_seq_id,
            mode,
            query_span_start_0based,
            query_span_end_0based,
            query_span_bp,
            reference_span_start_0based,
            reference_span_end_0based,
            reference_span_bp,
            word_size,
            step_bp,
            max_mismatches,
            tile_bp,
            query_windows,
            reference_windows,
            estimated_pair_evaluations,
        })
    }

    pub(super) fn build_dotplot_overlay_estimated_pair_evaluations(&self) -> Result<usize, String> {
        if !self.dotplot_ui.overlay_enabled {
            return Err("Dotplot overlay mode is disabled".to_string());
        }
        let choices = self.dotplot_overlay_choices();
        if choices.is_empty() {
            return Err(
                "Overlay queries are only available for the active locus reference".to_string(),
            );
        }
        let selected_feature_ids = self
            .dotplot_ui
            .overlay_transcript_feature_ids
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        if selected_feature_ids.is_empty() {
            return Err("Select at least one transcript for dotplot overlay".to_string());
        }
        let word_size =
            Self::parse_positive_usize_text(&self.dotplot_ui.word_size, "dotplot word")?;
        let step_bp = Self::parse_positive_usize_text(&self.dotplot_ui.step_bp, "dotplot step")?;
        let reference_seq_id = self.dotplot_ui.reference_seq_id.trim();
        if reference_seq_id.is_empty() {
            return Err("Overlay mode requires reference_seq_id".to_string());
        }
        let Some(engine) = self.engine.as_ref() else {
            return Err("No engine attached while resolving overlay diagnostics".to_string());
        };
        let guard = engine
            .read()
            .map_err(|_| "Engine lock poisoned while resolving overlay diagnostics".to_string())?;
        let reference_len = guard
            .state()
            .sequences
            .get(reference_seq_id)
            .map(|dna| dna.len())
            .ok_or_else(|| format!("Reference sequence '{reference_seq_id}' not found"))?;
        let reference_span_start_0based = Self::parse_optional_usize_text(
            &self.dotplot_ui.reference_span_start_0based,
            "dotplot ref_start",
        )?;
        let reference_span_end_0based = Self::parse_optional_usize_text(
            &self.dotplot_ui.reference_span_end_0based,
            "dotplot ref_end",
        )?;
        let (reference_start, reference_end) = Self::resolve_dotplot_span_bounds_for_status(
            reference_len,
            reference_span_start_0based,
            reference_span_end_0based,
            "Reference",
        )?;
        let reference_span_bp = reference_end.saturating_sub(reference_start);
        let reference_windows =
            Self::dotplot_window_count(reference_span_bp, word_size.max(1), step_bp.max(1));
        let mut estimated_pair_evaluations = 0usize;
        for choice in choices
            .iter()
            .filter(|choice| selected_feature_ids.contains(&choice.transcript_feature_id))
        {
            let query_len = guard
                .state()
                .sequences
                .get(&choice.derived_seq_id)
                .map(|dna| dna.len())
                .unwrap_or(choice.estimated_length_bp);
            if query_len < word_size {
                return Err(format!(
                    "Overlay query '{}' is shorter than word size {}",
                    choice.label, word_size
                ));
            }
            let query_windows =
                Self::dotplot_window_count(query_len, word_size.max(1), step_bp.max(1));
            estimated_pair_evaluations = estimated_pair_evaluations
                .saturating_add(query_windows.saturating_mul(reference_windows));
        }
        Ok(estimated_pair_evaluations)
    }

    pub(super) fn resolve_dotplot_overlay_query_specs_for_compute(
        &mut self,
    ) -> Result<Vec<DotplotOverlayQuerySpec>, String> {
        if !self.dotplot_ui.overlay_enabled {
            return Err("Dotplot overlay mode is disabled".to_string());
        }
        let Some(reference_view) = self.dotplot_overlay_reference_view() else {
            return Err(
                "Overlay queries are only available when the active locus DNA is the reference"
                    .to_string(),
            );
        };
        let choices = self.dotplot_overlay_choices();
        let selected_feature_ids = self
            .dotplot_ui
            .overlay_transcript_feature_ids
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        let selected_choices = choices
            .into_iter()
            .filter(|choice| selected_feature_ids.contains(&choice.transcript_feature_id))
            .collect::<Vec<_>>();
        if selected_choices.is_empty() {
            return Err("Select at least one transcript for dotplot overlay".to_string());
        }
        let Some(engine) = self.engine.as_ref() else {
            return Err("No engine attached while resolving dotplot overlay".to_string());
        };
        let mut resolved_seq_ids: BTreeMap<usize, String> = BTreeMap::new();
        let mut missing_feature_ids: Vec<usize> = vec![];
        {
            let guard = engine
                .read()
                .map_err(|_| "Engine lock poisoned while resolving dotplot overlay".to_string())?;
            for choice in &selected_choices {
                if guard.state().sequences.contains_key(&choice.derived_seq_id) {
                    resolved_seq_ids
                        .insert(choice.transcript_feature_id, choice.derived_seq_id.clone());
                } else {
                    missing_feature_ids.push(choice.transcript_feature_id);
                }
            }
        }
        if !missing_feature_ids.is_empty() {
            let result = self
                .derive_transcript_sequences_with_compact_feedback(
                    reference_view.seq_id.clone(),
                    missing_feature_ids.clone(),
                    None,
                    "selected dotplot overlay transcripts",
                )
                .ok_or_else(|| {
                    if self.op_status.trim().is_empty() {
                        "Could not derive overlay transcript sequences".to_string()
                    } else {
                        self.op_status.clone()
                    }
                })?;
            if result.created_seq_ids.len() != missing_feature_ids.len() {
                return Err(format!(
                    "Overlay transcript derivation created {} sequence(s) for {} requested transcript(s)",
                    result.created_seq_ids.len(),
                    missing_feature_ids.len()
                ));
            }
            for (feature_id, seq_id) in missing_feature_ids
                .iter()
                .copied()
                .zip(result.created_seq_ids.iter().cloned())
            {
                resolved_seq_ids.insert(feature_id, seq_id);
            }
        }
        selected_choices
            .into_iter()
            .map(|choice| {
                let seq_id = resolved_seq_ids
                    .get(&choice.transcript_feature_id)
                    .cloned()
                    .ok_or_else(|| {
                        format!(
                            "No derived sequence id was resolved for transcript '{}'",
                            choice.label
                        )
                    })?;
                Ok(DotplotOverlayQuerySpec {
                    seq_id,
                    label: choice.label,
                    transcript_feature_id: Some(choice.transcript_feature_id),
                    query_anchor_0based: None,
                    query_anchor_label: None,
                    span_start_0based: None,
                    span_end_0based: None,
                    mode: choice.mode,
                    color_rgb: Some(choice.color_rgb),
                })
            })
            .collect()
    }

    pub(super) fn maybe_auto_compute_dotplot(&mut self) {
        if !self.dotplot_auto_compute_dirty {
            return;
        }
        let Some(due_at) = self.dotplot_auto_compute_due_at else {
            return;
        };
        if Instant::now() < due_at {
            return;
        }
        let auto_compute_result = if self.dotplot_ui.overlay_enabled {
            self.build_dotplot_overlay_estimated_pair_evaluations()
                .map(|evals| ("overlay".to_string(), evals))
        } else {
            self.build_dotplot_compute_diagnostics().map(|diag| {
                (
                    diag.mode.as_str().to_string(),
                    diag.estimated_pair_evaluations,
                )
            })
        };
        match auto_compute_result {
            Ok((_label, estimated_pair_evaluations))
                if estimated_pair_evaluations <= MAX_DOTPLOT_PAIR_EVALUATIONS =>
            {
                self.compute_primary_dotplot();
            }
            Ok((label, estimated_pair_evaluations)) => {
                self.dotplot_last_compute_status = format!(
                    "Parameters changed, but auto-compute is paused for {label} because the estimate is {} pair evaluations (> {}). Press `Compute dotplot`.",
                    estimated_pair_evaluations, MAX_DOTPLOT_PAIR_EVALUATIONS
                );
                self.dotplot_auto_compute_due_at = None;
            }
            Err(message) => {
                self.dotplot_last_compute_status = format!(
                    "Parameters changed. Auto-compute is waiting for valid inputs: {message}"
                );
                self.dotplot_auto_compute_due_at = None;
            }
        }
    }

    pub(super) fn bounded_center_window(
        sequence_len: usize,
        center_0based: usize,
        half_window_bp: usize,
    ) -> Option<(usize, usize)> {
        if sequence_len == 0 {
            return None;
        }
        let half_window_bp = half_window_bp.max(1);
        let center_0based = center_0based.min(sequence_len.saturating_sub(1));
        let target_span = half_window_bp
            .saturating_mul(2)
            .saturating_add(1)
            .min(sequence_len);
        let mut start = center_0based.saturating_sub(half_window_bp);
        let mut end = center_0based
            .saturating_add(half_window_bp)
            .saturating_add(1)
            .min(sequence_len);
        let current_span = end.saturating_sub(start);
        if current_span < target_span {
            let deficit = target_span - current_span;
            let shift_left = deficit.min(start);
            start = start.saturating_sub(shift_left);
            let remaining = deficit.saturating_sub(shift_left);
            end = end.saturating_add(remaining).min(sequence_len);
            let second_span = end.saturating_sub(start);
            if second_span < target_span {
                let second_deficit = target_span - second_span;
                start = start.saturating_sub(second_deficit.min(start));
            }
        }
        if end <= start {
            end = (start + 1).min(sequence_len);
        }
        Some((start, end))
    }

    pub(super) fn default_dotplot_span_for_view(
        &self,
        half_window_bp: usize,
    ) -> Option<(usize, usize)> {
        let sequence_len = self.current_dotplot_query_length();
        if sequence_len == 0 {
            return None;
        }
        let center_0based = if self.using_dotplot_query_override() {
            sequence_len / 2
        } else if self.is_circular() {
            sequence_len / 2
        } else {
            let (start_bp, span_bp, _) = self.current_linear_viewport();
            if span_bp == 0 {
                sequence_len / 2
            } else {
                start_bp
                    .saturating_add(span_bp / 2)
                    .min(sequence_len.saturating_sub(1))
            }
        };
        Self::bounded_center_window(sequence_len, center_0based, half_window_bp)
    }

    pub(super) fn ensure_dotplot_cache_current(&mut self) {
        let Some(seq_id) = self.current_dotplot_owner_seq_id() else {
            self.invalidate_dotplot_cache();
            return;
        };
        let Some(engine) = self.engine.as_ref() else {
            self.invalidate_dotplot_cache();
            return;
        };
        let mut selected_view_id = self.dotplot_ui.dotplot_id.trim().to_string();
        let mut selected_track_id = self.dotplot_ui.flex_track_id.trim().to_string();
        let Ok(guard) = engine.read() else {
            return;
        };
        if selected_view_id.is_empty()
            && let Some(row) = guard.list_dotplot_views(Some(seq_id.as_str())).last()
        {
            selected_view_id = row.dotplot_id.clone();
        }
        if selected_track_id.is_empty()
            && let Some(row) = guard.list_flexibility_tracks(Some(seq_id.as_str())).last()
        {
            selected_track_id = row.track_id.clone();
        }
        if self.dotplot_cache_seq_id == seq_id
            && self.dotplot_cache_view_id == selected_view_id
            && self.dotplot_cache_track_id == selected_track_id
        {
            return;
        }
        let view = if selected_view_id.is_empty() {
            None
        } else {
            guard
                .get_dotplot_view(selected_view_id.as_str())
                .ok()
                .filter(|row| row.owner_seq_id == seq_id)
        };
        let track = if selected_track_id.is_empty() {
            None
        } else {
            guard
                .get_flexibility_track(selected_track_id.as_str())
                .ok()
                .filter(|row| row.seq_id == seq_id)
        };
        drop(guard);

        let mut changed = false;
        if self.dotplot_ui.dotplot_id.trim().is_empty() && !selected_view_id.is_empty() {
            self.dotplot_ui.dotplot_id = selected_view_id.clone();
            changed = true;
        }
        if self.dotplot_ui.flex_track_id.trim().is_empty() && !selected_track_id.is_empty() {
            self.dotplot_ui.flex_track_id = selected_track_id.clone();
            changed = true;
        }

        self.dotplot_cached_view = view;
        self.dotplot_cached_flex_track = track;
        self.dotplot_cache_seq_id = seq_id;
        self.dotplot_cache_view_id = selected_view_id;
        self.dotplot_cache_track_id = selected_track_id;
        if changed {
            self.save_engine_ops_state();
        }
    }

    pub(super) fn compute_primary_dotplot(&mut self) {
        let requires_reference = Self::dotplot_mode_requires_reference(self.dotplot_ui.mode);
        let overlay_enabled = self.dotplot_ui.overlay_enabled && requires_reference;
        let word_size = match Self::parse_positive_usize_text(
            &self.dotplot_ui.word_size,
            "dotplot word_size",
        ) {
            Ok(v) => v,
            Err(e) => {
                self.op_status = e;
                return;
            }
        };
        let step_bp =
            match Self::parse_positive_usize_text(&self.dotplot_ui.step_bp, "dotplot step_bp") {
                Ok(v) => v,
                Err(e) => {
                    self.op_status = e;
                    return;
                }
            };
        let max_mismatches = match Self::parse_optional_usize_text(
            &self.dotplot_ui.max_mismatches,
            "dotplot max_mismatches",
        ) {
            Ok(Some(v)) => v,
            Ok(None) => 0,
            Err(e) => {
                self.op_status = e;
                return;
            }
        };
        let tile_bp =
            match Self::parse_optional_usize_text(&self.dotplot_ui.tile_bp, "dotplot tile_bp") {
                Ok(v) => v.filter(|value| *value > 0),
                Err(e) => {
                    self.op_status = e;
                    return;
                }
            };
        if requires_reference {
            self.maybe_normalize_dotplot_reference_from_current_input();
        }
        let query_seq_id = self.current_dotplot_query_seq_id();
        if !overlay_enabled && query_seq_id.is_none() {
            self.op_status = "No active sequence selected for dotplot computation".to_string();
            return;
        }
        let reference_seq_id = if requires_reference {
            let value = self.dotplot_ui.reference_seq_id.trim();
            if value.is_empty() {
                self.op_status = format!(
                    "Dotplot mode '{}' requires a reference sequence id",
                    self.dotplot_ui.mode.as_str()
                );
                return;
            }
            Some(value.to_string())
        } else {
            None
        };
        let reference_span_start_0based = if requires_reference {
            match Self::parse_optional_usize_text(
                &self.dotplot_ui.reference_span_start_0based,
                "dotplot reference_span_start_0based",
            ) {
                Ok(value) => value,
                Err(e) => {
                    self.op_status = e;
                    return;
                }
            }
        } else {
            None
        };
        let reference_span_end_0based = if requires_reference {
            match Self::parse_optional_usize_text(
                &self.dotplot_ui.reference_span_end_0based,
                "dotplot reference_span_end_0based",
            ) {
                Ok(value) => value,
                Err(e) => {
                    self.op_status = e;
                    return;
                }
            }
        } else {
            None
        };
        let auto_fit_reference_span_requested = requires_reference
            && reference_span_start_0based.is_none()
            && reference_span_end_0based.is_none();
        let store_as = if self.dotplot_ui.dotplot_id.trim().is_empty() {
            if overlay_enabled {
                let owner = reference_seq_id
                    .as_deref()
                    .unwrap_or_else(|| self.seq_id.as_deref().unwrap_or("dotplot"));
                format!("{owner}.dotplot.overlay")
            } else {
                format!("{}.dotplot", query_seq_id.as_deref().unwrap_or("dotplot"))
            }
        } else {
            self.dotplot_ui.dotplot_id.trim().to_string()
        };
        self.dotplot_ui.dotplot_id = store_as.clone();
        if overlay_enabled {
            let estimated_pair_evaluations =
                self.build_dotplot_overlay_estimated_pair_evaluations().ok();
            if let Some(reference_id) = reference_seq_id.as_deref() {
                let selected_count = self.dotplot_ui.overlay_transcript_feature_ids.len();
                self.dotplot_last_compute_status = format!(
                    "Requested overlay compute: reference={} series={} word={} step={} mismatches={} tile_bp={} pair_evals≈{}",
                    reference_id,
                    selected_count,
                    word_size,
                    step_bp,
                    max_mismatches,
                    tile_bp
                        .map(|value| value.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    estimated_pair_evaluations
                        .map(|value| value.to_string())
                        .unwrap_or_else(|| "?".to_string())
                );
            }
            let queries = match self.resolve_dotplot_overlay_query_specs_for_compute() {
                Ok(queries) => queries,
                Err(message) => {
                    self.op_status = message;
                    self.clear_dotplot_auto_compute();
                    return;
                }
            };
            let owner_seq_id = match reference_seq_id.clone() {
                Some(value) => value,
                None => {
                    self.op_status = "Overlay mode requires reference_seq_id".to_string();
                    self.clear_dotplot_auto_compute();
                    return;
                }
            };
            self.apply_operation_with_feedback(Operation::ComputeDotplotOverlay {
                owner_seq_id,
                reference_seq_id: reference_seq_id.clone().unwrap_or_default(),
                reference_span_start_0based,
                reference_span_end_0based,
                queries,
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                store_as: Some(store_as.clone()),
            });
            self.invalidate_dotplot_cache();
            self.ensure_dotplot_cache_current();
            if let Some(view) = self.dotplot_cached_view.as_ref() {
                let total_point_count = view
                    .query_series
                    .iter()
                    .map(|series| series.point_count)
                    .sum::<usize>();
                let series_labels = view
                    .query_series
                    .iter()
                    .take(4)
                    .map(|series| series.label.clone())
                    .collect::<Vec<_>>();
                self.dotplot_last_compute_status = format!(
                    "Computed overlay dotplot '{}' owner={} reference={} [{}..{}] series={} total_points={} queries={}",
                    view.dotplot_id,
                    view.owner_seq_id,
                    view.reference_seq_id
                        .as_deref()
                        .unwrap_or(view.owner_seq_id.as_str()),
                    view.reference_span_start_0based.saturating_add(1),
                    view.reference_span_end_0based,
                    view.series_count.max(view.query_series.len()),
                    total_point_count,
                    series_labels.join(", ")
                );
            }
            self.clear_dotplot_auto_compute();
            self.save_engine_ops_state();
            return;
        }

        let Some(seq_id) = query_seq_id else {
            self.op_status = "No active sequence selected for dotplot computation".to_string();
            self.clear_dotplot_auto_compute();
            return;
        };
        let mut diagnostics_snapshot = self.build_dotplot_compute_diagnostics().ok();
        let half_window_bp = match self.resolve_dotplot_half_window_bp("dotplot half_window_bp") {
            Ok(v) => v,
            Err(e) => {
                self.op_status = e;
                self.clear_dotplot_auto_compute();
                return;
            }
        };
        let Some((span_start_0based, span_end_0based)) =
            self.default_dotplot_span_for_view(half_window_bp)
        else {
            self.op_status = "Active sequence is empty; dotplot span unavailable".to_string();
            self.clear_dotplot_auto_compute();
            return;
        };
        if span_end_0based.saturating_sub(span_start_0based) < word_size {
            self.op_status = format!(
                "Dotplot span {}..{} is smaller than word_size {}",
                span_start_0based, span_end_0based, word_size
            );
            self.clear_dotplot_auto_compute();
            return;
        }
        if let Some(diag) = diagnostics_snapshot.as_ref() {
            let reference_label = diag
                .reference_seq_id
                .as_deref()
                .unwrap_or(diag.seq_id.as_str());
            self.dotplot_last_compute_status = format!(
                "Requested compute: mode={} query={} [{}..{}] reference={} [{}..{}] word={} step={} mismatches={} tile_bp={} windows(query={},ref={}) pair_evals≈{}",
                diag.mode.as_str(),
                diag.seq_id,
                diag.query_span_start_0based.saturating_add(1),
                diag.query_span_end_0based,
                reference_label,
                diag.reference_span_start_0based.saturating_add(1),
                diag.reference_span_end_0based,
                diag.word_size,
                diag.step_bp,
                diag.max_mismatches,
                diag.tile_bp
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                diag.query_windows,
                diag.reference_windows,
                diag.estimated_pair_evaluations
            );
        } else {
            self.dotplot_last_compute_status.clear();
        }
        self.apply_operation_with_feedback(Operation::ComputeDotplot {
            seq_id: seq_id.clone(),
            reference_seq_id: reference_seq_id.clone(),
            span_start_0based: Some(span_start_0based),
            span_end_0based: Some(span_end_0based),
            reference_span_start_0based,
            reference_span_end_0based,
            mode: self.dotplot_ui.mode,
            word_size,
            step_bp,
            max_mismatches,
            tile_bp,
            store_as: Some(store_as.clone()),
        });
        self.invalidate_dotplot_cache();
        self.ensure_dotplot_cache_current();
        let mut auto_fit_applied: Option<(usize, usize)> = None;
        if let Some(view) = self.dotplot_cached_view.clone()
            && let Some((fit_start, fit_end)) = Self::auto_fit_reference_span_for_view_if_requested(
                auto_fit_reference_span_requested,
                &view,
                DOTPLOT_HIT_ENVELOPE_PADDING_BP,
            )
        {
            self.dotplot_ui.reference_span_start_0based = fit_start.to_string();
            self.dotplot_ui.reference_span_end_0based = fit_end.to_string();
            self.apply_operation_with_feedback(Operation::ComputeDotplot {
                seq_id,
                reference_seq_id,
                span_start_0based: Some(span_start_0based),
                span_end_0based: Some(span_end_0based),
                reference_span_start_0based: Some(fit_start),
                reference_span_end_0based: Some(fit_end),
                mode: self.dotplot_ui.mode,
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                store_as: Some(store_as),
            });
            self.invalidate_dotplot_cache();
            self.ensure_dotplot_cache_current();
            auto_fit_applied = Some((fit_start, fit_end));
        }
        if let Some(view) = self.dotplot_cached_view.as_ref() {
            let reference_seq_label = view
                .reference_seq_id
                .as_deref()
                .unwrap_or(view.seq_id.as_str());
            let query_span_bp = view.span_end_0based.saturating_sub(view.span_start_0based);
            let reference_span_bp = view
                .reference_span_end_0based
                .saturating_sub(view.reference_span_start_0based);
            let query_windows = Self::dotplot_window_count(
                query_span_bp,
                view.word_size.max(1),
                view.step_bp.max(1),
            );
            let reference_windows = Self::dotplot_window_count(
                reference_span_bp,
                view.word_size.max(1),
                view.step_bp.max(1),
            );
            let estimated_pair_evaluations = query_windows.saturating_mul(reference_windows);
            let estimated_hit_fraction = if estimated_pair_evaluations == 0 {
                0.0
            } else {
                (view.point_count as f64 / estimated_pair_evaluations as f64) * 100.0
            };
            self.dotplot_last_compute_status = format!(
                "Computed dotplot '{}' mode={} query={} [{}..{}] reference={} [{}..{}] points={} windows(query={},ref={}) pair_evals≈{} hit_fraction≈{:.6}%",
                view.dotplot_id,
                view.mode.as_str(),
                view.seq_id,
                view.span_start_0based.saturating_add(1),
                view.span_end_0based,
                reference_seq_label,
                view.reference_span_start_0based.saturating_add(1),
                view.reference_span_end_0based,
                view.point_count,
                query_windows,
                reference_windows,
                estimated_pair_evaluations,
                estimated_hit_fraction
            );
            if let Some((fit_start, fit_end)) = auto_fit_applied {
                self.dotplot_last_compute_status = format!(
                    "Auto-fit reference span to hit envelope [{}..{}] (padding={} bp), then recomputed. {}",
                    fit_start.saturating_add(1),
                    fit_end,
                    DOTPLOT_HIT_ENVELOPE_PADDING_BP,
                    self.dotplot_last_compute_status
                );
            }
        } else if diagnostics_snapshot.is_none() {
            diagnostics_snapshot = self.build_dotplot_compute_diagnostics().ok();
        }
        if self.dotplot_cached_view.is_none()
            && let Some(diag) = diagnostics_snapshot
        {
            let reference_label = diag
                .reference_seq_id
                .as_deref()
                .unwrap_or(diag.seq_id.as_str());
            self.dotplot_last_compute_status = format!(
                "Compute finished but no cached payload was loaded for dotplot_id='{}'. Last request: mode={} query={} [{}..{}] reference={} [{}..{}] pair_evals≈{}",
                self.dotplot_ui.dotplot_id.trim(),
                diag.mode.as_str(),
                diag.seq_id,
                diag.query_span_start_0based.saturating_add(1),
                diag.query_span_end_0based,
                reference_label,
                diag.reference_span_start_0based.saturating_add(1),
                diag.reference_span_end_0based,
                diag.estimated_pair_evaluations
            );
        }
        self.clear_dotplot_auto_compute();
        self.save_engine_ops_state();
    }

    pub(super) fn compute_primary_flexibility_track(&mut self) {
        let Some(seq_id) = self.seq_id.clone() else {
            self.op_status = "No active sequence selected for flexibility computation".to_string();
            return;
        };
        let half_window_bp = match self.resolve_dotplot_half_window_bp("flexibility half_window_bp")
        {
            Ok(v) => v,
            Err(e) => {
                self.op_status = e;
                return;
            }
        };
        let bin_bp = match Self::parse_positive_usize_text(
            &self.dotplot_ui.flex_bin_bp,
            "flexibility bin_bp",
        ) {
            Ok(v) => v,
            Err(e) => {
                self.op_status = e;
                return;
            }
        };
        let smoothing_bp = match Self::parse_optional_usize_text(
            &self.dotplot_ui.flex_smoothing_bp,
            "flexibility smoothing_bp",
        ) {
            Ok(v) => v.filter(|value| *value > 0),
            Err(e) => {
                self.op_status = e;
                return;
            }
        };
        let Some((span_start_0based, span_end_0based)) =
            self.default_dotplot_span_for_view(half_window_bp)
        else {
            self.op_status = "Active sequence is empty; flexibility span unavailable".to_string();
            return;
        };
        let store_as = if self.dotplot_ui.flex_track_id.trim().is_empty() {
            format!("{seq_id}.flex")
        } else {
            self.dotplot_ui.flex_track_id.trim().to_string()
        };
        self.dotplot_ui.flex_track_id = store_as.clone();
        self.apply_operation_with_feedback(Operation::ComputeFlexibilityTrack {
            seq_id,
            span_start_0based: Some(span_start_0based),
            span_end_0based: Some(span_end_0based),
            model: self.dotplot_ui.flex_model,
            bin_bp,
            smoothing_bp,
            store_as: Some(store_as),
        });
        self.invalidate_dotplot_cache();
        self.ensure_dotplot_cache_current();
        self.save_engine_ops_state();
    }

    pub(super) fn dotplot_crosshair_selection_bounds(
        x_0based: usize,
        y_0based: usize,
        sequence_length: usize,
    ) -> Option<(usize, usize)> {
        if sequence_length == 0 {
            return None;
        }
        let max_idx = sequence_length.saturating_sub(1);
        let x = x_0based.min(max_idx);
        let y = y_0based.min(max_idx);
        let from = x.min(y);
        let mut to = x.max(y).saturating_add(1).min(sequence_length);
        if to <= from {
            to = (from + 1).min(sequence_length);
        }
        Some((from, to))
    }

    pub(super) fn dotplot_crosshair_selection_bounds_for_mode(
        mode: DotplotMode,
        x_0based: usize,
        y_0based: usize,
        sequence_length: usize,
    ) -> Option<(usize, usize)> {
        if sequence_length == 0 {
            return None;
        }
        match mode {
            DotplotMode::SelfForward | DotplotMode::SelfReverseComplement => {
                Self::dotplot_crosshair_selection_bounds(x_0based, y_0based, sequence_length)
            }
            DotplotMode::PairForward | DotplotMode::PairReverseComplement => {
                let max_idx = sequence_length.saturating_sub(1);
                let from = x_0based.min(max_idx);
                let to = from.saturating_add(1).min(sequence_length);
                Some((from, to))
            }
        }
    }

    pub(super) fn sync_selection_to_dotplot_crosshair(
        &mut self,
        x_0based: usize,
        y_0based: usize,
        request_scroll: bool,
        mode: DotplotMode,
    ) {
        let sequence_length = self.dna.read().ok().map(|dna| dna.len()).unwrap_or(0);
        let Some((from, to)) = Self::dotplot_crosshair_selection_bounds_for_mode(
            mode,
            x_0based,
            y_0based,
            sequence_length,
        ) else {
            return;
        };
        if let Ok(mut display) = self.dna_display.write() {
            display.select(Selection::new(from, to, sequence_length));
        }
        if request_scroll {
            self.map_sequence.request_scroll_to_selection();
        }
    }

    pub(super) fn render_dotplot_status_ui(&self, ui: &mut egui::Ui) {
        ui.separator();
        ui.label(egui::RichText::new("Dotplot status").strong())
            .on_hover_text(
                "Deterministic compute diagnostics for current controls and most recently loaded payload",
            );

        if self.dotplot_ui.overlay_enabled {
            let overlay_choices = self.dotplot_overlay_choices();
            match self.build_dotplot_overlay_estimated_pair_evaluations() {
                Ok(pair_evals) => {
                    let reference_label = self.dotplot_ui.reference_seq_id.trim();
                    ui.label(
                        egui::RichText::new(format!(
                            "request mode=overlay reference={} series={} | word={} step={} mismatches={} tile_bp={} | pair_evals≈{}",
                            if reference_label.is_empty() {
                                "<unset>"
                            } else {
                                reference_label
                            },
                            self.dotplot_ui.overlay_transcript_feature_ids.len(),
                            self.dotplot_ui.word_size.trim(),
                            self.dotplot_ui.step_bp.trim(),
                            self.dotplot_ui.max_mismatches.trim(),
                            if self.dotplot_ui.tile_bp.trim().is_empty() {
                                "-"
                            } else {
                                self.dotplot_ui.tile_bp.trim()
                            },
                            pair_evals
                        ))
                        .monospace()
                        .size(self.feature_details_font_size()),
                    );
                    if !overlay_choices.is_empty() {
                        let labels = overlay_choices
                            .iter()
                            .filter(|choice| {
                                self.dotplot_ui
                                    .overlay_transcript_feature_ids
                                    .contains(&choice.transcript_feature_id)
                            })
                            .map(|choice| {
                                format!("{}({} bp)", choice.label, choice.estimated_length_bp)
                            })
                            .collect::<Vec<_>>();
                        if !labels.is_empty() {
                            ui.label(
                                egui::RichText::new(format!(
                                    "selected isoforms: {}",
                                    labels.join(", ")
                                ))
                                .monospace()
                                .size(self.feature_details_font_size()),
                            );
                        }
                    }
                    if pair_evals > MAX_DOTPLOT_PAIR_EVALUATIONS {
                        ui.colored_label(
                            egui::Color32::from_rgb(180, 83, 9),
                            format!(
                                "Overlay request exceeds the cheap auto-compute budget ({} > {}). The workspace will wait for an explicit compute.",
                                pair_evals, MAX_DOTPLOT_PAIR_EVALUATIONS
                            ),
                        );
                    }
                }
                Err(message) => {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 83, 9),
                        format!("Overlay request diagnostics unavailable: {message}"),
                    );
                }
            }
        } else {
            match self.build_dotplot_compute_diagnostics() {
                Ok(diag) => {
                    let reference_label = diag
                        .reference_seq_id
                        .as_deref()
                        .unwrap_or(diag.seq_id.as_str());
                    ui.label(
                    egui::RichText::new(format!(
                        "request mode={} query={} [{}..{}] reference={} [{}..{}] | word={} step={} mismatches={} tile_bp={}",
                        diag.mode.as_str(),
                        diag.seq_id,
                        diag.query_span_start_0based.saturating_add(1),
                        diag.query_span_end_0based,
                        reference_label,
                        diag.reference_span_start_0based.saturating_add(1),
                        diag.reference_span_end_0based,
                        diag.word_size,
                        diag.step_bp,
                        diag.max_mismatches,
                        diag.tile_bp
                            .map(|value| value.to_string())
                            .unwrap_or_else(|| "-".to_string())
                    ))
                    .monospace()
                    .size(self.feature_details_font_size()),
                );
                    ui.label(
                        egui::RichText::new(format!(
                            "estimated windows query={} reference={} => pair_evals≈{}",
                            diag.query_windows,
                            diag.reference_windows,
                            diag.estimated_pair_evaluations
                        ))
                        .monospace()
                        .size(self.feature_details_font_size()),
                    );
                    if diag.estimated_pair_evaluations == 0 {
                        ui.colored_label(
                        egui::Color32::from_rgb(180, 83, 9),
                        "Estimated pair evaluations are zero; reduce word size or increase span.",
                    );
                    } else if diag.max_mismatches > 0
                        && diag.estimated_pair_evaluations
                            > MAX_DOTPLOT_PAIR_EVALUATIONS.saturating_mul(9) / 10
                    {
                        ui.colored_label(
                        egui::Color32::from_rgb(180, 83, 9),
                        format!(
                            "Request is near compute guardrail ({} of {} pair evaluations). If this fails, increase step_bp or reduce query/reference span.",
                            diag.estimated_pair_evaluations, MAX_DOTPLOT_PAIR_EVALUATIONS
                        ),
                    );
                    } else if diag.max_mismatches == 0
                        && diag.estimated_pair_evaluations > MAX_DOTPLOT_PAIR_EVALUATIONS
                    {
                        ui.small(
                        "Large exact-seed request detected. Engine uses indexed exact matching for mismatches=0, so this can still complete without brute-force pair loops.",
                    );
                    }
                    if Self::dotplot_mode_requires_reference(diag.mode)
                        && diag.reference_span_bp > diag.query_span_bp.saturating_mul(6)
                        && diag.word_size >= 10
                        && diag.step_bp >= 10
                        && diag.max_mismatches == 0
                    {
                        ui.colored_label(
                        egui::Color32::from_rgb(180, 83, 9),
                        "Current settings are very strict for a wide reference span. For cDNA-vs-genomic controls, try smaller word/step or allow mismatches.",
                    );
                    }
                }
                Err(message) => {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 83, 9),
                        format!("Dotplot request diagnostics unavailable: {message}"),
                    );
                }
            }
        }

        if let Some(view) = self.dotplot_cached_view.as_ref() {
            let reference_seq_label = view
                .reference_seq_id
                .as_deref()
                .unwrap_or(view.seq_id.as_str());
            let overlay_mode = Self::dotplot_view_is_overlay(view);
            let reference_span_bp = view
                .reference_span_end_0based
                .saturating_sub(view.reference_span_start_0based);
            let reference_windows = Self::dotplot_window_count(
                reference_span_bp,
                view.word_size.max(1),
                view.step_bp.max(1),
            );
            let query_windows = if overlay_mode {
                view.query_series
                    .iter()
                    .map(|series| {
                        let query_span_bp = series
                            .span_end_0based
                            .saturating_sub(series.span_start_0based);
                        Self::dotplot_window_count(
                            query_span_bp,
                            view.word_size.max(1),
                            view.step_bp.max(1),
                        )
                    })
                    .sum()
            } else {
                let query_span_bp = view.span_end_0based.saturating_sub(view.span_start_0based);
                Self::dotplot_window_count(
                    query_span_bp,
                    view.word_size.max(1),
                    view.step_bp.max(1),
                )
            };
            let estimated_pair_evaluations = query_windows.saturating_mul(reference_windows);
            let point_count = Self::dotplot_view_total_point_count(view);
            let estimated_hit_fraction = if estimated_pair_evaluations == 0 {
                0.0
            } else {
                (point_count as f64 / estimated_pair_evaluations as f64) * 100.0
            };
            if overlay_mode {
                let annotation_count = view
                    .reference_annotation
                    .as_ref()
                    .map(|track| track.interval_count)
                    .unwrap_or(0);
                let labels = view
                    .query_series
                    .iter()
                    .map(|series| series.label.clone())
                    .collect::<Vec<_>>();
                ui.label(
                    egui::RichText::new(format!(
                        "loaded payload '{}' generated_at={} overlay series={} total_points={} hit_fraction≈{:.6}% owner={} reference={} [{}..{}] merged_exons={} | {}",
                        view.dotplot_id,
                        view.generated_at_unix_ms,
                        view.series_count.max(view.query_series.len()),
                        point_count,
                        estimated_hit_fraction,
                        view.owner_seq_id,
                        reference_seq_label,
                        view.reference_span_start_0based.saturating_add(1),
                        view.reference_span_end_0based,
                        annotation_count,
                        labels.join(", ")
                    ))
                    .monospace()
                    .size(self.feature_details_font_size()),
                );
            } else {
                let non_empty_box_bins = view
                    .boxplot_bins
                    .iter()
                    .filter(|bin| bin.hit_count > 0)
                    .count();
                ui.label(
                    egui::RichText::new(format!(
                        "loaded payload '{}' generated_at={} mode={} points={} hit_fraction≈{:.6}% box_bins={}/{} | query={} [{}..{}] reference={} [{}..{}]",
                        view.dotplot_id,
                        view.generated_at_unix_ms,
                        view.mode.as_str(),
                        point_count,
                        estimated_hit_fraction,
                        non_empty_box_bins,
                        view.boxplot_bin_count,
                        view.seq_id,
                        view.span_start_0based.saturating_add(1),
                        view.span_end_0based,
                        reference_seq_label,
                        view.reference_span_start_0based.saturating_add(1),
                        view.reference_span_end_0based
                    ))
                    .monospace()
                    .size(self.feature_details_font_size()),
                );
            }
            let same_reference = reference_seq_label == view.seq_id.as_str();
            if matches!(view.mode, DotplotMode::SelfForward) {
                ui.small(
                    "Self-forward compares a sequence to itself. A main diagonal is expected identity; off-diagonal signal indicates repeated motifs.",
                );
            }
            if point_count == 0 {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 24, 93),
                    "Loaded payload has zero seed hits. Try smaller word size, higher mismatches, or a narrower reference span around the expected locus.",
                );
                if same_reference
                    && matches!(
                        view.mode,
                        DotplotMode::PairReverseComplement | DotplotMode::SelfReverseComplement
                    )
                {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 83, 9),
                        "Reverse-complement self-pair mode is an inverted-repeat scan. Zero hits can be valid with sparse sampling/strict seeds; try smaller step and word for sensitivity.",
                    );
                }
            } else if !overlay_mode
                && Self::dotplot_mode_requires_reference(view.mode)
                && point_count <= DOTPLOT_SPARSE_POINT_HINT_THRESHOLD
            {
                if matches!(view.mode, DotplotMode::PairForward) {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 83, 9),
                        "Sparse pair_forward hits can indicate strand/orientation mismatch (common for cDNA vs genomic). Try pair_reverse_complement.",
                    );
                } else {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 83, 9),
                        "Pair reverse-complement is still sparse. Increase sensitivity (smaller word/step, allow mismatches) or restrict reference span.",
                    );
                    if let Some(forward_hits) = self.comparable_pair_forward_point_count(view)
                        && forward_hits > 0
                    {
                        ui.colored_label(
                            egui::Color32::from_rgb(180, 83, 9),
                            format!(
                                "Pair-forward has {forward_hits} hits for identical spans/seed settings. This typically means query and reference are already oriented in the same biological direction for this extraction."
                            ),
                        );
                    }
                }
                if let Some((hit_start, hit_end)) = Self::dotplot_reference_hit_envelope(view, 0) {
                    let left_gap = hit_start.saturating_sub(view.reference_span_start_0based);
                    let right_gap = view.reference_span_end_0based.saturating_sub(hit_end);
                    if left_gap <= DOTPLOT_HIT_ENVELOPE_PADDING_BP / 2
                        || right_gap <= DOTPLOT_HIT_ENVELOPE_PADDING_BP / 2
                    {
                        ui.colored_label(
                            egui::Color32::from_rgb(180, 83, 9),
                            "Detected hits are near reference-span edge. If you are using explicit ref_start/ref_end, click 'Fit ref span to hits' and recompute for a less edge-compressed map.",
                        );
                    }
                }
            }
        } else {
            ui.small("No dotplot payload is currently loaded.");
        }

        if !self.dotplot_last_compute_status.trim().is_empty() {
            ui.label(
                egui::RichText::new(format!("compute log: {}", self.dotplot_last_compute_status))
                    .monospace()
                    .size(self.feature_details_font_size()),
            );
        }

        if !self.op_status.trim().is_empty() {
            if let Some(first_line) = self.op_status.lines().next() {
                ui.label(
                    egui::RichText::new(format!("last operation: {first_line}"))
                        .monospace()
                        .size(self.feature_details_font_size()),
                );
            }
            if let Some(warnings_line) = self
                .op_status
                .lines()
                .find(|line| line.starts_with("warnings:") && !line.ends_with(" -"))
            {
                ui.label(
                    egui::RichText::new(warnings_line)
                        .monospace()
                        .size(self.feature_details_font_size()),
                );
            }
            if let Some(messages_line) = self
                .op_status
                .lines()
                .find(|line| line.starts_with("messages:") && !line.ends_with(" -"))
            {
                ui.label(
                    egui::RichText::new(messages_line)
                        .monospace()
                        .size(self.feature_details_font_size()),
                );
            }
        }
    }

    pub(super) fn render_dotplot_density_ui(
        &mut self,
        ui: &mut egui::Ui,
        view: &DotplotView,
        flex_track: Option<&FlexibilityTrack>,
        density_threshold: f32,
        intensity_gain: f32,
    ) {
        #[derive(Clone)]
        struct OverlaySeriesRenderData {
            label: String,
            seq_id: String,
            color: egui::Color32,
            span_start_0based: usize,
            span_end_0based: usize,
            shift_bp: usize,
            cells: HashMap<(i32, i32), usize>,
            visible_cells: HashMap<(i32, i32), egui::Color32>,
        }

        let overlay_mode = Self::dotplot_view_is_overlay(view);
        let has_reference_annotation = view
            .reference_annotation
            .as_ref()
            .is_some_and(|track| !track.intervals.is_empty());
        let desired_w = ui.available_width().max(740.0);
        let desired_h = if flex_track.is_some() { 470.0 } else { 390.0 };
        let (canvas_rect, response) =
            ui.allocate_exact_size(Vec2::new(desired_w, desired_h), egui::Sense::hover());
        let painter = ui.painter_at(canvas_rect);
        painter.rect_filled(canvas_rect, 6.0, egui::Color32::from_rgb(248, 250, 252));
        painter.rect_stroke(
            canvas_rect,
            6.0,
            egui::Stroke::new(1.0, egui::Color32::from_rgb(203, 213, 225)),
            egui::StrokeKind::Inside,
        );

        let top_margin = if overlay_mode { 48.0 } else { 26.0 };
        let left_margin = if has_reference_annotation { 78.0 } else { 56.0 };
        let right_margin = 16.0;
        let bottom_margin = 20.0;
        let flex_height = if flex_track.is_some() { 108.0 } else { 0.0 };
        let gap = if flex_track.is_some() { 10.0 } else { 0.0 };
        let dotplot_rect = egui::Rect::from_min_max(
            egui::pos2(
                canvas_rect.left() + left_margin,
                canvas_rect.top() + top_margin,
            ),
            egui::pos2(
                canvas_rect.right() - right_margin,
                canvas_rect.bottom() - bottom_margin - flex_height - gap,
            ),
        );
        painter.rect_filled(dotplot_rect, 4.0, egui::Color32::from_rgb(255, 255, 255));
        painter.rect_stroke(
            dotplot_rect,
            4.0,
            egui::Stroke::new(1.0, egui::Color32::from_rgb(203, 213, 225)),
            egui::StrokeKind::Inside,
        );

        for tick_idx in 1..10 {
            let frac = tick_idx as f32 / 10.0;
            let x = dotplot_rect.left() + frac * dotplot_rect.width();
            let y = dotplot_rect.top() + frac * dotplot_rect.height();
            painter.line_segment(
                [
                    egui::pos2(x, dotplot_rect.top()),
                    egui::pos2(x, dotplot_rect.bottom()),
                ],
                egui::Stroke::new(0.5, egui::Color32::from_rgb(241, 245, 249)),
            );
            painter.line_segment(
                [
                    egui::pos2(dotplot_rect.left(), y),
                    egui::pos2(dotplot_rect.right(), y),
                ],
                egui::Stroke::new(0.5, egui::Color32::from_rgb(241, 245, 249)),
            );
        }
        if overlay_mode {
            let overlay_x_axis_mode = self.dotplot_ui.overlay_x_axis_mode;
            let reference_span = view
                .reference_span_end_0based
                .saturating_sub(view.reference_span_start_0based)
                .max(1);
            let reference_span_max = reference_span.saturating_sub(1).max(1);
            let uses_shifted_anchor_axis = matches!(
                overlay_x_axis_mode,
                DotplotOverlayXAxisMode::SharedExonAnchor | DotplotOverlayXAxisMode::QueryAnchorBp
            );
            let resolved_anchor_series = match overlay_x_axis_mode {
                DotplotOverlayXAxisMode::SharedExonAnchor => self
                    .dotplot_ui
                    .overlay_anchor_exon
                    .as_ref()
                    .map(|exon| view.resolve_overlay_anchor_series(exon))
                    .unwrap_or_default(),
                DotplotOverlayXAxisMode::QueryAnchorBp => view.resolve_query_anchor_series(),
                _ => vec![],
            };
            let rendered_series_source = if uses_shifted_anchor_axis {
                resolved_anchor_series
                    .iter()
                    .filter_map(|resolved| {
                        view.query_series
                            .get(resolved.series_index)
                            .map(|series| (series, resolved.shift_bp))
                    })
                    .collect::<Vec<_>>()
            } else {
                view.query_series
                    .iter()
                    .map(|series| (series, 0usize))
                    .collect()
            };
            let max_query_span = rendered_series_source
                .iter()
                .map(|(series, _)| {
                    series
                        .span_end_0based
                        .saturating_sub(series.span_start_0based)
                        .max(1)
                })
                .max()
                .unwrap_or(1);
            let average_query_span = rendered_series_source
                .iter()
                .map(|(series, _)| {
                    series
                        .span_end_0based
                        .saturating_sub(series.span_start_0based)
                        .max(1)
                })
                .sum::<usize>()
                .checked_div(rendered_series_source.len().max(1))
                .unwrap_or(1)
                .max(1);
            let plotted_query_span = if uses_shifted_anchor_axis {
                resolved_anchor_series
                    .iter()
                    .map(|resolved| resolved.plotted_span_end_0based)
                    .max()
                    .unwrap_or(1)
                    .max(1)
            } else {
                overlay_x_axis_mode.plot_query_span_bp(max_query_span, average_query_span)
            };
            let plotted_query_span_max = plotted_query_span.saturating_sub(1).max(1);
            let (x_axis_start_label, x_axis_end_label) =
                overlay_x_axis_mode.axis_edge_labels(plotted_query_span);
            let anchor_status_message = match overlay_x_axis_mode {
                DotplotOverlayXAxisMode::SharedExonAnchor => {
                    if self.dotplot_ui.overlay_anchor_exon.is_none() {
                        Some("No shared exon anchor selected.".to_string())
                    } else if rendered_series_source.len() < 2 {
                        Some(
                            "Selected shared exon is not present in at least two plotted transcripts."
                                .to_string(),
                        )
                    } else {
                        None
                    }
                }
                DotplotOverlayXAxisMode::QueryAnchorBp => {
                    let anchored_series_count = view
                        .query_series
                        .iter()
                        .filter(|series| series.query_anchor_0based.is_some())
                        .count();
                    if anchored_series_count == 0 {
                        Some(
                            "This overlay payload has no explicit manual/domain query anchors."
                                .to_string(),
                        )
                    } else if rendered_series_source.len() < 2 {
                        Some(
                            "At least two plotted query series need manual/domain anchors."
                                .to_string(),
                        )
                    } else {
                        None
                    }
                }
                _ => None,
            };
            let reference_seq_label = view
                .reference_seq_id
                .as_deref()
                .unwrap_or(view.seq_id.as_str());
            let cols = dotplot_rect.width().max(2.0).round() as i32;
            let rows = dotplot_rect.height().max(2.0).round() as i32;
            let series_sample_cap =
                (DOTPLOT_RENDER_MAX_POINTS / rendered_series_source.len().max(1)).max(1);
            let density_threshold = density_threshold.clamp(0.0, 0.99);
            let intensity_gain = intensity_gain.clamp(0.1, 16.0);
            let mut rendered_series: Vec<OverlaySeriesRenderData> = vec![];
            let mut total_rendered_cells = 0usize;

            for (series, shift_bp) in &rendered_series_source {
                let sample_stride = (series.points.len() / series_sample_cap).max(1);
                let mut cells: HashMap<(i32, i32), usize> = HashMap::new();
                for point in series.points.iter().step_by(sample_stride) {
                    let y_local = point
                        .y_0based
                        .saturating_sub(view.reference_span_start_0based)
                        .min(reference_span_max);
                    let x_frac = if uses_shifted_anchor_axis {
                        let x_local = point.x_0based.saturating_sub(series.span_start_0based).min(
                            series
                                .span_end_0based
                                .saturating_sub(series.span_start_0based)
                                .saturating_sub(1)
                                .max(1),
                        );
                        (shift_bp.saturating_add(x_local) as f32 / plotted_query_span_max as f32)
                            .clamp(0.0, 1.0)
                    } else {
                        overlay_x_axis_mode.point_fraction(
                            point.x_0based,
                            series.span_start_0based,
                            series.span_end_0based,
                            max_query_span,
                        )
                    };
                    let y_frac = (y_local as f32 / reference_span_max as f32).clamp(0.0, 1.0);
                    let x_cell = ((x_frac * (cols - 1) as f32).round() as i32).clamp(0, cols - 1);
                    let y_cell = ((y_frac * (rows - 1) as f32).round() as i32).clamp(0, rows - 1);
                    let entry = cells.entry((x_cell, y_cell)).or_insert(0);
                    *entry = entry.saturating_add(1);
                }
                let max_cell_count = cells.values().copied().max().unwrap_or(1).max(1) as f32;
                let mut visible_cells: HashMap<(i32, i32), egui::Color32> = HashMap::new();
                for ((x_cell, y_cell), count) in &cells {
                    let density_raw = (*count as f32 / max_cell_count).clamp(0.0, 1.0);
                    if density_raw < density_threshold {
                        continue;
                    }
                    let normalized = if density_threshold <= 0.0 {
                        density_raw
                    } else {
                        ((density_raw - density_threshold) / (1.0 - density_threshold))
                            .clamp(0.0, 1.0)
                    };
                    let density = (normalized * intensity_gain).clamp(0.0, 1.0).sqrt();
                    let alpha = (72.0 + 178.0 * density).round() as u8;
                    let color = egui::Color32::from_rgba_unmultiplied(
                        series.color_rgb[0],
                        series.color_rgb[1],
                        series.color_rgb[2],
                        alpha,
                    );
                    let x0 =
                        dotplot_rect.left() + (*x_cell as f32 / cols as f32) * dotplot_rect.width();
                    let x1 = dotplot_rect.left()
                        + ((*x_cell + 1) as f32 / cols as f32) * dotplot_rect.width();
                    let y0 =
                        dotplot_rect.top() + (*y_cell as f32 / rows as f32) * dotplot_rect.height();
                    let y1 = dotplot_rect.top()
                        + ((*y_cell + 1) as f32 / rows as f32) * dotplot_rect.height();
                    painter.rect_filled(
                        egui::Rect::from_min_max(egui::pos2(x0, y0), egui::pos2(x1, y1)),
                        0.0,
                        color,
                    );
                    visible_cells.insert((*x_cell, *y_cell), color);
                    total_rendered_cells = total_rendered_cells.saturating_add(1);
                }
                if !visible_cells.is_empty()
                    && visible_cells.len()
                        <= DOTPLOT_CONNECT_DIAGONALS_MAX_CELLS
                            .checked_div(view.query_series.len().max(1))
                            .unwrap_or(DOTPLOT_CONNECT_DIAGONALS_MAX_CELLS)
                            .max(1)
                {
                    let cell_w = dotplot_rect.width() / cols as f32;
                    let cell_h = dotplot_rect.height() / rows as f32;
                    let stroke_width = (cell_w.min(cell_h) * 0.45).clamp(0.8, 2.0);
                    for ((x_cell, y_cell), color) in &visible_cells {
                        for dy in [-1, 0, 1] {
                            let neighbor = (*x_cell + 1, *y_cell + dy);
                            if visible_cells.contains_key(&neighbor) {
                                let x0 = dotplot_rect.left() + (*x_cell as f32 + 0.5) * cell_w;
                                let y0 = dotplot_rect.top() + (*y_cell as f32 + 0.5) * cell_h;
                                let x1 = dotplot_rect.left() + (neighbor.0 as f32 + 0.5) * cell_w;
                                let y1 = dotplot_rect.top() + (neighbor.1 as f32 + 0.5) * cell_h;
                                painter.line_segment(
                                    [egui::pos2(x0, y0), egui::pos2(x1, y1)],
                                    egui::Stroke::new(stroke_width, color.gamma_multiply(0.5)),
                                );
                            }
                        }
                    }
                }
                rendered_series.push(OverlaySeriesRenderData {
                    label: series.label.clone(),
                    seq_id: series.seq_id.clone(),
                    color: egui::Color32::from_rgb(
                        series.color_rgb[0],
                        series.color_rgb[1],
                        series.color_rgb[2],
                    ),
                    span_start_0based: series.span_start_0based,
                    span_end_0based: series.span_end_0based,
                    shift_bp: *shift_bp,
                    cells,
                    visible_cells,
                });
            }

            if total_rendered_cells == 0 {
                painter.text(
                    dotplot_rect.center(),
                    egui::Align2::CENTER_CENTER,
                    anchor_status_message.unwrap_or_else(|| {
                        format!(
                            "No visible cells (threshold={:.2}, points={}, series={})",
                            density_threshold,
                            rendered_series_source
                                .iter()
                                .map(|(series, _)| series.point_count)
                                .sum::<usize>(),
                            rendered_series_source.len()
                        )
                    }),
                    egui::FontId::monospace(11.0),
                    egui::Color32::from_rgb(100, 116, 139),
                );
            }

            if let Some(track) = view
                .reference_annotation
                .as_ref()
                .filter(|track| !track.intervals.is_empty())
            {
                let annotation_rect = egui::Rect::from_min_max(
                    egui::pos2(dotplot_rect.left() - 16.0, dotplot_rect.top()),
                    egui::pos2(dotplot_rect.left() - 6.0, dotplot_rect.bottom()),
                );
                painter.rect_filled(annotation_rect, 2.0, egui::Color32::from_rgb(248, 250, 252));
                painter.rect_stroke(
                    annotation_rect,
                    2.0,
                    egui::Stroke::new(1.0, egui::Color32::from_rgb(203, 213, 225)),
                    egui::StrokeKind::Inside,
                );
                let reference_span_f32 = reference_span.max(1) as f32;
                for interval in &track.intervals {
                    let local_start = interval
                        .start_0based
                        .saturating_sub(view.reference_span_start_0based)
                        .min(reference_span);
                    let local_end = interval
                        .end_0based_exclusive
                        .saturating_sub(view.reference_span_start_0based)
                        .min(reference_span);
                    let y0 = annotation_rect.top()
                        + (local_start as f32 / reference_span_f32) * annotation_rect.height();
                    let y1 = annotation_rect.top()
                        + (local_end as f32 / reference_span_f32) * annotation_rect.height();
                    painter.rect_filled(
                        egui::Rect::from_min_max(
                            egui::pos2(annotation_rect.left() + 1.0, y0),
                            egui::pos2(annotation_rect.right() - 1.0, y1.max(y0 + 1.0)),
                        ),
                        1.0,
                        egui::Color32::from_rgb(34, 197, 94),
                    );
                }
                painter.text(
                    egui::pos2(annotation_rect.center().x, canvas_rect.top() + 4.0),
                    egui::Align2::CENTER_TOP,
                    track.label.as_str(),
                    egui::FontId::monospace(9.0),
                    egui::Color32::from_rgb(51, 65, 85),
                );
            }

            let mut legend_cursor = egui::pos2(dotplot_rect.left(), canvas_rect.top() + 4.0);
            for series in &rendered_series {
                let swatch_rect = egui::Rect::from_min_max(
                    legend_cursor,
                    egui::pos2(legend_cursor.x + 10.0, legend_cursor.y + 10.0),
                );
                painter.rect_filled(swatch_rect, 1.0, series.color);
                let text_pos = egui::pos2(swatch_rect.right() + 4.0, legend_cursor.y - 1.0);
                painter.text(
                    text_pos,
                    egui::Align2::LEFT_TOP,
                    series.label.as_str(),
                    egui::FontId::monospace(9.5),
                    egui::Color32::from_rgb(51, 65, 85),
                );
                let estimated_width = 22.0 + series.label.len() as f32 * 6.1;
                if legend_cursor.x + estimated_width > dotplot_rect.right() - 120.0 {
                    legend_cursor.x = dotplot_rect.left();
                    legend_cursor.y += 12.0;
                } else {
                    legend_cursor.x += estimated_width;
                }
            }

            painter.text(
                egui::pos2(dotplot_rect.left(), canvas_rect.top() + 20.0),
                egui::Align2::LEFT_TOP,
                overlay_x_axis_mode.axis_label(),
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(51, 65, 85),
            );
            if overlay_x_axis_mode == DotplotOverlayXAxisMode::SharedExonAnchor
                && let Some(anchor) = self.dotplot_ui.overlay_anchor_exon.as_ref()
            {
                painter.text(
                    egui::pos2(dotplot_rect.left() + 220.0, canvas_rect.top() + 20.0),
                    egui::Align2::LEFT_TOP,
                    format!("anchor exon {}", anchor.token()),
                    egui::FontId::monospace(10.0),
                    egui::Color32::from_rgb(71, 85, 105),
                );
            } else if overlay_x_axis_mode == DotplotOverlayXAxisMode::QueryAnchorBp
                && let Some(anchor_label) = view.query_anchor_label()
            {
                painter.text(
                    egui::pos2(dotplot_rect.left() + 220.0, canvas_rect.top() + 20.0),
                    egui::Align2::LEFT_TOP,
                    format!("anchor {anchor_label}"),
                    egui::FontId::monospace(10.0),
                    egui::Color32::from_rgb(71, 85, 105),
                );
            }
            painter.text(
                egui::pos2(dotplot_rect.right(), canvas_rect.top() + 20.0),
                egui::Align2::RIGHT_TOP,
                format!("y: {reference_seq_label}"),
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(51, 65, 85),
            );
            painter.text(
                egui::pos2(dotplot_rect.left(), dotplot_rect.bottom() + 2.0),
                egui::Align2::LEFT_TOP,
                x_axis_start_label.as_str(),
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(71, 85, 105),
            );
            painter.text(
                egui::pos2(dotplot_rect.right(), dotplot_rect.bottom() + 2.0),
                egui::Align2::RIGHT_TOP,
                x_axis_end_label.as_str(),
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(71, 85, 105),
            );
            painter.text(
                egui::pos2(dotplot_rect.left() - 4.0, dotplot_rect.top()),
                egui::Align2::RIGHT_TOP,
                format!("{}", view.reference_span_start_0based.saturating_add(1)),
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(71, 85, 105),
            );
            painter.text(
                egui::pos2(dotplot_rect.left() - 4.0, dotplot_rect.bottom()),
                egui::Align2::RIGHT_BOTTOM,
                format!("{}", view.reference_span_end_0based),
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(71, 85, 105),
            );
            painter.text(
                egui::pos2(dotplot_rect.left() + 2.0, dotplot_rect.top() + 2.0),
                egui::Align2::LEFT_TOP,
                "y",
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(100, 116, 139),
            );
            painter.text(
                egui::pos2(dotplot_rect.right() - 2.0, dotplot_rect.bottom() - 2.0),
                egui::Align2::RIGHT_BOTTOM,
                "x",
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(100, 116, 139),
            );

            self.dotplot_hover_crosshair_bp = None;
            if response.hovered()
                && let Some(pointer) = response.hover_pos()
                && dotplot_rect.contains(pointer)
            {
                let fx = ((pointer.x - dotplot_rect.left()) / dotplot_rect.width()).clamp(0.0, 1.0);
                let fy = ((pointer.y - dotplot_rect.top()) / dotplot_rect.height()).clamp(0.0, 1.0);
                let y_bp = view.reference_span_start_0based
                    + (fy * reference_span_max as f32).round() as usize;
                let x_cell = ((fx * (cols - 1) as f32).round() as i32).clamp(0, cols - 1);
                let y_cell = ((fy * (rows - 1) as f32).round() as i32).clamp(0, rows - 1);
                let pointer_x = dotplot_rect.left() + fx * dotplot_rect.width();
                let pointer_y = dotplot_rect.top() + fy * dotplot_rect.height();
                painter.line_segment(
                    [
                        egui::pos2(pointer_x, dotplot_rect.top()),
                        egui::pos2(pointer_x, dotplot_rect.bottom()),
                    ],
                    egui::Stroke::new(1.0, egui::Color32::from_rgb(15, 23, 42)),
                );
                painter.line_segment(
                    [
                        egui::pos2(dotplot_rect.left(), pointer_y),
                        egui::pos2(dotplot_rect.right(), pointer_y),
                    ],
                    egui::Stroke::new(1.0, egui::Color32::from_rgb(15, 23, 42)),
                );
                painter.circle_filled(
                    egui::pos2(pointer_x, pointer_y),
                    2.0,
                    egui::Color32::from_rgb(15, 23, 42),
                );
                painter.text(
                    egui::pos2(dotplot_rect.left() + 4.0, dotplot_rect.top() + 14.0),
                    egui::Align2::LEFT_TOP,
                    format!(
                        "hover reference({})={} | {} isoforms",
                        reference_seq_label,
                        y_bp.saturating_add(1),
                        rendered_series.len()
                    ),
                    egui::FontId::monospace(10.0),
                    egui::Color32::from_rgb(15, 23, 42),
                );
                response.clone().on_hover_ui_at_pointer(|ui| {
                    ui.label(
                        egui::RichText::new(format!(
                            "reference({})={}",
                            reference_seq_label,
                            y_bp.saturating_add(1)
                        ))
                        .monospace()
                        .size(self.feature_details_font_size()),
                    );
                    for series in &rendered_series {
                        let cell_density = series.cells.get(&(x_cell, y_cell)).copied().unwrap_or(0);
                        let rendered = series.visible_cells.contains_key(&(x_cell, y_cell));
                        let query_label = if uses_shifted_anchor_axis {
                            let global_bp =
                                (fx * plotted_query_span_max as f32).round() as usize;
                            let series_span_max = series
                                .span_end_0based
                                .saturating_sub(series.span_start_0based)
                                .saturating_sub(1)
                                .max(1);
                            if global_bp < series.shift_bp
                                || global_bp
                                    > series.shift_bp.saturating_add(series_span_max)
                            {
                                "x=outside".to_string()
                            } else {
                                format!(
                                    "x={}",
                                    series
                                        .span_start_0based
                                        .saturating_add(global_bp - series.shift_bp)
                                        .saturating_add(1)
                                )
                            }
                        } else if let Some(query_bp) = overlay_x_axis_mode.query_coordinate_at_fraction(
                            fx,
                            series.span_start_0based,
                            series.span_end_0based,
                            max_query_span,
                        ) {
                            format!("x={}", query_bp.saturating_add(1))
                        } else if matches!(overlay_x_axis_mode, DotplotOverlayXAxisMode::PercentLength)
                        {
                            "x=n/a".to_string()
                        } else {
                            "x=outside".to_string()
                        };
                        ui.label(
                            egui::RichText::new(format!(
                                "{} ({}) {} cell={}{}",
                                series.label,
                                series.seq_id,
                                query_label,
                                cell_density,
                                if rendered { " visible" } else { "" }
                            ))
                            .monospace()
                            .color(series.color)
                            .size(self.feature_details_font_size()),
                        );
                    }
                    ui.small(
                        "Overlay mode uses one shared reference axis plus per-isoform x readouts according to the selected overlay x-axis layout. Selection sync and locked crosshair stay disabled.",
                    );
                });
            }
            if response.clicked_by(egui::PointerButton::Secondary) {
                self.dotplot_locked_crosshair_bp = None;
                self.dotplot_hover_crosshair_bp = None;
            }

            if let Some(track) = flex_track {
                let flex_rect = egui::Rect::from_min_max(
                    egui::pos2(dotplot_rect.left(), dotplot_rect.bottom() + gap),
                    egui::pos2(dotplot_rect.right(), canvas_rect.bottom() - bottom_margin),
                );
                painter.rect_filled(flex_rect, 4.0, egui::Color32::from_rgb(255, 255, 255));
                painter.rect_stroke(
                    flex_rect,
                    4.0,
                    egui::Stroke::new(1.0, egui::Color32::from_rgb(203, 213, 225)),
                    egui::StrokeKind::Inside,
                );
                let score_span = (track.max_score - track.min_score).abs().max(1e-12);
                let span_bp = track
                    .span_end_0based
                    .saturating_sub(track.span_start_0based)
                    .max(1) as f32;
                let mut points: Vec<egui::Pos2> = vec![];
                for bin in &track.bins {
                    let x_frac = (bin.start_0based.saturating_sub(track.span_start_0based) as f32
                        / span_bp)
                        .clamp(0.0, 1.0);
                    let y_frac =
                        ((bin.score - track.min_score) / score_span).clamp(0.0, 1.0) as f32;
                    points.push(egui::pos2(
                        flex_rect.left() + x_frac * flex_rect.width(),
                        flex_rect.bottom() - y_frac * flex_rect.height(),
                    ));
                }
                for segment in points.windows(2) {
                    let [left, right] = [segment[0], segment[1]];
                    painter.line_segment(
                        [left, right],
                        egui::Stroke::new(1.5, egui::Color32::from_rgb(22, 101, 52)),
                    );
                }
            }
            return;
        }
        let is_self_mode = matches!(
            view.mode,
            DotplotMode::SelfForward | DotplotMode::SelfReverseComplement
        );
        let selection_sync_enabled = self.dotplot_selection_sync_enabled_for_view(view);

        let query_span = view
            .span_end_0based
            .saturating_sub(view.span_start_0based)
            .max(1);
        let query_span_max = query_span.saturating_sub(1).max(1);
        let reference_span = view
            .reference_span_end_0based
            .saturating_sub(view.reference_span_start_0based)
            .max(1);
        let reference_span_max = reference_span.saturating_sub(1).max(1);
        let reference_seq_label = view
            .reference_seq_id
            .as_deref()
            .unwrap_or(view.seq_id.as_str());
        let cols = dotplot_rect.width().max(2.0).round() as i32;
        let rows = dotplot_rect.height().max(2.0).round() as i32;
        let sample_stride = (view.points.len() / DOTPLOT_RENDER_MAX_POINTS).max(1);
        let mut cells: HashMap<(i32, i32), (usize, usize)> = HashMap::new();
        for point in view.points.iter().step_by(sample_stride) {
            let x_local = point
                .x_0based
                .saturating_sub(view.span_start_0based)
                .min(query_span_max);
            let y_local = point
                .y_0based
                .saturating_sub(view.reference_span_start_0based)
                .min(reference_span_max);
            let x_frac = (x_local as f32 / query_span_max as f32).clamp(0.0, 1.0);
            let y_frac = (y_local as f32 / reference_span_max as f32).clamp(0.0, 1.0);
            let x_cell = ((x_frac * (cols - 1) as f32).round() as i32).clamp(0, cols - 1);
            let y_cell = ((y_frac * (rows - 1) as f32).round() as i32).clamp(0, rows - 1);
            let entry = cells
                .entry((x_cell, y_cell))
                .or_insert((0, point.mismatches));
            entry.0 = entry.0.saturating_add(1);
            entry.1 = entry.1.min(point.mismatches);
        }
        let max_cell_count = cells.values().map(|(count, _)| *count).max().unwrap_or(1) as f32;
        let density_threshold = density_threshold.clamp(0.0, 0.99);
        let intensity_gain = intensity_gain.clamp(0.1, 16.0);
        let mut rendered_cells = 0usize;
        let mut visible_cells: HashMap<(i32, i32), egui::Color32> = HashMap::new();
        for ((x_cell, y_cell), (count, min_mismatch)) in &cells {
            let density_raw = (*count as f32 / max_cell_count).clamp(0.0, 1.0);
            if density_raw < density_threshold {
                continue;
            }
            let normalized = if density_threshold <= 0.0 {
                density_raw
            } else {
                ((density_raw - density_threshold) / (1.0 - density_threshold)).clamp(0.0, 1.0)
            };
            let density = (normalized * intensity_gain).clamp(0.0, 1.0).sqrt();
            let mismatch_frac = if view.max_mismatches == 0 {
                0.0
            } else {
                (*min_mismatch as f32 / view.max_mismatches as f32).clamp(0.0, 1.0)
            };
            let rgb = Self::mix_rgb([29, 78, 216], [180, 83, 9], mismatch_frac);
            let alpha = (90.0 + 165.0 * density).round() as u8;
            let color = egui::Color32::from_rgba_unmultiplied(rgb[0], rgb[1], rgb[2], alpha);
            let x0 = dotplot_rect.left() + (*x_cell as f32 / cols as f32) * dotplot_rect.width();
            let x1 =
                dotplot_rect.left() + ((*x_cell + 1) as f32 / cols as f32) * dotplot_rect.width();
            let y0 = dotplot_rect.top() + (*y_cell as f32 / rows as f32) * dotplot_rect.height();
            let y1 =
                dotplot_rect.top() + ((*y_cell + 1) as f32 / rows as f32) * dotplot_rect.height();
            painter.rect_filled(
                egui::Rect::from_min_max(egui::pos2(x0, y0), egui::pos2(x1, y1)),
                0.0,
                color,
            );
            visible_cells.insert((*x_cell, *y_cell), color);
            rendered_cells = rendered_cells.saturating_add(1);
        }
        if !visible_cells.is_empty() && visible_cells.len() <= DOTPLOT_CONNECT_DIAGONALS_MAX_CELLS {
            let cell_w = dotplot_rect.width() / cols as f32;
            let cell_h = dotplot_rect.height() / rows as f32;
            let stroke_width = (cell_w.min(cell_h) * 0.45).clamp(0.8, 2.0);
            for ((x_cell, y_cell), color) in &visible_cells {
                for dy in [-1, 0, 1] {
                    let neighbor = (*x_cell + 1, *y_cell + dy);
                    if visible_cells.contains_key(&neighbor) {
                        let x0 = dotplot_rect.left() + (*x_cell as f32 + 0.5) * cell_w;
                        let y0 = dotplot_rect.top() + (*y_cell as f32 + 0.5) * cell_h;
                        let x1 = dotplot_rect.left() + (neighbor.0 as f32 + 0.5) * cell_w;
                        let y1 = dotplot_rect.top() + (neighbor.1 as f32 + 0.5) * cell_h;
                        painter.line_segment(
                            [egui::pos2(x0, y0), egui::pos2(x1, y1)],
                            egui::Stroke::new(stroke_width, color.gamma_multiply(0.55)),
                        );
                    }
                }
            }
        }
        if rendered_cells == 0 {
            painter.text(
                dotplot_rect.center(),
                egui::Align2::CENTER_CENTER,
                format!(
                    "No visible cells (threshold={:.2}, points={})",
                    density_threshold, view.point_count
                ),
                egui::FontId::monospace(11.0),
                egui::Color32::from_rgb(100, 116, 139),
            );
        }

        painter.text(
            egui::pos2(dotplot_rect.left(), canvas_rect.top() + 4.0),
            egui::Align2::LEFT_TOP,
            format!("x: {}", view.seq_id),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(51, 65, 85),
        );
        painter.text(
            egui::pos2(dotplot_rect.right(), canvas_rect.top() + 4.0),
            egui::Align2::RIGHT_TOP,
            format!("y: {reference_seq_label}"),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(51, 65, 85),
        );
        painter.text(
            egui::pos2(dotplot_rect.left(), dotplot_rect.bottom() + 2.0),
            egui::Align2::LEFT_TOP,
            format!("{}", view.span_start_0based.saturating_add(1)),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(71, 85, 105),
        );
        painter.text(
            egui::pos2(dotplot_rect.right(), dotplot_rect.bottom() + 2.0),
            egui::Align2::RIGHT_TOP,
            format!("{}", view.span_end_0based),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(71, 85, 105),
        );
        painter.text(
            egui::pos2(dotplot_rect.left() - 4.0, dotplot_rect.top()),
            egui::Align2::RIGHT_TOP,
            format!("{}", view.reference_span_start_0based.saturating_add(1)),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(71, 85, 105),
        );
        painter.text(
            egui::pos2(dotplot_rect.left() - 4.0, dotplot_rect.bottom()),
            egui::Align2::RIGHT_BOTTOM,
            format!("{}", view.reference_span_end_0based),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(71, 85, 105),
        );
        painter.text(
            egui::pos2(dotplot_rect.left() + 2.0, dotplot_rect.top() + 2.0),
            egui::Align2::LEFT_TOP,
            "y",
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(100, 116, 139),
        );
        painter.text(
            egui::pos2(dotplot_rect.right() - 2.0, dotplot_rect.bottom() - 2.0),
            egui::Align2::RIGHT_BOTTOM,
            "x",
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(100, 116, 139),
        );

        let mut hovered_crosshair_bp: Option<(usize, usize)> = None;
        if response.hovered()
            && let Some(pointer) = response.hover_pos()
            && dotplot_rect.contains(pointer)
        {
            let fx = ((pointer.x - dotplot_rect.left()) / dotplot_rect.width()).clamp(0.0, 1.0);
            let fy = ((pointer.y - dotplot_rect.top()) / dotplot_rect.height()).clamp(0.0, 1.0);
            let x_bp = view.span_start_0based + (fx * query_span_max as f32).round() as usize;
            let y_bp = view.reference_span_start_0based
                + (fy * reference_span_max as f32).round() as usize;
            hovered_crosshair_bp = Some((x_bp, y_bp));
            let x_cell = ((fx * (cols - 1) as f32).round() as i32).clamp(0, cols - 1);
            let y_cell = ((fy * (rows - 1) as f32).round() as i32).clamp(0, rows - 1);
            let cell_payload = cells.get(&(x_cell, y_cell)).cloned().unwrap_or((0, 0));
            let click_hint = if !selection_sync_enabled {
                "Click to lock crosshair (sequence selection sync disabled for query override)"
            } else if is_self_mode {
                "Click to lock crosshair + sync sequence selection interval"
            } else {
                "Click to lock crosshair + sync query selection from x-axis"
            };
            response.clone().on_hover_ui_at_pointer(|ui| {
                ui.label(
                    egui::RichText::new(format!(
                        "x({})={} y({})={} | cell density={} min-mismatches={}",
                        view.seq_id,
                        x_bp.saturating_add(1),
                        reference_seq_label,
                        y_bp.saturating_add(1),
                        cell_payload.0,
                        cell_payload.1
                    ))
                    .monospace()
                    .size(self.feature_details_font_size()),
                );
                ui.small(click_hint);
            });
        }
        self.dotplot_hover_crosshair_bp = hovered_crosshair_bp;
        if response.clicked_by(egui::PointerButton::Primary)
            && let Some(pointer) = response.interact_pointer_pos()
            && dotplot_rect.contains(pointer)
        {
            let fx = ((pointer.x - dotplot_rect.left()) / dotplot_rect.width()).clamp(0.0, 1.0);
            let fy = ((pointer.y - dotplot_rect.top()) / dotplot_rect.height()).clamp(0.0, 1.0);
            let x_bp = view.span_start_0based + (fx * query_span_max as f32).round() as usize;
            let y_bp = view.reference_span_start_0based
                + (fy * reference_span_max as f32).round() as usize;
            self.dotplot_locked_crosshair_bp = Some((x_bp, y_bp));
            if selection_sync_enabled {
                self.sync_selection_to_dotplot_crosshair(x_bp, y_bp, true, view.mode);
            }
        }
        if response.clicked_by(egui::PointerButton::Secondary) {
            self.dotplot_locked_crosshair_bp = None;
            self.dotplot_hover_crosshair_bp = None;
        }

        let active_crosshair = self
            .dotplot_locked_crosshair_bp
            .or(self.dotplot_hover_crosshair_bp);
        if let Some((x_bp, y_bp)) = active_crosshair {
            let x_local = x_bp
                .saturating_sub(view.span_start_0based)
                .min(query_span_max) as f32;
            let y_local = y_bp
                .saturating_sub(view.reference_span_start_0based)
                .min(reference_span_max) as f32;
            let x = dotplot_rect.left() + (x_local / query_span_max as f32) * dotplot_rect.width();
            let y =
                dotplot_rect.top() + (y_local / reference_span_max as f32) * dotplot_rect.height();
            let is_locked = self.dotplot_locked_crosshair_bp == Some((x_bp, y_bp));
            let crosshair_color = if is_locked {
                egui::Color32::from_rgb(190, 24, 93)
            } else {
                egui::Color32::from_rgb(15, 23, 42)
            };
            painter.line_segment(
                [
                    egui::pos2(x, dotplot_rect.top()),
                    egui::pos2(x, dotplot_rect.bottom()),
                ],
                egui::Stroke::new(if is_locked { 1.6 } else { 1.1 }, crosshair_color),
            );
            painter.line_segment(
                [
                    egui::pos2(dotplot_rect.left(), y),
                    egui::pos2(dotplot_rect.right(), y),
                ],
                egui::Stroke::new(if is_locked { 1.6 } else { 1.1 }, crosshair_color),
            );
            painter.circle_filled(
                egui::pos2(x, y),
                if is_locked { 2.8 } else { 2.0 },
                crosshair_color,
            );
            let crosshair_text = if is_self_mode {
                let dx = x_bp.max(y_bp).saturating_sub(x_bp.min(y_bp));
                format!(
                    "{} x={} y={} Δ={} bp",
                    if is_locked { "locked" } else { "hover" },
                    x_bp.saturating_add(1),
                    y_bp.saturating_add(1),
                    dx
                )
            } else {
                format!(
                    "{} x({})={} y({})={}",
                    if is_locked { "locked" } else { "hover" },
                    view.seq_id,
                    x_bp.saturating_add(1),
                    reference_seq_label,
                    y_bp.saturating_add(1)
                )
            };
            painter.text(
                egui::pos2(dotplot_rect.left() + 4.0, dotplot_rect.top() + 14.0),
                egui::Align2::LEFT_TOP,
                crosshair_text,
                egui::FontId::monospace(10.0),
                crosshair_color,
            );
        }

        if let Some(track) = flex_track {
            let flex_rect = egui::Rect::from_min_max(
                egui::pos2(dotplot_rect.left(), dotplot_rect.bottom() + gap),
                egui::pos2(dotplot_rect.right(), canvas_rect.bottom() - bottom_margin),
            );
            painter.rect_filled(flex_rect, 4.0, egui::Color32::from_rgb(255, 255, 255));
            painter.rect_stroke(
                flex_rect,
                4.0,
                egui::Stroke::new(1.0, egui::Color32::from_rgb(203, 213, 225)),
                egui::StrokeKind::Inside,
            );
            let score_span = (track.max_score - track.min_score).abs().max(1e-12);
            let span_bp = track
                .span_end_0based
                .saturating_sub(track.span_start_0based)
                .max(1) as f32;
            let mut points: Vec<egui::Pos2> = vec![];
            for bin in &track.bins {
                let x_frac = (bin.start_0based.saturating_sub(track.span_start_0based) as f32
                    / span_bp)
                    .clamp(0.0, 1.0);
                let y_frac = ((bin.score - track.min_score) / score_span).clamp(0.0, 1.0) as f32;
                points.push(egui::pos2(
                    flex_rect.left() + x_frac * flex_rect.width(),
                    flex_rect.bottom() - y_frac * flex_rect.height(),
                ));
            }
            for segment in points.windows(2) {
                let [left, right] = [segment[0], segment[1]];
                painter.line_segment(
                    [left, right],
                    egui::Stroke::new(1.5, egui::Color32::from_rgb(22, 101, 52)),
                );
            }
            painter.text(
                egui::pos2(flex_rect.left() + 4.0, flex_rect.top() + 4.0),
                egui::Align2::LEFT_TOP,
                format!(
                    "{} [{}..{}] min={:.3} max={:.3}",
                    track.model.as_str(),
                    track.span_start_0based.saturating_add(1),
                    track.span_end_0based,
                    track.min_score,
                    track.max_score
                ),
                egui::FontId::monospace(9.5),
                egui::Color32::from_rgb(15, 23, 42),
            );
        }
    }

    pub(super) fn render_dotplot_boxplot_summary_ui(&self, ui: &mut egui::Ui, view: &DotplotView) {
        if view.boxplot_bins.is_empty() {
            ui.small(
                "No boxplot summary stored in this payload. Recompute dotplot to generate bin summaries.",
            );
            return;
        }
        ui.label(egui::RichText::new("Reference distribution by query bins (boxplot)").strong())
            .on_hover_text(
                "Per-query-bin distribution of reference hit positions: whiskers=min/max, box=Q1..Q3, center line=median.",
            );
        let desired_w = ui.available_width().max(740.0);
        let desired_h = 220.0;
        let (canvas_rect, response) =
            ui.allocate_exact_size(Vec2::new(desired_w, desired_h), egui::Sense::hover());
        let painter = ui.painter_at(canvas_rect);
        painter.rect_filled(canvas_rect, 6.0, egui::Color32::from_rgb(248, 250, 252));
        painter.rect_stroke(
            canvas_rect,
            6.0,
            egui::Stroke::new(1.0, egui::Color32::from_rgb(203, 213, 225)),
            egui::StrokeKind::Inside,
        );

        let top_margin = 22.0;
        let left_margin = 56.0;
        let right_margin = 16.0;
        let bottom_margin = 20.0;
        let plot_rect = egui::Rect::from_min_max(
            egui::pos2(
                canvas_rect.left() + left_margin,
                canvas_rect.top() + top_margin,
            ),
            egui::pos2(
                canvas_rect.right() - right_margin,
                canvas_rect.bottom() - bottom_margin,
            ),
        );
        painter.rect_filled(plot_rect, 4.0, egui::Color32::from_rgb(255, 255, 255));
        painter.rect_stroke(
            plot_rect,
            4.0,
            egui::Stroke::new(1.0, egui::Color32::from_rgb(203, 213, 225)),
            egui::StrokeKind::Inside,
        );
        for tick_idx in 1..10 {
            let frac = tick_idx as f32 / 10.0;
            let x = plot_rect.left() + frac * plot_rect.width();
            let y = plot_rect.top() + frac * plot_rect.height();
            painter.line_segment(
                [
                    egui::pos2(x, plot_rect.top()),
                    egui::pos2(x, plot_rect.bottom()),
                ],
                egui::Stroke::new(0.5, egui::Color32::from_rgb(241, 245, 249)),
            );
            painter.line_segment(
                [
                    egui::pos2(plot_rect.left(), y),
                    egui::pos2(plot_rect.right(), y),
                ],
                egui::Stroke::new(0.5, egui::Color32::from_rgb(241, 245, 249)),
            );
        }

        let query_span = view
            .span_end_0based
            .saturating_sub(view.span_start_0based)
            .max(1);
        let query_span_max = query_span.saturating_sub(1).max(1);
        let ref_span = view
            .reference_span_end_0based
            .saturating_sub(view.reference_span_start_0based)
            .max(1);
        let ref_span_max = ref_span.saturating_sub(1).max(1);
        let max_hits = view
            .boxplot_bins
            .iter()
            .map(|bin| bin.hit_count)
            .max()
            .unwrap_or(1)
            .max(1) as f32;

        let y_from_ref = |value_0based: usize| -> f32 {
            let local = value_0based
                .saturating_sub(view.reference_span_start_0based)
                .min(ref_span_max);
            plot_rect.top() + (local as f32 / ref_span_max as f32) * plot_rect.height()
        };
        for bin in &view.boxplot_bins {
            if bin.hit_count == 0 {
                continue;
            }
            let bin_start_local = bin
                .query_start_0based
                .saturating_sub(view.span_start_0based)
                .min(query_span_max);
            let bin_end_local = bin
                .query_end_0based_exclusive
                .saturating_sub(view.span_start_0based)
                .min(query_span);
            let x0 = plot_rect.left()
                + (bin_start_local as f32 / query_span_max as f32) * plot_rect.width();
            let x1_raw =
                plot_rect.left() + (bin_end_local as f32 / query_span as f32) * plot_rect.width();
            let x1 = (x1_raw.max(x0 + 1.0)).min(plot_rect.right());
            let x_center = (x0 + x1) * 0.5;
            let density = (bin.hit_count as f32 / max_hits).clamp(0.0, 1.0);
            let fill_alpha = (60.0 + 120.0 * density).round() as u8;
            if let (Some(q1), Some(q3)) = (bin.q1_reference_0based, bin.q3_reference_0based) {
                let y1 = y_from_ref(q1);
                let y3 = y_from_ref(q3);
                let y_top = y1.min(y3);
                let y_bottom = y1.max(y3).max(y_top + 1.0);
                painter.rect_filled(
                    egui::Rect::from_min_max(egui::pos2(x0, y_top), egui::pos2(x1, y_bottom)),
                    0.5,
                    egui::Color32::from_rgba_unmultiplied(37, 99, 235, fill_alpha),
                );
                painter.rect_stroke(
                    egui::Rect::from_min_max(egui::pos2(x0, y_top), egui::pos2(x1, y_bottom)),
                    0.5,
                    egui::Stroke::new(0.9, egui::Color32::from_rgb(29, 78, 216)),
                    egui::StrokeKind::Inside,
                );
            }
            if let (Some(min_ref), Some(max_ref)) =
                (bin.min_reference_0based, bin.max_reference_0based)
            {
                let y_min = y_from_ref(min_ref);
                let y_max = y_from_ref(max_ref);
                painter.line_segment(
                    [egui::pos2(x_center, y_min), egui::pos2(x_center, y_max)],
                    egui::Stroke::new(0.9, egui::Color32::from_rgb(15, 23, 42)),
                );
            }
            if let Some(median) = bin.median_reference_0based {
                let y_median = y_from_ref(median);
                painter.line_segment(
                    [egui::pos2(x0, y_median), egui::pos2(x1, y_median)],
                    egui::Stroke::new(1.2, egui::Color32::from_rgb(190, 24, 93)),
                );
            }
        }

        painter.text(
            egui::pos2(plot_rect.left(), canvas_rect.top() + 4.0),
            egui::Align2::LEFT_TOP,
            format!(
                "query bins={} (non-empty={})",
                view.boxplot_bin_count,
                view.boxplot_bins
                    .iter()
                    .filter(|bin| bin.hit_count > 0)
                    .count()
            ),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(51, 65, 85),
        );
        painter.text(
            egui::pos2(plot_rect.left(), plot_rect.bottom() + 2.0),
            egui::Align2::LEFT_TOP,
            format!("{}", view.span_start_0based.saturating_add(1)),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(71, 85, 105),
        );
        painter.text(
            egui::pos2(plot_rect.right(), plot_rect.bottom() + 2.0),
            egui::Align2::RIGHT_TOP,
            format!("{}", view.span_end_0based),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(71, 85, 105),
        );
        painter.text(
            egui::pos2(plot_rect.left() - 4.0, plot_rect.top()),
            egui::Align2::RIGHT_TOP,
            format!("{}", view.reference_span_start_0based.saturating_add(1)),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(71, 85, 105),
        );
        painter.text(
            egui::pos2(plot_rect.left() - 4.0, plot_rect.bottom()),
            egui::Align2::RIGHT_BOTTOM,
            format!("{}", view.reference_span_end_0based),
            egui::FontId::monospace(10.0),
            egui::Color32::from_rgb(71, 85, 105),
        );
        if response.hovered()
            && let Some(pointer) = response.hover_pos()
            && plot_rect.contains(pointer)
        {
            let fx = ((pointer.x - plot_rect.left()) / plot_rect.width()).clamp(0.0, 1.0);
            let query_0based =
                view.span_start_0based + (fx * query_span_max as f32).round() as usize;
            if let Some(bin) = view.boxplot_bins.iter().find(|bin| {
                query_0based >= bin.query_start_0based
                    && query_0based < bin.query_end_0based_exclusive
            }) {
                response.clone().on_hover_ui_at_pointer(|ui| {
                    ui.label(
                        egui::RichText::new(format!(
                            "query [{}..{}] hits={}",
                            bin.query_start_0based.saturating_add(1),
                            bin.query_end_0based_exclusive,
                            bin.hit_count
                        ))
                        .monospace()
                        .size(self.feature_details_font_size()),
                    );
                    let summary = if bin.hit_count == 0 {
                        "No hits in this bin".to_string()
                    } else {
                        format!(
                            "ref min={} q1={} median={} q3={} max={}",
                            bin.min_reference_0based
                                .map(|v| v.saturating_add(1).to_string())
                                .unwrap_or_else(|| "-".to_string()),
                            bin.q1_reference_0based
                                .map(|v| v.saturating_add(1).to_string())
                                .unwrap_or_else(|| "-".to_string()),
                            bin.median_reference_0based
                                .map(|v| v.saturating_add(1).to_string())
                                .unwrap_or_else(|| "-".to_string()),
                            bin.q3_reference_0based
                                .map(|v| v.saturating_add(1).to_string())
                                .unwrap_or_else(|| "-".to_string()),
                            bin.max_reference_0based
                                .map(|v| v.saturating_add(1).to_string())
                                .unwrap_or_else(|| "-".to_string()),
                        )
                    };
                    ui.label(
                        egui::RichText::new(summary)
                            .monospace()
                            .size(self.feature_details_font_size()),
                    );
                });
            }
        }
    }

    pub(super) fn render_primary_dotplot_map_ui(&mut self, ui: &mut egui::Ui) {
        let panel_width = ui.available_width().max(600.0);
        self.last_linear_map_width_px = panel_width;
        self.ensure_dotplot_cache_current();
        let mut save_state = false;
        let dotplot_ids = self.dotplot_ids_for_active_sequence();
        let track_ids = self.flexibility_track_ids_for_active_sequence();
        let reference_seq_ids = self.dotplot_reference_sequence_ids();
        if Self::dotplot_mode_requires_reference(self.dotplot_ui.mode)
            && self.dotplot_ui.reference_seq_id.trim().is_empty()
            && !reference_seq_ids.is_empty()
        {
            let current_id = self.seq_id.as_deref();
            let default_reference = reference_seq_ids
                .iter()
                .find(|id| Some(id.as_str()) != current_id)
                .cloned()
                .or_else(|| reference_seq_ids.first().cloned());
            if let Some(reference_id) = default_reference {
                self.dotplot_ui.reference_seq_id = reference_id;
                save_state = true;
            }
        }
        if self.normalize_dotplot_half_window_default_if_needed() {
            save_state = true;
        }

        egui::Frame::NONE
            .fill(egui::Color32::from_gray(249))
            .stroke(egui::Stroke::new(1.0, egui::Color32::from_gray(206)))
            .show(ui, |ui| {
                ui.set_width(panel_width);
                ui.vertical(|ui| {
                    ui.label(egui::RichText::new("Dotplot map").strong())
                        .on_hover_text(
                            "Compact dotplot launcher for this sequence. Open the dedicated Dotplot window for full controls and rendering.",
                        );
                    ui.small(
                        "This in-sequence panel is intentionally compact. Use 'Open Dotplot Window' for full parameter editing and plot inspection.",
                    );

                    ui.horizontal_wrapped(|ui| {
                        let open_label = if self.show_dotplot_window {
                            "Focus Dotplot Window"
                        } else {
                            "Open Dotplot Window"
                        };
                        if ui
                            .button(open_label)
                            .on_hover_text(
                                "Open the standalone dotplot workspace window for this sequence",
                            )
                            .clicked()
                        {
                            self.open_dotplot_window();
                        }
                        if ui
                            .button("Compute dotplot")
                            .on_hover_text(
                                "Run ComputeDotplot using current compact controls and store to dotplot_id",
                            )
                            .clicked()
                        {
                            self.compute_primary_dotplot();
                        }
                        if ui
                            .button("Compute flexibility")
                            .on_hover_text(
                                "Run ComputeFlexibilityTrack using current compact controls and store to flex_track_id",
                            )
                            .clicked()
                        {
                            self.compute_primary_flexibility_track();
                        }
                        if ui
                            .button("Export Dotplot SVG...")
                            .on_hover_text(
                                "Export the currently loaded dotplot density view as SVG. The default filename includes mode, spans, and compute/display parameters.",
                            )
                            .clicked()
                        {
                            self.export_active_dotplot_svg();
                        }
                    });

                    ui.horizontal_wrapped(|ui| {
                        ui.label("mode").on_hover_text(
                            "Dotplot comparison mode. Pair modes require reference_seq_id.",
                        );
                        let mode_before = self.dotplot_ui.mode;
                        let mode_combo = egui::ComboBox::from_id_salt(format!(
                            "dotplot_mode_compact_{}",
                            self.seq_id.as_deref().unwrap_or("<none>")
                        ))
                        .selected_text(self.dotplot_ui.mode.as_str())
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.dotplot_ui.mode,
                                DotplotMode::SelfForward,
                                DotplotMode::SelfForward.as_str(),
                            )
                            .on_hover_text("Compare sequence against itself in forward orientation");
                            ui.selectable_value(
                                &mut self.dotplot_ui.mode,
                                DotplotMode::SelfReverseComplement,
                                DotplotMode::SelfReverseComplement.as_str(),
                            )
                            .on_hover_text(
                                "Compare sequence against its reverse complement (inverted-repeat discovery)",
                            );
                            ui.selectable_value(
                                &mut self.dotplot_ui.mode,
                                DotplotMode::PairForward,
                                DotplotMode::PairForward.as_str(),
                            )
                            .on_hover_text("Compare active query sequence against reference sequence");
                            ui.selectable_value(
                                &mut self.dotplot_ui.mode,
                                DotplotMode::PairReverseComplement,
                                DotplotMode::PairReverseComplement.as_str(),
                            )
                            .on_hover_text(
                                "Compare active query sequence against reverse-complemented reference orientation",
                            );
                        });
                        mode_combo.response.on_hover_text(
                            "Choose self vs pair mode and forward vs reverse-complement matching",
                        );
                        if self.dotplot_ui.mode != mode_before {
                            save_state = true;
                        }

                        ui.label("half_window_bp").on_hover_text(
                            "Half-width of default bounded compute span around current viewport center (default fills to the larger sequence length)",
                        );
                        if Self::add_small_uint_text_edit(
                            ui,
                            &mut self.dotplot_ui.half_window_bp,
                            5,
                        )
                            .on_hover_text(
                                "Default auto-fills full-span view from the larger query/reference sequence length. Enter a value to bound around center.",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        ui.label("word").on_hover_text(
                            "k-mer seed length used for candidate matches (larger=faster, lower sensitivity)",
                        );
                        if Self::add_small_uint_text_edit(ui, &mut self.dotplot_ui.word_size, 2)
                            .on_hover_text("Minimum exact seed length in bases")
                            .changed()
                        {
                            save_state = true;
                        }
                        ui.label("step")
                            .on_hover_text("Sampling step in bp along query/reference spans");
                        if Self::add_small_uint_text_edit(ui, &mut self.dotplot_ui.step_bp, 4)
                            .on_hover_text(
                                "Higher step reduces compute cost but may skip small local structures",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        ui.label("mismatches")
                            .on_hover_text("Allowed mismatches per seed comparison");
                        if Self::add_small_uint_text_edit(
                            ui,
                            &mut self.dotplot_ui.max_mismatches,
                            2,
                        )
                            .on_hover_text(
                                "Default is 0 (exact seeds). Higher values increase tolerance.",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                    });

                    if Self::dotplot_mode_requires_reference(self.dotplot_ui.mode) {
                        ui.horizontal_wrapped(|ui| {
                            ui.label("reference_seq_id").on_hover_text(
                                "Reference sequence used on y-axis in pair modes",
                            );
                            if ui
                                .text_edit_singleline(&mut self.dotplot_ui.reference_seq_id)
                                .on_hover_text(
                                    "Sequence ID from project state; required in pair modes",
                                )
                                .changed()
                            {
                                self.maybe_normalize_dotplot_reference_from_current_input();
                                save_state = true;
                            }
                            if !reference_seq_ids.is_empty() {
                                let selected_reference =
                                    if self.dotplot_ui.reference_seq_id.trim().is_empty() {
                                        "<select>".to_string()
                                    } else {
                                        self.dotplot_ui.reference_seq_id.clone()
                                    };
                                let reference_combo = egui::ComboBox::from_id_salt(format!(
                                    "dotplot_reference_select_compact_{}",
                                    self.seq_id.as_deref().unwrap_or("<none>")
                                ))
                                .selected_text(selected_reference)
                                .show_ui(ui, |ui| {
                                    for id in &reference_seq_ids {
                                        if ui
                                            .selectable_label(
                                                self.dotplot_ui.reference_seq_id == *id,
                                                id,
                                            )
                                            .clicked()
                                        {
                                            let _ = self.set_dotplot_reference_to_sequence(id);
                                            save_state = true;
                                        }
                                    }
                                });
                                reference_combo.response.on_hover_text(
                                    "Choose reference sequence from currently loaded project sequences",
                                );
                            }
                            if ui
                                .small_button("Refit ref span")
                                .on_hover_text(
                                    "Manually refit hidden ref_start/ref_end to loaded hit envelope (+padding). Default behavior auto-fits once when ref_start/ref_end are empty.",
                                )
                                .clicked()
                                && self.fit_dotplot_reference_span_to_loaded_hits()
                            {
                                save_state = true;
                            }
                            self.render_dotplot_reference_quick_actions(
                                ui,
                                "compact",
                                &mut save_state,
                            );
                        });
                    }

                    ui.horizontal_wrapped(|ui| {
                        ui.label("dotplot_id").on_hover_text(
                            "Storage key for dotplot payload in project metadata",
                        );
                        if ui
                            .text_edit_singleline(&mut self.dotplot_ui.dotplot_id)
                            .on_hover_text(
                                "Set explicit ID to overwrite/reload deterministic dotplot payload",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        if !dotplot_ids.is_empty() {
                            let selected = if self.dotplot_ui.dotplot_id.trim().is_empty() {
                                "<select>".to_string()
                            } else {
                                self.dotplot_ui.dotplot_id.clone()
                            };
                            let dotplot_combo = egui::ComboBox::from_id_salt(format!(
                                "dotplot_select_compact_{}",
                                self.seq_id.as_deref().unwrap_or("<none>")
                            ))
                            .selected_text(selected)
                            .show_ui(ui, |ui| {
                                for id in &dotplot_ids {
                                    if ui
                                        .selectable_label(self.dotplot_ui.dotplot_id == *id, id)
                                        .clicked()
                                    {
                                        self.dotplot_ui.dotplot_id = id.clone();
                                        self.invalidate_dotplot_cache();
                                        save_state = true;
                                    }
                                }
                            });
                            dotplot_combo.response.on_hover_text(
                                "Pick from stored dotplot payload IDs for this query sequence",
                            );
                        }
                        if ui
                            .small_button("Load")
                            .on_hover_text("Load selected dotplot from metadata store")
                            .clicked()
                        {
                            self.invalidate_dotplot_cache();
                        }
                    });

                    ui.horizontal_wrapped(|ui| {
                        ui.checkbox(
                            &mut self.dotplot_ui.show_flexibility_track,
                            "Show flexibility track",
                        )
                        .on_hover_text(
                            "Show an auxiliary composition track below the dotplot. The track is computed from the sequence span and can help reveal local A/T-rich or skewed segments.",
                        );
                        ui.label("flex_track_id").on_hover_text(
                            "Persistent metadata key for flexibility payloads. Compute writes under this ID; Load reuses it later (GUI, shell, and exports).",
                        );
                        if ui
                            .text_edit_singleline(&mut self.dotplot_ui.flex_track_id)
                            .on_hover_text(
                                "Set explicit ID for flexibility track load/compute target",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        if !track_ids.is_empty() {
                            let selected = if self.dotplot_ui.flex_track_id.trim().is_empty() {
                                "<select>".to_string()
                            } else {
                                self.dotplot_ui.flex_track_id.clone()
                            };
                            let flex_combo = egui::ComboBox::from_id_salt(format!(
                                "flex_select_compact_{}",
                                self.seq_id.as_deref().unwrap_or("<none>")
                            ))
                            .selected_text(selected)
                            .show_ui(ui, |ui| {
                                for id in &track_ids {
                                    if ui
                                        .selectable_label(self.dotplot_ui.flex_track_id == *id, id)
                                        .clicked()
                                    {
                                        self.dotplot_ui.flex_track_id = id.clone();
                                        self.invalidate_dotplot_cache();
                                        save_state = true;
                                    }
                                }
                            });
                            flex_combo.response.on_hover_text(
                                "Pick from stored flexibility-track IDs for this query sequence",
                            );
                        }
                    });

                    ui.horizontal_wrapped(|ui| {
                        ui.label("display threshold").on_hover_text(
                            "Dotplot cell-density sensitivity threshold (higher hides faint/sparse cells)",
                        );
                        if ui
                            .add(
                                egui::Slider::new(
                                    &mut self.dotplot_ui.display_density_threshold,
                                    0.0..=0.95,
                                )
                                .show_value(true)
                                .trailing_fill(true),
                            )
                            .on_hover_text(
                                "Only render cells with normalized density >= threshold",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        ui.label("intensity gain").on_hover_text(
                            "Amplify visible density intensity (Dotter-style contrast control)",
                        );
                        if ui
                            .add(
                                egui::Slider::new(
                                    &mut self.dotplot_ui.display_intensity_gain,
                                    0.25..=8.0,
                                )
                                .logarithmic(true)
                                .show_value(true)
                                .trailing_fill(true),
                            )
                            .on_hover_text(
                                "Multiplier applied to normalized cell density before alpha encoding",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        if ui
                            .button("Auto contrast")
                            .on_hover_text(
                                "Estimate threshold/gain from currently loaded payload density distribution",
                            )
                            .clicked()
                        {
                            self.apply_auto_dotplot_display_contrast();
                        }
                    });

                    self.ensure_dotplot_cache_current();
                    match self.dotplot_cached_view.as_ref() {
                        Some(view) => {
                            let reference_seq_label =
                                view.reference_seq_id.as_deref().unwrap_or(view.seq_id.as_str());
                            ui.label(
                                egui::RichText::new(format!(
                                    "selected dotplot '{}' query={} [{}..{}] reference={} [{}..{}] | mode={} | points={}",
                                    view.dotplot_id,
                                    view.seq_id,
                                    view.span_start_0based.saturating_add(1),
                                    view.span_end_0based,
                                    reference_seq_label,
                                    view.reference_span_start_0based.saturating_add(1),
                                    view.reference_span_end_0based,
                                    view.mode.as_str(),
                                    view.point_count
                                ))
                                .monospace()
                                .size(self.feature_details_font_size()),
                            )
                            .on_hover_text(
                                "Current selected payload summary. Open Dotplot Window for full visualization and crosshair inspection.",
                            );
                        }
                        None => {
                            ui.small(
                                "No selected dotplot payload yet. Compute one here or open Dotplot Window for full controls.",
                            );
                        }
                    }
                    self.render_dotplot_status_ui(ui);
                });
            });

        if save_state {
            self.save_engine_ops_state();
        }
    }

    pub(super) fn render_dotplot_workspace_ui(&mut self, ui: &mut egui::Ui) {
        let panel_width = ui.available_width().max(600.0);
        self.last_linear_map_width_px = panel_width;
        self.ensure_dotplot_cache_current();
        let mut save_state = false;
        let mut compute_inputs_changed = false;
        let dotplot_ids = self.dotplot_ids_for_active_sequence();
        let track_ids = self.flexibility_track_ids_for_active_sequence();
        let reference_seq_ids = self.dotplot_reference_sequence_ids();
        if Self::dotplot_mode_requires_reference(self.dotplot_ui.mode)
            && self.dotplot_ui.reference_seq_id.trim().is_empty()
            && !reference_seq_ids.is_empty()
        {
            let current_id = self.seq_id.as_deref();
            let default_reference = reference_seq_ids
                .iter()
                .find(|id| Some(id.as_str()) != current_id)
                .cloned()
                .or_else(|| reference_seq_ids.first().cloned());
            if let Some(reference_id) = default_reference {
                self.dotplot_ui.reference_seq_id = reference_id;
                save_state = true;
            }
        }
        if self.normalize_dotplot_half_window_default_if_needed() {
            save_state = true;
            compute_inputs_changed = true;
        }
        if self.sync_dotplot_overlay_selection() {
            save_state = true;
        }
        if self.sync_dotplot_overlay_anchor_selection() {
            save_state = true;
        }

        egui::Frame::NONE
            .fill(egui::Color32::from_gray(249))
            .stroke(egui::Stroke::new(1.0, egui::Color32::from_gray(206)))
            .show(ui, |ui| {
                ui.set_width(panel_width);
                ui.vertical(|ui| {
                    ui.label(egui::RichText::new("Dotplot workspace").strong()).on_hover_text(
                        "Promoter-oriented dotplot (self or pairwise) and flexibility track view for a bounded local span",
                    );
                    ui.label(
                        egui::RichText::new(
                            "Default compute span auto-fills to full view using the larger query/reference sequence length.",
                        )
                        .size(self.feature_details_font_size()),
                    );

                    ui.horizontal_wrapped(|ui| {
                        ui.label("Mode").on_hover_text(
                            "Dotplot comparison mode. Pair modes require a reference sequence.",
                        );
                        let mode_before = self.dotplot_ui.mode;
                        let mode_combo = egui::ComboBox::from_id_salt(format!(
                            "dotplot_mode_{}",
                            self.seq_id.as_deref().unwrap_or("<none>")
                        ))
                        .selected_text(self.dotplot_ui.mode.as_str())
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.dotplot_ui.mode,
                                DotplotMode::SelfForward,
                                DotplotMode::SelfForward.as_str(),
                            )
                            .on_hover_text("Self-vs-self comparison in forward orientation");
                            ui.selectable_value(
                                &mut self.dotplot_ui.mode,
                                DotplotMode::SelfReverseComplement,
                                DotplotMode::SelfReverseComplement.as_str(),
                            )
                            .on_hover_text("Self-vs-reverse-complement (inverted repeat detection)");
                            ui.selectable_value(
                                &mut self.dotplot_ui.mode,
                                DotplotMode::PairForward,
                                DotplotMode::PairForward.as_str(),
                            )
                            .on_hover_text("Query sequence versus reference sequence");
                            ui.selectable_value(
                                &mut self.dotplot_ui.mode,
                                DotplotMode::PairReverseComplement,
                                DotplotMode::PairReverseComplement.as_str(),
                            )
                            .on_hover_text(
                                "Query sequence versus reverse-complemented reference orientation",
                            );
                        });
                        mode_combo.response.on_hover_text(
                            "Select self/pair and forward/reverse-complement comparison behavior",
                        );
                        if self.dotplot_ui.mode != mode_before {
                            if !Self::dotplot_mode_requires_reference(self.dotplot_ui.mode) {
                                self.dotplot_ui.overlay_enabled = false;
                            }
                            save_state = true;
                            compute_inputs_changed = true;
                        }
                        ui.label("half_window_bp").on_hover_text(
                            "Half-width of default query span around current viewport center (default fills to larger query/reference sequence length)",
                        );
                        if Self::add_small_uint_text_edit(
                            ui,
                            &mut self.dotplot_ui.half_window_bp,
                            5,
                        )
                            .on_hover_text(
                                "Default auto-fills full span from larger sequence length; enter explicit value to bound around viewport center",
                            )
                            .changed()
                        {
                            save_state = true;
                            compute_inputs_changed = true;
                        }
                        ui.label("word")
                            .on_hover_text("Seed/k-mer length used for candidate matching");
                        if Self::add_small_uint_text_edit(ui, &mut self.dotplot_ui.word_size, 2)
                            .on_hover_text("Higher values are faster but less sensitive")
                            .changed()
                        {
                            save_state = true;
                            compute_inputs_changed = true;
                        }
                        ui.label("step")
                            .on_hover_text("Sampling step (bp) along query/reference spans");
                        if Self::add_small_uint_text_edit(ui, &mut self.dotplot_ui.step_bp, 4)
                            .on_hover_text("Higher values reduce work but can miss local structures")
                            .changed()
                        {
                            save_state = true;
                            compute_inputs_changed = true;
                        }
                        ui.label("mismatches")
                            .on_hover_text("Allowed mismatches per seed comparison");
                        if Self::add_small_uint_text_edit(
                            ui,
                            &mut self.dotplot_ui.max_mismatches,
                            2,
                        )
                            .on_hover_text("Default 0 = exact seed matches only")
                            .changed()
                        {
                            save_state = true;
                            compute_inputs_changed = true;
                        }
                        ui.label("tile_bp").on_hover_text(
                            "Optional rendering/export hint for tile width in bp (analysis metadata)",
                        );
                        if Self::add_small_uint_text_edit(ui, &mut self.dotplot_ui.tile_bp, 5)
                            .on_hover_text("Leave empty for automatic/default behavior")
                            .changed()
                        {
                            save_state = true;
                            compute_inputs_changed = true;
                        }
                        if ui
                            .button("Compute dotplot")
                            .on_hover_text("Run ComputeDotplot with current parameters")
                            .clicked()
                        {
                            self.compute_primary_dotplot();
                        }
                        if ui
                            .button("Export Dotplot SVG...")
                            .on_hover_text(
                                "Export the currently loaded dotplot density view as SVG. The default filename includes mode, spans, and compute/display parameters.",
                            )
                            .clicked()
                        {
                            self.export_active_dotplot_svg();
                        }
                    });

                    ui.horizontal_wrapped(|ui| {
                        ui.label("dotplot_id")
                            .on_hover_text("Storage key for dotplot payload in metadata");
                        if ui
                            .text_edit_singleline(&mut self.dotplot_ui.dotplot_id)
                            .on_hover_text(
                                "Use stable IDs to overwrite/load deterministic dotplot payloads",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        if !dotplot_ids.is_empty() {
                            let selected = if self.dotplot_ui.dotplot_id.trim().is_empty() {
                                "<select>".to_string()
                            } else {
                                self.dotplot_ui.dotplot_id.clone()
                            };
                            let dotplot_combo = egui::ComboBox::from_id_salt(format!(
                                "dotplot_select_{}",
                                self.seq_id.as_deref().unwrap_or("<none>")
                            ))
                            .selected_text(selected)
                            .show_ui(ui, |ui| {
                                for id in &dotplot_ids {
                                    if ui
                                        .selectable_label(self.dotplot_ui.dotplot_id == *id, id)
                                        .clicked()
                                    {
                                        self.dotplot_ui.dotplot_id = id.clone();
                                        self.invalidate_dotplot_cache();
                                        save_state = true;
                                    }
                                }
                            });
                            dotplot_combo.response.on_hover_text(
                                "Pick a stored dotplot payload for this query sequence",
                            );
                        }
                        if ui
                            .small_button("Load")
                            .on_hover_text("Load selected dotplot_id from engine metadata")
                            .clicked()
                        {
                            self.invalidate_dotplot_cache();
                        }
                    });
                    if Self::dotplot_mode_requires_reference(self.dotplot_ui.mode) {
                        ui.horizontal_wrapped(|ui| {
                            ui.label("reference_seq_id").on_hover_text(
                                "Reference sequence ID used for y-axis in pair modes",
                            );
                            if ui
                                .text_edit_singleline(&mut self.dotplot_ui.reference_seq_id)
                                .on_hover_text(
                                    "Required in pair modes; pick from currently loaded sequences",
                                )
                                .changed()
                            {
                                self.maybe_normalize_dotplot_reference_from_current_input();
                                save_state = true;
                                compute_inputs_changed = true;
                            }
                            if !reference_seq_ids.is_empty() {
                                let selected_reference =
                                    if self.dotplot_ui.reference_seq_id.trim().is_empty() {
                                        "<select>".to_string()
                                    } else {
                                        self.dotplot_ui.reference_seq_id.clone()
                                    };
                                let reference_combo = egui::ComboBox::from_id_salt(format!(
                                    "dotplot_reference_select_{}",
                                    self.seq_id.as_deref().unwrap_or("<none>")
                                ))
                                .selected_text(selected_reference)
                                .show_ui(ui, |ui| {
                                    for id in &reference_seq_ids {
                                        if ui
                                            .selectable_label(
                                                self.dotplot_ui.reference_seq_id == *id,
                                                id,
                                            )
                                            .clicked()
                                        {
                                            let _ = self.set_dotplot_reference_to_sequence(id);
                                            save_state = true;
                                            compute_inputs_changed = true;
                                        }
                                    }
                                });
                                reference_combo.response.on_hover_text(
                                    "Choose reference sequence from project sequence IDs",
                                );
                            }
                            ui.label("ref_start")
                                .on_hover_text("Optional 0-based inclusive reference span start");
                            if Self::add_small_uint_text_edit(
                                ui,
                                &mut self.dotplot_ui.reference_span_start_0based,
                                5,
                            )
                                .on_hover_text(
                                    "Leave empty for default auto-fit flow: initial full reference pass, then automatic hit-envelope fitting and recompute.",
                                )
                                .changed()
                            {
                                save_state = true;
                                compute_inputs_changed = true;
                            }
                            ui.label("ref_end")
                                .on_hover_text("Optional 0-based exclusive reference span end");
                            if Self::add_small_uint_text_edit(
                                ui,
                                &mut self.dotplot_ui.reference_span_end_0based,
                                5,
                            )
                                .on_hover_text(
                                    "Leave empty for default auto-fit flow: initial full reference pass, then automatic hit-envelope fitting and recompute.",
                                )
                                .changed()
                            {
                                save_state = true;
                                compute_inputs_changed = true;
                            }
                            if ui
                                .small_button("Fit ref span to hits")
                                .on_hover_text(
                                    "Manually refit ref_start/ref_end to loaded hit envelope (+padding). Default behavior auto-fits once when both fields are empty.",
                                )
                                .clicked()
                                && self.fit_dotplot_reference_span_to_loaded_hits()
                            {
                                save_state = true;
                            }
                            self.render_dotplot_reference_quick_actions(
                                ui,
                                "full",
                                &mut save_state,
                            );
                            ui.small("Leave ref_start/ref_end empty for automatic first-pass full-span compute followed by automatic hit-envelope refit.");
                        });
                        let overlay_choices = self.dotplot_overlay_choices();
                        ui.horizontal_wrapped(|ui| {
                            let overlay_enabled = !overlay_choices.is_empty();
                            if ui
                                .add_enabled(
                                    overlay_enabled,
                                    egui::Checkbox::new(
                                        &mut self.dotplot_ui.overlay_enabled,
                                        "Overlay transcript isoforms",
                                    ),
                                )
                                .on_hover_text(
                                    "Reference-centered overlay: compare multiple annotated transcript isoforms against the same reference sequence in one plot.",
                                )
                                .changed()
                            {
                                if self.dotplot_ui.overlay_enabled
                                    && self.dotplot_ui.overlay_transcript_feature_ids.is_empty()
                                {
                                    self.dotplot_ui.overlay_transcript_feature_ids = overlay_choices
                                        .iter()
                                        .map(|choice| choice.transcript_feature_id)
                                        .collect();
                                }
                                save_state = true;
                                compute_inputs_changed = true;
                            }
                            if !overlay_enabled {
                                ui.small(
                                    if self.using_dotplot_query_override() {
                                        "Overlay is unavailable while the workspace is bound to an RNA-read query override."
                                    } else {
                                        "Overlay becomes available when the current reference is the active locus DNA with transcript annotation."
                                    },
                                );
                            } else {
                                ui.small(
                                    "Selected isoforms reuse one shared reference span; colors distinguish the overlaid transcript queries.",
                                );
                            }
                        });
                        let loaded_overlay_view = self
                            .dotplot_cached_view
                            .as_ref()
                            .map(Self::dotplot_view_is_overlay)
                            .unwrap_or(false);
                        if (self.dotplot_ui.overlay_enabled && !overlay_choices.is_empty())
                            || loaded_overlay_view
                        {
                            let overlay_anchor_choices = self.dotplot_overlay_anchor_choices();
                            ui.horizontal_wrapped(|ui| {
                                ui.label("overlay x-axis").on_hover_text(
                                    "Choose how overlaid transcript queries share the x-axis: normalize each transcript to percent length, align by one end in bp, anchor one exon shared across transcripts, or align explicit manual/domain anchors carried in the overlay payload.",
                                );
                                let before = self.dotplot_ui.overlay_x_axis_mode;
                                egui::ComboBox::from_id_salt(format!(
                                    "dotplot_overlay_x_axis_{}",
                                    self.seq_id.as_deref().unwrap_or("<none>")
                                ))
                                .selected_text(self.dotplot_ui.overlay_x_axis_mode.ui_label())
                                .show_ui(ui, |ui| {
                                    ui.selectable_value(
                                        &mut self.dotplot_ui.overlay_x_axis_mode,
                                        DotplotOverlayXAxisMode::PercentLength,
                                        DotplotOverlayXAxisMode::PercentLength.ui_label(),
                                    )
                                    .on_hover_text(
                                        "Normalize each transcript independently from 0% to 100%. This is best for comparing shared relative exon structure.",
                                    );
                                    ui.selectable_value(
                                        &mut self.dotplot_ui.overlay_x_axis_mode,
                                        DotplotOverlayXAxisMode::LeftAlignedBp,
                                        DotplotOverlayXAxisMode::LeftAlignedBp.ui_label(),
                                    )
                                    .on_hover_text(
                                        "Use base-pair coordinates and align transcript starts. This makes different 3' splice endings stand out on the right.",
                                    );
                                    ui.selectable_value(
                                        &mut self.dotplot_ui.overlay_x_axis_mode,
                                        DotplotOverlayXAxisMode::RightAlignedBp,
                                        DotplotOverlayXAxisMode::RightAlignedBp.ui_label(),
                                    )
                                    .on_hover_text(
                                        "Use base-pair coordinates and align transcript ends. This makes different 5' splice starts stand out on the left.",
                                    );
                                    ui.selectable_value(
                                        &mut self.dotplot_ui.overlay_x_axis_mode,
                                        DotplotOverlayXAxisMode::SharedExonAnchor,
                                        DotplotOverlayXAxisMode::SharedExonAnchor.ui_label(),
                                    )
                                    .on_hover_text(
                                        "Choose one exon shared by multiple selected transcripts and align that exon start to one common x position.",
                                    );
                                    ui.selectable_value(
                                        &mut self.dotplot_ui.overlay_x_axis_mode,
                                        DotplotOverlayXAxisMode::QueryAnchorBp,
                                        DotplotOverlayXAxisMode::QueryAnchorBp.ui_label(),
                                    )
                                    .on_hover_text(
                                        "Use explicit per-query anchors stored in the overlay payload. This is intended for cross-family or domain-anchored comparisons where the plotted series do not share one genomic exon identity.",
                                    );
                                });
                                if self.dotplot_ui.overlay_x_axis_mode != before {
                                    save_state = true;
                                }
                                ui.small(
                                    "This is a display/export choice only; it does not recompute the dotplot seed matches.",
                                );
                            });
                            if self.dotplot_ui.overlay_x_axis_mode
                                == DotplotOverlayXAxisMode::SharedExonAnchor
                            {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label("anchor exon").on_hover_text(
                                        "Select the genomic exon that should land at one shared x position across all plotted transcripts that contain it.",
                                    );
                                    let selected_text = self
                                        .dotplot_ui
                                        .overlay_anchor_exon
                                        .as_ref()
                                        .map(|exon| exon.token())
                                        .unwrap_or_else(|| "<select exon>".to_string());
                                    egui::ComboBox::from_id_salt(format!(
                                        "dotplot_overlay_anchor_exon_{}",
                                        self.seq_id.as_deref().unwrap_or("<none>")
                                    ))
                                    .selected_text(selected_text)
                                    .show_ui(ui, |ui| {
                                        for (exon, support_count) in &overlay_anchor_choices {
                                            let label = format!(
                                                "{} ({} transcripts)",
                                                exon.token(),
                                                support_count
                                            );
                                            if ui
                                                .selectable_label(
                                                    self.dotplot_ui
                                                        .overlay_anchor_exon
                                                        .as_ref()
                                                        == Some(exon),
                                                    label,
                                                )
                                                .on_hover_text(
                                                    "Only selected transcripts containing this exon remain visible in the anchored overlay.",
                                                )
                                                .clicked()
                                            {
                                                self.dotplot_ui.overlay_anchor_exon =
                                                    Some(exon.clone());
                                                save_state = true;
                                            }
                                        }
                                    });
                                    if ui
                                        .small_button("Clear")
                                        .on_hover_text("Clear the shared-exon anchor selection")
                                        .clicked()
                                    {
                                        self.dotplot_ui.overlay_anchor_exon = None;
                                        save_state = true;
                                    }
                                    if overlay_anchor_choices.is_empty() {
                                        ui.small(
                                            "No exon is shared by at least two of the currently selected overlay transcripts.",
                                        );
                                    } else if self.dotplot_ui.overlay_anchor_exon.is_none() {
                                        ui.small(
                                            "Choose one shared exon to activate exon-anchored overlay rendering.",
                                        );
                                    }
                                });
                            }
                            if self.dotplot_ui.overlay_enabled && !overlay_choices.is_empty() {
                                ui.horizontal_wrapped(|ui| {
                                    for choice in &overlay_choices {
                                        let mut selected = self
                                            .dotplot_ui
                                            .overlay_transcript_feature_ids
                                            .contains(&choice.transcript_feature_id);
                                        let swatch = egui::Color32::from_rgb(
                                            choice.color_rgb[0],
                                            choice.color_rgb[1],
                                            choice.color_rgb[2],
                                        );
                                        let label = format!(
                                            "{} ({} bp)",
                                            choice.label, choice.estimated_length_bp
                                        );
                                        if ui
                                            .checkbox(
                                                &mut selected,
                                                egui::RichText::new(label).color(swatch),
                                            )
                                            .on_hover_text(format!(
                                                "Use transcript '{}' as an overlaid query against the shared reference span.",
                                                choice.transcript_id
                                            ))
                                            .changed()
                                        {
                                            if selected {
                                                self.dotplot_ui
                                                    .overlay_transcript_feature_ids
                                                    .push(choice.transcript_feature_id);
                                            } else {
                                                self.dotplot_ui
                                                    .overlay_transcript_feature_ids
                                                    .retain(|feature_id| {
                                                        *feature_id != choice.transcript_feature_id
                                                    });
                                            }
                                            self.dotplot_ui
                                                .overlay_transcript_feature_ids
                                                .sort_unstable();
                                            self.dotplot_ui
                                                .overlay_transcript_feature_ids
                                                .dedup();
                                            save_state = true;
                                            compute_inputs_changed = true;
                                        }
                                    }
                                });
                            }
                        }
                    }

                    ui.horizontal_wrapped(|ui| {
                        if ui
                            .checkbox(
                                &mut self.dotplot_ui.show_flexibility_track,
                                "Show flexibility track",
                            )
                            .on_hover_text(
                                "Render an auxiliary composition profile below the dotplot. This is visual context only; it does not modify dotplot matching.",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        ui.label("model")
                            .on_hover_text(
                                "Scoring model for flexibility profile: at_richness = (A+T)/len, at_skew = (A-T)/(A+T).",
                            );
                        let flex_model_before = self.dotplot_ui.flex_model;
                        let flex_model_combo = egui::ComboBox::from_id_salt(format!(
                            "flex_model_{}",
                            self.seq_id.as_deref().unwrap_or("<none>")
                        ))
                        .selected_text(self.dotplot_ui.flex_model.as_str())
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.dotplot_ui.flex_model,
                                FlexibilityModel::AtRichness,
                                FlexibilityModel::AtRichness.as_str(),
                            )
                            .on_hover_text(
                                "AT-richness per bin: fraction of A/T bases in each bin (0..1)",
                            );
                            ui.selectable_value(
                                &mut self.dotplot_ui.flex_model,
                                FlexibilityModel::AtSkew,
                                FlexibilityModel::AtSkew.as_str(),
                            )
                            .on_hover_text(
                                "AT-skew per bin: (A-T)/(A+T), highlighting strand composition bias",
                            );
                        });
                        flex_model_combo
                            .response
                            .on_hover_text("Select flexibility scoring model");
                        if self.dotplot_ui.flex_model != flex_model_before {
                            save_state = true;
                        }
                        ui.label("bin_bp")
                            .on_hover_text(
                                "Bin width for flexibility scoring. Smaller bins increase detail, larger bins emphasize broad trends.",
                            );
                        if Self::add_small_uint_text_edit(ui, &mut self.dotplot_ui.flex_bin_bp, 4)
                            .on_hover_text("Smaller bins give finer but noisier track detail")
                            .changed()
                        {
                            save_state = true;
                        }
                        ui.label("smoothing_bp")
                            .on_hover_text(
                                "Optional moving-average smoothing window (bp). Leave empty for raw per-bin profile.",
                            );
                        if Self::add_small_uint_text_edit(
                            ui,
                            &mut self.dotplot_ui.flex_smoothing_bp,
                            5,
                        )
                            .on_hover_text("Leave empty for no smoothing")
                            .changed()
                        {
                            save_state = true;
                        }
                        if ui
                            .button("Compute flexibility")
                            .on_hover_text(
                                "Run ComputeFlexibilityTrack on bounded local span via engine",
                            )
                            .clicked()
                        {
                            self.compute_primary_flexibility_track();
                        }
                    });

                    ui.horizontal_wrapped(|ui| {
                        ui.label("flex_track_id")
                            .on_hover_text(
                                "Persistent metadata key for flexibility payload. Compute writes under this ID; Load reuses the same payload later in GUI/shell/export flows.",
                            );
                        if ui
                            .text_edit_singleline(&mut self.dotplot_ui.flex_track_id)
                            .on_hover_text(
                                "Use stable IDs for reproducible workflows and later retrieval (for example before SVG export or in CLI review).",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        if !track_ids.is_empty() {
                            let selected = if self.dotplot_ui.flex_track_id.trim().is_empty() {
                                "<select>".to_string()
                            } else {
                                self.dotplot_ui.flex_track_id.clone()
                            };
                            let flex_combo = egui::ComboBox::from_id_salt(format!(
                                "flex_select_{}",
                                self.seq_id.as_deref().unwrap_or("<none>")
                            ))
                            .selected_text(selected)
                            .show_ui(ui, |ui| {
                                for id in &track_ids {
                                    if ui
                                        .selectable_label(self.dotplot_ui.flex_track_id == *id, id)
                                        .clicked()
                                    {
                                        self.dotplot_ui.flex_track_id = id.clone();
                                        self.invalidate_dotplot_cache();
                                        save_state = true;
                                    }
                                }
                            });
                            flex_combo.response.on_hover_text(
                                "Pick a stored flexibility payload for this query sequence",
                            );
                        }
                        if ui
                            .small_button("Load")
                            .on_hover_text("Load selected flexibility track from engine metadata")
                            .clicked()
                        {
                            self.invalidate_dotplot_cache();
                        }
                    });

                    ui.horizontal_wrapped(|ui| {
                        ui.label("display threshold").on_hover_text(
                            "Dotplot cell-density sensitivity threshold (hide sparse/faint cells below this)",
                        );
                        if ui
                            .add(
                                egui::Slider::new(
                                    &mut self.dotplot_ui.display_density_threshold,
                                    0.0..=0.95,
                                )
                                .show_value(true)
                                .trailing_fill(true),
                            )
                            .on_hover_text(
                                "Only cells with normalized density >= threshold are rendered",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        ui.label("intensity gain").on_hover_text(
                            "Amplify visible dot intensity after thresholding (Dotter-style contrast)",
                        );
                        if ui
                            .add(
                                egui::Slider::new(
                                    &mut self.dotplot_ui.display_intensity_gain,
                                    0.25..=8.0,
                                )
                                .logarithmic(true)
                                .show_value(true)
                                .trailing_fill(true),
                            )
                            .on_hover_text(
                                "Gain multiplier applied to normalized density before alpha mapping",
                            )
                            .changed()
                        {
                            save_state = true;
                        }
                        if ui
                            .button("Auto contrast")
                            .on_hover_text(
                                "Estimate threshold/gain from currently loaded payload density distribution",
                            )
                            .clicked()
                        {
                            self.apply_auto_dotplot_display_contrast();
                        }
                    });

                    ui.horizontal_wrapped(|ui| {
                        if self.using_dotplot_query_override() {
                            ui.label(
                                egui::RichText::new(format!(
                                    "Query override: {}",
                                    self.current_dotplot_query_label()
                                ))
                                .monospace()
                                .size(self.feature_details_font_size()),
                            )
                            .on_hover_text(
                                "This workspace is currently bound to a materialized RNA-read query sequence instead of the active genomic sequence window.",
                            );
                            if ui
                                .small_button("Use active sequence as query")
                                .on_hover_text(
                                    "Clear the RNA-read query override and return the workspace to the active sequence window.",
                                )
                                .clicked()
                            {
                                self.clear_dotplot_query_override();
                                save_state = true;
                            }
                            ui.separator();
                        }
                        let overlay_loaded = self
                            .dotplot_cached_view
                            .as_ref()
                            .map(Self::dotplot_view_is_overlay)
                            .unwrap_or(false);
                        if overlay_loaded {
                            ui.small(
                                "Overlay hover reports one shared reference coordinate plus one x/query coordinate per isoform using the selected overlay x-axis layout. Query-side selection sync and locked crosshair are disabled in overlay mode.",
                            );
                        } else if let Some((x_bp, y_bp)) = self.dotplot_locked_crosshair_bp {
                            let active_mode = self
                                .dotplot_cached_view
                                .as_ref()
                                .map(|view| view.mode)
                                .unwrap_or(self.dotplot_ui.mode);
                            let reference_label = self
                                .dotplot_cached_view
                                .as_ref()
                                .and_then(|view| view.reference_seq_id.clone())
                                .unwrap_or_else(|| {
                                    self.seq_id
                                        .clone()
                                        .unwrap_or_else(|| "reference".to_string())
                                });
                            ui.label(
                                egui::RichText::new(if Self::dotplot_mode_requires_reference(
                                    active_mode,
                                ) {
                                    format!(
                                        "locked crosshair x(query)={} y(reference:{})={}",
                                        x_bp.saturating_add(1),
                                        reference_label,
                                        y_bp.saturating_add(1)
                                    )
                                } else {
                                    format!(
                                        "locked crosshair x={} y={} Δ={} bp",
                                        x_bp.saturating_add(1),
                                        y_bp.saturating_add(1),
                                        x_bp.max(y_bp).saturating_sub(x_bp.min(y_bp))
                                    )
                                })
                                .monospace()
                                .size(self.feature_details_font_size()),
                            );
                            if ui
                                .add_enabled(
                                    self.dotplot_cached_view
                                        .as_ref()
                                        .map(|view| self.dotplot_selection_sync_enabled_for_view(view))
                                        .unwrap_or(!self.using_dotplot_query_override()),
                                    egui::Button::new("Re-sync selection"),
                                )
                                .on_hover_text(
                                    if self.using_dotplot_query_override() {
                                        "Selection sync is disabled while the dotplot query is an RNA-read override rather than the active genomic sequence."
                                    } else {
                                        "Apply locked crosshair selection to shared sequence selection (pair mode uses query/x axis)"
                                    },
                                )
                                .clicked()
                            {
                                self.sync_selection_to_dotplot_crosshair(
                                    x_bp,
                                    y_bp,
                                    true,
                                    active_mode,
                                );
                            }
                            if ui
                                .small_button("Clear crosshair")
                                .on_hover_text("Clear locked crosshair anchor")
                                .clicked()
                            {
                                self.dotplot_locked_crosshair_bp = None;
                            }
                        } else {
                            ui.small(
                                if self.using_dotplot_query_override() {
                                    "Hover inside dotplot for live crosshair. Click to lock the crosshair; sequence-selection sync is disabled for RNA-read query overrides. Right-click to clear."
                                } else {
                                    "Hover inside dotplot for live crosshair. Click to lock and sync selection from the x/query axis. Right-click to clear."
                                },
                            );
                        }
                    });

                    self.ensure_dotplot_cache_current();
                    let cached_view = self.dotplot_cached_view.clone();
                    let cached_track = self.dotplot_cached_flex_track.clone();
                    match cached_view.as_ref() {
                        Some(view) => {
                            let reference_seq_label =
                                view.reference_seq_id.as_deref().unwrap_or(view.seq_id.as_str());
                            let total_point_count = Self::dotplot_view_total_point_count(view);
                            ui.label(
                                egui::RichText::new(format!(
                                    "dotplot '{}' owner={} primary_query={} [{}..{}] reference={} [{}..{}] | series={}{} | mode={} | word={} step={} mismatches={} | points={}",
                                    view.dotplot_id,
                                    view.owner_seq_id,
                                    view.seq_id,
                                    view.span_start_0based.saturating_add(1),
                                    view.span_end_0based,
                                    reference_seq_label,
                                    view.reference_span_start_0based.saturating_add(1),
                                    view.reference_span_end_0based,
                                    view.series_count.max(view.query_series.len()),
                                    if Self::dotplot_view_is_overlay(view) {
                                        format!(
                                            " overlay_x={}",
                                            self.dotplot_ui.overlay_x_axis_mode.as_str()
                                        )
                                    } else {
                                        String::new()
                                    },
                                    view.mode.as_str(),
                                    view.word_size,
                                    view.step_bp,
                                    view.max_mismatches,
                                    total_point_count
                                ))
                                .monospace()
                                .size(self.feature_details_font_size()),
                            );
                            let overlay_mode = Self::dotplot_view_is_overlay(view);
                            let flex_track = if self.dotplot_ui.show_flexibility_track && !overlay_mode {
                                cached_track.as_ref()
                            } else {
                                None
                            };
                            if overlay_mode && self.dotplot_ui.show_flexibility_track {
                                ui.small(
                                    "Flexibility track is hidden in overlay mode because the x-axis can represent multiple transcript coordinate systems rather than one shared query span.",
                                );
                            }
                            self.render_dotplot_density_ui(
                                ui,
                                view,
                                flex_track,
                                self.dotplot_ui.display_density_threshold,
                                self.dotplot_ui.display_intensity_gain,
                            );
                            ui.add_space(6.0);
                            if !overlay_mode {
                                self.render_dotplot_boxplot_summary_ui(ui, view);
                            }
                        }
                        None => {
                            ui.add_space(6.0);
                            ui.label(
                                "No dotplot payload selected. Compute one with current bounded span, or select an existing dotplot_id.",
                            );
                        }
                    }
                    self.render_dotplot_status_ui(ui);
                });
            });
        if save_state {
            self.save_engine_ops_state();
        }
        if compute_inputs_changed {
            self.schedule_dotplot_auto_compute();
        }
        self.maybe_auto_compute_dotplot();
    }

    pub(super) fn render_primary_splicing_map_ui(&mut self, ui: &mut egui::Ui) {
        let panel_width = ui.available_width().max(600.0);
        self.last_linear_map_width_px = panel_width;
        egui::Frame::NONE
            .fill(egui::Color32::from_gray(249))
            .stroke(egui::Stroke::new(1.0, egui::Color32::from_gray(206)))
            .show(ui, |ui| {
                ui.set_width(panel_width);
                ui.vertical(|ui| {
                    ui.label(
                        egui::RichText::new(self.primary_map_mode.label())
                        .strong(),
                    )
                    .on_hover_text(
                        "Primary splicing map reuses the same transcript/exon lane geometry as the Splicing Expert window",
                    );
                    if let Some(view) = self.current_splicing_expert_view_for_primary_map() {
                        ui.label(
                            egui::RichText::new(format!(
                                "{}:{}..{}  transcripts={}  exons={}  strand {}",
                                view.seq_id,
                                view.region_start_1based,
                                view.region_end_1based,
                                view.transcript_count,
                                view.unique_exon_count,
                                view.strand
                            ))
                            .monospace()
                            .size(self.feature_details_font_size()),
                        );
                        if let Some(feature_id) = self.render_splicing_lane_canvas_ui(
                            ui,
                            &view,
                            SplicingLaneCanvasStyle::primary_map(),
                            true,
                            "primary_map",
                        ) {
                            self.focus_feature(feature_id);
                        }
                    } else {
                        ui.add_space(6.0);
                        ui.label(
                            "Select an mRNA, ncRNA, misc_RNA, transcript, exon, gene, or CDS feature in the feature tree to render splicing lanes.",
                        );
                    }
                });
            });
    }

    pub(super) fn render_splicing_expert_view_ui(
        &mut self,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
        id_namespace: &str,
    ) {
        ui.push_id(
            (
                "splicing_expert_view",
                id_namespace,
                view.seq_id.as_str(),
                view.target_feature_id,
            ),
            |ui| {
        ui.label(
            egui::RichText::new(format!(
                "{}  {}:{}..{}  strand {}  transcripts={}  exons={}",
                view.group_label,
                view.seq_id,
                view.region_start_1based,
                view.region_end_1based,
                view.strand,
                view.transcript_count,
                view.unique_exon_count
            ))
            .monospace()
            .size(self.feature_details_font_size()),
        );
        ui.add_space(2.0);
        let _ =
            self.render_splicing_lane_canvas_ui(
                ui,
                view,
                SplicingLaneCanvasStyle::expert(),
                false,
                "splicing_expert",
            );
        ui.horizontal_wrapped(|ui| {
            ui.label(
                egui::RichText::new("CDS flank phase colors (left/right exon edges):")
                    .size(9.0)
                    .color(egui::Color32::from_rgb(71, 85, 105)),
            )
            .on_hover_text(
                "Flank colors indicate CDS phase at exon entry/exit points when cds_ranges_1based is available",
            );
            for phase in 0_u8..=2 {
                ui.label(
                    egui::RichText::new(format!("{phase}"))
                        .monospace()
                        .size(8.5)
                        .background_color(Self::splicing_cds_phase_color(phase))
                        .color(egui::Color32::WHITE),
                )
                .on_hover_text(format!(
                    "CDS phase {phase} at exon flank (0=blue, 1=amber, 2=rose)"
                ));
            }
            ui.label(
                egui::RichText::new("only shown when CDS ranges are available")
                    .size(8.5)
                    .color(egui::Color32::from_rgb(100, 116, 139)),
            );
        });
        if !view.boundaries.is_empty() {
            let motif_rows = Self::splicing_boundary_motif_rows(view);
            let mut motif_class_counts = std::collections::BTreeMap::<String, usize>::new();
            for row in &motif_rows {
                *motif_class_counts.entry(row.motif_class.clone()).or_default() += 1;
            }
            ui.separator();
            ui.label(
                egui::RichText::new("Splice-site motifs")
                    .monospace()
                    .size(self.feature_details_font_size()),
            );
            ui.horizontal_wrapped(|ui| {
                ui.label(
                    egui::RichText::new(
                        "Known donor/acceptor motif classes are highlighted on the lane canvas and annotated here:",
                    )
                    .size(9.0)
                    .color(egui::Color32::from_rgb(71, 85, 105)),
                );
                for (class, count) in motif_class_counts {
                    let color = Self::splicing_boundary_motif_class_color(class.as_str());
                    ui.label(
                        egui::RichText::new(format!(
                            "{} x{}",
                            Self::splicing_boundary_motif_class_label(class.as_str()),
                            count
                        ))
                        .size(8.5)
                        .color(color),
                    );
                }
            });
            egui::Grid::new("splicing_boundary_motif_grid")
                .striped(true)
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("transcript").monospace().strong());
                    ui.label(egui::RichText::new("donor").monospace().strong());
                    ui.label(egui::RichText::new("acceptor").monospace().strong());
                    ui.label(egui::RichText::new("pair").monospace().strong());
                    ui.label(egui::RichText::new("class").monospace().strong());
                    ui.label(egui::RichText::new("note").strong());
                    ui.end_row();
                    for row in &motif_rows {
                        let class_color =
                            Self::splicing_boundary_motif_class_color(row.motif_class.as_str());
                        ui.label(
                            egui::RichText::new(format!(
                                "n-{} {}",
                                row.transcript_feature_id, row.transcript_id
                            ))
                            .monospace()
                            .size(9.0),
                        );
                        ui.label(
                            egui::RichText::new(format!(
                                "{} {}",
                                row.donor_position_1based, row.donor_motif_2bp
                            ))
                            .monospace()
                            .size(9.0)
                            .color(class_color),
                        );
                        ui.label(
                            egui::RichText::new(format!(
                                "{} {}",
                                row.acceptor_position_1based, row.acceptor_motif_2bp
                            ))
                            .monospace()
                            .size(9.0)
                            .color(class_color),
                        );
                        ui.label(
                            egui::RichText::new(row.paired_motif_signature.as_str())
                                .monospace()
                                .size(9.0)
                                .color(class_color),
                        );
                        ui.label(
                            egui::RichText::new(Self::splicing_boundary_motif_class_label(
                                row.motif_class.as_str(),
                            ))
                            .size(9.0)
                            .color(class_color),
                        );
                        ui.label(
                            egui::RichText::new(row.annotation.as_str())
                                .size(9.0)
                                .color(egui::Color32::from_rgb(71, 85, 105)),
                        );
                        ui.end_row();
                    }
                });
        }
        if !view.intron_signals.is_empty() {
            let signal_rows = Self::splicing_intron_signal_rows(view);
            let mut selected_signal_key = self.splicing_selected_intron_signal_key(view);
            ui.separator();
            ui.label(
                egui::RichText::new("Acceptor-proximal intron signals")
                    .monospace()
                    .size(self.feature_details_font_size()),
            );
            ui.label(
                egui::RichText::new(
                    "These are conservative heuristics, not a splice predictor: branchpoint-like adenines are searched in the usual 18-40 nt acceptor-proximal window and pyrimidine-rich tracts are summarized near the acceptor.",
                )
                .size(9.0)
                .color(egui::Color32::from_rgb(71, 85, 105)),
            );
            egui::Grid::new("splicing_intron_signal_grid")
                .striped(true)
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("transcript").monospace().strong());
                    ui.label(egui::RichText::new("intron").monospace().strong());
                    ui.label(egui::RichText::new("branchpoint-like").monospace().strong());
                    ui.label(egui::RichText::new("polyY tract").monospace().strong());
                    ui.label(egui::RichText::new("note").strong());
                    ui.end_row();
                    for row in &signal_rows {
                        let row_key = Self::splicing_intron_signal_key_for_row(row);
                        let is_selected = selected_signal_key
                            .as_ref()
                            .map(|selected| *selected == row_key)
                            .unwrap_or(false);
                        if ui
                            .selectable_label(
                                is_selected,
                                egui::RichText::new(format!(
                                    "n-{} {}",
                                    row.transcript_feature_id, row.transcript_id
                                ))
                                .monospace()
                                .size(9.0),
                            )
                            .on_hover_text(
                                "Click to focus this intron signal and highlight the intron in the lane view.",
                            )
                            .clicked()
                        {
                            selected_signal_key = Some(row_key.clone());
                        }
                        if ui
                            .selectable_label(
                                is_selected,
                                egui::RichText::new(format!(
                                    "{}..{} ({} bp)",
                                    row.donor_position_1based,
                                    row.acceptor_position_1based,
                                    row.intron_length_bp
                                ))
                                .monospace()
                                .size(9.0),
                            )
                            .clicked()
                        {
                            selected_signal_key = Some(row_key.clone());
                        }
                        let branchpoint_text = row
                            .branchpoint_position_1based
                            .map(|position| format!("{position} {}", row.branchpoint_motif))
                            .unwrap_or_else(|| "none".to_string());
                        if ui
                            .selectable_label(
                                is_selected,
                                egui::RichText::new(branchpoint_text)
                                    .monospace()
                                    .size(9.0)
                                    .color(egui::Color32::from_rgb(124, 58, 237)),
                            )
                            .clicked()
                        {
                            selected_signal_key = Some(row_key.clone());
                        }
                        let tract_text = match (
                            row.polypyrimidine_start_1based,
                            row.polypyrimidine_end_1based,
                        ) {
                            (Some(start), Some(end)) => {
                                format!("{start}..{end} ({}%)", row.polypyrimidine_fraction_percent)
                            }
                            _ => "none".to_string(),
                        };
                        if ui
                            .selectable_label(
                                is_selected,
                                egui::RichText::new(tract_text)
                                    .monospace()
                                    .size(9.0)
                                    .color(egui::Color32::from_rgb(37, 99, 235)),
                            )
                            .clicked()
                        {
                            selected_signal_key = Some(row_key.clone());
                        }
                        ui.vertical(|ui| {
                            ui.label(
                                egui::RichText::new(row.branchpoint_annotation.as_str())
                                    .size(9.0)
                                    .color(egui::Color32::from_rgb(71, 85, 105)),
                            );
                            ui.label(
                                egui::RichText::new(row.polypyrimidine_annotation.as_str())
                                    .size(9.0)
                                    .color(egui::Color32::from_rgb(100, 116, 139)),
                            );
                        });
                        ui.end_row();
                    }
                });
            self.splicing_expert_selected_intron_signal_key = selected_signal_key.clone();
            if let Some(selected_row) = self.splicing_selected_intron_signal_row(view) {
                ui.add_space(4.0);
                egui::Frame::group(ui.style())
                    .fill(egui::Color32::from_rgba_premultiplied(248, 250, 252, 220))
                    .show(ui, |ui| {
                        ui.label(
                            egui::RichText::new(format!(
                                "Selected intron signal: n-{} {}  {}..{}",
                                selected_row.transcript_feature_id,
                                selected_row.transcript_id,
                                selected_row.donor_position_1based,
                                selected_row.acceptor_position_1based
                            ))
                            .monospace()
                            .strong()
                            .size(9.5),
                        );
                        ui.label(
                            egui::RichText::new(
                                Self::splicing_intron_signal_description_text(&selected_row),
                            )
                            .size(9.0)
                            .color(egui::Color32::from_rgb(51, 65, 85)),
                        );
                        ui.label(
                            egui::RichText::new(selected_row.branchpoint_annotation.as_str())
                                .size(9.0)
                                .color(egui::Color32::from_rgb(109, 40, 217)),
                        );
                        ui.label(
                            egui::RichText::new(selected_row.polypyrimidine_annotation.as_str())
                                .size(9.0)
                                .color(egui::Color32::from_rgb(29, 78, 216)),
                        );
                    });
            }
        }

        if !view.matrix_rows.is_empty() && !view.unique_exons.is_empty() {
            let row_count = view.matrix_rows.len();
            let exon_count = view.unique_exons.len();
            let matrix_cell_count = row_count.saturating_mul(exon_count);
            let matrix_heavy = Self::splicing_matrix_should_default_collapsed(row_count, exon_count);
            ui.separator();
            if matrix_heavy {
                ui.small(
                    egui::RichText::new(format!(
                        "Large transcript x exon matrix ({}x{} = {} cells); collapsed by default to reduce idle CPU usage.",
                        row_count, exon_count, matrix_cell_count
                    ))
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
            }
            egui::CollapsingHeader::new(
                egui::RichText::new(format!(
                    "Transcript x exon matrix ({} x {})",
                    row_count, exon_count
                ))
                .monospace()
                .size(self.feature_details_font_size()),
            )
            .id_salt((
                "splicing_expert_matrix_section",
                id_namespace,
                view.seq_id.as_str(),
                view.target_feature_id,
            ))
            .default_open(!matrix_heavy)
            .show(ui, |ui| {
                let transcript_total = view.transcript_count.max(1);
                egui::ScrollArea::horizontal().show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    egui::Grid::new("splicing_expert_matrix")
                        .striped(true)
                        .show(ui, |ui| {
                            ui.label(egui::RichText::new("transcript").monospace().strong());
                            for exon in &view.unique_exons {
                                let mod3 = Self::splicing_exon_length_mod3(exon);
                                let (bg, fg) = Self::splicing_exon_mod3_colors(mod3);
                                ui.label(
                                    egui::RichText::new(format!(
                                        "{}..{}",
                                        exon.start_1based, exon.end_1based
                                    ))
                                    .monospace()
                                    .size(9.0)
                                    .background_color(bg)
                                    .color(fg),
                                )
                                .on_hover_text(format!(
                                    "len={} bp, len%3={mod3} (heuristic frame cue)",
                                    Self::splicing_exon_length_bp(exon)
                                ));
                            }
                            ui.end_row();
                            ui.label(egui::RichText::new("len%3").monospace().size(8.5).strong())
                                .on_hover_text(
                                    "Genomic exon length modulo 3 (heuristic reading-frame cue)",
                                );
                            for exon in &view.unique_exons {
                                let mod3 = Self::splicing_exon_length_mod3(exon);
                                let (bg, fg) = Self::splicing_exon_mod3_colors(mod3);
                                ui.label(
                                    egui::RichText::new(format!("{mod3}"))
                                        .monospace()
                                        .size(8.5)
                                        .background_color(bg)
                                        .color(fg),
                                )
                                .on_hover_text(format!(
                                    "len={} bp, len%3={mod3} (heuristic frame cue)",
                                    Self::splicing_exon_length_bp(exon)
                                ));
                            }
                            ui.end_row();
                            ui.label(
                                egui::RichText::new("support")
                                    .monospace()
                                    .size(9.0)
                                    .strong(),
                            );
                            for exon in &view.unique_exons {
                                let support_ratio = Self::support_ratio_percent(
                                    exon.support_transcript_count,
                                    transcript_total,
                                ) as f32
                                    / 100.0;
                                let support_label = Self::format_support_fraction(
                                    exon.support_transcript_count,
                                    transcript_total,
                                );
                                let label = if exon.constitutive {
                                    format!("{support_label} const")
                                } else {
                                    support_label
                                };
                                ui.label(egui::RichText::new(label).monospace().size(8.5).color(
                                    if exon.constitutive {
                                        let rgb = Self::mix_rgb(
                                            [59, 130, 246],
                                            [30, 64, 175],
                                            support_ratio,
                                        );
                                        egui::Color32::from_rgb(rgb[0], rgb[1], rgb[2])
                                    } else {
                                        let rgb = Self::mix_rgb(
                                            [107, 114, 128],
                                            [30, 64, 175],
                                            support_ratio,
                                        );
                                        egui::Color32::from_rgb(rgb[0], rgb[1], rgb[2])
                                    },
                                ));
                            }
                            ui.end_row();
                            for row in &view.matrix_rows {
                                ui.label(
                                    egui::RichText::new(format!(
                                        "n-{} {}",
                                        row.transcript_feature_id, row.transcript_id
                                    ))
                                    .monospace()
                                    .size(9.0),
                                );
                                for (col_idx, present) in row.exon_presence.iter().enumerate() {
                                    let support_count = view
                                        .unique_exons
                                        .get(col_idx)
                                        .map(|exon| exon.support_transcript_count)
                                        .unwrap_or(0);
                                    let support_ratio =
                                        Self::support_ratio_percent(support_count, transcript_total)
                                            as f32
                                            / 100.0;
                                    let (bg, fg) =
                                        Self::splicing_matrix_cell_colors(*present, support_ratio);
                                    let marker = if *present { "X" } else { "·" };
                                    ui.label(
                                        egui::RichText::new(marker)
                                            .monospace()
                                            .size(9.0)
                                            .color(fg)
                                            .background_color(bg),
                                    )
                                    .on_hover_text(format!(
                                        "Support {}",
                                        Self::format_support_fraction(
                                            support_count,
                                            transcript_total
                                        )
                                    ));
                                }
                                ui.end_row();
                            }
                        });
                });
            });
        }

        if view.unique_exons.len() >= 2 {
            let exon_count = view.unique_exons.len();
            let transition_cell_count = exon_count.saturating_mul(exon_count);
            let transition_heavy = Self::splicing_transition_should_default_collapsed(exon_count);
            ui.separator();
            if transition_heavy {
                ui.small(
                    egui::RichText::new(format!(
                        "Large exon transition matrix ({}x{} = {} cells); collapsed by default to reduce idle CPU usage.",
                        exon_count, exon_count, transition_cell_count
                    ))
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
            }
            egui::CollapsingHeader::new(
                egui::RichText::new(format!(
                    "Exon -> exon transition matrix ({} x {})",
                    exon_count, exon_count
                ))
                .monospace()
                .size(self.feature_details_font_size()),
            )
            .id_salt((
                "splicing_transition_matrix_section",
                id_namespace,
                view.seq_id.as_str(),
                view.target_feature_id,
            ))
            .default_open(!transition_heavy)
            .show(ui, |ui| {
                let transitions = compute_splicing_exon_transition_matrix(view);
                let transcript_total = view.transcript_count.max(1);
                ui.label(
                    egui::RichText::new(
                        "Cell color encodes transition support frequency; exon header color encodes exon genomic length %3 (heuristic frame cue).",
                    )
                    .size(9.0)
                    .color(egui::Color32::from_rgb(71, 85, 105)),
                );
                egui::ScrollArea::horizontal().show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    egui::Grid::new("splicing_exon_transition_grid")
                        .striped(true)
                        .show(ui, |ui| {
                            ui.label(
                                egui::RichText::new("from \\ to")
                                    .monospace()
                                    .size(9.0)
                                    .strong(),
                            );
                            for (to_idx, exon) in view.unique_exons.iter().enumerate() {
                                let mod3 = Self::splicing_exon_length_mod3(exon);
                                let (bg, fg) = Self::splicing_exon_mod3_colors(mod3);
                                ui.label(
                                    egui::RichText::new(format!("E{}", to_idx + 1))
                                        .monospace()
                                        .size(9.0)
                                        .background_color(bg)
                                        .color(fg),
                                )
                                .on_hover_text(format!(
                                    "{}..{} (len={} bp, len%3={mod3})",
                                    exon.start_1based,
                                    exon.end_1based,
                                    Self::splicing_exon_length_bp(exon),
                                ));
                            }
                            ui.end_row();
                            for (from_idx, from_exon) in view.unique_exons.iter().enumerate() {
                                let mod3 = Self::splicing_exon_length_mod3(from_exon);
                                let (bg, fg) = Self::splicing_exon_mod3_colors(mod3);
                                ui.label(
                                    egui::RichText::new(format!("E{}", from_idx + 1))
                                        .monospace()
                                        .size(9.0)
                                        .background_color(bg)
                                        .color(fg),
                                )
                                .on_hover_text(format!(
                                    "{}..{} (len={} bp, len%3={mod3})",
                                    from_exon.start_1based,
                                    from_exon.end_1based,
                                    Self::splicing_exon_length_bp(from_exon),
                                ));
                                for to_idx in 0..view.unique_exons.len() {
                                    let support_count = transitions
                                        .counts
                                        .get(from_idx)
                                        .and_then(|row| row.get(to_idx))
                                        .copied()
                                        .unwrap_or(0);
                                    let support_ratio =
                                        Self::support_ratio_percent(support_count, transcript_total)
                                            as f32
                                            / 100.0;
                                    let (bg, fg) =
                                        Self::splicing_transition_cell_colors(support_ratio);
                                    let marker = if support_count > 0 {
                                        support_count.to_string()
                                    } else {
                                        "·".to_string()
                                    };
                                    let participants = transitions
                                        .transcript_feature_ids
                                        .get(from_idx)
                                        .and_then(|row| row.get(to_idx))
                                        .cloned()
                                        .unwrap_or_default();
                                    let mut labels = participants
                                        .iter()
                                        .map(|feature_id| format!("n-{feature_id}"))
                                        .collect::<Vec<_>>();
                                    if labels.len() > 12 {
                                        let hidden = labels.len() - 12;
                                        labels.truncate(12);
                                        labels.push(format!("+{hidden} more"));
                                    }
                                    ui.label(
                                        egui::RichText::new(marker)
                                            .monospace()
                                            .size(9.0)
                                            .background_color(bg)
                                            .color(fg),
                                    )
                                    .on_hover_text(format!(
                                        "E{} -> E{} support {}{}",
                                        from_idx + 1,
                                        to_idx + 1,
                                        Self::format_support_fraction(
                                            support_count,
                                            transcript_total
                                        ),
                                        if labels.is_empty() {
                                            "".to_string()
                                        } else {
                                            format!("\nTranscripts: {}", labels.join(", "))
                                        }
                                    ));
                                }
                                ui.end_row();
                            }
                        });
                });
            });
        }

        if !view.junctions.is_empty() {
            ui.separator();
            ui.label(
                egui::RichText::new("Junction transition support")
                    .monospace()
                    .size(self.feature_details_font_size()),
            );
            let transcript_total = view.transcript_count.max(1);
            let mut junctions = view.junctions.clone();
            junctions.sort_by(|left, right| {
                right
                    .support_transcript_count
                    .cmp(&left.support_transcript_count)
                    .then_with(|| left.donor_1based.cmp(&right.donor_1based))
                    .then_with(|| left.acceptor_1based.cmp(&right.acceptor_1based))
            });
            egui::ScrollArea::horizontal().show(ui, |ui| {
                scroll_input_policy::apply_scrollarea_keyboard_navigation(
                    ui,
                    scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                );
                egui::Grid::new("splicing_junction_support_grid")
                    .striped(true)
                    .show(ui, |ui| {
                        ui.label(egui::RichText::new("donor").monospace().strong());
                        ui.label(egui::RichText::new("acceptor").monospace().strong());
                        ui.label(egui::RichText::new("span").monospace().strong());
                        ui.label(egui::RichText::new("support").monospace().strong());
                        ui.label(egui::RichText::new("transcripts").monospace().strong());
                        ui.end_row();
                        for junction in &junctions {
                            ui.label(
                                egui::RichText::new(junction.donor_1based.to_string())
                                    .monospace()
                                    .size(9.0),
                            );
                            ui.label(
                                egui::RichText::new(junction.acceptor_1based.to_string())
                                    .monospace()
                                    .size(9.0),
                            );
                            let span_bp = junction
                                .acceptor_1based
                                .max(junction.donor_1based)
                                .saturating_sub(
                                    junction.acceptor_1based.min(junction.donor_1based),
                                );
                            ui.label(
                                egui::RichText::new(span_bp.to_string())
                                    .monospace()
                                    .size(9.0),
                            );
                            ui.label(
                                egui::RichText::new(Self::format_support_fraction(
                                    junction.support_transcript_count,
                                    transcript_total,
                                ))
                                .monospace()
                                .size(9.0),
                            );
                            let mut transcript_labels = junction
                                .transcript_feature_ids
                                .iter()
                                .map(|feature_id| format!("n-{feature_id}"))
                                .collect::<Vec<_>>();
                            if transcript_labels.len() > 10 {
                                let hidden = transcript_labels.len() - 10;
                                transcript_labels.truncate(10);
                                transcript_labels.push(format!("+{hidden} more"));
                            }
                            ui.label(
                                egui::RichText::new(transcript_labels.join(", "))
                                    .monospace()
                                    .size(9.0),
                            );
                            ui.end_row();
                        }
                    });
            });
        }

        if !view.events.is_empty() {
            ui.separator();
            ui.label(
                egui::RichText::new("Splicing event summary")
                    .monospace()
                    .size(self.feature_details_font_size()),
            );
            for event in &view.events {
                ui.label(
                    egui::RichText::new(format!("{}: {}", event.event_type, event.count))
                        .monospace()
                        .size(self.feature_details_font_size()),
                );
                for detail in event.details.iter().take(3) {
                    ui.label(
                        egui::RichText::new(format!("  - {detail}"))
                            .monospace()
                            .size(9.0),
                    );
                }
            }
        }
            },
        );
    }

    pub(super) fn render_isoform_architecture_expert_view_ui(
        &self,
        ui: &mut egui::Ui,
        view: &IsoformArchitectureExpertView,
    ) {
        let font_size = self.feature_details_font_size();
        ui.label(
            egui::RichText::new(format!(
                "{} panel '{}' on {} ({}..{})",
                view.gene_symbol,
                view.panel_id,
                view.seq_id,
                view.region_start_1based,
                view.region_end_1based
            ))
            .monospace()
            .size(font_size + 1.0),
        );
        if let Some(source) = view
            .panel_source
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
        {
            ui.label(egui::RichText::new(format!("source: {source}")).size(font_size));
        }
        if view
            .panel_source
            .as_deref()
            .map(str::trim)
            .map(|value| {
                value.starts_with("UniProt projection ")
                    || value.starts_with("Transcript-native protein")
            })
            .unwrap_or(false)
        {
            ui.label(
                egui::RichText::new(
                    "Protein Expert: transcript-native translation is primary; optional external protein opinions are treated as comparison overlays. CONFLICT stays hidden by default; the top panel keeps coordinate-true exon/intron context; the lower panel reuses shared genomic exon-family columns before handing off to isoform-local protein lanes; exon colors follow genomic exon family/position; membrane topology features render in a dedicated lower band.",
                )
                .size(font_size)
                .color(egui::Color32::from_rgb(71, 85, 105)),
            );
        }
        if !view.warnings.is_empty() {
            ui.label(
                egui::RichText::new(format!("warnings: {}", view.warnings.len()))
                    .size(font_size)
                    .color(egui::Color32::from_rgb(180, 83, 9)),
            );
            for warning in view.warnings.iter().take(6) {
                ui.label(
                    egui::RichText::new(format!(" - {warning}"))
                        .size(font_size)
                        .color(egui::Color32::from_rgb(180, 83, 9)),
                );
            }
        }

        ui.separator();
        ui.label(
            egui::RichText::new("Transcript/exon architecture")
                .strong()
                .size(font_size),
        );
        egui::Grid::new("isoform_transcript_grid")
            .num_columns(5)
            .striped(true)
            .show(ui, |ui| {
                ui.label(egui::RichText::new("isoform").strong().size(font_size));
                ui.label(egui::RichText::new("transcript").strong().size(font_size));
                ui.label(egui::RichText::new("strand").strong().size(font_size));
                ui.label(egui::RichText::new("mapped").strong().size(font_size));
                ui.label(egui::RichText::new("exons").strong().size(font_size));
                ui.end_row();
                for lane in &view.transcript_lanes {
                    ui.label(egui::RichText::new(&lane.label).size(font_size));
                    ui.label(
                        egui::RichText::new(lane.transcript_id.as_deref().unwrap_or("-"))
                            .monospace()
                            .size(font_size),
                    );
                    ui.label(egui::RichText::new(&lane.strand).size(font_size));
                    ui.label(
                        egui::RichText::new(if lane.mapped { "yes" } else { "no" })
                            .size(font_size)
                            .color(if lane.mapped {
                                egui::Color32::from_rgb(22, 101, 52)
                            } else {
                                egui::Color32::from_rgb(180, 83, 9)
                            }),
                    );
                    ui.label(
                        egui::RichText::new(
                            lane.exons
                                .iter()
                                .map(|exon| format!("{}..{}", exon.start_1based, exon.end_1based))
                                .collect::<Vec<_>>()
                                .join(", "),
                        )
                        .monospace()
                        .size(font_size),
                    );
                    ui.end_row();
                }
            });

        ui.separator();
        ui.label(
            egui::RichText::new("Protein/domain architecture")
                .strong()
                .size(font_size),
        );
        egui::Grid::new("isoform_protein_grid")
            .num_columns(10)
            .striped(true)
            .show(ui, |ui| {
                ui.label(egui::RichText::new("isoform").strong().size(font_size));
                ui.label(egui::RichText::new("derived aa").strong().size(font_size));
                ui.label(
                    egui::RichText::new("reference span")
                        .strong()
                        .size(font_size),
                );
                ui.label(
                    egui::RichText::new("translation table/source")
                        .strong()
                        .size(font_size),
                );
                ui.label(
                    egui::RichText::new("speed profile/source")
                        .strong()
                        .size(font_size),
                );
                ui.label(egui::RichText::new("mode").strong().size(font_size));
                ui.label(egui::RichText::new("external opinion").strong().size(font_size));
                ui.label(egui::RichText::new("status").strong().size(font_size));
                ui.label(egui::RichText::new("mismatch summary").strong().size(font_size));
                ui.label(egui::RichText::new("domains").strong().size(font_size));
                ui.end_row();
                for lane in &view.protein_lanes {
                    let comparison = lane.comparison.as_ref();
                    let derived = comparison.and_then(|comparison| comparison.derived.as_ref());
                    let external = comparison.and_then(|comparison| comparison.external_opinion.as_ref());
                    let status = comparison.map(|comparison| comparison.status.as_str()).unwrap_or("-");
                    let status_color = match comparison.map(|comparison| comparison.status) {
                        Some(gentle_protocol::TranscriptProteinComparisonStatus::DerivedOnly)
                        | Some(
                            gentle_protocol::TranscriptProteinComparisonStatus::ConsistentWithExternalOpinion,
                        ) => egui::Color32::from_rgb(22, 101, 52),
                        Some(
                            gentle_protocol::TranscriptProteinComparisonStatus::LowConfidenceExternalOpinion,
                        ) => egui::Color32::from_rgb(180, 83, 9),
                        Some(
                            gentle_protocol::TranscriptProteinComparisonStatus::ExternalOpinionOnly,
                        )
                        | Some(gentle_protocol::TranscriptProteinComparisonStatus::NoTranscriptCds) => {
                            egui::Color32::from_rgb(148, 163, 184)
                        }
                        None => egui::Color32::from_rgb(71, 85, 105),
                    };
                    ui.label(egui::RichText::new(&lane.label).size(font_size));
                    ui.label(
                        egui::RichText::new(
                            derived
                                .map(|derivation| derivation.protein_length_aa.to_string())
                                .or_else(|| {
                                    lane.expected_length_aa
                                        .filter(|value| *value > 0)
                                        .map(|value| value.to_string())
                                })
                                .unwrap_or_else(|| "-".to_string()),
                        )
                        .monospace()
                        .size(font_size),
                    );
                    ui.label(
                        egui::RichText::new(
                            match (lane.reference_start_aa, lane.reference_end_aa) {
                                (Some(start), Some(end)) => format!("{start}..{end}"),
                                (Some(start), None) => format!("{start}.."),
                                (None, Some(end)) => format!("..{end}"),
                                (None, None) => "-".to_string(),
                            },
                        )
                        .monospace()
                        .size(font_size),
                    );
                    ui.label(
                        egui::RichText::new(
                            derived
                                .map(|derivation| {
                                    format!(
                                        "{} ({})",
                                        derivation.translation_table_label,
                                        derivation.translation_table_source.as_str()
                                    )
                                })
                                .unwrap_or_else(|| "-".to_string()),
                        )
                        .monospace()
                        .size(font_size),
                    );
                    ui.label(
                        egui::RichText::new(
                            derived
                                .map(|derivation| {
                                    let profile = derivation
                                        .translation_speed_profile_hint
                                        .as_deref()
                                        .unwrap_or("-");
                                    let source = derivation
                                        .translation_speed_profile_source
                                        .map(|source| source.as_str())
                                        .unwrap_or("unspecified");
                                    let reference = derivation
                                        .translation_speed_reference_species
                                        .as_deref()
                                        .unwrap_or("-");
                                    format!("{profile} ({source}, ref={reference})")
                                })
                                .unwrap_or_else(|| "-".to_string()),
                        )
                        .monospace()
                        .size(font_size),
                    );
                    ui.label(
                        egui::RichText::new(
                            derived
                                .map(|derivation| derivation.derivation_mode.as_str().to_string())
                                .unwrap_or_else(|| "-".to_string()),
                        )
                        .monospace()
                        .size(font_size),
                    );
                    ui.label(
                        egui::RichText::new(
                            external
                                .map(|external| external.source_label.clone())
                                .unwrap_or_else(|| "none".to_string()),
                        )
                        .size(font_size),
                    );
                    ui.label(
                        egui::RichText::new(
                            status.replace('_', " "),
                        )
                        .size(font_size)
                        .color(status_color),
                    );
                    ui.label(
                        egui::RichText::new(
                            comparison
                                .map(|comparison| {
                                    if comparison.mismatch_reasons.is_empty() {
                                        "-".to_string()
                                    } else {
                                        comparison.mismatch_reasons.join(" | ")
                                    }
                                })
                                .unwrap_or_else(|| "-".to_string()),
                        )
                        .monospace()
                        .size(font_size),
                    );
                    ui.label(
                        egui::RichText::new(
                            lane.domains
                                .iter()
                                .map(|domain| {
                                    format!(
                                        "{}:{}..{}",
                                        domain.name, domain.start_aa, domain.end_aa
                                    )
                                })
                                .collect::<Vec<_>>()
                                .join(", "),
                        )
                        .monospace()
                        .size(font_size),
                    );
                    ui.end_row();
                }
            });
    }

    pub(super) fn render_feature_expert_view_ui(
        &mut self,
        ui: &mut egui::Ui,
        view: &FeatureExpertView,
    ) {
        match view {
            FeatureExpertView::Tfbs(tfbs) => self.render_tfbs_expert_view_ui(ui, tfbs),
            FeatureExpertView::RestrictionSite(re) => {
                self.render_restriction_expert_view_ui(ui, re)
            }
            FeatureExpertView::Splicing(splicing) => {
                self.render_splicing_expert_view_ui(ui, splicing, "feature_expert_panel")
            }
            FeatureExpertView::IsoformArchitecture(isoform) => {
                self.render_isoform_architecture_expert_view_ui(ui, isoform)
            }
        }
    }

    pub(super) fn sync_splicing_expert_window_state(&mut self) {
        if let Some(FeatureExpertView::Splicing(view)) = self.description_cache_expert_view.clone()
        {
            if !self.show_splicing_expert_window {
                return;
            }
            let needs_refresh = self.splicing_expert_window_feature_id
                != Some(view.target_feature_id)
                || self
                    .splicing_expert_window_view
                    .as_ref()
                    .map(|cached| {
                        cached.seq_id != view.seq_id
                            || cached.region_start_1based != view.region_start_1based
                            || cached.region_end_1based != view.region_end_1based
                            || cached.transcript_count != view.transcript_count
                    })
                    .unwrap_or(true);
            if needs_refresh {
                self.splicing_expert_window_feature_id = Some(view.target_feature_id);
                self.splicing_expert_window_view = Some(Arc::new(view.clone()));
                self.splicing_expert_window_pending_initial_render = true;
                self.log_splicing_expert_status(&view, "window state refreshed", true);
            }
        }
    }

    pub(super) fn splicing_expert_viewport_id(seq_id: &str, feature_id: usize) -> egui::ViewportId {
        egui::ViewportId::from_hash_of(("splicing_expert_viewport", seq_id, feature_id))
    }

    pub(super) fn splicing_expert_window_default_size() -> Vec2 {
        Vec2::new(1040.0, 760.0)
    }

    pub(super) fn splicing_expert_window_min_size() -> Vec2 {
        Vec2::new(760.0, 520.0)
    }

    pub(super) fn splicing_expert_window_content_min_size() -> Vec2 {
        Vec2::new(1120.0, 700.0)
    }

    pub(super) fn splicing_expert_window_title(view: &SplicingExpertView) -> String {
        format!("Splicing Expert - {} ({})", view.group_label, view.seq_id)
    }

    pub(super) fn rna_read_mapping_viewport_id(
        seq_id: &str,
        feature_id: usize,
    ) -> egui::ViewportId {
        egui::ViewportId::from_hash_of(("rna_read_mapping_viewport", seq_id, feature_id))
    }

    pub(super) fn rna_read_mapping_window_default_size() -> Vec2 {
        Vec2::new(1120.0, 820.0)
    }

    pub(super) fn rna_read_mapping_window_min_size() -> Vec2 {
        Vec2::new(860.0, 620.0)
    }

    pub(super) fn rna_read_mapping_window_content_min_size() -> Vec2 {
        Vec2::new(1080.0, 760.0)
    }

    pub(super) fn rna_read_mapping_window_title(view: &SplicingExpertView) -> String {
        format!("RNA-read Mapping - {} ({})", view.group_label, view.seq_id)
    }

    pub(super) fn stale_auxiliary_window_title_layer_id(title: &str) -> egui::LayerId {
        egui::LayerId::new(egui::Order::Middle, egui::Id::new(title.to_string()))
    }

    pub(super) fn reset_auxiliary_window_areas_if_legacy_title_layer_visible(
        ctx: &egui::Context,
        title: &str,
    ) -> bool {
        let stale_title_layer = Self::stale_auxiliary_window_title_layer_id(title);
        if ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer)) {
            ctx.memory_mut(|mem| mem.reset_areas());
            true
        } else {
            false
        }
    }

    pub(super) fn rna_read_mapping_embedded_window_id(view: &SplicingExpertView) -> egui::Id {
        egui::Id::new(format!(
            "rna_read_mapping_window_embedded_{}_{}",
            view.seq_id, view.target_feature_id
        ))
    }

    pub(super) fn render_rna_read_mapping_embedded_window_shell(
        &mut self,
        ctx: &egui::Context,
        title: &str,
        view: &SplicingExpertView,
        default_size: Vec2,
        content_min_size: Vec2,
        pending_initial_render: bool,
    ) {
        let mut open = self.show_rna_read_mapping_window;
        Self::reset_auxiliary_window_areas_if_legacy_title_layer_visible(ctx, title);
        egui::Window::new(title)
            .id(Self::rna_read_mapping_embedded_window_id(view))
            .open(&mut open)
            .resizable(true)
            .default_size(default_size)
            .show(ctx, |ui| {
                let backdrop_settings = current_window_backdrop_settings();
                paint_window_backdrop(ui, WindowBackdropKind::Splicing, &backdrop_settings);
                egui::ScrollArea::both()
                    .id_salt(format!(
                        "rna_read_mapping_scroll_embedded_{}_{}",
                        view.seq_id, view.target_feature_id
                    ))
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        scroll_input_policy::apply_scrollarea_keyboard_navigation(
                            ui,
                            scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                        );
                        ui.set_min_size(content_min_size);
                        ui.horizontal_wrapped(|ui| {
                            ui.label(
                                egui::RichText::new("Workspace guide [?]")
                                    .size(9.0)
                                    .color(egui::Color32::from_rgb(71, 85, 105)),
                            )
                            .on_hover_text(Self::splicing_nanopore_cdna_panel_help_text());
                            ui.label(
                                egui::RichText::new(
                                    "This dedicated workspace owns RNA-read mapping controls, workflow staging, and report exports for the current splicing locus.",
                                )
                                .size(9.0)
                                .color(egui::Color32::from_rgb(100, 116, 139)),
                            );
                        });
                        ui.horizontal_wrapped(|ui| {
                            if ui
                                .button("Show in Splicing Expert")
                                .on_hover_text(
                                    "Jump back to the Splicing Expert with the current mapping report selected.",
                                )
                                .clicked()
                            {
                                self.show_splicing_expert_for_rna_read_mapping_view();
                            }
                            ui.small(
                                egui::RichText::new(format!(
                                    "Locus: {} | feature n-{}",
                                    view.seq_id, view.target_feature_id
                                ))
                                .color(egui::Color32::from_rgb(100, 116, 139)),
                            );
                        });
                        self.render_rna_read_mapping_window_body(
                            ctx,
                            ui,
                            view,
                            pending_initial_render,
                        );
                    });
            });
        self.show_rna_read_mapping_window = open;
    }

    pub(super) fn render_rna_read_mapping_window_body(
        &mut self,
        ctx: &egui::Context,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
        pending_initial_render: bool,
    ) {
        ui.horizontal_wrapped(|ui| {
            ui.label(
                egui::RichText::new("Workspace guide [?]")
                    .size(9.0)
                    .color(egui::Color32::from_rgb(71, 85, 105)),
            )
            .on_hover_text(Self::splicing_nanopore_cdna_panel_help_text());
            ui.label(
                egui::RichText::new(
                    "This dedicated workspace owns RNA-read mapping controls, workflow staging, and report exports for the current splicing locus.",
                )
                .size(9.0)
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
        });
        ui.horizontal_wrapped(|ui| {
            if ui
                .button("Show in Splicing Expert")
                .on_hover_text(
                    "Jump back to the Splicing Expert with the current mapping report selected.",
                )
                .clicked()
            {
                self.show_splicing_expert_for_rna_read_mapping_view();
            }
            ui.small(
                egui::RichText::new(format!(
                    "Locus: {} | feature n-{}",
                    view.seq_id, view.target_feature_id
                ))
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
        });
        if pending_initial_render {
            self.log_rna_read_mapping_status(view, "first frame deferred", true);
            ui.add_space(6.0);
            ui.label(
                egui::RichText::new(
                    "Preparing the RNA-read Mapping workspace. Detailed controls will appear on the next repaint.",
                )
                .size(10.0)
                .color(egui::Color32::from_rgb(71, 85, 105)),
            );
            ctx.request_repaint();
            self.rna_read_mapping_window_pending_initial_render = false;
            self.log_rna_read_mapping_status(view, "first frame queued repaint", false);
            return;
        }
        self.log_rna_read_mapping_status(view, "rendering workspace controls", false);
        self.render_rna_read_mapping_workspace_controls(ui, view);
        self.log_rna_read_mapping_status(view, "render complete", false);
    }

    pub(super) fn open_rna_read_mapping_workspace_for_view(&mut self, view: &SplicingExpertView) {
        let desired_report_id = self
            .selected_rna_read_evidence_report_id()
            .filter(|value| !value.trim().is_empty())
            .or_else(|| self.latest_rna_read_report_id_for_splicing_view(view));
        self.log_rna_read_mapping_status(
            view,
            "open requested",
            self.rna_read_mapping_window_pending_initial_render,
        );
        self.rna_read_mapping_window_feature_id = Some(view.target_feature_id);
        self.rna_read_mapping_window_view = Some(Arc::new(view.clone()));
        self.rna_read_mapping_window_pending_initial_render = true;
        self.rna_read_mapping_window_focus_requested = true;
        self.show_rna_read_mapping_window = true;
        self.rna_read_mapping_status.clear();
        self.log_rna_read_mapping_status(view, "window state stored", true);
        if self.rna_reads_ui.report_id.trim().is_empty()
            && let Some(report_id) = desired_report_id
        {
            self.rna_reads_ui.report_id = report_id;
        }
    }

    pub(super) fn show_splicing_expert_for_rna_read_mapping_view(&mut self) {
        let Some(view) = self.rna_read_mapping_window_view.clone() else {
            return;
        };
        self.open_splicing_expert_window_for_view(view.as_ref());
        let report_id = self.rna_reads_ui.report_id.trim();
        if !report_id.is_empty() {
            self.rna_read_evidence_ui.selected_report_id = report_id.to_string();
        }
    }

    pub(super) fn open_isoform_expert_window_for_view(
        &mut self,
        panel_id: &str,
        view: &IsoformArchitectureExpertView,
    ) {
        self.isoform_expert_window_panel_id = Some(panel_id.to_string());
        self.isoform_expert_window_view = Some(Arc::new(view.clone()));
        self.show_isoform_expert_window = true;
    }

    pub(super) fn open_uniprot_projection_expert(
        &mut self,
        projection_id: &str,
        protein_feature_filter: gentle_protocol::ProteinFeatureFilter,
    ) -> Result<IsoformArchitectureExpertView, String> {
        let projection_id = projection_id.trim();
        if projection_id.is_empty() {
            return Err("UniProt projection_id must not be empty".to_string());
        }
        match self.inspect_feature_expert_target(&FeatureExpertTarget::UniprotProjection {
            projection_id: projection_id.to_string(),
            protein_feature_filter,
        }) {
            Ok(FeatureExpertView::IsoformArchitecture(view)) => {
                self.open_isoform_expert_window_for_view(&view.panel_id, &view);
                self.op_status =
                    format!("Opened Protein Expert for UniProt projection '{projection_id}'");
                Ok(view)
            }
            Ok(other) => Err(format!(
                "Unexpected expert-view payload for UniProt projection: {other:?}"
            )),
            Err(err) => Err(err),
        }
    }

    pub(super) fn open_transcript_protein_expert(
        &mut self,
        transcript_id_filter: Option<&str>,
    ) -> Result<IsoformArchitectureExpertView, String> {
        match self.inspect_feature_expert_target(&FeatureExpertTarget::ProteinComparison {
            transcript_id_filter: transcript_id_filter
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            protein_feature_filter: Default::default(),
            external_source: None,
            external_entry_id: None,
        }) {
            Ok(FeatureExpertView::IsoformArchitecture(view)) => {
                self.open_isoform_expert_window_for_view(&view.panel_id, &view);
                self.op_status = if let Some(transcript_id_filter) = transcript_id_filter
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    format!("Opened transcript-native Protein Expert for '{transcript_id_filter}'")
                } else {
                    "Opened transcript-native Protein Expert".to_string()
                };
                Ok(view)
            }
            Ok(other) => Err(format!(
                "Unexpected expert-view payload for transcript protein comparison: {other:?}"
            )),
            Err(err) => Err(err),
        }
    }

    pub(super) fn open_ensembl_entry_protein_expert(
        &mut self,
        transcript_id_filter: Option<&str>,
        entry_id: &str,
        protein_feature_filter: gentle_protocol::ProteinFeatureFilter,
    ) -> Result<IsoformArchitectureExpertView, String> {
        let entry_id = entry_id.trim();
        if entry_id.is_empty() {
            return Err("Ensembl protein entry_id must not be empty".to_string());
        }
        match self.inspect_feature_expert_target(&FeatureExpertTarget::ProteinComparison {
            transcript_id_filter: transcript_id_filter
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            protein_feature_filter,
            external_source: Some(gentle_protocol::ProteinExternalOpinionSource::Ensembl),
            external_entry_id: Some(entry_id.to_string()),
        }) {
            Ok(FeatureExpertView::IsoformArchitecture(view)) => {
                self.open_isoform_expert_window_for_view(&view.panel_id, &view);
                self.op_status = if let Some(transcript_id_filter) = transcript_id_filter
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    format!(
                        "Opened Ensembl Protein Expert for '{}' / '{}'",
                        entry_id, transcript_id_filter
                    )
                } else {
                    format!("Opened Ensembl Protein Expert for '{}'", entry_id)
                };
                Ok(view)
            }
            Ok(other) => Err(format!(
                "Unexpected expert-view payload for Ensembl protein comparison: {other:?}"
            )),
            Err(err) => Err(err),
        }
    }

    pub(super) fn isoform_expert_window_title(
        panel_id: &str,
        seq_id: &str,
        view: &IsoformArchitectureExpertView,
    ) -> String {
        if view
            .panel_source
            .as_deref()
            .map(str::trim)
            .map(|value| {
                value.starts_with("UniProt projection ")
                    || value.starts_with("Transcript-native protein")
            })
            .unwrap_or(false)
        {
            format!("Protein Expert - {} ({seq_id})", view.gene_symbol)
        } else {
            format!("Isoform Expert - {panel_id} ({seq_id})")
        }
    }

    pub(super) fn dotplot_viewport_id(seq_id: &str) -> egui::ViewportId {
        egui::ViewportId::from_hash_of(("dotplot_viewport", seq_id))
    }

    pub(super) fn dotplot_window_default_size(has_loaded_payload: bool) -> Vec2 {
        Vec2::new(1240.0, if has_loaded_payload { 820.0 } else { 620.0 })
    }

    pub(super) fn dotplot_window_min_size(has_loaded_payload: bool) -> Vec2 {
        Vec2::new(900.0, if has_loaded_payload { 560.0 } else { 420.0 })
    }

    pub(super) fn dotplot_window_content_min_size(has_loaded_payload: bool) -> Vec2 {
        Vec2::new(1240.0, if has_loaded_payload { 760.0 } else { 0.0 })
    }

    pub(super) fn dotplot_window_title(seq_id: &str) -> String {
        format!("Dotplot - {seq_id}")
    }

    pub(super) fn open_dotplot_window(&mut self) {
        self.ensure_dotplot_cache_current();
        self.show_dotplot_window = true;
    }

    pub fn collect_open_auxiliary_window_entries(&self) -> Vec<(egui::ViewportId, String, String)> {
        let mut entries = Vec::new();
        if self.show_dotplot_window {
            if let Some(viewport_seq_id) = self.dotplot_window_identity_seq_id() {
                let query_label = self.current_dotplot_query_label();
                entries.push((
                    Self::dotplot_viewport_id(&viewport_seq_id),
                    Self::dotplot_window_title(&query_label),
                    format!("Dotplot workspace for '{query_label}'"),
                ));
            }
        }
        if self.show_splicing_expert_window {
            if let Some(view) = self.splicing_expert_window_view.as_ref() {
                entries.push((
                    Self::splicing_expert_viewport_id(&view.seq_id, view.target_feature_id),
                    Self::splicing_expert_window_title(view),
                    format!(
                        "Splicing expert for feature n-{} on '{}'",
                        view.target_feature_id, view.seq_id
                    ),
                ));
            }
        }
        if self.show_rna_read_mapping_window {
            if let Some(view) = self.rna_read_mapping_window_view.as_ref() {
                entries.push((
                    Self::rna_read_mapping_viewport_id(&view.seq_id, view.target_feature_id),
                    Self::rna_read_mapping_window_title(view),
                    format!(
                        "RNA-read mapping workspace for feature n-{} on '{}'",
                        view.target_feature_id, view.seq_id
                    ),
                ));
            }
        }
        if self.show_isoform_expert_window {
            if let Some(view) = self.isoform_expert_window_view.as_ref() {
                let panel_id = self
                    .isoform_expert_window_panel_id
                    .as_deref()
                    .unwrap_or(view.panel_id.as_str());
                entries.push((
                    Self::isoform_expert_viewport_id(&view.seq_id, panel_id),
                    Self::isoform_expert_window_title(panel_id, &view.seq_id, view),
                    if view
                        .panel_source
                        .as_deref()
                        .map(str::trim)
                        .map(|value| {
                            value.starts_with("UniProt projection ")
                                || value.starts_with("Transcript-native protein")
                        })
                        .unwrap_or(false)
                    {
                        format!("Protein Expert '{}' on '{}'", panel_id, view.seq_id)
                    } else {
                        format!(
                            "Isoform architecture panel '{panel_id}' on '{}'",
                            view.seq_id
                        )
                    },
                ));
            }
        }
        entries
    }

    pub(crate) fn embedded_auxiliary_window_layer_id(
        &self,
        viewport_id: egui::ViewportId,
    ) -> Option<egui::LayerId> {
        if self.show_dotplot_window {
            if let Some(viewport_seq_id) = self.dotplot_window_identity_seq_id() {
                if viewport_id == Self::dotplot_viewport_id(&viewport_seq_id) {
                    return Some(egui::LayerId::new(
                        egui::Order::Middle,
                        egui::Id::new(format!("dotplot_window_embedded_{viewport_seq_id}")),
                    ));
                }
            }
        }
        if self.show_splicing_expert_window {
            if let Some(view) = self.splicing_expert_window_view.as_ref() {
                if viewport_id
                    == Self::splicing_expert_viewport_id(&view.seq_id, view.target_feature_id)
                {
                    return Some(egui::LayerId::new(
                        egui::Order::Middle,
                        egui::Id::new(format!(
                            "splicing_expert_window_embedded_{}_{}",
                            view.seq_id, view.target_feature_id
                        )),
                    ));
                }
            }
        }
        if self.show_rna_read_mapping_window {
            if let Some(view) = self.rna_read_mapping_window_view.as_ref() {
                if viewport_id
                    == Self::rna_read_mapping_viewport_id(&view.seq_id, view.target_feature_id)
                {
                    return Some(egui::LayerId::new(
                        egui::Order::Middle,
                        Self::rna_read_mapping_embedded_window_id(view),
                    ));
                }
            }
        }
        if self.show_isoform_expert_window {
            if let Some(view) = self.isoform_expert_window_view.as_ref() {
                let panel_id = self
                    .isoform_expert_window_panel_id
                    .as_deref()
                    .unwrap_or(view.panel_id.as_str());
                if viewport_id == Self::isoform_expert_viewport_id(&view.seq_id, panel_id) {
                    return Some(egui::LayerId::new(
                        egui::Order::Middle,
                        egui::Id::new(format!(
                            "isoform_expert_window_embedded_{}_{}",
                            view.seq_id, panel_id
                        )),
                    ));
                }
            }
        }
        None
    }

    pub(super) fn render_dotplot_window(&mut self, ctx: &egui::Context) {
        if !self.show_dotplot_window {
            return;
        }
        let Some(viewport_seq_id) = self.dotplot_window_identity_seq_id() else {
            self.show_dotplot_window = false;
            return;
        };
        let title = Self::dotplot_window_title(&self.current_dotplot_query_label());
        let viewport_id = Self::dotplot_viewport_id(&viewport_seq_id);
        let has_loaded_payload = self.dotplot_cached_view.is_some()
            || (self.dotplot_ui.show_flexibility_track && self.dotplot_cached_flex_track.is_some());
        let default_size = Self::dotplot_window_default_size(has_loaded_payload);
        let min_size = Self::dotplot_window_min_size(has_loaded_payload);
        let content_min_size = Self::dotplot_window_content_min_size(has_loaded_payload);
        let builder = egui::ViewportBuilder::default()
            .with_title(title.clone())
            .with_inner_size([default_size.x, default_size.y])
            .with_min_inner_size([min_size.x, min_size.y]);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            if class == egui::ViewportClass::EmbeddedWindow {
                Self::reset_auxiliary_window_areas_if_legacy_title_layer_visible(
                    ctx,
                    title.as_str(),
                );
                let mut open = self.show_dotplot_window;
                egui::Window::new(title.clone())
                    .id(egui::Id::new(format!(
                        "dotplot_window_embedded_{}",
                        viewport_seq_id
                    )))
                    .open(&mut open)
                    .resizable(true)
                    .default_size(default_size)
                    .show(ctx, |ui| {
                        let backdrop_settings = current_window_backdrop_settings();
                        paint_window_backdrop(ui, WindowBackdropKind::Sequence, &backdrop_settings);
                        egui::ScrollArea::both()
                            .id_salt(format!(
                                "dotplot_window_scroll_embedded_{}",
                                viewport_seq_id
                            ))
                            .auto_shrink([false, false])
                            .show(ui, |ui| {
                                scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                    ui,
                                    scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                                );
                                ui.set_min_size(content_min_size);
                                self.render_dotplot_workspace_ui(ui);
                            });
                    });
                self.show_dotplot_window = open;
                return;
            }

            crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                let backdrop_settings = current_window_backdrop_settings();
                paint_window_backdrop(ui, WindowBackdropKind::Sequence, &backdrop_settings);
                egui::ScrollArea::both()
                    .id_salt(format!(
                        "dotplot_window_scroll_viewport_{}",
                        viewport_seq_id
                    ))
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        scroll_input_policy::apply_scrollarea_keyboard_navigation(
                            ui,
                            scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                        );
                        ui.set_min_size(content_min_size);
                        self.render_dotplot_workspace_ui(ui);
                    });
            });

            if crate::app::GENtleApp::viewport_close_requested_or_shortcut(ctx) {
                self.show_dotplot_window = false;
            }
        });
    }
}
