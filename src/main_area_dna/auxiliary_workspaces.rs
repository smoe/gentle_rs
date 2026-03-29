//! Dotplot, splicing, and auxiliary workspace support for `MainAreaDna`.
//!
//! This submodule groups the heavyweight secondary workspaces and their shared
//! dotplot support so the main sequence-window file can shrink toward clearer
//! GUI slices ahead of a later crate split.

use super::*;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct DotplotOpsUiState {
    pub(super) half_window_bp: String,
    pub(super) word_size: String,
    pub(super) step_bp: String,
    pub(super) max_mismatches: String,
    pub(super) tile_bp: String,
    pub(super) mode: DotplotMode,
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

impl MainAreaDna {
    pub(super) fn dotplot_window_identity_seq_id(&self) -> Option<String> {
        self.seq_id
            .clone()
            .or_else(|| self.current_dotplot_query_seq_id())
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

    pub(super) fn clear_dotplot_query_override(&mut self) {
        self.dotplot_query_override_seq_id.clear();
        self.dotplot_query_override_source_label.clear();
        self.invalidate_dotplot_cache();
    }

    pub(super) fn dotplot_selection_sync_enabled_for_view(&self, view: &DotplotView) -> bool {
        self.seq_id.as_deref() == Some(view.seq_id.as_str())
    }

    pub(super) fn dotplot_ids_for_active_sequence(&self) -> Vec<String> {
        let Some(seq_id) = self.current_dotplot_query_seq_id() else {
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
        let Some(seq_id) = self.current_dotplot_query_seq_id() else {
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
        if view.points.is_empty()
            || view.reference_span_end_0based <= view.reference_span_start_0based
        {
            return None;
        }
        let mut min_hit = usize::MAX;
        let mut max_hit = 0usize;
        for point in &view.points {
            min_hit = min_hit.min(point.y_0based);
            max_hit = max_hit.max(point.y_0based);
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
    ) -> Option<(f32, f32, String)> {
        if view.point_count == 0 || view.points.is_empty() {
            return None;
        }
        let query_span = view
            .span_end_0based
            .saturating_sub(view.span_start_0based)
            .max(1);
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
        let mut threshold = if view.point_count >= 800 {
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
        if view.point_count < 300 {
            gain = gain.max(1.5);
        }

        let summary = format!(
            "Auto contrast: threshold={:.3}, gain={:.3}, occupancy={:.2}%, q90={:.3}, points={}, cells={}/{}",
            threshold,
            gain,
            occupancy * 100.0,
            q90,
            view.point_count,
            cells.len(),
            total_cells
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
        let Some((threshold, gain, summary)) = Self::recommend_dotplot_display_from_view(view)
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
        let Some(seq_id) = self.current_dotplot_query_seq_id() else {
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
                .filter(|row| row.seq_id == seq_id)
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
        let Some(seq_id) = self.current_dotplot_query_seq_id() else {
            self.op_status = "No active sequence selected for dotplot computation".to_string();
            return;
        };
        let mut diagnostics_snapshot = self.build_dotplot_compute_diagnostics().ok();
        let half_window_bp = match self.resolve_dotplot_half_window_bp("dotplot half_window_bp") {
            Ok(v) => v,
            Err(e) => {
                self.op_status = e;
                return;
            }
        };
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
        let Some((span_start_0based, span_end_0based)) =
            self.default_dotplot_span_for_view(half_window_bp)
        else {
            self.op_status = "Active sequence is empty; dotplot span unavailable".to_string();
            return;
        };
        if span_end_0based.saturating_sub(span_start_0based) < word_size {
            self.op_status = format!(
                "Dotplot span {}..{} is smaller than word_size {}",
                span_start_0based, span_end_0based, word_size
            );
            return;
        }
        let store_as = if self.dotplot_ui.dotplot_id.trim().is_empty() {
            format!("{seq_id}.dotplot")
        } else {
            self.dotplot_ui.dotplot_id.trim().to_string()
        };
        let requires_reference = Self::dotplot_mode_requires_reference(self.dotplot_ui.mode);
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
        self.dotplot_ui.dotplot_id = store_as.clone();
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
                        diag.query_windows, diag.reference_windows, diag.estimated_pair_evaluations
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
                    view.point_count,
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
            let same_reference = reference_seq_label == view.seq_id.as_str();
            if matches!(view.mode, DotplotMode::SelfForward) {
                ui.small(
                    "Self-forward compares a sequence to itself. A main diagonal is expected identity; off-diagonal signal indicates repeated motifs.",
                );
            }
            if view.point_count == 0 {
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
            } else if Self::dotplot_mode_requires_reference(view.mode)
                && view.point_count <= DOTPLOT_SPARSE_POINT_HINT_THRESHOLD
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

        let top_margin = 26.0;
        let left_margin = 56.0;
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
                                            self.dotplot_ui.reference_seq_id = id.clone();
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
                            save_state = true;
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
                        }
                        ui.label("word")
                            .on_hover_text("Seed/k-mer length used for candidate matching");
                        if Self::add_small_uint_text_edit(ui, &mut self.dotplot_ui.word_size, 2)
                            .on_hover_text("Higher values are faster but less sensitive")
                            .changed()
                        {
                            save_state = true;
                        }
                        ui.label("step")
                            .on_hover_text("Sampling step (bp) along query/reference spans");
                        if Self::add_small_uint_text_edit(ui, &mut self.dotplot_ui.step_bp, 4)
                            .on_hover_text("Higher values reduce work but can miss local structures")
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
                            .on_hover_text("Default 0 = exact seed matches only")
                            .changed()
                        {
                            save_state = true;
                        }
                        ui.label("tile_bp").on_hover_text(
                            "Optional rendering/export hint for tile width in bp (analysis metadata)",
                        );
                        if Self::add_small_uint_text_edit(ui, &mut self.dotplot_ui.tile_bp, 5)
                            .on_hover_text("Leave empty for automatic/default behavior")
                            .changed()
                        {
                            save_state = true;
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
                                            self.dotplot_ui.reference_seq_id = id.clone();
                                            save_state = true;
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
                            ui.small("Leave ref_start/ref_end empty for automatic first-pass full-span compute followed by automatic hit-envelope refit.");
                        });
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
                        if let Some((x_bp, y_bp)) = self.dotplot_locked_crosshair_bp {
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
                            ui.label(
                                egui::RichText::new(format!(
                                    "dotplot '{}' query={} [{}..{}] reference={} [{}..{}] | mode={} | word={} step={} mismatches={} | points={}",
                                    view.dotplot_id,
                                    view.seq_id,
                                    view.span_start_0based.saturating_add(1),
                                    view.span_end_0based,
                                    reference_seq_label,
                                    view.reference_span_start_0based.saturating_add(1),
                                    view.reference_span_end_0based,
                                    view.mode.as_str(),
                                    view.word_size,
                                    view.step_bp,
                                    view.max_mismatches,
                                    view.point_count
                                ))
                                .monospace()
                                .size(self.feature_details_font_size()),
                            );
                            let flex_track = if self.dotplot_ui.show_flexibility_track {
                                cached_track.as_ref()
                            } else {
                                None
                            };
                            self.render_dotplot_density_ui(
                                ui,
                                view,
                                flex_track,
                                self.dotplot_ui.display_density_threshold,
                                self.dotplot_ui.display_intensity_gain,
                            );
                            ui.add_space(6.0);
                            self.render_dotplot_boxplot_summary_ui(ui, view);
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
        &self,
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
            .num_columns(4)
            .striped(true)
            .show(ui, |ui| {
                ui.label(egui::RichText::new("isoform").strong().size(font_size));
                ui.label(egui::RichText::new("expected aa").strong().size(font_size));
                ui.label(egui::RichText::new("TA class").strong().size(font_size));
                ui.label(egui::RichText::new("domains").strong().size(font_size));
                ui.end_row();
                for lane in &view.protein_lanes {
                    ui.label(egui::RichText::new(&lane.label).size(font_size));
                    ui.label(
                        egui::RichText::new(
                            lane.expected_length_aa
                                .map(|v| v.to_string())
                                .unwrap_or_else(|| "-".to_string()),
                        )
                        .monospace()
                        .size(font_size),
                    );
                    ui.label(
                        egui::RichText::new(lane.transactivation_class.as_deref().unwrap_or("-"))
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
        &self,
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
        if let Some(FeatureExpertView::Splicing(view)) = self.description_cache_expert_view.as_ref()
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
    ) {
        let mut open = self.show_rna_read_mapping_window;
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
                        self.render_rna_read_mapping_workspace_controls(ui, view);
                    });
            });
        self.show_rna_read_mapping_window = open;
    }

    pub(super) fn open_rna_read_mapping_workspace_for_view(&mut self, view: &SplicingExpertView) {
        let desired_report_id = self
            .selected_rna_read_evidence_report_id()
            .filter(|value| !value.trim().is_empty())
            .or_else(|| self.latest_rna_read_report_id_for_splicing_view(view));
        self.rna_read_mapping_window_feature_id = Some(view.target_feature_id);
        self.rna_read_mapping_window_view = Some(Arc::new(view.clone()));
        self.show_rna_read_mapping_window = true;
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

    pub(super) fn isoform_expert_window_title(panel_id: &str, seq_id: &str) -> String {
        format!("Isoform Expert - {panel_id} ({seq_id})")
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
                    Self::isoform_expert_window_title(panel_id, &view.seq_id),
                    format!(
                        "Isoform architecture panel '{panel_id}' on '{}'",
                        view.seq_id
                    ),
                ));
            }
        }
        entries
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
