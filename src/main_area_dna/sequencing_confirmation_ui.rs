//! Sequencing-confirmation UI support for `MainAreaDna`.
//!
//! This submodule owns the construct-read confirmation form state, report
//! helpers, chromatogram review helpers, and specialist rendering entrypoint.

use super::*;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct SequencingConfirmationUiState {
    pub(super) baseline_seq_id_text: String,
    pub(super) read_seq_ids_text: String,
    pub(super) trace_ids_text: String,
    pub(super) trace_import_path: String,
    pub(super) trace_import_id: String,
    pub(super) trace_import_associate_with_expected_seq: bool,
    pub(super) trace_import_add_to_run: bool,
    pub(super) report_id: String,
    pub(super) selected_report_id: String,
    pub(super) selected_target_id: String,
    pub(super) selected_evidence_id: String,
    pub(super) selected_trace_id: String,
    pub(super) selected_variant_id: String,
    pub(super) selected_review_focus_kind: Option<SequencingConfirmationReviewFocusKind>,
    pub(super) selected_gap_start_0based: Option<usize>,
    pub(super) selected_gap_end_0based_exclusive: Option<usize>,
    pub(super) review_unresolved_first: bool,
    pub(super) primer_seq_ids_text: String,
    pub(super) primer_min_3prime_anneal_bp: String,
    pub(super) primer_predicted_read_length_bp: String,
    #[serde(skip)]
    pub(super) primer_overlay_report: Option<SequencingPrimerOverlayReport>,
    pub(super) chromatogram_focus_mode: SequencingChromatogramFocusMode,
    pub(super) chromatogram_trace_base_index: usize,
    pub(super) chromatogram_window_bp: String,
    pub(super) junction_positions_0based: String,
    pub(super) junction_flank_bp: String,
    pub(super) include_full_span_target: bool,
    pub(super) allow_reverse_complement: bool,
    pub(super) alignment_mode: PairwiseAlignmentMode,
    pub(super) match_score: String,
    pub(super) mismatch_score: String,
    pub(super) gap_open: String,
    pub(super) gap_extend: String,
    pub(super) min_identity_fraction: String,
    pub(super) min_target_coverage_fraction: String,
}

impl Default for SequencingConfirmationUiState {
    fn default() -> Self {
        Self {
            baseline_seq_id_text: String::new(),
            read_seq_ids_text: String::new(),
            trace_ids_text: String::new(),
            trace_import_path: String::new(),
            trace_import_id: String::new(),
            trace_import_associate_with_expected_seq: true,
            trace_import_add_to_run: true,
            report_id: "seq_confirm_gui".to_string(),
            selected_report_id: String::new(),
            selected_target_id: String::new(),
            selected_evidence_id: String::new(),
            selected_trace_id: String::new(),
            selected_variant_id: String::new(),
            selected_review_focus_kind: None,
            selected_gap_start_0based: None,
            selected_gap_end_0based_exclusive: None,
            review_unresolved_first: false,
            primer_seq_ids_text: String::new(),
            primer_min_3prime_anneal_bp: "18".to_string(),
            primer_predicted_read_length_bp: "800".to_string(),
            primer_overlay_report: None,
            chromatogram_focus_mode: SequencingChromatogramFocusMode::VariantLocus,
            chromatogram_trace_base_index: 0,
            chromatogram_window_bp: "24".to_string(),
            junction_positions_0based: String::new(),
            junction_flank_bp: "12".to_string(),
            include_full_span_target: true,
            allow_reverse_complement: true,
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: "2".to_string(),
            mismatch_score: "-3".to_string(),
            gap_open: "-5".to_string(),
            gap_extend: "-1".to_string(),
            min_identity_fraction: "0.80".to_string(),
            min_target_coverage_fraction: "1.0".to_string(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub(super) enum SequencingChromatogramFocusMode {
    #[default]
    VariantLocus,
    TraceBaseBrowser,
}

impl SequencingChromatogramFocusMode {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::VariantLocus => "Variant locus",
            Self::TraceBaseBrowser => "Trace base browser",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(super) enum SequencingConfirmationReviewFocusKind {
    Target,
    Evidence,
    Variant,
    CoverageGap,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum SequencingConfirmationOverviewSelection {
    Target(String),
    Evidence(String),
    Variant(String),
    CoverageGap(usize, usize),
}

impl MainAreaDna {
    pub(super) fn sequencing_confirmation_default_report_id(&self) -> String {
        let seq_component = self
            .seq_id
            .as_deref()
            .map(|value| Self::sanitize_id_component(value, "seq"))
            .unwrap_or_else(|| "seq".to_string());
        format!("{seq_component}_seq_confirm")
    }

    pub(super) fn parse_sequencing_confirmation_junction_positions(
        text: &str,
    ) -> Result<Vec<usize>, String> {
        let trimmed = text.trim();
        if trimmed.is_empty() {
            return Ok(vec![]);
        }
        let mut out = vec![];
        let mut seen = BTreeSet::new();
        for token in trimmed.split(|c| matches!(c, ',' | ';' | '\n' | '\r')) {
            let token = token.trim();
            if token.is_empty() {
                continue;
            }
            let value = token.parse::<usize>().map_err(|_| {
                format!(
                    "Invalid sequencing-confirmation junction position '{}': expected 0-based integer",
                    token
                )
            })?;
            if seen.insert(value) {
                out.push(value);
            }
        }
        Ok(out)
    }

    pub(super) fn build_sequencing_confirmation_targets(
        &self,
    ) -> Result<Vec<SequencingConfirmationTargetSpec>, String> {
        let ui = &self.sequencing_confirmation_ui;
        let seq_len = self
            .dna
            .read()
            .map_err(|_| "DNA lock poisoned while preparing sequencing confirmation".to_string())?
            .len();
        if seq_len == 0 {
            return Err(
                "Sequencing confirmation requires a non-empty expected construct sequence"
                    .to_string(),
            );
        }
        let junction_flank = Self::parse_positive_usize_text(
            &ui.junction_flank_bp,
            "sequencing_confirmation.junction_flank_bp",
        )?;
        let junction_positions =
            Self::parse_sequencing_confirmation_junction_positions(&ui.junction_positions_0based)?;
        let mut targets = vec![];
        if ui.include_full_span_target {
            targets.push(SequencingConfirmationTargetSpec {
                target_id: "full_span".to_string(),
                label: "Full construct span".to_string(),
                kind: SequencingConfirmationTargetKind::FullSpan,
                start_0based: 0,
                end_0based_exclusive: seq_len,
                junction_left_end_0based: None,
                expected_bases: None,
                baseline_bases: None,
                required: true,
            });
        }
        for (idx, left_end) in junction_positions.iter().copied().enumerate() {
            if left_end > seq_len {
                return Err(format!(
                    "Sequencing-confirmation junction position {} is outside expected sequence length {}",
                    left_end, seq_len
                ));
            }
            let start_0based = left_end.saturating_sub(junction_flank);
            let end_0based_exclusive = left_end.saturating_add(junction_flank).min(seq_len);
            if end_0based_exclusive <= start_0based {
                return Err(format!(
                    "Sequencing-confirmation junction {} yields empty support interval {}..{}",
                    left_end, start_0based, end_0based_exclusive
                ));
            }
            targets.push(SequencingConfirmationTargetSpec {
                target_id: format!("junction_{}", idx + 1),
                label: format!("Junction @ {left_end}"),
                kind: SequencingConfirmationTargetKind::Junction,
                start_0based,
                end_0based_exclusive,
                junction_left_end_0based: Some(left_end),
                expected_bases: None,
                baseline_bases: None,
                required: true,
            });
        }
        if targets.is_empty() {
            return Err(
                "Enable full-span confirmation or add at least one junction breakpoint".to_string(),
            );
        }
        Ok(targets)
    }

    pub(super) fn build_confirm_construct_reads_operation(
        &self,
    ) -> Result<(String, Operation), String> {
        let expected_seq_id = self
            .seq_id
            .clone()
            .ok_or_else(|| "No active expected construct sequence".to_string())?;
        let ui = &self.sequencing_confirmation_ui;
        let read_seq_ids = Self::parse_ids(&ui.read_seq_ids_text);
        let trace_ids = Self::parse_ids(&ui.trace_ids_text);
        if read_seq_ids.is_empty() && trace_ids.is_empty() {
            return Err(
                "Provide at least one read sequence ID and/or one imported trace ID before running confirmation"
                    .to_string(),
            );
        }
        let report_id = if ui.report_id.trim().is_empty() {
            self.sequencing_confirmation_default_report_id()
        } else {
            Self::sanitize_id_component(&ui.report_id, "seq_confirm_gui")
        };
        let targets = self.build_sequencing_confirmation_targets()?;
        Ok((
            report_id.clone(),
            Operation::ConfirmConstructReads {
                expected_seq_id,
                baseline_seq_id: {
                    let trimmed = ui.baseline_seq_id_text.trim();
                    if trimmed.is_empty() {
                        None
                    } else {
                        Some(trimmed.to_string())
                    }
                },
                read_seq_ids,
                trace_ids,
                targets,
                alignment_mode: ui.alignment_mode,
                match_score: Self::parse_required_i32_text(
                    &ui.match_score,
                    "sequencing_confirmation.match_score",
                )?,
                mismatch_score: Self::parse_required_i32_text(
                    &ui.mismatch_score,
                    "sequencing_confirmation.mismatch_score",
                )?,
                gap_open: Self::parse_required_i32_text(
                    &ui.gap_open,
                    "sequencing_confirmation.gap_open",
                )?,
                gap_extend: Self::parse_required_i32_text(
                    &ui.gap_extend,
                    "sequencing_confirmation.gap_extend",
                )?,
                min_identity_fraction: Self::parse_fraction_0_to_1(
                    &ui.min_identity_fraction,
                    "sequencing_confirmation.min_identity_fraction",
                )?,
                min_target_coverage_fraction: Self::parse_fraction_0_to_1(
                    &ui.min_target_coverage_fraction,
                    "sequencing_confirmation.min_target_coverage_fraction",
                )?,
                allow_reverse_complement: ui.allow_reverse_complement,
                report_id: Some(report_id),
            },
        ))
    }

    pub(super) fn build_import_sequencing_trace_operation(&self) -> Result<Operation, String> {
        let ui = &self.sequencing_confirmation_ui;
        let path = ui.trace_import_path.trim();
        if path.is_empty() {
            return Err(
                "Choose a local ABI/AB1/SCF trace file before importing sequencing evidence"
                    .to_string(),
            );
        }
        let seq_id = if ui.trace_import_associate_with_expected_seq {
            Some(
                self.seq_id
                    .clone()
                    .ok_or_else(|| "No active expected construct sequence".to_string())?,
            )
        } else {
            None
        };
        let trace_id =
            (!ui.trace_import_id.trim().is_empty()).then(|| ui.trace_import_id.trim().to_string());
        Ok(Operation::ImportSequencingTrace {
            path: path.to_string(),
            trace_id,
            seq_id,
        })
    }

    pub(super) fn build_suggest_sequencing_primers_operation(&self) -> Result<Operation, String> {
        let expected_seq_id = self
            .seq_id
            .clone()
            .ok_or_else(|| "No active expected construct sequence".to_string())?;
        let ui = &self.sequencing_confirmation_ui;
        let primer_seq_ids = Self::parse_ids(&ui.primer_seq_ids_text);
        let confirmation_report_id = {
            let trimmed = ui.selected_report_id.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        };
        if primer_seq_ids.is_empty() && confirmation_report_id.is_none() {
            return Err(
                "Provide at least one primer sequence ID or select a saved confirmation report before suggesting sequencing-primer overlays"
                    .to_string(),
            );
        }
        Ok(Operation::SuggestSequencingPrimers {
            expected_seq_id,
            primer_seq_ids,
            confirmation_report_id,
            min_3prime_anneal_bp: Self::parse_positive_usize_text(
                &ui.primer_min_3prime_anneal_bp,
                "sequencing_confirmation.primer_min_3prime_anneal_bp",
            )?,
            predicted_read_length_bp: Self::parse_positive_usize_text(
                &ui.primer_predicted_read_length_bp,
                "sequencing_confirmation.primer_predicted_read_length_bp",
            )?,
        })
    }

    pub(super) fn sequencing_confirmation_overlay_report_for_expected_seq(
        &self,
        expected_seq_id: &str,
    ) -> Option<&SequencingPrimerOverlayReport> {
        let expected_seq_id = expected_seq_id.trim();
        if expected_seq_id.is_empty() {
            return None;
        }
        self.sequencing_confirmation_ui
            .primer_overlay_report
            .as_ref()
            .filter(|overlay| {
                overlay
                    .expected_seq_id
                    .eq_ignore_ascii_case(expected_seq_id)
            })
    }

    pub(super) fn sequencing_confirmation_report_matched_primer_overlay_report(
        &self,
        expected_seq_id: &str,
        report_id: &str,
    ) -> Option<&SequencingPrimerOverlayReport> {
        let report_id = report_id.trim();
        self.sequencing_confirmation_overlay_report_for_expected_seq(expected_seq_id)
            .filter(|overlay| {
                overlay
                    .confirmation_report_id
                    .as_deref()
                    .is_none_or(|id| id.eq_ignore_ascii_case(report_id))
            })
    }

    pub(super) fn build_sequencing_confirmation_unresolved_summary_markdown(
        &self,
        report_id: &str,
    ) -> Result<(SequencingConfirmationReport, String), EngineError> {
        let Some(engine) = self.engine.clone() else {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No engine attached".to_string(),
            });
        };
        let report = engine
            .read()
            .expect("Engine lock poisoned")
            .get_sequencing_confirmation_report(report_id)?;
        let overlay = self.sequencing_confirmation_report_matched_primer_overlay_report(
            report.expected_seq_id.as_str(),
            report.report_id.as_str(),
        );
        let text =
            GentleEngine::format_sequencing_confirmation_unresolved_summary_markdown_with_overlay(
                &report, overlay,
            );
        Ok((report, text))
    }

    pub(super) fn sequencing_trace_summaries(&self) -> Vec<SequencingTraceSummary> {
        let expected_seq_id = self.seq_id.as_deref().unwrap_or("");
        let Some(engine) = self.engine.as_ref() else {
            return vec![];
        };
        let mut rows = engine
            .read()
            .ok()
            .map(|guard| guard.list_sequencing_traces(None))
            .unwrap_or_default();
        rows.sort_by(|left, right| {
            let left_associated = left.seq_id.as_deref() == Some(expected_seq_id);
            let right_associated = right.seq_id.as_deref() == Some(expected_seq_id);
            right_associated
                .cmp(&left_associated)
                .then(right.imported_at_unix_ms.cmp(&left.imported_at_unix_ms))
                .then(
                    left.trace_id
                        .to_ascii_lowercase()
                        .cmp(&right.trace_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub(super) fn sequencing_trace_record(&self, trace_id: &str) -> Option<SequencingTraceRecord> {
        let trace_id = trace_id.trim();
        if trace_id.is_empty() {
            return None;
        }
        let engine = self.engine.as_ref()?;
        engine
            .read()
            .ok()
            .and_then(|guard| guard.get_sequencing_trace(trace_id).ok())
    }

    pub(super) fn handle_imported_sequencing_trace_result(&mut self, result: &OpResult) {
        let Some(import_report) = result.sequencing_trace_import_report.as_ref() else {
            return;
        };
        self.sequencing_confirmation_ui.selected_trace_id = import_report.trace_id.clone();
        if self.sequencing_confirmation_ui.trace_import_add_to_run {
            Self::append_unique_csv_token(
                &mut self.sequencing_confirmation_ui.trace_ids_text,
                import_report.trace_id.clone(),
            );
        }
        self.save_engine_ops_state();
    }

    pub(super) fn sequencing_trace_called_base_preview(
        text: &str,
        line_bp: usize,
        max_lines: usize,
    ) -> String {
        if text.is_empty() || line_bp == 0 || max_lines == 0 {
            return "<no called bases>".to_string();
        }
        let chars = text.chars().collect::<Vec<_>>();
        let mut out = String::new();
        let max_bp = line_bp.saturating_mul(max_lines);
        for (line_idx, chunk) in chars.chunks(line_bp).take(max_lines).enumerate() {
            if line_idx > 0 {
                out.push('\n');
            }
            for ch in chunk {
                out.push(*ch);
            }
        }
        if chars.len() > max_bp {
            out.push_str("\n...");
        }
        out
    }

    pub(super) fn sequencing_trace_confidence_summary(values: &[u8]) -> String {
        if values.is_empty() {
            return "confidences: none".to_string();
        }
        let mut min_value = u8::MAX;
        let mut max_value = u8::MIN;
        let mut sum = 0u64;
        for value in values {
            min_value = min_value.min(*value);
            max_value = max_value.max(*value);
            sum += *value as u64;
        }
        let mean = sum as f64 / values.len() as f64;
        format!(
            "confidences: n={} min={} mean={:.1} max={}",
            values.len(),
            min_value,
            mean,
            max_value
        )
    }

    pub(super) fn sequencing_trace_peak_summary(values: &[u32]) -> String {
        if values.is_empty() {
            return "peaks: none".to_string();
        }
        let preview = values
            .iter()
            .take(8)
            .map(|value| value.to_string())
            .collect::<Vec<_>>();
        let suffix = if values.len() > preview.len() {
            ", ..."
        } else {
            ""
        };
        format!(
            "peaks: n={} [{}{}]",
            values.len(),
            preview.join(", "),
            suffix
        )
    }

    pub(super) fn sequencing_confirmation_variant_color(
        classification: SequencingConfirmationVariantClassification,
    ) -> egui::Color32 {
        match classification {
            SequencingConfirmationVariantClassification::ExpectedMatch => {
                egui::Color32::from_rgb(46, 125, 50)
            }
            SequencingConfirmationVariantClassification::IntendedEditConfirmed => {
                egui::Color32::from_rgb(25, 118, 210)
            }
            SequencingConfirmationVariantClassification::ReferenceReversion => {
                egui::Color32::from_rgb(183, 28, 28)
            }
            SequencingConfirmationVariantClassification::UnexpectedDifference => {
                egui::Color32::from_rgb(230, 81, 0)
            }
            SequencingConfirmationVariantClassification::LowConfidenceOrAmbiguous => {
                egui::Color32::from_rgb(180, 83, 9)
            }
            SequencingConfirmationVariantClassification::InsufficientEvidence => {
                egui::Color32::from_gray(120)
            }
        }
    }

    pub(super) fn sequencing_trace_sample_window(
        trace: &SequencingTraceRecord,
        variant: Option<&SequencingConfirmationVariantRow>,
        trace_base_index: Option<usize>,
        flank_bp: usize,
    ) -> Option<(u32, u32)> {
        let max_point = trace
            .channel_data
            .iter()
            .filter_map(|row| row.points.len().checked_sub(1))
            .max()
            .map(|value| value as u32)?;
        if max_point == 0 {
            return Some((0, 1));
        }
        if let Some(base_index) = trace_base_index {
            return Self::sequencing_trace_sample_window_for_base_index(
                trace, base_index, flank_bp, max_point,
            );
        }
        let Some(variant) = variant else {
            return Some((0, max_point.min(400)));
        };
        let Some(center) = variant.peak_center else {
            return Some((0, max_point.min(400)));
        };
        if trace.peak_locations.is_empty() {
            let start = center.saturating_sub((flank_bp as u32).saturating_mul(8));
            let end = (center + (flank_bp as u32).saturating_mul(8)).min(max_point);
            return Some((start, end.max(start + 1)));
        }
        let nearest_peak_idx = trace
            .peak_locations
            .iter()
            .enumerate()
            .min_by_key(|(_, value)| value.abs_diff(center))
            .map(|(idx, _)| idx)
            .unwrap_or(0);
        let start_idx = nearest_peak_idx.saturating_sub(flank_bp);
        let end_idx =
            (nearest_peak_idx + flank_bp).min(trace.peak_locations.len().saturating_sub(1));
        let start = trace
            .peak_locations
            .get(start_idx)
            .copied()
            .unwrap_or(0)
            .saturating_sub(18);
        let end = trace
            .peak_locations
            .get(end_idx)
            .copied()
            .unwrap_or(max_point)
            .saturating_add(18)
            .min(max_point);
        Some((start, end.max(start + 1)))
    }

    pub(super) fn sequencing_trace_browser_base_count(trace: &SequencingTraceRecord) -> usize {
        let called_count = trace.called_bases.len();
        let peak_count = trace.peak_locations.len();
        match (called_count, peak_count) {
            (0, 0) => 0,
            (0, peak_count) => peak_count,
            (called_count, 0) => called_count,
            (called_count, peak_count) => called_count.min(peak_count),
        }
    }

    pub(super) fn sequencing_trace_clamp_base_index(
        trace: &SequencingTraceRecord,
        requested: usize,
    ) -> usize {
        let count = Self::sequencing_trace_browser_base_count(trace);
        if count == 0 {
            0
        } else {
            requested.min(count.saturating_sub(1))
        }
    }

    pub(super) fn sequencing_trace_sample_window_for_base_index(
        trace: &SequencingTraceRecord,
        base_index: usize,
        flank_bp: usize,
        max_point: u32,
    ) -> Option<(u32, u32)> {
        if !trace.peak_locations.is_empty() {
            let peak_index = Self::sequencing_trace_clamp_base_index(trace, base_index);
            let start_idx = peak_index.saturating_sub(flank_bp);
            let end_idx = (peak_index + flank_bp).min(trace.peak_locations.len().saturating_sub(1));
            let start = trace
                .peak_locations
                .get(start_idx)
                .copied()
                .unwrap_or(0)
                .saturating_sub(18);
            let end = trace
                .peak_locations
                .get(end_idx)
                .copied()
                .unwrap_or(max_point)
                .saturating_add(18)
                .min(max_point);
            return Some((start, end.max(start + 1)));
        }
        let base_count = trace.called_bases.len();
        if base_count == 0 {
            return Some((0, max_point.min(400)));
        }
        let clamped = base_index.min(base_count.saturating_sub(1));
        let center = if base_count <= 1 {
            max_point / 2
        } else {
            (((clamped as f64) / ((base_count - 1) as f64)) * (max_point as f64)).round() as u32
        };
        let span = (flank_bp as u32).saturating_mul(8);
        let start = center.saturating_sub(span);
        let end = (center + span).min(max_point);
        Some((start, end.max(start + 1)))
    }

    pub(super) fn sequencing_trace_base_index_for_variant(
        trace: &SequencingTraceRecord,
        variant: &SequencingConfirmationVariantRow,
    ) -> Option<usize> {
        let center = variant.peak_center?;
        if !trace.peak_locations.is_empty() {
            return trace
                .peak_locations
                .iter()
                .enumerate()
                .min_by_key(|(_, value)| value.abs_diff(center))
                .map(|(idx, _)| idx);
        }
        let base_count = trace.called_bases.len();
        let max_point = trace
            .channel_data
            .iter()
            .filter_map(|row| row.points.len().checked_sub(1))
            .max()
            .map(|value| value as u32)?;
        if base_count == 0 || max_point == 0 {
            return None;
        }
        let ratio = (center as f64 / max_point as f64).clamp(0.0, 1.0);
        Some(((ratio * (base_count.saturating_sub(1) as f64)).round() as usize).min(base_count - 1))
    }

    pub(super) fn sequencing_trace_selected_base_summary(
        trace: &SequencingTraceRecord,
        base_index: usize,
    ) -> Option<String> {
        let count = Self::sequencing_trace_browser_base_count(trace);
        if count == 0 {
            return None;
        }
        let base_index = Self::sequencing_trace_clamp_base_index(trace, base_index);
        let base = trace
            .called_bases
            .as_bytes()
            .get(base_index)
            .copied()
            .map(char::from)
            .unwrap_or('?');
        let confidence = trace
            .called_base_confidence_values
            .get(base_index)
            .map(|value| value.to_string())
            .unwrap_or_else(|| "-".to_string());
        let peak = trace
            .peak_locations
            .get(base_index)
            .map(|value| value.to_string())
            .unwrap_or_else(|| "-".to_string());
        let clip = match (
            trace.clip_start_base_index,
            trace.clip_end_base_index_exclusive,
        ) {
            (Some(start), Some(end)) => {
                if base_index >= start && base_index < end {
                    "inside clip"
                } else {
                    "outside clip"
                }
            }
            (Some(start), None) => {
                if base_index >= start {
                    "inside clip"
                } else {
                    "outside clip"
                }
            }
            _ => "clip unknown",
        };
        Some(format!(
            "base {} / {} = '{}' | qv={} | peak={} | {}",
            base_index + 1,
            count,
            base,
            confidence,
            peak,
            clip
        ))
    }

    pub(super) fn sequencing_trace_curve_unavailable_message(
        trace: &SequencingTraceRecord,
    ) -> Option<&'static str> {
        if trace.channel_data.is_empty() {
            Some("curve data unavailable; re-import this trace to inspect chromatogram curves")
        } else {
            None
        }
    }

    pub(super) fn render_sequencing_trace_chromatogram(
        ui: &mut egui::Ui,
        trace: &SequencingTraceRecord,
        variant: Option<&SequencingConfirmationVariantRow>,
        trace_base_index: Option<usize>,
        flank_bp: usize,
    ) {
        if let Some(message) = Self::sequencing_trace_curve_unavailable_message(trace) {
            ui.small(message);
            return;
        }
        let Some((sample_start, sample_end)) =
            Self::sequencing_trace_sample_window(trace, variant, trace_base_index, flank_bp)
        else {
            ui.small("Trace did not include enough raw curve points to render a chromatogram.");
            return;
        };
        let desired_size = egui::vec2(ui.available_width().max(320.0), 240.0);
        let (rect, _response) = ui.allocate_exact_size(desired_size, egui::Sense::hover());
        let painter = ui.painter_at(rect);
        painter.rect_stroke(
            rect,
            4.0,
            egui::Stroke::new(1.0, egui::Color32::from_gray(110)),
            egui::StrokeKind::Outside,
        );
        let mut max_y = 0u32;
        for channel in &trace.channel_data {
            for (idx, value) in channel.points.iter().enumerate() {
                let x = idx as u32;
                if x >= sample_start && x <= sample_end {
                    max_y = max_y.max(*value);
                }
            }
        }
        if max_y == 0 {
            max_y = 1;
        }
        let palette = [
            ("A", egui::Color32::from_rgb(46, 125, 50)),
            ("C", egui::Color32::from_rgb(25, 118, 210)),
            ("G", egui::Color32::from_rgb(230, 81, 0)),
            ("T", egui::Color32::from_rgb(183, 28, 28)),
        ];
        for channel in &trace.channel_data {
            let color = palette
                .iter()
                .find(|(label, _)| channel.channel.eq_ignore_ascii_case(label))
                .map(|(_, color)| *color)
                .unwrap_or(egui::Color32::LIGHT_GRAY);
            let mut points = Vec::new();
            for (idx, value) in channel.points.iter().enumerate() {
                let x = idx as u32;
                if x < sample_start || x > sample_end {
                    continue;
                }
                let x_t = (x.saturating_sub(sample_start)) as f32
                    / (sample_end.saturating_sub(sample_start).max(1)) as f32;
                let y_t = *value as f32 / max_y as f32;
                points.push(egui::pos2(
                    egui::lerp(rect.left()..=rect.right(), x_t),
                    egui::lerp((rect.bottom() - 16.0)..=rect.top() + 8.0, y_t),
                ));
            }
            if points.len() >= 2 {
                painter.add(egui::Shape::line(points, egui::Stroke::new(1.5, color)));
            }
        }
        let peak_count = trace.called_bases.len().min(trace.peak_locations.len());
        for idx in 0..peak_count {
            let peak = trace.peak_locations[idx];
            if peak < sample_start || peak > sample_end {
                continue;
            }
            let x_t = (peak.saturating_sub(sample_start)) as f32
                / (sample_end.saturating_sub(sample_start).max(1)) as f32;
            let x = egui::lerp(rect.left()..=rect.right(), x_t);
            painter.line_segment(
                [
                    egui::pos2(x, rect.bottom() - 14.0),
                    egui::pos2(x, rect.bottom() - 6.0),
                ],
                egui::Stroke::new(1.0, egui::Color32::from_gray(120)),
            );
            if peak_count <= 64 {
                let base = trace.called_bases.as_bytes()[idx] as char;
                let color = palette
                    .iter()
                    .find(|(label, _)| label.starts_with(base))
                    .map(|(_, color)| *color)
                    .unwrap_or(egui::Color32::WHITE);
                painter.text(
                    egui::pos2(x, rect.bottom() - 4.0),
                    egui::Align2::CENTER_BOTTOM,
                    base,
                    egui::TextStyle::Small.resolve(ui.style()),
                    color,
                );
            }
        }
        if let Some(variant) = variant {
            if let Some(center) = variant.peak_center
                && center >= sample_start
                && center <= sample_end
            {
                let x_t = (center.saturating_sub(sample_start)) as f32
                    / (sample_end.saturating_sub(sample_start).max(1)) as f32;
                let x = egui::lerp(rect.left()..=rect.right(), x_t);
                painter.line_segment(
                    [
                        egui::pos2(x, rect.top() + 4.0),
                        egui::pos2(x, rect.bottom() - 18.0),
                    ],
                    egui::Stroke::new(
                        2.0,
                        Self::sequencing_confirmation_variant_color(variant.classification),
                    ),
                );
            }
            painter.text(
                rect.left_top() + egui::vec2(8.0, 6.0),
                egui::Align2::LEFT_TOP,
                format!(
                    "{} [{}] expected='{}' observed='{}'{}",
                    variant.label,
                    variant.classification.as_str(),
                    variant.expected_bases,
                    if variant.observed_bases.is_empty() {
                        "-"
                    } else {
                        variant.observed_bases.as_str()
                    },
                    variant
                        .baseline_bases
                        .as_deref()
                        .map(|baseline| format!(" baseline='{}'", baseline))
                        .unwrap_or_default()
                ),
                egui::TextStyle::Small.resolve(ui.style()),
                Self::sequencing_confirmation_variant_color(variant.classification),
            );
        } else if let Some(base_index) = trace_base_index {
            let base_index = Self::sequencing_trace_clamp_base_index(trace, base_index);
            if let Some(center) = trace.peak_locations.get(base_index).copied()
                && center >= sample_start
                && center <= sample_end
            {
                let x_t = (center.saturating_sub(sample_start)) as f32
                    / (sample_end.saturating_sub(sample_start).max(1)) as f32;
                let x = egui::lerp(rect.left()..=rect.right(), x_t);
                painter.line_segment(
                    [
                        egui::pos2(x, rect.top() + 4.0),
                        egui::pos2(x, rect.bottom() - 18.0),
                    ],
                    egui::Stroke::new(2.0, egui::Color32::from_rgb(25, 118, 210)),
                );
            }
            if let Some(summary) = Self::sequencing_trace_selected_base_summary(trace, base_index) {
                painter.text(
                    rect.left_top() + egui::vec2(8.0, 6.0),
                    egui::Align2::LEFT_TOP,
                    summary,
                    egui::TextStyle::Small.resolve(ui.style()),
                    egui::Color32::from_rgb(25, 118, 210),
                );
            }
        }
    }

    pub(super) fn sequencing_confirmation_report_summaries(
        &self,
    ) -> Vec<SequencingConfirmationReportSummary> {
        let Some(engine) = self.engine.as_ref() else {
            return vec![];
        };
        let expected_seq_id = self.seq_id.as_deref();
        engine
            .read()
            .ok()
            .map(|guard| guard.list_sequencing_confirmation_reports(expected_seq_id))
            .unwrap_or_default()
    }

    pub(super) fn list_sequencing_confirmation_reports(&mut self) {
        let summaries = self.sequencing_confirmation_report_summaries();
        if summaries.is_empty() {
            self.op_status =
                "No persisted sequencing-confirmation reports for this sequence".to_string();
            return;
        }
        if self
            .sequencing_confirmation_ui
            .selected_report_id
            .trim()
            .is_empty()
            || !summaries.iter().any(|row| {
                row.report_id
                    .eq_ignore_ascii_case(self.sequencing_confirmation_ui.selected_report_id.trim())
            })
        {
            if let Some(row) = summaries.last() {
                self.sequencing_confirmation_ui.selected_report_id = row.report_id.clone();
                self.save_engine_ops_state();
            }
        }
        let preview_ids = summaries
            .iter()
            .take(8)
            .map(|row| row.report_id.clone())
            .collect::<Vec<_>>();
        let suffix = if summaries.len() > preview_ids.len() {
            ", ..."
        } else {
            ""
        };
        self.op_status = format!(
            "Sequencing-confirmation reports: {} total [{}{}]",
            summaries.len(),
            preview_ids.join(", "),
            suffix
        );
    }

    pub(super) fn show_sequencing_confirmation_report(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Sequencing-confirmation report_id is empty".to_string();
            return;
        }
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let report = match engine
            .read()
            .expect("Engine lock poisoned")
            .get_sequencing_confirmation_report(report_id)
        {
            Ok(report) => report,
            Err(err) => {
                self.op_status = format!(
                    "Could not load sequencing-confirmation report '{report_id}': {}",
                    err.message
                );
                return;
            }
        };
        self.sequencing_confirmation_ui.selected_report_id = report.report_id.clone();
        self.sequencing_confirmation_ui.selected_evidence_id =
            Self::sequencing_confirmation_selected_evidence_id(
                &report,
                &self.sequencing_confirmation_ui.selected_evidence_id,
            )
            .unwrap_or_default();
        let confirmed = report
            .targets
            .iter()
            .filter(|row| row.status == SequencingConfirmationStatus::Confirmed)
            .count();
        let contradicted = report
            .targets
            .iter()
            .filter(|row| row.status == SequencingConfirmationStatus::Contradicted)
            .count();
        let insufficient = report
            .targets
            .iter()
            .filter(|row| row.status == SequencingConfirmationStatus::InsufficientEvidence)
            .count();
        self.op_status = format!(
            "Sequencing-confirmation report '{}' expected='{}' baseline='{}' status={} reads={} traces={} evidence={} targets={} variants={} (confirmed={}, contradicted={}, insufficient={})",
            report.report_id,
            report.expected_seq_id,
            report.baseline_seq_id.as_deref().unwrap_or("-"),
            report.overall_status.as_str(),
            report.read_seq_ids.len(),
            report.trace_ids.len(),
            report.reads.len(),
            report.targets.len(),
            report.variants.len(),
            confirmed,
            contradicted,
            insufficient
        );
        self.save_engine_ops_state();
    }

    pub(super) fn export_sequencing_confirmation_report_dialog(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Sequencing-confirmation report_id is empty".to_string();
            return;
        }
        let default_name = format!("{report_id}.seq_confirm_report.json");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("JSON", &["json"])
            .save_file();
        let Some(path) = path else {
            self.op_status = "Sequencing-confirmation report export canceled".to_string();
            return;
        };
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let path_text = path.to_string_lossy().to_string();
        let exported = engine
            .read()
            .expect("Engine lock poisoned")
            .export_sequencing_confirmation_report(report_id, &path_text);
        match exported {
            Ok(report) => {
                self.op_status = format!(
                    "Exported sequencing-confirmation report '{}' to {}",
                    report.report_id, path_text
                );
            }
            Err(err) => {
                self.op_status = format!(
                    "Could not export sequencing-confirmation report '{report_id}': {}",
                    err.message
                );
            }
        }
    }

    pub(super) fn export_sequencing_confirmation_support_tsv_dialog(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Sequencing-confirmation report_id is empty".to_string();
            return;
        }
        let default_name = format!("{report_id}.seq_confirm_support.tsv");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("TSV", &["tsv", "txt"])
            .save_file();
        let Some(path) = path else {
            self.op_status = "Sequencing-confirmation TSV export canceled".to_string();
            return;
        };
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let path_text = path.to_string_lossy().to_string();
        let exported = engine
            .read()
            .expect("Engine lock poisoned")
            .export_sequencing_confirmation_support_tsv(report_id, &path_text);
        match exported {
            Ok(report) => {
                self.op_status = format!(
                    "Exported sequencing-confirmation support TSV '{}' to {}",
                    report.report_id, path_text
                );
            }
            Err(err) => {
                self.op_status = format!(
                    "Could not export sequencing-confirmation support TSV '{report_id}': {}",
                    err.message
                );
            }
        }
    }

    pub(super) fn copy_sequencing_confirmation_unresolved_summary(
        &mut self,
        report_id: &str,
        ctx: &egui::Context,
    ) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Sequencing-confirmation report_id is empty".to_string();
            return;
        }
        match self.build_sequencing_confirmation_unresolved_summary_markdown(report_id) {
            Ok((report, text)) => {
                ctx.copy_text(text);
                self.op_status = format!(
                    "Copied unresolved sequencing-confirmation summary for '{}'",
                    report.report_id
                );
            }
            Err(err) => {
                self.op_status = format!(
                    "Could not build unresolved sequencing-confirmation summary '{report_id}': {}",
                    err.message
                );
            }
        }
    }

    pub(super) fn export_sequencing_confirmation_unresolved_summary_dialog(
        &mut self,
        report_id: &str,
    ) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Sequencing-confirmation report_id is empty".to_string();
            return;
        }
        let default_name = format!("{report_id}.seq_confirm_unresolved.md");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("Markdown", &["md", "txt"])
            .save_file();
        let Some(path) = path else {
            self.op_status =
                "Sequencing-confirmation unresolved-summary export canceled".to_string();
            return;
        };
        let path_text = path.to_string_lossy().to_string();
        match self.build_sequencing_confirmation_unresolved_summary_markdown(report_id) {
            Ok((report, text)) => match fs::write(&path, text) {
                Ok(()) => {
                    self.op_status = format!(
                        "Exported unresolved sequencing-confirmation summary '{}' to {}",
                        report.report_id, path_text
                    );
                }
                Err(err) => {
                    self.op_status = format!(
                        "Could not write unresolved sequencing-confirmation summary '{}' to {}: {}",
                        report.report_id, path_text, err
                    );
                }
            },
            Err(err) => {
                self.op_status = format!(
                    "Could not export unresolved sequencing-confirmation summary '{report_id}': {}",
                    err.message
                );
            }
        }
    }

    pub(super) fn sequencing_confirmation_status_color(
        status: SequencingConfirmationStatus,
    ) -> egui::Color32 {
        match status {
            SequencingConfirmationStatus::Confirmed => egui::Color32::from_rgb(46, 125, 50),
            SequencingConfirmationStatus::Contradicted => egui::Color32::from_rgb(183, 28, 28),
            SequencingConfirmationStatus::InsufficientEvidence => {
                egui::Color32::from_rgb(180, 83, 9)
            }
        }
    }

    pub(super) fn sequencing_confirmation_selected_read<'a>(
        report: &'a SequencingConfirmationReport,
        selected_evidence_id: &str,
    ) -> Option<&'a SequencingConfirmationReadResult> {
        let normalized_id = selected_evidence_id.trim();
        if normalized_id.is_empty() {
            return report.reads.first();
        }
        report
            .reads
            .iter()
            .find(|row| row.evidence_id.eq_ignore_ascii_case(normalized_id))
            .or_else(|| report.reads.first())
    }

    pub(super) fn sequencing_confirmation_selected_evidence_id(
        report: &SequencingConfirmationReport,
        selected_evidence_id: &str,
    ) -> Option<String> {
        Self::sequencing_confirmation_selected_read(report, selected_evidence_id)
            .map(|row| row.evidence_id.clone())
    }

    pub(super) fn sequencing_confirmation_target_priority(
        status: SequencingConfirmationStatus,
    ) -> usize {
        match status {
            SequencingConfirmationStatus::Contradicted => 0,
            SequencingConfirmationStatus::InsufficientEvidence => 1,
            SequencingConfirmationStatus::Confirmed => 2,
        }
    }

    pub(super) fn sequencing_confirmation_variant_priority(
        classification: SequencingConfirmationVariantClassification,
    ) -> usize {
        match classification {
            SequencingConfirmationVariantClassification::ReferenceReversion
            | SequencingConfirmationVariantClassification::UnexpectedDifference => 0,
            SequencingConfirmationVariantClassification::LowConfidenceOrAmbiguous
            | SequencingConfirmationVariantClassification::InsufficientEvidence => 1,
            SequencingConfirmationVariantClassification::IntendedEditConfirmed
            | SequencingConfirmationVariantClassification::ExpectedMatch => 2,
        }
    }

    pub(super) fn sequencing_confirmation_target_review_rows<'a>(
        report: &'a SequencingConfirmationReport,
        review_unresolved_first: bool,
    ) -> Vec<&'a SequencingConfirmationTargetResult> {
        let mut rows = report.targets.iter().collect::<Vec<_>>();
        if review_unresolved_first {
            rows.sort_by(|a, b| {
                Self::sequencing_confirmation_target_priority(a.status)
                    .cmp(&Self::sequencing_confirmation_target_priority(b.status))
                    .then(a.start_0based.cmp(&b.start_0based))
                    .then(a.end_0based_exclusive.cmp(&b.end_0based_exclusive))
                    .then(a.label.cmp(&b.label))
            });
        }
        rows
    }

    pub(super) fn sequencing_confirmation_variant_review_rows<'a>(
        report: &'a SequencingConfirmationReport,
        review_unresolved_first: bool,
    ) -> Vec<&'a SequencingConfirmationVariantRow> {
        let mut rows = report.variants.iter().collect::<Vec<_>>();
        if review_unresolved_first {
            rows.sort_by(|a, b| {
                Self::sequencing_confirmation_variant_priority(a.classification)
                    .cmp(&Self::sequencing_confirmation_variant_priority(
                        b.classification,
                    ))
                    .then(a.start_0based.cmp(&b.start_0based))
                    .then(a.end_0based_exclusive.cmp(&b.end_0based_exclusive))
                    .then(a.label.cmp(&b.label))
            });
        }
        rows
    }

    pub(super) fn sequencing_confirmation_selected_target<'a>(
        report: &'a SequencingConfirmationReport,
        selected_target_id: &str,
        review_unresolved_first: bool,
    ) -> Option<&'a SequencingConfirmationTargetResult> {
        let normalized_id = selected_target_id.trim();
        if !normalized_id.is_empty()
            && let Some(row) = report
                .targets
                .iter()
                .find(|row| row.target_id.eq_ignore_ascii_case(normalized_id))
        {
            return Some(row);
        }
        Self::sequencing_confirmation_target_review_rows(report, review_unresolved_first)
            .into_iter()
            .next()
    }

    pub(super) fn sequencing_confirmation_selected_target_id(
        report: &SequencingConfirmationReport,
        selected_target_id: &str,
        review_unresolved_first: bool,
    ) -> Option<String> {
        Self::sequencing_confirmation_selected_target(
            report,
            selected_target_id,
            review_unresolved_first,
        )
        .map(|row| row.target_id.clone())
    }

    pub(super) fn sequencing_confirmation_report_sequence_length(
        report: &SequencingConfirmationReport,
        fallback_len: usize,
    ) -> usize {
        fallback_len.max(
            report
                .targets
                .iter()
                .map(|row| row.end_0based_exclusive)
                .chain(report.variants.iter().map(|row| row.end_0based_exclusive))
                .chain(
                    report
                        .reads
                        .iter()
                        .map(|row| row.best_alignment.aligned_target_end_0based_exclusive),
                )
                .max()
                .unwrap_or(0),
        )
    }

    pub(super) fn sequencing_confirmation_merged_evidence_spans(
        report: &SequencingConfirmationReport,
        sequence_length: usize,
    ) -> Vec<(usize, usize)> {
        if sequence_length == 0 {
            return Vec::new();
        }
        let mut spans = report
            .reads
            .iter()
            .filter_map(|row| {
                let start = row
                    .best_alignment
                    .aligned_target_start_0based
                    .min(sequence_length);
                let end = row
                    .best_alignment
                    .aligned_target_end_0based_exclusive
                    .min(sequence_length);
                (end > start).then_some((start, end))
            })
            .collect::<Vec<_>>();
        spans.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut merged: Vec<(usize, usize)> = Vec::new();
        for (start, end) in spans {
            if let Some(last) = merged.last_mut()
                && start <= last.1
            {
                last.1 = last.1.max(end);
            } else {
                merged.push((start, end));
            }
        }
        merged
    }

    pub(super) fn sequencing_confirmation_uncovered_spans(
        report: &SequencingConfirmationReport,
        sequence_length: usize,
    ) -> Vec<(usize, usize)> {
        if sequence_length == 0 {
            return Vec::new();
        }
        let merged = Self::sequencing_confirmation_merged_evidence_spans(report, sequence_length);
        if merged.is_empty() {
            return vec![(0, sequence_length)];
        }
        let mut uncovered = Vec::new();
        let mut cursor = 0usize;
        for (start, end) in merged {
            if start > cursor {
                uncovered.push((cursor, start));
            }
            cursor = cursor.max(end);
        }
        if cursor < sequence_length {
            uncovered.push((cursor, sequence_length));
        }
        uncovered
    }

    pub(super) fn sequencing_confirmation_target_intersects_span(
        target: &SequencingConfirmationTargetResult,
        start_0based: usize,
        end_0based_exclusive: usize,
    ) -> bool {
        target.start_0based < end_0based_exclusive && start_0based < target.end_0based_exclusive
    }

    pub(super) fn sequencing_confirmation_select_target_for_gap(
        report: &SequencingConfirmationReport,
        gap_start: usize,
        gap_end: usize,
    ) -> Option<String> {
        report
            .targets
            .iter()
            .filter(|row| {
                Self::sequencing_confirmation_target_intersects_span(row, gap_start, gap_end)
            })
            .min_by(|a, b| {
                Self::sequencing_confirmation_target_priority(a.status)
                    .cmp(&Self::sequencing_confirmation_target_priority(b.status))
                    .then(a.start_0based.cmp(&b.start_0based))
            })
            .map(|row| row.target_id.clone())
            .or_else(|| {
                report
                    .targets
                    .iter()
                    .min_by_key(|row| row.start_0based.abs_diff(gap_start))
                    .map(|row| row.target_id.clone())
            })
    }

    pub(super) fn sequencing_confirmation_selected_coverage_gap(
        report: &SequencingConfirmationReport,
        sequence_length: usize,
        selected_gap_start_0based: Option<usize>,
        selected_gap_end_0based_exclusive: Option<usize>,
    ) -> Option<(usize, usize)> {
        let selected = selected_gap_start_0based.zip(selected_gap_end_0based_exclusive)?;
        Self::sequencing_confirmation_uncovered_spans(report, sequence_length)
            .into_iter()
            .find(|row| *row == selected)
    }

    pub(super) fn sequencing_confirmation_span_overlap_bp(
        start_a: usize,
        end_a: usize,
        start_b: usize,
        end_b: usize,
    ) -> usize {
        let overlap_start = start_a.max(start_b);
        let overlap_end = end_a.min(end_b);
        overlap_end.saturating_sub(overlap_start)
    }

    pub(super) fn sequencing_confirmation_gap_center_0based(
        gap_start_0based: usize,
        gap_end_0based_exclusive: usize,
    ) -> usize {
        gap_start_0based + gap_end_0based_exclusive.saturating_sub(gap_start_0based) / 2
    }

    pub(super) fn sequencing_confirmation_gap_flanking_reads<'a>(
        report: &'a SequencingConfirmationReport,
        gap_start_0based: usize,
        gap_end_0based_exclusive: usize,
    ) -> (
        Option<&'a SequencingConfirmationReadResult>,
        Option<&'a SequencingConfirmationReadResult>,
    ) {
        let left = report
            .reads
            .iter()
            .filter(|row| {
                row.best_alignment.aligned_target_end_0based_exclusive <= gap_start_0based
            })
            .max_by(|a, b| {
                a.best_alignment
                    .aligned_target_end_0based_exclusive
                    .cmp(&b.best_alignment.aligned_target_end_0based_exclusive)
                    .then_with(|| {
                        a.best_alignment
                            .identity_fraction
                            .partial_cmp(&b.best_alignment.identity_fraction)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .then(a.evidence_id.cmp(&b.evidence_id))
            });
        let right = report
            .reads
            .iter()
            .filter(|row| {
                row.best_alignment.aligned_target_start_0based >= gap_end_0based_exclusive
            })
            .min_by(|a, b| {
                a.best_alignment
                    .aligned_target_start_0based
                    .cmp(&b.best_alignment.aligned_target_start_0based)
                    .then_with(|| {
                        b.best_alignment
                            .identity_fraction
                            .partial_cmp(&a.best_alignment.identity_fraction)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .then(a.evidence_id.cmp(&b.evidence_id))
            });
        (left, right)
    }

    pub(super) fn sequencing_confirmation_gap_primer_suggestions<'a>(
        overlay_report: &'a SequencingPrimerOverlayReport,
        gap_start_0based: usize,
        gap_end_0based_exclusive: usize,
    ) -> Vec<&'a SequencingPrimerOverlaySuggestion> {
        let gap_center = Self::sequencing_confirmation_gap_center_0based(
            gap_start_0based,
            gap_end_0based_exclusive,
        );
        let mut rows = overlay_report.suggestions.iter().collect::<Vec<_>>();
        rows.sort_by(|a, b| {
            let overlap_a = Self::sequencing_confirmation_span_overlap_bp(
                a.predicted_read_span_start_0based,
                a.predicted_read_span_end_0based_exclusive,
                gap_start_0based,
                gap_end_0based_exclusive,
            );
            let overlap_b = Self::sequencing_confirmation_span_overlap_bp(
                b.predicted_read_span_start_0based,
                b.predicted_read_span_end_0based_exclusive,
                gap_start_0based,
                gap_end_0based_exclusive,
            );
            overlap_b
                .cmp(&overlap_a)
                .then(
                    a.three_prime_position_0based
                        .abs_diff(gap_center)
                        .cmp(&b.three_prime_position_0based.abs_diff(gap_center)),
                )
                .then(a.primer_seq_id.cmp(&b.primer_seq_id))
                .then(a.orientation.as_str().cmp(b.orientation.as_str()))
        });
        rows
    }

    pub(super) fn sequencing_confirmation_gap_primer_proposals<'a>(
        overlay_report: &'a SequencingPrimerOverlayReport,
        gap_start_0based: usize,
        gap_end_0based_exclusive: usize,
    ) -> Vec<&'a SequencingPrimerProposalRow> {
        let gap_center = Self::sequencing_confirmation_gap_center_0based(
            gap_start_0based,
            gap_end_0based_exclusive,
        );
        let mut rows = overlay_report.proposals.iter().collect::<Vec<_>>();
        rows.sort_by(|a, b| {
            let overlap_a = Self::sequencing_confirmation_span_overlap_bp(
                a.predicted_read_span_start_0based,
                a.predicted_read_span_end_0based_exclusive,
                gap_start_0based,
                gap_end_0based_exclusive,
            );
            let overlap_b = Self::sequencing_confirmation_span_overlap_bp(
                b.predicted_read_span_start_0based,
                b.predicted_read_span_end_0based_exclusive,
                gap_start_0based,
                gap_end_0based_exclusive,
            );
            overlap_b
                .cmp(&overlap_a)
                .then(
                    a.three_prime_position_0based
                        .abs_diff(gap_center)
                        .cmp(&b.three_prime_position_0based.abs_diff(gap_center)),
                )
                .then(a.three_prime_distance_bp.cmp(&b.three_prime_distance_bp))
                .then(a.problem_id.cmp(&b.problem_id))
        });
        rows
    }

    pub(super) fn sequencing_confirmation_overview_x(
        rect: egui::Rect,
        sequence_length: usize,
        bp: usize,
    ) -> f32 {
        let denominator = sequence_length.max(1) as f32;
        let fraction = bp.min(sequence_length) as f32 / denominator;
        egui::lerp(rect.left()..=rect.right(), fraction)
    }

    pub(super) fn sequencing_confirmation_overview_span_rect(
        rect: egui::Rect,
        sequence_length: usize,
        start_0based: usize,
        end_0based_exclusive: usize,
    ) -> egui::Rect {
        let left = Self::sequencing_confirmation_overview_x(rect, sequence_length, start_0based);
        let right =
            Self::sequencing_confirmation_overview_x(rect, sequence_length, end_0based_exclusive);
        egui::Rect::from_min_max(
            egui::pos2(left.min(right), rect.top()),
            egui::pos2((left.max(right)).max(left.min(right) + 2.0), rect.bottom()),
        )
    }

    pub(super) fn sequencing_confirmation_evidence_color(
        read: &SequencingConfirmationReadResult,
        selected: bool,
    ) -> egui::Color32 {
        let base = if !read.usable {
            egui::Color32::from_rgb(180, 83, 9)
        } else if read.trace_id.is_some() {
            egui::Color32::from_rgb(25, 118, 210)
        } else {
            egui::Color32::from_rgb(46, 125, 50)
        };
        if selected {
            base
        } else {
            base.gamma_multiply(0.7)
        }
    }

    pub(super) fn sequencing_confirmation_sync_evidence_selection(
        &mut self,
        report: &SequencingConfirmationReport,
        evidence_id: &str,
    ) {
        let Some(selected) = Self::sequencing_confirmation_selected_read(report, evidence_id)
        else {
            self.sequencing_confirmation_ui.selected_evidence_id.clear();
            return;
        };
        self.sequencing_confirmation_ui.selected_review_focus_kind =
            Some(SequencingConfirmationReviewFocusKind::Evidence);
        self.sequencing_confirmation_ui.selected_gap_start_0based = None;
        self.sequencing_confirmation_ui
            .selected_gap_end_0based_exclusive = None;
        self.sequencing_confirmation_ui.selected_evidence_id = selected.evidence_id.clone();
        if let Some(trace_id) = selected.trace_id.as_deref() {
            self.sequencing_confirmation_ui.selected_trace_id = trace_id.to_string();
        }
        if let Some(target_id) = selected
            .contradicted_target_ids
            .first()
            .or(selected.confirmed_target_ids.first())
            .or(selected.covered_target_ids.first())
        {
            self.sequencing_confirmation_ui.selected_target_id = target_id.clone();
        }
        if let Some(variant) = report.variants.iter().find(|row| {
            row.evidence_id
                .eq_ignore_ascii_case(selected.evidence_id.as_str())
        }) {
            self.sequencing_confirmation_ui.selected_variant_id = variant.variant_id.clone();
            if let Some(trace_id) = variant.trace_id.as_deref() {
                self.sequencing_confirmation_ui.selected_trace_id = trace_id.to_string();
            }
        }
    }

    pub(super) fn sequencing_confirmation_sync_variant_selection(
        &mut self,
        report: &SequencingConfirmationReport,
        variant_id: &str,
    ) {
        let Some(variant) = report
            .variants
            .iter()
            .find(|row| row.variant_id.eq_ignore_ascii_case(variant_id.trim()))
        else {
            self.sequencing_confirmation_ui.selected_variant_id.clear();
            return;
        };
        self.sequencing_confirmation_ui.selected_review_focus_kind =
            Some(SequencingConfirmationReviewFocusKind::Variant);
        self.sequencing_confirmation_ui.selected_gap_start_0based = None;
        self.sequencing_confirmation_ui
            .selected_gap_end_0based_exclusive = None;
        self.sequencing_confirmation_ui.selected_variant_id = variant.variant_id.clone();
        if let Some(target_id) = variant.target_id.as_deref() {
            self.sequencing_confirmation_ui.selected_target_id = target_id.to_string();
        }
        if !variant.evidence_id.trim().is_empty() {
            self.sequencing_confirmation_ui.selected_evidence_id = variant.evidence_id.clone();
        }
        if let Some(trace_id) = variant.trace_id.as_deref() {
            self.sequencing_confirmation_ui.selected_trace_id = trace_id.to_string();
        } else if !variant.evidence_id.trim().is_empty() {
            self.sequencing_confirmation_sync_evidence_selection(report, &variant.evidence_id);
            self.sequencing_confirmation_ui.selected_variant_id = variant.variant_id.clone();
            self.sequencing_confirmation_ui.selected_review_focus_kind =
                Some(SequencingConfirmationReviewFocusKind::Variant);
        }
    }

    pub(super) fn sequencing_confirmation_sync_target_selection(
        &mut self,
        report: &SequencingConfirmationReport,
        target_id: &str,
    ) {
        let Some(target) = report
            .targets
            .iter()
            .find(|row| row.target_id.eq_ignore_ascii_case(target_id.trim()))
        else {
            self.sequencing_confirmation_ui.selected_target_id.clear();
            return;
        };
        self.sequencing_confirmation_ui.selected_review_focus_kind =
            Some(SequencingConfirmationReviewFocusKind::Target);
        self.sequencing_confirmation_ui.selected_gap_start_0based = None;
        self.sequencing_confirmation_ui
            .selected_gap_end_0based_exclusive = None;
        self.sequencing_confirmation_ui.selected_target_id = target.target_id.clone();
        let preferred_variant = report.variants.iter().find(|row| {
            row.target_id
                .as_deref()
                .is_some_and(|id| id.eq_ignore_ascii_case(target.target_id.as_str()))
        });
        let preferred_evidence_id = target
            .contradicting_read_ids
            .first()
            .or(target.support_read_ids.first())
            .cloned()
            .or_else(|| {
                report
                    .reads
                    .iter()
                    .find(|row| {
                        row.covered_target_ids
                            .iter()
                            .any(|id| id.eq_ignore_ascii_case(target.target_id.as_str()))
                    })
                    .map(|row| row.evidence_id.clone())
            });
        if let Some(evidence_id) = preferred_evidence_id.as_deref() {
            self.sequencing_confirmation_sync_evidence_selection(report, evidence_id);
            self.sequencing_confirmation_ui.selected_target_id = target.target_id.clone();
            self.sequencing_confirmation_ui.selected_review_focus_kind =
                Some(SequencingConfirmationReviewFocusKind::Target);
        }
        if let Some(variant) = preferred_variant {
            self.sequencing_confirmation_ui.selected_variant_id = variant.variant_id.clone();
            if let Some(trace_id) = variant.trace_id.as_deref() {
                self.sequencing_confirmation_ui.selected_trace_id = trace_id.to_string();
            }
        }
    }

    pub(super) fn sequencing_confirmation_apply_overview_selection(
        &mut self,
        report: &SequencingConfirmationReport,
        selection: SequencingConfirmationOverviewSelection,
    ) -> bool {
        match selection {
            SequencingConfirmationOverviewSelection::Target(target_id) => {
                self.sequencing_confirmation_sync_target_selection(report, &target_id);
                true
            }
            SequencingConfirmationOverviewSelection::Evidence(evidence_id) => {
                self.sequencing_confirmation_sync_evidence_selection(report, &evidence_id);
                true
            }
            SequencingConfirmationOverviewSelection::Variant(variant_id) => {
                self.sequencing_confirmation_sync_variant_selection(report, &variant_id);
                true
            }
            SequencingConfirmationOverviewSelection::CoverageGap(
                start_0based,
                end_0based_exclusive,
            ) => {
                self.sequencing_confirmation_ui.selected_review_focus_kind =
                    Some(SequencingConfirmationReviewFocusKind::CoverageGap);
                self.sequencing_confirmation_ui.selected_gap_start_0based = Some(start_0based);
                self.sequencing_confirmation_ui
                    .selected_gap_end_0based_exclusive = Some(end_0based_exclusive);
                if let Some(target_id) = Self::sequencing_confirmation_select_target_for_gap(
                    report,
                    start_0based,
                    end_0based_exclusive,
                ) {
                    self.sequencing_confirmation_sync_target_selection(report, &target_id);
                    self.sequencing_confirmation_ui.selected_gap_start_0based = Some(start_0based);
                    self.sequencing_confirmation_ui
                        .selected_gap_end_0based_exclusive = Some(end_0based_exclusive);
                    true
                } else {
                    true
                }
            }
        }
    }

    pub(super) fn sequencing_confirmation_unresolved_review_queue(
        report: &SequencingConfirmationReport,
        sequence_length: usize,
    ) -> Vec<SequencingConfirmationOverviewSelection> {
        let mut queue = Self::sequencing_confirmation_target_review_rows(report, true)
            .into_iter()
            .filter(|row| row.status != SequencingConfirmationStatus::Confirmed)
            .map(|row| SequencingConfirmationOverviewSelection::Target(row.target_id.clone()))
            .collect::<Vec<_>>();
        queue.extend(
            Self::sequencing_confirmation_variant_review_rows(report, true)
                .into_iter()
                .filter(|row| row.status != SequencingConfirmationStatus::Confirmed)
                .map(|row| {
                    SequencingConfirmationOverviewSelection::Variant(row.variant_id.clone())
                }),
        );
        queue.extend(
            Self::sequencing_confirmation_uncovered_spans(report, sequence_length)
                .into_iter()
                .map(|(start, end)| {
                    SequencingConfirmationOverviewSelection::CoverageGap(start, end)
                }),
        );
        queue
    }

    pub(super) fn sequencing_confirmation_current_unresolved_focus(
        report: &SequencingConfirmationReport,
        selected_target_id: &str,
        selected_variant_id: &str,
        selected_gap: Option<(usize, usize)>,
        selected_focus_kind: Option<SequencingConfirmationReviewFocusKind>,
    ) -> Option<SequencingConfirmationOverviewSelection> {
        let unresolved_target = report
            .targets
            .iter()
            .find(|row| {
                row.status != SequencingConfirmationStatus::Confirmed
                    && row
                        .target_id
                        .eq_ignore_ascii_case(selected_target_id.trim())
            })
            .map(|row| SequencingConfirmationOverviewSelection::Target(row.target_id.clone()));
        let unresolved_variant = report
            .variants
            .iter()
            .find(|row| {
                row.status != SequencingConfirmationStatus::Confirmed
                    && row
                        .variant_id
                        .eq_ignore_ascii_case(selected_variant_id.trim())
            })
            .map(|row| SequencingConfirmationOverviewSelection::Variant(row.variant_id.clone()));
        let unresolved_gap = selected_gap
            .map(|(start, end)| SequencingConfirmationOverviewSelection::CoverageGap(start, end));
        match selected_focus_kind {
            Some(SequencingConfirmationReviewFocusKind::CoverageGap) => unresolved_gap,
            Some(SequencingConfirmationReviewFocusKind::Variant) => {
                unresolved_variant.or(unresolved_target).or(unresolved_gap)
            }
            Some(SequencingConfirmationReviewFocusKind::Target) => {
                unresolved_target.or(unresolved_variant).or(unresolved_gap)
            }
            Some(SequencingConfirmationReviewFocusKind::Evidence) => {
                unresolved_variant.or(unresolved_target).or(unresolved_gap)
            }
            None => unresolved_gap.or(unresolved_variant).or(unresolved_target),
        }
    }

    pub(super) fn sequencing_confirmation_unresolved_focus_label(
        report: &SequencingConfirmationReport,
        focus: &SequencingConfirmationOverviewSelection,
    ) -> String {
        match focus {
            SequencingConfirmationOverviewSelection::Target(target_id) => report
                .targets
                .iter()
                .find(|row| row.target_id.eq_ignore_ascii_case(target_id.trim()))
                .map(|row| format!("Target: {}", row.label))
                .unwrap_or_else(|| format!("Target: {target_id}")),
            SequencingConfirmationOverviewSelection::Variant(variant_id) => report
                .variants
                .iter()
                .find(|row| row.variant_id.eq_ignore_ascii_case(variant_id.trim()))
                .map(|row| format!("Variant: {}", row.label))
                .unwrap_or_else(|| format!("Variant: {variant_id}")),
            SequencingConfirmationOverviewSelection::CoverageGap(start, end) => {
                format!("Coverage gap: {start}..{end}")
            }
            SequencingConfirmationOverviewSelection::Evidence(evidence_id) => {
                format!("Evidence: {evidence_id}")
            }
        }
    }

    pub(super) fn sequencing_confirmation_step_unresolved_focus(
        &mut self,
        report: &SequencingConfirmationReport,
        sequence_length: usize,
        direction: isize,
    ) -> bool {
        let queue = Self::sequencing_confirmation_unresolved_review_queue(report, sequence_length);
        if queue.is_empty() {
            return false;
        }
        let current_focus = Self::sequencing_confirmation_current_unresolved_focus(
            report,
            &self.sequencing_confirmation_ui.selected_target_id,
            &self.sequencing_confirmation_ui.selected_variant_id,
            Self::sequencing_confirmation_selected_coverage_gap(
                report,
                sequence_length,
                self.sequencing_confirmation_ui.selected_gap_start_0based,
                self.sequencing_confirmation_ui
                    .selected_gap_end_0based_exclusive,
            ),
            self.sequencing_confirmation_ui.selected_review_focus_kind,
        );
        let target_index = match current_focus
            .as_ref()
            .and_then(|focus| queue.iter().position(|row| row == focus))
        {
            Some(index) if direction < 0 => index.saturating_sub(1),
            Some(index) if direction > 0 => (index + 1).min(queue.len().saturating_sub(1)),
            Some(index) => index,
            None if direction < 0 => queue.len().saturating_sub(1),
            None => 0,
        };
        self.sequencing_confirmation_apply_overview_selection(report, queue[target_index].clone())
    }

    pub(super) fn render_sequencing_confirmation_construct_overview(
        ui: &mut egui::Ui,
        report: &SequencingConfirmationReport,
        sequence_length: usize,
        selected_target_id: &str,
        selected_evidence_id: &str,
        selected_variant_id: &str,
        selected_gap: Option<(usize, usize)>,
    ) -> Option<SequencingConfirmationOverviewSelection> {
        if sequence_length == 0 {
            ui.small(
                "Construct overview becomes available once the expected construct length is known.",
            );
            return None;
        }
        let mut selection = None;
        ui.small(
            "Click a target span, evidence span, coverage gap, or variant marker to sync the detailed review panes below.",
        );
        ui.horizontal_wrapped(|ui| {
            ui.small(format!("1"));
            ui.separator();
            ui.small(format!("{} bp", sequence_length));
            if sequence_length > 1 {
                ui.separator();
                ui.small(format!("midpoint {}", sequence_length / 2));
            }
        });

        let label_width = 108.0;

        ui.horizontal(|ui| {
            ui.add_sized([label_width, 18.0], egui::Label::new("Targets"));
            let desired_size = egui::vec2(ui.available_width().max(260.0), 20.0);
            let (rect, response) = ui.allocate_exact_size(desired_size, egui::Sense::click());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 3.0, egui::Color32::from_rgb(248, 250, 252));
            painter.rect_stroke(
                rect,
                3.0,
                egui::Stroke::new(1.0, egui::Color32::from_gray(180)),
                egui::StrokeKind::Outside,
            );
            for target in &report.targets {
                let span_rect = Self::sequencing_confirmation_overview_span_rect(
                    rect,
                    sequence_length,
                    target.start_0based,
                    target.end_0based_exclusive,
                )
                .shrink2(egui::vec2(0.0, 3.0));
                let selected = target
                    .target_id
                    .eq_ignore_ascii_case(selected_target_id.trim());
                painter.rect_filled(
                    span_rect,
                    2.0,
                    Self::sequencing_confirmation_status_color(target.status)
                        .gamma_multiply(if selected { 1.0 } else { 0.75 }),
                );
                if selected {
                    painter.rect_stroke(
                        span_rect.expand(1.0),
                        2.0,
                        egui::Stroke::new(1.5, egui::Color32::BLACK),
                        egui::StrokeKind::Outside,
                    );
                }
            }
            if response.clicked()
                && let Some(pointer) = response.interact_pointer_pos()
            {
                for target in report.targets.iter().rev() {
                    let span_rect = Self::sequencing_confirmation_overview_span_rect(
                        rect,
                        sequence_length,
                        target.start_0based,
                        target.end_0based_exclusive,
                    )
                    .shrink2(egui::vec2(0.0, 3.0));
                    if span_rect.contains(pointer) {
                        selection = Some(SequencingConfirmationOverviewSelection::Target(
                            target.target_id.clone(),
                        ));
                        break;
                    }
                }
            }
        });

        ui.horizontal(|ui| {
            ui.add_sized([label_width, 18.0], egui::Label::new("Evidence spans"));
            let desired_size = egui::vec2(ui.available_width().max(260.0), 24.0);
            let (rect, response) = ui.allocate_exact_size(desired_size, egui::Sense::click());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 3.0, egui::Color32::from_rgb(248, 250, 252));
            painter.rect_stroke(
                rect,
                3.0,
                egui::Stroke::new(1.0, egui::Color32::from_gray(180)),
                egui::StrokeKind::Outside,
            );
            for (idx, read) in report.reads.iter().enumerate() {
                let span_rect = Self::sequencing_confirmation_overview_span_rect(
                    rect,
                    sequence_length,
                    read.best_alignment.aligned_target_start_0based,
                    read.best_alignment.aligned_target_end_0based_exclusive,
                );
                let vertical_offset = (idx % 3) as f32 * 5.0;
                let bar_rect = egui::Rect::from_min_max(
                    egui::pos2(span_rect.left(), rect.top() + 3.0 + vertical_offset),
                    egui::pos2(span_rect.right(), rect.top() + 7.0 + vertical_offset),
                );
                let selected = read
                    .evidence_id
                    .eq_ignore_ascii_case(selected_evidence_id.trim());
                painter.rect_filled(
                    bar_rect,
                    1.0,
                    Self::sequencing_confirmation_evidence_color(read, selected),
                );
                if selected {
                    painter.rect_stroke(
                        bar_rect.expand(1.0),
                        1.0,
                        egui::Stroke::new(1.0, egui::Color32::BLACK),
                        egui::StrokeKind::Outside,
                    );
                }
            }
            if response.clicked()
                && let Some(pointer) = response.interact_pointer_pos()
            {
                let mut best_match: Option<(&SequencingConfirmationReadResult, f32)> = None;
                for read in &report.reads {
                    let span_rect = Self::sequencing_confirmation_overview_span_rect(
                        rect,
                        sequence_length,
                        read.best_alignment.aligned_target_start_0based,
                        read.best_alignment.aligned_target_end_0based_exclusive,
                    );
                    if span_rect.left() <= pointer.x && pointer.x <= span_rect.right() {
                        let midpoint = (span_rect.left() + span_rect.right()) / 2.0;
                        let distance = (pointer.x - midpoint).abs();
                        match best_match {
                            Some((_, best_distance)) if best_distance <= distance => {}
                            _ => best_match = Some((read, distance)),
                        }
                    }
                }
                if let Some((read, _)) = best_match {
                    selection = Some(SequencingConfirmationOverviewSelection::Evidence(
                        read.evidence_id.clone(),
                    ));
                }
            }
        });

        ui.horizontal(|ui| {
            ui.add_sized([label_width, 18.0], egui::Label::new("Coverage gaps"));
            let desired_size = egui::vec2(ui.available_width().max(260.0), 20.0);
            let (rect, response) = ui.allocate_exact_size(desired_size, egui::Sense::click());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 3.0, egui::Color32::from_rgb(235, 248, 239));
            painter.rect_stroke(
                rect,
                3.0,
                egui::Stroke::new(1.0, egui::Color32::from_gray(180)),
                egui::StrokeKind::Outside,
            );
            let uncovered_spans =
                Self::sequencing_confirmation_uncovered_spans(report, sequence_length);
            for (start, end) in &uncovered_spans {
                let gap_rect = Self::sequencing_confirmation_overview_span_rect(
                    rect,
                    sequence_length,
                    *start,
                    *end,
                )
                .shrink2(egui::vec2(0.0, 4.0));
                painter.rect_filled(
                    gap_rect,
                    1.0,
                    egui::Color32::from_rgb(180, 83, 9).gamma_multiply(0.9),
                );
                if selected_gap == Some((*start, *end)) {
                    painter.rect_stroke(
                        gap_rect.expand(1.0),
                        1.0,
                        egui::Stroke::new(1.5, egui::Color32::BLACK),
                        egui::StrokeKind::Outside,
                    );
                }
            }
            if uncovered_spans.is_empty() {
                painter.text(
                    rect.center(),
                    egui::Align2::CENTER_CENTER,
                    "full evidence coverage",
                    egui::TextStyle::Small.resolve(ui.style()),
                    egui::Color32::from_rgb(21, 128, 61),
                );
            }
            if response.clicked()
                && let Some(pointer) = response.interact_pointer_pos()
            {
                for (start, end) in uncovered_spans.iter().rev() {
                    let gap_rect = Self::sequencing_confirmation_overview_span_rect(
                        rect,
                        sequence_length,
                        *start,
                        *end,
                    )
                    .shrink2(egui::vec2(0.0, 4.0));
                    if gap_rect.contains(pointer) {
                        selection = Some(SequencingConfirmationOverviewSelection::CoverageGap(
                            *start, *end,
                        ));
                        break;
                    }
                }
            }
        });

        ui.horizontal(|ui| {
            ui.add_sized([label_width, 18.0], egui::Label::new("Variants"));
            let desired_size = egui::vec2(ui.available_width().max(260.0), 24.0);
            let (rect, response) = ui.allocate_exact_size(desired_size, egui::Sense::click());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 3.0, egui::Color32::from_rgb(248, 250, 252));
            painter.rect_stroke(
                rect,
                3.0,
                egui::Stroke::new(1.0, egui::Color32::from_gray(180)),
                egui::StrokeKind::Outside,
            );
            for variant in &report.variants {
                let center = Self::sequencing_confirmation_overview_x(
                    rect,
                    sequence_length,
                    variant.start_0based,
                );
                let selected = variant
                    .variant_id
                    .eq_ignore_ascii_case(selected_variant_id.trim());
                painter.line_segment(
                    [
                        egui::pos2(center, rect.top() + 4.0),
                        egui::pos2(center, rect.bottom() - 4.0),
                    ],
                    egui::Stroke::new(
                        if selected { 3.0 } else { 2.0 },
                        Self::sequencing_confirmation_variant_color(variant.classification),
                    ),
                );
                if selected {
                    painter.circle_filled(
                        egui::pos2(center, rect.center().y),
                        3.5,
                        egui::Color32::BLACK,
                    );
                }
            }
            if response.clicked()
                && let Some(pointer) = response.interact_pointer_pos()
            {
                let mut best_match: Option<(&SequencingConfirmationVariantRow, f32)> = None;
                for variant in &report.variants {
                    let center = Self::sequencing_confirmation_overview_x(
                        rect,
                        sequence_length,
                        variant.start_0based,
                    );
                    let distance = (pointer.x - center).abs();
                    if distance <= 8.0 {
                        match best_match {
                            Some((_, best_distance)) if best_distance <= distance => {}
                            _ => best_match = Some((variant, distance)),
                        }
                    }
                }
                if let Some((variant, _)) = best_match {
                    selection = Some(SequencingConfirmationOverviewSelection::Variant(
                        variant.variant_id.clone(),
                    ));
                }
            }
        });

        ui.horizontal(|ui| {
            ui.add_sized(
                [label_width, 18.0],
                egui::Label::new("Selected discrepancies"),
            );
            let desired_size = egui::vec2(ui.available_width().max(260.0), 20.0);
            let (rect, _response) = ui.allocate_exact_size(desired_size, egui::Sense::hover());
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 3.0, egui::Color32::from_rgb(248, 250, 252));
            painter.rect_stroke(
                rect,
                3.0,
                egui::Stroke::new(1.0, egui::Color32::from_gray(180)),
                egui::StrokeKind::Outside,
            );
            if let Some(read) =
                Self::sequencing_confirmation_selected_read(report, selected_evidence_id)
            {
                for discrepancy in &read.discrepancies {
                    let span_rect = Self::sequencing_confirmation_overview_span_rect(
                        rect,
                        sequence_length,
                        discrepancy.target_start_0based,
                        discrepancy.target_end_0based_exclusive,
                    )
                    .shrink2(egui::vec2(0.0, 4.0));
                    let color = match discrepancy.kind {
                        crate::engine::SequencingConfirmationDiscrepancyKind::Mismatch => {
                            egui::Color32::from_rgb(230, 81, 0)
                        }
                        crate::engine::SequencingConfirmationDiscrepancyKind::Insertion => {
                            egui::Color32::from_rgb(183, 28, 28)
                        }
                        crate::engine::SequencingConfirmationDiscrepancyKind::Deletion => {
                            egui::Color32::from_rgb(180, 83, 9)
                        }
                    };
                    painter.rect_filled(span_rect, 1.0, color);
                }
            }
        });

        selection
    }

    pub fn render_sequencing_confirmation_specialist(
        &mut self,
        ui: &mut egui::Ui,
        _ctx: &egui::Context,
    ) {
        let expected_seq_id = self.seq_id.clone().unwrap_or_default();
        let expected_title = self.window_title();
        let selection = self.current_selection_range_0based();
        let project_sequence_ids = self.project_sequence_ids();
        let mut report_summaries = self.sequencing_confirmation_report_summaries();
        let latest_report_id = report_summaries.last().map(|row| row.report_id.clone());
        if self.sequencing_confirmation_ui.report_id.trim().is_empty() {
            self.sequencing_confirmation_ui.report_id =
                self.sequencing_confirmation_default_report_id();
        }
        if self
            .sequencing_confirmation_ui
            .selected_report_id
            .trim()
            .is_empty()
            && let Some(report_id) = latest_report_id.clone()
        {
            self.sequencing_confirmation_ui.selected_report_id = report_id;
        }
        report_summaries = self.sequencing_confirmation_report_summaries();
        let selected_report = self
            .sequencing_confirmation_ui
            .selected_report_id
            .trim()
            .to_string();
        let selected_report = if selected_report.is_empty() {
            None
        } else {
            self.engine
                .as_ref()
                .and_then(|engine| engine.read().ok())
                .and_then(|guard| {
                    guard
                        .get_sequencing_confirmation_report(&selected_report)
                        .ok()
                })
        };
        let selected_variant_id_missing = self
            .sequencing_confirmation_ui
            .selected_variant_id
            .trim()
            .is_empty()
            || !selected_report.as_ref().is_some_and(|report| {
                report.variants.iter().any(|row| {
                    row.variant_id.eq_ignore_ascii_case(
                        self.sequencing_confirmation_ui.selected_variant_id.trim(),
                    )
                })
            });
        if selected_variant_id_missing
            && let Some(row) = selected_report
                .as_ref()
                .and_then(|report| report.variants.first())
        {
            self.sequencing_confirmation_ui.selected_variant_id = row.variant_id.clone();
            if let Some(trace_id) = row.trace_id.as_deref() {
                self.sequencing_confirmation_ui.selected_trace_id = trace_id.to_string();
            }
        }
        let selected_variant = selected_report.as_ref().and_then(|report| {
            report.variants.iter().find(|row| {
                row.variant_id.eq_ignore_ascii_case(
                    self.sequencing_confirmation_ui.selected_variant_id.trim(),
                )
            })
        });
        let selected_target_id_missing = self
            .sequencing_confirmation_ui
            .selected_target_id
            .trim()
            .is_empty()
            || !selected_report.as_ref().is_some_and(|report| {
                report.targets.iter().any(|row| {
                    row.target_id.eq_ignore_ascii_case(
                        self.sequencing_confirmation_ui.selected_target_id.trim(),
                    )
                })
            });
        if selected_target_id_missing {
            if let Some(target_id) = selected_report.as_ref().and_then(|report| {
                Self::sequencing_confirmation_selected_target_id(
                    report,
                    self.sequencing_confirmation_ui.selected_target_id.trim(),
                    self.sequencing_confirmation_ui.review_unresolved_first,
                )
            }) {
                self.sequencing_confirmation_ui.selected_target_id = target_id;
            } else {
                self.sequencing_confirmation_ui.selected_target_id.clear();
            }
        }
        let selected_target = selected_report.as_ref().and_then(|report| {
            Self::sequencing_confirmation_selected_target(
                report,
                &self.sequencing_confirmation_ui.selected_target_id,
                self.sequencing_confirmation_ui.review_unresolved_first,
            )
        });
        let selected_evidence_id_missing = self
            .sequencing_confirmation_ui
            .selected_evidence_id
            .trim()
            .is_empty()
            || !selected_report.as_ref().is_some_and(|report| {
                report.reads.iter().any(|row| {
                    row.evidence_id.eq_ignore_ascii_case(
                        self.sequencing_confirmation_ui.selected_evidence_id.trim(),
                    )
                })
            });
        if selected_evidence_id_missing {
            if let Some(evidence_id) = selected_report.as_ref().and_then(|report| {
                Self::sequencing_confirmation_selected_evidence_id(
                    report,
                    self.sequencing_confirmation_ui.selected_evidence_id.trim(),
                )
            }) {
                self.sequencing_confirmation_ui.selected_evidence_id = evidence_id;
            } else {
                self.sequencing_confirmation_ui.selected_evidence_id.clear();
            }
        }
        let selected_read = selected_report.as_ref().and_then(|report| {
            Self::sequencing_confirmation_selected_read(
                report,
                &self.sequencing_confirmation_ui.selected_evidence_id,
            )
        });
        let trace_summaries = self.sequencing_trace_summaries();
        let selected_trace_id_missing = self
            .sequencing_confirmation_ui
            .selected_trace_id
            .trim()
            .is_empty()
            || !trace_summaries.iter().any(|row| {
                row.trace_id
                    .eq_ignore_ascii_case(self.sequencing_confirmation_ui.selected_trace_id.trim())
            });
        if selected_trace_id_missing {
            if let Some(trace_id) = selected_report
                .as_ref()
                .and_then(|report| report.trace_ids.first())
                .cloned()
            {
                self.sequencing_confirmation_ui.selected_trace_id = trace_id;
            } else if let Some(row) = trace_summaries.first() {
                self.sequencing_confirmation_ui.selected_trace_id = row.trace_id.clone();
            }
        }
        let selected_trace =
            self.sequencing_trace_record(self.sequencing_confirmation_ui.selected_trace_id.trim());
        let selected_variant_trace = selected_variant
            .and_then(|row| row.trace_id.as_deref())
            .and_then(|trace_id| self.sequencing_trace_record(trace_id));
        let expected_sequence_length = self.dna.read().map(|dna| dna.len()).unwrap_or(0);
        let report_sequence_length = selected_report
            .as_ref()
            .map(|report| {
                Self::sequencing_confirmation_report_sequence_length(
                    report,
                    expected_sequence_length,
                )
            })
            .unwrap_or(expected_sequence_length);
        let selected_gap = selected_report.as_ref().and_then(|report| {
            Self::sequencing_confirmation_selected_coverage_gap(
                report,
                report_sequence_length,
                self.sequencing_confirmation_ui.selected_gap_start_0based,
                self.sequencing_confirmation_ui
                    .selected_gap_end_0based_exclusive,
            )
        });
        if let Some((start_0based, end_0based_exclusive)) = selected_gap {
            self.sequencing_confirmation_ui.selected_gap_start_0based = Some(start_0based);
            self.sequencing_confirmation_ui
                .selected_gap_end_0based_exclusive = Some(end_0based_exclusive);
        } else {
            self.sequencing_confirmation_ui.selected_gap_start_0based = None;
            self.sequencing_confirmation_ui
                .selected_gap_end_0based_exclusive = None;
            if self.sequencing_confirmation_ui.selected_review_focus_kind
                == Some(SequencingConfirmationReviewFocusKind::CoverageGap)
            {
                self.sequencing_confirmation_ui.selected_review_focus_kind = None;
            }
        }
        if let Some(trace) = selected_trace.as_ref() {
            let clamped = Self::sequencing_trace_clamp_base_index(
                trace,
                self.sequencing_confirmation_ui
                    .chromatogram_trace_base_index,
            );
            if clamped
                != self
                    .sequencing_confirmation_ui
                    .chromatogram_trace_base_index
            {
                self.sequencing_confirmation_ui
                    .chromatogram_trace_base_index = clamped;
            }
        } else {
            self.sequencing_confirmation_ui
                .chromatogram_trace_base_index = 0;
        }

        ui.label(
            "Construct-confirmation specialist for already-loaded sequencing reads and imported trace evidence. Target entry uses the same shared engine contract as `seq-confirm run`.",
        );
        ui.small(
            "Current GUI scope covers called reads plus imported ABI/AB1/SCF trace IDs through the same shared report store, including baseline-aware variant classification and chromatogram curve review for trace-backed loci.",
        );
        ui.separator();
        ui.columns(2, |columns| {
            columns[0].heading("Inputs + Run");
            columns[0].small(format!("Expected construct: {expected_title}"));
            columns[0].small(format!("Expected sequence ID: {expected_seq_id}"));
            columns[0].label("Baseline/reference sequence ID");
            let baseline_changed = columns[0]
                .add(
                    egui::TextEdit::singleline(
                        &mut self.sequencing_confirmation_ui.baseline_seq_id_text,
                    )
                    .desired_width(320.0),
                )
                .on_hover_text(
                    "Optional baseline sequence used to classify intended edits versus reference reversions. Leave empty for expected-only review.",
                )
                .changed();
            if baseline_changed {
                self.save_engine_ops_state();
            }
            columns[0].add_space(6.0);
            columns[0].label("Read sequence IDs");
            let read_ids_changed = columns[0]
                .add(
                    egui::TextEdit::singleline(
                        &mut self.sequencing_confirmation_ui.read_seq_ids_text,
                    )
                    .desired_width(320.0),
                )
                .on_hover_text(
                    "Comma-separated read sequence IDs already loaded in the current project state.",
                )
                .changed();
            if read_ids_changed {
                self.save_engine_ops_state();
            }
            columns[0].label("Imported trace IDs");
            let trace_ids_changed = columns[0]
                .add(
                    egui::TextEdit::singleline(
                        &mut self.sequencing_confirmation_ui.trace_ids_text,
                    )
                    .desired_width(320.0),
                )
                .on_hover_text(
                    "Comma-separated imported sequencing-trace IDs already present in the project evidence store.",
                )
                .changed();
            if trace_ids_changed {
                self.save_engine_ops_state();
            }
            columns[0].group(|ui| {
                ui.label(egui::RichText::new("Raw Trace Import").strong());
                ui.small(
                    "Import ABI/AB1/SCF files directly into the shared sequencing-trace evidence store used by GUI, CLI, and shell review.",
                );
                ui.horizontal_wrapped(|ui| {
                    ui.label("trace file");
                    if ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.sequencing_confirmation_ui.trace_import_path,
                            )
                            .desired_width(240.0),
                        )
                        .on_hover_text(
                            "Local ABI/AB1/SCF file to import into the project evidence store.",
                        )
                        .changed()
                    {
                        self.save_engine_ops_state();
                    }
                    if ui
                        .button("Browse...")
                        .on_hover_text(
                            "Choose a local ABI/AB1/SCF chromatogram file and fill the import path.",
                        )
                        .clicked()
                        && let Some(path) = rfd::FileDialog::new()
                            .add_filter("Sequencing traces", &["ab1", "abi", "scf"])
                            .pick_file()
                    {
                        self.sequencing_confirmation_ui.trace_import_path =
                            path.display().to_string();
                        self.save_engine_ops_state();
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("trace id");
                    if ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.sequencing_confirmation_ui.trace_import_id,
                            )
                            .desired_width(180.0),
                        )
                        .on_hover_text(
                            "Optional stable trace id. Leave empty to auto-derive one from the filename.",
                        )
                        .changed()
                    {
                        self.save_engine_ops_state();
                    }
                    let associate_changed = ui
                        .checkbox(
                            &mut self
                                .sequencing_confirmation_ui
                                .trace_import_associate_with_expected_seq,
                            "associate with expected construct",
                        )
                        .on_hover_text(
                            "Record the active expected construct sequence ID on the imported trace metadata.",
                        )
                        .changed();
                    let auto_add_changed = ui
                        .checkbox(
                            &mut self.sequencing_confirmation_ui.trace_import_add_to_run,
                            "add imported trace to current run",
                        )
                        .on_hover_text(
                            "Append the imported trace ID to the current trace input list after a successful import.",
                        )
                        .changed();
                    if associate_changed || auto_add_changed {
                        self.save_engine_ops_state();
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    if ui
                        .button("Import trace")
                        .on_hover_text(
                            "Run ImportSequencingTrace through the shared engine path, then focus the imported trace in this specialist.",
                        )
                        .clicked()
                    {
                        match self.build_import_sequencing_trace_operation() {
                            Ok(op) => {
                                if let Some(result) = self.apply_operation_with_feedback_and_result(op)
                                {
                                    self.handle_imported_sequencing_trace_result(&result);
                                }
                            }
                            Err(err) => {
                                self.op_status = err;
                            }
                        }
                    }
                    ui.small(
                        "Imported traces stay separate from confirmation verdicts until you explicitly include their IDs in a run.",
                    );
                });
            });
            let alternate_ids = project_sequence_ids
                .iter()
                .filter(|seq_id| *seq_id != &expected_seq_id)
                .take(10)
                .cloned()
                .collect::<Vec<_>>();
            if !alternate_ids.is_empty() {
                columns[0].small(format!(
                    "Available project sequence IDs: {}{}",
                    alternate_ids.join(", "),
                    if project_sequence_ids.len() > alternate_ids.len() {
                        ", ..."
                    } else {
                        ""
                    }
                ));
            }
            if trace_summaries.is_empty() {
                columns[0].small(
                    "No imported sequencing traces are stored yet. Use the Raw Trace Import box above or `seq-trace import ...` through CLI/shared shell, then select them here by trace ID.",
                );
            } else {
                let preview_ids = trace_summaries
                    .iter()
                    .take(8)
                    .map(|row| row.trace_id.clone())
                    .collect::<Vec<_>>();
                let suffix = if trace_summaries.len() > preview_ids.len() {
                    ", ..."
                } else {
                    ""
                };
                columns[0].small(format!(
                    "Imported traces available: {} total [{}{}]",
                    trace_summaries.len(),
                    preview_ids.join(", "),
                    suffix
                ));
            }
            columns[0].separator();
            columns[0].label("Targets");
            let include_full_span_changed = columns[0]
                .checkbox(
                    &mut self.sequencing_confirmation_ui.include_full_span_target,
                    "Include full construct span target",
                )
                .on_hover_text(
                    "Keep the default full-span confirmation target in addition to any explicit junction checkpoints.",
                )
                .changed();
            if include_full_span_changed {
                self.save_engine_ops_state();
            }
            columns[0].horizontal_wrapped(|ui| {
                ui.label("Junction breakpoints (0-based)");
                let changed = ui
                    .add(
                        egui::TextEdit::singleline(
                            &mut self.sequencing_confirmation_ui.junction_positions_0based,
                        )
                        .desired_width(180.0),
                    )
                    .on_hover_text(
                        "Comma-separated left-end positions for explicit junction targets; mirrors `seq-confirm run --junction`.",
                    )
                    .changed();
                ui.label("flank");
                let flank_changed = ui
                    .add(
                        egui::TextEdit::singleline(
                            &mut self.sequencing_confirmation_ui.junction_flank_bp,
                        )
                        .desired_width(56.0),
                    )
                    .on_hover_text(
                        "Number of bp to include on each side of a junction breakpoint.",
                    )
                    .changed();
                if changed || flank_changed {
                    self.save_engine_ops_state();
                }
            });
            columns[0].horizontal_wrapped(|ui| {
                if let Some((start, end_exclusive)) = selection {
                    let mid = start + (end_exclusive.saturating_sub(start) / 2);
                    ui.small(format!(
                        "Active selection: {start}..{end_exclusive} (midpoint {mid})"
                    ));
                    if ui
                        .small_button("Add start")
                        .on_hover_text(
                            "Append the current selection start as a junction breakpoint.",
                        )
                        .clicked()
                    {
                        Self::append_unique_csv_token(
                            &mut self.sequencing_confirmation_ui.junction_positions_0based,
                            start.to_string(),
                        );
                        self.save_engine_ops_state();
                    }
                    if ui
                        .small_button("Add midpoint")
                        .on_hover_text(
                            "Append the current selection midpoint as a junction breakpoint.",
                        )
                        .clicked()
                    {
                        Self::append_unique_csv_token(
                            &mut self.sequencing_confirmation_ui.junction_positions_0based,
                            mid.to_string(),
                        );
                        self.save_engine_ops_state();
                    }
                    if ui
                        .small_button("Add end")
                        .on_hover_text(
                            "Append the current selection end as a junction breakpoint.",
                        )
                        .clicked()
                    {
                        Self::append_unique_csv_token(
                            &mut self.sequencing_confirmation_ui.junction_positions_0based,
                            end_exclusive.to_string(),
                        );
                        self.save_engine_ops_state();
                    }
                } else {
                    ui.small(
                        "No active selection. Select a region in the sequence window if you want one-click junction seeding.",
                    );
                }
            });
            columns[0].separator();
            columns[0].label("Alignment + thresholds");
            columns[0].horizontal(|ui| {
                ui.label("mode");
                let mut changed = false;
                egui::ComboBox::from_id_salt(("seq_confirm_alignment_mode", self.panel_scope_key()))
                    .selected_text(self.sequencing_confirmation_ui.alignment_mode.as_str())
                    .show_ui(ui, |ui| {
                        changed |= ui
                            .selectable_value(
                                &mut self.sequencing_confirmation_ui.alignment_mode,
                                PairwiseAlignmentMode::Local,
                                "local",
                            )
                            .changed();
                        changed |= ui
                            .selectable_value(
                                &mut self.sequencing_confirmation_ui.alignment_mode,
                                PairwiseAlignmentMode::Global,
                                "global",
                            )
                            .changed();
                    });
                if changed {
                    self.save_engine_ops_state();
                }
                let rc_changed = ui
                    .checkbox(
                        &mut self.sequencing_confirmation_ui.allow_reverse_complement,
                        "allow reverse complement",
                    )
                    .on_hover_text(
                        "Try the reverse-complement orientation and keep the better usable alignment.",
                    )
                    .changed();
                if rc_changed {
                    self.save_engine_ops_state();
                }
            });
            egui::Grid::new(("seq_confirm_controls", self.panel_scope_key()))
                .num_columns(4)
                .striped(true)
                .show(&mut columns[0], |ui| {
                    ui.label("match");
                    if ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.sequencing_confirmation_ui.match_score,
                            )
                            .desired_width(60.0),
                        )
                        .changed()
                    {
                        self.save_engine_ops_state();
                    }
                    ui.label("mismatch");
                    if ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.sequencing_confirmation_ui.mismatch_score,
                            )
                            .desired_width(60.0),
                        )
                        .changed()
                    {
                        self.save_engine_ops_state();
                    }
                    ui.end_row();
                    ui.label("gap open");
                    if ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.sequencing_confirmation_ui.gap_open,
                            )
                            .desired_width(60.0),
                        )
                        .changed()
                    {
                        self.save_engine_ops_state();
                    }
                    ui.label("gap extend");
                    if ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.sequencing_confirmation_ui.gap_extend,
                            )
                            .desired_width(60.0),
                        )
                        .changed()
                    {
                        self.save_engine_ops_state();
                    }
                    ui.end_row();
                    ui.label("min identity");
                    if ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.sequencing_confirmation_ui.min_identity_fraction,
                            )
                            .desired_width(72.0),
                        )
                        .changed()
                    {
                        self.save_engine_ops_state();
                    }
                    ui.label("min target coverage");
                    if ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.sequencing_confirmation_ui.min_target_coverage_fraction,
                            )
                            .desired_width(72.0),
                        )
                        .changed()
                    {
                        self.save_engine_ops_state();
                    }
                    ui.end_row();
                    ui.label("report id");
                    if ui
                        .add(
                            egui::TextEdit::singleline(&mut self.sequencing_confirmation_ui.report_id)
                                .desired_width(180.0),
                        )
                        .changed()
                    {
                        self.save_engine_ops_state();
                    }
                    ui.end_row();
                });
            columns[0].horizontal_wrapped(|ui| {
                if ui
                    .button("Run confirmation")
                    .on_hover_text(
                        "Run ConfirmConstructReads with the current read/trace inputs and refresh the selected report.",
                    )
                    .clicked()
                {
                    match self.build_confirm_construct_reads_operation() {
                        Ok((report_id, op)) => {
                            self.sequencing_confirmation_ui.report_id = report_id.clone();
                            if self.apply_operation_with_feedback_and_result(op).is_some() {
                                self.sequencing_confirmation_ui.selected_report_id = report_id.clone();
                                self.show_sequencing_confirmation_report(&report_id);
                            }
                            self.save_engine_ops_state();
                        }
                        Err(err) => {
                            self.op_status = err;
                        }
                    }
                }
                if ui
                    .button("Refresh reports")
                    .on_hover_text(
                        "List persisted sequencing-confirmation reports for this expected construct sequence.",
                    )
                    .clicked()
                {
                    self.list_sequencing_confirmation_reports();
                }
                if ui
                    .button("Show selected")
                    .on_hover_text("Refresh status from the currently selected persisted report.")
                    .clicked()
                {
                    let report_id = self
                        .sequencing_confirmation_ui
                        .selected_report_id
                        .clone();
                    self.show_sequencing_confirmation_report(&report_id);
                }
            });
            columns[0].separator();
            columns[0].heading("Sequencing-Primer Overlays");
            columns[0].small(
                "Suggest read-coverage overlays from already-loaded primer sequences and, when unresolved loci still lack a good hit, propose fresh sequencing primers on the same shared engine path.",
            );
            columns[0].label("Primer sequence IDs");
            if columns[0]
                .add(
                    egui::TextEdit::singleline(
                        &mut self.sequencing_confirmation_ui.primer_seq_ids_text,
                    )
                    .desired_width(320.0),
                )
                .on_hover_text(
                    "Comma-separated primer sequence IDs already loaded in the current project state.",
                )
                .changed()
            {
                self.save_engine_ops_state();
            }
            columns[0].small(
                "Leave primer IDs empty if you only want de novo proposals for loci flagged in the selected saved report.",
            );
            columns[0].horizontal_wrapped(|ui| {
                ui.label("min 3' anneal bp");
                if ui
                    .add(
                        egui::TextEdit::singleline(
                            &mut self.sequencing_confirmation_ui.primer_min_3prime_anneal_bp,
                        )
                        .desired_width(56.0),
                    )
                    .changed()
                {
                    self.save_engine_ops_state();
                }
                ui.label("predicted read length bp");
                if ui
                    .add(
                        egui::TextEdit::singleline(
                            &mut self.sequencing_confirmation_ui.primer_predicted_read_length_bp,
                        )
                        .desired_width(72.0),
                    )
                    .changed()
                {
                    self.save_engine_ops_state();
                }
            });
            columns[0].small(format!(
                "Coverage annotation source: {}",
                selected_report
                    .as_ref()
                    .map(|report| report.report_id.as_str())
                    .or_else(|| {
                        let trimmed = self.sequencing_confirmation_ui.selected_report_id.trim();
                        (!trimmed.is_empty()).then_some(trimmed)
                    })
                    .unwrap_or("<none>")
            ));
            columns[0].horizontal_wrapped(|ui| {
                if ui
                    .button("Suggest primers")
                    .on_hover_text(
                        "Run SuggestSequencingPrimers for the active expected construct and annotate coverage against the selected saved report when available.",
                    )
                    .clicked()
                {
                    match self.build_suggest_sequencing_primers_operation() {
                        Ok(op) => match self.apply_operation_with_feedback_and_result(op) {
                            Some(result) => {
                                self.sequencing_confirmation_ui.primer_overlay_report =
                                    result.sequencing_primer_overlay_report.clone();
                            }
                            None => {
                                self.sequencing_confirmation_ui.primer_overlay_report = None;
                            }
                        },
                        Err(err) => {
                            self.op_status = err;
                        }
                    }
                }
                if ui
                    .button("Clear overlay")
                    .on_hover_text(
                        "Discard the currently displayed primer-overlay suggestion report from this window.",
                    )
                    .clicked()
                {
                    self.sequencing_confirmation_ui.primer_overlay_report = None;
                }
            });
            if let Some(report) = self.sequencing_confirmation_ui.primer_overlay_report.as_ref() {
                columns[0].group(|ui| {
                    ui.horizontal_wrapped(|ui| {
                        ui.label(egui::RichText::new("Current overlay report").strong());
                        ui.separator();
                        ui.small(format!("primers={}", report.primer_seq_ids.len()));
                        ui.separator();
                        ui.small(format!("suggestions={}", report.suggestion_count));
                        ui.separator();
                        ui.small(format!("guidance={}", report.problem_guidance_count));
                        ui.separator();
                        ui.small(format!("proposals={}", report.proposal_count));
                        ui.separator();
                        ui.small(format!(
                            "min_3prime_anneal_bp={}",
                            report.min_3prime_anneal_bp
                        ));
                        ui.separator();
                        ui.small(format!(
                            "predicted_read_length_bp={}",
                            report.predicted_read_length_bp
                        ));
                        ui.separator();
                        ui.small(format!(
                            "report={}",
                            report.confirmation_report_id.as_deref().unwrap_or("-")
                        ));
                    });
                    if !report.problem_guidance.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Recommended Next Primers").strong());
                        egui::Grid::new((
                            "seq_confirm_primer_guidance_grid",
                            self.panel_scope_key(),
                        ))
                        .num_columns(7)
                        .striped(true)
                        .show(ui, |ui| {
                            ui.strong("Problem");
                            ui.strong("Kind");
                            ui.strong("State");
                            ui.strong("Recommended primer");
                            ui.strong("Orientation");
                            ui.strong("3' dist");
                            ui.strong("Reason");
                            ui.end_row();
                            for row in &report.problem_guidance {
                                ui.label(&row.problem_label);
                                ui.small(row.problem_kind.as_str());
                                ui.small(&row.problem_summary);
                                ui.small(
                                    row.recommended_primer_seq_id
                                        .as_deref()
                                        .unwrap_or("-"),
                                );
                                ui.small(
                                    row.recommended_orientation
                                        .map(|value| value.as_str())
                                        .unwrap_or("-"),
                                );
                                ui.small(
                                    row.recommended_three_prime_distance_bp
                                        .map(|value| value.to_string())
                                        .unwrap_or_else(|| "-".to_string()),
                                );
                                ui.small(&row.reason);
                                ui.end_row();
                            }
                        });
                    }
                    if !report.proposals.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Proposed New Primers").strong());
                        egui::Grid::new((
                            "seq_confirm_primer_proposals_grid",
                            self.panel_scope_key(),
                        ))
                        .num_columns(8)
                        .striped(true)
                        .show(ui, |ui| {
                            ui.strong("Problem");
                            ui.strong("Orientation");
                            ui.strong("Sequence");
                            ui.strong("3' dist");
                            ui.strong("Read span");
                            ui.strong("Tm");
                            ui.strong("GC");
                            ui.strong("Reason");
                            ui.end_row();
                            for row in &report.proposals {
                                ui.label(&row.problem_label);
                                ui.small(row.orientation.as_str());
                                ui.small(egui::RichText::new(&row.primer_sequence).monospace());
                                ui.small(row.three_prime_distance_bp.to_string());
                                ui.small(format!(
                                    "{}..{}",
                                    row.predicted_read_span_start_0based,
                                    row.predicted_read_span_end_0based_exclusive
                                ));
                                ui.small(format!("{:.1}", row.tm_c));
                                ui.small(format!("{:.2}", row.gc_fraction));
                                ui.small(&row.reason);
                                ui.end_row();
                            }
                        });
                    }
                    if report.suggestions.is_empty() {
                        ui.small(
                            "No existing primer overlays matched the current expected construct with the requested exact 3' anneal length.",
                        );
                    } else {
                        egui::ScrollArea::vertical()
                            .id_salt(("seq_confirm_primer_overlay_scroll", self.panel_scope_key()))
                            .max_height(220.0)
                            .show(ui, |ui| {
                                egui::Grid::new((
                                    "seq_confirm_primer_overlay_grid",
                                    self.panel_scope_key(),
                                ))
                                .num_columns(8)
                                .striped(true)
                                .show(ui, |ui| {
                                    ui.strong("Primer");
                                    ui.strong("Orientation");
                                    ui.strong("Anneal");
                                    ui.strong("Read span");
                                    ui.strong("Targets");
                                    ui.strong("Problem targets");
                                    ui.strong("Variants");
                                    ui.strong("Problem variants");
                                    ui.end_row();
                                    for row in &report.suggestions {
                                        ui.label(&row.primer_seq_id);
                                        ui.small(row.orientation.as_str());
                                        ui.small(format!(
                                            "{}..{}",
                                            row.anneal_start_0based, row.anneal_end_0based_exclusive
                                        ));
                                        ui.small(format!(
                                            "{}..{}",
                                            row.predicted_read_span_start_0based,
                                            row.predicted_read_span_end_0based_exclusive
                                        ));
                                        ui.small(if row.covered_target_ids.is_empty() {
                                            "-".to_string()
                                        } else {
                                            row.covered_target_ids.join(", ")
                                        });
                                        ui.small(if row.covered_problem_target_ids.is_empty() {
                                            "-".to_string()
                                        } else {
                                            row.covered_problem_target_ids.join(", ")
                                        });
                                        ui.small(if row.covered_variant_ids.is_empty() {
                                            "-".to_string()
                                        } else {
                                            row.covered_variant_ids.join(", ")
                                        });
                                        ui.small(if row.covered_problem_variant_ids.is_empty() {
                                            "-".to_string()
                                        } else {
                                            row.covered_problem_variant_ids.join(", ")
                                        });
                                        ui.end_row();
                                    }
                                });
                            });
                    }
                });
            }

            columns[1].heading("Imported Trace Review");
            if trace_summaries.is_empty() {
                columns[1].small(
                    "No imported sequencing traces are stored in this project yet.",
                );
            } else {
                columns[1].horizontal_wrapped(|ui| {
                    ui.label("Selected trace");
                    let mut trace_selection_changed = false;
                    egui::ComboBox::from_id_salt((
                        "seq_confirm_selected_trace",
                        self.panel_scope_key(),
                    ))
                    .selected_text(
                        trace_summaries
                            .iter()
                            .find(|row| {
                                row.trace_id.eq_ignore_ascii_case(
                                    self.sequencing_confirmation_ui.selected_trace_id.trim(),
                                )
                            })
                            .map(|row| {
                                format!(
                                    "{} | {} | bases={} | seq={}",
                                    row.trace_id,
                                    row.format.as_str(),
                                    row.called_base_count,
                                    row.seq_id.as_deref().unwrap_or("-")
                                )
                            })
                            .unwrap_or_else(|| "<select trace>".to_string()),
                    )
                    .show_ui(ui, |ui| {
                        for row in &trace_summaries {
                            trace_selection_changed |= ui
                                .selectable_value(
                                    &mut self.sequencing_confirmation_ui.selected_trace_id,
                                    row.trace_id.clone(),
                                    format!(
                                        "{} | {} | bases={} | seq={}",
                                        row.trace_id,
                                        row.format.as_str(),
                                        row.called_base_count,
                                        row.seq_id.as_deref().unwrap_or("-")
                                    ),
                                )
                                .changed();
                        }
                    });
                    if trace_selection_changed {
                        self.save_engine_ops_state();
                    }
                    if ui
                        .small_button("Add selected trace to run")
                        .on_hover_text(
                            "Append the selected imported trace ID to the current confirmation input list.",
                        )
                        .clicked()
                    {
                        Self::append_unique_csv_token(
                            &mut self.sequencing_confirmation_ui.trace_ids_text,
                            self.sequencing_confirmation_ui.selected_trace_id.clone(),
                        );
                        self.save_engine_ops_state();
                    }
                    if ui
                        .small_button("Use selected only")
                        .on_hover_text(
                            "Replace the current imported trace input list with the selected trace ID.",
                        )
                        .clicked()
                    {
                        self.sequencing_confirmation_ui.trace_ids_text = self
                            .sequencing_confirmation_ui
                            .selected_trace_id
                            .trim()
                            .to_string();
                        self.save_engine_ops_state();
                    }
                });
                if let Some(trace) = selected_trace.as_ref() {
                    columns[1].group(|ui| {
                        ui.horizontal_wrapped(|ui| {
                            ui.label(egui::RichText::new(&trace.trace_id).strong());
                            ui.separator();
                            ui.small(trace.format.as_str());
                            ui.separator();
                            ui.small(format!(
                                "seq={}",
                                trace.seq_id.as_deref().unwrap_or("-")
                            ));
                            ui.separator();
                            ui.small(format!("bases={}", trace.called_bases.len()));
                            ui.separator();
                            ui.small(format!(
                                "confidences={}",
                                trace.called_base_confidence_values.len()
                            ));
                            ui.separator();
                            ui.small(format!("peaks={}", trace.peak_locations.len()));
                            ui.separator();
                            ui.small(format!("channels={}", trace.channel_summaries.len()));
                            ui.separator();
                            ui.small(format!(
                                "curves={}",
                                if trace.channel_data.is_empty() { "no" } else { "yes" }
                            ));
                        });
                        ui.small(format!("source={}", trace.source_path));
                        ui.small(format!("imported_at_unix_ms={}", trace.imported_at_unix_ms));
                        ui.small(format!(
                            "sample={} well={} run={} machine={} ({})",
                            trace.sample_name.as_deref().unwrap_or("-"),
                            trace.sample_well.as_deref().unwrap_or("-"),
                            trace.run_name.as_deref().unwrap_or("-"),
                            trace.machine_name.as_deref().unwrap_or("-"),
                            trace.machine_model.as_deref().unwrap_or("-")
                        ));
                        if !trace.channel_summaries.is_empty() {
                            let channel_text = trace
                                .channel_summaries
                                .iter()
                                .map(|row| {
                                    format!(
                                        "{}:{}:{}",
                                        row.trace_set, row.channel, row.point_count
                                    )
                                })
                                .collect::<Vec<_>>()
                                .join(", ");
                            ui.small(format!("channels: {channel_text}"));
                        }
                        ui.small(Self::sequencing_trace_confidence_summary(
                            &trace.called_base_confidence_values,
                        ));
                        ui.small(Self::sequencing_trace_peak_summary(&trace.peak_locations));
                        if let Some(start) = trace.clip_start_base_index {
                            ui.small(format!(
                                "clip window: {}..{}",
                                start,
                                trace.clip_end_base_index_exclusive
                                    .map(|value| value.to_string())
                                    .unwrap_or_else(|| "?".to_string())
                            ));
                        }
                        ui.label("Called bases preview");
                        let mut preview = Self::sequencing_trace_called_base_preview(
                            &trace.called_bases,
                            48,
                            4,
                        );
                        ui.add_enabled(
                            false,
                            egui::TextEdit::multiline(&mut preview)
                                .desired_rows(4)
                                .desired_width(f32::INFINITY),
                        );
                        if let Some(comments) = trace.comments_text.as_deref() {
                            ui.collapsing("Trace comments", |ui| {
                                let mut comments_text = comments.to_string();
                                ui.add_enabled(
                                    false,
                                    egui::TextEdit::multiline(&mut comments_text)
                                        .desired_rows(4)
                                        .desired_width(f32::INFINITY),
                                );
                            });
                        }
                    });
                }
            }
            columns[1].separator();
            columns[1].heading("Chromatogram");
            let flank_bp = Self::parse_positive_usize_text(
                &self.sequencing_confirmation_ui.chromatogram_window_bp,
                "chromatogram flank bp",
            )
            .unwrap_or(24);
            let variant_focus_available =
                selected_variant.is_some() && selected_variant_trace.is_some();
            let trace_browser_available = selected_trace.is_some();
            if variant_focus_available || trace_browser_available {
                let mut focus_changed = false;
                columns[1].horizontal_wrapped(|ui| {
                    ui.label("Focus mode");
                    if variant_focus_available {
                        focus_changed |= ui
                            .selectable_value(
                                &mut self.sequencing_confirmation_ui.chromatogram_focus_mode,
                                SequencingChromatogramFocusMode::VariantLocus,
                                SequencingChromatogramFocusMode::VariantLocus.label(),
                            )
                            .changed();
                    }
                    if trace_browser_available {
                        focus_changed |= ui
                            .selectable_value(
                                &mut self.sequencing_confirmation_ui.chromatogram_focus_mode,
                                SequencingChromatogramFocusMode::TraceBaseBrowser,
                                SequencingChromatogramFocusMode::TraceBaseBrowser.label(),
                            )
                            .changed();
                    }
                    ui.label("flank bp");
                    focus_changed |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.sequencing_confirmation_ui.chromatogram_window_bp,
                            )
                            .desired_width(56.0),
                        )
                        .changed();
                });
                if focus_changed {
                    self.save_engine_ops_state();
                }
            }
            let effective_chromatogram_focus = if self
                .sequencing_confirmation_ui
                .chromatogram_focus_mode
                == SequencingChromatogramFocusMode::VariantLocus
                && variant_focus_available
            {
                SequencingChromatogramFocusMode::VariantLocus
            } else {
                SequencingChromatogramFocusMode::TraceBaseBrowser
            };
            if effective_chromatogram_focus == SequencingChromatogramFocusMode::TraceBaseBrowser {
                if let Some(trace) = selected_trace.as_ref() {
                    let browser_base_count = Self::sequencing_trace_browser_base_count(trace);
                    if browser_base_count > 0 {
                        columns[1].horizontal_wrapped(|ui| {
                            let mut changed = false;
                            if ui
                                .small_button("First")
                                .on_hover_text("Jump to the first called base in the selected trace.")
                                .clicked()
                            {
                                self.sequencing_confirmation_ui.chromatogram_trace_base_index = 0;
                                changed = true;
                            }
                            if ui
                                .small_button("Prev")
                                .on_hover_text("Move one called base to the left in the selected trace.")
                                .clicked()
                            {
                                self.sequencing_confirmation_ui.chromatogram_trace_base_index =
                                    self.sequencing_confirmation_ui
                                        .chromatogram_trace_base_index
                                        .saturating_sub(1);
                                changed = true;
                            }
                            if ui
                                .small_button("Next")
                                .on_hover_text("Move one called base to the right in the selected trace.")
                                .clicked()
                            {
                                self.sequencing_confirmation_ui.chromatogram_trace_base_index = (
                                    self.sequencing_confirmation_ui.chromatogram_trace_base_index
                                        + 1
                                )
                                .min(browser_base_count.saturating_sub(1));
                                changed = true;
                            }
                            if ui
                                .small_button("Last")
                                .on_hover_text("Jump to the last called base in the selected trace.")
                                .clicked()
                            {
                                self.sequencing_confirmation_ui.chromatogram_trace_base_index =
                                    browser_base_count.saturating_sub(1);
                                changed = true;
                            }
                            if let Some(start) = trace.clip_start_base_index
                                && ui
                                    .small_button("Clip start")
                                    .on_hover_text("Jump to the first base inside the trace clip window.")
                                    .clicked()
                            {
                                self.sequencing_confirmation_ui.chromatogram_trace_base_index =
                                    Self::sequencing_trace_clamp_base_index(trace, start);
                                changed = true;
                            }
                            if let Some(end) = trace.clip_end_base_index_exclusive
                                && end > 0
                                && ui
                                    .small_button("Clip end")
                                    .on_hover_text("Jump to the last base inside the trace clip window.")
                                    .clicked()
                            {
                                self.sequencing_confirmation_ui.chromatogram_trace_base_index =
                                    Self::sequencing_trace_clamp_base_index(trace, end - 1);
                                changed = true;
                            }
                            if let Some(variant) = selected_variant
                                && variant.trace_id.as_deref().is_some_and(|trace_id| {
                                    trace_id.eq_ignore_ascii_case(trace.trace_id.as_str())
                                })
                                && let Some(variant_base_index) =
                                    Self::sequencing_trace_base_index_for_variant(trace, variant)
                                && ui
                                    .small_button("Jump to variant")
                                    .on_hover_text(
                                        "Center the trace-base browser on the currently selected variant locus.",
                                    )
                                    .clicked()
                            {
                                self.sequencing_confirmation_ui.chromatogram_trace_base_index =
                                    Self::sequencing_trace_clamp_base_index(
                                        trace,
                                        variant_base_index,
                                    );
                                changed = true;
                            }
                            let mut display_index = self
                                .sequencing_confirmation_ui
                                .chromatogram_trace_base_index
                                .saturating_add(1)
                                .min(browser_base_count);
                            changed |= ui
                                .add(
                                    egui::Slider::new(&mut display_index, 1..=browser_base_count)
                                        .text("called base"),
                                )
                                .on_hover_text(
                                    "Browse the chromatogram by called-base index across the selected trace.",
                                )
                                .changed();
                            if changed {
                                self.sequencing_confirmation_ui.chromatogram_trace_base_index =
                                    display_index.saturating_sub(1);
                                self.save_engine_ops_state();
                            }
                        });
                        if let Some(summary) = Self::sequencing_trace_selected_base_summary(
                            trace,
                            self.sequencing_confirmation_ui.chromatogram_trace_base_index,
                        ) {
                            columns[1].small(summary);
                        }
                    }
                }
            } else if !variant_focus_available && trace_browser_available {
                columns[1].small(
                    "Selected variant is not trace-backed, so the chromatogram falls back to the selected-trace base browser.",
                );
            }
            if let Some(report) = selected_report.as_ref() {
                if report.variants.is_empty() {
                    columns[1].small(
                        "No variant-focused loci are recorded in this report yet. Add a baseline sequence or expected-edit targets to get intended-edit/reversion classification.",
                    );
                } else {
                    columns[1].horizontal_wrapped(|ui| {
                        ui.label("Variant focus");
                        let mut selection_changed = false;
                        egui::ComboBox::from_id_salt((
                            "seq_confirm_selected_variant",
                            self.panel_scope_key(),
                        ))
                        .selected_text(
                            report
                                .variants
                                .iter()
                                .find(|row| {
                                    row.variant_id.eq_ignore_ascii_case(
                                        self.sequencing_confirmation_ui.selected_variant_id.trim(),
                                    )
                                })
                                .map(|row| {
                                    format!(
                                        "{} | {} | {}",
                                        row.label,
                                        row.classification.as_str(),
                                        row.trace_id.as_deref().unwrap_or("-")
                                    )
                                })
                                .unwrap_or_else(|| "<select variant>".to_string()),
                        )
                        .show_ui(ui, |ui| {
                            for row in &report.variants {
                                selection_changed |= ui
                                    .selectable_value(
                                        &mut self.sequencing_confirmation_ui.selected_variant_id,
                                        row.variant_id.clone(),
                                        format!(
                                            "{} | {} | {}",
                                            row.label,
                                            row.classification.as_str(),
                                            row.trace_id.as_deref().unwrap_or("-")
                                        ),
                                    )
                                    .changed();
                            }
                        });
                        if selection_changed {
                            let selected_variant_id =
                                self.sequencing_confirmation_ui.selected_variant_id.clone();
                            self.sequencing_confirmation_sync_variant_selection(
                                report,
                                &selected_variant_id,
                            );
                            self.save_engine_ops_state();
                        }
                    });
                    if effective_chromatogram_focus == SequencingChromatogramFocusMode::VariantLocus
                        && let Some(variant) = selected_variant
                    {
                        columns[1].group(|ui| {
                            ui.horizontal_wrapped(|ui| {
                                ui.label(egui::RichText::new(&variant.label).strong());
                                ui.separator();
                                ui.colored_label(
                                    Self::sequencing_confirmation_variant_color(
                                        variant.classification,
                                    ),
                                    variant.classification.as_str(),
                                );
                                ui.separator();
                                ui.small(format!(
                                    "expected='{}' observed='{}'",
                                    variant.expected_bases,
                                    if variant.observed_bases.is_empty() {
                                        "-"
                                    } else {
                                        variant.observed_bases.as_str()
                                    }
                                ));
                                if let Some(baseline) = variant.baseline_bases.as_deref() {
                                    ui.separator();
                                    ui.small(format!("baseline='{}'", baseline));
                                }
                            });
                            ui.small(format!(
                                "evidence={} trace={} confidence={} mean={} peak_center={}",
                                if variant.evidence_id.is_empty() {
                                    "-"
                                } else {
                                    variant.evidence_id.as_str()
                                },
                                variant.trace_id.as_deref().unwrap_or("-"),
                                variant.confidence_count,
                                variant
                                    .confidence_mean
                                    .map(|value| format!("{value:.1}"))
                                    .unwrap_or_else(|| "-".to_string()),
                                variant
                                    .peak_center
                                    .map(|value| value.to_string())
                                    .unwrap_or_else(|| "-".to_string())
                            ));
                            ui.small(&variant.reason);
                            if let Some(trace) = selected_variant_trace.as_ref() {
                                Self::render_sequencing_trace_chromatogram(
                                    ui,
                                    trace,
                                    Some(variant),
                                    None,
                                    flank_bp,
                                );
                            } else if variant.trace_id.is_some() {
                                ui.small(
                                    "Trace-backed variant selected, but the imported trace record could not be loaded.",
                                );
                            } else {
                                ui.small(
                                    "This variant row is not trace-backed, so no chromatogram curve is available.",
                                );
                            }
                        });
                    } else if let Some(trace) = selected_trace.as_ref() {
                        Self::render_sequencing_trace_chromatogram(
                            &mut columns[1],
                            trace,
                            None,
                            Some(self.sequencing_confirmation_ui.chromatogram_trace_base_index),
                            flank_bp,
                        );
                    }
                    columns[1].collapsing("Variant rows", |ui| {
                        let variant_rows = Self::sequencing_confirmation_variant_review_rows(
                            report,
                            self.sequencing_confirmation_ui.review_unresolved_first,
                        );
                        egui::Grid::new(("seq_confirm_variants_grid", self.panel_scope_key()))
                            .num_columns(6)
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Variant");
                                ui.strong("Class");
                                ui.strong("Observed");
                                ui.strong("Baseline");
                                ui.strong("Trace");
                                ui.strong("Reason");
                                ui.end_row();
                                for row in variant_rows {
                                    let selected = row.variant_id.eq_ignore_ascii_case(
                                        self.sequencing_confirmation_ui.selected_variant_id.trim(),
                                    );
                                    if ui
                                        .selectable_label(selected, &row.label)
                                        .on_hover_text(
                                            "Select this variant locus for chromatogram review.",
                                        )
                                        .clicked()
                                    {
                                        self.sequencing_confirmation_sync_variant_selection(
                                            report,
                                            &row.variant_id,
                                        );
                                        self.save_engine_ops_state();
                                    }
                                    ui.colored_label(
                                        Self::sequencing_confirmation_variant_color(
                                            row.classification,
                                        ),
                                        row.classification.as_str(),
                                    );
                                    ui.small(if row.observed_bases.is_empty() {
                                        "-"
                                    } else {
                                        row.observed_bases.as_str()
                                    });
                                    ui.small(row.baseline_bases.as_deref().unwrap_or("-"));
                                    ui.small(row.trace_id.as_deref().unwrap_or("-"));
                                    ui.small(&row.reason);
                                    ui.end_row();
                                }
                            });
                    });
                }
            } else if let Some(trace) = selected_trace.as_ref() {
                Self::render_sequencing_trace_chromatogram(
                    &mut columns[1],
                    trace,
                    None,
                    Some(self.sequencing_confirmation_ui.chromatogram_trace_base_index),
                    flank_bp,
                );
            } else {
                columns[1].small(
                    "Select a persisted report or imported trace to inspect chromatogram curves.",
                );
            }
            columns[1].separator();
            columns[1].heading("Saved Reports + Evidence");
            if report_summaries.is_empty() {
                columns[1].small(
                    "No sequencing-confirmation reports are stored for this sequence yet. Run confirmation to create one.",
                );
            } else {
                columns[1].horizontal_wrapped(|ui| {
                ui.label("Saved report");
                let mut selection_changed = false;
                egui::ComboBox::from_id_salt((
                    "seq_confirm_saved_report",
                    self.panel_scope_key(),
                ))
                .selected_text(
                    report_summaries
                        .iter()
                        .find(|row| {
                            row.report_id.eq_ignore_ascii_case(
                                self.sequencing_confirmation_ui.selected_report_id.trim(),
                            )
                        })
                        .map(|row| {
                            format!(
                                "{} | {} | reads={}",
                                row.report_id,
                                row.overall_status.as_str(),
                                row.read_count
                            )
                        })
                        .unwrap_or_else(|| "<select report>".to_string()),
                )
                .show_ui(ui, |ui| {
                    for row in &report_summaries {
                        selection_changed |= ui
                            .selectable_value(
                                &mut self.sequencing_confirmation_ui.selected_report_id,
                                row.report_id.clone(),
                                format!(
                                    "{} | {} | reads={} targets={}",
                                    row.report_id,
                                    row.overall_status.as_str(),
                                    row.read_count,
                                    row.target_count
                                ),
                            )
                            .changed();
                    }
                });
                if selection_changed {
                    self.save_engine_ops_state();
                }
                if ui
                    .button("Export JSON...")
                    .on_hover_text(
                        "Export the selected sequencing-confirmation report as pretty-printed JSON.",
                    )
                    .clicked()
                {
                    let report_id = self
                        .sequencing_confirmation_ui
                        .selected_report_id
                        .clone();
                    self.export_sequencing_confirmation_report_dialog(&report_id);
                }
                if ui
                    .button("Export TSV...")
                    .on_hover_text(
                        "Export one target-support row per checkpoint from the selected report.",
                    )
                    .clicked()
                {
                    let report_id = self
                        .sequencing_confirmation_ui
                        .selected_report_id
                        .clone();
                    self.export_sequencing_confirmation_support_tsv_dialog(&report_id);
                }
                if ui
                    .button("Copy summary")
                    .on_hover_text(
                        "Copy a compact Markdown summary of unresolved targets, variants, and coverage gaps from the selected report.",
                    )
                    .clicked()
                {
                    let report_id = self
                        .sequencing_confirmation_ui
                        .selected_report_id
                        .clone();
                    self.copy_sequencing_confirmation_unresolved_summary(
                        &report_id,
                        ui.ctx(),
                    );
                }
                if ui
                    .button("Export summary...")
                    .on_hover_text(
                        "Export a compact Markdown summary of unresolved targets, variants, and coverage gaps from the selected report.",
                    )
                    .clicked()
                {
                    let report_id = self
                        .sequencing_confirmation_ui
                        .selected_report_id
                        .clone();
                    self.export_sequencing_confirmation_unresolved_summary_dialog(&report_id);
                }
                });
                if let Some(report) = selected_report.as_ref() {
                    columns[1].group(|ui| {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Overall");
                        ui.colored_label(
                            Self::sequencing_confirmation_status_color(report.overall_status),
                            report.overall_status.as_str(),
                        );
                        ui.separator();
                        ui.small(format!("report_id={}", report.report_id));
                        ui.separator();
                        ui.small(format!(
                            "reads={} traces={} evidence={}",
                            report.read_seq_ids.len(),
                            report.trace_ids.len(),
                            report.reads.len()
                        ));
                        ui.separator();
                        ui.small(format!("targets={}", report.targets.len()));
                        ui.separator();
                        ui.small(format!("variants={}", report.variants.len()));
                        if let Some(baseline) = report.baseline_seq_id.as_deref() {
                            ui.separator();
                            ui.small(format!("baseline={baseline}"));
                        }
                        if !report.warnings.is_empty() {
                            ui.separator();
                            ui.small(format!("warnings={}", report.warnings.len()));
                        }
                    });
                    ui.small(format!(
                        "mode={} rc_allowed={} identity>={:.2} target_coverage>={:.2}",
                        report.alignment_mode.as_str(),
                        report.allow_reverse_complement,
                        report.min_identity_fraction,
                        report.min_target_coverage_fraction
                    ));
                    });
                    if !report.warnings.is_empty() {
                        columns[1].collapsing("Warnings", |ui| {
                        for warning in &report.warnings {
                            ui.small(warning);
                        }
                        });
                    }
                    columns[1].separator();
                    columns[1].label(egui::RichText::new("Construct overview").strong());
                    let review_mode_changed = columns[1]
                        .checkbox(
                            &mut self.sequencing_confirmation_ui.review_unresolved_first,
                            "Prioritize unresolved targets and variants",
                        )
                        .on_hover_text(
                            "Sort target and variant tables so contradicted or insufficient loci stay at the top of the inspection queue.",
                        )
                        .changed();
                    if review_mode_changed {
                        self.save_engine_ops_state();
                    }
                    if let Some(selection) = Self::render_sequencing_confirmation_construct_overview(
                        &mut columns[1],
                        report,
                        report_sequence_length,
                        &self.sequencing_confirmation_ui.selected_target_id,
                        &self.sequencing_confirmation_ui.selected_evidence_id,
                        &self.sequencing_confirmation_ui.selected_variant_id,
                        selected_gap,
                    ) {
                        if self
                            .sequencing_confirmation_apply_overview_selection(report, selection)
                        {
                            self.save_engine_ops_state();
                        }
                    }
                    let unresolved_queue = Self::sequencing_confirmation_unresolved_review_queue(
                        report,
                        report_sequence_length,
                    );
                    let current_unresolved_focus =
                        Self::sequencing_confirmation_current_unresolved_focus(
                            report,
                            &self.sequencing_confirmation_ui.selected_target_id,
                            &self.sequencing_confirmation_ui.selected_variant_id,
                            selected_gap,
                            self.sequencing_confirmation_ui.selected_review_focus_kind,
                        );
                    let current_unresolved_index = current_unresolved_focus
                        .as_ref()
                        .and_then(|focus| unresolved_queue.iter().position(|row| row == focus));
                    columns[1].horizontal_wrapped(|ui| {
                        ui.label("Unresolved review");
                        if unresolved_queue.is_empty() {
                            ui.small("all tracked loci currently look resolved");
                        } else {
                            if let Some(index) = current_unresolved_index {
                                ui.small(format!("{}/{}", index + 1, unresolved_queue.len()));
                                if let Some(focus) = unresolved_queue.get(index) {
                                    ui.separator();
                                    ui.small(Self::sequencing_confirmation_unresolved_focus_label(
                                        report, focus,
                                    ));
                                }
                            } else {
                                ui.small(format!("0/{}", unresolved_queue.len()));
                                ui.separator();
                                ui.small("Choose Next unresolved to start walking the queue.");
                            }
                            let prev_enabled =
                                current_unresolved_index.is_none_or(|index| index > 0);
                            if ui
                                .add_enabled(prev_enabled, egui::Button::new("Prev unresolved"))
                                .on_hover_text(
                                    "Move the saved-report review focus to the previous unresolved target, variant, or coverage gap.",
                                )
                                .clicked()
                                && self.sequencing_confirmation_step_unresolved_focus(
                                    report,
                                    report_sequence_length,
                                    -1,
                                )
                            {
                                self.save_engine_ops_state();
                            }
                            let next_enabled = current_unresolved_index
                                .is_none_or(|index| index + 1 < unresolved_queue.len());
                            if ui
                                .add_enabled(next_enabled, egui::Button::new("Next unresolved"))
                                .on_hover_text(
                                    "Move the saved-report review focus to the next unresolved target, variant, or coverage gap.",
                                )
                                .clicked()
                                && self.sequencing_confirmation_step_unresolved_focus(
                                    report,
                                    report_sequence_length,
                                    1,
                                )
                            {
                                self.save_engine_ops_state();
                            }
                        }
                    });
                    let unresolved_targets = report
                        .targets
                        .iter()
                        .filter(|row| row.status != SequencingConfirmationStatus::Confirmed)
                        .collect::<Vec<_>>();
                    let unresolved_variants = report
                        .variants
                        .iter()
                        .filter(|row| row.status != SequencingConfirmationStatus::Confirmed)
                        .collect::<Vec<_>>();
                    if !unresolved_targets.is_empty() || !unresolved_variants.is_empty() {
                        columns[1].collapsing("Review queue", |ui| {
                            if !unresolved_targets.is_empty() {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label("Targets");
                                    for row in Self::sequencing_confirmation_target_review_rows(
                                        report,
                                        true,
                                    )
                                    .into_iter()
                                    .filter(|row| row.status != SequencingConfirmationStatus::Confirmed)
                                    .take(8)
                                    {
                                        if ui
                                            .small_button(&row.label)
                                            .on_hover_text(&row.reason)
                                            .clicked()
                                        {
                                            self.sequencing_confirmation_sync_target_selection(
                                                report,
                                                &row.target_id,
                                            );
                                            self.save_engine_ops_state();
                                        }
                                    }
                                });
                            }
                            if !unresolved_variants.is_empty() {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label("Variants");
                                    for row in Self::sequencing_confirmation_variant_review_rows(
                                        report,
                                        true,
                                    )
                                    .into_iter()
                                    .filter(|row| row.status != SequencingConfirmationStatus::Confirmed)
                                    .take(8)
                                    {
                                        if ui
                                            .small_button(&row.label)
                                            .on_hover_text(&row.reason)
                                            .clicked()
                                        {
                                            self.sequencing_confirmation_sync_variant_selection(
                                                report,
                                                &row.variant_id,
                                            );
                                            self.save_engine_ops_state();
                                        }
                                    }
                                });
                            }
                        });
                    }
                    if let Some(target) = selected_target {
                        columns[1].group(|ui| {
                            ui.horizontal_wrapped(|ui| {
                                ui.label(egui::RichText::new(&target.label).strong());
                                ui.separator();
                                ui.small(target.kind.as_str());
                                ui.separator();
                                ui.colored_label(
                                    Self::sequencing_confirmation_status_color(target.status),
                                    target.status.as_str(),
                                );
                                ui.separator();
                                ui.small(format!(
                                    "{}..{} ({}/{})",
                                    target.start_0based,
                                    target.end_0based_exclusive,
                                    target.covered_bp,
                                    target.target_length_bp
                                ));
                            });
                            ui.small(format!(
                                "support={} contradiction={} required={}",
                                target.support_read_ids.len(),
                                target.contradicting_read_ids.len(),
                                target.required
                            ));
                            ui.small(&target.reason);
                        });
                    }
                    if let Some((gap_start, gap_end)) = selected_gap {
                        let primer_overlay_report = self
                            .sequencing_confirmation_overlay_report_for_expected_seq(
                                expected_seq_id.as_str(),
                            );
                        let report_matched_primer_overlay_report = self
                            .sequencing_confirmation_report_matched_primer_overlay_report(
                                expected_seq_id.as_str(),
                                report.report_id.as_str(),
                            );
                        let overlapping_targets = report
                            .targets
                            .iter()
                            .filter(|row| {
                                Self::sequencing_confirmation_target_intersects_span(
                                    row, gap_start, gap_end,
                                )
                            })
                            .collect::<Vec<_>>();
                        let (left_flank, right_flank) =
                            Self::sequencing_confirmation_gap_flanking_reads(
                                report, gap_start, gap_end,
                            );
                        let gap_suggestions = primer_overlay_report
                            .map(|overlay| {
                                Self::sequencing_confirmation_gap_primer_suggestions(
                                    overlay, gap_start, gap_end,
                                )
                            })
                            .unwrap_or_default();
                        let gap_proposals = report_matched_primer_overlay_report
                            .map(|overlay| {
                                Self::sequencing_confirmation_gap_primer_proposals(
                                    overlay, gap_start, gap_end,
                                )
                            })
                            .unwrap_or_default();
                        let target_guidance_rows = report_matched_primer_overlay_report
                            .map(|overlay| {
                                overlay
                                    .problem_guidance
                                    .iter()
                                    .filter(|row| {
                                        row.problem_kind == SequencingPrimerProblemKind::Target
                                            && overlapping_targets.iter().any(|target| {
                                                row.problem_id.eq_ignore_ascii_case(
                                                    target.target_id.as_str(),
                                                )
                                            })
                                    })
                                    .collect::<Vec<_>>()
                            })
                            .unwrap_or_default();
                        columns[1].group(|ui| {
                            ui.horizontal_wrapped(|ui| {
                                ui.label(egui::RichText::new("Selected coverage gap").strong());
                                ui.separator();
                                ui.colored_label(
                                    egui::Color32::from_rgb(180, 83, 9),
                                    format!(
                                        "{}..{} ({} bp)",
                                        gap_start,
                                        gap_end,
                                        gap_end.saturating_sub(gap_start)
                                    ),
                                );
                            });
                            if overlapping_targets.is_empty() {
                                ui.small(
                                    "No explicit confirmation target overlaps this unsupported region yet.",
                                );
                            } else {
                                ui.small(format!(
                                    "Overlapping targets: {}",
                                    overlapping_targets
                                        .iter()
                                        .map(|row| row.label.as_str())
                                        .collect::<Vec<_>>()
                                        .join(", ")
                                ));
                            }
                            if let Some(target_id) = Self::sequencing_confirmation_select_target_for_gap(
                                report,
                                gap_start,
                                gap_end,
                            ) {
                                ui.small(format!(
                                    "Closest review focus: {target_id}"
                                ));
                            }
                            ui.small(
                                "This interval currently lacks evidence-span coverage in the saved report.",
                            );
                            match (left_flank, right_flank) {
                                (Some(left), Some(right)) => {
                                    ui.small(format!(
                                        "Nearest covered edges: {} ends at {} ({} bp before gap), {} starts at {} ({} bp after gap).",
                                        left.evidence_id,
                                        left.best_alignment.aligned_target_end_0based_exclusive,
                                        gap_start.saturating_sub(
                                            left.best_alignment.aligned_target_end_0based_exclusive
                                        ),
                                        right.evidence_id,
                                        right.best_alignment.aligned_target_start_0based,
                                        right
                                            .best_alignment
                                            .aligned_target_start_0based
                                            .saturating_sub(gap_end)
                                    ));
                                }
                                (Some(left), None) => {
                                    ui.small(format!(
                                        "Nearest covered edge: {} ends at {} ({} bp before gap); no right-flanking evidence row is retained.",
                                        left.evidence_id,
                                        left.best_alignment.aligned_target_end_0based_exclusive,
                                        gap_start.saturating_sub(
                                            left.best_alignment.aligned_target_end_0based_exclusive
                                        ),
                                    ));
                                }
                                (None, Some(right)) => {
                                    ui.small(format!(
                                        "Nearest covered edge: {} starts at {} ({} bp after gap); no left-flanking evidence row is retained.",
                                        right.evidence_id,
                                        right.best_alignment.aligned_target_start_0based,
                                        right
                                            .best_alignment
                                            .aligned_target_start_0based
                                            .saturating_sub(gap_end),
                                    ));
                                }
                                (None, None) => {
                                    ui.small(
                                        "No retained evidence rows flank this unsupported interval yet.",
                                    );
                                }
                            }
                            if let Some(overlay) = primer_overlay_report {
                                ui.separator();
                                ui.small(format!(
                                    "Primer overlay source: primers={} report={}",
                                    overlay.primer_seq_ids.len(),
                                    overlay.confirmation_report_id.as_deref().unwrap_or("-")
                                ));
                                if !target_guidance_rows.is_empty() {
                                    ui.label(
                                        egui::RichText::new("Target-linked primer guidance")
                                            .strong(),
                                    );
                                    egui::Grid::new((
                                        "seq_confirm_gap_guidance_grid",
                                        self.panel_scope_key(),
                                    ))
                                    .num_columns(5)
                                    .striped(true)
                                    .show(ui, |ui| {
                                        ui.strong("Target");
                                        ui.strong("Primer");
                                        ui.strong("Orientation");
                                        ui.strong("3' dist");
                                        ui.strong("Reason");
                                        ui.end_row();
                                        for row in &target_guidance_rows {
                                            ui.small(&row.problem_label);
                                            ui.small(
                                                row.recommended_primer_seq_id
                                                    .as_deref()
                                                    .unwrap_or("-"),
                                            );
                                            ui.small(
                                                row.recommended_orientation
                                                    .map(|value| value.as_str())
                                                    .unwrap_or("-"),
                                            );
                                            ui.small(
                                                row.recommended_three_prime_distance_bp
                                                    .map(|value| value.to_string())
                                                    .unwrap_or_else(|| "-".to_string()),
                                            );
                                            ui.small(&row.reason);
                                            ui.end_row();
                                        }
                                    });
                                }
                                let helpful_suggestions = gap_suggestions
                                    .iter()
                                    .filter(|row| {
                                        Self::sequencing_confirmation_span_overlap_bp(
                                            row.predicted_read_span_start_0based,
                                            row.predicted_read_span_end_0based_exclusive,
                                            gap_start,
                                            gap_end,
                                        ) > 0
                                    })
                                    .take(4)
                                    .cloned()
                                    .collect::<Vec<_>>();
                                let closest_suggestions = if helpful_suggestions.is_empty() {
                                    gap_suggestions.iter().take(4).cloned().collect::<Vec<_>>()
                                } else {
                                    helpful_suggestions
                                };
                                if !closest_suggestions.is_empty() {
                                    ui.separator();
                                    ui.label(
                                        egui::RichText::new(
                                            if closest_suggestions.iter().any(|row| {
                                                Self::sequencing_confirmation_span_overlap_bp(
                                                    row.predicted_read_span_start_0based,
                                                    row.predicted_read_span_end_0based_exclusive,
                                                    gap_start,
                                                    gap_end,
                                                ) > 0
                                            }) {
                                                "Existing primers that could clarify this gap"
                                            } else {
                                                "Closest existing primer overlays"
                                            },
                                        )
                                        .strong(),
                                    );
                                    egui::Grid::new((
                                        "seq_confirm_gap_suggestions_grid",
                                        self.panel_scope_key(),
                                    ))
                                    .num_columns(6)
                                    .striped(true)
                                    .show(ui, |ui| {
                                        ui.strong("Primer");
                                        ui.strong("Orientation");
                                        ui.strong("Read span");
                                        ui.strong("Gap overlap");
                                        ui.strong("3' dist");
                                        ui.strong("Problem loci");
                                        ui.end_row();
                                        for row in &closest_suggestions {
                                            let overlap_bp =
                                                Self::sequencing_confirmation_span_overlap_bp(
                                                    row.predicted_read_span_start_0based,
                                                    row.predicted_read_span_end_0based_exclusive,
                                                    gap_start,
                                                    gap_end,
                                                );
                                            ui.small(&row.primer_seq_id);
                                            ui.small(row.orientation.as_str());
                                            ui.small(format!(
                                                "{}..{}",
                                                row.predicted_read_span_start_0based,
                                                row.predicted_read_span_end_0based_exclusive
                                            ));
                                            ui.small(format!("{overlap_bp} bp"));
                                            ui.small(
                                                row.three_prime_position_0based
                                                    .abs_diff(Self::sequencing_confirmation_gap_center_0based(gap_start, gap_end))
                                                    .to_string(),
                                            );
                                            ui.small(format!(
                                                "targets={} variants={}",
                                                row.covered_problem_target_ids.len(),
                                                row.covered_problem_variant_ids.len()
                                            ));
                                            ui.end_row();
                                        }
                                    });
                                }
                                let helpful_proposals = gap_proposals
                                    .iter()
                                    .filter(|row| {
                                        Self::sequencing_confirmation_span_overlap_bp(
                                            row.predicted_read_span_start_0based,
                                            row.predicted_read_span_end_0based_exclusive,
                                            gap_start,
                                            gap_end,
                                        ) > 0
                                    })
                                    .take(4)
                                    .cloned()
                                    .collect::<Vec<_>>();
                                if !helpful_proposals.is_empty() {
                                    ui.separator();
                                    ui.label(
                                        egui::RichText::new("Fresh primer proposals for this gap")
                                            .strong(),
                                    );
                                    egui::Grid::new((
                                        "seq_confirm_gap_proposals_grid",
                                        self.panel_scope_key(),
                                    ))
                                    .num_columns(7)
                                    .striped(true)
                                    .show(ui, |ui| {
                                        ui.strong("Problem");
                                        ui.strong("Orientation");
                                        ui.strong("Sequence");
                                        ui.strong("Read span");
                                        ui.strong("Gap overlap");
                                        ui.strong("Tm");
                                        ui.strong("Reason");
                                        ui.end_row();
                                        for row in &helpful_proposals {
                                            let overlap_bp =
                                                Self::sequencing_confirmation_span_overlap_bp(
                                                    row.predicted_read_span_start_0based,
                                                    row.predicted_read_span_end_0based_exclusive,
                                                    gap_start,
                                                    gap_end,
                                                );
                                            ui.small(&row.problem_label);
                                            ui.small(row.orientation.as_str());
                                            ui.small(
                                                egui::RichText::new(&row.primer_sequence)
                                                    .monospace(),
                                            );
                                            ui.small(format!(
                                                "{}..{}",
                                                row.predicted_read_span_start_0based,
                                                row.predicted_read_span_end_0based_exclusive
                                            ));
                                            ui.small(format!("{overlap_bp} bp"));
                                            ui.small(format!("{:.1}", row.tm_c));
                                            ui.small(&row.reason);
                                            ui.end_row();
                                        }
                                    });
                                } else if report_matched_primer_overlay_report.is_some()
                                    && target_guidance_rows.is_empty()
                                    && closest_suggestions.is_empty()
                                {
                                    ui.small(
                                        "The current primer overlay report does not yet offer an existing or proposed primer that reaches this gap.",
                                    );
                                }
                                if primer_overlay_report
                                    .and_then(|overlay| overlay.confirmation_report_id.as_deref())
                                    .is_some_and(|id| !id.eq_ignore_ascii_case(report.report_id.as_str()))
                                {
                                    ui.small(
                                        "The loaded primer overlay report belongs to a different saved confirmation report, so report-linked guidance/proposals are hidden here until you re-run `Suggest primers` for the current report.",
                                    );
                                }
                            } else {
                                ui.small(
                                    "No sequencing-primer overlay report is loaded for this construct. Run `Suggest primers` to see which existing or proposed primers could clarify this gap.",
                                );
                            }
                        });
                    }
                    columns[1].separator();
                    columns[1].label(egui::RichText::new("Targets").strong());
                    egui::ScrollArea::vertical()
                        .id_salt(("seq_confirm_targets_scroll", self.panel_scope_key()))
                        .max_height(220.0)
                        .show(&mut columns[1], |ui| {
                        let target_rows = Self::sequencing_confirmation_target_review_rows(
                            report,
                            self.sequencing_confirmation_ui.review_unresolved_first,
                        );
                        egui::Grid::new(("seq_confirm_targets_grid", self.panel_scope_key()))
                            .num_columns(6)
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Target");
                                ui.strong("Kind");
                                ui.strong("Status");
                                ui.strong("Span");
                                ui.strong("Support");
                                ui.strong("Reason");
                                ui.end_row();
                                for target in target_rows {
                                    let selected = target.target_id.eq_ignore_ascii_case(
                                        self.sequencing_confirmation_ui.selected_target_id.trim(),
                                    );
                                    if ui
                                        .selectable_label(selected, &target.label)
                                        .on_hover_text(
                                            "Select this target as the current construct-review focus.",
                                        )
                                        .clicked()
                                    {
                                        self.sequencing_confirmation_sync_target_selection(
                                            report,
                                            &target.target_id,
                                        );
                                        self.save_engine_ops_state();
                                    }
                                    ui.small(target.kind.as_str());
                                    ui.colored_label(
                                        Self::sequencing_confirmation_status_color(target.status),
                                        target.status.as_str(),
                                    );
                                    ui.small(format!(
                                        "{}..{} ({}/{})",
                                        target.start_0based,
                                        target.end_0based_exclusive,
                                        target.covered_bp,
                                        target.target_length_bp
                                    ));
                                    ui.small(format!(
                                        "+{} / -{}",
                                        target.support_read_ids.len(),
                                        target.contradicting_read_ids.len()
                                    ));
                                    ui.small(&target.reason);
                                    ui.end_row();
                                }
                            });
                        });
                    columns[1].separator();
                    columns[1].label(egui::RichText::new("Reads").strong());
                    if report.reads.is_empty() {
                        columns[1].small(
                            "This report does not contain any read/trace evidence rows yet.",
                        );
                    } else {
                        columns[1].horizontal_wrapped(|ui| {
                            ui.label("Evidence focus");
                            let mut selection_changed = false;
                            egui::ComboBox::from_id_salt((
                                "seq_confirm_selected_evidence",
                                self.panel_scope_key(),
                            ))
                            .selected_text(
                                selected_read
                                    .map(|row| {
                                        format!(
                                            "{} | {} | {}",
                                            row.evidence_id,
                                            row.evidence_kind.as_str(),
                                            row.linked_seq_id
                                                .as_deref()
                                                .unwrap_or(row.read_seq_id.as_str())
                                        )
                                    })
                                    .unwrap_or_else(|| "<select evidence>".to_string()),
                            )
                            .show_ui(ui, |ui| {
                                for row in &report.reads {
                                    selection_changed |= ui
                                        .selectable_value(
                                            &mut self.sequencing_confirmation_ui.selected_evidence_id,
                                            row.evidence_id.clone(),
                                            format!(
                                                "{} | {} | {}",
                                                row.evidence_id,
                                                row.evidence_kind.as_str(),
                                                row.linked_seq_id
                                                    .as_deref()
                                                    .unwrap_or(row.read_seq_id.as_str())
                                            ),
                                        )
                                        .changed();
                                }
                            });
                            let selected_index = selected_read.and_then(|selected| {
                                report.reads.iter().position(|row| {
                                    row.evidence_id.eq_ignore_ascii_case(selected.evidence_id.as_str())
                                })
                            });
                            if let Some(index) = selected_index {
                                if ui
                                    .small_button("Prev")
                                    .on_hover_text("Select the previous evidence row for alignment/discrepancy review.")
                                    .clicked()
                                    && index > 0
                                {
                                    self.sequencing_confirmation_sync_evidence_selection(
                                        report,
                                        &report.reads[index - 1].evidence_id,
                                    );
                                    selection_changed = true;
                                }
                                if ui
                                    .small_button("Next")
                                    .on_hover_text("Select the next evidence row for alignment/discrepancy review.")
                                    .clicked()
                                    && index + 1 < report.reads.len()
                                {
                                    self.sequencing_confirmation_sync_evidence_selection(
                                        report,
                                        &report.reads[index + 1].evidence_id,
                                    );
                                    selection_changed = true;
                                }
                            }
                            if selection_changed {
                                self.save_engine_ops_state();
                            }
                        });
                        egui::ScrollArea::vertical()
                            .id_salt(("seq_confirm_reads_scroll", self.panel_scope_key()))
                            .max_height(280.0)
                            .show(&mut columns[1], |ui| {
                            egui::Grid::new(("seq_confirm_reads_grid", self.panel_scope_key()))
                                .num_columns(10)
                                .striped(true)
                                .show(ui, |ui| {
                                    ui.strong("Evidence");
                                    ui.strong("Kind");
                                    ui.strong("Linked seq");
                                    ui.strong("Orientation");
                                    ui.strong("Usable");
                                    ui.strong("Identity");
                                    ui.strong("Coverage");
                                    ui.strong("Targets");
                                    ui.strong("Discrepancies");
                                    ui.strong("Review");
                                    ui.end_row();
                                    for read in &report.reads {
                                        let is_selected = read.evidence_id.eq_ignore_ascii_case(
                                            self.sequencing_confirmation_ui.selected_evidence_id.trim(),
                                        );
                                        if ui
                                            .selectable_label(is_selected, &read.evidence_id)
                                            .on_hover_text(
                                                "Select this evidence row for the alignment snapshot below.",
                                            )
                                            .clicked()
                                        {
                                            self.sequencing_confirmation_sync_evidence_selection(
                                                report,
                                                &read.evidence_id,
                                            );
                                            self.save_engine_ops_state();
                                        }
                                        ui.small(read.evidence_kind.as_str());
                                        ui.small(
                                            read.linked_seq_id
                                                .as_deref()
                                                .unwrap_or(read.read_seq_id.as_str()),
                                        );
                                        ui.small(match read.orientation {
                                            SequencingReadOrientation::Forward => "forward",
                                            SequencingReadOrientation::ReverseComplement => "reverse_complement",
                                        });
                                        ui.small(if read.usable { "yes" } else { "no" });
                                        ui.small(format!("{:.3}", read.best_alignment.identity_fraction));
                                        ui.small(format!(
                                            "q={:.3} / t={:.3}",
                                            read.best_alignment.query_coverage_fraction,
                                            read.best_alignment.target_coverage_fraction
                                        ));
                                        ui.small(format!(
                                            "+{} / -{}",
                                            read.confirmed_target_ids.len(),
                                            read.contradicted_target_ids.len()
                                        ));
                                        ui.small(format!("{}", read.discrepancies.len()));
                                        if let Some(trace_id) = read.trace_id.as_deref() {
                                            if ui
                                                .small_button("Inspect trace")
                                                .on_hover_text(
                                                    "Load the imported trace record behind this evidence row into the trace review pane.",
                                                )
                                                .clicked()
                                            {
                                                self.sequencing_confirmation_sync_evidence_selection(
                                                    report,
                                                    &read.evidence_id,
                                                );
                                                self.sequencing_confirmation_ui.selected_trace_id =
                                                    trace_id.to_string();
                                                self.save_engine_ops_state();
                                            }
                                        } else {
                                            ui.small("-");
                                        }
                                        ui.end_row();
                                    }
                                });
                            });
                    }
                    if let Some(best_read) = selected_read {
                        columns[1].separator();
                        columns[1].collapsing("Selected evidence alignment snapshot", |ui| {
                        let alignment: &SequenceAlignmentReport = &best_read.best_alignment;
                        ui.small(format!(
                            "evidence={} kind={} linked_seq={}",
                            best_read.evidence_id,
                            best_read.evidence_kind.as_str(),
                            best_read
                                .linked_seq_id
                                .as_deref()
                                .unwrap_or(best_read.read_seq_id.as_str())
                        ));
                        ui.small(format!(
                            "{} -> {} | mode={} score={} cigar={}",
                            alignment.query_seq_id,
                            alignment.target_seq_id,
                            alignment.mode.as_str(),
                            alignment.score,
                            alignment.cigar
                        ));
                        ui.small(format!(
                            "query {}..{} aligned {}..{}",
                            alignment.query_span_start_0based,
                            alignment.query_span_end_0based,
                            alignment.aligned_query_start_0based,
                            alignment.aligned_query_end_0based_exclusive
                        ));
                        ui.small(format!(
                            "target {}..{} aligned {}..{}",
                            alignment.target_span_start_0based,
                            alignment.target_span_end_0based,
                            alignment.aligned_target_start_0based,
                            alignment.aligned_target_end_0based_exclusive
                        ));
                        if let Some(trace_id) = best_read.trace_id.as_deref() {
                            if ui
                                .small_button("Inspect this trace")
                                .on_hover_text(
                                    "Load the imported trace record behind the selected evidence row into the trace review pane.",
                                )
                                .clicked()
                            {
                                self.sequencing_confirmation_ui.selected_trace_id =
                                    trace_id.to_string();
                                self.save_engine_ops_state();
                            }
                        }
                        if !best_read.discrepancies.is_empty() {
                            ui.separator();
                            ui.label(egui::RichText::new("Discrepancies").strong());
                            egui::Grid::new((
                                "seq_confirm_selected_evidence_discrepancies",
                                self.panel_scope_key(),
                            ))
                            .num_columns(5)
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Kind");
                                ui.strong("Query span");
                                ui.strong("Target span");
                                ui.strong("Observed");
                                ui.strong("Expected");
                                ui.end_row();
                                for discrepancy in &best_read.discrepancies {
                                    let discrepancy: &SequencingConfirmationDiscrepancy = discrepancy;
                                    ui.small(discrepancy.kind.as_str());
                                    ui.small(format!(
                                        "{}..{}",
                                        discrepancy.query_start_0based,
                                        discrepancy.query_end_0based_exclusive
                                    ));
                                    ui.small(format!(
                                        "{}..{}",
                                        discrepancy.target_start_0based,
                                        discrepancy.target_end_0based_exclusive
                                    ));
                                    ui.small(if discrepancy.query_bases.is_empty() {
                                        "-"
                                    } else {
                                        discrepancy.query_bases.as_str()
                                    });
                                    ui.small(if discrepancy.target_bases.is_empty() {
                                        "-"
                                    } else {
                                        discrepancy.target_bases.as_str()
                                    });
                                    ui.end_row();
                                }
                            });
                        }
                        });
                    }
                } else {
                    columns[1].small(
                        "Selected report is not loaded. Choose a saved report or run a new confirmation pass.",
                    );
                }
            }
        });
    }
}
