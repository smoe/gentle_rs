//! Primer, qPCR, and PCR Designer UI support for `MainAreaDna`.
//!
//! This submodule owns the primer/qPCR designer state records, report helpers,
//! and specialist rendering entrypoint while keeping the sequence window itself
//! as the orchestration parent.

use super::*;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct PrimerSideConstraintUiState {
    pub(super) min_length: String,
    pub(super) max_length: String,
    pub(super) location_0based: String,
    pub(super) start_0based: String,
    pub(super) end_0based: String,
    pub(super) min_tm_c: String,
    pub(super) max_tm_c: String,
    pub(super) min_gc_fraction: String,
    pub(super) max_gc_fraction: String,
    pub(super) max_anneal_hits: String,
    pub(super) non_annealing_5prime_tail: String,
    pub(super) fixed_5prime: String,
    pub(super) fixed_3prime: String,
    pub(super) required_motifs: String,
    pub(super) forbidden_motifs: String,
    pub(super) locked_positions: String,
}

impl Default for PrimerSideConstraintUiState {
    fn default() -> Self {
        Self {
            min_length: "20".to_string(),
            max_length: "30".to_string(),
            location_0based: String::new(),
            start_0based: String::new(),
            end_0based: String::new(),
            min_tm_c: "55.0".to_string(),
            max_tm_c: "68.0".to_string(),
            min_gc_fraction: "0.35".to_string(),
            max_gc_fraction: "0.70".to_string(),
            max_anneal_hits: "1".to_string(),
            non_annealing_5prime_tail: String::new(),
            fixed_5prime: String::new(),
            fixed_3prime: String::new(),
            required_motifs: String::new(),
            forbidden_motifs: String::new(),
            locked_positions: String::new(),
        }
    }
}

impl PrimerSideConstraintUiState {
    pub(super) fn qpcr_primer_default() -> Self {
        Self {
            min_length: "18".to_string(),
            max_length: "24".to_string(),
            min_tm_c: "55.0".to_string(),
            max_tm_c: "65.0".to_string(),
            max_gc_fraction: "0.75".to_string(),
            ..Self::default()
        }
    }

    pub(super) fn qpcr_probe_default() -> Self {
        Self {
            min_length: "20".to_string(),
            max_length: "30".to_string(),
            min_tm_c: "63.0".to_string(),
            max_tm_c: "72.0".to_string(),
            max_gc_fraction: "0.80".to_string(),
            ..Self::default()
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct PrimerPairConstraintUiState {
    pub(super) require_roi_flanking: bool,
    pub(super) required_amplicon_motifs: String,
    pub(super) forbidden_amplicon_motifs: String,
    pub(super) fixed_amplicon_start_0based: String,
    pub(super) fixed_amplicon_end_0based_exclusive: String,
}

impl Default for PrimerPairConstraintUiState {
    fn default() -> Self {
        Self {
            require_roi_flanking: false,
            required_amplicon_motifs: String::new(),
            forbidden_amplicon_motifs: String::new(),
            fixed_amplicon_start_0based: String::new(),
            fixed_amplicon_end_0based_exclusive: String::new(),
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub(super) struct PcrQueuedRegionUiState {
    pub(super) template: String,
    pub(super) source_label: String,
    pub(super) start_0based: usize,
    pub(super) end_0based_exclusive: usize,
}

impl Default for PcrQueuedRegionUiState {
    fn default() -> Self {
        Self {
            template: String::new(),
            source_label: String::new(),
            start_0based: 0,
            end_0based_exclusive: 0,
        }
    }
}

impl PcrQueuedRegionUiState {
    pub(super) fn span_len_bp(&self) -> usize {
        self.end_0based_exclusive.saturating_sub(self.start_0based)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub(super) struct PcrBatchResultRowUiState {
    pub(super) region_index_1based: usize,
    pub(super) template: String,
    pub(super) source_label: String,
    pub(super) start_0based: usize,
    pub(super) end_0based_exclusive: usize,
    pub(super) report_id: String,
    pub(super) report_pair_count: Option<usize>,
    pub(super) report_backend_requested: Option<String>,
    pub(super) report_backend_used: Option<String>,
    pub(super) report_primer3_explain: Option<String>,
    pub(super) report_note: Option<String>,
    pub(super) report_error: Option<String>,
    pub(super) copy_requested: bool,
    pub(super) copy_seq_id: Option<String>,
    pub(super) copy_error: Option<String>,
}

impl Default for PcrBatchResultRowUiState {
    fn default() -> Self {
        Self {
            region_index_1based: 0,
            template: String::new(),
            source_label: String::new(),
            start_0based: 0,
            end_0based_exclusive: 0,
            report_id: String::new(),
            report_pair_count: None,
            report_backend_requested: None,
            report_backend_used: None,
            report_primer3_explain: None,
            report_note: None,
            report_error: None,
            copy_requested: false,
            copy_seq_id: None,
            copy_error: None,
        }
    }
}

impl PcrBatchResultRowUiState {
    pub(super) fn span_len_bp(&self) -> usize {
        self.end_0based_exclusive.saturating_sub(self.start_0based)
    }

    pub(super) fn report_status_label(&self) -> &'static str {
        if self.report_error.is_some() {
            "report:failed"
        } else if matches!(self.report_pair_count, Some(0)) {
            "report:empty"
        } else {
            "report:ok"
        }
    }

    pub(super) fn copy_status_label(&self) -> &'static str {
        if !self.copy_requested {
            "copy:off"
        } else if self.copy_error.is_some() {
            "copy:failed"
        } else if self.copy_seq_id.is_some() {
            "copy:ok"
        } else {
            "copy:none"
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub(super) enum PcrPaintRole {
    #[default]
    Roi,
    UpstreamPrimerWindow,
    DownstreamPrimerWindow,
}

impl PcrPaintRole {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Roi => "ROI",
            Self::UpstreamPrimerWindow => "Upstream primer window",
            Self::DownstreamPrimerWindow => "Downstream primer window",
        }
    }

    pub(super) fn short_label(self) -> &'static str {
        match self {
            Self::Roi => "ROI",
            Self::UpstreamPrimerWindow => "Upstream",
            Self::DownstreamPrimerWindow => "Downstream",
        }
    }

    pub(super) fn color(self) -> egui::Color32 {
        match self {
            Self::Roi => egui::Color32::from_rgb(68, 153, 88),
            Self::UpstreamPrimerWindow => egui::Color32::from_rgb(198, 72, 78),
            Self::DownstreamPrimerWindow => egui::Color32::from_rgb(72, 124, 186),
        }
    }

    pub(super) fn source_label(self) -> &'static str {
        match self {
            Self::Roi => "painted ROI",
            Self::UpstreamPrimerWindow => "painted upstream window",
            Self::DownstreamPrimerWindow => "painted downstream window",
        }
    }

    pub(super) fn hover_text(self) -> &'static str {
        match self {
            Self::Roi => "Green role. Drag on linear map to paint one PCR ROI interval.",
            Self::UpstreamPrimerWindow => {
                "Red role. Drag on linear map to paint one upstream primer-window interval."
            }
            Self::DownstreamPrimerWindow => {
                "Blue role. Drag on linear map to paint one downstream primer-window interval."
            }
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
#[serde(default)]
pub(super) struct PcrPaintIntervalsUiState {
    pub(super) roi: Option<(usize, usize)>,
    pub(super) upstream_window: Option<(usize, usize)>,
    pub(super) downstream_window: Option<(usize, usize)>,
}

impl PcrPaintIntervalsUiState {
    pub(super) fn get(&self, role: PcrPaintRole) -> Option<(usize, usize)> {
        match role {
            PcrPaintRole::Roi => self.roi,
            PcrPaintRole::UpstreamPrimerWindow => self.upstream_window,
            PcrPaintRole::DownstreamPrimerWindow => self.downstream_window,
        }
    }

    pub(super) fn set(&mut self, role: PcrPaintRole, interval: Option<(usize, usize)>) {
        match role {
            PcrPaintRole::Roi => self.roi = interval,
            PcrPaintRole::UpstreamPrimerWindow => self.upstream_window = interval,
            PcrPaintRole::DownstreamPrimerWindow => self.downstream_window = interval,
        }
    }

    pub(super) fn clear_all(&mut self) {
        self.roi = None;
        self.upstream_window = None;
        self.downstream_window = None;
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct PrimerDesignOpsUiState {
    pub(super) roi_start_0based: String,
    pub(super) roi_end_0based: String,
    pub(super) forward: PrimerSideConstraintUiState,
    pub(super) reverse: PrimerSideConstraintUiState,
    pub(super) pair_constraints: PrimerPairConstraintUiState,
    pub(super) min_amplicon_bp: String,
    pub(super) max_amplicon_bp: String,
    pub(super) max_tm_delta_c: String,
    pub(super) max_pairs: String,
    pub(super) report_id: String,
    pub(super) specificity_target_genome_id: String,
    pub(super) specificity_pair_rank_1based: String,
    pub(super) specificity_max_target_amplicon_bp: String,
    pub(super) specificity_max_hits_per_primer: String,
    pub(super) restriction_cloning: RestrictionCloningPcrHandoffUiState,
}

impl Default for PrimerDesignOpsUiState {
    fn default() -> Self {
        Self {
            roi_start_0based: "0".to_string(),
            roi_end_0based: "0".to_string(),
            forward: PrimerSideConstraintUiState::default(),
            reverse: PrimerSideConstraintUiState::default(),
            pair_constraints: PrimerPairConstraintUiState::default(),
            min_amplicon_bp: "120".to_string(),
            max_amplicon_bp: "1200".to_string(),
            max_tm_delta_c: "2.0".to_string(),
            max_pairs: "200".to_string(),
            report_id: "primer_report_gui".to_string(),
            specificity_target_genome_id: String::new(),
            specificity_pair_rank_1based: "1".to_string(),
            specificity_max_target_amplicon_bp: "4000".to_string(),
            specificity_max_hits_per_primer: "500".to_string(),
            restriction_cloning: RestrictionCloningPcrHandoffUiState::default(),
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct RestrictionCloningPcrHandoffUiState {
    pub(super) selected_pair_rank_1based: String,
    pub(super) destination_vector_seq_id: String,
    pub(super) mode: RestrictionCloningPcrHandoffMode,
    pub(super) forward_enzyme: String,
    pub(super) reverse_enzyme: String,
    pub(super) forward_leader_5prime: String,
    pub(super) reverse_leader_5prime: String,
    pub(super) selected_saved_report_id: String,
}

impl Default for RestrictionCloningPcrHandoffUiState {
    fn default() -> Self {
        Self {
            selected_pair_rank_1based: "1".to_string(),
            destination_vector_seq_id: String::new(),
            mode: RestrictionCloningPcrHandoffMode::SingleSite,
            forward_enzyme: String::new(),
            reverse_enzyme: String::new(),
            forward_leader_5prime: String::new(),
            reverse_leader_5prime: String::new(),
            selected_saved_report_id: String::new(),
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct QpcrDesignOpsUiState {
    pub(super) transcript_targeting: QpcrTranscriptTargetingUiState,
    pub(super) roi_start_0based: String,
    pub(super) roi_end_0based: String,
    pub(super) forward: PrimerSideConstraintUiState,
    pub(super) reverse: PrimerSideConstraintUiState,
    pub(super) probe: PrimerSideConstraintUiState,
    pub(super) pair_constraints: PrimerPairConstraintUiState,
    pub(super) min_amplicon_bp: String,
    pub(super) max_amplicon_bp: String,
    pub(super) max_tm_delta_c: String,
    pub(super) max_probe_tm_delta_c: String,
    pub(super) max_assays: String,
    pub(super) report_id: String,
    pub(super) selected_assay_rank_1based: String,
}

impl Default for QpcrDesignOpsUiState {
    fn default() -> Self {
        Self {
            transcript_targeting: QpcrTranscriptTargetingUiState::default(),
            roi_start_0based: "0".to_string(),
            roi_end_0based: "0".to_string(),
            forward: PrimerSideConstraintUiState::qpcr_primer_default(),
            reverse: PrimerSideConstraintUiState::qpcr_primer_default(),
            probe: PrimerSideConstraintUiState::qpcr_probe_default(),
            pair_constraints: PrimerPairConstraintUiState::default(),
            min_amplicon_bp: "80".to_string(),
            max_amplicon_bp: "200".to_string(),
            max_tm_delta_c: "3.0".to_string(),
            max_probe_tm_delta_c: "10.0".to_string(),
            max_assays: "200".to_string(),
            report_id: "qpcr_report_gui".to_string(),
            selected_assay_rank_1based: "1".to_string(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub(super) enum QpcrTranscriptIntentUiMode {
    #[default]
    Genomic,
    SharedAcrossTranscripts,
    SpecificTranscript,
}

impl QpcrTranscriptIntentUiMode {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Genomic => "Genomic",
            Self::SharedAcrossTranscripts => "Shared across transcripts",
            Self::SpecificTranscript => "Specific transcript",
        }
    }

    pub(super) fn summary_label(self) -> &'static str {
        match self {
            Self::Genomic => "genomic / no transcript constraint",
            Self::SharedAcrossTranscripts => "shared exon or exon-chain across transcripts",
            Self::SpecificTranscript => "specific transcript",
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct QpcrTranscriptTargetingUiState {
    pub(super) assay_intent: QpcrTranscriptIntentUiMode,
    pub(super) specificity_evidence: QpcrTranscriptSpecificityEvidence,
    pub(super) source_seq_id: String,
    pub(super) source_feature_id: Option<usize>,
    pub(super) group_label: String,
    pub(super) transcript_count: usize,
    pub(super) strand: String,
    pub(super) selected_transcript_feature_id: Option<usize>,
    pub(super) selected_transcript_id: String,
    pub(super) selected_transcript_label: String,
}

impl Default for QpcrTranscriptTargetingUiState {
    fn default() -> Self {
        Self {
            assay_intent: QpcrTranscriptIntentUiMode::Genomic,
            specificity_evidence: QpcrTranscriptSpecificityEvidence::EitherPreferJunction,
            source_seq_id: String::new(),
            source_feature_id: None,
            group_label: String::new(),
            transcript_count: 0,
            strand: String::new(),
            selected_transcript_feature_id: None,
            selected_transcript_id: String::new(),
            selected_transcript_label: String::new(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct QpcrCartoonPreviewGeometry {
    pub(super) source_label: String,
    pub(super) roi_left_bp: usize,
    pub(super) probe_window_left_margin_bp: usize,
    pub(super) probe_site_bp: usize,
    pub(super) probe_window_right_margin_bp: usize,
    pub(super) roi_right_bp: usize,
    pub(super) forward_window_margin_bp: usize,
    pub(super) forward_primer_site_bp: usize,
    pub(super) reverse_window_margin_bp: usize,
    pub(super) reverse_primer_site_bp: usize,
}

impl QpcrCartoonPreviewGeometry {
    pub(super) fn total_roi_bp(&self) -> usize {
        self.roi_left_bp
            .saturating_add(self.probe_window_left_margin_bp)
            .saturating_add(self.probe_site_bp)
            .saturating_add(self.probe_window_right_margin_bp)
            .saturating_add(self.roi_right_bp)
    }

    pub(super) fn bindings_feature_override_count(&self) -> usize {
        pcr_assay_qpcr_geometry_bindings(
            self.roi_left_bp,
            self.probe_window_left_margin_bp,
            self.probe_site_bp,
            self.probe_window_right_margin_bp,
            self.roi_right_bp,
            self.forward_window_margin_bp,
            self.forward_primer_site_bp,
            self.reverse_window_margin_bp,
            self.reverse_primer_site_bp,
        )
        .feature_overrides
        .len()
    }

    pub(super) fn short_summary(&self) -> String {
        format!(
            "source={} | roi={} bp | F primer={} bp (+{} bp window margin) | probe={} bp (+{}/{} bp window margins) | R primer={} bp (+{} bp window margin)",
            self.source_label,
            self.total_roi_bp(),
            self.forward_primer_site_bp,
            self.forward_window_margin_bp,
            self.probe_site_bp,
            self.probe_window_left_margin_bp,
            self.probe_window_right_margin_bp,
            self.reverse_primer_site_bp,
            self.reverse_window_margin_bp
        )
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct QpcrSplicingContextSummary {
    pub(super) assay_class_label: String,
    pub(super) explanation: String,
    pub(super) covered_junction_labels: Vec<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub(super) enum PcrDesignerMode {
    #[default]
    PrimerPairs,
    QpcrAssays,
}

impl PcrDesignerMode {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::PrimerPairs => "Pair PCR",
            Self::QpcrAssays => "qPCR",
        }
    }

    pub(super) fn helper_text(self) -> &'static str {
        match self {
            Self::PrimerPairs => {
                "Pair-PCR mode focuses on ROI-flanking primer-pair search, simple-PCR starter flow, and queued multi-region batch runs."
            }
            Self::QpcrAssays => {
                "qPCR mode reuses the same template and ROI context, then adds probe-side constraints and assay scoring on top of retained primer pairs."
            }
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct PrimerDesignBatchSpec {
    pub(super) forward: PrimerDesignSideConstraint,
    pub(super) reverse: PrimerDesignSideConstraint,
    pub(super) pair_constraints: PrimerDesignPairConstraint,
    pub(super) min_amplicon_bp: usize,
    pub(super) max_amplicon_bp: usize,
    pub(super) max_tm_delta_c: Option<f64>,
    pub(super) max_pairs: Option<usize>,
}

#[derive(Clone, Debug)]
pub(super) struct PreparedPrimerDesignBatchInputs {
    pub(super) spec: PrimerDesignBatchSpec,
    pub(super) report_base: String,
    pub(super) queued_regions: Vec<PcrQueuedRegionUiState>,
    pub(super) create_copies: bool,
}

impl MainAreaDna {
    pub(super) fn parse_required_f64_text(raw: &str, field_name: &str) -> Result<f64, String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(format!("Invalid {field_name}: expected a number"));
        }
        trimmed
            .parse::<f64>()
            .map_err(|_| format!("Invalid {field_name}: expected a number"))
    }

    pub(super) fn parse_required_i32_text(raw: &str, field_name: &str) -> Result<i32, String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(format!("Invalid {field_name}: expected an integer"));
        }
        trimmed
            .parse::<i32>()
            .map_err(|_| format!("Invalid {field_name}: expected an integer"))
    }

    pub(super) fn parse_fraction_0_to_1(raw: &str, field_name: &str) -> Result<f64, String> {
        let value = Self::parse_required_f64_text(raw, field_name)?;
        if !(0.0..=1.0).contains(&value) {
            return Err(format!("Invalid {field_name}: expected within 0.0..=1.0"));
        }
        Ok(value)
    }

    pub(super) fn parse_primer_locked_positions(
        raw: &str,
        field_name: &str,
    ) -> Result<Vec<PrimerDesignBaseLock>, String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Ok(vec![]);
        }
        let mut seen_offsets: HashSet<usize> = HashSet::new();
        let mut locks: Vec<PrimerDesignBaseLock> = vec![];
        for token in trimmed.split(',') {
            let entry = token.trim();
            if entry.is_empty() {
                continue;
            }
            let (offset_raw, base_raw) = entry.split_once(':').ok_or_else(|| {
                format!("Invalid {field_name} entry '{entry}' (expected offset:base)")
            })?;
            let offset = offset_raw.trim().parse::<usize>().map_err(|_| {
                format!(
                    "Invalid {field_name} offset '{}' (expected integer)",
                    offset_raw.trim()
                )
            })?;
            if !seen_offsets.insert(offset) {
                return Err(format!(
                    "Invalid {field_name}: duplicate offset {} is not allowed",
                    offset
                ));
            }
            let base = base_raw.trim().to_string();
            if base.is_empty() {
                return Err(format!(
                    "Invalid {field_name} entry '{entry}': base must not be empty"
                ));
            }
            locks.push(PrimerDesignBaseLock {
                offset_0based: offset,
                base,
            });
        }
        Ok(locks)
    }

    pub(super) fn parse_primer_side_constraint_ui(
        side: &PrimerSideConstraintUiState,
        label: &str,
    ) -> Result<PrimerDesignSideConstraint, String> {
        let min_length =
            Self::parse_positive_usize_text(&side.min_length, &format!("{label}.min_length"))?;
        let max_length =
            Self::parse_positive_usize_text(&side.max_length, &format!("{label}.max_length"))?;
        if min_length > max_length {
            return Err(format!(
                "Invalid {label}: min_length ({min_length}) must be <= max_length ({max_length})"
            ));
        }
        let max_anneal_hits = Self::parse_positive_usize_text(
            &side.max_anneal_hits,
            &format!("{label}.max_anneal_hits"),
        )?;
        let min_tm_c = Self::parse_required_f64_text(&side.min_tm_c, &format!("{label}.min_tm_c"))?;
        let max_tm_c = Self::parse_required_f64_text(&side.max_tm_c, &format!("{label}.max_tm_c"))?;
        let min_gc_fraction = Self::parse_required_f64_text(
            &side.min_gc_fraction,
            &format!("{label}.min_gc_fraction"),
        )?;
        let max_gc_fraction = Self::parse_required_f64_text(
            &side.max_gc_fraction,
            &format!("{label}.max_gc_fraction"),
        )?;
        let non_annealing_5prime_tail = if side.non_annealing_5prime_tail.trim().is_empty() {
            None
        } else {
            Some(side.non_annealing_5prime_tail.trim().to_string())
        };
        let fixed_5prime = if side.fixed_5prime.trim().is_empty() {
            None
        } else {
            Some(side.fixed_5prime.trim().to_string())
        };
        let fixed_3prime = if side.fixed_3prime.trim().is_empty() {
            None
        } else {
            Some(side.fixed_3prime.trim().to_string())
        };
        Ok(PrimerDesignSideConstraint {
            min_length,
            max_length,
            location_0based: Self::parse_optional_usize_text(
                &side.location_0based,
                &format!("{label}.location_0based"),
            )?,
            start_0based: Self::parse_optional_usize_text(
                &side.start_0based,
                &format!("{label}.start_0based"),
            )?,
            end_0based: Self::parse_optional_usize_text(
                &side.end_0based,
                &format!("{label}.end_0based"),
            )?,
            min_tm_c,
            max_tm_c,
            min_gc_fraction,
            max_gc_fraction,
            max_anneal_hits,
            non_annealing_5prime_tail,
            fixed_5prime,
            fixed_3prime,
            required_motifs: Self::parse_ids(&side.required_motifs),
            forbidden_motifs: Self::parse_ids(&side.forbidden_motifs),
            locked_positions: Self::parse_primer_locked_positions(
                &side.locked_positions,
                &format!("{label}.locked_positions"),
            )?,
        })
    }

    pub(super) fn parse_primer_pair_constraint_ui(
        pair: &PrimerPairConstraintUiState,
    ) -> Result<PrimerDesignPairConstraint, String> {
        Ok(PrimerDesignPairConstraint {
            require_roi_flanking: pair.require_roi_flanking,
            required_amplicon_motifs: Self::parse_ids(&pair.required_amplicon_motifs),
            forbidden_amplicon_motifs: Self::parse_ids(&pair.forbidden_amplicon_motifs),
            fixed_amplicon_start_0based: Self::parse_optional_usize_text(
                &pair.fixed_amplicon_start_0based,
                "pair_constraints.fixed_amplicon_start_0based",
            )?,
            fixed_amplicon_end_0based_exclusive: Self::parse_optional_usize_text(
                &pair.fixed_amplicon_end_0based_exclusive,
                "pair_constraints.fixed_amplicon_end_0based_exclusive",
            )?,
        })
    }

    pub(super) fn sanitize_id_component(raw: &str, fallback: &str) -> String {
        let mut out = String::new();
        for ch in raw.trim().chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                out.push(ch);
            } else {
                out.push('_');
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            fallback.to_string()
        } else {
            trimmed.to_string()
        }
    }

    pub(super) fn build_primer_design_batch_spec(&self) -> Result<PrimerDesignBatchSpec, String> {
        let ui = &self.primer_design_ui;
        let max_pairs = Self::parse_optional_usize_text(&ui.max_pairs, "max_pairs")?;
        if matches!(max_pairs, Some(0)) {
            return Err("Invalid max_pairs: expected >= 1 when provided".to_string());
        }
        let min_amplicon_bp =
            Self::parse_positive_usize_text(&ui.min_amplicon_bp, "min_amplicon_bp")?;
        let max_amplicon_bp =
            Self::parse_positive_usize_text(&ui.max_amplicon_bp, "max_amplicon_bp")?;
        if min_amplicon_bp > max_amplicon_bp {
            return Err(format!(
                "Invalid amplicon window: min_amplicon_bp ({min_amplicon_bp}) must be <= max_amplicon_bp ({max_amplicon_bp})"
            ));
        }
        Ok(PrimerDesignBatchSpec {
            forward: Self::parse_primer_side_constraint_ui(&ui.forward, "forward")?,
            reverse: Self::parse_primer_side_constraint_ui(&ui.reverse, "reverse")?,
            pair_constraints: Self::parse_primer_pair_constraint_ui(&ui.pair_constraints)?,
            min_amplicon_bp,
            max_amplicon_bp,
            max_tm_delta_c: Self::parse_optional_f64_text(&ui.max_tm_delta_c, "max_tm_delta_c")?,
            max_pairs,
        })
    }

    pub(super) fn validate_pair_pcr_min_amplicon_against_roi(
        min_amplicon_bp: usize,
        roi_start_0based: usize,
        roi_end_0based: usize,
        context: &str,
    ) -> Result<(), String> {
        let roi_len_bp = roi_end_0based.saturating_sub(roi_start_0based);
        if min_amplicon_bp > roi_len_bp {
            return Err(format!(
                "Invalid amplicon window for {context}: min_amplicon_bp ({min_amplicon_bp}) must be <= ROI length ({roi_len_bp})"
            ));
        }
        Ok(())
    }

    pub(super) fn prepare_primer_pair_design_batch_inputs(
        &self,
    ) -> Result<PreparedPrimerDesignBatchInputs, String> {
        if self.pcr_queued_regions_ui.is_empty() {
            return Err("PCR region queue is empty; add at least one region first".to_string());
        }
        let spec = self.build_primer_design_batch_spec()?;
        if let Some((idx, row)) = self
            .pcr_queued_regions_ui
            .iter()
            .enumerate()
            .find(|(_, row)| spec.min_amplicon_bp > row.span_len_bp())
        {
            return Err(format!(
                "Queued PCR region #{:02} '{}' is {} bp long, so min_amplicon_bp ({}) must be <= ROI length before running batch primer design",
                idx + 1,
                row.source_label,
                row.span_len_bp(),
                spec.min_amplicon_bp
            ));
        }
        let fallback_template = self
            .pcr_queued_regions_ui
            .first()
            .map(|row| row.template.clone())
            .unwrap_or_else(|| "template".to_string());
        let report_base = if self.primer_design_ui.report_id.trim().is_empty() {
            Self::sanitize_id_component(
                &format!("{fallback_template}_primer_report"),
                "primer_report_gui",
            )
        } else {
            Self::sanitize_id_component(&self.primer_design_ui.report_id, "primer_report_gui")
        };
        Ok(PreparedPrimerDesignBatchInputs {
            spec,
            report_base,
            queued_regions: self.pcr_queued_regions_ui.clone(),
            create_copies: self.pcr_batch_create_extract_copies,
        })
    }

    pub(super) fn execute_primer_pair_design_batch(
        engine: &mut GentleEngine,
        queued_regions: &[PcrQueuedRegionUiState],
        spec: &PrimerDesignBatchSpec,
        report_base: &str,
        create_copies: bool,
        mut on_progress: impl FnMut(PrimerDesignBatchProgress),
    ) -> PrimerDesignBatchOutcome {
        let mut batch_rows: Vec<PcrBatchResultRowUiState> =
            Vec::with_capacity(queued_regions.len());
        let mut successful_report_ids: Vec<String> = vec![];
        let mut failed_report_rows: Vec<String> = vec![];
        let mut created_copy_ids: Vec<String> = vec![];
        let mut failed_copy_rows: Vec<String> = vec![];
        let mut created_seq_ids: Vec<String> = vec![];

        for (index, region) in queued_regions.iter().enumerate() {
            on_progress(PrimerDesignBatchProgress {
                region_index_1based: index + 1,
                region_count: queued_regions.len(),
                template: region.template.clone(),
                source_label: region.source_label.clone(),
                start_0based: region.start_0based,
                end_0based_exclusive: region.end_0based_exclusive,
            });
            let report_id = format!("{report_base}_r{:02}", index + 1);
            let region_label = format!(
                "r{:02} {} [{}..{} {}]",
                index + 1,
                region.template,
                region.start_0based,
                region.end_0based_exclusive,
                region.source_label
            );
            let mut row = PcrBatchResultRowUiState {
                region_index_1based: index + 1,
                template: region.template.clone(),
                source_label: region.source_label.clone(),
                start_0based: region.start_0based,
                end_0based_exclusive: region.end_0based_exclusive,
                report_id: report_id.clone(),
                report_pair_count: None,
                report_backend_requested: None,
                report_backend_used: None,
                report_primer3_explain: None,
                report_note: None,
                report_error: None,
                copy_requested: create_copies,
                copy_seq_id: None,
                copy_error: None,
            };
            if create_copies {
                let template_component = Self::sanitize_id_component(&region.template, "template");
                let copy_output_id = format!("{}_pcr_roi_{}", template_component, index + 1);
                match engine.apply(Operation::ExtractRegion {
                    input: region.template.clone(),
                    from: region.start_0based,
                    to: region.end_0based_exclusive,
                    output_id: Some(copy_output_id),
                }) {
                    Ok(result) => {
                        if result.created_seq_ids.is_empty() {
                            let message = format!(
                                "ExtractRegion for {region_label} did not create a sequence"
                            );
                            row.copy_error = Some(message.clone());
                            failed_copy_rows.push(message);
                        } else {
                            if let Some(copy_seq_id) = result.created_seq_ids.first() {
                                row.copy_seq_id = Some(copy_seq_id.clone());
                            }
                            created_copy_ids.extend(result.created_seq_ids.clone());
                            created_seq_ids.extend(result.created_seq_ids);
                        }
                    }
                    Err(err) => {
                        let message = format!("{region_label}: {}", err.message);
                        row.copy_error = Some(message.clone());
                        failed_copy_rows.push(message);
                    }
                }
            }
            match engine.apply(Operation::DesignPrimerPairs {
                template: region.template.clone(),
                roi_start_0based: region.start_0based,
                roi_end_0based: region.end_0based_exclusive,
                forward: spec.forward.clone(),
                reverse: spec.reverse.clone(),
                pair_constraints: spec.pair_constraints.clone(),
                min_amplicon_bp: spec.min_amplicon_bp,
                max_amplicon_bp: spec.max_amplicon_bp,
                max_tm_delta_c: spec.max_tm_delta_c,
                max_pairs: spec.max_pairs,
                report_id: Some(report_id.clone()),
            }) {
                Ok(_) => {
                    if let Ok(report) = engine.get_primer_design_report(&report_id) {
                        row.report_pair_count = Some(report.pair_count);
                        row.report_backend_requested = Some(report.backend.requested.clone());
                        row.report_backend_used = Some(report.backend.used.clone());
                        row.report_primer3_explain = report.backend.primer3_explain.clone();
                        if report.pair_count == 0 {
                            let mut detail = format!(
                                "no accepted primer pairs (backend {}->{})",
                                report.backend.requested, report.backend.used
                            );
                            if let Some(explain) = report
                                .backend
                                .primer3_explain
                                .as_deref()
                                .map(str::trim)
                                .filter(|raw| !raw.is_empty())
                            {
                                detail.push_str(&format!("; {explain}"));
                            }
                            row.report_note = Some(detail);
                        }
                    }
                    successful_report_ids.push(report_id);
                }
                Err(err) => {
                    let message = format!("{region_label}: {}", err.message);
                    row.report_error = Some(message.clone());
                    failed_report_rows.push(message);
                }
            }
            batch_rows.push(row);
        }

        PrimerDesignBatchOutcome {
            queued_region_count: queued_regions.len(),
            create_copies,
            batch_rows,
            successful_report_ids,
            failed_report_rows,
            created_copy_ids,
            failed_copy_rows,
            created_seq_ids,
        }
    }

    pub(super) fn apply_primer_pair_design_batch_outcome(
        &mut self,
        outcome: PrimerDesignBatchOutcome,
        elapsed_ms: u128,
    ) {
        self.pcr_batch_results_ui = outcome.batch_rows;
        self.save_engine_ops_state();

        if !outcome.created_seq_ids.is_empty() {
            self.last_created_seq_ids = outcome.created_seq_ids.clone();
            self.export_pool_inputs_text = outcome.created_seq_ids.join(", ");
        }

        let mut status_lines = vec![format!(
            "PCR primer batch finished in {} ms: {} queued region(s), {} succeeded, {} failed",
            elapsed_ms,
            outcome.queued_region_count,
            outcome.successful_report_ids.len(),
            outcome.failed_report_rows.len()
        )];
        if !outcome.successful_report_ids.is_empty() {
            status_lines.push(format!(
                "Report IDs: {}",
                outcome.successful_report_ids.join(", ")
            ));
        }
        let zero_pair_reports = self
            .pcr_batch_results_ui
            .iter()
            .filter(|row| matches!(row.report_pair_count, Some(0)))
            .count();
        if zero_pair_reports > 0 {
            status_lines.push(format!(
                "Reports with no accepted primer pairs: {}",
                zero_pair_reports
            ));
        }
        if outcome.create_copies {
            status_lines.push(format!(
                "Extracted region copies: {} succeeded, {} failed",
                outcome.created_copy_ids.len(),
                outcome.failed_copy_rows.len()
            ));
            if !outcome.created_copy_ids.is_empty() {
                status_lines.push(format!("Copy IDs: {}", outcome.created_copy_ids.join(", ")));
            }
        }
        if !outcome.failed_report_rows.is_empty() {
            status_lines.push(format!(
                "Primer design failures: {}",
                outcome.failed_report_rows.join(" | ")
            ));
        }
        if !outcome.failed_copy_rows.is_empty() {
            status_lines.push(format!(
                "Copy failures: {}",
                outcome.failed_copy_rows.join(" | ")
            ));
        }
        self.op_status = status_lines.join("\n");
        self.op_error_popup = None;
    }

    #[cfg(test)]
    pub(super) fn run_queued_primer_pair_design_batch(&mut self) {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let prepared = match self.prepare_primer_pair_design_batch_inputs() {
            Ok(value) => value,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        let started = Instant::now();
        let mut guard = match engine.write() {
            Ok(guard) => guard,
            Err(_) => {
                self.op_status =
                    "Engine lock poisoned while running queued PCR primer design".to_string();
                return;
            }
        };
        let outcome = Self::execute_primer_pair_design_batch(
            &mut guard,
            &prepared.queued_regions,
            &prepared.spec,
            &prepared.report_base,
            prepared.create_copies,
            |_| {},
        );
        drop(guard);
        self.apply_primer_pair_design_batch_outcome(outcome, started.elapsed().as_millis());
    }

    pub(super) fn build_design_primer_pairs_operation(
        &self,
        template: &str,
    ) -> Result<Operation, String> {
        let ui = &self.primer_design_ui;
        let max_pairs = Self::parse_optional_usize_text(&ui.max_pairs, "max_pairs")?;
        if matches!(max_pairs, Some(0)) {
            return Err("Invalid max_pairs: expected >= 1 when provided".to_string());
        }
        let (roi_start_0based, roi_end_0based) = self.resolve_roi_range_inputs_0based(
            &ui.roi_start_0based,
            &ui.roi_end_0based,
            "primer_design",
        )?;
        let min_amplicon_bp =
            Self::parse_positive_usize_text(&ui.min_amplicon_bp, "min_amplicon_bp")?;
        Self::validate_pair_pcr_min_amplicon_against_roi(
            min_amplicon_bp,
            roi_start_0based,
            roi_end_0based,
            "pair-PCR",
        )?;
        let max_amplicon_bp =
            Self::parse_positive_usize_text(&ui.max_amplicon_bp, "max_amplicon_bp")?;
        if min_amplicon_bp > max_amplicon_bp {
            return Err(format!(
                "Invalid amplicon window: min_amplicon_bp ({min_amplicon_bp}) must be <= max_amplicon_bp ({max_amplicon_bp})"
            ));
        }
        Ok(Operation::DesignPrimerPairs {
            template: template.to_string(),
            roi_start_0based,
            roi_end_0based,
            forward: Self::parse_primer_side_constraint_ui(&ui.forward, "forward")?,
            reverse: Self::parse_primer_side_constraint_ui(&ui.reverse, "reverse")?,
            pair_constraints: Self::parse_primer_pair_constraint_ui(&ui.pair_constraints)?,
            min_amplicon_bp,
            max_amplicon_bp,
            max_tm_delta_c: Self::parse_optional_f64_text(&ui.max_tm_delta_c, "max_tm_delta_c")?,
            max_pairs,
            report_id: if ui.report_id.trim().is_empty() {
                None
            } else {
                Some(ui.report_id.trim().to_string())
            },
        })
    }

    pub(super) fn build_design_qpcr_operation(&self, template: &str) -> Result<Operation, String> {
        let ui = &self.qpcr_design_ui;
        let max_assays = Self::parse_optional_usize_text(&ui.max_assays, "max_assays")?;
        if matches!(max_assays, Some(0)) {
            return Err("Invalid max_assays: expected >= 1 when provided".to_string());
        }
        let (roi_start_0based, roi_end_0based) = self.resolve_roi_range_inputs_0based(
            &ui.roi_start_0based,
            &ui.roi_end_0based,
            "qpcr_design",
        )?;
        Ok(Operation::DesignQpcrAssays {
            template: template.to_string(),
            roi_start_0based,
            roi_end_0based,
            forward: Self::parse_primer_side_constraint_ui(&ui.forward, "forward")?,
            reverse: Self::parse_primer_side_constraint_ui(&ui.reverse, "reverse")?,
            probe: Self::parse_primer_side_constraint_ui(&ui.probe, "probe")?,
            pair_constraints: Self::parse_primer_pair_constraint_ui(&ui.pair_constraints)?,
            min_amplicon_bp: Self::parse_positive_usize_text(
                &ui.min_amplicon_bp,
                "min_amplicon_bp",
            )?,
            max_amplicon_bp: Self::parse_positive_usize_text(
                &ui.max_amplicon_bp,
                "max_amplicon_bp",
            )?,
            max_tm_delta_c: Self::parse_optional_f64_text(&ui.max_tm_delta_c, "max_tm_delta_c")?,
            max_probe_tm_delta_c: Self::parse_optional_f64_text(
                &ui.max_probe_tm_delta_c,
                "max_probe_tm_delta_c",
            )?,
            max_assays,
            transcript_targeting: self.build_qpcr_transcript_targeting_from_ui(template)?,
            report_id: if ui.report_id.trim().is_empty() {
                None
            } else {
                Some(ui.report_id.trim().to_string())
            },
        })
    }

    pub(super) fn sync_primer_backend_controls_from_engine(&mut self) {
        let Some(engine) = self.engine.clone() else {
            return;
        };
        let Ok(guard) = engine.read() else {
            return;
        };
        self.primer_backend = guard.state().parameters.primer_design_backend;
        let raw = guard.state().parameters.primer3_executable.trim();
        self.primer3_executable = if raw.is_empty() {
            "primer3_core".to_string()
        } else {
            raw.to_string()
        };
    }

    pub(super) fn primer_backend_help_summary(&self) -> &'static str {
        match self.primer_backend {
            PrimerDesignBackend::Auto => {
                "Auto prefers Primer3 when the executable is reachable and otherwise falls back to GENtle's built-in internal backend. Saved reports record the actual requested->used backend."
            }
            PrimerDesignBackend::Internal => {
                "Internal uses GENtle's built-in primer designer only and does not require Primer3."
            }
            PrimerDesignBackend::Primer3 => {
                "Primer3 requires a reachable primer3 executable. When it is unavailable, design requests fail instead of falling back automatically."
            }
        }
    }

    pub(super) fn apply_primer_backend_settings(&mut self) {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let executable = {
            let trimmed = self.primer3_executable.trim();
            if trimmed.is_empty() {
                "primer3_core".to_string()
            } else {
                trimmed.to_string()
            }
        };
        let started = Instant::now();
        let apply_result = {
            let Ok(mut guard) = engine.write() else {
                self.op_status =
                    "Engine lock poisoned while applying primer backend settings".to_string();
                return;
            };
            match guard.apply(Operation::SetParameter {
                name: "primer_design_backend".to_string(),
                value: serde_json::json!(self.primer_backend.as_str()),
            }) {
                Ok(_) => guard.apply(Operation::SetParameter {
                    name: "primer3_executable".to_string(),
                    value: serde_json::json!(executable.clone()),
                }),
                Err(err) => Err(err),
            }
        };
        match apply_result {
            Ok(_) => {
                self.primer3_executable = executable.clone();
                let elapsed_ms = started.elapsed().as_millis();
                self.op_status = format!(
                    "Applied primer backend settings in {} ms (backend={}, primer3_executable='{}')",
                    elapsed_ms,
                    self.primer_backend.as_str(),
                    executable
                );
                self.op_error_popup = None;
            }
            Err(err) => {
                self.op_status =
                    format!("Could not apply primer backend settings: {}", err.message);
                self.op_error_popup = Some(err.message);
            }
        }
    }

    pub(super) fn run_primer3_preflight_probe(&mut self) {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let executable_override = self.primer3_executable.trim().to_string();
        let report = {
            let Ok(guard) = engine.read() else {
                self.op_status = "Engine lock poisoned while running Primer3 preflight".to_string();
                return;
            };
            guard.primer3_preflight_report(
                Some(self.primer_backend),
                Some(executable_override.as_str()),
            )
        };
        let detail = report
            .detail
            .as_deref()
            .or(report.version.as_deref())
            .or(report.error.as_deref())
            .unwrap_or("n/a");
        let configured_display = report
            .configured_executable
            .as_deref()
            .unwrap_or("<default: primer3_core>");
        let effective_display = report.executable.as_str();
        let resolved_display = report.resolved_path.as_deref().unwrap_or("-");
        let working_directory_display = report.working_directory.as_deref().unwrap_or("-");
        let target_display = report.resolved_path.as_deref().unwrap_or(effective_display);
        self.primer3_preflight_status = format!(
            "Primer3 preflight: reachable={} version_probe_ok={} backend={} configured='{}' effective='{}' resolved='{}' cwd='{}' status={} version='{}' detail='{}' probe={} ms",
            report.reachable,
            report.version_probe_ok,
            report.backend,
            configured_display,
            effective_display,
            resolved_display,
            working_directory_display,
            report
                .status_code
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string()),
            report.version.as_deref().unwrap_or("-"),
            detail,
            report.probe_time_ms
        );
        if report.reachable && report.version_probe_ok {
            self.op_status = format!(
                "Primer3 preflight OK (backend={}, executable='{}')",
                report.backend, target_display
            );
            self.op_error_popup = None;
        } else if report.reachable {
            self.op_status = format!(
                "Primer3 preflight warning: version probe returned status {} for '{}'",
                report
                    .status_code
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "signal".to_string()),
                target_display
            );
        } else {
            self.op_status = format!(
                "Primer3 preflight failed for '{}': {}",
                target_display,
                report.error.as_deref().unwrap_or("unreachable")
            );
        }
    }

    pub(super) fn list_primer_design_reports(&mut self) {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let reports = engine
            .read()
            .expect("Engine lock poisoned")
            .list_primer_design_reports();
        if reports.is_empty() {
            self.op_status = "No persisted primer-design reports".to_string();
            return;
        }
        let preview_ids = reports
            .iter()
            .take(8)
            .map(|row| row.report_id.clone())
            .collect::<Vec<_>>();
        let suffix = if reports.len() > preview_ids.len() {
            ", ..."
        } else {
            ""
        };
        self.op_status = format!(
            "Primer reports: {} total [{}{}]",
            reports.len(),
            preview_ids.join(", "),
            suffix
        );
    }

    pub(super) fn show_primer_design_report(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Primer report_id is empty".to_string();
            return;
        }
        let report = match self.load_primer_design_report(report_id) {
            Ok(report) => report,
            Err(message) => {
                self.op_status = message;
                return;
            }
        };
        self.pcr_designer_mode = PcrDesignerMode::PrimerPairs;
        let rejection = &report.rejection_summary;
        let rejection_summary = format!(
            "rejections(window/gc_tm/non_unique/amplicon/primer/pair/eval_skip)={}/{}/{}/{}/{}/{}/{}",
            rejection.out_of_window,
            rejection.gc_or_tm_out_of_bounds,
            rejection.non_unique_anneal,
            rejection.amplicon_or_roi_failure,
            rejection.primer_constraint_failure,
            rejection.pair_constraint_failure,
            rejection.pair_evaluation_limit_skipped
        );
        if report.pair_count == 0 {
            let mut detail_chunks = vec!["NO ACCEPTED PRIMER PAIRS".to_string(), rejection_summary];
            if let Some(explain) = report
                .backend
                .primer3_explain
                .as_deref()
                .map(str::trim)
                .filter(|raw| !raw.is_empty())
            {
                detail_chunks.push(format!("primer3_explain: {explain}"));
            }
            if report.backend.primer3_request_boulder_io.is_some() {
                detail_chunks
                    .push("primer3_request: available (use `Export Primer3 input...`)".to_string());
            }
            self.op_status = format!(
                "Primer report '{}' template='{}' roi={}..{} pairs={} backend={}->{}; {}",
                report.report_id,
                report.template,
                report.roi_start_0based,
                report.roi_end_0based,
                report.pair_count,
                report.backend.requested,
                report.backend.used,
                detail_chunks.join(" | ")
            );
            return;
        }
        let top = report.pairs.first().map(|pair| {
            let geometry = report.pair_core_geometry(pair);
            format!(
                "top amplicon={}..{} len={} left_to_core={} right_to_core={} flanks_core={} score={:.3} 3p_clamp(F/R)={}/{} dimer(3p,max)={}/{}",
                pair.amplicon_start_0based,
                pair.amplicon_end_0based_exclusive,
                pair.amplicon_length_bp,
                geometry.left_label(),
                geometry.right_label(),
                geometry.flanks_core_cleanly(),
                pair.score,
                pair.forward.three_prime_gc_clamp,
                pair.reverse.three_prime_gc_clamp,
                pair.primer_pair_3prime_complementary_run_bp,
                pair.primer_pair_complementary_run_bp
            )
        });
        self.op_status = format!(
            "Primer report '{}' template='{}' roi={}..{} pairs={} backend={}->{}; {} | {}",
            report.report_id,
            report.template,
            report.roi_start_0based,
            report.roi_end_0based,
            report.pair_count,
            report.backend.requested,
            report.backend.used,
            top.unwrap_or_else(|| "top pair unavailable".to_string()),
            rejection_summary
        );
    }

    pub(super) fn confirm_primer_specificity_for_report(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Primer report_id is empty".to_string();
            return;
        }
        let target_genome_id = self.primer_design_ui.specificity_target_genome_id.trim();
        if target_genome_id.is_empty() {
            self.op_status = "Primer specificity requires a prepared target genome id".to_string();
            return;
        }
        let pair_rank = match self
            .primer_design_ui
            .specificity_pair_rank_1based
            .trim()
            .parse::<usize>()
        {
            Ok(rank) if rank > 0 => rank,
            _ => {
                self.op_status = "Primer specificity pair rank must be >= 1".to_string();
                return;
            }
        };
        let max_target_amplicon_bp = match self
            .primer_design_ui
            .specificity_max_target_amplicon_bp
            .trim()
            .parse::<usize>()
        {
            Ok(value) if value > 0 => value,
            _ => {
                self.op_status = "Primer specificity max amplicon must be >= 1".to_string();
                return;
            }
        };
        let max_hits_per_primer = match self
            .primer_design_ui
            .specificity_max_hits_per_primer
            .trim()
            .parse::<usize>()
        {
            Ok(value) if value > 0 => value,
            _ => {
                self.op_status = "Primer specificity max hits per primer must be >= 1".to_string();
                return;
            }
        };
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let policy = PrimerSpecificityPolicy {
            specificity_target_genome_id: Some(target_genome_id.to_string()),
            max_target_amplicon_bp,
            max_hits_per_primer,
            ..PrimerSpecificityPolicy::default()
        };
        let report = match engine
            .read()
            .expect("Engine lock poisoned")
            .assess_primer_pair_specificity(
                Some(report_id),
                Some(pair_rank),
                None,
                None,
                None,
                target_genome_id,
                policy,
                None,
                None,
            ) {
            Ok(report) => report,
            Err(err) => {
                self.op_status = format!(
                    "Primer specificity failed for report '{}' rank {} against '{}': {}",
                    report_id, pair_rank, target_genome_id, err.message
                );
                return;
            }
        };
        self.op_status = format!(
            "Primer specificity {} for '{}' rank {} against '{}': intended={} unintended={} failing_unintended={} primer_hits={} accepted_hits={} warnings={}",
            report.summary.status,
            report_id,
            pair_rank,
            target_genome_id,
            report.summary.intended_amplicon_count,
            report.summary.unintended_amplicon_count,
            report.summary.failing_unintended_amplicon_count,
            report.summary.primer_hit_count,
            report.summary.accepted_primer_hit_count,
            report.warnings.len()
        );
    }

    pub(super) fn load_primer_design_report(
        &self,
        report_id: &str,
    ) -> Result<PrimerDesignReport, String> {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            return Err("Primer report_id is empty".to_string());
        }
        let Some(engine) = self.engine.clone() else {
            return Err("No engine attached".to_string());
        };
        engine
            .read()
            .expect("Engine lock poisoned")
            .get_primer_design_report(report_id)
            .map_err(|err| {
                format!(
                    "Could not load primer report '{report_id}': {}",
                    err.message
                )
            })
    }

    pub(super) fn load_qpcr_design_report(
        &self,
        report_id: &str,
    ) -> Result<QpcrDesignReport, String> {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            return Err("qPCR report_id is empty".to_string());
        }
        let Some(engine) = self.engine.clone() else {
            return Err("No engine attached".to_string());
        };
        engine
            .read()
            .expect("Engine lock poisoned")
            .get_qpcr_design_report(report_id)
            .map_err(|err| format!("Could not load qPCR report '{report_id}': {}", err.message))
    }

    pub(super) fn selected_qpcr_assay_rank_1based_for_report(
        &mut self,
        report: &QpcrDesignReport,
    ) -> Option<usize> {
        let first_rank = report.assays.first().map(|assay| assay.rank)?;
        let selected_rank = self
            .qpcr_design_ui
            .selected_assay_rank_1based
            .trim()
            .parse::<usize>()
            .ok()
            .and_then(|rank| {
                report
                    .assays
                    .iter()
                    .any(|assay| assay.rank == rank)
                    .then_some(rank)
            })
            .unwrap_or(first_rank);
        let normalized = selected_rank.to_string();
        if self.qpcr_design_ui.selected_assay_rank_1based.trim() != normalized {
            self.qpcr_design_ui.selected_assay_rank_1based = normalized;
        }
        Some(selected_rank)
    }

    pub(super) fn qpcr_assay_by_rank(
        report: &QpcrDesignReport,
        rank_1based: usize,
    ) -> Option<&crate::engine::QpcrAssayRecord> {
        report.assays.iter().find(|assay| assay.rank == rank_1based)
    }

    pub(super) fn qpcr_preview_geometry_from_assay(
        report: &QpcrDesignReport,
        assay: &crate::engine::QpcrAssayRecord,
    ) -> Option<QpcrCartoonPreviewGeometry> {
        let forward_window_margin_bp = report
            .roi_start_0based
            .saturating_sub(assay.forward.end_0based_exclusive)
            .max(1);
        let reverse_window_margin_bp = assay
            .reverse
            .start_0based
            .saturating_sub(report.roi_end_0based)
            .max(1);
        let roi_left_bp = assay
            .probe
            .start_0based
            .saturating_sub(report.roi_start_0based)
            .max(1);
        let roi_right_bp = report
            .roi_end_0based
            .saturating_sub(assay.probe.end_0based_exclusive)
            .max(1);
        let probe_window_margin_seed = assay.probe.length_bp.saturating_div(3).max(1);
        Some(QpcrCartoonPreviewGeometry {
            source_label: format!("saved report {} assay #{}", report.report_id, assay.rank),
            roi_left_bp,
            probe_window_left_margin_bp: probe_window_margin_seed,
            probe_site_bp: assay.probe.length_bp.max(1),
            probe_window_right_margin_bp: probe_window_margin_seed,
            roi_right_bp,
            forward_window_margin_bp,
            forward_primer_site_bp: assay.forward.length_bp.max(1),
            reverse_window_margin_bp,
            reverse_primer_site_bp: assay.reverse.length_bp.max(1),
        })
    }

    pub(super) fn relevant_qpcr_splicing_view_for_template(
        &mut self,
        template: &str,
    ) -> Option<SplicingExpertView> {
        if let Some(view) = self.current_splicing_expert_view_for_primary_map()
            && view.seq_id == template
        {
            return Some(view);
        }
        if let Some(view) = self
            .splicing_expert_window_view
            .as_ref()
            .filter(|view| view.seq_id == template)
        {
            return Some((**view).clone());
        }
        self.rna_read_mapping_window_view
            .as_ref()
            .filter(|view| view.seq_id == template)
            .map(|view| (**view).clone())
    }

    pub(super) fn qpcr_splicing_context_summary(
        view: &SplicingExpertView,
        assay: &crate::engine::QpcrAssayRecord,
    ) -> Option<QpcrSplicingContextSummary> {
        let merged_exons = Self::merged_splicing_exon_ranges(view);
        if merged_exons.is_empty() {
            return None;
        }
        let amplicon_start_1based = assay.amplicon_start_0based.saturating_add(1);
        let amplicon_end_1based = assay.amplicon_end_0based_exclusive;
        let probe_start_1based = assay.probe.start_0based.saturating_add(1);
        let probe_end_1based = assay.probe.end_0based_exclusive;
        let amplicon_segments = Self::project_genomic_interval_to_exonic(
            &merged_exons,
            amplicon_start_1based,
            amplicon_end_1based,
        );
        let probe_segments = Self::project_genomic_interval_to_exonic(
            &merged_exons,
            probe_start_1based,
            probe_end_1based,
        );
        let covered_junction_labels = view
            .junctions
            .iter()
            .filter(|junction| {
                amplicon_start_1based <= junction.donor_1based
                    && amplicon_end_1based >= junction.acceptor_1based
            })
            .map(|junction| format!("{}→{}", junction.donor_1based, junction.acceptor_1based))
            .collect::<Vec<_>>();
        let probe_crosses_junction = view.junctions.iter().any(|junction| {
            probe_start_1based <= junction.donor_1based
                && probe_end_1based >= junction.acceptor_1based
        });
        let assay_class_label = if probe_crosses_junction {
            "junction-crossing probe"
        } else if !covered_junction_labels.is_empty() {
            "junction-spanning amplicon"
        } else if amplicon_segments.len() > 1 {
            "multi-exon context"
        } else if amplicon_segments.len() == 1 {
            "single-exon context"
        } else {
            "splicing-region context"
        };
        let explanation = if probe_crosses_junction {
            format!(
                "{}: probe overlaps a modeled exon junction in splicing group '{}'.",
                assay_class_label, view.group_label
            )
        } else if !covered_junction_labels.is_empty() {
            format!(
                "{}: amplicon covers {} modeled junction(s) in group '{}'.",
                assay_class_label,
                covered_junction_labels.len(),
                view.group_label
            )
        } else if amplicon_segments.len() > 1 {
            format!(
                "{}: amplicon touches {} exon segments in group '{}'.",
                assay_class_label,
                amplicon_segments.len(),
                view.group_label
            )
        } else if amplicon_segments.len() == 1 {
            format!(
                "{}: amplicon stays inside one merged exon segment in group '{}'.",
                assay_class_label, view.group_label
            )
        } else if probe_segments.is_empty() {
            format!(
                "{}: assay stays inside the broader splicing ROI for '{}' but outside merged exon bodies.",
                assay_class_label, view.group_label
            )
        } else {
            format!(
                "{}: assay remains inside the saved splicing ROI for group '{}'.",
                assay_class_label, view.group_label
            )
        };
        Some(QpcrSplicingContextSummary {
            assay_class_label: assay_class_label.to_string(),
            explanation,
            covered_junction_labels,
        })
    }

    pub(super) fn qpcr_persisted_context_summary(
        assay: &crate::engine::QpcrAssayRecord,
    ) -> Option<QpcrSplicingContextSummary> {
        let context = assay.transcript_context.as_ref()?;
        let mut explanation = context.explanation.clone();
        explanation.push_str(&format!(
            " Supported transcripts: {} ({:.0}% of considered set). Design transcript: {}.",
            context.support_transcript_count,
            context.support_transcript_fraction * 100.0,
            context.design_transcript_id
        ));
        if let Some(distinguishing_primer) = context.transcript_distinguishing_primer.as_deref() {
            explanation.push_str(&format!(
                " Distinguishing primer(s): {distinguishing_primer}."
            ));
        }
        if let Some(realized_specificity_evidence) = context.realized_specificity_evidence {
            explanation.push_str(&format!(
                " Realized specificity evidence: {}.",
                Self::qpcr_specificity_evidence_label(realized_specificity_evidence)
            ));
        }
        if !context.genomic_carryover_risk.trim().is_empty() {
            explanation.push_str(&format!(
                " Genomic-DNA carryover risk: {}.",
                context.genomic_carryover_risk
            ));
        }
        Some(QpcrSplicingContextSummary {
            assay_class_label: context.assay_class_label.clone(),
            explanation,
            covered_junction_labels: context.covered_junction_labels.clone(),
        })
    }

    pub(super) fn qpcr_report_targeting_summary(report: &QpcrDesignReport) -> Option<String> {
        let targeting = report.transcript_targeting.as_ref()?;
        let mut chunks = vec![format!(
            "requested={}",
            match targeting.mode {
                QpcrTranscriptTargetingMode::SharedGene => "shared across transcripts",
                QpcrTranscriptTargetingMode::DistinguishTranscript => "specific transcript",
            }
        )];
        if let Some(transcript_id) = targeting
            .transcript_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            chunks.push(format!("transcript={transcript_id}"));
        }
        if let Some(specificity_evidence) = targeting.specificity_evidence {
            chunks.push(format!(
                "requested evidence={}",
                Self::qpcr_specificity_evidence_label(specificity_evidence)
            ));
        }
        if let Some(result) = report.transcript_targeting_result.as_ref() {
            chunks.push(format!(
                "support={}/{} ({:.0}%)",
                result.selected_support_transcript_count,
                result.transcript_count_considered,
                result.selected_support_transcript_fraction * 100.0
            ));
            if let Some(transcript_label) = result
                .transcript_label
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                chunks.push(format!("label={transcript_label}"));
            }
            if let Some(realized_specificity_evidence) = result.realized_specificity_evidence {
                chunks.push(format!(
                    "realized evidence={}",
                    Self::qpcr_specificity_evidence_label(realized_specificity_evidence)
                ));
            }
            if result.used_shared_support_fallback {
                chunks.push("used shared-support fallback".to_string());
            }
        }
        Some(chunks.join(" | "))
    }

    pub(super) fn sync_qpcr_transcript_targeting_ui_from_report(
        &mut self,
        report: &QpcrDesignReport,
    ) {
        let targeting_ui = &mut self.qpcr_design_ui.transcript_targeting;
        let Some(targeting) = report.transcript_targeting.as_ref() else {
            *targeting_ui = QpcrTranscriptTargetingUiState::default();
            return;
        };
        targeting_ui.assay_intent = match targeting.mode {
            QpcrTranscriptTargetingMode::SharedGene => {
                QpcrTranscriptIntentUiMode::SharedAcrossTranscripts
            }
            QpcrTranscriptTargetingMode::DistinguishTranscript => {
                QpcrTranscriptIntentUiMode::SpecificTranscript
            }
        };
        targeting_ui.specificity_evidence = targeting
            .specificity_evidence
            .unwrap_or(QpcrTranscriptSpecificityEvidence::JunctionOnly);
        targeting_ui.source_seq_id = report.template.clone();
        targeting_ui.source_feature_id = Some(targeting.source_feature_id);
        targeting_ui.selected_transcript_feature_id = None;
        targeting_ui.selected_transcript_id = targeting.transcript_id.clone().unwrap_or_default();
        targeting_ui.selected_transcript_label = String::new();
        targeting_ui.group_label.clear();
        targeting_ui.transcript_count = 0;
        targeting_ui.strand.clear();
        if let Some(result) = report.transcript_targeting_result.as_ref() {
            targeting_ui.group_label = result.group_label.clone();
            targeting_ui.transcript_count = result.transcript_count_considered;
            targeting_ui.strand = result.strand.clone();
            if targeting_ui.selected_transcript_id.trim().is_empty() {
                targeting_ui.selected_transcript_id =
                    result.transcript_id.clone().unwrap_or_default();
            }
            targeting_ui.selected_transcript_label =
                result.transcript_label.clone().unwrap_or_default();
        }
    }

    pub(super) fn qpcr_preview_geometry_from_ui(
        &self,
    ) -> Result<QpcrCartoonPreviewGeometry, String> {
        let (roi_start_0based, roi_end_0based) = self.resolve_roi_range_inputs_0based(
            &self.qpcr_design_ui.roi_start_0based,
            &self.qpcr_design_ui.roi_end_0based,
            "qpcr_design",
        )?;
        let roi_bp = roi_end_0based.saturating_sub(roi_start_0based).max(1);
        let forward_primer_site_bp = self
            .qpcr_design_ui
            .forward
            .min_length
            .trim()
            .parse::<usize>()
            .ok()
            .unwrap_or(20)
            .max(1);
        let reverse_primer_site_bp = self
            .qpcr_design_ui
            .reverse
            .min_length
            .trim()
            .parse::<usize>()
            .ok()
            .unwrap_or(20)
            .max(1);
        let probe_site_bp = self
            .qpcr_design_ui
            .probe
            .min_length
            .trim()
            .parse::<usize>()
            .ok()
            .unwrap_or(20)
            .max(1);
        let probe_window_left_margin_bp = (probe_site_bp / 4).max(1);
        let probe_window_right_margin_bp = (probe_site_bp / 4).max(1);
        let inner_roi_bp = roi_bp
            .saturating_sub(probe_site_bp)
            .saturating_sub(probe_window_left_margin_bp)
            .saturating_sub(probe_window_right_margin_bp)
            .max(2);
        let roi_left_bp = (inner_roi_bp / 2).max(1);
        let roi_right_bp = inner_roi_bp.saturating_sub(roi_left_bp).max(1);
        let forward_window_margin_bp = (forward_primer_site_bp / 2).max(1);
        let reverse_window_margin_bp = (reverse_primer_site_bp / 2).max(1);
        Ok(QpcrCartoonPreviewGeometry {
            source_label: "current qPCR ROI + constraint defaults".to_string(),
            roi_left_bp,
            probe_window_left_margin_bp,
            probe_site_bp,
            probe_window_right_margin_bp,
            roi_right_bp,
            forward_window_margin_bp,
            forward_primer_site_bp,
            reverse_window_margin_bp,
            reverse_primer_site_bp,
        })
    }

    pub(super) fn load_restriction_cloning_pcr_handoff_report(
        &self,
        report_id: &str,
    ) -> Result<RestrictionCloningPcrHandoffReport, String> {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            return Err("Restriction-cloning handoff report_id is empty".to_string());
        }
        let Some(engine) = self.engine.clone() else {
            return Err("No engine attached".to_string());
        };
        engine
            .read()
            .expect("Engine lock poisoned")
            .get_restriction_cloning_pcr_handoff(report_id)
            .map_err(|err| {
                format!(
                    "Could not load restriction-cloning handoff report '{report_id}': {}",
                    err.message
                )
            })
    }

    pub(super) fn list_restriction_cloning_pcr_handoffs(&mut self) {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let reports = engine
            .read()
            .expect("Engine lock poisoned")
            .list_restriction_cloning_pcr_handoffs();
        if reports.is_empty() {
            self.op_status = "No persisted restriction-cloning PCR handoff reports".to_string();
            return;
        }
        let preview_ids = reports
            .iter()
            .take(8)
            .map(|row| row.report_id.clone())
            .collect::<Vec<_>>();
        let suffix = if reports.len() > preview_ids.len() {
            ", ..."
        } else {
            ""
        };
        self.op_status = format!(
            "Restriction-cloning PCR handoffs: {} total [{}{}]",
            reports.len(),
            preview_ids.join(", "),
            suffix
        );
    }

    pub(super) fn show_restriction_cloning_pcr_handoff_report(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Restriction-cloning handoff report_id is empty".to_string();
            return;
        }
        let report = match self.load_restriction_cloning_pcr_handoff_report(report_id) {
            Ok(report) => report,
            Err(message) => {
                self.op_status = message;
                return;
            }
        };
        self.primer_design_ui
            .restriction_cloning
            .selected_saved_report_id = report.report_id.clone();
        let warning_count = report.compatibility.warnings.len();
        let blocking_count = report.compatibility.blocking_errors.len();
        let workflow_hint_count =
            usize::from(report.workflow_hints.pcr_advanced_operation.is_some())
                + usize::from(report.workflow_hints.staged_workflow.is_some())
                + usize::from(report.workflow_hints.insert_digest_operation.is_some())
                + usize::from(report.workflow_hints.vector_digest_operation.is_some())
                + usize::from(report.workflow_hints.ligation_operation_snippet.is_some());
        self.op_status = format!(
            "Restriction-cloning handoff '{}' template='{}' pair=#{} vector='{}' mode={} enzymes={}/{} status={} warnings={} blocking={} workflow_hints={}",
            report.report_id,
            report.template,
            report.pair_rank.max(report.pair_index + 1),
            report.destination_vector_seq_id,
            report.mode.as_str(),
            report.forward_enzyme,
            report.reverse_enzyme,
            report.compatibility.status,
            warning_count,
            blocking_count,
            workflow_hint_count
        );
    }

    pub(super) fn export_restriction_cloning_pcr_handoff_report_dialog(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Restriction-cloning handoff report_id is empty".to_string();
            return;
        }
        let default_name = format!("{report_id}.restriction_cloning_pcr_handoff.json");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .save_file();
        let Some(path) = path else {
            self.op_status = "Restriction-cloning handoff export canceled".to_string();
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
            .export_restriction_cloning_pcr_handoff(report_id, &path_text);
        match exported {
            Ok(report) => {
                self.op_status = format!(
                    "Exported restriction-cloning handoff '{}' to {}",
                    report.report_id, path_text
                );
            }
            Err(err) => {
                self.op_status = format!(
                    "Could not export restriction-cloning handoff '{report_id}': {}",
                    err.message
                );
            }
        }
    }

    pub(super) fn restriction_cloning_sequence_ids(&self) -> Vec<String> {
        let Some(engine) = self.engine.clone() else {
            return vec![];
        };
        let mut ids = engine
            .read()
            .ok()
            .map(|guard| guard.state().sequences.keys().cloned().collect::<Vec<_>>())
            .unwrap_or_default();
        ids.sort_unstable();
        let active = self.seq_id.as_deref().unwrap_or_default();
        ids.sort_by(|left, right| {
            let left_active = left == active;
            let right_active = right == active;
            left_active
                .cmp(&right_active)
                .reverse()
                .then(left.cmp(right))
        });
        ids
    }

    pub(super) fn restriction_cloning_vector_enzyme_suggestions(
        &self,
        seq_id: &str,
    ) -> Result<RestrictionCloningVectorEnzymeSuggestions, String> {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() {
            return Err("Restriction-cloning destination vector is empty".to_string());
        }
        let Some(engine) = self.engine.clone() else {
            return Err("No engine attached".to_string());
        };
        engine
            .read()
            .expect("Engine lock poisoned")
            .restriction_cloning_vector_enzyme_suggestions(seq_id)
            .map_err(|err| err.message)
    }

    pub(super) fn apply_restriction_cloning_single_site_recommendation(&mut self, enzyme: &str) {
        let enzyme = enzyme.trim();
        if enzyme.is_empty() {
            return;
        }
        self.primer_design_ui.restriction_cloning.mode =
            RestrictionCloningPcrHandoffMode::SingleSite;
        self.primer_design_ui.restriction_cloning.forward_enzyme = enzyme.to_string();
        self.primer_design_ui.restriction_cloning.reverse_enzyme = enzyme.to_string();
        self.primer_design_ui
            .restriction_cloning
            .selected_saved_report_id
            .clear();
        self.save_engine_ops_state();
        let vector = self
            .primer_design_ui
            .restriction_cloning
            .destination_vector_seq_id
            .trim();
        if vector.is_empty() {
            self.op_status = format!(
                "Applied restriction-cloning single-site recommendation '{}'",
                enzyme
            );
        } else {
            self.op_status = format!(
                "Applied restriction-cloning single-site recommendation '{}' for vector '{}'",
                enzyme, vector
            );
        }
    }

    pub(super) fn apply_restriction_cloning_directed_pair_recommendation(
        &mut self,
        forward_enzyme: &str,
        reverse_enzyme: &str,
        order_source: &str,
    ) {
        let forward_enzyme = forward_enzyme.trim();
        let reverse_enzyme = reverse_enzyme.trim();
        if forward_enzyme.is_empty() || reverse_enzyme.is_empty() {
            return;
        }
        self.primer_design_ui.restriction_cloning.mode =
            RestrictionCloningPcrHandoffMode::DirectedPair;
        self.primer_design_ui.restriction_cloning.forward_enzyme = forward_enzyme.to_string();
        self.primer_design_ui.restriction_cloning.reverse_enzyme = reverse_enzyme.to_string();
        self.primer_design_ui
            .restriction_cloning
            .selected_saved_report_id
            .clear();
        self.save_engine_ops_state();
        let vector = self
            .primer_design_ui
            .restriction_cloning
            .destination_vector_seq_id
            .trim();
        let source_label = if order_source.trim().is_empty() {
            "engine ordering".to_string()
        } else {
            order_source.trim().to_string()
        };
        if vector.is_empty() {
            self.op_status = format!(
                "Applied restriction-cloning directed-pair recommendation '{} -> {}' ({})",
                forward_enzyme, reverse_enzyme, source_label
            );
        } else {
            self.op_status = format!(
                "Applied restriction-cloning directed-pair recommendation '{} -> {}' for vector '{}' ({})",
                forward_enzyme, reverse_enzyme, vector, source_label
            );
        }
    }

    pub(super) fn latest_restriction_cloning_handoff_report_id_matching_ui(
        &self,
    ) -> Option<String> {
        let Some(engine) = self.engine.clone() else {
            return None;
        };
        let selected_rank = self
            .primer_design_ui
            .restriction_cloning
            .selected_pair_rank_1based
            .trim()
            .parse::<usize>()
            .ok()
            .and_then(|value| value.checked_sub(1));
        let mode = self.primer_design_ui.restriction_cloning.mode.as_str();
        let vector = self
            .primer_design_ui
            .restriction_cloning
            .destination_vector_seq_id
            .trim()
            .to_string();
        let forward = self
            .primer_design_ui
            .restriction_cloning
            .forward_enzyme
            .trim()
            .to_string();
        let reverse = if self.primer_design_ui.restriction_cloning.mode
            == RestrictionCloningPcrHandoffMode::SingleSite
        {
            if self
                .primer_design_ui
                .restriction_cloning
                .reverse_enzyme
                .trim()
                .is_empty()
            {
                forward.clone()
            } else {
                self.primer_design_ui
                    .restriction_cloning
                    .reverse_enzyme
                    .trim()
                    .to_string()
            }
        } else {
            self.primer_design_ui
                .restriction_cloning
                .reverse_enzyme
                .trim()
                .to_string()
        };
        engine
            .read()
            .ok()
            .map(|guard| {
                guard
                    .list_restriction_cloning_pcr_handoffs()
                    .into_iter()
                    .filter(|summary| {
                        summary.template == self.seq_id.clone().unwrap_or_default()
                            && summary.primer_report_id == self.primer_design_ui.report_id
                            && selected_rank
                                .is_none_or(|pair_index| summary.pair_index == pair_index)
                            && (vector.is_empty()
                                || summary
                                    .destination_vector_seq_id
                                    .eq_ignore_ascii_case(&vector))
                            && (forward.is_empty()
                                || summary.forward_enzyme.eq_ignore_ascii_case(&forward))
                            && (reverse.is_empty()
                                || summary.reverse_enzyme.eq_ignore_ascii_case(&reverse))
                            && summary.mode == mode
                    })
                    .max_by(|left, right| {
                        left.generated_at_unix_ms
                            .cmp(&right.generated_at_unix_ms)
                            .then(left.report_id.cmp(&right.report_id))
                    })
                    .map(|summary| summary.report_id)
            })
            .flatten()
    }

    pub(super) fn seed_restriction_cloning_pcr_handoff_request(
        &self,
    ) -> Result<RestrictionCloningPcrHandoffSeedRequest, String> {
        let ui = &self.primer_design_ui.restriction_cloning;
        let pair_rank = Self::parse_positive_usize_text(
            &ui.selected_pair_rank_1based,
            "restriction_cloning.selected_pair_rank_1based",
        )?;
        let vector_seq_id = ui.destination_vector_seq_id.trim();
        if vector_seq_id.is_empty() {
            return Err(
                "Choose a destination vector sequence before creating a restriction-tail handoff"
                    .to_string(),
            );
        }
        let primer_report_id = self.primer_design_ui.report_id.trim();
        if primer_report_id.is_empty() {
            return Err(
                "Primer report_id is empty; create or select a saved primer report first"
                    .to_string(),
            );
        }
        let Some(engine) = self.engine.clone() else {
            return Err("No engine attached".to_string());
        };
        let forward_enzyme =
            (!ui.forward_enzyme.trim().is_empty()).then_some(ui.forward_enzyme.trim());
        let reverse_enzyme =
            (!ui.reverse_enzyme.trim().is_empty()).then_some(ui.reverse_enzyme.trim());
        engine
            .read()
            .map_err(|_| {
                "Engine lock poisoned while seeding restriction-cloning handoff".to_string()
            })?
            .seed_restriction_cloning_pcr_handoff_request(
                primer_report_id,
                vector_seq_id,
                Some(pair_rank),
                ui.mode,
                forward_enzyme,
                reverse_enzyme,
                (!ui.forward_leader_5prime.trim().is_empty())
                    .then_some(ui.forward_leader_5prime.trim()),
                (!ui.reverse_leader_5prime.trim().is_empty())
                    .then_some(ui.reverse_leader_5prime.trim()),
            )
            .map_err(|err| err.message)
    }

    pub fn focus_restriction_cloning_handoff_report(&mut self, report_id: &str) {
        let normalized_id = report_id.trim();
        if normalized_id.is_empty() {
            self.op_status =
                "Could not open restriction-cloning handoff: report_id is empty".to_string();
            return;
        }
        let report = match self.load_restriction_cloning_pcr_handoff_report(normalized_id) {
            Ok(report) => report,
            Err(message) => {
                self.op_status = message;
                return;
            }
        };
        self.primer_design_ui.report_id = report.primer_report_id.clone();
        self.primer_design_ui
            .restriction_cloning
            .selected_pair_rank_1based = report.pair_rank.max(report.pair_index + 1).to_string();
        self.primer_design_ui
            .restriction_cloning
            .destination_vector_seq_id = report.destination_vector_seq_id.clone();
        self.primer_design_ui.restriction_cloning.mode = report.mode;
        self.primer_design_ui.restriction_cloning.forward_enzyme = report.forward_enzyme.clone();
        self.primer_design_ui.restriction_cloning.reverse_enzyme = report.reverse_enzyme.clone();
        self.primer_design_ui
            .restriction_cloning
            .forward_leader_5prime = report.forward_leader_5prime.clone();
        self.primer_design_ui
            .restriction_cloning
            .reverse_leader_5prime = report.reverse_leader_5prime.clone();
        self.primer_design_ui
            .restriction_cloning
            .selected_saved_report_id = report.report_id.clone();
        self.show_engine_ops = true;
        self.show_restriction_cloning_pcr_handoff_report(normalized_id);
        self.save_engine_ops_state();
    }

    pub(super) fn render_primer_design_report_preview(&mut self, ui: &mut egui::Ui) {
        let report_id = self.primer_design_ui.report_id.trim().to_string();
        if report_id.is_empty() {
            ui.small("Current report_id is empty.");
            return;
        }
        let report = match self.load_primer_design_report(&report_id) {
            Ok(report) => report,
            Err(message) => {
                ui.small(message);
                return;
            }
        };
        ui.group(|ui| {
            ui.label("Primer report preview");
            ui.small(format!(
                "report={} template={} core={}..{} (len {} bp) pairs={} backend={}->{}",
                report.report_id,
                report.template,
                report.roi_start_0based,
                report.roi_end_0based,
                report.roi_end_0based.saturating_sub(report.roi_start_0based),
                report.pair_count,
                report.backend.requested,
                report.backend.used
            ));
            if report.pairs.is_empty() {
                ui.colored_label(
                    egui::Color32::from_rgb(180, 83, 9),
                    "No accepted primer pairs in this saved report.",
                );
                ui.small(format!(
                    "Rejections: window={} gc/tm={} non-unique={} amplicon/ROI={} primer={} pair={} eval-skip={}",
                    report.rejection_summary.out_of_window,
                    report.rejection_summary.gc_or_tm_out_of_bounds,
                    report.rejection_summary.non_unique_anneal,
                    report.rejection_summary.amplicon_or_roi_failure,
                    report.rejection_summary.primer_constraint_failure,
                    report.rejection_summary.pair_constraint_failure,
                    report.rejection_summary.pair_evaluation_limit_skipped
                ));
                return;
            }
            ui.small(
                "Top saved pairs, shown in the simple-PCR vocabulary of distance from the core ROI.",
            );
            egui::Grid::new("primer_report_preview_grid")
                .striped(true)
                .num_columns(8)
                .show(ui, |ui| {
                    ui.strong("#");
                    ui.strong("amplicon");
                    ui.strong("len");
                    ui.strong("left->core");
                    ui.strong("right->core");
                    ui.strong("flanks");
                    ui.strong("ΔTm");
                    ui.strong("score");
                    ui.end_row();
                    for pair in report.pairs.iter().take(5) {
                        let geometry = report.pair_core_geometry(pair);
                        ui.monospace(format!("{}", pair.rank));
                        ui.monospace(format!(
                            "{}..{}",
                            pair.amplicon_start_0based, pair.amplicon_end_0based_exclusive
                        ));
                        ui.monospace(format!("{}", pair.amplicon_length_bp));
                        ui.monospace(geometry.left_label());
                        ui.monospace(geometry.right_label());
                        if geometry.flanks_core_cleanly() {
                            ui.colored_label(
                                egui::Color32::from_rgb(46, 125, 50),
                                "yes",
                            );
                        } else {
                            ui.colored_label(
                                egui::Color32::from_rgb(180, 83, 9),
                                "overlap",
                            );
                        }
                        ui.monospace(format!("{:.1}", pair.tm_delta_c));
                        ui.monospace(format!("{:.1}", pair.score));
                        ui.end_row();
                    }
                });
            if report.pairs.len() > 5 {
                ui.small(format!(
                    "Showing 5 of {} accepted primer pairs. Use `Export report_id...` for the full saved report.",
                    report.pairs.len()
                ));
            }
        });
    }

    pub(super) fn render_qpcr_design_report_preview(&mut self, ui: &mut egui::Ui) {
        let report_id = self.qpcr_design_ui.report_id.trim().to_string();
        if report_id.is_empty() {
            ui.small("Current qPCR report_id is empty.");
            return;
        }
        let report = match self.load_qpcr_design_report(&report_id) {
            Ok(report) => report,
            Err(message) => {
                ui.small(message);
                return;
            }
        };
        let selected_rank = self.selected_qpcr_assay_rank_1based_for_report(&report);
        let selected_assay = selected_rank
            .and_then(|rank| Self::qpcr_assay_by_rank(&report, rank))
            .or_else(|| report.assays.first());
        let splicing_view = self.relevant_qpcr_splicing_view_for_template(report.template.as_str());
        let show_context_column = splicing_view.is_some()
            || report
                .assays
                .iter()
                .any(|assay| assay.transcript_context.is_some());
        let prior_selected_rank = self.qpcr_design_ui.selected_assay_rank_1based.clone();
        ui.group(|ui| {
            ui.label("qPCR report preview");
            ui.small(format!(
                "report={} template={} roi={}..{} (len {} bp) assays={} backend={}->{}",
                report.report_id,
                report.template,
                report.roi_start_0based,
                report.roi_end_0based,
                report.roi_end_0based.saturating_sub(report.roi_start_0based),
                report.assay_count,
                report.backend.requested,
                report.backend.used
            ));
            if let Some(targeting_summary) = Self::qpcr_report_targeting_summary(&report) {
                ui.small(format!("Transcript targeting: {targeting_summary}"));
            }
            if report.assays.is_empty() {
                ui.colored_label(
                    egui::Color32::from_rgb(180, 83, 9),
                    "No accepted qPCR assays in this saved report.",
                );
                ui.small(format!(
                    "Rejections: primer_pair={} probe_window={} probe_gc/tm={} probe_non_unique={} assay={}",
                    report.rejection_summary.primer_pair.pair_constraint_failure,
                    report.rejection_summary.probe_out_of_window,
                    report.rejection_summary.probe_gc_or_tm_out_of_bounds,
                    report.rejection_summary.probe_non_unique_anneal,
                    report.rejection_summary.probe_or_assay_failure
                ));
                return;
            }
            ui.small(
                "Saved assays with forward / reverse / probe geometry. Click a row to drive the live qPCR preview and inspect transcript/splicing context when available.",
            );
            egui::Grid::new("qpcr_report_preview_grid")
                .striped(true)
                .num_columns(if show_context_column { 10 } else { 9 })
                .show(ui, |ui| {
                    ui.strong("#");
                    ui.strong("amplicon");
                    ui.strong("len");
                    ui.strong("forward");
                    ui.strong("probe");
                    ui.strong("reverse");
                    ui.strong("ΔTm primers");
                    ui.strong("ΔTm probe");
                    ui.strong("score");
                    if show_context_column {
                        ui.strong("context");
                    }
                    ui.end_row();
                    for assay in report.assays.iter().take(8) {
                        let selected = Some(assay.rank) == selected_rank;
                        let context_summary = Self::qpcr_persisted_context_summary(assay)
                            .or_else(|| {
                                splicing_view.as_ref().and_then(|view| {
                                    Self::qpcr_splicing_context_summary(view, assay)
                                })
                            });
                        if ui
                            .selectable_label(selected, format!("{}", assay.rank))
                            .on_hover_text("Select this assay and drive the live qPCR preview from it.")
                            .clicked()
                        {
                            self.qpcr_design_ui.selected_assay_rank_1based = assay.rank.to_string();
                        }
                        if ui
                            .selectable_label(
                                selected,
                                format!(
                                    "{}..{}",
                                    assay.amplicon_start_0based, assay.amplicon_end_0based_exclusive
                                ),
                            )
                            .clicked()
                        {
                            self.qpcr_design_ui.selected_assay_rank_1based = assay.rank.to_string();
                        }
                        if ui
                            .selectable_label(selected, format!("{}", assay.amplicon_length_bp))
                            .clicked()
                        {
                            self.qpcr_design_ui.selected_assay_rank_1based = assay.rank.to_string();
                        }
                        if ui
                            .selectable_label(
                                selected,
                                format!(
                                    "{}..{} ({} bp)",
                                    assay.forward.start_0based,
                                    assay.forward.end_0based_exclusive,
                                    assay.forward.length_bp
                                ),
                            )
                            .clicked()
                        {
                            self.qpcr_design_ui.selected_assay_rank_1based = assay.rank.to_string();
                        }
                        if ui
                            .selectable_label(
                                selected,
                                format!(
                                    "{}..{} ({} bp)",
                                    assay.probe.start_0based,
                                    assay.probe.end_0based_exclusive,
                                    assay.probe.length_bp
                                ),
                            )
                            .clicked()
                        {
                            self.qpcr_design_ui.selected_assay_rank_1based = assay.rank.to_string();
                        }
                        if ui
                            .selectable_label(
                                selected,
                                format!(
                                    "{}..{} ({} bp)",
                                    assay.reverse.start_0based,
                                    assay.reverse.end_0based_exclusive,
                                    assay.reverse.length_bp
                                ),
                            )
                            .clicked()
                        {
                            self.qpcr_design_ui.selected_assay_rank_1based = assay.rank.to_string();
                        }
                        if ui
                            .selectable_label(selected, format!("{:.1}", assay.primer_tm_delta_c))
                            .clicked()
                        {
                            self.qpcr_design_ui.selected_assay_rank_1based = assay.rank.to_string();
                        }
                        if ui
                            .selectable_label(selected, format!("{:.1}", assay.probe_tm_delta_c))
                            .clicked()
                        {
                            self.qpcr_design_ui.selected_assay_rank_1based = assay.rank.to_string();
                        }
                        if ui
                            .selectable_label(selected, format!("{:.1}", assay.score))
                            .clicked()
                        {
                            self.qpcr_design_ui.selected_assay_rank_1based = assay.rank.to_string();
                        }
                        if let Some(summary) = context_summary {
                            if ui
                                .selectable_label(selected, summary.assay_class_label)
                                .on_hover_text(summary.explanation)
                                .clicked()
                            {
                                self.qpcr_design_ui.selected_assay_rank_1based =
                                    assay.rank.to_string();
                            }
                        } else if show_context_column {
                            ui.small("-");
                        }
                        ui.end_row();
                    }
                });
            if let Some(selected_assay) = selected_assay {
                ui.separator();
                ui.label(format!("Selected assay #{}", selected_assay.rank));
                ui.small(format!(
                    "amplicon={}..{} ({} bp) | forward={}..{} | probe={}..{} | reverse={}..{} | score={:.1}",
                    selected_assay.amplicon_start_0based,
                    selected_assay.amplicon_end_0based_exclusive,
                    selected_assay.amplicon_length_bp,
                    selected_assay.forward.start_0based,
                    selected_assay.forward.end_0based_exclusive,
                    selected_assay.probe.start_0based,
                    selected_assay.probe.end_0based_exclusive,
                    selected_assay.reverse.start_0based,
                    selected_assay.reverse.end_0based_exclusive,
                    selected_assay.score
                ));
                ui.small(format!(
                    "Current best-assay rationale: highest retained score among assays that keep the probe inside the amplicon and satisfy both primer/probe ΔTm gates."
                ));
                if let Some(summary) = Self::qpcr_persisted_context_summary(selected_assay)
                    .or_else(|| {
                        splicing_view.as_ref().and_then(|view| {
                            Self::qpcr_splicing_context_summary(view, selected_assay)
                        })
                    })
                {
                    ui.small(format!(
                        "Context [{}]: {}",
                        summary.assay_class_label, summary.explanation
                    ));
                    if !summary.covered_junction_labels.is_empty() {
                        ui.small(format!(
                            "Covered modeled junctions: {}",
                            summary.covered_junction_labels.join(", ")
                        ));
                    }
                }
                ui.small(format!(
                    "Selected assay rule flags: roi={} amplicon={} primer ΔT={} probe-inside={} probe ΔT={}",
                    selected_assay.rule_flags.roi_covered,
                    selected_assay.rule_flags.amplicon_size_in_range,
                    selected_assay.rule_flags.primer_tm_delta_in_range,
                    selected_assay.rule_flags.probe_inside_amplicon,
                    selected_assay.rule_flags.probe_tm_delta_in_range
                ));
            }
            if report.assays.len() > 8 {
                ui.small(format!(
                    "Showing 8 of {} accepted qPCR assays. Use `Export report_id...` for the full saved report.",
                    report.assays.len()
                ));
            }
        });
        if self.qpcr_design_ui.selected_assay_rank_1based != prior_selected_rank {
            self.save_engine_ops_state();
        }
    }

    pub(super) fn render_live_qpcr_cartoon_preview(&mut self, ui: &mut egui::Ui) {
        ui.group(|ui| {
            ui.label("Live qPCR cartoon geometry");
            let geometry = self
                .load_qpcr_design_report(&self.qpcr_design_ui.report_id)
                .ok()
                .and_then(|report| {
                    let selected_rank = self.selected_qpcr_assay_rank_1based_for_report(&report);
                    selected_rank
                        .and_then(|rank| Self::qpcr_assay_by_rank(&report, rank))
                        .or_else(|| report.assays.first())
                        .and_then(|assay| Self::qpcr_preview_geometry_from_assay(&report, assay))
                })
                .or_else(|| self.qpcr_preview_geometry_from_ui().ok());
            let Some(geometry) = geometry else {
                ui.small(
                    "Enter a valid qPCR ROI or select a saved qPCR report to derive a live assay-cartoon preview.",
                );
                return;
            };
            ui.monospace(geometry.short_summary());
            ui.small(format!(
                "deterministic bindings ready ({} feature overrides)",
                geometry.bindings_feature_override_count()
            ));
            match pcr_assay_qpcr_spec_with_geometry(
                geometry.roi_left_bp,
                geometry.probe_window_left_margin_bp,
                geometry.probe_site_bp,
                geometry.probe_window_right_margin_bp,
                geometry.roi_right_bp,
                geometry.forward_window_margin_bp,
                geometry.forward_primer_site_bp,
                geometry.reverse_window_margin_bp,
                geometry.reverse_primer_site_bp,
            ) {
                Ok(spec) => {
                    let svg = render_protocol_cartoon_spec_svg(&spec);
                    ui.small(format!(
                        "qPCR cartoon spec resolves with {} event(s); svg ready ({} chars).",
                        spec.events.len(),
                        svg.len()
                    ));
                }
                Err(err) => {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!("Could not resolve qPCR cartoon: {err}"),
                    );
                }
            }
        });
    }

    pub(super) fn render_restriction_cloning_handoff_section(
        &mut self,
        ui: &mut egui::Ui,
        template: &str,
    ) {
        let report = match self.load_primer_design_report(&self.primer_design_ui.report_id) {
            Ok(report) => report,
            Err(message) => {
                ui.group(|ui| {
                    ui.label("Restriction-site cloning handoff");
                    ui.small(message);
                    ui.small(
                        "Create or select a saved primer report first, then extend one accepted pair with restriction-site tails for cloning preflight.",
                    );
                });
                return;
            }
        };
        let pair_count = report.pairs.len();
        let pair_preview = report
            .pairs
            .iter()
            .take(8)
            .map(|pair| {
                format!(
                    "#{} {}..{}",
                    pair.rank, pair.amplicon_start_0based, pair.amplicon_end_0based_exclusive
                )
            })
            .collect::<Vec<_>>();
        let sequence_ids = self.restriction_cloning_sequence_ids();
        let vector_options = sequence_ids
            .iter()
            .filter(|seq_id| seq_id.as_str() != template)
            .cloned()
            .collect::<Vec<_>>();
        if self
            .primer_design_ui
            .restriction_cloning
            .destination_vector_seq_id
            .trim()
            .is_empty()
            && let Some(first_vector) = vector_options.first()
        {
            self.primer_design_ui
                .restriction_cloning
                .destination_vector_seq_id = first_vector.clone();
        }
        let suggestion_result = if self
            .primer_design_ui
            .restriction_cloning
            .destination_vector_seq_id
            .trim()
            .is_empty()
        {
            None
        } else {
            Some(
                self.restriction_cloning_vector_enzyme_suggestions(
                    &self
                        .primer_design_ui
                        .restriction_cloning
                        .destination_vector_seq_id,
                ),
            )
        };
        let mut suggested_enzymes = vec![];
        let mut suggestion_note = None::<String>;
        let mut recommended_single_site = vec![];
        let mut recommended_directed_pairs = vec![];
        if let Some(result) = suggestion_result {
            match result {
                Ok(suggestions) => {
                    suggested_enzymes.extend(suggestions.selected_mcs.iter().cloned());
                    for enzyme in &suggestions.other_unique {
                        if !suggested_enzymes.iter().any(|value| value == enzyme) {
                            suggested_enzymes.push(enzyme.clone());
                        }
                    }
                    recommended_single_site = suggestions.recommended_single_site.clone();
                    recommended_directed_pairs = suggestions.recommended_directed_pairs.clone();
                    if !suggestions.missing_mcs.is_empty() {
                        suggestion_note = Some(format!(
                            "MCS annotations mention unavailable/non-unique sites: {}",
                            suggestions.missing_mcs.join(", ")
                        ));
                    }
                }
                Err(err) => suggestion_note = Some(err),
            }
        }
        if !self
            .primer_design_ui
            .restriction_cloning
            .forward_enzyme
            .trim()
            .is_empty()
            && !suggested_enzymes.iter().any(|value| {
                value.eq_ignore_ascii_case(
                    self.primer_design_ui
                        .restriction_cloning
                        .forward_enzyme
                        .trim(),
                )
            })
        {
            suggested_enzymes.push(
                self.primer_design_ui
                    .restriction_cloning
                    .forward_enzyme
                    .trim()
                    .to_string(),
            );
        }
        if !self
            .primer_design_ui
            .restriction_cloning
            .reverse_enzyme
            .trim()
            .is_empty()
            && !suggested_enzymes.iter().any(|value| {
                value.eq_ignore_ascii_case(
                    self.primer_design_ui
                        .restriction_cloning
                        .reverse_enzyme
                        .trim(),
                )
            })
        {
            suggested_enzymes.push(
                self.primer_design_ui
                    .restriction_cloning
                    .reverse_enzyme
                    .trim()
                    .to_string(),
            );
        }
        suggested_enzymes.sort_unstable();
        suggested_enzymes.dedup();

        let selected_report_id = self
            .primer_design_ui
            .restriction_cloning
            .selected_saved_report_id
            .trim()
            .to_string();
        let selected_report = if selected_report_id.is_empty() {
            None
        } else {
            self.load_restriction_cloning_pcr_handoff_report(&selected_report_id)
                .ok()
        };

        ui.group(|ui| {
            ui.label("Restriction-site cloning handoff");
            ui.small(
                "Post-design cloning preflight: extend one accepted primer pair with restriction-site tails, create new primer-vial artifacts plus a tailed amplicon, and stage digest/ligation hints without running cloning automatically.",
            );
            ui.small(format!(
                "Saved primer report '{}' has {} accepted pair(s) [{}{}]",
                report.report_id,
                pair_count,
                pair_preview.join(", "),
                if pair_count > pair_preview.len() {
                    ", ..."
                } else {
                    ""
                }
            ));
            egui::Grid::new("restriction_cloning_handoff_grid")
                .num_columns(4)
                .spacing([12.0, 6.0])
                .show(ui, |ui| {
                    ui.label("pair rank");
                    ui.add(
                        egui::TextEdit::singleline(
                            &mut self
                                .primer_design_ui
                                .restriction_cloning
                                .selected_pair_rank_1based,
                        )
                        .desired_width(80.0),
                    )
                    .on_hover_text("Accepted primer-pair rank from the saved primer report (1-based).");
                    ui.label("destination vector");
                    ui.add(
                        egui::TextEdit::singleline(
                            &mut self
                                .primer_design_ui
                                .restriction_cloning
                                .destination_vector_seq_id,
                        )
                        .desired_width(220.0),
                    )
                    .on_hover_text("Sequence ID of the destination cloning vector whose restriction sites constrain this handoff.");
                    ui.end_row();
                    ui.label("mode");
                    egui::ComboBox::from_id_salt("restriction_cloning_mode")
                        .selected_text(
                            self.primer_design_ui
                                .restriction_cloning
                                .mode
                                .as_str(),
                        )
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.primer_design_ui.restriction_cloning.mode,
                                RestrictionCloningPcrHandoffMode::SingleSite,
                                "single_site",
                            );
                            ui.selectable_value(
                                &mut self.primer_design_ui.restriction_cloning.mode,
                                RestrictionCloningPcrHandoffMode::DirectedPair,
                                "directed_pair",
                            );
                        });
                    ui.label("saved handoff");
                    ui.add(
                        egui::TextEdit::singleline(
                            &mut self
                                .primer_design_ui
                                .restriction_cloning
                                .selected_saved_report_id,
                        )
                        .desired_width(220.0),
                    )
                    .on_hover_text("Optional persisted restriction-cloning handoff report_id to reopen/export.");
                    ui.end_row();
                    ui.label("forward enzyme");
                    ui.add(
                        egui::TextEdit::singleline(
                            &mut self.primer_design_ui.restriction_cloning.forward_enzyme,
                        )
                        .desired_width(120.0),
                    )
                    .on_hover_text("Restriction enzyme on the forward-primer side (must be a usable vector cutter).");
                    ui.label("reverse enzyme");
                    ui.add_enabled(
                        self.primer_design_ui.restriction_cloning.mode
                            == RestrictionCloningPcrHandoffMode::DirectedPair,
                        egui::TextEdit::singleline(
                            &mut self.primer_design_ui.restriction_cloning.reverse_enzyme,
                        )
                        .desired_width(120.0),
                    )
                    .on_hover_text("Restriction enzyme on the reverse-primer side. In single_site mode the same enzyme is reused.");
                    ui.end_row();
                    ui.label("forward 5' leader");
                    ui.add(
                        egui::TextEdit::singleline(
                            &mut self
                                .primer_design_ui
                                .restriction_cloning
                                .forward_leader_5prime,
                        )
                        .desired_width(120.0),
                    )
                    .on_hover_text("Optional additional 5' bases placed before the forward restriction site.");
                    ui.label("reverse 5' leader");
                    ui.add(
                        egui::TextEdit::singleline(
                            &mut self
                                .primer_design_ui
                                .restriction_cloning
                                .reverse_leader_5prime,
                        )
                        .desired_width(120.0),
                    )
                    .on_hover_text("Optional additional 5' bases placed before the reverse restriction site.");
                    ui.end_row();
                });

            if self.primer_design_ui.restriction_cloning.mode
                == RestrictionCloningPcrHandoffMode::SingleSite
            {
                self.primer_design_ui.restriction_cloning.reverse_enzyme =
                    self.primer_design_ui.restriction_cloning.forward_enzyme.clone();
            }

            if !vector_options.is_empty() {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Vector shortcuts");
                    egui::ComboBox::from_id_salt("restriction_cloning_vector_combo")
                        .selected_text(
                            if self
                                .primer_design_ui
                                .restriction_cloning
                                .destination_vector_seq_id
                                .trim()
                                .is_empty()
                            {
                                "choose destination vector"
                            } else {
                                self.primer_design_ui
                                    .restriction_cloning
                                    .destination_vector_seq_id
                                    .as_str()
                            },
                        )
                        .show_ui(ui, |ui| {
                            for seq_id in &vector_options {
                                ui.selectable_value(
                                    &mut self
                                        .primer_design_ui
                                        .restriction_cloning
                                        .destination_vector_seq_id,
                                    seq_id.clone(),
                                    seq_id,
                                );
                            }
                        });
                });
            }
            if !recommended_single_site.is_empty() {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Recommended single-site");
                    for suggestion in &recommended_single_site {
                        let button_label = format!(
                            "{} @ {}",
                            suggestion.enzyme, suggestion.cut_position_0based
                        );
                        if ui
                            .button(&button_label)
                            .on_hover_text(
                                "Use this MCS-aware single-cutter recommendation for both primer tails",
                            )
                            .clicked()
                        {
                            self.apply_restriction_cloning_single_site_recommendation(
                                &suggestion.enzyme,
                            );
                        }
                    }
                });
            }
            if !recommended_directed_pairs.is_empty() {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Recommended directed pairs");
                    for suggestion in &recommended_directed_pairs {
                        let button_label = format!(
                            "{} -> {}",
                            suggestion.forward_enzyme, suggestion.reverse_enzyme
                        );
                        let hover = format!(
                            "Use directed insertion order from {} (cuts at {} and {})",
                            suggestion.order_source,
                            suggestion.forward_cut_position_0based,
                            suggestion.reverse_cut_position_0based
                        );
                        if ui
                            .button(&button_label)
                            .on_hover_text(hover)
                            .clicked()
                        {
                            self.apply_restriction_cloning_directed_pair_recommendation(
                                &suggestion.forward_enzyme,
                                &suggestion.reverse_enzyme,
                                &suggestion.order_source,
                            );
                        }
                    }
                });
            }
            if !suggested_enzymes.is_empty() {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Suggested enzymes");
                    egui::ComboBox::from_id_salt("restriction_cloning_forward_combo")
                        .selected_text(
                            if self
                                .primer_design_ui
                                .restriction_cloning
                                .forward_enzyme
                                .trim()
                                .is_empty()
                            {
                                "choose forward enzyme"
                            } else {
                                self.primer_design_ui
                                    .restriction_cloning
                                    .forward_enzyme
                                    .as_str()
                            },
                        )
                        .show_ui(ui, |ui| {
                            for enzyme in &suggested_enzymes {
                                ui.selectable_value(
                                    &mut self.primer_design_ui.restriction_cloning.forward_enzyme,
                                    enzyme.clone(),
                                    enzyme,
                                );
                            }
                        });
                    if self.primer_design_ui.restriction_cloning.mode
                        == RestrictionCloningPcrHandoffMode::DirectedPair
                    {
                        egui::ComboBox::from_id_salt("restriction_cloning_reverse_combo")
                            .selected_text(
                                if self
                                    .primer_design_ui
                                    .restriction_cloning
                                    .reverse_enzyme
                                    .trim()
                                    .is_empty()
                                {
                                    "choose reverse enzyme"
                                } else {
                                    self.primer_design_ui
                                        .restriction_cloning
                                        .reverse_enzyme
                                        .as_str()
                                },
                            )
                            .show_ui(ui, |ui| {
                                for enzyme in &suggested_enzymes {
                                    ui.selectable_value(
                                        &mut self
                                            .primer_design_ui
                                            .restriction_cloning
                                            .reverse_enzyme,
                                        enzyme.clone(),
                                        enzyme,
                                    );
                                }
                            });
                    }
                });
            }
            if let Some(note) = suggestion_note {
                ui.small(note);
            } else if !suggested_enzymes.is_empty() {
                ui.small(format!(
                    "Vector suggestions prefer annotated MCS cutters first, then fall back to other unique cutters: {}",
                    suggested_enzymes.join(", ")
                ));
            }
            if !recommended_directed_pairs.is_empty() {
                ui.small(
                    "Directed-pair recommendations follow MCS order when available, otherwise unique-cut position order on the vector.",
                );
            }

            ui.horizontal_wrapped(|ui| {
                if ui
                    .button("Create restriction-tail handoff")
                    .on_hover_text(
                        "Create extended primer-vial artifacts plus a tailed amplicon and digest/ligation preflight hints for the selected primer pair",
                    )
                    .clicked()
                {
                    match self.seed_restriction_cloning_pcr_handoff_request() {
                        Ok(seed_request) => {
                            let started = Instant::now();
                            let apply_result =
                                self.engine
                                    .clone()
                                    .ok_or_else(|| "No engine attached".to_string())
                                    .and_then(|engine| {
                                        engine
                                            .write()
                                            .map_err(|_| {
                                                "Engine lock poisoned while preparing restriction-cloning handoff".to_string()
                                            })?
                                            .apply(seed_request.operation.clone())
                                            .map_err(|err| err.message)
                                    });
                            match apply_result {
                                Ok(result) => {
                                    self.handle_operation_success(result, started);
                                    if let Some(report_id) =
                                        self.latest_restriction_cloning_handoff_report_id_matching_ui()
                                    {
                                        self.focus_restriction_cloning_handoff_report(&report_id);
                                    } else {
                                        self.save_engine_ops_state();
                                    }
                                }
                                Err(message) => {
                                    self.op_status =
                                        format!("Restriction-cloning handoff failed: {message}");
                                    self.op_error_popup = Some(message);
                                }
                            }
                        }
                        Err(err) => self.op_status = err,
                    }
                }
                if ui
                    .button("List saved handoffs")
                    .on_hover_text("List persisted restriction-cloning PCR handoff report IDs")
                    .clicked()
                {
                    self.list_restriction_cloning_pcr_handoffs();
                }
                let show_enabled = !self
                    .primer_design_ui
                    .restriction_cloning
                    .selected_saved_report_id
                    .trim()
                    .is_empty();
                if ui
                    .add_enabled(show_enabled, egui::Button::new("Show saved handoff"))
                    .on_hover_text("Show summary for the selected saved handoff report")
                    .clicked()
                {
                    let report_id = self
                        .primer_design_ui
                        .restriction_cloning
                        .selected_saved_report_id
                        .clone();
                    self.show_restriction_cloning_pcr_handoff_report(&report_id);
                }
                if ui
                    .add_enabled(show_enabled, egui::Button::new("Export saved handoff..."))
                    .on_hover_text("Export the selected saved handoff report to JSON")
                    .clicked()
                {
                    let report_id = self
                        .primer_design_ui
                        .restriction_cloning
                        .selected_saved_report_id
                        .clone();
                    self.export_restriction_cloning_pcr_handoff_report_dialog(&report_id);
                }
            });

            if let Some(saved) = selected_report {
                ui.separator();
                ui.label("Saved handoff preview");
                ui.small(format!(
                    "report={} pair=#{} vector={} mode={} status={}",
                    saved.report_id,
                    saved.pair_rank.max(saved.pair_index + 1),
                    saved.destination_vector_seq_id,
                    saved.mode.as_str(),
                    saved.compatibility.status
                ));
                egui::Grid::new("restriction_cloning_saved_detail_grid")
                    .num_columns(2)
                    .spacing([10.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Extended forward");
                        ui.monospace(&saved.extended_forward.sequence);
                        ui.end_row();
                        ui.label("Extended reverse");
                        ui.monospace(&saved.extended_reverse.sequence);
                        ui.end_row();
                        ui.label("Tailed amplicon ends");
                        ui.monospace(format!(
                            "5' {} | 3' {} | len {} bp",
                            saved.predicted_tailed_amplicon_5prime,
                            saved.predicted_tailed_amplicon_3prime,
                            saved.predicted_tailed_amplicon_length_bp
                        ));
                        ui.end_row();
                        ui.label("Created artifacts");
                        ui.monospace(format!(
                            "fwd={} | rev={} | amplicon={}",
                            saved.extended_forward_seq_id,
                            saved.extended_reverse_seq_id,
                            saved.tailed_amplicon_seq_id
                        ));
                        ui.end_row();
                        ui.label("Digest geometry");
                        ui.monospace(format!(
                            "forward={} reverse={} order={}",
                            saved.compatibility.forward_end_geometry,
                            saved.compatibility.reverse_end_geometry,
                            saved.compatibility.order_source
                        ));
                        ui.end_row();
                    });
                if !saved.compatibility.blocking_errors.is_empty() {
                    ui.colored_label(
                        egui::Color32::from_rgb(183, 28, 28),
                        format!(
                            "Blocking errors: {}",
                            saved.compatibility.blocking_errors.join(" | ")
                        ),
                    );
                }
                if !saved.compatibility.warnings.is_empty() {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 83, 9),
                        format!("Warnings: {}", saved.compatibility.warnings.join(" | ")),
                    );
                }
                ui.small(format!(
                    "Suggested next steps: staged workflow={} PCR advanced={} insert digest={} vector digest={} ligation snippet={}",
                    saved.workflow_hints.staged_workflow.is_some(),
                    saved.workflow_hints.pcr_advanced_operation.is_some(),
                    saved.workflow_hints.insert_digest_operation.is_some(),
                    saved.workflow_hints.vector_digest_operation.is_some(),
                    saved.workflow_hints.ligation_operation_snippet.is_some()
                ));
            }
        });
    }

    pub(super) fn export_primer_design_report_dialog(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Primer report_id is empty".to_string();
            return;
        }
        let default_name = format!("{report_id}.primer_report.json");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .save_file();
        let Some(path) = path else {
            self.op_status = "Primer report export canceled".to_string();
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
            .export_primer_design_report(report_id, &path_text);
        match exported {
            Ok(report) => {
                self.op_status = format!(
                    "Exported primer report '{}' (pairs={}) to {}",
                    report.report_id, report.pair_count, path_text
                );
            }
            Err(err) => {
                self.op_status = format!(
                    "Could not export primer report '{report_id}': {}",
                    err.message
                );
            }
        }
    }

    pub(super) fn export_primer3_request_dialog(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "Primer report_id is empty".to_string();
            return;
        }
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let report = match engine
            .read()
            .expect("Engine lock poisoned")
            .get_primer_design_report(report_id)
        {
            Ok(report) => report,
            Err(err) => {
                self.op_status = format!(
                    "Could not load primer report '{report_id}': {}",
                    err.message
                );
                return;
            }
        };
        let Some(request_payload) = report.backend.primer3_request_boulder_io.as_deref() else {
            self.op_status = format!(
                "Primer3 request payload is unavailable for report '{}' (backend used='{}')",
                report.report_id, report.backend.used
            );
            return;
        };
        let default_name = format!("{report_id}.primer3.boulderio.txt");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .save_file();
        let Some(path) = path else {
            self.op_status = "Primer3 request export canceled".to_string();
            return;
        };
        let path_text = path.to_string_lossy().to_string();
        match std::fs::write(&path, request_payload) {
            Ok(()) => {
                let executable = report
                    .backend
                    .primer3_executable
                    .as_deref()
                    .filter(|raw| !raw.trim().is_empty())
                    .unwrap_or("primer3_core");
                self.op_status = format!(
                    "Exported Primer3 request for report '{}' to {} (rerun: {} < '{}')",
                    report.report_id, path_text, executable, path_text
                );
            }
            Err(err) => {
                self.op_status = format!(
                    "Could not export Primer3 request for report '{}': {}",
                    report.report_id, err
                );
            }
        }
    }

    pub(super) fn render_protocol_cartoon_preview_svg(
        preview: &ProtocolCartoonPreviewTelemetry,
    ) -> Result<String, String> {
        let protocol = ProtocolCartoonKind::parse_id(&preview.protocol).ok_or_else(|| {
            format!(
                "Unknown protocol-cartoon id '{}' in preview payload",
                preview.protocol
            )
        })?;
        let template = protocol_cartoon_template_for_kind(&protocol);
        let spec = resolve_protocol_cartoon_template_with_bindings(&template, &preview.bindings)
            .map_err(|e| {
                format!(
                    "Could not resolve protocol-cartoon preview '{}' with bindings: {e}",
                    preview.protocol
                )
            })?;
        Ok(render_protocol_cartoon_spec_svg(&spec))
    }

    pub(super) fn default_protocol_cartoon_preview_file_name(&self) -> String {
        let protocol = self
            .last_protocol_cartoon_preview
            .as_ref()
            .map(|preview| preview.protocol.trim())
            .filter(|raw| !raw.is_empty())
            .unwrap_or("protocol_cartoon");
        let protocol_token = Self::sanitize_export_name_component(protocol, "protocol_cartoon");
        let op_token =
            Self::sanitize_export_name_component(&self.last_protocol_cartoon_preview_op_id, "op");
        format!("{protocol_token}.{op_token}.preview.svg")
    }

    pub(super) fn export_protocol_cartoon_preview_svg_dialog(&mut self) {
        let Some(preview) = self.last_protocol_cartoon_preview.as_ref() else {
            self.op_status =
                "No protocol-cartoon preview is available yet; run an OE mutagenesis operation first"
                    .to_string();
            return;
        };
        let default_name = self.default_protocol_cartoon_preview_file_name();
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.op_status = "Protocol-cartoon preview export canceled".to_string();
            return;
        };
        let svg = match Self::render_protocol_cartoon_preview_svg(preview) {
            Ok(svg) => svg,
            Err(err) => {
                self.op_status = err;
                return;
            }
        };
        match fs::write(&path, svg) {
            Ok(()) => {
                self.op_status = format!(
                    "Exported protocol-cartoon preview '{}' (op={}) to {}",
                    preview.protocol,
                    self.last_protocol_cartoon_preview_op_id,
                    path.display()
                );
            }
            Err(err) => {
                self.op_status = format!(
                    "Could not write protocol-cartoon preview SVG to '{}': {}",
                    path.display(),
                    err
                );
            }
        }
    }

    pub(super) fn export_builtin_protocol_cartoon_svg_dialog(
        &mut self,
        protocol: ProtocolCartoonKind,
        fallback_name: &str,
    ) {
        let default_name = format!(
            "{}.svg",
            Self::sanitize_export_name_component(protocol.id(), fallback_name)
        );
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("SVG", &["svg"])
            .save_file();
        let Some(path) = path else {
            self.op_status = format!("{} SVG export canceled", protocol.id());
            return;
        };
        let svg = render_protocol_cartoon_svg(&protocol);
        match fs::write(&path, svg) {
            Ok(()) => {
                self.op_status = format!(
                    "Exported built-in protocol cartoon '{}' to {}",
                    protocol.id(),
                    path.display()
                );
            }
            Err(err) => {
                self.op_status = format!(
                    "Could not export protocol cartoon '{}' to '{}': {}",
                    protocol.id(),
                    path.display(),
                    err
                );
            }
        }
    }

    pub(super) fn list_qpcr_design_reports(&mut self) {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let reports = engine
            .read()
            .expect("Engine lock poisoned")
            .list_qpcr_design_reports();
        if reports.is_empty() {
            self.op_status = "No persisted qPCR design reports".to_string();
            return;
        }
        let preview_ids = reports
            .iter()
            .take(8)
            .map(|row| row.report_id.clone())
            .collect::<Vec<_>>();
        let suffix = if reports.len() > preview_ids.len() {
            ", ..."
        } else {
            ""
        };
        self.op_status = format!(
            "qPCR reports: {} total [{}{}]",
            reports.len(),
            preview_ids.join(", "),
            suffix
        );
    }

    pub(super) fn show_qpcr_design_report(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "qPCR report_id is empty".to_string();
            return;
        }
        let report = match self.load_qpcr_design_report(report_id) {
            Ok(report) => report,
            Err(message) => {
                self.op_status = message;
                return;
            }
        };
        self.pcr_designer_mode = PcrDesignerMode::QpcrAssays;
        let selected = self
            .selected_qpcr_assay_rank_1based_for_report(&report)
            .and_then(|rank| Self::qpcr_assay_by_rank(&report, rank))
            .or_else(|| report.assays.first())
            .map(|assay| {
                format!(
                    "selected assay #{} amplicon={}..{} len={} score={:.3}",
                    assay.rank,
                    assay.amplicon_start_0based,
                    assay.amplicon_end_0based_exclusive,
                    assay.amplicon_length_bp,
                    assay.score
                )
            })
            .unwrap_or_else(|| "no accepted qPCR assays".to_string());
        self.sync_qpcr_transcript_targeting_ui_from_report(&report);
        let targeting_summary = Self::qpcr_report_targeting_summary(&report)
            .map(|summary| format!(" | {summary}"))
            .unwrap_or_default();
        self.op_status = format!(
            "qPCR report '{}' template='{}' roi={}..{} assays={} backend={}->{}; {}{}",
            report.report_id,
            report.template,
            report.roi_start_0based,
            report.roi_end_0based,
            report.assay_count,
            report.backend.requested,
            report.backend.used,
            selected,
            targeting_summary
        );
    }

    pub(super) fn export_qpcr_design_report_dialog(&mut self, report_id: &str) {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            self.op_status = "qPCR report_id is empty".to_string();
            return;
        }
        let default_name = format!("{report_id}.qpcr_report.json");
        let path = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .save_file();
        let Some(path) = path else {
            self.op_status = "qPCR report export canceled".to_string();
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
            .export_qpcr_design_report(report_id, &path_text);
        match exported {
            Ok(report) => {
                self.op_status = format!(
                    "Exported qPCR report '{}' (assays={}) to {}",
                    report.report_id, report.assay_count, path_text
                );
            }
            Err(err) => {
                self.op_status = format!(
                    "Could not export qPCR report '{report_id}': {}",
                    err.message
                );
            }
        }
    }

    pub(super) fn render_primer_side_constraint_editor(
        ui: &mut egui::Ui,
        id_prefix: &str,
        label: &str,
        side: &mut PrimerSideConstraintUiState,
    ) {
        const MOTIF_HELP: &str = "Comma-separated motifs using IUPAC DNA codes (A,C,G,T,R,Y,S,W,K,M,B,D,H,V,N), for example `ATG,GGNCC,TTTRAA`.";
        ui.group(|ui| {
            ui.label(label);
            egui::Grid::new(format!("{id_prefix}_grid_basic"))
                .num_columns(8)
                .show(ui, |ui| {
                    ui.label("len min")
                        .on_hover_text("Minimum annealing length in bp for this primer side.");
                    ui.add(egui::TextEdit::singleline(&mut side.min_length).desired_width(56.0));
                    ui.label("len max")
                        .on_hover_text("Maximum annealing length in bp for this primer side.");
                    ui.add(egui::TextEdit::singleline(&mut side.max_length).desired_width(56.0));
                    ui.label("max anneal hits").on_hover_text(
                        "Maximum allowed genomic/template annealing-hit count for the annealing segment.",
                    );
                    ui.add(
                        egui::TextEdit::singleline(&mut side.max_anneal_hits).desired_width(70.0),
                    );
                    ui.end_row();

                    ui.label("location").on_hover_text(
                        "Optional exact annealing anchor (0-based) for this primer side.",
                    );
                    ui.add(egui::TextEdit::singleline(&mut side.location_0based).desired_width(76.0))
                        .on_hover_text("Exact 0-based annealing position when a fixed anchor is required.");
                    ui.label("start").on_hover_text(
                        "Optional minimum annealing start position (0-based).",
                    );
                    ui.add(egui::TextEdit::singleline(&mut side.start_0based).desired_width(76.0))
                        .on_hover_text("Lower bound for annealing start (0-based).");
                    ui.label("end").on_hover_text(
                        "Optional maximum annealing end position (0-based, end-exclusive).",
                    );
                    ui.add(egui::TextEdit::singleline(&mut side.end_0based).desired_width(76.0))
                        .on_hover_text("Upper bound for annealing end (0-based, end-exclusive).");
                    ui.end_row();

                    ui.label("Tm min")
                        .on_hover_text("Minimum annealing melting temperature (°C).");
                    ui.add(egui::TextEdit::singleline(&mut side.min_tm_c).desired_width(64.0));
                    ui.label("Tm max")
                        .on_hover_text("Maximum annealing melting temperature (°C).");
                    ui.add(egui::TextEdit::singleline(&mut side.max_tm_c).desired_width(64.0));
                    ui.label("GC min")
                        .on_hover_text("Minimum GC fraction for annealing segment (0.0..1.0).");
                    ui.add(egui::TextEdit::singleline(&mut side.min_gc_fraction).desired_width(64.0));
                    ui.label("GC max")
                        .on_hover_text("Maximum GC fraction for annealing segment (0.0..1.0).");
                    ui.add(egui::TextEdit::singleline(&mut side.max_gc_fraction).desired_width(64.0));
                    ui.end_row();
            });
            ui.horizontal(|ui| {
                ui.label("5' tail (non-annealing)").on_hover_text(
                    "Optional 5' tail sequence added to the oligo but excluded from annealing constraints.",
                );
                ui.add(
                    egui::TextEdit::singleline(&mut side.non_annealing_5prime_tail)
                        .desired_width(180.0)
                        .hint_text("e.g. GAATTC"),
                )
                .on_hover_text(
                    "Optional 5' primer tail that is added to the oligo but excluded from anneal Tm/GC/hit calculations",
                );
            });
            ui.horizontal(|ui| {
                ui.label("fixed 5'")
                    .on_hover_text("Prefix motif that must appear at primer 5' end (IUPAC supported).");
                ui.add(egui::TextEdit::singleline(&mut side.fixed_5prime).desired_width(120.0))
                    .on_hover_text(MOTIF_HELP);
                ui.label("fixed 3'")
                    .on_hover_text("Suffix motif that must appear at primer 3' end (IUPAC supported).");
                ui.add(egui::TextEdit::singleline(&mut side.fixed_3prime).desired_width(120.0))
                    .on_hover_text(MOTIF_HELP);
            });
            ui.horizontal(|ui| {
                ui.label("required motifs")
                    .on_hover_text("Motifs that must be present in the annealing segment.");
                ui.add(egui::TextEdit::singleline(&mut side.required_motifs).desired_width(260.0))
                    .on_hover_text(MOTIF_HELP);
            });
            ui.horizontal(|ui| {
                ui.label("forbidden motifs")
                    .on_hover_text("Motifs that must not be present in the annealing segment.");
                ui.add(egui::TextEdit::singleline(&mut side.forbidden_motifs).desired_width(260.0))
                    .on_hover_text(MOTIF_HELP);
            });
            ui.horizontal(|ui| {
                ui.label("locked positions").on_hover_text(
                    "Comma-separated `offset:base` constraints inside annealing segment.",
                );
                ui.add(
                    egui::TextEdit::singleline(&mut side.locked_positions)
                        .desired_width(300.0)
                        .hint_text("offset:base,offset:base"),
                )
                .on_hover_text(
                    "Format: `offset:base` (0-based offsets; IUPAC base token allowed), for example `0:A,3:R`.",
                );
            });
            ui.small(
                "Motif fields are comma-separated IUPAC strings. Locked positions format: offset:base (e.g. 0:A,3:R). Length/Tm/GC/hits use the annealing segment; 5' tail contributes to full-oligo sequence and dimer diagnostics.",
            );
        });
    }

    pub(super) fn render_tm_delta_label(
        ui: &mut egui::Ui,
        prefix: &str,
    ) -> egui::InnerResponse<()> {
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 0.0;
            ui.label(prefix);
            ui.label(egui::RichText::new("m").small());
        })
    }

    pub(super) fn render_primer_pair_constraint_editor(
        ui: &mut egui::Ui,
        id_prefix: &str,
        pair: &mut PrimerPairConstraintUiState,
    ) {
        const AMPLICON_MOTIF_HELP: &str = "Comma-separated motifs using IUPAC DNA codes (for example `ATG,GGNCC,TTTRAA`) evaluated on the amplicon sequence.";
        ui.group(|ui| {
            ui.label("Pair constraints");
            ui.checkbox(&mut pair.require_roi_flanking, "require ROI flanking")
                .on_hover_text(
                    "Require forward/reverse primers to flank the ROI, so the amplicon spans across the selected ROI boundaries.",
                );
            egui::Grid::new(format!("{id_prefix}_grid"))
                .num_columns(4)
                .show(ui, |ui| {
                    ui.label("fixed amplicon start")
                        .on_hover_text(
                            "Optional fixed amplicon start coordinate (0-based). Large genomic coordinates are expected here.",
                        );
                    ui.add(
                        egui::TextEdit::singleline(&mut pair.fixed_amplicon_start_0based)
                            .desired_width(168.0),
                    )
                    .on_hover_text("Fix the amplicon start to this 0-based coordinate.");
                    ui.label("fixed amplicon end (exclusive)").on_hover_text(
                        "Optional fixed amplicon end coordinate (0-based, end-exclusive). Large genomic coordinates are expected here.",
                    );
                    ui.add(
                        egui::TextEdit::singleline(&mut pair.fixed_amplicon_end_0based_exclusive)
                            .desired_width(168.0),
                    )
                    .on_hover_text("Fix the amplicon end to this 0-based end-exclusive coordinate.");
                    ui.end_row();
                });
            ui.horizontal(|ui| {
                ui.label("required amplicon motifs")
                    .on_hover_text("Motifs that must appear in the amplicon sequence.");
                ui.add(
                    egui::TextEdit::singleline(&mut pair.required_amplicon_motifs)
                        .desired_width(300.0),
                )
                .on_hover_text(AMPLICON_MOTIF_HELP);
            });
            ui.horizontal(|ui| {
                ui.label("forbidden amplicon motifs")
                    .on_hover_text("Motifs that must not appear in the amplicon sequence.");
                ui.add(
                    egui::TextEdit::singleline(&mut pair.forbidden_amplicon_motifs)
                        .desired_width(300.0),
                )
                .on_hover_text(AMPLICON_MOTIF_HELP);
            });
        });
    }

    pub(super) fn render_live_pair_pcr_preview_from_paint(&mut self, ui: &mut egui::Ui) {
        ui.group(|ui| {
            ui.label("Live pair-PCR cartoon geometry (paint-first)");
            self.render_pcr_paint_interval_summary(ui);
            let Some((roi_bp, upstream_window_bp, downstream_window_bp)) =
                self.painted_pair_geometry_bp()
            else {
                ui.small(
                    "Paint ROI (green), upstream window (red), and downstream window (blue) on the linear map to activate this live preview.",
                );
                return;
            };
            ui.monospace(format!(
                "roi={} bp | upstream={} bp | downstream={} bp",
                roi_bp, upstream_window_bp, downstream_window_bp
            ));
            let bindings =
                pcr_assay_pair_geometry_bindings(roi_bp, upstream_window_bp, downstream_window_bp);
            ui.small(format!(
                "deterministic bindings ready ({} feature overrides)",
                bindings.feature_overrides.len()
            ));
            match pcr_assay_pair_spec_with_geometry(roi_bp, upstream_window_bp, downstream_window_bp)
            {
                Ok(_spec) => {
                    ui.small("Pair-PCR cartoon spec resolves with current painted geometry.");
                }
                Err(err) => {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!("Could not resolve pair-PCR cartoon: {err}"),
                    );
                }
            }
        });
    }

    pub(super) fn render_primer_design_ops(
        &mut self,
        ui: &mut egui::Ui,
        include_qpcr_section: bool,
    ) {
        let template = self.seq_id.clone().unwrap_or_default();
        if template.trim().is_empty() {
            ui.small("No active sequence selected.");
            return;
        }
        if !include_qpcr_section {
            self.pcr_designer_mode = PcrDesignerMode::PrimerPairs;
        }
        let primer_task_running = self.primer_design_task.is_some();
        ui.small(format!(
            "Template: {} (reports are persisted in project metadata)",
            template
        ));
        ui.group(|ui| {
            ui.label("Protocol cartoon preview (last operation)");
            if let Some(preview) = self.last_protocol_cartoon_preview.as_ref() {
                let op_label = if self.last_protocol_cartoon_preview_op_id.trim().is_empty() {
                    "-".to_string()
                } else {
                    self.last_protocol_cartoon_preview_op_id.clone()
                };
                ui.monospace(format!(
                    "protocol={} | op={} | flank={} bp | overlap={} bp | insert={} bp",
                    preview.protocol, op_label, preview.flank_bp, preview.overlap_bp, preview.insert_bp
                ));
                ui.horizontal(|ui| {
                    if ui
                        .button("Render + Export SVG...")
                        .on_hover_text(
                            "Render this protocol-cartoon preview from bound template geometry and export as SVG",
                        )
                        .clicked()
                    {
                        self.export_protocol_cartoon_preview_svg_dialog();
                    }
                    if ui
                        .button("Clear preview")
                        .on_hover_text("Clear currently cached protocol-cartoon preview payload")
                        .clicked()
                    {
                        self.last_protocol_cartoon_preview = None;
                        self.last_protocol_cartoon_preview_op_id.clear();
                        self.op_status = "Cleared protocol-cartoon preview payload".to_string();
                    }
                });
            } else {
                ui.small(
                    "No preview payload captured yet. OE mutagenesis insertion/replacement operations attach one automatically.",
                );
            }
        });
        ui.group(|ui| {
            ui.label("Primer backend and Primer3 preflight");
            ui.horizontal(|ui| {
                ui.label("backend").on_hover_text(
                    "Choose how primer designs are computed: auto prefers Primer3 when available, internal always uses GENtle's built-in backend, and primer3 requires the external executable.",
                );
                egui::ComboBox::from_id_salt("engine_ops_primer_backend")
                    .selected_text(self.primer_backend.as_str())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.primer_backend,
                            PrimerDesignBackend::Auto,
                            "auto",
                        );
                        ui.selectable_value(
                            &mut self.primer_backend,
                            PrimerDesignBackend::Internal,
                            "internal",
                        );
                        ui.selectable_value(
                            &mut self.primer_backend,
                            PrimerDesignBackend::Primer3,
                            "primer3",
                        );
                    })
                    .response
                    .on_hover_text(
                        "Select auto, internal, or primer3. Reports record both the requested backend and the backend actually used.",
                    );
                ui.label("primer3 executable").on_hover_text(
                    "Path to the Primer3 executable used when backend=primer3 or when auto tries Primer3 first.",
                );
                ui.add(
                    egui::TextEdit::singleline(&mut self.primer3_executable)
                        .desired_width(260.0)
                        .hint_text("primer3_core"),
                )
                .on_hover_text(
                    "Executable path used by primer3 backend (`primer3_core` when empty)",
                );
            });
            ui.horizontal(|ui| {
                if ui
                    .button("Apply Primer Backend")
                    .on_hover_text("Persist backend + primer3 executable into engine parameters")
                    .clicked()
                {
                    self.apply_primer_backend_settings();
                }
                if ui
                    .button("Probe Primer3")
                    .on_hover_text(
                        "Run `primer3 --version` preflight for the configured executable",
                    )
                    .clicked()
                {
                    self.run_primer3_preflight_probe();
                }
                if ui
                    .button("Reload from Engine")
                    .on_hover_text(
                        "Reload backend/executable controls from current engine parameters",
                    )
                    .clicked()
                {
                    self.sync_primer_backend_controls_from_engine();
                }
            });
            ui.small(self.primer_backend_help_summary());
            if !self.primer3_preflight_status.trim().is_empty() {
                ui.monospace(&self.primer3_preflight_status).on_hover_text(
                    "Latest explicit Primer3 probe result. When backend=auto and Primer3 is unavailable, GENtle falls back to the internal backend and records that fallback in saved reports.",
                );
            }
        });
        if include_qpcr_section {
            ui.group(|ui| {
                ui.label("PCR Designer mode");
                ui.horizontal(|ui| {
                    for mode in [PcrDesignerMode::PrimerPairs, PcrDesignerMode::QpcrAssays] {
                        let selected = self.pcr_designer_mode == mode;
                        if ui
                            .selectable_label(selected, mode.label())
                            .on_hover_text(mode.helper_text())
                            .clicked()
                        {
                            self.pcr_designer_mode = mode;
                            self.save_engine_ops_state();
                        }
                    }
                });
                ui.small(self.pcr_designer_mode.helper_text());
            });
        }

        if self.pcr_designer_mode == PcrDesignerMode::PrimerPairs {
            egui::CollapsingHeader::new("Design primer pairs")
                .default_open(true)
                .show(ui, |ui| {
                egui::Grid::new("primer_pairs_overview_grid")
                    .num_columns(4)
                    .spacing([12.0, 6.0])
                    .show(ui, |ui| {
                        ui.label("ROI start").on_hover_text(
                            "PCR region-of-interest start coordinate (0-based). Supports formulas: `=CDS.start+10` or range form `=CDS.start+10 .. CDS.end-500`.",
                        );
                        ui.add(
                            egui::TextEdit::singleline(&mut self.primer_design_ui.roi_start_0based)
                                .desired_width(168.0),
                        )
                        .on_hover_text(
                            "0-based ROI start. Accepts `=feature.boundary(+/-offset)` formulas. If range form is used here (`=left .. right`), ROI end field is ignored for this resolve.",
                        );
                        ui.label("ROI end").on_hover_text(
                            "PCR region-of-interest end coordinate (0-based, end-exclusive). Supports formulas: `=CDS.end-500`.",
                        );
                        ui.add(
                            egui::TextEdit::singleline(&mut self.primer_design_ui.roi_end_0based)
                                .desired_width(168.0),
                        )
                        .on_hover_text(
                            "0-based end-exclusive ROI end. Supports `=feature.boundary(+/-offset)` formulas.",
                        );
                        ui.end_row();
                        ui.label("min amplicon")
                            .on_hover_text("Minimum amplicon length in bp.");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.primer_design_ui.min_amplicon_bp)
                                .desired_width(92.0),
                        )
                        .on_hover_text("Lower amplicon length bound.");
                        ui.label("max amplicon")
                            .on_hover_text("Maximum amplicon length in bp.");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.primer_design_ui.max_amplicon_bp)
                                .desired_width(92.0),
                        )
                        .on_hover_text("Upper amplicon length bound.");
                        ui.end_row();
                        Self::render_tm_delta_label(ui, "max ΔT")
                            .response
                            .on_hover_text("Maximum allowed ΔTₘ (°C) between primer sides.");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.primer_design_ui.max_tm_delta_c)
                                .desired_width(92.0),
                        )
                        .on_hover_text("Optional ΔTₘ constraint for primer-pair ranking.");
                        ui.label("max pairs")
                            .on_hover_text("Maximum number of accepted primer pairs to return.");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.primer_design_ui.max_pairs)
                                .desired_width(92.0),
                        )
                        .on_hover_text("Upper limit for returned primer-pair candidates.");
                        ui.end_row();
                        ui.label("report_id").on_hover_text(
                            "Persisted report identifier. Batch mode derives deterministic suffixes (`_r01`, `_r02`, ...).",
                        );
                        ui.add(
                            egui::TextEdit::singleline(&mut self.primer_design_ui.report_id)
                                .desired_width(240.0),
                        )
                        .on_hover_text(
                            "Optional report id stem. Empty value auto-derives from template and ROI.",
                        );
                        ui.end_row();
                    });
                ui.horizontal_wrapped(|ui| {
                    if ui
                        .button("Apply ROI formula")
                        .on_hover_text(
                            "Resolve current ROI fields (supports `=` formulas) into numeric 0-based coordinates",
                        )
                        .clicked()
                    {
                        self.resolve_and_apply_primer_roi_fields();
                    }
                    ui.small(
                        "Formula syntax: `=KIND.start+N`, `=KIND.end-N`, optional occurrence `KIND[2]`, optional label filter `KIND[label=TP73]`, and range form `=left .. right` (or `=left to right`).",
                    );
                });
                ui.group(|ui| {
                    ui.label("Simple PCR starter");
                    ui.small(
                        "For a simple PCR, think in three inputs: core ROI, maximum primer distance from that core, and maximum amplicon length.",
                    );
                    ui.horizontal_wrapped(|ui| {
                        let selection_roi = self.current_selection_range_0based();
                        let seed_simple_response = ui.add_enabled(
                            selection_roi.is_some(),
                            egui::Button::new("Use current selection as simple PCR core"),
                        );
                        let seed_simple_response = if selection_roi.is_some() {
                            seed_simple_response.on_hover_text(
                                "Seed ROI from current selection, require ROI flanking, set min amplicon to core length, apply flank windows, and open PCR Designer if needed",
                            )
                        } else {
                            seed_simple_response.on_hover_text(
                                "Requires a non-empty current selection on the sequence map",
                            )
                        };
                        if seed_simple_response.clicked() {
                            self.open_simple_pcr_designer_from_current_selection();
                        }
                        ui.label("max primer distance from core").on_hover_text(
                            "Maximum allowed flanking search width on each side of the core ROI. This maps onto forward/reverse side start/end windows.",
                        );
                        ui.add(
                            egui::TextEdit::singleline(
                                &mut self.simple_pcr_max_primer_distance_bp,
                            )
                            .desired_width(92.0),
                        )
                        .on_hover_text(
                            "Width of the upstream and downstream primer-search windows measured from the ROI edges.",
                        );
                        if ui
                            .button("Apply simple flank windows")
                            .on_hover_text(
                                "Set forward side window to [ROI start - D .. ROI start] and reverse side window to [ROI end .. ROI end + D], while keeping ROI flanking enabled",
                            )
                            .clicked()
                        {
                            if let Err(message) = self.apply_simple_pcr_flank_windows_from_roi() {
                                self.op_status = message;
                            }
                        }
                    });
                    ui.small(
                        "Then set `max amplicon` below to the longest acceptable product. The helper keeps `min amplicon` at least the core length.",
                    );
                });
                ui.group(|ui| {
                    ui.label("PCR region queue (batch source)");
                    let copy_toggle = ui.checkbox(
                        &mut self.pcr_batch_create_extract_copies,
                        "Also create extracted region copies",
                    )
                    .on_hover_text(
                        "Create one extracted sequence artifact per queued ROI before/alongside primer design.",
                    );
                    if copy_toggle.changed() {
                        self.save_engine_ops_state();
                    }
                    ui.small(
                        "Queue entries are genomic/template ROI specs (`template + start + end`), not Primer3 jobs.",
                    );
                    if self.pcr_queued_regions_ui.is_empty() {
                        ui.small(
                            "Queue regions from DNA-window `PCR ROI` actions or `Queue current ROI spec`, then run batch primer design.",
                        );
                    } else {
                        let mut use_roi_idx: Option<usize> = None;
                        let mut remove_idx: Option<usize> = None;
                        egui::Grid::new("primer_pairs_region_queue_grid")
                            .striped(true)
                            .num_columns(8)
                            .show(ui, |ui| {
                                ui.strong("#");
                                ui.strong("source");
                                ui.strong("start");
                                ui.strong("end");
                                ui.strong("len");
                                ui.strong("template");
                                ui.strong("use");
                                ui.strong("remove");
                                ui.end_row();
                                for (idx, row) in self.pcr_queued_regions_ui.iter().enumerate() {
                                    ui.monospace(format!("{:02}", idx + 1));
                                    ui.label(&row.source_label);
                                    ui.monospace(format!("{}", row.start_0based));
                                    ui.monospace(format!("{}", row.end_0based_exclusive));
                                    ui.monospace(format!("{}", row.span_len_bp()));
                                    ui.monospace(&row.template);
                                    if ui
                                        .small_button("Use ROI")
                                        .on_hover_text(
                                            "Copy this queued region back into the active PCR ROI start/end fields above",
                                        )
                                        .clicked()
                                    {
                                        use_roi_idx = Some(idx);
                                    }
                                    if ui
                                        .small_button("Remove")
                                        .on_hover_text(
                                            "Remove this queued PCR region from the batch source table",
                                        )
                                        .clicked()
                                    {
                                        remove_idx = Some(idx);
                                    }
                                    ui.end_row();
                                }
                            });
                        if let Some(idx) = use_roi_idx {
                            if let Some(row) = self.pcr_queued_regions_ui.get(idx).cloned() {
                                self.set_primer_design_roi_fields_0based(
                                    row.start_0based,
                                    row.end_0based_exclusive,
                                );
                                self.show_engine_ops = true;
                                self.op_status = format!(
                                    "Set active PCR ROI from queued region #{:02}: {}..{}",
                                    idx + 1,
                                    row.start_0based,
                                    row.end_0based_exclusive
                                );
                                self.save_engine_ops_state();
                            }
                        }
                        if let Some(idx) = remove_idx {
                            self.remove_pcr_queued_region(idx);
                        }
                    }
                    ui.horizontal(|ui| {
                        if ui
                            .button("Queue current ROI spec")
                            .on_hover_text(
                                "Queue the ROI start/end fields above as one PCR region spec for batch processing",
                            )
                            .clicked()
                        {
                            self.queue_current_primer_roi_fields_for_pcr();
                        }
                        let run_batch_response = ui.add_enabled(
                            !self.pcr_queued_regions_ui.is_empty() && !primer_task_running,
                            egui::Button::new("Design Primer Pairs for queued regions"),
                        );
                        let run_batch_response = run_batch_response.on_hover_text(
                            "Run one DesignPrimerPairs operation per queued region using current primer constraints",
                        );
                        if run_batch_response.clicked() {
                            self.start_queued_primer_pair_design_batch();
                        }
                        if ui
                            .button("Clear queue")
                            .on_hover_text("Remove all queued PCR regions")
                            .clicked()
                        {
                            self.clear_pcr_region_queue();
                        }
                    });
                });
                ui.group(|ui| {
                    ui.label("PCR batch results");
                    if self.pcr_batch_results_ui.is_empty() {
                        ui.small(
                            "No batch run recorded yet. Run `Design Primer Pairs for queued regions` to populate this table.",
                        );
                    } else {
                        let mut show_report_id: Option<String> = None;
                        let mut export_report_id: Option<String> = None;
                        let mut open_seq_id: Option<String> = None;
                        egui::Grid::new("primer_pairs_batch_results_grid")
                            .striped(true)
                            .num_columns(11)
                            .show(ui, |ui| {
                                ui.strong("#");
                                ui.strong("status");
                                ui.strong("source");
                                ui.strong("region");
                                ui.strong("len");
                                ui.strong("report_id");
                                ui.strong("copy_id");
                                ui.strong("Show");
                                ui.strong("Export");
                                ui.strong("Open sequence");
                                ui.strong("detail");
                                ui.end_row();
                                for row in &self.pcr_batch_results_ui {
                                    ui.monospace(format!("{:02}", row.region_index_1based));
                                    ui.monospace(format!(
                                        "{} | {}",
                                        row.report_status_label(),
                                        row.copy_status_label()
                                    ));
                                    ui.label(&row.source_label);
                                    ui.monospace(format!(
                                        "{}:{}..{}",
                                        row.template, row.start_0based, row.end_0based_exclusive
                                    ));
                                    ui.monospace(format!("{}", row.span_len_bp()));
                                    ui.monospace(&row.report_id);
                                    if let Some(copy_seq_id) = &row.copy_seq_id {
                                        ui.monospace(copy_seq_id);
                                    } else if row.copy_requested {
                                        ui.label("-");
                                    } else {
                                        ui.label("n/a");
                                    }
                                    let report_ok = row.report_error.is_none();
                                    if ui
                                        .add_enabled(report_ok, egui::Button::new("Show"))
                                        .on_hover_text("Show summary for this report_id")
                                        .clicked()
                                    {
                                        show_report_id = Some(row.report_id.clone());
                                    }
                                    if ui
                                        .add_enabled(report_ok, egui::Button::new("Export"))
                                        .on_hover_text("Export this report_id to JSON")
                                        .clicked()
                                    {
                                        export_report_id = Some(row.report_id.clone());
                                    }
                                    let open_target = row
                                        .copy_seq_id
                                        .clone()
                                        .unwrap_or_else(|| row.template.clone());
                                    if ui
                                        .button("Open")
                                        .on_hover_text(
                                            "Open extracted copy when available; otherwise open template sequence",
                                        )
                                        .clicked()
                                    {
                                        open_seq_id = Some(open_target);
                                    }
                                    let detail = row
                                        .report_error
                                        .clone()
                                        .or_else(|| row.report_note.clone())
                                        .or_else(|| row.copy_error.clone())
                                        .unwrap_or_else(|| {
                                            row.report_pair_count
                                                .map(|count| format!("pairs={count}"))
                                                .unwrap_or_else(|| "ok".to_string())
                                        });
                                    ui.label(detail);
                                    ui.end_row();
                                }
                            });
                        ui.horizontal(|ui| {
                            if ui
                                .button("Clear results")
                                .on_hover_text("Remove all rows from PCR batch results table")
                                .clicked()
                            {
                                self.clear_pcr_batch_results();
                            }
                        });
                        if let Some(report_id) = show_report_id {
                            self.primer_design_ui.report_id = report_id.clone();
                            self.show_primer_design_report(&report_id);
                        }
                        if let Some(report_id) = export_report_id {
                            self.primer_design_ui.report_id = report_id.clone();
                            self.export_primer_design_report_dialog(&report_id);
                        }
                        if let Some(seq_id) = open_seq_id
                            && let Err(err) = self.open_sequence_by_id(&seq_id)
                        {
                            self.op_status = err;
                        }
                    }
                });
                Self::render_primer_pair_constraint_editor(
                    ui,
                    "primer_pairs_pair_constraints",
                    &mut self.primer_design_ui.pair_constraints,
                );
                self.render_live_pair_pcr_preview_from_paint(ui);
                Self::render_primer_side_constraint_editor(
                    ui,
                    "primer_pairs_forward",
                    "Forward side",
                    &mut self.primer_design_ui.forward,
                );
                Self::render_primer_side_constraint_editor(
                    ui,
                    "primer_pairs_reverse",
                    "Reverse side",
                    &mut self.primer_design_ui.reverse,
                );
                let design_button = if primer_task_running {
                    "Design Primer Pairs (running...)"
                } else {
                    "Design Primer Pairs"
                };
                if ui
                    .add_enabled(!primer_task_running, egui::Button::new(design_button))
                    .on_hover_text("Run DesignPrimerPairs with current GUI constraints")
                    .clicked()
                {
                    match self.build_design_primer_pairs_operation(&template) {
                        Ok(op) => self.start_primer_design_operation(op, "Primer-pair design"),
                        Err(err) => self.op_status = err,
                    }
                }
                if let Some(summary) = self
                    .active_primer_design_progress_summary_for_kind("primer_pairs")
                {
                    ui.small(format!("Progress: {summary}"));
                }
                ui.horizontal(|ui| {
                    if ui
                        .button("List Primer Reports")
                        .on_hover_text("List persisted primer-design report IDs")
                        .clicked()
                    {
                        self.list_primer_design_reports();
                    }
                    if ui
                        .button("Show report_id")
                        .on_hover_text("Show summary for report_id from this form")
                        .clicked()
                    {
                        let report_id = self.primer_design_ui.report_id.clone();
                        self.show_primer_design_report(&report_id);
                    }
                    if ui
                        .button("Export report_id...")
                        .on_hover_text("Export report_id from this form to JSON")
                        .clicked()
                    {
                        let report_id = self.primer_design_ui.report_id.clone();
                        self.export_primer_design_report_dialog(&report_id);
                    }
                    if ui
                        .button("Export Primer3 input...")
                        .on_hover_text(
                            "Export stored Primer3 Boulder-IO request payload for this report_id (when backend used Primer3)",
                        )
                        .clicked()
                    {
                        let report_id = self.primer_design_ui.report_id.clone();
                        self.export_primer3_request_dialog(&report_id);
                    }
                });
                ui.group(|ui| {
                    ui.label("Local specificity confirmation");
                    ui.small(
                        "Confirm one saved primer-pair rank against a prepared reference genome with local BLAST. GENtle uses only annealing segments for BLAST and keeps 5' tails as provenance.",
                    );
                    ui.horizontal_wrapped(|ui| {
                        ui.label("target genome").on_hover_text(
                            "Prepared reference genome id, for example `GRCh38.p14`.",
                        );
                        ui.add(
                            egui::TextEdit::singleline(
                                &mut self.primer_design_ui.specificity_target_genome_id,
                            )
                            .desired_width(160.0),
                        )
                        .on_hover_text("Genome must already be prepared with a BLAST index.");
                        ui.label("pair rank").on_hover_text(
                            "Saved primer-pair rank from the report, 1-based.",
                        );
                        ui.add(
                            egui::TextEdit::singleline(
                                &mut self.primer_design_ui.specificity_pair_rank_1based,
                            )
                            .desired_width(52.0),
                        );
                        ui.label("max product").on_hover_text(
                            "Maximum genomic product length retained as a candidate amplicon.",
                        );
                        ui.add(
                            egui::TextEdit::singleline(
                                &mut self.primer_design_ui.specificity_max_target_amplicon_bp,
                            )
                            .desired_width(72.0),
                        );
                        ui.label("max hits").on_hover_text(
                            "Maximum BLAST hits retained per primer query.",
                        );
                        ui.add(
                            egui::TextEdit::singleline(
                                &mut self.primer_design_ui.specificity_max_hits_per_primer,
                            )
                            .desired_width(62.0),
                        );
                        if ui
                            .button("Confirm specificity")
                            .on_hover_text(
                                "Run AssessPrimerPairSpecificity for report_id + pair rank against the prepared genome.",
                            )
                            .clicked()
                        {
                            let report_id = self.primer_design_ui.report_id.clone();
                            self.confirm_primer_specificity_for_report(&report_id);
                        }
                    });
                });
                self.render_primer_design_report_preview(ui);
                self.render_restriction_cloning_handoff_section(ui, &template);
            });
        }

        if include_qpcr_section && self.pcr_designer_mode == PcrDesignerMode::QpcrAssays {
            egui::CollapsingHeader::new("Design qPCR assays")
                .default_open(true)
                .show(ui, |ui| {
                ui.group(|ui| {
                    ui.label("qPCR protocol cartoon");
                    ui.small(
                        "Use the built-in probe-bearing qPCR strip to explain or promote the assay workflow from the same deterministic protocol-cartoon family.",
                    );
                    ui.horizontal_wrapped(|ui| {
                        if ui
                            .button("Export qPCR Cartoon SVG...")
                            .on_hover_text(
                                "Export the built-in qPCR assay protocol cartoon (`pcr.assay.qpcr`) as SVG.",
                            )
                            .clicked()
                        {
                            self.export_builtin_protocol_cartoon_svg_dialog(
                                ProtocolCartoonKind::PcrAssayQpcr,
                                "qpcr_assay_protocol_cartoon",
                            );
                        }
                        ui.monospace("protocol-cartoon render-svg pcr.assay.qpcr OUTPUT.svg");
                    });
                });
                ui.small(
                    "qPCR builds on the same template/ROI context as pair-PCR, then adds probe-side filtering and assay ranking. Pair-PCR queue/batch tools remain in the Pair PCR mode for now.",
                );
                let current_splicing_view =
                    self.relevant_qpcr_splicing_view_for_template(&template);
                let current_context_summary = self.qpcr_transcript_context_summary(&template);
                let transcript_targeting_error =
                    self.validate_qpcr_transcript_targeting_ui(&template).err();
                ui.group(|ui| {
                    ui.label("Transcript targeting");
                    let mut intent_changed = false;
                    ui.horizontal_wrapped(|ui| {
                        intent_changed |= ui
                            .radio_value(
                                &mut self.qpcr_design_ui.transcript_targeting.assay_intent,
                                QpcrTranscriptIntentUiMode::Genomic,
                                QpcrTranscriptIntentUiMode::Genomic.label(),
                            )
                            .changed();
                        intent_changed |= ui
                            .radio_value(
                                &mut self.qpcr_design_ui.transcript_targeting.assay_intent,
                                QpcrTranscriptIntentUiMode::SharedAcrossTranscripts,
                                QpcrTranscriptIntentUiMode::SharedAcrossTranscripts.label(),
                            )
                            .changed();
                        intent_changed |= ui
                            .radio_value(
                                &mut self.qpcr_design_ui.transcript_targeting.assay_intent,
                                QpcrTranscriptIntentUiMode::SpecificTranscript,
                                QpcrTranscriptIntentUiMode::SpecificTranscript.label(),
                            )
                            .changed();
                    });
                    if intent_changed {
                        self.save_engine_ops_state();
                    }
                    let assay_intent = self.qpcr_design_ui.transcript_targeting.assay_intent;
                    ui.small(format!(
                        "Assay intent: {}.",
                        assay_intent.summary_label()
                    ));
                    match current_context_summary.as_deref() {
                        Some(summary) => ui.small(format!("Transcript context: {summary}")),
                        None if assay_intent == QpcrTranscriptIntentUiMode::Genomic => ui.small(
                            "Transcript context: optional. Genomic qPCR keeps the existing direct PCR Designer flow.",
                        ),
                        None => ui.small(
                            "Transcript context: none attached yet for this sequence. Seed from a splicing overview before running transcript-aware qPCR.",
                        ),
                    };
                    if assay_intent != QpcrTranscriptIntentUiMode::Genomic {
                        ui.horizontal_wrapped(|ui| {
                            if let Some(view) = current_splicing_view.as_ref() {
                                let button_label = if current_context_summary.is_some() {
                                    "Refresh from current Splicing Expert context"
                                } else {
                                    "Use current Splicing Expert context"
                                };
                                if ui
                                    .button(button_label)
                                    .on_hover_text(
                                        "Attach the current Splicing Expert group and transcript selection to this qPCR request",
                                    )
                                    .clicked()
                                {
                                    match self.apply_qpcr_transcript_targeting_context_from_view(
                                        view,
                                        assay_intent,
                                    ) {
                                        Ok(()) => {
                                            self.op_status = format!(
                                                "Attached transcript-aware qPCR context from splicing group '{}'",
                                                view.group_label
                                            );
                                            self.save_engine_ops_state();
                                        }
                                        Err(err) => self.op_status = err,
                                    }
                                }
                            }
                            if ui
                                .button("Open / Focus Splicing Expert")
                                .on_hover_text(
                                    "Open or focus Splicing Expert so transcript-aware qPCR can reuse its group and transcript selection",
                                )
                                .clicked()
                            {
                                self.focus_or_open_qpcr_splicing_expert(ui.ctx(), &template);
                            }
                        });
                    }
                    match assay_intent {
                        QpcrTranscriptIntentUiMode::Genomic => {}
                        QpcrTranscriptIntentUiMode::SharedAcrossTranscripts => {
                            ui.small(
                                "GENtle will prefer exon or exon-chain context shared across the whole group. If no fully shared span exists, it falls back to the largest supported transcript subset and reports that explicitly.",
                            );
                        }
                        QpcrTranscriptIntentUiMode::SpecificTranscript => {
                            let transcript_label = if self
                                .qpcr_design_ui
                                .transcript_targeting
                                .selected_transcript_id
                                .trim()
                                .is_empty()
                            {
                                "<select in Splicing Expert>".to_string()
                            } else {
                                let transcript_id = self
                                    .qpcr_design_ui
                                    .transcript_targeting
                                    .selected_transcript_id
                                    .trim()
                                    .to_string();
                                let transcript_label = self
                                    .qpcr_design_ui
                                    .transcript_targeting
                                    .selected_transcript_label
                                    .trim()
                                    .to_string();
                                if transcript_label.is_empty()
                                    || transcript_label.eq_ignore_ascii_case(&transcript_id)
                                {
                                    transcript_id
                                } else {
                                    format!("{transcript_id} ({transcript_label})")
                                }
                            };
                            ui.small(format!("Selected transcript: {transcript_label}"));
                            ui.horizontal_wrapped(|ui| {
                                ui.label("Specificity evidence").on_hover_text(
                                    "Require a junction-spanning primer, a transcript-unique exon/exon-chain primer, or let GENtle prefer junction evidence when both are possible.",
                                );
                                egui::ComboBox::from_id_salt("qpcr_transcript_specificity_evidence")
                                    .selected_text(Self::qpcr_specificity_evidence_label(
                                        self.qpcr_design_ui
                                            .transcript_targeting
                                            .specificity_evidence,
                                    ))
                                    .show_ui(ui, |ui| {
                                        for evidence in [
                                            QpcrTranscriptSpecificityEvidence::JunctionOnly,
                                            QpcrTranscriptSpecificityEvidence::UniqueExonOrChain,
                                            QpcrTranscriptSpecificityEvidence::EitherPreferJunction,
                                        ] {
                                            ui.selectable_value(
                                                &mut self
                                                    .qpcr_design_ui
                                                    .transcript_targeting
                                                    .specificity_evidence,
                                                evidence,
                                                Self::qpcr_specificity_evidence_label(evidence),
                                            );
                                        }
                                    });
                            });
                        }
                    }
                    if let Some(err) = transcript_targeting_error.as_deref() {
                        ui.colored_label(
                            egui::Color32::from_rgb(180, 83, 9),
                            err,
                        );
                    }
                });
                    ui.horizontal(|ui| {
                        ui.label("ROI start").on_hover_text(
                            "PCR region-of-interest start coordinate (0-based). Supports formulas: `=CDS.start+10` or range form `=CDS.start+10 .. CDS.end-500`.",
                        );
                        ui.add(
                            egui::TextEdit::singleline(&mut self.qpcr_design_ui.roi_start_0based)
                            .desired_width(168.0),
                        )
                        .on_hover_text(
                            "0-based ROI start. Accepts `=feature.boundary(+/-offset)` formulas.",
                        );
                        ui.label("ROI end").on_hover_text(
                            "PCR region-of-interest end coordinate (0-based, end-exclusive). Supports formulas: `=CDS.end-500`.",
                        );
                        ui.add(
                            egui::TextEdit::singleline(&mut self.qpcr_design_ui.roi_end_0based)
                            .desired_width(168.0),
                        )
                        .on_hover_text(
                            "0-based end-exclusive ROI end. Supports `=feature.boundary(+/-offset)` formulas.",
                        );
                        ui.label("min amplicon")
                            .on_hover_text("Minimum amplicon length in bp.");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.qpcr_design_ui.min_amplicon_bp)
                            .desired_width(92.0),
                    )
                    .on_hover_text("Lower amplicon length bound.");
                    ui.label("max amplicon")
                        .on_hover_text("Maximum amplicon length in bp.");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.qpcr_design_ui.max_amplicon_bp)
                            .desired_width(92.0),
                        )
                        .on_hover_text("Upper amplicon length bound.");
                    });
                    ui.horizontal_wrapped(|ui| {
                        if ui
                            .button("Apply ROI formula")
                            .on_hover_text(
                                "Resolve current qPCR ROI fields (supports `=` formulas) into numeric 0-based coordinates",
                            )
                            .clicked()
                        {
                            self.resolve_and_apply_qpcr_roi_fields();
                        }
                        ui.small(
                            "Formula syntax mirrors primer-pair ROI: `=KIND.start+N`, optional `KIND[2]` or `KIND[label=TP73]`, range form `=left .. right`.",
                        );
                    });
                    ui.horizontal(|ui| {
                    Self::render_tm_delta_label(ui, "max primer ΔT")
                        .response
                        .on_hover_text("Maximum allowed ΔTₘ (°C) between qPCR primer sides.");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.qpcr_design_ui.max_tm_delta_c)
                            .desired_width(92.0),
                    )
                    .on_hover_text("Optional ΔTₘ constraint for qPCR primer-side ranking.");
                    Self::render_tm_delta_label(ui, "max probe ΔT")
                        .response
                        .on_hover_text("Maximum allowed ΔTₘ (°C) between probe-side annealing temperatures.");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.qpcr_design_ui.max_probe_tm_delta_c)
                            .desired_width(92.0),
                    )
                    .on_hover_text("Optional ΔTₘ constraint for probe-side ranking.");
                    ui.label("max assays")
                        .on_hover_text("Maximum number of accepted qPCR assays to return.");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.qpcr_design_ui.max_assays)
                            .desired_width(80.0),
                    )
                    .on_hover_text("Upper limit for returned qPCR assay candidates.");
                    ui.label("report_id").on_hover_text(
                        "Persisted report identifier for qPCR. Batch mode is not used here, so no `_rNN` suffixing applies automatically.",
                    );
                    ui.add(
                        egui::TextEdit::singleline(&mut self.qpcr_design_ui.report_id)
                            .desired_width(180.0),
                    )
                    .on_hover_text(
                        "Optional report id stem. Empty value auto-derives from template and ROI.",
                    );
                });
                Self::render_primer_pair_constraint_editor(
                    ui,
                    "qpcr_pair_constraints",
                    &mut self.qpcr_design_ui.pair_constraints,
                );
                Self::render_primer_side_constraint_editor(
                    ui,
                    "qpcr_forward",
                    "Forward side",
                    &mut self.qpcr_design_ui.forward,
                );
                Self::render_primer_side_constraint_editor(
                    ui,
                    "qpcr_reverse",
                    "Reverse side",
                    &mut self.qpcr_design_ui.reverse,
                );
                Self::render_primer_side_constraint_editor(
                    ui,
                    "qpcr_probe",
                    "Probe side",
                    &mut self.qpcr_design_ui.probe,
                );
                self.render_live_qpcr_cartoon_preview(ui);
                let qpcr_button = if primer_task_running {
                    "Design qPCR Assays (running...)"
                } else {
                    "Design qPCR Assays"
                };
                if ui
                    .add_enabled(
                        !primer_task_running && transcript_targeting_error.is_none(),
                        egui::Button::new(qpcr_button),
                    )
                    .on_hover_text("Run DesignQpcrAssays with current GUI constraints")
                    .clicked()
                {
                    match self.build_design_qpcr_operation(&template) {
                        Ok(op) => self.start_primer_design_operation(op, "qPCR design"),
                        Err(err) => self.op_status = err,
                    }
                }
                if let Some(summary) = self
                    .active_primer_design_progress_summary_for_kind("qpcr_assays")
                {
                    ui.small(format!("Progress: {summary}"));
                }
                ui.horizontal(|ui| {
                    if ui
                        .button("List qPCR Reports")
                        .on_hover_text("List persisted qPCR report IDs")
                        .clicked()
                    {
                        self.list_qpcr_design_reports();
                    }
                    if ui
                        .button("Show report_id")
                        .on_hover_text("Show summary for qPCR report_id from this form")
                        .clicked()
                    {
                        let report_id = self.qpcr_design_ui.report_id.clone();
                        self.show_qpcr_design_report(&report_id);
                    }
                    if ui
                        .button("Export report_id...")
                        .on_hover_text("Export qPCR report_id from this form to JSON")
                        .clicked()
                    {
                        let report_id = self.qpcr_design_ui.report_id.clone();
                        self.export_qpcr_design_report_dialog(&report_id);
                    }
                });
                self.render_qpcr_design_report_preview(ui);
                });
        }
    }

    pub(super) fn render_pcr_queue_compact_table(&mut self, ui: &mut egui::Ui) {
        ui.label("PCR region queue");
        if self.pcr_queued_regions_ui.is_empty() {
            ui.small("Queue is empty.");
            return;
        }
        let mut remove_idx: Option<usize> = None;
        egui::Grid::new("pcr_designer_queue_compact")
            .striped(true)
            .num_columns(6)
            .show(ui, |ui| {
                ui.strong("#");
                ui.strong("source");
                ui.strong("start");
                ui.strong("end");
                ui.strong("len");
                ui.strong("remove");
                ui.end_row();
                for (idx, row) in self.pcr_queued_regions_ui.iter().enumerate() {
                    ui.monospace(format!("{:02}", idx + 1));
                    ui.label(&row.source_label);
                    ui.monospace(format!("{}", row.start_0based));
                    ui.monospace(format!("{}", row.end_0based_exclusive));
                    ui.monospace(format!("{}", row.span_len_bp()));
                    if ui.small_button("x").clicked() {
                        remove_idx = Some(idx);
                    }
                    ui.end_row();
                }
            });
        if let Some(idx) = remove_idx {
            self.remove_pcr_queued_region(idx);
        }
    }

    pub(super) fn render_pcr_designer_map_panel(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        let width = ui.available_width().max(320.0);
        let response = ui.add_sized(egui::Vec2::new(width, 220.0), self.map_dna.to_owned());
        if response.rect.width().is_finite() && response.rect.width() > 0.0 {
            self.last_linear_map_width_px = response.rect.width();
        }
        if !self.is_circular() {
            self.handle_linear_map_drag_for_pcr_paint(&response, ctx);
            if self.linear_pan_drag_origin_bp.is_some() {
                ui.output_mut(|o| o.cursor_icon = egui::CursorIcon::Grabbing);
            }
        } else {
            ui.colored_label(
                egui::Color32::from_rgb(170, 120, 50),
                "Paint-first PCR interaction is available in linear map mode. Use existing ROI fallback actions for circular map mode.",
            );
        }
        self.draw_pcr_paint_overlays(ui, &response);
        self.render_pcr_post_drag_actions(ui);
        response.context_menu(|ui| {
            if self.render_selection_simple_pcr_context_action(ui) {
                ui.close();
                return;
            }
            ui.label("No selection-specific PCR action available");
        });
    }

    pub fn render_pcr_designer_specialist(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        ui.label(
            "Paint-first PCR/qPCR designer. Drag on the linear DNA map with fixed semantic colors: ROI (green), upstream window (red), downstream window (blue).",
        );
        self.render_pcr_paint_role_controls(ui, true);
        ui.small(
            "Mouse shortcuts: drag paints selected role; Shift+drag on ROI queues immediately; Option/Alt+drag pans.",
        );
        ui.separator();
        ui.columns(2, |columns| {
            columns[0].heading(match self.pcr_designer_mode {
                PcrDesignerMode::PrimerPairs => "Paint + Queue",
                PcrDesignerMode::QpcrAssays => "Paint + ROI",
            });
            self.render_selection_formula_inline_controls(&mut columns[0], 280.0);
            let selection_roi = self.current_selection_range_0based();
            columns[0].horizontal_wrapped(|ui| {
                let set_clicked = ui
                    .add_enabled(
                        selection_roi.is_some(),
                        egui::Button::new("Set ROI from selection"),
                    )
                    .on_hover_text("Copy current map/text selection into pair-PCR ROI form fields")
                    .clicked();
                if set_clicked {
                    if let Err(err) = self
                        .set_primer_design_roi_from_current_selection("current sequence selection")
                    {
                        self.op_status = err;
                    }
                }
                if self.pcr_designer_mode == PcrDesignerMode::PrimerPairs {
                    let queue_clicked = ui
                        .add_enabled(
                            selection_roi.is_some(),
                            egui::Button::new("Queue selection"),
                        )
                        .on_hover_text("Queue current map/text selection as one PCR region")
                        .clicked();
                    if queue_clicked {
                        self.queue_current_selection_for_pcr();
                    }
                }
            });
            if let Some((start, end_exclusive)) = selection_roi {
                columns[0].small(format!(
                    "Active selection: {start}..{end_exclusive} (len {} bp)",
                    end_exclusive.saturating_sub(start)
                ));
            } else {
                columns[0].small("No active selection. Use formula or drag-select on map.");
            }
            if self.pcr_designer_mode == PcrDesignerMode::QpcrAssays {
                columns[0].small(
                    "qPCR mode reuses the same ROI painting and selection workflow, but queued batch-region actions remain pair-PCR-only in this first pass.",
                );
            }
            columns[0].separator();
            self.render_pcr_designer_map_panel(&mut columns[0], ctx);
            columns[0].separator();
            self.render_pcr_paint_interval_summary(&mut columns[0]);
            columns[0].horizontal(|ui| {
                if self.pcr_designer_mode == PcrDesignerMode::PrimerPairs {
                    if ui
                        .button("Queue painted ROI")
                        .on_hover_text("Queue the currently painted ROI interval")
                        .clicked()
                    {
                        if let Some((start, end_exclusive)) = self.pcr_paint_intervals.roi {
                            self.queue_painted_roi_from_interval(
                                start,
                                end_exclusive,
                                "painted ROI",
                            );
                        } else {
                            self.op_status =
                                "Paint an ROI interval first (green role) before queueing"
                                    .to_string();
                        }
                    }
                }
                if ui
                    .button("Set form ROI from paint")
                    .on_hover_text("Copy painted ROI interval into the active PCR/qPCR ROI form fields")
                    .clicked()
                {
                    if let Some((start, end_exclusive)) = self.pcr_paint_intervals.roi {
                        self.set_primer_design_roi_fields_0based(start, end_exclusive);
                        self.op_status = format!(
                            "Set PCR Designer ROI fields from painted interval: {}..{}",
                            start, end_exclusive
                        );
                        self.save_engine_ops_state();
                    } else {
                        self.op_status =
                            "Paint an ROI interval first (green role) before setting form ROI"
                                .to_string();
                    }
                }
            });
            if self.pcr_designer_mode == PcrDesignerMode::PrimerPairs {
                self.render_pcr_queue_compact_table(&mut columns[0]);
            }

            columns[1].heading(match self.pcr_designer_mode {
                PcrDesignerMode::PrimerPairs => "Pair-PCR Constraints + Run",
                PcrDesignerMode::QpcrAssays => "qPCR Constraints + Run",
            });
            egui::ScrollArea::vertical()
                .id_salt("pcr_designer_pair_scroll")
                .show(&mut columns[1], |ui| {
                    self.render_primer_design_ops(ui, true);
                });
        });
    }
}
