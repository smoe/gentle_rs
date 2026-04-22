//! RNA-read UI state, caching, and read-inspection helpers for `MainAreaDna`.
//!
//! Keeping the RNA-read mapping/report support layer in its own submodule makes
//! later crate separation easier: the sequence window can depend on one
//! coherent mapping/report slice instead of a long interleaved block inside the
//! monolithic `main_area_dna.rs` file.
//!
//! Look here for:
//! - persisted RNA-read workspace UI state/defaults
//! - cached report/summary wrappers reused by the DNA window
//! - selection/filter/export helper types for `Mapped cDNA -> Read effects`

use super::*;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct RnaReadInterpretOpsUiState {
    pub(super) input_path: String,
    pub(super) report_id: String,
    #[serde(default = "super::default_true")]
    pub(super) report_id_auto_sync: bool,
    #[serde(default)]
    pub(super) report_id_auto_sync_explicit: bool,
    pub(super) scope: SplicingScopePreset,
    pub(super) profile: RnaReadInterpretationProfile,
    pub(super) input_format: RnaReadInputFormat,
    #[serde(default)]
    pub(super) origin_mode: RnaReadOriginMode,
    #[serde(default)]
    pub(super) target_gene_ids: String,
    #[serde(default)]
    pub(super) roi_seed_capture_enabled: bool,
    #[serde(default)]
    pub(super) report_mode: RnaReadReportMode,
    #[serde(default)]
    pub(super) checkpoint_path: String,
    #[serde(default = "super::default_rna_checkpoint_every_reads_text")]
    pub(super) checkpoint_every_reads: String,
    #[serde(default)]
    pub(super) resume_from_checkpoint: bool,
    #[serde(default = "super::default_true")]
    pub(super) cdna_poly_t_flip_enabled: bool,
    #[serde(default = "super::default_poly_t_prefix_min_bp_text")]
    pub(super) poly_t_prefix_min_bp: String,
    pub(super) kmer_len: String,
    #[serde(default = "super::default_rna_seed_stride_bp_text")]
    pub(super) seed_stride_bp: String,
    pub(super) min_seed_hit_fraction: String,
    pub(super) min_weighted_seed_hit_fraction: String,
    pub(super) min_unique_matched_kmers: String,
    pub(super) max_median_transcript_gap: String,
    #[serde(default = "super::default_min_chain_consistency_fraction_text")]
    pub(super) min_chain_consistency_fraction: String,
    pub(super) min_confirmed_exon_transitions: String,
    pub(super) min_transition_support_fraction: String,
    pub(super) align_band_width_bp: String,
    pub(super) align_min_identity_fraction: String,
    pub(super) align_max_secondary_mappings: String,
    #[serde(default = "super::default_rna_align_selection")]
    pub(super) align_phase_selection: RnaReadHitSelection,
    pub(super) show_advanced: bool,
}

impl Default for RnaReadInterpretOpsUiState {
    fn default() -> Self {
        Self {
            input_path: String::new(),
            report_id: String::new(),
            report_id_auto_sync: true,
            report_id_auto_sync_explicit: true,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_format: RnaReadInputFormat::Fasta,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: String::new(),
            roi_seed_capture_enabled: false,
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: String::new(),
            checkpoint_every_reads: "10000".to_string(),
            resume_from_checkpoint: false,
            cdna_poly_t_flip_enabled: true,
            poly_t_prefix_min_bp: "18".to_string(),
            kmer_len: "10".to_string(),
            seed_stride_bp: "1".to_string(),
            min_seed_hit_fraction: "0.30".to_string(),
            min_weighted_seed_hit_fraction: "0.05".to_string(),
            min_unique_matched_kmers: "12".to_string(),
            max_median_transcript_gap: "4.0".to_string(),
            min_chain_consistency_fraction: "0.40".to_string(),
            min_confirmed_exon_transitions: "1".to_string(),
            min_transition_support_fraction: "0.05".to_string(),
            align_band_width_bp: "24".to_string(),
            align_min_identity_fraction: "0.60".to_string(),
            align_max_secondary_mappings: "0".to_string(),
            align_phase_selection: RnaReadHitSelection::SeedPassed,
            show_advanced: false,
        }
    }
}

impl RnaReadInterpretOpsUiState {
    pub(super) fn normalize_loaded_state(&mut self) {
        // Legacy saved UI state predates the explicit auto/manual toggle. When
        // that marker is absent, treat non-empty report IDs as manual overrides
        // so reopening an older workspace does not silently rename reports.
        if !self.report_id_auto_sync_explicit && !self.report_id.trim().is_empty() {
            self.report_id_auto_sync = false;
        }
        self.report_id_auto_sync_explicit = true;
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct RnaReadEvidenceUiState {
    pub(super) selected_report_id: String,
    #[serde(default = "super::default_true")]
    pub(super) score_density_use_log_scale: bool,
    #[serde(default)]
    pub(super) score_density_variant: RnaReadScoreDensityVariant,
}

impl Default for RnaReadEvidenceUiState {
    fn default() -> Self {
        Self {
            selected_report_id: String::new(),
            score_density_use_log_scale: true,
            score_density_variant: RnaReadScoreDensityVariant::AllScored,
        }
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub(super) enum RnaReadConcatemerSubsetMode {
    #[default]
    AllMatchingRows,
    SelectedRowsOnly,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub(super) struct RnaReadConcatemerUiState {
    pub(super) selection: RnaReadHitSelection,
    pub(super) subset_mode: RnaReadConcatemerSubsetMode,
    pub(super) limit: String,
    pub(super) adapter_fasta_path: String,
    pub(super) transcript_fasta_paths_text: String,
    pub(super) transcript_index_paths_text: String,
    pub(super) show_advanced: bool,
    pub(super) internal_homopolymer_min_bp: String,
    pub(super) end_margin_bp: String,
    pub(super) max_primary_query_coverage_fraction: String,
    pub(super) min_secondary_identity_fraction: String,
    pub(super) max_secondary_query_overlap_fraction: String,
    pub(super) adapter_min_match_bp: String,
    pub(super) fragment_min_bp: String,
    pub(super) fragment_max_parts: String,
    pub(super) fragment_min_identity_fraction: String,
    pub(super) fragment_min_query_coverage_fraction: String,
}

impl Default for RnaReadConcatemerUiState {
    fn default() -> Self {
        let defaults = RnaReadConcatemerInspectionSettings::default();
        Self {
            selection: RnaReadHitSelection::Aligned,
            subset_mode: RnaReadConcatemerSubsetMode::AllMatchingRows,
            limit: "50".to_string(),
            adapter_fasta_path: String::new(),
            transcript_fasta_paths_text: String::new(),
            transcript_index_paths_text: String::new(),
            show_advanced: false,
            internal_homopolymer_min_bp: defaults.internal_homopolymer_min_bp.to_string(),
            end_margin_bp: defaults.end_margin_bp.to_string(),
            max_primary_query_coverage_fraction: format!(
                "{:.2}",
                defaults.max_primary_query_coverage_fraction
            ),
            min_secondary_identity_fraction: format!(
                "{:.2}",
                defaults.min_secondary_identity_fraction
            ),
            max_secondary_query_overlap_fraction: format!(
                "{:.2}",
                defaults.max_secondary_query_overlap_fraction
            ),
            adapter_min_match_bp: defaults.adapter_min_match_bp.to_string(),
            fragment_min_bp: defaults.fragment_min_bp.to_string(),
            fragment_max_parts: defaults.fragment_max_parts.to_string(),
            fragment_min_identity_fraction: format!(
                "{:.2}",
                defaults.fragment_min_identity_fraction
            ),
            fragment_min_query_coverage_fraction: format!(
                "{:.2}",
                defaults.fragment_min_query_coverage_fraction
            ),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum RnaReadEvidenceSourceTab {
    ReportedTranscript,
    ThresholdedCdna,
    MappedCdna,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum RnaReadMappedCdnaSubview {
    ReadEffects,
    AggregateSupport,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum RnaReadAlignmentEffectFilter {
    AllAligned,
    ConfirmedOnly,
    DisagreementOnly,
    ReassignedOnly,
    NoPhase1Only,
    SelectedOnly,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum RnaReadAlignmentEffectSortKey {
    Rank,
    Identity,
    Coverage,
    Score,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum RnaReadSelectedExportKind {
    Fasta,
    AlignmentsTsv,
    ExonPathsTsv,
    ExonAbundanceTsv,
}

impl RnaReadSelectedExportKind {
    pub(super) fn menu_label(self) -> &'static str {
        match self {
            Self::Fasta => "FASTA...",
            Self::AlignmentsTsv => "Alignments TSV...",
            Self::ExonPathsTsv => "Exon paths TSV...",
            Self::ExonAbundanceTsv => "Exon abundance TSV...",
        }
    }

    pub(super) fn missing_report_message(self) -> &'static str {
        match self {
            Self::Fasta => "Set a Report ID first to export selected reads",
            Self::AlignmentsTsv => "Set a Report ID first to export selected aligned rows",
            Self::ExonPathsTsv => "Set a Report ID first to export selected exon paths",
            Self::ExonAbundanceTsv => "Set a Report ID first to export selected exon abundance",
        }
    }

    pub(super) fn empty_selection_message(self) -> &'static str {
        match self {
            Self::Fasta => "Select one or more reads first to export selected FASTA",
            Self::AlignmentsTsv => "Select one or more aligned rows first to export them as TSV",
            Self::ExonPathsTsv => {
                "Select one or more aligned rows first to export their exon paths"
            }
            Self::ExonAbundanceTsv => {
                "Select one or more aligned rows first to export their exon abundance"
            }
        }
    }

    pub(super) fn default_file_name(self, report_id: &str) -> String {
        match self {
            Self::Fasta => format!("{report_id}_selected_reads.fa"),
            Self::AlignmentsTsv => format!("{report_id}_selected_alignments.tsv"),
            Self::ExonPathsTsv => format!("{report_id}_selected_exon_paths.tsv"),
            Self::ExonAbundanceTsv => format!("{report_id}_selected_exon_abundance.tsv"),
        }
    }

    pub(super) fn configure_dialog(self, default_name: &str) -> rfd::FileDialog {
        let dialog = rfd::FileDialog::new().set_file_name(default_name);
        match self {
            Self::Fasta => dialog.add_filter("FASTA", &["fa", "fasta"]),
            Self::AlignmentsTsv | Self::ExonPathsTsv | Self::ExonAbundanceTsv => {
                dialog.add_filter("TSV", &["tsv", "txt"])
            }
        }
    }

    pub(super) fn build_operation(
        self,
        report_id: String,
        path: String,
        selected_record_indices: Vec<usize>,
        subset_spec: Option<String>,
    ) -> Operation {
        match self {
            Self::Fasta => Operation::ExportRnaReadHitsFasta {
                report_id,
                path,
                selection: RnaReadHitSelection::Aligned,
                selected_record_indices,
                subset_spec,
            },
            Self::AlignmentsTsv => Operation::ExportRnaReadAlignmentsTsv {
                report_id,
                path,
                selection: RnaReadHitSelection::Aligned,
                limit: None,
                selected_record_indices,
                subset_spec,
            },
            Self::ExonPathsTsv => Operation::ExportRnaReadExonPathsTsv {
                report_id,
                path,
                selection: RnaReadHitSelection::Aligned,
                selected_record_indices,
                subset_spec,
            },
            Self::ExonAbundanceTsv => Operation::ExportRnaReadExonAbundanceTsv {
                report_id,
                path,
                selection: RnaReadHitSelection::Aligned,
                selected_record_indices,
                subset_spec,
            },
        }
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct RnaReadAlignmentEffectSummary {
    pub(super) aligned_rows: usize,
    pub(super) confirmed_assignments: usize,
    pub(super) reassigned_transcripts: usize,
    pub(super) aligned_without_phase1_assignment: usize,
    pub(super) seed_passed_but_unaligned: usize,
}

#[derive(Clone, Debug)]
pub(super) struct CachedRnaReadReport {
    pub(super) report_id: String,
    pub(super) report: Arc<RnaReadInterpretationReport>,
}

#[derive(Clone, Debug)]
pub(super) struct CachedRnaReadReportSummaries {
    pub(super) seq_id: String,
    pub(super) target_feature_id: usize,
    pub(super) summaries: Arc<Vec<RnaReadInterpretationReportSummary>>,
}

#[derive(Clone, Debug)]
pub(super) struct CachedRnaReadProgress {
    pub(super) report_id: String,
    pub(super) progress: Arc<RnaReadInterpretProgress>,
}

#[derive(Clone, Debug)]
pub(super) struct CachedRnaReadGeneSupportSummary {
    pub(super) cache_key: String,
    pub(super) result: Result<Arc<RnaReadGeneSupportSummary>, String>,
}

#[derive(Clone, Debug)]
pub(super) struct CachedRnaReadAlignmentInspection {
    pub(super) cache_key: String,
    pub(super) result: Result<Arc<RnaReadAlignmentInspection>, String>,
}

#[derive(Clone, Debug)]
pub(super) struct CachedRnaReadAlignmentDetail {
    pub(super) cache_key: String,
    pub(super) result: Result<Arc<RnaReadPairwiseAlignmentDetail>, String>,
}

#[derive(Clone, Debug)]
pub(super) struct CachedRnaReadAlignmentDisplay {
    pub(super) cache_key: String,
    pub(super) result: Result<Arc<RnaReadAlignmentDisplay>, String>,
}

#[derive(Clone, Debug)]
pub(super) struct CachedRnaReadConcatemerInspection {
    pub(super) cache_key: String,
    pub(super) result: Result<Arc<RnaReadConcatemerInspection>, String>,
}

impl MainAreaDna {
    pub(super) fn selected_rna_read_evidence_report_id(&self) -> Option<String> {
        let report_id = self.rna_read_evidence_ui.selected_report_id.trim();
        if report_id.is_empty() {
            None
        } else {
            Some(report_id.to_string())
        }
    }

    pub(super) fn current_rna_read_mapping_workspace_report_id(&self) -> Option<String> {
        let report_id = self.rna_reads_ui.report_id.trim();
        if report_id.is_empty() {
            None
        } else {
            Some(report_id.to_string())
        }
    }

    pub(super) fn sync_rna_read_evidence_selection_to_mapping_report(&mut self) {
        let Some(report_id) = self.current_rna_read_mapping_workspace_report_id() else {
            return;
        };
        if !self
            .rna_read_evidence_ui
            .selected_report_id
            .eq_ignore_ascii_case(&report_id)
        {
            self.rna_read_evidence_ui.selected_report_id = report_id;
        }
    }

    pub(super) fn invalidate_rna_read_report_display_cache(&mut self) {
        self.cached_saved_rna_read_report = None;
        self.cached_rna_read_report_summaries = None;
        self.cached_saved_rna_read_progress = None;
        self.cached_rna_read_gene_support_summary = None;
        self.cached_rna_read_alignment_inspections.clear();
        self.cached_rna_read_alignment_detail = None;
        self.cached_rna_read_alignment_display = None;
        self.cached_rna_read_concatemer_inspection = None;
        self.rna_read_alignment_detail_visible_key = None;
    }

    pub(super) fn rna_read_alignment_inspection_cache_key(
        report_id: &str,
        limit: usize,
        subset_spec: Option<&RnaReadAlignmentInspectionSubsetSpec>,
    ) -> String {
        let subset_spec_json = subset_spec
            .and_then(|spec| serde_json::to_string(spec).ok())
            .unwrap_or_else(|| "null".to_string());
        format!(
            "report={report_id}|limit={}|subset={subset_spec_json}",
            limit.max(1)
        )
    }

    pub(super) fn rna_read_gene_support_summary_cache_key(
        report_id: &str,
        gene_ids: &[String],
        selected_record_indices: &[usize],
        complete_rule: RnaReadGeneSupportCompleteRule,
    ) -> String {
        let gene_ids_json = serde_json::to_string(gene_ids).unwrap_or_else(|_| "[]".to_string());
        let selected_json =
            serde_json::to_string(selected_record_indices).unwrap_or_else(|_| "[]".to_string());
        format!(
            "report={report_id}|genes={gene_ids_json}|selected={selected_json}|complete_rule={}",
            complete_rule.as_str()
        )
    }

    pub(super) fn rna_read_concatemer_inspection_cache_key(
        report_id: &str,
        selection: RnaReadHitSelection,
        limit: usize,
        selected_record_indices: &[usize],
        settings: &RnaReadConcatemerInspectionSettings,
    ) -> String {
        let settings_json = serde_json::to_string(settings).unwrap_or_else(|_| "null".to_string());
        let selected_record_indices_json =
            serde_json::to_string(selected_record_indices).unwrap_or_else(|_| "[]".to_string());
        format!(
            "report={report_id}|selection={}|limit={}|selected={selected_record_indices_json}|settings={settings_json}",
            selection.as_str(),
            limit.max(1)
        )
    }

    pub(super) fn get_saved_rna_read_report_by_id(
        &mut self,
        report_id: &str,
    ) -> Option<Arc<RnaReadInterpretationReport>> {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            return None;
        }
        if let Some(cached) = self.cached_saved_rna_read_report.as_ref()
            && cached.report_id.eq_ignore_ascii_case(report_id)
        {
            return Some(cached.report.clone());
        }
        let report = {
            let engine = self.engine.as_ref()?;
            let guard = engine.read().ok()?;
            Arc::new(guard.get_rna_read_report(report_id).ok()?)
        };
        self.cached_saved_rna_read_report = Some(CachedRnaReadReport {
            report_id: report_id.to_string(),
            report: report.clone(),
        });
        if self
            .cached_saved_rna_read_progress
            .as_ref()
            .is_some_and(|cached| !cached.report_id.eq_ignore_ascii_case(report_id))
        {
            self.cached_saved_rna_read_progress = None;
        }
        Some(report)
    }

    pub(super) fn matching_rna_read_report_summaries_for_splicing_view(
        &mut self,
        view: &SplicingExpertView,
    ) -> Arc<Vec<RnaReadInterpretationReportSummary>> {
        if let Some(cached) = self.cached_rna_read_report_summaries.as_ref()
            && cached.seq_id == view.seq_id
            && cached.target_feature_id == view.target_feature_id
        {
            return cached.summaries.clone();
        }
        let mut rows = {
            let Some(engine) = self.engine.as_ref() else {
                return Arc::new(vec![]);
            };
            let Ok(guard) = engine.read() else {
                return Arc::new(vec![]);
            };
            guard.list_rna_read_reports(Some(&view.seq_id))
        };
        rows.retain(|row| row.seed_feature_id == view.target_feature_id);
        rows.sort_by(|left, right| {
            right
                .generated_at_unix_ms
                .cmp(&left.generated_at_unix_ms)
                .then_with(|| left.report_id.cmp(&right.report_id))
        });
        let rows = Arc::new(rows);
        self.cached_rna_read_report_summaries = Some(CachedRnaReadReportSummaries {
            seq_id: view.seq_id.clone(),
            target_feature_id: view.target_feature_id,
            summaries: rows.clone(),
        });
        rows
    }

    pub(super) fn latest_matching_rna_read_report_id(
        summaries: &[RnaReadInterpretationReportSummary],
    ) -> Option<String> {
        summaries.first().map(|row| row.report_id.clone())
    }

    pub(super) fn latest_rna_read_report_id_for_splicing_view(
        &mut self,
        view: &SplicingExpertView,
    ) -> Option<String> {
        let summaries = self.matching_rna_read_report_summaries_for_splicing_view(view);
        Self::latest_matching_rna_read_report_id(&summaries)
    }

    pub(super) fn ensure_selected_rna_read_evidence_report_for_view(
        &mut self,
        view: &SplicingExpertView,
    ) -> Arc<Vec<RnaReadInterpretationReportSummary>> {
        let summaries = self.matching_rna_read_report_summaries_for_splicing_view(view);
        let selected = self.rna_read_evidence_ui.selected_report_id.trim();
        let selected_matches = !selected.is_empty()
            && summaries
                .iter()
                .any(|row| row.report_id.eq_ignore_ascii_case(selected));
        if selected_matches {
            return summaries;
        }
        self.rna_read_evidence_ui.selected_report_id =
            Self::latest_matching_rna_read_report_id(&summaries).unwrap_or_default();
        summaries
    }

    pub(super) fn current_saved_rna_read_report(
        &mut self,
    ) -> Option<Arc<RnaReadInterpretationReport>> {
        self.selected_rna_read_evidence_report_id()
            .and_then(|report_id| self.get_saved_rna_read_report_by_id(&report_id))
    }

    pub(super) fn saved_rna_read_gene_support_summary_for_report_id(
        &mut self,
        report_id: &str,
        gene_ids: Vec<String>,
        selected_record_indices: Vec<usize>,
        complete_rule: RnaReadGeneSupportCompleteRule,
    ) -> Result<Arc<RnaReadGeneSupportSummary>, String> {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            return Err("Select a Report first before summarizing target-gene support".to_string());
        }
        let cache_key = Self::rna_read_gene_support_summary_cache_key(
            report_id,
            &gene_ids,
            &selected_record_indices,
            complete_rule,
        );
        if let Some(cached) = self.cached_rna_read_gene_support_summary.as_ref()
            && cached.cache_key == cache_key
        {
            return cached.result.as_ref().map(Arc::clone).map_err(Clone::clone);
        }
        let result = {
            let Some(engine) = &self.engine else {
                return Err("No engine attached".to_string());
            };
            match engine.read() {
                Ok(guard) => guard
                    .summarize_rna_read_gene_support(
                        report_id,
                        &gene_ids,
                        &selected_record_indices,
                        complete_rule,
                    )
                    .map(Arc::new)
                    .map_err(|e| e.message),
                Err(_) => {
                    Err("Engine lock poisoned while summarizing RNA-read gene support".to_string())
                }
            }
        };
        self.cached_rna_read_gene_support_summary = Some(CachedRnaReadGeneSupportSummary {
            cache_key,
            result: result.clone(),
        });
        result
    }

    pub(super) fn saved_rna_read_concatemer_inspection_for_report_id(
        &mut self,
        report_id: &str,
        selection: RnaReadHitSelection,
        limit: usize,
        selected_record_indices: Vec<usize>,
        settings: RnaReadConcatemerInspectionSettings,
    ) -> Result<Arc<RnaReadConcatemerInspection>, String> {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            return Err("Select a Report first before inspecting concatemer suspicion".to_string());
        }
        let cache_key = Self::rna_read_concatemer_inspection_cache_key(
            report_id,
            selection,
            limit,
            &selected_record_indices,
            &settings,
        );
        if let Some(cached) = self.cached_rna_read_concatemer_inspection.as_ref()
            && cached.cache_key == cache_key
        {
            return cached.result.as_ref().map(Arc::clone).map_err(Clone::clone);
        }
        let result = {
            let Some(engine) = &self.engine else {
                return Err("No engine attached".to_string());
            };
            match engine.read() {
                Ok(guard) => guard
                    .inspect_rna_read_concatemers(
                        report_id,
                        selection,
                        limit,
                        selected_record_indices,
                        settings,
                    )
                    .map(Arc::new)
                    .map_err(|e| e.message),
                Err(_) => Err(
                    "Engine lock poisoned while inspecting RNA-read concatemer suspicion"
                        .to_string(),
                ),
            }
        };
        self.cached_rna_read_concatemer_inspection = Some(CachedRnaReadConcatemerInspection {
            cache_key,
            result: result.clone(),
        });
        result
    }

    pub(super) fn saved_rna_read_alignment_inspection_for_report_id(
        &mut self,
        report_id: &str,
        limit: usize,
        subset_spec: Option<RnaReadAlignmentInspectionSubsetSpec>,
    ) -> Result<Arc<RnaReadAlignmentInspection>, String> {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            return Err("Select a Report first before inspecting aligned reads".to_string());
        }
        let cache_key =
            Self::rna_read_alignment_inspection_cache_key(report_id, limit, subset_spec.as_ref());
        if let Some(cached) = self
            .cached_rna_read_alignment_inspections
            .iter()
            .find(|row| row.cache_key == cache_key)
        {
            return cached.result.as_ref().map(Arc::clone).map_err(Clone::clone);
        }
        let result = {
            let Some(engine) = &self.engine else {
                return Err("No engine attached".to_string());
            };
            engine
                .read()
                .map_err(|_| {
                    "Engine lock poisoned while inspecting RNA-read alignments".to_string()
                })?
                .inspect_rna_read_alignments_with_subset(
                    &report_id,
                    RnaReadHitSelection::All,
                    limit.max(1),
                    subset_spec,
                )
                .map(Arc::new)
                .map_err(|error| error.message)
        };
        self.cached_rna_read_alignment_inspections
            .push(CachedRnaReadAlignmentInspection {
                cache_key,
                result: result.clone(),
            });
        if self.cached_rna_read_alignment_inspections.len() > 8 {
            self.cached_rna_read_alignment_inspections.remove(0);
        }
        result
    }

    pub(super) fn current_saved_rna_read_alignment_inspection(
        &mut self,
        limit: usize,
        subset_spec: Option<RnaReadAlignmentInspectionSubsetSpec>,
    ) -> Result<Arc<RnaReadAlignmentInspection>, String> {
        let report_id = self
            .selected_rna_read_evidence_report_id()
            .ok_or_else(|| "Select a Report first before inspecting aligned reads".to_string())?;
        self.saved_rna_read_alignment_inspection_for_report_id(&report_id, limit, subset_spec)
    }

    pub(super) fn rna_read_alignment_detail_cache_key(
        report_id: &str,
        record_index: usize,
    ) -> String {
        format!("{}|{}", report_id.trim(), record_index)
    }

    pub(super) fn rna_read_alignment_display_cache_key(
        report_id: &str,
        record_index: usize,
    ) -> String {
        format!("{}|{}", report_id.trim(), record_index)
    }

    pub(super) fn saved_rna_read_alignment_detail_for_report_id(
        &mut self,
        report_id: &str,
        record_index: usize,
    ) -> Result<Arc<RnaReadPairwiseAlignmentDetail>, String> {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            return Err("Select a Report first before inspecting alignment detail".to_string());
        }
        let cache_key = Self::rna_read_alignment_detail_cache_key(report_id, record_index);
        if let Some(cached) = self
            .cached_rna_read_alignment_detail
            .as_ref()
            .filter(|cached| cached.cache_key == cache_key)
        {
            return cached.result.as_ref().map(Arc::clone).map_err(Clone::clone);
        }
        let result = {
            let Some(engine) = &self.engine else {
                return Err("No engine attached".to_string());
            };
            engine
                .read()
                .map_err(|_| {
                    "Engine lock poisoned while inspecting RNA-read alignment detail".to_string()
                })?
                .inspect_rna_read_alignment_detail(report_id, record_index)
                .map(Arc::new)
                .map_err(|error| error.message)
        };
        self.cached_rna_read_alignment_detail = Some(CachedRnaReadAlignmentDetail {
            cache_key,
            result: result.clone(),
        });
        result
    }

    pub(super) fn saved_rna_read_alignment_display_for_report_id(
        &mut self,
        report_id: &str,
        record_index: usize,
    ) -> Result<Arc<RnaReadAlignmentDisplay>, String> {
        let report_id = report_id.trim();
        if report_id.is_empty() {
            return Err(
                "Select a Report first before inspecting phase-2 pairwise alignment".to_string(),
            );
        }
        let cache_key = Self::rna_read_alignment_display_cache_key(report_id, record_index);
        if let Some(cached) = self
            .cached_rna_read_alignment_display
            .as_ref()
            .filter(|cached| cached.cache_key == cache_key)
        {
            return cached.result.as_ref().map(Arc::clone).map_err(Clone::clone);
        }
        let result = {
            let Some(engine) = &self.engine else {
                return Err(
                    "No engine attached while building RNA-read alignment detail".to_string(),
                );
            };
            engine
                .read()
                .map_err(|_| {
                    "Engine lock poisoned while building RNA-read alignment detail".to_string()
                })?
                .build_rna_read_alignment_display(report_id, record_index)
                .map(Arc::new)
                .map_err(|error| error.message)
        };
        self.cached_rna_read_alignment_display = Some(CachedRnaReadAlignmentDisplay {
            cache_key,
            result: result.clone(),
        });
        result
    }

    pub(super) fn selected_highlighted_rna_report_hit<'a>(
        &self,
        report: &'a RnaReadInterpretationReport,
    ) -> Option<&'a RnaReadInterpretationHit> {
        let selected_index = self.rna_seed_highlight_record_index?;
        report
            .hits
            .iter()
            .find(|hit| hit.record_index == selected_index)
    }

    pub(super) fn rna_read_top_hit_preview_from_hit(
        hit: &RnaReadInterpretationHit,
    ) -> RnaReadTopHitPreview {
        let mapping = hit.best_mapping.as_ref();
        let sequence_preview = if hit.sequence.len() > 120 {
            format!("{}...", &hit.sequence[..120])
        } else {
            hit.sequence.clone()
        };
        RnaReadTopHitPreview {
            record_index: hit.record_index,
            header_id: hit.header_id.clone(),
            seed_hit_fraction: hit.seed_hit_fraction,
            weighted_seed_hit_fraction: hit.weighted_seed_hit_fraction,
            weighted_matched_kmers: hit.weighted_matched_kmers,
            seed_chain_transcript_id: hit.seed_chain_transcript_id.clone(),
            seed_chain_support_kmers: hit.seed_chain_support_kmers,
            seed_chain_support_fraction: hit.seed_chain_support_fraction,
            seed_median_transcript_gap: hit.seed_median_transcript_gap,
            seed_transcript_gap_count: hit.seed_transcript_gap_count,
            matched_kmers: hit.matched_kmers,
            tested_kmers: hit.tested_kmers,
            passed_seed_filter: hit.passed_seed_filter,
            reverse_complement_applied: hit.reverse_complement_applied,
            selected_strand: hit.strand_diagnostics.selected_strand.clone(),
            competing_opposite_strand: hit.strand_diagnostics.competing_opposite_strand,
            ambiguous_strand_tie: hit.strand_diagnostics.ambiguous_near_tie,
            origin_class: hit.origin_class,
            origin_reason: hit.origin_reason.clone(),
            origin_confidence: hit.origin_confidence,
            strand_confidence: hit.strand_confidence,
            origin_candidates: hit.origin_candidates.clone(),
            msa_eligible: hit.msa_eligible,
            msa_eligibility_reason: hit.msa_eligibility_reason.clone(),
            aligned: mapping.is_some(),
            best_alignment_mode: mapping
                .map(|row| row.alignment_mode.as_str().to_string())
                .unwrap_or_default(),
            best_alignment_transcript_id: mapping
                .map(|row| row.transcript_id.clone())
                .unwrap_or_default(),
            best_alignment_transcript_label: mapping
                .map(|row| row.transcript_label.clone())
                .unwrap_or_default(),
            best_alignment_strand: mapping.map(|row| row.strand.clone()).unwrap_or_default(),
            best_alignment_target_start_1based: mapping
                .map(|row| row.target_start_1based)
                .unwrap_or_default(),
            best_alignment_target_end_1based: mapping
                .map(|row| row.target_end_1based)
                .unwrap_or_default(),
            best_alignment_identity_fraction: mapping
                .map(|row| row.identity_fraction)
                .unwrap_or_default(),
            best_alignment_query_coverage_fraction: mapping
                .map(|row| row.query_coverage_fraction)
                .unwrap_or_default(),
            best_alignment_score: mapping.map(|row| row.score).unwrap_or_default(),
            secondary_mapping_count: hit.secondary_mappings.len(),
            read_length_bp: hit.read_length_bp,
            sequence: hit.sequence.clone(),
            sequence_preview,
        }
    }

    pub(super) fn mapped_support_frequencies_from_inspection(
        inspection: &RnaReadAlignmentInspection,
    ) -> (
        Vec<RnaReadExonSupportFrequency>,
        Vec<RnaReadJunctionSupportFrequency>,
    ) {
        let mut exon_counts = BTreeMap::<(usize, usize), usize>::new();
        let mut junction_counts = BTreeMap::<(usize, usize), usize>::new();
        let denominator = inspection.aligned_count.max(1) as f64;
        for row in &inspection.rows {
            let mut seen_exons = BTreeSet::<(usize, usize)>::new();
            for exon in &row.mapped_exon_support {
                if seen_exons.insert((exon.start_1based, exon.end_1based)) {
                    *exon_counts
                        .entry((exon.start_1based, exon.end_1based))
                        .or_insert(0) += 1;
                }
            }
            let mut seen_junctions = BTreeSet::<(usize, usize)>::new();
            for junction in &row.mapped_junction_support {
                if seen_junctions.insert((junction.donor_1based, junction.acceptor_1based)) {
                    *junction_counts
                        .entry((junction.donor_1based, junction.acceptor_1based))
                        .or_insert(0) += 1;
                }
            }
        }
        let exon_rows = exon_counts
            .into_iter()
            .map(
                |((start_1based, end_1based), support_read_count)| RnaReadExonSupportFrequency {
                    start_1based,
                    end_1based,
                    support_read_count,
                    support_fraction: support_read_count as f64 / denominator,
                },
            )
            .collect::<Vec<_>>();
        let junction_rows = junction_counts
            .into_iter()
            .map(|((donor_1based, acceptor_1based), support_read_count)| {
                RnaReadJunctionSupportFrequency {
                    donor_1based,
                    acceptor_1based,
                    support_read_count,
                    support_fraction: support_read_count as f64 / denominator,
                }
            })
            .collect::<Vec<_>>();
        (exon_rows, junction_rows)
    }

    pub(super) fn synthesize_rna_read_progress_from_report(
        &mut self,
        report: &RnaReadInterpretationReport,
    ) -> Arc<RnaReadInterpretProgress> {
        if let Some(cached) = self.cached_saved_rna_read_progress.as_ref()
            && cached.report_id.eq_ignore_ascii_case(&report.report_id)
        {
            return cached.progress.clone();
        }
        let mut lengths = report
            .hits
            .iter()
            .map(|hit| hit.read_length_bp)
            .collect::<Vec<_>>();
        lengths.sort_unstable();
        let total_reads = report.read_count_total.max(report.hits.len());
        let total_bases = report
            .hits
            .iter()
            .map(|hit| hit.read_length_bp as u64)
            .sum::<u64>();
        let mean_read_length_bp = if report.hits.is_empty() {
            0.0
        } else {
            total_bases as f64 / report.hits.len() as f64
        };
        let median_read_length_bp = if lengths.is_empty() {
            0
        } else {
            lengths[lengths.len() / 2]
        };
        let p95_read_length_bp = if lengths.is_empty() {
            0
        } else {
            let idx = ((lengths.len() - 1) as f64 * 0.95).round() as usize;
            lengths[idx.min(lengths.len() - 1)]
        };
        let top_hits_preview = report
            .hits
            .iter()
            .take(256)
            .map(Self::rna_read_top_hit_preview_from_hit)
            .collect::<Vec<_>>();
        let has_saved_alignments = report.read_count_aligned > 0
            || report.hits.iter().any(|hit| hit.best_mapping.is_some());
        let inspection = if has_saved_alignments {
            self.saved_rna_read_alignment_inspection_for_report_id(
                &report.report_id,
                report.hits.len().max(1),
                None,
            )
            .ok()
        } else {
            None
        };
        let (mapped_exon_support_frequencies, mapped_junction_support_frequencies) = inspection
            .as_ref()
            .map(|inspection| Self::mapped_support_frequencies_from_inspection(inspection.as_ref()))
            .unwrap_or_else(|| (vec![], vec![]));
        let reads_with_transition_support = report
            .hits
            .iter()
            .filter(|hit| hit.exon_transitions_confirmed > 0)
            .count();
        let transition_confirmations = report
            .hits
            .iter()
            .map(|hit| hit.exon_transitions_confirmed)
            .sum::<usize>();
        let tested_kmers = report
            .hits
            .iter()
            .map(|hit| hit.tested_kmers)
            .sum::<usize>();
        let matched_kmers = report
            .hits
            .iter()
            .map(|hit| hit.matched_kmers)
            .sum::<usize>();
        let progress = Arc::new(RnaReadInterpretProgress {
            seq_id: report.seq_id.clone(),
            reads_processed: total_reads,
            reads_total: total_reads,
            read_bases_processed: total_bases,
            mean_read_length_bp,
            median_read_length_bp,
            p95_read_length_bp,
            input_bytes_processed: 0,
            input_bytes_total: 0,
            seed_passed: report.read_count_seed_passed,
            aligned: report.read_count_aligned,
            tested_kmers,
            matched_kmers,
            seed_compute_ms: 0.0,
            align_compute_ms: 0.0,
            io_read_ms: 0.0,
            fasta_parse_ms: 0.0,
            normalize_compute_ms: 0.0,
            inference_compute_ms: 0.0,
            progress_emit_ms: 0.0,
            update_every_reads: 0,
            done: true,
            bins: vec![],
            score_density_bins: report.score_density_bins.clone(),
            seed_pass_score_density_bins: report.seed_pass_score_density_bins.clone(),
            top_hits_preview,
            transition_support_rows: report.transition_support_rows.clone(),
            isoform_support_rows: report.isoform_support_rows.clone(),
            mapped_exon_support_frequencies,
            mapped_junction_support_frequencies,
            mapped_isoform_support_rows: report.mapped_isoform_support_rows.clone(),
            reads_with_transition_support,
            transition_confirmations,
            junction_crossing_seed_bits_indexed: 0,
            origin_class_counts: report.origin_class_counts.clone(),
        });
        self.cached_saved_rna_read_progress = Some(CachedRnaReadProgress {
            report_id: report.report_id.clone(),
            progress: progress.clone(),
        });
        progress
    }

    pub(super) fn active_rna_read_task_matches_splicing_view(
        &self,
        view: &SplicingExpertView,
    ) -> bool {
        self.rna_read_task
            .as_ref()
            .map(|task| {
                task.seq_id == view.seq_id && task.seed_feature_id == view.target_feature_id
            })
            .unwrap_or(false)
    }

    pub(super) fn current_rna_read_evidence_progress_for_view(
        &mut self,
        view: &SplicingExpertView,
        report: Option<&RnaReadInterpretationReport>,
    ) -> Option<Arc<RnaReadInterpretProgress>> {
        let selected_report_matches_active = self
            .rna_read_task
            .as_ref()
            .and_then(|task| task.report_id_hint.as_deref())
            .zip(self.selected_rna_read_evidence_report_id())
            .map(|(task_report_id, selected_report_id)| {
                task_report_id.eq_ignore_ascii_case(selected_report_id.as_str())
            })
            .unwrap_or(false);
        if self.active_rna_read_task_matches_splicing_view(view)
            && selected_report_matches_active
            && let Some(progress) = self.rna_read_progress.clone()
        {
            return Some(Arc::new(progress));
        }
        report.map(|row| self.synthesize_rna_read_progress_from_report(row))
    }

    pub(super) fn current_rna_read_mapping_progress_for_view(
        &mut self,
        view: &SplicingExpertView,
    ) -> Option<Arc<RnaReadInterpretProgress>> {
        if self.active_rna_read_task_matches_splicing_view(view)
            && let Some(progress) = self.rna_read_progress.clone()
        {
            return Some(Arc::new(progress));
        }
        self.current_rna_read_mapping_workspace_report_id()
            .and_then(|report_id| self.get_saved_rna_read_report_by_id(&report_id))
            .map(|report| self.synthesize_rna_read_progress_from_report(report.as_ref()))
    }

    pub(super) fn summarize_rna_read_alignment_effects(
        report: &RnaReadInterpretationReport,
        inspection: &RnaReadAlignmentInspection,
    ) -> RnaReadAlignmentEffectSummary {
        let mut summary = RnaReadAlignmentEffectSummary {
            aligned_rows: inspection.aligned_count,
            ..RnaReadAlignmentEffectSummary::default()
        };
        for row in &inspection.rows {
            match row.alignment_effect {
                RnaReadAlignmentEffect::ConfirmedAssignment => {
                    summary.confirmed_assignments = summary.confirmed_assignments.saturating_add(1);
                }
                RnaReadAlignmentEffect::ReassignedTranscript => {
                    summary.reassigned_transcripts =
                        summary.reassigned_transcripts.saturating_add(1);
                }
                RnaReadAlignmentEffect::AlignedWithoutPhase1Assignment => {
                    summary.aligned_without_phase1_assignment =
                        summary.aligned_without_phase1_assignment.saturating_add(1);
                }
            }
        }
        summary.seed_passed_but_unaligned = report
            .hits
            .iter()
            .filter(|hit| hit.passed_seed_filter && hit.best_mapping.is_none())
            .count();
        summary
    }

    pub(super) fn rna_read_alignment_effect_label(effect: RnaReadAlignmentEffect) -> &'static str {
        match effect {
            RnaReadAlignmentEffect::ConfirmedAssignment => "confirmed",
            RnaReadAlignmentEffect::ReassignedTranscript => "reassigned",
            RnaReadAlignmentEffect::AlignedWithoutPhase1Assignment => "aligned (no phase-1 tx)",
        }
    }

    pub(super) fn rna_read_alignment_effect_filter_label(
        filter: RnaReadAlignmentEffectFilter,
    ) -> &'static str {
        match filter {
            RnaReadAlignmentEffectFilter::AllAligned => "all aligned",
            RnaReadAlignmentEffectFilter::ConfirmedOnly => "confirmed only",
            RnaReadAlignmentEffectFilter::DisagreementOnly => "disagreement only",
            RnaReadAlignmentEffectFilter::ReassignedOnly => "reassigned only",
            RnaReadAlignmentEffectFilter::NoPhase1Only => "no phase-1 tx",
            RnaReadAlignmentEffectFilter::SelectedOnly => "selected only",
        }
    }

    pub(super) fn rna_read_alignment_effect_sort_key_label(
        sort_key: RnaReadAlignmentEffectSortKey,
    ) -> &'static str {
        match sort_key {
            RnaReadAlignmentEffectSortKey::Rank => "rank",
            RnaReadAlignmentEffectSortKey::Identity => "identity",
            RnaReadAlignmentEffectSortKey::Coverage => "coverage",
            RnaReadAlignmentEffectSortKey::Score => "score",
        }
    }

    pub(super) fn engine_rna_read_alignment_effect_filter(
        filter: RnaReadAlignmentEffectFilter,
    ) -> RnaReadAlignmentInspectionEffectFilter {
        match filter {
            RnaReadAlignmentEffectFilter::AllAligned => {
                RnaReadAlignmentInspectionEffectFilter::AllAligned
            }
            RnaReadAlignmentEffectFilter::ConfirmedOnly => {
                RnaReadAlignmentInspectionEffectFilter::ConfirmedOnly
            }
            RnaReadAlignmentEffectFilter::DisagreementOnly => {
                RnaReadAlignmentInspectionEffectFilter::DisagreementOnly
            }
            RnaReadAlignmentEffectFilter::ReassignedOnly => {
                RnaReadAlignmentInspectionEffectFilter::ReassignedOnly
            }
            RnaReadAlignmentEffectFilter::NoPhase1Only => {
                RnaReadAlignmentInspectionEffectFilter::NoPhase1Only
            }
            RnaReadAlignmentEffectFilter::SelectedOnly => {
                RnaReadAlignmentInspectionEffectFilter::SelectedOnly
            }
        }
    }

    pub(super) fn engine_rna_read_alignment_sort_key(
        sort_key: RnaReadAlignmentEffectSortKey,
    ) -> RnaReadAlignmentInspectionSortKey {
        match sort_key {
            RnaReadAlignmentEffectSortKey::Rank => RnaReadAlignmentInspectionSortKey::Rank,
            RnaReadAlignmentEffectSortKey::Identity => RnaReadAlignmentInspectionSortKey::Identity,
            RnaReadAlignmentEffectSortKey::Coverage => RnaReadAlignmentInspectionSortKey::Coverage,
            RnaReadAlignmentEffectSortKey::Score => RnaReadAlignmentInspectionSortKey::Score,
        }
    }

    pub(super) fn current_rna_read_alignment_subset_spec(
        &mut self,
    ) -> RnaReadAlignmentInspectionSubsetSpec {
        let score_bin_count = self
            .current_saved_rna_read_report()
            .map(|report| report.score_density_bins.len().max(40))
            .unwrap_or(40);
        let score_density_seed_filter_override = self
            .current_rna_read_score_density_seed_filter_override(
                self.rna_read_evidence_ui.score_density_variant,
            )
            .ok()
            .flatten();
        Self::rna_read_alignment_effect_subset_struct(
            self.rna_read_alignment_effect_filter,
            &self.rna_read_alignment_effect_search,
            self.rna_read_alignment_effect_sort_key,
            &self.rna_seed_selected_record_indices,
            self.rna_read_evidence_ui.score_density_variant,
            score_density_seed_filter_override,
            self.rna_read_alignment_effect_score_bin_index,
            score_bin_count,
        )
    }

    pub(super) fn current_rna_read_score_density_seed_filter_override(
        &self,
        variant: RnaReadScoreDensityVariant,
    ) -> Result<Option<RnaReadSeedFilterConfig>, String> {
        match variant {
            RnaReadScoreDensityVariant::RetainedReplayCurrentControls => {
                self.parse_rna_seed_filter_from_ui().map(Some)
            }
            _ => Ok(None),
        }
    }

    pub(super) fn rna_read_alignment_effect_subset_struct(
        filter: RnaReadAlignmentEffectFilter,
        search: &str,
        sort_key: RnaReadAlignmentEffectSortKey,
        selected_record_indices: &BTreeSet<usize>,
        score_density_variant: RnaReadScoreDensityVariant,
        score_density_seed_filter_override: Option<RnaReadSeedFilterConfig>,
        score_bin_index: Option<usize>,
        score_bin_count: usize,
    ) -> RnaReadAlignmentInspectionSubsetSpec {
        RnaReadAlignmentInspectionSubsetSpec {
            effect_filter: Self::engine_rna_read_alignment_effect_filter(filter),
            sort_key: Self::engine_rna_read_alignment_sort_key(sort_key),
            search: search.trim().to_string(),
            selected_record_indices: selected_record_indices.iter().copied().collect(),
            score_density_variant,
            score_density_seed_filter_override,
            score_bin_index,
            score_bin_count: score_bin_count.max(1),
        }
    }

    #[cfg(test)]
    pub(super) fn rna_read_alignment_effect_matches_filter(
        row: &RnaReadAlignmentInspectionRow,
        filter: RnaReadAlignmentEffectFilter,
        selected_record_indices: &BTreeSet<usize>,
    ) -> bool {
        match filter {
            RnaReadAlignmentEffectFilter::AllAligned => true,
            RnaReadAlignmentEffectFilter::ConfirmedOnly => {
                row.alignment_effect == RnaReadAlignmentEffect::ConfirmedAssignment
            }
            RnaReadAlignmentEffectFilter::DisagreementOnly => {
                row.alignment_effect != RnaReadAlignmentEffect::ConfirmedAssignment
            }
            RnaReadAlignmentEffectFilter::ReassignedOnly => {
                row.alignment_effect == RnaReadAlignmentEffect::ReassignedTranscript
            }
            RnaReadAlignmentEffectFilter::NoPhase1Only => {
                row.alignment_effect == RnaReadAlignmentEffect::AlignedWithoutPhase1Assignment
            }
            RnaReadAlignmentEffectFilter::SelectedOnly => {
                selected_record_indices.contains(&row.record_index)
            }
        }
    }

    #[cfg(test)]
    pub(super) fn rna_read_alignment_effect_matches_search(
        row: &RnaReadAlignmentInspectionRow,
        search: &str,
    ) -> bool {
        let needle = search.trim().to_ascii_lowercase();
        if needle.is_empty() {
            return true;
        }
        let effect_label = Self::rna_read_alignment_effect_label(row.alignment_effect);
        let record_index_label = format!("#{}", row.record_index + 1);
        let rank_label = row.rank.to_string();
        [
            row.header_id.as_str(),
            row.phase1_primary_transcript_id.as_str(),
            row.seed_chain_transcript_id.as_str(),
            row.exon_path_transcript_id.as_str(),
            row.exon_path.as_str(),
            row.transcript_id.as_str(),
            row.transcript_label.as_str(),
            row.strand.as_str(),
            row.selected_strand.as_str(),
            row.origin_class.as_str(),
            effect_label,
            record_index_label.as_str(),
            rank_label.as_str(),
        ]
        .iter()
        .any(|field| field.to_ascii_lowercase().contains(&needle))
    }

    #[cfg(test)]
    pub(super) fn compare_rna_read_alignment_effect_rows(
        left: &RnaReadAlignmentInspectionRow,
        right: &RnaReadAlignmentInspectionRow,
        sort_key: RnaReadAlignmentEffectSortKey,
    ) -> Ordering {
        match sort_key {
            RnaReadAlignmentEffectSortKey::Rank => left
                .rank
                .cmp(&right.rank)
                .then_with(|| left.record_index.cmp(&right.record_index)),
            RnaReadAlignmentEffectSortKey::Identity => right
                .identity_fraction
                .partial_cmp(&left.identity_fraction)
                .unwrap_or(Ordering::Equal)
                .then_with(|| {
                    right
                        .query_coverage_fraction
                        .partial_cmp(&left.query_coverage_fraction)
                        .unwrap_or(Ordering::Equal)
                })
                .then_with(|| right.score.cmp(&left.score))
                .then_with(|| left.rank.cmp(&right.rank)),
            RnaReadAlignmentEffectSortKey::Coverage => right
                .query_coverage_fraction
                .partial_cmp(&left.query_coverage_fraction)
                .unwrap_or(Ordering::Equal)
                .then_with(|| {
                    right
                        .identity_fraction
                        .partial_cmp(&left.identity_fraction)
                        .unwrap_or(Ordering::Equal)
                })
                .then_with(|| right.score.cmp(&left.score))
                .then_with(|| left.rank.cmp(&right.rank)),
            RnaReadAlignmentEffectSortKey::Score => right
                .score
                .cmp(&left.score)
                .then_with(|| {
                    right
                        .identity_fraction
                        .partial_cmp(&left.identity_fraction)
                        .unwrap_or(Ordering::Equal)
                })
                .then_with(|| {
                    right
                        .query_coverage_fraction
                        .partial_cmp(&left.query_coverage_fraction)
                        .unwrap_or(Ordering::Equal)
                })
                .then_with(|| left.rank.cmp(&right.rank)),
        }
    }

    #[cfg(test)]
    pub(super) fn collect_visible_rna_read_alignment_effect_rows<'a>(
        inspection: &'a RnaReadAlignmentInspection,
        filter: RnaReadAlignmentEffectFilter,
        search: &str,
        selected_record_indices: &BTreeSet<usize>,
        sort_key: RnaReadAlignmentEffectSortKey,
    ) -> Vec<&'a RnaReadAlignmentInspectionRow> {
        let mut rows = inspection
            .rows
            .iter()
            .filter(|row| {
                Self::rna_read_alignment_effect_matches_filter(row, filter, selected_record_indices)
                    && Self::rna_read_alignment_effect_matches_search(row, search)
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            Self::compare_rna_read_alignment_effect_rows(left, right, sort_key)
        });
        rows
    }

    pub(super) fn collect_rna_read_alignment_effect_record_indices(
        rows: &[&RnaReadAlignmentInspectionRow],
    ) -> Vec<usize> {
        rows.iter().map(|row| row.record_index).collect::<Vec<_>>()
    }

    pub(super) fn rna_read_alignment_effect_subset_spec(
        filter: RnaReadAlignmentEffectFilter,
        search: &str,
        sort_key: RnaReadAlignmentEffectSortKey,
        score_density_variant: RnaReadScoreDensityVariant,
        score_bin_index: Option<usize>,
        score_bin_count: usize,
    ) -> String {
        let search = search.trim();
        let score_bin = score_bin_index
            .map(|idx| Self::format_rna_read_score_bin_spec(idx, score_bin_count))
            .unwrap_or_else(|| "<none>".to_string());
        format!(
            "filter={} | histogram={} | score_bin={} | sort={} | search={}",
            Self::rna_read_alignment_effect_filter_label(filter),
            Self::rna_read_score_density_variant_label(score_density_variant),
            score_bin,
            Self::rna_read_alignment_effect_sort_key_label(sort_key),
            if search.is_empty() { "<none>" } else { search }
        )
    }

    pub(super) fn rna_read_score_density_variant_label(
        variant: RnaReadScoreDensityVariant,
    ) -> &'static str {
        match variant {
            RnaReadScoreDensityVariant::AllScored => "all scored",
            RnaReadScoreDensityVariant::CompositeSeedGate => "composite gate",
            RnaReadScoreDensityVariant::RetainedReplayCurrentControls => "retained replay",
        }
    }

    pub(super) fn rna_read_score_density_variant_note(
        variant: RnaReadScoreDensityVariant,
    ) -> &'static str {
        match variant {
            RnaReadScoreDensityVariant::AllScored => {
                "Histogram over all scored reads from phase 1."
            }
            RnaReadScoreDensityVariant::CompositeSeedGate => {
                "Histogram over reads that passed the full composite seed gate during phase 1."
            }
            RnaReadScoreDensityVariant::RetainedReplayCurrentControls => {
                "Instant replay over retained saved-report rows using the current seed controls. Reads never retained by the original run are not revisited."
            }
        }
    }

    pub(super) fn focus_rna_read_alignment_effect_record_indices(
        &mut self,
        record_indices: Vec<usize>,
        source_label: &str,
    ) {
        self.rna_read_statistics_tab = RnaReadEvidenceSourceTab::MappedCdna;
        self.rna_read_mapped_cdna_subview = RnaReadMappedCdnaSubview::ReadEffects;
        self.rna_read_alignment_effect_filter = RnaReadAlignmentEffectFilter::SelectedOnly;
        self.rna_read_alignment_effect_sort_key = RnaReadAlignmentEffectSortKey::Score;
        self.rna_read_alignment_effect_search.clear();
        self.rna_read_alignment_effect_score_bin_index = None;
        self.rna_seed_selected_record_indices = record_indices.iter().copied().collect();
        self.rna_seed_highlight_record_index = record_indices.first().copied();
        self.op_status = format!(
            "Focused read effects on {} contributing read(s) from {source_label}",
            self.rna_seed_selected_record_indices.len()
        );
    }

    pub(super) fn clear_rna_read_alignment_score_bin_focus(&mut self) {
        self.rna_read_alignment_effect_score_bin_index = None;
    }

    pub(super) fn rna_read_score_density_bin_bounds(
        bin_index: usize,
        bin_count: usize,
    ) -> (f64, f64) {
        let bin_count = bin_count.max(1);
        let left = bin_index.min(bin_count.saturating_sub(1)) as f64 / bin_count as f64;
        let right = (bin_index.min(bin_count.saturating_sub(1)) + 1) as f64 / bin_count as f64;
        (left, right)
    }

    pub(super) fn format_rna_read_score_bin_spec(bin_index: usize, bin_count: usize) -> String {
        let (left, right) = Self::rna_read_score_density_bin_bounds(bin_index, bin_count);
        format!(
            "{}/{} [{left:.3},{right:.3})",
            bin_index.min(bin_count.saturating_sub(1)),
            bin_count.max(1)
        )
    }

    pub(super) fn rna_read_score_density_bins_for_progress(
        progress: &RnaReadInterpretProgress,
        variant: RnaReadScoreDensityVariant,
    ) -> Vec<u64> {
        match variant {
            RnaReadScoreDensityVariant::AllScored => progress.score_density_bins.clone(),
            RnaReadScoreDensityVariant::CompositeSeedGate => {
                progress.seed_pass_score_density_bins.clone()
            }
            RnaReadScoreDensityVariant::RetainedReplayCurrentControls => vec![],
        }
    }

    pub(super) fn collect_rna_read_mapped_exon_contributors(
        inspection: &RnaReadAlignmentInspection,
    ) -> HashMap<(usize, usize), Vec<usize>> {
        let mut contributors = HashMap::<(usize, usize), Vec<usize>>::new();
        for row in &inspection.rows {
            for exon in &row.mapped_exon_support {
                let entry = contributors
                    .entry((exon.start_1based, exon.end_1based))
                    .or_default();
                if !entry.contains(&row.record_index) {
                    entry.push(row.record_index);
                }
            }
        }
        contributors
    }

    pub(super) fn collect_rna_read_mapped_junction_contributors(
        inspection: &RnaReadAlignmentInspection,
    ) -> HashMap<(usize, usize), Vec<usize>> {
        let mut contributors = HashMap::<(usize, usize), Vec<usize>>::new();
        for row in &inspection.rows {
            for junction in &row.mapped_junction_support {
                let entry = contributors
                    .entry((junction.donor_1based, junction.acceptor_1based))
                    .or_default();
                if !entry.contains(&row.record_index) {
                    entry.push(row.record_index);
                }
            }
        }
        contributors
    }

    pub(super) fn collect_rna_read_mapped_isoform_contributors(
        inspection: &RnaReadAlignmentInspection,
    ) -> HashMap<String, Vec<usize>> {
        let mut contributors = HashMap::<String, Vec<usize>>::new();
        for row in &inspection.rows {
            let transcript_id = row.transcript_id.trim();
            if transcript_id.is_empty() {
                continue;
            }
            let entry = contributors.entry(transcript_id.to_string()).or_default();
            if !entry.contains(&row.record_index) {
                entry.push(row.record_index);
            }
        }
        contributors
    }

    pub(super) fn format_rna_read_contributor_hover_text(
        inspection: &RnaReadAlignmentInspection,
        record_indices: &[usize],
    ) -> String {
        let row_by_record_index = inspection
            .rows
            .iter()
            .map(|row| (row.record_index, row))
            .collect::<HashMap<_, _>>();
        record_indices
            .iter()
            .filter_map(|record_index| {
                row_by_record_index.get(record_index).map(|row| {
                    format!(
                        "#{} {} | {} | {}",
                        row.record_index + 1,
                        row.header_id,
                        Self::rna_read_alignment_effect_label(row.alignment_effect),
                        Self::compact_rna_read_transcript_label(
                            &row.transcript_id,
                            &row.transcript_label,
                        )
                    )
                })
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    pub(super) fn compact_rna_read_transcript_label(
        transcript_id: &str,
        transcript_label: &str,
    ) -> String {
        let id = transcript_id.trim();
        let label = transcript_label.trim();
        if id.is_empty() && label.is_empty() {
            "none".to_string()
        } else if id.is_empty() {
            label.to_string()
        } else if label.is_empty() || label.eq_ignore_ascii_case(id) {
            id.to_string()
        } else {
            format!("{label} ({id})")
        }
    }

    pub(super) fn format_rna_read_mapped_exon_support_compact(
        row: &RnaReadAlignmentInspectionRow,
    ) -> String {
        if row.mapped_exon_support.is_empty() {
            "none".to_string()
        } else {
            row.mapped_exon_support
                .iter()
                .map(|entry| format!("{}..{}", entry.start_1based, entry.end_1based))
                .collect::<Vec<_>>()
                .join(", ")
        }
    }

    pub(super) fn format_rna_read_mapped_junction_support_compact(
        row: &RnaReadAlignmentInspectionRow,
    ) -> String {
        if row.mapped_junction_support.is_empty() {
            "none".to_string()
        } else {
            row.mapped_junction_support
                .iter()
                .map(|entry| format!("{}->{}", entry.donor_1based, entry.acceptor_1based))
                .collect::<Vec<_>>()
                .join(", ")
        }
    }

    pub(super) fn rna_read_score_density_bin_index(
        seed_hit_fraction: f64,
        bin_count: usize,
    ) -> usize {
        if bin_count <= 1 {
            return 0;
        }
        let clamped = seed_hit_fraction.clamp(0.0, 1.0);
        let scaled = (clamped * bin_count as f64).floor() as usize;
        scaled.min(bin_count.saturating_sub(1))
    }

    pub(super) fn select_rna_read_report_max_score_record_indices(
        report: &RnaReadInterpretationReport,
    ) -> Vec<usize> {
        let Some(max_score) = report
            .hits
            .iter()
            .map(|hit| hit.seed_hit_fraction)
            .max_by(|left, right| left.partial_cmp(right).unwrap_or(Ordering::Equal))
        else {
            return vec![];
        };
        report
            .hits
            .iter()
            .filter(|hit| (hit.seed_hit_fraction - max_score).abs() <= f64::EPSILON)
            .map(|hit| hit.record_index)
            .collect::<Vec<_>>()
    }

    pub(super) fn select_rna_read_report_rightmost_score_bin_record_indices(
        report: &RnaReadInterpretationReport,
        variant: RnaReadScoreDensityVariant,
        seed_filter_override: Option<&RnaReadSeedFilterConfig>,
    ) -> Vec<usize> {
        if report.hits.is_empty() {
            return vec![];
        }
        let bin_count = report.score_density_bins.len().max(40);
        let mut rightmost_idx = None;
        for hit in &report.hits {
            if !GentleEngine::rna_read_hit_matches_score_density_variant(
                hit,
                variant,
                seed_filter_override,
            ) {
                continue;
            }
            let idx = Self::rna_read_score_density_bin_index(hit.seed_hit_fraction, bin_count);
            rightmost_idx = Some(rightmost_idx.map_or(idx, |current: usize| current.max(idx)));
        }
        let Some(target_idx) = rightmost_idx else {
            return vec![];
        };
        report
            .hits
            .iter()
            .filter(|hit| {
                GentleEngine::rna_read_hit_matches_score_density_variant(
                    hit,
                    variant,
                    seed_filter_override,
                ) && Self::rna_read_score_density_bin_index(hit.seed_hit_fraction, bin_count)
                    == target_idx
            })
            .map(|hit| hit.record_index)
            .collect::<Vec<_>>()
    }

    pub(super) fn select_rna_read_report_score_bin_record_indices(
        report: &RnaReadInterpretationReport,
        target_idx: usize,
        bin_count: usize,
        variant: RnaReadScoreDensityVariant,
        seed_filter_override: Option<&RnaReadSeedFilterConfig>,
    ) -> Vec<usize> {
        if report.hits.is_empty() {
            return vec![];
        }
        let bin_count = bin_count.max(1);
        let target_idx = target_idx.min(bin_count.saturating_sub(1));
        report
            .hits
            .iter()
            .filter(|hit| {
                GentleEngine::rna_read_hit_matches_score_density_variant(
                    hit,
                    variant,
                    seed_filter_override,
                ) && Self::rna_read_score_density_bin_index(hit.seed_hit_fraction, bin_count)
                    == target_idx
            })
            .map(|hit| hit.record_index)
            .collect::<Vec<_>>()
    }

    pub(super) fn collect_rna_read_top_hit_previews_for_score_bin(
        report: &RnaReadInterpretationReport,
        target_idx: usize,
        bin_count: usize,
        variant: RnaReadScoreDensityVariant,
        seed_filter_override: Option<&RnaReadSeedFilterConfig>,
    ) -> Vec<RnaReadTopHitPreview> {
        let selected_record_indices = Self::select_rna_read_report_score_bin_record_indices(
            report,
            target_idx,
            bin_count,
            variant,
            seed_filter_override,
        )
        .into_iter()
        .collect::<BTreeSet<_>>();
        let mut rows = report
            .hits
            .iter()
            .filter(|hit| selected_record_indices.contains(&hit.record_index))
            .map(Self::rna_read_top_hit_preview_from_hit)
            .collect::<Vec<_>>();
        Self::sort_rna_read_top_hit_previews_by_phase1_score(&mut rows);
        rows
    }

    pub(super) fn format_rna_read_alignment_selection_summary(
        report: &RnaReadInterpretationReport,
        selection: RnaReadHitSelection,
    ) -> String {
        let retained_total = report.hits.len();
        let seed_passed_retained = report
            .hits
            .iter()
            .filter(|hit| hit.passed_seed_filter)
            .count();
        let raw_min_hit_retained = report
            .hits
            .iter()
            .filter(|hit| {
                hit.seed_hit_fraction + f64::EPSILON >= report.seed_filter.min_seed_hit_fraction
            })
            .count();
        let aligned_retained = report
            .hits
            .iter()
            .filter(|hit| hit.best_mapping.is_some())
            .count();
        match selection {
            RnaReadHitSelection::All => format!(
                "Phase-2 selection='all' aligns all retained saved-report rows ({retained_total}). Of those, {seed_passed_retained} currently pass the composite seed gate recorded for this report, and {raw_min_hit_retained} are at or above raw min_hit={:.2}. The red histogram line marks only that raw min_hit component, not the full composite gate.",
                report.seed_filter.min_seed_hit_fraction
            ),
            RnaReadHitSelection::SeedPassed => {
                if seed_passed_retained > 0 {
                    format!(
                        "Phase-2 selection='seed_passed' aligns the retained saved-report rows that currently pass the composite seed gate ({seed_passed_retained}). That count comes from `passed_seed_filter=yes` in the saved report, so it is usually stricter than the raw red min_hit line alone. If that subset ever becomes empty, the engine falls back to retained rows at or above raw min_hit={:.2} ({raw_min_hit_retained}), and failing that, to the single best retained row.",
                        report.seed_filter.min_seed_hit_fraction
                    )
                } else {
                    format!(
                        "Phase-2 selection='seed_passed' currently has no retained saved-report rows with `passed_seed_filter=yes`. The engine will therefore fall back to retained rows at or above raw min_hit={:.2} ({raw_min_hit_retained}), and if that is still empty, to the single best retained row.",
                        report.seed_filter.min_seed_hit_fraction
                    )
                }
            }
            RnaReadHitSelection::Aligned => format!(
                "Phase-2 selection='already_aligned' aligns only retained rows that already have a stored phase-2 mapping from an earlier pass ({aligned_retained})."
            ),
        }
    }

    pub(super) fn sort_rna_read_top_hit_previews_by_phase1_score(
        rows: &mut [RnaReadTopHitPreview],
    ) {
        rows.sort_by(|left, right| {
            right
                .seed_hit_fraction
                .partial_cmp(&left.seed_hit_fraction)
                .unwrap_or(Ordering::Equal)
                .then_with(|| {
                    right
                        .weighted_seed_hit_fraction
                        .partial_cmp(&left.weighted_seed_hit_fraction)
                        .unwrap_or(Ordering::Equal)
                })
                .then(right.matched_kmers.cmp(&left.matched_kmers))
                .then(right.tested_kmers.cmp(&left.tested_kmers))
                .then(left.record_index.cmp(&right.record_index))
        });
    }

    pub(super) fn selected_rna_report_hits<'a>(
        report: &'a RnaReadInterpretationReport,
        selected_record_indices: &BTreeSet<usize>,
    ) -> Vec<&'a RnaReadInterpretationHit> {
        report
            .hits
            .iter()
            .filter(|hit| selected_record_indices.contains(&hit.record_index))
            .collect::<Vec<_>>()
    }

    pub(super) fn format_rna_read_hit_fasta_entry(hit: &RnaReadInterpretationHit) -> String {
        let header_id = hit.header_id.trim().replace(['\n', '\r', '\t'], " ");
        let gap_median = if hit.seed_transcript_gap_count == 0 {
            "na".to_string()
        } else {
            format!("{:.2}", hit.seed_median_transcript_gap)
        };
        let mapping_summary = if let Some(best) = &hit.best_mapping {
            format!(
                "{}_tx={}_strand={}_t={}-{}_id={:.3}_cov={:.3}_score={}_sec={}",
                best.alignment_mode.as_str(),
                if best.transcript_id.is_empty() {
                    "unknown"
                } else {
                    best.transcript_id.as_str()
                },
                if best.strand.is_empty() {
                    "na"
                } else {
                    best.strand.as_str()
                },
                best.target_start_1based,
                best.target_end_1based,
                best.identity_fraction,
                best.query_coverage_fraction,
                best.score,
                hit.secondary_mappings.len()
            )
        } else {
            "none".to_string()
        };
        format!(
            ">{header_id} record_index={} score={:.3} wscore={:.4} wsupport={:.2} gap_med={} gap_n={} chain={:.2}/{} tx={} class={} origin_conf={:.3} strand_conf={:.3} strand={} opp={} ambig={} matched/tested={}/{} pass={} rc={} msa_eligible={} msa_reason={} align={} len={}\n{}",
            hit.record_index + 1,
            hit.seed_hit_fraction,
            hit.weighted_seed_hit_fraction,
            hit.weighted_matched_kmers,
            gap_median,
            hit.seed_transcript_gap_count,
            hit.seed_chain_support_fraction,
            hit.seed_chain_support_kmers,
            if hit.seed_chain_transcript_id.is_empty() {
                "none"
            } else {
                hit.seed_chain_transcript_id.as_str()
            },
            hit.origin_class.as_str(),
            hit.origin_confidence,
            hit.strand_confidence,
            if hit.strand_diagnostics.selected_strand.is_empty() {
                "na"
            } else {
                hit.strand_diagnostics.selected_strand.as_str()
            },
            hit.strand_diagnostics.competing_opposite_strand,
            hit.strand_diagnostics.ambiguous_near_tie,
            hit.matched_kmers,
            hit.tested_kmers,
            hit.passed_seed_filter,
            hit.reverse_complement_applied,
            hit.msa_eligible,
            hit.msa_eligibility_reason.trim(),
            mapping_summary,
            hit.read_length_bp,
            hit.sequence,
        )
    }

    pub(super) fn copy_rna_report_hits_as_fasta(
        &mut self,
        ui: &egui::Ui,
        hits: &[&RnaReadInterpretationHit],
        source_label: &str,
    ) {
        if hits.is_empty() {
            self.op_status = format!("No reads available to copy from {source_label}");
            return;
        }
        let fasta = hits
            .iter()
            .map(|hit| Self::format_rna_read_hit_fasta_entry(hit))
            .collect::<Vec<_>>()
            .join("\n");
        ui.ctx().copy_text(fasta);
        self.op_status = format!(
            "Copied {} read sequence(s) as FASTA from {source_label}",
            hits.len()
        );
    }

    pub(super) fn materialize_rna_read_report_hits(
        &mut self,
        selected_record_indices: Vec<usize>,
        source_label: &str,
    ) {
        if selected_record_indices.is_empty() {
            self.op_status = format!("No reads selected to materialize from {source_label}");
            return;
        }
        let report_id = self.rna_reads_ui.report_id.trim().to_string();
        if report_id.is_empty() {
            self.op_status = "Load/save a Report ID before materializing RNA-read hits".to_string();
            return;
        }
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let mut guard = match engine.write() {
            Ok(guard) => guard,
            Err(_) => {
                self.op_status =
                    "Engine lock poisoned while materializing RNA-read hit sequences".to_string();
                return;
            }
        };
        let output_prefix = Some(format!(
            "{}_reads",
            Self::sanitize_workflow_run_id_component(&report_id)
        ));
        let result = match guard.apply(Operation::MaterializeRnaReadHitSequences {
            report_id: report_id.clone(),
            selection: RnaReadHitSelection::All,
            selected_record_indices,
            output_prefix,
        }) {
            Ok(result) => result,
            Err(error) => {
                self.op_status = error.message.clone();
                self.op_error_popup = Some(error.message);
                return;
            }
        };
        self.last_created_seq_ids = result.created_seq_ids.clone();
        self.export_pool_inputs_text = self.last_created_seq_ids.join(", ");
        self.op_status = format!(
            "Materialized {} RNA-read hit sequence(s) from {source_label}",
            self.last_created_seq_ids.len()
        );
    }

    pub(super) fn materialize_selected_rna_read_report_hits(&mut self) {
        let selected_indices = self
            .rna_seed_selected_record_indices
            .iter()
            .copied()
            .collect::<Vec<_>>();
        self.materialize_rna_read_report_hits(selected_indices, "selected report reads");
    }

    pub(super) fn materialize_highlighted_rna_read_report_hit(&mut self) {
        let Some(record_index) = self.rna_seed_highlight_record_index else {
            self.op_status = "Highlight one read before materializing it".to_string();
            return;
        };
        self.materialize_rna_read_report_hits(vec![record_index], "highlighted report read");
    }
}
