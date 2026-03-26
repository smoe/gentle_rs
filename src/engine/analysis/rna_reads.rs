//! RNA-read report storage, parsing, scoring, and export helper routines.
//!
//! The RNA-read pipeline is large enough to deserve its own implementation
//! slice, but it still extends the shared engine contract rather than creating
//! a separate subsystem.

use super::*;

#[derive(Debug, Clone)]
struct RnaReadAlignmentScatterPoint {
    rank: RnaReadRetentionRank,
    record_index: usize,
    header_id: String,
    transcript_id: String,
    identity_fraction: f64,
    query_coverage_fraction: f64,
    score: isize,
}

impl GentleEngine {
    pub(super) fn read_rna_read_report_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> RnaReadReportStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<RnaReadReportStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = RNA_READ_REPORTS_SCHEMA.to_string();
        }
        store
    }

    pub(super) fn read_rna_read_report_store(&self) -> RnaReadReportStore {
        Self::read_rna_read_report_store_from_metadata(
            self.state.metadata.get(RNA_READ_REPORTS_METADATA_KEY),
        )
    }

    pub(super) fn write_rna_read_report_store(
        &mut self,
        mut store: RnaReadReportStore,
    ) -> Result<(), EngineError> {
        if store.reports.is_empty() {
            self.state.metadata.remove(RNA_READ_REPORTS_METADATA_KEY);
            return Ok(());
        }
        store.schema = RNA_READ_REPORTS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize RNA-read report metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(RNA_READ_REPORTS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub(super) fn upsert_rna_read_report(
        &mut self,
        report: RnaReadInterpretationReport,
    ) -> Result<(), EngineError> {
        let mut store = self.read_rna_read_report_store();
        store.reports.insert(report.report_id.clone(), report);
        self.write_rna_read_report_store(store)
    }

    pub(super) fn normalize_rna_read_report_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "report_id cannot be empty".to_string(),
            });
        }
        let mut out = String::with_capacity(trimmed.len());
        for ch in trimmed.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                out.push(ch);
            } else {
                out.push('_');
            }
        }
        if out.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "report_id must contain at least one ASCII letter, digit, '-', '_' or '.'"
                    .to_string(),
            });
        }
        Ok(out)
    }

    pub(super) fn normalize_rna_read_checkpoint_path(raw: Option<&str>) -> Option<String> {
        raw.map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
    }

    pub(super) fn read_rna_read_interpret_checkpoint(
        path: &str,
    ) -> Result<RnaReadInterpretCheckpoint, EngineError> {
        let text = std::fs::read_to_string(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read RNA-read checkpoint '{}': {e}", path),
        })?;
        let mut checkpoint =
            serde_json::from_str::<RnaReadInterpretCheckpoint>(&text).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not parse RNA-read checkpoint '{}': {e}", path),
            })?;
        if checkpoint.schema.trim().is_empty() {
            checkpoint.schema = RNA_READ_CHECKPOINT_SCHEMA.to_string();
        }
        if checkpoint.schema != RNA_READ_CHECKPOINT_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported RNA-read checkpoint schema '{}' in '{}'",
                    checkpoint.schema, path
                ),
            });
        }
        Ok(checkpoint)
    }

    pub(super) fn write_rna_read_interpret_checkpoint(
        path: &str,
        checkpoint: &RnaReadInterpretCheckpoint,
    ) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(checkpoint).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize RNA-read checkpoint for '{}': {e}",
                path
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write RNA-read checkpoint '{}': {e}", path),
        })
    }

    pub fn list_rna_read_reports(
        &self,
        seq_id_filter: Option<&str>,
    ) -> Vec<RnaReadInterpretationReportSummary> {
        let filter = seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let mut rows = self
            .read_rna_read_report_store()
            .reports
            .values()
            .filter(|report| {
                filter.as_ref().is_none_or(|needle| {
                    report
                        .seq_id
                        .to_ascii_lowercase()
                        .eq_ignore_ascii_case(needle)
                })
            })
            .map(|report| RnaReadInterpretationReportSummary {
                report_id: report.report_id.clone(),
                report_mode: report.report_mode,
                seq_id: report.seq_id.clone(),
                generated_at_unix_ms: report.generated_at_unix_ms,
                profile: report.profile,
                input_path: report.input_path.clone(),
                input_format: report.input_format,
                seed_feature_id: report.seed_feature_id,
                scope: report.scope,
                origin_mode: report.origin_mode,
                target_gene_count: report.target_gene_ids.len(),
                roi_seed_capture_enabled: report.roi_seed_capture_enabled,
                read_count_total: report.read_count_total,
                read_count_seed_passed: report.read_count_seed_passed,
                read_count_aligned: report.read_count_aligned,
                retained_count_msa_eligible: report.retained_count_msa_eligible,
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.seq_id
                .to_ascii_lowercase()
                .cmp(&right.seq_id.to_ascii_lowercase())
                .then(left.generated_at_unix_ms.cmp(&right.generated_at_unix_ms))
                .then(
                    left.report_id
                        .to_ascii_lowercase()
                        .cmp(&right.report_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub(crate) fn format_rna_read_report_summary_row(
        row: &RnaReadInterpretationReportSummary,
    ) -> String {
        let seed_pass_pct = if row.read_count_total == 0 {
            0.0
        } else {
            (row.read_count_seed_passed as f64 / row.read_count_total as f64) * 100.0
        };
        format!(
            "{} seq={} mode={} origin={} targets={} roi_capture={} reads={} seed_passed={} ({:.2}%) aligned={} msa_eligible(retained)={}",
            row.report_id,
            row.seq_id,
            row.report_mode.as_str(),
            row.origin_mode.as_str(),
            row.target_gene_count,
            row.roi_seed_capture_enabled,
            row.read_count_total,
            row.read_count_seed_passed,
            seed_pass_pct,
            row.read_count_aligned,
            row.retained_count_msa_eligible
        )
    }

    pub(crate) fn format_rna_read_report_detail_summary(
        report: &RnaReadInterpretationReport,
    ) -> String {
        let seed_pass_pct = if report.read_count_total == 0 {
            0.0
        } else {
            (report.read_count_seed_passed as f64 / report.read_count_total as f64) * 100.0
        };
        format!(
            "{} seq={} profile={} mode={} origin={} targets={} roi_capture={} reads={} seed_passed={} ({:.2}%) aligned={} msa_eligible(retained)={}",
            report.report_id,
            report.seq_id,
            report.profile.as_str(),
            report.report_mode.as_str(),
            report.origin_mode.as_str(),
            report.target_gene_ids.len(),
            report.roi_seed_capture_enabled,
            report.read_count_total,
            report.read_count_seed_passed,
            seed_pass_pct,
            report.read_count_aligned,
            report.retained_count_msa_eligible
        )
    }

    pub fn get_rna_read_report(
        &self,
        report_id: &str,
    ) -> Result<RnaReadInterpretationReport, EngineError> {
        let report_id = Self::normalize_rna_read_report_id(report_id)?;
        self.read_rna_read_report_store()
            .reports
            .get(report_id.as_str())
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("RNA-read report '{}' not found", report_id),
            })
    }

    pub fn export_rna_read_report(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<RnaReadInterpretationReport, EngineError> {
        let report = self.get_rna_read_report(report_id)?;
        let text = serde_json::to_string_pretty(&report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize RNA-read report '{}': {e}",
                report.report_id
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write RNA-read report to '{path}': {e}"),
        })?;
        Ok(report)
    }

    pub fn export_rna_read_hits_fasta(
        &self,
        report_id: &str,
        path: &str,
        selection: RnaReadHitSelection,
        selected_record_indices: &[usize],
        subset_spec: Option<&str>,
    ) -> Result<usize, EngineError> {
        let report = self.get_rna_read_report(report_id)?;
        let mut writer = BufWriter::new(File::create(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not create RNA-read FASTA output '{}': {e}", path),
        })?);
        let explicit_record_filter = selected_record_indices
            .iter()
            .copied()
            .collect::<HashSet<_>>();
        let mut written = 0usize;
        for hit in &report.hits {
            if !Self::include_rna_read_hit_by_selection_and_indices(
                hit,
                selection,
                &explicit_record_filter,
            ) {
                continue;
            }
            // Report payload keeps detailed mapping metrics while FASTA headers provide compact
            // human-readable context for quick triage.
            let mapping_header = if let Some(best) = &hit.best_mapping {
                format!(
                    "best={} mode={} strand={} target={}..{} identity={:.3} coverage={:.3}",
                    best.transcript_id,
                    best.alignment_mode.as_str(),
                    best.strand,
                    best.target_start_1based,
                    best.target_end_1based,
                    best.identity_fraction,
                    best.query_coverage_fraction
                )
            } else {
                "best=none".to_string()
            };
            let seed_gap_median = if hit.seed_transcript_gap_count == 0 {
                "na".to_string()
            } else {
                format!("{:.2}", hit.seed_median_transcript_gap)
            };
            let header = format!(
                "{} record_index={} byte_offset={} seed_hit_fraction={:.3} seed_gap_median={} seed_gap_count={} seed_chain_support={} seed_chain_frac={:.3} seed_chain_tx={} origin_class={} origin_conf={:.3} strand_conf={:.3} perfect={} rc_applied={} msa_eligible={} msa_reason={} exon_path_tx={} exon_path={} exon_transitions={}/{} subset_spec={} {}",
                hit.header_id,
                hit.record_index,
                hit.source_byte_offset,
                hit.seed_hit_fraction,
                seed_gap_median,
                hit.seed_transcript_gap_count,
                hit.seed_chain_support_kmers,
                hit.seed_chain_support_fraction,
                if hit.seed_chain_transcript_id.is_empty() {
                    "none"
                } else {
                    hit.seed_chain_transcript_id.as_str()
                },
                hit.origin_class.as_str(),
                hit.origin_confidence,
                hit.strand_confidence,
                hit.perfect_seed_match,
                hit.reverse_complement_applied,
                hit.msa_eligible,
                Self::sanitize_tsv_cell(&hit.msa_eligibility_reason),
                if hit.exon_path_transcript_id.is_empty() {
                    "none"
                } else {
                    hit.exon_path_transcript_id.as_str()
                },
                if hit.exon_path.is_empty() {
                    "none"
                } else {
                    hit.exon_path.as_str()
                },
                hit.exon_transitions_confirmed,
                hit.exon_transitions_total,
                Self::format_subset_spec_for_metadata(subset_spec),
                mapping_header
            );
            writeln!(writer, ">{header}").map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write FASTA header to '{}': {e}", path),
            })?;
            writeln!(writer, "{}", hit.sequence).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write FASTA sequence to '{}': {e}", path),
            })?;
            written += 1;
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not flush FASTA output '{}': {e}", path),
        })?;
        Ok(written)
    }

    pub(super) fn include_rna_read_hit_by_selection(
        hit: &RnaReadInterpretationHit,
        selection: RnaReadHitSelection,
    ) -> bool {
        match selection {
            RnaReadHitSelection::All => true,
            RnaReadHitSelection::SeedPassed => hit.passed_seed_filter,
            RnaReadHitSelection::Aligned => hit.best_mapping.is_some(),
        }
    }

    fn include_rna_read_hit_by_selection_and_indices(
        hit: &RnaReadInterpretationHit,
        selection: RnaReadHitSelection,
        explicit_record_filter: &HashSet<usize>,
    ) -> bool {
        if !explicit_record_filter.is_empty() {
            explicit_record_filter.contains(&hit.record_index)
        } else {
            Self::include_rna_read_hit_by_selection(hit, selection)
        }
    }

    fn primary_phase1_transcript_id(hit: &RnaReadInterpretationHit) -> &str {
        if !hit.exon_path_transcript_id.trim().is_empty() {
            hit.exon_path_transcript_id.as_str()
        } else {
            hit.seed_chain_transcript_id.as_str()
        }
    }

    fn alignment_effect_for_hit(
        hit: &RnaReadInterpretationHit,
        mapping: &RnaReadMappingHit,
    ) -> RnaReadAlignmentEffect {
        let phase1_transcript_id = Self::primary_phase1_transcript_id(hit).trim();
        if phase1_transcript_id.is_empty() {
            RnaReadAlignmentEffect::AlignedWithoutPhase1Assignment
        } else if phase1_transcript_id == mapping.transcript_id.trim() {
            RnaReadAlignmentEffect::ConfirmedAssignment
        } else {
            RnaReadAlignmentEffect::ReassignedTranscript
        }
    }

    fn normalize_rna_read_alignment_inspection_subset_spec(
        subset_spec: Option<RnaReadAlignmentInspectionSubsetSpec>,
    ) -> RnaReadAlignmentInspectionSubsetSpec {
        let mut subset_spec = subset_spec.unwrap_or_default();
        subset_spec.search = subset_spec.search.trim().to_string();
        subset_spec.selected_record_indices.sort_unstable();
        subset_spec.selected_record_indices.dedup();
        subset_spec
    }

    fn rna_read_alignment_effect_search_label(effect: RnaReadAlignmentEffect) -> &'static str {
        match effect {
            RnaReadAlignmentEffect::ConfirmedAssignment => "confirmed",
            RnaReadAlignmentEffect::ReassignedTranscript => "reassigned",
            RnaReadAlignmentEffect::AlignedWithoutPhase1Assignment => "aligned_without_phase1",
        }
    }

    fn rna_read_alignment_inspection_matches_filter(
        row: &RnaReadAlignmentInspectionRow,
        effect_filter: RnaReadAlignmentInspectionEffectFilter,
        selected_record_indices: &[usize],
    ) -> bool {
        match effect_filter {
            RnaReadAlignmentInspectionEffectFilter::AllAligned => true,
            RnaReadAlignmentInspectionEffectFilter::ConfirmedOnly => {
                row.alignment_effect == RnaReadAlignmentEffect::ConfirmedAssignment
            }
            RnaReadAlignmentInspectionEffectFilter::DisagreementOnly => {
                row.alignment_effect != RnaReadAlignmentEffect::ConfirmedAssignment
            }
            RnaReadAlignmentInspectionEffectFilter::ReassignedOnly => {
                row.alignment_effect == RnaReadAlignmentEffect::ReassignedTranscript
            }
            RnaReadAlignmentInspectionEffectFilter::NoPhase1Only => {
                row.alignment_effect == RnaReadAlignmentEffect::AlignedWithoutPhase1Assignment
            }
            RnaReadAlignmentInspectionEffectFilter::SelectedOnly => selected_record_indices
                .binary_search(&row.record_index)
                .is_ok(),
        }
    }

    fn rna_read_alignment_inspection_matches_search(
        row: &RnaReadAlignmentInspectionRow,
        search: &str,
    ) -> bool {
        let needle = search.trim().to_ascii_lowercase();
        if needle.is_empty() {
            return true;
        }
        let effect_label = Self::rna_read_alignment_effect_search_label(row.alignment_effect);
        let effect_schema = row.alignment_effect.as_str();
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
            effect_schema,
            record_index_label.as_str(),
            rank_label.as_str(),
        ]
        .iter()
        .any(|field| field.to_ascii_lowercase().contains(&needle))
    }

    fn compare_rna_read_alignment_inspection_rows(
        left: &RnaReadAlignmentInspectionRow,
        right: &RnaReadAlignmentInspectionRow,
        sort_key: RnaReadAlignmentInspectionSortKey,
    ) -> Ordering {
        match sort_key {
            RnaReadAlignmentInspectionSortKey::Rank => left
                .rank
                .cmp(&right.rank)
                .then_with(|| left.record_index.cmp(&right.record_index)),
            RnaReadAlignmentInspectionSortKey::Identity => right
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
            RnaReadAlignmentInspectionSortKey::Coverage => right
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
            RnaReadAlignmentInspectionSortKey::Score => right
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

    pub fn inspect_rna_read_alignments(
        &self,
        report_id: &str,
        selection: RnaReadHitSelection,
        limit: usize,
    ) -> Result<RnaReadAlignmentInspection, EngineError> {
        self.inspect_rna_read_alignments_with_subset(report_id, selection, limit, None)
    }

    pub fn inspect_rna_read_alignments_with_subset(
        &self,
        report_id: &str,
        selection: RnaReadHitSelection,
        limit: usize,
        subset_spec: Option<RnaReadAlignmentInspectionSubsetSpec>,
    ) -> Result<RnaReadAlignmentInspection, EngineError> {
        if limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "RNA-read alignment inspection requires --limit >= 1".to_string(),
            });
        }
        let subset_spec = Self::normalize_rna_read_alignment_inspection_subset_spec(subset_spec);
        let report = self.get_rna_read_report(report_id)?;
        let splicing_view = self
            .build_splicing_expert_view(&report.seq_id, report.seed_feature_id, report.scope)
            .ok();
        let mut ranked_rows = report
            .hits
            .iter()
            .filter(|hit| Self::include_rna_read_hit_by_selection(hit, selection))
            .filter_map(|hit| {
                let mapping = hit.best_mapping.as_ref()?;
                let rank = Self::rna_read_retention_rank(hit);
                let (mapped_exon_support, mapped_junction_support) = splicing_view
                    .as_ref()
                    .map(|splicing| {
                        Self::collect_mapped_support_attribution_rows(mapping, splicing)
                    })
                    .unwrap_or_default();
                Some((
                    rank,
                    RnaReadAlignmentInspectionRow {
                        rank: 0,
                        record_index: hit.record_index,
                        header_id: hit.header_id.clone(),
                        phase1_primary_transcript_id: Self::primary_phase1_transcript_id(hit)
                            .to_string(),
                        seed_chain_transcript_id: hit.seed_chain_transcript_id.clone(),
                        exon_path_transcript_id: hit.exon_path_transcript_id.clone(),
                        exon_path: hit.exon_path.clone(),
                        exon_transitions_confirmed: hit.exon_transitions_confirmed,
                        exon_transitions_total: hit.exon_transitions_total,
                        selected_strand: hit.strand_diagnostics.selected_strand.clone(),
                        reverse_complement_applied: hit.reverse_complement_applied,
                        alignment_effect: Self::alignment_effect_for_hit(hit, mapping),
                        transcript_id: mapping.transcript_id.clone(),
                        transcript_label: mapping.transcript_label.clone(),
                        strand: mapping.strand.clone(),
                        alignment_mode: mapping.alignment_mode,
                        target_start_1based: mapping.target_start_1based,
                        target_end_1based: mapping.target_end_1based,
                        score: mapping.score,
                        identity_fraction: mapping.identity_fraction,
                        query_coverage_fraction: mapping.query_coverage_fraction,
                        secondary_mapping_count: hit.secondary_mappings.len(),
                        seed_hit_fraction: hit.seed_hit_fraction,
                        weighted_seed_hit_fraction: hit.weighted_seed_hit_fraction,
                        passed_seed_filter: hit.passed_seed_filter,
                        msa_eligible: hit.msa_eligible,
                        origin_class: hit.origin_class,
                        mapped_exon_support,
                        mapped_junction_support,
                    },
                ))
            })
            .collect::<Vec<_>>();
        ranked_rows.sort_by(|left, right| right.0.cmp(&left.0));
        let aligned_count = ranked_rows.len();
        let mut rows = ranked_rows
            .into_iter()
            .enumerate()
            .map(|(idx, (_, mut row))| {
                row.rank = idx + 1;
                row
            })
            .collect::<Vec<_>>();
        rows.retain(|row| {
            Self::rna_read_alignment_inspection_matches_filter(
                row,
                subset_spec.effect_filter,
                &subset_spec.selected_record_indices,
            ) && Self::rna_read_alignment_inspection_matches_search(row, &subset_spec.search)
        });
        rows.sort_by(|left, right| {
            Self::compare_rna_read_alignment_inspection_rows(left, right, subset_spec.sort_key)
        });
        let subset_match_count = rows.len();
        rows.truncate(limit);
        Ok(RnaReadAlignmentInspection {
            schema: RNA_READ_ALIGNMENT_INSPECTION_SCHEMA.to_string(),
            report_id: report.report_id,
            seq_id: report.seq_id,
            selection,
            row_count: rows.len(),
            aligned_count,
            subset_match_count,
            limit,
            subset_spec,
            align_min_identity_fraction: report.align_config.min_identity_fraction,
            max_secondary_mappings: report.align_config.max_secondary_mappings,
            rows,
        })
    }

    pub fn export_rna_read_alignments_tsv(
        &self,
        report_id: &str,
        path: &str,
        selection: RnaReadHitSelection,
        limit: Option<usize>,
        selected_record_indices: &[usize],
        subset_spec: Option<&str>,
    ) -> Result<RnaReadAlignmentTsvExport, EngineError> {
        let path = path.trim();
        if path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "RNA-read alignment TSV export requires non-empty path".to_string(),
            });
        }
        if let Some(value) = limit {
            if value == 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "RNA-read alignment TSV export requires --limit >= 1".to_string(),
                });
            }
        }
        let report = self.get_rna_read_report(report_id)?;
        let mut inspection = self.inspect_rna_read_alignments(report_id, selection, usize::MAX)?;
        let explicit_record_filter = selected_record_indices
            .iter()
            .copied()
            .collect::<HashSet<_>>();
        if !explicit_record_filter.is_empty() {
            inspection
                .rows
                .retain(|row| explicit_record_filter.contains(&row.record_index));
            inspection.row_count = inspection.rows.len();
            inspection.aligned_count = inspection.row_count;
        }
        if let Some(value) = limit {
            inspection.rows.truncate(value);
            inspection.row_count = inspection.rows.len();
        }
        let file = File::create(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create RNA-read alignment TSV export '{}': {e}",
                path
            ),
        })?;
        let mut writer = BufWriter::new(file);
        for line in Self::rna_read_alignment_tsv_metadata_lines(
            &report,
            &inspection,
            selection,
            limit,
            selected_record_indices,
            subset_spec,
        ) {
            writeln!(writer, "{line}").map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not write RNA-read alignment TSV metadata to '{}': {e}",
                    path
                ),
            })?;
        }
        writeln!(
            writer,
            "report_id\tseq_id\trank\trecord_index\theader_id\tphase1_primary_transcript_id\tseed_chain_transcript_id\texon_path_transcript_id\texon_path\texon_transitions_confirmed\texon_transitions_total\tselected_strand\treverse_complement_applied\talignment_effect\ttranscript_id\ttranscript_label\tstrand\talignment_mode\ttarget_start_1based\ttarget_end_1based\tscore\tidentity_fraction\tquery_coverage_fraction\tsecondary_mapping_count\tseed_hit_fraction\tweighted_seed_hit_fraction\tpassed_seed_filter\tmsa_eligible\torigin_class\tmapped_exon_support\tmapped_junction_support"
        )
        .map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write RNA-read alignment TSV header to '{}': {e}",
                path
            ),
        })?;
        for row in &inspection.rows {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}",
                Self::sanitize_tsv_cell(&inspection.report_id),
                Self::sanitize_tsv_cell(&inspection.seq_id),
                row.rank,
                row.record_index,
                Self::sanitize_tsv_cell(&row.header_id),
                Self::sanitize_tsv_cell(&row.phase1_primary_transcript_id),
                Self::sanitize_tsv_cell(&row.seed_chain_transcript_id),
                Self::sanitize_tsv_cell(&row.exon_path_transcript_id),
                Self::sanitize_tsv_cell(&row.exon_path),
                row.exon_transitions_confirmed,
                row.exon_transitions_total,
                Self::sanitize_tsv_cell(&row.selected_strand),
                row.reverse_complement_applied,
                row.alignment_effect.as_str(),
                Self::sanitize_tsv_cell(&row.transcript_id),
                Self::sanitize_tsv_cell(&row.transcript_label),
                Self::sanitize_tsv_cell(&row.strand),
                row.alignment_mode.as_str(),
                row.target_start_1based,
                row.target_end_1based,
                row.score,
                row.identity_fraction,
                row.query_coverage_fraction,
                row.secondary_mapping_count,
                row.seed_hit_fraction,
                row.weighted_seed_hit_fraction,
                row.passed_seed_filter,
                row.msa_eligible,
                row.origin_class.as_str(),
                Self::sanitize_tsv_cell(&Self::format_mapped_exon_support_for_tsv(
                    &row.mapped_exon_support
                )),
                Self::sanitize_tsv_cell(&Self::format_mapped_junction_support_for_tsv(
                    &row.mapped_junction_support
                )),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not write RNA-read alignment TSV row to '{}': {e}",
                    path
                ),
            })?;
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not flush RNA-read alignment TSV '{}': {e}", path),
        })?;
        Ok(RnaReadAlignmentTsvExport {
            schema: RNA_READ_ALIGNMENT_TSV_EXPORT_SCHEMA.to_string(),
            path: path.to_string(),
            report_id: inspection.report_id,
            selection,
            row_count: inspection.row_count,
            aligned_count: inspection.aligned_count,
            limit,
        })
    }

    pub(super) fn alignment_scatter_color_hex(
        score: isize,
        min_score: isize,
        max_score: isize,
    ) -> String {
        let t = if max_score <= min_score {
            0.5
        } else {
            ((score - min_score) as f64 / (max_score - min_score) as f64).clamp(0.0, 1.0)
        };
        let low = (37u8, 99u8, 235u8);
        let high = (220u8, 38u8, 38u8);
        let mix = |a: u8, b: u8| -> u8 { (a as f64 + (b as f64 - a as f64) * t).round() as u8 };
        let r = mix(low.0, high.0);
        let g = mix(low.1, high.1);
        let b = mix(low.2, high.2);
        format!("#{r:02x}{g:02x}{b:02x}")
    }

    fn render_rna_read_alignment_dotplot_svg_text(
        report: &RnaReadInterpretationReport,
        selection: RnaReadHitSelection,
        points: &[RnaReadAlignmentScatterPoint],
        total_points: usize,
        max_points: usize,
        min_score: isize,
        max_score: isize,
    ) -> String {
        let width = 960.0f64;
        let height = 620.0f64;
        let margin_left = 72.0f64;
        let margin_right = 28.0f64;
        let margin_top = 58.0f64;
        let margin_bottom = 72.0f64;
        let chart_left = margin_left;
        let chart_right = width - margin_right;
        let chart_top = margin_top;
        let chart_bottom = height - margin_bottom;
        let chart_width = (chart_right - chart_left).max(1.0);
        let chart_height = (chart_bottom - chart_top).max(1.0);
        let identity_threshold = report.align_config.min_identity_fraction.clamp(0.0, 1.0);

        let mut svg = String::new();
        svg.push_str(&format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width:.0}\" height=\"{height:.0}\" viewBox=\"0 0 {width:.0} {height:.0}\">\n"
        ));
        svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n");
        svg.push_str(&format!(
            "<text x=\"{x:.1}\" y=\"24\" font-family=\"monospace\" font-size=\"15\" fill=\"#111827\">{title}</text>\n",
            x = chart_left,
            title = Self::escape_svg_text(&format!(
                "RNA-read alignment dotplot (coverage vs identity): {}",
                report.report_id
            ))
        ));
        svg.push_str(&format!(
            "<text x=\"{x:.1}\" y=\"40\" font-family=\"monospace\" font-size=\"11\" fill=\"#4b5563\">{subtitle}</text>\n",
            x = chart_left,
            subtitle = Self::escape_svg_text(&format!(
                "selection={} rendered_points={} total_points={} max_points={} min_identity={:.3} score_range=[{}, {}]",
                selection.as_str(),
                points.len(),
                total_points,
                max_points,
                identity_threshold,
                min_score,
                max_score
            ))
        ));

        for tick in [0.0f64, 0.25, 0.50, 0.75, 1.0] {
            let x = chart_left + tick * chart_width;
            let y = chart_bottom - tick * chart_height;
            svg.push_str(&format!(
                "<line x1=\"{x:.2}\" y1=\"{top:.2}\" x2=\"{x:.2}\" y2=\"{bottom:.2}\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n",
                top = chart_top,
                bottom = chart_bottom
            ));
            svg.push_str(&format!(
                "<line x1=\"{left:.2}\" y1=\"{y:.2}\" x2=\"{right:.2}\" y2=\"{y:.2}\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n",
                left = chart_left,
                right = chart_right
            ));
            svg.push_str(&format!(
                "<text x=\"{x:.2}\" y=\"{label_y:.2}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"10\" fill=\"#6b7280\">{tick:.2}</text>\n",
                label_y = chart_bottom + 16.0
            ));
            svg.push_str(&format!(
                "<text x=\"{label_x:.2}\" y=\"{y_plus:.2}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#6b7280\">{tick:.2}</text>\n",
                label_x = chart_left - 8.0,
                y_plus = y + 3.0
            ));
        }

        svg.push_str(&format!(
            "<rect x=\"{x:.2}\" y=\"{y:.2}\" width=\"{w:.2}\" height=\"{h:.2}\" fill=\"none\" stroke=\"#9ca3af\" stroke-width=\"1.2\"/>\n",
            x = chart_left,
            y = chart_top,
            w = chart_width,
            h = chart_height
        ));
        svg.push_str(&format!(
            "<line x1=\"{left:.2}\" y1=\"{y:.2}\" x2=\"{right:.2}\" y2=\"{y:.2}\" stroke=\"#ef4444\" stroke-width=\"1\" stroke-dasharray=\"4 3\"/>\n",
            left = chart_left,
            right = chart_right,
            y = chart_bottom - identity_threshold * chart_height
        ));
        svg.push_str(&format!(
            "<text x=\"{x:.2}\" y=\"{y:.2}\" font-family=\"monospace\" font-size=\"10\" fill=\"#b91c1c\">min identity={:.3}</text>\n",
            identity_threshold,
            x = chart_right - 180.0,
            y = chart_bottom - identity_threshold * chart_height - 6.0
        ));
        svg.push_str(&format!(
            "<text x=\"{x:.2}\" y=\"{y:.2}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"11\" fill=\"#111827\">query coverage fraction</text>\n",
            x = (chart_left + chart_right) * 0.5,
            y = height - 24.0
        ));
        svg.push_str(&format!(
            "<text x=\"18\" y=\"{y:.2}\" transform=\"rotate(-90 18 {y:.2})\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"11\" fill=\"#111827\">identity fraction</text>\n",
            y = (chart_top + chart_bottom) * 0.5
        ));

        for point in points {
            let coverage = point.query_coverage_fraction.clamp(0.0, 1.0);
            let identity = point.identity_fraction.clamp(0.0, 1.0);
            let x = chart_left + coverage * chart_width;
            let y = chart_bottom - identity * chart_height;
            let color = Self::alignment_scatter_color_hex(point.score, min_score, max_score);
            svg.push_str(&format!(
                "<circle cx=\"{x:.2}\" cy=\"{y:.2}\" r=\"2.3\" fill=\"{color}\" fill-opacity=\"0.72\" stroke=\"#0f172a\" stroke-opacity=\"0.20\" stroke-width=\"0.35\"><title>{tip}</title></circle>\n",
                tip = Self::escape_svg_text(&format!(
                    "record={} tx={} score={} identity={:.4} coverage={:.4} header={}",
                    point.record_index,
                    point.transcript_id,
                    point.score,
                    point.identity_fraction,
                    point.query_coverage_fraction,
                    point.header_id
                ))
            ));
        }

        svg.push_str("</svg>\n");
        svg
    }

    pub fn export_rna_read_alignment_dotplot_svg(
        &self,
        report_id: &str,
        path: &str,
        selection: RnaReadHitSelection,
        max_points: usize,
    ) -> Result<RnaReadAlignmentDotplotSvgExport, EngineError> {
        let report = self.get_rna_read_report(report_id)?;
        let path = path.trim();
        if path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "RNA-read alignment dotplot export requires non-empty path".to_string(),
            });
        }
        if max_points == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "RNA-read alignment dotplot export requires max_points >= 1".to_string(),
            });
        }

        let mut points = report
            .hits
            .iter()
            .filter(|hit| Self::include_rna_read_hit_by_selection(hit, selection))
            .filter_map(|hit| {
                let mapping = hit.best_mapping.as_ref()?;
                Some(RnaReadAlignmentScatterPoint {
                    rank: Self::rna_read_retention_rank(hit),
                    record_index: hit.record_index,
                    header_id: hit.header_id.clone(),
                    transcript_id: mapping.transcript_id.clone(),
                    identity_fraction: mapping.identity_fraction,
                    query_coverage_fraction: mapping.query_coverage_fraction,
                    score: mapping.score,
                })
            })
            .collect::<Vec<_>>();
        points.sort_by(|left, right| right.rank.cmp(&left.rank));
        let point_count = points.len();
        if points.len() > max_points {
            points.truncate(max_points);
        }
        let rendered_point_count = points.len();
        let min_score = points.iter().map(|point| point.score).min().unwrap_or(0);
        let max_score = points.iter().map(|point| point.score).max().unwrap_or(0);
        let svg = Self::render_rna_read_alignment_dotplot_svg_text(
            &report,
            selection,
            &points,
            point_count,
            max_points,
            min_score,
            max_score,
        );
        std::fs::write(path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write RNA-read alignment dotplot SVG '{}': {e}",
                path
            ),
        })?;
        Ok(RnaReadAlignmentDotplotSvgExport {
            schema: RNA_READ_ALIGNMENT_DOTPLOT_SVG_EXPORT_SCHEMA.to_string(),
            path: path.to_string(),
            report_id: report.report_id,
            selection,
            point_count,
            rendered_point_count,
            max_points,
            min_score,
            max_score,
        })
    }

    pub(super) fn parse_exon_path(path: &str) -> (Vec<usize>, Vec<char>) {
        if path.trim().is_empty() {
            return (vec![], vec![]);
        }
        let mut ordinals = Vec::<usize>::new();
        let mut transitions = Vec::<char>::new();
        let mut current = String::new();
        for ch in path.chars() {
            if ch.is_ascii_digit() {
                current.push(ch);
                continue;
            }
            if (ch == ':' || ch == '-') && !current.is_empty() {
                if let Ok(value) = current.parse::<usize>() {
                    ordinals.push(value);
                    transitions.push(ch);
                }
                current.clear();
            }
        }
        if !current.is_empty() {
            if let Ok(value) = current.parse::<usize>() {
                ordinals.push(value);
            }
        }
        while transitions.len() > ordinals.len().saturating_sub(1) {
            let _ = transitions.pop();
        }
        (ordinals, transitions)
    }

    pub fn export_rna_read_exon_paths_tsv(
        &self,
        report_id: &str,
        path: &str,
        selection: RnaReadHitSelection,
        selected_record_indices: &[usize],
        subset_spec: Option<&str>,
    ) -> Result<RnaReadExonPathsExport, EngineError> {
        let report = self.get_rna_read_report(report_id)?;
        let path = path.trim();
        if path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "RNA-read exon-path export requires non-empty path".to_string(),
            });
        }
        let file = File::create(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not create RNA-read exon-path export '{}': {e}", path),
        })?;
        let mut writer = BufWriter::new(file);
        for line in Self::rna_read_tsv_common_metadata_lines(
            &report,
            selection,
            selected_record_indices,
            subset_spec,
        ) {
            writeln!(writer, "{line}").map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not write RNA-read exon-path export metadata to '{}': {e}",
                    path
                ),
            })?;
        }
        writeln!(
            writer,
            "report_id\tseq_id\trecord_index\theader_id\tsource_byte_offset\tread_length_bp\tseed_hit_fraction\tweighted_seed_hit_fraction\tweighted_matched_kmers\tseed_chain_transcript_id\tseed_chain_support_kmers\tseed_chain_support_fraction\tseed_median_transcript_gap\tseed_transcript_gap_count\tmatched_kmers\ttested_kmers\tpassed_seed_filter\treverse_complement_applied\torigin_class\torigin_reason\torigin_confidence\tstrand_confidence\tmsa_eligible\tmsa_eligibility_reason\texon_path_transcript_id\texon_path\texon_transitions_confirmed\texon_transitions_total\tbest_transcript_id\tbest_alignment_mode\tbest_strand\tbest_target_start_1based\tbest_target_end_1based\tbest_identity_fraction\tbest_query_coverage_fraction"
        )
        .map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write RNA-read exon-path export header to '{}': {e}",
                path
            ),
        })?;
        let explicit_record_filter = selected_record_indices
            .iter()
            .copied()
            .collect::<HashSet<_>>();
        let mut row_count = 0usize;
        for hit in &report.hits {
            if !Self::include_rna_read_hit_by_selection_and_indices(
                hit,
                selection,
                &explicit_record_filter,
            ) {
                continue;
            }
            let (
                best_transcript_id,
                best_alignment_mode,
                best_strand,
                best_target_start_1based,
                best_target_end_1based,
                best_identity_fraction,
                best_query_coverage_fraction,
            ) = if let Some(best) = &hit.best_mapping {
                (
                    best.transcript_id.as_str(),
                    best.alignment_mode.as_str(),
                    best.strand.as_str(),
                    best.target_start_1based.to_string(),
                    best.target_end_1based.to_string(),
                    format!("{:.6}", best.identity_fraction),
                    format!("{:.6}", best.query_coverage_fraction),
                )
            } else {
                (
                    "none",
                    "none",
                    "none",
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                )
            };
            let row = [
                Self::sanitize_tsv_cell(&report.report_id),
                Self::sanitize_tsv_cell(&report.seq_id),
                hit.record_index.to_string(),
                Self::sanitize_tsv_cell(&hit.header_id),
                hit.source_byte_offset.to_string(),
                hit.read_length_bp.to_string(),
                format!("{:.6}", hit.seed_hit_fraction),
                format!("{:.6}", hit.weighted_seed_hit_fraction),
                format!("{:.6}", hit.weighted_matched_kmers),
                Self::sanitize_tsv_cell(&hit.seed_chain_transcript_id),
                hit.seed_chain_support_kmers.to_string(),
                format!("{:.6}", hit.seed_chain_support_fraction),
                format!("{:.6}", hit.seed_median_transcript_gap),
                hit.seed_transcript_gap_count.to_string(),
                hit.matched_kmers.to_string(),
                hit.tested_kmers.to_string(),
                hit.passed_seed_filter.to_string(),
                hit.reverse_complement_applied.to_string(),
                Self::sanitize_tsv_cell(hit.origin_class.as_str()),
                Self::sanitize_tsv_cell(&hit.origin_reason),
                format!("{:.6}", hit.origin_confidence),
                format!("{:.6}", hit.strand_confidence),
                hit.msa_eligible.to_string(),
                Self::sanitize_tsv_cell(&hit.msa_eligibility_reason),
                Self::sanitize_tsv_cell(&hit.exon_path_transcript_id),
                Self::sanitize_tsv_cell(&hit.exon_path),
                hit.exon_transitions_confirmed.to_string(),
                hit.exon_transitions_total.to_string(),
                Self::sanitize_tsv_cell(best_transcript_id),
                Self::sanitize_tsv_cell(best_alignment_mode),
                Self::sanitize_tsv_cell(best_strand),
                best_target_start_1based,
                best_target_end_1based,
                best_identity_fraction,
                best_query_coverage_fraction,
            ]
            .join("\t");
            writeln!(writer, "{row}").map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write RNA-read exon-path row to '{}': {e}", path),
            })?;
            row_count = row_count.saturating_add(1);
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not flush RNA-read exon-path export '{}': {e}", path),
        })?;
        Ok(RnaReadExonPathsExport {
            schema: RNA_READ_EXON_PATHS_EXPORT_SCHEMA.to_string(),
            path: path.to_string(),
            report_id: report.report_id,
            selection,
            row_count,
        })
    }

    pub fn export_rna_read_exon_abundance_tsv(
        &self,
        report_id: &str,
        path: &str,
        selection: RnaReadHitSelection,
        selected_record_indices: &[usize],
        subset_spec: Option<&str>,
    ) -> Result<RnaReadExonAbundanceExport, EngineError> {
        let report = self.get_rna_read_report(report_id)?;
        let path = path.trim();
        if path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "RNA-read exon-abundance export requires non-empty path".to_string(),
            });
        }
        let file = File::create(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create RNA-read exon-abundance export '{}': {e}",
                path
            ),
        })?;
        let mut writer = BufWriter::new(file);
        for line in Self::rna_read_tsv_common_metadata_lines(
            &report,
            selection,
            selected_record_indices,
            subset_spec,
        ) {
            writeln!(writer, "{line}").map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not write RNA-read exon-abundance metadata to '{}': {e}",
                    path
                ),
            })?;
        }
        writeln!(
            writer,
            "report_id\tseq_id\tselection\trow_kind\texon_ordinal\tfrom_exon_ordinal\tto_exon_ordinal\tsupport_read_count\tsupport_fraction\tconfirmed_read_count\tconfirmed_fraction\tunconfirmed_read_count\tunconfirmed_fraction"
        )
        .map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write RNA-read exon-abundance header to '{}': {e}",
                path
            ),
        })?;
        let mut exon_counts = BTreeMap::<usize, usize>::new();
        let mut transition_counts = BTreeMap::<(usize, usize), (usize, usize)>::new();
        let explicit_record_filter = selected_record_indices
            .iter()
            .copied()
            .collect::<HashSet<_>>();
        let mut selected_read_count = 0usize;
        for hit in &report.hits {
            if !Self::include_rna_read_hit_by_selection_and_indices(
                hit,
                selection,
                &explicit_record_filter,
            ) {
                continue;
            }
            selected_read_count = selected_read_count.saturating_add(1);
            let (ordinals, transitions) = Self::parse_exon_path(&hit.exon_path);
            if ordinals.is_empty() {
                continue;
            }
            let mut touched = BTreeSet::<usize>::new();
            for exon_ordinal in ordinals.iter().copied() {
                touched.insert(exon_ordinal);
            }
            for exon_ordinal in touched {
                *exon_counts.entry(exon_ordinal).or_insert(0) += 1;
            }
            for (idx, transition) in transitions.iter().enumerate() {
                if idx + 1 >= ordinals.len() {
                    break;
                }
                let pair = (ordinals[idx], ordinals[idx + 1]);
                let entry = transition_counts.entry(pair).or_insert((0, 0));
                entry.1 = entry.1.saturating_add(1);
                if *transition == ':' {
                    entry.0 = entry.0.saturating_add(1);
                }
            }
        }
        let denominator = selected_read_count.max(1) as f64;
        for (exon_ordinal, count) in &exon_counts {
            let support_fraction = *count as f64 / denominator;
            writeln!(
                writer,
                "{}\t{}\t{}\texon\t{}\t\t\t{}\t{:.6}\t\t\t\t",
                Self::sanitize_tsv_cell(&report.report_id),
                Self::sanitize_tsv_cell(&report.seq_id),
                selection.as_str(),
                exon_ordinal,
                count,
                support_fraction,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write exon-abundance row to '{}': {e}", path),
            })?;
        }
        for ((from_exon, to_exon), (confirmed, total)) in &transition_counts {
            let support_fraction = *total as f64 / denominator;
            let confirmed_fraction = *confirmed as f64 / denominator;
            let unconfirmed = total.saturating_sub(*confirmed);
            let unconfirmed_fraction = unconfirmed as f64 / denominator;
            writeln!(
                writer,
                "{}\t{}\t{}\ttransition\t\t{}\t{}\t{}\t{:.6}\t{}\t{:.6}\t{}\t{:.6}",
                Self::sanitize_tsv_cell(&report.report_id),
                Self::sanitize_tsv_cell(&report.seq_id),
                selection.as_str(),
                from_exon,
                to_exon,
                total,
                support_fraction,
                confirmed,
                confirmed_fraction,
                unconfirmed,
                unconfirmed_fraction,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not write transition-abundance row to '{}': {e}",
                    path
                ),
            })?;
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not flush RNA-read exon-abundance export '{}': {e}",
                path
            ),
        })?;
        Ok(RnaReadExonAbundanceExport {
            schema: RNA_READ_EXON_ABUNDANCE_EXPORT_SCHEMA.to_string(),
            path: path.to_string(),
            report_id: report.report_id,
            selection,
            selected_read_count,
            exon_row_count: exon_counts.len(),
            transition_row_count: transition_counts.len(),
        })
    }

    pub(super) fn derive_score_density_bins_from_report(
        report: &RnaReadInterpretationReport,
    ) -> (Vec<u64>, bool) {
        if !report.score_density_bins.is_empty() {
            return (report.score_density_bins.clone(), false);
        }
        let mut bins = vec![0u64; RNA_READ_SCORE_DENSITY_BIN_COUNT];
        for hit in &report.hits {
            let idx = Self::score_density_bin_index(hit.seed_hit_fraction);
            if let Some(bucket) = bins.get_mut(idx) {
                *bucket = bucket.saturating_add(1);
            }
        }
        (bins, true)
    }

    pub(super) fn escape_svg_text(raw: &str) -> String {
        raw.replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
            .replace('"', "&quot;")
            .replace('\'', "&apos;")
    }

    fn format_score_density_count_compact(count: u64) -> String {
        if count >= 1_000_000 {
            let value = count as f64 / 1_000_000.0;
            if value >= 10.0 {
                format!("{value:.0}M")
            } else {
                format!("{value:.1}M")
            }
        } else if count >= 1_000 {
            let value = count as f64 / 1_000.0;
            if value >= 10.0 {
                format!("{value:.0}k")
            } else {
                format!("{value:.1}k")
            }
        } else {
            count.to_string()
        }
    }

    pub(super) fn render_rna_read_score_density_svg_text(
        report: &RnaReadInterpretationReport,
        bins: &[u64],
        scale: RnaReadScoreDensityScale,
    ) -> String {
        let width = 960.0f64;
        let height = 300.0f64;
        let margin_left = 56.0f64;
        let margin_right = 24.0f64;
        let margin_top = 66.0f64;
        let margin_bottom = 40.0f64;
        let axis_bottom = height - margin_bottom;
        let plot_top = margin_top + 16.0;
        let plot_height = (axis_bottom - plot_top).max(1.0);
        let chart_left = margin_left;
        let chart_right = width - margin_right;
        let chart_width = (chart_right - chart_left).max(1.0);
        let max_count = bins.iter().copied().max().unwrap_or(0);
        let max_scaled = if scale == RnaReadScoreDensityScale::Log {
            (max_count as f64 + 1.0).ln().max(1.0)
        } else {
            (max_count.max(1)) as f64
        };
        let threshold = report.seed_filter.min_seed_hit_fraction.clamp(0.0, 1.0);
        let threshold_x = chart_left + chart_width * threshold;
        let bar_count = bins.len().max(1);
        let bin_width = chart_width / bar_count as f64;
        let scale_text = if scale == RnaReadScoreDensityScale::Log {
            "log(1+count)"
        } else {
            "linear count"
        };
        let title = format!("RNA-read seed-hit score density ({scale_text})");
        let subtitle = format!(
            "report={} | reads={} seed-passed={} aligned={}",
            report.report_id,
            report.read_count_total,
            report.read_count_seed_passed,
            report.read_count_aligned
        );
        let provenance = format!(
            "profile={} mode={} scope={} origin={} | {}",
            report.profile.as_str(),
            report.report_mode.as_str(),
            report.scope.as_str(),
            report.origin_mode.as_str(),
            Self::rna_read_seed_filter_summary(&report.seed_filter),
        );
        let bins_source = if report.score_density_bins.is_empty() {
            "score-density bins derived from retained hits"
        } else {
            "score-density bins stored in report"
        };

        let mut svg = String::new();
        svg.push_str(&format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width:.0}\" height=\"{height:.0}\" viewBox=\"0 0 {width:.0} {height:.0}\">\n"
        ));
        svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n");
        svg.push_str(&format!(
            "<text x=\"{x:.1}\" y=\"22\" font-family=\"monospace\" font-size=\"14\" fill=\"#111827\">{title}</text>\n",
            x = chart_left,
            title = Self::escape_svg_text(&title)
        ));
        svg.push_str(&format!(
            "<text x=\"{x:.1}\" y=\"36\" font-family=\"monospace\" font-size=\"11\" fill=\"#4b5563\">{subtitle}</text>\n",
            x = chart_left,
            subtitle = Self::escape_svg_text(&subtitle)
        ));
        svg.push_str(&format!(
            "<text x=\"{x:.1}\" y=\"50\" font-family=\"monospace\" font-size=\"10\" fill=\"#4b5563\">{provenance}</text>\n",
            x = chart_left,
            provenance = Self::escape_svg_text(&provenance)
        ));
        svg.push_str(&format!(
            "<text x=\"{x:.1}\" y=\"62\" font-family=\"monospace\" font-size=\"10\" fill=\"#6b7280\">{bins_source}</text>\n",
            x = chart_left,
            bins_source = Self::escape_svg_text(bins_source)
        ));
        svg.push_str(&format!(
            "<rect x=\"{x:.1}\" y=\"{y:.1}\" width=\"{w:.1}\" height=\"{h:.1}\" fill=\"none\" stroke=\"#9ca3af\" stroke-width=\"1\"/>\n",
            x = chart_left,
            y = plot_top,
            w = chart_width,
            h = axis_bottom - plot_top
        ));
        for (idx, count) in bins.iter().enumerate() {
            if *count == 0 {
                continue;
            }
            let x0 = chart_left + idx as f64 * bin_width + 0.6;
            let x1 = if idx + 1 == bins.len() {
                chart_right - 0.6
            } else {
                chart_left + (idx + 1) as f64 * bin_width - 0.6
            };
            if x1 <= x0 {
                continue;
            }
            let scaled = if scale == RnaReadScoreDensityScale::Log {
                (*count as f64 + 1.0).ln()
            } else {
                *count as f64
            };
            let bar_height = (scaled / max_scaled) * plot_height;
            let y = axis_bottom - bar_height;
            svg.push_str(&format!(
                "<rect x=\"{x:.2}\" y=\"{y:.2}\" width=\"{w:.2}\" height=\"{h:.2}\" fill=\"#10b981\"/>\n",
                x = x0,
                y = y,
                w = (x1 - x0).max(0.2),
                h = (axis_bottom - y).max(0.2)
            ));
            if (x1 - x0) >= 22.0 {
                let label =
                    Self::escape_svg_text(&Self::format_score_density_count_compact(*count));
                if bar_height >= 12.0 {
                    svg.push_str(&format!(
                        "<text x=\"{x:.2}\" y=\"{y:.2}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"8\" fill=\"#071811\">{label}</text>\n",
                        x = (x0 + x1) * 0.5,
                        y = (y + 8.0).max(plot_top + 8.0),
                    ));
                } else {
                    svg.push_str(&format!(
                        "<text x=\"{x:.2}\" y=\"{y:.2}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"8\" fill=\"#121212\">{label}</text>\n",
                        x = (x0 + x1) * 0.5,
                        y = (y - 2.0).max(plot_top + 8.0),
                    ));
                }
            }
        }
        svg.push_str(&format!(
            "<line x1=\"{x1:.1}\" y1=\"{y:.1}\" x2=\"{x2:.1}\" y2=\"{y:.1}\" stroke=\"#6b7280\" stroke-width=\"1\"/>\n",
            x1 = chart_left,
            x2 = chart_right,
            y = axis_bottom
        ));
        svg.push_str(&format!(
            "<line x1=\"{x:.2}\" y1=\"{y1:.1}\" x2=\"{x:.2}\" y2=\"{y2:.1}\" stroke=\"#dc2626\" stroke-width=\"1\"/>\n",
            x = threshold_x,
            y1 = plot_top,
            y2 = axis_bottom
        ));
        svg.push_str(&format!(
            "<text x=\"{x:.2}\" y=\"{y:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#dc2626\">{label}</text>\n",
            x = threshold_x + 2.0,
            y = plot_top + 10.0,
            label = format!("{threshold:.2}")
        ));
        svg.push_str(&format!(
            "<text x=\"{x:.1}\" y=\"{y:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#6b7280\">0.0</text>\n",
            x = chart_left,
            y = axis_bottom + 14.0
        ));
        svg.push_str(&format!(
            "<text x=\"{x:.1}\" y=\"{y:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10\" fill=\"#6b7280\">1.0</text>\n",
            x = chart_right,
            y = axis_bottom + 14.0
        ));
        svg.push_str(&format!(
            "<text x=\"{x:.1}\" y=\"{y:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#4b5563\">bins={bins} max_count={max_count} total_scored={total_scored}</text>\n",
            x = chart_left,
            y = height - 8.0,
            bins = bins.len(),
            max_count = max_count,
            total_scored = bins.iter().copied().sum::<u64>()
        ));
        svg.push_str("</svg>\n");
        svg
    }

    pub fn export_rna_read_score_density_svg(
        &self,
        report_id: &str,
        path: &str,
        scale: RnaReadScoreDensityScale,
    ) -> Result<RnaReadScoreDensitySvgExport, EngineError> {
        let report = self.get_rna_read_report(report_id)?;
        let path = path.trim();
        if path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "RNA-read score-density SVG export requires non-empty path".to_string(),
            });
        }
        let (bins, derived_from_report_hits_only) =
            Self::derive_score_density_bins_from_report(&report);
        let svg = Self::render_rna_read_score_density_svg_text(&report, &bins, scale);
        std::fs::write(path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write RNA-read score-density SVG to '{path}': {e}"),
        })?;
        Ok(RnaReadScoreDensitySvgExport {
            schema: RNA_READ_SCORE_DENSITY_SVG_EXPORT_SCHEMA.to_string(),
            path: path.to_string(),
            report_id: report.report_id,
            scale,
            bin_count: bins.len(),
            max_bin_count: bins.iter().copied().max().unwrap_or(0),
            total_scored_reads: bins.iter().copied().sum::<u64>(),
            derived_from_report_hits_only,
        })
    }

    pub(super) fn sanitize_tsv_cell(raw: &str) -> String {
        raw.replace(['\t', '\n', '\r'], " ")
    }

    fn ordered_window_overlap_summary(window_len: usize, step_bp: usize) -> String {
        let safe_step = step_bp.max(1);
        let overlap_bp = window_len.saturating_sub(safe_step);
        let windows_per_base = (window_len.saturating_add(safe_step).saturating_sub(1)) / safe_step;
        if safe_step == 1 {
            format!(
                "adjacent windows overlap by {overlap_bp} bp; each interior base participates in {windows_per_base} consecutive ordered windows (dense sliding)"
            )
        } else if safe_step < window_len {
            format!(
                "adjacent windows overlap by {overlap_bp} bp; each interior base participates in up to {windows_per_base} consecutive ordered windows (subsampled sliding)"
            )
        } else if safe_step == window_len {
            "adjacent windows do not overlap; each interior base participates in at most 1 ordered window (edge-touching sampling)".to_string()
        } else {
            format!(
                "adjacent windows do not overlap and can leave up to {} bp unsampled between starts; each interior base participates in at most 1 ordered window (sparse sampling)",
                safe_step - window_len
            )
        }
    }

    fn rna_read_seed_filter_summary(seed_filter: &RnaReadSeedFilterConfig) -> String {
        format!(
            "seed_filter: k={} stride={} min_seed_hit_fraction={:.2} min_weighted_seed_hit_fraction={:.2} min_unique_matched_kmers={} max_median_transcript_gap={:.2} min_chain_consistency_fraction={:.2} min_confirmed_exon_transitions={} min_transition_support_fraction={:.2} poly_t_flip={} poly_t_prefix_min_bp={} | {}",
            seed_filter.kmer_len,
            seed_filter.seed_stride_bp,
            seed_filter.min_seed_hit_fraction,
            seed_filter.min_weighted_seed_hit_fraction,
            seed_filter.min_unique_matched_kmers,
            seed_filter.max_median_transcript_gap,
            seed_filter.min_chain_consistency_fraction,
            seed_filter.min_confirmed_exon_transitions,
            seed_filter.min_transition_support_fraction,
            seed_filter.cdna_poly_t_flip_enabled,
            seed_filter.poly_t_prefix_min_bp,
            Self::ordered_window_overlap_summary(
                seed_filter.kmer_len.max(1),
                seed_filter.seed_stride_bp.max(1),
            ),
        )
    }

    fn format_selected_record_indices_for_metadata(selected_record_indices: &[usize]) -> String {
        if selected_record_indices.is_empty() {
            "none".to_string()
        } else {
            selected_record_indices
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(",")
        }
    }

    fn format_subset_spec_for_metadata(subset_spec: Option<&str>) -> String {
        subset_spec
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(Self::sanitize_tsv_cell)
            .unwrap_or_else(|| "none".to_string())
    }

    fn rna_read_tsv_common_metadata_lines(
        report: &RnaReadInterpretationReport,
        selection: RnaReadHitSelection,
        selected_record_indices: &[usize],
        subset_spec: Option<&str>,
    ) -> Vec<String> {
        let target_gene_ids = if report.target_gene_ids.is_empty() {
            "none".to_string()
        } else {
            report.target_gene_ids.join(",")
        };
        vec![
            format!(
                "# report_id={} seq_id={} selection={} selected_record_indices={} subset_spec={}",
                Self::sanitize_tsv_cell(&report.report_id),
                Self::sanitize_tsv_cell(&report.seq_id),
                selection.as_str(),
                Self::format_selected_record_indices_for_metadata(selected_record_indices),
                Self::format_subset_spec_for_metadata(subset_spec),
            ),
            format!(
                "# profile={} report_mode={} input_format={} scope={} origin_mode={} roi_seed_capture_enabled={} target_gene_ids={}",
                report.profile.as_str(),
                report.report_mode.as_str(),
                report.input_format.as_str(),
                report.scope.as_str(),
                report.origin_mode.as_str(),
                report.roi_seed_capture_enabled,
                Self::sanitize_tsv_cell(&target_gene_ids),
            ),
            format!(
                "# {}",
                Self::rna_read_seed_filter_summary(&report.seed_filter)
            ),
        ]
    }

    fn rna_read_alignment_tsv_metadata_lines(
        report: &RnaReadInterpretationReport,
        inspection: &RnaReadAlignmentInspection,
        selection: RnaReadHitSelection,
        limit: Option<usize>,
        selected_record_indices: &[usize],
        subset_spec: Option<&str>,
    ) -> Vec<String> {
        let limit_text = limit
            .map(|value| value.to_string())
            .unwrap_or_else(|| "all".to_string());
        let mut lines = Self::rna_read_tsv_common_metadata_lines(
            report,
            selection,
            selected_record_indices,
            subset_spec,
        );
        lines[0] = format!(
            "# report_id={} seq_id={} selection={} selected_record_indices={} subset_spec={} limit={} row_count={} aligned_count={}",
            Self::sanitize_tsv_cell(&inspection.report_id),
            Self::sanitize_tsv_cell(&inspection.seq_id),
            selection.as_str(),
            Self::format_selected_record_indices_for_metadata(selected_record_indices),
            Self::format_subset_spec_for_metadata(subset_spec),
            limit_text,
            inspection.row_count,
            inspection.aligned_count,
        );
        lines.push(format!(
            "# align_config: min_identity_fraction={:.2} max_secondary_mappings={}",
            report.align_config.min_identity_fraction, report.align_config.max_secondary_mappings,
        ));
        lines
    }

    pub fn export_rna_read_sample_sheet(
        &self,
        path: &str,
        seq_id_filter: Option<&str>,
        report_ids: &[String],
        append: bool,
    ) -> Result<RnaReadSampleSheetExport, EngineError> {
        let path = path.trim();
        if path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "RNA-read sample sheet export requires non-empty path".to_string(),
            });
        }
        let seq_filter = seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_ascii_lowercase);
        let store = self.read_rna_read_report_store();
        let mut selected = if report_ids.is_empty() {
            store.reports.values().cloned().collect::<Vec<_>>()
        } else {
            let mut out = Vec::<RnaReadInterpretationReport>::with_capacity(report_ids.len());
            for raw_id in report_ids {
                let normalized = Self::normalize_rna_read_report_id(raw_id)?;
                let Some(report) = store.reports.get(normalized.as_str()) else {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("RNA-read report '{}' not found", normalized),
                    });
                };
                out.push(report.clone());
            }
            out
        };
        if let Some(filter) = seq_filter.as_deref() {
            selected.retain(|report| report.seq_id.eq_ignore_ascii_case(filter));
        }
        if selected.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: "No RNA-read reports matched the requested sample-sheet selection"
                    .to_string(),
            });
        }
        selected.sort_by(|left, right| {
            left.seq_id
                .to_ascii_lowercase()
                .cmp(&right.seq_id.to_ascii_lowercase())
                .then(left.generated_at_unix_ms.cmp(&right.generated_at_unix_ms))
                .then(
                    left.report_id
                        .to_ascii_lowercase()
                        .cmp(&right.report_id.to_ascii_lowercase()),
                )
        });

        let existing_nonempty = append
            && std::fs::metadata(path)
                .map(|meta| meta.len() > 0)
                .unwrap_or(false);
        let file = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(!append)
            .append(append)
            .open(path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not open RNA-read sample sheet '{}': {e}", path),
            })?;
        let mut writer = BufWriter::new(file);
        if !existing_nonempty {
            writeln!(
                writer,
                "sample_id\tsample_name\tsample_description\treport_id\tseq_id\tseed_feature_id\tgenerated_at_unix_ms\tinput_path\tprofile\tscope\treport_mode\torigin_mode\ttarget_gene_count\ttarget_gene_ids_json\troi_seed_capture_enabled\tread_count_total\tread_count_seed_passed\tread_count_aligned\tseed_pass_fraction\taligned_fraction\texon_support_frequencies_json\tjunction_support_frequencies_json\torigin_class_counts_json"
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write sample-sheet header to '{}': {e}", path),
            })?;
        }

        for report in &selected {
            let sample_name = Path::new(&report.input_path)
                .file_stem()
                .map(|stem| stem.to_string_lossy().to_string())
                .filter(|value| !value.trim().is_empty())
                .unwrap_or_else(|| report.report_id.clone());
            let sample_description = format!(
                "input={} profile={} scope={}",
                report.input_path,
                report.profile.as_str(),
                report.scope.as_str()
            );
            let seed_pass_fraction = if report.read_count_total == 0 {
                0.0
            } else {
                report.read_count_seed_passed as f64 / report.read_count_total as f64
            };
            let aligned_fraction = if report.read_count_total == 0 {
                0.0
            } else {
                report.read_count_aligned as f64 / report.read_count_total as f64
            };
            let exon_json =
                serde_json::to_string(&report.exon_support_frequencies).map_err(|e| {
                    EngineError {
                        code: ErrorCode::Internal,
                        message: format!(
                            "Could not serialize exon support frequencies for report '{}': {e}",
                            report.report_id
                        ),
                    }
                })?;
            let junction_json = serde_json::to_string(&report.junction_support_frequencies)
                .map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Could not serialize junction support frequencies for report '{}': {e}",
                        report.report_id
                    ),
                })?;
            let target_gene_ids_json =
                serde_json::to_string(&report.target_gene_ids).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Could not serialize target_gene_ids for report '{}': {e}",
                        report.report_id
                    ),
                })?;
            let origin_class_counts_json = serde_json::to_string(&report.origin_class_counts)
                .map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Could not serialize origin_class_counts for report '{}': {e}",
                        report.report_id
                    ),
                })?;
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}",
                Self::sanitize_tsv_cell(&report.report_id),
                Self::sanitize_tsv_cell(&sample_name),
                Self::sanitize_tsv_cell(&sample_description),
                Self::sanitize_tsv_cell(&report.report_id),
                Self::sanitize_tsv_cell(&report.seq_id),
                report.seed_feature_id,
                report.generated_at_unix_ms,
                Self::sanitize_tsv_cell(&report.input_path),
                report.profile.as_str(),
                report.scope.as_str(),
                report.report_mode.as_str(),
                report.origin_mode.as_str(),
                report.target_gene_ids.len(),
                Self::sanitize_tsv_cell(&target_gene_ids_json),
                report.roi_seed_capture_enabled,
                report.read_count_total,
                report.read_count_seed_passed,
                report.read_count_aligned,
                seed_pass_fraction,
                aligned_fraction,
                Self::sanitize_tsv_cell(&exon_json),
                Self::sanitize_tsv_cell(&junction_json),
                Self::sanitize_tsv_cell(&origin_class_counts_json),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not write RNA-read sample sheet row to '{}': {e}",
                    path
                ),
            })?;
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not flush RNA-read sample sheet '{}': {e}", path),
        })?;
        Ok(RnaReadSampleSheetExport {
            schema: RNA_READ_SAMPLE_SHEET_EXPORT_SCHEMA.to_string(),
            path: path.to_string(),
            report_count: selected.len(),
            appended: append,
        })
    }

    pub fn collect_rna_seed_hash_catalog(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
        scope: SplicingScopePreset,
        seed_filter: &RnaReadSeedFilterConfig,
    ) -> Result<Vec<RnaSeedHashCatalogEntry>, EngineError> {
        let templates =
            self.collect_rna_seed_templates(seq_id, seed_feature_id, scope, seed_filter)?;
        Ok(Self::collect_rna_seed_hash_catalog_rows(
            &templates,
            seed_filter.kmer_len,
        ))
    }

    pub fn collect_rna_seed_hash_template_audit(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
        scope: SplicingScopePreset,
        seed_filter: &RnaReadSeedFilterConfig,
    ) -> Result<Vec<RnaSeedHashTemplateAuditEntry>, EngineError> {
        let templates =
            self.collect_rna_seed_templates(seq_id, seed_feature_id, scope, seed_filter)?;
        Ok(Self::collect_rna_seed_hash_template_audit_rows(&templates))
    }

    pub(super) fn collect_rna_seed_templates(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
        scope: SplicingScopePreset,
        seed_filter: &RnaReadSeedFilterConfig,
    ) -> Result<Vec<SplicingTranscriptTemplate>, EngineError> {
        if seed_filter.kmer_len == 0 || seed_filter.kmer_len > 16 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "kmer_len must be within 1..=16".to_string(),
            });
        }
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let splicing = self.build_splicing_expert_view(seq_id, seed_feature_id, scope)?;
        let templates = splicing
            .transcripts
            .iter()
            .map(|lane| Self::make_transcript_template(dna, lane, seed_filter.kmer_len))
            .filter(|template| !template.sequence.is_empty())
            .collect::<Vec<_>>();
        if templates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No transcript templates available for splicing scope '{}' on '{}'",
                    scope.as_str(),
                    seq_id
                ),
            });
        }
        Ok(templates)
    }

    pub fn export_rna_seed_hash_catalog(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
        scope: SplicingScopePreset,
        seed_filter: &RnaReadSeedFilterConfig,
        path: &str,
    ) -> Result<(usize, usize), EngineError> {
        let rows =
            self.collect_rna_seed_hash_catalog(seq_id, seed_feature_id, scope, seed_filter)?;

        let path = path.trim();
        if path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "RNA seed-hash catalog export requires non-empty path".to_string(),
            });
        }
        let file = File::create(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not create RNA seed-hash catalog '{}': {e}", path),
        })?;
        let mut writer = BufWriter::new(file);
        writeln!(
            writer,
            "seed_bits\tseed_bits_hex\tkmer_sequence\ttranscript_feature_id\ttranscript_id\ttranscript_label\tstrand\ttemplate_offset_0based\tgenomic_pos_1based"
        )
        .map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write RNA seed-hash catalog header to '{}': {e}", path),
        })?;

        let mut unique_hashes = HashSet::<u32>::new();
        for row in &rows {
            unique_hashes.insert(row.seed_bits);
            writeln!(
                writer,
                "{}\t{:08X}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.seed_bits,
                row.seed_bits,
                Self::sanitize_tsv_cell(&row.kmer_sequence),
                row.transcript_feature_id,
                Self::sanitize_tsv_cell(&row.transcript_id),
                Self::sanitize_tsv_cell(&row.transcript_label),
                Self::sanitize_tsv_cell(&row.strand),
                row.template_offset_0based,
                row.genomic_pos_1based,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write RNA seed-hash row to '{}': {e}", path),
            })?;
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not flush RNA seed-hash catalog '{}': {e}", path),
        })?;
        Ok((rows.len(), unique_hashes.len()))
    }

    pub(super) fn collect_rna_seed_hash_catalog_rows(
        templates: &[SplicingTranscriptTemplate],
        kmer_len: usize,
    ) -> Vec<RnaSeedHashCatalogEntry> {
        if kmer_len == 0 {
            return Vec::new();
        }
        let mut seen = HashSet::<(u32, usize, bool)>::new();
        let mut rows = Vec::<RnaSeedHashCatalogEntry>::new();
        for template in templates {
            let strand_minus = template.strand.trim() == "-";
            if template.sequence.len() < kmer_len {
                continue;
            }
            for start in 0..=template.sequence.len() - kmer_len {
                let window = &template.sequence[start..start + kmer_len];
                let Some(seed_bits) = Self::encode_kmer_bits(window) else {
                    continue;
                };
                let genomic_pos_1based = template
                    .genomic_positions_1based
                    .get(start)
                    .copied()
                    .unwrap_or(1);
                if !seen.insert((seed_bits, genomic_pos_1based, strand_minus)) {
                    continue;
                }
                rows.push(RnaSeedHashCatalogEntry {
                    seed_bits,
                    kmer_sequence: String::from_utf8_lossy(window).to_string(),
                    transcript_feature_id: template.transcript_feature_id,
                    transcript_id: template.transcript_id.clone(),
                    transcript_label: template.transcript_label.clone(),
                    strand: template.strand.clone(),
                    template_offset_0based: start,
                    genomic_pos_1based,
                });
            }
        }
        rows.sort_by(|left, right| {
            left.genomic_pos_1based
                .cmp(&right.genomic_pos_1based)
                .then_with(|| left.strand.cmp(&right.strand))
                .then_with(|| left.transcript_feature_id.cmp(&right.transcript_feature_id))
                .then_with(|| {
                    left.template_offset_0based
                        .cmp(&right.template_offset_0based)
                })
                .then_with(|| left.seed_bits.cmp(&right.seed_bits))
        });
        rows
    }

    pub(super) fn collect_rna_seed_hash_template_audit_rows(
        templates: &[SplicingTranscriptTemplate],
    ) -> Vec<RnaSeedHashTemplateAuditEntry> {
        let mut rows = templates
            .iter()
            .map(|template| RnaSeedHashTemplateAuditEntry {
                transcript_feature_id: template.transcript_feature_id,
                transcript_id: template.transcript_id.clone(),
                transcript_label: template.transcript_label.clone(),
                strand: template.strand.clone(),
                template_sequence: String::from_utf8_lossy(&template.sequence).to_string(),
                template_length_bp: template.sequence.len(),
                template_first_genomic_pos_1based: template
                    .genomic_positions_1based
                    .first()
                    .copied()
                    .unwrap_or(0),
                template_last_genomic_pos_1based: template
                    .genomic_positions_1based
                    .last()
                    .copied()
                    .unwrap_or(0),
                reverse_complemented_from_genome: template.strand.trim() == "-",
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.strand
                .cmp(&right.strand)
                .then_with(|| left.transcript_id.cmp(&right.transcript_id))
                .then_with(|| left.transcript_feature_id.cmp(&right.transcript_feature_id))
        });
        rows
    }

    pub(super) fn visit_fasta_records_with_offsets(
        path: &str,
        on_record: &mut dyn FnMut(
            ParsedFastaReadRecord,
            FastaVisitProgress,
        ) -> Result<(), EngineError>,
    ) -> Result<FastaVisitProgress, EngineError> {
        let file = File::open(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not open FASTA input '{}': {e}", path),
        })?;
        let input_bytes_total = file.metadata().map(|meta| meta.len()).unwrap_or(0);
        let source_bytes_read = Arc::new(AtomicU64::new(0));
        let counted_file = CountingReader::new(file, Arc::clone(&source_bytes_read));
        let lower = path.to_ascii_lowercase();
        let mut reader: Box<dyn BufRead> = if lower.ends_with(".gz") {
            let decoder = MultiGzDecoder::new(BufReader::new(counted_file));
            Box::new(BufReader::new(decoder))
        } else {
            Box::new(BufReader::new(counted_file))
        };
        let mut processed = 0usize;
        let mut current_header: Option<String> = None;
        let mut current_header_offset: usize = 0usize;
        let mut current_sequence = String::new();
        let mut offset = 0usize;
        let mut line = String::new();
        let io_read_ms = Cell::new(0.0f64);
        let record_parse_ms = Cell::new(0.0f64);
        let mut flush_record = |header: Option<String>,
                                header_offset: usize,
                                sequence: &mut String|
         -> Result<(), EngineError> {
            let Some(header_text) = header else {
                return Ok(());
            };
            let parse_started = Instant::now();
            let normalized = sequence
                .chars()
                .filter(|ch| !ch.is_whitespace())
                .map(|ch| {
                    let upper = ch.to_ascii_uppercase();
                    if upper == 'U' { 'T' } else { upper }
                })
                .collect::<String>();
            sequence.clear();
            if normalized.is_empty() {
                return Ok(());
            }
            let record_index = processed;
            let header_id = header_text
                .trim()
                .split_ascii_whitespace()
                .next()
                .filter(|v| !v.is_empty())
                .map(str::to_string)
                .unwrap_or_else(|| format!("record_{}", record_index + 1));
            processed = processed.saturating_add(1);
            let input_bytes_processed = if input_bytes_total > 0 {
                source_bytes_read
                    .load(AtomicOrdering::Relaxed)
                    .min(input_bytes_total)
            } else {
                source_bytes_read.load(AtomicOrdering::Relaxed)
            };
            record_parse_ms
                .set(record_parse_ms.get() + parse_started.elapsed().as_secs_f64() * 1000.0);
            on_record(
                ParsedFastaReadRecord {
                    record_index,
                    source_byte_offset: header_offset,
                    header_id,
                    sequence: normalized.as_bytes().to_vec(),
                },
                FastaVisitProgress {
                    records_processed: processed,
                    input_bytes_processed,
                    input_bytes_total,
                    io_read_ms: io_read_ms.get(),
                    record_parse_ms: record_parse_ms.get(),
                },
            )?;
            Ok(())
        };

        loop {
            line.clear();
            let read_started = Instant::now();
            let bytes = reader.read_line(&mut line).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read FASTA input '{}': {e}", path),
            })?;
            io_read_ms.set(io_read_ms.get() + read_started.elapsed().as_secs_f64() * 1000.0);
            if bytes == 0 {
                break;
            }
            let line_start_offset = offset;
            offset = offset.saturating_add(bytes);
            if line.starts_with('>') {
                flush_record(
                    current_header.take(),
                    current_header_offset,
                    &mut current_sequence,
                )?;
                current_header = Some(line[1..].trim().to_string());
                current_header_offset = line_start_offset;
                continue;
            }
            if current_header.is_none() {
                continue;
            }
            current_sequence.push_str(line.trim());
        }
        flush_record(
            current_header.take(),
            current_header_offset,
            &mut current_sequence,
        )?;
        if processed == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("No FASTA records found in '{}'", path),
            });
        }
        let input_bytes_processed = if input_bytes_total > 0 {
            input_bytes_total
        } else {
            source_bytes_read.load(AtomicOrdering::Relaxed)
        };
        Ok(FastaVisitProgress {
            records_processed: processed,
            input_bytes_processed,
            input_bytes_total,
            io_read_ms: io_read_ms.get(),
            record_parse_ms: record_parse_ms.get(),
        })
    }

    #[cfg(test)]
    pub(super) fn parse_fasta_records_with_offsets_with_progress(
        path: &str,
        on_record_progress: &mut dyn FnMut(usize) -> bool,
    ) -> Result<Vec<ParsedFastaReadRecord>, EngineError> {
        let mut out = Vec::<ParsedFastaReadRecord>::new();
        let mut processed = 0usize;
        Self::visit_fasta_records_with_offsets(path, &mut |record, progress| {
            out.push(record);
            processed = progress.records_processed;
            if (processed % RNA_READ_PROGRESS_UPDATE_EVERY_READS == 0 || processed <= 3)
                && !on_record_progress(processed)
            {
                return Err(EngineError {
                    code: ErrorCode::Internal,
                    message: "RNA-read FASTA parsing cancelled during progress reporting"
                        .to_string(),
                });
            }
            Ok(())
        })?;
        if !on_record_progress(processed) {
            return Err(EngineError {
                code: ErrorCode::Internal,
                message: "RNA-read FASTA parsing cancelled at completion".to_string(),
            });
        }
        Ok(out)
    }

    #[cfg(test)]
    pub(super) fn parse_fasta_records_with_offsets(
        path: &str,
    ) -> Result<Vec<ParsedFastaReadRecord>, EngineError> {
        let mut noop = |_processed: usize| true;
        Self::parse_fasta_records_with_offsets_with_progress(path, &mut noop)
    }

    pub(super) fn normalize_nucleotide_base(base: u8) -> u8 {
        match base.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        }
    }

    pub(super) fn complement_nucleotide_base(base: u8) -> u8 {
        match Self::normalize_nucleotide_base(base) {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'N' => b'N',
            _ => b'N',
        }
    }

    pub(super) fn reverse_complement_sequence(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| Self::complement_nucleotide_base(*base))
            .collect::<Vec<_>>()
    }

    // Nanopore cDNA tails are often T-rich but can contain a few interruptions
    // (e.g. one or two mismatching bases) near the read head.
    pub(super) fn has_poly_t_prefix(sequence: &[u8], min_prefix_bp: usize) -> bool {
        const POLY_T_HEAD_MAX_SCAN_BP: usize = 96;
        const POLY_T_HEAD_WINDOW_PADDING_BP: usize = 8;
        const POLY_T_HEAD_MAX_WINDOW_START_BP: usize = 24;
        const POLY_T_HEAD_MIN_WINDOW_FRACTION: f64 = 0.60;
        const POLY_T_HEAD_MIN_CONTIGUOUS_RUN_BP: usize = 8;

        let required = min_prefix_bp.max(1);
        if sequence.len() < required {
            return false;
        }

        let scan_len = sequence
            .len()
            .min((required.saturating_add(POLY_T_HEAD_WINDOW_PADDING_BP)).max(32))
            .min(POLY_T_HEAD_MAX_SCAN_BP);
        if scan_len < required {
            return false;
        }

        let mut head = Vec::<u8>::with_capacity(scan_len);
        for base in sequence.iter().take(scan_len) {
            head.push(Self::normalize_nucleotide_base(*base));
        }

        let mut leading_t = 0usize;
        for base in &head {
            if *base == b'T' {
                leading_t = leading_t.saturating_add(1);
            } else {
                break;
            }
        }
        if leading_t >= required {
            return true;
        }

        let mut longest_t_run = 0usize;
        let mut current_t_run = 0usize;
        for base in &head {
            if *base == b'T' {
                current_t_run = current_t_run.saturating_add(1);
                longest_t_run = longest_t_run.max(current_t_run);
            } else {
                current_t_run = 0;
            }
        }
        if longest_t_run < POLY_T_HEAD_MIN_CONTIGUOUS_RUN_BP {
            return false;
        }

        let window_len = required
            .saturating_add(POLY_T_HEAD_WINDOW_PADDING_BP)
            .max(required)
            .min(scan_len);
        if window_len == 0 {
            return false;
        }
        let max_start = scan_len
            .saturating_sub(window_len)
            .min(POLY_T_HEAD_MAX_WINDOW_START_BP);
        for start in 0..=max_start {
            let end = start.saturating_add(window_len);
            let t_count = head[start..end]
                .iter()
                .filter(|base| **base == b'T')
                .count();
            if t_count < required {
                continue;
            }
            let fraction = t_count as f64 / window_len as f64;
            if fraction >= POLY_T_HEAD_MIN_WINDOW_FRACTION {
                return true;
            }
        }
        false
    }

    pub(super) fn normalize_rna_read_sequence_for_scoring(
        sequence: &[u8],
        seed_filter: &RnaReadSeedFilterConfig,
    ) -> (Vec<u8>, bool) {
        if seed_filter.cdna_poly_t_flip_enabled
            && Self::has_poly_t_prefix(sequence, seed_filter.poly_t_prefix_min_bp)
        {
            (Self::reverse_complement_sequence(sequence), true)
        } else {
            (
                sequence
                    .iter()
                    .map(|base| Self::normalize_nucleotide_base(*base))
                    .collect::<Vec<_>>(),
                false,
            )
        }
    }

    pub(super) fn base_to_2bit(base: u8) -> Option<u32> {
        match Self::normalize_nucleotide_base(base) {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        }
    }

    pub(super) fn encode_kmer_bits(window: &[u8]) -> Option<u32> {
        let mut bits = 0u32;
        for base in window {
            let value = Self::base_to_2bit(*base)?;
            bits = (bits << 2) | value;
        }
        Some(bits)
    }

    pub(super) fn sampled_read_windows(
        read_len: usize,
        _seed_filter: &RnaReadSeedFilterConfig,
    ) -> Vec<(usize, usize)> {
        if read_len == 0 {
            return vec![];
        }
        vec![(0, read_len)]
    }

    pub(super) fn build_rna_read_seed_histogram_bins(
        seq_len_bp: usize,
    ) -> Vec<RnaReadSeedHistogramBin> {
        if seq_len_bp == 0 {
            return vec![];
        }
        let bin_count = seq_len_bp.min(RNA_READ_PROGRESS_MAX_HISTOGRAM_BINS).max(1);
        let bin_size_bp = seq_len_bp.div_ceil(bin_count);
        (0..bin_count)
            .map(|idx| {
                let start_1based = idx.saturating_mul(bin_size_bp).saturating_add(1);
                let end_1based = ((idx + 1).saturating_mul(bin_size_bp)).min(seq_len_bp);
                RnaReadSeedHistogramBin {
                    start_1based,
                    end_1based,
                    confirmed_plus: 0,
                    confirmed_minus: 0,
                }
            })
            .collect()
    }

    pub(super) fn build_rna_read_seed_histogram_index(
        templates: &[SplicingTranscriptTemplate],
        seq_len_bp: usize,
        bins: &[RnaReadSeedHistogramBin],
    ) -> HashMap<u32, Vec<SeedHistogramWeight>> {
        if seq_len_bp == 0 || bins.is_empty() {
            return HashMap::new();
        }
        let bin_size_bp = bins[0]
            .end_1based
            .saturating_sub(bins[0].start_1based)
            .saturating_add(1)
            .max(1);
        let mut expanded_positions: HashMap<u32, HashSet<(usize, bool)>> = HashMap::new();
        for template in templates {
            let strand_minus = template.strand.trim() == "-";
            for (bits, starts) in &template.kmer_positions {
                let entry = expanded_positions.entry(*bits).or_default();
                for start in starts {
                    let Some(pos_1based) = template.genomic_positions_1based.get(*start).copied()
                    else {
                        continue;
                    };
                    let pos_0based = pos_1based
                        .saturating_sub(1)
                        .min(seq_len_bp.saturating_sub(1));
                    entry.insert((pos_0based, strand_minus));
                }
            }
        }
        expanded_positions
            .into_iter()
            .map(|(bits, by_pos)| {
                let mut by_bin = HashSet::<(usize, bool)>::new();
                for (pos_0based, strand_minus) in by_pos {
                    let bin_index = (pos_0based / bin_size_bp).min(bins.len().saturating_sub(1));
                    by_bin.insert((bin_index, strand_minus));
                }
                let mut weights = by_bin
                    .into_iter()
                    .map(|(bin_index, strand_minus)| SeedHistogramWeight {
                        bin_index,
                        strand_minus,
                    })
                    .collect::<Vec<_>>();
                weights.sort_by(|left, right| {
                    left.bin_index
                        .cmp(&right.bin_index)
                        .then(left.strand_minus.cmp(&right.strand_minus))
                });
                (bits, weights)
            })
            .collect()
    }

    pub(super) fn build_seed_template_position_index(
        templates: &[SplicingTranscriptTemplate],
    ) -> HashMap<u32, Vec<SeedTemplatePosition>> {
        let mut index = HashMap::<u32, HashSet<(usize, usize)>>::new();
        for (template_idx, template) in templates.iter().enumerate() {
            for (bits, starts) in &template.kmer_positions {
                let entry = index.entry(*bits).or_default();
                for start in starts {
                    entry.insert((template_idx, *start));
                }
            }
        }
        index
            .into_iter()
            .map(|(bits, rows)| {
                let mut positions = rows
                    .into_iter()
                    .map(|(template_idx, template_pos)| SeedTemplatePosition {
                        template_idx,
                        template_pos,
                    })
                    .collect::<Vec<_>>();
                positions.sort_by(|left, right| {
                    left.template_idx
                        .cmp(&right.template_idx)
                        .then(left.template_pos.cmp(&right.template_pos))
                });
                (bits, positions)
            })
            .collect()
    }

    pub(super) fn count_seed_hits_in_window_with_histogram(
        sequence: &[u8],
        read_start_offset: usize,
        kmer_len: usize,
        seed_stride_bp: usize,
        seed_index: &HashSet<u32>,
        histogram_index: &HashMap<u32, Vec<SeedHistogramWeight>>,
        bins: &mut [RnaReadSeedHistogramBin],
        _seed_occurrence_counts: &HashMap<u32, usize>,
    ) -> (usize, usize, HashSet<u32>, Vec<SeedMatchObservation>) {
        if kmer_len == 0 || sequence.len() < kmer_len {
            return (0, 0, HashSet::new(), vec![]);
        }
        let mut tested = 0usize;
        let mut matched = 0usize;
        let mut matched_seed_bits = HashSet::<u32>::new();
        let mut matched_seed_observations = Vec::<SeedMatchObservation>::new();
        let mut plus_bins_touched = HashSet::<usize>::new();
        let mut minus_bins_touched = HashSet::<usize>::new();
        let stride = seed_stride_bp.max(1);
        for start in (0..=sequence.len() - kmer_len).step_by(stride) {
            let window = &sequence[start..start + kmer_len];
            let Some(bits) = Self::encode_kmer_bits(window) else {
                continue;
            };
            tested += 1;
            if seed_index.contains(&bits) {
                matched += 1;
                matched_seed_bits.insert(bits);
                matched_seed_observations.push(SeedMatchObservation {
                    read_start: read_start_offset.saturating_add(start),
                    bits,
                });
            }
        }
        for bits in &matched_seed_bits {
            if let Some(weights) = histogram_index.get(&bits) {
                for weight in weights {
                    if weight.strand_minus {
                        minus_bins_touched.insert(weight.bin_index);
                    } else {
                        plus_bins_touched.insert(weight.bin_index);
                    }
                }
            }
        }
        for bin_index in plus_bins_touched {
            if let Some(bin) = bins.get_mut(bin_index) {
                bin.confirmed_plus = bin.confirmed_plus.saturating_add(1);
            }
        }
        for bin_index in minus_bins_touched {
            if let Some(bin) = bins.get_mut(bin_index) {
                bin.confirmed_minus = bin.confirmed_minus.saturating_add(1);
            }
        }
        (
            tested,
            matched,
            matched_seed_bits,
            matched_seed_observations,
        )
    }

    pub(super) fn seed_hit_metrics(
        tested_kmers: usize,
        matched_kmers: usize,
        min_seed_hit_fraction: f64,
    ) -> (f64, bool, bool) {
        let seed_hit_fraction = if tested_kmers == 0 {
            0.0
        } else {
            matched_kmers as f64 / tested_kmers as f64
        };
        let perfect_seed_match = tested_kmers > 0 && tested_kmers == matched_kmers;
        let passed_seed_filter = seed_hit_fraction + f64::EPSILON >= min_seed_hit_fraction;
        (seed_hit_fraction, perfect_seed_match, passed_seed_filter)
    }

    pub(super) fn weighted_seed_support_from_occurrences(
        matched_seed_bits: &HashSet<u32>,
        seed_occurrence_counts: &HashMap<u32, usize>,
    ) -> f64 {
        matched_seed_bits
            .iter()
            .map(|bits| {
                let occurrence_count = seed_occurrence_counts
                    .get(bits)
                    .copied()
                    .unwrap_or(1)
                    .max(1);
                1.0 / occurrence_count as f64
            })
            .sum::<f64>()
    }

    pub(super) fn score_density_bin_index(seed_hit_fraction: f64) -> usize {
        if RNA_READ_SCORE_DENSITY_BIN_COUNT <= 1 {
            return 0;
        }
        let clamped = seed_hit_fraction.clamp(0.0, 1.0);
        let scaled = (clamped * RNA_READ_SCORE_DENSITY_BIN_COUNT as f64).floor() as usize;
        scaled.min(RNA_READ_SCORE_DENSITY_BIN_COUNT - 1)
    }

    pub(super) fn update_read_length_counts(length_counts: &mut Vec<u64>, read_length_bp: usize) {
        if length_counts.len() <= read_length_bp {
            length_counts.resize(read_length_bp.saturating_add(1), 0);
        }
        if let Some(bucket) = length_counts.get_mut(read_length_bp) {
            *bucket = bucket.saturating_add(1);
        }
    }

    pub(super) fn quantile_read_length_from_counts(
        length_counts: &[u64],
        total_reads: usize,
        quantile: f64,
    ) -> usize {
        if total_reads == 0 || length_counts.is_empty() {
            return 0;
        }
        let q = quantile.clamp(0.0, 1.0);
        let rank = ((total_reads as f64 * q).ceil() as usize).clamp(1, total_reads);
        let mut cumulative = 0usize;
        for (len_bp, count) in length_counts.iter().enumerate() {
            cumulative = cumulative.saturating_add((*count).min(usize::MAX as u64) as usize);
            if cumulative >= rank {
                return len_bp;
            }
        }
        length_counts.len().saturating_sub(1)
    }

    pub(super) fn summarize_read_lengths(
        length_counts: &[u64],
        total_reads: usize,
        total_bases: u64,
    ) -> (f64, usize, usize) {
        if total_reads == 0 {
            return (0.0, 0, 0);
        }
        let mean = total_bases as f64 / total_reads as f64;
        let median = Self::quantile_read_length_from_counts(length_counts, total_reads, 0.50);
        let p95 = Self::quantile_read_length_from_counts(length_counts, total_reads, 0.95);
        (mean, median, p95)
    }

    pub(super) fn seed_filter_passes(
        seed_hit_fraction: f64,
        weighted_seed_hit_fraction: f64,
        tested_kmers: usize,
        unique_matched_kmers: usize,
        chain_consistency_fraction: f64,
        median_transcript_gap: f64,
        transcript_gap_count: usize,
        confirmed_transitions: usize,
        total_transitions: usize,
        seed_filter: &RnaReadSeedFilterConfig,
    ) -> bool {
        let raw_pass = seed_hit_fraction + f64::EPSILON >= seed_filter.min_seed_hit_fraction;
        let weighted_pass =
            weighted_seed_hit_fraction + f64::EPSILON >= seed_filter.min_weighted_seed_hit_fraction;
        let unique_required = seed_filter
            .min_unique_matched_kmers
            .min(tested_kmers.max(1));
        let unique_pass = unique_matched_kmers >= unique_required;
        let chain_pass =
            chain_consistency_fraction + f64::EPSILON >= seed_filter.min_chain_consistency_fraction;
        let gap_pass = if transcript_gap_count == 0 {
            tested_kmers <= 1
        } else {
            median_transcript_gap <= seed_filter.max_median_transcript_gap + f64::EPSILON
        };
        let transition_pass =
            Self::transition_gate_passes(confirmed_transitions, total_transitions, seed_filter);
        raw_pass && weighted_pass && unique_pass && chain_pass && gap_pass && transition_pass
    }

    pub(super) fn compute_msa_eligibility(
        hit: &RnaReadInterpretationHit,
        seed_filter: &RnaReadSeedFilterConfig,
    ) -> (bool, String) {
        if !hit.passed_seed_filter {
            return (false, "seed gate failed".to_string());
        }
        let Some(mapping) = &hit.best_mapping else {
            return (false, "no transcript alignment".to_string());
        };
        if mapping.query_coverage_fraction + f64::EPSILON < 0.60 {
            return (
                false,
                format!(
                    "low query coverage {:.3} < 0.600",
                    mapping.query_coverage_fraction
                ),
            );
        }
        if hit.seed_transcript_gap_count > 0
            && hit.seed_median_transcript_gap > seed_filter.max_median_transcript_gap + f64::EPSILON
        {
            return (
                false,
                format!(
                    "seed gap median {:.2} > {:.2}",
                    hit.seed_median_transcript_gap, seed_filter.max_median_transcript_gap
                ),
            );
        }
        if !Self::transition_gate_passes(
            hit.exon_transitions_confirmed,
            hit.exon_transitions_total,
            seed_filter,
        ) {
            return (
                false,
                format!(
                    "transition gate failed ({}/{})",
                    hit.exon_transitions_confirmed, hit.exon_transitions_total
                ),
            );
        }
        (
            true,
            format!(
                "seed+alignment support (identity {:.3}, coverage {:.3})",
                mapping.identity_fraction, mapping.query_coverage_fraction
            ),
        )
    }

    pub(super) fn classify_rna_read_origin(
        seed_hit_fraction: f64,
        weighted_seed_hit_fraction: f64,
        passed_seed_filter: bool,
        spacing_metrics: &SeedChainSpacingMetrics,
        path_inference: &ReadExonPathInference,
        target_feature_strand: &str,
        seed_filter: &RnaReadSeedFilterConfig,
    ) -> RnaReadOriginClassification {
        let selected_strand = path_inference.strand.trim();
        let target_strand = target_feature_strand.trim();
        let has_target_strand = !target_strand.is_empty();
        let has_selected_strand = !selected_strand.is_empty();
        let strand_matches_target =
            has_target_strand && has_selected_strand && selected_strand == target_strand;
        let transition_fraction = if path_inference.total_transitions == 0 {
            0.0
        } else {
            path_inference.confirmed_transitions as f64 / path_inference.total_transitions as f64
        };
        let has_local_block = seed_hit_fraction + f64::EPSILON >= seed_filter.min_seed_hit_fraction
            || weighted_seed_hit_fraction + f64::EPSILON
                >= seed_filter.min_weighted_seed_hit_fraction;
        let coherent_chain = spacing_metrics.support_fraction + f64::EPSILON
            >= seed_filter.min_chain_consistency_fraction;
        let strand_confidence = if !has_selected_strand {
            0.0
        } else if path_inference.strand_diagnostics.ambiguous_near_tie {
            0.35
        } else if path_inference.strand_diagnostics.competing_opposite_strand {
            0.70
        } else {
            1.0
        };
        let origin_confidence = (weighted_seed_hit_fraction.clamp(0.0, 1.0) * 0.45
            + spacing_metrics.support_fraction.clamp(0.0, 1.0) * 0.35
            + transition_fraction.clamp(0.0, 1.0) * 0.20)
            .clamp(0.0, 1.0);
        let (origin_class, reason) = if passed_seed_filter && coherent_chain {
            (
                RnaReadOriginClass::TargetCoherent,
                "passed strict seed gate with coherent transcript-chain support".to_string(),
            )
        } else if has_local_block && coherent_chain {
            (
                RnaReadOriginClass::TargetPartialLocalBlock,
                "strong local seed support with coherent chain, but strict gate not fully met"
                    .to_string(),
            )
        } else if has_local_block && has_target_strand && has_selected_strand {
            if strand_matches_target {
                (
                    RnaReadOriginClass::RoiSameStrandLocalBlock,
                    "local ROI support on target-feature strand".to_string(),
                )
            } else {
                (
                    RnaReadOriginClass::RoiReverseStrandLocalBlock,
                    "local ROI support on reverse strand relative to target feature".to_string(),
                )
            }
        } else if path_inference.strand_diagnostics.ambiguous_near_tie
            || path_inference.strand_diagnostics.competing_opposite_strand
        {
            (
                RnaReadOriginClass::TpFamilyAmbiguous,
                "cross-strand competition/ambiguity suggests family-like or shared-domain seed support"
                    .to_string(),
            )
        } else {
            (
                RnaReadOriginClass::BackgroundLikely,
                "insufficient coherent local support under current seed/coherence model"
                    .to_string(),
            )
        };
        RnaReadOriginClassification {
            origin_class,
            reason,
            origin_confidence,
            strand_confidence,
        }
    }

    pub(super) fn build_rna_read_origin_candidates(
        path_inference: &ReadExonPathInference,
        spacing_metrics: &SeedChainSpacingMetrics,
        transcript_models_by_id: &HashMap<String, TranscriptExonPathModel>,
    ) -> Vec<RnaReadOriginCandidateContribution> {
        let mut rows = Vec::<RnaReadOriginCandidateContribution>::new();
        let mut seen = HashSet::<String>::new();
        let mut push_row = |candidate_role: &str,
                            transcript_id: &str,
                            strand: &str,
                            transition_hits: usize,
                            exon_hits: usize| {
            let tid = transcript_id.trim();
            if tid.is_empty() || !seen.insert(tid.to_string()) {
                return;
            }
            rows.push(RnaReadOriginCandidateContribution {
                candidate_role: candidate_role.to_string(),
                transcript_id: tid.to_string(),
                strand: strand.to_string(),
                transition_hits,
                exon_hits,
            });
        };
        push_row(
            "selected_path",
            &path_inference.transcript_id,
            &path_inference.strand,
            path_inference.confirmed_transitions,
            path_inference.strand_diagnostics.selected_exon_hits,
        );
        push_row(
            "plus_best",
            &path_inference.strand_diagnostics.plus_best_transcript_id,
            "+",
            path_inference.strand_diagnostics.plus_best_transition_hits,
            path_inference.strand_diagnostics.plus_best_exon_hits,
        );
        push_row(
            "minus_best",
            &path_inference.strand_diagnostics.minus_best_transcript_id,
            "-",
            path_inference.strand_diagnostics.minus_best_transition_hits,
            path_inference.strand_diagnostics.minus_best_exon_hits,
        );
        if !spacing_metrics.transcript_id.trim().is_empty() {
            let seed_chain_strand = transcript_models_by_id
                .get(spacing_metrics.transcript_id.as_str())
                .map(|model| model.strand.as_str())
                .unwrap_or("");
            push_row(
                "seed_chain_best",
                spacing_metrics.transcript_id.as_str(),
                seed_chain_strand,
                0,
                spacing_metrics.support_kmers,
            );
        }
        rows
    }

    pub(super) fn build_read_seed_support_sets(
        matched_seed_bits: &HashSet<u32>,
        seed_to_exons: &HashMap<u32, Vec<usize>>,
        seed_to_transitions: &HashMap<u32, Vec<(usize, usize)>>,
    ) -> (HashSet<usize>, HashSet<(usize, usize)>) {
        let mut read_supported_exons = HashSet::<usize>::new();
        let mut read_supported_transitions = HashSet::<(usize, usize)>::new();
        for bits in matched_seed_bits {
            if let Some(exons) = seed_to_exons.get(bits) {
                read_supported_exons.extend(exons.iter().copied());
            }
            if let Some(transitions) = seed_to_transitions.get(bits) {
                read_supported_transitions.extend(transitions.iter().copied());
            }
        }
        (read_supported_exons, read_supported_transitions)
    }

    pub(super) fn compute_seed_chain_spacing_metrics(
        matched_seed_observations: &[SeedMatchObservation],
        seed_template_positions: &HashMap<u32, Vec<SeedTemplatePosition>>,
        templates: &[SplicingTranscriptTemplate],
    ) -> SeedChainSpacingMetrics {
        let mut key_counts = HashMap::<(usize, isize), usize>::new();
        for observation in matched_seed_observations {
            let Some(positions) = seed_template_positions.get(&observation.bits) else {
                continue;
            };
            for pos in positions
                .iter()
                .take(RNA_READ_SEED_CHAIN_MAX_CANDIDATES_PER_BIT)
            {
                let delta = pos.template_pos as isize - observation.read_start as isize;
                *key_counts.entry((pos.template_idx, delta)).or_insert(0) += 1;
            }
        }
        let Some(((best_template_idx, best_delta), best_support)) =
            key_counts.into_iter().max_by(|left, right| {
                left.1
                    .cmp(&right.1)
                    .then_with(|| right.0.0.cmp(&left.0.0))
                    .then_with(|| right.0.1.cmp(&left.0.1))
            })
        else {
            return SeedChainSpacingMetrics {
                median_transcript_gap: -1.0,
                ..SeedChainSpacingMetrics::default()
            };
        };
        let support_fraction = if matched_seed_observations.is_empty() {
            0.0
        } else {
            best_support as f64 / matched_seed_observations.len() as f64
        };
        let mut transcript_positions = Vec::<usize>::new();
        for observation in matched_seed_observations {
            let Some(positions) = seed_template_positions.get(&observation.bits) else {
                continue;
            };
            if let Some(pos) = positions
                .iter()
                .take(RNA_READ_SEED_CHAIN_MAX_CANDIDATES_PER_BIT)
                .find(|pos| {
                    pos.template_idx == best_template_idx
                        && (pos.template_pos as isize - observation.read_start as isize)
                            == best_delta
                })
            {
                transcript_positions.push(pos.template_pos);
            }
        }
        transcript_positions.sort_unstable();
        transcript_positions.dedup();
        let mut gaps = transcript_positions
            .windows(2)
            .map(|pair| pair[1].saturating_sub(pair[0]))
            .collect::<Vec<_>>();
        let transcript_gap_count = gaps.len();
        let median_transcript_gap = if gaps.is_empty() {
            -1.0
        } else {
            gaps.sort_unstable();
            let mid = gaps.len() / 2;
            if gaps.len() % 2 == 0 {
                (gaps[mid - 1] as f64 + gaps[mid] as f64) * 0.5
            } else {
                gaps[mid] as f64
            }
        };
        SeedChainSpacingMetrics {
            transcript_id: templates
                .get(best_template_idx)
                .map(|template| template.transcript_id.clone())
                .unwrap_or_default(),
            support_kmers: best_support,
            support_fraction,
            median_transcript_gap,
            transcript_gap_count,
        }
    }

    pub(super) fn rna_read_retention_rank(hit: &RnaReadInterpretationHit) -> RnaReadRetentionRank {
        let weighted_support_milli = (hit.weighted_matched_kmers.max(0.0) * 1000.0).round() as u64;
        let weighted_seed_hit_ppm =
            (hit.weighted_seed_hit_fraction.clamp(0.0, 1.0) * 1_000_000.0).round() as u32;
        let seed_hit_ppm = (hit.seed_hit_fraction.clamp(0.0, 1.0) * 1_000_000.0).round() as u32;
        let (has_alignment, alignment_identity_ppm, alignment_query_coverage_ppm, alignment_score) =
            if let Some(mapping) = hit.best_mapping.as_ref() {
                (
                    true,
                    (mapping.identity_fraction.clamp(0.0, 1.0) * 1_000_000.0).round() as u32,
                    (mapping.query_coverage_fraction.clamp(0.0, 1.0) * 1_000_000.0).round() as u32,
                    mapping.score as i64,
                )
            } else {
                (false, 0, 0, i64::MIN)
            };
        RnaReadRetentionRank {
            passed_seed_filter: hit.passed_seed_filter,
            has_alignment,
            alignment_identity_ppm,
            alignment_query_coverage_ppm,
            alignment_score,
            weighted_support_milli,
            weighted_seed_hit_ppm,
            seed_hit_ppm,
            matched_kmers: hit.matched_kmers,
            tested_kmers: hit.tested_kmers,
            read_length_bp: hit.read_length_bp,
            record_index: hit.record_index,
        }
    }

    pub(super) fn sort_rna_read_hits_by_retention_rank(hits: &mut [RnaReadInterpretationHit]) {
        hits.sort_by(|left, right| {
            let left_rank = Self::rna_read_retention_rank(left);
            let right_rank = Self::rna_read_retention_rank(right);
            right_rank.cmp(&left_rank)
        });
    }

    pub(super) fn should_emit_rna_read_progress(
        reads_processed: usize,
        elapsed_since_last_emit: Duration,
    ) -> bool {
        reads_processed <= 3
            || reads_processed % RNA_READ_PROGRESS_UPDATE_EVERY_READS == 0
            || elapsed_since_last_emit >= RNA_READ_PROGRESS_UPDATE_MAX_INTERVAL
    }

    pub(super) fn retain_top_rna_read_hit(
        retained_hits: &mut BinaryHeap<RetainedRnaReadHit>,
        hit: RnaReadInterpretationHit,
    ) {
        let rank = Self::rna_read_retention_rank(&hit);
        if retained_hits.len() < RNA_READ_RETAINED_HITS_MAX {
            retained_hits.push(RetainedRnaReadHit { rank, hit });
            return;
        }
        let should_replace = retained_hits
            .peek()
            .map(|worst| rank > worst.rank)
            .unwrap_or(true);
        if should_replace {
            let _ = retained_hits.pop();
            retained_hits.push(RetainedRnaReadHit { rank, hit });
        }
    }

    pub(super) fn retained_rna_read_hit_heap_from_rows(
        rows: &[RnaReadInterpretationHit],
    ) -> BinaryHeap<RetainedRnaReadHit> {
        let mut heap = BinaryHeap::<RetainedRnaReadHit>::new();
        for hit in rows.iter().cloned() {
            Self::retain_top_rna_read_hit(&mut heap, hit);
        }
        heap
    }

    pub(super) fn make_rna_read_top_hit_preview(
        hit: &RnaReadInterpretationHit,
    ) -> RnaReadTopHitPreview {
        let sequence_preview = hit.sequence.chars().take(80).collect::<String>();
        let (
            aligned,
            best_alignment_mode,
            best_alignment_transcript_id,
            best_alignment_transcript_label,
            best_alignment_strand,
            best_alignment_target_start_1based,
            best_alignment_target_end_1based,
            best_alignment_identity_fraction,
            best_alignment_query_coverage_fraction,
            best_alignment_score,
        ) = if let Some(best) = &hit.best_mapping {
            (
                true,
                best.alignment_mode.as_str().to_string(),
                best.transcript_id.clone(),
                best.transcript_label.clone(),
                best.strand.clone(),
                best.target_start_1based,
                best.target_end_1based,
                best.identity_fraction,
                best.query_coverage_fraction,
                best.score,
            )
        } else {
            (
                false,
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                0,
                0,
                0.0,
                0.0,
                0,
            )
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
            aligned,
            best_alignment_mode,
            best_alignment_transcript_id,
            best_alignment_transcript_label,
            best_alignment_strand,
            best_alignment_target_start_1based,
            best_alignment_target_end_1based,
            best_alignment_identity_fraction,
            best_alignment_query_coverage_fraction,
            best_alignment_score,
            secondary_mapping_count: hit.secondary_mappings.len(),
            read_length_bp: hit.read_length_bp,
            sequence: hit.sequence.clone(),
            sequence_preview,
        }
    }

    pub(super) fn retain_top_rna_read_preview_hit(
        retained_hits: &mut BinaryHeap<RetainedRnaReadPreviewHit>,
        hit: &RnaReadInterpretationHit,
    ) {
        let rank = Self::rna_read_retention_rank(hit);
        let preview = Self::make_rna_read_top_hit_preview(hit);
        if retained_hits.len() < RNA_READ_PROGRESS_TOP_HITS_PREVIEW_MAX {
            retained_hits.push(RetainedRnaReadPreviewHit { rank, preview });
            return;
        }
        let should_replace = retained_hits
            .peek()
            .map(|worst| rank > worst.rank)
            .unwrap_or(true);
        if should_replace {
            let _ = retained_hits.pop();
            retained_hits.push(RetainedRnaReadPreviewHit { rank, preview });
        }
    }

    pub(super) fn retained_rna_read_preview_heap_from_rows(
        rows: &[RnaReadInterpretationHit],
    ) -> BinaryHeap<RetainedRnaReadPreviewHit> {
        let mut heap = BinaryHeap::<RetainedRnaReadPreviewHit>::new();
        for hit in rows {
            Self::retain_top_rna_read_preview_hit(&mut heap, hit);
        }
        heap
    }

    pub(super) fn apply_rna_read_report_mode_to_hits(
        report_mode: RnaReadReportMode,
        hits: Vec<RnaReadInterpretationHit>,
    ) -> Vec<RnaReadInterpretationHit> {
        match report_mode {
            RnaReadReportMode::Full => hits,
            RnaReadReportMode::SeedPassedOnly => hits
                .into_iter()
                .filter(|hit| hit.passed_seed_filter)
                .collect::<Vec<_>>(),
        }
    }

    pub(super) fn collect_rna_read_top_hit_previews(
        retained_hits: &BinaryHeap<RetainedRnaReadPreviewHit>,
    ) -> Vec<RnaReadTopHitPreview> {
        let mut rows = retained_hits.iter().cloned().collect::<Vec<_>>();
        rows.sort_by(|left, right| right.rank.cmp(&left.rank));
        rows.into_iter().map(|row| row.preview).collect::<Vec<_>>()
    }

    fn mapped_support_indices_for_mapping_by_genomic_span(
        mapping: &RnaReadMappingHit,
        splicing: &SplicingExpertView,
    ) -> (Vec<usize>, Vec<usize>) {
        let span_start = mapping.target_start_1based.min(mapping.target_end_1based);
        let span_end = mapping.target_start_1based.max(mapping.target_end_1based);
        let mut exon_indices = Vec::<usize>::new();
        for (idx, exon) in splicing.unique_exons.iter().enumerate() {
            if span_start <= exon.end_1based && span_end >= exon.start_1based {
                exon_indices.push(idx);
            }
        }
        let mut junction_indices = Vec::<usize>::new();
        for (idx, junction) in splicing.junctions.iter().enumerate() {
            let donor = junction.donor_1based.min(junction.acceptor_1based);
            let acceptor = junction.donor_1based.max(junction.acceptor_1based);
            if span_start <= donor && span_end >= acceptor {
                junction_indices.push(idx);
            }
        }
        (exon_indices, junction_indices)
    }

    fn mapped_support_indices_for_mapping_by_transcript_offsets(
        mapping: &RnaReadMappingHit,
        splicing: &SplicingExpertView,
    ) -> Option<(Vec<usize>, Vec<usize>)> {
        if mapping.target_end_offset_0based_exclusive <= mapping.target_start_offset_0based {
            return None;
        }
        let Some(transcript) = splicing
            .transcripts
            .iter()
            .find(|row| row.transcript_feature_id == mapping.transcript_feature_id)
        else {
            return None;
        };
        if transcript.exons.is_empty() {
            return None;
        }
        let exon_index = splicing
            .unique_exons
            .iter()
            .enumerate()
            .map(|(idx, exon)| ((exon.start_1based, exon.end_1based), idx))
            .collect::<HashMap<_, _>>();
        let junction_index = splicing
            .junctions
            .iter()
            .enumerate()
            .map(|(idx, junction)| ((junction.donor_1based, junction.acceptor_1based), idx))
            .collect::<HashMap<_, _>>();
        let ordered_exons = if transcript.strand.trim() == "-" {
            transcript.exons.iter().rev().collect::<Vec<_>>()
        } else {
            transcript.exons.iter().collect::<Vec<_>>()
        };
        let aligned_start = mapping.target_start_offset_0based;
        let aligned_end = mapping.target_end_offset_0based_exclusive;
        let mut template_cursor = 0usize;
        let mut exon_indices = BTreeSet::<usize>::new();
        let mut junction_indices = BTreeSet::<usize>::new();
        for (idx, exon) in ordered_exons.iter().enumerate() {
            let exon_start_1based = exon.start_1based.min(exon.end_1based);
            let exon_end_1based = exon.start_1based.max(exon.end_1based);
            let exon_len = exon_end_1based
                .saturating_sub(exon_start_1based)
                .saturating_add(1);
            if exon_len == 0 {
                continue;
            }
            let exon_offset_start = template_cursor;
            let exon_offset_end = template_cursor.saturating_add(exon_len);
            if aligned_start < exon_offset_end && aligned_end > exon_offset_start {
                if let Some(exon_idx) = exon_index.get(&(exon_start_1based, exon_end_1based)) {
                    exon_indices.insert(*exon_idx);
                }
            }
            if idx + 1 < ordered_exons.len() {
                let boundary_offset = exon_offset_end;
                if aligned_start < boundary_offset && aligned_end > boundary_offset {
                    let next_exon = ordered_exons[idx + 1];
                    let next_start_1based = next_exon.start_1based.min(next_exon.end_1based);
                    let next_end_1based = next_exon.start_1based.max(next_exon.end_1based);
                    let (donor_1based, acceptor_1based) = if exon_start_1based <= next_start_1based
                    {
                        (exon_end_1based, next_start_1based)
                    } else {
                        (next_end_1based, exon_start_1based)
                    };
                    if let Some(junction_idx) = junction_index.get(&(donor_1based, acceptor_1based))
                    {
                        junction_indices.insert(*junction_idx);
                    }
                }
            }
            template_cursor = exon_offset_end;
            if template_cursor >= aligned_end {
                break;
            }
        }
        Some((
            exon_indices.into_iter().collect::<Vec<_>>(),
            junction_indices.into_iter().collect::<Vec<_>>(),
        ))
    }

    fn mapped_support_indices_for_mapping(
        mapping: &RnaReadMappingHit,
        splicing: &SplicingExpertView,
    ) -> (Vec<usize>, Vec<usize>) {
        Self::mapped_support_indices_for_mapping_by_transcript_offsets(mapping, splicing)
            .unwrap_or_else(|| {
                Self::mapped_support_indices_for_mapping_by_genomic_span(mapping, splicing)
            })
    }

    fn collect_mapped_support_attribution_rows(
        mapping: &RnaReadMappingHit,
        splicing: &SplicingExpertView,
    ) -> (
        Vec<RnaReadMappedSupportExonAttribution>,
        Vec<RnaReadMappedSupportJunctionAttribution>,
    ) {
        let (exon_indices, junction_indices) =
            Self::mapped_support_indices_for_mapping(mapping, splicing);
        let mapped_exon_support = exon_indices
            .into_iter()
            .filter_map(|idx| {
                let exon = splicing.unique_exons.get(idx)?;
                Some(RnaReadMappedSupportExonAttribution {
                    start_1based: exon.start_1based,
                    end_1based: exon.end_1based,
                })
            })
            .collect::<Vec<_>>();
        let mapped_junction_support = junction_indices
            .into_iter()
            .filter_map(|idx| {
                let junction = splicing.junctions.get(idx)?;
                Some(RnaReadMappedSupportJunctionAttribution {
                    donor_1based: junction.donor_1based,
                    acceptor_1based: junction.acceptor_1based,
                })
            })
            .collect::<Vec<_>>();
        (mapped_exon_support, mapped_junction_support)
    }

    fn format_mapped_exon_support_for_tsv(rows: &[RnaReadMappedSupportExonAttribution]) -> String {
        rows.iter()
            .map(|row| format!("{}..{}", row.start_1based, row.end_1based))
            .collect::<Vec<_>>()
            .join(";")
    }

    fn format_mapped_junction_support_for_tsv(
        rows: &[RnaReadMappedSupportJunctionAttribution],
    ) -> String {
        rows.iter()
            .map(|row| format!("{}->{}", row.donor_1based, row.acceptor_1based))
            .collect::<Vec<_>>()
            .join(";")
    }

    pub(super) fn accumulate_support_counts_for_mapping(
        mapping: &RnaReadMappingHit,
        splicing: &SplicingExpertView,
        exon_counts: &mut [usize],
        junction_counts: &mut [usize],
    ) {
        let (exon_indices, junction_indices) =
            Self::mapped_support_indices_for_mapping(mapping, splicing);
        for idx in exon_indices {
            if let Some(count) = exon_counts.get_mut(idx) {
                *count = count.saturating_add(1);
            }
        }
        for idx in junction_indices {
            if let Some(count) = junction_counts.get_mut(idx) {
                *count = count.saturating_add(1);
            }
        }
    }

    pub(super) fn build_rna_read_support_frequencies_from_counts(
        splicing: &SplicingExpertView,
        exon_counts: &[usize],
        junction_counts: &[usize],
        aligned_reads: usize,
    ) -> (
        Vec<RnaReadExonSupportFrequency>,
        Vec<RnaReadJunctionSupportFrequency>,
    ) {
        let denom = aligned_reads as f64;
        let exon = splicing
            .unique_exons
            .iter()
            .enumerate()
            .map(|(idx, exon)| {
                let support_read_count = exon_counts[idx];
                let support_fraction = if aligned_reads == 0 {
                    0.0
                } else {
                    support_read_count as f64 / denom
                };
                RnaReadExonSupportFrequency {
                    start_1based: exon.start_1based,
                    end_1based: exon.end_1based,
                    support_read_count,
                    support_fraction,
                }
            })
            .collect::<Vec<_>>();
        let junction = splicing
            .junctions
            .iter()
            .enumerate()
            .map(|(idx, j)| {
                let support_read_count = junction_counts[idx];
                let support_fraction = if aligned_reads == 0 {
                    0.0
                } else {
                    support_read_count as f64 / denom
                };
                RnaReadJunctionSupportFrequency {
                    donor_1based: j.donor_1based,
                    acceptor_1based: j.acceptor_1based,
                    support_read_count,
                    support_fraction,
                }
            })
            .collect::<Vec<_>>();
        (exon, junction)
    }

    pub(super) fn make_transcript_template(
        dna: &DNAsequence,
        lane: &SplicingTranscriptLane,
        kmer_len: usize,
    ) -> SplicingTranscriptTemplate {
        let forward = dna.forward_bytes();
        let is_reverse = lane.strand.trim() == "-";
        let mut sequence = Vec::<u8>::new();
        let mut genomic_positions_1based = Vec::<usize>::new();

        if is_reverse {
            for exon in lane.exons.iter().rev() {
                for pos_1based in (exon.start_1based..=exon.end_1based).rev() {
                    let idx = pos_1based.saturating_sub(1);
                    if idx >= forward.len() {
                        continue;
                    }
                    sequence.push(Self::complement_nucleotide_base(forward[idx]));
                    genomic_positions_1based.push(pos_1based);
                }
            }
        } else {
            for exon in &lane.exons {
                for pos_1based in exon.start_1based..=exon.end_1based {
                    let idx = pos_1based.saturating_sub(1);
                    if idx >= forward.len() {
                        continue;
                    }
                    sequence.push(Self::normalize_nucleotide_base(forward[idx]));
                    genomic_positions_1based.push(pos_1based);
                }
            }
        }

        let mut kmer_positions: HashMap<u32, Vec<usize>> = HashMap::new();
        if kmer_len > 0 && sequence.len() >= kmer_len {
            for start in 0..=sequence.len() - kmer_len {
                let window = &sequence[start..start + kmer_len];
                let Some(bits) = Self::encode_kmer_bits(window) else {
                    continue;
                };
                kmer_positions.entry(bits).or_default().push(start);
            }
        }

        SplicingTranscriptTemplate {
            transcript_feature_id: lane.transcript_feature_id,
            transcript_id: lane.transcript_id.clone(),
            transcript_label: lane.label.clone(),
            strand: lane.strand.clone(),
            sequence,
            genomic_positions_1based,
            kmer_positions,
        }
    }

    pub(super) fn build_exon_position_ordinal_map(
        exon_summaries: &[SplicingExonSummary],
    ) -> HashMap<usize, (usize, usize, usize)> {
        let mut map = HashMap::<usize, (usize, usize, usize)>::new();
        for (idx, exon) in exon_summaries.iter().enumerate() {
            let ordinal = idx + 1;
            for pos in exon.start_1based..=exon.end_1based {
                map.insert(pos, (ordinal, exon.start_1based, exon.end_1based));
            }
        }
        map
    }

    pub(super) fn collect_seed_support_exon_summaries(
        transcript_lanes: &[SplicingTranscriptLane],
    ) -> Vec<SplicingExonSummary> {
        let mut exon_support = HashMap::<(usize, usize), HashSet<usize>>::new();
        for lane in transcript_lanes {
            for exon in &lane.exons {
                exon_support
                    .entry((exon.start_1based, exon.end_1based))
                    .or_default()
                    .insert(lane.transcript_feature_id);
            }
        }
        let transcript_count = transcript_lanes.len().max(1);
        let mut exons = exon_support
            .into_iter()
            .map(
                |((start_1based, end_1based), supporting_transcripts)| SplicingExonSummary {
                    start_1based,
                    end_1based,
                    support_transcript_count: supporting_transcripts.len(),
                    constitutive: supporting_transcripts.len() == transcript_count,
                },
            )
            .collect::<Vec<_>>();
        exons.sort_by(|left, right| {
            left.start_1based
                .cmp(&right.start_1based)
                .then(left.end_1based.cmp(&right.end_1based))
        });
        exons
    }

    pub(super) fn build_seed_support_indexes(
        exon_summaries: &[SplicingExonSummary],
        templates: &[SplicingTranscriptTemplate],
        kmer_len: usize,
    ) -> (
        HashMap<u32, Vec<usize>>,
        HashMap<u32, Vec<(usize, usize)>>,
        Vec<TranscriptExonPathModel>,
        Vec<RnaReadTransitionSupportRow>,
    ) {
        if kmer_len == 0 {
            return (HashMap::new(), HashMap::new(), vec![], vec![]);
        }
        let exon_position_ordinal = Self::build_exon_position_ordinal_map(exon_summaries);
        let mut seed_to_exons = HashMap::<u32, HashSet<usize>>::new();
        let mut seed_to_transitions = HashMap::<u32, HashSet<(usize, usize)>>::new();
        let mut transcript_models = Vec::<TranscriptExonPathModel>::new();
        let mut transition_catalog = BTreeSet::<(usize, usize)>::new();

        for template in templates {
            let mut transcript_exons = Vec::<usize>::new();
            for pos in &template.genomic_positions_1based {
                let Some((ordinal, _, _)) = exon_position_ordinal.get(pos).copied() else {
                    continue;
                };
                if transcript_exons.last().copied() != Some(ordinal) {
                    transcript_exons.push(ordinal);
                }
            }
            if transcript_exons.is_empty() {
                continue;
            }
            let transcript_transitions = transcript_exons
                .windows(2)
                .map(|pair| (pair[0], pair[1]))
                .collect::<Vec<_>>();
            for transition in &transcript_transitions {
                transition_catalog.insert(*transition);
            }
            transcript_models.push(TranscriptExonPathModel {
                transcript_feature_id: template.transcript_feature_id,
                transcript_id: template.transcript_id.clone(),
                transcript_label: template.transcript_label.clone(),
                strand: template.strand.clone(),
                exon_ordinals: transcript_exons.clone(),
                transitions: transcript_transitions.clone(),
            });

            if template.sequence.len() < kmer_len {
                continue;
            }
            for start in 0..=template.sequence.len() - kmer_len {
                let window = &template.sequence[start..start + kmer_len];
                let Some(bits) = Self::encode_kmer_bits(window) else {
                    continue;
                };
                let mut touched_exons = Vec::<usize>::new();
                for offset in start..start + kmer_len {
                    let Some(pos_1based) = template.genomic_positions_1based.get(offset).copied()
                    else {
                        continue;
                    };
                    let Some((ordinal, _, _)) = exon_position_ordinal.get(&pos_1based).copied()
                    else {
                        continue;
                    };
                    if touched_exons.last().copied() != Some(ordinal) {
                        touched_exons.push(ordinal);
                    }
                }
                if touched_exons.is_empty() {
                    continue;
                }
                seed_to_exons
                    .entry(bits)
                    .or_default()
                    .extend(touched_exons.iter().copied());
                if touched_exons.len() >= 2 {
                    let row = seed_to_transitions.entry(bits).or_default();
                    for pair in touched_exons.windows(2) {
                        row.insert((pair[0], pair[1]));
                    }
                }
            }
        }

        let mut seed_to_exons_vec = seed_to_exons
            .into_iter()
            .map(|(bits, ordinals)| {
                let mut rows = ordinals.into_iter().collect::<Vec<_>>();
                rows.sort_unstable();
                (bits, rows)
            })
            .collect::<HashMap<_, _>>();
        for rows in seed_to_exons_vec.values_mut() {
            rows.dedup();
        }

        let mut seed_to_transitions_vec = seed_to_transitions
            .into_iter()
            .map(|(bits, transitions)| {
                let mut rows = transitions.into_iter().collect::<Vec<_>>();
                rows.sort_unstable_by(|left, right| {
                    left.0.cmp(&right.0).then(left.1.cmp(&right.1))
                });
                (bits, rows)
            })
            .collect::<HashMap<_, _>>();
        for rows in seed_to_transitions_vec.values_mut() {
            rows.dedup();
        }

        transcript_models.sort_by(|left, right| {
            left.transcript_feature_id
                .cmp(&right.transcript_feature_id)
                .then_with(|| left.transcript_id.cmp(&right.transcript_id))
        });

        let transition_support_rows = transition_catalog
            .into_iter()
            .map(|(from_ord, to_ord)| {
                let from_exon = exon_summaries.get(from_ord.saturating_sub(1));
                let to_exon = exon_summaries.get(to_ord.saturating_sub(1));
                RnaReadTransitionSupportRow {
                    from_exon_ordinal: from_ord,
                    to_exon_ordinal: to_ord,
                    from_start_1based: from_exon.map(|row| row.start_1based).unwrap_or(0),
                    from_end_1based: from_exon.map(|row| row.end_1based).unwrap_or(0),
                    to_start_1based: to_exon.map(|row| row.start_1based).unwrap_or(0),
                    to_end_1based: to_exon.map(|row| row.end_1based).unwrap_or(0),
                    support_read_count: 0,
                }
            })
            .collect::<Vec<_>>();

        (
            seed_to_exons_vec,
            seed_to_transitions_vec,
            transcript_models,
            transition_support_rows,
        )
    }

    pub(super) fn collect_sparse_target_transcript_lanes(
        dna: &DNAsequence,
        seed_feature_id: usize,
        scope: SplicingScopePreset,
        target_feature_strand: &str,
        target_gene_ids: &[String],
        already_present_feature_ids: &HashSet<usize>,
    ) -> Result<(Vec<SplicingTranscriptLane>, Vec<String>, Vec<String>), EngineError> {
        if target_gene_ids.is_empty() {
            return Ok((vec![], vec![], vec![]));
        }
        let features = dna.features();
        let target_feature = features.get(seed_feature_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!("Feature id '{}' was not found in sequence", seed_feature_id),
        })?;
        let target_is_reverse = feature_is_reverse(target_feature);
        let restrict_to_target_strand = scope.restrict_to_target_strand();

        let mut requested = target_gene_ids
            .iter()
            .map(|raw| raw.trim())
            .filter(|raw| !raw.is_empty())
            .map(|raw| raw.to_string())
            .collect::<Vec<_>>();
        requested.sort_by(|left, right| left.to_ascii_lowercase().cmp(&right.to_ascii_lowercase()));
        requested.dedup_by(|left, right| left.eq_ignore_ascii_case(right));
        let requested_lower = requested
            .iter()
            .map(|gene| gene.to_ascii_lowercase())
            .collect::<HashSet<_>>();

        let mut matched_requested = HashSet::<String>::new();
        let mut lanes = Vec::<SplicingTranscriptLane>::new();
        for (idx, feature) in features.iter().enumerate() {
            if !Self::is_splicing_transcript_feature(feature) {
                continue;
            }
            let group = Self::splicing_group_label(feature, idx);
            if !requested_lower.contains(&group.to_ascii_lowercase()) {
                continue;
            }
            matched_requested.insert(group.to_ascii_lowercase());
            if already_present_feature_ids.contains(&idx) {
                continue;
            }
            let is_reverse = feature_is_reverse(feature);
            if restrict_to_target_strand && is_reverse != target_is_reverse {
                continue;
            }
            let mut exon_ranges = vec![];
            collect_location_ranges_usize(&feature.location, &mut exon_ranges);
            if exon_ranges.is_empty() {
                let (from, to) = feature.location.find_bounds().map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not parse transcript range: {e}"),
                })?;
                if from >= 0 && to >= 0 {
                    exon_ranges.push((from as usize, to as usize));
                }
            }
            exon_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            exon_ranges.retain(|(start, end)| end > start);
            if exon_ranges.is_empty() {
                continue;
            }
            let mut introns = Vec::<(usize, usize)>::new();
            for pair in exon_ranges.windows(2) {
                let left = pair[0];
                let right = pair[1];
                if right.0 > left.1 {
                    introns.push((left.1, right.0));
                }
            }
            let exon_cds_phases =
                Self::exon_cds_phases_for_transcript(feature, &exon_ranges, is_reverse);
            lanes.push(SplicingTranscriptLane {
                transcript_feature_id: idx,
                transcript_id: Self::feature_transcript_id(feature, idx),
                label: Self::feature_display_label(feature, idx),
                strand: if is_reverse {
                    "-".to_string()
                } else {
                    "+".to_string()
                },
                exons: Self::range_vec_to_splicing(exon_ranges),
                exon_cds_phases,
                introns: Self::range_vec_to_splicing(introns),
                has_target_feature: idx == seed_feature_id,
            });
        }
        lanes.sort_by(|left, right| {
            left.transcript_id
                .cmp(&right.transcript_id)
                .then(left.transcript_feature_id.cmp(&right.transcript_feature_id))
        });
        let mut matched_gene_ids = requested
            .iter()
            .filter(|gene| matched_requested.contains(&gene.to_ascii_lowercase()))
            .cloned()
            .collect::<Vec<_>>();
        let mut missing_gene_ids = requested
            .iter()
            .filter(|gene| !matched_requested.contains(&gene.to_ascii_lowercase()))
            .cloned()
            .collect::<Vec<_>>();
        matched_gene_ids
            .sort_by(|left, right| left.to_ascii_lowercase().cmp(&right.to_ascii_lowercase()));
        missing_gene_ids
            .sort_by(|left, right| left.to_ascii_lowercase().cmp(&right.to_ascii_lowercase()));
        if restrict_to_target_strand && target_feature_strand.trim().is_empty() {
            missing_gene_ids = requested;
            lanes.clear();
            matched_gene_ids.clear();
        }
        Ok((lanes, matched_gene_ids, missing_gene_ids))
    }

    pub(super) fn build_isoform_support_accumulators(
        transcript_models: &[TranscriptExonPathModel],
    ) -> BTreeMap<String, IsoformSupportAccumulator> {
        let mut out = BTreeMap::<String, IsoformSupportAccumulator>::new();
        for model in transcript_models {
            out.insert(
                model.transcript_id.clone(),
                IsoformSupportAccumulator {
                    transcript_feature_id: model.transcript_feature_id,
                    transcript_id: model.transcript_id.clone(),
                    transcript_label: model.transcript_label.clone(),
                    strand: model.strand.clone(),
                    exon_count: model.exon_ordinals.len(),
                    expected_transition_count: model.transitions.len(),
                    ..IsoformSupportAccumulator::default()
                },
            );
        }
        out
    }

    pub(super) fn update_isoform_support_accumulator(
        accumulators: &mut BTreeMap<String, IsoformSupportAccumulator>,
        assigned_transcript_id: &str,
        assigned_strand: &str,
        passed_seed_filter: bool,
        seed_median_gap: f64,
        transition_confirmed: usize,
        transition_total: usize,
        seed_hit_fraction: f64,
        weighted_seed_hit_fraction: f64,
        chain_preferred_strand: &str,
        competing_opposite_strand: bool,
        ambiguous_strand_tie: bool,
        read_supported_transitions: &HashSet<(usize, usize)>,
        transcript_models_by_id: &HashMap<String, TranscriptExonPathModel>,
    ) {
        let Some(row) = accumulators.get_mut(assigned_transcript_id) else {
            return;
        };
        row.reads_assigned = row.reads_assigned.saturating_add(1);
        if passed_seed_filter {
            row.reads_seed_passed = row.reads_seed_passed.saturating_add(1);
        }
        if seed_median_gap.is_finite() && seed_median_gap >= 0.0 {
            row.seed_gap_sum += seed_median_gap;
            row.seed_gap_count = row.seed_gap_count.saturating_add(1);
        }
        if passed_seed_filter && transition_total > 0 {
            row.confirmed_transition_fraction_sum +=
                transition_confirmed as f64 / transition_total as f64;
            row.confirmed_transition_fraction_count =
                row.confirmed_transition_fraction_count.saturating_add(1);
        }
        row.best_seed_hit_fraction = row.best_seed_hit_fraction.max(seed_hit_fraction);
        row.best_weighted_seed_hit_fraction = row
            .best_weighted_seed_hit_fraction
            .max(weighted_seed_hit_fraction);
        if !chain_preferred_strand.is_empty() && chain_preferred_strand == assigned_strand {
            row.reads_chain_same_strand = row.reads_chain_same_strand.saturating_add(1);
        }
        if competing_opposite_strand {
            row.reads_with_opposite_strand_competition =
                row.reads_with_opposite_strand_competition.saturating_add(1);
        }
        if ambiguous_strand_tie {
            row.reads_ambiguous_strand_ties = row.reads_ambiguous_strand_ties.saturating_add(1);
        }
        if passed_seed_filter {
            if let Some(model) = transcript_models_by_id.get(assigned_transcript_id) {
                for transition in &model.transitions {
                    if read_supported_transitions.contains(transition) {
                        row.transition_rows_supported.insert(*transition);
                    }
                }
            }
        }
    }

    pub(super) fn collect_isoform_support_rows(
        accumulators: &BTreeMap<String, IsoformSupportAccumulator>,
    ) -> Vec<RnaReadIsoformSupportRow> {
        let mut rows = accumulators
            .values()
            .map(|row| {
                let transition_rows_supported = row.transition_rows_supported.len();
                let transition_rows_supported_fraction = if row.expected_transition_count == 0 {
                    0.0
                } else {
                    transition_rows_supported as f64 / row.expected_transition_count as f64
                };
                let mean_seed_median_gap = if row.seed_gap_count == 0 {
                    -1.0
                } else {
                    row.seed_gap_sum / row.seed_gap_count as f64
                };
                let mean_confirmed_transition_fraction =
                    if row.confirmed_transition_fraction_count == 0 {
                        0.0
                    } else {
                        row.confirmed_transition_fraction_sum
                            / row.confirmed_transition_fraction_count as f64
                    };
                RnaReadIsoformSupportRow {
                    transcript_feature_id: row.transcript_feature_id,
                    transcript_id: row.transcript_id.clone(),
                    transcript_label: row.transcript_label.clone(),
                    strand: row.strand.clone(),
                    exon_count: row.exon_count,
                    expected_transition_count: row.expected_transition_count,
                    reads_assigned: row.reads_assigned,
                    reads_seed_passed: row.reads_seed_passed,
                    transition_rows_supported,
                    transition_rows_supported_fraction,
                    mean_seed_median_gap,
                    mean_confirmed_transition_fraction,
                    best_seed_hit_fraction: row.best_seed_hit_fraction,
                    best_weighted_seed_hit_fraction: row.best_weighted_seed_hit_fraction,
                    reads_chain_same_strand: row.reads_chain_same_strand,
                    reads_with_opposite_strand_competition: row
                        .reads_with_opposite_strand_competition,
                    reads_ambiguous_strand_ties: row.reads_ambiguous_strand_ties,
                }
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            right
                .reads_seed_passed
                .cmp(&left.reads_seed_passed)
                .then(right.reads_assigned.cmp(&left.reads_assigned))
                .then_with(|| {
                    right
                        .mean_confirmed_transition_fraction
                        .partial_cmp(&left.mean_confirmed_transition_fraction)
                        .unwrap_or(Ordering::Equal)
                })
                .then_with(|| {
                    right
                        .best_weighted_seed_hit_fraction
                        .partial_cmp(&left.best_weighted_seed_hit_fraction)
                        .unwrap_or(Ordering::Equal)
                })
                .then_with(|| left.transcript_id.cmp(&right.transcript_id))
        });
        rows
    }

    pub(super) fn collect_mapped_isoform_support_rows(
        hits: &[RnaReadInterpretationHit],
    ) -> Vec<RnaReadMappedIsoformSupportRow> {
        #[derive(Default)]
        struct Accumulator {
            transcript_feature_id: usize,
            transcript_id: String,
            transcript_label: String,
            strand: String,
            aligned_read_count: usize,
            msa_eligible_read_count: usize,
            identity_sum: f64,
            query_coverage_sum: f64,
            best_alignment_score: isize,
            secondary_mapping_total: usize,
        }

        let mut accumulators = BTreeMap::<String, Accumulator>::new();
        for hit in hits {
            let Some(mapping) = hit.best_mapping.as_ref() else {
                continue;
            };
            let row = accumulators
                .entry(mapping.transcript_id.clone())
                .or_insert_with(|| Accumulator {
                    transcript_feature_id: mapping.transcript_feature_id,
                    transcript_id: mapping.transcript_id.clone(),
                    transcript_label: mapping.transcript_label.clone(),
                    strand: mapping.strand.clone(),
                    best_alignment_score: isize::MIN,
                    ..Accumulator::default()
                });
            row.aligned_read_count = row.aligned_read_count.saturating_add(1);
            if hit.msa_eligible {
                row.msa_eligible_read_count = row.msa_eligible_read_count.saturating_add(1);
            }
            row.identity_sum += mapping.identity_fraction;
            row.query_coverage_sum += mapping.query_coverage_fraction;
            row.best_alignment_score = row.best_alignment_score.max(mapping.score);
            row.secondary_mapping_total = row
                .secondary_mapping_total
                .saturating_add(hit.secondary_mappings.len());
        }

        let mut rows = accumulators
            .into_values()
            .map(|row| {
                let denom = row.aligned_read_count.max(1) as f64;
                RnaReadMappedIsoformSupportRow {
                    transcript_feature_id: row.transcript_feature_id,
                    transcript_id: row.transcript_id,
                    transcript_label: row.transcript_label,
                    strand: row.strand,
                    aligned_read_count: row.aligned_read_count,
                    msa_eligible_read_count: row.msa_eligible_read_count,
                    mean_identity_fraction: row.identity_sum / denom,
                    mean_query_coverage_fraction: row.query_coverage_sum / denom,
                    best_alignment_score: if row.best_alignment_score == isize::MIN {
                        0
                    } else {
                        row.best_alignment_score
                    },
                    secondary_mapping_total: row.secondary_mapping_total,
                }
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            right
                .aligned_read_count
                .cmp(&left.aligned_read_count)
                .then(
                    right
                        .msa_eligible_read_count
                        .cmp(&left.msa_eligible_read_count),
                )
                .then_with(|| {
                    right
                        .mean_identity_fraction
                        .partial_cmp(&left.mean_identity_fraction)
                        .unwrap_or(Ordering::Equal)
                })
                .then_with(|| {
                    right
                        .mean_query_coverage_fraction
                        .partial_cmp(&left.mean_query_coverage_fraction)
                        .unwrap_or(Ordering::Equal)
                })
                .then(right.best_alignment_score.cmp(&left.best_alignment_score))
                .then_with(|| left.transcript_id.cmp(&right.transcript_id))
        });
        rows
    }

    pub(super) fn transition_gate_passes(
        confirmed_transitions: usize,
        total_transitions: usize,
        seed_filter: &RnaReadSeedFilterConfig,
    ) -> bool {
        if total_transitions == 0 {
            return true;
        }
        let count_pass = confirmed_transitions >= seed_filter.min_confirmed_exon_transitions;
        let fraction_pass = (confirmed_transitions as f64 / total_transitions as f64)
            + f64::EPSILON
            >= seed_filter.min_transition_support_fraction;
        count_pass && fraction_pass
    }

    pub(super) fn prefer_transcript_score(
        candidate: TranscriptSupportScore<'_>,
        current: TranscriptSupportScore<'_>,
        chain_preferred_strand: Option<&str>,
    ) -> bool {
        if candidate.transition_hits != current.transition_hits {
            return candidate.transition_hits > current.transition_hits;
        }
        if candidate.exon_hits != current.exon_hits {
            return candidate.exon_hits > current.exon_hits;
        }
        if candidate.model.strand != current.model.strand {
            if let Some(preferred) = chain_preferred_strand {
                let candidate_pref = candidate.model.strand == preferred;
                let current_pref = current.model.strand == preferred;
                if candidate_pref != current_pref {
                    return candidate_pref;
                }
            }
        }
        candidate.model.transcript_feature_id < current.model.transcript_feature_id
    }

    pub(super) fn infer_read_exon_path(
        transcript_models: &[TranscriptExonPathModel],
        supported_exons: &HashSet<usize>,
        supported_transitions: &HashSet<(usize, usize)>,
        seed_chain_transcript_id: &str,
    ) -> ReadExonPathInference {
        let chain_preferred_strand = transcript_models
            .iter()
            .find(|model| model.transcript_id == seed_chain_transcript_id)
            .map(|model| model.strand.as_str());
        let mut best: Option<TranscriptSupportScore<'_>> = None;
        let mut best_plus: Option<TranscriptSupportScore<'_>> = None;
        let mut best_minus: Option<TranscriptSupportScore<'_>> = None;
        for model in transcript_models {
            if model.exon_ordinals.is_empty() {
                continue;
            }
            let exon_hits = model
                .exon_ordinals
                .iter()
                .filter(|ordinal| supported_exons.contains(ordinal))
                .count();
            let transition_hits = model
                .transitions
                .iter()
                .filter(|pair| supported_transitions.contains(pair))
                .count();
            let candidate = TranscriptSupportScore {
                model,
                exon_hits,
                transition_hits,
            };
            let should_replace = match best {
                None => true,
                Some(current) => {
                    Self::prefer_transcript_score(candidate, current, chain_preferred_strand)
                }
            };
            if should_replace {
                best = Some(candidate);
            }

            if model.strand.trim() == "-" {
                let should_replace_minus = match best_minus {
                    None => true,
                    Some(current) => {
                        Self::prefer_transcript_score(candidate, current, chain_preferred_strand)
                    }
                };
                if should_replace_minus {
                    best_minus = Some(candidate);
                }
            } else {
                let should_replace_plus = match best_plus {
                    None => true,
                    Some(current) => {
                        Self::prefer_transcript_score(candidate, current, chain_preferred_strand)
                    }
                };
                if should_replace_plus {
                    best_plus = Some(candidate);
                }
            }
        }
        let Some(selected) = best else {
            return ReadExonPathInference {
                strand_diagnostics: RnaReadStrandAssignmentDiagnostics {
                    selected_reason: "no transcript model matched".to_string(),
                    chain_preferred_strand: chain_preferred_strand.unwrap_or("").to_string(),
                    ..RnaReadStrandAssignmentDiagnostics::default()
                },
                ..ReadExonPathInference::default()
            };
        };
        let model = selected.model;
        let confirmed_transitions = selected.transition_hits;
        if model.exon_ordinals.is_empty() {
            return ReadExonPathInference {
                strand_diagnostics: RnaReadStrandAssignmentDiagnostics {
                    selected_reason: "selected transcript had no exons".to_string(),
                    chain_preferred_strand: chain_preferred_strand.unwrap_or("").to_string(),
                    ..RnaReadStrandAssignmentDiagnostics::default()
                },
                ..ReadExonPathInference::default()
            };
        }
        let mut path = model.exon_ordinals[0].to_string();
        for pair in model.exon_ordinals.windows(2) {
            let edge = (pair[0], pair[1]);
            let sep = if supported_transitions.contains(&edge) {
                ':'
            } else {
                '-'
            };
            path.push(sep);
            path.push_str(pair[1].to_string().as_str());
        }
        let plus_best = best_plus;
        let minus_best = best_minus;
        let opposite_best = if model.strand.trim() == "-" {
            plus_best
        } else {
            minus_best
        };
        let competing_opposite_strand = opposite_best
            .map(|row| row.transition_hits > 0 || row.exon_hits > 0)
            .unwrap_or(false);
        let ambiguous_near_tie = opposite_best
            .map(|row| {
                row.transition_hits == selected.transition_hits
                    && row.exon_hits == selected.exon_hits
                    && (row.transition_hits > 0 || row.exon_hits > 0)
            })
            .unwrap_or(false);
        let selected_reason = if let Some(opposite) = opposite_best {
            if selected.transition_hits != opposite.transition_hits {
                format!(
                    "selected higher transition support {}>{}",
                    selected.transition_hits, opposite.transition_hits
                )
            } else if selected.exon_hits != opposite.exon_hits {
                format!(
                    "selected higher exon support {}>{}",
                    selected.exon_hits, opposite.exon_hits
                )
            } else if let Some(preferred) = chain_preferred_strand {
                if model.strand == preferred && opposite.model.strand != preferred {
                    format!("tie resolved by chain-preferred strand '{}'", preferred)
                } else {
                    format!(
                        "tie resolved by transcript feature id {}<={}",
                        model.transcript_feature_id, opposite.model.transcript_feature_id
                    )
                }
            } else {
                format!(
                    "tie resolved by transcript feature id {}<={}",
                    model.transcript_feature_id, opposite.model.transcript_feature_id
                )
            }
        } else {
            format!("only strand '{}' has indexed support", model.strand)
        };
        ReadExonPathInference {
            path,
            confirmed_transitions,
            total_transitions: model.transitions.len(),
            transcript_id: model.transcript_id.clone(),
            strand: model.strand.clone(),
            strand_diagnostics: RnaReadStrandAssignmentDiagnostics {
                selected_strand: model.strand.clone(),
                selected_reason,
                selected_transition_hits: selected.transition_hits,
                selected_exon_hits: selected.exon_hits,
                plus_best_transcript_id: plus_best
                    .map(|row| row.model.transcript_id.clone())
                    .unwrap_or_default(),
                plus_best_transition_hits: plus_best.map(|row| row.transition_hits).unwrap_or(0),
                plus_best_exon_hits: plus_best.map(|row| row.exon_hits).unwrap_or(0),
                minus_best_transcript_id: minus_best
                    .map(|row| row.model.transcript_id.clone())
                    .unwrap_or_default(),
                minus_best_transition_hits: minus_best.map(|row| row.transition_hits).unwrap_or(0),
                minus_best_exon_hits: minus_best.map(|row| row.exon_hits).unwrap_or(0),
                competing_opposite_strand,
                ambiguous_near_tie,
                chain_preferred_strand: chain_preferred_strand.unwrap_or("").to_string(),
            },
        }
    }

    fn align_read_to_template_dense_with_mode(
        read: &[u8],
        template: &SplicingTranscriptTemplate,
        config: &RnaReadAlignConfig,
        mode: RnaReadAlignmentMode,
    ) -> Option<RnaReadMappingHit> {
        let mut aligner = bio::alignment::pairwise::Aligner::new(-5, -1, &|a: u8, b: u8| {
            if a.eq_ignore_ascii_case(&b) {
                2i32
            } else {
                -3i32
            }
        });
        let alignment = match mode {
            RnaReadAlignmentMode::Local => aligner.local(read, &template.sequence),
            RnaReadAlignmentMode::Semiglobal => aligner.semiglobal(read, &template.sequence),
        };
        Self::rna_mapping_hit_from_alignment(read, template, config, mode, alignment)
    }

    fn effective_banded_k(seed_kmer_len: usize, read_len: usize, template_len: usize) -> usize {
        let min_len = read_len.min(template_len);
        if min_len <= 2 {
            return min_len.max(1);
        }
        seed_kmer_len.clamp(3, min_len)
    }

    fn align_read_to_template_banded_with_mode(
        read: &[u8],
        template: &SplicingTranscriptTemplate,
        config: &RnaReadAlignConfig,
        seed_kmer_len: usize,
        mode: RnaReadAlignmentMode,
    ) -> Option<RnaReadMappingHit> {
        let k = Self::effective_banded_k(seed_kmer_len, read.len(), template.sequence.len());
        let w = config.band_width_bp.max(1);
        let mut aligner = bio::alignment::pairwise::banded::Aligner::new(
            -5,
            -1,
            |a: u8, b: u8| {
                if a.eq_ignore_ascii_case(&b) {
                    2i32
                } else {
                    -3i32
                }
            },
            k,
            w,
        );
        let alignment = match mode {
            RnaReadAlignmentMode::Local => aligner.local(read, &template.sequence),
            RnaReadAlignmentMode::Semiglobal => aligner.semiglobal(read, &template.sequence),
        };
        Self::rna_mapping_hit_from_alignment(read, template, config, mode, alignment)
    }

    fn rna_mapping_hit_from_alignment(
        read: &[u8],
        template: &SplicingTranscriptTemplate,
        config: &RnaReadAlignConfig,
        mode: RnaReadAlignmentMode,
        alignment: bio::alignment::Alignment,
    ) -> Option<RnaReadMappingHit> {
        if alignment.operations.is_empty() {
            return None;
        }
        let query_span = alignment.xend.saturating_sub(alignment.xstart);
        if query_span == 0 {
            return None;
        }
        if alignment.yend <= alignment.ystart {
            return None;
        }
        let mut matches = 0usize;
        let mut mismatches = 0usize;
        let mut indels = 0usize;
        for op in &alignment.operations {
            match op {
                bio::alignment::AlignmentOperation::Match => matches += 1,
                bio::alignment::AlignmentOperation::Subst => mismatches += 1,
                bio::alignment::AlignmentOperation::Del => indels += 1,
                bio::alignment::AlignmentOperation::Ins => indels += 1,
                bio::alignment::AlignmentOperation::Xclip(_)
                | bio::alignment::AlignmentOperation::Yclip(_) => {}
            }
        }
        let aligned_columns = matches + mismatches + indels;
        if aligned_columns == 0 {
            return None;
        }
        let identity_fraction = matches as f64 / aligned_columns as f64;
        if identity_fraction + f64::EPSILON < config.min_identity_fraction {
            return None;
        }
        let query_coverage_fraction = if read.is_empty() {
            0.0
        } else {
            query_span as f64 / read.len() as f64
        };
        let target_start_1based = template
            .genomic_positions_1based
            .get(alignment.ystart)
            .copied()
            .unwrap_or(1);
        let target_end_1based = template
            .genomic_positions_1based
            .get(alignment.yend.saturating_sub(1))
            .copied()
            .unwrap_or(target_start_1based);
        let (target_start_1based, target_end_1based) = if target_start_1based <= target_end_1based {
            (target_start_1based, target_end_1based)
        } else {
            (target_end_1based, target_start_1based)
        };
        Some(RnaReadMappingHit {
            alignment_mode: mode,
            transcript_feature_id: template.transcript_feature_id,
            transcript_id: template.transcript_id.clone(),
            transcript_label: template.transcript_label.clone(),
            strand: template.strand.clone(),
            query_start_0based: alignment.xstart,
            query_end_0based_exclusive: alignment.xend,
            target_start_1based,
            target_end_1based,
            target_start_offset_0based: alignment.ystart,
            target_end_offset_0based_exclusive: alignment.yend,
            matches,
            mismatches,
            score: alignment.score as isize,
            identity_fraction,
            query_coverage_fraction,
        })
    }

    pub(super) fn align_read_to_template(
        read: &[u8],
        template: &SplicingTranscriptTemplate,
        config: &RnaReadAlignConfig,
        seed_kmer_len: usize,
    ) -> Option<RnaReadMappingHit> {
        if read.is_empty() || template.sequence.is_empty() {
            return None;
        }
        let semiglobal = Self::align_read_to_template_banded_with_mode(
            read,
            template,
            config,
            seed_kmer_len,
            RnaReadAlignmentMode::Semiglobal,
        )
        .or_else(|| {
            Self::align_read_to_template_dense_with_mode(
                read,
                template,
                config,
                RnaReadAlignmentMode::Semiglobal,
            )
        });
        let local = Self::align_read_to_template_banded_with_mode(
            read,
            template,
            config,
            seed_kmer_len,
            RnaReadAlignmentMode::Local,
        )
        .or_else(|| {
            Self::align_read_to_template_dense_with_mode(
                read,
                template,
                config,
                RnaReadAlignmentMode::Local,
            )
        });
        match (semiglobal, local) {
            (Some(left), Some(right)) => {
                if Self::compare_mapping_quality(&left, &right) != Ordering::Less {
                    Some(left)
                } else {
                    Some(right)
                }
            }
            (Some(row), None) | (None, Some(row)) => Some(row),
            (None, None) => None,
        }
    }

    fn compare_mapping_quality(left: &RnaReadMappingHit, right: &RnaReadMappingHit) -> Ordering {
        left.query_coverage_fraction
            .partial_cmp(&right.query_coverage_fraction)
            .unwrap_or(Ordering::Equal)
            .then_with(|| {
                left.identity_fraction
                    .partial_cmp(&right.identity_fraction)
                    .unwrap_or(Ordering::Equal)
            })
            .then(left.matches.cmp(&right.matches))
            .then(right.mismatches.cmp(&left.mismatches))
            .then(left.score.cmp(&right.score))
    }

    pub(super) fn align_read_to_templates(
        read: &[u8],
        templates: &[SplicingTranscriptTemplate],
        config: &RnaReadAlignConfig,
        seed_kmer_len: usize,
    ) -> (Option<RnaReadMappingHit>, Vec<RnaReadMappingHit>) {
        let mut mappings = templates
            .iter()
            .filter_map(|template| {
                Self::align_read_to_template(read, template, config, seed_kmer_len)
            })
            .collect::<Vec<_>>();
        mappings.sort_by(|left, right| {
            right
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
                .then(left.transcript_feature_id.cmp(&right.transcript_feature_id))
                .then(left.strand.cmp(&right.strand))
                .then(left.target_start_1based.cmp(&right.target_start_1based))
        });
        if mappings.is_empty() {
            return (None, vec![]);
        }
        let best = mappings.first().cloned();
        let secondary = mappings
            .into_iter()
            .skip(1)
            .take(config.max_secondary_mappings)
            .collect::<Vec<_>>();
        (best, secondary)
    }

    #[allow(dead_code)]
    pub(super) fn compute_rna_read_support_frequencies(
        hits: &[RnaReadInterpretationHit],
        splicing: &SplicingExpertView,
    ) -> (
        Vec<RnaReadExonSupportFrequency>,
        Vec<RnaReadJunctionSupportFrequency>,
    ) {
        let mut exon_counts = vec![0usize; splicing.unique_exons.len()];
        let mut junction_counts = vec![0usize; splicing.junctions.len()];
        let mut aligned_reads = 0usize;
        for hit in hits {
            let Some(mapping) = &hit.best_mapping else {
                continue;
            };
            aligned_reads = aligned_reads.saturating_add(1);
            Self::accumulate_support_counts_for_mapping(
                mapping,
                splicing,
                &mut exon_counts,
                &mut junction_counts,
            );
        }
        Self::build_rna_read_support_frequencies_from_counts(
            splicing,
            &exon_counts,
            &junction_counts,
            aligned_reads,
        )
    }

    pub(super) fn push_rna_read_report_result_message(
        &mut self,
        report: RnaReadInterpretationReport,
        result: &mut OpResult,
    ) -> Result<(), EngineError> {
        let origin_summary = if report.origin_class_counts.is_empty() {
            "origin_classes=none".to_string()
        } else {
            let parts = report
                .origin_class_counts
                .iter()
                .map(|(class, count)| format!("{class}:{count}"))
                .collect::<Vec<_>>();
            format!("origin_classes={}", parts.join(","))
        };
        let replaced = self
            .read_rna_read_report_store()
            .reports
            .contains_key(report.report_id.as_str());
        self.upsert_rna_read_report(report.clone())?;
        result.messages.push(format!(
            "{} RNA-read report '{}' (profile={}, origin_mode={}, report_mode={}, targets={}, reads={}, seed_passed={}, aligned={}, msa_eligible(retained)={}, {})",
            if replaced { "Updated" } else { "Created" },
            report.report_id,
            report.profile.as_str(),
            report.origin_mode.as_str(),
            report.report_mode.as_str(),
            report.target_gene_ids.len(),
            report.read_count_total,
            report.read_count_seed_passed,
            report.read_count_aligned,
            report.retained_count_msa_eligible,
            origin_summary
        ));
        Ok(())
    }

    pub fn compute_rna_read_report_with_progress(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
        profile: RnaReadInterpretationProfile,
        input_path: &str,
        input_format: RnaReadInputFormat,
        scope: SplicingScopePreset,
        origin_mode: RnaReadOriginMode,
        target_gene_ids: &[String],
        roi_seed_capture_enabled: bool,
        seed_filter: &RnaReadSeedFilterConfig,
        align_config: &RnaReadAlignConfig,
        report_id: Option<&str>,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<RnaReadInterpretationReport, EngineError> {
        let mut keep_running = || true;
        self.compute_rna_read_report_with_progress_and_cancel(
            seq_id,
            seed_feature_id,
            profile,
            input_path,
            input_format,
            scope,
            origin_mode,
            target_gene_ids,
            roi_seed_capture_enabled,
            seed_filter,
            align_config,
            report_id,
            on_progress,
            &mut keep_running,
        )
    }

    pub fn compute_rna_read_report_with_progress_and_cancel(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
        profile: RnaReadInterpretationProfile,
        input_path: &str,
        input_format: RnaReadInputFormat,
        scope: SplicingScopePreset,
        origin_mode: RnaReadOriginMode,
        target_gene_ids: &[String],
        roi_seed_capture_enabled: bool,
        seed_filter: &RnaReadSeedFilterConfig,
        align_config: &RnaReadAlignConfig,
        report_id: Option<&str>,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
        should_continue: &mut dyn FnMut() -> bool,
    ) -> Result<RnaReadInterpretationReport, EngineError> {
        self.interpret_rna_reads_report_with_progress(
            seq_id,
            seed_feature_id,
            profile,
            input_path,
            input_format,
            scope,
            origin_mode,
            target_gene_ids,
            roi_seed_capture_enabled,
            seed_filter,
            align_config,
            report_id,
            &RnaReadInterpretOptions::default(),
            on_progress,
            should_continue,
        )
    }

    pub fn compute_rna_read_report_with_runtime_options_and_progress_and_cancel(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
        profile: RnaReadInterpretationProfile,
        input_path: &str,
        input_format: RnaReadInputFormat,
        scope: SplicingScopePreset,
        origin_mode: RnaReadOriginMode,
        target_gene_ids: &[String],
        roi_seed_capture_enabled: bool,
        seed_filter: &RnaReadSeedFilterConfig,
        align_config: &RnaReadAlignConfig,
        report_id: Option<&str>,
        report_mode: RnaReadReportMode,
        checkpoint_path: Option<&str>,
        checkpoint_every_reads: usize,
        resume_from_checkpoint: bool,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
        should_continue: &mut dyn FnMut() -> bool,
    ) -> Result<RnaReadInterpretationReport, EngineError> {
        let options = RnaReadInterpretOptions {
            report_mode,
            checkpoint_path: checkpoint_path.map(|raw| raw.to_string()),
            checkpoint_every_reads,
            resume_from_checkpoint,
        };
        self.compute_rna_read_report_with_options_and_progress_and_cancel(
            seq_id,
            seed_feature_id,
            profile,
            input_path,
            input_format,
            scope,
            origin_mode,
            target_gene_ids,
            roi_seed_capture_enabled,
            seed_filter,
            align_config,
            report_id,
            &options,
            on_progress,
            should_continue,
        )
    }

    pub fn align_rna_read_report_with_progress(
        &self,
        report_id: &str,
        selection: RnaReadHitSelection,
        align_config_override: Option<RnaReadAlignConfig>,
        selected_record_indices: &[usize],
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<RnaReadInterpretationReport, EngineError> {
        let mut keep_running = || true;
        self.align_rna_read_report_with_progress_and_cancel(
            report_id,
            selection,
            align_config_override,
            selected_record_indices,
            on_progress,
            &mut keep_running,
        )
    }

    pub fn align_rna_read_report_with_progress_and_cancel(
        &self,
        report_id: &str,
        selection: RnaReadHitSelection,
        align_config_override: Option<RnaReadAlignConfig>,
        selected_record_indices: &[usize],
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
        should_continue: &mut dyn FnMut() -> bool,
    ) -> Result<RnaReadInterpretationReport, EngineError> {
        let mut report = self.get_rna_read_report(report_id)?;
        if !should_continue() {
            return Err(EngineError {
                code: ErrorCode::Internal,
                message: "RNA-read alignment phase cancelled before start".to_string(),
            });
        }
        let align_config = align_config_override.unwrap_or_else(|| report.align_config.clone());
        if align_config.band_width_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "align_band_width_bp must be > 0".to_string(),
            });
        }
        if !(0.0..=1.0).contains(&align_config.min_identity_fraction) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "align_min_identity_fraction must be within 0.0..=1.0".to_string(),
            });
        }
        let dna = self
            .state
            .sequences
            .get(&report.seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", report.seq_id),
            })?;
        let splicing =
            self.build_splicing_expert_view(&report.seq_id, report.seed_feature_id, report.scope)?;
        let mut transcript_lanes = splicing.transcripts.clone();
        if matches!(report.origin_mode, RnaReadOriginMode::MultiGeneSparse) {
            let existing_feature_ids = transcript_lanes
                .iter()
                .map(|lane| lane.transcript_feature_id)
                .collect::<HashSet<_>>();
            let (extra_lanes, _matched_genes, _missing_genes) =
                Self::collect_sparse_target_transcript_lanes(
                    dna,
                    report.seed_feature_id,
                    report.scope,
                    splicing.strand.as_str(),
                    &report.target_gene_ids,
                    &existing_feature_ids,
                )?;
            if !extra_lanes.is_empty() {
                transcript_lanes.extend(extra_lanes);
                transcript_lanes.sort_by(|left, right| {
                    left.transcript_id
                        .cmp(&right.transcript_id)
                        .then(left.transcript_feature_id.cmp(&right.transcript_feature_id))
                });
            }
        }
        let templates = transcript_lanes
            .iter()
            .map(|lane| Self::make_transcript_template(dna, lane, report.seed_filter.kmer_len))
            .filter(|template| !template.sequence.is_empty())
            .collect::<Vec<_>>();
        if templates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No transcript templates available for splicing scope '{}' on '{}'",
                    report.scope.as_str(),
                    report.seq_id
                ),
            });
        }
        let mut seed_index: HashSet<u32> = HashSet::new();
        for template in &templates {
            seed_index.extend(template.kmer_positions.keys().copied());
        }
        if seed_index.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: "No directional k-mer seeds could be generated from transcript templates"
                    .to_string(),
            });
        }
        let seed_catalog_rows =
            Self::collect_rna_seed_hash_catalog_rows(&templates, report.seed_filter.kmer_len);
        let mut seed_occurrence_counts = HashMap::<u32, usize>::new();
        for row in &seed_catalog_rows {
            *seed_occurrence_counts.entry(row.seed_bits).or_insert(0) += 1;
        }
        let seed_support_exons = Self::collect_seed_support_exon_summaries(&transcript_lanes);
        let (
            seed_to_exons,
            seed_to_transitions,
            transcript_exon_models,
            mut transition_support_rows,
        ) = Self::build_seed_support_indexes(
            &seed_support_exons,
            &templates,
            report.seed_filter.kmer_len,
        );
        let transcript_models_by_id = transcript_exon_models
            .iter()
            .map(|model| (model.transcript_id.clone(), model.clone()))
            .collect::<HashMap<_, _>>();
        let transition_row_index = transition_support_rows
            .iter()
            .enumerate()
            .map(|(idx, row)| ((row.from_exon_ordinal, row.to_exon_ordinal), idx))
            .collect::<HashMap<_, _>>();
        let seed_template_positions = Self::build_seed_template_position_index(&templates);
        let mut bins = Self::build_rna_read_seed_histogram_bins(dna.len());
        let histogram_index =
            Self::build_rna_read_seed_histogram_index(&templates, dna.len(), &bins);
        let explicit_record_filter = selected_record_indices
            .iter()
            .copied()
            .collect::<HashSet<_>>();
        let selected_indices = report
            .hits
            .iter()
            .enumerate()
            .filter_map(|(idx, hit)| {
                if !explicit_record_filter.is_empty() {
                    explicit_record_filter
                        .contains(&hit.record_index)
                        .then_some(idx)
                } else {
                    Self::include_rna_read_hit_by_selection(hit, selection).then_some(idx)
                }
            })
            .collect::<Vec<_>>();
        let selected_total = selected_indices.len();
        let mut cumulative_read_bases_processed = 0u64;
        let mut read_length_counts = vec![0u64; 1];
        let mut cumulative_tested_kmers = 0usize;
        let mut cumulative_matched_kmers = 0usize;
        let mut cumulative_seed_compute_ms = 0.0f64;
        let mut cumulative_align_compute_ms = 0.0f64;
        let mut cumulative_normalize_compute_ms = 0.0f64;
        let mut cumulative_inference_compute_ms = 0.0f64;
        let mut cumulative_progress_emit_ms = 0.0f64;
        let mut reads_processed = 0usize;
        let mut seed_passed = 0usize;
        let mut aligned = 0usize;
        let mut support_aligned_reads = 0usize;
        let mut support_exon_counts = vec![0usize; splicing.unique_exons.len()];
        let mut support_junction_counts = vec![0usize; splicing.junctions.len()];
        let mut reads_with_transition_support = 0usize;
        let mut transition_confirmations = 0usize;
        let mut isoform_support_accumulators =
            Self::build_isoform_support_accumulators(&transcript_exon_models);
        let junction_crossing_seed_bits_indexed = seed_to_transitions.len();
        let mut progress_top_hits = Self::retained_rna_read_preview_heap_from_rows(&report.hits);
        if !on_progress(OperationProgress::RnaReadInterpret(
            RnaReadInterpretProgress {
                seq_id: report.seq_id.clone(),
                reads_processed: 0,
                reads_total: selected_total,
                read_bases_processed: 0,
                mean_read_length_bp: 0.0,
                median_read_length_bp: 0,
                p95_read_length_bp: 0,
                input_bytes_processed: 0,
                input_bytes_total: 0,
                seed_passed: 0,
                aligned: 0,
                tested_kmers: 0,
                matched_kmers: 0,
                seed_compute_ms: 0.0,
                align_compute_ms: 0.0,
                io_read_ms: 0.0,
                fasta_parse_ms: 0.0,
                normalize_compute_ms: 0.0,
                inference_compute_ms: 0.0,
                progress_emit_ms: 0.0,
                update_every_reads: RNA_READ_PROGRESS_UPDATE_EVERY_READS,
                done: false,
                bins: bins.clone(),
                score_density_bins: report.score_density_bins.clone(),
                top_hits_preview: Self::collect_rna_read_top_hit_previews(&progress_top_hits),
                transition_support_rows: transition_support_rows.clone(),
                isoform_support_rows: Self::collect_isoform_support_rows(
                    &isoform_support_accumulators,
                ),
                mapped_exon_support_frequencies: vec![],
                mapped_junction_support_frequencies: vec![],
                mapped_isoform_support_rows: vec![],
                reads_with_transition_support,
                transition_confirmations,
                junction_crossing_seed_bits_indexed,
                origin_class_counts: BTreeMap::new(),
            },
        )) {
            return Err(EngineError {
                code: ErrorCode::Internal,
                message: "RNA-read alignment phase cancelled during progress reporting".to_string(),
            });
        }
        let mut last_progress_emit_at = Instant::now();
        for idx in selected_indices {
            if !should_continue() {
                return Err(EngineError {
                    code: ErrorCode::Internal,
                    message: "RNA-read alignment phase cancelled during processing".to_string(),
                });
            }
            let Some(hit) = report.hits.get_mut(idx) else {
                continue;
            };
            let normalize_started = Instant::now();
            let normalized_sequence = hit
                .sequence
                .as_bytes()
                .iter()
                .map(|base| Self::normalize_nucleotide_base(*base))
                .collect::<Vec<_>>();
            cumulative_normalize_compute_ms += normalize_started.elapsed().as_secs_f64() * 1000.0;
            cumulative_read_bases_processed =
                cumulative_read_bases_processed.saturating_add(normalized_sequence.len() as u64);
            Self::update_read_length_counts(&mut read_length_counts, normalized_sequence.len());
            let seed_started = Instant::now();
            let windows =
                Self::sampled_read_windows(normalized_sequence.len(), &report.seed_filter);
            let mut tested_kmers = 0usize;
            let mut matched_kmers = 0usize;
            let mut matched_seed_bits = HashSet::<u32>::new();
            let mut matched_seed_observations = Vec::<SeedMatchObservation>::new();
            for (start, end) in windows {
                let (tested, matched, matched_bits, matched_observations) =
                    Self::count_seed_hits_in_window_with_histogram(
                        &normalized_sequence[start..end],
                        start,
                        report.seed_filter.kmer_len,
                        report.seed_filter.seed_stride_bp,
                        &seed_index,
                        &histogram_index,
                        &mut bins,
                        &seed_occurrence_counts,
                    );
                tested_kmers = tested_kmers.saturating_add(tested);
                matched_kmers = matched_kmers.saturating_add(matched);
                matched_seed_bits.extend(matched_bits);
                matched_seed_observations.extend(matched_observations);
            }
            cumulative_seed_compute_ms += seed_started.elapsed().as_secs_f64() * 1000.0;
            cumulative_tested_kmers = cumulative_tested_kmers.saturating_add(tested_kmers);
            cumulative_matched_kmers = cumulative_matched_kmers.saturating_add(matched_kmers);
            let weighted_matched_kmers = Self::weighted_seed_support_from_occurrences(
                &matched_seed_bits,
                &seed_occurrence_counts,
            );
            let weighted_seed_hit_fraction = if tested_kmers == 0 {
                0.0
            } else {
                weighted_matched_kmers / tested_kmers as f64
            };
            let inference_started = Instant::now();
            let (read_supported_exons, read_supported_transitions) =
                Self::build_read_seed_support_sets(
                    &matched_seed_bits,
                    &seed_to_exons,
                    &seed_to_transitions,
                );
            let spacing_metrics = Self::compute_seed_chain_spacing_metrics(
                &matched_seed_observations,
                &seed_template_positions,
                &templates,
            );
            let (seed_hit_fraction, perfect_seed_match, _raw_seed_pass) = Self::seed_hit_metrics(
                tested_kmers,
                matched_kmers,
                report.seed_filter.min_seed_hit_fraction,
            );
            let path_inference = Self::infer_read_exon_path(
                &transcript_exon_models,
                &read_supported_exons,
                &read_supported_transitions,
                &spacing_metrics.transcript_id,
            );
            let passed_seed_filter = Self::seed_filter_passes(
                seed_hit_fraction,
                weighted_seed_hit_fraction,
                tested_kmers,
                matched_seed_bits.len(),
                spacing_metrics.support_fraction,
                spacing_metrics.median_transcript_gap,
                spacing_metrics.transcript_gap_count,
                path_inference.confirmed_transitions,
                path_inference.total_transitions,
                &report.seed_filter,
            );
            let target_feature_strand = splicing.strand.as_str();
            let origin_classification = Self::classify_rna_read_origin(
                seed_hit_fraction,
                weighted_seed_hit_fraction,
                passed_seed_filter,
                &spacing_metrics,
                &path_inference,
                target_feature_strand,
                &report.seed_filter,
            );
            let origin_candidates = Self::build_rna_read_origin_candidates(
                &path_inference,
                &spacing_metrics,
                &transcript_models_by_id,
            );
            if passed_seed_filter && !read_supported_transitions.is_empty() {
                reads_with_transition_support = reads_with_transition_support.saturating_add(1);
                transition_confirmations =
                    transition_confirmations.saturating_add(read_supported_transitions.len());
                for transition in &read_supported_transitions {
                    let Some(row_idx) = transition_row_index.get(transition).copied() else {
                        continue;
                    };
                    if let Some(row) = transition_support_rows.get_mut(row_idx) {
                        row.support_read_count = row.support_read_count.saturating_add(1);
                    }
                }
            }
            if !path_inference.transcript_id.is_empty() {
                Self::update_isoform_support_accumulator(
                    &mut isoform_support_accumulators,
                    &path_inference.transcript_id,
                    &path_inference.strand,
                    passed_seed_filter,
                    spacing_metrics.median_transcript_gap,
                    path_inference.confirmed_transitions,
                    path_inference.total_transitions,
                    seed_hit_fraction,
                    weighted_seed_hit_fraction,
                    &path_inference.strand_diagnostics.chain_preferred_strand,
                    path_inference.strand_diagnostics.competing_opposite_strand,
                    path_inference.strand_diagnostics.ambiguous_near_tie,
                    &read_supported_transitions,
                    &transcript_models_by_id,
                );
            }
            cumulative_inference_compute_ms += inference_started.elapsed().as_secs_f64() * 1000.0;
            let align_started = Instant::now();
            let (best_mapping, secondary_mappings) = if passed_seed_filter {
                Self::align_read_to_templates(
                    &normalized_sequence,
                    &templates,
                    &align_config,
                    report.seed_filter.kmer_len,
                )
            } else {
                (None, vec![])
            };
            cumulative_align_compute_ms += align_started.elapsed().as_secs_f64() * 1000.0;
            if passed_seed_filter {
                seed_passed = seed_passed.saturating_add(1);
            }
            if let Some(mapping) = &best_mapping {
                aligned = aligned.saturating_add(1);
                support_aligned_reads = support_aligned_reads.saturating_add(1);
                Self::accumulate_support_counts_for_mapping(
                    mapping,
                    &splicing,
                    &mut support_exon_counts,
                    &mut support_junction_counts,
                );
            }
            hit.read_length_bp = normalized_sequence.len();
            hit.tested_kmers = tested_kmers;
            hit.matched_kmers = matched_kmers;
            hit.seed_hit_fraction = seed_hit_fraction;
            hit.weighted_seed_hit_fraction = weighted_seed_hit_fraction;
            hit.weighted_matched_kmers = weighted_matched_kmers;
            hit.seed_chain_transcript_id = spacing_metrics.transcript_id;
            hit.seed_chain_support_kmers = spacing_metrics.support_kmers;
            hit.seed_chain_support_fraction = spacing_metrics.support_fraction;
            hit.seed_median_transcript_gap = spacing_metrics.median_transcript_gap;
            hit.seed_transcript_gap_count = spacing_metrics.transcript_gap_count;
            hit.exon_path_transcript_id = path_inference.transcript_id;
            hit.exon_path = path_inference.path;
            hit.exon_transitions_confirmed = path_inference.confirmed_transitions;
            hit.exon_transitions_total = path_inference.total_transitions;
            hit.strand_diagnostics = path_inference.strand_diagnostics;
            hit.origin_class = origin_classification.origin_class;
            hit.origin_reason = origin_classification.reason;
            hit.origin_confidence = origin_classification.origin_confidence;
            hit.strand_confidence = origin_classification.strand_confidence;
            hit.origin_candidates = origin_candidates;
            hit.perfect_seed_match = perfect_seed_match;
            hit.passed_seed_filter = passed_seed_filter;
            hit.best_mapping = best_mapping;
            hit.secondary_mappings = secondary_mappings;
            let (msa_eligible, msa_reason) =
                Self::compute_msa_eligibility(hit, &report.seed_filter);
            hit.msa_eligible = msa_eligible;
            hit.msa_eligibility_reason = msa_reason;
            reads_processed = reads_processed.saturating_add(1);
            if Self::should_emit_rna_read_progress(reads_processed, last_progress_emit_at.elapsed())
            {
                let (mean_len, median_len, p95_len) = Self::summarize_read_lengths(
                    &read_length_counts,
                    reads_processed,
                    cumulative_read_bases_processed,
                );
                progress_top_hits = Self::retained_rna_read_preview_heap_from_rows(&report.hits);
                let emit_started = Instant::now();
                let mut origin_class_counts = BTreeMap::<String, usize>::new();
                for row in &report.hits {
                    *origin_class_counts
                        .entry(row.origin_class.as_str().to_string())
                        .or_insert(0) += 1;
                }
                let (mapped_exon_support_frequencies, mapped_junction_support_frequencies) =
                    Self::build_rna_read_support_frequencies_from_counts(
                        &splicing,
                        &support_exon_counts,
                        &support_junction_counts,
                        support_aligned_reads,
                    );
                let mapped_isoform_support_rows =
                    Self::collect_mapped_isoform_support_rows(&report.hits);
                if !on_progress(OperationProgress::RnaReadInterpret(
                    RnaReadInterpretProgress {
                        seq_id: report.seq_id.clone(),
                        reads_processed,
                        reads_total: selected_total,
                        read_bases_processed: cumulative_read_bases_processed,
                        mean_read_length_bp: mean_len,
                        median_read_length_bp: median_len,
                        p95_read_length_bp: p95_len,
                        input_bytes_processed: 0,
                        input_bytes_total: 0,
                        seed_passed,
                        aligned,
                        tested_kmers: cumulative_tested_kmers,
                        matched_kmers: cumulative_matched_kmers,
                        seed_compute_ms: cumulative_seed_compute_ms,
                        align_compute_ms: cumulative_align_compute_ms,
                        io_read_ms: 0.0,
                        fasta_parse_ms: 0.0,
                        normalize_compute_ms: cumulative_normalize_compute_ms,
                        inference_compute_ms: cumulative_inference_compute_ms,
                        progress_emit_ms: cumulative_progress_emit_ms,
                        update_every_reads: RNA_READ_PROGRESS_UPDATE_EVERY_READS,
                        done: false,
                        bins: bins.clone(),
                        score_density_bins: report.score_density_bins.clone(),
                        top_hits_preview: Self::collect_rna_read_top_hit_previews(
                            &progress_top_hits,
                        ),
                        transition_support_rows: transition_support_rows.clone(),
                        isoform_support_rows: Self::collect_isoform_support_rows(
                            &isoform_support_accumulators,
                        ),
                        mapped_exon_support_frequencies,
                        mapped_junction_support_frequencies,
                        mapped_isoform_support_rows,
                        reads_with_transition_support,
                        transition_confirmations,
                        junction_crossing_seed_bits_indexed,
                        origin_class_counts,
                    },
                )) {
                    return Err(EngineError {
                        code: ErrorCode::Internal,
                        message: "RNA-read alignment phase cancelled during progress reporting"
                            .to_string(),
                    });
                }
                cumulative_progress_emit_ms += emit_started.elapsed().as_secs_f64() * 1000.0;
                last_progress_emit_at = Instant::now();
            }
            if reads_processed % RNA_READ_COOPERATIVE_YIELD_EVERY_READS == 0 {
                std::thread::yield_now();
            }
        }
        Self::sort_rna_read_hits_by_retention_rank(&mut report.hits);
        let mut origin_class_counts = BTreeMap::<String, usize>::new();
        for row in &report.hits {
            *origin_class_counts
                .entry(row.origin_class.as_str().to_string())
                .or_insert(0) += 1;
        }
        let (exon_support_frequencies, junction_support_frequencies) =
            Self::build_rna_read_support_frequencies_from_counts(
                &splicing,
                &support_exon_counts,
                &support_junction_counts,
                support_aligned_reads,
            );
        report.align_config = align_config;
        report.generated_at_unix_ms = Self::now_unix_ms();
        report.read_count_aligned = aligned;
        report.retained_count_msa_eligible =
            report.hits.iter().filter(|hit| hit.msa_eligible).count();
        report.exon_support_frequencies = exon_support_frequencies;
        report.junction_support_frequencies = junction_support_frequencies;
        report.transition_support_rows = transition_support_rows.clone();
        report.isoform_support_rows =
            Self::collect_isoform_support_rows(&isoform_support_accumulators);
        report.mapped_isoform_support_rows =
            Self::collect_mapped_isoform_support_rows(&report.hits);
        report.origin_class_counts = origin_class_counts.clone();
        report
            .warnings
            .retain(|warning| !warning.contains("alignment is deferred to a later pass"));
        report.warnings.push(format!(
            "alignment phase completed over retained-hit selection '{}'{} (selected={}, aligned={}, msa_eligible={}); retained hits were re-ranked by alignment-aware retention rank",
            selection.as_str(),
            if explicit_record_filter.is_empty() {
                String::new()
            } else {
                format!(" with explicit_record_indices={}", explicit_record_filter.len())
            },
            selected_total,
            report.read_count_aligned,
            report.retained_count_msa_eligible,
        ));
        if selected_total == 0 {
            if explicit_record_filter.is_empty() {
                report.warnings.push(format!(
                    "alignment phase selection '{}' matched no retained hits",
                    selection.as_str()
                ));
            } else {
                report.warnings.push(format!(
                    "alignment phase explicit_record_indices={} matched no retained hits",
                    explicit_record_filter.len()
                ));
            }
        }
        if report.align_config.max_secondary_mappings == 0 {
            report.warnings.push(
                "alignment phase ran with max_secondary_mappings=0; only best mapping is retained"
                    .to_string(),
            );
        }
        progress_top_hits = Self::retained_rna_read_preview_heap_from_rows(&report.hits);
        let (mean_len, median_len, p95_len) = Self::summarize_read_lengths(
            &read_length_counts,
            reads_processed,
            cumulative_read_bases_processed,
        );
        let mapped_isoform_support_rows = Self::collect_mapped_isoform_support_rows(&report.hits);
        if !on_progress(OperationProgress::RnaReadInterpret(
            RnaReadInterpretProgress {
                seq_id: report.seq_id.clone(),
                reads_processed,
                reads_total: selected_total,
                read_bases_processed: cumulative_read_bases_processed,
                mean_read_length_bp: mean_len,
                median_read_length_bp: median_len,
                p95_read_length_bp: p95_len,
                input_bytes_processed: 0,
                input_bytes_total: 0,
                seed_passed,
                aligned,
                tested_kmers: cumulative_tested_kmers,
                matched_kmers: cumulative_matched_kmers,
                seed_compute_ms: cumulative_seed_compute_ms,
                align_compute_ms: cumulative_align_compute_ms,
                io_read_ms: 0.0,
                fasta_parse_ms: 0.0,
                normalize_compute_ms: cumulative_normalize_compute_ms,
                inference_compute_ms: cumulative_inference_compute_ms,
                progress_emit_ms: cumulative_progress_emit_ms,
                update_every_reads: RNA_READ_PROGRESS_UPDATE_EVERY_READS,
                done: true,
                bins,
                score_density_bins: report.score_density_bins.clone(),
                top_hits_preview: Self::collect_rna_read_top_hit_previews(&progress_top_hits),
                transition_support_rows,
                isoform_support_rows: report.isoform_support_rows.clone(),
                mapped_exon_support_frequencies: report.exon_support_frequencies.clone(),
                mapped_junction_support_frequencies: report.junction_support_frequencies.clone(),
                mapped_isoform_support_rows,
                reads_with_transition_support,
                transition_confirmations,
                junction_crossing_seed_bits_indexed,
                origin_class_counts,
            },
        )) {
            return Err(EngineError {
                code: ErrorCode::Internal,
                message: "RNA-read alignment phase cancelled at completion".to_string(),
            });
        }
        Ok(report)
    }

    pub(super) fn compute_rna_read_report_with_options_and_progress_and_cancel(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
        profile: RnaReadInterpretationProfile,
        input_path: &str,
        input_format: RnaReadInputFormat,
        scope: SplicingScopePreset,
        origin_mode: RnaReadOriginMode,
        target_gene_ids: &[String],
        roi_seed_capture_enabled: bool,
        seed_filter: &RnaReadSeedFilterConfig,
        align_config: &RnaReadAlignConfig,
        report_id: Option<&str>,
        options: &RnaReadInterpretOptions,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
        should_continue: &mut dyn FnMut() -> bool,
    ) -> Result<RnaReadInterpretationReport, EngineError> {
        self.interpret_rna_reads_report_with_progress(
            seq_id,
            seed_feature_id,
            profile,
            input_path,
            input_format,
            scope,
            origin_mode,
            target_gene_ids,
            roi_seed_capture_enabled,
            seed_filter,
            align_config,
            report_id,
            options,
            on_progress,
            should_continue,
        )
    }

    pub fn commit_rna_read_report(
        &mut self,
        report: RnaReadInterpretationReport,
    ) -> Result<OpResult, EngineError> {
        let op = Operation::InterpretRnaReads {
            seq_id: report.seq_id.clone(),
            seed_feature_id: report.seed_feature_id,
            profile: report.profile,
            input_path: report.input_path.clone(),
            input_format: report.input_format,
            scope: report.scope,
            origin_mode: report.origin_mode,
            target_gene_ids: report.target_gene_ids.clone(),
            roi_seed_capture_enabled: report.roi_seed_capture_enabled,
            seed_filter: report.seed_filter.clone(),
            align_config: report.align_config.clone(),
            report_id: Some(report.report_id.clone()),
            report_mode: report.report_mode,
            checkpoint_path: report.checkpoint_path.clone(),
            checkpoint_every_reads: report.checkpoint_every_reads,
            resume_from_checkpoint: report.resumed_from_checkpoint,
        };
        let run_id = "interactive".to_string();
        let checkpoint = self.maybe_capture_checkpoint(&op);
        let mut result = OpResult {
            op_id: self.next_op_id(),
            created_seq_ids: vec![],
            changed_seq_ids: vec![],
            warnings: vec![],
            messages: vec![],
            genome_annotation_projection: None,
            sequence_alignment: None,
        };
        self.push_rna_read_report_result_message(report, &mut result)?;
        self.journal.push(OperationRecord {
            run_id,
            op,
            result: result.clone(),
        });
        if let Some(checkpoint) = checkpoint {
            self.push_undo_checkpoint(checkpoint);
        }
        Ok(result)
    }

    pub fn commit_aligned_rna_read_report(
        &mut self,
        report: RnaReadInterpretationReport,
        selection: RnaReadHitSelection,
        align_config_override: Option<RnaReadAlignConfig>,
        selected_record_indices: Vec<usize>,
    ) -> Result<OpResult, EngineError> {
        let op = Operation::AlignRnaReadReport {
            report_id: report.report_id.clone(),
            selection,
            align_config_override,
            selected_record_indices,
        };
        let run_id = "interactive".to_string();
        let checkpoint = self.maybe_capture_checkpoint(&op);
        let mut result = OpResult {
            op_id: self.next_op_id(),
            created_seq_ids: vec![],
            changed_seq_ids: vec![],
            warnings: vec![],
            messages: vec![],
            genome_annotation_projection: None,
            sequence_alignment: None,
        };
        self.push_rna_read_report_result_message(report, &mut result)?;
        self.journal.push(OperationRecord {
            run_id,
            op,
            result: result.clone(),
        });
        if let Some(checkpoint) = checkpoint {
            self.push_undo_checkpoint(checkpoint);
        }
        Ok(result)
    }

    #[allow(dead_code)]
    pub(super) fn interpret_rna_reads_report(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
        profile: RnaReadInterpretationProfile,
        input_path: &str,
        input_format: RnaReadInputFormat,
        scope: SplicingScopePreset,
        origin_mode: RnaReadOriginMode,
        target_gene_ids: &[String],
        roi_seed_capture_enabled: bool,
        seed_filter: &RnaReadSeedFilterConfig,
        align_config: &RnaReadAlignConfig,
        report_id: Option<&str>,
    ) -> Result<RnaReadInterpretationReport, EngineError> {
        let mut on_progress = |_progress: OperationProgress| true;
        let mut keep_running = || true;
        self.interpret_rna_reads_report_with_progress(
            seq_id,
            seed_feature_id,
            profile,
            input_path,
            input_format,
            scope,
            origin_mode,
            target_gene_ids,
            roi_seed_capture_enabled,
            seed_filter,
            align_config,
            report_id,
            &RnaReadInterpretOptions::default(),
            &mut on_progress,
            &mut keep_running,
        )
    }

    pub(super) fn interpret_rna_reads_report_with_progress(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
        profile: RnaReadInterpretationProfile,
        input_path: &str,
        input_format: RnaReadInputFormat,
        scope: SplicingScopePreset,
        origin_mode: RnaReadOriginMode,
        target_gene_ids: &[String],
        roi_seed_capture_enabled: bool,
        seed_filter: &RnaReadSeedFilterConfig,
        align_config: &RnaReadAlignConfig,
        report_id: Option<&str>,
        options: &RnaReadInterpretOptions,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
        should_continue: &mut dyn FnMut() -> bool,
    ) -> Result<RnaReadInterpretationReport, EngineError> {
        let mut normalized_target_gene_ids = target_gene_ids
            .iter()
            .map(|id| id.trim())
            .filter(|id| !id.is_empty())
            .map(|id| id.to_string())
            .collect::<Vec<_>>();
        normalized_target_gene_ids
            .sort_by(|left, right| left.to_ascii_lowercase().cmp(&right.to_ascii_lowercase()));
        normalized_target_gene_ids.dedup_by(|left, right| left.eq_ignore_ascii_case(right));
        let mut origin_mode_warnings = Vec::<String>::new();
        if matches!(origin_mode, RnaReadOriginMode::SingleGene)
            && !normalized_target_gene_ids.is_empty()
        {
            origin_mode_warnings.push(
                "target_gene_ids were provided but origin_mode=single_gene; extra target genes are ignored unless origin_mode=multi_gene_sparse"
                    .to_string(),
            );
        }
        if roi_seed_capture_enabled {
            origin_mode_warnings.push(
                "roi_seed_capture_enabled=true requested; ROI capture layer is planned but not active in phase-1 scoring"
                    .to_string(),
            );
        }
        if !matches!(profile, RnaReadInterpretationProfile::NanoporeCdnaV1) {
            return Err(EngineError {
                code: ErrorCode::Unsupported,
                message: format!(
                    "Profile '{}' is not implemented yet in phase 1",
                    profile.as_str()
                ),
            });
        }
        if !matches!(input_format, RnaReadInputFormat::Fasta) {
            return Err(EngineError {
                code: ErrorCode::Unsupported,
                message: format!(
                    "Input format '{}' is not supported yet in phase 1",
                    input_format.as_str()
                ),
            });
        }
        if input_path.trim().to_ascii_lowercase().ends_with(".sra") {
            return Err(EngineError {
                code: ErrorCode::Unsupported,
                message: "Direct .sra input is not supported; convert externally to FASTA first"
                    .to_string(),
            });
        }
        if seed_filter.kmer_len == 0 || seed_filter.kmer_len > 16 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "kmer_len must be within 1..=16".to_string(),
            });
        }
        if !(0.0..=1.0).contains(&seed_filter.min_seed_hit_fraction) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "min_seed_hit_fraction must be within 0.0..=1.0".to_string(),
            });
        }
        if !(0.0..=1.0).contains(&seed_filter.min_weighted_seed_hit_fraction) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "min_weighted_seed_hit_fraction must be within 0.0..=1.0".to_string(),
            });
        }
        if !seed_filter.max_median_transcript_gap.is_finite()
            || seed_filter.max_median_transcript_gap < 1.0
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "max_median_transcript_gap must be >= 1.0".to_string(),
            });
        }
        if !(0.0..=1.0).contains(&seed_filter.min_chain_consistency_fraction) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "min_chain_consistency_fraction must be within 0.0..=1.0".to_string(),
            });
        }
        if !(0.0..=1.0).contains(&seed_filter.min_transition_support_fraction) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "min_transition_support_fraction must be within 0.0..=1.0".to_string(),
            });
        }

        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let splicing = self.build_splicing_expert_view(seq_id, seed_feature_id, scope)?;
        let target_feature_strand = splicing.strand.clone();
        let mut transcript_lanes = splicing.transcripts.clone();
        if matches!(origin_mode, RnaReadOriginMode::MultiGeneSparse) {
            let existing_feature_ids = transcript_lanes
                .iter()
                .map(|lane| lane.transcript_feature_id)
                .collect::<HashSet<_>>();
            let (extra_lanes, matched_genes, missing_genes) =
                Self::collect_sparse_target_transcript_lanes(
                    dna,
                    seed_feature_id,
                    scope,
                    target_feature_strand.as_str(),
                    &normalized_target_gene_ids,
                    &existing_feature_ids,
                )?;
            if !extra_lanes.is_empty() {
                let added = extra_lanes.len();
                transcript_lanes.extend(extra_lanes);
                transcript_lanes.sort_by(|left, right| {
                    left.transcript_id
                        .cmp(&right.transcript_id)
                        .then(left.transcript_feature_id.cmp(&right.transcript_feature_id))
                });
                origin_mode_warnings.push(format!(
                    "origin_mode=multi_gene_sparse active: added {} transcript lane(s) from target_gene_ids; total indexed transcript lanes={}",
                    added,
                    transcript_lanes.len()
                ));
            } else if normalized_target_gene_ids.is_empty() {
                origin_mode_warnings.push(
                    "origin_mode=multi_gene_sparse active with empty target_gene_ids; using baseline transcript scope only"
                        .to_string(),
                );
            } else {
                origin_mode_warnings.push(
                    "origin_mode=multi_gene_sparse active: no additional transcript lanes were added from target_gene_ids (baseline scope already covered requested genes or none matched)"
                        .to_string(),
                );
            }
            if !matched_genes.is_empty() {
                origin_mode_warnings.push(format!(
                    "multi_gene_sparse matched target genes: {}",
                    matched_genes.join(",")
                ));
            }
            if !missing_genes.is_empty() {
                origin_mode_warnings.push(format!(
                    "multi_gene_sparse target genes not found in local annotation: {}",
                    missing_genes.join(",")
                ));
            }
        }
        let templates = transcript_lanes
            .iter()
            .map(|lane| Self::make_transcript_template(dna, lane, seed_filter.kmer_len))
            .filter(|template| !template.sequence.is_empty())
            .collect::<Vec<_>>();
        if templates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No transcript templates available for splicing scope '{}' on '{}'",
                    scope.as_str(),
                    seq_id
                ),
            });
        }
        let mut seed_index: HashSet<u32> = HashSet::new();
        for template in &templates {
            seed_index.extend(template.kmer_positions.keys().copied());
        }
        if seed_index.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: "No directional k-mer seeds could be generated from transcript templates"
                    .to_string(),
            });
        }
        let seed_catalog_rows =
            Self::collect_rna_seed_hash_catalog_rows(&templates, seed_filter.kmer_len);
        let mut seed_occurrence_counts = HashMap::<u32, usize>::new();
        for row in &seed_catalog_rows {
            *seed_occurrence_counts.entry(row.seed_bits).or_insert(0) += 1;
        }
        let seed_support_exons = Self::collect_seed_support_exon_summaries(&transcript_lanes);
        let (
            seed_to_exons,
            seed_to_transitions,
            transcript_exon_models,
            mut transition_support_rows,
        ) = Self::build_seed_support_indexes(&seed_support_exons, &templates, seed_filter.kmer_len);
        let transcript_models_by_id = transcript_exon_models
            .iter()
            .map(|model| (model.transcript_id.clone(), model.clone()))
            .collect::<HashMap<_, _>>();
        let mut isoform_support_accumulators =
            Self::build_isoform_support_accumulators(&transcript_exon_models);
        let mut isoform_support_rows: Vec<RnaReadIsoformSupportRow>;
        let transition_row_index = transition_support_rows
            .iter()
            .enumerate()
            .map(|(idx, row)| ((row.from_exon_ordinal, row.to_exon_ordinal), idx))
            .collect::<HashMap<_, _>>();
        let junction_crossing_seed_bits_indexed = seed_to_transitions.len();
        let mut bins = Self::build_rna_read_seed_histogram_bins(dna.len());
        let histogram_index =
            Self::build_rna_read_seed_histogram_index(&templates, dna.len(), &bins);
        let seed_template_positions = Self::build_seed_template_position_index(&templates);
        let report_mode = options.report_mode;
        let checkpoint_path =
            Self::normalize_rna_read_checkpoint_path(options.checkpoint_path.as_deref());
        if options.resume_from_checkpoint && checkpoint_path.is_none() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "resume_from_checkpoint=true requires checkpoint_path to be set"
                    .to_string(),
            });
        }
        let checkpoint_every_reads = options.checkpoint_every_reads.max(1);
        let loaded_checkpoint = if options.resume_from_checkpoint {
            let path = checkpoint_path.as_deref().ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: "resume_from_checkpoint=true requires checkpoint_path to be set"
                    .to_string(),
            })?;
            let checkpoint = Self::read_rna_read_interpret_checkpoint(path)?;
            if checkpoint.seq_id != seq_id
                || checkpoint.seed_feature_id != seed_feature_id
                || checkpoint.profile != profile
                || checkpoint.input_path != input_path
                || checkpoint.input_format != input_format
                || checkpoint.scope != scope
                || checkpoint.origin_mode != origin_mode
                || checkpoint.target_gene_ids != normalized_target_gene_ids
                || checkpoint.roi_seed_capture_enabled != roi_seed_capture_enabled
                || checkpoint.seed_filter != *seed_filter
                || checkpoint.align_config != *align_config
                || checkpoint.report_mode != report_mode
            {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "RNA-read checkpoint '{}' does not match requested interpret parameters",
                        path
                    ),
                });
            }
            Some(checkpoint)
        } else {
            None
        };
        let report_id = match report_id {
            Some(raw) => Self::normalize_rna_read_report_id(raw)?,
            None => {
                if let Some(checkpoint) = loaded_checkpoint.as_ref() {
                    Self::normalize_rna_read_report_id(&checkpoint.report_id)?
                } else {
                    Self::normalize_rna_read_report_id(&format!(
                        "rna_reads_{}_{}",
                        seq_id,
                        Self::now_unix_ms()
                    ))?
                }
            }
        };
        if let Some(checkpoint) = loaded_checkpoint.as_ref() {
            let checkpoint_report_id = Self::normalize_rna_read_report_id(&checkpoint.report_id)?;
            if checkpoint_report_id != report_id {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Checkpoint report_id '{}' does not match requested report_id '{}'",
                        checkpoint_report_id, report_id
                    ),
                });
            }
        }
        let mut score_density_bins = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.score_density_bins.clone())
            .unwrap_or_else(|| vec![0u64; RNA_READ_SCORE_DENSITY_BIN_COUNT]);
        if score_density_bins.len() != RNA_READ_SCORE_DENSITY_BIN_COUNT {
            score_density_bins.resize(RNA_READ_SCORE_DENSITY_BIN_COUNT, 0);
        }
        if let Some(checkpoint) = loaded_checkpoint.as_ref() {
            if !checkpoint.bins.is_empty() {
                bins = checkpoint.bins.clone();
            }
            if !checkpoint.transition_support_rows.is_empty()
                && checkpoint.transition_support_rows.len() == transition_support_rows.len()
            {
                transition_support_rows = checkpoint.transition_support_rows.clone();
            }
            if !checkpoint.isoform_support_accumulators.is_empty() {
                isoform_support_accumulators = checkpoint.isoform_support_accumulators.clone();
            }
        }
        isoform_support_rows = Self::collect_isoform_support_rows(&isoform_support_accumulators);
        let mut input_bytes_processed = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.input_bytes_processed)
            .unwrap_or(0);
        let mut input_bytes_total = std::fs::metadata(input_path)
            .map(|meta| meta.len())
            .unwrap_or(0);
        if let Some(checkpoint) = loaded_checkpoint.as_ref() {
            input_bytes_total = input_bytes_total.max(checkpoint.input_bytes_total);
        }
        let mut retained_hits = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| Self::retained_rna_read_hit_heap_from_rows(&checkpoint.retained_hits))
            .unwrap_or_default();
        let mut progress_top_hits = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| {
                Self::retained_rna_read_preview_heap_from_rows(&checkpoint.retained_hits)
            })
            .unwrap_or_default();
        let mut seed_passed = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.read_count_seed_passed)
            .unwrap_or(0);
        let mut aligned = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.read_count_aligned)
            .unwrap_or(0);
        let mut reads_processed = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.reads_processed)
            .unwrap_or(0);
        let mut cumulative_tested_kmers = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.cumulative_tested_kmers)
            .unwrap_or(0);
        let mut cumulative_matched_kmers = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.cumulative_matched_kmers)
            .unwrap_or(0);
        let mut cumulative_seed_compute_ms = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.cumulative_seed_compute_ms)
            .unwrap_or(0.0);
        let mut cumulative_align_compute_ms = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.cumulative_align_compute_ms)
            .unwrap_or(0.0);
        let resume_io_read_ms_base = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.cumulative_io_read_ms)
            .unwrap_or(0.0);
        let resume_fasta_parse_ms_base = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.cumulative_fasta_parse_ms)
            .unwrap_or(0.0);
        let mut cumulative_io_read_ms = resume_io_read_ms_base;
        let mut cumulative_fasta_parse_ms = resume_fasta_parse_ms_base;
        let mut cumulative_normalize_compute_ms = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.cumulative_normalize_compute_ms)
            .unwrap_or(0.0);
        let mut cumulative_inference_compute_ms = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.cumulative_inference_compute_ms)
            .unwrap_or(0.0);
        let mut cumulative_progress_emit_ms = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.cumulative_progress_emit_ms)
            .unwrap_or(0.0);
        let mut cumulative_read_bases_processed = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.cumulative_read_bases_processed)
            .unwrap_or(0);
        let mut read_length_counts = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.read_length_counts.clone())
            .unwrap_or_else(|| vec![0u64; 1]);
        if read_length_counts.is_empty() {
            read_length_counts.push(0);
        }
        let mut support_aligned_reads = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.support_aligned_reads)
            .unwrap_or(0);
        let mut support_exon_counts = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.support_exon_counts.clone())
            .unwrap_or_else(|| vec![0usize; splicing.unique_exons.len()]);
        if support_exon_counts.len() != splicing.unique_exons.len() {
            support_exon_counts.resize(splicing.unique_exons.len(), 0);
        }
        let mut support_junction_counts = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.support_junction_counts.clone())
            .unwrap_or_else(|| vec![0usize; splicing.junctions.len()]);
        if support_junction_counts.len() != splicing.junctions.len() {
            support_junction_counts.resize(splicing.junctions.len(), 0);
        }
        let mut reads_with_transition_support = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.reads_with_transition_support)
            .unwrap_or(0);
        let mut transition_confirmations = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.transition_confirmations)
            .unwrap_or(0);
        let mut origin_class_counts = loaded_checkpoint
            .as_ref()
            .map(|checkpoint| checkpoint.origin_class_counts.clone())
            .unwrap_or_default();
        let resume_skip_records = reads_processed;
        let resumed_from_checkpoint = loaded_checkpoint.is_some();
        if resumed_from_checkpoint {
            let checkpoint_origin = checkpoint_path.as_deref().unwrap_or("<unknown>");
            origin_mode_warnings.push(format!(
                "resumed from checkpoint '{}' at record {}",
                checkpoint_origin, resume_skip_records
            ));
        }
        if checkpoint_path.is_some() {
            origin_mode_warnings.push(format!(
                "checkpoint snapshots enabled every {} reads",
                checkpoint_every_reads
            ));
        }
        if !on_progress(OperationProgress::RnaReadInterpret(
            RnaReadInterpretProgress {
                seq_id: seq_id.to_string(),
                reads_processed,
                reads_total: 0,
                read_bases_processed: cumulative_read_bases_processed,
                mean_read_length_bp: 0.0,
                median_read_length_bp: 0,
                p95_read_length_bp: 0,
                input_bytes_processed,
                input_bytes_total,
                seed_passed,
                aligned,
                tested_kmers: cumulative_tested_kmers,
                matched_kmers: cumulative_matched_kmers,
                seed_compute_ms: cumulative_seed_compute_ms,
                align_compute_ms: cumulative_align_compute_ms,
                io_read_ms: cumulative_io_read_ms,
                fasta_parse_ms: cumulative_fasta_parse_ms,
                normalize_compute_ms: cumulative_normalize_compute_ms,
                inference_compute_ms: cumulative_inference_compute_ms,
                progress_emit_ms: cumulative_progress_emit_ms,
                update_every_reads: RNA_READ_PROGRESS_UPDATE_EVERY_READS,
                done: false,
                bins: bins.clone(),
                score_density_bins: score_density_bins.clone(),
                top_hits_preview: Self::collect_rna_read_top_hit_previews(&progress_top_hits),
                transition_support_rows: transition_support_rows.clone(),
                isoform_support_rows: isoform_support_rows.clone(),
                mapped_exon_support_frequencies: vec![],
                mapped_junction_support_frequencies: vec![],
                mapped_isoform_support_rows: vec![],
                reads_with_transition_support,
                transition_confirmations,
                junction_crossing_seed_bits_indexed,
                origin_class_counts: origin_class_counts.clone(),
            },
        )) {
            return Err(EngineError {
                code: ErrorCode::Internal,
                message: "RNA-read interpretation cancelled before FASTA scan started".to_string(),
            });
        }
        let mut last_progress_emit_at = Instant::now();
        if !should_continue() {
            return Err(EngineError {
                code: ErrorCode::Internal,
                message: "RNA-read interpretation cancelled before processing started".to_string(),
            });
        }
        let alignment_enabled = !matches!(profile, RnaReadInterpretationProfile::NanoporeCdnaV1)
            && align_config.max_secondary_mappings > 0;
        let final_visit_progress =
            Self::visit_fasta_records_with_offsets(input_path, &mut |record, visit_progress| {
                if !should_continue() {
                    return Err(EngineError {
                        code: ErrorCode::Internal,
                        message: "RNA-read interpretation cancelled during FASTA scan".to_string(),
                    });
                }
                input_bytes_processed =
                    input_bytes_processed.max(visit_progress.input_bytes_processed);
                input_bytes_total = input_bytes_total.max(visit_progress.input_bytes_total);
                cumulative_io_read_ms = resume_io_read_ms_base + visit_progress.io_read_ms;
                cumulative_fasta_parse_ms =
                    resume_fasta_parse_ms_base + visit_progress.record_parse_ms;
                if record.record_index < resume_skip_records {
                    return Ok(());
                }
                let normalize_started = Instant::now();
                let (normalized_sequence, reverse_complement_applied) =
                    Self::normalize_rna_read_sequence_for_scoring(&record.sequence, seed_filter);
                cumulative_normalize_compute_ms +=
                    normalize_started.elapsed().as_secs_f64() * 1000.0;
                cumulative_read_bases_processed = cumulative_read_bases_processed
                    .saturating_add(normalized_sequence.len() as u64);
                Self::update_read_length_counts(&mut read_length_counts, normalized_sequence.len());
                let windows = Self::sampled_read_windows(normalized_sequence.len(), seed_filter);
                let mut tested_kmers = 0usize;
                let mut matched_kmers = 0usize;
                let mut matched_seed_bits = HashSet::<u32>::new();
                let mut matched_seed_observations = Vec::<SeedMatchObservation>::new();
                let seed_started = Instant::now();
                for (start, end) in windows {
                    if !should_continue() {
                        return Err(EngineError {
                            code: ErrorCode::Internal,
                            message: "RNA-read interpretation cancelled during seed scanning"
                                .to_string(),
                        });
                    }
                    let (tested, matched, matched_bits, matched_observations) =
                        Self::count_seed_hits_in_window_with_histogram(
                            &normalized_sequence[start..end],
                            start,
                            seed_filter.kmer_len,
                            seed_filter.seed_stride_bp,
                            &seed_index,
                            &histogram_index,
                            &mut bins,
                            &seed_occurrence_counts,
                        );
                    tested_kmers = tested_kmers.saturating_add(tested);
                    matched_kmers = matched_kmers.saturating_add(matched);
                    matched_seed_bits.extend(matched_bits);
                    matched_seed_observations.extend(matched_observations);
                }
                cumulative_seed_compute_ms += seed_started.elapsed().as_secs_f64() * 1000.0;
                cumulative_tested_kmers = cumulative_tested_kmers.saturating_add(tested_kmers);
                cumulative_matched_kmers = cumulative_matched_kmers.saturating_add(matched_kmers);
                let weighted_matched_kmers = Self::weighted_seed_support_from_occurrences(
                    &matched_seed_bits,
                    &seed_occurrence_counts,
                );
                let weighted_seed_hit_fraction = if tested_kmers == 0 {
                    0.0
                } else {
                    weighted_matched_kmers / tested_kmers as f64
                };
                let inference_started = Instant::now();
                let should_parallelize_infer = matched_seed_bits.len()
                    >= RNA_READ_INFER_PARALLEL_MIN_MATCHED_BITS
                    || matched_seed_observations.len()
                        >= RNA_READ_INFER_PARALLEL_MIN_MATCHED_OBSERVATIONS;
                let ((read_supported_exons, read_supported_transitions), spacing_metrics) =
                    if should_parallelize_infer {
                        join(
                            || {
                                Self::build_read_seed_support_sets(
                                    &matched_seed_bits,
                                    &seed_to_exons,
                                    &seed_to_transitions,
                                )
                            },
                            || {
                                Self::compute_seed_chain_spacing_metrics(
                                    &matched_seed_observations,
                                    &seed_template_positions,
                                    &templates,
                                )
                            },
                        )
                    } else {
                        (
                            Self::build_read_seed_support_sets(
                                &matched_seed_bits,
                                &seed_to_exons,
                                &seed_to_transitions,
                            ),
                            Self::compute_seed_chain_spacing_metrics(
                                &matched_seed_observations,
                                &seed_template_positions,
                                &templates,
                            ),
                        )
                    };
                let (seed_hit_fraction, perfect_seed_match, _raw_passed_seed_filter) =
                    Self::seed_hit_metrics(
                        tested_kmers,
                        matched_kmers,
                        seed_filter.min_seed_hit_fraction,
                    );
                let path_inference = Self::infer_read_exon_path(
                    &transcript_exon_models,
                    &read_supported_exons,
                    &read_supported_transitions,
                    &spacing_metrics.transcript_id,
                );
                let passed_seed_filter = Self::seed_filter_passes(
                    seed_hit_fraction,
                    weighted_seed_hit_fraction,
                    tested_kmers,
                    matched_seed_bits.len(),
                    spacing_metrics.support_fraction,
                    spacing_metrics.median_transcript_gap,
                    spacing_metrics.transcript_gap_count,
                    path_inference.confirmed_transitions,
                    path_inference.total_transitions,
                    seed_filter,
                );
                let origin_classification = Self::classify_rna_read_origin(
                    seed_hit_fraction,
                    weighted_seed_hit_fraction,
                    passed_seed_filter,
                    &spacing_metrics,
                    &path_inference,
                    target_feature_strand.as_str(),
                    seed_filter,
                );
                let origin_candidates = Self::build_rna_read_origin_candidates(
                    &path_inference,
                    &spacing_metrics,
                    &transcript_models_by_id,
                );
                *origin_class_counts
                    .entry(origin_classification.origin_class.as_str().to_string())
                    .or_insert(0) += 1;
                if passed_seed_filter && !read_supported_transitions.is_empty() {
                    reads_with_transition_support = reads_with_transition_support.saturating_add(1);
                    transition_confirmations =
                        transition_confirmations.saturating_add(read_supported_transitions.len());
                    for transition in &read_supported_transitions {
                        let Some(row_idx) = transition_row_index.get(transition).copied() else {
                            continue;
                        };
                        if let Some(row) = transition_support_rows.get_mut(row_idx) {
                            row.support_read_count = row.support_read_count.saturating_add(1);
                        }
                    }
                }
                if !path_inference.transcript_id.is_empty() {
                    Self::update_isoform_support_accumulator(
                        &mut isoform_support_accumulators,
                        &path_inference.transcript_id,
                        &path_inference.strand,
                        passed_seed_filter,
                        spacing_metrics.median_transcript_gap,
                        path_inference.confirmed_transitions,
                        path_inference.total_transitions,
                        seed_hit_fraction,
                        weighted_seed_hit_fraction,
                        &path_inference.strand_diagnostics.chain_preferred_strand,
                        path_inference.strand_diagnostics.competing_opposite_strand,
                        path_inference.strand_diagnostics.ambiguous_near_tie,
                        &read_supported_transitions,
                        &transcript_models_by_id,
                    );
                }
                cumulative_inference_compute_ms +=
                    inference_started.elapsed().as_secs_f64() * 1000.0;
                let density_idx = Self::score_density_bin_index(seed_hit_fraction);
                if let Some(bin) = score_density_bins.get_mut(density_idx) {
                    *bin = bin.saturating_add(1);
                }
                let align_started = Instant::now();
                let (best_mapping, secondary_mappings) = if passed_seed_filter && alignment_enabled
                {
                    if !should_continue() {
                        return Err(EngineError {
                            code: ErrorCode::Internal,
                            message: "RNA-read interpretation cancelled before alignment"
                                .to_string(),
                        });
                    }
                    Self::align_read_to_templates(
                        &normalized_sequence,
                        &templates,
                        align_config,
                        seed_filter.kmer_len,
                    )
                } else {
                    (None, vec![])
                };
                if passed_seed_filter && alignment_enabled {
                    cumulative_align_compute_ms += align_started.elapsed().as_secs_f64() * 1000.0;
                }
                if passed_seed_filter {
                    seed_passed += 1;
                }
                if let Some(mapping) = &best_mapping {
                    aligned += 1;
                    support_aligned_reads = support_aligned_reads.saturating_add(1);
                    Self::accumulate_support_counts_for_mapping(
                        mapping,
                        &splicing,
                        &mut support_exon_counts,
                        &mut support_junction_counts,
                    );
                }
                let hit = RnaReadInterpretationHit {
                    record_index: record.record_index,
                    source_byte_offset: record.source_byte_offset,
                    header_id: record.header_id,
                    sequence: String::from_utf8_lossy(&normalized_sequence).to_string(),
                    read_length_bp: normalized_sequence.len(),
                    tested_kmers,
                    matched_kmers,
                    seed_hit_fraction,
                    weighted_seed_hit_fraction,
                    weighted_matched_kmers,
                    seed_chain_transcript_id: spacing_metrics.transcript_id,
                    seed_chain_support_kmers: spacing_metrics.support_kmers,
                    seed_chain_support_fraction: spacing_metrics.support_fraction,
                    seed_median_transcript_gap: spacing_metrics.median_transcript_gap,
                    seed_transcript_gap_count: spacing_metrics.transcript_gap_count,
                    exon_path_transcript_id: path_inference.transcript_id,
                    exon_path: path_inference.path,
                    exon_transitions_confirmed: path_inference.confirmed_transitions,
                    exon_transitions_total: path_inference.total_transitions,
                    reverse_complement_applied,
                    strand_diagnostics: path_inference.strand_diagnostics,
                    origin_class: origin_classification.origin_class,
                    origin_reason: origin_classification.reason,
                    origin_confidence: origin_classification.origin_confidence,
                    strand_confidence: origin_classification.strand_confidence,
                    origin_candidates,
                    perfect_seed_match,
                    passed_seed_filter,
                    msa_eligible: false,
                    msa_eligibility_reason: String::new(),
                    best_mapping,
                    secondary_mappings,
                };
                let mut hit = hit;
                let (msa_eligible, msa_reason) = Self::compute_msa_eligibility(&hit, seed_filter);
                hit.msa_eligible = msa_eligible;
                hit.msa_eligibility_reason = msa_reason;
                Self::retain_top_rna_read_preview_hit(&mut progress_top_hits, &hit);
                Self::retain_top_rna_read_hit(&mut retained_hits, hit);
                reads_processed = reads_processed.saturating_add(1);
                let should_emit = Self::should_emit_rna_read_progress(
                    reads_processed,
                    last_progress_emit_at.elapsed(),
                );
                if should_emit {
                    isoform_support_rows =
                        Self::collect_isoform_support_rows(&isoform_support_accumulators);
                    let (mean_len, median_len, p95_len) = Self::summarize_read_lengths(
                        &read_length_counts,
                        reads_processed,
                        cumulative_read_bases_processed,
                    );
                    let emit_started = Instant::now();
                    let progress_event =
                        OperationProgress::RnaReadInterpret(RnaReadInterpretProgress {
                            seq_id: seq_id.to_string(),
                            reads_processed,
                            reads_total: 0,
                            read_bases_processed: cumulative_read_bases_processed,
                            mean_read_length_bp: mean_len,
                            median_read_length_bp: median_len,
                            p95_read_length_bp: p95_len,
                            input_bytes_processed,
                            input_bytes_total,
                            seed_passed,
                            aligned,
                            tested_kmers: cumulative_tested_kmers,
                            matched_kmers: cumulative_matched_kmers,
                            seed_compute_ms: cumulative_seed_compute_ms,
                            align_compute_ms: cumulative_align_compute_ms,
                            io_read_ms: cumulative_io_read_ms,
                            fasta_parse_ms: cumulative_fasta_parse_ms,
                            normalize_compute_ms: cumulative_normalize_compute_ms,
                            inference_compute_ms: cumulative_inference_compute_ms,
                            progress_emit_ms: cumulative_progress_emit_ms,
                            update_every_reads: RNA_READ_PROGRESS_UPDATE_EVERY_READS,
                            done: false,
                            bins: bins.clone(),
                            score_density_bins: score_density_bins.clone(),
                            top_hits_preview: Self::collect_rna_read_top_hit_previews(
                                &progress_top_hits,
                            ),
                            transition_support_rows: transition_support_rows.clone(),
                            isoform_support_rows: isoform_support_rows.clone(),
                            mapped_exon_support_frequencies: vec![],
                            mapped_junction_support_frequencies: vec![],
                            mapped_isoform_support_rows: vec![],
                            reads_with_transition_support,
                            transition_confirmations,
                            junction_crossing_seed_bits_indexed,
                            origin_class_counts: origin_class_counts.clone(),
                        });
                    if !on_progress(progress_event) {
                        return Err(EngineError {
                            code: ErrorCode::Internal,
                            message: "RNA-read interpretation cancelled during progress reporting"
                                .to_string(),
                        });
                    }
                    cumulative_progress_emit_ms += emit_started.elapsed().as_secs_f64() * 1000.0;
                    last_progress_emit_at = Instant::now();
                }
                // Cooperative scheduler yield to keep GUI/event-loop responsiveness while
                // high-throughput mapping saturates available CPU.
                if reads_processed % RNA_READ_COOPERATIVE_YIELD_EVERY_READS == 0 {
                    std::thread::yield_now();
                }
                if let Some(path) = checkpoint_path.as_deref() {
                    if reads_processed > 0 && reads_processed % checkpoint_every_reads == 0 {
                        let mut retained_rows = retained_hits.iter().cloned().collect::<Vec<_>>();
                        retained_rows.sort_by(|left, right| right.rank.cmp(&left.rank));
                        let checkpoint = RnaReadInterpretCheckpoint {
                            schema: RNA_READ_CHECKPOINT_SCHEMA.to_string(),
                            created_at_unix_ms: Self::now_unix_ms(),
                            report_id: report_id.clone(),
                            report_mode,
                            seq_id: seq_id.to_string(),
                            seed_feature_id,
                            profile,
                            input_path: input_path.to_string(),
                            input_format,
                            scope,
                            origin_mode,
                            target_gene_ids: normalized_target_gene_ids.clone(),
                            roi_seed_capture_enabled,
                            seed_filter: seed_filter.clone(),
                            align_config: align_config.clone(),
                            checkpoint_path: checkpoint_path.clone(),
                            checkpoint_every_reads,
                            reads_processed,
                            read_count_seed_passed: seed_passed,
                            read_count_aligned: aligned,
                            input_bytes_processed,
                            input_bytes_total,
                            cumulative_tested_kmers,
                            cumulative_matched_kmers,
                            cumulative_seed_compute_ms,
                            cumulative_align_compute_ms,
                            cumulative_io_read_ms,
                            cumulative_fasta_parse_ms,
                            cumulative_normalize_compute_ms,
                            cumulative_inference_compute_ms,
                            cumulative_progress_emit_ms,
                            cumulative_read_bases_processed,
                            read_length_counts: read_length_counts.clone(),
                            support_aligned_reads,
                            support_exon_counts: support_exon_counts.clone(),
                            support_junction_counts: support_junction_counts.clone(),
                            reads_with_transition_support,
                            transition_confirmations,
                            transition_support_rows: transition_support_rows.clone(),
                            isoform_support_accumulators: isoform_support_accumulators.clone(),
                            origin_class_counts: origin_class_counts.clone(),
                            bins: bins.clone(),
                            score_density_bins: score_density_bins.clone(),
                            retained_hits: retained_rows.into_iter().map(|row| row.hit).collect(),
                        };
                        Self::write_rna_read_interpret_checkpoint(path, &checkpoint)?;
                    }
                }
                Ok(())
            })?;
        input_bytes_processed = final_visit_progress.input_bytes_processed;
        input_bytes_total = final_visit_progress.input_bytes_total;
        let reads_total = final_visit_progress.records_processed;
        let mut ranked_hits = retained_hits.into_vec();
        ranked_hits.sort_by(|left, right| right.rank.cmp(&left.rank));
        let retained_hits = ranked_hits
            .into_iter()
            .map(|row| row.hit)
            .collect::<Vec<_>>();
        let retained_hit_count = retained_hits.len();
        let hits = Self::apply_rna_read_report_mode_to_hits(report_mode, retained_hits.clone());
        let report_hit_count = hits.len();
        let retained_count_msa_eligible = hits.iter().filter(|hit| hit.msa_eligible).count();
        isoform_support_rows = Self::collect_isoform_support_rows(&isoform_support_accumulators);
        let (mapped_exon_support_frequencies, mapped_junction_support_frequencies) =
            Self::build_rna_read_support_frequencies_from_counts(
                &splicing,
                &support_exon_counts,
                &support_junction_counts,
                support_aligned_reads,
            );
        let mapped_isoform_support_rows = Self::collect_mapped_isoform_support_rows(&hits);
        let (mean_len, median_len, p95_len) = Self::summarize_read_lengths(
            &read_length_counts,
            reads_total,
            cumulative_read_bases_processed,
        );
        let final_emit_started = Instant::now();
        let final_progress = OperationProgress::RnaReadInterpret(RnaReadInterpretProgress {
            seq_id: seq_id.to_string(),
            reads_processed: reads_total,
            reads_total,
            read_bases_processed: cumulative_read_bases_processed,
            mean_read_length_bp: mean_len,
            median_read_length_bp: median_len,
            p95_read_length_bp: p95_len,
            input_bytes_processed,
            input_bytes_total,
            seed_passed,
            aligned,
            tested_kmers: cumulative_tested_kmers,
            matched_kmers: cumulative_matched_kmers,
            seed_compute_ms: cumulative_seed_compute_ms,
            align_compute_ms: cumulative_align_compute_ms,
            io_read_ms: cumulative_io_read_ms,
            fasta_parse_ms: cumulative_fasta_parse_ms,
            normalize_compute_ms: cumulative_normalize_compute_ms,
            inference_compute_ms: cumulative_inference_compute_ms,
            progress_emit_ms: cumulative_progress_emit_ms,
            update_every_reads: RNA_READ_PROGRESS_UPDATE_EVERY_READS,
            done: true,
            bins: bins.clone(),
            score_density_bins: score_density_bins.clone(),
            top_hits_preview: Self::collect_rna_read_top_hit_previews(&progress_top_hits),
            transition_support_rows: transition_support_rows.clone(),
            isoform_support_rows: isoform_support_rows.clone(),
            mapped_exon_support_frequencies: mapped_exon_support_frequencies.clone(),
            mapped_junction_support_frequencies: mapped_junction_support_frequencies.clone(),
            mapped_isoform_support_rows: mapped_isoform_support_rows.clone(),
            reads_with_transition_support,
            transition_confirmations,
            junction_crossing_seed_bits_indexed,
            origin_class_counts: origin_class_counts.clone(),
        });
        if !on_progress(final_progress) {
            return Err(EngineError {
                code: ErrorCode::Internal,
                message: "RNA-read interpretation cancelled at completion".to_string(),
            });
        }
        let _final_emit_ms = final_emit_started.elapsed().as_secs_f64() * 1000.0;
        if let Some(path) = checkpoint_path.as_deref() {
            let checkpoint = RnaReadInterpretCheckpoint {
                schema: RNA_READ_CHECKPOINT_SCHEMA.to_string(),
                created_at_unix_ms: Self::now_unix_ms(),
                report_id: report_id.clone(),
                report_mode,
                seq_id: seq_id.to_string(),
                seed_feature_id,
                profile,
                input_path: input_path.to_string(),
                input_format,
                scope,
                origin_mode,
                target_gene_ids: normalized_target_gene_ids.clone(),
                roi_seed_capture_enabled,
                seed_filter: seed_filter.clone(),
                align_config: align_config.clone(),
                checkpoint_path: checkpoint_path.clone(),
                checkpoint_every_reads,
                reads_processed: reads_total,
                read_count_seed_passed: seed_passed,
                read_count_aligned: aligned,
                input_bytes_processed,
                input_bytes_total,
                cumulative_tested_kmers,
                cumulative_matched_kmers,
                cumulative_seed_compute_ms,
                cumulative_align_compute_ms,
                cumulative_io_read_ms,
                cumulative_fasta_parse_ms,
                cumulative_normalize_compute_ms,
                cumulative_inference_compute_ms,
                cumulative_progress_emit_ms,
                cumulative_read_bases_processed,
                read_length_counts: read_length_counts.clone(),
                support_aligned_reads,
                support_exon_counts: support_exon_counts.clone(),
                support_junction_counts: support_junction_counts.clone(),
                reads_with_transition_support,
                transition_confirmations,
                transition_support_rows: transition_support_rows.clone(),
                isoform_support_accumulators: isoform_support_accumulators.clone(),
                origin_class_counts: origin_class_counts.clone(),
                bins: bins.clone(),
                score_density_bins: score_density_bins.clone(),
                retained_hits: retained_hits.clone(),
            };
            Self::write_rna_read_interpret_checkpoint(path, &checkpoint)?;
        }

        let mut warnings = origin_mode_warnings;
        if !alignment_enabled {
            warnings.push(
                "phase-1 profile runs seed filtering only; alignment is deferred to a later pass"
                    .to_string(),
            );
        }
        warnings.push(
            "phase-1 seed filtering hashes full read span for every read; short/long window knobs are compatibility fields and currently have no effect"
                .to_string(),
        );
        warnings.push(format!(
            "seed-pass gate: raw>={:.3} and weighted>={:.3} and unique_matched_kmers>={} and chain_consistency>={:.2} and median_transcript_gap<={:.2} and confirmed_transitions>={} and transition_fraction>={:.2}",
            seed_filter.min_seed_hit_fraction,
            seed_filter.min_weighted_seed_hit_fraction,
            seed_filter.min_unique_matched_kmers,
            seed_filter.min_chain_consistency_fraction,
            seed_filter.max_median_transcript_gap,
            seed_filter.min_confirmed_exon_transitions,
            seed_filter.min_transition_support_fraction,
        ));
        if reads_total > retained_hit_count {
            warnings.push(format!(
                "retained_top_hits={} out of total_reads={} (scored by retention rank; alignment quality contributes when mappings are available)",
                retained_hit_count, reads_total
            ));
        }
        if report_mode == RnaReadReportMode::SeedPassedOnly {
            warnings.push(format!(
                "report_mode=seed_passed_only kept {} of {} retained hits",
                report_hit_count, retained_hit_count
            ));
        }

        Ok(RnaReadInterpretationReport {
            schema: RNA_READ_REPORT_SCHEMA.to_string(),
            report_id,
            report_mode,
            seq_id: seq_id.to_string(),
            seed_feature_id,
            generated_at_unix_ms: Self::now_unix_ms(),
            profile,
            input_path: input_path.to_string(),
            input_format,
            scope,
            origin_mode,
            target_gene_ids: normalized_target_gene_ids,
            roi_seed_capture_enabled,
            checkpoint_path,
            checkpoint_every_reads,
            resumed_from_checkpoint,
            seed_filter: seed_filter.clone(),
            align_config: align_config.clone(),
            read_count_total: reads_total,
            read_count_seed_passed: seed_passed,
            read_count_aligned: aligned,
            retained_count_msa_eligible,
            warnings,
            hits,
            exon_support_frequencies: mapped_exon_support_frequencies,
            junction_support_frequencies: mapped_junction_support_frequencies,
            transition_support_rows,
            isoform_support_rows,
            mapped_isoform_support_rows,
            origin_class_counts,
            score_density_bins,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{dna_sequence::DNAsequence, feature_expert::SplicingRange};

    fn range(start_1based: usize, end_1based: usize) -> SplicingRange {
        SplicingRange {
            start_1based,
            end_1based,
        }
    }

    #[test]
    fn make_transcript_template_reverse_strand_uses_reverse_complemented_exon_concat() {
        let dna = DNAsequence::from_sequence("AACCGGTT").expect("sequence");
        let lane = SplicingTranscriptLane {
            transcript_feature_id: 73,
            transcript_id: "tp73_tx".to_string(),
            label: "TP73 transcript".to_string(),
            strand: "-".to_string(),
            exons: vec![range(1, 2), range(5, 6)],
            exon_cds_phases: vec![],
            introns: vec![range(3, 4)],
            has_target_feature: true,
        };

        let template = GentleEngine::make_transcript_template(&dna, &lane, 2);

        assert_eq!(String::from_utf8_lossy(&template.sequence), "CCTT");
        assert_eq!(template.genomic_positions_1based, vec![6, 5, 2, 1]);

        let rows = GentleEngine::collect_rna_seed_hash_catalog_rows(&[template], 2);
        assert!(rows.iter().any(|row| {
            row.strand == "-"
                && row.kmer_sequence == "CC"
                && row.template_offset_0based == 0
                && row.genomic_pos_1based == 6
        }));
    }

    #[test]
    fn seed_hash_template_audit_rows_expose_exact_template_sequence_and_orientation() {
        let template = SplicingTranscriptTemplate {
            transcript_feature_id: 73,
            transcript_id: "tp73_tx".to_string(),
            transcript_label: "TP73 transcript".to_string(),
            strand: "-".to_string(),
            sequence: b"CCTT".to_vec(),
            genomic_positions_1based: vec![6, 5, 2, 1],
            kmer_positions: HashMap::new(),
        };

        let rows = GentleEngine::collect_rna_seed_hash_template_audit_rows(&[template]);
        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert_eq!(row.transcript_feature_id, 73);
        assert_eq!(row.transcript_id, "tp73_tx");
        assert_eq!(row.strand, "-");
        assert_eq!(row.template_sequence, "CCTT");
        assert_eq!(row.template_length_bp, 4);
        assert_eq!(row.template_first_genomic_pos_1based, 6);
        assert_eq!(row.template_last_genomic_pos_1based, 1);
        assert!(row.reverse_complemented_from_genome);
    }
}
