//! Sequencing-confirmation storage and called-read construct checks.
//!
//! Phase 1 focuses on deterministic confirmation from already-called read
//! sequences. Raw ABI/AB1/SCF intake and trace-aware enrichment follow later.

use super::*;

#[derive(Debug, Clone)]
pub(super) struct ComputedPairwiseAlignment {
    pub report: SequenceAlignmentReport,
    pub operations: Vec<bio::alignment::AlignmentOperation>,
    query_span_bases: Vec<u8>,
    target_span_bases: Vec<u8>,
}

#[derive(Debug, Clone, Default)]
struct TargetEvidenceStats {
    covered_bp: usize,
    contradicted: bool,
}

impl GentleEngine {
    pub(super) fn read_sequencing_confirmation_report_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> SequencingConfirmationReportStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<SequencingConfirmationReportStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = SEQUENCING_CONFIRMATION_REPORTS_SCHEMA.to_string();
        }
        store
    }

    pub(super) fn read_sequencing_confirmation_report_store(
        &self,
    ) -> SequencingConfirmationReportStore {
        Self::read_sequencing_confirmation_report_store_from_metadata(
            self.state
                .metadata
                .get(SEQUENCING_CONFIRMATION_REPORTS_METADATA_KEY),
        )
    }

    pub(super) fn write_sequencing_confirmation_report_store(
        &mut self,
        mut store: SequencingConfirmationReportStore,
    ) -> Result<(), EngineError> {
        if store.reports.is_empty() {
            self.state
                .metadata
                .remove(SEQUENCING_CONFIRMATION_REPORTS_METADATA_KEY);
            return Ok(());
        }
        store.schema = SEQUENCING_CONFIRMATION_REPORTS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize sequencing-confirmation metadata: {e}"),
        })?;
        self.state.metadata.insert(
            SEQUENCING_CONFIRMATION_REPORTS_METADATA_KEY.to_string(),
            value,
        );
        Ok(())
    }

    pub(super) fn upsert_sequencing_confirmation_report(
        &mut self,
        report: SequencingConfirmationReport,
    ) -> Result<(), EngineError> {
        let mut store = self.read_sequencing_confirmation_report_store();
        store.reports.insert(report.report_id.clone(), report);
        self.write_sequencing_confirmation_report_store(store)
    }

    pub(super) fn normalize_sequencing_confirmation_report_id(
        raw: &str,
    ) -> Result<String, EngineError> {
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

    fn unique_sequencing_confirmation_report_id(&self, base: &str) -> String {
        let normalized = Self::normalize_sequencing_confirmation_report_id(base)
            .unwrap_or_else(|_| "sequencing_confirmation".to_string());
        let store = self.read_sequencing_confirmation_report_store();
        if !store.reports.contains_key(&normalized) {
            return normalized;
        }
        let mut counter = 2usize;
        loop {
            let candidate = format!("{normalized}_{counter}");
            if !store.reports.contains_key(&candidate) {
                return candidate;
            }
            counter += 1;
        }
    }

    fn default_sequencing_confirmation_targets(
        expected_len: usize,
    ) -> Vec<SequencingConfirmationTargetSpec> {
        vec![SequencingConfirmationTargetSpec {
            target_id: "full_span".to_string(),
            label: "Full construct span".to_string(),
            kind: SequencingConfirmationTargetKind::FullSpan,
            start_0based: 0,
            end_0based_exclusive: expected_len,
            junction_left_end_0based: None,
            required: true,
        }]
    }

    fn normalize_sequencing_confirmation_targets(
        expected_len: usize,
        targets: &[SequencingConfirmationTargetSpec],
    ) -> Result<Vec<SequencingConfirmationTargetSpec>, EngineError> {
        let source_targets = if targets.is_empty() {
            Self::default_sequencing_confirmation_targets(expected_len)
        } else {
            targets.to_vec()
        };
        let mut normalized = Vec::with_capacity(source_targets.len());
        for (idx, target) in source_targets.into_iter().enumerate() {
            if target.end_0based_exclusive <= target.start_0based {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Sequencing-confirmation target {} has an empty interval {}..{}",
                        idx + 1,
                        target.start_0based,
                        target.end_0based_exclusive
                    ),
                });
            }
            if target.end_0based_exclusive > expected_len {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Sequencing-confirmation target {} interval {}..{} exceeds expected sequence length {}",
                        idx + 1,
                        target.start_0based,
                        target.end_0based_exclusive,
                        expected_len
                    ),
                });
            }
            let target_id = if target.target_id.trim().is_empty() {
                format!("{}_{}", target.kind.as_str(), idx + 1)
            } else {
                Self::normalize_sequencing_confirmation_report_id(target.target_id.as_str())?
            };
            let label = if target.label.trim().is_empty() {
                match target.kind {
                    SequencingConfirmationTargetKind::FullSpan => "Full construct span".to_string(),
                    SequencingConfirmationTargetKind::Junction => {
                        let left_end = target.junction_left_end_0based.unwrap_or(
                            target.start_0based.saturating_add(
                                (target.end_0based_exclusive - target.start_0based) / 2,
                            ),
                        );
                        format!("Junction @ {left_end}")
                    }
                    _ => format!(
                        "{} {}..{}",
                        target.kind.as_str(),
                        target.start_0based,
                        target.end_0based_exclusive
                    ),
                }
            } else {
                target.label.trim().to_string()
            };
            let junction_left_end_0based = match target.kind {
                SequencingConfirmationTargetKind::Junction => {
                    let left_end = target.junction_left_end_0based.ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Sequencing-confirmation target '{}' requires junction_left_end_0based",
                            target_id
                        ),
                    })?;
                    if left_end < target.start_0based || left_end > target.end_0based_exclusive {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Sequencing-confirmation target '{}' has junction_left_end_0based={} outside interval {}..{}",
                                target_id,
                                left_end,
                                target.start_0based,
                                target.end_0based_exclusive
                            ),
                        });
                    }
                    Some(left_end)
                }
                _ => target.junction_left_end_0based,
            };
            normalized.push(SequencingConfirmationTargetSpec {
                target_id,
                label,
                kind: target.kind,
                start_0based: target.start_0based,
                end_0based_exclusive: target.end_0based_exclusive,
                junction_left_end_0based,
                required: target.required,
            });
        }
        Ok(normalized)
    }

    pub(super) fn compute_pairwise_alignment_report(
        query_seq_id: &str,
        query_text: &str,
        query_span_start_0based: Option<usize>,
        query_span_end_0based: Option<usize>,
        target_seq_id: &str,
        target_text: &str,
        target_span_start_0based: Option<usize>,
        target_span_end_0based: Option<usize>,
        mode: PairwiseAlignmentMode,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Result<ComputedPairwiseAlignment, EngineError> {
        if match_score <= 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Pairwise alignment requires match_score > 0".to_string(),
            });
        }
        if gap_open > 0 || gap_extend > 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Pairwise alignment requires gap_open and gap_extend to be <= 0"
                    .to_string(),
            });
        }
        let query_text = query_text.to_ascii_uppercase();
        let target_text = target_text.to_ascii_uppercase();
        let query_bytes = query_text.as_bytes();
        let target_bytes = target_text.as_bytes();
        let (query_span_start_0based, query_span_end_0based) = Self::resolve_analysis_span(
            query_bytes.len(),
            query_span_start_0based,
            query_span_end_0based,
        )?;
        let (target_span_start_0based, target_span_end_0based) = Self::resolve_analysis_span(
            target_bytes.len(),
            target_span_start_0based,
            target_span_end_0based,
        )?;
        let query_span = &query_bytes[query_span_start_0based..query_span_end_0based];
        let target_span = &target_bytes[target_span_start_0based..target_span_end_0based];

        let score = |a: u8, b: u8| {
            if a.eq_ignore_ascii_case(&b) {
                match_score
            } else {
                mismatch_score
            }
        };
        let mut aligner = bio::alignment::pairwise::Aligner::new(gap_open, gap_extend, &score);
        let alignment = match mode {
            PairwiseAlignmentMode::Global => aligner.global(query_span, target_span),
            PairwiseAlignmentMode::Local => aligner.local(query_span, target_span),
        };

        let mut matches = 0usize;
        let mut mismatches = 0usize;
        let mut insertions = 0usize;
        let mut deletions = 0usize;
        let mut cigar = String::new();
        let mut run_len = 0usize;
        let mut run_code = 'M';
        for op in &alignment.operations {
            let code = match op {
                bio::alignment::AlignmentOperation::Match => {
                    matches += 1;
                    Some('=')
                }
                bio::alignment::AlignmentOperation::Subst => {
                    mismatches += 1;
                    Some('X')
                }
                bio::alignment::AlignmentOperation::Ins => {
                    insertions += 1;
                    Some('I')
                }
                bio::alignment::AlignmentOperation::Del => {
                    deletions += 1;
                    Some('D')
                }
                bio::alignment::AlignmentOperation::Xclip(_) => Some('S'),
                bio::alignment::AlignmentOperation::Yclip(_) => None,
            };
            let Some(code) = code else {
                continue;
            };
            if run_len == 0 {
                run_code = code;
                run_len = 1;
            } else if run_code == code {
                run_len += 1;
            } else {
                cigar.push_str(&format!("{run_len}{run_code}"));
                run_code = code;
                run_len = 1;
            }
        }
        if run_len > 0 {
            cigar.push_str(&format!("{run_len}{run_code}"));
        }
        if cigar.is_empty() {
            cigar = "*".to_string();
        }

        let aligned_columns = matches + mismatches + insertions + deletions;
        let query_covered = alignment.xend.saturating_sub(alignment.xstart);
        let target_covered = alignment.yend.saturating_sub(alignment.ystart);
        let identity_fraction = if aligned_columns == 0 {
            0.0
        } else {
            matches as f64 / aligned_columns as f64
        };
        let query_coverage_fraction = if query_span.is_empty() {
            0.0
        } else {
            query_covered as f64 / query_span.len() as f64
        };
        let target_coverage_fraction = if target_span.is_empty() {
            0.0
        } else {
            target_covered as f64 / target_span.len() as f64
        };
        let report = SequenceAlignmentReport {
            schema: SEQUENCE_ALIGNMENT_REPORT_SCHEMA.to_string(),
            mode,
            query_seq_id: query_seq_id.to_string(),
            target_seq_id: target_seq_id.to_string(),
            query_span_start_0based,
            query_span_end_0based,
            target_span_start_0based,
            target_span_end_0based,
            aligned_query_start_0based: query_span_start_0based + alignment.xstart,
            aligned_query_end_0based_exclusive: query_span_start_0based + alignment.xend,
            aligned_target_start_0based: target_span_start_0based + alignment.ystart,
            aligned_target_end_0based_exclusive: target_span_start_0based + alignment.yend,
            score: alignment.score,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            aligned_columns,
            matches,
            mismatches,
            insertions,
            deletions,
            identity_fraction,
            query_coverage_fraction,
            target_coverage_fraction,
            cigar,
        };
        Ok(ComputedPairwiseAlignment {
            report,
            operations: alignment.operations,
            query_span_bases: query_span.to_vec(),
            target_span_bases: target_span.to_vec(),
        })
    }

    fn extract_sequencing_confirmation_discrepancies(
        alignment: &ComputedPairwiseAlignment,
    ) -> Vec<SequencingConfirmationDiscrepancy> {
        let mut out = vec![];
        let mut query_pos = alignment
            .report
            .aligned_query_start_0based
            .saturating_sub(alignment.report.query_span_start_0based);
        let mut target_pos = alignment
            .report
            .aligned_target_start_0based
            .saturating_sub(alignment.report.target_span_start_0based);
        for op in &alignment.operations {
            match op {
                bio::alignment::AlignmentOperation::Match => {
                    query_pos += 1;
                    target_pos += 1;
                }
                bio::alignment::AlignmentOperation::Subst => {
                    let query_base = alignment
                        .query_span_bases
                        .get(query_pos)
                        .copied()
                        .unwrap_or(b'N') as char;
                    let target_base = alignment
                        .target_span_bases
                        .get(target_pos)
                        .copied()
                        .unwrap_or(b'N') as char;
                    out.push(SequencingConfirmationDiscrepancy {
                        kind: SequencingConfirmationDiscrepancyKind::Mismatch,
                        query_start_0based: alignment.report.query_span_start_0based + query_pos,
                        query_end_0based_exclusive: alignment.report.query_span_start_0based
                            + query_pos
                            + 1,
                        target_start_0based: alignment.report.target_span_start_0based + target_pos,
                        target_end_0based_exclusive: alignment.report.target_span_start_0based
                            + target_pos
                            + 1,
                        query_bases: query_base.to_string(),
                        target_bases: target_base.to_string(),
                    });
                    query_pos += 1;
                    target_pos += 1;
                }
                bio::alignment::AlignmentOperation::Ins => {
                    let query_base = alignment
                        .query_span_bases
                        .get(query_pos)
                        .copied()
                        .unwrap_or(b'N') as char;
                    out.push(SequencingConfirmationDiscrepancy {
                        kind: SequencingConfirmationDiscrepancyKind::Insertion,
                        query_start_0based: alignment.report.query_span_start_0based + query_pos,
                        query_end_0based_exclusive: alignment.report.query_span_start_0based
                            + query_pos
                            + 1,
                        target_start_0based: alignment.report.target_span_start_0based + target_pos,
                        target_end_0based_exclusive: alignment.report.target_span_start_0based
                            + target_pos,
                        query_bases: query_base.to_string(),
                        target_bases: String::new(),
                    });
                    query_pos += 1;
                }
                bio::alignment::AlignmentOperation::Del => {
                    let target_base = alignment
                        .target_span_bases
                        .get(target_pos)
                        .copied()
                        .unwrap_or(b'N') as char;
                    out.push(SequencingConfirmationDiscrepancy {
                        kind: SequencingConfirmationDiscrepancyKind::Deletion,
                        query_start_0based: alignment.report.query_span_start_0based + query_pos,
                        query_end_0based_exclusive: alignment.report.query_span_start_0based
                            + query_pos,
                        target_start_0based: alignment.report.target_span_start_0based + target_pos,
                        target_end_0based_exclusive: alignment.report.target_span_start_0based
                            + target_pos
                            + 1,
                        query_bases: String::new(),
                        target_bases: target_base.to_string(),
                    });
                    target_pos += 1;
                }
                bio::alignment::AlignmentOperation::Xclip(_) => {}
                bio::alignment::AlignmentOperation::Yclip(_) => {}
            }
        }
        out
    }

    fn collect_target_evidence_stats(
        alignment: &ComputedPairwiseAlignment,
        target: &SequencingConfirmationTargetSpec,
    ) -> TargetEvidenceStats {
        let mut stats = TargetEvidenceStats::default();
        let mut target_pos = alignment
            .report
            .aligned_target_start_0based
            .saturating_sub(alignment.report.target_span_start_0based);
        for op in &alignment.operations {
            match op {
                bio::alignment::AlignmentOperation::Match => {
                    let target_abs = alignment.report.target_span_start_0based + target_pos;
                    if target_abs >= target.start_0based && target_abs < target.end_0based_exclusive
                    {
                        stats.covered_bp += 1;
                    }
                    target_pos += 1;
                }
                bio::alignment::AlignmentOperation::Subst => {
                    let target_abs = alignment.report.target_span_start_0based + target_pos;
                    if target_abs >= target.start_0based && target_abs < target.end_0based_exclusive
                    {
                        stats.covered_bp += 1;
                        stats.contradicted = true;
                    }
                    target_pos += 1;
                }
                bio::alignment::AlignmentOperation::Ins => {
                    let target_abs = alignment.report.target_span_start_0based + target_pos;
                    if target_abs >= target.start_0based
                        && target_abs <= target.end_0based_exclusive
                    {
                        stats.contradicted = true;
                    }
                }
                bio::alignment::AlignmentOperation::Del => {
                    let target_abs = alignment.report.target_span_start_0based + target_pos;
                    if target_abs >= target.start_0based && target_abs < target.end_0based_exclusive
                    {
                        stats.covered_bp += 1;
                        stats.contradicted = true;
                    }
                    target_pos += 1;
                }
                bio::alignment::AlignmentOperation::Xclip(_) => {}
                bio::alignment::AlignmentOperation::Yclip(_) => {}
            }
        }
        stats
    }

    pub(super) fn confirm_construct_reads(
        &mut self,
        expected_seq_id: &str,
        read_seq_ids: &[String],
        targets: &[SequencingConfirmationTargetSpec],
        alignment_mode: PairwiseAlignmentMode,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
        min_identity_fraction: f64,
        min_target_coverage_fraction: f64,
        allow_reverse_complement: bool,
        report_id: Option<&str>,
    ) -> Result<SequencingConfirmationReport, EngineError> {
        if read_seq_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ConfirmConstructReads requires at least one read sequence".to_string(),
            });
        }
        if !(0.0..=1.0).contains(&min_identity_fraction) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "ConfirmConstructReads requires min_identity_fraction in the range 0.0..=1.0"
                        .to_string(),
            });
        }
        if !(0.0..=1.0).contains(&min_target_coverage_fraction)
            || min_target_coverage_fraction == 0.0
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ConfirmConstructReads requires min_target_coverage_fraction in the range (0.0, 1.0]"
                    .to_string(),
            });
        }
        let expected_dna =
            self.state
                .sequences
                .get(expected_seq_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Sequence '{expected_seq_id}' not found"),
                })?;
        let expected_text = expected_dna.get_forward_string().to_ascii_uppercase();
        if expected_text.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "ConfirmConstructReads cannot use empty expected sequence '{}'",
                    expected_seq_id
                ),
            });
        }
        let normalized_targets =
            Self::normalize_sequencing_confirmation_targets(expected_text.len(), targets)?;
        let report_id = match report_id {
            Some(raw) => {
                let normalized = Self::normalize_sequencing_confirmation_report_id(raw)?;
                self.unique_sequencing_confirmation_report_id(&normalized)
            }
            None => self.unique_sequencing_confirmation_report_id(&format!(
                "{expected_seq_id}_seq_confirm"
            )),
        };

        let mut reads = vec![];
        let mut warnings = vec![];
        let mut target_results = normalized_targets
            .iter()
            .map(|target| SequencingConfirmationTargetResult {
                target_id: target.target_id.clone(),
                label: target.label.clone(),
                kind: target.kind,
                start_0based: target.start_0based,
                end_0based_exclusive: target.end_0based_exclusive,
                junction_left_end_0based: target.junction_left_end_0based,
                required: target.required,
                status: SequencingConfirmationStatus::InsufficientEvidence,
                covered_bp: 0,
                target_length_bp: target
                    .end_0based_exclusive
                    .saturating_sub(target.start_0based),
                support_read_ids: vec![],
                contradicting_read_ids: vec![],
                reason: String::new(),
            })
            .collect::<Vec<_>>();

        for read_seq_id in read_seq_ids {
            let read_dna = self
                .state
                .sequences
                .get(read_seq_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Sequence '{read_seq_id}' not found"),
                })?;
            let read_text = read_dna.get_forward_string().to_ascii_uppercase();
            if read_text.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "ConfirmConstructReads cannot use empty read sequence '{}'",
                        read_seq_id
                    ),
                });
            }
            let forward_alignment = Self::compute_pairwise_alignment_report(
                read_seq_id,
                read_text.as_str(),
                None,
                None,
                expected_seq_id,
                expected_text.as_str(),
                None,
                None,
                alignment_mode,
                match_score,
                mismatch_score,
                gap_open,
                gap_extend,
            )?;
            let reverse_alignment = if allow_reverse_complement {
                let reverse_text = Self::reverse_complement(&read_text);
                Some(Self::compute_pairwise_alignment_report(
                    read_seq_id,
                    reverse_text.as_str(),
                    None,
                    None,
                    expected_seq_id,
                    expected_text.as_str(),
                    None,
                    None,
                    alignment_mode,
                    match_score,
                    mismatch_score,
                    gap_open,
                    gap_extend,
                )?)
            } else {
                None
            };
            let (orientation, alignment) = if let Some(reverse_alignment) = reverse_alignment {
                let forward_score = (
                    (forward_alignment.report.query_coverage_fraction * 1_000_000.0) as i64,
                    (forward_alignment.report.identity_fraction * 1_000_000.0) as i64,
                    i64::from(forward_alignment.report.score),
                    forward_alignment.report.aligned_columns as i64,
                );
                let reverse_score = (
                    (reverse_alignment.report.query_coverage_fraction * 1_000_000.0) as i64,
                    (reverse_alignment.report.identity_fraction * 1_000_000.0) as i64,
                    i64::from(reverse_alignment.report.score),
                    reverse_alignment.report.aligned_columns as i64,
                );
                if reverse_score > forward_score {
                    (
                        SequencingReadOrientation::ReverseComplement,
                        reverse_alignment,
                    )
                } else {
                    (SequencingReadOrientation::Forward, forward_alignment)
                }
            } else {
                (SequencingReadOrientation::Forward, forward_alignment)
            };
            let discrepancies = Self::extract_sequencing_confirmation_discrepancies(&alignment);
            let mut covered_target_ids = vec![];
            let mut confirmed_target_ids = vec![];
            let mut contradicted_target_ids = vec![];
            for (target_result, target_spec) in
                target_results.iter_mut().zip(normalized_targets.iter())
            {
                let stats = Self::collect_target_evidence_stats(&alignment, target_spec);
                if stats.covered_bp > target_result.covered_bp {
                    target_result.covered_bp = stats.covered_bp;
                }
                if stats.covered_bp > 0 {
                    covered_target_ids.push(target_spec.target_id.clone());
                }
                let target_len = target_result.target_length_bp.max(1);
                let coverage_fraction = stats.covered_bp as f64 / target_len as f64;
                if stats.contradicted {
                    target_result
                        .contradicting_read_ids
                        .push(read_seq_id.clone());
                    contradicted_target_ids.push(target_spec.target_id.clone());
                } else if coverage_fraction >= min_target_coverage_fraction
                    && alignment.report.identity_fraction >= min_identity_fraction
                {
                    target_result.support_read_ids.push(read_seq_id.clone());
                    confirmed_target_ids.push(target_spec.target_id.clone());
                }
            }
            if alignment.report.identity_fraction < min_identity_fraction {
                warnings.push(format!(
                    "Read '{}' best alignment identity {:.3} is below min_identity_fraction {:.3}",
                    read_seq_id, alignment.report.identity_fraction, min_identity_fraction
                ));
            }
            reads.push(SequencingConfirmationReadResult {
                read_seq_id: read_seq_id.clone(),
                orientation,
                usable: !covered_target_ids.is_empty(),
                best_alignment: alignment.report,
                covered_target_ids,
                confirmed_target_ids,
                contradicted_target_ids,
                discrepancies,
            });
        }

        for target in &mut target_results {
            target.status = if !target.contradicting_read_ids.is_empty() {
                SequencingConfirmationStatus::Contradicted
            } else if !target.support_read_ids.is_empty() {
                SequencingConfirmationStatus::Confirmed
            } else {
                SequencingConfirmationStatus::InsufficientEvidence
            };
            let coverage_pct = if target.target_length_bp == 0 {
                0.0
            } else {
                (target.covered_bp as f64 / target.target_length_bp as f64) * 100.0
            };
            target.reason = match target.status {
                SequencingConfirmationStatus::Confirmed => format!(
                    "{} supporting read(s); best target coverage {:.1}%",
                    target.support_read_ids.len(),
                    coverage_pct
                ),
                SequencingConfirmationStatus::Contradicted => format!(
                    "{} contradicting read(s) overlap this target",
                    target.contradicting_read_ids.len()
                ),
                SequencingConfirmationStatus::InsufficientEvidence => format!(
                    "No read met target coverage >= {:.1}% without contradiction (best coverage {:.1}%)",
                    min_target_coverage_fraction * 100.0,
                    coverage_pct
                ),
            };
        }

        let overall_status = if target_results.iter().any(|target| {
            target.required && target.status == SequencingConfirmationStatus::Contradicted
        }) {
            SequencingConfirmationStatus::Contradicted
        } else if target_results
            .iter()
            .filter(|target| target.required)
            .all(|target| target.status == SequencingConfirmationStatus::Confirmed)
        {
            SequencingConfirmationStatus::Confirmed
        } else {
            SequencingConfirmationStatus::InsufficientEvidence
        };

        let report = SequencingConfirmationReport {
            schema: SEQUENCING_CONFIRMATION_REPORT_SCHEMA.to_string(),
            report_id,
            expected_seq_id: expected_seq_id.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            overall_status,
            alignment_mode,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            min_identity_fraction,
            min_target_coverage_fraction,
            allow_reverse_complement,
            read_seq_ids: read_seq_ids.to_vec(),
            target_count: target_results.len(),
            reads,
            targets: target_results,
            warnings,
        };
        self.upsert_sequencing_confirmation_report(report.clone())?;
        Ok(report)
    }

    pub fn list_sequencing_confirmation_reports(
        &self,
        expected_seq_id_filter: Option<&str>,
    ) -> Vec<SequencingConfirmationReportSummary> {
        let filter = expected_seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let mut rows = self
            .read_sequencing_confirmation_report_store()
            .reports
            .values()
            .filter(|report| {
                filter.as_ref().is_none_or(|needle| {
                    report
                        .expected_seq_id
                        .to_ascii_lowercase()
                        .eq_ignore_ascii_case(needle)
                })
            })
            .map(|report| SequencingConfirmationReportSummary {
                report_id: report.report_id.clone(),
                expected_seq_id: report.expected_seq_id.clone(),
                generated_at_unix_ms: report.generated_at_unix_ms,
                overall_status: report.overall_status,
                read_count: report.reads.len(),
                target_count: report.targets.len(),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.expected_seq_id
                .to_ascii_lowercase()
                .cmp(&right.expected_seq_id.to_ascii_lowercase())
                .then(left.generated_at_unix_ms.cmp(&right.generated_at_unix_ms))
                .then(
                    left.report_id
                        .to_ascii_lowercase()
                        .cmp(&right.report_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub(crate) fn format_sequencing_confirmation_report_summary_row(
        row: &SequencingConfirmationReportSummary,
    ) -> String {
        format!(
            "{} expected={} status={} reads={} targets={}",
            row.report_id,
            row.expected_seq_id,
            row.overall_status.as_str(),
            row.read_count,
            row.target_count
        )
    }

    pub(crate) fn format_sequencing_confirmation_report_detail_summary(
        report: &SequencingConfirmationReport,
    ) -> String {
        format!(
            "{} expected={} status={} reads={} targets={} mode={} rc_allowed={} identity>={:.2} coverage>={:.2}",
            report.report_id,
            report.expected_seq_id,
            report.overall_status.as_str(),
            report.reads.len(),
            report.targets.len(),
            report.alignment_mode.as_str(),
            report.allow_reverse_complement,
            report.min_identity_fraction,
            report.min_target_coverage_fraction
        )
    }

    pub fn get_sequencing_confirmation_report(
        &self,
        report_id: &str,
    ) -> Result<SequencingConfirmationReport, EngineError> {
        let report_id = Self::normalize_sequencing_confirmation_report_id(report_id)?;
        self.read_sequencing_confirmation_report_store()
            .reports
            .get(report_id.as_str())
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequencing-confirmation report '{}' not found", report_id),
            })
    }

    pub fn export_sequencing_confirmation_report(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<SequencingConfirmationReport, EngineError> {
        let report = self.get_sequencing_confirmation_report(report_id)?;
        let text = serde_json::to_string_pretty(&report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize sequencing-confirmation report '{}': {e}",
                report.report_id
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write sequencing-confirmation report to '{path}': {e}"),
        })?;
        Ok(report)
    }

    pub fn export_sequencing_confirmation_support_tsv(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<SequencingConfirmationReport, EngineError> {
        let report = self.get_sequencing_confirmation_report(report_id)?;
        let mut writer = BufWriter::new(File::create(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create sequencing-confirmation TSV output '{}': {e}",
                path
            ),
        })?);
        writeln!(
            writer,
            "report_id\texpected_seq_id\toverall_status\ttarget_id\tlabel\tkind\trequired\tstatus\tcovered_bp\ttarget_length_bp\tsupport_read_ids\tcontradicting_read_ids\treason"
        )
        .map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write sequencing-confirmation TSV header '{}': {e}",
                path
            ),
        })?;
        for target in &report.targets {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                report.report_id,
                report.expected_seq_id,
                report.overall_status.as_str(),
                target.target_id,
                Self::sanitize_tsv_cell(&target.label),
                target.kind.as_str(),
                target.required,
                target.status.as_str(),
                target.covered_bp,
                target.target_length_bp,
                Self::sanitize_tsv_cell(&target.support_read_ids.join(",")),
                Self::sanitize_tsv_cell(&target.contradicting_read_ids.join(",")),
                Self::sanitize_tsv_cell(&target.reason),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not write sequencing-confirmation TSV row '{}': {e}",
                    path
                ),
            })?;
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not flush sequencing-confirmation TSV output '{}': {e}",
                path
            ),
        })?;
        Ok(report)
    }
}
