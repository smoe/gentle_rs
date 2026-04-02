//! Sequencing-confirmation storage and construct checks from read or trace evidence.
//!
//! The baseline workflow started with deterministic confirmation from already-
//! called read sequences. Imported ABI/AB1/SCF traces now participate as
//! evidence too by reusing their stored called bases while keeping the raw
//! trace store separate from confirmation reports.
//!
//! Look here for:
//! - sequencing-confirmation report-store helpers
//! - target-by-target evidence aggregation from pairwise alignments
//! - JSON/TSV-facing construct-confirmation exports and verdict logic

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

#[derive(Debug, Clone, Default)]
struct VariantConfidenceSummary {
    min: Option<u8>,
    max: Option<u8>,
    mean: Option<f64>,
    count: usize,
}

#[derive(Debug, Clone)]
struct VariantLocus {
    variant_id: String,
    label: String,
    target_id: Option<String>,
    start_0based: usize,
    end_0based_exclusive: usize,
    expected_bases: String,
    baseline_bases: Option<String>,
}

#[derive(Debug, Clone)]
struct VariantObservation {
    classification: SequencingConfirmationVariantClassification,
    status: SequencingConfirmationStatus,
    observed_bases: String,
    confidence: VariantConfidenceSummary,
    peak_center: Option<u32>,
    reason: String,
}

#[derive(Debug, Clone)]
struct SequencingConfirmationEvidenceInput {
    evidence_kind: SequencingConfirmationEvidenceKind,
    evidence_id: String,
    display_read_seq_id: String,
    trace_id: Option<String>,
    linked_seq_id: Option<String>,
    called_bases: String,
    confidence_values: Option<Vec<u8>>,
    peak_locations: Option<Vec<u32>>,
}

#[derive(Debug, Clone)]
struct SequencingPrimerProposalCandidate {
    orientation: SequencingPrimerOrientation,
    primer: PrimerDesignPrimerRecord,
    anneal_sequence: String,
    three_prime_position_0based: usize,
    predicted_read_span_start_0based: usize,
    predicted_read_span_end_0based_exclusive: usize,
    three_prime_distance_bp: usize,
    score: f64,
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

    pub(crate) fn sequencing_confirmation_reports_from_state(
        state: &ProjectState,
    ) -> Vec<SequencingConfirmationReport> {
        let mut reports = Self::read_sequencing_confirmation_report_store_from_metadata(
            state
                .metadata
                .get(SEQUENCING_CONFIRMATION_REPORTS_METADATA_KEY),
        )
        .reports
        .into_values()
        .collect::<Vec<_>>();
        reports.sort_by(|left, right| {
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
        reports
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
            expected_bases: None,
            baseline_bases: None,
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
            if target.end_0based_exclusive < target.start_0based
                || (target.kind != SequencingConfirmationTargetKind::ExpectedEdit
                    && target.end_0based_exclusive == target.start_0based)
            {
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
                expected_bases: target.expected_bases.clone(),
                baseline_bases: target.baseline_bases.clone(),
                required: target.required,
            });
        }
        Ok(normalized)
    }

    fn sequencing_trace_variant_low_confidence_threshold() -> u8 {
        20
    }

    fn merge_sequencing_confirmation_discrepancies(
        discrepancies: &[SequencingConfirmationDiscrepancy],
    ) -> Vec<SequencingConfirmationDiscrepancy> {
        let mut out: Vec<SequencingConfirmationDiscrepancy> =
            Vec::with_capacity(discrepancies.len());
        for row in discrepancies {
            let Some(last) = out.last_mut() else {
                out.push(row.clone());
                continue;
            };
            let can_merge = last.kind == row.kind
                && last.query_end_0based_exclusive == row.query_start_0based
                && last.target_end_0based_exclusive == row.target_start_0based;
            if can_merge {
                last.query_end_0based_exclusive = row.query_end_0based_exclusive;
                last.target_end_0based_exclusive = row.target_end_0based_exclusive;
                last.query_bases.push_str(&row.query_bases);
                last.target_bases.push_str(&row.target_bases);
            } else {
                out.push(row.clone());
            }
        }
        out
    }

    fn infer_expected_edit_targets(
        expected_seq_id: &str,
        expected_text: &str,
        baseline_seq_id: &str,
        baseline_text: &str,
    ) -> Result<Vec<SequencingConfirmationTargetSpec>, EngineError> {
        let alignment = Self::compute_pairwise_alignment_report(
            baseline_seq_id,
            baseline_text,
            None,
            None,
            expected_seq_id,
            expected_text,
            None,
            None,
            PairwiseAlignmentMode::Global,
            2,
            -3,
            -5,
            -1,
        )?;
        let discrepancies = Self::extract_sequencing_confirmation_discrepancies(&alignment);
        let merged = Self::merge_sequencing_confirmation_discrepancies(&discrepancies);
        Ok(merged
            .into_iter()
            .enumerate()
            .map(|(idx, row)| {
                let expected_bases = row.target_bases.to_ascii_uppercase();
                let baseline_bases = row.query_bases.to_ascii_uppercase();
                let label = if expected_bases.is_empty() {
                    format!(
                        "Expected deletion @ {} (remove '{}')",
                        row.target_start_0based, baseline_bases
                    )
                } else if baseline_bases.is_empty() {
                    format!(
                        "Expected insertion @ {} ('{}')",
                        row.target_start_0based, expected_bases
                    )
                } else {
                    format!(
                        "Expected edit @ {}..{} ({} -> {})",
                        row.target_start_0based,
                        row.target_end_0based_exclusive,
                        baseline_bases,
                        expected_bases
                    )
                };
                SequencingConfirmationTargetSpec {
                    target_id: format!("expected_edit_{}", idx + 1),
                    label,
                    kind: SequencingConfirmationTargetKind::ExpectedEdit,
                    start_0based: row.target_start_0based,
                    end_0based_exclusive: row.target_end_0based_exclusive,
                    junction_left_end_0based: None,
                    expected_bases: Some(expected_bases),
                    baseline_bases: Some(baseline_bases),
                    required: true,
                }
            })
            .collect())
    }

    fn collect_variant_loci(
        expected_text: &str,
        targets: &[SequencingConfirmationTargetSpec],
    ) -> Vec<VariantLocus> {
        targets
            .iter()
            .filter(|target| target.kind == SequencingConfirmationTargetKind::ExpectedEdit)
            .map(|target| VariantLocus {
                variant_id: target.target_id.clone(),
                label: target.label.clone(),
                target_id: Some(target.target_id.clone()),
                start_0based: target.start_0based,
                end_0based_exclusive: target.end_0based_exclusive,
                expected_bases: target
                    .expected_bases
                    .clone()
                    .unwrap_or_else(|| {
                        expected_text[target.start_0based..target.end_0based_exclusive].to_string()
                    })
                    .to_ascii_uppercase(),
                baseline_bases: target
                    .baseline_bases
                    .clone()
                    .map(|value| value.to_ascii_uppercase()),
            })
            .collect()
    }

    fn sequencing_primer_spans_overlap(
        left_start: usize,
        left_end_exclusive: usize,
        right_start: usize,
        right_end_exclusive: usize,
    ) -> bool {
        if left_start == left_end_exclusive {
            return left_start >= right_start && left_start < right_end_exclusive;
        }
        if right_start == right_end_exclusive {
            return right_start >= left_start && right_start < left_end_exclusive;
        }
        left_start < right_end_exclusive && right_start < left_end_exclusive
    }

    fn sequencing_primer_problem_center(
        locus_start_0based: usize,
        locus_end_0based_exclusive: usize,
    ) -> usize {
        if locus_end_0based_exclusive > locus_start_0based {
            locus_start_0based + ((locus_end_0based_exclusive - locus_start_0based) / 2)
        } else {
            locus_start_0based
        }
    }

    fn sequencing_primer_three_prime_distance_bp(
        suggestion: &SequencingPrimerOverlaySuggestion,
        locus_center_0based: usize,
    ) -> (bool, usize) {
        match suggestion.orientation {
            SequencingPrimerOrientation::ForwardRead => (
                locus_center_0based >= suggestion.three_prime_position_0based,
                locus_center_0based.abs_diff(suggestion.three_prime_position_0based),
            ),
            SequencingPrimerOrientation::ReverseRead => (
                locus_center_0based <= suggestion.three_prime_position_0based,
                locus_center_0based.abs_diff(suggestion.three_prime_position_0based),
            ),
        }
    }

    fn sequencing_primer_guidance_penalty(
        in_read_direction: bool,
        three_prime_distance_bp: usize,
    ) -> usize {
        const MIN_PREFERRED_DISTANCE_BP: usize = 40;
        const MAX_PREFERRED_DISTANCE_BP: usize = 500;
        let direction_penalty = if in_read_direction { 0 } else { 1000 };
        let distance_penalty = if three_prime_distance_bp < MIN_PREFERRED_DISTANCE_BP {
            MIN_PREFERRED_DISTANCE_BP - three_prime_distance_bp
        } else if three_prime_distance_bp > MAX_PREFERRED_DISTANCE_BP {
            three_prime_distance_bp - MAX_PREFERRED_DISTANCE_BP
        } else {
            0
        };
        direction_penalty + distance_penalty
    }

    fn sequencing_primer_guidance_reason(
        candidate_count: usize,
        in_read_direction: bool,
        three_prime_distance_bp: usize,
        problem_target_count: usize,
        problem_variant_count: usize,
    ) -> String {
        let total_problem_count = problem_target_count + problem_variant_count;
        let coverage_text = if total_problem_count > 1 {
            format!(
                " The same predicted read span also covers {} unresolved loci in total.",
                total_problem_count
            )
        } else {
            String::new()
        };
        if in_read_direction {
            if candidate_count > 1 {
                format!(
                    "Ranked best among {candidate_count} covering primer hits because the locus sits {} bp from the primer 3' end.{}",
                    three_prime_distance_bp, coverage_text
                )
            } else {
                format!(
                    "Only covering primer hit; the locus sits {} bp from the primer 3' end.{}",
                    three_prime_distance_bp, coverage_text
                )
            }
        } else if candidate_count > 1 {
            format!(
                "Ranks best among {candidate_count} covering primer hits, but the locus is not downstream of the primer 3' end in read direction (distance {} bp).{}",
                three_prime_distance_bp, coverage_text
            )
        } else {
            format!(
                "Only covering primer hit, but the locus is not downstream of the primer 3' end in read direction (distance {} bp).{}",
                three_prime_distance_bp, coverage_text
            )
        }
    }

    fn sequencing_primer_proposal_distance_window(
        predicted_read_length_bp: usize,
    ) -> (usize, usize, usize) {
        let min_distance_bp = predicted_read_length_bp.saturating_sub(1).min(40);
        let max_distance_bp = predicted_read_length_bp
            .saturating_sub(1)
            .min(500)
            .max(min_distance_bp);
        let target_distance_bp = 120usize.min(max_distance_bp).max(min_distance_bp);
        (min_distance_bp, max_distance_bp, target_distance_bp)
    }

    fn sequencing_primer_guidance_is_good_enough(
        row: &SequencingPrimerProblemGuidanceRow,
        predicted_read_length_bp: usize,
    ) -> bool {
        let (min_distance_bp, max_distance_bp, _) =
            Self::sequencing_primer_proposal_distance_window(predicted_read_length_bp);
        matches!(
            (
                row.recommended_primer_seq_id.as_deref(),
                row.recommended_in_read_direction,
                row.recommended_three_prime_distance_bp,
            ),
            (Some(_), Some(true), Some(distance_bp))
                if distance_bp >= min_distance_bp && distance_bp <= max_distance_bp
        )
    }

    fn sequencing_primer_proposal_score(
        primer: &PrimerDesignPrimerRecord,
        three_prime_distance_bp: usize,
        target_distance_bp: usize,
    ) -> f64 {
        let distance_penalty = three_prime_distance_bp.abs_diff(target_distance_bp) as f64 * 2.0;
        let tm_penalty = (primer.tm_c - 60.0).abs() * 12.0;
        let gc_penalty = (primer.gc_fraction - 0.5).abs() * 80.0;
        let length_penalty = Self::preferred_primer_length_penalty(primer.anneal_length_bp) * 6.0;
        let hit_penalty = primer.anneal_hits.saturating_sub(1) as f64 * 25.0;
        let homopolymer_penalty = if primer.longest_homopolymer_run_bp
            > PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP
        {
            (primer.longest_homopolymer_run_bp - PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP) as f64
                * 10.0
        } else {
            0.0
        };
        let self_complementary_penalty = if primer.self_complementary_run_bp
            > PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP
        {
            (primer.self_complementary_run_bp - PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP)
                as f64
                * 12.0
        } else {
            0.0
        };
        let gc_clamp_bonus = if primer.three_prime_gc_clamp {
            8.0
        } else {
            -8.0
        };
        1000.0
            - distance_penalty
            - tm_penalty
            - gc_penalty
            - length_penalty
            - hit_penalty
            - homopolymer_penalty
            - self_complementary_penalty
            + gc_clamp_bonus
    }

    fn build_sequencing_primer_proposal_candidates(
        expected_text: &str,
        locus_start_0based: usize,
        locus_end_0based_exclusive: usize,
        predicted_read_length_bp: usize,
    ) -> Vec<SequencingPrimerProposalCandidate> {
        let template = expected_text.as_bytes();
        if template.is_empty() || predicted_read_length_bp == 0 {
            return vec![];
        }
        let locus_center_0based =
            Self::sequencing_primer_problem_center(locus_start_0based, locus_end_0based_exclusive)
                .min(template.len().saturating_sub(1));
        let default_side = PrimerDesignSideConstraint::default();
        let (min_distance_bp, max_distance_bp, target_distance_bp) =
            Self::sequencing_primer_proposal_distance_window(predicted_read_length_bp);
        let mut ret = Vec::new();

        for length_bp in default_side.min_length..=default_side.max_length {
            if length_bp > template.len() {
                break;
            }
            for distance_bp in min_distance_bp..=max_distance_bp {
                let Some(anneal_end_0based_exclusive) = locus_center_0based
                    .checked_add(1)
                    .and_then(|value| value.checked_sub(distance_bp))
                else {
                    continue;
                };
                if anneal_end_0based_exclusive == 0 || anneal_end_0based_exclusive > template.len()
                {
                    continue;
                }
                let Some(anneal_start_0based) = anneal_end_0based_exclusive.checked_sub(length_bp)
                else {
                    continue;
                };
                let binding_window = &template[anneal_start_0based..anneal_end_0based_exclusive];
                let anneal_hits = Self::find_all_subsequences(template, binding_window).len();
                if anneal_hits != 1 {
                    continue;
                }
                let anneal_sequence =
                    String::from_utf8(binding_window.to_vec()).unwrap_or_default();
                let tm_c = Self::estimate_primer_tm_c(binding_window);
                let gc_fraction = Self::sequence_gc_fraction(binding_window).unwrap_or(0.0);
                if gc_fraction < default_side.min_gc_fraction
                    || gc_fraction > default_side.max_gc_fraction
                    || tm_c < default_side.min_tm_c
                    || tm_c > default_side.max_tm_c
                {
                    continue;
                }
                let predicted_read_span_start_0based = anneal_start_0based;
                let predicted_read_span_end_0based_exclusive =
                    (anneal_start_0based + predicted_read_length_bp).min(template.len());
                if locus_center_0based < predicted_read_span_start_0based
                    || locus_center_0based >= predicted_read_span_end_0based_exclusive
                {
                    continue;
                }
                let primer = Self::annotate_primer_record_heuristics(
                    PrimerDesignPrimerRecord {
                        sequence: anneal_sequence.clone(),
                        start_0based: anneal_start_0based,
                        end_0based_exclusive: anneal_end_0based_exclusive,
                        tm_c,
                        gc_fraction,
                        anneal_hits,
                        ..PrimerDesignPrimerRecord::default()
                    },
                    length_bp,
                );
                let score = Self::sequencing_primer_proposal_score(
                    &primer,
                    distance_bp,
                    target_distance_bp,
                );
                ret.push(SequencingPrimerProposalCandidate {
                    orientation: SequencingPrimerOrientation::ForwardRead,
                    primer,
                    anneal_sequence,
                    three_prime_position_0based: anneal_end_0based_exclusive.saturating_sub(1),
                    predicted_read_span_start_0based,
                    predicted_read_span_end_0based_exclusive,
                    three_prime_distance_bp: distance_bp,
                    score,
                });
            }

            for distance_bp in min_distance_bp..=max_distance_bp {
                let anneal_start_0based = locus_center_0based.saturating_add(distance_bp);
                let Some(anneal_end_0based_exclusive) = anneal_start_0based.checked_add(length_bp)
                else {
                    continue;
                };
                if anneal_end_0based_exclusive > template.len() {
                    continue;
                }
                let binding_window = &template[anneal_start_0based..anneal_end_0based_exclusive];
                let anneal_hits = Self::find_all_subsequences(template, binding_window).len();
                if anneal_hits != 1 {
                    continue;
                }
                let anneal_bytes = Self::reverse_complement_bytes(binding_window);
                let anneal_sequence = String::from_utf8(anneal_bytes).unwrap_or_default();
                let tm_c = Self::estimate_primer_tm_c(anneal_sequence.as_bytes());
                let gc_fraction =
                    Self::sequence_gc_fraction(anneal_sequence.as_bytes()).unwrap_or(0.0);
                if gc_fraction < default_side.min_gc_fraction
                    || gc_fraction > default_side.max_gc_fraction
                    || tm_c < default_side.min_tm_c
                    || tm_c > default_side.max_tm_c
                {
                    continue;
                }
                let predicted_read_span_end_0based_exclusive = anneal_end_0based_exclusive;
                let predicted_read_span_start_0based =
                    anneal_end_0based_exclusive.saturating_sub(predicted_read_length_bp);
                if locus_center_0based < predicted_read_span_start_0based
                    || locus_center_0based >= predicted_read_span_end_0based_exclusive
                {
                    continue;
                }
                let primer = Self::annotate_primer_record_heuristics(
                    PrimerDesignPrimerRecord {
                        sequence: anneal_sequence.clone(),
                        start_0based: anneal_start_0based,
                        end_0based_exclusive: anneal_end_0based_exclusive,
                        tm_c,
                        gc_fraction,
                        anneal_hits,
                        ..PrimerDesignPrimerRecord::default()
                    },
                    length_bp,
                );
                let score = Self::sequencing_primer_proposal_score(
                    &primer,
                    distance_bp,
                    target_distance_bp,
                );
                ret.push(SequencingPrimerProposalCandidate {
                    orientation: SequencingPrimerOrientation::ReverseRead,
                    primer,
                    anneal_sequence,
                    three_prime_position_0based: anneal_start_0based,
                    predicted_read_span_start_0based,
                    predicted_read_span_end_0based_exclusive,
                    three_prime_distance_bp: distance_bp,
                    score,
                });
            }
        }

        ret.sort_by(|left, right| {
            right
                .score
                .total_cmp(&left.score)
                .then(
                    left.three_prime_distance_bp
                        .cmp(&right.three_prime_distance_bp),
                )
                .then(left.orientation.as_str().cmp(right.orientation.as_str()))
                .then(left.primer.start_0based.cmp(&right.primer.start_0based))
                .then(left.primer.sequence.cmp(&right.primer.sequence))
        });
        ret.dedup_by(|left, right| {
            left.orientation == right.orientation
                && left.primer.start_0based == right.primer.start_0based
                && left.primer.end_0based_exclusive == right.primer.end_0based_exclusive
                && left.primer.sequence == right.primer.sequence
        });
        ret
    }

    fn sequencing_primer_proposal_reason(
        problem_summary: &str,
        replaced_existing_primer_seq_id: Option<&str>,
        candidate_count: usize,
        three_prime_distance_bp: usize,
    ) -> String {
        match replaced_existing_primer_seq_id {
            Some(existing) => format!(
                "Existing primer '{}' covered this {} locus but not within the preferred sequencing window, so the best of {candidate_count} fresh proposals was chosen {} bp from the new primer 3' end.",
                existing, problem_summary, three_prime_distance_bp
            ),
            None => format!(
                "No existing primer covered this {} locus, so the best of {candidate_count} fresh proposals was chosen {} bp from the new primer 3' end.",
                problem_summary, three_prime_distance_bp
            ),
        }
    }

    pub fn suggest_sequencing_primers(
        &self,
        expected_seq_id: &str,
        primer_seq_ids: &[String],
        confirmation_report_id: Option<&str>,
        min_3prime_anneal_bp: usize,
        predicted_read_length_bp: usize,
    ) -> Result<SequencingPrimerOverlayReport, EngineError> {
        if primer_seq_ids.is_empty() && confirmation_report_id.is_none() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "SuggestSequencingPrimers requires at least one primer sequence id or a saved sequencing-confirmation report"
                        .to_string(),
            });
        }
        if min_3prime_anneal_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SuggestSequencingPrimers requires min_3prime_anneal_bp >= 1".to_string(),
            });
        }
        if predicted_read_length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SuggestSequencingPrimers requires predicted_read_length_bp >= 1"
                    .to_string(),
            });
        }
        let expected_seq_id = expected_seq_id.trim();
        let expected = self
            .state
            .sequences
            .get(expected_seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Expected sequence '{}' not found", expected_seq_id),
            })?;
        let expected_text = expected.get_forward_string().to_ascii_uppercase();
        if expected_text.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "SuggestSequencingPrimers cannot scan empty expected sequence '{}'",
                    expected_seq_id
                ),
            });
        }
        let confirmation_report = if let Some(report_id) = confirmation_report_id {
            let report = self.get_sequencing_confirmation_report(report_id)?;
            if !report.expected_seq_id.eq_ignore_ascii_case(expected_seq_id) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Sequencing-confirmation report '{}' belongs to expected sequence '{}' rather than '{}'",
                        report.report_id, report.expected_seq_id, expected_seq_id
                    ),
                });
            }
            Some(report)
        } else {
            None
        };
        let template = expected_text.as_bytes();
        let mut suggestions = Vec::new();
        for primer_seq_id in primer_seq_ids {
            let primer_seq_id = primer_seq_id.trim();
            let primer = self
                .state
                .sequences
                .get(primer_seq_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Primer sequence '{}' not found", primer_seq_id),
                })?;
            let primer_sequence = primer.get_forward_string().to_ascii_uppercase();
            if primer_sequence.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Primer sequence '{}' is empty and cannot be scanned",
                        primer_seq_id
                    ),
                });
            }
            if primer_sequence.len() < min_3prime_anneal_bp {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Primer sequence '{}' is shorter than min_3prime_anneal_bp={}",
                        primer_seq_id, min_3prime_anneal_bp
                    ),
                });
            }
            let primer_label = primer
                .name()
                .clone()
                .filter(|value| !value.trim().is_empty())
                .unwrap_or_else(|| primer_seq_id.to_string());
            let anneal_sequence =
                primer_sequence[primer_sequence.len() - min_3prime_anneal_bp..].to_string();
            let forward_hits = Self::find_anneal_sites(
                template,
                anneal_sequence.as_bytes(),
                0,
                min_3prime_anneal_bp,
                true,
            );
            for anneal_start_0based in forward_hits {
                let anneal_end_0based_exclusive = anneal_start_0based + min_3prime_anneal_bp;
                let predicted_read_span_start_0based = anneal_start_0based;
                let predicted_read_span_end_0based_exclusive =
                    (anneal_start_0based + predicted_read_length_bp).min(template.len());
                let mut covered_target_ids = Vec::new();
                let mut covered_problem_target_ids = Vec::new();
                let mut covered_variant_ids = Vec::new();
                let mut covered_problem_variant_ids = Vec::new();
                if let Some(report) = confirmation_report.as_ref() {
                    for target in &report.targets {
                        if Self::sequencing_primer_spans_overlap(
                            predicted_read_span_start_0based,
                            predicted_read_span_end_0based_exclusive,
                            target.start_0based,
                            target.end_0based_exclusive,
                        ) {
                            covered_target_ids.push(target.target_id.clone());
                            if target.status != SequencingConfirmationStatus::Confirmed {
                                covered_problem_target_ids.push(target.target_id.clone());
                            }
                        }
                    }
                    for variant in &report.variants {
                        if Self::sequencing_primer_spans_overlap(
                            predicted_read_span_start_0based,
                            predicted_read_span_end_0based_exclusive,
                            variant.start_0based,
                            variant.end_0based_exclusive,
                        ) {
                            covered_variant_ids.push(variant.variant_id.clone());
                            if !matches!(
                                variant.classification,
                                SequencingConfirmationVariantClassification::ExpectedMatch
                                    | SequencingConfirmationVariantClassification::IntendedEditConfirmed
                            ) {
                                covered_problem_variant_ids.push(variant.variant_id.clone());
                            }
                        }
                    }
                }
                suggestions.push(SequencingPrimerOverlaySuggestion {
                    primer_seq_id: primer_seq_id.to_string(),
                    primer_label: primer_label.clone(),
                    primer_sequence: primer_sequence.clone(),
                    orientation: SequencingPrimerOrientation::ForwardRead,
                    anneal_sequence: anneal_sequence.clone(),
                    anneal_start_0based,
                    anneal_end_0based_exclusive,
                    three_prime_position_0based: anneal_end_0based_exclusive.saturating_sub(1),
                    predicted_read_span_start_0based,
                    predicted_read_span_end_0based_exclusive,
                    covered_target_ids,
                    covered_problem_target_ids,
                    covered_variant_ids,
                    covered_problem_variant_ids,
                });
            }
            let reverse_binding = Self::reverse_complement(&anneal_sequence);
            let reverse_hits = Self::find_anneal_sites(
                template,
                reverse_binding.as_bytes(),
                0,
                min_3prime_anneal_bp,
                false,
            );
            for anneal_start_0based in reverse_hits {
                let anneal_end_0based_exclusive = anneal_start_0based + min_3prime_anneal_bp;
                let predicted_read_span_end_0based_exclusive = anneal_end_0based_exclusive;
                let predicted_read_span_start_0based =
                    anneal_end_0based_exclusive.saturating_sub(predicted_read_length_bp);
                let mut covered_target_ids = Vec::new();
                let mut covered_problem_target_ids = Vec::new();
                let mut covered_variant_ids = Vec::new();
                let mut covered_problem_variant_ids = Vec::new();
                if let Some(report) = confirmation_report.as_ref() {
                    for target in &report.targets {
                        if Self::sequencing_primer_spans_overlap(
                            predicted_read_span_start_0based,
                            predicted_read_span_end_0based_exclusive,
                            target.start_0based,
                            target.end_0based_exclusive,
                        ) {
                            covered_target_ids.push(target.target_id.clone());
                            if target.status != SequencingConfirmationStatus::Confirmed {
                                covered_problem_target_ids.push(target.target_id.clone());
                            }
                        }
                    }
                    for variant in &report.variants {
                        if Self::sequencing_primer_spans_overlap(
                            predicted_read_span_start_0based,
                            predicted_read_span_end_0based_exclusive,
                            variant.start_0based,
                            variant.end_0based_exclusive,
                        ) {
                            covered_variant_ids.push(variant.variant_id.clone());
                            if !matches!(
                                variant.classification,
                                SequencingConfirmationVariantClassification::ExpectedMatch
                                    | SequencingConfirmationVariantClassification::IntendedEditConfirmed
                            ) {
                                covered_problem_variant_ids.push(variant.variant_id.clone());
                            }
                        }
                    }
                }
                suggestions.push(SequencingPrimerOverlaySuggestion {
                    primer_seq_id: primer_seq_id.to_string(),
                    primer_label: primer_label.clone(),
                    primer_sequence: primer_sequence.clone(),
                    orientation: SequencingPrimerOrientation::ReverseRead,
                    anneal_sequence: anneal_sequence.clone(),
                    anneal_start_0based,
                    anneal_end_0based_exclusive,
                    three_prime_position_0based: anneal_start_0based,
                    predicted_read_span_start_0based,
                    predicted_read_span_end_0based_exclusive,
                    covered_target_ids,
                    covered_problem_target_ids,
                    covered_variant_ids,
                    covered_problem_variant_ids,
                });
            }
        }
        suggestions.sort_by(|left, right| {
            left.primer_seq_id
                .to_ascii_lowercase()
                .cmp(&right.primer_seq_id.to_ascii_lowercase())
                .then(left.orientation.as_str().cmp(right.orientation.as_str()))
                .then(left.anneal_start_0based.cmp(&right.anneal_start_0based))
                .then(
                    left.predicted_read_span_start_0based
                        .cmp(&right.predicted_read_span_start_0based),
                )
        });
        let mut problem_guidance = Vec::new();
        let mut proposals = Vec::new();
        if let Some(report) = confirmation_report.as_ref() {
            for target in &report.targets {
                if target.status == SequencingConfirmationStatus::Confirmed {
                    continue;
                }
                let locus_center_0based = Self::sequencing_primer_problem_center(
                    target.start_0based,
                    target.end_0based_exclusive,
                );
                let mut covering_candidates = suggestions
                    .iter()
                    .filter(|row| {
                        row.covered_problem_target_ids
                            .iter()
                            .any(|id| id.eq_ignore_ascii_case(&target.target_id))
                    })
                    .collect::<Vec<_>>();
                let candidate_count = covering_candidates.len();
                let recommendation = if covering_candidates.is_empty() {
                    SequencingPrimerProblemGuidanceRow {
                        problem_id: target.target_id.clone(),
                        problem_kind: SequencingPrimerProblemKind::Target,
                        problem_label: target.label.clone(),
                        problem_summary: target.status.as_str().to_string(),
                        recommended_primer_seq_id: None,
                        recommended_primer_label: None,
                        recommended_orientation: None,
                        recommended_read_span_start_0based: None,
                        recommended_read_span_end_0based_exclusive: None,
                        recommended_three_prime_distance_bp: None,
                        recommended_in_read_direction: None,
                        recommended_problem_target_count: 0,
                        recommended_problem_variant_count: 0,
                        candidate_count,
                        reason:
                            "No existing primer overlay covered this unresolved confirmation target."
                                .to_string(),
                    }
                } else {
                    covering_candidates.sort_by(|left, right| {
                        let (left_in_direction, left_distance) =
                            Self::sequencing_primer_three_prime_distance_bp(
                                left,
                                locus_center_0based,
                            );
                        let (right_in_direction, right_distance) =
                            Self::sequencing_primer_three_prime_distance_bp(
                                right,
                                locus_center_0based,
                            );
                        Self::sequencing_primer_guidance_penalty(left_in_direction, left_distance)
                            .cmp(&Self::sequencing_primer_guidance_penalty(
                                right_in_direction,
                                right_distance,
                            ))
                            .then(
                                (right.covered_problem_target_ids.len()
                                    + right.covered_problem_variant_ids.len())
                                .cmp(
                                    &(left.covered_problem_target_ids.len()
                                        + left.covered_problem_variant_ids.len()),
                                ),
                            )
                            .then(left_distance.cmp(&right_distance))
                            .then(
                                left.primer_seq_id
                                    .to_ascii_lowercase()
                                    .cmp(&right.primer_seq_id.to_ascii_lowercase()),
                            )
                            .then(left.orientation.as_str().cmp(right.orientation.as_str()))
                            .then(left.anneal_start_0based.cmp(&right.anneal_start_0based))
                    });
                    let recommended = covering_candidates[0];
                    let (in_read_direction, three_prime_distance_bp) =
                        Self::sequencing_primer_three_prime_distance_bp(
                            recommended,
                            locus_center_0based,
                        );
                    SequencingPrimerProblemGuidanceRow {
                        problem_id: target.target_id.clone(),
                        problem_kind: SequencingPrimerProblemKind::Target,
                        problem_label: target.label.clone(),
                        problem_summary: target.status.as_str().to_string(),
                        recommended_primer_seq_id: Some(recommended.primer_seq_id.clone()),
                        recommended_primer_label: Some(recommended.primer_label.clone()),
                        recommended_orientation: Some(recommended.orientation),
                        recommended_read_span_start_0based: Some(
                            recommended.predicted_read_span_start_0based,
                        ),
                        recommended_read_span_end_0based_exclusive: Some(
                            recommended.predicted_read_span_end_0based_exclusive,
                        ),
                        recommended_three_prime_distance_bp: Some(three_prime_distance_bp),
                        recommended_in_read_direction: Some(in_read_direction),
                        recommended_problem_target_count: recommended
                            .covered_problem_target_ids
                            .len(),
                        recommended_problem_variant_count: recommended
                            .covered_problem_variant_ids
                            .len(),
                        candidate_count,
                        reason: Self::sequencing_primer_guidance_reason(
                            candidate_count,
                            in_read_direction,
                            three_prime_distance_bp,
                            recommended.covered_problem_target_ids.len(),
                            recommended.covered_problem_variant_ids.len(),
                        ),
                    }
                };
                let needs_proposal = !Self::sequencing_primer_guidance_is_good_enough(
                    &recommendation,
                    predicted_read_length_bp,
                );
                if needs_proposal {
                    let proposal_candidates = Self::build_sequencing_primer_proposal_candidates(
                        &expected_text,
                        target.start_0based,
                        target.end_0based_exclusive,
                        predicted_read_length_bp,
                    );
                    if let Some(candidate) = proposal_candidates.first() {
                        proposals.push(SequencingPrimerProposalRow {
                            proposal_id: format!("{}_proposal", target.target_id),
                            problem_id: target.target_id.clone(),
                            problem_kind: SequencingPrimerProblemKind::Target,
                            problem_label: target.label.clone(),
                            problem_summary: target.status.as_str().to_string(),
                            orientation: candidate.orientation,
                            primer_sequence: candidate.primer.sequence.clone(),
                            anneal_sequence: candidate.anneal_sequence.clone(),
                            anneal_start_0based: candidate.primer.start_0based,
                            anneal_end_0based_exclusive: candidate.primer.end_0based_exclusive,
                            three_prime_position_0based: candidate.three_prime_position_0based,
                            predicted_read_span_start_0based: candidate
                                .predicted_read_span_start_0based,
                            predicted_read_span_end_0based_exclusive: candidate
                                .predicted_read_span_end_0based_exclusive,
                            three_prime_distance_bp: candidate.three_prime_distance_bp,
                            tm_c: candidate.primer.tm_c,
                            gc_fraction: candidate.primer.gc_fraction,
                            anneal_hits: candidate.primer.anneal_hits,
                            three_prime_gc_clamp: candidate.primer.three_prime_gc_clamp,
                            longest_homopolymer_run_bp: candidate.primer.longest_homopolymer_run_bp,
                            self_complementary_run_bp: candidate.primer.self_complementary_run_bp,
                            reason: Self::sequencing_primer_proposal_reason(
                                target.status.as_str(),
                                recommendation.recommended_primer_seq_id.as_deref(),
                                proposal_candidates.len(),
                                candidate.three_prime_distance_bp,
                            ),
                        });
                    }
                }
                problem_guidance.push(recommendation);
            }
            for variant in &report.variants {
                if matches!(
                    variant.classification,
                    SequencingConfirmationVariantClassification::ExpectedMatch
                        | SequencingConfirmationVariantClassification::IntendedEditConfirmed
                ) {
                    continue;
                }
                let locus_center_0based = Self::sequencing_primer_problem_center(
                    variant.start_0based,
                    variant.end_0based_exclusive,
                );
                let mut covering_candidates = suggestions
                    .iter()
                    .filter(|row| {
                        row.covered_problem_variant_ids
                            .iter()
                            .any(|id| id.eq_ignore_ascii_case(&variant.variant_id))
                    })
                    .collect::<Vec<_>>();
                let candidate_count = covering_candidates.len();
                let recommendation = if covering_candidates.is_empty() {
                    SequencingPrimerProblemGuidanceRow {
                        problem_id: variant.variant_id.clone(),
                        problem_kind: SequencingPrimerProblemKind::Variant,
                        problem_label: variant.label.clone(),
                        problem_summary: variant.classification.as_str().to_string(),
                        recommended_primer_seq_id: None,
                        recommended_primer_label: None,
                        recommended_orientation: None,
                        recommended_read_span_start_0based: None,
                        recommended_read_span_end_0based_exclusive: None,
                        recommended_three_prime_distance_bp: None,
                        recommended_in_read_direction: None,
                        recommended_problem_target_count: 0,
                        recommended_problem_variant_count: 0,
                        candidate_count,
                        reason: "No existing primer overlay covered this unresolved variant locus."
                            .to_string(),
                    }
                } else {
                    covering_candidates.sort_by(|left, right| {
                        let (left_in_direction, left_distance) =
                            Self::sequencing_primer_three_prime_distance_bp(
                                left,
                                locus_center_0based,
                            );
                        let (right_in_direction, right_distance) =
                            Self::sequencing_primer_three_prime_distance_bp(
                                right,
                                locus_center_0based,
                            );
                        Self::sequencing_primer_guidance_penalty(left_in_direction, left_distance)
                            .cmp(&Self::sequencing_primer_guidance_penalty(
                                right_in_direction,
                                right_distance,
                            ))
                            .then(
                                (right.covered_problem_target_ids.len()
                                    + right.covered_problem_variant_ids.len())
                                .cmp(
                                    &(left.covered_problem_target_ids.len()
                                        + left.covered_problem_variant_ids.len()),
                                ),
                            )
                            .then(left_distance.cmp(&right_distance))
                            .then(
                                left.primer_seq_id
                                    .to_ascii_lowercase()
                                    .cmp(&right.primer_seq_id.to_ascii_lowercase()),
                            )
                            .then(left.orientation.as_str().cmp(right.orientation.as_str()))
                            .then(left.anneal_start_0based.cmp(&right.anneal_start_0based))
                    });
                    let recommended = covering_candidates[0];
                    let (in_read_direction, three_prime_distance_bp) =
                        Self::sequencing_primer_three_prime_distance_bp(
                            recommended,
                            locus_center_0based,
                        );
                    SequencingPrimerProblemGuidanceRow {
                        problem_id: variant.variant_id.clone(),
                        problem_kind: SequencingPrimerProblemKind::Variant,
                        problem_label: variant.label.clone(),
                        problem_summary: variant.classification.as_str().to_string(),
                        recommended_primer_seq_id: Some(recommended.primer_seq_id.clone()),
                        recommended_primer_label: Some(recommended.primer_label.clone()),
                        recommended_orientation: Some(recommended.orientation),
                        recommended_read_span_start_0based: Some(
                            recommended.predicted_read_span_start_0based,
                        ),
                        recommended_read_span_end_0based_exclusive: Some(
                            recommended.predicted_read_span_end_0based_exclusive,
                        ),
                        recommended_three_prime_distance_bp: Some(three_prime_distance_bp),
                        recommended_in_read_direction: Some(in_read_direction),
                        recommended_problem_target_count: recommended
                            .covered_problem_target_ids
                            .len(),
                        recommended_problem_variant_count: recommended
                            .covered_problem_variant_ids
                            .len(),
                        candidate_count,
                        reason: Self::sequencing_primer_guidance_reason(
                            candidate_count,
                            in_read_direction,
                            three_prime_distance_bp,
                            recommended.covered_problem_target_ids.len(),
                            recommended.covered_problem_variant_ids.len(),
                        ),
                    }
                };
                let needs_proposal = !Self::sequencing_primer_guidance_is_good_enough(
                    &recommendation,
                    predicted_read_length_bp,
                );
                if needs_proposal {
                    let proposal_candidates = Self::build_sequencing_primer_proposal_candidates(
                        &expected_text,
                        variant.start_0based,
                        variant.end_0based_exclusive,
                        predicted_read_length_bp,
                    );
                    if let Some(candidate) = proposal_candidates.first() {
                        proposals.push(SequencingPrimerProposalRow {
                            proposal_id: format!("{}_proposal", variant.variant_id),
                            problem_id: variant.variant_id.clone(),
                            problem_kind: SequencingPrimerProblemKind::Variant,
                            problem_label: variant.label.clone(),
                            problem_summary: variant.classification.as_str().to_string(),
                            orientation: candidate.orientation,
                            primer_sequence: candidate.primer.sequence.clone(),
                            anneal_sequence: candidate.anneal_sequence.clone(),
                            anneal_start_0based: candidate.primer.start_0based,
                            anneal_end_0based_exclusive: candidate.primer.end_0based_exclusive,
                            three_prime_position_0based: candidate.three_prime_position_0based,
                            predicted_read_span_start_0based: candidate
                                .predicted_read_span_start_0based,
                            predicted_read_span_end_0based_exclusive: candidate
                                .predicted_read_span_end_0based_exclusive,
                            three_prime_distance_bp: candidate.three_prime_distance_bp,
                            tm_c: candidate.primer.tm_c,
                            gc_fraction: candidate.primer.gc_fraction,
                            anneal_hits: candidate.primer.anneal_hits,
                            three_prime_gc_clamp: candidate.primer.three_prime_gc_clamp,
                            longest_homopolymer_run_bp: candidate.primer.longest_homopolymer_run_bp,
                            self_complementary_run_bp: candidate.primer.self_complementary_run_bp,
                            reason: Self::sequencing_primer_proposal_reason(
                                variant.classification.as_str(),
                                recommendation.recommended_primer_seq_id.as_deref(),
                                proposal_candidates.len(),
                                candidate.three_prime_distance_bp,
                            ),
                        });
                    }
                }
                problem_guidance.push(recommendation);
            }
        }
        Ok(SequencingPrimerOverlayReport {
            schema: SEQUENCING_PRIMER_OVERLAY_REPORT_SCHEMA.to_string(),
            expected_seq_id: expected_seq_id.to_string(),
            confirmation_report_id: confirmation_report
                .as_ref()
                .map(|row| row.report_id.clone()),
            min_3prime_anneal_bp,
            predicted_read_length_bp,
            primer_seq_ids: primer_seq_ids
                .iter()
                .map(|row| row.trim().to_string())
                .collect(),
            suggestion_count: suggestions.len(),
            suggestions,
            problem_guidance_count: problem_guidance.len(),
            problem_guidance,
            proposal_count: proposals.len(),
            proposals,
        })
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

    fn summarize_variant_confidence(values: &[u8]) -> VariantConfidenceSummary {
        if values.is_empty() {
            return VariantConfidenceSummary::default();
        }
        let min = values.iter().copied().min();
        let max = values.iter().copied().max();
        let sum: usize = values.iter().map(|value| *value as usize).sum();
        VariantConfidenceSummary {
            min,
            max,
            mean: Some(sum as f64 / values.len() as f64),
            count: values.len(),
        }
    }

    fn variant_is_low_confidence_or_ambiguous(
        observed_bases: &str,
        confidence: &VariantConfidenceSummary,
    ) -> bool {
        if observed_bases
            .chars()
            .any(|base| !matches!(base, 'A' | 'C' | 'G' | 'T'))
        {
            return true;
        }
        confidence
            .min
            .is_some_and(|value| value < Self::sequencing_trace_variant_low_confidence_threshold())
    }

    fn extract_variant_observed_allele(
        alignment: &ComputedPairwiseAlignment,
        variant: &VariantLocus,
    ) -> (bool, String, Vec<usize>, Option<usize>) {
        let mut query_pos = alignment
            .report
            .aligned_query_start_0based
            .saturating_sub(alignment.report.query_span_start_0based);
        let mut target_pos = alignment
            .report
            .aligned_target_start_0based
            .saturating_sub(alignment.report.target_span_start_0based);
        let mut observed = String::new();
        let mut query_indices = Vec::new();
        let mut covered_positions = 0usize;
        let mut zero_length_anchor_query_index = None;
        let variant_len = variant
            .end_0based_exclusive
            .saturating_sub(variant.start_0based);
        for op in &alignment.operations {
            match op {
                bio::alignment::AlignmentOperation::Match
                | bio::alignment::AlignmentOperation::Subst => {
                    let target_abs = alignment.report.target_span_start_0based + target_pos;
                    if target_abs >= variant.start_0based
                        && target_abs < variant.end_0based_exclusive
                    {
                        covered_positions += 1;
                        let query_base = alignment
                            .query_span_bases
                            .get(query_pos)
                            .copied()
                            .unwrap_or(b'N') as char;
                        observed.push(query_base);
                        query_indices.push(query_pos);
                    }
                    if variant_len == 0
                        && zero_length_anchor_query_index.is_none()
                        && target_abs >= variant.start_0based
                    {
                        zero_length_anchor_query_index = Some(query_pos);
                    }
                    query_pos += 1;
                    target_pos += 1;
                }
                bio::alignment::AlignmentOperation::Del => {
                    let target_abs = alignment.report.target_span_start_0based + target_pos;
                    if target_abs >= variant.start_0based
                        && target_abs < variant.end_0based_exclusive
                    {
                        covered_positions += 1;
                    }
                    if variant_len == 0
                        && zero_length_anchor_query_index.is_none()
                        && target_abs >= variant.start_0based
                    {
                        zero_length_anchor_query_index = Some(query_pos);
                    }
                    target_pos += 1;
                }
                bio::alignment::AlignmentOperation::Ins => {
                    let target_abs = alignment.report.target_span_start_0based + target_pos;
                    let in_window = if variant_len == 0 {
                        target_abs == variant.start_0based
                    } else {
                        target_abs >= variant.start_0based
                            && target_abs <= variant.end_0based_exclusive
                    };
                    if in_window {
                        let query_base = alignment
                            .query_span_bases
                            .get(query_pos)
                            .copied()
                            .unwrap_or(b'N') as char;
                        observed.push(query_base);
                        query_indices.push(query_pos);
                    }
                    if variant_len == 0 && zero_length_anchor_query_index.is_none() {
                        zero_length_anchor_query_index = Some(query_pos);
                    }
                    query_pos += 1;
                }
                bio::alignment::AlignmentOperation::Xclip(_) => {}
                bio::alignment::AlignmentOperation::Yclip(_) => {}
            }
        }
        let covered = if variant_len == 0 {
            alignment.report.aligned_target_start_0based <= variant.start_0based
                && alignment.report.aligned_target_end_0based_exclusive >= variant.start_0based
        } else {
            covered_positions >= variant_len
        };
        (
            covered,
            observed,
            query_indices,
            zero_length_anchor_query_index,
        )
    }

    fn classify_variant_observation(
        alignment: &ComputedPairwiseAlignment,
        evidence: &SequencingConfirmationEvidenceInput,
        variant: &VariantLocus,
    ) -> VariantObservation {
        let (covered, observed_bases, query_indices, zero_length_anchor_query_index) =
            Self::extract_variant_observed_allele(alignment, variant);
        if !covered {
            return VariantObservation {
                classification: SequencingConfirmationVariantClassification::InsufficientEvidence,
                status: SequencingConfirmationStatus::InsufficientEvidence,
                observed_bases: String::new(),
                confidence: VariantConfidenceSummary::default(),
                peak_center: None,
                reason: "Variant locus is not fully covered by this evidence row".to_string(),
            };
        }
        let confidence_values = evidence
            .confidence_values
            .as_ref()
            .map(|values| {
                query_indices
                    .iter()
                    .filter_map(|idx| values.get(*idx).copied())
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let confidence = Self::summarize_variant_confidence(&confidence_values);
        let peak_center = evidence.peak_locations.as_ref().and_then(|peaks| {
            let positions = if !query_indices.is_empty() {
                query_indices
                    .iter()
                    .filter_map(|idx| peaks.get(*idx).copied())
                    .collect::<Vec<_>>()
            } else if let Some(anchor_idx) = zero_length_anchor_query_index {
                peaks
                    .get(anchor_idx)
                    .copied()
                    .into_iter()
                    .collect::<Vec<_>>()
            } else {
                Vec::new()
            };
            if positions.is_empty() {
                None
            } else {
                Some(positions.iter().sum::<u32>() / positions.len() as u32)
            }
        });
        if Self::variant_is_low_confidence_or_ambiguous(&observed_bases, &confidence) {
            return VariantObservation {
                classification:
                    SequencingConfirmationVariantClassification::LowConfidenceOrAmbiguous,
                status: SequencingConfirmationStatus::InsufficientEvidence,
                observed_bases,
                confidence,
                peak_center,
                reason: format!(
                    "Observed locus is low-confidence or ambiguous (threshold < {})",
                    Self::sequencing_trace_variant_low_confidence_threshold()
                ),
            };
        }
        let expected_bases = variant.expected_bases.to_ascii_uppercase();
        let baseline_bases = variant
            .baseline_bases
            .as_ref()
            .map(|value| value.to_ascii_uppercase());
        let classification = if observed_bases == expected_bases {
            if baseline_bases
                .as_deref()
                .is_some_and(|value| value != expected_bases)
            {
                SequencingConfirmationVariantClassification::IntendedEditConfirmed
            } else {
                SequencingConfirmationVariantClassification::ExpectedMatch
            }
        } else if baseline_bases
            .as_deref()
            .is_some_and(|value| value == observed_bases)
        {
            SequencingConfirmationVariantClassification::ReferenceReversion
        } else {
            SequencingConfirmationVariantClassification::UnexpectedDifference
        };
        let status = match classification {
            SequencingConfirmationVariantClassification::ExpectedMatch
            | SequencingConfirmationVariantClassification::IntendedEditConfirmed => {
                SequencingConfirmationStatus::Confirmed
            }
            SequencingConfirmationVariantClassification::ReferenceReversion
            | SequencingConfirmationVariantClassification::UnexpectedDifference => {
                SequencingConfirmationStatus::Contradicted
            }
            SequencingConfirmationVariantClassification::LowConfidenceOrAmbiguous
            | SequencingConfirmationVariantClassification::InsufficientEvidence => {
                SequencingConfirmationStatus::InsufficientEvidence
            }
        };
        let reason = match classification {
            SequencingConfirmationVariantClassification::ExpectedMatch => {
                "Observed allele matches the expected construct".to_string()
            }
            SequencingConfirmationVariantClassification::IntendedEditConfirmed => format!(
                "Observed allele matches the intended edit relative to baseline '{}'",
                baseline_bases.unwrap_or_default()
            ),
            SequencingConfirmationVariantClassification::ReferenceReversion => format!(
                "Observed allele matches the baseline/reference allele '{}'",
                baseline_bases.unwrap_or_default()
            ),
            SequencingConfirmationVariantClassification::UnexpectedDifference => {
                "Observed allele matches neither expected construct nor baseline".to_string()
            }
            SequencingConfirmationVariantClassification::LowConfidenceOrAmbiguous => {
                "Observed locus is low-confidence or ambiguous".to_string()
            }
            SequencingConfirmationVariantClassification::InsufficientEvidence => {
                "Variant locus lacks enough evidence".to_string()
            }
        };
        VariantObservation {
            classification,
            status,
            observed_bases,
            confidence,
            peak_center,
            reason,
        }
    }

    fn choose_best_variant_observation<'a>(
        observations: &'a [(SequencingConfirmationReadResult, VariantObservation)],
    ) -> Option<&'a (SequencingConfirmationReadResult, VariantObservation)> {
        fn rank(observation: &VariantObservation) -> (i32, usize, i32) {
            let priority = match observation.classification {
                SequencingConfirmationVariantClassification::ReferenceReversion
                | SequencingConfirmationVariantClassification::UnexpectedDifference => 3,
                SequencingConfirmationVariantClassification::IntendedEditConfirmed
                | SequencingConfirmationVariantClassification::ExpectedMatch => 2,
                SequencingConfirmationVariantClassification::LowConfidenceOrAmbiguous => 1,
                SequencingConfirmationVariantClassification::InsufficientEvidence => 0,
            };
            let mean_scaled = (observation.confidence.mean.unwrap_or(0.0) * 1000.0) as i32;
            (priority, observation.confidence.count, mean_scaled)
        }
        observations
            .iter()
            .max_by_key(|(_, observation)| rank(observation))
    }

    pub(super) fn confirm_construct_reads(
        &mut self,
        expected_seq_id: &str,
        baseline_seq_id: Option<&str>,
        read_seq_ids: &[String],
        trace_ids: &[String],
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
        if read_seq_ids.is_empty() && trace_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "ConfirmConstructReads requires at least one read sequence or sequencing trace"
                        .to_string(),
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
        let normalized_baseline_seq_id = baseline_seq_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let mut requested_targets = targets.to_vec();
        if let Some(baseline_seq_id) = normalized_baseline_seq_id.as_deref() {
            let baseline_dna =
                self.state
                    .sequences
                    .get(baseline_seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Baseline sequence '{}' not found", baseline_seq_id),
                    })?;
            let baseline_text = baseline_dna.get_forward_string().to_ascii_uppercase();
            if baseline_text.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "ConfirmConstructReads cannot use empty baseline sequence '{}'",
                        baseline_seq_id
                    ),
                });
            }
            requested_targets.extend(Self::infer_expected_edit_targets(
                expected_seq_id,
                expected_text.as_str(),
                baseline_seq_id,
                baseline_text.as_str(),
            )?);
        }
        let normalized_targets = Self::normalize_sequencing_confirmation_targets(
            expected_text.len(),
            &requested_targets,
        )?;
        let variant_loci = Self::collect_variant_loci(expected_text.as_str(), &normalized_targets);
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
                expected_bases: target.expected_bases.clone(),
                baseline_bases: target.baseline_bases.clone(),
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

        let mut evidence_inputs = vec![];
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
            evidence_inputs.push(SequencingConfirmationEvidenceInput {
                evidence_kind: SequencingConfirmationEvidenceKind::Sequence,
                evidence_id: read_seq_id.clone(),
                display_read_seq_id: read_seq_id.clone(),
                trace_id: None,
                linked_seq_id: Some(read_seq_id.clone()),
                called_bases: read_text,
                confidence_values: None,
                peak_locations: None,
            });
        }
        for trace_id in trace_ids {
            let trace = self.get_sequencing_trace(trace_id)?;
            let called_bases = trace.called_bases.trim().to_ascii_uppercase();
            if called_bases.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "ConfirmConstructReads cannot use sequencing trace '{}' because it has no called bases",
                        trace.trace_id
                    ),
                });
            }
            evidence_inputs.push(SequencingConfirmationEvidenceInput {
                evidence_kind: SequencingConfirmationEvidenceKind::Trace,
                evidence_id: trace.trace_id.clone(),
                display_read_seq_id: trace
                    .seq_id
                    .clone()
                    .unwrap_or_else(|| trace.trace_id.clone()),
                trace_id: Some(trace.trace_id.clone()),
                linked_seq_id: trace.seq_id.clone(),
                called_bases,
                confidence_values: Some(trace.called_base_confidence_values.clone()),
                peak_locations: Some(trace.peak_locations.clone()),
            });
        }

        let mut variant_observations: BTreeMap<
            String,
            Vec<(SequencingConfirmationReadResult, VariantObservation)>,
        > = BTreeMap::new();
        for evidence in evidence_inputs {
            let forward_alignment = Self::compute_pairwise_alignment_report(
                evidence.evidence_id.as_str(),
                evidence.called_bases.as_str(),
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
                let reverse_text = Self::reverse_complement(&evidence.called_bases);
                Some(Self::compute_pairwise_alignment_report(
                    evidence.evidence_id.as_str(),
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
            let (orientation, alignment, oriented_confidence_values, oriented_peak_locations) =
                if let Some(reverse_alignment) = reverse_alignment {
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
                        let reversed_confidences = evidence
                            .confidence_values
                            .as_ref()
                            .map(|values| values.iter().rev().copied().collect::<Vec<_>>());
                        let reversed_peaks = evidence
                            .peak_locations
                            .as_ref()
                            .map(|values| values.iter().rev().copied().collect::<Vec<_>>());
                        (
                            SequencingReadOrientation::ReverseComplement,
                            reverse_alignment,
                            reversed_confidences,
                            reversed_peaks,
                        )
                    } else {
                        (
                            SequencingReadOrientation::Forward,
                            forward_alignment,
                            evidence.confidence_values.clone(),
                            evidence.peak_locations.clone(),
                        )
                    }
                } else {
                    (
                        SequencingReadOrientation::Forward,
                        forward_alignment,
                        evidence.confidence_values.clone(),
                        evidence.peak_locations.clone(),
                    )
                };
            let discrepancies = Self::merge_sequencing_confirmation_discrepancies(
                &Self::extract_sequencing_confirmation_discrepancies(&alignment),
            );
            let mut covered_target_ids = vec![];
            let mut confirmed_target_ids = vec![];
            let mut contradicted_target_ids = vec![];
            for (target_result, target_spec) in
                target_results.iter_mut().zip(normalized_targets.iter())
            {
                if target_spec.kind == SequencingConfirmationTargetKind::ExpectedEdit {
                    continue;
                }
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
                        .push(evidence.evidence_id.clone());
                    contradicted_target_ids.push(target_spec.target_id.clone());
                } else if coverage_fraction >= min_target_coverage_fraction
                    && alignment.report.identity_fraction >= min_identity_fraction
                {
                    target_result
                        .support_read_ids
                        .push(evidence.evidence_id.clone());
                    confirmed_target_ids.push(target_spec.target_id.clone());
                }
            }
            if alignment.report.identity_fraction < min_identity_fraction {
                warnings.push(format!(
                    "{} '{}' best alignment identity {:.3} is below min_identity_fraction {:.3}",
                    evidence.evidence_kind.as_str(),
                    evidence.evidence_id,
                    alignment.report.identity_fraction,
                    min_identity_fraction
                ));
            }
            let mut read_result = SequencingConfirmationReadResult {
                evidence_kind: evidence.evidence_kind,
                evidence_id: evidence.evidence_id.clone(),
                read_seq_id: evidence.display_read_seq_id,
                trace_id: evidence.trace_id,
                linked_seq_id: evidence.linked_seq_id,
                orientation,
                usable: !covered_target_ids.is_empty(),
                best_alignment: alignment.report.clone(),
                covered_target_ids,
                confirmed_target_ids,
                contradicted_target_ids,
                discrepancies,
            };
            let mut variant_contributed = false;
            for variant in &variant_loci {
                let observed = Self::classify_variant_observation(
                    &alignment,
                    &SequencingConfirmationEvidenceInput {
                        evidence_kind: evidence.evidence_kind,
                        evidence_id: evidence.evidence_id.clone(),
                        display_read_seq_id: read_result.read_seq_id.clone(),
                        trace_id: read_result.trace_id.clone(),
                        linked_seq_id: read_result.linked_seq_id.clone(),
                        called_bases: evidence.called_bases.clone(),
                        confidence_values: oriented_confidence_values.clone(),
                        peak_locations: oriented_peak_locations.clone(),
                    },
                    variant,
                );
                if observed.classification
                    != SequencingConfirmationVariantClassification::InsufficientEvidence
                {
                    variant_contributed = true;
                }
                variant_observations
                    .entry(variant.variant_id.clone())
                    .or_default()
                    .push((read_result.clone(), observed));
            }
            if variant_contributed {
                read_result.usable = true;
            }
            reads.push(read_result);
        }

        for target in &mut target_results {
            if target.kind == SequencingConfirmationTargetKind::ExpectedEdit {
                continue;
            }
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

        let mut variant_rows = Vec::with_capacity(variant_loci.len());
        for variant in &variant_loci {
            let observations = variant_observations
                .remove(&variant.variant_id)
                .unwrap_or_default();
            let supporting_ids = observations
                .iter()
                .filter(|(_, observation)| {
                    observation.status == SequencingConfirmationStatus::Confirmed
                })
                .map(|(read, _)| read.evidence_id.clone())
                .collect::<Vec<_>>();
            let contradicting_ids = observations
                .iter()
                .filter(|(_, observation)| {
                    observation.status == SequencingConfirmationStatus::Contradicted
                })
                .map(|(read, _)| read.evidence_id.clone())
                .collect::<Vec<_>>();
            let best = Self::choose_best_variant_observation(&observations);
            let (status, row) = if let Some((read, observation)) = best {
                let status = if !contradicting_ids.is_empty() {
                    SequencingConfirmationStatus::Contradicted
                } else if !supporting_ids.is_empty() {
                    SequencingConfirmationStatus::Confirmed
                } else {
                    SequencingConfirmationStatus::InsufficientEvidence
                };
                (
                    status,
                    SequencingConfirmationVariantRow {
                        variant_id: variant.variant_id.clone(),
                        label: variant.label.clone(),
                        target_id: variant.target_id.clone(),
                        start_0based: variant.start_0based,
                        end_0based_exclusive: variant.end_0based_exclusive,
                        expected_bases: variant.expected_bases.clone(),
                        baseline_bases: variant.baseline_bases.clone(),
                        observed_bases: observation.observed_bases.clone(),
                        classification: observation.classification,
                        status,
                        evidence_kind: read.evidence_kind,
                        evidence_id: read.evidence_id.clone(),
                        trace_id: read.trace_id.clone(),
                        read_seq_id: read.read_seq_id.clone(),
                        linked_seq_id: read.linked_seq_id.clone(),
                        confidence_min: observation.confidence.min,
                        confidence_max: observation.confidence.max,
                        confidence_mean: observation.confidence.mean,
                        confidence_count: observation.confidence.count,
                        peak_center: observation.peak_center,
                        reason: observation.reason.clone(),
                    },
                )
            } else {
                (
                    SequencingConfirmationStatus::InsufficientEvidence,
                    SequencingConfirmationVariantRow {
                        variant_id: variant.variant_id.clone(),
                        label: variant.label.clone(),
                        target_id: variant.target_id.clone(),
                        start_0based: variant.start_0based,
                        end_0based_exclusive: variant.end_0based_exclusive,
                        expected_bases: variant.expected_bases.clone(),
                        baseline_bases: variant.baseline_bases.clone(),
                        observed_bases: String::new(),
                        classification:
                            SequencingConfirmationVariantClassification::InsufficientEvidence,
                        status: SequencingConfirmationStatus::InsufficientEvidence,
                        evidence_kind: SequencingConfirmationEvidenceKind::Sequence,
                        evidence_id: String::new(),
                        trace_id: None,
                        read_seq_id: String::new(),
                        linked_seq_id: None,
                        confidence_min: None,
                        confidence_max: None,
                        confidence_mean: None,
                        confidence_count: 0,
                        peak_center: None,
                        reason: "No evidence row covered this expected-edit locus".to_string(),
                    },
                )
            };
            if let Some(target_id) = variant.target_id.as_deref()
                && let Some(target) = target_results
                    .iter_mut()
                    .find(|target| target.target_id == target_id)
            {
                target.expected_bases = Some(variant.expected_bases.clone());
                target.baseline_bases = variant.baseline_bases.clone();
                target.support_read_ids = supporting_ids.clone();
                target.contradicting_read_ids = contradicting_ids.clone();
                target.status = status;
                target.reason = row.reason.clone();
                target.covered_bp = if status == SequencingConfirmationStatus::InsufficientEvidence
                {
                    0
                } else {
                    target.target_length_bp.max(1)
                };
            }
            variant_rows.push(row);
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
            baseline_seq_id: normalized_baseline_seq_id,
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
            trace_ids: trace_ids.to_vec(),
            target_count: target_results.len(),
            reads,
            targets: target_results,
            variants: variant_rows,
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
            "{} expected={} baseline={} status={} reads={} traces={} evidence={} targets={} variants={} mode={} rc_allowed={} identity>={:.2} coverage>={:.2}",
            report.report_id,
            report.expected_seq_id,
            report.baseline_seq_id.as_deref().unwrap_or("-"),
            report.overall_status.as_str(),
            report.read_seq_ids.len(),
            report.trace_ids.len(),
            report.reads.len(),
            report.targets.len(),
            report.variants.len(),
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
