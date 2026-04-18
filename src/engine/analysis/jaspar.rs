//! Shared JASPAR entry presentation helpers.
//!
//! This module owns the deterministic "show me what this motif likes and what
//! random DNA does" report used by shell/CLI/agent resource inspection paths.
//!
//! Look here for:
//! - per-entry max/min scoring sequence derivation
//! - deterministic pseudorandom background scoring summaries
//! - JSON export helpers for JASPAR presentation reports

use super::*;

impl GentleEngine {
    fn jaspar_extreme_sequence(score_matrix: &[[f64; 4]], pick_max: bool) -> String {
        const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
        let mut out = String::with_capacity(score_matrix.len());
        for column in score_matrix {
            let mut best_idx = 0usize;
            let mut best_score = column[0];
            for (idx, score) in column.iter().copied().enumerate().skip(1) {
                let better = if pick_max {
                    score > best_score
                } else {
                    score < best_score
                };
                if better {
                    best_idx = idx;
                    best_score = score;
                }
            }
            out.push(BASES[best_idx]);
        }
        out
    }

    fn deterministic_random_dna_bytes(length_bp: usize, seed: u64) -> Vec<u8> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let mut state = if seed == 0 {
            DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED
        } else {
            seed
        };
        let mut out = Vec::with_capacity(length_bp);
        for _ in 0..length_bp {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            out.push(BASES[((state >> 62) & 0b11) as usize]);
        }
        out
    }

    fn summarize_jaspar_score_distribution(scores: &[f64]) -> JasparScoreDistributionSummary {
        fn percentile(sorted_scores: &[f64], fraction: f64) -> f64 {
            if sorted_scores.is_empty() {
                return 0.0;
            }
            let capped = fraction.clamp(0.0, 1.0);
            let idx = ((sorted_scores.len().saturating_sub(1) as f64) * capped).round() as usize;
            sorted_scores[idx.min(sorted_scores.len().saturating_sub(1))]
        }

        if scores.is_empty() {
            return JasparScoreDistributionSummary::default();
        }

        let sample_count = scores.len();
        let mut sorted_scores = scores.to_vec();
        sorted_scores
            .sort_by(|left, right| left.partial_cmp(right).unwrap_or(std::cmp::Ordering::Equal));

        let sum = scores.iter().sum::<f64>();
        let mean_score = sum / sample_count as f64;
        let variance = scores
            .iter()
            .map(|score| {
                let delta = *score - mean_score;
                delta * delta
            })
            .sum::<f64>()
            / sample_count as f64;

        JasparScoreDistributionSummary {
            sample_count,
            min_score: *sorted_scores.first().unwrap_or(&0.0),
            max_score: *sorted_scores.last().unwrap_or(&0.0),
            mean_score,
            stddev_score: variance.sqrt(),
            p01_score: percentile(&sorted_scores, 0.01),
            p05_score: percentile(&sorted_scores, 0.05),
            p25_score: percentile(&sorted_scores, 0.25),
            p50_score: percentile(&sorted_scores, 0.50),
            p75_score: percentile(&sorted_scores, 0.75),
            p95_score: percentile(&sorted_scores, 0.95),
            p99_score: percentile(&sorted_scores, 0.99),
            positive_fraction: scores.iter().filter(|score| **score > 0.0).count() as f64
                / sample_count as f64,
        }
    }

    pub(crate) fn summarize_jaspar_entries(
        &self,
        motifs: &[String],
        random_sequence_length_bp: usize,
        random_seed: u64,
    ) -> Result<JasparEntryPresentationReport, EngineError> {
        if random_sequence_length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeJasparEntries requires random_sequence_length_bp >= 1"
                    .to_string(),
            });
        }

        let requested_motifs = if motifs.is_empty() {
            tf_motifs::all_motif_ids()
        } else {
            motifs
                .iter()
                .map(|motif| motif.trim().to_string())
                .filter(|motif| !motif.is_empty())
                .collect::<Vec<_>>()
        };
        if requested_motifs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeJasparEntries requires at least one motif or a non-empty local JASPAR registry".to_string(),
            });
        }

        let random_background =
            Self::deterministic_random_dna_bytes(random_sequence_length_bp, random_seed);

        let mut seen_ids = BTreeSet::new();
        let mut rows = vec![];
        for token in &requested_motifs {
            let motif = tf_motifs::resolve_motif_definition(token).ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "JASPAR entry '{}' was not found in the local motif registry",
                    token
                ),
            })?;
            if !seen_ids.insert(motif.id.clone()) {
                continue;
            }

            let (llr_matrix, true_log_odds_matrix) =
                Self::prepare_scoring_matrices(&motif.matrix_counts);
            if llr_matrix.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::Internal,
                    message: format!("JASPAR entry '{}' has no usable scoring columns", motif.id),
                });
            }
            if random_background.len() < llr_matrix.len() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "random_sequence_length_bp ({}) must be >= motif length ({}) for '{}'",
                        random_background.len(),
                        llr_matrix.len(),
                        motif.id
                    ),
                });
            }

            let maximizing_sequence = Self::jaspar_extreme_sequence(&llr_matrix, true);
            let minimizing_sequence = Self::jaspar_extreme_sequence(&llr_matrix, false);
            let maximizing_bytes = maximizing_sequence.as_bytes();
            let minimizing_bytes = minimizing_sequence.as_bytes();
            let maximizing_llr_bits =
                Self::score_matrix_window(maximizing_bytes, &llr_matrix).unwrap_or(0.0);
            let minimizing_llr_bits =
                Self::score_matrix_window(minimizing_bytes, &llr_matrix).unwrap_or(0.0);
            let maximizing_true_log_odds_bits =
                Self::score_matrix_window(maximizing_bytes, &true_log_odds_matrix).unwrap_or(0.0);
            let minimizing_true_log_odds_bits =
                Self::score_matrix_window(minimizing_bytes, &true_log_odds_matrix).unwrap_or(0.0);

            let hits = Self::scan_tf_scores(
                &random_background,
                &llr_matrix,
                &true_log_odds_matrix,
                |_, _| {},
            );
            let llr_scores = hits.iter().map(|hit| hit.2).collect::<Vec<_>>();
            let true_log_odds_scores = hits.iter().map(|hit| hit.4).collect::<Vec<_>>();
            let mut sorted_llr_scores = llr_scores.clone();
            sorted_llr_scores.sort_by(|left, right| {
                left.partial_cmp(right).unwrap_or(std::cmp::Ordering::Equal)
            });
            let mut sorted_true_log_odds_scores = true_log_odds_scores.clone();
            sorted_true_log_odds_scores.sort_by(|left, right| {
                left.partial_cmp(right).unwrap_or(std::cmp::Ordering::Equal)
            });

            rows.push(JasparEntryPresentationRow {
                motif_id: motif.id,
                motif_name: motif.name,
                consensus_iupac: motif.consensus_iupac,
                motif_length_bp: llr_matrix.len(),
                maximizing_sequence,
                minimizing_sequence,
                maximizing_llr_bits,
                maximizing_llr_quantile: Self::empirical_quantile(
                    &sorted_llr_scores,
                    maximizing_llr_bits,
                ),
                minimizing_llr_bits,
                minimizing_llr_quantile: Self::empirical_quantile(
                    &sorted_llr_scores,
                    minimizing_llr_bits,
                ),
                maximizing_true_log_odds_bits,
                maximizing_true_log_odds_quantile: Self::empirical_quantile(
                    &sorted_true_log_odds_scores,
                    maximizing_true_log_odds_bits,
                ),
                minimizing_true_log_odds_bits,
                minimizing_true_log_odds_quantile: Self::empirical_quantile(
                    &sorted_true_log_odds_scores,
                    minimizing_true_log_odds_bits,
                ),
                llr_bits_distribution: Self::summarize_jaspar_score_distribution(&llr_scores),
                true_log_odds_bits_distribution: Self::summarize_jaspar_score_distribution(
                    &true_log_odds_scores,
                ),
            });
        }

        rows.sort_by(|left, right| {
            left.motif_id
                .to_ascii_uppercase()
                .cmp(&right.motif_id.to_ascii_uppercase())
                .then_with(|| {
                    left.motif_name
                        .as_deref()
                        .unwrap_or("")
                        .to_ascii_uppercase()
                        .cmp(
                            &right
                                .motif_name
                                .as_deref()
                                .unwrap_or("")
                                .to_ascii_uppercase(),
                        )
                })
        });

        Ok(JasparEntryPresentationReport {
            schema: JASPAR_ENTRY_PRESENTATION_REPORT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            requested_motifs,
            registry_entry_count: tf_motifs::all_motif_ids().len(),
            resolved_entry_count: rows.len(),
            random_sequence_length_bp,
            random_seed,
            background_model: "uniform_acgt_lcg".to_string(),
            rows,
        })
    }

    pub(crate) fn write_jaspar_entry_presentation_report_json(
        &self,
        report: &JasparEntryPresentationReport,
        path: &str,
    ) -> Result<(), EngineError> {
        self.write_pretty_json_file(report, path, "JASPAR entry presentation report")
    }
}
