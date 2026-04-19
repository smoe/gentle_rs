//! Shared motif-statistics helpers for TFBS/JASPAR-style scoring.
//!
//! This module owns the generic matrix math that sits below the current
//! adapter-facing motif products:
//! - JASPAR entry presentation / statistics
//! - TFBS score tracks and TFBS feature expert scoring columns
//!
//! The intent is to keep this as the reusable "motif math" layer without
//! collapsing higher-level products together. In particular, future ATtRACT
//! PWM-backed scoring can reuse these helpers while still returning a distinct
//! splice-aware RBP evidence payload.

use super::*;

#[derive(Debug, Clone)]
pub(super) struct ModeledTfbsScoreDistribution {
    pub quantum_bits: f64,
    pub theoretical_min_score: f64,
    pub theoretical_max_score: f64,
    cumulative_bins: Vec<(i32, f64, f64)>,
}

impl ModeledTfbsScoreDistribution {
    fn cumulative_probability_at_or_below_score(&self, score: f64) -> f64 {
        if self.cumulative_bins.is_empty() {
            return 0.0;
        }
        let threshold = score + 0.5 * self.quantum_bits;
        let partition_idx = self
            .cumulative_bins
            .partition_point(|(bin, _, _)| (*bin as f64 * self.quantum_bits) <= threshold);
        if partition_idx == 0 {
            0.0
        } else {
            self.cumulative_bins[partition_idx - 1].2.clamp(0.0, 1.0)
        }
    }

    fn cumulative_probability_below_score(&self, score: f64) -> f64 {
        if self.cumulative_bins.is_empty() {
            return 0.0;
        }
        let threshold = score - 0.5 * self.quantum_bits;
        let partition_idx = self
            .cumulative_bins
            .partition_point(|(bin, _, _)| (*bin as f64 * self.quantum_bits) < threshold);
        if partition_idx == 0 {
            0.0
        } else {
            self.cumulative_bins[partition_idx - 1].2.clamp(0.0, 1.0)
        }
    }

    pub fn modeled_quantile(&self, score: f64) -> f64 {
        let below = self.cumulative_probability_below_score(score);
        let at_or_below = self.cumulative_probability_at_or_below_score(score);
        (below + 0.5 * (at_or_below - below)).clamp(0.0, 1.0)
    }

    pub fn modeled_tail_probability(&self, score: f64) -> f64 {
        (1.0 - self.cumulative_probability_below_score(score)).clamp(0.0, 1.0)
    }

    pub fn modeled_tail_log10(&self, score: f64) -> f64 {
        let tail = self.modeled_tail_probability(score).max(1e-300);
        -tail.log10()
    }

    pub fn score_at_quantile(&self, quantile: f64) -> f64 {
        if self.cumulative_bins.is_empty() {
            return 0.0;
        }
        let target = quantile.clamp(0.0, 1.0);
        let idx = self
            .cumulative_bins
            .partition_point(|(_, _, cumulative_probability)| *cumulative_probability < target)
            .min(self.cumulative_bins.len().saturating_sub(1));
        self.cumulative_bins[idx].0 as f64 * self.quantum_bits
    }
}

impl GentleEngine {
    pub(super) const TFBS_MODELED_SCORE_QUANTUM_BITS: f64 = 1e-3;

    pub(super) fn smooth_probability_matrix(matrix_counts: &[[f64; 4]]) -> Vec<[f64; 4]> {
        if matrix_counts.is_empty() {
            return vec![];
        }
        let max_col_sum = matrix_counts
            .iter()
            .map(|c| c.iter().sum::<f64>())
            .fold(0.0_f64, f64::max);
        let baseline = max_col_sum.max(1.0);

        let mut out = Vec::with_capacity(matrix_counts.len());
        for col in matrix_counts {
            let mut adjusted = *col;
            let col_sum = adjusted.iter().sum::<f64>();

            // Missing observations are distributed uniformly to match the
            // highest-supported column count in this motif.
            if baseline > col_sum {
                let add = (baseline - col_sum) / 4.0;
                for v in &mut adjusted {
                    *v += add;
                }
            }

            let epsilon = baseline * 1e-9;
            for v in &mut adjusted {
                *v += epsilon;
            }

            let total = adjusted.iter().sum::<f64>();
            let mut p_col = [0.0_f64; 4];
            for i in 0..4 {
                let p = (adjusted[i] / total).clamp(f64::MIN_POSITIVE, 1.0 - f64::EPSILON);
                p_col[i] = p;
            }
            out.push(p_col);
        }
        out
    }

    pub(super) fn prepare_scoring_matrices(
        matrix_counts: &[[f64; 4]],
    ) -> (Vec<[f64; 4]>, Vec<[f64; 4]>) {
        let probabilities = Self::smooth_probability_matrix(matrix_counts);
        let background = [0.25_f64, 0.25_f64, 0.25_f64, 0.25_f64];
        let mut llr = Vec::with_capacity(probabilities.len());
        let mut true_log_odds = Vec::with_capacity(probabilities.len());

        for col in probabilities {
            let mut llr_col = [0.0_f64; 4];
            let mut lor_col = [0.0_f64; 4];
            for i in 0..4 {
                let p = col[i];
                let q = background[i];
                llr_col[i] = (p / q).log2();
                let odds_p = p / (1.0 - p);
                let odds_q = q / (1.0 - q);
                lor_col[i] = (odds_p / odds_q).log2();
            }
            llr.push(llr_col);
            true_log_odds.push(lor_col);
        }
        (llr, true_log_odds)
    }

    pub(super) fn base_to_idx(base: u8) -> Option<usize> {
        match base.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        }
    }

    pub(super) fn score_matrix_window(window: &[u8], score_matrix: &[[f64; 4]]) -> Option<f64> {
        if window.len() != score_matrix.len() {
            return None;
        }
        let mut score = 0.0_f64;
        for (idx, base) in window.iter().enumerate() {
            let b = Self::base_to_idx(*base)?;
            score += score_matrix[idx][b];
        }
        Some(score)
    }

    pub(super) fn empirical_quantile(sorted_scores: &[f64], score: f64) -> f64 {
        if sorted_scores.is_empty() {
            return 0.0;
        }
        let mut lo = 0usize;
        let mut hi = sorted_scores.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            if sorted_scores[mid] <= score {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo as f64 / sorted_scores.len() as f64
    }

    pub(super) fn motif_score_theoretical_bounds(score_matrix: &[[f64; 4]]) -> Option<(f64, f64)> {
        if score_matrix.is_empty() {
            return None;
        }
        let theoretical_min_score = score_matrix
            .iter()
            .map(|column| column.iter().copied().fold(f64::INFINITY, f64::min))
            .sum::<f64>();
        let theoretical_max_score = score_matrix
            .iter()
            .map(|column| column.iter().copied().fold(f64::NEG_INFINITY, f64::max))
            .sum::<f64>();
        Some((theoretical_min_score, theoretical_max_score))
    }

    pub(super) fn modeled_tfbs_score_distribution(
        score_matrix: &[[f64; 4]],
    ) -> Option<ModeledTfbsScoreDistribution> {
        if score_matrix.is_empty() {
            return None;
        }
        let (theoretical_min_score, theoretical_max_score) =
            Self::motif_score_theoretical_bounds(score_matrix)?;
        let quantum_bits = Self::TFBS_MODELED_SCORE_QUANTUM_BITS;
        let quantize = |score: f64| -> i32 { (score / quantum_bits).round() as i32 };
        let mut support = std::collections::HashMap::<i32, f64>::from([(0_i32, 1.0_f64)]);
        for column in score_matrix {
            let quantized_column = column.map(quantize);
            let mut next = std::collections::HashMap::<i32, f64>::with_capacity(
                support.len().saturating_mul(quantized_column.len()),
            );
            for (partial_score, probability) in &support {
                for quantized_score in quantized_column {
                    *next
                        .entry(partial_score.saturating_add(quantized_score))
                        .or_insert(0.0) += probability * 0.25;
                }
            }
            support = next;
        }
        let mut bins = support.into_iter().collect::<Vec<_>>();
        bins.sort_by_key(|(score_bin, _)| *score_bin);
        let mut cumulative_probability = 0.0_f64;
        let cumulative_bins = bins
            .into_iter()
            .map(|(score_bin, probability)| {
                cumulative_probability += probability;
                (score_bin, probability, cumulative_probability)
            })
            .collect::<Vec<_>>();
        Some(ModeledTfbsScoreDistribution {
            quantum_bits,
            theoretical_min_score,
            theoretical_max_score,
            cumulative_bins,
        })
    }

    pub(super) fn scan_tf_scores(
        sequence: &[u8],
        llr_matrix: &[[f64; 4]],
        true_log_odds_matrix: &[[f64; 4]],
        mut on_progress: impl FnMut(usize, usize),
    ) -> Vec<(usize, bool, f64, f64, f64, f64)> {
        Self::scan_tf_scores_with_topology(
            sequence,
            llr_matrix,
            true_log_odds_matrix,
            InlineSequenceTopology::Linear,
            |scanned_steps, total_steps| on_progress(scanned_steps, total_steps),
        )
    }

    pub(super) fn scan_tf_scores_with_topology(
        sequence: &[u8],
        llr_matrix: &[[f64; 4]],
        true_log_odds_matrix: &[[f64; 4]],
        topology: InlineSequenceTopology,
        mut on_progress: impl FnMut(usize, usize),
    ) -> Vec<(usize, bool, f64, f64, f64, f64)> {
        if llr_matrix.is_empty()
            || sequence.len() < llr_matrix.len()
            || llr_matrix.len() != true_log_odds_matrix.len()
        {
            return vec![];
        }
        let mut raw_hits = Vec::new();
        let mut all_llr_scores = Vec::new();
        let mut all_true_log_odds_scores = Vec::new();
        let len = llr_matrix.len();
        let circular_windows = matches!(topology, InlineSequenceTopology::Circular);
        let windows = if circular_windows {
            sequence.len()
        } else {
            sequence.len().saturating_sub(len).saturating_add(1)
        };
        let total_steps = windows.saturating_mul(2);
        let progress_stride = (total_steps / 200).max(1);
        let mut scanned_steps = 0usize;
        on_progress(scanned_steps, total_steps);
        let circular_sequence = if circular_windows && len > 1 {
            let mut bytes = Vec::with_capacity(sequence.len() + len - 1);
            bytes.extend_from_slice(sequence);
            bytes.extend_from_slice(&sequence[..len - 1]);
            Some(bytes)
        } else {
            None
        };
        let source_bytes = circular_sequence.as_deref().unwrap_or(sequence);
        for start in 0..windows {
            let window = &source_bytes[start..start + len];
            if let (Some(llr), Some(true_log_odds)) = (
                Self::score_matrix_window(window, llr_matrix),
                Self::score_matrix_window(window, true_log_odds_matrix),
            ) {
                all_llr_scores.push(llr);
                all_true_log_odds_scores.push(true_log_odds);
                raw_hits.push((start, false, llr, true_log_odds));
            }
            scanned_steps = scanned_steps.saturating_add(1);
            if scanned_steps % progress_stride == 0 || scanned_steps == total_steps {
                on_progress(scanned_steps, total_steps);
            }
            let rc_window = Self::reverse_complement_bytes(window);
            if let (Some(llr), Some(true_log_odds)) = (
                Self::score_matrix_window(&rc_window, llr_matrix),
                Self::score_matrix_window(&rc_window, true_log_odds_matrix),
            ) {
                all_llr_scores.push(llr);
                all_true_log_odds_scores.push(true_log_odds);
                raw_hits.push((start, true, llr, true_log_odds));
            }
            scanned_steps = scanned_steps.saturating_add(1);
            if scanned_steps % progress_stride == 0 || scanned_steps == total_steps {
                on_progress(scanned_steps, total_steps);
            }
        }
        if scanned_steps != total_steps {
            on_progress(total_steps, total_steps);
        }
        all_llr_scores.sort_by(|a, b| a.total_cmp(b));
        all_true_log_odds_scores.sort_by(|a, b| a.total_cmp(b));
        raw_hits
            .into_iter()
            .map(|(start, reverse, llr_bits, true_log_odds_bits)| {
                (
                    start,
                    reverse,
                    llr_bits,
                    Self::empirical_quantile(&all_llr_scores, llr_bits),
                    true_log_odds_bits,
                    Self::empirical_quantile(&all_true_log_odds_scores, true_log_odds_bits),
                )
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modeled_tfbs_score_distribution_tracks_bounds_and_nonzero_tail() {
        let matrix_counts = GentleEngine::matrix_from_iupac("GGGGCGGGG");
        let (llr_matrix, _true_log_odds_matrix) =
            GentleEngine::prepare_scoring_matrices(&matrix_counts);
        let modeled = GentleEngine::modeled_tfbs_score_distribution(&llr_matrix)
            .expect("modeled distribution");
        let maximizing_sequence = b"GGGGCGGGG";
        let maximizing_score =
            GentleEngine::score_matrix_window(maximizing_sequence, &llr_matrix).expect("score");

        assert!(modeled.theoretical_max_score >= maximizing_score);
        assert!(modeled.theoretical_min_score <= maximizing_score);
        assert!(modeled.modeled_quantile(maximizing_score) < 1.0);
        assert!(modeled.modeled_tail_probability(maximizing_score) > 0.0);
        assert!(modeled.modeled_tail_probability(maximizing_score) < 1.0);
        assert!(modeled.score_at_quantile(0.99) <= modeled.theoretical_max_score);
    }
}
