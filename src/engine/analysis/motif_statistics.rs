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
    rounding_error_bits: f64,
    survival_log_probabilities: Vec<f64>,
}

impl ModeledTfbsScoreDistribution {
    fn cumulative_probability_at_or_below_score(&self, score: f64) -> f64 {
        if self.cumulative_bins.is_empty() {
            return 0.0;
        }
        let threshold = score + self.rounding_error_bits;
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
        let threshold = score - self.rounding_error_bits;
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
        self.tail_log_probability(score).exp()
    }

    pub fn modeled_tail_log10(&self, score: f64) -> f64 {
        -self.tail_log_probability(score) / std::f64::consts::LN_10
    }

    fn tail_log_probability(&self, score: f64) -> f64 {
        // If raw S >= observed, rounded Q >= observed - sum(column errors).
        // Including that whole bin gives an upper bound on the inclusive tail,
        // not an overstatement of significance. Never subtract a CDF from one.
        let threshold = score - self.rounding_error_bits;
        let index = self
            .cumulative_bins
            .partition_point(|(bin, _, _)| *bin as f64 * self.quantum_bits < threshold);
        self.survival_log_probabilities[index].min(0.0)
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

fn log_add_probability(a: f64, b: f64) -> f64 {
    let hi = a.max(b);
    if hi == f64::NEG_INFINITY {
        hi
    } else {
        hi + (a.min(b) - hi).exp().ln_1p()
    }
}

impl GentleEngine {
    pub(super) fn tfbs_cancelled_error(context: &str) -> EngineError {
        EngineError {
            code: ErrorCode::Internal,
            message: format!("TFBS scoring cancelled during {context}"),
            cause_chain: vec![],
        }
    }

    pub(super) const TFBS_MODELED_SCORE_QUANTUM_BITS: f64 = 1e-3;
    pub(super) const TFBS_MODELED_TAIL_METHOD: &str =
        "uniform_iid_quantized_conservative_survival_v2";

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
        if score_matrix.is_empty()
            || score_matrix
                .iter()
                .flatten()
                .any(|score| !score.is_finite())
        {
            return None;
        }
        let (theoretical_min_score, theoretical_max_score) =
            Self::motif_score_theoretical_bounds(score_matrix)?;
        let quantum_bits = Self::TFBS_MODELED_SCORE_QUANTUM_BITS;
        let quantize = |score: f64| -> i32 { (score / quantum_bits).round() as i32 };
        let rounding_error_bits = score_matrix
            .iter()
            .map(|column| {
                column
                    .iter()
                    .map(|score| (score - quantize(*score) as f64 * quantum_bits).abs())
                    .fold(0.0_f64, f64::max)
            })
            .sum::<f64>();
        // Also cover floating-point accumulation in raw sums, bounds and Q*q.
        let magnitude = score_matrix
            .iter()
            .map(|column| {
                column
                    .iter()
                    .map(|score| score.abs())
                    .fold(0.0_f64, f64::max)
            })
            .sum::<f64>()
            .max(1.0);
        let rounding_error_bits =
            rounding_error_bits + 8.0 * f64::EPSILON * score_matrix.len() as f64 * magnitude;
        // Short-motif masses are dyadic and safely representable. Long motifs
        // use log masses so even a unique >500-base maximum cannot underflow.
        let logarithmic = score_matrix.len() > 500;
        let mut support = vec![(0_i32, if logarithmic { 0.0 } else { 1.0 })];
        for column in score_matrix {
            let quantized_column = column.map(quantize);
            let mut next = std::collections::HashMap::<i32, f64>::with_capacity(
                support.len().saturating_mul(quantized_column.len()),
            );
            for (partial_score, probability) in &support {
                for quantized_score in quantized_column {
                    let entry = next
                        .entry(partial_score.checked_add(quantized_score)?)
                        .or_insert(if logarithmic { f64::NEG_INFINITY } else { 0.0 });
                    if logarithmic {
                        *entry = log_add_probability(*entry, probability - 4.0_f64.ln());
                    } else {
                        *entry += probability * 0.25;
                    }
                }
            }
            support = next.into_iter().collect();
            support.sort_by_key(|(bin, _)| *bin);
        }
        let mut survival_log_probabilities = vec![f64::NEG_INFINITY; support.len() + 1];
        for (index, (_, mass)) in support.iter().enumerate().rev() {
            survival_log_probabilities[index] = log_add_probability(
                if logarithmic { *mass } else { mass.ln() },
                survival_log_probabilities[index + 1],
            );
        }
        let mut cumulative_probability = 0.0_f64;
        let cumulative_bins = support
            .into_iter()
            .map(|(score_bin, mass)| {
                let probability = if logarithmic { mass.exp() } else { mass };
                cumulative_probability += probability;
                (score_bin, probability, cumulative_probability)
            })
            .collect::<Vec<_>>();
        Some(ModeledTfbsScoreDistribution {
            quantum_bits,
            theoretical_min_score,
            theoretical_max_score,
            cumulative_bins,
            rounding_error_bits,
            survival_log_probabilities,
        })
    }

    pub(super) fn scan_tf_scores(
        sequence: &[u8],
        llr_matrix: &[[f64; 4]],
        true_log_odds_matrix: &[[f64; 4]],
        mut on_progress: impl FnMut(usize, usize),
    ) -> Vec<(usize, bool, f64, f64, f64, f64)> {
        Self::scan_tf_scores_with_topology_and_cancel(
            sequence,
            llr_matrix,
            true_log_odds_matrix,
            InlineSequenceTopology::Linear,
            |scanned_steps, total_steps| {
                on_progress(scanned_steps, total_steps);
                true
            },
        )
        .unwrap_or_default()
    }

    pub(super) fn scan_tf_scores_with_topology(
        sequence: &[u8],
        llr_matrix: &[[f64; 4]],
        true_log_odds_matrix: &[[f64; 4]],
        topology: InlineSequenceTopology,
        mut on_progress: impl FnMut(usize, usize),
    ) -> Vec<(usize, bool, f64, f64, f64, f64)> {
        Self::scan_tf_scores_with_topology_and_cancel(
            sequence,
            llr_matrix,
            true_log_odds_matrix,
            topology,
            |scanned_steps, total_steps| {
                on_progress(scanned_steps, total_steps);
                true
            },
        )
        .unwrap_or_default()
    }

    pub(super) fn scan_tf_scores_with_topology_and_cancel(
        sequence: &[u8],
        llr_matrix: &[[f64; 4]],
        true_log_odds_matrix: &[[f64; 4]],
        topology: InlineSequenceTopology,
        mut on_progress: impl FnMut(usize, usize) -> bool,
    ) -> Result<Vec<(usize, bool, f64, f64, f64, f64)>, EngineError> {
        if llr_matrix.is_empty()
            || sequence.len() < llr_matrix.len()
            || llr_matrix.len() != true_log_odds_matrix.len()
        {
            return Ok(vec![]);
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
        if !on_progress(scanned_steps, total_steps) {
            return Err(Self::tfbs_cancelled_error("scan setup"));
        }
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
            if (scanned_steps.is_multiple_of(progress_stride) || scanned_steps == total_steps)
                && !on_progress(scanned_steps, total_steps)
            {
                return Err(Self::tfbs_cancelled_error("forward-strand scan"));
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
            if (scanned_steps.is_multiple_of(progress_stride) || scanned_steps == total_steps)
                && !on_progress(scanned_steps, total_steps)
            {
                return Err(Self::tfbs_cancelled_error("reverse-strand scan"));
            }
        }
        if scanned_steps != total_steps && !on_progress(total_steps, total_steps) {
            return Err(Self::tfbs_cancelled_error("scan completion"));
        }
        all_llr_scores.sort_by(|a, b| a.total_cmp(b));
        all_true_log_odds_scores.sort_by(|a, b| a.total_cmp(b));
        Ok(raw_hits
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
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Hand-crafted score matrices, not empirical PFMs. Exhaustively recreate
    // their uniform-background oracle by summing all 4^L possible words.
    fn enumerate_scores(matrix: &[[f64; 4]]) -> Vec<f64> {
        let mut scores = vec![0.0];
        for column in matrix {
            scores = scores
                .into_iter()
                .flat_map(|score| column.map(|v| score + v))
                .collect();
        }
        scores.sort_by(f64::total_cmp);
        scores
    }

    #[test]
    fn inclusive_survival_matches_exhaustive_short_motifs_and_tied_maxima() {
        for matrix in [
            vec![[0.0, 1.0, 1.0, -1.0]; 4],
            vec![[0.00049, 0.80049, 0.90049, 1.00049]; 3],
            vec![[0.0; 4]; 5],
        ] {
            let scores = enumerate_scores(&matrix);
            let model = GentleEngine::modeled_tfbs_score_distribution(&matrix).unwrap();
            let mut previous = 1.0_f64;
            for &score in &scores {
                // Include arithmetic ties to within the same floating-point guard.
                let exact = scores.iter().filter(|s| **s >= score - 1e-12).count() as f64
                    / scores.len() as f64;
                let tail = model.modeled_tail_probability(score);
                assert!((tail - exact).abs() < 1e-12, "{score}: {tail} != {exact}");
                assert!(tail <= previous + 1e-15);
                previous = tail;
                assert!(
                    model.modeled_tail_log10(score)
                        <= matrix.len() as f64 * 4.0_f64.log10() + 1e-12
                );
            }
        }
    }

    #[test]
    fn quantized_survival_is_conservative_when_distinct_raw_scores_share_bins() {
        let matrix = vec![[0.0001, 0.0002, 0.0003, 0.0004]; 4];
        let scores = enumerate_scores(&matrix);
        let model = GentleEngine::modeled_tfbs_score_distribution(&matrix).unwrap();
        for &score in &scores {
            let exact = scores.iter().filter(|s| **s >= score).count() as f64 / scores.len() as f64;
            assert!(model.modeled_tail_probability(score) + 1e-14 >= exact);
        }
        assert_eq!(model.modeled_tail_log10(*scores.last().unwrap()), 0.0);
    }

    #[test]
    fn long_motif_rounding_and_tiny_survival_do_not_inflate_significance() {
        for length in [16, 80, 500, 600] {
            let matrix = vec![[1.00049, 0.0, 0.0, 0.0]; length];
            let score = GentleEngine::score_matrix_window(&vec![b'A'; length], &matrix).unwrap();
            let model = GentleEngine::modeled_tfbs_score_distribution(&matrix).unwrap();
            let expected = length as f64 * 4.0_f64.log10();
            assert!((model.modeled_tail_log10(score) - expected).abs() < 1e-9);
            if length <= 500 {
                let expected_tail = 4.0_f64.powi(-(length as i32));
                assert!((model.modeled_tail_probability(score) / expected_tail - 1.0).abs() < 1e-9);
            }
            assert_eq!(
                GentleEngine::score_matrix_window(&vec![b'N'; length], &matrix),
                None
            );
            assert_eq!(GentleEngine::score_matrix_window(b"R", &matrix), None);
        }
    }

    #[test]
    fn tp73_ma0861_2_maximum_has_inclusive_four_to_minus_sixteen_tail() {
        // Exact bundled JASPAR PFM; no active/user registry or network lookup.
        let db = crate::tf_motifs::TfMotifDb::from_json_for_test(include_str!(
            "../../../assets/jaspar.motifs.json"
        ))
        .unwrap();
        let motif = db.resolve_exact_full_pfm("MA0861.2").unwrap();
        let (matrix, _) = GentleEngine::prepare_scoring_matrices(&motif.matrix_counts);
        let score = GentleEngine::score_matrix_window(b"ACATGTCTGGACATGT", &matrix).unwrap();
        assert!((score - 19.543680326691).abs() < 1e-11);
        let model = GentleEngine::modeled_tfbs_score_distribution(&matrix).unwrap();
        assert!((model.modeled_tail_probability(score) / 4.0_f64.powi(-16) - 1.0).abs() < 1e-12);
        assert!((model.modeled_tail_log10(score) - 9.632959861247).abs() < 1e-11);
    }

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
