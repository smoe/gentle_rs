//! Shared promoter-design helpers that extend beyond one specific variant.
//!
//! This module currently owns continuous TF motif score tracks used by the
//! GUI-side Promoter design expert and by headless JSON export paths.

use super::motif_statistics::ModeledTfbsScoreDistribution;
use super::*;
use crate::feature_location::{feature_is_reverse, feature_ranges_sorted_i64};

impl GentleEngine {
    const TFBS_BACKGROUND_TAIL_SHOW_QUANTILE: f64 = 0.95;
    const TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP: usize = 25;

    fn tfbs_background_tail_log10(tail_probability: f64, modeled_quantile: f64) -> f64 {
        if modeled_quantile < Self::TFBS_BACKGROUND_TAIL_SHOW_QUANTILE {
            return 0.0;
        }
        -tail_probability.max(1e-300).log10()
    }

    fn tfbs_modeled_background_quantile(
        score: f64,
        modeled_distribution: Option<&ModeledTfbsScoreDistribution>,
        empirical_background_sorted_scores: &[f64],
    ) -> f64 {
        modeled_distribution
            .map(|distribution| distribution.modeled_quantile(score))
            .unwrap_or_else(|| Self::empirical_quantile(empirical_background_sorted_scores, score))
    }

    fn tfbs_modeled_background_tail_probability(
        score: f64,
        modeled_distribution: Option<&ModeledTfbsScoreDistribution>,
        empirical_background_sorted_scores: &[f64],
    ) -> f64 {
        modeled_distribution
            .map(|distribution| distribution.modeled_tail_probability(score))
            .unwrap_or_else(|| {
                let quantile = Self::empirical_quantile(empirical_background_sorted_scores, score);
                (1.0 - quantile)
                    .max(1.0 / empirical_background_sorted_scores.len().max(1) as f64)
                    .clamp(0.0, 1.0)
            })
    }

    fn tfbs_score_track_presented_value(
        llr_bits: f64,
        llr_quantile: f64,
        true_log_odds_bits: f64,
        true_log_odds_quantile: f64,
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        llr_modeled_distribution: Option<&ModeledTfbsScoreDistribution>,
        true_log_odds_modeled_distribution: Option<&ModeledTfbsScoreDistribution>,
        llr_background_sorted_scores: &[f64],
        true_log_odds_background_sorted_scores: &[f64],
    ) -> f64 {
        match score_kind {
            TfbsScoreTrackValueKind::LlrBits => {
                Self::promoter_design_clip_score(llr_bits, clip_negative)
            }
            TfbsScoreTrackValueKind::LlrQuantile => llr_quantile,
            TfbsScoreTrackValueKind::LlrBackgroundQuantile => {
                let quantile = Self::tfbs_modeled_background_quantile(
                    llr_bits,
                    llr_modeled_distribution,
                    llr_background_sorted_scores,
                );
                if quantile >= Self::TFBS_BACKGROUND_TAIL_SHOW_QUANTILE {
                    quantile
                } else {
                    0.0
                }
            }
            TfbsScoreTrackValueKind::LlrBackgroundTailLog10 => {
                let quantile = Self::tfbs_modeled_background_quantile(
                    llr_bits,
                    llr_modeled_distribution,
                    llr_background_sorted_scores,
                );
                let tail_probability = Self::tfbs_modeled_background_tail_probability(
                    llr_bits,
                    llr_modeled_distribution,
                    llr_background_sorted_scores,
                );
                Self::tfbs_background_tail_log10(tail_probability, quantile)
            }
            TfbsScoreTrackValueKind::TrueLogOddsBits => {
                Self::promoter_design_clip_score(true_log_odds_bits, clip_negative)
            }
            TfbsScoreTrackValueKind::TrueLogOddsQuantile => true_log_odds_quantile,
            TfbsScoreTrackValueKind::TrueLogOddsBackgroundQuantile => {
                let quantile = Self::tfbs_modeled_background_quantile(
                    true_log_odds_bits,
                    true_log_odds_modeled_distribution,
                    true_log_odds_background_sorted_scores,
                );
                if quantile >= Self::TFBS_BACKGROUND_TAIL_SHOW_QUANTILE {
                    quantile
                } else {
                    0.0
                }
            }
            TfbsScoreTrackValueKind::TrueLogOddsBackgroundTailLog10 => {
                let quantile = Self::tfbs_modeled_background_quantile(
                    true_log_odds_bits,
                    true_log_odds_modeled_distribution,
                    true_log_odds_background_sorted_scores,
                );
                let tail_probability = Self::tfbs_modeled_background_tail_probability(
                    true_log_odds_bits,
                    true_log_odds_modeled_distribution,
                    true_log_odds_background_sorted_scores,
                );
                Self::tfbs_background_tail_log10(tail_probability, quantile)
            }
        }
    }

    fn promoter_design_clip_score(score: f64, clip_negative: bool) -> f64 {
        if clip_negative && score.is_sign_negative() {
            0.0
        } else if score.is_finite() {
            score
        } else {
            0.0
        }
    }

    fn summarize_tfbs_score_track_tss_markers_from_provenance(
        &self,
        seq_id: &str,
        start_0based: usize,
        end_0based_exclusive: usize,
    ) -> Vec<TfbsScoreTrackTssMarker> {
        let Some(provenance) = self.latest_genome_extraction_provenance_for_seq(seq_id) else {
            return vec![];
        };
        let (Some(interval_start_1based), Some(interval_end_1based), Some(tss_1based)) = (
            provenance.start_1based,
            provenance.end_1based,
            provenance.tss_1based,
        ) else {
            return vec![];
        };
        if tss_1based < interval_start_1based || tss_1based > interval_end_1based {
            return vec![];
        }
        let position_0based = tss_1based.saturating_sub(interval_start_1based);
        if position_0based.saturating_add(1) < start_0based
            || position_0based > end_0based_exclusive
        {
            return vec![];
        }
        let label = provenance
            .transcript_id
            .clone()
            .or(provenance.gene_name.clone())
            .or(provenance.gene_query.clone())
            .unwrap_or_else(|| "genome_promoter_slice".to_string());
        vec![TfbsScoreTrackTssMarker {
            feature_id: usize::MAX,
            feature_kind: "genome_promoter_slice".to_string(),
            label,
            position_0based,
            is_reverse: matches!(provenance.strand, Some('-')),
        }]
    }

    fn summarize_tfbs_score_track_tss_markers(
        &self,
        seq_id: &str,
        dna: &DNAsequence,
        start_0based: usize,
        end_0based_exclusive: usize,
    ) -> Vec<TfbsScoreTrackTssMarker> {
        let mut markers = dna
            .features()
            .iter()
            .enumerate()
            .filter_map(|(feature_id, feature)| {
                let feature_kind = feature.kind.to_string();
                let normalized_kind = feature_kind.to_ascii_uppercase();
                if !matches!(normalized_kind.as_str(), "GENE" | "MRNA" | "PROMOTER") {
                    return None;
                }
                let ranges = feature_ranges_sorted_i64(feature);
                let (first_start, last_end) = match (ranges.first(), ranges.last()) {
                    (Some(first), Some(last)) => (first.0, last.1),
                    _ => return None,
                };
                let is_reverse = feature_is_reverse(feature);
                let tss_i64 = if is_reverse {
                    last_end.saturating_sub(1)
                } else {
                    first_start
                };
                let position_0based = usize::try_from(tss_i64).ok()?;
                if position_0based.saturating_add(1) < start_0based
                    || position_0based > end_0based_exclusive
                {
                    return None;
                }
                let label = Self::first_nonempty_feature_qualifier(
                    feature,
                    &["transcript_id", "label", "name", "standard_name", "gene"],
                )
                .unwrap_or_else(|| Self::feature_display_label(feature, feature_id));
                Some(TfbsScoreTrackTssMarker {
                    feature_id,
                    feature_kind,
                    label,
                    position_0based,
                    is_reverse,
                })
            })
            .collect::<Vec<_>>();
        markers.sort_by(|left, right| {
            left.position_0based
                .cmp(&right.position_0based)
                .then(left.feature_kind.cmp(&right.feature_kind))
                .then(left.label.cmp(&right.label))
        });
        markers.dedup_by(|left, right| {
            left.position_0based == right.position_0based
                && left.is_reverse == right.is_reverse
                && left.label == right.label
        });
        if markers.is_empty() {
            markers.extend(self.summarize_tfbs_score_track_tss_markers_from_provenance(
                seq_id,
                start_0based,
                end_0based_exclusive,
            ));
        }
        markers
    }

    fn tfbs_track_display_signal(track: &TfbsScoreTrackRow) -> Vec<f64> {
        track
            .forward_scores
            .iter()
            .copied()
            .zip(track.reverse_scores.iter().copied())
            .map(|(forward, reverse)| forward.max(reverse))
            .collect()
    }

    fn smooth_tfbs_track_signal(signal: &[f64], window_bp: usize) -> Vec<f64> {
        if signal.is_empty() || window_bp <= 1 {
            return signal.to_vec();
        }
        let radius = window_bp / 2;
        let mut prefix = Vec::with_capacity(signal.len() + 1);
        prefix.push(0.0);
        for value in signal {
            prefix.push(prefix.last().copied().unwrap_or(0.0) + value);
        }
        let mut smoothed = Vec::with_capacity(signal.len());
        for idx in 0..signal.len() {
            let start = idx.saturating_sub(radius);
            let end = (idx + radius + 1).min(signal.len());
            let sum = prefix[end] - prefix[start];
            smoothed.push(sum / (end - start).max(1) as f64);
        }
        smoothed
    }

    fn pearson_correlation(left: &[f64], right: &[f64]) -> f64 {
        let len = left.len().min(right.len());
        if len < 2 {
            return 0.0;
        }
        let left = &left[..len];
        let right = &right[..len];
        let left_mean = left.iter().sum::<f64>() / len as f64;
        let right_mean = right.iter().sum::<f64>() / len as f64;
        let mut covariance = 0.0;
        let mut left_variance = 0.0;
        let mut right_variance = 0.0;
        for (left_value, right_value) in left.iter().zip(right.iter()) {
            let left_centered = *left_value - left_mean;
            let right_centered = *right_value - right_mean;
            covariance += left_centered * right_centered;
            left_variance += left_centered * left_centered;
            right_variance += right_centered * right_centered;
        }
        if left_variance <= f64::EPSILON || right_variance <= f64::EPSILON {
            0.0
        } else {
            (covariance / (left_variance.sqrt() * right_variance.sqrt())).clamp(-1.0, 1.0)
        }
    }

    fn summarize_tfbs_track_primary_peak_offset_bp(
        left: &TfbsScoreTrackRow,
        right: &TfbsScoreTrackRow,
    ) -> Option<i64> {
        let left_peak = left
            .top_peaks
            .first()
            .map(|peak| peak.start_0based)
            .or(left.max_position_0based)?;
        let right_peak = right
            .top_peaks
            .first()
            .map(|peak| peak.start_0based)
            .or(right.max_position_0based)?;
        Some(right_peak as i64 - left_peak as i64)
    }

    fn summarize_tfbs_score_track_correlation_summary(
        tracks: &[TfbsScoreTrackRow],
    ) -> Option<TfbsScoreTrackCorrelationSummary> {
        if tracks.len() < 2 {
            return None;
        }
        let mut rows = vec![];
        for left_idx in 0..tracks.len() {
            let left = &tracks[left_idx];
            let left_signal = Self::tfbs_track_display_signal(left);
            if left_signal.len() < 2 {
                continue;
            }
            let left_smoothed = Self::smooth_tfbs_track_signal(
                &left_signal,
                Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
            );
            for right in tracks.iter().skip(left_idx + 1) {
                let right_signal = Self::tfbs_track_display_signal(right);
                let overlap_window_count = left_signal.len().min(right_signal.len());
                if overlap_window_count < 2 {
                    continue;
                }
                let right_smoothed = Self::smooth_tfbs_track_signal(
                    &right_signal,
                    Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
                );
                rows.push(TfbsScoreTrackCorrelationRow {
                    left_tf_id: left.tf_id.clone(),
                    left_tf_name: left.tf_name.clone(),
                    right_tf_id: right.tf_id.clone(),
                    right_tf_name: right.tf_name.clone(),
                    overlap_window_count,
                    raw_pearson: Self::pearson_correlation(
                        &left_signal[..overlap_window_count],
                        &right_signal[..overlap_window_count],
                    ),
                    smoothed_pearson: Self::pearson_correlation(
                        &left_smoothed[..overlap_window_count],
                        &right_smoothed[..overlap_window_count],
                    ),
                    signed_primary_peak_offset_bp:
                        Self::summarize_tfbs_track_primary_peak_offset_bp(left, right),
                });
            }
        }
        if rows.is_empty() {
            return None;
        }
        rows.sort_by(|left, right| {
            right
                .smoothed_pearson
                .abs()
                .total_cmp(&left.smoothed_pearson.abs())
                .then(right.raw_pearson.abs().total_cmp(&left.raw_pearson.abs()))
                .then(left.left_tf_id.cmp(&right.left_tf_id))
                .then(left.right_tf_id.cmp(&right.right_tf_id))
        });
        Some(TfbsScoreTrackCorrelationSummary {
            signal_source: "max(forward_score, reverse_score)".to_string(),
            smoothing_method: "centered_boxcar".to_string(),
            smoothing_window_bp: Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
            pair_count: rows.len(),
            rows,
        })
    }

    fn collect_tfbs_score_track_background_scores(
        random_background: &[u8],
        llr_matrix: &[[f64; 4]],
        true_log_odds_matrix: &[[f64; 4]],
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        llr_modeled_distribution: Option<&ModeledTfbsScoreDistribution>,
        true_log_odds_modeled_distribution: Option<&ModeledTfbsScoreDistribution>,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let hits = Self::scan_tf_scores(
            random_background,
            llr_matrix,
            true_log_odds_matrix,
            |_, _| {},
        );
        let mut llr_background_scores = hits.iter().map(|row| row.2).collect::<Vec<_>>();
        let mut true_log_odds_background_scores = hits.iter().map(|row| row.4).collect::<Vec<_>>();
        llr_background_scores.sort_by(|left, right| left.total_cmp(right));
        true_log_odds_background_scores.sort_by(|left, right| left.total_cmp(right));
        let displayed_scores = hits
            .into_iter()
            .map(
                |(_, _, llr_bits, llr_quantile, true_log_odds_bits, true_log_odds_quantile)| {
                    Self::tfbs_score_track_presented_value(
                        llr_bits,
                        llr_quantile,
                        true_log_odds_bits,
                        true_log_odds_quantile,
                        score_kind,
                        clip_negative,
                        llr_modeled_distribution,
                        true_log_odds_modeled_distribution,
                        &llr_background_scores,
                        &true_log_odds_background_scores,
                    )
                },
            )
            .collect::<Vec<_>>();
        let underlying_background_scores = if score_kind.uses_llr_background_bits() {
            llr_background_scores.clone()
        } else {
            true_log_odds_background_scores.clone()
        };
        (
            displayed_scores,
            underlying_background_scores,
            llr_background_scores,
            true_log_odds_background_scores,
        )
    }

    fn summarize_tfbs_score_track_normalization_reference(
        underlying_background_scores: &[f64],
        observed_peak_underlying_score: f64,
        modeled_distribution: Option<&ModeledTfbsScoreDistribution>,
    ) -> Option<TfbsScoreTrackNormalizationReference> {
        if underlying_background_scores.is_empty() {
            return None;
        }
        let distribution = Self::summarize_jaspar_score_distribution(underlying_background_scores);
        let mut sorted_scores = underlying_background_scores.to_vec();
        sorted_scores.sort_by(|left, right| left.total_cmp(right));
        let modeled_quantile = modeled_distribution
            .map(|distribution| distribution.modeled_quantile(observed_peak_underlying_score))
            .unwrap_or_else(|| {
                Self::empirical_quantile(&sorted_scores, observed_peak_underlying_score)
            });
        let modeled_tail_probability = modeled_distribution
            .map(|distribution| {
                distribution.modeled_tail_probability(observed_peak_underlying_score)
            })
            .unwrap_or_else(|| {
                (1.0 - modeled_quantile)
                    .max(1.0_f64 / sorted_scores.len().max(1) as f64)
                    .clamp(0.0, 1.0)
            });
        let modeled_tail_log10 = modeled_distribution
            .map(|distribution| distribution.modeled_tail_log10(observed_peak_underlying_score))
            .unwrap_or_else(|| -modeled_tail_probability.max(1e-300).log10());
        let theoretical_min_score = modeled_distribution
            .map(|distribution| distribution.theoretical_min_score)
            .unwrap_or_else(|| distribution.min_score.min(observed_peak_underlying_score));
        let theoretical_max_score = modeled_distribution
            .map(|distribution| distribution.theoretical_max_score)
            .unwrap_or_else(|| distribution.max_score.max(observed_peak_underlying_score));
        Some(TfbsScoreTrackNormalizationReference {
            background_model: "uniform_random_dna".to_string(),
            chance_model: "quantized_iid_uniform_window_dp".to_string(),
            random_sequence_length_bp: DEFAULT_TFBS_SCORE_TRACK_RANDOM_SEQUENCE_LENGTH_BP,
            random_seed: DEFAULT_TFBS_SCORE_TRACK_RANDOM_SEED,
            sample_count: distribution.sample_count,
            mean_score: distribution.mean_score,
            stddev_score: distribution.stddev_score,
            p95_score: modeled_distribution
                .map(|distribution| distribution.score_at_quantile(0.95))
                .unwrap_or(distribution.p95_score),
            p99_score: modeled_distribution
                .map(|distribution| distribution.score_at_quantile(0.99))
                .unwrap_or(distribution.p99_score),
            positive_fraction: distribution.positive_fraction,
            observed_peak_empirical_quantile: Self::empirical_quantile(
                &sorted_scores,
                observed_peak_underlying_score,
            ),
            observed_peak_modeled_quantile: modeled_quantile,
            observed_peak_modeled_tail_probability: modeled_tail_probability,
            observed_peak_modeled_tail_log10: modeled_tail_log10,
            observed_peak_delta_from_p95: observed_peak_underlying_score
                - modeled_distribution
                    .map(|distribution| distribution.score_at_quantile(0.95))
                    .unwrap_or(distribution.p95_score),
            observed_peak_delta_from_p99: observed_peak_underlying_score
                - modeled_distribution
                    .map(|distribution| distribution.score_at_quantile(0.99))
                    .unwrap_or(distribution.p99_score),
            theoretical_min_score,
            theoretical_max_score,
        })
    }

    fn summarize_tfbs_score_track_top_peaks(
        track_start_0based: usize,
        motif_length_bp: usize,
        forward_scores: &[f64],
        reverse_scores: &[f64],
        displayed_background_scores: &[f64],
    ) -> Vec<TfbsScoreTrackPeak> {
        fn is_local_peak(scores: &[f64], idx: usize) -> bool {
            let score = scores[idx];
            let left = idx.checked_sub(1).and_then(|i| scores.get(i)).copied();
            let right = scores.get(idx + 1).copied();
            left.is_none_or(|neighbor| score >= neighbor)
                && right.is_none_or(|neighbor| score >= neighbor)
        }

        #[derive(Clone)]
        struct PeakCandidate {
            start_0based: usize,
            end_0based_exclusive: usize,
            is_reverse: bool,
            score: f64,
            empirical_quantile: f64,
            delta_from_p99: f64,
        }

        let mut candidates = vec![];
        let mut background_sorted_scores = displayed_background_scores.to_vec();
        background_sorted_scores.sort_by(|left, right| left.total_cmp(right));
        let displayed_distribution =
            Self::summarize_jaspar_score_distribution(displayed_background_scores);
        let displayed_p99 = displayed_distribution.p99_score;
        let mut push_candidates = |scores: &[f64], is_reverse: bool| {
            for (idx, score) in scores.iter().copied().enumerate() {
                if !score.is_finite() || !is_local_peak(scores, idx) {
                    continue;
                }
                let start_0based = track_start_0based + idx;
                let empirical_quantile = if background_sorted_scores.is_empty() {
                    0.0
                } else {
                    Self::empirical_quantile(&background_sorted_scores, score)
                };
                let delta_from_p99 = score - displayed_p99;
                candidates.push(PeakCandidate {
                    start_0based,
                    end_0based_exclusive: start_0based + motif_length_bp,
                    is_reverse,
                    score,
                    empirical_quantile,
                    delta_from_p99,
                });
            }
        };
        push_candidates(forward_scores, false);
        push_candidates(reverse_scores, true);

        if candidates.is_empty() {
            return vec![];
        }
        candidates.sort_by(|left, right| {
            right
                .score
                .total_cmp(&left.score)
                .then(left.start_0based.cmp(&right.start_0based))
                .then(left.is_reverse.cmp(&right.is_reverse))
        });

        let threshold = displayed_p99;
        let mut selected = vec![];
        for candidate in &candidates {
            if candidate.score <= threshold {
                continue;
            }
            if selected.iter().any(|picked: &PeakCandidate| {
                picked.start_0based.abs_diff(candidate.start_0based) < motif_length_bp
            }) {
                continue;
            }
            selected.push(candidate.clone());
            if selected.len() == 3 {
                break;
            }
        }
        if selected.is_empty() {
            selected.push(candidates[0].clone());
        }
        selected
            .into_iter()
            .enumerate()
            .map(|(rank_idx, peak)| TfbsScoreTrackPeak {
                rank: rank_idx + 1,
                start_0based: peak.start_0based,
                end_0based_exclusive: peak.end_0based_exclusive,
                is_reverse: peak.is_reverse,
                score: peak.score,
                empirical_quantile: peak.empirical_quantile,
                delta_from_p99: peak.delta_from_p99,
            })
            .collect()
    }

    pub(crate) fn summarize_tfbs_score_tracks(
        &self,
        target: SequenceScanTarget,
        motifs: &[String],
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
    ) -> Result<TfbsScoreTrackReport, EngineError> {
        if motifs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeTfbsScoreTracks requires at least one motif".to_string(),
            });
        }
        let tss_source_seq_id = match &target {
            SequenceScanTarget::SeqId { seq_id, .. } => Some(seq_id.clone()),
            SequenceScanTarget::InlineSequence { .. } => None,
        };
        let (
            target_kind,
            target_label,
            source_sequence_length_bp,
            scan_start_0based,
            scan_end_0based_exclusive,
            scan_topology,
            scan_dna,
        ) = self.resolve_sequence_scan_target(&target, "SummarizeTfbsScoreTracks")?;

        let sequence = scan_dna.get_forward_string();
        let view = sequence.as_bytes();
        let random_background = Self::deterministic_random_dna_bytes(
            DEFAULT_TFBS_SCORE_TRACK_RANDOM_SEQUENCE_LENGTH_BP,
            DEFAULT_TFBS_SCORE_TRACK_RANDOM_SEED,
        );

        let mut tracks = vec![];
        let mut global_max_score = 0.0_f64;
        for motif in motifs {
            let (tf_id, tf_name, _consensus, matrix_counts) =
                Self::resolve_tf_motif_for_scoring(motif)?;
            let (llr_matrix, true_log_odds_matrix) = Self::prepare_scoring_matrices(&matrix_counts);
            let llr_modeled_distribution = Self::modeled_tfbs_score_distribution(&llr_matrix);
            let true_log_odds_modeled_distribution =
                Self::modeled_tfbs_score_distribution(&true_log_odds_matrix);
            let motif_length_bp = llr_matrix.len();
            let scored_window_count = if motif_length_bp == 0 || view.len() < motif_length_bp {
                0
            } else if matches!(scan_topology, InlineSequenceTopology::Circular) {
                view.len()
            } else {
                view.len() - motif_length_bp + 1
            };
            let mut forward_scores = Vec::with_capacity(scored_window_count);
            let mut reverse_scores = Vec::with_capacity(scored_window_count);
            let mut max_score = 0.0_f64;
            let mut max_underlying_score = f64::NEG_INFINITY;
            let mut max_position_0based = None;
            let (
                displayed_background_scores,
                underlying_background_scores,
                llr_background_sorted_scores,
                true_log_odds_background_sorted_scores,
            ) = Self::collect_tfbs_score_track_background_scores(
                &random_background,
                &llr_matrix,
                &true_log_odds_matrix,
                score_kind,
                clip_negative,
                llr_modeled_distribution.as_ref(),
                true_log_odds_modeled_distribution.as_ref(),
            );

            if scored_window_count > 0 {
                let hits = Self::scan_tf_scores_with_topology(
                    view,
                    &llr_matrix,
                    &true_log_odds_matrix,
                    scan_topology,
                    |_, _| {},
                );
                forward_scores.resize(scored_window_count, 0.0);
                reverse_scores.resize(scored_window_count, 0.0);
                for (
                    offset,
                    reverse,
                    llr_bits,
                    llr_quantile,
                    true_log_odds_bits,
                    true_log_odds_quantile,
                ) in hits
                {
                    let score = Self::tfbs_score_track_presented_value(
                        llr_bits,
                        llr_quantile,
                        true_log_odds_bits,
                        true_log_odds_quantile,
                        score_kind,
                        clip_negative,
                        llr_modeled_distribution.as_ref(),
                        true_log_odds_modeled_distribution.as_ref(),
                        &llr_background_sorted_scores,
                        &true_log_odds_background_sorted_scores,
                    );
                    if reverse {
                        reverse_scores[offset] = score;
                    } else {
                        forward_scores[offset] = score;
                    }
                    if score > max_score {
                        max_score = score;
                        max_position_0based = Some(scan_start_0based + offset);
                    }
                    let underlying_score = if score_kind.uses_llr_background_bits() {
                        llr_bits
                    } else {
                        true_log_odds_bits
                    };
                    if underlying_score > max_underlying_score {
                        max_underlying_score = underlying_score;
                    }
                }
            }
            if !max_underlying_score.is_finite() {
                max_underlying_score = 0.0;
            }
            let normalization_reference = Self::summarize_tfbs_score_track_normalization_reference(
                &underlying_background_scores,
                max_underlying_score,
                if score_kind.uses_llr_background_bits() {
                    llr_modeled_distribution.as_ref()
                } else {
                    true_log_odds_modeled_distribution.as_ref()
                },
            );
            let top_peaks = Self::summarize_tfbs_score_track_top_peaks(
                scan_start_0based,
                motif_length_bp,
                &forward_scores,
                &reverse_scores,
                &displayed_background_scores,
            );
            global_max_score = global_max_score.max(max_score);
            tracks.push(TfbsScoreTrackRow {
                tf_id,
                tf_name: tf_name.or_else(|| Some(motif.clone())),
                motif_length_bp,
                track_start_0based: scan_start_0based,
                scored_window_count,
                max_score,
                max_position_0based,
                normalization_reference,
                top_peaks,
                forward_scores,
                reverse_scores,
            });
        }

        let tss_markers = if let Some(seq_id) = tss_source_seq_id {
            self.state
                .sequences
                .get(&seq_id)
                .map(|dna| {
                    Self::summarize_tfbs_score_track_tss_markers(
                        dna,
                        scan_start_0based,
                        scan_end_0based_exclusive,
                    )
                })
                .unwrap_or_default()
        } else {
            vec![]
        };

        Ok(TfbsScoreTrackReport {
            schema: TFBS_SCORE_TRACK_REPORT_SCHEMA.to_string(),
            target_kind,
            seq_id: target_label.clone(),
            target_label,
            source_sequence_length_bp,
            sequence_length_bp: source_sequence_length_bp,
            scan_topology,
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            score_kind,
            view_start_0based: scan_start_0based,
            view_end_0based_exclusive: scan_end_0based_exclusive,
            clip_negative,
            motifs_requested: motifs.to_vec(),
            global_max_score,
            tss_markers,
            correlation_summary: Self::summarize_tfbs_score_track_correlation_summary(&tracks),
            tracks,
        })
    }

    pub(crate) fn scan_tfbs_hits(
        &self,
        target: SequenceScanTarget,
        motifs: &[String],
        min_llr_bits: Option<f64>,
        min_llr_quantile: Option<f64>,
        per_tf_thresholds: &[TfThresholdOverride],
        max_hits: Option<usize>,
        op_id: Option<&str>,
        run_id: Option<&str>,
    ) -> Result<TfbsHitScanReport, EngineError> {
        let effective_motifs = if motifs.len() == 1
            && matches!(motifs[0].trim().to_ascii_uppercase().as_str(), "ALL" | "*")
        {
            crate::tf_motifs::all_motif_ids()
        } else {
            motifs.to_vec()
        };
        if effective_motifs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ScanTfbsHits requires at least one motif".to_string(),
            });
        }

        let default_min_llr_quantile = min_llr_quantile.unwrap_or(0.0);
        Self::validate_tf_thresholds(default_min_llr_quantile)?;

        let mut override_map: std::collections::HashMap<String, (Option<f64>, Option<f64>)> =
            std::collections::HashMap::new();
        for override_row in per_tf_thresholds {
            let key = override_row.tf.trim().to_ascii_uppercase();
            if key.is_empty() {
                continue;
            }
            if let Some(q) = override_row.min_llr_quantile {
                Self::validate_tf_thresholds(q)?;
            }
            override_map.insert(
                key,
                (override_row.min_llr_bits, override_row.min_llr_quantile),
            );
        }

        let (
            target_kind,
            target_label,
            source_sequence_length_bp,
            scan_start_0based,
            scan_end_0based_exclusive,
            scan_topology,
            scan_dna,
        ) = self.resolve_sequence_scan_target(&target, "ScanTfbsHits")?;
        let scan_sequence = scan_dna.get_forward_string();
        let scan_bytes = scan_sequence.as_bytes();

        let effective_max_hits = match max_hits {
            Some(0) => None,
            Some(value) => Some(value),
            None => None,
        };
        let default_min_llr_bits = min_llr_bits.unwrap_or(f64::NEG_INFINITY);
        let mut report = TfbsHitScanReport {
            schema: TFBS_HIT_SCAN_REPORT_SCHEMA.to_string(),
            target_kind,
            target_label,
            source_sequence_length_bp,
            scan_start_0based,
            scan_end_0based_exclusive,
            scan_length_bp: scan_dna.len(),
            scan_topology,
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: op_id.map(str::to_string),
            run_id: run_id.map(str::to_string),
            motifs_requested: motifs.to_vec(),
            motifs_scanned: vec![],
            default_min_llr_bits: min_llr_bits,
            default_min_llr_quantile: min_llr_quantile,
            per_tf_thresholds: per_tf_thresholds.to_vec(),
            max_hits: effective_max_hits,
            truncated_at_max_hits: false,
            matched_hit_count: 0,
            path: None,
            rows: vec![],
        };

        'motif_loop: for token in effective_motifs {
            let token_key = token.trim().to_ascii_uppercase();
            if token_key.is_empty() {
                continue;
            }
            let (tf_id, tf_name, consensus, matrix_counts) =
                Self::resolve_tf_motif_for_scoring(&token)?;
            let (llr_matrix, true_log_odds_matrix) = Self::prepare_scoring_matrices(&matrix_counts);
            if llr_matrix.is_empty() || llr_matrix.len() > scan_bytes.len() {
                continue;
            }

            let mut effective_bits = default_min_llr_bits;
            let mut effective_quantile = default_min_llr_quantile;
            let id_key = tf_id.to_ascii_uppercase();
            let name_key = tf_name
                .as_ref()
                .map(|name| name.trim().to_ascii_uppercase())
                .unwrap_or_default();
            for key in [token_key.as_str(), id_key.as_str(), name_key.as_str()] {
                if key.is_empty() {
                    continue;
                }
                if let Some((bits, quantile)) = override_map.get(key) {
                    if let Some(value) = bits {
                        effective_bits = *value;
                    }
                    if let Some(value) = quantile {
                        effective_quantile = *value;
                    }
                    break;
                }
            }

            report.motifs_scanned.push(tf_id.clone());
            for (
                start,
                reverse,
                llr_bits,
                llr_quantile,
                true_log_odds_bits,
                true_log_odds_quantile,
            ) in Self::scan_tf_scores_with_topology(
                scan_bytes,
                &llr_matrix,
                &true_log_odds_matrix,
                report.scan_topology,
                |_, _| {},
            ) {
                if llr_bits < effective_bits || llr_quantile < effective_quantile {
                    continue;
                }
                let raw_end = start + llr_matrix.len();
                let wraps_origin = raw_end > report.scan_length_bp;
                let end = if wraps_origin {
                    raw_end - report.scan_length_bp
                } else {
                    raw_end
                };
                let matched_sequence = if wraps_origin {
                    format!("{}{}", &scan_sequence[start..], &scan_sequence[..end])
                } else {
                    scan_sequence[start..end].to_string()
                };
                report.rows.push(TfbsHitScanRow {
                    tf_id: tf_id.clone(),
                    tf_name: tf_name.clone(),
                    motif_consensus_iupac: consensus.clone(),
                    motif_length_bp: llr_matrix.len(),
                    match_start_0based: start,
                    match_end_0based_exclusive: end,
                    source_match_start_0based: report.scan_start_0based + start,
                    source_match_end_0based_exclusive: report.scan_start_0based + end,
                    wraps_origin,
                    forward_strand: !reverse,
                    matched_sequence,
                    llr_bits,
                    llr_quantile,
                    true_log_odds_bits,
                    true_log_odds_quantile,
                });
                if let Some(limit) = effective_max_hits
                    && report.rows.len() >= limit
                {
                    report.truncated_at_max_hits = true;
                    break 'motif_loop;
                }
            }
        }

        report.matched_hit_count = report.rows.len();
        Ok(report)
    }

    pub(crate) fn write_tfbs_score_track_report_json(
        &self,
        report: &TfbsScoreTrackReport,
        path: &str,
    ) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize TFBS score-track report '{}' for '{}': {e}",
                report.seq_id, path
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write TFBS score-track report to '{path}': {e}"),
        })
    }
}
