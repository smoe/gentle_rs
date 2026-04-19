//! Shared promoter-design helpers that extend beyond one specific variant.
//!
//! This module currently owns continuous TF motif score tracks used by the
//! GUI-side Promoter design expert and by headless JSON export paths.

use super::*;
use crate::feature_location::{feature_is_reverse, feature_ranges_sorted_i64};

impl GentleEngine {
    const TFBS_BACKGROUND_TAIL_SHOW_QUANTILE: f64 = 0.95;

    fn tfbs_background_tail_log10(quantile: f64, sample_count: usize) -> f64 {
        if quantile < Self::TFBS_BACKGROUND_TAIL_SHOW_QUANTILE {
            return 0.0;
        }
        let min_tail = (1.0 / sample_count.max(1) as f64).max(1e-12);
        let tail = (1.0 - quantile).max(min_tail);
        -tail.log10()
    }

    fn tfbs_score_track_presented_value(
        llr_bits: f64,
        llr_quantile: f64,
        true_log_odds_bits: f64,
        true_log_odds_quantile: f64,
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        llr_background_sorted_scores: &[f64],
        true_log_odds_background_sorted_scores: &[f64],
    ) -> f64 {
        match score_kind {
            TfbsScoreTrackValueKind::LlrBits => {
                Self::promoter_design_clip_score(llr_bits, clip_negative)
            }
            TfbsScoreTrackValueKind::LlrQuantile => llr_quantile,
            TfbsScoreTrackValueKind::LlrBackgroundQuantile => {
                let quantile = Self::empirical_quantile(llr_background_sorted_scores, llr_bits);
                if quantile >= Self::TFBS_BACKGROUND_TAIL_SHOW_QUANTILE {
                    quantile
                } else {
                    0.0
                }
            }
            TfbsScoreTrackValueKind::LlrBackgroundTailLog10 => {
                let quantile = Self::empirical_quantile(llr_background_sorted_scores, llr_bits);
                Self::tfbs_background_tail_log10(quantile, llr_background_sorted_scores.len())
            }
            TfbsScoreTrackValueKind::TrueLogOddsBits => {
                Self::promoter_design_clip_score(true_log_odds_bits, clip_negative)
            }
            TfbsScoreTrackValueKind::TrueLogOddsQuantile => true_log_odds_quantile,
            TfbsScoreTrackValueKind::TrueLogOddsBackgroundQuantile => {
                let quantile = Self::empirical_quantile(
                    true_log_odds_background_sorted_scores,
                    true_log_odds_bits,
                );
                if quantile >= Self::TFBS_BACKGROUND_TAIL_SHOW_QUANTILE {
                    quantile
                } else {
                    0.0
                }
            }
            TfbsScoreTrackValueKind::TrueLogOddsBackgroundTailLog10 => {
                let quantile = Self::empirical_quantile(
                    true_log_odds_background_sorted_scores,
                    true_log_odds_bits,
                );
                Self::tfbs_background_tail_log10(
                    quantile,
                    true_log_odds_background_sorted_scores.len(),
                )
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

    fn summarize_tfbs_score_track_tss_markers(
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
        markers
    }

    fn collect_tfbs_score_track_background_scores(
        random_background: &[u8],
        llr_matrix: &[[f64; 4]],
        true_log_odds_matrix: &[[f64; 4]],
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
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
                        &llr_background_scores,
                        &true_log_odds_background_scores,
                    )
                },
            )
            .collect::<Vec<_>>();
        (
            displayed_scores,
            llr_background_scores,
            true_log_odds_background_scores,
        )
    }

    fn summarize_tfbs_score_track_normalization_reference(
        background_scores: &[f64],
        observed_peak_score: f64,
    ) -> Option<TfbsScoreTrackNormalizationReference> {
        if background_scores.is_empty() {
            return None;
        }
        let distribution = Self::summarize_jaspar_score_distribution(background_scores);
        let mut sorted_scores = background_scores.to_vec();
        sorted_scores.sort_by(|left, right| left.total_cmp(right));
        Some(TfbsScoreTrackNormalizationReference {
            background_model: "uniform_random_dna".to_string(),
            random_sequence_length_bp: DEFAULT_JASPAR_PRESENTATION_RANDOM_SEQUENCE_LENGTH_BP,
            random_seed: DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED,
            sample_count: distribution.sample_count,
            mean_score: distribution.mean_score,
            stddev_score: distribution.stddev_score,
            p95_score: distribution.p95_score,
            p99_score: distribution.p99_score,
            positive_fraction: distribution.positive_fraction,
            observed_peak_empirical_quantile: Self::empirical_quantile(
                &sorted_scores,
                observed_peak_score,
            ),
            observed_peak_delta_from_p95: observed_peak_score - distribution.p95_score,
            observed_peak_delta_from_p99: observed_peak_score - distribution.p99_score,
        })
    }

    fn summarize_tfbs_score_track_top_peaks(
        track_start_0based: usize,
        motif_length_bp: usize,
        forward_scores: &[f64],
        reverse_scores: &[f64],
        background_scores: &[f64],
        normalization_reference: Option<&TfbsScoreTrackNormalizationReference>,
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
        let mut background_sorted_scores = background_scores.to_vec();
        background_sorted_scores.sort_by(|left, right| left.total_cmp(right));
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
                let delta_from_p99 = normalization_reference
                    .map(|normalization| score - normalization.p99_score)
                    .unwrap_or(0.0);
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

        let threshold = normalization_reference
            .map(|normalization| normalization.p99_score)
            .unwrap_or(0.0);
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
        seq_id: &str,
        motifs: &[String],
        start_0based: usize,
        end_0based_exclusive: usize,
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
    ) -> Result<TfbsScoreTrackReport, EngineError> {
        if motifs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeTfbsScoreTracks requires at least one motif".to_string(),
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
        if start_0based >= end_0based_exclusive {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "SummarizeTfbsScoreTracks requires start < end (got {}..{})",
                    start_0based, end_0based_exclusive
                ),
            });
        }
        if end_0based_exclusive > dna.len() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "SummarizeTfbsScoreTracks range {}..{} exceeds sequence length {}",
                    start_0based,
                    end_0based_exclusive,
                    dna.len()
                ),
            });
        }

        let sequence = dna.get_forward_string();
        let bytes = sequence.as_bytes();
        let view = &bytes[start_0based..end_0based_exclusive];
        let random_background = Self::deterministic_random_dna_bytes(
            DEFAULT_JASPAR_PRESENTATION_RANDOM_SEQUENCE_LENGTH_BP,
            DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED,
        );

        let mut tracks = vec![];
        let mut global_max_score = 0.0_f64;
        for motif in motifs {
            let (tf_id, tf_name, _consensus, matrix_counts) =
                Self::resolve_tf_motif_for_scoring(motif)?;
            let (llr_matrix, true_log_odds_matrix) = Self::prepare_scoring_matrices(&matrix_counts);
            let motif_length_bp = llr_matrix.len();
            let scored_window_count = if motif_length_bp == 0 || view.len() < motif_length_bp {
                0
            } else {
                view.len() - motif_length_bp + 1
            };
            let mut forward_scores = Vec::with_capacity(scored_window_count);
            let mut reverse_scores = Vec::with_capacity(scored_window_count);
            let mut max_score = 0.0_f64;
            let mut max_position_0based = None;
            let (
                background_scores,
                llr_background_sorted_scores,
                true_log_odds_background_sorted_scores,
            ) = Self::collect_tfbs_score_track_background_scores(
                &random_background,
                &llr_matrix,
                &true_log_odds_matrix,
                score_kind,
                clip_negative,
            );

            if scored_window_count > 0 {
                let hits =
                    Self::scan_tf_scores(view, &llr_matrix, &true_log_odds_matrix, |_, _| {});
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
                        &llr_background_sorted_scores,
                        &true_log_odds_background_sorted_scores,
                    );
                    if reverse {
                        reverse_scores[offset] = score;
                    } else {
                        forward_scores[offset] = score;
                    }
                }
                for offset in 0..scored_window_count {
                    let forward_score = forward_scores[offset];
                    let reverse_score = reverse_scores[offset];
                    let local_max = forward_score.max(reverse_score);
                    if local_max > max_score {
                        max_score = local_max;
                        max_position_0based = Some(start_0based + offset);
                    }
                }
            }
            let normalization_reference = Self::summarize_tfbs_score_track_normalization_reference(
                &background_scores,
                max_score,
            );
            let top_peaks = Self::summarize_tfbs_score_track_top_peaks(
                start_0based,
                motif_length_bp,
                &forward_scores,
                &reverse_scores,
                &background_scores,
                normalization_reference.as_ref(),
            );
            global_max_score = global_max_score.max(max_score);
            tracks.push(TfbsScoreTrackRow {
                tf_id,
                tf_name,
                motif_length_bp,
                track_start_0based: start_0based,
                scored_window_count,
                max_score,
                max_position_0based,
                normalization_reference,
                top_peaks,
                forward_scores,
                reverse_scores,
            });
        }

        Ok(TfbsScoreTrackReport {
            schema: TFBS_SCORE_TRACK_REPORT_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            sequence_length_bp: dna.len(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            score_kind,
            view_start_0based: start_0based,
            view_end_0based_exclusive: end_0based_exclusive,
            clip_negative,
            motifs_requested: motifs.to_vec(),
            global_max_score,
            tss_markers: Self::summarize_tfbs_score_track_tss_markers(
                dna,
                start_0based,
                end_0based_exclusive,
            ),
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
