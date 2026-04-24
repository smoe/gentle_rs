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

    fn emit_tfbs_score_track_progress(
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
        seq_id: &str,
        motif_id: &str,
        motif_index: usize,
        motif_count: usize,
        stage_index: usize,
        stage_count: usize,
        stage_label: &str,
        detail: &str,
        scanned_steps: usize,
        total_steps: usize,
    ) -> bool {
        let stage_fraction = if total_steps == 0 {
            1.0
        } else {
            (scanned_steps as f64 / total_steps as f64).clamp(0.0, 1.0)
        };
        let motif_fraction = if stage_count == 0 {
            1.0
        } else {
            ((stage_index as f64 - 1.0) + stage_fraction) / stage_count as f64
        }
        .clamp(0.0, 1.0);
        let total_fraction = if motif_count == 0 || stage_count == 0 {
            1.0
        } else {
            (((motif_index.saturating_sub(1) * stage_count) as f64)
                + (stage_index as f64 - 1.0)
                + stage_fraction)
                / (motif_count * stage_count) as f64
        }
        .clamp(0.0, 1.0);
        on_progress(OperationProgress::Tfbs(TfbsProgress {
            seq_id: seq_id.to_string(),
            motif_id: motif_id.to_string(),
            motif_index,
            motif_count,
            scanned_steps,
            total_steps,
            motif_percent: motif_fraction * 100.0,
            total_percent: total_fraction * 100.0,
            task_kind: Some("score_tracks".to_string()),
            stage_label: Some(stage_label.to_string()),
            detail: (!detail.trim().is_empty()).then(|| detail.to_string()),
            stage_percent: Some(stage_fraction * 100.0),
        }))
    }

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

    fn summarize_tfbs_score_track_overlay_tracks(
        &self,
        dna: &DNAsequence,
        start_0based: usize,
        end_0based_exclusive: usize,
    ) -> Vec<TfbsScoreTrackOverlayTrack> {
        let mut grouped = std::collections::BTreeMap::<String, TfbsScoreTrackOverlayTrack>::new();
        for feature in dna.features() {
            if !Self::is_generated_genome_bed_feature(feature) {
                continue;
            }
            let source_kind = Self::feature_qualifier_text(feature, "gentle_track_source")
                .unwrap_or_else(|| "BED".to_string());
            let track_name =
                Self::first_nonempty_feature_qualifier(feature, &["gentle_track_name", "name"])
                    .unwrap_or_else(|| "BED track".to_string());
            let source_file_name = Self::feature_qualifier_text(feature, "gentle_track_file")
                .and_then(|value| {
                    std::path::Path::new(value.trim())
                        .file_name()
                        .and_then(|name| name.to_str())
                        .map(|name| name.to_string())
                });
            let display_label = match source_file_name.as_deref() {
                Some(file_name) if file_name != track_name => {
                    format!("{track_name} ({file_name})")
                }
                _ => track_name.clone(),
            };
            let group_key = match source_file_name.as_deref() {
                Some(file_name) => format!("{track_name}\u{1f}{file_name}"),
                None => track_name.clone(),
            };
            let interval_label = Self::feature_qualifier_text(feature, "label");
            let score = Self::feature_qualifier_text(feature, "score")
                .and_then(|value| value.trim().parse::<f64>().ok());
            let strand = Self::feature_qualifier_text(feature, "bed_strand")
                .or_else(|| Self::feature_qualifier_text(feature, "strand"));

            let mut ranges = vec![];
            collect_location_ranges_usize(&feature.location, &mut ranges);
            for (range_start, range_end) in ranges {
                if range_end <= range_start {
                    continue;
                }
                let clipped_start = range_start.max(start_0based);
                let clipped_end = range_end.min(end_0based_exclusive);
                if clipped_end <= clipped_start {
                    continue;
                }
                let track = grouped.entry(group_key.clone()).or_insert_with(|| {
                    TfbsScoreTrackOverlayTrack {
                        source_kind: source_kind.clone(),
                        track_name: track_name.clone(),
                        display_label: display_label.clone(),
                        source_file_name: source_file_name.clone(),
                        interval_count: 0,
                        max_score: None,
                        intervals: vec![],
                    }
                });
                track.intervals.push(TfbsScoreTrackOverlayInterval {
                    start_0based: clipped_start,
                    end_0based_exclusive: clipped_end,
                    label: interval_label.clone(),
                    score,
                    strand: strand.clone(),
                });
                track.interval_count += 1;
                if let Some(score) = score {
                    let current_max = track.max_score.unwrap_or(score);
                    track.max_score = Some(current_max.max(score));
                }
            }
        }

        let mut tracks = grouped.into_values().collect::<Vec<_>>();
        for track in &mut tracks {
            track.intervals.sort_by(|left, right| {
                left.start_0based
                    .cmp(&right.start_0based)
                    .then(left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
                    .then(left.label.cmp(&right.label))
            });
        }
        tracks.sort_by(|left, right| {
            left.display_label
                .cmp(&right.display_label)
                .then(left.track_name.cmp(&right.track_name))
        });
        tracks
    }

    fn tfbs_track_display_signal(
        track: &TfbsScoreTrackRow,
        signal_source: TfbsScoreTrackCorrelationSignalSource,
    ) -> Vec<f64> {
        match signal_source {
            TfbsScoreTrackCorrelationSignalSource::MaxStrands => track
                .forward_scores
                .iter()
                .copied()
                .zip(track.reverse_scores.iter().copied())
                .map(|(forward, reverse)| forward.max(reverse))
                .collect(),
            TfbsScoreTrackCorrelationSignalSource::ForwardOnly => track.forward_scores.clone(),
            TfbsScoreTrackCorrelationSignalSource::ReverseOnly => track.reverse_scores.clone(),
        }
    }

    fn tfbs_track_strand_signal(
        track: &TfbsScoreTrackRow,
        strand: TfbsScoreTrackStrandComponent,
    ) -> &[f64] {
        match strand {
            TfbsScoreTrackStrandComponent::Forward => &track.forward_scores,
            TfbsScoreTrackStrandComponent::Reverse => &track.reverse_scores,
        }
    }

    fn primary_signal_peak_position(signal: &[f64], track_start_0based: usize) -> Option<usize> {
        signal
            .iter()
            .copied()
            .enumerate()
            .max_by(|left, right| left.1.total_cmp(&right.1).then(left.0.cmp(&right.0)))
            .map(|(idx, _)| track_start_0based + idx)
    }

    fn summarize_tfbs_track_directional_summary(
        track: &TfbsScoreTrackRow,
    ) -> Option<TfbsScoreTrackDirectionalSummary> {
        let overlap_window_count = track.forward_scores.len().min(track.reverse_scores.len());
        if overlap_window_count < 2 {
            return None;
        }
        let forward_signal = &track.forward_scores[..overlap_window_count];
        let reverse_signal = &track.reverse_scores[..overlap_window_count];
        let forward_smoothed = Self::smooth_tfbs_track_signal(
            forward_signal,
            Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
        );
        let reverse_smoothed = Self::smooth_tfbs_track_signal(
            reverse_signal,
            Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
        );
        let forward_primary_peak_position_0based =
            Self::primary_signal_peak_position(forward_signal, track.track_start_0based);
        let reverse_primary_peak_position_0based =
            Self::primary_signal_peak_position(reverse_signal, track.track_start_0based);
        Some(TfbsScoreTrackDirectionalSummary {
            forward_primary_peak_position_0based,
            reverse_primary_peak_position_0based,
            signed_primary_peak_offset_bp: match (
                forward_primary_peak_position_0based,
                reverse_primary_peak_position_0based,
            ) {
                (Some(forward), Some(reverse)) => Some(reverse as i64 - forward as i64),
                _ => None,
            },
            raw_pearson: Self::pearson_correlation(forward_signal, reverse_signal),
            smoothed_pearson: Self::pearson_correlation(&forward_smoothed, &reverse_smoothed),
            raw_spearman: Self::spearman_correlation(forward_signal, reverse_signal),
            smoothed_spearman: Self::spearman_correlation(&forward_smoothed, &reverse_smoothed),
        })
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

    fn rank_signal_for_spearman(signal: &[f64]) -> Vec<f64> {
        if signal.is_empty() {
            return vec![];
        }
        let mut indexed = signal.iter().copied().enumerate().collect::<Vec<_>>();
        indexed.sort_by(|left, right| left.1.total_cmp(&right.1).then(left.0.cmp(&right.0)));
        let mut ranks = vec![0.0; signal.len()];
        let mut idx = 0usize;
        while idx < indexed.len() {
            let start = idx;
            let value = indexed[idx].1;
            idx += 1;
            while idx < indexed.len() && indexed[idx].1 == value {
                idx += 1;
            }
            let average_rank = ((start + 1) as f64 + idx as f64) / 2.0;
            for entry in &indexed[start..idx] {
                ranks[entry.0] = average_rank;
            }
        }
        ranks
    }

    fn spearman_correlation(left: &[f64], right: &[f64]) -> f64 {
        let len = left.len().min(right.len());
        if len < 2 {
            return 0.0;
        }
        let left_ranks = Self::rank_signal_for_spearman(&left[..len]);
        let right_ranks = Self::rank_signal_for_spearman(&right[..len]);
        Self::pearson_correlation(&left_ranks, &right_ranks)
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

    fn summarize_tfbs_score_track_correlation_summary_for_source(
        tracks: &[TfbsScoreTrackRow],
        signal_source: TfbsScoreTrackCorrelationSignalSource,
    ) -> Option<TfbsScoreTrackCorrelationSummary> {
        if tracks.len() < 2 {
            return None;
        }
        let mut rows = vec![];
        for left_idx in 0..tracks.len() {
            let left = &tracks[left_idx];
            let left_signal = Self::tfbs_track_display_signal(left, signal_source);
            if left_signal.len() < 2 {
                continue;
            }
            let left_smoothed = Self::smooth_tfbs_track_signal(
                &left_signal,
                Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
            );
            for right in tracks.iter().skip(left_idx + 1) {
                let right_signal = Self::tfbs_track_display_signal(right, signal_source);
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
                    raw_spearman: Self::spearman_correlation(
                        &left_signal[..overlap_window_count],
                        &right_signal[..overlap_window_count],
                    ),
                    smoothed_spearman: Self::spearman_correlation(
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
            signal_source,
            smoothing_method: "centered_boxcar".to_string(),
            smoothing_window_bp: Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
            pair_count: rows.len(),
            rows,
        })
    }

    fn summarize_tfbs_score_track_correlation_summaries(
        tracks: &[TfbsScoreTrackRow],
    ) -> Vec<TfbsScoreTrackCorrelationSummary> {
        [
            TfbsScoreTrackCorrelationSignalSource::MaxStrands,
            TfbsScoreTrackCorrelationSignalSource::ForwardOnly,
            TfbsScoreTrackCorrelationSignalSource::ReverseOnly,
        ]
        .into_iter()
        .filter_map(|signal_source| {
            Self::summarize_tfbs_score_track_correlation_summary_for_source(tracks, signal_source)
        })
        .collect()
    }

    fn summarize_tfbs_score_track_cross_strand_correlation_summary(
        tracks: &[TfbsScoreTrackRow],
    ) -> Option<TfbsScoreTrackCrossStrandCorrelationSummary> {
        if tracks.len() < 2 {
            return None;
        }
        let mut rows = vec![];
        for left_idx in 0..tracks.len() {
            let left = &tracks[left_idx];
            for right in tracks.iter().skip(left_idx + 1) {
                let mut cells = vec![];
                for (left_strand, right_strand) in [
                    (
                        TfbsScoreTrackStrandComponent::Forward,
                        TfbsScoreTrackStrandComponent::Forward,
                    ),
                    (
                        TfbsScoreTrackStrandComponent::Forward,
                        TfbsScoreTrackStrandComponent::Reverse,
                    ),
                    (
                        TfbsScoreTrackStrandComponent::Reverse,
                        TfbsScoreTrackStrandComponent::Forward,
                    ),
                    (
                        TfbsScoreTrackStrandComponent::Reverse,
                        TfbsScoreTrackStrandComponent::Reverse,
                    ),
                ] {
                    let left_signal = Self::tfbs_track_strand_signal(left, left_strand);
                    let right_signal = Self::tfbs_track_strand_signal(right, right_strand);
                    let overlap_window_count = left_signal.len().min(right_signal.len());
                    if overlap_window_count < 2 {
                        continue;
                    }
                    let left_signal = &left_signal[..overlap_window_count];
                    let right_signal = &right_signal[..overlap_window_count];
                    let left_smoothed = Self::smooth_tfbs_track_signal(
                        left_signal,
                        Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
                    );
                    let right_smoothed = Self::smooth_tfbs_track_signal(
                        right_signal,
                        Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
                    );
                    let left_peak =
                        Self::primary_signal_peak_position(left_signal, left.track_start_0based);
                    let right_peak =
                        Self::primary_signal_peak_position(right_signal, right.track_start_0based);
                    cells.push(TfbsScoreTrackCrossStrandCorrelationCell {
                        left_strand,
                        right_strand,
                        overlap_window_count,
                        raw_pearson: Self::pearson_correlation(left_signal, right_signal),
                        smoothed_pearson: Self::pearson_correlation(
                            &left_smoothed,
                            &right_smoothed,
                        ),
                        raw_spearman: Self::spearman_correlation(left_signal, right_signal),
                        smoothed_spearman: Self::spearman_correlation(
                            &left_smoothed,
                            &right_smoothed,
                        ),
                        signed_primary_peak_offset_bp: match (left_peak, right_peak) {
                            (Some(left_peak), Some(right_peak)) => {
                                Some(right_peak as i64 - left_peak as i64)
                            }
                            _ => None,
                        },
                    });
                }
                if cells.is_empty() {
                    continue;
                }
                rows.push(TfbsScoreTrackCrossStrandCorrelationRow {
                    left_tf_id: left.tf_id.clone(),
                    left_tf_name: left.tf_name.clone(),
                    right_tf_id: right.tf_id.clone(),
                    right_tf_name: right.tf_name.clone(),
                    cells,
                });
            }
        }
        if rows.is_empty() {
            return None;
        }
        rows.sort_by(|left, right| {
            let left_smoothed = left
                .cells
                .iter()
                .map(|cell| cell.smoothed_spearman.abs())
                .fold(0.0_f64, f64::max);
            let right_smoothed = right
                .cells
                .iter()
                .map(|cell| cell.smoothed_spearman.abs())
                .fold(0.0_f64, f64::max);
            let left_raw = left
                .cells
                .iter()
                .map(|cell| cell.raw_spearman.abs())
                .fold(0.0_f64, f64::max);
            let right_raw = right
                .cells
                .iter()
                .map(|cell| cell.raw_spearman.abs())
                .fold(0.0_f64, f64::max);
            right_smoothed
                .total_cmp(&left_smoothed)
                .then(right_raw.total_cmp(&left_raw))
                .then(left.left_tf_id.cmp(&right.left_tf_id))
                .then(left.right_tf_id.cmp(&right.right_tf_id))
        });
        Some(TfbsScoreTrackCrossStrandCorrelationSummary {
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
        mut on_progress: impl FnMut(usize, usize) -> bool,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), EngineError> {
        let hits = Self::scan_tf_scores_with_topology_and_cancel(
            random_background,
            llr_matrix,
            true_log_odds_matrix,
            InlineSequenceTopology::Linear,
            |scanned_steps, total_steps| on_progress(scanned_steps, total_steps),
        )?;
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
        Ok((
            displayed_scores,
            underlying_background_scores,
            llr_background_scores,
            true_log_odds_background_scores,
        ))
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

    fn summarize_tfbs_score_tracks_internal(
        &self,
        target: SequenceScanTarget,
        motifs: &[String],
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        include_correlation_summary: bool,
    ) -> Result<TfbsScoreTrackReport, EngineError> {
        self.summarize_tfbs_score_tracks_with_progress(
            target,
            motifs,
            score_kind,
            clip_negative,
            include_correlation_summary,
            &mut |_| true,
        )
    }

    pub(crate) fn summarize_tfbs_score_tracks_with_progress(
        &self,
        target: SequenceScanTarget,
        motifs: &[String],
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        include_correlation_summary: bool,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<TfbsScoreTrackReport, EngineError> {
        if motifs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeTfbsScoreTracks requires at least one motif".to_string(),
            });
        }
        let effective_motifs = Self::expand_tf_query_tokens(motifs)?;
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
        let motif_count = effective_motifs.len();
        const SCORE_TRACK_STAGE_COUNT: usize = 2;
        for (motif_idx, motif) in effective_motifs.iter().enumerate() {
            let (tf_id, tf_name, _consensus, matrix_counts) =
                Self::resolve_tf_motif_for_scoring(motif)?;
            let (llr_matrix, true_log_odds_matrix) = Self::prepare_scoring_matrices(&matrix_counts);
            let llr_modeled_distribution = Self::modeled_tfbs_score_distribution(&llr_matrix);
            let true_log_odds_modeled_distribution =
                Self::modeled_tfbs_score_distribution(&true_log_odds_matrix);
            let motif_length_bp = llr_matrix.len();
            let motif_logo_columns = Self::jaspar_expert_columns(&matrix_counts);
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
                |scanned_steps, total_steps| {
                    Self::emit_tfbs_score_track_progress(
                        on_progress,
                        &target_label,
                        &tf_id,
                        motif_idx + 1,
                        motif_count,
                        1,
                        SCORE_TRACK_STAGE_COUNT,
                        "background calibration",
                        &format!(
                            "{} bp deterministic random DNA",
                            DEFAULT_TFBS_SCORE_TRACK_RANDOM_SEQUENCE_LENGTH_BP
                        ),
                        scanned_steps,
                        total_steps,
                    )
                },
            )?;

            if scored_window_count > 0 {
                let hits = Self::scan_tf_scores_with_topology_and_cancel(
                    view,
                    &llr_matrix,
                    &true_log_odds_matrix,
                    scan_topology,
                    |scanned_steps, total_steps| {
                        Self::emit_tfbs_score_track_progress(
                            on_progress,
                            &target_label,
                            &tf_id,
                            motif_idx + 1,
                            motif_count,
                            2,
                            SCORE_TRACK_STAGE_COUNT,
                            "target scan",
                            &format!(
                                "{}..{} ({})",
                                scan_start_0based,
                                scan_end_0based_exclusive,
                                score_kind.as_str()
                            ),
                            scanned_steps,
                            total_steps,
                        )
                    },
                )?;
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
            } else {
                if !Self::emit_tfbs_score_track_progress(
                    on_progress,
                    &target_label,
                    &tf_id,
                    motif_idx + 1,
                    motif_count,
                    2,
                    SCORE_TRACK_STAGE_COUNT,
                    "target scan",
                    &format!(
                        "{}..{} ({})",
                        scan_start_0based,
                        scan_end_0based_exclusive,
                        score_kind.as_str()
                    ),
                    1,
                    1,
                ) {
                    return Err(Self::tfbs_cancelled_error("target scan"));
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
            let mut track = TfbsScoreTrackRow {
                tf_id,
                tf_name: tf_name.or_else(|| Some(motif.clone())),
                motif_length_bp,
                motif_logo_columns,
                track_start_0based: scan_start_0based,
                scored_window_count,
                max_score,
                max_position_0based,
                normalization_reference,
                directional_summary: None,
                top_peaks,
                forward_scores,
                reverse_scores,
            };
            track.directional_summary = Self::summarize_tfbs_track_directional_summary(&track);
            tracks.push(track);
        }

        let (tss_markers, overlay_tracks) = if let Some(seq_id) = tss_source_seq_id {
            self.state
                .sequences
                .get(&seq_id)
                .map(|dna| {
                    (
                        self.summarize_tfbs_score_track_tss_markers(
                            &seq_id,
                            dna,
                            scan_start_0based,
                            scan_end_0based_exclusive,
                        ),
                        self.summarize_tfbs_score_track_overlay_tracks(
                            dna,
                            scan_start_0based,
                            scan_end_0based_exclusive,
                        ),
                    )
                })
                .unwrap_or_else(|| (vec![], vec![]))
        } else {
            (vec![], vec![])
        };

        let correlation_summaries = if include_correlation_summary {
            Self::summarize_tfbs_score_track_correlation_summaries(&tracks)
        } else {
            vec![]
        };
        let cross_strand_correlation_summary = if include_correlation_summary {
            Self::summarize_tfbs_score_track_cross_strand_correlation_summary(&tracks)
        } else {
            None
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
            overlay_tracks,
            correlation_summary: correlation_summaries
                .iter()
                .find(|summary| {
                    summary.signal_source == TfbsScoreTrackCorrelationSignalSource::MaxStrands
                })
                .cloned(),
            correlation_summaries,
            cross_strand_correlation_summary,
            tracks,
        })
    }

    pub(crate) fn summarize_tfbs_score_tracks(
        &self,
        target: SequenceScanTarget,
        motifs: &[String],
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
    ) -> Result<TfbsScoreTrackReport, EngineError> {
        self.summarize_tfbs_score_tracks_internal(target, motifs, score_kind, clip_negative, true)
    }

    fn tfbs_track_similarity_metric_value(
        row: &TfbsTrackSimilarityRow,
        metric: TfbsTrackSimilarityRankingMetric,
    ) -> f64 {
        match metric {
            TfbsTrackSimilarityRankingMetric::RawPearson => row.raw_pearson,
            TfbsTrackSimilarityRankingMetric::SmoothedPearson => row.smoothed_pearson,
            TfbsTrackSimilarityRankingMetric::RawSpearman => row.raw_spearman,
            TfbsTrackSimilarityRankingMetric::SmoothedSpearman => row.smoothed_spearman,
        }
    }

    fn summarize_tfbs_track_similarity_row(
        anchor: &TfbsScoreTrackRow,
        candidate: &TfbsScoreTrackRow,
        remote_summary: Option<JasparCatalogRemoteSummary>,
    ) -> Option<TfbsTrackSimilarityRow> {
        let anchor_signal = Self::tfbs_track_display_signal(
            anchor,
            TfbsScoreTrackCorrelationSignalSource::MaxStrands,
        );
        let candidate_signal = Self::tfbs_track_display_signal(
            candidate,
            TfbsScoreTrackCorrelationSignalSource::MaxStrands,
        );
        let overlap_window_count = anchor_signal.len().min(candidate_signal.len());
        if overlap_window_count < 2 {
            return None;
        }
        let anchor_smoothed = Self::smooth_tfbs_track_signal(
            &anchor_signal,
            Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
        );
        let candidate_smoothed = Self::smooth_tfbs_track_signal(
            &candidate_signal,
            Self::TFBS_TRACK_CORRELATION_SMOOTHING_WINDOW_BP,
        );
        Some(TfbsTrackSimilarityRow {
            candidate_tf_id: candidate.tf_id.clone(),
            candidate_tf_name: candidate.tf_name.clone(),
            overlap_window_count,
            raw_pearson: Self::pearson_correlation(
                &anchor_signal[..overlap_window_count],
                &candidate_signal[..overlap_window_count],
            ),
            smoothed_pearson: Self::pearson_correlation(
                &anchor_smoothed[..overlap_window_count],
                &candidate_smoothed[..overlap_window_count],
            ),
            raw_spearman: Self::spearman_correlation(
                &anchor_signal[..overlap_window_count],
                &candidate_signal[..overlap_window_count],
            ),
            smoothed_spearman: Self::spearman_correlation(
                &anchor_smoothed[..overlap_window_count],
                &candidate_smoothed[..overlap_window_count],
            ),
            signed_primary_peak_offset_bp: Self::summarize_tfbs_track_primary_peak_offset_bp(
                anchor, candidate,
            ),
            remote_summary,
        })
    }

    fn jaspar_species_filter_matches(
        row: &JasparRemoteMetadataSnapshotRow,
        filters: &[String],
    ) -> bool {
        if filters.is_empty() {
            return true;
        }
        filters.iter().any(|filter| {
            let filter = filter.trim().to_ascii_lowercase();
            if filter.is_empty() {
                return false;
            }
            row.remote_metadata
                .species_assignments
                .iter()
                .any(|assignment| {
                    assignment
                        .scientific_name
                        .to_ascii_lowercase()
                        .contains(&filter)
                        || assignment
                            .common_name
                            .as_deref()
                            .unwrap_or("")
                            .to_ascii_lowercase()
                            .contains(&filter)
                        || assignment
                            .tax_id
                            .as_deref()
                            .unwrap_or("")
                            .to_ascii_lowercase()
                            .contains(&filter)
                })
        })
    }

    pub(crate) fn summarize_tfbs_track_similarity(
        &self,
        target: SequenceScanTarget,
        anchor_motif: &str,
        candidate_motifs: &[String],
        ranking_metric: TfbsTrackSimilarityRankingMetric,
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        species_filters: &[String],
        include_remote_metadata: bool,
        limit: Option<usize>,
    ) -> Result<TfbsTrackSimilarityReport, EngineError> {
        self.summarize_tfbs_track_similarity_with_snapshot_path(
            target,
            anchor_motif,
            candidate_motifs,
            ranking_metric,
            score_kind,
            clip_negative,
            species_filters,
            include_remote_metadata,
            limit,
            None,
        )
    }

    pub(crate) fn summarize_tfbs_track_similarity_with_snapshot_path(
        &self,
        target: SequenceScanTarget,
        anchor_motif: &str,
        candidate_motifs: &[String],
        ranking_metric: TfbsTrackSimilarityRankingMetric,
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        species_filters: &[String],
        include_remote_metadata: bool,
        limit: Option<usize>,
        snapshot_path: Option<&str>,
    ) -> Result<TfbsTrackSimilarityReport, EngineError> {
        let anchor_requested = anchor_motif.trim().to_string();
        if anchor_requested.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeTfbsTrackSimilarity requires non-empty anchor_motif".to_string(),
            });
        }
        let (anchor_tf_id, anchor_tf_name_resolved, _anchor_consensus, _anchor_matrix_counts) =
            Self::resolve_tf_motif_for_scoring(&anchor_requested)?;
        let candidate_scope = if candidate_motifs.is_empty()
            || (candidate_motifs.len() == 1
                && matches!(
                    candidate_motifs[0].trim().to_ascii_uppercase().as_str(),
                    "ALL" | "*"
                )) {
            "all_registry".to_string()
        } else {
            "explicit".to_string()
        };
        let requested_candidate_tokens = if candidate_scope == "all_registry" {
            crate::tf_motifs::all_motif_ids()
        } else {
            Self::expand_tf_query_tokens(candidate_motifs)?
        };
        if requested_candidate_tokens.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "SummarizeTfbsTrackSimilarity requires at least one candidate motif or ALL"
                        .to_string(),
            });
        }

        let normalized_species_filters = species_filters
            .iter()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();

        let needs_remote_snapshot =
            include_remote_metadata || !normalized_species_filters.is_empty();
        let mut warnings = vec![];
        let cached_remote_rows = if needs_remote_snapshot {
            match Self::load_jaspar_remote_metadata_snapshot_rows(snapshot_path) {
                Ok(rows) => rows,
                Err(err) => {
                    warnings.push(err);
                    std::collections::BTreeMap::new()
                }
            }
        } else {
            std::collections::BTreeMap::new()
        };

        let mut seen_candidate_ids = std::collections::BTreeSet::new();
        let mut candidate_ids = vec![];
        let mut species_excluded_count = 0usize;
        let mut species_metadata_missing_count = 0usize;
        for token in requested_candidate_tokens.iter() {
            let token = token.trim();
            if token.is_empty() {
                continue;
            }
            let (candidate_tf_id, _candidate_tf_name, _consensus, _matrix_counts) =
                Self::resolve_tf_motif_for_scoring(token)?;
            if candidate_tf_id == anchor_tf_id {
                continue;
            }
            if !normalized_species_filters.is_empty() {
                match cached_remote_rows.get(&candidate_tf_id) {
                    Some(row) => {
                        if !Self::jaspar_species_filter_matches(row, &normalized_species_filters) {
                            species_excluded_count = species_excluded_count.saturating_add(1);
                            continue;
                        }
                    }
                    None => {
                        species_metadata_missing_count =
                            species_metadata_missing_count.saturating_add(1);
                        continue;
                    }
                }
            }
            if seen_candidate_ids.insert(candidate_tf_id.clone()) {
                candidate_ids.push(candidate_tf_id);
            }
        }
        if !normalized_species_filters.is_empty() {
            if species_excluded_count > 0 {
                warnings.push(format!(
                    "Species filter excluded {} candidate motif(s) from the similarity ranking",
                    species_excluded_count
                ));
            }
            if species_metadata_missing_count > 0 {
                warnings.push(format!(
                    "Species filter skipped {} candidate motif(s) without cached JASPAR remote metadata",
                    species_metadata_missing_count
                ));
            }
        }

        let mut motifs_to_score = Vec::with_capacity(candidate_ids.len().saturating_add(1));
        motifs_to_score.push(anchor_tf_id.clone());
        motifs_to_score.extend(candidate_ids.iter().cloned());
        let base_report = self.summarize_tfbs_score_tracks_internal(
            target,
            &motifs_to_score,
            score_kind,
            clip_negative,
            false,
        )?;
        let anchor_track = base_report
            .tracks
            .iter()
            .find(|track| track.tf_id == anchor_tf_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Could not recover anchor motif '{}' from score-track summary",
                    anchor_tf_id
                ),
            })?;

        let mut rows = base_report
            .tracks
            .iter()
            .filter(|track| track.tf_id != anchor_tf_id)
            .filter_map(|track| {
                Self::summarize_tfbs_track_similarity_row(
                    anchor_track,
                    track,
                    if include_remote_metadata {
                        cached_remote_rows
                            .get(&track.tf_id)
                            .and_then(|row| row.remote_summary.clone())
                            .or_else(|| {
                                cached_remote_rows.get(&track.tf_id).map(|row| {
                                    Self::jaspar_remote_metadata_summary(&row.remote_metadata)
                                })
                            })
                    } else {
                        None
                    },
                )
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            Self::tfbs_track_similarity_metric_value(right, ranking_metric)
                .total_cmp(&Self::tfbs_track_similarity_metric_value(
                    left,
                    ranking_metric,
                ))
                .then(right.smoothed_spearman.total_cmp(&left.smoothed_spearman))
                .then(right.smoothed_pearson.total_cmp(&left.smoothed_pearson))
                .then_with(|| {
                    left.signed_primary_peak_offset_bp
                        .unwrap_or(i64::MAX)
                        .abs()
                        .cmp(
                            &right
                                .signed_primary_peak_offset_bp
                                .unwrap_or(i64::MAX)
                                .abs(),
                        )
                })
                .then(left.candidate_tf_id.cmp(&right.candidate_tf_id))
        });
        let scanned_candidate_count = rows.len();
        if let Some(limit) = limit {
            if rows.len() > limit {
                rows.truncate(limit);
            }
        }

        Ok(TfbsTrackSimilarityReport {
            schema: TFBS_TRACK_SIMILARITY_REPORT_SCHEMA.to_string(),
            target_kind: base_report.target_kind,
            target_label: base_report.target_label.clone(),
            seq_id: base_report.seq_id,
            source_sequence_length_bp: base_report.source_sequence_length_bp,
            scan_topology: base_report.scan_topology,
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            score_kind,
            view_start_0based: base_report.view_start_0based,
            view_end_0based_exclusive: base_report.view_end_0based_exclusive,
            clip_negative,
            anchor_requested,
            anchor_tf_id,
            anchor_tf_name: anchor_track.tf_name.clone().or(anchor_tf_name_resolved),
            candidate_scope,
            candidates_requested: candidate_motifs.to_vec(),
            species_filters: normalized_species_filters,
            include_remote_metadata,
            ranking_metric,
            scanned_candidate_count,
            returned_candidate_count: rows.len(),
            warnings,
            rows,
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
        let effective_motifs = Self::expand_tf_query_tokens(motifs)?;

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

    fn promoter_aligned_sequence(sequence: &str, strand: Option<char>) -> String {
        if strand == Some('-') {
            Self::reverse_complement(sequence)
        } else {
            sequence.to_string()
        }
    }

    fn promoter_oriented_tss_position_0based(
        sequence_length_bp: usize,
        extract_start_1based: usize,
        tss_1based: usize,
        strand: Option<char>,
    ) -> Option<usize> {
        if extract_start_1based == 0 || tss_1based < extract_start_1based {
            return None;
        }
        let local_0based = tss_1based.saturating_sub(extract_start_1based);
        if local_0based >= sequence_length_bp {
            return None;
        }
        Some(if strand == Some('-') {
            sequence_length_bp
                .saturating_sub(1)
                .saturating_sub(local_0based)
        } else {
            local_0based
        })
    }

    fn promoter_local_position_to_genomic_1based(
        strand: Option<char>,
        promoter_start_1based: usize,
        promoter_end_1based: usize,
        promoter_length_bp: usize,
        local_0based: usize,
    ) -> Option<usize> {
        if promoter_length_bp == 0 || local_0based >= promoter_length_bp {
            return None;
        }
        Some(if strand == Some('-') {
            promoter_end_1based.saturating_sub(local_0based)
        } else {
            promoter_start_1based.saturating_add(local_0based)
        })
    }

    fn tfbs_track_positive_support_fraction(track: &TfbsScoreTrackRow) -> f64 {
        if track.scored_window_count == 0 {
            return 0.0;
        }
        let positive_windows = track
            .forward_scores
            .iter()
            .zip(track.reverse_scores.iter())
            .filter(|(forward, reverse)| **forward > 0.0 || **reverse > 0.0)
            .count();
        positive_windows as f64 / track.scored_window_count as f64
    }

    pub(crate) fn summarize_multi_gene_promoter_tfbs(
        &self,
        genome_id: &str,
        genes: &[PromoterTfbsGeneQuery],
        motifs: &[String],
        upstream_bp: usize,
        downstream_bp: usize,
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<MultiGenePromoterTfbsReport, EngineError> {
        if genes.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeMultiGenePromoterTfbs requires at least one gene query"
                    .to_string(),
            });
        }
        if motifs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeMultiGenePromoterTfbs requires at least one TF motif query"
                    .to_string(),
            });
        }

        let effective_catalog_path =
            catalog_path.unwrap_or(crate::genomes::default_catalog_discovery_token(false));
        let (catalog, _) = Self::open_reference_genome_catalog(Some(effective_catalog_path))?;
        let mut gene_reports = vec![];
        let mut summary_rows = vec![];
        let mut warnings = vec![];

        for gene in genes {
            let resolved = Self::resolve_genome_promoter_slice_request(
                &catalog,
                genome_id,
                &gene.gene_query,
                gene.occurrence,
                gene.transcript_id.as_deref(),
                upstream_bp,
                downstream_bp,
                cache_dir,
            )?;
            warnings.extend(resolved.warnings.clone());

            let strand = resolved
                .selected_gene
                .strand
                .or(resolved.selected_transcript.strand)
                .unwrap_or('+');
            let promoter_sequence = catalog
                .get_sequence_region_with_cache(
                    genome_id,
                    &resolved.selected_transcript.chromosome,
                    resolved.extract_start_1based,
                    resolved.extract_end_1based,
                    cache_dir,
                )
                .map_err(|e| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Could not load promoter slice {}:{}-{} from '{}': {}",
                        resolved.selected_transcript.chromosome,
                        resolved.extract_start_1based,
                        resolved.extract_end_1based,
                        genome_id,
                        e
                    ),
                })?;
            let oriented_sequence =
                Self::promoter_aligned_sequence(&promoter_sequence, Some(strand));
            let promoter_length_bp = oriented_sequence.len();
            let tss_position_0based = Self::promoter_oriented_tss_position_0based(
                promoter_length_bp,
                resolved.extract_start_1based,
                resolved.tss_1based,
                Some(strand),
            )
            .ok_or_else(|| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Could not place TSS for gene '{}' into promoter-oriented coordinates",
                    Self::genome_gene_display_label(&resolved.selected_gene)
                ),
            })?;
            let desired_promoter_length =
                upstream_bp.saturating_add(downstream_bp).saturating_add(1);
            if promoter_length_bp != desired_promoter_length {
                warnings.push(format!(
                    "Promoter slice for '{}' is {} bp instead of the requested {} bp, likely because the genomic interval clipped at a contig boundary.",
                    Self::genome_gene_display_label(&resolved.selected_gene),
                    promoter_length_bp,
                    desired_promoter_length
                ));
            }

            let display_label = gene
                .display_label
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string)
                .or_else(|| resolved.selected_gene.gene_name.clone())
                .unwrap_or_else(|| resolved.query.clone());
            let mut tfbs_score_tracks = self.summarize_tfbs_score_tracks(
                SequenceScanTarget::InlineSequence {
                    sequence_text: oriented_sequence,
                    topology: InlineSequenceTopology::Linear,
                    id_hint: Some(display_label.clone()),
                    span_start_0based: None,
                    span_end_0based_exclusive: None,
                },
                motifs,
                score_kind,
                clip_negative,
            )?;
            tfbs_score_tracks.tss_markers = vec![TfbsScoreTrackTssMarker {
                feature_id: usize::MAX,
                feature_kind: "genome_promoter_slice".to_string(),
                label: resolved.selected_transcript.transcript_id.clone(),
                position_0based: tss_position_0based,
                is_reverse: false,
            }];

            for track in &tfbs_score_tracks.tracks {
                let peak_position_0based = track.max_position_0based;
                summary_rows.push(MultiGenePromoterTfbsSummaryRow {
                    gene_label: display_label.clone(),
                    gene_query: resolved.query.clone(),
                    transcript_id: resolved.selected_transcript.transcript_id.clone(),
                    tf_id: track.tf_id.clone(),
                    tf_name: track.tf_name.clone(),
                    max_score: track.max_score,
                    peak_position_0based,
                    peak_position_promoter_relative_bp: peak_position_0based
                        .map(|value| value as i64 - tss_position_0based as i64),
                    peak_genomic_position_1based: peak_position_0based.and_then(|value| {
                        Self::promoter_local_position_to_genomic_1based(
                            Some(strand),
                            resolved.extract_start_1based,
                            resolved.extract_end_1based,
                            promoter_length_bp,
                            value,
                        )
                    }),
                    positive_fraction: Self::tfbs_track_positive_support_fraction(track),
                });
            }

            gene_reports.push(MultiGenePromoterTfbsGeneReport {
                gene_query: resolved.query.clone(),
                occurrence: resolved.occurrence,
                transcript_id_requested: gene.transcript_id.clone(),
                display_label,
                gene_id: resolved.selected_gene.gene_id.clone(),
                gene_name: resolved.selected_gene.gene_name.clone(),
                transcript_id: resolved.selected_transcript.transcript_id.clone(),
                chromosome: resolved.selected_transcript.chromosome.clone(),
                strand: strand.to_string(),
                promoter_start_1based: resolved.extract_start_1based,
                promoter_end_1based: resolved.extract_end_1based,
                promoter_length_bp,
                tss_1based: resolved.tss_1based,
                sequence_orientation: "transcription_aligned".to_string(),
                used_fuzzy_gene_match: resolved.used_fuzzy_gene_match,
                tfbs_score_tracks,
            });
        }

        summary_rows.sort_by(|left, right| {
            left.gene_label
                .cmp(&right.gene_label)
                .then(left.transcript_id.cmp(&right.transcript_id))
                .then(left.tf_id.cmp(&right.tf_id))
        });

        Ok(MultiGenePromoterTfbsReport {
            schema: MULTI_GENE_PROMOTER_TFBS_REPORT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            genome_id: genome_id.to_string(),
            upstream_bp,
            downstream_bp,
            score_kind,
            clip_negative,
            motifs_requested: motifs.to_vec(),
            gene_queries_requested: genes.to_vec(),
            returned_gene_count: gene_reports.len(),
            genes: gene_reports,
            summary_rows,
            warnings,
        })
    }
}
