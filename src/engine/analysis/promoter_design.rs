//! Shared promoter-design helpers that extend beyond one specific variant.
//!
//! This module currently owns continuous TF motif score tracks used by the
//! GUI-side Promoter design expert and by headless JSON export paths.

use super::*;
use crate::feature_location::{feature_is_reverse, feature_ranges_sorted_i64};

impl GentleEngine {
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
                    let raw_score = match score_kind {
                        TfbsScoreTrackValueKind::LlrBits => llr_bits,
                        TfbsScoreTrackValueKind::LlrQuantile => llr_quantile,
                        TfbsScoreTrackValueKind::TrueLogOddsBits => true_log_odds_bits,
                        TfbsScoreTrackValueKind::TrueLogOddsQuantile => true_log_odds_quantile,
                    };
                    let score = if score_kind.supports_negative_values() {
                        Self::promoter_design_clip_score(raw_score, clip_negative)
                    } else {
                        raw_score
                    };
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
            global_max_score = global_max_score.max(max_score);
            tracks.push(TfbsScoreTrackRow {
                tf_id,
                tf_name,
                motif_length_bp,
                track_start_0based: start_0based,
                scored_window_count,
                max_score,
                max_position_0based,
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
