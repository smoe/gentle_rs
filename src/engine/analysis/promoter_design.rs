//! Shared promoter-design helpers that extend beyond one specific variant.
//!
//! This module currently owns continuous TF motif score tracks used by the
//! GUI-side Promoter design expert and by headless JSON export paths.

use super::*;

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

    pub(crate) fn summarize_tfbs_score_tracks(
        &self,
        seq_id: &str,
        motifs: &[String],
        start_0based: usize,
        end_0based_exclusive: usize,
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
            let (llr_matrix, _true_log_odds_matrix) =
                Self::prepare_scoring_matrices(&matrix_counts);
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
                for offset in 0..scored_window_count {
                    let window = &view[offset..offset + motif_length_bp];
                    let forward_raw = Self::score_matrix_window(window, &llr_matrix).unwrap_or(0.0);
                    let reverse_window = Self::reverse_complement_bytes(window);
                    let reverse_raw =
                        Self::score_matrix_window(&reverse_window, &llr_matrix).unwrap_or(0.0);
                    let forward_score =
                        Self::promoter_design_clip_score(forward_raw, clip_negative);
                    let reverse_score =
                        Self::promoter_design_clip_score(reverse_raw, clip_negative);
                    let local_max = forward_score.max(reverse_score);
                    if local_max > max_score {
                        max_score = local_max;
                        max_position_0based = Some(start_0based + offset);
                    }
                    forward_scores.push(forward_score);
                    reverse_scores.push(reverse_score);
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
            score_kind: "llr_bits".to_string(),
            view_start_0based: start_0based,
            view_end_0based_exclusive: end_0based_exclusive,
            clip_negative,
            motifs_requested: motifs.to_vec(),
            global_max_score,
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
