//! Sequence-bound presentation of explicit, shared-engine local scoring.
//! No imported evidence is replaced or reinterpreted by these lanes.

use super::*;
use crate::engine::{
    GentleEngine, InlineSequenceTopology, OperationProgress, SequenceScanTarget,
    TfbsScoreTrackMatrixBinding, TfbsScoreTrackReport, TfbsScoreTrackValueKind,
};
use serde::Deserialize;

/// Independent TSS settings; never borrowed from a different TFBS panel.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub struct TssLocalScoreRequest {
    pub matrix_ids: Vec<String>,
    pub score_kind: TfbsScoreTrackValueKind,
    pub clip_negative: bool,
}

/// Retained with the view/export; the score report remains the scientific result.
#[derive(Clone, Debug, Serialize)]
pub struct TssLocalScoreAttachment {
    pub request: TssLocalScoreRequest,
    pub source_binding_sha256: String,
    pub cache_key_sha256: String,
    pub report_sha256: String,
    pub producer_revision: String,
    pub provenance: crate::engine::TfbsScoreTrackProvenance,
}

fn digest<T: Serialize>(value: &T) -> Result<String, String> {
    serde_json::to_vec(value)
        .map(|bytes| sha256_hex_bytes(&bytes))
        .map_err(|e| e.to_string())
}

impl TssLocalScoreRequest {
    /// A bounded convenience surface, not a new scorer or a limit on headless analysis.
    pub fn validate_budget(&self, length: usize) -> Result<(), String> {
        if self.matrix_ids.is_empty() || self.matrix_ids.len() > 32 {
            return Err("Local TSS scoring requires 1..32 explicit matrix IDs (no ALL, aliases or IUPAC fallback)".into());
        }
        if length == 0 || length > 50_000 || length * self.matrix_ids.len() > 1_000_000 {
            return Err("Local TSS scoring budget: at most 50,000 bp and 1,000,000 matrix/base combinations; use a smaller annotated window or the shared headless scorer".into());
        }
        let mut seen = std::collections::BTreeSet::new();
        if self
            .matrix_ids
            .iter()
            .any(|id| id.trim() != id || id.is_empty() || !seen.insert(id))
        {
            return Err("Local matrix IDs must be nonempty, exact and unique".into());
        }
        Ok(())
    }
}

impl TssSequenceView {
    /// Reference, locus, orientation and bases are independent from attached evidence.
    pub fn local_score_source_binding(&self) -> Result<String, String> {
        digest(&(
            &self.promoter_id,
            &self.assembly,
            &self.genome_id,
            &self.annotation_release,
            &self.geometry,
            &self.sequence_sha256,
        ))
    }

    /// Uses precisely the already oriented, validated DNA, never a project SeqId lookup.
    pub fn local_score_target(&self, sequence: &str) -> Result<SequenceScanTarget, String> {
        let sequence = sequence.to_ascii_uppercase();
        if self.geometry.length().ok_or("Invalid TSS geometry")? != sequence.len()
            || sha256_hex_bytes(sequence.as_bytes()) != self.sequence_sha256
        {
            return Err("Local scoring sequence no longer matches the validated TSS window".into());
        }
        Ok(SequenceScanTarget::InlineSequence {
            sequence_text: sequence,
            topology: InlineSequenceTopology::Linear,
            id_hint: Some(self.promoter_id.clone()),
            span_start_0based: None,
            span_end_0based_exclusive: None,
        })
    }

    /// Worker-only orchestration of the existing engine scorer, with exact PFM admission.
    pub fn compute_local_scores(
        &self,
        sequence: &str,
        request: &TssLocalScoreRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<TfbsScoreTrackReport, String> {
        request.validate_budget(sequence.len())?;
        let target = self.local_score_target(sequence)?;
        let registry = crate::tf_motifs::snapshot_db();
        let mut expected = Vec::new();
        for id in &request.matrix_ids {
            let matrix = registry.resolve_exact_full_pfm(id)?;
            if matrix.matrix_counts.len() > 64 {
                return Err(format!(
                    "Matrix {id} exceeds the local 64-column scoring budget"
                ));
            }
            expected.push(TfbsScoreTrackMatrixBinding {
                matrix_id: matrix.id.clone(),
                matrix_name: matrix.name.clone(),
                matrix_sha256: digest(&(&matrix.id, &matrix.name, &matrix.matrix_counts))?,
            });
        }
        // The shared scorer records the actual matrices used. Reject a concurrent
        // registry change instead of attaching provenance from the admission snapshot.
        let report = GentleEngine::default()
            .summarize_tfbs_score_tracks_with_progress(
                target,
                &request.matrix_ids,
                request.score_kind,
                request.clip_negative,
                false,
                progress,
            )
            .map_err(|e| e.to_string())?;
        if report.scoring_provenance.as_ref().map(|p| &p.matrices) != Some(&expected) {
            return Err("Matrix registry changed during local scoring; rerun explicitly".into());
        }
        Ok(report)
    }

    /// Project shared results only after checking their exact input and array shapes.
    pub fn with_local_scores(
        &self,
        request: &TssLocalScoreRequest,
        report: &TfbsScoreTrackReport,
    ) -> Result<Self, String> {
        let length = self.geometry.length().ok_or("Invalid TSS geometry")?;
        request.validate_budget(length)?;
        let provenance = report
            .scoring_provenance
            .as_ref()
            .ok_or("Local scoring provenance unavailable in legacy report")?;
        if report.schema != crate::engine::TFBS_SCORE_TRACK_REPORT_SCHEMA
            || report.target_kind != "inline_sequence"
            || report.target_label != self.promoter_id
            || report.scan_topology != InlineSequenceTopology::Linear
            || report.source_sequence_length_bp != length
            || report.view_start_0based != 0
            || report.view_end_0based_exclusive != length
            || provenance.sequence_sha256 != self.sequence_sha256
            || report.score_kind != request.score_kind
            || report.clip_negative != request.clip_negative
            || report.motifs_requested != request.matrix_ids
            || report.tracks.len() != request.matrix_ids.len()
            || provenance.matrices.len() != report.tracks.len()
        {
            return Err("Local score report does not match the TSS InlineSequence request".into());
        }
        let mut view = self.clone();
        view.clear_local_scores();
        for ((track, binding), id) in report
            .tracks
            .iter()
            .zip(&provenance.matrices)
            .zip(&request.matrix_ids)
        {
            let count = if track.motif_length_bp > length {
                0
            } else {
                length + 1 - track.motif_length_bp
            };
            let validity = track
                .score_validity
                .as_ref()
                .ok_or("Local score validity unavailable")?;
            if &track.tf_id != id
                || binding.matrix_id != *id
                || track.track_start_0based != 0
                || binding.matrix_sha256.len() != 64
                || !binding.matrix_sha256.bytes().all(|c| c.is_ascii_hexdigit())
                || track.motif_length_bp == 0
                || track.scored_window_count != count
                || track.forward_scores.len() != count
                || track.reverse_scores.len() != count
                || validity.forward.len() != count
                || validity.reverse.len() != count
            {
                return Err("Local score matrix identity or array shape is inconsistent".into());
            }
            let forward: Vec<_> = (0..count).map(|i| track.score_at(i, false)).collect();
            let reverse: Vec<_> = (0..count).map(|i| track.score_at(i, true)).collect();
            let valid = forward.iter().chain(&reverse).flatten().count();
            let (min, max) = forward
                .iter()
                .chain(&reverse)
                .flatten()
                .fold((0.0_f64, 0.0_f64), |(lo, hi), &v| (lo.min(v), hi.max(v)));
            view.lanes.push(TssViewLane {
                kind: TssLaneKind::LocalScoreTrace,
                id: format!("local/{id}/{}", report.score_kind.as_str()),
                label: format!("Locally computed | {} | {id}", track.tf_name.as_deref().unwrap_or(id)),
                details: format!("GENtle shared InlineSequence scorer; matrix SHA-256 {}; sequence SHA-256 {}; scorer {}. Motif-window starts in displayed orientation, not measured binding or luciferase activity. Independent from saved report/imported evidence.", binding.matrix_sha256, provenance.sequence_sha256, provenance.scorer),
                units: report.score_kind.as_str().into(),
                state: format!("{valid}/{} evaluated strand-windows; {}", 2 * count,
                    if report.clip_negative { "negative scores clipped" } else { "signed scores retained" }),
                features: vec![],
                trace: Some(TssViewTrace { motif_length_bp: track.motif_length_bp, forward, reverse,
                    clip_negative: report.clip_negative, range_is_fallback: min == max }),
                scale_min: min,
                scale_max: if min == max { min + 1.0 } else { max },
            });
        }
        let source_binding_sha256 = self.local_score_source_binding()?;
        view.local_scoring = Some(TssLocalScoreAttachment {
            request: request.clone(),
            cache_key_sha256: digest(&(&source_binding_sha256, request, provenance))?,
            source_binding_sha256,
            report_sha256: digest(report)?,
            producer_revision: option_env!("GENTLE_SOURCE_REVISION")
                .unwrap_or("unknown")
                .into(),
            provenance: provenance.clone(),
        });
        Ok(view)
    }

    pub fn clear_local_scores(&mut self) {
        self.lanes
            .retain(|l| l.kind != TssLaneKind::LocalScoreTrace);
        self.local_scoring = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Synthetic nine-base annotated window; no private sequence or experimental evidence.
    fn fixture(minus: bool) -> (TssSequenceView, String, TssLocalScoreRequest) {
        let mut seq = gb_io::seq::Seq::empty();
        seq.seq = b"AAAAAAANA".to_vec();
        seq.comments = vec![
            format!(
                "GENtle promoter_id=toy; sequence_sha256={}",
                sha256_hex_bytes(&seq.seq)
            ),
            format!(
                "Reference=toy; assembly=toy1; chromosome=1; genomic=100..108; genomic_strand={}; local_axis=transcript_5prime_to_3prime; TSS_local_1based=5",
                if minus { "-" } else { "+" }
            ),
        ];
        seq.features = vec![gb_io::seq::Feature {
            kind: "misc_feature".into(),
            location: gb_io::seq::Location::single(4),
            qualifiers: vec![
                ("label".into(), Some("Annotated TSS candidate".into())),
                (
                    "note".into(),
                    Some("Genomic 104; annotation-derived".into()),
                ),
            ],
        }];
        let dna = DNAsequence::from_genbank_seq(seq);
        (
            TssSequenceView::from_dna(&dna).unwrap(),
            dna.get_forward_string(),
            TssLocalScoreRequest {
                matrix_ids: vec!["MA0004.1".into()],
                score_kind: TfbsScoreTrackValueKind::LlrBits,
                clip_negative: true,
            },
        )
    }

    #[test]
    fn tss_local_scoring_binds_inline_sequence_matrices_gaps_zero_and_svg() {
        let _guard = crate::tf_motifs::test_registry_lock().lock().unwrap();
        for minus in [false, true] {
            let (view, sequence, request) = fixture(minus);
            let before = serde_json::to_value(&view).unwrap();
            let mut events = 0;
            let report = view
                .compute_local_scores(&sequence, &request, &mut |_| {
                    events += 1;
                    true
                })
                .unwrap();
            assert!(events > 0);
            assert_eq!(serde_json::to_value(&view).unwrap(), before);
            assert_eq!(report.target_kind, "inline_sequence");
            assert!(report.seq_id == view.promoter_id);
            let rendered = view.with_local_scores(&request, &report).unwrap();
            let lane = rendered
                .lanes
                .iter()
                .find(|l| l.kind == TssLaneKind::LocalScoreTrace)
                .unwrap();
            let trace = lane.trace.as_ref().unwrap();
            assert_eq!(trace.forward[0], Some(0.0));
            assert_eq!(trace.forward[2], None);
            assert!(lane.label.starts_with("Locally computed"));
            assert_eq!(lane.units, "llr_bits");
            assert_eq!(rendered.geometry, view.geometry);
            let p = report.scoring_provenance.as_ref().unwrap();
            let matrix = crate::tf_motifs::snapshot_db()
                .resolve_exact_full_pfm("MA0004.1")
                .unwrap();
            assert_eq!(
                p.matrices[0].matrix_sha256,
                digest(&(&matrix.id, &matrix.name, &matrix.matrix_counts)).unwrap()
            );
            assert_eq!(p.sequence_sha256, view.sequence_sha256);
            let options = super::super::TssViewSvgOptions {
                start_0based: 0,
                end_0based_exclusive: 9,
                lane_indices: (0..rendered.lanes.len()).collect(),
                width_px: 1200,
                print_size_mm: None,
            };
            let svg = render_tss_view_svg(&rendered, &options).unwrap();
            assert!(svg.contains("Locally computed"));
            assert!(svg.contains(&rendered.local_scoring.as_ref().unwrap().report_sha256));
            assert!(svg.contains(&p.matrices[0].matrix_sha256));
            assert!(svg.contains("local_scoring"));
            assert!(matches!(
                view.local_score_target(&sequence).unwrap(),
                SequenceScanTarget::InlineSequence {
                    topology: InlineSequenceTopology::Linear,
                    ..
                }
            ));

            let mut wrong = report.clone();
            wrong.scoring_provenance.as_mut().unwrap().sequence_sha256 = "0".repeat(64);
            assert!(view.with_local_scores(&request, &wrong).is_err());
            wrong = report.clone();
            wrong.tracks[0].score_validity = None;
            assert!(view.with_local_scores(&request, &wrong).is_err());
            wrong = report.clone();
            wrong.scoring_provenance = None;
            assert!(view.with_local_scores(&request, &wrong).is_err());
            wrong = report.clone();
            wrong.tracks[0].forward_scores.pop();
            assert!(view.with_local_scores(&request, &wrong).is_err());
            let changed = TssLocalScoreRequest {
                score_kind: TfbsScoreTrackValueKind::TrueLogOddsBits,
                ..request.clone()
            };
            assert!(view.with_local_scores(&changed, &report).is_err());
            assert!(view.local_score_target("AAAAAAAAA").is_err());
            let mut cleared = rendered;
            cleared.clear_local_scores();
            assert_eq!(serde_json::to_value(&cleared).unwrap(), before);
        }
    }

    #[test]
    fn tss_local_scoring_explicit_budgets_exact_ids_and_cancellation() {
        let _guard = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let (view, sequence, request) = fixture(false);
        assert!(request.validate_budget(50_001).is_err());
        let duplicate = TssLocalScoreRequest {
            matrix_ids: vec!["MA0004.1".into(); 2],
            ..request.clone()
        };
        assert!(duplicate.validate_budget(9).is_err());
        let error = view
            .compute_local_scores(&sequence, &request, &mut |_| false)
            .unwrap_err();
        assert!(error.to_lowercase().contains("cancel"));
        for token in ["ALL", "ACGT", "Arnt", "MA99999.999"] {
            let invalid = TssLocalScoreRequest {
                matrix_ids: vec![token.into()],
                ..request.clone()
            };
            assert!(
                view.compute_local_scores(&sequence, &invalid, &mut |_| panic!(
                    "must reject before scoring"
                ))
                .is_err()
            );
        }
    }

    #[test]
    fn tss_local_scoring_report_and_local_attachment_are_independent() {
        let _guard = crate::tf_motifs::test_registry_lock().lock().unwrap();
        let (dna, attached) = crate::tss_sequence_view::profile_fixture(false);
        let view = TssSequenceView::from_dna(&dna).unwrap();
        let request = TssLocalScoreRequest {
            matrix_ids: vec!["MA0004.1".into()],
            score_kind: TfbsScoreTrackValueKind::TrueLogOddsBits,
            clip_negative: false,
        };
        let report = view
            .compute_local_scores(&dna.get_forward_string(), &request, &mut |_| true)
            .unwrap();
        let local = view.with_local_scores(&request, &report).unwrap();
        let mut both = local.with_profile(&attached).unwrap();
        assert_eq!(
            both.local_scoring.as_ref().unwrap().report_sha256,
            digest(&report).unwrap()
        );
        assert!(
            both.lanes
                .iter()
                .any(|l| l.kind == TssLaneKind::ScoreTrace && l.units == "llr_bits")
        );
        assert!(
            both.lanes
                .iter()
                .any(|l| l.kind == TssLaneKind::LocalScoreTrace && l.units == "true_log_odds_bits")
        );
        both.clear_profile();
        assert_eq!(
            serde_json::to_value(&both).unwrap(),
            serde_json::to_value(&local).unwrap()
        );
        let mut both = view
            .with_profile(&attached)
            .unwrap()
            .with_local_scores(&request, &report)
            .unwrap();
        both.clear_local_scores();
        assert_eq!(
            serde_json::to_value(&both).unwrap(),
            serde_json::to_value(view.with_profile(&attached).unwrap()).unwrap()
        );
        let mut relocated = view.clone();
        relocated.assembly = "other-assembly".into();
        assert_ne!(
            view.local_score_source_binding().unwrap(),
            relocated.local_score_source_binding().unwrap()
        );
    }
}
