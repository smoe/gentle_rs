//! Strict TSS documents compose the existing motif scorer without alias expansion.
//!
//! Registry definitions and backgrounds are frozen once per run. Display clipping
//! never changes the raw arrays used for within-factor matrix comparisons.

use super::*;
use crate::digest_utils::sha256_hex_bytes;
use gentle_protocol::tss_profiles::*;

#[path = "tss_profiles_context.rs"]
mod context;

#[cfg(test)]
#[path = "tss_profiles_tests.rs"]
mod tests;

fn invalid(message: impl Into<String>) -> EngineError {
    EngineError::invalid_input(message)
}

fn peaks(
    values: &[Option<f64>],
    limit: usize,
    motif_length: usize,
) -> (Option<TssPeak>, Vec<TssPeak>) {
    let mut ranked = values
        .iter()
        .enumerate()
        .filter_map(|(i, score)| {
            score.map(|score| TssPeak {
                local_start_0based: i,
                score,
            })
        })
        .collect::<Vec<_>>();
    ranked.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then(a.local_start_0based.cmp(&b.local_start_0based))
    });
    let maximum = ranked.first().cloned();
    let mut selected: Vec<TssPeak> = vec![];
    for candidate in ranked {
        if selected.len() == limit {
            break;
        }
        let i = candidate.local_start_0based;
        if i.checked_sub(1)
            .and_then(|j| values[j])
            .is_some_and(|v| v > candidate.score)
            || values
                .get(i + 1)
                .copied()
                .flatten()
                .is_some_and(|v| v > candidate.score)
            || selected
                .iter()
                .any(|p| p.local_start_0based.abs_diff(i) < motif_length)
        {
            continue;
        }
        selected.push(candidate);
    }
    (maximum, selected)
}

impl GentleEngine {
    pub(crate) fn export_tss_profile_report(
        report: &TssProfileReport,
        request: &ExportTssProfilesRequest,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<(Option<TssProfileReport>, TssProfileReceipt), EngineError> {
        let mut enriched = None;
        if let Some(path) = &request.context_manifest {
            let mut report = report.clone();
            let label = report.panel_resolution.panel.label.clone();
            context::attach(&mut report, Path::new(path), &mut || {
                Self::emit_tfbs_score_track_progress(
                    on_progress,
                    &label,
                    "",
                    1,
                    1,
                    1,
                    1,
                    "TSS context",
                    "verifying source hashes and projecting locus evidence",
                    0,
                    1,
                )
            })?;
            enriched = Some(report);
        }
        let report = enriched.as_ref().unwrap_or(report);
        let mut resolved = request.clone();
        resolved.context_manifest = None;
        let receipt = crate::tss_profile_export::export_tss_profiles_with_cancel(
            report,
            &resolved,
            &mut || {
                Self::emit_tfbs_score_track_progress(
                    on_progress,
                    &report.panel_resolution.panel.label,
                    "",
                    1,
                    1,
                    1,
                    1,
                    "document export",
                    "validating, rendering and staging a complete output set",
                    0,
                    1,
                )
            },
        )?;
        Ok((enriched, receipt))
    }

    /// Resolve an immutable panel using exact IDs; no runtime alias/family routing.
    pub(crate) fn resolve_tss_panel(
        panel: JasparTargetPanel,
        panel_bytes: &[u8],
        registry: &tf_motifs::TfMotifDb,
    ) -> Result<TssPanelResolution, EngineError> {
        gentle_engine::tss_profiles::validate_panel(&panel)?;
        if panel.factors.len() > 256 {
            return Err(invalid(
                "TSS profile documents support at most 256 exact matrices per request",
            ));
        }
        let _: TfbsScoreTrackValueKind = serde_json::from_value(json!(panel.score_kind))
            .map_err(|e| invalid(format!("Panel score_kind: {e}")))?;
        let mut matrices = vec![];
        for specification in &panel.factors {
            let motif = registry
                .resolve_exact_full_pfm(&specification.source_id)
                .map_err(invalid)?;
            if motif.matrix_counts.len() > 64
                || motif.consensus_iupac.len() != motif.matrix_counts.len()
            {
                return Err(invalid(format!(
                    "Matrix {}: scoring supports at most 64 columns and requires consistent consensus/PFM shape",
                    motif.id
                )));
            }
            if motif.name.as_deref() != Some(specification.factor_id.as_str()) {
                return Err(invalid(format!(
                    "Matrix {}: factor_id '{}' differs from exact registry name {:?}",
                    specification.source_id, specification.factor_id, motif.name
                )));
            }
            // Bind identity, counts and declared name, not merely a consensus string.
            let bytes = serde_json::to_vec(&(&motif.id, &motif.name, &motif.matrix_counts))
                .map_err(|e| invalid(e.to_string()))?;
            matrices.push(ResolvedTssMatrix {
                specification: specification.clone(),
                version: motif.id.rsplit('.').next().unwrap_or_default().into(),
                consensus: motif.consensus_iupac,
                matrix_counts: motif.matrix_counts,
                matrix_sha256: sha256_hex_bytes(&bytes),
            });
        }
        let (registry_sources, registry_source_url) = registry.strict_source_bindings();
        Ok(TssPanelResolution {
            panel,
            matrices,
            registry_sources,
            registry_source_url,
            panel_sha256: sha256_hex_bytes(panel_bytes),
        })
    }

    pub(crate) fn compute_tss_profiles(
        &self,
        request: &ComputeTssProfilesRequest,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<TssProfileReport, EngineError> {
        let bundle = crate::tss_fasta_bundle::read_bundle(request)?;
        let panel_file = crate::tss_fasta_bundle::open_regular_input(Path::new(&request.panel))?;
        let mut panel_bytes = Vec::new();
        panel_file
            .take(2 * 1024 * 1024 + 1)
            .read_to_end(&mut panel_bytes)
            .map_err(|e| invalid(format!("Panel read: {e}")))?;
        if panel_bytes.len() > 2 * 1024 * 1024 {
            return Err(invalid("Panel exceeds 2 MiB"));
        }
        let panel = crate::tfbs_track_panel::parse_tss_panel(&panel_bytes)?;
        let registry = tf_motifs::snapshot_db();
        let resolution = Self::resolve_tss_panel(panel, &panel_bytes, &registry)?;
        let mut report = self.compute_verified_tss_profiles(bundle, resolution, on_progress)?;
        let executable = std::env::current_exe()
            .map_err(|e| invalid(format!("Cannot identify profile producer executable: {e}")))?;
        report.producer_executable_sha256 = Some(
            crate::digest_utils::sha256_file_hex(&executable)
                .map_err(|e| invalid(format!("Cannot bind profile producer executable: {e}")))?,
        );
        Ok(report)
    }

    fn compute_verified_tss_profiles(
        &self,
        mut bundle: crate::tss_fasta_bundle::VerifiedTssBundle,
        resolution: TssPanelResolution,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<TssProfileReport, EngineError> {
        if bundle.records.len() > 4096 {
            return Err(invalid(
                "TSS profile documents support at most 4096 windows; split the input bundle",
            ));
        }
        let points = bundle
            .records
            .iter()
            .try_fold(0usize, |n, (_, seq, _)| {
                n.checked_add(
                    seq.len()
                        .checked_mul(resolution.matrices.len())?
                        .checked_mul(2)?,
                )
            })
            .ok_or_else(|| invalid("Profile size overflows"))?;
        if points > 10_000_000 {
            return Err(invalid(
                "Profile request exceeds 10 million strand/position scores; split the input bundle",
            ));
        }
        bundle.records.sort_by(|a, b| {
            a.0.gene_id
                .cmp(&b.0.gene_id)
                .then(b.2.cmp(&a.2))
                .then(a.0.geometry.chromosome.cmp(&b.0.geometry.chromosome))
                .then(a.0.geometry.start_1based.cmp(&b.0.geometry.start_1based))
                .then(a.0.promoter_id.cmp(&b.0.promoter_id))
        });
        let score_kind: TfbsScoreTrackValueKind =
            serde_json::from_value(json!(resolution.panel.score_kind))
                .map_err(|e| invalid(e.to_string()))?;
        let mut windows = bundle
            .records
            .iter()
            .map(|(record, _, selected)| TssProfileWindow {
                detail_context: None,
                record: record.clone(),
                selected: *selected,
                selection_evidence: bundle.selection_evidence.get(&record.promoter_id).cloned(),
                tracks: vec![],
                comparisons: vec![],
            })
            .collect::<Vec<_>>();
        let random = Self::deterministic_random_dna_bytes(
            DEFAULT_TFBS_SCORE_TRACK_RANDOM_SEQUENCE_LENGTH_BP,
            DEFAULT_TFBS_SCORE_TRACK_RANDOM_SEED,
        );
        let stages = windows.len() + 1;
        for (matrix_index, matrix) in resolution.matrices.iter().enumerate() {
            let accession = &matrix.specification.source_id;
            let emit = |callback: &mut dyn FnMut(OperationProgress) -> bool,
                        stage: usize,
                        count: usize,
                        total: usize| {
                Self::emit_tfbs_score_track_progress(
                    callback,
                    &resolution.panel.label,
                    accession,
                    matrix_index + 1,
                    resolution.matrices.len(),
                    stage,
                    stages,
                    if stage == 1 {
                        "background calibration"
                    } else {
                        "TSS scan"
                    },
                    "frozen accession; unclipped computational scores",
                    count,
                    total,
                )
            };
            if !emit(on_progress, 1, 0, 1) {
                return Err(Self::tfbs_cancelled_error("background calibration"));
            }
            let (llr, lor) = Self::prepare_scoring_matrices(&matrix.matrix_counts);
            let llr_model = Self::modeled_tfbs_score_distribution(&llr);
            let lor_model = Self::modeled_tfbs_score_distribution(&lor);
            let (_, background, llr_sorted, lor_sorted) =
                Self::collect_tfbs_score_track_background_scores(
                    &random,
                    &llr,
                    &lor,
                    score_kind,
                    false,
                    llr_model.as_ref(),
                    lor_model.as_ref(),
                    |n, t| emit(on_progress, 1, n, t),
                )?;
            for (window_index, ((_, sequence, _), window)) in
                bundle.records.iter().zip(windows.iter_mut()).enumerate()
            {
                let length = matrix.matrix_counts.len();
                let count = sequence
                    .len()
                    .checked_sub(length)
                    .map(|n| n + 1)
                    .unwrap_or(0);
                let mut forward = vec![None; count];
                let mut reverse = vec![None; count];
                if !emit(on_progress, window_index + 2, 0, 1) {
                    return Err(Self::tfbs_cancelled_error("TSS scan"));
                }
                let hits = Self::scan_tf_scores_with_topology_and_cancel(
                    sequence.as_bytes(),
                    &llr,
                    &lor,
                    InlineSequenceTopology::Linear,
                    |n, t| emit(on_progress, window_index + 2, n, t),
                )?;
                let mut max_underlying: Option<f64> = None;
                for (offset, is_reverse, bits, quantile, odds, odds_quantile) in hits {
                    let value = Self::tfbs_score_track_presented_value(
                        bits,
                        quantile,
                        odds,
                        odds_quantile,
                        score_kind,
                        false,
                        llr_model.as_ref(),
                        lor_model.as_ref(),
                        &llr_sorted,
                        &lor_sorted,
                    );
                    if !value.is_finite() {
                        return Err(invalid(format!("Matrix {accession}: nonfinite score")));
                    }
                    if is_reverse {
                        reverse[offset] = Some(value);
                    } else {
                        forward[offset] = Some(value);
                    }
                    let underlying = if score_kind.uses_llr_background_bits() {
                        bits
                    } else {
                        odds
                    };
                    max_underlying = Some(max_underlying.map_or(underlying, |v| v.max(underlying)));
                }
                let normalization_reference = match max_underlying {
                    Some(peak) => serde_json::to_value(
                        Self::summarize_tfbs_score_track_normalization_reference(
                            &background,
                            peak,
                            if score_kind.uses_llr_background_bits() {
                                llr_model.as_ref()
                            } else {
                                lor_model.as_ref()
                            },
                        ),
                    )
                    .map_err(|e| invalid(e.to_string()))?,
                    None => serde_json::Value::Null,
                };
                let (forward_maximum, forward_peaks) =
                    peaks(&forward, resolution.panel.top_hit_count, length);
                let (reverse_maximum, reverse_peaks) =
                    peaks(&reverse, resolution.panel.top_hit_count, length);
                window.tracks.push(TssProfileTrack {
                    accession: accession.clone(),
                    motif_length_bp: length,
                    forward_scores: forward,
                    reverse_scores: reverse,
                    forward_maximum,
                    reverse_maximum,
                    forward_peaks,
                    reverse_peaks,
                    normalization_reference,
                });
            }
        }
        for window in &mut windows {
            window.comparisons = Self::tss_matrix_comparisons(window, &resolution);
        }
        let score_policy = BTreeMap::from([
            ("background".into(), "uniform_iid_ACGT_0.25".into()),
            (
                "pseudocount".into(),
                "pad_columns_to_max_total_then_add_max_total_times_1e-9_per_base".into(),
            ),
            (
                "modeled_quantum_bits".into(),
                Self::TFBS_MODELED_SCORE_QUANTUM_BITS.to_string(),
            ),
            ("background_tail_show_quantile".into(), "0.95".into()),
            (
                "background_seed".into(),
                DEFAULT_TFBS_SCORE_TRACK_RANDOM_SEED.to_string(),
            ),
            (
                "background_length_bp".into(),
                DEFAULT_TFBS_SCORE_TRACK_RANDOM_SEQUENCE_LENGTH_BP.to_string(),
            ),
            (
                "correlation_signal".into(),
                "unclipped_declared_score_kind; no_smoothing; common_valid_window_starts".into(),
            ),
            (
                "ambiguous_windows".into(),
                "unavailable_null_not_zero".into(),
            ),
        ]);
        bundle.inputs.push(TssInputBinding {
            role: "panel".into(),
            name: "panel.json".into(),
            sha256: resolution.panel_sha256.clone(),
        });
        Ok(TssProfileReport {
            schema: REPORT_SCHEMA.into(),
            reference: bundle.reference,
            panel_resolution: resolution,
            inputs: bundle.inputs,
            windows,
            source: Some(bundle.source),
            producer_executable_sha256: None,
            producer_revision: option_env!("GENTLE_SOURCE_REVISION")
                .unwrap_or("unknown")
                .into(),
            lockfile_sha256: sha256_hex_bytes(include_bytes!("../../../Cargo.lock")),
            score_policy,
            verification: "bundle_consistency_verified; prepared_reference_not_assessed".into(),
            warnings: bundle.warnings,
            non_claims: NON_CLAIMS.into(),
        })
    }

    fn tss_matrix_comparisons(
        window: &TssProfileWindow,
        resolution: &TssPanelResolution,
    ) -> Vec<TssMatrixComparison> {
        // Pearson is scale invariant. Normalize here so the legacy helper's
        // absolute variance cutoff cannot mislabel a small nonconstant signal.
        let normalized = |values: &[f64]| {
            let magnitude = values.iter().map(|v| v.abs()).fold(0.0, f64::max);
            let scaled = values
                .iter()
                .map(|value| value / magnitude)
                .collect::<Vec<_>>();
            let min = scaled.iter().copied().fold(f64::INFINITY, f64::min);
            let max = scaled.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            scaled
                .iter()
                .map(|value| (value - min) / (max - min))
                .collect::<Vec<_>>()
        };
        let mut rows = vec![];
        for (i, left) in window.tracks.iter().enumerate() {
            for (j, right) in window.tracks.iter().enumerate().skip(i + 1) {
                let factor = &resolution.matrices[i].specification.factor_id;
                if *factor != resolution.matrices[j].specification.factor_id {
                    continue;
                }
                for (strand, a, b) in [
                    (TssStrand::Plus, &left.forward_scores, &right.forward_scores),
                    (
                        TssStrand::Minus,
                        &left.reverse_scores,
                        &right.reverse_scores,
                    ),
                ] {
                    let paired = a
                        .iter()
                        .zip(b)
                        .filter_map(|(a, b)| Some(((*a)?, (*b)?)))
                        .collect::<Vec<_>>();
                    let av = paired.iter().map(|p| p.0).collect::<Vec<_>>();
                    let bv = paired.iter().map(|p| p.1).collect::<Vec<_>>();
                    let undefined_reason = if paired.len() < 2 {
                        Some("insufficient_common_valid_windows".into())
                    } else if av.iter().all(|v| *v == av[0]) || bv.iter().all(|v| *v == bv[0]) {
                        Some("constant_signal".into())
                    } else {
                        None
                    };
                    rows.push(TssMatrixComparison {factor_id:factor.clone(),left_accession:left.accession.clone(),
                        right_accession:right.accession.clone(),local_strand:strand,paired_window_count:paired.len(),
                        excluded_window_count:a.len().max(b.len())-paired.len(),
                        pearson:undefined_reason.is_none().then(||Self::pearson_correlation(&normalized(&av),&normalized(&bv))),
                        spearman:undefined_reason.is_none().then(||Self::spearman_correlation(&av,&bv)),
                        undefined_reason,method:"Pearson and Spearman (average ranks for ties); same local strand; unclipped declared score kind; no smoothing; common valid window starts".into()});
                }
            }
        }
        rows
    }
}
