//! Deterministic structural screening for potential cryptic splice removal.
//!
//! This screen reports sequence architecture and heuristic support. It does
//! not claim that a candidate is used by the spliceosome. Optional statistical
//! splice-site models can populate the reserved score fields without changing
//! structural candidates into observed junction evidence.

use super::*;
use crate::{AMINO_ACIDS, amino_acids::STOP_CODON};
use gentle_protocol::{
    CRYPTIC_SPLICING_EVIDENCE_OVERLAY_SCHEMA, CRYPTIC_SPLICING_PROTEIN_PROJECTION_SCHEMA,
    CRYPTIC_SPLICING_SCREEN_SCHEMA, CrypticSplicingBranchpointSignal, CrypticSplicingBudgetSummary,
    CrypticSplicingCandidateEvidenceOverlayRow, CrypticSplicingCandidateProteinProjectionRow,
    CrypticSplicingCandidateRow, CrypticSplicingCdsConsequence,
    CrypticSplicingEvidenceBindingStatus, CrypticSplicingEvidenceClass,
    CrypticSplicingEvidenceOverlayReport, CrypticSplicingEvidenceOverlayRequest,
    CrypticSplicingGenomicProvenanceRow, CrypticSplicingModelPolicy,
    CrypticSplicingModelProvenance, CrypticSplicingModelStatus,
    CrypticSplicingPolypyrimidineSignal, CrypticSplicingProjectedProteinFeatureRow,
    CrypticSplicingProteinProjectionReport, CrypticSplicingProteinProjectionRequest,
    CrypticSplicingProteinProjectionStatus, CrypticSplicingRnaEvidenceStatus,
    CrypticSplicingRnaJunctionMatch, CrypticSplicingRnaReportBinding, CrypticSplicingScreenRequest,
    CrypticSplicingScreenView, CrypticSplicingSignalStatus, CrypticSplicingSiteKind,
    CrypticSplicingSiteRow, CrypticSplicingSpan, CrypticSplicingStrand,
};

#[derive(Clone)]
struct CrypticCdsContext {
    feature_id: usize,
    strand_matches: bool,
    source_positions_0based: Vec<usize>,
    coding_sequence: Vec<u8>,
    translation_table: usize,
}

#[derive(Clone, Copy, Debug)]
struct CrypticRankedPair {
    donor_index: usize,
    acceptor_index: usize,
    score: Option<f64>,
}

impl PartialEq for CrypticRankedPair {
    fn eq(&self, other: &Self) -> bool {
        self.donor_index == other.donor_index
            && self.acceptor_index == other.acceptor_index
            && self.score.map(f64::to_bits) == other.score.map(f64::to_bits)
    }
}

impl Eq for CrypticRankedPair {}

impl PartialOrd for CrypticRankedPair {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CrypticRankedPair {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // BinaryHeap keeps the greatest value at its root. Here "greater"
        // deliberately means worse so the current worst top-K row is replaced
        // in O(log K) time. Ascending sort consequently yields best-first rows.
        match (self.score, other.score) {
            (Some(left), Some(right)) => right.total_cmp(&left),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        }
        .then_with(|| self.donor_index.cmp(&other.donor_index))
        .then_with(|| self.acceptor_index.cmp(&other.acceptor_index))
    }
}

impl GentleEngine {
    fn cryptic_local_to_source_1based(
        request: &CrypticSplicingScreenRequest,
        local_position_0based: usize,
    ) -> usize {
        match request.strand {
            CrypticSplicingStrand::Forward => request.start_1based + local_position_0based,
            CrypticSplicingStrand::Reverse => {
                request.end_1based.saturating_sub(local_position_0based)
            }
        }
    }

    fn cryptic_source_to_local_1based(
        request: &CrypticSplicingScreenRequest,
        source_position_1based: usize,
    ) -> usize {
        match request.strand {
            CrypticSplicingStrand::Forward => source_position_1based
                .saturating_sub(request.start_1based)
                .saturating_add(1),
            CrypticSplicingStrand::Reverse => request
                .end_1based
                .saturating_sub(source_position_1based)
                .saturating_add(1),
        }
    }

    fn cryptic_site_row(
        request: &CrypticSplicingScreenRequest,
        source_digest: &str,
        oriented_sequence: &str,
        local_position_0based: usize,
        kind: CrypticSplicingSiteKind,
    ) -> CrypticSplicingSiteRow {
        let (kind_label, motif, window_range) = match kind {
            CrypticSplicingSiteKind::Donor => (
                "donor",
                "GT",
                local_position_0based
                    .checked_sub(3)
                    .map(|start| (start, local_position_0based.saturating_add(6))),
            ),
            CrypticSplicingSiteKind::Acceptor => (
                "acceptor",
                "AG",
                local_position_0based
                    .checked_sub(18)
                    .map(|start| (start, local_position_0based.saturating_add(5))),
            ),
        };
        let model_window = window_range
            .filter(|(start, end)| *start < *end && *end <= oriented_sequence.len())
            .map(|(start, end)| oriented_sequence[start..end].to_string());
        let not_evaluable_reason = model_window.is_none().then(|| {
            format!(
                "The {kind_label} is too close to the scanned-span boundary for the reserved fixed-width MaxEnt window."
            )
        });
        let scanned_position_1based = local_position_0based + 1;
        let identity = format!(
            "{source_digest}|{}|{kind_label}|{scanned_position_1based}|{motif}",
            request.strand.as_str()
        );
        CrypticSplicingSiteRow {
            site_id: short_sha256_id("cryptic_site", &identity),
            evidence_class: CrypticSplicingEvidenceClass::StructuralCandidate,
            kind,
            motif_2bp: motif.to_string(),
            scanned_position_1based,
            source_position_1based: Self::cryptic_local_to_source_1based(
                request,
                local_position_0based,
            ),
            model_window,
            model_status: CrypticSplicingModelStatus::Absent,
            maxent_score: None,
            not_evaluable_reason,
        }
    }

    fn cryptic_branchpoint_signal(
        request: &CrypticSplicingScreenRequest,
        donor_0based: usize,
        oriented_candidate: &str,
    ) -> CrypticSplicingBranchpointSignal {
        if oriented_candidate.len() < 20 {
            return CrypticSplicingBranchpointSignal {
                status: CrypticSplicingSignalStatus::NotEvaluable,
                not_evaluable_reason: Some(
                    "The candidate is shorter than the 20 nt minimum needed for the acceptor-proximal branchpoint heuristic window."
                        .to_string(),
                ),
                ..CrypticSplicingBranchpointSignal::default()
            };
        }
        match Self::detect_branchpoint_like_site(oriented_candidate) {
            Some((offset_0based, motif, score, annotation)) => {
                let local_position_0based = donor_0based + offset_0based;
                CrypticSplicingBranchpointSignal {
                    status: CrypticSplicingSignalStatus::Detected,
                    scanned_position_1based: Some(local_position_0based + 1),
                    source_position_1based: Some(Self::cryptic_local_to_source_1based(
                        request,
                        local_position_0based,
                    )),
                    motif: Some(motif),
                    heuristic_score: Some(score),
                    annotation: Some(annotation.to_string()),
                    not_evaluable_reason: None,
                }
            }
            None => CrypticSplicingBranchpointSignal {
                status: CrypticSplicingSignalStatus::NotDetected,
                annotation: Some(
                    "No branchpoint-like adenine was detected by the bounded acceptor-proximal heuristic."
                        .to_string(),
                ),
                ..CrypticSplicingBranchpointSignal::default()
            },
        }
    }

    fn cryptic_polypyrimidine_signal(
        request: &CrypticSplicingScreenRequest,
        donor_0based: usize,
        oriented_candidate: &str,
    ) -> CrypticSplicingPolypyrimidineSignal {
        if oriented_candidate.len() < 12 {
            return CrypticSplicingPolypyrimidineSignal {
                status: CrypticSplicingSignalStatus::NotEvaluable,
                not_evaluable_reason: Some(
                    "The candidate is shorter than the 12 nt minimum needed for the polypyrimidine-tract heuristic."
                        .to_string(),
                ),
                ..CrypticSplicingPolypyrimidineSignal::default()
            };
        }
        match Self::detect_polypyrimidine_tract(oriented_candidate) {
            Some((start_0based, end_0based_exclusive, fraction, annotation)) => {
                let global_start_0based = donor_0based + start_0based;
                let global_end_0based = donor_0based + end_0based_exclusive - 1;
                let source_left = Self::cryptic_local_to_source_1based(
                    request,
                    global_start_0based,
                );
                let source_right =
                    Self::cryptic_local_to_source_1based(request, global_end_0based);
                CrypticSplicingPolypyrimidineSignal {
                    status: CrypticSplicingSignalStatus::Detected,
                    scanned_start_1based: Some(global_start_0based + 1),
                    scanned_end_1based: Some(global_end_0based + 1),
                    source_start_1based: Some(source_left.min(source_right)),
                    source_end_1based: Some(source_left.max(source_right)),
                    pyrimidine_fraction: Some(fraction),
                    annotation: Some(annotation.to_string()),
                    not_evaluable_reason: None,
                }
            }
            None => CrypticSplicingPolypyrimidineSignal {
                status: CrypticSplicingSignalStatus::NotDetected,
                annotation: Some(
                    "No acceptor-proximal polypyrimidine-rich tract met the current density heuristic."
                        .to_string(),
                ),
                ..CrypticSplicingPolypyrimidineSignal::default()
            },
        }
    }

    fn cryptic_boundary_class(
        request: &CrypticSplicingScreenRequest,
        donor_source_1based: usize,
        acceptor_source_1based: usize,
    ) -> String {
        let Some(insert) = request.insert_span.as_ref() else {
            return "insert_span_not_specified".to_string();
        };
        let compartment = |position: usize| {
            if position >= insert.start_1based && position <= insert.end_1based {
                "insert"
            } else {
                "vector"
            }
        };
        format!(
            "{}_to_{}",
            compartment(donor_source_1based),
            compartment(acceptor_source_1based)
        )
    }

    fn cryptic_cds_context(
        dna: &DNAsequence,
        request: &CrypticSplicingScreenRequest,
    ) -> Result<Option<CrypticCdsContext>, EngineError> {
        let Some(feature_id) = request.cds_feature_id else {
            return Ok(None);
        };
        let feature = dna.features().get(feature_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "CDS feature id '{}' was not found in sequence '{}'",
                feature_id, request.seq_id
            ),
            cause_chain: vec![],
        })?;
        if !feature.kind.eq_ignore_ascii_case("CDS") {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Feature id '{}' in sequence '{}' is '{}' rather than CDS",
                    feature_id, request.seq_id, feature.kind
                ),
                cause_chain: vec![],
            });
        }
        let feature_reverse = feature_is_reverse(feature);
        let request_reverse = request.strand == CrypticSplicingStrand::Reverse;
        let mut ranges = vec![];
        collect_location_ranges_usize(&feature.location, &mut ranges);
        ranges.sort_unstable_by_key(|range| range.0);
        if feature_reverse {
            ranges.reverse();
        }
        let mut source_positions_0based = vec![];
        let mut coding_sequence = vec![];
        for (start, end) in ranges {
            if end <= start || end > dna.len() {
                continue;
            }
            if feature_reverse {
                source_positions_0based.extend((start..end).rev());
                coding_sequence.extend(Self::reverse_complement_bytes(
                    &dna.forward_bytes()[start..end],
                ));
            } else {
                source_positions_0based.extend(start..end);
                coding_sequence.extend_from_slice(&dna.forward_bytes()[start..end]);
            }
        }
        let translation_table = Self::feature_qualifier_text(feature, "transl_table")
            .and_then(|raw| raw.parse::<usize>().ok())
            .unwrap_or(1);
        Ok(Some(CrypticCdsContext {
            feature_id,
            strand_matches: feature_reverse == request_reverse,
            source_positions_0based,
            coding_sequence,
            translation_table,
        }))
    }

    fn cryptic_translated_length(sequence: &[u8], translation_table: usize) -> usize {
        let mut amino_acid_count = 0usize;
        for codon in sequence.chunks_exact(3) {
            let aa = AMINO_ACIDS.codon2aa([codon[0], codon[1], codon[2]], Some(translation_table));
            if aa == STOP_CODON {
                break;
            }
            amino_acid_count += 1;
        }
        amino_acid_count
    }

    fn cryptic_cds_consequence(
        context: Option<&CrypticCdsContext>,
        removed_source_start_1based: usize,
        removed_source_end_1based: usize,
    ) -> Option<CrypticSplicingCdsConsequence> {
        let context = context?;
        if !context.strand_matches {
            return Some(CrypticSplicingCdsConsequence {
                status: "strand_mismatch_not_evaluable".to_string(),
                cds_feature_id: Some(context.feature_id),
                interpretation: "The requested scan strand and CDS feature strand differ; coding consequences were not inferred."
                    .to_string(),
                ..CrypticSplicingCdsConsequence::default()
            });
        }
        let removed_start_0based = removed_source_start_1based.saturating_sub(1);
        let removed_end_0based_exclusive = removed_source_end_1based;
        let removed_indices = context
            .source_positions_0based
            .iter()
            .enumerate()
            .filter_map(|(coding_index, source_position)| {
                (*source_position >= removed_start_0based
                    && *source_position < removed_end_0based_exclusive)
                    .then_some(coding_index)
            })
            .collect::<Vec<_>>();
        if removed_indices.is_empty() {
            return Some(CrypticSplicingCdsConsequence {
                status: "no_coding_overlap".to_string(),
                cds_feature_id: Some(context.feature_id),
                interpretation: "The candidate removal does not overlap the selected CDS."
                    .to_string(),
                ..CrypticSplicingCdsConsequence::default()
            });
        }
        let first = removed_indices[0];
        let last = *removed_indices.last().unwrap_or(&first);
        let removed_coding_bp = removed_indices.len();
        let removed_set = removed_indices.iter().copied().collect::<HashSet<_>>();
        let altered_coding_sequence = context
            .coding_sequence
            .iter()
            .enumerate()
            .filter_map(|(index, base)| (!removed_set.contains(&index)).then_some(*base))
            .collect::<Vec<_>>();
        let native_terminal_codon_is_stop = context
            .coding_sequence
            .chunks_exact(3)
            .next_back()
            .map(|codon| {
                AMINO_ACIDS.codon2aa(
                    [codon[0], codon[1], codon[2]],
                    Some(context.translation_table),
                ) == STOP_CODON
            })
            .unwrap_or(false);
        let terminal_codon_start = context.coding_sequence.len().saturating_sub(3);
        let native_stop_removed = native_terminal_codon_is_stop.then_some(
            removed_indices
                .iter()
                .any(|index| *index >= terminal_codon_start),
        );
        let frame_delta = removed_coding_bp % 3;
        let interpretation = if frame_delta == 0 {
            "The candidate removes coding bases in-frame; this is a sequence consequence, not evidence that splicing occurs."
        } else {
            "The candidate removes a non-multiple of three coding bases and would alter the downstream reading frame if used."
        };
        Some(CrypticSplicingCdsConsequence {
            status: if frame_delta == 0 {
                "coding_in_frame_removal"
            } else {
                "coding_frameshift"
            }
            .to_string(),
            cds_feature_id: Some(context.feature_id),
            removed_coding_bp,
            frame_delta_mod_3: Some(frame_delta),
            affected_coding_interval_1based: Some(CrypticSplicingSpan {
                start_1based: first + 1,
                end_1based: last + 1,
            }),
            affected_aa_interval: None,
            native_stop_removed,
            first_altered_aa_position: Some(first / 3 + 1),
            predicted_protein_length_aa: Some(Self::cryptic_translated_length(
                &altered_coding_sequence,
                context.translation_table,
            )),
            interpretation: interpretation.to_string(),
        })
    }

    fn cryptic_genomic_provenance(
        dna: &DNAsequence,
        request: &CrypticSplicingScreenRequest,
    ) -> Vec<CrypticSplicingGenomicProvenanceRow> {
        let scan_start_0based = request.start_1based.saturating_sub(1);
        let scan_end_0based_exclusive = request.end_1based;
        let mut rows = vec![];
        for feature in dna.features() {
            if !feature.kind.eq_ignore_ascii_case("exon")
                || Self::feature_qualifier_text(feature, "synthetic_origin").as_deref()
                    != Some("mrna_transcript_derived")
            {
                continue;
            }
            let mut ranges = vec![];
            collect_location_ranges_usize(&feature.location, &mut ranges);
            for (feature_start, feature_end) in ranges {
                let Some((overlap_start, overlap_end)) = Self::range_intersection_0based(
                    (feature_start, feature_end),
                    (scan_start_0based, scan_end_0based_exclusive),
                ) else {
                    continue;
                };
                let source_left = overlap_start + 1;
                let source_right = overlap_end;
                let local_left = Self::cryptic_source_to_local_1based(request, source_left);
                let local_right = Self::cryptic_source_to_local_1based(request, source_right);
                rows.push(CrypticSplicingGenomicProvenanceRow {
                    scanned_start_1based: local_left.min(local_right),
                    scanned_end_1based: local_left.max(local_right),
                    source_seq_id: Self::feature_qualifier_text(feature, "source_seq_id"),
                    source_feature_id: Self::feature_qualifier_text(feature, "source_feature_id")
                        .and_then(|raw| raw.parse::<usize>().ok()),
                    genomic_start_1based: Self::feature_qualifier_text(
                        feature,
                        "genomic_start_1based",
                    )
                    .and_then(|raw| raw.parse::<usize>().ok()),
                    genomic_end_1based: Self::feature_qualifier_text(feature, "genomic_end_1based")
                        .and_then(|raw| raw.parse::<usize>().ok()),
                    strand: Self::feature_qualifier_text(feature, "strand"),
                });
            }
        }
        rows.sort_by(|left, right| {
            left.scanned_start_1based
                .cmp(&right.scanned_start_1based)
                .then(left.scanned_end_1based.cmp(&right.scanned_end_1based))
        });
        rows
    }

    pub fn inspect_cryptic_splicing_screen(
        &self,
        request: &CrypticSplicingScreenRequest,
    ) -> Result<CrypticSplicingScreenView, EngineError> {
        let dna = self
            .state
            .sequences
            .get(&request.seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", request.seq_id),
                cause_chain: vec![],
            })?;
        if request.start_1based == 0 || request.end_1based == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Cryptic-splicing scan coordinates are 1-based and must be greater than zero."
                        .to_string(),
                cause_chain: vec![],
            });
        }
        if request.start_1based > request.end_1based {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Cryptic-splicing v1 does not support a circular span that wraps across the sequence origin; materialize a linear span first."
                    .to_string(),
                cause_chain: vec![],
            });
        }
        if request.end_1based > dna.len() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Cryptic-splicing scan end {} exceeds sequence '{}' length {}",
                    request.end_1based,
                    request.seq_id,
                    dna.len()
                ),
                cause_chain: vec![],
            });
        }
        if request.min_pseudo_intron_bp < 4
            || request.max_pseudo_intron_bp < request.min_pseudo_intron_bp
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Cryptic-splicing pseudo-intron bounds require min >= 4 and max >= min."
                    .to_string(),
                cause_chain: vec![],
            });
        }
        if request.max_candidate_pairs == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Cryptic-splicing max_candidate_pairs must be greater than zero."
                    .to_string(),
                cause_chain: vec![],
            });
        }
        if request.model_policy == CrypticSplicingModelPolicy::StructuralOnly
            && request.expected_model_fingerprint_sha256.is_some()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "expected_model_fingerprint_sha256 requires model_policy=use_active_maxent"
                        .to_string(),
                cause_chain: vec![],
            });
        }
        if let Some(insert) = request.insert_span.as_ref()
            && (insert.start_1based == 0
                || insert.end_1based < insert.start_1based
                || insert.end_1based > dna.len())
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Cryptic-splicing insert_span must be a non-wrapping 1-based inclusive span on the source sequence."
                    .to_string(),
                cause_chain: vec![],
            });
        }

        let raw = &dna.forward_bytes()[request.start_1based.saturating_sub(1)..request.end_1based];
        let mut normalized = raw
            .iter()
            .map(|base| match base.to_ascii_uppercase() {
                b'U' => b'T',
                other => other,
            })
            .collect::<Vec<_>>();
        if request.strand == CrypticSplicingStrand::Reverse {
            normalized = Self::reverse_complement_bytes(&normalized);
        }
        let oriented_sequence = String::from_utf8_lossy(&normalized).to_string();
        let source_digest = sha256_prefixed_bytes(&normalized);
        let request_sha256 = serde_json::to_vec(request)
            .map(|bytes| sha256_prefixed_bytes(&bytes))
            .map_err(|error| EngineError {
                code: ErrorCode::Internal,
                message: format!("Could not fingerprint cryptic-splicing request: {error}"),
                cause_chain: vec![],
            })?;

        let mut model_provenance = CrypticSplicingModelProvenance {
            status: CrypticSplicingModelStatus::Absent,
            note: "The request selected the dependency-free structural-only policy; no optional statistical splice-site model was applied."
                .to_string(),
            ..CrypticSplicingModelProvenance::default()
        };
        let mut active_model = None;
        if request.model_policy == CrypticSplicingModelPolicy::UseActiveMaxent {
            let status = crate::maxent_splicing::active_model_status();
            model_provenance.resource_path = Some(status.path.clone());
            match status.model {
                Some(model) => {
                    let snapshot = model.snapshot();
                    let fingerprint = snapshot.model_fingerprint_sha256.clone();
                    if request
                        .expected_model_fingerprint_sha256
                        .as_deref()
                        .is_some_and(|expected| expected != fingerprint)
                    {
                        model_provenance.status = CrypticSplicingModelStatus::Failed;
                        model_provenance.resource_sha256 = Some(fingerprint.clone());
                        model_provenance.note = format!(
                            "The active MaxEnt model fingerprint '{}' does not match the request lock '{}'; structural results are retained without model scores.",
                            fingerprint,
                            request
                                .expected_model_fingerprint_sha256
                                .as_deref()
                                .unwrap_or_default()
                        );
                    } else {
                        model_provenance = CrypticSplicingModelProvenance {
                            status: CrypticSplicingModelStatus::Present,
                            model_kind: Some("maxentscan_native_yeo_burge_2004".to_string()),
                            resource_path: Some(status.path),
                            resource_sha256: Some(fingerprint),
                            source_url: snapshot.source_url.clone(),
                            retrieved_on: snapshot.retrieved_on.clone(),
                            redistribution_status: Some(
                                snapshot.redistribution_status.clone(),
                            ),
                            table_sha256: snapshot.table_sha256.clone(),
                            scoring_implementation: Some(
                                "GENtle native reproduction of MaxEntScan score5.pl/score3.pl"
                                    .to_string(),
                            ),
                            note: "Applied the active user-supplied MaxEntScan probability tables. Scores prioritize sequence-model compatibility; they do not establish observed splicing."
                                .to_string(),
                        };
                        active_model = Some(model);
                    }
                }
                None if status.exists => {
                    model_provenance.status = CrypticSplicingModelStatus::Failed;
                    model_provenance.note = format!(
                        "The requested MaxEnt resource is present but invalid: {}. Structural results are retained without model scores.",
                        status
                            .error
                            .unwrap_or_else(|| "unknown validation error".to_string())
                    );
                }
                None => {
                    model_provenance.note = format!(
                        "No user-supplied MaxEnt snapshot is available at '{}'; structural results are retained and model scores are unavailable.",
                        status.path
                    );
                }
            }
        }
        let effective_input_sha256 = sha256_prefixed_bytes(
            format!(
                "{}|{}|{}|{}",
                request_sha256,
                source_digest,
                match model_provenance.status {
                    CrypticSplicingModelStatus::Absent => "absent",
                    CrypticSplicingModelStatus::Present => "present",
                    CrypticSplicingModelStatus::Failed => "failed",
                },
                model_provenance.resource_sha256.as_deref().unwrap_or("-")
            )
            .as_bytes(),
        );

        let donor_positions = normalized
            .windows(2)
            .enumerate()
            .filter_map(|(index, motif)| (motif == b"GT").then_some(index))
            .collect::<Vec<_>>();
        let acceptor_positions = normalized
            .windows(2)
            .enumerate()
            .filter_map(|(index, motif)| (motif == b"AG").then_some(index))
            .collect::<Vec<_>>();
        let mut donor_sites = donor_positions
            .iter()
            .map(|position| {
                Self::cryptic_site_row(
                    request,
                    &source_digest,
                    &oriented_sequence,
                    *position,
                    CrypticSplicingSiteKind::Donor,
                )
            })
            .collect::<Vec<_>>();
        let mut acceptor_sites = acceptor_positions
            .iter()
            .map(|position| {
                Self::cryptic_site_row(
                    request,
                    &source_digest,
                    &oriented_sequence,
                    *position,
                    CrypticSplicingSiteKind::Acceptor,
                )
            })
            .collect::<Vec<_>>();
        if let Some(model) = active_model.as_ref() {
            for site in &mut donor_sites {
                match site
                    .model_window
                    .as_deref()
                    .ok_or_else(|| {
                        site.not_evaluable_reason
                            .clone()
                            .unwrap_or_else(|| "MaxEnt donor window is unavailable".to_string())
                    })
                    .and_then(|window| model.score_donor(window))
                {
                    Ok(score) => {
                        site.evidence_class = CrypticSplicingEvidenceClass::ModelScored;
                        site.model_status = CrypticSplicingModelStatus::Present;
                        site.maxent_score = Some(score);
                        site.not_evaluable_reason = None;
                    }
                    Err(error) => {
                        site.model_status = CrypticSplicingModelStatus::Failed;
                        site.not_evaluable_reason = Some(error);
                    }
                }
            }
            for site in &mut acceptor_sites {
                match site
                    .model_window
                    .as_deref()
                    .ok_or_else(|| {
                        site.not_evaluable_reason
                            .clone()
                            .unwrap_or_else(|| "MaxEnt acceptor window is unavailable".to_string())
                    })
                    .and_then(|window| model.score_acceptor(window))
                {
                    Ok(score) => {
                        site.evidence_class = CrypticSplicingEvidenceClass::ModelScored;
                        site.model_status = CrypticSplicingModelStatus::Present;
                        site.maxent_score = Some(score);
                        site.not_evaluable_reason = None;
                    }
                    Err(error) => {
                        site.model_status = CrypticSplicingModelStatus::Failed;
                        site.not_evaluable_reason = Some(error);
                    }
                }
            }
        }
        let donor_ids = donor_sites
            .iter()
            .map(|site| site.site_id.clone())
            .collect::<Vec<_>>();
        let acceptor_ids = acceptor_sites
            .iter()
            .map(|site| site.site_id.clone())
            .collect::<Vec<_>>();

        let model_ranking_active = active_model.is_some();
        let mut admissible_pair_count = 0usize;
        let mut evaluated_pair_count = 0usize;
        let mut positional_pair_indices = Vec::with_capacity(request.max_candidate_pairs);
        let mut ranked_pairs = std::collections::BinaryHeap::with_capacity(
            request.max_candidate_pairs.saturating_add(1),
        );
        for (donor_index, donor_position) in donor_positions.iter().enumerate() {
            let min_acceptor = donor_position
                .saturating_add(request.min_pseudo_intron_bp)
                .saturating_sub(2);
            let max_acceptor = donor_position
                .saturating_add(request.max_pseudo_intron_bp)
                .saturating_sub(2);
            let first = acceptor_positions.partition_point(|position| *position < min_acceptor);
            let after_last =
                acceptor_positions.partition_point(|position| *position <= max_acceptor);
            admissible_pair_count = admissible_pair_count.saturating_add(after_last - first);
            if model_ranking_active {
                for acceptor_index in first..after_last {
                    evaluated_pair_count = evaluated_pair_count.saturating_add(1);
                    let ranked = CrypticRankedPair {
                        donor_index,
                        acceptor_index,
                        score: donor_sites[donor_index]
                            .maxent_score
                            .zip(acceptor_sites[acceptor_index].maxent_score)
                            .map(|(donor, acceptor)| donor + acceptor),
                    };
                    if ranked_pairs.len() < request.max_candidate_pairs {
                        ranked_pairs.push(ranked);
                    } else if ranked_pairs
                        .peek()
                        .is_some_and(|current_worst| ranked < *current_worst)
                    {
                        ranked_pairs.pop();
                        ranked_pairs.push(ranked);
                    }
                }
            } else if positional_pair_indices.len() < request.max_candidate_pairs {
                let remaining = request.max_candidate_pairs - positional_pair_indices.len();
                positional_pair_indices.extend(
                    (first..after_last)
                        .take(remaining)
                        .map(|acceptor_index| (donor_index, acceptor_index)),
                );
            }
        }
        let pair_indices = if model_ranking_active {
            let mut selected = ranked_pairs.into_vec();
            selected.sort();
            selected
                .into_iter()
                .map(|row| (row.donor_index, row.acceptor_index))
                .collect::<Vec<_>>()
        } else {
            evaluated_pair_count = positional_pair_indices.len();
            positional_pair_indices
        };

        let cds_context = Self::cryptic_cds_context(dna, request)?;
        let mut candidates = Vec::with_capacity(pair_indices.len());
        for (donor_index, acceptor_index) in pair_indices {
            let donor_0based = donor_positions[donor_index];
            let acceptor_0based = acceptor_positions[acceptor_index];
            let pseudo_intron_end_exclusive = acceptor_0based + 2;
            let oriented_candidate = &oriented_sequence[donor_0based..pseudo_intron_end_exclusive];
            let donor_source_1based = Self::cryptic_local_to_source_1based(request, donor_0based);
            let acceptor_source_1based =
                Self::cryptic_local_to_source_1based(request, acceptor_0based);
            let terminal_source_1based =
                Self::cryptic_local_to_source_1based(request, pseudo_intron_end_exclusive - 1);
            let removed_source_start_1based = donor_source_1based
                .min(acceptor_source_1based)
                .min(terminal_source_1based);
            let removed_source_end_1based = donor_source_1based
                .max(acceptor_source_1based)
                .max(terminal_source_1based);
            let (canonical_pair, paired_motif_signature, motif_class, _annotation) =
                Self::classify_splice_boundary_pair("GT", "AG");
            debug_assert!(canonical_pair);
            let identity = format!(
                "{source_digest}|{}|{}|{}|{}",
                request.strand.as_str(),
                donor_0based + 1,
                acceptor_0based + 1,
                pseudo_intron_end_exclusive - donor_0based
            );
            candidates.push(CrypticSplicingCandidateRow {
                candidate_id: short_sha256_id("cryptic_pair", &identity),
                evidence_class: CrypticSplicingEvidenceClass::StructuralCandidate,
                donor_site_id: donor_ids[donor_index].clone(),
                acceptor_site_id: acceptor_ids[acceptor_index].clone(),
                donor_scanned_position_1based: donor_0based + 1,
                acceptor_scanned_position_1based: acceptor_0based + 1,
                donor_source_position_1based: donor_source_1based,
                acceptor_source_position_1based: acceptor_source_1based,
                removed_source_start_1based,
                removed_source_end_1based,
                pseudo_intron_length_bp: pseudo_intron_end_exclusive - donor_0based,
                paired_motif_signature,
                motif_class: motif_class.to_string(),
                boundary_class: Self::cryptic_boundary_class(
                    request,
                    donor_source_1based,
                    acceptor_source_1based,
                ),
                branchpoint: Self::cryptic_branchpoint_signal(
                    request,
                    donor_0based,
                    oriented_candidate,
                ),
                polypyrimidine_tract: Self::cryptic_polypyrimidine_signal(
                    request,
                    donor_0based,
                    oriented_candidate,
                ),
                model_status: match (
                    donor_sites[donor_index].model_status,
                    acceptor_sites[acceptor_index].model_status,
                ) {
                    (CrypticSplicingModelStatus::Present, CrypticSplicingModelStatus::Present) => {
                        CrypticSplicingModelStatus::Present
                    }
                    (CrypticSplicingModelStatus::Absent, _)
                    | (_, CrypticSplicingModelStatus::Absent) => CrypticSplicingModelStatus::Absent,
                    _ => CrypticSplicingModelStatus::Failed,
                },
                donor_maxent_score: donor_sites[donor_index].maxent_score,
                acceptor_maxent_score: acceptor_sites[acceptor_index].maxent_score,
                prioritization_heuristic: donor_sites[donor_index]
                    .maxent_score
                    .zip(acceptor_sites[acceptor_index].maxent_score)
                    .map(|(donor, acceptor)| donor + acceptor),
                cds_consequence: Self::cryptic_cds_consequence(
                    cds_context.as_ref(),
                    removed_source_start_1based,
                    removed_source_end_1based,
                ),
                warnings: vec![],
            });
        }
        if model_ranking_active {
            for candidate in &mut candidates {
                if candidate.model_status == CrypticSplicingModelStatus::Present {
                    candidate.evidence_class = CrypticSplicingEvidenceClass::ModelScored;
                }
            }
        }
        let truncated = admissible_pair_count > candidates.len();
        let mut warnings = vec![
            "Structural GT-AG candidates are sequence hypotheses, not evidence that splicing occurs. Validate prioritized candidates with an appropriate splice-site model and experimental RNA evidence."
                .to_string(),
            model_provenance.note.clone(),
        ];
        let truncation_rule = truncated.then(|| {
            if model_ranking_active {
                format!(
                    "Reported the top {} pairs after evaluating all {} admissible pairs by donor_maxent_score + acceptor_maxent_score; unscored pairs rank after scored pairs and ties use donor then acceptor position.",
                    candidates.len(), admissible_pair_count
                )
            } else {
                format!(
                    "Reported the first {} admissible pairs in deterministic donor-position then acceptor-position order.",
                    candidates.len()
                )
            }
        });
        if let Some(rule) = truncation_rule.as_ref() {
            warnings.push(format!(
                "Candidate output was truncated: {admissible_pair_count} pairs met the length bounds. {rule}"
            ));
        }
        if request.cds_feature_id.is_none() {
            warnings.push(
                "No CDS feature was selected; coding and protein consequences are unavailable."
                    .to_string(),
            );
        }
        let reported_pair_count = candidates.len();

        Ok(CrypticSplicingScreenView {
            schema: CRYPTIC_SPLICING_SCREEN_SCHEMA.to_string(),
            op_id: None,
            run_id: None,
            request: request.clone(),
            request_sha256,
            effective_input_sha256,
            source_digest,
            source_length_bp: dna.len(),
            scanned_sequence_length_bp: normalized.len(),
            coordinate_space: "scanned_span".to_string(),
            coordinate_convention: "Request and source coordinates are 1-based inclusive. Scanned positions are 1-based in requested strand orientation; reverse-strand positions increase from source end toward source start."
                .to_string(),
            genomic_provenance: Self::cryptic_genomic_provenance(dna, request),
            model: model_provenance,
            donor_sites,
            acceptor_sites,
            candidates,
            budget: CrypticSplicingBudgetSummary {
                donor_site_count: donor_positions.len(),
                acceptor_site_count: acceptor_positions.len(),
                admissible_pair_count,
                evaluated_pair_count,
                reported_pair_count,
                ranking_complete: model_ranking_active || !truncated,
                truncated,
                truncation_rule,
            },
            warnings,
        })
    }

    pub fn inspect_cryptic_splicing_evidence_overlay(
        &self,
        request: &CrypticSplicingEvidenceOverlayRequest,
    ) -> Result<CrypticSplicingEvidenceOverlayReport, EngineError> {
        let screen = self.inspect_cryptic_splicing_screen(&request.screen_request)?;
        let dna = self
            .state
            .sequences
            .get(&request.screen_request.seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", request.screen_request.seq_id),
                cause_chain: vec![],
            })?;
        let source_sequence_sha256 = sha256_prefixed_bytes(dna.forward_bytes());
        let request_sha256 = serde_json::to_vec(request)
            .map(|bytes| sha256_prefixed_bytes(&bytes))
            .map_err(|error| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Could not fingerprint cryptic-splicing evidence request: {error}"
                ),
                cause_chain: vec![],
            })?;
        let mut report_ids = request
            .rna_read_report_ids
            .iter()
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>();
        report_ids.sort();
        report_ids.dedup();
        let mut bindings = Vec::with_capacity(report_ids.len());
        let mut bound_reports = Vec::new();
        let mut warnings = vec![
            "RNA junction support is an observed-evidence layer over structural candidates; lack of retained support is not evidence that a candidate is biologically absent."
                .to_string(),
        ];
        for report_id in report_ids {
            let report = self.get_rna_read_report(&report_id)?;
            let report_sha256 = serde_json::to_vec(&report)
                .map(|bytes| sha256_prefixed_bytes(&bytes))
                .map_err(|error| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Could not fingerprint RNA-read report '{}': {error}",
                        report.report_id
                    ),
                    cause_chain: vec![],
                })?;
            let binding_error = if report.seq_id != request.screen_request.seq_id {
                Some(format!(
                    "RNA report seq_id '{}' does not match screen seq_id '{}'",
                    report.seq_id, request.screen_request.seq_id
                ))
            } else if report.source_sequence_sha256.as_deref()
                != Some(source_sequence_sha256.as_str())
            {
                Some(match report.source_sequence_sha256.as_deref() {
                    Some(actual) => format!(
                        "RNA report source digest '{}' does not match current sequence digest '{}'",
                        actual, source_sequence_sha256
                    ),
                    None => "Legacy RNA report has no source-sequence digest; re-run it before coordinate-bound evidence projection."
                        .to_string(),
                })
            } else if report.coordinate_space
                != gentle_protocol::RNA_READ_LOADED_SEQUENCE_COORDINATE_SPACE
                || report.coordinate_convention
                    != gentle_protocol::RNA_READ_LOADED_SEQUENCE_COORDINATE_CONVENTION
            {
                Some(format!(
                    "RNA report coordinate binding '{}/{}' is incompatible with the expected '{}/{}'",
                    report.coordinate_space,
                    report.coordinate_convention,
                    gentle_protocol::RNA_READ_LOADED_SEQUENCE_COORDINATE_SPACE,
                    gentle_protocol::RNA_READ_LOADED_SEQUENCE_COORDINATE_CONVENTION
                ))
            } else {
                None
            };
            let status = if binding_error.is_none() {
                CrypticSplicingEvidenceBindingStatus::Bound
            } else {
                CrypticSplicingEvidenceBindingStatus::NotAssessable
            };
            let note = binding_error.unwrap_or_else(|| {
                "RNA report is bound to the current sequence digest and coordinate convention."
                    .to_string()
            });
            if status == CrypticSplicingEvidenceBindingStatus::Bound {
                bound_reports.push(report.clone());
            } else {
                warnings.push(format!("RNA report '{}': {note}", report.report_id));
            }
            bindings.push(CrypticSplicingRnaReportBinding {
                report_id: report.report_id,
                report_sha256,
                seq_id: report.seq_id,
                source_sequence_sha256: report.source_sequence_sha256,
                coordinate_space: report.coordinate_space,
                coordinate_convention: report.coordinate_convention,
                status,
                note,
            });
        }

        let candidates = screen
            .candidates
            .iter()
            .map(|candidate| {
                let mut matches = Vec::new();
                for report in &bound_reports {
                    for junction in &report.junction_support_frequencies {
                        if junction.support_read_count == 0 {
                            continue;
                        }
                        let reported_boundary_low_1based =
                            junction.donor_1based.min(junction.acceptor_1based);
                        let reported_boundary_high_1based =
                            junction.donor_1based.max(junction.acceptor_1based);
                        let candidate_boundary_low_1based = candidate
                            .donor_source_position_1based
                            .min(candidate.acceptor_source_position_1based);
                        let candidate_boundary_high_1based = candidate
                            .donor_source_position_1based
                            .max(candidate.acceptor_source_position_1based);
                        let boundary_low_distance_bp = reported_boundary_low_1based
                            .abs_diff(candidate_boundary_low_1based);
                        let boundary_high_distance_bp = reported_boundary_high_1based
                            .abs_diff(candidate_boundary_high_1based);
                        if boundary_low_distance_bp <= request.nearby_tolerance_bp
                            && boundary_high_distance_bp <= request.nearby_tolerance_bp
                        {
                            matches.push(CrypticSplicingRnaJunctionMatch {
                                report_id: report.report_id.clone(),
                                reported_boundary_low_1based,
                                reported_boundary_high_1based,
                                boundary_low_distance_bp,
                                boundary_high_distance_bp,
                                support_read_count: junction.support_read_count,
                                support_fraction: junction.support_fraction,
                                exact: boundary_low_distance_bp == 0
                                    && boundary_high_distance_bp == 0,
                            });
                        }
                    }
                }
                matches.sort_by(|left, right| {
                    (!left.exact)
                        .cmp(&(!right.exact))
                        .then(
                            (left.boundary_low_distance_bp + left.boundary_high_distance_bp).cmp(
                                &(right.boundary_low_distance_bp
                                    + right.boundary_high_distance_bp),
                            ),
                        )
                        .then(right.support_read_count.cmp(&left.support_read_count))
                        .then(left.report_id.cmp(&right.report_id))
                });
                let status = if matches.iter().any(|row| row.exact) {
                    CrypticSplicingRnaEvidenceStatus::ExactObservedJunction
                } else if !matches.is_empty() {
                    CrypticSplicingRnaEvidenceStatus::NearbyObservedJunction
                } else if bound_reports.is_empty() {
                    CrypticSplicingRnaEvidenceStatus::NotAssessable
                } else {
                    CrypticSplicingRnaEvidenceStatus::NoRetainedJunctionSupport
                };
                let interpretation = match status {
                    CrypticSplicingRnaEvidenceStatus::ExactObservedJunction => {
                        "At least one bound RNA-read report retained junction support at these exact donor/acceptor coordinates."
                    }
                    CrypticSplicingRnaEvidenceStatus::NearbyObservedJunction => {
                        "A bound RNA-read report retained a nearby junction within the configured coordinate tolerance; inspect alignment uncertainty before interpreting it."
                    }
                    CrypticSplicingRnaEvidenceStatus::NoRetainedJunctionSupport => {
                        "No positive retained junction row matched this candidate in the bound reports; this is not evidence of biological absence."
                    }
                    CrypticSplicingRnaEvidenceStatus::NotAssessable => {
                        "No RNA-read report could be bound to the current sequence and coordinate convention."
                    }
                }
                .to_string();
                CrypticSplicingCandidateEvidenceOverlayRow {
                    candidate_id: candidate.candidate_id.clone(),
                    status,
                    matches,
                    interpretation,
                }
            })
            .collect::<Vec<_>>();
        let report_id = short_sha256_id(
            "cryptic_overlay",
            &format!(
                "{}|{}|{}",
                request_sha256, screen.effective_input_sha256, source_sequence_sha256
            ),
        );
        Ok(CrypticSplicingEvidenceOverlayReport {
            schema: CRYPTIC_SPLICING_EVIDENCE_OVERLAY_SCHEMA.to_string(),
            report_id,
            op_id: None,
            run_id: None,
            request_sha256,
            screen_request_sha256: screen.request_sha256.clone(),
            screen_effective_input_sha256: screen.effective_input_sha256.clone(),
            screen_source_digest: screen.source_digest.clone(),
            seq_id: request.screen_request.seq_id.clone(),
            source_sequence_sha256,
            coordinate_space: gentle_protocol::RNA_READ_LOADED_SEQUENCE_COORDINATE_SPACE
                .to_string(),
            coordinate_convention: gentle_protocol::RNA_READ_LOADED_SEQUENCE_COORDINATE_CONVENTION
                .to_string(),
            screen,
            rna_report_bindings: bindings,
            candidates,
            warnings,
        })
    }

    pub fn inspect_cryptic_splicing_protein_projection(
        &self,
        request: &CrypticSplicingProteinProjectionRequest,
    ) -> Result<CrypticSplicingProteinProjectionReport, EngineError> {
        let screen = self.inspect_cryptic_splicing_screen(&request.screen_request)?;
        let projection = self.get_uniprot_genome_projection(&request.uniprot_projection_id)?;
        let dna = self
            .state
            .sequences
            .get(&request.screen_request.seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", request.screen_request.seq_id),
                cause_chain: vec![],
            })?;
        let source_sequence_sha256 = sha256_prefixed_bytes(dna.forward_bytes());
        let request_sha256 = serde_json::to_vec(request)
            .map(|bytes| sha256_prefixed_bytes(&bytes))
            .map_err(|error| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Could not fingerprint cryptic-splicing protein projection request: {error}"
                ),
                cause_chain: vec![],
            })?;
        let uniprot_projection_sha256 = serde_json::to_vec(&projection)
            .map(|bytes| sha256_prefixed_bytes(&bytes))
            .map_err(|error| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Could not fingerprint UniProt projection '{}': {error}",
                    projection.projection_id
                ),
                cause_chain: vec![],
            })?;
        let binding_error = if projection.seq_id != request.screen_request.seq_id {
            Some(format!(
                "UniProt projection seq_id '{}' does not match screen seq_id '{}'",
                projection.seq_id, request.screen_request.seq_id
            ))
        } else if projection.source_sequence_sha256.as_deref()
            != Some(source_sequence_sha256.as_str())
        {
            Some(match projection.source_sequence_sha256.as_deref() {
                Some(actual) => format!(
                    "UniProt projection source digest '{}' does not match current sequence digest '{}'",
                    actual, source_sequence_sha256
                ),
                None => "Legacy UniProt projection has no source-sequence digest; re-project it before candidate-to-protein evidence mapping."
                    .to_string(),
            })
        } else {
            None
        };
        let binding_status = if binding_error.is_none() {
            CrypticSplicingEvidenceBindingStatus::Bound
        } else {
            CrypticSplicingEvidenceBindingStatus::NotAssessable
        };
        let binding_note = binding_error.unwrap_or_else(|| {
            "UniProt projection is bound to the current loaded-sequence digest.".to_string()
        });
        let candidates = screen
            .candidates
            .iter()
            .map(|candidate| {
                let mut overlapping_features = Vec::new();
                if binding_status == CrypticSplicingEvidenceBindingStatus::Bound {
                    for transcript in &projection.transcript_projections {
                        for feature in &transcript.feature_projections {
                            for segment in &feature.genomic_segments {
                                let segment_start = segment
                                    .genomic_start_1based
                                    .min(segment.genomic_end_1based);
                                let segment_end = segment
                                    .genomic_start_1based
                                    .max(segment.genomic_end_1based);
                                let overlap_start = segment_start
                                    .max(candidate.removed_source_start_1based);
                                let overlap_end =
                                    segment_end.min(candidate.removed_source_end_1based);
                                if overlap_start <= overlap_end {
                                    overlapping_features.push(
                                        CrypticSplicingProjectedProteinFeatureRow {
                                            transcript_id: transcript.transcript_id.clone(),
                                            feature_key: feature.feature_key.clone(),
                                            feature_note: feature.feature_note.clone(),
                                            aa_start: feature.aa_start,
                                            aa_end: feature.aa_end,
                                            genomic_overlap_start_1based: overlap_start,
                                            genomic_overlap_end_1based: overlap_end,
                                        },
                                    );
                                }
                            }
                        }
                    }
                }
                overlapping_features.sort_by(|left, right| {
                    left.transcript_id
                        .cmp(&right.transcript_id)
                        .then(left.feature_key.cmp(&right.feature_key))
                        .then(left.aa_start.cmp(&right.aa_start))
                        .then(
                            left.genomic_overlap_start_1based
                                .cmp(&right.genomic_overlap_start_1based),
                        )
                });
                overlapping_features.dedup();
                let status = if binding_status
                    == CrypticSplicingEvidenceBindingStatus::NotAssessable
                {
                    CrypticSplicingProteinProjectionStatus::NotAssessable
                } else if overlapping_features.is_empty() {
                    CrypticSplicingProteinProjectionStatus::NoProjectedFeatureOverlap
                } else {
                    CrypticSplicingProteinProjectionStatus::OverlapsProjectedFeature
                };
                let interpretation = match status {
                    CrypticSplicingProteinProjectionStatus::OverlapsProjectedFeature => {
                        "The candidate removal overlaps one or more feature intervals in the bound UniProt projection; this is a coordinate consequence, not evidence that splicing occurs."
                    }
                    CrypticSplicingProteinProjectionStatus::NoProjectedFeatureOverlap => {
                        "No annotated UniProt feature interval overlaps the candidate removal in this projection; coding or protein consequences may still exist outside annotated features."
                    }
                    CrypticSplicingProteinProjectionStatus::NotAssessable => {
                        "The UniProt projection could not be bound to the current sequence digest."
                    }
                }
                .to_string();
                CrypticSplicingCandidateProteinProjectionRow {
                    candidate_id: candidate.candidate_id.clone(),
                    status,
                    overlapping_features,
                    interpretation,
                }
            })
            .collect();
        let warnings = vec![
            "Projected UniProt feature overlap is a separate annotation-evidence layer and does not establish cryptic splice usage or a complete altered open reading frame."
                .to_string(),
        ];
        let report_id = short_sha256_id(
            "cryptic_protein",
            &format!(
                "{}|{}|{}|{}",
                request_sha256,
                screen.effective_input_sha256,
                uniprot_projection_sha256,
                source_sequence_sha256
            ),
        );
        Ok(CrypticSplicingProteinProjectionReport {
            schema: CRYPTIC_SPLICING_PROTEIN_PROJECTION_SCHEMA.to_string(),
            report_id,
            op_id: None,
            run_id: None,
            request_sha256,
            screen_request_sha256: screen.request_sha256.clone(),
            screen_effective_input_sha256: screen.effective_input_sha256.clone(),
            screen_source_digest: screen.source_digest.clone(),
            seq_id: request.screen_request.seq_id.clone(),
            source_sequence_sha256,
            screen,
            uniprot_projection_id: projection.projection_id,
            uniprot_projection_sha256,
            uniprot_entry_id: projection.entry_id,
            binding_status,
            binding_note,
            candidates,
            warnings,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gb_io::seq::{Feature, Location};

    static MODEL_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

    struct RestoreDefaultModel;

    impl Drop for RestoreDefaultModel {
        fn drop(&mut self) {
            crate::maxent_splicing::reload_from_path(None);
        }
    }

    fn base4_sequence(mut index: usize, len: usize) -> String {
        let mut bases = vec![b'A'; len];
        for position in (0..len).rev() {
            bases[position] = match index % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            index /= 4;
        }
        String::from_utf8(bases).expect("base4 DNA")
    }

    fn synthetic_maxent_snapshot(scale: f64) -> crate::maxent_splicing::MaxEntSpliceModelSnapshot {
        let mut snapshot = crate::maxent_splicing::MaxEntSpliceModelSnapshot {
            schema: crate::maxent_splicing::MAXENT_SPLICE_MODEL_SCHEMA.to_string(),
            source: "synthetic cryptic-splicing test model".to_string(),
            redistribution_status: "synthetic_no_restriction".to_string(),
            donor_sequences: (0..16_384).map(|index| base4_sequence(index, 7)).collect(),
            donor_scores: vec![scale; 16_384],
            acceptor_tables: [16_384, 16_384, 16_384, 16_384, 16_384, 64, 256, 64, 256]
                .into_iter()
                .map(|len| vec![scale; len])
                .collect(),
            ..crate::maxent_splicing::MaxEntSpliceModelSnapshot::default()
        };
        snapshot.finalize_fingerprint().expect("fingerprint");
        snapshot
    }

    fn engine_with_sequence(seq_id: &str, sequence: &str, with_cds: bool) -> GentleEngine {
        let mut dna = DNAsequence::from_sequence(sequence).expect("synthetic DNA");
        if with_cds {
            dna.features_mut().push(Feature {
                kind: "CDS".into(),
                location: Location::simple_range(0, sequence.len() as i64),
                qualifiers: vec![("transl_table".into(), Some("1".to_string()))],
            });
        }
        let mut state = ProjectState::default();
        state.sequences.insert(seq_id.to_string(), dna);
        GentleEngine::from_state(state)
    }

    fn structural_sequence() -> String {
        "AAAGTCCCCCTACTAACCCCCCCCCCCCCCCCCCCAGAAA".to_string()
    }

    fn request(seq_id: &str, len: usize) -> CrypticSplicingScreenRequest {
        CrypticSplicingScreenRequest {
            seq_id: seq_id.to_string(),
            start_1based: 1,
            end_1based: len,
            min_pseudo_intron_bp: 20,
            max_pseudo_intron_bp: 200,
            max_candidate_pairs: 100,
            ..CrypticSplicingScreenRequest::default()
        }
    }

    #[test]
    fn cryptic_splicing_structural_screen_keeps_missing_model_explicit() {
        let sequence = structural_sequence();
        let engine = engine_with_sequence("cassette", &sequence, true);
        let mut screen_request = request("cassette", sequence.len());
        screen_request.cds_feature_id = Some(0);
        let view = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("structural screen");

        assert_eq!(view.schema, CRYPTIC_SPLICING_SCREEN_SCHEMA);
        assert_eq!(view.model.status, CrypticSplicingModelStatus::Absent);
        assert!(!view.candidates.is_empty());
        let candidate = &view.candidates[0];
        assert_eq!(candidate.model_status, CrypticSplicingModelStatus::Absent);
        assert_eq!(candidate.donor_maxent_score, None);
        assert_eq!(candidate.acceptor_maxent_score, None);
        assert_eq!(
            candidate.branchpoint.status,
            CrypticSplicingSignalStatus::Detected
        );
        assert_eq!(
            candidate.polypyrimidine_tract.status,
            CrypticSplicingSignalStatus::Detected
        );
        let consequence = candidate
            .cds_consequence
            .as_ref()
            .expect("selected CDS consequence");
        assert!(consequence.removed_coding_bp > 0);
        assert!(consequence.first_altered_aa_position.is_some());
    }

    #[test]
    fn cryptic_splicing_cds_consequence_distinguishes_in_frame_and_frameshift_removals() {
        let coding_sequence = b"ATGAAACCCGGGTAA".to_vec();
        let context = CrypticCdsContext {
            feature_id: 7,
            strand_matches: true,
            source_positions_0based: (0..coding_sequence.len()).collect(),
            coding_sequence,
            translation_table: 1,
        };

        let in_frame = GentleEngine::cryptic_cds_consequence(Some(&context), 4, 6)
            .expect("in-frame consequence");
        assert_eq!(in_frame.status, "coding_in_frame_removal");
        assert_eq!(in_frame.removed_coding_bp, 3);
        assert_eq!(in_frame.frame_delta_mod_3, Some(0));

        let frameshift = GentleEngine::cryptic_cds_consequence(Some(&context), 4, 7)
            .expect("frameshift consequence");
        assert_eq!(frameshift.status, "coding_frameshift");
        assert_eq!(frameshift.removed_coding_bp, 4);
        assert_eq!(frameshift.frame_delta_mod_3, Some(1));
    }

    #[test]
    fn cryptic_splicing_boundary_class_reports_donor_and_acceptor_compartments() {
        let without_insert = request("cassette", 40);
        assert_eq!(
            GentleEngine::cryptic_boundary_class(&without_insert, 12, 24),
            "insert_span_not_specified"
        );

        let mut with_insert = without_insert;
        with_insert.insert_span = Some(CrypticSplicingSpan {
            start_1based: 10,
            end_1based: 30,
        });
        assert_eq!(
            GentleEngine::cryptic_boundary_class(&with_insert, 12, 24),
            "insert_to_insert"
        );
        assert_eq!(
            GentleEngine::cryptic_boundary_class(&with_insert, 2, 8),
            "vector_to_vector"
        );
        assert_eq!(
            GentleEngine::cryptic_boundary_class(&with_insert, 12, 35),
            "insert_to_vector"
        );
        assert_eq!(
            GentleEngine::cryptic_boundary_class(&with_insert, 2, 12),
            "vector_to_insert"
        );
    }

    #[test]
    fn cryptic_splicing_reverse_screen_uses_requested_orientation() {
        let forward = structural_sequence();
        let reverse = String::from_utf8(GentleEngine::reverse_complement_bytes(forward.as_bytes()))
            .expect("reverse complement text");
        let engine = engine_with_sequence("reverse_cassette", &reverse, false);
        let mut screen_request = request("reverse_cassette", reverse.len());
        screen_request.strand = CrypticSplicingStrand::Reverse;
        let view = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("reverse structural screen");

        assert!(!view.candidates.is_empty());
        assert_eq!(
            view.candidates[0].donor_scanned_position_1based,
            forward.find("GT").expect("donor") + 1
        );
        assert!(
            view.candidates[0].donor_source_position_1based
                > view.candidates[0].acceptor_source_position_1based
        );
    }

    #[test]
    fn cryptic_splicing_candidate_budget_is_deterministic_and_reported() {
        let sequence = "GTAAAAAAAAAAAAAAAAAAAAAG".repeat(12);
        let engine = engine_with_sequence("repeat", &sequence, false);
        let mut screen_request = request("repeat", sequence.len());
        screen_request.max_pseudo_intron_bp = sequence.len();
        screen_request.max_candidate_pairs = 3;
        let first = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("first screen");
        let second = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("second screen");

        assert!(first.budget.truncated);
        assert_eq!(first.candidates.len(), 3);
        assert_eq!(first, second);
        assert!(first.budget.admissible_pair_count > first.candidates.len());
        assert_eq!(first.budget.evaluated_pair_count, 3);
        assert_eq!(first.budget.reported_pair_count, 3);
        assert!(!first.budget.ranking_complete);
        assert!(first.budget.truncation_rule.is_some());
    }

    #[test]
    fn cryptic_splicing_rejects_implicit_circular_wraparound() {
        let engine = engine_with_sequence("circle", "AAAGTCCCCAGAAA", false);
        let error = engine
            .inspect_cryptic_splicing_screen(&CrypticSplicingScreenRequest {
                seq_id: "circle".to_string(),
                start_1based: 12,
                end_1based: 4,
                min_pseudo_intron_bp: 4,
                max_pseudo_intron_bp: 20,
                max_candidate_pairs: 10,
                ..CrypticSplicingScreenRequest::default()
            })
            .expect_err("wraparound must be explicit in a later schema");
        assert!(error.message.contains("does not support a circular span"));
    }

    #[test]
    fn cryptic_splicing_short_candidate_is_not_misreported_as_signal_absence() {
        let sequence = "AAAGTAAAAAGAAA";
        let engine = engine_with_sequence("short", sequence, false);
        let view = engine
            .inspect_cryptic_splicing_screen(&CrypticSplicingScreenRequest {
                seq_id: "short".to_string(),
                start_1based: 1,
                end_1based: sequence.len(),
                min_pseudo_intron_bp: 4,
                max_pseudo_intron_bp: 20,
                max_candidate_pairs: 10,
                ..CrypticSplicingScreenRequest::default()
            })
            .expect("short screen");
        let candidate = view.candidates.first().expect("short candidate");
        assert_eq!(
            candidate.branchpoint.status,
            CrypticSplicingSignalStatus::NotEvaluable
        );
        assert!(candidate.branchpoint.not_evaluable_reason.is_some());
    }

    #[test]
    fn cryptic_splicing_model_scores_and_effective_identity_follow_the_active_snapshot() {
        let _lock = MODEL_TEST_LOCK.lock().expect("model test lock");
        let _restore = RestoreDefaultModel;
        let temp = tempfile::tempdir().expect("tempdir");
        let first_path = temp.path().join("model-1.json");
        let second_path = temp.path().join("model-2.json");
        std::fs::write(
            &first_path,
            serde_json::to_vec(&synthetic_maxent_snapshot(1.0)).expect("first JSON"),
        )
        .expect("write first model");
        std::fs::write(
            &second_path,
            serde_json::to_vec(&synthetic_maxent_snapshot(2.0)).expect("second JSON"),
        )
        .expect("write second model");

        let sequence = structural_sequence();
        let engine = engine_with_sequence("cassette", &sequence, false);
        let mut screen_request = request("cassette", sequence.len());
        screen_request.model_policy = CrypticSplicingModelPolicy::UseActiveMaxent;

        crate::maxent_splicing::reload_from_path(Some(first_path.to_string_lossy().as_ref()));
        let first = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("first model screen");
        assert_eq!(first.model.status, CrypticSplicingModelStatus::Present);
        assert!(first.candidates.iter().any(|row| {
            row.donor_maxent_score.is_some() && row.acceptor_maxent_score.is_some()
        }));

        crate::maxent_splicing::reload_from_path(Some(second_path.to_string_lossy().as_ref()));
        let second = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("replacement model screen");
        assert_ne!(first.model.resource_sha256, second.model.resource_sha256);
        assert_ne!(first.effective_input_sha256, second.effective_input_sha256);

        screen_request.expected_model_fingerprint_sha256 = first.model.resource_sha256.clone();
        let mismatch = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("fingerprint mismatch retains structure");
        assert_eq!(mismatch.model.status, CrypticSplicingModelStatus::Failed);
        assert!(!mismatch.candidates.is_empty());
        assert!(
            mismatch
                .candidates
                .iter()
                .all(|row| row.donor_maxent_score.is_none())
        );
    }

    #[test]
    fn cryptic_splicing_model_budget_selects_global_top_pair() {
        let _lock = MODEL_TEST_LOCK.lock().expect("model test lock");
        let _restore = RestoreDefaultModel;
        let temp = tempfile::tempdir().expect("tempdir");
        let model_path = temp.path().join("ranked-model.json");
        let mut snapshot = synthetic_maxent_snapshot(1.0);
        let preferred_index = snapshot
            .donor_sequences
            .iter()
            .position(|sequence| sequence == "CCCCCCC")
            .expect("preferred donor table row");
        snapshot.donor_scores[preferred_index] = 100.0;
        snapshot
            .finalize_fingerprint()
            .expect("updated fingerprint");
        std::fs::write(
            &model_path,
            serde_json::to_vec(&snapshot).expect("model JSON"),
        )
        .expect("write ranked model");
        crate::maxent_splicing::reload_from_path(Some(model_path.to_string_lossy().as_ref()));

        let sequence = format!(
            "AAAGTAAAA{}CCCGTCCCC{}AGAAA",
            "A".repeat(20),
            "A".repeat(25)
        );
        let engine = engine_with_sequence("ranked", &sequence, false);
        let mut screen_request = request("ranked", sequence.len());
        screen_request.min_pseudo_intron_bp = 4;
        screen_request.max_pseudo_intron_bp = sequence.len();
        screen_request.max_candidate_pairs = 1;
        screen_request.model_policy = CrypticSplicingModelPolicy::UseActiveMaxent;
        let report = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("ranked screen");

        assert!(report.budget.truncated);
        assert!(report.budget.ranking_complete);
        assert_eq!(
            report.budget.evaluated_pair_count,
            report.budget.admissible_pair_count
        );
        assert_eq!(report.budget.reported_pair_count, 1);
        assert_eq!(
            report.candidates[0].donor_scanned_position_1based,
            sequence.find("CCCGTCCCC").expect("preferred donor window") + 4
        );
        assert!(
            report
                .budget
                .truncation_rule
                .as_deref()
                .is_some_and(|rule| rule.contains("top 1 pairs"))
        );
    }

    #[test]
    fn cryptic_splicing_rna_overlay_binds_exact_junctions_and_refuses_legacy_coordinates() {
        let sequence = structural_sequence();
        let mut engine = engine_with_sequence("cassette", &sequence, false);
        let screen_request = request("cassette", sequence.len());
        let screen = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("screen");
        let candidate = screen.candidates.first().expect("candidate");
        let source_sequence_sha256 = sha256_prefixed_bytes(sequence.as_bytes());
        let bound_report = RnaReadInterpretationReport {
            schema: RNA_READ_REPORT_SCHEMA.to_string(),
            report_id: "bound_rna".to_string(),
            seq_id: "cassette".to_string(),
            source_sequence_sha256: Some(source_sequence_sha256),
            coordinate_space: gentle_protocol::RNA_READ_LOADED_SEQUENCE_COORDINATE_SPACE
                .to_string(),
            coordinate_convention: gentle_protocol::RNA_READ_LOADED_SEQUENCE_COORDINATE_CONVENTION
                .to_string(),
            junction_support_frequencies: vec![RnaReadJunctionSupportFrequency {
                donor_1based: candidate
                    .donor_source_position_1based
                    .min(candidate.acceptor_source_position_1based),
                acceptor_1based: candidate
                    .donor_source_position_1based
                    .max(candidate.acceptor_source_position_1based),
                support_read_count: 9,
                support_fraction: 0.75,
            }],
            ..RnaReadInterpretationReport::default()
        };
        engine
            .upsert_rna_read_report(bound_report.clone())
            .expect("store bound RNA report");
        let overlay = engine
            .inspect_cryptic_splicing_evidence_overlay(&CrypticSplicingEvidenceOverlayRequest {
                screen_request: screen_request.clone(),
                rna_read_report_ids: vec!["bound_rna".to_string()],
                nearby_tolerance_bp: 2,
            })
            .expect("bound overlay");
        assert_eq!(
            overlay.rna_report_bindings[0].status,
            CrypticSplicingEvidenceBindingStatus::Bound
        );
        assert_eq!(
            overlay.candidates[0].status,
            CrypticSplicingRnaEvidenceStatus::ExactObservedJunction
        );

        let mut legacy_value = serde_json::to_value(bound_report).expect("RNA report JSON");
        let legacy_object = legacy_value.as_object_mut().expect("RNA report object");
        legacy_object.insert(
            "report_id".to_string(),
            serde_json::Value::String("legacy_rna".to_string()),
        );
        legacy_object.remove("source_sequence_sha256");
        legacy_object.remove("coordinate_space");
        legacy_object.remove("coordinate_convention");
        let legacy_report: RnaReadInterpretationReport =
            serde_json::from_value(legacy_value).expect("legacy RNA report still deserializes");
        assert_eq!(legacy_report.source_sequence_sha256, None);
        assert!(legacy_report.coordinate_space.is_empty());
        engine
            .upsert_rna_read_report(legacy_report)
            .expect("store legacy RNA report");
        let legacy = engine
            .inspect_cryptic_splicing_evidence_overlay(&CrypticSplicingEvidenceOverlayRequest {
                screen_request,
                rna_read_report_ids: vec!["legacy_rna".to_string()],
                nearby_tolerance_bp: 2,
            })
            .expect("legacy overlay remains inspectable");
        assert_eq!(
            legacy.rna_report_bindings[0].status,
            CrypticSplicingEvidenceBindingStatus::NotAssessable
        );
        assert!(
            legacy
                .candidates
                .iter()
                .all(|row| { row.status == CrypticSplicingRnaEvidenceStatus::NotAssessable })
        );
    }

    #[test]
    fn cryptic_splicing_rna_overlay_matches_reverse_candidates_by_normalized_geometry() {
        let forward = structural_sequence();
        let reverse = String::from_utf8(GentleEngine::reverse_complement_bytes(forward.as_bytes()))
            .expect("reverse complement text");
        let mut engine = engine_with_sequence("reverse_cassette", &reverse, false);
        let mut screen_request = request("reverse_cassette", reverse.len());
        screen_request.strand = CrypticSplicingStrand::Reverse;
        let candidate = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("reverse screen")
            .candidates
            .into_iter()
            .next()
            .expect("reverse candidate");
        assert!(candidate.donor_source_position_1based > candidate.acceptor_source_position_1based);
        engine
            .upsert_rna_read_report(RnaReadInterpretationReport {
                schema: RNA_READ_REPORT_SCHEMA.to_string(),
                report_id: "reverse_rna".to_string(),
                seq_id: "reverse_cassette".to_string(),
                source_sequence_sha256: Some(sha256_prefixed_bytes(reverse.as_bytes())),
                coordinate_space: gentle_protocol::RNA_READ_LOADED_SEQUENCE_COORDINATE_SPACE
                    .to_string(),
                coordinate_convention:
                    gentle_protocol::RNA_READ_LOADED_SEQUENCE_COORDINATE_CONVENTION.to_string(),
                junction_support_frequencies: vec![RnaReadJunctionSupportFrequency {
                    donor_1based: candidate.acceptor_source_position_1based,
                    acceptor_1based: candidate.donor_source_position_1based,
                    support_read_count: 4,
                    support_fraction: 0.5,
                }],
                ..RnaReadInterpretationReport::default()
            })
            .expect("store reverse RNA report");

        let overlay = engine
            .inspect_cryptic_splicing_evidence_overlay(&CrypticSplicingEvidenceOverlayRequest {
                screen_request,
                rna_read_report_ids: vec!["reverse_rna".to_string()],
                nearby_tolerance_bp: 0,
            })
            .expect("reverse overlay");
        assert_eq!(
            overlay.candidates[0].status,
            CrypticSplicingRnaEvidenceStatus::ExactObservedJunction
        );
        assert!(overlay.candidates[0].matches[0].exact);
    }

    #[test]
    fn cryptic_splicing_projects_only_digest_bound_uniprot_features() {
        let sequence = structural_sequence();
        let mut engine = engine_with_sequence("cassette", &sequence, false);
        let screen_request = request("cassette", sequence.len());
        let candidate = engine
            .inspect_cryptic_splicing_screen(&screen_request)
            .expect("screen")
            .candidates
            .into_iter()
            .next()
            .expect("candidate");
        let source_sequence_sha256 = sha256_prefixed_bytes(sequence.as_bytes());
        engine
            .upsert_uniprot_projection(crate::uniprot::UniprotGenomeProjection {
                schema: UNIPROT_GENOME_PROJECTION_SCHEMA.to_string(),
                projection_id: "PTEST@cassette".to_string(),
                entry_id: "PTEST".to_string(),
                seq_id: "cassette".to_string(),
                source_sequence_sha256: Some(source_sequence_sha256),
                transcript_projections: vec![crate::uniprot::UniprotTranscriptProjection {
                    transcript_id: "TX1".to_string(),
                    feature_projections: vec![crate::uniprot::UniprotFeatureProjection {
                        feature_key: "DNA_BIND".to_string(),
                        feature_note: Some("synthetic DNA-binding region".to_string()),
                        aa_start: 10,
                        aa_end: 20,
                        genomic_segments: vec![crate::uniprot::UniprotAaGenomicSegment {
                            aa_start: 10,
                            aa_end: 20,
                            genomic_start_1based: candidate.removed_source_start_1based,
                            genomic_end_1based: candidate.removed_source_end_1based,
                            strand: "+".to_string(),
                        }],
                    }],
                    ..crate::uniprot::UniprotTranscriptProjection::default()
                }],
                ..crate::uniprot::UniprotGenomeProjection::default()
            })
            .expect("store projection");
        let report = engine
            .inspect_cryptic_splicing_protein_projection(&CrypticSplicingProteinProjectionRequest {
                screen_request: screen_request.clone(),
                uniprot_projection_id: "PTEST@cassette".to_string(),
            })
            .expect("protein projection");
        assert_eq!(
            report.binding_status,
            CrypticSplicingEvidenceBindingStatus::Bound
        );
        assert_eq!(
            report.candidates[0].status,
            CrypticSplicingProteinProjectionStatus::OverlapsProjectedFeature
        );

        engine
            .upsert_uniprot_projection(crate::uniprot::UniprotGenomeProjection {
                schema: UNIPROT_GENOME_PROJECTION_SCHEMA.to_string(),
                projection_id: "PLEGACY@cassette".to_string(),
                entry_id: "PLEGACY".to_string(),
                seq_id: "cassette".to_string(),
                ..crate::uniprot::UniprotGenomeProjection::default()
            })
            .expect("store legacy projection");
        let legacy = engine
            .inspect_cryptic_splicing_protein_projection(&CrypticSplicingProteinProjectionRequest {
                screen_request,
                uniprot_projection_id: "PLEGACY@cassette".to_string(),
            })
            .expect("legacy projection remains inspectable");
        assert_eq!(
            legacy.binding_status,
            CrypticSplicingEvidenceBindingStatus::NotAssessable
        );
        assert!(
            legacy
                .candidates
                .iter()
                .all(|row| { row.status == CrypticSplicingProteinProjectionStatus::NotAssessable })
        );
    }
}
