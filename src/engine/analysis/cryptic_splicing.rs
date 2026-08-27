//! Deterministic structural screening for potential cryptic splice removal.
//!
//! This screen reports sequence architecture and heuristic support. It does
//! not claim that a candidate is used by the spliceosome. Optional statistical
//! splice-site models can populate the reserved score fields in a later slice.

use super::*;
use crate::{AMINO_ACIDS, amino_acids::STOP_CODON};
use gentle_protocol::{
    CRYPTIC_SPLICING_SCREEN_SCHEMA, CrypticSplicingBranchpointSignal, CrypticSplicingBudgetSummary,
    CrypticSplicingCandidateRow, CrypticSplicingCdsConsequence, CrypticSplicingEvidenceClass,
    CrypticSplicingGenomicProvenanceRow, CrypticSplicingModelProvenance,
    CrypticSplicingModelStatus, CrypticSplicingPolypyrimidineSignal, CrypticSplicingScreenRequest,
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
        let donor_sites = donor_positions
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
        let acceptor_sites = acceptor_positions
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
        let donor_ids = donor_sites
            .iter()
            .map(|site| site.site_id.clone())
            .collect::<Vec<_>>();
        let acceptor_ids = acceptor_sites
            .iter()
            .map(|site| site.site_id.clone())
            .collect::<Vec<_>>();

        let mut admissible_pair_count = 0usize;
        let mut pair_indices = Vec::with_capacity(request.max_candidate_pairs);
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
            if pair_indices.len() < request.max_candidate_pairs {
                let remaining = request.max_candidate_pairs - pair_indices.len();
                pair_indices.extend(
                    (first..after_last)
                        .take(remaining)
                        .map(|acceptor_index| (donor_index, acceptor_index)),
                );
            }
        }

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
                model_status: CrypticSplicingModelStatus::Absent,
                donor_maxent_score: None,
                acceptor_maxent_score: None,
                prioritization_heuristic: None,
                cds_consequence: Self::cryptic_cds_consequence(
                    cds_context.as_ref(),
                    removed_source_start_1based,
                    removed_source_end_1based,
                ),
                warnings: vec![],
            });
        }
        let truncated = admissible_pair_count > candidates.len();
        let mut warnings = vec![
            "Structural GT-AG candidates are sequence hypotheses, not evidence that splicing occurs. Validate prioritized candidates with an appropriate splice-site model and experimental RNA evidence."
                .to_string(),
            "No optional MaxEnt model resource was applied; donor and acceptor model scores are unavailable rather than zero."
                .to_string(),
        ];
        let truncation_rule = truncated.then(|| {
            format!(
                "Reported the first {} admissible pairs in deterministic donor-position then acceptor-position order.",
                candidates.len()
            )
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

        Ok(CrypticSplicingScreenView {
            schema: CRYPTIC_SPLICING_SCREEN_SCHEMA.to_string(),
            request: request.clone(),
            request_sha256,
            source_digest,
            source_length_bp: dna.len(),
            scanned_sequence_length_bp: normalized.len(),
            coordinate_space: "scanned_span".to_string(),
            coordinate_convention: "Request and source coordinates are 1-based inclusive. Scanned positions are 1-based in requested strand orientation; reverse-strand positions increase from source end toward source start."
                .to_string(),
            genomic_provenance: Self::cryptic_genomic_provenance(dna, request),
            model: CrypticSplicingModelProvenance {
                status: CrypticSplicingModelStatus::Absent,
                model_kind: None,
                resource_path: None,
                resource_sha256: None,
                note: "No optional statistical splice-site model was applied; v1 reports the dependency-free structural screen."
                    .to_string(),
            },
            donor_sites,
            acceptor_sites,
            candidates,
            budget: CrypticSplicingBudgetSummary {
                donor_site_count: donor_positions.len(),
                acceptor_site_count: acceptor_positions.len(),
                admissible_pair_count,
                evaluated_pair_count: admissible_pair_count.min(request.max_candidate_pairs),
                truncated,
                truncation_rule,
            },
            warnings,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gb_io::seq::{Feature, Location};

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
}
