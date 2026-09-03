//! Transcript-aware, read-only promoter-reporter architecture comparison.
//!
//! The analysis composes promoter TSS classes, the existing fail-closed CDS
//! boundary audit, optional motif scans, and optional CUT&RUN support. It does
//! not create reporter constructs or infer promoter usage from occupancy.

use super::*;
use std::path::Path;

impl GentleEngine {
    fn promoter_architecture_transcript_id_stem(raw: &str) -> &str {
        raw.trim()
            .split_once('.')
            .map_or(raw.trim(), |(stem, _)| stem)
    }

    fn promoter_architecture_transcript_ids_match(left: &str, right: &str) -> bool {
        left.trim().eq_ignore_ascii_case(right.trim())
            || Self::promoter_architecture_transcript_id_stem(left)
                .eq_ignore_ascii_case(Self::promoter_architecture_transcript_id_stem(right))
    }

    fn promoter_architecture_genomic_strand(
        local_strand: &str,
        anchor: Option<&SequenceGenomeAnchorSummary>,
    ) -> String {
        match (local_strand.trim(), anchor.and_then(|value| value.strand)) {
            ("+", Some('-')) => "-",
            ("-", Some('-')) => "+",
            ("+", _) => "+",
            ("-", _) => "-",
            _ => "?",
        }
        .to_string()
    }

    fn promoter_architecture_explicit_tss_class_id(
        class_id: &str,
        local_strand: &str,
        genomic_strand: &str,
    ) -> String {
        let local = if local_strand == "-" { "minus" } else { "plus" };
        let genomic = if genomic_strand == "-" {
            "minus"
        } else {
            "plus"
        };
        let legacy_marker = format!("_{local}_tss_");
        let explicit_marker = format!("_local_{local}_genomic_{genomic}_tss_");
        if class_id.contains(&legacy_marker) {
            class_id.replacen(&legacy_marker, &explicit_marker, 1)
        } else {
            format!("{class_id}_local_{local}_genomic_{genomic}")
        }
    }

    fn promoter_architecture_local_point_to_genomic(
        anchor: Option<&SequenceGenomeAnchorSummary>,
        local_0based: usize,
    ) -> Option<usize> {
        let anchor = anchor?;
        if anchor.strand == Some('-') {
            anchor.end_1based.checked_sub(local_0based)
        } else {
            anchor.start_1based.checked_add(local_0based)
        }
    }

    fn promoter_architecture_local_range_to_genomic(
        anchor: Option<&SequenceGenomeAnchorSummary>,
        start_0based: usize,
        end_0based_exclusive: usize,
    ) -> Option<(usize, usize)> {
        let anchor = anchor?;
        if end_0based_exclusive <= start_0based {
            return None;
        }
        if anchor.strand == Some('-') {
            Some((
                anchor
                    .end_1based
                    .checked_sub(end_0based_exclusive.saturating_sub(1))?,
                anchor.end_1based.checked_sub(start_0based)?,
            ))
        } else {
            Some((
                anchor.start_1based.checked_add(start_0based)?,
                anchor
                    .start_1based
                    .checked_add(end_0based_exclusive.saturating_sub(1))?,
            ))
        }
    }

    fn promoter_architecture_feature_ranges(
        source: &DNAsequence,
        kind: &str,
        transcript_id: &str,
        is_reverse: bool,
    ) -> Vec<(usize, usize)> {
        let mut ranges = vec![];
        for feature in source.features().iter().filter(|feature| {
            feature.kind.to_string().eq_ignore_ascii_case(kind)
                && feature_is_reverse(feature) == is_reverse
                && feature.qualifier_values("transcript_id").any(|value| {
                    Self::promoter_architecture_transcript_ids_match(value, transcript_id)
                })
        }) {
            collect_location_ranges_usize(&feature.location, &mut ranges);
        }
        ranges.retain(|(start, end)| end > start);
        ranges.sort_unstable();
        ranges.dedup();
        if is_reverse {
            ranges.sort_unstable_by(|left, right| {
                right.1.cmp(&left.1).then_with(|| right.0.cmp(&left.0))
            });
        }
        ranges
    }

    fn promoter_architecture_support_metadata(
        transcript: &gb_io::seq::Feature,
    ) -> BTreeMap<String, Vec<String>> {
        let mut metadata = BTreeMap::new();
        for key in [
            "transcript_support_level",
            "tag",
            "biotype",
            "transcript_biotype",
            "source",
            "logic_name",
            "havana_transcript",
            "ccds",
        ] {
            let mut values = transcript
                .qualifier_values(key)
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string)
                .collect::<Vec<_>>();
            values.sort_by_key(|value| value.to_ascii_lowercase());
            values.dedup_by(|left, right| left.eq_ignore_ascii_case(right));
            if !values.is_empty() {
                metadata.insert(key.to_string(), values);
            }
        }
        metadata
    }

    fn promoter_architecture_transcript_audit(
        &self,
        source: &DNAsequence,
        record: &PromoterWindowRecord,
        class_id: &str,
        anchor: Option<&SequenceGenomeAnchorSummary>,
    ) -> Result<PromoterReporterTranscriptArchitectureAudit, EngineError> {
        let transcript_feature_id = record.transcript_feature_id.ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Transcript '{}' has no stable feature id for reporter architecture comparison",
                record.transcript_id
            ),
            cause_chain: vec![],
        })?;
        let transcript = source
            .features()
            .get(transcript_feature_id)
            .filter(|feature| Self::is_mrna_feature(feature))
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Transcript '{}' feature {} is not an mRNA annotation",
                    record.transcript_id,
                    transcript_feature_id.saturating_add(1)
                ),
                cause_chain: vec![],
            })?;
        let dummy_candidate = PromoterReporterFragmentCandidate {
            candidate_id: format!("{}_boundary_audit", record.transcript_id),
            transcript_id: record.transcript_id.clone(),
            strand: record.strand.clone(),
            variant_start_0based: record.tss_local_0based,
            variant_end_0based_exclusive: record.tss_local_0based.saturating_add(1),
            start_0based: 0,
            end_0based_exclusive: source.len(),
            ..PromoterReporterFragmentCandidate::default()
        };
        let (_, _, boundary) = Self::promoter_reporter_panel_extended_geometry(
            source,
            &dummy_candidate,
            &PromoterReporterPanelExtendedBoundaryPolicy {
                kind: PromoterReporterPanelExtendedBoundaryKind::CanonicalCdsStartCodon,
                transcript_id: record.transcript_id.clone(),
            },
        )?;
        let is_reverse = record.strand == "-";
        let mut exon_ranges = Self::promoter_architecture_feature_ranges(
            source,
            "exon",
            &record.transcript_id,
            is_reverse,
        );
        if exon_ranges.is_empty() {
            collect_location_ranges_usize(&transcript.location, &mut exon_ranges);
            exon_ranges.retain(|(start, end)| end > start);
            exon_ranges.sort_unstable();
            exon_ranges.dedup();
            if is_reverse {
                exon_ranges.sort_unstable_by(|left, right| right.1.cmp(&left.1));
            }
        }
        let cds_ranges = Self::promoter_architecture_feature_ranges(
            source,
            "CDS",
            &record.transcript_id,
            is_reverse,
        );
        let cds_start_codon_genomic_ranges_1based = boundary
            .cds_start_codon_source_ranges_0based
            .iter()
            .filter_map(|(start, end)| {
                Self::promoter_architecture_local_range_to_genomic(anchor, *start, *end)
            })
            .collect::<Vec<_>>();
        let five_prime_utr_introns = boundary
            .five_prime_utr_intron_ranges_0based
            .iter()
            .enumerate()
            .map(
                |(index, (start, end))| PromoterReporterArchitectureIntronAudit {
                    start_0based: *start,
                    end_0based_exclusive: *end,
                    length_bp: end.saturating_sub(*start),
                    transcript_order: index + 1,
                },
            )
            .collect::<Vec<_>>();
        let tss_to_atg_through_span_bp = if is_reverse {
            record
                .tss_local_0based
                .saturating_sub(boundary.cds_start_codon_source_start_0based)
                .saturating_add(1)
        } else {
            boundary
                .cds_start_codon_source_end_0based_exclusive
                .saturating_sub(record.tss_local_0based)
        };
        let local_strand = record.strand.clone();
        let genomic_strand = Self::promoter_architecture_genomic_strand(&local_strand, anchor);
        Ok(PromoterReporterTranscriptArchitectureAudit {
            transcript_id: record.transcript_id.clone(),
            transcript_label: record.transcript_label.clone(),
            transcript_feature_id,
            gene_label: record.gene_label.clone(),
            strand: local_strand.clone(),
            local_strand,
            genomic_strand,
            tss_class_id: class_id.to_string(),
            tss_local_0based: record.tss_local_0based,
            tss_genomic_1based: Self::promoter_architecture_local_point_to_genomic(
                anchor,
                record.tss_local_0based,
            ),
            cds_start_codon_source_ranges_0based: boundary.cds_start_codon_source_ranges_0based,
            cds_start_codon_genomic_ranges_1based,
            cds_start_codon_5prime_to_3prime: boundary.cds_start_codon_5prime_to_3prime,
            tss_to_atg_through_span_bp,
            transcript_exon_ranges_0based: exon_ranges,
            cds_ranges_0based: cds_ranges,
            five_prime_utr_exon_ranges_0based: boundary.five_prime_utr_exon_ranges_0based,
            five_prime_utr_spliced_length_bp: boundary.five_prime_utr_spliced_length_bp,
            five_prime_utr_introns,
            five_prime_utr_total_intron_bp: boundary.five_prime_utr_total_intron_bp,
            five_prime_utr_max_intron_bp: boundary.five_prime_utr_max_intron_bp,
            kozak_context_5prime_to_3prime: boundary.kozak_context_5prime_to_3prime,
            kozak_context_class: boundary.kozak_context_class,
            support_metadata: Self::promoter_architecture_support_metadata(transcript),
            warnings: boundary.warnings,
        })
    }

    fn promoter_architecture_segment(
        architecture_id: &str,
        index: usize,
        kind: &str,
        material: &str,
        range: Option<(usize, usize)>,
        length_bp: usize,
        description: &str,
    ) -> PromoterReporterArchitectureSegment {
        PromoterReporterArchitectureSegment {
            segment_id: format!("{architecture_id}_segment_{index}"),
            kind: kind.to_string(),
            material: material.to_string(),
            source_start_0based: range.map(|(start, _)| start),
            source_end_0based_exclusive: range.map(|(_, end)| end),
            length_bp,
            transcript_orientation: "5prime_to_3prime".to_string(),
            description: description.to_string(),
        }
    }

    fn promoter_architecture_source_bounds(
        source_len: usize,
        audit: &PromoterReporterTranscriptArchitectureAudit,
        upstream_bp: usize,
    ) -> (usize, usize) {
        if audit.strand == "-" {
            (
                audit.tss_local_0based.saturating_add(1),
                audit
                    .tss_local_0based
                    .saturating_add(upstream_bp)
                    .saturating_add(1)
                    .min(source_len),
            )
        } else {
            (
                audit.tss_local_0based.saturating_sub(upstream_bp),
                audit.tss_local_0based,
            )
        }
    }

    fn promoter_architecture_stable_id(
        class_id: &str,
        transcript_id: &str,
        kind: PromoterReporterArchitectureKind,
        policy_suffix: &str,
    ) -> String {
        let identity = format!(
            "{class_id}|{transcript_id}|{}|{policy_suffix}",
            kind.as_str()
        );
        let digest = crate::digest_utils::sha256_hex_str(&identity);
        format!(
            "arch_{}_{}_{}_{}",
            Self::normalize_id_token(class_id),
            Self::normalize_id_token(Self::promoter_architecture_transcript_id_stem(
                transcript_id
            )),
            kind.as_str(),
            &digest[..12]
        )
    }

    fn promoter_architecture_cds_prefix_segments(
        architecture_id: &str,
        audit: &PromoterReporterTranscriptArchitectureAudit,
        retained_bp: usize,
        first_segment_index: usize,
    ) -> Result<Vec<PromoterReporterArchitectureSegment>, EngineError> {
        let mut remaining = retained_bp;
        let mut segments = vec![];
        for (range_index, (start, end)) in audit.cds_ranges_0based.iter().copied().enumerate() {
            let take = remaining.min(end.saturating_sub(start));
            if take == 0 {
                continue;
            }
            let selected = if audit.strand == "-" {
                (end - take, end)
            } else {
                (start, start + take)
            };
            segments.push(Self::promoter_architecture_segment(
                architecture_id,
                first_segment_index + range_index,
                "endogenous_cds_prefix",
                "spliced_genomic_exon",
                Some(selected),
                take,
                "Transcript-ordered endogenous CDS prefix retained before luciferase",
            ));
            remaining -= take;
            if remaining == 0 {
                break;
            }
        }
        if remaining != 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Fusion request for transcript '{}' retains {} CDS bp, but only {} annotated coding bp were available",
                    audit.transcript_id,
                    retained_bp,
                    retained_bp.saturating_sub(remaining)
                ),
                cause_chain: vec![],
            });
        }
        Ok(segments)
    }

    #[allow(clippy::too_many_arguments)]
    fn promoter_architecture_row(
        source_len: usize,
        class: &PromoterReporterTssClass,
        audit: &PromoterReporterTranscriptArchitectureAudit,
        kind: PromoterReporterArchitectureKind,
        promoter_upstream_bp: usize,
        tss_proximal_downstream_bp: usize,
        fusion: Option<&PromoterReporterFusionRequest>,
    ) -> Result<PromoterReporterArchitectureRow, EngineError> {
        let policy_suffix = fusion.map_or("default".to_string(), |request| {
            let junction = request
                .junction_sequence_5prime_to_3prime
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_ascii_uppercase)
                .unwrap_or_else(|| "none".to_string());
            format!(
                "{}|{}|{:?}|{}",
                request.mode.as_str(),
                request.retained_endogenous_cds_bp,
                request.vector_luciferase_atg_removed,
                junction
            )
        });
        let architecture_id = Self::promoter_architecture_stable_id(
            &class.class_id,
            &audit.transcript_id,
            kind,
            &policy_suffix,
        );
        let promoter_range =
            Self::promoter_architecture_source_bounds(source_len, audit, promoter_upstream_bp);
        let mut segments = vec![];
        let junction: PromoterReporterArchitectureJunctionAudit;
        let mut major_confounds = vec![];
        let mut cloning_feasibility_notes = vec![];
        let (modeled_insert_length_bp, testable_question) = match kind {
            PromoterReporterArchitectureKind::TssProximalTranscriptional => {
                let range = if audit.strand == "-" {
                    (
                        audit
                            .tss_local_0based
                            .saturating_sub(tss_proximal_downstream_bp),
                        promoter_range.1,
                    )
                } else {
                    (
                        promoter_range.0,
                        audit
                            .tss_local_0based
                            .saturating_add(tss_proximal_downstream_bp)
                            .saturating_add(1)
                            .min(source_len),
                    )
                };
                segments.push(Self::promoter_architecture_segment(
                    &architecture_id,
                    1,
                    "tss_proximal_genomic_promoter",
                    "genomic_dna",
                    Some(range),
                    range.1.saturating_sub(range.0),
                    "Genomic TSS-proximal fragment; the vector luciferase CDS and ATG remain intact",
                ));
                junction = PromoterReporterArchitectureJunctionAudit {
                    endogenous_atg_disposition:
                        PromoterReporterAtgDisposition::EndogenousExcludedVectorLuciferaseRetained,
                    reporter_atg_source: "vector_luciferase_atg".to_string(),
                    vector_luciferase_atg_removed: Some(false),
                    frame_status: "not_applicable_transcriptional_reporter".to_string(),
                    audit_status: "explicit".to_string(),
                    ..PromoterReporterArchitectureJunctionAudit::default()
                };
                major_confounds.push(
                    "The chosen downstream TSS distance may include leader sequence or nearby regulatory elements"
                        .to_string(),
                );
                (
                    range.1.saturating_sub(range.0),
                    "Compare promoter-proximal transcriptional output while leaving luciferase translation under the vector's own ATG".to_string(),
                )
            }
            PromoterReporterArchitectureKind::Spliced5utrLucAtgReplacement => {
                segments.push(Self::promoter_architecture_segment(
                    &architecture_id,
                    1,
                    "genomic_promoter_upstream_of_tss",
                    "genomic_dna",
                    Some(promoter_range),
                    promoter_range.1.saturating_sub(promoter_range.0),
                    "Genomic promoter ending at the annotated TSS boundary",
                ));
                for (index, range) in audit
                    .five_prime_utr_exon_ranges_0based
                    .iter()
                    .copied()
                    .enumerate()
                {
                    segments.push(Self::promoter_architecture_segment(
                        &architecture_id,
                        index + 2,
                        "spliced_5utr_exon",
                        "spliced_cdna",
                        Some(range),
                        range.1.saturating_sub(range.0),
                        "Transcript-oriented mature 5' UTR exon segment",
                    ));
                }
                segments.push(Self::promoter_architecture_segment(
                    &architecture_id,
                    segments.len() + 1,
                    "luciferase_start_codon",
                    "reporter_sequence",
                    None,
                    3,
                    "Luciferase ATG placed at the endogenous translation-start position",
                ));
                junction = PromoterReporterArchitectureJunctionAudit {
                    endogenous_atg_disposition:
                        PromoterReporterAtgDisposition::EndogenousReplacedByLuciferase,
                    reporter_atg_source: "luciferase_cds_atg".to_string(),
                    vector_luciferase_atg_removed: Some(false),
                    frame_status: "not_applicable_atg_replacement".to_string(),
                    audit_status: "explicit".to_string(),
                    ..PromoterReporterArchitectureJunctionAudit::default()
                };
                major_confounds.extend([
                    "The construct is a synthetic genomic-promoter/cDNA-leader hybrid".to_string(),
                    "Leader-dependent translation and RNA stability contribute to reporter output"
                        .to_string(),
                ]);
                let length = promoter_range.1.saturating_sub(promoter_range.0)
                    + audit.five_prime_utr_spliced_length_bp
                    + 3;
                (
                    length,
                    "Test effects of the mature transcript leader while excluding its genomic introns".to_string(),
                )
            }
            PromoterReporterArchitectureKind::Genomic5utrLucAtgReplacement => {
                let codon_outer_start = audit
                    .cds_start_codon_source_ranges_0based
                    .iter()
                    .map(|(start, _)| *start)
                    .min()
                    .unwrap_or_default();
                let codon_outer_end = audit
                    .cds_start_codon_source_ranges_0based
                    .iter()
                    .map(|(_, end)| *end)
                    .max()
                    .unwrap_or_default();
                let range = if audit.strand == "-" {
                    (codon_outer_end, promoter_range.1)
                } else {
                    (promoter_range.0, codon_outer_start)
                };
                if range.1 <= range.0 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Transcript '{}' genomic 5' UTR reporter range is empty",
                            audit.transcript_id
                        ),
                        cause_chain: vec![],
                    });
                }
                segments.push(Self::promoter_architecture_segment(
                    &architecture_id,
                    1,
                    "genomic_promoter_and_5utr",
                    "genomic_dna",
                    Some(range),
                    range.1 - range.0,
                    "Contiguous genomic promoter and 5' leader ending immediately before the endogenous ATG",
                ));
                segments.push(Self::promoter_architecture_segment(
                    &architecture_id,
                    2,
                    "luciferase_start_codon",
                    "reporter_sequence",
                    None,
                    3,
                    "Luciferase ATG replaces the excluded endogenous ATG",
                ));
                junction = PromoterReporterArchitectureJunctionAudit {
                    endogenous_atg_disposition:
                        PromoterReporterAtgDisposition::EndogenousReplacedByLuciferase,
                    reporter_atg_source: "luciferase_cds_atg".to_string(),
                    vector_luciferase_atg_removed: Some(false),
                    frame_status: "not_applicable_atg_replacement".to_string(),
                    audit_status: "explicit".to_string(),
                    ..PromoterReporterArchitectureJunctionAudit::default()
                };
                major_confounds.extend([
                    "Cis sequence, spacing, splicing, and RNA processing are changed together"
                        .to_string(),
                    "Reporter output cannot be interpreted as pure promoter activity".to_string(),
                ]);
                if range.1 - range.0 > 5_000 {
                    cloning_feasibility_notes.push(
                        "Long genomic insert exceeds 5 kb and belongs in a higher-complexity second tier"
                            .to_string(),
                    );
                }
                (
                    range.1 - range.0 + 3,
                    "Test the combined contribution of genomic leader sequence, introns, spacing, and RNA processing".to_string(),
                )
            }
            PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion => {
                let fusion = fusion.ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Transcript '{}' requires an explicit fusion request for endogenous_atg_retained_fusion",
                        audit.transcript_id
                    ),
                    cause_chain: vec![],
                })?;
                if fusion.retained_endogenous_cds_bp < 3 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Fusion request for transcript '{}' must retain at least the three-base endogenous ATG",
                            audit.transcript_id
                        ),
                        cause_chain: vec![],
                    });
                }
                segments.push(Self::promoter_architecture_segment(
                    &architecture_id,
                    1,
                    "genomic_promoter_upstream_of_tss",
                    "genomic_dna",
                    Some(promoter_range),
                    promoter_range.1.saturating_sub(promoter_range.0),
                    "Genomic promoter ending at the annotated TSS boundary",
                ));
                for (index, range) in audit
                    .five_prime_utr_exon_ranges_0based
                    .iter()
                    .copied()
                    .enumerate()
                {
                    segments.push(Self::promoter_architecture_segment(
                        &architecture_id,
                        index + 2,
                        "spliced_5utr_exon",
                        "spliced_cdna",
                        Some(range),
                        range.1.saturating_sub(range.0),
                        "Transcript-oriented mature 5' UTR exon segment",
                    ));
                }
                let first_cds_segment = segments.len() + 1;
                segments.extend(Self::promoter_architecture_cds_prefix_segments(
                    &architecture_id,
                    audit,
                    fusion.retained_endogenous_cds_bp,
                    first_cds_segment,
                )?);
                let mut audit_warnings = vec![];
                let (disposition, frame_status, audit_status) = match fusion.mode {
                    PromoterReporterFusionMode::InFrameTranslationFusion => {
                        if fusion.retained_endogenous_cds_bp % 3 != 0 {
                            audit_warnings.push(format!(
                                "Retained endogenous CDS length {} is not divisible by three",
                                fusion.retained_endogenous_cds_bp
                            ));
                        }
                        if fusion.vector_luciferase_atg_removed != Some(true) {
                            audit_warnings.push(
                                "An in-frame fusion requires explicit removal of the vector luciferase ATG"
                                    .to_string(),
                            );
                        }
                        let pass = audit_warnings.is_empty()
                            && fusion.junction_sequence_5prime_to_3prime.is_some();
                        (
                            PromoterReporterAtgDisposition::EndogenousRetainedVectorLuciferaseRemoved,
                            if fusion.retained_endogenous_cds_bp % 3 == 0 {
                                "in_frame_by_retained_length"
                            } else {
                                "out_of_frame_by_retained_length"
                            },
                            if pass { "pass" } else { "incomplete_or_failed" },
                        )
                    }
                    PromoterReporterFusionMode::ExtraUpstreamAtg => {
                        audit_warnings.push(
                            "This is not an audited translational fusion: an endogenous upstream ATG remains before the intact vector luciferase ATG"
                                .to_string(),
                        );
                        (
                            PromoterReporterAtgDisposition::EndogenousRetainedExtraUpstreamOfVectorLuciferase,
                            "not_a_single_orf_fusion",
                            "warning",
                        )
                    }
                };
                junction = PromoterReporterArchitectureJunctionAudit {
                    endogenous_atg_disposition: disposition,
                    reporter_atg_source: "endogenous_atg_then_luciferase_cds".to_string(),
                    vector_luciferase_atg_removed: fusion.vector_luciferase_atg_removed,
                    frame_status: frame_status.to_string(),
                    junction_sequence_5prime_to_3prime: fusion
                        .junction_sequence_5prime_to_3prime
                        .clone(),
                    audit_status: audit_status.to_string(),
                    warnings: audit_warnings,
                };
                major_confounds.extend([
                    "The retained endogenous coding peptide can alter luciferase stability or localization".to_string(),
                    "Translation-initiation and fusion-frame effects are inseparable from transcriptional effects".to_string(),
                ]);
                let length = promoter_range.1.saturating_sub(promoter_range.0)
                    + audit.five_prime_utr_spliced_length_bp
                    + fusion.retained_endogenous_cds_bp;
                (
                    length,
                    "Explicitly test a caller-defined endogenous translation-start/fusion architecture".to_string(),
                )
            }
        };
        Ok(PromoterReporterArchitectureRow {
            architecture_id,
            architecture_kind: kind,
            tss_class_id: class.class_id.clone(),
            representative_transcript_id: class.representative_transcript_id.clone(),
            tss_class_transcript_ids: class.transcript_ids.clone(),
            transcript_id: audit.transcript_id.clone(),
            strand: audit.strand.clone(),
            local_strand: if audit.local_strand.is_empty() {
                audit.strand.clone()
            } else {
                audit.local_strand.clone()
            },
            genomic_strand: audit.genomic_strand.clone(),
            tss_local_0based: audit.tss_local_0based,
            tss_genomic_1based: audit.tss_genomic_1based,
            cds_start_codon_source_ranges_0based: audit
                .cds_start_codon_source_ranges_0based
                .clone(),
            modeled_insert_length_bp,
            modeled_length_scope:
                "regulatory_or_leader_insert_through_start_boundary; excludes luciferase coding body"
                    .to_string(),
            five_prime_utr_spliced_length_bp: audit.five_prime_utr_spliced_length_bp,
            five_prime_utr_total_intron_bp: audit.five_prime_utr_total_intron_bp,
            five_prime_utr_max_intron_bp: audit.five_prime_utr_max_intron_bp,
            five_prime_utr_introns: audit.five_prime_utr_introns.clone(),
            segments,
            junction,
            transcript_warnings: audit.warnings.clone(),
            testable_question,
            major_confounds,
            cloning_feasibility_notes,
            non_claims: vec![
                "This sequence architecture does not establish promoter usage in the assayed cells"
                    .to_string(),
                "This sequence architecture does not predict luciferase output".to_string(),
            ],
        })
    }

    fn promoter_architecture_portable_source(
        raw: Option<String>,
        include_local_paths: bool,
    ) -> Option<String> {
        raw.and_then(|value| {
            let trimmed = value.trim();
            if trimmed.is_empty() {
                return None;
            }
            if include_local_paths
                || trimmed.starts_with("http://")
                || trimmed.starts_with("https://")
                || trimmed.starts_with("ftp://")
            {
                Some(trimmed.to_string())
            } else {
                Path::new(trimmed)
                    .file_name()
                    .and_then(|name| name.to_str())
                    .map(str::to_string)
                    .or_else(|| Some("local_source_redacted".to_string()))
            }
        })
    }

    fn promoter_architecture_portable_source_error(
        message: &str,
        include_local_paths: bool,
    ) -> String {
        if include_local_paths {
            message.to_string()
        } else {
            "Source could not be resolved or prepared; local path details are redacted.".to_string()
        }
    }

    fn promoter_architecture_source_provenance(
        &self,
        request: &PromoterReporterArchitectureComparisonRequest,
        source: &DNAsequence,
    ) -> PromoterReporterArchitectureSourceProvenance {
        let extraction = self.latest_genome_extraction_provenance_for_seq(&request.seq_id);
        PromoterReporterArchitectureSourceProvenance {
            source_sequence_sha256: crate::digest_utils::sha256_hex_str(
                &source.get_forward_string().to_ascii_uppercase(),
            ),
            genome_anchor: self.sequence_genome_anchor_summary(&request.seq_id).ok(),
            extraction_operation: extraction.as_ref().map(|row| row.operation.clone()),
            sequence_source_type: extraction
                .as_ref()
                .and_then(|row| row.sequence_source_type.clone()),
            annotation_source_type: extraction
                .as_ref()
                .and_then(|row| row.annotation_source_type.clone()),
            sequence_source: Self::promoter_architecture_portable_source(
                extraction
                    .as_ref()
                    .and_then(|row| row.sequence_source.clone()),
                request.include_local_source_paths,
            ),
            annotation_source: Self::promoter_architecture_portable_source(
                extraction
                    .as_ref()
                    .and_then(|row| row.annotation_source.clone()),
                request.include_local_source_paths,
            ),
            sequence_sha1: extraction
                .as_ref()
                .and_then(|row| row.sequence_sha1.clone()),
            annotation_sha1: extraction
                .as_ref()
                .and_then(|row| row.annotation_sha1.clone()),
            local_paths_included: request.include_local_source_paths,
        }
    }

    fn promoter_architecture_existing_construct(
        &self,
        input: Option<&PromoterReporterExistingConstructInput>,
    ) -> Result<PromoterReporterExistingConstructStatus, EngineError> {
        let Some(input) = input else {
            return Ok(PromoterReporterExistingConstructStatus {
                status: "existing_construct_unverified".to_string(),
                label: "existing laboratory construct".to_string(),
                detail: "No insert/vector sequence or boundary coordinates were supplied; the existing construct cannot be assigned to a TSS or fusion architecture."
                    .to_string(),
                ..PromoterReporterExistingConstructStatus::default()
            });
        };
        let label = if input.label.trim().is_empty() {
            "existing laboratory construct".to_string()
        } else {
            input.label.trim().to_string()
        };
        let Some(insert_seq_id) = input
            .insert_seq_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        else {
            return Ok(PromoterReporterExistingConstructStatus {
                status: "existing_construct_unverified".to_string(),
                label,
                vector_seq_id: input.vector_seq_id.clone(),
                source_start_0based: input.source_start_0based,
                source_end_0based_exclusive: input.source_end_0based_exclusive,
                detail: "Construct metadata was supplied without an insert seq_id; sequence and boundaries remain unverified."
                    .to_string(),
                ..PromoterReporterExistingConstructStatus::default()
            });
        };
        let insert = self
            .state
            .sequences
            .get(insert_seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Existing construct insert sequence '{}' not found",
                    insert_seq_id
                ),
                cause_chain: vec![],
            })?;
        let text = insert.get_forward_string();
        let (start, end) = match (input.source_start_0based, input.source_end_0based_exclusive) {
            (Some(start), Some(end)) if end > start && end <= text.len() => (start, end),
            (None, None) => (0, text.len()),
            _ => {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Existing construct insert '{}' has invalid or incomplete source boundaries",
                        insert_seq_id
                    ),
                    cause_chain: vec![],
                });
            }
        };
        if let Some(vector_seq_id) = input
            .vector_seq_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            if !self.state.sequences.contains_key(vector_seq_id) {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Existing construct vector sequence '{}' not found",
                        vector_seq_id
                    ),
                    cause_chain: vec![],
                });
            }
        }
        Ok(PromoterReporterExistingConstructStatus {
            status: "sequence_bound_material_identity_unverified".to_string(),
            label,
            insert_seq_id: Some(insert_seq_id.to_string()),
            vector_seq_id: input.vector_seq_id.clone(),
            insert_sequence_sha256: Some(crate::digest_utils::sha256_hex_str(
                &text[start..end].to_ascii_uppercase(),
            )),
            source_start_0based: Some(start),
            source_end_0based_exclusive: Some(end),
            detail: "The project sequence and selected boundaries are hash-bound, but this does not verify the physical laboratory construct or exact vector stock."
                .to_string(),
        })
    }

    fn promoter_architecture_cutrun_role(
        source_id: &str,
        target_factor: Option<&str>,
        sample_label: Option<&str>,
    ) -> String {
        let joined = format!(
            "{} {} {}",
            source_id,
            target_factor.unwrap_or_default(),
            sample_label.unwrap_or_default()
        )
        .to_ascii_lowercase();
        if joined.contains("h3k4me3") {
            "h3k4me3_control".to_string()
        } else if joined.contains("igg") {
            "igg_control".to_string()
        } else if joined.contains("input") || joined.contains("genome_control") {
            "input_control".to_string()
        } else if joined.contains("negative_control") || joined.contains("neg_") {
            "negative_control".to_string()
        } else if joined.contains("tp73") || joined.contains("p73") {
            "p73_antibody".to_string()
        } else if joined.contains("positive_control") || joined.contains("pos_") {
            "unspecified_positive_control".to_string()
        } else {
            "unspecified".to_string()
        }
    }

    fn promoter_architecture_cutrun_evidence(
        &self,
        request: &PromoterReporterArchitectureComparisonRequest,
    ) -> PromoterReporterCutRunEvidence {
        if request.cutrun_dataset_ids.is_empty() && request.cutrun_read_report_ids.is_empty() {
            return PromoterReporterCutRunEvidence {
                evaluation_state: "unevaluated_no_sources_selected".to_string(),
                quantitative_comparison_status: "not_requested".to_string(),
                interpretation: "No CUT&RUN source was selected. Occupancy is unknown and does not affect TSS classification."
                    .to_string(),
                warnings: vec![
                    "Select prepared CUT&RUN datasets or saved ROI read reports to inspect occupancy; no raw reads are acquired automatically."
                        .to_string(),
                ],
                ..PromoterReporterCutRunEvidence::default()
            };
        }
        let mut lanes = vec![];
        let mut warnings = vec![];
        let mut usable_dataset_ids = vec![];
        let mut usable_read_report_ids = vec![];
        for dataset_id in &request.cutrun_dataset_ids {
            match self.show_cutrun_dataset_status(
                dataset_id,
                request.cutrun_catalog_path.as_deref(),
                request.cutrun_cache_dir.as_deref(),
            ) {
                Ok(status) => {
                    let has_evidence = status.peaks.prepared || status.signal.prepared;
                    if has_evidence {
                        usable_dataset_ids.push(status.dataset_id.clone());
                    }
                    let manifest_sha256 = status.manifest.as_ref().and_then(|manifest| {
                        serde_json::to_vec(manifest)
                            .ok()
                            .map(|bytes| crate::digest_utils::sha256_hex_bytes(&bytes))
                    });
                    let role = Self::promoter_architecture_cutrun_role(
                        &status.dataset_id,
                        status.target_factor.as_deref(),
                        status.sample_label.as_deref(),
                    );
                    let mut lane_warnings = status.warnings.clone();
                    if role == "unspecified_positive_control" {
                        lane_warnings.push(
                            "A generic 'pos' or positive-control label is insufficient; declare H3K4me3 explicitly when that is the assay target."
                                .to_string(),
                        );
                    }
                    lanes.push(PromoterReporterCutRunLane {
                        source_kind: "dataset".to_string(),
                        source_id: status.dataset_id,
                        state: if has_evidence {
                            PromoterReporterCutRunLaneState::Evaluated
                        } else if status.prepared {
                            PromoterReporterCutRunLaneState::PreparedNoCompatibleEvidence
                        } else {
                            PromoterReporterCutRunLaneState::NotPrepared
                        },
                        role,
                        target_factor: status.target_factor,
                        sample_label: status.sample_label,
                        condition: status.condition,
                        tissue_or_cell_type: status.tissue_or_cell_type,
                        source_accession: status.source_accession,
                        manifest_sha256,
                        local_source_path: request
                            .include_local_source_paths
                            .then_some(status.manifest_path)
                            .flatten(),
                        support_window_count: 0,
                        warnings: lane_warnings,
                    });
                }
                Err(error) => {
                    let detail = Self::promoter_architecture_portable_source_error(
                        &error.message,
                        request.include_local_source_paths,
                    );
                    warnings.push(format!(
                        "CUT&RUN dataset '{}' was not evaluated: {}",
                        dataset_id, detail
                    ));
                    lanes.push(PromoterReporterCutRunLane {
                        source_kind: "dataset".to_string(),
                        source_id: dataset_id.clone(),
                        state: PromoterReporterCutRunLaneState::Unevaluated,
                        role: "unresolved".to_string(),
                        warnings: vec![detail],
                        ..PromoterReporterCutRunLane::default()
                    });
                }
            }
        }
        for report_id in &request.cutrun_read_report_ids {
            match self.get_cutrun_read_report(report_id) {
                Ok(report) => {
                    usable_read_report_ids.push(report.report_id.clone());
                    lanes.push(PromoterReporterCutRunLane {
                        source_kind: "read_report".to_string(),
                        source_id: report.report_id.clone(),
                        state: PromoterReporterCutRunLaneState::Evaluated,
                        role: Self::promoter_architecture_cutrun_role(
                            &report.report_id,
                            report.target_factor.as_deref(),
                            None,
                        ),
                        target_factor: report.target_factor.clone(),
                        source_accession: report.dataset_id.clone(),
                        local_source_path: request
                            .include_local_source_paths
                            .then_some(report.input_r1_path.clone()),
                        support_window_count: report.support_clusters.len(),
                        warnings: report.warnings.clone(),
                        ..PromoterReporterCutRunLane::default()
                    });
                }
                Err(error) => {
                    let detail = Self::promoter_architecture_portable_source_error(
                        &error.message,
                        request.include_local_source_paths,
                    );
                    warnings.push(format!(
                        "CUT&RUN read report '{}' was not evaluated: {}",
                        report_id, detail
                    ));
                    lanes.push(PromoterReporterCutRunLane {
                        source_kind: "read_report".to_string(),
                        source_id: report_id.clone(),
                        state: PromoterReporterCutRunLaneState::Unevaluated,
                        role: "unresolved".to_string(),
                        warnings: vec![detail],
                        ..PromoterReporterCutRunLane::default()
                    });
                }
            }
        }
        let regulatory_support =
            if usable_dataset_ids.is_empty() && usable_read_report_ids.is_empty() {
                None
            } else {
                match self.inspect_cutrun_regulatory_support_with_catalog(
                    &request.seq_id,
                    &usable_dataset_ids,
                    &usable_read_report_ids,
                    request.cutrun_catalog_path.as_deref(),
                    request.cutrun_cache_dir.as_deref(),
                    None,
                    None,
                    150,
                    &request.cutrun_species_filters,
                ) {
                    Ok(mut report) => {
                        if !request.include_local_source_paths {
                            report.catalog_path = None;
                            report.cache_dir = None;
                        }
                        for lane in &mut lanes {
                            lane.support_window_count = report
                                .support_windows
                                .iter()
                                .filter(|window| {
                                    window
                                        .contributing_dataset_ids
                                        .iter()
                                        .any(|id| id.eq_ignore_ascii_case(&lane.source_id))
                                        || window
                                            .contributing_read_report_ids
                                            .iter()
                                            .any(|id| id.eq_ignore_ascii_case(&lane.source_id))
                                })
                                .count()
                                .max(lane.support_window_count);
                        }
                        warnings.extend(report.warnings.iter().cloned());
                        Some(report)
                    }
                    Err(error) => {
                        let detail = Self::promoter_architecture_portable_source_error(
                            &error.message,
                            request.include_local_source_paths,
                        );
                        warnings.push(format!(
                            "CUT&RUN regulatory support remained unevaluated: {}",
                            detail
                        ));
                        None
                    }
                }
            };
        let evaluated_lanes = lanes
            .iter()
            .filter(|lane| lane.state == PromoterReporterCutRunLaneState::Evaluated)
            .count();
        PromoterReporterCutRunEvidence {
            evaluation_state: if evaluated_lanes > 0 {
                "evaluated_qualitatively"
            } else {
                "not_prepared_or_unevaluated"
            }
            .to_string(),
            quantitative_comparison_status: if evaluated_lanes > 1 {
                "not_comparable_without_explicit_normalization_and_compatible_controls"
            } else {
                "qualitative_only"
            }
            .to_string(),
            lanes,
            regulatory_support,
            interpretation: "CUT&RUN overlap is occupancy/enrichment support only. It is not evidence of TSS usage, promoter activity, direct regulation, biochemical affinity, or luciferase output."
                .to_string(),
            warnings,
        }
    }

    fn promoter_architecture_staged_panel(
        classes: &[PromoterReporterTssClass],
        architectures: &[PromoterReporterArchitectureRow],
        existing: &PromoterReporterExistingConstructStatus,
    ) -> Vec<PromoterReporterPanelRecommendation> {
        let mut rows = vec![PromoterReporterPanelRecommendation {
            recommendation_id: "stage0_existing_construct".to_string(),
            stage: "preflight".to_string(),
            role: "existing_construct_baseline".to_string(),
            rationale: if existing.status == "existing_construct_unverified" {
                "Retain the laboratory construct as a named baseline, but do not assign it to a TSS or fusion model before sequence and boundaries are verified."
            } else {
                "Use the sequence-bound laboratory construct as a baseline after physical identity and vector/insert junctions are independently confirmed."
            }
            .to_string(),
            required_before_materialization:
                "Supply and verify insert sequence, vector identity, and both physical fusion boundaries"
                    .to_string(),
            ..PromoterReporterPanelRecommendation::default()
        }];
        let mut ordered_classes = classes.iter().collect::<Vec<_>>();
        ordered_classes.sort_by(|left, right| {
            right
                .transcript_ids
                .len()
                .cmp(&left.transcript_ids.len())
                .then_with(|| left.min_tss_local_0based.cmp(&right.min_tss_local_0based))
                .then_with(|| left.class_id.cmp(&right.class_id))
        });
        for (class_index, class) in ordered_classes.into_iter().enumerate() {
            let mut class_rows = architectures
                .iter()
                .filter(|row| {
                    row.tss_class_id == class.class_id
                        && Self::promoter_architecture_transcript_ids_match(
                            &row.transcript_id,
                            &class.representative_transcript_id,
                        )
                        && row.architecture_kind
                            != PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion
                })
                .collect::<Vec<_>>();
            class_rows.sort_by_key(|row| match row.architecture_kind {
                PromoterReporterArchitectureKind::TssProximalTranscriptional => 0,
                PromoterReporterArchitectureKind::Spliced5utrLucAtgReplacement => 1,
                PromoterReporterArchitectureKind::Genomic5utrLucAtgReplacement => 2,
                PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion => 3,
            });
            for architecture in class_rows {
                let second_tier = architecture.architecture_kind
                    == PromoterReporterArchitectureKind::Genomic5utrLucAtgReplacement
                    && architecture.modeled_insert_length_bp > 5_000;
                rows.push(PromoterReporterPanelRecommendation {
                    recommendation_id: format!(
                        "stage{}_{}",
                        if class_index == 0 && !second_tier { 1 } else { 2 },
                        architecture.architecture_id
                    ),
                    stage: if class_index == 0 && !second_tier {
                        "first_tier"
                    } else {
                        "second_tier"
                    }
                    .to_string(),
                    role: architecture.architecture_kind.as_str().to_string(),
                    architecture_id: Some(architecture.architecture_id.clone()),
                    rationale: if second_tier {
                        "Retain as a higher-complexity long genomic-leader construct after the shorter main-TSS panel is interpretable."
                    } else if class_index == 0 {
                        "Representative architecture from the largest annotation-supported TSS class."
                    } else {
                        "Alternative TSS-class counterpart; annotation supports the model but does not establish frequent use."
                    }
                    .to_string(),
                    required_before_materialization:
                        "Confirm exact vector, insert boundary, reporter ATG policy, and cloning route"
                            .to_string(),
                });
            }
        }
        rows.extend([
            PromoterReporterPanelRecommendation {
                recommendation_id: "control_promoterless".to_string(),
                stage: "matched_controls".to_string(),
                role: "promoterless_control".to_string(),
                rationale: "Measure vector/background reporter signal without the tested promoter insert."
                    .to_string(),
                required_before_materialization: "Use the same reporter backbone and assay context"
                    .to_string(),
                ..PromoterReporterPanelRecommendation::default()
            },
            PromoterReporterPanelRecommendation {
                recommendation_id: "control_constitutive".to_string(),
                stage: "matched_controls".to_string(),
                role: "constitutive_control".to_string(),
                rationale: "Provide an assay-competence and normalization reference without treating it as a SERPINE1 promoter model."
                    .to_string(),
                required_before_materialization:
                    "Choose and document one sequence-verified constitutive promoter control"
                        .to_string(),
                ..PromoterReporterPanelRecommendation::default()
            },
        ]);
        rows
    }

    pub(crate) fn compare_promoter_reporter_architectures(
        &self,
        mut request: PromoterReporterArchitectureComparisonRequest,
    ) -> Result<PromoterReporterArchitectureComparisonReport, EngineError> {
        if request.schema.trim().is_empty() {
            request.schema = PROMOTER_REPORTER_ARCHITECTURE_REQUEST_SCHEMA.to_string();
        }
        if request.schema != PROMOTER_REPORTER_ARCHITECTURE_REQUEST_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Promoter-reporter architecture request schema '{}' is unsupported; expected '{}'",
                    request.schema, PROMOTER_REPORTER_ARCHITECTURE_REQUEST_SCHEMA
                ),
                cause_chain: vec![],
            });
        }
        request.seq_id = request.seq_id.trim().to_string();
        if request.seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Promoter-reporter architecture comparison requires a seq_id".to_string(),
                cause_chain: vec![],
            });
        }
        let source = self
            .state
            .sequences
            .get(&request.seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", request.seq_id),
                cause_chain: vec![],
            })?;
        let mut requested_architectures = request.architectures.clone();
        if requested_architectures.is_empty() {
            requested_architectures =
                PromoterReporterArchitectureComparisonRequest::default().architectures;
        }
        let mut seen_architectures = BTreeSet::new();
        requested_architectures.retain(|kind| seen_architectures.insert(kind.as_str().to_string()));
        if requested_architectures
            .contains(&PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion)
            && request.fusion_requests.is_empty()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "endogenous_atg_retained_fusion requires at least one explicit fusion request"
                        .to_string(),
                cause_chain: vec![],
            });
        }
        if !request.fusion_requests.is_empty()
            && !requested_architectures
                .contains(&PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion)
        {
            requested_architectures
                .push(PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion);
        }
        let selected_ids = request
            .transcript_ids
            .iter()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();
        let mut transcript_windows = self.derive_promoter_window_records(
            source,
            request.gene_label.as_deref(),
            None,
            request.promoter_upstream_bp,
            request.tss_proximal_downstream_bp,
            PromoterWindowCollapseMode::Transcript,
        );
        if !selected_ids.is_empty() {
            let unresolved = selected_ids
                .iter()
                .filter(|selected| {
                    !transcript_windows.iter().any(|record| {
                        Self::promoter_architecture_transcript_ids_match(
                            selected,
                            &record.transcript_id,
                        )
                    })
                })
                .cloned()
                .collect::<Vec<_>>();
            if !unresolved.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Promoter-reporter architecture request could not resolve selected transcript(s): {}",
                        unresolved.join(", ")
                    ),
                    cause_chain: vec![],
                });
            }
            transcript_windows.retain(|record| {
                selected_ids.iter().any(|selected| {
                    Self::promoter_architecture_transcript_ids_match(
                        selected,
                        &record.transcript_id,
                    )
                })
            });
        }
        if transcript_windows.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No transcript annotations matched promoter-reporter architecture request on '{}'",
                    request.seq_id
                ),
                cause_chain: vec![],
            });
        }
        let class_records = Self::collapse_promoter_window_records_by_tss_cluster(
            transcript_windows.clone(),
            request.tss_cluster_tolerance_bp,
        );
        let anchor = self.sequence_genome_anchor_summary(&request.seq_id).ok();
        let mut tss_classes = vec![];
        for class_record in class_records {
            let class_id = class_record.promoter_class_id.clone().unwrap_or_else(|| {
                format!(
                    "{}_tss_{}",
                    Self::normalize_id_token(&class_record.transcript_id),
                    class_record.tss_local_0based.saturating_add(1)
                )
            });
            let local_strand = class_record.strand.clone();
            let genomic_strand =
                Self::promoter_architecture_genomic_strand(&local_strand, anchor.as_ref());
            let class_id = Self::promoter_architecture_explicit_tss_class_id(
                &class_id,
                &local_strand,
                &genomic_strand,
            );
            let member_windows = transcript_windows
                .iter()
                .filter(|record| {
                    class_record.transcript_ids.iter().any(|id| {
                        Self::promoter_architecture_transcript_ids_match(id, &record.transcript_id)
                    })
                })
                .collect::<Vec<_>>();
            let class_evidence = request
                .initiation_evidence
                .iter()
                .filter(|evidence| {
                    class_record.transcript_ids.iter().any(|id| {
                        Self::promoter_architecture_transcript_ids_match(
                            id,
                            &evidence.transcript_id,
                        )
                    })
                })
                .cloned()
                .collect::<Vec<_>>();
            tss_classes.push(PromoterReporterTssClass {
                class_id,
                gene_label: class_record.gene_label.clone(),
                strand: local_strand.clone(),
                local_strand,
                genomic_strand,
                representative_transcript_id: class_record.transcript_id.clone(),
                representative_tss_local_0based: class_record.tss_local_0based,
                representative_tss_genomic_1based:
                    Self::promoter_architecture_local_point_to_genomic(
                        anchor.as_ref(),
                        class_record.tss_local_0based,
                    ),
                min_tss_local_0based: member_windows
                    .iter()
                    .map(|row| row.tss_local_0based)
                    .min()
                    .unwrap_or(class_record.tss_local_0based),
                max_tss_local_0based: member_windows
                    .iter()
                    .map(|row| row.tss_local_0based)
                    .max()
                    .unwrap_or(class_record.tss_local_0based),
                transcript_ids: class_record.transcript_ids.clone(),
                grouping_reason: class_record.grouping_reason.unwrap_or_else(|| {
                    "same_gene_and_strand_with_overlapping_first_exon_and_tss_within_tolerance"
                        .to_string()
                }),
                transcription_initiation_evidence_state: if class_evidence.is_empty() {
                    "annotation_only_usage_unevaluated"
                } else {
                    "caller_supplied_initiation_evidence"
                }
                .to_string(),
                initiation_evidence: class_evidence,
                recommended_usage_evidence: vec![
                    "CAGE or RAMPAGE in the relevant cell context".to_string(),
                    "Long-read 5' transcript evidence".to_string(),
                    "5' RACE or an equivalent initiation-site assay".to_string(),
                ],
            });
        }
        tss_classes.sort_by(|left, right| {
            left.min_tss_local_0based
                .cmp(&right.min_tss_local_0based)
                .then_with(|| left.class_id.cmp(&right.class_id))
        });
        let mut transcripts = vec![];
        for record in &transcript_windows {
            let class = tss_classes
                .iter()
                .find(|class| {
                    class.transcript_ids.iter().any(|id| {
                        Self::promoter_architecture_transcript_ids_match(id, &record.transcript_id)
                    })
                })
                .ok_or_else(|| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Transcript '{}' was not assigned to a TSS class",
                        record.transcript_id
                    ),
                    cause_chain: vec![],
                })?;
            transcripts.push(self.promoter_architecture_transcript_audit(
                source,
                record,
                &class.class_id,
                anchor.as_ref(),
            )?);
        }
        transcripts.sort_by(|left, right| {
            left.tss_local_0based
                .cmp(&right.tss_local_0based)
                .then_with(|| left.transcript_id.cmp(&right.transcript_id))
        });
        let unresolved_fusion_targets = request
            .fusion_requests
            .iter()
            .filter(|fusion| {
                !transcripts.iter().any(|audit| {
                    Self::promoter_architecture_transcript_ids_match(
                        &fusion.transcript_id,
                        &audit.transcript_id,
                    )
                })
            })
            .map(|fusion| fusion.transcript_id.trim().to_string())
            .collect::<BTreeSet<_>>();
        if !unresolved_fusion_targets.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Promoter-reporter fusion request could not resolve transcript target(s): {}",
                    unresolved_fusion_targets
                        .into_iter()
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
                cause_chain: vec![],
            });
        }
        let mut architectures = vec![];
        for audit in &transcripts {
            let class = tss_classes
                .iter()
                .find(|class| class.class_id == audit.tss_class_id)
                .expect("validated TSS class assignment");
            for kind in &requested_architectures {
                if *kind == PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion {
                    for fusion in request.fusion_requests.iter().filter(|fusion| {
                        Self::promoter_architecture_transcript_ids_match(
                            &fusion.transcript_id,
                            &audit.transcript_id,
                        )
                    }) {
                        architectures.push(Self::promoter_architecture_row(
                            source.len(),
                            class,
                            audit,
                            *kind,
                            request.promoter_upstream_bp,
                            request.tss_proximal_downstream_bp,
                            Some(fusion),
                        )?);
                    }
                } else {
                    architectures.push(Self::promoter_architecture_row(
                        source.len(),
                        class,
                        audit,
                        *kind,
                        request.promoter_upstream_bp,
                        request.tss_proximal_downstream_bp,
                        None,
                    )?);
                }
            }
        }
        architectures.sort_by(|left, right| {
            left.tss_class_id
                .cmp(&right.tss_class_id)
                .then_with(|| left.transcript_id.cmp(&right.transcript_id))
                .then_with(|| {
                    left.architecture_kind
                        .as_str()
                        .cmp(right.architecture_kind.as_str())
                })
                .then_with(|| left.architecture_id.cmp(&right.architecture_id))
        });
        architectures.dedup_by(|left, right| left.architecture_id == right.architecture_id);
        let common_cds_local = transcripts
            .iter()
            .filter_map(|audit| {
                if audit.strand == "-" {
                    audit
                        .cds_start_codon_source_ranges_0based
                        .iter()
                        .map(|(_, end)| end.saturating_sub(1))
                        .max()
                } else {
                    audit
                        .cds_start_codon_source_ranges_0based
                        .iter()
                        .map(|(start, _)| *start)
                        .min()
                }
            })
            .collect::<BTreeSet<_>>();
        let common_cds_start_local_0based =
            (common_cds_local.len() == 1).then(|| *common_cds_local.iter().next().unwrap());
        let common_cds_start_genomic_1based = common_cds_start_local_0based.and_then(|local| {
            Self::promoter_architecture_local_point_to_genomic(anchor.as_ref(), local)
        });
        let existing_construct =
            self.promoter_architecture_existing_construct(request.existing_construct.as_ref())?;
        let staged_panel = Self::promoter_architecture_staged_panel(
            &tss_classes,
            &architectures,
            &existing_construct,
        );
        let cutrun_evidence = self.promoter_architecture_cutrun_evidence(&request);
        let theoretical_motif_hits = if request.motif_tokens.is_empty() {
            None
        } else {
            Some(self.scan_tfbs_hits(
                SequenceScanTarget::SeqId {
                    seq_id: request.seq_id.clone(),
                    span_start_0based: None,
                    span_end_0based_exclusive: None,
                },
                &request.motif_tokens,
                None,
                request.motif_min_llr_quantile,
                &[],
                Some(2_000),
                None,
                None,
            )?)
        };
        let locus_evidence = request
            .locus_evidence_request
            .as_ref()
            .map(|locus_request| {
                let mut locus_request = locus_request.clone();
                if !request.include_local_source_paths {
                    locus_request.include_local_source_paths = false;
                }
                self.build_gene_locus_evidence_display_report(&request.seq_id, &locus_request)
            })
            .transpose()?;
        let mut warnings = vec![];
        if tss_classes.len() > 1 && request.initiation_evidence.is_empty() {
            warnings.push(
                "Multiple annotation-derived TSS classes were found, but no cell-context transcription-initiation evidence was supplied. Alternative promoter usage remains unevaluated."
                    .to_string(),
            );
        }
        if existing_construct.status == "existing_construct_unverified" {
            warnings.push(
                "The existing laboratory construct is unverified and was not assigned to any architecture row."
                    .to_string(),
            );
        }
        warnings.extend(cutrun_evidence.warnings.iter().cloned());
        let source_provenance = self.promoter_architecture_source_provenance(&request, source);
        Ok(PromoterReporterArchitectureComparisonReport {
            schema: PROMOTER_REPORTER_ARCHITECTURE_REPORT_SCHEMA.to_string(),
            request_schema: request.schema,
            seq_id: request.seq_id,
            gene_label: request.gene_label,
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            source_provenance,
            promoter_upstream_bp: request.promoter_upstream_bp,
            tss_proximal_downstream_bp: request.tss_proximal_downstream_bp,
            tss_cluster_tolerance_bp: request.tss_cluster_tolerance_bp,
            tss_classes,
            transcripts,
            architectures,
            common_cds_start_local_0based,
            common_cds_start_genomic_1based,
            existing_construct,
            staged_panel,
            cutrun_evidence,
            theoretical_motif_hits,
            locus_evidence,
            warnings,
            non_claims: vec![
                "Transcript annotation alone does not establish which TSS is used in the assayed cells."
                    .to_string(),
                "CUT&RUN overlap is occupancy/enrichment evidence, not proof of TSS usage, promoter activity, direct regulation, biochemical affinity, or luciferase output."
                    .to_string(),
                "Reporter geometry does not verify the identity of a physical laboratory construct."
                    .to_string(),
            ],
            json_path: None,
            svg_path: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn feature(
        kind: &str,
        ranges: &[(usize, usize)],
        transcript_id: &str,
        gene: &str,
        reverse: bool,
    ) -> gb_io::seq::Feature {
        let parts = ranges
            .iter()
            .map(|(start, end)| gb_io::seq::Location::simple_range(*start as i64, *end as i64))
            .collect::<Vec<_>>();
        let location = if parts.len() == 1 {
            parts.into_iter().next().expect("one feature range")
        } else {
            gb_io::seq::Location::Join(parts)
        };
        gb_io::seq::Feature {
            kind: kind.to_string().into(),
            location: if reverse {
                gb_io::seq::Location::Complement(Box::new(location))
            } else {
                location
            },
            qualifiers: vec![
                ("transcript_id".into(), Some(transcript_id.to_string())),
                ("gene".into(), Some(gene.to_string())),
                ("label".into(), Some(transcript_id.to_string())),
                ("transcript_support_level".into(), Some("1".to_string())),
            ],
        }
    }

    fn plus_locus_engine() -> GentleEngine {
        let mut bases = vec![b'C'; 800];
        bases[360..363].copy_from_slice(b"ATG");
        let mut dna = DNAsequence::from_sequence(&String::from_utf8(bases).unwrap())
            .expect("synthetic plus locus");
        for (transcript, exons) in [
            ("SYN1-201", vec![(100, 180), (300, 500)]),
            ("SYN1-202", vec![(105, 500)]),
            ("SYN1-203", vec![(20, 80), (120, 160), (300, 500)]),
        ] {
            dna.features_mut()
                .push(feature("mRNA", &exons, transcript, "SYN1", false));
            for exon in &exons {
                dna.features_mut()
                    .push(feature("exon", &[*exon], transcript, "SYN1", false));
            }
            dna.features_mut()
                .push(feature("CDS", &[(360, 500)], transcript, "SYN1", false));
        }
        GentleEngine::prepare_sequence(&mut dna);
        let mut state = ProjectState::default();
        state.sequences.insert("synthetic_plus".to_string(), dna);
        GentleEngine::from_state(state)
    }

    fn minus_locus_engine() -> GentleEngine {
        let mut bases = vec![b'C'; 800];
        bases[497..500].copy_from_slice(b"CAT");
        let mut dna = DNAsequence::from_sequence(&String::from_utf8(bases).unwrap())
            .expect("synthetic minus locus");
        let exons = vec![(400, 550), (600, 700)];
        dna.features_mut()
            .push(feature("mRNA", &exons, "SYNREV-201", "SYNREV", true));
        for exon in &exons {
            dna.features_mut()
                .push(feature("exon", &[*exon], "SYNREV-201", "SYNREV", true));
        }
        dna.features_mut()
            .push(feature("CDS", &[(400, 500)], "SYNREV-201", "SYNREV", true));
        GentleEngine::prepare_sequence(&mut dna);
        let mut state = ProjectState::default();
        state.sequences.insert("synthetic_minus".to_string(), dna);
        GentleEngine::from_state(state)
    }

    fn serpine1_fixture_engine() -> GentleEngine {
        let geometry: serde_json::Value = serde_json::from_str(include_str!(
            "../../../test_files/fixtures/promoter_reporter_architecture/ensembl_116_serpine1_geometry.json"
        ))
        .expect("SERPINE1 geometry fixture");
        let fasta = include_str!(
            "../../../test_files/fixtures/promoter_reporter_architecture/ensembl_116_serpine1_region.fa"
        );
        let sequence = fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .collect::<String>();
        assert_eq!(
            crate::digest_utils::sha256_hex_str(&sequence),
            geometry["provenance"]["sequence_response_sha256"]
                .as_str()
                .expect("sequence response SHA-256")
        );
        let region_start = geometry["fixture_region"]["start_1based"]
            .as_u64()
            .expect("fixture region start") as usize;
        let region_end = geometry["fixture_region"]["end_1based"]
            .as_u64()
            .expect("fixture region end") as usize;
        let mut dna = DNAsequence::from_sequence(&sequence).expect("SERPINE1 fixture sequence");
        let transcripts = geometry["transcripts"]
            .as_array()
            .expect("SERPINE1 transcripts");
        assert_eq!(transcripts.len(), 15);
        let mut panel_isoforms = Vec::with_capacity(transcripts.len());
        for transcript in transcripts {
            let transcript_id = transcript["id"].as_str().expect("versioned transcript id");
            let display_name = transcript["display_name"]
                .as_str()
                .expect("transcript display name");
            let exons = transcript["exons"]
                .as_array()
                .expect("transcript exons")
                .iter()
                .map(|exon| {
                    let start = exon["start_1based"].as_u64().expect("exon start") as usize;
                    let end = exon["end_1based"].as_u64().expect("exon end") as usize;
                    (start - region_start, end - region_start + 1)
                })
                .collect::<Vec<_>>();
            let mut mrna = feature("mRNA", &exons, transcript_id, "SERPINE1", false);
            mrna.qualifiers.extend([
                ("gene_id".into(), Some("ENSG00000106366.11".to_string())),
                ("product".into(), Some(display_name.to_string())),
                (
                    "biotype".into(),
                    transcript["biotype"].as_str().map(str::to_string),
                ),
            ]);
            if transcript["is_canonical"].as_u64() == Some(1) {
                mrna.qualifiers
                    .push(("tag".into(), Some("Ensembl_canonical".to_string())));
            }
            dna.features_mut().push(mrna);
            for exon in &exons {
                dna.features_mut().push(feature(
                    "exon",
                    &[*exon],
                    transcript_id,
                    "SERPINE1",
                    false,
                ));
            }
            let translation_start = transcript["translation_start_1based"]
                .as_u64()
                .expect("translation start") as usize;
            let translation_end = transcript["translation_end_1based"]
                .as_u64()
                .expect("translation end") as usize;
            let cds_ranges = transcript["exons"]
                .as_array()
                .expect("transcript exons")
                .iter()
                .filter_map(|exon| {
                    let exon_start = exon["start_1based"].as_u64()? as usize;
                    let exon_end = exon["end_1based"].as_u64()? as usize;
                    let start = exon_start.max(translation_start);
                    let end = exon_end.min(translation_end);
                    (end >= start).then_some((start - region_start, end - region_start + 1))
                })
                .collect::<Vec<_>>();
            dna.features_mut().push(feature(
                "CDS",
                &cds_ranges,
                transcript_id,
                "SERPINE1",
                false,
            ));
            panel_isoforms.push(IsoformPanelIsoformSpec {
                isoform_id: transcript_id.to_string(),
                label: Some(display_name.to_string()),
                transcript_ids: vec![transcript_id.to_string()],
                annotation_transcript_id: Some(transcript_id.to_string()),
                transcript_strand: Some("+".to_string()),
                transcript_biotype: transcript["biotype"].as_str().map(str::to_string),
                annotation_is_canonical: Some(transcript["is_canonical"].as_u64() == Some(1)),
                exon_ranges_genomic_1based: transcript["exons"]
                    .as_array()
                    .expect("transcript exons")
                    .iter()
                    .map(|exon| IsoformPanelAnnotationRange {
                        start_1based: exon["start_1based"].as_u64().expect("exon start") as usize,
                        end_1based: exon["end_1based"].as_u64().expect("exon end") as usize,
                        exon_id: exon["id"].as_str().map(str::to_string),
                        ..Default::default()
                    })
                    .collect(),
                cds_ranges_genomic_1based: cds_ranges
                    .iter()
                    .map(|(start, end)| IsoformPanelAnnotationRange {
                        start_1based: region_start + start,
                        end_1based: region_start + end - 1,
                        ..Default::default()
                    })
                    .collect(),
                coding_start_status: Some("annotation_translation_boundary".to_string()),
                coding_completeness_status: Some("not_reported_by_ensembl_lookup".to_string()),
                annotation_provenance: Some(
                    "Pinned Ensembl 116 expanded SERPINE1 lookup projection".to_string(),
                ),
                ..Default::default()
            });
        }
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "track".into(),
            location: gb_io::seq::Location::simple_range(6_900, 7_300),
            qualifiers: vec![
                (
                    "label".into(),
                    Some("synthetic SERPINE1 occupancy".to_string()),
                ),
                ("score".into(), Some("9.5".to_string())),
                ("gentle_track_source".into(), Some("BED".to_string())),
                (
                    "gentle_generated".into(),
                    Some("genome_bed_track".to_string()),
                ),
                (
                    "gentle_track_name".into(),
                    Some("Synthetic SERPINE1 CUT&RUN".to_string()),
                ),
                (
                    "gentle_track_file".into(),
                    Some("synthetic://serpine1-occupancy".to_string()),
                ),
            ],
        });
        GentleEngine::prepare_sequence(&mut dna);
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("serpine1_ensembl116".to_string(), dna);
        state.metadata.insert(
            "provenance".to_string(),
            serde_json::json!({
                "genome_extractions": [{
                    "seq_id": "serpine1_ensembl116",
                    "recorded_at_unix_ms": 1,
                    "operation": "pinned_ensembl_rest_fixture",
                    "genome_id": "GRCh38 Ensembl 116",
                    "catalog_path": "test_files/fixtures/promoter_reporter_architecture",
                    "cache_dir": null,
                    "chromosome": "7",
                    "start_1based": region_start,
                    "end_1based": region_end,
                    "gene_query": "ENSG00000106366",
                    "occurrence": 1,
                    "gene_extract_mode": "gene_with_flanks",
                    "transcript_id": null,
                    "tss_1based": null,
                    "promoter_upstream_bp": 1000,
                    "promoter_downstream_bp": 0,
                    "gene_id": "ENSG00000106366.11",
                    "gene_name": "SERPINE1",
                    "strand": "+",
                    "anchor_strand": "+",
                    "anchor_verified": true,
                    "sequence_source_type": "pinned_ensembl_rest_region",
                    "annotation_source_type": "pinned_ensembl_rest_lookup_projection",
                    "sequence_source": geometry["provenance"]["sequence_url"],
                    "annotation_source": geometry["provenance"]["lookup_url"],
                    "sequence_sha1": "19d7e55abec1a990384b008cf6e05b32ae4bc022",
                    "annotation_sha1": "3780869ba4fd35e654ad27a5b9be332a3011d232"
                }]
            }),
        );
        let mut engine = GentleEngine::from_state(state);
        engine
            .upsert_isoform_panel_record(IsoformPanelRecord {
                seq_id: "serpine1_ensembl116".to_string(),
                panel_id: "serpine1_ensembl116_panel".to_string(),
                imported_at_unix_ms: 1,
                source_path:
                    "test_files/fixtures/promoter_reporter_architecture/ensembl_116_serpine1_geometry.json"
                        .to_string(),
                strict: false,
                resource: IsoformPanelResource {
                    schema: ISOFORM_PANEL_RESOURCE_SCHEMA.to_string(),
                    panel_id: "serpine1_ensembl116_panel".to_string(),
                    gene_symbol: "SERPINE1".to_string(),
                    assembly: Some("GRCh38 Ensembl 116".to_string()),
                    source: Some("Pinned Ensembl 116 acceptance fixture".to_string()),
                    isoforms: panel_isoforms,
                    ..Default::default()
                },
            })
            .expect("register pinned SERPINE1 panel");
        engine
    }

    #[test]
    fn legacy_architecture_payloads_default_to_no_locus_evidence() {
        let request: PromoterReporterArchitectureComparisonRequest =
            serde_json::from_value(serde_json::json!({
                "schema": PROMOTER_REPORTER_ARCHITECTURE_REQUEST_SCHEMA,
                "seq_id": "legacy_locus"
            }))
            .expect("legacy architecture request");
        assert!(request.locus_evidence_request.is_none());
        assert!(
            !serde_json::to_value(&request)
                .expect("serialize request")
                .as_object()
                .expect("request object")
                .contains_key("locus_evidence_request")
        );

        let report: PromoterReporterArchitectureComparisonReport =
            serde_json::from_value(serde_json::json!({
                "schema": PROMOTER_REPORTER_ARCHITECTURE_REPORT_SCHEMA,
                "request_schema": PROMOTER_REPORTER_ARCHITECTURE_REQUEST_SCHEMA,
                "seq_id": "legacy_locus"
            }))
            .expect("legacy architecture report");
        assert!(report.locus_evidence.is_none());
    }

    #[test]
    fn comparison_separates_tss_classes_and_genomic_from_spliced_leaders() {
        let report = plus_locus_engine()
            .compare_promoter_reporter_architectures(
                PromoterReporterArchitectureComparisonRequest {
                    seq_id: "synthetic_plus".to_string(),
                    gene_label: Some("SYN1".to_string()),
                    promoter_upstream_bp: 20,
                    tss_proximal_downstream_bp: 10,
                    ..Default::default()
                },
            )
            .expect("compare plus-strand architectures");

        assert_eq!(report.tss_classes.len(), 2);
        assert_eq!(report.transcripts.len(), 3);
        assert_eq!(report.architectures.len(), 9);
        assert!(report.transcripts.iter().all(|row| {
            row.strand == "+" && row.local_strand == "+" && row.genomic_strand == "+"
        }));
        assert!(
            report
                .tss_classes
                .iter()
                .all(|class| { class.class_id.contains("_local_plus_genomic_plus_tss_") })
        );
        assert_eq!(report.common_cds_start_local_0based, Some(360));
        let intronless = report
            .transcripts
            .iter()
            .find(|row| row.transcript_id == "SYN1-202")
            .expect("intronless transcript");
        assert_eq!(intronless.five_prime_utr_total_intron_bp, 0);
        let long_leader = report
            .transcripts
            .iter()
            .find(|row| row.transcript_id == "SYN1-203")
            .expect("multi-intron transcript");
        assert_eq!(
            long_leader
                .five_prime_utr_introns
                .iter()
                .map(|row| (row.start_0based, row.end_0based_exclusive))
                .collect::<Vec<_>>(),
            vec![(80, 120), (160, 300)]
        );
        assert_eq!(
            report.cutrun_evidence.evaluation_state,
            "unevaluated_no_sources_selected"
        );
        assert!(report.warnings.iter().any(|warning| {
            warning.contains("no cell-context transcription-initiation evidence")
        }));
    }

    #[test]
    fn ensembl_116_serpine1_fixture_matches_observed_representative_geometry() {
        let mut report = serpine1_fixture_engine()
            .compare_promoter_reporter_architectures(
                PromoterReporterArchitectureComparisonRequest {
                    seq_id: "serpine1_ensembl116".to_string(),
                    gene_label: Some("SERPINE1".to_string()),
                    transcript_ids: vec![
                        "ENST00000223095.5".to_string(),
                        "ENST00000870828.2".to_string(),
                        "ENST00000950058.1".to_string(),
                    ],
                    promoter_upstream_bp: 1_000,
                    tss_proximal_downstream_bp: 200,
                    tss_cluster_tolerance_bp: 600,
                    locus_evidence_request: Some(GeneLocusEvidenceDisplayRequest {
                        isoform_evidence: GeneIsoformEvidenceRequest {
                            panel_id: "serpine1_ensembl116_panel".to_string(),
                            annotation_release: Some("Ensembl 116".to_string()),
                            ..Default::default()
                        },
                        upstream_bp: 1_000,
                        downstream_bp: 500,
                        occupancy_layout: GeneLocusOccupancyLayout {
                            schema: GENE_LOCUS_OCCUPANCY_LAYOUT_SCHEMA.to_string(),
                            groups: vec![GeneLocusOccupancyGroupRequest {
                                group_id: "synthetic_serpine1_occupancy".to_string(),
                                label: "Synthetic SERPINE1 evidence".to_string(),
                                scale_mode: GeneLocusOccupancyScaleMode::Independent,
                                lanes: vec![GeneLocusOccupancyLaneRequest {
                                    track_name: "Synthetic SERPINE1 CUT&RUN".to_string(),
                                    source_id: Some(
                                        "synthetic:serpine1:cutrun:acceptance".to_string(),
                                    ),
                                    source_assembly: Some("GRCh38 Ensembl 116".to_string()),
                                    assay: Some("CUT&RUN".to_string()),
                                    factor: Some("SERPINE1".to_string()),
                                    role: GeneLocusOccupancyLaneRole::Experimental,
                                    ..Default::default()
                                }],
                                ..Default::default()
                            }],
                        },
                        scale_bar: GeneLocusScaleBarPolicy {
                            mode: GeneLocusScaleBarMode::Auto,
                            ..Default::default()
                        },
                        ..Default::default()
                    }),
                    ..Default::default()
                },
            )
            .expect("compare pinned Ensembl 116 SERPINE1 representatives");

        let tss = report.tss_classes[0].representative_tss_local_0based;
        let tss_genomic = report.tss_classes[0]
            .representative_tss_genomic_1based
            .expect("SERPINE1 genomic TSS");
        let tss_class_id = report.tss_classes[0].class_id.clone();
        let architecture_id = report.architectures[0].architecture_id.clone();
        let source = crate::ensembl_regulation::source_descriptor(
            crate::ensembl_regulation::ENSEMBL_REGULATION_2026_08_GRCH38_SOURCE_ID,
        )
        .expect("pinned Ensembl Regulation source");
        report
            .locus_evidence
            .as_mut()
            .expect("SERPINE1 locus evidence")
            .ensembl_regulation = Some(GeneLocusEnsemblRegulationEvidence {
            schema: GENE_LOCUS_ENSEMBL_REGULATION_EVIDENCE_SCHEMA.to_string(),
            availability: GeneLocusEnsemblRegulationAvailability::Available,
            requested_source_id: source.source_id.clone(),
            source: Some(source.clone()),
            source_binding: Some(GeneLocusEnsemblRegulationSourceBinding {
                source: source.clone(),
                overlap_report_id: "synthetic_serpine1_ensembl_overlap".to_string(),
                index_sha256: "sha256:synthetic-serpine1-index".to_string(),
                intervals_sha256: "sha256:synthetic-serpine1-intervals".to_string(),
                content_identity_verified: true,
                assembly_match_status: "assembly_and_species_matched_source_release_pinned"
                    .to_string(),
                requested_feature_types: vec!["promoter".to_string()],
                max_rows: 10,
                matched_feature_count: 1,
                returned_feature_count: 1,
                ..Default::default()
            }),
            rows: vec![GeneLocusEnsemblRegulationFeatureRow {
                source_id: source.source_id.clone(),
                provider: source.provider.clone(),
                annotation_release: source.annotation_release.clone(),
                annotation_api_version: source.annotation_api_version.clone(),
                pipeline_version: source.pipeline_version.clone(),
                assembly_name: source.assembly_name.clone(),
                assembly_accession: source.assembly_accession.clone(),
                feature_id: "ENSR_TEST_SERPINE1_PROMOTER".to_string(),
                feature_type: "promoter".to_string(),
                core_genomic_start_1based: tss_genomic.saturating_sub(40),
                core_genomic_end_1based: tss_genomic.saturating_add(40),
                displayed_genomic_start_1based: tss_genomic.saturating_sub(40),
                displayed_genomic_end_1based: tss_genomic.saturating_add(40),
                displayed_local_start_1based: tss.saturating_sub(40).saturating_add(1),
                displayed_local_end_1based: tss.saturating_add(41),
                local_strand: ".".to_string(),
                genomic_strand: ".".to_string(),
                locus_relation: "contained_in_displayed_locus".to_string(),
                associated_gene_ids: vec!["ENSG00000106366".to_string()],
                associated_gene_names: vec!["SERPINE1".to_string()],
                canonical_feature_url: crate::ensembl_regulation::feature_page_url(
                    &source,
                    "ENSR_TEST_SERPINE1_PROMOTER",
                )
                .expect("synthetic provider URL"),
                feature_url_template: source.feature_page_url_template.clone(),
                overlapping_tss_class_ids: vec![tss_class_id.clone()],
                nearest_tss_class_ids: vec![tss_class_id],
                signed_distance_to_nearest_tss_bp: Some(0),
                absolute_distance_to_nearest_tss_bp: Some(0),
                signed_distance_basis: "negative=upstream; positive=downstream; zero=contains TSS"
                    .to_string(),
                overlapping_reporter_architecture_ids: vec![architecture_id],
                relationship_statement: "Synthetic geometric overlap for renderer acceptance; provider-associated genes remain annotations, not causal assignments."
                    .to_string(),
                ..Default::default()
            }],
            evidence_statement: "Synthetic Ensembl Regulation overlap for public SERPINE1 renderer acceptance."
                .to_string(),
            non_claims: vec![
                "This synthetic overlap does not establish activity or causal regulation."
                    .to_string(),
            ],
            provider_link_note: "Link opens the live provider page; plotted evidence is synthetic and release-bound for testing."
                .to_string(),
            ..Default::default()
        });

        assert_eq!(report.tss_classes.len(), 2);
        assert_eq!(report.transcripts.len(), 3);
        assert_eq!(report.common_cds_start_genomic_1based, Some(101_128_394));
        for (transcript_id, span, spliced_utr, introns, max_intron) in [
            ("ENST00000223095.5", 1_293, 142, 1_148, 1_148),
            ("ENST00000870828.2", 1_796, 753, 1_040, 852),
            ("ENST00000950058.1", 6_239, 1_007, 5_229, 3_122),
        ] {
            let audit = report
                .transcripts
                .iter()
                .find(|row| row.transcript_id == transcript_id)
                .unwrap_or_else(|| panic!("missing transcript audit {transcript_id}"));
            assert_eq!(audit.tss_to_atg_through_span_bp, span, "{transcript_id}");
            assert_eq!(
                audit.five_prime_utr_spliced_length_bp, spliced_utr,
                "{transcript_id}"
            );
            assert_eq!(
                audit.five_prime_utr_total_intron_bp, introns,
                "{transcript_id}"
            );
            assert_eq!(
                audit.five_prime_utr_max_intron_bp, max_intron,
                "{transcript_id}"
            );
        }
        let serpine1_202 = report
            .transcripts
            .iter()
            .find(|row| row.transcript_id == "ENST00000870828.2")
            .expect("SERPINE1-202 transcript audit");
        assert_eq!(
            serpine1_202
                .five_prime_utr_introns
                .iter()
                .map(|intron| intron.length_bp)
                .collect::<Vec<_>>(),
            vec![852, 188]
        );
        assert_eq!(
            report
                .source_provenance
                .genome_anchor
                .as_ref()
                .map(|row| row.genome_id.as_str()),
            Some("GRCh38 Ensembl 116")
        );
        assert_eq!(
            report.source_provenance.sequence_sha1.as_deref(),
            Some("19d7e55abec1a990384b008cf6e05b32ae4bc022")
        );
        assert_eq!(
            report.source_provenance.annotation_sha1.as_deref(),
            Some("3780869ba4fd35e654ad27a5b9be332a3011d232")
        );
        let locus_evidence = report
            .locus_evidence
            .as_ref()
            .expect("SERPINE1 reporter comparison should carry the composed locus evidence");
        let occupancy_lane = &locus_evidence.occupancy_groups[0].lanes[0];
        assert_eq!(occupancy_lane.state, GeneLocusOccupancyLaneState::Available);
        assert_eq!(occupancy_lane.assay.as_deref(), Some("CUT&RUN"));
        let combined_svg =
            crate::render_promoter_reporter_architecture::render_promoter_reporter_architecture_svg(
                &report,
            );
        assert!(
            combined_svg
                .contains("data-gentle-overlay-id=\"promoter_reporter_architecture_comparison\"")
        );
        assert!(
            combined_svg
                .contains("data-gentle-occupancy-source=\"synthetic:serpine1:cutrun:acceptance\"")
        );
        assert!(combined_svg.contains("data-gentle-transcript=\"ENST00000223095.5\""));
        assert!(combined_svg.contains("Ensembl-annotated regulatory regions"));
        assert!(combined_svg.contains("ENSR_TEST_SERPINE1_PROMOTER"));
        let linked =
            crate::render_promoter_reporter_architecture::render_promoter_reporter_architecture_with_links(
                &report,
            );
        assert_eq!(linked.uri_links.len(), 2);
        if let Some(output_dir) = std::env::var_os("GENTLE_SERPINE1_REPORTER_ARTIFACT_DIR") {
            let output_dir = std::path::PathBuf::from(output_dir);
            std::fs::create_dir_all(&output_dir).expect("create SERPINE1 artifact directory");
            let mut artifact_report = report.clone();
            artifact_report.generated_at_unix_ms = 0;
            artifact_report.op_id = None;
            artifact_report.run_id = None;
            artifact_report.json_path = Some("serpine1_architecture_comparison.json".to_string());
            artifact_report.svg_path = Some("serpine1_architecture_comparison.svg".to_string());
            std::fs::write(
                output_dir.join("serpine1_architecture_comparison.json"),
                serde_json::to_string_pretty(&artifact_report)
                    .expect("serialize SERPINE1 architecture artifact"),
            )
            .expect("write SERPINE1 architecture JSON artifact");
            std::fs::write(
                output_dir.join("serpine1_architecture_comparison.svg"),
                crate::render_promoter_reporter_architecture::render_promoter_reporter_architecture_svg(
                    &artifact_report,
                ),
            )
            .expect("write SERPINE1 architecture SVG artifact");
        }
    }

    #[test]
    fn comparison_is_strand_aware_and_keeps_atg_policy_explicit() {
        let report = minus_locus_engine()
            .compare_promoter_reporter_architectures(
                PromoterReporterArchitectureComparisonRequest {
                    seq_id: "synthetic_minus".to_string(),
                    gene_label: Some("SYNREV".to_string()),
                    promoter_upstream_bp: 50,
                    tss_proximal_downstream_bp: 10,
                    ..Default::default()
                },
            )
            .expect("compare minus-strand architectures");

        assert_eq!(report.common_cds_start_local_0based, Some(499));
        assert_eq!(report.transcripts[0].strand, "-");
        assert_eq!(report.transcripts[0].local_strand, "-");
        assert_eq!(report.transcripts[0].genomic_strand, "-");
        assert_eq!(
            report.transcripts[0].cds_start_codon_5prime_to_3prime,
            "ATG"
        );
        let proximal = report
            .architectures
            .iter()
            .find(|row| {
                row.architecture_kind
                    == PromoterReporterArchitectureKind::TssProximalTranscriptional
            })
            .expect("minus-strand proximal row");
        assert_eq!(
            proximal.junction.endogenous_atg_disposition,
            PromoterReporterAtgDisposition::EndogenousExcludedVectorLuciferaseRetained
        );
    }

    #[test]
    fn retained_atg_fusions_require_and_audit_explicit_semantics() {
        let engine = plus_locus_engine();
        let missing = engine.compare_promoter_reporter_architectures(
            PromoterReporterArchitectureComparisonRequest {
                seq_id: "synthetic_plus".to_string(),
                transcript_ids: vec!["SYN1-201".to_string()],
                architectures: vec![PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion],
                ..Default::default()
            },
        );
        assert!(missing.is_err());

        let report = engine
            .compare_promoter_reporter_architectures(
                PromoterReporterArchitectureComparisonRequest {
                    seq_id: "synthetic_plus".to_string(),
                    transcript_ids: vec!["SYN1-201".to_string()],
                    architectures: vec![
                        PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion,
                    ],
                    fusion_requests: vec![PromoterReporterFusionRequest {
                        transcript_id: "SYN1-201".to_string(),
                        mode: PromoterReporterFusionMode::InFrameTranslationFusion,
                        retained_endogenous_cds_bp: 6,
                        vector_luciferase_atg_removed: Some(true),
                        junction_sequence_5prime_to_3prime: Some("ATGGCCGGT".to_string()),
                    }],
                    ..Default::default()
                },
            )
            .expect("audited in-frame fusion");
        assert_eq!(report.architectures.len(), 1);
        assert_eq!(report.architectures[0].junction.audit_status, "pass");
        assert_eq!(
            report.architectures[0].junction.endogenous_atg_disposition,
            PromoterReporterAtgDisposition::EndogenousRetainedVectorLuciferaseRemoved
        );
    }

    #[test]
    fn comparison_rejects_unresolved_targets_and_distinguishes_fusion_junctions() {
        let engine = plus_locus_engine();
        let unresolved = engine.compare_promoter_reporter_architectures(
            PromoterReporterArchitectureComparisonRequest {
                seq_id: "synthetic_plus".to_string(),
                transcript_ids: vec!["SYN1-201".to_string(), "missing-transcript".to_string()],
                ..Default::default()
            },
        );
        assert!(
            unresolved
                .expect_err("unresolved selected transcript must fail closed")
                .message
                .contains("missing-transcript")
        );

        let fusion = |junction: &str| PromoterReporterFusionRequest {
            transcript_id: "SYN1-201".to_string(),
            mode: PromoterReporterFusionMode::InFrameTranslationFusion,
            retained_endogenous_cds_bp: 6,
            vector_luciferase_atg_removed: Some(true),
            junction_sequence_5prime_to_3prime: Some(junction.to_string()),
        };
        let report = engine
            .compare_promoter_reporter_architectures(
                PromoterReporterArchitectureComparisonRequest {
                    seq_id: "synthetic_plus".to_string(),
                    transcript_ids: vec!["SYN1-201".to_string()],
                    architectures: vec![
                        PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion,
                    ],
                    fusion_requests: vec![
                        fusion("ATGGCCGGT"),
                        fusion("ATGGCCAAA"),
                        fusion("ATGGCCGGT"),
                    ],
                    ..Default::default()
                },
            )
            .expect("compare two distinct fusion junctions");
        assert_eq!(report.architectures.len(), 2);
        assert_ne!(
            report.architectures[0].architecture_id,
            report.architectures[1].architecture_id
        );

        let unresolved_fusion = engine.compare_promoter_reporter_architectures(
            PromoterReporterArchitectureComparisonRequest {
                seq_id: "synthetic_plus".to_string(),
                transcript_ids: vec!["SYN1-201".to_string()],
                fusion_requests: vec![PromoterReporterFusionRequest {
                    transcript_id: "missing-fusion-target".to_string(),
                    mode: PromoterReporterFusionMode::InFrameTranslationFusion,
                    retained_endogenous_cds_bp: 6,
                    vector_luciferase_atg_removed: Some(true),
                    junction_sequence_5prime_to_3prime: Some("ATGGCCGGT".to_string()),
                }],
                ..Default::default()
            },
        );
        assert!(
            unresolved_fusion
                .expect_err("unresolved fusion target must fail closed")
                .message
                .contains("missing-fusion-target")
        );
    }

    #[test]
    fn operation_exports_portable_json_svg_and_normalized_analysis_is_deterministic() {
        fn normalized_json(
            report: &PromoterReporterArchitectureComparisonReport,
        ) -> serde_json::Value {
            let mut value = serde_json::to_value(report).expect("serialize comparison report");
            let object = value.as_object_mut().expect("comparison report object");
            object.remove("generated_at_unix_ms");
            object.remove("op_id");
            object.remove("run_id");
            object.remove("json_path");
            object.remove("svg_path");
            value
        }

        let request = PromoterReporterArchitectureComparisonRequest {
            seq_id: "synthetic_plus".to_string(),
            gene_label: Some("SYN1".to_string()),
            transcript_ids: vec!["SYN1-201".to_string(), "SYN1-203".to_string()],
            promoter_upstream_bp: 20,
            tss_proximal_downstream_bp: 10,
            ..Default::default()
        };
        let mut engine = plus_locus_engine();
        let first = engine
            .compare_promoter_reporter_architectures(request.clone())
            .expect("first deterministic comparison");
        assert!(
            first.locus_evidence.is_none(),
            "legacy requests must keep the original renderer contract"
        );
        let second = engine
            .compare_promoter_reporter_architectures(request.clone())
            .expect("second deterministic comparison");
        assert_eq!(normalized_json(&first), normalized_json(&second));
        assert_eq!(
            crate::render_promoter_reporter_architecture::render_promoter_reporter_architecture_svg(
                &first
            ),
            crate::render_promoter_reporter_architecture::render_promoter_reporter_architecture_svg(
                &second
            )
        );

        let temp = tempfile::tempdir().expect("architecture export tempdir");
        let json_path = temp.path().join("comparison.json");
        let svg_path = temp.path().join("comparison.svg");
        let result = engine
            .apply(Operation::ComparePromoterReporterArchitectures {
                request: Box::new(request),
                path: Some(json_path.to_string_lossy().into_owned()),
                svg_path: Some(svg_path.to_string_lossy().into_owned()),
            })
            .expect("execute architecture comparison operation");
        assert!(result.promoter_reporter_architecture_comparison.is_some());
        let json = std::fs::read_to_string(json_path).expect("read exported comparison JSON");
        let svg = std::fs::read_to_string(svg_path).expect("read exported comparison SVG");
        assert!(json.contains(PROMOTER_REPORTER_ARCHITECTURE_REPORT_SCHEMA));
        assert!(svg.contains("Proposed reporter architectures"));
        assert!(svg.contains("Occupancy does not establish promoter use"));
    }

    #[test]
    fn unresolved_cutrun_source_remains_unevaluated_and_paths_stay_private() {
        let report = plus_locus_engine()
            .compare_promoter_reporter_architectures(
                PromoterReporterArchitectureComparisonRequest {
                    seq_id: "synthetic_plus".to_string(),
                    gene_label: Some("SYN1".to_string()),
                    transcript_ids: vec!["SYN1-201".to_string()],
                    cutrun_dataset_ids: vec!["missing_rostock_lane".to_string()],
                    cutrun_catalog_path: Some("/private/study/cutrun.json".to_string()),
                    cutrun_cache_dir: Some("/private/study/cache".to_string()),
                    include_local_source_paths: false,
                    ..Default::default()
                },
            )
            .expect("typed unavailable CUT&RUN evidence");

        assert_eq!(
            report.cutrun_evidence.evaluation_state,
            "not_prepared_or_unevaluated"
        );
        assert_eq!(
            report.cutrun_evidence.lanes[0].state,
            PromoterReporterCutRunLaneState::Unevaluated
        );
        assert!(report.cutrun_evidence.lanes[0].local_source_path.is_none());
        let json = serde_json::to_string(&report).expect("serialize portable report");
        assert!(!json.contains("/private/study"));
        assert!(json.contains("occupancy/enrichment evidence"));
    }
}
