//! Shared deterministic analyses for multi-promoter reporter panels.
//!
//! This module composes existing vector-site ranking and restriction-site
//! scanning. It does not introduce a second enzyme catalog or cloning model.

use super::*;
use crate::gibson_planning::{
    GIBSON_ASSEMBLY_PLAN_SCHEMA, GibsonAssemblyPlan, GibsonPlanAssemblyMember,
    GibsonPlanDestination, GibsonPlanEndStrategy, GibsonPlanFragment, GibsonPlanJunction,
    GibsonPlanOpening, GibsonPlanOverlapPartition, GibsonPlanProduct, GibsonPlanValidationPolicy,
};
use std::fs;

#[derive(Debug)]
struct PromoterReporterPanelSimulation {
    detached: DetachedEngineExecution,
    members: Vec<PromoterReporterPanelMemberProposal>,
    cloning_strategy: PromoterReporterPanelCloningStrategyReport,
    workflow: Vec<PromoterReporterPanelWorkflowStep>,
    products: Vec<PromoterReporterPanelProductProposal>,
    created_seq_ids: Vec<String>,
}

#[derive(Debug)]
struct PromoterReporterPanelResolvedMember {
    request: PromoterReporterPanelMemberRequest,
    candidate_set_binding: PromoterReporterPanelSourceBinding,
    candidate: PromoterReporterFragmentCandidate,
    source_seq_id: String,
    source_sequence_sha256: String,
    member_id: String,
    label: String,
    motif_start_in_fragment_0based: usize,
    motif_end_in_fragment_0based_exclusive: usize,
    motif_forward_strand: bool,
    wild_type_seq_id: String,
    mutant_seq_id: String,
    extended_boundary_audit: Option<PromoterReporterPanelExtendedBoundaryAudit>,
}

impl GentleEngine {
    fn promoter_reporter_cds_entry_phase_offset(
        transcript: &gb_io::seq::Feature,
        entry_cds_features: &[&gb_io::seq::Feature],
        transcript_id: &str,
    ) -> Result<usize, EngineError> {
        let mut offsets = BTreeSet::new();
        let mut evidence = Vec::new();
        for feature in std::iter::once(transcript).chain(entry_cds_features.iter().copied()) {
            for (key, value) in &feature.qualifiers {
                let key = key.as_ref();
                if !matches!(key, "codon_start" | "phase") {
                    continue;
                }
                let raw = value.as_deref().map(str::trim).ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Extended transcript '{}' has /{} without a value at the CDS entry",
                        transcript_id, key
                    ),
                    cause_chain: vec![],
                })?;
                let parsed = raw.parse::<usize>().map_err(|_| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Extended transcript '{}' has invalid /{}='{}' at the CDS entry",
                        transcript_id, key, raw
                    ),
                    cause_chain: vec![],
                })?;
                let offset = match key {
                    "codon_start" if (1..=3).contains(&parsed) => parsed - 1,
                    "phase" if parsed <= 2 => parsed,
                    _ => {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Extended transcript '{}' has invalid /{}='{}' at the CDS entry",
                                transcript_id, key, raw
                            ),
                            cause_chain: vec![],
                        });
                    }
                };
                offsets.insert(offset);
                evidence.push(format!("{key}={raw}"));
            }
        }
        if offsets.len() > 1 {
            evidence.sort();
            evidence.dedup();
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Extended transcript '{}' has conflicting CDS-entry phase annotations: {}",
                    transcript_id,
                    evidence.join(", ")
                ),
                cause_chain: vec![],
            });
        }
        Ok(offsets.into_iter().next().unwrap_or_default())
    }

    fn promoter_reporter_location_has_partial_cds_entry(
        location: &gb_io::seq::Location,
        is_reverse: bool,
        cds_entry: usize,
    ) -> bool {
        use gb_io::seq::Location;

        match location {
            Location::Range((start, before), (end, after)) => {
                (!is_reverse && *start == cds_entry as i64 && before.0)
                    || (is_reverse && *end == cds_entry as i64 && after.0)
            }
            Location::Complement(inner) => {
                Self::promoter_reporter_location_has_partial_cds_entry(inner, is_reverse, cds_entry)
            }
            Location::Join(parts)
            | Location::Order(parts)
            | Location::Bond(parts)
            | Location::OneOf(parts) => parts.iter().any(|part| {
                Self::promoter_reporter_location_has_partial_cds_entry(part, is_reverse, cds_entry)
            }),
            Location::External(_, Some(inner)) => {
                Self::promoter_reporter_location_has_partial_cds_entry(inner, is_reverse, cds_entry)
            }
            Location::Between(_, _) | Location::External(_, None) | Location::Gap(_) => false,
        }
    }

    fn promoter_reporter_panel_extended_geometry(
        source: &DNAsequence,
        candidate: &PromoterReporterFragmentCandidate,
        policy: &PromoterReporterPanelExtendedBoundaryPolicy,
    ) -> Result<(usize, usize, PromoterReporterPanelExtendedBoundaryAudit), EngineError> {
        let transcript_id = policy.transcript_id.trim();
        if transcript_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Extended promoter fragments require a non-empty transcript_id"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        let transcript_matches = source
            .features()
            .iter()
            .filter(|feature| Self::is_mrna_feature(feature))
            .filter(|feature| {
                feature
                    .qualifier_values("transcript_id")
                    .any(|value| value.trim().eq_ignore_ascii_case(transcript_id))
            })
            .collect::<Vec<_>>();
        if transcript_matches.len() != 1 {
            return Err(EngineError {
                code: if transcript_matches.is_empty() {
                    ErrorCode::NotFound
                } else {
                    ErrorCode::InvalidInput
                },
                message: format!(
                    "Extended promoter fragment transcript '{}' resolved to {} transcript annotations; exactly one is required",
                    transcript_id,
                    transcript_matches.len()
                ),
                cause_chain: vec![],
            });
        }
        let transcript = transcript_matches[0];
        let is_reverse = feature_is_reverse(transcript);
        let candidate_reverse = matches!(
            candidate.strand.trim().to_ascii_lowercase().as_str(),
            "-" | "reverse" | "minus"
        );
        if is_reverse != candidate_reverse {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Extended transcript '{}' strand does not match candidate '{}'",
                    transcript_id, candidate.candidate_id
                ),
                cause_chain: vec![],
            });
        }

        let cds_features = source
            .features()
            .iter()
            .filter(|feature| {
                feature.kind.to_string().eq_ignore_ascii_case("CDS")
                    && feature_is_reverse(feature) == is_reverse
                    && feature
                        .qualifier_values("transcript_id")
                        .any(|value| value.trim().eq_ignore_ascii_case(transcript_id))
            })
            .collect::<Vec<_>>();
        let mut cds_ranges = vec![];
        for feature in &cds_features {
            collect_location_ranges_usize(&feature.location, &mut cds_ranges);
        }
        cds_ranges.retain(|(start, end)| end > start);
        cds_ranges.sort_unstable();
        cds_ranges.dedup();
        if cds_ranges.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Extended transcript '{}' has no unambiguous CDS start annotation",
                    transcript_id
                ),
                cause_chain: vec![],
            });
        }
        if is_reverse {
            cds_ranges.sort_unstable_by(|left, right| {
                right.1.cmp(&left.1).then_with(|| right.0.cmp(&left.0))
            });
        }
        let cds_entry = if is_reverse {
            cds_ranges[0].1
        } else {
            cds_ranges[0].0
        };
        let entry_cds_features = cds_features
            .iter()
            .copied()
            .filter(|feature| {
                let mut ranges = vec![];
                collect_location_ranges_usize(&feature.location, &mut ranges);
                ranges.iter().any(|(start, end)| {
                    if is_reverse {
                        *end == cds_entry
                    } else {
                        *start == cds_entry
                    }
                })
            })
            .collect::<Vec<_>>();
        if entry_cds_features.iter().any(|feature| {
            Self::promoter_reporter_location_has_partial_cds_entry(
                &feature.location,
                is_reverse,
                cds_entry,
            )
        }) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Extended transcript '{}' has a 5'-partial CDS location and cannot define a canonical start codon",
                    transcript_id
                ),
                cause_chain: vec![],
            });
        }
        let cds_entry_phase_offset = Self::promoter_reporter_cds_entry_phase_offset(
            transcript,
            &entry_cds_features,
            transcript_id,
        )?;
        if cds_entry_phase_offset != 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Extended transcript '{}' has CDS-entry phase offset {}; this is a 5'-partial CDS and cannot define a canonical start codon",
                    transcript_id, cds_entry_phase_offset
                ),
                cause_chain: vec![],
            });
        }
        let source_text = source.get_forward_string();
        let mut remaining_codon_bases = 3usize;
        let mut codon_ranges = Vec::new();
        let mut cds_start_codon = String::with_capacity(3);
        for &(start, end) in &cds_ranges {
            if end > source_text.len() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Extended transcript '{}' CDS start codon is outside the loaded sequence",
                        transcript_id
                    ),
                    cause_chain: vec![],
                });
            }
            let take = remaining_codon_bases.min(end.saturating_sub(start));
            if take == 0 {
                continue;
            }
            let range = if is_reverse {
                (end - take, end)
            } else {
                (start, start + take)
            };
            let segment = &source_text[range.0..range.1];
            if is_reverse {
                cds_start_codon.push_str(&Self::reverse_complement(segment));
            } else {
                cds_start_codon.push_str(&segment.to_ascii_uppercase());
            }
            codon_ranges.push(range);
            remaining_codon_bases -= take;
            if remaining_codon_bases == 0 {
                break;
            }
        }
        if remaining_codon_bases != 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Extended transcript '{}' has fewer than three annotated coding bases at its CDS start",
                    transcript_id
                ),
                cause_chain: vec![],
            });
        }
        if cds_start_codon != "ATG" {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Extended transcript '{}' CDS entry resolves to '{}', not the canonical ATG start codon",
                    transcript_id, cds_start_codon
                ),
                cause_chain: vec![],
            });
        }
        let codon_start = codon_ranges
            .iter()
            .map(|(start, _)| *start)
            .min()
            .unwrap_or_default();
        let codon_end = codon_ranges
            .iter()
            .map(|(_, end)| *end)
            .max()
            .unwrap_or_default();

        let mut exon_ranges = vec![];
        for feature in source.features().iter().filter(|feature| {
            feature.kind.to_string().eq_ignore_ascii_case("exon")
                && feature_is_reverse(feature) == is_reverse
                && feature
                    .qualifier_values("transcript_id")
                    .any(|value| value.trim().eq_ignore_ascii_case(transcript_id))
        }) {
            collect_location_ranges_usize(&feature.location, &mut exon_ranges);
        }
        if exon_ranges.is_empty() {
            collect_location_ranges_usize(&transcript.location, &mut exon_ranges);
        }
        exon_ranges.retain(|(start, end)| end > start);
        exon_ranges.sort_unstable();
        exon_ranges.dedup();
        let mut utr_ranges = exon_ranges
            .iter()
            .filter_map(|(start, end)| {
                if is_reverse {
                    let utr_start = (*start).max(cds_entry);
                    (*end > utr_start).then_some((utr_start, *end))
                } else {
                    let utr_end = (*end).min(cds_entry);
                    (utr_end > *start).then_some((*start, utr_end))
                }
            })
            .collect::<Vec<_>>();
        if is_reverse {
            utr_ranges.sort_unstable_by(|left, right| right.0.cmp(&left.0));
        } else {
            utr_ranges.sort_unstable();
        }
        let mut utr_sequence = String::new();
        for (start, end) in &utr_ranges {
            let segment = &source_text[*start..*end];
            if is_reverse {
                utr_sequence.push_str(&Self::reverse_complement(segment));
            } else {
                utr_sequence.push_str(segment);
            }
        }
        let utr_upper = utr_sequence.to_ascii_uppercase();
        let upstream_atgs = (0..utr_upper.len().saturating_sub(2))
            .filter(|offset| &utr_upper[*offset..*offset + 3] == "ATG")
            .collect::<Vec<_>>();
        let upstream_orfs = upstream_atgs
            .iter()
            .copied()
            .filter(|start| {
                ((*start + 3)..utr_upper.len().saturating_sub(2))
                    .step_by(3)
                    .any(|offset| matches!(&utr_upper[offset..offset + 3], "TAA" | "TAG" | "TGA"))
            })
            .collect::<Vec<_>>();
        let mut warnings = vec![];
        if !upstream_atgs.is_empty() {
            warnings.push(PromoterReporterPanelExtendedWarning {
                kind: PromoterReporterPanelExtendedWarningKind::UpstreamAtg,
                transcript_id: transcript_id.to_string(),
                detail: format!(
                    "The spliced 5' UTR contains {} upstream ATG start candidate(s)",
                    upstream_atgs.len()
                ),
                transcript_positions_0based: upstream_atgs,
            });
        }
        if !upstream_orfs.is_empty() {
            warnings.push(PromoterReporterPanelExtendedWarning {
                kind: PromoterReporterPanelExtendedWarningKind::UpstreamOrf,
                transcript_id: transcript_id.to_string(),
                detail: format!(
                    "The spliced 5' UTR contains {} upstream ATG(s) with an in-frame stop before the annotated CDS",
                    upstream_orfs.len()
                ),
                transcript_positions_0based: upstream_orfs,
            });
        }
        let mut genomic_utr_ranges = utr_ranges.clone();
        genomic_utr_ranges.sort_unstable();
        let five_prime_utr_intron_count = genomic_utr_ranges
            .windows(2)
            .filter(|pair| pair[1].0 > pair[0].1)
            .count();
        if five_prime_utr_intron_count > 0 {
            warnings.push(PromoterReporterPanelExtendedWarning {
                kind: PromoterReporterPanelExtendedWarningKind::FivePrimeUtrIntron,
                transcript_id: transcript_id.to_string(),
                detail: format!(
                    "The genomic extended fragment contains {} intron(s) between 5' UTR exons",
                    five_prime_utr_intron_count
                ),
                transcript_positions_0based: vec![],
            });
        }
        let (fragment_start, fragment_end) = if is_reverse {
            (codon_start, candidate.end_0based_exclusive)
        } else {
            (candidate.start_0based, codon_end)
        };
        if fragment_start >= fragment_end
            || candidate.variant_start_0based < fragment_start
            || candidate.variant_end_0based_exclusive > fragment_end
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Extended transcript '{}' CDS boundary does not contain candidate '{}' motif anchor",
                    transcript_id, candidate.candidate_id
                ),
                cause_chain: vec![],
            });
        }
        Ok((
            fragment_start,
            fragment_end,
            PromoterReporterPanelExtendedBoundaryAudit {
                policy: policy.clone(),
                strand: if is_reverse { "-" } else { "+" }.to_string(),
                cds_start_codon_source_ranges_0based: codon_ranges,
                cds_start_codon_source_start_0based: codon_start,
                cds_start_codon_source_end_0based_exclusive: codon_end,
                cds_start_codon_5prime_to_3prime: cds_start_codon,
                five_prime_utr_spliced_length_bp: utr_upper.len(),
                five_prime_utr_exon_ranges_0based: utr_ranges,
                warnings,
            },
        ))
    }

    /// Select one directional restriction pair that is usable across every
    /// panel insert, or return an explicit Gibson fallback.
    pub fn restriction_cloning_panel_strategy(
        &self,
        destination_vector_seq_id: &str,
        insert_seq_ids: &[String],
    ) -> Result<PromoterReporterPanelCloningStrategyReport, EngineError> {
        let destination_vector_seq_id = destination_vector_seq_id.trim();
        if destination_vector_seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel cloning strategy requires a destination vector seq_id".to_string(),
                cause_chain: vec![],
            });
        }
        if insert_seq_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel cloning strategy requires at least one insert seq_id".to_string(),
                cause_chain: vec![],
            });
        }

        let mut normalized_insert_ids = Vec::with_capacity(insert_seq_ids.len());
        let mut seen = HashSet::new();
        for raw in insert_seq_ids {
            let insert_seq_id = raw.trim();
            if insert_seq_id.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "Panel cloning strategy insert seq_ids must be non-empty".to_string(),
                    cause_chain: vec![],
                });
            }
            if insert_seq_id == destination_vector_seq_id {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Destination vector '{}' cannot also be a panel insert",
                        destination_vector_seq_id
                    ),
                    cause_chain: vec![],
                });
            }
            if !seen.insert(insert_seq_id.to_string()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Duplicate panel insert seq_id '{}'", insert_seq_id),
                    cause_chain: vec![],
                });
            }
            if !self.state.sequences.contains_key(insert_seq_id) {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Panel insert sequence '{}' not found", insert_seq_id),
                    cause_chain: vec![],
                });
            }
            normalized_insert_ids.push(insert_seq_id.to_string());
        }

        let vector_suggestions =
            self.restriction_cloning_vector_enzyme_suggestions(destination_vector_seq_id)?;
        let candidate_enzymes = vector_suggestions
            .recommended_directed_pairs
            .iter()
            .flat_map(|pair| [&pair.forward_enzyme, &pair.reverse_enzyme])
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();

        let mut insert_site_summaries = Vec::with_capacity(normalized_insert_ids.len());
        for insert_seq_id in &normalized_insert_ids {
            let insert_length_bp = self
                .state
                .sequences
                .get(insert_seq_id)
                .map(DNAsequence::len)
                .unwrap_or_default();
            let mut physical_sites = BTreeSet::new();
            if !candidate_enzymes.is_empty() {
                let scan = self.find_restriction_sites(
                    SequenceScanTarget::SeqId {
                        seq_id: insert_seq_id.clone(),
                        span_start_0based: None,
                        span_end_0based_exclusive: None,
                    },
                    &candidate_enzymes,
                    None,
                    false,
                    None,
                    None,
                )?;
                for row in scan.rows {
                    physical_sites.insert((
                        row.enzyme_name,
                        row.recognition_start_0based,
                        row.recognition_end_0based_exclusive,
                    ));
                }
            }
            let mut site_count_by_enzyme = candidate_enzymes
                .iter()
                .cloned()
                .map(|enzyme| (enzyme, 0usize))
                .collect::<BTreeMap<_, _>>();
            for (enzyme, _, _) in physical_sites {
                *site_count_by_enzyme.entry(enzyme).or_default() += 1;
            }
            insert_site_summaries.push(PromoterReporterPanelInsertSiteSummary {
                insert_seq_id: insert_seq_id.clone(),
                insert_length_bp,
                site_count_by_enzyme,
            });
        }

        let mut pair_evaluations =
            Vec::with_capacity(vector_suggestions.recommended_directed_pairs.len());
        let mut selected_pair = None;
        for (pair_index, pair) in vector_suggestions
            .recommended_directed_pairs
            .iter()
            .enumerate()
        {
            let mut blockers = vec![];
            for insert in &insert_site_summaries {
                for enzyme in [&pair.forward_enzyme, &pair.reverse_enzyme] {
                    let site_count = insert
                        .site_count_by_enzyme
                        .get(enzyme)
                        .copied()
                        .unwrap_or_default();
                    if site_count > 0 {
                        blockers.push(PromoterReporterPanelCloningBlocker {
                            insert_seq_id: insert.insert_seq_id.clone(),
                            enzyme: enzyme.clone(),
                            site_count,
                        });
                    }
                }
            }
            let compatible = blockers.is_empty();
            if compatible && selected_pair.is_none() {
                selected_pair = Some(pair.clone());
            }
            pair_evaluations.push(PromoterReporterPanelPairEvaluation {
                rank_1based: pair_index + 1,
                pair: pair.clone(),
                compatible,
                blockers,
            });
        }

        let (strategy, fallback_reason) = if selected_pair.is_some() {
            (
                PromoterReporterPanelCloningStrategy::DirectionalRestriction,
                None,
            )
        } else if vector_suggestions.recommended_directed_pairs.is_empty() {
            (
                PromoterReporterPanelCloningStrategy::Gibson,
                Some(format!(
                    "Destination vector '{}' does not expose two ordered unique restriction sites; use Gibson assembly",
                    destination_vector_seq_id
                )),
            )
        } else {
            (
                PromoterReporterPanelCloningStrategy::Gibson,
                Some(
                    "Every vector-ranked directional restriction pair is blocked by at least one internal site across the panel; use Gibson assembly"
                        .to_string(),
                ),
            )
        };

        Ok(PromoterReporterPanelCloningStrategyReport {
            schema: PROMOTER_REPORTER_PANEL_CLONING_STRATEGY_SCHEMA.to_string(),
            destination_vector_seq_id: destination_vector_seq_id.to_string(),
            insert_seq_ids: normalized_insert_ids,
            vector_suggestions,
            insert_site_summaries,
            pair_evaluations,
            strategy,
            selected_pair,
            fallback_reason,
            warnings: vec![
                "Restriction compatibility is a sequence-derived planning screen; verify enzyme availability and reaction conditions before use."
                    .to_string(),
            ],
        })
    }

    /// Apply the fixed p53-family core-disruption rule and audit its local PWM
    /// and selected restriction-site consequences.
    pub fn design_p53_family_motif_disruption(
        &self,
        sequence_text: &str,
        motif_start_0based: usize,
        motif_end_0based_exclusive: usize,
        motif_forward_strand: bool,
        selected_enzymes: &[String],
    ) -> Result<P53FamilyMotifDisruptionReport, EngineError> {
        const STATED_RULE: &str = "For each oriented 10-bp p53-family half-site, replace position 4 C with A and position 7 G with T (1-based within the half-site)";
        const P53_FAMILY_MOTIFS: [&str; 3] = ["MA0861.1", "MA0106.3", "MA0525.2"];

        let normalized_sequence = sequence_text
            .chars()
            .filter(|base| !base.is_ascii_whitespace())
            .map(|base| base.to_ascii_uppercase())
            .collect::<String>();
        if normalized_sequence.is_empty()
            || normalized_sequence
                .bytes()
                .any(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "p53-family motif disruption requires a non-empty A/C/G/T sequence"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        if motif_start_0based >= motif_end_0based_exclusive
            || motif_end_0based_exclusive > normalized_sequence.len()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "p53-family motif interval {}..{} is outside a {} bp sequence",
                    motif_start_0based,
                    motif_end_0based_exclusive,
                    normalized_sequence.len()
                ),
                cause_chain: vec![],
            });
        }
        let motif_length = motif_end_0based_exclusive - motif_start_0based;
        if !matches!(motif_length, 10 | 20) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "p53-family motif disruption v1 requires one or two contiguous 10-bp half-sites, found {} bp",
                    motif_length
                ),
                cause_chain: vec![],
            });
        }

        let original_motif_source =
            normalized_sequence[motif_start_0based..motif_end_0based_exclusive].to_string();
        let original_motif_oriented = if motif_forward_strand {
            original_motif_source.clone()
        } else {
            Self::reverse_complement(&original_motif_source)
        };
        let oriented_bytes = original_motif_oriented.as_bytes();
        for half_site_index in 0..(motif_length / 10) {
            let offset = half_site_index * 10;
            if oriented_bytes[offset + 3] != b'C' || oriented_bytes[offset + 6] != b'G' {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "p53-family half-site {} does not match the stated-rule core: expected C at oriented position {} and G at position {}",
                        half_site_index + 1,
                        offset + 4,
                        offset + 7
                    ),
                    cause_chain: vec![],
                });
            }
        }

        let mut mutant_bytes = normalized_sequence.as_bytes().to_vec();
        let mut changes = vec![];
        for half_site_index in 0..(motif_length / 10) {
            for (within_half_0based, oriented_replacement) in [(3usize, b'A'), (6usize, b'T')] {
                let oriented_position_0based = half_site_index * 10 + within_half_0based;
                let source_position_0based = if motif_forward_strand {
                    motif_start_0based + oriented_position_0based
                } else {
                    motif_end_0based_exclusive - 1 - oriented_position_0based
                };
                let source_replacement = if motif_forward_strand {
                    oriented_replacement
                } else {
                    IupacCode::letter_complement(oriented_replacement).to_ascii_uppercase()
                };
                let original_base = mutant_bytes[source_position_0based];
                mutant_bytes[source_position_0based] = source_replacement;
                changes.push(P53FamilyMotifBaseChange {
                    source_position_0based,
                    oriented_position_0based,
                    half_site_index_1based: half_site_index + 1,
                    original_base: original_base as char,
                    replacement_base: source_replacement as char,
                    rule: if within_half_0based == 3 {
                        "oriented_core_C_to_A".to_string()
                    } else {
                        "oriented_core_G_to_T".to_string()
                    },
                });
            }
        }
        changes.sort_by_key(|change| change.source_position_0based);
        let mutant_sequence = String::from_utf8(mutant_bytes).map_err(|error| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not encode p53-family motif mutant: {error}"),
            cause_chain: vec![],
        })?;
        let mutant_motif_source =
            mutant_sequence[motif_start_0based..motif_end_0based_exclusive].to_string();
        let mutant_motif_oriented = if motif_forward_strand {
            mutant_motif_source.clone()
        } else {
            Self::reverse_complement(&mutant_motif_source)
        };
        let changed_positions = changes
            .iter()
            .map(|change| change.source_position_0based)
            .collect::<BTreeSet<_>>();

        let motif_tokens = P53_FAMILY_MOTIFS
            .iter()
            .map(|motif| motif.to_string())
            .collect::<Vec<_>>();
        let wild_type_pwm = self.scan_tfbs_hits(
            SequenceScanTarget::InlineSequence {
                sequence_text: normalized_sequence.clone(),
                topology: InlineSequenceTopology::Linear,
                id_hint: Some("p53_family_wild_type".to_string()),
                span_start_0based: None,
                span_end_0based_exclusive: None,
            },
            &motif_tokens,
            None,
            Some(0.0),
            &[],
            None,
            None,
            None,
        )?;
        let mutant_pwm = self.scan_tfbs_hits(
            SequenceScanTarget::InlineSequence {
                sequence_text: mutant_sequence.clone(),
                topology: InlineSequenceTopology::Linear,
                id_hint: Some("p53_family_mutant".to_string()),
                span_start_0based: None,
                span_end_0based_exclusive: None,
            },
            &motif_tokens,
            None,
            Some(0.0),
            &[],
            None,
            None,
            None,
        )?;
        let mut pwm_audits = vec![];
        for tf_id in P53_FAMILY_MOTIFS {
            let wild_type_scores = Self::p53_family_affected_pwm_scores(
                &wild_type_pwm.rows,
                tf_id,
                &changed_positions,
            );
            let mutant_scores =
                Self::p53_family_affected_pwm_scores(&mutant_pwm.rows, tf_id, &changed_positions);
            let keys = wild_type_scores
                .keys()
                .chain(mutant_scores.keys())
                .copied()
                .collect::<BTreeSet<_>>();
            let stronger_mutant_window_count = keys
                .iter()
                .filter(|key| {
                    mutant_scores.get(key).copied().unwrap_or(f64::NEG_INFINITY)
                        > wild_type_scores
                            .get(key)
                            .copied()
                            .unwrap_or(f64::NEG_INFINITY)
                            + 1e-9
                })
                .count();
            let wild_type_max_llr_bits = Self::finite_max(wild_type_scores.values().copied());
            let mutant_max_llr_bits = Self::finite_max(mutant_scores.values().copied());
            let max_llr_delta_bits = wild_type_max_llr_bits
                .zip(mutant_max_llr_bits)
                .map(|(wild_type, mutant)| mutant - wild_type);
            let tf_name = wild_type_pwm
                .rows
                .iter()
                .chain(mutant_pwm.rows.iter())
                .find(|row| row.tf_id == tf_id)
                .and_then(|row| row.tf_name.clone());
            pwm_audits.push(P53FamilyMotifPwmAudit {
                tf_id: tf_id.to_string(),
                tf_name,
                windows_compared: keys.len(),
                wild_type_max_llr_bits,
                mutant_max_llr_bits,
                max_llr_delta_bits,
                stronger_mutant_window_count,
            });
        }
        // Shifted windows can become less unfavorable even when the strongest
        // local p53-family match is weakened. Treat only an increased local
        // maximum as a stronger competing motif; retain per-window increases
        // in the audit rows as diagnostics.
        let no_stronger_p53_family_hit_near_edit = pwm_audits.iter().all(|audit| {
            audit
                .max_llr_delta_bits
                .is_none_or(|delta_bits| delta_bits <= 1e-9)
        });

        let (restriction_site_audits, creates_selected_restriction_site) = self
            .p53_family_restriction_site_audits(
                &normalized_sequence,
                &mutant_sequence,
                selected_enzymes,
            )?;

        Ok(P53FamilyMotifDisruptionReport {
            schema: P53_FAMILY_MOTIF_DISRUPTION_SCHEMA.to_string(),
            stated_rule: STATED_RULE.to_string(),
            source_sequence_sha256: sha256_prefixed_str(&normalized_sequence),
            motif_start_0based,
            motif_end_0based_exclusive,
            motif_forward_strand,
            original_motif_source,
            original_motif_oriented,
            mutant_motif_source,
            mutant_motif_oriented,
            length_preserved: mutant_sequence.len() == normalized_sequence.len(),
            mutant_sequence,
            changes,
            pwm_audits,
            no_stronger_p53_family_hit_near_edit,
            restriction_site_audits,
            creates_selected_restriction_site,
            warnings: vec![
                "This is a sequence-only disruption proposal, not evidence of lost binding, occupancy, transcriptional effect, or reporter activity."
                    .to_string(),
                "PWM comparisons are post-edit audits of the stated rule; they do not choose or optimize the substitutions."
                    .to_string(),
            ],
        })
    }

    fn p53_family_affected_pwm_scores(
        rows: &[TfbsHitScanRow],
        tf_id: &str,
        changed_positions: &BTreeSet<usize>,
    ) -> BTreeMap<(usize, bool), f64> {
        rows.iter()
            .filter(|row| row.tf_id == tf_id && !row.wraps_origin)
            .filter(|row| {
                changed_positions.iter().any(|position| {
                    row.match_start_0based <= *position
                        && *position < row.match_end_0based_exclusive
                })
            })
            .map(|row| ((row.match_start_0based, row.forward_strand), row.llr_bits))
            .collect()
    }

    fn finite_max(values: impl Iterator<Item = f64>) -> Option<f64> {
        values
            .filter(|value| value.is_finite())
            .reduce(|left, right| left.max(right))
    }

    fn p53_family_restriction_site_audits(
        &self,
        wild_type_sequence: &str,
        mutant_sequence: &str,
        selected_enzymes: &[String],
    ) -> Result<(Vec<P53FamilyRestrictionSiteAudit>, bool), EngineError> {
        let mut enzymes = vec![];
        let mut seen = HashSet::new();
        for raw in selected_enzymes {
            let trimmed = raw.trim();
            if !trimmed.is_empty() && seen.insert(trimmed.to_ascii_uppercase()) {
                enzymes.push(trimmed.to_string());
            }
        }
        if enzymes.is_empty() {
            return Ok((vec![], false));
        }
        let scan = |sequence_text: &str, id_hint: &str| {
            self.find_restriction_sites(
                SequenceScanTarget::InlineSequence {
                    sequence_text: sequence_text.to_string(),
                    topology: InlineSequenceTopology::Linear,
                    id_hint: Some(id_hint.to_string()),
                    span_start_0based: None,
                    span_end_0based_exclusive: None,
                },
                &enzymes,
                None,
                false,
                None,
                None,
            )
        };
        let wild_type = scan(wild_type_sequence, "p53_family_wild_type")?;
        let mutant = scan(mutant_sequence, "p53_family_mutant")?;
        let to_sites = |report: RestrictionSiteScanReport| {
            let mut by_enzyme: BTreeMap<String, BTreeSet<P53FamilyRestrictionSiteInterval>> =
                report
                    .enzymes_scanned
                    .iter()
                    .cloned()
                    .map(|enzyme| (enzyme, BTreeSet::new()))
                    .collect();
            for row in report.rows {
                by_enzyme.entry(row.enzyme_name).or_default().insert(
                    P53FamilyRestrictionSiteInterval {
                        start_0based: row.recognition_start_0based,
                        end_0based_exclusive: row.recognition_end_0based_exclusive,
                    },
                );
            }
            by_enzyme
        };
        let wild_type_sites = to_sites(wild_type);
        let mutant_sites = to_sites(mutant);
        let enzyme_names = wild_type_sites
            .keys()
            .chain(mutant_sites.keys())
            .cloned()
            .collect::<BTreeSet<_>>();
        let mut audits = vec![];
        for enzyme in enzyme_names {
            let before = wild_type_sites.get(&enzyme).cloned().unwrap_or_default();
            let after = mutant_sites.get(&enzyme).cloned().unwrap_or_default();
            audits.push(P53FamilyRestrictionSiteAudit {
                enzyme,
                wild_type_site_count: before.len(),
                mutant_site_count: after.len(),
                newly_created_sites: after.difference(&before).cloned().collect(),
                lost_sites: before.difference(&after).cloned().collect(),
            });
        }
        let creates_selected_restriction_site = audits
            .iter()
            .any(|audit| !audit.newly_created_sites.is_empty());
        Ok((audits, creates_selected_restriction_site))
    }
}

impl GentleEngine {
    /// Build a read-only, content-addressed proposal for a complete promoter
    /// reporter panel. All sequence construction happens on a detached clone.
    pub fn plan_promoter_reporter_panel(
        &self,
        request: PromoterReporterPanelRequest,
    ) -> Result<PromoterReporterPanelProposal, EngineError> {
        self.build_promoter_reporter_panel_proposal(request)
            .map(|(proposal, _)| proposal)
    }

    fn build_promoter_reporter_panel_proposal(
        &self,
        request: PromoterReporterPanelRequest,
    ) -> Result<
        (
            PromoterReporterPanelProposal,
            PromoterReporterPanelSimulation,
        ),
        EngineError,
    > {
        let request = Self::normalize_promoter_reporter_panel_request(request)?;
        let baseline_state_sha256 = self.promoter_reporter_panel_value_sha256(
            self.snapshot(),
            "promoter-reporter panel project baseline",
        )?;
        let request_sha256 =
            self.promoter_reporter_panel_value_sha256(&request, "panel request")?;

        let backbone = self.resolve_reporter_backbone(
            &request.vector_seq_id,
            None,
            true,
            Some(&request.vector_catalog_id),
            request.helper_catalog_path.as_deref(),
        )?;
        let vector_validation = backbone.validation.ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: "Panel planning requires catalog-owned exact vector validation".to_string(),
            cause_chain: vec![],
        })?;
        if vector_validation.status != ReporterVectorValidationStatus::Verified {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Panel planning requires a verified exact vector; '{}' is {:?}",
                    request.vector_seq_id, vector_validation.status
                ),
                cause_chain: vec![],
            });
        }
        let vector = self
            .state
            .sequences
            .get(&request.vector_seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Verified panel vector '{}' is not loaded in project state",
                    request.vector_seq_id
                ),
                cause_chain: vec![],
            })?;
        let vector_sequence_sha256 = sha256_prefixed_str(&vector.get_forward_string());

        let (helper_catalog, helper_catalog_origin) =
            Self::open_helper_genome_catalog(request.helper_catalog_path.as_deref())?;
        let helper_catalog_path =
            Self::canonical_existing_path(&helper_catalog_origin, "helper-vector catalog")?;
        let helper_catalog_bytes = fs::read(&helper_catalog_path).map_err(|error| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not read helper-vector catalog '{}': {error}",
                helper_catalog_path.display()
            ),
            cause_chain: vec![],
        })?;
        let helper_catalog_binding = PromoterReporterPanelSourceBinding {
            path: helper_catalog_path.to_string_lossy().to_string(),
            sha256: sha256_prefixed_bytes(&helper_catalog_bytes),
        };
        let expectation = helper_catalog
            .helper_vector_sequence_expectation(&request.vector_catalog_id)
            .map_err(|message| EngineError {
                code: ErrorCode::InvalidInput,
                message,
                cause_chain: vec![],
            })?
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Helper-vector '{}' has no exact sequence expectation",
                    request.vector_catalog_id
                ),
                cause_chain: vec![],
            })?;
        let mcs_expectation = expectation
            .required_features
            .iter()
            .find(|feature| feature.id.eq_ignore_ascii_case("multiple_cloning_region"))
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Helper-vector '{}' does not declare a multiple_cloning_region feature",
                    request.vector_catalog_id
                ),
                cause_chain: vec![],
            })?;
        let (mcs_start_1based, mcs_end_1based, _) =
            Self::find_reporter_vector_feature(vector, mcs_expectation).ok_or_else(|| {
                EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Verified vector '{}' has no resolved multiple-cloning-region annotation",
                        request.vector_seq_id
                    ),
                    cause_chain: vec![],
                }
            })?;
        let mcs_start_0based = mcs_start_1based.saturating_sub(1);
        let mcs_end_0based_exclusive = mcs_end_1based;

        let resolved_members = self.resolve_promoter_reporter_panel_members(&request)?;
        let mut simulation = self.simulate_promoter_reporter_panel(
            &request,
            &resolved_members,
            mcs_start_0based,
            mcs_end_0based_exclusive,
        )?;
        let artifacts = Self::promoter_reporter_panel_artifacts(&request, &simulation.products);
        Self::append_promoter_reporter_panel_export_steps(&mut simulation.workflow, &artifacts)?;
        let workflow_sha256 =
            self.promoter_reporter_panel_value_sha256(&simulation.workflow, "panel workflow")?;
        let mut nonclaims = vec![
            "Sequence motifs are candidate regulatory evidence, not proof of p53-family occupancy or promoter function."
                .to_string(),
            "In-silico construct and primer models require wet-lab validation; GENtle does not order oligos or execute experiments."
                .to_string(),
        ];
        nonclaims.extend(request.scientific_caveats.iter().cloned());
        let mut proposal = PromoterReporterPanelProposal {
            schema: PROMOTER_REPORTER_PANEL_PROPOSAL_SCHEMA.to_string(),
            request,
            request_sha256,
            baseline_state_sha256,
            helper_catalog: helper_catalog_binding,
            vector_sequence_sha256,
            vector_validation,
            mcs_start_0based,
            mcs_end_0based_exclusive,
            members: simulation.members.clone(),
            cloning_strategy: simulation.cloning_strategy.clone(),
            workflow: simulation.workflow.clone(),
            workflow_sha256,
            products: simulation.products.clone(),
            artifacts,
            approval_required: true,
            warnings: simulation
                .cloning_strategy
                .warnings
                .iter()
                .cloned()
                .chain(
                    simulation
                        .members
                        .iter()
                        .flat_map(|member| member.warnings.clone()),
                )
                .chain(
                    simulation
                        .products
                        .iter()
                        .flat_map(|product| product.final_product_audit.warnings.clone()),
                )
                .collect(),
            nonclaims,
            ..PromoterReporterPanelProposal::default()
        };
        proposal.proposal_digest = Self::promoter_reporter_panel_proposal_digest(&proposal)?;
        proposal.proposal_id = format!(
            "promoter_panel_{}",
            &proposal
                .proposal_digest
                .strip_prefix("sha256:")
                .unwrap_or(&proposal.proposal_digest)[..16]
        );
        Ok((proposal, simulation))
    }

    /// Materialize only an unchanged, digest-approved panel proposal. The
    /// exact proposal is recomputed before files or project state are changed.
    pub fn materialize_promoter_reporter_panel(
        &mut self,
        proposal: PromoterReporterPanelProposal,
        approval_digest: &str,
    ) -> Result<PromoterReporterPanelReceipt, EngineError> {
        if proposal.schema != PROMOTER_REPORTER_PANEL_PROPOSAL_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported promoter-reporter panel proposal schema '{}'",
                    proposal.schema
                ),
                cause_chain: vec![],
            });
        }
        let embedded_digest = Self::promoter_reporter_panel_proposal_digest(&proposal)?;
        if proposal.proposal_digest != embedded_digest {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel proposal content does not match its embedded proposal_digest"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        if approval_digest.trim() != proposal.proposal_digest {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Panel materialization requires exact approval digest '{}'",
                    proposal.proposal_digest
                ),
                cause_chain: vec![],
            });
        }

        let (fresh_proposal, mut simulation) =
            self.build_promoter_reporter_panel_proposal(proposal.request.clone())?;
        if fresh_proposal.proposal_digest != proposal.proposal_digest {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel proposal is stale because project state or a bound input changed; plan again and approve the new digest"
                    .to_string(),
                cause_chain: vec![],
            });
        }

        let mut artifact_paths = BTreeSet::new();
        for artifact in &fresh_proposal.artifacts {
            if !artifact_paths.insert(artifact.path.clone()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Panel proposal repeats artifact path '{}'", artifact.path),
                    cause_chain: vec![],
                });
            }
            if Path::new(&artifact.path).exists() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Panel materialization refuses to overwrite existing artifact '{}'",
                        artifact.path
                    ),
                    cause_chain: vec![],
                });
            }
        }
        fs::create_dir_all(&fresh_proposal.request.output_dir).map_err(|error| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create panel output directory '{}': {error}",
                fresh_proposal.request.output_dir
            ),
            cause_chain: vec![],
        })?;

        let export_result = Self::write_promoter_reporter_panel_artifacts(
            simulation.detached.engine_mut(),
            &fresh_proposal,
        );
        if let Err(error) = export_result {
            Self::remove_promoter_reporter_panel_artifacts(&fresh_proposal.artifacts);
            return Err(error);
        }

        let created_seq_ids = simulation.created_seq_ids.clone();
        if let Err(error) = self.commit_detached_execution(&mut simulation.detached) {
            Self::remove_promoter_reporter_panel_artifacts(&fresh_proposal.artifacts);
            return Err(error);
        }
        let completed_state_sha256 = self.promoter_reporter_panel_value_sha256(
            self.snapshot(),
            "materialized promoter-reporter panel state",
        )?;
        let manifest_path = fresh_proposal
            .artifacts
            .iter()
            .find(|artifact| artifact.artifact_kind == "construct_primer_manifest_tsv")
            .map(|artifact| artifact.path.clone())
            .unwrap_or_default();
        Ok(PromoterReporterPanelReceipt {
            schema: PROMOTER_REPORTER_PANEL_RECEIPT_SCHEMA.to_string(),
            proposal_id: fresh_proposal.proposal_id,
            approved_proposal_digest: fresh_proposal.proposal_digest,
            completed_state_sha256,
            created_seq_ids,
            product_seq_ids: fresh_proposal
                .products
                .iter()
                .map(|product| product.product_seq_id.clone())
                .collect(),
            artifact_paths: fresh_proposal
                .artifacts
                .iter()
                .map(|artifact| artifact.path.clone())
                .collect(),
            manifest_path,
            warnings: fresh_proposal.warnings,
        })
    }

    fn promoter_reporter_panel_artifacts(
        request: &PromoterReporterPanelRequest,
        products: &[PromoterReporterPanelProductProposal],
    ) -> Vec<PromoterReporterPanelArtifactProposal> {
        let output_dir = PathBuf::from(&request.output_dir);
        let mut artifacts = vec![];
        for product in products {
            let file_stem = Self::promoter_reporter_panel_slug(&product.product_seq_id);
            artifacts.push(PromoterReporterPanelArtifactProposal {
                artifact_kind: "product_genbank".to_string(),
                product_seq_id: Some(product.product_seq_id.clone()),
                path: output_dir
                    .join(format!("{file_stem}.gb"))
                    .to_string_lossy()
                    .to_string(),
            });
            artifacts.push(PromoterReporterPanelArtifactProposal {
                artifact_kind: "product_svg".to_string(),
                product_seq_id: Some(product.product_seq_id.clone()),
                path: output_dir
                    .join(format!("{file_stem}.svg"))
                    .to_string_lossy()
                    .to_string(),
            });
        }
        artifacts.push(PromoterReporterPanelArtifactProposal {
            artifact_kind: "construct_primer_manifest_tsv".to_string(),
            product_seq_id: None,
            path: output_dir
                .join(format!(
                    "{}.construct_primer_manifest.tsv",
                    request.panel_id
                ))
                .to_string_lossy()
                .to_string(),
        });
        artifacts
    }

    fn append_promoter_reporter_panel_export_steps(
        workflow: &mut Vec<PromoterReporterPanelWorkflowStep>,
        artifacts: &[PromoterReporterPanelArtifactProposal],
    ) -> Result<(), EngineError> {
        for artifact in artifacts {
            let operation = match artifact.artifact_kind.as_str() {
                "product_genbank" => serde_json::to_value(Operation::SaveFile {
                    seq_id: artifact.product_seq_id.clone().unwrap_or_default(),
                    path: artifact.path.clone(),
                    format: ExportFormat::GenBank,
                }),
                "product_svg" => serde_json::to_value(Operation::RenderSequenceSvg {
                    seq_id: artifact.product_seq_id.clone().unwrap_or_default(),
                    mode: RenderSvgMode::Circular,
                    path: artifact.path.clone(),
                }),
                "construct_primer_manifest_tsv" => Ok(json!({
                    "kind": "write_promoter_reporter_panel_construct_primer_manifest",
                    "path": artifact.path,
                })),
                other => {
                    return Err(EngineError {
                        code: ErrorCode::Internal,
                        message: format!("Unsupported planned panel artifact kind '{other}'"),
                        cause_chain: vec![],
                    });
                }
            }
            .map_err(|error| EngineError {
                code: ErrorCode::Internal,
                message: format!("Could not serialize panel export operation: {error}"),
                cause_chain: vec![],
            })?;
            Self::record_promoter_reporter_panel_action(
                workflow,
                "export_approved_panel_artifact",
                None,
                operation,
                vec![],
            )?;
        }
        Ok(())
    }

    fn write_promoter_reporter_panel_artifacts(
        engine: &mut GentleEngine,
        proposal: &PromoterReporterPanelProposal,
    ) -> Result<(), EngineError> {
        for artifact in &proposal.artifacts {
            match artifact.artifact_kind.as_str() {
                "product_genbank" => {
                    engine.apply(Operation::SaveFile {
                        seq_id: artifact.product_seq_id.clone().unwrap_or_default(),
                        path: artifact.path.clone(),
                        format: ExportFormat::GenBank,
                    })?;
                }
                "product_svg" => {
                    engine.apply(Operation::RenderSequenceSvg {
                        seq_id: artifact.product_seq_id.clone().unwrap_or_default(),
                        mode: RenderSvgMode::Circular,
                        path: artifact.path.clone(),
                    })?;
                }
                "construct_primer_manifest_tsv" => {
                    let text = Self::promoter_reporter_panel_manifest_tsv(proposal)?;
                    fs::write(&artifact.path, text).map_err(|error| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not write panel construct/primer manifest '{}': {error}",
                            artifact.path
                        ),
                        cause_chain: vec![],
                    })?;
                }
                other => {
                    return Err(EngineError {
                        code: ErrorCode::Internal,
                        message: format!("Unsupported panel artifact kind '{other}'"),
                        cause_chain: vec![],
                    });
                }
            }
        }
        Ok(())
    }

    fn promoter_reporter_panel_manifest_tsv(
        proposal: &PromoterReporterPanelProposal,
    ) -> Result<String, EngineError> {
        let member_by_id = proposal
            .members
            .iter()
            .map(|member| (member.member_id.as_str(), member))
            .collect::<BTreeMap<_, _>>();
        let mut text = String::from(
            "panel_id\tmember_id\tallele\tfragment_role\tgene_label\tpromoter_class_id\tinsert_seq_id\tproduct_seq_id\tassembly_model\tjunction_validation\tfinal_restriction_site_count\tunexpected_restriction_site_count_excess\tlength_bp\tsequence_sha256\tprimer_seq_ids\tprimer_sequences_5prime_to_3prime\tgenbank_path\tsvg_path\n",
        );
        for product in &proposal.products {
            let member = member_by_id
                .get(product.member_id.as_str())
                .ok_or_else(|| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Panel product '{}' refers to missing member '{}'",
                        product.product_seq_id, product.member_id
                    ),
                    cause_chain: vec![],
                })?;
            let artifact_path = |kind: &str| {
                proposal
                    .artifacts
                    .iter()
                    .find(|artifact| {
                        artifact.artifact_kind == kind
                            && artifact.product_seq_id.as_deref()
                                == Some(product.product_seq_id.as_str())
                    })
                    .map(|artifact| artifact.path.as_str())
                    .unwrap_or_default()
            };
            let cells = [
                proposal.request.panel_id.clone(),
                product.member_id.clone(),
                product.allele.clone(),
                match member.fragment_role {
                    PromoterReporterPanelFragmentRole::Core => "core".to_string(),
                    PromoterReporterPanelFragmentRole::Extended => "extended".to_string(),
                },
                member.gene_label.clone().unwrap_or_default(),
                member.promoter_class_id.clone().unwrap_or_default(),
                product.insert_seq_id.clone(),
                product.product_seq_id.clone(),
                product.assembly_model.clone(),
                product.final_product_audit.junction_validation.clone(),
                product.final_product_audit.matched_site_count.to_string(),
                product
                    .final_product_audit
                    .unexpected_count_excess_by_enzyme
                    .iter()
                    .map(|(enzyme, count)| format!("{enzyme}:+{count}"))
                    .collect::<Vec<_>>()
                    .join(";"),
                product.length_bp.to_string(),
                product.sequence_sha256.clone(),
                product.primer_seq_ids.join(";"),
                product.primer_sequences_5prime_to_3prime.join(";"),
                artifact_path("product_genbank").to_string(),
                artifact_path("product_svg").to_string(),
            ];
            text.push_str(
                &cells
                    .iter()
                    .map(|cell| Self::promoter_reporter_panel_tsv_cell(cell))
                    .collect::<Vec<_>>()
                    .join("\t"),
            );
            text.push('\n');
        }
        Ok(text)
    }

    fn promoter_reporter_panel_tsv_cell(value: &str) -> String {
        value.replace('\t', " ").replace(['\r', '\n'], " ")
    }

    fn remove_promoter_reporter_panel_artifacts(
        artifacts: &[PromoterReporterPanelArtifactProposal],
    ) {
        for artifact in artifacts {
            let _ = fs::remove_file(&artifact.path);
        }
    }

    fn normalize_promoter_reporter_panel_request(
        mut request: PromoterReporterPanelRequest,
    ) -> Result<PromoterReporterPanelRequest, EngineError> {
        if request.schema.trim().is_empty() {
            request.schema = PROMOTER_REPORTER_PANEL_REQUEST_SCHEMA.to_string();
        }
        if request.schema != PROMOTER_REPORTER_PANEL_REQUEST_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported promoter-reporter panel request schema '{}'",
                    request.schema
                ),
                cause_chain: vec![],
            });
        }
        request.panel_id = Self::promoter_reporter_panel_slug(&request.panel_id);
        if request.panel_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel request requires a non-empty panel_id".to_string(),
                cause_chain: vec![],
            });
        }
        request.vector_seq_id = request.vector_seq_id.trim().to_string();
        request.vector_catalog_id = request.vector_catalog_id.trim().to_string();
        if request.vector_seq_id.is_empty() || request.vector_catalog_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel request requires vector_seq_id and vector_catalog_id".to_string(),
                cause_chain: vec![],
            });
        }
        if request.members.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel request requires at least one selected promoter candidate"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        if request.mutation_policy == PromoterReporterPanelMutationPolicy::Unspecified {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel request requires an explicit mutation_policy; promoter-reporter panel v1 currently supports p53_family_core_disruption_v1"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        request.scientific_caveats = request
            .scientific_caveats
            .into_iter()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if let Some(path) = request.helper_catalog_path.as_deref() {
            request.helper_catalog_path = Some(
                Self::canonical_existing_path(path, "helper-vector catalog")?
                    .to_string_lossy()
                    .to_string(),
            );
        }
        for member in &mut request.members {
            member.candidate_set_path = Self::canonical_existing_path(
                &member.candidate_set_path,
                "promoter candidate set",
            )?
            .to_string_lossy()
            .to_string();
            member.candidate_id = member
                .candidate_id
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string);
            member.label = member
                .label
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string);
        }
        let output = PathBuf::from(request.output_dir.trim());
        if output.as_os_str().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Panel request requires output_dir".to_string(),
                cause_chain: vec![],
            });
        }
        request.output_dir = if output.is_absolute() {
            output
        } else {
            env::current_dir()
                .map_err(|error| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not resolve current directory: {error}"),
                    cause_chain: vec![],
                })?
                .join(output)
        }
        .to_string_lossy()
        .to_string();
        Ok(request)
    }

    fn canonical_existing_path(raw: &str, label: &str) -> Result<PathBuf, EngineError> {
        fs::canonicalize(raw.trim()).map_err(|error| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not resolve {label} '{}': {error}", raw.trim()),
            cause_chain: vec![],
        })
    }

    fn promoter_reporter_panel_slug(raw: &str) -> String {
        raw.trim()
            .to_ascii_lowercase()
            .chars()
            .map(|character| {
                if character.is_ascii_alphanumeric() {
                    character
                } else {
                    '_'
                }
            })
            .collect::<String>()
            .split('_')
            .filter(|part| !part.is_empty())
            .collect::<Vec<_>>()
            .join("_")
    }

    fn promoter_reporter_panel_value_sha256<T: Serialize>(
        &self,
        value: &T,
        label: &str,
    ) -> Result<String, EngineError> {
        Self::construct_reasoning_canonical_json(value, label)
            .map(|canonical| sha256_prefixed_str(&canonical))
    }

    fn promoter_reporter_panel_proposal_digest(
        proposal: &PromoterReporterPanelProposal,
    ) -> Result<String, EngineError> {
        let mut basis = proposal.clone();
        basis.proposal_id.clear();
        basis.proposal_digest.clear();
        Self::construct_reasoning_canonical_json(&basis, "promoter-reporter panel proposal")
            .map(|canonical| sha256_prefixed_str(&canonical))
    }

    fn resolve_promoter_reporter_panel_members(
        &self,
        request: &PromoterReporterPanelRequest,
    ) -> Result<Vec<PromoterReporterPanelResolvedMember>, EngineError> {
        let mut resolved = vec![];
        let mut seen_ids = HashSet::new();
        for (index, member_request) in request.members.iter().enumerate() {
            let bytes =
                fs::read(&member_request.candidate_set_path).map_err(|error| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not read promoter candidate set '{}': {error}",
                        member_request.candidate_set_path
                    ),
                    cause_chain: vec![],
                })?;
            let candidate_set =
                Self::read_promoter_reporter_candidate_set(&member_request.candidate_set_path)?;
            if candidate_set.schema != PROMOTER_REPORTER_CANDIDATES_SCHEMA {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate set '{}' uses unsupported schema '{}'",
                        member_request.candidate_set_path, candidate_set.schema
                    ),
                    cause_chain: vec![],
                });
            }
            let candidate_id = member_request
                .candidate_id
                .as_deref()
                .unwrap_or(&candidate_set.recommended_candidate_id);
            let mut candidate = candidate_set
                .candidates
                .iter()
                .find(|candidate| candidate.candidate_id == candidate_id)
                .cloned()
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Candidate '{}' was not found in '{}'",
                        candidate_id, member_request.candidate_set_path
                    ),
                    cause_chain: vec![],
                })?;
            let source = self
                .state
                .sequences
                .get(&candidate_set.seq_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Candidate source sequence '{}' is not loaded",
                        candidate_set.seq_id
                    ),
                    cause_chain: vec![],
                })?;
            if source.len() != candidate_set.sequence_length_bp
                || candidate.end_0based_exclusive > source.len()
                || candidate.start_0based >= candidate.end_0based_exclusive
            {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate '{}' geometry no longer matches loaded sequence '{}'",
                        candidate.candidate_id, candidate_set.seq_id
                    ),
                    cause_chain: vec![],
                });
            }
            let anchor = candidate
                .anchor
                .clone()
                .or_else(|| candidate_set.anchor.clone())
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate '{}' has no resolved motif anchor",
                        candidate.candidate_id
                    ),
                    cause_chain: vec![],
                })?;
            if anchor.kind != PromoterReporterAnchorKind::MotifHit
                || anchor.start_0based < candidate.start_0based
                || anchor.end_0based_exclusive > candidate.end_0based_exclusive
            {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate '{}' must contain one resolved motif_hit anchor",
                        candidate.candidate_id
                    ),
                    cause_chain: vec![],
                });
            }
            let extended_boundary_audit = match member_request.fragment_role {
                PromoterReporterPanelFragmentRole::Core => {
                    if member_request.extended_boundary.is_some() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Core panel member '{}' must not provide extended_boundary",
                                candidate.candidate_id
                            ),
                            cause_chain: vec![],
                        });
                    }
                    None
                }
                PromoterReporterPanelFragmentRole::Extended => {
                    let policy = member_request.extended_boundary.as_ref().ok_or_else(|| {
                        EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Extended panel member '{}' requires an explicit extended_boundary transcript policy",
                                candidate.candidate_id
                            ),
                            cause_chain: vec![],
                        }
                    })?;
                    let (start, end, audit) = Self::promoter_reporter_panel_extended_geometry(
                        source, &candidate, policy,
                    )?;
                    let length = end.saturating_sub(start);
                    if length > candidate_set.fragment_policy.max_fragment_length_bp {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Extended panel member '{}' is {} bp, exceeding the candidate-set maximum of {} bp",
                                candidate.candidate_id,
                                length,
                                candidate_set.fragment_policy.max_fragment_length_bp
                            ),
                            cause_chain: vec![],
                        });
                    }
                    candidate.start_0based = start;
                    candidate.end_0based_exclusive = end;
                    candidate.length_bp = length;
                    Some(audit)
                }
            };
            let base_label = member_request
                .label
                .clone()
                .or_else(|| candidate.gene_label.clone())
                .unwrap_or_else(|| candidate.candidate_id.clone());
            let member_id = format!(
                "{}_{:02}_{}",
                request.panel_id,
                index + 1,
                Self::promoter_reporter_panel_slug(&base_label)
            );
            if !seen_ids.insert(member_id.clone()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Panel member id collision for '{}'", member_id),
                    cause_chain: vec![],
                });
            }
            let wild_type_seq_id = format!("{}_wt", member_id);
            let mutant_seq_id = format!("{}_mut", member_id);
            for seq_id in [&wild_type_seq_id, &mutant_seq_id] {
                if self.state.sequences.contains_key(seq_id) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Planned sequence id '{}' already exists; choose another panel_id",
                            seq_id
                        ),
                        cause_chain: vec![],
                    });
                }
            }
            let motif_forward_strand = !anchor
                .strand
                .as_deref()
                .map(|strand| {
                    matches!(
                        strand.trim().to_ascii_lowercase().as_str(),
                        "-" | "reverse" | "minus"
                    )
                })
                .unwrap_or(false);
            let candidate_start_0based = candidate.start_0based;
            resolved.push(PromoterReporterPanelResolvedMember {
                request: member_request.clone(),
                candidate_set_binding: PromoterReporterPanelSourceBinding {
                    path: member_request.candidate_set_path.clone(),
                    sha256: sha256_prefixed_bytes(&bytes),
                },
                candidate,
                source_seq_id: candidate_set.seq_id.clone(),
                source_sequence_sha256: sha256_prefixed_str(&source.get_forward_string()),
                member_id,
                label: base_label,
                motif_start_in_fragment_0based: anchor.start_0based - candidate_start_0based,
                motif_end_in_fragment_0based_exclusive: anchor.end_0based_exclusive
                    - candidate_start_0based,
                motif_forward_strand,
                wild_type_seq_id,
                mutant_seq_id,
                extended_boundary_audit,
            });
        }
        Ok(resolved)
    }

    fn simulate_promoter_reporter_panel(
        &self,
        request: &PromoterReporterPanelRequest,
        resolved_members: &[PromoterReporterPanelResolvedMember],
        mcs_start_0based: usize,
        mcs_end_0based_exclusive: usize,
    ) -> Result<PromoterReporterPanelSimulation, EngineError> {
        let mut detached = self.fork_detached_execution();
        let mut workflow = vec![];
        let mut created_seq_ids = vec![];
        let mut members = vec![];

        for resolved in resolved_members {
            let extract = Operation::ExtractRegion {
                input: resolved.source_seq_id.clone(),
                from: resolved.candidate.start_0based,
                to: resolved.candidate.end_0based_exclusive,
                output_id: Some(resolved.wild_type_seq_id.clone()),
            };
            let extract_result = Self::apply_promoter_reporter_panel_step(
                &mut detached,
                &mut workflow,
                "extract_promoter_fragment",
                Some(&resolved.member_id),
                extract,
            )?;
            Self::extend_unique(&mut created_seq_ids, &extract_result.created_seq_ids);
            if extract_result.created_seq_ids.as_slice() != [resolved.wild_type_seq_id.as_str()] {
                return Err(EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Panel extraction for '{}' did not create exact id '{}'",
                        resolved.member_id, resolved.wild_type_seq_id
                    ),
                    cause_chain: vec![],
                });
            }
            let wild_type_sequence = detached
                .engine()
                .state
                .sequences
                .get(&resolved.wild_type_seq_id)
                .map(DNAsequence::get_forward_string)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Extracted panel fragment '{}' disappeared",
                        resolved.wild_type_seq_id
                    ),
                    cause_chain: vec![],
                })?;
            let mutation = detached.engine().design_promoter_reporter_panel_mutation(
                request.mutation_policy,
                &wild_type_sequence,
                resolved.motif_start_in_fragment_0based,
                resolved.motif_end_in_fragment_0based_exclusive,
                resolved.motif_forward_strand,
                &[],
            )?;
            let branch_result = Self::apply_promoter_reporter_panel_step(
                &mut detached,
                &mut workflow,
                "branch_wild_type_for_mutation",
                Some(&resolved.member_id),
                Operation::Branch {
                    input: resolved.wild_type_seq_id.clone(),
                    output_id: Some(resolved.mutant_seq_id.clone()),
                },
            )?;
            Self::extend_unique(&mut created_seq_ids, &branch_result.created_seq_ids);
            Self::replace_panel_mutant_sequence(
                detached.engine_mut(),
                &resolved.wild_type_seq_id,
                &resolved.mutant_seq_id,
                &mutation,
            )?;
            Self::record_promoter_reporter_panel_action(
                &mut workflow,
                "apply_stated_p53_family_rule",
                Some(&resolved.member_id),
                json!({"stated_p53_family_motif_disruption": mutation}),
                vec![resolved.mutant_seq_id.clone()],
            )?;

            members.push(PromoterReporterPanelMemberProposal {
                member_id: resolved.member_id.clone(),
                label: resolved.label.clone(),
                fragment_role: resolved.request.fragment_role,
                candidate_set: resolved.candidate_set_binding.clone(),
                candidate_id: resolved.candidate.candidate_id.clone(),
                source_seq_id: resolved.source_seq_id.clone(),
                source_sequence_sha256: resolved.source_sequence_sha256.clone(),
                gene_label: resolved.candidate.gene_label.clone(),
                promoter_class_id: resolved.candidate.promoter_class_id.clone(),
                transcript_ids: if resolved.candidate.transcript_ids.is_empty() {
                    vec![resolved.candidate.transcript_id.clone()]
                } else {
                    resolved.candidate.transcript_ids.clone()
                },
                strand: resolved.candidate.strand.clone(),
                fragment_start_0based: resolved.candidate.start_0based,
                fragment_end_0based_exclusive: resolved.candidate.end_0based_exclusive,
                fragment_length_bp: resolved.candidate.length_bp,
                extended_boundary_audit: resolved.extended_boundary_audit.clone(),
                motif_start_in_fragment_0based: resolved.motif_start_in_fragment_0based,
                motif_end_in_fragment_0based_exclusive: resolved
                    .motif_end_in_fragment_0based_exclusive,
                motif_forward_strand: resolved.motif_forward_strand,
                wild_type_seq_id: resolved.wild_type_seq_id.clone(),
                wild_type_sha256: sha256_prefixed_str(&wild_type_sequence),
                mutant_seq_id: resolved.mutant_seq_id.clone(),
                mutant_sha256: sha256_prefixed_str(&mutation.mutant_sequence),
                mutation,
                warnings: resolved
                    .extended_boundary_audit
                    .as_ref()
                    .map(|audit| {
                        audit
                            .warnings
                            .iter()
                            .map(|warning| warning.detail.clone())
                            .collect()
                    })
                    .unwrap_or_default(),
            });
        }

        let insert_seq_ids = members
            .iter()
            .flat_map(|member| {
                [
                    member.wild_type_seq_id.clone(),
                    member.mutant_seq_id.clone(),
                ]
            })
            .collect::<Vec<_>>();
        let cloning_strategy = detached
            .engine()
            .restriction_cloning_panel_strategy(&request.vector_seq_id, &insert_seq_ids)?;
        let selected_enzymes = cloning_strategy
            .selected_pair
            .as_ref()
            .map(|pair| vec![pair.forward_enzyme.clone(), pair.reverse_enzyme.clone()])
            .unwrap_or_default();
        for member in &mut members {
            let wild_type = detached
                .engine()
                .state
                .sequences
                .get(&member.wild_type_seq_id)
                .map(DNAsequence::get_forward_string)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Panel insert '{}' disappeared", member.wild_type_seq_id),
                    cause_chain: vec![],
                })?;
            let audited = detached.engine().design_promoter_reporter_panel_mutation(
                request.mutation_policy,
                &wild_type,
                member.motif_start_in_fragment_0based,
                member.motif_end_in_fragment_0based_exclusive,
                member.motif_forward_strand,
                &selected_enzymes,
            )?;
            if audited.mutant_sequence != member.mutation.mutant_sequence {
                return Err(EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Restriction audit changed the stated mutant for '{}'",
                        member.member_id
                    ),
                    cause_chain: vec![],
                });
            }
            if !audited.no_stronger_p53_family_hit_near_edit {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Stated motif edit for '{}' creates a stronger local p53-family PWM window",
                        member.member_id
                    ),
                    cause_chain: vec![],
                });
            }
            if audited.creates_selected_restriction_site {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Stated motif edit for '{}' creates a selected cloning-enzyme site",
                        member.member_id
                    ),
                    cause_chain: vec![],
                });
            }
            member.mutation = audited;
        }

        let products = match cloning_strategy.strategy {
            PromoterReporterPanelCloningStrategy::DirectionalRestriction => self
                .simulate_directional_promoter_reporter_products(
                    &mut detached,
                    request,
                    &members,
                    &cloning_strategy,
                    &mut workflow,
                    &mut created_seq_ids,
                )?,
            PromoterReporterPanelCloningStrategy::Gibson => self
                .simulate_gibson_promoter_reporter_products(
                    &mut detached,
                    request,
                    &members,
                    mcs_start_0based,
                    mcs_end_0based_exclusive,
                    &mut workflow,
                    &mut created_seq_ids,
                )?,
        };

        Ok(PromoterReporterPanelSimulation {
            detached,
            members,
            cloning_strategy,
            workflow,
            products,
            created_seq_ids,
        })
    }

    fn design_promoter_reporter_panel_mutation(
        &self,
        policy: PromoterReporterPanelMutationPolicy,
        source_sequence: &str,
        motif_start_0based: usize,
        motif_end_0based_exclusive: usize,
        motif_forward_strand: bool,
        selected_enzymes: &[String],
    ) -> Result<P53FamilyMotifDisruptionReport, EngineError> {
        match policy {
            PromoterReporterPanelMutationPolicy::P53FamilyCoreDisruptionV1 => self
                .design_p53_family_motif_disruption(
                    source_sequence,
                    motif_start_0based,
                    motif_end_0based_exclusive,
                    motif_forward_strand,
                    selected_enzymes,
                ),
            PromoterReporterPanelMutationPolicy::Unspecified => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Promoter-reporter panel mutation policy is unspecified".to_string(),
                cause_chain: vec![],
            }),
        }
    }

    fn replace_panel_mutant_sequence(
        engine: &mut GentleEngine,
        wild_type_seq_id: &str,
        mutant_seq_id: &str,
        mutation: &P53FamilyMotifDisruptionReport,
    ) -> Result<(), EngineError> {
        let wild_type = engine
            .state
            .sequences
            .get(wild_type_seq_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::Internal,
                message: format!("Wild-type fragment '{}' disappeared", wild_type_seq_id),
                cause_chain: vec![],
            })?;
        let mut mutant =
            DNAsequence::from_sequence(&mutation.mutant_sequence).map_err(|error| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not materialize panel mutant: {error}"),
                cause_chain: vec![],
            })?;
        *mutant.features_mut() = wild_type.features().clone();
        mutant.features_mut().push(gb_io::seq::Feature {
            kind: "misc_feature".into(),
            location: gb_io::seq::Location::simple_range(
                mutation
                    .motif_start_0based
                    .try_into()
                    .map_err(|_| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Motif start does not fit GenBank coordinates".to_string(),
                        cause_chain: vec![],
                    })?,
                mutation
                    .motif_end_0based_exclusive
                    .try_into()
                    .map_err(|_| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Motif end does not fit GenBank coordinates".to_string(),
                        cause_chain: vec![],
                    })?,
            ),
            qualifiers: vec![
                (
                    "label".into(),
                    Some("stated-rule p53-family motif mutant".to_string()),
                ),
                (
                    "note".into(),
                    Some(
                        "sequence-only disruption proposal; functional effect not established"
                            .to_string(),
                    ),
                ),
                (
                    "gentle_mutation_schema".into(),
                    Some(P53_FAMILY_MOTIF_DISRUPTION_SCHEMA.to_string()),
                ),
            ],
        });
        mutant.set_name(format!("Mutant branch of {wild_type_seq_id}"));
        mutant.set_circular(false);
        Self::prepare_sequence(&mut mutant);
        engine
            .state_mut()
            .sequences
            .insert(mutant_seq_id.to_string(), mutant);
        Ok(())
    }

    fn apply_promoter_reporter_panel_step(
        detached: &mut DetachedEngineExecution,
        workflow: &mut Vec<PromoterReporterPanelWorkflowStep>,
        stage: &str,
        member_id: Option<&str>,
        operation: Operation,
    ) -> Result<OpResult, EngineError> {
        let operation_value = serde_json::to_value(&operation).map_err(|error| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize panel workflow operation: {error}"),
            cause_chain: vec![],
        })?;
        let result = detached.engine_mut().apply(operation)?;
        Self::record_promoter_reporter_panel_action(
            workflow,
            stage,
            member_id,
            operation_value,
            result.created_seq_ids.clone(),
        )?;
        Ok(result)
    }

    fn record_promoter_reporter_panel_action(
        workflow: &mut Vec<PromoterReporterPanelWorkflowStep>,
        stage: &str,
        member_id: Option<&str>,
        operation: serde_json::Value,
        expected_created_seq_ids: Vec<String>,
    ) -> Result<(), EngineError> {
        let canonical = Self::construct_reasoning_canonical_json(
            &operation,
            "promoter-reporter panel workflow step",
        )?;
        workflow.push(PromoterReporterPanelWorkflowStep {
            ordinal: workflow.len() + 1,
            stage: stage.to_string(),
            member_id: member_id.map(str::to_string),
            operation,
            operation_sha256: sha256_prefixed_str(&canonical),
            expected_created_seq_ids,
        });
        Ok(())
    }

    fn extend_unique(target: &mut Vec<String>, values: &[String]) {
        for value in values {
            if !target.contains(value) {
                target.push(value.clone());
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn simulate_directional_promoter_reporter_products(
        &self,
        detached: &mut DetachedEngineExecution,
        request: &PromoterReporterPanelRequest,
        members: &[PromoterReporterPanelMemberProposal],
        cloning_strategy: &PromoterReporterPanelCloningStrategyReport,
        workflow: &mut Vec<PromoterReporterPanelWorkflowStep>,
        created_seq_ids: &mut Vec<String>,
    ) -> Result<Vec<PromoterReporterPanelProductProposal>, EngineError> {
        let pair = cloning_strategy
            .selected_pair
            .as_ref()
            .ok_or_else(|| EngineError {
                code: ErrorCode::Internal,
                message: "Directional panel strategy has no selected enzyme pair".to_string(),
                cause_chain: vec![],
            })?;
        let enzymes = vec![pair.forward_enzyme.clone(), pair.reverse_enzyme.clone()];
        let vector_digest = Self::apply_promoter_reporter_panel_step(
            detached,
            workflow,
            "digest_destination_vector",
            None,
            Operation::Digest {
                input: request.vector_seq_id.clone(),
                enzymes: enzymes.clone(),
                output_prefix: Some(format!("{}_vector_digest", request.panel_id)),
            },
        )?;
        Self::extend_unique(created_seq_ids, &vector_digest.created_seq_ids);
        let vector_backbone_seq_id = Self::largest_created_sequence(
            detached.engine(),
            &vector_digest.created_seq_ids,
            "digested destination-vector backbone",
        )?;

        let mut products = vec![];
        for member in members {
            for (allele, insert_seq_id) in [
                ("wild_type", member.wild_type_seq_id.as_str()),
                ("mutant", member.mutant_seq_id.as_str()),
            ] {
                let insert_length = detached
                    .engine()
                    .state
                    .sequences
                    .get(insert_seq_id)
                    .map(DNAsequence::len)
                    .unwrap_or_default();
                const TERMINAL_PRIMER_WINDOW_BP: usize = 25;
                if insert_length <= TERMINAL_PRIMER_WINDOW_BP * 2 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Panel insert '{}' is too short for deterministic cloning-primer design",
                            insert_seq_id
                        ),
                        cause_chain: vec![],
                    });
                }
                let report_id = format!("{}_{}_primers", member.member_id, allele);
                let design = Operation::DesignPrimerPairs {
                    template: insert_seq_id.to_string(),
                    roi_start_0based: TERMINAL_PRIMER_WINDOW_BP,
                    roi_end_0based: insert_length - TERMINAL_PRIMER_WINDOW_BP,
                    forward: PrimerDesignSideConstraint {
                        min_length: 18,
                        max_length: 25,
                        location_0based: Some(0),
                        min_tm_c: 0.0,
                        max_tm_c: 100.0,
                        min_gc_fraction: 0.0,
                        max_gc_fraction: 1.0,
                        max_anneal_hits: 1_000,
                        ..PrimerDesignSideConstraint::default()
                    },
                    reverse: PrimerDesignSideConstraint {
                        min_length: 18,
                        max_length: 25,
                        start_0based: Some(insert_length - TERMINAL_PRIMER_WINDOW_BP),
                        end_0based: Some(insert_length),
                        min_tm_c: 0.0,
                        max_tm_c: 100.0,
                        min_gc_fraction: 0.0,
                        max_gc_fraction: 1.0,
                        max_anneal_hits: 1_000,
                        ..PrimerDesignSideConstraint::default()
                    },
                    pair_constraints: PrimerDesignPairConstraint {
                        fixed_amplicon_start_0based: Some(0),
                        fixed_amplicon_end_0based_exclusive: Some(insert_length),
                        ..PrimerDesignPairConstraint::default()
                    },
                    min_amplicon_bp: insert_length.saturating_sub(5),
                    max_amplicon_bp: insert_length.saturating_add(5),
                    max_tm_delta_c: Some(100.0),
                    max_pairs: Some(10),
                    report_id: Some(report_id.clone()),
                };
                Self::apply_promoter_reporter_panel_step(
                    detached,
                    workflow,
                    "design_cloning_primers",
                    Some(&member.member_id),
                    design,
                )?;
                let handoff_result = Self::apply_promoter_reporter_panel_step(
                    detached,
                    workflow,
                    "prepare_directional_cloning_handoff",
                    Some(&member.member_id),
                    Operation::PrepareRestrictionCloningPcrHandoff {
                        template: insert_seq_id.to_string(),
                        primer_report_id: report_id.clone(),
                        pair_index: 0,
                        destination_vector_seq_id: request.vector_seq_id.clone(),
                        mode: RestrictionCloningPcrHandoffMode::DirectedPair,
                        forward_enzyme: pair.forward_enzyme.clone(),
                        reverse_enzyme: Some(pair.reverse_enzyme.clone()),
                        forward_leader_5prime: Some("GCGC".to_string()),
                        reverse_leader_5prime: Some("GCGC".to_string()),
                    },
                )?;
                Self::extend_unique(created_seq_ids, &handoff_result.created_seq_ids);
                let handoff = detached
                    .engine()
                    .list_restriction_cloning_pcr_handoffs()
                    .into_iter()
                    .find(|summary| {
                        summary.primer_report_id == report_id && summary.template == insert_seq_id
                    })
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::Internal,
                        message: format!(
                            "Restriction handoff for '{}' was not persisted",
                            insert_seq_id
                        ),
                        cause_chain: vec![],
                    })?;
                let handoff = detached
                    .engine()
                    .get_restriction_cloning_pcr_handoff(&handoff.report_id)?;
                let primer_seq_ids = vec![
                    handoff.extended_forward_seq_id.clone(),
                    handoff.extended_reverse_seq_id.clone(),
                ];
                let primer_sequences = primer_seq_ids
                    .iter()
                    .map(|seq_id| {
                        detached
                            .engine()
                            .state
                            .sequences
                            .get(seq_id)
                            .map(DNAsequence::get_forward_string)
                            .ok_or_else(|| EngineError {
                                code: ErrorCode::Internal,
                                message: format!(
                                    "Cloning primer sequence '{}' disappeared",
                                    seq_id
                                ),
                                cause_chain: vec![],
                            })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let insert_digest = Self::apply_promoter_reporter_panel_step(
                    detached,
                    workflow,
                    "digest_tailed_promoter_amplicon",
                    Some(&member.member_id),
                    Operation::Digest {
                        input: handoff.tailed_amplicon_seq_id.clone(),
                        enzymes: enzymes.clone(),
                        output_prefix: Some(format!(
                            "{}_{}_insert_digest",
                            member.member_id, allele
                        )),
                    },
                )?;
                Self::extend_unique(created_seq_ids, &insert_digest.created_seq_ids);
                let digested_insert_seq_id = Self::largest_created_sequence(
                    detached.engine(),
                    &insert_digest.created_seq_ids,
                    "digested promoter insert",
                )?;
                let product_seq_id = format!("{}_{}_pgl", member.member_id, allele);
                let ligation = Self::apply_promoter_reporter_panel_step(
                    detached,
                    workflow,
                    "ligate_directional_reporter_construct",
                    Some(&member.member_id),
                    Operation::Ligation {
                        inputs: vec![vector_backbone_seq_id.clone(), digested_insert_seq_id],
                        circularize_if_possible: true,
                        output_id: Some(product_seq_id.clone()),
                        protocol: LigationProtocol::Sticky,
                        output_prefix: None,
                        unique: Some(true),
                    },
                )?;
                Self::extend_unique(created_seq_ids, &ligation.created_seq_ids);
                let final_scan = Self::apply_promoter_reporter_panel_step(
                    detached,
                    workflow,
                    "scan_final_reporter_construct",
                    Some(&member.member_id),
                    Operation::FindRestrictionSites {
                        target: SequenceScanTarget::SeqId {
                            seq_id: product_seq_id.clone(),
                            span_start_0based: None,
                            span_end_0based_exclusive: None,
                        },
                        enzymes: Self::promoter_reporter_panel_all_enzyme_names(),
                        max_sites_per_enzyme: None,
                        include_cut_geometry: false,
                        path: None,
                    },
                )?
                .restriction_site_scan
                .ok_or_else(|| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Final restriction scan returned no report for '{}'",
                        product_seq_id
                    ),
                    cause_chain: vec![],
                })?;
                products.push(Self::promoter_reporter_panel_product(
                    detached.engine(),
                    &member.member_id,
                    allele,
                    insert_seq_id,
                    &product_seq_id,
                    &request.vector_seq_id,
                    &handoff.tailed_amplicon_seq_id,
                    "directional_restriction_digest_ligation",
                    "validated_by_digest_end_compatibility_and_sticky_ligation",
                    primer_seq_ids,
                    primer_sequences,
                    &final_scan,
                )?);
            }
        }
        Ok(products)
    }

    #[allow(clippy::too_many_arguments)]
    fn simulate_gibson_promoter_reporter_products(
        &self,
        detached: &mut DetachedEngineExecution,
        request: &PromoterReporterPanelRequest,
        members: &[PromoterReporterPanelMemberProposal],
        mcs_start_0based: usize,
        mcs_end_0based_exclusive: usize,
        workflow: &mut Vec<PromoterReporterPanelWorkflowStep>,
        created_seq_ids: &mut Vec<String>,
    ) -> Result<Vec<PromoterReporterPanelProductProposal>, EngineError> {
        let mut products = vec![];
        for member in members {
            for (allele, insert_seq_id) in [
                ("wild_type", member.wild_type_seq_id.as_str()),
                ("mutant", member.mutant_seq_id.as_str()),
            ] {
                let product_seq_id = format!("{}_{}_pgl", member.member_id, allele);
                let plan = Self::single_insert_panel_gibson_plan(
                    request,
                    &member.member_id,
                    allele,
                    insert_seq_id,
                    &product_seq_id,
                    mcs_start_0based,
                    mcs_end_0based_exclusive,
                );
                let operation = Operation::ApplyGibsonAssemblyPlan {
                    plan_json: serde_json::to_string(&plan).map_err(|error| EngineError {
                        code: ErrorCode::Internal,
                        message: format!("Could not serialize panel Gibson plan: {error}"),
                        cause_chain: vec![],
                    })?,
                };
                let result = Self::apply_promoter_reporter_panel_step(
                    detached,
                    workflow,
                    "apply_gibson_reporter_construct",
                    Some(&member.member_id),
                    operation,
                )?;
                Self::extend_unique(created_seq_ids, &result.created_seq_ids);
                let mut primer_seq_ids = result
                    .created_seq_ids
                    .iter()
                    .filter(|seq_id| seq_id.to_ascii_lowercase().contains("primer"))
                    .cloned()
                    .collect::<Vec<_>>();
                primer_seq_ids.sort();
                let primer_sequences = primer_seq_ids
                    .iter()
                    .map(|seq_id| {
                        detached
                            .engine()
                            .state
                            .sequences
                            .get(seq_id)
                            .map(DNAsequence::get_forward_string)
                            .ok_or_else(|| EngineError {
                                code: ErrorCode::Internal,
                                message: format!("Gibson primer '{}' disappeared", seq_id),
                                cause_chain: vec![],
                            })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let final_scan = Self::apply_promoter_reporter_panel_step(
                    detached,
                    workflow,
                    "scan_final_reporter_construct",
                    Some(&member.member_id),
                    Operation::FindRestrictionSites {
                        target: SequenceScanTarget::SeqId {
                            seq_id: product_seq_id.clone(),
                            span_start_0based: None,
                            span_end_0based_exclusive: None,
                        },
                        enzymes: Self::promoter_reporter_panel_all_enzyme_names(),
                        max_sites_per_enzyme: None,
                        include_cut_geometry: false,
                        path: None,
                    },
                )?
                .restriction_site_scan
                .ok_or_else(|| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Final restriction scan returned no report for '{}'",
                        product_seq_id
                    ),
                    cause_chain: vec![],
                })?;
                products.push(Self::promoter_reporter_panel_product(
                    detached.engine(),
                    &member.member_id,
                    allele,
                    insert_seq_id,
                    &product_seq_id,
                    &request.vector_seq_id,
                    insert_seq_id,
                    "gibson_assembly",
                    "validated_by_apply_gibson_assembly_plan_junction_policy",
                    primer_seq_ids,
                    primer_sequences,
                    &final_scan,
                )?);
            }
        }
        Ok(products)
    }

    #[allow(clippy::too_many_arguments)]
    fn single_insert_panel_gibson_plan(
        request: &PromoterReporterPanelRequest,
        member_id: &str,
        allele: &str,
        insert_seq_id: &str,
        product_seq_id: &str,
        mcs_start_0based: usize,
        mcs_end_0based_exclusive: usize,
    ) -> GibsonAssemblyPlan {
        let fragment_id = format!("{}_{}_insert", member_id, allele);
        let junction_left = format!("{}_{}_left", member_id, allele);
        let junction_right = format!("{}_{}_right", member_id, allele);
        let destination_left = format!("{}_{}_dest_left", member_id, allele);
        let destination_right = format!("{}_{}_dest_right", member_id, allele);
        let mut validation_policy = GibsonPlanValidationPolicy::default();
        validation_policy.require_unambiguous_destination_opening = true;
        validation_policy.require_distinct_terminal_junctions = true;
        validation_policy.adjacency_overlap_mismatch = "error".to_string();
        validation_policy.design_targets.overlap_bp_min = 18;
        validation_policy.design_targets.overlap_bp_max = 30;
        validation_policy.design_targets.minimum_overlap_tm_celsius = 45.0;
        validation_policy
            .design_targets
            .priming_segment_tm_min_celsius = 45.0;
        validation_policy
            .design_targets
            .priming_segment_tm_max_celsius = 75.0;
        validation_policy.design_targets.max_anneal_hits = 1_000;
        validation_policy.uniqueness_checks.destination_context = "warn".to_string();
        validation_policy.uniqueness_checks.participating_fragments = "warn".to_string();
        GibsonAssemblyPlan {
            schema: GIBSON_ASSEMBLY_PLAN_SCHEMA.to_string(),
            id: format!("{}_{}_gibson", member_id, allele),
            title: format!("{} {} promoter reporter", member_id, allele),
            summary: "Sequence-model one promoter insert in the validated reporter-vector MCS"
                .to_string(),
            destination: GibsonPlanDestination {
                seq_id: request.vector_seq_id.clone(),
                topology_before_opening: "circular".to_string(),
                opening: GibsonPlanOpening {
                    mode: "defined_site".to_string(),
                    label: "validated multiple cloning region".to_string(),
                    start_0based: Some(mcs_start_0based),
                    end_0based_exclusive: Some(mcs_end_0based_exclusive),
                    left_end_id: destination_left.clone(),
                    right_end_id: destination_right.clone(),
                    uniqueness_requirement: "must_be_unambiguous".to_string(),
                },
            },
            product: GibsonPlanProduct {
                topology: "circular".to_string(),
                output_id_hint: product_seq_id.to_string(),
            },
            fragments: vec![GibsonPlanFragment {
                id: fragment_id.clone(),
                seq_id: insert_seq_id.to_string(),
                role: "promoter_insert".to_string(),
                orientation: "forward".to_string(),
                left_end_strategy: Some(GibsonPlanEndStrategy {
                    mode: "primer_added_overlap".to_string(),
                    target_junction_id: junction_left.clone(),
                }),
                right_end_strategy: Some(GibsonPlanEndStrategy {
                    mode: "primer_added_overlap".to_string(),
                    target_junction_id: junction_right.clone(),
                }),
                source_span_1based: None,
            }],
            assembly_order: vec![
                GibsonPlanAssemblyMember {
                    kind: "destination_end".to_string(),
                    id: destination_left.clone(),
                },
                GibsonPlanAssemblyMember {
                    kind: "fragment".to_string(),
                    id: fragment_id.clone(),
                },
                GibsonPlanAssemblyMember {
                    kind: "destination_end".to_string(),
                    id: destination_right.clone(),
                },
            ],
            junctions: vec![
                GibsonPlanJunction {
                    id: junction_left.clone(),
                    left_member: GibsonPlanAssemblyMember {
                        kind: "destination_end".to_string(),
                        id: destination_left,
                    },
                    right_member: GibsonPlanAssemblyMember {
                        kind: "fragment".to_string(),
                        id: fragment_id.clone(),
                    },
                    required_overlap_bp: Some(20),
                    overlap_partition: Some(GibsonPlanOverlapPartition {
                        left_member_bp: 20,
                        right_member_bp: 0,
                    }),
                    overlap_source: "derive_from_destination_left_flank".to_string(),
                    distinct_from: vec![junction_right.clone()],
                },
                GibsonPlanJunction {
                    id: junction_right.clone(),
                    left_member: GibsonPlanAssemblyMember {
                        kind: "fragment".to_string(),
                        id: fragment_id,
                    },
                    right_member: GibsonPlanAssemblyMember {
                        kind: "destination_end".to_string(),
                        id: destination_right,
                    },
                    required_overlap_bp: Some(20),
                    overlap_partition: Some(GibsonPlanOverlapPartition {
                        left_member_bp: 0,
                        right_member_bp: 20,
                    }),
                    overlap_source: "derive_from_destination_right_flank".to_string(),
                    distinct_from: vec![junction_left],
                },
            ],
            validation_policy,
            derived_design: None,
        }
    }

    fn largest_created_sequence(
        engine: &GentleEngine,
        seq_ids: &[String],
        label: &str,
    ) -> Result<String, EngineError> {
        seq_ids
            .iter()
            .filter_map(|seq_id| {
                engine
                    .state
                    .sequences
                    .get(seq_id)
                    .map(|dna| (dna.len(), seq_id))
            })
            .max_by(|left, right| left.0.cmp(&right.0).then_with(|| right.1.cmp(left.1)))
            .map(|(_, seq_id)| seq_id.clone())
            .ok_or_else(|| EngineError {
                code: ErrorCode::Internal,
                message: format!("No sequence was created for {label}"),
                cause_chain: vec![],
            })
    }

    fn promoter_reporter_panel_product(
        engine: &GentleEngine,
        member_id: &str,
        allele: &str,
        insert_seq_id: &str,
        product_seq_id: &str,
        vector_seq_id: &str,
        restriction_comparison_insert_seq_id: &str,
        assembly_model: &str,
        junction_validation: &str,
        primer_seq_ids: Vec<String>,
        primer_sequences_5prime_to_3prime: Vec<String>,
        final_restriction_site_scan: &RestrictionSiteScanReport,
    ) -> Result<PromoterReporterPanelProductProposal, EngineError> {
        let product = engine
            .state
            .sequences
            .get(product_seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Expected panel product '{}' was not created",
                    product_seq_id
                ),
                cause_chain: vec![],
            })?;
        if !product.is_circular() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Panel product '{}' is not circular", product_seq_id),
                cause_chain: vec![],
            });
        }
        let luc2_annotation_preserved = product.features().iter().any(|feature| {
            feature.qualifier_values("label").any(|value| {
                value.to_ascii_lowercase().contains("luc2")
                    || value.to_ascii_lowercase().contains("luciferase")
            }) || feature.qualifier_values("gene").any(|value| {
                value.to_ascii_lowercase().contains("luc2")
                    || value.to_ascii_lowercase().contains("luciferase")
            }) || feature.qualifier_values("product").any(|value| {
                value.to_ascii_lowercase().contains("luc2")
                    || value.to_ascii_lowercase().contains("luciferase")
            })
        });
        if !luc2_annotation_preserved {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Panel product '{}' lost the luc2/luciferase annotation",
                    product_seq_id
                ),
                cause_chain: vec![],
            });
        }
        let final_product_audit = Self::promoter_reporter_panel_final_product_audit(
            engine,
            vector_seq_id,
            restriction_comparison_insert_seq_id,
            product_seq_id,
            junction_validation,
            final_restriction_site_scan,
        )?;
        Ok(PromoterReporterPanelProductProposal {
            member_id: member_id.to_string(),
            allele: allele.to_string(),
            insert_seq_id: insert_seq_id.to_string(),
            product_seq_id: product_seq_id.to_string(),
            sequence_sha256: sha256_prefixed_str(&product.get_forward_string()),
            length_bp: product.len(),
            circular: true,
            luc2_annotation_preserved,
            assembly_model: assembly_model.to_string(),
            primer_seq_ids,
            primer_sequences_5prime_to_3prime,
            final_product_audit,
        })
    }

    fn promoter_reporter_panel_restriction_counts(
        report: &RestrictionSiteScanReport,
    ) -> BTreeMap<String, usize> {
        let mut counts = BTreeMap::new();
        for row in &report.rows {
            *counts.entry(row.enzyme_name.clone()).or_insert(0) += 1;
        }
        counts
    }

    fn promoter_reporter_panel_all_enzyme_names() -> Vec<String> {
        let mut enzyme_names = crate::enzymes::active_restriction_enzymes()
            .into_iter()
            .map(|enzyme| enzyme.name)
            .collect::<Vec<_>>();
        enzyme_names.sort();
        enzyme_names.dedup();
        enzyme_names
    }

    fn promoter_reporter_panel_final_product_audit(
        engine: &GentleEngine,
        vector_seq_id: &str,
        insert_seq_id: &str,
        product_seq_id: &str,
        junction_validation: &str,
        product_scan: &RestrictionSiteScanReport,
    ) -> Result<PromoterReporterPanelFinalProductAudit, EngineError> {
        let mut enzyme_names = product_scan.enzymes_scanned.clone();
        enzyme_names.sort();
        enzyme_names.dedup();
        let scan = |seq_id: &str| {
            engine.find_restriction_sites(
                SequenceScanTarget::SeqId {
                    seq_id: seq_id.to_string(),
                    span_start_0based: None,
                    span_end_0based_exclusive: None,
                },
                &enzyme_names,
                None,
                false,
                None,
                None,
            )
        };
        let vector_scan = scan(vector_seq_id)?;
        let insert_scan = scan(insert_seq_id)?;
        let vector_counts = Self::promoter_reporter_panel_restriction_counts(&vector_scan);
        let insert_counts = Self::promoter_reporter_panel_restriction_counts(&insert_scan);
        let product_counts = Self::promoter_reporter_panel_restriction_counts(product_scan);
        let mut unexpected_count_excess_by_enzyme = BTreeMap::new();
        for (enzyme, product_count) in &product_counts {
            let input_count = vector_counts
                .get(enzyme)
                .copied()
                .unwrap_or_default()
                .saturating_add(insert_counts.get(enzyme).copied().unwrap_or_default());
            if *product_count > input_count {
                unexpected_count_excess_by_enzyme
                    .insert(enzyme.clone(), product_count - input_count);
            }
        }
        let warnings = if unexpected_count_excess_by_enzyme.is_empty() {
            vec![]
        } else {
            vec![format!(
                "Final product '{}' has restriction-site count excess versus vector plus insert for {}; review assembly junctions before experimental use",
                product_seq_id,
                unexpected_count_excess_by_enzyme
                    .iter()
                    .map(|(enzyme, count)| format!("{enzyme} (+{count})"))
                    .collect::<Vec<_>>()
                    .join(", ")
            )]
        };
        let mut enzymes_scanned = product_scan.enzymes_scanned.clone();
        enzymes_scanned.sort();
        enzymes_scanned.dedup();
        Ok(PromoterReporterPanelFinalProductAudit {
            junction_validation: junction_validation.to_string(),
            restriction_scan_schema: product_scan.schema.clone(),
            comparison_vector_seq_id: vector_seq_id.to_string(),
            comparison_insert_seq_id: insert_seq_id.to_string(),
            enzymes_scanned,
            matched_site_count: product_scan.matched_site_count,
            site_count_by_enzyme: product_counts,
            sites: product_scan
                .rows
                .iter()
                .map(|row| PromoterReporterPanelRestrictionSiteRow {
                    enzyme_name: row.enzyme_name.clone(),
                    recognition_start_0based: row.recognition_start_0based,
                    recognition_end_0based_exclusive: row.recognition_end_0based_exclusive,
                    forward_strand: row.forward_strand,
                })
                .collect(),
            unexpected_count_excess_by_enzyme,
            comparison_basis:
                "conservative final-product count excess versus vector plus assembly insert"
                    .to_string(),
            warnings,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::{TempDir, tempdir};

    fn promoter_reporter_panel_materialization_fixture()
    -> (TempDir, GentleEngine, PromoterReporterPanelRequest) {
        let temp = tempdir().expect("panel test directory");
        let helper_catalog_path = temp.path().join("helper_vectors.json");
        let candidate_set_path = temp.path().join("candidate_set.json");
        let output_dir = temp.path().join("materialized_panel");
        let fixture_path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("test_files/fixtures/reporter_vectors/synthetic_mcs_backbone.gb");

        let helper_catalog = serde_json::json!({
            "Synthetic panel vector": {
                "description": "Repository-owned synthetic promoter-reporter test vector",
                "sequence_local": fixture_path.to_string_lossy(),
                "annotations_local": fixture_path.to_string_lossy(),
                "usable_as_empty_backbone": true,
                "helper_kind": "plasmid_vector",
                "sequence_expectation": {
                    "schema": crate::genomes::HELPER_VECTOR_SEQUENCE_EXPECTATION_SCHEMA,
                    "provider": "GENtle tests",
                    "product_name": "synthetic MCS backbone",
                    "catalog_number": "SYNTH-MCS-1",
                    "accession_version": "GENTLE_SYNTHETIC_MCS.1",
                    "expected_length_bp": 240,
                    "expected_topology": "circular",
                    "required_features": [
                        {
                            "id": "multiple_cloning_region",
                            "feature_kinds": ["misc_feature"],
                            "qualifier_terms": ["multiple cloning site region"],
                            "expected_start_1based": 1,
                            "expected_end_1based": 70
                        },
                        {
                            "id": "luc2",
                            "feature_kinds": ["CDS"],
                            "qualifier_terms": ["luciferase luc2 marker"],
                            "expected_start_1based": 100,
                            "expected_end_1based": 180
                        }
                    ],
                    "restriction_site_equivalences": [],
                    "provenance": [{
                        "source_id": "synthetic-test-fixture",
                        "source_url": "test_files/fixtures/reporter_vectors/synthetic_mcs_backbone.gb",
                        "asserted_on": "2026-08-28",
                        "note": "Repository-owned deterministic fixture."
                    }]
                }
            }
        });
        fs::write(
            &helper_catalog_path,
            serde_json::to_vec_pretty(&helper_catalog).expect("helper catalog JSON"),
        )
        .expect("write helper catalog");

        let mut vector = crate::dna_sequence::load_from_file(&fixture_path.to_string_lossy())
            .expect("load synthetic panel vector");
        GentleEngine::prepare_sequence(&mut vector);
        let prefix = "ATCGTACGATGCTAGCTACGTTAGCGATCGTACCTGACTGATCGTAGCTAGC";
        let spacer = "TCGATGACCTAGTCGATCGTAGCATGCTACGATCGT";
        let motif = "GGGCATGCCCGGGCATGCCC";
        let suffix = "CTAGTCGATGCTAGCATCGATCGTACGATGCTAGCTACGATCG";
        let source_sequence = format!("{prefix}{spacer}{motif}{suffix}");
        let motif_start = prefix.len() + spacer.len();
        let anchor = PromoterReporterAnchor {
            kind: PromoterReporterAnchorKind::MotifHit,
            label: "synthetic p53-family response element".to_string(),
            start_0based: motif_start,
            end_0based_exclusive: motif_start + motif.len(),
            strand: Some("+".to_string()),
            motif_id: Some("synthetic_p53_family_re".to_string()),
            evidence_kind: "sequence_motif_candidate".to_string(),
            interpretation_tags: vec!["candidate_not_occupancy_proof".to_string()],
            ..PromoterReporterAnchor::default()
        };
        let candidate = PromoterReporterFragmentCandidate {
            candidate_id: "tp73_core".to_string(),
            gene_label: Some("TP73".to_string()),
            transcript_id: "ENST_SYNTHETIC_TP73".to_string(),
            transcript_label: "synthetic TP73 promoter transcript".to_string(),
            strand: "+".to_string(),
            tss_local_0based: 0,
            promoter_class_id: Some("tp73_core_promoter".to_string()),
            transcript_count: Some(1),
            transcript_ids: vec!["ENST_SYNTHETIC_TP73".to_string()],
            grouping_reason: Some("synthetic deterministic panel test".to_string()),
            anchor: Some(anchor.clone()),
            variant_start_0based: motif_start,
            variant_end_0based_exclusive: motif_start + motif.len(),
            start_0based: 0,
            end_0based_exclusive: source_sequence.len(),
            length_bp: source_sequence.len(),
            rank: 1,
            recommended: true,
            rationale: "Synthetic promoter fragment containing a candidate p53-family motif."
                .to_string(),
            ..PromoterReporterFragmentCandidate::default()
        };
        let candidate_set = PromoterReporterCandidateSet {
            schema: PROMOTER_REPORTER_CANDIDATES_SCHEMA.to_string(),
            seq_id: "panel_source".to_string(),
            sequence_length_bp: source_sequence.len(),
            generated_at_unix_ms: 0,
            anchor: Some(anchor),
            chosen_gene_label: Some("TP73".to_string()),
            chosen_transcript_id: Some("ENST_SYNTHETIC_TP73".to_string()),
            chosen_transcript_label: Some("synthetic TP73 promoter transcript".to_string()),
            transcript_ambiguity_status: "resolved_synthetic_fixture".to_string(),
            recommended_candidate_id: candidate.candidate_id.clone(),
            candidates: vec![candidate],
            ..PromoterReporterCandidateSet::default()
        };
        fs::write(
            &candidate_set_path,
            serde_json::to_vec_pretty(&candidate_set).expect("candidate set JSON"),
        )
        .expect("write candidate set");

        let mut source = DNAsequence::from_sequence(&source_sequence).expect("source sequence");
        GentleEngine::prepare_sequence(&mut source);
        let mut state = ProjectState::default();
        state
            .sequences
            .insert("synthetic_vector".to_string(), vector);
        state.sequences.insert("panel_source".to_string(), source);
        let request = PromoterReporterPanelRequest {
            schema: PROMOTER_REPORTER_PANEL_REQUEST_SCHEMA.to_string(),
            panel_id: "TP73 synthetic panel".to_string(),
            vector_seq_id: "synthetic_vector".to_string(),
            vector_catalog_id: "Synthetic panel vector".to_string(),
            helper_catalog_path: Some(helper_catalog_path.to_string_lossy().to_string()),
            mutation_policy: PromoterReporterPanelMutationPolicy::P53FamilyCoreDisruptionV1,
            scientific_caveats: vec![
                "Synthetic test fixture; no functional promoter claim is made.".to_string(),
            ],
            members: vec![PromoterReporterPanelMemberRequest {
                candidate_set_path: candidate_set_path.to_string_lossy().to_string(),
                candidate_id: Some("tp73_core".to_string()),
                fragment_role: PromoterReporterPanelFragmentRole::Core,
                extended_boundary: None,
                label: Some("TP73 core".to_string()),
            }],
            output_dir: output_dir.to_string_lossy().to_string(),
        };
        (temp, GentleEngine::from_state(state), request)
    }

    fn promoter_reporter_panel_extended_materialization_fixture()
    -> (TempDir, GentleEngine, PromoterReporterPanelRequest) {
        let (temp, mut engine, mut request) = promoter_reporter_panel_materialization_fixture();
        let candidate_set_path = Path::new(&request.members[0].candidate_set_path);
        let mut candidate_set: PromoterReporterCandidateSet = serde_json::from_slice(
            &fs::read(candidate_set_path).expect("extended candidate-set bytes"),
        )
        .expect("extended candidate set");

        let mut bases = engine.snapshot().sequences["panel_source"]
            .get_forward_string()
            .into_bytes();
        bases.extend_from_slice(
            b"AGTCTGACCGTATCGGATCAGTGCATACGTCAGGCTATGACCTGATCGTACGATGCTAGTCAGTCCGATGCTACGTAGCATGTCAGC",
        );
        bases.resize(240, b'A');
        bases[170..172].copy_from_slice(b"AT");
        bases[200] = b'G';
        let source_text = String::from_utf8(bases).expect("extended source DNA");
        let mut source = DNAsequence::from_sequence(&source_text).expect("extended panel source");
        source.features_mut().extend([
            extended_boundary_test_feature("mRNA", 0, 240, "ENST_SYNTHETIC_TP73", false),
            extended_boundary_test_feature("exon", 0, 172, "ENST_SYNTHETIC_TP73", false),
            extended_boundary_test_feature("exon", 200, 240, "ENST_SYNTHETIC_TP73", false),
            extended_boundary_test_feature("CDS", 170, 172, "ENST_SYNTHETIC_TP73", false),
            extended_boundary_test_feature("CDS", 200, 240, "ENST_SYNTHETIC_TP73", false),
        ]);
        GentleEngine::prepare_sequence(&mut source);
        engine
            .state_mut()
            .sequences
            .insert("panel_source".to_string(), source);

        candidate_set.sequence_length_bp = source_text.len();
        let candidate = &mut candidate_set.candidates[0];
        candidate.end_0based_exclusive = candidate.variant_end_0based_exclusive + 10;
        candidate.length_bp = candidate.end_0based_exclusive - candidate.start_0based;
        fs::write(
            candidate_set_path,
            serde_json::to_vec_pretty(&candidate_set).expect("extended candidate-set JSON"),
        )
        .expect("write extended candidate set");

        request.members[0].fragment_role = PromoterReporterPanelFragmentRole::Extended;
        request.members[0].extended_boundary = Some(PromoterReporterPanelExtendedBoundaryPolicy {
            kind: PromoterReporterPanelExtendedBoundaryKind::CanonicalCdsStartCodon,
            transcript_id: "ENST_SYNTHETIC_TP73".to_string(),
        });
        (temp, engine, request)
    }

    fn promoter_reporter_panel_vector(sequence: &str) -> DNAsequence {
        let mut dna = DNAsequence::from_sequence(sequence).expect("vector sequence");
        *dna.restriction_enzymes_mut() = crate::enzymes::active_restriction_enzymes();
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "misc_feature".into(),
            location: gb_io::seq::Location::simple_range(
                0,
                sequence.len().try_into().expect("sequence length"),
            ),
            qualifiers: vec![
                ("label".into(), Some("MCS".to_string())),
                ("note".into(), Some("Multiple cloning site".to_string())),
                (
                    "mcs_expected_sites".into(),
                    Some("EcoRI,HindIII,KpnI".to_string()),
                ),
            ],
        });
        dna.update_computed_features();
        dna
    }

    fn promoter_reporter_panel_engine(insert_a: &str, insert_b: &str) -> GentleEngine {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "vector".to_string(),
            promoter_reporter_panel_vector("AAAAGAATTCGGGGGAAGCTTTTTTGGTACCAAAA"),
        );
        state.sequences.insert(
            "insert_a".to_string(),
            DNAsequence::from_sequence(insert_a).expect("insert A"),
        );
        state.sequences.insert(
            "insert_b".to_string(),
            DNAsequence::from_sequence(insert_b).expect("insert B"),
        );
        GentleEngine::from_state(state)
    }

    #[test]
    fn promoter_reporter_panel_selects_first_pair_clear_across_every_insert() {
        let engine = promoter_reporter_panel_engine(
            "ACGTACGTACGTGAATTCACGTACGTACGT",
            "TGCATGCATGCATGCATGCATGCA",
        );
        let report = engine
            .restriction_cloning_panel_strategy(
                "vector",
                &["insert_a".to_string(), "insert_b".to_string()],
            )
            .expect("panel strategy");

        assert_eq!(
            report.strategy,
            PromoterReporterPanelCloningStrategy::DirectionalRestriction
        );
        let selected = report.selected_pair.expect("selected pair");
        assert_eq!(selected.forward_enzyme, "HindIII");
        assert_eq!(selected.reverse_enzyme, "KpnI");
        assert_eq!(report.pair_evaluations.len(), 3);
        assert!(!report.pair_evaluations[0].compatible);
        assert!(!report.pair_evaluations[1].compatible);
        assert!(report.pair_evaluations[2].compatible);
        assert_eq!(
            report.insert_site_summaries[0]
                .site_count_by_enzyme
                .get("EcoRI"),
            Some(&1)
        );
        assert_eq!(
            report.insert_site_summaries[1]
                .site_count_by_enzyme
                .get("EcoRI"),
            Some(&0)
        );
    }

    #[test]
    fn promoter_reporter_panel_falls_back_to_gibson_when_all_pairs_are_blocked() {
        let engine = promoter_reporter_panel_engine(
            "ACGTGAATTCACGTAAGCTTACGTGGTACCACGT",
            "TGCATGCATGCATGCATGCA",
        );
        let report = engine
            .restriction_cloning_panel_strategy(
                "vector",
                &["insert_a".to_string(), "insert_b".to_string()],
            )
            .expect("panel strategy");

        assert_eq!(
            report.strategy,
            PromoterReporterPanelCloningStrategy::Gibson
        );
        assert!(report.selected_pair.is_none());
        assert!(report.pair_evaluations.iter().all(|row| !row.compatible));
        assert!(
            report
                .fallback_reason
                .as_deref()
                .is_some_and(|reason| reason.contains("internal site"))
        );
    }

    #[test]
    fn promoter_reporter_panel_p53_rule_is_length_preserving_and_audited() {
        let engine = GentleEngine::default();
        let report = engine
            .design_p53_family_motif_disruption(
                "GGGCATGCCCGGGCATGCCC",
                0,
                20,
                true,
                &["KpnI".to_string(), "HindIII".to_string()],
            )
            .expect("p53-family disruption");

        assert_eq!(report.mutant_motif_oriented, "GGGAATTCCCGGGAATTCCC");
        assert_eq!(report.changes.len(), 4);
        assert!(report.length_preserved);
        assert_eq!(report.pwm_audits.len(), 3);
        assert!(report.no_stronger_p53_family_hit_near_edit);
        assert!(
            report
                .pwm_audits
                .iter()
                .any(|audit| audit.stronger_mutant_window_count > 0),
            "shifted submaximal windows remain diagnostic even when every local maximum decreases"
        );
        assert!(
            report
                .pwm_audits
                .iter()
                .all(|audit| { audit.max_llr_delta_bits.is_none_or(|delta| delta <= 1.0e-9) })
        );
        assert!(!report.creates_selected_restriction_site);
        assert!(
            report
                .restriction_site_audits
                .iter()
                .all(|audit| audit.newly_created_sites.is_empty())
        );
    }

    #[test]
    fn promoter_reporter_panel_p53_rule_maps_reverse_strand_changes_to_source() {
        let engine = GentleEngine::default();
        let source = GentleEngine::reverse_complement("GGGCATGCCCGGGCATGCCC");
        let report = engine
            .design_p53_family_motif_disruption(&source, 0, 20, false, &[])
            .expect("reverse p53-family disruption");

        assert_eq!(report.original_motif_oriented, "GGGCATGCCCGGGCATGCCC");
        assert_eq!(report.mutant_motif_oriented, "GGGAATTCCCGGGAATTCCC");
        assert_eq!(
            report.mutant_sequence,
            GentleEngine::reverse_complement("GGGAATTCCCGGGAATTCCC")
        );
        assert_eq!(
            report
                .changes
                .iter()
                .map(|change| change.source_position_0based)
                .collect::<Vec<_>>(),
            vec![3, 6, 13, 16]
        );
    }

    #[test]
    fn promoter_reporter_panel_p53_rule_rejects_noncanonical_core() {
        let engine = GentleEngine::default();
        let error = engine
            .design_p53_family_motif_disruption("GGGAATGCCC", 0, 10, true, &[])
            .expect_err("noncanonical core must fail");
        assert_eq!(error.code, ErrorCode::InvalidInput);
        assert!(error.message.contains("expected C"));
    }

    #[test]
    fn promoter_reporter_panel_rejects_an_unspecified_mutation_policy() {
        let request = PromoterReporterPanelRequest {
            schema: PROMOTER_REPORTER_PANEL_REQUEST_SCHEMA.to_string(),
            panel_id: "unscoped_panel".to_string(),
            vector_seq_id: "vector".to_string(),
            vector_catalog_id: "catalog_vector".to_string(),
            members: vec![PromoterReporterPanelMemberRequest {
                candidate_set_path: "unused.json".to_string(),
                ..PromoterReporterPanelMemberRequest::default()
            }],
            output_dir: "unused-output".to_string(),
            ..PromoterReporterPanelRequest::default()
        };
        let error = GentleEngine::default()
            .plan_promoter_reporter_panel(request)
            .expect_err("an omitted mutation policy must fail before sequence work");
        assert_eq!(error.code, ErrorCode::InvalidInput);
        assert!(error.message.contains("explicit mutation_policy"));
        assert!(error.message.contains("p53_family_core_disruption_v1"));
    }

    fn run_promoter_reporter_panel_end_to_end_test(name: &str, test: fn()) {
        std::thread::Builder::new()
            .name(name.to_string())
            .stack_size(16 * 1024 * 1024)
            .spawn(test)
            .expect("spawn promoter panel end-to-end test")
            .join()
            .expect("promoter panel end-to-end test thread");
    }

    #[test]
    fn promoter_reporter_panel_proposal_is_deterministic_read_only_and_exactly_approved() {
        run_promoter_reporter_panel_end_to_end_test(
            "promoter-panel-proposal-test",
            promoter_reporter_panel_proposal_runs_on_expanded_stack,
        );
    }

    fn promoter_reporter_panel_proposal_runs_on_expanded_stack() {
        let (_temp, mut engine, request) = promoter_reporter_panel_materialization_fixture();
        let baseline = engine
            .promoter_reporter_panel_value_sha256(engine.snapshot(), "test baseline")
            .expect("baseline digest");

        let first = engine
            .plan_promoter_reporter_panel(request.clone())
            .expect("first panel proposal");
        let second = engine
            .plan_promoter_reporter_panel(request)
            .expect("second panel proposal");

        assert_eq!(first.proposal_digest, second.proposal_digest);
        assert_eq!(
            first.request.mutation_policy,
            PromoterReporterPanelMutationPolicy::P53FamilyCoreDisruptionV1
        );
        assert!(
            first
                .nonclaims
                .iter()
                .any(|value| value.contains("Synthetic test fixture"))
        );
        assert!(
            first
                .nonclaims
                .iter()
                .all(|value| !value.contains("TGFB1") && !value.contains("SERPINE1"))
        );
        assert_eq!(first.workflow_sha256, second.workflow_sha256);
        assert_eq!(first.products.len(), 2);
        assert_eq!(
            first.cloning_strategy.strategy,
            PromoterReporterPanelCloningStrategy::DirectionalRestriction
        );
        assert!(first.products.iter().all(|product| {
            product.circular
                && product.luc2_annotation_preserved
                && !product.primer_seq_ids.is_empty()
                && product.primer_seq_ids.len() == product.primer_sequences_5prime_to_3prime.len()
                && !product.final_product_audit.junction_validation.is_empty()
                && product.final_product_audit.enzymes_scanned.len() > 100
                && product.final_product_audit.matched_site_count
                    == product.final_product_audit.sites.len()
        }));
        assert_eq!(
            baseline,
            engine
                .promoter_reporter_panel_value_sha256(engine.snapshot(), "test baseline")
                .expect("unchanged baseline digest")
        );

        let error = engine
            .materialize_promoter_reporter_panel(first.clone(), "sha256:not-approved")
            .expect_err("wrong digest must fail");
        assert!(error.message.contains("exact approval digest"));
        assert!(!Path::new(&first.request.output_dir).exists());
        assert_eq!(
            baseline,
            engine
                .promoter_reporter_panel_value_sha256(engine.snapshot(), "test baseline")
                .expect("state after rejected approval")
        );
    }

    #[test]
    fn promoter_reporter_panel_rejects_stale_bound_input_without_side_effects() {
        run_promoter_reporter_panel_end_to_end_test(
            "promoter-panel-stale-input-test",
            promoter_reporter_panel_stale_input_runs_on_expanded_stack,
        );
    }

    fn promoter_reporter_panel_stale_input_runs_on_expanded_stack() {
        let (_temp, mut engine, request) = promoter_reporter_panel_materialization_fixture();
        let proposal = engine
            .plan_promoter_reporter_panel(request)
            .expect("panel proposal");
        let candidate_path = proposal.members[0].candidate_set.path.clone();
        let mut candidate_bytes = fs::read(&candidate_path).expect("candidate bytes");
        candidate_bytes.push(b'\n');
        fs::write(&candidate_path, candidate_bytes).expect("touch candidate binding");
        let baseline = engine
            .promoter_reporter_panel_value_sha256(engine.snapshot(), "test stale baseline")
            .expect("stale baseline digest");

        let error = engine
            .materialize_promoter_reporter_panel(
                proposal.clone(),
                proposal.proposal_digest.as_str(),
            )
            .expect_err("stale candidate source must fail");

        assert!(error.message.contains("stale"));
        assert!(!Path::new(&proposal.request.output_dir).exists());
        assert_eq!(
            baseline,
            engine
                .promoter_reporter_panel_value_sha256(engine.snapshot(), "test stale baseline")
                .expect("state after stale rejection")
        );
    }

    #[test]
    fn promoter_reporter_panel_materializes_exact_products_and_manifest() {
        run_promoter_reporter_panel_end_to_end_test(
            "promoter-panel-materialization-test",
            promoter_reporter_panel_materialization_runs_on_expanded_stack,
        );
    }

    fn promoter_reporter_panel_materialization_runs_on_expanded_stack() {
        let (_temp, mut engine, request) = promoter_reporter_panel_materialization_fixture();
        let proposal = engine
            .plan_promoter_reporter_panel(request)
            .expect("panel proposal");
        let receipt = engine
            .materialize_promoter_reporter_panel(
                proposal.clone(),
                proposal.proposal_digest.as_str(),
            )
            .expect("approved materialization");

        assert_eq!(receipt.approved_proposal_digest, proposal.proposal_digest);
        assert_eq!(receipt.product_seq_ids.len(), 2);
        assert_eq!(receipt.artifact_paths.len(), 5);
        assert!(
            receipt
                .artifact_paths
                .iter()
                .all(|path| Path::new(path).is_file())
        );
        assert_eq!(
            receipt.completed_state_sha256,
            engine
                .promoter_reporter_panel_value_sha256(engine.snapshot(), "completed state")
                .expect("completed state digest")
        );
        for product in &proposal.products {
            let sequence = engine
                .snapshot()
                .sequences
                .get(&product.product_seq_id)
                .expect("materialized product");
            assert!(sequence.is_circular());
            assert!(sequence.features().iter().any(|feature| {
                feature.kind.eq_ignore_ascii_case("CDS")
                    && feature.qualifiers.iter().any(|(key, value)| {
                        matches!(
                            key.to_ascii_lowercase().as_str(),
                            "label" | "gene" | "product"
                        ) && value
                            .as_deref()
                            .is_some_and(|value| value.to_ascii_lowercase().contains("luc2"))
                    })
            }));
        }
        let manifest = fs::read_to_string(&receipt.manifest_path).expect("panel manifest");
        assert!(manifest.contains("primer_seq_ids\tprimer_sequences_5prime_to_3prime"));
        assert!(manifest.contains("junction_validation\tfinal_restriction_site_count"));
        assert!(manifest.contains("tp73_synthetic_panel_01_tp73_core"));
    }

    #[test]
    fn promoter_reporter_split_codon_extended_member_plans_and_materializes_exact_insert() {
        run_promoter_reporter_panel_end_to_end_test(
            "promoter-panel-split-codon-test",
            promoter_reporter_split_codon_extended_member_runs_on_expanded_stack,
        );
    }

    fn promoter_reporter_split_codon_extended_member_runs_on_expanded_stack() {
        let (_temp, mut engine, request) =
            promoter_reporter_panel_extended_materialization_fixture();
        let proposal = engine
            .plan_promoter_reporter_panel(request)
            .expect("split-codon extended proposal");
        let member = &proposal.members[0];
        let audit = member
            .extended_boundary_audit
            .as_ref()
            .expect("extended boundary audit");
        assert_eq!(member.fragment_end_0based_exclusive, 201);
        assert_eq!(audit.cds_start_codon_5prime_to_3prime, "ATG");
        assert_eq!(
            audit.cds_start_codon_source_ranges_0based,
            vec![(170, 172), (200, 201)]
        );

        engine
            .materialize_promoter_reporter_panel(
                proposal.clone(),
                proposal.proposal_digest.as_str(),
            )
            .expect("materialize split-codon extended panel");
        let insert = &engine.snapshot().sequences[&member.wild_type_seq_id];
        assert_eq!(insert.len(), 201);
        assert!(insert.get_forward_string().ends_with('G'));
    }

    fn extended_boundary_test_feature(
        kind: &str,
        start: usize,
        end: usize,
        transcript_id: &str,
        reverse: bool,
    ) -> gb_io::seq::Feature {
        let location = gb_io::seq::Location::simple_range(start as i64, end as i64);
        gb_io::seq::Feature {
            kind: kind.to_string().into(),
            location: if reverse {
                gb_io::seq::Location::Complement(Box::new(location))
            } else {
                location
            },
            qualifiers: vec![("transcript_id".into(), Some(transcript_id.to_string()))],
        }
    }

    fn promoter_reporter_cds_entry_validation_fixture(
        codon: &str,
        qualifier: Option<(&str, &str)>,
    ) -> (
        DNAsequence,
        PromoterReporterFragmentCandidate,
        PromoterReporterPanelExtendedBoundaryPolicy,
    ) {
        let mut bases = vec![b'C'; 120];
        bases[80..83].copy_from_slice(codon.as_bytes());
        let mut source = DNAsequence::from_sequence(&String::from_utf8(bases).unwrap())
            .expect("CDS-entry validation source");
        let mut cds = extended_boundary_test_feature("CDS", 80, 100, "ENST_PHASE_TEST", false);
        if let Some((key, value)) = qualifier {
            cds.qualifiers
                .push((key.to_string().into(), Some(value.to_string())));
        }
        source.features_mut().extend([
            extended_boundary_test_feature("mRNA", 0, 100, "ENST_PHASE_TEST", false),
            extended_boundary_test_feature("exon", 0, 100, "ENST_PHASE_TEST", false),
            cds,
        ]);
        (
            source,
            PromoterReporterFragmentCandidate {
                candidate_id: "phase_test_core".to_string(),
                strand: "+".to_string(),
                variant_start_0based: 20,
                variant_end_0based_exclusive: 25,
                start_0based: 0,
                end_0based_exclusive: 40,
                length_bp: 40,
                ..PromoterReporterFragmentCandidate::default()
            },
            PromoterReporterPanelExtendedBoundaryPolicy {
                kind: PromoterReporterPanelExtendedBoundaryKind::CanonicalCdsStartCodon,
                transcript_id: "ENST_PHASE_TEST".to_string(),
            },
        )
    }

    #[test]
    fn promoter_reporter_extended_boundary_rejects_partial_cds_entry_phase() {
        for (key, value) in [
            ("codon_start", "2"),
            ("codon_start", "3"),
            ("phase", "1"),
            ("phase", "2"),
        ] {
            let (source, candidate, policy) =
                promoter_reporter_cds_entry_validation_fixture("ATG", Some((key, value)));
            let error = GentleEngine::promoter_reporter_panel_extended_geometry(
                &source, &candidate, &policy,
            )
            .expect_err("partial CDS entry must fail closed");
            assert_eq!(error.code, ErrorCode::InvalidInput);
            assert!(
                error.message.contains("5'-partial CDS"),
                "{}",
                error.message
            );
        }
    }

    #[test]
    fn promoter_reporter_extended_boundary_rejects_invalid_or_conflicting_phase() {
        let (source, candidate, policy) =
            promoter_reporter_cds_entry_validation_fixture("ATG", Some(("codon_start", "4")));
        let error =
            GentleEngine::promoter_reporter_panel_extended_geometry(&source, &candidate, &policy)
                .expect_err("invalid codon_start must fail closed");
        assert!(error.message.contains("invalid /codon_start='4'"));

        let (mut source, candidate, policy) =
            promoter_reporter_cds_entry_validation_fixture("ATG", Some(("phase", "1")));
        source.features_mut()[0]
            .qualifiers
            .push(("codon_start".into(), Some("1".to_string())));
        let error =
            GentleEngine::promoter_reporter_panel_extended_geometry(&source, &candidate, &policy)
                .expect_err("conflicting phase annotations must fail closed");
        assert!(error.message.contains("conflicting CDS-entry phase"));
    }

    #[test]
    fn promoter_reporter_extended_boundary_rejects_noncanonical_start_triplet() {
        let (source, candidate, policy) =
            promoter_reporter_cds_entry_validation_fixture("GTG", Some(("codon_start", "1")));
        let error =
            GentleEngine::promoter_reporter_panel_extended_geometry(&source, &candidate, &policy)
                .expect_err("noncanonical start triplet must fail closed");
        assert_eq!(error.code, ErrorCode::InvalidInput);
        assert!(error.message.contains("'GTG', not the canonical ATG"));
    }

    #[test]
    fn promoter_reporter_extended_boundary_rejects_partial_location_marker() {
        let (mut source, candidate, policy) =
            promoter_reporter_cds_entry_validation_fixture("ATG", None);
        source.features_mut()[2].location = gb_io::seq::Location::Range(
            (80, gb_io::seq::Before(true)),
            (100, gb_io::seq::After(false)),
        );
        let error =
            GentleEngine::promoter_reporter_panel_extended_geometry(&source, &candidate, &policy)
                .expect_err("partial CDS location must fail closed");
        assert!(error.message.contains("5'-partial CDS location"));
    }

    #[test]
    fn promoter_reporter_extended_boundary_uses_spliced_utr_and_reports_uorf_and_intron() {
        let mut bases = vec![b'C'; 120];
        bases[5..14].copy_from_slice(b"ATGAAATAA");
        bases[80..83].copy_from_slice(b"ATG");
        let mut source = DNAsequence::from_sequence(&String::from_utf8(bases).unwrap())
            .expect("extended source");
        source.features_mut().extend([
            extended_boundary_test_feature("mRNA", 0, 100, "ENST_TGFB1", false),
            extended_boundary_test_feature("exon", 0, 30, "ENST_TGFB1", false),
            extended_boundary_test_feature("exon", 50, 100, "ENST_TGFB1", false),
            extended_boundary_test_feature("CDS", 80, 100, "ENST_TGFB1", false),
        ]);
        let candidate = PromoterReporterFragmentCandidate {
            candidate_id: "tgfb1_core".to_string(),
            strand: "+".to_string(),
            variant_start_0based: 20,
            variant_end_0based_exclusive: 25,
            start_0based: 0,
            end_0based_exclusive: 40,
            length_bp: 40,
            ..PromoterReporterFragmentCandidate::default()
        };
        let policy = PromoterReporterPanelExtendedBoundaryPolicy {
            kind: PromoterReporterPanelExtendedBoundaryKind::CanonicalCdsStartCodon,
            transcript_id: "ENST_TGFB1".to_string(),
        };

        let (start, end, audit) =
            GentleEngine::promoter_reporter_panel_extended_geometry(&source, &candidate, &policy)
                .expect("extended boundary");

        assert_eq!((start, end), (0, 83));
        assert_eq!(
            audit.five_prime_utr_exon_ranges_0based,
            vec![(0, 30), (50, 80)]
        );
        assert_eq!(audit.five_prime_utr_spliced_length_bp, 60);
        assert_eq!(audit.cds_start_codon_5prime_to_3prime, "ATG");
        assert_eq!(audit.cds_start_codon_source_ranges_0based, vec![(80, 83)]);
        assert!(audit.warnings.iter().any(|warning| {
            warning.kind == PromoterReporterPanelExtendedWarningKind::UpstreamAtg
        }));
        assert!(audit.warnings.iter().any(|warning| {
            warning.kind == PromoterReporterPanelExtendedWarningKind::UpstreamOrf
        }));
        assert!(audit.warnings.iter().any(|warning| {
            warning.kind == PromoterReporterPanelExtendedWarningKind::FivePrimeUtrIntron
        }));
        let round_trip: PromoterReporterPanelExtendedBoundaryAudit = serde_json::from_str(
            &serde_json::to_string(&audit).expect("serialize extended boundary audit"),
        )
        .expect("deserialize extended boundary audit");
        assert_eq!(round_trip, audit);
    }

    #[test]
    fn promoter_reporter_extended_boundary_is_strand_aware() {
        let mut bases = vec![b'G'; 120];
        bases[37..40].copy_from_slice(b"CAT");
        let mut source = DNAsequence::from_sequence(&String::from_utf8(bases).unwrap())
            .expect("reverse extended source");
        source.features_mut().extend([
            extended_boundary_test_feature("mRNA", 20, 120, "ENST_CD44", true),
            extended_boundary_test_feature("exon", 20, 70, "ENST_CD44", true),
            extended_boundary_test_feature("exon", 90, 120, "ENST_CD44", true),
            extended_boundary_test_feature("CDS", 20, 40, "ENST_CD44", true),
        ]);
        let candidate = PromoterReporterFragmentCandidate {
            candidate_id: "cd44_core".to_string(),
            strand: "-".to_string(),
            variant_start_0based: 100,
            variant_end_0based_exclusive: 105,
            start_0based: 80,
            end_0based_exclusive: 120,
            length_bp: 40,
            ..PromoterReporterFragmentCandidate::default()
        };
        let policy = PromoterReporterPanelExtendedBoundaryPolicy {
            kind: PromoterReporterPanelExtendedBoundaryKind::CanonicalCdsStartCodon,
            transcript_id: "ENST_CD44".to_string(),
        };

        let (start, end, audit) =
            GentleEngine::promoter_reporter_panel_extended_geometry(&source, &candidate, &policy)
                .expect("reverse extended boundary");

        assert_eq!((start, end), (37, 120));
        assert_eq!(audit.strand, "-");
        assert_eq!(audit.cds_start_codon_source_start_0based, 37);
        assert_eq!(audit.cds_start_codon_source_end_0based_exclusive, 40);
        assert_eq!(audit.cds_start_codon_5prime_to_3prime, "ATG");
        assert_eq!(audit.cds_start_codon_source_ranges_0based, vec![(37, 40)]);
        assert_eq!(
            audit.five_prime_utr_exon_ranges_0based,
            vec![(90, 120), (40, 70)]
        );
    }

    #[test]
    fn promoter_reporter_extended_boundary_resolves_forward_split_start_codon() {
        let mut bases = vec![b'C'; 130];
        bases[80..82].copy_from_slice(b"AT");
        bases[100] = b'G';
        let mut source = DNAsequence::from_sequence(&String::from_utf8(bases).unwrap())
            .expect("forward split-codon source");
        source.features_mut().extend([
            extended_boundary_test_feature("mRNA", 0, 130, "ENST_SPLIT_FWD", false),
            extended_boundary_test_feature("exon", 0, 82, "ENST_SPLIT_FWD", false),
            extended_boundary_test_feature("exon", 100, 130, "ENST_SPLIT_FWD", false),
            extended_boundary_test_feature("CDS", 80, 82, "ENST_SPLIT_FWD", false),
            extended_boundary_test_feature("CDS", 100, 130, "ENST_SPLIT_FWD", false),
        ]);
        let candidate = PromoterReporterFragmentCandidate {
            candidate_id: "split_forward".to_string(),
            strand: "+".to_string(),
            variant_start_0based: 20,
            variant_end_0based_exclusive: 25,
            start_0based: 0,
            end_0based_exclusive: 40,
            length_bp: 40,
            ..PromoterReporterFragmentCandidate::default()
        };
        let policy = PromoterReporterPanelExtendedBoundaryPolicy {
            kind: PromoterReporterPanelExtendedBoundaryKind::CanonicalCdsStartCodon,
            transcript_id: "ENST_SPLIT_FWD".to_string(),
        };

        let (start, end, audit) =
            GentleEngine::promoter_reporter_panel_extended_geometry(&source, &candidate, &policy)
                .expect("forward split-codon boundary");

        assert_eq!((start, end), (0, 101));
        assert_eq!(audit.cds_start_codon_5prime_to_3prime, "ATG");
        assert_eq!(
            audit.cds_start_codon_source_ranges_0based,
            vec![(80, 82), (100, 101)]
        );
        assert_eq!(audit.cds_start_codon_source_start_0based, 80);
        assert_eq!(audit.cds_start_codon_source_end_0based_exclusive, 101);
    }

    #[test]
    fn promoter_reporter_extended_boundary_resolves_reverse_split_start_codon() {
        let mut bases = vec![b'C'; 160];
        bases[118..120].copy_from_slice(b"AT");
        bases[99] = b'C';
        let mut source = DNAsequence::from_sequence(&String::from_utf8(bases).unwrap())
            .expect("reverse split-codon source");
        source.features_mut().extend([
            extended_boundary_test_feature("mRNA", 80, 160, "ENST_SPLIT_REV", true),
            extended_boundary_test_feature("exon", 80, 100, "ENST_SPLIT_REV", true),
            extended_boundary_test_feature("exon", 118, 160, "ENST_SPLIT_REV", true),
            extended_boundary_test_feature("CDS", 80, 100, "ENST_SPLIT_REV", true),
            extended_boundary_test_feature("CDS", 118, 120, "ENST_SPLIT_REV", true),
        ]);
        let candidate = PromoterReporterFragmentCandidate {
            candidate_id: "split_reverse".to_string(),
            strand: "-".to_string(),
            variant_start_0based: 135,
            variant_end_0based_exclusive: 140,
            start_0based: 130,
            end_0based_exclusive: 155,
            length_bp: 25,
            ..PromoterReporterFragmentCandidate::default()
        };
        let policy = PromoterReporterPanelExtendedBoundaryPolicy {
            kind: PromoterReporterPanelExtendedBoundaryKind::CanonicalCdsStartCodon,
            transcript_id: "ENST_SPLIT_REV".to_string(),
        };

        let (start, end, audit) =
            GentleEngine::promoter_reporter_panel_extended_geometry(&source, &candidate, &policy)
                .expect("reverse split-codon boundary");

        assert_eq!((start, end), (99, 155));
        assert_eq!(audit.cds_start_codon_5prime_to_3prime, "ATG");
        assert_eq!(
            audit.cds_start_codon_source_ranges_0based,
            vec![(118, 120), (99, 100)]
        );
        assert_eq!(audit.cds_start_codon_source_start_0based, 99);
        assert_eq!(audit.cds_start_codon_source_end_0based_exclusive, 120);
    }

    #[test]
    fn promoter_reporter_extended_boundary_requires_exact_transcript_cds_annotation() {
        let source = DNAsequence::from_sequence(&"A".repeat(100)).expect("source");
        let candidate = PromoterReporterFragmentCandidate {
            candidate_id: "serpine1_core".to_string(),
            strand: "+".to_string(),
            variant_start_0based: 10,
            variant_end_0based_exclusive: 20,
            start_0based: 0,
            end_0based_exclusive: 30,
            length_bp: 30,
            ..PromoterReporterFragmentCandidate::default()
        };
        let error = GentleEngine::promoter_reporter_panel_extended_geometry(
            &source,
            &candidate,
            &PromoterReporterPanelExtendedBoundaryPolicy {
                kind: PromoterReporterPanelExtendedBoundaryKind::CanonicalCdsStartCodon,
                transcript_id: "ENST_SERPINE1".to_string(),
            },
        )
        .expect_err("missing transcript annotation must fail");
        assert_eq!(error.code, ErrorCode::NotFound);
        assert!(error.message.contains("exactly one is required"));
    }
}
