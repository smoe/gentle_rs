//! Hash-bound reuse of typed locus evidence; citations alone never imply assessment.

use super::*;
use std::io::Read;

impl GentleEngine {
    pub(super) fn regulatory_fragment_external_dimension(
        &self,
        request: &RegulatoryFragmentPanelRequest,
        fragments: &[ResolvedRegulatoryFragment],
        kind: RegulatoryFragmentEvidenceDimensionKind,
    ) -> Result<RegulatoryFragmentEvidenceDimension, EngineError> {
        let mut observations = vec![];
        let mut all_assessed = true;
        let mut truncated = false;
        for binding in request
            .evidence_bindings
            .iter()
            .filter(|b| b.dimension == kind)
        {
            let Some(path) = &binding.report_path else {
                all_assessed = false;
                continue;
            };
            let invalid = |detail: String| {
                Self::regulatory_fragment_error(
                    "external_evidence_binding_invalid",
                    [&binding.report_id],
                    detail,
                )
            };
            // Hash and parse the same bounded byte buffer, avoiding file replacement races.
            let mut bytes = Vec::new();
            fs::File::open(path)
                .map_err(|e| invalid(e.to_string()))?
                .take(128 * 1024 * 1024 + 1)
                .read_to_end(&mut bytes)
                .map_err(|e| invalid(e.to_string()))?;
            if bytes.len() > 128 * 1024 * 1024
                || crate::digest_utils::sha256_prefixed_bytes(&bytes) != binding.report_sha256
            {
                return Err(invalid(
                    "Locus document exceeds the size limit or its bytes differ from report_sha256"
                        .into(),
                ));
            }
            let document =
                crate::locus_report::LocusDocument::from_json(&bytes).map_err(invalid)?;
            let locus = document.locus();
            if locus.panel_id != binding.report_id {
                return Err(invalid(
                    "report_id must identify the nested locus panel_id".into(),
                ));
            }
            let dna = self
                .state
                .sequences
                .get(&locus.seq_id)
                .ok_or_else(|| invalid("Source sequence is not loaded".into()))?;
            let anchor = self.sequence_genome_anchor_summary(&locus.seq_id)?;
            crate::locus_report::verify_live_binding(
                locus,
                &locus.seq_id,
                &crate::locus_report::sequence_binding(dna, Some(&anchor)),
            )
            .map_err(invalid)?;
            let matching = fragments
                .iter()
                .filter(|f| f.projection.seq_id == locus.seq_id)
                .collect::<Vec<_>>();
            if matching.is_empty()
                || matching.iter().any(|f| {
                    locus.isoform_evidence.assembly
                        != f.binding.region.interval.reference.assembly_name
                        || locus.isoform_evidence.annotation_release.as_deref()
                            != Some(f.binding.reference_release.as_str())
                })
            {
                return Err(invalid("Locus sequence, assembly or annotation release does not match the bound fragments".into()));
            }
            let overlaps = |start: usize, end: usize| -> Result<Vec<String>, EngineError> {
                if start == 0 || start > end || end > dna.len() {
                    return Err(invalid(
                        "Source evidence has invalid local coordinates".into(),
                    ));
                }
                Ok(matching
                    .iter()
                    .filter(|f| {
                        ((start - 1) as u64) < f.projection.local_end_0based_exclusive
                            && end as u64 > f.projection.local_start_0based
                    })
                    .map(|f| f.binding.fragment_id.clone())
                    .collect())
            };
            let mut rows = vec![];
            match kind {
                RegulatoryFragmentEvidenceDimensionKind::EnsemblRegulatoryOverlap => {
                    if let Some(source) = &locus.ensembl_regulation {
                        truncated |= source.source_binding.as_ref().is_some_and(|b| b.truncated);
                        let available = source.availability
                            == gp::GeneLocusEnsemblRegulationAvailability::Available;
                        if available
                            && !source.source_binding.as_ref().is_some_and(|s| {
                                s.content_identity_verified
                                    && Self::regulatory_fragment_valid_sha256(&s.index_sha256)
                                    && Self::regulatory_fragment_valid_sha256(&s.intervals_sha256)
                                    && s.source.assembly_name == locus.isoform_evidence.assembly
                            })
                        {
                            return Err(invalid(
                                "Available Ensembl evidence requires verified source hashes".into(),
                            ));
                        }
                        for row in &source.rows {
                            overlaps(
                                row.displayed_local_start_1based,
                                row.displayed_local_end_1based,
                            )?;
                            if row.assembly_name != locus.isoform_evidence.assembly {
                                return Err(invalid(
                                    "Ensembl feature assembly disagrees with its locus".into(),
                                ));
                            }
                        }
                        if binding.row_id.is_none() {
                            let ids = overlaps(
                                locus.locus_local_start_1based,
                                locus.locus_local_end_1based,
                            )?;
                            rows.push((
                                source.requested_source_id.clone(),
                                ids,
                                available,
                                RegulatoryFragmentExternalEvidence::EnsemblRegulation(
                                    source.clone(),
                                ),
                            ));
                        } else {
                            for row in &source.rows {
                                if binding.row_id.as_deref() != Some(row.feature_id.as_str()) {
                                    continue;
                                }
                                let ids = overlaps(
                                    row.displayed_local_start_1based,
                                    row.displayed_local_end_1based,
                                )?;
                                let mut selected = source.clone();
                                selected.rows = vec![row.clone()];
                                rows.push((
                                    row.feature_id.clone(),
                                    ids,
                                    available,
                                    RegulatoryFragmentExternalEvidence::EnsemblRegulation(selected),
                                ));
                            }
                        }
                    }
                }
                RegulatoryFragmentEvidenceDimensionKind::TfbsModelScoreContext => {
                    for track in &locus.regulatory_score_tracks {
                        if binding
                            .row_id
                            .as_deref()
                            .is_some_and(|id| id != track.track_id)
                        {
                            continue;
                        }
                        let available = track.state == gp::GeneLocusRegulatoryScoreState::Available;
                        if available
                            && (track.input_sequence_id != locus.seq_id
                                || track.input_sequence_sha256.trim_start_matches("sha256:")
                                    != locus
                                        .sequence_binding
                                        .as_ref()
                                        .map(|b| b.sequence_sha256.trim_start_matches("sha256:"))
                                        .unwrap_or("")
                                || track.assembly != locus.isoform_evidence.assembly
                                || track.chromosome != anchor.chromosome
                                || track.anchor_start_1based != anchor.start_1based
                                || track.anchor_end_1based != anchor.end_1based)
                        {
                            return Err(invalid(
                                "Score track identity disagrees with its locus document".into(),
                            ));
                        }
                        for values in [&track.forward_scores, &track.reverse_scores] {
                            if values.is_empty() {
                                continue;
                            }
                            let end = (values.len() - 1)
                                .checked_mul(track.stride_bp)
                                .and_then(|n| n.checked_add(track.track_start_0based))
                                .and_then(|n| n.checked_add(track.window_length_bp));
                            if track.stride_bp == 0
                                || track.window_length_bp == 0
                                || end.is_none_or(|end| end > dna.len())
                                || values.iter().any(|v| !v.is_finite())
                            {
                                return Err(invalid(
                                    "Invalid score vector coordinates or values".into(),
                                ));
                            }
                        }
                        for site in &track.sites {
                            let start = site
                                .local_start_0based
                                .checked_add(1)
                                .ok_or_else(|| invalid("Score site coordinate overflow".into()))?;
                            overlaps(start, site.local_end_0based_exclusive)?;
                        }
                        let ids =
                            overlaps(locus.locus_local_start_1based, locus.locus_local_end_1based)?;
                        rows.push((
                            track.track_id.clone(),
                            ids,
                            available,
                            RegulatoryFragmentExternalEvidence::RegulatoryScore(track.clone()),
                        ));
                    }
                }
                RegulatoryFragmentEvidenceDimensionKind::CutrunAndChromatinContext => {
                    for lane in locus.occupancy_groups.iter().flat_map(|g| &g.lanes) {
                        if binding
                            .row_id
                            .as_deref()
                            .is_some_and(|id| id != lane.lane.lane_id)
                        {
                            continue;
                        }
                        truncated |= lane.lane.interval_count > lane.lane.intervals.len();
                        let available = lane.state == gp::GeneLocusOccupancyLaneState::Available;
                        if available
                            && (!lane
                                .source_sha256
                                .as_deref()
                                .is_some_and(Self::regulatory_fragment_valid_sha256)
                                || lane.source_assembly.as_deref()
                                    != Some(locus.isoform_evidence.assembly.as_str()))
                        {
                            return Err(invalid("Available occupancy evidence requires source hash and matching assembly".into()));
                        }
                        for interval in &lane.lane.intervals {
                            overlaps(interval.local_start_1based, interval.local_end_1based)?;
                        }
                        let ids =
                            overlaps(locus.locus_local_start_1based, locus.locus_local_end_1based)?;
                        rows.push((
                            lane.lane.lane_id.clone(),
                            ids,
                            available,
                            RegulatoryFragmentExternalEvidence::Occupancy(lane.clone()),
                        ));
                    }
                }
                _ => {
                    return Err(invalid(
                        "Content resolution supports only the three external evidence dimensions"
                            .into(),
                    ));
                }
            }
            if rows.is_empty() {
                if binding.row_id.is_some() {
                    return Err(invalid(
                        "Requested evidence row is absent from the typed report".into(),
                    ));
                }
                all_assessed = false;
            }
            rows.sort_by(|a, b| a.0.cmp(&b.0));
            for (row_id, fragment_ids, available, evidence) in rows {
                all_assessed &= available;
                if observations.len() == request.policy.max_evidence_observations_per_dimension {
                    truncated = true;
                    continue;
                }
                observations.push(
                    RegulatoryFragmentEvidenceObservation::ExternalLocusContext {
                        report_id: binding.report_id.clone(),
                        report_sha256: binding.report_sha256.clone(),
                        row_id,
                        fragment_ids,
                        state: if available {
                            RegulatoryFragmentEvidenceState::Evaluated
                        } else {
                            RegulatoryFragmentEvidenceState::Unavailable
                        },
                        evidence: Box::new(evidence),
                    },
                );
            }
        }
        let state = if !observations.is_empty() && all_assessed {
            RegulatoryFragmentEvidenceState::Evaluated
        } else {
            RegulatoryFragmentEvidenceState::NotEvaluated
        };
        Self::regulatory_fragment_evidence_dimension(
            request,
            kind,
            state,
            if observations.is_empty() {
                vec![]
            } else {
                vec!["hash_bound_locus_document_v1".into()]
            },
            observations,
            truncated,
            "Original typed locus evidence, bound to current sequence/anchor, assembly and gene-annotation release. Row geometry and source provenance remain authoritative; fragment_ids identify the assessed local context, not activity, causality or sufficiency. Missing or citation-only inputs are not evaluated, never a pass.",
        )
    }
}
