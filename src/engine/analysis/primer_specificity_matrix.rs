//! Read-only panel specificity history, preserving receipt scopes and readiness.

use super::*;

impl GentleEngine {
    /// Inspect retained receipt snapshots without resource probes or store writes.
    pub fn inspect_transcript_assay_panel_specificity_matrix(
        &self,
        panel_report_id: &str,
    ) -> Result<TranscriptAssayPanelSpecificityMatrix, EngineError> {
        let panel = self.get_transcript_assay_panel_report(panel_report_id)?;
        let design_digest = Self::primer_specificity_selection_panel_digest(&panel)?;
        let store = self.read_primer_design_store();
        let mut assays = Vec::new();
        let mut current_pairs = BTreeMap::new();
        for assay in &panel.selected_assays {
            let forward = Self::primer_specificity_input_from_record(
                PrimerSpecificityPrimerRole::Forward,
                &assay.primer_pair.forward,
            )?;
            let reverse = Self::primer_specificity_input_from_record(
                PrimerSpecificityPrimerRole::Reverse,
                &assay.primer_pair.reverse,
            )?;
            let binding = Self::primer_specificity_pair_binding(&forward, &reverse)?;
            if current_pairs
                .insert(assay.assay_id.clone(), binding.clone())
                .is_some()
            {
                return Err(EngineError::invalid_input(
                    "Duplicate selected assay identity",
                ));
            }
            assays.push(TranscriptAssaySpecificityMatrixAssay {
                assay_id: assay.assay_id.clone(),
                assay_rank: assay.rank,
                pair_binding_sha256: binding,
            });
        }
        assays.sort_by(|a, b| {
            a.assay_rank
                .cmp(&b.assay_rank)
                .then(a.assay_id.cmp(&b.assay_id))
        });
        let retained = store
            .primer_specificity_reference_selections
            .iter()
            .filter(|(_, selection)| selection.panel_report_id == panel.report_id)
            .collect::<Vec<_>>();
        let superseded = retained
            .iter()
            .filter_map(|(_, s)| s.replaces_selection_id.clone())
            .collect::<BTreeSet<_>>();
        let mut selections = Vec::new();
        let mut assessed = BTreeSet::new();
        let mut known_pairs = current_pairs.values().cloned().collect::<BTreeSet<_>>();
        for (key, selection) in retained {
            let matches_design = selection.panel_digest == design_digest;
            let mut diagnostics = Vec::new();
            let acceptance = &selection.acceptance;
            let expected_id = short_sha256_id(
                "specificity_selection",
                &serde_json::to_string(&json!({
                    "panel_digest":selection.panel_digest, "reference":selection.reference,
                    "policy":selection.policy, "reports":selection.report_ids_by_assay,
                    "pairs":selection.pair_bindings_by_assay,
                    "report_hashes":selection.report_content_sha256_by_assay,
                    "acceptance_id":acceptance.acceptance_id,
                }))
                .map_err(|e| EngineError::internal(e.to_string()))?,
            );
            if selection.schema != "gentle.primer_specificity_reference_selection.v1"
                || key != &selection.selection_id
                || selection.selection_id != expected_id
                || selection.panel_digest.is_empty()
            {
                diagnostics.push("selection_identity_invalid".into());
            }
            if acceptance.panel_report_id != panel.report_id
                || !matches!(
                    acceptance.status,
                    TranscriptAssayPanelSpecificityAcceptanceStatus::Pass
                        | TranscriptAssayPanelSpecificityAcceptanceStatus::SpecificityFail
                )
                || acceptance.accepted
                    != (acceptance.status == TranscriptAssayPanelSpecificityAcceptanceStatus::Pass)
                || serde_json::to_value(&acceptance.policy).ok()
                    != serde_json::to_value(&selection.policy).ok()
                || acceptance.expected_assay_count != acceptance.assessments.len()
                || acceptance.assessed_assay_count != acceptance.assessments.len()
                || (acceptance.status == TranscriptAssayPanelSpecificityAcceptanceStatus::Pass
                    && acceptance
                        .assessments
                        .iter()
                        .any(|a| a.report.summary.status != "pass"))
                || (acceptance.status
                    == TranscriptAssayPanelSpecificityAcceptanceStatus::SpecificityFail
                    && !acceptance
                        .assessments
                        .iter()
                        .any(|a| a.report.summary.status == "fail"))
            {
                diagnostics.push("selection_acceptance_incomplete_or_inconsistent".into());
            }
            let recorded_ids = acceptance
                .assessments
                .iter()
                .map(|a| a.assay_id.clone())
                .collect::<BTreeSet<_>>();
            if recorded_ids.len() != acceptance.assessments.len()
                || recorded_ids != selection.report_ids_by_assay.keys().cloned().collect()
                || recorded_ids
                    != selection
                        .report_content_sha256_by_assay
                        .keys()
                        .cloned()
                        .collect()
                || recorded_ids != selection.pair_bindings_by_assay.keys().cloned().collect()
            {
                diagnostics.push("selection_assay_inventory_inconsistent".into());
            }
            let mut ids = recorded_ids;
            ids.extend(current_pairs.keys().cloned());
            let mut cells = Vec::new();
            for assay_id in ids {
                let assessment = acceptance
                    .assessments
                    .iter()
                    .find(|a| a.assay_id == assay_id);
                let mut issues = diagnostics.clone();
                let report_id = selection.report_ids_by_assay.get(&assay_id).cloned();
                let report_hash = selection
                    .report_content_sha256_by_assay
                    .get(&assay_id)
                    .cloned();
                let mut recorded_status = None;
                let mut verdict = PrimerSpecificityMultiVerdict::Incomplete;
                let pair = selection.pair_bindings_by_assay.get(&assay_id);
                let matches_assay =
                    matches_design && pair.is_some() && pair == current_pairs.get(&assay_id);
                if let Some(assessment) = assessment {
                    let report = &assessment.report;
                    recorded_status = Some(report.summary.status.clone());
                    let dimension = match selection.reference.index_kind {
                        BlastDatabaseIndexKind::GenomicDna => &report.genomic_specificity,
                        BlastDatabaseIndexKind::TranscriptomeCdna => {
                            &report.transcriptome_specificity
                        }
                    };
                    let expected_status = match report.summary.status.as_str() {
                        "pass" => "external_blast_pass",
                        "fail" => "external_blast_specificity_fail",
                        "not_assessed" => "external_blast_not_assessed",
                        _ => "external_blast_incomplete",
                    };
                    if dimension.status != report.summary.status
                        || assessment.status != expected_status
                        || (report.summary.status == "pass" && !report.summary.specificity_pass)
                        || (report.summary.status == "fail" && report.summary.specificity_pass)
                    {
                        issues.push("retained_verdict_inconsistent".into());
                    }
                    let hash = sha256_prefixed_bytes(
                        &serde_json::to_vec(report)
                            .map_err(|e| EngineError::internal(e.to_string()))?,
                    );
                    if report_id.as_deref() != Some(report.report_id.as_str())
                        || report_hash.as_ref() != Some(&hash)
                    {
                        issues.push("retained_report_hash_or_identity_mismatch".into());
                    }
                    if report
                        .blast_database
                        .as_ref()
                        .and_then(Self::primer_specificity_reference_identity)
                        .as_ref()
                        != Some(&selection.reference)
                        || report.intended_target.reference_binding.as_ref()
                            != Some(&selection.reference)
                        || report.target_kind != selection.reference.index_kind.as_str()
                        || report.target_genome_id != selection.reference.genome_id
                        || serde_json::to_value(&report.policy).ok()
                            != serde_json::to_value(&selection.policy).ok()
                    {
                        issues.push("retained_reference_or_policy_mismatch".into());
                    }
                    let report_pair = Self::primer_specificity_report_pair_binding(report)?;
                    if pair.is_none()
                        || report_pair.as_ref() != pair
                        || (matches_design && !matches_assay)
                    {
                        issues.push("retained_pair_or_current_assay_mismatch".into());
                    }
                    if !report.search_completeness.complete
                        || report.raw_detail_artifacts.len() != 2
                        || report
                            .raw_detail_artifacts
                            .iter()
                            .any(|a| a.checksum.as_deref().is_none_or(str::is_empty))
                    {
                        issues.push("retained_search_evidence_incomplete".into());
                    }
                    if issues.is_empty() {
                        verdict = match report.summary.status.as_str() {
                            "pass" if report.summary.specificity_pass => {
                                PrimerSpecificityMultiVerdict::Pass
                            }
                            "fail" => PrimerSpecificityMultiVerdict::Fail,
                            _ => PrimerSpecificityMultiVerdict::Incomplete,
                        };
                        if let Some(pair) = pair {
                            known_pairs.insert(pair.clone());
                        }
                    }
                } else {
                    issues.push("assay_not_in_retained_selection".into());
                }
                let verified = issues.is_empty();
                if verified
                    && matches_assay
                    && matches!(
                        verdict,
                        PrimerSpecificityMultiVerdict::Pass | PrimerSpecificityMultiVerdict::Fail
                    )
                {
                    assessed.insert((
                        selection.reference.index_kind.as_str().to_string(),
                        assay_id.clone(),
                    ));
                }
                cells.push(TranscriptAssaySpecificityMatrixCell {
                    assay_id,
                    report_id,
                    report_content_sha256: report_hash,
                    recorded_status,
                    verdict,
                    retained_binding_verified: verified,
                    matches_current_assay: matches_assay,
                    diagnostics: issues,
                });
            }
            let active_key = format!(
                "{}:{}",
                panel.report_id,
                selection.reference.index_kind.as_str()
            );
            selections.push(TranscriptAssaySpecificityMatrixSelection {
                selection_id: selection.selection_id.clone(),
                panel_design_digest: selection.panel_digest.clone(),
                reference: selection.reference.clone(),
                policy: selection.policy.clone(),
                acceptance_id: acceptance.acceptance_id.clone(),
                recorded_acceptance_status: acceptance.status,
                is_active: store
                    .active_primer_specificity_reference_selections
                    .get(&active_key)
                    == Some(key),
                is_superseded: superseded.contains(&selection.selection_id),
                replaces_selection_id: selection.replaces_selection_id.clone(),
                matches_current_panel_design: matches_design,
                cells,
                diagnostics,
            });
        }
        let missing = |kind: BlastDatabaseIndexKind| {
            assays
                .iter()
                .filter(|a| !assessed.contains(&(kind.as_str().to_string(), a.assay_id.clone())))
                .map(|a| a.assay_id.clone())
                .collect()
        };
        let unlinked_pair_summaries = store.primer_specificity_multi_summaries.values()
            .filter(|s| known_pairs.contains(&s.pair_binding_sha256))
            .map(|s| TranscriptAssaySpecificityMatrixUnlinkedSummary {
                summary_id:s.summary_id.clone(), pair_binding_sha256:s.pair_binding_sha256.clone(),
                same_oligo_current_assay_ids:assays.iter().filter(|a| a.pair_binding_sha256 == s.pair_binding_sha256)
                    .map(|a| a.assay_id.clone()).collect(),
                reason:"Same full oligos and tail interpretation only; no verified panel/assay/context association. Not counted as panel evidence.".into(),
            }).collect();
        Ok(TranscriptAssayPanelSpecificityMatrix {
            schema:"gentle.transcript_assay_panel_specificity_matrix.v1".into(),
            panel_report_id:panel.report_id, current_panel_design_digest:design_digest,
            unassessed_genomic_assay_ids:missing(BlastDatabaseIndexKind::GenomicDna),
            unassessed_transcriptome_assay_ids:missing(BlastDatabaseIndexKind::TranscriptomeCdna),
            assays, selections, unlinked_pair_summaries,
            nonclaims:vec![
                "Retained binding verification is not authentication of external execution or current resource currency; no databases were probed.".into(),
                "Matching the saved panel design does not revalidate present source DNA or reference currency.".into(),
                "Each receipt keeps its own reference, policy and panel-design scope. No aggregate panel pass or required/optional status is inferred.".into(),
                "Active selection metadata is displayed only; this inspection does not change readiness, establish isoform discrimination or authorize ordering.".into(),
            ],
        })
    }
}

#[cfg(test)]
#[path = "primer_specificity_matrix_tests.rs"]
mod tests;
