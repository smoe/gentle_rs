//! Exact reference admission and companion selections for specificity gates.

use super::*;

const REFERENCE_SELECTION_SCHEMA: &str = "gentle.primer_specificity_reference_selection.v1";

impl GentleEngine {
    pub(super) fn primer_specificity_reference_identity(
        database: &crate::genomes::BlastDatabaseInspectionReport,
    ) -> Option<PrimerSpecificityReferenceIdentity> {
        let fingerprint = database.content_fingerprint.as_ref()?.trim();
        if database.validation_status != "valid"
            || fingerprint.is_empty()
            || database.fingerprint_algorithm.trim().is_empty()
            || database.subject_annotation_status == "stale"
            || (database.subject_annotation_fingerprint.is_some()
                && database.subject_annotation_status != "ready")
        {
            return None;
        }
        Some(PrimerSpecificityReferenceIdentity {
            genome_id: database.source_genome_id.clone(),
            index_kind: database.index_kind,
            assembly: database.source_assembly.clone(),
            release: database.source_release.clone(),
            content_fingerprint: fingerprint.to_string(),
            fingerprint_algorithm: database.fingerprint_algorithm.clone(),
            subject_annotation_fingerprint: database.subject_annotation_fingerprint.clone(),
            subject_annotation_fingerprint_algorithm: database
                .subject_annotation_fingerprint_algorithm
                .clone(),
        })
    }

    pub(super) fn primer_specificity_source_reference_from_anchor(
        &self,
        anchor: &SequenceGenomeAnchorSummary,
    ) -> Option<PrimerSpecificitySourceReference> {
        if anchor.anchor_verified != Some(true) {
            return None;
        }
        let provenance = self.latest_genome_extraction_provenance_for_seq(&anchor.seq_id)?;
        if provenance.genome_id != anchor.genome_id {
            return None;
        }
        let (catalog, _) =
            Self::open_reference_genome_catalog(Some(&provenance.catalog_path)).ok()?;
        let (genome_id, entry) = catalog.exact_catalog_entry(&anchor.genome_id).ok()?;
        Some(PrimerSpecificitySourceReference {
            genome_id,
            assembly: entry
                .reference_name
                .clone()
                .or_else(|| entry.ncbi_assembly_accession.clone())
                .or_else(|| entry.ncbi_assembly_name.clone()),
            release: entry.reference_release.clone().or_else(|| {
                entry
                    .ensembl_template
                    .as_ref()
                    .map(|template| format!("{} release {}", template.provider, template.release))
            }),
            sequence_sha1: provenance.sequence_sha1,
            annotation_sha1: provenance.annotation_sha1,
        })
    }

    /// Admit coordinates before subject alias normalization. Unknown geometry
    /// must not turn matches on another assembly into intended products.
    pub(super) fn primer_specificity_bind_intended_target(
        target: &PrimerSpecificityIntendedTarget,
        database: Option<&crate::genomes::BlastDatabaseInspectionReport>,
    ) -> PrimerSpecificityIntendedTarget {
        let identity = database.and_then(Self::primer_specificity_reference_identity);
        let explicit_match = target
            .reference_binding
            .as_ref()
            .is_some_and(|binding| identity.as_ref() == Some(binding));
        let source_match = target.source_reference.as_ref().is_some_and(|source| {
            identity.as_ref().is_some_and(|reference| {
                source.genome_id == reference.genome_id
                    && source
                        .assembly
                        .as_ref()
                        .is_some_and(|value| !value.trim().is_empty())
                    && source
                        .release
                        .as_ref()
                        .is_some_and(|value| !value.trim().is_empty())
                    && source.assembly == reference.assembly
                    && source.release == reference.release
            })
        });
        if explicit_match || (target.reference_binding.is_none() && source_match) {
            let mut bound = target.clone();
            bound.reference_binding = identity;
            return bound;
        }
        let mut unknown = PrimerSpecificityIntendedTarget {
            source_reference: target.source_reference.clone(),
            source: target.source.clone(),
            warnings: target.warnings.clone(),
            ..PrimerSpecificityIntendedTarget::default()
        };
        unknown.warnings.push(
            "Intended-target reference binding is missing or incompatible; coordinates and transcript identifiers were not transferred to the searched reference. Status is not_assessed, not biological absence."
                .to_string(),
        );
        unknown
    }

    pub(super) fn primer_specificity_pair_binding(
        forward: &PrimerSpecificityInputPrimer,
        reverse: &PrimerSpecificityInputPrimer,
    ) -> Result<String, EngineError> {
        let full_pair = Self::primer_specificity_pair_content_sha256(
            &forward.full_sequence,
            &reverse.full_sequence,
        )?;
        serde_json::to_vec(&json!({
            "full_pair": full_pair, "forward": forward, "reverse": reverse,
        }))
        .map(|bytes| sha256_prefixed_bytes(&bytes))
        .map_err(|error| EngineError::internal(format!("Could not bind primer pair: {error}")))
    }

    fn primer_specificity_selection_key(panel_id: &str, kind: BlastDatabaseIndexKind) -> String {
        format!("{panel_id}:{}", kind.as_str())
    }

    pub(super) fn primer_specificity_selection_panel_digest(
        panel: &TranscriptAssayPanelReport,
    ) -> Result<String, EngineError> {
        // The communication summary is refreshed by QA, not a design edit.
        let mut design = panel.clone();
        for assay in &mut design.selected_assays {
            assay.primer_pair_summary = PrimerPairCommunicationSummary::default();
        }
        Self::transcript_assay_panel_specificity_digest(&design)
    }

    pub(super) fn record_primer_specificity_reference_selection(
        store: &mut PrimerDesignStore,
        panel: &TranscriptAssayPanelReport,
        acceptance: &TranscriptAssayPanelSpecificityAcceptance,
    ) -> Result<Option<String>, EngineError> {
        if !matches!(
            acceptance.status,
            TranscriptAssayPanelSpecificityAcceptanceStatus::Pass
                | TranscriptAssayPanelSpecificityAcceptanceStatus::SpecificityFail
        ) || acceptance.assessments.len() != panel.selected_assays.len()
        {
            return Ok(None);
        }
        let mut reference = None;
        let mut report_ids = BTreeMap::new();
        let mut report_hashes = BTreeMap::new();
        let mut pairs = BTreeMap::new();
        for assessment in &acceptance.assessments {
            let report = &assessment.report;
            let Some(identity) = report
                .blast_database
                .as_ref()
                .and_then(Self::primer_specificity_reference_identity)
            else {
                return Ok(None);
            };
            if reference
                .as_ref()
                .is_some_and(|expected| expected != &identity)
                || report.intended_target.reference_binding.as_ref() != Some(&identity)
                || !report.search_completeness.complete
                || report.raw_detail_artifacts.len() != 2
                || report
                    .raw_detail_artifacts
                    .iter()
                    .any(|artifact| artifact.checksum.is_none())
            {
                return Ok(None);
            }
            let Some(assay) = panel
                .selected_assays
                .iter()
                .find(|assay| assay.assay_id == assessment.assay_id)
            else {
                return Ok(None);
            };
            let forward = Self::primer_specificity_input_from_record(
                PrimerSpecificityPrimerRole::Forward,
                &assay.primer_pair.forward,
            )?;
            let reverse = Self::primer_specificity_input_from_record(
                PrimerSpecificityPrimerRole::Reverse,
                &assay.primer_pair.reverse,
            )?;
            let expected_pair = Self::primer_specificity_pair_binding(&forward, &reverse)?;
            let Some(report_pair) = Self::primer_specificity_report_pair_binding(report)? else {
                return Ok(None);
            };
            if expected_pair != report_pair {
                return Ok(None);
            }
            report_ids.insert(assay.assay_id.clone(), report.report_id.clone());
            report_hashes.insert(
                assay.assay_id.clone(),
                sha256_prefixed_bytes(
                    &serde_json::to_vec(report)
                        .map_err(|error| EngineError::internal(error.to_string()))?,
                ),
            );
            pairs.insert(assay.assay_id.clone(), expected_pair);
            reference = Some(identity);
        }
        let Some(reference) = reference else {
            return Ok(None);
        };
        let key = Self::primer_specificity_selection_key(&panel.report_id, reference.index_kind);
        let previous = store
            .active_primer_specificity_reference_selections
            .get(&key)
            .cloned();
        let panel_digest = Self::primer_specificity_selection_panel_digest(panel)?;
        let selection_id = short_sha256_id(
            "specificity_selection",
            &serde_json::to_string(&json!({
                "panel_digest": panel_digest, "reference": reference,
                "policy": acceptance.policy, "reports": report_ids, "pairs": pairs,
                "report_hashes": report_hashes,
                "acceptance_id": acceptance.acceptance_id,
            }))
            .map_err(|error| EngineError::internal(error.to_string()))?,
        );
        let selection = PrimerSpecificityReferenceSelection {
            schema: REFERENCE_SELECTION_SCHEMA.to_string(),
            selection_id: selection_id.clone(),
            panel_report_id: panel.report_id.clone(),
            panel_digest,
            reference,
            policy: acceptance.policy.clone(),
            report_ids_by_assay: report_ids,
            report_content_sha256_by_assay: report_hashes,
            pair_bindings_by_assay: pairs,
            acceptance: acceptance.clone(),
            replaces_selection_id: previous.clone().filter(|id| id != &selection_id),
        };
        for assessment in &acceptance.assessments {
            store.primer_specificity_reports.insert(
                assessment.report.report_id.clone(),
                assessment.report.clone(),
            );
        }
        store
            .primer_specificity_reference_selections
            .insert(selection_id.clone(), selection);
        store
            .active_primer_specificity_reference_selections
            .insert(key, selection_id.clone());
        Ok(previous.filter(|id| id != &selection_id))
    }

    pub(super) fn primer_specificity_report_pair_binding(
        report: &PrimerSpecificityReport,
    ) -> Result<Option<String>, EngineError> {
        let forward = report
            .primers
            .iter()
            .find(|primer| primer.role == PrimerSpecificityPrimerRole::Forward);
        let reverse = report
            .primers
            .iter()
            .find(|primer| primer.role == PrimerSpecificityPrimerRole::Reverse);
        match (forward, reverse) {
            (Some(forward), Some(reverse)) if report.primers.len() == 2 => {
                Self::primer_specificity_pair_binding(forward, reverse).map(Some)
            }
            _ => Ok(None),
        }
    }

    pub(super) fn selected_primer_specificity_report(
        &self,
        panel: &TranscriptAssayPanelReport,
        assay_id: &str,
        kind: BlastDatabaseIndexKind,
    ) -> Option<PrimerSpecificityReport> {
        let store = self.read_primer_design_store();
        let key = Self::primer_specificity_selection_key(&panel.report_id, kind);
        let selection_id = store
            .active_primer_specificity_reference_selections
            .get(&key)?;
        let selection = store
            .primer_specificity_reference_selections
            .get(selection_id)?;
        if selection.schema != REFERENCE_SELECTION_SCHEMA
            || selection.panel_digest
                != Self::primer_specificity_selection_panel_digest(panel).ok()?
        {
            return None;
        }
        let report_id = selection.report_ids_by_assay.get(assay_id)?;
        let report = &selection
            .acceptance
            .assessments
            .iter()
            .find(|row| row.assay_id == assay_id && &row.report.report_id == report_id)?
            .report;
        let report_hash = sha256_prefixed_bytes(&serde_json::to_vec(report).ok()?);
        if selection.report_content_sha256_by_assay.get(assay_id) != Some(&report_hash) {
            return None;
        }
        let reference = report
            .blast_database
            .as_ref()
            .and_then(Self::primer_specificity_reference_identity)?;
        if reference != selection.reference
            || report.target_kind != kind.as_str()
            || report.intended_target.reference_binding.as_ref() != Some(&reference)
            || serde_json::to_value(&report.policy).ok()?
                != serde_json::to_value(&selection.policy).ok()?
        {
            return None;
        }
        let assay = panel
            .selected_assays
            .iter()
            .find(|assay| assay.assay_id == assay_id)?;
        let forward = Self::primer_specificity_input_from_record(
            PrimerSpecificityPrimerRole::Forward,
            &assay.primer_pair.forward,
        )
        .ok()?;
        let reverse = Self::primer_specificity_input_from_record(
            PrimerSpecificityPrimerRole::Reverse,
            &assay.primer_pair.reverse,
        )
        .ok()?;
        let expected_pair = Self::primer_specificity_pair_binding(&forward, &reverse).ok()?;
        if selection.pair_bindings_by_assay.get(assay_id) != Some(&expected_pair)
            || Self::primer_specificity_report_pair_binding(report).ok()?? != expected_pair
        {
            return None;
        }
        Some(report.clone())
    }

    /// Explicit handoff construction rechecks resources; display never does.
    pub(super) fn primer_specificity_selection_freshness(
        &self,
        panel: &TranscriptAssayPanelReport,
    ) -> BTreeMap<String, bool> {
        let mut freshness = BTreeMap::new();
        let Some(assay) = panel.selected_assays.first() else {
            return freshness;
        };
        for kind in [
            BlastDatabaseIndexKind::GenomicDna,
            BlastDatabaseIndexKind::TranscriptomeCdna,
        ] {
            let Some(report) =
                self.selected_primer_specificity_report(panel, &assay.assay_id, kind)
            else {
                continue;
            };
            let expected = report
                .blast_database
                .as_ref()
                .and_then(Self::primer_specificity_reference_identity);
            let current = Self::open_reference_genome_catalog(report.catalog_path.as_deref())
                .ok()
                .and_then(|(catalog, _)| {
                    catalog
                        .inspect_blast_database(
                            &report.target_genome_id,
                            report.cache_dir.as_deref(),
                        )
                        .ok()
                        .flatten()
                })
                .as_ref()
                .and_then(Self::primer_specificity_reference_identity);
            freshness.insert(
                kind.as_str().to_string(),
                expected.is_some() && expected == current,
            );
        }
        freshness
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Hand-crafted reference/pair records, not genes or production evidence.
    fn database(
        id: &str,
        kind: BlastDatabaseIndexKind,
    ) -> crate::genomes::BlastDatabaseInspectionReport {
        crate::genomes::BlastDatabaseInspectionReport {
            source_genome_id: id.to_string(),
            index_kind: kind,
            source_assembly: Some("synthetic-assembly".to_string()),
            source_release: Some("synthetic-release-1".to_string()),
            content_fingerprint: Some(format!("sha256:{id}")),
            fingerprint_algorithm: "synthetic-sha256".to_string(),
            validation_status: "valid".to_string(),
            ..Default::default()
        }
    }

    fn panel() -> TranscriptAssayPanelReport {
        let primer = |sequence: &str| PrimerDesignPrimerRecord {
            sequence: sequence.to_string(),
            length_bp: sequence.len(),
            anneal_length_bp: sequence.len(),
            ..Default::default()
        };
        TranscriptAssayPanelReport {
            report_id: "synthetic_panel".to_string(),
            selected_assays: vec![TranscriptAssayPanelAssay {
                assay_id: "synthetic_assay".to_string(),
                rank: 1,
                primer_pair: PrimerDesignPairRecord {
                    forward: primer("ACGTACGT"),
                    reverse: primer("TGCATGCA"),
                    ..Default::default()
                },
                ..Default::default()
            }],
            ..Default::default()
        }
    }

    fn acceptance(
        panel: &TranscriptAssayPanelReport,
        id: &str,
        kind: BlastDatabaseIndexKind,
        pass: bool,
    ) -> TranscriptAssayPanelSpecificityAcceptance {
        let database = database(id, kind);
        let reference = GentleEngine::primer_specificity_reference_identity(&database).unwrap();
        let assay = &panel.selected_assays[0];
        let forward = GentleEngine::primer_specificity_input_from_record(
            PrimerSpecificityPrimerRole::Forward,
            &assay.primer_pair.forward,
        )
        .unwrap();
        let reverse = GentleEngine::primer_specificity_input_from_record(
            PrimerSpecificityPrimerRole::Reverse,
            &assay.primer_pair.reverse,
        )
        .unwrap();
        let dimension = PrimerSpecificityTargetAssessment {
            status: if pass { "pass" } else { "fail" }.to_string(),
            ..Default::default()
        };
        let report = PrimerSpecificityReport {
            report_id: format!("synthetic_{id}"),
            target_kind: kind.as_str().to_string(),
            target_genome_id: id.to_string(),
            blast_database: Some(database),
            intended_target: PrimerSpecificityIntendedTarget {
                reference_binding: Some(reference),
                model: PrimerSpecificityIntendedTargetModel::GenomicInterval,
                ..Default::default()
            },
            primers: vec![forward, reverse],
            search_completeness: PrimerSpecificitySearchCompleteness {
                complete: true,
                ..Default::default()
            },
            genomic_specificity: dimension.clone(),
            transcriptome_specificity: dimension,
            raw_detail_artifacts: vec![
                ComputationalArtifactExternalInput {
                    checksum: Some("sha256:forward".to_string()),
                    ..Default::default()
                },
                ComputationalArtifactExternalInput {
                    checksum: Some("sha256:reverse".to_string()),
                    ..Default::default()
                },
            ],
            ..Default::default()
        };
        TranscriptAssayPanelSpecificityAcceptance {
            acceptance_id: format!("synthetic_acceptance_{id}"),
            status: if pass {
                TranscriptAssayPanelSpecificityAcceptanceStatus::Pass
            } else {
                TranscriptAssayPanelSpecificityAcceptanceStatus::SpecificityFail
            },
            assessments: vec![TranscriptAssayGenomicSpecificityAssessment {
                assay_id: assay.assay_id.clone(),
                assay_rank: 1,
                report,
                ..Default::default()
            }],
            ..Default::default()
        }
    }

    #[test]
    fn specificity_reference_identity_rejects_stale_annotation_with_planned_hash() {
        let mut database = database("reference_a", BlastDatabaseIndexKind::GenomicDna);
        database.subject_annotation_fingerprint = Some("sha256:planned-annotation".into());
        database.subject_annotation_fingerprint_algorithm = Some("sha256".into());
        database.subject_annotation_status = "ready".into();
        assert!(GentleEngine::primer_specificity_reference_identity(&database).is_some());
        for status in ["stale", "unavailable", "invalid"] {
            database.subject_annotation_status = status.into();
            assert!(GentleEngine::primer_specificity_reference_identity(&database).is_none());
        }
    }

    #[test]
    fn specificity_reference_binding_rejects_other_assembly_before_aliasing() {
        let database = database("reference_a", BlastDatabaseIndexKind::GenomicDna);
        let target = PrimerSpecificityIntendedTarget {
            model: PrimerSpecificityIntendedTargetModel::GenomicInterval,
            subject_id: Some("chr1".to_string()),
            source_reference: Some(PrimerSpecificitySourceReference {
                genome_id: "reference_a".to_string(),
                assembly: Some("different-assembly".to_string()),
                release: database.source_release.clone(),
                ..Default::default()
            }),
            ..Default::default()
        };
        let refused =
            GentleEngine::primer_specificity_bind_intended_target(&target, Some(&database));
        assert_eq!(refused.model, PrimerSpecificityIntendedTargetModel::Unknown);
        assert!(refused.subject_id.is_none());
        let mut explicit = target;
        explicit.reference_binding = GentleEngine::primer_specificity_reference_identity(&database);
        explicit.source = "caller_provided_intended_target".to_string();
        let admitted =
            GentleEngine::primer_specificity_bind_intended_target(&explicit, Some(&database));
        assert_eq!(admitted.subject_id.as_deref(), Some("chr1"));
        let mut changed = database;
        changed.content_fingerprint = Some("sha256:replacement".to_string());
        assert_eq!(
            GentleEngine::primer_specificity_bind_intended_target(&explicit, Some(&changed)).model,
            PrimerSpecificityIntendedTargetModel::Unknown
        );
        explicit.reference_binding = None;
        explicit.source_reference = None;
        assert!(
            GentleEngine::primer_specificity_bind_intended_target(&explicit, Some(&changed))
                .reference_binding
                .is_none()
        );
    }

    #[test]
    fn specificity_reference_selection_is_pinned_and_replacement_retains_history() {
        let mut engine = GentleEngine::default();
        let mut panel = panel();
        let failed = acceptance(
            &panel,
            "reference_a",
            BlastDatabaseIndexKind::GenomicDna,
            false,
        );
        let mut newer = acceptance(
            &panel,
            "reference_b",
            BlastDatabaseIndexKind::GenomicDna,
            true,
        );
        newer.assessments[0].report.generated_at_unix_ms = 100;
        let mut store = PrimerDesignStore::default();
        GentleEngine::record_primer_specificity_reference_selection(&mut store, &panel, &failed)
            .unwrap();
        panel.genomic_specificity_assessments = newer.assessments.clone();
        engine.write_primer_design_store(store).unwrap();
        assert_eq!(
            engine
                .selected_primer_specificity_report(
                    &panel,
                    "synthetic_assay",
                    BlastDatabaseIndexKind::GenomicDna
                )
                .unwrap()
                .genomic_specificity
                .status,
            "fail"
        );
        let mut store = engine.read_primer_design_store();
        assert!(
            GentleEngine::record_primer_specificity_reference_selection(&mut store, &panel, &newer)
                .unwrap()
                .is_some()
        );
        assert_eq!(store.primer_specificity_reference_selections.len(), 2);
        let transcript = acceptance(
            &panel,
            "transcript_ref",
            BlastDatabaseIndexKind::TranscriptomeCdna,
            true,
        );
        GentleEngine::record_primer_specificity_reference_selection(
            &mut store,
            &panel,
            &transcript,
        )
        .unwrap();
        let before = store.active_primer_specificity_reference_selections.clone();
        newer.status = TranscriptAssayPanelSpecificityAcceptanceStatus::Incomplete;
        GentleEngine::record_primer_specificity_reference_selection(&mut store, &panel, &newer)
            .unwrap();
        assert_eq!(before, store.active_primer_specificity_reference_selections);
        engine.write_primer_design_store(store).unwrap();
        let state = serde_json::from_slice(&serde_json::to_vec(engine.state()).unwrap()).unwrap();
        let restored = GentleEngine::from_state(state);
        for kind in [
            BlastDatabaseIndexKind::GenomicDna,
            BlastDatabaseIndexKind::TranscriptomeCdna,
        ] {
            assert!(
                restored
                    .selected_primer_specificity_report(&panel, "synthetic_assay", kind)
                    .is_some()
            );
        }
        panel.selected_assays[0].primer_pair.forward.sequence = "ACGTACGA".to_string();
        assert!(
            restored
                .selected_primer_specificity_report(
                    &panel,
                    "synthetic_assay",
                    BlastDatabaseIndexKind::GenomicDna
                )
                .is_none()
        );
    }

    #[test]
    fn specificity_reference_legacy_and_changed_tail_do_not_claim_pass() {
        let engine = GentleEngine::default();
        assert!(
            engine
                .selected_primer_specificity_report(
                    &panel(),
                    "synthetic_assay",
                    BlastDatabaseIndexKind::GenomicDna
                )
                .is_none()
        );
        let forward = GentleEngine::primer_specificity_input_from_record(
            PrimerSpecificityPrimerRole::Forward,
            &panel().selected_assays[0].primer_pair.forward,
        )
        .unwrap();
        let mut tailed = forward.clone();
        tailed.annealing_sequence.remove(0);
        tailed.non_annealing_5prime_tail = "A".to_string();
        tailed.non_annealing_5prime_tail_bp = 1;
        tailed.annealing_length_bp -= 1;
        assert_ne!(
            GentleEngine::primer_specificity_pair_binding(&forward, &forward).unwrap(),
            GentleEngine::primer_specificity_pair_binding(&tailed, &forward).unwrap()
        );
    }
}
