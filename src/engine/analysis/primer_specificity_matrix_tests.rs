//! Inline synthetic receipt fixtures for panel inspection; no genes or resources.
//! Recreated by `fixture`/`acceptance`, used only by these shared-engine tests.
//! Stored IDs follow production normalization; adapter requests retain punctuation.

use super::*;

const RAW_PANEL_REPORT_ID: &str = "synthetic panel ' inspection";

fn fixture() -> (GentleEngine, TranscriptAssayPanelReport) {
    let primer = |sequence: &str| PrimerDesignPrimerRecord {
        sequence: sequence.into(),
        length_bp: sequence.len(),
        anneal_length_bp: sequence.len(),
        ..Default::default()
    };
    // Identical oligos in two different assay contexts must remain two identities.
    let panel = TranscriptAssayPanelReport {
        report_id: GentleEngine::normalize_primer_design_report_id(RAW_PANEL_REPORT_ID).unwrap(),
        source_seq_id: "synthetic_source".into(),
        selected_assays: ["assay_b", "assay_a"]
            .into_iter()
            .enumerate()
            .map(|(i, id)| TranscriptAssayPanelAssay {
                assay_id: id.into(),
                rank: 2 - i,
                design_transcript_id: format!("synthetic_tx_{id}"),
                primer_pair: PrimerDesignPairRecord {
                    forward: primer("ACGTACGT"),
                    reverse: primer("TGCATGCA"),
                    ..Default::default()
                },
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    };
    let mut engine = GentleEngine::default();
    let mut store = PrimerDesignStore::default();
    store
        .transcript_assay_panels
        .insert(panel.report_id.clone(), panel.clone());
    engine.write_primer_design_store(store).unwrap();
    (engine, panel)
}

#[test]
fn specificity_matrix_raw_report_id_resolves_canonical_fixture_without_mutation() {
    let (engine, panel) = fixture();
    assert_eq!(panel.report_id, "synthetic_panel___inspection");
    assert_ne!(panel.report_id, RAW_PANEL_REPORT_ID);
    let before = serde_json::to_value(engine.snapshot()).unwrap();
    for requested_id in [RAW_PANEL_REPORT_ID, panel.report_id.as_str()] {
        assert_eq!(
            engine
                .get_transcript_assay_panel_report(requested_id)
                .unwrap()
                .report_id,
            panel.report_id
        );
        let matrix = engine
            .inspect_transcript_assay_panel_specificity_matrix(requested_id)
            .unwrap();
        assert_eq!(matrix.panel_report_id, panel.report_id);
        assert_eq!(matrix.assays.len(), panel.selected_assays.len());
        assert_eq!(matrix.unassessed_genomic_assay_ids.len(), 2);
        assert_eq!(matrix.unassessed_transcriptome_assay_ids.len(), 2);
    }
    assert_eq!(before, serde_json::to_value(engine.snapshot()).unwrap());
}

fn acceptance(
    panel: &TranscriptAssayPanelReport,
    id: &str,
    kind: BlastDatabaseIndexKind,
    fail: bool,
) -> TranscriptAssayPanelSpecificityAcceptance {
    let policy = PrimerSpecificityPolicy::default();
    let assessments = panel
        .selected_assays
        .iter()
        .map(|assay| {
            let pass = !fail || assay.assay_id != "assay_a";
            let status = if pass { "pass" } else { "fail" };
            let database = crate::genomes::BlastDatabaseInspectionReport {
                source_genome_id: id.into(),
                index_kind: kind,
                source_assembly: Some("synthetic_assembly".into()),
                source_release: Some("synthetic_release".into()),
                content_fingerprint: Some(format!("sha256:{id}")),
                fingerprint_algorithm: "synthetic_sha256".into(),
                validation_status: "valid".into(),
                // A nonexistent prefix proves inspection does not require installed resources.
                ..Default::default()
            };
            let reference = GentleEngine::primer_specificity_reference_identity(&database).unwrap();
            TranscriptAssayGenomicSpecificityAssessment {
                assay_id: assay.assay_id.clone(),
                assay_rank: assay.rank,
                status: if pass {
                    "external_blast_pass"
                } else {
                    "external_blast_specificity_fail"
                }
                .into(),
                report: PrimerSpecificityReport {
                    report_id: format!("synthetic_{id}_{}", assay.assay_id),
                    target_genome_id: id.into(),
                    target_kind: kind.as_str().into(),
                    policy: policy.clone(),
                    blast_database: Some(database),
                    intended_target: PrimerSpecificityIntendedTarget {
                        reference_binding: Some(reference),
                        model: PrimerSpecificityIntendedTargetModel::TranscriptSet,
                        ..Default::default()
                    },
                    primers: vec![
                        GentleEngine::primer_specificity_input_from_record(
                            PrimerSpecificityPrimerRole::Forward,
                            &assay.primer_pair.forward,
                        )
                        .unwrap(),
                        GentleEngine::primer_specificity_input_from_record(
                            PrimerSpecificityPrimerRole::Reverse,
                            &assay.primer_pair.reverse,
                        )
                        .unwrap(),
                    ],
                    summary: PrimerSpecificitySummary {
                        status: status.into(),
                        specificity_pass: pass,
                        ..Default::default()
                    },
                    genomic_specificity: PrimerSpecificityTargetAssessment {
                        status: status.into(),
                        ..Default::default()
                    },
                    transcriptome_specificity: PrimerSpecificityTargetAssessment {
                        status: status.into(),
                        ..Default::default()
                    },
                    search_completeness: PrimerSpecificitySearchCompleteness {
                        complete: true,
                        ..Default::default()
                    },
                    raw_detail_artifacts: ["forward", "reverse"]
                        .map(|role| ComputationalArtifactExternalInput {
                            checksum: Some(format!("sha256:synthetic_{role}")),
                            ..Default::default()
                        })
                        .into(),
                    ..Default::default()
                },
            }
        })
        .collect::<Vec<_>>();
    TranscriptAssayPanelSpecificityAcceptance {
        panel_report_id: panel.report_id.clone(),
        panel_digest: GentleEngine::transcript_assay_panel_specificity_digest(panel).unwrap(),
        acceptance_id: format!("synthetic_acceptance_{id}"),
        policy,
        status: if fail {
            TranscriptAssayPanelSpecificityAcceptanceStatus::SpecificityFail
        } else {
            TranscriptAssayPanelSpecificityAcceptanceStatus::Pass
        },
        accepted: !fail,
        expected_assay_count: assessments.len(),
        assessed_assay_count: assessments.len(),
        assessments,
        ..Default::default()
    }
}

fn retain(
    engine: &mut GentleEngine,
    panel: &TranscriptAssayPanelReport,
    acceptance: &TranscriptAssayPanelSpecificityAcceptance,
) {
    let mut store = engine.read_primer_design_store();
    GentleEngine::record_primer_specificity_reference_selection(&mut store, panel, acceptance)
        .unwrap();
    engine.write_primer_design_store(store).unwrap();
}

#[test]
fn specificity_matrix_preserves_failed_history_policies_and_active_dimensions_without_probes() {
    let (mut engine, panel) = fixture();
    let failed = acceptance(
        &panel,
        "reference_a",
        BlastDatabaseIndexKind::GenomicDna,
        true,
    );
    retain(&mut engine, &panel, &failed);
    let mut newer = acceptance(
        &panel,
        "reference_b",
        BlastDatabaseIndexKind::GenomicDna,
        false,
    );
    newer.policy.max_3prime_mismatches += 1;
    for a in &mut newer.assessments {
        a.report.policy = newer.policy.clone();
    }
    retain(&mut engine, &panel, &newer);
    retain(
        &mut engine,
        &panel,
        &acceptance(
            &panel,
            "cdna_reference",
            BlastDatabaseIndexKind::TranscriptomeCdna,
            false,
        ),
    );
    // Later global-store records must not replace receipt-owned history.
    let mut store = engine.read_primer_design_store();
    for report in store.primer_specificity_reports.values_mut() {
        report.summary.status = "not_assessed".into();
        report.summary.specificity_pass = false;
    }
    engine.write_primer_design_store(store).unwrap();
    let before = serde_json::to_value(engine.snapshot()).unwrap();
    let matrix = engine
        .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
        .unwrap();
    assert_eq!(
        matrix
            .assays
            .iter()
            .map(|a| a.assay_id.as_str())
            .collect::<Vec<_>>(),
        ["assay_a", "assay_b"]
    );
    assert_eq!(matrix.selections.len(), 3);
    let a = matrix
        .selections
        .iter()
        .find(|s| s.reference.genome_id == "reference_a")
        .unwrap();
    let b = matrix
        .selections
        .iter()
        .find(|s| s.reference.genome_id == "reference_b")
        .unwrap();
    assert!(!a.is_active && a.is_superseded && b.is_active);
    assert_ne!(
        serde_json::to_value(&a.policy).unwrap(),
        serde_json::to_value(&b.policy).unwrap()
    );
    assert_eq!(a.cells[0].verdict, PrimerSpecificityMultiVerdict::Fail);
    assert_eq!(b.cells[0].verdict, PrimerSpecificityMultiVerdict::Pass);
    assert!(matrix.unassessed_genomic_assay_ids.is_empty());
    assert!(matrix.unassessed_transcriptome_assay_ids.is_empty());
    let json = serde_json::to_value(&matrix).unwrap();
    assert!(json.get("aggregate").is_none() && json.get("order_ready").is_none());
    let mut forged = json.clone();
    forged["aggregate_pass"] = json!(true);
    assert!(serde_json::from_value::<TranscriptAssayPanelSpecificityMatrix>(forged).is_err());
    assert_eq!(
        json,
        serde_json::to_value(
            engine
                .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
                .unwrap()
        )
        .unwrap()
    );
    assert_eq!(before, serde_json::to_value(engine.snapshot()).unwrap());
}

#[test]
fn specificity_matrix_history_does_not_become_evidence_for_an_edited_panel() {
    let (mut engine, mut panel) = fixture();
    retain(
        &mut engine,
        &panel,
        &acceptance(
            &panel,
            "old_reference",
            BlastDatabaseIndexKind::GenomicDna,
            false,
        ),
    );
    panel.selected_assays[0].assay_id = "new_assay".into();
    panel.selected_assays[0].primer_pair.forward.sequence = "ACGTACGA".into();
    let mut store = engine.read_primer_design_store();
    store
        .transcript_assay_panels
        .insert(panel.report_id.clone(), panel.clone());
    engine.write_primer_design_store(store).unwrap();
    let matrix = engine
        .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
        .unwrap();
    let column = &matrix.selections[0];
    assert!(!column.matches_current_panel_design);
    assert!(column.cells.iter().any(|c| c.assay_id == "assay_b"
        && c.retained_binding_verified
        && c.verdict == PrimerSpecificityMultiVerdict::Pass
        && !c.matches_current_assay));
    assert!(column.cells.iter().any(
        |c| c.assay_id == "new_assay" && c.verdict == PrimerSpecificityMultiVerdict::Incomplete
    ));
    assert_eq!(matrix.unassessed_genomic_assay_ids.len(), 2);
    assert_eq!(matrix.unassessed_transcriptome_assay_ids.len(), 2);
}

#[test]
fn specificity_matrix_same_reference_label_keeps_distinct_fingerprint_history() {
    let (mut engine, panel) = fixture();
    let old = acceptance(
        &panel,
        "same_reference",
        BlastDatabaseIndexKind::GenomicDna,
        true,
    );
    retain(&mut engine, &panel, &old);
    let mut replacement = acceptance(
        &panel,
        "same_reference",
        BlastDatabaseIndexKind::GenomicDna,
        false,
    );
    replacement.acceptance_id.push_str("_replacement");
    for assessment in &mut replacement.assessments {
        let report = &mut assessment.report;
        let database = report.blast_database.as_mut().unwrap();
        database.content_fingerprint = Some("sha256:replacement_content".into());
        report.intended_target.reference_binding =
            GentleEngine::primer_specificity_reference_identity(database);
    }
    // Report IDs deliberately collide in the global store; snapshots must not.
    retain(&mut engine, &panel, &replacement);
    let matrix = engine
        .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
        .unwrap();
    assert_eq!(matrix.selections.len(), 2);
    assert_ne!(
        matrix.selections[0].reference.content_fingerprint,
        matrix.selections[1].reference.content_fingerprint
    );
    let prior = matrix.selections.iter().find(|s| s.is_superseded).unwrap();
    assert_eq!(prior.cells[0].verdict, PrimerSpecificityMultiVerdict::Fail);
    let active = matrix.selections.iter().find(|s| s.is_active).unwrap();
    assert_eq!(active.cells[0].verdict, PrimerSpecificityMultiVerdict::Pass);
}

#[test]
fn specificity_matrix_tampered_receipts_and_inconsistent_verdicts_never_pass() {
    let (mut engine, panel) = fixture();
    retain(
        &mut engine,
        &panel,
        &acceptance(
            &panel,
            "reference",
            BlastDatabaseIndexKind::GenomicDna,
            false,
        ),
    );
    let baseline = engine.read_primer_design_store();
    for field in [
        "hash",
        "reference",
        "policy",
        "pair",
        "missing",
        "duplicate",
        "incomplete",
        "verdict",
    ] {
        let mut store = baseline.clone();
        let selection = store
            .primer_specificity_reference_selections
            .values_mut()
            .next()
            .unwrap();
        match field {
            "hash" => {
                selection
                    .report_content_sha256_by_assay
                    .insert("assay_a".into(), "sha256:tampered".into());
            }
            "reference" => selection.reference.assembly = Some("other_assembly".into()),
            "policy" => selection.policy.max_3prime_mismatches += 1,
            "pair" => {
                selection
                    .pair_bindings_by_assay
                    .insert("assay_a".into(), "sha256:other_tail".into());
            }
            "missing" => selection
                .acceptance
                .assessments
                .retain(|a| a.assay_id != "assay_a"),
            "duplicate" => selection
                .acceptance
                .assessments
                .push(selection.acceptance.assessments[0].clone()),
            "incomplete" => {
                selection.acceptance.status =
                    TranscriptAssayPanelSpecificityAcceptanceStatus::Incomplete
            }
            "verdict" => {
                selection
                    .acceptance
                    .assessments
                    .iter_mut()
                    .find(|a| a.assay_id == "assay_a")
                    .unwrap()
                    .report
                    .genomic_specificity
                    .status = "fail".into()
            }
            _ => unreachable!(),
        }
        engine.write_primer_design_store(store).unwrap();
        let before = serde_json::to_value(engine.snapshot()).unwrap();
        let matrix = engine
            .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
            .unwrap();
        let cell = matrix.selections[0]
            .cells
            .iter()
            .find(|c| c.assay_id == "assay_a")
            .unwrap();
        assert_eq!(
            cell.verdict,
            PrimerSpecificityMultiVerdict::Incomplete,
            "{field}"
        );
        assert!(!cell.retained_binding_verified, "{field}");
        assert!(
            matrix
                .unassessed_genomic_assay_ids
                .contains(&"assay_a".into())
        );
        assert_eq!(before, serde_json::to_value(engine.snapshot()).unwrap());
    }
}

#[test]
fn specificity_matrix_legacy_and_same_oligo_summaries_remain_unassessed() {
    let (mut engine, panel) = fixture();
    let mut store = engine.read_primer_design_store();
    store
        .transcript_assay_panels
        .get_mut(&panel.report_id)
        .unwrap()
        .genomic_specificity_assessments =
        acceptance(&panel, "legacy", BlastDatabaseIndexKind::GenomicDna, false).assessments;
    let empty = engine
        .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
        .unwrap();
    let pair = empty.assays[0].pair_binding_sha256.clone();
    store.primer_specificity_multi_summaries.insert(
        "unlinked".into(),
        PrimerSpecificityMultiSummary {
            schema: PRIMER_SPECIFICITY_MULTI_SUMMARY_SCHEMA.into(),
            summary_id: "unlinked".into(),
            content_sha256: "sha256:synthetic".into(),
            handoff_id: "synthetic_handoff".into(),
            pair_binding_sha256: pair,
            source_snapshot_sha256: None,
            request: PrimerSpecificityMultiRequest {
                schema: PRIMER_SPECIFICITY_MULTI_REQUEST_SCHEMA.into(),
                pair: PrimerSpecificityMultiPair::SavedPair {
                    primer_report_id: "different_design_context".into(),
                    pair_rank: Some(1),
                    pair_index: None,
                },
                policy: Default::default(),
                references: vec![],
                catalog_path: None,
                cache_dir: None,
            },
            references: vec![],
            genomic: PrimerSpecificityMultiVerdict::Pass,
            transcriptome: PrimerSpecificityMultiVerdict::Pass,
            execution_complete: true,
            execution_manifest: PrimerSpecificityMultiExecutionManifest {
                schema: PRIMER_SPECIFICITY_MULTI_MANIFEST_SCHEMA.into(),
                handoff_id: "synthetic_handoff".into(),
                handoff_content_sha256: "sha256:synthetic".into(),
                pair_binding_sha256: String::new(),
                commands: vec![],
            },
            nonclaims: vec![],
        },
    );
    engine.write_primer_design_store(store).unwrap();
    let matrix = engine
        .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
        .unwrap();
    assert!(matrix.selections.is_empty());
    assert_eq!(matrix.unlinked_pair_summaries.len(), 1);
    assert_eq!(
        matrix.unlinked_pair_summaries[0].same_oligo_current_assay_ids,
        ["assay_a", "assay_b"]
    );
    assert_eq!(matrix.unassessed_genomic_assay_ids, ["assay_a", "assay_b"]);
    assert_eq!(
        matrix.unassessed_transcriptome_assay_ids,
        ["assay_a", "assay_b"]
    );
}

#[test]
fn specificity_matrix_rehashed_inconsistent_verdict_is_not_a_pass() {
    let (mut engine, panel) = fixture();
    let mut bound = acceptance(
        &panel,
        "reference",
        BlastDatabaseIndexKind::GenomicDna,
        false,
    );
    bound
        .assessments
        .iter_mut()
        .find(|a| a.assay_id == "assay_a")
        .unwrap()
        .report
        .genomic_specificity
        .status = "fail".into();
    // Recompute every selection/report binding via the existing recorder. This
    // tests verdict consistency independently of stale hashes, not authenticity.
    retain(&mut engine, &panel, &bound);
    let matrix = engine
        .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
        .unwrap();
    let cell = matrix.selections[0]
        .cells
        .iter()
        .find(|c| c.assay_id == "assay_a")
        .unwrap();
    assert_eq!(cell.recorded_status.as_deref(), Some("pass"));
    assert_eq!(cell.verdict, PrimerSpecificityMultiVerdict::Incomplete);
    assert!(
        cell.diagnostics
            .iter()
            .any(|s| s == "retained_verdict_inconsistent")
    );
    assert!(
        !cell
            .diagnostics
            .iter()
            .any(|s| s == "retained_report_hash_or_identity_mismatch")
    );
}

#[test]
fn specificity_matrix_shared_shell_raw_operation_and_reload_agree() {
    let (mut engine, panel) = fixture();
    retain(
        &mut engine,
        &panel,
        &acceptance(
            &panel,
            "reference",
            BlastDatabaseIndexKind::GenomicDna,
            true,
        ),
    );
    let before = serde_json::to_value(engine.snapshot()).unwrap();
    let expected = serde_json::to_value(
        engine
            .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
            .unwrap(),
    )
    .unwrap();
    let tutorial =
        include_str!("../../../docs/tutorial/04-09_reference_bound_primer_specificity.md");
    let template = tutorial
        .lines()
        .filter_map(|line| line.strip_prefix("gentle_cli --state study.gentle.json "))
        .find(|line| line.starts_with("primers transcript-assay-specificity-matrix "))
        .expect("tutorial must exercise the real matrix route");
    assert_eq!(template.matches("YOUR_PANEL_REPORT_ID").count(), 1);
    let command = crate::engine_shell::parse_shell_line(&template.replace(
        "YOUR_PANEL_REPORT_ID",
        &crate::engine_shell::shell_quote(RAW_PANEL_REPORT_ID),
    ))
    .unwrap();
    assert!(!command.is_state_mutating());
    let result = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
    assert!(!result.state_changed);
    assert_eq!(
        result.output["result"]["transcript_assay_specificity_matrix"],
        expected
    );
    let raw = crate::engine_shell::ShellCommand::Op {
        payload: serde_json::to_string(&Operation::InspectTranscriptAssayPanelSpecificityMatrix {
            panel_report_id: RAW_PANEL_REPORT_ID.into(),
        })
        .unwrap(),
    };
    assert!(!raw.is_state_mutating());
    let result = crate::engine_shell::execute_shell_command(&mut engine, &raw).unwrap();
    assert_eq!(
        result.output["result"]["transcript_assay_specificity_matrix"],
        expected
    );
    assert_eq!(before, serde_json::to_value(engine.snapshot()).unwrap());
    let restored = GentleEngine::from_state(serde_json::from_value(before).unwrap());
    assert_eq!(
        serde_json::to_value(
            restored
                .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
                .unwrap()
        )
        .unwrap(),
        expected
    );
    let capabilities = crate::engine_shell::execute_shell_command(
        &mut engine,
        &crate::engine_shell::parse_shell_line("introspect capabilities").unwrap(),
    )
    .unwrap();
    assert!(
        capabilities
            .output
            .to_string()
            .contains("transcript-assay-specificity-matrix")
    );
}

#[test]
fn specificity_matrix_mcp_reads_same_payload_without_persisting_state() {
    let (mut engine, panel) = fixture();
    retain(
        &mut engine,
        &panel,
        &acceptance(
            &panel,
            "reference",
            BlastDatabaseIndexKind::GenomicDna,
            false,
        ),
    );
    let expected = serde_json::to_value(
        engine
            .inspect_transcript_assay_panel_specificity_matrix(&panel.report_id)
            .unwrap(),
    )
    .unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("state ' with spaces.json");
    engine.state().save_to_path(path.to_str().unwrap()).unwrap();
    let before = std::fs::read(&path).unwrap();
    let output = crate::mcp_server::mcp_tool_call_for_capability_surface_tests(
        path.to_str().unwrap(),
        "transcript_assay_specificity_matrix",
        json!({"panel_report_id":RAW_PANEL_REPORT_ID}),
    );
    assert_eq!(output["isError"], false);
    let content = output["content"].as_array().unwrap();
    let payload: serde_json::Value =
        serde_json::from_str(content.iter().find_map(|c| c["text"].as_str()).unwrap()).unwrap();
    assert_eq!(
        payload["result"]["transcript_assay_specificity_matrix"],
        expected
    );
    assert_eq!(before, std::fs::read(path).unwrap());
}
