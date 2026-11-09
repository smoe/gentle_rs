//! Synthetic inline fixtures for external-report orchestration, not biological evidence.
//! Unix integration reuses GENtle's invented transcript panel and mock local BLAST tools.

use super::*;

fn request() -> TranscriptAssayPanelMultiReferenceRequest {
    serde_json::from_value(json!({
        "schema": PANEL_MULTI_REQUEST_SCHEMA, "panel_report_id":"panel_B",
        "references":[
            {"target_genome_id":"synthetic-a","expected_index_kind":"genomic_dna","required":true},
            {"target_genome_id":"synthetic-b","expected_index_kind":"transcriptome_cdna","required":true},
            {"target_genome_id":"synthetic-optional","expected_index_kind":"genomic_dna","required":false}
        ]
    })).unwrap()
}

fn bound() -> TranscriptAssayPanelMultiReferenceHandoff {
    let request = request();
    let pair = PrimerSpecificityHandoff {
        commands: ["F", "R"]
            .into_iter()
            .enumerate()
            .map(|(i, id)| PrimerSpecificityHandoffCommand {
                command_id: id.into(),
                role: if i == 0 {
                    PrimerSpecificityPrimerRole::Forward
                } else {
                    PrimerSpecificityPrimerRole::Reverse
                },
                output_tsv_path: format!("{id}.tsv"),
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    };
    TranscriptAssayPanelMultiReferenceHandoff {
        schema: PANEL_MULTI_HANDOFF_SCHEMA.into(),
        handoff_id: "h".into(),
        content_sha256: "science".into(),
        bundle_sha256: "bundle".into(),
        panel_report_id: request.panel_report_id.clone(),
        request,
        panel_digest: "design".into(),
        source_snapshot_sha256: "source".into(),
        effective_limits: LIMITS,
        expected_command_count: 2,
        references: vec![PanelSpecificityChild {
            requested_genome_id: "g".into(),
            resolved_genome_id: "g".into(),
            required: true,
            expected_index_kind: BlastDatabaseIndexKind::GenomicDna,
            availability: PrimerSpecificityMultiAvailability::Prepared,
            diagnostic: None,
            reference: None,
            child: Some(TranscriptAssayPanelSpecificityHandoff {
                assays: vec![TranscriptAssayPanelSpecificityHandoffAssay {
                    handoff: pair,
                    ..Default::default()
                }],
                ..Default::default()
            }),
        }],
        manifest_path: "execution_manifest.json".into(),
        nonclaims: nonclaims(),
    }
}

#[test]
fn transcript_assay_panel_multi_admission_and_smaller_effective_limits_fail_closed() {
    let mut r = request();
    assert!(request_refs(&r, LIMITS).is_ok());
    let original = r.references.clone();
    r.references.clear();
    assert!(request_refs(&r, LIMITS).is_err());
    r.references = vec![original[0].clone(); LIMITS.max_references + 1];
    assert!(request_refs(&r, LIMITS).is_err());
    r.references = vec![original[0].clone(), original[0].clone()];
    assert!(request_refs(&r, LIMITS).is_err());
    r.references = vec![original[2].clone()];
    assert!(request_refs(&r, LIMITS).is_err());
    r.references = original;
    r.references[0].target_genome_id = "*".into();
    assert!(request_refs(&r, LIMITS).is_err());
    r = request();
    r.policy.specificity_target_genome_id = Some("hidden".into());
    assert!(request_refs(&r, LIMITS).is_err());
    let mut value = serde_json::to_value(request()).unwrap();
    value["references"][0]["intended_target"] = json!({});
    assert!(serde_json::from_value::<TranscriptAssayPanelMultiReferenceRequest>(value).is_err());
    assert!(LIMITS.commands(0, 1).is_err());
    assert!(LIMITS.commands(LIMITS.max_assays + 1, 1).is_err());
    assert!(
        LIMITS
            .commands(1, LIMITS.max_prepared_references + 1)
            .is_err()
    );
    let mut tight = LIMITS;
    tight.max_commands = 3;
    assert!(tight.validate().is_ok());
    assert!(tight.commands(1, 2).is_err());
    tight.max_total_output_bytes += 1;
    assert!(tight.validate().is_err());
}

#[test]
fn transcript_assay_panel_multi_manifest_rejects_omitted_duplicate_unknown_and_cross_bundle_rows() {
    let h = bound();
    let m = panel_manifest(&h).unwrap();
    assert!(manifest_rows(&h, &m).is_ok());
    let mut fake_pending = m.clone();
    fake_pending.commands[0].output_size_bytes = Some(0);
    assert!(manifest_rows(&h, &fake_pending).is_err());
    for change in 0..6 {
        let mut wrong = m.clone();
        match change {
            0 => {
                wrong.commands.pop();
            }
            1 => wrong.commands[1] = wrong.commands[0].clone(),
            2 => wrong.commands[1].command_id = "unknown".into(),
            3 => wrong.commands[1].output_path = "other.tsv".into(),
            4 => wrong.commands[1].command_sha256 = "other".into(),
            _ => wrong.bundle_sha256 = "other".into(),
        }
        assert!(manifest_rows(&h, &wrong).is_err());
    }
    let mut duplicate = h.clone();
    duplicate.references[0].child.as_mut().unwrap().assays[0]
        .handoff
        .commands[1]
        .command_id = "F".into();
    assert!(panel_manifest(&duplicate).is_err());
}

#[test]
fn transcript_assay_panel_multi_budget_refusal_precedes_file_access_and_parsing() {
    let mut h = bound();
    h.effective_limits.max_total_output_bytes = 9;
    let mut m = panel_manifest(&h).unwrap();
    for c in &mut m.commands {
        c.state = PrimerSpecificityMultiExecutionState::Completed;
        c.exit_code = Some(0);
        c.output_sha256 = Some(sha256_prefixed_bytes(b""));
        c.output_size_bytes = Some(5);
    }
    assert!(
        manifest_rows(&h, &m)
            .unwrap_err()
            .message
            .contains("total declared")
    );
    m.commands[1].output_size_bytes = None;
    assert!(
        manifest_rows(&h, &m)
            .unwrap_err()
            .message
            .contains("size unknown")
    );
    let root = tempfile::tempdir().unwrap();
    let root = root.path().canonicalize().unwrap();
    let mut pair = h.references[0]
        .child
        .take()
        .unwrap()
        .assays
        .remove(0)
        .handoff;
    let mut primers = Vec::new();
    for c in &mut pair.commands {
        c.query_label = c.command_id.clone();
        c.output_tsv_path = root
            .join(format!("{}.tsv", c.command_id))
            .to_string_lossy()
            .into_owned();
        c.query_fasta_path = root
            .join(format!("{}.fa", c.command_id))
            .to_string_lossy()
            .into_owned();
        fs::write(&c.query_fasta_path, format!(">{}\nA\n", c.query_label)).unwrap();
        primers.push(PrimerSpecificityInputPrimer {
            role: c.role,
            annealing_sequence: "A".into(),
            ..Default::default()
        });
    }
    let text = "bad row\n".repeat(LIMITS.max_parsed_rows + 1);
    fs::write(&pair.commands[0].output_tsv_path, &text).unwrap();
    let c = &mut m.commands[0];
    c.output_size_bytes = Some(text.len() as u64);
    c.output_sha256 = Some(sha256_prefixed_bytes(text.as_bytes()));
    let runs = m
        .commands
        .iter()
        .map(|c| (c.command_id.clone(), c))
        .collect();
    let mut budget = SpecificityOutputBudget {
        per_output: LIMITS.max_bytes_per_output,
        bytes_left: LIMITS.max_total_output_bytes,
        rows_left: LIMITS.max_parsed_rows,
    };
    let mut complete = true;
    let error = GentleEngine::specificity_multi_read_outputs(
        &pair,
        &primers,
        &runs,
        &mut vec![],
        &mut complete,
        Some(&mut budget),
    )
    .unwrap_err();
    assert!(error.message.contains("row limit"), "{error:?}");
    pair.commands[1].role = pair.commands[0].role;
    assert!(
        GentleEngine::specificity_multi_read_outputs(
            &pair,
            &primers,
            &runs,
            &mut vec![],
            &mut complete,
            None
        )
        .unwrap_err()
        .message
        .contains("role")
    );
    pair.commands.clear();
    assert!(
        GentleEngine::specificity_multi_read_outputs(
            &pair,
            &primers,
            &runs,
            &mut vec![],
            &mut complete,
            None
        )
        .is_err()
    );
}

#[test]
fn transcript_assay_panel_multi_verdicts_keep_required_failure_and_optional_diagnostics_distinct() {
    use PrimerSpecificityMultiVerdict::{Fail, Incomplete, Pass};
    let row = |required, verdict| PrimerSpecificityMultiSummaryReference {
        genome_id: "synthetic".into(),
        required,
        index_kind: BlastDatabaseIndexKind::GenomicDna,
        reference: None,
        verdict,
        applicability: "fixture".into(),
        diagnostics: vec![],
        report: None,
    };
    let mut cells = vec![row(true, Pass), row(false, Fail)];
    let verdict = |cells: &[PrimerSpecificityMultiSummaryReference], complete| {
        GentleEngine::specificity_multi_dimension(
            cells,
            BlastDatabaseIndexKind::GenomicDna,
            complete,
        )
    };
    assert_eq!(verdict(&cells, true), Pass);
    assert_eq!(verdict(&cells, false), Incomplete);
    cells[0].verdict = Fail;
    assert_eq!(verdict(&cells, false), Fail);
    assert_eq!(cells[1].verdict, Fail);
}

#[test]
fn transcript_assay_panel_multi_destination_rollback_leaves_existing_evidence_alone() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path().canonicalize().unwrap();
    let dir = root.join("panel ' \u{03b2}");
    assert!(
        with_fresh_specificity_bundle(dir.to_str().unwrap(), |owned| {
            fs::create_dir(owned.join("reference_00")).unwrap();
            Err::<(), _>(invalid("second child refused"))
        })
        .is_err()
    );
    assert!(!dir.exists());
    fs::create_dir(&dir).unwrap();
    fs::write(dir.join("foreign"), b"retained").unwrap();
    assert!(with_fresh_specificity_bundle(dir.to_str().unwrap(), |_| Ok(())).is_err());
    assert_eq!(fs::read(dir.join("foreign")).unwrap(), b"retained");
    // Joining onto a Windows verbatim path can normalize away the test's `..`.
    let raw = format!(
        "{}{separator}child{separator}..{separator}escape",
        root.to_str().unwrap(),
        separator = std::path::MAIN_SEPARATOR
    );
    assert!(std::path::Path::new(&raw).components().any(|component| {
        component == std::path::Component::ParentDir || component.as_os_str() == ".."
    }));
    let error = with_fresh_specificity_bundle(&raw, |_| Ok(())).unwrap_err();
    assert_eq!(error.message, "Unsafe multi-reference output path");
    assert!(!root.join("escape").exists());
    assert_eq!(fs::read(dir.join("foreign")).unwrap(), b"retained");
}

#[cfg(unix)]
#[test]
fn transcript_assay_panel_multi_bind_import_replay_partial_and_atomic_refusal() {
    use super::super::tests::{
        transcript_assay_panel_operation, transcript_assay_panel_relaxed_side,
        transcript_qpcr_panel_test_engine,
    };
    use crate::engine::primer_specificity_multi::tests::synthetic_specificity_resources;
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path().canonicalize().unwrap();
    let (catalog, _, _tools) = synthetic_specificity_resources(&root);
    let catalog_name = catalog.to_str().unwrap();
    let mut engine = transcript_qpcr_panel_test_engine();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "synthetic-a".into(),
            catalog_path: Some(catalog_name.into()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    let catalog = crate::genomes::GenomeCatalog::from_json_file(catalog_name).unwrap();
    catalog
        .prepare_blast_sequence_resource("synthetic-b", None)
        .unwrap();
    let mut design = transcript_assay_panel_operation(
        TranscriptAssayCoveragePolicy::BestEffort,
        transcript_assay_panel_relaxed_side(),
        1,
        "panel_B",
    );
    if let Operation::DesignTranscriptAssayPanel { objective, .. } = &mut design {
        *objective = TranscriptAssayPanelObjective::OnePerClass;
    }
    let mut panel = *engine
        .apply(design)
        .unwrap()
        .transcript_assay_panel
        .unwrap();
    assert!(panel.selected_assays.len() >= 2);
    panel.selected_assays.truncate(2);
    // A separate assay with the same oligos still needs its own command receipts.
    let mut third = panel.selected_assays[0].clone();
    third.rank = panel.selected_assays.iter().map(|a| a.rank).max().unwrap() + 1;
    third.assay_id = "synthetic_repeat_assay".into();
    panel.selected_assays.push(third);
    let panel_id = panel.report_id.clone();
    let mut store = engine.read_primer_design_store();
    store
        .transcript_assay_panels
        .insert(panel_id.clone(), panel.clone());
    let identity = GentleEngine::primer_specificity_reference_identity(
        &catalog
            .inspect_blast_database("synthetic-a", None)
            .unwrap()
            .unwrap(),
    )
    .unwrap();
    store.primer_specificity_panel_sources.insert(
        panel_id.clone(),
        PrimerSpecificityPanelSource {
            panel_design_digest: GentleEngine::primer_specificity_selection_panel_digest(&panel)
                .unwrap(),
            source_reference: PrimerSpecificitySourceReference {
                genome_id: identity.genome_id,
                assembly: identity.assembly,
                release: identity.release,
                ..Default::default()
            },
        },
    );
    engine.write_primer_design_store(store).unwrap();
    let untouched = |e: &GentleEngine| {
        let s = e.read_primer_design_store();
        json!([
            s.transcript_assay_panels,
            s.primer_specificity_panel_sources,
            s.primer_specificity_reference_selections,
            s.active_primer_specificity_reference_selections
        ])
    };
    let original_panel = untouched(&engine);
    let mut req = request();
    req.panel_report_id = panel_id.clone();
    req.catalog_path = Some(catalog_name.into());
    req.policy.full_alignment.mode = PrimerSpecificityFullAlignmentMode::Disabled;
    let dir = root.join("panel ' \u{03b2}");
    let hpath = dir.join("handoff.json");
    let hp = hpath.to_str().unwrap();
    let before = json!(engine.snapshot());
    let mut bad = req.clone();
    bad.references[2].required = true;
    assert!(
        engine
            .prepare_panel_multi_specificity(bad, dir.to_str().unwrap())
            .is_err()
    );
    let mut bad = req.clone();
    bad.references[0].expected_index_kind = BlastDatabaseIndexKind::TranscriptomeCdna;
    assert!(
        engine
            .prepare_panel_multi_specificity(bad, dir.to_str().unwrap())
            .is_err()
    );
    let mut bad = req.clone();
    bad.expected_panel_digest = Some("stale".into());
    assert!(
        engine
            .prepare_panel_multi_specificity(bad, dir.to_str().unwrap())
            .is_err()
    );
    let mut bad = req.clone();
    bad.panel_report_id = "no such panel".into();
    assert!(
        engine
            .prepare_panel_multi_specificity(bad, dir.to_str().unwrap())
            .is_err()
    );
    assert!(!dir.exists());
    let h = engine
        .prepare_panel_multi_specificity(req.clone(), dir.to_str().unwrap())
        .unwrap();
    let h2 = engine
        .prepare_panel_multi_specificity(req.clone(), root.join("second").to_str().unwrap())
        .unwrap();
    assert_eq!(h.handoff_id, h2.handoff_id);
    assert_eq!(h.content_sha256, h2.content_sha256);
    assert_ne!(h.bundle_sha256, h2.bundle_sha256);
    assert_eq!(h.expected_command_count, 12);
    assert_eq!(h.references.len(), 3);
    assert!(h.references[2].child.is_none());
    assert!(h.references[2].diagnostic.is_some());
    assert_eq!(json!(engine.snapshot()), before);
    for r in &h.references[..2] {
        let child = r.child.as_ref().unwrap();
        assert_eq!(child.assays.len(), 3);
        let direct = engine
            .prepare_transcript_assay_panel_specificity_handoff(
                &panel_id,
                &r.resolved_genome_id,
                req.policy.clone(),
                Some(catalog_name),
                None,
                root.join(format!("direct_{}", r.resolved_genome_id))
                    .to_str()
                    .unwrap(),
            )
            .unwrap();
        for (a, d) in child.assays.iter().zip(direct.assays) {
            assert_eq!(a.handoff.handoff_id, d.handoff.handoff_id);
            assert_eq!(a.primer_pair_digest, d.primer_pair_digest);
            assert!(multi_equal(
                &a.handoff.intended_target,
                &d.handoff.intended_target
            ));
            if r.expected_index_kind == BlastDatabaseIndexKind::GenomicDna {
                assert_eq!(a.handoff.intended_target.reference_binding, r.reference);
            } else {
                assert!(
                    a.handoff.intended_target.reference_binding.is_none(),
                    "no coordinate transfer to a different reference ID"
                );
            }
        }
    }
    let mut m = panel_manifest(&h).unwrap();
    let mp = h.manifest_path.as_str();
    for c in &mut m.commands {
        fs::write(&c.output_path, b"").unwrap();
        c.state = PrimerSpecificityMultiExecutionState::Completed;
        c.exit_code = Some(0);
        c.output_size_bytes = Some(0);
        c.output_sha256 = Some(sha256_prefixed_bytes(b""));
    }
    specificity_multi_write_json(Path::new(mp), &m).unwrap();
    let summary = engine.import_panel_multi_specificity(hp, mp).unwrap();
    assert!(summary.execution_complete);
    assert_eq!(summary.assays.len(), 3);
    assert_eq!(
        summary.retained_command_ids.len(),
        12,
        "zero bytes with a verified receipt are not missing output"
    );
    for a in &summary.assays {
        assert_eq!(a.cells.len(), 3);
        assert!(a.cells[..2].iter().all(|r| r.report.is_some()));
        assert_eq!(a.transcriptome, PrimerSpecificityMultiVerdict::Incomplete);
        assert_eq!(a.cells[2].applicability, "unavailable_at_preparation");
    }
    assert_eq!(untouched(&engine), original_panel);
    let mut reopened =
        GentleEngine::from_state(serde_json::from_value(json!(engine.snapshot())).unwrap());
    assert_eq!(reopened.list_panel_multi_specificity().len(), 1);
    let stable = json!(reopened.snapshot());
    assert_eq!(
        json!(
            reopened
                .get_panel_multi_specificity(&summary.summary_id)
                .unwrap()
        ),
        json!(summary)
    );
    assert_eq!(
        json!(reopened.import_panel_multi_specificity(hp, mp).unwrap()),
        json!(summary)
    );
    assert_eq!(json!(reopened.snapshot()), stable);

    let mut partial = m.clone();
    for (i, state) in [
        PrimerSpecificityMultiExecutionState::Cancelled,
        PrimerSpecificityMultiExecutionState::Failed,
    ]
    .into_iter()
    .enumerate()
    {
        let c = &mut partial.commands[i];
        c.state = state;
        c.exit_code = None;
        c.output_size_bytes = None;
        c.output_sha256 = None;
    }
    specificity_multi_write_json(Path::new(mp), &partial).unwrap();
    let partial = engine.import_panel_multi_specificity(hp, mp).unwrap();
    assert!(!partial.execution_complete);
    assert_eq!(partial.retained_command_ids.len(), 10);
    assert!(!partial.assays[0].execution_complete);
    assert!(partial.assays[1].execution_complete);
    assert_eq!(untouched(&engine), original_panel);
    let stable = json!(engine.snapshot());
    let cancelled_query = &h.references[0].child.as_ref().unwrap().assays[0]
        .handoff
        .commands[0]
        .query_fasta_path;
    let original_query = fs::read(cancelled_query).unwrap();
    fs::write(cancelled_query, b"changed pending query").unwrap();
    assert!(
        engine
            .import_panel_multi_specificity(hp, mp)
            .unwrap_err()
            .message
            .contains("query_sequence_changed")
    );
    fs::write(cancelled_query, original_query).unwrap();
    assert_eq!(json!(engine.snapshot()), stable);
    let mut budget = m.clone();
    budget.commands[0].output_size_bytes = Some(LIMITS.max_total_output_bytes + 1);
    fs::remove_file(&budget.commands[0].output_path).unwrap();
    specificity_multi_write_json(Path::new(mp), &budget).unwrap();
    assert!(
        engine
            .import_panel_multi_specificity(hp, mp)
            .unwrap_err()
            .message
            .contains("byte limit")
    );
    specificity_multi_write_json(Path::new(mp), &m).unwrap();
    fs::write(&m.commands[0].output_path, b"tampered").unwrap();
    assert!(engine.import_panel_multi_specificity(hp, mp).is_err());
    fs::write(&m.commands[0].output_path, b"").unwrap();
    let query = &h.references[0].child.as_ref().unwrap().assays[0]
        .handoff
        .commands[0]
        .query_fasta_path;
    let saved = fs::read(query).unwrap();
    fs::write(query, b">changed\nA\n").unwrap();
    assert!(engine.import_panel_multi_specificity(hp, mp).is_err());
    fs::write(query, saved).unwrap();
    let mut wrong = h.clone();
    wrong.references[0]
        .child
        .as_mut()
        .unwrap()
        .policy
        .max_3prime_mismatches += 1;
    wrong.bundle_sha256 = bundle_digest(&wrong).unwrap();
    specificity_multi_write_json(&hpath, &wrong).unwrap();
    assert!(engine.import_panel_multi_specificity(hp, mp).is_err());
    specificity_multi_write_json(&hpath, &h).unwrap();
    let child = h.references[1].child.as_ref().unwrap();
    for file in [&child.handoff_path, &child.assays[2].handoff.handoff_path] {
        let saved = fs::read(file).unwrap();
        fs::write(file, b"{}").unwrap();
        assert!(
            engine
                .import_panel_multi_specificity(hp, mp)
                .unwrap_err()
                .message
                .contains("child file")
        );
        fs::write(file, saved).unwrap();
    }
    assert_eq!(
        json!(engine.snapshot()),
        stable,
        "late failures must leave no partial write"
    );

    for mutation in 0..3 {
        let mut changed = GentleEngine::from_state(serde_json::from_value(stable.clone()).unwrap());
        if mutation == 2 {
            changed.state_mut().sequences.remove(&panel.source_seq_id);
        } else {
            let mut s = changed.read_primer_design_store();
            let p = s.transcript_assay_panels.get_mut(&panel_id).unwrap();
            if mutation == 0 {
                p.generated_at_unix_ms += 1;
            } else {
                p.selected_assays[0]
                    .primer_pair
                    .forward
                    .non_annealing_5prime_tail_bp += 1;
            }
            changed.write_primer_design_store(s).unwrap();
        }
        let changed_state = json!(changed.snapshot());
        assert!(changed.import_panel_multi_specificity(hp, mp).is_err());
        assert_eq!(json!(changed.snapshot()), changed_state);
    }

    // A separately bound cDNA source can retain transcript products; never transfer
    // its coordinates back to the genomic sibling just because the release matches.
    let id = GentleEngine::primer_specificity_reference_identity(
        &catalog
            .inspect_blast_database("synthetic-b", None)
            .unwrap()
            .unwrap(),
    )
    .unwrap();
    let mut changed = GentleEngine::from_state(serde_json::from_value(stable.clone()).unwrap());
    let mut s = changed.read_primer_design_store();
    s.primer_specificity_panel_sources
        .get_mut(&panel_id)
        .unwrap()
        .source_reference = PrimerSpecificitySourceReference {
        genome_id: id.genome_id.clone(),
        assembly: id.assembly.clone(),
        release: id.release.clone(),
        ..Default::default()
    };
    changed.write_primer_design_store_snapshot(s).unwrap();
    let cdna = changed
        .prepare_panel_multi_specificity(req, root.join("cdna source").to_str().unwrap())
        .unwrap();
    for a in &cdna.references[1].child.as_ref().unwrap().assays {
        assert_eq!(
            a.handoff.intended_target.reference_binding.as_ref(),
            Some(&id)
        );
        assert!(
            a.handoff
                .intended_target
                .expected_products
                .iter()
                .any(|p| p.target_space == "transcriptome_cdna")
        );
    }
    for a in &cdna.references[0].child.as_ref().unwrap().assays {
        assert!(a.handoff.intended_target.reference_binding.is_none());
    }
    let prefix = &h.references[0].child.as_ref().unwrap().assays[0]
        .handoff
        .blast_db_prefix;
    fs::write(format!("{prefix}.nsq"), b"replaced-index").unwrap();
    let stale = engine.import_panel_multi_specificity(hp, mp).unwrap();
    assert!(
        stale
            .assays
            .iter()
            .all(|a| a.genomic == PrimerSpecificityMultiVerdict::Incomplete
                && a.cells[0].applicability == "stale_or_unavailable_reference")
    );
    assert_eq!(untouched(&engine), original_panel);
}
