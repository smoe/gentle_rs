//! Frozen, synthetic pre-profile records exercise the actual stored identity seams.
//! See test_files/fixtures/primer_specificity_identity/README.md for provenance.

use super::*;
use serde::de::DeserializeOwned;

const FIXTURE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/test_files/fixtures/primer_specificity_identity"
);

fn fixture() -> Value {
    let root = Path::new(FIXTURE_DIR);
    let bytes = fs::read(root.join("records.json")).unwrap();
    let sums = fs::read_to_string(root.join("SHA256SUMS")).unwrap();
    assert_eq!(
        sums.trim(),
        format!(
            "{}  records.json",
            sha256_prefixed_bytes(&bytes).trim_start_matches("sha256:")
        )
    );
    serde_json::from_slice(&bytes).unwrap()
}

fn stored<T: DeserializeOwned + Serialize>(records: &Value, key: &str) -> T {
    let value = &records[key];
    let restored: T = serde_json::from_value(value.clone()).unwrap();
    assert_eq!(
        serde_json::to_value(&restored).unwrap(),
        *value,
        "{key} serialization drift"
    );
    restored
}

#[test]
fn stored_primer_identity_policy_and_single_handoff_survive_reopening() {
    let records = fixture();
    let policy: PrimerSpecificityPolicy = stored(&records, "policy");
    assert_eq!(
        sha256_prefixed_bytes(&serde_json::to_vec(&policy).unwrap()),
        records["policy_sha256"]
    );
    let handoff: PrimerSpecificityHandoff = stored(&records, "single_handoff");
    assert_eq!(
        GentleEngine::primer_specificity_handoff_id_from_record(&handoff).unwrap(),
        handoff.handoff_id
    );
    let mut changed = handoff;
    changed.policy.max_3prime_mismatches += 1;
    assert_ne!(
        GentleEngine::primer_specificity_handoff_id_from_record(&changed).unwrap(),
        changed.handoff_id
    );
}

#[test]
fn stored_primer_identity_multi_handoff_digest_survives_reopening() {
    let records = fixture();
    let handoff: PrimerSpecificityMultiHandoff = stored(&records, "multi_handoff");
    assert_eq!(
        handoff.nonclaims,
        MULTI_NONCLAIMS
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
    );
    let digest = GentleEngine::specificity_multi_handoff_digest(&handoff).unwrap();
    assert_eq!(digest, handoff.content_sha256);
    assert_eq!(
        handoff.handoff_id,
        format!("primer_multi_{}", digest.trim_start_matches("sha256:"))
    );
    let child = handoff.references[0].child.as_ref().unwrap();
    assert_eq!(
        GentleEngine::primer_specificity_handoff_id_from_record(child).unwrap(),
        child.handoff_id
    );
    let mut changed = handoff;
    changed.references[0].child.as_mut().unwrap().commands[0]
        .args
        .push("changed".into());
    assert_ne!(
        GentleEngine::specificity_multi_handoff_digest(&changed).unwrap(),
        changed.content_sha256
    );
}

#[test]
fn stored_primer_identity_summary_and_report_hashes_survive_reopening() {
    let records = fixture();
    let report: PrimerSpecificityReport = stored(&records, "report");
    assert_eq!(
        sha256_prefixed_bytes(&serde_json::to_vec(&report).unwrap()),
        records["report_sha256"]
    );
    assert!(report.forward_hits[0].full_alignment.is_some());
    assert!(report.reverse_hits[0].full_alignment.is_none());
    assert!(!report.amplicons.is_empty());
    let summary: PrimerSpecificityMultiSummary = stored(&records, "summary");
    assert_eq!(
        GentleEngine::specificity_multi_summary_digest(&summary).unwrap(),
        summary.content_sha256
    );
    assert_eq!(
        summary.summary_id,
        format!(
            "primer_multi_summary_{}",
            summary.content_sha256.trim_start_matches("sha256:")
        )
    );
    let mut changed = summary;
    changed.references[0].report.as_mut().unwrap().forward_hits[0]
        .full_alignment
        .as_mut()
        .unwrap()
        .three_prime_mismatches += 1;
    assert_ne!(
        GentleEngine::specificity_multi_summary_digest(&changed).unwrap(),
        changed.content_sha256
    );
}

#[test]
fn stored_primer_identity_panel_selection_remains_selected() {
    let records = fixture();
    let panel: TranscriptAssayPanelReport = stored(&records, "panel");
    // The store envelope may grow unrelated collections; only bound records
    // and their real lookup/hash seams are the compatibility contract here.
    let store: PrimerDesignStore = serde_json::from_value(records["store"].clone()).unwrap();
    assert_eq!(
        GentleEngine::primer_specificity_selection_panel_digest(&panel).unwrap(),
        records["panel_digest"]
    );
    let selection = store
        .primer_specificity_reference_selections
        .values()
        .next()
        .unwrap();
    assert_eq!(selection.selection_id, records["selection_id"]);
    assert_eq!(
        selection.report_content_sha256_by_assay["synthetic_assay"],
        records["report_sha256"]
    );
    let mut engine = GentleEngine::new();
    engine.write_primer_design_store(store).unwrap();
    let restored = engine
        .selected_primer_specificity_report(
            &panel,
            "synthetic_assay",
            BlastDatabaseIndexKind::TranscriptomeCdna,
        )
        .expect("saved report must not silently lose its active selection");
    assert_eq!(serde_json::to_value(restored).unwrap(), records["report"]);
    let summary = engine
        .get_primer_pair_multi_reference_specificity_summary(
            records["summary"]["summary_id"].as_str().unwrap(),
        )
        .unwrap();
    assert_eq!(
        GentleEngine::specificity_multi_summary_digest(&summary).unwrap(),
        summary.content_sha256
    );

    // Exercise the hash guard itself, not merely equality with a default object.
    let mut changed = engine.read_primer_design_store();
    changed
        .primer_specificity_reference_selections
        .values_mut()
        .next()
        .unwrap()
        .acceptance
        .assessments[0]
        .report
        .forward_hits[0]
        .effective_mismatches += 1;
    engine.write_primer_design_store(changed).unwrap();
    assert!(
        engine
            .selected_primer_specificity_report(
                &panel,
                "synthetic_assay",
                BlastDatabaseIndexKind::TranscriptomeCdna
            )
            .is_none()
    );
}

#[test]
fn stored_primer_identity_selection_can_be_recomputed_from_saved_acceptance() {
    let records = fixture();
    let panel: TranscriptAssayPanelReport = stored(&records, "panel");
    let original: PrimerDesignStore = serde_json::from_value(records["store"].clone()).unwrap();
    let selection = original
        .primer_specificity_reference_selections
        .values()
        .next()
        .unwrap();
    let mut recomputed = PrimerDesignStore::default();
    GentleEngine::record_primer_specificity_reference_selection(
        &mut recomputed,
        &panel,
        &selection.acceptance,
    )
    .unwrap();
    assert_eq!(
        serde_json::to_value(&recomputed.primer_specificity_reference_selections).unwrap(),
        serde_json::to_value(&original.primer_specificity_reference_selections).unwrap()
    );
    assert_eq!(
        recomputed.active_primer_specificity_reference_selections,
        original.active_primer_specificity_reference_selections
    );
}

// Capture uses only synthetic literals and baseline serializers/identity helpers.
// Never replace the golden records to accommodate a changed serializer.
#[test]
#[ignore = "explicit baseline capture only; requires a new output directory"]
fn capture_primer_identity_baseline() {
    let output = std::env::var_os("GENTLE_PRIMER_IDENTITY_CAPTURE_DIR")
        .expect("set a NEW capture directory");
    let output = Path::new(&output);
    fs::create_dir(output).expect("capture must not overwrite existing records");
    let records = baseline_records();
    let mut bytes = serde_json::to_vec_pretty(&records).unwrap();
    bytes.push(b'\n');
    fs::write(output.join("records.json"), &bytes).unwrap();
    fs::write(
        output.join("SHA256SUMS"),
        format!(
            "{}  records.json\n",
            sha256_prefixed_bytes(&bytes).trim_start_matches("sha256:")
        ),
    )
    .unwrap();
}

fn baseline_records() -> Value {
    let policy = PrimerSpecificityPolicy::default();
    let outputs = [
        "forward\tsynthetic_transcript.1\t100\t20\t0\t0\t1\t20\t1\t20\t0.001\t40\t100\t20\tACGTACGTACGTACGTACGT\tACGTACGTACGTACGTACGT\n",
        "reverse\tsynthetic_transcript.1\t100\t20\t0\t0\t1\t20\t100\t81\t0.001\t40\t100\n",
    ];
    let pair = PrimerDesignPairRecord {
        forward: PrimerDesignPrimerRecord {
            sequence: "ACGTACGTACGTACGTACGT".into(),
            length_bp: 20,
            anneal_length_bp: 20,
            ..Default::default()
        },
        reverse: PrimerDesignPrimerRecord {
            sequence: "TGCATGCATGCATGCATGCA".into(),
            length_bp: 20,
            anneal_length_bp: 20,
            ..Default::default()
        },
        ..Default::default()
    };
    let primers = vec![
        GentleEngine::primer_specificity_input_from_record(
            PrimerSpecificityPrimerRole::Forward,
            &pair.forward,
        )
        .unwrap(),
        GentleEngine::primer_specificity_input_from_record(
            PrimerSpecificityPrimerRole::Reverse,
            &pair.reverse,
        )
        .unwrap(),
    ];
    let database = crate::genomes::BlastDatabaseInspectionReport {
        source_genome_id: "synthetic_cdna".into(),
        index_kind: BlastDatabaseIndexKind::TranscriptomeCdna,
        source_release: Some("synthetic_release".into()),
        content_fingerprint: Some(format!("sha256:{}", "a".repeat(64))),
        fingerprint_algorithm: "synthetic_sha256".into(),
        validation_status: "valid".into(),
        sequence_count: Some(2),
        ..Default::default()
    };
    let reference = GentleEngine::primer_specificity_reference_identity(&database).unwrap();
    let intended = PrimerSpecificityIntendedTarget {
        reference_binding: Some(reference.clone()),
        model: PrimerSpecificityIntendedTargetModel::TranscriptSet,
        expected_products: vec![PrimerSpecificityExpectedProduct {
            target_space: "transcriptome_cdna".into(),
            subject_id: "synthetic_transcript.1".into(),
            expected_product_range: Some(PrimerSpecificitySubjectRange {
                start_1based: 1,
                end_1based: 100,
            }),
            source_transcript_id: Some("synthetic_transcript.1".into()),
        }],
        source: "synthetic_fixture".into(),
        ..Default::default()
    };
    let complete = PrimerSpecificitySearchCompleteness {
        complete: true,
        status: "complete".into(),
        database_sequence_count: Some(2),
        required_max_target_seqs: Some(3),
        observed_min_max_target_seqs: Some(3),
        command_count: 2,
        reason: "synthetic subject-coverage evidence only".into(),
    };
    let mut child = PrimerSpecificityHandoff {
        schema: PRIMER_SPECIFICITY_HANDOFF_SCHEMA.into(),
        requested_target_genome_id: "synthetic_cdna".into(),
        resolved_target_genome_id: "synthetic_cdna".into(),
        policy: policy.clone(),
        primers: primers.clone(),
        blast_database: Some(database.clone()),
        intended_target: intended.clone(),
        search_completeness: complete.clone(),
        completion_policy: "all_commands_success".into(),
        commands: primers
            .iter()
            .map(|primer| PrimerSpecificityHandoffCommand {
                command_id: format!("synthetic:{}", primer.role.as_str()),
                role: primer.role,
                query_label: primer.role.as_str().into(),
                query_length_bp: 20,
                query_fasta_path: format!("synthetic/{}.fa", primer.role.as_str()),
                output_tsv_path: format!("synthetic/{}.tsv", primer.role.as_str()),
                program: "not_executed_blastn".into(),
                args: vec![
                    "-evalue".into(),
                    "1000".into(),
                    "-max_target_seqs".into(),
                    "3".into(),
                ],
                success_exit_codes: vec![0],
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    };
    child.handoff_id = GentleEngine::primer_specificity_handoff_id_from_record(&child).unwrap();
    let request = PrimerSpecificityMultiRequest {
        schema: PRIMER_SPECIFICITY_MULTI_REQUEST_SCHEMA.into(),
        pair: PrimerSpecificityMultiPair::ExplicitPair {
            forward: primers[0].clone(),
            reverse: primers[1].clone(),
        },
        policy: policy.clone(),
        references: vec![PrimerSpecificityMultiReference {
            target_genome_id: "synthetic_cdna".into(),
            expected_index_kind: BlastDatabaseIndexKind::TranscriptomeCdna,
            required: true,
            intended_target: Some(intended.clone()),
        }],
        catalog_path: None,
        cache_dir: None,
    };
    let mut multi = PrimerSpecificityMultiHandoff {
        schema: PRIMER_SPECIFICITY_MULTI_HANDOFF_SCHEMA.into(),
        handoff_id: String::new(),
        content_sha256: String::new(),
        request: request.clone(),
        pair_binding_sha256: GentleEngine::primer_specificity_pair_binding(
            &primers[0],
            &primers[1],
        )
        .unwrap(),
        source_snapshot_sha256: None,
        primers: primers.clone(),
        references: vec![PrimerSpecificityMultiHandoffReference {
            requested_genome_id: "synthetic_cdna".into(),
            resolved_genome_id: "synthetic_cdna".into(),
            required: true,
            expected_index_kind: BlastDatabaseIndexKind::TranscriptomeCdna,
            availability: PrimerSpecificityMultiAvailability::Prepared,
            diagnostic: None,
            reference: Some(reference.clone()),
            child: Some(child.clone()),
        }],
        manifest_path: "synthetic/execution_manifest.json".into(),
        nonclaims: MULTI_NONCLAIMS.iter().map(|s| s.to_string()).collect(),
    };
    multi.content_sha256 = GentleEngine::specificity_multi_handoff_digest(&multi).unwrap();
    multi.handoff_id = format!(
        "primer_multi_{}",
        multi.content_sha256.trim_start_matches("sha256:")
    );
    let alignment = PrimerSpecificityFullAlignment {
        status: "complete".into(),
        algorithm: "blast_full_query_aligned_strings_v1".into(),
        match_score: 2,
        mismatch_score: -3,
        gap_open_score: -5,
        gap_extend_score: -2,
        score: 40,
        subject_window_start_1based: 1,
        subject_window_end_1based: 20,
        aligned_subject_start_1based: 1,
        aligned_subject_end_1based: 20,
        strand: "+".into(),
        query_alignment: primers[0].annealing_sequence.clone(),
        relation_alignment: "|".repeat(20),
        subject_alignment: primers[0].annealing_sequence.clone(),
        cigar: "20M".into(),
        matches: 20,
        three_prime_window_bp: 5,
        subject_window_sha256: format!("sha256:{}", "b".repeat(64)),
        ..Default::default()
    };
    let forward = PrimerSpecificityPrimerHit {
        role: PrimerSpecificityPrimerRole::Forward,
        subject_id: "synthetic_transcript.1".into(),
        identity_percent: 100.0,
        alignment_length_bp: 20,
        query_start_1based: 1,
        query_end_1based: 20,
        subject_start_1based: 1,
        subject_end_1based: 20,
        subject_min_1based: 1,
        subject_max_1based: 20,
        strand: "+".into(),
        evalue: 0.001,
        bit_score: 40.0,
        query_coverage_fraction: 1.0,
        three_prime_window_bp: 5,
        binding_start_1based: Some(1),
        binding_end_1based: Some(20),
        full_alignment: Some(alignment),
        accepted_by_policy: true,
        ..Default::default()
    };
    let reverse = PrimerSpecificityPrimerHit {
        role: PrimerSpecificityPrimerRole::Reverse,
        subject_start_1based: 100,
        subject_end_1based: 81,
        subject_min_1based: 81,
        subject_max_1based: 100,
        strand: "-".into(),
        binding_start_1based: None,
        binding_end_1based: None,
        full_alignment: None,
        ..forward.clone()
    };
    let report = PrimerSpecificityReport {
        schema: PRIMER_SPECIFICITY_REPORT_SCHEMA.into(),
        report_id: "synthetic_saved_specificity_report".into(),
        target_kind: "transcriptome_cdna".into(),
        target_genome_id: "synthetic_cdna".into(),
        blast_database: Some(database),
        policy: policy.clone(),
        intended_target: intended,
        primers: primers.clone(),
        forward_hits: vec![forward],
        reverse_hits: vec![reverse],
        amplicons: vec![PrimerSpecificityAmplicon {
            kind: PrimerSpecificityAmpliconKind::ForwardReverse,
            subject_id: "synthetic_transcript.1".into(),
            left_role: PrimerSpecificityPrimerRole::Forward,
            right_role: PrimerSpecificityPrimerRole::Reverse,
            start_1based: 1,
            end_1based: 100,
            length_bp: 100,
            terminal_policy_pass: true,
            within_readiness_amplicon_range: true,
            intended: true,
            intended_reason: Some("synthetic_fixture".into()),
            ..Default::default()
        }],
        raw_detail_artifacts: outputs
            .map(|output| ComputationalArtifactExternalInput {
                checksum: Some(sha256_prefixed_bytes(output.as_bytes())),
                ..Default::default()
            })
            .into(),
        compaction: PrimerSpecificityReportCompaction {
            detail_mode: policy.report_detail_mode,
            raw_forward_hit_count: 1,
            retained_forward_hit_count: 1,
            raw_reverse_hit_count: 1,
            retained_reverse_hit_count: 1,
            raw_amplicon_count: 1,
            retained_amplicon_count: 1,
            pairing_candidate_comparison_count: 1,
        },
        search_completeness: complete,
        summary: PrimerSpecificitySummary {
            status: "pass".into(),
            specificity_pass: true,
            primer_hit_count: 2,
            accepted_primer_hit_count: 2,
            amplicon_count: 1,
            intended_amplicon_count: 1,
            ..Default::default()
        },
        transcriptome_specificity: PrimerSpecificityTargetAssessment {
            target_space: "transcriptome_cdna".into(),
            status: "pass".into(),
            intended_target_model: PrimerSpecificityIntendedTargetModel::TranscriptSet,
            intended_product_observed: true,
            expected_intended_product_count: 1,
            observed_intended_product_count: 1,
            compatible_product_count: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut summary = PrimerSpecificityMultiSummary {
        schema: PRIMER_SPECIFICITY_MULTI_SUMMARY_SCHEMA.into(),
        summary_id: String::new(),
        content_sha256: String::new(),
        handoff_id: multi.handoff_id.clone(),
        pair_binding_sha256: multi.pair_binding_sha256.clone(),
        source_snapshot_sha256: None,
        request,
        execution_manifest: PrimerSpecificityMultiExecutionManifest {
            schema: PRIMER_SPECIFICITY_MULTI_MANIFEST_SCHEMA.into(),
            handoff_id: multi.handoff_id.clone(),
            handoff_content_sha256: multi.content_sha256.clone(),
            pair_binding_sha256: multi.pair_binding_sha256.clone(),
            commands: child
                .commands
                .iter()
                .zip(outputs)
                .map(|(command, output)| PrimerSpecificityMultiExecutionCommand {
                    command_id: command.command_id.clone(),
                    command_sha256: GentleEngine::specificity_multi_command_digest(command)
                        .unwrap(),
                    output_path: command.output_tsv_path.clone(),
                    state: PrimerSpecificityMultiExecutionState::Completed,
                    exit_code: Some(0),
                    output_size_bytes: Some(output.len() as u64),
                    output_sha256: Some(sha256_prefixed_bytes(output.as_bytes())),
                })
                .collect(),
        },
        execution_complete: true,
        references: vec![PrimerSpecificityMultiSummaryReference {
            genome_id: "synthetic_cdna".into(),
            required: true,
            index_kind: BlastDatabaseIndexKind::TranscriptomeCdna,
            reference: Some(reference),
            verdict: PrimerSpecificityMultiVerdict::Pass,
            applicability: "synthetic_fixture_only".into(),
            diagnostics: vec![],
            report: Some(report.clone()),
        }],
        genomic: PrimerSpecificityMultiVerdict::NotRequested,
        transcriptome: PrimerSpecificityMultiVerdict::Pass,
        nonclaims: MULTI_NONCLAIMS.iter().map(|s| s.to_string()).collect(),
    };
    summary.content_sha256 = GentleEngine::specificity_multi_summary_digest(&summary).unwrap();
    summary.summary_id = format!(
        "primer_multi_summary_{}",
        summary.content_sha256.trim_start_matches("sha256:")
    );
    let panel = TranscriptAssayPanelReport {
        report_id: "synthetic_panel".into(),
        selected_assays: vec![TranscriptAssayPanelAssay {
            assay_id: "synthetic_assay".into(),
            rank: 1,
            primer_pair: pair,
            ..Default::default()
        }],
        ..Default::default()
    };
    let acceptance = TranscriptAssayPanelSpecificityAcceptance {
        acceptance_id: "synthetic_acceptance".into(),
        panel_report_id: panel.report_id.clone(),
        panel_digest: GentleEngine::transcript_assay_panel_specificity_digest(&panel).unwrap(),
        policy: policy.clone(),
        status: TranscriptAssayPanelSpecificityAcceptanceStatus::Pass,
        accepted: true,
        expected_assay_count: 1,
        assessed_assay_count: 1,
        assessments: vec![TranscriptAssayGenomicSpecificityAssessment {
            assay_id: "synthetic_assay".into(),
            assay_rank: 1,
            status: "external_blast_pass".into(),
            report: report.clone(),
            ..Default::default()
        }],
        ..Default::default()
    };
    let mut store = PrimerDesignStore::default();
    GentleEngine::record_primer_specificity_reference_selection(&mut store, &panel, &acceptance)
        .unwrap();
    store
        .primer_specificity_multi_summaries
        .insert(summary.summary_id.clone(), summary.clone());
    let selection_id = store
        .primer_specificity_reference_selections
        .keys()
        .next()
        .unwrap()
        .clone();
    json!({
        "origin": "synthetic records captured from b62d21c971c207704b96396523eec0ce61659806 plus test-only harness; no external execution or biological acceptance",
        "policy_sha256": sha256_prefixed_bytes(&serde_json::to_vec(&policy).unwrap()),
        "report_sha256": sha256_prefixed_bytes(&serde_json::to_vec(&report).unwrap()),
        "panel_digest": GentleEngine::primer_specificity_selection_panel_digest(&panel).unwrap(),
        "selection_id": selection_id,
        "synthetic_outputs": outputs,
        "policy": policy, "single_handoff": child, "multi_handoff": multi,
        "report": report, "summary": summary, "panel": panel, "store": store
    })
}
