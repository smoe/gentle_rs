//! Hand-written synthetic oligos/templates/FASTQ, generated only in temporary directories.
//! No order workbook, real gene, private read or kit sequence is used. These tests exercise
//! shared intake, matching, publication and inspection; they do not establish assay validity.

use super::*;

const HANDLE: &str = "TCAGAGCTGACA";
const CORE: &str = "ACGTCAGTACCGTTAGCCTAGG";
const PAYLOAD: &str = "GATCCATGAGCTAGTCGACATCGTAGCTACGATGCTAGCATCG";

fn primers(orientation: PrimerCoreOrientation) -> RnaPrimerPanel {
    let sequence = format!("{HANDLE}{CORE}");
    let mut panel = RnaPrimerPanel {
        schema: RNA_PRIMER_PANEL_SCHEMA.into(),
        content_sha256: String::new(),
        primers: vec![RnaLibraryPrimer {
            primer_id: "synthetic_a".into(),
            assay_label: "Artificial trace".into(),
            sequence: sequence.clone(),
            components: vec![
                RnaPrimerComponent {
                    component_id: "h".into(),
                    role: "common_handle".into(),
                    start: 0,
                    end: HANDLE.len(),
                },
                RnaPrimerComponent {
                    component_id: "c".into(),
                    role: "target_core".into(),
                    start: HANDLE.len(),
                    end: sequence.len(),
                },
            ],
            core_orientation: orientation,
            expected_targets: vec![RnaPrimerTarget {
                locus_id: "synthetic_locus".into(),
                transcript_feature_id: Some(0),
            }],
            variant: "invented".into(),
            source: "Hand-written synthetic regression, not a physical oligo".into(),
            source_sha256: None,
        }],
        policy: RnaPrimerTracePolicy::default(),
    };
    GentleEngine::seal_rna_primer_panel(&mut panel).unwrap();
    panel
}
fn reference(orientation: PrimerCoreOrientation) -> RnaReadPanel {
    let sequence = if orientation == PrimerCoreOrientation::Sense5prime {
        format!("GGATCA{CORE}{PAYLOAD}TGCCTA")
    } else {
        format!(
            "TGCCTA{}{}GGATCA",
            String::from_utf8(matching::rc(PAYLOAD.as_bytes())).unwrap(),
            String::from_utf8(matching::rc(CORE.as_bytes())).unwrap()
        )
    };
    let mut panel:RnaReadPanel=serde_json::from_value(serde_json::json!({"schema":RNA_PANEL_SCHEMA,"content_sha256":"",
        "loci":[{"source":{"locus_id":"synthetic_locus","seq_id":"synthetic_record","source_feature_id":0,"source":"synthetic","assembly":"synthetic","annotation_release":"synthetic","gene_id":"synthetic_gene","gene_symbol":"Artificial","selected_transcript_ids":[]},
            "source_record_sha256":sha256_prefixed_str("hand-written synthetic record"),"genome_anchor":null,"annotated_transcript_count":1,"annotated_transcript_ids":["synthetic_tx"]}],
        "templates":[{"member_id":"synthetic_locus:0","locus_id":"synthetic_locus","transcript_feature_id":0,"transcript_id":"synthetic_tx","strand":"+","sequence":sequence,
            "sequence_sha256":sha256_prefixed_str(&sequence),"junction_offsets":[],"source_exons":[[0,sequence.len()]]}],
        "gene_count":1,"locus_count":1,"selected_transcript_count":1,"annotated_transcript_count":1,"overlapping_source_pairs":[]})).unwrap();
    seal_panel(&mut panel).unwrap();
    validate_panel(&panel).unwrap();
    panel
}
fn record(sequence: &[u8]) -> input::FastqRecord {
    input::FastqRecord {
        read_id: "synthetic_read".into(),
        sequence: sequence.to_vec(),
        quality: vec![b'I'; sequence.len()],
        logical_offset: 0,
        compressed_bytes_read: 0,
        compressed_bytes_total: 0,
    }
}
fn classify(
    panel: &RnaPrimerPanel,
    reference: Option<&RnaReadPanel>,
    sequence: &[u8],
) -> RnaPrimerTraceRow {
    let b = RnaPrimerTraceBudgets::default();
    let sites = reference
        .map(|r| matching::sites(panel, r, &b, &mut || false).unwrap().0)
        .unwrap_or_default();
    matching::row(
        panel,
        reference,
        &sites,
        &record(sequence),
        0,
        &b,
        &mut || false,
    )
    .unwrap()
}
fn request(root: &Path, panel: RnaPrimerPanel) -> TracePrimerReadsRequest {
    let sequence = format!("GGGTTA{}{PAYLOAD}", panel.primers[0].sequence);
    let bytes = format!("@invented\n{sequence}\n+\n{}\n", "I".repeat(sequence.len()));
    let path = root.join("input with spaces.fastq");
    fs::write(&path, bytes.as_bytes()).unwrap();
    TracePrimerReadsRequest {
        schema: RNA_PRIMER_TRACE_REQUEST_SCHEMA.into(),
        primer_panel: RnaPrimerPanelBinding::Inline {
            panel: Box::new(panel),
        },
        reference_panel: None,
        files: vec![RnaReadPanelInput {
            sample_id: "sample".into(),
            barcode: "artificial_barcode".into(),
            priming: "synthetic".into(),
            path: path.display().to_string(),
            declared_complete: true,
            expected_md5: None,
            expected_sha256: Some(sha256_prefixed_bytes(bytes.as_bytes())),
            checksum_manifest_key: None,
        }],
        md5_manifest_path: None,
        output_dir: root.join("fresh output").display().to_string(),
        export_payload_fasta: false,
        budgets: RnaPrimerTraceBudgets::default(),
        cancel_path: None,
    }
}

#[test]
fn primer_trace_four_orientation_geometries_keep_the_boundary_adjacent_probe() {
    for orientation in [
        PrimerCoreOrientation::Sense5prime,
        PrimerCoreOrientation::Antisense3prime,
    ] {
        let panel = primers(orientation);
        let reference = reference(orientation);
        let forward = format!("GGGTTA{}{PAYLOAD}", panel.primers[0].sequence).into_bytes();
        for reverse in [false, true] {
            let read = if reverse {
                matching::rc(&forward)
            } else {
                forward.clone()
            };
            let row = classify(&panel, Some(&reference), &read);
            assert_eq!(row.anchor_class, "full_exact");
            assert_eq!(row.locus_relation, "expected_target_contiguous");
            let hit = row.anchor.as_ref().unwrap();
            assert_eq!(
                hit.strand,
                if reverse {
                    PrimerReadStrand::ReverseComplement
                } else {
                    PrimerReadStrand::Forward
                }
            );
            assert_eq!(
                row.payload,
                Some(if reverse {
                    (0, PAYLOAD.len())
                } else {
                    (6 + HANDLE.len() + CORE.len(), forward.len())
                })
            );
            assert_eq!(
                row.payload_transcript_sense,
                Some(orientation == PrimerCoreOrientation::Sense5prime)
            );
            assert!(!row.artifact_flags.iter().any(|f| f.contains("antisense")));
        }
    }
}

#[test]
fn primer_trace_substitutions_insertions_deletions_and_core_limit() {
    let panel = primers(PrimerCoreOrientation::Sense5prime);
    let reference = reference(PrimerCoreOrientation::Sense5prime);
    for edits in 1..=2 {
        let mut oligo = panel.primers[0].sequence.clone().into_bytes();
        for i in 0..edits {
            oligo[i] = if oligo[i] == b'A' { b'C' } else { b'A' };
        }
        let read = [b"GGGTTA".as_slice(), &oligo, PAYLOAD.as_bytes()].concat();
        let row = classify(&panel, Some(&reference), &read);
        assert_eq!(row.anchor_class, "full_approximate");
        assert_eq!(row.anchor.as_ref().unwrap().edits, edits);
        assert_eq!(row.locus_relation, "expected_target_contiguous");
    }
    for insertion in [false, true] {
        let mut oligo = panel.primers[0].sequence.clone().into_bytes();
        if insertion {
            oligo.insert(HANDLE.len() + 4, b'C');
        } else {
            oligo.remove(HANDLE.len() + 4);
        }
        let read = [b"GGGTTA".as_slice(), &oligo, PAYLOAD.as_bytes()].concat();
        let row = classify(&panel, Some(&reference), &read);
        assert_eq!(row.anchor_class, "full_approximate");
        assert_eq!(
            row.payload.unwrap().0,
            6 + oligo.len(),
            "observed span, not nominal oligo length"
        );
        assert_eq!(row.anchor.unwrap().core_edits, Some(1));
        assert_eq!(row.locus_relation, "expected_target_contiguous");
    }
    let mut read = format!("{}{}", panel.primers[0].sequence, PAYLOAD).into_bytes();
    for offset in [2, 7] {
        read[HANDLE.len() + offset] = b'A';
    }
    let row = classify(&panel, Some(&reference), &read);
    assert_ne!(row.locus_relation, "expected_target_contiguous");
    assert!(row.reasons.contains(&"core_edit_limit".into()));
}

#[test]
fn primer_trace_handle_no_anchor_short_and_quality_are_not_target_calls() {
    let panel = primers(PrimerCoreOrientation::Sense5prime);
    let row = classify(&panel, None, b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
    assert_eq!(row.anchor_class, "no_primer");
    assert!(row.anchor_position.is_none());
    assert!(row.primer_multiplicity.is_none());
    assert!(row.payload.is_none());
    for reverse in [false, true] {
        let sequence = format!("{HANDLE}TTTTTTTTTTTT{PAYLOAD}").into_bytes();
        let sequence = if reverse {
            matching::rc(&sequence)
        } else {
            sequence
        };
        let row = classify(&panel, None, &sequence);
        assert_eq!(row.anchor_class, "handle_only");
        assert!(row.anchor.unwrap().primer_id.is_none());
        assert!(
            row.artifact_flags
                .contains(&"handle_followed_by_poly_t".into())
        );
        assert!(row.payload.is_none());
    }
    let sequence = format!("{}AC", panel.primers[0].sequence);
    let mut r = record(sequence.as_bytes());
    let n = r.quality.len();
    r.quality[n - 2] = b'!';
    r.quality[n - 1] = b'I';
    let row = matching::row(
        &panel,
        None,
        &matching::ReferencePreparation::default(),
        &r,
        0,
        &RnaPrimerTraceBudgets::default(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(row.payload_mean_phred, Some(20.));
    assert_eq!(row.payload_fraction_below_q10, Some(0.5));
    assert!(row.artifact_flags.contains(&"short_payload".into()));
}

#[test]
fn primer_trace_multiple_anchors_do_not_become_unchecked_molecules() {
    let panel = primers(PrimerCoreOrientation::Sense5prime);
    let read = format!(
        "{}{}{}{}",
        panel.primers[0].sequence, PAYLOAD, panel.primers[0].sequence, PAYLOAD
    );
    let row = classify(&panel, None, read.as_bytes());
    assert_eq!(
        row.primer_multiplicity.as_deref(),
        Some("repeated_same_primer")
    );
    assert!(
        row.reasons
            .contains(&"multiple_anchor_payload_not_segmented".into())
    );
    let reverse = String::from_utf8(matching::rc(panel.primers[0].sequence.as_bytes())).unwrap();
    let read = format!("{}{}{}", panel.primers[0].sequence, PAYLOAD, reverse);
    assert_eq!(
        classify(&panel, None, read.as_bytes())
            .primer_multiplicity
            .as_deref(),
        Some("conflicting_strands")
    );
    let mut p = panel.clone();
    p.policy.end_window_nt = 64;
    GentleEngine::seal_rna_primer_panel(&mut p).unwrap();
    let read = format!("{}{}{PAYLOAD}", "G".repeat(100), p.primers[0].sequence);
    assert_eq!(
        classify(&p, None, read.as_bytes())
            .anchor_position
            .as_deref(),
        Some("internal")
    );
}

#[test]
fn primer_trace_reference_ambiguity_is_retained_and_bad_flanks_stay_unknown() {
    let panel = primers(PrimerCoreOrientation::Sense5prime);
    let mut reference = reference(PrimerCoreOrientation::Sense5prime);
    let mut other = reference.templates[0].clone();
    other.member_id = "synthetic_locus:1".into();
    other.transcript_feature_id = 1;
    reference.templates.push(other);
    let mut panel = panel;
    panel.primers[0].expected_targets[0].transcript_feature_id = None;
    GentleEngine::seal_rna_primer_panel(&mut panel).unwrap();
    let sequence = format!("{}{PAYLOAD}", panel.primers[0].sequence);
    let row = classify(&panel, Some(&reference), sequence.as_bytes());
    assert!(
        row.template_matches
            .iter()
            .any(|m| m.member_id == "synthetic_locus:0")
    );
    assert!(
        row.template_matches
            .iter()
            .any(|m| m.member_id == "synthetic_locus:1")
    );
    for t in &mut reference.templates {
        t.sequence = format!("{CORE}AC");
    }
    let row = classify(&panel, Some(&reference), sequence.as_bytes());
    assert_eq!(row.locus_relation, "unresolved");
    assert!(
        !row.artifact_flags
            .contains(&"primer_boundary_divergence".into())
    );
    let sequence = format!("{}{}", panel.primers[0].sequence, "N".repeat(45));
    assert!(
        classify(&panel, Some(&reference), sequence.as_bytes())
            .reasons
            .contains(&"ambiguous_payload_probe".into())
    );
}

#[test]
fn primer_trace_policy_binding_budgets_and_cancel_are_explicit() {
    let mut panel = primers(PrimerCoreOrientation::Sense5prime);
    panel.policy.max_full_edits = 3;
    assert!(validate(&panel).is_err());
    GentleEngine::seal_rna_primer_panel(&mut panel).unwrap();
    panel.primers[0].components[0].end += 1;
    assert!(GentleEngine::seal_rna_primer_panel(&mut panel).is_err());
    let root = tempfile::tempdir().unwrap();
    let mut request = request(root.path(), primers(PrimerCoreOrientation::Sense5prime));
    request.cancel_path = Some(root.path().join("CANCEL").display().to_string());
    GentleEngine::cancel_primer_trace(&request).unwrap();
    let engine = GentleEngine::new();
    let report = engine
        .trace_primer_reads(request.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(report.status, "cancelled");
    assert!(report.files_completed.is_empty());
    assert_eq!(report.inputs_not_completed.len(), 1);
    assert!(report.summary.by_sample_barcode.is_empty());
    let binding = GentleEngine::primer_trace_result_binding(&request, &report).unwrap();
    assert_eq!(binding.outcome, "cancelled");
}

#[test]
fn primer_trace_publication_inspection_tamper_and_digest_failure() {
    let root = tempfile::tempdir().unwrap();
    let mut request = request(root.path(), primers(PrimerCoreOrientation::Sense5prime));
    request.export_payload_fasta = true;
    let engine = GentleEngine::new();
    let before = serde_json::to_value(engine.state()).unwrap();
    let report = engine
        .trace_primer_reads(request.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(report.status, "completed");
    assert_eq!(report.files_completed[0].counts.reads, 1);
    assert_eq!(
        report.files_completed[0]
            .counts
            .anchor_start_histogram
            .get("6"),
        Some(&1)
    );
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    let binding = GentleEngine::primer_trace_result_binding(&request, &report).unwrap();
    let inspected = GentleEngine::inspect_primer_trace(
        &binding.report_path,
        Some(&binding),
        Some("file_0000"),
        1,
    )
    .unwrap();
    assert_eq!(inspected.rows[0].read_id, "invented");
    assert!(!inspected.artifact_content_verified);
    for artifact in &report.files_completed[0].artifacts {
        let bytes = fs::read(
            Path::new(&request.output_dir)
                .join("artifacts")
                .join(&artifact.path),
        )
        .unwrap();
        assert_eq!(artifact.sha256, sha256_prefixed_bytes(&bytes));
        assert_eq!(artifact.bytes, bytes.len() as u64);
    }
    let mut wrong = binding.clone();
    wrong.report_content_sha256 = sha256_prefixed_str("wrong");
    assert!(
        GentleEngine::inspect_primer_trace(&binding.report_path, Some(&wrong), None, 0).is_err()
    );
    request.output_dir = root.path().join("digest failure").display().to_string();
    request.files[0].expected_sha256 = Some(sha256_prefixed_str("wrong"));
    assert!(
        engine
            .trace_primer_reads(request.clone(), &mut |_| true)
            .is_err()
    );
    assert!(!Path::new(&request.output_dir).exists());
    request.files[0].expected_sha256 = Some(sha256_prefixed_bytes(
        &fs::read(&request.files[0].path).unwrap(),
    ));
    request.budgets.max_sidecar_bytes = 1;
    assert!(
        engine
            .trace_primer_reads(request.clone(), &mut |_| true)
            .is_err()
    );
    assert!(!Path::new(&request.output_dir).exists());
}

#[test]
fn primer_trace_shared_fastq_keeps_quality_and_raw_lf_crlf_hashes() {
    let root = tempfile::tempdir().unwrap();
    for newline in ["\n", "\r\n"] {
        let text = ["@synthetic", "ACGT", "+", "!I!I", ""].join(newline);
        let path = root.path().join("read.fastq");
        fs::write(&path, &text).unwrap();
        let digest = input::visit_fastq(&path, &mut |r| {
            assert_eq!(r.quality, b"!I!I");
            assert_eq!(r.sequence, b"ACGT");
            Ok(true)
        })
        .unwrap()
        .unwrap();
        assert_eq!(digest.sha256, sha256_prefixed_bytes(text.as_bytes()));
    }
}

#[test]
fn primer_trace_core_only_limits_and_primer_order_do_not_hide_competitors() {
    let mut panel = primers(PrimerCoreOrientation::Sense5prime);
    for errors in [1, 2] {
        let mut core = CORE.as_bytes().to_vec();
        core[2] = b'A';
        if errors == 2 {
            core[7] = b'A';
        }
        let row = classify(&panel, None, &[core, PAYLOAD.as_bytes().to_vec()].concat());
        if errors == 1 {
            assert_eq!(row.anchor_class, "core_approximate");
            assert_eq!(row.anchor.unwrap().core_edits, Some(1));
        } else {
            assert_eq!(row.anchor_class, "no_primer");
        }
    }
    let mut other = panel.primers[0].clone();
    other.primer_id = "synthetic_b".into();
    let mut oligo = other.sequence.into_bytes();
    oligo[HANDLE.len() + 2] = b'A';
    oligo[HANDLE.len() + 7] = b'A';
    other.sequence = String::from_utf8(oligo).unwrap();
    panel.primers.push(other);
    GentleEngine::seal_rna_primer_panel(&mut panel).unwrap();
    let mut read = panel.primers[0].sequence.as_bytes().to_vec();
    read[HANDLE.len() + 2] = b'A';
    read.extend_from_slice(PAYLOAD.as_bytes());
    let first = classify(&panel, None, &read);
    assert_eq!(
        first.anchor.as_ref().unwrap().primer_id.as_deref(),
        Some("synthetic_a")
    );
    assert_eq!(
        first.primer_multiplicity.as_deref(),
        Some("multiple_primers")
    );
    panel.primers.reverse();
    GentleEngine::seal_rna_primer_panel(&mut panel).unwrap();
    assert_eq!(
        encode(&first).unwrap(),
        encode(&classify(&panel, None, &read)).unwrap()
    );
}

#[test]
fn primer_trace_other_template_probe_is_not_contiguous_core_evidence() {
    let panel = primers(PrimerCoreOrientation::Sense5prime);
    let mut refs = reference(PrimerCoreOrientation::Sense5prime);
    refs.templates[0].sequence = format!("{CORE}{}", "C".repeat(45));
    let mut other = refs.templates[0].clone();
    other.member_id = "synthetic_locus:1".into();
    other.transcript_feature_id = 1;
    other.sequence = format!("{}{PAYLOAD}", "A".repeat(22));
    refs.templates.push(other);
    let sequence = format!("{}{PAYLOAD}", panel.primers[0].sequence);
    let row = classify(&panel, Some(&refs), sequence.as_bytes());
    assert_eq!(row.locus_relation, "other_panel_template");
    assert!(
        row.template_matches
            .iter()
            .any(|m| m.member_id == "synthetic_locus:1"
                && !m.core_and_probe_contiguous
                && m.core_3prime_mismatches == Some(8))
    );
    assert!(
        row.artifact_flags
            .contains(&"primer_boundary_divergence".into())
    );
    let mut changed = PAYLOAD.as_bytes().to_vec();
    changed[0] = b'C';
    let ordinary_reference = reference(PrimerCoreOrientation::Sense5prime);
    let read = [panel.primers[0].sequence.as_bytes(), changed.as_slice()].concat();
    let row = classify(&panel, Some(&ordinary_reference), &read);
    assert_eq!(row.locus_relation, "expected_target_contiguous");
    assert!(
        !row.artifact_flags
            .contains(&"primer_boundary_divergence".into())
    );
}

#[test]
fn primer_trace_bound_reference_does_not_change_primer_observation_fields() {
    let panel = primers(PrimerCoreOrientation::Antisense3prime);
    let reference = reference(PrimerCoreOrientation::Antisense3prime);
    let read = format!("GGGTTA{}{PAYLOAD}", panel.primers[0].sequence);
    let project = |row: RnaPrimerTraceRow| {
        (
            row.anchor_class,
            row.anchor_position,
            row.primer_multiplicity,
            row.anchor,
            row.hits,
            row.leader,
            row.payload,
            row.payload_transcript_sense,
            row.payload_mean_phred,
            row.payload_fraction_below_q10,
        )
    };
    assert_eq!(
        project(classify(&panel, None, read.as_bytes())),
        project(classify(&panel, Some(&reference), read.as_bytes()))
    );
}

#[test]
fn primer_trace_marker_fault_preserves_verified_units_and_reservation_fails_closed() {
    let root = tempfile::tempdir().unwrap();
    let mut request = request(root.path(), primers(PrimerCoreOrientation::Sense5prime));
    let marker = root.path().join("marker");
    request.cancel_path = Some(marker.display().to_string());
    let engine = GentleEngine::new();
    let report = engine
        .trace_primer_reads(request.clone(), &mut |_| {
            fs::write(&marker, "invalid marker").unwrap();
            true
        })
        .unwrap();
    assert_eq!(report.status, "control_error");
    assert_eq!(report.files_completed.len(), 1);
    assert_eq!(report.files_completed[0].counts.reads, 1);
    assert!(report.inputs_not_completed.is_empty());
    assert!(report.control_error.is_some());
    fs::remove_file(&marker).unwrap();
    request.cancel_path = None;
    request.output_dir = root
        .path()
        .join("occupied during publication")
        .display()
        .to_string();
    let output = PathBuf::from(&request.output_dir);
    let error = engine
        .trace_primer_reads(request.clone(), &mut |_| {
            if !output.exists() {
                fs::create_dir(&output).unwrap();
                fs::write(output.join("sentinel"), "owner data").unwrap();
            }
            true
        })
        .unwrap_err();
    assert!(error.message.contains("output_reserve"));
    assert_eq!(
        fs::read_to_string(output.join("sentinel")).unwrap(),
        "owner data"
    );
    assert!(!output.join("report.json").exists());
}

#[test]
fn primer_trace_honest_redigest_cannot_authorize_traversal_or_fabricated_counts() {
    let root = tempfile::tempdir().unwrap();
    let request = request(root.path(), primers(PrimerCoreOrientation::Sense5prime));
    let mut report = GentleEngine::new()
        .trace_primer_reads(request, &mut |_| true)
        .unwrap();
    report.files_completed[0].file_id = "..".into();
    report.content_sha256.clear();
    report.content_sha256 = digest(&report).unwrap();
    assert!(validate_report(&report).is_err());
    report.files_completed[0].file_id = "file_0000".into();
    report.files_completed[0].counts.reads += 1;
    report.content_sha256.clear();
    report.content_sha256 = digest(&report).unwrap();
    assert!(validate_report(&report).is_err());
}

#[test]
fn primer_trace_draft_sealing_is_explicit_and_never_repairs_a_bound_hash() {
    let panel = primers(PrimerCoreOrientation::Sense5prime);
    let expected = panel.content_sha256.clone();
    let mut draft = panel.clone();
    draft.content_sha256.clear();
    let binding = RnaPrimerPanelBinding::Inline {
        panel: Box::new(draft),
    };
    assert!(GentleEngine::validate_rna_primer_panel(&binding, false).is_err());
    let sealed = GentleEngine::validate_rna_primer_panel(&binding, true).unwrap();
    assert_eq!(sealed.content_sha256, expected);
    let mut changed = panel;
    changed.content_sha256 = sha256_prefixed_str("wrong");
    let binding = RnaPrimerPanelBinding::Inline {
        panel: Box::new(changed),
    };
    assert!(GentleEngine::validate_rna_primer_panel(&binding, false).is_err());
    assert!(GentleEngine::validate_rna_primer_panel(&binding, true).is_err());
}

#[test]
fn primer_trace_engine_operation_binds_reference_before_streaming_and_preserves_state() {
    let root = tempfile::tempdir().unwrap();
    let mut request = request(root.path(), primers(PrimerCoreOrientation::Sense5prime));
    let reference = reference(PrimerCoreOrientation::Sense5prime);
    let path = root.path().join("reference.json");
    fs::write(&path, serde_json::to_vec_pretty(&reference).unwrap()).unwrap();
    request.reference_panel = Some(RnaPrimerReferenceBinding {
        path: path.display().to_string(),
        expected_sha256: reference.content_sha256,
    });
    let mut engine = GentleEngine::new();
    let before = serde_json::to_value(engine.state()).unwrap();
    let result = engine
        .apply(Operation::TracePrimerReads {
            request: Box::new(request.clone()),
        })
        .unwrap();
    assert!(result.created_seq_ids.is_empty());
    assert!(result.changed_seq_ids.is_empty());
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    let report = result.rna_primer_trace.unwrap();
    let binding = GentleEngine::primer_trace_result_binding(&request, &report).unwrap();
    let inspection = GentleEngine::inspect_primer_trace(
        &binding.report_path,
        Some(&binding),
        Some("file_0000"),
        1,
    )
    .unwrap();
    assert_eq!(
        inspection.rows[0].locus_relation,
        "expected_target_contiguous"
    );
    assert!(
        inspection.rows[0]
            .template_matches
            .iter()
            .all(|m| m.core_and_probe_contiguous)
    );
    request.reference_panel.as_mut().unwrap().expected_sha256 = sha256_prefixed_str("wrong");
    request.output_dir = root.path().join("wrong binding").display().to_string();
    assert!(
        engine
            .trace_primer_reads(request.clone(), &mut |_| true)
            .is_err()
    );
    assert!(!Path::new(&request.output_dir).exists());
}

#[test]
fn primer_trace_three_full_edits_are_not_relabelled_as_a_full_match() {
    let panel = primers(PrimerCoreOrientation::Sense5prime);
    let mut oligo = panel.primers[0].sequence.as_bytes().to_vec();
    for i in [1, 5, 9] {
        oligo[i] = if oligo[i] == b'A' { b'C' } else { b'A' };
    }
    oligo.extend_from_slice(PAYLOAD.as_bytes());
    let row = classify(&panel, None, &oligo);
    assert!(!row.hits.iter().any(|h| h.anchor_class.starts_with("full_")));
    assert_eq!(
        row.anchor_class, "core_exact",
        "a retained lower-level observation is not a full hit"
    );
}

#[test]
fn primer_trace_malformed_panels_and_read_work_budgets_fail_without_publication() {
    let mut panel = primers(PrimerCoreOrientation::Sense5prime);
    panel.primers[0].sequence.replace_range(0..1, "N");
    assert!(GentleEngine::seal_rna_primer_panel(&mut panel).is_err());
    let mut panel = primers(PrimerCoreOrientation::Sense5prime);
    let mut duplicate = panel.primers[0].clone();
    duplicate.primer_id = "different_id".into();
    panel.primers.push(duplicate);
    assert!(GentleEngine::seal_rna_primer_panel(&mut panel).is_err());
    for read_limit in [true, false] {
        let root = tempfile::tempdir().unwrap();
        let mut request = request(root.path(), primers(PrimerCoreOrientation::Sense5prime));
        if read_limit {
            request.budgets.max_read_nt = 1;
        } else {
            request.budgets.max_scan_cells = 1;
        }
        assert!(
            GentleEngine::new()
                .trace_primer_reads(request.clone(), &mut |_| true)
                .is_err()
        );
        assert!(!Path::new(&request.output_dir).exists());
    }
}
