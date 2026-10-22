//! Inline synthetic loci/FASTQ, recreated from fixed arithmetic and literal bytes.
//! No human gene sequences, private downloads, basecaller or biological acceptance.

use super::*;
use flate2::{Compression, write::GzEncoder};
use md5::{Digest, Md5};
use serde_json::Value;

pub(super) fn bases(seed: u64, length: usize) -> String {
    let mut value = seed;
    (0..length)
        .map(|_| {
            value = value.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(value >> 32) as usize % 4] as char
        })
        .collect()
}

pub(super) fn starter(count: usize, directory: &Path) -> (GentleEngine, BuildRnaReadPanelRequest) {
    let mut state = ProjectState::default();
    let mut sources = vec![];
    for i in 0..count {
        let seq_id = format!("source_{i}");
        let gene = format!("SYN_{i}");
        let mut dna = DNAsequence::from_sequence(&bases(100 + i as u64, 150)).unwrap();
        for (id, end) in [("reused_tx", 100), ("reused_alt", 150)] {
            dna.features_mut().push(gb_io::seq::Feature {
                kind: "mRNA".into(),
                location: gb_io::seq::Location::Join(vec![
                    gb_io::seq::Location::simple_range(0, 40),
                    gb_io::seq::Location::simple_range(end - 40, end),
                ]),
                qualifiers: vec![
                    ("gene".into(), Some(gene.clone())),
                    ("transcript_id".into(), Some(id.into())),
                ],
            });
        }
        state.sequences.insert(seq_id.clone(), dna);
        sources.push(RnaReadPanelSource {
            locus_id: format!("locus_{i}"),
            seq_id,
            source_feature_id: 0,
            source: "hand_crafted".into(),
            assembly: "synthetic".into(),
            annotation_release: "v1".into(),
            gene_id: gene.clone(),
            gene_symbol: gene,
            selected_transcript_ids: vec![],
        });
    }
    (
        GentleEngine::from_state(state),
        BuildRnaReadPanelRequest {
            schema: RNA_PANEL_BUILD_SCHEMA.into(),
            sources,
            output_path: directory.join("panel.json").display().to_string(),
        },
    )
}

pub(super) fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut writer = GzEncoder::new(Vec::new(), Compression::fast());
    writer.write_all(bytes).unwrap();
    writer.finish().unwrap()
}
pub(super) fn fastq(sequence: &str, quality: u8) -> Vec<u8> {
    format!(
        "@read same-header-may-recur\n{sequence}\n+\n{}\n",
        (quality as char).to_string().repeat(sequence.len())
    )
    .into_bytes()
}
pub(super) fn request(
    panel: &RnaReadPanel,
    directory: &Path,
    sequence: &str,
) -> ScreenRnaReadPanelRequest {
    let path = directory.join("reads.fastq.gz");
    let bytes = gzip(&fastq(sequence, b'>'));
    fs::write(&path, &bytes).unwrap();
    write_json(&directory.join("panel-input.json"), panel).unwrap();
    ScreenRnaReadPanelRequest {
        schema: RNA_PANEL_SCREEN_REQUEST_SCHEMA.into(),
        panel_path: directory.join("panel-input.json").display().to_string(),
        files: vec![RnaReadPanelInput {
            sample_id: "sample_a".into(),
            barcode: "barcode01".into(),
            priming: "synthetic_custom".into(),
            path: path.display().to_string(),
            declared_complete: true,
            expected_md5: Some(input::hex(&Md5::digest(&bytes))),
            expected_sha256: None,
            checksum_manifest_key: None,
        }],
        md5_manifest_path: None,
        output_dir: directory.join("run").display().to_string(),
        policy: RnaReadPanelPolicy {
            kmer_len: 8,
            min_chain_bases: 16,
            min_junction_flank_bp: 3,
            ..Default::default()
        },
        cancel_path: None,
    }
}

#[test]
fn rna_panel_independent_sources_are_namespaced_and_build_is_read_only() {
    let td = tempfile::tempdir().unwrap();
    let (mut engine, req) = starter(5, td.path());
    let before = serde_json::to_value(engine.state()).unwrap();
    let panel = engine.build_rna_read_panel(&req).unwrap();
    assert_eq!(
        (
            panel.locus_count,
            panel.gene_count,
            panel.selected_transcript_count,
            panel.annotated_transcript_count
        ),
        (5, 5, 10, 10)
    );
    assert_eq!(
        panel
            .templates
            .iter()
            .map(|t| &t.member_id)
            .collect::<BTreeSet<_>>()
            .len(),
        10
    );
    assert_eq!(
        panel
            .templates
            .iter()
            .map(|t| &t.transcript_id)
            .collect::<BTreeSet<_>>()
            .len(),
        2
    );
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    let mut req2 = req.clone();
    req2.output_path = td.path().join("second.json").display().to_string();
    use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};
    let line = format!(
        "rna-reads panel-build {}",
        quote_shell_arg(&serde_json::to_string(&req2).unwrap())
    );
    let parsed = parse_shell_line(&line).unwrap();
    assert!(
        parsed.is_state_mutating(),
        "artifact writes require agent confirmation"
    );
    let out = execute_shell_command(&mut engine, &parsed).unwrap();
    assert!(!out.state_changed);
    assert_eq!(
        out.output["result"]["rna_read_panel"]["content_sha256"],
        panel.content_sha256
    );
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    assert!(engine.build_rna_read_panel(&req).is_err(), "no overwrite");
}

#[test]
fn rna_panel_shared_prefix_retains_isoforms_and_roundtrips_reports_and_artifacts() {
    let td = tempfile::tempdir().unwrap();
    let (engine, req) = starter(2, td.path());
    let panel = engine.build_rna_read_panel(&req).unwrap();
    let mut run = request(&panel, td.path(), &panel.templates[0].sequence[..40]);
    run.policy.write_locus_fastas = true;
    let before = serde_json::to_value(engine.state()).unwrap();
    let report = engine
        .screen_rna_read_panel(run.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(report.status, "completed");
    assert_eq!(report.files_completed[0].reads, 1);
    let receipt = &report.files_completed[0];
    assert_eq!(receipt.classes.get("single_locus_candidate"), Some(&1));
    let equivalence: Vec<String> =
        serde_json::from_str(receipt.equivalence_classes.keys().next().unwrap()).unwrap();
    assert_eq!(equivalence, vec!["locus_0:0", "locus_0:1"]);
    for artifact in &receipt.artifacts {
        let bytes = fs::read(Path::new(&run.output_dir).join(&artifact.path)).unwrap();
        assert_eq!(artifact.sha256, sha256_prefixed_bytes(&bytes));
        assert_eq!(artifact.bytes, bytes.len() as u64);
        if artifact.locus_id.as_deref() == Some("locus_0") {
            let fasta = String::from_utf8(bytes).unwrap();
            assert!(fasta.contains("ambiguous=true"));
            assert!(fasta.contains(&receipt.observed_md5));
        }
    }
    let retained: RnaReadPanelScreenReport = load(
        &Path::new(&run.output_dir)
            .join("report.json")
            .display()
            .to_string(),
    )
    .unwrap();
    assert_eq!(encode(&retained).unwrap(), encode(&report).unwrap());
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    assert!(engine.screen_rna_read_panel(run, &mut |_| true).is_err());
}

#[test]
fn rna_panel_twenty_loci_two_hundred_shared_isoforms_have_no_sixty_four_posting_bias() {
    let td = tempfile::tempdir().unwrap();
    let (engine, req) = starter(20, td.path());
    let mut panel = engine.build_rna_read_panel(&req).unwrap();
    let prototype = panel.templates[0].clone();
    panel.templates.clear();
    for locus in &mut panel.loci {
        locus.annotated_transcript_count = 10;
        locus.annotated_transcript_ids = vec![prototype.transcript_id.clone(); 10];
        for i in 0..10 {
            let mut t = prototype.clone();
            t.locus_id = locus.source.locus_id.clone();
            t.member_id = format!("{}:{i}", t.locus_id);
            t.transcript_feature_id = i;
            panel.templates.push(t);
        }
    }
    panel.selected_transcript_count = 200;
    panel.annotated_transcript_count = 200;
    seal_panel(&mut panel).unwrap();
    validate_panel(&panel).unwrap();
    let policy = RnaReadPanelPolicy {
        kmer_len: 8,
        min_chain_bases: 16,
        min_junction_flank_bp: 3,
        ..Default::default()
    };
    let (index, _) = build_index(&panel, &policy).unwrap();
    assert!(index.values().any(|postings| postings.len() >= 200));
    let (class, candidates) = nominate(
        &panel,
        &index,
        &policy,
        prototype.sequence.as_bytes(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(class, "multi_locus_ambiguous");
    assert_eq!(
        candidates
            .iter()
            .map(|c| &c.member_id)
            .collect::<BTreeSet<_>>()
            .len(),
        200
    );
    panel.templates.reverse();
    let (reversed, _) = build_index(&panel, &policy).unwrap();
    let (_, again) = nominate(
        &panel,
        &reversed,
        &policy,
        prototype.sequence.as_bytes(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(encode(&candidates).unwrap(), encode(&again).unwrap());
    let mut tiny = policy.clone();
    tiny.max_index_postings = 1;
    assert!(build_index(&panel, &tiny).is_err());
    tiny = policy.clone();
    tiny.max_anchors_per_read = 1;
    assert!(
        nominate(
            &panel,
            &index,
            &tiny,
            prototype.sequence.as_bytes(),
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn rna_panel_junction_flanks_and_both_orientations_are_explicit() {
    let td = tempfile::tempdir().unwrap();
    let (engine, req) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&req).unwrap();
    let mut policy = RnaReadPanelPolicy {
        kmer_len: 8,
        min_chain_bases: 8,
        min_junction_flank_bp: 3,
        ..Default::default()
    };
    let (index, _) = build_index(&panel, &policy).unwrap();
    let t = &panel.templates[0];
    let j = t.junction_offsets[0];
    let (_, weak) = nominate(
        &panel,
        &index,
        &policy,
        &t.sequence.as_bytes()[j - 2..j + 6],
        &mut || false,
    )
    .unwrap();
    assert!(
        weak.iter()
            .filter(|c| c.member_id == t.member_id)
            .all(|c| c.supported_junction_offsets.is_empty())
    );
    policy.min_chain_bases = 16;
    let rc = GentleEngine::reverse_complement_bytes(t.sequence.as_bytes());
    let (_, strong) = nominate(&panel, &index, &policy, &rc, &mut || false).unwrap();
    assert!(strong.iter().any(|c| c.member_id == t.member_id
        && c.read_orientation == "reverse_complement"
        && c.supported_junction_offsets == vec![j]));
    assert!(
        nominate(&panel, &index, &policy, &rc, &mut || true)
            .unwrap_err()
            .message
            .contains("cancelled")
    );
}

#[test]
fn rna_panel_fastq_validates_crlf_quality_markers_lengths_and_gzip_members() {
    let td = tempfile::tempdir().unwrap();
    let path = td.path().join("no-extension");
    let raw = fastq("ACGTACGT", b'>');
    let mut members = gzip(&raw);
    members.extend(gzip(
        &String::from_utf8(raw.clone())
            .unwrap()
            .replace('\n', "\r\n")
            .into_bytes(),
    ));
    fs::write(&path, &members).unwrap();
    let mut offsets = vec![];
    let digest = visit_fastq(&path, &mut |record| {
        offsets.push(record.logical_offset);
        Ok(true)
    })
    .unwrap()
    .unwrap();
    assert_eq!(offsets, vec![0, raw.len() as u64]);
    assert_eq!(digest.md5, input::hex(&Md5::digest(&members)));
    assert_eq!(digest.sha256, sha256_prefixed_bytes(&members));
    for invalid in [
        b">fasta\nACGT\n".as_slice(),
        b"@r\nACGT\n-\nIIII\n",
        b"@r\nACGT\n+\nIII\n",
        b"@r\nACGT\n+\n",
        b"@r\nACGT\n+x\nIIII\n",
    ] {
        fs::write(&path, invalid).unwrap();
        assert!(visit_fastq(&path, &mut |_| Ok(true)).is_err());
    }
    members.truncate(members.len() - 5);
    fs::write(&path, &members).unwrap();
    assert!(
        visit_fastq(&path, &mut |_| Ok(true)).is_err(),
        "truncated gzip must not pass"
    );
    fs::write(&path, &raw).unwrap();
    assert!(visit_fastq(&path, &mut |_| Ok(false)).unwrap().is_none());
}

#[test]
fn rna_panel_detects_growing_input_and_bad_transfer_digest_without_publication() {
    let td = tempfile::tempdir().unwrap();
    let path = td.path().join("growing.fastq");
    let raw = fastq("ACGTACGT", b'I');
    fs::write(&path, &raw).unwrap();
    let mut changed = false;
    assert!(
        visit_fastq(&path, &mut |_| {
            if !changed {
                fs::OpenOptions::new()
                    .append(true)
                    .open(&path)
                    .unwrap()
                    .write_all(&raw)
                    .unwrap();
                changed = true;
            }
            Ok(true)
        })
        .unwrap_err()
        .message
        .contains("changed_during_read")
    );
    let (engine, req) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&req).unwrap();
    let mut run = request(&panel, td.path(), &panel.templates[0].sequence);
    let mut sha_only = run.clone();
    sha_only.files[0].expected_md5 = None;
    sha_only.files[0].expected_sha256 = Some(sha256_prefixed_bytes(
        &fs::read(&run.files[0].path).unwrap(),
    ));
    run.files[0].expected_md5 = Some("0".repeat(32));
    let output = run.output_dir.clone();
    assert!(
        engine
            .screen_rna_read_panel(run, &mut |_| true)
            .unwrap_err()
            .message
            .contains("checksum_mismatch")
    );
    assert!(!Path::new(&output).exists());
    assert_eq!(
        engine
            .screen_rna_read_panel(sha_only, &mut |_| true)
            .unwrap()
            .status,
        "completed"
    );
}

#[test]
fn rna_panel_supplier_manifest_is_rename_tolerant_but_not_ambiguous() {
    let td = tempfile::tempdir().unwrap();
    let (engine, req) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&req).unwrap();
    let mut run = request(&panel, td.path(), &panel.templates[0].sequence);
    let manifest = td.path().join("md5_checksums.txt");
    let hash = run.files[0].expected_md5.take().unwrap();
    fs::write(
        &manifest,
        format!("{hash}  /supplier/original/barcode01/reads.fastq.gz\r\n"),
    )
    .unwrap();
    run.md5_manifest_path = Some(manifest.display().to_string());
    assert_eq!(
        expected_md5(
            &run.files[0],
            Some(&load_md5_manifest(&manifest.display().to_string()).unwrap())
        )
        .unwrap(),
        Some(hash.clone())
    );
    fs::write(
        &manifest,
        format!("{hash}  /a/reads.fastq.gz\n{hash}  /b/reads.fastq.gz\n"),
    )
    .unwrap();
    let error = validate_request(&mut run).unwrap_err();
    assert!(error.message.contains("/a/reads.fastq.gz"));
    assert!(error.message.contains("/b/reads.fastq.gz"));
    run.files[0].checksum_manifest_key = Some("/a/reads.fastq.gz".into());
    validate_request(&mut run).unwrap();
}

#[test]
fn rna_panel_cancelled_prefix_can_continue_at_file_boundaries_without_false_success() {
    let td = tempfile::tempdir().unwrap();
    let (engine, req) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&req).unwrap();
    let mut run = request(&panel, td.path(), &panel.templates[0].sequence);
    let mut second = run.files[0].clone();
    second.path = td.path().join("second.fastq.gz").display().to_string();
    fs::copy(&run.files[0].path, &second.path).unwrap();
    run.files.push(second);
    let report = engine
        .screen_rna_read_panel(
            run.clone(),
            &mut |p| !matches!(p,OperationProgress::Workflow{completed,..} if completed==1),
        )
        .unwrap();
    assert_eq!(report.status, "cancelled");
    assert_eq!(report.files_completed.len(), 1);
    assert_eq!(report.files_not_completed.len(), 1);
    let feedback = serde_json::json!({"result":{"rna_read_panel_screen":report}});
    assert_eq!(
        crate::agent_feedback::AgentExecutionStatus::from_shell_output(&feedback),
        crate::agent_feedback::AgentExecutionStatus::Cancelled
    );
    run.files.remove(0);
    run.output_dir = td.path().join("continuation").display().to_string();
    let continued = engine.screen_rna_read_panel(run, &mut |_| true).unwrap();
    assert_eq!(continued.status, "completed");
    assert_eq!(
        continued.files_completed[0].equivalence_classes,
        report.files_completed[0].equivalence_classes
    );
}

#[test]
fn rna_panel_binding_tampering_and_missing_selection_fail_closed() {
    let td = tempfile::tempdir().unwrap();
    let (engine, mut req) = starter(1, td.path());
    let mut panel = engine.build_rna_read_panel(&req).unwrap();
    panel.templates[0].sequence.replace_range(0..1, "N");
    assert!(validate_panel(&panel).is_err());
    panel.templates[0].sequence_sha256 = sha256_prefixed_str(&panel.templates[0].sequence);
    panel.templates[0].junction_offsets[0] += 1;
    seal_panel(&mut panel).unwrap();
    assert!(validate_panel(&panel).is_err());
    req.output_path = td.path().join("missing.json").display().to_string();
    req.sources[0].selected_transcript_ids = vec!["NOT_PRESENT".into()];
    assert!(engine.build_rna_read_panel(&req).is_err());
    assert!(
        serde_json::from_str::<RnaReadPanelPolicy>(r#"{"silently_discard_postings":true}"#)
            .is_err()
    );
}

#[test]
fn rna_panel_fastq_cannot_be_misread_as_legacy_fasta() {
    let td = tempfile::tempdir().unwrap();
    for (name, bytes) in [
        ("reads.fasta", fastq("ACGTACGT", b'>')),
        ("reads.fasta.gz", gzip(&fastq("ACGTACGT", b'>'))),
    ] {
        let path = td.path().join(name);
        fs::write(&path, bytes).unwrap();
        let mut records = 0;
        let result = GentleEngine::visit_fasta_records_with_offsets(
            &path.display().to_string(),
            &mut |_, _| {
                records += 1;
                Ok(())
            },
        );
        assert!(result.unwrap_err().message.contains("FASTA requires"));
        assert_eq!(records, 0, "Phred > is quality, not a FASTA header");
    }
}

#[test]
fn rna_panel_source_overlap_is_disclosed_and_screening_refuses_double_counting() {
    let td = tempfile::tempdir().unwrap();
    let (engine, mut req) = starter(1, td.path());
    let mut second = req.sources[0].clone();
    second.locus_id = "second_view".into();
    second.source_feature_id = 1;
    req.sources.push(second);
    let panel = engine.build_rna_read_panel(&req).unwrap();
    assert_eq!(panel.overlapping_source_pairs.len(), 1);
    let run = request(&panel, td.path(), &panel.templates[0].sequence);
    let out = run.output_dir.clone();
    assert!(
        engine
            .screen_rna_read_panel(run, &mut |_| true)
            .unwrap_err()
            .message
            .contains("overlapping_sources")
    );
    assert!(!Path::new(&out).exists());
}

#[test]
fn rna_panel_reverse_locus_and_disjoint_multi_locus_segments_remain_candidates() {
    let td = tempfile::tempdir().unwrap();
    let (mut engine, req) = starter(2, td.path());
    for feature in engine
        .state
        .sequences
        .get_mut("source_1")
        .unwrap()
        .features_mut()
    {
        feature.location = gb_io::seq::Location::Complement(Box::new(feature.location.clone()));
    }
    let panel = engine.build_rna_read_panel(&req).unwrap();
    let policy = RnaReadPanelPolicy {
        kmer_len: 12,
        min_chain_bases: 24,
        min_junction_flank_bp: 4,
        ..Default::default()
    };
    let (index, _) = build_index(&panel, &policy).unwrap();
    let a = panel
        .templates
        .iter()
        .find(|t| t.locus_id == "locus_0")
        .unwrap();
    let b = panel
        .templates
        .iter()
        .find(|t| t.locus_id == "locus_1")
        .unwrap();
    assert_eq!(b.strand, "-");
    let joined = format!("{}{}{}", a.sequence, "N".repeat(40), b.sequence);
    let (class, candidates) =
        nominate(&panel, &index, &policy, joined.as_bytes(), &mut || false).unwrap();
    assert_eq!(class, "multi_segment_candidate");
    assert!(
        candidates
            .iter()
            .any(|c| c.member_id == b.member_id && c.read_orientation == "as_read")
    );
    // Hand-introduced substitution and insertion need not improve a score;
    // they must not turn a deterministic seed screen into an exclusive call.
    let mut noisy = a.sequence.clone();
    noisy.replace_range(16..17, "N");
    noisy.insert_str(32, "NN");
    let first = nominate(&panel, &index, &policy, noisy.as_bytes(), &mut || false).unwrap();
    let second = nominate(&panel, &index, &policy, noisy.as_bytes(), &mut || false).unwrap();
    assert_eq!(encode(&first).unwrap(), encode(&second).unwrap());
    assert!(first.1.iter().any(|c| c.member_id == a.member_id));
}

#[test]
fn rna_panel_tutorial_is_executable_through_shared_agent_commands() {
    use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};
    let td = tempfile::tempdir().unwrap();
    let run = td.path().join("tutorial with spaces");
    let mut engine = GentleEngine::new();
    let before = serde_json::to_value(engine.state()).unwrap();
    let demo = format!(
        "rna-reads panel-demo {}",
        quote_shell_arg(&run.display().to_string())
    );
    let command = parse_shell_line(&demo).unwrap();
    assert!(
        command.is_state_mutating(),
        "external writes require approval"
    );
    execute_shell_command(&mut engine, &command).unwrap();
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    let directory = run.join("artifacts");
    let state: ProjectState = load(&directory.join("starter.json").display().to_string()).unwrap();
    engine = GentleEngine::from_state(state);
    let before = serde_json::to_value(engine.state()).unwrap();
    for (command, name) in [
        ("panel-build", "build.json"),
        ("panel-screen", "screen.json"),
        ("panel-export-summary", "summary.json"),
    ] {
        let line = format!(
            "rna-reads {command} {}",
            quote_shell_arg(&format!("@{}", directory.join(name).display()))
        );
        let out = execute_shell_command(&mut engine, &parse_shell_line(&line).unwrap()).unwrap();
        assert!(!out.state_changed);
    }
    let report_path = run.join("screen/report.json").display().to_string();
    let line = format!(
        "rna-reads panel-show {} --file-id file_0000 --rows 4",
        quote_shell_arg(&report_path)
    );
    assert!(!parse_shell_line(&line).unwrap().is_state_mutating());
    let out = execute_shell_command(&mut engine, &parse_shell_line(&line).unwrap()).unwrap();
    let report = &out.output["result"]["rna_read_panel_screen"];
    assert_eq!(report["status"], "completed");
    assert_eq!(report["files_completed"][0]["reads"], 4);
    assert_eq!(report["isoform_assignment_performed"], false);
    assert_eq!(report["alignment_confirmation_performed"], false);
    let inspected = &out.output["result"]["rna_read_panel_inspection"];
    assert_eq!(inspected["artifact_content_verified"], false);
    assert_eq!(
        inspected["rows"][0]["candidates"].as_array().unwrap().len(),
        2
    );
    assert_eq!(inspected["rows"][3]["class"], "unmatched");
    let mut request: ScreenRnaReadPanelRequest =
        load(&directory.join("screen.json").display().to_string()).unwrap();
    request.output_dir = run.join("counterfactual").display().to_string();
    let direct = engine
        .screen_rna_read_panel(request, &mut |_| true)
        .unwrap();
    assert_eq!(
        report["files_completed"],
        serde_json::to_value(direct.files_completed).unwrap()
    );
    let cancel = format!(
        "rna-reads panel-cancel {}",
        quote_shell_arg(&format!("@{}", directory.join("screen.json").display()))
    );
    let command = parse_shell_line(&cancel).unwrap();
    assert!(command.is_state_mutating());
    execute_shell_command(&mut engine, &command).unwrap();
    let marker = fs::read(run.join("CANCEL")).unwrap();
    execute_shell_command(&mut engine, &command).unwrap();
    assert_eq!(fs::read(run.join("CANCEL")).unwrap(), marker);
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    let mut altered: Value = load(&report_path).unwrap();
    altered["status"] = "cancelled".into();
    fs::write(&report_path, serde_json::to_vec(&altered).unwrap()).unwrap();
    assert!(
        engine
            .inspect_rna_read_panel(&report_path, None, 0)
            .is_err()
    );
}

#[test]
fn rna_panel_named_controls_do_not_invalidate_a_running_detached_screen() {
    use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};
    let td = tempfile::tempdir().unwrap();
    let (mut engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let mut run = request(&panel, td.path(), &panel.templates[0].sequence);
    run.cancel_path = Some(td.path().join("CANCEL").display().to_string());
    let mut detached = engine.fork_detached_execution();
    let before = engine.journal_len();
    let show = format!(
        "rna-reads panel-show {}",
        quote_shell_arg(&build.output_path)
    );
    execute_shell_command(&mut engine, &parse_shell_line(&show).unwrap()).unwrap();
    let cancel = format!(
        "rna-reads panel-cancel {}",
        quote_shell_arg(&serde_json::to_string(&run).unwrap())
    );
    let out = execute_shell_command(&mut engine, &parse_shell_line(&cancel).unwrap()).unwrap();
    assert_eq!(out.output["cancel_requested"], true);
    assert_eq!(out.output["termination_confirmed"], false);
    assert_eq!(engine.journal_len(), before);
    let report_path = Path::new(&run.output_dir).join("report.json");
    let result = detached
        .engine_mut()
        .apply(Operation::ScreenRnaReadPanel {
            request: Box::new(run),
        })
        .unwrap();
    assert_eq!(result.rna_read_panel_screen.unwrap().status, "cancelled");
    engine.commit_detached_execution(&mut detached).unwrap();
    let show = format!(
        "rna-reads panel-show {}",
        quote_shell_arg(&report_path.display().to_string())
    );
    let out = execute_shell_command(&mut engine, &parse_shell_line(&show).unwrap()).unwrap();
    assert_eq!(
        out.output["result"]["rna_read_panel_screen"]["status"],
        "cancelled"
    );
    assert_eq!(
        crate::agent_feedback::AgentExecutionStatus::from_shell_output(&out.output),
        crate::agent_feedback::AgentExecutionStatus::Completed,
        "inspection completed, analysis cancelled"
    );
}
