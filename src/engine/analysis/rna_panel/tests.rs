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
            expected_source_record_sha256: None,
        });
    }
    (
        GentleEngine::from_state(state),
        BuildRnaReadPanelRequest {
            gene_set: None,
            register_in_project: None,
            schema: RNA_PANEL_BUILD_SCHEMA.into(),
            sources,
            output_path: directory.join("panel.json").display().to_string(),
            base_panel_path: None,
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
        expected_panel_sha256: None,
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
fn rna_panel_cross_source_identical_dna_without_anchors_is_flagged_before_streaming() {
    let td = tempfile::tempdir().unwrap();
    let (mut engine, mut build) = starter(2, td.path());
    let sequence = engine.state().sequences["source_0"].get_forward_string();
    let mut copied = DNAsequence::from_sequence(&sequence).unwrap();
    copied
        .features_mut()
        .extend(engine.state().sequences["source_1"].features().clone());
    engine
        .state_mut()
        .sequences
        .insert("source_1".into(), copied);
    // Deliberately distinct source/gene IDs and feature annotations, with no genome anchors.
    build.sources[0].source = "synthetic_ensembl".into();
    build.sources[0].gene_id = "SYN_ENSEMBL_ID".into();
    build.sources[1].source = "synthetic_refseq".into();
    build.sources[1].gene_id = "SYN_REFSEQ_ID".into();
    let panel = engine.build_rna_read_panel(&build).unwrap();
    assert!(panel.loci.iter().all(|l| l.genome_anchor.is_none()));
    assert_ne!(
        panel.loci[0].source_record_sha256,
        panel.loci[1].source_record_sha256
    );
    assert_eq!(
        panel.loci[0].source_sequence_sha256,
        panel.loci[1].source_sequence_sha256
    );
    assert_eq!(
        panel.overlapping_source_pairs,
        vec![("locus_0".into(), "locus_1".into())]
    );
    let req = request(&panel, td.path(), &panel.templates[0].sequence);
    let mut progress = 0;
    let error = engine
        .screen_rna_read_panel(req.clone(), &mut |_| {
            progress += 1;
            true
        })
        .unwrap_err();
    assert!(error.message.starts_with("panel_overlapping_sources:"));
    assert_eq!(progress, 0);
    assert!(!Path::new(&req.output_dir).exists());
}

#[test]
fn rna_panel_legacy_locus_sequence_hash_is_absent_not_fabricated_and_stays_digest_bound() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(2, td.path());
    let mut panel = engine.build_rna_read_panel(&build).unwrap();
    for locus in &mut panel.loci {
        locus.source_sequence_sha256 = None;
    }
    seal_panel(&mut panel).unwrap();
    let bytes = encode(&panel).unwrap();
    assert!(
        !String::from_utf8(bytes.clone())
            .unwrap()
            .contains("source_sequence_sha256")
    );
    let legacy: RnaReadPanel = serde_json::from_slice(&bytes).unwrap();
    validate_panel(&legacy).unwrap();
    assert_eq!(encode(&legacy).unwrap(), bytes);
    let mut tampered = legacy;
    tampered.loci[0].source_sequence_sha256 = Some("sha256:invalid".into());
    seal_panel(&mut tampered).unwrap();
    assert!(validate_panel(&tampered).is_err());
}

#[test]
fn rna_panel_restricted_selection_can_hide_shared_competitors_but_never_claims_resolution() {
    let td = tempfile::tempdir().unwrap();
    let (engine, mut build) = starter(1, td.path());
    let full = engine.build_rna_read_panel(&build).unwrap();
    let shared = &full.templates[0].sequence[..40];
    let full_request = request(&full, td.path(), shared);
    let full_report = engine
        .screen_rna_read_panel(full_request, &mut |_| true)
        .unwrap();
    build.sources[0].selected_transcript_ids = vec!["reused_tx".into()];
    build.output_path = td.path().join("restricted.json").display().to_string();
    let restricted = engine.build_rna_read_panel(&build).unwrap();
    let restricted_inputs = td.path().join("restricted-inputs");
    fs::create_dir(&restricted_inputs).unwrap();
    let mut restricted_request = request(&restricted, &restricted_inputs, shared);
    restricted_request.output_dir = td.path().join("restricted-run").display().to_string();
    let report = engine
        .screen_rna_read_panel(restricted_request, &mut |_| true)
        .unwrap();
    for (report, expected_members) in [(&full_report, 2), (&report, 1)] {
        let key: Value = serde_json::from_str(
            report.files_completed[0]
                .equivalence_classes
                .keys()
                .next()
                .unwrap(),
        )
        .unwrap();
        assert_eq!(
            key["member_ids"].as_array().unwrap().len(),
            expected_members
        );
        assert!(!report.isoform_assignment_performed);
        assert!(
            report
                .nonclaims
                .iter()
                .any(|s| s.contains("removes unselected competitors"))
        );
    }
    assert_eq!(restricted.annotated_transcript_count, 2);
    assert_eq!(restricted.selected_transcript_count, 1);
}

#[test]
fn rna_panel_shared_isoforms_can_exceed_the_default_anchor_budget_without_publishing() {
    let td = tempfile::tempdir().unwrap();
    let (mut engine, build) = starter(1, td.path());
    let sequence = bases(981, 4096);
    let mut dna = DNAsequence::from_sequence(&sequence).unwrap();
    for number in 0..64 {
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: gb_io::seq::Location::simple_range(0, 3000 + number),
            qualifiers: vec![
                ("gene".into(), Some("SYN_0".into())),
                ("transcript_id".into(), Some(format!("shared_{number}"))),
            ],
        });
    }
    engine.state_mut().sequences.insert("source_0".into(), dna);
    let panel = engine.build_rna_read_panel(&build).unwrap();
    assert_eq!(panel.selected_transcript_count, 64);
    let mut req = request(&panel, td.path(), &sequence[..3000]);
    req.policy = RnaReadPanelPolicy::default();
    let error = engine
        .screen_rna_read_panel(req.clone(), &mut |_| true)
        .unwrap_err();
    assert!(error.message.contains("panel_anchor_budget:"));
    assert!(error.message.contains("limit 100000"));
    assert!(error.message.contains("file_id=file_0000"));
    assert!(!Path::new(&req.output_dir).exists());
}

#[test]
fn rna_panel_key_byte_budget_is_run_cumulative_and_overflow_safe_without_large_allocations() {
    // Synthetic debits represent new keys in four different files, not real read classifications.
    let mut used = 0;
    for _file in 0..4 {
        used = admit_run_key_bytes(used, 4_000_000).unwrap();
    }
    assert_eq!(used, 16_000_000);
    assert_eq!(admit_run_key_bytes(used, 1), None);
    assert_eq!(admit_run_key_bytes(usize::MAX, 1), None);
}

#[test]
fn rna_panel_whole_gzip_member_loss_is_caught_by_the_expected_transfer_digest() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let req_sequence = &panel.templates[0].sequence;
    let mut req = request(&panel, td.path(), req_sequence);
    let first = fs::read(&req.files[0].path).unwrap();
    let mut complete = first.clone();
    complete.extend(gzip(&fastq(req_sequence, b'I')));
    req.files[0].expected_md5 = Some(input::hex(&Md5::digest(&complete)));
    // A valid first gzip member remains; decompression itself cannot detect the lost second one.
    fs::write(&req.files[0].path, &first).unwrap();
    let error = engine
        .screen_rna_read_panel(req.clone(), &mut |_| true)
        .unwrap_err();
    assert!(error.message.starts_with("panel_checksum_mismatch:"));
    assert!(!Path::new(&req.output_dir).exists());
}

#[test]
fn rna_panel_extension_matches_full_rebuild_and_rejects_stale_sources() {
    // Synthetic source records generated by starter; no private reads or biological claims.
    let td = tempfile::tempdir().unwrap();
    let (mut engine, all) = starter(20, td.path());
    let mut first = all.clone();
    first.sources.truncate(10);
    let base = engine.build_rna_read_panel(&first).unwrap();
    let original = fs::read(&first.output_path).unwrap();
    let before = serde_json::to_value(engine.state()).unwrap();
    let mut extension = all.clone();
    extension.sources = all.sources[10..].to_vec();
    extension.base_panel_path = Some(first.output_path.clone());
    extension.output_path = td.path().join("extended panel.json").display().to_string();
    let command = crate::engine_shell::parse_shell_line(&format!(
        "rna-reads panel-build {}",
        crate::engine_shell::quote_shell_arg(&serde_json::to_string(&extension).unwrap())
    ))
    .unwrap();
    assert!(
        command.is_state_mutating(),
        "extension writes require approval"
    );
    let out = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
    let extended: RnaReadPanel =
        serde_json::from_value(out.output["result"]["rna_read_panel"].clone()).unwrap();
    let mut complete = all.clone();
    complete.output_path = td.path().join("full.json").display().to_string();
    let full = engine.build_rna_read_panel(&complete).unwrap();
    assert_eq!(
        serde_json::to_value(&extended.loci).unwrap(),
        serde_json::to_value(&full.loci).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&extended.templates).unwrap(),
        serde_json::to_value(&full.templates).unwrap()
    );
    assert_eq!(
        extended.parent_panel_sha256.as_deref(),
        Some(base.content_sha256.as_str())
    );
    assert_ne!(
        extended.content_sha256, full.content_sha256,
        "extension provenance is identity-bound"
    );
    assert_eq!(
        (extended.gene_count, extended.selected_transcript_count),
        (20, 40)
    );
    assert_ne!(base.content_sha256, extended.content_sha256);
    assert_eq!(fs::read(&first.output_path).unwrap(), original);
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    let legacy_request = serde_json::to_value(&all).unwrap();
    for key in ["base_panel_path", "gene_set", "register_in_project"] {
        assert!(!legacy_request.as_object().unwrap().contains_key(key));
    }
    assert!(
        legacy_request["sources"][0]
            .get("expected_source_record_sha256")
            .is_none()
    );
    let legacy_panel = serde_json::to_value(&base).unwrap();
    for key in ["gene_set_resolution", "parent_panel_sha256"] {
        assert!(legacy_panel.get(key).is_none());
    }

    let mut duplicate = extension.clone();
    duplicate.sources = first.sources[..1].to_vec();
    duplicate.output_path = td.path().join("duplicate.json").display().to_string();
    assert!(engine.build_rna_read_panel(&duplicate).is_err());
    assert!(!Path::new(&duplicate.output_path).exists());
    let mut crlf = extension.clone();
    let crlf_path = td.path().join("base CRLF.json");
    fs::write(
        &crlf_path,
        String::from_utf8(original.clone())
            .unwrap()
            .replace('\n', "\r\n"),
    )
    .unwrap();
    crlf.base_panel_path = Some(crlf_path.display().to_string());
    crlf.output_path = td.path().join("from CRLF.json").display().to_string();
    assert_eq!(
        engine.build_rna_read_panel(&crlf).unwrap().content_sha256,
        extended.content_sha256
    );

    engine
        .state_mut()
        .sequences
        .get_mut("source_0")
        .unwrap()
        .features_mut()[0]
        .qualifiers
        .push(("note".into(), Some("edited source".into())));
    extension.output_path = td.path().join("stale.json").display().to_string();
    assert!(
        engine
            .build_rna_read_panel(&extension)
            .unwrap_err()
            .message
            .contains("panel_base_source_changed")
    );
    assert!(!Path::new(&extension.output_path).exists());
    assert_eq!(fs::read(&first.output_path).unwrap(), original);
}

#[test]
fn rna_panel_source_picker_matches_read_only_shell_and_preserves_unknowns() {
    let td = tempfile::tempdir().unwrap();
    let (mut engine, _) = starter(20, td.path());
    let before = serde_json::to_value(&engine).unwrap();
    let choices = engine.rna_read_panel_source_choices().unwrap();
    assert_eq!(
        choices.len(),
        20,
        "one row per exact transcript group, not per isoform"
    );
    assert!(
        choices
            .iter()
            .all(|c| c.transcript_ids.len() == 2 && c.diagnostic.is_none())
    );
    assert!(choices.iter().all(|c| c.source.assembly.is_empty()
        && c.source.annotation_release.is_empty()
        && c.source.gene_id.is_empty()
        && c.source.selected_transcript_ids.is_empty()));
    let command =
        crate::engine_shell::parse_shell_line("rna-reads panel-build --list-sources").unwrap();
    assert!(!command.is_state_mutating());
    assert!(crate::command_execution::CommandExecutionService::manages(
        &command
    ));
    let actual = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
    assert!(!actual.state_changed);
    assert_eq!(
        actual.output["sources"],
        serde_json::to_value(&choices).unwrap()
    );
    assert_eq!(serde_json::to_value(&engine).unwrap(), before);
    let mut stale = choices[0].source.clone();
    stale.source = "synthetic".into();
    stale.assembly = "synthetic".into();
    stale.annotation_release = "v1".into();
    stale.gene_id = "SYN_0".into();
    let stale_path = td.path().join("stale-picker.json");
    let dna = engine.state_mut().sequences.get_mut("source_0").unwrap();
    dna.features_mut().clear();
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(0, 150),
        qualifiers: vec![("gene".into(), Some("NO_TRANSCRIPTS".into()))],
    });
    let request = BuildRnaReadPanelRequest {
        schema: RNA_PANEL_BUILD_SCHEMA.into(),
        sources: vec![stale],
        output_path: stale_path.display().to_string(),
        base_panel_path: None,
        gene_set: None,
        register_in_project: None,
    };
    assert!(
        engine
            .build_rna_read_panel(&request)
            .unwrap_err()
            .message
            .contains("panel_source_changed")
    );
    assert!(!stale_path.exists());
    assert!(
        engine
            .rna_read_panel_source_choices()
            .unwrap()
            .iter()
            .any(|c| c.source.seq_id == "source_0" && c.diagnostic.is_some())
    );
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
    let equivalence: RnaReadPanelNominationSetKey =
        serde_json::from_str(receipt.equivalence_classes.keys().next().unwrap()).unwrap();
    assert_eq!(equivalence.member_ids, vec!["locus_0:0", "locus_0:1"]);
    assert_eq!(equivalence.class, "single_locus_candidate");
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
        b"@r\n\n+\n\n",
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
fn rna_panel_same_locus_opposite_orientations_remain_candidates_not_a_fusion_call() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let template = &panel.templates[0];
    let rc = GentleEngine::reverse_complement_bytes(template.sequence.as_bytes());
    let read = [template.sequence.as_bytes(), rc.as_slice()].concat();
    let policy = RnaReadPanelPolicy {
        kmer_len: 8,
        min_chain_bases: 16,
        ..Default::default()
    };
    let (index, _) = build_index(&panel, &policy).unwrap();
    let (class, candidates) = nominate(&panel, &index, &policy, &read, &mut || false).unwrap();
    assert_eq!(class, "single_locus_candidate");
    let orientations = candidates
        .iter()
        .filter(|c| c.member_id == template.member_id)
        .map(|c| c.read_orientation.as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        orientations,
        BTreeSet::from(["as_read", "reverse_complement"])
    );
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
    let mut run = request(&panel, td.path(), &joined);
    run.policy = policy.clone();
    let report = engine
        .screen_rna_read_panel(run.clone(), &mut |_| true)
        .unwrap();
    let key: RnaReadPanelNominationSetKey = serde_json::from_str(
        report.files_completed[0]
            .equivalence_classes
            .keys()
            .next()
            .unwrap(),
    )
    .unwrap();
    assert_eq!(key.class, "multi_segment_candidate");
    let summary = engine
        .export_rna_read_panel_summary(&ExportRnaReadPanelSummaryRequest {
            schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
            report_paths: vec![
                Path::new(&run.output_dir)
                    .join("report.json")
                    .display()
                    .to_string(),
            ],
            output_dir: td.path().join("segment summary").display().to_string(),
        })
        .unwrap();
    let set = &summary.co_nomination_sets.as_ref().unwrap()[0];
    assert_eq!(set.class, key.class);
    assert_eq!(set.member_ids, key.member_ids);
    assert_eq!(set.reads, 1);
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
    let junction_read = inspected["rows"][1]["candidates"].as_array().unwrap();
    let a = junction_read
        .iter()
        .find(|c| c["member_id"] == "locus_0:0" && c["read_orientation"] == "as_read")
        .unwrap();
    let b = junction_read
        .iter()
        .find(|c| c["member_id"] == "locus_0:1" && c["read_orientation"] == "as_read")
        .unwrap();
    assert_eq!(a["supported_junction_offsets"], serde_json::json!([40]));
    assert_eq!(b["supported_junction_offsets"], serde_json::json!([]));
    // demo_1 spans A's junction, but shared-exon seeds still nominate B. This is not whole-read fit.
    let verify = format!(
        "rna-reads panel-show {} --verify-artifacts",
        quote_shell_arg(&run.join("summary/report.json").display().to_string())
    );
    let command = parse_shell_line(&verify).unwrap();
    assert!(!command.is_state_mutating());
    let verified = execute_shell_command(&mut engine, &command).unwrap();
    assert_eq!(
        verified.output["result"]["rna_read_panel_verification"]["original_inputs_reverified"],
        false
    );
    assert_eq!(
        verified.output["result"]["rna_read_panel_verification"]["artifacts"]
            .as_array()
            .unwrap()
            .len(),
        3
    );
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

#[test]
fn rna_panel_cancellation_is_request_bound_and_marker_checks_are_throttled() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let mut run = request(&panel, td.path(), &panel.templates[0].sequence);
    let marker = td.path().join("CANCEL");
    run.cancel_path = Some(marker.display().to_string());
    engine.cancel_rna_read_panel(&run).unwrap();
    let bytes = fs::read(&marker).unwrap();
    engine.cancel_rna_read_panel(&run).unwrap();
    assert_eq!(fs::read(&marker).unwrap(), bytes);
    let mut other = run.clone();
    other.output_dir = td.path().join("other").display().to_string();
    assert!(
        engine
            .screen_rna_read_panel(other.clone(), &mut |_| true)
            .unwrap_err()
            .message
            .contains("marker_mismatch")
    );
    assert!(!Path::new(&other.output_dir).exists());
    assert!(engine.cancel_rna_read_panel(&other).is_err());
    assert_eq!(fs::read(&marker).unwrap(), bytes);
    let cancelled = engine.screen_rna_read_panel(run, &mut |_| true).unwrap();
    assert_eq!(cancelled.status, "cancelled");
    assert!(
        cancelled.files_completed.is_empty(),
        "matching pre-start cancellation is supported"
    );
    assert_eq!(cancelled.inputs_not_completed.as_ref().unwrap().len(), 1);
    let start = std::time::Instant::now();
    let mut gate = control::PollGate::default();
    let mut checks = 0;
    for micros in 0..1_000_000 {
        if gate.due(start + std::time::Duration::from_micros(micros), false) {
            checks += 1;
        }
    }
    assert_eq!(
        checks, 10,
        "marker I/O and progress only follow admitted polls, never each k-mer"
    );
    assert!(gate.due(start + std::time::Duration::from_secs(1), true));
    fs::write(&marker, b"unrelated old marker").unwrap();
    assert!(engine.cancel_rna_read_panel(&other).is_err());
}

#[test]
fn rna_panel_budget_errors_identify_the_record_and_publish_no_partial_results() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let mut run = request(&panel, td.path(), &panel.templates[0].sequence);
    run.policy.max_anchors_per_read = 1;
    let output = run.output_dir.clone();
    let error = engine
        .screen_rna_read_panel(run, &mut |_| true)
        .unwrap_err();
    for part in [
        "budget",
        "file_id=file_0000",
        "record_index=0",
        "read_id=read",
        "no result directory published",
    ] {
        assert!(error.message.contains(part), "{}: {part}", error.message);
    }
    assert!(!Path::new(&output).exists());
}

#[test]
fn rna_panel_mid_run_marker_errors_preserve_verified_units_with_explicit_outcomes() {
    use RnaReadPanelControlErrorPhase::{FileBoundary, Finalization, Streaming};
    // All faults are injected by deterministic progress callbacks, never sleeps or private data.
    for (trigger, fault, expected_phase, verified) in [
        ("boundary", "foreign", Some(FileBoundary), 1),
        ("streaming", "foreign", Some(Streaming), 1),
        ("partial", "foreign", Some(Streaming), 1),
        ("finalization", "foreign", Some(Finalization), 1),
        ("empty", "foreign", Some(Streaming), 0),
        ("streaming", "malformed", Some(Streaming), 1),
        ("boundary", "directory", Some(FileBoundary), 1),
        ("streaming", "matching", None, 1),
    ] {
        let td = tempfile::tempdir().unwrap();
        let (engine, build) = starter(1, td.path());
        let panel = engine.build_rna_read_panel(&build).unwrap();
        let mut run = request(&panel, td.path(), &panel.templates[0].sequence);
        run.policy.write_locus_fastas = true;
        let marker = td.path().join("CANCEL");
        run.cancel_path = Some(marker.display().to_string());
        if trigger != "finalization" {
            let mut second = run.files[0].clone();
            second.sample_id = "sample_b".into();
            second.barcode = "barcode02".into();
            second.path = td.path().join("second.fastq.gz").display().to_string();
            let mut raw = fastq(&panel.templates[0].sequence, b'I');
            if trigger == "partial" {
                raw.extend(fastq(&panel.templates[0].sequence, b'I'));
            }
            let bytes = gzip(&raw);
            fs::write(&second.path, &bytes).unwrap();
            second.expected_md5 = Some(input::hex(&Md5::digest(&bytes)));
            run.files.push(second);
        }
        let before = encode(engine.state()).unwrap();
        let mut injected = false;
        let mut callback = |event: OperationProgress| {
            let inject_here = match &event {
                OperationProgress::RnaReadPanel(p) => {
                    (matches!(trigger, "boundary" | "finalization")
                        && p.file_id == "file_0000"
                        && p.phase == "file_checksum_verified")
                        || (trigger == "partial"
                            && p.file_id == "file_0001"
                            && p.phase == "streaming_unverified"
                            && p.processed_reads == 1)
                }
                OperationProgress::Workflow { completed, .. } => {
                    (trigger == "streaming" && *completed == 1)
                        || (trigger == "empty" && *completed == 0)
                }
                _ => false,
            };
            if !injected && inject_here {
                match fault {
                    "foreign" => fs::write(&marker, format!("sha256:{}", "f".repeat(64))).unwrap(),
                    "malformed" => fs::write(&marker, [0xff]).unwrap(),
                    "directory" => fs::create_dir(&marker).unwrap(),
                    "matching" => {
                        engine.cancel_rna_read_panel(&run).unwrap();
                    }
                    _ => unreachable!(),
                }
                injected = true;
                return true;
            }
            // A later callback cancellation must not relabel an already-observed control error.
            !(injected && matches!(event, OperationProgress::Workflow { .. }))
        };
        let report = if trigger == "partial" {
            let start = std::time::Instant::now();
            let mut tick = 0;
            screen_with_clock(run.clone(), &mut callback, || {
                tick += 1;
                start + std::time::Duration::from_millis(tick * 100)
            })
        } else {
            engine.screen_rna_read_panel(run.clone(), &mut callback)
        }
        .unwrap();
        assert!(injected);
        assert_eq!(report.files_completed.len(), verified, "{trigger}/{fault}");
        assert_eq!(report.files_not_completed.len(), run.files.len() - verified);
        assert_eq!(
            report.inputs_not_completed.as_ref().unwrap().len(),
            run.files.len() - verified
        );
        if let Some(phase) = expected_phase {
            assert_eq!(report.status, "control_error");
            let error = report.control_error.as_ref().unwrap();
            assert_eq!(error.phase, phase);
            assert_eq!(error.code, "panel_control_marker_error");
            assert_eq!(error.marker_path, marker.display().to_string());
            assert!(error.message.starts_with("panel_cancel_marker"));
            assert_eq!(
                error.file_id,
                if phase == Finalization {
                    None
                } else {
                    Some(format!("file_{verified:04}"))
                }
            );
        } else {
            assert_eq!(report.status, "cancelled");
            assert!(report.control_error.is_none());
        }
        assert_eq!(encode(engine.state()).unwrap(), before);
        let output = Path::new(&run.output_dir);
        for file in &report.files_completed {
            for artifact in &file.artifacts {
                let bytes = fs::read(output.join(&artifact.path)).unwrap();
                assert_eq!(sha256_prefixed_bytes(&bytes), artifact.sha256);
                assert_eq!(bytes.len() as u64, artifact.bytes);
            }
        }
        assert!(
            !output
                .join(format!("artifacts/file_{verified:04}"))
                .exists(),
            "unfinished units must be discarded"
        );
        let path = output.join("report.json");
        let shown = engine
            .inspect_rna_read_panel(&path.to_string_lossy(), None, 0)
            .unwrap()
            .1
            .unwrap();
        assert_eq!(encode(&shown).unwrap(), encode(&report).unwrap());
        let mut inspector = engine.clone_without_history();
        let journal = inspector.journal_len();
        let command = crate::engine_shell::parse_shell_line(&format!(
            "rna-reads panel-show {}",
            crate::engine_shell::quote_shell_arg(&path.to_string_lossy())
        ))
        .unwrap();
        assert!(!command.is_state_mutating());
        let inspection =
            crate::engine_shell::execute_shell_command(&mut inspector, &command).unwrap();
        assert_eq!(
            inspection.output["result"]["rna_read_panel_screen"],
            serde_json::to_value(&report).unwrap()
        );
        assert_eq!(
            crate::agent_feedback::AgentExecutionStatus::from_shell_output(&inspection.output),
            crate::agent_feedback::AgentExecutionStatus::Completed,
            "inspection completed; this does not change the saved analysis outcome"
        );
        assert_eq!(inspector.journal_len(), journal);
        let feedback = serde_json::json!({"result":{"rna_read_panel_screen": &report}});
        assert_eq!(
            crate::agent_feedback::AgentExecutionStatus::from_shell_output(&feedback),
            if expected_phase.is_some() {
                crate::agent_feedback::AgentExecutionStatus::Failed
            } else {
                crate::agent_feedback::AgentExecutionStatus::Cancelled
            }
        );
        let summary = engine
            .export_rna_read_panel_summary(&ExportRnaReadPanelSummaryRequest {
                schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
                report_paths: vec![path.display().to_string()],
                output_dir: td.path().join("summary").display().to_string(),
            })
            .unwrap();
        assert_eq!(summary.completeness, "partial");
        assert_eq!(summary.source_runs[0].status, report.status);
        assert_eq!(
            encode(&summary.source_runs[0].control_error).unwrap(),
            encode(&report.control_error).unwrap()
        );
        assert_eq!(summary.source_files.len(), verified);
        assert_eq!(
            summary.samples.iter().map(|s| s.reads).sum::<u64>(),
            verified as u64
        );
        let coverage = summary.sample_coverage.as_ref().unwrap();
        if verified == 0 {
            assert!(
                coverage
                    .iter()
                    .all(|s| s.status == "unassessed" && s.reads.is_none())
            );
        }
        if run.files.len() > verified {
            assert_eq!(coverage.last().unwrap().reads, None);
        }
    }
}

#[test]
fn rna_panel_control_repair_keeps_checksum_and_resource_failures_hard() {
    for failure in ["checksum", "resource"] {
        let td = tempfile::tempdir().unwrap();
        let (engine, build) = starter(1, td.path());
        let panel = engine.build_rna_read_panel(&build).unwrap();
        let mut run = request(&panel, td.path(), &panel.templates[0].sequence[..12]);
        let mut second = run.files[0].clone();
        second.path = td.path().join("second.fastq.gz").display().to_string();
        let bytes = gzip(&fastq(&panel.templates[0].sequence, b'I'));
        fs::write(&second.path, &bytes).unwrap();
        second.expected_md5 = Some(if failure == "checksum" {
            "0".repeat(32)
        } else {
            input::hex(&Md5::digest(&bytes))
        });
        run.files.push(second);
        if failure == "resource" {
            run.policy.max_anchors_per_read = 64;
        }
        let mut verified_first = false;
        let error = engine
            .screen_rna_read_panel(run.clone(), &mut |event| {
                if let OperationProgress::RnaReadPanel(p) = event {
                    verified_first |=
                        p.file_id == "file_0000" && p.phase == "file_checksum_verified";
                }
                true
            })
            .unwrap_err();
        assert!(verified_first);
        assert!(
            error.message.contains(if failure == "checksum" {
                "checksum_mismatch"
            } else {
                "anchor_budget"
            }),
            "{}",
            error.message
        );
        assert!(
            !Path::new(&run.output_dir).exists(),
            "only control errors have the new retention policy"
        );
    }
}

#[test]
fn rna_panel_control_error_binding_cannot_be_relabelled_as_success_or_legacy() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let run = request(&panel, td.path(), &panel.templates[0].sequence);
    let report = engine.screen_rna_read_panel(run, &mut |_| true).unwrap();
    assert!(
        serde_json::to_value(&report)
            .unwrap()
            .get("control_error")
            .is_none(),
        "historical digest encoding stays unchanged when absent"
    );
    let mut invalid = report.clone();
    invalid.request.cancel_path = Some("synthetic_marker".into());
    invalid.control_error = Some(control::diagnostic(
        bad("panel_cancel_marker_mismatch: synthetic"),
        Some("synthetic_marker"),
        RnaReadPanelControlErrorPhase::Finalization,
        None,
    ));
    let seal = |r: &mut RnaReadPanelScreenReport| {
        r.content_sha256.clear();
        r.content_sha256 = sha256_prefixed_bytes(&encode(r).unwrap());
    };
    for status in ["completed", "cancelled"] {
        invalid.status = status.into();
        seal(&mut invalid);
        assert!(validate_screen_report_binding(&invalid).is_err());
    }
    invalid.status = "control_error".into();
    seal(&mut invalid);
    validate_screen_report_binding(&invalid).unwrap();
    invalid.schema = RNA_PANEL_SCREEN_LEGACY_SCHEMA.into();
    seal(&mut invalid);
    assert!(validate_screen_report_binding(&invalid).is_err());
    invalid.schema = RNA_PANEL_SCREEN_SCHEMA.into();
    invalid.control_error.as_mut().unwrap().marker_path = "another_request_marker".into();
    seal(&mut invalid);
    assert!(validate_screen_report_binding(&invalid).is_err());
    invalid.control_error = None;
    seal(&mut invalid);
    assert!(validate_screen_report_binding(&invalid).is_err());
}

#[test]
fn rna_panel_progress_separates_decoder_read_ahead_from_verified_file_units() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let run = request(&panel, td.path(), &panel.templates[0].sequence);
    let mut observations = vec![];
    let report = engine
        .screen_rna_read_panel(run, &mut |event| {
            if let OperationProgress::RnaReadPanel(progress) = event {
                observations.push(progress);
            }
            true
        })
        .unwrap();
    assert!(observations.len() >= 3);
    assert!(
        observations
            .windows(2)
            .all(|p| p[0].compressed_bytes_read <= p[1].compressed_bytes_read
                && p[0].processed_reads <= p[1].processed_reads)
    );
    assert!(
        observations[..observations.len() - 1]
            .iter()
            .all(|p| p.phase == "streaming_unverified")
    );
    let last = observations.last().unwrap();
    assert_eq!(last.phase, "file_checksum_verified");
    assert_eq!(last.processed_reads, report.files_completed[0].reads);
    assert_eq!(
        last.compressed_bytes_read,
        report.files_completed[0].compressed_bytes
    );
    assert_eq!(last.compressed_bytes_total, last.compressed_bytes_read);
}
