//! Synthetic deterministic work-profile regressions for the RNA panel engine.
//!
//! Provenance: all loci and reads are hand-crafted literals or recreated with
//! `tests::{starter, bases, fastq, gzip}` and fixed seeds. Temporary sidecars
//! exercise nomination/report contracts only; no private or biological data,
//! external reference, throughput claim or scientific acceptance is involved.

use super::tests::{bases, fastq, gzip, request, starter};
use super::*;
use serde_json::json;
use std::time::Instant;

const EXPECTED_ANCHOR_BOUNDS: [u64; 20] = [
    0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072,
    262144,
];

#[test]
fn rna_panel_work_total_overflow_is_unknown_not_a_new_resource_failure() {
    // Synthetic near-limit counter state, not a huge file or changed read budget.
    let mut profile = work_profile::file_profile();
    profile.total_anchor_postings = Some(u64::MAX);
    for record in 0..2 {
        work_profile::record_read(
            &mut profile,
            &work_profile::ReadWork {
                anchors: 1,
                template_chain_evaluations: 1,
            },
            record,
            "synthetic-overflow-counter",
        )
        .unwrap();
        assert_eq!(profile.total_anchor_postings, None);
    }
    assert_eq!(profile.anchor_histogram_counts[1], 2);
}

#[test]
fn rna_panel_work_profiles_reject_rehashed_counter_and_peak_inconsistency() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let run = request(&panel, td.path(), &panel.templates[0].sequence);
    let report = screen(run, &mut |_| true).unwrap();
    validate_screen_report_binding(&report).unwrap();
    for change in 0..4 {
        let mut changed = report.clone();
        if change == 0 {
            changed
                .index_work_profile
                .as_mut()
                .unwrap()
                .words_above_posting_thresholds
                .remove(&8);
        } else {
            let profile = changed.files_completed[0].work_profile.as_mut().unwrap();
            match change {
                1 => profile.anchor_histogram_counts[0] += 1,
                2 => profile.peak_anchor_postings.as_mut().unwrap().read_id = "x".repeat(257),
                3 => profile.equivalence_key_bytes_added += 1,
                _ => unreachable!(),
            }
        }
        changed.content_sha256.clear();
        changed.content_sha256 = sha256_prefixed_bytes(&encode(&changed).unwrap());
        assert!(
            validate_screen_report_binding(&changed)
                .unwrap_err()
                .message
                .contains("work_profile_binding")
        );
    }
}

fn replace_synthetic_input(run: &mut ScreenRnaReadPanelRequest, bytes: &[u8]) {
    fs::write(&run.files[0].path, bytes).unwrap();
    run.files[0].expected_md5 = None;
    run.files[0].expected_sha256 = Some(sha256_prefixed_bytes(bytes));
}

#[test]
fn rna_panel_work_key_bytes_debit_repeated_keys_independently_for_each_verified_file() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let sequence = &panel.templates[0].sequence;
    let mut run = request(&panel, td.path(), sequence);
    let bytes = gzip(&fastq(sequence, b'I'));
    let path = td.path().join("second.fastq.gz");
    fs::write(&path, &bytes).unwrap();
    let mut second = run.files[0].clone();
    second.path = path.display().to_string();
    second.expected_md5 = None;
    second.expected_sha256 = Some(sha256_prefixed_bytes(&bytes));
    run.files.push(second);
    let report = screen(run, &mut |_| true).unwrap();
    validate_screen_report_binding(&report).unwrap();
    assert_eq!(report.files_completed.len(), 2);
    assert_eq!(
        report.files_completed[0].equivalence_classes,
        report.files_completed[1].equivalence_classes
    );
    let key_bytes = report.files_completed[0]
        .equivalence_classes
        .keys()
        .map(String::len)
        .sum::<usize>();
    assert!(key_bytes > 0);
    assert!(
        report.files_completed.iter().all(|f| f
            .work_profile
            .as_ref()
            .unwrap()
            .equivalence_key_bytes_added
            == key_bytes)
    );
    assert_eq!(
        report
            .files_completed
            .iter()
            .map(|f| f.work_profile.as_ref().unwrap().equivalence_key_bytes_added)
            .sum::<usize>(),
        2 * key_bytes
    );
}

#[test]
fn rna_panel_index_work_fanout_counts_all_four_strict_thresholds() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(2, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    assert_eq!(panel.templates.len(), 4);
    assert_eq!(panel.templates[0].locus_id, panel.templates[1].locus_id);
    assert_ne!(panel.templates[0].locus_id, panel.templates[2].locus_id);

    // One synthetic word at and immediately above each threshold. Positions
    // repeat within each template; fanout counts distinct templates/loci only.
    let mut index = Index::new();
    for (word, count) in [1, 2, 8, 9, 64, 65, 256, 257].into_iter().enumerate() {
        index.insert(
            word as u32,
            (0..count).map(|posting| (posting % 4, 0)).collect(),
        );
    }
    assert_eq!(
        work_profile::index_profile(&panel, &index),
        RnaReadPanelIndexWorkProfile {
            distinct_words: 8,
            max_postings_per_word: 257,
            max_templates_per_word: 4,
            max_loci_per_word: 2,
            words_above_posting_thresholds: BTreeMap::from([(1, 7), (8, 5), (64, 3), (256, 1)]),
        }
    );
}

#[test]
fn rna_panel_file_work_histogram_totals_and_peak_ties_are_hand_computed() {
    let mut profile = work_profile::file_profile();
    let anchors = [0, 1, 2, 3, 4, 5, 8, 9, 250_000, 250_000];
    for (record, anchors) in anchors.into_iter().enumerate() {
        let read_id = if record == 8 { "z-first" } else { "a-later" };
        work_profile::record_read(
            &mut profile,
            &work_profile::ReadWork {
                anchors,
                template_chain_evaluations: record,
            },
            record as u64,
            read_id,
        )
        .unwrap();
    }
    assert_eq!(
        profile,
        RnaReadPanelFileWorkProfile {
            peak_anchor_postings: Some(RnaReadPanelReadWorkPeak {
                anchors: 250_000,
                record_index: 8,
                read_id: "z-first".into(),
                read_id_truncated: false,
                read_id_sha256: None,
            }),
            anchor_histogram_upper_bounds: EXPECTED_ANCHOR_BOUNDS.to_vec(),
            anchor_histogram_counts: vec![
                1, 1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2
            ],
            total_anchor_postings: Some(500_032),
            peak_template_chain_evaluations: 9,
            equivalence_key_bytes_added: 0,
        }
    );
    assert_eq!(profile.anchor_histogram_counts.iter().sum::<u64>(), 10);
}

#[test]
fn rna_panel_peak_read_id_is_bounded_by_unicode_scalars_with_full_id_hash() {
    let mut profile = work_profile::file_profile();
    let prefix = format!("{}\u{1d11e}", "\u{03b1}".repeat(255));
    assert_eq!(prefix.chars().count(), 256);
    assert!(
        prefix.len() > 256,
        "the boundary must not count UTF-8 bytes"
    );
    work_profile::record_read(
        &mut profile,
        &work_profile::ReadWork {
            anchors: 1,
            template_chain_evaluations: 1,
        },
        0,
        &prefix,
    )
    .unwrap();
    let exact = profile.peak_anchor_postings.as_ref().unwrap();
    assert_eq!(exact.read_id, prefix);
    assert!(!exact.read_id_truncated);
    assert_eq!(exact.read_id_sha256, None);

    let full_id = format!("{prefix}synthetic-tail");
    work_profile::record_read(
        &mut profile,
        &work_profile::ReadWork {
            anchors: 2,
            template_chain_evaluations: 1,
        },
        1,
        &full_id,
    )
    .unwrap();
    let peak = profile.peak_anchor_postings.as_ref().unwrap().clone();
    assert_eq!(peak.read_id, prefix);
    assert_eq!(peak.read_id.chars().count(), 256);
    assert!(peak.read_id_truncated);
    assert_eq!(
        peak.read_id_sha256,
        Some(sha256_prefixed_bytes(full_id.as_bytes()))
    );
    assert_ne!(
        peak.read_id_sha256,
        Some(sha256_prefixed_bytes(prefix.as_bytes()))
    );

    work_profile::record_read(
        &mut profile,
        &work_profile::ReadWork {
            anchors: 2,
            template_chain_evaluations: 2,
        },
        2,
        "a-tied-later-id",
    )
    .unwrap();
    assert_eq!(profile.peak_anchor_postings.as_ref(), Some(&peak));
}

#[test]
fn rna_panel_read_work_counts_postings_before_dedup_and_rejected_chains() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(2, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    // A nine-base read has two A8 windows and its reverse has two T8 windows.
    // The duplicate posting is deliberately synthetic: 2*4 + 2*2 = 12 hits,
    // but only two seeded templates per orientation invoke LCSk++ (four total).
    let index = Index::from([
        (
            GentleEngine::encode_kmer_bits(b"AAAAAAAA").unwrap(),
            vec![(0, 0), (0, 0), (1, 0), (1, 1)],
        ),
        (
            GentleEngine::encode_kmer_bits(b"TTTTTTTT").unwrap(),
            vec![(0, 0), (2, 0)],
        ),
    ]);
    let mut policy = RnaReadPanelPolicy {
        kmer_len: 8,
        min_chain_bases: 8,
        min_junction_flank_bp: 3,
        ..Default::default()
    };
    let mut work = work_profile::ReadWork::default();
    let (class, candidates) = nominate_with_work(
        &panel,
        &index,
        &policy,
        b"AAAAAAAAA",
        &mut || false,
        Some(&mut work),
    )
    .unwrap();
    assert_eq!(class, "multi_locus_ambiguous");
    assert_eq!(candidates.len(), 4);
    assert_eq!(work.anchors, 12);
    assert_eq!(work.template_chain_evaluations, 4);

    policy.min_chain_bases = 10;
    let (class, candidates) = nominate_with_work(
        &panel,
        &index,
        &policy,
        b"AAAAAAAAA",
        &mut || false,
        Some(&mut work),
    )
    .unwrap();
    assert_eq!(class, "partial");
    assert!(candidates.is_empty());
    assert_eq!(work.anchors, 12);
    assert_eq!(work.template_chain_evaluations, 4);

    let (class, candidates) = nominate_with_work(
        &panel,
        &index,
        &policy,
        b"NNNNNNNNN",
        &mut || false,
        Some(&mut work),
    )
    .unwrap();
    assert_eq!(class, "unmatched");
    assert!(candidates.is_empty());
    assert_eq!(work.anchors, 0);
    assert_eq!(work.template_chain_evaluations, 0);
}

#[test]
fn rna_panel_optional_read_work_preserves_nominations_and_histogram_after_template_reordering() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(3, td.path());
    let mut panel = engine.build_rna_read_panel(&build).unwrap();
    let templates = panel.templates.clone();
    let reads = vec![
        templates[0].sequence.as_bytes().to_vec(),
        templates[0].sequence.as_bytes()[..40].to_vec(),
        GentleEngine::reverse_complement_bytes(templates[0].sequence.as_bytes()),
        format!("{}NNNNNNNN{}", templates[0].sequence, templates[2].sequence).into_bytes(),
        b"NNNNNNNNNNNNNNNN".to_vec(),
        templates[0].sequence.as_bytes()[..8].to_vec(),
        bases(8_751, 90).into_bytes(),
    ];
    let policy = RnaReadPanelPolicy {
        kmer_len: 8,
        min_chain_bases: 16,
        min_junction_flank_bp: 3,
        ..Default::default()
    };
    let mut expected = None;
    for order in 0..3 {
        panel.templates = templates.clone();
        match order {
            1 => panel.templates.reverse(),
            2 => panel.templates.rotate_left(1),
            _ => {}
        }
        let (index, postings) = build_index(&panel, &policy).unwrap();
        let mut profile = work_profile::file_profile();
        let mut nominations = Vec::new();
        let mut per_read_work = Vec::new();
        for (record, read) in reads.iter().enumerate() {
            let plain = nominate(&panel, &index, &policy, read, &mut || false).unwrap();
            let without_work =
                nominate_with_work(&panel, &index, &policy, read, &mut || false, None).unwrap();
            let mut work = work_profile::ReadWork::default();
            let with_work = nominate_with_work(
                &panel,
                &index,
                &policy,
                read,
                &mut || false,
                Some(&mut work),
            )
            .unwrap();
            assert_eq!(encode(&plain).unwrap(), encode(&without_work).unwrap());
            assert_eq!(encode(&plain).unwrap(), encode(&with_work).unwrap());
            work_profile::record_read(
                &mut profile,
                &work,
                record as u64,
                &format!("synthetic_{record}"),
            )
            .unwrap();
            per_read_work.push((work.anchors, work.template_chain_evaluations));
            nominations.push(plain);
        }
        assert_eq!(nominations[4].0, "unmatched");
        assert_eq!(nominations[5].0, "partial");
        assert!(profile.total_anchor_postings.is_some_and(|n| n > 0));
        assert_eq!(
            profile.anchor_histogram_counts.iter().sum::<u64>(),
            reads.len() as u64
        );
        let observed = (
            postings,
            work_profile::index_profile(&panel, &index),
            profile,
            per_read_work,
            encode(&nominations).unwrap(),
        );
        if let Some(expected) = &expected {
            assert_eq!(&observed, expected, "template order {order}");
        } else {
            expected = Some(observed);
        }
    }
}

#[test]
fn rna_panel_gzip_sidecars_and_counts_are_byte_exact_with_or_without_work_profiles() {
    for crlf in [false, true] {
        let td = tempfile::tempdir().unwrap();
        let (engine, build) = starter(2, td.path());
        let panel = engine.build_rna_read_panel(&build).unwrap();
        let shared = &panel.templates[0].sequence[..40];
        let mut run = request(&panel, td.path(), shared);
        run.policy.write_locus_fastas = true;
        let mut first_member = fastq(shared, b'>');
        first_member.extend(fastq(shared, b'>'));
        let mut second_member = fastq(&panel.templates[0].sequence, b'>');
        second_member.extend(fastq(&panel.templates[2].sequence, b'>'));
        second_member.extend(fastq(
            &format!(
                "{}NNNNNNNN{}",
                panel.templates[0].sequence, panel.templates[2].sequence
            ),
            b'>',
        ));
        second_member.extend(fastq("NNNNNNNNNNNNNNNN", b'>'));
        second_member.extend(fastq(&panel.templates[0].sequence[..8], b'>'));
        if crlf {
            first_member = String::from_utf8(first_member)
                .unwrap()
                .replace('\n', "\r\n")
                .into_bytes();
            second_member = String::from_utf8(second_member)
                .unwrap()
                .replace('\n', "\r\n")
                .into_bytes();
        }
        let mut compressed = gzip(&first_member);
        compressed.extend(gzip(&second_member));
        replace_synthetic_input(&mut run, &compressed);
        let mut off_request = run.clone();
        off_request.output_dir = td.path().join("profiling-off").display().to_string();
        run.output_dir = td.path().join("profiling-on").display().to_string();
        assert_ne!(run.output_dir, off_request.output_dir);
        assert!(!Path::new(&run.output_dir).exists());
        assert!(!Path::new(&off_request.output_dir).exists());
        assert_eq!(
            encode(&run.files).unwrap(),
            encode(&off_request.files).unwrap()
        );

        let off = screen_with_clock_and_profiling(off_request, &mut |_| true, Instant::now, false)
            .unwrap();
        let on = screen_with_clock_and_profiling(run, &mut |_| true, Instant::now, true).unwrap();
        validate_screen_report_binding(&off).unwrap();
        validate_screen_report_binding(&on).unwrap();
        assert_eq!(off.status, "completed");
        assert_eq!(on.status, off.status);
        assert_eq!(on.panel_sha256, off.panel_sha256);
        assert_eq!(on.indexed_postings, off.indexed_postings);
        assert_eq!(on.nonclaims, off.nonclaims);
        assert!(!on.isoform_assignment_performed && !off.isoform_assignment_performed);
        assert!(!on.alignment_confirmation_performed && !off.alignment_confirmation_performed);
        assert!(off.index_work_profile.is_none());
        assert!(on.index_work_profile.is_some());
        assert_eq!(on.files_completed.len(), 1);
        assert_eq!(off.files_completed.len(), 1);
        let without = &off.files_completed[0];
        let with = &on.files_completed[0];
        assert!(without.work_profile.is_none());
        assert_eq!(with.reads, 7);
        assert_eq!(with.classes, without.classes);
        assert_eq!(with.equivalence_classes, without.equivalence_classes);
        assert!(with.equivalence_classes.values().any(|count| *count > 1));
        assert_eq!(with.artifacts.len(), panel.loci.len() + 1);
        assert_eq!(
            encode(&with.artifacts).unwrap(),
            encode(&without.artifacts).unwrap()
        );
        let mut tsv_count = 0;
        let mut fasta_count = 0;
        for artifact in &with.artifacts {
            let on_bytes =
                fs::read(Path::new(&on.request.output_dir).join(&artifact.path)).unwrap();
            let off_bytes =
                fs::read(Path::new(&off.request.output_dir).join(&artifact.path)).unwrap();
            assert_eq!(on_bytes, off_bytes, "{}; crlf={crlf}", artifact.path);
            assert_eq!(artifact.sha256, sha256_prefixed_bytes(&on_bytes));
            assert_eq!(artifact.bytes, on_bytes.len() as u64);
            if artifact.path.ends_with(".tsv") {
                tsv_count += 1;
            } else {
                assert!(artifact.path.ends_with(".fasta"));
                assert!(
                    on_bytes.starts_with(b">"),
                    "both loci must have FASTA support"
                );
                fasta_count += 1;
            }
        }
        assert_eq!(tsv_count, 1);
        assert_eq!(fasta_count, 2);
        let profile = with.work_profile.as_ref().unwrap();
        assert_eq!(
            profile.anchor_histogram_upper_bounds,
            EXPECTED_ANCHOR_BOUNDS
        );
        assert_eq!(
            profile.anchor_histogram_counts.iter().sum::<u64>(),
            with.reads
        );
        assert_eq!(
            profile.equivalence_key_bytes_added,
            with.equivalence_classes
                .keys()
                .map(String::len)
                .sum::<usize>()
        );
        assert!(profile.total_anchor_postings.is_some_and(|n| n > 0));
        assert!(profile.peak_template_chain_evaluations >= 2);
        let mut without_profile = with.clone();
        without_profile.work_profile = None;
        assert_eq!(encode(&without_profile).unwrap(), encode(without).unwrap());
    }
}

#[test]
fn rna_panel_verified_empty_gzip_has_an_explicit_zero_work_profile_by_default() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = starter(1, td.path());
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let mut run = request(&panel, td.path(), "NNNNNNNN");
    let compressed = gzip(b"");
    replace_synthetic_input(&mut run, &compressed);
    let report = screen(run, &mut |_| true).unwrap();
    validate_screen_report_binding(&report).unwrap();
    assert_eq!(report.status, "completed");
    assert!(report.index_work_profile.is_some());
    assert_eq!(report.files_completed.len(), 1);
    let file = &report.files_completed[0];
    assert_eq!(file.reads, 0);
    assert!(file.classes.is_empty());
    assert!(file.equivalence_classes.is_empty());
    assert_eq!(file.observed_sha256, sha256_prefixed_bytes(&compressed));
    assert_eq!(file.compressed_bytes, compressed.len() as u64);
    assert_eq!(
        file.work_profile,
        Some(RnaReadPanelFileWorkProfile {
            peak_anchor_postings: None,
            anchor_histogram_upper_bounds: EXPECTED_ANCHOR_BOUNDS.to_vec(),
            anchor_histogram_counts: vec![0; 20],
            total_anchor_postings: Some(0),
            peak_template_chain_evaluations: 0,
            equivalence_key_bytes_added: 0,
        })
    );
}

#[test]
fn rna_panel_legacy_v2_without_work_fields_roundtrips_bytes_and_content_digest() {
    // Independent hand-authored wire fixture: deliberately do not serialize a
    // current report and strip fields, which could hide absent-field regressions.
    let input = json!({
        "sample_id": "synthetic_sample",
        "barcode": "synthetic_barcode",
        "priming": "hand_crafted",
        "path": "synthetic.fastq.gz",
        "declared_complete": true,
        "expected_md5": null,
        "expected_sha256": "sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
        "checksum_manifest_key": null
    });
    let mut legacy = json!({
        "schema": "gentle.rna_read_panel_screen.v2",
        "content_sha256": "",
        "status": "completed",
        "request": {
            "schema": "gentle.rna_read_panel_screen_request.v1",
            "panel_path": "synthetic-panel.json",
            "files": [input.clone()],
            "md5_manifest_path": null,
            "output_dir": "synthetic-output",
            "policy": {
                "kmer_len": 8,
                "min_chain_bases": 16,
                "min_junction_flank_bp": 3,
                "max_index_postings": 2000000,
                "max_anchors_per_read": 100000,
                "max_equivalence_classes": 4096,
                "write_locus_fastas": false
            },
            "cancel_path": null
        },
        "panel_sha256": "sha256:bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
        "indexed_postings": 1,
        "isoform_assignment_performed": false,
        "alignment_confirmation_performed": false,
        "files_completed": [{
            "file_id": "file_0000",
            "input": input,
            "observed_md5": "cccccccccccccccccccccccccccccccc",
            "observed_sha256": "sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            "compressed_bytes": 100,
            "reads": 1,
            "classes": {"unmatched": 1},
            "equivalence_classes": {"{\"class\":\"unmatched\",\"member_ids\":[]}": 1},
            "artifacts": [{
                "path": "artifacts/file_0000/nominations.tsv",
                "sha256": "sha256:dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd",
                "bytes": 200
            }]
        }],
        "files_not_completed": [],
        "nonclaims": ["Synthetic wire compatibility only; no biological evidence."]
    });
    let digest = sha256_prefixed_bytes(&encode(&legacy).unwrap());
    legacy["content_sha256"] = json!(digest);
    let bytes = encode(&legacy).unwrap();
    let restored: RnaReadPanelScreenReport = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(restored.schema, RNA_PANEL_SCREEN_SCHEMA);
    assert!(restored.index_work_profile.is_none());
    assert!(restored.files_completed[0].work_profile.is_none());
    assert_eq!(restored.content_sha256, digest);
    validate_screen_report_binding(&restored).unwrap();
    assert_eq!(encode(&restored).unwrap(), bytes);
    let roundtrip = serde_json::to_value(&restored).unwrap();
    assert!(roundtrip.get("index_work_profile").is_none());
    assert!(
        roundtrip["files_completed"][0]
            .get("work_profile")
            .is_none()
    );
    let mut basis = restored.clone();
    basis.content_sha256.clear();
    assert_eq!(sha256_prefixed_bytes(&encode(&basis).unwrap()), digest);
}
