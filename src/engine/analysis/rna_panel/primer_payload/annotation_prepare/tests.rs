//! Inline hand-crafted GTF/dictionary fixtures; recreated here in isolated tempdirs.
//! Used only by S2.1b file-boundary tests; no real annotation, genome or read data.

use super::*;
use serde_json::json;

#[path = "adapter_tests.rs"]
mod adapters;

impl GentleEngine {
    /// Inline synthetic annotation inputs, absent from non-test builds.
    pub(crate) fn observed_annotation_test_fixture(
        root: &Path,
    ) -> PrepareObservedTranscriptAnnotationRequest {
        fixture(root)
    }
}

const GTF: &str = concat!(
    "# synthetic snapshot\n",
    "chrToy\ttest\ttranscript\t2\t14\t.\t-\t.\tgene_id \"gene_minus\"; transcript_id \"tx.7\"; transcript_version \"7\";\n",
    "chrToy\ttest\texon\t2\t5\t.\t-\t.\tgene_id \"gene_minus\"; transcript_id \"tx.7\"; transcript_version \"7\";\n",
    "chrToy\ttest\texon\t10\t14\t.\t-\t.\tgene_id \"gene_minus\"; transcript_id \"tx.7\"; transcript_version \"7\";\n",
);

fn bound(path: &Path) -> PrimerPayloadFile {
    let raw = fs::read(path).unwrap();
    PrimerPayloadFile {
        path: path.display().to_string(),
        bytes: raw.len() as u64,
        sha256: sha256_prefixed_bytes(&raw),
    }
}

fn fixture(root: &Path) -> PrepareObservedTranscriptAnnotationRequest {
    let gtf = root.join("a gene's source.gtf");
    fs::write(&gtf, GTF).unwrap();
    let dictionary = root.join("declared reference.json");
    fs::write(&dictionary, serde_json::to_vec(&json!({
        "schema": REFERENCE_CONTIG_DICTIONARY_SCHEMA, "reference_id":"toy_ref",
        "assembly":"toy_assembly", "reference_fasta_sha256":format!("sha256:{}", "4".repeat(64)),
        "contigs":[{"name":"chrToy","length":1000},{"name":"competitor","length":800}]
    })).unwrap()).unwrap();
    PrepareObservedTranscriptAnnotationRequest {
        schema: OBSERVED_ANNOTATION_REQUEST_SCHEMA.into(),
        format: ObservedAnnotationFormat::Gtf,
        gtf: bound(&gtf),
        dictionary: bound(&dictionary),
        reference_id: "toy_ref".into(),
        reference: crate::observed_transcript_structures::ReferenceBinding {
            assembly: "toy_assembly".into(),
            fasta_sha256: format!("sha256:{}", "4".repeat(64)),
        },
        selector: crate::observed_transcript_structures::ReferenceGene {
            source_id: "synthetic_source".into(),
            gene_id: "gene_minus".into(),
        },
        source: ObservedAnnotationSourceDeclaration {
            scope: ObservedAnnotationInputScope::ProducerExcerpt,
            provider: Some("synthetic provider".into()),
            declared_provider_release: Some("DECLARED-not-verified".into()),
            declared_derivation: Some("inline hand-crafted fixture".into()),
            source_uri: None,
            notice: None,
        },
        policy: ObservedAnnotationPolicy {
            transcript_feature_rows: ObservedAnnotationTranscriptRows::Required,
            minimum_reference_gap_bases: 1,
        },
        budgets: ObservedAnnotationBudgets {
            max_input_bytes: 65536,
            max_dictionary_bytes: 8192,
            max_line_bytes: 2048,
            max_lines: 100,
            max_dictionary_entries: 10,
            max_indexed_genes: 10,
            max_indexed_transcripts: 10,
            max_exon_rows: 100,
            max_retained_bytes: 8 * 1024 * 1024,
            max_work: 100000,
            max_diagnostics: 10,
            max_output_bytes: 1024 * 1024,
        },
        output_dir: root.join("prepared gene's evidence").display().to_string(),
        cancel_path: None,
    }
}

fn run(r: PrepareObservedTranscriptAnnotationRequest) -> ObservedTranscriptAnnotationPreparation {
    GentleEngine::new()
        .prepare_observed_transcript_annotation(r, &mut |_| true)
        .unwrap()
}

fn no_receipt(r: PrepareObservedTranscriptAnnotationRequest) -> EngineError {
    let out = PathBuf::from(&r.output_dir);
    let e = GentleEngine::new()
        .prepare_observed_transcript_annotation(r, &mut |_| true)
        .unwrap_err();
    assert!(!out.exists(), "failed input left public output: {e:?}");
    e
}

#[test]
fn observed_annotation_file_success_is_bound_read_only_and_report_last() {
    let root = tempfile::tempdir().unwrap();
    let r = fixture(root.path());
    let engine = GentleEngine::new();
    let before = serde_json::to_value(engine.state()).unwrap();
    let a = engine
        .prepare_observed_transcript_annotation(r.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(before, serde_json::to_value(engine.state()).unwrap());
    let report = &a.report;
    assert_eq!(report.status, ObservedAnnotationExecutionStatus::Completed);
    assert!(report.dictionary_verified && report.gtf_verified);
    let analysis = report.analysis.as_ref().unwrap();
    assert_eq!(analysis.outcome, ObservedAnnotationOutcome::Completed);
    assert_eq!(analysis.dictionary_entries, 2);
    assert_eq!(analysis.counts.lines, 4);
    assert_eq!(analysis.counts.exon_rows, 2);
    assert_eq!(analysis.coordinate_authority, "declared");
    assert_eq!(analysis.raw_sha256, r.gtf.sha256);
    let artifact = report.annotation.as_ref().unwrap();
    let path = Path::new(&r.output_dir).join(&artifact.path);
    assert_eq!(bound(&path).sha256, artifact.sha256);
    let annotation: PrimerPayloadAnnotation =
        serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    assert_eq!(annotation.transcripts[0].genomic_strand, Some('-'));
    assert_eq!(
        annotation.transcripts[0].genomic_exons,
        vec![(1, 5), (9, 14)]
    );
    assert_eq!(annotation.release, format!("snapshot-{}", r.gtf.sha256));
    let file = bound(Path::new(&a.report_file.path));
    assert_eq!(file.sha256, a.report_file.sha256);
    let saved: ObservedTranscriptAnnotationReport =
        serde_json::from_slice(&fs::read(&file.path).unwrap()).unwrap();
    assert_eq!(saved.content_sha256, report.content_sha256);
    let mut basis = saved;
    basis.content_sha256.clear();
    assert_eq!(
        semantic_hash(&basis, r.budgets.max_retained_bytes).unwrap(),
        report.content_sha256
    );
    assert_ne!(file.sha256, report.content_sha256);
    assert_eq!(
        report.request_sha256,
        semantic_hash(&r, r.budgets.max_retained_bytes).unwrap()
    );
    // A second publication cannot adopt or clobber the first.
    assert!(
        engine
            .prepare_observed_transcript_annotation(r, &mut |_| true)
            .is_err()
    );
    assert_eq!(bound(Path::new(&file.path)).sha256, file.sha256);
}

#[test]
fn observed_annotation_lf_crlf_final_line_and_policies_bind_distinct_receipts() {
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    let first = run(r.clone());
    let mut raw = GTF.replace('\n', "\r\n");
    raw.truncate(raw.len() - 2);
    fs::write(&r.gtf.path, raw).unwrap();
    r.gtf = bound(Path::new(&r.gtf.path));
    r.output_dir = root.path().join("CRLF output").display().to_string();
    let second = run(r.clone());
    let annotation = |p: &ObservedTranscriptAnnotationPreparation| {
        let path = Path::new(&p.report.request.output_dir)
            .join(&p.report.annotation.as_ref().unwrap().path);
        serde_json::from_slice::<PrimerPayloadAnnotation>(&fs::read(path).unwrap()).unwrap()
    };
    assert_eq!(
        serde_json::to_value(annotation(&first).transcripts).unwrap(),
        serde_json::to_value(annotation(&second).transcripts).unwrap()
    );
    assert_ne!(annotation(&first).release, annotation(&second).release);
    r.policy.minimum_reference_gap_bases = 3;
    r.output_dir = root
        .path()
        .join("other admitted policy")
        .display()
        .to_string();
    let third = run(r);
    assert_eq!(
        second.report.annotation.as_ref().unwrap().sha256,
        third.report.annotation.as_ref().unwrap().sha256
    );
    assert_ne!(second.report.content_sha256, third.report.content_sha256);
    // The policy itself affects the semantic request, independently of destination.
    let mut same_path = third.report.request.clone();
    same_path.output_dir = second.report.request.output_dir.clone();
    assert_ne!(
        semantic_hash(&same_path, same_path.budgets.max_retained_bytes).unwrap(),
        second.report.request_sha256
    );
}

#[test]
fn observed_annotation_missing_and_refused_genes_publish_receipt_only() {
    let root = tempfile::tempdir().unwrap();
    let original = fixture(root.path());
    for (i, outcome) in [
        ObservedAnnotationOutcome::MissingGene,
        ObservedAnnotationOutcome::ReferenceGapPolicyRefused,
        ObservedAnnotationOutcome::IdenticalDuplicateRow,
        ObservedAnnotationOutcome::TranscriptModelConflict,
    ]
    .into_iter()
    .enumerate()
    {
        let mut r = original.clone();
        fs::write(&r.gtf.path, GTF).unwrap();
        match i {
            0 => r.selector.gene_id = "GENE_minus".into(),
            1 => r.policy.minimum_reference_gap_bases = 5,
            2 => fs::write(&r.gtf.path, format!("{GTF}{}", GTF.lines().last().unwrap())).unwrap(),
            _ => fs::write(&r.gtf.path, GTF.replace("\ttranscript\t", "\tCDS\t")).unwrap(),
        }
        r.gtf = bound(Path::new(&r.gtf.path));
        r.output_dir = root
            .path()
            .join(format!("refusal {i}"))
            .display()
            .to_string();
        let p = run(r.clone());
        assert_eq!(p.report.analysis.as_ref().unwrap().outcome, outcome);
        assert!(p.report.annotation.is_none());
        assert_eq!(
            fs::read_dir(Path::new(&r.output_dir).join("artifacts"))
                .unwrap()
                .count(),
            0
        );
    }
}

#[test]
fn observed_annotation_dictionary_is_strict_full_scope_and_metadata_bound() {
    let root = tempfile::tempdir().unwrap();
    let original = fixture(root.path());
    let raw = fs::read_to_string(&original.dictionary.path).unwrap();
    for bad_json in [
        raw.replacen("{", "{\"schema\":\"duplicate\",", 1),
        raw.replacen("{", "{\"unknown\":1,", 1),
        raw.replace("competitor", "chrToy"),
        raw.replace("toy_assembly", "other_assembly"),
        raw.replace("1000", "0"),
        raw.replace("competitor", "bad contig"),
        raw.replace("1000", "18446744073709551616"),
    ] {
        let mut r = original.clone();
        fs::write(&r.dictionary.path, bad_json).unwrap();
        r.dictionary = bound(Path::new(&r.dictionary.path));
        no_receipt(r);
    }
    fs::write(&original.dictionary.path, raw).unwrap();
    let mut r = original;
    fs::write(&r.gtf.path, format!("{GTF}unknown\ttest\texon\t1\t4\t.\t+\t.\tgene_id \"other\"; transcript_id \"other_t\";\n")).unwrap();
    r.gtf = bound(Path::new(&r.gtf.path));
    assert!(no_receipt(r).message.contains("contig"));
}

#[test]
fn observed_annotation_hash_parse_and_budget_failures_never_publish() {
    let root = tempfile::tempdir().unwrap();
    let original = fixture(root.path());
    let mut r = original.clone();
    r.gtf.sha256 = format!("sha256:{}", "0".repeat(64));
    assert!(no_receipt(r).message.contains("hash"));
    for i in 0..10 {
        let mut r = original.clone();
        match i {
            0 => r.budgets.max_input_bytes = r.gtf.bytes + r.dictionary.bytes - 1,
            1 => r.budgets.max_dictionary_bytes = r.dictionary.bytes as usize - 1,
            2 => r.budgets.max_line_bytes = 20,
            3 => r.budgets.max_lines = 3,
            4 => r.budgets.max_dictionary_entries = 1,
            5 => r.budgets.max_exon_rows = 1,
            6 => r.budgets.max_retained_bytes = 1000,
            7 => r.budgets.max_work = 1,
            8 => r.budgets.max_output_bytes = 10,
            _ => r.schema = "unknown".into(),
        }
        no_receipt(r);
    }
    fs::write(&original.gtf.path, format!("{GTF}malformed\n")).unwrap();
    let mut r = original;
    r.gtf = bound(Path::new(&r.gtf.path));
    no_receipt(r);
}

#[test]
fn observed_annotation_early_and_late_control_have_distinct_analysis_states() {
    let root = tempfile::tempdir().unwrap();
    let original = fixture(root.path());
    for late in [false, true] {
        let mut r = original.clone();
        r.output_dir = root
            .path()
            .join(format!("cancel {late}"))
            .display()
            .to_string();
        let p = GentleEngine::new().prepare_observed_transcript_annotation(r, &mut |progress| {
            !matches!(progress, OperationProgress::Workflow { completed, .. } if !late || completed == 3)
        }).unwrap();
        assert_eq!(
            p.report.status,
            ObservedAnnotationExecutionStatus::Cancelled
        );
        assert_eq!(p.report.analysis.is_some(), late);
        assert_eq!(p.report.annotation.is_some(), late);
        assert_eq!(p.report.gtf_verified, late);
    }
    for late in [false, true] {
        let mut r = original.clone();
        let marker = root.path().join(format!("bad marker {late}"));
        r.cancel_path = Some(marker.display().to_string());
        r.output_dir = root
            .path()
            .join(format!("control failure {late}"))
            .display()
            .to_string();
        if !late {
            fs::write(&marker, b"wrong identity").unwrap();
        }
        let p = GentleEngine::new()
            .prepare_observed_transcript_annotation(r, &mut |progress| {
                if matches!(progress, OperationProgress::Workflow { completed: 3, .. }) {
                    fs::write(&marker, b"wrong identity").unwrap();
                }
                true
            })
            .unwrap();
        assert_eq!(
            p.report.status,
            ObservedAnnotationExecutionStatus::ControlError
        );
        assert!(
            p.report
                .control_error
                .as_ref()
                .unwrap()
                .contains("mismatch")
        );
        assert_eq!(p.report.analysis.is_some(), late);
    }
}

#[test]
fn observed_annotation_marker_is_request_bound_and_cancels_before_reading() {
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    let marker = root.path().join("request cancel marker");
    r.cancel_path = Some(marker.display().to_string());
    fs::write(
        &marker,
        semantic_hash(&r, r.budgets.max_retained_bytes).unwrap(),
    )
    .unwrap();
    let p = run(r);
    assert_eq!(
        p.report.status,
        ObservedAnnotationExecutionStatus::Cancelled
    );
    assert!(p.report.analysis.is_none() && p.report.annotation.is_none());
    assert!(!p.report.dictionary_verified && !p.report.gtf_verified);
}

#[test]
fn observed_annotation_raw_paths_and_hardlink_aliases_are_refused() {
    let root = tempfile::tempdir().unwrap();
    let original = fixture(root.path());
    for suffix in ["/./result", "/../result", "//result", "/result/"] {
        let mut r = original.clone();
        r.output_dir = format!("{}{suffix}", root.path().display());
        assert!(publication::Paths::new(&r).is_err());
    }
    let mut r = original.clone();
    r.cancel_path = Some(r.gtf.path.clone());
    no_receipt(r);
    let mut r = original.clone();
    r.cancel_path = Some(format!("{}/marker", r.output_dir));
    no_receipt(r);
    let alias = root.path().join("hardlink marker");
    fs::hard_link(&original.gtf.path, &alias).unwrap();
    let mut r = original.clone();
    r.cancel_path = Some(alias.display().to_string());
    no_receipt(r);
    let mut r = original;
    r.dictionary.path = alias.display().to_string();
    no_receipt(r);
}

#[cfg(unix)]
#[test]
fn observed_annotation_symlinks_and_replaced_input_identity_are_refused() {
    use std::os::unix::fs::symlink;
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    let link = root.path().join("linked input.gtf");
    symlink(&r.gtf.path, &link).unwrap();
    let mut linked = r.clone();
    linked.gtf.path = link.display().to_string();
    no_receipt(linked);
    let paths = publication::Paths::new(&r).unwrap();
    let replacement = root.path().join("replacement.gtf");
    fs::write(&replacement, GTF).unwrap();
    fs::rename(&replacement, &r.gtf.path).unwrap();
    assert!(paths.reader(&r.gtf, r.budgets.max_input_bytes).is_err());
    // A symlinked ancestor resolving to an input cannot conceal a marker alias.
    let parent_link = root.path().join("parent alias");
    symlink(root.path(), &parent_link).unwrap();
    r.cancel_path = Some(
        parent_link
            .join(Path::new(&r.gtf.path).file_name().unwrap())
            .display()
            .to_string(),
    );
    no_receipt(r);
}

#[test]
fn observed_annotation_publication_visibility_and_fault_cleanup_are_owned() {
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    let p = run(r.clone());
    let artifact_raw =
        fs::read(Path::new(&r.output_dir).join("artifacts/annotation.json")).unwrap();
    for fail in [None, Some(0), Some(1), Some(2), Some(3)] {
        r.output_dir = root
            .path()
            .join(format!("fault {fail:?}"))
            .display()
            .to_string();
        let paths = publication::Paths::new(&r).unwrap();
        let stage = staging(&r.output_dir).unwrap();
        fs::write(stage.path().join("annotation.json"), &artifact_raw).unwrap();
        let out = PathBuf::from(&r.output_dir);
        let mut report = p.report.clone();
        report.request = r.clone();
        let mut observations = Vec::new();
        let mut index = 0;
        let result = publication::publish(stage, &paths, &report, &mut |phase| {
            observations.push(format!("{phase:?}"));
            assert!(!out.join("report.json").exists());
            if matches!(
                phase,
                publication::Phase::ArtifactsVisible | publication::Phase::ReportReady
            ) {
                assert!(out.join("artifacts/annotation.json").exists());
            }
            let fail_here = fail == Some(index);
            index += 1;
            if fail_here {
                Err(bad("synthetic publication fault"))
            } else {
                Ok(())
            }
        });
        assert_eq!(observations.len(), fail.map_or(4, |n| n + 1));
        if fail.is_some() {
            assert!(result.is_err());
            assert!(!out.exists());
        } else {
            assert!(result.is_ok());
            assert!(out.join("report.json").exists());
        }
    }
}

#[test]
fn observed_annotation_publication_preserves_replaced_artifacts() {
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    let p = run(r.clone());
    let raw = fs::read(Path::new(&r.output_dir).join("artifacts/annotation.json")).unwrap();
    for with_annotation in [false, true] {
        r.output_dir = root
            .path()
            .join(format!("replacement {with_annotation}"))
            .display()
            .to_string();
        let mut report = p.report.clone();
        report.request = r.clone();
        if !with_annotation {
            report.annotation = None;
        }
        let paths = publication::Paths::new(&r).unwrap();
        let stage = staging(&r.output_dir).unwrap();
        if with_annotation {
            fs::write(stage.path().join("annotation.json"), &raw).unwrap();
        }
        let out = PathBuf::from(&r.output_dir);
        let moved = root
            .path()
            .join(format!("original artifacts {with_annotation}"));
        let e = publication::publish(stage, &paths, &report, &mut |phase| {
            if matches!(phase, publication::Phase::ArtifactsVisible) {
                fs::rename(out.join("artifacts"), &moved).unwrap();
                fs::create_dir(out.join("artifacts")).unwrap();
                fs::write(
                    out.join("artifacts/annotation.json"),
                    b"unowned replacement",
                )
                .unwrap();
                return Err(bad("synthetic replacement fault"));
            }
            Ok(())
        })
        .unwrap_err();
        assert!(e.message.contains("observed_annotation_orphaned_output"));
        assert_eq!(
            fs::read(out.join("artifacts/annotation.json")).unwrap(),
            b"unowned replacement"
        );
        assert!(!out.join("report.json").exists());
        if with_annotation {
            assert_eq!(fs::read(moved.join("annotation.json")).unwrap(), raw);
        } else {
            assert_eq!(fs::read_dir(moved).unwrap().count(), 0);
        }
    }
}

#[test]
fn observed_annotation_publication_failed_reacquire_discloses_orphan() {
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    let mut p = run(r.clone());
    r.output_dir = root.path().join("reacquire failure").display().to_string();
    p.report.request = r.clone();
    let paths = publication::Paths::new(&r).unwrap();
    let stage = staging(&r.output_dir).unwrap();
    let original = stage.path().join("annotation.json");
    fs::write(&original, b"synthetic annotation").unwrap();
    let e = publication::publish(stage, &paths, &p.report, &mut |phase| {
        if matches!(phase, publication::Phase::Reserved) {
            fs::remove_file(&original).unwrap();
        }
        Ok(())
    })
    .unwrap_err();
    assert!(e.message.contains("observed_annotation_orphaned_output"));
    assert!(Path::new(&r.output_dir).join("artifacts").is_dir());
    assert!(!Path::new(&r.output_dir).join("report.json").exists());
}

#[test]
fn observed_annotation_publication_discloses_unowned_orphan_without_deleting_it() {
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    let mut p = run(r.clone());
    r.output_dir = root.path().join("orphan guard").display().to_string();
    p.report.request = r.clone();
    p.report.annotation = None;
    let paths = publication::Paths::new(&r).unwrap();
    let stage = staging(&r.output_dir).unwrap();
    let out = PathBuf::from(&r.output_dir);
    let e = publication::publish(stage, &paths, &p.report, &mut |phase| {
        if matches!(phase, publication::Phase::ArtifactsVisible) {
            fs::write(out.join("unowned evidence"), b"preserve").unwrap();
            return Err(bad("synthetic fault"));
        }
        Ok(())
    })
    .unwrap_err();
    assert!(e.message.contains("orphaned_output"));
    assert_eq!(fs::read(out.join("unowned evidence")).unwrap(), b"preserve");
    assert!(!out.join("report.json").exists());
}

#[test]
fn observed_annotation_new_wire_contract_rejects_unknown_and_duplicate_fields() {
    let root = tempfile::tempdir().unwrap();
    let r = fixture(root.path());
    let raw = serde_json::to_string(&r).unwrap();
    let bad = raw.replacen("{", "{\"schema\":\"duplicate\",", 1);
    assert!(serde_json::from_str::<PrepareObservedTranscriptAnnotationRequest>(&bad).is_err());
    let mut value = serde_json::to_value(&r).unwrap();
    value["policy"]["guess_gene"] = json!(true);
    assert!(serde_json::from_value::<PrepareObservedTranscriptAnnotationRequest>(value).is_err());
    let value = serde_json::to_value(run(r)).unwrap();
    let roundtrip: ObservedTranscriptAnnotationPreparation = serde_json::from_value(value).unwrap();
    assert!(roundtrip.report.analysis.is_some());
}

#[test]
fn observed_annotation_verified_input_without_finished_analysis_has_no_counts() {
    let root = tempfile::tempdir().unwrap();
    let original = fixture(root.path());
    for stop_at in [1, 2] {
        let mut r = original.clone();
        r.output_dir = root
            .path()
            .join(format!("partial {stop_at}"))
            .display()
            .to_string();
        let p = GentleEngine::new().prepare_observed_transcript_annotation(r, &mut |progress| {
            !matches!(progress, OperationProgress::Workflow { completed, .. } if completed == stop_at)
        }).unwrap();
        assert_eq!(
            p.report.status,
            ObservedAnnotationExecutionStatus::Cancelled
        );
        assert!(p.report.dictionary_verified);
        assert_eq!(p.report.gtf_verified, stop_at == 2);
        assert!(p.report.analysis.is_none() && p.report.annotation.is_none());
        let value = serde_json::to_value(p.report).unwrap();
        assert_eq!(value["analysis"], serde_json::Value::Null);
    }
    let mut r = original;
    let marker = root.path().join("late invalid directory marker");
    r.cancel_path = Some(marker.display().to_string());
    r.output_dir = root
        .path()
        .join("late directory failure")
        .display()
        .to_string();
    let p = GentleEngine::new()
        .prepare_observed_transcript_annotation(r, &mut |progress| {
            if matches!(progress, OperationProgress::Workflow { completed: 3, .. }) {
                fs::create_dir(&marker).unwrap();
            }
            true
        })
        .unwrap();
    assert_eq!(
        p.report.status,
        ObservedAnnotationExecutionStatus::ControlError
    );
    assert!(p.report.analysis.is_some() && p.report.annotation.is_some());
}

#[test]
fn observed_annotation_changed_input_during_work_never_publishes() {
    let root = tempfile::tempdir().unwrap();
    let r = fixture(root.path());
    let source = r.gtf.path.clone();
    let out = PathBuf::from(&r.output_dir);
    let e = GentleEngine::new()
        .prepare_observed_transcript_annotation(r, &mut |progress| {
            if matches!(progress, OperationProgress::Workflow { completed: 1, .. }) {
                fs::write(&source, b"changed input\n").unwrap();
            }
            true
        })
        .unwrap_err();
    assert!(e.message.contains("file_size"));
    assert!(!out.exists());
}

#[test]
fn observed_annotation_publication_cannot_clobber_racing_report_or_old_stage_name() {
    let root = tempfile::tempdir().unwrap();
    let mut r = fixture(root.path());
    let mut p = run(r.clone());
    for collide in [false, true] {
        r.output_dir = root
            .path()
            .join(format!("racing report {collide}"))
            .display()
            .to_string();
        p.report.request = r.clone();
        p.report.annotation = None;
        let paths = publication::Paths::new(&r).unwrap();
        let stage = staging(&r.output_dir).unwrap();
        let old_stage = stage.path().to_owned();
        let out = PathBuf::from(&r.output_dir);
        let result = publication::publish(stage, &paths, &p.report, &mut |phase| {
            if matches!(phase, publication::Phase::ArtifactsVisible) {
                fs::create_dir(&old_stage).unwrap();
                fs::write(old_stage.join("unowned"), b"preserve old name").unwrap();
                if collide {
                    fs::write(out.join("report.json"), b"unowned report").unwrap();
                }
            }
            Ok(())
        });
        assert_eq!(
            fs::read(old_stage.join("unowned")).unwrap(),
            b"preserve old name"
        );
        if collide {
            assert!(result.unwrap_err().message.contains("orphaned_output"));
            assert_eq!(
                fs::read(out.join("report.json")).unwrap(),
                b"unowned report"
            );
        } else {
            assert!(result.is_ok());
        }
    }
}
