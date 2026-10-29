//! Inline synthetic S1/Phase 2 fixtures; no downloaded RNA/reference or committed artifacts.

use super::*;

fn fixture(root: &Path) -> (GentleEngine, AssessObservedSpliceJunctionsRequest) {
    let (engine, req) = GentleEngine::observed_transcripts_test_fixture(root);
    let report = engine
        .summarize_observed_transcript_structures(req, &mut |_| true)
        .unwrap();
    let binding =
        GentleEngine::observed_transcript_structures_result_binding(&report.request, &report)
            .unwrap();
    let file = super::super::tests::bound_file(Path::new(&binding.report_path));
    (
        engine,
        AssessObservedSpliceJunctionsRequest {
            schema: OBSERVED_JUNCTION_REQUEST_SCHEMA.into(),
            source: ObservedJunctionSource::S1Report {
                result: binding,
                report: file,
            },
            reference_sequence: None,
            sample_manifest: None,
            near_window_bp: 5,
            recurrence_policy: None,
            budgets: ObservedJunctionBudgets {
                max_input_bytes: 128_000_000,
                max_reference_bytes: 128_000_000,
                max_manifest_bytes: 1_000_000,
                max_distinct_junctions: 1000,
                max_annotated_junctions_compared: 10000,
                max_near_hits_per_junction: 1000,
                max_total_work: 1_000_000_000,
                max_output_bytes: 128_000_000,
            },
            output_dir: root
                .join("junction assessment's result")
                .display()
                .to_string(),
            cancel_path: None,
        },
    )
}

fn data(r: &ObservedSpliceJunctionAssessmentReport) -> ObservedJunctionAssessmentData {
    serde_json::from_slice(
        &fs::read(Path::new(&r.request.output_dir).join("artifacts/assessment.json")).unwrap(),
    )
    .unwrap()
}

#[test]
fn junction_assessment_saved_source_is_distinct_read_only_and_inspectable() {
    let root = tempfile::tempdir().unwrap();
    let (engine, req) = fixture(root.path());
    let before = encode(engine.state()).unwrap();
    let journal = engine.journal_len();
    let r = engine
        .assess_observed_splice_junctions(req, &mut |_| true)
        .unwrap();
    assert_eq!(r.status, "completed");
    assert!(r.counts.is_some());
    let d = data(&r);
    assert_eq!(d.junctions.len(), 1);
    assert!(d.junctions[0].sources.iter().any(|s| s.relation == "exact"));
    assert_eq!(
        d.junctions[0].motif.motif_class,
        ObservedMotifClass::NotAssessable
    );
    assert!(d.junctions[0].recurrence.is_none());
    assert!(d.chains.iter().all(|c| c.rna_end_status == "not_verified"));
    assert_eq!(encode(engine.state()).unwrap(), before);
    assert_eq!(engine.journal_len(), journal);
    let binding = GentleEngine::observed_splice_junction_result_binding(&r.request, &r).unwrap();
    let inspected = GentleEngine::inspect_observed_splice_junctions(
        &binding.report_path,
        Some(&binding),
        "junction_assessment.tsv",
        1,
        true,
    )
    .unwrap();
    assert!(inspected.artifacts_verified);
    assert_eq!(inspected.rows.len(), 1);
}

#[test]
fn junction_assessment_wire_rejects_defaults_unknown_fields_and_missing_budgets() {
    let root = tempfile::tempdir().unwrap();
    let (_, req) = fixture(root.path());
    let mut v = serde_json::to_value(&req).unwrap();
    v["near_window_bp"] = json!(0);
    let r: AssessObservedSpliceJunctionsRequest = serde_json::from_value(v).unwrap();
    assert!(
        GentleEngine::new()
            .assess_observed_splice_junctions(r, &mut |_| true)
            .is_err()
    );
    let mut v = serde_json::to_value(&req).unwrap();
    v["invented_default_gene"] = json!("bad");
    assert!(serde_json::from_value::<AssessObservedSpliceJunctionsRequest>(v).is_err());
    let mut v = serde_json::to_value(req).unwrap();
    v["budgets"]
        .as_object_mut()
        .unwrap()
        .remove("max_total_work");
    assert!(serde_json::from_value::<AssessObservedSpliceJunctionsRequest>(v).is_err());
}

fn model(exons: &[(usize, usize)], strand: char) -> PrimerPayloadTranscript {
    PrimerPayloadTranscript {
        transcript_id: "synthetic_model".into(),
        gene_ids: vec!["synthetic_gene".into()],
        target_name: "synthetic_contig".into(),
        target_length: 1000,
        genomic_contig: Some("synthetic_contig".into()),
        genomic_strand: Some(strand),
        genomic_exons: exons.to_vec(),
    }
}

#[test]
fn junction_assessment_exact_and_every_one_to_five_base_shift_are_not_snapped() {
    let m = model(&[(10, 100), (200, 260)], '+');
    for a in -5i64..=5 {
        for z in -5i64..=5 {
            let skip = ((100 + a) as usize, (200 + z) as usize);
            let hits = near_hits_for_models(skip, std::slice::from_ref(&m), 5);
            assert_eq!(hits.len(), 1);
            assert_eq!(hits[0].start_offset, -a);
            assert_eq!(hits[0].end_offset, -z);
            let expected = if a == 0 && z == 0 {
                "exact"
            } else if a == 0 {
                "shared_start_only"
            } else if z == 0 {
                "shared_end_only"
            } else {
                "near"
            };
            assert_eq!(near_relation(&hits), expected);
        }
    }
    assert!(near_hits_for_models((106, 206), &[m], 5).is_empty());
}

#[test]
fn junction_assessment_motif_orientation_softmask_and_ambiguity_are_explicit() {
    for (first, last, expected) in [
        ("gt", "ag", ObservedMotifClass::GtAg),
        ("GC", "AG", ObservedMotifClass::GcAg),
        ("AT", "AC", ObservedMotifClass::AtAc),
        ("AA", "AA", ObservedMotifClass::Other),
    ] {
        let plus = motif_from_dinucleotides(Some(first), Some(last), Some('+'), "verified");
        assert_eq!(plus.motif_class, expected);
        assert_eq!(plus.reference_forward_first.as_deref(), Some(first));
    }
    let minus = motif_from_dinucleotides(Some("ct"), Some("ac"), Some('-'), "declared");
    assert_eq!(minus.motif_class, ObservedMotifClass::GtAg);
    assert_eq!(
        minus.opposite_orientation_motif_class,
        ObservedMotifClass::Other
    );
    for (first, last, strand) in [
        (Some("NN"), Some("AG"), Some('+')),
        (None, Some("AG"), Some('+')),
        (Some("GT"), Some("AG"), None),
    ] {
        assert_eq!(
            motif_from_dinucleotides(first, last, strand, "verified").motif_class,
            ObservedMotifClass::NotAssessable
        );
    }
}

#[test]
fn junction_assessment_geometry_labels_are_nonexclusive_and_strand_correct() {
    let plus = model(&[(10, 100), (200, 250), (300, 350), (400, 500)], '+');
    for (skip, label) in [
        ((100, 210), "alternative_acceptor"),
        ((90, 200), "alternative_donor"),
        ((100, 300), "exon_skip"),
        ((90, 200), "alternative_first_junction"),
        ((350, 450), "alternative_terminal_acceptor"),
        ((110, 220), "unannotated_both_ends"),
    ] {
        assert!(
            geometry_events(skip, &plus).contains(&label.to_string()),
            "{skip:?}: {label}"
        );
    }
    let minus = model(&plus.genomic_exons, '-');
    assert!(geometry_events((100, 210), &minus).contains(&"alternative_donor".to_string()));
    assert!(
        geometry_events((110, 200), &minus).contains(&"alternative_terminal_acceptor".to_string())
    );
    assert!(
        geometry_events((350, 450), &minus).contains(&"alternative_first_junction".to_string())
    );
}

#[test]
fn junction_assessment_file_budgets_fail_without_final_report() {
    for field in ["max_input_bytes", "max_total_work", "max_output_bytes"] {
        let root = tempfile::tempdir().unwrap();
        let (e, r) = fixture(root.path());
        let mut v = serde_json::to_value(r).unwrap();
        v["budgets"][field] = json!(1);
        let r: AssessObservedSpliceJunctionsRequest = serde_json::from_value(v).unwrap();
        let out = PathBuf::from(&r.output_dir);
        assert!(
            e.assess_observed_splice_junctions(r, &mut |_| true)
                .is_err(),
            "{field}"
        );
        assert!(!out.join("report.json").exists());
    }
}

#[test]
fn junction_assessment_manifest_hierarchy_is_functional_not_guessed() {
    let root = tempfile::tempdir().unwrap();
    let (e, mut r) = fixture(root.path());
    let context = ObservedSampleContext {
        sample_id: "invented_sample".into(),
        barcode: "invented_barcode".into(),
        priming: "invented".into(),
        file_id: "file_0000".into(),
    };
    let row = ObservedSampleBinding {
        context,
        library_id: "L".into(),
        cdna_prep_id: "P1".into(),
        biological_sample_id: "S".into(),
        cell_line: None,
    };
    let mut other = row.clone();
    other.context.file_id = "unused_file".into();
    other.cdna_prep_id = "P2".into();
    let p = root.path().join("sample hierarchy.json");
    fs::write(
        &p,
        serde_json::to_vec(&ObservedSampleManifest {
            schema: OBSERVED_SAMPLE_MANIFEST_SCHEMA.into(),
            rows: vec![row, other],
        })
        .unwrap(),
    )
    .unwrap();
    r.sample_manifest = Some(super::super::tests::bound_file(&p));
    assert!(
        e.assess_observed_splice_junctions(r, &mut |_| true)
            .unwrap_err()
            .message
            .contains("hierarchy")
    );
}

#[test]
fn junction_assessment_cancellation_and_no_overwrite_preserve_content() {
    let root = tempfile::tempdir().unwrap();
    let (e, mut r) = fixture(root.path());
    r.cancel_path = Some(root.path().join("cancel marker").display().to_string());
    GentleEngine::cancel_observed_splice_junctions(&r).unwrap();
    let before = encode(e.state()).unwrap();
    let n = e.journal_len();
    let result = e
        .assess_observed_splice_junctions(r.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(result.status, "cancelled");
    assert!(result.counts.is_none());
    assert!(
        e.assess_observed_splice_junctions(r, &mut |_| true)
            .is_err()
    );
    assert_eq!(before, encode(e.state()).unwrap());
    assert_eq!(n, e.journal_len());
}

#[test]
fn junction_assessment_source_hash_mismatch_is_not_empty_or_novel() {
    let root = tempfile::tempdir().unwrap();
    let (e, mut r) = fixture(root.path());
    let ObservedJunctionSource::S1Report { report, .. } = &mut r.source;
    report.sha256 = format!("sha256:{}", "0".repeat(64));
    assert!(
        e.assess_observed_splice_junctions(r, &mut |_| true)
            .is_err()
    );
}

#[test]
fn junction_assessment_shared_shell_uses_one_typed_request_and_lifecycle() {
    use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};
    let root = tempfile::tempdir().unwrap();
    let (mut e, r) = fixture(root.path());
    let path = root.path().join("request's path.json");
    fs::write(
        &path,
        serde_json::to_string_pretty(&r)
            .unwrap()
            .replace('\n', "\r\n"),
    )
    .unwrap();
    let cmd = parse_shell_line(&format!(
        "rna-reads splice-junctions-assess {}",
        quote_shell_arg(&format!("@{}", path.display()))
    ))
    .unwrap();
    let before = encode(e.state()).unwrap();
    let run = execute_shell_command(&mut e, &cmd).unwrap();
    assert!(!run.state_changed);
    assert_eq!(before, encode(e.state()).unwrap());
    assert_eq!(
        run.output["result"]["observed_splice_junctions"]["status"],
        "completed"
    );
}

#[test]
fn junction_assessment_declared_excerpt_lf_crlf_edges_and_softmask() {
    for (i, bases, declared_len, expected) in [
        (0, "gtcag", 5, ObservedMotifClass::Other),
        (1, "GTCA", 4, ObservedMotifClass::NotAssessable),
        (2, "GNcAG", 5, ObservedMotifClass::NotAssessable),
    ] {
        let root = tempfile::tempdir().unwrap();
        let (e, mut r) = fixture(root.path());
        let ObservedJunctionSource::S1Report { report, .. } = &r.source;
        let saved: ObservedTranscriptStructuresReport =
            serde_json::from_slice(&fs::read(&report.path).unwrap()).unwrap();
        let d: ObservedTranscriptStructuresData = serde_json::from_slice(
            &fs::read(Path::new(&saved.request.output_dir).join("artifacts/structures.json"))
                .unwrap(),
        )
        .unwrap();
        let family = &d.analysis.unwrap().families[0];
        let path = root.path().join("excerpt.fa");
        fs::write(&path, format!(">{}\n{bases}\n", family.contig)).unwrap();
        r.reference_sequence = Some(ObservedJunctionReference {
            file: super::super::tests::bound_file(&path),
            format: ObservedJunctionReferenceFormat::FastaPlainV1,
            authority: ObservedJunctionReferenceAuthority::DeclaredExcerpt {
                contig: family.contig.clone(),
                offset_0based: family.skips[0].0,
                length: declared_len,
            },
        });
        let result = e
            .assess_observed_splice_junctions(r.clone(), &mut |_| true)
            .unwrap();
        let one = data(&result);
        assert_eq!(one.junctions[0].motif.motif_class, expected, "case {i}");
        assert_eq!(one.junctions[0].motif.reference_authority, "declared");
        let lf_hash = r.reference_sequence.as_ref().unwrap().file.sha256.clone();
        fs::write(&path, format!(">{}\r\n{bases}\r\n", family.contig)).unwrap();
        r.reference_sequence.as_mut().unwrap().file = super::super::tests::bound_file(&path);
        assert_ne!(lf_hash, r.reference_sequence.as_ref().unwrap().file.sha256);
        r.output_dir = root.path().join("crlf result").display().to_string();
        let two = data(
            &e.assess_observed_splice_junctions(r, &mut |_| true)
                .unwrap(),
        );
        assert_eq!(
            serde_json::to_value(&one.junctions[0].motif).unwrap(),
            serde_json::to_value(&two.junctions[0].motif).unwrap()
        );
    }
}

#[test]
fn junction_assessment_full_assembly_authority_requires_exact_raw_hash() {
    let root = tempfile::tempdir().unwrap();
    let (e, mut r) = fixture(root.path());
    let ObservedJunctionSource::S1Report { report, .. } = &r.source;
    let saved: ObservedTranscriptStructuresReport =
        serde_json::from_slice(&fs::read(&report.path).unwrap()).unwrap();
    r.reference_sequence = Some(ObservedJunctionReference {
        file: saved.request.annotation_sources[0].reference.fasta.clone(),
        format: ObservedJunctionReferenceFormat::FastaPlainV1,
        authority: ObservedJunctionReferenceAuthority::FullAssembly,
    });
    let d = data(
        &e.assess_observed_splice_junctions(r.clone(), &mut |_| true)
            .unwrap(),
    );
    assert_eq!(d.junctions[0].motif.reference_authority, "verified");
    let newpath = root.path().join("wrong assembly.fa");
    fs::write(&newpath, b">other\nACGT\n").unwrap();
    r.reference_sequence.as_mut().unwrap().file = super::super::tests::bound_file(&newpath);
    r.output_dir = root.path().join("refused reference").display().to_string();
    assert!(
        e.assess_observed_splice_junctions(r, &mut |_| true)
            .is_err()
    );
}

#[test]
fn junction_assessment_missing_catalog_remains_unassessed() {
    let root = tempfile::tempdir().unwrap();
    let (e, r) = fixture(root.path());
    let ObservedJunctionSource::S1Report { report, .. } = &r.source;
    let saved: ObservedTranscriptStructuresReport =
        serde_json::from_slice(&fs::read(&report.path).unwrap()).unwrap();
    fs::remove_file(
        &saved.request.annotation_sources[0]
            .reference
            .annotation
            .path,
    )
    .unwrap();
    let d = data(
        &e.assess_observed_splice_junctions(r, &mut |_| true)
            .unwrap(),
    );
    assert!(
        d.junctions
            .iter()
            .all(|j| j.sources.iter().all(|s| s.relation == "not_assessable"))
    );
}
