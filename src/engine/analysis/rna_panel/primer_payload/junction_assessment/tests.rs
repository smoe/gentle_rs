//! Inline synthetic S1/Phase 2 fixtures; no downloaded RNA/reference or committed artifacts.

use super::*;

#[path = "acceptance_tests.rs"]
mod acceptance;

impl GentleEngine {
    pub(crate) fn observed_junctions_test_fixture(
        root: &Path,
    ) -> (Self, AssessObservedSpliceJunctionsRequest) {
        fixture(root)
    }
}

fn loaded(r: &AssessObservedSpliceJunctionsRequest) -> source::Loaded {
    let mut m = Meter {
        b: &r.budgets,
        input: 0,
        work: 0,
        comparisons: 0,
    };
    let mut progress = |_| true;
    let mut ctl = Control::new(None, "synthetic".into(), 1, &mut progress).unwrap();
    source::load(r, &mut m, &mut ctl).unwrap()
}

fn index(
    s: &s1::Summary,
    b: &ObservedJunctionBudgets,
) -> BTreeMap<(String, s1::Interval), Vec<ObservedJunctionRecordSupport>> {
    let mut m = Meter {
        b,
        input: 0,
        work: 0,
        comparisons: 0,
    };
    let mut progress = |_| true;
    let mut ctl = Control::new(None, "synthetic".into(), 1, &mut progress).unwrap();
    support::record_index(s, &BTreeMap::new(), &mut m, &mut ctl).unwrap()
}

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

#[test]
fn junction_assessment_records_are_deduplicated_and_never_bridge_missing_skips() {
    let td = tempfile::tempdir().unwrap();
    let (_, req) = fixture(td.path());
    let mut s = loaded(&req).data.analysis.unwrap();
    let one = s.observations[0].clone();
    let r = s
        .records
        .iter_mut()
        .find(|r| r.record == one.record)
        .unwrap();
    r.hypotheses = vec![
        vec![one.observation_id.clone()],
        vec![one.observation_id.clone()],
    ];
    r.search_declared_complete = true;
    r.all_query_views_assessed = true;
    r.definite_chains = r.possible_chains.clone();
    let a = index(&s, &req.budgets);
    assert_eq!(a.len(), 1);
    assert_eq!(a.values().next().unwrap().len(), s.records.len());
    assert_eq!(
        support::counts(a.values().next().unwrap())
            .definite
            .qualified_only,
        1
    );
    s.records[0].search_declared_complete = false;
    let b = index(&s, &req.budgets);
    let counts = support::counts(b.values().next().unwrap());
    assert_eq!(counts.definite, ObservedQualifiedCounts::default());
    assert!(counts.possible_only.qualified_only > 0);
    // Two separate physical records with no explicit N never manufacture a junction.
    for f in &mut s.families {
        f.skips.clear();
    }
    for o in &mut s.observations {
        o.skips.clear();
        o.flank_qualification.clear();
    }
    assert!(index(&s, &req.budgets).is_empty());
}

fn synthetic_support(
    n: usize,
    certainty: &str,
    qualification: &str,
) -> (Vec<ObservedJunctionRecordSupport>, source::Manifest) {
    let mut records = Vec::new();
    let mut manifest = BTreeMap::new();
    for i in 0..n {
        let context = ObservedSampleContext {
            sample_id: format!("context{i}"),
            barcode: String::new(),
            priming: "synthetic".into(),
            file_id: format!("file{i}"),
        };
        manifest.insert(
            context.clone(),
            ObservedSampleBinding {
                context: context.clone(),
                library_id: format!("L{i}"),
                cdna_prep_id: "shared_prep".into(),
                biological_sample_id: "one_biological_sample".into(),
                cell_line: None,
            },
        );
        records.push(ObservedJunctionRecordSupport {
            record: s1::ReadRecord {
                source_sha256: format!("sha256:{}", "a".repeat(64)),
                record_index: i as u64,
                read_id: format!("r{i}"),
            },
            certainty: certainty.into(),
            qualification: qualification.into(),
            contexts: vec![context],
        });
    }
    (records, manifest)
}

#[test]
fn junction_assessment_recurrence_needs_qualified_explicit_independent_units() {
    let policy = ObservedRecurrencePolicy {
        unit: ObservedRecurrenceUnit::Library,
        min_units: 2,
        certainty: ObservedRecurrenceCertainty::Definite,
        require_qualified_flanks: true,
    };
    let (mut rows, manifest) = synthetic_support(3, "definite", "qualified_only");
    assert_eq!(
        support::recurrence(&policy, &rows, &manifest, false).status,
        "recurrent_candidate"
    );
    let mut prep_policy = policy.clone();
    prep_policy.unit = ObservedRecurrenceUnit::CdnaPrep;
    assert_eq!(
        support::recurrence(&prep_policy, &rows, &manifest, false).status,
        "not_recurrent_under_policy"
    );
    for q in ["insufficient_only", "mixed", "unassessed_only"] {
        for r in &mut rows {
            r.qualification = q.into();
        }
        assert_ne!(
            support::recurrence(&policy, &rows, &manifest, false).status,
            "recurrent_candidate"
        );
    }
    for r in &mut rows {
        r.qualification = "qualified_only".into();
        r.certainty = "possible".into();
    }
    assert_eq!(
        support::recurrence(&policy, &rows, &manifest, false).status,
        "not_assessable"
    );
    let mut possible = policy.clone();
    possible.certainty = ObservedRecurrenceCertainty::DefiniteOrPossible;
    assert_eq!(
        support::recurrence(&possible, &rows, &manifest, false).status,
        "recurrent_candidate"
    );
    assert_eq!(
        support::recurrence(&possible, &rows, &BTreeMap::new(), false).status,
        "not_assessable"
    );
    let other_context = rows[1].contexts[0].clone();
    rows[0].contexts.push(other_context);
    assert_ne!(
        support::recurrence(&possible, &rows[..1], &manifest, false).status,
        "recurrent_candidate"
    );
    let grouped = support::groups(&rows, &manifest);
    let biological = grouped
        .iter()
        .find(|g| g.dimension == "biological_sample")
        .unwrap();
    assert_eq!(biological.counts.possible_only.qualified_only, 3);
}

#[test]
fn junction_assessment_sources_stay_separate_novel_combinations_and_strands() {
    let a = model(&[(10, 100), (200, 250)], '+');
    let mut b = model(&[(10, 150), (300, 350)], '+');
    b.transcript_id = "second".into();
    let c = source::Catalog {
        source_id: "ensembl_like".into(),
        gene_id: Some("G".into()),
        reason: "resolved".into(),
        models: vec![a.clone(), b.clone()],
    };
    let row = source::compare_source("synthetic_contig", (100, 300), &c, 5);
    assert_eq!(row.relation, "none_within_window");
    assert!(row.events.contains(&"novel_combination".into()));
    let mut d = source::Catalog {
        source_id: "refseq_like".into(),
        gene_id: Some("G".into()),
        reason: "resolved".into(),
        models: vec![model(&[(10, 100), (300, 350)], '-')],
    };
    let exact = source::compare_source("synthetic_contig", (100, 300), &d, 5);
    assert_eq!(exact.relation, "exact");
    assert_eq!(assigned_orientation(&[row.clone(), exact]), None);
    d.reason = "ambiguous_gene".into();
    assert_eq!(
        source::compare_source("synthetic_contig", (100, 300), &d, 5).relation,
        "not_assessable"
    );
    d.reason = "missing_gene".into();
    assert_eq!(
        source::compare_source("synthetic_contig", (100, 300), &d, 5).relation,
        "not_assessable"
    );
    let mut reordered = source::Catalog {
        models: vec![b, a],
        ..c
    };
    let first = source::compare_source("synthetic_contig", (101, 201), &reordered, 5);
    reordered.models.reverse();
    let second = source::compare_source("synthetic_contig", (101, 201), &reordered, 5);
    assert_eq!(
        encode(&first.near_hits).unwrap(),
        encode(&second.near_hits).unwrap()
    );
    assert_eq!(encode(&first).unwrap(), encode(&second).unwrap());
    assert_eq!(row.source_id, "ensembl_like");
}

#[test]
fn junction_assessment_full_identity_survives_pinned_short_prefix_collision() {
    // Offline bounded search of start=1..300000, end=start+17 using canonical JSON
    // SHA-256 over (domain, reference, contig, interval). No biological reference.
    let reference = s1::ReferenceBinding {
        assembly: "synthetic_collision_v1".into(),
        fasta_sha256: format!("sha256:{}", "a".repeat(64)),
    };
    let a = junction_id(&reference, "synthetic_contig", (57505, 57522)).unwrap();
    let b = junction_id(&reference, "synthetic_contig", (62716, 62733)).unwrap();
    assert_eq!(
        a,
        "sha256:69d1a3a75dbce6798c73f0fd4170720390de938a11daaaf377fbd489a0ced9a1"
    );
    assert_eq!(
        b,
        "sha256:69d1a3a7aa953d3a96811e15d911271e7d703378a3099f5f318f2778711dc73b"
    );
    assert_eq!(&a[7..15], &b[7..15]);
    assert_ne!(a, b);
}

#[test]
fn junction_assessment_unknown_contexts_and_zero_gene_coherent_denominators() {
    let td = tempfile::tempdir().unwrap();
    let (_, req) = fixture(td.path());
    let l = loaded(&req);
    let mut s = l.data.analysis.unwrap();
    let mut contexts = l.contexts;
    let mut manifest = BTreeMap::new();
    for (record, cs) in &contexts {
        for c in cs {
            manifest.insert(
                c.clone(),
                ObservedSampleBinding {
                    context: c.clone(),
                    library_id: "L".into(),
                    cdna_prep_id: "P".into(),
                    biological_sample_id: "S".into(),
                    cell_line: None,
                },
            );
        }
        assert!(s.records.iter().any(|r| &r.record == record));
    }
    for o in &mut s.observations {
        for c in &mut o.comparisons {
            for m in &mut c.models {
                m.relation = s1::ModelRelation::OutsideModel;
            }
        }
    }
    let d = support::metadata(&s, &contexts, &manifest, &l.report, &l.data.unavailable);
    assert_eq!(d.libraries.len(), 1);
    assert!(d.libraries[0].total_records > 0);
    assert_eq!(d.libraries[0].gene_coherent_records, 0);
    contexts.clear();
    let d = support::metadata(&s, &contexts, &manifest, &l.report, &l.data.unavailable);
    assert!(!d.unassessed.is_empty());
    assert!(!d.unused.is_empty());
    assert!(d.libraries.is_empty());
}

#[test]
fn junction_assessment_publication_is_report_last_and_discloses_unowned_orphans() {
    for orphan in [false, true] {
        let td = tempfile::tempdir().unwrap();
        let (e, req) = fixture(td.path());
        let report = e
            .assess_observed_splice_junctions(req, &mut |_| true)
            .unwrap();
        let d = data(&report);
        let mut r = report.clone();
        r.request.output_dir = td.path().join("faulted output").display().to_string();
        let stage = staging(&r.request.output_dir).unwrap();
        r.artifacts = publication::artifacts(stage.path(), &d, &r).unwrap();
        let out = PathBuf::from(&r.request.output_dir);
        let err = publication::publish(stage, &r, &mut |phase| {
            assert!(!out.join("report.json").exists());
            if matches!(phase, publication::Phase::ReportReady) {
                assert!(out.join("artifacts/assessment.json").exists());
                if orphan {
                    fs::write(out.join("unowned sentinel"), "preserve").unwrap();
                }
                return Err(bad("injected_publication_fault"));
            }
            Ok(())
        })
        .unwrap_err();
        assert_eq!(out.exists(), orphan);
        assert_eq!(
            err.message.contains("junction_assessment_orphaned_output"),
            orphan
        );
    }
}

#[test]
fn junction_assessment_s1_tampering_and_noncompleted_sources_do_not_become_counts() {
    for unavailable in [false, true] {
        let td = tempfile::tempdir().unwrap();
        let (e, mut r) = fixture(td.path());
        let ObservedJunctionSource::S1Report {
            result,
            report: file,
        } = &mut r.source;
        let mut saved: ObservedTranscriptStructuresReport =
            serde_json::from_slice(&fs::read(&file.path).unwrap()).unwrap();
        if unavailable {
            saved.status = "cancelled".into();
            saved.counts = None;
            saved.content_sha256.clear();
            saved.content_sha256 = hash(&saved).unwrap();
            fs::write(&file.path, serde_json::to_vec_pretty(&saved).unwrap()).unwrap();
            *file = super::super::tests::bound_file(Path::new(&file.path));
            *result =
                GentleEngine::observed_transcript_structures_result_binding(&saved.request, &saved)
                    .unwrap();
            let report = e
                .assess_observed_splice_junctions(r, &mut |_| true)
                .unwrap();
            assert_eq!(report.status, "unavailable");
            assert!(report.counts.is_none());
            assert!(data(&report).junctions.is_empty());
        } else {
            let path = Path::new(&saved.request.output_dir).join("artifacts/structures.json");
            fs::write(path, "{}").unwrap();
            assert!(
                e.assess_observed_splice_junctions(r, &mut |_| true)
                    .is_err()
            );
        }
    }
}

#[test]
fn junction_assessment_every_budget_accepts_boundary_and_rejects_one_over() {
    let td = tempfile::tempdir().unwrap();
    let (e, base) = acceptance::cd274_like(td.path(), false);
    let report = e
        .assess_observed_splice_junctions(base.clone(), &mut |_| true)
        .unwrap();
    let d = data(&report);
    let c = report.counts.as_ref().unwrap();
    let max_hits = d
        .junctions
        .iter()
        .map(|j| j.sources.iter().map(|s| s.near_hits.len()).sum::<usize>())
        .max()
        .unwrap();
    let boundaries = [
        ("max_input_bytes", c.charged_input_bytes),
        ("max_total_work", c.charged_work),
        ("max_distinct_junctions", c.distinct_junctions as u64),
        (
            "max_annotated_junctions_compared",
            c.annotated_comparisons as u64,
        ),
        ("max_near_hits_per_junction", max_hits as u64),
        (
            "max_reference_bytes",
            base.reference_sequence.as_ref().unwrap().file.bytes,
        ),
        (
            "max_manifest_bytes",
            base.sample_manifest.as_ref().unwrap().bytes,
        ),
    ];
    for (i, (name, limit)) in boundaries.into_iter().enumerate() {
        assert!(
            limit > 1,
            "{name} must exercise a positive one-over boundary"
        );
        for exact in [true, false] {
            let mut v = serde_json::to_value(&base).unwrap();
            v["budgets"][name] = json!(if exact { limit } else { limit - 1 });
            v["output_dir"] = json!(
                td.path()
                    .join(format!("budget-{i}-{exact}"))
                    .display()
                    .to_string()
            );
            let r: AssessObservedSpliceJunctionsRequest = serde_json::from_value(v).unwrap();
            let path = PathBuf::from(&r.output_dir);
            let result = e.assess_observed_splice_junctions(r, &mut |_| true);
            assert_eq!(result.is_ok(), exact, "{name} {exact}: {result:?}");
            assert_eq!(path.join("report.json").exists(), exact);
        }
    }
    // Output covers the final report as well as all five artifacts, not just rows.
    for exact in [true, false] {
        let mut r = report.clone();
        r.request.output_dir = td
            .path()
            .join(format!("output-bound-{exact}"))
            .display()
            .to_string();
        let stage = staging(&r.request.output_dir).unwrap();
        r.artifacts = publication::artifacts(stage.path(), &d, &r).unwrap();
        let artifacts = r.artifacts.iter().map(|a| a.bytes).sum::<u64>();
        loop {
            let total = artifacts + serde_json::to_vec_pretty(&r).unwrap().len() as u64;
            if r.request.budgets.max_output_bytes == total {
                break;
            }
            r.request.budgets.max_output_bytes = total;
        }
        if !exact {
            r.request.budgets.max_output_bytes -= 1;
        }
        let out = PathBuf::from(&r.request.output_dir);
        assert_eq!(
            publication::publish(stage, &r, &mut |_| Ok(())).is_ok(),
            exact
        );
        assert_eq!(out.join("report.json").exists(), exact);
    }
}

#[test]
fn junction_assessment_raw_lifecycle_and_inspection_refusal_preserve_content_and_history() {
    use crate::engine::operation_lifecycle::LifecycleRequest;
    let td = tempfile::tempdir().unwrap();
    let (_, r) = fixture(td.path());
    let op = Operation::AssessObservedSpliceJunctions {
        request: Box::new(r.clone()),
    };
    let state = crate::engine::lifecycle_slice3_tests::state();
    let mut e = GentleEngine::from_state_for(state.clone(), LifecycleRequest::Operation(&op));
    let before = encode(e.state()).unwrap();
    let history = e.history_summary();
    let result = e.apply(op).unwrap();
    let report = result.observed_splice_junctions.unwrap();
    assert_eq!(encode(e.state()).unwrap(), before);
    let h = e.history_summary();
    assert_eq!(
        (h.undo_count, h.redo_count),
        (history.undo_count, history.redo_count)
    );
    assert_eq!(h.operation_log_count, history.operation_log_count + 1);
    let before_log = e.journal_len();
    let cmd = crate::engine_shell::parse_shell_line(&format!(
        "rna-reads splice-junctions-show {} --rows 1 --verify-artifacts",
        crate::engine_shell::quote_shell_arg(
            &Path::new(&r.output_dir)
                .join("report.json")
                .display()
                .to_string()
        )
    ))
    .unwrap();
    assert!(
        !crate::engine_shell::execute_shell_command(&mut e, &cmd)
            .unwrap()
            .state_changed
    );
    assert_eq!(e.journal_len(), before_log);
    assert!(
        e.apply(Operation::AssessObservedSpliceJunctions {
            request: Box::new(r)
        })
        .is_err()
    );
    assert_eq!(e.journal_len(), before_log);
    assert_eq!(encode(e.state()).unwrap(), before);
    assert_eq!(report.status, "completed");
}

#[test]
fn junction_assessment_agent_outcomes_do_not_turn_unavailable_into_success() {
    use crate::agent_feedback::AgentExecutionStatus as Status;
    for (outcome, expected) in [
        ("completed", Status::Completed),
        ("unavailable", Status::Partial),
        ("cancelled", Status::Cancelled),
        ("control_error", Status::Failed),
    ] {
        let output = json!({"schema":"gentle.command_observed_junction_result.v1","observed_junctions_result":{"schema":OBSERVED_JUNCTION_BINDING_SCHEMA,"outcome":outcome}});
        assert_eq!(Status::from_shell_output(&output), expected);
        let output = json!({"result":{"observed_splice_junctions":{"schema":OBSERVED_JUNCTION_REPORT_SCHEMA,"status":outcome}}});
        assert_eq!(Status::from_shell_output(&output), expected);
    }
}

#[cfg(feature = "js-interface")]
#[test]
fn junction_assessment_js_generic_operation_runtime() {
    let td = tempfile::tempdir().unwrap();
    let (_, req) = fixture(td.path());
    let op = Operation::AssessObservedSpliceJunctions {
        request: Box::new(req),
    };
    let mut js = crate::js_interface::JavaScriptInterface::new();
    js.run_checked(format!("const result=apply_operation({},{}); if(result.result.observed_splice_junctions.status!=='completed') throw Error('missing assessment');",serde_json::to_string(GentleEngine::new().state()).unwrap(),serde_json::to_string(&op).unwrap())).unwrap();
}

#[cfg(feature = "lua-interface")]
#[test]
fn junction_assessment_lua_generic_operation_runtime() {
    use mlua::LuaSerdeExt;
    let td = tempfile::tempdir().unwrap();
    let (_, req) = fixture(td.path());
    let op = Operation::AssessObservedSpliceJunctions {
        request: Box::new(req),
    };
    let lua = crate::lua_interface::LuaInterface::new();
    lua.register_rust_functions().unwrap();
    lua.lua()
        .globals()
        .set(
            "input",
            lua.lua().to_value(GentleEngine::new().state()).unwrap(),
        )
        .unwrap();
    lua.lua()
        .globals()
        .set("request", serde_json::to_string(&op).unwrap())
        .unwrap();
    lua.lua().load("local result=apply_operation(input,request); assert(result.result.observed_splice_junctions.status=='completed')").exec().unwrap();
}
