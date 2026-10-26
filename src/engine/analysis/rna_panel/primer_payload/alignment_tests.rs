//! Synthetic PAF/annotation oracles and temporary artifacts; no aligner or real sequencing input.

use super::tests::{PAYLOAD, limits, prepare, source};
use super::*;

fn file(path: &Path) -> PrimerPayloadFile {
    let bytes = fs::read(path).unwrap();
    PrimerPayloadFile {
        path: path.display().to_string(),
        sha256: sha256_prefixed_bytes(&bytes),
        bytes: bytes.len() as u64,
    }
}

fn request(
    root: &Path,
    trace: RnaPrimerTraceResultBinding,
    prepared: &PrimerPayloadQueries,
) -> SummarizePrimerPayloadAlignmentsRequest {
    let fasta = root.join("synthetic reference.fa");
    fs::write(
        &fasta,
        format!(
            ">toy_expected\n{}\n>toy_foreign\nAGTCCGATCA{PAYLOAD}\n",
            "C".repeat(PAYLOAD.len() + 10)
        ),
    )
    .unwrap();
    let reference_file = file(&fasta);
    let annotation = PrimerPayloadAnnotation {
        schema: PRIMER_PAYLOAD_ANNOTATION_SCHEMA.into(),
        reference_id: "synthetic_reference".into(),
        assembly: "invented_assembly".into(),
        release: "invented_release".into(),
        reference_fasta_sha256: reference_file.sha256.clone(),
        coordinate_convention: "zero_based_half_open".into(),
        transcripts: [
            ("toy_expected", "expected_gene"),
            ("toy_foreign", "foreign_gene"),
        ]
        .into_iter()
        .map(|(tx, gene)| PrimerPayloadTranscript {
            transcript_id: tx.into(),
            gene_ids: vec![gene.into()],
            target_name: tx.into(),
            target_length: PAYLOAD.len() + 10,
            genomic_contig: None,
            genomic_strand: None,
            genomic_exons: Vec::new(),
        })
        .collect(),
    };
    let annotation_path = root.join("annotation.json");
    fs::write(&annotation_path, serde_json::to_vec(&annotation).unwrap()).unwrap();
    let paf = root.join("query grouped.paf");
    fs::write(&paf,prepared.queries.iter().map(|q|format!("{}\t{}\t0\t{}\t+\ttoy_foreign\t{}\t10\t{}\t{}\t{}\t60\tcg:Z:{}M\tcs:Z::{}\tNM:i:0\ttp:A:P\n",q.query_id,q.query_length,q.query_length,PAYLOAD.len()+10,q.query_length+10,q.query_length,q.query_length,q.query_length,q.query_length)).collect::<String>()).unwrap();
    SummarizePrimerPayloadAlignmentsRequest {
        schema: PRIMER_PAYLOAD_REQUEST_SCHEMA.into(),
        trace,
        queries: file(&Path::new(&prepared.request.output_dir).join("report.json")),
        query_manifest_sha256: prepared.content_sha256.clone(),
        references: vec![PrimerPayloadReference {
            reference_id: annotation.reference_id.clone(),
            scope: PrimerPayloadReferenceScope::Transcriptome,
            assembly: annotation.assembly,
            release: annotation.release,
            fasta: reference_file,
            annotation: file(&annotation_path),
        }],
        expected_targets: vec![PrimerPayloadExpectedTarget {
            primer_id: "invented_primer".into(),
            reference_id: annotation.reference_id.clone(),
            phase1_locus_id: "invented_expected_locus".into(),
            phase1_transcript_feature_id: None,
            gene_ids: vec!["expected_gene".into()],
            transcript_ids: vec!["toy_expected".into()],
        }],
        batches: vec![PrimerPayloadAlignmentBatch {
            batch_id: "synthetic_batch".into(),
            reference_id: annotation.reference_id,
            paf: file(&paf),
            query_manifest_sha256: prepared.content_sha256.clone(),
            producer: "literal PAF fixture".into(),
            producer_version: "1".into(),
            arguments: Vec::new(),
            producer_exit_code: 0,
            declared_complete: true,
            competitor_search_limitations:
                "Invented rows; no aligner/search and no completeness proof".into(),
            accounting: prepared
                .queries
                .iter()
                .map(|q| (q.query_id.clone(), PrimerPayloadQueryAccounting::Mapped))
                .collect(),
        }],
        policy: PrimerPayloadPolicy {
            strict_alignment_detail: true,
            min_segment_bases: 12,
            min_identity_per_mille: 800,
            min_mapq: 1,
            dominant_fraction_per_mille: 800,
            min_cross_locus_segment_bases: 12,
            breakpoint_tolerance_bases: 2,
            terminal_seed_lengths: vec![6, 8, 10, 12],
        },
        budgets: limits(),
        output_dir: root.join("new payload evidence").display().to_string(),
        cancel_path: None,
    }
}

#[test]
fn primer_payload_alignment_round_trip_keeps_phase1_bytes_and_separates_partner_edges() {
    let td = tempfile::tempdir().unwrap();
    let (engine, trace) = source(td.path(), true, PAYLOAD);
    let before = fs::read(&trace.report_path).unwrap();
    let state = encode(engine.state()).unwrap();
    let journal = engine.journal_len();
    let prepared = engine
        .prepare_primer_payload_alignment_queries(prepare(trace.clone(), td.path()), &mut |_| true)
        .unwrap();
    let req = request(td.path(), trace.clone(), &prepared);
    let report = engine
        .summarize_primer_payload_alignments(req.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(report.status, "completed");
    assert_eq!(report.evidence_rows, 1);
    assert_eq!(report.profiles[0].supporting_reads, 1);
    assert_eq!(report.profiles[0].partner_gene_ids, vec!["foreign_gene"]);
    assert_eq!(report.breakpoint_clusters[0].minimum_edge, 10);
    let binding = GentleEngine::primer_payload_alignment_result_binding(&req, &report).unwrap();
    let inspect =
        GentleEngine::inspect_primer_payload(&binding.report_path, Some(&binding), 1, true)
            .unwrap();
    assert!(inspect.artifact_content_verified);
    assert_eq!(inspect.rows[0]["classification"], "off_target_dominant");
    assert_eq!(inspect.rows[0]["expected_target_bases"], 0);
    assert_eq!(inspect.rows[0]["cross_locus_payload_evidence"], false);
    let seeds = inspect.rows[0]["terminal_seeds"]
        .as_object()
        .unwrap()
        .values()
        .next()
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(seeds[0]["mismatches"], 0);
    assert!(seeds[3]["mismatches"].is_null());
    assert_eq!(fs::read(&trace.report_path).unwrap(), before);
    assert_eq!(encode(engine.state()).unwrap(), state);
    assert_eq!(engine.journal_len(), journal);
    assert!(
        engine
            .summarize_primer_payload_alignments(req, &mut |_| true)
            .is_err()
    );
}

#[test]
fn primer_payload_alignment_refuses_omitted_accounting_and_hash_mutation_without_publication() {
    let td = tempfile::tempdir().unwrap();
    let (engine, trace) = source(td.path(), true, PAYLOAD);
    let prepared = engine
        .prepare_primer_payload_alignment_queries(prepare(trace.clone(), td.path()), &mut |_| true)
        .unwrap();
    let mut req = request(td.path(), trace, &prepared);
    let original = req.batches[0].accounting.clone();
    req.batches[0].accounting.clear();
    assert!(
        engine
            .summarize_primer_payload_alignments(req.clone(), &mut |_| true)
            .unwrap_err()
            .message
            .contains("accounting_coverage")
    );
    assert!(!Path::new(&req.output_dir).exists());
    req.batches[0].accounting = original;
    let path = &req.batches[0].paf.path;
    let mut bytes = fs::read(path).unwrap();
    bytes[0] = b'x';
    fs::write(path, bytes).unwrap();
    assert!(
        engine
            .summarize_primer_payload_alignments(req.clone(), &mut |_| true)
            .is_err()
    );
    assert!(!Path::new(&req.output_dir).exists());
}

#[test]
fn primer_payload_explicit_unmapped_is_not_missing_or_cross_locus() {
    let td = tempfile::tempdir().unwrap();
    let (engine, trace) = source(td.path(), true, PAYLOAD);
    let prepared = engine
        .prepare_primer_payload_alignment_queries(prepare(trace.clone(), td.path()), &mut |_| true)
        .unwrap();
    let mut req = request(td.path(), trace, &prepared);
    fs::write(&req.batches[0].paf.path, b"").unwrap();
    req.batches[0].paf = file(Path::new(&req.batches[0].paf.path));
    req.batches[0]
        .accounting
        .values_mut()
        .for_each(|x| *x = PrimerPayloadQueryAccounting::ProducerReportedUnmapped);
    let report = engine
        .summarize_primer_payload_alignments(req.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(report.profiles[0].unmapped_reads, 1);
    assert_eq!(report.profiles[0].supporting_reads, 0);
    assert!(report.breakpoint_clusters.is_empty());
    let inspected = GentleEngine::inspect_primer_payload(
        &Path::new(&req.output_dir)
            .join("report.json")
            .display()
            .to_string(),
        None,
        1,
        false,
    )
    .unwrap();
    assert!(!inspected.artifact_content_verified);
    assert_eq!(
        inspected.rows[0]["classification"],
        "producer_reported_unmapped"
    );
}

fn q() -> PrimerPayloadQuery {
    PrimerPayloadQuery {
        query_id: "q_synthetic".into(),
        anchor_id: "synthetic_anchor".into(),
        file_id: "invented_file".into(),
        record_index: 0,
        read_id: "duplicate read name".into(),
        sample_id: "sample".into(),
        barcode: "arbitrary".into(),
        priming: "unknown".into(),
        primer_ids: vec!["invented_primer".into()],
        anchor: RnaPrimerHit {
            primer_id: Some("invented_primer".into()),
            anchor_class: "full_exact".into(),
            strand: PrimerReadStrand::Forward,
            start: 0,
            end: 10,
            edits: 0,
            core_start: Some(0),
            core_end: Some(10),
            core_edits: Some(0),
            terminal: true,
        },
        raw_read_length: 60,
        raw_payload_interval: (10, 60),
        excluded_query_intervals: Vec::new(),
        query_length: 50,
        usable_bases: 50,
        sequence_sha256: sha256_prefixed_str(&"A".repeat(50)),
    }
}

#[test]
fn primer_payload_dominance_and_alternative_chains_do_not_invent_cross_locus_reads() {
    for case in [
        "expected",
        "foreign",
        "mapq_zero",
        "mapq_missing",
        "overlapping_alternatives",
        "same_bounds_distinct_indels",
        "two_ordered_loci",
    ] {
        let td = tempfile::tempdir().unwrap();
        let (engine, trace) = source(td.path(), true, PAYLOAD);
        let prepared = engine
            .prepare_primer_payload_alignment_queries(
                prepare(trace.clone(), td.path()),
                &mut |_| true,
            )
            .unwrap();
        let mut req = request(td.path(), trace, &prepared);
        let query = &prepared.queries[0];
        let n = query.query_length;
        let line = |target: &str, a: usize, z: usize, mapq: usize| {
            format!(
                "{}\t{n}\t{a}\t{z}\t+\t{target}\t{}\t10\t{}\t{}\t{}\t{mapq}\tcg:Z:{}M\n",
                query.query_id,
                PAYLOAD.len() + 10,
                10 + z - a,
                z - a,
                z - a,
                z - a
            )
        };
        let rows = match case {
            "expected" => line("toy_expected", 0, n, 60),
            "foreign" => line("toy_foreign", 0, n, 60),
            "mapq_zero" => line("toy_expected", 0, n, 0),
            "mapq_missing" => line("toy_expected", 0, n, 255),
            "overlapping_alternatives" => {
                line("toy_expected", 0, n, 60) + &line("toy_foreign", 0, n, 60)
            }
            "same_bounds_distinct_indels" => {
                let base = format!(
                    "{}\t{n}\t0\t{n}\t+\ttoy_expected\t{}\t9\t{}\t{n}\t{}\t60\t",
                    query.query_id,
                    n + 10,
                    n + 10,
                    n + 1
                );
                format!("{base}cg:Z:1D{n}M\n{base}cg:Z:{}M1D{}M\n", n / 2, n - n / 2)
            }
            _ => line("toy_expected", 0, n / 2, 60) + &line("toy_foreign", n / 2, n, 60),
        };
        fs::write(&req.batches[0].paf.path, rows).unwrap();
        req.batches[0].paf = file(Path::new(&req.batches[0].paf.path));
        engine
            .summarize_primer_payload_alignments(req.clone(), &mut |_| true)
            .unwrap();
        let inspected = GentleEngine::inspect_primer_payload(
            &Path::new(&req.output_dir)
                .join("report.json")
                .display()
                .to_string(),
            None,
            1,
            true,
        )
        .unwrap();
        let row = &inspected.rows[0];
        let expected = match case {
            "expected" => "expected_dominant",
            "foreign" => "off_target_dominant",
            "two_ordered_loci" => "mixed",
            _ => "ambiguous",
        };
        assert_eq!(row["classification"], expected, "{case}");
        assert_eq!(
            row["cross_locus_payload_evidence"],
            case == "two_ordered_loci",
            "{case}"
        );
        assert!(row["expected_target_bases"].as_u64().unwrap() <= n as u64);
    }
}

#[test]
fn primer_payload_profiles_deduplicate_anchors_and_clusters_cannot_transitively_bridge() {
    let td = tempfile::tempdir().unwrap();
    let (engine, trace_binding) = source(td.path(), true, PAYLOAD);
    let trace = *GentleEngine::inspect_primer_trace(
        &trace_binding.report_path,
        Some(&trace_binding),
        None,
        0,
    )
    .unwrap()
    .report;
    let prepared = engine
        .prepare_primer_payload_alignment_queries(
            prepare(trace_binding.clone(), td.path()),
            &mut |_| true,
        )
        .unwrap();
    let req = request(td.path(), trace_binding, &prepared);
    engine
        .summarize_primer_payload_alignments(req.clone(), &mut |_| true)
        .unwrap();
    let inspected = GentleEngine::inspect_primer_payload(
        &Path::new(&req.output_dir)
            .join("report.json")
            .display()
            .to_string(),
        None,
        1,
        false,
    )
    .unwrap();
    let original: PrimerPayloadReadEvidence =
        serde_json::from_value(inspected.rows[0].clone()).unwrap();
    let mut rows = (0..3)
        .map(|i| {
            let mut r = original.clone();
            r.record_index = i;
            r.primer_edge_alternatives[0].coordinate = 10 + i as usize * 2;
            r
        })
        .collect::<Vec<_>>();
    rows.push(rows[0].clone());
    let (profiles, clusters) =
        evidence::aggregate(&rows, &prepared, &trace, &["synthetic_reference".into()], 2);
    assert_eq!(profiles[0].supporting_reads, 3);
    assert_eq!(profiles[0].eligible_reads, 3);
    assert_eq!(profiles[0].anchor_observations, 4);
    assert_eq!(clusters.len(), 2);
    assert_eq!(
        (clusters[0].minimum_edge, clusters[0].maximum_edge),
        (10, 12)
    );
    assert_eq!(clusters[1].minimum_edge, 14);
    rows[0].primer_ids = vec!["a_different_primer_same_gene".into()];
    let (profiles, _) =
        evidence::aggregate(&rows, &prepared, &trace, &["synthetic_reference".into()], 2);
    assert_eq!(profiles.len(), 2);
}

#[test]
fn primer_payload_paf_cg_cs_reverse_indels_and_contradictions() {
    let query = q();
    let b = limits();
    let row = "q_synthetic\t50\t2\t17\t-\tcontig\t1000\t50\t66\t14\t17\t255\tcg:Z:5M1I4M2D5M\tcs:Z::5+a:4-cc:5\tNM:i:3\ttp:A:S";
    let (a, length) = paf::row(row, &query, true, &b).unwrap().unwrap();
    assert_eq!(length, 1000);
    assert_eq!(a.segment.mapq, None);
    assert_eq!(a.segment.reported_role.as_deref(), Some("S"));
    assert_eq!(a.blocks[0].query, (2, 7));
    assert_eq!(a.blocks.last().unwrap().query, (12, 17));
    for changed in [
        row.replace("1I", "2I"),
        row.replace("NM:i:3", "NM:i:2"),
        row.replace("cs:Z:", "cs:i:"),
        row.replace("\t2\t17\t", "\t2\t51\t"),
        format!("{row}\tNM:i:3"),
    ] {
        assert!(paf::row(&changed, &query, true, &b).is_err(), "{changed}");
    }
    assert_eq!(
        paf::difference("=acg*at+tt-cc~gt12ag:2", &b).unwrap().len(),
        6
    );
    for cs in ["", ":0", "*a", "~gtxag", "+", ":2?"] {
        assert!(paf::difference(cs, &b).is_err(), "{cs}");
    }
}

#[test]
fn primer_payload_masked_union_never_counts_a_primer_or_duplicates_a_base() {
    assert_eq!(
        paf::subtract_masks(vec![(0, 20), (10, 30), (0, 20)], &[(8, 12), (20, 25)]),
        vec![(0, 8), (12, 20), (25, 30)]
    );
    let b = paf::Block {
        query: (3, 13),
        target: (20, 30),
    };
    assert_eq!(reference::project(&b, (22, 25), '+'), Some((5, 8)));
    assert_eq!(reference::project(&b, (22, 25), '-'), Some((8, 11)));
}

#[test]
fn primer_payload_coarse_alignment_cannot_claim_exact_geometry() {
    let text = "q_synthetic\t50\t0\t50\t+\ttarget\t100\t10\t60\t50\t50\t0";
    assert!(paf::row(text, &q(), true, &limits()).is_err());
    let (row, _) = paf::row(text, &q(), false, &limits()).unwrap().unwrap();
    assert!(!row.segment.alignment_detail_verified);
    assert!(row.blocks.is_empty());
    assert_eq!(row.segment.mapq, Some(0));
}

#[test]
fn primer_payload_genomic_splice_citations_collapse_without_chromosome_aliases() {
    let td = tempfile::tempdir().unwrap();
    let (engine, trace) = source(td.path(), true, PAYLOAD);
    let prepared = engine
        .prepare_primer_payload_alignment_queries(prepare(trace.clone(), td.path()), &mut |_| true)
        .unwrap();
    let mut req = request(td.path(), trace, &prepared);
    let reference = &mut req.references[0];
    reference.scope = PrimerPayloadReferenceScope::Genome;
    let annotation = PrimerPayloadAnnotation {
        schema: PRIMER_PAYLOAD_ANNOTATION_SCHEMA.into(),
        reference_id: reference.reference_id.clone(),
        assembly: reference.assembly.clone(),
        release: reference.release.clone(),
        reference_fasta_sha256: reference.fasta.sha256.clone(),
        coordinate_convention: "zero_based_half_open".into(),
        transcripts: ["isoform_a", "isoform_b"]
            .into_iter()
            .map(|id| PrimerPayloadTranscript {
                transcript_id: id.into(),
                gene_ids: vec!["shared_gene".into()],
                target_name: "exact_contig".into(),
                target_length: 1000,
                genomic_contig: Some("exact_contig".into()),
                genomic_strand: Some('-'),
                genomic_exons: vec![(10, 20), (30, 40)],
            })
            .collect(),
    };
    fs::write(
        &reference.annotation.path,
        serde_json::to_vec(&annotation).unwrap(),
    )
    .unwrap();
    reference.annotation = file(Path::new(&reference.annotation.path));
    let cat = reference::catalog(reference, &limits()).unwrap();
    let query = q();
    let (mut a, _) = paf::row(
        "q_synthetic\t50\t0\t20\t-\texact_contig\t1000\t10\t40\t20\t20\t60\tcg:Z:10M10N10M",
        &query,
        true,
        &limits(),
    )
    .unwrap()
    .unwrap();
    reference::annotate(&mut a, reference, &cat, &mut 0, &limits()).unwrap();
    assert_eq!(a.segment.splice_edges.len(), 1);
    assert_eq!(a.segment.splice_edges[0].transcript_ids.len(), 2);
    assert_eq!(
        (
            a.segment.splice_edges[0].donor,
            a.segment.splice_edges[0].acceptor
        ),
        (30, 20)
    );
}

#[test]
fn primer_payload_alignment_prestart_cancel_and_small_budgets_fail_closed() {
    let td = tempfile::tempdir().unwrap();
    let (engine, trace) = source(td.path(), true, PAYLOAD);
    let prepared = engine
        .prepare_primer_payload_alignment_queries(prepare(trace.clone(), td.path()), &mut |_| true)
        .unwrap();
    let mut req = request(td.path(), trace, &prepared);
    req.cancel_path = Some(td.path().join("cancel").display().to_string());
    GentleEngine::cancel_primer_payload(&PrimerPayloadControlRequest::Alignments(req.clone()))
        .unwrap();
    let report = engine
        .summarize_primer_payload_alignments(req.clone(), &mut |_| true)
        .unwrap();
    assert_eq!(report.status, "cancelled");
    assert!(report.batches_completed.is_empty());
    assert_eq!(report.batches_not_completed.len(), 1);
    req.output_dir = td.path().join("budget refused").display().to_string();
    req.cancel_path = None;
    req.budgets.max_output_bytes = 10;
    assert!(
        engine
            .summarize_primer_payload_alignments(req.clone(), &mut |_| true)
            .is_err()
    );
    assert!(!Path::new(&req.output_dir).exists());
}
