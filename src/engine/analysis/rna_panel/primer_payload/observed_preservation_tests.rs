//! S2.0a independent v1 publication replay for the parent's inline synthetic S1 fixture.
//!
//! Run the `s2a_` library-test filter. This freezes field/row recipes, not a second
//! interpreter: expected geometry uses the separately fingerprinted S0 kernel.
//! Actual temporary paths and every dependent digest are retained, never scrubbed
//! or normalized, so the same oracle works for Unix and Windows path spellings.

use super::*;

#[derive(Serialize)]
struct RequestV1<'a> {
    schema: &'a str,
    source: &'a PrimerPayloadResultBinding,
    source_report: &'a PrimerPayloadFile,
    batch_id: &'a str,
    reference_gene: &'a kernel::ReferenceGene,
    comparison_genes: &'a [kernel::ReferenceGene],
    annotation_sources: &'a [ObservedTranscriptAnnotation],
    required_flank_bases: usize,
    budgets: &'a PrimerPayloadBudgets,
    kernel_budgets: &'a kernel::Budgets,
    output_dir: &'a str,
    cancel_path: &'a Option<String>,
}

#[derive(Serialize)]
struct DataV1<'a> {
    schema: &'static str,
    request_sha256: &'a str,
    analysis: Option<&'a kernel::Summary>,
    unavailable: &'a [ObservedTranscriptUnavailable],
}

#[derive(Serialize)]
struct CountsV1 {
    records: usize,
    observations: usize,
    chain_families: usize,
    unavailable_views: usize,
}

#[derive(Serialize)]
struct ReportV1<'a> {
    schema: &'static str,
    content_sha256: String,
    request_sha256: &'a str,
    request: &'a RequestV1<'a>,
    status: &'static str,
    control_error: Option<String>,
    source_status: &'a str,
    source_trace_status: &'a str,
    source_query_status: &'a str,
    cohort: &'a PrimerPayloadCohort,
    source_batches_completed: &'a [String],
    source_batches_not_completed: &'a [String],
    source_inputs_missing: usize,
    source_query_units_missing: usize,
    source_excluded_anchor_count: usize,
    counts: Option<CountsV1>,
    artifacts: Vec<RnaReadPanelArtifact>,
    nonclaims: [&'static str; 4],
}

fn canonical_sha<T: Serialize>(value: &T) -> String {
    sha256_prefixed_bytes(&serde_json::to_vec(&serde_json::to_value(value).unwrap()).unwrap())
}

fn json_cell_table<T: Serialize>(rows: &[T]) -> Vec<u8> {
    // Frozen CSV v1 encoding, deliberately independent of the production writer.
    let mut bytes = b"row_json\n".to_vec();
    for row in rows {
        let cell = serde_json::to_string(row).unwrap();
        bytes.extend_from_slice(b"\"");
        bytes.extend_from_slice(cell.replace('"', "\"\"").as_bytes());
        bytes.extend_from_slice(b"\"\n");
    }
    bytes
}

fn expected_analysis(
    r: &SummarizeObservedTranscriptStructuresRequest,
    source: &PrimerPayloadAlignmentReport,
    q: &PrimerPayloadQuery,
) -> kernel::Summary {
    let trace = GentleEngine::inspect_primer_trace(
        &source.request.trace.report_path,
        Some(&source.request.trace),
        None,
        0,
    )
    .unwrap();
    let file = trace
        .report
        .files_completed
        .iter()
        .find(|f| f.file_id == q.file_id)
        .unwrap();
    let mut rows = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(Path::new(&source.request.output_dir).join("artifacts/read_anchor_evidence.tsv"))
        .unwrap();
    let row = rows.records().next().unwrap().unwrap();
    let evidence: PrimerPayloadReadEvidence = serde_json::from_str(&row[0]).unwrap();
    let reference = &r.annotation_sources[0].reference;
    let annotation: PrimerPayloadAnnotation = read_json(&reference.annotation, &r.budgets).unwrap();
    let lengths = annotation
        .transcripts
        .iter()
        .map(|t| (t.target_name.clone(), t.target_length))
        .collect::<BTreeMap<_, _>>();
    assert!(q.excluded_query_intervals.is_empty());
    let input = kernel::Input {
        reference: kernel::ReferenceBinding {
            assembly: reference.assembly.clone(),
            fasta_sha256: reference.fasta.sha256.clone(),
        },
        reference_gene: r.reference_gene.clone(),
        comparison_genes: r.comparison_genes.clone(),
        annotation_sources: vec![kernel::AnnotationSource {
            source_id: r.annotation_sources[0].source_id.clone(),
            provider: r.annotation_sources[0].provider.clone(),
            annotation_sha256: reference.annotation.sha256.clone(),
            annotation,
        }],
        queries: vec![kernel::QueryEvidence {
            record: kernel::ReadRecord {
                source_sha256: file.observed_sha256.clone(),
                record_index: q.record_index,
                read_id: q.read_id.clone(),
            },
            alignment_sha256: source.request.batches[0].paf.sha256.clone(),
            query_id: q.query_id.clone(),
            query_length: q.query_length,
            scope: PrimerPayloadReferenceScope::Genome,
            target_lengths: evidence
                .segments
                .iter()
                .map(|s| (s.target_name.clone(), lengths[&s.target_name]))
                .collect(),
            masks: Some(vec![]),
            ends: [
                vec![kernel::EndConstraint::PrimerBounded],
                vec![kernel::EndConstraint::ReadEndUnverified],
            ],
            // The synthetic producer explicitly disclaims search-completeness proof.
            search_complete: false,
            segments: evidence.segments,
        }],
        unassessed_records: vec![],
        required_flank_bases: r.required_flank_bases,
        budgets: r.kernel_budgets.clone(),
    };
    kernel::summarize(&input).unwrap()
}

#[test]
fn observed_transcripts_s2a_saved_v1_report_data_and_all_tables_keep_exact_bytes() {
    for completed in [true, false] {
        let td = tempfile::tempdir().unwrap();
        let (engine, r) = fixture_mode(td.path(), completed);
        let source: PrimerPayloadAlignmentReport = read_json(&r.source_report, &limits()).unwrap();
        let queries: PrimerPayloadQueries = read_json(&source.request.queries, &limits()).unwrap();
        assert_eq!(queries.queries.len(), 1);
        let q = &queries.queries[0];
        let request = RequestV1 {
            schema: "gentle.observed_transcript_structures_request.v1",
            source: &r.source,
            source_report: &r.source_report,
            batch_id: "synthetic_batch",
            reference_gene: &r.reference_gene,
            comparison_genes: &r.comparison_genes,
            annotation_sources: &r.annotation_sources,
            required_flank_bases: 2,
            budgets: &r.budgets,
            kernel_budgets: &r.kernel_budgets,
            output_dir: &r.output_dir,
            cancel_path: &None,
        };
        assert_eq!(
            serde_json::to_vec(&r).unwrap(),
            serde_json::to_vec(&request).unwrap()
        );
        let request_json = serde_json::to_value(&r).unwrap();
        for field in [
            "schema",
            "source",
            "source_report",
            "batch_id",
            "reference_gene",
            "comparison_genes",
            "annotation_sources",
            "required_flank_bases",
            "budgets",
            "kernel_budgets",
            "output_dir",
        ] {
            let mut missing = request_json.clone();
            missing.as_object_mut().unwrap().remove(field).unwrap();
            assert!(
                serde_json::from_value::<SummarizeObservedTranscriptStructuresRequest>(missing)
                    .is_err(),
                "request default added for {field}"
            );
        }
        let mut unknown = request_json;
        unknown["not_a_v1_request_field"] = json!(true);
        assert!(
            serde_json::from_value::<SummarizeObservedTranscriptStructuresRequest>(unknown)
                .is_err()
        );
        let identity = canonical_sha(&request);
        let analysis = completed.then(|| expected_analysis(&r, &source, q));
        let unavailable = if completed {
            vec![]
        } else {
            vec![ObservedTranscriptUnavailable {
                file_id: q.file_id.clone(),
                record_index: q.record_index,
                read_id: q.read_id.clone(),
                query_id: Some(q.query_id.clone()),
                reason: "batch_not_completed".into(),
            }]
        };
        let expected_data = DataV1 {
            schema: "gentle.observed_transcript_structures_data.v1",
            request_sha256: &identity,
            analysis: analysis.as_ref(),
            unavailable: &unavailable,
        };
        let mut expected_files = vec![(
            "structures.json",
            serde_json::to_vec(&expected_data).unwrap(),
        )];
        let (mut reads, mut families, mut comparisons, mut junctions) =
            (vec![], vec![], vec![], vec![]);
        if let Some(s) = &analysis {
            assert_eq!(
                (s.records.len(), s.observations.len(), s.families.len()),
                (1, 1, 1)
            );
            let record = &s.records[0];
            let observation = &s.observations[0];
            reads.push(json!({"support":record,"queries":[q],"observations":[observation]}));
            families.push(json!({"family":s.families[0],"records":[record],"support_by_context":[{
                "context":{"sample_id":q.sample_id,"barcode":q.barcode,"priming":q.priming,"primer_id":"invented_primer","file_id":q.file_id},
                "source_records":1,"records":[record.record]
            }],"groups_are_nonadditive":true}));
            for comparison in &observation.comparisons {
                comparisons.push(json!({"observation_id":observation.observation_id,"record":record.record,"comparison":comparison}));
            }
            assert_eq!(observation.skips.len(), 1);
            let skip = observation.skips[0];
            let comparison = &observation.comparisons[0];
            junctions.push(json!({"contig":observation.contig,"skip":skip,"record":record.record,
                "observation_id":observation.observation_id,"flank_qualification":observation.flank_qualification[0],
                "on_best_path":true,"source_comparisons":[{"source_id":comparison.source_id,
                "gene_id":comparison.gene_id,"resolution":comparison.resolution,"junction":comparison.junctions.as_ref().unwrap()[0]}],
                "note":"Physical skip observation, not independent molecule or natural RNA strand evidence"}));
        }
        for (name, rows) in [
            ("read_structures.tsv", reads),
            ("structure_catalog.tsv", families),
            ("reference_compatibility.tsv", comparisons),
            ("junction_catalog.tsv", junctions),
        ] {
            expected_files.push((name, json_cell_table(&rows)));
        }
        expected_files.push(("unavailable.tsv", json_cell_table(&unavailable)));
        let artifacts = expected_files
            .iter()
            .map(|(name, bytes)| RnaReadPanelArtifact {
                path: (*name).into(),
                sha256: sha256_prefixed_bytes(bytes),
                bytes: bytes.len() as u64,
                locus_id: None,
            })
            .collect();
        let mut expected_report = ReportV1 {
            schema: "gentle.observed_transcript_structures_report.v1",
            content_sha256: String::new(),
            request_sha256: &identity,
            request: &request,
            status: if completed {
                "completed"
            } else {
                "unavailable"
            },
            control_error: None,
            source_status: &source.status,
            source_trace_status: &source.source_trace_status,
            source_query_status: &source.source_query_status,
            cohort: &source.cohort,
            source_batches_completed: &source.batches_completed,
            source_batches_not_completed: &source.batches_not_completed,
            source_inputs_missing: source.source_inputs_missing,
            source_query_units_missing: source.source_query_units_missing,
            source_excluded_anchor_count: source.source_excluded_anchor_count,
            counts: completed.then_some(CountsV1 {
                records: 1,
                observations: 1,
                chain_families: 1,
                unavailable_views: 0,
            }),
            artifacts,
            nonclaims: [
                "Internal content consistency is not authenticated provenance or an independently repeated alignment search.",
                "Support counts source records, not independent molecules, expression, exclusive isoforms or proteins.",
                "Novel means relative to the explicitly bound annotation snapshots; query orientation does not establish RNA strand.",
                "Missing views and unknown target lengths remain unavailable; payload ends are not natural RNA ends. Original FASTQ, PAF and reference FASTA are not reopened.",
            ],
        };
        expected_report.content_sha256 = canonical_sha(&expected_report);
        let report = engine
            .summarize_observed_transcript_structures(r.clone(), &mut |_| true)
            .unwrap();
        assert_eq!(
            serde_json::to_vec(&report).unwrap(),
            serde_json::to_vec(&expected_report).unwrap()
        );
        let report_json = serde_json::to_value(&report).unwrap();
        for field in [
            "schema",
            "content_sha256",
            "request_sha256",
            "request",
            "status",
            "source_status",
            "source_trace_status",
            "source_query_status",
            "cohort",
            "source_batches_completed",
            "source_batches_not_completed",
            "source_inputs_missing",
            "source_query_units_missing",
            "source_excluded_anchor_count",
            "artifacts",
            "nonclaims",
        ] {
            let mut missing = report_json.clone();
            missing.as_object_mut().unwrap().remove(field).unwrap();
            assert!(
                serde_json::from_value::<ObservedTranscriptStructuresReport>(missing).is_err(),
                "report default added for {field}"
            );
        }
        let mut unknown = report_json;
        unknown["not_a_v1_report_field"] = json!(true);
        assert!(serde_json::from_value::<ObservedTranscriptStructuresReport>(unknown).is_err());
        assert_eq!(
            fs::read(Path::new(&r.output_dir).join("report.json")).unwrap(),
            serde_json::to_vec_pretty(&expected_report).unwrap()
        );
        for (name, bytes) in expected_files {
            assert_eq!(
                fs::read(Path::new(&r.output_dir).join("artifacts").join(name)).unwrap(),
                bytes,
                "{name}, completed={completed}"
            );
        }
    }
}

#[test]
fn observed_transcripts_s2a_request_and_source_error_strings_stay_exact() {
    let td = tempfile::tempdir().unwrap();
    let (engine, request) = fixture(td.path());
    for (case, expected) in [
        (0, "observed_transcripts_request"),
        (1, "observed_transcripts_source_binding"),
        (2, "primer_payload_file_hash_or_change"),
    ] {
        let mut r = request.clone();
        match case {
            0 => r.schema = "not-a-supported-schema".into(),
            1 => r.source.report_content_sha256 = format!("sha256:{}", "0".repeat(64)),
            _ => r.source_report.sha256 = format!("sha256:{}", "0".repeat(64)),
        }
        let output = PathBuf::from(&r.output_dir);
        assert_eq!(
            engine
                .summarize_observed_transcript_structures(r, &mut |_| true)
                .unwrap_err()
                .message,
            expected,
            "case {case}"
        );
        assert!(!output.exists());
    }
}
