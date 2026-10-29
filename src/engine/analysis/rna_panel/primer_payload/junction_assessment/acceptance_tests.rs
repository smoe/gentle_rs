//! Generated CD274-like geometry only, not patient/cell-line evidence or an assembly snapshot.
use super::*;

pub(super) fn cd274_like(
    root: &Path,
    incomplete: bool,
) -> (GentleEngine, AssessObservedSpliceJunctionsRequest) {
    use super::super::super::tests::{PAYLOAD, bound_file, limits, prepare, source};
    let (e, trace) = source(root, true, PAYLOAD);
    let old: RnaPrimerTraceReport =
        serde_json::from_slice(&fs::read(&trace.report_path).unwrap()).unwrap();
    let mut trace_request = old.request;
    trace_request.files.clear();
    for (i, total) in [4, 2, 4, 1].into_iter().enumerate() {
        let path = root.join(format!("library{i}.fastq"));
        let sequence = format!("AGTCCGATCA{PAYLOAD}");
        let bytes = (0..total)
            .map(|n| {
                format!(
                    "@library{i}_read{n}\n{sequence}\n+\n{}\n",
                    "I".repeat(sequence.len())
                )
            })
            .collect::<String>();
        fs::write(&path, &bytes).unwrap();
        trace_request.files.push(serde_json::from_value(json!({"sample_id":format!("synthetic_context{i}"),"barcode":format!("barcode{i}"),"priming":"synthetic","path":path,"declared_complete":true,"expected_sha256":sha256_prefixed_bytes(bytes.as_bytes())})).unwrap());
    }
    trace_request.output_dir = root.join("four library trace").display().to_string();
    let trace_report = e
        .trace_primer_reads(trace_request.clone(), &mut |_| true)
        .unwrap();
    let binding = GentleEngine::primer_trace_result_binding(&trace_request, &trace_report).unwrap();
    let mut prep = prepare(binding.clone(), root);
    prep.file_ids = (0..4).map(|i| format!("file_{i:04}")).collect();
    let queries = e
        .prepare_primer_payload_alignment_queries(prep, &mut |_| true)
        .unwrap();
    let mut phase2 = super::super::super::alignment_tests::request(root, binding, &queries);
    let (start, end, len) = (5_466_826usize, 5_468_192usize, 5_468_300usize);
    // Reproducible LCG: seed 0x4a31, x = x*6364136223846793005+1 (u64),
    // bases selected by the high two bits. Only the declared motif is planted.
    let mut seed = 0x4a31u64;
    let mut bases = (0..len)
        .map(|_| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(seed >> 62) as usize]
        })
        .collect::<Vec<_>>();
    bases[start..start + 2].copy_from_slice(b"GT");
    bases[end - 2..end].copy_from_slice(b"AG");
    let mut fasta = b">synthetic_cd274_contig\n".to_vec();
    for line in bases.chunks(80) {
        fasta.extend_from_slice(line);
        fasta.push(b'\n');
    }
    let reference = &mut phase2.references[0];
    fs::write(&reference.fasta.path, fasta).unwrap();
    reference.fasta = bound_file(Path::new(&reference.fasta.path));
    reference.scope = PrimerPayloadReferenceScope::Genome;
    reference.assembly = "synthetic_cd274_like_v1".into();
    let mut ann = PrimerPayloadAnnotation {
        schema: PRIMER_PAYLOAD_ANNOTATION_SCHEMA.into(),
        reference_id: reference.reference_id.clone(),
        assembly: reference.assembly.clone(),
        release: "synthetic_ensembl_like".into(),
        reference_fasta_sha256: reference.fasta.sha256.clone(),
        coordinate_convention: "zero_based_half_open".into(),
        transcripts: vec![PrimerPayloadTranscript {
            transcript_id: "synthetic_ensembl_tx".into(),
            gene_ids: vec!["synthetic_cd274_gene".into()],
            target_name: "synthetic_cd274_contig".into(),
            target_length: len,
            genomic_contig: Some("synthetic_cd274_contig".into()),
            genomic_strand: Some('+'),
            genomic_exons: vec![(start - 30, start), (start + 600, start + 650)],
        }],
    };
    reference.release = ann.release.clone();
    fs::write(
        &reference.annotation.path,
        serde_json::to_vec(&ann).unwrap(),
    )
    .unwrap();
    reference.annotation = bound_file(Path::new(&reference.annotation.path));
    let primary = reference.clone();
    phase2.expected_targets[0].gene_ids = vec!["synthetic_cd274_gene".into()];
    phase2.expected_targets[0].transcript_ids = vec!["synthetic_ensembl_tx".into()];
    let mut seen = BTreeMap::<String, usize>::new();
    let paf=queries.queries.iter().map(|q| {
        let f=q.file_id.trim_start_matches("file_").parse::<usize>().unwrap();
        let n=seen.entry(q.file_id.clone()).or_default(); let support=*n < [3,1,3,0][f]; *n+=1;
        let accept=if support {end} else {start+600}; let half=q.query_length/2;
        format!("{}\t{}\t0\t{}\t+\tsynthetic_cd274_contig\t{len}\t{}\t{}\t{}\t{}\t60\tcg:Z:{half}M{}N{}M\ttp:A:P\n",q.query_id,q.query_length,q.query_length,start-half,accept+q.query_length-half,q.query_length,q.query_length,accept-start,q.query_length-half)
    }).collect::<String>();
    fs::write(&phase2.batches[0].paf.path, paf).unwrap();
    phase2.batches[0].paf = bound_file(Path::new(&phase2.batches[0].paf.path));
    phase2.batches[0].competitor_search_limitations = if incomplete {
        "synthetic incomplete search"
    } else {
        "synthetic exhaustive intention; frozen producer still requires nonempty limitations"
    }
    .into();
    let p = e
        .summarize_primer_payload_alignments(phase2, &mut |_| true)
        .unwrap();
    let pb = GentleEngine::primer_payload_alignment_result_binding(&p.request, &p).unwrap();
    let mut second = primary.clone();
    second.reference_id = "synthetic_refseq_reference".into();
    second.release = "synthetic_refseq_like".into();
    ann.reference_id = second.reference_id.clone();
    ann.release = second.release.clone();
    ann.transcripts[0].transcript_id = "synthetic_refseq_tx".into();
    ann.transcripts[0].genomic_exons[1] = (start + 601, start + 651);
    let path = root.join("refseq-like.json");
    fs::write(&path, serde_json::to_vec(&ann).unwrap()).unwrap();
    second.annotation = bound_file(&path);
    let s1req = SummarizeObservedTranscriptStructuresRequest {
        schema: "gentle.observed_transcript_structures_request.v1".into(),
        source_report: bound_file(Path::new(&pb.report_path)),
        source: pb,
        batch_id: "synthetic_batch".into(),
        reference_gene: s1::ReferenceGene {
            source_id: "ensembl_like".into(),
            gene_id: "synthetic_cd274_gene".into(),
        },
        comparison_genes: vec![s1::ReferenceGene {
            source_id: "refseq_like".into(),
            gene_id: "synthetic_cd274_gene".into(),
        }],
        annotation_sources: vec![
            ObservedTranscriptAnnotation {
                source_id: "ensembl_like".into(),
                provider: Some("synthetic".into()),
                reference: primary.clone(),
            },
            ObservedTranscriptAnnotation {
                source_id: "refseq_like".into(),
                provider: Some("synthetic".into()),
                reference: second,
            },
        ],
        required_flank_bases: 2,
        budgets: limits(),
        kernel_budgets: s1::Budgets {
            max_queries: 100,
            max_records: 100,
            max_segments_per_query: 100,
            max_operations: 10000,
            max_skips_per_chain: 100,
            max_chain_families: 100,
            max_alternatives_per_record: 100,
            max_annotation_rows: 100,
            max_model_comparisons_per_family: 100,
            max_total_work: 1_000_000,
        },
        output_dir: root.join("saved structures").display().to_string(),
        cancel_path: None,
    };
    let s = e
        .summarize_observed_transcript_structures(s1req, &mut |_| true)
        .unwrap();
    let sb = GentleEngine::observed_transcript_structures_result_binding(&s.request, &s).unwrap();
    let manifest = ObservedSampleManifest {
        schema: OBSERVED_SAMPLE_MANIFEST_SCHEMA.into(),
        rows: (0..4)
            .map(|i| ObservedSampleBinding {
                context: ObservedSampleContext {
                    sample_id: format!("synthetic_context{i}"),
                    barcode: format!("barcode{i}"),
                    priming: "synthetic".into(),
                    file_id: format!("file_{i:04}"),
                },
                library_id: format!("L{i}"),
                cdna_prep_id: format!("prep{}_sample{}", i / 2, i % 2),
                biological_sample_id: format!("sample{}", i % 2),
                cell_line: None,
            })
            .collect(),
    };
    let path = root.join("sample manifest.json");
    fs::write(&path, serde_json::to_vec(&manifest).unwrap()).unwrap();
    let req = AssessObservedSpliceJunctionsRequest {
        schema: OBSERVED_JUNCTION_REQUEST_SCHEMA.into(),
        source: ObservedJunctionSource::S1Report {
            report: bound_file(Path::new(&sb.report_path)),
            result: sb,
        },
        reference_sequence: Some(ObservedJunctionReference {
            file: primary.fasta,
            format: ObservedJunctionReferenceFormat::FastaPlainV1,
            authority: ObservedJunctionReferenceAuthority::FullAssembly,
        }),
        sample_manifest: Some(bound_file(&path)),
        near_window_bp: 5,
        recurrence_policy: Some(ObservedRecurrencePolicy {
            unit: ObservedRecurrenceUnit::Library,
            min_units: 2,
            certainty: ObservedRecurrenceCertainty::Definite,
            require_qualified_flanks: true,
        }),
        budgets: ObservedJunctionBudgets {
            max_input_bytes: 128_000_000,
            max_reference_bytes: 16_000_000,
            max_manifest_bytes: 1_000_000,
            max_distinct_junctions: 100,
            max_annotated_junctions_compared: 10000,
            max_near_hits_per_junction: 1000,
            max_total_work: 1_000_000_000,
            max_output_bytes: 128_000_000,
        },
        output_dir: root.join("junctions").display().to_string(),
        cancel_path: None,
    };
    (e, req)
}

#[test]
fn junction_assessment_synthetic_cd274_like_3130_not_real_rna_evidence() {
    for incomplete in [false, true] {
        let td = tempfile::tempdir().unwrap();
        let (e, r) = cd274_like(td.path(), incomplete);
        let request = r.clone();
        let report = e
            .assess_observed_splice_junctions(r, &mut |_| true)
            .unwrap();
        let data = data(&report);
        let j = data
            .junctions
            .iter()
            .find(|j| j.skipped_interval == (5_466_826, 5_468_192))
            .unwrap();
        assert_eq!(j.donor_exon_last_base_1based, Some(5_466_826));
        assert_eq!(j.acceptor_exon_first_base_1based, Some(5_468_193));
        assert_eq!(j.intron_first_base_1based, 5_466_827);
        assert_eq!(j.intron_last_base_1based, 5_468_192);
        assert_eq!(j.motif.motif_class, ObservedMotifClass::GtAg);
        assert_eq!(j.motif.reference_authority, "verified");
        assert_eq!(j.sources.len(), 2);
        assert!(j.sources.iter().all(|s| s.relation != "exact"));
        assert!(
            j.sources
                .iter()
                .any(|s| s.events.contains(&"alternative_terminal_acceptor".into()))
        );
        for (i, n) in [3, 1, 3, 0].into_iter().enumerate() {
            let support = j
                .grouped_support
                .iter()
                .find(|g| g.dimension == "library" && g.id == format!("L{i}"));
            let observed = support
                .map(|g| g.counts.possible_only.qualified_only)
                .unwrap_or(0);
            assert_eq!(observed, n, "L{i} incomplete={incomplete}");
            assert!(
                data.libraries
                    .iter()
                    .find(|l| l.library_id == format!("L{i}"))
                    .unwrap()
                    .gene_coherent_records
                    > n
            );
        }
        assert_eq!(j.support.definite.qualified_only, 0);
        assert_eq!(j.support.possible_only.qualified_only, 7);
        assert_eq!(j.recurrence.as_ref().unwrap().status, "not_assessable");
        if !incomplete {
            // Pure synthetic contract test only: the frozen Phase 2 producer cannot
            // emit search-complete S1 today. Never republish/promote its saved evidence.
            let l = loaded(&request);
            let mut summary = l.data.analysis.unwrap();
            for record in &mut summary.records {
                record.search_declared_complete = true;
                record.all_query_views_assessed = true;
                record.definite_chains = record.possible_chains.clone();
            }
            let mut meter = Meter {
                b: &request.budgets,
                input: 0,
                work: 0,
                comparisons: 0,
            };
            let mut progress = |_| true;
            let mut ctl = Control::new(None, "synthetic".into(), 1, &mut progress).unwrap();
            let manifest = source::manifest(&request, &mut meter, &mut ctl).unwrap();
            let index = support::record_index(&summary, &l.contexts, &mut meter, &mut ctl).unwrap();
            let records = &index[&(j.contig.clone(), j.skipped_interval)];
            assert_eq!(support::counts(records).definite.qualified_only, 7);
            assert_eq!(
                support::recurrence(
                    request.recurrence_policy.as_ref().unwrap(),
                    records,
                    &manifest,
                    false
                )
                .status,
                "recurrent_candidate"
            );
            let groups = support::groups(records, &manifest);
            for (i, n) in [3, 1, 3, 0].into_iter().enumerate() {
                assert_eq!(
                    groups
                        .iter()
                        .find(|g| g.dimension == "library" && g.id == format!("L{i}"))
                        .map(|g| g.counts.definite.qualified_only)
                        .unwrap_or(0),
                    n
                );
            }
        }
        assert!(serde_json::to_value(j).unwrap().get("protein").is_none());
    }
}
