//! S1: verify saved Phase 2 bytes, interpret one genomic batch, publish separate evidence.

use super::*;
use crate::observed_transcript_structures as kernel;
use serde_json::{Value, json};

const REQUEST: &str = "gentle.observed_transcript_structures_request.v1";
const REPORT: &str = "gentle.observed_transcript_structures_report.v1";
const BINDING: &str = "gentle.observed_transcript_structures_result_binding.v1";
const TABLES: [&str; 5] = [
    "read_structures.tsv",
    "structure_catalog.tsv",
    "reference_compatibility.tsv",
    "junction_catalog.tsv",
    "unavailable.tsv",
];

struct BoundedOutput {
    writer: HashedWriter,
    left: u64,
}
impl Write for BoundedOutput {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() as u64 > self.left {
            return Err(std::io::Error::other("observed_transcripts_output_budget"));
        }
        let n = self.writer.write(bytes)?;
        self.left -= n as u64;
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }
}

fn output_table<T: Serialize>(
    root: &Path,
    name: &str,
    values: &[T],
    left: &mut u64,
    max_line: usize,
) -> Result<RnaReadPanelArtifact, EngineError> {
    let writer = BoundedOutput {
        writer: HashedWriter::create(&root.join(name))?,
        left: *left,
    };
    let mut csv = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);
    csv.write_record(["row_json"])
        .map_err(|e| bad(e.to_string()))?;
    for value in values {
        let text = serde_json::to_string(value).map_err(|e| bad(e.to_string()))?;
        if text.len() > max_line {
            return Err(bad("observed_transcripts_output_row_budget"));
        }
        csv.write_record([text]).map_err(|e| bad(e.to_string()))?;
    }
    csv.flush().map_err(|e| bad(e.to_string()))?;
    let writer = csv.into_inner().map_err(|e| bad(e.to_string()))?;
    *left = writer.left;
    writer.writer.receipt(name.into())
}

fn push(
    rows: &mut Vec<Value>,
    value: Value,
    memory: &mut u64,
    b: &PrimerPayloadBudgets,
) -> Result<(), EngineError> {
    let size = serde_json::to_vec(&value)
        .map_err(|e| bad(e.to_string()))?
        .len();
    if size > b.max_line_bytes {
        return Err(bad("observed_transcripts_output_row_budget"));
    }
    *memory = memory
        .checked_add(size as u64)
        .filter(|n| *n <= b.max_output_bytes)
        .ok_or_else(|| bad("observed_transcripts_output_memory_budget"))?;
    rows.push(value);
    Ok(())
}

fn support_groups(
    records: &[&kernel::RecordSupport],
    metadata: &BTreeMap<kernel::ReadRecord, Vec<PrimerPayloadQuery>>,
) -> Result<Vec<Value>, EngineError> {
    let mut groups = BTreeMap::<String, BTreeSet<kernel::ReadRecord>>::new();
    for record in records {
        for q in metadata.get(&record.record).into_iter().flatten() {
            for primer in &q.primer_ids {
                let context=serde_json::to_string(&json!({"sample_id":q.sample_id,"barcode":q.barcode,"priming":q.priming,"primer_id":primer,"file_id":q.file_id})).map_err(|e|bad(e.to_string()))?;
                groups
                    .entry(context)
                    .or_default()
                    .insert(record.record.clone());
            }
        }
    }
    groups.into_iter().map(|(context,records)|Ok(json!({"context":serde_json::from_str::<Value>(&context).map_err(|e|bad(e.to_string()))?,"source_records":records.len(),"records":records}))).collect()
}

fn admit(r: &SummarizeObservedTranscriptStructuresRequest) -> Result<(), EngineError> {
    budgets(&r.budgets)?;
    r.kernel_budgets.validate(r.required_flank_bases)?;
    if r.schema != REQUEST
        || r.batch_id.is_empty()
        || r.annotation_sources.is_empty()
        || r.annotation_sources.len() > 32
        || r.budgets.max_input_bytes > 128_000_000
        || r.budgets.max_output_bytes > 128_000_000
        || r.budgets.max_queries > 4096
        || r.budgets.max_records > 4096
        || r.source.kind != "alignments"
        || r.reference_gene.source_id.is_empty()
        || r.reference_gene.gene_id.is_empty()
        || Path::new(&r.source_report.path) != Path::new(&r.source.report_path)
        || !is_hash(&r.source_report.sha256)
    {
        return Err(bad("observed_transcripts_request"));
    }
    Ok(())
}

fn report_binding(
    r: &ObservedTranscriptStructuresReport,
) -> Result<ObservedTranscriptStructuresResultBinding, EngineError> {
    admit(&r.request)?;
    let mut basis = r.clone();
    basis.content_sha256.clear();
    let names = r
        .artifacts
        .iter()
        .map(|a| a.path.as_str())
        .collect::<BTreeSet<_>>();
    let expected = TABLES
        .into_iter()
        .chain(["structures.json"])
        .collect::<BTreeSet<_>>();
    if r.schema != REPORT
        || hash(&basis)? != r.content_sha256
        || hash(&r.request)? != r.request_sha256
        || !(binding::outcome(&r.status, &r.control_error)
            || (r.status == "unavailable" && r.control_error.is_none()))
        || (r.status == "completed" && r.counts.is_none())
        || names != expected
        || r.artifacts.len() != expected.len()
        || r.artifacts
            .iter()
            .any(|a| !is_hash(&a.sha256) || a.bytes > r.request.budgets.max_output_bytes)
        || (r.status == "unavailable" && r.counts.is_some())
    {
        return Err(bad("observed_transcripts_report_binding"));
    }
    Ok(ObservedTranscriptStructuresResultBinding {
        schema: BINDING.into(),
        request_sha256: r.request_sha256.clone(),
        report_content_sha256: r.content_sha256.clone(),
        report_path: Path::new(&r.request.output_dir)
            .join("report.json")
            .display()
            .to_string(),
        outcome: r.status.clone(),
    })
}

/// A global actual-input budget, including source/manifest/retained trace dependencies.
fn add_bytes(n: &mut u64, amount: u64, b: &PrimerPayloadBudgets) -> Result<(), EngineError> {
    *n = n
        .checked_add(amount)
        .filter(|n| *n <= b.max_input_bytes)
        .ok_or_else(|| bad("observed_transcripts_input_budget"))?;
    Ok(())
}

fn consume(
    r: &SummarizeObservedTranscriptStructuresRequest,
    source: &PrimerPayloadAlignmentReport,
    ctl: &mut Control<'_>,
) -> Result<
    (
        kernel::Input,
        Vec<ObservedTranscriptUnavailable>,
        Vec<PrimerPayloadQuery>,
    ),
    EngineError,
> {
    let b = &r.budgets;
    let batch = source
        .request
        .batches
        .iter()
        .find(|x| x.batch_id == r.batch_id)
        .ok_or_else(|| bad("observed_transcripts_batch"))?;
    let selected = source
        .request
        .references
        .iter()
        .find(|x| x.reference_id == batch.reference_id)
        .ok_or_else(|| bad("observed_transcripts_reference"))?;
    if selected.scope != PrimerPayloadReferenceScope::Genome {
        return Err(bad("observed_transcripts_genomic_batch_required"));
    }
    let trace = prepare::trace_source(&source.request.trace, b)?;
    let mut additional = 0;
    add_bytes(&mut additional, r.source_report.bytes, b)?;
    for a in &source.artifacts {
        add_bytes(&mut additional, a.bytes, b)?;
    }
    for a in &r.annotation_sources {
        add_bytes(&mut additional, a.reference.annotation.bytes, b)?;
    }
    let prepared = binding::queries(
        &source.request.queries,
        &trace,
        &source.request.trace,
        b,
        additional,
    )?;
    if prepared.content_sha256 != source.request.query_manifest_sha256
        || source.source_trace_status != trace.status
        || source.source_inputs_missing != trace.inputs_not_completed.len()
        || source.source_query_status != prepared.status
        || source.source_query_units_missing != prepared.files_not_completed.len()
        || source.source_excluded_anchor_count != prepared.unassessed.len()
    {
        return Err(bad("observed_transcripts_source_accounting"));
    }
    let mut annotations = Vec::new();
    let mut lengths = BTreeMap::new();
    for a in &r.annotation_sources {
        ctl.check()?;
        if a.reference.scope != PrimerPayloadReferenceScope::Genome
            || a.reference.assembly != selected.assembly
            || a.reference.fasta.sha256 != selected.fasta.sha256
        {
            return Err(bad("observed_transcripts_annotation_reference"));
        }
        let catalog = reference::catalog(&a.reference, b)?;
        for (target, length) in catalog.lengths {
            if lengths.insert(target, length).is_some_and(|n| n != length) {
                return Err(bad("observed_transcripts_target_length_conflict"));
            }
        }
        annotations.push(kernel::AnnotationSource {
            source_id: a.source_id.clone(),
            provider: a.provider.clone(),
            annotation_sha256: a.reference.annotation.sha256.clone(),
            annotation: catalog.annotation,
        });
    }
    let root = Path::new(&r.source_report.path)
        .parent()
        .unwrap_or(Path::new("."))
        .join("artifacts");
    let queries = prepared
        .queries
        .iter()
        .map(|q| (q.query_id.as_str(), q))
        .collect::<BTreeMap<_, _>>();
    let batches = source
        .request
        .batches
        .iter()
        .map(|x| (x.batch_id.as_str(), x))
        .collect::<BTreeMap<_, _>>();
    let mut seen = BTreeSet::new();
    let mut evidence = Vec::new();
    let mut memory = 0u64;
    // Verify every published artifact, including rows from nonselected batches, to EOF.
    for a in &source.artifacts {
        let file = artifact_file(&root, a)?;
        let input = VerifiedReader::open(&file, b.max_input_bytes)?;
        if a.path != "read_anchor_evidence.tsv" {
            let mut input = input;
            std::io::copy(&mut input, &mut std::io::sink()).map_err(|e| bad(e.to_string()))?;
            input.finish()?;
            continue;
        }
        let mut csv = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(input);
        if csv
            .headers()
            .map_err(|e| bad(e.to_string()))?
            .iter()
            .collect::<Vec<_>>()
            != ["evidence_json"]
        {
            return Err(bad("observed_transcripts_evidence_columns"));
        }
        for row in csv.records() {
            ctl.check()?;
            let row = row.map_err(|e| bad(e.to_string()))?;
            if row.len() != 1 || row[0].len() > b.max_line_bytes || seen.len() >= b.max_records {
                return Err(bad("observed_transcripts_evidence_budget"));
            }
            let e: PrimerPayloadReadEvidence =
                serde_json::from_str(&row[0]).map_err(|e| bad(e.to_string()))?;
            let q = queries
                .get(e.query_id.as_str())
                .ok_or_else(|| bad("observed_transcripts_foreign_query"))?;
            let batch = batches
                .get(e.batch_id.as_str())
                .ok_or_else(|| bad("observed_transcripts_foreign_batch"))?;
            let accounting = batch
                .accounting
                .get(&e.query_id)
                .ok_or_else(|| bad("observed_transcripts_accounting"))?;
            if !source.batches_completed.contains(&e.batch_id)
                || !seen.insert((e.batch_id.clone(), e.query_id.clone()))
                || e.reference_id != batch.reference_id
                || e.file_id != q.file_id
                || e.record_index != q.record_index
                || e.read_id != q.read_id
                || e.primer_ids != q.primer_ids
                || e.sample_id != q.sample_id
                || e.barcode != q.barcode
                || encode(&e.priming)? != encode(&q.priming)?
                || e.usable_bases != q.usable_bases
                || (matches!(accounting, PrimerPayloadQueryAccounting::Mapped)
                    == e.segments.is_empty())
            {
                return Err(bad("observed_transcripts_evidence_binding"));
            }
            if e.batch_id == r.batch_id {
                add_bytes(&mut memory, row[0].len() as u64, b)?;
                if memory > 128_000_000 {
                    return Err(bad("observed_transcripts_memory_budget"));
                }
                evidence.push(e);
            }
        }
        csv.into_inner().finish()?;
    }
    if seen.len() != source.evidence_rows
        || source.batches_completed.iter().any(|id| {
            prepared
                .queries
                .iter()
                .any(|q| !seen.contains(&(id.clone(), q.query_id.clone())))
        })
    {
        return Err(bad("observed_transcripts_evidence_coverage"));
    }
    let mut unavailable = prepared
        .unassessed
        .iter()
        .map(|u| ObservedTranscriptUnavailable {
            file_id: u.file_id.clone(),
            record_index: u.record_index,
            read_id: u.read_id.clone(),
            query_id: None,
            reason: u.reason.clone(),
        })
        .collect::<Vec<_>>();
    let mut input = kernel::Input {
        reference: kernel::ReferenceBinding {
            assembly: selected.assembly.clone(),
            fasta_sha256: selected.fasta.sha256.clone(),
        },
        reference_gene: r.reference_gene.clone(),
        comparison_genes: r.comparison_genes.clone(),
        annotation_sources: annotations,
        queries: Vec::new(),
        unassessed_records: Vec::new(),
        required_flank_bases: r.required_flank_bases,
        budgets: r.kernel_budgets.clone(),
    };
    for e in evidence {
        let q = queries[e.query_id.as_str()];
        let file = trace
            .files_completed
            .iter()
            .find(|f| f.file_id == q.file_id)
            .ok_or_else(|| bad("observed_transcripts_record_source"))?;
        let unknown_length = e
            .segments
            .iter()
            .any(|s| !lengths.contains_key(&s.target_name));
        let missing_view = unavailable
            .iter()
            .any(|u| u.file_id == q.file_id && u.record_index == q.record_index);
        let unassessed = matches!(
            batch.accounting[&e.query_id],
            PrimerPayloadQueryAccounting::Unassessed
        );
        if unknown_length || unassessed {
            unavailable.push(ObservedTranscriptUnavailable {
                file_id: q.file_id.clone(),
                record_index: q.record_index,
                read_id: q.read_id.clone(),
                query_id: Some(q.query_id.clone()),
                reason: if unknown_length {
                    "target_length_unavailable"
                } else {
                    "producer_unassessed"
                }
                .into(),
            });
        }
        let mut ends = [
            vec![kernel::EndConstraint::PrimerBounded],
            vec![kernel::EndConstraint::ReadEndUnverified],
        ];
        for &(a, z) in &q.excluded_query_intervals {
            if a == 0 {
                ends[0].push(kernel::EndConstraint::Masked);
            }
            if z == q.query_length {
                ends[1].push(kernel::EndConstraint::Masked);
            }
        }
        input.queries.push(kernel::QueryEvidence {
            record: kernel::ReadRecord {
                source_sha256: file.observed_sha256.clone(),
                record_index: q.record_index,
                read_id: q.read_id.clone(),
            },
            alignment_sha256: batch.paf.sha256.clone(),
            query_id: q.query_id.clone(),
            query_length: q.query_length,
            scope: PrimerPayloadReferenceScope::Genome,
            target_lengths: e
                .segments
                .iter()
                .filter_map(|s| {
                    lengths
                        .get(&s.target_name)
                        .map(|n| (s.target_name.clone(), *n))
                })
                .collect(),
            masks: Some(q.excluded_query_intervals.clone()),
            ends,
            search_complete: batch.declared_complete
                && batch.competitor_search_limitations.is_empty()
                && !missing_view
                && !unknown_length
                && !unassessed,
            // Dropping only an unknown competitor would falsely strengthen the remaining hypothesis.
            segments: if unknown_length {
                Vec::new()
            } else {
                e.segments
            },
        });
    }
    if !source.batches_completed.contains(&r.batch_id) {
        for q in &prepared.queries {
            unavailable.push(ObservedTranscriptUnavailable {
                file_id: q.file_id.clone(),
                record_index: q.record_index,
                read_id: q.read_id.clone(),
                query_id: Some(q.query_id.clone()),
                reason: "batch_not_completed".into(),
            });
        }
    }
    for u in &unavailable {
        let file = trace
            .files_completed
            .iter()
            .find(|f| f.file_id == u.file_id)
            .ok_or_else(|| bad("observed_transcripts_unavailable_source"))?;
        input.unassessed_records.push(kernel::ReadRecord {
            source_sha256: file.observed_sha256.clone(),
            record_index: u.record_index,
            read_id: u.read_id.clone(),
        });
    }
    input.unassessed_records.sort();
    input.unassessed_records.dedup();
    Ok((input, unavailable, prepared.queries))
}

impl GentleEngine {
    /// Consume saved bindings only; no producer invocation or original-read/reference replay.
    pub fn summarize_observed_transcript_structures(
        &self,
        request: SummarizeObservedTranscriptStructuresRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<ObservedTranscriptStructuresReport, EngineError> {
        admit(&request)?;
        let identity = hash(&request)?;
        let mut ctl = Control::new(request.cancel_path.clone(), identity.clone(), 1, progress)?;
        let source: PrimerPayloadAlignmentReport =
            read_json(&request.source_report, &request.budgets)?;
        let bound = Self::primer_payload_alignment_result_binding(&source.request, &source)?;
        if encode(&bound)? != encode(&request.source)? {
            return Err(bad("observed_transcripts_source_binding"));
        }
        let stage = staging(&request.output_dir)?;
        let mut data = ObservedTranscriptStructuresData {
            schema: "gentle.observed_transcript_structures_data.v1".into(),
            request_sha256: identity.clone(),
            analysis: None,
            unavailable: Vec::new(),
        };
        let mut metadata = BTreeMap::<kernel::ReadRecord, Vec<PrimerPayloadQuery>>::new();
        let work = (|| {
            ctl.check()?;
            let (input, unavailable, queries) = consume(&request, &source, &mut ctl)?;
            data.unavailable = unavailable;
            let by_id = queries
                .iter()
                .map(|q| (q.query_id.as_str(), q))
                .collect::<BTreeMap<_, _>>();
            for e in &input.queries {
                metadata
                    .entry(e.record.clone())
                    .or_default()
                    .push(by_id[e.query_id.as_str()].clone());
            }
            if source.batches_completed.contains(&request.batch_id) {
                data.analysis = Some(kernel::summarize_with_control(&input, &mut || ctl.check())?);
                ctl.completed = 1;
            }
            Ok::<_, EngineError>(())
        })();
        match work {
            Err(e) if e.message == INTERRUPTED => {}
            Err(e) => return Err(e),
            Ok(()) => {}
        }
        ctl.stop(true);
        let mut left = request.budgets.max_output_bytes;
        let mut output = BoundedOutput {
            writer: HashedWriter::create(&stage.path().join("structures.json"))?,
            left,
        };
        serde_json::to_writer(&mut output, &data).map_err(|e| bad(e.to_string()))?;
        left = output.left;
        let mut artifacts = vec![output.writer.receipt("structures.json".into())?];
        let mut memory = 0;
        let mut reads = Vec::new();
        let mut families = Vec::new();
        let mut comparisons = Vec::new();
        let mut junctions = Vec::new();
        if let Some(s) = &data.analysis {
            for record in &s.records {
                let queries = metadata
                    .get(&record.record)
                    .map(Vec::as_slice)
                    .unwrap_or(&[]);
                push(
                    &mut reads,
                    json!({"support":record,"queries":queries,"observations":s.observations.iter().filter(|o| o.record==record.record).collect::<Vec<_>>() }),
                    &mut memory,
                    &request.budgets,
                )?;
            }
            for family in &s.families {
                let records = s
                    .records
                    .iter()
                    .filter(|r| r.possible_chains.contains(&family.chain_id))
                    .collect::<Vec<_>>();
                push(
                    &mut families,
                    json!({"family":family,"records":records,"support_by_context":support_groups(&records,&metadata)?,"groups_are_nonadditive":true}),
                    &mut memory,
                    &request.budgets,
                )?;
            }
            for o in &s.observations {
                for c in &o.comparisons {
                    push(
                        &mut comparisons,
                        json!({"observation_id":o.observation_id,"record":o.record,"comparison":c}),
                        &mut memory,
                        &request.budgets,
                    )?;
                }
                for (skip, qualification) in o.skips.iter().zip(&o.flank_qualification) {
                    push(
                        &mut junctions,
                        json!({"contig":o.contig,"skip":skip,"record":o.record,"observation_id":o.observation_id,"flank_qualification":qualification,
                        "on_best_path":s.records.iter().filter(|r| r.record==o.record).any(|r| r.hypotheses.iter().any(|p| p.contains(&o.observation_id))),
                        "source_comparisons":o.comparisons.iter().map(|c|json!({"source_id":c.source_id,"gene_id":c.gene_id,"resolution":c.resolution,"junction":c.junctions.as_ref().and_then(|j|j.iter().find(|j|j.skipped_interval==*skip))})).collect::<Vec<_>>(),
                        "note":"Physical skip observation, not independent molecule or natural RNA strand evidence"}),
                        &mut memory,
                        &request.budgets,
                    )?;
                }
            }
        }
        for (name, rows) in TABLES[..4]
            .iter()
            .zip([reads, families, comparisons, junctions])
        {
            artifacts.push(output_table(
                stage.path(),
                name,
                &rows,
                &mut left,
                request.budgets.max_line_bytes,
            )?);
        }
        artifacts.push(output_table(
            stage.path(),
            TABLES[4],
            &data.unavailable,
            &mut left,
            request.budgets.max_line_bytes,
        )?);
        let mut report = ObservedTranscriptStructuresReport {
            schema:REPORT.into(),content_sha256:String::new(),request_sha256:identity,request:request.clone(),
            status:if ctl.status()=="completed" && data.analysis.is_none() {"unavailable"} else {ctl.status()}.into(),control_error:ctl.error.clone(),
            source_status:source.status.clone(),source_trace_status:source.source_trace_status.clone(),source_query_status:source.source_query_status.clone(),
            cohort:source.cohort.clone(),source_batches_completed:source.batches_completed.clone(),source_batches_not_completed:source.batches_not_completed.clone(),
            source_inputs_missing:source.source_inputs_missing,source_query_units_missing:source.source_query_units_missing,source_excluded_anchor_count:source.source_excluded_anchor_count,
            counts:data.analysis.as_ref().map(|s| ObservedTranscriptCounts { records:s.records.len(),observations:s.observations.len(),chain_families:s.families.len(),unavailable_views:data.unavailable.len() }),artifacts,
            nonclaims:vec!["Internal content consistency is not authenticated provenance or an independently repeated alignment search.".into(),
                "Support counts source records, not independent molecules, expression, exclusive isoforms or proteins.".into(),
                "Novel means relative to the explicitly bound annotation snapshots; query orientation does not establish RNA strand.".into(),
                "Missing views and unknown target lengths remain unavailable; payload ends are not natural RNA ends. Original FASTQ, PAF and reference FASTA are not reopened.".into()],
        };
        report.content_sha256 = hash(&report)?;
        report_binding(&report)?;
        publish(stage, &request.output_dir, &report, &request.budgets)?;
        Ok(report)
    }

    pub fn observed_transcript_structures_result_binding(
        request: &SummarizeObservedTranscriptStructuresRequest,
        report: &ObservedTranscriptStructuresReport,
    ) -> Result<ObservedTranscriptStructuresResultBinding, EngineError> {
        if hash(request)? != report.request_sha256 {
            return Err(bad("observed_transcripts_request_binding"));
        }
        report_binding(report)
    }
}

impl GentleEngine {
    /// Inspect saved S1 output only. Row selection is bounded and verifies its table at EOF.
    pub fn inspect_observed_transcript_structures(
        path: &str,
        expected: Option<&ObservedTranscriptStructuresResultBinding>,
        table: &str,
        rows: usize,
        verify_artifacts: bool,
    ) -> Result<ObservedTranscriptStructuresInspection, EngineError> {
        if rows > 1000 || !TABLES.contains(&table) {
            return Err(bad("observed_transcripts_inspection_bounds"));
        }
        regular(Path::new(path))?;
        let report: ObservedTranscriptStructuresReport = load(path)?;
        let binding = report_binding(&report)?;
        let actual = fs::canonicalize(path).map_err(|e| bad(e.to_string()))?;
        if fs::canonicalize(&binding.report_path).ok().as_ref() != Some(&actual)
            || expected.is_some_and(|e| encode(e).ok() != encode(&binding).ok())
        {
            return Err(bad("observed_transcripts_inspection_binding"));
        }
        let root = Path::new(path)
            .parent()
            .unwrap_or(Path::new("."))
            .join("artifacts");
        let mut selected = Vec::new();
        let mut size = 0u64;
        let mut response = 0usize;
        for a in &report.artifacts {
            if !verify_artifacts && !(rows > 0 && a.path == table) {
                continue;
            }
            size = size
                .checked_add(a.bytes)
                .ok_or_else(|| bad("observed_transcripts_inspection_byte_budget"))?;
            if size > report.request.budgets.max_output_bytes {
                return Err(bad("observed_transcripts_inspection_byte_budget"));
            }
            let mut reader = VerifiedReader::open(
                &artifact_file(&root, a)?,
                report.request.budgets.max_output_bytes,
            )?;
            if a.path == table && rows > 0 {
                let mut csv = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .from_reader(reader);
                if csv
                    .headers()
                    .map_err(|e| bad(e.to_string()))?
                    .iter()
                    .collect::<Vec<_>>()
                    != ["row_json"]
                {
                    return Err(bad("observed_transcripts_table_columns"));
                }
                for row in csv.records() {
                    let row = row.map_err(|e| bad(e.to_string()))?;
                    if row.len() != 1 || row[0].len() > report.request.budgets.max_line_bytes {
                        return Err(bad("observed_transcripts_row_budget"));
                    }
                    if selected.len() < rows {
                        response += row[0].len();
                        if response > 8_000_000 {
                            return Err(bad("observed_transcripts_response_budget"));
                        }
                        selected
                            .push(serde_json::from_str(&row[0]).map_err(|e| bad(e.to_string()))?);
                    }
                }
                reader = csv.into_inner();
            } else {
                std::io::copy(&mut reader, &mut std::io::sink()).map_err(|e| bad(e.to_string()))?;
            }
            reader.finish()?;
        }
        Ok(ObservedTranscriptStructuresInspection {
            report,
            artifacts_verified: verify_artifacts,
            table: table.into(),
            rows: selected,
        })
    }

    /// Request-bound control marker; a successful write is not confirmed worker termination.
    pub fn cancel_observed_transcript_structures(
        r: &SummarizeObservedTranscriptStructuresRequest,
    ) -> Result<String, EngineError> {
        admit(r)?;
        let identity = hash(r)?;
        let path = r
            .cancel_path
            .as_deref()
            .ok_or_else(|| bad("observed_transcripts_cancel_path_required"))?;
        if control::marker_requested(Some(path), std::slice::from_ref(&identity))? {
            return Ok("Cancellation already requested; termination not confirmed".into());
        }
        let parent = Path::new(path)
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let mut file = tempfile::NamedTempFile::new_in(parent).map_err(|e| bad(e.to_string()))?;
        write!(file, "{identity}").map_err(|e| bad(e.to_string()))?;
        file.as_file().sync_all().map_err(|e| bad(e.to_string()))?;
        file.persist_noclobber(path)
            .map_err(|e| bad(e.to_string()))?;
        Ok("Cancellation requested for this exact request; termination not confirmed".into())
    }
}
