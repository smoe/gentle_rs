//! S2.1b: verify one bounded GTF snapshot and publish a separate, report-last receipt.

use super::*;
use crate::observed_transcript_structures::annotation as kernel;

#[path = "annotation_prepare/publication.rs"]
mod publication;
#[cfg(test)]
#[path = "annotation_prepare/tests.rs"]
mod tests;

const CAP: u64 = 128 * 1024 * 1024;

fn text(s: &str, cap: usize) -> bool {
    !s.trim().is_empty() && s.len() <= cap && !s.chars().any(char::is_control)
}

fn limits(r: &PrepareObservedTranscriptAnnotationRequest) -> Result<kernel::Limits, EngineError> {
    let b = &r.budgets;
    if r.schema != OBSERVED_ANNOTATION_REQUEST_SCHEMA
        || !text(&r.reference_id, 4096)
        || !text(&r.reference.assembly, 512)
        || !is_hash(&r.reference.fasta_sha256)
        || [&r.selector.source_id, &r.selector.gene_id]
            .iter()
            .any(|s| !text(s, 512) || s.chars().any(char::is_whitespace))
        || [
            &r.source.provider,
            &r.source.declared_provider_release,
            &r.source.declared_derivation,
            &r.source.source_uri,
            &r.source.notice,
        ]
        .iter()
        .any(|s| s.as_ref().is_some_and(|s| !text(s, 4096)))
        || !(1..=CAP).contains(&b.max_input_bytes)
        || !(1..=CAP).contains(&b.max_output_bytes)
        || !(1..=8 * 1024 * 1024).contains(&b.max_dictionary_bytes)
        || !(1..=1024 * 1024).contains(&b.max_line_bytes)
        || !(1..=1_000_000).contains(&b.max_lines)
        || !(1..=100_000).contains(&b.max_dictionary_entries)
        || !(1..=100_000).contains(&b.max_indexed_genes)
        || !(1..=100_000).contains(&b.max_indexed_transcripts)
        || !(1..=1_000_000).contains(&b.max_exon_rows)
        || !(1..=CAP as usize).contains(&b.max_retained_bytes)
        || !(1..=1_000_000_000).contains(&b.max_work)
        || b.max_diagnostics > 1024
        || r.policy.minimum_reference_gap_bases == 0
        || [&r.gtf, &r.dictionary].iter().any(|f| !is_hash(&f.sha256))
        || r.gtf
            .bytes
            .checked_add(r.dictionary.bytes)
            .is_none_or(|n| n > b.max_input_bytes)
        || r.dictionary.bytes > b.max_dictionary_bytes as u64
    {
        return Err(bad("observed_annotation_admission"));
    }
    Ok(kernel::Limits {
        max_input_bytes: usize::try_from(r.gtf.bytes)
            .map_err(|_| bad("observed_annotation_integer_range"))?
            .max(1),
        max_line_bytes: b.max_line_bytes,
        max_lines: b.max_lines,
        max_dictionary_entries: b.max_dictionary_entries,
        max_indexed_genes: b.max_indexed_genes,
        max_indexed_transcripts: b.max_indexed_transcripts,
        max_exon_rows: b.max_exon_rows,
        max_retained_bytes: b.max_retained_bytes,
        max_work: b.max_work,
        max_diagnostics: b.max_diagnostics,
    })
}

/// Conservative JSON allocation admission, not allocator/RSS telemetry.
/// Charge raw JSON at 64x plus parser scratch before allocating/deserializing it.
fn json_charge(bytes: u64, memory: usize) -> Result<usize, EngineError> {
    usize::try_from(bytes)
        .ok()
        .and_then(|n| n.checked_mul(64))
        .and_then(|n| n.checked_add(65536))
        .filter(|n| *n <= memory)
        .ok_or_else(|| bad("observed_annotation_json_memory_budget"))
}

struct Counter {
    bytes: u64,
    limit: u64,
}
impl Write for Counter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len() as u64)
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| std::io::Error::other("observed_annotation_json_budget"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn semantic_hash<T: Serialize>(value: &T, memory: usize) -> Result<String, EngineError> {
    let mut counter = Counter {
        bytes: 0,
        limit: memory.saturating_sub(65536) as u64 / 64,
    };
    serde_json::to_writer(&mut counter, value).map_err(|e| bad(e.to_string()))?;
    json_charge(counter.bytes, memory)?;
    crate::digest_utils::canonical_json_string(value)
        .map(|s| sha256_prefixed_bytes(s.as_bytes()))
        .map_err(|e| bad(e.to_string()))
}

fn dictionary(
    r: &PrepareObservedTranscriptAnnotationRequest,
    paths: &publication::Paths,
) -> Result<(kernel::Dictionary, usize), EngineError> {
    json_charge(r.dictionary.bytes, r.budgets.max_retained_bytes)?;
    let mut reader = paths.reader(&r.dictionary, r.budgets.max_dictionary_bytes as u64)?;
    // Bound allocation by verified declared bytes; never use an unbounded read_to_end.
    let mut raw = vec![0; r.dictionary.bytes as usize];
    reader
        .read_exact(&mut raw)
        .map_err(|e| bad(e.to_string()))?;
    let mut extra = [0];
    if reader.read(&mut extra).map_err(|e| bad(e.to_string()))? != 0 {
        return Err(bad("observed_annotation_dictionary_growth"));
    }
    reader.finish()?;
    paths.unchanged(&r.dictionary.path)?;
    let d: ReferenceContigDictionary = serde_json::from_slice(&raw)
        .map_err(|e| bad(format!("observed_annotation_dictionary_json: {e}")))?;
    if d.schema != REFERENCE_CONTIG_DICTIONARY_SCHEMA
        || d.reference_id != r.reference_id
        || d.assembly != r.reference.assembly
        || d.reference_fasta_sha256 != r.reference.fasta_sha256
        || d.contigs.is_empty()
        || d.contigs.len() > r.budgets.max_dictionary_entries
    {
        return Err(bad("observed_annotation_dictionary_binding"));
    }
    let retained = d.contigs.iter().try_fold(4096usize, |n, c| {
        n.checked_add(c.name.len())
            .and_then(|n| n.checked_add(512))
            .ok_or_else(|| bad("observed_annotation_dictionary_memory"))
    })?;
    let d = kernel::Dictionary {
        reference_id: d.reference_id,
        reference: r.reference.clone(),
        contigs: d
            .contigs
            .into_iter()
            .map(|c| kernel::Contig {
                name: c.name,
                length: c.length,
            })
            .collect(),
    };
    Ok((d, retained))
}

fn raw_line(reader: &mut impl BufRead, limit: usize) -> Result<Option<Vec<u8>>, EngineError> {
    let mut raw = Vec::new();
    let n = reader
        .take(limit as u64 + 1)
        .read_until(b'\n', &mut raw)
        .map_err(|e| bad(e.to_string()))?;
    if raw.len() > limit {
        return Err(bad("observed_annotation_line_budget"));
    }
    Ok((n != 0).then_some(raw))
}

fn analysis(
    p: kernel::Preparation,
    dictionary_entries: usize,
    retained: usize,
    dictionary_bytes: usize,
) -> ObservedAnnotationAnalysis {
    ObservedAnnotationAnalysis {
        grammar: p.grammar.into(),
        outcome: p.outcome,
        coordinate_authority: p.coordinate_authority.into(),
        snapshot_release: format!("snapshot-{}", p.raw_sha256),
        raw_sha256: p.raw_sha256,
        raw_bytes: p.raw_bytes,
        dictionary_entries,
        counts: p.counts,
        selected_reasons: p.selected_reasons,
        unrelated_reasons: p.unrelated_reasons,
        diagnostics: p.diagnostics,
        diagnostic_count: p.diagnostic_count,
        omitted_diagnostics: p.omitted_diagnostics,
        indexed_genes: p.indexed_genes,
        indexed_transcripts: p.indexed_transcripts,
        gene_version: p.gene_version,
        transcript_versions: p.transcript_versions,
        included_contigs: p.included_contigs,
        work: p.work + dictionary_bytes,
        retained_peak_bytes: p.retained_peak_bytes + retained,
    }
}

fn prepare(
    r: &PrepareObservedTranscriptAnnotationRequest,
    paths: &publication::Paths,
    limits: &mut kernel::Limits,
    control: &mut Control<'_>,
    dictionary_verified: &mut bool,
    gtf_verified: &mut bool,
) -> Result<(kernel::Preparation, usize, usize), EngineError> {
    control.check()?;
    if r.dictionary.bytes as usize >= limits.max_work {
        return Err(bad("observed_annotation_work_budget"));
    }
    limits.max_work -= r.dictionary.bytes as usize;
    let (dictionary, retained) = dictionary(r, paths)?;
    let retained = retained
        .checked_add(65536)
        .ok_or_else(|| bad("observed_annotation_dictionary_memory"))?;
    *dictionary_verified = true;
    limits.max_retained_bytes = limits
        .max_retained_bytes
        .checked_sub(retained)
        .filter(|n| *n > 0)
        .ok_or_else(|| bad("observed_annotation_dictionary_memory"))?;
    let context = kernel::Context {
        reference_id: r.reference_id.clone(),
        reference: r.reference.clone(),
        selector: r.selector.clone(),
        provider: r.source.provider.clone(),
        declared_provider_release: r.source.declared_provider_release.clone(),
        scope: r.source.scope,
        policy: r.policy.clone(),
    };
    control.completed = 1;
    if control.stop(true) {
        return Err(bad(INTERRUPTED));
    }
    let mut consumer =
        kernel::Consumer::new(&context, &dictionary, limits, &mut || control.check())?;
    let mut reader =
        BufReader::with_capacity(8192, paths.reader(&r.gtf, r.budgets.max_input_bytes)?);
    loop {
        control.check()?;
        let Some(raw) = raw_line(&mut reader, r.budgets.max_line_bytes)? else {
            break;
        };
        consumer.push_line(&raw, &mut || control.check())?;
    }
    reader.into_inner().finish()?;
    paths.unchanged(&r.gtf.path)?;
    *gtf_verified = true;
    control.completed = 2;
    if control.stop(true) {
        return Err(bad(INTERRUPTED));
    }
    let p = consumer.finish(&mut || control.check())?;
    if p.raw_sha256 != r.gtf.sha256 || p.raw_bytes as u64 != r.gtf.bytes {
        return Err(bad("observed_annotation_consumer_binding"));
    }
    Ok((p, dictionary.contigs.len(), retained))
}

impl GentleEngine {
    /// Direct S2.1b file boundary only. No Operation, shell/MCP route or project mutation.
    /// Callers must authorize the external write attempt, even for a refusal/cancellation.
    pub fn prepare_observed_transcript_annotation(
        &self,
        request: PrepareObservedTranscriptAnnotationRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<ObservedTranscriptAnnotationPreparation, EngineError> {
        let mut limits = limits(&request)?;
        let paths = publication::Paths::new(&request)?;
        let identity = semantic_hash(&request, request.budgets.max_retained_bytes)?;
        // Unlike legacy Control::new(Some(...)), a pre-start malformed marker is
        // recorded as a control failure, not confused with invalid data or cancellation.
        let mut control = Control::new(None, identity.clone(), 3, progress)?;
        control.marker = request.cancel_path.clone();
        let mut report = ObservedTranscriptAnnotationReport {
            schema: OBSERVED_ANNOTATION_REPORT_SCHEMA.into(), content_sha256: String::new(),
            request_sha256: identity, request: request.clone(),
            status: ObservedAnnotationExecutionStatus::Completed, control_error: None,
            dictionary_verified: false, gtf_verified: false, analysis: None, annotation: None,
            nonclaims: vec![
                "EOF verifies supplied bytes, not provider completeness, authentic provenance or the declared FASTA.".into(),
                "Gene absence is exact-ID absence in this snapshot, not biological absence; no isoform/protein identity is inferred.".into(),
                "Retained-memory charges are conservative representation bounds, not allocator/RSS measurements.".into(),
                "Legacy annotation readers do not enforce this preparation receipt or its policies.".into(),
            ],
        };
        let stage = staging(&request.output_dir)?;
        let prepared = if control.stop(true) {
            Err(bad(INTERRUPTED))
        } else {
            prepare(
                &request,
                &paths,
                &mut limits,
                &mut control,
                &mut report.dictionary_verified,
                &mut report.gtf_verified,
            )
        };
        match prepared {
            Ok((mut p, dictionary_entries, retained)) => {
                report
                    .nonclaims
                    .extend(p.nonclaims.iter().map(|s| (*s).into()));
                // Touch provenance fields so this bridge cannot silently forget its context.
                if p.selector.source_id != request.selector.source_id
                    || p.selector.gene_id != request.selector.gene_id
                    || p.scope != request.source.scope
                    || p.policy != request.policy
                    || p.provider != request.source.provider
                    || p.declared_provider_release != request.source.declared_provider_release
                {
                    return Err(bad("observed_annotation_context_binding"));
                }
                if let Some(annotation) = p.annotation.take() {
                    report.annotation = Some(publication::write_json(
                        &stage.path().join("annotation.json"),
                        "artifacts/annotation.json",
                        &annotation,
                        request.budgets.max_output_bytes,
                    )?);
                }
                report.analysis = Some(analysis(
                    p,
                    dictionary_entries,
                    retained,
                    request.dictionary.bytes as usize,
                ));
                control.completed = 3;
                control.stop(true);
            }
            Err(e) if e.message == INTERRUPTED && control.stopped => {}
            Err(e) => return Err(e),
        }
        report.status = if control.error.is_some() {
            ObservedAnnotationExecutionStatus::ControlError
        } else if control.stopped {
            ObservedAnnotationExecutionStatus::Cancelled
        } else {
            ObservedAnnotationExecutionStatus::Completed
        };
        report.control_error = control.error;
        report.content_sha256 = semantic_hash(&report, request.budgets.max_retained_bytes)?;
        let report_file = publication::publish(stage, &paths, &report, &mut |_| Ok(()))?;
        Ok(ObservedTranscriptAnnotationPreparation {
            report,
            report_file,
        })
    }
}
