//! Verified source units and opaque anchor identities; source sequences are never rewritten.

use super::*;

pub(super) fn admit(
    request: &PreparePrimerPayloadAlignmentQueriesRequest,
) -> Result<(), EngineError> {
    budgets(&request.budgets)?;
    if request.schema != PRIMER_PAYLOAD_PREPARE_SCHEMA
        || !safe_id(&request.cohort.cohort_id)
        || request.cohort.label.trim().is_empty()
        || request.cohort.label.len() > 512
        || request.cohort.selection_provenance.trim().is_empty()
        || request.cohort.selection_provenance.len() > 4096
        || request.file_ids.is_empty()
        || request.file_ids.len() > 4096
        || request.min_payload_bases == 0
        || request.min_payload_bases > request.budgets.max_query_bases
        || (request.source_mode == PrimerPayloadSourceMode::RetainedPrimaryFasta
            && request.anchor_selection != PrimerPayloadAnchorSelection::Primary)
    {
        return Err(bad("primer_payload_prepare_admission"));
    }
    Ok(())
}

pub(super) fn trace_source(
    binding: &RnaPrimerTraceResultBinding,
    b: &PrimerPayloadBudgets,
) -> Result<RnaPrimerTraceReport, EngineError> {
    if regular(Path::new(&binding.report_path))?.len() > b.max_input_bytes.min(128_000_000) {
        return Err(bad("primer_payload_trace_report_budget"));
    }
    let report =
        *GentleEngine::inspect_primer_trace(&binding.report_path, Some(binding), None, 0)?.report;
    if report
        .artifacts
        .iter()
        .try_fold(regular(Path::new(&binding.report_path))?.len(), |n, a| {
            n.checked_add(a.bytes)
        })
        .is_none_or(|n| n > b.max_input_bytes)
    {
        return Err(bad("primer_payload_trace_metadata_budget"));
    }
    let parent = Path::new(&binding.report_path)
        .parent()
        .unwrap_or(Path::new("."));
    let root = parent.join("artifacts");
    let meta = fs::symlink_metadata(&root).map_err(|e| bad(e.to_string()))?;
    if !meta.is_dir() || meta.file_type().is_symlink() {
        return Err(bad("primer_payload_trace_artifacts_directory"));
    }
    for artifact in &report.artifacts {
        let bytes = verified_bytes(&artifact_file(&root, artifact)?, b)?;
        let value: serde_json::Value =
            serde_json::from_slice(&bytes).map_err(|e| bad(e.to_string()))?;
        match artifact.path.as_str() {
            "primer_panel.json" if encode(&value)? == encode(&report.primer_panel)? => {}
            "summary.json" if encode(&value)? == encode(&report.summary)? => {}
            "reference_panel.json" => {
                let reference: RnaReadPanel =
                    serde_json::from_value(value).map_err(|e| bad(e.to_string()))?;
                validate_panel(&reference)?;
                if report.reference_panel_sha256.as_ref() != Some(&reference.content_sha256) {
                    return Err(bad("primer_payload_trace_reference_binding"));
                }
            }
            _ => return Err(bad("primer_payload_trace_metadata_binding")),
        }
    }
    Ok(report)
}

pub(super) fn source_rows(
    root: &Path,
    source: &RnaPrimerTraceFile,
    panel: &RnaPrimerPanel,
    b: &PrimerPayloadBudgets,
) -> Result<Vec<RnaPrimerTraceRow>, EngineError> {
    let artifact = source
        .artifacts
        .iter()
        .find(|a| a.path == format!("{}/reads.tsv", source.file_id))
        .ok_or_else(|| bad("primer_payload_trace_rows_missing"))?;
    let bytes = verified_bytes(&artifact_file(root, artifact)?, b)?;
    let mut csv = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(bytes.as_slice());
    if csv
        .headers()
        .map_err(|e| bad(e.to_string()))?
        .iter()
        .collect::<Vec<_>>()
        != primer_trace::HEADER
    {
        return Err(bad("primer_payload_trace_columns"));
    }
    let ids = panel
        .primers
        .iter()
        .map(|p| p.primer_id.as_str())
        .collect::<BTreeSet<_>>();
    let mut rows = Vec::new();
    let mut counts = RnaPrimerTraceCounts::default();
    for record in csv.records() {
        let record = record.map_err(|e| bad(e.to_string()))?;
        if record.iter().map(str::len).sum::<usize>() > b.max_line_bytes
            || rows.len() >= b.max_records
        {
            return Err(bad("primer_payload_trace_row_budget"));
        }
        let row = primer_trace::parse_row(&record)?;
        if row.record_index != rows.len() as u64
            || row.read_id.is_empty()
            || row.read_length > input::MAX_READ_BASES
            || row.hits.len()
                != record[16]
                    .parse::<usize>()
                    .map_err(|_| bad("primer_payload_trace_hit_count"))?
            || row.hits.len() > 1024
            || row.hits.iter().any(|h| {
                h.start >= h.end
                    || h.end > row.read_length
                    || h.primer_id
                        .as_ref()
                        .is_some_and(|p| !ids.contains(p.as_str()))
                    || h.core_start
                        .zip(h.core_end)
                        .is_some_and(|(a, z)| a < h.start || a >= z || z > h.end)
                    || h.core_start.is_some() != h.core_end.is_some()
                    || (h.primer_id.is_some() && h.core_start.is_none())
            })
            || row.anchor.as_ref().is_some_and(|a| !row.hits.contains(a))
        {
            return Err(bad("primer_payload_trace_row_binding"));
        }
        primer_trace::count(&mut counts, &row, row.read_length);
        rows.push(row);
    }
    if encode(&counts)? != encode(&source.counts)? {
        return Err(bad("primer_payload_trace_count_binding"));
    }
    Ok(rows)
}

pub(super) fn retained_fasta(
    root: &Path,
    source: &RnaPrimerTraceFile,
    b: &PrimerPayloadBudgets,
) -> Result<BTreeMap<u64, Vec<u8>>, EngineError> {
    let artifact = source
        .artifacts
        .iter()
        .find(|a| a.path == format!("{}/payload.fasta", source.file_id))
        .ok_or_else(|| bad("primer_payload_retained_fasta_missing_use_explicit_replay"))?;
    let file = artifact_file(root, artifact)?;
    let mut reader = BufReader::new(VerifiedReader::open(
        &file,
        b.max_input_bytes.min(128_000_000),
    )?);
    let mut result = BTreeMap::new();
    let mut current = None;
    while let Some(text) = line(&mut reader, b.max_line_bytes)? {
        if let Some(header) = text.strip_prefix('>') {
            let (id, index) = header
                .split_ascii_whitespace()
                .next()
                .unwrap_or("")
                .rsplit_once(':')
                .ok_or_else(|| bad("primer_payload_fasta_identity"))?;
            let index = index
                .parse::<u64>()
                .map_err(|_| bad("primer_payload_fasta_identity"))?;
            if id != source.file_id
                || result.len() >= b.max_records
                || result.insert(index, Vec::new()).is_some()
            {
                return Err(bad("primer_payload_fasta_identity"));
            }
            current = Some(index);
        } else {
            let index = current.ok_or_else(|| bad("primer_payload_fasta_header"))?;
            let sequence = result.get_mut(&index).expect("inserted header");
            if !text.bytes().all(|c| b"ACGTRYSWKMBDHVN".contains(&c))
                || sequence
                    .len()
                    .checked_add(text.len())
                    .is_none_or(|n| n > b.max_query_bases)
            {
                return Err(bad("primer_payload_fasta_bases_or_budget"));
            }
            sequence.extend_from_slice(text.as_bytes());
        }
    }
    reader.into_inner().finish()?;
    Ok(result)
}

pub(super) struct PayloadGeometry {
    pub interval: (usize, usize),
    pub excluded: Vec<(usize, usize)>,
    pub usable: usize,
}

/// Clip observed hits before subtraction, using the same transform for writing and verification.
pub(super) fn payload_geometry(
    length: usize,
    anchor: &RnaPrimerHit,
    hits: &[RnaPrimerHit],
) -> Result<PayloadGeometry, EngineError> {
    if anchor.primer_id.is_none()
        || !hits.contains(anchor)
        || hits.iter().any(|h| h.start >= h.end || h.end > length)
    {
        return Err(bad("primer_payload_anchor_geometry"));
    }
    let (a, z) = if anchor.strand == PrimerReadStrand::Forward {
        (anchor.end, length)
    } else {
        (0, anchor.start)
    };
    let excluded = union(
        hits.iter()
            .filter_map(|other| {
                let s = a.max(other.start);
                let e = z.min(other.end);
                if s >= e {
                    return None;
                }
                Some(if anchor.strand == PrimerReadStrand::Forward {
                    (s - a, e - a)
                } else {
                    (z - e, z - s)
                })
            })
            .collect(),
    );
    let usable = (z - a)
        .checked_sub(span_bases(&excluded))
        .ok_or_else(|| bad("primer_payload_anchor_geometry"))?;
    Ok(PayloadGeometry {
        interval: (a, z),
        excluded,
        usable,
    })
}

fn add_row(
    row: &RnaPrimerTraceRow,
    source: &RnaPrimerTraceFile,
    request: &PreparePrimerPayloadAlignmentQueriesRequest,
    raw: Option<&[u8]>,
    retained: Option<&[u8]>,
    queries: &mut Vec<PrimerPayloadQuery>,
    unassessed: &mut Vec<PrimerPayloadUnassessed>,
    writer: &mut impl Write,
) -> Result<(), EngineError> {
    let hits = match request.anchor_selection {
        PrimerPayloadAnchorSelection::Primary => row.anchor.iter().collect::<Vec<_>>(),
        PrimerPayloadAnchorSelection::All => {
            row.hits.iter().filter(|h| h.primer_id.is_some()).collect()
        }
    };
    let mut groups = BTreeMap::new();
    for hit in hits.into_iter().filter(|h| h.primer_id.is_some()) {
        groups
            .entry((hit.start, hit.end, hit.strand, hit.core_start, hit.core_end))
            .or_insert(hit);
    }
    if groups.is_empty() {
        unassessed.push(PrimerPayloadUnassessed {
            file_id: source.file_id.clone(),
            record_index: row.record_index,
            read_id: row.read_id.clone(),
            anchor_id: None,
            primer_ids: Vec::new(),
            reason: "no_attributed_anchor".into(),
        });
        if unassessed.len() > request.budgets.max_records {
            return Err(bad("primer_payload_query_budget"));
        }
    }
    for ((start, end, strand, core_start, core_end), hit) in groups {
        let primer_ids = row
            .hits
            .iter()
            .filter(|h| {
                (h.start, h.end, h.strand, h.core_start, h.core_end)
                    == (start, end, strand, core_start, core_end)
            })
            .filter_map(|h| h.primer_id.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let anchor_id = hash(&(
            &request.trace.report_content_sha256,
            &source.file_id,
            row.record_index,
            &row.read_id,
            (start, end, strand, core_start, core_end),
            &primer_ids,
        ))?;
        let PayloadGeometry {
            interval: (a, z),
            excluded,
            usable,
        } = payload_geometry(row.read_length, hit, &row.hits)?;
        if usable < request.min_payload_bases
            || (!excluded.is_empty()
                && request.internal_primer_policy == PrimerPayloadInternalPrimerPolicy::Unassessed)
        {
            unassessed.push(PrimerPayloadUnassessed {
                file_id: source.file_id.clone(),
                record_index: row.record_index,
                read_id: row.read_id.clone(),
                anchor_id: Some(anchor_id),
                primer_ids,
                reason: if usable < request.min_payload_bases {
                    "insufficient_primer_free_payload"
                } else {
                    "additional_primer_interval_unassessed"
                }
                .into(),
            });
            if unassessed.len() > request.budgets.max_records {
                return Err(bad("primer_payload_query_budget"));
            }
            continue;
        }
        let mut sequence = if let Some(raw) = raw {
            if raw.len() != row.read_length {
                return Err(bad("primer_payload_read_length_binding"));
            }
            if strand == PrimerReadStrand::Forward {
                raw[a..z].to_vec()
            } else {
                GentleEngine::reverse_complement_bytes(&raw[a..z])
            }
        } else {
            if row.payload != Some((a, z)) {
                return Err(bad("primer_payload_primary_interval_binding"));
            }
            retained
                .ok_or_else(|| bad("primer_payload_retained_query_missing"))?
                .to_vec()
        };
        if sequence.len() != z - a || sequence.len() > request.budgets.max_query_bases {
            return Err(bad("primer_payload_query_length_or_budget"));
        }
        for &(s, e) in &excluded {
            sequence[s..e].fill(b'N');
        }
        let sequence_sha256 = sha256_prefixed_bytes(&sequence);
        let query_id = format!(
            "q_{}",
            hash(&(&anchor_id, &sequence_sha256))?.trim_start_matches("sha256:")
        );
        writeln!(
            writer,
            ">{query_id}\n{}",
            String::from_utf8_lossy(&sequence)
        )
        .map_err(|e| bad(e.to_string()))?;
        queries.push(PrimerPayloadQuery {
            query_id,
            anchor_id,
            file_id: source.file_id.clone(),
            record_index: row.record_index,
            read_id: row.read_id.clone(),
            sample_id: source.input.sample_id.clone(),
            barcode: source.input.barcode.clone(),
            priming: source.input.priming.clone(),
            primer_ids,
            anchor: hit.clone(),
            raw_read_length: row.read_length,
            raw_payload_interval: (a, z),
            excluded_query_intervals: excluded,
            query_length: sequence.len(),
            usable_bases: usable,
            sequence_sha256,
        });
        if queries.len() > request.budgets.max_queries
            || unassessed.len() > request.budgets.max_records
        {
            return Err(bad("primer_payload_query_budget"));
        }
    }
    Ok(())
}

impl GentleEngine {
    /// Verified sidecars or explicitly approved FASTQ replay, never implicit source discovery.
    pub fn prepare_primer_payload_alignment_queries(
        &self,
        request: PreparePrimerPayloadAlignmentQueriesRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<PrimerPayloadQueries, EngineError> {
        admit(&request)?;
        let identity = hash(&request)?;
        let mut ctl = Control::new(
            request.cancel_path.clone(),
            identity.clone(),
            request.file_ids.len(),
            progress,
        )?;
        let source = trace_source(&request.trace, &request.budgets)?;
        let ids = request.file_ids.iter().cloned().collect::<BTreeSet<_>>();
        let available = source
            .files_completed
            .iter()
            .map(|f| f.file_id.clone())
            .collect::<BTreeSet<_>>();
        if ids.len() != request.file_ids.len()
            || !ids.is_subset(&available)
            || (request.cohort.all_verified_units && ids != available)
        {
            return Err(bad("primer_payload_source_selection"));
        }
        let input_bytes = source
            .artifacts
            .iter()
            .chain(
                source
                    .files_completed
                    .iter()
                    .filter(|f| ids.contains(&f.file_id))
                    .flat_map(|f| &f.artifacts),
            )
            .try_fold(
                regular(Path::new(&request.trace.report_path))?.len(),
                |sum, a| sum.checked_add(a.bytes),
            )
            .ok_or_else(|| bad("primer_payload_input_budget"))?;
        let replay_bytes = if request.source_mode == PrimerPayloadSourceMode::ReplayVerifiedFastq {
            source
                .files_completed
                .iter()
                .filter(|f| ids.contains(&f.file_id))
                .try_fold(0u64, |sum, f| sum.checked_add(f.compressed_bytes))
                .ok_or_else(|| bad("primer_payload_input_budget"))?
        } else {
            0
        };
        if input_bytes
            .checked_add(replay_bytes)
            .is_none_or(|n| n > request.budgets.max_input_bytes)
        {
            return Err(bad("primer_payload_input_budget"));
        }
        let stage = staging(&request.output_dir)?;
        let mut fasta = HashedWriter::create(&stage.path().join("queries.fasta"))?;
        let root = Path::new(&request.trace.report_path)
            .parent()
            .unwrap_or(Path::new("."))
            .join("artifacts");
        let mut queries = Vec::new();
        let mut unassessed = Vec::new();
        let mut completed = Vec::new();
        let mut records = 0usize;
        for id in &ids {
            if ctl.stop(true) {
                break;
            }
            let unit = source
                .files_completed
                .iter()
                .find(|f| &f.file_id == id)
                .expect("admitted ID");
            let rows = source_rows(&root, unit, &source.primer_panel, &request.budgets)?;
            records = records
                .checked_add(rows.len())
                .filter(|n| *n <= request.budgets.max_records)
                .ok_or_else(|| bad("primer_payload_record_budget"))?;
            let mut unit_queries = Vec::new();
            let mut unit_unassessed = Vec::new();
            let mut spool =
                tempfile::NamedTempFile::new_in(stage.path()).map_err(|e| bad(e.to_string()))?;
            let result = match request.source_mode {
                PrimerPayloadSourceMode::RetainedPrimaryFasta => {
                    let mut retained = retained_fasta(&root, unit, &request.budgets)?;
                    let mut result = Ok(());
                    for (index, row) in rows.iter().enumerate() {
                        if ctl.stop(index == 0) {
                            result = Err(bad(INTERRUPTED));
                            break;
                        }
                        let sequence = retained.remove(&row.record_index);
                        if let Err(e) = add_row(
                            row,
                            unit,
                            &request,
                            None,
                            sequence.as_deref(),
                            &mut unit_queries,
                            &mut unit_unassessed,
                            &mut spool,
                        ) {
                            result = Err(e);
                            break;
                        }
                        if spool
                            .as_file()
                            .metadata()
                            .map_err(|e| bad(e.to_string()))?
                            .len()
                            > request.budgets.max_output_bytes
                        {
                            return Err(bad("primer_payload_output_budget"));
                        }
                    }
                    if result.is_ok() && !retained.is_empty() {
                        return Err(bad("primer_payload_foreign_fasta_query"));
                    }
                    result
                }
                PrimerPayloadSourceMode::ReplayVerifiedFastq => {
                    regular(Path::new(&unit.input.path))?;
                    let mut index = 0usize;
                    let visited = visit_fastq(Path::new(&unit.input.path), &mut |record| {
                        if ctl.stop(index == 0) {
                            return Ok(false);
                        }
                        let row = rows
                            .get(index)
                            .ok_or_else(|| bad("primer_payload_unexpected_source_record"))?;
                        if row.read_id != record.read_id || row.read_length != record.sequence.len()
                        {
                            return Err(bad("primer_payload_source_record_binding"));
                        }
                        add_row(
                            row,
                            unit,
                            &request,
                            Some(&record.sequence),
                            None,
                            &mut unit_queries,
                            &mut unit_unassessed,
                            &mut spool,
                        )?;
                        if spool
                            .as_file()
                            .metadata()
                            .map_err(|e| bad(e.to_string()))?
                            .len()
                            > request.budgets.max_output_bytes
                        {
                            return Err(bad("primer_payload_output_budget"));
                        }
                        index += 1;
                        Ok(true)
                    })?;
                    if let Some(digest) = visited {
                        if index != rows.len()
                            || digest.sha256 != unit.observed_sha256
                            || digest.md5 != unit.observed_md5
                            || digest.bytes != unit.compressed_bytes
                        {
                            return Err(bad("primer_payload_original_input_binding"));
                        }
                        Ok(())
                    } else {
                        Err(bad(INTERRUPTED))
                    }
                }
            };
            match result {
                Ok(()) => {}
                Err(e) if e.message == INTERRUPTED => break,
                Err(e) => return Err(e),
            }
            if queries
                .len()
                .checked_add(unit_queries.len())
                .is_none_or(|n| n > request.budgets.max_queries)
                || unassessed
                    .len()
                    .checked_add(unit_unassessed.len())
                    .is_none_or(|n| n > request.budgets.max_records)
            {
                return Err(bad("primer_payload_query_budget"));
            }
            spool
                .as_file_mut()
                .flush()
                .map_err(|e| bad(e.to_string()))?;
            let mut input = fs::File::open(spool.path()).map_err(|e| bad(e.to_string()))?;
            std::io::copy(&mut input, &mut fasta).map_err(|e| bad(e.to_string()))?;
            drop(input);
            drop(spool);
            if fs::metadata(stage.path().join("queries.fasta"))
                .map_err(|e| bad(e.to_string()))?
                .len()
                > request.budgets.max_output_bytes
            {
                return Err(bad("primer_payload_output_budget"));
            }
            queries.extend(unit_queries);
            unassessed.extend(unit_unassessed);
            completed.push(id.clone());
            ctl.completed = completed.len();
        }
        ctl.stop(true);
        let artifact = fasta.receipt("queries.fasta".into())?;
        let mut report = PrimerPayloadQueries {
            schema: PRIMER_PAYLOAD_QUERIES_SCHEMA.into(),
            content_sha256: String::new(),
            request_sha256: identity,
            request: request.clone(),
            status: ctl.status().into(),
            control_error: ctl.error.clone(),
            source_trace_status: source.status,
            source_inputs_missing: source.inputs_not_completed.len(),
            original_fastq_reverified: request.source_mode
                == PrimerPayloadSourceMode::ReplayVerifiedFastq
                && !completed.is_empty(),
            primer_panel: source.primer_panel,
            files_not_completed: ids
                .into_iter()
                .filter(|id| !completed.contains(id))
                .collect(),
            files_completed: completed,
            queries,
            unassessed,
            artifacts: vec![artifact],
            nonclaims: nonclaims(),
        };
        report.content_sha256 = hash(&report)?;
        publish(stage, &request.output_dir, &report, &request.budgets)?;
        Ok(report)
    }
}
