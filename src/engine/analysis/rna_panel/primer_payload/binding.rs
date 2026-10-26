//! Content-bound prepared queries and bounded inspection; never reopens original FASTQ implicitly.

use super::*;

pub(super) fn outcome(status: &str, error: &Option<String>) -> bool {
    error.as_ref().is_none_or(|s| !s.trim().is_empty())
        && matches!(
            (status, error.is_some()),
            ("completed" | "cancelled", false) | ("control_error", true)
        )
}

pub(super) fn validate_queries(report: &PrimerPayloadQueries) -> Result<(), EngineError> {
    prepare::admit(&report.request)?;
    let mut basis = report.clone();
    basis.content_sha256.clear();
    let requested = report.request.file_ids.iter().collect::<BTreeSet<_>>();
    let done = report.files_completed.iter().collect::<BTreeSet<_>>();
    let missing = report.files_not_completed.iter().collect::<BTreeSet<_>>();
    if report.schema != PRIMER_PAYLOAD_QUERIES_SCHEMA
        || hash(&basis)? != report.content_sha256
        || report.request_sha256 != hash(&report.request)?
        || report.request.schema != PRIMER_PAYLOAD_PREPARE_SCHEMA
        || !outcome(&report.status, &report.control_error)
        || requested.len() != report.request.file_ids.len()
        || done.len() != report.files_completed.len()
        || missing.len() != report.files_not_completed.len()
        || !done.is_disjoint(&missing)
        || done.union(&missing).copied().collect::<BTreeSet<_>>() != requested
        || (report.status == "completed" && !missing.is_empty())
        || report.artifacts.len() != 1
        || report.artifacts[0].path != "queries.fasta"
        || report.queries.len() > report.request.budgets.max_queries
        || report.unassessed.len() > report.request.budgets.max_records
        || report.original_fastq_reverified
            != (report.request.source_mode == PrimerPayloadSourceMode::ReplayVerifiedFastq
                && !done.is_empty())
    {
        return Err(bad("primer_payload_query_manifest_binding"));
    }
    budgets(&report.request.budgets)?;
    Ok(())
}

pub(super) fn queries(
    file: &PrimerPayloadFile,
    trace: &RnaPrimerTraceReport,
    expected_trace: &RnaPrimerTraceResultBinding,
    b: &PrimerPayloadBudgets,
    additional_input_bytes: u64,
) -> Result<PrimerPayloadQueries, EngineError> {
    let report: PrimerPayloadQueries = read_json(file, b)?;
    validate_queries(&report)?;
    if encode(&report.request.trace)? != encode(expected_trace)? {
        return Err(bad("primer_payload_query_manifest_source_binding"));
    }
    let requested = report.request.file_ids.iter().collect::<BTreeSet<_>>();
    let available = trace
        .files_completed
        .iter()
        .map(|u| &u.file_id)
        .collect::<BTreeSet<_>>();
    if !requested.is_subset(&available)
        || (report.request.cohort.all_verified_units && requested != available)
    {
        return Err(bad("primer_payload_source_selection"));
    }
    let report_bytes = regular(Path::new(&report.request.trace.report_path))?.len();
    let mut bytes = additional_input_bytes
        .checked_add(file.bytes)
        .and_then(|n| n.checked_add(report_bytes))
        .ok_or_else(|| bad("primer_payload_input_budget"))?;
    for n in report
        .artifacts
        .iter()
        .map(|a| a.bytes)
        .chain(trace.artifacts.iter().map(|a| a.bytes))
        .chain(
            trace
                .files_completed
                .iter()
                .filter(|u| report.files_completed.contains(&u.file_id))
                .flat_map(|u| u.artifacts.iter().map(|a| a.bytes)),
        )
    {
        bytes = bytes
            .checked_add(n)
            .filter(|n| *n <= b.max_input_bytes)
            .ok_or_else(|| bad("primer_payload_input_budget"))?;
    }
    if bytes > b.max_input_bytes {
        return Err(bad("primer_payload_input_budget"));
    }
    if report.request.trace.report_content_sha256 != trace.content_sha256
        || report.source_trace_status != trace.status
        || report.source_inputs_missing != trace.inputs_not_completed.len()
        || encode(&report.primer_panel)? != encode(&trace.primer_panel)?
        || report.queries.len() > b.max_queries
    {
        return Err(bad("primer_payload_query_source_binding"));
    }
    let expected = Path::new(&report.request.output_dir).join("report.json");
    if fs::canonicalize(&expected).ok() != fs::canonicalize(&file.path).ok() {
        return Err(bad("primer_payload_query_location_binding"));
    }
    let root = Path::new(&file.path)
        .parent()
        .unwrap_or(Path::new("."))
        .join("artifacts");
    let bytes = verified_bytes(&artifact_file(&root, &report.artifacts[0])?, b)?;
    let mut reader = BufReader::new(bytes.as_slice());
    let mut seqs = BTreeMap::<String, Vec<u8>>::new();
    let mut current = None;
    while let Some(text) = line(&mut reader, b.max_line_bytes)? {
        if let Some(id) = text.strip_prefix('>') {
            if !id.starts_with("q_")
                || !valid_digest(&id[2..], 64)
                || seqs.len() >= b.max_queries
                || seqs.insert(id.into(), Vec::new()).is_some()
            {
                return Err(bad("primer_payload_query_fasta_id"));
            }
            current = Some(id.to_string());
        } else {
            let id = current
                .as_ref()
                .ok_or_else(|| bad("primer_payload_query_fasta_header"))?;
            let seq = seqs.get_mut(id).expect("inserted");
            if !text.bytes().all(|x| b"ACGTRYSWKMBDHVN".contains(&x))
                || seq
                    .len()
                    .checked_add(text.len())
                    .is_none_or(|n| n > b.max_query_bases)
            {
                return Err(bad("primer_payload_query_fasta_bases"));
            }
            seq.extend(text.as_bytes());
        }
    }
    let trace_root = Path::new(&report.request.trace.report_path)
        .parent()
        .unwrap_or(Path::new("."))
        .join("artifacts");
    let mut seen = BTreeSet::new();
    let mut anchors = BTreeSet::new();
    let mut seen_rows = BTreeMap::new();
    let mut retained = BTreeMap::new();
    for id in &report.files_completed {
        let unit = trace
            .files_completed
            .iter()
            .find(|f| &f.file_id == id)
            .ok_or_else(|| bad("primer_payload_query_unit"))?;
        let rows = prepare::source_rows(&trace_root, unit, &trace.primer_panel, b)?;
        if report.request.source_mode == PrimerPayloadSourceMode::RetainedPrimaryFasta {
            retained.insert(id.clone(), prepare::retained_fasta(&trace_root, unit, b)?);
        }
        seen_rows.insert(id.clone(), rows);
    }
    for q in &report.queries {
        let unit = trace
            .files_completed
            .iter()
            .find(|f| f.file_id == q.file_id)
            .ok_or_else(|| bad("primer_payload_query_unit"))?;
        let row = seen_rows
            .get(&q.file_id)
            .and_then(|rows| {
                usize::try_from(q.record_index)
                    .ok()
                    .and_then(|i| rows.get(i))
            })
            .ok_or_else(|| bad("primer_payload_query_source_row"))?;
        let h = &q.anchor;
        let interval = if h.strand == PrimerReadStrand::Forward {
            (h.end, row.read_length)
        } else {
            (0, h.start)
        };
        let mask = union(
            row.hits
                .iter()
                .filter_map(|x| {
                    let (a, z) = (interval.0.max(x.start), interval.1.min(x.end));
                    (a < z).then_some(if h.strand == PrimerReadStrand::Forward {
                        (a - interval.0, z - interval.0)
                    } else {
                        (interval.1 - z, interval.1 - a)
                    })
                })
                .collect(),
        );
        let ids = row
            .hits
            .iter()
            .filter(|x| {
                (x.start, x.end, x.strand, x.core_start, x.core_end)
                    == (h.start, h.end, h.strand, h.core_start, h.core_end)
            })
            .filter_map(|x| x.primer_id.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let anchor = hash(&(
            &trace.content_sha256,
            &q.file_id,
            q.record_index,
            &q.read_id,
            (h.start, h.end, h.strand, h.core_start, h.core_end),
            &ids,
        ))?;
        let seq = seqs
            .get(&q.query_id)
            .ok_or_else(|| bad("primer_payload_query_sequence_missing"))?;
        if !seen.insert(&q.query_id)
            || !anchors.insert(&q.anchor_id)
            || q.anchor_id != anchor
            || q.query_id
                != format!(
                    "q_{}",
                    hash(&(&anchor, &q.sequence_sha256))?.trim_start_matches("sha256:")
                )
            || q.primer_ids != ids
            || ids.is_empty()
            || !row.hits.contains(h)
            || q.read_id != row.read_id
            || q.raw_read_length != row.read_length
            || q.raw_payload_interval != interval
            || q.excluded_query_intervals != mask
            || q.query_length != interval.1 - interval.0
            || seq.len() != q.query_length
            || q.sequence_sha256 != sha256_prefixed_bytes(seq)
            || q.usable_bases != q.query_length - span_bases(&mask)
            || q.usable_bases < report.request.min_payload_bases
            || q.sample_id != unit.input.sample_id
            || q.barcode != unit.input.barcode
            || q.priming != unit.input.priming
            || mask
                .iter()
                .any(|(a, z)| seq[*a..*z].iter().any(|c| *c != b'N'))
            || (report.request.anchor_selection == PrimerPayloadAnchorSelection::Primary
                && row.anchor.as_ref() != Some(h))
        {
            return Err(bad("primer_payload_query_geometry_or_identity"));
        }
        if let Some(unit) = retained.get(&q.file_id) {
            let mut original = unit
                .get(&q.record_index)
                .ok_or_else(|| bad("primer_payload_retained_query_missing"))?
                .clone();
            if original.len() != q.query_length {
                return Err(bad("primer_payload_retained_query_length"));
            }
            for &(a, z) in &mask {
                original[a..z].fill(b'N');
            }
            if &original != seq {
                return Err(bad("primer_payload_retained_query_bases"));
            }
        }
    }
    if seen.len() != seqs.len() {
        return Err(bad("primer_payload_query_fasta_coverage"));
    }
    for missing in &report.unassessed {
        let row = seen_rows
            .get(&missing.file_id)
            .and_then(|r| {
                usize::try_from(missing.record_index)
                    .ok()
                    .and_then(|i| r.get(i))
            })
            .ok_or_else(|| bad("primer_payload_unassessed_row"))?;
        if missing.read_id != row.read_id
            || missing.reason.trim().is_empty()
            || missing
                .primer_ids
                .iter()
                .any(|id| !row.hits.iter().any(|h| h.primer_id.as_ref() == Some(id)))
        {
            return Err(bad("primer_payload_unassessed_binding"));
        }
    }
    // A rehashed manifest cannot silently drop a source row/anchor or turn it into an off-target zero.
    let mut wanted = BTreeSet::new();
    let mut no_anchor = BTreeSet::new();
    for (file, rows) in &seen_rows {
        for row in rows {
            let hits = if report.request.anchor_selection == PrimerPayloadAnchorSelection::Primary {
                row.anchor.iter().collect::<Vec<_>>()
            } else {
                row.hits.iter().filter(|h| h.primer_id.is_some()).collect()
            };
            let groups = hits
                .into_iter()
                .filter(|h| h.primer_id.is_some())
                .map(|h| (h.start, h.end, h.strand, h.core_start, h.core_end))
                .collect::<BTreeSet<_>>();
            if groups.is_empty() {
                no_anchor.insert((file.clone(), row.record_index));
            }
            for geometry in groups {
                let ids = row
                    .hits
                    .iter()
                    .filter(|h| (h.start, h.end, h.strand, h.core_start, h.core_end) == geometry)
                    .filter_map(|h| h.primer_id.clone())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect::<Vec<_>>();
                wanted.insert(hash(&(
                    &trace.content_sha256,
                    file,
                    row.record_index,
                    &row.read_id,
                    geometry,
                    &ids,
                ))?);
            }
        }
    }
    let mut accounted = report
        .queries
        .iter()
        .map(|q| q.anchor_id.clone())
        .collect::<BTreeSet<_>>();
    let mut missing_no_anchor = BTreeSet::new();
    for missing in &report.unassessed {
        if let Some(id) = &missing.anchor_id {
            if !accounted.insert(id.clone()) {
                return Err(bad("primer_payload_duplicate_anchor_accounting"));
            }
        } else if !missing_no_anchor.insert((missing.file_id.clone(), missing.record_index))
            || missing.reason != "no_attributed_anchor"
        {
            return Err(bad("primer_payload_missing_anchor_accounting"));
        }
    }
    if accounted != wanted || no_anchor != missing_no_anchor {
        return Err(bad("primer_payload_anchor_accounting_coverage"));
    }
    Ok(report)
}

impl GentleEngine {
    pub fn primer_payload_query_result_binding(
        request: &PreparePrimerPayloadAlignmentQueriesRequest,
        report: &PrimerPayloadQueries,
    ) -> Result<PrimerPayloadResultBinding, EngineError> {
        validate_queries(report)?;
        if hash(request)? != report.request_sha256 {
            return Err(bad("primer_payload_request_binding"));
        }
        Ok(PrimerPayloadResultBinding {
            schema: PRIMER_PAYLOAD_BINDING_SCHEMA.into(),
            kind: "queries".into(),
            request_sha256: report.request_sha256.clone(),
            report_content_sha256: report.content_sha256.clone(),
            report_path: Path::new(&request.output_dir)
                .join("report.json")
                .display()
                .to_string(),
            outcome: report.status.clone(),
        })
    }
}
