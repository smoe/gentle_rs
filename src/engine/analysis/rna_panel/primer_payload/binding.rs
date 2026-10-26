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

fn validate_anchor_accounting(
    report: &PrimerPayloadQueries,
    by_unit: &BTreeMap<String, Vec<RnaPrimerTraceRow>>,
) -> Result<(), EngineError> {
    let mut queries = BTreeMap::new();
    let mut excluded = BTreeMap::new();
    let mut unattributed = BTreeMap::new();
    for q in &report.queries {
        if queries.insert(q.anchor_id.as_str(), q).is_some() {
            return Err(bad("primer_payload_duplicate_anchor_accounting"));
        }
    }
    for missing in &report.unassessed {
        if let Some(id) = missing.anchor_id.as_deref() {
            if queries.contains_key(id) || excluded.insert(id, missing).is_some() {
                return Err(bad("primer_payload_duplicate_anchor_accounting"));
            }
        } else if unattributed
            .insert((missing.file_id.as_str(), missing.record_index), missing)
            .is_some()
        {
            return Err(bad("primer_payload_duplicate_anchor_accounting"));
        }
    }
    // Match each source observation to exactly one policy-derived disposition; never trust a reason label.
    for (file, rows) in by_unit {
        for row in rows {
            let anchors = prepare::selected_anchors(row, report.request.anchor_selection);
            if anchors.is_empty() {
                let missing = unattributed
                    .remove(&(file.as_str(), row.record_index))
                    .ok_or_else(|| bad("primer_payload_anchor_accounting_coverage"))?;
                if missing.read_id != row.read_id
                    || !missing.primer_ids.is_empty()
                    || missing.reason != "no_attributed_anchor"
                {
                    return Err(bad("primer_payload_exclusion_binding"));
                }
            }
            for hit in anchors {
                let primers = prepare::anchor_primers(row, hit);
                let id = prepare::anchor_identity(
                    &report.request.trace.report_content_sha256,
                    file,
                    row,
                    hit,
                    &primers,
                )?;
                let geometry = prepare::payload_geometry(row.read_length, hit, &row.hits)?;
                if let Some(reason) = prepare::exclusion_reason(&report.request, &geometry) {
                    let missing = excluded
                        .remove(id.as_str())
                        .ok_or_else(|| bad("primer_payload_exclusion_policy_binding"))?;
                    if missing.file_id != *file
                        || missing.record_index != row.record_index
                        || missing.read_id != row.read_id
                        || missing.primer_ids != primers
                        || missing.reason != reason
                    {
                        return Err(bad("primer_payload_exclusion_binding"));
                    }
                } else {
                    if excluded.contains_key(id.as_str()) {
                        return Err(bad("primer_payload_exclusion_policy_binding"));
                    }
                    if queries.remove(id.as_str()).is_none() {
                        return Err(bad("primer_payload_anchor_accounting_coverage"));
                    }
                }
            }
        }
    }
    if !queries.is_empty() || !excluded.is_empty() || !unattributed.is_empty() {
        return Err(bad("primer_payload_anchor_accounting_coverage"));
    }
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
    {
        return Err(bad("primer_payload_query_source_binding"));
    }
    if report.queries.len() > b.max_queries || report.unassessed.len() > b.max_records {
        return Err(bad("primer_payload_query_import_budget"));
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
    let mut records = 0usize;
    for id in &report.files_completed {
        let unit = trace
            .files_completed
            .iter()
            .find(|f| &f.file_id == id)
            .ok_or_else(|| bad("primer_payload_query_unit"))?;
        let rows = prepare::source_rows(&trace_root, unit, &trace.primer_panel, b)?;
        records = records
            .checked_add(rows.len())
            .filter(|n| *n <= b.max_records && *n <= report.request.budgets.max_records)
            .ok_or_else(|| bad("primer_payload_trace_record_budget"))?;
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
        let prepare::PayloadGeometry {
            interval,
            excluded: mask,
            usable,
        } = prepare::payload_geometry(row.read_length, h, &row.hits)?;
        let ids = prepare::anchor_primers(row, h);
        let anchor = prepare::anchor_identity(&trace.content_sha256, &q.file_id, row, h, &ids)?;
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
            || q.usable_bases != usable
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
    validate_anchor_accounting(&report, &seen_rows)?;
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
