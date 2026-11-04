//! Joined saved-evidence projection. No original FASTQ, PAF, reference or annotation replay.

use super::*;

fn admit(r: &InspectPrimerPayloadJunctionsRequest) -> Result<(), EngineError> {
    budgets(&r.budgets)?;
    let b = &r.budgets;
    if r.schema != PRIMER_JUNCTION_INSPECTION_REQUEST_SCHEMA
        || r.source.kind != "alignments"
        || Path::new(&r.source_report.path) != Path::new(&r.source.report_path)
        || !is_hash(&r.source_report.sha256)
        || r.rows > 1000
        || r.offset > 4096
        || b.max_input_bytes > 128_000_000
        || b.max_output_bytes > 8_000_000
        || b.max_records > 4096
        || b.max_queries > 4096
        || b.max_line_bytes > 1_000_000
        || b.max_alignments_per_query > 128
        || r.file_id
            .iter()
            .chain(r.primer_id.iter())
            .any(|s| !safe_id(s))
    {
        return Err(bad("primer_junction_inspection_admission"));
    }
    Ok(())
}

fn read_table<T: serde::de::DeserializeOwned>(
    f: &PrimerPayloadFile,
    column: &str,
    b: &PrimerPayloadBudgets,
    ctl: &mut Control<'_>,
) -> Result<Vec<T>, EngineError> {
    let reader = VerifiedReader::open(f, b.max_input_bytes)?;
    let mut csv = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    if csv
        .headers()
        .map_err(|e| bad(e.to_string()))?
        .iter()
        .collect::<Vec<_>>()
        != [column]
    {
        return Err(bad("primer_junction_inspection_columns"));
    }
    let mut out = Vec::new();
    for row in csv.records() {
        ctl.check()?;
        let row = row.map_err(|e| bad(e.to_string()))?;
        if row.len() != 1 || row[0].len() > b.max_line_bytes || out.len() >= b.max_records {
            return Err(bad("primer_junction_inspection_table_budget"));
        }
        out.push(serde_json::from_str(&row[0]).map_err(|e| bad(e.to_string()))?);
    }
    csv.into_inner().finish()?;
    Ok(out)
}

fn validate_evidence(
    e: &PrimerPayloadReadEvidence,
    q: &PrimerPayloadQuery,
    batch: &PrimerPayloadAlignmentBatch,
    source: &PrimerPayloadAlignmentReport,
    b: &PrimerPayloadBudgets,
) -> Result<(), EngineError> {
    let accounting = batch
        .accounting
        .get(&e.query_id)
        .ok_or_else(|| bad("primer_junction_inspection_accounting"))?;
    if e.reference_id != batch.reference_id
        || e.file_id != q.file_id
        || e.record_index != q.record_index
        || e.read_id != q.read_id
        || e.primer_ids != q.primer_ids
        || e.sample_id != q.sample_id
        || e.barcode != q.barcode
        || e.priming != q.priming
        || e.usable_bases != q.usable_bases
        || (matches!(accounting, PrimerPayloadQueryAccounting::Mapped) == e.segments.is_empty())
        || e.expected_target_bases.is_some_and(|n| n > q.usable_bases)
        || e.expected_target_fraction_per_mille
            != evidence::target_fraction(e.expected_target_bases, q.usable_bases)
        || e.segments.len() > b.max_alignments_per_query
        || e.internal_query_boundaries.iter().any(|edge| {
            edge.left_query_end > edge.right_query_start
                || edge.right_query_start > q.query_length
                || !matches!(edge.left.strand, '+' | '-')
                || !matches!(edge.right.strand, '+' | '-')
        })
        || e.primer_edge_alternatives.len() > b.max_alignments_per_query
        || !matches!(
            e.classification.as_str(),
            "unassessed"
                | "producer_reported_unmapped"
                | "unresolved"
                | "ambiguous"
                | "expected_dominant"
                | "off_target_dominant"
                | "mixed"
        )
        || (q.primer_ids.len() != 1
            && matches!(
                e.classification.as_str(),
                "expected_dominant" | "off_target_dominant" | "mixed"
            ))
        || (*accounting == PrimerPayloadQueryAccounting::Unassessed
            && e.classification != "unassessed")
        || (*accounting == PrimerPayloadQueryAccounting::ProducerReportedUnmapped
            && e.classification != "producer_reported_unmapped")
        || (*accounting == PrimerPayloadQueryAccounting::Mapped
            && matches!(
                e.classification.as_str(),
                "unassessed" | "producer_reported_unmapped"
            ))
        || e.segments.iter().any(|s| {
            s.query_interval.0 >= s.query_interval.1
                || s.query_interval.1 > q.query_length
                || s.target_interval.0 >= s.target_interval.1
                || !matches!(s.strand, '+' | '-')
                || s.reported_identity_per_mille > 1000
                || s.aligned_query_intervals
                    .iter()
                    .any(|&(a, z)| a >= z || z > q.query_length)
        })
        || e.primer_edge_alternatives
            .iter()
            .any(|edge| !matches!(edge.strand, '+' | '-'))
    {
        return Err(bad("primer_junction_inspection_evidence_binding"));
    }
    let keys = e
        .primer_edge_alternatives
        .iter()
        .flat_map(|edge| {
            e.primer_ids.iter().map(move |id| {
                (
                    id.clone(),
                    edge.target_name.clone(),
                    edge.coordinate,
                    edge.strand,
                )
            })
        })
        .collect::<BTreeSet<_>>();
    let mut seen = BTreeSet::new();
    for (key, seeds) in &e.terminal_seeds {
        let key: (String, String, usize, char) =
            serde_json::from_str(key).map_err(|_| bad("primer_junction_inspection_seed_key"))?;
        if !keys.contains(&key)
            || !seen.insert(key)
            || seeds.iter().map(|s| s.length).collect::<Vec<_>>()
                != source.request.policy.terminal_seed_lengths
            || seeds.iter().any(|s| {
                s.mismatches.is_some_and(|n| n > s.length)
                    || (s.mismatches.is_some() == s.reason.is_some())
                    || s.reason.as_ref().is_some_and(|r| r.trim().is_empty())
            })
        {
            return Err(bad("primer_junction_inspection_seed_binding"));
        }
    }
    if seen != keys {
        return Err(bad("primer_junction_inspection_seed_coverage"));
    }
    Ok(())
}

impl GentleEngine {
    /// Verify and join saved evidence before pagination; no source mutation or implicit biological reassignment.
    pub fn inspect_primer_payload_junctions(
        &self,
        request: InspectPrimerPayloadJunctionsRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<PrimerPayloadJunctionInspection, EngineError> {
        admit(&request)?;
        let b = &request.budgets;
        let identity = hash(&request)?;
        let mut ctl = Control::new(None, identity.clone(), 1, progress)?;
        if ctl.stop(true) {
            return Err(bad("primer_junction_inspection_cancelled"));
        }
        let source: PrimerPayloadAlignmentReport = read_json(&request.source_report, b)?;
        let actual = Self::primer_payload_alignment_result_binding(&source.request, &source)?;
        if encode(&actual)? != encode(&request.source)?
            || source.evidence_rows > b.max_records
            || source.request.references.len() > 16
            || source.request.batches.len() > 16
            || source.profiles.len() > b.max_records
            || source.breakpoint_clusters.len() > b.max_records
        {
            return Err(bad("primer_junction_inspection_source_binding_or_budget"));
        }
        let canonical =
            fs::canonicalize(&request.source_report.path).map_err(|e| bad(e.to_string()))?;
        if fs::canonicalize(&actual.report_path).ok().as_ref() != Some(&canonical) {
            return Err(bad("primer_junction_inspection_source_location"));
        }
        let additional = source
            .artifacts
            .iter()
            .try_fold(request.source_report.bytes, |n, a| {
                n.checked_add(a.bytes).filter(|n| *n <= b.max_input_bytes)
            })
            .ok_or_else(|| bad("primer_junction_inspection_input_budget"))?;
        ctl.check()?;
        let reserved = additional
            .checked_add(source.request.queries.bytes)
            .filter(|n| *n < b.max_input_bytes)
            .ok_or_else(|| bad("primer_junction_inspection_input_budget"))?;
        let mut trace_budget = b.clone();
        trace_budget.max_input_bytes = b.max_input_bytes - reserved;
        let trace = prepare::trace_source(&source.request.trace, &trace_budget)?;
        ctl.check()?;
        let (prepared, trace_rows, query_bases) = binding::queries_with_rows(
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
            || encode(&source.cohort)? != encode(&prepared.request.cohort)?
            || request
                .file_id
                .as_ref()
                .is_some_and(|id| !prepared.request.file_ids.contains(id))
            || request.primer_id.as_ref().is_some_and(|id| {
                !trace
                    .primer_panel
                    .primers
                    .iter()
                    .any(|p| &p.primer_id == id)
            })
        {
            return Err(bad(
                "primer_junction_inspection_source_accounting_or_filter",
            ));
        }
        let units = trace
            .files_completed
            .iter()
            .map(|u| (u.file_id.as_str(), u))
            .collect::<BTreeMap<_, _>>();
        let mut hashes = BTreeSet::new();
        let mut cohort_reads = 0usize;
        for (id, rows) in &trace_rows {
            let unit = units
                .get(id.as_str())
                .ok_or_else(|| bad("primer_junction_inspection_unit"))?;
            if !hashes.insert(&unit.observed_sha256) {
                return Err(bad("primer_junction_inspection_duplicate_source_bytes"));
            }
            cohort_reads = cohort_reads
                .checked_add(rows.len())
                .filter(|n| *n <= b.max_records)
                .ok_or_else(|| bad("primer_junction_inspection_record_budget"))?;
        }
        let root = Path::new(&request.source_report.path)
            .parent()
            .unwrap_or(Path::new("."))
            .join("artifacts");
        let find = |name: &str| {
            source
                .artifacts
                .iter()
                .find(|a| a.path == name)
                .ok_or_else(|| bad("primer_junction_inspection_artifact"))
        };
        let all: Vec<PrimerPayloadReadEvidence> = read_table(
            &artifact_file(&root, find("read_anchor_evidence.tsv")?)?,
            "evidence_json",
            b,
            &mut ctl,
        )?;
        let saved_profiles: Vec<PrimerPayloadPartnerProfile> = read_table(
            &artifact_file(&root, find("partner_profiles.tsv")?)?,
            "partner_profile_json",
            b,
            &mut ctl,
        )?;
        let saved_clusters: Vec<PrimerPayloadBreakpointCluster> = read_table(
            &artifact_file(&root, find("breakpoint_clusters.tsv")?)?,
            "breakpoint_cluster_json",
            b,
            &mut ctl,
        )?;
        let queries = prepared
            .queries
            .iter()
            .map(|q| (q.query_id.as_str(), q))
            .collect::<BTreeMap<_, _>>();
        let batches = source
            .request
            .batches
            .iter()
            .map(|b| (b.batch_id.as_str(), b))
            .collect::<BTreeMap<_, _>>();
        let mut seen = BTreeSet::new();
        for e in &all {
            ctl.check()?;
            let q = queries
                .get(e.query_id.as_str())
                .ok_or_else(|| bad("primer_junction_inspection_foreign_query"))?;
            let batch = batches
                .get(e.batch_id.as_str())
                .ok_or_else(|| bad("primer_junction_inspection_foreign_batch"))?;
            if !source.batches_completed.contains(&e.batch_id)
                || !seen.insert((&e.batch_id, &e.query_id))
            {
                return Err(bad(
                    "primer_junction_inspection_duplicate_or_incomplete_batch",
                ));
            }
            validate_evidence(e, q, batch, &source, b)?;
        }
        if seen.len() != source.evidence_rows
            || source.batches_completed.iter().any(|id| {
                prepared
                    .queries
                    .iter()
                    .any(|q| !seen.contains(&(id, &q.query_id)))
            })
        {
            return Err(bad("primer_junction_inspection_evidence_population"));
        }
        ctl.check()?;
        let refs = source
            .request
            .batches
            .iter()
            .filter(|batch| source.batches_completed.contains(&batch.batch_id))
            .map(|batch| batch.reference_id.clone())
            .collect::<Vec<_>>();
        let (profiles, clusters) = evidence::aggregate(
            &all,
            &prepared,
            &trace,
            &refs,
            source.request.policy.breakpoint_tolerance_bases,
        );
        if encode(&profiles)? != encode(&source.profiles)?
            || encode(&clusters)? != encode(&source.breakpoint_clusters)?
            || encode(&profiles)? != encode(&saved_profiles)?
            || encode(&clusters)? != encode(&saved_clusters)?
        {
            return Err(bad(
                "primer_junction_inspection_profile_or_cluster_accounting",
            ));
        }
        let mut by_read = BTreeMap::<(&str, u64), Vec<&PrimerPayloadReadEvidence>>::new();
        let mut by_query = BTreeMap::<(&str, u64), Vec<&PrimerPayloadQuery>>::new();
        let mut by_missing = BTreeMap::<(&str, u64), Vec<&PrimerPayloadUnassessed>>::new();
        for e in &all {
            by_read
                .entry((&e.file_id, e.record_index))
                .or_default()
                .push(e);
        }
        for q in &prepared.queries {
            by_query
                .entry((&q.file_id, q.record_index))
                .or_default()
                .push(q);
        }
        for u in &prepared.unassessed {
            by_missing
                .entry((&u.file_id, u.record_index))
                .or_default()
                .push(u);
        }
        let mut selected = Vec::new();
        let mut matched = 0usize;
        for id in &prepared.files_completed {
            let unit = units[id.as_str()];
            for row in &trace_rows[id] {
                ctl.check()?;
                if request.file_id.as_ref().is_some_and(|f| f != id)
                    || request
                        .primer_id
                        .as_ref()
                        .is_some_and(|p| !row.hits.iter().any(|h| h.primer_id.as_ref() == Some(p)))
                {
                    continue;
                }
                let index = matched;
                matched += 1;
                if index < request.offset || selected.len() >= request.rows {
                    continue;
                }
                let key = (id.as_str(), row.record_index);
                selected.push(PrimerPayloadJunctionRow {
                    file_id: id.clone(),
                    source_sha256: unit.observed_sha256.clone(),
                    sample_id: unit.input.sample_id.clone(),
                    barcode: unit.input.barcode.clone(),
                    priming: unit.input.priming.clone(),
                    trace: row.clone(),
                    queries: by_query
                        .get(&key)
                        .into_iter()
                        .flatten()
                        .map(|q| (**q).clone())
                        .collect(),
                    payload_previews: by_query
                        .get(&key)
                        .into_iter()
                        .flatten()
                        .map(|q| {
                            let bases = &query_bases[&q.query_id];
                            (
                                q.query_id.clone(),
                                String::from_utf8_lossy(&bases[..bases.len().min(64)]).into_owned(),
                            )
                        })
                        .collect(),
                    evidence: by_read
                        .get(&key)
                        .into_iter()
                        .flatten()
                        .map(|e| (**e).clone())
                        .collect(),
                    unassessed: by_missing
                        .get(&key)
                        .into_iter()
                        .flatten()
                        .map(|u| (**u).clone())
                        .collect(),
                });
            }
        }
        let mut report = PrimerPayloadJunctionInspection { schema:PRIMER_JUNCTION_INSPECTION_SCHEMA.into(),content_sha256:String::new(),
            request_sha256:identity, request:request.clone(), source_status:source.status, source_control_error:source.control_error,
            trace_status:trace.status, trace_control_error:trace.control_error, query_status:prepared.status,query_control_error:prepared.control_error,
            cohort:source.cohort,primer_panel:trace.primer_panel, trace_inputs_not_completed:trace.inputs_not_completed,
            query_files_not_completed:prepared.files_not_completed,batches_not_completed:source.batches_not_completed,
            cohort_distinct_reads:cohort_reads,matched_trace_rows:matched,rows_truncated:matched>request.offset.saturating_add(selected.len()),rows:selected,
            profiles,breakpoint_clusters:clusters,saved_artifacts_verified:true, original_reads_reverified:false, reference_bases_reverified:false,
            warnings:vec![
                "Joined saved trace/query/payload evidence only. Original FASTQ, original reference/annotation and PAF were not reopened; no new trimming, alignment or junction mechanism was inferred.".into(),
                "Trace intervals are original-read half-open coordinates. Query bases/interval masks were already normalized into the observed primer direction; do not reverse them again or infer transcript orientation.".into(),
                "Each query anchor is a geometry-group representative. Full/core classes and edit counts belong to the individual trace hits; a core-only hit does not prove that a full synthetic primer was observed or trimmed.".into(),
                "Payload previews contain at most 64 saved query bases in their existing masked primer orientation, not original-read or reference bases. Terminal seed mismatch counts/reasons are saved evidence; literal reference flank letters are unavailable and have not been synthesized or reverified.".into(),
                "Distinct-read keys are (file_id,record_index), not read names. Duplicate source-byte hashes in the verified selected cohort are refused. Anchor observations, ambiguous equivalence sets and overlapping denominator buckets remain separate.".into(),
                "Row filters/pagination change only the trace preview. Profiles/clusters retain full saved cohort/reference/primer-set scope; never sum alternative partner groups or reference scopes as independent molecules.".into(),
                "Successful inspection does not complete cancelled/control-error sources or fill missing units. Producer/search completeness is not authenticated or experimentally validated. No RT mechanism, expression, fusion, ordering or candidate redesign claim.".into(),
            ] };
        if encode(&report)?.len() as u64 > b.max_output_bytes {
            return Err(bad("primer_junction_inspection_response_budget"));
        }
        ctl.completed = 1;
        if ctl.stop(true) {
            return Err(bad("primer_junction_inspection_cancelled"));
        }
        report.content_sha256 = hash(&report)?;
        if encode(&report)?.len() as u64 > b.max_output_bytes {
            return Err(bad("primer_junction_inspection_response_budget"));
        }
        Ok(report)
    }
}

#[cfg(test)]
#[path = "junction_inspector_tests.rs"]
mod tests;
