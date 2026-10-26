//! Verified alignment batches publish as new artifacts, without project or Phase 1 mutation.

use super::*;

fn admission(request: &SummarizePrimerPayloadAlignmentsRequest) -> Result<(), EngineError> {
    budgets(&request.budgets)?;
    let p = &request.policy;
    if request.schema != PRIMER_PAYLOAD_REQUEST_SCHEMA
        || !is_hash(&request.query_manifest_sha256)
        || request.references.is_empty()
        || request.references.len() > 64
        || request.batches.len() != request.references.len()
        || request.expected_targets.len() > 100_000
        || p.min_segment_bases == 0
        || p.min_segment_bases > request.budgets.max_query_bases
        || p.min_identity_per_mille > 1000
        || p.min_mapq > 254
        || !(1..=1000).contains(&p.dominant_fraction_per_mille)
        || p.min_cross_locus_segment_bases == 0
        || p.min_cross_locus_segment_bases > request.budgets.max_query_bases
        || p.breakpoint_tolerance_bases > request.budgets.max_query_bases
        || p.terminal_seed_lengths.len() > 64
        || p.terminal_seed_lengths
            .iter()
            .any(|n| !(1..=1024).contains(n))
        || !p.terminal_seed_lengths.windows(2).all(|w| w[0] < w[1])
    {
        return Err(bad("primer_payload_summary_admission"));
    }
    let mut ids = BTreeSet::new();
    for reference in &request.references {
        if !safe_id(&reference.reference_id)
            || !ids.insert(&reference.reference_id)
            || reference.assembly.trim().is_empty()
            || reference.release.trim().is_empty()
        {
            return Err(bad("primer_payload_reference_identity"));
        }
    }
    let mut batches = BTreeSet::new();
    let mut refs = BTreeSet::new();
    for batch in &request.batches {
        if !safe_id(&batch.batch_id)
            || !batches.insert(&batch.batch_id)
            || !refs.insert(&batch.reference_id)
            || !ids.contains(&batch.reference_id)
            || batch.query_manifest_sha256 != request.query_manifest_sha256
            || batch.producer.trim().is_empty()
            || batch.producer_version.trim().is_empty()
            || batch.arguments.len() > 4096
            || batch
                .arguments
                .iter()
                .any(|a| a.len() > 4096 || a.contains('\0'))
            || batch.producer_exit_code != 0
            || batch.competitor_search_limitations.trim().is_empty()
        {
            return Err(bad("primer_payload_batch_identity_or_provenance"));
        }
    }
    Ok(())
}

fn crosswalk(
    request: &SummarizePrimerPayloadAlignmentsRequest,
    panel: &RnaPrimerPanel,
    catalogs: &BTreeMap<String, reference::Catalog>,
) -> Result<(), EngineError> {
    let mut seen = BTreeSet::new();
    let mut comparisons = 0usize;
    for expected in &request.expected_targets {
        let primer = panel
            .primers
            .iter()
            .find(|p| p.primer_id == expected.primer_id)
            .ok_or_else(|| bad("primer_payload_crosswalk_primer"))?;
        let cat = catalogs
            .get(&expected.reference_id)
            .ok_or_else(|| bad("primer_payload_crosswalk_reference"))?;
        comparisons = comparisons
            .checked_add(
                cat.annotation
                    .transcripts
                    .len()
                    .checked_mul(expected.gene_ids.len() + expected.transcript_ids.len())
                    .ok_or_else(|| bad("primer_payload_crosswalk_comparison_budget"))?,
            )
            .filter(|n| *n <= request.budgets.max_chain_comparisons)
            .ok_or_else(|| bad("primer_payload_crosswalk_comparison_budget"))?;
        if !primer.expected_targets.iter().any(|t| {
            t.locus_id == expected.phase1_locus_id
                && t.transcript_feature_id == expected.phase1_transcript_feature_id
        }) || !seen.insert((
            &expected.primer_id,
            &expected.reference_id,
            &expected.phase1_locus_id,
            expected.phase1_transcript_feature_id,
        )) || (expected.gene_ids.is_empty() && expected.transcript_ids.is_empty())
            || !expected.gene_ids.windows(2).all(|w| w[0] < w[1])
            || !expected.transcript_ids.windows(2).all(|w| w[0] < w[1])
            || expected.gene_ids.iter().any(|g| {
                !cat.annotation
                    .transcripts
                    .iter()
                    .any(|t| t.gene_ids.contains(g))
            })
            || expected.transcript_ids.iter().any(|id| {
                !cat.annotation.transcripts.iter().any(|t| {
                    &t.transcript_id == id
                        && (expected.gene_ids.is_empty()
                            || t.gene_ids.iter().any(|g| expected.gene_ids.contains(g)))
                })
            })
        {
            return Err(bad("primer_payload_crosswalk_source_identity"));
        }
    }
    Ok(())
}

fn batch_rows(
    batch: &PrimerPayloadAlignmentBatch,
    reference: &PrimerPayloadReference,
    cat: &reference::Catalog,
    prepared: &PrimerPayloadQueries,
    request: &SummarizePrimerPayloadAlignmentsRequest,
    ctl: &mut Control<'_>,
    work: &mut usize,
) -> Result<(Vec<PrimerPayloadReadEvidence>, BTreeMap<String, usize>), EngineError> {
    let queries = prepared
        .queries
        .iter()
        .map(|q| (q.query_id.as_str(), q))
        .collect::<BTreeMap<_, _>>();
    if queries.len() != batch.accounting.len()
        || batch
            .accounting
            .keys()
            .any(|id| !queries.contains_key(id.as_str()))
    {
        return Err(bad("primer_payload_query_accounting_coverage"));
    }
    let mut reader = BufReader::new(VerifiedReader::open(
        &batch.paf,
        request.budgets.max_input_bytes,
    )?);
    let mut current = String::new();
    let mut alignments = Vec::new();
    let mut seen = BTreeSet::new();
    let mut lengths = cat.lengths.clone();
    let mut rows = Vec::new();
    let mut count = 0usize;
    let mut no_hit = false;
    let mut row_bytes = 0usize;
    let mut group_bytes = 0usize;
    let mut first_line = true;
    let finish = |id: &str,
                  alignments: &mut Vec<paf::Alignment>,
                  no_hit: bool,
                  rows: &mut Vec<PrimerPayloadReadEvidence>,
                  work: &mut usize,
                  row_bytes: &mut usize|
     -> Result<(), EngineError> {
        let accounting = batch
            .accounting
            .get(id)
            .ok_or_else(|| bad("primer_payload_query_accounting_missing"))?;
        if (alignments.is_empty() && *accounting == PrimerPayloadQueryAccounting::Mapped)
            || (!alignments.is_empty() && *accounting != PrimerPayloadQueryAccounting::Mapped)
            || (no_hit
                && (!alignments.is_empty()
                    || *accounting != PrimerPayloadQueryAccounting::ProducerReportedUnmapped))
        {
            return Err(bad("primer_payload_query_accounting_conflict"));
        }
        let row = evidence::analyze(
            queries[id],
            std::mem::take(alignments),
            batch,
            reference,
            cat,
            request,
            work,
        )?;
        *row_bytes = row_bytes
            .checked_add(encode(&row)?.len())
            .filter(|n| *n as u64 <= request.budgets.max_output_bytes.min(128_000_000))
            .ok_or_else(|| bad("primer_payload_evidence_memory_budget"))?;
        rows.push(row);
        Ok(())
    };
    while let Some(text) = line(&mut reader, request.budgets.max_line_bytes)? {
        if first_line && ctl.stop(true) {
            return Err(bad(INTERRUPTED));
        }
        first_line = false;
        ctl.check()?;
        let id = text.split('\t').next().unwrap_or("");
        let query = queries
            .get(id)
            .ok_or_else(|| bad("primer_payload_foreign_query"))?;
        if id != current {
            if !current.is_empty() {
                finish(
                    &current,
                    &mut alignments,
                    no_hit,
                    &mut rows,
                    work,
                    &mut row_bytes,
                )?;
            }
            if !seen.insert(id.to_string()) {
                return Err(bad("primer_payload_paf_requires_query_grouping"));
            }
            current = id.into();
            count = 0;
            group_bytes = 0;
            no_hit = false;
        }
        count += 1;
        if count > request.budgets.max_alignments_per_query {
            return Err(bad("primer_payload_alignment_count_budget"));
        }
        match paf::row(
            &text,
            query,
            request.policy.strict_alignment_detail,
            &request.budgets,
        )? {
            None => {
                if no_hit {
                    return Err(bad("primer_payload_duplicate_no_hit"));
                }
                no_hit = true;
            }
            Some((mut alignment, length)) => {
                if lengths
                    .insert(alignment.segment.target_name.clone(), length)
                    .is_some_and(|old| old != length)
                {
                    return Err(bad("primer_payload_reference_target_length_conflict"));
                }
                reference::annotate(&mut alignment, reference, cat, work, &request.budgets)?;
                let estimate = encode(&alignment.segment)?
                    .len()
                    .checked_add(
                        alignment
                            .blocks
                            .len()
                            .saturating_mul(std::mem::size_of::<paf::Block>()),
                    )
                    .and_then(|n| {
                        n.checked_add(
                            alignment
                                .skips
                                .len()
                                .saturating_mul(2 * std::mem::size_of::<usize>()),
                        )
                    })
                    .ok_or_else(|| bad("primer_payload_group_memory_budget"))?;
                group_bytes = group_bytes
                    .checked_add(estimate)
                    .filter(|n| *n as u64 <= request.budgets.max_output_bytes.min(128_000_000))
                    .ok_or_else(|| bad("primer_payload_group_memory_budget"))?;
                alignments.push(alignment);
            }
        }
    }
    reader.into_inner().finish()?;
    if !current.is_empty() {
        finish(
            &current,
            &mut alignments,
            no_hit,
            &mut rows,
            work,
            &mut row_bytes,
        )?;
    }
    for id in queries.keys().filter(|id| !seen.contains(**id)) {
        finish(id, &mut Vec::new(), false, &mut rows, work, &mut row_bytes)?;
    }
    rows.sort_by(|a, z| a.query_id.cmp(&z.query_id));
    Ok((rows, lengths))
}

fn table<T: Serialize>(
    path: &Path,
    values: &[T],
    header: &str,
) -> Result<RnaReadPanelArtifact, EngineError> {
    let writer = HashedWriter::create(path)?;
    let mut csv = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);
    csv.write_record([header]).map_err(|e| bad(e.to_string()))?;
    for value in values {
        csv.write_record([serde_json::to_string(value).map_err(|e| bad(e.to_string()))?])
            .map_err(|e| bad(e.to_string()))?;
    }
    csv.flush().map_err(|e| bad(e.to_string()))?;
    csv.into_inner()
        .map_err(|e| bad(e.to_string()))?
        .receipt(path.file_name().unwrap().to_string_lossy().into_owned())
}

impl GentleEngine {
    /// One explicit batch per reference, verified at EOF; no aligner or reference discovery.
    pub fn summarize_primer_payload_alignments(
        &self,
        request: SummarizePrimerPayloadAlignmentsRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<PrimerPayloadAlignmentReport, EngineError> {
        admission(&request)?;
        let identity = hash(&request)?;
        let mut ctl = Control::new(
            request.cancel_path.clone(),
            identity.clone(),
            request.batches.len(),
            progress,
        )?;
        let trace = prepare::trace_source(&request.trace, &request.budgets)?;
        let additional_input_bytes = request
            .references
            .iter()
            .flat_map(|r| [r.fasta.bytes, r.annotation.bytes])
            .chain(request.batches.iter().map(|b| b.paf.bytes))
            .try_fold(0u64, |n, x| n.checked_add(x))
            .ok_or_else(|| bad("primer_payload_input_budget"))?;
        let prepared = binding::queries(
            &request.queries,
            &trace,
            &request.trace,
            &request.budgets,
            additional_input_bytes,
        )?;
        if prepared.content_sha256 != request.query_manifest_sha256
            || encode(&prepared.request.trace)? != encode(&request.trace)?
        {
            return Err(bad("primer_payload_query_manifest_source_binding"));
        }
        let mut input_bytes = request.queries.bytes;
        for n in request
            .references
            .iter()
            .flat_map(|r| [r.fasta.bytes, r.annotation.bytes])
            .chain(request.batches.iter().map(|b| b.paf.bytes))
            .chain(prepared.artifacts.iter().map(|a| a.bytes))
            .chain(trace.artifacts.iter().map(|a| a.bytes))
            .chain(
                trace
                    .files_completed
                    .iter()
                    .filter(|u| prepared.files_completed.contains(&u.file_id))
                    .flat_map(|u| u.artifacts.iter().map(|a| a.bytes)),
            )
        {
            input_bytes = input_bytes
                .checked_add(n)
                .filter(|n| *n <= request.budgets.max_input_bytes)
                .ok_or_else(|| bad("primer_payload_input_budget"))?;
        }
        let mut catalogs = BTreeMap::new();
        for r in &request.references {
            catalogs.insert(
                r.reference_id.clone(),
                reference::catalog(r, &request.budgets)?,
            );
        }
        crosswalk(&request, &trace.primer_panel, &catalogs)?;
        let stage = staging(&request.output_dir)?;
        let mut all_rows = Vec::new();
        let mut completed = Vec::new();
        let mut work = 0usize;
        let mut memory = 0usize;
        let mut batches = request.batches.iter().collect::<Vec<_>>();
        batches.sort_by(|a, z| a.batch_id.cmp(&z.batch_id));
        for batch in batches {
            if ctl.stop(true) {
                break;
            }
            let reference = request
                .references
                .iter()
                .find(|r| r.reference_id == batch.reference_id)
                .expect("admitted reference");
            let cat = &catalogs[&reference.reference_id];
            let result = (|| {
                let (mut rows, lengths) = batch_rows(
                    batch, reference, cat, &prepared, &request, &mut ctl, &mut work,
                )?;
                let ranges =
                    evidence::context_ranges(&rows, &trace.primer_panel, &lengths, &request.policy);
                let contexts = reference::verify_fasta(
                    reference,
                    &lengths,
                    &ranges,
                    &request.budgets,
                    &mut ctl,
                )?;
                evidence::seeds(&mut rows, &trace.primer_panel, &contexts, &request.policy)?;
                Ok::<_, EngineError>(rows)
            })();
            let rows = match result {
                Ok(rows) => rows,
                Err(e) if e.message == INTERRUPTED => break,
                Err(e) => return Err(e),
            };
            if all_rows.len() + rows.len() > request.budgets.max_records {
                return Err(bad("primer_payload_evidence_row_budget"));
            }
            for row in &rows {
                memory = memory
                    .checked_add(encode(row)?.len())
                    .filter(|n| *n as u64 <= request.budgets.max_output_bytes.min(128_000_000))
                    .ok_or_else(|| bad("primer_payload_evidence_memory_budget"))?;
            }
            all_rows.extend(rows);
            completed.push(batch.batch_id.clone());
            ctl.completed = completed.len();
        }
        ctl.stop(true);
        let refs = request
            .batches
            .iter()
            .filter(|b| completed.contains(&b.batch_id))
            .map(|b| b.reference_id.clone())
            .collect::<Vec<_>>();
        let (profiles, clusters) = evidence::aggregate(
            &all_rows,
            &prepared,
            &trace,
            &refs,
            request.policy.breakpoint_tolerance_bases,
        );
        let artifacts = vec![
            table(
                &stage.path().join("read_anchor_evidence.tsv"),
                &all_rows,
                "evidence_json",
            )?,
            table(
                &stage.path().join("partner_profiles.tsv"),
                &profiles,
                "partner_profile_json",
            )?,
            table(
                &stage.path().join("breakpoint_clusters.tsv"),
                &clusters,
                "breakpoint_cluster_json",
            )?,
        ];
        let mut report = PrimerPayloadAlignmentReport {
            schema: PRIMER_PAYLOAD_REPORT_SCHEMA.into(),
            content_sha256: String::new(),
            request_sha256: identity,
            request: request.clone(),
            status: ctl.status().into(),
            control_error: ctl.error.clone(),
            cohort: prepared.request.cohort,
            source_trace_status: trace.status,
            source_inputs_missing: trace.inputs_not_completed.len(),
            source_query_status: prepared.status,
            source_query_units_missing: prepared.files_not_completed.len(),
            source_excluded_anchor_count: prepared.unassessed.len(),
            batches_not_completed: request
                .batches
                .iter()
                .filter(|b| !completed.contains(&b.batch_id))
                .map(|b| b.batch_id.clone())
                .collect(),
            batches_completed: completed,
            evidence_rows: all_rows.len(),
            profiles,
            breakpoint_clusters: clusters,
            artifacts,
            nonclaims: nonclaims(),
        };
        report.nonclaims.push("Alignment detail verifies coordinate consumption, not independently reproduced base matching or search completeness. The source query-preparation witness is content-bound, not authenticated.".into());
        report.nonclaims.push("Profiles count distinct reads per key; alternative partner groups, anchors and reference scopes are not additive. Clusters use sorted maximum-span tolerance and the minimum observed edge representative.".into());
        report.content_sha256 = hash(&report)?;
        publish(stage, &request.output_dir, &report, &request.budgets)?;
        Ok(report)
    }
}

pub(super) fn validate_report(report: &PrimerPayloadAlignmentReport) -> Result<(), EngineError> {
    admission(&report.request)?;
    let mut basis = report.clone();
    basis.content_sha256.clear();
    let done = report.batches_completed.iter().collect::<BTreeSet<_>>();
    let missing = report.batches_not_completed.iter().collect::<BTreeSet<_>>();
    let ids = report
        .request
        .batches
        .iter()
        .map(|b| &b.batch_id)
        .collect::<BTreeSet<_>>();
    if report.schema != PRIMER_PAYLOAD_REPORT_SCHEMA
        || hash(&basis)? != report.content_sha256
        || report.request_sha256 != hash(&report.request)?
        || !binding::outcome(&report.status, &report.control_error)
        || done.len() != report.batches_completed.len()
        || missing.len() != report.batches_not_completed.len()
        || !done.is_disjoint(&missing)
        || done.union(&missing).copied().collect::<BTreeSet<_>>() != ids
        || (report.status == "completed" && !missing.is_empty())
        || report
            .artifacts
            .iter()
            .map(|a| a.path.as_str())
            .collect::<BTreeSet<_>>()
            != BTreeSet::from([
                "read_anchor_evidence.tsv",
                "partner_profiles.tsv",
                "breakpoint_clusters.tsv",
            ])
        || report.artifacts.len() != 3
    {
        return Err(bad("primer_payload_report_binding"));
    }
    Ok(())
}
