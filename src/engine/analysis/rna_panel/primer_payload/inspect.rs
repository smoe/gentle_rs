//! Bounded summaries and explicit sidecar hash verification, never implicit source replay.

use super::*;
use serde_json::{Value, json};

impl GentleEngine {
    pub fn primer_payload_alignment_result_binding(
        request: &SummarizePrimerPayloadAlignmentsRequest,
        report: &PrimerPayloadAlignmentReport,
    ) -> Result<PrimerPayloadResultBinding, EngineError> {
        summarize::validate_report(report)?;
        if hash(request)? != report.request_sha256 {
            return Err(bad("primer_payload_request_binding"));
        }
        Ok(PrimerPayloadResultBinding {
            schema: PRIMER_PAYLOAD_BINDING_SCHEMA.into(),
            kind: "alignments".into(),
            request_sha256: report.request_sha256.clone(),
            report_content_sha256: report.content_sha256.clone(),
            report_path: Path::new(&request.output_dir)
                .join("report.json")
                .display()
                .to_string(),
            outcome: report.status.clone(),
        })
    }

    /// The bounded summary does not echo a potentially large query/accounting population to agents.
    pub fn inspect_primer_payload(
        path: &str,
        expected: Option<&PrimerPayloadResultBinding>,
        rows: usize,
        verify_artifacts: bool,
    ) -> Result<PrimerPayloadInspection, EngineError> {
        if rows > 1000 {
            return Err(bad("primer_payload_inspection_row_budget"));
        }
        regular(Path::new(path))?;
        let value: Value = load(path)?;
        let (binding, artifacts, b, summary, mut selected) = match value
            .get("schema")
            .and_then(Value::as_str)
        {
            Some(PRIMER_PAYLOAD_QUERIES_SCHEMA) => {
                let report: PrimerPayloadQueries =
                    serde_json::from_value(value).map_err(|e| bad(e.to_string()))?;
                let binding = Self::primer_payload_query_result_binding(&report.request, &report)?;
                let summary = json!({"schema":report.schema,"content_sha256":report.content_sha256,"request_sha256":report.request_sha256,"status":report.status,"control_error":report.control_error,
                    "cohort":report.request.cohort,"source_trace_status":report.source_trace_status,"source_inputs_missing":report.source_inputs_missing,"source_mode":report.request.source_mode,"anchor_selection":report.request.anchor_selection,
                    "files_completed":report.files_completed,"files_not_completed":report.files_not_completed,"query_count":report.queries.len(),"excluded_anchor_count":report.unassessed.len(),"artifacts":report.artifacts,"nonclaims":report.nonclaims});
                let selected = report
                    .queries
                    .iter()
                    .take(rows)
                    .map(|q| serde_json::to_value(q).map_err(|e| bad(e.to_string())))
                    .collect::<Result<Vec<_>, _>>()?;
                (
                    binding,
                    report.artifacts,
                    report.request.budgets,
                    summary,
                    selected,
                )
            }
            Some(PRIMER_PAYLOAD_REPORT_SCHEMA) => {
                let report: PrimerPayloadAlignmentReport =
                    serde_json::from_value(value).map_err(|e| bad(e.to_string()))?;
                let binding =
                    Self::primer_payload_alignment_result_binding(&report.request, &report)?;
                let summary = json!({"schema":report.schema,"content_sha256":report.content_sha256,"request_sha256":report.request_sha256,"status":report.status,"control_error":report.control_error,
                    "cohort":report.cohort,"source_trace_status":report.source_trace_status,"source_inputs_missing":report.source_inputs_missing,"batches_completed":report.batches_completed,"batches_not_completed":report.batches_not_completed,
                    "source_query_status":report.source_query_status,"source_query_units_missing":report.source_query_units_missing,"source_excluded_anchor_count":report.source_excluded_anchor_count,
                    "evidence_rows":report.evidence_rows,"profile_count":report.profiles.len(),"cluster_count":report.breakpoint_clusters.len(),"artifacts":report.artifacts,"nonclaims":report.nonclaims});
                (
                    binding,
                    report.artifacts,
                    report.request.budgets,
                    summary,
                    Vec::new(),
                )
            }
            _ => return Err(bad("primer_payload_inspection_schema")),
        };
        let actual = fs::canonicalize(path).map_err(|e| bad(e.to_string()))?;
        if fs::canonicalize(&binding.report_path).ok().as_ref() != Some(&actual)
            || expected.is_some_and(|e| encode(e).ok() != encode(&binding).ok())
        {
            return Err(bad("primer_payload_inspection_result_binding"));
        }
        let root = Path::new(path)
            .parent()
            .unwrap_or(Path::new("."))
            .join("artifacts");
        if verify_artifacts {
            let mut total = 0u64;
            for artifact in &artifacts {
                total = total
                    .checked_add(artifact.bytes)
                    .filter(|n| *n <= b.max_output_bytes)
                    .ok_or_else(|| bad("primer_payload_inspection_byte_budget"))?;
                let mut input =
                    VerifiedReader::open(&artifact_file(&root, artifact)?, b.max_output_bytes)?;
                std::io::copy(&mut input, &mut std::io::sink()).map_err(|e| bad(e.to_string()))?;
                input.finish()?;
            }
        }
        if binding.kind == "alignments" && rows > 0 {
            let artifact = artifacts
                .iter()
                .find(|a| a.path == "read_anchor_evidence.tsv")
                .expect("validated artifact");
            let file = artifact_file(&root, artifact)?;
            regular(Path::new(&file.path))?;
            let input = fs::File::open(&file.path).map_err(|e| bad(e.to_string()))?;
            let mut csv = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .from_reader(input.take(8_000_001));
            if csv
                .headers()
                .map_err(|e| bad(e.to_string()))?
                .iter()
                .collect::<Vec<_>>()
                != vec!["evidence_json"]
            {
                return Err(bad("primer_payload_inspection_columns"));
            }
            for record in csv.records().take(rows) {
                let record = record.map_err(|e| bad(e.to_string()))?;
                if record.len() != 1 || record[0].len() > b.max_line_bytes {
                    return Err(bad("primer_payload_inspection_row_budget"));
                }
                let row: PrimerPayloadReadEvidence =
                    serde_json::from_str(&record[0]).map_err(|e| bad(e.to_string()))?;
                selected.push(serde_json::to_value(row).map_err(|e| bad(e.to_string()))?);
            }
        }
        if encode(&selected)?.len() > 8_000_000 {
            return Err(bad("primer_payload_inspection_response_budget"));
        }
        Ok(PrimerPayloadInspection {
            report: summary,
            rows: selected,
            artifact_content_verified: verify_artifacts,
        })
    }

    pub fn cancel_primer_payload(
        request: &PrimerPayloadControlRequest,
    ) -> Result<String, EngineError> {
        let (identity, path) = match request {
            PrimerPayloadControlRequest::Queries(r)
                if r.schema == PRIMER_PAYLOAD_PREPARE_SCHEMA =>
            {
                (hash(r)?, r.cancel_path.as_deref())
            }
            PrimerPayloadControlRequest::Alignments(r)
                if r.schema == PRIMER_PAYLOAD_REQUEST_SCHEMA =>
            {
                (hash(r)?, r.cancel_path.as_deref())
            }
            _ => return Err(bad("primer_payload_cancel_schema")),
        };
        let path = path.ok_or_else(|| bad("primer_payload_cancel_path_required"))?;
        if control::marker_requested(Some(path), std::slice::from_ref(&identity))? {
            return Ok("Cancellation already requested; termination not confirmed".into());
        }
        let parent = Path::new(path)
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let mut marker = tempfile::NamedTempFile::new_in(parent).map_err(|e| bad(e.to_string()))?;
        write!(marker, "{identity}").map_err(|e| bad(e.to_string()))?;
        marker
            .as_file()
            .sync_all()
            .map_err(|e| bad(e.to_string()))?;
        marker
            .persist_noclobber(path)
            .map_err(|e| bad(e.to_string()))?;
        Ok("Cancellation requested for this exact request; termination not confirmed".into())
    }
}
