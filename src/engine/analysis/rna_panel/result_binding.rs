//! Shared submitted/effective-request and published-result bindings; no live project mutation.

use super::*;

const SCHEMA: &str = "gentle.rna_read_panel_result_binding.v1";

fn digest<T: Serialize>(value: &T) -> Result<String, EngineError> {
    Ok(sha256_prefixed_bytes(&encode(value)?))
}

fn validate_summary(summary: &RnaReadPanelSummary) -> Result<(), EngineError> {
    let mut basis = summary.clone();
    basis.content_sha256.clear();
    if summary.schema != RNA_PANEL_SUMMARY_SCHEMA
        || digest(&basis)? != summary.content_sha256
        || !matches!(summary.completeness.as_str(), "complete" | "partial")
    {
        return Err(bad("panel_result_summary_binding"));
    }
    Ok(())
}

impl GentleEngine {
    /// Bind a screen returned by the shared engine, before any host cancellation/stale gate.
    /// Only MD5 resolution/normalization may differ from the submitted request.
    pub fn rna_read_panel_screen_result_binding(
        submitted: &ScreenRnaReadPanelRequest,
        report: &RnaReadPanelScreenReport,
    ) -> Result<RnaReadPanelResultBinding, EngineError> {
        validate_screen_report_binding(report)?;
        let submitted_hash = digest(submitted)?;
        if report.submitted_request_sha256.as_deref() != Some(submitted_hash.as_str())
            || submitted.files.len() != report.request.files.len()
            || submitted
                .expected_panel_sha256
                .as_ref()
                .is_some_and(|h| h != &report.panel_sha256)
        {
            return Err(bad("panel_result_submitted_request_binding"));
        }
        let mut effective = submitted.clone();
        for (input, resolved) in effective.files.iter_mut().zip(&report.request.files) {
            if input.expected_md5.as_ref().is_some_and(|h| {
                resolved.expected_md5.as_deref() != Some(h.to_ascii_lowercase().as_str())
            }) {
                return Err(bad("panel_result_explicit_checksum_binding"));
            }
            input.expected_md5 = resolved.expected_md5.clone();
        }
        if encode(&effective)? != encode(&report.request)? {
            return Err(bad("panel_result_effective_request_binding"));
        }
        Ok(RnaReadPanelResultBinding {
            schema: SCHEMA.into(),
            kind: RnaReadPanelResultKind::Screen,
            submitted_request_sha256: submitted_hash,
            effective_request_sha256: digest(&report.request)?,
            panel_sha256: report.panel_sha256.clone(),
            report_content_sha256: report.content_sha256.clone(),
            report_path: Path::new(&submitted.output_dir)
                .join("report.json")
                .display()
                .to_string(),
            outcome: report.status.clone(),
        })
    }

    /// Summary requests are not normalized; source-report hashes remain in the summary itself.
    pub fn rna_read_panel_summary_result_binding(
        submitted: &ExportRnaReadPanelSummaryRequest,
        summary: &RnaReadPanelSummary,
    ) -> Result<RnaReadPanelResultBinding, EngineError> {
        validate_summary(summary)?;
        let hash = digest(submitted)?;
        let mut requested = submitted.report_paths.clone();
        let mut retained = summary
            .source_runs
            .iter()
            .map(|r| r.report_path.clone())
            .collect::<Vec<_>>();
        requested.sort();
        retained.sort();
        if summary.submitted_request_sha256.as_deref() != Some(hash.as_str())
            || requested != retained
        {
            return Err(bad("panel_result_summary_request_binding"));
        }
        Ok(RnaReadPanelResultBinding {
            schema: SCHEMA.into(),
            kind: RnaReadPanelResultKind::Summary,
            submitted_request_sha256: hash.clone(),
            effective_request_sha256: hash,
            panel_sha256: summary.panel_sha256.clone(),
            report_content_sha256: summary.content_sha256.clone(),
            report_path: Path::new(&submitted.output_dir)
                .join("report.json")
                .display()
                .to_string(),
            outcome: summary.completeness.clone(),
        })
    }

    /// Check an already inspected result against a retained receipt, without rereading a manifest,
    /// original reads, or live source records. A filename alone never establishes this relation.
    pub fn validate_rna_read_panel_result_binding(
        binding: &RnaReadPanelResultBinding,
        path: &str,
        report: Option<&RnaReadPanelScreenReport>,
        summary: Option<&RnaReadPanelSummary>,
    ) -> Result<(), EngineError> {
        if binding.schema != SCHEMA || binding.report_path != path {
            return Err(bad("panel_result_location_binding"));
        }
        let agrees = match (binding.kind, report, summary) {
            (RnaReadPanelResultKind::Screen, Some(report), None) => {
                validate_screen_report_binding(report)?;
                report.submitted_request_sha256.as_ref() == Some(&binding.submitted_request_sha256)
                    && digest(&report.request)? == binding.effective_request_sha256
                    && report.panel_sha256 == binding.panel_sha256
                    && report.content_sha256 == binding.report_content_sha256
                    && report.status == binding.outcome
            }
            (RnaReadPanelResultKind::Summary, None, Some(summary)) => {
                validate_summary(summary)?;
                summary.submitted_request_sha256.as_ref() == Some(&binding.submitted_request_sha256)
                    && binding.submitted_request_sha256 == binding.effective_request_sha256
                    && summary.panel_sha256 == binding.panel_sha256
                    && summary.content_sha256 == binding.report_content_sha256
                    && summary.completeness == binding.outcome
            }
            _ => false,
        };
        if !agrees {
            return Err(bad("panel_result_content_binding"));
        }
        Ok(())
    }
}
