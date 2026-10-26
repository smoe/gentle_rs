//! Named RNA artifact-only and independent observation/control policies, never inferred from output JSON.

use super::*;
use crate::engine::{Operation, RnaReadPanelResultBinding};

pub(super) enum Policy {
    Project,
    Screen(Box<crate::engine::ScreenRnaReadPanelRequest>),
    Summary(Box<crate::engine::ExportRnaReadPanelSummaryRequest>),
    PrimerTrace(Box<crate::engine::TracePrimerReadsRequest>),
    Observe { expensive: bool },
    Control,
}

impl Policy {
    pub(super) fn lane(&self) -> Lane {
        match self {
            Self::Control => Lane::Control,
            Self::Observe { expensive: false } => Lane::Observation,
            _ => Lane::Execution,
        }
    }
    pub(super) fn needs_project(&self) -> bool {
        matches!(self, Self::Project)
    }
}

fn classify(op: Operation) -> Policy {
    match op {
        Operation::TracePrimerReads { request } => Policy::PrimerTrace(request),
        Operation::ValidateRnaPrimerPanel { .. } | Operation::InspectPrimerTrace { .. } => Policy::Observe { expensive: true },
        Operation::CancelPrimerTrace { .. } => Policy::Control,
        Operation::ScreenRnaReadPanel { request } => Policy::Screen(request),
        Operation::ExportRnaReadPanelSummary { request } => Policy::Summary(request),
        Operation::InspectRnaReadPanel {
            verify_artifacts, ..
        } => Policy::Observe {
            expensive: verify_artifacts,
        },
        Operation::CancelRnaReadPanel { .. } => Policy::Control,
        _ => Policy::Project,
    }
}

#[cfg(test)]
mod primer_trace_tests {
    use super::*;

    #[test]
    fn primer_trace_lifecycle_is_explicit_not_inferred_from_read_only() {
        let request:crate::engine::TracePrimerReadsRequest=serde_json::from_value(serde_json::json!({
            "schema":"gentle.rna_primer_trace_request.v1","primer_panel":{"kind":"file","path":"missing.json","expected_sha256":"sha256:synthetic"},
            "files":[],"output_dir":"synthetic-output"})).unwrap();
        let policy=classify(Operation::TracePrimerReads {request:Box::new(request.clone())});
        assert!(!policy.needs_project());assert!(matches!(policy,Policy::PrimerTrace(_)));
        let policy=classify(Operation::CancelPrimerTrace {request:Box::new(request)});
        assert!(!policy.needs_project());assert!(matches!(policy.lane(),Lane::Control));
        let policy=classify(Operation::InspectPrimerTrace {path:"missing.json".into(),expected_result:None,file_id:None,rows:0});
        assert!(!policy.needs_project());assert!(matches!(policy,Policy::Observe {expensive:true}));
        assert!(publication(policy,&serde_json::json!({})).is_ok());
    }
}

/// Admission examines only an already parsed inline payload, never reads an @file.
/// Native/agent callers pass the frozen shared-parser result through submit_parsed.
pub(super) fn admission(command: &ShellCommand) -> Policy {
    match command {
        ShellCommand::RnaReadPanelLoadDraft { .. } => Policy::Observe { expensive: false },
        ShellCommand::RnaReadPanelScreenDraft {
            output_path: None, ..
        } => Policy::Observe { expensive: false },
        ShellCommand::Op { payload } => serde_json::from_str::<Operation>(payload)
            .map(classify)
            .unwrap_or(Policy::Project),
        _ => Policy::Project,
    }
}

/// Resolve an Op payload once on the worker, then execute exactly those captured bytes.
pub(super) fn freeze(command: &mut ShellCommand) -> Result<Policy, String> {
    if matches!(
        command,
        ShellCommand::RnaReadPanelLoadDraft { .. }
            | ShellCommand::RnaReadPanelScreenDraft {
                output_path: None,
                ..
            }
    ) {
        return Ok(Policy::Observe { expensive: false });
    }
    let ShellCommand::Op { payload } = command else {
        return Ok(Policy::Project);
    };
    let op = crate::engine_shell::resolve_operation_payload(payload)?;
    *payload = serde_json::to_string(&op).map_err(|e| e.to_string())?;
    Ok(classify(op))
}

pub(super) fn publication(
    policy: Policy,
    output: &serde_json::Value,
) -> Result<Publication, String> {
    if let Policy::PrimerTrace(request) = &policy {
        let report = serde_json::from_value(output.pointer("/result/rna_primer_trace")
            .ok_or("Primer trace returned no typed report")?.clone()).map_err(|e| e.to_string())?;
        let binding = GentleEngine::primer_trace_result_binding(request, &report).map_err(|e| e.to_string())?;
        return Ok(Publication::PrimerTrace(binding));
    }
    let binding: RnaReadPanelResultBinding = match policy {
        Policy::Project => return Ok(Publication::Project),
        Policy::Observe { .. } | Policy::Control => return Ok(Publication::Observation),
        Policy::PrimerTrace(_) => unreachable!("handled above"),
        Policy::Screen(request) => {
            let report = serde_json::from_value(
                output
                    .pointer("/result/rna_read_panel_screen")
                    .ok_or("Screen returned no typed report")?
                    .clone(),
            )
            .map_err(|e| e.to_string())?;
            GentleEngine::rna_read_panel_screen_result_binding(&request, &report)
                .map_err(|e| e.to_string())?
        }
        Policy::Summary(request) => {
            let report = serde_json::from_value(
                output
                    .pointer("/result/rna_read_panel_summary")
                    .ok_or("Summary export returned no typed report")?
                    .clone(),
            )
            .map_err(|e| e.to_string())?;
            GentleEngine::rna_read_panel_summary_result_binding(&request, &report)
                .map_err(|e| e.to_string())?
        }
    };
    Ok(Publication::Artifact(binding))
}
