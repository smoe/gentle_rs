//! Named RNA artifact-only and independent observation/control policies, never inferred from output JSON.

use super::*;
use crate::engine::{Operation, RnaReadPanelResultBinding};

pub(super) enum Policy {
    Project,
    Screen(Box<crate::engine::ScreenRnaReadPanelRequest>),
    Summary(Box<crate::engine::ExportRnaReadPanelSummaryRequest>),
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

/// Admission examines only an already parsed inline payload, never reads an @file.
/// Native/agent callers pass the frozen shared-parser result through submit_parsed.
pub(super) fn admission(command: &ShellCommand) -> Policy {
    match command {
        ShellCommand::Op { payload } => serde_json::from_str::<Operation>(payload)
            .map(classify)
            .unwrap_or(Policy::Project),
        _ => Policy::Project,
    }
}

/// Resolve an Op payload once on the worker, then execute exactly those captured bytes.
pub(super) fn freeze(command: &mut ShellCommand) -> Result<Policy, String> {
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
    let binding: RnaReadPanelResultBinding = match policy {
        Policy::Project => return Ok(Publication::Project),
        Policy::Observe { .. } | Policy::Control => return Ok(Publication::Observation),
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
