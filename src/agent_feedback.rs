//! Privacy-filtered execution evidence, independent of project persistence and approvals.

use crate::{digest_utils::sha256_prefixed_str, engine::GentleEngine};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::{
    collections::BTreeSet,
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};

pub const AGENT_EXECUTION_FEEDBACK_SCHEMA: &str = "gentle.agent_execution_feedback.v1";
pub const AGENT_EXECUTION_RECEIPT_LIMIT: usize = 100;

pub(crate) fn agent_context_id_is_valid(value: &str) -> bool {
    value
        .strip_prefix("sha256:")
        .is_some_and(|hex| hex.len() == 64 && hex.bytes().all(|b| b.is_ascii_hexdigit()))
}

/// Process/session identity, not a scientific content fingerprint or an approval.
pub fn new_agent_context_id() -> String {
    static NEXT: AtomicU64 = AtomicU64::new(1);
    sha256_prefixed_str(&format!(
        "{}:{}:{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum AgentExecutionStatus {
    #[default]
    NotRun,
    Blocked,
    Failed,
    Partial,
    /// The dispatcher returned, but the requested UI/job effect is not yet verified.
    Dispatched,
    Running,
    Cancelled,
    /// Command completion only; never an assertion that scientific gates passed.
    Completed,
}

impl AgentExecutionStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NotRun => "not_run",
            Self::Blocked => "blocked",
            Self::Failed => "failed",
            Self::Partial => "partial",
            Self::Dispatched => "dispatched",
            Self::Running => "running",
            Self::Cancelled => "cancelled",
            Self::Completed => "completed",
        }
    }

    /// Interpret only engine-owned job envelopes, never free-text messages.
    pub fn from_shell_output(output: &Value) -> Self {
        match output.get("schema").and_then(Value::as_str) {
            Some(
                "gentle.blast_async_start.v1"
                | "gentle.blast_async_start.v2"
                | "gentle.blast_async_status.v1"
                | "gentle.blast_async_cancel.v1",
            ) => match output.pointer("/job/state").and_then(Value::as_str) {
                Some("completed")
                    if output
                        .pointer("/job/result_available")
                        .and_then(Value::as_bool)
                        == Some(true) =>
                {
                    Self::Completed
                }
                Some("failed") => Self::Failed,
                Some("cancelled") => Self::Cancelled,
                Some("running" | "queued") => Self::Running,
                _ => Self::Dispatched,
            },
            _ if output.get("applied").and_then(Value::as_bool) == Some(false)
                && output.get("ui_intent").is_some() =>
            {
                Self::Dispatched
            }
            _ => Self::Completed,
        }
    }

    pub fn after_sequence_import(self, imported: bool) -> Self {
        if imported { self } else { Self::Partial }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct AgentExecutionRevision {
    pub execution: u64,
    pub mutation: u64,
    pub structural: u64,
}

impl AgentExecutionRevision {
    pub fn capture(engine: &GentleEngine) -> Self {
        Self {
            execution: engine.execution_revision(),
            mutation: engine.mutation_revision(),
            structural: engine.structural_revision(),
        }
    }
}

/// No raw command, error, result body, path, sequence, or model-authored text.
/// Exact bytes remain in the caller's local execution log; hashes correlate them.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct AgentExecutionReceipt {
    pub receipt_id: String,
    pub session_id: String,
    pub turn_id: Option<String>,
    pub suggestion_index: Option<usize>,
    pub command_sha256: String,
    pub output_sha256: Option<String>,
    pub job_id_sha256: Option<String>,
    pub error_sha256: Option<String>,
    pub status: AgentExecutionStatus,
    pub before: Option<AgentExecutionRevision>,
    pub after: Option<AgentExecutionRevision>,
}

impl AgentExecutionReceipt {
    pub fn new(
        command: &str,
        status: AgentExecutionStatus,
        before: Option<AgentExecutionRevision>,
        after: Option<AgentExecutionRevision>,
    ) -> Self {
        Self {
            receipt_id: new_agent_context_id(),
            session_id: String::new(),
            turn_id: None,
            suggestion_index: None,
            command_sha256: sha256_prefixed_str(command),
            output_sha256: None,
            job_id_sha256: None,
            error_sha256: None,
            status,
            before,
            after,
        }
    }

    pub fn bind_output(&mut self, output: &Value) {
        self.job_id_sha256 = output
            .pointer("/job/job_id")
            .and_then(Value::as_str)
            .map(sha256_prefixed_str);
        self.output_sha256 = serde_json::to_string(output)
            .ok()
            .map(|text| sha256_prefixed_str(&text));
    }

    pub fn bind_error(&mut self, error: &str) {
        self.error_sha256 = Some(sha256_prefixed_str(error));
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum AgentFeedbackApplicability {
    SameStructuralRevision,
    RecheckRequired,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentExecutionFeedbackRow {
    pub receipt: AgentExecutionReceipt,
    pub applicability: AgentFeedbackApplicability,
}

/// Session-only projection. Structural equality does not revalidate external files or QA.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentExecutionFeedback {
    pub schema: String,
    pub session_id: String,
    pub current_revision: Option<AgentExecutionRevision>,
    pub omitted_receipt_count: usize,
    pub rows: Vec<AgentExecutionFeedbackRow>,
}

impl AgentExecutionFeedback {
    /// Validate the bounded path-free wire projection before provider transport.
    pub fn validate(&self) -> Result<(), String> {
        if self.schema != AGENT_EXECUTION_FEEDBACK_SCHEMA
            || !agent_context_id_is_valid(&self.session_id)
            || self.rows.len() > AGENT_EXECUTION_RECEIPT_LIMIT
        {
            return Err(
                "Invalid execution feedback schema, session identity or receipt limit".to_string(),
            );
        }
        for row in &self.rows {
            let r = &row.receipt;
            if r.session_id != self.session_id
                || [&r.receipt_id, &r.command_sha256]
                    .into_iter()
                    .any(|id| !agent_context_id_is_valid(id))
                || [
                    &r.turn_id,
                    &r.output_sha256,
                    &r.error_sha256,
                    &r.job_id_sha256,
                ]
                .into_iter()
                .flatten()
                .any(|id| !agent_context_id_is_valid(id))
            {
                return Err(
                    "Execution feedback identities must be scoped SHA-256 values".to_string(),
                );
            }
            let expected = if r
                .after
                .zip(self.current_revision)
                .is_some_and(|(after, now)| after.structural == now.structural)
            {
                AgentFeedbackApplicability::SameStructuralRevision
            } else {
                AgentFeedbackApplicability::RecheckRequired
            };
            if row.applicability != expected {
                return Err(
                    "Execution feedback applicability disagrees with its revision binding"
                        .to_string(),
                );
            }
        }
        Ok(())
    }

    /// Project only receipts attributable to a suggestion in the visible conversation.
    /// Unattributed local commands remain in the host log, not in model context.
    pub fn project(
        session_id: &str,
        current: Option<AgentExecutionRevision>,
        receipts: &[AgentExecutionReceipt],
        turn_ids: &BTreeSet<String>,
    ) -> Self {
        let mut rows = receipts
            .iter()
            .rev()
            .filter(|receipt| receipt.session_id == session_id)
            .filter(|receipt| {
                receipt
                    .turn_id
                    .as_ref()
                    .is_some_and(|id| turn_ids.contains(id))
                    && receipt.suggestion_index.is_some_and(|index| index > 0)
            })
            .take(AGENT_EXECUTION_RECEIPT_LIMIT)
            .map(|receipt| AgentExecutionFeedbackRow {
                receipt: receipt.clone(),
                applicability: if receipt
                    .after
                    .zip(current)
                    .is_some_and(|(after, now)| after.structural == now.structural)
                {
                    AgentFeedbackApplicability::SameStructuralRevision
                } else {
                    AgentFeedbackApplicability::RecheckRequired
                },
            })
            .collect::<Vec<_>>();
        rows.reverse();
        Self {
            schema: AGENT_EXECUTION_FEEDBACK_SCHEMA.to_string(),
            session_id: session_id.to_string(),
            current_revision: current,
            omitted_receipt_count: receipts.len().saturating_sub(rows.len()),
            rows,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn agent_feedback_async_start_is_not_completion() {
        for schema in ["gentle.blast_async_start.v1", "gentle.blast_async_start.v2"] {
            for state in ["queued", "running", "future_state"] {
                let output = serde_json::json!({"schema":schema, "job":{"state":state}});
                assert_ne!(
                    AgentExecutionStatus::from_shell_output(&output),
                    AgentExecutionStatus::Completed
                );
            }
        }
        for (state, expected) in [
            ("completed", AgentExecutionStatus::Completed),
            ("failed", AgentExecutionStatus::Failed),
            ("cancelled", AgentExecutionStatus::Cancelled),
        ] {
            let output = serde_json::json!({"schema":"gentle.blast_async_status.v1", "job":{"state":state,"result_available":true}});
            assert_eq!(AgentExecutionStatus::from_shell_output(&output), expected);
        }
        assert_eq!(
            AgentExecutionStatus::from_shell_output(
                &serde_json::json!({"ui_intent":{},"applied":false})
            ),
            AgentExecutionStatus::Dispatched
        );
    }

    #[test]
    fn agent_feedback_is_bounded_private_and_session_scoped() {
        let mut receipt = AgentExecutionReceipt::new(
            "import /private/secret.fa --token PASSWORD ACGTACGT",
            AgentExecutionStatus::Completed,
            None,
            None,
        );
        let session_id = new_agent_context_id();
        receipt.session_id = session_id.clone();
        let turn_id = new_agent_context_id();
        receipt.turn_id = Some(turn_id.clone());
        receipt.suggestion_index = Some(1);
        receipt
            .bind_output(&serde_json::json!({"sequence":"ACGTACGT", "path":"/private/secret.fa"}));
        receipt.bind_error("PASSWORD");
        let mut receipts = vec![receipt; 103];
        receipts[0].session_id = "previous_project".into();
        receipts[1].turn_id = Some("expired_turn".into());
        let feedback = AgentExecutionFeedback::project(
            &session_id,
            None,
            &receipts,
            &BTreeSet::from([turn_id]),
        );
        feedback.validate().expect("valid private projection");
        assert_eq!(feedback.rows.len(), 100);
        assert_eq!(feedback.omitted_receipt_count, 3);
        let json = serde_json::to_string(&feedback).expect("feedback");
        for private in [
            "ACGTACGT",
            "/private/",
            "PASSWORD",
            "previous_project",
            "expired_turn",
        ] {
            assert!(!json.contains(private));
        }
        assert!(
            feedback
                .rows
                .iter()
                .all(|r| r.applicability == AgentFeedbackApplicability::RecheckRequired)
        );
    }

    #[test]
    fn agent_feedback_omits_unattributable_receipts_without_erasing_local_history() {
        let session_id = new_agent_context_id();
        let turn_id = new_agent_context_id();
        let mut bound = AgentExecutionReceipt::new(
            "state-summary",
            AgentExecutionStatus::Completed,
            None,
            None,
        );
        bound.session_id = session_id.clone();
        bound.turn_id = Some(turn_id.clone());
        bound.suggestion_index = Some(1);
        let visible_turns = BTreeSet::from([turn_id.clone()]);
        let mut receipts = vec![bound.clone()];
        for (turn, index) in [
            (None, None),
            (None, Some(1)),
            (Some(turn_id.clone()), None),
            (Some(turn_id), Some(0)),
            (Some(new_agent_context_id()), Some(1)),
        ] {
            let mut receipt = bound.clone();
            receipt.turn_id = turn;
            receipt.suggestion_index = index;
            // Later unbound observations must not exhaust the projection's row budget.
            receipts.extend(std::iter::repeat_n(receipt, AGENT_EXECUTION_RECEIPT_LIMIT));
        }
        let before = serde_json::to_value(&receipts).expect("local log");
        let feedback =
            AgentExecutionFeedback::project(&session_id, None, &receipts, &visible_turns);
        feedback.validate().expect("valid projection");
        assert_eq!(feedback.rows.len(), 1);
        assert_eq!(feedback.rows[0].receipt.receipt_id, bound.receipt_id);
        assert_eq!(feedback.omitted_receipt_count, receipts.len() - 1);
        assert_eq!(
            before,
            serde_json::to_value(&receipts).expect("retained local log")
        );
        let no_conversation =
            AgentExecutionFeedback::project(&session_id, None, &receipts, &BTreeSet::new());
        assert!(no_conversation.rows.is_empty());
        assert_eq!(no_conversation.omitted_receipt_count, receipts.len());
    }

    #[test]
    fn agent_feedback_keeps_historical_completion_separate_from_applicability() {
        let mut engine = GentleEngine::from_state(crate::engine::ProjectState::default());
        let revision = AgentExecutionRevision::capture(&engine);
        let mut receipt = AgentExecutionReceipt::new(
            "state-summary",
            AgentExecutionStatus::Completed,
            Some(revision),
            Some(revision),
        );
        receipt.session_id = "session".into();
        let turn_id = new_agent_context_id();
        receipt.turn_id = Some(turn_id.clone());
        receipt.suggestion_index = Some(1);
        let visible_turns = BTreeSet::from([turn_id]);
        engine.display_state_mut().show_features = false;
        let feedback = AgentExecutionFeedback::project(
            "session",
            Some(AgentExecutionRevision::capture(&engine)),
            &[receipt.clone()],
            &visible_turns,
        );
        assert_eq!(
            feedback.rows[0].applicability,
            AgentFeedbackApplicability::SameStructuralRevision
        );
        engine.state_mut();
        let feedback = AgentExecutionFeedback::project(
            "session",
            Some(AgentExecutionRevision::capture(&engine)),
            &[receipt],
            &visible_turns,
        );
        assert_eq!(
            feedback.rows[0].applicability,
            AgentFeedbackApplicability::RecheckRequired
        );
        assert_eq!(
            feedback.rows[0].receipt.status,
            AgentExecutionStatus::Completed
        );
    }

    #[test]
    fn agent_feedback_rejects_unscoped_or_raw_wire_values() {
        let session_id = new_agent_context_id();
        let mut receipt = AgentExecutionReceipt::new(
            "state-summary",
            AgentExecutionStatus::Completed,
            None,
            None,
        );
        receipt.session_id = session_id.clone();
        let turn_id = new_agent_context_id();
        receipt.turn_id = Some(turn_id.clone());
        receipt.suggestion_index = Some(1);
        let feedback = AgentExecutionFeedback::project(
            &session_id,
            None,
            &[receipt],
            &BTreeSet::from([turn_id]),
        );
        feedback.validate().expect("valid receipt");
        let mut invalid = feedback.clone();
        invalid.rows[0].receipt.error_sha256 = Some("/private/result.fa".into());
        assert!(invalid.validate().is_err());
        let mut invalid = feedback.clone();
        invalid.rows[0].receipt.session_id = new_agent_context_id();
        assert!(invalid.validate().is_err());
        let mut invalid = feedback;
        invalid.rows[0].applicability = AgentFeedbackApplicability::SameStructuralRevision;
        assert!(invalid.validate().is_err());
    }
}
