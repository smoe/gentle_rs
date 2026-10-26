//! Request-aware conversation approval for the explicitly migrated sequence-design preview.
//!
//! This is only the approval facet of the proposed lifecycle policy. It does not
//! change construction, reconciliation, history, persistence or destination safety.

use super::Operation;

/// Operation-derived approval only; candidate and host/tool gates remain additive.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum OperationApprovalPolicy {
    /// This concrete request adds no confirmation, not permission to bypass other gates.
    NoAdditionalConfirmation,
    /// Require review before dispatch, even if search or the subsequent write fails.
    ConfirmArtifactWriteAttempt,
}

impl OperationApprovalPolicy {
    pub(crate) fn requires_confirmation(self) -> bool {
        matches!(self, Self::ConfirmArtifactWriteAttempt)
    }
}

/// Select only from the typed request, without I/O or registry/model inference.
/// `None` leaves an unmigrated operation's existing conservative rules intact.
pub(crate) fn operation_approval_policy(op: &Operation) -> Option<OperationApprovalPolicy> {
    match op {
        Operation::PlanDnaSequenceDesign { path: Some(_), .. } => {
            Some(OperationApprovalPolicy::ConfirmArtifactWriteAttempt)
        }
        Operation::PlanDnaSequenceDesign { path: None, .. } => {
            Some(OperationApprovalPolicy::NoAdditionalConfirmation)
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequence_design_approval_policy_is_request_aware_and_legacy_conservative() {
        // Existing explicitly synthetic tutorial input; never a natural assay target.
        let request = serde_json::from_str(include_str!(
            "../../docs/tutorial/inputs/synthetic_sequence_design.json"
        ))
        .unwrap();
        let mut plan = Operation::PlanDnaSequenceDesign {
            request: Box::new(request),
            path: None,
        };
        assert_eq!(
            operation_approval_policy(&plan),
            Some(OperationApprovalPolicy::NoAdditionalConfirmation)
        );
        for path in ["preview.json", "", "  ", "invalid\0path"] {
            if let Operation::PlanDnaSequenceDesign { path: output, .. } = &mut plan {
                *output = Some(path.into());
            }
            assert_eq!(
                operation_approval_policy(&plan),
                Some(OperationApprovalPolicy::ConfirmArtifactWriteAttempt)
            );
        }
        let save = Operation::SaveFile {
            seq_id: "synthetic".into(),
            path: "preview.json".into(),
            format: super::super::ExportFormat::GenBank,
        };
        assert_eq!(operation_approval_policy(&save), None);
        assert_eq!(
            operation_approval_policy(&Operation::SetParameter {
                name: "synthetic".into(),
                value: serde_json::json!(1),
            }),
            None
        );
    }
}
