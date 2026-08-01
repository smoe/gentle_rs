//! Registry-backed presentation model for collection operation launchers.
//!
//! The GUI projects canonical engine-operation policies and keeps execution
//! adapters explicit. Policy support and GUI readiness are deliberately
//! separate so a missing frontend adapter is never presented as a biological
//! rejection.

use gentle_protocol::{
    CapabilitySource, CollectionLiftRejectionReason, CollectionLiftSupport, CollectionLiftingMode,
    CollectionSubjectKind, CollectionSubjectLiftPolicy, capability_descriptor,
    collection_lift_policy_registry,
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) enum CollectionLauncherAdapter {
    #[default]
    PrimerSpecificity,
    RestrictionScan,
    TfbsScan,
    Digest,
    PromoterCohort,
}

impl CollectionLauncherAdapter {
    pub(super) fn capability_name(self) -> &'static str {
        match self {
            Self::PrimerSpecificity => "AssessPrimerPairSpecificity",
            Self::RestrictionScan => "FindRestrictionSites",
            Self::TfbsScan => "ScanTfbsHits",
            Self::Digest => "Digest",
            Self::PromoterCohort => "BuildGeneSetPromoterCohort",
        }
    }

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::PrimerSpecificity => "Primer specificity",
            Self::RestrictionScan => "Restriction scan",
            Self::TfbsScan => "TFBS hit scan",
            Self::Digest => "Restriction digest",
            Self::PromoterCohort => "Promoter cohort",
        }
    }

    fn from_capability_name(name: &str) -> Option<Self> {
        match name {
            "AssessPrimerPairSpecificity" => Some(Self::PrimerSpecificity),
            "FindRestrictionSites" => Some(Self::RestrictionScan),
            "ScanTfbsHits" => Some(Self::TfbsScan),
            "Digest" => Some(Self::Digest),
            "BuildGeneSetPromoterCohort" => Some(Self::PromoterCohort),
            _ => None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct CollectionLauncherRow {
    pub(super) capability_name: String,
    pub(super) title: String,
    pub(super) description: String,
    pub(super) policy: CollectionSubjectLiftPolicy,
    pub(super) adapter: Option<CollectionLauncherAdapter>,
}

impl CollectionLauncherRow {
    pub(super) fn lifting_mode_label(&self) -> &'static str {
        match self.policy.support.mode() {
            Some(CollectionLiftingMode::Map) => "Map",
            Some(CollectionLiftingMode::Combine) => "Combine",
            Some(CollectionLiftingMode::Compare) => "Compare",
            Some(CollectionLiftingMode::Arrange) => "Arrange",
            Some(CollectionLiftingMode::Derive) => "Derive",
            None => "-",
        }
    }

    pub(super) fn result_payload_kind(&self) -> &str {
        match &self.policy.support {
            CollectionLiftSupport::Supported {
                result_payload_kind,
                ..
            } => result_payload_kind,
            CollectionLiftSupport::Rejected { .. } => "-",
        }
    }

    pub(super) fn baseline_readiness(&self) -> CollectionLauncherReadiness {
        match &self.policy.support {
            CollectionLiftSupport::Supported { .. } => match self.adapter {
                Some(_) => CollectionLauncherReadiness::Ready,
                None => CollectionLauncherReadiness::AdapterUnavailable {
                    detail: format!(
                        "'{}' has a supported collection policy but no typed GUI adapter",
                        self.capability_name
                    ),
                },
            },
            CollectionLiftSupport::Rejected { reason, detail } => match reason {
                CollectionLiftRejectionReason::RequiresMaterializedMembers => {
                    CollectionLauncherReadiness::RequiresMaterialization {
                        reason: *reason,
                        detail: detail.clone(),
                    }
                }
                CollectionLiftRejectionReason::RequiresPhysicalPool => {
                    CollectionLauncherReadiness::RequiresPhysicalPool {
                        reason: *reason,
                        detail: detail.clone(),
                    }
                }
                _ => CollectionLauncherReadiness::PolicyRejected {
                    reason: *reason,
                    detail: detail.clone(),
                },
            },
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum CollectionLauncherReadiness {
    Ready,
    NeedsInput {
        detail: String,
    },
    NeedsBindings {
        detail: String,
    },
    RequiresMaterialization {
        reason: CollectionLiftRejectionReason,
        detail: String,
    },
    RequiresPhysicalPool {
        reason: CollectionLiftRejectionReason,
        detail: String,
    },
    PolicyRejected {
        reason: CollectionLiftRejectionReason,
        detail: String,
    },
    AdapterUnavailable {
        detail: String,
    },
}

impl CollectionLauncherReadiness {
    pub(super) fn is_ready(&self) -> bool {
        matches!(self, Self::Ready)
    }

    pub(super) fn label(&self) -> &'static str {
        match self {
            Self::Ready => "Ready",
            Self::NeedsInput { .. } => "Needs input",
            Self::NeedsBindings { .. } => "Needs bindings",
            Self::RequiresMaterialization { .. } => "Requires materialization",
            Self::RequiresPhysicalPool { .. } => "Requires physical pool",
            Self::PolicyRejected { .. } => "Unsupported",
            Self::AdapterUnavailable { .. } => "GUI adapter unavailable",
        }
    }

    pub(super) fn detail(&self) -> Option<&str> {
        match self {
            Self::Ready => None,
            Self::NeedsInput { detail }
            | Self::NeedsBindings { detail }
            | Self::RequiresMaterialization { detail, .. }
            | Self::RequiresPhysicalPool { detail, .. }
            | Self::PolicyRejected { detail, .. }
            | Self::AdapterUnavailable { detail } => Some(detail),
        }
    }

    pub(super) fn rejection_reason(&self) -> Option<CollectionLiftRejectionReason> {
        match self {
            Self::RequiresMaterialization { reason, .. }
            | Self::RequiresPhysicalPool { reason, .. }
            | Self::PolicyRejected { reason, .. } => Some(*reason),
            _ => None,
        }
    }
}

pub(super) fn collection_launcher_rows(
    subject_kind: CollectionSubjectKind,
) -> Vec<CollectionLauncherRow> {
    let mut rows = collection_lift_policy_registry()
        .iter()
        .filter(|entry| entry.source == CapabilitySource::EngineOperation)
        .filter_map(|entry| {
            let policy = entry
                .subjects
                .iter()
                .find(|policy| policy.subject_kind == subject_kind)?;
            let descriptor = capability_descriptor(entry.source, &entry.name)?;
            Some(CollectionLauncherRow {
                capability_name: entry.name.clone(),
                title: descriptor
                    .title
                    .clone()
                    .unwrap_or_else(|| entry.name.clone()),
                description: descriptor.description.clone(),
                policy: policy.clone(),
                adapter: CollectionLauncherAdapter::from_capability_name(&entry.name),
            })
        })
        .collect::<Vec<_>>();
    rows.sort_by(|left, right| {
        left.title
            .to_ascii_lowercase()
            .cmp(&right.title.to_ascii_lowercase())
            .then_with(|| left.capability_name.cmp(&right.capability_name))
    });
    rows
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collection_operation_launcher_projects_every_canonical_gene_set_policy_once() {
        let rows = collection_launcher_rows(CollectionSubjectKind::GeneSetResolution);
        let projected = rows
            .iter()
            .map(|row| (row.capability_name.clone(), row.policy.clone()))
            .collect::<Vec<_>>();
        let mut expected = collection_lift_policy_registry()
            .iter()
            .filter(|entry| entry.source == CapabilitySource::EngineOperation)
            .filter_map(|entry| {
                entry
                    .subjects
                    .iter()
                    .find(|policy| policy.subject_kind == CollectionSubjectKind::GeneSetResolution)
                    .map(|policy| (entry.name.clone(), policy.clone()))
            })
            .collect::<Vec<_>>();
        expected.sort_by(|left, right| {
            let left_title = capability_descriptor(CapabilitySource::EngineOperation, &left.0)
                .and_then(|descriptor| descriptor.title.as_deref())
                .unwrap_or(&left.0)
                .to_ascii_lowercase();
            let right_title = capability_descriptor(CapabilitySource::EngineOperation, &right.0)
                .and_then(|descriptor| descriptor.title.as_deref())
                .unwrap_or(&right.0)
                .to_ascii_lowercase();
            left_title
                .cmp(&right_title)
                .then_with(|| left.0.cmp(&right.0))
        });
        assert_eq!(projected, expected);
        assert_eq!(
            rows.iter()
                .map(|row| row.capability_name.as_str())
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            rows.len(),
            "glossary aliases must not create duplicate launcher rows"
        );
    }

    #[test]
    fn collection_operation_launcher_engine_projection_covers_all_gene_set_policy_semantics() {
        let engine_policies = collection_lift_policy_registry()
            .iter()
            .filter(|entry| entry.source == CapabilitySource::EngineOperation)
            .filter_map(|entry| {
                entry
                    .subjects
                    .iter()
                    .find(|policy| policy.subject_kind == CollectionSubjectKind::GeneSetResolution)
            })
            .collect::<Vec<_>>();
        let uncovered = collection_lift_policy_registry()
            .iter()
            .filter(|entry| entry.source != CapabilitySource::EngineOperation)
            .filter_map(|entry| {
                let policy = entry.subjects.iter().find(|policy| {
                    policy.subject_kind == CollectionSubjectKind::GeneSetResolution
                })?;
                (!engine_policies.contains(&policy))
                    .then(|| format!("{}:{}", entry.source.as_str(), entry.name))
            })
            .collect::<Vec<_>>();
        assert!(
            uncovered.is_empty(),
            "these gene-set policies have no semantically identical canonical engine-operation row: {uncovered:?}"
        );
    }

    #[test]
    fn collection_operation_launcher_preserves_typed_policy_rejections() {
        let rows = collection_launcher_rows(CollectionSubjectKind::GeneSetResolution);
        let export_pool = rows
            .iter()
            .find(|row| row.capability_name == "ExportPool")
            .expect("ExportPool policy row");
        let readiness = export_pool.baseline_readiness();
        assert_eq!(
            readiness.rejection_reason(),
            Some(CollectionLiftRejectionReason::RequiresPhysicalPool)
        );
        assert!(
            readiness
                .detail()
                .is_some_and(|detail| detail.contains("logical gene set"))
        );
    }

    #[test]
    fn collection_operation_launcher_distinguishes_missing_adapter_from_policy_rejection() {
        let row = CollectionLauncherRow {
            capability_name: "FutureCollectionOperation".to_string(),
            title: "Future collection operation".to_string(),
            description: String::new(),
            policy: CollectionSubjectLiftPolicy {
                subject_kind: CollectionSubjectKind::GeneSetResolution,
                support: CollectionLiftSupport::Supported {
                    mode: CollectionLiftingMode::Map,
                    result_payload_kind: "gentle.future.v1".to_string(),
                },
                ..CollectionSubjectLiftPolicy::default()
            },
            adapter: None,
        };
        assert!(matches!(
            row.baseline_readiness(),
            CollectionLauncherReadiness::AdapterUnavailable { .. }
        ));
    }

    #[test]
    fn collection_digest_adapter_projects_supported_logical_and_rejected_physical_subjects() {
        for subject_kind in [
            CollectionSubjectKind::ProjectSequences,
            CollectionSubjectKind::GeneSetResolution,
        ] {
            let row = collection_launcher_rows(subject_kind)
                .into_iter()
                .find(|row| row.adapter == Some(CollectionLauncherAdapter::Digest))
                .expect("digest launcher row");
            assert!(row.baseline_readiness().is_ready());
            assert_eq!(
                row.policy.context_requirement,
                gentle_protocol::CollectionContextRequirement::ContextAgnostic
            );
            assert_eq!(row.result_payload_kind(), "gentle.collection_digest.v1");
        }

        let container = collection_launcher_rows(CollectionSubjectKind::Container)
            .into_iter()
            .find(|row| row.capability_name == "Digest")
            .expect("physical-container digest policy row");
        assert_eq!(
            container.baseline_readiness().rejection_reason(),
            Some(CollectionLiftRejectionReason::UnsupportedSubjectKind)
        );
        assert!(
            container
                .baseline_readiness()
                .detail()
                .is_some_and(|detail| detail.contains("DigestContainer"))
        );
    }

    #[test]
    fn container_collection_projection_exposes_pool_export_and_existing_gel_combine_policies() {
        let rows = collection_launcher_rows(CollectionSubjectKind::Container);
        for (capability_name, payload_kind) in [
            ("ExportPoolCollection", "gentle.collection_pool_export.v1"),
            ("RenderPoolGelSvg", "gentle.pool_gel_svg.v1"),
        ] {
            let row = rows
                .iter()
                .find(|row| row.capability_name == capability_name)
                .unwrap_or_else(|| panic!("missing container policy row for {capability_name}"));
            assert!(matches!(
                row.policy.support,
                CollectionLiftSupport::Supported {
                    mode: CollectionLiftingMode::Combine,
                    ..
                }
            ));
            assert_eq!(row.result_payload_kind(), payload_kind);
            assert_eq!(
                row.policy.context_requirement,
                gentle_protocol::CollectionContextRequirement::ContextAgnostic
            );
        }
    }
}
