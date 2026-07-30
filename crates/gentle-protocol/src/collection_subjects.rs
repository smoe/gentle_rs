//! Portable collection-subject and operation-lifting contracts.
//!
//! Collections retain their biological semantics: a logical set is not a
//! physical pool, and an arrangement has meaningful member order. Capability
//! policies declare those distinctions before an adapter exposes a collection
//! action.

use crate::{
    ArrangementId, CapabilitySource, ContainerId, EngineError, GeneSetProvenanceRow, SeqId,
};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::collections::BTreeSet;

/// Schema for one operation applied to a collection subject.
pub const COLLECTION_OPERATION_REPORT_SCHEMA: &str = "gentle.collection_operation.v1";
/// Schema for the curated capability-to-collection policy registry.
pub const COLLECTION_LIFT_POLICY_REGISTRY_SCHEMA: &str =
    "gentle.collection_lift_policy_registry.v1";
/// Canonical membership fingerprint algorithm carried in collection reports.
///
/// The SHA-256 input is the UTF-8 encoding returned by
/// [`canonical_collection_membership_json`]. Set-like subjects sort and
/// de-duplicate stable member ids; arrangements retain numeric position order.
pub const COLLECTION_MEMBERSHIP_FINGERPRINT_ALGORITHM: &str =
    "sha256_canonical_collection_members_v1";

/// Collection classes whose semantics affect lifting and membership identity.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
pub enum CollectionSubjectKind {
    #[default]
    ProjectSequences,
    Container,
    Arrangement,
    GeneSetResolution,
}

/// Portable reference to an engine-owned collection.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "source_kind", rename_all = "snake_case")]
pub enum CollectionSubjectRef {
    ProjectSequences {
        #[serde(default)]
        seq_ids: Vec<SeqId>,
    },
    Container {
        container_id: ContainerId,
    },
    Arrangement {
        arrangement_id: ArrangementId,
    },
    GeneSetResolution {
        report_id: String,
    },
}

impl Default for CollectionSubjectRef {
    fn default() -> Self {
        Self::ProjectSequences { seq_ids: vec![] }
    }
}

impl CollectionSubjectRef {
    pub fn kind(&self) -> CollectionSubjectKind {
        match self {
            Self::ProjectSequences { .. } => CollectionSubjectKind::ProjectSequences,
            Self::Container { .. } => CollectionSubjectKind::Container,
            Self::Arrangement { .. } => CollectionSubjectKind::Arrangement,
            Self::GeneSetResolution { .. } => CollectionSubjectKind::GeneSetResolution,
        }
    }
}

/// One source or derived member participating in a collection operation.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct CollectionMemberRef {
    /// Stable identity within the collection. Gene-set members use `dedup_key`.
    pub stable_member_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seq_id: Option<SeqId>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub container_id: Option<ContainerId>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_symbol: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    /// Required for semantically ordered subjects such as arrangements.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ordering_index: Option<usize>,
    /// Derived members point to their source member through this field.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parent_member_id: Option<String>,
    #[serde(default)]
    pub source_provenance: Vec<GeneSetProvenanceRow>,
}

/// How an operation consumes or transforms a collection.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CollectionLiftingMode {
    #[default]
    Map,
    Combine,
    Compare,
    Arrange,
    Derive,
}

/// Stable outcome vocabulary for one source or derived member.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CollectionMemberOutcome {
    Succeeded,
    Failed,
    Skipped,
    #[default]
    NotAssessed,
}

/// Typed reason why a capability cannot consume a collection subject.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum CollectionLiftRejectionReason {
    SingleSequenceOnly,
    IncompatibleMemberTypes,
    RequiresPhysicalPool,
    RequiresMaterializedMembers,
    UnsupportedSubjectKind,
}

/// Supported lifting behavior or an explicit typed rejection.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum CollectionLiftSupport {
    Supported {
        mode: CollectionLiftingMode,
        result_payload_kind: String,
    },
    Rejected {
        reason: CollectionLiftRejectionReason,
        detail: String,
    },
}

impl Default for CollectionLiftSupport {
    fn default() -> Self {
        Self::Rejected {
            reason: CollectionLiftRejectionReason::UnsupportedSubjectKind,
            detail: String::new(),
        }
    }
}

impl CollectionLiftSupport {
    pub fn mode(&self) -> Option<CollectionLiftingMode> {
        match self {
            Self::Supported { mode, .. } => Some(*mode),
            Self::Rejected { .. } => None,
        }
    }
}

/// Policy for one subject kind on one capability.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct CollectionSubjectLiftPolicy {
    pub subject_kind: CollectionSubjectKind,
    pub support: CollectionLiftSupport,
}

/// Curated collection policy attached to a named capability.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct CollectionCapabilityLiftPolicy {
    pub source: CapabilitySource,
    pub name: String,
    #[serde(default)]
    pub subjects: Vec<CollectionSubjectLiftPolicy>,
}

/// Root payload loaded from `docs/collection_lift_policies.json`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct CollectionLiftPolicyRegistry {
    pub schema: String,
    pub policies: Vec<CollectionCapabilityLiftPolicy>,
}

/// Status and produced-report links for one collection member.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct CollectionMemberStatusRow {
    pub member: CollectionMemberRef,
    pub outcome: CollectionMemberOutcome,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error: Option<EngineError>,
    #[serde(default)]
    pub produced_report_ids: Vec<String>,
}

/// Portable summary of one lifted collection operation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct CollectionOperationReport {
    pub schema: String,
    pub report_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub generated_at_unix_ms: u128,
    pub capability_source: CapabilitySource,
    pub capability_name: String,
    pub collection_subject: CollectionSubjectRef,
    pub lifting_mode: CollectionLiftingMode,
    pub lift_policy: CollectionSubjectLiftPolicy,
    pub fingerprint_algorithm: String,
    pub collection_membership_fingerprint_sha256: String,
    pub dry_run: bool,
    pub applied: bool,
    #[serde(default)]
    pub per_member_status: Vec<CollectionMemberStatusRow>,
    #[serde(default)]
    pub aggregate_warnings: Vec<String>,
    #[serde(default)]
    pub provenance: Vec<GeneSetProvenanceRow>,
}

impl Default for CollectionOperationReport {
    fn default() -> Self {
        Self {
            schema: String::new(),
            report_id: String::new(),
            op_id: None,
            run_id: None,
            generated_at_unix_ms: 0,
            capability_source: CapabilitySource::EngineOperation,
            capability_name: String::new(),
            collection_subject: CollectionSubjectRef::default(),
            lifting_mode: CollectionLiftingMode::default(),
            lift_policy: CollectionSubjectLiftPolicy::default(),
            fingerprint_algorithm: String::new(),
            collection_membership_fingerprint_sha256: String::new(),
            dry_run: false,
            applied: false,
            per_member_status: vec![],
            aggregate_warnings: vec![],
            provenance: vec![],
        }
    }
}

/// Produce the canonical JSON bytes that are hashed for a membership lock.
///
/// Logical sets, project selections, and containers are set-like for this
/// identity and therefore sort and de-duplicate stable ids. Arrangements keep
/// duplicate members and sort by numeric `ordering_index`, avoiding the
/// `_10`-before-`_2` ambiguity of lexicographic position strings.
pub fn canonical_collection_membership_json(
    subject_kind: CollectionSubjectKind,
    members: &[CollectionMemberRef],
) -> String {
    if subject_kind == CollectionSubjectKind::Arrangement {
        let mut ordered = members
            .iter()
            .map(|member| (member.ordering_index, member.stable_member_id.clone()))
            .collect::<Vec<_>>();
        ordered.sort_by(|left, right| {
            left.0
                .is_none()
                .cmp(&right.0.is_none())
                .then(left.0.cmp(&right.0))
                .then(left.1.cmp(&right.1))
        });
        return json!({
            "subject_kind": subject_kind,
            "order_semantics": "ordered",
            "members": ordered
                .into_iter()
                .map(|(ordering_index, stable_member_id)| json!({
                    "ordering_index": ordering_index,
                    "stable_member_id": stable_member_id,
                }))
                .collect::<Vec<_>>(),
        })
        .to_string();
    }

    let member_ids = members
        .iter()
        .map(|member| member.stable_member_id.clone())
        .collect::<BTreeSet<_>>();
    json!({
        "subject_kind": subject_kind,
        "order_semantics": "set",
        "members": member_ids,
    })
    .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_membership_is_set_like_except_for_arrangements() {
        let members = vec![
            CollectionMemberRef {
                stable_member_id: "member_10".to_string(),
                ordering_index: Some(10),
                ..CollectionMemberRef::default()
            },
            CollectionMemberRef {
                stable_member_id: "member_2".to_string(),
                ordering_index: Some(2),
                ..CollectionMemberRef::default()
            },
            CollectionMemberRef {
                stable_member_id: "member_2".to_string(),
                ordering_index: Some(11),
                ..CollectionMemberRef::default()
            },
        ];

        let logical = canonical_collection_membership_json(
            CollectionSubjectKind::GeneSetResolution,
            &members,
        );
        assert_eq!(
            logical,
            r#"{"members":["member_10","member_2"],"order_semantics":"set","subject_kind":"gene_set_resolution"}"#
        );

        let arrangement =
            canonical_collection_membership_json(CollectionSubjectKind::Arrangement, &members);
        assert_eq!(
            arrangement,
            r#"{"members":[{"ordering_index":2,"stable_member_id":"member_2"},{"ordering_index":10,"stable_member_id":"member_10"},{"ordering_index":11,"stable_member_id":"member_2"}],"order_semantics":"ordered","subject_kind":"arrangement"}"#
        );
    }

    #[test]
    fn collection_operation_report_round_trips_with_derived_member() {
        let report = CollectionOperationReport {
            schema: COLLECTION_OPERATION_REPORT_SCHEMA.to_string(),
            report_id: "collection_operation:op_1".to_string(),
            capability_source: CapabilitySource::EngineOperation,
            capability_name: "BuildGeneSetPromoterCohort".to_string(),
            collection_subject: CollectionSubjectRef::GeneSetResolution {
                report_id: "resolution:op_0".to_string(),
            },
            lifting_mode: CollectionLiftingMode::Derive,
            lift_policy: CollectionSubjectLiftPolicy {
                subject_kind: CollectionSubjectKind::GeneSetResolution,
                support: CollectionLiftSupport::Supported {
                    mode: CollectionLiftingMode::Derive,
                    result_payload_kind: "gentle.gene_set_promoter_cohort.v1".to_string(),
                },
            },
            fingerprint_algorithm: COLLECTION_MEMBERSHIP_FINGERPRINT_ALGORITHM.to_string(),
            collection_membership_fingerprint_sha256: "sha256:test".to_string(),
            applied: true,
            per_member_status: vec![CollectionMemberStatusRow {
                member: CollectionMemberRef {
                    stable_member_id: "promoter:gene_id:ensg1".to_string(),
                    parent_member_id: Some("gene_id:ensg1".to_string()),
                    ..CollectionMemberRef::default()
                },
                outcome: CollectionMemberOutcome::Succeeded,
                produced_report_ids: vec!["promoter_cohort:op_1".to_string()],
                ..CollectionMemberStatusRow::default()
            }],
            ..CollectionOperationReport::default()
        };

        let encoded = serde_json::to_string(&report).expect("serialize collection report");
        let decoded: CollectionOperationReport =
            serde_json::from_str(&encoded).expect("deserialize collection report");
        assert_eq!(decoded, report);
    }
}
