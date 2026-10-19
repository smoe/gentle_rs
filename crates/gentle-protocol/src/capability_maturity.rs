//! Portable capability evidence, independent of availability and permission.
//!
//! Discovery is always unassessed. A scoped assessment compares exact runtime,
//! input/resource and capability-definition bindings for one evidence kind.
//! Neither registry reachability nor a content hash authenticates an audit.

use crate::{CapabilityDescriptor, CapabilitySource, capability_registry};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::LazyLock;

pub const CAPABILITY_MATURITY_SCHEMA: &str = "gentle.clawbio_capability_maturity.v1";

/// Claimed maturity within one attestation, not a global software badge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MaturityLevel {
    Unassessed,
    Experimental,
    ReviewRequired,
    Verified,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MaturityEvidenceKind {
    Execution,
    NativeGui,
    BiologicalReview,
    Performance,
}

/// A retained audit outcome, separate from an operation's scientific result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MaturityVerdict {
    Pass,
    Fail,
    Incomplete,
    Withdrawn,
}

/// Whether exact registry references exist for a route's individual plan step.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MaturityStepScope {
    WrapperMetadata,
    RuntimeCommand,
    Workflow,
    RuntimeUnregistered,
}

/// A reference to an existing registry identity, not another operation catalog.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaturityRegistryRef {
    pub source: CapabilitySource,
    pub name: String,
}

/// Ordered step metadata; approval is still owned by the original plan/executor.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaturityStepBinding {
    pub step_index: usize,
    pub request_mode: String,
    pub scope: MaturityStepScope,
    pub registry_refs: Vec<MaturityRegistryRef>,
}

/// All fields are exact; no wildcard platform, inferred resource or feature set.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaturityEvidenceScope {
    pub source_sha: String,
    pub tool_version: String,
    pub binary_sha256: String,
    pub platform: String,
    pub architecture: String,
    pub build_profile: String,
    pub features: Vec<String>,
    pub input_class: String,
    pub resources: BTreeMap<String, String>,
}

/// Public repository-relative bytes; a digest is a binding, not authentication.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaturitySourceFile {
    pub path: String,
    pub sha256: String,
}

/// A retained public attestation; its assertions and non-claims stay inspectable.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CapabilityMaturityEvidence {
    pub evidence_id: String,
    pub capability_id: String,
    pub definition_sha256: String,
    pub scope: MaturityEvidenceScope,
    pub kind: MaturityEvidenceKind,
    pub verdict: MaturityVerdict,
    pub maturity: MaturityLevel,
    pub artifact: MaturitySourceFile,
    pub assertions: Vec<String>,
    pub limitations: Vec<String>,
}

/// One stable skill/intent identity with its current definition and historical evidence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CapabilityMaturityRecord {
    pub capability_id: String,
    pub skill: String,
    pub skill_version: String,
    pub intent_id: String,
    pub route_sha256: String,
    pub definition_sha256: String,
    pub steps: Vec<MaturityStepBinding>,
    pub maturity: MaturityLevel,
    pub evidence: Vec<CapabilityMaturityEvidence>,
}

/// Generated projection. No field authorizes execution or reports local readiness.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CapabilityMaturityCatalog {
    pub schema: String,
    pub runtime_binding: String,
    pub availability: String,
    pub policy: String,
    pub source_files: Vec<MaturitySourceFile>,
    pub capabilities: Vec<CapabilityMaturityRecord>,
}

/// Conservative result for the requested evidence kind; never execution permission.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MaturityAssessmentStatus {
    Unassessed,
    Experimental,
    ReviewRequired,
    Verified,
    Fail,
    Incomplete,
    Withdrawn,
}

/// Applies only to `kind` within the exact supplied scope, never all GENtle use.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CapabilityMaturityAssessment {
    pub capability_id: String,
    pub kind: MaturityEvidenceKind,
    pub scope: MaturityEvidenceScope,
    pub status: MaturityAssessmentStatus,
    pub evidence_ids: Vec<String>,
    pub reason: String,
}

fn hex(raw: &str, length: usize) -> bool {
    raw.len() == length
        && raw
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

fn hash(raw: &str) -> bool {
    raw.strip_prefix("sha256:").is_some_and(|raw| hex(raw, 64))
}

fn text(raw: &str) -> bool {
    !raw.trim().is_empty()
}

fn known_label(raw: &str) -> bool {
    text(raw)
        && !raw.contains('*')
        && !matches!(
            raw.trim().to_ascii_lowercase().as_str(),
            "unknown"
                | "unavailable"
                | "not_assessed"
                | "not_checked"
                | "unspecified"
                | "any"
                | "all"
        )
}

fn relative(raw: &str) -> bool {
    text(raw) && !raw.contains(['\\', ':']) && raw.split('/').all(|p| !matches!(p, "" | "." | ".."))
}

fn source_valid(source: &MaturitySourceFile) -> bool {
    relative(&source.path) && hash(&source.sha256)
}

impl MaturityEvidenceScope {
    pub fn validate(&self) -> Result<(), String> {
        let labels = [
            &self.tool_version,
            &self.platform,
            &self.architecture,
            &self.build_profile,
            &self.input_class,
        ];
        if !hex(&self.source_sha, 40)
            || !hash(&self.binary_sha256)
            || labels.iter().any(|s| !known_label(s))
            || self.features.iter().any(|s| !text(s) || s.contains('*'))
            || self.features.windows(2).any(|w| w[0] >= w[1])
            || self
                .resources
                .iter()
                .any(|(name, sha)| !text(name) || name.contains('*') || !hash(sha))
        {
            return Err("invalid_maturity_evidence_scope".into());
        }
        Ok(())
    }
}

impl CapabilityMaturityCatalog {
    pub fn validate(&self) -> Result<(), String> {
        if self.schema != CAPABILITY_MATURITY_SCHEMA
            || self.runtime_binding != "not_checked"
            || self.availability != "not_checked"
            || self.policy != "existing_execution_checks"
            || self.source_files.is_empty()
            || self.capabilities.is_empty()
            || self.source_files.iter().any(|s| !source_valid(s))
            || self.source_files.windows(2).any(|w| w[0].path >= w[1].path)
            || self
                .capabilities
                .windows(2)
                .any(|w| w[0].capability_id >= w[1].capability_id)
        {
            return Err("invalid_maturity_catalog".into());
        }
        let mut evidence_ids = BTreeSet::new();
        for row in &self.capabilities {
            if !text(&row.skill)
                || !text(&row.skill_version)
                || !text(&row.intent_id)
                || row.capability_id != format!("{}/{}", row.skill, row.intent_id)
                || !hash(&row.route_sha256)
                || !hash(&row.definition_sha256)
                || row.maturity != MaturityLevel::Unassessed
                || row.steps.is_empty()
            {
                return Err("invalid_maturity_capability".into());
            }
            for (index, step) in row.steps.iter().enumerate() {
                let refs = step
                    .registry_refs
                    .iter()
                    .map(|r| (r.source, &r.name))
                    .collect::<BTreeSet<_>>();
                if step.step_index != index
                    || !text(&step.request_mode)
                    || step.registry_refs.iter().any(|r| !text(&r.name))
                    || refs.len() != step.registry_refs.len()
                    || (refs.is_empty()
                        && !matches!(
                            step.scope,
                            MaturityStepScope::WrapperMetadata
                                | MaturityStepScope::RuntimeUnregistered
                        ))
                {
                    return Err("invalid_maturity_step_binding".into());
                }
            }
            for record in &row.evidence {
                record.scope.validate()?;
                if record.capability_id != row.capability_id
                    || !hash(&record.definition_sha256)
                    || !text(&record.evidence_id)
                    || !evidence_ids.insert(&record.evidence_id)
                    || !source_valid(&record.artifact)
                    || record.assertions.is_empty()
                    || record.limitations.is_empty()
                    || record
                        .assertions
                        .iter()
                        .chain(&record.limitations)
                        .any(|s| !text(s))
                    || (record.verdict != MaturityVerdict::Pass
                        && record.maturity != MaturityLevel::Unassessed)
                {
                    return Err("invalid_maturity_evidence".into());
                }
            }
        }
        Ok(())
    }

    /// Does not probe binaries, resources, projects or external services.
    pub fn assess(
        &self,
        capability_id: &str,
        scope: &MaturityEvidenceScope,
        kind: MaturityEvidenceKind,
        registry: &[CapabilityDescriptor],
    ) -> Result<CapabilityMaturityAssessment, String> {
        self.validate()?;
        scope.validate()?;
        let mut assessment = CapabilityMaturityAssessment {
            capability_id: capability_id.into(),
            kind,
            scope: scope.clone(),
            status: MaturityAssessmentStatus::Unassessed,
            evidence_ids: vec![],
            reason: "no_matching_scoped_evidence".into(),
        };
        let Some(row) = self
            .capabilities
            .iter()
            .find(|row| row.capability_id == capability_id)
        else {
            assessment.reason = "unknown_capability".into();
            return Ok(assessment);
        };
        if row.steps.iter().any(|step| {
            step.scope == MaturityStepScope::RuntimeUnregistered
                || step.registry_refs.iter().any(|reference| {
                    registry
                        .iter()
                        .filter(|item| {
                            item.source == reference.source && item.name == reference.name
                        })
                        .count()
                        != 1
                })
        }) {
            assessment.reason = "registry_binding_missing_or_ambiguous".into();
            return Ok(assessment);
        }
        let matches = row
            .evidence
            .iter()
            .filter(|e| {
                e.scope == *scope && e.kind == kind && e.definition_sha256 == row.definition_sha256
            })
            .collect::<Vec<_>>();
        assessment.evidence_ids = matches.iter().map(|e| e.evidence_id.clone()).collect();
        assessment.evidence_ids.sort();
        if matches.is_empty() {
            return Ok(assessment);
        }
        assessment.reason = "exact_scoped_attestations_only_not_permission_or_availability".into();
        assessment.status = if matches
            .iter()
            .any(|e| e.verdict == MaturityVerdict::Withdrawn)
        {
            MaturityAssessmentStatus::Withdrawn
        } else if matches.iter().any(|e| e.verdict == MaturityVerdict::Fail) {
            MaturityAssessmentStatus::Fail
        } else if matches
            .iter()
            .any(|e| e.verdict == MaturityVerdict::Incomplete)
        {
            MaturityAssessmentStatus::Incomplete
        } else if matches
            .iter()
            .any(|e| e.maturity == MaturityLevel::Unassessed)
        {
            MaturityAssessmentStatus::Unassessed
        } else if matches
            .iter()
            .any(|e| e.maturity == MaturityLevel::Experimental)
        {
            MaturityAssessmentStatus::Experimental
        } else if matches
            .iter()
            .any(|e| e.maturity == MaturityLevel::ReviewRequired)
        {
            MaturityAssessmentStatus::ReviewRequired
        } else {
            MaturityAssessmentStatus::Verified
        };
        Ok(assessment)
    }
}

static BUNDLED: LazyLock<Result<CapabilityMaturityCatalog, String>> = LazyLock::new(|| {
    let catalog: CapabilityMaturityCatalog =
        serde_json::from_str(include_str!("../../../assets/capability_maturity.json"))
            .map_err(|error| format!("invalid_bundled_maturity_catalog: {error}"))?;
    catalog.validate()?;
    Ok(catalog)
});

/// Shared read-only discovery, including diagnostics rather than inferred coverage.
pub fn bundled_capability_maturity() -> serde_json::Value {
    match BUNDLED.as_ref() {
        Err(error) => serde_json::json!({"schema": CAPABILITY_MATURITY_SCHEMA,
            "status": "invalid", "diagnostic": error, "maturity": "unassessed"}),
        Ok(catalog) => {
            let mut diagnostics = vec![];
            for row in &catalog.capabilities {
                for step in &row.steps {
                    if step.scope == MaturityStepScope::RuntimeUnregistered {
                        diagnostics.push(serde_json::json!({"capability_id": row.capability_id,
                            "step_index": step.step_index, "status": "unregistered"}));
                    }
                    for reference in &step.registry_refs {
                        let count = capability_registry()
                            .iter()
                            .filter(|item| {
                                item.source == reference.source && item.name == reference.name
                            })
                            .count();
                        if count != 1 {
                            diagnostics.push(
                                serde_json::json!({"capability_id": row.capability_id,
                                "step_index": step.step_index, "reference": reference,
                                "status": if count == 0 { "missing" } else { "ambiguous" }}),
                            );
                        }
                    }
                }
            }
            serde_json::json!({"status": "available", "catalog": catalog,
                "registry_diagnostics": diagnostics})
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // All added evidence below is hand-crafted synthetic metadata, never an audit.
    fn fixture() -> (CapabilityMaturityCatalog, MaturityEvidenceScope, String) {
        let mut catalog = BUNDLED.as_ref().expect("bundled projection").clone();
        let row = catalog
            .capabilities
            .iter_mut()
            .find(|r| r.intent_id == "primer_preflight")
            .unwrap();
        let id = row.capability_id.clone();
        let scope = MaturityEvidenceScope {
            source_sha: "a".repeat(40),
            tool_version: "synthetic".into(),
            binary_sha256: format!("sha256:{}", "b".repeat(64)),
            platform: "linux".into(),
            architecture: "x86_64".into(),
            build_profile: "synthetic".into(),
            features: vec![],
            input_class: "synthetic".into(),
            resources: BTreeMap::new(),
        };
        row.evidence.push(CapabilityMaturityEvidence {
            evidence_id: "synthetic-execution".into(),
            capability_id: id.clone(),
            definition_sha256: row.definition_sha256.clone(),
            scope: scope.clone(),
            kind: MaturityEvidenceKind::Execution,
            verdict: MaturityVerdict::Pass,
            maturity: MaturityLevel::Verified,
            artifact: MaturitySourceFile {
                path: "synthetic.json".into(),
                sha256: format!("sha256:{}", "c".repeat(64)),
            },
            assertions: vec!["synthetic execution only".into()],
            limitations: vec!["not biological or native acceptance".into()],
        });
        (catalog, scope, id)
    }

    #[test]
    fn capability_maturity_bundled_routes_are_unassessed() {
        let catalog = BUNDLED.as_ref().expect("bundled catalog");
        assert_eq!(catalog.capabilities.len(), 97);
        assert!(
            catalog
                .capabilities
                .iter()
                .all(|r| r.maturity == MaturityLevel::Unassessed && r.evidence.is_empty())
        );
    }

    #[test]
    fn capability_maturity_exact_scope_never_upgrades_other_kinds_or_platforms() {
        let (catalog, scope, id) = fixture();
        let assess = |scope: &MaturityEvidenceScope, kind| {
            catalog
                .assess(&id, scope, kind, capability_registry())
                .unwrap()
                .status
        };
        assert_eq!(
            assess(&scope, MaturityEvidenceKind::Execution),
            MaturityAssessmentStatus::Verified
        );
        assert_eq!(
            assess(&scope, MaturityEvidenceKind::BiologicalReview),
            MaturityAssessmentStatus::Unassessed
        );
        assert_eq!(
            assess(&scope, MaturityEvidenceKind::NativeGui),
            MaturityAssessmentStatus::Unassessed
        );
        let mut scopes = vec![];
        let mut changed = scope.clone();
        changed.platform = "windows".into();
        scopes.push(changed);
        let mut changed = scope.clone();
        changed.source_sha = "d".repeat(40);
        scopes.push(changed);
        let mut changed = scope.clone();
        changed.binary_sha256 = format!("sha256:{}", "e".repeat(64));
        scopes.push(changed);
        let mut changed = scope.clone();
        changed.input_class = "real_data".into();
        scopes.push(changed);
        let mut changed = scope.clone();
        changed
            .resources
            .insert("genome".into(), format!("sha256:{}", "f".repeat(64)));
        scopes.push(changed);
        for changed in scopes {
            assert_eq!(
                assess(&changed, MaturityEvidenceKind::Execution),
                MaturityAssessmentStatus::Unassessed
            );
        }
    }

    #[test]
    fn capability_maturity_conflicts_and_ambiguous_registry_fail_closed() {
        let (mut catalog, scope, id) = fixture();
        let row = catalog
            .capabilities
            .iter_mut()
            .find(|r| r.capability_id == id)
            .unwrap();
        let mut failure = row.evidence[0].clone();
        failure.evidence_id = "synthetic-failure".into();
        failure.verdict = MaturityVerdict::Fail;
        failure.maturity = MaturityLevel::Unassessed;
        row.evidence.push(failure);
        assert_eq!(
            catalog
                .assess(
                    &id,
                    &scope,
                    MaturityEvidenceKind::Execution,
                    capability_registry()
                )
                .unwrap()
                .status,
            MaturityAssessmentStatus::Fail
        );
        let mut registry = capability_registry().to_vec();
        let reference = &catalog
            .capabilities
            .iter()
            .find(|r| r.capability_id == id)
            .unwrap()
            .steps[0]
            .registry_refs[0];
        registry.push(
            registry
                .iter()
                .find(|r| r.source == reference.source && r.name == reference.name)
                .unwrap()
                .clone(),
        );
        assert_eq!(
            catalog
                .assess(&id, &scope, MaturityEvidenceKind::Execution, &registry)
                .unwrap()
                .status,
            MaturityAssessmentStatus::Unassessed
        );
    }

    #[test]
    fn capability_maturity_tampering_and_unknown_fields_are_rejected() {
        let (mut catalog, scope, id) = fixture();
        let row = catalog
            .capabilities
            .iter_mut()
            .find(|r| r.capability_id == id)
            .unwrap();
        row.definition_sha256 = format!("sha256:{}", "f".repeat(64));
        assert_eq!(
            catalog
                .assess(
                    &id,
                    &scope,
                    MaturityEvidenceKind::Execution,
                    capability_registry()
                )
                .unwrap()
                .status,
            MaturityAssessmentStatus::Unassessed
        );
        let mut value = serde_json::to_value(&catalog).unwrap();
        value["auto_execute"] = serde_json::json!(true);
        assert!(serde_json::from_value::<CapabilityMaturityCatalog>(value).is_err());
    }

    #[test]
    fn capability_maturity_negative_verdicts_and_lower_maturity_are_not_overridden() {
        for (verdict, maturity, expected) in [
            (
                MaturityVerdict::Incomplete,
                MaturityLevel::Unassessed,
                MaturityAssessmentStatus::Incomplete,
            ),
            (
                MaturityVerdict::Withdrawn,
                MaturityLevel::Unassessed,
                MaturityAssessmentStatus::Withdrawn,
            ),
            (
                MaturityVerdict::Pass,
                MaturityLevel::Experimental,
                MaturityAssessmentStatus::Experimental,
            ),
            (
                MaturityVerdict::Pass,
                MaturityLevel::ReviewRequired,
                MaturityAssessmentStatus::ReviewRequired,
            ),
        ] {
            let (mut catalog, scope, id) = fixture();
            let row = catalog
                .capabilities
                .iter_mut()
                .find(|r| r.capability_id == id)
                .unwrap();
            let mut other = row.evidence[0].clone();
            other.evidence_id = "synthetic-other".into();
            other.verdict = verdict;
            other.maturity = maturity;
            row.evidence.push(other);
            assert_eq!(
                catalog
                    .assess(
                        &id,
                        &scope,
                        MaturityEvidenceKind::Execution,
                        capability_registry()
                    )
                    .unwrap()
                    .status,
                expected
            );
        }
    }

    #[test]
    fn capability_maturity_scope_and_step_coverage_cannot_be_unknown_or_partial() {
        let (catalog, scope, id) = fixture();
        for field in ["tool_version", "architecture", "build_profile"] {
            let mut value = serde_json::to_value(&scope).unwrap();
            value[field] = serde_json::json!("not_checked");
            let invalid: MaturityEvidenceScope = serde_json::from_value(value).unwrap();
            assert!(invalid.validate().is_err());
        }
        for field in ["tool_version", "architecture", "build_profile", "features"] {
            let mut value = serde_json::to_value(&scope).unwrap();
            value[field] = if field == "features" {
                serde_json::json!(["desktop-gui"])
            } else {
                serde_json::json!("different")
            };
            let changed: MaturityEvidenceScope = serde_json::from_value(value).unwrap();
            assert_eq!(
                catalog
                    .assess(
                        &id,
                        &changed,
                        MaturityEvidenceKind::Execution,
                        capability_registry()
                    )
                    .unwrap()
                    .status,
                MaturityAssessmentStatus::Unassessed
            );
        }
        let mut value = serde_json::to_value(&scope).unwrap();
        value.as_object_mut().unwrap().remove("binary_sha256");
        assert!(serde_json::from_value::<MaturityEvidenceScope>(value).is_err());
        assert_eq!(
            catalog
                .assess(&id, &scope, MaturityEvidenceKind::Execution, &[])
                .unwrap()
                .status,
            MaturityAssessmentStatus::Unassessed
        );
        let mut catalog = catalog;
        let row = catalog
            .capabilities
            .iter_mut()
            .find(|r| r.capability_id == id)
            .unwrap();
        row.steps[0].scope = MaturityStepScope::RuntimeUnregistered;
        assert_eq!(
            catalog
                .assess(
                    &id,
                    &scope,
                    MaturityEvidenceKind::Execution,
                    capability_registry()
                )
                .unwrap()
                .status,
            MaturityAssessmentStatus::Unassessed
        );
    }
}
