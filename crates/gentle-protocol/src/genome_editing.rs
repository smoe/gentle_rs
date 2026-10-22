//! Read-only, reference-bound SNV intent contracts, not editor or efficacy models.
//!
//! Requests compose existing reference/context records but never accept a
//! caller-certified projection. Engine validation supplies the checked mapping.

use crate::{
    BiologicalContext, GenomicRegionLocalProjection, GenomicRegionReference, GenomicRegionStrand,
};
use serde::{Deserialize, Serialize};

pub const GENOME_EDIT_INTENT_REQUEST_SCHEMA: &str = "gentle.genome_edit_intent_request.v1";
pub const GENOME_EDIT_INTENT_REPORT_SCHEMA: &str = "gentle.genome_edit_intent_report.v1";
/// Bound both returned sequence copies and validation work, not a performance promise.
pub const MAX_GENOME_EDIT_WINDOW_BP: usize = 100_000;
pub const MAX_GENOME_EDIT_ID_BYTES: usize = 256;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum GenomeEditTarget {
    LoadedSequence { seq_id: String },
    InlineSequence { sequence: String },
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum GenomeEditReference {
    Genomic {
        reference: GenomicRegionReference,
    },
    /// A caller-declared design identity, never a fabricated biological accession.
    Synthetic {
        namespace: String,
        reference_id: String,
    },
}

/// Source-window geometry is 1-based inclusive; intent coordinates are not.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct GenomeEditAnchor {
    pub start_1based: u64,
    pub end_1based: u64,
    pub strand: GenomicRegionStrand,
}

/// One genomic-forward replacement. The first engine slice admits one SNV only.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct GenomeEditIntent {
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
    pub deleted: String,
    pub inserted: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct GenomeEditIntentRequest {
    pub schema: String,
    pub target: GenomeEditTarget,
    pub reference: GenomeEditReference,
    pub anchor: GenomeEditAnchor,
    /// SHA-256 of uppercase source DNA in its loaded/supplied orientation.
    pub expected_sequence_sha256: String,
    /// Explicit design-basis label (for example `reference`), not a measured genotype.
    pub design_haplotype: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub context: Option<BiologicalContext>,
    pub intent: GenomeEditIntent,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomeEditReferenceVerification {
    SyntheticDeclared,
    InlineCallerDeclared,
    /// Matched retained project provenance; no fresh reference-provider query.
    RetainedAnchorMatched,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomeEditAnnotationStatus {
    AnnotationNotSupplied,
    /// Source/release labels are bound, but no annotations were fetched or evaluated.
    SuppliedNotEvaluated,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct GenomeEditIntentReport {
    pub schema: String,
    pub request: GenomeEditIntentRequest,
    pub source_sequence: String,
    pub source_sequence_sha256: String,
    pub genomic_forward_sequence: String,
    pub genomic_forward_sha256: String,
    pub desired_genomic_forward_sequence: String,
    pub desired_genomic_forward_sha256: String,
    pub desired_local_sequence: String,
    pub desired_local_sha256: String,
    pub local_start_0based: u64,
    pub local_end_0based_exclusive: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub local_projection: Option<GenomicRegionLocalProjection>,
    pub reference_verification: GenomeEditReferenceVerification,
    pub annotation_status: GenomeEditAnnotationStatus,
    pub nonclaims: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    // Hand-written parser fixture; no real locus, guide or biological profile.
    fn request() -> serde_json::Value {
        json!({
            "schema": GENOME_EDIT_INTENT_REQUEST_SCHEMA,
            "target": {"kind":"inline_sequence", "sequence":"ACGTAC"},
            "reference": {"kind":"synthetic", "namespace":"synthetic:snv-test", "reference_id":"toy"},
            "anchor": {"start_1based":1001, "end_1based":1006, "strand":"plus"},
            "expected_sequence_sha256": format!("sha256:{}", "0".repeat(64)),
            "design_haplotype":"reference",
            "intent":{"start_0based":1001,"end_0based_exclusive":1002,"deleted":"C","inserted":"T"}
        })
    }

    #[test]
    fn genome_edit_intent_requires_fields_and_refuses_caller_projection_or_status() {
        assert!(serde_json::from_str::<GenomeEditIntentRequest>("{}").is_err());
        let base = request();
        let parsed: GenomeEditIntentRequest = serde_json::from_value(base.clone()).unwrap();
        assert_eq!(serde_json::to_value(parsed).unwrap(), base);
        for field in [
            "schema",
            "target",
            "reference",
            "anchor",
            "expected_sequence_sha256",
            "design_haplotype",
            "intent",
        ] {
            let mut missing = base.clone();
            missing.as_object_mut().unwrap().remove(field);
            assert!(
                serde_json::from_value::<GenomeEditIntentRequest>(missing).is_err(),
                "{field}"
            );
        }
        for field in ["local_projection", "status", "path", "editor"] {
            let mut extra = base.clone();
            extra[field] = json!({});
            assert!(
                serde_json::from_value::<GenomeEditIntentRequest>(extra).is_err(),
                "{field}"
            );
        }
        for field in ["target", "reference", "anchor", "intent"] {
            let mut extra = base.clone();
            extra[field]["status"] = json!("current");
            assert!(
                serde_json::from_value::<GenomeEditIntentRequest>(extra).is_err(),
                "{field}"
            );
        }
    }

    #[test]
    fn genome_edit_intent_registry_declares_read_only_inline_shell_reachability() {
        let descriptor = crate::capability_descriptor(
            crate::CapabilitySource::EngineOperation,
            "EvaluateGenomeEditIntent",
        )
        .unwrap();
        assert_eq!(descriptor.mutating, crate::CapabilityMutation::ReadOnly);
        assert_eq!(descriptor.inline_operand_ok, Some(true));
        for adapter in [
            crate::CapabilityAdapter::Gui,
            crate::CapabilityAdapter::Cli,
            crate::CapabilityAdapter::Mcp,
            crate::CapabilityAdapter::Js,
            crate::CapabilityAdapter::Lua,
        ] {
            assert_eq!(
                descriptor.surfacing_for_adapter(adapter),
                crate::AdapterSurfacing::ShellPassthrough
            );
        }
    }
}
