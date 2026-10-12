//! Portable synthetic-CDS motif-removal requests, preview reports and exact approvals.

use serde::{Deserialize, Serialize};

pub const REQUEST_SCHEMA: &str = "gentle.dna_sequence_design_request.v1";
pub const REPORT_SCHEMA: &str = "gentle.dna_sequence_design_report.v1";
pub const RECEIPT_SCHEMA: &str = "gentle.dna_sequence_design_receipt.v1";

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignInterval {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DesignStrand {
    Forward,
    Reverse,
    Both,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AvoidDesignMotif {
    pub pattern: String,
    pub strand: DesignStrand,
}

/// Inclusive basis-point GC limits on the complete declared CDS, not flanks.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignGcBounds {
    pub min_basis_points: u16,
    pub max_basis_points: u16,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignGcMeasurement {
    pub gc_bases: usize,
    pub satisfies_bounds: bool,
}

/// Integer-only facts. Reachable extrema are a relaxation, not a feasibility proof.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignGcAssessment {
    pub denominator_bases: usize,
    pub minimum_gc_bases: usize,
    pub maximum_gc_bases: usize,
    pub reachable_minimum_gc_bases: usize,
    pub reachable_maximum_gc_bases: usize,
    pub input: DesignGcMeasurement,
    pub output: Option<DesignGcMeasurement>,
}

/// Inclusive GC bounds on all full windows inside the CDS, with a fixed 1-bp step.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignWindowGcBounds {
    pub window_bp: usize,
    pub min_basis_points: u16,
    pub max_basis_points: u16,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignWindowGcMeasurement {
    pub interval: DesignInterval,
    pub gc_bases: usize,
    pub satisfies_bounds: bool,
}

/// Every complete window in ascending local order. Missing output is not zero.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignWindowGcAssessment {
    pub window_bp: usize,
    pub window_count: usize,
    pub minimum_gc_bases: usize,
    pub maximum_gc_bases: usize,
    pub input: Vec<DesignWindowGcMeasurement>,
    pub output: Option<Vec<DesignWindowGcMeasurement>>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum DnaDesignTarget {
    InlineSequence { sequence: String },
    LoadedSequence { seq_id: String },
}

/// Explicit solver policy; omitted defaults retain legacy report/digest bytes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DnaDesignSearchStrategy {
    #[default]
    FullEnumeration,
    ConflictDirected,
}

impl DnaDesignSearchStrategy {
    pub fn is_full_enumeration(&self) -> bool {
        *self == Self::FullEnumeration
    }
}

/// Only explicitly declared synthetic coding inserts are admitted; annotation labels
/// never infer permission or CDS geometry. Protein excludes the terminal stop.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DnaSequenceDesignRequest {
    pub schema: String,
    pub target: DnaDesignTarget,
    pub purpose: String,
    pub cds: DesignInterval,
    pub protein_sequence: String,
    pub genetic_code: usize,
    #[serde(default)]
    pub protected_intervals: Vec<DesignInterval>,
    pub avoid_motifs: Vec<AvoidDesignMotif>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gc_content: Option<DesignGcBounds>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gc_window: Option<DesignWindowGcBounds>,
    pub max_evaluations: u64,
    #[serde(
        default,
        skip_serializing_if = "DnaDesignSearchStrategy::is_full_enumeration"
    )]
    pub search_strategy: DnaDesignSearchStrategy,
    pub output_seq_id: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignCodon {
    pub codon: String,
    pub residue: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DnaDesignStatus {
    Invalid,
    Unsupported,
    Feasible,
    SearchExhausted,
    ProvenInfeasible,
    Cancelled,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignMotifMatch {
    pub motif_index: usize,
    pub interval: DesignInterval,
    pub strand: DesignStrand,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DesignNucleotideEdit {
    pub position_0based: usize,
    pub before: String,
    pub after: String,
}

/// A feasible preview carries an exact approval digest. No other outcome does.
/// This report is not an order or functional/experimental approval.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DnaSequenceDesignReport {
    pub schema: String,
    pub request: DnaSequenceDesignRequest,
    pub source_sequence: String,
    pub source_sha256: String,
    pub source_features_sha256: Option<String>,
    pub omitted_source_feature_count: usize,
    pub genetic_code_mapping: Vec<DesignCodon>,
    pub genetic_code_mapping_sha256: String,
    pub algorithm: String,
    pub status: DnaDesignStatus,
    pub reason: String,
    pub output_sequence: Option<String>,
    pub output_sha256: Option<String>,
    pub edits: Vec<DesignNucleotideEdit>,
    pub initial_matches: Vec<DesignMotifMatch>,
    pub evaluated_candidates: u64,
    pub effective_evaluation_budget: u64,
    pub search_space: Option<u64>,
    /// Preview search claims; apply does not authenticate or re-prove them.
    pub optimization_complete: bool,
    /// Minimum-edit claim from the preview search, not apply-time verification.
    pub minimum_edits_proven: bool,
    pub approval_digest: Option<String>,
    pub nonclaims: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gc_content: Option<DesignGcAssessment>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gc_window: Option<DesignWindowGcAssessment>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DnaSequenceDesignReceipt {
    pub schema: String,
    pub approval_digest: String,
    pub created_seq_id: String,
    pub output_sha256: String,
    /// Apply rechecked the DNA constraints and exact nucleotide edits.
    /// Missing legacy fields do not establish verification.
    #[serde(default)]
    pub output_constraints_verified: bool,
    /// A content-bound preview is not authenticated search provenance.
    /// Currently always false: apply does not rerun either search strategy.
    #[serde(default)]
    pub search_claims_verified: bool,
    pub nonclaims: Vec<String>,
}

/// Submitted preview content, never authenticated or re-proven by apply.
/// The tag is mandatory and deliberately has no verified variant.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "verification", rename_all = "snake_case", deny_unknown_fields)]
pub enum DnaSequenceDesignSubmittedProposal {
    UnverifiedPortablePreview { report: DnaSequenceDesignReport },
}

/// Persisted `dna_sequence_design:OUTPUT_ID` metadata, separate from approval bytes.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DnaSequenceDesignMaterializationRecord {
    pub submitted_proposal: DnaSequenceDesignSubmittedProposal,
    pub receipt: DnaSequenceDesignReceipt,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequence_design_legacy_receipt_does_not_imply_verified_claims() {
        // Hand-crafted legacy receipt shape, not retained scientific evidence.
        let legacy = serde_json::json!({
            "schema": RECEIPT_SCHEMA,
            "approval_digest": "sha256:synthetic-example",
            "created_seq_id": "synthetic_without_ecori",
            "output_sha256": "sha256:synthetic-output",
            "nonclaims": []
        });
        let mut receipt: DnaSequenceDesignReceipt = serde_json::from_value(legacy).unwrap();
        assert!(!receipt.output_constraints_verified);
        assert!(!receipt.search_claims_verified);
        receipt.output_constraints_verified = true;
        let roundtrip: DnaSequenceDesignReceipt =
            serde_json::from_value(serde_json::to_value(&receipt).unwrap()).unwrap();
        assert_eq!(roundtrip, receipt);
        assert!(!roundtrip.search_claims_verified);
    }

    #[test]
    fn sequence_design_receipt_rejects_unknown_verification_fields() {
        // Synthetic receipt, recreated literally; not a scientific result.
        let receipt = serde_json::json!({
            "schema": RECEIPT_SCHEMA,
            "approval_digest": "sha256:synthetic-example",
            "created_seq_id": "synthetic_without_ecori",
            "output_sha256": "sha256:synthetic-output",
            "output_constraints_verified": true,
            "search_claims_verified": false,
            "nonclaims": []
        });
        let decoded: DnaSequenceDesignReceipt = serde_json::from_value(receipt.clone()).unwrap();
        assert_eq!(serde_json::to_value(decoded).unwrap(), receipt);
        let mut unknown = receipt;
        unknown["minimum_edits_verified"] = true.into();
        assert!(serde_json::from_value::<DnaSequenceDesignReceipt>(unknown).is_err());
    }

    #[test]
    fn sequence_design_strategy_preserves_legacy_request_bytes_and_rejects_unknown_modes() {
        // Literal synthetic MEF insert; recreates the documented request, not a gene.
        let legacy = serde_json::json!({
            "schema": REQUEST_SCHEMA,
            "target": {"kind": "inline_sequence", "sequence": "ATGGAATTCTAA"},
            "purpose": "synthetic_coding_insert",
            "cds": {"start_0based": 0, "end_0based_exclusive": 12},
            "protein_sequence": "MEF", "genetic_code": 1,
            "protected_intervals": [],
            "avoid_motifs": [{"pattern": "GAATTC", "strand": "both"}],
            "max_evaluations": 4096, "output_seq_id": "synthetic_without_ecori"
        });
        let request: DnaSequenceDesignRequest = serde_json::from_value(legacy.clone()).unwrap();
        assert_eq!(
            request.search_strategy,
            DnaDesignSearchStrategy::FullEnumeration
        );
        assert_eq!(serde_json::to_value(&request).unwrap(), legacy);
        assert_eq!(
            serde_json::to_string(&request).unwrap(),
            concat!(
                "{\"schema\":\"gentle.dna_sequence_design_request.v1\",",
                "\"target\":{\"kind\":\"inline_sequence\",\"sequence\":\"ATGGAATTCTAA\"},",
                "\"purpose\":\"synthetic_coding_insert\",",
                "\"cds\":{\"start_0based\":0,\"end_0based_exclusive\":12},",
                "\"protein_sequence\":\"MEF\",\"genetic_code\":1,",
                "\"protected_intervals\":[],",
                "\"avoid_motifs\":[{\"pattern\":\"GAATTC\",\"strand\":\"both\"}],",
                "\"max_evaluations\":4096,\"output_seq_id\":\"synthetic_without_ecori\"}"
            )
        );
        let mut explicit = legacy.clone();
        explicit["search_strategy"] = "full_enumeration".into();
        let request: DnaSequenceDesignRequest = serde_json::from_value(explicit.clone()).unwrap();
        assert_eq!(serde_json::to_value(request).unwrap(), legacy);
        explicit["search_strategy"] = "conflict_directed".into();
        let request: DnaSequenceDesignRequest = serde_json::from_value(explicit.clone()).unwrap();
        assert_eq!(
            request.search_strategy,
            DnaDesignSearchStrategy::ConflictDirected
        );
        assert_eq!(serde_json::to_value(request).unwrap(), explicit);
        explicit["search_strategy"] = "automatic".into();
        assert!(serde_json::from_value::<DnaSequenceDesignRequest>(explicit).is_err());
    }

    #[test]
    fn sequence_design_gc_roundtrip_is_additive_and_rejects_unknown_gc_fields() {
        // Synthetic MFK request, recreated literally; no natural gene is represented.
        let wire = serde_json::json!({
            "schema": REQUEST_SCHEMA,
            "target": {"kind": "inline_sequence", "sequence": "ATGTTTAAATAA"},
            "purpose": "synthetic_coding_insert",
            "cds": {"start_0based": 0, "end_0based_exclusive": 12},
            "protein_sequence": "MFK", "genetic_code": 1,
            "protected_intervals": [], "avoid_motifs": [],
            "gc_content": {"min_basis_points": 1666, "max_basis_points": 1667},
            "max_evaluations": 4096, "output_seq_id": "synthetic_gc"
        });
        let request: DnaSequenceDesignRequest = serde_json::from_value(wire.clone()).unwrap();
        assert_eq!(serde_json::to_value(request).unwrap(), wire);
        let mut unknown = wire.clone();
        unknown["gc_content"]["window_bp"] = 6.into();
        assert!(serde_json::from_value::<DnaSequenceDesignRequest>(unknown).is_err());
        let mut fractional = wire;
        fractional["gc_content"]["min_basis_points"] = serde_json::json!(16.66);
        assert!(serde_json::from_value::<DnaSequenceDesignRequest>(fractional).is_err());
    }

    #[test]
    fn sequence_design_window_gc_roundtrip_preserves_explicit_policy_and_rejects_stride() {
        // Synthetic MFK request with every 6-bp CDS window, recreated literally.
        let wire = serde_json::json!({
            "schema": REQUEST_SCHEMA,
            "target": {"kind": "inline_sequence", "sequence": "ATGTTTAAATAA"},
            "purpose": "synthetic_coding_insert",
            "cds": {"start_0based": 0, "end_0based_exclusive": 12},
            "protein_sequence": "MFK", "genetic_code": 1,
            "protected_intervals": [], "avoid_motifs": [],
            "gc_window": {"window_bp": 6, "min_basis_points": 1666, "max_basis_points": 1667},
            "max_evaluations": 4096, "output_seq_id": "synthetic_window_gc"
        });
        let request: DnaSequenceDesignRequest = serde_json::from_value(wire.clone()).unwrap();
        assert!(request.gc_content.is_none());
        assert_eq!(serde_json::to_value(request).unwrap(), wire);
        for field in ["stride", "include_partial", "strand"] {
            let mut unknown = wire.clone();
            unknown["gc_window"][field] = 1.into();
            assert!(serde_json::from_value::<DnaSequenceDesignRequest>(unknown).is_err());
        }
        let mut fractional = wire;
        fractional["gc_window"]["min_basis_points"] = serde_json::json!(16.66);
        assert!(serde_json::from_value::<DnaSequenceDesignRequest>(fractional).is_err());
    }
}
