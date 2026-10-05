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

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum DnaDesignTarget {
    InlineSequence { sequence: String },
    LoadedSequence { seq_id: String },
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
    pub max_evaluations: u64,
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
    pub optimization_complete: bool,
    pub minimum_edits_proven: bool,
    pub approval_digest: Option<String>,
    pub nonclaims: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DnaSequenceDesignReceipt {
    pub schema: String,
    pub approval_digest: String,
    pub created_seq_id: String,
    pub output_sha256: String,
    pub nonclaims: Vec<String>,
}
