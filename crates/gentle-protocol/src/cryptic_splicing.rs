//! Portable contracts for deterministic cryptic-splicing candidate screens.
//!
//! These records deliberately distinguish sequence-structural candidates from
//! annotated introns, model scores, and observed RNA junctions. The v1 screen
//! is computed on demand and does not persist candidate rows in project state.

use serde::{Deserialize, Serialize};

pub const CRYPTIC_SPLICING_SCREEN_SCHEMA: &str = "gentle.cryptic_splicing_screen.v1";

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CrypticSplicingStrand {
    #[default]
    Forward,
    Reverse,
}

impl CrypticSplicingStrand {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CrypticSplicingEvidenceClass {
    Annotated,
    #[default]
    StructuralCandidate,
    ModelScored,
    ObservedJunction,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CrypticSplicingModelStatus {
    #[default]
    Absent,
    Present,
    Failed,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CrypticSplicingSignalStatus {
    Detected,
    NotDetected,
    #[default]
    NotEvaluable,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CrypticSplicingSiteKind {
    #[default]
    Donor,
    Acceptor,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingSpan {
    pub start_1based: usize,
    pub end_1based: usize,
}

fn default_min_pseudo_intron_bp() -> usize {
    50
}

fn default_max_pseudo_intron_bp() -> usize {
    5_000
}

fn default_max_candidate_pairs() -> usize {
    1_000
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingScreenRequest {
    pub seq_id: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: CrypticSplicingStrand,
    pub insert_span: Option<CrypticSplicingSpan>,
    /// Zero-based feature index of a CDS on `seq_id`.
    pub cds_feature_id: Option<usize>,
    #[serde(default = "default_min_pseudo_intron_bp")]
    pub min_pseudo_intron_bp: usize,
    #[serde(default = "default_max_pseudo_intron_bp")]
    pub max_pseudo_intron_bp: usize,
    #[serde(default = "default_max_candidate_pairs")]
    pub max_candidate_pairs: usize,
}

impl Default for CrypticSplicingScreenRequest {
    fn default() -> Self {
        Self {
            seq_id: String::new(),
            start_1based: 1,
            end_1based: 0,
            strand: CrypticSplicingStrand::Forward,
            insert_span: None,
            cds_feature_id: None,
            min_pseudo_intron_bp: default_min_pseudo_intron_bp(),
            max_pseudo_intron_bp: default_max_pseudo_intron_bp(),
            max_candidate_pairs: default_max_candidate_pairs(),
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingGenomicProvenanceRow {
    pub scanned_start_1based: usize,
    pub scanned_end_1based: usize,
    pub source_seq_id: Option<String>,
    pub source_feature_id: Option<usize>,
    pub genomic_start_1based: Option<usize>,
    pub genomic_end_1based: Option<usize>,
    pub strand: Option<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingSiteRow {
    pub site_id: String,
    pub evidence_class: CrypticSplicingEvidenceClass,
    pub kind: CrypticSplicingSiteKind,
    pub motif_2bp: String,
    pub scanned_position_1based: usize,
    pub source_position_1based: usize,
    pub model_window: Option<String>,
    pub model_status: CrypticSplicingModelStatus,
    pub maxent_score: Option<f64>,
    pub not_evaluable_reason: Option<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingBranchpointSignal {
    pub status: CrypticSplicingSignalStatus,
    pub scanned_position_1based: Option<usize>,
    pub source_position_1based: Option<usize>,
    pub motif: Option<String>,
    pub heuristic_score: Option<f32>,
    pub annotation: Option<String>,
    pub not_evaluable_reason: Option<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingPolypyrimidineSignal {
    pub status: CrypticSplicingSignalStatus,
    pub scanned_start_1based: Option<usize>,
    pub scanned_end_1based: Option<usize>,
    pub source_start_1based: Option<usize>,
    pub source_end_1based: Option<usize>,
    pub pyrimidine_fraction: Option<f32>,
    pub annotation: Option<String>,
    pub not_evaluable_reason: Option<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingCdsConsequence {
    pub status: String,
    pub cds_feature_id: Option<usize>,
    pub removed_coding_bp: usize,
    pub frame_delta_mod_3: Option<usize>,
    pub affected_coding_interval_1based: Option<CrypticSplicingSpan>,
    /// Reserved for a future projection against a specific protein model.
    pub affected_aa_interval: Option<CrypticSplicingSpan>,
    pub native_stop_removed: Option<bool>,
    pub first_altered_aa_position: Option<usize>,
    pub predicted_protein_length_aa: Option<usize>,
    pub interpretation: String,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingCandidateRow {
    pub candidate_id: String,
    pub evidence_class: CrypticSplicingEvidenceClass,
    pub donor_site_id: String,
    pub acceptor_site_id: String,
    pub donor_scanned_position_1based: usize,
    pub acceptor_scanned_position_1based: usize,
    pub donor_source_position_1based: usize,
    pub acceptor_source_position_1based: usize,
    pub removed_source_start_1based: usize,
    pub removed_source_end_1based: usize,
    pub pseudo_intron_length_bp: usize,
    pub paired_motif_signature: String,
    pub motif_class: String,
    pub boundary_class: String,
    pub branchpoint: CrypticSplicingBranchpointSignal,
    pub polypyrimidine_tract: CrypticSplicingPolypyrimidineSignal,
    pub model_status: CrypticSplicingModelStatus,
    pub donor_maxent_score: Option<f64>,
    pub acceptor_maxent_score: Option<f64>,
    pub prioritization_heuristic: Option<f64>,
    pub cds_consequence: Option<CrypticSplicingCdsConsequence>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingBudgetSummary {
    pub donor_site_count: usize,
    pub acceptor_site_count: usize,
    pub admissible_pair_count: usize,
    pub evaluated_pair_count: usize,
    pub truncated: bool,
    pub truncation_rule: Option<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingModelProvenance {
    pub status: CrypticSplicingModelStatus,
    pub model_kind: Option<String>,
    pub resource_path: Option<String>,
    pub resource_sha256: Option<String>,
    pub note: String,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CrypticSplicingScreenView {
    pub schema: String,
    pub request: CrypticSplicingScreenRequest,
    pub request_sha256: String,
    pub source_digest: String,
    pub source_length_bp: usize,
    pub scanned_sequence_length_bp: usize,
    pub coordinate_space: String,
    pub coordinate_convention: String,
    pub genomic_provenance: Vec<CrypticSplicingGenomicProvenanceRow>,
    pub model: CrypticSplicingModelProvenance,
    pub donor_sites: Vec<CrypticSplicingSiteRow>,
    pub acceptor_sites: Vec<CrypticSplicingSiteRow>,
    pub candidates: Vec<CrypticSplicingCandidateRow>,
    pub budget: CrypticSplicingBudgetSummary,
    pub warnings: Vec<String>,
}
