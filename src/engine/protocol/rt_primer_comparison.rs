//! Isoform-aware RT candidate comparison; annotation context and compatibility are independent evidence.

use super::{
    AssessRtPrimerCandidatesRequest, RtPrimerAssessmentReport, TranscriptCaptureSource,
    TranscriptCaptureSourceBinding,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const RT_PRIMER_COMPARISON_REQUEST_SCHEMA: &str = "gentle.rt_primer_comparison_request.v1";
pub const RT_PRIMER_COMPARISON_SCHEMA: &str = "gentle.rt_primer_comparison.v1";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RtPrimerComparisonSource {
    pub source: TranscriptCaptureSource,
    /// Optional guard against a changed loaded DNA/annotation record, not just changed bases.
    pub expected_source_record_sha256: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RtPrimerComparisonBudgets {
    pub max_transcript_base_comparisons: u64,
    pub max_bindings: usize,
}
impl Default for RtPrimerComparisonBudgets {
    fn default() -> Self {
        Self {
            max_transcript_base_comparisons: 100_000_000,
            max_bindings: 10_000,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CompareRtPrimerCandidatesRequest {
    pub schema: String,
    /// Composed exactly once through the existing prospective assessment operation.
    pub assessment: AssessRtPrimerCandidatesRequest,
    pub sources: Vec<RtPrimerComparisonSource>,
    #[serde(default)]
    pub budgets: RtPrimerComparisonBudgets,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum RtPrimerBindingRegion {
    Cds,
    FivePrimeUtr,
    ThreePrimeUtr,
    Mixed,
    Unknown,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RtPrimerSitePresence {
    Present,
    Absent,
    Unresolved,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerComparisonMember {
    pub member_id: String,
    pub source_index: usize,
    pub seq_id: String,
    pub transcript_id: String,
    pub transcript_feature_id: usize,
    pub strand: String,
    pub cdna_length_bp: usize,
    pub cdna_sha256: String,
    /// Annotated CDS union in biological cDNA coordinates; empty means unknown, not all-UTR.
    pub cds_ranges_0based: Vec<(usize, usize)>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerClassifiedBinding {
    pub transcript_start_0based: usize,
    pub transcript_end_0based_exclusive: usize,
    pub source_ranges_0based: Vec<(usize, usize)>,
    pub region: RtPrimerBindingRegion,
    pub cds_bases: usize,
    pub five_prime_utr_bases: usize,
    pub three_prime_utr_bases: usize,
    pub unknown_bases: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerIsoformAssessment {
    pub primer_id: String,
    pub member_id: String,
    pub presence: RtPrimerSitePresence,
    /// Compatible ambiguous windows, distinct from proven canonical exact sites.
    pub possible_ambiguous_sites: usize,
    pub bindings: Vec<RtPrimerClassifiedBinding>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerComparisonMetrics {
    pub length_bp: usize,
    pub gc_fraction: f64,
    pub tm_c: f64,
    pub self_complementary_run_bp: usize,
    pub self_3prime_complementary_run_bp: usize,
    pub longest_homopolymer_run_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerCandidateComparison {
    pub primer_id: String,
    pub target_tail: RtPrimerComparisonMetrics,
    pub full_oligo: RtPrimerComparisonMetrics,
    pub total_members: usize,
    pub members_with_exact_site: usize,
    pub members_with_possible_only_site: usize,
    pub members_without_site: usize,
    /// One member per region regardless of repeat-site count. Different regions may overlap.
    pub region_member_counts: BTreeMap<RtPrimerBindingRegion, usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerComparisonReport {
    pub schema: String,
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: CompareRtPrimerCandidatesRequest,
    pub assessment: Box<RtPrimerAssessmentReport>,
    pub sources: Vec<TranscriptCaptureSourceBinding>,
    pub members: Vec<RtPrimerComparisonMember>,
    pub isoforms: Vec<RtPrimerIsoformAssessment>,
    /// Input candidate order; no implicit ranking or order-readiness decision.
    pub candidates: Vec<RtPrimerCandidateComparison>,
    pub transcript_base_comparisons: u64,
    pub tm_model: String,
    pub secondary_structure_status: String,
    pub warnings: Vec<String>,
}
