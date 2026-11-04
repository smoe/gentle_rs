//! Constrained RT candidate generation for explicit transcript groups; output remains advisory.

use super::*;

pub const RT_PRIMER_REDESIGN_REQUEST_SCHEMA: &str = "gentle.rt_primer_redesign_request.v1";
pub const RT_PRIMER_REDESIGN_SCHEMA: &str = "gentle.rt_primer_redesign.v1";

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Required group-wide binding rule; UTR mode also excludes compatible off-region alternatives.
pub enum RtPrimerRedesignMode {
    SharedTranscripts,
    ThreePrimeUtr,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct RtPrimerRedesignMember {
    pub seq_id: String,
    pub transcript_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RtPrimerRedesignAnchor {
    pub member: RtPrimerRedesignMember,
    /// Search interval in mature biological 5′–3′ cDNA, including possible exon junctions.
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Inclusive tail ranges and tail/full-oligo diagnostic ceilings; no experimental approval.
pub struct RtPrimerRedesignConstraints {
    pub min_tail_length_bp: usize,
    pub max_tail_length_bp: usize,
    pub min_tail_gc_fraction: f64,
    pub max_tail_gc_fraction: f64,
    pub min_tail_tm_c: f64,
    pub max_tail_tm_c: f64,
    /// Diagnostics constrain the complete handle-bearing oligo as well as its tail.
    pub max_homopolymer_run_bp: usize,
    pub max_self_complementary_run_bp: usize,
    pub max_self_3prime_complementary_run_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RtPrimerRedesignSeedLimit {
    pub reference_id: String,
    pub seed_length: usize,
    pub max_mismatches: usize,
    /// All physical windows, including intended sites; no inferred intended-site exclusions.
    pub max_compatible_windows: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RtPrimerRedesignBudgets {
    pub max_enumerated_windows: usize,
    /// Covers canonical intake, heuristic work and group-member matching, before comparison.
    pub max_generation_work: u64,
    /// First qualifying distinct tails in start-coordinate then length order; never a best-score claim.
    pub max_candidates: usize,
}
impl Default for RtPrimerRedesignBudgets {
    fn default() -> Self {
        Self {
            max_enumerated_windows: 50_000,
            max_generation_work: 100_000_000,
            max_candidates: 16,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Explicit group, anchor search scope and background inputs for one bounded advisory design.
pub struct RedesignRtPrimersRequest {
    pub schema: String,
    pub group_id: String,
    pub mode: RtPrimerRedesignMode,
    /// Each source uses an explicit transcript universe. May include comparison-only members.
    pub sources: Vec<RtPrimerComparisonSource>,
    pub required_members: Vec<RtPrimerRedesignMember>,
    pub anchor: RtPrimerRedesignAnchor,
    pub common_handle: String,
    pub constraints: RtPrimerRedesignConstraints,
    pub references: Vec<RtPrimerBackground>,
    #[serde(default)]
    pub seed_limits: Vec<RtPrimerRedesignSeedLimit>,
    #[serde(default)]
    pub budgets: RtPrimerRedesignBudgets,
    #[serde(default)]
    pub assessment_budgets: RtPrimerAssessmentBudgets,
    #[serde(default)]
    pub comparison_budgets: RtPrimerComparisonBudgets,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// A declared constraint outcome; unresolved evidence must not become a clean pass.
pub enum RtPrimerRedesignConstraintStatus {
    Passed,
    Failed,
    Unresolved,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerRedesignSeedCheck {
    pub limit: RtPrimerRedesignSeedLimit,
    pub compatible_windows: Option<u64>,
    pub ambiguous_windows: Option<u64>,
    pub status: RtPrimerRedesignConstraintStatus,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerRedesignCandidate {
    pub primer_id: String,
    pub target_tail: String,
    pub anchor_start_0based: usize,
    pub anchor_end_0based_exclusive: usize,
    /// Only declared seed limits, never an experimental or order-readiness decision.
    pub seed_constraint_status: RtPrimerRedesignConstraintStatus,
    pub seed_checks: Vec<RtPrimerRedesignSeedCheck>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Exhaustive admitted generation counts and a capped shortlist with unchanged comparison evidence.
pub struct RtPrimerRedesignReport {
    pub schema: String,
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: RedesignRtPrimersRequest,
    pub sources: Vec<TranscriptCaptureSourceBinding>,
    pub members: Vec<RtPrimerComparisonMember>,
    pub required_member_ids: Vec<String>,
    pub enumerated_windows: usize,
    pub skipped_ambiguous_windows: usize,
    pub duplicate_tail_windows: usize,
    pub metric_rejected_tails: usize,
    pub group_rejected_tails: usize,
    pub qualifying_distinct_tails: usize,
    pub generation_work: u64,
    /// Additional qualifying tails were not background-assessed, never silently rejected.
    pub shortlist_truncated: bool,
    pub candidates: Vec<RtPrimerRedesignCandidate>,
    /// None only for an exhaustive empty generation result; backgrounds then remain unopened.
    pub comparison: Option<Box<RtPrimerComparisonReport>>,
    pub warnings: Vec<String>,
}
