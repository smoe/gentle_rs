//! Gene-agnostic transcript capture discovery, not an approved assay or order.

use super::{
    SequenceGenomeAnchorSummary, TranscriptAssayCdnaSynthesis, TranscriptAssayCoveragePolicy,
    TranscriptAssayCoverageResolution, TranscriptAssayCoverageUniverse,
};
use serde::{Deserialize, Serialize};

pub const TRANSCRIPT_CAPTURE_REQUEST_SCHEMA: &str = "gentle.transcript_capture_pool_request.v1";
pub const TRANSCRIPT_CAPTURE_REPORT_SCHEMA: &str = "gentle.transcript_capture_pool.v1";

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum TranscriptCaptureRole {
    /// Order the transcript-sense sequence; it binds antisense first-strand cDNA.
    SenseForward,
    /// Order the reverse complement; it binds RNA or sense cDNA.
    AntisenseReverse,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum TranscriptCaptureWindow {
    /// No automatic CDS or transcript-relative fallback when CDS is unknown.
    FivePrimeUtr,
    TranscriptRange {
        start_0based: usize,
        end_0based_exclusive: usize,
    },
    TerminalExonStart {
        search_window_bp: usize,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranscriptCaptureSource {
    pub seq_id: String,
    pub source_feature_id: usize,
    #[serde(default)]
    pub coverage_universe: TranscriptAssayCoverageUniverse,
    #[serde(default)]
    pub annotation_release: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranscriptCaptureTarget {
    pub target_id: String,
    pub sources: Vec<TranscriptCaptureSource>,
    pub role: TranscriptCaptureRole,
    pub window: TranscriptCaptureWindow,
    pub max_primers: usize,
    /// Equal nonempty keys permit, but do not require, one oligo across targets.
    #[serde(default)]
    pub sharing_group: Option<String>,
    pub stage_ids: Vec<String>,
    #[serde(default)]
    pub tail_5prime: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranscriptCaptureFixedOligo {
    pub oligo_id: String,
    /// Complete ordered sequence, including any adapter; IUPAC is supported.
    pub full_oligo_5_to_3: String,
    pub stage_ids: Vec<String>,
    pub provenance: String,
    pub reorder_same_sequence: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranscriptCaptureTmRange {
    pub min_c: f64,
    pub max_c: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct TranscriptCaptureSearchPolicy {
    pub min_length_bp: usize,
    pub max_length_bp: usize,
    pub max_candidates_per_target: usize,
    pub beam_width: usize,
    /// Optional filter under the explicitly reported shared Tm model, not Ta.
    pub tm_range: Option<TranscriptCaptureTmRange>,
    /// Report exact A-runs of at least this length as possible internal priming.
    pub internal_a_run_min_bp: usize,
}

impl Default for TranscriptCaptureSearchPolicy {
    fn default() -> Self {
        Self {
            min_length_bp: 20,
            max_length_bp: 24,
            max_candidates_per_target: 12,
            beam_width: 128,
            tm_range: None,
            internal_a_run_min_bp: 12,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranscriptCapturePoolRequest {
    pub schema: String,
    pub report_id: String,
    pub targets: Vec<TranscriptCaptureTarget>,
    #[serde(default)]
    pub fixed_oligos: Vec<TranscriptCaptureFixedOligo>,
    #[serde(default)]
    pub coverage_policy: TranscriptAssayCoveragePolicy,
    pub cdna_synthesis: TranscriptAssayCdnaSynthesis,
    #[serde(default)]
    pub search: TranscriptCaptureSearchPolicy,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptCaptureSourceBinding {
    pub target_id: String,
    pub source: TranscriptCaptureSource,
    /// Hash of the complete loaded DNA record, including annotation.
    pub source_record_sha256: String,
    pub genome_anchor: Option<SequenceGenomeAnchorSummary>,
    pub coverage: TranscriptAssayCoverageResolution,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptCaptureMember {
    pub member_id: String,
    pub target_id: String,
    pub seq_id: String,
    pub transcript_id: String,
    pub transcript_feature_id: usize,
    pub source_strand: String,
    pub cdna_length_bp: usize,
    pub cdna_sha256: String,
    pub five_prime_utr_length_bp: Option<usize>,
    pub search_range_0based: Option<(usize, usize)>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptCaptureBinding {
    pub member_id: String,
    pub transcript_start_0based: usize,
    pub transcript_end_0based_exclusive: usize,
    pub source_ranges_0based: Vec<(usize, usize)>,
    pub within_requested_window: bool,
    /// False for an occurrence in a target that did not authorize sharing.
    pub permitted_target: bool,
    pub upstream_bases_omitted: usize,
    pub downstream_bases_omitted: usize,
    /// Annotation-derived retained interval, excluding synthetic tails/poly(A).
    pub retained_length_bp: usize,
    pub retained_sequence_sha256: String,
    /// False means identical ambiguity symbols are not proof of equivalence.
    pub retained_sequence_canonical: bool,
    pub internal_a_runs_0based: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptCaptureCandidate {
    pub candidate_id: String,
    pub role: TranscriptCaptureRole,
    pub permitted_target_ids: Vec<String>,
    pub stage_ids: Vec<String>,
    pub annealing_5_to_3: String,
    pub full_oligo_5_to_3: String,
    pub tm_c: f64,
    pub self_3prime_run_bp: usize,
    pub self_complementary_run_bp: usize,
    pub homopolymer_run_bp: usize,
    pub max_fixed_3prime_run_bp: usize,
    pub max_fixed_complementary_run_bp: usize,
    pub covered_member_ids: Vec<String>,
    pub bindings: Vec<TranscriptCaptureBinding>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptCaptureInteraction {
    pub left_oligo_id: String,
    pub right_oligo_id: String,
    pub shared_stage_ids: Vec<String>,
    pub max_complementary_run_bp: usize,
    pub max_3prime_complementary_run_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptCaptureEquivalence {
    pub retained_sequence_sha256: String,
    pub retained_length_bp: usize,
    /// Candidate/member/position instances; one transcript can yield several.
    pub binding_instances: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptCapturePoolReport {
    pub schema: String,
    pub report_id: String,
    pub op_id: String,
    pub run_id: String,
    pub request_sha256: String,
    pub request: TranscriptCapturePoolRequest,
    pub sources: Vec<TranscriptCaptureSourceBinding>,
    pub members: Vec<TranscriptCaptureMember>,
    pub candidates: Vec<TranscriptCaptureCandidate>,
    /// Best found within the declared bounds; never approval or an order.
    pub proposed_candidate_ids: Vec<String>,
    pub uncovered_member_ids: Vec<String>,
    pub coverage_satisfied: bool,
    pub distinct_candidates_evaluated: usize,
    pub ambiguous_windows_skipped: usize,
    pub candidate_retention_truncated: bool,
    pub beam_truncated: bool,
    pub pool_states_evaluated: usize,
    pub interactions: Vec<TranscriptCaptureInteraction>,
    pub captured_equivalence_groups: Vec<TranscriptCaptureEquivalence>,
    pub retained_length_range_bp: Option<(usize, usize)>,
    pub tm_model: String,
    pub specificity_status: String,
    pub warnings: Vec<String>,
}
