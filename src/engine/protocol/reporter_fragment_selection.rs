//! Read-only evidence-guided ROI proposals and bounded, reversible boundary choices.

use super::{
    PromoterReporterPanelCloningStrategyReport, PromoterReporterPanelExtendedBoundaryAudit,
    ReporterVectorValidationReport,
};
use gentle_protocol::{GenomicRegionInterval, GenomicRegionOfInterest, GenomicRegionSet};
use serde::{Deserialize, Serialize};

pub const FRAGMENT_SELECTION_REQUEST_SCHEMA: &str = "gentle.reporter_fragment_selection_request.v1";
pub const FRAGMENT_SELECTION_REPORT_SCHEMA: &str = "gentle.reporter_fragment_selection.v1";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FragmentSelectionLocus {
    pub path: String,
    pub sha256: String,
    pub panel_id: String,
    pub annotation_release: String,
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FragmentSelectionPurpose {
    #[default]
    EndogenousPromoter,
    ResponseElement,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct FragmentSelectionAnchor {
    pub transcript_id: String,
    /// Empty means all annotation/model/called-peak seeds in this anchor's envelope.
    pub seed_evidence_ids: Vec<String>,
    pub required_evidence_ids: Vec<String>,
    /// Optional exact prior ROI, for revisiting a human-selected boundary.
    pub seed_region: Option<GenomicRegionOfInterest>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct FragmentSelectionPolicy {
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub preferred_length_bp: usize,
    pub maximum_length_bp: usize,
    pub flank_bp: usize,
    /// The outer adjustment limit beyond the initial search envelope.
    pub maximum_extension_bp: usize,
    pub maximum_boundary_shift_bp: usize,
    /// Explicit provisional context when no promoter annotation covers the TSS.
    pub fallback_promoter_upstream_bp: usize,
    pub fallback_promoter_downstream_bp: usize,
    pub maximum_candidates: usize,
    pub maximum_evidence_intervals: usize,
    pub purpose: FragmentSelectionPurpose,
}

impl Default for FragmentSelectionPolicy {
    fn default() -> Self {
        Self {
            upstream_bp: 700,
            downstream_bp: 300,
            preferred_length_bp: 700,
            maximum_length_bp: 5_000,
            flank_bp: 20,
            maximum_extension_bp: 200,
            maximum_boundary_shift_bp: 200,
            fallback_promoter_upstream_bp: 100,
            fallback_promoter_downstream_bp: 50,
            maximum_candidates: 128,
            maximum_evidence_intervals: 10_000,
            purpose: FragmentSelectionPurpose::EndogenousPromoter,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FragmentBoundaryReason {
    #[default]
    HumanReview,
    TfbsContext,
    CutrunContext,
    RestrictionSite,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FragmentBoundaryAdjustment {
    pub adjustment_id: String,
    pub transcript_id: String,
    /// Seed evidence ID, or the exact seed_region.region_id.
    pub seed_id: String,
    /// Positive extends outward, negative shortens, always in transcript orientation.
    pub upstream_delta_bp: i64,
    pub downstream_delta_bp: i64,
    pub reason: FragmentBoundaryReason,
    pub explanation: String,
    #[serde(default)]
    pub evidence_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FragmentSelectionVector {
    pub seq_id: String,
    pub catalog_id: String,
    #[serde(default)]
    pub helper_catalog_path: Option<String>,
    /// Suggest padding-only trims around internal MCS-compatible recognition sites.
    #[serde(default)]
    pub suggest_restriction_adjustments: bool,
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FragmentSignalMissingPolicy {
    #[default]
    RequireComplete,
    MissingAsZero,
}

/// A declared descriptive comparison on each FIXED search envelope, never trial inserts.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FragmentSignalComparison {
    pub comparison_id: String,
    pub sample_lane_id: String,
    pub control_lane_id: String,
    pub cell_line: String,
    pub sample_replicate_id: String,
    pub control_replicate_id: String,
    /// Caller-declared compatible units, retained as a declaration, not inferred.
    pub units: String,
    pub minimum_mean_difference: f64,
    pub missing_policy: FragmentSignalMissingPolicy,
}

/// Caller-declared replicate, bound to one available source lane in the locus report.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FragmentPeakLaneBinding {
    pub lane_id: String,
    pub source_sha256: String,
    pub replicate_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case", deny_unknown_fields)]
pub enum FragmentPeakControl {
    Matched { control: FragmentPeakLaneBinding },
    NotUsed { reason: String },
}

/// Explicitly declared BED peak calls, never inferred from generic coverage or labels.
/// Caller settings and replicate identities are provenance declarations, not verified QA.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FragmentCalledPeakSource {
    pub source_id: String,
    pub path: String,
    pub sha256: String,
    pub assembly: String,
    pub chromosome: String,
    pub caller: String,
    pub caller_version: String,
    pub parameters: std::collections::BTreeMap<String, String>,
    pub cell_line: String,
    pub sample: FragmentPeakLaneBinding,
    pub control: FragmentPeakControl,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FragmentSelectionRequest {
    pub schema: String,
    pub locus: FragmentSelectionLocus,
    pub anchors: Vec<FragmentSelectionAnchor>,
    #[serde(default)]
    pub policy: FragmentSelectionPolicy,
    #[serde(default)]
    pub adjustments: Vec<FragmentBoundaryAdjustment>,
    #[serde(default)]
    pub signal_comparisons: Vec<FragmentSignalComparison>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub called_peak_sources: Vec<FragmentCalledPeakSource>,
    #[serde(default)]
    pub vector: Option<FragmentSelectionVector>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FragmentEvidenceKind {
    Annotation,
    ModelSite,
    RawCoverage,
    CalledPeak,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentSelectionEvidence {
    pub evidence_id: String,
    pub kind: FragmentEvidenceKind,
    pub label: String,
    pub interval: GenomicRegionInterval,
    pub source_sha256: String,
    pub source_id: String,
    pub available: bool,
    pub may_seed_boundary: bool,
    pub score: Option<f64>,
    pub statement: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentFixedComparison {
    pub comparison: FragmentSignalComparison,
    pub transcript_id: String,
    pub measurement_window: GenomicRegionInterval,
    pub sample_source_sha256: Option<String>,
    pub control_source_sha256: Option<String>,
    pub sample_mean: Option<f64>,
    pub control_mean: Option<f64>,
    pub mean_difference: Option<f64>,
    pub passes_descriptive_rule: Option<bool>,
    pub status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentSelectionAnchorResult {
    pub transcript_id: String,
    pub tss: GenomicRegionInterval,
    pub search_envelope: GenomicRegionInterval,
    pub opposite_strand_transcript_ids: Vec<String>,
    pub gene_structure: Option<PromoterReporterPanelExtendedBoundaryAudit>,
    pub findings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentRestrictionSite {
    pub enzyme: String,
    /// Recognition footprint, NOT an inferred ligation junction.
    pub interval: GenomicRegionInterval,
    pub local_top_cut_0based: Option<usize>,
    pub local_bottom_cut_0based: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentCandidateRanking {
    pub preserves_required_context: bool,
    pub descriptive_enrichment_supported: bool,
    pub retains_annotation: bool,
    pub retains_model_site: bool,
    /// Presence only; shares sample provenance with coverage, not an independent vote.
    #[serde(default)]
    pub retains_called_peak: bool,
    pub bisected_feature_count: usize,
    pub fits_preferred_length: bool,
    pub length_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReporterFragmentCandidate {
    pub candidate_id: String,
    pub transcript_id: String,
    pub member_transcript_ids: Vec<String>,
    pub seed_id: String,
    pub variant: String,
    pub parent_candidate_id: Option<String>,
    pub region: GenomicRegionOfInterest,
    pub sequence_5prime_to_3prime: String,
    pub sequence_sha256: String,
    pub retained_evidence_ids: Vec<String>,
    pub excluded_evidence_ids: Vec<String>,
    pub bisected_evidence_ids: Vec<String>,
    pub ranking: FragmentCandidateRanking,
    pub blockers: Vec<String>,
    pub reasons: Vec<String>,
    pub cloning: Option<PromoterReporterPanelCloningStrategyReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentSelectionVectorResult {
    pub validation: ReporterVectorValidationReport,
    pub sequence_sha256: String,
    pub catalog_sha256: String,
    pub mcs_start_0based: usize,
    pub mcs_end_0based_exclusive: usize,
    pub source_sites: Vec<FragmentRestrictionSite>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentSelectionReport {
    pub schema: String,
    pub request: FragmentSelectionRequest,
    pub request_sha256: String,
    pub proposal_sha256: String,
    pub source_sequence_sha256: String,
    pub evidence: Vec<FragmentSelectionEvidence>,
    pub anchors: Vec<FragmentSelectionAnchorResult>,
    pub fixed_comparisons: Vec<FragmentFixedComparison>,
    pub candidates: Vec<ReporterFragmentCandidate>,
    /// Eligible candidates only. Import explicitly through the existing region operation.
    pub proposed_region_set: GenomicRegionSet,
    pub vector: Option<FragmentSelectionVectorResult>,
    pub findings: Vec<String>,
    pub non_claims: Vec<String>,
}
