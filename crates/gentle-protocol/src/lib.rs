//! Shared machine-readable GENtle contracts.
//!
//! This crate is intentionally scaffolded first in the workspace split so
//! stable cross-adapter payloads can move here before execution or GUI code.
//! The first extracted slice is intentionally small: stable identifier aliases,
//! shared analysis enums, and the portable engine error payload.

use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, error::Error, fmt};

/// Stable identifier for one sequence entry stored in project state.
pub type SeqId = String;
/// Stable identifier for one executed operation journal row.
pub type OpId = String;
/// Caller-supplied identifier that groups operations into one workflow/run.
pub type RunId = String;
/// Stable identifier for one lineage graph node.
pub type NodeId = String;
/// Stable identifier for one wet-lab-style container record.
pub type ContainerId = String;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Dotplot comparison mode for dotplot computation/export.
pub enum DotplotMode {
    #[default]
    SelfForward,
    SelfReverseComplement,
    PairForward,
    PairReverseComplement,
}

impl DotplotMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SelfForward => "self_forward",
            Self::SelfReverseComplement => "self_reverse_complement",
            Self::PairForward => "pair_forward",
            Self::PairReverseComplement => "pair_reverse_complement",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Pairwise alignment mode for sequence and confirmation alignments.
pub enum PairwiseAlignmentMode {
    #[default]
    Global,
    Local,
}

impl PairwiseAlignmentMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Global => "global",
            Self::Local => "local",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Flexibility score model for flexibility-track computation.
pub enum FlexibilityModel {
    #[default]
    AtRichness,
    AtSkew,
}

impl FlexibilityModel {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AtRichness => "at_richness",
            Self::AtSkew => "at_skew",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadInputFormat {
    #[default]
    Fasta,
}

impl RnaReadInputFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Fasta => "fasta",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadInterpretationProfile {
    #[default]
    NanoporeCdnaV1,
    ShortReadV1,
    TransposonV1,
}

impl RnaReadInterpretationProfile {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NanoporeCdnaV1 => "nanopore_cdna_v1",
            Self::ShortReadV1 => "short_read_v1",
            Self::TransposonV1 => "transposon_v1",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadReportMode {
    #[default]
    Full,
    SeedPassedOnly,
}

impl RnaReadReportMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Full => "full",
            Self::SeedPassedOnly => "seed_passed_only",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadOriginMode {
    #[default]
    SingleGene,
    MultiGeneSparse,
}

impl RnaReadOriginMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SingleGene => "single_gene",
            Self::MultiGeneSparse => "multi_gene_sparse",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadHitSelection {
    All,
    SeedPassed,
    #[default]
    Aligned,
}

impl RnaReadHitSelection {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::SeedPassed => "seed_passed",
            Self::Aligned => "aligned",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadScoreDensityScale {
    Linear,
    #[default]
    Log,
}

impl RnaReadScoreDensityScale {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::Log => "log",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadScoreDensityVariant {
    #[default]
    AllScored,
    CompositeSeedGate,
    RetainedReplayCurrentControls,
}

impl RnaReadScoreDensityVariant {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AllScored => "all_scored",
            Self::CompositeSeedGate => "composite_seed_gate",
            Self::RetainedReplayCurrentControls => "retained_replay_current_controls",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadAlignmentMode {
    #[default]
    Local,
    Semiglobal,
}

impl RnaReadAlignmentMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Local => "local",
            Self::Semiglobal => "semiglobal",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadOriginClass {
    TargetCoherent,
    TargetPartialLocalBlock,
    RoiSameStrandLocalBlock,
    RoiReverseStrandLocalBlock,
    TpFamilyAmbiguous,
    #[default]
    BackgroundLikely,
}

impl RnaReadOriginClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::TargetCoherent => "target_coherent",
            Self::TargetPartialLocalBlock => "target_partial_local_block",
            Self::RoiSameStrandLocalBlock => "roi_same_strand_local_block",
            Self::RoiReverseStrandLocalBlock => "roi_reverse_strand_local_block",
            Self::TpFamilyAmbiguous => "tp_family_ambiguous",
            Self::BackgroundLikely => "background_likely",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadAlignmentBackend {
    #[default]
    Banded,
    DenseFallback,
}

impl RnaReadAlignmentBackend {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Banded => "banded",
            Self::DenseFallback => "dense_fallback",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadAlignmentEffect {
    #[default]
    ConfirmedAssignment,
    ReassignedTranscript,
    AlignedWithoutPhase1Assignment,
}

impl RnaReadAlignmentEffect {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ConfirmedAssignment => "confirmed_assignment",
            Self::ReassignedTranscript => "reassigned_transcript",
            Self::AlignedWithoutPhase1Assignment => "aligned_without_phase1_assignment",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadAlignmentInspectionEffectFilter {
    #[default]
    AllAligned,
    ConfirmedOnly,
    DisagreementOnly,
    ReassignedOnly,
    NoPhase1Only,
    SelectedOnly,
}

impl RnaReadAlignmentInspectionEffectFilter {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AllAligned => "all_aligned",
            Self::ConfirmedOnly => "confirmed_only",
            Self::DisagreementOnly => "disagreement_only",
            Self::ReassignedOnly => "reassigned_only",
            Self::NoPhase1Only => "no_phase1_only",
            Self::SelectedOnly => "selected_only",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadAlignmentInspectionSortKey {
    #[default]
    Rank,
    Identity,
    Coverage,
    Score,
}

impl RnaReadAlignmentInspectionSortKey {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Rank => "rank",
            Self::Identity => "identity",
            Self::Coverage => "coverage",
            Self::Score => "score",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
/// Stable high-level error class for engine operation failures.
pub enum ErrorCode {
    InvalidInput,
    NotFound,
    Unsupported,
    Io,
    Internal,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Shared structured error payload returned by engine operations/adapters.
pub struct EngineError {
    pub code: ErrorCode,
    pub message: String,
}

impl fmt::Display for EngineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}: {}", self.code, self.message)
    }
}

impl Error for EngineError {}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotMatchPoint {
    pub x_0based: usize,
    pub y_0based: usize,
    pub mismatches: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotBoxplotBin {
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub hit_count: usize,
    pub min_reference_0based: Option<usize>,
    pub q1_reference_0based: Option<usize>,
    pub median_reference_0based: Option<usize>,
    pub q3_reference_0based: Option<usize>,
    pub max_reference_0based: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotQuerySeries {
    pub series_id: String,
    pub seq_id: String,
    pub label: String,
    pub color_rgb: [u8; 3],
    #[serde(default)]
    pub mode: DotplotMode,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub point_count: usize,
    pub points: Vec<DotplotMatchPoint>,
    pub boxplot_bin_count: usize,
    pub boxplot_bins: Vec<DotplotBoxplotBin>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotReferenceAnnotationInterval {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub label: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotReferenceAnnotationTrack {
    pub seq_id: String,
    pub label: String,
    pub interval_count: usize,
    pub intervals: Vec<DotplotReferenceAnnotationInterval>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotOverlayQuerySpec {
    pub seq_id: String,
    pub label: String,
    #[serde(default)]
    pub span_start_0based: Option<usize>,
    #[serde(default)]
    pub span_end_0based: Option<usize>,
    #[serde(default)]
    pub mode: DotplotMode,
    #[serde(default)]
    pub color_rgb: Option<[u8; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotView {
    pub schema: String,
    pub dotplot_id: String,
    pub owner_seq_id: String,
    pub seq_id: String,
    pub reference_seq_id: Option<String>,
    pub generated_at_unix_ms: u128,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub reference_span_start_0based: usize,
    pub reference_span_end_0based: usize,
    pub mode: DotplotMode,
    pub word_size: usize,
    pub step_bp: usize,
    pub max_mismatches: usize,
    pub tile_bp: Option<usize>,
    pub point_count: usize,
    pub points: Vec<DotplotMatchPoint>,
    pub boxplot_bin_count: usize,
    pub boxplot_bins: Vec<DotplotBoxplotBin>,
    pub series_count: usize,
    pub query_series: Vec<DotplotQuerySeries>,
    #[serde(default)]
    pub reference_annotation: Option<DotplotReferenceAnnotationTrack>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DotplotViewSummary {
    pub dotplot_id: String,
    pub owner_seq_id: String,
    pub seq_id: String,
    pub reference_seq_id: Option<String>,
    pub generated_at_unix_ms: u128,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub reference_span_start_0based: usize,
    pub reference_span_end_0based: usize,
    pub mode: DotplotMode,
    pub word_size: usize,
    pub step_bp: usize,
    pub max_mismatches: usize,
    pub point_count: usize,
    pub series_count: usize,
}

impl DotplotView {
    pub fn normalize_v3_defaults(&mut self) {
        if self.owner_seq_id.trim().is_empty() {
            self.owner_seq_id = self.seq_id.clone();
        }
        if self.query_series.is_empty() {
            let label = if self.seq_id.trim().is_empty() {
                "<query>".to_string()
            } else {
                self.seq_id.clone()
            };
            self.query_series.push(DotplotQuerySeries {
                series_id: if self.dotplot_id.trim().is_empty() {
                    "series_1".to_string()
                } else {
                    format!("{}_series_1", self.dotplot_id)
                },
                seq_id: self.seq_id.clone(),
                label,
                color_rgb: [29, 78, 216],
                mode: self.mode,
                span_start_0based: self.span_start_0based,
                span_end_0based: self.span_end_0based,
                point_count: if self.point_count == 0 {
                    self.points.len()
                } else {
                    self.point_count
                },
                points: self.points.clone(),
                boxplot_bin_count: if self.boxplot_bin_count == 0 {
                    self.boxplot_bins.len()
                } else {
                    self.boxplot_bin_count
                },
                boxplot_bins: self.boxplot_bins.clone(),
            });
        }
        self.series_count = self.query_series.len();
        if let Some(primary) = self.query_series.first() {
            if self.seq_id.trim().is_empty() {
                self.seq_id = primary.seq_id.clone();
            }
            self.mode = primary.mode;
            self.span_start_0based = primary.span_start_0based;
            self.span_end_0based = primary.span_end_0based;
            self.point_count = primary.point_count.max(primary.points.len());
            self.points = primary.points.clone();
            self.boxplot_bin_count = primary.boxplot_bin_count.max(primary.boxplot_bins.len());
            self.boxplot_bins = primary.boxplot_bins.clone();
        }
        if let Some(annotation) = self.reference_annotation.as_mut() {
            annotation.interval_count = annotation.intervals.len();
        }
    }

    pub fn primary_series(&self) -> Option<&DotplotQuerySeries> {
        self.query_series.first()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceAlignmentReport {
    pub schema: String,
    pub mode: PairwiseAlignmentMode,
    pub query_seq_id: String,
    pub target_seq_id: String,
    pub query_span_start_0based: usize,
    pub query_span_end_0based: usize,
    pub target_span_start_0based: usize,
    pub target_span_end_0based: usize,
    pub aligned_query_start_0based: usize,
    pub aligned_query_end_0based_exclusive: usize,
    pub aligned_target_start_0based: usize,
    pub aligned_target_end_0based_exclusive: usize,
    pub score: i32,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub aligned_columns: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
    pub target_coverage_fraction: f64,
    pub cigar: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequenceFeatureRangeRelation {
    #[default]
    Overlap,
    Within,
    Contains,
}

impl SequenceFeatureRangeRelation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Overlap => "overlap",
            Self::Within => "within",
            Self::Contains => "contains",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequenceFeatureStrandFilter {
    #[default]
    Any,
    Forward,
    Reverse,
}

impl SequenceFeatureStrandFilter {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Any => "any",
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequenceFeatureSortBy {
    FeatureId,
    #[default]
    Start,
    End,
    Kind,
    Length,
}

impl SequenceFeatureSortBy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FeatureId => "feature_id",
            Self::Start => "start",
            Self::End => "end",
            Self::Kind => "kind",
            Self::Length => "length",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQualifierFilter {
    pub key: String,
    pub value_contains: Option<String>,
    pub value_regex: Option<String>,
    pub case_sensitive: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQuery {
    pub seq_id: SeqId,
    pub include_source: bool,
    pub include_qualifiers: bool,
    pub kind_in: Vec<String>,
    pub kind_not_in: Vec<String>,
    pub start_0based: Option<usize>,
    pub end_0based_exclusive: Option<usize>,
    pub range_relation: SequenceFeatureRangeRelation,
    pub strand: SequenceFeatureStrandFilter,
    pub label_contains: Option<String>,
    pub label_regex: Option<String>,
    pub qualifier_filters: Vec<SequenceFeatureQualifierFilter>,
    pub min_len_bp: Option<usize>,
    pub max_len_bp: Option<usize>,
    pub limit: Option<usize>,
    pub offset: usize,
    pub sort_by: SequenceFeatureSortBy,
    pub descending: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQueryRow {
    pub feature_id: usize,
    pub kind: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub length_bp: usize,
    pub strand: String,
    pub label: String,
    pub labels: Vec<String>,
    pub qualifiers: BTreeMap<String, Vec<String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQueryResult {
    pub schema: String,
    pub seq_id: SeqId,
    pub sequence_length_bp: usize,
    pub total_feature_count: usize,
    pub matched_count: usize,
    pub returned_count: usize,
    pub offset: usize,
    pub limit: usize,
    pub range_relation: String,
    pub strand_filter: String,
    pub sort_by: String,
    pub descending: bool,
    pub query: SequenceFeatureQuery,
    pub rows: Vec<SequenceFeatureQueryRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Progress payload emitted by TFBS annotation operations.
pub struct TfbsProgress {
    pub seq_id: String,
    pub motif_id: String,
    pub motif_index: usize,
    pub motif_count: usize,
    pub scanned_steps: usize,
    pub total_steps: usize,
    pub motif_percent: f64,
    pub total_percent: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Progress payload emitted by genome track import operations.
pub struct GenomeTrackImportProgress {
    pub seq_id: String,
    pub source: String,
    pub path: String,
    pub parsed_records: usize,
    pub imported_features: usize,
    pub skipped_records: usize,
    pub done: bool,
}

/// Engine capability snapshot used by adapters for discovery/negotiation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Capabilities {
    pub protocol_version: String,
    pub supported_operations: Vec<String>,
    pub supported_export_formats: Vec<String>,
    pub deterministic_operation_log: bool,
}
