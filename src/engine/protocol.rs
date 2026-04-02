//! Stable public analysis/report contracts extracted from the monolithic
//! engine.
//!
//! This is the narrow, serialization-friendly layer that GUI/CLI/JS/Lua/Python
//! adapters should lean on when they only need record shapes rather than full
//! engine execution.
//!
//! Look here for:
//! - persisted report/result payloads such as RNA-read, dotplot, planning, and
//!   sequencing-confirmation records
//! - small enums that appear in JSON-facing operation contracts
//! - state-summary structs that should remain slower-changing than
//!   `src/engine.rs`

pub use gentle_protocol::{
    Capabilities, DotplotBoxplotBin, DotplotMatchPoint, DotplotMode, DotplotOverlayQuerySpec,
    DotplotQuerySeries, DotplotReferenceAnnotationInterval, DotplotReferenceAnnotationTrack,
    DotplotView, DotplotViewSummary, EngineError, ErrorCode, FlexibilityModel,
    GenomeTrackImportProgress, PairwiseAlignmentMode, RnaReadAlignConfig, RnaReadAlignmentBackend,
    RnaReadAlignmentDisplay, RnaReadAlignmentDotplotSvgExport, RnaReadAlignmentEffect,
    RnaReadAlignmentInspection, RnaReadAlignmentInspectionEffectFilter,
    RnaReadAlignmentInspectionRow, RnaReadAlignmentInspectionSortKey,
    RnaReadAlignmentInspectionSubsetSpec, RnaReadAlignmentMode, RnaReadAlignmentTsvExport,
    RnaReadExonAbundanceExport, RnaReadExonPathsExport, RnaReadExonSupportFrequency,
    RnaReadHitSelection, RnaReadInputFormat, RnaReadInterpretProgress, RnaReadInterpretationHit,
    RnaReadInterpretationProfile, RnaReadInterpretationReport, RnaReadInterpretationReportSummary,
    RnaReadIsoformSupportRow, RnaReadJunctionSupportFrequency, RnaReadMappedIsoformSupportRow,
    RnaReadMappedSupportExonAttribution, RnaReadMappedSupportJunctionAttribution,
    RnaReadMappingHit, RnaReadOriginCandidateContribution, RnaReadOriginClass, RnaReadOriginMode,
    RnaReadPairwiseAlignmentDetail, RnaReadReportMode, RnaReadSampleSheetExport,
    RnaReadScoreDensityScale, RnaReadScoreDensitySvgExport, RnaReadScoreDensityVariant,
    RnaReadSeedFilterConfig, RnaReadSeedHistogramBin, RnaReadStrandAssignmentDiagnostics,
    RnaReadTopHitPreview, RnaReadTransitionSupportRow, RnaSeedHashCatalogEntry,
    RnaSeedHashTemplateAuditEntry, SequenceAlignmentReport, SequenceFeatureQualifierFilter,
    SequenceFeatureQuery, SequenceFeatureQueryResult, SequenceFeatureQueryRow,
    SequenceFeatureRangeRelation, SequenceFeatureSortBy, SequenceFeatureStrandFilter,
    SequencingPrimerOrientation, SequencingPrimerOverlayReport, SequencingPrimerOverlaySuggestion,
    SequencingPrimerProblemGuidanceRow, SequencingPrimerProblemKind, SequencingPrimerProposalRow,
    SequencingTraceChannelData, SequencingTraceChannelSummary, SequencingTraceFormat,
    SequencingTraceImportReport, SequencingTraceRecord, SequencingTraceSummary, TfbsProgress,
};
use serde::{Deserialize, Serialize};

use crate::enzymes::default_preferred_restriction_enzyme_names;

use super::{
    OpId, Operation, PrepareGenomeProgress, ProtocolCartoonTemplateBindings, RunId, SeqId,
    SequencingConfirmationReport,
};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadGeneSupportCompleteRule {
    #[default]
    Near,
    Strict,
    Exact,
}

impl RnaReadGeneSupportCompleteRule {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Near => "near",
            Self::Strict => "strict",
            Self::Exact => "exact",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneExonSupportRow {
    pub gene_id: String,
    pub exon_ordinal: usize,
    pub start_1based: usize,
    pub end_1based: usize,
    pub support_read_count: usize,
    pub support_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneExonPairSupportRow {
    pub gene_id: String,
    pub from_exon_ordinal: usize,
    pub from_start_1based: usize,
    pub from_end_1based: usize,
    pub to_exon_ordinal: usize,
    pub to_start_1based: usize,
    pub to_end_1based: usize,
    pub support_read_count: usize,
    pub support_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportCohortSummary {
    pub read_count: usize,
    pub exon_support: Vec<RnaReadGeneExonSupportRow>,
    pub exon_pair_support: Vec<RnaReadGeneExonPairSupportRow>,
    pub direct_transition_support: Vec<RnaReadGeneExonPairSupportRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportSummary {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub requested_gene_ids: Vec<String>,
    pub matched_gene_ids: Vec<String>,
    pub missing_gene_ids: Vec<String>,
    pub selected_record_indices: Vec<usize>,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
    pub aligned_base_count: usize,
    pub accepted_target_count: usize,
    pub fragment_count: usize,
    pub complete_count: usize,
    pub complete_strict_count: usize,
    pub complete_exact_count: usize,
    pub all_target: RnaReadGeneSupportCohortSummary,
    pub fragments: RnaReadGeneSupportCohortSummary,
    pub complete: RnaReadGeneSupportCohortSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Named visibility targets controlled through `Operation::SetDisplayVisibility`.
///
/// Adapters should treat these as the canonical shared display toggles.
pub enum DisplayTarget {
    SequencePanel,
    MapPanel,
    Features,
    CdsFeatures,
    GeneFeatures,
    MrnaFeatures,
    Tfbs,
    RestrictionEnzymes,
    GcContents,
    OpenReadingFrames,
    MethylationSites,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Hash, Default)]
#[serde(rename_all = "snake_case")]
/// Rendering mode selector for linear DNA letter layout.
///
/// `AutoAdaptive` is the default policy and should be preferred for parity
/// unless the caller explicitly requests one fixed layout mode.
pub enum LinearSequenceLetterLayoutMode {
    #[default]
    AutoAdaptive,
    StandardLinear,
    #[serde(alias = "continuous")]
    ContinuousHelical,
    #[serde(
        alias = "condensed_10_row",
        alias = "condensed10row",
        alias = "condensed"
    )]
    Condensed10Row,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Hash, Default)]
#[serde(rename_all = "snake_case")]
pub enum RestrictionEnzymeDisplayMode {
    PreferredOnly,
    #[default]
    PreferredAndUnique,
    UniqueOnly,
    AllInView,
}

impl RestrictionEnzymeDisplayMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::PreferredOnly => "Preferred only",
            Self::PreferredAndUnique => "Preferred + unique",
            Self::UniqueOnly => "Unique only",
            Self::AllInView => "All in view",
        }
    }

    pub fn short_label(self) -> &'static str {
        match self {
            Self::PreferredOnly => "Pref",
            Self::PreferredAndUnique => "Pref+uniq",
            Self::UniqueOnly => "Unique",
            Self::AllInView => "All",
        }
    }

    pub fn count_label(self) -> &'static str {
        match self {
            Self::PreferredOnly => "preferred cutters",
            Self::PreferredAndUnique => "preferred/unique cutters",
            Self::UniqueOnly => "unique cutters",
            Self::AllInView => "cut sites",
        }
    }

    pub fn empty_state_label(self) -> &'static str {
        match self {
            Self::PreferredOnly => "No preferred cutters in view.",
            Self::PreferredAndUnique => "No preferred or unique cutters in view.",
            Self::UniqueOnly => "No unique cutters in view.",
            Self::AllInView => "No restriction cut sites in view.",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
/// Project-level display settings persisted in `ProjectState`.
///
/// These settings are consumed by GUI renderers and export paths, and are kept
/// in engine state so behavior is adapter-equivalent.
pub struct DisplaySettings {
    pub show_sequence_panel: bool,
    #[serde(default = "DisplaySettings::default_sequence_panel_max_text_length_bp")]
    pub sequence_panel_max_text_length_bp: usize,
    pub auto_hide_sequence_panel_when_linear_bases_visible: bool,
    pub show_map_panel: bool,
    pub show_features: bool,
    pub show_cds_features: bool,
    pub show_gene_features: bool,
    pub show_mrna_features: bool,
    pub show_tfbs: bool,
    pub regulatory_tracks_near_baseline: bool,
    pub regulatory_feature_max_view_span_bp: usize,
    pub tfbs_display_use_llr_bits: bool,
    pub tfbs_display_min_llr_bits: f64,
    pub tfbs_display_use_llr_quantile: bool,
    pub tfbs_display_min_llr_quantile: f64,
    pub tfbs_display_use_true_log_odds_bits: bool,
    pub tfbs_display_min_true_log_odds_bits: f64,
    pub tfbs_display_use_true_log_odds_quantile: bool,
    pub tfbs_display_min_true_log_odds_quantile: f64,
    pub vcf_display_show_snp: bool,
    pub vcf_display_show_ins: bool,
    pub vcf_display_show_del: bool,
    pub vcf_display_show_sv: bool,
    pub vcf_display_show_other: bool,
    pub vcf_display_pass_only: bool,
    pub vcf_display_use_min_qual: bool,
    pub vcf_display_min_qual: f64,
    pub vcf_display_use_max_qual: bool,
    pub vcf_display_max_qual: f64,
    #[serde(default)]
    pub vcf_display_required_info_keys: Vec<String>,
    pub show_restriction_enzymes: bool,
    #[serde(default)]
    pub restriction_enzyme_display_mode: RestrictionEnzymeDisplayMode,
    #[serde(default = "DisplaySettings::default_preferred_restriction_enzymes")]
    pub preferred_restriction_enzymes: Vec<String>,
    pub show_gc_contents: bool,
    #[serde(default = "DisplaySettings::default_gc_content_bin_size_bp")]
    pub gc_content_bin_size_bp: usize,
    pub show_open_reading_frames: bool,
    pub show_methylation_sites: bool,
    pub linear_view_start_bp: usize,
    pub linear_view_span_bp: usize,
    pub linear_view_vertical_offset_px: f32,
    pub linear_sequence_base_text_max_view_span_bp: usize,
    pub linear_sequence_helical_letters_enabled: bool,
    pub linear_sequence_helical_max_view_span_bp: usize,
    #[serde(default = "DisplaySettings::default_linear_sequence_condensed_max_view_span_bp")]
    pub linear_sequence_condensed_max_view_span_bp: usize,
    #[serde(default)]
    pub linear_sequence_letter_layout_mode: LinearSequenceLetterLayoutMode,
    pub linear_sequence_helical_phase_offset_bp: usize,
    pub linear_show_double_strand_bases: bool,
    #[serde(default = "DisplaySettings::default_linear_helical_parallel_strands")]
    pub linear_helical_parallel_strands: bool,
    pub linear_hide_backbone_when_sequence_bases_visible: bool,
    pub linear_reverse_strand_use_upside_down_letters: bool,
    #[serde(default = "DisplaySettings::default_reverse_strand_visual_opacity")]
    pub reverse_strand_visual_opacity: f32,
    pub feature_details_font_size: f32,
    pub linear_external_feature_label_font_size: f32,
    pub linear_external_feature_label_background_opacity: f32,
}

impl DisplaySettings {
    pub const fn default_sequence_panel_max_text_length_bp() -> usize {
        200_000
    }

    pub const fn default_gc_content_bin_size_bp() -> usize {
        100
    }

    pub const fn default_linear_sequence_condensed_max_view_span_bp() -> usize {
        1500
    }

    pub const fn default_linear_helical_parallel_strands() -> bool {
        true
    }

    pub const fn default_reverse_strand_visual_opacity() -> f32 {
        0.55
    }

    pub fn default_preferred_restriction_enzymes() -> Vec<String> {
        default_preferred_restriction_enzyme_names()
    }
}

impl Default for DisplaySettings {
    fn default() -> Self {
        Self {
            show_sequence_panel: true,
            sequence_panel_max_text_length_bp: Self::default_sequence_panel_max_text_length_bp(),
            auto_hide_sequence_panel_when_linear_bases_visible: false,
            show_map_panel: true,
            show_features: true,
            show_cds_features: true,
            show_gene_features: true,
            show_mrna_features: true,
            show_tfbs: false,
            regulatory_tracks_near_baseline: false,
            regulatory_feature_max_view_span_bp: 50_000,
            tfbs_display_use_llr_bits: true,
            tfbs_display_min_llr_bits: 0.0,
            tfbs_display_use_llr_quantile: true,
            tfbs_display_min_llr_quantile: 0.95,
            tfbs_display_use_true_log_odds_bits: false,
            tfbs_display_min_true_log_odds_bits: 0.0,
            tfbs_display_use_true_log_odds_quantile: false,
            tfbs_display_min_true_log_odds_quantile: 0.95,
            vcf_display_show_snp: true,
            vcf_display_show_ins: true,
            vcf_display_show_del: true,
            vcf_display_show_sv: true,
            vcf_display_show_other: true,
            vcf_display_pass_only: false,
            vcf_display_use_min_qual: false,
            vcf_display_min_qual: 0.0,
            vcf_display_use_max_qual: false,
            vcf_display_max_qual: 0.0,
            vcf_display_required_info_keys: vec![],
            show_restriction_enzymes: true,
            restriction_enzyme_display_mode: RestrictionEnzymeDisplayMode::default(),
            preferred_restriction_enzymes: Self::default_preferred_restriction_enzymes(),
            show_gc_contents: true,
            gc_content_bin_size_bp: Self::default_gc_content_bin_size_bp(),
            show_open_reading_frames: false,
            show_methylation_sites: false,
            linear_view_start_bp: 0,
            linear_view_span_bp: 0,
            linear_view_vertical_offset_px: 0.0,
            linear_sequence_base_text_max_view_span_bp: 500,
            linear_sequence_helical_letters_enabled: true,
            linear_sequence_helical_max_view_span_bp: 2000,
            linear_sequence_condensed_max_view_span_bp:
                Self::default_linear_sequence_condensed_max_view_span_bp(),
            linear_sequence_letter_layout_mode: LinearSequenceLetterLayoutMode::AutoAdaptive,
            linear_sequence_helical_phase_offset_bp: 0,
            linear_show_double_strand_bases: true,
            linear_helical_parallel_strands: Self::default_linear_helical_parallel_strands(),
            linear_hide_backbone_when_sequence_bases_visible: false,
            linear_reverse_strand_use_upside_down_letters: true,
            reverse_strand_visual_opacity: Self::default_reverse_strand_visual_opacity(),
            feature_details_font_size: 9.0,
            linear_external_feature_label_font_size: 11.0,
            linear_external_feature_label_background_opacity: 0.9,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// One deterministic workflow run: ordered operations with a caller-supplied
/// `run_id`.
///
/// Operations are applied sequentially. The current workflow executor is not
/// transactional: if a later step fails, earlier successful steps remain in
/// state and in the operation journal.
pub struct Workflow {
    pub run_id: RunId,
    pub ops: Vec<Operation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Canonical result payload returned after one operation completes.
///
/// `created_seq_ids` and `changed_seq_ids` are the stable adapter-facing hint
/// for which sequence windows/views may need refresh after an operation.
pub struct OpResult {
    pub op_id: OpId,
    pub created_seq_ids: Vec<SeqId>,
    pub changed_seq_ids: Vec<SeqId>,
    pub warnings: Vec<String>,
    pub messages: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub protocol_cartoon_preview: Option<ProtocolCartoonPreviewTelemetry>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_annotation_projection: Option<GenomeAnnotationProjectionTelemetry>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequence_alignment: Option<SequenceAlignmentReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequencing_confirmation_report: Option<SequencingConfirmationReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequencing_trace_import_report: Option<SequencingTraceImportReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequencing_trace_record: Option<SequencingTraceRecord>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequencing_trace_summaries: Option<Vec<SequencingTraceSummary>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequencing_primer_overlay_report: Option<SequencingPrimerOverlayReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rna_read_gene_support_summary: Option<RnaReadGeneSupportSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Optional protocol-cartoon preview payload emitted by operations that can
/// project a deterministic mechanism strip from operation geometry.
pub struct ProtocolCartoonPreviewTelemetry {
    pub protocol: String,
    pub flank_bp: usize,
    pub overlap_bp: usize,
    pub insert_bp: usize,
    pub bindings: ProtocolCartoonTemplateBindings,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Structured annotation projection telemetry emitted by genomic extraction
/// operations (`ExtractGenomeRegion`, `ExtractGenomeGene`).
pub struct GenomeAnnotationProjectionTelemetry {
    pub requested_scope: String,
    pub effective_scope: String,
    pub max_features_cap: Option<usize>,
    pub candidate_feature_count: usize,
    pub attached_feature_count: usize,
    pub dropped_feature_count: usize,
    pub genes_attached: usize,
    pub transcripts_attached: usize,
    pub exons_attached: usize,
    pub cds_attached: usize,
    pub fallback_applied: bool,
    pub fallback_reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DbSnpFetchStage {
    ValidateInput,
    InspectPreparedGenome,
    ContactServer,
    WaitResponse,
    ParseResponse,
    ResolvePlacement,
    ExtractRegion,
    AttachVariantMarker,
}

impl DbSnpFetchStage {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::ValidateInput => "validate_input",
            Self::InspectPreparedGenome => "inspect_prepared_genome",
            Self::ContactServer => "contact_server",
            Self::WaitResponse => "wait_response",
            Self::ParseResponse => "parse_response",
            Self::ResolvePlacement => "resolve_placement",
            Self::ExtractRegion => "extract_region",
            Self::AttachVariantMarker => "attach_variant_marker",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Progress snapshot for shared dbSNP locus resolution and extraction.
pub struct DbSnpFetchProgress {
    pub rs_id: String,
    pub genome_id: String,
    pub stage: DbSnpFetchStage,
    pub detail: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Union of long-running operation progress events.
pub enum OperationProgress {
    Tfbs(TfbsProgress),
    GenomePrepare(PrepareGenomeProgress),
    GenomeTrackImport(GenomeTrackImportProgress),
    DbSnpFetch(DbSnpFetchProgress),
    RnaReadInterpret(RnaReadInterpretProgress),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Immutable operation journal row.
pub struct OperationRecord {
    pub run_id: RunId,
    pub op: Operation,
    pub result: OpResult,
}

/// Compact sequence row used by state-summary style adapter surfaces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineSequenceSummary {
    pub id: String,
    pub name: Option<String>,
    pub length: usize,
    pub circular: bool,
}

/// Compact container row used by shell/CLI inspection surfaces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineContainerSummary {
    pub id: String,
    pub kind: String,
    pub member_count: usize,
    pub members: Vec<String>,
}

/// Compact arrangement row used by shell/CLI inspection surfaces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineArrangementSummary {
    pub id: String,
    pub mode: String,
    pub lane_count: usize,
    pub lane_container_ids: Vec<String>,
    pub ladders: Vec<String>,
}

/// Machine-readable snapshot of top-level engine state counts and summaries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineStateSummary {
    pub sequence_count: usize,
    pub sequences: Vec<EngineSequenceSummary>,
    pub container_count: usize,
    pub containers: Vec<EngineContainerSummary>,
    pub arrangement_count: usize,
    pub arrangements: Vec<EngineArrangementSummary>,
    pub display: DisplaySettings,
}
