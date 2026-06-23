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

use std::collections::HashMap;

pub use gentle_protocol::{
    ANNOTATION_CANDIDATE_SCHEMA, ANNOTATION_CANDIDATE_SUMMARY_SCHEMA,
    ANNOTATION_CANDIDATE_WRITEBACK_SCHEMA, AdapterCaptureProtectionMode, AdapterCaptureStyle,
    AdapterRestrictionCapturePlan, AnchorBoundary, AnchorDirection, AnnotationCandidate,
    AnnotationCandidateSummary, AnnotationCandidateWriteback, AttractPwmMappingPolicy,
    AttractRegionClass, AttractSpeciesMatchMode, AttractSplicingEvidenceHitRow,
    AttractSplicingEvidencePolicySummary, AttractSplicingEvidenceSettings,
    AttractSplicingEvidenceSummaryRow, AttractSplicingEvidenceView, CONSTRUCT_CANDIDATE_SCHEMA,
    CONSTRUCT_OBJECTIVE_SCHEMA, CONSTRUCT_REASONING_GRAPH_SCHEMA, CONSTRUCT_REASONING_STORE_SCHEMA,
    CandidateFeatureBoundaryMode, CandidateFeatureGeometryMode, CandidateFeatureStrandRelation,
    CandidateMacroTemplateParam, CandidateObjectiveDirection, CandidateObjectiveSpec,
    CandidateSetOperator, CandidateTieBreakPolicy, CandidateWeightedObjectiveTerm, Capabilities,
    CdnaAssayMaterializedProductRow, CdnaAssayProductGelBandRow, CdnaAssayProductMaterialization,
    CdnaAssayTranscriptMapCoordinateMode, CdnaAssayTranscriptOrder, ConstructCandidate,
    ConstructObjective, ConstructReasoningGraph, ConstructReasoningStore, ConstructRole,
    ContainerId, ContainerKind, CutRunAlignConfig, CutRunCatalogEntry, CutRunCatalogListEntry,
    CutRunCoverageKind, CutRunDatasetListReport, CutRunDatasetProjectionReport,
    CutRunDatasetStatus, CutRunFragmentSpan, CutRunInputFormat,
    CutRunMotifAbsentOccupancyInterpretation, CutRunMotifAbsentSupportWindow,
    CutRunMotifContextHit, CutRunMotifContextScope, CutRunMotifContextSummaryRow,
    CutRunPreparedAssetManifest, CutRunPreparedAssetStatus, CutRunPreparedManifest,
    CutRunReadCoverageExport, CutRunReadLayout, CutRunReadOrientation, CutRunReadPlacement,
    CutRunReadReport, CutRunReadReportStore, CutRunReadReportSummary, CutRunReadUnitRow,
    CutRunReadUnitStatus, CutRunRegulatoryEvidenceSourceKind, CutRunRegulatoryEvidenceSourceRef,
    CutRunRegulatorySupportReport, CutRunRegulatoryTfbsConfirmationStatus, CutRunRegulatoryTfbsRow,
    CutRunSeedFilterConfig, CutRunSupportCluster, CutRunSupportStrength, CutRunSupportWindowRecord,
    DESIGN_DECISION_NODE_SCHEMA, DESIGN_EVIDENCE_SCHEMA, DESIGN_FACT_SCHEMA, DecisionMethod,
    DesignDecisionNode, DesignEvidence, DesignFact, DotplotBoxplotBin, DotplotMatchPoint,
    DotplotMode, DotplotOverlayAnchorExon, DotplotOverlayAnchorExonRef,
    DotplotOverlayAnchorSeriesSupport, DotplotOverlayQuerySpec, DotplotOverlayResolvedAnchorSeries,
    DotplotOverlayXAxisMode, DotplotQuerySeries, DotplotReferenceAnnotationInterval,
    DotplotReferenceAnnotationTrack, DotplotView, DotplotViewSummary,
    EXON_SKIP_MATERIALIZATION_SCHEMA, EXON_SKIP_SELECTION_PLAN_SCHEMA, EditableStatus, EngineError,
    ErrorCode, EvidenceClass, EvidenceScope, ExonSkipCandidateExon, ExonSkipMaterializationReport,
    ExonSkipReturnKind, ExonSkipReturnPayload, ExonSkipSelectionCriterion, ExonSkipSelectionPlan,
    FeatureBedCoordinateMode, FlexibilityModel, GENE_SET_CO_REGULATED_CACHE_SCHEMA,
    GENE_SET_CUTRUN_REGULATORY_SUPPORT_SCHEMA, GENE_SET_DIRECT_LIST_CACHE_SCHEMA,
    GENE_SET_ONTOLOGY_ASSIGNMENT_CACHE_SCHEMA, GENE_SET_PROMOTER_COHORT_SCHEMA,
    GENE_SET_RESOLUTION_SCHEMA, GeneSetCoRegulatedProducerMetadata, GeneSetCohortRelationship,
    GeneSetCohortRelationshipFlag, GeneSetCutRunEvaluationState, GeneSetCutRunMemberSupport,
    GeneSetCutRunRegulatorySupportReport, GeneSetCutRunSupportAggregate, GeneSetProducerFilter,
    GeneSetProducerKind, GeneSetProducerProvenance, GeneSetProducerQueryMetadata,
    GeneSetPromoterCohortReport, GeneSetPromoterWindow, GeneSetProvenanceRow,
    GeneSetRandomProvenance, GeneSetRequest, GeneSetResolutionReport,
    GeneSetResolutionReviewStatus, GeneSetResolvedMember, GeneSetUnresolvedMember,
    GenomeAnchorSide, GenomeAnnotationScope, GenomeGeneExtractMode, GenomeTrackImportProgress,
    GenomeTrackSource, GenomeTrackSubscription, HOST_PROFILE_CATALOG_SCHEMA,
    HelperConstructProfile, HostLifecycleRole, HostProfileCatalog, HostProfileRecord,
    HostRouteStep, ORTHOLOG_PROMOTER_COHORT_SCHEMA, ORTHOLOG_PROMOTER_COMPARISON_SCHEMA,
    ORTHOLOG_RESOURCE_SCHEMA, OrthologAmbiguityPolicy, OrthologCutRunSupportRow,
    OrthologCutRunSupportStatus, OrthologExpressionAssignment, OrthologMappingRow,
    OrthologPairwiseTfbsSimilarity, OrthologPromoterCohortReport, OrthologPromoterCohortRequest,
    OrthologPromoterComparisonReport, OrthologPromoterRole, OrthologPromoterRow, OrthologResource,
    OrthologSequenceSimilarityRow, OrthologSpeciesAlias, OrthologTfbsPeakSummary,
    OrthologTfbsSummaryRow, OrthologUnresolvedRow, PairwiseAlignmentMode, PortBindingStatus,
    PreparedCacheCleanupMode, PreparedCacheCleanupRequest, PrimerDesignBackend,
    PrimerSpecificityCheckMode, PrimerSpecificityPolicy, ProteinResidueGenomicCoordinateBase,
    ProteinResidueGenomicCoordinateMatch, ProteinResidueGenomicCoordinateReport,
    ProteinToDnaHandoffCandidate, ProteinToDnaHandoffCoverage, ProteinToDnaHandoffRankingGoal,
    ProteinToDnaHandoffStrategy, ProtocolCartoonKind, QpcrTranscriptSpecificityEvidence,
    QpcrTranscriptTargetingMode, READ_ACQUISITION_REPORT_SCHEMA, REPORTER_CATALOG_REPORT_SCHEMA,
    REPORTER_CATALOG_SCHEMA, REPORTER_CONSTRUCT_HANDOFF_SCHEMA, REPORTER_CORPUS_EXPORT_SCHEMA,
    REPORTER_RECOMMENDATION_SCHEMA, RNA_READ_ALIGNMENT_DISPLAY_BATCH_SCHEMA,
    RNA_READ_BATCH_MAP_REPORT_SCHEMA, RNA_READ_GENE_SCREEN_SUMMARY_SCHEMA,
    RNA_READ_TRANSCRIPT_CATALOG_INDEX_SCHEMA, ReadAcquisitionAnalysisFormat,
    ReadAcquisitionCommandProvenance, ReadAcquisitionManifestRow, ReadAcquisitionOutputPath,
    ReadAcquisitionReadLayout, ReadAcquisitionReadLengthStats, ReadAcquisitionReport,
    ReadAcquisitionRunReport, RenderSvgMode, ReporterAnnotatedRecord, ReporterBackboneResolution,
    ReporterBackboneResolutionStatus, ReporterCatalog, ReporterCatalogReport,
    ReporterComputedAnnotation, ReporterConstraints, ReporterConstructHandoffCommand,
    ReporterConstructHandoffPlan, ReporterConstructHandoffProvenance, ReporterConstructPortBinding,
    ReporterConstructSelectedFragment, ReporterConstructSelectedReporter, ReporterCorpusExport,
    ReporterCorpusExportFormat, ReporterPreferenceWeights, ReporterQuarantinedRecord,
    ReporterRecommendation, ReporterRecommendationResult, ReporterRecord,
    ReporterRejectedCandidate, ReporterSourceRef, ReporterSpectralProfile,
    RestrictionCloningPcrHandoffMode, RnaReadAlignConfig, RnaReadAlignmentBackend,
    RnaReadAlignmentDisplay, RnaReadAlignmentDotplotSvgExport, RnaReadAlignmentEffect,
    RnaReadAlignmentInspection, RnaReadAlignmentInspectionEffectFilter,
    RnaReadAlignmentInspectionRow, RnaReadAlignmentInspectionSortKey,
    RnaReadAlignmentInspectionSubsetSpec, RnaReadAlignmentMode, RnaReadAlignmentTsvExport,
    RnaReadBatchConcatemerPartnerRow, RnaReadBatchIsoformSupportRow, RnaReadBatchMapReport,
    RnaReadBatchMapSampleRow, RnaReadBatchMapSampleStatus, RnaReadBatchMapSraPreparationRow,
    RnaReadConcatemerAdapterHit, RnaReadConcatemerFragmentOrigin, RnaReadConcatemerInspection,
    RnaReadConcatemerInspectionSettings, RnaReadConcatemerPartnerGeneSummary,
    RnaReadConcatemerPartnerTranscriptSummary, RnaReadConcatemerSuspicionLevel,
    RnaReadConcatemerSuspicionRow, RnaReadExonAbundanceExport, RnaReadExonPathsExport,
    RnaReadExonSupportFrequency, RnaReadHitSelection, RnaReadInputFormat, RnaReadInterpretProgress,
    RnaReadInterpretationHit, RnaReadInterpretationProfile, RnaReadInterpretationReport,
    RnaReadInterpretationReportSummary, RnaReadIsoformPreflightControlSummary,
    RnaReadIsoformPreflightReport, RnaReadIsoformPreflightScore,
    RnaReadIsoformPreflightThresholdRecommendation, RnaReadIsoformSupportRow,
    RnaReadIsoformTriageBin, RnaReadIsoformTriageTsvExport, RnaReadJunctionSupportFrequency,
    RnaReadMappedIsoformSupportRow, RnaReadMappedSupportExonAttribution,
    RnaReadMappedSupportJunctionAttribution, RnaReadMappingHit, RnaReadOriginCandidateContribution,
    RnaReadOriginClass, RnaReadOriginMode, RnaReadPairwiseAlignmentDetail, RnaReadReportMode,
    RnaReadSampleSheetExport, RnaReadScoreDensityScale, RnaReadScoreDensitySvgExport,
    RnaReadScoreDensityVariant, RnaReadSeedFilterConfig, RnaReadSeedHistogramBin,
    RnaReadStrandAssignmentDiagnostics, RnaReadTopHitPreview, RnaReadTranscriptCatalogIndex,
    RnaReadTranscriptCatalogTemplateRecord, RnaReadTransitionSupportRow, RnaSeedHashCatalogEntry,
    RnaSeedHashTemplateAuditEntry, SequenceAlignmentReport, SequenceAnchor,
    SequenceFeatureBedExportReport, SequenceFeatureQualifierFilter, SequenceFeatureQuery,
    SequenceFeatureQueryResult, SequenceFeatureQueryRow, SequenceFeatureRangeRelation,
    SequenceFeatureSortBy, SequenceFeatureStrandFilter, SequencingPrimerOrientation,
    SequencingPrimerOverlayReport, SequencingPrimerOverlaySuggestion,
    SequencingPrimerProblemGuidanceRow, SequencingPrimerProblemKind, SequencingPrimerProposalRow,
    SequencingTraceChannelData, SequencingTraceChannelSummary, SequencingTraceFormat,
    SequencingTraceImportReport, SequencingTraceRecord, SequencingTraceSummary,
    SharedAssetActivityStatus, SplicingScopePreset, TfThresholdOverride, TfbsProgress,
    TranscriptProteinDerivation, TranscriptProteinDerivationMode,
    TranscriptProteinTranslationTableSource, TranslationSpeedMark, TranslationSpeedProfile,
    TranslationSpeedProfileSource, UniprotFeatureCodingDnaExonPair,
    UniprotFeatureCodingDnaExonSpan, UniprotFeatureCodingDnaMatch,
    UniprotFeatureCodingDnaQueryMode, UniprotFeatureCodingDnaQueryReport,
    UniprotFeatureCodingDnaSegment,
};
use serde::{Deserialize, Deserializer, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;

use crate::enzymes::default_preferred_restriction_enzyme_names;
use crate::genomes::BlastExternalBinaryPreflightReport;

use super::{
    CLONING_MACRO_TEMPLATE_SCHEMA, OpId, Operation, PrepareGenomeProgress,
    ProtocolCartoonTemplateBindings, RunId, SeqId,
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadGeneSupportAuditStatus {
    #[default]
    Unaligned,
    AlignedOtherGene,
    AcceptedFragment,
    AcceptedComplete,
}

impl RnaReadGeneSupportAuditStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Unaligned => "unaligned",
            Self::AlignedOtherGene => "aligned_other_gene",
            Self::AcceptedFragment => "accepted_fragment",
            Self::AcceptedComplete => "accepted_complete",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadGeneSupportAuditCohortFilter {
    #[default]
    All,
    Accepted,
    Fragment,
    Complete,
    Rejected,
}

impl RnaReadGeneSupportAuditCohortFilter {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::Accepted => "accepted",
            Self::Fragment => "fragment",
            Self::Complete => "complete",
            Self::Rejected => "rejected",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentDisplayBatchEntry {
    pub record_index: usize,
    pub header_id: String,
    pub alignment: RnaReadAlignmentDisplay,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentDisplayBatchSkippedRecord {
    pub record_index: usize,
    pub header_id: String,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentDisplayBatch {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub seed_feature_id: usize,
    pub requested_gene_ids: Vec<String>,
    pub matched_gene_ids: Vec<String>,
    pub missing_gene_ids: Vec<String>,
    pub selection_mode: String,
    #[serde(default)]
    pub cohort_filter: RnaReadGeneSupportAuditCohortFilter,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
    pub selected_record_indices: Vec<usize>,
    pub limit: Option<usize>,
    pub entry_count: usize,
    pub skipped_record_indices: Vec<usize>,
    pub skipped_records: Vec<RnaReadAlignmentDisplayBatchSkippedRecord>,
    pub entries: Vec<RnaReadAlignmentDisplayBatchEntry>,
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
pub struct RnaReadLengthDistributionSummary {
    pub sample_count: usize,
    pub mean_length_bp: f64,
    pub min_length_bp: usize,
    pub q25_length_bp: usize,
    pub median_length_bp: usize,
    pub q75_length_bp: usize,
    pub max_length_bp: usize,
    pub p95_length_bp: usize,
    pub length_counts: Vec<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadFractionDistributionSummary {
    pub sample_count: usize,
    pub mean_fraction: f64,
    pub median_fraction: f64,
    pub p95_fraction: f64,
    pub bin_counts: Vec<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportSummary {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub source_report_generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_report_op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_report_run_id: Option<String>,
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
    pub evaluated_read_lengths: RnaReadLengthDistributionSummary,
    pub accepted_target_read_lengths: RnaReadLengthDistributionSummary,
    pub accepted_target_fragment_lengths: RnaReadLengthDistributionSummary,
    pub accepted_target_query_coverage: RnaReadFractionDistributionSummary,
    pub all_target: RnaReadGeneSupportCohortSummary,
    pub fragments: RnaReadGeneSupportCohortSummary,
    pub complete: RnaReadGeneSupportCohortSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadTargetQualityComparisonEntry {
    pub entry_id: String,
    pub comparison_label: String,
    pub gentle_version: String,
    pub report_id: String,
    pub seq_id: String,
    pub report_generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub report_op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub report_run_id: Option<String>,
    #[serde(default)]
    pub requested_gene_ids: Vec<String>,
    #[serde(default)]
    pub matched_gene_ids: Vec<String>,
    #[serde(default)]
    pub missing_gene_ids: Vec<String>,
    #[serde(default)]
    pub report_target_gene_ids: Vec<String>,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
    #[serde(default)]
    pub profile: RnaReadInterpretationProfile,
    #[serde(default)]
    pub scope: SplicingScopePreset,
    #[serde(default)]
    pub origin_mode: RnaReadOriginMode,
    #[serde(default)]
    pub report_mode: RnaReadReportMode,
    #[serde(default)]
    pub seed_filter: RnaReadSeedFilterConfig,
    #[serde(default)]
    pub align_config: RnaReadAlignConfig,
    #[serde(default)]
    pub all_read_lengths: RnaReadLengthDistributionSummary,
    pub summary: RnaReadGeneSupportSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadTargetQualityComparisonBundle {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default)]
    pub entries: Vec<RnaReadTargetQualityComparisonEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadTargetQualityExport {
    pub schema: String,
    pub requested_path: String,
    pub written_path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bundle_path: Option<String>,
    pub format: String,
    pub report_id: String,
    #[serde(default)]
    pub requested_gene_ids: Vec<String>,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
    pub entry_count: usize,
    #[serde(default)]
    pub appended_to_existing_bundle: bool,
    #[serde(default)]
    pub reused_existing_entry_slot: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportAuditPair {
    pub from_exon_ordinal: usize,
    pub from_start_1based: usize,
    pub from_end_1based: usize,
    pub to_exon_ordinal: usize,
    pub to_start_1based: usize,
    pub to_end_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportAuditRow {
    pub record_index: usize,
    pub header_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_label: Option<String>,
    #[serde(default)]
    pub status: RnaReadGeneSupportAuditStatus,
    pub status_reason: String,
    #[serde(default)]
    pub full_length_exact: bool,
    #[serde(default)]
    pub full_length_near: bool,
    #[serde(default)]
    pub full_length_strict: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub full_length_class: Option<String>,
    #[serde(default)]
    pub mapped_exon_ordinals: Vec<usize>,
    #[serde(default)]
    pub exon_pairs: Vec<RnaReadGeneSupportAuditPair>,
    #[serde(default)]
    pub direct_transition_pairs: Vec<RnaReadGeneSupportAuditPair>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub score: Option<isize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub identity_fraction: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_coverage_fraction: Option<f64>,
    #[serde(default)]
    pub passed_seed_filter: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportAudit {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub source_report_generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_report_op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_report_run_id: Option<String>,
    pub requested_gene_ids: Vec<String>,
    pub matched_gene_ids: Vec<String>,
    pub missing_gene_ids: Vec<String>,
    pub selected_record_indices: Vec<usize>,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
    #[serde(default)]
    pub cohort_filter: RnaReadGeneSupportAuditCohortFilter,
    pub evaluated_row_count: usize,
    pub row_count: usize,
    pub accepted_target_record_indices: Vec<usize>,
    pub fragment_record_indices: Vec<usize>,
    pub complete_record_indices: Vec<usize>,
    pub complete_strict_record_indices: Vec<usize>,
    pub complete_exact_record_indices: Vec<usize>,
    pub rows: Vec<RnaReadGeneSupportAuditRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
/// Request parameters for summarizing TFBS density in one focus window relative
/// to a wider context window on the same sequence.
///
/// The current engine implementation counts TFBS features by overlap with the
/// requested spans. When `context_*` is omitted, the full sequence is used as
/// the wider comparison window.
pub struct TfbsRegionSummaryRequest {
    pub seq_id: String,
    pub focus_start_0based: usize,
    pub focus_end_0based_exclusive: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub context_start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub context_end_0based_exclusive: Option<usize>,
    #[serde(default = "default_tfbs_region_summary_min_focus_occurrences")]
    pub min_focus_occurrences: usize,
    #[serde(default)]
    pub min_context_occurrences: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub limit: Option<usize>,
}

fn default_tfbs_region_summary_min_focus_occurrences() -> usize {
    1
}

impl Default for TfbsRegionSummaryRequest {
    fn default() -> Self {
        Self {
            seq_id: String::new(),
            focus_start_0based: 0,
            focus_end_0based_exclusive: 0,
            context_start_0based: None,
            context_end_0based_exclusive: None,
            min_focus_occurrences: default_tfbs_region_summary_min_focus_occurrences(),
            min_context_occurrences: 0,
            limit: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One grouped TFBS summary row for a transcription factor across a focus and
/// context span.
pub struct TfbsRegionSummaryRow {
    pub tf_name: String,
    pub motif_ids: Vec<String>,
    pub focus_occurrences: usize,
    pub context_occurrences: usize,
    pub outside_focus_occurrences: usize,
    pub focus_density_per_kb: f64,
    pub context_density_per_kb: f64,
    pub outside_focus_density_per_kb: f64,
    pub focus_share_of_context_occurrences: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub focus_vs_context_density_ratio: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub focus_vs_outside_density_ratio: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable grouped TFBS summary for one focus window and one wider context.
pub struct TfbsRegionSummary {
    pub schema: String,
    pub seq_id: String,
    pub sequence_length_bp: usize,
    pub focus_start_0based: usize,
    pub focus_end_0based_exclusive: usize,
    pub context_start_0based: usize,
    pub context_end_0based_exclusive: usize,
    pub focus_width_bp: usize,
    pub context_width_bp: usize,
    pub outside_focus_width_bp: usize,
    pub total_tfbs_feature_count: usize,
    pub focus_hit_count: usize,
    pub context_hit_count: usize,
    pub matched_tf_count: usize,
    pub returned_tf_count: usize,
    pub min_focus_occurrences: usize,
    pub min_context_occurrences: usize,
    pub limit: usize,
    pub rows: Vec<TfbsRegionSummaryRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Deterministic background-normalization reference for one TF score track.
///
/// The reference distribution is computed with the same score family and the
/// same clipping semantics as the exported/displayed track itself so
/// `observed_peak_delta_from_p99` stays interpretable even when negative raw
/// scores are hidden from the presentation.
pub struct TfbsScoreTrackNormalizationReference {
    pub background_model: String,
    #[serde(default)]
    pub chance_model: String,
    pub random_sequence_length_bp: usize,
    pub random_seed: u64,
    pub sample_count: usize,
    pub mean_score: f64,
    pub stddev_score: f64,
    pub p95_score: f64,
    pub p99_score: f64,
    pub positive_fraction: f64,
    pub observed_peak_empirical_quantile: f64,
    #[serde(default)]
    pub observed_peak_modeled_quantile: f64,
    #[serde(default)]
    pub observed_peak_modeled_tail_probability: f64,
    #[serde(default)]
    pub observed_peak_modeled_tail_log10: f64,
    pub observed_peak_delta_from_p95: f64,
    pub observed_peak_delta_from_p99: f64,
    #[serde(default)]
    pub theoretical_min_score: f64,
    #[serde(default)]
    pub theoretical_max_score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One highlighted motif window from a continuous TF score track.
pub struct TfbsScoreTrackPeak {
    pub rank: usize,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub is_reverse: bool,
    pub score: f64,
    pub empirical_quantile: f64,
    pub delta_from_p99: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One pairwise synchrony estimate between two TF score tracks.
///
/// `raw_pearson` works on the exact displayed per-position vectors, while
/// `smoothed_pearson` uses the same vectors after centered boxcar smoothing.
/// The smoothed value is generally the better "are these peaks in the same
/// neighborhood?" measure for promoter interpretation.
pub struct TfbsScoreTrackCorrelationRow {
    pub left_tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub left_tf_name: Option<String>,
    pub right_tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub right_tf_name: Option<String>,
    pub overlap_window_count: usize,
    pub raw_pearson: f64,
    pub smoothed_pearson: f64,
    pub raw_spearman: f64,
    pub smoothed_spearman: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signed_primary_peak_offset_bp: Option<i64>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// One strand-specific TFBS score-track component used in cross-strand
/// correlation summaries.
pub enum TfbsScoreTrackStrandComponent {
    #[default]
    Forward,
    Reverse,
}

impl TfbsScoreTrackStrandComponent {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        }
    }

    pub fn short_label(self) -> &'static str {
        match self {
            Self::Forward => "F",
            Self::Reverse => "R",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Which strand-handling rule is used before comparing TFBS score tracks.
pub enum TfbsScoreTrackCorrelationSignalSource {
    #[default]
    MaxStrands,
    ForwardOnly,
    ReverseOnly,
}

impl TfbsScoreTrackCorrelationSignalSource {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::MaxStrands => "max_strands",
            Self::ForwardOnly => "forward_only",
            Self::ReverseOnly => "reverse_only",
        }
    }

    pub fn display_label(self) -> &'static str {
        match self {
            Self::MaxStrands => "Max(forward, reverse)",
            Self::ForwardOnly => "Forward strand only",
            Self::ReverseOnly => "Reverse strand only",
        }
    }

    pub fn summary_label(self) -> &'static str {
        match self {
            Self::MaxStrands => "max(forward_score, reverse_score)",
            Self::ForwardOnly => "forward_score",
            Self::ReverseOnly => "reverse_score",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Which pairwise statistic is used when presenting TFBS score-track
/// synchrony.
pub enum TfbsScoreTrackCorrelationMetric {
    #[default]
    Pearson,
    Spearman,
}

impl TfbsScoreTrackCorrelationMetric {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pearson => "pearson",
            Self::Spearman => "spearman",
        }
    }

    pub fn display_label(self, smoothed: bool) -> &'static str {
        match (self, smoothed) {
            (Self::Pearson, true) => "Smoothed Pearson r",
            (Self::Pearson, false) => "Raw Pearson r",
            (Self::Spearman, true) => "Smoothed Spearman rho",
            (Self::Spearman, false) => "Raw Spearman rho",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Shared correlation-sidecar for one TFBS score-track report.
pub struct TfbsScoreTrackCorrelationSummary {
    #[serde(default)]
    pub signal_source: TfbsScoreTrackCorrelationSignalSource,
    pub smoothing_method: String,
    pub smoothing_window_bp: usize,
    pub pair_count: usize,
    pub rows: Vec<TfbsScoreTrackCorrelationRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One explicit strand-pair synchrony estimate between two TF score tracks.
pub struct TfbsScoreTrackCrossStrandCorrelationCell {
    pub left_strand: TfbsScoreTrackStrandComponent,
    pub right_strand: TfbsScoreTrackStrandComponent,
    pub overlap_window_count: usize,
    pub raw_pearson: f64,
    pub smoothed_pearson: f64,
    pub raw_spearman: f64,
    pub smoothed_spearman: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signed_primary_peak_offset_bp: Option<i64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One TF-pair block containing all four `F-F / F-R / R-F / R-R` synchrony
/// estimates.
pub struct TfbsScoreTrackCrossStrandCorrelationRow {
    pub left_tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub left_tf_name: Option<String>,
    pub right_tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub right_tf_name: Option<String>,
    #[serde(default)]
    pub cells: Vec<TfbsScoreTrackCrossStrandCorrelationCell>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Shared cross-strand TFBS synchrony summary where each TF-pair carries all
/// four strand pairings in one row so renderers can either show grouped 2x2
/// TF-pair blocks or expand the matrix into explicit `F, R` curve axes.
pub struct TfbsScoreTrackCrossStrandCorrelationSummary {
    pub smoothing_method: String,
    pub smoothing_window_bp: usize,
    pub pair_count: usize,
    #[serde(default)]
    pub rows: Vec<TfbsScoreTrackCrossStrandCorrelationRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// How one motif's forward and reverse strand score tracks relate to each
/// other over the selected span.
pub struct TfbsScoreTrackDirectionalSummary {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub forward_primary_peak_position_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reverse_primary_peak_position_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signed_primary_peak_offset_bp: Option<i64>,
    pub raw_pearson: f64,
    pub smoothed_pearson: f64,
    pub raw_spearman: f64,
    pub smoothed_spearman: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Which anchor-vs-candidate statistic is used when ranking TFBS score-track
/// similarity for one selected DNA span.
pub enum TfbsTrackSimilarityRankingMetric {
    RawPearson,
    SmoothedPearson,
    RawSpearman,
    #[default]
    SmoothedSpearman,
}

impl TfbsTrackSimilarityRankingMetric {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::RawPearson => "raw_pearson",
            Self::SmoothedPearson => "smoothed_pearson",
            Self::RawSpearman => "raw_spearman",
            Self::SmoothedSpearman => "smoothed_spearman",
        }
    }

    pub fn display_label(self) -> &'static str {
        match self {
            Self::RawPearson => "Raw Pearson r",
            Self::SmoothedPearson => "Smoothed Pearson r",
            Self::RawSpearman => "Raw Spearman rho",
            Self::SmoothedSpearman => "Smoothed Spearman rho",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One candidate TF ranked against one anchor TF over the same displayed
/// score-track span.
pub struct TfbsTrackSimilarityRow {
    pub candidate_tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub candidate_tf_name: Option<String>,
    pub overlap_window_count: usize,
    pub raw_pearson: f64,
    pub smoothed_pearson: f64,
    pub raw_spearman: f64,
    pub smoothed_spearman: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signed_primary_peak_offset_bp: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_summary: Option<JasparCatalogRemoteSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Shared anchor-vs-candidate TFBS score-track ranking over one selected DNA
/// span.
pub struct TfbsTrackSimilarityReport {
    pub schema: String,
    #[serde(default)]
    pub target_kind: String,
    #[serde(default)]
    pub target_label: String,
    pub seq_id: String,
    pub source_sequence_length_bp: usize,
    #[serde(default)]
    pub scan_topology: InlineSequenceTopology,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub score_kind: TfbsScoreTrackValueKind,
    pub view_start_0based: usize,
    pub view_end_0based_exclusive: usize,
    pub clip_negative: bool,
    pub anchor_requested: String,
    pub anchor_tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchor_tf_name: Option<String>,
    #[serde(default)]
    pub candidate_scope: String,
    #[serde(default)]
    pub candidates_requested: Vec<String>,
    #[serde(default)]
    pub species_filters: Vec<String>,
    #[serde(default)]
    pub include_remote_metadata: bool,
    #[serde(default)]
    pub ranking_metric: TfbsTrackSimilarityRankingMetric,
    pub scanned_candidate_count: usize,
    pub returned_candidate_count: usize,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub rows: Vec<TfbsTrackSimilarityRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One continuous TF motif score track over a requested DNA span.
///
/// `forward_scores[i]` and `reverse_scores[i]` correspond to the motif window
/// that starts at `track_start_0based + i`. Scores may be clipped to `0.0`
/// when the parent report requests positive-only display.
pub struct TfbsScoreTrackRow {
    pub tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_name: Option<String>,
    pub motif_length_bp: usize,
    #[serde(default)]
    pub motif_logo_columns: Vec<JasparExpertColumn>,
    pub track_start_0based: usize,
    pub scored_window_count: usize,
    pub max_score: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_position_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub normalization_reference: Option<TfbsScoreTrackNormalizationReference>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub directional_summary: Option<TfbsScoreTrackDirectionalSummary>,
    #[serde(default)]
    pub top_peaks: Vec<TfbsScoreTrackPeak>,
    #[serde(default)]
    pub forward_scores: Vec<f64>,
    #[serde(default)]
    pub reverse_scores: Vec<f64>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Which per-window TF motif score is carried by a continuous score-track
/// report.
pub enum TfbsScoreTrackValueKind {
    #[default]
    LlrBits,
    LlrQuantile,
    LlrBackgroundQuantile,
    LlrBackgroundTailLog10,
    TrueLogOddsBits,
    TrueLogOddsQuantile,
    TrueLogOddsBackgroundQuantile,
    TrueLogOddsBackgroundTailLog10,
}

impl TfbsScoreTrackValueKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::LlrBits => "llr_bits",
            Self::LlrQuantile => "llr_quantile",
            Self::LlrBackgroundQuantile => "llr_background_quantile",
            Self::LlrBackgroundTailLog10 => "llr_background_tail_log10",
            Self::TrueLogOddsBits => "true_log_odds_bits",
            Self::TrueLogOddsQuantile => "true_log_odds_quantile",
            Self::TrueLogOddsBackgroundQuantile => "true_log_odds_background_quantile",
            Self::TrueLogOddsBackgroundTailLog10 => "true_log_odds_background_tail_log10",
        }
    }

    pub fn supports_negative_values(self) -> bool {
        matches!(self, Self::LlrBits | Self::TrueLogOddsBits)
    }

    pub fn uses_llr_background_bits(self) -> bool {
        matches!(
            self,
            Self::LlrBits
                | Self::LlrQuantile
                | Self::LlrBackgroundQuantile
                | Self::LlrBackgroundTailLog10
        )
    }

    pub fn uses_true_log_odds_background_bits(self) -> bool {
        matches!(
            self,
            Self::TrueLogOddsBits
                | Self::TrueLogOddsQuantile
                | Self::TrueLogOddsBackgroundQuantile
                | Self::TrueLogOddsBackgroundTailLog10
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One transcription-start marker relevant to one TFBS score-track span.
pub struct TfbsScoreTrackTssMarker {
    pub feature_id: usize,
    pub feature_kind: String,
    pub label: String,
    pub position_0based: usize,
    pub is_reverse: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One external interval projected into the selected TFBS score-track span.
pub struct TfbsScoreTrackOverlayInterval {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub score: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One external interval-track lane rendered underneath the motif score rows.
pub struct TfbsScoreTrackOverlayTrack {
    pub source_kind: String,
    pub track_name: String,
    pub display_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_file_name: Option<String>,
    pub interval_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_score: Option<f64>,
    #[serde(default)]
    pub intervals: Vec<TfbsScoreTrackOverlayInterval>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable per-position TF motif score tracks for promoter-design review.
pub struct TfbsScoreTrackReport {
    pub schema: String,
    #[serde(default)]
    pub target_kind: String,
    #[serde(default)]
    pub target_label: String,
    pub seq_id: String,
    pub source_sequence_length_bp: usize,
    pub sequence_length_bp: usize,
    #[serde(default)]
    pub scan_topology: InlineSequenceTopology,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default)]
    pub score_kind: TfbsScoreTrackValueKind,
    pub view_start_0based: usize,
    pub view_end_0based_exclusive: usize,
    pub clip_negative: bool,
    #[serde(default)]
    pub motifs_requested: Vec<String>,
    pub global_max_score: f64,
    #[serde(default)]
    pub tss_markers: Vec<TfbsScoreTrackTssMarker>,
    #[serde(default)]
    pub overlay_tracks: Vec<TfbsScoreTrackOverlayTrack>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub correlation_summary: Option<TfbsScoreTrackCorrelationSummary>,
    #[serde(default)]
    pub correlation_summaries: Vec<TfbsScoreTrackCorrelationSummary>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cross_strand_correlation_summary: Option<TfbsScoreTrackCrossStrandCorrelationSummary>,
    #[serde(default)]
    pub tracks: Vec<TfbsScoreTrackRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// One requested gene token for a multi-gene promoter TFBS analysis.
pub struct PromoterTfbsGeneQuery {
    pub gene_query: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub occurrence: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub display_label: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One promoter-resolved gene entry inside a multi-gene TFBS analysis report.
pub struct MultiGenePromoterTfbsGeneReport {
    pub gene_query: String,
    pub occurrence: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_requested: Option<String>,
    pub display_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_name: Option<String>,
    pub transcript_id: String,
    pub chromosome: String,
    pub strand: String,
    pub promoter_start_1based: usize,
    pub promoter_end_1based: usize,
    pub promoter_length_bp: usize,
    pub tss_1based: usize,
    pub sequence_orientation: String,
    pub used_fuzzy_gene_match: bool,
    pub tfbs_score_tracks: TfbsScoreTrackReport,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One per-gene/per-motif summary row derived from a multi-gene promoter TFBS
/// analysis.
pub struct MultiGenePromoterTfbsSummaryRow {
    pub gene_label: String,
    pub gene_query: String,
    pub transcript_id: String,
    pub tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_name: Option<String>,
    pub max_score: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_position_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_position_promoter_relative_bp: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_genomic_position_1based: Option<usize>,
    pub positive_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable multi-gene promoter TFBS analysis report with one promoter-aligned
/// score-track panel per requested gene and one compact comparison matrix.
pub struct MultiGenePromoterTfbsReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub genome_id: String,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub score_kind: TfbsScoreTrackValueKind,
    pub clip_negative: bool,
    #[serde(default)]
    pub motifs_requested: Vec<String>,
    #[serde(default)]
    pub gene_queries_requested: Vec<PromoterTfbsGeneQuery>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_set: Option<GeneSetRequest>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_set_resolution: Option<GeneSetResolutionReport>,
    pub returned_gene_count: usize,
    #[serde(default)]
    pub genes: Vec<MultiGenePromoterTfbsGeneReport>,
    #[serde(default)]
    pub summary_rows: Vec<MultiGenePromoterTfbsSummaryRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// First-slice cohort relationship for promoter comparison.
pub enum PromoterCohortKind {
    #[default]
    Manual,
    CoRegulated,
    AntiCoRegulated,
}

impl PromoterCohortKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Manual => "manual",
            Self::CoRegulated => "co_regulated",
            Self::AntiCoRegulated => "anti_co_regulated",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One promoter window resolved for a promoter cohort comparison.
pub struct PromoterCohortResolvedWindow {
    pub gene_label: String,
    pub gene_query: String,
    pub transcript_id: String,
    pub chromosome: String,
    pub strand: String,
    pub promoter_start_1based: usize,
    pub promoter_end_1based: usize,
    pub promoter_length_bp: usize,
    pub tss_1based: usize,
    pub tss_position_0based: usize,
    pub sequence_orientation: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One pairwise promoter similarity row derived from shared TFBS score tracks.
pub struct PromoterCohortPairwiseSimilarity {
    pub left_gene_label: String,
    pub right_gene_label: String,
    pub shared_motif_count: usize,
    pub mean_raw_pearson: f64,
    pub mean_smoothed_spearman: f64,
    #[serde(default)]
    pub motif_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Non-blocking flag comparing available promoter evidence with the declared
/// cohort relationship expectation.
pub struct PromoterCohortRelationshipFlag {
    pub flag_kind: String,
    pub evidence_kind: String,
    #[serde(default)]
    pub gene_labels: Vec<String>,
    pub detail: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Summary of a motif peak that is shared across or specific to cohort members.
pub struct PromoterCohortTfbsPeakSummary {
    pub tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_name: Option<String>,
    pub promoter_count: usize,
    #[serde(default)]
    pub gene_labels: Vec<String>,
    pub max_score: f64,
    #[serde(default)]
    pub peak_positions_promoter_relative_bp: Vec<i64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Expression metadata echoed into a promoter cohort comparison.
pub struct PromoterCohortExpressionAssociation {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    pub value: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    pub source: String,
    pub association_note: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Engine-owned first-slice promoter cohort comparison.
pub struct PromoterCohortComparisonReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub cohort_label: String,
    pub cohort_kind: PromoterCohortKind,
    pub genome_id: String,
    #[serde(default)]
    pub source_seq_ids: Vec<String>,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub score_kind: TfbsScoreTrackValueKind,
    pub clip_negative: bool,
    #[serde(default)]
    pub motifs_requested: Vec<String>,
    #[serde(default)]
    pub gene_queries_requested: Vec<PromoterTfbsGeneQuery>,
    pub resolved_promoter_count: usize,
    #[serde(default)]
    pub resolved_promoter_windows: Vec<PromoterCohortResolvedWindow>,
    #[serde(default)]
    pub tfbs_score_track_summaries: Vec<MultiGenePromoterTfbsSummaryRow>,
    #[serde(default)]
    pub pairwise_similarity: Vec<PromoterCohortPairwiseSimilarity>,
    #[serde(default)]
    pub relationship_flags: Vec<PromoterCohortRelationshipFlag>,
    #[serde(default)]
    pub shared_tfbs_peaks: Vec<PromoterCohortTfbsPeakSummary>,
    #[serde(default)]
    pub cohort_specific_tfbs_peaks: Vec<PromoterCohortTfbsPeakSummary>,
    #[serde(default)]
    pub expression_associations: Vec<PromoterCohortExpressionAssociation>,
    #[serde(default)]
    pub cutrun_dataset_ids: Vec<String>,
    #[serde(default)]
    pub cutrun_read_report_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub multi_gene_promoter_tfbs: Option<MultiGenePromoterTfbsReport>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Filter over UCSC RepeatMasker (`rmsk`) annotations.
pub struct RepeatAnnotationFilter {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub rep_names: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub rep_classes: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub rep_families: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub normalized_aliases: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub span_start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub span_end_0based_exclusive: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One RepeatMasker annotation row normalized from a UCSC `rmsk` table.
pub struct RepeatAnnotationRecord {
    pub annotation_id: String,
    pub genome_id: String,
    pub chromosome: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: String,
    pub rep_name: String,
    pub rep_class: String,
    pub rep_family: String,
    pub normalized_alias: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub score: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub milli_div: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_line_number: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable result of filtering one RepeatMasker source file.
pub struct RepeatAnnotationQueryReport {
    pub schema: String,
    pub genome_id: String,
    pub rmsk_path: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub filter: RepeatAnnotationFilter,
    pub parsed_row_count: usize,
    pub matched_row_count: usize,
    pub returned_row_count: usize,
    pub malformed_line_count: usize,
    #[serde(default)]
    pub malformed_examples: Vec<String>,
    #[serde(default)]
    pub rows: Vec<RepeatAnnotationRecord>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One UCSC `rmsk` interval projected onto a GENtle sequence anchor.
pub struct SequenceRepeatOverlapRow {
    pub repeat: RepeatAnnotationRecord,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub local_strand: String,
    pub genomic_start_0based: usize,
    pub genomic_end_0based_exclusive: usize,
    pub overlap_bp: usize,
    pub clipped: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Aggregate repeat-overlap metrics for one class/family/alias bucket.
pub struct RepeatOverlapSummaryRow {
    pub rep_class: String,
    pub rep_family: String,
    pub normalized_alias: String,
    pub repeat_count: usize,
    pub overlap_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Repeat-overlap lookup over one genome-anchored sequence and prepared rmsk index.
pub struct SequenceRepeatOverlapReport {
    pub schema: String,
    pub seq_id: String,
    pub rmsk_index_path: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub genome_anchor: Option<SequenceGenomeAnchorSummary>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_end_0based_exclusive: Option<usize>,
    pub query_length_bp: usize,
    pub matched_repeat_count: usize,
    pub returned_repeat_count: usize,
    pub total_overlap_bp: usize,
    pub covered_query_bp: usize,
    pub coverage_fraction: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nearest_repeat_distance_bp: Option<usize>,
    #[serde(default)]
    pub class_summaries: Vec<RepeatOverlapSummaryRow>,
    #[serde(default)]
    pub rows: Vec<SequenceRepeatOverlapRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Result of attaching UCSC rmsk overlaps as ordinary sequence features.
pub struct RepeatFeatureMaterializationReport {
    pub schema: String,
    pub seq_id: String,
    pub rmsk_index_path: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub genome_anchor: Option<SequenceGenomeAnchorSummary>,
    pub matched_repeat_count: usize,
    pub added_feature_count: usize,
    pub skipped_existing_count: usize,
    pub removed_existing_count: usize,
    #[serde(default)]
    pub feature_ids: Vec<usize>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Coordinate frame used to derive one repeat-environment window.
pub enum RepeatEnvironmentGeometryMode {
    #[default]
    RepeatMidpoint,
    Transcript5utrStart,
    Pol2PromoterUpstream,
    CdsStopContext,
}

impl RepeatEnvironmentGeometryMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::RepeatMidpoint => "repeat_midpoint",
            Self::Transcript5utrStart => "transcript_5utr_start",
            Self::Pol2PromoterUpstream => "pol2_promoter_upstream",
            Self::CdsStopContext => "cds_stop_context",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One transcript/gene projection attached to a repeat-cohort row.
pub struct RepeatTranscriptContext {
    pub gene_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    pub transcript_id: String,
    pub chromosome: String,
    pub strand: String,
    pub transcript_start_1based: usize,
    pub transcript_end_1based: usize,
    pub overlaps_repeat: bool,
    pub repeat_relation: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tss_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inferred_5utr_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inferred_5utr_end_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cds_stop_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signed_tss_distance_bp: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signed_stop_distance_bp: Option<i64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Availability and selected coordinates for one geometry mode on one row.
pub struct RepeatEnvironmentGeometryWindow {
    pub mode: RepeatEnvironmentGeometryMode,
    pub available: bool,
    pub reason: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub end_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchor_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One repeat-locus environment row with transcript-aware geometry metadata.
pub struct RepeatEnvironmentCohortRow {
    pub row_id: String,
    pub repeat: RepeatAnnotationRecord,
    pub selected_geometry: RepeatEnvironmentGeometryMode,
    pub selected_window: RepeatEnvironmentGeometryWindow,
    #[serde(default)]
    pub geometry_windows: Vec<RepeatEnvironmentGeometryWindow>,
    #[serde(default)]
    pub transcript_contexts: Vec<RepeatTranscriptContext>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub primary_gene_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub primary_transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rna_rank_score: Option<f64>,
    #[serde(default)]
    pub rna_support_labels: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Repeat-family cohort with multiple inspectable genomic coordinate frames.
pub struct RepeatEnvironmentCohortReport {
    pub schema: String,
    pub genome_id: String,
    pub rmsk_path: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub filter: RepeatAnnotationFilter,
    pub selected_geometry: RepeatEnvironmentGeometryMode,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub parsed_repeat_count: usize,
    pub matched_repeat_count: usize,
    pub returned_row_count: usize,
    #[serde(default)]
    pub rows: Vec<RepeatEnvironmentCohortRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Per-window TFBS score-track report inside a generic cohort comparison.
pub struct WindowCohortTfbsWindowReport {
    pub row_id: String,
    pub label: String,
    pub chromosome: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub anchor_1based: Option<usize>,
    pub geometry_mode: RepeatEnvironmentGeometryMode,
    pub tfbs_score_tracks: TfbsScoreTrackReport,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Compact per-window/per-motif summary row for a generic cohort TFBS report.
pub struct WindowCohortTfbsSummaryRow {
    pub row_id: String,
    pub label: String,
    pub tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_name: Option<String>,
    pub max_score: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_position_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_genomic_position_1based: Option<usize>,
    pub positive_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Generic multi-window TFBS comparison report for stored genomic cohorts.
pub struct WindowCohortTfbsReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub genome_id: String,
    pub cohort_schema: String,
    pub geometry_mode: RepeatEnvironmentGeometryMode,
    pub score_kind: TfbsScoreTrackValueKind,
    pub clip_negative: bool,
    #[serde(default)]
    pub motifs_requested: Vec<String>,
    pub returned_window_count: usize,
    #[serde(default)]
    pub windows: Vec<WindowCohortTfbsWindowReport>,
    #[serde(default)]
    pub summary_rows: Vec<WindowCohortTfbsSummaryRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Default,
)]
#[serde(rename_all = "snake_case")]
/// Topology hint for inline sequence operands used by state-optional scans.
pub enum InlineSequenceTopology {
    #[default]
    Linear,
    Circular,
}

impl InlineSequenceTopology {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::Circular => "circular",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "kind", rename_all = "snake_case")]
/// Shared operand for sequence inspections that may run against a stored
/// sequence or an inline ASCII DNA payload.
pub enum SequenceScanTarget {
    SeqId {
        seq_id: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        span_start_0based: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        span_end_0based_exclusive: Option<usize>,
    },
    InlineSequence {
        sequence_text: String,
        #[serde(default)]
        topology: InlineSequenceTopology,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        id_hint: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        span_start_0based: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        span_end_0based_exclusive: Option<usize>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One restriction-site hit returned by a state-optional restriction scan.
pub struct RestrictionSiteScanHit {
    pub enzyme_name: String,
    pub recognition_sequence: String,
    pub recognition_start_0based: usize,
    pub recognition_end_0based_exclusive: usize,
    pub source_recognition_start_0based: usize,
    pub source_recognition_end_0based_exclusive: usize,
    pub recognition_length_bp: usize,
    pub forward_strand: bool,
    pub end_geometry: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub forward_cut_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reverse_cut_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub opening_start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub opening_end_0based_exclusive: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_forward_cut_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_reverse_cut_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_opening_start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_opening_end_0based_exclusive: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable result payload for direct restriction-site inspection on either a
/// stored `seq_id` or inline ASCII DNA input.
pub struct RestrictionSiteScanReport {
    pub schema: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub report_id: String,
    pub target_kind: String,
    pub target_label: String,
    pub source_sequence_length_bp: usize,
    pub scan_start_0based: usize,
    pub scan_end_0based_exclusive: usize,
    pub scan_length_bp: usize,
    #[serde(default)]
    pub scan_topology: InlineSequenceTopology,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default)]
    pub enzyme_filters: Vec<String>,
    #[serde(default)]
    pub enzymes_scanned: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_sites_per_enzyme: Option<usize>,
    pub include_cut_geometry: bool,
    pub matched_site_count: usize,
    #[serde(default)]
    pub skipped_enzyme_names_due_to_max_sites: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    #[serde(default)]
    pub rows: Vec<RestrictionSiteScanHit>,
}

pub const PROJECT_FACT_GRAPH_SCHEMA: &str = "gentle.project_fact_graph.v1";
pub const FACT_EXPRESSION_SCHEMA: &str = "gentle.fact_expression.v1";
pub const FACT_EVALUATION_SCHEMA: &str = "gentle.fact_evaluation.v1";

fn project_fact_graph_schema() -> String {
    PROJECT_FACT_GRAPH_SCHEMA.to_string()
}

fn fact_evaluation_schema() -> String {
    FACT_EVALUATION_SCHEMA.to_string()
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[serde(rename_all = "snake_case")]
/// Top-level project object class addressed by a fact atom.
pub enum FactSubjectKind {
    Sequence,
    Report,
    Container,
    Ui,
    Other,
}

impl FactSubjectKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Sequence => "sequence",
            Self::Report => "report",
            Self::Container => "container",
            Self::Ui => "ui",
            Self::Other => "other",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// World assumption used when evaluating a known project fact type.
pub enum ProjectFactWorld {
    ClosedWorld,
    OpenWorld,
}

impl ProjectFactWorld {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ClosedWorld => "closed_world",
            Self::OpenWorld => "open_world",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
/// Registry entry for one controlled project-fact vocabulary item.
pub struct ProjectFactTypeSpec {
    pub name: &'static str,
    pub world: ProjectFactWorld,
    pub requires_basis: bool,
    pub subject_kind: FactSubjectKind,
    pub description: &'static str,
}

const PROJECT_FACT_TYPE_SPECS: &[ProjectFactTypeSpec] = &[
    ProjectFactTypeSpec {
        name: "sequence.exists",
        world: ProjectFactWorld::ClosedWorld,
        requires_basis: false,
        subject_kind: FactSubjectKind::Sequence,
        description: "A sequence with this id is loaded in the current project state.",
    },
    ProjectFactTypeSpec {
        name: "sequence.kind",
        world: ProjectFactWorld::ClosedWorld,
        requires_basis: false,
        subject_kind: FactSubjectKind::Sequence,
        description: "Loaded sequence molecule class; use equals \"dna\", \"rna\", or \"protein\".",
    },
    ProjectFactTypeSpec {
        name: "sequence.length",
        world: ProjectFactWorld::ClosedWorld,
        requires_basis: false,
        subject_kind: FactSubjectKind::Sequence,
        description: "Loaded sequence length in residues/bases; use compare for numeric predicates.",
    },
    ProjectFactTypeSpec {
        name: "sequence.circular",
        world: ProjectFactWorld::ClosedWorld,
        requires_basis: false,
        subject_kind: FactSubjectKind::Sequence,
        description: "Whether the loaded sequence is circular.",
    },
    ProjectFactTypeSpec {
        name: "report.exists",
        world: ProjectFactWorld::OpenWorld,
        requires_basis: true,
        subject_kind: FactSubjectKind::Report,
        description: "An explicit proof/report artifact is available as evidence.",
    },
    ProjectFactTypeSpec {
        name: "restriction_site.present",
        world: ProjectFactWorld::OpenWorld,
        requires_basis: true,
        subject_kind: FactSubjectKind::Sequence,
        description: "A restriction-site scan report proves at least one matching enzyme site in the requested range.",
    },
    ProjectFactTypeSpec {
        name: "restriction_site.absent",
        world: ProjectFactWorld::OpenWorld,
        requires_basis: true,
        subject_kind: FactSubjectKind::Sequence,
        description: "A covering restriction-site scan report proves zero matching enzyme sites in the requested range.",
    },
    ProjectFactTypeSpec {
        name: "ui.host_available",
        world: ProjectFactWorld::ClosedWorld,
        requires_basis: false,
        subject_kind: FactSubjectKind::Ui,
        description: "The current adapter can ask the user to perform GUI-hosted actions such as picking a file.",
    },
];

pub fn project_fact_type_specs() -> &'static [ProjectFactTypeSpec] {
    PROJECT_FACT_TYPE_SPECS
}

pub fn project_fact_type_spec(name: &str) -> Option<&'static ProjectFactTypeSpec> {
    PROJECT_FACT_TYPE_SPECS
        .iter()
        .find(|spec| spec.name == name)
}

pub fn project_fact_registry_prompt_block() -> String {
    let mut lines = vec![
        format!("Known project fact vocabulary ({FACT_EXPRESSION_SCHEMA}):"),
        "Use only these fact names in precondition_expr/expected_effects unless you are deliberately asking GENtle to treat a future fact as unknown.".to_string(),
    ];
    for spec in PROJECT_FACT_TYPE_SPECS {
        lines.push(format!(
            "- {}: subject_kind={}; world={}; requires_basis={}; {}",
            spec.name,
            spec.subject_kind.as_str(),
            spec.world.as_str(),
            spec.requires_basis,
            spec.description
        ));
    }
    lines.join("\n")
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
/// Stable subject reference shared by projected facts and fact expressions.
pub struct FactSubject {
    pub kind: FactSubjectKind,
    pub id: String,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[serde(tag = "kind", rename_all = "snake_case")]
/// Coordinate scope covered or required by a project fact.
///
/// Deserialization also accepts the legacy string `"whole_sequence"` as sugar
/// for `{"kind":"whole_sequence"}`.
pub enum FactRange {
    WholeSequence,
    Span {
        start_0based: usize,
        end_0based_exclusive: usize,
        #[serde(default)]
        topology: InlineSequenceTopology,
    },
}

impl Default for FactRange {
    fn default() -> Self {
        Self::WholeSequence
    }
}

impl<'de> Deserialize<'de> for FactRange {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        if let Some(text) = value.as_str()
            && text == "whole_sequence"
        {
            return Ok(Self::WholeSequence);
        }
        #[derive(Deserialize)]
        #[serde(tag = "kind", rename_all = "snake_case")]
        enum Helper {
            WholeSequence,
            Span {
                start_0based: usize,
                end_0based_exclusive: usize,
                #[serde(default)]
                topology: InlineSequenceTopology,
            },
        }
        match Helper::deserialize(value).map_err(serde::de::Error::custom)? {
            Helper::WholeSequence => Ok(Self::WholeSequence),
            Helper::Span {
                start_0based,
                end_0based_exclusive,
                topology,
            } => Ok(Self::Span {
                start_0based,
                end_0based_exclusive,
                topology,
            }),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
/// Evidence pointer that turns an open-world observation into a proof-backed
/// fact.
pub struct FactBasis {
    pub report_id: String,
    pub report_kind: String,
    pub evidence_class: EvidenceClass,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
}

impl Default for FactBasis {
    fn default() -> Self {
        Self {
            report_id: String::new(),
            report_kind: String::new(),
            evidence_class: EvidenceClass::HardFact,
            op_id: None,
            run_id: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
/// One ground fact projected from project state or explicit evidence.
pub struct ProjectFact {
    pub fact: String,
    pub subject: FactSubject,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub enzyme: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub range: Option<FactRange>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub value: Option<Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub basis: Option<FactBasis>,
}

impl Default for ProjectFact {
    fn default() -> Self {
        Self {
            fact: String::new(),
            subject: FactSubject {
                kind: FactSubjectKind::Other,
                id: String::new(),
            },
            enzyme: None,
            range: None,
            value: None,
            basis: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(default)]
/// Deterministic projection of current project state plus supplied proof
/// evidence into controlled-vocabulary facts.
pub struct ProjectFactGraph {
    #[serde(default = "project_fact_graph_schema")]
    pub schema: String,
    pub facts: Vec<ProjectFact>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(default)]
/// Atomic fact pattern used inside a fact expression.
pub struct FactAtom {
    pub fact: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subject: Option<FactSubject>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub enzyme: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub range: Option<FactRange>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub equals: Option<Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub compare: Option<FactComparison>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
/// Numeric comparison predicate for scalar fact values.
pub struct FactComparison {
    pub op: String,
    pub value: Value,
}

impl Default for FactComparison {
    fn default() -> Self {
        Self {
            op: "eq".to_string(),
            value: Value::Null,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(untagged)]
/// Boolean fact-expression tree. Evaluation uses three-valued Kleene logic.
pub enum FactExpression {
    All { all: Vec<FactExpression> },
    Any { any: Vec<FactExpression> },
    Not { not: Box<FactExpression> },
    Atom(FactAtom),
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Three-valued fact-expression result.
pub enum FactTruth {
    #[default]
    Unknown,
    Satisfied,
    Unsatisfied,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(default)]
/// Deterministic result of evaluating one fact expression against a fact graph.
pub struct FactEvaluationResult {
    #[serde(default = "fact_evaluation_schema")]
    pub schema: String,
    pub truth: FactTruth,
    pub unmet_atoms: Vec<FactAtom>,
    pub unknown_atoms: Vec<FactAtom>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(default)]
/// One TFBS/JASPAR motif hit returned by a non-mutating direct scan on either a
/// stored `seq_id` or inline ASCII DNA input.
pub struct TfbsHitScanRow {
    pub tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_name: Option<String>,
    pub motif_consensus_iupac: String,
    pub motif_length_bp: usize,
    pub match_start_0based: usize,
    pub match_end_0based_exclusive: usize,
    pub source_match_start_0based: usize,
    pub source_match_end_0based_exclusive: usize,
    #[serde(default)]
    pub wraps_origin: bool,
    pub forward_strand: bool,
    pub matched_sequence: String,
    pub llr_bits: f64,
    pub llr_quantile: f64,
    pub true_log_odds_bits: f64,
    pub true_log_odds_quantile: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(default)]
/// Portable result payload for state-optional TFBS/JASPAR hit scans on either a
/// stored `seq_id` or inline ASCII DNA input.
pub struct TfbsHitScanReport {
    pub schema: String,
    pub target_kind: String,
    pub target_label: String,
    pub source_sequence_length_bp: usize,
    pub scan_start_0based: usize,
    pub scan_end_0based_exclusive: usize,
    pub scan_length_bp: usize,
    #[serde(default)]
    pub scan_topology: InlineSequenceTopology,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default)]
    pub motifs_requested: Vec<String>,
    #[serde(default)]
    pub motifs_scanned: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub default_min_llr_bits: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub default_min_llr_quantile: Option<f64>,
    #[serde(default)]
    pub per_tf_thresholds: Vec<TfThresholdOverride>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_hits: Option<usize>,
    pub truncated_at_max_hits: bool,
    pub matched_hit_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    #[serde(default)]
    pub rows: Vec<TfbsHitScanRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Deterministic score-distribution summary for one JASPAR motif over one
/// pseudorandom DNA background.
pub struct JasparScoreDistributionSummary {
    pub sample_count: usize,
    pub min_score: f64,
    pub max_score: f64,
    pub mean_score: f64,
    pub stddev_score: f64,
    pub p01_score: f64,
    pub p05_score: f64,
    pub p25_score: f64,
    pub p50_score: f64,
    pub p75_score: f64,
    pub p95_score: f64,
    pub p99_score: f64,
    pub positive_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One JASPAR entry presented with score-extreme sequences and one
/// deterministic random-background score summary.
pub struct JasparEntryPresentationRow {
    pub motif_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub motif_name: Option<String>,
    pub consensus_iupac: String,
    pub motif_length_bp: usize,
    pub maximizing_sequence: String,
    pub minimizing_sequence: String,
    pub maximizing_llr_bits: f64,
    pub maximizing_llr_quantile: f64,
    pub minimizing_llr_bits: f64,
    pub minimizing_llr_quantile: f64,
    pub maximizing_true_log_odds_bits: f64,
    pub maximizing_true_log_odds_quantile: f64,
    pub minimizing_true_log_odds_bits: f64,
    pub minimizing_true_log_odds_quantile: f64,
    pub llr_bits_distribution: JasparScoreDistributionSummary,
    pub true_log_odds_bits_distribution: JasparScoreDistributionSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable registry-wide presentation of JASPAR entries that shows the
/// deterministic best/worst-scoring motif strings and one random-background
/// score expectation summary per entry.
pub struct JasparEntryPresentationReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default)]
    pub requested_motifs: Vec<String>,
    pub registry_entry_count: usize,
    pub resolved_entry_count: usize,
    pub random_sequence_length_bp: usize,
    pub random_seed: u64,
    pub background_model: String,
    #[serde(default)]
    pub rows: Vec<JasparEntryPresentationRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One top-ranking motif row inside a registry-wide JASPAR benchmark summary.
pub struct JasparRegistryBenchmarkTopRow {
    pub motif_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub motif_name: Option<String>,
    pub value: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One score-family aggregate summary across a registry-wide JASPAR benchmark.
pub struct JasparRegistryScoreFamilySummary {
    pub score_kind: TfbsScoreTrackValueKind,
    pub label: String,
    pub motif_count: usize,
    pub global_min_score: f64,
    pub global_max_score: f64,
    pub mean_of_mean_scores: f64,
    pub mean_of_stddev_scores: f64,
    pub median_of_p50_scores: f64,
    pub mean_positive_fraction: f64,
    #[serde(default)]
    pub top_max_score_rows: Vec<JasparRegistryBenchmarkTopRow>,
    #[serde(default)]
    pub top_positive_fraction_rows: Vec<JasparRegistryBenchmarkTopRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Registry-wide deterministic JASPAR benchmark snapshot over one shared
/// pseudorandom DNA background.
pub struct JasparRegistryBenchmarkReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub registry_entry_count: usize,
    pub benchmarked_entry_count: usize,
    pub random_sequence_length_bp: usize,
    pub random_seed: u64,
    pub background_model: String,
    #[serde(default)]
    pub score_family_summaries: Vec<JasparRegistryScoreFamilySummary>,
    #[serde(default)]
    pub rows: Vec<JasparEntryPresentationRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Compact remote JASPAR metadata summary suitable for catalog/list views.
pub struct JasparCatalogRemoteSummary {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub collection: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tax_group: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_class: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_family: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub data_type: Option<String>,
    pub species_count: usize,
    #[serde(default)]
    pub species_preview: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One row inside the shared JASPAR catalog report.
pub struct JasparCatalogRow {
    pub motif_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub motif_name: Option<String>,
    pub consensus_iupac: String,
    pub motif_length_bp: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_summary: Option<JasparCatalogRemoteSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable catalog/list view of local JASPAR entries with optional remote
/// metadata summaries for the returned subset.
pub struct JasparCatalogReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub filter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub limit: Option<usize>,
    pub include_remote_metadata: bool,
    pub registry_entry_count: usize,
    pub returned_entry_count: usize,
    #[serde(default)]
    pub rows: Vec<JasparCatalogRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One motif matched while resolving a user-facing TF query or TF-group token.
pub struct TfQueryResolutionMatch {
    pub motif_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub motif_name: Option<String>,
    pub consensus_iupac: String,
    pub motif_length_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Resolution result for one user-facing TF query token.
pub struct TfQueryResolutionEntry {
    pub query: String,
    pub normalized_query: String,
    pub resolution_kind: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default)]
    pub aliases: Vec<String>,
    #[serde(default)]
    pub matches: Vec<TfQueryResolutionMatch>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unresolved_reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable report for resolving user-facing TF/group/family queries into local
/// motif-registry entries.
pub struct TfQueryResolutionReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default)]
    pub queries: Vec<String>,
    pub returned_query_count: usize,
    pub matched_motif_count: usize,
    #[serde(default)]
    pub entries: Vec<TfQueryResolutionEntry>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One species assignment reported by optional remote JASPAR metadata.
pub struct JasparSpeciesAssignment {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tax_id: Option<String>,
    pub scientific_name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub common_name: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Optional JASPAR REST metadata enrichment for one local motif entry.
pub struct JasparRemoteMetadata {
    pub source_url: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub collection: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tax_group: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_class: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_family: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub data_type: Option<String>,
    #[serde(default)]
    pub species_assignments: Vec<JasparSpeciesAssignment>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One persisted JASPAR remote-metadata row keyed by local motif id.
pub struct JasparRemoteMetadataSnapshotRow {
    pub motif_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub motif_name: Option<String>,
    pub consensus_iupac: String,
    pub motif_length_bp: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_summary: Option<JasparCatalogRemoteSummary>,
    pub remote_metadata: JasparRemoteMetadata,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Persisted reusable JASPAR remote-metadata snapshot for catalog/expert
/// enrichment across sessions.
pub struct JasparRemoteMetadataSnapshot {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub filter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub limit: Option<usize>,
    #[serde(default)]
    pub requested_motifs: Vec<String>,
    pub registry_entry_count: usize,
    pub fetched_entry_count: usize,
    pub persisted_entry_count: usize,
    pub source: String,
    #[serde(default)]
    pub rows: Vec<JasparRemoteMetadataSnapshotRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One matrix column expanded for GUI/CLI expert inspection and simple logo
/// rendering.
pub struct JasparExpertColumn {
    pub position_1based: usize,
    pub total_count: f64,
    pub dominant_base: String,
    pub a_count: f64,
    pub c_count: f64,
    pub g_count: f64,
    pub t_count: f64,
    pub a_fraction: f64,
    pub c_fraction: f64,
    pub g_fraction: f64,
    pub t_fraction: f64,
    pub information_content_bits: f64,
    pub a_logo_bits: f64,
    pub c_logo_bits: f64,
    pub g_logo_bits: f64,
    pub t_logo_bits: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One histogram bin for a JASPAR score-family distribution.
pub struct JasparScoreDistributionBin {
    pub start_score: f64,
    pub end_score: f64,
    pub count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// One score-family panel inside the JASPAR entry expert view.
pub struct JasparScoreDistributionPanel {
    pub score_kind: TfbsScoreTrackValueKind,
    pub label: String,
    pub maximizing_sequence: String,
    pub maximizing_score: f64,
    pub maximizing_quantile: f64,
    pub minimizing_sequence: String,
    pub minimizing_score: f64,
    pub minimizing_quantile: f64,
    pub distribution: JasparScoreDistributionSummary,
    #[serde(default)]
    pub histogram_bins: Vec<JasparScoreDistributionBin>,
}

impl Default for JasparScoreDistributionPanel {
    fn default() -> Self {
        Self {
            score_kind: TfbsScoreTrackValueKind::LlrBits,
            label: String::new(),
            maximizing_sequence: String::new(),
            maximizing_score: 0.0,
            maximizing_quantile: 0.0,
            minimizing_sequence: String::new(),
            minimizing_score: 0.0,
            minimizing_quantile: 0.0,
            distribution: JasparScoreDistributionSummary::default(),
            histogram_bins: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Detailed expert-oriented view for one local JASPAR entry, including the
/// count matrix, a simple logo payload, multiple score families, and optional
/// remote species metadata.
pub struct JasparEntryExpertView {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub motif_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub motif_name: Option<String>,
    pub consensus_iupac: String,
    pub motif_length_bp: usize,
    pub registry_entry_count: usize,
    pub requested_token: String,
    pub random_sequence_length_bp: usize,
    pub random_seed: u64,
    pub background_model: String,
    pub include_remote_metadata: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_metadata: Option<JasparRemoteMetadata>,
    #[serde(default)]
    pub columns: Vec<JasparExpertColumn>,
    #[serde(default)]
    pub score_panels: Vec<JasparScoreDistributionPanel>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One visible feature-class lane summarized in a sequence-context view.
pub struct SequenceContextVisibleClass {
    pub class_id: String,
    #[serde(default)]
    pub feature_kinds: Vec<String>,
    pub matched_count: usize,
    pub returned_count: usize,
    #[serde(default)]
    pub prominent_labels: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One feature row surfaced in a compact sequence-context summary.
pub struct SequenceContextFeatureRow {
    pub feature_id: usize,
    pub kind: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub length_bp: usize,
    pub strand: String,
    pub label: String,
    #[serde(default)]
    pub labels: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_end_1based: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable summary of one DNA-sequence viewer context.
///
/// This report is intentionally chat-friendly and bundle-friendly: it exposes
/// the current viewport, the active visible feature classes, a compact feature
/// table, and short summary lines that ClawBio or other automation layers can
/// relay without having to infer biology from a raw SVG alone.
pub struct SequenceContextViewReport {
    pub schema: String,
    pub seq_id: String,
    pub sequence_length_bp: usize,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub mode: String,
    pub coordinate_mode: String,
    pub viewport_start_0based: usize,
    pub viewport_end_0based_exclusive: usize,
    pub viewport_span_bp: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_anchor: Option<SequenceGenomeAnchorSummary>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_view_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_view_end_1based: Option<usize>,
    #[serde(default)]
    pub visible_classes: Vec<SequenceContextVisibleClass>,
    pub matched_feature_count: usize,
    pub returned_feature_count: usize,
    pub limit: usize,
    #[serde(default)]
    pub rows: Vec<SequenceContextFeatureRow>,
    #[serde(default)]
    pub summary_lines: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One exported artifact within a sequence-context bundle, with deterministic
/// presentation metadata for chat/report consumers.
pub struct SequenceContextBundleArtifact {
    pub artifact_id: String,
    pub path: String,
    pub media_type: String,
    pub artifact_kind: String,
    pub caption: String,
    pub recommended_use: String,
    pub presentation_rank: usize,
    pub is_best_first_artifact: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable deterministic artifact manifest for one DNA-sequence view bundle.
///
/// This packages the current sequence-view SVG with the shared compact
/// sequence-context summary and an optional coordinate-bearing BED companion so
/// bundle-oriented callers can consume one operation result instead of
/// assembling multiple exports manually.
pub struct SequenceContextBundleExport {
    pub schema: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub output_dir: String,
    pub svg_path: String,
    pub summary_json_path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub summary_text_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub feature_bed_path: Option<String>,
    pub bundle_json_path: String,
    pub include_text_summary: bool,
    pub include_feature_bed: bool,
    pub include_restriction_sites: bool,
    #[serde(default)]
    pub restriction_enzymes: Vec<String>,
    #[serde(default)]
    pub artifacts: Vec<SequenceContextBundleArtifact>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub best_first_artifact_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub best_first_artifact_path: Option<String>,
    pub sequence_context_view: SequenceContextViewReport,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub feature_bed_export: Option<SequenceFeatureBedExportReport>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// How transcript-derived promoter windows are collapsed before annotation or
/// downstream reporting.
pub enum PromoterWindowCollapseMode {
    #[default]
    Transcript,
    Gene,
}

impl PromoterWindowCollapseMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Transcript => "transcript",
            Self::Gene => "gene",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Which allele should be materialized from one single-nucleotide variant.
pub enum VariantAlleleChoice {
    #[default]
    Reference,
    Alternate,
}

impl VariantAlleleChoice {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Reference => "reference",
            Self::Alternate => "alternate",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One promoter window derived from transcript TSS geometry.
pub struct PromoterWindowRecord {
    pub gene_label: Option<String>,
    pub gene_id: Option<String>,
    pub transcript_id: String,
    pub transcript_label: String,
    pub transcript_feature_id: Option<usize>,
    pub transcript_count: usize,
    pub transcript_ids: Vec<String>,
    pub transcript_labels: Vec<String>,
    pub strand: String,
    pub tss_local_0based: usize,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub source: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One grouped transcript-derived promoter interpretation row used for
/// alternative-promoter comparison in Promoter design.
pub struct AlternativePromoterComparisonRow {
    pub label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub representative_transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub representative_transcript_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub representative_transcript_feature_id: Option<usize>,
    pub transcript_count: usize,
    #[serde(default)]
    pub transcript_ids: Vec<String>,
    #[serde(default)]
    pub transcript_labels: Vec<String>,
    pub strand: String,
    pub representative_tss_local_0based: usize,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    #[serde(default)]
    pub source: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable grouped transcript-derived promoter comparison for one locus.
pub struct AlternativePromoterComparisonReport {
    pub schema: String,
    pub seq_id: String,
    pub sequence_length_bp: usize,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label_filter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_filter: Option<String>,
    pub promoter_upstream_bp: usize,
    pub promoter_downstream_bp: usize,
    pub transcript_window_count: usize,
    pub collapsed_window_count: usize,
    #[serde(default)]
    pub rows: Vec<AlternativePromoterComparisonRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One evidence item attached to a promoter candidate in the conservative
/// promoter-evidence ledger.
///
/// `kind` is intentionally string-valued so orthology, co-regulation,
/// anti-co-regulation, CUT&RUN, and future evidence layers can be added
/// without revving readers that only need to preserve unknown evidence rows.
pub struct PromoterEvidenceItem {
    pub evidence_id: String,
    pub kind: String,
    pub source: String,
    pub summary: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub end_0based_exclusive: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub overlap_bp: Option<usize>,
    #[serde(default)]
    pub metrics: BTreeMap<String, f64>,
    #[serde(default)]
    pub attributes: BTreeMap<String, String>,
    #[serde(default)]
    pub provenance_refs: Vec<String>,
    #[serde(default)]
    pub interpretation_tags: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One promoter candidate row with a transparent list of evidence items.
pub struct PromoterEvidenceMatrixRow {
    pub row_id: String,
    pub display_rank: usize,
    pub label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub representative_transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub representative_transcript_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub representative_transcript_feature_id: Option<usize>,
    pub transcript_count: usize,
    #[serde(default)]
    pub transcript_ids: Vec<String>,
    #[serde(default)]
    pub transcript_labels: Vec<String>,
    pub strand: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub representative_tss_local_0based: Option<usize>,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub source: String,
    pub evidence_count: usize,
    #[serde(default)]
    pub evidence_kind_counts: BTreeMap<String, usize>,
    #[serde(default)]
    pub evidence: Vec<PromoterEvidenceItem>,
    #[serde(default)]
    pub interpretation_tags: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Engine-owned promoter evidence matrix for one sequence/locus.
///
/// This report is an evidence ledger, not a final biological verdict. Rows are
/// deterministically ordered for review, and every future scoring/ranking layer
/// should keep its method and provenance explicit.
pub struct PromoterEvidenceMatrixReport {
    pub schema: String,
    pub seq_id: String,
    pub sequence_length_bp: usize,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label_filter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_filter: Option<String>,
    pub promoter_upstream_bp: usize,
    pub promoter_downstream_bp: usize,
    pub transcript_window_count: usize,
    pub promoter_candidate_count: usize,
    pub evidence_item_count: usize,
    #[serde(default)]
    pub evidence_kinds_observed: Vec<String>,
    pub ranking_mode: String,
    pub ranking_basis: String,
    #[serde(default)]
    pub rows: Vec<PromoterEvidenceMatrixRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Compact evidence signature used to compare promoter candidates across
/// transcript/isoform-derived promoter windows.
pub struct IsoformPromoterEvidenceSignature {
    pub signature_id: String,
    pub kind: String,
    pub source: String,
    pub label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub end_0based_exclusive: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
    pub summary: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One promoter group in an isoform-aware comparison for a single gene.
pub struct IsoformPromoterComparisonGroup {
    pub group_id: String,
    pub label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    pub strand: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub representative_tss_local_0based: Option<usize>,
    pub transcript_count: usize,
    #[serde(default)]
    pub transcript_ids: Vec<String>,
    #[serde(default)]
    pub transcript_labels: Vec<String>,
    #[serde(default)]
    pub evidence_kind_counts: BTreeMap<String, usize>,
    #[serde(default)]
    pub evidence_signatures: Vec<IsoformPromoterEvidenceSignature>,
    #[serde(default)]
    pub common_evidence_signatures: Vec<IsoformPromoterEvidenceSignature>,
    #[serde(default)]
    pub unique_evidence_signatures: Vec<IsoformPromoterEvidenceSignature>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One evidence signature that differs between isoform-derived promoter groups.
pub struct IsoformPromoterDifferentialEvidence {
    pub signature: IsoformPromoterEvidenceSignature,
    #[serde(default)]
    pub present_group_ids: Vec<String>,
    #[serde(default)]
    pub absent_group_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Engine-owned comparison of common and differential promoter-region evidence
/// between isoforms of one gene.
pub struct IsoformPromoterComparisonReport {
    pub schema: String,
    pub seq_id: String,
    pub sequence_length_bp: usize,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label_filter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_filter: Option<String>,
    pub promoter_upstream_bp: usize,
    pub promoter_downstream_bp: usize,
    pub transcript_window_count: usize,
    pub promoter_group_count: usize,
    pub comparison_evidence_item_count: usize,
    #[serde(default)]
    pub comparison_evidence_kinds_observed: Vec<String>,
    #[serde(default)]
    pub groups: Vec<IsoformPromoterComparisonGroup>,
    #[serde(default)]
    pub common_evidence_signatures: Vec<IsoformPromoterEvidenceSignature>,
    #[serde(default)]
    pub differential_evidence: Vec<IsoformPromoterDifferentialEvidence>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// User- or workflow-supplied expression evidence row that can be attached to
/// promoter candidates without making GENtle own a particular RNA-seq pipeline.
pub struct PromoterExpressionEvidenceInput {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    pub value: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub artifact_path: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One expression evidence row after it has been assigned to a promoter group.
pub struct PromoterExpressionEvidenceRecord {
    pub evidence_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    pub value: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    pub source: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub artifact_path: Option<String>,
    #[serde(default)]
    pub matched_by: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Expression evidence summarized for one promoter candidate.
pub struct PromoterExpressionEvidenceRow {
    pub promoter_group_id: String,
    pub label: String,
    #[serde(default)]
    pub transcript_ids: Vec<String>,
    #[serde(default)]
    pub transcript_labels: Vec<String>,
    pub expression_record_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mean_value: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_value: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    #[serde(default)]
    pub sample_ids: Vec<String>,
    #[serde(default)]
    pub conditions: Vec<String>,
    #[serde(default)]
    pub records: Vec<PromoterExpressionEvidenceRecord>,
    pub interpretation: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Engine-owned bridge between promoter candidates and expression-level
/// evidence supplied by workflows, RNA-read reports, or external artifacts.
pub struct PromoterExpressionEvidenceReport {
    pub schema: String,
    pub seq_id: String,
    pub sequence_length_bp: usize,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label_filter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_filter: Option<String>,
    pub promoter_upstream_bp: usize,
    pub promoter_downstream_bp: usize,
    pub promoter_group_count: usize,
    pub supplied_expression_record_count: usize,
    pub assigned_expression_record_count: usize,
    pub expression_source_label: String,
    #[serde(default)]
    pub rows: Vec<PromoterExpressionEvidenceRow>,
    #[serde(default)]
    pub unassigned_expression_records: Vec<PromoterExpressionEvidenceInput>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One component artifact in a promoter-design handoff index.
pub struct PromoterArtifactManifestEntry {
    pub artifact_id: String,
    pub artifact_kind: String,
    pub path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub schema_hint: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recommended_use: Option<String>,
    #[serde(default)]
    pub required: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One resolved component artifact row in a promoter-design handoff index.
pub struct PromoterArtifactManifestResolvedEntry {
    pub artifact_id: String,
    pub artifact_kind: String,
    pub path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub schema_hint: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recommended_use: Option<String>,
    pub required: bool,
    pub status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Lightweight index of promoter-design component artifacts. This is not a
/// narrative dossier; downstream tools can choose their own ordering/story.
pub struct PromoterArtifactManifestReport {
    pub schema: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label_filter: Option<String>,
    pub artifact_count: usize,
    pub present_artifact_count: usize,
    pub missing_required_artifact_count: usize,
    #[serde(default)]
    pub artifacts: Vec<PromoterArtifactManifestResolvedEntry>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One overlapping annotation/evidence row surfaced in a variant-promoter
/// context report.
pub struct VariantPromoterContextEvidenceRow {
    pub role: String,
    pub kind: String,
    pub label: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub strand: Option<String>,
    pub source: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable promoter-context summary for one variant on one sequence.
pub struct VariantPromoterContextReport {
    pub schema: String,
    pub seq_id: String,
    pub sequence_length_bp: usize,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub variant_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub variant_feature_id: Option<usize>,
    pub variant_start_0based: usize,
    pub variant_end_0based_exclusive: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub variant_class: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_ref: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_alt: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_anchor: Option<SequenceGenomeAnchorSummary>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chosen_gene_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chosen_transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chosen_transcript_label: Option<String>,
    pub transcript_ambiguity_status: String,
    pub promoter_upstream_bp: usize,
    pub promoter_downstream_bp: usize,
    pub promoter_overlap: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signed_tss_distance_bp: Option<isize>,
    #[serde(default)]
    pub overlapping_gene_labels: Vec<String>,
    #[serde(default)]
    pub overlapping_transcript_labels: Vec<String>,
    #[serde(default)]
    pub overlapping_promoter_labels: Vec<String>,
    #[serde(default)]
    pub overlapping_tfbs_labels: Vec<String>,
    #[serde(default)]
    pub overlapping_evidence: Vec<VariantPromoterContextEvidenceRow>,
    #[serde(default)]
    pub promoter_windows_considered: Vec<PromoterWindowRecord>,
    #[serde(default)]
    pub effect_tags: Vec<String>,
    #[serde(default)]
    pub suggested_assay_ids: Vec<String>,
    pub tfbs_focus_half_window_bp: usize,
    pub tfbs_near_variant_status: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tfbs_region_summary: Option<TfbsRegionSummary>,
    pub rationale: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One deterministic promoter-reporter fragment candidate derived from one
/// transcript/TSS and one variant.
pub struct PromoterReporterFragmentCandidate {
    pub candidate_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label: Option<String>,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub tss_local_0based: usize,
    pub variant_start_0based: usize,
    pub variant_end_0based_exclusive: usize,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub length_bp: usize,
    pub retain_downstream_from_tss_bp: usize,
    pub retain_upstream_beyond_variant_bp: usize,
    pub promoter_overlap: bool,
    pub signed_tss_distance_bp: isize,
    pub rank: usize,
    pub recommended: bool,
    pub rationale: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable candidate set for promoter-fragment suggestions that can be turned
/// into luciferase reporter inserts.
pub struct PromoterReporterCandidateSet {
    pub schema: String,
    pub seq_id: String,
    pub sequence_length_bp: usize,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub variant_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chosen_gene_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chosen_transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chosen_transcript_label: Option<String>,
    pub transcript_ambiguity_status: String,
    pub retain_downstream_from_tss_bp: usize,
    pub retain_upstream_beyond_variant_bp: usize,
    pub max_candidates: usize,
    pub recommended_candidate_id: String,
    #[serde(default)]
    pub suggested_assay_ids: Vec<String>,
    #[serde(default)]
    pub candidates: Vec<PromoterReporterFragmentCandidate>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequencingConfirmationStatus {
    Confirmed,
    Contradicted,
    #[default]
    InsufficientEvidence,
}

impl SequencingConfirmationStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Confirmed => "confirmed",
            Self::Contradicted => "contradicted",
            Self::InsufficientEvidence => "insufficient_evidence",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Orientation chosen for the best alignment of one confirmation read.
pub enum SequencingReadOrientation {
    #[default]
    Forward,
    ReverseComplement,
}

impl SequencingReadOrientation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::ReverseComplement => "reverse_complement",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Origin of one sequencing-confirmation evidence row.
pub enum SequencingConfirmationEvidenceKind {
    #[default]
    Sequence,
    Trace,
}

impl SequencingConfirmationEvidenceKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Sequence => "sequence",
            Self::Trace => "trace",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Biology-facing target classes supported by sequencing confirmation v1.
pub enum SequencingConfirmationTargetKind {
    #[default]
    FullSpan,
    Junction,
    FeaturePresence,
    ExpectedEdit,
    RestrictionSite,
}

impl SequencingConfirmationTargetKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FullSpan => "full_span",
            Self::Junction => "junction",
            Self::FeaturePresence => "feature_presence",
            Self::ExpectedEdit => "expected_edit",
            Self::RestrictionSite => "restriction_site",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
/// One requested construct-confirmation target on the expected sequence.
///
/// Coordinates are 0-based half-open on the expected construct. Junction
/// targets also set `junction_left_end_0based` to mark the break between the
/// two halves that must both be supported.
pub struct SequencingConfirmationTargetSpec {
    pub target_id: String,
    pub label: String,
    pub kind: SequencingConfirmationTargetKind,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub junction_left_end_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_bases: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub baseline_bases: Option<String>,
    pub required: bool,
}

impl Default for SequencingConfirmationTargetSpec {
    fn default() -> Self {
        Self {
            target_id: String::new(),
            label: String::new(),
            kind: SequencingConfirmationTargetKind::FullSpan,
            start_0based: 0,
            end_0based_exclusive: 0,
            junction_left_end_0based: None,
            expected_bases: None,
            baseline_bases: None,
            required: true,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Alignment-level discrepancy class recorded for one confirmation read.
pub enum SequencingConfirmationDiscrepancyKind {
    #[default]
    Mismatch,
    Insertion,
    Deletion,
}

impl SequencingConfirmationDiscrepancyKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Mismatch => "mismatch",
            Self::Insertion => "insertion",
            Self::Deletion => "deletion",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One mismatch/indel segment extracted from the best read alignment.
pub struct SequencingConfirmationDiscrepancy {
    pub kind: SequencingConfirmationDiscrepancyKind,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub target_start_0based: usize,
    pub target_end_0based_exclusive: usize,
    pub query_bases: String,
    pub target_bases: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequencingConfirmationVariantClassification {
    #[default]
    ExpectedMatch,
    IntendedEditConfirmed,
    ReferenceReversion,
    UnexpectedDifference,
    LowConfidenceOrAmbiguous,
    InsufficientEvidence,
}

impl SequencingConfirmationVariantClassification {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ExpectedMatch => "expected_match",
            Self::IntendedEditConfirmed => "intended_edit_confirmed",
            Self::ReferenceReversion => "reference_reversion",
            Self::UnexpectedDifference => "unexpected_difference",
            Self::LowConfidenceOrAmbiguous => "low_confidence_or_ambiguous",
            Self::InsufficientEvidence => "insufficient_evidence",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequencingConfirmationVariantRow {
    pub variant_id: String,
    pub label: String,
    pub target_id: Option<String>,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub expected_bases: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub baseline_bases: Option<String>,
    pub observed_bases: String,
    pub classification: SequencingConfirmationVariantClassification,
    pub status: SequencingConfirmationStatus,
    pub evidence_kind: SequencingConfirmationEvidenceKind,
    pub evidence_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub trace_id: Option<String>,
    pub read_seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub linked_seq_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence_min: Option<u8>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence_max: Option<u8>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence_mean: Option<f64>,
    pub confidence_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_center: Option<u32>,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Construct-level verdict for one requested confirmation target.
pub struct SequencingConfirmationTargetResult {
    pub target_id: String,
    pub label: String,
    pub kind: SequencingConfirmationTargetKind,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub junction_left_end_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_bases: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub baseline_bases: Option<String>,
    pub required: bool,
    pub status: SequencingConfirmationStatus,
    pub covered_bp: usize,
    pub target_length_bp: usize,
    pub support_read_ids: Vec<String>,
    pub contradicting_read_ids: Vec<String>,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Best-alignment and discrepancy summary for one input sequencing evidence row.
pub struct SequencingConfirmationReadResult {
    pub evidence_kind: SequencingConfirmationEvidenceKind,
    pub evidence_id: String,
    pub read_seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub trace_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub linked_seq_id: Option<String>,
    pub orientation: SequencingReadOrientation,
    pub usable: bool,
    pub best_alignment: SequenceAlignmentReport,
    pub covered_target_ids: Vec<String>,
    pub confirmed_target_ids: Vec<String>,
    pub contradicted_target_ids: Vec<String>,
    pub discrepancies: Vec<SequencingConfirmationDiscrepancy>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Stored construct-confirmation report shared by shell, GUI, and exports.
pub struct SequencingConfirmationReport {
    pub schema: String,
    pub report_id: String,
    pub expected_seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub baseline_seq_id: Option<String>,
    pub generated_at_unix_ms: u128,
    pub overall_status: SequencingConfirmationStatus,
    pub alignment_mode: PairwiseAlignmentMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub min_identity_fraction: f64,
    pub min_target_coverage_fraction: f64,
    pub allow_reverse_complement: bool,
    pub read_seq_ids: Vec<String>,
    #[serde(default)]
    pub trace_ids: Vec<String>,
    pub target_count: usize,
    pub reads: Vec<SequencingConfirmationReadResult>,
    pub targets: Vec<SequencingConfirmationTargetResult>,
    #[serde(default)]
    pub variants: Vec<SequencingConfirmationVariantRow>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Lightweight listing row for one sequencing-confirmation report.
pub struct SequencingConfirmationReportSummary {
    pub report_id: String,
    pub expected_seq_id: String,
    pub generated_at_unix_ms: u128,
    pub overall_status: SequencingConfirmationStatus,
    pub read_count: usize,
    pub target_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct FlexibilityBinScore {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct FlexibilityTrack {
    pub schema: String,
    pub track_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub model: FlexibilityModel,
    pub bin_bp: usize,
    pub smoothing_bp: Option<usize>,
    pub min_score: f64,
    pub max_score: f64,
    pub bins: Vec<FlexibilityBinScore>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlexibilityTrackSummary {
    pub track_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub model: FlexibilityModel,
    pub bin_bp: usize,
    pub smoothing_bp: Option<usize>,
    pub bin_count: usize,
    pub min_score: f64,
    pub max_score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderBandInfo {
    pub length_bp: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub relative_strength: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderInfo {
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub loading_hint: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_bp: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_bp: Option<usize>,
    pub band_count: usize,
    pub bands: Vec<DnaLadderBandInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderCatalog {
    pub schema: String,
    pub ladder_count: usize,
    pub ladders: Vec<DnaLadderInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderExportReport {
    pub path: String,
    pub ladder_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderBandInfo {
    pub length_nt: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub relative_strength: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderInfo {
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub loading_hint: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_nt: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_nt: Option<usize>,
    pub band_count: usize,
    pub bands: Vec<RnaLadderBandInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderCatalog {
    pub schema: String,
    pub ladder_count: usize,
    pub ladders: Vec<RnaLadderInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderExportReport {
    pub path: String,
    pub ladder_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteaseCatalogEntry {
    pub schema: String,
    pub name: String,
    pub cleavage_pattern: String,
    pub cut_offset: isize,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub aliases: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub category: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specificity: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cleavage_side: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub typical_applications: Vec<String>,
    #[serde(default)]
    pub sequence_specific: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteaseCatalogReport {
    pub schema: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub filter: Option<String>,
    pub total_protease_count: usize,
    pub returned_protease_count: usize,
    pub proteases: Vec<ProteaseCatalogEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteaseCatalogExportReport {
    pub path: String,
    pub returned_protease_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteaseDigestProteaseSummary {
    pub name: String,
    pub cleavage_pattern: String,
    pub cut_offset: isize,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub aliases: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteaseCleavageSite {
    pub protease_name: String,
    pub cleavage_boundary_0based: usize,
    pub cleavage_after_aa_1based: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upstream_residue: Option<char>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub downstream_residue: Option<char>,
    pub context: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteaseDigestPeptide {
    pub peptide_index: usize,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub source_start_aa_1based: usize,
    pub source_end_aa_1based: usize,
    pub length_aa: usize,
    pub sequence: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_left_cleavage_boundary_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_right_cleavage_boundary_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub created_seq_id: Option<SeqId>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteaseDigestReport {
    pub schema: String,
    pub source_seq_id: SeqId,
    pub source_length_aa: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_protein_derivation_mode: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_translation_table: Option<String>,
    pub requested_proteases: Vec<String>,
    pub resolved_proteases: Vec<ProteaseDigestProteaseSummary>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub missing_proteases: Vec<String>,
    pub min_length_aa: usize,
    pub materialized: bool,
    pub cleavage_site_count: usize,
    pub peptide_count: usize,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub created_seq_ids: Vec<SeqId>,
    pub sites: Vec<ProteaseCleavageSite>,
    pub peptides: Vec<ProteaseDigestPeptide>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeExtractionProvenance {
    pub seq_id: SeqId,
    pub recorded_at_unix_ms: u128,
    pub operation: String,
    pub genome_id: String,
    pub catalog_path: String,
    pub cache_dir: Option<String>,
    pub chromosome: Option<String>,
    pub start_1based: Option<usize>,
    pub end_1based: Option<usize>,
    pub gene_query: Option<String>,
    pub occurrence: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_extract_mode: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tss_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_upstream_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_downstream_bp: Option<usize>,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub strand: Option<char>,
    pub anchor_strand: Option<char>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchor_verified: Option<bool>,
    pub sequence_source_type: Option<String>,
    pub annotation_source_type: Option<String>,
    pub sequence_source: Option<String>,
    pub annotation_source: Option<String>,
    pub sequence_sha1: Option<String>,
    pub annotation_sha1: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceGenomeAnchorSummary {
    pub seq_id: String,
    pub genome_id: String,
    pub chromosome: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: Option<char>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchor_verified: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceAnchorPreparedGenomeOptionsSummary {
    pub seq_id: String,
    pub requested_genome_id: String,
    pub requested_catalog_key: String,
    pub requested_family: Option<String>,
    pub exact_prepared: bool,
    pub compatible_prepared_options: Vec<String>,
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
    RepeatFeatures,
    ArrayFeatures,
    ConstructReasoningOverlay,
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
    #[serde(default)]
    pub show_linear_sequence_panel: bool,
    #[serde(default = "DisplaySettings::default_sequence_panel_max_text_length_bp")]
    pub sequence_panel_max_text_length_bp: usize,
    pub auto_hide_sequence_panel_when_linear_bases_visible: bool,
    pub show_map_panel: bool,
    pub show_features: bool,
    pub show_cds_features: bool,
    pub show_gene_features: bool,
    pub show_mrna_features: bool,
    #[serde(default = "DisplaySettings::default_show_repeat_features")]
    pub show_repeat_features: bool,
    #[serde(default = "DisplaySettings::default_show_array_features")]
    pub show_array_features: bool,
    #[serde(default = "DisplaySettings::default_show_construct_reasoning_overlay")]
    pub show_construct_reasoning_overlay: bool,
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
    pub const fn default_show_repeat_features() -> bool {
        true
    }

    pub const fn default_show_array_features() -> bool {
        true
    }

    pub const fn default_show_construct_reasoning_overlay() -> bool {
        true
    }

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
            show_linear_sequence_panel: false,
            sequence_panel_max_text_length_bp: Self::default_sequence_panel_max_text_length_bp(),
            auto_hide_sequence_panel_when_linear_bases_visible: false,
            show_map_panel: true,
            show_features: true,
            show_cds_features: true,
            show_gene_features: true,
            show_mrna_features: true,
            show_repeat_features: Self::default_show_repeat_features(),
            show_array_features: Self::default_show_array_features(),
            show_construct_reasoning_overlay: Self::default_show_construct_reasoning_overlay(),
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
    pub protein_derivation_report: Option<ProteinDerivationReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reverse_translation_report: Option<ReverseTranslationReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub protease_digest_report: Option<ProteaseDigestReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub protein_residue_genomic_coordinates: Option<ProteinResidueGenomicCoordinateReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exon_skip_selection_plan: Option<ExonSkipSelectionPlan>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exon_skip_materialization: Option<ExonSkipMaterializationReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cdna_assay_test_report: Option<Box<CdnaAssayTestReport>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cdna_assay_product_materialization: Option<Box<CdnaAssayProductMaterialization>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_qpcr_panel: Option<Box<TranscriptQpcrPanelReport>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub primer_specificity_report: Option<Box<PrimerSpecificityReport>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub construct_reasoning_graph: Option<Box<ConstructReasoningGraph>>,
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
    pub cutrun_dataset_list: Option<CutRunDatasetListReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cutrun_dataset_status: Option<CutRunDatasetStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cutrun_dataset_projection: Option<CutRunDatasetProjectionReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cutrun_read_report: Option<CutRunReadReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cutrun_read_report_summaries: Option<Vec<CutRunReadReportSummary>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cutrun_read_coverage_export: Option<CutRunReadCoverageExport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cutrun_regulatory_support: Option<CutRunRegulatorySupportReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_set_resolution: Option<GeneSetResolutionReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_set_promoter_cohort: Option<GeneSetPromoterCohortReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_set_cutrun_regulatory_support: Option<GeneSetCutRunRegulatorySupportReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ortholog_promoter_cohort: Option<OrthologPromoterCohortReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ortholog_promoter_comparison: Option<OrthologPromoterComparisonReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub read_acquisition_report: Option<ReadAcquisitionReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub microarray_projection: Option<MicroarrayProjectionReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub probe_region_evidence_interpretation: Option<ProbeRegionEvidenceInterpretationReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_coordinate_projection: Option<GenomeCoordinateProjectionReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rna_read_gene_support_summary: Option<RnaReadGeneSupportSummary>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rna_read_gene_support_audit: Option<RnaReadGeneSupportAudit>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rna_read_target_quality_export: Option<RnaReadTargetQualityExport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rna_read_batch_map_report: Option<RnaReadBatchMapReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rna_read_isoform_preflight: Option<RnaReadIsoformPreflightReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tfbs_region_summary: Option<TfbsRegionSummary>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tfbs_score_tracks: Option<TfbsScoreTrackReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tfbs_track_similarity: Option<TfbsTrackSimilarityReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub multi_gene_promoter_tfbs: Option<MultiGenePromoterTfbsReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_cohort_comparison: Option<PromoterCohortComparisonReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub repeat_annotation_query: Option<RepeatAnnotationQueryReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequence_repeat_overlaps: Option<SequenceRepeatOverlapReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub repeat_feature_materialization: Option<RepeatFeatureMaterializationReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub repeat_environment_cohort: Option<RepeatEnvironmentCohortReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub window_cohort_tfbs: Option<WindowCohortTfbsReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tfbs_hit_scan: Option<TfbsHitScanReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub restriction_site_scan: Option<RestrictionSiteScanReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub jaspar_remote_metadata_snapshot: Option<JasparRemoteMetadataSnapshot>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub jaspar_catalog_report: Option<JasparCatalogReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_query_resolution_report: Option<TfQueryResolutionReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub jaspar_entry_expert_view: Option<JasparEntryExpertView>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub jaspar_registry_benchmark: Option<JasparRegistryBenchmarkReport>,
    pub jaspar_entry_presentation: Option<JasparEntryPresentationReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequence_context_view: Option<SequenceContextViewReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequence_context_bundle: Option<SequenceContextBundleExport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub variant_promoter_context: Option<VariantPromoterContextReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub alternative_promoter_comparison: Option<AlternativePromoterComparisonReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_evidence_matrix: Option<PromoterEvidenceMatrixReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub isoform_promoter_comparison: Option<IsoformPromoterComparisonReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_expression_evidence: Option<PromoterExpressionEvidenceReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_artifact_manifest: Option<PromoterArtifactManifestReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_reporter_candidates: Option<PromoterReporterCandidateSet>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reporter_catalog: Option<ReporterCatalogReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reporter_recommendation: Option<ReporterRecommendationResult>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reporter_corpus_export: Option<ReporterCorpusExport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reporter_construct_handoff: Option<ReporterConstructHandoffPlan>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uniprot_projection_audit: Option<Box<UniprotProjectionAuditReport>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uniprot_projection_audit_parity: Option<Box<UniprotProjectionAuditParityReport>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lab_assistant_instructions: Option<Box<LabAssistantInstructionsExport>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
/// Compact manifest for genome-anchored microarray projection tracks.
///
/// The manifest is intentionally small: it records the coordinate contract and
/// points at per-contrast TSVs so projection can stream only intervals
/// overlapping the current sequence anchor.
pub struct MicroarrayTrackManifest {
    pub schema: String,
    #[serde(alias = "dataset_id")]
    pub dataset: String,
    pub platform: String,
    #[serde(alias = "normalization_method")]
    pub normalization: String,
    pub coordinate_system: String,
    pub supported_genome_ids: Vec<String>,
    pub contrast_order: Vec<String>,
    #[serde(alias = "tracks")]
    pub contrasts: Vec<MicroarrayTrackContrast>,
    #[serde(default, alias = "projection_maps")]
    pub coordinate_projections: Vec<GenomeCoordinateProjectionSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_path: Option<String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// One per-contrast projected TSV entry in a microarray track manifest.
pub struct MicroarrayTrackContrast {
    #[serde(alias = "name")]
    pub contrast: String,
    pub level: String,
    #[serde(alias = "tsv_path")]
    pub path: String,
    pub row_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Request contract for chromosome-ordered Affymetrix probe/probeset planning.
pub struct ProbeRegionRequest {
    pub cel_paths: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dataset: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata_path: Option<String>,
    pub genes: Vec<String>,
    pub loci: Vec<String>,
    pub transcript_cluster_ids: Vec<String>,
    pub probeset_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub platform: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_library_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub condition_column: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_column: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub block_column: Option<String>,
    pub paired_by_replicate_suffix: bool,
    pub plot: bool,
    pub normalization: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub output_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cache_dir: Option<String>,
    pub dry_run: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Selector summary separated from the original request for GUI-friendly display.
pub struct ProbeRegionSelectorPlan {
    pub genes: Vec<String>,
    pub loci: Vec<String>,
    pub transcript_cluster_ids: Vec<String>,
    pub probeset_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Filesystem status for one CEL, metadata, annotation, output, or cache path.
pub struct ProbeRegionFileStatus {
    pub path: String,
    pub role: String,
    pub exists: bool,
    pub is_file: bool,
    pub is_dir: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size_bytes: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub modified_unix_seconds: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Resolved annotation/library source for a planned Affymetrix probe-region run.
pub struct ProbeRegionAnnotationSourcePlan {
    pub path: Option<ProbeRegionFileStatus>,
    pub source_kind: String,
    pub usable: bool,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub vendor_support_files: Vec<ProbeRegionFileStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub required_r_package: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Platform/backend plan derived from user input and known Affymetrix mappings.
pub struct ProbeRegionPlatformPlan {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub requested: Option<String>,
    pub normalized: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub registry_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub family: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub species: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub genome_builds: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub backend_kinds: Vec<String>,
    pub backend_hint: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bioconductor_package: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cdf_package: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_package: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub staging_directory: Option<String>,
    pub confidence: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Dependency check row for non-downloading preflight of local analysis tools.
pub struct ProbeRegionDependencyCheck {
    pub name: String,
    pub kind: String,
    pub required: bool,
    pub status: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Per-condition sample count parsed from optional CEL sample metadata.
pub struct ProbeRegionConditionSummary {
    pub condition: String,
    pub sample_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Metadata table preview used for sample grouping and default contrast planning.
pub struct ProbeRegionMetadataPlan {
    pub status: String,
    pub delimiter: String,
    pub columns: Vec<String>,
    pub row_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_column: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub condition_column: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub block_column: Option<String>,
    pub sample_count: usize,
    pub conditions: Vec<ProbeRegionConditionSummary>,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Planned default comparison between two metadata-defined conditions.
pub struct ProbeRegionContrastPlan {
    pub contrast: String,
    pub numerator_condition: String,
    pub denominator_condition: String,
    pub numerator_sample_count: usize,
    pub denominator_sample_count: usize,
    pub status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Candidate execution backend and the local inputs it still needs.
pub struct ProbeRegionBackendCandidate {
    pub backend: String,
    pub status: String,
    pub required_inputs: Vec<String>,
    pub missing: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub helper_script: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub suggested_command: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Stage-one preflight plan for a future Affymetrix probe-region execution.
pub struct ProbeRegionPlan {
    pub schema: String,
    pub stage: String,
    pub implementation_status: String,
    pub input_mode: String,
    pub request: ProbeRegionRequest,
    pub selectors: ProbeRegionSelectorPlan,
    pub cel_files: Vec<ProbeRegionFileStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata: Option<ProbeRegionFileStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata_plan: Option<ProbeRegionMetadataPlan>,
    pub annotation_source: ProbeRegionAnnotationSourcePlan,
    pub platform: ProbeRegionPlatformPlan,
    pub dependencies: Vec<ProbeRegionDependencyCheck>,
    pub backend_candidates: Vec<ProbeRegionBackendCandidate>,
    pub contrasts: Vec<ProbeRegionContrastPlan>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub output_dir_status: Option<ProbeRegionFileStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cache_dir_status: Option<ProbeRegionFileStatus>,
    pub planned_outputs: Vec<String>,
    pub cache_compatibility_keys: Vec<String>,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
    pub preflight_ok: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Explicit external-backend execution report for a persisted probe-region plan.
pub struct ProbeRegionBackendRunReport {
    pub schema: String,
    pub plan_path: String,
    pub backend: String,
    pub command: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub output_dir: Option<String>,
    pub allow_external_execution: bool,
    pub executed: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exit_code: Option<i32>,
    pub stdout: String,
    pub stderr: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inspection: Option<ProbeRegionOutputInspection>,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// One chromosome-ordered preview row from a completed probe-region output table.
pub struct ProbeRegionOutputPreviewRow {
    pub chromosome: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stop_1based: Option<usize>,
    pub probeset_or_region_id: String,
    pub transcript_cluster_id: String,
    pub gene_symbol: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Read-only inspection of a completed `probe_regions_oligo.R` output folder.
pub struct ProbeRegionOutputInspection {
    pub schema: String,
    pub output_dir: String,
    pub usable: bool,
    pub region_table: ProbeRegionFileStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub probe_table: Option<ProbeRegionFileStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_table: Option<ProbeRegionFileStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub normalized_matrix_manifest: Option<ProbeRegionFileStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provenance: Option<ProbeRegionFileStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub backend: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub platform: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub platform_package: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub normalization: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordinate_system: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_build: Option<String>,
    #[serde(default, alias = "projection_maps")]
    pub coordinate_projections: Vec<GenomeCoordinateProjectionSpec>,
    pub projection_ready: bool,
    pub projection_blockers: Vec<String>,
    pub target_levels: Vec<String>,
    pub artifact_paths: Vec<String>,
    pub row_count: usize,
    pub column_count: usize,
    pub probe_row_count: usize,
    pub probe_parent_feature_count: usize,
    pub feature_count: usize,
    pub transcript_cluster_count: usize,
    pub chromosome_count: usize,
    pub chromosomes: Vec<String>,
    pub gene_symbols: Vec<String>,
    pub sample_columns: Vec<String>,
    pub condition_summary_columns: Vec<String>,
    pub logfc_columns: Vec<String>,
    pub preview_rows: Vec<ProbeRegionOutputPreviewRow>,
    pub required_columns_missing: Vec<String>,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Deterministic SVG export summary for inspected probe-region output.
pub struct ProbeRegionOutputSvgExport {
    pub schema: String,
    pub output_dir: String,
    pub svg_path: String,
    pub row_count: usize,
    pub intensity_track_count: usize,
    pub logfc_track_count: usize,
    pub chromosome_count: usize,
    pub projection_ready: bool,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Deterministic SVG export summary for a probe-region evidence interpretation report.
pub struct ProbeRegionEvidenceSvgExport {
    pub schema: String,
    pub report_path: String,
    pub svg_path: String,
    pub evidence_row_count: usize,
    pub transcript_count: usize,
    pub parent_feature_count: usize,
    pub junction_span_count: usize,
    pub ambiguity_tags: Vec<String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Deterministic import report for explicit APT output plus annotation tables.
pub struct ProbeRegionAptImportReport {
    pub schema: String,
    pub apt_summary_path: String,
    pub annotation_path: String,
    pub output_dir: String,
    pub platform: String,
    pub normalization: String,
    pub coordinate_system: String,
    pub genome_build: String,
    pub summary_row_count: usize,
    pub annotation_row_count: usize,
    pub written_row_count: usize,
    pub missing_annotation_count: usize,
    pub skipped_invalid_count: usize,
    pub probe_row_count: usize,
    pub sample_columns: Vec<String>,
    pub probe_intensity_path: Option<String>,
    pub probe_intensity_source: Option<String>,
    pub probe_intensity_sample_columns: Vec<String>,
    pub missing_probe_intensity_count: usize,
    pub metadata_path: Option<String>,
    pub condition_column: Option<String>,
    pub sample_column: Option<String>,
    pub condition_columns: Vec<String>,
    pub logfc_columns: Vec<String>,
    pub warnings: Vec<String>,
    pub inspection: ProbeRegionOutputInspection,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
/// Ambiguity-preserving interpretation of projected probe-region evidence.
///
/// This report deliberately stops short of isoform support claims. It compares
/// already-projected array features with transcript/exon geometry and records
/// compatibility, constraints, and unresolved ambiguity for downstream review.
pub struct ProbeRegionEvidenceInterpretationReport {
    pub schema: String,
    pub seq_id: SeqId,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label: Option<String>,
    pub level: String,
    pub coordinate_frame: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordinate_system: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordinate_chromosome: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_abs_logfc: Option<f64>,
    pub array_feature_count: usize,
    pub transcript_count: usize,
    pub evidence_rows: Vec<ProbeRegionEvidenceMappingRow>,
    pub transcript_rows: Vec<ProbeRegionEvidenceTranscriptRow>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
/// One projected array feature interpreted against transcript/exon geometry.
pub struct ProbeRegionEvidenceMappingRow {
    pub evidence_id: String,
    pub level: String,
    pub feature_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parent_feature_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub intensity_source: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub end_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub logfc: Option<f64>,
    pub overlapping_transcript_ids: Vec<String>,
    pub overlapping_exon_count: usize,
    pub transcript_mappings: Vec<ProbeRegionEvidenceTranscriptMapping>,
    pub mapping_status: String,
    pub ambiguity_tags: Vec<String>,
    pub relationship: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
/// Coordinate-normalized geometry mapping for one projected array evidence row.
pub struct ProbeRegionEvidenceTranscriptMapping {
    pub transcript_id: String,
    pub coordinate_frame: String,
    pub mapping_kind: String,
    pub geometry_score: f64,
    pub geometry_score_class: String,
    pub score_basis: Vec<String>,
    pub exon_ordinals: Vec<usize>,
    pub exon_ranges_1based: Vec<String>,
    pub local_exon_ranges_1based: Vec<String>,
    pub junction_spans: Vec<ProbeRegionEvidenceJunctionSpan>,
    pub overlap_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// One exon-exon boundary spanned by a projected array evidence interval.
///
/// The `genomic_*` fields are the report coordinate-frame values used for
/// rendering. When the report frame is local, the optional `local_*` fields
/// preserve the same source coordinates explicitly.
pub struct ProbeRegionEvidenceJunctionSpan {
    pub from_exon_ordinal: usize,
    pub to_exon_ordinal: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub local_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub local_end_1based: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
/// Per-transcript view of projected array evidence geometry.
pub struct ProbeRegionEvidenceTranscriptRow {
    pub transcript_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
    pub exon_count: usize,
    pub compatible_evidence_count: usize,
    pub constraining_evidence_count: usize,
    pub shared_evidence_count: usize,
    pub unique_evidence_count: usize,
    pub unmapped_evidence_count: usize,
    pub compatible_geometry_score: f64,
    pub shared_geometry_score: f64,
    pub unique_geometry_score: f64,
    pub constraining_geometry_score: f64,
    pub review_status: String,
    pub relationship_summary: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Describes an explicit coordinate-projection map between two genome builds.
pub struct GenomeCoordinateProjectionSpec {
    pub source_genome_id: String,
    pub target_genome_id: String,
    pub method: String,
    pub path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// One projected interval emitted by the shared genome-coordinate projection path.
pub struct GenomeCoordinateProjectionInterval {
    pub source_chrom: String,
    pub source_start_1based: usize,
    pub source_end_1based: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_strand: Option<String>,
    pub target_chrom: String,
    pub target_start_1based: usize,
    pub target_end_1based: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_strand: Option<String>,
    pub method: String,
    pub status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Structured result emitted by genome-coordinate interval projection.
pub struct GenomeCoordinateProjectionReport {
    pub schema: String,
    pub source_genome_id: String,
    pub target_genome_id: String,
    pub projection_path: String,
    pub method: String,
    pub input_chrom: String,
    pub input_start_1based: usize,
    pub input_end_1based: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub input_strand: Option<String>,
    pub mapped: bool,
    #[serde(default)]
    pub intervals: Vec<GenomeCoordinateProjectionInterval>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
/// Structured result emitted by `ProjectMicroarrayTrack`.
pub struct MicroarrayProjectionReport {
    pub schema: String,
    pub seq_id: SeqId,
    pub manifest_path: String,
    pub dataset: String,
    pub platform: String,
    pub normalization: String,
    pub coordinate_system: String,
    pub coordinate_projection_used: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordinate_projection_method: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordinate_projection_path: Option<String>,
    pub anchor_genome_id: String,
    pub anchor_chromosome: String,
    pub anchor_start_1based: usize,
    pub anchor_end_1based: usize,
    pub anchor_strand: String,
    pub requested_contrasts: Vec<String>,
    pub projected_contrasts: Vec<String>,
    pub level: String,
    pub parsed_rows: usize,
    pub imported_features: usize,
    pub skipped_rows: usize,
    pub skipped_invalid: usize,
    pub skipped_wrong_chromosome: usize,
    pub skipped_non_overlap: usize,
    pub skipped_filter: usize,
    pub skipped_projection_unmapped: usize,
    pub truncated_at_limit: bool,
    pub warnings: Vec<String>,
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
/// Progress snapshot for shared primer-pair / qPCR design operations.
pub struct PrimerDesignProgress {
    pub seq_id: String,
    pub design_kind: String,
    pub backend_requested: String,
    pub backend_used: String,
    pub stage: String,
    pub detail: String,
    pub roi_start_0based: usize,
    pub roi_end_0based_exclusive: usize,
    pub forward_candidate_count: Option<usize>,
    pub reverse_candidate_count: Option<usize>,
    pub probe_candidate_count: Option<usize>,
    pub pair_candidate_combinations: Option<usize>,
    pub pair_evaluated: Option<usize>,
    pub pair_evaluation_limit: Option<usize>,
    pub pair_evaluation_limited: Option<bool>,
    pub accepted_pair_count: Option<usize>,
    pub assay_candidate_combinations: Option<usize>,
    pub assays_evaluated: Option<usize>,
    pub accepted_assay_count: Option<usize>,
    pub max_output: usize,
    pub done: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Union of long-running operation progress events.
pub enum OperationProgress {
    PrimerDesign(PrimerDesignProgress),
    Tfbs(TfbsProgress),
    GenomePrepare(PrepareGenomeProgress),
    GenomeTrackImport(GenomeTrackImportProgress),
    DbSnpFetch(DbSnpFetchProgress),
    ReadAcquisition(SharedAssetActivityStatus),
    RnaReadInterpret(RnaReadInterpretProgress),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Immutable operation journal row.
pub struct OperationRecord {
    pub run_id: RunId,
    pub op: Operation,
    pub result: OpResult,
}

/// Compact one-row summary for a history transition exposed to adapters.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EngineHistoryTransitionSummary {
    pub op_id: OpId,
    pub run_id: RunId,
    pub operation: String,
}

/// Session-local undo/redo availability summary.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EngineHistorySummary {
    pub schema: String,
    pub undo_count: usize,
    pub redo_count: usize,
    pub history_limit: usize,
    pub operation_log_count: usize,
    pub next_undo: Option<EngineHistoryTransitionSummary>,
    pub next_redo: Option<EngineHistoryTransitionSummary>,
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
    pub declared_contents_exclusive: bool,
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

#[derive(Debug, Clone, Serialize, Deserialize)]
/// BLAST hit reduced to the fields needed for feature import/overlay pipelines.
pub struct BlastHitFeatureInput {
    pub subject_id: String,
    pub query_start_1based: usize,
    pub query_end_1based: usize,
    pub subject_start_1based: usize,
    pub subject_end_1based: usize,
    pub identity_percent: f64,
    pub bit_score: f64,
    pub evalue: f64,
    pub query_coverage_percent: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Provenance bundle describing exactly how a BLAST search was invoked.
///
/// This is the record to inspect when reproducing a search or explaining which
/// executable/options/catalog paths produced a report.
pub struct BlastInvocationProvenance {
    pub genome_id: String,
    pub query_label: String,
    pub query_length: usize,
    pub max_hits: usize,
    pub task: String,
    pub blastn_executable: String,
    pub blast_db_prefix: String,
    pub command: Vec<String>,
    pub command_line: String,
    pub catalog_path: Option<String>,
    pub cache_dir: Option<String>,
    #[serde(default)]
    pub options_override_json: Option<serde_json::Value>,
    #[serde(default)]
    pub effective_options_json: Option<serde_json::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
/// Optional threshold overrides for genome BLAST post-filtering.
pub struct BlastThresholdOptions {
    pub max_evalue: Option<f64>,
    pub min_identity_percent: Option<f64>,
    pub min_query_coverage_percent: Option<f64>,
    pub min_alignment_length_bp: Option<usize>,
    pub min_bit_score: Option<f64>,
    pub unique_best_hit: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
/// Request-time BLAST execution options before defaults/project overrides merge.
pub struct BlastRunOptions {
    pub task: Option<String>,
    pub max_hits: Option<usize>,
    #[serde(default)]
    pub thresholds: BlastThresholdOptions,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Resolved BLAST options after defaults/project/request layering.
pub struct BlastResolvedOptions {
    pub task: String,
    pub max_hits: usize,
    #[serde(default)]
    pub thresholds: BlastThresholdOptions,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
/// Prepared-genome fallback policy when an anchor references a non-prepared id.
pub enum GenomeAnchorPreparedFallbackPolicy {
    Off,
    SingleCompatible,
    AlwaysExplicit,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GenomeTrackSyncReport {
    pub subscriptions_considered: usize,
    pub target_sequences: usize,
    pub applied_imports: usize,
    pub failed_imports: usize,
    pub warnings_count: usize,
    pub errors: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LigationProtocol {
    Sticky,
    Blunt,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ExportFormat {
    GenBank,
    Fasta,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PrimerLibraryMode {
    Enumerate,
    Sample,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PcrPrimerSpec {
    pub sequence: String,
    pub anneal_len: Option<usize>,
    pub max_mismatches: Option<usize>,
    pub require_3prime_exact_bases: Option<usize>,
    pub library_mode: Option<PrimerLibraryMode>,
    pub max_variants: Option<usize>,
    pub sample_seed: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpMutationSpec {
    pub zero_based_position: usize,
    pub reference: String,
    pub alternate: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct PrimerDesignBaseLock {
    pub offset_0based: usize,
    pub base: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PrimerDesignPairConstraint {
    pub require_roi_flanking: bool,
    #[serde(default)]
    pub required_amplicon_motifs: Vec<String>,
    #[serde(default)]
    pub forbidden_amplicon_motifs: Vec<String>,
    pub fixed_amplicon_start_0based: Option<usize>,
    pub fixed_amplicon_end_0based_exclusive: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct PrimerDesignSideConstraint {
    pub min_length: usize,
    pub max_length: usize,
    pub location_0based: Option<usize>,
    pub start_0based: Option<usize>,
    pub end_0based: Option<usize>,
    pub min_tm_c: f64,
    pub max_tm_c: f64,
    pub min_gc_fraction: f64,
    pub max_gc_fraction: f64,
    pub max_anneal_hits: usize,
    pub non_annealing_5prime_tail: Option<String>,
    pub fixed_5prime: Option<String>,
    pub fixed_3prime: Option<String>,
    #[serde(default)]
    pub required_motifs: Vec<String>,
    #[serde(default)]
    pub forbidden_motifs: Vec<String>,
    #[serde(default)]
    pub locked_positions: Vec<PrimerDesignBaseLock>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignPrimerRecord {
    pub sequence: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub length_bp: usize,
    pub anneal_length_bp: usize,
    pub non_annealing_5prime_tail_bp: usize,
    pub tm_c: f64,
    pub gc_fraction: f64,
    pub anneal_hits: usize,
    pub three_prime_base: String,
    pub three_prime_gc_clamp: bool,
    pub longest_homopolymer_run_bp: usize,
    pub self_complementary_run_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignPairRuleFlags {
    pub roi_covered: bool,
    pub amplicon_size_in_range: bool,
    pub tm_delta_in_range: bool,
    pub forward_secondary_structure_ok: bool,
    pub reverse_secondary_structure_ok: bool,
    pub primer_pair_dimer_risk_low: bool,
    pub forward_three_prime_gc_clamp: bool,
    pub reverse_three_prime_gc_clamp: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignPairRecord {
    pub rank: usize,
    pub score: f64,
    pub forward: PrimerDesignPrimerRecord,
    pub reverse: PrimerDesignPrimerRecord,
    pub amplicon_start_0based: usize,
    pub amplicon_end_0based_exclusive: usize,
    pub amplicon_length_bp: usize,
    pub tm_delta_c: f64,
    pub primer_pair_complementary_run_bp: usize,
    pub primer_pair_3prime_complementary_run_bp: usize,
    pub rule_flags: PrimerDesignPairRuleFlags,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Explicit exon-junction placement preference for transcript-aware
/// PCR/qPCR designs.
pub enum PrimerExonJunctionPolicy {
    #[default]
    NoPreference,
    MustSpan,
    MustNotSpan,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Explicit intron-separation preference for transcript-aware PCR/qPCR
/// designs.
pub enum PrimerIntronSeparationPolicy {
    #[default]
    NoPreference,
    MustSeparateByIntron,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Primer role in a specificity-confirmation run.
#[derive(Default)]
pub enum PrimerSpecificityPrimerRole {
    #[default]
    Forward,
    Reverse,
}

impl PrimerSpecificityPrimerRole {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Pairing class inspected for Primer-BLAST-style warnings.
#[derive(Default)]
pub enum PrimerSpecificityAmpliconKind {
    #[default]
    ForwardReverse,
    ForwardForward,
    ReverseReverse,
}

impl PrimerSpecificityAmpliconKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ForwardReverse => "forward_reverse",
            Self::ForwardForward => "forward_forward",
            Self::ReverseReverse => "reverse_reverse",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Primer sequence bundle used for local specificity confirmation. BLAST uses
/// only the annealing sequence; 5' tails are kept as provenance.
pub struct PrimerSpecificityInputPrimer {
    pub role: PrimerSpecificityPrimerRole,
    pub full_sequence: String,
    pub annealing_sequence: String,
    pub annealing_length_bp: usize,
    pub non_annealing_5prime_tail: String,
    pub non_annealing_5prime_tail_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One primer BLAST hit normalized into genomic coordinates and primer-end
/// compatibility metrics.
pub struct PrimerSpecificityPrimerHit {
    pub hit_index: usize,
    pub role: PrimerSpecificityPrimerRole,
    pub subject_id: String,
    pub identity_percent: f64,
    pub alignment_length_bp: usize,
    pub mismatches: usize,
    pub gap_opens: usize,
    pub query_start_1based: usize,
    pub query_end_1based: usize,
    pub subject_start_1based: usize,
    pub subject_end_1based: usize,
    pub subject_min_1based: usize,
    pub subject_max_1based: usize,
    pub strand: String,
    pub evalue: f64,
    pub bit_score: f64,
    pub query_coverage_fraction: f64,
    pub three_prime_window_bp: usize,
    pub three_prime_mismatches: usize,
    pub accepted_by_policy: bool,
    #[serde(default)]
    pub rejection_reasons: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Candidate product implied by compatible primer-hit orientations.
pub struct PrimerSpecificityAmplicon {
    pub kind: PrimerSpecificityAmpliconKind,
    pub subject_id: String,
    pub left_role: PrimerSpecificityPrimerRole,
    pub left_hit_index: usize,
    pub right_role: PrimerSpecificityPrimerRole,
    pub right_hit_index: usize,
    pub start_1based: usize,
    pub end_1based: usize,
    pub length_bp: usize,
    pub combined_mismatches: usize,
    pub max_three_prime_mismatches: usize,
    pub terminal_policy_pass: bool,
    pub intended: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub intended_reason: Option<String>,
    pub specificity_failure: bool,
    #[serde(default)]
    pub failure_reasons: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Summary badge for a local primer specificity report.
pub struct PrimerSpecificitySummary {
    pub specificity_pass: bool,
    pub status: String,
    pub primer_hit_count: usize,
    pub accepted_primer_hit_count: usize,
    pub amplicon_count: usize,
    pub intended_amplicon_count: usize,
    pub unintended_amplicon_count: usize,
    pub failing_unintended_amplicon_count: usize,
    pub summary: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Engine-owned local Primer-BLAST-parity report for one primer pair against a
/// prepared genome BLAST database.
pub struct PrimerSpecificityReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub primer_report_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pair_rank: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pair_index: Option<usize>,
    pub target_kind: String,
    pub target_genome_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub catalog_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cache_dir: Option<String>,
    #[serde(default)]
    pub policy: PrimerSpecificityPolicy,
    #[serde(default)]
    pub primers: Vec<PrimerSpecificityInputPrimer>,
    #[serde(default)]
    pub blast_preflight: BlastExternalBinaryPreflightReport,
    #[serde(default)]
    pub blast_runs: Vec<BlastInvocationProvenance>,
    #[serde(default)]
    pub forward_hits: Vec<PrimerSpecificityPrimerHit>,
    #[serde(default)]
    pub reverse_hits: Vec<PrimerSpecificityPrimerHit>,
    #[serde(default)]
    pub amplicons: Vec<PrimerSpecificityAmplicon>,
    #[serde(default)]
    pub summary: PrimerSpecificitySummary,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Deterministic exact-run QC for one assay oligo.
pub struct OligoQcOligoRecord {
    pub label: String,
    pub role: String,
    pub sequence: String,
    pub length_bp: usize,
    pub gc_fraction: f64,
    pub three_prime_base: String,
    pub three_prime_gc_clamp: bool,
    pub longest_homopolymer_run_bp: usize,
    pub self_complementary_run_bp: usize,
    pub self_3prime_complementary_run_bp: usize,
    pub status: String,
    pub summary: String,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Deterministic exact-run QC for a pair of assay oligos.
pub struct OligoQcInteractionRecord {
    pub left_label: String,
    pub right_label: String,
    pub left_role: String,
    pub right_role: String,
    pub max_complementary_run_bp: usize,
    pub left_3prime_complementary_run_bp: usize,
    pub right_3prime_complementary_run_bp: usize,
    pub max_3prime_complementary_run_bp: usize,
    pub primer_3prime_extension_risk: bool,
    pub status: String,
    pub summary: String,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Shared oligo-interaction QC report for supplied PCR/qPCR assay oligos.
pub struct OligoQcReport {
    pub schema: String,
    pub assay_kind: String,
    pub method: String,
    pub method_reference: String,
    pub status: String,
    pub summary: String,
    pub oligo_count: usize,
    pub interaction_count: usize,
    #[serde(default)]
    pub oligos: Vec<OligoQcOligoRecord>,
    #[serde(default)]
    pub interactions: Vec<OligoQcInteractionRecord>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
/// Distance/overlap summary for one primer pair relative to the required core
/// ROI of a simple-PCR design.
pub struct PrimerPairCoreGeometry {
    pub left_distance_from_core_bp: usize,
    pub right_distance_from_core_bp: usize,
    pub left_overlap_into_core_bp: usize,
    pub right_overlap_into_core_bp: usize,
}

impl PrimerPairCoreGeometry {
    pub fn flanks_core_cleanly(&self) -> bool {
        self.left_overlap_into_core_bp == 0 && self.right_overlap_into_core_bp == 0
    }

    pub fn side_label(distance_bp: usize, overlap_bp: usize) -> String {
        if overlap_bp > 0 {
            format!("overlap {overlap_bp} bp")
        } else {
            format!("{distance_bp} bp")
        }
    }

    pub fn left_label(&self) -> String {
        Self::side_label(
            self.left_distance_from_core_bp,
            self.left_overlap_into_core_bp,
        )
    }

    pub fn right_label(&self) -> String {
        Self::side_label(
            self.right_distance_from_core_bp,
            self.right_overlap_into_core_bp,
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct PrimerInsertionIntent {
    pub requested_forward_3prime_end_0based_exclusive: usize,
    pub requested_reverse_3prime_start_0based: usize,
    pub forward_extension_5prime: String,
    pub reverse_extension_5prime: String,
    pub forward_window_start_0based: usize,
    pub forward_window_end_0based_exclusive: usize,
    pub reverse_window_start_0based: usize,
    pub reverse_window_end_0based_exclusive: usize,
    pub max_anchor_shift_bp: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerInsertionPairCompensation {
    pub rank: usize,
    pub forward_anchor_shift_bp: isize,
    pub reverse_anchor_shift_bp: isize,
    pub within_shift_budget: bool,
    pub compensable: bool,
    pub forward_compensation_5prime: String,
    pub reverse_compensation_5prime: String,
    pub compensated_forward_5prime_tail: String,
    pub compensated_reverse_5prime_tail: String,
    pub compensated_forward_sequence: String,
    pub compensated_reverse_sequence: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerInsertionContextReport {
    pub requested_forward_3prime_end_0based_exclusive: usize,
    pub requested_reverse_3prime_start_0based: usize,
    pub forward_extension_5prime: String,
    pub reverse_extension_5prime: String,
    pub forward_window_start_0based: usize,
    pub forward_window_end_0based_exclusive: usize,
    pub reverse_window_start_0based: usize,
    pub reverse_window_end_0based_exclusive: usize,
    pub max_anchor_shift_bp: usize,
    pub uncompensable_pair_count: usize,
    pub out_of_shift_budget_pair_count: usize,
    #[serde(default)]
    pub pairs: Vec<PrimerInsertionPairCompensation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct OverlapExtensionMutagenesisConstraints {
    pub overlap_bp: usize,
    pub outer_forward: PrimerDesignSideConstraint,
    pub outer_reverse: PrimerDesignSideConstraint,
    pub inner_forward: PrimerDesignSideConstraint,
    pub inner_reverse: PrimerDesignSideConstraint,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignRejectionSummary {
    pub out_of_window: usize,
    pub gc_or_tm_out_of_bounds: usize,
    pub non_unique_anneal: usize,
    pub amplicon_or_roi_failure: usize,
    pub primer_constraint_failure: usize,
    pub pair_constraint_failure: usize,
    pub pair_evaluation_limit_skipped: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignReport {
    pub schema: String,
    pub report_id: String,
    pub template: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub roi_start_0based: usize,
    pub roi_end_0based: usize,
    pub forward: PrimerDesignSideConstraint,
    pub reverse: PrimerDesignSideConstraint,
    #[serde(default)]
    pub pair_constraints: PrimerDesignPairConstraint,
    pub min_amplicon_bp: usize,
    pub max_amplicon_bp: usize,
    pub max_tm_delta_c: f64,
    pub max_pairs: usize,
    pub pair_count: usize,
    #[serde(default)]
    pub pairs: Vec<PrimerDesignPairRecord>,
    #[serde(default)]
    pub rejection_summary: PrimerDesignRejectionSummary,
    #[serde(default)]
    pub backend: PrimerDesignBackendInfo,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub insertion_context: Option<PrimerInsertionContextReport>,
}

impl PrimerDesignReport {
    pub fn pair_core_geometry(&self, pair: &PrimerDesignPairRecord) -> PrimerPairCoreGeometry {
        PrimerPairCoreGeometry {
            left_distance_from_core_bp: self
                .roi_start_0based
                .saturating_sub(pair.forward.end_0based_exclusive),
            right_distance_from_core_bp: pair
                .reverse
                .start_0based
                .saturating_sub(self.roi_end_0based),
            left_overlap_into_core_bp: pair
                .forward
                .end_0based_exclusive
                .saturating_sub(self.roi_start_0based),
            right_overlap_into_core_bp: self
                .roi_end_0based
                .saturating_sub(pair.reverse.start_0based),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignBackendInfo {
    pub requested: String,
    pub used: String,
    pub fallback_reason: Option<String>,
    pub primer3_executable: Option<String>,
    pub primer3_version: Option<String>,
    pub primer3_explain: Option<String>,
    pub primer3_request_boulder_io: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct Primer3PreflightReport {
    pub backend: String,
    pub configured_executable: Option<String>,
    pub used_default_executable: bool,
    pub executable: String,
    pub resolved_path: Option<String>,
    pub working_directory: Option<String>,
    pub reachable: bool,
    pub version_probe_ok: bool,
    pub status_code: Option<i32>,
    pub version: Option<String>,
    pub detail: Option<String>,
    pub error: Option<String>,
    pub probe_time_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignReportSummary {
    pub report_id: String,
    pub template: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub roi_start_0based: usize,
    pub roi_end_0based: usize,
    pub pair_count: usize,
    #[serde(default)]
    pub backend_used: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Source trace for one oligo-order line item.
pub struct OligoOrderLineProvenance {
    pub source_kind: String,
    pub report_id: String,
    pub report_schema: String,
    pub template: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pair_rank: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assay_rank: Option<usize>,
    pub role: String,
    #[serde(default)]
    pub source_coordinates_0based: Vec<SequenceRange0Based>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One preserved procurement row in an oligo order form.
pub struct OligoOrderLineItem {
    pub line_id: String,
    pub line_no: usize,
    pub name: String,
    pub role: String,
    pub sequence_5_to_3: String,
    pub length_nt: usize,
    #[serde(default)]
    pub modifications: Vec<String>,
    pub scale: String,
    pub purification: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub notes: Option<String>,
    pub provenance: OligoOrderLineProvenance,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Exact procurement duplicate group. Creation never collapses these rows.
pub struct OligoOrderDuplicateGroup {
    pub group_id: String,
    #[serde(default)]
    pub line_ids: Vec<String>,
    pub sequence_5_to_3: String,
    #[serde(default)]
    pub modifications: Vec<String>,
    pub scale: String,
    pub purification: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Same-sequence reuse group across different procurement settings.
pub struct OligoOrderSequenceReuseGroup {
    pub group_id: String,
    #[serde(default)]
    pub line_ids: Vec<String>,
    pub sequence_5_to_3: String,
    #[serde(default)]
    pub procurement_tuple_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Explicit duplicate-review status for one order form.
pub struct OligoOrderDuplicateReview {
    pub status: String,
    pub default_action: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reviewer: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reviewed_at_unix_ms: Option<u128>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// First-class, reviewable batch artifact for planned oligo procurement.
pub struct OligoOrderForm {
    pub schema: String,
    pub form_id: String,
    #[serde(default)]
    pub target_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_note: Option<String>,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
    #[serde(default)]
    pub line_items: Vec<OligoOrderLineItem>,
    #[serde(default)]
    pub duplicate_groups: Vec<OligoOrderDuplicateGroup>,
    #[serde(default)]
    pub sequence_reuse_groups: Vec<OligoOrderSequenceReuseGroup>,
    pub duplicate_review: OligoOrderDuplicateReview,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Compact list row for persisted oligo order forms.
pub struct OligoOrderFormSummary {
    pub form_id: String,
    #[serde(default)]
    pub target_label: String,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
    pub line_count: usize,
    pub duplicate_group_count: usize,
    pub sequence_reuse_group_count: usize,
    pub duplicate_review_status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Generic JSON creation payload for a persisted oligo order form.
pub struct OligoOrderFormCreateRequest {
    pub schema: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub form_id: Option<String>,
    #[serde(default)]
    pub target_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_note: Option<String>,
    #[serde(default)]
    pub line_items: Vec<OligoOrderLineItem>,
    #[serde(default)]
    pub scale: String,
    #[serde(default)]
    pub purification: String,
    #[serde(default)]
    pub modifications: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RestrictionCloningSingleSiteSuggestion {
    pub enzyme: String,
    pub cut_position_0based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RestrictionCloningDirectedPairSuggestion {
    pub order_source: String,
    pub forward_enzyme: String,
    pub reverse_enzyme: String,
    pub forward_cut_position_0based: usize,
    pub reverse_cut_position_0based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RestrictionCloningVectorEnzymeSuggestions {
    pub seq_id: String,
    #[serde(default)]
    pub selected_mcs: Vec<String>,
    #[serde(default)]
    pub other_unique: Vec<String>,
    #[serde(default)]
    pub missing_mcs: Vec<String>,
    #[serde(default)]
    pub recommended_single_site: Vec<RestrictionCloningSingleSiteSuggestion>,
    #[serde(default)]
    pub recommended_directed_pairs: Vec<RestrictionCloningDirectedPairSuggestion>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RestrictionCloningPcrHandoffSeedRequest {
    pub schema: String,
    pub primer_report_id: String,
    pub template: String,
    pub destination_vector_seq_id: String,
    pub pair_index: usize,
    pub pair_rank: usize,
    pub selected_pair: PrimerDesignPairRecord,
    pub selected_pair_core_geometry: PrimerPairCoreGeometry,
    pub mode: RestrictionCloningPcrHandoffMode,
    pub forward_enzyme: String,
    pub reverse_enzyme: String,
    pub forward_leader_5prime: String,
    pub reverse_leader_5prime: String,
    pub selection_source: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub suggestion_order_source: Option<String>,
    pub vector_suggestions: RestrictionCloningVectorEnzymeSuggestions,
    pub operation: Operation,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RestrictionCloningPcrDigestCompatibility {
    pub status: String,
    pub order_source: String,
    #[serde(default)]
    pub vector_site_count_by_enzyme: BTreeMap<String, usize>,
    #[serde(default)]
    pub insert_site_count_by_enzyme: BTreeMap<String, usize>,
    #[serde(default)]
    pub expected_insert_site_count_by_enzyme: BTreeMap<String, usize>,
    #[serde(default)]
    pub vector_cut_position_0based_by_enzyme: BTreeMap<String, usize>,
    pub vector_sites_unique_ok: bool,
    pub insert_sites_unique_to_tails_ok: bool,
    pub directed_order_ok: bool,
    pub termini_compatible: bool,
    pub forward_end_geometry: String,
    pub reverse_end_geometry: String,
    #[serde(default)]
    pub blocking_errors: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RestrictionCloningPcrWorkflowHints {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub staged_workflow: Option<Workflow>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pcr_advanced_operation: Option<Operation>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub insert_digest_operation: Option<Operation>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vector_digest_operation: Option<Operation>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ligation_operation_snippet: Option<serde_json::Value>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub insert_fragment_hint: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vector_fragment_hint: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RestrictionCloningPcrHandoffReport {
    pub schema: String,
    pub report_id: String,
    pub template: String,
    pub primer_report_id: String,
    pub pair_index: usize,
    pub pair_rank: usize,
    pub destination_vector_seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub mode: RestrictionCloningPcrHandoffMode,
    pub forward_enzyme: String,
    pub reverse_enzyme: String,
    pub forward_leader_5prime: String,
    pub reverse_leader_5prime: String,
    pub original_forward: PrimerDesignPrimerRecord,
    pub original_reverse: PrimerDesignPrimerRecord,
    pub extended_forward: PrimerDesignPrimerRecord,
    pub extended_reverse: PrimerDesignPrimerRecord,
    pub extended_forward_seq_id: String,
    pub extended_reverse_seq_id: String,
    pub tailed_amplicon_seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub handoff_container_id: Option<String>,
    pub predicted_tailed_amplicon_length_bp: usize,
    pub predicted_tailed_amplicon_5prime: String,
    pub predicted_tailed_amplicon_3prime: String,
    pub extended_pair_complementary_run_bp: usize,
    pub extended_pair_3prime_complementary_run_bp: usize,
    #[serde(default)]
    pub compatibility: RestrictionCloningPcrDigestCompatibility,
    #[serde(default)]
    pub workflow_hints: RestrictionCloningPcrWorkflowHints,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RestrictionCloningPcrHandoffReportSummary {
    pub report_id: String,
    pub template: String,
    pub primer_report_id: String,
    pub pair_index: usize,
    pub pair_rank: usize,
    pub destination_vector_seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub mode: String,
    pub forward_enzyme: String,
    pub reverse_enzyme: String,
    pub compatibility_status: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum UniprotPeptideComparisonMode {
    DirectCompare,
    GlobalAlignment,
    #[default]
    Unavailable,
}

impl UniprotPeptideComparisonMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::DirectCompare => "direct_compare",
            Self::GlobalAlignment => "global_alignment",
            Self::Unavailable => "unavailable",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum UniprotVariantFeatureEvidenceStatus {
    Supported,
    Contradicted,
    Unaligned,
    #[default]
    OutsideComparableCoverage,
}

impl UniprotVariantFeatureEvidenceStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Supported => "supported",
            Self::Contradicted => "contradicted",
            Self::Unaligned => "unaligned",
            Self::OutsideComparableCoverage => "outside_comparable_coverage",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum UniprotProjectionAuditRowStatus {
    #[default]
    Consistent,
    Warning,
    Mismatch,
    MissingEvidence,
}

impl UniprotProjectionAuditRowStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Consistent => "consistent",
            Self::Warning => "warning",
            Self::Mismatch => "mismatch",
            Self::MissingEvidence => "missing_evidence",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotEnsemblLinkedXref {
    pub transcript_id: Option<String>,
    pub protein_id: Option<String>,
    pub gene_id: Option<String>,
    pub isoform_id: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotEnsemblLinkResolutionRow {
    pub transcript_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    #[serde(default)]
    pub matched_xrefs: Vec<UniprotEnsemblLinkedXref>,
    #[serde(default)]
    pub normalized_xref_transcript_ids: Vec<String>,
    #[serde(default)]
    pub normalized_xref_protein_ids: Vec<String>,
    #[serde(default)]
    pub normalized_xref_gene_ids: Vec<String>,
    pub status: String,
    #[serde(default)]
    pub diagnostics: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotEnsemblLinkReport {
    pub schema: String,
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_filter: Option<String>,
    #[serde(default)]
    pub rows: Vec<UniprotEnsemblLinkResolutionRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct UniprotTranscriptExonContribution {
    pub ordinal: usize,
    pub exon_start_1based: usize,
    pub exon_end_1based: usize,
    pub exon_nt: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coding_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coding_end_1based: Option<usize>,
    pub coding_nt: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotProjectionTranscriptAccountingRow {
    pub transcript_id: String,
    pub transcript_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    pub strand: String,
    #[serde(default)]
    pub contributing_exons: Vec<UniprotTranscriptExonContribution>,
    pub contributing_exon_nt_sum: usize,
    pub untranslated_5prime_nt: usize,
    pub untranslated_3prime_nt: usize,
    pub translated_nt: usize,
    pub translated_nt_divisible_by_3: bool,
    pub expected_aa_count: usize,
    pub derived_protein_length_aa: usize,
    pub derived_protein_sequence: String,
    pub uniprot_aa_count: usize,
    pub init_met_declared: bool,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotProjectionTranscriptAccountingReport {
    pub schema: String,
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_filter: Option<String>,
    #[serde(default)]
    pub rows: Vec<UniprotProjectionTranscriptAccountingRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct UniprotEnsemblExonBoundaryDifference {
    pub ordinal: usize,
    pub side: String,
    pub current_coordinate_1based: usize,
    pub ensembl_coordinate_1based: usize,
    pub note: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotEnsemblExonCompareRow {
    pub transcript_id: String,
    pub transcript_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_entry_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_transcript_id: Option<String>,
    #[serde(default)]
    pub current_contributing_exons: Vec<UniprotTranscriptExonContribution>,
    #[serde(default)]
    pub ensembl_contributing_exons: Vec<UniprotTranscriptExonContribution>,
    #[serde(default)]
    pub matched_exon_ordinals: Vec<usize>,
    #[serde(default)]
    pub missing_in_ensembl_ordinals: Vec<usize>,
    #[serde(default)]
    pub excess_in_ensembl_ordinals: Vec<usize>,
    #[serde(default)]
    pub boundary_differences: Vec<UniprotEnsemblExonBoundaryDifference>,
    pub status: String,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotEnsemblExonCompareReport {
    pub schema: String,
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_filter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_entry_id: Option<String>,
    #[serde(default)]
    pub rows: Vec<UniprotEnsemblExonCompareRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct UniprotPeptideDirectMismatch {
    pub position_1based: usize,
    pub gentle_residue: String,
    pub uniprot_residue: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_residue: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotVariantFeatureEvidenceRow {
    pub feature_key: String,
    pub feature_note: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub start_aa: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub end_aa: Option<usize>,
    pub status: UniprotVariantFeatureEvidenceStatus,
    pub detail: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotEnsemblPeptideCompareRow {
    pub transcript_id: String,
    pub transcript_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_entry_id: Option<String>,
    pub derived_protein_length_aa: usize,
    pub uniprot_protein_length_aa: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_protein_length_aa: Option<usize>,
    pub comparison_mode: UniprotPeptideComparisonMode,
    pub init_met_declared: bool,
    #[serde(default)]
    pub direct_mismatches: Vec<UniprotPeptideDirectMismatch>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub global_alignment: Option<SequenceAlignmentReport>,
    #[serde(default)]
    pub variant_feature_evidence: Vec<UniprotVariantFeatureEvidenceRow>,
    pub status: String,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotEnsemblPeptideCompareReport {
    pub schema: String,
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_filter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_entry_id: Option<String>,
    #[serde(default)]
    pub rows: Vec<UniprotEnsemblPeptideCompareRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotProjectionAuditEmailDraft {
    pub subject: String,
    pub body: String,
    #[serde(default)]
    pub failing_transcript_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotProjectionAuditRow {
    pub transcript_id: String,
    pub transcript_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    pub status: UniprotProjectionAuditRowStatus,
    #[serde(default)]
    pub mismatch_reasons: Vec<String>,
    pub link_resolution: UniprotEnsemblLinkResolutionRow,
    pub accounting: UniprotProjectionTranscriptAccountingRow,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_exon_compare: Option<UniprotEnsemblExonCompareRow>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peptide_compare: Option<UniprotEnsemblPeptideCompareRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotProjectionAuditReport {
    pub schema: String,
    pub report_id: String,
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_filter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_entry_id: Option<String>,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default)]
    pub rows: Vec<UniprotProjectionAuditRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub maintainer_email_draft: Option<UniprotProjectionAuditEmailDraft>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotProjectionAuditReportSummary {
    pub report_id: String,
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_entry_id: Option<String>,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub transcript_count: usize,
    pub failing_transcript_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotProjectionAuditParityRow {
    pub transcript_id: String,
    pub direct_status: UniprotProjectionAuditRowStatus,
    pub composed_status: UniprotProjectionAuditRowStatus,
    pub statuses_match: bool,
    pub accounting_match: bool,
    pub mismatch_reason_match: bool,
    pub comparison_mode_match: bool,
    #[serde(default)]
    pub differences: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotProjectionAuditParityReport {
    pub schema: String,
    pub report_id: String,
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_filter: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_entry_id: Option<String>,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub direct_report: UniprotProjectionAuditReport,
    pub composed_report: UniprotProjectionAuditReport,
    #[serde(default)]
    pub rows: Vec<UniprotProjectionAuditParityRow>,
    pub email_draft_transcripts_match: bool,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotProjectionAuditParityReportSummary {
    pub report_id: String,
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_entry_id: Option<String>,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub transcript_count: usize,
    pub divergent_transcript_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct QpcrAssayRuleFlags {
    pub roi_covered: bool,
    pub amplicon_size_in_range: bool,
    pub primer_tm_delta_in_range: bool,
    pub probe_inside_amplicon: bool,
    pub probe_tm_delta_in_range: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Exact source-template interval used by a transcript-aware primer or amplicon.
pub struct SequenceRange0Based {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One primer/probe binding site detected on a transcript-derived cDNA
/// template.
pub struct CdnaAssayPrimerHit {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub binding_sequence: String,
    pub oligo_binding_strand: String,
    pub mismatch_count: usize,
    #[serde(default)]
    pub source_ranges_0based: Vec<SequenceRange0Based>,
    #[serde(default)]
    pub covered_junction_labels: Vec<String>,
    pub spans_junction: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One PCR/qPCR amplicon detected on a transcript-derived cDNA template.
pub struct CdnaAssayProduct {
    pub amplicon_start_0based: usize,
    pub amplicon_end_0based_exclusive: usize,
    pub amplicon_length_bp: usize,
    pub forward_hit_index: usize,
    pub reverse_hit_index: usize,
    #[serde(default)]
    pub probe_hit_indices: Vec<usize>,
    #[serde(default)]
    pub source_ranges_0based: Vec<SequenceRange0Based>,
    #[serde(default)]
    pub covered_junction_labels: Vec<String>,
    pub spans_junction: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_equivalent_start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_equivalent_end_0based_exclusive: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_equivalent_length_bp: Option<usize>,
    #[serde(default)]
    pub genomic_carryover_risk: String,
    #[serde(default)]
    pub genomic_carryover_rationale: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One exon segment in transcript-local cDNA coordinates, preserving the source
/// interval that identifies the exon across isoforms.
pub struct CdnaAssayTranscriptExonSegment {
    pub exon_ordinal: usize,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub source_start_0based: usize,
    pub source_end_0based_exclusive: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Per-transcript result row for a cDNA PCR/qPCR assay test.
pub struct CdnaAssayTranscriptResult {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_path: Option<String>,
    pub strand: String,
    pub cdna_length_bp: usize,
    pub status: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exon_segments: Vec<CdnaAssayTranscriptExonSegment>,
    #[serde(default)]
    pub forward_hits: Vec<CdnaAssayPrimerHit>,
    #[serde(default)]
    pub reverse_hits: Vec<CdnaAssayPrimerHit>,
    #[serde(default)]
    pub probe_hits: Vec<CdnaAssayPrimerHit>,
    #[serde(default)]
    pub products: Vec<CdnaAssayProduct>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Report-level construct/oligo length summary for a tested cDNA assay.
pub struct CdnaAssayConstructLengthSummary {
    pub forward_primer_length_bp: usize,
    pub reverse_primer_length_bp: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub probe_length_bp: Option<usize>,
    #[serde(default)]
    pub detected_amplicon_lengths_bp: Vec<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_detected_amplicon_length_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_detected_amplicon_length_bp: Option<usize>,
    pub summary: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Report-level summary of residual genomic-DNA carryover risk for a cDNA
/// PCR/qPCR assay.
pub struct CdnaAssayGenomicCarryoverRiskSummary {
    pub risk_level: String,
    pub summary: String,
    pub product_count: usize,
    pub low_risk_product_count: usize,
    pub medium_risk_product_count: usize,
    pub high_risk_product_count: usize,
    pub unknown_risk_product_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_genomic_equivalent_length_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_genomic_equivalent_length_bp: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable SVG map showing where a tested PCR/qPCR assay is functional across
/// transcript cDNA templates.
pub struct CdnaAssayTranscriptMap {
    pub schema: String,
    pub artifact_id: String,
    pub media_type: String,
    pub title: String,
    pub summary: String,
    pub width_px: usize,
    pub height_px: usize,
    pub row_count: usize,
    pub shown_transcript_count: usize,
    pub omitted_transcript_count: usize,
    pub column_count: usize,
    pub rows_per_column: usize,
    pub transcript_order: CdnaAssayTranscriptOrder,
    pub coordinate_mode: CdnaAssayTranscriptMapCoordinateMode,
    pub product_count: usize,
    pub svg: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Deterministic cDNA PCR/qPCR assay-test report shared by shell/CLI/agents.
pub struct CdnaAssayTestReport {
    pub schema: String,
    pub assay_kind: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub template_source_kind: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_paths: Vec<String>,
    pub source_seq_id: String,
    pub source_feature_id: usize,
    pub group_label: String,
    pub strand: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub requested_transcript_id: Option<String>,
    pub forward_primer: String,
    pub reverse_primer: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub probe: Option<String>,
    #[serde(default)]
    pub construct_lengths: CdnaAssayConstructLengthSummary,
    #[serde(default)]
    pub genomic_carryover_risk: CdnaAssayGenomicCarryoverRiskSummary,
    #[serde(default)]
    pub oligo_qc: OligoQcReport,
    pub max_mismatches: usize,
    pub require_3prime_exact_bases: usize,
    pub min_amplicon_bp: usize,
    pub max_amplicon_bp: usize,
    pub transcript_count: usize,
    pub detected_transcript_count: usize,
    pub product_count: usize,
    pub transcript_order: CdnaAssayTranscriptOrder,
    pub transcript_map_coordinate_mode: CdnaAssayTranscriptMapCoordinateMode,
    pub overall_status: String,
    pub summary: String,
    #[serde(default)]
    pub transcript_results: Vec<CdnaAssayTranscriptResult>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_map: Option<CdnaAssayTranscriptMap>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One source-reference interval that supports a qPCR panel oligo.
///
/// `source_*` coordinates are local to the GENtle source sequence. The start is
/// 0-based and the end is exclusive; the paired 1-based fields are inclusive
/// display coordinates for the same interval. `genomic_*` coordinates are
/// 1-based inclusive absolute reference coordinates when a genome anchor is
/// available.
pub struct TranscriptQpcrPanelSourceRange {
    pub source_start_0based: usize,
    pub source_end_0based_exclusive: usize,
    pub source_start_1based: usize,
    pub source_end_1based_inclusive: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_end_1based: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Primer/probe row used by transcript qPCR panel summaries.
///
/// `primer.start_0based` and `primer.end_0based_exclusive` use the
/// `binding_coordinate_space` named here: either the source reference sequence
/// or the selected transcript cDNA template. Source/genomic coordinates are
/// repeated in `source_ranges` for adapter-neutral rendering.
pub struct TranscriptQpcrPanelOligoRecord {
    pub role: String,
    pub primer: PrimerDesignPrimerRecord,
    pub binding_coordinate_space: String,
    pub source_ranges: Vec<TranscriptQpcrPanelSourceRange>,
    pub source_range_label: String,
    pub reference_strand: String,
    pub spans_junction: bool,
    #[serde(default)]
    pub covered_junction_labels: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Per-transcript characteristic-primer row for a transcript qPCR panel.
pub struct TranscriptQpcrPanelTranscriptRow {
    /// Zero-based feature index into the source sequence feature table.
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub cdna_length_bp: usize,
    pub characteristic_status: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub characteristic_forward: Option<TranscriptQpcrPanelOligoRecord>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub realized_specificity_evidence: Option<QpcrTranscriptSpecificityEvidence>,
    pub exact_target_hit_count: usize,
    pub exact_competitor_hit_count: usize,
    #[serde(default)]
    pub exact_hit_transcript_ids: Vec<String>,
    #[serde(default)]
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Shared report for a qPCR panel table: common assay components plus
/// transcript-characteristic forward primers.
pub struct TranscriptQpcrPanelReport {
    pub schema: String,
    pub source_seq_id: String,
    /// Zero-based feature index into the source sequence feature table.
    pub source_feature_id: usize,
    pub group_label: String,
    pub strand: String,
    pub transcript_count: usize,
    pub shared_qpcr_report_id: String,
    pub shared_assay_rank: usize,
    pub shared_forward: TranscriptQpcrPanelOligoRecord,
    pub shared_reverse: TranscriptQpcrPanelOligoRecord,
    pub shared_probe: TranscriptQpcrPanelOligoRecord,
    #[serde(default)]
    pub transcript_rows: Vec<TranscriptQpcrPanelTranscriptRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Optional transcript-aware targeting request for one qPCR design operation.
pub struct QpcrTranscriptTargeting {
    pub source_feature_id: usize,
    #[serde(default)]
    pub mode: QpcrTranscriptTargetingMode,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specificity_evidence: Option<QpcrTranscriptSpecificityEvidence>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Persisted transcript-context summary for one accepted qPCR assay.
pub struct QpcrTranscriptAssayContext {
    pub assay_class_label: String,
    pub explanation: String,
    pub probe_placement: String,
    pub design_transcript_feature_id: usize,
    pub design_transcript_id: String,
    pub design_transcript_label: String,
    pub support_transcript_count: usize,
    pub support_transcript_fraction: f64,
    #[serde(default)]
    pub supported_transcript_ids: Vec<String>,
    #[serde(default)]
    pub forward_source_ranges_0based: Vec<SequenceRange0Based>,
    #[serde(default)]
    pub reverse_source_ranges_0based: Vec<SequenceRange0Based>,
    #[serde(default)]
    pub probe_source_ranges_0based: Vec<SequenceRange0Based>,
    #[serde(default)]
    pub amplicon_source_ranges_0based: Vec<SequenceRange0Based>,
    #[serde(default)]
    pub covered_junction_labels: Vec<String>,
    pub forward_spans_junction: bool,
    pub reverse_spans_junction: bool,
    pub probe_spans_junction: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_equivalent_start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_equivalent_end_0based_exclusive: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genomic_equivalent_length_bp: Option<usize>,
    #[serde(default)]
    pub genomic_carryover_risk: String,
    #[serde(default)]
    pub genomic_carryover_rationale: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_distinguishing_primer: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub realized_specificity_evidence: Option<QpcrTranscriptSpecificityEvidence>,
    pub satisfies_requested_targeting: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Report-level transcript-aware outcome summary shared across adapters.
pub struct QpcrTranscriptTargetingResult {
    pub source_feature_id: usize,
    #[serde(default)]
    pub mode: QpcrTranscriptTargetingMode,
    pub group_label: String,
    pub strand: String,
    pub transcript_count_considered: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub realized_specificity_evidence: Option<QpcrTranscriptSpecificityEvidence>,
    pub selected_support_transcript_count: usize,
    pub selected_support_transcript_fraction: f64,
    pub used_shared_support_fallback: bool,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct QpcrAssayRecord {
    pub rank: usize,
    pub score: f64,
    pub forward: PrimerDesignPrimerRecord,
    pub reverse: PrimerDesignPrimerRecord,
    pub probe: PrimerDesignPrimerRecord,
    pub amplicon_start_0based: usize,
    pub amplicon_end_0based_exclusive: usize,
    pub amplicon_length_bp: usize,
    pub primer_tm_delta_c: f64,
    pub probe_tm_delta_c: f64,
    pub rule_flags: QpcrAssayRuleFlags,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_context: Option<QpcrTranscriptAssayContext>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct QpcrDesignRejectionSummary {
    pub primer_pair: PrimerDesignRejectionSummary,
    pub probe_out_of_window: usize,
    pub probe_gc_or_tm_out_of_bounds: usize,
    pub probe_non_unique_anneal: usize,
    pub probe_or_assay_failure: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct QpcrDesignReport {
    pub schema: String,
    pub report_id: String,
    pub template: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub roi_start_0based: usize,
    pub roi_end_0based: usize,
    pub forward: PrimerDesignSideConstraint,
    pub reverse: PrimerDesignSideConstraint,
    pub probe: PrimerDesignSideConstraint,
    #[serde(default)]
    pub pair_constraints: PrimerDesignPairConstraint,
    pub min_amplicon_bp: usize,
    pub max_amplicon_bp: usize,
    pub max_tm_delta_c: f64,
    pub max_probe_tm_delta_c: f64,
    pub max_assays: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_targeting: Option<QpcrTranscriptTargeting>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_targeting_result: Option<QpcrTranscriptTargetingResult>,
    pub assay_count: usize,
    #[serde(default)]
    pub best_assay_probe_placement: String,
    #[serde(default)]
    pub best_assay_summary: String,
    #[serde(default)]
    pub assays: Vec<QpcrAssayRecord>,
    #[serde(default)]
    pub rejection_summary: QpcrDesignRejectionSummary,
    #[serde(default)]
    pub backend: PrimerDesignBackendInfo,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct QpcrDesignReportSummary {
    pub report_id: String,
    pub template: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub roi_start_0based: usize,
    pub roi_end_0based: usize,
    pub assay_count: usize,
    #[serde(default)]
    pub best_assay_probe_placement: String,
    #[serde(default)]
    pub best_assay_summary: String,
    #[serde(default)]
    pub backend_used: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One persisted transcript-native protein derivation row.
///
/// The report keeps both the created protein sequence id and the shared
/// transcript/CDS derivation summary so lineage, GUI reopen paths, and future
/// shell/agent surfaces can inspect the same deterministic translation
/// decision without re-deriving biology locally.
pub struct ProteinDerivationReportRow {
    pub protein_seq_id: String,
    pub transcript_feature_id: usize,
    pub protein_name: String,
    pub derivation: TranscriptProteinDerivation,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Persisted transcript-native protein derivation artifact.
///
/// This is the computational provenance companion to the derived protein
/// sequence nodes created by `DeriveProteinSequences`: the sequence entries
/// remain the first-class biological products, while this report captures the
/// transcript selection, coding-span resolution mode, and report-level
/// operation/run linkage needed for lineage-visible audit/reopen flows.
pub struct ProteinDerivationReport {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default)]
    pub requested_feature_ids: Vec<usize>,
    #[serde(default)]
    pub selected_feature_ids: Vec<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub feature_query: Option<SequenceFeatureQuery>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scope: Option<SplicingScopePreset>,
    pub effective_output_prefix: String,
    pub derived_count: usize,
    #[serde(default)]
    pub rows: Vec<ProteinDerivationReportRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Compact lineage/listing summary for one persisted protein-derivation report.
pub struct ProteinDerivationReportSummary {
    pub report_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub effective_output_prefix: String,
    pub derived_count: usize,
    pub derivation_mode_summary: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Persisted reverse-translation artifact for one protein-to-coding-DNA run.
///
/// The created coding sequence remains the first-class biological product in
/// project state. This report is the computational provenance companion that
/// captures the translation table, codon-bias profile resolution, optional Tm
/// steering, and operation/run linkage needed for lineage-visible audit and
/// reopen behavior.
pub struct ReverseTranslationReport {
    pub schema: String,
    pub report_id: String,
    pub protein_seq_id: String,
    pub coding_seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub requested_output_id: Option<String>,
    pub effective_output_id: String,
    pub protein_length_aa: usize,
    pub coding_length_bp: usize,
    pub translation_table: usize,
    pub translation_table_label: String,
    pub translation_table_source: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub translation_context_organism: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub translation_context_organelle: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub requested_speed_profile: Option<TranslationSpeedProfile>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resolved_speed_profile: Option<TranslationSpeedProfile>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resolved_speed_profile_source: Option<TranslationSpeedProfileSource>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub translation_speed_reference_species: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub speed_mark: Option<TranslationSpeedMark>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_anneal_tm_c: Option<f64>,
    pub anneal_window_bp: usize,
    pub preferred_synonymous_choice_count: usize,
    pub alternative_synonymous_choice_count: usize,
    pub fallback_unknown_codon_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gc_fraction: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub realized_anneal_tm_c: Option<f64>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Compact lineage/listing summary for one persisted reverse-translation report.
pub struct ReverseTranslationReportSummary {
    pub report_id: String,
    pub protein_seq_id: String,
    pub coding_seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub protein_length_aa: usize,
    pub coding_length_bp: usize,
    pub translation_table: usize,
    pub speed_profile_summary: String,
    pub diagnostics_summary: String,
}
/// Compact lineage/listing summary for one persisted construct-reasoning graph.
///
/// This keeps the graph itself as the canonical portable reasoning artifact,
/// while giving GUI/CLI lineage surfaces a stable, cheap-to-list record with
/// the counts and objective labels needed to present it as a first-class
/// computational contribution.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ConstructReasoningGraphSummary {
    pub graph_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub objective_id: String,
    pub objective_title: String,
    pub objective_goal: String,
    pub evidence_count: usize,
    pub decision_count: usize,
    pub candidate_count: usize,
    #[serde(default)]
    pub contains_protein_to_dna_handoff: bool,
    #[serde(default)]
    pub protein_to_dna_handoff_candidate_count: usize,
    #[serde(default)]
    pub protein_to_dna_source_protein_seq_ids: Vec<String>,
    #[serde(default)]
    pub summary_lines: Vec<String>,
    #[serde(default)]
    pub warning_lines: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GuideU6TerminatorWindow {
    SpacerOnly,
    #[default]
    SpacerPlusTail,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GuideOligoExportFormat {
    CsvTable,
    PlateCsv,
    Fasta,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GuideOligoPlateFormat {
    #[default]
    Plate96,
    Plate384,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideSet {
    pub guide_set_id: String,
    pub guides: Vec<GuideCandidate>,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideCandidate {
    pub guide_id: String,
    pub seq_id: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub strand: String,
    pub protospacer: String,
    pub pam: String,
    pub nuclease: String,
    pub cut_offset_from_protospacer_start: usize,
    pub rank: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuideSetSummary {
    pub guide_set_id: String,
    pub guide_count: usize,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct GuidePracticalFilterConfig {
    pub gc_min: Option<f64>,
    pub gc_max: Option<f64>,
    pub max_homopolymer_run: Option<usize>,
    #[serde(default)]
    pub max_homopolymer_run_per_base: HashMap<String, usize>,
    pub reject_ambiguous_bases: bool,
    pub avoid_u6_terminator_tttt: bool,
    pub u6_terminator_window: GuideU6TerminatorWindow,
    pub max_dinucleotide_repeat_units: Option<usize>,
    #[serde(default)]
    pub forbidden_motifs: Vec<String>,
    pub required_5prime_base: Option<String>,
    pub allow_5prime_g_extension: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideFilterReason {
    pub code: String,
    pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuidePracticalFilterResult {
    pub guide_id: String,
    pub passed: bool,
    #[serde(default)]
    pub reasons: Vec<GuideFilterReason>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub metrics: HashMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuidePracticalFilterReport {
    pub guide_set_id: String,
    pub generated_at_unix_ms: u128,
    pub config: GuidePracticalFilterConfig,
    pub passed_count: usize,
    pub rejected_count: usize,
    #[serde(default)]
    pub results: Vec<GuidePracticalFilterResult>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoTemplateSpec {
    pub template_id: String,
    pub description: String,
    pub forward_prefix: String,
    pub forward_suffix: String,
    pub reverse_prefix: String,
    pub reverse_suffix: String,
    pub reverse_uses_reverse_complement_of_spacer: bool,
    pub uppercase_output: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoRecord {
    pub guide_id: String,
    pub rank: Option<usize>,
    pub forward_oligo: String,
    pub reverse_oligo: String,
    #[serde(default)]
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoSet {
    pub oligo_set_id: String,
    pub guide_set_id: String,
    pub generated_at_unix_ms: u128,
    pub template: GuideOligoTemplateSpec,
    pub apply_5prime_g_extension: bool,
    #[serde(default)]
    pub records: Vec<GuideOligoRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoExportReport {
    pub guide_set_id: String,
    pub oligo_set_id: String,
    pub format: String,
    pub path: String,
    pub exported_records: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideProtocolExportReport {
    pub guide_set_id: String,
    pub oligo_set_id: String,
    pub path: String,
    pub guide_count: usize,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Which stored planning profile layer a read/write command refers to.
///
/// The effective merge order is `global -> confirmed_agent_overlay ->
/// project_override`.
pub enum PlanningProfileScope {
    #[default]
    Global,
    ProjectOverride,
    ConfirmedAgentOverlay,
    Effective,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Review state for an advisory planning suggestion.
pub enum PlanningSuggestionStatus {
    #[default]
    Pending,
    Accepted,
    Rejected,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Availability/cost note for one inventory item referenced by planning.
pub struct PlanningInventoryItem {
    pub available: bool,
    pub unit_cost: Option<f64>,
    pub procurement_business_days: Option<f64>,
    pub note: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Availability/queue note for one machine or platform used in planning.
pub struct PlanningMachineAvailability {
    pub available: bool,
    pub queue_business_days: f64,
    pub run_cost_per_hour: Option<f64>,
    pub note: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// User- or agent-supplied planning profile describing local capabilities.
///
/// This is the main durable contract for inventory, procurement latency, and
/// machine availability. Start here when tracing `planning profile ...`.
pub struct PlanningProfile {
    pub schema: String,
    pub profile_id: Option<String>,
    pub currency: Option<String>,
    pub procurement_business_days_default: f64,
    pub capabilities: Vec<String>,
    pub inventory: HashMap<String, PlanningInventoryItem>,
    pub machine_availability: HashMap<String, PlanningMachineAvailability>,
    pub notes: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Optimization weights and hard guardrails used by planning estimates.
pub struct PlanningObjective {
    pub schema: String,
    pub biological_intent: Option<String>,
    pub weight_time: f64,
    pub weight_cost: f64,
    pub weight_local_fit: f64,
    pub max_cost: Option<f64>,
    pub max_time_hours: Option<f64>,
    pub required_capabilities: Vec<String>,
    pub helper_profile_id: Option<String>,
    pub preferred_routine_families: Vec<String>,
    pub enforce_guardrails: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Deterministic estimate payload produced by planning evaluation.
pub struct PlanningEstimate {
    pub schema: String,
    pub estimated_time_hours: f64,
    pub estimated_cost: f64,
    pub local_fit_score: f64,
    pub composite_meta_score: f64,
    pub passes_guardrails: bool,
    pub guardrail_failures: Vec<String>,
    pub explanation: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Read-only cloning strategy and helper/vector planning consultation.
pub struct PlanningCloningConsultation {
    pub schema: String,
    pub profile_scope: String,
    pub seq_id: Option<String>,
    pub biological_intent: String,
    pub objective_summary: serde_json::Value,
    pub local_profile_summary: serde_json::Value,
    pub available_helper_vectors: Vec<PlanningCloningHelperVectorSummary>,
    pub available_host_profiles: Vec<PlanningCloningHostProfileSummary>,
    pub strategy_candidates: Vec<PlanningCloningStrategyCandidate>,
    pub vector_candidates: Vec<PlanningCloningVectorCandidate>,
    pub missing_questions: Vec<PlanningCloningMissingQuestion>,
    pub local_constraints: Vec<PlanningCloningLocalConstraint>,
    pub warnings: Vec<String>,
    pub suggested_next_actions: Vec<PlanningCloningSuggestedNextAction>,
    pub suggested_sync_payload: Option<serde_json::Value>,
    pub text_report: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Best routine candidate for one cloning routine family.
pub struct PlanningCloningStrategyCandidate {
    pub rank: usize,
    pub family: String,
    pub routine_id: String,
    pub title: String,
    pub status: String,
    pub summary: Option<String>,
    pub estimate: PlanningEstimate,
    pub rationale: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Helper/vector catalog row exposed to cloning consultation clients.
pub struct PlanningCloningHelperVectorSummary {
    pub helper_id: String,
    pub description: Option<String>,
    pub summary: Option<String>,
    pub helper_kind: Option<String>,
    pub host_systems: Vec<String>,
    pub structured_tags: Vec<String>,
    pub offered_functions: Vec<String>,
    pub routine_family_hints: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Host/strain catalog row exposed to cloning consultation clients.
pub struct PlanningCloningHostProfileSummary {
    pub profile_id: String,
    pub species: String,
    pub strain: String,
    pub aliases: Vec<String>,
    pub genotype_tags: Vec<String>,
    pub phenotype_tags: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Ranked helper/vector candidate based only on structured catalog fields.
pub struct PlanningCloningVectorCandidate {
    pub rank: usize,
    pub helper_id: String,
    pub score: f64,
    pub helper_kind: Option<String>,
    pub host_systems: Vec<String>,
    pub routine_family_hints: Vec<String>,
    pub rationale: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Stable follow-up question for underspecified cloning planning requests.
pub struct PlanningCloningMissingQuestion {
    pub question_id: String,
    pub prompt: String,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Local infrastructure note that affects strategy/vector planning confidence.
pub struct PlanningCloningLocalConstraint {
    pub constraint_id: String,
    pub status: String,
    pub summary: String,
    pub details: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Safe next action suggested by a cloning consultation.
pub struct PlanningCloningSuggestedNextAction {
    pub action_id: String,
    pub label: String,
    pub shell_line: String,
    pub rationale: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Read-only protein-expression planning handoff for underspecified yield requests.
pub struct ProteinExpressionHandoffReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub status: String,
    pub biological_intent: String,
    pub product_definition: ProteinExpressionProductDefinition,
    pub product_readiness: ProteinExpressionProductReadiness,
    pub sequence_context: Option<ProteinExpressionSequenceContext>,
    pub cds_assessment: ProteinExpressionCdsAssessment,
    pub tag_assessment: ProteinExpressionTagAssessment,
    pub host_chassis_candidates: Vec<ProteinExpressionHostChassisCandidate>,
    pub vector_route_candidates: Vec<ProteinExpressionVectorRouteCandidate>,
    pub missing_questions: Vec<PlanningCloningMissingQuestion>,
    pub service_handoff_candidates: Vec<ProteinExpressionServiceHandoffCandidate>,
    pub warnings: Vec<String>,
    pub suggested_next_actions: Vec<PlanningCloningSuggestedNextAction>,
    pub text_report: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Product context supplied to a protein-expression handoff.
pub struct ProteinExpressionProductDefinition {
    pub seq_id: Option<String>,
    pub sequence_present: bool,
    pub sequence_name: Option<String>,
    pub length_bp: Option<usize>,
    pub feature_count: Option<usize>,
    pub readiness: ProteinExpressionProductReadiness,
    pub product_metric: String,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Readiness summary for deciding whether a product can enter expression review.
pub struct ProteinExpressionProductReadiness {
    pub status: String,
    pub usable_sequence_context: bool,
    pub usable_cds_context: bool,
    pub translation_possible: bool,
    pub review_gate: String,
    pub blockers: Vec<String>,
    pub review_notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Sequence-level context consumed by a read-only protein-expression handoff.
pub struct ProteinExpressionSequenceContext {
    pub seq_id: Option<String>,
    pub sequence_name: Option<String>,
    pub molecule_type: Option<String>,
    pub nucleotide_length: Option<usize>,
    pub protein_length_aa: Option<usize>,
    pub feature_count: usize,
    pub gc_percent: Option<f64>,
    pub gc_window_bp: Option<usize>,
    pub gc_min_percent: Option<f64>,
    pub gc_max_percent: Option<f64>,
    pub ambiguous_base_count: usize,
    pub ambiguous_bases: Vec<String>,
    pub annotation_summaries: Vec<ProteinExpressionFeatureSummary>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Compact annotation row relevant to protein-expression handoff review.
pub struct ProteinExpressionFeatureSummary {
    pub feature_id: usize,
    pub kind: String,
    pub label: Option<String>,
    pub start_0based: Option<usize>,
    pub end_0based_exclusive: Option<usize>,
    pub strand: Option<String>,
    pub length_nt: Option<usize>,
    pub selected_qualifiers: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Read-only CDS/ORF sanity assessment for a protein-expression product.
pub struct ProteinExpressionCdsAssessment {
    pub context_source: String,
    pub plausible_cds: bool,
    pub translation_possible: bool,
    pub nucleotide_length: Option<usize>,
    pub protein_length_aa: Option<usize>,
    pub codon_count: Option<usize>,
    pub translation_table: Option<usize>,
    pub start_codon: Option<String>,
    pub stop_codon: Option<String>,
    pub starts_with_atg: Option<bool>,
    pub has_terminal_stop: Option<bool>,
    pub has_internal_stops: bool,
    pub internal_stop_count: usize,
    pub ambiguous_codon_count: usize,
    pub length_multiple_of_three: Option<bool>,
    pub annotated_cds_features: Vec<ProteinExpressionFeatureSummary>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Read-only tag/signalling annotation context for protein-expression review.
pub struct ProteinExpressionTagAssessment {
    pub annotated_tag_count: usize,
    pub annotated_signal_peptide_count: usize,
    pub annotated_tags: Vec<ProteinExpressionFeatureSummary>,
    pub tag_policy_known: bool,
    pub missing_inputs: Vec<String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Candidate expression chassis row.
pub struct ProteinExpressionHostChassisCandidate {
    pub rank: usize,
    pub chassis_id: String,
    pub label: String,
    pub status: String,
    pub rationale: Vec<String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Candidate construct/vector route row.
pub struct ProteinExpressionVectorRouteCandidate {
    pub rank: usize,
    pub route_id: String,
    pub label: String,
    pub status: String,
    pub rationale: Vec<String>,
    pub missing_inputs: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Provider-neutral service handoff preview row.
pub struct ProteinExpressionServiceHandoffCandidate {
    pub provider: String,
    pub service_kind: String,
    pub status: String,
    pub example_request_path: String,
    pub draft_request_preview: serde_json::Value,
    pub shell_line: String,
    pub rationale: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Advisory pull/push suggestion awaiting explicit user acceptance or rejection.
pub struct PlanningSuggestion {
    pub schema: String,
    pub suggestion_id: String,
    pub status: PlanningSuggestionStatus,
    pub direction: String,
    pub source: String,
    pub confidence: Option<f64>,
    pub snapshot_id: Option<String>,
    pub message: Option<String>,
    pub profile_patch: Option<PlanningProfile>,
    pub objective_patch: Option<PlanningObjective>,
    pub diff: serde_json::Value,
    pub created_at_unix_ms: u128,
    pub resolved_at_unix_ms: Option<u128>,
    pub rejection_reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// High-level sync status for the planning suggestion channel.
pub struct PlanningSyncStatus {
    pub schema: String,
    pub pending_suggestion_count: usize,
    pub last_pull_at_unix_ms: Option<u128>,
    pub last_push_at_unix_ms: Option<u128>,
    pub last_source: Option<String>,
    pub last_snapshot_id: Option<String>,
    pub last_error: Option<String>,
}

/// Stored workflow-macro template shared by shell, GUI helpers, and future MCP.
///
/// A workflow macro expands to one workflow script plus typed parameters and
/// declared input/output ports for lineage/protocol reporting.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct WorkflowMacroTemplate {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameters: Vec<WorkflowMacroTemplateParam>,
    #[serde(default)]
    pub input_ports: Vec<WorkflowMacroTemplatePort>,
    #[serde(default)]
    pub output_ports: Vec<WorkflowMacroTemplatePort>,
    #[serde(default = "default_cloning_macro_template_schema")]
    pub template_schema: String,
    pub script: String,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

fn default_cloning_macro_template_schema() -> String {
    CLONING_MACRO_TEMPLATE_SCHEMA.to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One typed parameter exposed by a workflow macro template.
pub struct WorkflowMacroTemplateParam {
    pub name: String,
    pub default_value: Option<String>,
    pub required: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Declared input or output port on a workflow macro template.
pub struct WorkflowMacroTemplatePort {
    pub port_id: String,
    pub kind: String,
    pub required: bool,
    pub cardinality: String,
    pub description: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Listing row for one workflow macro template.
pub struct WorkflowMacroTemplateSummary {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameter_count: usize,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Stored candidate-macro template that expands to candidate-shell commands.
pub struct CandidateMacroTemplate {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameters: Vec<CandidateMacroTemplateParam>,
    pub script: String,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Listing row for one candidate macro template.
pub struct CandidateMacroTemplateSummary {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameter_count: usize,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
/// Engine-owned helper/objective-derived routine preference context.
pub struct RoutinePreferenceContextRecord {
    pub helper_profile_id: Option<String>,
    pub construct_reasoning_seq_id: Option<String>,
    pub helper_resolution_status: String,
    pub explicit_preferred_routine_families: Vec<String>,
    pub helper_derived_preferred_routine_families: Vec<String>,
    pub construct_strategy_derived_preferred_routine_families: Vec<String>,
    pub variant_derived_preferred_routine_families: Vec<String>,
    pub effective_preferred_routine_families: Vec<String>,
    pub helper_offered_functions: Vec<String>,
    pub helper_component_labels: Vec<String>,
    pub variant_effect_tags: Vec<String>,
    pub variant_suggested_assay_ids: Vec<String>,
    pub rationale: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
/// Candidate routine score snapshot captured in routine-decision traces.
pub struct RoutineDecisionTraceCandidateScore {
    pub routine_id: String,
    pub routine_title: Option<String>,
    pub routine_family: String,
    pub passes_guardrails: bool,
    pub estimated_time_hours: Option<f64>,
    pub estimated_cost: Option<f64>,
    pub local_fit_score: Option<f64>,
    pub composite_meta_score: Option<f64>,
    pub routine_family_alignment_bonus: Option<f64>,
    pub routine_family_alignment_sources: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
/// Suggested workflow/candidate macro template aligned to current planning context.
pub struct MacroTemplateSuggestion {
    pub macro_kind: String,
    pub template_name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub score: f64,
    pub matched_routine_families: Vec<String>,
    pub matched_terms: Vec<String>,
    pub rationale: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelDomainSpec {
    pub name: String,
    pub start_aa: usize,
    pub end_aa: usize,
    #[serde(default)]
    pub color_hex: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelCurationInfo {
    #[serde(default)]
    pub source_kind: Option<String>,
    #[serde(default)]
    pub source_label: Option<String>,
    #[serde(default)]
    pub evidence: Vec<String>,
    #[serde(default)]
    pub validation_tags: Vec<String>,
    #[serde(default)]
    pub public_database_status: Option<String>,
    #[serde(default)]
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
/// One upstream or local source record that supports a curated isoform panel.
pub struct IsoformPanelEvidenceRecord {
    pub evidence_id: String,
    pub source_type: String,
    pub accession: String,
    #[serde(default)]
    pub version: Option<String>,
    #[serde(default)]
    pub label: Option<String>,
    #[serde(default)]
    pub description: Option<String>,
    #[serde(default)]
    pub url: Option<String>,
    #[serde(default)]
    pub sequence_path: Option<String>,
    #[serde(default)]
    pub sequence_length_bp: Option<usize>,
    #[serde(default)]
    pub sequence_sha256: Option<String>,
    #[serde(default)]
    pub cds_start_1based: Option<usize>,
    #[serde(default)]
    pub cds_end_1based: Option<usize>,
    #[serde(default)]
    pub retrieved_on: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
/// Small key/value metric attached to one curated isoform-panel evaluation row.
pub struct IsoformPanelEvaluationMetric {
    pub key: String,
    pub value: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
/// One comparison row inside a curated isoform-panel evaluation.
pub struct IsoformPanelEvaluationRow {
    #[serde(default)]
    pub isoform_id: Option<String>,
    #[serde(default)]
    pub evidence_id: Option<String>,
    #[serde(default)]
    pub compared_to: Option<String>,
    pub status: String,
    pub summary: String,
    #[serde(default)]
    pub metrics: Vec<IsoformPanelEvaluationMetric>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
/// A stored curation/comparison note for a panel, such as ENA-vs-Ensembl checks.
pub struct IsoformPanelEvaluationRecord {
    pub evaluation_id: String,
    #[serde(default)]
    pub title: Option<String>,
    #[serde(default)]
    pub created_on: Option<String>,
    #[serde(default)]
    pub method: Option<String>,
    pub status: String,
    pub summary: String,
    #[serde(default)]
    pub source_evidence_ids: Vec<String>,
    #[serde(default)]
    pub isoform_ids: Vec<String>,
    #[serde(default)]
    pub rows: Vec<IsoformPanelEvaluationRow>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum IsoformTranscriptGeometryMode {
    #[default]
    Exon,
    Cds,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelIsoformSpec {
    pub isoform_id: String,
    #[serde(default)]
    pub label: Option<String>,
    #[serde(default)]
    pub transcript_ids: Vec<String>,
    #[serde(default)]
    pub transactivation_class: Option<String>,
    #[serde(default)]
    pub expected_length_aa: Option<usize>,
    #[serde(default)]
    pub reference_start_aa: Option<usize>,
    #[serde(default)]
    pub reference_end_aa: Option<usize>,
    #[serde(default)]
    pub curation: Option<IsoformPanelCurationInfo>,
    #[serde(default)]
    pub domains: Vec<IsoformPanelDomainSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelResource {
    pub schema: String,
    pub panel_id: String,
    pub gene_symbol: String,
    #[serde(default)]
    pub transcript_geometry_mode: IsoformTranscriptGeometryMode,
    #[serde(default)]
    pub assembly: Option<String>,
    #[serde(default)]
    pub source: Option<String>,
    #[serde(default)]
    pub notes: Option<String>,
    #[serde(default)]
    pub curation: Option<IsoformPanelCurationInfo>,
    #[serde(default)]
    pub evidence: Vec<IsoformPanelEvidenceRecord>,
    #[serde(default)]
    pub evaluations: Vec<IsoformPanelEvaluationRecord>,
    #[serde(default)]
    pub isoforms: Vec<IsoformPanelIsoformSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelValidationIssue {
    pub severity: String,
    pub code: String,
    pub message: String,
    pub isoform_id: Option<String>,
    pub transcript_probe: Option<String>,
    pub domain_name: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelValidationIsoformSummary {
    pub isoform_id: String,
    pub label: String,
    pub transcript_probe_count: usize,
    pub domain_count: usize,
    pub expected_length_aa: Option<usize>,
    pub max_domain_end_aa: Option<usize>,
    #[serde(default)]
    pub curation_source_kind: Option<String>,
    #[serde(default)]
    pub validation_tags: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelValidationReport {
    pub schema: String,
    pub path: String,
    pub panel_id: String,
    pub gene_symbol: String,
    pub assembly: Option<String>,
    pub isoform_count: usize,
    pub transcript_probe_count: usize,
    pub unique_transcript_probe_count: usize,
    pub domain_count: usize,
    #[serde(default)]
    pub curation_source_kind: Option<String>,
    #[serde(default)]
    pub curated_isoform_count: usize,
    #[serde(default)]
    pub evidence_count: usize,
    #[serde(default)]
    pub evaluation_count: usize,
    #[serde(default)]
    pub evaluation_row_count: usize,
    pub issue_count: usize,
    pub status: String,
    pub isoforms: Vec<IsoformPanelValidationIsoformSummary>,
    pub issues: Vec<IsoformPanelValidationIssue>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessRunBundleOperationInputSummary {
    pub op_id: String,
    pub run_id: String,
    pub operation: String,
    pub record_index: usize,
    #[serde(default)]
    pub sequence_ids: Vec<String>,
    #[serde(default)]
    pub container_ids: Vec<String>,
    #[serde(default)]
    pub arrangement_ids: Vec<String>,
    #[serde(default)]
    pub candidate_set_ids: Vec<String>,
    #[serde(default)]
    pub guide_set_ids: Vec<String>,
    #[serde(default)]
    pub genome_ids: Vec<String>,
    #[serde(default)]
    pub file_paths: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ProcessRunBundleInputs {
    pub root_sequence_ids: Vec<String>,
    pub referenced_sequence_ids: Vec<String>,
    pub referenced_container_ids: Vec<String>,
    pub referenced_arrangement_ids: Vec<String>,
    pub referenced_candidate_set_ids: Vec<String>,
    pub referenced_guide_set_ids: Vec<String>,
    pub referenced_genome_ids: Vec<String>,
    pub file_inputs: Vec<String>,
    pub operation_inputs: Vec<ProcessRunBundleOperationInputSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessRunBundleParameterOverride {
    pub op_id: String,
    pub run_id: String,
    pub record_index: usize,
    pub name: String,
    pub value: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ProcessRunBundleOutputs {
    pub created_seq_ids: Vec<String>,
    pub changed_seq_ids: Vec<String>,
    pub final_sequences: Vec<EngineSequenceSummary>,
    pub created_container_ids: Vec<String>,
    pub created_arrangement_ids: Vec<String>,
    pub exported_paths: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct LabAssistantMaterialRow {
    pub material_id: String,
    pub display_name: String,
    pub kind: String,
    pub source: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub length_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub topology: Option<String>,
    #[serde(default)]
    pub members: Vec<String>,
    #[serde(default)]
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct LabAssistantInstructionSection {
    pub heading: String,
    pub steps: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum LabAssistantInstructionsFormat {
    #[default]
    Markdown,
    Odt,
    Docx,
}

impl LabAssistantInstructionsFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Markdown => "markdown",
            Self::Odt => "odt",
            Self::Docx => "docx",
        }
    }

    pub fn extension(self) -> &'static str {
        match self {
            Self::Markdown => "md",
            Self::Odt => "odt",
            Self::Docx => "docx",
        }
    }

    pub fn from_token(value: &str) -> Option<Self> {
        match value.trim().to_ascii_lowercase().as_str() {
            "markdown" | "md" => Some(Self::Markdown),
            "odt" | "opendocument" | "open_document" => Some(Self::Odt),
            "docx" | "word" => Some(Self::Docx),
            _ => None,
        }
    }

    pub fn infer_from_path(path: &str) -> Self {
        let extension = std::path::Path::new(path)
            .extension()
            .and_then(|value| value.to_str())
            .unwrap_or_default();
        Self::from_token(extension).unwrap_or(Self::Markdown)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct LabAssistantVisualRow {
    pub visual_id: String,
    pub label: String,
    pub format: String,
    pub source: String,
    pub width_px: usize,
    pub height_px: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct LabAssistantInstructionsExport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub title: String,
    pub audience: String,
    pub output_path: String,
    pub output_format: LabAssistantInstructionsFormat,
    pub run_id_filter: Option<String>,
    pub selected_record_count: usize,
    pub material_rows: Vec<LabAssistantMaterialRow>,
    pub step_sections: Vec<LabAssistantInstructionSection>,
    pub embedded_visuals: Vec<LabAssistantVisualRow>,
    pub checkpoint_lines: Vec<String>,
    pub safety_lines: Vec<String>,
    pub record_keeping_lines: Vec<String>,
    pub warning_lines: Vec<String>,
    pub summary_lines: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTraceComparison {
    pub left_routine_id: String,
    pub right_routine_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTraceDisambiguationQuestion {
    pub question_id: String,
    pub question_text: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTraceDisambiguationAnswer {
    pub question_id: String,
    pub answer_text: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTracePreflightSnapshot {
    pub can_execute: bool,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
    pub contract_source: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTraceExportEvent {
    pub run_bundle_path: String,
    pub exported_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTrace {
    pub schema: String,
    pub trace_id: String,
    pub source: String,
    pub status: String,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
    pub goal_text: String,
    pub query_text: String,
    pub candidate_routine_ids: Vec<String>,
    pub routine_preference_context: Option<RoutinePreferenceContextRecord>,
    pub candidate_planning_scores: Vec<RoutineDecisionTraceCandidateScore>,
    pub selected_routine_id: Option<String>,
    pub selected_routine_title: Option<String>,
    pub selected_routine_family: Option<String>,
    pub macro_suggestions: Vec<MacroTemplateSuggestion>,
    pub alternatives_presented: Vec<String>,
    pub comparisons: Vec<RoutineDecisionTraceComparison>,
    pub disambiguation_questions_presented: Vec<RoutineDecisionTraceDisambiguationQuestion>,
    pub disambiguation_answers: Vec<RoutineDecisionTraceDisambiguationAnswer>,
    pub bindings_snapshot: BTreeMap<String, String>,
    pub preflight_history: Vec<RoutineDecisionTracePreflightSnapshot>,
    pub preflight_snapshot: Option<RoutineDecisionTracePreflightSnapshot>,
    pub execution_attempted: bool,
    pub execution_success: Option<bool>,
    pub transactional: Option<bool>,
    pub macro_instance_id: Option<String>,
    pub emitted_operation_ids: Vec<String>,
    pub execution_error: Option<String>,
    pub export_events: Vec<RoutineDecisionTraceExportEvent>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTraceStore {
    pub schema: String,
    pub traces: Vec<RoutineDecisionTrace>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
#[derive(Default)]
pub struct ProcessRunBundleConstructReasoningSummary {
    pub seq_id: String,
    pub graph_id: String,
    pub objective_id: String,
    pub objective_title: String,
    pub objective_goal: String,
    pub fact_types: Vec<String>,
    pub decision_types: Vec<String>,
    pub fact_statuses: BTreeMap<String, String>,
    pub host_profile_ids: Vec<String>,
    pub helper_profile_id: Option<String>,
    pub medium_conditions: Vec<String>,
    pub growth_condition_signals: Vec<String>,
    pub supported_selection_rule_ids: Vec<String>,
    pub variant_effect_tags: Vec<String>,
    pub suggested_variant_assay_ids: Vec<String>,
    pub summary_lines: Vec<String>,
    pub warning_lines: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ProcessRunBundleConstructReasoningExport {
    pub seq_ids_considered: Vec<String>,
    pub summaries: Vec<ProcessRunBundleConstructReasoningSummary>,
    pub graphs: Vec<ConstructReasoningGraph>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ProcessRunBundleExport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default)]
    pub run_id_filter: Option<String>,
    pub selected_record_count: usize,
    pub inputs: ProcessRunBundleInputs,
    #[serde(default)]
    pub parameter_overrides: Vec<ProcessRunBundleParameterOverride>,
    #[serde(default)]
    pub decision_traces: Vec<RoutineDecisionTrace>,
    pub operation_log: Vec<OperationRecord>,
    pub outputs: ProcessRunBundleOutputs,
    pub parameter_snapshot: serde_json::Value,
    pub construct_reasoning: ProcessRunBundleConstructReasoningExport,
}

impl Default for ProcessRunBundleExport {
    fn default() -> Self {
        Self {
            schema: String::new(),
            generated_at_unix_ms: 0,
            run_id_filter: None,
            selected_record_count: 0,
            inputs: ProcessRunBundleInputs::default(),
            parameter_overrides: vec![],
            decision_traces: vec![],
            operation_log: vec![],
            outputs: ProcessRunBundleOutputs::default(),
            parameter_snapshot: serde_json::Value::Null,
            construct_reasoning: ProcessRunBundleConstructReasoningExport::default(),
        }
    }
}
