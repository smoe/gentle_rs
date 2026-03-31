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
    SequenceFeatureRangeRelation, SequenceFeatureSortBy, SequenceFeatureStrandFilter, TfbsProgress,
};
use serde::{Deserialize, Serialize};

use super::{
    DisplaySettings, OpId, Operation, PrepareGenomeProgress, ProtocolCartoonTemplateBindings,
    RunId, SeqId, SequencingConfirmationReport,
};

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
/// Union of long-running operation progress events.
pub enum OperationProgress {
    Tfbs(TfbsProgress),
    GenomePrepare(PrepareGenomeProgress),
    GenomeTrackImport(GenomeTrackImportProgress),
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
