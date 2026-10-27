//! Bound Phase 2 consumption and portable, independently inspectable structural evidence.

use super::*;
use crate::observed_transcript_structures as structures;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedTranscriptAnnotation {
    pub source_id: String,
    pub provider: Option<String>,
    pub reference: PrimerPayloadReference,
}

/// One genomic batch and explicit stable gene IDs; neither a symbol nor a default gene.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SummarizeObservedTranscriptStructuresRequest {
    pub schema: String,
    pub source: PrimerPayloadResultBinding,
    pub source_report: PrimerPayloadFile,
    pub batch_id: String,
    pub reference_gene: structures::ReferenceGene,
    pub comparison_genes: Vec<structures::ReferenceGene>,
    pub annotation_sources: Vec<ObservedTranscriptAnnotation>,
    pub required_flank_bases: usize,
    pub budgets: PrimerPayloadBudgets,
    pub kernel_budgets: structures::Budgets,
    pub output_dir: String,
    pub cancel_path: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedTranscriptUnavailable {
    pub file_id: String,
    pub record_index: u64,
    pub read_id: String,
    pub query_id: Option<String>,
    pub reason: String,
}

/// Counts are absent when the atomic analysis unit did not complete, never measured zeros.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedTranscriptCounts {
    pub records: usize,
    pub observations: usize,
    pub chain_families: usize,
    pub unavailable_views: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedTranscriptStructuresReport {
    pub schema: String,
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: SummarizeObservedTranscriptStructuresRequest,
    pub status: String,
    pub control_error: Option<String>,
    pub source_status: String,
    pub source_trace_status: String,
    pub source_query_status: String,
    pub cohort: PrimerPayloadCohort,
    pub source_batches_completed: Vec<String>,
    pub source_batches_not_completed: Vec<String>,
    pub source_inputs_missing: usize,
    pub source_query_units_missing: usize,
    pub source_excluded_anchor_count: usize,
    pub counts: Option<ObservedTranscriptCounts>,
    pub artifacts: Vec<RnaReadPanelArtifact>,
    pub nonclaims: Vec<String>,
}

/// Full analysis lives in a hashed artifact, not in an unbounded adapter response.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedTranscriptStructuresData {
    pub schema: String,
    pub request_sha256: String,
    pub analysis: Option<structures::Summary>,
    pub unavailable: Vec<ObservedTranscriptUnavailable>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedTranscriptStructuresResultBinding {
    pub schema: String,
    pub request_sha256: String,
    pub report_content_sha256: String,
    pub report_path: String,
    pub outcome: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedTranscriptStructuresInspection {
    pub report: ObservedTranscriptStructuresReport,
    pub artifacts_verified: bool,
    /// Bounded JSON rows from the requested saved table. No source re-analysis.
    pub table: String,
    pub rows: Vec<serde_json::Value>,
}
