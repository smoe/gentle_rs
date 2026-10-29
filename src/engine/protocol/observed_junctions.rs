//! J1 saved-S1 assessment contracts; no S1 mutation, translation or biological defaults.

use super::*;
use crate::observed_transcript_structures as s1;

pub const OBSERVED_JUNCTION_REQUEST_SCHEMA: &str =
    "gentle.observed_splice_junction_assessment_request.v1";
pub const OBSERVED_JUNCTION_REPORT_SCHEMA: &str =
    "gentle.observed_splice_junction_assessment_report.v1";
pub const OBSERVED_JUNCTION_BINDING_SCHEMA: &str =
    "gentle.observed_splice_junction_assessment_result_binding.v1";
pub const OBSERVED_JUNCTION_DATA_SCHEMA: &str =
    "gentle.observed_splice_junction_assessment_data.v1";
pub const OBSERVED_SAMPLE_MANIFEST_SCHEMA: &str = "gentle.observed_sample_manifest.v1";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ObservedJunctionSource {
    S1Report {
        result: ObservedTranscriptStructuresResultBinding,
        report: PrimerPayloadFile,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ObservedJunctionReferenceAuthority {
    FullAssembly,
    DeclaredExcerpt {
        contig: String,
        offset_0based: usize,
        length: usize,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedJunctionReferenceFormat {
    FastaPlainV1,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedJunctionReference {
    pub file: PrimerPayloadFile,
    pub format: ObservedJunctionReferenceFormat,
    pub authority: ObservedJunctionReferenceAuthority,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSampleContext {
    pub sample_id: String,
    pub barcode: String,
    pub priming: String,
    pub file_id: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSampleBinding {
    pub context: ObservedSampleContext,
    pub library_id: String,
    pub cdna_prep_id: String,
    pub biological_sample_id: String,
    pub cell_line: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSampleManifest {
    pub schema: String,
    pub rows: Vec<ObservedSampleBinding>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedRecurrenceUnit {
    Library,
    CdnaPrep,
    BiologicalSample,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedRecurrenceCertainty {
    Definite,
    DefiniteOrPossible,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedRecurrencePolicy {
    pub unit: ObservedRecurrenceUnit,
    pub min_units: usize,
    pub certainty: ObservedRecurrenceCertainty,
    pub require_qualified_flanks: bool,
}

/// Explicit ceilings, never truncation, assay thresholds or throughput promises.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedJunctionBudgets {
    pub max_input_bytes: u64,
    pub max_reference_bytes: u64,
    pub max_manifest_bytes: u64,
    pub max_distinct_junctions: usize,
    pub max_annotated_junctions_compared: usize,
    pub max_near_hits_per_junction: usize,
    pub max_total_work: u64,
    pub max_output_bytes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AssessObservedSpliceJunctionsRequest {
    pub schema: String,
    pub source: ObservedJunctionSource,
    pub reference_sequence: Option<ObservedJunctionReference>,
    pub sample_manifest: Option<PrimerPayloadFile>,
    pub near_window_bp: usize,
    pub recurrence_policy: Option<ObservedRecurrencePolicy>,
    pub budgets: ObservedJunctionBudgets,
    pub output_dir: String,
    pub cancel_path: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedMotifClass {
    GtAg,
    GcAg,
    AtAc,
    Other,
    NotAssessable,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedJunctionMotif {
    pub reference_authority: String,
    pub reference_forward_first: Option<String>,
    pub reference_forward_last: Option<String>,
    pub reference_forward_motif_class: ObservedMotifClass,
    pub reference_reverse_motif_class: ObservedMotifClass,
    pub motif_class: ObservedMotifClass,
    pub opposite_orientation_motif_class: ObservedMotifClass,
    pub reason: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedAnnotatedJunctionHit {
    pub skipped_interval: s1::Interval,
    pub start_offset: i64,
    pub end_offset: i64,
    pub model_ids: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedModelJunctionEvents {
    pub transcript_id: String,
    pub labels: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedJunctionSourceAssessment {
    pub source_id: String,
    pub gene_id: Option<String>,
    pub relation: String,
    pub reason: String,
    pub model_strand: Option<char>,
    pub shared_end: Option<String>,
    pub near_hits: Vec<ObservedAnnotatedJunctionHit>,
    pub events: Vec<String>,
    pub model_events: Vec<ObservedModelJunctionEvents>,
}

/// Mutually exclusive record partitions within one certainty class.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedQualifiedCounts {
    pub qualified_only: usize,
    pub insufficient_only: usize,
    pub unassessed_only: usize,
    pub mixed: usize,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedJunctionSupportCounts {
    pub definite: ObservedQualifiedCounts,
    pub possible_only: ObservedQualifiedCounts,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedJunctionRecordSupport {
    pub record: s1::ReadRecord,
    pub certainty: String,
    pub qualification: String,
    pub contexts: Vec<ObservedSampleContext>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedGroupedJunctionSupport {
    pub dimension: String,
    pub id: String,
    pub counts: ObservedJunctionSupportCounts,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedJunctionRecurrence {
    pub status: String,
    pub reason: String,
    pub independent_units: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSpliceJunctionAssessment {
    pub junction_id: String,
    pub display_id: String,
    pub contig: String,
    pub skipped_interval: s1::Interval,
    pub intron_first_base_1based: usize,
    pub intron_last_base_1based: usize,
    pub donor_exon_last_base_1based: Option<usize>,
    pub acceptor_exon_first_base_1based: Option<usize>,
    pub orientation: Option<char>,
    pub orientation_reason: String,
    pub motif: ObservedJunctionMotif,
    pub sources: Vec<ObservedJunctionSourceAssessment>,
    pub support: ObservedJunctionSupportCounts,
    pub records: Vec<ObservedJunctionRecordSupport>,
    pub grouped_support: Vec<ObservedGroupedJunctionSupport>,
    pub recurrence: Option<ObservedJunctionRecurrence>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedChainAssessment {
    pub observation_id: String,
    pub record: s1::ReadRecord,
    pub chain_id: Option<String>,
    pub comparisons: Vec<s1::SourceComparison>,
    pub end_constraints: [Vec<s1::EndConstraint>; 2],
    pub alignment_clipped: [bool; 2],
    pub rna_end_status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedLibraryDenominator {
    pub library_id: String,
    pub total_records: usize,
    pub gene_coherent_records: usize,
    pub excluded_records: usize,
    pub unassessed_records: usize,
    pub context_coverage_complete: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedJunctionUnassessedContext {
    pub record: s1::ReadRecord,
    pub context: Option<ObservedSampleContext>,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedJunctionAssessmentData {
    pub schema: String,
    pub junctions: Vec<ObservedSpliceJunctionAssessment>,
    pub chains: Vec<ObservedChainAssessment>,
    pub libraries: Vec<ObservedLibraryDenominator>,
    pub unassessed_contexts: Vec<ObservedJunctionUnassessedContext>,
    pub unused_manifest_rows: Vec<ObservedSampleBinding>,
    pub source_unavailable: Vec<ObservedTranscriptUnavailable>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedJunctionAssessmentCounts {
    pub source_records: usize,
    pub distinct_junctions: usize,
    pub unassessed_contexts: usize,
    pub annotated_comparisons: usize,
    pub charged_input_bytes: u64,
    pub charged_work: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSpliceJunctionAssessmentReport {
    pub schema: String,
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: AssessObservedSpliceJunctionsRequest,
    pub status: String,
    pub control_error: Option<String>,
    pub source_status: String,
    pub counts: Option<ObservedJunctionAssessmentCounts>,
    pub artifacts: Vec<RnaReadPanelArtifact>,
    pub nonclaims: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSpliceJunctionResultBinding {
    pub schema: String,
    pub request_sha256: String,
    pub report_content_sha256: String,
    pub report_path: String,
    pub outcome: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedSpliceJunctionInspection {
    pub report: ObservedSpliceJunctionAssessmentReport,
    pub artifacts_verified: bool,
    pub table: String,
    pub rows: Vec<serde_json::Value>,
}
