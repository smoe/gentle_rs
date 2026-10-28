//! S2.1b: strict annotation-preparation contracts, not operation/adapter registration.

use super::*;
use crate::observed_transcript_structures::{ReferenceBinding, ReferenceGene};

// One semantic vocabulary shared by the pure consumer and the versioned receipt.
pub use crate::observed_transcript_structures::annotation::{
    Counts as ObservedAnnotationCounts, Diagnostic as ObservedAnnotationDiagnostic,
    Entity as ObservedAnnotationEntity, InputScope as ObservedAnnotationInputScope,
    Outcome as ObservedAnnotationOutcome, Policy as ObservedAnnotationPolicy,
    Reason as ObservedAnnotationReason, TranscriptRows as ObservedAnnotationTranscriptRows,
};

pub const REFERENCE_CONTIG_DICTIONARY_SCHEMA: &str = "gentle.reference_contig_dictionary.v1";
pub const OBSERVED_ANNOTATION_REQUEST_SCHEMA: &str =
    "gentle.observed_transcript_annotation_request.v1";
pub const OBSERVED_ANNOTATION_REPORT_SCHEMA: &str =
    "gentle.observed_transcript_annotation_report.v1";

/// Ordered exact-name dictionary; its hash does not verify the declared FASTA.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceContigDictionary {
    pub schema: String,
    pub reference_id: String,
    pub assembly: String,
    pub reference_fasta_sha256: String,
    pub contigs: Vec<ReferenceContigDictionaryEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceContigDictionaryEntry {
    pub name: String,
    pub length: u64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedAnnotationFormat {
    Gtf,
}

/// Declared provenance is retained verbatim, not authenticated or inferred from a header.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedAnnotationSourceDeclaration {
    pub scope: ObservedAnnotationInputScope,
    pub provider: Option<String>,
    pub declared_provider_release: Option<String>,
    pub declared_derivation: Option<String>,
    pub source_uri: Option<String>,
    pub notice: Option<String>,
}

/// Explicit preparation-only ceilings, not throughput promises or biological defaults.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedAnnotationBudgets {
    /// Total bound GTF plus dictionary bytes.
    pub max_input_bytes: u64,
    pub max_dictionary_bytes: usize,
    pub max_line_bytes: usize,
    pub max_lines: usize,
    pub max_dictionary_entries: usize,
    pub max_indexed_genes: usize,
    pub max_indexed_transcripts: usize,
    pub max_exon_rows: usize,
    pub max_retained_bytes: usize,
    pub max_work: usize,
    pub max_diagnostics: usize,
    pub max_output_bytes: u64,
}

/// One configurable exact gene in one source. All filesystem writes are external effects.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrepareObservedTranscriptAnnotationRequest {
    pub schema: String,
    pub format: ObservedAnnotationFormat,
    pub gtf: PrimerPayloadFile,
    pub dictionary: PrimerPayloadFile,
    pub reference_id: String,
    pub reference: ReferenceBinding,
    pub selector: ReferenceGene,
    pub source: ObservedAnnotationSourceDeclaration,
    pub policy: ObservedAnnotationPolicy,
    pub budgets: ObservedAnnotationBudgets,
    pub output_dir: String,
    pub cancel_path: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ObservedAnnotationExecutionStatus {
    Completed,
    Cancelled,
    ControlError,
}

/// Present only after the whole supplied input is verified and analysis completes.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedAnnotationAnalysis {
    pub grammar: String,
    pub outcome: ObservedAnnotationOutcome,
    pub coordinate_authority: String,
    pub raw_sha256: String,
    pub raw_bytes: usize,
    pub dictionary_entries: usize,
    pub snapshot_release: String,
    pub counts: ObservedAnnotationCounts,
    pub selected_reasons: BTreeMap<ObservedAnnotationReason, usize>,
    pub unrelated_reasons: BTreeMap<ObservedAnnotationReason, usize>,
    pub diagnostics: Vec<ObservedAnnotationDiagnostic>,
    pub diagnostic_count: usize,
    pub omitted_diagnostics: usize,
    pub indexed_genes: usize,
    pub indexed_transcripts: usize,
    pub gene_version: Option<String>,
    pub transcript_versions: BTreeMap<String, Option<String>>,
    pub included_contigs: Vec<String>,
    pub work: usize,
    pub retained_peak_bytes: usize,
}

/// Content digest uses canonical JSON with only `content_sha256` omitted.
/// Raw report bytes have a separate file binding returned by the direct engine method.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedTranscriptAnnotationReport {
    pub schema: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: PrepareObservedTranscriptAnnotationRequest,
    pub status: ObservedAnnotationExecutionStatus,
    pub control_error: Option<String>,
    pub dictionary_verified: bool,
    pub gtf_verified: bool,
    pub analysis: Option<ObservedAnnotationAnalysis>,
    pub annotation: Option<RnaReadPanelArtifact>,
    pub nonclaims: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservedTranscriptAnnotationPreparation {
    pub report: ObservedTranscriptAnnotationReport,
    pub report_file: PrimerPayloadFile,
}
