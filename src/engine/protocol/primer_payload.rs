//! Explicit, reference-bound payload evidence; no aligner, assay or molecule identity inference.

use super::{RnaPrimerHit, RnaPrimerPanel, RnaPrimerTraceResultBinding, RnaReadPanelArtifact};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const PRIMER_PAYLOAD_PREPARE_SCHEMA: &str = "gentle.primer_payload_prepare_request.v1";
pub const PRIMER_PAYLOAD_QUERIES_SCHEMA: &str = "gentle.primer_payload_queries.v1";
pub const PRIMER_PAYLOAD_REQUEST_SCHEMA: &str = "gentle.primer_payload_alignment_request.v1";
pub const PRIMER_PAYLOAD_REPORT_SCHEMA: &str = "gentle.primer_payload_alignment_report.v1";
pub const PRIMER_PAYLOAD_BINDING_SCHEMA: &str = "gentle.primer_payload_result_binding.v1";
pub const PRIMER_PAYLOAD_ANNOTATION_SCHEMA: &str = "gentle.primer_payload_annotation.v1";

/// Raw file bytes, not a semantic JSON digest or an authenticated producer signature.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadFile {
    pub path: String,
    pub sha256: String,
    pub bytes: u64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerPayloadSourceMode {
    RetainedPrimaryFasta,
    ReplayVerifiedFastq,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerPayloadAnchorSelection {
    Primary,
    All,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerPayloadInternalPrimerPolicy {
    Mask,
    Unassessed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadCohort {
    pub cohort_id: String,
    pub label: String,
    /// True means all verified trace units, not all originally requested inputs.
    pub all_verified_units: bool,
    pub selection_provenance: String,
}

/// Explicit resource ceilings, shared by preparation and interpretation; never truncation limits.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadBudgets {
    pub max_input_bytes: u64,
    pub max_output_bytes: u64,
    pub max_records: usize,
    pub max_queries: usize,
    pub max_query_bases: usize,
    pub max_line_bytes: usize,
    pub max_alignments_per_query: usize,
    pub max_cigar_operations: usize,
    pub max_annotation_rows: usize,
    pub max_reference_context_bases: usize,
    pub max_chain_comparisons: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreparePrimerPayloadAlignmentQueriesRequest {
    pub schema: String,
    pub trace: RnaPrimerTraceResultBinding,
    pub file_ids: Vec<String>,
    pub cohort: PrimerPayloadCohort,
    pub source_mode: PrimerPayloadSourceMode,
    pub anchor_selection: PrimerPayloadAnchorSelection,
    pub internal_primer_policy: PrimerPayloadInternalPrimerPolicy,
    pub min_payload_bases: usize,
    pub budgets: PrimerPayloadBudgets,
    pub output_dir: String,
    pub cancel_path: Option<String>,
}

/// All intervals are half-open. Query bases are normalized into the observed primer direction.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadQuery {
    pub query_id: String,
    pub anchor_id: String,
    pub file_id: String,
    pub record_index: u64,
    pub read_id: String,
    pub sample_id: String,
    pub barcode: String,
    pub priming: String,
    pub primer_ids: Vec<String>,
    pub anchor: RnaPrimerHit,
    pub raw_read_length: usize,
    pub raw_payload_interval: (usize, usize),
    pub excluded_query_intervals: Vec<(usize, usize)>,
    pub query_length: usize,
    pub usable_bases: usize,
    pub sequence_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadUnassessed {
    pub file_id: String,
    pub record_index: u64,
    pub read_id: String,
    pub anchor_id: Option<String>,
    pub primer_ids: Vec<String>,
    pub reason: String,
}

/// Only fully verified source units may contribute. Missing units retain explicit coverage.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadQueries {
    pub schema: String,
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: PreparePrimerPayloadAlignmentQueriesRequest,
    pub status: String,
    pub control_error: Option<String>,
    pub source_trace_status: String,
    pub source_inputs_missing: usize,
    pub original_fastq_reverified: bool,
    pub primer_panel: RnaPrimerPanel,
    pub files_completed: Vec<String>,
    pub files_not_completed: Vec<String>,
    pub queries: Vec<PrimerPayloadQuery>,
    pub unassessed: Vec<PrimerPayloadUnassessed>,
    pub artifacts: Vec<RnaReadPanelArtifact>,
    pub nonclaims: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum PrimerPayloadReferenceScope {
    Transcriptome,
    Genome,
}

/// Normalized annotations explicitly associate reference targets and stable gene/transcript IDs.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadAnnotation {
    pub schema: String,
    pub reference_id: String,
    pub assembly: String,
    pub release: String,
    pub reference_fasta_sha256: String,
    pub coordinate_convention: String,
    pub transcripts: Vec<PrimerPayloadTranscript>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadTranscript {
    pub transcript_id: String,
    pub gene_ids: Vec<String>,
    pub target_name: String,
    pub target_length: usize,
    pub genomic_contig: Option<String>,
    pub genomic_strand: Option<char>,
    /// Ascending, nonoverlapping genomic exons; required for genomic assertions.
    pub genomic_exons: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadReference {
    pub reference_id: String,
    pub scope: PrimerPayloadReferenceScope,
    pub assembly: String,
    pub release: String,
    pub fasta: PrimerPayloadFile,
    pub annotation: PrimerPayloadFile,
}

/// No implicit equivalence between a Phase 1 locus name and a new reference gene ID.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadExpectedTarget {
    pub primer_id: String,
    pub reference_id: String,
    pub phase1_locus_id: String,
    pub phase1_transcript_feature_id: Option<usize>,
    pub gene_ids: Vec<String>,
    pub transcript_ids: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PrimerPayloadQueryAccounting {
    Mapped,
    ProducerReportedUnmapped,
    Unassessed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadAlignmentBatch {
    pub batch_id: String,
    pub reference_id: String,
    pub paf: PrimerPayloadFile,
    pub query_manifest_sha256: String,
    pub producer: String,
    pub producer_version: String,
    pub arguments: Vec<String>,
    pub producer_exit_code: i32,
    pub declared_complete: bool,
    pub competitor_search_limitations: String,
    /// Every selected query is accounted for; missing entries cannot become unmapped zeroes.
    pub accounting: BTreeMap<String, PrimerPayloadQueryAccounting>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadPolicy {
    pub strict_alignment_detail: bool,
    pub min_segment_bases: usize,
    pub min_identity_per_mille: usize,
    pub min_mapq: u8,
    pub dominant_fraction_per_mille: usize,
    pub min_cross_locus_segment_bases: usize,
    pub breakpoint_tolerance_bases: usize,
    pub terminal_seed_lengths: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SummarizePrimerPayloadAlignmentsRequest {
    pub schema: String,
    pub trace: RnaPrimerTraceResultBinding,
    pub queries: PrimerPayloadFile,
    pub query_manifest_sha256: String,
    pub references: Vec<PrimerPayloadReference>,
    pub expected_targets: Vec<PrimerPayloadExpectedTarget>,
    pub batches: Vec<PrimerPayloadAlignmentBatch>,
    pub policy: PrimerPayloadPolicy,
    pub budgets: PrimerPayloadBudgets,
    pub output_dir: String,
    pub cancel_path: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadSegment {
    pub query_interval: (usize, usize),
    pub target_name: String,
    pub target_interval: (usize, usize),
    pub strand: char,
    pub mapq: Option<u8>,
    pub reported_role: Option<String>,
    pub reported_identity_per_mille: usize,
    pub edits: Option<usize>,
    pub alignment_detail_verified: bool,
    /// Validated operation lengths in increasing reference order, not a reproduced search.
    pub alignment_operations: Vec<(char, usize)>,
    pub aligned_query_intervals: Vec<(usize, usize)>,
    pub gene_ids: Vec<String>,
    pub transcript_ids: Vec<String>,
    pub splice_edges: Vec<PrimerPayloadSpliceEdge>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadSpliceEdge {
    pub contig: String,
    pub donor: usize,
    pub acceptor: usize,
    pub strand: char,
    pub transcript_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadBoundary {
    pub left_query_end: usize,
    pub right_query_start: usize,
    pub left: PrimerPayloadEdge,
    pub right: PrimerPayloadEdge,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadEdge {
    pub target_name: String,
    pub coordinate: usize,
    pub strand: char,
    pub gene_ids: Vec<String>,
    pub transcript_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadSeed {
    pub length: usize,
    pub mismatches: Option<usize>,
    pub reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadReadEvidence {
    pub query_id: String,
    pub batch_id: String,
    pub reference_id: String,
    pub primer_ids: Vec<String>,
    pub file_id: String,
    pub record_index: u64,
    pub read_id: String,
    pub sample_id: String,
    pub barcode: String,
    pub priming: String,
    pub usable_bases: usize,
    pub expected_target_bases: Option<usize>,
    pub expected_target_fraction_per_mille: Option<usize>,
    pub classification: String,
    pub partner_gene_ids: Vec<String>,
    pub segments: Vec<PrimerPayloadSegment>,
    pub primer_edge_alternatives: Vec<PrimerPayloadEdge>,
    pub internal_query_boundaries: Vec<PrimerPayloadBoundary>,
    pub terminal_seeds: BTreeMap<String, Vec<PrimerPayloadSeed>>,
    pub cross_locus_payload_evidence: bool,
    pub reasons: Vec<String>,
    pub artifact_flags: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadPartnerProfile {
    pub reference_id: String,
    pub sample_id: String,
    pub barcode: String,
    pub priming: String,
    pub primer_ids: Vec<String>,
    pub partner_gene_ids: Vec<String>,
    pub supporting_reads: usize,
    pub eligible_reads: usize,
    pub assessed_reads: usize,
    pub ambiguous_reads: usize,
    pub unmapped_reads: usize,
    pub unassessed_reads: usize,
    pub excluded_reads: usize,
    pub anchor_observations: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadBreakpointCluster {
    pub reference_id: String,
    pub sample_id: String,
    pub barcode: String,
    pub priming: String,
    pub primer_ids: Vec<String>,
    pub target_name: String,
    pub strand: char,
    pub edge_type: String,
    pub partner_gene_ids: Vec<String>,
    pub minimum_edge: usize,
    pub maximum_edge: usize,
    pub representative_edge: usize,
    pub supporting_reads: usize,
    pub eligible_reads: usize,
    pub assessed_reads: usize,
    pub ambiguous_reads: usize,
    pub unmapped_reads: usize,
    pub unassessed_reads: usize,
    pub excluded_reads: usize,
}

/// Large read evidence lives in hashed TSV, not an unbounded agent JSON response.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadAlignmentReport {
    pub schema: String,
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: SummarizePrimerPayloadAlignmentsRequest,
    pub status: String,
    pub control_error: Option<String>,
    pub cohort: PrimerPayloadCohort,
    pub source_trace_status: String,
    pub source_inputs_missing: usize,
    pub source_query_status: String,
    pub source_query_units_missing: usize,
    pub source_excluded_anchor_count: usize,
    pub batches_completed: Vec<String>,
    pub batches_not_completed: Vec<String>,
    pub evidence_rows: usize,
    pub profiles: Vec<PrimerPayloadPartnerProfile>,
    pub breakpoint_clusters: Vec<PrimerPayloadBreakpointCluster>,
    pub artifacts: Vec<RnaReadPanelArtifact>,
    pub nonclaims: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadResultBinding {
    pub schema: String,
    pub kind: String,
    pub request_sha256: String,
    pub report_content_sha256: String,
    pub report_path: String,
    pub outcome: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerPayloadInspection {
    pub report: serde_json::Value,
    pub rows: Vec<serde_json::Value>,
    pub artifact_content_verified: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    content = "request",
    rename_all = "snake_case",
    deny_unknown_fields
)]
pub enum PrimerPayloadControlRequest {
    Queries(PreparePrimerPayloadAlignmentQueriesRequest),
    Alignments(SummarizePrimerPayloadAlignmentsRequest),
}
