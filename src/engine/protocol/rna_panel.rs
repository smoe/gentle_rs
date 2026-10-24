//! Portable multi-locus seed-screen contracts; co-nomination is not whole-read compatibility.

use super::{GeneSetResolutionReport, SequenceGenomeAnchorSummary};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const RNA_PANEL_SCHEMA: &str = "gentle.rna_read_panel.v1";
pub const RNA_PANEL_BUILD_SCHEMA: &str = "gentle.rna_read_panel_build_request.v1";
pub const RNA_PANEL_SCREEN_LEGACY_SCHEMA: &str = "gentle.rna_read_panel_screen.v1";
pub const RNA_PANEL_SCREEN_SCHEMA: &str = "gentle.rna_read_panel_screen.v2";
pub const RNA_PANEL_SCREEN_REQUEST_SCHEMA: &str = "gentle.rna_read_panel_screen_request.v1";
pub const RNA_PANEL_SUMMARY_REQUEST_SCHEMA: &str = "gentle.rna_read_panel_summary_request.v1";
pub const RNA_PANEL_SUMMARY_LEGACY_SCHEMA: &str = "gentle.rna_read_panel_summary.v1";
pub const RNA_PANEL_SUMMARY_SCHEMA: &str = "gentle.rna_read_panel_summary.v2";

/// V2 count keys preserve segment classification, unlike v1 member-only keys.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelNominationSetKey {
    pub class: String,
    pub member_ids: Vec<String>,
}

/// Raw compressed bytes consumed can include decoder read-ahead; they are not verified EOF.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelProgress {
    pub file_id: String,
    pub sample_id: String,
    pub file_ordinal: usize,
    pub files_total: usize,
    pub processed_reads: u64,
    pub compressed_bytes_read: u64,
    pub compressed_bytes_total: u64,
    pub phase: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// One explicitly declared source group; labels are provenance, not authenticated identity.
pub struct RnaReadPanelSource {
    pub locus_id: String,
    pub seq_id: String,
    pub source_feature_id: usize,
    pub source: String,
    pub assembly: String,
    pub annotation_release: String,
    pub gene_id: String,
    pub gene_symbol: String,
    #[serde(default)]
    pub selected_transcript_ids: Vec<String>,
    /// Optional optimistic source binding, supplied by shared source discovery.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_source_record_sha256: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Build from loaded loci into a new portable file without mutating those loci.
pub struct BuildRnaReadPanelRequest {
    pub schema: String,
    pub sources: Vec<RnaReadPanelSource>,
    pub output_path: String,
    /// Extend this validated panel with new loci; recheck its loaded source hashes.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub base_panel_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_set: Option<RnaReadPanelGeneSetBinding>,
    /// Explicit project-graph persistence through the BuildRnaReadPanel operation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub register_in_project: Option<bool>,
}

/// Exact reference to the existing gene-set operand, not a second gene-list model.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelGeneSetBinding {
    pub resolution_id: String,
    pub resolution_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaReadPanelGeneSetChoice {
    pub binding: RnaReadPanelGeneSetBinding,
    pub resolution: GeneSetResolutionReport,
}

/// Project-owned immutable panel snapshot; external file existence is not implied.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaReadPanelProjectEntry {
    pub panel: RnaReadPanel,
    pub output_path: String,
    pub op_id: String,
    pub created_at_unix_ms: u128,
}

/// Editable source defaults from loaded annotations, not authenticated provider identity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaReadPanelSourceChoice {
    pub source: RnaReadPanelSource,
    pub transcript_ids: Vec<String>,
    /// Geometry failures remain visible rather than becoming an empty biological result.
    pub diagnostic: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelLocus {
    pub source: RnaReadPanelSource,
    /// Semantic loaded DNA + ordered features, not raw source-file bytes or mutable caches.
    pub source_record_sha256: String,
    pub genome_anchor: Option<SequenceGenomeAnchorSummary>,
    pub annotated_transcript_count: usize,
    pub annotated_transcript_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelTemplate {
    pub member_id: String,
    pub locus_id: String,
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub strand: String,
    pub sequence: String,
    pub sequence_sha256: String,
    /// Exon starts in transcript-oriented 0-based coordinates, excluding zero.
    pub junction_offsets: Vec<usize>,
    /// Source-sequence half-open exon intervals, in transcript order.
    pub source_exons: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Independent source namespaces and mature templates, bound by a canonical content digest.
pub struct RnaReadPanel {
    pub schema: String,
    pub content_sha256: String,
    pub loci: Vec<RnaReadPanelLocus>,
    pub templates: Vec<RnaReadPanelTemplate>,
    pub locus_count: usize,
    pub gene_count: usize,
    pub selected_transcript_count: usize,
    pub annotated_transcript_count: usize,
    /// Potential duplicate source views are disclosed, never silently counted as independent loci.
    pub overlapping_source_pairs: Vec<(String, String)>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_set_resolution: Option<Box<GeneSetResolutionReport>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parent_panel_sha256: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelInput {
    pub sample_id: String,
    pub barcode: String,
    pub priming: String,
    pub path: String,
    /// Caller attests copying this individual file has finished. Digest still verifies at EOF.
    pub declared_complete: bool,
    #[serde(default)]
    pub expected_md5: Option<String>,
    #[serde(default)]
    pub expected_sha256: Option<String>,
    /// Exact supplier manifest path when a basename is not unique.
    #[serde(default)]
    pub checksum_manifest_key: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
/// Uncalibrated nomination thresholds and hard budgets; overflow never drops associations.
pub struct RnaReadPanelPolicy {
    pub kmer_len: usize,
    pub min_chain_bases: u32,
    pub min_junction_flank_bp: usize,
    pub max_index_postings: usize,
    pub max_anchors_per_read: usize,
    pub max_equivalence_classes: usize,
    pub write_locus_fastas: bool,
}

impl Default for RnaReadPanelPolicy {
    fn default() -> Self {
        Self {
            kmer_len: 15,
            min_chain_bases: 45,
            min_junction_flank_bp: 5,
            max_index_postings: 2_000_000,
            max_anchors_per_read: 100_000,
            max_equivalence_classes: 4096,
            write_locus_fastas: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Explicit completed files and sample bindings; compressed checksums verify during streaming.
pub struct ScreenRnaReadPanelRequest {
    pub schema: String,
    pub panel_path: String,
    /// Bind the reviewed panel, not merely its mutable filesystem location.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_panel_sha256: Option<String>,
    pub files: Vec<RnaReadPanelInput>,
    #[serde(default)]
    pub md5_manifest_path: Option<String>,
    pub output_dir: String,
    #[serde(default)]
    pub policy: RnaReadPanelPolicy,
    /// Explicit cooperative cancellation; may be created externally while the screen runs.
    #[serde(default)]
    pub cancel_path: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Qualifying seed-chain evidence, never an exclusive transcript assignment.
pub struct RnaReadPanelCandidate {
    pub member_id: String,
    pub locus_id: String,
    pub transcript_id: String,
    pub read_orientation: String,
    pub chain_bases: u32,
    /// Original read coordinates, irrespective of orientation.
    pub read_interval_0based: (usize, usize),
    pub supported_junction_offsets: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelArtifact {
    pub path: String,
    pub sha256: String,
    pub bytes: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub locus_id: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelFileReceipt {
    pub file_id: String,
    pub input: RnaReadPanelInput,
    pub observed_md5: String,
    pub observed_sha256: String,
    pub compressed_bytes: u64,
    pub reads: u64,
    pub classes: BTreeMap<String, u64>,
    /// V2 keys encode NominationSetKey; v1 keys encode only member IDs. Neither is an isoform call.
    pub equivalence_classes: BTreeMap<String, u64>,
    pub artifacts: Vec<RnaReadPanelArtifact>,
}

/// Where a marker failure was observed, not a request to cancel.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadPanelControlErrorPhase {
    FileBoundary,
    Streaming,
    Finalization,
}

/// Post-start marker failure; retaining verified units does not authorize cancellation or prove completion.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelControlError {
    pub code: String,
    pub message: String,
    pub marker_path: String,
    pub phase: RnaReadPanelControlErrorPhase,
    pub file_id: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Verified units survive cancellation/control errors; neither outcome means a complete analysis.
pub struct RnaReadPanelScreenReport {
    pub schema: String,
    pub content_sha256: String,
    /// Canonical submitted request, before manifest resolution/normalization.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub submitted_request_sha256: Option<String>,
    pub status: String,
    pub request: ScreenRnaReadPanelRequest,
    pub panel_sha256: String,
    pub indexed_postings: usize,
    pub isoform_assignment_performed: bool,
    pub alignment_confirmation_performed: bool,
    pub files_completed: Vec<RnaReadPanelFileReceipt>,
    pub files_not_completed: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inputs_not_completed: Option<Vec<RnaReadPanelInput>>,
    /// Absent on historical v1/v2 reports and on successful or requested-cancel outcomes.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub control_error: Option<RnaReadPanelControlError>,
    pub nonclaims: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelInspection {
    pub schema: String,
    pub report_sha256: String,
    pub file_id: String,
    pub expected_artifact_sha256: String,
    /// A bounded prefix inspection does not re-hash the entire potentially huge TSV.
    pub artifact_content_verified: bool,
    pub rows: Vec<RnaReadPanelNomination>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelNomination {
    pub sample_id: String,
    pub barcode: String,
    pub file_id: String,
    pub expected_input_digest: String,
    pub record_index: u64,
    pub read_id: String,
    pub decompressed_record_offset: u64,
    pub class: String,
    pub candidates: Vec<RnaReadPanelCandidate>,
}

/// Consolidate explicit screen reports without reading original FASTQ or nomination sidecars.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExportRnaReadPanelSummaryRequest {
    pub schema: String,
    pub report_paths: Vec<String>,
    pub output_dir: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelSummaryRun {
    pub report_path: String,
    pub report_sha256: String,
    pub status: String,
    pub completed_files: usize,
    pub files_not_completed: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inputs_not_completed: Option<Vec<RnaReadPanelInput>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub control_error: Option<RnaReadPanelControlError>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelSummaryFile {
    pub report_sha256: String,
    pub file_id: String,
    pub sample_id: String,
    pub input_path: String,
    pub observed_sha256: String,
    pub observed_md5: String,
    pub reads: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelSummarySample {
    pub sample_id: String,
    pub barcode: String,
    pub priming: String,
    pub verified_files: usize,
    pub reads: u64,
    pub single_candidate_reads: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub co_nominated_reads: Option<u64>,
    /// Retained only when reading historical v1 summaries; not whole-read compatibility.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ambiguous_reads: Option<u64>,
    pub no_candidate_reads: u64,
    pub classes: BTreeMap<String, u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelSummaryTranscript {
    pub sample_id: String,
    pub gene_id: String,
    pub gene_symbol: String,
    pub locus_id: String,
    pub member_id: String,
    pub transcript_id: String,
    pub single_candidate_reads: u64,
    /// Overlapping partial seed support; do not sum this column as read totals.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub co_nominated_reads: Option<u64>,
    /// Historical v1 name only; does not assert compatibility with the entire read.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ambiguous_compatible_reads: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelSummaryAmbiguousSet {
    pub sample_id: String,
    pub member_ids: Vec<String>,
    pub reads: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelSummaryNominationSet {
    pub sample_id: String,
    pub class: String,
    pub member_ids: Vec<String>,
    pub reads: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelSampleCoverage {
    pub sample_id: String,
    pub barcode: String,
    pub priming: String,
    pub requested_files: usize,
    pub verified_files: usize,
    pub status: String,
    /// None means no verified unit, not zero observed reads.
    pub reads: Option<u64>,
    pub inputs_not_completed: Vec<RnaReadPanelInput>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelArtifactVerification {
    pub schema: String,
    pub summary_sha256: String,
    pub artifacts: Vec<RnaReadPanelArtifact>,
    pub verified_bytes: u64,
    pub original_inputs_reverified: bool,
    pub scope: String,
}

/// Read support, not molecular abundance, exclusive isoform assignment or expression.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelSummary {
    pub schema: String,
    pub content_sha256: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub submitted_request_sha256: Option<String>,
    pub panel_sha256: String,
    pub policy: RnaReadPanelPolicy,
    /// Complete only for the explicitly supplied reports; never whole-run or transcriptome coverage.
    pub completeness: String,
    pub source_runs: Vec<RnaReadPanelSummaryRun>,
    pub source_files: Vec<RnaReadPanelSummaryFile>,
    pub samples: Vec<RnaReadPanelSummarySample>,
    pub transcripts: Vec<RnaReadPanelSummaryTranscript>,
    /// Historical v1 rows have no recoverable class/member association.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ambiguous_sets: Option<Vec<RnaReadPanelSummaryAmbiguousSet>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub co_nomination_sets: Option<Vec<RnaReadPanelSummaryNominationSet>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_coverage: Option<Vec<RnaReadPanelSampleCoverage>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub class_member_counts_complete: Option<bool>,
    pub original_inputs_reverified: bool,
    pub nomination_sidecars_reverified: bool,
    pub artifacts: Vec<RnaReadPanelArtifact>,
    pub nonclaims: Vec<String>,
}

/// Portable binding for reopening an exact file-only result, not an authenticated signature.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelResultBinding {
    pub schema: String,
    pub kind: RnaReadPanelResultKind,
    pub submitted_request_sha256: String,
    pub effective_request_sha256: String,
    pub panel_sha256: String,
    pub report_content_sha256: String,
    pub report_path: String,
    /// Screen: completed/cancelled/control_error. Summary: complete/partial source coverage.
    pub outcome: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadPanelResultKind {
    Screen,
    Summary,
}
