//! Portable multi-locus seed-screen contracts; compatibility is not isoform identification.

use super::SequenceGenomeAnchorSummary;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const RNA_PANEL_SCHEMA: &str = "gentle.rna_read_panel.v1";
pub const RNA_PANEL_BUILD_SCHEMA: &str = "gentle.rna_read_panel_build_request.v1";
pub const RNA_PANEL_SCREEN_SCHEMA: &str = "gentle.rna_read_panel_screen.v1";
pub const RNA_PANEL_SCREEN_REQUEST_SCHEMA: &str = "gentle.rna_read_panel_screen_request.v1";
pub const RNA_PANEL_SUMMARY_REQUEST_SCHEMA: &str = "gentle.rna_read_panel_summary_request.v1";
pub const RNA_PANEL_SUMMARY_SCHEMA: &str = "gentle.rna_read_panel_summary.v1";

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
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Build from loaded loci into a new portable file without mutating those loci.
pub struct BuildRnaReadPanelRequest {
    pub schema: String,
    pub sources: Vec<RnaReadPanelSource>,
    pub output_path: String,
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
    /// Each read increments one member-set class, not every member's exclusive count.
    pub equivalence_classes: BTreeMap<String, u64>,
    pub artifacts: Vec<RnaReadPanelArtifact>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Verified file units and ambiguity counts; cancelled never means a complete run.
pub struct RnaReadPanelScreenReport {
    pub schema: String,
    pub content_sha256: String,
    pub status: String,
    pub request: ScreenRnaReadPanelRequest,
    pub panel_sha256: String,
    pub indexed_postings: usize,
    pub isoform_assignment_performed: bool,
    pub alignment_confirmation_performed: bool,
    pub files_completed: Vec<RnaReadPanelFileReceipt>,
    pub files_not_completed: Vec<String>,
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
    pub ambiguous_reads: u64,
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
    /// The same ambiguous read can support several rows; do not sum this column as read totals.
    pub ambiguous_compatible_reads: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelSummaryAmbiguousSet {
    pub sample_id: String,
    pub member_ids: Vec<String>,
    pub reads: u64,
}

/// Read support, not molecular abundance, exclusive isoform assignment or expression.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaReadPanelSummary {
    pub schema: String,
    pub content_sha256: String,
    pub panel_sha256: String,
    pub policy: RnaReadPanelPolicy,
    /// Complete only for the explicitly supplied reports; never whole-run or transcriptome coverage.
    pub completeness: String,
    pub source_runs: Vec<RnaReadPanelSummaryRun>,
    pub source_files: Vec<RnaReadPanelSummaryFile>,
    pub samples: Vec<RnaReadPanelSummarySample>,
    pub transcripts: Vec<RnaReadPanelSummaryTranscript>,
    pub ambiguous_sets: Vec<RnaReadPanelSummaryAmbiguousSet>,
    pub original_inputs_reverified: bool,
    pub nomination_sidecars_reverified: bool,
    pub artifacts: Vec<RnaReadPanelArtifact>,
    pub nonclaims: Vec<String>,
}
