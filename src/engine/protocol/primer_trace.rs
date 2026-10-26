//! Primer-anchored read observations; local concordance is not transcript assignment.

use super::{RnaReadPanelArtifact, RnaReadPanelControlError, RnaReadPanelInput};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const RNA_PRIMER_PANEL_SCHEMA: &str = "gentle.rna_primer_panel.v1";
pub const RNA_PRIMER_TRACE_REQUEST_SCHEMA: &str = "gentle.rna_primer_trace_request.v1";
pub const RNA_PRIMER_TRACE_SCHEMA: &str = "gentle.rna_primer_trace_report.v1";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PrimerCoreOrientation {
    #[serde(rename = "sense_5prime")]
    Sense5prime,
    #[serde(rename = "antisense_3prime")]
    Antisense3prime,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PrimerReadStrand {
    Forward,
    ReverseComplement,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaPrimerComponent {
    pub component_id: String,
    pub role: String,
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaPrimerTarget {
    pub locus_id: String,
    pub transcript_feature_id: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Caller-supplied ordered oligo, never reconstructed from a gene name or proposed core.
pub struct RnaLibraryPrimer {
    pub primer_id: String,
    pub assay_label: String,
    pub sequence: String,
    pub components: Vec<RnaPrimerComponent>,
    pub core_orientation: PrimerCoreOrientation,
    pub expected_targets: Vec<RnaPrimerTarget>,
    pub variant: String,
    pub source: String,
    pub source_sha256: Option<String>,
}

/// Uncalibrated matching policy; every field participates in panel identity.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RnaPrimerTracePolicy {
    pub max_full_edits: u8,
    pub max_core_edits_for_attribution: u8,
    pub end_window_nt: usize,
    pub check_handle_followed_by_poly_t: bool,
    pub poly_t_min_run: usize,
    pub poly_t_max_mismatches: usize,
    pub min_payload_nt: usize,
    pub concordance_probe_nt: usize,
    pub max_probe_edits: usize,
    pub concordance_seed_nt: usize,
}
impl Default for RnaPrimerTracePolicy {
    fn default() -> Self {
        Self {
            max_full_edits: 2,
            max_core_edits_for_attribution: 1,
            end_window_nt: 256,
            check_handle_followed_by_poly_t: true,
            poly_t_min_run: 12,
            poly_t_max_mismatches: 1,
            min_payload_nt: 40,
            concordance_probe_nt: 30,
            max_probe_edits: 2,
            concordance_seed_nt: 8,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Semantic content binding includes ordered definitions, provenance and the complete policy.
pub struct RnaPrimerPanel {
    pub schema: String,
    pub content_sha256: String,
    pub primers: Vec<RnaLibraryPrimer>,
    pub policy: RnaPrimerTracePolicy,
}

/// Exactly one file binding or inline panel; a file binding always pins its semantic digest.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum RnaPrimerPanelBinding {
    Inline {
        panel: Box<RnaPrimerPanel>,
    },
    File {
        path: String,
        expected_sha256: String,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaPrimerReferenceBinding {
    pub path: String,
    pub expected_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RnaPrimerTraceBudgets {
    pub max_read_nt: usize,
    pub max_hits_per_read: usize,
    pub max_reference_sites: usize,
    pub max_index_postings: usize,
    pub max_scan_cells: usize,
    pub max_sidecar_bytes: u64,
}
impl Default for RnaPrimerTraceBudgets {
    fn default() -> Self {
        Self {
            max_read_nt: 200_000,
            max_hits_per_read: 256,
            max_reference_sites: 4096,
            max_index_postings: 2_000_000,
            max_scan_cells: 100_000_000,
            max_sidecar_bytes: 1_000_000_000,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Explicit completed inputs and a fresh external result directory; no project mutation or discovery.
pub struct TracePrimerReadsRequest {
    pub schema: String,
    pub primer_panel: RnaPrimerPanelBinding,
    pub reference_panel: Option<RnaPrimerReferenceBinding>,
    pub files: Vec<RnaReadPanelInput>,
    pub md5_manifest_path: Option<String>,
    pub output_dir: String,
    #[serde(default)]
    pub export_payload_fasta: bool,
    #[serde(default)]
    pub budgets: RnaPrimerTraceBudgets,
    #[serde(default)]
    pub cancel_path: Option<String>,
}

/// Coordinates are always original-read half-open intervals, never trimmed read coordinates.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaPrimerHit {
    pub primer_id: Option<String>,
    pub anchor_class: String,
    pub strand: PrimerReadStrand,
    pub start: usize,
    pub end: usize,
    pub edits: usize,
    pub core_start: Option<usize>,
    pub core_end: Option<usize>,
    pub core_edits: Option<usize>,
    pub terminal: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaPrimerTemplateMatch {
    pub member_id: String,
    pub template_start: usize,
    pub template_end: usize,
    pub expected_target: bool,
    pub core_and_probe_contiguous: bool,
    pub payload_direction_matches_declared: bool,
    pub probe_edits: usize,
    /// Exact matching prefix before the first unequal base, not an aligned-read length.
    pub concordant_nt: usize,
    pub core_3prime_mismatches: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Primary-anchor interpretation plus retained alternatives; contains no read sequence.
pub struct RnaPrimerTraceRow {
    pub record_index: u64,
    pub read_id: String,
    pub read_length: usize,
    pub anchor_class: String,
    pub anchor_position: Option<String>,
    pub primer_multiplicity: Option<String>,
    pub locus_relation: String,
    pub evidence_tier: String,
    pub anchor: Option<RnaPrimerHit>,
    pub hits: Vec<RnaPrimerHit>,
    pub leader: Option<(usize, usize)>,
    pub payload: Option<(usize, usize)>,
    pub leader_length: Option<usize>,
    pub payload_length: Option<usize>,
    pub payload_transcript_sense: Option<bool>,
    pub payload_mean_phred: Option<f64>,
    pub payload_fraction_below_q10: Option<f64>,
    pub template_matches: Vec<RnaPrimerTemplateMatch>,
    pub artifact_flags: Vec<String>,
    pub reasons: Vec<String>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaPrimerTraceCounts {
    pub reads: u64,
    pub anchor_classes: BTreeMap<String, u64>,
    pub locus_relations: BTreeMap<String, u64>,
    /// All retained full hits, not just the primary anchor; a read may contribute more than once.
    pub full_hits_by_primer_class_strand: BTreeMap<String, u64>,
    pub anchor_start_histogram: BTreeMap<String, u64>,
    pub remaining_tail_histogram: BTreeMap<String, u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaPrimerTraceFile {
    pub file_id: String,
    pub input: RnaReadPanelInput,
    pub observed_sha256: String,
    pub observed_md5: String,
    pub compressed_bytes: u64,
    pub counts: RnaPrimerTraceCounts,
    pub artifacts: Vec<RnaReadPanelArtifact>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaPrimerTraceSummary {
    pub schema: String,
    /// Exact positions 0..256, then fixed right-closed doubling buckets through 2,000,000.
    pub histogram_policy: String,
    pub by_sample_barcode: BTreeMap<String, RnaPrimerTraceCounts>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Only EOF-verified file units contribute counts; missing units are explicit, never numerical zero.
pub struct RnaPrimerTraceReport {
    pub schema: String,
    pub content_sha256: String,
    pub algorithm: String,
    pub submitted_request_sha256: String,
    pub status: String,
    pub request: TracePrimerReadsRequest,
    pub primer_panel: RnaPrimerPanel,
    pub reference_panel_sha256: Option<String>,
    pub artifacts: Vec<RnaReadPanelArtifact>,
    pub files_completed: Vec<RnaPrimerTraceFile>,
    pub inputs_not_completed: Vec<RnaReadPanelInput>,
    pub control_error: Option<RnaReadPanelControlError>,
    pub summary: RnaPrimerTraceSummary,
    pub warnings: Vec<String>,
    pub nonclaims: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RnaPrimerTraceInspection {
    pub report: Box<RnaPrimerTraceReport>,
    pub rows: Vec<RnaPrimerTraceRow>,
    pub artifact_content_verified: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Exact content/location/outcome receipt, not authenticated provenance or a worker-recovery handle.
pub struct RnaPrimerTraceResultBinding {
    pub schema: String,
    pub submitted_request_sha256: String,
    pub effective_request_sha256: String,
    pub primer_panel_sha256: String,
    pub reference_panel_sha256: Option<String>,
    pub report_content_sha256: String,
    pub report_path: String,
    pub outcome: String,
}
