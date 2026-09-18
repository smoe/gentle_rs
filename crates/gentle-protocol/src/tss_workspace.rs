//! Preview-bound derivation of exact annotated transcript starts into project sequences.

use serde::{Deserialize, Serialize};

/// Registry discovery is not member validation or a biological coverage claim.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TssCollectionListReport {
    pub schema: String,
    pub collections: Vec<TssCollectionListEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TssCollectionListEntry {
    pub collection_id: String,
    pub source_seq_id: Option<String>,
    pub gene_query: Option<String>,
    pub window_count: Option<usize>,
    pub record_status: TssCollectionRecordStatus,
    pub validation_status: TssCollectionValidationStatus,
    pub diagnostic: Option<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TssCollectionRecordStatus {
    Readable,
    Legacy,
    Invalid,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TssCollectionValidationStatus {
    NotChecked,
}

use crate::{collection_subjects::CollectionSubjectRef, genomic_regions::GenomicRegionInterval};

fn upstream() -> usize {
    500
}
fn downstream() -> usize {
    200
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TssInventoryRequest {
    pub seq_id: String,
    pub gene_query: String,
    /// Explicit namespace; existing unrelated outputs are never replaced.
    pub collection_id: String,
    #[serde(default = "upstream")]
    pub upstream_bp: usize,
    #[serde(default = "downstream")]
    pub downstream_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TssWindowAvailability {
    Available,
    MissingFlanks,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TssInventoryRow {
    pub tss_id: String,
    pub output_seq_id: String,
    pub gene_id: Option<String>,
    pub gene_label: Option<String>,
    pub annotation_source: String,
    pub transcript_ids: Vec<String>,
    pub transcript_feature_ids: Vec<usize>,
    pub tss_local_0based: usize,
    pub local_strand: String,
    /// One-base genomic interval, preserving genomic versus local strand.
    pub genomic_tss: GenomicRegionInterval,
    pub window_local_start_0based: Option<usize>,
    pub window_local_end_0based_exclusive: Option<usize>,
    pub availability: TssWindowAvailability,
    pub explanation: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TssInventoryReport {
    pub schema: String,
    /// Empty in legacy payloads whose snapshots included nondeterministic caches.
    #[serde(default)]
    pub snapshot_algorithm: String,
    pub request: TssInventoryRequest,
    pub source_snapshot_sha256: String,
    /// Includes source, geometry, transcript membership and output namespace.
    pub approval_sha256: String,
    pub rows: Vec<TssInventoryRow>,
    #[serde(default)]
    pub excluded_transcripts: Vec<TssExcludedTranscript>,
    /// Locus-level annotations with no gene linkage, not exclusions from the requested gene.
    /// Omitted when empty to preserve serialization of existing bound collections.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub unassigned_transcripts: Vec<TssExcludedTranscript>,
    pub warnings: Vec<String>,
}

/// An annotation that must not be promoted to an exact transcript start.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TssExcludedTranscript {
    pub feature_id: usize,
    pub transcript_id: String,
    pub reason: TssTranscriptExclusionReason,
    pub explanation: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TssTranscriptExclusionReason {
    MissingGeneLink,
    UncertainFivePrimeEnd,
    TruncatedFivePrimeEnd,
    InvalidGenomicBounds,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TssMaterializeRequest {
    pub inventory: TssInventoryRequest,
    pub expected_approval_sha256: String,
    /// Explicit, nonempty selection from the preview. No implicit all-project fallback.
    pub selected_tss_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssCollectionMember {
    pub tss: TssInventoryRow,
    pub sequence_sha256: String,
    pub record_snapshot_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssCollectionReport {
    pub schema: String,
    pub collection_id: String,
    pub inventory: TssInventoryReport,
    pub members: Vec<TssCollectionMember>,
    pub lifting_mode: crate::collection_subjects::CollectionLiftingMode,
    pub collection_membership_fingerprint_sha256: String,
    /// Map operations resolve this typed reference and revalidate member snapshots.
    pub subject: CollectionSubjectRef,
}
