//! Preview-bound derivation of exact annotated transcript starts into project sequences.

use serde::{Deserialize, Serialize};

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
    pub request: TssInventoryRequest,
    pub source_snapshot_sha256: String,
    /// Includes source, geometry, transcript membership and output namespace.
    pub approval_sha256: String,
    pub rows: Vec<TssInventoryRow>,
    pub warnings: Vec<String>,
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
    /// Existing map operations consume this explicit sequence subject.
    pub subject: CollectionSubjectRef,
}
