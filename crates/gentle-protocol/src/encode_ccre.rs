//! Portable contracts for local ENCODE SCREEN candidate cis-regulatory elements.
//!
//! SCREEN Registry files are optional external evidence. GENtle keeps the
//! original BED payload outside project state and uses a compact, content-bound
//! index for explicit overlap queries and feature materialization.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const ENCODE_CCRE_SOURCE_CATALOG_SCHEMA: &str = "gentle.encode_ccre_source_catalog.v1";
pub const ENCODE_CCRE_INTERVAL_INDEX_SCHEMA: &str = "gentle.encode_ccre_interval_index.v1";
pub const ENCODE_CCRE_OVERLAP_SCHEMA: &str = "gentle.encode_ccre_overlap.v1";
pub const ENCODE_CCRE_MATERIALIZATION_SCHEMA: &str = "gentle.encode_ccre_materialization.v1";
pub const ENCODE_CCRE_INSTALL_REPORT_SCHEMA: &str = "gentle.encode_ccre_install_report.v1";

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EncodeCcreSourceDescriptor {
    pub source_id: String,
    pub provider: String,
    pub registry_version: String,
    pub species_scientific_name: String,
    pub taxon_id: u32,
    pub assembly_name: String,
    pub assembly_aliases: Vec<String>,
    pub subset_id: String,
    pub subset_label: String,
    pub primary_classes: Vec<String>,
    pub dhs_accession_prefix: String,
    pub ccre_accession_prefix: String,
    pub source_url: String,
    pub download_page_url: String,
    pub publication_url: String,
    pub coordinate_system: String,
    pub field_order: Vec<String>,
    pub scope_note: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EncodeCcreSourceCatalog {
    pub schema: String,
    pub registry_version: String,
    pub sources: Vec<EncodeCcreSourceDescriptor>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EncodeCcreBinCheckpoint {
    pub bin_index: u64,
    pub scan_start_byte: u64,
    pub scan_start_line_number: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EncodeCcreChromosomeIndex {
    pub chromosome: String,
    pub row_count: u64,
    pub first_row_byte: u64,
    pub end_byte_exclusive: u64,
    pub bins: Vec<EncodeCcreBinCheckpoint>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EncodeCcreIntervalIndex {
    pub schema: String,
    pub source: EncodeCcreSourceDescriptor,
    pub bed_path_hint: String,
    pub bed_file_name: String,
    pub bed_sha256: String,
    pub bed_size_bytes: u64,
    pub row_count: u64,
    pub chromosome_count: usize,
    pub class_counts: BTreeMap<String, u64>,
    pub bin_size_bp: u64,
    pub chromosomes: BTreeMap<String, EncodeCcreChromosomeIndex>,
    pub chromosome_aliases: BTreeMap<String, String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EncodeCcreInterval {
    pub chromosome: String,
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
    pub dhs_accession: String,
    pub ccre_accession: String,
    pub ccre_class: String,
    pub source_line_number: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EncodeCcreOverlapRow {
    pub interval: EncodeCcreInterval,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub genomic_start_0based: u64,
    pub genomic_end_0based_exclusive: u64,
    pub overlap_bp: usize,
    pub clipped: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EncodeCcreClassSummary {
    pub ccre_class: String,
    pub matched_count: usize,
    pub overlap_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EncodeCcreGenomeAnchor {
    pub seq_id: String,
    pub genome_id: String,
    pub chromosome: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: Option<char>,
    pub anchor_verified: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EncodeCcreOverlapReport {
    pub schema: String,
    pub report_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub seq_id: String,
    pub index_path: String,
    pub resolved_bed_path: String,
    pub index_sha256: String,
    pub source: EncodeCcreSourceDescriptor,
    pub bed_sha256: String,
    pub content_identity_verified: bool,
    pub assembly_match_status: String,
    pub genome_anchor: Option<EncodeCcreGenomeAnchor>,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub query_length_bp: usize,
    pub requested_classes: Vec<String>,
    pub max_rows: usize,
    pub matched_ccre_count: usize,
    pub returned_ccre_count: usize,
    pub truncated: bool,
    pub class_summaries: Vec<EncodeCcreClassSummary>,
    pub rows: Vec<EncodeCcreOverlapRow>,
    pub evidence_statement: String,
    pub non_claims: Vec<String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EncodeCcreMaterializationReport {
    pub schema: String,
    pub report_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub seq_id: String,
    pub index_path: String,
    pub resolved_bed_path: String,
    pub source: EncodeCcreSourceDescriptor,
    pub bed_sha256: String,
    pub requested_classes: Vec<String>,
    pub matched_ccre_count: usize,
    pub added_feature_count: usize,
    pub skipped_existing_count: usize,
    pub removed_existing_count: usize,
    pub feature_ids: Vec<usize>,
    pub evidence_statement: String,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EncodeCcreInstallReport {
    pub schema: String,
    pub report_id: String,
    pub source: EncodeCcreSourceDescriptor,
    pub input: String,
    pub bed_output: String,
    pub index_output: String,
    pub bed_sha256: String,
    pub row_count: u64,
    pub class_counts: BTreeMap<String, u64>,
    pub warnings: Vec<String>,
}
