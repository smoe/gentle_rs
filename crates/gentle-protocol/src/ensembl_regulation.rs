//! Portable contracts for optional Ensembl Regulation annotations.
//!
//! Large provider snapshots remain external to project state. GENtle records
//! their release, assembly, content identity, and derived interval evidence so
//! every adapter can distinguish annotation from activity and primary signal.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const ENSEMBL_REGULATION_SOURCE_CATALOG_SCHEMA: &str =
    "gentle.ensembl_regulation_source_catalog.v1";
pub const ENSEMBL_REGULATION_INTERVAL_INDEX_SCHEMA: &str =
    "gentle.ensembl_regulation_interval_index.v1";
pub const ENSEMBL_REGULATION_OVERLAP_SCHEMA: &str = "gentle.ensembl_regulation_overlap.v1";
pub const ENSEMBL_REGULATION_MATERIALIZATION_SCHEMA: &str =
    "gentle.ensembl_regulation_materialization.v1";
pub const ENSEMBL_REGULATION_INSTALL_REPORT_SCHEMA: &str =
    "gentle.ensembl_regulation_install_report.v1";

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EnsemblRegulationSourceDescriptor {
    pub source_id: String,
    pub provider: String,
    pub annotation_release: String,
    pub annotation_api_version: String,
    pub pipeline_version: String,
    pub gene_annotation_source: String,
    pub gene_annotation_release: String,
    pub species_scientific_name: String,
    pub taxon_id: u32,
    pub assembly_name: String,
    pub assembly_aliases: Vec<String>,
    pub assembly_accession: String,
    pub feature_types: Vec<String>,
    pub annotation_api_url: String,
    pub regulatory_activity_url: Option<String>,
    pub data_access_url: String,
    pub primary_analysis_url: String,
    pub coordinate_system: String,
    pub scope_note: String,
    pub activity_note: String,
    pub primary_signal_note: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EnsemblRegulationSourceCatalog {
    pub schema: String,
    pub sources: Vec<EnsemblRegulationSourceDescriptor>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EnsemblRegulationBinCheckpoint {
    pub bin_index: u64,
    pub scan_start_byte: u64,
    pub scan_start_line_number: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EnsemblRegulationChromosomeIndex {
    pub chromosome: String,
    pub row_count: u64,
    pub first_row_byte: u64,
    pub end_byte_exclusive: u64,
    pub bins: Vec<EnsemblRegulationBinCheckpoint>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EnsemblRegulationIntervalIndex {
    pub schema: String,
    pub source: EnsemblRegulationSourceDescriptor,
    pub intervals_path_hint: String,
    pub intervals_file_name: String,
    pub intervals_sha256: String,
    pub intervals_size_bytes: u64,
    pub row_count: u64,
    pub chromosome_count: usize,
    pub feature_type_counts: BTreeMap<String, u64>,
    pub bin_size_bp: u64,
    pub chromosomes: BTreeMap<String, EnsemblRegulationChromosomeIndex>,
    pub chromosome_aliases: BTreeMap<String, String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EnsemblRegulationInterval {
    pub chromosome: String,
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
    pub feature_id: String,
    pub feature_type: String,
    pub strand: Option<char>,
    pub extended_start_0based: Option<u64>,
    pub extended_end_0based_exclusive: Option<u64>,
    pub associated_gene_ids: Vec<String>,
    pub associated_gene_names: Vec<String>,
    pub source_line_number: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EnsemblRegulationOverlapRow {
    pub interval: EnsemblRegulationInterval,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub genomic_start_0based: u64,
    pub genomic_end_0based_exclusive: u64,
    pub overlap_bp: usize,
    pub clipped: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EnsemblRegulationTypeSummary {
    pub feature_type: String,
    pub matched_count: usize,
    pub overlap_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct EnsemblRegulationGenomeAnchor {
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
pub struct EnsemblRegulationOverlapReport {
    pub schema: String,
    pub report_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub seq_id: String,
    pub index_path: String,
    pub resolved_intervals_path: String,
    pub index_sha256: String,
    pub source: EnsemblRegulationSourceDescriptor,
    pub intervals_sha256: String,
    pub content_identity_verified: bool,
    pub assembly_match_status: String,
    pub genome_anchor: Option<EnsemblRegulationGenomeAnchor>,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub query_length_bp: usize,
    pub requested_feature_types: Vec<String>,
    pub max_rows: usize,
    pub matched_feature_count: usize,
    pub returned_feature_count: usize,
    pub truncated: bool,
    pub type_summaries: Vec<EnsemblRegulationTypeSummary>,
    pub rows: Vec<EnsemblRegulationOverlapRow>,
    pub evidence_statement: String,
    pub non_claims: Vec<String>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EnsemblRegulationMaterializationReport {
    pub schema: String,
    pub report_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub seq_id: String,
    pub index_path: String,
    pub resolved_intervals_path: String,
    pub source: EnsemblRegulationSourceDescriptor,
    pub intervals_sha256: String,
    pub requested_feature_types: Vec<String>,
    pub matched_feature_count: usize,
    pub added_feature_count: usize,
    pub skipped_existing_count: usize,
    pub removed_existing_count: usize,
    pub feature_ids: Vec<usize>,
    pub evidence_statement: String,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EnsemblRegulationInstallReport {
    pub schema: String,
    pub report_id: String,
    pub source: EnsemblRegulationSourceDescriptor,
    pub input: String,
    pub intervals_output: String,
    pub index_output: String,
    pub intervals_sha256: String,
    pub row_count: u64,
    pub feature_type_counts: BTreeMap<String, u64>,
    pub page_count: usize,
    pub fetched_release: String,
    pub warnings: Vec<String>,
}
