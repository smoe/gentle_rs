//! Portable contracts for querying precomputed whole-genome motif evidence.
//!
//! The provider is deliberately optional. A request may resolve to an
//! unavailable report when no package or DuckDB executable is configured;
//! existing local motif scoring and every other GENtle workflow remain usable.

use serde::{Deserialize, Serialize};

pub const GENOMIC_MOTIF_EVIDENCE_SCHEMA: &str = "gentle.genomic_motif_evidence.v1";
pub const DEFAULT_GENOMIC_MOTIF_EVIDENCE_MAX_ROWS: usize = 10_000;
pub const DEFAULT_GENOMIC_MOTIF_EVIDENCE_MAX_PAYLOAD_FILES: usize = 256;
pub const DEFAULT_GENOMIC_MOTIF_EVIDENCE_TIMEOUT_SECONDS: u64 = 30;
pub const MAX_GENOMIC_MOTIF_EVIDENCE_QUERY_MOTIFS: usize = 64;

fn default_max_rows() -> usize {
    DEFAULT_GENOMIC_MOTIF_EVIDENCE_MAX_ROWS
}

fn default_max_payload_files() -> usize {
    DEFAULT_GENOMIC_MOTIF_EVIDENCE_MAX_PAYLOAD_FILES
}

fn default_timeout_seconds() -> u64 {
    DEFAULT_GENOMIC_MOTIF_EVIDENCE_TIMEOUT_SECONDS
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicMotifEvidenceAvailability {
    Available,
    #[default]
    PackageNotConfigured,
    PackageMissing,
    DuckdbUnavailable,
    InvalidPackage,
    IncompatiblePackage,
    QueryFailed,
}

impl GenomicMotifEvidenceAvailability {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Available => "available",
            Self::PackageNotConfigured => "package_not_configured",
            Self::PackageMissing => "package_missing",
            Self::DuckdbUnavailable => "duckdb_unavailable",
            Self::InvalidPackage => "invalid_package",
            Self::IncompatiblePackage => "incompatible_package",
            Self::QueryFailed => "query_failed",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicMotifEvidenceCompatibilityStatus {
    #[default]
    NotAssessed,
    ContigGeometryMatchedOnly,
    ContigGeometryMismatch,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicMotifEvidenceCoverageStatus {
    #[default]
    NotAssessed,
    CompleteForPackageRetention,
    CompleteForRequestedThreshold,
    IncompleteBelowStorageFloor,
    DensityLimitedAtSource,
    MotifNotInPackage,
    TruncatedAtMaxRows,
}

impl GenomicMotifEvidenceCoverageStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NotAssessed => "not_assessed",
            Self::CompleteForPackageRetention => "complete_for_package_retention",
            Self::CompleteForRequestedThreshold => "complete_for_requested_threshold",
            Self::IncompleteBelowStorageFloor => "incomplete_below_storage_floor",
            Self::DensityLimitedAtSource => "density_limited_at_source",
            Self::MotifNotInPackage => "motif_not_in_package",
            Self::TruncatedAtMaxRows => "truncated_at_max_rows",
        }
    }
}

impl GenomicMotifEvidenceCompatibilityStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NotAssessed => "not_assessed",
            Self::ContigGeometryMatchedOnly => "contig_geometry_matched_only",
            Self::ContigGeometryMismatch => "contig_geometry_mismatch",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicMotifEvidenceInterval {
    pub interval_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub chromosome: String,
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "target_kind", rename_all = "snake_case")]
pub enum GenomicMotifEvidenceTarget {
    AnchoredSequence {
        seq_id: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        span_start_0based: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        span_end_0based_exclusive: Option<usize>,
    },
    GenomicIntervals {
        intervals: Vec<GenomicMotifEvidenceInterval>,
    },
}

impl Default for GenomicMotifEvidenceTarget {
    fn default() -> Self {
        Self::GenomicIntervals { intervals: vec![] }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct GenomicMotifEvidenceRequest {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub package_root: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub database_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub duckdb_executable: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expected_genome_id: Option<String>,
    pub target: GenomicMotifEvidenceTarget,
    pub motif_ids: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub minimum_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub minimum_pwm_relative_score: Option<f64>,
    #[serde(default = "default_max_rows")]
    pub max_rows: usize,
    #[serde(default = "default_max_payload_files")]
    pub max_payload_files: usize,
    #[serde(default = "default_timeout_seconds")]
    pub timeout_seconds: u64,
}

impl Default for GenomicMotifEvidenceRequest {
    fn default() -> Self {
        Self {
            package_root: None,
            database_path: None,
            duckdb_executable: None,
            expected_genome_id: None,
            target: GenomicMotifEvidenceTarget::default(),
            motif_ids: vec![],
            minimum_score: None,
            minimum_pwm_relative_score: None,
            max_rows: DEFAULT_GENOMIC_MOTIF_EVIDENCE_MAX_ROWS,
            max_payload_files: DEFAULT_GENOMIC_MOTIF_EVIDENCE_MAX_PAYLOAD_FILES,
            timeout_seconds: DEFAULT_GENOMIC_MOTIF_EVIDENCE_TIMEOUT_SECONDS,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicMotifEvidenceProviderProvenance {
    pub provider_kind: String,
    pub package_root: String,
    pub manifest_path: String,
    pub manifest_sha256: String,
    pub declared_content_fingerprint_sha256: String,
    pub manifest_schema_version: u64,
    pub database_path: String,
    pub duckdb_executable: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub duckdb_version: Option<String>,
    pub run_id: String,
    pub genome_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assembly_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assembly_accession: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ensembl_release: Option<String>,
    pub motif_set_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub jaspar_version: Option<String>,
    pub score_mode: String,
    pub pseudocount: f64,
    pub pseudocount_scheme: String,
    pub background_model_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub minimum_pwm_relative_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub maximum_pwm_relative_score: Option<f64>,
    pub coordinate_mode: String,
    pub n_policy: String,
    pub matched_sequence_policy: String,
    pub selected_payloads: Vec<GenomicMotifEvidencePayloadProvenance>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicMotifEvidencePayloadProvenance {
    pub task_id: String,
    pub output_relative_path: String,
    pub declared_sha256: String,
    pub emitted_hits: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicMotifEvidenceMotifCoverage {
    pub motif_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub motif_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub threshold_set_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub informative_threshold: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_minimum_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_minimum_pwm_relative_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub requested_minimum_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub requested_minimum_pwm_relative_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub density_limited: Option<bool>,
    pub status: GenomicMotifEvidenceCoverageStatus,
    pub returned_hit_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicMotifEvidenceResolvedRegion {
    pub interval_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub requested_chromosome: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub resolved_chromosome: Option<String>,
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_seq_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_start_0based: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_end_0based_exclusive: Option<usize>,
    pub compatibility_status: GenomicMotifEvidenceCompatibilityStatus,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub package_contig_length_bp: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicMotifEvidenceHit {
    pub interval_id: String,
    pub chromosome: String,
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
    pub motif_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub motif_name: Option<String>,
    pub strand: String,
    pub score: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pwm_relative_score: Option<f64>,
    pub score_mode: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub minimum_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub matched_sequence: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_seq_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_start_0based: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_end_0based_exclusive: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_forward_strand: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicMotifEvidenceReport {
    pub schema: String,
    pub report_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub availability: GenomicMotifEvidenceAvailability,
    pub request: GenomicMotifEvidenceRequest,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub provider: Option<GenomicMotifEvidenceProviderProvenance>,
    pub regions: Vec<GenomicMotifEvidenceResolvedRegion>,
    pub motif_coverage: Vec<GenomicMotifEvidenceMotifCoverage>,
    pub selected_payload_file_count: usize,
    pub selected_payload_emitted_hit_count: u64,
    pub matched_hit_count: usize,
    pub returned_hit_count: usize,
    pub truncated: bool,
    pub query_complete: bool,
    pub hits: Vec<GenomicMotifEvidenceHit>,
    pub warnings: Vec<String>,
}
