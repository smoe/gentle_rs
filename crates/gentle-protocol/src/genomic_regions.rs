//! Portable, assembly-bound genomic regions of interest.
//!
//! The contracts in this module keep canonical genomic geometry separate from
//! sequence-local projections and from the evidence explaining why a region
//! was selected. JSON is lossless; BED is intentionally accompanied by a
//! content-bound manifest when it is used for round-tripping.

use crate::{
    CutRunRegulatorySupportReport, EnsemblRegulationOverlapReport,
    GeneLocusEnsemblRegulationFeatureRow, GeneLocusEnsemblRegulationSourceBinding,
};
use serde::{Deserialize, Serialize};

pub const GENOMIC_REGION_OF_INTEREST_SCHEMA: &str = "gentle.genomic_region_of_interest.v1";
pub const GENOMIC_REGION_SET_SCHEMA: &str = "gentle.genomic_region_set.v1";
pub const GENOMIC_REGION_STORE_SCHEMA: &str = "gentle.genomic_region_store.v1";
pub const GENOMIC_REGION_BED_MANIFEST_SCHEMA: &str = "gentle.genomic_region_bed_manifest.v1";
pub const GENOMIC_REGION_OPERATION_REPORT_SCHEMA: &str =
    "gentle.genomic_region_operation_report.v1";

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionCoordinateConvention {
    #[default]
    ZeroBasedHalfOpen,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionStrand {
    Plus,
    Minus,
    #[default]
    Unstranded,
}

impl GenomicRegionStrand {
    pub fn bed_value(self) -> &'static str {
        match self {
            Self::Plus => "+",
            Self::Minus => "-",
            Self::Unstranded => ".",
        }
    }

    pub fn human_value(self) -> &'static str {
        match self {
            Self::Plus => "+",
            Self::Minus => "-",
            Self::Unstranded => "unstranded",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionReference {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub species_scientific_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub taxon_id: Option<u32>,
    pub assembly_name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assembly_accession: Option<String>,
    pub contig_name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub contig_accession: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub contig_aliases: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionInterval {
    pub reference: GenomicRegionReference,
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
    pub strand: GenomicRegionStrand,
    pub coordinate_convention: GenomicRegionCoordinateConvention,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionEvidenceAvailability {
    #[default]
    Available,
    Unavailable,
    Stale,
    Unverified,
}

impl GenomicRegionEvidenceAvailability {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Available => "available",
            Self::Unavailable => "unavailable",
            Self::Stale => "stale",
            Self::Unverified => "unverified",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionPurpose {
    CandidateCisRegulatoryRegion,
    OccupancyRegion,
    PromoterRegion,
    ReporterCandidate,
    #[default]
    Other,
}

impl GenomicRegionPurpose {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CandidateCisRegulatoryRegion => "candidate_cis_regulatory_region",
            Self::OccupancyRegion => "occupancy_region",
            Self::PromoterRegion => "promoter_region",
            Self::ReporterCandidate => "reporter_candidate",
            Self::Other => "other",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionSelectionMethod {
    #[default]
    ManualSpan,
    ExistingInterval,
    CutrunSupportWindow,
    EnsemblRegulatoryFeature,
    ProviderAnnotation,
    Imported,
    Derived,
}

impl GenomicRegionSelectionMethod {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ManualSpan => "manual_span",
            Self::ExistingInterval => "existing_interval",
            Self::CutrunSupportWindow => "cutrun_support_window",
            Self::EnsemblRegulatoryFeature => "ensembl_regulatory_feature",
            Self::ProviderAnnotation => "provider_annotation",
            Self::Imported => "imported",
            Self::Derived => "derived",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionLocalProjectionStatus {
    #[default]
    Current,
    SequenceUnavailable,
    SequenceDigestMismatch,
    AnchorMismatch,
    OutsideSequence,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionLocalProjection {
    pub seq_id: String,
    pub sequence_sha256: String,
    pub source_genome_id: String,
    pub anchor_start_1based: u64,
    pub anchor_end_1based: u64,
    pub anchor_strand: GenomicRegionStrand,
    pub local_start_0based: u64,
    pub local_end_0based_exclusive: u64,
    pub local_strand: GenomicRegionStrand,
    pub status: GenomicRegionLocalProjectionStatus,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionEvidenceReference {
    pub evidence_id: String,
    pub source_kind: String,
    pub source_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_release: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_sha256: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub report_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub feature_or_window_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub feature_type: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub provider_url: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub associated_gene_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub associated_gene_names: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assay: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub factor: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub support_strength: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_signal_value: Option<f64>,
    pub availability: GenomicRegionEvidenceAvailability,
    /// Conservative evidence statement retained verbatim from the producer.
    pub evidence_statement: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub non_claims: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionDerivationMethod {
    Union,
    Intersection,
    #[default]
    Hull,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionParentReference {
    pub region_id: String,
    pub identity_sha256: String,
    pub content_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionDerivation {
    pub method: GenomicRegionDerivationMethod,
    pub parents: Vec<GenomicRegionParentReference>,
    pub rule: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionOfInterest {
    pub schema: String,
    pub region_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    pub interval: GenomicRegionInterval,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub local_projection: Option<GenomicRegionLocalProjection>,
    pub purpose: GenomicRegionPurpose,
    pub selection_method: GenomicRegionSelectionMethod,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub evidence: Vec<GenomicRegionEvidenceReference>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub derivation: Option<GenomicRegionDerivation>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub created_at_unix_ms: Option<u128>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
    pub identity_sha256: String,
    pub content_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionSet {
    pub schema: String,
    pub set_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    pub regions: Vec<GenomicRegionOfInterest>,
    pub content_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionStore {
    pub schema: String,
    pub sets: Vec<GenomicRegionSet>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionCollisionPolicy {
    #[default]
    Reject,
    Replace,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionCreateRequest {
    pub set_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub set_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub region_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    pub interval: GenomicRegionInterval,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub local_projection: Option<GenomicRegionLocalProjection>,
    pub purpose: GenomicRegionPurpose,
    pub selection_method: GenomicRegionSelectionMethod,
    pub evidence: Vec<GenomicRegionEvidenceReference>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub created_at_unix_ms: Option<u128>,
    pub notes: Vec<String>,
    pub collision_policy: GenomicRegionCollisionPolicy,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionEnsemblIntervalKind {
    #[default]
    Core,
    Extended,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "source_kind", rename_all = "snake_case")]
pub enum GenomicRegionCaptureSource {
    SequenceSelection {
        seq_id: String,
        local_start_0based: u64,
        local_end_0based_exclusive: u64,
        #[serde(default)]
        strand: GenomicRegionStrand,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        reference_override: Option<GenomicRegionReference>,
    },
    SequenceFeature {
        seq_id: String,
        feature_id: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        reference_override: Option<GenomicRegionReference>,
    },
    CutrunSupportWindow {
        report: Box<CutRunRegulatorySupportReport>,
        window_id: String,
        reference: GenomicRegionReference,
    },
    EnsemblRegulatoryFeature {
        report: Box<EnsemblRegulationOverlapReport>,
        feature_id: String,
        #[serde(default)]
        interval_kind: GenomicRegionEnsemblIntervalKind,
    },
    /// Capture an already-normalized row carried by a gene-locus report without
    /// requiring the original provider overlap report to remain in memory.
    GeneLocusEnsemblRegulatoryFeature {
        row: GeneLocusEnsemblRegulationFeatureRow,
        source_binding: GeneLocusEnsemblRegulationSourceBinding,
        reference: GenomicRegionReference,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        seq_id: Option<String>,
        #[serde(default)]
        interval_kind: GenomicRegionEnsemblIntervalKind,
        evidence_statement: String,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        non_claims: Vec<String>,
    },
    ProviderAnnotation {
        interval: GenomicRegionInterval,
        evidence: GenomicRegionEvidenceReference,
    },
}

impl Default for GenomicRegionCaptureSource {
    fn default() -> Self {
        Self::SequenceSelection {
            seq_id: String::new(),
            local_start_0based: 0,
            local_end_0based_exclusive: 0,
            strand: GenomicRegionStrand::Unstranded,
            reference_override: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GenomicRegionCaptureRequest {
    pub set_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub set_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub region_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    pub purpose: GenomicRegionPurpose,
    pub source: GenomicRegionCaptureSource,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub created_at_unix_ms: Option<u128>,
    pub notes: Vec<String>,
    pub collision_policy: GenomicRegionCollisionPolicy,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionUpdateRequest {
    pub set_id: String,
    pub region_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionDeriveRequest {
    pub set_id: String,
    pub parent_region_ids: Vec<String>,
    pub method: GenomicRegionDerivationMethod,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub region_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    pub purpose: GenomicRegionPurpose,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub created_at_unix_ms: Option<u128>,
    pub notes: Vec<String>,
    pub collision_policy: GenomicRegionCollisionPolicy,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionListRequest {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub set_id: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionInspectRequest {
    pub set_id: String,
    pub region_id: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionImportFormat {
    #[default]
    Json,
    Bed,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionImportRequest {
    pub path: String,
    pub format: GenomicRegionImportFormat,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub manifest_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bare_bed_reference: Option<GenomicRegionReference>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub set_id_override: Option<String>,
    pub collision_policy: GenomicRegionCollisionPolicy,
    #[serde(default = "default_genomic_region_import_max_bytes")]
    pub max_bytes: u64,
    #[serde(default = "default_genomic_region_import_max_rows")]
    pub max_rows: usize,
}

impl Default for GenomicRegionImportRequest {
    fn default() -> Self {
        Self {
            path: String::new(),
            format: GenomicRegionImportFormat::Json,
            manifest_path: None,
            bare_bed_reference: None,
            set_id_override: None,
            collision_policy: GenomicRegionCollisionPolicy::Reject,
            max_bytes: default_genomic_region_import_max_bytes(),
            max_rows: default_genomic_region_import_max_rows(),
        }
    }
}

pub const fn default_genomic_region_import_max_bytes() -> u64 {
    10 * 1024 * 1024
}

pub const fn default_genomic_region_import_max_rows() -> usize {
    100_000
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionExportRequest {
    pub set_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub json_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bed_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub manifest_path: Option<String>,
    pub include_local_paths: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionBedManifest {
    pub schema: String,
    pub bed_sha256: String,
    pub bed_file_name: String,
    pub coordinate_convention: String,
    pub columns: Vec<String>,
    pub region_set: GenomicRegionSet,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionSetSummary {
    pub set_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub region_count: usize,
    pub content_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionOperationReport {
    pub schema: String,
    pub action: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub set: Option<GenomicRegionSet>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub region: Option<GenomicRegionOfInterest>,
    pub sets: Vec<GenomicRegionSetSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub human_copy: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bed_row: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub canonical_roi_json: Option<String>,
    pub written_artifacts: Vec<String>,
    pub warnings: Vec<String>,
}
