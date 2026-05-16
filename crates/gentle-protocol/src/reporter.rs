//! Portable reporter-catalog and recommender contracts.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Built-in/local reporter catalog schema.
pub const REPORTER_CATALOG_SCHEMA: &str = "gentle.reporter_catalog.v1";
/// Annotated reporter catalog report schema.
pub const REPORTER_CATALOG_REPORT_SCHEMA: &str = "gentle.reporter_catalog_report.v1";
/// Constraint-based reporter recommendation report schema.
pub const REPORTER_RECOMMENDATION_SCHEMA: &str = "gentle.reporter_recommendation.v1";
/// Agent/local-AI reporter corpus export schema.
pub const REPORTER_CORPUS_EXPORT_SCHEMA: &str = "gentle.reporter_corpus_export.v1";

/// Supported corpus export shapes for local AI retrieval/training prep.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum ReporterCorpusExportFormat {
    #[default]
    Json,
    Jsonl,
}

impl ReporterCorpusExportFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Json => "json",
            Self::Jsonl => "jsonl",
        }
    }
}

/// One catalog/source reference backing a reporter record.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterSourceRef {
    pub source_id: String,
    pub accession: String,
    pub url: String,
    pub retrieved_at: String,
    pub license_status: String,
    pub license_note: String,
}

/// Compact spectral and practical-imaging metadata.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterSpectralProfile {
    pub color: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub excitation_nm: Option<u16>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub emission_nm: Option<u16>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub brightness: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub maturation_minutes: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pka: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub oligomerization: Option<String>,
}

/// One curated reporter candidate in a local catalog.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterRecord {
    pub id: String,
    pub name: String,
    pub aliases: Vec<String>,
    pub reporter_class: String,
    pub sequence: String,
    pub sequence_sha1: String,
    pub source_refs: Vec<ReporterSourceRef>,
    pub license_status: String,
    pub provenance_note: String,
    pub colors: Vec<String>,
    pub assay_modes: Vec<String>,
    pub substrate_required: bool,
    pub compatible_hosts: Vec<String>,
    pub fusion_compatibility: Vec<String>,
    pub characterization_confidence: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub catalog_use_count: Option<usize>,
    pub safety_scope: String,
    pub spectral: ReporterSpectralProfile,
    pub notes: Vec<String>,
}

/// Local reporter catalog file.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterCatalog {
    pub schema: String,
    pub curated_at: String,
    pub sources: Vec<ReporterSourceRef>,
    pub records: Vec<ReporterRecord>,
    pub notes: Vec<String>,
}

/// Deterministic annotations computed by GENtle from one reporter sequence.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterComputedAnnotation {
    pub length_bp: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gc_fraction: Option<f64>,
    pub starts_with_atg: bool,
    pub ends_with_stop: bool,
    pub multiple_of_three: bool,
    pub likely_complete_cds: bool,
    pub checksum_ok: bool,
    #[serde(default)]
    pub forbidden_motif_hits: Vec<String>,
}

/// Reporter record plus deterministic GENtle-computed annotations.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterAnnotatedRecord {
    pub record: ReporterRecord,
    pub annotation: ReporterComputedAnnotation,
    pub warnings: Vec<String>,
}

/// Quarantined catalog row that did not satisfy V1 provenance/safety gates.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterQuarantinedRecord {
    pub id: String,
    pub name: String,
    pub reasons: Vec<String>,
}

/// Annotated local catalog report.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterCatalogReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub catalog_path: String,
    pub record_count: usize,
    pub active_record_count: usize,
    pub quarantined_record_count: usize,
    pub records: Vec<ReporterAnnotatedRecord>,
    pub quarantined_records: Vec<ReporterQuarantinedRecord>,
    pub warnings: Vec<String>,
}

/// Optional weights for deterministic soft ranking.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct ReporterPreferenceWeights {
    pub characterization_confidence: f64,
    pub host_match: f64,
    pub assay_match: f64,
    pub spectral_match: f64,
    pub brightness: f64,
    pub short_sequence: f64,
    pub complete_cds: f64,
}

impl Default for ReporterPreferenceWeights {
    fn default() -> Self {
        Self {
            characterization_confidence: 1.0,
            host_match: 1.0,
            assay_match: 1.0,
            spectral_match: 1.0,
            brightness: 1.0,
            short_sequence: 1.0,
            complete_cds: 1.0,
        }
    }
}

/// User/agent constraints for reporter selection.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterConstraints {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub intended_assay: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chassis: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub live_assay: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub desired_color: Option<String>,
    #[serde(default)]
    pub allowed_reporter_classes: Vec<String>,
    #[serde(default)]
    pub available_excitation_nm: Vec<u16>,
    #[serde(default)]
    pub available_emission_nm: Vec<u16>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fusion_mode: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_coding_length_bp: Option<usize>,
    #[serde(default)]
    pub forbidden_motifs: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub substrate_allowed: Option<bool>,
    #[serde(default)]
    pub preference_weights: ReporterPreferenceWeights,
}

/// One accepted ranked reporter candidate.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterRecommendation {
    pub rank: usize,
    pub reporter_id: String,
    pub name: String,
    pub score: f64,
    pub score_components: BTreeMap<String, f64>,
    pub rationale: Vec<String>,
    pub warnings: Vec<String>,
    pub record: ReporterAnnotatedRecord,
}

/// One rejected reporter candidate with explicit hard-constraint reasons.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterRejectedCandidate {
    pub reporter_id: String,
    pub name: String,
    pub reasons: Vec<String>,
}

/// Deterministic reporter recommendation result.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterRecommendationResult {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub catalog_path: String,
    pub constraints: ReporterConstraints,
    pub considered_candidate_count: usize,
    pub recommended_candidate_count: usize,
    pub rejected_candidate_count: usize,
    pub recommendations: Vec<ReporterRecommendation>,
    pub rejected_candidates: Vec<ReporterRejectedCandidate>,
    pub warnings: Vec<String>,
}

/// Exported annotated corpus for retrieval or local training pipelines.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct ReporterCorpusExport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub catalog_path: String,
    pub format: ReporterCorpusExportFormat,
    pub record_count: usize,
    pub records: Vec<ReporterAnnotatedRecord>,
    pub warnings: Vec<String>,
}
