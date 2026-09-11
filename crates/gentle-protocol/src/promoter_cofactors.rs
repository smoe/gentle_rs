//! Read-only queries over reduced, provenance-bound promoter cofactor packages.
//! Association statistics belong to anchor cohorts, never to individual sites.

use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;

pub const PROMOTER_COFACTOR_SCHEMA: &str = "gentle.promoter_cofactor_query.v1";

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CofactorQuery {
    #[default]
    Inspect,
    Rankings,
    Candidates,
    Anchors,
    AnchorDetail,
    Promoters,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CofactorRanking {
    #[default]
    TaEnriched,
    DnEnriched,
    TaDepleted,
    DnDepleted,
    IsoformDifference,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CofactorRegion {
    pub chromosome: String,
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PromoterCofactorRequest {
    pub package_path: String,
    /// Exact assembly expectation; no implicit aliasing or liftover.
    pub assembly: String,
    pub query: CofactorQuery,
    pub motif: Option<String>,
    /// Exact Ensembl gene ID in promoter ownership, not the cofactor TF name.
    pub gene_id: Option<String>,
    pub region: Option<CofactorRegion>,
    /// Package-local key, meaningful only with the returned package digest.
    pub anchor_id: Option<u64>,
    pub distance_band: Option<String>,
    pub source_species: Option<String>,
    pub ranking: CofactorRanking,
    /// Filters the corresponding precomputed BH q-value, without refitting.
    pub max_q_value: Option<f64>,
    /// Presence only; arbitrary-threshold occurrence counts are unavailable.
    pub presence_threshold: f64,
    pub max_rows: usize,
    pub timeout_seconds: u64,
    pub duckdb_executable: Option<String>,
}

impl Default for PromoterCofactorRequest {
    fn default() -> Self {
        Self {
            package_path: String::new(),
            assembly: "GRCh38".into(),
            query: CofactorQuery::Inspect,
            motif: None,
            gene_id: None,
            region: None,
            anchor_id: None,
            distance_band: None,
            source_species: None,
            ranking: CofactorRanking::TaEnriched,
            max_q_value: None,
            presence_threshold: 0.0,
            max_rows: 200,
            timeout_seconds: 30,
            duckdb_executable: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CofactorAvailability {
    Available,
    PackageMissing,
    InvalidPackage,
    AssemblyMismatch,
    UnsupportedCoverage,
    RuntimeUnavailable,
    QueryFailed,
    RowLimitExceeded,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CofactorCandidate {
    pub gene: String,
    pub group: String,
    pub motif_ids: Vec<String>,
    pub status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CofactorPackageCoverage {
    pub assembly: String,
    pub taxon_id: u64,
    pub chromosomes: Vec<String>,
    pub detailed_motif_ids: Vec<String>,
    pub distance_bands: Vec<String>,
    pub source_score_floor: f64,
    pub positive_threshold: f64,
    pub score_configuration: BTreeMap<String, Value>,
    pub retention: String,
    pub complete_genome_scan: bool,
    pub h3k4me3_model_effects: String,
    pub requested_candidates: Vec<CofactorCandidate>,
}

/// Original estimates and uncertainty are copied, not recomputed by GENtle.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CofactorRankingRow {
    pub motif_id: String,
    pub motif_name: String,
    pub distance_band: String,
    pub anchors_total: u64,
    pub anchors_positive: u64,
    pub positive_anchor_fraction: Option<f64>,
    pub ta_adjusted_odds_ratio: Option<f64>,
    pub dn_adjusted_odds_ratio: Option<f64>,
    pub ta_confidence_interval_95_lower: Option<f64>,
    pub ta_confidence_interval_95_upper: Option<f64>,
    pub dn_confidence_interval_95_lower: Option<f64>,
    pub dn_confidence_interval_95_upper: Option<f64>,
    pub ta_q_value_bh_tax_group: Option<f64>,
    pub dn_q_value_bh_tax_group: Option<f64>,
    pub ta_vs_dn_odds_ratio_ratio: Option<f64>,
    pub confidence_interval_95_lower: Option<f64>,
    pub confidence_interval_95_upper: Option<f64>,
    pub q_value_bh_tax_group: Option<f64>,
    pub ta_evaluation_status: String,
    pub dn_evaluation_status: String,
    pub evaluation_status: String,
    pub source_species: String,
    pub class_support_flag: String,
    /// Preserve series-level diagnostics and future annotation columns verbatim.
    #[serde(flatten)]
    pub source_fields: BTreeMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CofactorAnchor {
    pub anchor_id: u64,
    pub chrom: String,
    pub anchor_start: u64,
    pub anchor_end: u64,
    pub anchor_score: f64,
    /// Per-series support/depth fields; never merge samples or infer raw tracks.
    #[serde(flatten)]
    pub source_fields: BTreeMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CofactorDetail {
    pub anchor_id: u64,
    pub motif_id: String,
    pub distance_band: String,
    pub hit_start: Option<u64>,
    pub hit_end: Option<u64>,
    pub best_score: Option<f64>,
    pub plus_score: Option<f64>,
    pub minus_score: Option<f64>,
    pub best_strand: Option<String>,
    pub interval_distance_bp: Option<i64>,
    pub genomic_side: Option<String>,
    pub n_source_loci: u64,
    pub n_score_zero_loci: u64,
    pub present_at_requested_threshold: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CofactorPromoter {
    pub regulatory_feature_id: String,
    pub chrom: String,
    pub extended_start: u64,
    pub extended_end: u64,
    pub gene_links: Vec<CofactorGeneLink>,
    #[serde(flatten)]
    pub source_fields: BTreeMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CofactorGeneLink {
    pub gene_id: String,
    pub link_source: String,
    pub annotation_release: String,
}

/// Membership is separate from physical observations and gene ownership.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CofactorPromoterMembership {
    pub anchor_id: u64,
    pub regulatory_feature_id: String,
    #[serde(flatten)]
    pub source_fields: BTreeMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PromoterCofactorReport {
    pub schema: String,
    pub report_id: String,
    pub request: PromoterCofactorRequest,
    pub availability: CofactorAvailability,
    pub diagnostic: Option<String>,
    pub package_manifest_sha256: Option<String>,
    pub completion_sha256: Option<String>,
    pub verified_file_sha256: BTreeMap<String, String>,
    pub duckdb_version: Option<String>,
    pub coverage: Option<CofactorPackageCoverage>,
    pub rankings: Vec<CofactorRankingRow>,
    pub more_rankings_available: bool,
    pub anchors: Vec<CofactorAnchor>,
    pub details: Vec<CofactorDetail>,
    pub promoters: Vec<CofactorPromoter>,
    pub memberships: Vec<CofactorPromoterMembership>,
    pub non_claims: Vec<String>,
}
