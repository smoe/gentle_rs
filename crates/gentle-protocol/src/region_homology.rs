//! Query-referenced homology screens for portable genomic regions.
//!
//! These contracts keep observed sequence similarity separate from asserted
//! orthology and from reporter-module interpretation. A target insertion is
//! retained as provenance but never adds a column to the displayed query
//! coordinate system.

use crate::{GenomicRegionOfInterest, GenomicRegionReference, GenomicRegionStrand};
use serde::{Deserialize, Serialize};

pub const GENOMIC_REGION_HOMOLOGY_SCREEN_SCHEMA: &str = "gentle.genomic_region_homology_screen.v1";
pub const PROMOTER_MODULE_ASSESSMENT_SCHEMA: &str = "gentle.promoter_module_assessment.v1";
pub const GENOMIC_REGION_HOMOLOGY_PROJECTION_VERSION: &str = "query_projection_v1";
pub const PROMOTER_SIMILARITY_MATRIX_SCHEMA: &str = "gentle.promoter_similarity_matrix.v1";

pub const fn default_homology_min_identity_percent() -> f64 {
    70.0
}

pub const fn default_homology_min_alignment_length_bp() -> usize {
    24
}

pub const fn default_homology_max_evalue() -> f64 {
    1.0e-3
}

pub const fn default_homology_max_chain_gap_bp() -> usize {
    200
}

pub const fn default_homology_max_loci_per_target() -> usize {
    50
}

pub const fn default_homology_max_hsps_per_target() -> usize {
    100_000
}

pub const fn default_homology_min_conserved_block_bp() -> usize {
    12
}

pub const fn default_promoter_similarity_upstream_bp() -> usize {
    2_000
}

pub const fn default_promoter_similarity_downstream_bp() -> usize {
    200
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct PromoterSimilarityMatrixPolicy {
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    /// Maximum displayed promoter rows. Matching rows beyond this bound remain
    /// counted and are reported as omitted, never silently treated as absent.
    pub max_rows: usize,
}

impl Default for PromoterSimilarityMatrixPolicy {
    fn default() -> Self {
        Self {
            upstream_bp: default_promoter_similarity_upstream_bp(),
            downstream_bp: default_promoter_similarity_downstream_bp(),
            max_rows: 500,
        }
    }
}

pub const fn default_promoter_module_max_partner_gap_bp() -> usize {
    500
}

pub const fn default_promoter_module_max_same_genome_coverage_percent() -> f64 {
    80.0
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionHomologyTargetRole {
    SameGenome,
    ExpectedOrtholog,
    #[default]
    CrossSpeciesUnassigned,
}

impl GenomicRegionHomologyTargetRole {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SameGenome => "same_genome",
            Self::ExpectedOrtholog => "expected_ortholog",
            Self::CrossSpeciesUnassigned => "cross_species_unassigned",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionHomologyTargetStatus {
    #[default]
    Available,
    Unavailable,
    SearchOutputTooBroad,
    NoAcceptedSimilarity,
}

impl GenomicRegionHomologyTargetStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Available => "available",
            Self::Unavailable => "unavailable",
            Self::SearchOutputTooBroad => "search_output_too_broad",
            Self::NoAcceptedSimilarity => "no_accepted_similarity",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionHomologyLocusClass {
    Query,
    ExpectedOrtholog,
    CrossSpeciesUnassigned,
    SameGenomeSelf,
    #[default]
    SameGenomeNonself,
}

impl GenomicRegionHomologyLocusClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Query => "query",
            Self::ExpectedOrtholog => "expected_ortholog",
            Self::CrossSpeciesUnassigned => "cross_species_unassigned",
            Self::SameGenomeSelf => "same_genome_self",
            Self::SameGenomeNonself => "same_genome_nonself",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum GenomicRegionHomologySupportClass {
    ExpectedOrtholog,
    CrossSpeciesUnassigned,
    #[default]
    SameGenomeNonself,
}

impl GenomicRegionHomologySupportClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ExpectedOrtholog => "expected_ortholog",
            Self::CrossSpeciesUnassigned => "cross_species_unassigned",
            Self::SameGenomeNonself => "same_genome_nonself",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionHomologyExpectedLocus {
    pub expected_locus_id: String,
    pub reference: GenomicRegionReference,
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
    pub strand: GenomicRegionStrand,
    pub evidence_id: String,
    pub source_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_sha256: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionHomologyTargetRequest {
    pub genome_id: String,
    pub required: bool,
    pub role: GenomicRegionHomologyTargetRole,
    pub expected_loci: Vec<GenomicRegionHomologyExpectedLocus>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct GenomicRegionHomologySearchPolicy {
    pub min_identity_percent: f64,
    pub min_alignment_length_bp: usize,
    pub max_evalue: f64,
    pub max_chain_gap_bp: usize,
    pub max_loci_per_target: usize,
    pub max_hsps_per_target: usize,
    pub min_conserved_block_bp: usize,
    /// When present, annotate same-genome similarity against transcript-derived
    /// promoter windows and retain an ordered block matrix in the report.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub promoter_similarity_matrix: Option<PromoterSimilarityMatrixPolicy>,
}

impl Default for GenomicRegionHomologySearchPolicy {
    fn default() -> Self {
        Self {
            min_identity_percent: default_homology_min_identity_percent(),
            min_alignment_length_bp: default_homology_min_alignment_length_bp(),
            max_evalue: default_homology_max_evalue(),
            max_chain_gap_bp: default_homology_max_chain_gap_bp(),
            max_loci_per_target: default_homology_max_loci_per_target(),
            max_hsps_per_target: default_homology_max_hsps_per_target(),
            min_conserved_block_bp: default_homology_min_conserved_block_bp(),
            promoter_similarity_matrix: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct PromoterSimilarityBlock {
    pub block_id: String,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    /// Target-promoter coordinates in transcriptional 5'-to-3' orientation.
    pub target_start_0based: u64,
    pub target_end_0based_exclusive: u64,
    /// One-based order among blocks in this target promoter.
    pub target_order: usize,
    pub strand: GenomicRegionStrand,
    pub identity_percent: f64,
    pub bit_score: f64,
    pub source_hsp_ids: Vec<String>,
    /// True when query order and target-promoter order differ from the preceding
    /// block. Such blocks must not be visually joined.
    pub order_break_before: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct PromoterSimilarityMatrixRow {
    pub row_id: String,
    pub target_genome_id: String,
    pub chromosome: String,
    pub promoter_start_0based: u64,
    pub promoter_end_0based_exclusive: u64,
    pub tss_1based: u64,
    pub strand: GenomicRegionStrand,
    pub gene_ids: Vec<String>,
    pub gene_names: Vec<String>,
    pub transcript_ids: Vec<String>,
    pub query_coverage_percent: f64,
    pub mean_identity_percent: f64,
    pub blocks: Vec<PromoterSimilarityBlock>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct PromoterSimilarityMatrix {
    pub schema: String,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub annotated_promoter_window_count: usize,
    pub distinct_gene_count: usize,
    pub distinct_transcript_count: usize,
    pub displayed_row_count: usize,
    pub omitted_row_count: usize,
    /// False when the underlying BLAST/HSP or retained-locus budget prevents a
    /// complete frequency interpretation.
    pub frequency_complete: bool,
    pub rows: Vec<PromoterSimilarityMatrixRow>,
    pub warnings: Vec<String>,
    pub non_claims: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionHomologyScreenRequest {
    pub set_id: String,
    pub region_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expected_region_content_sha256: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub query_genome_id: Option<String>,
    pub targets: Vec<GenomicRegionHomologyTargetRequest>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub catalog_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cache_dir: Option<String>,
    pub policy: GenomicRegionHomologySearchPolicy,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionHomologyEffectiveRequest {
    pub set_id: String,
    pub region_id: String,
    pub region_content_sha256: String,
    pub query_genome_id: String,
    pub targets: Vec<GenomicRegionHomologyTargetRequest>,
    pub catalog_origin: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cache_dir: Option<String>,
    pub policy: GenomicRegionHomologySearchPolicy,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionHomologyDatabaseBinding {
    pub genome_id: String,
    pub index_kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_assembly: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_release: Option<String>,
    pub masking: String,
    pub prefix: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub blast_database_version: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sequence_count: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_bases: Option<u64>,
    pub tool_executable: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tool_version: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub content_fingerprint: Option<String>,
    pub fingerprint_algorithm: String,
    pub validation_status: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionHomologyQueryBinding {
    pub set_id: String,
    pub region: GenomicRegionOfInterest,
    pub sequence: String,
    pub sequence_sha256: String,
    pub source_kind: String,
    pub source_resource_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_fingerprint: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionHomologyHsp {
    pub hsp_id: String,
    pub target_genome_id: String,
    pub subject_id_raw: String,
    pub subject_id: String,
    pub strand: GenomicRegionStrand,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub subject_start_0based: u64,
    pub subject_end_0based_exclusive: u64,
    pub identity_percent: f64,
    pub alignment_length_bp: usize,
    pub mismatches: usize,
    pub gap_opens: usize,
    pub evalue: f64,
    pub bit_score: f64,
    pub aligned_query: String,
    pub aligned_subject: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionHomologyOmittedInsertion {
    pub insertion_id: String,
    pub alignment_row_id: String,
    pub hsp_id: String,
    pub query_anchor_0based: usize,
    pub target_start_0based: u64,
    pub target_end_0based_exclusive: u64,
    pub length_bp: usize,
    pub target_sequence: String,
    pub strand: GenomicRegionStrand,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionHomologyLocus {
    pub locus_id: String,
    pub target_genome_id: String,
    pub target_role: GenomicRegionHomologyTargetRole,
    pub locus_class: GenomicRegionHomologyLocusClass,
    pub subject_id: String,
    pub strand: GenomicRegionStrand,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub subject_start_0based: u64,
    pub subject_end_0based_exclusive: u64,
    pub identity_percent: f64,
    pub query_coverage_percent: f64,
    pub bit_score: f64,
    pub source_hsp_ids: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub orthology_evidence_id: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionHomologyAlignmentConflict {
    pub query_position_0based: usize,
    pub retained_hsp_id: String,
    pub rejected_hsp_id: String,
    pub retained_symbol: String,
    pub rejected_symbol: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GenomicRegionHomologyAlignmentRow {
    pub row_id: String,
    pub locus_id: String,
    pub target_genome_id: String,
    pub locus_class: GenomicRegionHomologyLocusClass,
    pub subject_id: String,
    pub strand: GenomicRegionStrand,
    pub subject_start_0based: u64,
    pub subject_end_0based_exclusive: u64,
    /// One byte per query base: `.`, a substituted base, `-`, or a blank.
    pub query_projection: String,
    pub exact_match_count: usize,
    pub covered_query_base_count: usize,
    pub omitted_insertion_ids: Vec<String>,
    pub source_hsp_ids: Vec<String>,
    pub conflicts: Vec<GenomicRegionHomologyAlignmentConflict>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionHomologyTargetResult {
    pub target: GenomicRegionHomologyTargetRequest,
    pub status: GenomicRegionHomologyTargetStatus,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub database: Option<GenomicRegionHomologyDatabaseBinding>,
    pub raw_hsp_count: usize,
    pub accepted_hsp_count: usize,
    pub retained_locus_count: usize,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionHomologyConservedBlock {
    pub block_id: String,
    pub support_class: GenomicRegionHomologySupportClass,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub query_sequence: String,
    pub genomic_start_0based: u64,
    pub genomic_end_0based_exclusive: u64,
    pub genomic_strand: GenomicRegionStrand,
    pub supporting_row_ids: Vec<String>,
    pub supporting_genome_ids: Vec<String>,
    pub available_genome_ids: Vec<String>,
    pub unavailable_genome_ids: Vec<String>,
    pub support_fraction: f64,
    pub mean_identity_percent: f64,
    pub source_hsp_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GenomicRegionHomologyScreenReport {
    pub schema: String,
    pub projection_version: String,
    pub request_sha256: String,
    pub effective_request_sha256: String,
    pub content_sha256: String,
    pub query: GenomicRegionHomologyQueryBinding,
    pub effective_request: GenomicRegionHomologyEffectiveRequest,
    pub targets: Vec<GenomicRegionHomologyTargetResult>,
    pub hsps: Vec<GenomicRegionHomologyHsp>,
    pub loci: Vec<GenomicRegionHomologyLocus>,
    pub alignment_rows: Vec<GenomicRegionHomologyAlignmentRow>,
    pub omitted_insertions: Vec<GenomicRegionHomologyOmittedInsertion>,
    pub conserved_blocks: Vec<GenomicRegionHomologyConservedBlock>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub promoter_similarity_matrix: Option<PromoterSimilarityMatrix>,
    pub same_genome_nonself_locus_count: usize,
    pub same_genome_nonself_query_coverage_percent: f64,
    pub warnings: Vec<String>,
    pub non_claims: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct PromoterModuleEvidenceSpan {
    pub evidence_id: String,
    pub evidence_kind: String,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub required: bool,
    pub source_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_sha256: Option<String>,
    pub evidence_statement: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct PromoterModuleAssessmentRequest {
    pub homology_report: Box<GenomicRegionHomologyScreenReport>,
    pub selected_evidence_spans: Vec<PromoterModuleEvidenceSpan>,
    pub max_partner_gap_bp: usize,
    /// Allowed absolute difference between query and ortholog gaps; zero is exact.
    pub max_partner_gap_difference_bp: usize,
    pub max_same_genome_query_coverage_percent: f64,
}

impl Default for PromoterModuleAssessmentRequest {
    fn default() -> Self {
        Self {
            homology_report: Box::default(),
            selected_evidence_spans: vec![],
            max_partner_gap_bp: default_promoter_module_max_partner_gap_bp(),
            max_partner_gap_difference_bp: 0,
            max_same_genome_query_coverage_percent:
                default_promoter_module_max_same_genome_coverage_percent(),
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PromoterModuleHypothesisKind {
    StandaloneReporterCandidate,
    PairedContextCandidate,
    RepetitiveOrAmbiguous,
    #[default]
    InsufficientEvidence,
}

impl PromoterModuleHypothesisKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::StandaloneReporterCandidate => "standalone_reporter_candidate",
            Self::PairedContextCandidate => "paired_context_candidate",
            Self::RepetitiveOrAmbiguous => "repetitive_or_ambiguous",
            Self::InsufficientEvidence => "insufficient_evidence",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct PromoterModuleDecisionRule {
    pub rule_id: String,
    pub description: String,
    pub satisfied: bool,
    pub evidence_ids: Vec<String>,
    pub block_ids: Vec<String>,
    pub detail: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct PromoterModuleAlternativeFragment {
    pub fragment_id: String,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub block_ids: Vec<String>,
    pub covered_evidence_ids: Vec<String>,
    pub rationale: String,
}

/// One conserved block mapped through the winning HSPs onto a target locus.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct PromoterModuleTargetBlock {
    pub block_id: String,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub target_start_0based: u64,
    pub target_end_0based_exclusive: u64,
}

/// Evidence for (or against) ordered, spacing-compatible blocks in one ortholog.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct PromoterModulePartnerContext {
    pub row_id: String,
    pub locus_id: String,
    pub genome_id: String,
    pub subject_id: String,
    pub strand: GenomicRegionStrand,
    pub source_hsp_ids: Vec<String>,
    pub blocks: Vec<PromoterModuleTargetBlock>,
    pub query_gaps_bp: Vec<usize>,
    pub target_gaps_bp: Vec<Option<u64>>,
    pub passed: bool,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct PromoterModuleAssessmentReport {
    pub schema: String,
    pub request_sha256: String,
    pub content_sha256: String,
    pub homology_report_sha256: String,
    pub hypothesis: PromoterModuleHypothesisKind,
    pub selected_evidence_spans: Vec<PromoterModuleEvidenceSpan>,
    pub selected_block_ids: Vec<String>,
    pub decision_trace: Vec<PromoterModuleDecisionRule>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub partner_contexts: Vec<PromoterModulePartnerContext>,
    pub alternative_fragments: Vec<PromoterModuleAlternativeFragment>,
    pub suggested_validation: Vec<String>,
    pub warnings: Vec<String>,
    pub non_claims: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
}
