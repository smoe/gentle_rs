//! Stable public analysis/report contracts extracted from the monolithic
//! engine.
//!
//! This is the narrow, serialization-friendly layer that GUI/CLI/JS/Lua/Python
//! adapters should lean on when they only need record shapes rather than full
//! engine execution.
//!
//! Look here for:
//! - persisted report/result payloads such as RNA-read, dotplot, planning, and
//!   sequencing-confirmation records
//! - small enums that appear in JSON-facing operation contracts
//! - state-summary structs that should remain slower-changing than
//!   `src/engine.rs`

pub use gentle_protocol::{
    DotplotMode, EngineError, ErrorCode, FlexibilityModel, PairwiseAlignmentMode,
    RnaReadAlignmentBackend, RnaReadAlignmentEffect, RnaReadAlignmentInspectionEffectFilter,
    RnaReadAlignmentInspectionSortKey, RnaReadAlignmentMode, RnaReadHitSelection,
    RnaReadInputFormat, RnaReadInterpretationProfile, RnaReadOriginClass, RnaReadOriginMode,
    RnaReadReportMode, RnaReadScoreDensityScale, RnaReadScoreDensityVariant,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

use super::{
    DisplaySettings, OpId, Operation, PrepareGenomeProgress, ProtocolCartoonTemplateBindings,
    RunId, SeqId, SequencingConfirmationReport, SplicingScopePreset,
};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct RnaReadSeedFilterConfig {
    pub kmer_len: usize,
    #[serde(default = "super::default_rna_seed_stride_bp")]
    pub seed_stride_bp: usize,
    pub min_seed_hit_fraction: f64,
    #[serde(default = "super::default_min_weighted_seed_hit_fraction")]
    pub min_weighted_seed_hit_fraction: f64,
    #[serde(default = "super::default_min_unique_matched_kmers")]
    pub min_unique_matched_kmers: usize,
    #[serde(default = "super::default_max_median_transcript_gap")]
    pub max_median_transcript_gap: f64,
    #[serde(default = "super::default_min_chain_consistency_fraction")]
    pub min_chain_consistency_fraction: f64,
    #[serde(default = "super::default_min_confirmed_exon_transitions")]
    pub min_confirmed_exon_transitions: usize,
    #[serde(default = "super::default_min_transition_support_fraction")]
    pub min_transition_support_fraction: f64,
    #[serde(default = "super::default_true")]
    pub cdna_poly_t_flip_enabled: bool,
    #[serde(default = "super::default_poly_t_prefix_min_bp")]
    pub poly_t_prefix_min_bp: usize,
}

impl Default for RnaReadSeedFilterConfig {
    fn default() -> Self {
        Self {
            kmer_len: 10,
            seed_stride_bp: 1,
            min_seed_hit_fraction: 0.30,
            min_weighted_seed_hit_fraction: 0.05,
            min_unique_matched_kmers: 12,
            max_median_transcript_gap: 4.0,
            min_chain_consistency_fraction: 0.40,
            min_confirmed_exon_transitions: 1,
            min_transition_support_fraction: 0.05,
            cdna_poly_t_flip_enabled: true,
            poly_t_prefix_min_bp: 18,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct RnaReadAlignConfig {
    pub band_width_bp: usize,
    pub min_identity_fraction: f64,
    pub max_secondary_mappings: usize,
}

impl Default for RnaReadAlignConfig {
    fn default() -> Self {
        Self {
            band_width_bp: 24,
            min_identity_fraction: 0.60,
            max_secondary_mappings: 3,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadMappingHit {
    #[serde(default)]
    pub alignment_mode: RnaReadAlignmentMode,
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    #[serde(default)]
    pub query_reverse_complemented: bool,
    pub target_start_1based: usize,
    pub target_end_1based: usize,
    #[serde(default)]
    pub target_start_offset_0based: usize,
    #[serde(default)]
    pub target_end_offset_0based_exclusive: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub score: isize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentDisplay {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    #[serde(default)]
    pub alignment_mode: RnaReadAlignmentMode,
    #[serde(default)]
    pub query_reverse_complemented: bool,
    #[serde(default)]
    pub query_start_0based: usize,
    #[serde(default)]
    pub query_end_0based_exclusive: usize,
    #[serde(default)]
    pub target_start_1based: usize,
    #[serde(default)]
    pub target_end_1based: usize,
    #[serde(default)]
    pub target_start_offset_0based: usize,
    #[serde(default)]
    pub target_end_offset_0based_exclusive: usize,
    #[serde(default)]
    pub target_length_bp: usize,
    pub score: isize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
    #[serde(default)]
    pub target_coverage_fraction: f64,
    #[serde(default)]
    pub matches: usize,
    #[serde(default)]
    pub mismatches: usize,
    #[serde(default)]
    pub insertions: usize,
    #[serde(default)]
    pub deletions: usize,
    #[serde(default)]
    pub aligned_columns: usize,
    #[serde(default)]
    pub aligned_query: String,
    #[serde(default)]
    pub aligned_midline: String,
    #[serde(default)]
    pub aligned_target: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadStrandAssignmentDiagnostics {
    pub selected_strand: String,
    pub selected_reason: String,
    pub selected_transition_hits: usize,
    pub selected_exon_hits: usize,
    pub plus_best_transcript_id: String,
    pub plus_best_transition_hits: usize,
    pub plus_best_exon_hits: usize,
    pub minus_best_transcript_id: String,
    pub minus_best_transition_hits: usize,
    pub minus_best_exon_hits: usize,
    pub competing_opposite_strand: bool,
    pub ambiguous_near_tie: bool,
    pub chain_preferred_strand: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadOriginCandidateContribution {
    pub candidate_role: String,
    pub transcript_id: String,
    pub strand: String,
    pub transition_hits: usize,
    pub exon_hits: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadInterpretationHit {
    pub record_index: usize,
    pub source_byte_offset: usize,
    pub header_id: String,
    pub sequence: String,
    pub read_length_bp: usize,
    pub tested_kmers: usize,
    pub matched_kmers: usize,
    pub seed_hit_fraction: f64,
    #[serde(default)]
    pub weighted_seed_hit_fraction: f64,
    #[serde(default)]
    pub weighted_matched_kmers: f64,
    #[serde(default)]
    pub seed_chain_transcript_id: String,
    #[serde(default)]
    pub seed_chain_support_kmers: usize,
    #[serde(default)]
    pub seed_chain_support_fraction: f64,
    #[serde(default)]
    pub seed_median_transcript_gap: f64,
    #[serde(default)]
    pub seed_transcript_gap_count: usize,
    #[serde(default)]
    pub exon_path_transcript_id: String,
    #[serde(default)]
    pub exon_path: String,
    #[serde(default)]
    pub exon_transitions_confirmed: usize,
    #[serde(default)]
    pub exon_transitions_total: usize,
    #[serde(default)]
    pub reverse_complement_applied: bool,
    #[serde(default)]
    pub strand_diagnostics: RnaReadStrandAssignmentDiagnostics,
    #[serde(default)]
    pub origin_class: RnaReadOriginClass,
    #[serde(default)]
    pub origin_reason: String,
    #[serde(default)]
    pub origin_confidence: f64,
    #[serde(default)]
    pub strand_confidence: f64,
    #[serde(default)]
    pub origin_candidates: Vec<RnaReadOriginCandidateContribution>,
    pub perfect_seed_match: bool,
    pub passed_seed_filter: bool,
    #[serde(default)]
    pub msa_eligible: bool,
    #[serde(default)]
    pub msa_eligibility_reason: String,
    pub best_mapping: Option<RnaReadMappingHit>,
    pub secondary_mappings: Vec<RnaReadMappingHit>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadExonSupportFrequency {
    pub start_1based: usize,
    pub end_1based: usize,
    pub support_read_count: usize,
    pub support_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadJunctionSupportFrequency {
    pub donor_1based: usize,
    pub acceptor_1based: usize,
    pub support_read_count: usize,
    pub support_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadTransitionSupportRow {
    pub from_exon_ordinal: usize,
    pub to_exon_ordinal: usize,
    pub from_start_1based: usize,
    pub from_end_1based: usize,
    pub to_start_1based: usize,
    pub to_end_1based: usize,
    pub support_read_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadIsoformSupportRow {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub exon_count: usize,
    pub expected_transition_count: usize,
    pub reads_assigned: usize,
    pub reads_seed_passed: usize,
    pub transition_rows_supported: usize,
    pub transition_rows_supported_fraction: f64,
    pub mean_seed_median_gap: f64,
    pub mean_confirmed_transition_fraction: f64,
    pub best_seed_hit_fraction: f64,
    pub best_weighted_seed_hit_fraction: f64,
    #[serde(default)]
    pub reads_chain_same_strand: usize,
    #[serde(default)]
    pub reads_with_opposite_strand_competition: usize,
    #[serde(default)]
    pub reads_ambiguous_strand_ties: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadMappedIsoformSupportRow {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    #[serde(default)]
    pub aligned_read_count: usize,
    #[serde(default)]
    pub msa_eligible_read_count: usize,
    #[serde(default)]
    pub mean_identity_fraction: f64,
    #[serde(default)]
    pub mean_query_coverage_fraction: f64,
    #[serde(default)]
    pub best_alignment_score: isize,
    #[serde(default)]
    pub secondary_mapping_total: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadSampleSheetExport {
    pub schema: String,
    pub path: String,
    pub report_count: usize,
    pub appended: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadExonPathsExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub selection: RnaReadHitSelection,
    pub row_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadExonAbundanceExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub selection: RnaReadHitSelection,
    pub selected_read_count: usize,
    pub exon_row_count: usize,
    pub transition_row_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadScoreDensitySvgExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub scale: RnaReadScoreDensityScale,
    #[serde(default)]
    pub variant: RnaReadScoreDensityVariant,
    pub bin_count: usize,
    pub max_bin_count: u64,
    pub total_scored_reads: u64,
    pub derived_from_report_hits_only: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentDotplotSvgExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub selection: RnaReadHitSelection,
    pub point_count: usize,
    pub rendered_point_count: usize,
    pub max_points: usize,
    pub min_score: isize,
    pub max_score: isize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentTsvExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub selection: RnaReadHitSelection,
    pub row_count: usize,
    pub aligned_count: usize,
    pub limit: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentInspectionSubsetSpec {
    pub effect_filter: RnaReadAlignmentInspectionEffectFilter,
    pub sort_key: RnaReadAlignmentInspectionSortKey,
    pub search: String,
    pub selected_record_indices: Vec<usize>,
    #[serde(default)]
    pub score_density_variant: RnaReadScoreDensityVariant,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub score_density_seed_filter_override: Option<RnaReadSeedFilterConfig>,
    pub score_bin_index: Option<usize>,
    pub score_bin_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadMappedSupportExonAttribution {
    pub start_1based: usize,
    pub end_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadMappedSupportJunctionAttribution {
    pub donor_1based: usize,
    pub acceptor_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentInspectionRow {
    pub rank: usize,
    pub record_index: usize,
    pub header_id: String,
    #[serde(default)]
    pub phase1_primary_transcript_id: String,
    #[serde(default)]
    pub seed_chain_transcript_id: String,
    #[serde(default)]
    pub exon_path_transcript_id: String,
    #[serde(default)]
    pub exon_path: String,
    #[serde(default)]
    pub exon_transitions_confirmed: usize,
    #[serde(default)]
    pub exon_transitions_total: usize,
    #[serde(default)]
    pub selected_strand: String,
    #[serde(default)]
    pub reverse_complement_applied: bool,
    #[serde(default)]
    pub alignment_effect: RnaReadAlignmentEffect,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub alignment_mode: RnaReadAlignmentMode,
    #[serde(default)]
    pub target_start_1based: usize,
    #[serde(default)]
    pub target_end_1based: usize,
    #[serde(default)]
    pub target_length_bp: usize,
    pub score: isize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
    #[serde(default)]
    pub target_coverage_fraction: f64,
    #[serde(default)]
    pub secondary_mapping_count: usize,
    pub seed_hit_fraction: f64,
    pub weighted_seed_hit_fraction: f64,
    pub passed_seed_filter: bool,
    pub msa_eligible: bool,
    pub origin_class: RnaReadOriginClass,
    #[serde(default)]
    pub mapped_exon_support: Vec<RnaReadMappedSupportExonAttribution>,
    #[serde(default)]
    pub mapped_junction_support: Vec<RnaReadMappedSupportJunctionAttribution>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentInspection {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub selection: RnaReadHitSelection,
    pub row_count: usize,
    pub aligned_count: usize,
    pub subset_match_count: usize,
    pub limit: usize,
    pub subset_spec: RnaReadAlignmentInspectionSubsetSpec,
    pub align_min_identity_fraction: f64,
    pub max_secondary_mappings: usize,
    pub rows: Vec<RnaReadAlignmentInspectionRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadPairwiseAlignmentDetail {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub record_index: usize,
    pub header_id: String,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    #[serde(default)]
    pub alignment_mode: RnaReadAlignmentMode,
    #[serde(default)]
    pub backend: RnaReadAlignmentBackend,
    pub query_length_bp: usize,
    pub target_length_bp: usize,
    pub aligned_query_start_0based: usize,
    pub aligned_query_end_0based_exclusive: usize,
    pub aligned_target_start_offset_0based: usize,
    pub aligned_target_end_offset_0based_exclusive: usize,
    pub target_start_1based: usize,
    pub target_end_1based: usize,
    pub aligned_columns: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub score: isize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
    #[serde(default)]
    pub target_coverage_fraction: f64,
    pub cigar: String,
    pub aligned_query: String,
    pub aligned_relation: String,
    pub aligned_target: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaSeedHashCatalogEntry {
    pub seed_bits: u32,
    pub kmer_sequence: String,
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub template_offset_0based: usize,
    pub genomic_pos_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaSeedHashTemplateAuditEntry {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub template_sequence: String,
    pub template_length_bp: usize,
    pub template_first_genomic_pos_1based: usize,
    pub template_last_genomic_pos_1based: usize,
    pub reverse_complemented_from_genome: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadInterpretationReport {
    pub schema: String,
    pub report_id: String,
    #[serde(default)]
    pub report_mode: RnaReadReportMode,
    pub seq_id: String,
    pub seed_feature_id: usize,
    pub generated_at_unix_ms: u128,
    pub profile: RnaReadInterpretationProfile,
    pub input_path: String,
    pub input_format: RnaReadInputFormat,
    pub scope: SplicingScopePreset,
    #[serde(default)]
    pub origin_mode: RnaReadOriginMode,
    #[serde(default)]
    pub target_gene_ids: Vec<String>,
    #[serde(default)]
    pub roi_seed_capture_enabled: bool,
    #[serde(default)]
    pub checkpoint_path: Option<String>,
    #[serde(default = "super::default_rna_read_checkpoint_every_reads")]
    pub checkpoint_every_reads: usize,
    #[serde(default)]
    pub resumed_from_checkpoint: bool,
    pub seed_filter: RnaReadSeedFilterConfig,
    pub align_config: RnaReadAlignConfig,
    pub read_count_total: usize,
    pub read_count_seed_passed: usize,
    pub read_count_aligned: usize,
    #[serde(default)]
    pub retained_count_msa_eligible: usize,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub hits: Vec<RnaReadInterpretationHit>,
    #[serde(default)]
    pub exon_support_frequencies: Vec<RnaReadExonSupportFrequency>,
    #[serde(default)]
    pub junction_support_frequencies: Vec<RnaReadJunctionSupportFrequency>,
    #[serde(default)]
    pub transition_support_rows: Vec<RnaReadTransitionSupportRow>,
    #[serde(default)]
    pub isoform_support_rows: Vec<RnaReadIsoformSupportRow>,
    #[serde(default)]
    pub mapped_isoform_support_rows: Vec<RnaReadMappedIsoformSupportRow>,
    #[serde(default)]
    pub origin_class_counts: BTreeMap<String, usize>,
    #[serde(default)]
    pub score_density_bins: Vec<u64>,
    #[serde(default)]
    pub seed_pass_score_density_bins: Vec<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaReadInterpretationReportSummary {
    pub report_id: String,
    #[serde(default)]
    pub report_mode: RnaReadReportMode,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    pub profile: RnaReadInterpretationProfile,
    pub input_path: String,
    pub input_format: RnaReadInputFormat,
    pub seed_feature_id: usize,
    pub scope: SplicingScopePreset,
    #[serde(default)]
    pub origin_mode: RnaReadOriginMode,
    #[serde(default)]
    pub target_gene_count: usize,
    #[serde(default)]
    pub roi_seed_capture_enabled: bool,
    pub read_count_total: usize,
    pub read_count_seed_passed: usize,
    pub read_count_aligned: usize,
    #[serde(default)]
    pub retained_count_msa_eligible: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotMatchPoint {
    pub x_0based: usize,
    pub y_0based: usize,
    pub mismatches: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotBoxplotBin {
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub hit_count: usize,
    pub min_reference_0based: Option<usize>,
    pub q1_reference_0based: Option<usize>,
    pub median_reference_0based: Option<usize>,
    pub q3_reference_0based: Option<usize>,
    pub max_reference_0based: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotQuerySeries {
    pub series_id: String,
    pub seq_id: String,
    pub label: String,
    pub color_rgb: [u8; 3],
    #[serde(default)]
    pub mode: DotplotMode,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub point_count: usize,
    pub points: Vec<DotplotMatchPoint>,
    pub boxplot_bin_count: usize,
    pub boxplot_bins: Vec<DotplotBoxplotBin>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotReferenceAnnotationInterval {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub label: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotReferenceAnnotationTrack {
    pub seq_id: String,
    pub label: String,
    pub interval_count: usize,
    pub intervals: Vec<DotplotReferenceAnnotationInterval>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotOverlayQuerySpec {
    pub seq_id: String,
    pub label: String,
    #[serde(default)]
    pub span_start_0based: Option<usize>,
    #[serde(default)]
    pub span_end_0based: Option<usize>,
    #[serde(default)]
    pub mode: DotplotMode,
    #[serde(default)]
    pub color_rgb: Option<[u8; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotView {
    pub schema: String,
    pub dotplot_id: String,
    pub owner_seq_id: String,
    pub seq_id: String,
    pub reference_seq_id: Option<String>,
    pub generated_at_unix_ms: u128,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub reference_span_start_0based: usize,
    pub reference_span_end_0based: usize,
    pub mode: DotplotMode,
    pub word_size: usize,
    pub step_bp: usize,
    pub max_mismatches: usize,
    pub tile_bp: Option<usize>,
    pub point_count: usize,
    pub points: Vec<DotplotMatchPoint>,
    pub boxplot_bin_count: usize,
    pub boxplot_bins: Vec<DotplotBoxplotBin>,
    pub series_count: usize,
    pub query_series: Vec<DotplotQuerySeries>,
    #[serde(default)]
    pub reference_annotation: Option<DotplotReferenceAnnotationTrack>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DotplotViewSummary {
    pub dotplot_id: String,
    pub owner_seq_id: String,
    pub seq_id: String,
    pub reference_seq_id: Option<String>,
    pub generated_at_unix_ms: u128,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub reference_span_start_0based: usize,
    pub reference_span_end_0based: usize,
    pub mode: DotplotMode,
    pub word_size: usize,
    pub step_bp: usize,
    pub max_mismatches: usize,
    pub point_count: usize,
    pub series_count: usize,
}

impl DotplotView {
    pub fn normalize_v3_defaults(&mut self) {
        if self.owner_seq_id.trim().is_empty() {
            self.owner_seq_id = self.seq_id.clone();
        }
        if self.query_series.is_empty() {
            let label = if self.seq_id.trim().is_empty() {
                "<query>".to_string()
            } else {
                self.seq_id.clone()
            };
            self.query_series.push(DotplotQuerySeries {
                series_id: if self.dotplot_id.trim().is_empty() {
                    "series_1".to_string()
                } else {
                    format!("{}_series_1", self.dotplot_id)
                },
                seq_id: self.seq_id.clone(),
                label,
                color_rgb: [29, 78, 216],
                mode: self.mode,
                span_start_0based: self.span_start_0based,
                span_end_0based: self.span_end_0based,
                point_count: if self.point_count == 0 {
                    self.points.len()
                } else {
                    self.point_count
                },
                points: self.points.clone(),
                boxplot_bin_count: if self.boxplot_bin_count == 0 {
                    self.boxplot_bins.len()
                } else {
                    self.boxplot_bin_count
                },
                boxplot_bins: self.boxplot_bins.clone(),
            });
        }
        self.series_count = self.query_series.len();
        if let Some(primary) = self.query_series.first() {
            if self.seq_id.trim().is_empty() {
                self.seq_id = primary.seq_id.clone();
            }
            self.mode = primary.mode;
            self.span_start_0based = primary.span_start_0based;
            self.span_end_0based = primary.span_end_0based;
            self.point_count = primary.point_count.max(primary.points.len());
            self.points = primary.points.clone();
            self.boxplot_bin_count = primary.boxplot_bin_count.max(primary.boxplot_bins.len());
            self.boxplot_bins = primary.boxplot_bins.clone();
        }
        if let Some(annotation) = self.reference_annotation.as_mut() {
            annotation.interval_count = annotation.intervals.len();
        }
    }

    pub fn primary_series(&self) -> Option<&DotplotQuerySeries> {
        self.query_series.first()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceAlignmentReport {
    pub schema: String,
    pub mode: PairwiseAlignmentMode,
    pub query_seq_id: String,
    pub target_seq_id: String,
    pub query_span_start_0based: usize,
    pub query_span_end_0based: usize,
    pub target_span_start_0based: usize,
    pub target_span_end_0based: usize,
    pub aligned_query_start_0based: usize,
    pub aligned_query_end_0based_exclusive: usize,
    pub aligned_target_start_0based: usize,
    pub aligned_target_end_0based_exclusive: usize,
    pub score: i32,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub aligned_columns: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
    pub target_coverage_fraction: f64,
    pub cigar: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// One genome-position bin used for running RNA-read seed-confirmation statistics.
pub struct RnaReadSeedHistogramBin {
    pub start_1based: usize,
    pub end_1based: usize,
    pub confirmed_plus: u64,
    pub confirmed_minus: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
/// Lightweight top-hit row included in running RNA-read progress updates.
pub struct RnaReadTopHitPreview {
    pub record_index: usize,
    pub header_id: String,
    pub seed_hit_fraction: f64,
    pub weighted_seed_hit_fraction: f64,
    #[serde(default)]
    pub weighted_matched_kmers: f64,
    #[serde(default)]
    pub seed_chain_transcript_id: String,
    #[serde(default)]
    pub seed_chain_support_kmers: usize,
    #[serde(default)]
    pub seed_chain_support_fraction: f64,
    #[serde(default)]
    pub seed_median_transcript_gap: f64,
    #[serde(default)]
    pub seed_transcript_gap_count: usize,
    pub matched_kmers: usize,
    pub tested_kmers: usize,
    pub passed_seed_filter: bool,
    #[serde(default)]
    pub reverse_complement_applied: bool,
    #[serde(default)]
    pub selected_strand: String,
    #[serde(default)]
    pub competing_opposite_strand: bool,
    #[serde(default)]
    pub ambiguous_strand_tie: bool,
    #[serde(default)]
    pub origin_class: RnaReadOriginClass,
    #[serde(default)]
    pub origin_reason: String,
    #[serde(default)]
    pub origin_confidence: f64,
    #[serde(default)]
    pub strand_confidence: f64,
    #[serde(default)]
    pub origin_candidates: Vec<RnaReadOriginCandidateContribution>,
    #[serde(default)]
    pub msa_eligible: bool,
    #[serde(default)]
    pub msa_eligibility_reason: String,
    #[serde(default)]
    pub aligned: bool,
    #[serde(default)]
    pub best_alignment_mode: String,
    #[serde(default)]
    pub best_alignment_transcript_id: String,
    #[serde(default)]
    pub best_alignment_transcript_label: String,
    #[serde(default)]
    pub best_alignment_strand: String,
    #[serde(default)]
    pub best_alignment_target_start_1based: usize,
    #[serde(default)]
    pub best_alignment_target_end_1based: usize,
    #[serde(default)]
    pub best_alignment_identity_fraction: f64,
    #[serde(default)]
    pub best_alignment_query_coverage_fraction: f64,
    #[serde(default)]
    pub best_alignment_score: isize,
    #[serde(default)]
    pub secondary_mapping_count: usize,
    pub read_length_bp: usize,
    pub sequence: String,
    pub sequence_preview: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Progress payload emitted by RNA-read interpretation operations.
pub struct RnaReadInterpretProgress {
    pub seq_id: String,
    pub reads_processed: usize,
    pub reads_total: usize,
    #[serde(default)]
    pub read_bases_processed: u64,
    #[serde(default)]
    pub mean_read_length_bp: f64,
    #[serde(default)]
    pub median_read_length_bp: usize,
    #[serde(default)]
    pub p95_read_length_bp: usize,
    #[serde(default)]
    pub input_bytes_processed: u64,
    #[serde(default)]
    pub input_bytes_total: u64,
    pub seed_passed: usize,
    pub aligned: usize,
    pub tested_kmers: usize,
    pub matched_kmers: usize,
    #[serde(default)]
    pub seed_compute_ms: f64,
    #[serde(default)]
    pub align_compute_ms: f64,
    #[serde(default)]
    pub io_read_ms: f64,
    #[serde(default)]
    pub fasta_parse_ms: f64,
    #[serde(default)]
    pub normalize_compute_ms: f64,
    #[serde(default)]
    pub inference_compute_ms: f64,
    #[serde(default)]
    pub progress_emit_ms: f64,
    pub update_every_reads: usize,
    pub done: bool,
    pub bins: Vec<RnaReadSeedHistogramBin>,
    pub score_density_bins: Vec<u64>,
    #[serde(default)]
    pub seed_pass_score_density_bins: Vec<u64>,
    #[serde(default)]
    pub top_hits_preview: Vec<RnaReadTopHitPreview>,
    #[serde(default)]
    pub transition_support_rows: Vec<RnaReadTransitionSupportRow>,
    #[serde(default)]
    pub isoform_support_rows: Vec<RnaReadIsoformSupportRow>,
    #[serde(default)]
    pub mapped_exon_support_frequencies: Vec<RnaReadExonSupportFrequency>,
    #[serde(default)]
    pub mapped_junction_support_frequencies: Vec<RnaReadJunctionSupportFrequency>,
    #[serde(default)]
    pub mapped_isoform_support_rows: Vec<RnaReadMappedIsoformSupportRow>,
    #[serde(default)]
    pub reads_with_transition_support: usize,
    #[serde(default)]
    pub transition_confirmations: usize,
    #[serde(default)]
    pub junction_crossing_seed_bits_indexed: usize,
    #[serde(default)]
    pub origin_class_counts: BTreeMap<String, usize>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequenceFeatureRangeRelation {
    #[default]
    Overlap,
    Within,
    Contains,
}

impl SequenceFeatureRangeRelation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Overlap => "overlap",
            Self::Within => "within",
            Self::Contains => "contains",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequenceFeatureStrandFilter {
    #[default]
    Any,
    Forward,
    Reverse,
}

impl SequenceFeatureStrandFilter {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Any => "any",
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequenceFeatureSortBy {
    FeatureId,
    #[default]
    Start,
    End,
    Kind,
    Length,
}

impl SequenceFeatureSortBy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FeatureId => "feature_id",
            Self::Start => "start",
            Self::End => "end",
            Self::Kind => "kind",
            Self::Length => "length",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQualifierFilter {
    pub key: String,
    pub value_contains: Option<String>,
    pub value_regex: Option<String>,
    pub case_sensitive: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQuery {
    pub seq_id: SeqId,
    pub include_source: bool,
    pub include_qualifiers: bool,
    pub kind_in: Vec<String>,
    pub kind_not_in: Vec<String>,
    pub start_0based: Option<usize>,
    pub end_0based_exclusive: Option<usize>,
    pub range_relation: SequenceFeatureRangeRelation,
    pub strand: SequenceFeatureStrandFilter,
    pub label_contains: Option<String>,
    pub label_regex: Option<String>,
    pub qualifier_filters: Vec<SequenceFeatureQualifierFilter>,
    pub min_len_bp: Option<usize>,
    pub max_len_bp: Option<usize>,
    pub limit: Option<usize>,
    pub offset: usize,
    pub sort_by: SequenceFeatureSortBy,
    pub descending: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQueryRow {
    pub feature_id: usize,
    pub kind: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub length_bp: usize,
    pub strand: String,
    pub label: String,
    pub labels: Vec<String>,
    pub qualifiers: BTreeMap<String, Vec<String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQueryResult {
    pub schema: String,
    pub seq_id: SeqId,
    pub sequence_length_bp: usize,
    pub total_feature_count: usize,
    pub matched_count: usize,
    pub returned_count: usize,
    pub offset: usize,
    pub limit: usize,
    pub range_relation: String,
    pub strand_filter: String,
    pub sort_by: String,
    pub descending: bool,
    pub query: SequenceFeatureQuery,
    pub rows: Vec<SequenceFeatureQueryRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// One deterministic workflow run: ordered operations with a caller-supplied
/// `run_id`.
///
/// Operations are applied sequentially. The current workflow executor is not
/// transactional: if a later step fails, earlier successful steps remain in
/// state and in the operation journal.
pub struct Workflow {
    pub run_id: RunId,
    pub ops: Vec<Operation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Canonical result payload returned after one operation completes.
///
/// `created_seq_ids` and `changed_seq_ids` are the stable adapter-facing hint
/// for which sequence windows/views may need refresh after an operation.
pub struct OpResult {
    pub op_id: OpId,
    pub created_seq_ids: Vec<SeqId>,
    pub changed_seq_ids: Vec<SeqId>,
    pub warnings: Vec<String>,
    pub messages: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub protocol_cartoon_preview: Option<ProtocolCartoonPreviewTelemetry>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_annotation_projection: Option<GenomeAnnotationProjectionTelemetry>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequence_alignment: Option<SequenceAlignmentReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequencing_confirmation_report: Option<SequencingConfirmationReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Optional protocol-cartoon preview payload emitted by operations that can
/// project a deterministic mechanism strip from operation geometry.
pub struct ProtocolCartoonPreviewTelemetry {
    pub protocol: String,
    pub flank_bp: usize,
    pub overlap_bp: usize,
    pub insert_bp: usize,
    pub bindings: ProtocolCartoonTemplateBindings,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Structured annotation projection telemetry emitted by genomic extraction
/// operations (`ExtractGenomeRegion`, `ExtractGenomeGene`).
pub struct GenomeAnnotationProjectionTelemetry {
    pub requested_scope: String,
    pub effective_scope: String,
    pub max_features_cap: Option<usize>,
    pub candidate_feature_count: usize,
    pub attached_feature_count: usize,
    pub dropped_feature_count: usize,
    pub genes_attached: usize,
    pub transcripts_attached: usize,
    pub exons_attached: usize,
    pub cds_attached: usize,
    pub fallback_applied: bool,
    pub fallback_reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Progress payload emitted by TFBS annotation operations.
pub struct TfbsProgress {
    pub seq_id: String,
    pub motif_id: String,
    pub motif_index: usize,
    pub motif_count: usize,
    pub scanned_steps: usize,
    pub total_steps: usize,
    pub motif_percent: f64,
    pub total_percent: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Progress payload emitted by genome track import operations.
pub struct GenomeTrackImportProgress {
    pub seq_id: String,
    pub source: String,
    pub path: String,
    pub parsed_records: usize,
    pub imported_features: usize,
    pub skipped_records: usize,
    pub done: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Union of long-running operation progress events.
pub enum OperationProgress {
    Tfbs(TfbsProgress),
    GenomePrepare(PrepareGenomeProgress),
    GenomeTrackImport(GenomeTrackImportProgress),
    RnaReadInterpret(RnaReadInterpretProgress),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Immutable operation journal row.
pub struct OperationRecord {
    pub run_id: RunId,
    pub op: Operation,
    pub result: OpResult,
}

/// Engine capability snapshot used by adapters for discovery/negotiation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Capabilities {
    pub protocol_version: String,
    pub supported_operations: Vec<String>,
    pub supported_export_formats: Vec<String>,
    pub deterministic_operation_log: bool,
}

/// Compact sequence row used by state-summary style adapter surfaces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineSequenceSummary {
    pub id: String,
    pub name: Option<String>,
    pub length: usize,
    pub circular: bool,
}

/// Compact container row used by shell/CLI inspection surfaces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineContainerSummary {
    pub id: String,
    pub kind: String,
    pub member_count: usize,
    pub members: Vec<String>,
}

/// Compact arrangement row used by shell/CLI inspection surfaces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineArrangementSummary {
    pub id: String,
    pub mode: String,
    pub lane_count: usize,
    pub lane_container_ids: Vec<String>,
    pub ladders: Vec<String>,
}

/// Machine-readable snapshot of top-level engine state counts and summaries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineStateSummary {
    pub sequence_count: usize,
    pub sequences: Vec<EngineSequenceSummary>,
    pub container_count: usize,
    pub containers: Vec<EngineContainerSummary>,
    pub arrangement_count: usize,
    pub arrangements: Vec<EngineArrangementSummary>,
    pub display: DisplaySettings,
}
