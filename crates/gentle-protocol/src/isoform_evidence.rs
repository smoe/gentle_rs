//! Portable contracts for gene-level isoform evidence inspection.
//!
//! The inspector composes annotation geometry with optional experimental and
//! assay reports. It is deliberately an evidence ledger: missing evidence is
//! not evidence of absence, and probe overlap alone is never promoted to an
//! isoform-validation claim.

use crate::{
    EnsemblRegulationSourceDescriptor, GenomicRegionEvidenceAvailability, GenomicRegionPurpose,
    GenomicRegionSelectionMethod, GenomicRegionStrand, SplicingExpertView,
};
use serde::{Deserialize, Serialize};

/// Legacy machine-readable report produced before per-measurement evidence retention.
pub const GENE_ISOFORM_EVIDENCE_SCHEMA_V1: &str = "gentle.gene_isoform_evidence.v1";
/// Machine-readable report produced by the gene isoform evidence inspector.
pub const GENE_ISOFORM_EVIDENCE_SCHEMA: &str = "gentle.gene_isoform_evidence.v2";
/// Publication-oriented composition of isoform, occupancy, motif, and assay evidence.
pub const GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA: &str = "gentle.gene_locus_evidence_display.v1";
/// Declarative grouping and scaling metadata for projected occupancy tracks.
pub const GENE_LOCUS_OCCUPANCY_LAYOUT_SCHEMA: &str = "gentle.gene_locus_occupancy_layout.v1";
/// Portable offline input for externally computed, coordinate-bound TF scores.
pub const GENE_LOCUS_EXTERNAL_REGULATORY_SCORE_SCHEMA: &str =
    "gentle.gene_locus_external_regulatory_scores.v1";
/// Optional, source-bound Ensembl Regulation layer in a locus display report.
pub const GENE_LOCUS_ENSEMBL_REGULATION_EVIDENCE_SCHEMA: &str =
    "gentle.gene_locus_ensembl_regulation_evidence.v1";
/// Small offline resource for cDNA/EST support linked to exon/junction geometry.
pub const CDNA_EST_EVIDENCE_RESOURCE_SCHEMA: &str = "gentle.cdna_est_evidence_resource.v1";
/// Human-facing interpretation boundary shared by GUI and SVG renderers.
pub const GENE_ISOFORM_EVIDENCE_INSTRUCTION: &str = "Isoform evidence inspector: transcript models and coordinate geometry are annotation-derived; RNA reads, cDNA/EST records, array probes, expression values, projected occupancy tracks, and qPCR assays are shown as separate evidence layers. Missing evidence is unknown. Occupancy is locus-level evidence; spatial overlap alone neither identifies a regulated isoform nor establishes causality.";
/// Human-facing interpretation boundary for composed locus figures.
pub const GENE_LOCUS_EVIDENCE_DISPLAY_INSTRUCTION: &str = "Gene locus evidence display: transcript/CDS geometry, predicted TF binding scores, projected occupancy, chromatin context, and validation-assay candidates are aligned on one strand-aware axis. These are distinct evidence layers. CUT&RUN enrichment is occupancy evidence, not biochemical affinity; chromatin context does not establish TF binding; PWM or model scores are predictions, not observed binding. Shared visual position or scale does not establish isoform-specific regulation, causal binding, or cross-sample comparability.";

/// Whether a missing optional Ensembl Regulation snapshot should abort locus preparation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusEnsemblRegulationAvailabilityPolicy {
    /// Preserve a typed unavailable layer so reports remain honest and portable.
    #[default]
    RetainUnavailable,
    /// Refuse preparation when the requested source files are unavailable.
    Required,
}

/// Availability state of an explicitly requested Ensembl Regulation layer.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusEnsemblRegulationAvailability {
    Available,
    #[default]
    Unavailable,
}

/// Bounded, offline-after-install request for Ensembl Regulation locus evidence.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusEnsemblRegulationRequest {
    pub source_id: String,
    pub interval_index_path: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub intervals_path: Option<String>,
    /// Optional exact content lock for the selected sparse index JSON.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expected_index_sha256: Option<String>,
    /// Optional exact content lock for the normalized interval TSV.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expected_intervals_sha256: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub feature_types: Vec<String>,
    pub max_rows: usize,
    pub availability_policy: GeneLocusEnsemblRegulationAvailabilityPolicy,
}

impl Default for GeneLocusEnsemblRegulationRequest {
    fn default() -> Self {
        Self {
            source_id: String::new(),
            interval_index_path: String::new(),
            intervals_path: None,
            expected_index_sha256: None,
            expected_intervals_sha256: None,
            feature_types: vec![],
            max_rows: 10_000,
            availability_policy: GeneLocusEnsemblRegulationAvailabilityPolicy::default(),
        }
    }
}

/// Content and release identity retained for a locus-level Ensembl query.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusEnsemblRegulationSourceBinding {
    pub source: EnsemblRegulationSourceDescriptor,
    pub overlap_report_id: String,
    pub index_sha256: String,
    pub intervals_sha256: String,
    pub content_identity_verified: bool,
    pub assembly_match_status: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub interval_index_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub resolved_intervals_path: Option<String>,
    pub requested_feature_types: Vec<String>,
    pub max_rows: usize,
    pub matched_feature_count: usize,
    pub returned_feature_count: usize,
    pub truncated: bool,
}

/// One source-bound Ensembl regulatory feature projected onto a locus figure.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusEnsemblRegulationFeatureRow {
    pub source_id: String,
    pub provider: String,
    pub annotation_release: String,
    pub annotation_api_version: String,
    pub pipeline_version: String,
    pub assembly_name: String,
    pub assembly_accession: String,
    pub feature_id: String,
    pub feature_type: String,
    pub core_genomic_start_1based: usize,
    pub core_genomic_end_1based: usize,
    pub displayed_genomic_start_1based: usize,
    pub displayed_genomic_end_1based: usize,
    pub displayed_local_start_1based: usize,
    pub displayed_local_end_1based: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub extended_genomic_start_1based: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub extended_genomic_end_1based: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub extended_local_start_1based: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub extended_local_end_1based: Option<usize>,
    pub local_strand: String,
    pub genomic_strand: String,
    pub clipped: bool,
    pub locus_relation: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub associated_gene_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub associated_gene_names: Vec<String>,
    pub canonical_feature_url: String,
    pub feature_url_template: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub overlapping_tss_class_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub nearest_tss_class_ids: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub signed_distance_to_nearest_tss_bp: Option<i64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub absolute_distance_to_nearest_tss_bp: Option<u64>,
    pub signed_distance_basis: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub overlapping_reporter_architecture_ids: Vec<String>,
    pub relationship_statement: String,
}

/// Optional Ensembl Regulation layer carried by the portable locus report.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusEnsemblRegulationEvidence {
    pub schema: String,
    pub availability: GeneLocusEnsemblRegulationAvailability,
    pub requested_source_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source: Option<EnsemblRegulationSourceDescriptor>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unavailable_reason: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_binding: Option<GeneLocusEnsemblRegulationSourceBinding>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub rows: Vec<GeneLocusEnsemblRegulationFeatureRow>,
    pub evidence_statement: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub non_claims: Vec<String>,
    pub provider_link_note: String,
}

/// References to optional evidence sources composed by a pure-read inspection.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneIsoformEvidenceRequest {
    pub panel_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub annotation_release: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub rna_read_report_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub qpcr_report_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub probe_evidence_paths: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub cdna_est_resource_paths: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expression_tsv_path: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub occupancy_track_names: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusOccupancyScaleMode {
    /// One scale within each biological/sample group.
    #[default]
    SharedGroup,
    /// One scale across all groups; requires an explicit comparability rationale.
    SharedAll,
    /// Each lane uses its own maximum.
    Independent,
    /// Every lane in the group uses `fixed_abs_max_score`.
    Fixed,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusOccupancyLaneRole {
    #[default]
    Experimental,
    GfpControl,
    InputControl,
    IggControl,
    PositiveControl,
    NegativeControl,
    /// Promoter/chromatin-state context such as H3K4me3, not TF occupancy.
    ChromatinContext,
    Other,
}

/// Whether one requested occupancy/chromatin lane could be evaluated.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusOccupancyLaneState {
    #[default]
    Available,
    NotPrepared,
    NoCompatibleInterval,
    AssemblyMismatch,
}

impl GeneLocusOccupancyLaneState {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Available => "available",
            Self::NotPrepared => "not_prepared",
            Self::NoCompatibleInterval => "no_compatible_interval",
            Self::AssemblyMismatch => "assembly_mismatch",
        }
    }
}

/// One projected track and its explicitly declared figure role.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusOccupancyLaneRequest {
    pub track_name: String,
    /// Stable source identity independent of a local display name or path.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_sha256: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_assembly: Option<String>,
    /// Optional preparation outcome supplied by an engine-owned composer.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_state: Option<GeneLocusOccupancyLaneState>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub display_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub condition_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_line_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub batch_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assay: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mark: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub factor: Option<String>,
    pub role: GeneLocusOccupancyLaneRole,
}

/// One cell-line/sample group with an explicit scaling policy.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusOccupancyGroupRequest {
    pub group_id: String,
    pub label: String,
    pub scale_mode: GeneLocusOccupancyScaleMode,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fixed_abs_max_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cross_group_scale_justification: Option<String>,
    pub lanes: Vec<GeneLocusOccupancyLaneRequest>,
}

/// Offline layout resource used to keep lane grouping out of filename heuristics.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusOccupancyLayout {
    pub schema: String,
    pub groups: Vec<GeneLocusOccupancyGroupRequest>,
}

/// Whether and how a genomic scale bar is included in a locus figure.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusScaleBarMode {
    /// Backward-compatible default for reports created before scale bars existed.
    #[default]
    Hidden,
    /// Select a deterministic 1/2/5 x 10^n length from the displayed span.
    Auto,
    /// Use the positive `length_bp` supplied by the request.
    Fixed,
}

/// Requested genomic scale-bar policy.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusScaleBarPolicy {
    pub mode: GeneLocusScaleBarMode,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub length_bp: Option<usize>,
}

/// Resolved scale bar retained in the portable report.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusScaleBar {
    pub mode: GeneLocusScaleBarMode,
    pub length_bp: usize,
    pub label: String,
}

/// Source family for one predicted regulatory-score track.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusRegulatoryScoreProviderKind {
    #[default]
    JasparPwm,
    ExternalModelScores,
    Other,
}

/// Strand components retained from the provider output.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusRegulatoryScoreStrandPolicy {
    #[default]
    Both,
    ForwardOnly,
    ReverseOnly,
}

/// Display-scale policy for one predicted score track.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusRegulatoryScoreScaleMode {
    #[default]
    Independent,
    SharedGroup,
    Fixed,
}

/// Whether a normalized score track is available for rendering.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusRegulatoryScoreState {
    Available,
    #[default]
    NotAssessable,
}

/// Typed calibration status; descriptive prose is never parsed as policy.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusRegulatoryCalibrationState {
    #[default]
    Unspecified,
    MatrixSpecific,
    ProviderDeclaredUncalibrated,
    ProviderDeclaredCalibrated,
    /// Different score sources may share a scale only when an exact
    /// `calibration_id` and `calibration_sha256` bind them.
    CrossSourceCalibrated,
}

impl GeneLocusRegulatoryCalibrationState {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Unspecified => "unspecified",
            Self::MatrixSpecific => "matrix_specific",
            Self::ProviderDeclaredUncalibrated => "provider_declared_uncalibrated",
            Self::ProviderDeclaredCalibrated => "provider_declared_calibrated",
            Self::CrossSourceCalibrated => "cross_source_calibrated",
        }
    }
}

/// Explicit factor identity carried independently of display labels and paths.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusRegulatoryFactor {
    pub factor_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub factor_label: Option<String>,
}

/// Factors bound to one exact matrix/model source in a multi-source request.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusRegulatorySourceFactorBinding {
    pub source_id: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub factors: Vec<GeneLocusRegulatoryFactor>,
}

/// One independently configured predicted regulatory-score request.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct GeneLocusRegulatoryScoreTrackRequest {
    pub track_id: String,
    pub label: String,
    pub provider_kind: GeneLocusRegulatoryScoreProviderKind,
    /// JASPAR matrix/factor tokens in deterministic request order.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub source_ids: Vec<String>,
    /// Required for `external_model_scores`; ignored by the local JASPAR adapter.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_path: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub factors: Vec<GeneLocusRegulatoryFactor>,
    /// Preferred for requests resolving more than one matrix/model source.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub source_factor_bindings: Vec<GeneLocusRegulatorySourceFactorBinding>,
    pub score_kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub score_units: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub calibration_statement: Option<String>,
    pub calibration_state: GeneLocusRegulatoryCalibrationState,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub calibration_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub calibration_sha256: Option<String>,
    pub strand_policy: GeneLocusRegulatoryScoreStrandPolicy,
    pub clip_negative: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub display_threshold: Option<f64>,
    pub top_hit_count: usize,
    pub scale_mode: GeneLocusRegulatoryScoreScaleMode,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scale_group: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub shared_scale_justification: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fixed_scale_min: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fixed_scale_max: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub color_hint: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub line_style_hint: Option<String>,
}

impl Default for GeneLocusRegulatoryScoreTrackRequest {
    fn default() -> Self {
        Self {
            track_id: String::new(),
            label: String::new(),
            provider_kind: GeneLocusRegulatoryScoreProviderKind::JasparPwm,
            source_ids: vec![],
            source_path: None,
            factors: vec![],
            source_factor_bindings: vec![],
            score_kind: default_locus_motif_score_kind(),
            score_units: None,
            calibration_statement: None,
            calibration_state: GeneLocusRegulatoryCalibrationState::Unspecified,
            calibration_id: None,
            calibration_sha256: None,
            strand_policy: GeneLocusRegulatoryScoreStrandPolicy::Both,
            clip_negative: true,
            display_threshold: None,
            top_hit_count: default_locus_motif_top_hit_count(),
            scale_mode: GeneLocusRegulatoryScoreScaleMode::Independent,
            scale_group: None,
            shared_scale_justification: None,
            fixed_scale_min: None,
            fixed_scale_max: None,
            color_hint: None,
            line_style_hint: None,
        }
    }
}

fn default_locus_upstream_bp() -> usize {
    5_000
}

fn default_locus_downstream_bp() -> usize {
    1_000
}

fn default_locus_motif_score_kind() -> String {
    "llr_background_tail_log10".to_string()
}

fn default_locus_motif_top_hit_count() -> usize {
    5
}

/// Request for a publication-oriented locus composition over one imported panel.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct GeneLocusEvidenceDisplayRequest {
    pub isoform_evidence: GeneIsoformEvidenceRequest,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub probe_effect_table_paths: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub probe_effect_contrasts: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub probe_effect_coordinate_system: Option<String>,
    pub occupancy_layout: GeneLocusOccupancyLayout,
    pub motifs: Vec<String>,
    pub motif_score_kind: String,
    pub motif_clip_negative: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub motif_display_threshold: Option<f64>,
    pub motif_top_hit_count: usize,
    /// Additive provider-neutral requests. Legacy `motifs` are normalized into
    /// equivalent local-JASPAR requests by the engine.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub regulatory_score_tracks: Vec<GeneLocusRegulatoryScoreTrackRequest>,
    pub scale_bar: GeneLocusScaleBarPolicy,
    /// Optional persisted genomic-region sets projected onto this report's
    /// canonical locus axis.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub region_set_ids: Vec<String>,
    /// Preserve local source paths in the portable report. Disabled by default
    /// so reports can be shared without exposing workstation paths.
    pub include_local_source_paths: bool,
}

impl Default for GeneLocusEvidenceDisplayRequest {
    fn default() -> Self {
        Self {
            isoform_evidence: GeneIsoformEvidenceRequest::default(),
            upstream_bp: default_locus_upstream_bp(),
            downstream_bp: default_locus_downstream_bp(),
            probe_effect_table_paths: vec![],
            probe_effect_contrasts: vec![],
            probe_effect_coordinate_system: None,
            occupancy_layout: GeneLocusOccupancyLayout::default(),
            motifs: vec![],
            motif_score_kind: default_locus_motif_score_kind(),
            motif_clip_negative: true,
            motif_display_threshold: None,
            motif_top_hit_count: default_locus_motif_top_hit_count(),
            regulatory_score_tracks: vec![],
            scale_bar: GeneLocusScaleBarPolicy::default(),
            region_set_ids: vec![],
            include_local_source_paths: false,
        }
    }
}

/// One saved genomic ROI projected onto the canonical locus display axis.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusSavedRegionOverlayRow {
    pub set_id: String,
    pub set_content_sha256: String,
    pub region_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub purpose: GenomicRegionPurpose,
    pub selection_method: GenomicRegionSelectionMethod,
    pub evidence_availability: GenomicRegionEvidenceAvailability,
    pub local_start_1based: usize,
    pub local_end_1based: usize,
    pub genomic_start_0based: u64,
    pub genomic_end_0based_exclusive: u64,
    pub genomic_strand: GenomicRegionStrand,
    pub region_content_sha256: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub evidence_ids: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum CdnaEstEvidenceKind {
    Cdna,
    Est,
    CuratedTranscript,
    #[default]
    Other,
}

/// One independently sourced sequence-evidence record.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct CdnaEstEvidenceRecord {
    pub evidence_id: String,
    pub kind: CdnaEstEvidenceKind,
    pub source: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub accession: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assembly: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genomic_start_1based: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genomic_end_1based: Option<usize>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub exon_family_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub junction_ids: Vec<String>,
    #[serde(default = "default_support_count")]
    pub support_count: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub alignment_identity_fraction: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub alignment_coverage_fraction: Option<f64>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
}

fn default_support_count() -> usize {
    1
}

/// Offline cDNA/EST evidence with explicit provenance and geometry references.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct CdnaEstEvidenceResource {
    pub schema: String,
    pub resource_id: String,
    pub source: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub annotation_release: Option<String>,
    pub records: Vec<CdnaEstEvidenceRecord>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum IsoformEvidenceSourceKind {
    #[default]
    AnnotationModel,
    RnaRead,
    Cdna,
    Est,
    CuratedTranscript,
    OtherSequence,
    ArrayProbe,
    Expression,
    OccupancyTrack,
    QpcrAssay,
}

/// One projected interval retained on a gene-level occupancy lane.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneIsoformOccupancyInterval {
    pub interval_id: String,
    pub local_start_1based: usize,
    pub local_end_1based: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
}

/// One source-specific BED/BigWig lane aligned to the transcript models.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneIsoformOccupancyLane {
    pub lane_id: String,
    pub track_name: String,
    pub display_label: String,
    pub source_kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_path: Option<String>,
    pub interval_count: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_score: Option<f64>,
    pub intervals: Vec<GeneIsoformOccupancyInterval>,
}

/// One occupancy lane after applying its declared role and display scale.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusOccupancyLane {
    pub lane: GeneIsoformOccupancyLane,
    pub state: GeneLocusOccupancyLaneState,
    pub source_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_sha256: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_assembly: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub display_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub condition_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cell_line_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub batch_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assay: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mark: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub factor: Option<String>,
    pub role: GeneLocusOccupancyLaneRole,
    pub display_abs_max_score: f64,
}

/// One explicitly separated cell-line/sample block in the locus figure.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusOccupancyGroup {
    pub group_id: String,
    pub label: String,
    pub scale_mode: GeneLocusOccupancyScaleMode,
    pub group_abs_max_score: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cross_group_scale_justification: Option<String>,
    pub lanes: Vec<GeneLocusOccupancyLane>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusCodonKind {
    #[default]
    Start,
    Stop,
}

/// Direction in which loaded-sequence coordinates increase on the rendered
/// locus axis.
///
/// This is intentionally independent of the genomic gene strand: an imported
/// negative-strand Ensembl locus may already be reverse-complemented into
/// gene-oriented local coordinates.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusLocalAxisDirection {
    /// Compatibility behavior for reports written before the local-axis field
    /// existed. Renderers fall back to the historical genomic-strand rule.
    #[default]
    LegacyGeneStrandFallback,
    IncreasingLeftToRight,
    DecreasingLeftToRight,
}

impl GeneLocusLocalAxisDirection {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::LegacyGeneStrandFallback => "legacy_gene_strand_fallback",
            Self::IncreasingLeftToRight => "increasing_left_to_right",
            Self::DecreasingLeftToRight => "decreasing_left_to_right",
        }
    }
}

/// One evidence-backed translation boundary drawn on one transcript lane.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusCodonMarker {
    pub transcript_id: String,
    pub kind: GeneLocusCodonKind,
    pub local_position_1based: usize,
    pub genomic_position_1based: usize,
    /// Legacy local-sequence strand retained for backward compatibility.
    pub strand: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub local_strand: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub genomic_strand: String,
    pub basis: String,
}

/// Annotation-derived metrics shown beside one transcript model.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusTranscriptMetrics {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub spliced_exon_length_bp: usize,
    pub cds_length_bp: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expected_peptide_length_aa: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub protein_identity_sha256: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub predicted_molecular_weight_kda: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub protein_mass_basis: Option<String>,
    pub coding_status: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub biotype: Option<String>,
    pub cds_ranges_local_1based: Vec<(usize, usize)>,
    pub flags: Vec<String>,
}

/// One junction-spanning validation marker, deduplicated across transcript rows.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusAssayOverlay {
    pub junction_id: String,
    pub local_donor_1based: usize,
    pub local_acceptor_1based: usize,
    pub genomic_donor_1based: usize,
    pub genomic_acceptor_1based: usize,
    pub assay_ids: Vec<String>,
    pub family_ids: Vec<String>,
    pub transcript_ids: Vec<String>,
}

/// Array feature class retained by the locus probe-effect overlay.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneLocusProbeClass {
    Psr,
    Juc,
    #[default]
    Other,
}

impl GeneLocusProbeClass {
    /// Classify a probe row using explicit type metadata before stable feature-id
    /// and geometry fallbacks. The returned basis is persisted for inspection.
    pub fn classify(
        raw_type_or_level: &str,
        feature_id: &str,
        has_junction_edges: bool,
    ) -> (Self, String) {
        let raw = raw_type_or_level.trim().to_ascii_lowercase();
        if raw == "psr"
            || raw == "probeset"
            || raw == "exon"
            || raw.contains("->psr")
            || raw.contains("probe_selection_region")
        {
            return (Self::Psr, "explicit_probe_type_or_level".to_string());
        }
        if raw == "juc" || raw == "junction" || raw.contains("->juc") {
            return (Self::Juc, "explicit_probe_type_or_level".to_string());
        }
        let feature_upper = feature_id.trim().to_ascii_uppercase();
        if feature_upper.starts_with("PSR") {
            return (Self::Psr, "clariom_feature_id_prefix".to_string());
        }
        if feature_upper.starts_with("JUC") {
            return (Self::Juc, "clariom_feature_id_prefix".to_string());
        }
        if has_junction_edges {
            return (Self::Juc, "junction_edge_columns".to_string());
        }
        (Self::Other, "unclassified".to_string())
    }
}

/// One ordered contrast column exposed by a probe-effect source table.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusProbeEffectContrast {
    pub contrast_id: String,
    pub display_label: String,
    pub source_column: String,
}

/// One effect value associated with one probe geometry.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusProbeEffectValue {
    pub contrast_id: String,
    pub display_label: String,
    pub source_column: String,
    pub value: f64,
    pub unit: String,
}

/// One PSR interval or JUC span with its ordered contrast effects.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusProbeEffectOverlay {
    pub feature_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub parent_feature_id: Option<String>,
    pub probe_class: GeneLocusProbeClass,
    pub classification_basis: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub coordinate_system: String,
    pub chromosome: String,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub local_start_1based: usize,
    pub local_end_1based: usize,
    pub strand: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub transcript_cluster_id: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub exon_ids: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub junction_start_edge_1based: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub junction_stop_edge_1based: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pm_probe_count: Option<usize>,
    pub effects: Vec<GeneLocusProbeEffectValue>,
    pub source_path: String,
    pub source_row_1based: usize,
    pub mapping_status: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
}

/// One labeled high-scoring position on a continuous motif track.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusMotifHit {
    pub rank: usize,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    /// Legacy local-sequence strand retained for backward compatibility.
    pub strand: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub local_strand: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub genomic_strand: String,
    pub score: f64,
}

/// One existing TFBS score track projected into the locus composition contract.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusMotifTrack {
    pub motif_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub motif_name: Option<String>,
    pub motif_length_bp: usize,
    pub score_kind: String,
    pub track_start_0based: usize,
    pub max_score: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub display_threshold: Option<f64>,
    pub forward_scores: Vec<f64>,
    pub reverse_scores: Vec<f64>,
    pub top_hits: Vec<GeneLocusMotifHit>,
    pub provenance: String,
}

/// One labeled site call retained by a normalized regulatory-score track.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusRegulatoryScoreSite {
    pub site_id: String,
    pub rank: usize,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    /// Legacy strand field retained for backward compatibility. New reports
    /// must use the explicit local/genomic fields below.
    pub strand: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub local_strand: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub genomic_strand: String,
    pub score: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
}

/// Provider-neutral score-track result consumed directly by renderers.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusRegulatoryScoreTrack {
    pub track_id: String,
    pub label: String,
    pub provider_kind: GeneLocusRegulatoryScoreProviderKind,
    pub provider_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub provider_version: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub factors: Vec<GeneLocusRegulatoryFactor>,
    pub source_ids: Vec<String>,
    pub input_sequence_id: String,
    pub input_sequence_sha256: String,
    pub assembly: String,
    pub chromosome: String,
    pub anchor_start_1based: usize,
    pub anchor_end_1based: usize,
    pub coordinate_convention: String,
    pub window_length_bp: usize,
    pub stride_bp: usize,
    pub strand_policy: GeneLocusRegulatoryScoreStrandPolicy,
    pub score_kind: String,
    pub score_units: String,
    pub score_directionality: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub theoretical_min: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub theoretical_max: Option<f64>,
    pub calibration_status: String,
    pub calibration_state: GeneLocusRegulatoryCalibrationState,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub calibration_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub calibration_sha256: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub calibration_statement: Option<String>,
    pub track_start_0based: usize,
    pub forward_scores: Vec<f64>,
    pub reverse_scores: Vec<f64>,
    pub sites: Vec<GeneLocusRegulatoryScoreSite>,
    pub display_threshold: Option<f64>,
    pub scale_mode: GeneLocusRegulatoryScoreScaleMode,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scale_group: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub shared_scale_justification: Option<String>,
    pub display_scale_min: f64,
    pub display_scale_max: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub color_hint: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub line_style_hint: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_sha256: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub request_sha256: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_sha256: Option<String>,
    pub provenance: String,
    pub state: GeneLocusRegulatoryScoreState,
    pub warnings: Vec<String>,
}

/// One portable external score payload. It is deliberately data-only and does
/// not identify or invoke a networked model provider.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusExternalRegulatoryScoreResource {
    pub schema: String,
    pub track_id: String,
    pub label: String,
    pub model_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub model_version: Option<String>,
    pub factors: Vec<GeneLocusRegulatoryFactor>,
    pub input_sequence_id: String,
    pub input_sequence_sha256: String,
    pub assembly: String,
    pub chromosome: String,
    pub anchor_start_1based: usize,
    pub anchor_end_1based: usize,
    pub coordinate_convention: String,
    pub window_length_bp: usize,
    pub stride_bp: usize,
    pub strand_policy: GeneLocusRegulatoryScoreStrandPolicy,
    pub score_kind: String,
    pub score_units: String,
    pub score_directionality: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub theoretical_min: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub theoretical_max: Option<f64>,
    pub calibration_state: GeneLocusRegulatoryCalibrationState,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub calibration_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub calibration_sha256: Option<String>,
    pub calibration_statement: String,
    pub track_start_0based: usize,
    pub forward_scores: Vec<f64>,
    pub reverse_scores: Vec<f64>,
    pub sites: Vec<GeneLocusRegulatoryScoreSite>,
    pub request_schema_sha256: String,
    pub request_sha256: String,
    pub score_payload_sha256: String,
    pub provenance: String,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum IsoformEvidenceAssessmentStatus {
    Observed,
    Candidate,
    ConstraintOnly,
    #[default]
    NotEvaluated,
    Unknown,
}

impl IsoformEvidenceAssessmentStatus {
    /// Stable human-facing label shared by evidence-ledger frontends.
    pub const fn display_label(self) -> &'static str {
        match self {
            Self::Observed => "observed evidence",
            Self::Candidate => "candidate association",
            Self::ConstraintOnly => "design constraint",
            Self::NotEvaluated => "not evaluated",
            Self::Unknown => "unresolved evidence",
        }
    }

    /// Conservative interpretation boundary for the corresponding status.
    pub const fn interpretation(self) -> &'static str {
        match self {
            Self::Observed => {
                "A source record or measurement was attached to this geometry; this does not by itself establish causality or isoform-specific regulation."
            }
            Self::Candidate => {
                "A computed or projected association is available for review and is not a validation claim."
            }
            Self::ConstraintOnly => {
                "This evidence constrains assay or feature geometry but does not demonstrate biological support."
            }
            Self::NotEvaluated => {
                "The applicable analysis was not run or no comparable source was supplied."
            }
            Self::Unknown => {
                "A source was supplied, but its relationship to this geometry could not be resolved confidently."
            }
        }
    }
}

/// One evidence statement attached to one or more geometry objects or families.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneIsoformEvidenceItem {
    pub evidence_id: String,
    pub source_kind: IsoformEvidenceSourceKind,
    pub source_id: String,
    pub status: IsoformEvidenceAssessmentStatus,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub target_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub family_ids: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub support_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub support_fraction: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub value: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub probe_class: Option<GeneLocusProbeClass>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub probe_classification_basis: Option<String>,
    pub method: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
}

/// One source measurement retained without cross-condition or cross-unit collapse.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneIsoformEvidenceMeasurement {
    pub evidence_id: String,
    pub status: IsoformEvidenceAssessmentStatus,
    pub classification: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub value: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    pub explanation: String,
}

/// One independently interpretable scoring dimension.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneIsoformEvidenceComponent {
    pub status: IsoformEvidenceAssessmentStatus,
    pub classification: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub value: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub evidence_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub measurements: Vec<GeneIsoformEvidenceMeasurement>,
    pub explanation: String,
}

/// The four intentionally separate dimensions used for assay triage.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneIsoformEvidenceComponents {
    pub specificity: GeneIsoformEvidenceComponent,
    pub abundance: GeneIsoformEvidenceComponent,
    pub responsiveness: GeneIsoformEvidenceComponent,
    pub assayability: GeneIsoformEvidenceComponent,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneIsoformRecommendationTier {
    AssayReady,
    EvidencePrioritized,
    AnnotationCandidate,
    #[default]
    NotEvaluated,
}

/// Rule-based triage guidance. This is deliberately not a weighted evidence score.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneIsoformRecommendation {
    pub tier: GeneIsoformRecommendationTier,
    pub recommended_use: String,
    pub rule: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub evidence_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneIsoformFamilyRow {
    pub family_id: String,
    pub label: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub transcript_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub aliases: Vec<String>,
    pub mapped_transcript_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneIsoformTranscriptRow {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub label: String,
    pub strand: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub family_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub exon_family_ids_5_to_3: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub exon_family_ids_genomic_ascending: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneIsoformExonFamilyRow {
    pub exon_family_id: String,
    pub coordinate_frame: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub local_start_1based: usize,
    pub local_end_1based: usize,
    pub strand: String,
    pub annotation_model_count: usize,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub transcript_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub family_ids: Vec<String>,
    pub specificity_class: String,
    pub components: GeneIsoformEvidenceComponents,
    pub recommendation: GeneIsoformRecommendation,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneIsoformJunctionRow {
    pub junction_id: String,
    pub coordinate_frame: String,
    pub genomic_low_1based: usize,
    pub genomic_high_1based: usize,
    pub local_low_1based: usize,
    pub local_high_1based: usize,
    pub transcript_donor_1based: usize,
    pub transcript_acceptor_1based: usize,
    pub strand: String,
    pub annotation_model_count: usize,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub transcript_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub family_ids: Vec<String>,
    pub specificity_class: String,
    pub components: GeneIsoformEvidenceComponents,
    pub recommendation: GeneIsoformRecommendation,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneIsoformAssayCandidate {
    pub qpcr_report_id: String,
    pub assay_rank: usize,
    pub score: f64,
    pub assay_class_label: String,
    pub explanation: String,
    pub amplicon_length_bp: usize,
    pub forward_sequence_5_to_3: String,
    pub reverse_sequence_5_to_3: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub probe_sequence_5_to_3: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub target_junction_ids: Vec<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub supported_transcript_ids: Vec<String>,
    pub satisfies_requested_targeting: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneIsoformEvidenceProvenanceSource {
    pub source_kind: String,
    pub source_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub schema: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
}

/// Deterministic gene-level ledger joining annotation, evidence, and assays.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneIsoformEvidenceReport {
    pub schema: String,
    pub seq_id: String,
    pub gene_symbol: String,
    pub panel_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub annotation_release: Option<String>,
    pub assembly: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,
    pub coordinate_frame: String,
    pub gene_strand: String,
    pub region_start_1based: usize,
    pub region_end_1based: usize,
    pub instruction: String,
    pub splicing: Option<SplicingExpertView>,
    pub families: Vec<GeneIsoformFamilyRow>,
    pub transcripts: Vec<GeneIsoformTranscriptRow>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub transcript_metrics: Vec<GeneLocusTranscriptMetrics>,
    pub exon_families: Vec<GeneIsoformExonFamilyRow>,
    pub junctions: Vec<GeneIsoformJunctionRow>,
    pub occupancy_lanes: Vec<GeneIsoformOccupancyLane>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub occupancy_shared_abs_max_score: Option<f64>,
    pub evidence_items: Vec<GeneIsoformEvidenceItem>,
    pub assay_candidates: Vec<GeneIsoformAssayCandidate>,
    pub provenance: Vec<GeneIsoformEvidenceProvenanceSource>,
    pub warnings: Vec<String>,
}

/// Deterministic publication composition for one gene locus.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GeneLocusEvidenceDisplayReport {
    pub schema: String,
    pub seq_id: String,
    pub gene_symbol: String,
    pub panel_id: String,
    pub instruction: String,
    /// Genomic strand of the interpreted gene.
    pub gene_strand: String,
    /// Presentation direction of loaded-sequence coordinates. This prevents a
    /// gene-oriented negative-strand import from being mirrored a second time.
    pub local_axis_direction: GeneLocusLocalAxisDirection,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub locus_local_start_1based: usize,
    pub locus_local_end_1based: usize,
    pub locus_genomic_start_1based: usize,
    pub locus_genomic_end_1based: usize,
    pub axis_left_genomic_1based: usize,
    pub axis_right_genomic_1based: usize,
    pub isoform_evidence: GeneIsoformEvidenceReport,
    pub transcript_metrics: Vec<GeneLocusTranscriptMetrics>,
    pub codon_markers: Vec<GeneLocusCodonMarker>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub probe_effect_contrasts: Vec<GeneLocusProbeEffectContrast>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub probe_effect_overlays: Vec<GeneLocusProbeEffectOverlay>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub probe_effect_shared_abs_max: Option<f64>,
    pub occupancy_groups: Vec<GeneLocusOccupancyGroup>,
    pub motif_tracks: Vec<GeneLocusMotifTrack>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub regulatory_score_tracks: Vec<GeneLocusRegulatoryScoreTrack>,
    pub scale_bar: GeneLocusScaleBar,
    pub assay_overlays: Vec<GeneLocusAssayOverlay>,
    /// Present only when Ensembl Regulation evidence was explicitly requested.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ensembl_regulation: Option<GeneLocusEnsemblRegulationEvidence>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub saved_region_overlays: Vec<GeneLocusSavedRegionOverlayRow>,
    pub provenance: Vec<GeneIsoformEvidenceProvenanceSource>,
    pub warnings: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isoform_evidence_status_language_keeps_observation_distinct_from_validation() {
        assert_eq!(
            IsoformEvidenceAssessmentStatus::Observed.display_label(),
            "observed evidence"
        );
        assert!(
            IsoformEvidenceAssessmentStatus::Observed
                .interpretation()
                .contains("does not by itself establish causality")
        );
        assert!(
            IsoformEvidenceAssessmentStatus::ConstraintOnly
                .interpretation()
                .contains("does not demonstrate biological support")
        );
        assert_ne!(
            IsoformEvidenceAssessmentStatus::NotEvaluated.display_label(),
            IsoformEvidenceAssessmentStatus::Unknown.display_label()
        );
    }

    #[test]
    fn pre_probe_overlay_locus_payloads_keep_empty_defaults() {
        let request: GeneLocusEvidenceDisplayRequest = serde_json::from_value(serde_json::json!({
            "isoform_evidence": { "panel_id": "patz1" },
            "upstream_bp": 5000,
            "downstream_bp": 1000
        }))
        .expect("deserialize pre-overlay locus request");
        assert!(request.probe_effect_table_paths.is_empty());
        assert!(request.probe_effect_contrasts.is_empty());
        assert_eq!(request.probe_effect_coordinate_system, None);
        assert!(request.regulatory_score_tracks.is_empty());
        assert_eq!(request.scale_bar.mode, GeneLocusScaleBarMode::Hidden);
        assert_eq!(request.scale_bar.length_bp, None);
        assert!(!request.include_local_source_paths);

        let report: GeneLocusEvidenceDisplayReport = serde_json::from_value(serde_json::json!({
            "schema": GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA,
            "seq_id": "patz1",
            "gene_symbol": "PATZ1"
        }))
        .expect("deserialize pre-overlay locus report");
        assert!(report.probe_effect_contrasts.is_empty());
        assert!(report.probe_effect_overlays.is_empty());
        assert_eq!(report.probe_effect_shared_abs_max, None);
        assert!(report.regulatory_score_tracks.is_empty());
        assert_eq!(report.scale_bar.mode, GeneLocusScaleBarMode::Hidden);
        assert_eq!(report.scale_bar.length_bp, 0);
        assert_eq!(
            report.local_axis_direction,
            GeneLocusLocalAxisDirection::LegacyGeneStrandFallback
        );
        assert!(report.ensembl_regulation.is_none());
    }

    #[test]
    fn ensembl_regulation_locus_contract_round_trips_without_changing_report_schema() {
        let evidence = GeneLocusEnsemblRegulationEvidence {
            schema: GENE_LOCUS_ENSEMBL_REGULATION_EVIDENCE_SCHEMA.to_string(),
            availability: GeneLocusEnsemblRegulationAvailability::Available,
            requested_source_id: "ensembl_regulation_2026_08_grch38".to_string(),
            source: Some(EnsemblRegulationSourceDescriptor {
                source_id: "ensembl_regulation_2026_08_grch38".to_string(),
                provider: "Ensembl Regulation".to_string(),
                annotation_release: "2026-08".to_string(),
                browser_species_slug: "homo_sapiens".to_string(),
                feature_page_url_template:
                    "https://regulation.ensembl.org/regulatory_features/homo_sapiens/{FEATURE_ID}"
                        .to_string(),
                ..Default::default()
            }),
            rows: vec![GeneLocusEnsemblRegulationFeatureRow {
                source_id: "ensembl_regulation_2026_08_grch38".to_string(),
                feature_id: "ENSR1_958".to_string(),
                feature_type: "promoter".to_string(),
                core_genomic_start_1based: 100,
                core_genomic_end_1based: 180,
                displayed_local_start_1based: 11,
                displayed_local_end_1based: 91,
                canonical_feature_url:
                    "https://regulation.ensembl.org/regulatory_features/homo_sapiens/ENSR1_958"
                        .to_string(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let report = GeneLocusEvidenceDisplayReport {
            schema: GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA.to_string(),
            ensembl_regulation: Some(evidence),
            ..Default::default()
        };
        let value = serde_json::to_value(&report).expect("serialize locus report");
        assert_eq!(value["schema"], GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA);
        assert_eq!(
            value["ensembl_regulation"]["rows"][0]["feature_id"],
            "ENSR1_958"
        );
        let decoded: GeneLocusEvidenceDisplayReport =
            serde_json::from_value(value).expect("deserialize locus report");
        assert_eq!(
            decoded.ensembl_regulation.expect("Ensembl evidence").rows[0].feature_id,
            "ENSR1_958"
        );

        let request: GeneLocusEnsemblRegulationRequest =
            serde_json::from_value(serde_json::json!({
                "source_id": "ensembl_regulation_2026_08_grch38",
                "interval_index_path": "regulation.index.json"
            }))
            .expect("deserialize Ensembl locus request");
        assert_eq!(request.max_rows, 10_000);
        assert_eq!(
            request.availability_policy,
            GeneLocusEnsemblRegulationAvailabilityPolicy::RetainUnavailable
        );
    }

    #[test]
    fn locus_axis_direction_and_explicit_strands_round_trip() {
        let report = GeneLocusEvidenceDisplayReport {
            local_axis_direction: GeneLocusLocalAxisDirection::IncreasingLeftToRight,
            codon_markers: vec![GeneLocusCodonMarker {
                strand: "+".to_string(),
                local_strand: "+".to_string(),
                genomic_strand: "-".to_string(),
                ..Default::default()
            }],
            ..Default::default()
        };
        let encoded = serde_json::to_value(&report).expect("serialize explicit axis semantics");
        assert_eq!(encoded["local_axis_direction"], "increasing_left_to_right");
        assert_eq!(encoded["codon_markers"][0]["local_strand"], "+");
        assert_eq!(encoded["codon_markers"][0]["genomic_strand"], "-");
        let decoded: GeneLocusEvidenceDisplayReport =
            serde_json::from_value(encoded).expect("deserialize explicit axis semantics");
        assert_eq!(decoded.local_axis_direction, report.local_axis_direction);
        assert_eq!(decoded.codon_markers[0].genomic_strand, "-");
    }

    #[test]
    fn chromatin_context_lane_and_score_track_settings_round_trip() {
        let request: GeneLocusEvidenceDisplayRequest = serde_json::from_value(serde_json::json!({
            "isoform_evidence": { "panel_id": "demo" },
            "occupancy_layout": {
                "schema": GENE_LOCUS_OCCUPANCY_LAYOUT_SCHEMA,
                "groups": [{
                    "group_id": "chromatin",
                    "label": "Chromatin context",
                    "scale_mode": "independent",
                    "lanes": [{
                        "track_name": "h3k4me3",
                        "source_id": "synthetic:h3k4me3",
                        "source_sha256": "sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                        "source_assembly": "GRCh38",
                        "source_state": "not_prepared",
                        "condition_label": "untreated",
                        "cell_line_label": "Saos-2",
                        "batch_label": "batch-1",
                        "assay": "CUT&RUN",
                        "mark": "H3K4me3",
                        "role": "chromatin_context"
                    }]
                }]
            },
            "regulatory_score_tracks": [{
                "track_id": "tp73_pwm",
                "label": "TP73 PWM",
                "provider_kind": "jaspar_pwm",
                "source_ids": ["MA0861.1"],
                "factors": [{"factor_id": "TP73", "factor_label": "p73"}],
                "source_factor_bindings": [{
                    "source_id": "MA0861.1",
                    "factors": [{"factor_id": "TP73", "factor_label": "p73"}]
                }],
                "score_kind": "llr_bits",
                "calibration_state": "matrix_specific",
                "strand_policy": "both",
                "clip_negative": false,
                "display_threshold": 2.5,
                "top_hit_count": 7,
                "scale_mode": "independent"
            }],
            "scale_bar": {"mode": "fixed", "length_bp": 1000},
            "include_local_source_paths": true
        }))
        .expect("deserialize additive gene-locus request");

        let lane = &request.occupancy_layout.groups[0].lanes[0];
        assert_eq!(lane.role, GeneLocusOccupancyLaneRole::ChromatinContext);
        assert_eq!(lane.cell_line_label.as_deref(), Some("Saos-2"));
        assert_eq!(lane.batch_label.as_deref(), Some("batch-1"));
        assert_eq!(lane.source_id.as_deref(), Some("synthetic:h3k4me3"));
        assert_eq!(
            lane.source_state,
            Some(GeneLocusOccupancyLaneState::NotPrepared)
        );
        assert_eq!(lane.assay.as_deref(), Some("CUT&RUN"));
        assert_eq!(lane.mark.as_deref(), Some("H3K4me3"));
        assert_eq!(request.regulatory_score_tracks.len(), 1);
        assert_eq!(
            request.regulatory_score_tracks[0].provider_kind,
            GeneLocusRegulatoryScoreProviderKind::JasparPwm
        );
        assert_eq!(request.regulatory_score_tracks[0].top_hit_count, 7);
        assert_eq!(
            request.regulatory_score_tracks[0].calibration_state,
            GeneLocusRegulatoryCalibrationState::MatrixSpecific
        );
        assert_eq!(
            request.regulatory_score_tracks[0].source_factor_bindings[0].source_id,
            "MA0861.1"
        );
        assert_eq!(request.scale_bar.mode, GeneLocusScaleBarMode::Fixed);
        assert_eq!(request.scale_bar.length_bp, Some(1000));
        assert!(request.include_local_source_paths);

        let encoded = serde_json::to_value(&request).expect("serialize gene-locus request");
        assert_eq!(
            encoded["occupancy_layout"]["groups"][0]["lanes"][0]["role"],
            "chromatin_context"
        );
        assert_eq!(encoded["scale_bar"]["length_bp"], 1000);
        assert_eq!(encoded["include_local_source_paths"], true);
    }

    #[test]
    fn pre_availability_locus_lanes_remain_available_when_state_is_absent() {
        let lane: GeneLocusOccupancyLane = serde_json::from_value(serde_json::json!({
            "lane": {
                "lane_id": "OCC:legacy",
                "track_name": "legacy projected track",
                "intervals": [{
                    "interval_id": "legacy-1",
                    "local_start_1based": 1,
                    "local_end_1based": 4
                }]
            },
            "role": "experimental"
        }))
        .expect("deserialize pre-availability occupancy lane");
        assert_eq!(lane.state, GeneLocusOccupancyLaneState::Available);

        let score_request: GeneLocusRegulatoryScoreTrackRequest =
            serde_json::from_value(serde_json::json!({
                "track_id": "legacy_jaspar",
                "provider_kind": "jaspar_pwm",
                "source_ids": ["MA0861.1"]
            }))
            .expect("deserialize pre-calibration score request");
        assert_eq!(
            score_request.calibration_state,
            GeneLocusRegulatoryCalibrationState::Unspecified
        );
        assert!(score_request.source_factor_bindings.is_empty());
    }

    #[test]
    fn isoform_evidence_v1_payload_defaults_v2_additions() {
        let report: GeneIsoformEvidenceReport = serde_json::from_value(serde_json::json!({
            "schema": GENE_ISOFORM_EVIDENCE_SCHEMA_V1,
            "seq_id": "synthetic",
            "panel_id": "panel",
            "exon_families": [{
                "exon_family_id": "EXF:test",
                "components": {
                    "responsiveness": {
                        "status": "candidate",
                        "classification": "array_logfc_candidate",
                        "value": -1.2,
                        "unit": "logFC"
                    }
                }
            }]
        }))
        .expect("deserialize legacy isoform-evidence report");
        assert!(report.transcript_metrics.is_empty());
        assert!(
            report.exon_families[0]
                .components
                .responsiveness
                .measurements
                .is_empty()
        );
        assert_eq!(
            report.exon_families[0].recommendation.tier,
            GeneIsoformRecommendationTier::NotEvaluated
        );
    }

    #[test]
    fn probe_classification_prefers_explicit_metadata_then_stable_fallbacks() {
        assert_eq!(
            GeneLocusProbeClass::classify("probeset", "JUC_conflict", true),
            (
                GeneLocusProbeClass::Psr,
                "explicit_probe_type_or_level".to_string()
            )
        );
        assert_eq!(
            GeneLocusProbeClass::classify("", "JUC_demo", false).0,
            GeneLocusProbeClass::Juc
        );
        assert_eq!(
            GeneLocusProbeClass::classify("", "feature", true).0,
            GeneLocusProbeClass::Juc
        );
    }
}
