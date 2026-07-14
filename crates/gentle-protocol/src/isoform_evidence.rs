//! Portable contracts for gene-level isoform evidence inspection.
//!
//! The inspector composes annotation geometry with optional experimental and
//! assay reports. It is deliberately an evidence ledger: missing evidence is
//! not evidence of absence, and probe overlap alone is never promoted to an
//! isoform-validation claim.

use crate::SplicingExpertView;
use serde::{Deserialize, Serialize};

/// Machine-readable report produced by the gene isoform evidence inspector.
pub const GENE_ISOFORM_EVIDENCE_SCHEMA: &str = "gentle.gene_isoform_evidence.v1";
/// Publication-oriented composition of isoform, occupancy, motif, and assay evidence.
pub const GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA: &str = "gentle.gene_locus_evidence_display.v1";
/// Declarative grouping and scaling metadata for projected occupancy tracks.
pub const GENE_LOCUS_OCCUPANCY_LAYOUT_SCHEMA: &str = "gentle.gene_locus_occupancy_layout.v1";
/// Small offline resource for cDNA/EST support linked to exon/junction geometry.
pub const CDNA_EST_EVIDENCE_RESOURCE_SCHEMA: &str = "gentle.cdna_est_evidence_resource.v1";
/// Human-facing interpretation boundary shared by GUI and SVG renderers.
pub const GENE_ISOFORM_EVIDENCE_INSTRUCTION: &str = "Isoform evidence inspector: transcript models and coordinate geometry are annotation-derived; RNA reads, cDNA/EST records, array probes, expression values, projected occupancy tracks, and qPCR assays are shown as separate evidence layers. Missing evidence is unknown. Occupancy is locus-level evidence; spatial overlap alone neither identifies a regulated isoform nor establishes causality.";
/// Human-facing interpretation boundary for composed locus figures.
pub const GENE_LOCUS_EVIDENCE_DISPLAY_INSTRUCTION: &str = "Gene locus evidence display: transcript/CDS geometry, motif scores, projected occupancy, and validation-assay candidates are aligned on one strand-aware axis. These are distinct evidence layers. Shared visual position or scale does not establish isoform-specific regulation, causal binding, or cross-sample comparability.";

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
    Other,
}

/// One projected track and its explicitly declared figure role.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct GeneLocusOccupancyLaneRequest {
    pub track_name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub display_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub condition_label: Option<String>,
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
        }
    }
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
    #[serde(skip_serializing_if = "Option::is_none")]
    pub display_label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub condition_label: Option<String>,
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

/// One evidence-backed translation boundary drawn on one transcript lane.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusCodonMarker {
    pub transcript_id: String,
    pub kind: GeneLocusCodonKind,
    pub local_position_1based: usize,
    pub genomic_position_1based: usize,
    pub strand: String,
    pub basis: String,
}

/// Annotation-derived metrics shown beside one transcript model.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct GeneLocusTranscriptMetrics {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub spliced_exon_length_bp: usize,
    pub cds_length_bp: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expected_peptide_length_aa: Option<usize>,
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
    pub strand: String,
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
    pub method: String,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
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
    pub gene_strand: String,
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
    pub assay_overlays: Vec<GeneLocusAssayOverlay>,
    pub provenance: Vec<GeneIsoformEvidenceProvenanceSource>,
    pub warnings: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

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

        let report: GeneLocusEvidenceDisplayReport =
            serde_json::from_value(serde_json::json!({
                "schema": GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA,
                "seq_id": "patz1",
                "gene_symbol": "PATZ1"
            }))
            .expect("deserialize pre-overlay locus report");
        assert!(report.probe_effect_contrasts.is_empty());
        assert!(report.probe_effect_overlays.is_empty());
        assert_eq!(report.probe_effect_shared_abs_max, None);
    }
}
