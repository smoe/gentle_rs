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
/// Small offline resource for cDNA/EST support linked to exon/junction geometry.
pub const CDNA_EST_EVIDENCE_RESOURCE_SCHEMA: &str = "gentle.cdna_est_evidence_resource.v1";
/// Human-facing interpretation boundary shared by GUI and SVG renderers.
pub const GENE_ISOFORM_EVIDENCE_INSTRUCTION: &str = "Isoform evidence inspector: transcript models and coordinate geometry are annotation-derived; RNA reads, cDNA/EST records, array probes, expression values, projected occupancy tracks, and qPCR assays are shown as separate evidence layers. Missing evidence is unknown. Occupancy is locus-level evidence; spatial overlap alone neither identifies a regulated isoform nor establishes causality.";

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
