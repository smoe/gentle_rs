//! Source-coherent genomic structure presentation, not transcript or TSS selection.

use serde::{Deserialize, Serialize};

pub const SCHEMA: &str = "gentle.transcript_structure_presentation.v1";
pub const NON_CLAIMS: &str = "Source designations and exact coordinate agreement do not establish a preferred, consensus or biologically correct TSS. Structure equivalence is not mature-cDNA or protein equivalence. Omitted payload records are unassessed, not missing.";

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum TranscriptProvider {
    Ensembl,
    RefSeq,
}

/// Exact full-chain membership in the supplied annotations, never biological absence.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TranscriptSourceMembership {
    EnsemblOnly,
    RefSeqOnly,
    Both,
    EnsemblOtherUnassessed,
    RefSeqOtherUnassessed,
}

/// An exon-chain comparison keeps CDS alternatives and exact source records separate.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptSourceComparisonRow {
    pub exon_chain_id: String,
    pub membership: TranscriptSourceMembership,
    pub ensembl_transcript_ids: Vec<String>,
    pub refseq_transcript_ids: Vec<String>,
    pub member_record_ids: Vec<String>,
    pub structure_ids: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TranscriptAnnotationFormat {
    Gff3,
    EnsemblGeneEntry,
}

/// Explicit source identity. No source is inferred from an accession prefix.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TranscriptAnnotationSource {
    pub path: String,
    pub sha256: String,
    pub provider: TranscriptProvider,
    pub format: TranscriptAnnotationFormat,
    pub assembly: String,
    pub release: String,
    pub accession: String,
    /// Exact sequence identifier used in this annotation, not a guessed chr alias.
    pub chromosome: String,
    /// Hash of the loaded locus sequence to which this request is attached.
    pub locus_sequence_sha256: String,
    /// GFF3 gene IDs (including any provider prefix). Ensembl entry: its gene ID.
    pub gene_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptSourceBinding {
    pub source_id: String,
    pub provider: TranscriptProvider,
    pub assembly: String,
    pub release: String,
    pub accession: String,
    pub chromosome: String,
    pub annotation_sha256: String,
    pub locus_sequence_sha256: String,
}

/// Source fields are retained verbatim alongside their human-facing designation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct TranscriptDesignation {
    pub label: String,
    pub field: String,
    pub value: String,
}

/// Inclusive genomic coordinates; array order is transcription order.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct TranscriptInterval {
    pub start_1based: u64,
    pub end_1based: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptExon {
    pub interval: TranscriptInterval,
    pub source_exon_id: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct TranscriptCds {
    pub interval: TranscriptInterval,
    pub phase: Option<u8>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SourceTranscriptStructure {
    pub source_id: String,
    pub transcript_id: String,
    pub label: String,
    pub strand: i8,
    pub exons: Vec<TranscriptExon>,
    /// None = not supplied/assessed, Some([]) = explicitly annotated noncoding.
    pub cds: Option<Vec<TranscriptCds>>,
    pub designations: Vec<TranscriptDesignation>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct PhysicalExon {
    pub exon_id: String,
    pub interval: TranscriptInterval,
    pub strand: i8,
    pub member_record_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptStructureGroup {
    pub structure_id: String,
    pub exon_chain_id: String,
    pub cds_geometry_id: String,
    pub exon_ids: Vec<String>,
    pub cds: Option<Vec<TranscriptCds>>,
    pub member_record_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptPresentationRecord {
    pub record_id: String,
    pub content_sha256: String,
    pub structure: SourceTranscriptStructure,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SourceTssTick {
    pub provider: TranscriptProvider,
    pub genomic_position_1based: u64,
    pub strand: i8,
    pub member_record_ids: Vec<String>,
    /// Other provider has a start at this exact coordinate, within supplied sources.
    pub exact_cross_source_agreement: bool,
}

/// All supplied Ensembl x RefSeq starts, not a nearest/preferred-TSS assignment.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SourceTssDelta {
    pub ensembl_position_1based: u64,
    pub refseq_position_1based: u64,
    pub strand: i8,
    /// (RefSeq - Ensembl) * strand, in transcript orientation.
    pub transcript_oriented_delta_bp: i64,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptPayloadCoverage {
    pub requested_transcript_ids: Vec<String>,
    pub included_transcript_ids: Vec<String>,
    pub unassessed_transcript_ids: Vec<String>,
    pub statement: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TranscriptStructurePresentation {
    pub schema: String,
    pub content_sha256: String,
    pub assembly: String,
    pub chromosome: String,
    pub locus_sequence_sha256: String,
    pub sources: Vec<TranscriptSourceBinding>,
    pub records: Vec<TranscriptPresentationRecord>,
    pub physical_exons: Vec<PhysicalExon>,
    pub structure_groups: Vec<TranscriptStructureGroup>,
    pub tss_ticks: Vec<SourceTssTick>,
    pub tss_deltas: Vec<SourceTssDelta>,
    pub non_claims: String,
}
