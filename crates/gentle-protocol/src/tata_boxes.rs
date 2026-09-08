//! TATA evidence keeps source annotations, EPD promoter classifications, and
//! sequence-model predictions separate. All local intervals are zero-based,
//! half-open; a TSS is the position of its first transcribed base.

use serde::{Deserialize, Serialize};

pub const TATA_BOX_SCREEN_SCHEMA: &str = "gentle.tata_box_screen.v1";
pub const TATA_BOX_NON_CLAIM: &str = "Source annotation is not necessarily experimental validation. EPD TATA status is a FindM promoter classification, not an exact site. TBP matrix hits are sequence predictions, not evidence of occupancy or promoter activity; not detected is not biological absence.";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
pub struct TataBoxScreenRequest {
    pub seq_id: String,
    pub start_0based: usize,
    pub end_0based_exclusive: Option<usize>,
    pub include_annotations: bool,
    pub predict: bool,
    /// Explicit opt-in to scan both strands without restricting to known TSSs.
    pub scan_without_tss: bool,
    pub motif_id: String,
    pub minimum_llr_bits: f64,
    /// Inclusive signed distance of the motif's transcript-oriented first base.
    pub minimum_tss_offset: i32,
    pub maximum_tss_offset: i32,
    pub additional_tss: Vec<TataBoxTss>,
    pub epd: Option<TataBoxEpdSource>,
    pub max_rows: usize,
}

impl Default for TataBoxScreenRequest {
    fn default() -> Self {
        Self {
            seq_id: String::new(),
            start_0based: 0,
            end_0based_exclusive: None,
            include_annotations: true,
            predict: true,
            scan_without_tss: false,
            motif_id: "MA0108.3".into(),
            minimum_llr_bits: 6.0,
            minimum_tss_offset: -40,
            maximum_tss_offset: -15,
            additional_tss: vec![],
            epd: None,
            max_rows: 5000,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct TataBoxTss {
    pub id: String,
    pub position_0based: usize,
    pub reverse: bool,
    pub source: String,
}

/// EPD's published BED8 and promoter_motifs.txt, joined by promoter ID.
/// Files are local-only; reads never trigger a download. Identity can be pinned.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct TataBoxEpdSource {
    pub bed_path: String,
    pub motifs_path: String,
    pub assembly: String,
    pub taxon_id: u64,
    pub release: String,
    pub source_url: String,
    #[serde(default)]
    pub required: bool,
    #[serde(default)]
    pub expected_bed_sha256: Option<String>,
    #[serde(default)]
    pub expected_motifs_sha256: Option<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum TataBoxEvidenceKind {
    SourceAnnotation,
    EpdClassification,
    MotifPrediction,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TataBoxTssAssociation {
    pub tss_id: String,
    pub source: String,
    pub signed_distance_bp: i64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TataBoxEvidenceRow {
    pub row_id: String,
    pub evidence_kind: TataBoxEvidenceKind,
    pub label: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub reverse: bool,
    /// EPD rows mark the TSS only, never an inferred exact TATA interval.
    pub geometry_kind: String,
    pub source_feature_id: Option<usize>,
    pub source_qualifiers: Vec<(String, Option<String>)>,
    pub tata_positive: Option<bool>,
    pub llr_bits: Option<f64>,
    pub sequence_5prime_to_3prime: Option<String>,
    pub tss_associations: Vec<TataBoxTssAssociation>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TataBoxScreenReport {
    pub schema: String,
    pub report_id: String,
    /// Hash of this report with this field empty.
    pub content_sha256: String,
    pub request: TataBoxScreenRequest,
    pub sequence_sha256: String,
    pub annotation_sha256: String,
    pub anchor_sha256: String,
    pub genome_id: Option<String>,
    pub chromosome: Option<String>,
    pub genomic_start_1based: Option<usize>,
    pub genomic_end_1based: Option<usize>,
    pub genomic_reverse: Option<bool>,
    pub motif_id: Option<String>,
    pub matrix_sha256: Option<String>,
    pub score_policy: String,
    pub epd_status: String,
    pub epd_bed_sha256: Option<String>,
    pub epd_motifs_sha256: Option<String>,
    pub tss: Vec<TataBoxTss>,
    pub rows: Vec<TataBoxEvidenceRow>,
    pub scored_windows: usize,
    pub ambiguous_windows: usize,
    pub warnings: Vec<String>,
    pub non_claim: String,
}

/// Materialization recomputes the screen and requires the exact reviewed digest.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct TataBoxMaterializeRequest {
    pub screen: TataBoxScreenRequest,
    pub expected_report_sha256: String,
    pub row_ids: Vec<String>,
}
