//! Expert-view data contracts for feature-centric deep-inspection UIs.

use serde::{Deserialize, Serialize};

pub const TFBS_EXPERT_INSTRUCTION: &str = "TFBS expert view: each column is one PSSM position. Bar height is information content (2 - entropy in bits) from column base frequencies. Colored segments show A/C/G/T relative frequencies; the black polyline marks the matched base across positions.";

pub const RESTRICTION_EXPERT_INSTRUCTION: &str = "Restriction-site expert view: top strand is 5'->3', bottom strand is complementary 3'->5'. The vertical cut marker shows cleavage position for the selected enzyme/site.";

pub const SPLICING_EXPERT_INSTRUCTION: &str = "Splicing expert view: one lane per transcript on a shared genomic axis. Exon geometry is coordinate-true (labels never resize exon/intron footprints). Donor/acceptor splice boundaries are marked, junction arcs summarize support across transcripts, and the transcript-vs-exon matrix shows isoform differences.";

pub const ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION: &str = "Isoform architecture view: top panel shows transcript/exon or transcript/CDS structure on genomic coordinates with 5'->3' orientation left-to-right (strand-aware axis), bottom panel shows per-isoform protein-domain architecture on amino-acid coordinates. Row order is shared across both panels; CDS-to-protein guide lines indicate which coding segments contribute to which amino-acid spans.";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FeatureExpertTarget {
    TfbsFeature {
        feature_id: usize,
    },
    RestrictionSite {
        cut_pos_1based: usize,
        #[serde(default)]
        enzyme: Option<String>,
        #[serde(default)]
        recognition_start_1based: Option<usize>,
        #[serde(default)]
        recognition_end_1based: Option<usize>,
    },
    SplicingFeature {
        feature_id: usize,
    },
    IsoformArchitecture {
        panel_id: String,
    },
}

impl FeatureExpertTarget {
    pub fn describe(&self) -> String {
        match self {
            Self::TfbsFeature { feature_id } => format!("tfbs feature #{feature_id}"),
            Self::RestrictionSite {
                cut_pos_1based,
                enzyme,
                recognition_start_1based,
                recognition_end_1based,
            } => {
                let enzyme = enzyme
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("*");
                let mut out = format!("restriction cut@{cut_pos_1based} enzyme={enzyme}");
                if let (Some(start), Some(end)) = (recognition_start_1based, recognition_end_1based)
                {
                    out.push_str(&format!(" range={start}..{end}"));
                }
                out
            }
            Self::SplicingFeature { feature_id } => {
                format!("splicing feature #{feature_id}")
            }
            Self::IsoformArchitecture { panel_id } => {
                format!("isoform architecture panel '{panel_id}'")
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingRange {
    pub start_1based: usize,
    pub end_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingExonSummary {
    pub start_1based: usize,
    pub end_1based: usize,
    pub support_transcript_count: usize,
    pub constitutive: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingBoundaryMarker {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub side: String,
    pub position_1based: usize,
    pub motif_2bp: String,
    pub canonical: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingJunctionArc {
    pub donor_1based: usize,
    pub acceptor_1based: usize,
    pub support_transcript_count: usize,
    pub transcript_feature_ids: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingTranscriptLane {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub label: String,
    pub strand: String,
    pub exons: Vec<SplicingRange>,
    pub introns: Vec<SplicingRange>,
    pub has_target_feature: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingMatrixRow {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub label: String,
    pub exon_presence: Vec<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingEventSummary {
    pub event_type: String,
    pub count: usize,
    pub details: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingExpertView {
    pub seq_id: String,
    pub target_feature_id: usize,
    pub group_label: String,
    pub strand: String,
    pub region_start_1based: usize,
    pub region_end_1based: usize,
    pub transcript_count: usize,
    pub unique_exon_count: usize,
    pub instruction: String,
    pub transcripts: Vec<SplicingTranscriptLane>,
    pub unique_exons: Vec<SplicingExonSummary>,
    pub matrix_rows: Vec<SplicingMatrixRow>,
    pub boundaries: Vec<SplicingBoundaryMarker>,
    pub junctions: Vec<SplicingJunctionArc>,
    pub events: Vec<SplicingEventSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformArchitectureTranscriptLane {
    pub isoform_id: String,
    pub label: String,
    #[serde(default)]
    pub transcript_id: Option<String>,
    #[serde(default)]
    pub transcript_feature_id: Option<usize>,
    pub strand: String,
    #[serde(default)]
    pub transcript_exons: Vec<SplicingRange>,
    pub exons: Vec<SplicingRange>,
    pub introns: Vec<SplicingRange>,
    pub mapped: bool,
    #[serde(default)]
    pub transactivation_class: Option<String>,
    #[serde(default)]
    pub cds_to_protein_segments: Vec<IsoformArchitectureCdsAaSegment>,
    #[serde(default)]
    pub note: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformArchitectureCdsAaSegment {
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub aa_start: usize,
    pub aa_end: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformArchitectureProteinDomain {
    pub name: String,
    pub start_aa: usize,
    pub end_aa: usize,
    #[serde(default)]
    pub color_hex: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformArchitectureProteinLane {
    pub isoform_id: String,
    pub label: String,
    #[serde(default)]
    pub transcript_id: Option<String>,
    #[serde(default)]
    pub expected_length_aa: Option<usize>,
    #[serde(default)]
    pub reference_start_aa: Option<usize>,
    #[serde(default)]
    pub reference_end_aa: Option<usize>,
    pub domains: Vec<IsoformArchitectureProteinDomain>,
    #[serde(default)]
    pub transactivation_class: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformArchitectureExpertView {
    pub seq_id: String,
    pub panel_id: String,
    pub gene_symbol: String,
    #[serde(default = "default_isoform_transcript_geometry_mode")]
    pub transcript_geometry_mode: String,
    #[serde(default)]
    pub panel_source: Option<String>,
    pub region_start_1based: usize,
    pub region_end_1based: usize,
    pub instruction: String,
    pub transcript_lanes: Vec<IsoformArchitectureTranscriptLane>,
    pub protein_lanes: Vec<IsoformArchitectureProteinLane>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

fn default_isoform_transcript_geometry_mode() -> String {
    "exon".to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TfbsExpertColumn {
    pub index_1based: usize,
    pub counts: [f64; 4],
    pub frequencies: [f64; 4],
    pub information_content_bits: f64,
    pub match_base: Option<char>,
    pub match_frequency: Option<f64>,
    pub llr_bits: Option<f64>,
    pub true_log_odds_bits: Option<f64>,
    pub llr_rank_desc: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TfbsExpertView {
    pub seq_id: String,
    pub feature_id: usize,
    pub feature_label: String,
    pub tf_id: String,
    #[serde(default)]
    pub tf_name: Option<String>,
    pub strand: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub motif_length: usize,
    pub matched_sequence: String,
    #[serde(default)]
    pub llr_total_bits: Option<f64>,
    #[serde(default)]
    pub llr_quantile: Option<f64>,
    #[serde(default)]
    pub true_log_odds_total_bits: Option<f64>,
    #[serde(default)]
    pub true_log_odds_quantile: Option<f64>,
    pub instruction: String,
    pub columns: Vec<TfbsExpertColumn>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RestrictionSiteExpertView {
    pub seq_id: String,
    pub cut_pos_1based: usize,
    pub recognition_start_1based: usize,
    pub recognition_end_1based: usize,
    pub cut_index_0based: usize,
    pub number_of_cuts_for_enzyme: usize,
    #[serde(default)]
    pub selected_enzyme: Option<String>,
    pub enzyme_names: Vec<String>,
    #[serde(default)]
    pub recognition_iupac: Option<String>,
    pub site_sequence: String,
    pub site_sequence_complement: String,
    #[serde(default)]
    pub overlap_bp: Option<isize>,
    pub instruction: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", content = "data", rename_all = "snake_case")]
pub enum FeatureExpertView {
    Tfbs(TfbsExpertView),
    RestrictionSite(RestrictionSiteExpertView),
    Splicing(SplicingExpertView),
    IsoformArchitecture(IsoformArchitectureExpertView),
}

impl FeatureExpertView {
    pub fn instruction(&self) -> &str {
        match self {
            Self::Tfbs(_) => TFBS_EXPERT_INSTRUCTION,
            Self::RestrictionSite(_) => RESTRICTION_EXPERT_INSTRUCTION,
            Self::Splicing(_) => SPLICING_EXPERT_INSTRUCTION,
            Self::IsoformArchitecture(_) => ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION,
        }
    }
}
