//! Expert-view data contracts for feature-centric deep-inspection UIs.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

pub const TFBS_EXPERT_INSTRUCTION: &str = "TFBS expert view: each column is one PSSM position. Bar height is information content (2 - entropy in bits) from column base frequencies. Colored segments show A/C/G/T relative frequencies; the black polyline marks the matched base across positions.";

pub const RESTRICTION_EXPERT_INSTRUCTION: &str = "Restriction-site expert view: top strand is 5'->3', bottom strand is complementary 3'->5'. The vertical cut marker shows cleavage position for the selected enzyme/site.";

pub const SPLICING_EXPERT_INSTRUCTION: &str = "Splicing expert view: one lane per transcript on a shared genomic axis. Exon geometry is coordinate-true (labels never resize exon/intron footprints). Donor/acceptor splice boundaries are marked, junction arcs summarize support across transcripts, and the transcript-vs-exon matrix shows isoform differences.";

pub const ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION: &str = "Isoform architecture view: top panel shows transcript/exon or transcript/CDS structure on genomic coordinates with 5'->3' orientation left-to-right (strand-aware axis), bottom panel shows per-isoform protein-domain architecture on amino-acid coordinates. Row order is shared across both panels; CDS-to-protein guide lines indicate which coding segments contribute to which amino-acid spans.";

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SplicingScopePreset {
    #[default]
    AllOverlappingBothStrands,
    TargetGroupAnyStrand,
    AllOverlappingTargetStrand,
    TargetGroupTargetStrand,
}

impl SplicingScopePreset {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AllOverlappingBothStrands => "all_overlapping_both_strands",
            Self::TargetGroupAnyStrand => "target_group_any_strand",
            Self::AllOverlappingTargetStrand => "all_overlapping_target_strand",
            Self::TargetGroupTargetStrand => "target_group_target_strand",
        }
    }

    pub fn restrict_to_target_group(self) -> bool {
        matches!(
            self,
            Self::TargetGroupAnyStrand | Self::TargetGroupTargetStrand
        )
    }

    pub fn restrict_to_target_strand(self) -> bool {
        matches!(
            self,
            Self::AllOverlappingTargetStrand | Self::TargetGroupTargetStrand
        )
    }
}

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
        #[serde(default)]
        scope: SplicingScopePreset,
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
            Self::SplicingFeature { feature_id, scope } => {
                format!("splicing feature #{feature_id} scope={}", scope.as_str())
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
pub struct SplicingExonCdsPhase {
    pub start_1based: usize,
    pub end_1based: usize,
    #[serde(default)]
    pub left_cds_phase: Option<u8>,
    #[serde(default)]
    pub right_cds_phase: Option<u8>,
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
    #[serde(default)]
    pub exon_cds_phases: Vec<SplicingExonCdsPhase>,
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

#[derive(Debug, Clone)]
pub struct SplicingExonTransitionMatrix {
    pub counts: Vec<Vec<usize>>,
    pub transcript_feature_ids: Vec<Vec<Vec<usize>>>,
}

pub fn compute_splicing_exon_transition_matrix(
    view: &SplicingExpertView,
) -> SplicingExonTransitionMatrix {
    let exon_count = view.unique_exons.len();
    let mut counts = vec![vec![0usize; exon_count]; exon_count];
    let mut transcript_feature_ids = vec![vec![Vec::<usize>::new(); exon_count]; exon_count];
    if exon_count == 0 {
        return SplicingExonTransitionMatrix {
            counts,
            transcript_feature_ids,
        };
    }
    let exon_index: HashMap<(usize, usize), usize> = view
        .unique_exons
        .iter()
        .enumerate()
        .map(|(idx, exon)| ((exon.start_1based, exon.end_1based), idx))
        .collect();
    for transcript in &view.transcripts {
        let mut ordered_indices = transcript
            .exons
            .iter()
            .filter_map(|exon| {
                exon_index
                    .get(&(exon.start_1based, exon.end_1based))
                    .copied()
            })
            .collect::<Vec<_>>();
        if transcript.strand.trim() == "-" {
            ordered_indices.reverse();
        }
        for pair in ordered_indices.windows(2) {
            let from = pair[0];
            let to = pair[1];
            if from == to {
                continue;
            }
            counts[from][to] += 1;
            transcript_feature_ids[from][to].push(transcript.transcript_feature_id);
        }
    }
    for row in &mut transcript_feature_ids {
        for participants in row {
            participants.sort_unstable();
            participants.dedup();
        }
    }
    SplicingExonTransitionMatrix {
        counts,
        transcript_feature_ids,
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn range(start_1based: usize, end_1based: usize) -> SplicingRange {
        SplicingRange {
            start_1based,
            end_1based,
        }
    }

    #[test]
    fn compute_splicing_exon_transition_matrix_respects_transcript_strand() {
        let exon_a = SplicingExonSummary {
            start_1based: 100,
            end_1based: 150,
            support_transcript_count: 3,
            constitutive: true,
        };
        let exon_b = SplicingExonSummary {
            start_1based: 200,
            end_1based: 250,
            support_transcript_count: 2,
            constitutive: false,
        };
        let exon_c = SplicingExonSummary {
            start_1based: 300,
            end_1based: 330,
            support_transcript_count: 3,
            constitutive: true,
        };
        let view = SplicingExpertView {
            seq_id: "s".to_string(),
            target_feature_id: 1,
            group_label: "GENE1".to_string(),
            strand: "+".to_string(),
            region_start_1based: 100,
            region_end_1based: 330,
            transcript_count: 3,
            unique_exon_count: 3,
            instruction: "splicing".to_string(),
            transcripts: vec![
                SplicingTranscriptLane {
                    transcript_feature_id: 11,
                    transcript_id: "tx_plus_1".to_string(),
                    label: "tx_plus_1".to_string(),
                    strand: "+".to_string(),
                    exons: vec![range(100, 150), range(200, 250), range(300, 330)],
                    exon_cds_phases: vec![],
                    introns: vec![range(151, 199), range(251, 299)],
                    has_target_feature: true,
                },
                SplicingTranscriptLane {
                    transcript_feature_id: 12,
                    transcript_id: "tx_plus_2".to_string(),
                    label: "tx_plus_2".to_string(),
                    strand: "+".to_string(),
                    exons: vec![range(100, 150), range(300, 330)],
                    exon_cds_phases: vec![],
                    introns: vec![range(151, 299)],
                    has_target_feature: false,
                },
                SplicingTranscriptLane {
                    transcript_feature_id: 13,
                    transcript_id: "tx_minus".to_string(),
                    label: "tx_minus".to_string(),
                    strand: "-".to_string(),
                    exons: vec![range(100, 150), range(200, 250), range(300, 330)],
                    exon_cds_phases: vec![],
                    introns: vec![range(151, 199), range(251, 299)],
                    has_target_feature: false,
                },
            ],
            unique_exons: vec![exon_a, exon_b, exon_c],
            matrix_rows: vec![],
            boundaries: vec![],
            junctions: vec![],
            events: vec![],
        };
        let matrix = compute_splicing_exon_transition_matrix(&view);
        assert_eq!(matrix.counts.len(), 3);
        assert_eq!(matrix.counts[0][1], 1, "E1->E2 from tx_plus_1");
        assert_eq!(matrix.counts[1][2], 1, "E2->E3 from tx_plus_1");
        assert_eq!(matrix.counts[0][2], 1, "E1->E3 from tx_plus_2");
        assert_eq!(matrix.counts[2][1], 1, "E3->E2 from tx_minus");
        assert_eq!(matrix.counts[1][0], 1, "E2->E1 from tx_minus");
        assert_eq!(matrix.transcript_feature_ids[0][2], vec![12]);
        assert_eq!(matrix.transcript_feature_ids[2][1], vec![13]);
    }
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
    pub enzyme_cut_offset_0based: Option<isize>,
    #[serde(default)]
    pub overlap_bp: Option<isize>,
    #[serde(default)]
    pub enzyme_note: Option<String>,
    #[serde(default)]
    pub rebase_url: Option<String>,
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
