//! Expert-view helpers for feature-centric deep-inspection UIs.
//!
//! Most portable expert-view payloads now live in `gentle-protocol`; this
//! module re-exports them for compatibility and keeps only helper logic that
//! computes derived matrices from those shared records.

pub use gentle_protocol::{
    FeatureExpertTarget, FeatureExpertView, ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION,
    IsoformArchitectureCdsAaSegment, IsoformArchitectureExpertView,
    IsoformArchitectureProteinDomain, IsoformArchitectureProteinLane,
    IsoformArchitectureTranscriptLane, RESTRICTION_EXPERT_INSTRUCTION, RestrictionSiteExpertView,
    SPLICING_EXPERT_INSTRUCTION, SplicingBoundaryMarker, SplicingEventSummary,
    SplicingExonCdsPhase, SplicingExonSummary, SplicingExpertView, SplicingJunctionArc,
    SplicingMatrixRow, SplicingRange, SplicingScopePreset, SplicingTranscriptLane,
    TFBS_EXPERT_INSTRUCTION, TfbsExpertColumn, TfbsExpertView,
};
use std::collections::HashMap;

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
