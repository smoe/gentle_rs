//! Display-only grouping, exact feature-ID joins and shared splicing row geometry.
//!
//! No evidence is fetched or inferred here. Biological payload order and hashes
//! remain untouched; ambiguous legacy joins cannot establish evaluated absence.

use crate::{
    SplicingBoundaryMarker, SplicingExpertView, SplicingUniprotReferenceEvidence,
    SplicingUniprotReferenceStatus,
};
use std::collections::{BTreeMap, BTreeSet};

/// One unique oriented boundary/pair, retaining every source marker and transcript.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplicingBoundarySummary {
    pub marker_index: usize,
    pub marker_indices: Vec<usize>,
    pub transcript_feature_ids: Vec<usize>,
    pub strand: String,
}

/// Deduplicate presentation text, not biological records. Opposite strands stay separate.
pub fn boundary_summaries(view: &SplicingExpertView) -> Vec<SplicingBoundarySummary> {
    let mut summaries = BTreeMap::new();
    for (index, marker) in view.boundaries.iter().enumerate() {
        let strand = view
            .transcripts
            .iter()
            .find(|lane| lane.transcript_feature_id == marker.transcript_feature_id)
            .map(|lane| lane.strand.clone())
            .unwrap_or_else(|| "?".to_string());
        let key = (
            marker.position_1based,
            marker.partner_position_1based,
            marker.side.clone(),
            marker.motif_2bp.clone(),
            marker.paired_motif_signature.clone(),
            marker.motif_class.clone(),
            marker.canonical,
            marker.canonical_pair,
            strand.clone(),
        );
        let summary = summaries
            .entry(key)
            .or_insert_with(|| SplicingBoundarySummary {
                marker_index: index,
                marker_indices: Vec::new(),
                transcript_feature_ids: Vec::new(),
                strand,
            });
        summary.marker_indices.push(index);
        summary
            .transcript_feature_ids
            .push(marker.transcript_feature_id);
    }
    summaries
        .into_values()
        .map(|mut summary| {
            summary.transcript_feature_ids.sort_unstable();
            summary.transcript_feature_ids.dedup();
            summary
        })
        .collect()
}

/// Readable summary/hover wording shared by GUI and SVG, retaining every source marker.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplicingBoundaryPresentationRow {
    pub summary: SplicingBoundarySummary,
    pub label: String,
    pub hover_text: String,
    pub exceptional: bool,
}

/// A canonical dinucleotide can still belong to a non-canonical paired signature.
pub fn boundary_is_exceptional(marker: &SplicingBoundaryMarker) -> bool {
    !marker.canonical || !marker.canonical_pair
}

/// Build once per immutable view for cached GUI presentation or one-shot SVG export.
pub fn boundary_presentation_rows(
    view: &SplicingExpertView,
) -> Vec<SplicingBoundaryPresentationRow> {
    boundary_summaries(view)
        .into_iter()
        .map(|summary| {
            let marker = &view.boundaries[summary.marker_index];
            let exceptional = boundary_is_exceptional(marker);
            let label = format!(
                "{} / {} | {} | {}:{}{} | {} | {} transcripts",
                marker.position_1based,
                marker.partner_position_1based,
                summary.strand,
                marker.side,
                marker.motif_2bp,
                if exceptional { "*" } else { "" },
                marker.paired_motif_signature,
                summary.transcript_feature_ids.len()
            );
            let hover_text = summary
                .marker_indices
                .iter()
                .map(|&index| {
                    let marker = &view.boundaries[index];
                    format!(
                        "n-{} {}: {}",
                        marker.transcript_feature_id, marker.transcript_id, marker.annotation
                    )
                })
                .collect::<Vec<_>>()
                .join("\n");
            SplicingBoundaryPresentationRow {
                summary,
                label,
                hover_text,
                exceptional,
            }
        })
        .collect()
}

/// Diagnostics may concern gene-level evidence, not necessarily transcript ambiguity.
pub const REFERENCE_BADGE_LEGEND: &str = "UniProt badges: +N = additional exact-xref entries; ? = loaded-evidence diagnostic (see hover), not necessarily transcript ambiguity.";

/// One chart lane, joined to at most one matrix record by exact feature ID.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplicingPresentationLane {
    pub lane_index: usize,
    pub transcript_feature_id: usize,
    pub matrix_row_index: Option<usize>,
    pub status: SplicingUniprotReferenceStatus,
    pub badge: Option<String>,
}

/// A header or transcript row in the presentation-only matrix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SplicingMatrixDisplayRow {
    Header(String),
    Transcript {
        lane_index: Option<usize>,
        matrix_row_index: Option<usize>,
    },
}

/// One half-open vertical band. Header bands never resolve to a transcript.
#[derive(Debug, Clone, PartialEq)]
pub struct SplicingCanvasRow {
    pub top: f32,
    pub bottom: f32,
    pub lane_index: Option<usize>,
    pub group: Option<SplicingUniprotReferenceStatus>,
}

/// Shared painting and hit-testing coordinates, relative to the canvas origin.
#[derive(Debug, Clone, PartialEq)]
pub struct SplicingCanvasLayout {
    pub rows: Vec<SplicingCanvasRow>,
    pub height: f32,
}

impl SplicingCanvasLayout {
    pub fn lane_at_y(&self, y: f32) -> Option<usize> {
        if !y.is_finite() {
            return None;
        }
        self.rows
            .iter()
            .find(|row| row.top <= y && y < row.bottom)?
            .lane_index
    }

    pub fn centre_for_lane(&self, lane_index: usize) -> Option<f32> {
        self.rows
            .iter()
            .find(|row| row.lane_index == Some(lane_index))
            .map(|row| (row.top + row.bottom) * 0.5)
    }
}

/// Deterministic display order with all original records retained.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplicingPresentationLayout {
    pub lanes: Vec<SplicingPresentationLane>,
    pub orphan_matrix_rows: Vec<usize>,
    pub show_headers: bool,
    pub diagnostics: Vec<String>,
}

/// Local-evidence wording shared by GUI and SVG, never a global absence claim.
pub fn group_heading(status: SplicingUniprotReferenceStatus) -> &'static str {
    match status {
        SplicingUniprotReferenceStatus::Referenced => {
            "UniProt-referenced isoforms (loaded evidence)"
        }
        SplicingUniprotReferenceStatus::NotReferenced => {
            "Other transcripts (no exact xref in relevant loaded evidence)"
        }
        SplicingUniprotReferenceStatus::NotEvaluated => "UniProt status not evaluated",
    }
}

fn rank(status: SplicingUniprotReferenceStatus) -> u8 {
    match status {
        SplicingUniprotReferenceStatus::Referenced => 0,
        SplicingUniprotReferenceStatus::NotReferenced => 1,
        SplicingUniprotReferenceStatus::NotEvaluated => 2,
    }
}

/// Sorting only: the existing exact-xref matcher remains engine-owned.
fn transcript_sort_key(id: &str) -> String {
    let id = id.trim().to_ascii_uppercase();
    match id.rsplit_once('.') {
        Some((prefix, suffix))
            if !prefix.is_empty()
                && !suffix.is_empty()
                && suffix.bytes().all(|byte| byte.is_ascii_digit()) =>
        {
            prefix.to_string()
        }
        _ => id,
    }
}

/// A single compact badge; gene-only sources never count. Full tuples remain in hovers/JSON.
pub fn reference_badge(evidence: &SplicingUniprotReferenceEvidence) -> Option<String> {
    if evidence.status != SplicingUniprotReferenceStatus::Referenced {
        return None;
    }
    let mut matches = evidence
        .sources
        .iter()
        .filter(|source| !source.matched_transcript_xrefs.is_empty())
        .map(|source| {
            (
                source.accession.as_str(),
                source.entry_id.as_str(),
                source.reviewed,
            )
        })
        .collect::<Vec<_>>();
    matches.sort_by_key(|(accession, entry_id, _)| (*accession, *entry_id));
    let (accession, _, reviewed) = matches.first()?;
    let review = match reviewed {
        Some(true) => "reviewed",
        Some(false) => "unreviewed",
        None => "unknown",
    };
    let overflow = if matches.len() > 1 {
        format!(" +{}", matches.len() - 1)
    } else {
        String::new()
    };
    let diagnostic = if evidence.diagnostics.is_empty() {
        ""
    } else {
        " ?"
    };
    Some(format!(
        "UniProt {accession} [{review}]{overflow}{diagnostic}"
    ))
}

impl SplicingPresentationLayout {
    pub fn new(view: &SplicingExpertView) -> Self {
        let mut lane_counts = BTreeMap::<usize, usize>::new();
        let mut matrix_indices = BTreeMap::<usize, Vec<usize>>::new();
        for lane in &view.transcripts {
            *lane_counts.entry(lane.transcript_feature_id).or_default() += 1;
        }
        for (index, row) in view.matrix_rows.iter().enumerate() {
            matrix_indices
                .entry(row.transcript_feature_id)
                .or_default()
                .push(index);
        }
        let mut diagnostics = Vec::new();
        for (id, count) in &lane_counts {
            if *count > 1 {
                diagnostics.push(format!(
                    "Duplicate transcript feature ID {id}; evidence join not evaluated."
                ));
            }
        }
        for (id, indices) in &matrix_indices {
            if indices.len() > 1 {
                diagnostics.push(format!(
                    "Duplicate matrix feature ID {id}; evidence join not evaluated."
                ));
            }
        }
        let mut used_rows = BTreeSet::new();
        let mut lanes = view
            .transcripts
            .iter()
            .enumerate()
            .map(|(lane_index, lane)| {
                let matrix_row_index = matrix_indices
                    .get(&lane.transcript_feature_id)
                    .filter(|rows| rows.len() == 1 && lane_counts[&lane.transcript_feature_id] == 1)
                    .map(|rows| rows[0]);
                let evidence = matrix_row_index.map(|index| {
                    used_rows.insert(index);
                    &view.matrix_rows[index].uniprot_reference
                });
                if matrix_row_index.is_none() {
                    diagnostics.push(format!(
                        "No unambiguous matrix row for feature {}; status not evaluated.",
                        lane.transcript_feature_id
                    ));
                }
                SplicingPresentationLane {
                    lane_index,
                    transcript_feature_id: lane.transcript_feature_id,
                    matrix_row_index,
                    status: evidence.map(|e| e.status).unwrap_or_default(),
                    badge: evidence.and_then(reference_badge),
                }
            })
            .collect::<Vec<_>>();
        lanes.sort_by_key(|lane| {
            (
                rank(lane.status),
                !view.transcripts[lane.lane_index].has_target_feature,
                transcript_sort_key(&view.transcripts[lane.lane_index].transcript_id),
                lane.transcript_feature_id,
                lane.lane_index,
            )
        });
        let orphan_matrix_rows = (0..view.matrix_rows.len())
            .filter(|index| !used_rows.contains(index))
            .collect::<Vec<_>>();
        if !orphan_matrix_rows.is_empty() {
            diagnostics.push(format!(
                "{} matrix row(s) have no unambiguous lane; retained at the end.",
                orphan_matrix_rows.len()
            ));
        }
        Self {
            show_headers: lanes
                .iter()
                .any(|lane| lane.status != SplicingUniprotReferenceStatus::NotEvaluated),
            lanes,
            orphan_matrix_rows,
            diagnostics,
        }
    }

    pub fn canvas(&self, lane_height: f32, header_height: f32) -> SplicingCanvasLayout {
        let mut canvas = SplicingCanvasLayout {
            rows: Vec::new(),
            height: 0.0,
        };
        if !lane_height.is_finite()
            || lane_height <= 0.0
            || !header_height.is_finite()
            || header_height < 0.0
        {
            return canvas;
        }
        let mut previous = None;
        for lane in &self.lanes {
            if self.show_headers && previous != Some(lane.status) {
                canvas.rows.push(SplicingCanvasRow {
                    top: canvas.height,
                    bottom: canvas.height + header_height,
                    lane_index: None,
                    group: Some(lane.status),
                });
                canvas.height += header_height;
            }
            canvas.rows.push(SplicingCanvasRow {
                top: canvas.height,
                bottom: canvas.height + lane_height,
                lane_index: Some(lane.lane_index),
                group: None,
            });
            canvas.height += lane_height;
            previous = Some(lane.status);
        }
        canvas
    }

    pub fn matrix_display_rows(&self) -> Vec<SplicingMatrixDisplayRow> {
        let mut rows = Vec::new();
        let mut previous = None;
        for lane in &self.lanes {
            if self.show_headers && previous != Some(lane.status) {
                rows.push(SplicingMatrixDisplayRow::Header(
                    group_heading(lane.status).to_string(),
                ));
            }
            rows.push(SplicingMatrixDisplayRow::Transcript {
                lane_index: Some(lane.lane_index),
                matrix_row_index: lane.matrix_row_index,
            });
            previous = Some(lane.status);
        }
        if !self.orphan_matrix_rows.is_empty() {
            rows.push(SplicingMatrixDisplayRow::Header(
                "Saved matrix rows without an unambiguous transcript lane".to_string(),
            ));
            rows.extend(self.orphan_matrix_rows.iter().map(|index| {
                SplicingMatrixDisplayRow::Transcript {
                    lane_index: None,
                    matrix_row_index: Some(*index),
                }
            }));
        }
        rows
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    // Explicitly synthetic, in-memory annotation/evidence; no natural-gene fixture.
    fn view() -> SplicingExpertView {
        let lanes = [
            (30, "enst_synthetic_c.2", false),
            (10, "ENST_SYNTHETIC_A.7", false),
            (20, "ENST_SYNTHETIC_B.1", true),
        ]
        .into_iter()
        .map(|(id, transcript, target)| {
            json!({
                "transcript_feature_id":id,"transcript_id":transcript,"label":transcript,
                "strand":"+","exons":[],"introns":[],"has_target_feature":target
            })
        })
        .collect::<Vec<_>>();
        let matrices = [20, 30, 10].into_iter().map(|id| json!({
            "transcript_feature_id":id,"transcript_id":format!("row_{id}"),"label":"synthetic",
            "exon_presence":[]
        })).collect::<Vec<_>>();
        serde_json::from_value(json!({"seq_id":"synthetic","target_feature_id":20,
            "group_label":"synthetic","strand":"+","region_start_1based":1,"region_end_1based":100,
            "transcript_count":3,"unique_exon_count":0,"instruction":"synthetic",
            "transcripts":lanes,"matrix_rows":matrices,"unique_exons":[],"boundaries":[],
            "junctions":[],"events":[],"presentation_fingerprint_sha256":"retained-exactly"}))
        .unwrap()
    }

    #[test]
    fn splicing_boundary_summary_retains_oriented_exceptions_and_all_source_markers() {
        let mut view = view();
        for lane in &view.transcripts {
            view.boundaries.push(serde_json::from_value(json!({
                "transcript_feature_id":lane.transcript_feature_id,"transcript_id":lane.transcript_id,
                "side":"donor","position_1based":10,"partner_position_1based":50,
                "motif_2bp":"GT","paired_motif_signature":"GT-AG","canonical":true,"canonical_pair":true,
                "annotation":format!("unique source note {}",lane.transcript_feature_id)
            })).unwrap());
        }
        let grouped = boundary_summaries(&view);
        assert_eq!(grouped.len(), 1);
        assert_eq!(grouped[0].marker_indices, vec![0, 1, 2]);
        assert_eq!(grouped[0].transcript_feature_ids, vec![10, 20, 30]);
        view.transcripts[2].strand = "-".to_string();
        assert_eq!(boundary_summaries(&view).len(), 2);
        let mut exceptional = view.boundaries[0].clone();
        exceptional.canonical_pair = false;
        exceptional.paired_motif_signature = "GT-TT".to_string();
        view.boundaries.push(exceptional);
        let before = serde_json::to_value(&view).unwrap();
        let grouped = boundary_summaries(&view);
        assert_eq!(grouped.len(), 3);
        assert_eq!(
            grouped
                .iter()
                .map(|row| row.marker_indices.len())
                .sum::<usize>(),
            4
        );
        assert_eq!(serde_json::to_value(&view).unwrap(), before);
    }

    #[test]
    fn splicing_boundary_presentation_keeps_readable_sources_and_opposite_strands() {
        // Hand-crafted display-only records, including a deliberately long synthetic ID.
        let mut view = view();
        view.transcripts[0].transcript_id = format!("SYNTHETIC_{}", "LONG_ID_".repeat(30));
        view.transcripts[2].strand = "-".to_string();
        for lane in &view.transcripts {
            view.boundaries.push(serde_json::from_value(json!({
                "transcript_feature_id":lane.transcript_feature_id,"transcript_id":lane.transcript_id,
                "side":"donor","position_1based":10,"partner_position_1based":50,
                "motif_2bp":"GT","paired_motif_signature":"GT-AG","canonical":true,"canonical_pair":true,
                "annotation":format!("synthetic annotation {}",lane.transcript_feature_id)
            })).unwrap());
        }
        let before = serde_json::to_value(&view).unwrap();
        let rows = boundary_presentation_rows(&view);
        assert_eq!(rows.len(), 2);
        let plus = rows.iter().find(|row| row.summary.strand == "+").unwrap();
        assert_eq!(plus.label, "10 / 50 | + | donor:GT | GT-AG | 2 transcripts");
        assert_eq!(plus.summary.marker_indices, [0, 1]);
        assert_eq!(plus.summary.transcript_feature_ids, [10, 30]);
        assert_eq!(
            plus.hover_text,
            format!(
                "n-30 {}: synthetic annotation 30\nn-10 {}: synthetic annotation 10",
                view.transcripts[0].transcript_id, view.transcripts[1].transcript_id
            )
        );
        let minus = rows.iter().find(|row| row.summary.strand == "-").unwrap();
        assert_eq!(
            minus.label,
            "10 / 50 | - | donor:GT | GT-AG | 1 transcripts"
        );
        assert!(!plus.exceptional);
        assert!(!minus.exceptional);
        assert!(!plus.hover_text.contains("SplicingBoundaryMarker"));
        assert_eq!(serde_json::to_value(&view).unwrap(), before);
    }

    #[test]
    fn splicing_boundary_exception_includes_noncanonical_pairs() {
        // Synthetic boundary: its own dinucleotide and paired signature are independent facts.
        let mut view = view();
        let mut marker: SplicingBoundaryMarker = serde_json::from_value(json!({
            "transcript_feature_id":30,"transcript_id":"synthetic",
            "side":"acceptor","position_1based":50,"partner_position_1based":10,
            "motif_2bp":"AG","paired_motif_signature":"GC-AG","canonical":true,"canonical_pair":false
        })).unwrap();
        for (canonical, canonical_pair, exceptional) in [
            (true, true, false),
            (true, false, true),
            (false, true, true),
            (false, false, true),
        ] {
            marker.canonical = canonical;
            marker.canonical_pair = canonical_pair;
            assert_eq!(boundary_is_exceptional(&marker), exceptional);
            view.boundaries = vec![marker.clone()];
            let row = boundary_presentation_rows(&view).pop().unwrap();
            assert_eq!(row.exceptional, exceptional);
            assert_eq!(row.label.contains("AG*"), exceptional);
        }
    }

    #[test]
    fn splicing_layout_is_display_only_and_joins_shuffled_matrix_ids() {
        let mut view = view();
        view.matrix_rows[1].uniprot_reference.status = SplicingUniprotReferenceStatus::Referenced;
        view.matrix_rows[0].uniprot_reference.status =
            SplicingUniprotReferenceStatus::NotReferenced;
        let before = serde_json::to_value(&view).unwrap();
        let layout = view.uniprot_presentation_layout();
        assert_eq!(
            layout
                .lanes
                .iter()
                .map(|lane| lane.lane_index)
                .collect::<Vec<_>>(),
            [0, 2, 1]
        );
        assert_eq!(
            layout
                .lanes
                .iter()
                .map(|lane| lane.matrix_row_index)
                .collect::<Vec<_>>(),
            [Some(1), Some(0), Some(2)]
        );
        assert!(layout.show_headers);
        assert_eq!(serde_json::to_value(&view).unwrap(), before);
        assert_eq!(view.presentation_fingerprint_sha256, "retained-exactly");

        let mut tied = view.clone();
        for row in &mut tied.matrix_rows {
            row.uniprot_reference.status = SplicingUniprotReferenceStatus::Referenced;
        }
        tied.transcripts[0].transcript_id = "enst_synthetic_a.1".to_string();
        let before = serde_json::to_value(&tied).unwrap();
        assert_eq!(
            tied.uniprot_presentation_layout()
                .lanes
                .iter()
                .map(|lane| lane.transcript_feature_id)
                .collect::<Vec<_>>(),
            [20, 10, 30],
            "target first, then normalized-ID ties use feature ID, not version or payload order"
        );
        assert_eq!(serde_json::to_value(&tied).unwrap(), before);
    }

    #[test]
    fn splicing_layout_legacy_target_sort_and_header_hit_tests_are_explicit() {
        let mut view = view();
        let layout = view.uniprot_presentation_layout();
        assert!(!layout.show_headers);
        assert_eq!(
            layout
                .lanes
                .iter()
                .map(|lane| lane.transcript_feature_id)
                .collect::<Vec<_>>(),
            [20, 10, 30]
        );
        view.matrix_rows[0].uniprot_reference.status =
            SplicingUniprotReferenceStatus::NotReferenced;
        let layout = view.uniprot_presentation_layout();
        let canvas = layout.canvas(30.0, 20.0);
        assert_eq!(canvas.lane_at_y(0.0), None);
        assert_eq!(canvas.lane_at_y(19.999), None);
        assert_eq!(canvas.lane_at_y(20.0), Some(2));
        assert_eq!(canvas.centre_for_lane(2), Some(35.0));
        assert_eq!(canvas.lane_at_y(50.0), None);
        assert_eq!(canvas.lane_at_y(canvas.height), None);
        assert_eq!(canvas.lane_at_y(f32::NAN), None);
        assert!(layout.canvas(f32::NAN, 20.0).rows.is_empty());
    }

    #[test]
    fn splicing_layout_missing_duplicates_and_orphans_do_not_borrow_evidence() {
        let mut view = view();
        view.matrix_rows[1].transcript_feature_id = 20;
        view.matrix_rows[0].uniprot_reference.status = SplicingUniprotReferenceStatus::Referenced;
        let layout = view.uniprot_presentation_layout();
        assert!(
            layout
                .lanes
                .iter()
                .all(|lane| lane.status == SplicingUniprotReferenceStatus::NotEvaluated)
        );
        assert_eq!(layout.orphan_matrix_rows, [0, 1]);
        assert!(
            layout
                .diagnostics
                .iter()
                .any(|message| message.contains("Duplicate matrix"))
        );
        assert_eq!(layout.matrix_display_rows().len(), 6);
        view.transcripts.push(view.transcripts[1].clone());
        let layout = view.uniprot_presentation_layout();
        assert!(
            layout
                .diagnostics
                .iter()
                .any(|message| message.contains("Duplicate transcript"))
        );
        assert_eq!(layout.orphan_matrix_rows, [0, 1, 2]);
        assert_eq!(layout.lanes.len(), 4);
    }

    #[test]
    fn splicing_badge_counts_only_exact_sources_with_stable_accession_order() {
        use crate::{SplicingUniprotReferenceSource, UniprotEnsemblLinkedXref};
        let mut evidence = SplicingUniprotReferenceEvidence {
            status: SplicingUniprotReferenceStatus::Referenced,
            sources: vec![SplicingUniprotReferenceSource {
                accession: "GENE_ONLY".into(),
                matched_locus_gene_id: Some("ENSG_SYNTHETIC".into()),
                ..Default::default()
            }],
            diagnostics: vec![],
        };
        assert_eq!(reference_badge(&evidence), None);
        for accession in ["Z_SYNTHETIC", "A_SYNTHETIC"] {
            evidence.sources.push(SplicingUniprotReferenceSource {
                accession: accession.into(),
                reviewed: Some(false),
                matched_transcript_xrefs: vec![UniprotEnsemblLinkedXref {
                    transcript_id: Some("ENST_SYNTHETIC".into()),
                    ..Default::default()
                }],
                ..Default::default()
            });
        }
        evidence.diagnostics.push("Ambiguous relation".into());
        assert_eq!(
            reference_badge(&evidence).as_deref(),
            Some("UniProt A_SYNTHETIC [unreviewed] +1 ?")
        );
        evidence.sources.reverse();
        assert_eq!(
            reference_badge(&evidence).as_deref(),
            Some("UniProt A_SYNTHETIC [unreviewed] +1 ?")
        );
        evidence.diagnostics = vec!["Multiple explicit locus gene IDs".into()];
        assert_eq!(
            reference_badge(&evidence).as_deref(),
            Some("UniProt A_SYNTHETIC [unreviewed] +1 ?")
        );
        assert!(REFERENCE_BADGE_LEGEND.contains("? = loaded-evidence diagnostic"));
        assert!(REFERENCE_BADGE_LEGEND.contains("not necessarily transcript ambiguity"));
    }
}
