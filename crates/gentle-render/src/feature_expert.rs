//! Feature expert-view SVG renderer.

use gentle_protocol::{
    FeatureExpertView, IsoformArchitectureExpertView, RestrictionSiteExpertView,
    SplicingExonSummary, SplicingExpertView, TfbsExpertView,
};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use svg::Document;
use svg::node::element::path::Data;
use svg::node::element::{Circle, Line, Path, Rectangle, Text};

const W: f32 = 1200.0;
const H: f32 = 700.0;
const PROTEIN_DOMAIN_LABEL_MAX_CHARS: usize = 44;
const PROTEIN_DOMAIN_LABEL_FONT_SIZE: f32 = 8.0;
const PROTEIN_DOMAIN_LABEL_CHAR_W: f32 = 4.8;
const PROTEIN_DOMAIN_LABEL_BOX_X_PAD: f32 = 3.0;
const PROTEIN_DOMAIN_LABEL_BOX_Y_PAD: f32 = 2.0;
const PROTEIN_DOMAIN_LABEL_ROW_PITCH: f32 = 12.0;
const PROTEIN_DOMAIN_LABEL_ROW_GAP: f32 = 8.0;
const PROTEIN_DOMAIN_LABEL_DOMAIN_GAP: f32 = 7.0;
const PROTEIN_DOMAIN_HALF_HEIGHT: f32 = 7.0;
const PROTEIN_TOPOLOGY_ROW_HEIGHT: f32 = 10.0;
const PROTEIN_TOPOLOGY_ROW_GAP: f32 = 4.0;
const PROTEIN_TOPOLOGY_TRACK_GAP: f32 = 8.0;
const PROTEIN_TOPOLOGY_CHAR_W: f32 = 4.4;
const PROTEIN_TOPOLOGY_LABEL_FONT_SIZE: f32 = 7.0;
const PROTEIN_LANE_BOTTOM_PAD: f32 = 10.0;
const TRANSCRIPT_PRODUCT_LANE_GAP: f32 = 16.0;
const TRANSCRIPT_PRODUCT_TOP_PAD: f32 = 10.0;
const TRANSCRIPT_PRODUCT_CONNECTOR_GAP: f32 = 22.0;
const COMPRESSED_EXON_HALF_HEIGHT: f32 = 7.0;
const COMPRESSED_EXON_INTRON_GAP: f32 = 18.0;
const COMPRESSED_EXON_MIN_WIDTH: f32 = 28.0;
const COMPRESSED_EXON_MAX_WIDTH: f32 = 92.0;
const PROTEIN_CONTRIBUTION_HALF_HEIGHT: f32 = 4.5;

#[derive(Debug, Clone)]
struct ProteinDomainLabelPlacement {
    compact_label: String,
    box_left: f32,
    box_right: f32,
    row_index: usize,
}

#[derive(Debug, Clone)]
struct ProteinLabelLayout {
    block_height: f32,
    placements: Vec<ProteinDomainLabelPlacement>,
}

#[derive(Debug, Clone)]
struct ProteinIntervalRowLayout {
    row_count: usize,
    row_indices: Vec<usize>,
}

#[derive(Debug, Clone, Copy)]
struct ExonSegmentColor {
    fill: &'static str,
    stroke: &'static str,
}

#[derive(Debug, Clone)]
struct CompressedTranscriptExonBox {
    exon_key: (usize, usize),
    family_range: (usize, usize),
    x1: f32,
    x2: f32,
    fill: &'static str,
    stroke: &'static str,
}

#[derive(Debug, Clone)]
struct CompressedTranscriptBlock {
    segment_index: usize,
    exon_key: (usize, usize),
    family_range: (usize, usize),
    x1: f32,
    x2: f32,
    fill: &'static str,
    stroke: &'static str,
}

#[derive(Debug, Clone)]
struct CompressedTranscriptLaneLayout {
    exon_boxes: Vec<CompressedTranscriptExonBox>,
    coding_blocks: Vec<CompressedTranscriptBlock>,
}

#[derive(Debug, Clone)]
struct LocalProteinContributionSegment {
    segment_index: usize,
    exon_key: (usize, usize),
    family_range: (usize, usize),
    reference_aa_start: usize,
    reference_aa_end: usize,
    local_aa_start: usize,
    local_aa_end: usize,
    fill: &'static str,
    stroke: &'static str,
}

#[derive(Debug, Clone)]
struct LocalProteinRect {
    segment_index: usize,
    exon_key: (usize, usize),
    family_range: (usize, usize),
    x1: f32,
    x2: f32,
    fill: &'static str,
    stroke: &'static str,
}

#[derive(Debug, Clone)]
struct ExonFamilyRegistry {
    family_ranges: Vec<(usize, usize)>,
    family_index_by_key: BTreeMap<(usize, usize), usize>,
}

#[derive(Debug, Clone)]
struct ExonFamilyColumn {
    family_range: (usize, usize),
    x1: f32,
    x2: f32,
    fill: &'static str,
    stroke: &'static str,
}

#[derive(Debug, Clone)]
struct LocalProteinFeaturePiece {
    domain_index: usize,
    x1: f32,
    x2: f32,
}

#[derive(Debug, Clone)]
struct LocalProteinLaneLayout {
    local_length_aa: usize,
    contribution_rects: Vec<LocalProteinRect>,
    overlay_pieces: Vec<LocalProteinFeaturePiece>,
    topology_pieces: Vec<LocalProteinFeaturePiece>,
    overlay_labels: ProteinLabelLayout,
    overlay_label_domain_indices: Vec<usize>,
    topology_rows: ProteinIntervalRowLayout,
    lane_height: f32,
    reference_span_label: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ProteinFeatureTrack {
    Overlay,
    Topology,
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

fn nuc_color(base_idx: usize) -> &'static str {
    match base_idx {
        0 => "#2b8a3e", // A
        1 => "#1d4ed8", // C
        2 => "#d97706", // G
        3 => "#b91c1c", // T
        _ => "#555555",
    }
}

fn wrap_text(input: &str, max_chars: usize) -> Vec<String> {
    let mut out = Vec::new();
    let mut current = String::new();
    for token in input.split_whitespace() {
        let next_len = if current.is_empty() {
            token.len()
        } else {
            current.len() + 1 + token.len()
        };
        if !current.is_empty() && next_len > max_chars {
            out.push(current);
            current = token.to_string();
        } else {
            if !current.is_empty() {
                current.push(' ');
            }
            current.push_str(token);
        }
    }
    if !current.is_empty() {
        out.push(current);
    }
    if out.is_empty() {
        out.push(String::new());
    }
    out
}

fn compact_protein_domain_label(input: &str, max_chars: usize) -> String {
    if input.chars().count() <= max_chars {
        return input.to_string();
    }
    if max_chars <= 3 {
        return input.chars().take(max_chars).collect();
    }
    let keep = max_chars.saturating_sub(3);
    let mut out = String::new();
    for ch in input.chars().take(keep) {
        out.push(ch);
    }
    out.push_str("...");
    out
}

fn estimate_protein_domain_label_box_width(label: &str) -> f32 {
    label.chars().count() as f32 * PROTEIN_DOMAIN_LABEL_CHAR_W
        + PROTEIN_DOMAIN_LABEL_BOX_X_PAD * 2.0
}

fn protein_feature_key(label: &str) -> Option<&str> {
    let trimmed = label.trim();
    let prefix = trimmed.split(':').next().unwrap_or(trimmed).trim();
    (!prefix.is_empty()).then_some(prefix)
}

fn protein_feature_track(label: &str) -> ProteinFeatureTrack {
    match protein_feature_key(label)
        .unwrap_or("")
        .trim()
        .to_ascii_uppercase()
        .as_str()
    {
        "TOPO_DOM" | "TRANSMEM" | "INTRAMEM" | "SIGNAL" | "TRANSIT" => {
            ProteinFeatureTrack::Topology
        }
        _ => ProteinFeatureTrack::Overlay,
    }
}

fn topology_feature_fill(feature_label: &str) -> &'static str {
    match protein_feature_key(feature_label)
        .unwrap_or("")
        .trim()
        .to_ascii_uppercase()
        .as_str()
    {
        "SIGNAL" | "TRANSIT" => "#86efac",
        "TRANSMEM" | "INTRAMEM" => "#fca5a5",
        "TOPO_DOM" => "#bfdbfe",
        _ => "#cbd5e1",
    }
}

fn topology_feature_stroke(feature_label: &str) -> &'static str {
    match protein_feature_key(feature_label)
        .unwrap_or("")
        .trim()
        .to_ascii_uppercase()
        .as_str()
    {
        "SIGNAL" | "TRANSIT" => "#15803d",
        "TRANSMEM" | "INTRAMEM" => "#b91c1c",
        "TOPO_DOM" => "#1d4ed8",
        _ => "#475569",
    }
}

fn topology_feature_label(feature_label: &str) -> String {
    let trimmed = feature_label.trim();
    let note = trimmed
        .split_once(':')
        .map(|(_, note)| note.trim())
        .filter(|note| !note.is_empty());
    match protein_feature_key(trimmed)
        .unwrap_or("")
        .trim()
        .to_ascii_uppercase()
        .as_str()
    {
        "SIGNAL" => note.unwrap_or("signal peptide").to_string(),
        "TRANSIT" => note.unwrap_or("transit peptide").to_string(),
        "TRANSMEM" => note.unwrap_or("transmembrane").to_string(),
        "INTRAMEM" => note.unwrap_or("intramembrane").to_string(),
        "TOPO_DOM" => note.unwrap_or("topological domain").to_string(),
        _ => trimmed.to_string(),
    }
}

fn compact_topology_label(feature_label: &str, available_width: f32) -> Option<String> {
    let max_chars = ((available_width - 6.0) / PROTEIN_TOPOLOGY_CHAR_W).floor() as usize;
    if max_chars < 3 {
        return None;
    }
    Some(compact_protein_domain_label(
        &topology_feature_label(feature_label),
        max_chars,
    ))
}

fn layout_interval_rows(intervals: &[(f32, f32)]) -> ProteinIntervalRowLayout {
    #[derive(Debug, Clone)]
    struct PendingInterval {
        original_index: usize,
        left: f32,
        right: f32,
    }

    let mut pending = intervals
        .iter()
        .enumerate()
        .map(|(original_index, (left, right))| PendingInterval {
            original_index,
            left: *left,
            right: *right,
        })
        .collect::<Vec<_>>();
    pending.sort_by(|a, b| {
        a.left
            .partial_cmp(&b.left)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                a.right
                    .partial_cmp(&b.right)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    let mut row_right_edges: Vec<f32> = Vec::new();
    let mut row_indices = vec![0usize; intervals.len()];
    for interval in pending {
        let mut row_index = None;
        for (idx, row_right) in row_right_edges.iter_mut().enumerate() {
            if interval.left >= *row_right + PROTEIN_TOPOLOGY_ROW_GAP {
                *row_right = interval.right;
                row_index = Some(idx);
                break;
            }
        }
        let row_index = row_index.unwrap_or_else(|| {
            row_right_edges.push(interval.right);
            row_right_edges.len() - 1
        });
        row_indices[interval.original_index] = row_index;
    }
    ProteinIntervalRowLayout {
        row_count: row_right_edges.len(),
        row_indices,
    }
}

fn layout_protein_domain_labels(
    domains: &[(f32, f32, String)],
    left: f32,
    right: f32,
) -> ProteinLabelLayout {
    #[derive(Debug, Clone)]
    struct PendingLabel {
        domain_index: usize,
        compact_label: String,
        box_left: f32,
        box_right: f32,
    }

    let mut pending = domains
        .iter()
        .enumerate()
        .map(|(domain_index, (domain_left, domain_right, raw_label))| {
            let compact_label =
                compact_protein_domain_label(raw_label, PROTEIN_DOMAIN_LABEL_MAX_CHARS);
            let mut box_width = estimate_protein_domain_label_box_width(&compact_label);
            let available_width = (right - left).max(1.0);
            if box_width > available_width {
                box_width = available_width;
            }
            let domain_center = (*domain_left + *domain_right) * 0.5;
            let max_left = (right - box_width).max(left);
            let box_left = (domain_center - box_width * 0.5).clamp(left, max_left);
            PendingLabel {
                domain_index,
                compact_label,
                box_left,
                box_right: box_left + box_width,
            }
        })
        .collect::<Vec<_>>();
    pending.sort_by(|a, b| {
        a.box_left
            .partial_cmp(&b.box_left)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                b.box_right
                    .partial_cmp(&a.box_right)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    let mut row_right_edges: Vec<f32> = Vec::new();
    let mut placements = vec![None; pending.len()];
    for label in pending {
        let mut row_index = None;
        for (idx, row_right) in row_right_edges.iter_mut().enumerate() {
            if label.box_left >= *row_right + PROTEIN_DOMAIN_LABEL_ROW_GAP {
                *row_right = label.box_right;
                row_index = Some(idx);
                break;
            }
        }
        let row_index = row_index.unwrap_or_else(|| {
            row_right_edges.push(label.box_right);
            row_right_edges.len() - 1
        });
        placements[label.domain_index] = Some(ProteinDomainLabelPlacement {
            compact_label: label.compact_label,
            box_left: label.box_left,
            box_right: label.box_right,
            row_index,
        });
    }

    let label_rows = row_right_edges.len();
    let label_band_height = if label_rows == 0 {
        0.0
    } else {
        PROTEIN_DOMAIN_LABEL_FONT_SIZE
            + (label_rows.saturating_sub(1) as f32) * PROTEIN_DOMAIN_LABEL_ROW_PITCH
            + PROTEIN_DOMAIN_LABEL_BOX_Y_PAD
    };
    ProteinLabelLayout {
        block_height: label_band_height,
        placements: placements.into_iter().flatten().collect(),
    }
}

fn palette_exon_color(segment_index: usize) -> ExonSegmentColor {
    const PALETTE: [ExonSegmentColor; 10] = [
        ExonSegmentColor {
            fill: "#2563eb",
            stroke: "#1d4ed8",
        },
        ExonSegmentColor {
            fill: "#16a34a",
            stroke: "#15803d",
        },
        ExonSegmentColor {
            fill: "#ea580c",
            stroke: "#c2410c",
        },
        ExonSegmentColor {
            fill: "#7c3aed",
            stroke: "#6d28d9",
        },
        ExonSegmentColor {
            fill: "#db2777",
            stroke: "#be185d",
        },
        ExonSegmentColor {
            fill: "#0891b2",
            stroke: "#0e7490",
        },
        ExonSegmentColor {
            fill: "#65a30d",
            stroke: "#4d7c0f",
        },
        ExonSegmentColor {
            fill: "#dc2626",
            stroke: "#b91c1c",
        },
        ExonSegmentColor {
            fill: "#d97706",
            stroke: "#b45309",
        },
        ExonSegmentColor {
            fill: "#0f766e",
            stroke: "#115e59",
        },
    ];
    PALETTE[segment_index % PALETTE.len()]
}

fn genomic_interval_key(start_1based: usize, end_1based: usize) -> (usize, usize) {
    (start_1based.min(end_1based), start_1based.max(end_1based))
}

fn exon_identity_key(range: &gentle_protocol::SplicingRange) -> (usize, usize) {
    genomic_interval_key(range.start_1based, range.end_1based)
}

fn exon_identity_attr(exon_key: (usize, usize)) -> String {
    format!("{}..{}", exon_key.0, exon_key.1)
}

fn exon_family_attr(family_range: (usize, usize)) -> String {
    format!("{}..{}", family_range.0, family_range.1)
}

fn best_matching_exon_identity_key(
    ranges: &[gentle_protocol::SplicingRange],
    genomic_start_1based: usize,
    genomic_end_1based: usize,
) -> Option<(usize, usize)> {
    let target = genomic_interval_key(genomic_start_1based, genomic_end_1based);
    let mut containing = ranges
        .iter()
        .map(exon_identity_key)
        .filter(|(start, end)| *start <= target.0 && *end >= target.1)
        .collect::<Vec<_>>();
    containing.sort_by(|left, right| {
        let left_len = left.1.saturating_sub(left.0).saturating_add(1);
        let right_len = right.1.saturating_sub(right.0).saturating_add(1);
        left_len.cmp(&right_len).then(left.cmp(right))
    });
    if let Some(best) = containing.into_iter().next() {
        return Some(best);
    }

    let mut overlaps = ranges
        .iter()
        .map(exon_identity_key)
        .filter_map(|key| {
            let overlap_start = key.0.max(target.0);
            let overlap_end = key.1.min(target.1);
            (overlap_end >= overlap_start).then_some((key, overlap_end - overlap_start + 1))
        })
        .collect::<Vec<_>>();
    overlaps.sort_by(|left, right| {
        right
            .1
            .cmp(&left.1)
            .then_with(|| {
                let left_len = left.0.1.saturating_sub(left.0.0).saturating_add(1);
                let right_len = right.0.1.saturating_sub(right.0.0).saturating_add(1);
                left_len.cmp(&right_len)
            })
            .then(left.0.cmp(&right.0))
    });
    overlaps.into_iter().map(|(key, _)| key).next()
}

fn lane_exon_identity_key_for_coding_range(
    lane: &gentle_protocol::IsoformArchitectureTranscriptLane,
    genomic_start_1based: usize,
    genomic_end_1based: usize,
) -> (usize, usize) {
    best_matching_exon_identity_key(
        &lane.transcript_exons,
        genomic_start_1based,
        genomic_end_1based,
    )
    .or_else(|| {
        best_matching_exon_identity_key(&lane.exons, genomic_start_1based, genomic_end_1based)
    })
    .unwrap_or_else(|| genomic_interval_key(genomic_start_1based, genomic_end_1based))
}

fn build_exon_family_registry(view: &IsoformArchitectureExpertView) -> ExonFamilyRegistry {
    let mut keys = BTreeSet::new();
    for lane in &view.transcript_lanes {
        let exon_ranges = if !lane.transcript_exons.is_empty() {
            &lane.transcript_exons
        } else {
            &lane.exons
        };
        for exon in exon_ranges {
            keys.insert(exon_identity_key(exon));
        }
        for segment in &lane.cds_to_protein_segments {
            keys.insert(lane_exon_identity_key_for_coding_range(
                lane,
                segment.genomic_start_1based,
                segment.genomic_end_1based,
            ));
        }
    }
    let mut sorted_keys = keys.into_iter().collect::<Vec<_>>();
    sorted_keys.sort();

    let mut family_ranges: Vec<(usize, usize)> = Vec::new();
    let mut family_index_by_key: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for key in sorted_keys {
        let family_index = if let Some((_, family_end)) = family_ranges.last_mut() {
            if key.0 <= *family_end {
                *family_end = (*family_end).max(key.1);
                family_ranges.len() - 1
            } else {
                family_ranges.push(key);
                family_ranges.len() - 1
            }
        } else {
            family_ranges.push(key);
            0
        };
        family_index_by_key.insert(key, family_index);
    }

    ExonFamilyRegistry {
        family_ranges,
        family_index_by_key,
    }
}

fn exon_family_index_for_key(
    exon_key: (usize, usize),
    exon_family_registry: &ExonFamilyRegistry,
) -> usize {
    if let Some(index) = exon_family_registry
        .family_index_by_key
        .get(&exon_key)
        .copied()
    {
        return index;
    }
    exon_family_registry
        .family_ranges
        .iter()
        .enumerate()
        .find_map(|(index, family_range)| {
            let overlap_start = family_range.0.max(exon_key.0);
            let overlap_end = family_range.1.min(exon_key.1);
            (overlap_end >= overlap_start).then_some(index)
        })
        .unwrap_or(0)
}

fn exon_family_range_for_key(
    exon_key: (usize, usize),
    exon_family_registry: &ExonFamilyRegistry,
) -> (usize, usize) {
    exon_family_registry.family_ranges[exon_family_index_for_key(exon_key, exon_family_registry)]
}

fn exon_color_for_key(
    exon_key: (usize, usize),
    exon_family_registry: &ExonFamilyRegistry,
) -> ExonSegmentColor {
    palette_exon_color(exon_family_index_for_key(exon_key, exon_family_registry))
}

fn compressed_exon_width(exon_bp: usize) -> f32 {
    let bp = exon_bp.max(1) as f32;
    (18.0 + bp.sqrt() * 3.1).clamp(COMPRESSED_EXON_MIN_WIDTH, COMPRESSED_EXON_MAX_WIDTH)
}

fn layout_exon_family_columns(
    exon_family_registry: &ExonFamilyRegistry,
    left: f32,
    right: f32,
    reverse_display: bool,
) -> BTreeMap<usize, ExonFamilyColumn> {
    let mut ordered_families = exon_family_registry
        .family_ranges
        .iter()
        .copied()
        .enumerate()
        .collect::<Vec<_>>();
    if reverse_display {
        ordered_families.reverse();
    }

    let raw_widths = ordered_families
        .iter()
        .map(|(_, family_range)| {
            compressed_exon_width(
                family_range
                    .1
                    .saturating_sub(family_range.0)
                    .saturating_add(1),
            )
        })
        .collect::<Vec<_>>();
    let available_width = (right - left).max(1.0);
    let raw_total = raw_widths.iter().sum::<f32>()
        + COMPRESSED_EXON_INTRON_GAP * raw_widths.len().saturating_sub(1) as f32;
    let scale = if raw_total > available_width {
        (available_width / raw_total).clamp(0.18, 1.0)
    } else {
        1.0
    };

    let mut columns = BTreeMap::new();
    let mut cursor_x = left;
    for ((family_index, family_range), raw_width) in ordered_families.into_iter().zip(raw_widths) {
        let width = (raw_width * scale).max(8.0);
        let color = palette_exon_color(family_index);
        columns.insert(
            family_index,
            ExonFamilyColumn {
                family_range,
                x1: cursor_x,
                x2: cursor_x + width,
                fill: color.fill,
                stroke: color.stroke,
            },
        );
        cursor_x += width + COMPRESSED_EXON_INTRON_GAP * scale;
    }
    columns
}

fn build_local_protein_contribution_segments(
    lane: &gentle_protocol::IsoformArchitectureTranscriptLane,
    exon_family_registry: &ExonFamilyRegistry,
) -> Vec<LocalProteinContributionSegment> {
    let mut local_start_aa = 1usize;
    let mut out = Vec::new();
    for (segment_index, segment) in lane.cds_to_protein_segments.iter().enumerate() {
        let reference_start_aa = segment.aa_start.min(segment.aa_end);
        let reference_end_aa = segment.aa_start.max(segment.aa_end);
        if reference_end_aa < reference_start_aa {
            continue;
        }
        let segment_len_aa = reference_end_aa
            .saturating_sub(reference_start_aa)
            .saturating_add(1);
        if segment_len_aa == 0 {
            continue;
        }
        let local_end_aa = local_start_aa
            .saturating_add(segment_len_aa)
            .saturating_sub(1);
        let exon_key = lane_exon_identity_key_for_coding_range(
            lane,
            segment.genomic_start_1based,
            segment.genomic_end_1based,
        );
        let family_range = exon_family_range_for_key(exon_key, exon_family_registry);
        let color = exon_color_for_key(exon_key, exon_family_registry);
        out.push(LocalProteinContributionSegment {
            segment_index,
            exon_key,
            family_range,
            reference_aa_start: reference_start_aa,
            reference_aa_end: reference_end_aa,
            local_aa_start: local_start_aa,
            local_aa_end: local_end_aa,
            fill: color.fill,
            stroke: color.stroke,
        });
        local_start_aa = local_end_aa.saturating_add(1);
    }
    out
}

fn map_reference_interval_to_local_segments(
    interval_start_aa: usize,
    interval_end_aa: usize,
    contributions: &[LocalProteinContributionSegment],
    lane_reference_start_aa: Option<usize>,
    lane_reference_end_aa: Option<usize>,
) -> Vec<(usize, usize)> {
    let start_aa = interval_start_aa.min(interval_end_aa);
    let end_aa = interval_start_aa.max(interval_end_aa);
    let mut out = Vec::new();
    if !contributions.is_empty() {
        for segment in contributions {
            if segment.reference_aa_end < start_aa || segment.reference_aa_start > end_aa {
                continue;
            }
            let overlap_start = start_aa.max(segment.reference_aa_start);
            let overlap_end = end_aa.min(segment.reference_aa_end);
            if overlap_end < overlap_start {
                continue;
            }
            let local_start = segment
                .local_aa_start
                .saturating_add(overlap_start.saturating_sub(segment.reference_aa_start));
            let local_end = segment
                .local_aa_start
                .saturating_add(overlap_end.saturating_sub(segment.reference_aa_start));
            if local_end >= local_start {
                out.push((local_start, local_end));
            }
        }
        return out;
    }

    if let (Some(reference_start), Some(reference_end)) =
        (lane_reference_start_aa, lane_reference_end_aa)
    {
        let overlap_start = start_aa.max(reference_start.min(reference_end));
        let overlap_end = end_aa.min(reference_start.max(reference_end));
        if overlap_end >= overlap_start {
            let local_start = overlap_start.saturating_sub(reference_start.min(reference_end)) + 1;
            let local_end = overlap_end.saturating_sub(reference_start.min(reference_end)) + 1;
            out.push((local_start, local_end));
        }
        return out;
    }

    if end_aa >= start_aa {
        out.push((start_aa, end_aa));
    }
    out
}

fn layout_compressed_transcript_lane(
    lane: &gentle_protocol::IsoformArchitectureTranscriptLane,
    exon_family_columns: &BTreeMap<usize, ExonFamilyColumn>,
    exon_family_registry: &ExonFamilyRegistry,
) -> CompressedTranscriptLaneLayout {
    let transcript_exons = if !lane.transcript_exons.is_empty() {
        lane.transcript_exons.clone()
    } else {
        lane.exons.clone()
    };
    let mut lane_key_by_family = BTreeMap::new();
    for exon in &transcript_exons {
        let exon_key = exon_identity_key(exon);
        let family_index = exon_family_index_for_key(exon_key, exon_family_registry);
        lane_key_by_family.entry(family_index).or_insert(exon_key);
    }

    let mut exon_boxes = Vec::with_capacity(lane_key_by_family.len());
    for (family_index, exon_key) in lane_key_by_family {
        let Some(column) = exon_family_columns.get(&family_index) else {
            continue;
        };
        exon_boxes.push(CompressedTranscriptExonBox {
            exon_key,
            family_range: column.family_range,
            x1: column.x1,
            x2: column.x2,
            fill: column.fill,
            stroke: column.stroke,
        });
    }
    exon_boxes.sort_by(|left_box, right_box| {
        left_box
            .x1
            .partial_cmp(&right_box.x1)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let ordered_segments = lane
        .cds_to_protein_segments
        .iter()
        .enumerate()
        .map(|(segment_index, segment)| {
            (
                segment_index,
                segment.genomic_start_1based.min(segment.genomic_end_1based),
                segment.genomic_start_1based.max(segment.genomic_end_1based),
            )
        })
        .collect::<Vec<_>>();
    let mut coding_blocks = Vec::new();
    for (segment_index, overlap_start, overlap_end) in ordered_segments {
        let exon_key = lane_exon_identity_key_for_coding_range(lane, overlap_start, overlap_end);
        let family_index = exon_family_index_for_key(exon_key, exon_family_registry);
        let Some(column) = exon_family_columns.get(&family_index) else {
            continue;
        };
        let family_start = column.family_range.0;
        let family_end = column.family_range.1;
        let family_len = family_end
            .saturating_sub(family_start)
            .saturating_add(1)
            .max(1) as f32;
        let clipped_start = overlap_start.max(family_start).min(family_end);
        let clipped_end = overlap_end.min(family_end).max(family_start);
        let start_offset = if lane.strand.trim() == "-" {
            family_end.saturating_sub(clipped_end) as f32
        } else {
            clipped_start.saturating_sub(family_start) as f32
        };
        let end_offset_exclusive = if lane.strand.trim() == "-" {
            family_end.saturating_sub(clipped_start).saturating_add(1) as f32
        } else {
            clipped_end.saturating_sub(family_start).saturating_add(1) as f32
        };
        let x1 = column.x1 + (start_offset / family_len) * (column.x2 - column.x1).max(1.0);
        let x2 = column.x1 + (end_offset_exclusive / family_len) * (column.x2 - column.x1).max(1.0);
        coding_blocks.push(CompressedTranscriptBlock {
            segment_index,
            exon_key,
            family_range: column.family_range,
            x1: x1.min(x2),
            x2: x1.max(x2),
            fill: column.fill,
            stroke: column.stroke,
        });
    }

    CompressedTranscriptLaneLayout {
        exon_boxes,
        coding_blocks,
    }
}

fn layout_local_protein_lane(
    transcript_lane: &gentle_protocol::IsoformArchitectureTranscriptLane,
    protein_lane: &gentle_protocol::IsoformArchitectureProteinLane,
    left: f32,
    right: f32,
    exon_family_registry: &ExonFamilyRegistry,
) -> LocalProteinLaneLayout {
    let contributions =
        build_local_protein_contribution_segments(transcript_lane, exon_family_registry);
    let local_length_aa = contributions
        .last()
        .map(|segment| segment.local_aa_end)
        .or_else(|| protein_lane.expected_length_aa)
        .or_else(|| {
            match (
                protein_lane.reference_start_aa,
                protein_lane.reference_end_aa,
            ) {
                (Some(start), Some(end)) => Some(end.max(start).saturating_sub(start.min(end)) + 1),
                _ => None,
            }
        })
        .unwrap_or_else(|| {
            protein_lane
                .domains
                .iter()
                .map(|domain| domain.end_aa.max(domain.start_aa))
                .max()
                .unwrap_or(1)
        })
        .max(1);
    let local_span = local_length_aa as f32;
    let x_for_local_edge = |edge_1based_exclusive: usize| -> f32 {
        let clamped = edge_1based_exclusive.clamp(1, local_length_aa.saturating_add(1));
        left + ((clamped.saturating_sub(1)) as f32 / local_span) * (right - left).max(1.0)
    };
    let contribution_rects = contributions
        .iter()
        .map(|segment| LocalProteinRect {
            segment_index: segment.segment_index,
            exon_key: segment.exon_key,
            family_range: segment.family_range,
            x1: x_for_local_edge(segment.local_aa_start),
            x2: x_for_local_edge(segment.local_aa_end.saturating_add(1)),
            fill: segment.fill,
            stroke: segment.stroke,
        })
        .collect::<Vec<_>>();

    let mut overlay_unions = Vec::new();
    let mut overlay_pieces = Vec::new();
    let mut topology_pieces = Vec::new();
    for (domain_index, domain) in protein_lane.domains.iter().enumerate() {
        let local_ranges = map_reference_interval_to_local_segments(
            domain.start_aa,
            domain.end_aa.max(domain.start_aa),
            &contributions,
            protein_lane.reference_start_aa,
            protein_lane.reference_end_aa,
        );
        if local_ranges.is_empty() {
            continue;
        }
        let mut union_left = f32::MAX;
        let mut union_right = f32::MIN;
        for (local_start, local_end) in local_ranges {
            let x1 = x_for_local_edge(local_start);
            let x2 = x_for_local_edge(local_end.saturating_add(1));
            union_left = union_left.min(x1.min(x2));
            union_right = union_right.max(x1.max(x2));
            let piece = LocalProteinFeaturePiece {
                domain_index,
                x1: x1.min(x2),
                x2: x1.max(x2),
            };
            if protein_feature_track(&domain.name) == ProteinFeatureTrack::Topology {
                topology_pieces.push(piece);
            } else {
                overlay_pieces.push(piece);
            }
        }
        if union_right >= union_left {
            overlay_unions.push((union_left, union_right, domain.name.clone(), domain_index));
        }
    }

    let overlay_label_domain_indices = overlay_unions
        .iter()
        .map(|(_, _, _, domain_index)| *domain_index)
        .collect::<Vec<_>>();
    let overlay_labels = layout_protein_domain_labels(
        &overlay_unions
            .iter()
            .map(|(x1, x2, label, _)| (*x1, *x2, label.clone()))
            .collect::<Vec<_>>(),
        left,
        right,
    );
    let topology_rows = layout_interval_rows(
        &topology_pieces
            .iter()
            .map(|piece| (piece.x1, piece.x2))
            .collect::<Vec<_>>(),
    );
    let topology_block_height = if topology_rows.row_count == 0 {
        0.0
    } else {
        topology_rows.row_count as f32 * PROTEIN_TOPOLOGY_ROW_HEIGHT
            + (topology_rows.row_count.saturating_sub(1) as f32) * PROTEIN_TOPOLOGY_ROW_GAP
    };
    let lane_height = PROTEIN_DOMAIN_HALF_HEIGHT * 2.0
        + PROTEIN_CONTRIBUTION_HALF_HEIGHT * 2.0
        + if topology_block_height > 0.0 {
            PROTEIN_TOPOLOGY_TRACK_GAP + topology_block_height
        } else {
            0.0
        }
        + if overlay_labels.block_height > 0.0 {
            PROTEIN_DOMAIN_LABEL_DOMAIN_GAP + overlay_labels.block_height
        } else {
            0.0
        }
        + PROTEIN_LANE_BOTTOM_PAD;
    let reference_span_label = match (
        protein_lane.reference_start_aa,
        protein_lane.reference_end_aa,
    ) {
        (Some(start), Some(end)) => Some(format!("ref {}..{}", start.min(end), start.max(end))),
        _ => None,
    };

    LocalProteinLaneLayout {
        local_length_aa,
        contribution_rects,
        overlay_pieces,
        topology_pieces,
        overlay_labels,
        overlay_label_domain_indices,
        topology_rows,
        lane_height,
        reference_span_label,
    }
}

fn match_base_to_idx(base: Option<char>) -> Option<usize> {
    match base.map(|c| c.to_ascii_uppercase()) {
        Some('A') => Some(0),
        Some('C') => Some(1),
        Some('G') => Some(2),
        Some('T') => Some(3),
        _ => None,
    }
}

fn mix_rgb(from: [u8; 3], to: [u8; 3], t: f32) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    let lerp = |a: u8, b: u8| -> u8 {
        ((a as f32) + ((b as f32) - (a as f32)) * t)
            .round()
            .clamp(0.0, 255.0) as u8
    };
    [
        lerp(from[0], to[0]),
        lerp(from[1], to[1]),
        lerp(from[2], to[2]),
    ]
}

fn rgb_hex(rgb: [u8; 3]) -> String {
    format!("#{:02x}{:02x}{:02x}", rgb[0], rgb[1], rgb[2])
}

fn support_ratio_percent(support_count: usize, transcript_total: usize) -> f64 {
    if transcript_total == 0 {
        return 0.0;
    }
    support_count as f64 * 100.0 / transcript_total as f64
}

fn format_support_fraction(support_count: usize, transcript_total: usize) -> String {
    format!(
        "{}/{} ({:.1}%)",
        support_count,
        transcript_total,
        support_ratio_percent(support_count, transcript_total)
    )
}

fn splicing_matrix_cell_fill_hex(present: bool, support_ratio: f32) -> String {
    let ratio = support_ratio.clamp(0.0, 1.0);
    if present {
        rgb_hex(mix_rgb([191, 219, 254], [30, 64, 175], ratio))
    } else {
        rgb_hex(mix_rgb([248, 250, 252], [219, 234, 254], ratio))
    }
}

fn splicing_matrix_cell_text_hex(present: bool, support_ratio: f32) -> &'static str {
    if present && support_ratio >= 0.42 {
        "#ffffff"
    } else {
        "#0f172a"
    }
}

fn splicing_transition_cell_fill_hex(support_ratio: f32) -> String {
    let ratio = support_ratio.clamp(0.0, 1.0);
    if ratio > 0.0 {
        rgb_hex(mix_rgb([220, 252, 231], [22, 101, 52], ratio))
    } else {
        "#f8fafc".to_string()
    }
}

fn splicing_transition_cell_text_hex(support_ratio: f32) -> &'static str {
    if support_ratio >= 0.42 {
        "#ffffff"
    } else if support_ratio > 0.0 {
        "#14532d"
    } else {
        "#64748b"
    }
}

fn splicing_exon_length_bp(exon: &SplicingExonSummary) -> usize {
    exon.end_1based
        .max(exon.start_1based)
        .saturating_sub(exon.end_1based.min(exon.start_1based))
        + 1
}

fn splicing_exon_mod3(exon: &SplicingExonSummary) -> usize {
    splicing_exon_length_bp(exon) % 3
}

fn splicing_exon_mod3_fill_hex(mod3: usize) -> &'static str {
    match mod3 % 3 {
        0 => "#dbeafe",
        1 => "#fef3c7",
        _ => "#ffe4e6",
    }
}

fn splicing_exon_mod3_text_hex(mod3: usize) -> &'static str {
    match mod3 % 3 {
        0 => "#1e40af",
        1 => "#92400e",
        _ => "#9f1239",
    }
}

fn splicing_cds_phase_fill_hex(phase: u8) -> &'static str {
    match phase % 3 {
        0 => "#2563eb",
        1 => "#d97706",
        _ => "#e11d48",
    }
}

fn render_tfbs(view: &TfbsExpertView) -> String {
    let mut doc = Document::new()
        .set("viewBox", (0, 0, W, H))
        .set("width", W)
        .set("height", H)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", W)
                .set("height", H)
                .set("fill", "#ffffff"),
        );

    let left = 84.0;
    let right = W - 54.0;
    let chart_top = 120.0;
    let chart_bottom = H - 190.0;
    let chart_height = chart_bottom - chart_top;
    let baseline = chart_bottom;
    let n = view.columns.len().max(1);
    let available = (right - left).max(1.0);
    let step = (available / n as f32).max(8.0);
    let bar_width = (step * 0.72).clamp(5.0, 22.0);
    let scale = chart_height / 2.0;

    doc = doc.add(
        Text::new(format!(
            "TFBS expert: {} ({}) on {}:{}..{} ({})",
            view.tf_name.clone().unwrap_or_else(|| view.tf_id.clone()),
            view.tf_id,
            view.seq_id,
            view.start_1based,
            view.end_1based,
            view.strand
        ))
        .set("x", left)
        .set("y", 36)
        .set("font-family", "monospace")
        .set("font-size", 18)
        .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new(format!(
            "matched={} | motif_length={} | llr_bits={:.4} | true_log_odds_bits={:.4}",
            view.matched_sequence,
            view.motif_length,
            view.llr_total_bits.unwrap_or_default(),
            view.true_log_odds_total_bits.unwrap_or_default()
        ))
        .set("x", left)
        .set("y", 58)
        .set("font-family", "monospace")
        .set("font-size", 13)
        .set("fill", "#444444"),
    );

    doc = doc.add(
        Line::new()
            .set("x1", left)
            .set("y1", baseline)
            .set("x2", right)
            .set("y2", baseline)
            .set("stroke", "#333333")
            .set("stroke-width", 1),
    );
    for tick in [0.0_f32, 1.0, 2.0] {
        let y = baseline - tick * scale;
        doc = doc.add(
            Line::new()
                .set("x1", left - 6.0)
                .set("y1", y)
                .set("x2", right)
                .set("y2", y)
                .set("stroke", if tick == 0.0 { "#333333" } else { "#e5e7eb" })
                .set("stroke-width", 1),
        );
        doc = doc.add(
            Text::new(format!("{tick:.0}"))
                .set("x", left - 10.0)
                .set("y", y + 4.0)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#555555"),
        );
    }
    doc = doc.add(
        Text::new("information content (bits)")
            .set("x", left - 48.0)
            .set("y", chart_top - 10.0)
            .set("font-family", "monospace")
            .set("font-size", 11)
            .set("fill", "#555555"),
    );

    let mut match_points: Vec<(f32, f32)> = Vec::new();

    for (idx, col) in view.columns.iter().enumerate() {
        let x_center = left + step * (idx as f32 + 0.5);
        let x = x_center - bar_width / 2.0;
        let ic = col.information_content_bits.clamp(0.0, 2.0) as f32;
        let full_height = ic * scale;
        let mut y_cursor = baseline;
        let match_idx = match_base_to_idx(col.match_base);
        let mut match_y: Option<f32> = None;

        for base_idx in 0..4 {
            let frac = col.frequencies[base_idx].clamp(0.0, 1.0) as f32;
            let segment_h = (full_height * frac).max(0.0);
            if segment_h <= 0.0 {
                continue;
            }
            let y_top = y_cursor - segment_h;
            doc = doc.add(
                Rectangle::new()
                    .set("x", x)
                    .set("y", y_top)
                    .set("width", bar_width)
                    .set("height", segment_h)
                    .set("fill", nuc_color(base_idx)),
            );
            if match_idx == Some(base_idx) {
                match_y = Some(y_top + segment_h * 0.5);
            }
            y_cursor = y_top;
        }

        if let Some(y) = match_y {
            match_points.push((x_center, y));
        }

        if idx % 2 == 0 || n <= 24 {
            doc = doc.add(
                Text::new(format!("{}", col.index_1based))
                    .set("x", x_center)
                    .set("y", baseline + 16.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 10)
                    .set("fill", "#4b5563"),
            );
        }

        if let Some(base) = col.match_base {
            doc = doc.add(
                Text::new(base.to_ascii_uppercase().to_string())
                    .set("x", x_center)
                    .set("y", baseline + 32.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 11)
                    .set("fill", "#111111"),
            );
        }
    }

    if match_points.len() >= 2 {
        let mut data = Data::new().move_to((match_points[0].0, match_points[0].1));
        for (x, y) in match_points.iter().skip(1) {
            data = data.line_to((*x, *y));
        }
        doc = doc.add(
            Path::new()
                .set("d", data)
                .set("fill", "none")
                .set("stroke", "#111111")
                .set("stroke-width", 2.0),
        );
        for (x, y) in match_points {
            doc = doc.add(
                Circle::new()
                    .set("cx", x)
                    .set("cy", y)
                    .set("r", 2.3)
                    .set("fill", "#111111"),
            );
        }
    }

    let legend_items = [("A", 0usize), ("C", 1usize), ("G", 2usize), ("T", 3usize)];
    let mut lx = left;
    let ly = H - 132.0;
    for (label, idx) in legend_items {
        doc = doc.add(
            Rectangle::new()
                .set("x", lx)
                .set("y", ly - 9.0)
                .set("width", 12)
                .set("height", 12)
                .set("fill", nuc_color(idx)),
        );
        doc = doc.add(
            Text::new(label)
                .set("x", lx + 16.0)
                .set("y", ly + 1.0)
                .set("font-family", "monospace")
                .set("font-size", 12)
                .set("fill", "#111111"),
        );
        lx += 52.0;
    }
    doc = doc.add(
        Text::new("black line = matched sequence path across motif columns")
            .set("x", left + 230.0)
            .set("y", ly + 1.0)
            .set("font-family", "monospace")
            .set("font-size", 11)
            .set("fill", "#4b5563"),
    );

    for (line_idx, line) in wrap_text(&view.instruction, 125).into_iter().enumerate() {
        doc = doc.add(
            Text::new(line)
                .set("x", left)
                .set("y", H - 86.0 + line_idx as f32 * 14.0)
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#374151"),
        );
    }

    doc.to_string()
}

fn restriction_end_geometry_display_label(end_geometry: &str, overlap_bp: Option<isize>) -> String {
    match end_geometry {
        "5prime_overhang" => format!(
            "5' overhang ({} bp)",
            overlap_bp.unwrap_or(0).unsigned_abs()
        ),
        "3prime_overhang" => format!(
            "3' overhang ({} bp)",
            overlap_bp.unwrap_or(0).unsigned_abs()
        ),
        _ => match overlap_bp {
            Some(value) if value > 0 => format!("5' overhang ({} bp)", value as usize),
            Some(value) if value < 0 => format!("3' overhang ({} bp)", value.unsigned_abs()),
            _ => "blunt".to_string(),
        },
    }
}

fn restriction_end_geometry_color(end_geometry: &str, overlap_bp: Option<isize>) -> &'static str {
    match end_geometry {
        "5prime_overhang" => "#2563eb",
        "3prime_overhang" => "#d97706",
        _ => match overlap_bp {
            Some(value) if value > 0 => "#2563eb",
            Some(value) if value < 0 => "#d97706",
            _ => "#475569",
        },
    }
}

fn render_restriction(view: &RestrictionSiteExpertView) -> String {
    let mut doc = Document::new()
        .set("viewBox", (0, 0, W, H))
        .set("width", W)
        .set("height", H)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", W)
                .set("height", H)
                .set("fill", "#ffffff"),
        );

    let title = format!(
        "Restriction expert: {} @ {}:{}..{}",
        view.selected_enzyme
            .clone()
            .unwrap_or_else(|| view.enzyme_names.join(",")),
        view.seq_id,
        view.recognition_start_1based,
        view.recognition_end_1based
    );
    doc = doc.add(
        Text::new(title)
            .set("x", 90)
            .set("y", 42)
            .set("font-family", "monospace")
            .set("font-size", 18)
            .set("fill", "#111111"),
    );
    let geometry_label =
        restriction_end_geometry_display_label(&view.end_geometry, view.overlap_bp);
    let cut_color = restriction_end_geometry_color(&view.end_geometry, view.overlap_bp);
    let paired_cut_pos_1based = if view.paired_cut_pos_1based == 0 {
        view.cut_pos_1based
    } else {
        view.paired_cut_pos_1based
    };
    let paired_cut_index_0based =
        if view.paired_cut_pos_1based == 0 && view.paired_cut_index_0based == 0 {
            view.cut_index_0based
        } else {
            view.paired_cut_index_0based
        };
    doc = doc.add(
        Text::new(format!(
            "cuts={}|{} | geometry={} | cuts_for_enzyme={} | recognition_iupac={}",
            view.cut_pos_1based,
            paired_cut_pos_1based,
            geometry_label,
            view.number_of_cuts_for_enzyme,
            view.recognition_iupac
                .clone()
                .unwrap_or_else(|| "-".to_string())
        ))
        .set("x", 90)
        .set("y", 66)
        .set("font-family", "monospace")
        .set("font-size", 12)
        .set("fill", "#4b5563"),
    );
    if view.enzyme_cut_offset_0based.is_some() || view.overlap_bp.is_some() {
        let mut fields = Vec::new();
        if let Some(cut_offset) = view.enzyme_cut_offset_0based {
            fields.push(format!("enzyme_cut_offset_0based={cut_offset}"));
        }
        if let Some(overlap) = view.overlap_bp {
            fields.push(format!("overlap_bp={overlap}"));
        }
        if !fields.is_empty() {
            doc = doc.add(
                Text::new(fields.join(" | "))
                    .set("x", 90)
                    .set("y", 84)
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("fill", "#4b5563"),
            );
        }
    }
    if let Some(note) = &view.enzyme_note {
        for (line_idx, line) in wrap_text(&format!("note={note}"), 125)
            .into_iter()
            .enumerate()
        {
            doc = doc.add(
                Text::new(line)
                    .set("x", 90)
                    .set("y", 102.0 + line_idx as f32 * 14.0)
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("fill", "#4b5563"),
            );
        }
    }
    if let Some(url) = &view.rebase_url {
        doc = doc.add(
            Text::new(format!("rebase_url={url}"))
                .set("x", 90)
                .set("y", 120)
                .set("font-family", "monospace")
                .set("font-size", 12)
                .set("fill", "#2563eb"),
        );
    }

    let top = if view.site_sequence.is_empty() {
        view.recognition_iupac.clone().unwrap_or_default()
    } else {
        view.site_sequence.clone()
    };
    let bottom = if view.site_sequence_complement.is_empty() {
        String::new()
    } else {
        view.site_sequence_complement.clone()
    };

    let n = top.chars().count().max(1);
    let step = ((W - 220.0) / n as f32).clamp(16.0, 28.0);
    let total = step * n as f32;
    let start_x = (W - total) * 0.5;
    let top_y = 260.0;
    let bottom_y = 320.0;
    let rail_left = start_x - 18.0;
    let rail_right = start_x + total + 18.0;

    doc = doc.add(
        Line::new()
            .set("x1", rail_left)
            .set("y1", top_y - 10.0)
            .set("x2", rail_right)
            .set("y2", top_y - 10.0)
            .set("stroke", "#111111")
            .set("stroke-width", 2),
    );
    doc = doc.add(
        Line::new()
            .set("x1", rail_left)
            .set("y1", bottom_y + 10.0)
            .set("x2", rail_right)
            .set("y2", bottom_y + 10.0)
            .set("stroke", "#111111")
            .set("stroke-width", 2),
    );
    doc = doc.add(
        Text::new("5'")
            .set("x", rail_left - 32.0)
            .set("y", top_y - 6.0)
            .set("font-family", "monospace")
            .set("font-size", 14)
            .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new("3'")
            .set("x", rail_right + 10.0)
            .set("y", top_y - 6.0)
            .set("font-family", "monospace")
            .set("font-size", 14)
            .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new("3'")
            .set("x", rail_left - 32.0)
            .set("y", bottom_y + 14.0)
            .set("font-family", "monospace")
            .set("font-size", 14)
            .set("fill", "#111111"),
    );
    doc = doc.add(
        Text::new("5'")
            .set("x", rail_right + 10.0)
            .set("y", bottom_y + 14.0)
            .set("font-family", "monospace")
            .set("font-size", 14)
            .set("fill", "#111111"),
    );

    for (idx, ch) in top.chars().enumerate() {
        let x = start_x + step * (idx as f32 + 0.5);
        doc = doc.add(
            Text::new(ch.to_string())
                .set("x", x)
                .set("y", top_y)
                .set("text-anchor", "middle")
                .set("font-family", "monospace")
                .set("font-size", 22)
                .set("fill", "#0f172a"),
        );
    }
    for (idx, ch) in bottom.chars().enumerate() {
        let x = start_x + step * (idx as f32 + 0.5);
        doc = doc.add(
            Text::new(ch.to_string())
                .set("x", x)
                .set("y", bottom_y)
                .set("text-anchor", "middle")
                .set("font-family", "monospace")
                .set("font-size", 22)
                .set("fill", "#334155"),
        );
    }

    let top_cut_x = start_x + step * view.cut_index_0based as f32;
    let bottom_cut_x = start_x + step * paired_cut_index_0based as f32;
    if (top_cut_x - bottom_cut_x).abs() < f32::EPSILON {
        doc = doc.add(
            Line::new()
                .set("x1", top_cut_x)
                .set("y1", top_y - 42.0)
                .set("x2", top_cut_x)
                .set("y2", bottom_y + 32.0)
                .set("stroke", cut_color)
                .set("stroke-width", 2.5),
        );
    } else {
        doc = doc.add(
            Line::new()
                .set("x1", top_cut_x)
                .set("y1", top_y - 42.0)
                .set("x2", top_cut_x)
                .set("y2", top_y - 6.0)
                .set("stroke", cut_color)
                .set("stroke-width", 2.5),
        );
        doc = doc.add(
            Line::new()
                .set("x1", top_cut_x)
                .set("y1", top_y - 6.0)
                .set("x2", bottom_cut_x)
                .set("y2", bottom_y + 6.0)
                .set("stroke", cut_color)
                .set("stroke-width", 2.5),
        );
        doc = doc.add(
            Line::new()
                .set("x1", bottom_cut_x)
                .set("y1", bottom_y + 6.0)
                .set("x2", bottom_cut_x)
                .set("y2", bottom_y + 32.0)
                .set("stroke", cut_color)
                .set("stroke-width", 2.5),
        );
    }
    doc = doc.add(
        Text::new(geometry_label)
            .set("x", top_cut_x.max(bottom_cut_x) + 4.0)
            .set("y", top_y - 46.0)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", cut_color),
    );

    for (line_idx, line) in wrap_text(&view.instruction, 125).into_iter().enumerate() {
        doc = doc.add(
            Text::new(line)
                .set("x", 90)
                .set("y", H - 84.0 + line_idx as f32 * 14.0)
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#374151"),
        );
    }

    doc.to_string()
}

fn render_splicing(view: &SplicingExpertView) -> String {
    let lane_count = view.transcripts.len().max(1);
    let unique_exon_total = view.unique_exons.len();
    let exon_count = unique_exon_total.max(1);
    let transcript_total = view.transcript_count.max(1);
    let exon_transitions = compute_splicing_exon_transition_matrix(view);
    let lane_height = 30.0_f32;
    let lane_gap = 14.0_f32;
    let chart_top = 126.0_f32;
    let chart_height =
        lane_count as f32 * lane_height + (lane_count.saturating_sub(1) as f32) * lane_gap;
    let matrix_top = chart_top + chart_height + 78.0;
    let matrix_row_h = 17.0_f32;
    let matrix_support_row_h = 16.0_f32;
    let matrix_mod_row_h = 15.0_f32;
    let matrix_h =
        lane_count as f32 * matrix_row_h + matrix_support_row_h + matrix_mod_row_h + 38.0;
    let transition_top = matrix_top + matrix_h + 24.0;
    let transition_h = if unique_exon_total >= 2 {
        45.0 + unique_exon_total as f32 * 14.0
    } else {
        0.0
    };
    let junction_table_top = transition_top + transition_h + 24.0;
    let mut junction_rows = view.junctions.clone();
    junction_rows.sort_by(|left, right| {
        right
            .support_transcript_count
            .cmp(&left.support_transcript_count)
            .then_with(|| left.donor_1based.cmp(&right.donor_1based))
            .then_with(|| left.acceptor_1based.cmp(&right.acceptor_1based))
    });
    let junction_rows_rendered: Vec<_> = junction_rows.iter().take(24).cloned().collect();
    let junction_h = if junction_rows_rendered.is_empty() {
        0.0
    } else {
        let base = 34.0 + junction_rows_rendered.len() as f32 * 14.0;
        if view.junctions.len() > junction_rows_rendered.len() {
            base + 13.0
        } else {
            base
        }
    };
    let footer_top = junction_table_top + junction_h + 28.0;
    let dyn_h = (footer_top + 126.0).max(H);
    let left = 240.0_f32;
    let right = W - 56.0_f32;
    let map_width = (right - left).max(1.0);
    let region_start = view.region_start_1based.max(1);
    let region_end = view.region_end_1based.max(region_start);
    let region_span = (region_end - region_start + 1).max(1) as f32;
    let x_for = |pos_1based: usize| -> f32 {
        let clamped = pos_1based.clamp(region_start, region_end);
        let rel = (clamped.saturating_sub(region_start)) as f32 / region_span;
        left + rel * map_width
    };

    let exon_meta: HashMap<(usize, usize), (usize, bool)> = view
        .unique_exons
        .iter()
        .map(|exon| {
            (
                (exon.start_1based, exon.end_1based),
                (exon.support_transcript_count, exon.constitutive),
            )
        })
        .collect();

    let mut doc = Document::new()
        .set("viewBox", (0, 0, W, dyn_h))
        .set("width", W)
        .set("height", dyn_h)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", W)
                .set("height", dyn_h)
                .set("fill", "#ffffff"),
        );

    doc = doc
        .add(
            Text::new(format!(
                "Splicing expert: {} ({}) on {}:{}..{} ({})",
                view.group_label,
                view.transcript_count,
                view.seq_id,
                view.region_start_1based,
                view.region_end_1based,
                view.strand
            ))
            .set("x", 88)
            .set("y", 38)
            .set("font-family", "monospace")
            .set("font-size", 18)
            .set("fill", "#111111"),
        )
        .add(
            Text::new(format!(
                "target_feature={} | unique_exons={} | donor/acceptor markers shown",
                view.target_feature_id, view.unique_exon_count
            ))
            .set("x", 88)
            .set("y", 60)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#4b5563"),
        )
        .add(
            Text::new("CDS flank phase colors: 0=blue 1=amber 2=rose (drawn on exon left/right edges when CDS ranges exist)")
                .set("x", 88)
                .set("y", 76)
                .set("font-family", "monospace")
                .set("font-size", 9)
                .set("fill", "#64748b"),
        )
        .add(
            Line::new()
                .set("x1", left)
                .set("y1", chart_top - 16.0)
                .set("x2", right)
                .set("y2", chart_top - 16.0)
                .set("stroke", "#6b7280")
                .set("stroke-width", 1.0),
        )
        .add(
            Text::new(format!("{} bp", region_start))
                .set("x", left)
                .set("y", chart_top - 21.0)
                .set("text-anchor", "start")
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#4b5563"),
        )
        .add(
            Text::new(format!("{} bp", region_end))
                .set("x", right)
                .set("y", chart_top - 21.0)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#4b5563"),
        );

    for junction in &view.junctions {
        let x1 = x_for(junction.donor_1based);
        let x2 = x_for(junction.acceptor_1based);
        if (x2 - x1).abs() < 2.0 {
            continue;
        }
        let mid_x = (x1 + x2) * 0.5;
        let span = (x2 - x1).abs().max(1.0);
        let lift = (22.0 + span * 0.04).clamp(24.0, 84.0);
        let apex_y = chart_top - lift;
        let stroke_w = (1.0 + junction.support_transcript_count as f32 * 0.35).clamp(1.2, 4.4);
        let path = Data::new()
            .move_to((x1, chart_top - 2.0))
            .quadratic_curve_to((mid_x, apex_y, x2, chart_top - 2.0));
        doc = doc.add(
            Path::new()
                .set("d", path)
                .set("fill", "none")
                .set("stroke", "#64748b")
                .set("stroke-width", stroke_w),
        );
        let show_support_label =
            view.junctions.len() <= 48 || junction.support_transcript_count > 1;
        if show_support_label {
            doc = doc.add(
                Text::new(format!("{}", junction.support_transcript_count))
                    .set("x", mid_x)
                    .set("y", apex_y - 4.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 8)
                    .set("fill", "#1d4ed8"),
            );
        }
    }

    for (lane_idx, lane) in view.transcripts.iter().enumerate() {
        let y = chart_top + lane_idx as f32 * (lane_height + lane_gap) + lane_height * 0.5;
        if lane.has_target_feature {
            doc = doc.add(
                Rectangle::new()
                    .set("x", left - 4.0)
                    .set("y", y - lane_height * 0.5 - 3.0)
                    .set("width", map_width + 8.0)
                    .set("height", lane_height + 6.0)
                    .set("fill", "#fef9c3")
                    .set("fill-opacity", 0.35),
            );
        }
        doc = doc.add(
            Text::new(format!(
                "{} ({})",
                lane.transcript_id,
                if lane.strand.trim().is_empty() {
                    "?"
                } else {
                    lane.strand.as_str()
                }
            ))
            .set("x", left - 10.0)
            .set("y", y + 4.0)
            .set("text-anchor", "end")
            .set("font-family", "monospace")
            .set("font-size", 11)
            .set(
                "fill",
                if lane.has_target_feature {
                    "#854d0e"
                } else {
                    "#111827"
                },
            ),
        );

        for intron in &lane.introns {
            let x1 = x_for(intron.start_1based);
            let x2 = x_for(intron.end_1based);
            if x2 <= x1 {
                continue;
            }
            doc = doc.add(
                Line::new()
                    .set("x1", x1)
                    .set("y1", y)
                    .set("x2", x2)
                    .set("y2", y)
                    .set("stroke", "#6b7280")
                    .set("stroke-width", 1.2),
            );
        }

        for (exon_idx, exon) in lane.exons.iter().enumerate() {
            let x1 = x_for(exon.start_1based);
            let x2 = x_for(exon.end_1based);
            if x2 <= x1 {
                continue;
            }
            let (support, constitutive) = exon_meta
                .get(&(exon.start_1based, exon.end_1based))
                .copied()
                .unwrap_or((1, false));
            let fill = if constitutive { "#2563eb" } else { "#b45309" };
            let opacity = if constitutive {
                0.92
            } else if support <= 1 {
                0.65
            } else {
                0.78
            };
            doc = doc
                .add(
                    Rectangle::new()
                        .set("x", x1)
                        .set("y", y - 9.5)
                        .set("width", (x2 - x1).max(1.0))
                        .set("height", 19.0)
                        .set("fill", fill)
                        .set("fill-opacity", opacity)
                        .set("stroke", "#1f2937")
                        .set("stroke-width", 0.6),
                )
                .add(
                    Text::new(format!("{support}"))
                        .set("x", x1 + 2.0)
                        .set("y", y - 11.0)
                        .set("font-family", "monospace")
                        .set("font-size", 8)
                        .set("fill", "#334155"),
                );
            let phase_info = lane
                .exon_cds_phases
                .get(exon_idx)
                .filter(|phase| {
                    phase.start_1based == exon.start_1based && phase.end_1based == exon.end_1based
                })
                .or_else(|| {
                    lane.exon_cds_phases.iter().find(|phase| {
                        phase.start_1based == exon.start_1based
                            && phase.end_1based == exon.end_1based
                    })
                });
            if let Some(phase_info) = phase_info {
                let flank_w = ((x2 - x1) * 0.15).clamp(1.5, 3.5);
                if let Some(left_phase) = phase_info.left_cds_phase {
                    doc = doc.add(
                        Rectangle::new()
                            .set("x", x1)
                            .set("y", y - 9.5)
                            .set("width", flank_w)
                            .set("height", 19.0)
                            .set("fill", splicing_cds_phase_fill_hex(left_phase)),
                    );
                }
                if let Some(right_phase) = phase_info.right_cds_phase {
                    doc = doc.add(
                        Rectangle::new()
                            .set("x", (x2 - flank_w).max(x1))
                            .set("y", y - 9.5)
                            .set("width", flank_w)
                            .set("height", 19.0)
                            .set("fill", splicing_cds_phase_fill_hex(right_phase)),
                    );
                }
            }
        }
    }

    for marker in &view.boundaries {
        let lane_idx = view
            .transcripts
            .iter()
            .position(|lane| lane.transcript_feature_id == marker.transcript_feature_id);
        let Some(lane_idx) = lane_idx else {
            continue;
        };
        let y = chart_top + lane_idx as f32 * (lane_height + lane_gap) + lane_height * 0.5;
        let x = x_for(marker.position_1based);
        let (tick_top, tick_bottom, color, dy) = if marker.side.eq_ignore_ascii_case("donor") {
            (y - 18.0, y - 8.0, "#be123c", -20.0)
        } else {
            (y + 8.0, y + 18.0, "#0f766e", 28.0)
        };
        doc = doc
            .add(
                Line::new()
                    .set("x1", x)
                    .set("y1", tick_top)
                    .set("x2", x)
                    .set("y2", tick_bottom)
                    .set("stroke", color)
                    .set("stroke-width", 1.4),
            )
            .add(
                Text::new(format!(
                    "{}:{}{}",
                    marker.side,
                    marker.motif_2bp,
                    if marker.canonical { "" } else { "*" }
                ))
                .set("x", x + 2.0)
                .set("y", y + dy)
                .set("font-family", "monospace")
                .set("font-size", 8)
                .set("fill", color),
            );
    }

    let matrix_label_x = 88.0_f32;
    let matrix_left = 250.0_f32;
    let matrix_cell_w = ((right - matrix_left - 6.0) / exon_count as f32).clamp(7.0, 16.0);
    let matrix_cell_h = 12.0_f32;
    let matrix_support_y = matrix_top + 12.0;
    let matrix_mod_y = matrix_support_y + matrix_support_row_h;
    let matrix_rows_top = matrix_mod_y + matrix_mod_row_h - 1.0;
    doc = doc.add(
        Text::new("Transcript vs exon matrix")
            .set("x", matrix_label_x)
            .set("y", matrix_top - 14.0)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#111827"),
    );
    doc = doc.add(
        Text::new("cell color intensity encodes exon support frequency")
            .set("x", matrix_left)
            .set("y", matrix_top - 14.0)
            .set("font-family", "monospace")
            .set("font-size", 9)
            .set("fill", "#475569"),
    );
    for (col_idx, exon) in view.unique_exons.iter().enumerate() {
        let mod3 = splicing_exon_mod3(exon);
        let x = matrix_left + col_idx as f32 * matrix_cell_w;
        doc = doc
            .add(
                Rectangle::new()
                    .set("x", x)
                    .set("y", matrix_top - 10.5)
                    .set("width", (matrix_cell_w - 1.0).max(4.0))
                    .set("height", 9.5)
                    .set("fill", splicing_exon_mod3_fill_hex(mod3))
                    .set("stroke", "#cbd5e1")
                    .set("stroke-width", 0.4),
            )
            .add(
                Text::new(format!("E{}", col_idx + 1))
                    .set("x", x + (matrix_cell_w - 1.0) * 0.5)
                    .set("y", matrix_top - 2.2)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 7)
                    .set("fill", splicing_exon_mod3_text_hex(mod3)),
            );
    }
    doc = doc.add(
        Text::new("support")
            .set("x", matrix_label_x)
            .set("y", matrix_support_y)
            .set("font-family", "monospace")
            .set("font-size", 9)
            .set("fill", "#334155"),
    );
    doc = doc.add(
        Text::new("len%3")
            .set("x", matrix_label_x)
            .set("y", matrix_mod_y)
            .set("font-family", "monospace")
            .set("font-size", 9)
            .set("fill", "#334155"),
    );
    for (col_idx, exon) in view.unique_exons.iter().enumerate() {
        let x = matrix_left + col_idx as f32 * matrix_cell_w + 0.5;
        let support_ratio =
            support_ratio_percent(exon.support_transcript_count, transcript_total) as f32 / 100.0;
        let label = if exon.constitutive {
            format!(
                "{} const",
                format_support_fraction(exon.support_transcript_count, transcript_total)
            )
        } else {
            format_support_fraction(exon.support_transcript_count, transcript_total)
        };
        doc = doc.add(
            Text::new(label)
                .set("x", x)
                .set("y", matrix_support_y)
                .set("font-family", "monospace")
                .set("font-size", 7)
                .set(
                    "fill",
                    rgb_hex(mix_rgb([100, 116, 139], [30, 64, 175], support_ratio)),
                ),
        );
        let mod3 = splicing_exon_mod3(exon);
        doc = doc.add(
            Text::new(format!("{mod3}"))
                .set("x", x + (matrix_cell_w - 1.0) * 0.5)
                .set("y", matrix_mod_y)
                .set("text-anchor", "middle")
                .set("font-family", "monospace")
                .set("font-size", 8)
                .set("fill", splicing_exon_mod3_text_hex(mod3)),
        );
    }
    for (row_idx, row) in view.matrix_rows.iter().enumerate() {
        let y = matrix_rows_top + row_idx as f32 * matrix_row_h;
        doc = doc.add(
            Text::new(row.transcript_id.clone())
                .set("x", matrix_label_x)
                .set("y", y + 10.0)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#1f2937"),
        );
        for (col_idx, present) in row.exon_presence.iter().copied().enumerate() {
            let x = matrix_left + col_idx as f32 * matrix_cell_w;
            let support_count = view
                .unique_exons
                .get(col_idx)
                .map(|exon| exon.support_transcript_count)
                .unwrap_or(0);
            let support_ratio =
                support_ratio_percent(support_count, transcript_total) as f32 / 100.0;
            let fill = splicing_matrix_cell_fill_hex(present, support_ratio);
            doc = doc.add(
                Rectangle::new()
                    .set("x", x)
                    .set("y", y)
                    .set("width", matrix_cell_w - 1.0)
                    .set("height", matrix_cell_h)
                    .set("fill", fill)
                    .set("stroke", "#cbd5e1")
                    .set("stroke-width", 0.5),
            );
            doc = doc.add(
                Text::new(if present { "X" } else { "·" })
                    .set("x", x + (matrix_cell_w - 1.0) * 0.5)
                    .set("y", y + matrix_cell_h - 2.0)
                    .set("text-anchor", "middle")
                    .set("font-family", "monospace")
                    .set("font-size", 8)
                    .set(
                        "fill",
                        splicing_matrix_cell_text_hex(present, support_ratio),
                    ),
            );
        }
    }

    if unique_exon_total >= 2 {
        let transition_label_x = 88.0_f32;
        let transition_axis_w = 58.0_f32;
        let transition_matrix_left = 250.0_f32 + transition_axis_w;
        let transition_cell_w =
            ((right - transition_matrix_left - 6.0) / exon_count as f32).clamp(7.0, 16.0);
        let transition_cell_h = 11.0_f32;
        let transition_header_y = transition_top + 13.0;
        let transition_rows_top = transition_top + 22.0;
        doc = doc
            .add(
                Text::new("Exon -> exon transition matrix")
                    .set("x", transition_label_x)
                    .set("y", transition_top)
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("fill", "#111827"),
            )
            .add(
                Text::new(
                    "frequency color: transition support; E# color: exon genomic length %3 (heuristic frame cue)",
                )
                .set("x", transition_matrix_left - transition_axis_w)
                .set("y", transition_top + 12.5)
                .set("font-family", "monospace")
                .set("font-size", 8)
                .set("fill", "#475569"),
            );
        for (to_idx, exon) in view.unique_exons.iter().enumerate() {
            let mod3 = splicing_exon_mod3(exon);
            let x = transition_matrix_left + to_idx as f32 * transition_cell_w;
            doc = doc
                .add(
                    Rectangle::new()
                        .set("x", x)
                        .set("y", transition_header_y - 8.5)
                        .set("width", (transition_cell_w - 1.0).max(4.0))
                        .set("height", 9.5)
                        .set("fill", splicing_exon_mod3_fill_hex(mod3))
                        .set("stroke", "#cbd5e1")
                        .set("stroke-width", 0.4),
                )
                .add(
                    Text::new(format!("E{}", to_idx + 1))
                        .set("x", x + (transition_cell_w - 1.0) * 0.5)
                        .set("y", transition_header_y - 1.0)
                        .set("text-anchor", "middle")
                        .set("font-family", "monospace")
                        .set("font-size", 7)
                        .set("fill", splicing_exon_mod3_text_hex(mod3)),
                );
        }
        for (from_idx, exon) in view.unique_exons.iter().enumerate() {
            let y = transition_rows_top + from_idx as f32 * 14.0;
            let mod3 = splicing_exon_mod3(exon);
            doc = doc
                .add(
                    Rectangle::new()
                        .set("x", transition_matrix_left - transition_axis_w)
                        .set("y", y)
                        .set("width", transition_axis_w - 4.0)
                        .set("height", transition_cell_h)
                        .set("fill", splicing_exon_mod3_fill_hex(mod3))
                        .set("stroke", "#cbd5e1")
                        .set("stroke-width", 0.4),
                )
                .add(
                    Text::new(format!("E{}", from_idx + 1))
                        .set("x", transition_matrix_left - transition_axis_w + 5.0)
                        .set("y", y + transition_cell_h - 2.0)
                        .set("font-family", "monospace")
                        .set("font-size", 8)
                        .set("fill", splicing_exon_mod3_text_hex(mod3)),
                );
            for to_idx in 0..view.unique_exons.len() {
                let support_count = exon_transitions
                    .counts
                    .get(from_idx)
                    .and_then(|row| row.get(to_idx))
                    .copied()
                    .unwrap_or(0);
                let support_ratio =
                    support_ratio_percent(support_count, transcript_total) as f32 / 100.0;
                let fill = splicing_transition_cell_fill_hex(support_ratio);
                let x = transition_matrix_left + to_idx as f32 * transition_cell_w;
                doc = doc
                    .add(
                        Rectangle::new()
                            .set("x", x)
                            .set("y", y)
                            .set("width", transition_cell_w - 1.0)
                            .set("height", transition_cell_h)
                            .set("fill", fill)
                            .set("stroke", "#cbd5e1")
                            .set("stroke-width", 0.5),
                    )
                    .add(
                        Text::new(if support_count > 0 {
                            support_count.to_string()
                        } else {
                            "·".to_string()
                        })
                        .set("x", x + (transition_cell_w - 1.0) * 0.5)
                        .set("y", y + transition_cell_h - 2.0)
                        .set("text-anchor", "middle")
                        .set("font-family", "monospace")
                        .set("font-size", 8)
                        .set("fill", splicing_transition_cell_text_hex(support_ratio)),
                    );
            }
        }
    }

    if !junction_rows_rendered.is_empty() {
        let table_left = 88.0_f32;
        let donor_x = table_left;
        let acceptor_x = table_left + 88.0;
        let span_x = table_left + 176.0;
        let support_x = table_left + 246.0;
        let transcripts_x = table_left + 362.0;
        let mut y = junction_table_top;
        doc = doc
            .add(
                Text::new("Junction transition support")
                    .set("x", table_left)
                    .set("y", y)
                    .set("font-family", "monospace")
                    .set("font-size", 12)
                    .set("fill", "#111827"),
            )
            .add(
                Text::new("donor")
                    .set("x", donor_x)
                    .set("y", y + 13.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#334155"),
            )
            .add(
                Text::new("acceptor")
                    .set("x", acceptor_x)
                    .set("y", y + 13.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#334155"),
            )
            .add(
                Text::new("span")
                    .set("x", span_x)
                    .set("y", y + 13.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#334155"),
            )
            .add(
                Text::new("support")
                    .set("x", support_x)
                    .set("y", y + 13.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#334155"),
            )
            .add(
                Text::new("transcripts")
                    .set("x", transcripts_x)
                    .set("y", y + 13.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#334155"),
            );
        y += 26.0;

        for junction in &junction_rows_rendered {
            let span_bp = junction
                .acceptor_1based
                .max(junction.donor_1based)
                .saturating_sub(junction.acceptor_1based.min(junction.donor_1based));
            let mut transcript_labels = junction
                .transcript_feature_ids
                .iter()
                .map(|feature_id| format!("n-{feature_id}"))
                .collect::<Vec<_>>();
            if transcript_labels.len() > 10 {
                let hidden = transcript_labels.len() - 10;
                transcript_labels.truncate(10);
                transcript_labels.push(format!("+{hidden}"));
            }
            doc = doc
                .add(
                    Text::new(junction.donor_1based.to_string())
                        .set("x", donor_x)
                        .set("y", y)
                        .set("font-family", "monospace")
                        .set("font-size", 9)
                        .set("fill", "#1f2937"),
                )
                .add(
                    Text::new(junction.acceptor_1based.to_string())
                        .set("x", acceptor_x)
                        .set("y", y)
                        .set("font-family", "monospace")
                        .set("font-size", 9)
                        .set("fill", "#1f2937"),
                )
                .add(
                    Text::new(span_bp.to_string())
                        .set("x", span_x)
                        .set("y", y)
                        .set("font-family", "monospace")
                        .set("font-size", 9)
                        .set("fill", "#1f2937"),
                )
                .add(
                    Text::new(format_support_fraction(
                        junction.support_transcript_count,
                        transcript_total,
                    ))
                    .set("x", support_x)
                    .set("y", y)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#1f2937"),
                )
                .add(
                    Text::new(transcript_labels.join(", "))
                        .set("x", transcripts_x)
                        .set("y", y)
                        .set("font-family", "monospace")
                        .set("font-size", 8)
                        .set("fill", "#475569"),
                );
            y += 13.0;
        }
        if view.junctions.len() > junction_rows_rendered.len() {
            let hidden = view.junctions.len() - junction_rows_rendered.len();
            doc = doc.add(
                Text::new(format!("... {} additional transitions omitted", hidden))
                    .set("x", table_left)
                    .set("y", y)
                    .set("font-family", "monospace")
                    .set("font-size", 8)
                    .set("fill", "#6b7280"),
            );
        }
    }

    let mut event_y = footer_top;
    doc = doc.add(
        Text::new("Event summary")
            .set("x", 88.0)
            .set("y", event_y)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#111827"),
    );
    event_y += 14.0;
    for event in &view.events {
        doc = doc.add(
            Text::new(format!("{}: {}", event.event_type, event.count))
                .set("x", 88.0)
                .set("y", event_y)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#334155"),
        );
        event_y += 12.0;
        for detail in event.details.iter().take(2) {
            doc = doc.add(
                Text::new(format!("  - {detail}"))
                    .set("x", 88.0)
                    .set("y", event_y)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#475569"),
            );
            event_y += 11.0;
        }
    }

    for (line_idx, line) in wrap_text(&view.instruction, 140).into_iter().enumerate() {
        doc = doc.add(
            Text::new(line)
                .set("x", 88)
                .set("y", dyn_h - 64.0 + line_idx as f32 * 14.0)
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#374151"),
        );
    }

    doc.to_string()
}

fn render_isoform_architecture(view: &IsoformArchitectureExpertView) -> String {
    #[derive(Debug, Clone)]
    struct PairedLaneLayout {
        compressed: CompressedTranscriptLaneLayout,
        protein: LocalProteinLaneLayout,
    }

    let lane_count = view.transcript_lanes.len().max(1);
    let lane_h = 24.0_f32;
    let lane_gap = 10.0_f32;
    let top_header = 114.0_f32;
    let exon_chart_h =
        lane_count as f32 * lane_h + (lane_count.saturating_sub(1) as f32) * lane_gap;
    let label_x = 44.0_f32;
    let left = 340.0_f32;
    let right = W - 44.0_f32;
    let width = (right - left).max(1.0);

    let region_start = view.region_start_1based.max(1);
    let region_end = view.region_end_1based.max(region_start);
    let region_span = (region_end - region_start + 1).max(1) as f32;
    let dominant_strand_is_reverse = {
        let mut plus = 0usize;
        let mut minus = 0usize;
        for lane in view.transcript_lanes.iter().filter(|lane| lane.mapped) {
            match lane.strand.trim() {
                "+" => plus += 1,
                "-" => minus += 1,
                _ => {}
            }
        }
        if plus == 0 && minus == 0 {
            for lane in &view.transcript_lanes {
                match lane.strand.trim() {
                    "+" => plus += 1,
                    "-" => minus += 1,
                    _ => {}
                }
            }
        }
        minus > plus
    };
    let left_axis_bp = if dominant_strand_is_reverse {
        region_end
    } else {
        region_start
    };
    let right_axis_bp = if dominant_strand_is_reverse {
        region_start
    } else {
        region_end
    };
    let x_for_genomic = |pos_1based: usize| -> f32 {
        let clamped = pos_1based.clamp(region_start, region_end);
        let rel = if dominant_strand_is_reverse {
            (region_end.saturating_sub(clamped)) as f32 / region_span
        } else {
            (clamped.saturating_sub(region_start)) as f32 / region_span
        };
        left + rel * width
    };
    let exon_family_registry = build_exon_family_registry(view);
    let exon_family_columns = layout_exon_family_columns(
        &exon_family_registry,
        left,
        right,
        dominant_strand_is_reverse,
    );

    let paired_lane_layouts = view
        .transcript_lanes
        .iter()
        .enumerate()
        .map(|(idx, transcript_lane)| {
            let fallback_protein_lane = gentle_protocol::IsoformArchitectureProteinLane {
                isoform_id: transcript_lane.isoform_id.clone(),
                label: transcript_lane.label.clone(),
                transcript_id: transcript_lane.transcript_id.clone(),
                expected_length_aa: None,
                reference_start_aa: None,
                reference_end_aa: None,
                domains: vec![],
                transactivation_class: transcript_lane.transactivation_class.clone(),
                comparison: None,
            };
            let protein_lane = view
                .protein_lanes
                .get(idx)
                .unwrap_or(&fallback_protein_lane);
            PairedLaneLayout {
                compressed: layout_compressed_transcript_lane(
                    transcript_lane,
                    &exon_family_columns,
                    &exon_family_registry,
                ),
                protein: layout_local_protein_lane(
                    transcript_lane,
                    protein_lane,
                    left,
                    right,
                    &exon_family_registry,
                ),
            }
        })
        .collect::<Vec<_>>();
    let longest_local_product = paired_lane_layouts
        .iter()
        .map(|layout| layout.protein.local_length_aa)
        .max()
        .unwrap_or(1)
        .max(1);
    let paired_chart_top = top_header + exon_chart_h + 106.0;
    let paired_chart_h = if paired_lane_layouts.is_empty() {
        TRANSCRIPT_PRODUCT_TOP_PAD
            + COMPRESSED_EXON_HALF_HEIGHT * 2.0
            + TRANSCRIPT_PRODUCT_CONNECTOR_GAP
            + PROTEIN_DOMAIN_HALF_HEIGHT * 2.0
            + PROTEIN_LANE_BOTTOM_PAD
    } else {
        paired_lane_layouts
            .iter()
            .map(|layout| {
                TRANSCRIPT_PRODUCT_TOP_PAD
                    + COMPRESSED_EXON_HALF_HEIGHT * 2.0
                    + TRANSCRIPT_PRODUCT_CONNECTOR_GAP
                    + layout.protein.lane_height
            })
            .sum::<f32>()
            + TRANSCRIPT_PRODUCT_LANE_GAP * paired_lane_layouts.len().saturating_sub(1) as f32
    };
    let footer_top = paired_chart_top + paired_chart_h + 64.0;
    let dyn_h = (footer_top + 120.0).max(H + 180.0);
    let top_geometry_kind = if view
        .transcript_geometry_mode
        .trim()
        .eq_ignore_ascii_case("cds")
    {
        "cds"
    } else {
        "exon"
    };
    let top_geometry_title = if top_geometry_kind == "cds" {
        "A) coordinate-true transcript architecture (faint transcript exons, colored CDS blocks)"
    } else {
        "A) coordinate-true transcript / exon architecture"
    };

    let mut doc = Document::new()
        .set("viewBox", (0, 0, W, dyn_h))
        .set("width", W)
        .set("height", dyn_h)
        .add(
            Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", W)
                .set("height", dyn_h)
                .set("fill", "#ffffff"),
        );

    doc = doc
        .add(
            Text::new(format!(
                "Isoform architecture: {} panel '{}' on {}",
                view.gene_symbol, view.panel_id, view.seq_id
            ))
            .set("x", label_x)
            .set("y", 36)
            .set("font-family", "monospace")
            .set("font-size", 18)
            .set("fill", "#111111"),
        )
        .add(
            Text::new(format!(
                "genomic span {}..{} | isoforms={} | longest local product={} aa",
                view.region_start_1based,
                view.region_end_1based,
                view.transcript_lanes.len(),
                longest_local_product
            ))
            .set("x", label_x)
            .set("y", 58)
            .set("font-family", "monospace")
            .set("font-size", 12)
            .set("fill", "#4b5563"),
        )
        .add(
            Text::new(format!(
                "display orientation: transcript 5'->3' left-to-right (dominant strand {}) | top-panel geometry: {} | lower panel: exon-chain transcript + isoform-local protein axis",
                if dominant_strand_is_reverse { "-" } else { "+" },
                top_geometry_kind
            ))
            .set("x", label_x)
            .set("y", 74)
            .set("font-family", "monospace")
            .set("font-size", 11)
            .set("fill", "#4b5563"),
        )
        .add(
            Text::new(top_geometry_title)
                .set("x", label_x)
                .set("y", 88)
                .set("font-family", "monospace")
                .set("font-size", 13)
                .set("fill", "#111827"),
        )
        .add(
            Line::new()
                .set("x1", left)
                .set("y1", top_header - 16.0)
                .set("x2", right)
                .set("y2", top_header - 16.0)
                .set("stroke", "#6b7280")
                .set("stroke-width", 1.0),
        )
        .add(
            Text::new(format!("{} bp", left_axis_bp))
                .set("x", left)
                .set("y", top_header - 20.0)
                .set("text-anchor", "start")
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#4b5563"),
        )
        .add(
            Text::new(format!("{} bp", right_axis_bp))
                .set("x", right)
                .set("y", top_header - 20.0)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#4b5563"),
        );

    for (idx, lane) in view.transcript_lanes.iter().enumerate() {
        let y = top_header + idx as f32 * (lane_h + lane_gap) + lane_h * 0.5;
        let label = lane
            .transcript_id
            .as_deref()
            .map(|tx| format!("{} ({tx})", lane.label))
            .unwrap_or_else(|| lane.label.clone());
        doc = doc.add(
            Text::new(label)
                .set("x", left - 12.0)
                .set("y", y + 3.5)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", if lane.mapped { "#111827" } else { "#b45309" }),
        );
        for intron in &lane.introns {
            let xa = x_for_genomic(intron.start_1based);
            let xb = x_for_genomic(intron.end_1based);
            let (x1, x2) = if xa <= xb { (xa, xb) } else { (xb, xa) };
            if (x2 - x1).abs() < 0.5 {
                continue;
            }
            doc = doc.add(
                Line::new()
                    .set("x1", x1)
                    .set("y1", y)
                    .set("x2", x2)
                    .set("y2", y)
                    .set("stroke", "#6b7280")
                    .set("stroke-width", 1.2)
                    .set("data-track", "coordinate-intron"),
            );
        }
        if top_geometry_kind == "cds" && !lane.transcript_exons.is_empty() {
            for exon in &lane.transcript_exons {
                let exon_key = exon_identity_key(exon);
                let family_range = exon_family_range_for_key(exon_key, &exon_family_registry);
                let color = exon_color_for_key(exon_key, &exon_family_registry);
                let xa = x_for_genomic(exon.start_1based);
                let xb = x_for_genomic(exon.end_1based);
                let (x1, x2) = if xa <= xb { (xa, xb) } else { (xb, xa) };
                if (x2 - x1).abs() < 0.5 {
                    continue;
                }
                doc = doc.add(
                    Rectangle::new()
                        .set("x", x1)
                        .set("y", y - 7.5)
                        .set("width", (x2 - x1).max(1.0))
                        .set("height", 15.0)
                        .set("fill", color.fill)
                        .set("fill-opacity", 0.25)
                        .set("stroke", color.stroke)
                        .set("stroke-width", 0.35)
                        .set("data-track", "coordinate-transcript-exon")
                        .set("data-exon-key", exon_identity_attr(exon_key))
                        .set("data-exon-family", exon_family_attr(family_range)),
                );
            }
        }
        if !lane.cds_to_protein_segments.is_empty() {
            for (segment_index, segment) in lane.cds_to_protein_segments.iter().enumerate() {
                let exon_key = lane_exon_identity_key_for_coding_range(
                    lane,
                    segment.genomic_start_1based,
                    segment.genomic_end_1based,
                );
                let family_range = exon_family_range_for_key(exon_key, &exon_family_registry);
                let color = exon_color_for_key(exon_key, &exon_family_registry);
                let xa = x_for_genomic(segment.genomic_start_1based);
                let xb = x_for_genomic(segment.genomic_end_1based);
                let (x1, x2) = if xa <= xb { (xa, xb) } else { (xb, xa) };
                if (x2 - x1).abs() < 0.5 {
                    continue;
                }
                doc = doc.add(
                    Rectangle::new()
                        .set("x", x1)
                        .set("y", y - 7.5)
                        .set("width", (x2 - x1).max(1.0))
                        .set("height", 15.0)
                        .set("fill", color.fill)
                        .set("fill-opacity", if lane.mapped { 0.88 } else { 0.45 })
                        .set("stroke", color.stroke)
                        .set("stroke-width", 0.6)
                        .set("data-track", "coordinate-cds-block")
                        .set("data-segment-rank", segment_index + 1)
                        .set("data-exon-key", exon_identity_attr(exon_key))
                        .set("data-exon-family", exon_family_attr(family_range)),
                );
            }
        } else {
            for (segment_index, exon) in lane.exons.iter().enumerate() {
                let exon_key = exon_identity_key(exon);
                let family_range = exon_family_range_for_key(exon_key, &exon_family_registry);
                let color = exon_color_for_key(exon_key, &exon_family_registry);
                let xa = x_for_genomic(exon.start_1based);
                let xb = x_for_genomic(exon.end_1based);
                let (x1, x2) = if xa <= xb { (xa, xb) } else { (xb, xa) };
                if (x2 - x1).abs() < 0.5 {
                    continue;
                }
                doc = doc.add(
                    Rectangle::new()
                        .set("x", x1)
                        .set("y", y - 7.5)
                        .set("width", (x2 - x1).max(1.0))
                        .set("height", 15.0)
                        .set("fill", if lane.mapped { color.fill } else { "#f59e0b" })
                        .set("fill-opacity", if lane.mapped { 0.85 } else { 0.45 })
                        .set("stroke", if lane.mapped { color.stroke } else { "#92400e" })
                        .set("stroke-width", 0.5)
                        .set("data-track", "coordinate-cds-block")
                        .set("data-segment-rank", segment_index + 1)
                        .set("data-exon-key", exon_identity_attr(exon_key))
                        .set("data-exon-family", exon_family_attr(family_range)),
                );
            }
        }
        if let Some(tag) = lane.transactivation_class.as_deref() {
            doc = doc.add(
                Text::new(format!("TA={tag}"))
                    .set("x", right + 6.0)
                    .set("y", y + 3.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#374151"),
            );
        }
    }

    doc = doc
        .add(
            Text::new(
                "B) shared genomic-exon columns + isoform-local protein products (same colors = one genomic exon family and its translated peptide contribution)",
            )
            .set("x", label_x)
            .set("y", paired_chart_top - 34.0)
            .set("font-family", "monospace")
            .set("font-size", 13)
            .set("fill", "#111827"),
        )
        .add(
            Text::new("Each lower protein rail is isoform-local (1..length aa), while the top panel keeps true genomic exon and intron positions.")
                .set("x", label_x)
                .set("y", paired_chart_top - 18.0)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#475569"),
        )
        .add(
            Line::new()
                .set("x1", left)
                .set("y1", paired_chart_top - 12.0)
                .set("x2", right)
                .set("y2", paired_chart_top - 12.0)
                .set("stroke", "#6b7280")
                .set("stroke-width", 1.0),
        );

    let mut paired_lane_top = paired_chart_top;
    for (idx, transcript_lane) in view.transcript_lanes.iter().enumerate() {
        let Some(layout) = paired_lane_layouts.get(idx) else {
            continue;
        };
        let fallback_protein_lane = gentle_protocol::IsoformArchitectureProteinLane {
            isoform_id: transcript_lane.isoform_id.clone(),
            label: transcript_lane.label.clone(),
            transcript_id: transcript_lane.transcript_id.clone(),
            expected_length_aa: None,
            reference_start_aa: None,
            reference_end_aa: None,
            domains: vec![],
            transactivation_class: transcript_lane.transactivation_class.clone(),
            comparison: None,
        };
        let protein_lane = view
            .protein_lanes
            .get(idx)
            .unwrap_or(&fallback_protein_lane);
        let transcript_y =
            paired_lane_top + TRANSCRIPT_PRODUCT_TOP_PAD + COMPRESSED_EXON_HALF_HEIGHT;
        let protein_y = transcript_y
            + COMPRESSED_EXON_HALF_HEIGHT
            + TRANSCRIPT_PRODUCT_CONNECTOR_GAP
            + PROTEIN_DOMAIN_HALF_HEIGHT;
        let lane_label = transcript_lane
            .transcript_id
            .as_deref()
            .map(|tx| format!("{} ({tx})", transcript_lane.label))
            .unwrap_or_else(|| transcript_lane.label.clone());
        doc = doc.add(
            Text::new(lane_label)
                .set("x", left - 12.0)
                .set("y", transcript_y + 3.5)
                .set("text-anchor", "end")
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set(
                    "fill",
                    if transcript_lane.mapped {
                        "#111827"
                    } else {
                        "#b45309"
                    },
                ),
        );
        let exon_boxes = &layout.compressed.exon_boxes;
        for pair in exon_boxes.windows(2) {
            let left_exon_right = pair[0].x2;
            let right_exon_left = pair[1].x1;
            doc = doc.add(
                Line::new()
                    .set("x1", left_exon_right)
                    .set("y1", transcript_y)
                    .set("x2", right_exon_left)
                    .set("y2", transcript_y)
                    .set("stroke", "#94a3b8")
                    .set("stroke-width", 1.1)
                    .set("stroke-dasharray", "4 3")
                    .set("data-track", "compressed-intron"),
            );
        }
        for exon_box in exon_boxes {
            doc = doc.add(
                Rectangle::new()
                    .set("x", exon_box.x1)
                    .set("y", transcript_y - COMPRESSED_EXON_HALF_HEIGHT)
                    .set("width", (exon_box.x2 - exon_box.x1).max(1.0))
                    .set("height", COMPRESSED_EXON_HALF_HEIGHT * 2.0)
                    .set("fill", exon_box.fill)
                    .set("fill-opacity", 0.26)
                    .set("stroke", exon_box.stroke)
                    .set("stroke-width", 0.5)
                    .set("data-track", "compressed-transcript-exon")
                    .set("data-exon-key", exon_identity_attr(exon_box.exon_key))
                    .set("data-exon-family", exon_family_attr(exon_box.family_range)),
            );
        }
        let mut transcript_segment_bounds: BTreeMap<usize, (f32, f32, &'static str, &'static str)> =
            BTreeMap::new();
        for block in &layout.compressed.coding_blocks {
            doc = doc.add(
                Rectangle::new()
                    .set("x", block.x1)
                    .set("y", transcript_y - COMPRESSED_EXON_HALF_HEIGHT)
                    .set("width", (block.x2 - block.x1).max(1.0))
                    .set("height", COMPRESSED_EXON_HALF_HEIGHT * 2.0)
                    .set("fill", block.fill)
                    .set(
                        "fill-opacity",
                        if transcript_lane.mapped { 0.88 } else { 0.45 },
                    )
                    .set("stroke", block.stroke)
                    .set("stroke-width", 0.6)
                    .set("data-track", "compressed-cds-block")
                    .set("data-segment-rank", block.segment_index + 1)
                    .set("data-exon-key", exon_identity_attr(block.exon_key))
                    .set("data-exon-family", exon_family_attr(block.family_range)),
            );
            transcript_segment_bounds
                .entry(block.segment_index)
                .and_modify(|entry| {
                    entry.0 = entry.0.min(block.x1);
                    entry.1 = entry.1.max(block.x2);
                })
                .or_insert((block.x1, block.x2, block.fill, block.stroke));
        }

        for rect in &layout.protein.contribution_rects {
            if let Some((tx_left, tx_right, fill, _stroke)) =
                transcript_segment_bounds.get(&rect.segment_index)
            {
                let ribbon = Data::new()
                    .move_to((*tx_left, transcript_y + COMPRESSED_EXON_HALF_HEIGHT + 0.5))
                    .line_to((*tx_right, transcript_y + COMPRESSED_EXON_HALF_HEIGHT + 0.5))
                    .line_to((rect.x2, protein_y - PROTEIN_CONTRIBUTION_HALF_HEIGHT - 0.5))
                    .line_to((rect.x1, protein_y - PROTEIN_CONTRIBUTION_HALF_HEIGHT - 0.5))
                    .close();
                doc = doc.add(
                    Path::new()
                        .set("d", ribbon)
                        .set("fill", *fill)
                        .set("fill-opacity", 0.12)
                        .set("stroke", "none")
                        .set("data-track", "transcript-product-link")
                        .set("data-segment-rank", rect.segment_index + 1)
                        .set("data-exon-key", exon_identity_attr(rect.exon_key))
                        .set("data-exon-family", exon_family_attr(rect.family_range)),
                );
            }
        }

        doc = doc.add(
            Line::new()
                .set("x1", left)
                .set("y1", protein_y)
                .set("x2", right)
                .set("y2", protein_y)
                .set("stroke", "#94a3b8")
                .set("stroke-width", 2.0)
                .set("data-track", "protein-rail"),
        );
        for rect in &layout.protein.contribution_rects {
            doc = doc.add(
                Rectangle::new()
                    .set("x", rect.x1)
                    .set("y", protein_y - PROTEIN_CONTRIBUTION_HALF_HEIGHT)
                    .set("width", (rect.x2 - rect.x1).max(1.0))
                    .set("height", PROTEIN_CONTRIBUTION_HALF_HEIGHT * 2.0)
                    .set("rx", 1.5)
                    .set("ry", 1.5)
                    .set("fill", rect.fill)
                    .set("fill-opacity", 0.32)
                    .set("stroke", rect.stroke)
                    .set("stroke-width", 0.45)
                    .set("data-track", "protein-contribution")
                    .set("data-segment-rank", rect.segment_index + 1)
                    .set("data-exon-key", exon_identity_attr(rect.exon_key))
                    .set("data-exon-family", exon_family_attr(rect.family_range)),
            );
        }

        let topology_block_top =
            protein_y + PROTEIN_DOMAIN_HALF_HEIGHT + PROTEIN_TOPOLOGY_TRACK_GAP;
        let topology_block_height = if layout.protein.topology_rows.row_count == 0 {
            0.0
        } else {
            layout.protein.topology_rows.row_count as f32 * PROTEIN_TOPOLOGY_ROW_HEIGHT
                + (layout.protein.topology_rows.row_count.saturating_sub(1) as f32)
                    * PROTEIN_TOPOLOGY_ROW_GAP
        };
        for (piece_idx, piece) in layout.protein.topology_pieces.iter().enumerate() {
            let domain = &protein_lane.domains[piece.domain_index];
            let row_index = layout
                .protein
                .topology_rows
                .row_indices
                .get(piece_idx)
                .copied()
                .unwrap_or(0);
            let rect_y = topology_block_top
                + row_index as f32 * (PROTEIN_TOPOLOGY_ROW_HEIGHT + PROTEIN_TOPOLOGY_ROW_GAP);
            let rect_width = (piece.x2 - piece.x1).max(1.0);
            doc = doc.add(
                Rectangle::new()
                    .set("x", piece.x1)
                    .set("y", rect_y)
                    .set("width", rect_width)
                    .set("height", PROTEIN_TOPOLOGY_ROW_HEIGHT)
                    .set("rx", 2.0)
                    .set("ry", 2.0)
                    .set("fill", topology_feature_fill(&domain.name))
                    .set("fill-opacity", 0.9)
                    .set("stroke", topology_feature_stroke(&domain.name))
                    .set("stroke-width", 0.6)
                    .set("data-track", "protein-topology")
                    .set(
                        "data-feature-key",
                        protein_feature_key(&domain.name).unwrap_or("unknown"),
                    ),
            );
            if let Some(label) = compact_topology_label(&domain.name, rect_width) {
                doc = doc.add(
                    Text::new(label)
                        .set("x", (piece.x1 + piece.x2) * 0.5)
                        .set(
                            "y",
                            rect_y
                                + PROTEIN_TOPOLOGY_ROW_HEIGHT * 0.5
                                + PROTEIN_TOPOLOGY_LABEL_FONT_SIZE * 0.35,
                        )
                        .set("text-anchor", "middle")
                        .set("font-family", "monospace")
                        .set("font-size", PROTEIN_TOPOLOGY_LABEL_FONT_SIZE)
                        .set("fill", "#111827"),
                );
            }
        }

        let label_band_top = if topology_block_height > 0.0 {
            topology_block_top + topology_block_height + PROTEIN_DOMAIN_LABEL_DOMAIN_GAP
        } else {
            protein_y + PROTEIN_DOMAIN_HALF_HEIGHT + PROTEIN_DOMAIN_LABEL_DOMAIN_GAP
        };
        for piece in &layout.protein.overlay_pieces {
            let domain = &protein_lane.domains[piece.domain_index];
            let fill = domain.color_hex.as_deref().unwrap_or("#7c3aed");
            doc = doc.add(
                Rectangle::new()
                    .set("x", piece.x1)
                    .set("y", protein_y - PROTEIN_DOMAIN_HALF_HEIGHT)
                    .set("width", (piece.x2 - piece.x1).max(1.0))
                    .set("height", PROTEIN_DOMAIN_HALF_HEIGHT * 2.0)
                    .set("fill", fill)
                    .set("fill-opacity", 0.82)
                    .set("stroke", "#1f2937")
                    .set("stroke-width", 0.5)
                    .set("data-track", "protein-overlay")
                    .set(
                        "data-feature-key",
                        protein_feature_key(&domain.name).unwrap_or("unknown"),
                    ),
            );
        }
        for (label_idx, label) in layout.protein.overlay_labels.placements.iter().enumerate() {
            let Some(&domain_index) = layout.protein.overlay_label_domain_indices.get(label_idx)
            else {
                continue;
            };
            let anchor_domain = layout
                .protein
                .overlay_pieces
                .iter()
                .filter(|piece| piece.domain_index == domain_index)
                .fold(None::<(f32, f32)>, |acc, piece| match acc {
                    Some((left_edge, right_edge)) => {
                        Some((left_edge.min(piece.x1), right_edge.max(piece.x2)))
                    }
                    None => Some((piece.x1, piece.x2)),
                });
            let Some((anchor_left, anchor_right)) = anchor_domain else {
                continue;
            };
            let label_baseline_y = label_band_top
                + PROTEIN_DOMAIN_LABEL_FONT_SIZE
                + label.row_index as f32 * PROTEIN_DOMAIN_LABEL_ROW_PITCH;
            let label_box_top = label_baseline_y - PROTEIN_DOMAIN_LABEL_FONT_SIZE;
            let label_box_height =
                PROTEIN_DOMAIN_LABEL_FONT_SIZE + PROTEIN_DOMAIN_LABEL_BOX_Y_PAD * 2.0;
            let domain_center = (anchor_left + anchor_right) * 0.5;
            let leader_top = protein_y + PROTEIN_DOMAIN_HALF_HEIGHT + 1.0;
            let leader_bottom =
                (label_box_top - PROTEIN_DOMAIN_LABEL_BOX_Y_PAD - 1.0).max(leader_top);
            doc = doc
                .add(
                    Line::new()
                        .set("x1", domain_center)
                        .set("y1", leader_top)
                        .set("x2", domain_center)
                        .set("y2", leader_bottom)
                        .set("stroke", "#94a3b8")
                        .set("stroke-width", 0.8)
                        .set("stroke-opacity", 0.9),
                )
                .add(
                    Rectangle::new()
                        .set("x", label.box_left)
                        .set("y", label_box_top - PROTEIN_DOMAIN_LABEL_BOX_Y_PAD)
                        .set("width", (label.box_right - label.box_left).max(1.0))
                        .set("height", label_box_height)
                        .set("rx", 2.0)
                        .set("ry", 2.0)
                        .set("fill", "#ffffff")
                        .set("fill-opacity", 0.96)
                        .set("stroke", "#cbd5e1")
                        .set("stroke-width", 0.5),
                )
                .add(
                    Text::new(label.compact_label.clone())
                        .set("x", label.box_left + PROTEIN_DOMAIN_LABEL_BOX_X_PAD)
                        .set("y", label_baseline_y + PROTEIN_DOMAIN_LABEL_BOX_Y_PAD * 0.5)
                        .set("font-family", "monospace")
                        .set("font-size", PROTEIN_DOMAIN_LABEL_FONT_SIZE)
                        .set("fill", "#334155"),
                );
        }

        let protein_note = match layout.protein.reference_span_label.as_deref() {
            Some(reference_span) => {
                format!("{} aa | {}", layout.protein.local_length_aa, reference_span)
            }
            None => format!("{} aa", layout.protein.local_length_aa),
        };
        doc = doc.add(
            Text::new(protein_note)
                .set("x", right + 6.0)
                .set("y", protein_y + 3.0)
                .set("font-family", "monospace")
                .set("font-size", 9)
                .set("fill", "#374151"),
        );
        if let Some(comparison) = protein_lane.comparison.as_ref() {
            let status_label = match comparison.status {
                gentle_protocol::TranscriptProteinComparisonStatus::DerivedOnly => None,
                gentle_protocol::TranscriptProteinComparisonStatus::ConsistentWithExternalOpinion => {
                    None
                }
                gentle_protocol::TranscriptProteinComparisonStatus::LowConfidenceExternalOpinion => {
                    Some(("status=low_confidence_external_opinion", "#b45309"))
                }
                gentle_protocol::TranscriptProteinComparisonStatus::ExternalOpinionOnly => {
                    Some(("status=external_opinion_only", "#64748b"))
                }
                gentle_protocol::TranscriptProteinComparisonStatus::NoTranscriptCds => {
                    Some(("status=no_transcript_cds", "#64748b"))
                }
            };
            if let Some((status_label, color)) = status_label {
                doc = doc.add(
                    Text::new(status_label)
                        .set("x", right + 6.0)
                        .set("y", protein_y + 15.0)
                        .set("font-family", "monospace")
                        .set("font-size", 8.5)
                        .set("fill", color),
                );
            }
        }
        if let Some(tag) = protein_lane.transactivation_class.as_deref() {
            doc = doc.add(
                Text::new(format!("TA={tag}"))
                    .set("x", right + 6.0)
                    .set("y", transcript_y + 3.0)
                    .set("font-family", "monospace")
                    .set("font-size", 9)
                    .set("fill", "#374151"),
            );
        }

        paired_lane_top += TRANSCRIPT_PRODUCT_TOP_PAD
            + COMPRESSED_EXON_HALF_HEIGHT * 2.0
            + TRANSCRIPT_PRODUCT_CONNECTOR_GAP
            + layout.protein.lane_height
            + TRANSCRIPT_PRODUCT_LANE_GAP;
    }

    let mut y = footer_top;
    if let Some(source) = view
        .panel_source
        .as_deref()
        .filter(|value| !value.trim().is_empty())
    {
        doc = doc.add(
            Text::new(format!("panel source: {source}"))
                .set("x", label_x)
                .set("y", y)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#334155"),
        );
        y += 14.0;
    }
    for warning in view.warnings.iter().take(6) {
        doc = doc.add(
            Text::new(format!("warning: {warning}"))
                .set("x", label_x)
                .set("y", y)
                .set("font-family", "monospace")
                .set("font-size", 10)
                .set("fill", "#b45309"),
        );
        y += 12.0;
    }
    for (line_idx, line) in wrap_text(&view.instruction, 140).into_iter().enumerate() {
        doc = doc.add(
            Text::new(line)
                .set("x", label_x)
                .set("y", dyn_h - 64.0 + line_idx as f32 * 14.0)
                .set("font-family", "monospace")
                .set("font-size", 11)
                .set("fill", "#374151"),
        );
    }

    doc.to_string()
}

pub fn render_feature_expert_svg(view: &FeatureExpertView) -> String {
    match view {
        FeatureExpertView::Tfbs(tfbs) => render_tfbs(tfbs),
        FeatureExpertView::RestrictionSite(re) => render_restriction(re),
        FeatureExpertView::Splicing(splicing) => render_splicing(splicing),
        FeatureExpertView::IsoformArchitecture(isoform) => render_isoform_architecture(isoform),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle_protocol::{
        IsoformArchitectureCdsAaSegment, IsoformArchitectureProteinDomain,
        IsoformArchitectureProteinLane, IsoformArchitectureTranscriptLane, SplicingRange,
    };
    use std::collections::{BTreeMap, BTreeSet};

    fn isoform_test_view(transcript_strand: &str) -> IsoformArchitectureExpertView {
        IsoformArchitectureExpertView {
            seq_id: "tp53".to_string(),
            panel_id: "panel".to_string(),
            gene_symbol: "TP53".to_string(),
            transcript_geometry_mode: "cds".to_string(),
            panel_source: Some("test".to_string()),
            region_start_1based: 101,
            region_end_1based: 220,
            instruction: "test".to_string(),
            transcript_lanes: vec![IsoformArchitectureTranscriptLane {
                isoform_id: "i1".to_string(),
                label: "iso1".to_string(),
                transcript_id: Some("tx1".to_string()),
                transcript_feature_id: Some(1),
                strand: transcript_strand.to_string(),
                transcript_exons: vec![
                    SplicingRange {
                        start_1based: 105,
                        end_1based: 138,
                    },
                    SplicingRange {
                        start_1based: 148,
                        end_1based: 186,
                    },
                ],
                exons: vec![
                    SplicingRange {
                        start_1based: 110,
                        end_1based: 130,
                    },
                    SplicingRange {
                        start_1based: 150,
                        end_1based: 176,
                    },
                ],
                introns: vec![SplicingRange {
                    start_1based: 131,
                    end_1based: 149,
                }],
                mapped: true,
                transactivation_class: None,
                cds_to_protein_segments: vec![
                    IsoformArchitectureCdsAaSegment {
                        genomic_start_1based: 110,
                        genomic_end_1based: 130,
                        aa_start: 1,
                        aa_end: 21,
                    },
                    IsoformArchitectureCdsAaSegment {
                        genomic_start_1based: 150,
                        genomic_end_1based: 176,
                        aa_start: 22,
                        aa_end: 48,
                    },
                ],
                note: None,
            }],
            protein_lanes: vec![IsoformArchitectureProteinLane {
                isoform_id: "i1".to_string(),
                label: "iso1".to_string(),
                transcript_id: Some("tx1".to_string()),
                expected_length_aa: Some(48),
                reference_start_aa: Some(1),
                reference_end_aa: Some(48),
                domains: vec![IsoformArchitectureProteinDomain {
                    name: "dbd".to_string(),
                    start_aa: 8,
                    end_aa: 38,
                    color_hex: Some("#ff0000".to_string()),
                }],
                transactivation_class: None,
                comparison: None,
            }],
            warnings: vec![],
        }
    }

    fn extract_track_x_positions(svg: &str, track: &str) -> Vec<f32> {
        let mut xs = Vec::new();
        for part in svg.split("<rect").skip(1) {
            if !part.contains(track) {
                continue;
            }
            let Some(attr_start) = part.find("x=\"") else {
                continue;
            };
            let rest = &part[attr_start + 3..];
            let Some(attr_end) = rest.find('"') else {
                continue;
            };
            if let Ok(value) = rest[..attr_end].parse::<f32>() {
                xs.push(value);
            }
        }
        xs
    }

    fn track_count(svg: &str, track: &str) -> usize {
        svg.matches(track).count()
    }

    fn extract_fill_for_track_and_exon_key(
        svg: &str,
        element: &str,
        track: &str,
        exon_key: &str,
    ) -> Option<String> {
        for part in svg.split(element).skip(1) {
            if !part.contains(track) || !part.contains(exon_key) {
                continue;
            }
            let Some(attr_start) = part.find("fill=\"") else {
                continue;
            };
            let rest = &part[attr_start + 6..];
            let Some(attr_end) = rest.find('"') else {
                continue;
            };
            return Some(rest[..attr_end].to_string());
        }
        None
    }

    fn extract_fill_for_track_and_exon_family(
        svg: &str,
        element: &str,
        track: &str,
        exon_family: &str,
    ) -> Option<String> {
        for part in svg.split(element).skip(1) {
            if !part.contains(track) || !part.contains(exon_family) {
                continue;
            }
            let Some(attr_start) = part.find("fill=\"") else {
                continue;
            };
            let rest = &part[attr_start + 6..];
            let Some(attr_end) = rest.find('"') else {
                continue;
            };
            return Some(rest[..attr_end].to_string());
        }
        None
    }

    fn collect_track_xs_by_exon_family(
        svg: &str,
        track: &str,
    ) -> BTreeMap<String, BTreeSet<String>> {
        let mut out: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        for part in svg.split("<rect").skip(1) {
            if !part.contains(track) {
                continue;
            }
            let Some(x_attr_start) = part.find("x=\"") else {
                continue;
            };
            let x_rest = &part[x_attr_start + 3..];
            let Some(x_attr_end) = x_rest.find('"') else {
                continue;
            };
            let Ok(value) = x_rest[..x_attr_end].parse::<f32>() else {
                continue;
            };
            let Some(family_attr_start) = part.find("data-exon-family=\"") else {
                continue;
            };
            let family_rest = &part[family_attr_start + 18..];
            let Some(family_attr_end) = family_rest.find('"') else {
                continue;
            };
            out.entry(format!("{value:.2}"))
                .or_default()
                .insert(family_rest[..family_attr_end].to_string());
        }
        out
    }

    fn extract_first_protein_rail_y(svg: &str) -> Option<f32> {
        for part in svg.split("<line").skip(1) {
            if !part.contains("stroke=\"#94a3b8\"") || !part.contains("stroke-width=\"2\"") {
                continue;
            }
            let Some(attr_start) = part.find(" y1=\"") else {
                continue;
            };
            let rest = &part[attr_start + 5..];
            let Some(attr_end) = rest.find('"') else {
                continue;
            };
            if let Ok(value) = rest[..attr_end].parse::<f32>() {
                return Some(value);
            }
        }
        None
    }

    fn extract_first_label_box_y(svg: &str) -> Option<f32> {
        for part in svg.split("<rect").skip(1) {
            if !part.contains("fill=\"#ffffff\"") || !part.contains("stroke=\"#cbd5e1\"") {
                continue;
            }
            let Some(attr_start) = part.find(" y=\"") else {
                continue;
            };
            let rest = &part[attr_start + 4..];
            let Some(attr_end) = rest.find('"') else {
                continue;
            };
            if let Ok(value) = rest[..attr_end].parse::<f32>() {
                return Some(value);
            }
        }
        None
    }

    fn extract_first_topology_rect_y(svg: &str) -> Option<f32> {
        for part in svg.split("<rect").skip(1) {
            if !part.contains("data-track=\"protein-topology\"") {
                continue;
            }
            let Some(attr_start) = part.find(" y=\"") else {
                continue;
            };
            let rest = &part[attr_start + 4..];
            let Some(attr_end) = rest.find('"') else {
                continue;
            };
            if let Ok(value) = rest[..attr_end].parse::<f32>() {
                return Some(value);
            }
        }
        None
    }

    #[test]
    fn isoform_renderer_keeps_forward_strand_exon_order_left_to_right() {
        let svg = render_isoform_architecture(&isoform_test_view("+"));
        let exon_x = extract_track_x_positions(&svg, "data-track=\"coordinate-cds-block\"");
        assert_eq!(exon_x.len(), 2);
        assert!(exon_x[0] < exon_x[1]);
        assert!(svg.contains("dominant strand +"));
        assert!(svg.contains("101 bp"));
        assert!(svg.contains("220 bp"));
    }

    #[test]
    fn isoform_renderer_flips_reverse_strand_exons_but_keeps_them_visible() {
        let svg = render_isoform_architecture(&isoform_test_view("-"));
        let exon_x = extract_track_x_positions(&svg, "data-track=\"coordinate-cds-block\"");
        assert_eq!(exon_x.len(), 2);
        assert!(exon_x[0] > exon_x[1]);
        assert!(svg.contains("dominant strand -"));
        assert!(svg.contains("220 bp"));
        assert!(svg.contains("101 bp"));
    }

    #[test]
    fn isoform_renderer_draws_cds_to_protein_connector_guides() {
        let view = isoform_test_view("+");
        let svg = render_isoform_architecture(&view);
        assert_eq!(track_count(&svg, "data-track=\"compressed-cds-block\""), 2);
        assert_eq!(track_count(&svg, "data-track=\"protein-contribution\""), 2);
        assert_eq!(
            track_count(&svg, "data-track=\"transcript-product-link\""),
            2
        );
        assert!(svg.contains("isoform-local protein axis"));
    }

    #[test]
    fn isoform_renderer_keeps_same_genomic_exon_color_across_isoforms_and_panels() {
        let mut view = isoform_test_view("+");
        view.transcript_lanes
            .push(IsoformArchitectureTranscriptLane {
                isoform_id: "i2".to_string(),
                label: "iso2".to_string(),
                transcript_id: Some("tx2".to_string()),
                transcript_feature_id: Some(2),
                strand: "+".to_string(),
                transcript_exons: vec![
                    SplicingRange {
                        start_1based: 148,
                        end_1based: 186,
                    },
                    SplicingRange {
                        start_1based: 196,
                        end_1based: 214,
                    },
                ],
                exons: vec![
                    SplicingRange {
                        start_1based: 150,
                        end_1based: 176,
                    },
                    SplicingRange {
                        start_1based: 198,
                        end_1based: 210,
                    },
                ],
                introns: vec![SplicingRange {
                    start_1based: 177,
                    end_1based: 197,
                }],
                mapped: true,
                transactivation_class: None,
                cds_to_protein_segments: vec![
                    IsoformArchitectureCdsAaSegment {
                        genomic_start_1based: 150,
                        genomic_end_1based: 176,
                        aa_start: 1,
                        aa_end: 27,
                    },
                    IsoformArchitectureCdsAaSegment {
                        genomic_start_1based: 198,
                        genomic_end_1based: 210,
                        aa_start: 28,
                        aa_end: 40,
                    },
                ],
                note: None,
            });
        view.protein_lanes.push(IsoformArchitectureProteinLane {
            isoform_id: "i2".to_string(),
            label: "iso2".to_string(),
            transcript_id: Some("tx2".to_string()),
            expected_length_aa: Some(40),
            reference_start_aa: Some(1),
            reference_end_aa: Some(40),
            domains: vec![IsoformArchitectureProteinDomain {
                name: "dbd-short".to_string(),
                start_aa: 4,
                end_aa: 24,
                color_hex: Some("#00aa00".to_string()),
            }],
            transactivation_class: None,
            comparison: None,
        });

        let svg = render_isoform_architecture(&view);
        let exon_key = "data-exon-key=\"148..186\"";
        let coordinate_fill = extract_fill_for_track_and_exon_key(
            &svg,
            "<rect",
            "data-track=\"coordinate-cds-block\"",
            exon_key,
        )
        .expect("coordinate fill");
        let compressed_fill = extract_fill_for_track_and_exon_key(
            &svg,
            "<rect",
            "data-track=\"compressed-cds-block\"",
            exon_key,
        )
        .expect("compressed fill");
        let protein_fill = extract_fill_for_track_and_exon_key(
            &svg,
            "<rect",
            "data-track=\"protein-contribution\"",
            exon_key,
        )
        .expect("protein contribution fill");
        let link_fill = extract_fill_for_track_and_exon_key(
            &svg,
            "<path",
            "data-track=\"transcript-product-link\"",
            exon_key,
        )
        .expect("link fill");
        assert_eq!(coordinate_fill, compressed_fill);
        assert_eq!(compressed_fill, protein_fill);
        assert_eq!(protein_fill, link_fill);
    }

    #[test]
    fn isoform_renderer_uses_shared_genomic_exon_columns_in_lower_panel() {
        let mut view = isoform_test_view("+");
        view.transcript_lanes
            .push(IsoformArchitectureTranscriptLane {
                isoform_id: "i2".to_string(),
                label: "iso2".to_string(),
                transcript_id: Some("tx2".to_string()),
                transcript_feature_id: Some(2),
                strand: "+".to_string(),
                transcript_exons: vec![
                    SplicingRange {
                        start_1based: 148,
                        end_1based: 186,
                    },
                    SplicingRange {
                        start_1based: 196,
                        end_1based: 214,
                    },
                ],
                exons: vec![
                    SplicingRange {
                        start_1based: 150,
                        end_1based: 176,
                    },
                    SplicingRange {
                        start_1based: 198,
                        end_1based: 210,
                    },
                ],
                introns: vec![SplicingRange {
                    start_1based: 177,
                    end_1based: 197,
                }],
                mapped: true,
                transactivation_class: None,
                cds_to_protein_segments: vec![
                    IsoformArchitectureCdsAaSegment {
                        genomic_start_1based: 150,
                        genomic_end_1based: 176,
                        aa_start: 1,
                        aa_end: 27,
                    },
                    IsoformArchitectureCdsAaSegment {
                        genomic_start_1based: 198,
                        genomic_end_1based: 210,
                        aa_start: 28,
                        aa_end: 40,
                    },
                ],
                note: None,
            });
        view.protein_lanes.push(IsoformArchitectureProteinLane {
            isoform_id: "i2".to_string(),
            label: "iso2".to_string(),
            transcript_id: Some("tx2".to_string()),
            expected_length_aa: Some(40),
            reference_start_aa: Some(1),
            reference_end_aa: Some(40),
            domains: vec![IsoformArchitectureProteinDomain {
                name: "dbd-short".to_string(),
                start_aa: 4,
                end_aa: 24,
                color_hex: Some("#00aa00".to_string()),
            }],
            transactivation_class: None,
            comparison: None,
        });

        let svg = render_isoform_architecture(&view);
        let column_families =
            collect_track_xs_by_exon_family(&svg, "data-track=\"compressed-transcript-exon\"");
        assert!(
            column_families.values().all(|families| families.len() == 1),
            "expected each compressed-transcript column x to map to one genomic exon family: {column_families:?}"
        );

        let shared_family = "data-exon-family=\"148..186\"";
        let compressed_fill = extract_fill_for_track_and_exon_family(
            &svg,
            "<rect",
            "data-track=\"compressed-transcript-exon\"",
            shared_family,
        )
        .expect("compressed transcript family fill");
        let protein_fill = extract_fill_for_track_and_exon_family(
            &svg,
            "<rect",
            "data-track=\"protein-contribution\"",
            shared_family,
        )
        .expect("protein family fill");
        assert_eq!(compressed_fill, protein_fill);
    }

    #[test]
    fn isoform_renderer_splits_domains_across_skipped_reference_gaps_on_local_axis() {
        let mut view = isoform_test_view("+");
        view.transcript_lanes[0].cds_to_protein_segments = vec![
            IsoformArchitectureCdsAaSegment {
                genomic_start_1based: 110,
                genomic_end_1based: 130,
                aa_start: 1,
                aa_end: 21,
            },
            IsoformArchitectureCdsAaSegment {
                genomic_start_1based: 150,
                genomic_end_1based: 176,
                aa_start: 30,
                aa_end: 56,
            },
        ];
        view.protein_lanes[0].expected_length_aa = Some(48);
        view.protein_lanes[0].reference_start_aa = Some(1);
        view.protein_lanes[0].reference_end_aa = Some(56);
        view.protein_lanes[0].domains = vec![IsoformArchitectureProteinDomain {
            name: "DOMAIN: split by skipped exon".to_string(),
            start_aa: 18,
            end_aa: 34,
            color_hex: Some("#7c3aed".to_string()),
        }];

        let svg = render_isoform_architecture(&view);
        assert_eq!(track_count(&svg, "data-track=\"protein-overlay\""), 2);
        assert_eq!(track_count(&svg, "data-track=\"protein-contribution\""), 2);
    }

    #[test]
    fn isoform_renderer_places_protein_domain_labels_below_domain_lane() {
        let svg = render_isoform_architecture(&isoform_test_view("+"));
        let protein_y = extract_first_protein_rail_y(&svg).expect("protein rail y");
        let label_box_y = extract_first_label_box_y(&svg).expect("protein label box y");
        assert!(
            label_box_y > protein_y,
            "expected first label box below protein rail: label_box_y={label_box_y}, protein_y={protein_y}"
        );
        assert!(svg.contains("stroke=\"#cbd5e1\""));
    }

    #[test]
    fn isoform_renderer_staggers_dense_protein_domain_labels_across_rows() {
        let dense_layout = layout_protein_domain_labels(
            &[
                (
                    420.0,
                    450.0,
                    "label alpha domain segment with long annotation text".to_string(),
                ),
                (
                    420.0,
                    450.0,
                    "label beta overlapping segment with long annotation text".to_string(),
                ),
                (
                    420.0,
                    450.0,
                    "label gamma dense segment with long annotation text".to_string(),
                ),
                (
                    420.0,
                    450.0,
                    "label delta dense segment with long annotation text".to_string(),
                ),
            ],
            340.0,
            1156.0,
        );
        assert!(
            dense_layout
                .placements
                .iter()
                .map(|placement| placement.row_index)
                .collect::<std::collections::BTreeSet<_>>()
                .len()
                >= 2
        );

        let mut view = isoform_test_view("+");
        view.protein_lanes[0].domains = vec![
            IsoformArchitectureProteinDomain {
                name: "label alpha domain segment with long annotation text".to_string(),
                start_aa: 10,
                end_aa: 18,
                color_hex: Some("#ff0000".to_string()),
            },
            IsoformArchitectureProteinDomain {
                name: "label beta overlapping segment with long annotation text".to_string(),
                start_aa: 10,
                end_aa: 18,
                color_hex: Some("#00aa00".to_string()),
            },
            IsoformArchitectureProteinDomain {
                name: "label gamma dense segment with long annotation text".to_string(),
                start_aa: 10,
                end_aa: 18,
                color_hex: Some("#0000ff".to_string()),
            },
            IsoformArchitectureProteinDomain {
                name: "label delta dense segment with long annotation text".to_string(),
                start_aa: 10,
                end_aa: 18,
                color_hex: Some("#ffaa00".to_string()),
            },
        ];

        let svg = render_isoform_architecture(&view);
        assert!(svg.contains("stroke=\"#cbd5e1\""));
    }

    #[test]
    fn isoform_renderer_places_topology_features_in_dedicated_lower_band() {
        let mut view = isoform_test_view("+");
        view.protein_lanes[0].domains = vec![
            IsoformArchitectureProteinDomain {
                name: "SIGNAL: signal peptide".to_string(),
                start_aa: 1,
                end_aa: 14,
                color_hex: Some("#22c55e".to_string()),
            },
            IsoformArchitectureProteinDomain {
                name: "TOPO_DOM: luminal".to_string(),
                start_aa: 8,
                end_aa: 24,
                color_hex: Some("#60a5fa".to_string()),
            },
            IsoformArchitectureProteinDomain {
                name: "TRANSMEM: helix".to_string(),
                start_aa: 25,
                end_aa: 32,
                color_hex: Some("#ef4444".to_string()),
            },
            IsoformArchitectureProteinDomain {
                name: "DOMAIN: catalytic tail".to_string(),
                start_aa: 34,
                end_aa: 44,
                color_hex: Some("#7c3aed".to_string()),
            },
        ];

        let svg = render_isoform_architecture(&view);
        let protein_y = extract_first_protein_rail_y(&svg).expect("protein rail y");
        let topology_rect_y = extract_first_topology_rect_y(&svg).expect("topology rect y");
        let label_box_y = extract_first_label_box_y(&svg).expect("protein label box y");
        assert!(
            topology_rect_y > protein_y,
            "expected topology band below protein rail: topology_rect_y={topology_rect_y}, protein_y={protein_y}"
        );
        assert!(
            label_box_y > topology_rect_y,
            "expected overlay/domain labels below topology band: label_box_y={label_box_y}, topology_rect_y={topology_rect_y}"
        );
        assert!(svg.contains("data-track=\"protein-topology\""));
        assert!(svg.contains("data-feature-key=\"SIGNAL\""));
        assert!(svg.contains("data-feature-key=\"TOPO_DOM\""));
        assert!(svg.contains("data-feature-key=\"TRANSMEM\""));
    }
}
