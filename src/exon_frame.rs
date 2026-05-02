//! Shared exon length and coding-frame cue helpers.
//!
//! These helpers keep display-time exon frame cues and exon-skip planning on the
//! same vocabulary without tying headless engine code to GUI rendering.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExonLengthFrameCue {
    /// Exon span length in base pairs.
    pub length_bp: usize,
    /// `length_bp % 3`, useful as a coarse whole-exon skip frame cue.
    pub length_mod3: usize,
    /// Whether omitting the full exon preserves coding length modulo 3.
    pub frame_neutral_length: bool,
}

/// Coding phase observed at the transcript-facing exon boundaries.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ExonCdsPhaseCue {
    /// CDS phase at the exon start in genomic left-to-right coordinates.
    pub left_cds_phase: Option<u8>,
    /// CDS phase at the exon end in genomic left-to-right coordinates.
    pub right_cds_phase: Option<u8>,
}

/// Frame consequence from skipping only the coding portion of an exon.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExonCodingFrameCue {
    /// CDS-overlapping base pairs removed by skipping the exon.
    pub coding_skip_bp: usize,
    /// `coding_skip_bp % 3`.
    pub coding_skip_mod3: usize,
    /// Whether the coding portion removed by the skip preserves frame length.
    pub frame_neutral_coding_skip: bool,
    /// Coarse UTR/CDS composition bucket for the exon.
    pub coding_context: &'static str,
}

impl ExonLengthFrameCue {
    /// Build a frame cue from a known exon length.
    pub fn from_length(length_bp: usize) -> Self {
        let length_mod3 = length_bp % 3;
        Self {
            length_bp,
            length_mod3,
            frame_neutral_length: length_mod3 == 0,
        }
    }

    /// Build a frame cue from a 0-based half-open genomic range.
    pub fn from_range(start_0based: usize, end_0based_exclusive: usize) -> Self {
        Self::from_length(end_0based_exclusive.saturating_sub(start_0based))
    }

    /// Human-facing frame interpretation for whole-exon skipping.
    pub fn skip_frame_hint(self, coding_context: bool) -> &'static str {
        match (self.frame_neutral_length, coding_context) {
            (true, true) => {
                "length is a multiple of 3; whole-exon omission is length-frame-neutral"
            }
            (true, false) => {
                "length is a multiple of 3; whole-exon omission would be length-frame-neutral if coding"
            }
            (false, true) => "whole-exon omission changes coding-frame length",
            (false, false) => {
                "whole-exon omission would change coding-frame length if the exon is coding"
            }
        }
    }

    /// Human-facing note when full-exon and CDS-only modulo cues differ.
    pub fn coding_frame_note(self, coding_cue: ExonCodingFrameCue) -> Option<String> {
        (coding_cue.coding_skip_bp > 0 && self.length_mod3 != coding_cue.coding_skip_mod3).then(
            || {
                format!(
                    "Whole exon length modulo 3 is {}, but the CDS-overlap skip length modulo 3 is {}; use the coding value for frame consequence.",
                    self.length_mod3, coding_cue.coding_skip_mod3
                )
            },
        )
    }
}

impl ExonCodingFrameCue {
    /// Build coding-frame skip cues from a 0-based half-open exon range and CDS ranges.
    pub fn from_exon_and_cds(
        exon_range_0based: (usize, usize),
        cds_ranges_0based: &[(usize, usize)],
    ) -> Self {
        let exon_len = exon_range_0based.1.saturating_sub(exon_range_0based.0);
        let coding_skip_bp = coding_overlap_length_0based(exon_range_0based, cds_ranges_0based);
        let coding_skip_mod3 = coding_skip_bp % 3;
        Self {
            coding_skip_bp,
            coding_skip_mod3,
            frame_neutral_coding_skip: coding_skip_mod3 == 0,
            coding_context: coding_context_for_lengths(exon_len, coding_skip_bp),
        }
    }

    /// Warning emitted when a CDS-overlapping skip changes coding length modulo 3.
    pub fn cds_phase_warning(self) -> Option<String> {
        (self.coding_skip_bp > 0 && !self.frame_neutral_coding_skip).then(|| {
            format!(
                "Skipping this CDS-overlapping exon removes {} coding bp (not divisible by 3).",
                self.coding_skip_bp
            )
        })
    }

    /// Human-facing interpretation for coding portion of whole-exon skipping.
    pub fn coding_skip_hint(self) -> &'static str {
        match (self.coding_skip_bp, self.frame_neutral_coding_skip) {
            (0, _) => "exon has no CDS-overlapping bases in this transcript",
            (_, true) => "CDS-overlap skip length is frame-neutral",
            (_, false) => "CDS-overlap skip length changes the coding frame",
        }
    }
}

/// Return the intron length between two 0-based half-open exon ranges.
pub fn intron_length_between_exons_0based(
    left_exon_0based: (usize, usize),
    right_exon_0based: (usize, usize),
) -> usize {
    let left_end = left_exon_0based.1.min(right_exon_0based.1);
    let right_start = left_exon_0based.0.max(right_exon_0based.0);
    right_start.saturating_sub(left_end)
}

/// Derive exon boundary CDS phases from ordered exon ranges and CDS ranges.
///
/// Ranges are 0-based half-open genomic intervals. For reverse-strand
/// transcripts, CDS consumption proceeds from high coordinates to low
/// coordinates, while returned phases remain attached to genomic left/right
/// exon boundaries.
pub fn exon_cds_phase_cues(
    exon_ranges_0based: &[(usize, usize)],
    cds_ranges_0based: &[(usize, usize)],
    is_reverse: bool,
) -> Vec<ExonCdsPhaseCue> {
    let mut phases = vec![ExonCdsPhaseCue::default(); exon_ranges_0based.len()];
    if exon_ranges_0based.is_empty() || cds_ranges_0based.is_empty() {
        return phases;
    }

    let mut cds_ranges = cds_ranges_0based.to_vec();
    cds_ranges.sort_unstable_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
    cds_ranges.dedup();

    let mut consumed_cds_bp = 0usize;
    let exon_indices = if is_reverse {
        (0..exon_ranges_0based.len()).rev().collect::<Vec<_>>()
    } else {
        (0..exon_ranges_0based.len()).collect::<Vec<_>>()
    };

    for exon_idx in exon_indices {
        let exon = exon_ranges_0based[exon_idx];
        let mut coding_segments = cds_ranges
            .iter()
            .filter_map(|cds| range_intersection_0based(exon, *cds))
            .collect::<Vec<_>>();
        if coding_segments.is_empty() {
            continue;
        }
        if is_reverse {
            coding_segments
                .sort_unstable_by(|left, right| right.0.cmp(&left.0).then(right.1.cmp(&left.1)));
        } else {
            coding_segments
                .sort_unstable_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
        }

        for (seg_start, seg_end) in coding_segments {
            let seg_len = seg_end.saturating_sub(seg_start);
            if seg_len == 0 {
                continue;
            }
            let entry_phase = (consumed_cds_bp % 3) as u8;
            let exit_phase = ((consumed_cds_bp + seg_len - 1) % 3) as u8;
            if is_reverse {
                if seg_end == exon.1 && phases[exon_idx].right_cds_phase.is_none() {
                    phases[exon_idx].right_cds_phase = Some(entry_phase);
                }
                if seg_start == exon.0 {
                    phases[exon_idx].left_cds_phase = Some(exit_phase);
                }
            } else {
                if seg_start == exon.0 && phases[exon_idx].left_cds_phase.is_none() {
                    phases[exon_idx].left_cds_phase = Some(entry_phase);
                }
                if seg_end == exon.1 {
                    phases[exon_idx].right_cds_phase = Some(exit_phase);
                }
            }
            consumed_cds_bp += seg_len;
        }
    }

    phases
}

/// Pick the transcript-entry CDS phase from genomic left/right boundary phases.
pub fn transcript_entry_phase(
    left_cds_phase: Option<u8>,
    right_cds_phase: Option<u8>,
    is_reverse: bool,
) -> Option<u8> {
    if is_reverse {
        right_cds_phase
    } else {
        left_cds_phase
    }
}

/// Stable machine-readable bucket for a transcript-entry CDS phase.
pub fn phase_entry_kind(entry_cds_phase: Option<u8>) -> &'static str {
    match entry_cds_phase.map(|phase| phase % 3) {
        Some(0) => "codon_boundary",
        Some(1) => "split_codon_1",
        Some(2) => "split_codon_2",
        None => "unavailable",
        _ => "unavailable",
    }
}

/// Short GUI/CLI hint for the coding entry phase.
pub fn phase_entry_hint(entry_cds_phase: Option<u8>) -> &'static str {
    match entry_cds_phase.map(|phase| phase % 3) {
        Some(0) => "starts at a codon boundary/new amino acid",
        Some(1) => "starts one nucleotide into a split codon",
        Some(2) => "starts two nucleotides into a split codon",
        None => "coding entry phase unavailable",
        _ => "coding entry phase unavailable",
    }
}

fn range_intersection_0based(
    left: (usize, usize),
    right: (usize, usize),
) -> Option<(usize, usize)> {
    let start = left.0.max(right.0);
    let end = left.1.min(right.1);
    (end > start).then_some((start, end))
}

fn coding_overlap_length_0based(
    exon_range_0based: (usize, usize),
    cds_ranges_0based: &[(usize, usize)],
) -> usize {
    let mut intersections = cds_ranges_0based
        .iter()
        .filter_map(|cds| range_intersection_0based(exon_range_0based, *cds))
        .collect::<Vec<_>>();
    if intersections.is_empty() {
        return 0;
    }
    intersections.sort_unstable_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
    let mut total = 0usize;
    let mut current = intersections[0];
    for next in intersections.into_iter().skip(1) {
        if next.0 <= current.1 {
            current.1 = current.1.max(next.1);
        } else {
            total += current.1.saturating_sub(current.0);
            current = next;
        }
    }
    total + current.1.saturating_sub(current.0)
}

fn coding_context_for_lengths(exon_len: usize, coding_skip_bp: usize) -> &'static str {
    match (coding_skip_bp, exon_len) {
        (0, _) => "utr_only",
        (coding, exon) if coding >= exon => "cds_only",
        _ => "mixed_utr_cds",
    }
}
