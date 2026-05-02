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

    /// Warning emitted when a CDS-overlapping whole-exon skip changes length modulo 3.
    pub fn cds_phase_warning(self, cds_overlap: bool) -> Option<String> {
        (cds_overlap && !self.frame_neutral_length).then(|| {
            format!(
                "Skipping this CDS-overlapping exon changes coding length by {} bp (not divisible by 3).",
                self.length_bp
            )
        })
    }
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
