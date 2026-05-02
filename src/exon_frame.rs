//! Shared exon length and coding-frame cue helpers.
//!
//! These helpers keep display-time exon frame cues and exon-skip planning on the
//! same vocabulary without tying headless engine code to GUI rendering.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExonLengthFrameCue {
    pub length_bp: usize,
    pub length_mod3: usize,
    pub frame_neutral_length: bool,
}

impl ExonLengthFrameCue {
    pub fn from_length(length_bp: usize) -> Self {
        let length_mod3 = length_bp % 3;
        Self {
            length_bp,
            length_mod3,
            frame_neutral_length: length_mod3 == 0,
        }
    }

    pub fn from_range(start_0based: usize, end_0based_exclusive: usize) -> Self {
        Self::from_length(end_0based_exclusive.saturating_sub(start_0based))
    }

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

    pub fn cds_phase_warning(self, cds_overlap: bool) -> Option<String> {
        (cds_overlap && !self.frame_neutral_length).then(|| {
            format!(
                "Skipping this CDS-overlapping exon changes coding length by {} bp (not divisible by 3).",
                self.length_bp
            )
        })
    }
}

pub fn phase_entry_hint(left_cds_phase: Option<u8>) -> &'static str {
    match left_cds_phase.map(|phase| phase % 3) {
        Some(0) => "starts at a codon boundary/new amino acid",
        Some(1) => "starts one nucleotide into a split codon",
        Some(2) => "starts two nucleotides into a split codon",
        None => "coding entry phase unavailable",
        _ => "coding entry phase unavailable",
    }
}
