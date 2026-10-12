//! Open-reading-frame detection logic.

use crate::{amino_acids::AminoAcids, iupac_code::IupacCode};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

const MIN_ORF_LENGTH: i32 = 100;

/// ATG-to-first-stop prediction, not an annotated CDS or evidence of translation.
///
/// `from`/`to` are inclusive, 0-based bounds in reference-strand order, including
/// the stop codon. `from > to` denotes an origin crossing on either strand.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OpenReadingFrame {
    from: i32,
    to: i32,
    frame: i32,
}

impl OpenReadingFrame {
    pub fn new(from: i32, to: i32, frame: i32) -> Self {
        OpenReadingFrame { from, to, frame }
    }

    #[inline(always)]
    pub fn from(&self) -> i32 {
        self.from
    }

    #[inline(always)]
    pub fn to(&self) -> i32 {
        self.to
    }

    #[inline(always)]
    pub fn frame(&self) -> i32 {
        self.frame
    }

    #[inline(always)]
    pub fn is_reverse(&self) -> bool {
        self.frame < 0
    }

    /// Half-open reference bounds; the end exceeds `sequence_length` for a wrap.
    pub fn unrolled_bounds_0based(&self, sequence_length: usize) -> Option<(usize, usize)> {
        let from = usize::try_from(self.from).ok()?;
        let to = usize::try_from(self.to).ok()?;
        if from >= sequence_length || to >= sequence_length {
            return None;
        }
        let end = to.checked_add(1)?;
        Some((
            from,
            if from > to {
                end.checked_add(sequence_length)?
            } else {
                end
            },
        ))
    }

    /// One or two half-open spans, retaining the origin crossing without sorting.
    pub fn spans_0based(&self, sequence_length: usize) -> impl Iterator<Item = (usize, usize)> {
        let spans = match self.unrolled_bounds_0based(sequence_length) {
            Some((from, end)) if end > sequence_length => [
                Some((from, sequence_length)),
                Some((0, end - sequence_length)),
            ],
            Some(bounds) => [Some(bounds), None],
            None => [None, None],
        };
        spans.into_iter().flatten()
    }

    pub fn find_orfs(sequence: &[u8], is_circular: bool) -> Vec<OpenReadingFrame> {
        // The serialized coordinate fields remain i32 for compatibility.
        if sequence.len() < 3 || sequence.len() > i32::MAX as usize {
            return Vec::new();
        }
        [1, 2, 3, -1, -2, -3]
            .par_iter()
            .map(|offset| Self::add_orfs(sequence, is_circular, *offset))
            .collect::<Vec<_>>()
            .into_iter()
            .flatten()
            .collect()
    }

    fn add_orfs(sequence: &[u8], is_circular: bool, offset: i32) -> Vec<OpenReadingFrame> {
        let mut ret = vec![];
        let seq_len = sequence.len() as i64;
        let reverse = offset < 0;
        let direction = if reverse { -1 } else { 1 };
        // Retain the existing reference-anchored frame numbering on both strands.
        let phase = if reverse {
            (offset + 1).rem_euclid(3)
        } else {
            offset - 1
        } as i64;
        let mut start = if reverse {
            seq_len - 1 - (seq_len - 1 - phase).rem_euclid(3)
        } else {
            phase
        };
        while (0..seq_len).contains(&start) {
            if Self::get_codon(sequence, start, reverse, is_circular)
                .is_some_and(|codon| AminoAcids::is_start_codon(&codon))
            {
                // Include both start and stop, but never reuse a base on a second lap.
                for step in 1..seq_len / 3 {
                    let stop = start + direction * step * 3;
                    let Some(codon) = Self::get_codon(sequence, stop, reverse, is_circular) else {
                        break;
                    };
                    if AminoAcids::is_stop_codon(&codon) {
                        if step + 1 >= i64::from(MIN_ORF_LENGTH) {
                            let last = (stop + direction * 2).rem_euclid(seq_len) as i32;
                            let (from, to) = if reverse {
                                (last, start as i32)
                            } else {
                                (start as i32, last)
                            };
                            ret.push(Self::new(from, to, offset));
                        }
                        break;
                    }
                }
            }
            start += direction * 3;
        }
        ret
    }

    #[inline(always)]
    fn get_codon(sequence: &[u8], start: i64, reverse: bool, circular: bool) -> Option<[u8; 3]> {
        let seq_len = i64::try_from(sequence.len()).ok().filter(|len| *len > 0)?;
        let mut codon = [b'N'; 3];
        for (index, base) in codon.iter_mut().enumerate() {
            let pos = start
                + if reverse {
                    -(index as i64)
                } else {
                    index as i64
                };
            let pos = if circular {
                pos.rem_euclid(seq_len)
            } else {
                pos
            };
            let letter = *sequence.get(usize::try_from(pos).ok()?)?;
            *base = if reverse {
                IupacCode::from_letter(letter).complement().to_letter()
            } else {
                letter.to_ascii_uppercase()
            };
        }
        Some(codon)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_orfs_linear_forward_taa() {
        let mut sequence = "AAATG".to_string();
        sequence += &"AAA".repeat(105); // Filler
        sequence += "TAAGG";
        let orfs = OpenReadingFrame::find_orfs(sequence.as_bytes(), false);
        assert_eq!(orfs, vec![OpenReadingFrame::new(2, 322, 3)]);
    }

    #[test]
    fn test_find_orfs_linear_forward_tag() {
        let mut sequence = "AAATG".to_string();
        sequence += &"AAA".repeat(105); // Filler
        sequence += "TAGGG";
        let orfs = OpenReadingFrame::find_orfs(sequence.as_bytes(), false);
        assert_eq!(orfs, vec![OpenReadingFrame::new(2, 322, 3)]);
    }

    #[test]
    fn test_find_orfs_linear_forward_tga() {
        let mut sequence = "AAATG".to_string();
        sequence += &"AAA".repeat(105); // Filler
        sequence += "TGAGG";
        let orfs = OpenReadingFrame::find_orfs(sequence.as_bytes(), false);
        assert_eq!(orfs, vec![OpenReadingFrame::new(2, 322, 3)]);
    }

    #[test]
    fn test_find_orfs_linear_reverse() {
        let mut sequence = "AATTA".to_string();
        sequence += &"AAA".repeat(105); // Filler
        sequence += "CATGG";
        let orfs = OpenReadingFrame::find_orfs(sequence.as_bytes(), false);
        assert_eq!(orfs, vec![OpenReadingFrame::new(2, 322, -3)]);
    }

    #[test]
    fn test_find_orfs_circular_forward() {
        let mut sequence = "CCCCCCTAA".to_string();
        sequence += "GGGGGGATG";
        sequence += &"CCC".repeat(105);

        // Try linear
        let orfs = OpenReadingFrame::find_orfs(sequence.as_bytes(), false);
        assert!(orfs.is_empty());

        // Try circular
        let orfs = OpenReadingFrame::find_orfs(sequence.as_bytes(), true);
        assert_eq!(orfs, vec![OpenReadingFrame::new(15, 8, 1)]);
    }

    #[test]
    fn test_find_orfs_circular_reverse() {
        let sequence = format!("{}{}{}", "GGGCATAGGG", "GGGGGGTTAGGG", "CCC".repeat(105));
        assert!(OpenReadingFrame::find_orfs(sequence.as_bytes(), false).is_empty());
        assert_eq!(
            OpenReadingFrame::find_orfs(sequence.as_bytes(), true),
            vec![OpenReadingFrame::new(16, 5, -2)]
        );
    }

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                _ => b'N',
            })
            .collect()
    }

    #[test]
    fn circular_orfs_preserve_sequence_at_every_origin_and_strand() {
        // Synthetic single ORF; rotation covers every split-start/split-stop case
        // and all molecule-length residues modulo three, without external data.
        let coding = format!("ATG{}TAA", "CCC".repeat(105)).into_bytes();
        for extra in 0..3 {
            let mut template = coding.clone();
            template.extend(vec![b'C'; 6 + extra]);
            let len = template.len();
            for rotation in 0..len {
                let mut forward = template.clone();
                forward.rotate_left(rotation);
                let start = (len - rotation) % len;
                let stop = (start + coding.len() - 1) % len;
                for reverse in [false, true] {
                    let sequence = if reverse {
                        reverse_complement(&forward)
                    } else {
                        forward.clone()
                    };
                    let result = OpenReadingFrame::find_orfs(&sequence, true);
                    assert_eq!(
                        result.len(),
                        1,
                        "rotation={rotation}, reverse={reverse}, len={len}"
                    );
                    let orf = &result[0];
                    let bounds = if reverse {
                        (len - 1 - stop, len - 1 - start)
                    } else {
                        (start, stop)
                    };
                    assert_eq!((orf.from(), orf.to()), (bounds.0 as i32, bounds.1 as i32));
                    assert_eq!(orf.is_reverse(), reverse);
                    let extracted: Vec<u8> = orf
                        .spans_0based(len)
                        .flat_map(|(a, b)| sequence[a..b].iter().copied())
                        .collect();
                    assert_eq!(
                        if reverse {
                            reverse_complement(&extracted)
                        } else {
                            extracted
                        },
                        coding
                    );
                    let (from, end) = orf.unrolled_bounds_0based(len).expect("valid bounds");
                    assert_eq!(end - from, coding.len());
                }
            }
        }
    }

    #[test]
    fn orfs_at_linear_edges_and_exact_one_lap_are_complete() {
        let forward = format!("ATG{}TGA", "CCC".repeat(105)).into_bytes();
        for sequence in [forward.clone(), reverse_complement(&forward)] {
            for circular in [false, true] {
                let result = OpenReadingFrame::find_orfs(&sequence, circular);
                assert_eq!(result.len(), 1);
                assert_eq!(
                    (result[0].from(), result[0].to()),
                    (0, sequence.len() as i32 - 1)
                );
            }
        }
    }

    #[test]
    fn orfs_keep_first_stop_threshold_and_do_not_invent_missing_codons() {
        let short = format!("ATG{}TAACCCTAA", "CCC".repeat(97));
        let at_threshold = format!("ATG{}TAA", "CCC".repeat(98));
        assert!(OpenReadingFrame::find_orfs(short.as_bytes(), false).is_empty());
        assert_eq!(
            OpenReadingFrame::find_orfs(at_threshold.as_bytes(), false).len(),
            1
        );
        for sequence in [
            "".to_owned(),
            "A".to_owned(),
            "AT".to_owned(),
            format!("ATG{}TA", "CCC".repeat(105)),
            format!("ATN{}TAA", "CCC".repeat(105)),
        ] {
            assert!(OpenReadingFrame::find_orfs(sequence.as_bytes(), false).is_empty());
        }
        let no_stop = format!("ATG{}", "CCC".repeat(105));
        assert!(OpenReadingFrame::find_orfs(no_stop.as_bytes(), true).is_empty());
        assert_eq!(
            OpenReadingFrame::find_orfs(at_threshold.to_lowercase().as_bytes(), false).len(),
            1
        );
    }

    #[test]
    fn orf_spans_retain_last_base_and_reject_invalid_legacy_bounds() {
        assert_eq!(
            OpenReadingFrame::new(90, 9, -1)
                .spans_0based(100)
                .collect::<Vec<_>>(),
            vec![(90, 100), (0, 10)]
        );
        assert_eq!(
            OpenReadingFrame::new(0, 99, 1)
                .spans_0based(100)
                .collect::<Vec<_>>(),
            vec![(0, 100)]
        );
        for orf in [
            OpenReadingFrame::new(-1, 5, 1),
            OpenReadingFrame::new(10, 100, 1),
        ] {
            assert_eq!(orf.unrolled_bounds_0based(100), None);
        }
    }
}
