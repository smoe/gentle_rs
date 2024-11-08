use crate::FACILITY;
use rayon::prelude::*;

const MIN_ORF_LENGTH: i32 = 100;

#[derive(Clone, Debug, PartialEq)]
pub struct OpenReadingFrame {
    from: i32,
    to: i32,
    offset: i32,
}

impl OpenReadingFrame {
    pub fn new(from: i32, to: i32, offset: i32) -> Self {
        OpenReadingFrame { from, to, offset }
    }

    pub fn from(&self) -> i32 {
        self.from
    }

    pub fn to(&self) -> i32 {
        self.to
    }

    pub fn offset(&self) -> i32 {
        self.offset
    }

    pub fn find_orfs(sequence: &[u8], is_circular: bool) -> Vec<OpenReadingFrame> {
        [1, 2, 3, -1, -2, -3]
            .par_iter()
            .flat_map(|offset| Self::add_orfs(sequence, is_circular, *offset))
            .collect()
    }

    fn get_nucleotide(sequence: &[u8], pos: i32, complement: bool) -> char {
        match sequence.get(pos as usize) {
            Some(c) => {
                let c = *c as char;
                if complement {
                    FACILITY.complement(c)
                } else {
                    c
                }
            }
            None => 'N',
        }
    }

    fn add_orfs(sequence: &[u8], is_circular: bool, offset: i32) -> Vec<OpenReadingFrame> {
        let mut ret = vec![];
        let seq_len = sequence.len() as i32;
        let (direction, inisial_start_position, complement) = if offset > 0 {
            (1, offset - 1, false)
        } else {
            let mut b = 0;
            while b <= seq_len {
                b += 3;
            }
            b -= 3;
            b += offset + 1;
            (-1, b, true)
        };

        // Max ORF length, in AA => ~ 1/3 of the sequence length
        let max_amino_acids = (seq_len + 2) / 3;

        let mut start_codon_position = inisial_start_position;
        while start_codon_position + direction * 3 > 0
            && (start_codon_position + direction * 3) < seq_len
        {
            let codon = [
                Self::get_nucleotide(sequence, start_codon_position, complement),
                Self::get_nucleotide(sequence, start_codon_position + direction, complement),
                Self::get_nucleotide(sequence, start_codon_position + direction * 2, complement),
            ];

            if codon == ['A', 'T', 'G'] {
                let mut amino_acids = 0;
                let mut stop_codon_position = start_codon_position;

                while amino_acids <= max_amino_acids {
                    amino_acids += 1;

                    // Handle wrapping
                    if stop_codon_position < 0 {
                        if !is_circular {
                            break;
                        }
                        stop_codon_position = seq_len - stop_codon_position;
                    }
                    if stop_codon_position >= seq_len {
                        if !is_circular {
                            break;
                        }
                        stop_codon_position -= seq_len;
                    }

                    let codon = [
                        Self::get_nucleotide(sequence, stop_codon_position, complement),
                        Self::get_nucleotide(sequence, stop_codon_position + direction, complement),
                        Self::get_nucleotide(
                            sequence,
                            stop_codon_position + direction * 2,
                            complement,
                        ),
                    ];

                    if codon == ['T', 'A', 'A']
                        || codon == ['T', 'A', 'G']
                        || codon == ['T', 'G', 'A']
                    {
                        if amino_acids >= MIN_ORF_LENGTH {
                            let from = start_codon_position;
                            let to = stop_codon_position + direction * 2;

                            if from < to || direction == 1 {
                                ret.push(OpenReadingFrame::new(from, to, offset));
                            } else {
                                ret.push(OpenReadingFrame::new(to, from, offset));
                            }
                        }
                        break;
                    }
                    stop_codon_position += direction * 3;
                }
            }
            start_codon_position += direction * 3;
        }
        ret
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

    // TODO: Fix this test
    // #[test]
    // fn test_find_orfs_circular_reverse() {
    //     let mut sequence = "GGGCATAGGG".to_string();
    //     sequence += "GGGGGGTTAGGG";
    //     sequence += &"CCC".repeat(105);

    //     // Try linear
    //     let orfs = OpenReadingFrame::find_orfs(sequence.as_bytes(), false);
    //     assert!(orfs.is_empty());

    //     // Try circular
    //     let orfs = OpenReadingFrame::find_orfs(sequence.as_bytes(), true);
    //     // assert_eq!(orfs, vec![OpenReadingFrame::new(15, 8, 1)]);
    //     println!("{orfs:?} / {}", sequence.len());
    // }
}
