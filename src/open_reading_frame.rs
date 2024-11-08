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
        let (dir, b, complement) = if offset > 0 {
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

        let max_aa = (seq_len + 2) / 3; // Max ORF length, in AA => ~ 1/3 of the sequence length

        let mut a = b;
        while a + dir * 3 > 0 && (a + dir * 3) < seq_len {
            let codon = format!(
                "{}{}{}",
                Self::get_nucleotide(sequence, a, complement),
                Self::get_nucleotide(sequence, a + dir, complement),
                Self::get_nucleotide(sequence, a + dir * 2, complement)
            );

            if codon == "ATG" {
                // let mut cnt = (seq_len + 2) / 3; // Max ORF length, in AA => ~ 1/3 of the sequence length
                let mut aa = 0;
                let mut b = a;

                while aa <= max_aa {
                    // cnt -= 1;
                    aa += 1;

                    // Handle wrapping
                    if b < 0 {
                        if !is_circular {
                            break;
                        }
                        b = seq_len - b;
                    }
                    if b >= seq_len {
                        if !is_circular {
                            break;
                        }
                        b -= seq_len;
                    }

                    let codon = format!(
                        "{}{}{}",
                        Self::get_nucleotide(sequence, b, complement),
                        Self::get_nucleotide(sequence, b + dir, complement),
                        Self::get_nucleotide(sequence, b + dir * 2, complement)
                    );
                    // println!("{b}: {codon}");

                    if codon == "TAA" || codon == "TAG" || codon == "TGA" {
                        if aa >= MIN_ORF_LENGTH {
                            let from = a;
                            let to = b + dir * 2;

                            if from < to || dir == 1 {
                                ret.push(OpenReadingFrame::new(from, to, offset));
                            } else {
                                ret.push(OpenReadingFrame::new(to, from, offset));
                            }
                        }
                        break;
                    }
                    b += dir * 3;
                }
            }

            a += dir * 3;
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_orfs_linear_forward() {
        let mut sequence = "AAATG".to_string();
        sequence += &"AAA".repeat(105); // Filler
        sequence += "TAAGG";
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
