//! IUPAC nucleotide-code conversion and validation helpers.

const DNA_BITMASK_A: u8 = 1;
const DNA_BITMASK_C: u8 = 2;
const DNA_BITMASK_G: u8 = 4;
const DNA_BITMASK_T: u8 = 8;
const DNA_BITMASK_W: u8 = DNA_BITMASK_A | DNA_BITMASK_T;
const DNA_BITMASK_S: u8 = DNA_BITMASK_C | DNA_BITMASK_G;
const DNA_BITMASK_M: u8 = DNA_BITMASK_A | DNA_BITMASK_C;
const DNA_BITMASK_K: u8 = DNA_BITMASK_G | DNA_BITMASK_T;
const DNA_BITMASK_R: u8 = DNA_BITMASK_A | DNA_BITMASK_G;
const DNA_BITMASK_Y: u8 = DNA_BITMASK_C | DNA_BITMASK_T;
const DNA_BITMASK_B: u8 = DNA_BITMASK_C | DNA_BITMASK_G | DNA_BITMASK_T;
const DNA_BITMASK_D: u8 = DNA_BITMASK_A | DNA_BITMASK_G | DNA_BITMASK_T;
const DNA_BITMASK_H: u8 = DNA_BITMASK_A | DNA_BITMASK_C | DNA_BITMASK_T;
const DNA_BITMASK_V: u8 = DNA_BITMASK_A | DNA_BITMASK_C | DNA_BITMASK_G;
const DNA_BITMASK_N: u8 = DNA_BITMASK_A | DNA_BITMASK_C | DNA_BITMASK_G | DNA_BITMASK_T;

/// A bitmasked IUPAC code for DNA bases, eg DNA_BITMASK_A|DNA_BITMASK_C
#[derive(Debug, Copy, Clone, PartialEq, Hash)]
pub struct IupacCode(u8);

impl IupacCode {
    pub fn new(bitmask: u8) -> Self {
        Self(bitmask)
    }

    #[inline(always)]
    pub fn from_letter(letter: u8) -> Self {
        match letter.to_ascii_uppercase() {
            b'A' => Self(DNA_BITMASK_A),
            b'C' => Self(DNA_BITMASK_C),
            b'G' => Self(DNA_BITMASK_G),
            b'T' => Self(DNA_BITMASK_T),
            b'U' => Self(DNA_BITMASK_T),
            b'W' => Self(DNA_BITMASK_A | DNA_BITMASK_T),
            b'S' => Self(DNA_BITMASK_C | DNA_BITMASK_G),
            b'M' => Self(DNA_BITMASK_A | DNA_BITMASK_C),
            b'K' => Self(DNA_BITMASK_G | DNA_BITMASK_T),
            b'R' => Self(DNA_BITMASK_A | DNA_BITMASK_G),
            b'Y' => Self(DNA_BITMASK_C | DNA_BITMASK_T),
            b'B' => Self(DNA_BITMASK_C | DNA_BITMASK_G | DNA_BITMASK_T),
            b'D' => Self(DNA_BITMASK_A | DNA_BITMASK_G | DNA_BITMASK_T),
            b'H' => Self(DNA_BITMASK_A | DNA_BITMASK_C | DNA_BITMASK_T),
            b'V' => Self(DNA_BITMASK_A | DNA_BITMASK_C | DNA_BITMASK_G),
            b'N' => Self(DNA_BITMASK_N),
            _ => Self(0),
        }
    }

    #[inline(always)]
    pub fn to_letter(&self) -> u8 {
        match self.0 {
            DNA_BITMASK_A => b'A',
            DNA_BITMASK_C => b'C',
            DNA_BITMASK_G => b'G',
            DNA_BITMASK_T => b'T',
            DNA_BITMASK_N => b'N',
            DNA_BITMASK_W => b'W',
            DNA_BITMASK_S => b'S',
            DNA_BITMASK_M => b'M',
            DNA_BITMASK_K => b'K',
            DNA_BITMASK_R => b'R',
            DNA_BITMASK_Y => b'Y',
            DNA_BITMASK_B => b'B',
            DNA_BITMASK_D => b'D',
            DNA_BITMASK_H => b'H',
            DNA_BITMASK_V => b'V',
            _ => b' ',
        }
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.0 == 0
    }

    #[inline(always)]
    pub fn subset(self, other: Self) -> Self {
        Self(self.0 & other.0)
    }

    #[inline(always)]
    pub fn is_valid_letter(letter: u8) -> bool {
        matches!(
            letter,
            b'A' | b'C'
                | b'G'
                | b'T'
                | b'U'
                | b'W'
                | b'S'
                | b'M'
                | b'K'
                | b'R'
                | b'Y'
                | b'B'
                | b'D'
                | b'H'
                | b'V'
                | b'N'
                | b'a'
                | b'c'
                | b'g'
                | b't'
                | b'u'
                | b'w'
                | b's'
                | b'm'
                | b'k'
                | b'r'
                | b'y'
                | b'b'
                | b'd'
                | b'h'
                | b'v'
                | b'n'
        )
    }

    #[inline(always)]
    pub fn to_vec(&self) -> Vec<u8> {
        let mut ret = Vec::with_capacity(4);
        if self.0 & DNA_BITMASK_A != 0 {
            ret.push(b'A');
        }
        if self.0 & DNA_BITMASK_C != 0 {
            ret.push(b'C');
        }
        if self.0 & DNA_BITMASK_G != 0 {
            ret.push(b'G');
        }
        if self.0 & DNA_BITMASK_T != 0 {
            ret.push(b'T');
        }
        ret
    }

    #[inline(always)]
    pub fn letter_complement(letter: u8) -> u8 {
        match letter.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'U' => b'A',
            _ => b' ',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base2iupac() {
        assert!(
            !IupacCode::from_letter(b'V')
                .subset(IupacCode::from_letter(b'G'))
                .is_empty()
        );
        assert!(
            IupacCode::from_letter(b'H')
                .subset(IupacCode::from_letter(b'G'))
                .is_empty()
        );
        assert_eq!(IupacCode::from_letter(b'A'), IupacCode::new(DNA_BITMASK_A));
        assert_eq!(IupacCode::from_letter(b'C'), IupacCode::new(DNA_BITMASK_C));
        assert_eq!(IupacCode::from_letter(b'G'), IupacCode::new(DNA_BITMASK_G));
        assert_eq!(IupacCode::from_letter(b'T'), IupacCode::new(DNA_BITMASK_T));
        assert_eq!(IupacCode::from_letter(b'U'), IupacCode::new(DNA_BITMASK_T));
        assert_eq!(
            IupacCode::from_letter(b'M'),
            IupacCode::new(DNA_BITMASK_A | DNA_BITMASK_C)
        );
        assert_eq!(IupacCode::from_letter(b'X'), IupacCode::new(0));
    }

    #[test]
    fn test_split_iupac() {
        assert_eq!(IupacCode::new(DNA_BITMASK_A).to_vec(), vec![b'A']);
        assert_eq!(IupacCode::new(DNA_BITMASK_C).to_vec(), vec![b'C']);
        assert_eq!(IupacCode::new(DNA_BITMASK_G).to_vec(), vec![b'G']);
        assert_eq!(IupacCode::new(DNA_BITMASK_T).to_vec(), vec![b'T']);
        assert_eq!(
            IupacCode::new(DNA_BITMASK_A | DNA_BITMASK_C).to_vec(),
            vec![b'A', b'C']
        );
        assert_eq!(
            IupacCode::new(DNA_BITMASK_A | DNA_BITMASK_C | DNA_BITMASK_G).to_vec(),
            vec![b'A', b'C', b'G']
        );
        assert_eq!(
            IupacCode::new(DNA_BITMASK_A | DNA_BITMASK_C | DNA_BITMASK_G | DNA_BITMASK_T).to_vec(),
            vec![b'A', b'C', b'G', b'T']
        );
    }

    #[test]
    fn test_complement() {
        assert_eq!(IupacCode::letter_complement(b'A'), b'T');
        assert_eq!(IupacCode::letter_complement(b'C'), b'G');
        assert_eq!(IupacCode::letter_complement(b'G'), b'C');
        assert_eq!(IupacCode::letter_complement(b'T'), b'A');
        assert_eq!(IupacCode::letter_complement(b'U'), b'A');
        assert_eq!(IupacCode::letter_complement(b'X'), b' ');
        assert_eq!(IupacCode::letter_complement(b'a'), b'T');
    }

    #[test]
    fn test_to_letter() {
        assert_eq!(IupacCode::from_letter(b'A').to_letter(), b'A');
        assert_eq!(IupacCode::from_letter(b'C').to_letter(), b'C');
        assert_eq!(IupacCode::from_letter(b'G').to_letter(), b'G');
        assert_eq!(IupacCode::from_letter(b'T').to_letter(), b'T');
        assert_eq!(IupacCode::from_letter(b'U').to_letter(), b'T');
        assert_eq!(IupacCode::from_letter(b'M').to_letter(), b'M');
        assert_eq!(IupacCode::from_letter(b'X').to_letter(), b' ');
    }
}
