use crate::dna_sequence::DNAsequence;
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RestrictionEnzymeKey {
    pos: isize,
    cut_size: isize,
    number_of_cuts: usize,
    from: isize,
    to: isize,
}

impl RestrictionEnzymeKey {
    pub fn new(pos: isize, cut_size: isize, number_of_cuts: usize, from: isize, to: isize) -> Self {
        Self {
            pos,
            cut_size,
            number_of_cuts,
            from,
            to,
        }
    }

    pub fn number_of_cuts(&self) -> usize {
        self.number_of_cuts
    }

    pub fn cut_size(&self) -> isize {
        self.cut_size
    }

    pub fn pos(&self) -> isize {
        self.pos
    }

    pub fn from(&self) -> isize {
        self.from
    }

    pub fn to(&self) -> isize {
        self.to
    }
}

impl PartialOrd for RestrictionEnzymeKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RestrictionEnzymeKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.pos.cmp(&other.pos)
    }
}

impl fmt::Display for RestrictionEnzymeKey {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let j = serde_json::to_string(&self).unwrap();
        write!(f, "{}", j)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RestrictionEnzyme {
    pub name: String,
    pub sequence: String,
    pub note: Option<String>,
    pub cut: isize,
    pub overlap: isize,
    #[serde(skip_serializing, default)]
    is_palindromic: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RestrictionEnzymeSite {
    pub offset: isize,
    pub enzyme: RestrictionEnzyme,
    pub forward_strand: bool,
}

impl RestrictionEnzyme {
    pub fn check_palimdromic(&mut self) {
        self.is_palindromic = self.sequence == self.get_sequence_rc();
    }

    pub fn is_palindromic(&self) -> bool {
        self.is_palindromic
    }

    fn get_sequence_rc(&self) -> String {
        // TODO cache this?
        let rc = self.sequence.as_bytes();
        let rc = match std::str::from_utf8(rc) {
            Ok(rc) => rc,
            Err(_) => panic!("RestrictionEnzyme::check_palimdromic: non-utf8 char"),
        };
        rc.to_string()
    }

    pub fn get_sites(
        &self,
        seq: &DNAsequence,
        max_sites: Option<usize>,
    ) -> Vec<RestrictionEnzymeSite> {
        // TODO reverse-complement if required
        let mut ret = vec![];
        let recognition_len = self.sequence.len();
        let seq_len = if seq.is_circular() {
            seq.len()
        } else {
            seq.len() - recognition_len + 1
        };
        for start in 0..seq_len {
            let range = std::ops::Range {
                start,
                end: start + recognition_len,
            };
            let s = seq.get_range_safe(range); // Safe for circular
            if let Some(s) = s {
                let s = std::str::from_utf8(&s).unwrap().to_uppercase(); // TODO do this once?
                if s == self.sequence {
                    // TODO IUPAC
                    ret.push(RestrictionEnzymeSite {
                        offset: start as isize,
                        enzyme: self.to_owned(),
                        forward_strand: true,
                    });
                }
            }
        }
        if let Some(max) = max_sites {
            if max < ret.len() {
                return vec![];
            }
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence::DNAsequence;

    #[test]
    fn test_restriction_enzyme() {
        let mut re = RestrictionEnzyme {
            name: "EcoRI".to_string(),
            sequence: "GAATTC".to_string(),
            note: None,
            cut: 1,
            overlap: 1,
            is_palindromic: false,
        };
        re.check_palimdromic();
        assert!(re.is_palindromic());
        let seq = DNAsequence::from_sequence("GAATTC").unwrap();
        let sites = re.get_sites(&seq, None);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].offset, 0);
        assert!(sites[0].forward_strand);
    }

    #[test]
    fn test_restriction_enzyme_sites() {
        let mut re = RestrictionEnzyme {
            name: "EcoRI".to_string(),
            sequence: "GAATTC".to_string(),
            note: None,
            cut: 1,
            overlap: 1,
            is_palindromic: false,
        };
        re.check_palimdromic();
        assert!(re.is_palindromic());
        let seq = DNAsequence::from_sequence("GAATTCGAATTC").unwrap();
        let sites = re.get_sites(&seq, None);
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].offset, 0);
        assert!(sites[0].forward_strand);
        assert_eq!(sites[1].offset, 6);
        assert!(sites[1].forward_strand);
    }
}
