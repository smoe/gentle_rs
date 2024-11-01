use serde::{Deserialize, Serialize};

use crate::dna_sequence::DNAsequence;

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

#[derive(Clone, Debug)]
pub struct RestrictionEnzymeSite {
    pub offset: usize,
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
        let forward = seq.forward();
        let seq_len = if seq.is_circular() {
            forward.len()
        } else {
            forward.len() - recognition_len + 1
        };
        for start in 0..seq_len {
            let range = std::ops::Range {
                start,
                end: start + recognition_len,
            };
            let s = forward.get(range); // TODO circular
            if let Some(s) = s {
                let s = std::str::from_utf8(s).unwrap();
                if s == self.sequence {
                    // TODO IUPAC
                    ret.push(RestrictionEnzymeSite {
                        offset: start,
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
