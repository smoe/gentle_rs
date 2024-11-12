// Various constants and DNA-related functions

use crate::amino_acids::AminoAcids;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct DNAmarkerPart {
    pub length: f64,
    pub strength: Option<u64>,
}

#[derive(Clone, Debug)]
pub struct Facility {
    pub amino_acids: AminoAcids,
    // pub dna_iupac: [u8; 256],
    pub dna_iupac_complement: [u8; 256],
    pub dna_markers: HashMap<String, Vec<DNAmarkerPart>>,
}

impl Default for Facility {
    fn default() -> Self {
        Self::new()
    }
}

impl Facility {
    pub fn new() -> Self {
        Self {
            amino_acids: AminoAcids::load(),
            dna_iupac_complement: Self::initialize_dna_iupac_complement(),
            dna_markers: Self::initialize_dna_markers(),
        }
    }

    #[inline(always)]
    pub fn complement(&self, base: u8) -> u8 {
        self.dna_iupac_complement[base as usize]
    }

    #[inline(always)]
    pub fn get_dna_marker(&self, name: &str) -> Option<&Vec<DNAmarkerPart>> {
        self.dna_markers.get(name)
    }

    #[inline(always)]
    pub fn is_start_codon(&self, codon: &[u8; 3]) -> bool {
        *codon == [b'A', b'T', b'G']
    }

    #[inline(always)]
    pub fn is_stop_codon(&self, codon: &[u8; 3]) -> bool {
        *codon == [b'T', b'A', b'A'] || *codon == [b'T', b'A', b'G'] || *codon == [b'T', b'G', b'A']
    }

    fn initialize_dna_iupac_complement() -> [u8; 256] {
        let mut dna: [u8; 256] = [b' '; 256];
        dna['A' as usize] = b'T';
        dna['C' as usize] = b'G';
        dna['G' as usize] = b'C';
        dna['T' as usize] = b'A';
        dna['U' as usize] = b'A';
        // TODO IUPAC bases too?
        dna
    }

    fn initialize_dna_markers() -> HashMap<String, Vec<DNAmarkerPart>> {
        let mut ret = HashMap::new();
        let data = include_str!("../assets/dna_markers.json");
        let res: serde_json::Value = serde_json::from_str(data).expect("Invalid JSON");
        let map = res.as_object().expect("JSON is not an object");
        for (name, parts) in map.iter() {
            let parts: Vec<DNAmarkerPart> = parts
                .as_array()
                .expect("DNA marker part is not an array")
                .iter()
                .map(|p| p.as_array().expect("DNA marker subpart not an array"))
                .map(|p| DNAmarkerPart {
                    length: p[0].as_f64().expect("Not an f64"),
                    strength: match p.get(1) {
                        Some(s) => s.as_u64(),
                        None => None,
                    },
                })
                .collect();
            ret.insert(name.to_owned(), parts);
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    // use super::*;
    use crate::{iupac_code::IupacCode, FACILITY};

    // NOTE: amino_acids is tested in amino_acids.rs

    #[test]
    fn test_from_json_file() {
        let marker = FACILITY.dna_markers.get("GeneRuler Mix").unwrap();
        assert_eq!(marker[2].length, 8000.0);
        assert_eq!(marker[2].strength, Some(13));

        assert!(!IupacCode::from_letter(b'V')
            .subset(IupacCode::from_letter(b'G'))
            .is_empty());
        assert!(IupacCode::from_letter(b'H')
            .subset(IupacCode::from_letter(b'G'))
            .is_empty());
    }

    #[test]
    fn test_complement() {
        assert_eq!(FACILITY.complement(b'A'), b'T');
        assert_eq!(FACILITY.complement(b'C'), b'G');
        assert_eq!(FACILITY.complement(b'G'), b'C');
        assert_eq!(FACILITY.complement(b'T'), b'A');
        assert_eq!(FACILITY.complement(b'U'), b'A');
        assert_eq!(FACILITY.complement(b'X'), b' ');
    }

    #[test]
    fn test_get_dna_marker() {
        let marker = FACILITY.get_dna_marker("GeneRuler Mix").unwrap();
        assert_eq!(marker[2].length, 8000.0);
        assert_eq!(marker[2].strength, Some(13));
    }
}
