// Various constants and DNA-related functions

use crate::amino_acids::AminoAcids;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

const DNA_A: u8 = 1;
const DNA_C: u8 = 2;
const DNA_G: u8 = 4;
const DNA_T: u8 = 8;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct DNAmarkerPart {
    pub length: f64,
    pub strength: Option<u64>,
}

#[derive(Clone, Debug)]
pub struct Facility {
    pub amino_acids: AminoAcids,
    pub dna_iupac: [u8; 256],
    pub dna_iupac_complement: [char; 256],
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
            dna_iupac: Self::initialize_dna_iupac(),
            dna_iupac_complement: Self::initialize_dna_iupac_complement(),
            dna_markers: Self::initialize_dna_markers(),
        }
    }

    #[inline(always)]
    pub fn complement(&self, base: char) -> char {
        self.dna_iupac_complement[base as usize]
    }

    #[inline(always)]
    pub fn base2iupac(&self, base: char) -> u8 {
        self.dna_iupac
            .get(base.to_ascii_uppercase() as usize)
            .copied()
            .unwrap_or(0)
    }

    #[inline(always)]
    pub fn split_iupac(&self, iupac_code: u8) -> Vec<char> {
        let mut ret = Vec::with_capacity(4);
        if iupac_code & DNA_A != 0 {
            ret.push('A');
        }
        if iupac_code & DNA_C != 0 {
            ret.push('C');
        }
        if iupac_code & DNA_G != 0 {
            ret.push('G');
        }
        if iupac_code & DNA_T != 0 {
            ret.push('T');
        }
        ret
    }

    pub fn get_dna_marker(&self, name: &str) -> Option<&Vec<DNAmarkerPart>> {
        self.dna_markers.get(name)
    }

    fn initialize_dna_iupac() -> [u8; 256] {
        let mut dna: [u8; 256] = [0; 256];
        dna['A' as usize] = DNA_A;
        dna['C' as usize] = DNA_C;
        dna['G' as usize] = DNA_G;
        dna['T' as usize] = DNA_T;
        dna['U' as usize] = DNA_T;
        dna['W' as usize] = DNA_A | DNA_T;
        dna['S' as usize] = DNA_C | DNA_G;
        dna['M' as usize] = DNA_A | DNA_C;
        dna['K' as usize] = DNA_G | DNA_T;
        dna['R' as usize] = DNA_A | DNA_G;
        dna['Y' as usize] = DNA_C | DNA_T;
        dna['B' as usize] = DNA_C | DNA_G | DNA_T;
        dna['D' as usize] = DNA_A | DNA_G | DNA_T;
        dna['H' as usize] = DNA_A | DNA_C | DNA_T;
        dna['V' as usize] = DNA_A | DNA_C | DNA_G;
        dna['N' as usize] = DNA_A | DNA_C | DNA_G | DNA_T;
        dna
    }

    fn initialize_dna_iupac_complement() -> [char; 256] {
        let mut dna: [char; 256] = [' '; 256];
        dna['A' as usize] = 'T';
        dna['C' as usize] = 'G';
        dna['G' as usize] = 'C';
        dna['T' as usize] = 'A';
        dna['U' as usize] = 'A';
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
    use super::*;
    use crate::FACILITY;

    // NOTE: amino_acids is tested in amino_acids.rs

    #[test]
    fn test_from_json_file() {
        let marker = FACILITY.dna_markers.get("GeneRuler Mix").unwrap();
        assert_eq!(marker[2].length, 8000.0);
        assert_eq!(marker[2].strength, Some(13));

        assert!(FACILITY.dna_iupac['V' as usize] & FACILITY.dna_iupac['G' as usize] > 0);
        assert!(FACILITY.dna_iupac['H' as usize] & FACILITY.dna_iupac['G' as usize] == 0);
    }

    #[test]
    fn test_complement() {
        assert_eq!(FACILITY.complement('A'), 'T');
        assert_eq!(FACILITY.complement('C'), 'G');
        assert_eq!(FACILITY.complement('G'), 'C');
        assert_eq!(FACILITY.complement('T'), 'A');
        assert_eq!(FACILITY.complement('U'), 'A');
        assert_eq!(FACILITY.complement('X'), ' ');
    }

    #[test]
    fn test_base2iupac() {
        assert_eq!(FACILITY.base2iupac('A'), DNA_A);
        assert_eq!(FACILITY.base2iupac('C'), DNA_C);
        assert_eq!(FACILITY.base2iupac('G'), DNA_G);
        assert_eq!(FACILITY.base2iupac('T'), DNA_T);
        assert_eq!(FACILITY.base2iupac('U'), DNA_T);
        assert_eq!(FACILITY.base2iupac('X'), 0);
    }

    #[test]
    fn test_split_iupac() {
        assert_eq!(FACILITY.split_iupac(DNA_A), vec!['A']);
        assert_eq!(FACILITY.split_iupac(DNA_C), vec!['C']);
        assert_eq!(FACILITY.split_iupac(DNA_G), vec!['G']);
        assert_eq!(FACILITY.split_iupac(DNA_T), vec!['T']);
        assert_eq!(FACILITY.split_iupac(DNA_A | DNA_C), vec!['A', 'C']);
        assert_eq!(
            FACILITY.split_iupac(DNA_A | DNA_C | DNA_G),
            vec!['A', 'C', 'G']
        );
        assert_eq!(
            FACILITY.split_iupac(DNA_A | DNA_C | DNA_G | DNA_T),
            vec!['A', 'C', 'G', 'T']
        );
    }

    #[test]
    fn test_get_dna_marker() {
        let marker = FACILITY.get_dna_marker("GeneRuler Mix").unwrap();
        assert_eq!(marker[2].length, 8000.0);
        assert_eq!(marker[2].strength, Some(13));
    }
}
