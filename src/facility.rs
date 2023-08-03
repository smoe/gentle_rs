use std::{collections::HashMap, fs};
use serde::{Deserialize, Serialize};

use crate::amino_acids::AminoAcids;

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
    pub dna_markers: HashMap<String,Vec<DNAmarkerPart>>,
}

impl Facility {
    pub fn new() -> Self {
        Self {
            amino_acids: AminoAcids::load(),
            dna_iupac: Self::get_dna(),
            dna_markers: Self::get_dna_markers(),
        }
    }

    fn get_dna() -> [u8; 256] {
        let mut dna: [u8; 256] = [0; 256];
        dna['A' as usize] = DNA_A;
        dna['C' as usize] = DNA_C;
        dna['G' as usize] = DNA_G;
        dna['T' as usize] = DNA_T;
        dna['U' as usize] = DNA_T;
        dna['W' as usize] = DNA_A|DNA_T;
        dna['S' as usize] = DNA_C|DNA_G;
        dna['M' as usize] = DNA_A|DNA_C;
        dna['K' as usize] = DNA_G|DNA_T;
        dna['R' as usize] = DNA_A|DNA_G;
        dna['Y' as usize] = DNA_C|DNA_T;
        dna['B' as usize] = DNA_C|DNA_G|DNA_T;
        dna['D' as usize] = DNA_A|DNA_G|DNA_T;
        dna['H' as usize] = DNA_A|DNA_C|DNA_T;
        dna['V' as usize] = DNA_A|DNA_C|DNA_G;
        dna['N' as usize] = DNA_A|DNA_C|DNA_G|DNA_T;
        dna
    }

    fn get_dna_markers() -> HashMap<String,Vec<DNAmarkerPart>>{
        let mut ret = HashMap::new();
        let data = fs::read_to_string("assets/dna_markers.json").expect("File not found");
        let res: serde_json::Value = serde_json::from_str(&data).expect("Invalid JSON");
        let map = res.as_object().expect("JSON is not an object");
        for (name,parts) in map.iter() {
            let parts: Vec<DNAmarkerPart> = parts
                .as_array().expect("DNA marker part is not an array")
                .iter()
                .map(|p|p.as_array().expect("DNA marker subpart not an array"))
                .map(|p|{
                    DNAmarkerPart {
                        length: p[0].as_f64().expect("Not an f64") ,
                        strength: match p.get(1) {
                            Some(s) => s.as_u64(),
                            None => None,
                        }
                    }
                })
                .collect();
            ret.insert(name.to_owned(),parts);
        }
        ret
    }

}

#[cfg(test)]
mod tests {
    // use super::*;
    use crate::FACILITY;

    #[test]
    fn test_from_json_file() {
        let marker = FACILITY.dna_markers.get("GeneRuler Mix").unwrap();
        assert_eq!(marker[2].length,8000.0);
        assert_eq!(marker[2].strength,Some(13));

        assert!(FACILITY.dna_iupac['V' as usize]&FACILITY.dna_iupac['G' as usize]>0);
        assert!(FACILITY.dna_iupac['H' as usize]&FACILITY.dna_iupac['G' as usize]==0);
    }
}
