use csv::ReaderBuilder;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

use crate::FACILITY;

const DEFAULT_TRANSLATION_TABLE: usize = 1;
pub const UNKNOWN_CODON: char = '?';
pub const STOP_CODON: char = '|';

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct CodonTable {
    pub id: usize, // Not sure if needed
    pub sequence: String,
    pub organism: String,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AminoAcidHydrophobicity {
    pub kyle_doolittle: f32,
    pub hopp_woods: f32,
}

#[allow(non_snake_case)]
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AminoAcidAtoms {
    pub C: u8,
    pub H: u8,
    pub N: u8,
    pub O: u8,
    pub S: u8,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AminoAcid {
    pub aa: char,    // single-character amino acid code
    pub mw: f32,     // molecular weight
    pub pi: f32,     // isoelectric point
    pub tla: String, // Three Letter Acronym
    pub atoms: AminoAcidAtoms,
    pub halflife: Vec<isize>,
    pub chou_fasman: Vec<f32>,
    pub hydrophobicity: AminoAcidHydrophobicity,
    #[serde(skip_serializing, default)]
    pub species_codons: HashMap<String, String>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AminoAcids {
    pub aas: HashMap<char, AminoAcid>,
    pub codon_tables: HashMap<usize, CodonTable>,
}

impl AminoAcids {
    pub fn load() -> Self {
        let mut ret = Self::default();
        let data = include_str!("../assets/amino_acids.json");
        let res: serde_json::Value = serde_json::from_str(data).expect("Can not parse JSON");
        let arr = res.as_array().expect("JSON is not an array");
        for row in arr {
            let aa: AminoAcid = match serde_json::from_str(&row.to_string()) {
                Ok(aa) => aa,
                Err(e) => {
                    eprintln!("Bad restriction enzyme: {}: {e}", row);
                    continue;
                }
            };
            ret.aas.insert(aa.aa, aa);
        }
        ret.codon_tables = Self::get_codon_tables();
        ret.load_codon_catalog();
        ret
    }

    fn load_codon_catalog(&mut self) {
        let text = include_str!("../assets/codon_catalog.csv");
        let mut rdr = ReaderBuilder::new()
            .has_headers(false)
            .from_reader(text.as_bytes());
        let mut header = vec![];
        for result in rdr.records() {
            let record = result.expect("Bad CSV line");
            if header.is_empty() {
                header = record.iter().skip(2).map(|s| s.to_string()).collect();
                continue;
            }
            let mut record = record.iter();
            let letter = record.next().expect("Bad record").chars().next().unwrap();
            let record = record.skip(1); // TLA
            self.aas
                .get_mut(&letter)
                .expect("Wrong letter")
                .species_codons = record
                .zip(header.iter())
                .map(|(codon, species)| (species.to_string(), codon.to_string()))
                .collect();
        }
    }

    fn get_codon_tables() -> HashMap<usize, CodonTable> {
        let mut ret = HashMap::new();
        let data = include_str!("../assets/codon_tables.json");
        let res: serde_json::Value = serde_json::from_str(data).expect("Invalid JSON");
        let arr = res.as_array().expect("JSON not an array");
        for row in arr {
            let ct: CodonTable = match serde_json::from_str(&row.to_string()) {
                Ok(ct) => ct,
                Err(e) => {
                    eprintln!("Invalid codon table: {}: {e}", row);
                    continue;
                }
            };
            ret.insert(ct.id, ct);
        }
        ret
    }

    #[inline(always)]
    pub fn get(&self, aa: char) -> Option<&AminoAcid> {
        self.aas.get(&aa)
    }

    /// Translates a codon base to an index for the codon table.
    /// THESE INDICES ARE SPECIFIC FOR THE CODON TABLES AND NOT thE SAME AS IN `FACILITY`!
    /// #[inline(always)]
    fn acgt(c: u8) -> usize {
        match c.to_ascii_uppercase() {
            b'A' => 2,
            b'C' => 1,
            b'G' => 3,
            b'T' => 0,
            b'U' => 0,
            _ => 250, // Out-of-range
        }
    }

    #[inline(always)]
    fn base2bases(base: u8) -> Vec<u8> {
        FACILITY.split_iupac(FACILITY.base2iupac(base))
    }

    pub fn codon2aa(&self, codon: [u8; 3], translation_table: Option<usize>) -> char {
        let translation_table = translation_table.unwrap_or(DEFAULT_TRANSLATION_TABLE);
        let tt = match self.codon_tables.get(&translation_table) {
            Some(tt) => tt,
            None => return UNKNOWN_CODON,
        };

        // This will return an out-of-range error if a codon base is not A, C, G, or T
        let pos = Self::acgt(codon[0]) * 16 + Self::acgt(codon[1]) * 4 + Self::acgt(codon[2]);
        if pos < 64 {
            tt.sequence.chars().nth(pos).unwrap_or(UNKNOWN_CODON)
        } else {
            let result: HashSet<char> = Self::base2bases(codon[0])
                .into_iter()
                .cartesian_product(Self::base2bases(codon[1]))
                .cartesian_product(Self::base2bases(codon[2]))
                .map(|((c1, c2), c3)| Self::acgt(c1) * 16 + Self::acgt(c2) * 4 + Self::acgt(c3))
                .map(|pos| tt.sequence.chars().nth(pos).unwrap_or(UNKNOWN_CODON))
                .collect();
            if result.len() == 1 {
                // Only one possible amino acid
                *result.iter().next().unwrap()
            } else {
                UNKNOWN_CODON // More than one possible amino acid
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::FACILITY;

    use super::*;

    #[test]
    fn test_from_json_file() {
        let aas = AminoAcids::load();
        assert_eq!(aas.get('C').unwrap().atoms.S, 1);
        assert_eq!(
            aas.get('F')
                .unwrap()
                .species_codons
                .get("Chlamydomonas reinhardtii")
                .unwrap(),
            "TTC"
        );
    }

    #[test]
    fn test_translation_tables() {
        let aas = &FACILITY.amino_acids;

        // Using standard table
        assert_eq!(aas.codon2aa([b'G', b'G', b'T'], None), 'G');
        assert_eq!(aas.codon2aa([b'A', b'C', b'C'], None), 'T');

        // Difference from standard table
        assert_eq!(aas.codon2aa([b'T', b'A', b'A'], None), STOP_CODON);
        assert_eq!(aas.codon2aa([b'T', b'A', b'A'], Some(6)), 'Q');

        // SIUPAC codon
        assert_eq!(aas.codon2aa([b'G', b'C', b'N'], None), 'A');
        assert_eq!(aas.codon2aa([b'G', b'A', b'N'], None), UNKNOWN_CODON);
    }
}
