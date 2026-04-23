//! Amino-acid lookup tables and codon translation helpers.

use csv::ReaderBuilder;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

use crate::iupac_code::IupacCode;

const DEFAULT_TRANSLATION_TABLE: usize = 1;
pub const UNKNOWN_CODON: char = '?';
pub const STOP_CODON: char = '|';
const PROTEIN_PKA_N_TERM: f64 = 9.69;
const PROTEIN_PKA_C_TERM: f64 = 2.34;
const PROTEIN_PKA_SIDECHAIN_ACIDIC: &[(char, f64)] =
    &[('D', 3.86), ('E', 4.25), ('C', 8.33), ('Y', 10.07)];
const PROTEIN_PKA_SIDECHAIN_BASIC: &[(char, f64)] = &[('H', 6.00), ('K', 10.53), ('R', 12.48)];

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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AminoAcids {
    pub aas: HashMap<char, AminoAcid>,
    pub codon_tables: HashMap<usize, CodonTable>,
}

impl AminoAcids {
    #[inline(always)]
    pub fn is_start_codon(codon: &[u8; 3]) -> bool {
        *codon == [b'A', b'T', b'G']
    }

    #[inline(always)]
    pub fn is_stop_codon(codon: &[u8; 3]) -> bool {
        *codon == [b'T', b'A', b'A'] || *codon == [b'T', b'A', b'G'] || *codon == [b'T', b'G', b'A']
    }

    fn load_codon_catalog(&mut self) {
        let text = include_str!("../assets/codon_catalog.csv");
        let mut rdr = ReaderBuilder::new()
            .has_headers(false)
            .from_reader(text.as_bytes());
        let mut header: Vec<String> = Vec::new();
        for result in rdr.records() {
            let record: csv::StringRecord = result.expect("Bad CSV line");
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
    fn base2bases(letter: u8) -> Vec<u8> {
        IupacCode::from_letter(letter).to_vec()
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

    pub fn aa2codons(&self, aa: char, translation_table: Option<usize>) -> Vec<[u8; 3]> {
        let translation_table = translation_table.unwrap_or(DEFAULT_TRANSLATION_TABLE);
        let bases = [b'T', b'C', b'A', b'G'];
        let mut out = vec![];
        for b1 in bases {
            for b2 in bases {
                for b3 in bases {
                    let codon = [b1, b2, b3];
                    if self.codon2aa(codon, Some(translation_table)) == aa {
                        out.push(codon);
                    }
                }
            }
        }
        out
    }

    pub fn preferred_species_codon(&self, aa: char, species: &str) -> Option<String> {
        let species = species.trim();
        if species.is_empty() {
            return None;
        }
        self.aas
            .get(&aa)?
            .species_codons
            .get(species)
            .map(|codon| codon.to_ascii_uppercase())
            .filter(|codon| codon.len() == 3)
    }

    /// Estimate a protein's isoelectric point from its amino-acid sequence.
    ///
    /// The estimate uses a deterministic Henderson-Hasselbalch charge balance
    /// model with standard terminal and side-chain pKa values. It is intended
    /// for renderer placement and other reproducible UI summaries, not for
    /// wet-lab-grade pI calibration.
    pub fn protein_isoelectric_point(&self, sequence: &str) -> Option<f32> {
        let mut residues: HashMap<char, usize> = HashMap::new();
        let mut residue_count = 0usize;
        for aa in sequence.trim().chars() {
            let aa = aa.to_ascii_uppercase();
            if aa == STOP_CODON {
                continue;
            }
            if self.aas.contains_key(&aa) {
                residue_count += 1;
                *residues.entry(aa).or_insert(0) += 1;
            }
        }
        if residue_count == 0 {
            return None;
        }

        fn positive_fraction(ph: f64, pka: f64) -> f64 {
            1.0 / (1.0 + 10f64.powf(ph - pka))
        }
        fn negative_fraction(ph: f64, pka: f64) -> f64 {
            1.0 / (1.0 + 10f64.powf(pka - ph))
        }
        fn net_charge(residues: &HashMap<char, usize>, ph: f64) -> f64 {
            let mut charge = positive_fraction(ph, PROTEIN_PKA_N_TERM)
                - negative_fraction(ph, PROTEIN_PKA_C_TERM);
            for (aa, pka) in PROTEIN_PKA_SIDECHAIN_BASIC {
                let count = residues.get(aa).copied().unwrap_or(0) as f64;
                charge += count * positive_fraction(ph, *pka);
            }
            for (aa, pka) in PROTEIN_PKA_SIDECHAIN_ACIDIC {
                let count = residues.get(aa).copied().unwrap_or(0) as f64;
                charge -= count * negative_fraction(ph, *pka);
            }
            charge
        }

        let mut lo = 0.0f64;
        let mut hi = 14.0f64;
        for _ in 0..80 {
            let mid = (lo + hi) * 0.5;
            if net_charge(&residues, mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(((lo + hi) * 0.5) as f32)
    }
}

impl Default for AminoAcids {
    fn default() -> Self {
        let mut aas: HashMap<char, AminoAcid> = HashMap::new();
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
            aas.insert(aa.aa, aa);
        }
        let mut ret = Self {
            aas,
            codon_tables: Self::get_codon_tables(),
        };
        ret.load_codon_catalog();
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AMINO_ACIDS;

    #[test]
    fn test_from_json_file() {
        let aas = AminoAcids::default();
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
        let aas = &AMINO_ACIDS;

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

    #[test]
    fn protein_isoelectric_point_distinguishes_basic_and_acidic_sequences() {
        let aas = AminoAcids::default();
        let acidic = aas
            .protein_isoelectric_point("EEEEEEEE")
            .expect("acidic peptide pI");
        let basic = aas
            .protein_isoelectric_point("KKKKKKKK")
            .expect("basic peptide pI");
        assert!(basic > acidic);
        assert!(acidic < 5.0);
        assert!(basic > 9.0);
    }
}
