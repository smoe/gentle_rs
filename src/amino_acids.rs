use std::{fs::{self, File}, collections::HashMap};
use csv::ReaderBuilder;
use serde::{Deserialize, Serialize};

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
    pub aa: char,
    pub mw: f32,
    pub pi: f32,
    pub tla: String,
    pub atoms: AminoAcidAtoms,
    pub halflife: Vec<isize>,
    pub chou_fasman: Vec<f32>,
    pub hydrophobicity: AminoAcidHydrophobicity,
    #[serde(skip_serializing,default)]
    pub species_codons: HashMap<String,String>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AminoAcids {
    pub aas : HashMap<char,AminoAcid>,
    pub codon_tables: HashMap<usize,CodonTable>,
}

impl AminoAcids {
    pub fn load() -> Self {
        let mut ret = Self::default();
        let data = fs::read_to_string("assets/amino_acids.json").expect("File not found");
        let res: serde_json::Value = serde_json::from_str(&data).expect("Can not parse JSON");
        let arr = res.as_array().expect("JSON is not an array");
        for row in arr {
            let aa: AminoAcid = match serde_json::from_str(&row.to_string()) {
                Ok(aa) => aa,
                Err(e) => {
                    eprintln!("Bad restriction enzyme: {}: {e}",row.to_string());
                    continue
                }
            };
            ret.aas.insert(aa.aa,aa);
        }
        ret.codon_tables = Self::get_codon_tables();
        ret.load_codon_catalog();
        ret
    }

    fn load_codon_catalog(&mut self) {
        let file = File::open("assets/codon_catalog.csv").expect("Codon catalog file not found");
        let mut rdr = ReaderBuilder::new().has_headers(false).from_reader(file);
        let mut header = vec![];
        for result in rdr.records() {
            let record = result.expect("Bad CSV line");
            if header.is_empty() {
                header = record.iter().skip(2).map(|s|s.to_string()).collect();
                continue;
            }
            let mut record = record.iter();
            let letter = record.next().expect("Bad record").chars().next().unwrap();
            let record = record.skip(1); // TLA
            self.aas.get_mut(&letter).expect("Wrong letter").species_codons = record
                .zip(header.iter())
                .map(|(codon,species)|(species.to_string(),codon.to_string()))
                .collect();
        }
    }

    fn get_codon_tables() -> HashMap<usize,CodonTable> {
        let mut ret = HashMap::new();
        let data = fs::read_to_string("assets/codon_tables.json").expect("File not found");
        let res: serde_json::Value = serde_json::from_str(&data).expect("Invalid JSON");
        let arr = res.as_array().expect("JSON not an array");
        for row in arr {
            let ct: CodonTable = match serde_json::from_str(&row.to_string()) {
                Ok(ct) => ct,
                Err(e) => {
                    eprintln!("Invalid codon table: {}: {e}",row.to_string());
                    continue
                }
            };
            ret.insert(ct.id,ct);
        }
        ret
    }

    pub fn get(&self, aa: char) -> Option<&AminoAcid> {
        self.aas.get(&aa)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_json_file() {
        let aas = AminoAcids::load();
        assert_eq!(aas.get('C').unwrap().atoms.S,1);
        assert_eq!(aas.get('F').unwrap().species_codons.get("Chlamydomonas reinhardtii").unwrap(),"TTC");
    }
}