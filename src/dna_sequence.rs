// use bio::alphabets;
// use bio::data_structures::bwt::{bwt, less, Occ};
// use bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable};
// use bio::data_structures::suffix_array::suffix_array;

use std::{fs::File, collections::HashMap};

use bio::io::fasta;
use gb_io::seq::Seq;

use crate::{restriction_enzyme::{RestrictionEnzyme, RestrictionEnzymeSite}, error::GENtleError, FACILITY};

type DNAstring = Vec<u8>;

#[derive(Clone, Debug, Default)]
pub struct DNAoverhang {
    forward_3: DNAstring,
    forward_5: DNAstring,
    reverse_3: DNAstring,
    reverse_5: DNAstring,
}

impl DNAoverhang {
    pub fn is_blunt(&self) -> bool {
        self.forward_3.is_empty() &&
        self.forward_5.is_empty() &&
        self.reverse_3.is_empty() &&
        self.reverse_5.is_empty()
    }

    pub fn is_sticky(&self) -> bool {
        !self.is_blunt()
    }
}

#[derive(Clone, Debug, Default)]
pub struct Location {
    pub start: i64,
    pub stop: i64,
    pub is_reverse: bool,
}

impl Location {
    pub fn from_genbank_location(location: gb_io::seq::Location) -> Self {
        match location {
            gb_io::seq::Location::Range(from, to) => {
                Self {
                    start: from.0, // ignoring Before
                    stop: to.0, // ignoring After
                    is_reverse: false,
                }
            },
            gb_io::seq::Location::Between(_, _) => todo!(),
            gb_io::seq::Location::Complement(_) => todo!(),
            gb_io::seq::Location::Join(_) => todo!(),
            gb_io::seq::Location::Order(_) => todo!(),
            gb_io::seq::Location::Bond(_) => todo!(),
            gb_io::seq::Location::OneOf(_) => todo!(),
            gb_io::seq::Location::External(_, _) => todo!(),
            gb_io::seq::Location::Gap(_) => todo!(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub enum FeatureKind {
    #[default]
    Source,
    CDS,
    Gene,
    Regulatory,
    ORI,
    ProteinBinding,
    Misc,
    Other(String),
}

impl FeatureKind {
    pub fn from_string(s: String) -> Self {
        match s.to_ascii_uppercase().as_str() {
            "CDS" => Self::CDS,
            "GENE" => Self::Gene,
            "REGULATORY" => Self::Regulatory,
            "REP_ORIGIN" => Self::ORI,
            "SOURCE" => Self::Source,
            "PROTEIN_BIND" => Self::ProteinBinding,
            "MISC_FEATURE" => Self::Misc,
            _ => {
                eprintln!("Unknown feature kind {s}");
                Self::Other(s)
            },
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Feature {
    pub kind: FeatureKind,
    pub location: Location,
    pub kv: HashMap<String,String>,
}

impl Feature {
    pub fn from_genbank(f: gb_io::seq::Feature) -> Self {
        // println!("{:#?}",&f);
        Self {
            kind: FeatureKind::from_string(f.kind.to_string()),
            location: Location::from_genbank_location(f.location),
            kv: f.qualifiers.into_iter().map(|(k,v)|(k.to_string(),v.unwrap_or("".to_string()))).collect(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct DNAsequence {
    name: String,
    description: String,
    forward: DNAstring,
    overhang: DNAoverhang,
    features: Vec<Feature>,
    is_circular: bool,
}

impl DNAsequence {
    pub fn from_fasta_file(filename: &str) -> Result<Vec<DNAsequence>,GENtleError> {
        let file = File::open(filename)?;
        Ok(fasta::Reader::new(file)
            .records()
            .into_iter()
            .filter_map(|record|record.ok())
            .map(|record|DNAsequence::from_fasta_record(&record))
            .collect())
    }

    pub fn from_genbank_file(filename: &str) -> Result<Vec<DNAsequence>,GENtleError> {
        Ok(gb_io::reader::parse_file(filename)?
            .into_iter()
            .map(|seq|DNAsequence::from_genbank_seq(seq))
            .collect())
    }

    pub fn restriction_enzyme_sites(&self, restriction_enzymes: &Vec<RestrictionEnzyme>, max: Option<usize>) -> Vec<RestrictionEnzymeSite> {
        restriction_enzymes.iter()
            .flat_map(|re|re.get_sites(&self, max))
            .collect()
    }

    pub fn from_genbank_seq(seq: Seq) -> Self {
        // Not imported: date, len, molecule_type, division, definition, accession, version, source, dblink, keywords, references, contig
        Self {
            name: seq.name.unwrap_or(format!("Unnamed sequence of {} bp",seq.len.unwrap_or(0))),
            description: seq.comments.join("\n"),
            forward: Self::validate_dna_sequence(&seq.seq),
            overhang: DNAoverhang::default(),
            is_circular: seq.topology==gb_io::seq::Topology::Circular,
            features: seq.features.into_iter().map(|f|Feature::from_genbank(f)).collect()
        }
    }

    pub fn from_fasta_record(record: &bio::io::fasta::Record) -> Self {
        let seq = record.seq();
        let name = record.id().to_string();
        let desc = record.desc().unwrap_or("").to_string();
        Self {
            name,
            description: desc,
            forward: seq.to_owned(),
            overhang: DNAoverhang::default(),
            is_circular: false,
            features: vec![],
        }
    }

    pub fn features(&self) -> &Vec<Feature> {
        &self.features
    }

    pub fn name(&self) -> &String {
        &self.name
    }

    pub fn description(&self) -> &String {
        &self.description
    }

    pub fn get_forward(&self) -> &DNAstring {
        &self.forward
    }

    pub fn to_string(&self) -> String {
        String::from_utf8_lossy(&self.forward).into()
    }

    pub fn get_overhang(&self) -> &DNAoverhang {
        &self.overhang
    }

    pub fn is_circular(&self) -> bool {
        self.is_circular
    }

    pub fn set_circular(&mut self, is_circular: bool) {
        self.is_circular = is_circular;
        // TODO clear overhang if is_circular=true?
    }

    fn validate_dna_sequence(v: &[u8]) -> Vec<u8> {
        v.iter()
            .filter(|c|!c.is_ascii_whitespace())
            .map(|c| c.to_ascii_uppercase())
            .map(|c| if FACILITY.dna_iupac[c as usize]>0 { c } else { b'N' } )
            .collect()
    }
}

impl From<String> for DNAsequence {
    fn from(s: String) -> Self {
        let mut ret = DNAsequence::default();
        let allowed_chars = b"ACGTN"; // TODO IUPAC
        ret.forward = s.to_ascii_uppercase()
            .as_bytes()
            .iter()
            .filter(|c|!c.is_ascii_whitespace())
            .map(|c| if allowed_chars.contains(c) { *c } else { b'N' } )
            .collect();
        ret
    }
}

#[cfg(test)]
mod tests {
    use crate::enzymes::Enzymes;
    use super::*;

    #[test]
    fn test_pgex_3x_fasta() {
        let seq = DNAsequence::from_fasta_file("test_files/pGEX_3X.fa").unwrap();
        let seq = seq.get(0).unwrap();
        assert_eq!(seq.name,"U13852.1");

        let enzymes = Enzymes::from_json_file("assets/enzymes.json").unwrap();
        let all = seq.restriction_enzyme_sites(&enzymes.restriction_enzymes(),None);
        let max3 = seq.restriction_enzyme_sites(&enzymes.restriction_enzymes(),Some(3));
        assert!(all.len()>max3.len());
    }

    #[test]
    fn test_pgex_3x_genbank() {
        let dna = DNAsequence::from_genbank_file("test_files/pGEX-3X.gb").unwrap();
        let dna = dna.get(0).unwrap();
        assert_eq!(dna.name(),"XXU13852");
        assert_eq!(dna.features().len(),12);
    }
}
