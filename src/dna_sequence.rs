// use bio::alphabets;
// use bio::data_structures::bwt::{bwt, less, Occ};
// use bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable};
// use bio::data_structures::suffix_array::suffix_array;

use std::fs::File;

use bio::io::fasta;

use crate::{restriction_enzyme::{RestrictionEnzyme, RestrictionEnzymeSite}, error::GENtleError};

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
pub struct DNAsequence {
    name: String,
    description: String,
    forward: DNAstring,
    overhang: DNAoverhang,
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

    pub fn restriction_enzyme_sites(&self, restriction_enzymes: &Vec<RestrictionEnzyme>, max: Option<usize>) -> Vec<RestrictionEnzymeSite> {
        restriction_enzymes.iter()
            .flat_map(|re|re.get_sites(&self, max))
            .collect()
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
        }
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
    fn test_pgex_3x() {
        let seq = DNAsequence::from_fasta_file("test_files/pGEX_3X.fa").unwrap();
        let seq = seq.get(0).unwrap();
        assert_eq!(seq.name,"U13852.1");

        let enzymes = Enzymes::from_json_file("assets/enzymes.json").unwrap();
        let all = seq.restriction_enzyme_sites(&enzymes.restriction_enzymes(),None);
        let max3 = seq.restriction_enzyme_sites(&enzymes.restriction_enzymes(),Some(3));
        assert!(all.len()>max3.len());
    }
}
