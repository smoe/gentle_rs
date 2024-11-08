use crate::{
    methylation_sites::{MethylationMode, MethylationSites},
    open_reading_frame::OpenReadingFrame,
    restriction_enzyme::{RestrictionEnzyme, RestrictionEnzymeSite},
    FACILITY,
};
use anyhow::Result;
use bio::io::fasta;
use gb_io::seq::{Feature, Seq, Topology};
use rayon::prelude::*;
use std::{fmt, fs::File};

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
        self.forward_3.is_empty()
            && self.forward_5.is_empty()
            && self.reverse_3.is_empty()
            && self.reverse_5.is_empty()
    }

    pub fn is_sticky(&self) -> bool {
        !self.is_blunt()
    }
}

#[derive(Clone, Debug)]
pub struct DNAsequence {
    seq: Seq,
    overhang: DNAoverhang,
    restriction_enzymes: Vec<RestrictionEnzyme>,
    restriction_enzyme_sites: Vec<RestrictionEnzymeSite>,
    max_restriction_enzyme_sites: Option<usize>,
    open_reading_frames: Vec<OpenReadingFrame>,
    methylation_sites: MethylationSites,
}

impl DNAsequence {
    pub fn from_sequence(sequence: &str) -> Result<DNAsequence> {
        Ok(DNAsequence::from_u8(sequence.as_bytes()))
    }

    pub fn from_fasta_file(filename: &str) -> Result<Vec<DNAsequence>> {
        let file = File::open(filename)?;
        Ok(fasta::Reader::new(file)
            .records()
            .filter_map(|record| record.ok())
            .map(|record| DNAsequence::from_fasta_record(&record))
            .collect())
    }

    pub fn from_genbank_file(filename: &str) -> Result<Vec<DNAsequence>> {
        Ok(gb_io::reader::parse_file(filename)?
            .into_iter()
            .map(DNAsequence::from_genbank_seq)
            .collect())
    }

    pub fn calculate_restriction_enzyme_sites(
        &self,
        restriction_enzymes: &[RestrictionEnzyme],
        max: Option<usize>,
    ) -> Vec<RestrictionEnzymeSite> {
        restriction_enzymes
            .iter()
            .flat_map(|re| re.get_sites(self, max))
            .collect()
    }

    pub fn forward(&self) -> &Vec<u8> {
        &self.seq.seq
    }

    pub fn len(&self) -> usize {
        self.forward().len()
    }

    pub fn is_empty(&self) -> bool {
        self.forward().is_empty()
    }

    pub fn from_genbank_seq(seq: Seq) -> Self {
        Self {
            seq,
            overhang: DNAoverhang::default(),
            restriction_enzymes: vec![],
            restriction_enzyme_sites: vec![],
            max_restriction_enzyme_sites: Some(3), // TODO default?
            open_reading_frames: vec![],
            methylation_sites: MethylationSites::default(),
        }
    }

    pub fn from_fasta_record(record: &bio::io::fasta::Record) -> Self {
        let seq = record.seq();
        let name = record.id().to_string();
        let mut ret = Self::from_u8(seq);
        ret.seq.name = Some(name);
        if let Some(desc) = record.desc() {
            ret.seq.comments.push(desc.to_string())
        }
        ret
    }

    fn from_u8(s: &[u8]) -> Self {
        let seq = Seq {
            name: None,
            topology: Topology::Linear,
            date: None,
            len: Some(s.len()),
            molecule_type: None, // TODO dna?
            division: String::new(),
            definition: None,
            accession: None,
            version: None,
            source: None,
            dblink: None,
            keywords: None,
            references: vec![],
            comments: vec![],
            seq: s.to_vec(),
            contig: None,
            features: vec![],
        };

        Self {
            seq,
            overhang: DNAoverhang::default(),
            restriction_enzymes: vec![],
            restriction_enzyme_sites: vec![],
            max_restriction_enzyme_sites: Some(3), // TODO default?
            open_reading_frames: vec![],
            methylation_sites: MethylationSites::default(),
        }
    }

    pub fn set_max_restriction_enzyme_sites(&mut self, max_re_sites: Option<usize>) {
        self.max_restriction_enzyme_sites = max_re_sites;
    }

    pub fn restriction_enzymes(&self) -> &Vec<RestrictionEnzyme> {
        &self.restriction_enzymes
    }

    pub fn restriction_enzymes_mut(&mut self) -> &mut Vec<RestrictionEnzyme> {
        &mut self.restriction_enzymes
    }

    pub fn restriction_enzyme_sites(&self) -> &Vec<RestrictionEnzymeSite> {
        &self.restriction_enzyme_sites
    }

    pub fn methylation_sites(&self) -> &MethylationSites {
        &self.methylation_sites
    }

    fn update_restriction_enyzme_sites(&mut self) {
        self.restriction_enzyme_sites = self
            .restriction_enzymes
            .par_iter()
            .flat_map(|re| re.get_sites(self, self.max_restriction_enzyme_sites.to_owned()))
            .collect();
    }

    fn update_open_reading_frames(&mut self) {
        self.open_reading_frames = OpenReadingFrame::find_orfs(self.forward(), self.is_circular());
    }

    fn update_methylation_sites(&mut self) {
        let mode = MethylationMode::default(); // TODO FIXME as a option somewhere
        self.methylation_sites = MethylationSites::new_from_sequence(self.forward(), mode);
    }

    pub fn update_computed_features(&mut self) {
        self.update_restriction_enyzme_sites();
        self.update_open_reading_frames();
        self.update_methylation_sites();
        // TODO amino acids
        // TODO protease sites
    }

    pub fn features(&self) -> &Vec<Feature> {
        &self.seq.features
    }

    pub fn name(&self) -> &Option<String> {
        &self.seq.name
    }

    pub fn description(&self) -> &Vec<String> {
        &self.seq.comments
    }

    pub fn get_forward_string(&self) -> String {
        std::str::from_utf8(self.forward()).unwrap().to_string()
    }

    pub fn get_overhang(&self) -> &DNAoverhang {
        &self.overhang
    }

    pub fn is_circular(&self) -> bool {
        self.seq.topology == Topology::Circular
    }

    pub fn set_circular(&mut self, is_circular: bool) {
        self.seq.topology = match is_circular {
            true => Topology::Circular,
            false => Topology::Linear,
        };
        // TODO clear overhang if is_circular=true?
    }

    pub fn validate_dna_sequence(v: &[u8]) -> Vec<u8> {
        v.iter()
            .filter(|c| !c.is_ascii_whitespace())
            .map(|c| c.to_ascii_uppercase())
            .map(|c| {
                if FACILITY.dna_iupac[c as usize] > 0 {
                    c
                } else {
                    b'N'
                }
            })
            .collect()
    }
}

impl fmt::Display for DNAsequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(self.forward()))
    }
}

impl From<String> for DNAsequence {
    fn from(s: String) -> Self {
        DNAsequence::from_u8(s.as_bytes())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzymes::Enzymes;

    #[test]
    fn test_pgex_3x_fasta() {
        let seq = DNAsequence::from_fasta_file("test_files/pGEX_3X.fa").unwrap();
        let seq = seq.first().unwrap();
        assert_eq!(seq.name().clone().unwrap(), "U13852.1");

        let enzymes = Enzymes::new().unwrap();
        let all = seq.calculate_restriction_enzyme_sites(enzymes.restriction_enzymes(), None);
        let max3 = seq.calculate_restriction_enzyme_sites(enzymes.restriction_enzymes(), Some(3));
        assert!(all.len() > max3.len());
    }

    #[test]
    fn test_pgex_3x_genbank() {
        let dna = DNAsequence::from_genbank_file("test_files/pGEX-3X.gb").unwrap();
        let dna = dna.first().unwrap();
        assert_eq!(dna.name().clone().unwrap(), "XXU13852");
        assert_eq!(dna.features().len(), 12);
    }
}
