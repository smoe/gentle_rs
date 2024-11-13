use crate::{
    gc_contents::GcContents,
    iupac_code::IupacCode,
    methylation_sites::{MethylationMode, MethylationSites},
    open_reading_frame::OpenReadingFrame,
    restriction_enzyme::{RestrictionEnzyme, RestrictionEnzymeKey, RestrictionEnzymeSite},
};
use anyhow::Result;
use bio::io::fasta;
use gb_io::seq::{Feature, Seq, Topology};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_with::serde_as;
use std::{
    collections::HashMap,
    fmt,
    fs::File,
    ops::{Range, RangeInclusive},
};

type DNAstring = Vec<u8>;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct DNAoverhang {
    pub forward_3: DNAstring,
    pub forward_5: DNAstring,
    pub reverse_3: DNAstring,
    pub reverse_5: DNAstring,
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

#[serde_as]
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DNAsequence {
    seq: Seq,
    overhang: DNAoverhang,
    restriction_enzymes: Vec<RestrictionEnzyme>,
    restriction_enzyme_sites: Vec<RestrictionEnzymeSite>,
    #[serde_as(as = "Vec<(_, _)>")]
    restriction_enzyme_groups: HashMap<RestrictionEnzymeKey, Vec<String>>,
    max_restriction_enzyme_sites: Option<usize>,
    open_reading_frames: Vec<OpenReadingFrame>,
    methylation_sites: MethylationSites,
    methylation_mode: MethylationMode,
    gc_content: GcContents,
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

    pub fn write_genbank_file(&self, filename: &str) -> Result<()> {
        let file = File::create(filename)?;
        gb_io::writer::write(file, &self.seq)?;
        Ok(())
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

    #[inline(always)]
    pub fn get_base_safe(&self, i: usize) -> Option<u8> {
        let i = if self.is_circular() {
            i % self.len()
        } else {
            i
        };
        self.forward().get(i).copied()
    }

    #[inline(always)]
    pub fn get_base_or_n(&self, i: usize) -> u8 {
        let i = if self.is_circular() {
            i % self.len()
        } else {
            i
        };
        self.forward().get(i).unwrap_or(&b'N').to_owned()
    }

    pub fn get_inclusive_range_safe(&self, range: RangeInclusive<usize>) -> Option<DNAstring> {
        let start = *range.start();
        let end = *range.end() + 1;
        self.get_range_safe(start..end)
    }

    pub fn get_range_safe(&self, range: Range<usize>) -> Option<DNAstring> {
        let Range { start, end } = range;
        if start >= end {
            return None;
        }
        let start = if self.is_circular() {
            start % self.len()
        } else {
            start
        };
        let end = if self.is_circular() {
            (end - 1) % self.len()
        } else {
            end - 1
        };
        if start >= self.len() || end >= self.len() {
            return None;
        }
        if start > end {
            if self.is_circular() {
                Some(
                    self.forward()[start..]
                        .iter()
                        .chain(self.forward()[..=end].iter())
                        .copied()
                        .collect(),
                )
            } else {
                None
            }
        } else {
            Some(self.forward()[start..=end].to_vec())
        }
    }

    #[inline(always)]
    fn forward(&self) -> &Vec<u8> {
        &self.seq.seq
    }

    #[inline(always)]
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
            restriction_enzyme_groups: HashMap::new(),
            max_restriction_enzyme_sites: Some(3), // TODO default?
            open_reading_frames: vec![],
            methylation_sites: MethylationSites::default(),
            methylation_mode: MethylationMode::default(),
            gc_content: GcContents::default(),
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
            restriction_enzyme_groups: HashMap::new(),
            max_restriction_enzyme_sites: Some(3), // TODO default?
            open_reading_frames: vec![],
            methylation_sites: MethylationSites::default(),
            methylation_mode: MethylationMode::default(), // TODO default?
            gc_content: GcContents::default(),
        }
    }

    pub fn open_reading_frames(&self) -> &Vec<OpenReadingFrame> {
        &self.open_reading_frames
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

    pub fn methylation_mode(&self) -> MethylationMode {
        self.methylation_mode.to_owned()
    }

    pub fn set_methylation_mode(&mut self, mode: MethylationMode) {
        self.methylation_mode = mode;
    }

    pub fn gc_content(&self) -> &GcContents {
        &self.gc_content
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
        let mode = self.methylation_mode.to_owned();
        self.methylation_sites = MethylationSites::new_from_sequence(self.forward(), mode);
    }

    fn update_gc_content(&mut self) {
        self.gc_content = GcContents::new_from_sequence(self.forward());
    }

    pub fn update_computed_features(&mut self) {
        self.update_restriction_enyzme_sites();
        self.update_restriction_enzyme_groups();
        self.update_open_reading_frames();
        self.update_methylation_sites();
        self.update_gc_content();
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
            .map(|c| {
                if IupacCode::is_valid_letter(*c) {
                    c.to_ascii_uppercase()
                } else {
                    b'N'
                }
            })
            .collect()
    }

    pub fn restriction_enzyme_groups(&self) -> &HashMap<RestrictionEnzymeKey, Vec<String>> {
        &self.restriction_enzyme_groups
    }

    /// Draws restriction enzyme sites
    fn update_restriction_enzyme_groups(&mut self) {
        let mut name2cut_count = HashMap::new();
        self.restriction_enzyme_groups = HashMap::new();

        for re_site in self
            .restriction_enzyme_sites
            .iter()
            .filter(|site| site.forward_strand)
        {
            name2cut_count
                .entry(&re_site.enzyme.name)
                .and_modify(|c| *c += 1)
                .or_insert(1);
        }
        for re_site in self
            .restriction_enzyme_sites
            .iter()
            .filter(|site| site.forward_strand)
        {
            let pos = re_site.offset + re_site.enzyme.cut;
            let cut_size = re_site.enzyme.cut;
            let number_of_cuts = name2cut_count.get(&re_site.enzyme.name).unwrap();
            let from = re_site.offset;
            let to = from + re_site.enzyme.sequence.len() as isize;
            let key = RestrictionEnzymeKey::new(pos, cut_size, *number_of_cuts, from, to);
            self.restriction_enzyme_groups
                .entry(key)
                .or_default()
                .push(re_site.enzyme.name.to_owned());
        }
    }

    fn split_at_restriction_enzyme_site_circular(&self, site: &RestrictionEnzymeSite) -> Self {
        let pos1 = site.offset;
        let pos2 = site.offset + site.enzyme.overlap;
        let pos2 = pos2 % self.len() as isize;
        let right = pos1.max(pos2) + 1;

        // Rotate so that position 0 is now the sequence after the cut
        let seq = self.seq.set_origin(right as i64);

        // Cut off the overhanging part of the sequence, and keep it around
        let new_size = self.len() as i64 - site.enzyme.overlap.abs() as i64;
        let overhang = seq.seq[new_size as usize..self.len()].to_owned();
        let overhang_rc: Vec<u8> = overhang
            .iter()
            .map(|c| IupacCode::letter_complement(*c))
            .collect();
        let mut seq = seq.extract_range(0, new_size);

        // Sequence is now linear
        seq.topology = Topology::Linear;

        let mut ret = Self::from_u8(self.forward());
        ret.seq = seq;

        // Add overhangs
        if site.enzyme.overlap > 0 {
            ret.overhang.forward_5 = overhang;
            ret.overhang.reverse_3 = overhang_rc;
        } else {
            // TODO test this
            ret.overhang.forward_3 = overhang;
            ret.overhang.reverse_5 = overhang_rc;
        }

        ret
    }

    fn split_at_restriction_enzyme_site_linear(&self, site: &RestrictionEnzymeSite) -> Vec<Self> {
        let pos1 = site.offset + 1;
        let pos2 = site.offset + site.enzyme.overlap;
        let pos2 = pos2 % self.len() as isize;
        let left = pos1;
        let right = pos1.max(pos2) + 1;

        let overhang = self.seq.seq[pos1 as usize..(pos2 + 1) as usize].to_owned();
        let overhang_rc: Vec<u8> = overhang
            .iter()
            .map(|c| IupacCode::letter_complement(*c))
            .collect();

        let mut seq1 = Self::from_u8(self.forward());
        seq1.seq = self.seq.extract_range(0, left as i64);

        let mut seq2 = Self::from_u8(self.forward());
        seq2.seq = self.seq.extract_range(right as i64, self.len() as i64);

        // Add overhangs
        if site.enzyme.overlap > 0 {
            seq1.overhang.forward_3 = vec![];
            seq1.overhang.reverse_5 = overhang_rc;
            seq2.overhang.forward_5 = overhang;
            seq2.overhang.reverse_3 = vec![];
        } else {
            // TODO test this
            seq1.overhang.forward_3 = overhang;
            seq1.overhang.reverse_5 = vec![];
            seq2.overhang.forward_5 = vec![];
            seq2.overhang.reverse_3 = overhang_rc;
        }

        vec![seq1, seq2]
    }

    pub fn split_at_restriction_enzyme_site(&self, site: &RestrictionEnzymeSite) -> Vec<Self> {
        if self.is_circular() {
            vec![self.split_at_restriction_enzyme_site_circular(site)]
        } else {
            self.split_at_restriction_enzyme_site_linear(site)
        }
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
    fn test_split_at_restriction_enzyme_site_circular() {
        // Create circular test sequence
        let mut orig_seq = DNAsequence::from_sequence("ATGGATCCGC").unwrap();
        orig_seq.seq.topology = Topology::Circular;

        // Load BamHI enzyme
        let bam_hi = Enzymes::default()
            .restriction_enzymes()
            .iter()
            .find(|e| e.name == "BamHI")
            .unwrap()
            .to_owned();

        /*
        ATGGATCCGC => rotate
        CGCATGGATC => remove overlap
        CGCATG

        ATG|GATC CGC
        TAC CTAG|GCG

        GATC CGCATG
             GCGTAC CTAG
        */

        // Get restriction site
        let sites = bam_hi.get_sites(&orig_seq, None);
        assert_eq!(sites.len(), 1);
        let site = sites.first().unwrap();

        // Create new sequence from cut
        let new_seq = orig_seq.split_at_restriction_enzyme_site_circular(site);
        assert_eq!(new_seq.get_forward_string(), "CGCATG");
        assert_eq!(new_seq.overhang.forward_5, "GATC".as_bytes());
        assert_eq!(new_seq.overhang.forward_3, "".as_bytes());
        assert_eq!(new_seq.overhang.reverse_5, "".as_bytes());
        assert_eq!(new_seq.overhang.reverse_3, "CTAG".as_bytes());
    }

    #[test]
    fn test_split_at_restriction_enzyme_site_linear() {
        // Create circular test sequence
        let mut orig_seq = DNAsequence::from_sequence("ATGGATCCGC").unwrap();
        orig_seq.seq.topology = Topology::Linear;

        // Load BamHI enzyme
        let bam_hi = Enzymes::default()
            .restriction_enzymes()
            .iter()
            .find(|e| e.name == "BamHI")
            .unwrap()
            .to_owned();

        // Get restriction site
        let sites = bam_hi.get_sites(&orig_seq, None);
        assert_eq!(sites.len(), 1);
        let site = sites.first().unwrap();

        // Create new sequence from cut
        let seqs = orig_seq.split_at_restriction_enzyme_site_linear(site);
        assert_eq!(seqs[0].get_forward_string(), "ATG");
        assert_eq!(seqs[1].get_forward_string(), "CGC");
        assert_eq!(seqs[0].overhang.forward_3, "".as_bytes());
        assert_eq!(seqs[0].overhang.reverse_5, "CTAG".as_bytes());
        assert_eq!(seqs[1].overhang.forward_5, "GATC".as_bytes());
        assert_eq!(seqs[1].overhang.reverse_3, "".as_bytes());
    }

    #[test]
    fn test_pgex_3x_fasta() {
        let seq = DNAsequence::from_fasta_file("test_files/pGEX_3X.fa").unwrap();
        let seq = seq.first().unwrap();
        assert_eq!(seq.name().clone().unwrap(), "U13852.1");

        let enzymes = Enzymes::default();
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

    #[test]
    fn test_get_base_safe() {
        let mut dna = DNAsequence::from("ATGC".to_string());

        // linear
        dna.set_circular(false);
        assert_eq!(dna.get_base_safe(0), Some(b'A'));
        assert_eq!(dna.get_base_safe(1), Some(b'T'));
        assert_eq!(dna.get_base_safe(2), Some(b'G'));
        assert_eq!(dna.get_base_safe(3), Some(b'C'));
        assert_eq!(dna.get_base_safe(4), None);

        // circular
        dna.set_circular(true);
        assert_eq!(dna.get_base_safe(4), Some(b'A'));
    }

    #[test]
    fn test_get_base_or_n() {
        let mut dna = DNAsequence::from("ATGC".to_string());

        // linear
        dna.set_circular(false);
        assert_eq!(dna.get_base_or_n(0), b'A');
        assert_eq!(dna.get_base_or_n(1), b'T');
        assert_eq!(dna.get_base_or_n(2), b'G');
        assert_eq!(dna.get_base_or_n(3), b'C');
        assert_eq!(dna.get_base_or_n(4), b'N');

        // circular
        dna.set_circular(true);
        assert_eq!(dna.get_base_or_n(4), b'A');
    }

    #[test]
    fn test_get_range_safe() {
        let mut dna = DNAsequence::from("ATGC".to_string());

        // linear
        dna.set_circular(false);
        assert_eq!(dna.get_range_safe(0..4), Some("ATGC".as_bytes().to_vec()));
        assert_eq!(dna.get_range_safe(0..5), None);

        // circular
        dna.set_circular(true);
        assert_eq!(dna.get_range_safe(0..4), Some("ATGC".as_bytes().to_vec()));
        assert_eq!(dna.get_range_safe(4..8), Some("ATGC".as_bytes().to_vec()));
        assert_eq!(dna.get_range_safe(0..5), Some("A".as_bytes().to_vec())); // Converts to 0..1
        assert_eq!(dna.get_range_safe(1..5), Some("TGCA".as_bytes().to_vec())); // Wraps around 0 point
    }

    #[test]
    fn test_get_inclsive_range_safe() {
        let mut dna = DNAsequence::from("ATGC".to_string());

        // linear
        dna.set_circular(false);
        assert_eq!(
            dna.get_inclusive_range_safe(0..=3),
            Some("ATGC".as_bytes().to_vec())
        );
        assert_eq!(dna.get_inclusive_range_safe(0..=4), None);

        // circular
        dna.set_circular(true);
        assert_eq!(
            dna.get_inclusive_range_safe(0..=3),
            Some("ATGC".as_bytes().to_vec())
        );
        assert_eq!(
            dna.get_inclusive_range_safe(4..=7),
            Some("ATGC".as_bytes().to_vec())
        );
        assert_eq!(
            dna.get_inclusive_range_safe(0..=4),
            Some("A".as_bytes().to_vec())
        ); // Converts to 0..1
        assert_eq!(
            dna.get_inclusive_range_safe(1..=4),
            Some("TGCA".as_bytes().to_vec())
        ); // Wraps around 0 point
    }
}
