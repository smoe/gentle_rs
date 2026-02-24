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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SyntheticMoleculeType {
    DsDna,
    SsDna,
    Rna,
}

impl SyntheticMoleculeType {
    fn parse(raw: &str) -> Option<Self> {
        let normalized = raw
            .trim()
            .to_ascii_lowercase()
            .replace(['_', '-'], "")
            .replace(' ', "");
        match normalized.as_str() {
            "dsdna" | "dna" | "doublestrandeddna" | "double" | "ds" => Some(Self::DsDna),
            "ssdna" | "singlestrandeddna" | "single" | "ssdnaoligo" => Some(Self::SsDna),
            "rna" | "ssrna" | "singlestrandedrna" | "transcript" | "mrna" | "cdna" => {
                Some(Self::Rna)
            }
            _ => None,
        }
    }

    fn molecule_type_value(self) -> &'static str {
        match self {
            Self::DsDna => "dsDNA",
            Self::SsDna => "ssDNA",
            Self::Rna => "RNA",
        }
    }

    fn supports_overhangs(self) -> bool {
        matches!(self, Self::DsDna)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
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

impl fmt::Display for DNAoverhang {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let left = self.forward_5.len().max(self.reverse_3.len());
        let line1 = format!(
            "{: >left$} - ... - {}",
            String::from_utf8(self.forward_5.to_owned()).unwrap(),
            String::from_utf8(self.forward_3.to_owned()).unwrap(),
            left = left
        );
        let line2 = format!(
            "{: >left$} - ... - {}",
            String::from_utf8(self.reverse_3.to_owned()).unwrap(),
            String::from_utf8(self.reverse_5.to_owned()).unwrap(),
            left = left
        );
        write!(f, "{}\n{}", line1, line2)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SequenceEqualityError {
    BiotypeMismatch {
        left: Option<String>,
        right: Option<String>,
    },
    TopologyMismatch {
        left_circular: bool,
        right_circular: bool,
    },
    OverhangMismatch {
        left: DNAoverhang,
        right: DNAoverhang,
    },
    LengthMismatch {
        left: usize,
        right: usize,
    },
    BaseMismatch {
        zero_based_position: usize,
        left: u8,
        right: u8,
    },
}

impl fmt::Display for SequenceEqualityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SequenceEqualityError::BiotypeMismatch { left, right } => write!(
                f,
                "Molecule/biotype mismatch (left={:?}, right={:?})",
                left, right
            ),
            SequenceEqualityError::TopologyMismatch {
                left_circular,
                right_circular,
            } => write!(
                f,
                "Topology mismatch (left_circular={}, right_circular={})",
                left_circular, right_circular
            ),
            SequenceEqualityError::OverhangMismatch { left, right } => {
                write!(f, "Overhang mismatch (left='{}', right='{}')", left, right)
            }
            SequenceEqualityError::LengthMismatch { left, right } => {
                write!(f, "Length mismatch (left={}, right={})", left, right)
            }
            SequenceEqualityError::BaseMismatch {
                zero_based_position,
                left,
                right,
            } => write!(
                f,
                "Base mismatch at position {} (left='{}', right='{}')",
                zero_based_position, *left as char, *right as char
            ),
        }
    }
}

impl std::error::Error for SequenceEqualityError {}

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

    pub fn from_embl_file(filename: &str) -> Result<Vec<DNAsequence>> {
        let text = std::fs::read_to_string(filename)?;
        let parsed = parse_embl_records(&text)?;
        Ok(parsed.into_iter().map(DNAsequence::from_genbank_seq).collect())
    }

    pub fn write_genbank_file(&self, filename: &str) -> Result<()> {
        let file = File::create(filename)?;
        let mut seq = self.seq.clone();
        for feature in &mut seq.features {
            let location_text = feature.location.to_gb_format();
            feature.location = gb_io::seq::Location::from_gb_format(&location_text).map_err(|e| {
                anyhow::anyhow!(
                    "Could not canonicalize feature location '{location_text}' before GenBank write: {e}"
                )
            })?;
        }
        gb_io::writer::write(file, &seq)?;
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
            ret.seq.comments.push(desc.to_string());
            ret.apply_fasta_header_metadata(desc);
        } else {
            // FASTA records are treated as synthetic dsDNA by default.
            ret.seq.molecule_type = Some(SyntheticMoleculeType::DsDna.molecule_type_value().into());
        }
        ret
    }

    fn strip_metadata_value(raw: &str) -> String {
        raw.trim_matches(|c: char| {
            c.is_ascii_whitespace() || c == '"' || c == '\'' || c == ',' || c == ';'
        })
        .to_string()
    }

    fn parse_fasta_header_metadata(desc: &str) -> HashMap<String, String> {
        let mut out = HashMap::new();
        for raw in desc.split_whitespace() {
            let token = raw.trim_matches(|c: char| c == ',' || c == ';');
            let (key, value) = token
                .split_once('=')
                .or_else(|| token.split_once(':'))
                .unwrap_or(("", ""));
            if key.is_empty() {
                continue;
            }
            let value = Self::strip_metadata_value(value);
            if value.is_empty() {
                continue;
            }
            out.insert(key.to_ascii_lowercase(), value);
        }
        out
    }

    fn metadata_value<'a>(meta: &'a HashMap<String, String>, keys: &[&str]) -> Option<&'a str> {
        for key in keys {
            if let Some(value) = meta.get(&key.to_ascii_lowercase()) {
                return Some(value.as_str());
            }
        }
        None
    }

    fn parse_overhang_value(raw: &str) -> DNAstring {
        let trimmed = raw.trim();
        if trimmed.is_empty()
            || trimmed.eq_ignore_ascii_case("none")
            || trimmed.eq_ignore_ascii_case("blunt")
            || trimmed == "."
            || trimmed == "-"
        {
            return vec![];
        }
        Self::validate_dna_sequence(trimmed.as_bytes())
    }

    fn apply_fasta_header_metadata(&mut self, desc: &str) {
        let meta = Self::parse_fasta_header_metadata(desc);

        let molecule = Self::metadata_value(
            &meta,
            &[
                "molecule",
                "molecule_type",
                "mol",
                "mol_type",
                "type",
                "biotype",
            ],
        )
        .and_then(SyntheticMoleculeType::parse)
        .unwrap_or(SyntheticMoleculeType::DsDna);
        self.seq.molecule_type = Some(molecule.molecule_type_value().to_string());

        if matches!(molecule, SyntheticMoleculeType::Rna) {
            // Normalize RNA imports to U when users provide T in FASTA.
            for nt in &mut self.seq.seq {
                *nt = match nt.to_ascii_uppercase() {
                    b'T' => b'U',
                    other => other,
                };
            }
        }

        if let Some(topology_raw) = Self::metadata_value(&meta, &["topology"]) {
            self.seq.topology = if topology_raw.eq_ignore_ascii_case("circular") {
                Topology::Circular
            } else {
                Topology::Linear
            };
        }

        let mut had_overhang_field = false;
        let mut overhang = DNAoverhang::default();

        if let Some(v) = Self::metadata_value(&meta, &["forward_5", "f5", "oh_f5", "overhang_f5"]) {
            had_overhang_field = true;
            overhang.forward_5 = Self::parse_overhang_value(v);
        }
        if let Some(v) = Self::metadata_value(&meta, &["forward_3", "f3", "oh_f3", "overhang_f3"]) {
            had_overhang_field = true;
            overhang.forward_3 = Self::parse_overhang_value(v);
        }
        if let Some(v) = Self::metadata_value(&meta, &["reverse_5", "r5", "oh_r5", "overhang_r5"]) {
            had_overhang_field = true;
            overhang.reverse_5 = Self::parse_overhang_value(v);
        }
        if let Some(v) = Self::metadata_value(&meta, &["reverse_3", "r3", "oh_r3", "overhang_r3"]) {
            had_overhang_field = true;
            overhang.reverse_3 = Self::parse_overhang_value(v);
        }

        if had_overhang_field {
            if molecule.supports_overhangs() {
                if self.seq.topology == Topology::Circular {
                    self.seq.comments.push(
                        "Overhang metadata requested circular topology; forced linear topology"
                            .to_string(),
                    );
                    self.seq.topology = Topology::Linear;
                }
                self.overhang = overhang;
            } else {
                self.seq.comments.push(
                    "Ignored overhang metadata for non-double-stranded molecule type".to_string(),
                );
            }
        }
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

    pub fn features_mut(&mut self) -> &mut Vec<Feature> {
        &mut self.seq.features
    }

    pub fn name(&self) -> &Option<String> {
        &self.seq.name
    }

    pub fn molecule_type(&self) -> Option<&str> {
        self.seq.molecule_type.as_deref()
    }

    pub fn description(&self) -> &Vec<String> {
        &self.seq.comments
    }

    pub fn get_forward_string(&self) -> String {
        std::str::from_utf8(self.forward()).unwrap().to_string()
    }

    pub fn overhang(&self) -> &DNAoverhang {
        &self.overhang
    }

    pub fn assert_sequence_equality(&self, other: &Self) -> Result<(), SequenceEqualityError> {
        let left_biotype = self.molecule_type().map(ToString::to_string);
        let right_biotype = other.molecule_type().map(ToString::to_string);
        if left_biotype != right_biotype {
            return Err(SequenceEqualityError::BiotypeMismatch {
                left: left_biotype,
                right: right_biotype,
            });
        }

        if self.is_circular() != other.is_circular() {
            return Err(SequenceEqualityError::TopologyMismatch {
                left_circular: self.is_circular(),
                right_circular: other.is_circular(),
            });
        }

        if self.overhang != other.overhang {
            return Err(SequenceEqualityError::OverhangMismatch {
                left: self.overhang.clone(),
                right: other.overhang.clone(),
            });
        }

        let left = self.forward();
        let right = other.forward();

        if left.len() != right.len() {
            return Err(SequenceEqualityError::LengthMismatch {
                left: left.len(),
                right: right.len(),
            });
        }

        if let Some((zero_based_position, (left, right))) = left
            .iter()
            .zip(right.iter())
            .enumerate()
            .find(|(_, (l, r))| l != r)
        {
            return Err(SequenceEqualityError::BaseMismatch {
                zero_based_position,
                left: *left,
                right: *right,
            });
        }

        Ok(())
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
            ret.overhang.reverse_5 = overhang_rc;
        } else {
            ret.overhang.forward_3 = overhang;
            ret.overhang.reverse_3 = overhang_rc;
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
        seq1.overhang = self.overhang.clone();

        let mut seq2 = Self::from_u8(self.forward());
        seq2.seq = self.seq.extract_range(right as i64, self.len() as i64);
        seq2.overhang = self.overhang.clone();

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

    pub fn restriction_enzymes_full_digest(&self, enzymes: Vec<RestrictionEnzyme>) -> Vec<Self> {
        let mut ret = vec![self.to_owned()];
        for enzyme in &enzymes {
            loop {
                let mut found_one = false;
                let mut new_ret = vec![];
                for seq in ret.drain(..) {
                    if let Some(site) = enzyme.get_sites(&seq, None).first() {
                        let tmp = seq.split_at_restriction_enzyme_site(site);
                        new_ret.extend(tmp);
                        found_one = true;
                    } else {
                        new_ret.push(seq);
                    }
                }
                ret = new_ret;
                if !found_one {
                    break;
                }
            }
        }
        ret
    }
}

fn parse_embl_records(text: &str) -> Result<Vec<Seq>> {
    let mut records: Vec<Seq> = vec![];
    let mut current_lines: Vec<String> = vec![];
    for line in text.lines() {
        if line.trim() == "//" {
            if !current_lines.is_empty() {
                records.push(parse_embl_record(&current_lines)?);
                current_lines.clear();
            }
            continue;
        }
        current_lines.push(line.to_string());
    }
    if !current_lines.is_empty() {
        records.push(parse_embl_record(&current_lines)?);
    }
    if records.is_empty() {
        return Err(anyhow::anyhow!("Could not parse EMBL file: no records found"));
    }
    Ok(records)
}

fn parse_embl_record(lines: &[String]) -> Result<Seq> {
    #[derive(Debug, Default)]
    struct PendingEmblFeature {
        kind: Option<gb_io::seq::FeatureKind>,
        location_text: String,
        qualifiers: Vec<(gb_io::seq::QualifierKey, Option<String>)>,
        last_qualifier_index: Option<usize>,
    }

    impl PendingEmblFeature {
        fn push_qualifier_line(&mut self, trimmed: &str) {
            if let Some(raw_qualifier) = trimmed.strip_prefix('/') {
                let (qk, qv) = if let Some((key, value)) = raw_qualifier.split_once('=') {
                    (key.trim(), Some(value.trim().to_string()))
                } else {
                    (raw_qualifier.trim(), None)
                };
                self.qualifiers.push((qk.into(), qv));
                self.last_qualifier_index = Some(self.qualifiers.len().saturating_sub(1));
            } else if let Some(idx) = self.last_qualifier_index {
                if let Some(existing) = self.qualifiers.get_mut(idx).and_then(|(_, v)| v.as_mut()) {
                    if !existing.ends_with(' ') {
                        existing.push(' ');
                    }
                    existing.push_str(trimmed);
                }
            } else {
                self.location_text.push_str(trimmed);
            }
        }

        fn into_feature(self) -> Result<Feature> {
            let kind = self
                .kind
                .ok_or_else(|| anyhow::anyhow!("Encountered empty EMBL feature record"))?;
            let location_text = self.location_text.trim().to_string();
            let parsed_location =
                gb_io::seq::Location::from_gb_format(&location_text).map_err(|e| {
                    anyhow::anyhow!("Could not parse EMBL feature location '{location_text}': {e}")
                })?;
            let location = canonicalize_location(parsed_location).map_err(|e| {
                anyhow::anyhow!(
                    "Could not canonicalize EMBL feature location '{location_text}': {e}"
                )
            })?;
            Ok(Feature {
                kind,
                location,
                qualifiers: self.qualifiers,
            })
        }
    }

    let mut seq = Seq::empty();
    let mut accession: Option<String> = None;
    let mut version: Option<String> = None;
    let mut definition_lines: Vec<String> = vec![];
    let mut sequence_started = false;
    let mut current_feature: Option<PendingEmblFeature> = None;

    for raw_line in lines {
        let line = raw_line.trim_end_matches('\r');
        if sequence_started {
            let chunk: String = line
                .chars()
                .filter(|ch| ch.is_ascii_alphabetic())
                .map(|ch| ch.to_ascii_uppercase())
                .collect();
            if !chunk.is_empty() {
                seq.seq.extend_from_slice(chunk.as_bytes());
            }
            continue;
        }

        if let Some(raw_id) = line.strip_prefix("ID") {
            let id = raw_id.trim();
            if let Some(name) = id.split(';').next().map(str::trim).filter(|v| !v.is_empty()) {
                seq.name = Some(name.to_string());
            }
            let lower = id.to_ascii_lowercase();
            if lower.contains("circular") {
                seq.topology = Topology::Circular;
            } else if lower.contains("linear") {
                seq.topology = Topology::Linear;
            }
            continue;
        }
        if let Some(raw_de) = line.strip_prefix("DE") {
            let value = raw_de.trim();
            if !value.is_empty() {
                definition_lines.push(value.to_string());
            }
            continue;
        }
        if let Some(raw_ac) = line.strip_prefix("AC") {
            if accession.is_none() {
                accession = raw_ac
                    .split(';')
                    .map(str::trim)
                    .find(|part| !part.is_empty())
                    .map(str::to_string);
            }
            continue;
        }
        if let Some(raw_sv) = line.strip_prefix("SV") {
            if version.is_none() {
                version = raw_sv
                    .split(';')
                    .map(str::trim)
                    .find(|part| !part.is_empty())
                    .map(str::to_string);
            }
            continue;
        }
        if line.starts_with("SQ") {
            sequence_started = true;
            continue;
        }
        if !line.starts_with("FT") {
            continue;
        }

        let key = line.get(5..21).unwrap_or_default().trim();
        let value = line.get(21..).unwrap_or_default().trim_end();
        if !key.is_empty() {
            if let Some(feature) = current_feature.take() {
                seq.features.push(feature.into_feature()?);
            }
            current_feature = Some(PendingEmblFeature {
                kind: Some(gb_io::seq::FeatureKind::from(key)),
                location_text: value.trim().to_string(),
                qualifiers: vec![],
                last_qualifier_index: None,
            });
            continue;
        }

        let trimmed = value.trim_start();
        if trimmed.is_empty() {
            continue;
        }
        let Some(feature) = current_feature.as_mut() else {
            return Err(anyhow::anyhow!(
                "Encountered EMBL qualifier before any feature: '{trimmed}'"
            ));
        };
        feature.push_qualifier_line(trimmed);
    }

    if let Some(feature) = current_feature.take() {
        seq.features.push(feature.into_feature()?);
    }
    for feature in &mut seq.features {
        for (_, value) in &mut feature.qualifiers {
            if let Some(raw) = value.as_mut() {
                *raw = normalize_embl_qualifier_value(raw);
            }
        }
    }

    seq.definition = (!definition_lines.is_empty()).then_some(definition_lines.join(" "));
    seq.accession = accession.clone();
    seq.version = match (accession, version) {
        (_, Some(raw)) if raw.contains('.') => Some(raw),
        (Some(acc), Some(raw)) if raw.chars().all(|ch| ch.is_ascii_digit()) => {
            Some(format!("{acc}.{raw}"))
        }
        (_, Some(raw)) => Some(raw),
        _ => None,
    };
    seq.len = Some(seq.seq.len());
    if seq.seq.is_empty() {
        return Err(anyhow::anyhow!(
            "Could not parse EMBL record '{}': missing sequence data",
            seq.name.clone().unwrap_or_else(|| "<unnamed>".to_string())
        ));
    }
    Ok(seq)
}

fn canonicalize_location(location: gb_io::seq::Location) -> Result<gb_io::seq::Location> {
    let value = serde_json::to_value(&location)?;
    Ok(serde_json::from_value(value)?)
}

fn normalize_embl_qualifier_value(raw: &str) -> String {
    let mut value = raw.trim().to_string();
    if value.starts_with('"') && value.ends_with('"') && value.len() >= 2 {
        value = value[1..value.len() - 1].to_string();
    }
    value.replace("\"\"", "\"")
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
    use crate::{app::GENtleApp, enzymes::Enzymes};
    use std::io::Write;
    use tempfile::Builder;

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
        assert_eq!(new_seq.overhang.reverse_5, "CTAG".as_bytes());
        assert_eq!(new_seq.overhang.reverse_3, "".as_bytes());
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
    fn test_restriction_enzymes_full_digest() {
        let seq = GENtleApp::load_from_file("test_files/pGEX-3X.gb").unwrap();
        let enzymes = Enzymes::default();
        let res = enzymes.restriction_enzymes_by_name(&["BamHI", "EcoRI"]);
        assert_eq!(res.len(), 2);
        let seqs = seq.restriction_enzymes_full_digest(res);
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].len(), 6);
        assert_eq!(seqs[1].len(), 4938);
        assert_eq!(seqs[0].overhang.forward_5, "gatc".as_bytes());
        assert_eq!(seqs[0].overhang.reverse_5, "TTAA".as_bytes());
        assert_eq!(seqs[1].overhang.forward_5, "aatt".as_bytes());
        assert_eq!(seqs[1].overhang.reverse_5, "CTAG".as_bytes());
    }

    #[test]
    fn test_pgex_3x_fasta() {
        let seq = DNAsequence::from_fasta_file("test_files/pGEX_3X.fa").unwrap();
        let seq = seq.first().unwrap();
        assert_eq!(seq.name().clone().unwrap(), "U13852.1");
        assert_eq!(seq.molecule_type(), Some("dsDNA"));

        let enzymes = Enzymes::default();
        let all = seq.calculate_restriction_enzyme_sites(enzymes.restriction_enzymes(), None);
        let max3 = seq.calculate_restriction_enzyme_sites(enzymes.restriction_enzymes(), Some(3));
        assert!(all.len() > max3.len());
    }

    #[test]
    fn test_toy_small_cross_format_sequence_parity() {
        let fasta =
            DNAsequence::from_fasta_file("test_files/fixtures/import_parity/toy.small.fa").unwrap();
        let genbank =
            DNAsequence::from_genbank_file("test_files/fixtures/import_parity/toy.small.gb")
                .unwrap();
        let embl =
            DNAsequence::from_embl_file("test_files/fixtures/import_parity/toy.small.embl")
                .unwrap();
        let fasta = fasta.first().unwrap();
        let genbank = genbank.first().unwrap();
        let embl = embl.first().unwrap();
        assert_eq!(fasta.len(), 120);
        assert_eq!(genbank.len(), 120);
        assert_eq!(embl.len(), 120);
        assert_eq!(
            fasta.get_forward_string().to_ascii_uppercase(),
            genbank.get_forward_string().to_ascii_uppercase()
        );
        assert_eq!(
            genbank.get_forward_string().to_ascii_uppercase(),
            embl.get_forward_string().to_ascii_uppercase()
        );
        let gene_names: Vec<String> = genbank
            .features()
            .iter()
            .filter(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
            .filter_map(|feature| {
                feature
                    .qualifier_values("gene".into())
                    .next()
                    .map(|value| value.to_string())
            })
            .collect();
        assert!(gene_names.iter().any(|name| name == "toyA"));
        assert!(gene_names.iter().any(|name| name == "toyB"));
    }

    #[test]
    fn test_toy_multi_embl_parses_multiple_records() {
        let seqs =
            DNAsequence::from_embl_file("test_files/fixtures/import_parity/toy.multi.embl")
                .expect("parse multi-record EMBL fixture");
        assert_eq!(seqs.len(), 2);
        let first = &seqs[0];
        let second = &seqs[1];
        assert_eq!(first.len(), 24);
        assert_eq!(second.len(), 24);

        let first_gene = first
            .features()
            .iter()
            .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
            .expect("first record should contain a gene feature");
        assert_eq!(first_gene.location.to_gb_format(), "join(1..4,9..12)");

        let second_misc = second
            .features()
            .iter()
            .find(|feature| {
                feature
                    .kind
                    .to_string()
                    .eq_ignore_ascii_case("misc_feature")
            })
            .expect("second record should contain a misc_feature");
        let note = second_misc
            .qualifier_values("note".into())
            .next()
            .unwrap_or_default();
        assert!(note.contains("line one"));
        assert!(note.contains("line two"));
    }

    #[test]
    fn test_pgex_3x_genbank() {
        let dna = DNAsequence::from_genbank_file("test_files/pGEX-3X.gb").unwrap();
        let dna = dna.first().unwrap();
        assert_eq!(dna.name().clone().unwrap(), "XXU13852");
        assert_eq!(dna.features().len(), 12);
    }

    #[test]
    fn test_pgex_3x_genbank_embl_parity_sequence_and_feature_locations() {
        let genbank = DNAsequence::from_genbank_file("test_files/pGEX-3X.gb").unwrap();
        let embl = DNAsequence::from_embl_file("test_files/pGEX-3X.embl").unwrap();
        let genbank = genbank.first().unwrap();
        let embl = embl.first().unwrap();

        assert_eq!(genbank.len(), embl.len());
        assert_eq!(
            genbank.get_forward_string().to_ascii_uppercase(),
            embl.get_forward_string().to_ascii_uppercase()
        );

        let genbank_features: Vec<(String, (i64, i64))> = genbank
            .features()
            .iter()
            .map(|feature| {
                (
                    feature.kind.to_string().to_ascii_lowercase(),
                    feature
                        .location
                        .find_bounds()
                        .expect("GenBank feature location should be bounded"),
                )
            })
            .collect();
        let embl_features: Vec<(String, (i64, i64))> = embl
            .features()
            .iter()
            .map(|feature| {
                (
                    feature.kind.to_string().to_ascii_lowercase(),
                    feature
                        .location
                        .find_bounds()
                        .expect("EMBL feature location should be bounded"),
                )
            })
            .collect();

        let mut genbank_counts: HashMap<(String, (i64, i64)), usize> = HashMap::new();
        for feature in genbank_features {
            *genbank_counts.entry(feature).or_insert(0) += 1;
        }
        for feature in &embl_features {
            let count = genbank_counts
                .get_mut(feature)
                .expect("EMBL feature location should exist in GenBank");
            assert!(*count > 0);
            *count -= 1;
        }

        let mut missing_from_embl: Vec<(String, (i64, i64))> = genbank_counts
            .into_iter()
            .filter_map(|(feature, count)| (count > 0).then_some(feature))
            .collect();
        missing_from_embl.sort();
        assert_eq!(
            missing_from_embl,
            vec![
                ("gene".to_string(), (1289, 2220)),
                ("gene".to_string(), (3300, 4383))
            ]
        );
    }

    #[test]
    fn test_pgex_3x_embl_load_from_file() {
        let dna = GENtleApp::load_from_file("test_files/pGEX-3X.embl").unwrap();
        assert_eq!(dna.len(), 4952);
        assert!(!dna.features().is_empty());
    }

    #[test]
    fn test_embl_wrapped_feature_location_is_parsed() {
        let embl = "\
ID   TEST1; SV 1; linear; genomic DNA; STD; SYN; 40 BP.\n\
XX\n\
AC   TEST1;\n\
XX\n\
DE   Test wrapped location.\n\
XX\n\
FH   Key             Location/Qualifiers\n\
FH\n\
FT   gene            join(1..5,\n\
FT                   10..15)\n\
FT                   /gene=\"g1\"\n\
FT   misc_feature    complement(20..25)\n\
FT                   /note=\"line one\"\n\
FT                   \"line two\"\n\
SQ   Sequence 40 BP; 10 A; 10 C; 10 G; 10 T; 0 other;\n\
     acgtacgtac gtacgtacgt gtacgtacgt gtacgtacgt        40\n\
//\n";
        let parsed = parse_embl_records(embl).expect("parse EMBL");
        assert_eq!(parsed.len(), 1);
        let dna = DNAsequence::from_genbank_seq(parsed.into_iter().next().unwrap());
        assert_eq!(dna.len(), 40);

        let gene = dna
            .features()
            .iter()
            .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
            .expect("gene feature");
        assert_eq!(
            gene.location.find_bounds().expect("gene bounds"),
            (0, 15),
            "joined location should preserve wrapped continuation segment"
        );
        assert_eq!(
            gene.location.to_gb_format(),
            "join(1..5,10..15)",
            "canonical location formatting should preserve both joined intervals"
        );
        assert!(
            gene.qualifier_values("gene".into()).any(|value| value == "g1"),
            "gene qualifier should be parsed"
        );

        let misc = dna
            .features()
            .iter()
            .find(|feature| {
                feature
                    .kind
                    .to_string()
                    .eq_ignore_ascii_case("misc_feature")
            })
            .expect("misc_feature");
        let note = misc
            .qualifier_values("note".into())
            .next()
            .unwrap_or_default();
        assert!(
            note.contains("line one")
                && note.contains("line two")
                && !note.contains("\"\""),
            "wrapped qualifier text should be concatenated and normalized"
        );

        let mut tmp = Builder::new()
            .suffix(".embl")
            .tempfile()
            .expect("temp EMBL file");
        tmp.write_all(embl.as_bytes()).expect("write EMBL fixture");
        let dna_via_loader = GENtleApp::load_from_file(
            tmp.path()
                .to_str()
                .expect("temp EMBL path should be valid UTF-8"),
        )
        .expect("load_from_file should parse EMBL by extension");
        let gene_via_loader = dna_via_loader
            .features()
            .iter()
            .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
            .expect("gene feature through load_from_file");
        assert_eq!(
            gene_via_loader
                .location
                .find_bounds()
                .expect("gene bounds through load_from_file"),
            (0, 15)
        );
        assert_eq!(
            gene_via_loader.location.to_gb_format(),
            "join(1..5,10..15)",
            "load_from_file should keep wrapped EMBL join locations intact"
        );

        let gb_out = Builder::new()
            .suffix(".gb")
            .tempfile()
            .expect("temp GenBank output");
        let gb_out_path = gb_out.path().to_path_buf();
        drop(gb_out);
        dna_via_loader
            .write_genbank_file(
                gb_out_path
                    .to_str()
                    .expect("temp GenBank output path should be valid UTF-8"),
            )
            .expect("write wrapped EMBL feature to GenBank");
        let exported = std::fs::read_to_string(&gb_out_path).expect("read exported GenBank");
        assert!(
            exported.contains("join(1..5,10..15)") || exported.contains("10..15)"),
            "exported GenBank should retain full join location: {exported}"
        );
    }

    #[test]
    fn test_embl_wrapped_single_feature_genbank_export_keeps_join() {
        let embl = "\
ID   TEST1; SV 1; linear; genomic DNA; STD; SYN; 40 BP.\n\
XX\n\
AC   TEST1;\n\
XX\n\
DE   Test wrapped single feature.\n\
XX\n\
FH   Key             Location/Qualifiers\n\
FH\n\
FT   gene            join(1..5,\n\
FT                   10..15)\n\
FT                   /gene=\"g1\"\n\
SQ   Sequence 40 BP; 10 A; 10 C; 10 G; 10 T; 0 other;\n\
     acgtacgtac gtacgtacgt gtacgtacgt gtacgtacgt        40\n\
//\n";
        let mut tmp = Builder::new()
            .suffix(".embl")
            .tempfile()
            .expect("temp EMBL file");
        tmp.write_all(embl.as_bytes()).expect("write EMBL fixture");

        let dna = GENtleApp::load_from_file(
            tmp.path()
                .to_str()
                .expect("temp EMBL path should be valid UTF-8"),
        )
        .expect("load wrapped EMBL");
        let gene = dna
            .features()
            .iter()
            .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
            .expect("gene feature");
        assert_eq!(
            gene.location.to_gb_format(),
            "join(1..5,10..15)",
            "wrapped single feature should canonicalize to full join"
        );

        let gb_out = Builder::new()
            .suffix(".gb")
            .tempfile()
            .expect("temp GenBank output");
        let gb_out_path = gb_out.path().to_path_buf();
        drop(gb_out);
        dna.write_genbank_file(
            gb_out_path
                .to_str()
                .expect("temp GenBank output path should be valid UTF-8"),
        )
        .expect("write GenBank");
        let exported = std::fs::read_to_string(&gb_out_path).expect("read exported GenBank");
        assert!(
            exported.contains("join(1..5,10..15)") || exported.contains("10..15)"),
            "single-feature wrapped join should survive export: {exported}"
        );
    }

    #[test]
    fn test_genbank_regulatory_qualifiers_are_preserved() {
        let dna = DNAsequence::from_genbank_file("test_files/tp73.ncbi.gb").unwrap();
        let dna = dna.first().unwrap();
        let regulatory = dna
            .features()
            .iter()
            .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("regulatory"))
            .expect("expected at least one regulatory feature");
        assert!(
            regulatory
                .qualifier_values("regulatory_class".into())
                .any(|value| !value.trim().is_empty())
        );
        assert!(
            regulatory
                .qualifier_values("function".into())
                .any(|value| value.to_ascii_lowercase().contains("promoter"))
        );
        assert!(
            regulatory
                .qualifier_values("experiment".into())
                .any(|value| value.to_ascii_lowercase().contains("reporter gene assay"))
        );
        assert!(
            regulatory
                .qualifier_values("db_xref".into())
                .any(|value| value.contains("GeneID:"))
        );
    }

    #[test]
    fn test_fasta_header_sets_ssdna_molecule_type() {
        let record = fasta::Record::with_attrs("oligo_ss", Some("molecule=ssdna"), b"ATGCATGC");
        let dna = DNAsequence::from_fasta_record(&record);
        assert_eq!(dna.molecule_type(), Some("ssDNA"));
        assert!(dna.overhang().is_blunt());
    }

    #[test]
    fn test_fasta_header_sets_rna_and_normalizes_t_to_u() {
        let record = fasta::Record::with_attrs("oligo_rna", Some("molecule=rna"), b"AUGTT");
        let dna = DNAsequence::from_fasta_record(&record);
        assert_eq!(dna.molecule_type(), Some("RNA"));
        assert_eq!(dna.get_forward_string(), "AUGUU".to_string());
        assert!(dna.overhang().is_blunt());
    }

    #[test]
    fn test_fasta_header_overhangs_for_double_stranded_oligo() {
        let record = fasta::Record::with_attrs(
            "oligo_ds",
            Some("molecule=dsdna f5=gatc r5=ctag topology=linear"),
            b"ATGCATGC",
        );
        let dna = DNAsequence::from_fasta_record(&record);
        assert_eq!(dna.molecule_type(), Some("dsDNA"));
        assert_eq!(dna.overhang().forward_5, b"GATC".to_vec());
        assert_eq!(dna.overhang().reverse_5, b"CTAG".to_vec());
        assert_eq!(dna.overhang().forward_3, b"".to_vec());
        assert_eq!(dna.overhang().reverse_3, b"".to_vec());
    }

    #[test]
    fn test_sequence_equality_success_for_identical_synthetic_sequence() {
        let record = fasta::Record::with_attrs(
            "oligo_ds",
            Some("molecule=dsdna f5=gatc r5=ctag topology=linear"),
            b"ATGCATGC",
        );
        let left = DNAsequence::from_fasta_record(&record);
        let right = DNAsequence::from_fasta_record(&record);
        assert_eq!(left.assert_sequence_equality(&right), Ok(()));
    }

    #[test]
    fn test_sequence_equality_reports_biotype_mismatch() {
        let dsdna = DNAsequence::from_fasta_record(&fasta::Record::with_attrs(
            "ds",
            Some("molecule=dsdna"),
            b"ATGC",
        ));
        let ssdna = DNAsequence::from_fasta_record(&fasta::Record::with_attrs(
            "ss",
            Some("molecule=ssdna"),
            b"ATGC",
        ));

        let err = dsdna.assert_sequence_equality(&ssdna).unwrap_err();
        assert_eq!(
            err,
            SequenceEqualityError::BiotypeMismatch {
                left: Some("dsDNA".to_string()),
                right: Some("ssDNA".to_string())
            }
        );
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
