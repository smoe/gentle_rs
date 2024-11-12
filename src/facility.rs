// Various constants and DNA-related functions

use crate::{amino_acids::AminoAcids, dna_marker::DNAMarkers};

#[derive(Clone, Debug)]
pub struct Facility {
    pub amino_acids: AminoAcids,
    pub dna_markers: DNAMarkers,
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
            dna_markers: DNAMarkers::default(),
        }
    }

    #[inline(always)]
    pub fn dna_markers(&self) -> &DNAMarkers {
        &self.dna_markers
    }

    #[inline(always)]
    pub fn amino_acids(&self) -> &AminoAcids {
        &self.amino_acids
    }
}
