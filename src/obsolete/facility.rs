// Various constants and DNA-related functions

use crate::dna_marker::DNAMarkers;

#[derive(Clone, Debug)]
pub struct Facility {
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
            dna_markers: DNAMarkers::default(),
        }
    }

    #[inline(always)]
    pub fn dna_markers(&self) -> &DNAMarkers {
        &self.dna_markers
    }
}
