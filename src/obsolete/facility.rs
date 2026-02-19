// Various constants and DNA-related functions

use crate::dna_ladder::DNALadders;

#[derive(Clone, Debug)]
pub struct Facility {
    pub dna_ladders: DNALadders,
}

impl Default for Facility {
    fn default() -> Self {
        Self::new()
    }
}

impl Facility {
    pub fn new() -> Self {
        Self {
            dna_ladders: DNALadders::default(),
        }
    }

    #[inline(always)]
    pub fn dna_ladders(&self) -> &DNALadders {
        &self.dna_ladders
    }
}
