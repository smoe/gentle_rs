use lazy_static::lazy_static;
use dna_sequence::DNAsequence;

use crate::facility::Facility;

pub mod error;
pub mod facility;
pub mod restriction_enzyme;
pub mod protease;
pub mod enzymes;
pub mod dna_sequence;
pub mod amino_acids;

lazy_static! {
    pub static ref FACILITY: Facility = Facility::new();
}

fn main() {
    let dna = DNAsequence::from_genbank_file("test_files/pGEX-3X.gb").unwrap();
    let dna = dna.get(0).unwrap();
    println!("{}",dna.get_forward_string());
}
