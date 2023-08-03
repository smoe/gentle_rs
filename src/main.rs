use std::io::stdout;

use lazy_static::lazy_static;
use dna_sequence::DNAsequence;
use render_dna::RenderSVG;

use crate::facility::Facility;

pub mod error;
pub mod facility;
pub mod restriction_enzyme;
pub mod protease;
pub mod enzymes;
pub mod dna_sequence;
pub mod amino_acids;
pub mod render_dna;

lazy_static! {
    pub static ref FACILITY: Facility = Facility::new();
}

fn main() {
    let dna = DNAsequence::from_genbank_file("test_files/pGEX-3X.gb").unwrap();
    let dna = dna.get(0).unwrap();
    let r = RenderSVG::from_dna_sequence(dna);
    let _ = svg::write(stdout(),&r.document);
    // println!("{}",dna.get_forward_string());
}
