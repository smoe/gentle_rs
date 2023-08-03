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
    let enzymes = enzymes::Enzymes::from_json_file("assets/enzymes.json").unwrap();
    let res = enzymes.restriction_enzymes().into_iter().filter(|re|re.name=="EcoRI").cloned().collect();

    let dna: DNAsequence = "AAAGAATTCTT".to_string().into();
    let all = dna.restriction_enzyme_sites(&res, None);
    let max4 = dna.restriction_enzyme_sites(&res, Some(4));
    println!("{} / {}",all.len(),max4.len());
    // println!("Hello, world!");
}
