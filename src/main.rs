#![allow(non_snake_case)]
// import the prelude to get access to the `rsx!` macro and the `Scope` and `Element` types
use dioxus::prelude::*;

use std::io::stdout;

use lazy_static::lazy_static;
use dna_sequence::DNAsequence;
use render_dna::RenderCircularMapSVG;

use crate::facility::Facility;

pub mod error;
pub mod facility;
// pub mod location;
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
    dioxus_desktop::launch(App);
}

fn App(cx: Scope) -> Element {
    let dna = DNAsequence::from_genbank_file("test_files/pGEX-3X.gb").unwrap();
    let dna = dna.get(0).unwrap();
    let r = RenderCircularMapSVG::from_dna_sequence(dna);
    if false {
        let _ = svg::write(stdout(),&r.document);
    }
    let svg_str = r.document.to_string();


    cx.render(rsx! {
        div {
            color:"blue",
            "Hello, world!"
        }
        div {
            dangerous_inner_html: "{svg_str}",
        }
    })
}
