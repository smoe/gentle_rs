use amino_acids::AminoAcids;
use dna_marker::DNAMarkers;
use enzymes::Enzymes;
use lazy_static::lazy_static;
use translations::Translations;

pub mod amino_acids;
pub mod app;
pub mod dna_display;
pub mod dna_marker;
pub mod dna_sequence;
pub mod engine;
pub mod enzymes;
pub mod gc_contents;
pub mod genomes;
pub mod icons;
pub mod iupac_code;
pub mod js_interface;
pub mod lua_interface;
pub mod main_area_dna;
pub mod methylation_sites;
pub mod open_reading_frame;
pub mod protease;
pub mod pssm;
pub mod render_dna;
pub mod render_dna_circular;
pub mod render_dna_linear;
pub mod render_export;
pub mod render_sequence;
pub mod restriction_enzyme;
pub mod sequence_rows;
pub mod sequence_rows_blank;
pub mod sequence_rows_dna;
pub mod sequence_rows_restriction_enzymes;
pub mod translations;
pub mod window;
pub mod window_dna;

lazy_static! {
    // Interface translations
    pub static ref TRANSLATIONS: Translations = Translations::default();

    // Restriction enzymes and proteases
    pub static ref ENZYMES: Enzymes = Enzymes::default();

    // Amino acids
    pub static ref AMINO_ACIDS : AminoAcids = AminoAcids::default();

    // DNA markers
    pub static ref DNA_MARKERS: DNAMarkers = DNAMarkers::default();
}
