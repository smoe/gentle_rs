use enzymes::Enzymes;
use facility::Facility;
use lazy_static::lazy_static;
use translations::Translations;

pub mod amino_acids;
pub mod app;
pub mod dna_display;
pub mod dna_sequence;
pub mod enzymes;
pub mod facility;
pub mod gc_contents;
pub mod icons;
pub mod main_area_dna;
pub mod methylation_sites;
pub mod open_reading_frame;
pub mod protease;
pub mod render_dna;
pub mod render_dna_circular;
pub mod render_dna_linear;
pub mod render_sequence;
pub mod restriction_enzyme;
pub mod sequence_rows;
pub mod translations;
pub mod window;
pub mod window_dna;

lazy_static! {
    // IUPAC codes etc.
    pub static ref FACILITY: Facility = Facility::new();

    // Interface translations
    pub static ref TRANSLATIONS: Translations = Translations::new();

    // Restriction enzymes and proteases
    pub static ref ENZYMES: Enzymes = Enzymes::new().unwrap();
}
