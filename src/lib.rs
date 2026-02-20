use amino_acids::AminoAcids;
use dna_ladder::{default_dna_ladders, default_rna_ladders, LadderCatalog};
use enzymes::Enzymes;
use lazy_static::lazy_static;
use translations::Translations;

pub mod about;
pub mod agent_bridge;
pub mod amino_acids;
pub mod app;
pub mod dna_display;
pub mod dna_ladder;
pub mod dna_sequence;
pub mod engine;
pub mod engine_shell;
pub mod enzymes;
pub mod feature_location;
pub mod gc_contents;
pub mod genomes;
pub mod icons;
pub mod iupac_code;
pub mod js_interface;
pub mod lineage_export;
pub mod lua_interface;
pub mod main_area_dna;
pub mod methylation_sites;
pub mod open_reading_frame;
pub mod pool_gel;
pub mod protease;
pub mod pssm;
pub mod render_dna;
pub mod render_dna_circular;
pub mod render_dna_linear;
pub mod render_export;
pub mod render_sequence;
pub mod resource_sync;
pub mod restriction_enzyme;
pub mod rna_structure;
pub mod sequence_rows;
pub mod sequence_rows_blank;
pub mod sequence_rows_dna;
pub mod sequence_rows_restriction_enzymes;
pub mod shell_docs;
pub mod tf_motifs;
pub mod tool_overrides;
pub mod translations;
pub mod window;
pub mod window_dna;
pub mod workflow_examples;

lazy_static! {
    // Interface translations
    pub static ref TRANSLATIONS: Translations = Translations::default();

    // Restriction enzymes and proteases
    pub static ref ENZYMES: Enzymes = Enzymes::default();

    // Amino acids
    pub static ref AMINO_ACIDS : AminoAcids = AminoAcids::default();

    // DNA ladders
    pub static ref DNA_LADDERS: LadderCatalog = default_dna_ladders();

    // RNA ladders
    pub static ref RNA_LADDERS: LadderCatalog = default_rna_ladders();
}
