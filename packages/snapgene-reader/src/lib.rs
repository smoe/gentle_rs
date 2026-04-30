//! Read SnapGene `.dna` files into reusable Rust data structures.
//!
//! The crate exposes a raw, headless parser first and keeps any downstream
//! biological-model adapters behind optional Cargo features.

mod error;
mod model;
mod parser;
mod xml;

#[cfg(feature = "gb-io")]
mod adapter;

pub use error::SnapGeneError;
pub use model::{
    AdditionalSequenceProperties, DnaPacket, FeatureDirectionality, FeatureRecord, FeatureSegment,
    HybridizationParams, NotesRecord, PacketType, PrimerBindingSite, PrimerRecord, PrimerStrand,
    QualifierRecord, QualifierValue, RawPacket, SequenceTopology, SnapGeneFile,
};
pub use parser::{parse_bytes, parse_path, parse_reader};
