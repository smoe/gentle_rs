//! Compatibility shim for the IUPAC nucleotide-code helpers.
//!
//! The implementation now lives in `gentle-engine`; this module preserves
//! the existing root-crate import path while the engine crate extraction
//! continues.

pub use gentle_engine::iupac_code::*;
