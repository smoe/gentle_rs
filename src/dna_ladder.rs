//! Compatibility shim for DNA/RNA ladder catalogs and defaults.
//!
//! The implementation now lives in `gentle-protocol`; this module preserves
//! the existing root-crate import path while the workspace split is staged.

pub use gentle_protocol::dna_ladder::*;
