//! Compatibility re-export for the headless allele-aware RNA-read screen.

pub use gentle_engine::allele_hash_screen::*;

/// Bound saved-alignment/base-quality ingress for the shared allele-screen comparison kernel.
pub mod comparison;
pub use comparison::*;
