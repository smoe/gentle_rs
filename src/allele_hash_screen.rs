//! Compatibility re-export for the headless allele-aware RNA-read screen.

pub use gentle_engine::allele_hash_screen::*;

/// Bound saved-alignment/base-quality ingress for the shared allele-screen comparison kernel.
pub mod comparison;
pub use comparison::*;
/// Authenticated CDS geometry and translation before downstream allele interpretation.
pub mod coding_models;
pub use coding_models::*;
/// File-bound downstream interpretation and pure JSON/text projections.
pub mod interpretation;
