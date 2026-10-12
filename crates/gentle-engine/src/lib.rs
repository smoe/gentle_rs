//! Deterministic GENtle biology engine — first slice.
//!
//! Owns root-independent biological parsing and analysis helpers. The root
//! crate keeps compatibility re-exports while the workspace split proceeds.

/// Allele-aware RNA-read hash screen over transcript-coordinate variants.
pub mod allele_hash_screen;
pub mod iupac_code;
/// Standalone streaming RNA target-rescue screen.
pub mod target_rescue;
/// Deterministic source-coherent structure grouping and binding validation.
pub mod transcript_presentation;
/// Pure identity, panel-policy and coordinate checks for portable TSS profiles.
pub mod tss_profiles;
/// Stateless fixed-window geometry shared by engine operations and reporter envelopes.
pub mod tss_window_geometry;
/// UCSC RepeatMasker resource parsing, indexing, and overlap projection.
pub mod ucsc_rmsk;
/// UniProt/SWISS-PROT parsing and transcript/genome projection contracts.
pub mod uniprot;
