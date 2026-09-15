//! Deterministic GENtle biology engine — first slice.
//!
//! Currently owns the IUPAC nucleotide-code conversion and validation
//! helpers. Further engine modules will move here in subsequent PRs.

pub mod iupac_code;
/// Deterministic source-coherent structure grouping and binding validation.
pub mod transcript_presentation;
/// Pure identity, panel-policy and coordinate checks for portable TSS profiles.
pub mod tss_profiles;
/// Stateless fixed-window geometry shared by engine operations and reporter envelopes.
pub mod tss_window_geometry;
