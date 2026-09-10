//! Deterministic GENtle biology engine — first slice.
//!
//! Currently owns the IUPAC nucleotide-code conversion and validation
//! helpers. Further engine modules will move here in subsequent PRs.

pub mod iupac_code;
/// Pure identity, panel-policy and coordinate checks for portable TSS profiles.
pub mod tss_profiles;
