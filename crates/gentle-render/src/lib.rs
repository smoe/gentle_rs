//! Headless GENtle export and rendering helpers.
//!
//! This crate owns the reusable non-egui render/export paths that should stay
//! adapter-neutral: feature-expert SVG, lineage SVG, protocol cartoons, and
//! the shared gel-figure renderers used by workflow-driven demos.

mod feature_expert;
pub mod pool_gel;
pub mod protein_gel;
pub mod protocol_cartoon;

pub use feature_expert::{
    SplicingExonTransitionMatrix, compute_splicing_exon_transition_matrix,
    render_feature_expert_svg,
};
