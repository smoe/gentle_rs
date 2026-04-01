//! Headless GENtle export and rendering helpers.
//!
//! This crate now owns the feature-expert SVG renderer and related shared
//! helper logic, and remains the intended home for further reusable non-egui
//! render/export paths such as lineage SVG, protocol cartoons, gel exports,
//! and similar adapter-neutral figure generation.

mod feature_expert;
pub mod pool_gel;
pub mod protocol_cartoon;

pub use feature_expert::{
    SplicingExonTransitionMatrix, compute_splicing_exon_transition_matrix,
    render_feature_expert_svg,
};
