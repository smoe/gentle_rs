//! Compatibility shim for the headless feature expert-view SVG renderer.
//!
//! The implementation now lives in `gentle-render`; this module preserves the
//! existing root-crate import path while the workspace split is still in
//! progress.

pub use gentle_render::render_feature_expert_svg;
