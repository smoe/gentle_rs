//! Compatibility shim for virtual pool-gel layout and SVG export.
//!
//! The implementation now lives in `gentle-render`; this module preserves the
//! existing root-crate import path while the render crate extraction is staged.

pub use gentle_render::pool_gel::*;
