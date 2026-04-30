//! Compatibility shim for protocol-cartoon catalog and deterministic SVG rendering.
//!
//! The implementation now lives in `gentle-render`; this module preserves the
//! existing root-crate import path while the render crate extraction is staged.

pub use gentle_render::protocol_cartoon::*;
