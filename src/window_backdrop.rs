//! Compatibility shim for window backdrop configuration and rendering helpers.
//!
//! The implementation now lives in `gentle-gui`; this module preserves the
//! existing root-crate import path while the GUI crate extraction is staged.

pub use gentle_gui::window_backdrop::*;
