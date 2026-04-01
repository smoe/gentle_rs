//! Compatibility shim for embedded GUI icon/resource helpers.
//!
//! The implementation now lives in `gentle-gui`; this module preserves the
//! existing root-crate import path while the GUI crate extraction is staged.

pub use gentle_gui::icons::*;
