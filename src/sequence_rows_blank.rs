//! Blank-row renderer implementation.
//!
//! This module implements the `RowBlank` renderer for separator rows between
//! content rows in sequence views. It provides visual spacing and separation
//! without any actual content, acting as a simple spacer.
//!
//! ## Core Responsibilities
//!
//! - **Visual Separation**: Creates empty vertical space between content rows
//!   (DNA, restriction enzymes, etc.) to improve readability and visual hierarchy.
//!
//! - **Line Height Calculation**: Computes and stores the row's height based on
//!   the provided viewport dimensions. The height matches the viewport's Y dimension.
//!
//! - **Minimal Footprint**: Intentionally simple with no dependencies on DNA sequence
//!   or display state, as it carries no content to render.
//!
//! ## Design Principles
//!
//! - **Zero Blocks**: Reports `blocks() -> 0` since it has no content blocks to render.
//!   This is used by the orchestrator to skip block-based calculations.
//!
//! - **No State**: Contains only a `line_height` field, with no sequence or display
//!   references, making it lightweight and fast to clone.
//!
//! - **Default Constructible**: Implements `Default` for easy creation without
//!   parameters via `new()` or `Default::default()`.
//!
//! - **Performance**: All methods are marked `#[inline(always)]` as this is a
//!   hot path component that's created frequently.
//!
//! ## Invariants
//!
//! - Always reports 0 blocks (no content)
//! - Line height is always non-negative
//! - Has no side effects when painted (renders nothing)
//!
//! ## Failure Modes
//!
//! - None: This is a trivial component with no failure modes
//!
//! ## Relationship to Other Modules
//!
//! - **`sequence_rows`**: Defines the row enum that includes this renderer
//! - **`render_sequence`**: Orchestrates multiple rows including this one as separators
//!
//! ## Usage
//!
//! Separator rows are typically inserted by `RenderSequence` between content rows:
//! ```ignore
//! rows.push(SequenceRow::Dna(RowDna::new(...)));
//! rows.push(SequenceRow::Separator(RowBlank::new()));  // <- separator
//! rows.push(SequenceRow::RestrictionEnzymes(RowRestrictionEnzymes::new(...)));
//! ```

use eframe::egui::Vec2;

#[derive(Clone, Debug, Default)]
pub struct RowBlank {
    line_height: f32,
}

impl RowBlank {
    #[inline(always)]
    pub fn new() -> Self {
        Self::default()
    }

    #[inline(always)]
    pub fn blocks(&self) -> usize {
        0
    }

    pub fn line_height(&self) -> f32 {
        self.line_height
    }

    pub fn compute_line_height(&mut self, size: &Vec2) {
        self.line_height = size.y;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_row_blank() {
        let row = RowBlank::default();
        assert_eq!(row.blocks(), 0);
    }
}
