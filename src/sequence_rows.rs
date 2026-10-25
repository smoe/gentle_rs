//! Sequence-row abstraction shared by specialized row renderers.
//!
//! This module defines the `SequenceRow` enum that serves as the common type for all
//! sequence row implementations. It provides the abstraction layer between the
//! `RenderSequence` orchestrator and concrete row renderers, enabling polymorphic
//! row composition with a unified interface.
//!
//! ## Core Responsibilities
//!
//! - **Row Enum Definition**: Defines all supported row types as variants of
//!   `SequenceRow`, each carrying its specific renderer implementation.
//!
//! - **Uniform Interface**: Provides common methods that all row types implement:
//!   - `compute_line_height()`: Calculates the row's vertical height based on viewport
//!   - `paint()`: Renders the row content within its allocated rectangle
//!   - `get_line_height()`: Returns the current computed height
//!   - `get_row_type_name()`: Returns a human-readable type identifier
//!
//! - **Deferred Implementations**: Includes placeholder variants and methods for
//!   features that are explicitly deferred in the architecture (`AminoAcids`, `Features`,
//!   `Proteases`). These are marked with architecture contract references.
//!
//! - **Optimization Support**: The `#[inline(always)]` attribute on `paint()` ensures
//!   hot rendering paths have minimal call overhead.
//!
//! ## Design Principles
//!
//! - **Extensible Enum**: New row types can be added as new variants without changing
//!   the core interface. Each variant carries its own state.
//!
//! - **Trait-Based Dispatch**: Uses Rust's enum dispatch pattern rather than dynamic
//!   dispatch (trait objects) for maximum performance in the rendering hot path.
//!
//! - **Explicit Deferral**: Deferred row types (AminoAcids, Features, Proteases) have
//!   no-op or minimal implementations with clear architecture references, preventing
//!   accidental implementation before contracts are ready.
//!
//! - **Shared Dependencies**: All row types share common dependencies (painter,
//!   rectangle, display state) to ensure consistent rendering behavior.
//!
//! ## Invariants
//!
//! - Each row variant carries all state needed for its rendering
//! - The `paint()` method never panics (implementations handle edge cases)
//! - Deferred variants have minimal or no-op implementations
//!
//! ## Failure Modes
//!
//! - Deferred row types render as empty (no crash, no content)
//! - Invalid display state defaults to sensible behavior
//! - Zero-size viewport dimensions are handled gracefully
//!
//! ## Relationship to Other Modules
//!
//! - **`render_sequence`**: The orchestrator that composes multiple rows
//! - **`sequence_rows_dna`**: DNA base sequence row implementation
//! - **`sequence_rows_blank`**: Blank separator row implementation
//! - **`sequence_rows_restriction_enzymes`**: Restriction enzyme row implementation
//!
//! ## Architecture Notes
//!
//! Multiple variants are explicitly deferred per architecture contracts:
//! - `AminoAcids`: Deferred pending transcript/CDS-aware translation in engine
//! - `Features`: Deferred pending feature-row visualization contracts
//! - `Proteases`: Deferred pending protease digestion visualization contracts
//!
//! See `docs/architecture.md` for the amino-acid translation row contract details.

use crate::{
    dna_display::DnaDisplay, sequence_rows_blank::RowBlank, sequence_rows_dna::RowDna,
    sequence_rows_restriction_enzymes::RowRestrictionEnzymes,
};
use eframe::egui::{Painter, Rect, Vec2};
use std::sync::{Arc, RwLock};

#[derive(Clone, Debug)]
pub enum SequenceRow {
    Separator(RowBlank),
    Dna(RowDna),
    Features,
    AminoAcids,
    RestrictionEnzymes(RowRestrictionEnzymes),
    Proteases,
}

impl SequenceRow {
    pub fn compute_line_height(&mut self, size: &Vec2) {
        match self {
            Self::Separator(row) => {
                row.compute_line_height(size);
            }
            Self::Dna(row) => {
                row.compute_line_height(size);
            }
            Self::RestrictionEnzymes(row) => {
                row.compute_line_height(size);
            }
            Self::AminoAcids => {
                // Deferred by architecture contract:
                // docs/architecture.md -> "Amino-acid translation row contract (deferred)".
            }
            Self::Features | Self::Proteases => {}
        }
    }

    #[inline(always)]
    pub fn line_height(&self) -> f32 {
        match self {
            Self::Separator(row) => row.line_height(),
            Self::Dna(row) => row.line_height(),
            Self::RestrictionEnzymes(row) => row.line_height(),
            _ => 0.0,
        }
    }

    #[inline(always)]
    pub fn blocks(&self) -> usize {
        match self {
            Self::Separator(row) => row.blocks(),
            Self::Dna(row) => row.blocks(),
            Self::RestrictionEnzymes(row) => row.blocks(),
            _ => 0,
        }
    }

    pub fn layout(
        &mut self,
        _display: &Arc<RwLock<DnaDisplay>>,
        _row_num: usize,
        block_height: f32,
        block_offset: f32,
        area: &Rect,
    ) {
        match self {
            Self::Separator(_row) => {
                // Ignore
            }
            Self::Dna(row) => {
                row.layout(block_offset, block_height, area);
            }
            Self::RestrictionEnzymes(row) => {
                row.layout(block_offset, block_height, area);
            }
            Self::AminoAcids => {
                // Deferred by architecture contract:
                // docs/architecture.md -> "Amino-acid translation row contract (deferred)".
            }
            Self::Features | Self::Proteases => {}
        }
    }

    pub fn render(&self, row_num: usize, block_num: usize, painter: &Painter, rect: &Rect) {
        match self {
            Self::Separator(_row) => {
                // Ignore
            }
            Self::Dna(row) => {
                row.render(row_num, block_num, painter, rect);
            }
            Self::RestrictionEnzymes(row) => {
                row.render(row_num, block_num, painter, rect);
            }
            Self::AminoAcids => {
                // Deferred by architecture contract:
                // docs/architecture.md -> "Amino-acid translation row contract (deferred)".
            }
            Self::Features | Self::Proteases => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::SequenceRow;
    use crate::dna_display::DnaDisplay;
    use eframe::egui::{Pos2, Rect, Vec2};
    use std::sync::{Arc, RwLock};

    #[test]
    fn amino_acid_row_path_is_deferred_noop() {
        let mut row = SequenceRow::AminoAcids;
        row.compute_line_height(&Vec2::new(8.0, 12.0));
        assert_eq!(row.line_height(), 0.0);
        assert_eq!(row.blocks(), 0);
        let display = Arc::new(RwLock::new(DnaDisplay::default()));
        row.layout(
            &display,
            0,
            24.0,
            0.0,
            &Rect::from_min_size(Pos2::new(0.0, 0.0), Vec2::new(240.0, 80.0)),
        );
        assert_eq!(row.line_height(), 0.0);
        assert_eq!(row.blocks(), 0);
    }
}
