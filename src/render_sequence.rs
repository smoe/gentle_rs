//! Sequence export helpers and render-side formatting.
//!
//! This module provides the `RenderSequence` component for rendering DNA sequences
//! as text rows in the GUI. It orchestrates multiple sequence row types (DNA bases,
//! restriction enzyme sites, reverse complement, amino acid translations) into a
//! unified, scrollable view.
//!
//! ## Core Responsibilities
//!
//! - **Multi-Row Sequence Display**: Composes multiple `SequenceRow` variants into
//!   a vertical stack, enabling users to view the DNA sequence alongside derived
//!   information (restriction sites, reverse complement, translations).
//!
//! - **Row Type Management**: Creates and manages different row types:
//!   - `RowDna`: The primary DNA base sequence
//!   - `RowRestrictionEnzymes`: Restriction enzyme cut sites aligned to sequence
//!   - Reverse complement row: Shows the complementary strand
//!   - Amino acid translation rows: Deferred pending engine contracts (see architecture)
//!
//! - **Scroll and Selection Tracking**: Manages vertical scrolling through
//!   sequence blocks and tracks the selection anchor for multi-block selections.
//!
//! - **Block-Based Virtualization**: Only renders visible sequence blocks to
//!   maintain performance with long sequences. Blocks are sized consistently
//!   based on the current zoom level.
//!
//! - **Separator Rows**: Inserts blank separator rows between content rows for
//!   visual clarity.
//!
//! ## Design Principles
//!
//! - **Deferred Translation**: Amino acid translation rows are explicitly deferred
//!   in the architecture contract until transcript/CDS-aware translation is
//!   implemented in the engine. See `docs/architecture.md` for details.
//!
//! - **Dynamic Row Composition**: Rows are created based on `DnaDisplay` settings,
//!   allowing users to toggle which derived information is visible without
//!   reloading the sequence.
//!
//! - **Consistent Block Sizing**: All DNA sequence rows use the same block height
//!   and horizontal layout to ensure alignment across rows.
//!
//! - **Thread Safety**: Uses `Arc<RwLock<...>>` for shared display state; each row
//!   holds its own reference to the shared state.
//!
//! ## Invariants
//!
//! - The DNA base row is always present (cannot be toggled off)
//! - Separator rows are inserted between content rows, never at the start
//! - Block dimensions are consistent across all rows at the same zoom level
//! - Selection anchors are cleared when scrolling changes the visible region
//!
//! ## Failure Modes
//!
//! - Empty sequences render as empty rows without errors
//! - Missing display state defaults to sensible visibility settings
//! - Extremely long sequences use virtualization to maintain performance
//!
//! ## Relationship to Other Modules
//!
//! - **`dna_display`**: Provides the `DnaDisplay` configuration for row visibility
//! - **`sequence_rows`**: Defines the `SequenceRow` enum and common row traits
//! - **`sequence_rows_dna`**: Implements the DNA base row renderer
//! - **`sequence_rows_restriction_enzymes`**: Implements the restriction site row
//! - **`sequence_rows_blank`**: Implements the blank separator row
//!
//! ## Architecture Note
//!
//! Amino acid translation is explicitly deferred per:
//! ```text
//! docs/architecture.md -> "Amino-acid translation row contract (deferred)"
//! ```
//! The row creation logic includes a placeholder comment marking this as a
//! future extension point. Do not implement translation rows without first
//! implementing the engine-side transcript/CDS contracts.

use crate::{
    dna_display::{AminoAcidFrame, DnaDisplay},
    dna_sequence::DNAsequence,
    sequence_rows::*,
    sequence_rows_blank::RowBlank,
    sequence_rows_dna::RowDna,
    sequence_rows_restriction_enzymes::RowRestrictionEnzymes,
};
use eframe::egui::{self, Align2, Color32, FontFamily, FontId, Painter, Pos2, Rect, Sense, Vec2};
use std::sync::{Arc, RwLock};

#[derive(Clone, Debug)]
pub struct RenderSequence {
    area: Rect,
    display: Arc<RwLock<DnaDisplay>>,
    rows: Vec<SequenceRow>,
    block_height: f32,
    max_blocks: usize,
    last_selection_anchor: Option<usize>,
    pending_scroll_block: Option<usize>,
}

impl RenderSequence {
    pub fn new_single_sequence(
        dna: Arc<RwLock<DNAsequence>>,
        display: Arc<RwLock<DnaDisplay>>,
    ) -> Self {
        let (show_restriction_enzymes, show_reverse_complement, aa_frame) = display
            .read()
            .map(|d| {
                (
                    d.show_restriction_enzyme_sites(),
                    d.show_reverse_complement(),
                    d.aa_frame(),
                )
            })
            .unwrap_or((false, false, AminoAcidFrame::None));
        let mut rows = Vec::new();
        rows.push(SequenceRow::Dna(RowDna::new(dna.clone(), display.clone())));
        if show_restriction_enzymes {
            rows.push(SequenceRow::RestrictionEnzymes(RowRestrictionEnzymes::new(
                dna.clone(),
                display.clone(),
            )));
        }
        if show_reverse_complement {
            rows.push(SequenceRow::Dna(
                RowDna::new(dna.clone(), display.clone()).reverse_complement(),
            ));
        }
        if aa_frame != AminoAcidFrame::None {
            // Deferred by architecture contract:
            // docs/architecture.md -> "Amino-acid translation row contract (deferred)".
            // Keep sequence-row pipeline deterministic and panic-free until
            // transcript/CDS-aware translation contracts are implemented in engine.
        }
        rows.push(SequenceRow::Separator(RowBlank::new()));
        Self {
            area: Rect::NOTHING,
            display,
            rows,
            block_height: 0.0,
            max_blocks: 0,
            last_selection_anchor: None,
            pending_scroll_block: None,
        }
    }

    pub fn area(&self) -> &Rect {
        &self.area
    }

    fn is_rect_usable(rect: Rect) -> bool {
        rect.min.x.is_finite()
            && rect.min.y.is_finite()
            && rect.max.x.is_finite()
            && rect.max.y.is_finite()
            && rect.width() > 0.0
            && rect.height() > 0.0
    }

    fn layout_needs_recomputing(&mut self, ui: &mut egui::Ui) -> bool {
        let mut ret = false;

        // Recompute layout if area has changed
        let mut new_area = ui.available_rect_before_wrap();
        if !Self::is_rect_usable(new_area) {
            new_area = ui.max_rect();
        }
        if !Self::is_rect_usable(new_area) {
            return false;
        }
        if self.area != new_area {
            ret = true;
            self.area = new_area;
        }

        // Recompute layout if update flag is set
        let display_needs_update = self
            .display
            .read()
            .map(|d| d.update_layout().update_map_sequence())
            .unwrap_or(false);
        ret = ret || display_needs_update;

        ret
    }

    fn layout_was_updated(&self) {
        if let Ok(mut display) = self.display.write() {
            display.update_layout_mut().map_sequence_updated();
        }
    }

    fn block_for_position(&self, position: usize) -> Option<usize> {
        self.rows.iter().find_map(|row| match row {
            SequenceRow::Dna(row) => row.block_for_position(position),
            _ => None,
        })
    }

    fn update_scroll_target_from_selection(&mut self) {
        let selection_anchor = self
            .display
            .read()
            .ok()
            .and_then(|display| display.selection().map(|selection| selection.from()));
        if selection_anchor == self.last_selection_anchor {
            return;
        }
        self.last_selection_anchor = selection_anchor;
        self.pending_scroll_block = selection_anchor.and_then(|pos| self.block_for_position(pos));
    }

    pub fn request_scroll_to_selection(&mut self) {
        self.pending_scroll_block = self
            .display
            .read()
            .ok()
            .and_then(|display| display.selection().map(|selection| selection.from()))
            .and_then(|pos| self.block_for_position(pos));
    }

    pub fn font() -> FontId {
        FontId {
            size: 12.0,
            family: FontFamily::Monospace,
        }
    }

    fn get_size_of_a(painter: &Painter) -> Vec2 {
        let rect = painter.text(
            Pos2 { x: 0.0, y: 0.0 },
            Align2::LEFT_TOP,
            "A",
            Self::font(),
            Color32::WHITE,
        );
        rect.size()
    }

    fn layout(&mut self, ui: &mut egui::Ui) {
        let mut area = ui.available_rect_before_wrap();
        if !Self::is_rect_usable(area) {
            area = ui.max_rect();
        }
        if !Self::is_rect_usable(area) {
            self.area = Rect::NOTHING;
            self.max_blocks = 0;
            return;
        }
        self.area = area;
        let painter = ui.painter();
        let size = Self::get_size_of_a(painter);
        for row in self.rows.iter_mut() {
            row.compute_line_height(&size);
        }
        self.block_height = self.rows.iter().map(|row| row.line_height()).sum();
        let mut block_offset = 0.0;
        for (row_num, row) in self.rows.iter_mut().enumerate() {
            row.layout(
                &self.display,
                row_num,
                self.block_height,
                block_offset,
                &self.area,
            );
            block_offset += row.line_height();
        }
        self.max_blocks = self.rows.iter().map(|row| row.blocks()).max().unwrap_or(0);
    }

    pub fn render(&mut self, ui: &mut egui::Ui) {
        crate::gentle_gui_profile_scope!("RenderSequence::render");
        if self.layout_needs_recomputing(ui) {
            crate::gentle_gui_profile_scope!("RenderSequence::layout");
            self.layout(ui);
            self.layout_was_updated();
        }
        if !Self::is_rect_usable(self.area) || self.max_blocks == 0 || self.block_height <= 0.0 {
            return;
        }

        self.update_scroll_target_from_selection();
        let mut scroll_area = egui::ScrollArea::vertical();
        if let Some(target_block) = self.pending_scroll_block.take() {
            let target_offset = (target_block as f32 * self.block_height).max(0.0);
            scroll_area = scroll_area.vertical_scroll_offset(target_offset);
        }

        scroll_area.show_rows(ui, self.block_height, self.max_blocks, |ui, row_range| {
            for block_num in row_range {
                let size = Vec2::new(
                    self.area.width().clamp(1.0, 100_000.0),
                    self.block_height.clamp(1.0, 100_000.0),
                );
                let (response, painter) = ui.allocate_painter(size, Sense::hover());
                let rect = response.rect;
                for (row_num, row) in self.rows.iter_mut().enumerate() {
                    row.render(row_num, block_num, &painter, &rect);
                }
            }
        });
    }
}
