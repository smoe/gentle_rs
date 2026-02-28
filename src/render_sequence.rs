//! Sequence export helpers and render-side formatting.

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
        if self.layout_needs_recomputing(ui) {
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
