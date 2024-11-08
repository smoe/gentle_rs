use crate::{
    dna_display::{AminoAcidFrame, DnaDisplay},
    dna_sequence::DNAsequence,
    sequence_rows::*,
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
}

impl RenderSequence {
    pub fn new_single_sequence(
        dna: Arc<RwLock<DNAsequence>>,
        display: Arc<RwLock<DnaDisplay>>,
    ) -> Self {
        let mut rows = Vec::new();
        rows.push(SequenceRow::Dna(RowDna::new(dna.clone())));
        // if display.read().unwrap().show_re() {
        //     rows.push(SequenceRow::RestrictionEnzymes);
        // }
        if display.read().unwrap().show_rc() {
            rows.push(SequenceRow::Dna(
                RowDna::new(dna.clone()).reverse_complement(),
            ));
        }
        if display.read().unwrap().aa_frame() != AminoAcidFrame::None {
            rows.push(SequenceRow::AminoAcids);
        }
        rows.push(SequenceRow::Separator(RowBlank::new()));
        Self {
            area: Rect::NOTHING,
            display,
            rows,
            block_height: 0.0,
            max_blocks: 0,
        }
    }

    pub fn area(&self) -> &Rect {
        &self.area
    }

    fn layout_needs_recomputing(&mut self, ui: &mut egui::Ui) -> bool {
        let mut ret = false;

        // Recompute layout if area has changed
        let new_area = ui.available_rect_before_wrap();
        if self.area != new_area {
            ret = true;
            self.area = new_area;
        }

        // Recompute layout if update flag is set
        ret = ret
            || self
                .display
                .read()
                .unwrap()
                .update_layout()
                .update_map_sequence();

        ret
    }

    fn layout_was_updated(&self) {
        self.display
            .write()
            .unwrap()
            .update_layout_mut()
            .map_sequence_updated();
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
        self.area = ui.available_rect_before_wrap();
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
        self.max_blocks = self.rows.iter().map(|row| row.blocks()).max().unwrap();
    }

    pub fn render(&mut self, ui: &mut egui::Ui) {
        if self.layout_needs_recomputing(ui) {
            self.layout(ui);
            self.layout_was_updated();
        }

        egui::ScrollArea::vertical().show_rows(
            ui,
            self.block_height,
            self.max_blocks,
            |ui, row_range| {
                for block_num in row_range {
                    let size = Vec2 {
                        x: self.area.width(),
                        y: self.block_height,
                    };
                    let (response, painter) = ui.allocate_painter(size, Sense::hover());
                    let rect = response.rect;
                    for (row_num, row) in self.rows.iter_mut().enumerate() {
                        row.render(row_num, block_num, &painter, &rect);
                    }
                }
            },
        );
    }
}
