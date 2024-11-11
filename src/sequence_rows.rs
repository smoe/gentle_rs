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
            Self::Separator(ref mut row) => {
                row.compute_line_height(size);
            }
            Self::Dna(ref mut row) => {
                row.compute_line_height(size);
            }
            Self::RestrictionEnzymes(ref mut row) => {
                row.compute_line_height(size);
            }
            _ => {
                todo!();
            }
        }
    }

    #[inline(always)]
    pub fn line_height(&self) -> f32 {
        match self {
            Self::Separator(ref row) => row.line_height(),
            Self::Dna(ref row) => row.line_height(),
            Self::RestrictionEnzymes(ref row) => row.line_height(),
            _ => 0.0,
        }
    }

    #[inline(always)]
    pub fn blocks(&self) -> usize {
        match self {
            Self::Separator(ref row) => row.blocks(),
            Self::Dna(ref row) => row.blocks(),
            Self::RestrictionEnzymes(ref row) => row.blocks(),
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
            Self::Dna(ref mut row) => {
                row.layout(block_offset, block_height, area);
            }
            Self::RestrictionEnzymes(ref mut row) => {
                row.layout(block_offset, block_height, area);
            }
            _ => {
                todo!();
            }
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
            _ => {
                todo!();
            }
        }
    }
}
