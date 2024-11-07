use crate::{
    dna_display::DnaDisplay, dna_sequence::DNAsequence, render_sequence::RenderSequence, FACILITY,
};
use eframe::egui::{Align2, Color32, Painter, Pos2, Rect, Vec2};
use std::sync::{Arc, RwLock};

#[derive(Clone, Debug, Default)]
pub struct RowBlank {
    line_height: f32,
}

impl RowBlank {
    pub fn new() -> Self {
        Self::default()
    }

    fn blocks(&self) -> usize {
        0
    }
}

#[derive(Clone, Debug)]
pub struct RowDna {
    dna: Arc<RwLock<DNAsequence>>,
    blocks: usize,
    number_offset: f32,
    line_height: f32,
    char_width: f32,
    block_height: f32,
    block_offset: f32,
    batch_bases: usize,
    bases_per_line: usize,
    show_reverse_complement: bool,
    show_position: bool,
}

impl RowDna {
    pub fn new(dna: Arc<RwLock<DNAsequence>>) -> Self {
        Self {
            dna,
            blocks: 0,
            line_height: 0.0,
            char_width: 0.0,
            number_offset: 0.0,
            block_height: 0.0,
            block_offset: 0.0,
            batch_bases: 10,
            bases_per_line: 0,
            show_reverse_complement: false,
            show_position: true,
        }
    }

    pub fn reverse_complement(mut self) -> Self {
        self.show_reverse_complement = true;
        self.show_position = false;
        self
    }

    fn sequence_position_length(&self) -> usize {
        // TODO more elegant way to get the length of the sequence
        format!("{}", self.seq_len()).len()
    }

    fn seq_len(&self) -> usize {
        self.dna.read().unwrap().len()
    }

    pub fn layout(&mut self, block_offset: f32, block_height: f32, area: &Rect) {
        self.block_height = block_height;
        self.block_offset = block_offset;
        self.number_offset = self.sequence_position_length() as f32 * self.char_width;
        let block_width = area.width() - self.number_offset;
        let batches_per_line =
            (block_width / (self.char_width * (self.batch_bases + 1) as f32)) as usize;
        self.bases_per_line = batches_per_line * self.batch_bases;
        self.blocks = (self.seq_len() + self.bases_per_line - 1) / self.bases_per_line;
    }

    pub fn render(&self, _row_num: usize, block_num: usize, painter: &Painter, rect: &Rect) {
        let seq_offset = block_num * self.bases_per_line;
        let pos = Pos2 {
            x: rect.left() + self.number_offset,
            y: rect.top() + self.block_offset,
        };
        if self.show_position {
            painter.text(
                pos,
                Align2::RIGHT_TOP,
                format!("{}", seq_offset + 1),
                RenderSequence::font(),
                Color32::BLACK,
            );
        }
        let seq_end = (seq_offset + self.bases_per_line).min(self.seq_len());
        if let Some(seq) = self.dna.read().unwrap().forward().get(seq_offset..seq_end) {
            let y = rect.top() + self.block_offset;
            let mut x = pos.x + self.char_width * 2.0;
            seq.iter().enumerate().for_each(|(offset, base)| {
                let mut base = *base as char;
                base.make_ascii_uppercase();
                let base = if self.show_reverse_complement {
                    FACILITY.complement(base)
                } else {
                    base
                };

                painter.text(
                    Pos2 { x, y },
                    Align2::RIGHT_TOP,
                    base,
                    RenderSequence::font(),
                    if self.show_reverse_complement {
                        Color32::DARK_GRAY
                    } else {
                        Color32::BLACK
                    },
                );

                x += self.char_width;
                if (offset + 1) % self.batch_bases == 0 {
                    x += self.char_width;
                }
            });
        }
    }

    fn blocks(&self) -> usize {
        self.blocks
    }
}

#[derive(Clone, Debug)]
pub enum SequenceRow {
    Separator(RowBlank),
    Dna(RowDna),
    Features,
    AminoAcids,
    RestrictionEnzymes,
    Proteases,
}

impl SequenceRow {
    pub fn compute_line_height(&mut self, size: &Vec2) {
        match self {
            Self::Separator(ref mut row) => {
                row.line_height = size.y;
            }
            Self::Dna(ref mut row) => {
                row.char_width = size.x;
                row.line_height = size.y;
            }
            _ => {
                todo!();
            }
        }
    }

    pub fn line_height(&self) -> f32 {
        match self {
            Self::Separator(ref row) => row.line_height,
            Self::Dna(ref row) => row.line_height,
            _ => 0.0,
        }
    }

    pub fn blocks(&self) -> usize {
        match self {
            Self::Separator(ref row) => row.blocks(),
            Self::Dna(ref row) => row.blocks(),
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
            _ => {
                todo!();
            }
        }
    }
}
