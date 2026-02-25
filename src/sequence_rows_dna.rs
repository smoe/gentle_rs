//! DNA base-row renderer implementation.

use crate::{
    dna_display::DnaDisplay, dna_sequence::DNAsequence, iupac_code::IupacCode,
    render_sequence::RenderSequence,
};
use eframe::egui::{Align2, Color32, Painter, Pos2, Rect, Stroke, Vec2};
use std::sync::{Arc, RwLock};

const MAX_SEQUENCE_PANEL_BP: usize = 100_000;

#[derive(Clone, Debug)]
pub struct RowDna {
    dna: Arc<RwLock<DNAsequence>>,
    display: Arc<RwLock<DnaDisplay>>,
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
    pub fn new(dna: Arc<RwLock<DNAsequence>>, display: Arc<RwLock<DnaDisplay>>) -> Self {
        Self {
            dna,
            display,
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

    #[inline(always)]
    pub fn reverse_complement(mut self) -> Self {
        self.show_reverse_complement = true;
        self.show_position = false;
        self
    }

    #[inline(always)]
    fn sequence_position_length(&self) -> usize {
        // TODO more elegant way to get the length of the sequence
        format!("{}", self.seq_len()).len()
    }

    #[inline(always)]
    fn seq_len(&self) -> usize {
        self.dna.read().map(|d| d.len()).unwrap_or(0)
    }

    fn window(&self) -> (usize, usize) {
        let Ok(dna) = self.dna.read() else {
            return (0, 0);
        };
        let seq_len = dna.len();
        if seq_len == 0 {
            return (0, 0);
        }
        if dna.is_circular() {
            return (0, seq_len.min(MAX_SEQUENCE_PANEL_BP));
        }
        let Ok(display) = self.display.read() else {
            return (0, 0);
        };
        let mut span = display.linear_view_span_bp();
        if span == 0 || span > seq_len {
            span = seq_len;
        }
        span = span.min(MAX_SEQUENCE_PANEL_BP);
        let max_start = seq_len.saturating_sub(span);
        let start = display.linear_view_start_bp().min(max_start);
        (start, span)
    }

    pub fn layout(&mut self, block_offset: f32, block_height: f32, area: &Rect) {
        self.block_height = block_height;
        self.block_offset = block_offset;
        self.number_offset = self.sequence_position_length() as f32 * self.char_width;
        let block_width = area.width() - self.number_offset;
        let batches_per_line =
            (block_width / (self.char_width * (self.batch_bases + 1) as f32)) as usize;
        let batches_per_line = batches_per_line.max(1);
        self.bases_per_line = batches_per_line * self.batch_bases;
        let (_, span) = self.window();
        self.blocks = if span == 0 {
            0
        } else {
            (span + self.bases_per_line - 1) / self.bases_per_line
        };
    }

    pub fn render(&self, _row_num: usize, block_num: usize, painter: &Painter, rect: &Rect) {
        let (window_start, window_span) = self.window();
        if window_span == 0 || block_num >= self.blocks {
            return;
        }
        let seq_offset = window_start + block_num * self.bases_per_line;
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
        let selection = self.display.read().ok().and_then(|d| d.selection());
        let window_end_exclusive = window_start.saturating_add(window_span).min(self.seq_len());
        let seq_end_exclusive = (seq_offset + self.bases_per_line).min(window_end_exclusive);
        if seq_end_exclusive <= seq_offset {
            return;
        }
        let seq = self.dna.read().ok().and_then(|dna| {
            dna.get_inclusive_range_safe(seq_offset..=seq_end_exclusive.saturating_sub(1))
        });
        if let Some(seq) = seq {
            let y = rect.top() + self.block_offset;
            let mut x = pos.x + self.char_width * 2.0;
            seq.iter().enumerate().for_each(|(offset, base)| {
                let base = if self.show_reverse_complement {
                    IupacCode::letter_complement(*base)
                } else {
                    base.to_ascii_uppercase()
                } as char;

                // Show selection, if any, in primary sequence only
                if !self.show_reverse_complement {
                    if let Some(selection) = &selection {
                        let position = seq_offset + offset;
                        if selection.contains(position) {
                            painter.rect(
                                Rect::from_min_size(
                                    Pos2 {
                                        x: x - self.char_width,
                                        y,
                                    },
                                    Vec2 {
                                        x: self.char_width,
                                        y: self.line_height,
                                    },
                                ),
                                0.0,
                                Color32::LIGHT_GRAY,
                                Stroke::NONE,
                            );
                        }
                    }
                }

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

    #[inline(always)]
    pub fn blocks(&self) -> usize {
        self.blocks
    }

    pub fn block_for_position(&self, position: usize) -> Option<usize> {
        if self.bases_per_line == 0 || self.blocks == 0 {
            return None;
        }
        let (window_start, window_span) = self.window();
        if window_span == 0 {
            return None;
        }
        let window_end_exclusive = window_start.saturating_add(window_span);
        if position < window_start || position >= window_end_exclusive {
            return None;
        }
        Some((position - window_start) / self.bases_per_line)
    }

    pub fn compute_line_height(&mut self, size: &Vec2) {
        self.char_width = size.x;
        self.line_height = size.y;
    }

    pub fn line_height(&self) -> f32 {
        self.line_height
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence_rows::SequenceRow;

    #[test]
    fn test_row_dna() {
        let dna_display = Arc::new(RwLock::new(DnaDisplay::default()));
        let dna = Arc::new(RwLock::new(DNAsequence::from_sequence("ACGT").unwrap()));
        let row = RowDna::new(dna, dna_display.clone());
        let mut row = SequenceRow::Dna(row);
        row.compute_line_height(&Vec2::new(10.0, 10.0));
        row.layout(
            &dna_display,
            0,
            0.0,
            0.0,
            &Rect::from_min_size(Pos2::ZERO, Vec2::new(500.0, 100.0)),
        );
        match row {
            SequenceRow::Dna(inner) => {
                assert_eq!(inner.blocks, 1);
                assert_eq!(inner.bases_per_line, 40);
                assert_eq!(inner.sequence_position_length(), 1);
                assert_eq!(inner.seq_len(), 4);
            }
            _ => {
                panic!("Expected Dna row");
            }
        }
    }
}
