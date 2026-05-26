//! Restriction-enzyme sequence-row renderer implementation.

use crate::{
    dna_display::DnaDisplay, dna_sequence::DNAsequence, restriction_enzyme::RestrictionEnzymeKey,
};
use eframe::egui::{Align2, FontFamily, FontId, Painter, Pos2, Rect, Stroke, Vec2};
use std::sync::{Arc, RwLock};

#[derive(Clone, Debug)]
pub struct RowRestrictionEnzymes {
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
}

impl RowRestrictionEnzymes {
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
        }
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
        let sequence_panel_limit_bp = self
            .display
            .read()
            .map(|display| display.sequence_panel_max_text_length_bp())
            .unwrap_or(200_000);
        if sequence_panel_limit_bp == 0 {
            (0, seq_len)
        } else {
            (0, seq_len.min(sequence_panel_limit_bp))
        }
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
            span.div_ceil(self.bases_per_line)
        };
    }

    fn restriction_site_overlaps_block(
        key: &RestrictionEnzymeKey,
        start: usize,
        end_exclusive: usize,
    ) -> bool {
        if start >= end_exclusive || key.to() < 0 {
            return false;
        }
        let start = start as isize;
        let end_exclusive = end_exclusive as isize;
        key.to() >= start && key.from() < end_exclusive
    }

    fn x_for_sequence_position(&self, rect: &Rect, seq_offset: usize, position: usize) -> f32 {
        let offset = position.saturating_sub(seq_offset);
        let batch_bases = self.batch_bases.max(1);
        rect.left()
            + self.number_offset
            + self.char_width * 2.0
            + self.char_width * offset as f32
            + self.char_width * (offset / batch_bases) as f32
    }

    fn draw_recognition_span(
        &self,
        painter: &Painter,
        rect: &Rect,
        seq_offset: usize,
        from: usize,
        to_inclusive: usize,
        y: f32,
        stroke: Stroke,
    ) {
        if to_inclusive < from {
            return;
        }
        let batch_bases = self.batch_bases.max(1);
        let mut segment_start = from;
        while segment_start <= to_inclusive {
            let offset = segment_start.saturating_sub(seq_offset);
            let remaining_in_batch = batch_bases - (offset % batch_bases);
            let segment_end = to_inclusive.min(
                segment_start
                    .saturating_add(remaining_in_batch)
                    .saturating_sub(1),
            );
            let x1 =
                self.x_for_sequence_position(rect, seq_offset, segment_start) - self.char_width;
            let x2 = self.x_for_sequence_position(rect, seq_offset, segment_end);
            painter.line_segment([Pos2::new(x1, y), Pos2::new(x2, y)], stroke);
            if segment_end == usize::MAX {
                break;
            }
            segment_start = segment_end.saturating_add(1);
        }
    }

    pub fn render(&self, _row_num: usize, block_num: usize, painter: &Painter, rect: &Rect) {
        if self.line_height() <= 0.0 || self.char_width <= 0.0 {
            return;
        }
        let (window_start, window_span) = self.window();
        if window_span == 0 || block_num >= self.blocks {
            return;
        }
        let seq_offset = window_start + block_num * self.bases_per_line;
        let pos = Pos2 {
            x: rect.left() + self.number_offset,
            y: rect.top() + self.block_offset,
        };
        let (seq_end_exclusive, visible_groups) = {
            let Ok(dna) = self.dna.read() else {
                return;
            };
            let window_end_exclusive = window_start.saturating_add(window_span).min(dna.len());
            let seq_end_exclusive = (seq_offset + self.bases_per_line).min(window_end_exclusive);
            if seq_end_exclusive <= seq_offset {
                return;
            }
            let visible_groups = dna
                .restriction_enzyme_groups()
                .iter()
                .filter(|(re, _)| {
                    Self::restriction_site_overlaps_block(re, seq_offset, seq_end_exclusive)
                })
                .map(|(re, label)| (re.clone(), label.clone()))
                .collect::<Vec<_>>();
            (seq_end_exclusive, visible_groups)
        };
        let base_y = pos.y;
        for (re, label) in visible_groups {
            let start = seq_offset as isize;
            let end_exclusive = seq_end_exclusive as isize;
            let visible_from = re.from().max(start) as usize;
            let visible_to = re.to().min(end_exclusive.saturating_sub(1)) as usize;
            let y = base_y + re.number_of_cuts() as f32 - 1.0;
            let label_color = DnaDisplay::restriction_enzyme_group_color(re.number_of_cuts());
            let cut_color = DnaDisplay::restriction_enzyme_geometry_color(re.cut_geometry());
            self.draw_recognition_span(
                painter,
                rect,
                seq_offset,
                visible_from,
                visible_to,
                y,
                Stroke::new(1.0, label_color),
            );
            if re.to() >= start && re.to() < end_exclusive {
                let x = self.x_for_sequence_position(rect, seq_offset, re.to() as usize);
                painter.text(
                    Pos2::new(x, y),
                    Align2::RIGHT_TOP,
                    label.join(","),
                    FontId {
                        size: 9.0,
                        family: FontFamily::Proportional,
                    },
                    label_color,
                );
            }
            if re.pos() >= start && re.pos() < end_exclusive {
                let x = self.x_for_sequence_position(rect, seq_offset, re.pos() as usize);
                painter.line_segment(
                    [Pos2::new(x, y - self.line_height / 2.0), Pos2::new(x, y)],
                    Stroke::new(1.0, cut_color),
                );
            }
            if re.mate_pos() >= start && re.mate_pos() < end_exclusive {
                let x = self.x_for_sequence_position(rect, seq_offset, re.mate_pos() as usize);
                painter.line_segment(
                    [Pos2::new(x, y + self.line_height / 2.0), Pos2::new(x, y)],
                    Stroke::new(1.0, cut_color),
                );
            }
        }
    }

    #[inline(always)]
    pub fn blocks(&self) -> usize {
        if self.line_height() <= 0.0 {
            return 0;
        }
        self.blocks
    }

    pub fn compute_line_height(&mut self, size: &Vec2) {
        self.char_width = size.x;
        self.line_height = size.y;
    }

    pub fn line_height(&self) -> f32 {
        if self
            .display
            .read()
            .map(|d| d.show_reverse_complement())
            .unwrap_or(false)
        {
            self.line_height / 2.0
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hidden_restriction_row_reports_no_blocks() {
        let dna = Arc::new(RwLock::new(
            DNAsequence::from_sequence(&"A".repeat(5_000)).expect("sequence"),
        ));
        let display = Arc::new(RwLock::new(DnaDisplay::default()));
        display
            .write()
            .expect("display")
            .set_show_reverse_complement(false);
        let mut row = RowRestrictionEnzymes::new(dna, display);
        row.compute_line_height(&Vec2::new(10.0, 12.0));
        row.layout(
            0.0,
            12.0,
            &Rect::from_min_size(Pos2::ZERO, Vec2::new(500.0, 100.0)),
        );
        assert_eq!(row.line_height(), 0.0);
        assert_eq!(row.blocks(), 0);
    }

    #[test]
    fn restriction_site_overlap_uses_current_text_block() {
        let site = RestrictionEnzymeKey::new(100, 100, 0, 1, 95, 105);
        assert!(RowRestrictionEnzymes::restriction_site_overlaps_block(
            &site, 90, 96
        ));
        assert!(RowRestrictionEnzymes::restriction_site_overlaps_block(
            &site, 100, 110
        ));
        assert!(!RowRestrictionEnzymes::restriction_site_overlaps_block(
            &site, 0, 95
        ));
        assert!(!RowRestrictionEnzymes::restriction_site_overlaps_block(
            &site, 106, 120
        ));
    }
}
