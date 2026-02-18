use crate::{dna_display::DnaDisplay, dna_sequence::DNAsequence};
use eframe::egui::{Align2, FontFamily, FontId, Painter, Pos2, Rect, Stroke, Vec2};
use std::sync::{Arc, RwLock};

const MAX_SEQUENCE_PANEL_BP: usize = 100_000;

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
        let window_end_exclusive = window_start.saturating_add(window_span).min(self.seq_len());
        let seq_end_exclusive = (seq_offset + self.bases_per_line).min(window_end_exclusive);
        if seq_end_exclusive <= seq_offset {
            return;
        }
        let seq = self.dna.read().ok().and_then(|dna| {
            dna.get_inclusive_range_safe(seq_offset..=seq_end_exclusive.saturating_sub(1))
        });
        if let Some(seq) = seq {
            let groups = self
                .dna
                .read()
                .map(|dna| dna.restriction_enzyme_groups().clone())
                .unwrap_or_default();
            let y = rect.top() + self.block_offset;
            let mut x = pos.x + self.char_width * 2.0;
            seq.iter().enumerate().for_each(|(offset, _base)| {
                groups
                    .iter()
                    .filter(|(re, _label)| re.from() <= (seq_offset + offset) as isize)
                    .filter(|(re, _label)| re.to() >= (seq_offset + offset) as isize)
                    .for_each(|(re, label)| {
                        let pos = (seq_offset + offset) as isize;
                        let y = y + re.number_of_cuts() as f32 - 1.0;
                        let color = DnaDisplay::restriction_enzyme_group_color(re.number_of_cuts());
                        painter.line_segment(
                            [Pos2::new(x - self.char_width, y), Pos2::new(x, y)],
                            Stroke::new(1.0, color.to_owned()),
                        );
                        if pos == re.to() {
                            painter.text(
                                Pos2::new(x, y),
                                Align2::RIGHT_TOP,
                                label.join(",").to_owned(),
                                FontId {
                                    size: 9.0,
                                    family: FontFamily::Proportional,
                                },
                                color.to_owned(),
                            );
                        }
                        if pos == re.pos() {
                            painter.line_segment(
                                [Pos2::new(x, y - self.line_height / 2.0), Pos2::new(x, y)],
                                Stroke::new(1.0, color.to_owned()),
                            );
                        }
                        if pos == re.pos() + re.cut_size() {
                            painter.line_segment(
                                [Pos2::new(x, y + self.line_height / 2.0), Pos2::new(x, y)],
                                Stroke::new(1.0, color.to_owned()),
                            );
                        }
                    });
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
