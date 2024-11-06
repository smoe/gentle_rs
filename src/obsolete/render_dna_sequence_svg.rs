// THIS FILE IS CURRENTLY NOT IN USE

use std::io::{self, Write};
// use gb_io::seq::Feature;
// use svg::node::element::{path::Data, Path, SVG};

use crate::dna_sequence::DNAsequence;

const SVG_WIDTH: usize = 4000;
const FONT_SIZE_SMALL: i32 = 96;
const CHAR_WIDTH: usize = 60;
const CHAR_HEIGHT: usize = 70;
const CHARS_PER_BLOCK: usize = 10;

#[derive(Clone, Debug)]
pub enum SubRow<'a> {
    DNA(&'a DNAsequence),
    BlockNumber(&'a DNAsequence),
    Blank,
}

impl<'a> SubRow<'a> {
    fn get_left_offset(&self) -> usize {
        match self {
            SubRow::DNA(dna) => (format!("{}", dna.forward().len()).len() + 1) * CHAR_WIDTH,
            _ => 0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct RenderDnaSequenceSVG {
    pub document: svg::Document,
    number_of_sub_rows: usize,
    left_offset: usize,
}

impl RenderDnaSequenceSVG {
    pub fn write<T>(&self, target: T) -> io::Result<()>
    where
        T: Write,
    {
        svg::write(target, &self.document)
    }

    pub fn get_click_location(&self, x: usize, y: usize) -> (usize, usize) {
        let row_absolute = y / CHAR_HEIGHT;
        let row_blocks_before = row_absolute / self.number_of_sub_rows;
        let sub_row = row_absolute - row_blocks_before * self.number_of_sub_rows;
        (0, sub_row)
    }

    fn nominal_width(&self) -> usize {
        SVG_WIDTH - self.left_offset
    }

    fn linear2xy(&self, position: usize, current_sub_row: usize) -> (usize, usize) {
        let blocks_per_row = self.nominal_width() / (CHAR_WIDTH * (CHARS_PER_BLOCK + 1));
        let chars_per_row = CHARS_PER_BLOCK * blocks_per_row;
        let row = position / (chars_per_row);
        let position_in_row = position - row * chars_per_row;
        let blocks_before = position_in_row / CHARS_PER_BLOCK;
        let x = (position_in_row + blocks_before) * CHAR_WIDTH + self.left_offset;
        let y = CHAR_HEIGHT * (row * (self.number_of_sub_rows) + current_sub_row + 1);
        (x, y)
    }

    pub fn from_dna_sequence(dna: &DNAsequence) -> Self {
        let sub_rows = vec![
            // SubRow::BlockNumber(dna),
            SubRow::DNA(dna),
            SubRow::Blank,
        ];

        let mut ret = Self {
            document: svg::Document::new(),
            number_of_sub_rows: sub_rows.len(),
            left_offset: sub_rows
                .iter()
                .map(|sr| sr.get_left_offset())
                .max()
                .unwrap_or_default(),
        };

        let mut document = svg::Document::new();
        let mut max_y = 0;
        for (current_sub_row, sub_row) in sub_rows.iter().enumerate() {
            document = match sub_row {
                SubRow::DNA(dna) => ret.render_dna_row(document, dna, current_sub_row, &mut max_y),
                SubRow::BlockNumber(_dna) => {
                    ret.render_block_number(document, dna, current_sub_row, &mut max_y)
                }
                SubRow::Blank => document,
            };
        }

        document = document.set("viewBox", (0, 0, SVG_WIDTH, max_y));
        ret.document = document;
        ret
    }

    fn render_dna_row(
        &self,
        mut document: svg::Document,
        dna: &DNAsequence,
        current_sub_row: usize,
        max_y: &mut usize,
    ) -> svg::Document {
        let mut last_y = usize::MAX;
        let number_digits = format!("{}", dna.forward().len()).len();
        for (pos, c) in dna.forward().iter().enumerate() {
            let (x, y) = self.linear2xy(pos, current_sub_row);
            let lower_y = y + 70;
            if *max_y < lower_y {
                *max_y = lower_y;
            }
            if y != last_y {
                last_y = y;
                let mut numtext = format!("{}", pos + 1);
                while numtext.len() < number_digits {
                    numtext = format!("0{numtext}");
                }
                let text = svg::node::element::Text::new(numtext)
                    .set("x", 0)
                    .set("y", y)
                    .set("font-size", FONT_SIZE_SMALL)
                    .set("font-family", "Courier")
                    .set("text-anchor", "left")
                    .set("fill", "grey");
                document = document.add(text);
            }
            let text = svg::node::element::Text::new(format!("{}", (*c) as char).to_uppercase())
                .set("x", x)
                .set("y", y)
                // .set("pos", pos)
                .set("font-size", FONT_SIZE_SMALL)
                .set("font-family", "Courier")
                .set("text-anchor", "left")
                .set("fill", "grey");
            document = document.add(text);
        }
        document
    }

    fn render_block_number(
        &self,
        mut document: svg::Document,
        dna: &DNAsequence,
        current_sub_row: usize,
        max_y: &mut usize,
    ) -> svg::Document {
        for (pos, c) in dna.forward().iter().enumerate() {
            if pos % 10 != 0 {
                continue;
            }
            let (x, y) = self.linear2xy(pos, current_sub_row);
            let lower_y = y + 70;
            if *max_y < lower_y {
                *max_y = lower_y;
            }
            let text = svg::node::element::Text::new(format!("{}", pos + 1).to_uppercase())
                .set("x", x)
                .set("y", y)
                .set("pos", pos)
                .set("font-size", FONT_SIZE_SMALL)
                .set("font-family", "Courier")
                .set("text-anchor", "left")
                .set("fill", "grey");
            document = document.add(text);
        }
        document
    }
}
