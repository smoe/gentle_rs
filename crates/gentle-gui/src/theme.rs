//! Shared egui theme tokens for GUI-only presentation.
//!
//! Window-surface colors are decorative and identify workspace types. Science
//! colors encode biological or analysis semantics. Keep those roles separate so
//! a window tint never changes the meaning of a strand, feature, or analysis.

use egui::{self, Color32, Frame, Margin, Stroke};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WindowSurfacePalette {
    pub main: Color32,
    pub sequence: Color32,
    pub splicing: Color32,
    pub pool: Color32,
    pub configuration: Color32,
    pub help: Color32,
}

impl Default for WindowSurfacePalette {
    fn default() -> Self {
        Self {
            main: Color32::from_rgb(158, 108, 66),
            sequence: Color32::from_rgb(176, 116, 72),
            splicing: Color32::from_rgb(113, 135, 86),
            pool: Color32::from_rgb(148, 95, 58),
            configuration: Color32::from_rgb(169, 118, 79),
            help: Color32::from_rgb(141, 92, 56),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SciencePalette {
    pub forward_strand: Color32,
    pub reverse_strand: Color32,
    pub gene: Color32,
    pub cds: Color32,
    pub exon: Color32,
    pub sequence_node: Color32,
    pub pool_node: Color32,
    pub analysis_node: Color32,
    pub arrangement_node: Color32,
    pub macro_node: Color32,
    pub retrieval_node: Color32,
    pub node_group: Color32,
    pub warning: Color32,
    pub success: Color32,
}

impl Default for SciencePalette {
    fn default() -> Self {
        Self {
            forward_strand: Color32::from_rgb(66, 132, 196),
            reverse_strand: Color32::from_rgb(208, 121, 59),
            gene: Color32::from_rgb(86, 145, 105),
            cds: Color32::from_rgb(143, 91, 175),
            exon: Color32::from_rgb(206, 158, 55),
            sequence_node: Color32::from_rgb(90, 140, 210),
            pool_node: Color32::from_rgb(180, 120, 70),
            analysis_node: Color32::from_rgb(118, 136, 166),
            arrangement_node: Color32::from_rgb(108, 154, 122),
            macro_node: Color32::from_rgb(98, 98, 108),
            retrieval_node: Color32::from_rgb(188, 146, 48),
            node_group: Color32::from_rgb(145, 145, 145),
            warning: Color32::from_rgb(193, 108, 28),
            success: Color32::from_rgb(58, 142, 84),
        }
    }
}

pub const CANVAS_INNER_MARGIN_PX: i8 = 8;

pub fn canvas_fill(dark_mode: bool) -> Color32 {
    if dark_mode {
        Color32::from_rgb(22, 27, 34)
    } else {
        Color32::from_rgb(244, 247, 250)
    }
}

pub fn canvas_stroke(dark_mode: bool) -> Color32 {
    if dark_mode {
        Color32::from_rgb(78, 92, 112)
    } else {
        Color32::from_rgb(184, 194, 207)
    }
}

pub fn canvas_frame(dark_mode: bool) -> Frame {
    Frame::NONE
        .fill(canvas_fill(dark_mode))
        .stroke(Stroke::new(1.0_f32, canvas_stroke(dark_mode)))
        .inner_margin(Margin::same(CANVAS_INNER_MARGIN_PX))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn color_distance(left: Color32, right: Color32) -> u32 {
        let dr = left.r() as i32 - right.r() as i32;
        let dg = left.g() as i32 - right.g() as i32;
        let db = left.b() as i32 - right.b() as i32;
        (dr * dr + dg * dg + db * db) as u32
    }

    #[test]
    fn canvas_helpers_are_mode_aware() {
        assert_ne!(canvas_fill(false), canvas_fill(true));
        assert_ne!(canvas_stroke(false), canvas_stroke(true));
        assert!(color_distance(canvas_fill(false), canvas_stroke(false)) > 800);
        assert!(color_distance(canvas_fill(true), canvas_stroke(true)) > 800);
    }

    #[test]
    fn window_and_science_palettes_stay_semantically_separate() {
        let window = WindowSurfacePalette::default();
        let science = SciencePalette::default();

        assert_ne!(window.sequence, science.forward_strand);
        assert_ne!(window.sequence, science.reverse_strand);
        assert_ne!(window.main, science.analysis_node);
        assert_ne!(window.pool, science.pool_node);
    }
}
