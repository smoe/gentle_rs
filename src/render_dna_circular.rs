use eframe::egui::{self, Color32, Pos2, Rect, Rounding, Sense, Stroke, Vec2};

use crate::dna_sequence::DNAsequence;
use std::sync::{Arc, Mutex};

pub struct RenderDnaCircular {
    dna: Arc<Mutex<DNAsequence>>,
}

impl RenderDnaCircular {
    pub fn new(dna: Arc<Mutex<DNAsequence>>) -> Self {
        Self { dna }
    }

    pub fn render_dna_circular(&self, ui: &mut egui::Ui) {
        let area = ui.available_rect_before_wrap();

        let _name = self.dna.lock().expect("DNA lock poisoned").name().as_ref();
        let painter = ui.painter();
        // painter.rect_stroke(
        //     Rect {
        //         min: Pos2 {
        //             x: area.left() + 1.0,
        //             y: area.top() + 1.0,
        //         },
        //         max: Pos2 {
        //             x: area.right() - 1.0,
        //             y: area.bottom() - 1.0,
        //         },
        //     },
        //     Rounding::ZERO,
        //     Stroke {
        //         width: 1.0,
        //         color: Color32::RED,
        //     },
        // );

        let max_width = area.width();
        let max_height = area.height();
        let center = Pos2::new(
            (area.left() + area.right()) / 2.0,
            (area.top() + area.bottom()) / 2.0,
        );
        painter.circle_stroke(
            center,
            max_width.min(max_height) * 0.4,
            Stroke {
                width: 1.0,
                color: Color32::BLACK,
            },
        );
    }
}
