use eframe::egui::Vec2;

#[derive(Clone, Debug, Default)]
pub struct RowBlank {
    line_height: f32,
}

impl RowBlank {
    #[inline(always)]
    pub fn new() -> Self {
        Self::default()
    }

    #[inline(always)]
    pub fn blocks(&self) -> usize {
        0
    }

    pub fn line_height(&self) -> f32 {
        self.line_height
    }

    pub fn compute_line_height(&mut self, size: &Vec2) {
        self.line_height = size.y;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_row_blank() {
        let row = RowBlank::default();
        assert_eq!(row.blocks(), 0);
    }
}
