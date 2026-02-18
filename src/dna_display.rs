use eframe::egui::Color32;

#[derive(Debug, Clone, PartialEq)]
pub struct Selection {
    from: usize,
    to: usize,
    sequence_length: usize,
}

impl Selection {
    pub fn new(from: usize, to: usize, sequence_length: usize) -> Self {
        if to > sequence_length {
            Self {
                from: to % sequence_length,
                to: from,
                sequence_length,
            }
        } else {
            Self {
                from,
                to,
                sequence_length,
            }
        }
    }

    pub fn from(&self) -> usize {
        self.from
    }

    pub fn to(&self) -> usize {
        self.to
    }

    pub fn parts(&self) -> Vec<(usize, usize)> {
        if self.from < self.to {
            vec![(self.from, self.to)]
        } else {
            vec![(self.to, self.sequence_length), (0, self.from)]
        }
    }

    pub fn contains(&self, position: usize) -> bool {
        if self.from < self.to {
            position >= self.from && position < self.to
        } else {
            position >= self.from || position < self.to
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum AminoAcidLetters {
    Single,
    Three,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AminoAcidFrame {
    None,
    Automatic,
    Forward(u8),           // 1,2,3
    ReverseCompelment(u8), // 1,2,3
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TfbsDisplayCriteria {
    pub use_llr_bits: bool,
    pub min_llr_bits: f64,
    pub use_llr_quantile: bool,
    pub min_llr_quantile: f64,
    pub use_true_log_odds_bits: bool,
    pub min_true_log_odds_bits: f64,
    pub use_true_log_odds_quantile: bool,
    pub min_true_log_odds_quantile: f64,
}

impl Default for TfbsDisplayCriteria {
    fn default() -> Self {
        Self {
            use_llr_bits: true,
            min_llr_bits: 0.0,
            use_llr_quantile: true,
            min_llr_quantile: 0.95,
            use_true_log_odds_bits: false,
            min_true_log_odds_bits: 0.0,
            use_true_log_odds_quantile: false,
            min_true_log_odds_quantile: 0.95,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct UpdateLayoutParts {
    update_map_dna: bool,
    update_map_sequence: bool,
}

impl UpdateLayoutParts {
    pub fn update_all(&mut self) {
        self.update_map_dna = true;
        self.update_map_sequence = true;
    }

    pub fn update_map_dna(&self) -> bool {
        self.update_map_dna
    }

    pub fn update_map_sequence(&self) -> bool {
        self.update_map_sequence
    }

    pub fn map_dna_updated(&mut self) {
        self.update_map_dna = false;
    }

    pub fn map_sequence_updated(&mut self) {
        self.update_map_sequence = false;
    }
}

#[derive(Debug)]
pub struct DnaDisplay {
    show_restriction_enzymes: bool,
    show_reverse_complement: bool,
    show_open_reading_frames: bool,
    show_features: bool,
    show_tfbs: bool,
    tfbs_display_criteria: TfbsDisplayCriteria,
    show_gc_contents: bool,
    show_methylation_sites: bool,
    update_layout: UpdateLayoutParts,
    aa_letters: AminoAcidLetters,
    aa_frame: AminoAcidFrame,
    selection: Option<Selection>,
}

impl DnaDisplay {
    fn mark_layout_dirty(&mut self) {
        self.update_layout.update_all();
    }

    pub fn show_restriction_enzyme_sites(&self) -> bool {
        self.show_restriction_enzymes
    }

    pub fn toggle_show_restriction_enzyme_sites(&mut self) {
        self.show_restriction_enzymes = !self.show_restriction_enzymes;
        self.mark_layout_dirty();
    }

    pub fn set_show_restriction_enzyme_sites(&mut self, value: bool) {
        if self.show_restriction_enzymes != value {
            self.show_restriction_enzymes = value;
            self.mark_layout_dirty();
        }
    }

    pub fn toggle_show_features(&mut self) {
        self.show_features = !self.show_features;
        self.mark_layout_dirty();
    }

    pub fn set_show_features(&mut self, value: bool) {
        if self.show_features != value {
            self.show_features = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_features(&self) -> bool {
        self.show_features
    }

    pub fn show_tfbs(&self) -> bool {
        self.show_tfbs
    }

    pub fn toggle_show_tfbs(&mut self) {
        self.show_tfbs = !self.show_tfbs;
        self.mark_layout_dirty();
    }

    pub fn set_show_tfbs(&mut self, value: bool) {
        if self.show_tfbs != value {
            self.show_tfbs = value;
            self.mark_layout_dirty();
        }
    }

    pub fn tfbs_display_criteria(&self) -> TfbsDisplayCriteria {
        self.tfbs_display_criteria
    }

    pub fn set_tfbs_display_criteria(&mut self, criteria: TfbsDisplayCriteria) {
        if self.tfbs_display_criteria != criteria {
            self.tfbs_display_criteria = criteria;
            self.mark_layout_dirty();
        }
    }

    pub fn show_reverse_complement(&self) -> bool {
        self.show_reverse_complement
    }

    pub fn toggle_reverse_complement(&mut self) {
        self.show_reverse_complement = !self.show_reverse_complement;
        self.mark_layout_dirty();
    }

    pub fn set_show_reverse_complement(&mut self, value: bool) {
        if self.show_reverse_complement != value {
            self.show_reverse_complement = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_open_reading_frames(&self) -> bool {
        self.show_open_reading_frames
    }

    pub fn toggle_show_open_reading_frames(&mut self) {
        self.show_open_reading_frames = !self.show_open_reading_frames;
        self.mark_layout_dirty();
    }

    pub fn set_show_open_reading_frames(&mut self, value: bool) {
        if self.show_open_reading_frames != value {
            self.show_open_reading_frames = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_gc_contents(&self) -> bool {
        self.show_gc_contents
    }

    pub fn toggle_show_gc_contents(&mut self) {
        self.show_gc_contents = !self.show_gc_contents;
        self.mark_layout_dirty();
    }

    pub fn set_show_gc_contents(&mut self, value: bool) {
        if self.show_gc_contents != value {
            self.show_gc_contents = value;
            self.mark_layout_dirty();
        }
    }

    pub fn show_methylation_sites(&self) -> bool {
        self.show_methylation_sites
    }

    pub fn toggle_show_methylation_sites(&mut self) {
        self.show_methylation_sites = !self.show_methylation_sites;
        self.mark_layout_dirty();
    }

    pub fn set_show_methylation_sites(&mut self, value: bool) {
        if self.show_methylation_sites != value {
            self.show_methylation_sites = value;
            self.mark_layout_dirty();
        }
    }

    pub fn update_layout(&self) -> &UpdateLayoutParts {
        &self.update_layout
    }

    pub fn update_layout_mut(&mut self) -> &mut UpdateLayoutParts {
        &mut self.update_layout
    }

    pub fn aa_letters(&self) -> AminoAcidLetters {
        self.aa_letters.clone()
    }

    pub fn set_aa_letters(&mut self, aa_letters: AminoAcidLetters) {
        self.aa_letters = aa_letters;
    }

    pub fn aa_frame(&self) -> AminoAcidFrame {
        self.aa_frame.clone()
    }

    pub fn set_aa_frame(&mut self, aa_frame: AminoAcidFrame) {
        self.aa_frame = aa_frame;
    }

    pub fn selection(&self) -> Option<Selection> {
        self.selection.clone()
    }

    pub fn select(&mut self, selection: Selection) {
        self.selection = Some(selection);
    }

    pub fn deselect(&mut self) {
        self.selection = None;
    }

    pub fn restriction_enzyme_group_color(cuts: usize) -> Color32 {
        match cuts {
            1 => Color32::DARK_RED,
            2 => Color32::DARK_BLUE,
            3 => Color32::DARK_GREEN,
            _ => Color32::BLACK,
        }
    }
}

impl Default for DnaDisplay {
    fn default() -> Self {
        Self {
            show_restriction_enzymes: true,
            show_reverse_complement: true,
            show_open_reading_frames: true,
            show_features: true,
            show_tfbs: false,
            tfbs_display_criteria: TfbsDisplayCriteria::default(),
            show_gc_contents: true,
            show_methylation_sites: false,
            update_layout: UpdateLayoutParts::default(),
            aa_letters: AminoAcidLetters::Single,
            aa_frame: AminoAcidFrame::None,
            selection: None,
        }
    }
}
