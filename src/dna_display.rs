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
    show_gc_contents: bool,
    update_layout: UpdateLayoutParts,
    aa_letters: AminoAcidLetters,
    aa_frame: AminoAcidFrame,
}

impl DnaDisplay {
    pub fn show_re(&self) -> bool {
        self.show_restriction_enzymes
    }

    pub fn set_re(&mut self, show_re: bool) {
        self.show_restriction_enzymes = show_re;
    }

    pub fn show_rc(&self) -> bool {
        self.show_reverse_complement
    }

    pub fn set_rc(&mut self, show_rc: bool) {
        self.show_reverse_complement = show_rc;
    }

    pub fn show_orf(&self) -> bool {
        self.show_open_reading_frames
    }

    pub fn set_orf(&mut self, show_orf: bool) {
        self.show_open_reading_frames = show_orf;
    }

    pub fn show_gc(&self) -> bool {
        self.show_gc_contents
    }

    pub fn set_gc(&mut self, show_gc: bool) {
        self.show_gc_contents = show_gc;
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
}

impl Default for DnaDisplay {
    fn default() -> Self {
        Self {
            show_restriction_enzymes: true,
            show_reverse_complement: true,
            show_open_reading_frames: true,
            show_gc_contents: true,
            update_layout: UpdateLayoutParts::default(),
            aa_letters: AminoAcidLetters::Single,
            aa_frame: AminoAcidFrame::None,
        }
    }
}
