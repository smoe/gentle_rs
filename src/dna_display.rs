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

#[derive(Debug)]
pub struct DnaDisplay {
    show_re: bool, // Restriction enzymes
    show_rc: bool, // Reverse complement
    aa_letters: AminoAcidLetters,
    aa_frame: AminoAcidFrame,
}

impl DnaDisplay {
    pub fn show_re(&self) -> bool {
        self.show_re
    }

    pub fn set_re(&mut self, show_re: bool) {
        self.show_re = show_re;
    }

    pub fn show_rc(&self) -> bool {
        self.show_rc
    }

    pub fn set_rc(&mut self, show_rc: bool) {
        self.show_rc = show_rc;
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
            show_re: true,
            show_rc: true,
            aa_letters: AminoAcidLetters::Single,
            aa_frame: AminoAcidFrame::None,
        }
    }
}
