#[derive(Clone, Debug, Default)]
pub enum MethylationMode {
    #[default]
    DCM,
    DAM,
}

#[derive(Clone, Debug, Default)]
pub struct MethylationSites {
    sites: Vec<usize>,
    last_mode: MethylationMode,
}

impl MethylationSites {
    pub fn new_from_sequence(sequence: &[u8], mode: MethylationMode) -> Self {
        let mut ret = Self::default();
        for pos in 0..sequence.len() {
            match mode {
                MethylationMode::DAM => {
                    let bases = [
                        Self::get_nucleotide(sequence, pos),
                        Self::get_nucleotide(sequence, pos + 1),
                        Self::get_nucleotide(sequence, pos + 2),
                        Self::get_nucleotide(sequence, pos + 3),
                    ];
                    if bases == ['G', 'A', 'T', 'C'] {
                        ret.sites.push(pos + 1);
                    }
                }
                MethylationMode::DCM => {
                    let bases = [
                        Self::get_nucleotide(sequence, pos),
                        Self::get_nucleotide(sequence, pos + 1),
                        Self::get_nucleotide(sequence, pos + 2), // W=A|T
                        Self::get_nucleotide(sequence, pos + 3),
                        Self::get_nucleotide(sequence, pos + 4),
                    ];
                    if bases == ['C', 'C', 'A', 'G', 'G'] || bases == ['C', 'C', 'T', 'G', 'G'] {
                        ret.sites.push(pos + 1);
                    }
                }
            }
        }
        ret.last_mode = mode;
        ret
    }

    #[inline(always)]
    fn get_nucleotide(sequence: &[u8], pos: usize) -> char {
        sequence.get(pos).map_or('N', |c| *c as char)
    }

    #[inline(always)]
    pub fn sites(&self) -> &[usize] {
        &self.sites
    }

    #[inline(always)]
    pub fn last_mode(&self) -> MethylationMode {
        self.last_mode.to_owned()
    }
}
