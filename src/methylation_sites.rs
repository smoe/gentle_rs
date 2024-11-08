#[derive(Clone, Debug, PartialEq)]
pub struct MethylationMode {
    dcm: bool,
    dam: bool,
}

impl MethylationMode {
    pub fn dcm(&self) -> bool {
        self.dcm
    }

    pub fn set_dcm(&mut self, dcm: bool) {
        self.dcm = dcm;
    }

    pub fn dam(&self) -> bool {
        self.dam
    }

    pub fn set_dam(&mut self, dam: bool) {
        self.dam = dam;
    }
}

impl Default for MethylationMode {
    fn default() -> Self {
        Self {
            dcm: true,
            dam: true,
        }
    }
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
            if mode.dam() {
                let bases = [
                    Self::get_nucleotide(sequence, pos),
                    Self::get_nucleotide(sequence, pos + 1),
                    Self::get_nucleotide(sequence, pos + 2),
                    Self::get_nucleotide(sequence, pos + 3),
                ];
                if bases == ['G', 'A', 'T', 'C'] {
                    ret.sites.push(pos + 1);
                    continue; // No need to check DCM, site is already added
                }
            }

            if mode.dcm() {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_methylation_sites() {
        let sequence = b"CCAGGATCCAGG";
        let mode = MethylationMode {
            dcm: true,
            dam: true,
        };
        let sites = MethylationSites::new_from_sequence(sequence, mode.to_owned());
        assert_eq!(sites.sites(), &[1, 5, 8]);
        assert_eq!(sites.last_mode(), mode);
    }
}
