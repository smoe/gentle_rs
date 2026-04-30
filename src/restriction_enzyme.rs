//! Restriction-enzyme site model and cut geometry utilities.

use crate::dna_sequence::DNAsequence;
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RestrictionEndGeometry {
    Blunt,
    FivePrimeOverhang(usize),
    ThreePrimeOverhang(usize),
}

impl RestrictionEndGeometry {
    pub fn kind_label(self) -> &'static str {
        match self {
            Self::Blunt => "blunt",
            Self::FivePrimeOverhang(_) => "5prime_overhang",
            Self::ThreePrimeOverhang(_) => "3prime_overhang",
        }
    }

    pub fn display_label(self) -> String {
        match self {
            Self::Blunt => "blunt".to_string(),
            Self::FivePrimeOverhang(bp) => format!("5' overhang ({bp} bp)"),
            Self::ThreePrimeOverhang(bp) => format!("3' overhang ({bp} bp)"),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RestrictionEnzymeKey {
    pos: isize,
    #[serde(default)]
    mate_pos: Option<isize>,
    cut_size: isize,
    number_of_cuts: usize,
    from: isize,
    to: isize,
}

impl RestrictionEnzymeKey {
    pub fn new(
        pos: isize,
        mate_pos: isize,
        cut_size: isize,
        number_of_cuts: usize,
        from: isize,
        to: isize,
    ) -> Self {
        Self {
            pos,
            mate_pos: Some(mate_pos),
            cut_size,
            number_of_cuts,
            from,
            to,
        }
    }

    pub fn number_of_cuts(&self) -> usize {
        self.number_of_cuts
    }

    pub fn cut_size(&self) -> isize {
        self.cut_size
    }

    pub fn pos(&self) -> isize {
        self.pos
    }

    pub fn mate_pos(&self) -> isize {
        self.mate_pos.unwrap_or(self.pos)
    }

    pub fn from(&self) -> isize {
        self.from
    }

    pub fn to(&self) -> isize {
        self.to
    }

    pub fn cut_bounds(&self) -> (isize, isize) {
        let mate = self.mate_pos();
        if self.pos <= mate {
            (self.pos, mate)
        } else {
            (mate, self.pos)
        }
    }

    pub fn cut_geometry(&self) -> RestrictionEndGeometry {
        let mate = self.mate_pos();
        if mate == self.pos {
            RestrictionEndGeometry::Blunt
        } else if mate > self.pos {
            RestrictionEndGeometry::FivePrimeOverhang((mate - self.pos) as usize)
        } else {
            RestrictionEndGeometry::ThreePrimeOverhang((self.pos - mate) as usize)
        }
    }
}

impl PartialOrd for RestrictionEnzymeKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RestrictionEnzymeKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.cut_bounds(), self.pos).cmp(&(other.cut_bounds(), other.pos))
    }
}

impl fmt::Display for RestrictionEnzymeKey {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let j = serde_json::to_string(&self).unwrap();
        write!(f, "{}", j)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RestrictionEnzyme {
    pub name: String,
    pub sequence: String,
    pub note: Option<String>,
    pub cut: isize,
    pub overlap: isize,
    #[serde(skip_serializing, default)]
    is_palindromic: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RestrictionEnzymeSite {
    pub offset: isize,
    pub enzyme: RestrictionEnzyme,
    pub forward_strand: bool,
}

impl RestrictionEnzyme {
    pub fn check_palimdromic(&mut self) {
        self.is_palindromic = self.sequence == self.get_sequence_rc();
    }

    #[inline(always)]
    pub fn is_palindromic(&self) -> bool {
        self.is_palindromic
    }

    fn get_sequence_rc(&self) -> String {
        // TODO cache this?
        let rc = self.sequence.as_bytes();
        let rc = match std::str::from_utf8(rc) {
            Ok(rc) => rc,
            Err(_) => panic!("RestrictionEnzyme::check_palimdromic: non-utf8 char"),
        };
        rc.to_string()
    }

    pub fn get_sites(
        &self,
        seq: &DNAsequence,
        max_sites: Option<usize>,
    ) -> Vec<RestrictionEnzymeSite> {
        // TODO reverse-complement if required
        let mut ret = vec![];
        let recognition_len = self.sequence.len();
        let seq_len = if seq.is_circular() {
            seq.len()
        } else {
            if seq.len() < recognition_len {
                return ret;
            }
            seq.len() - recognition_len + 1
        };
        for start in 0..seq_len {
            let range = std::ops::Range {
                start,
                end: start + recognition_len,
            };
            let s = seq.get_range_safe(range); // Safe for circular
            if let Some(s) = s {
                let s = std::str::from_utf8(&s).unwrap().to_uppercase(); // TODO do this once?
                if s == self.sequence {
                    // TODO IUPAC
                    ret.push(RestrictionEnzymeSite {
                        offset: start as isize,
                        enzyme: self.to_owned(),
                        forward_strand: true,
                    });
                }
            }
        }
        if let Some(max) = max_sites {
            if max < ret.len() {
                return vec![];
            }
        }
        ret
    }

    /// Return the two recessed-end offsets inside the recognition sequence.
    ///
    /// GENtle uses this shared helper when a cloning workflow needs the
    /// double-stranded opening produced by a restriction digest rather than the
    /// whole recognition span. For sticky-end cutters this follows the stored
    /// `cut + overlap` geometry. For blunt cutters, many built-in catalogs only
    /// preserve bluntness and not the exact midpoint, so we fall back to the
    /// recognition midpoint as the truthful opening coordinate.
    pub fn strand_cut_offsets(&self) -> (isize, isize) {
        let forward_cut = if self.overlap == 0 {
            (self.sequence.len() / 2) as isize
        } else {
            self.cut
        };
        let reverse_cut = if self.overlap == 0 {
            forward_cut
        } else {
            forward_cut + self.overlap
        };
        (forward_cut, reverse_cut)
    }

    pub fn end_geometry(&self) -> RestrictionEndGeometry {
        if self.overlap == 0 {
            RestrictionEndGeometry::Blunt
        } else if self.overlap > 0 {
            RestrictionEndGeometry::FivePrimeOverhang(self.overlap as usize)
        } else {
            RestrictionEndGeometry::ThreePrimeOverhang(self.overlap.unsigned_abs())
        }
    }

    pub fn recessed_end_offsets(&self) -> (isize, isize) {
        let (forward_cut, reverse_cut) = self.strand_cut_offsets();
        if forward_cut <= reverse_cut {
            (forward_cut, reverse_cut)
        } else {
            (reverse_cut, forward_cut)
        }
    }
}

impl RestrictionEnzymeSite {
    pub fn recognition_bounds_0based(&self, seq_len: usize) -> Option<(usize, usize)> {
        let start = usize::try_from(self.offset).ok()?;
        let end = start.checked_add(self.enzyme.sequence.len())?;
        (end <= seq_len).then_some((start, end))
    }

    /// Return the zero-based opening window between the two recessed ends of
    /// the digested DNA arms.
    ///
    /// For sticky-end cutters this is a non-empty interval spanning the
    /// single-stranded overhang region between the recessed 3' termini. For a
    /// blunt cutter this is a zero-length cutpoint (`start == end`).
    pub fn recessed_opening_window_0based(&self, seq_len: usize) -> Option<(usize, usize)> {
        let (forward_cut, reverse_cut) = self.strand_cut_positions_0based(seq_len)?;
        Some((forward_cut.min(reverse_cut), forward_cut.max(reverse_cut)))
    }

    pub fn strand_cut_positions_0based(&self, seq_len: usize) -> Option<(usize, usize)> {
        let (recognition_start, recognition_end) = self.recognition_bounds_0based(seq_len)?;
        let (forward_offset, reverse_offset) = self.enzyme.strand_cut_offsets();
        let forward_cut = usize::try_from(self.offset.checked_add(forward_offset)?).ok()?;
        let reverse_cut = usize::try_from(self.offset.checked_add(reverse_offset)?).ok()?;
        (forward_cut >= recognition_start
            && forward_cut <= recognition_end
            && reverse_cut >= recognition_start
            && reverse_cut <= recognition_end)
            .then_some((forward_cut, reverse_cut))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence::DNAsequence;

    #[test]
    fn test_restriction_enzyme() {
        let mut re = RestrictionEnzyme {
            name: "EcoRI".to_string(),
            sequence: "GAATTC".to_string(),
            note: None,
            cut: 1,
            overlap: 1,
            is_palindromic: false,
        };
        re.check_palimdromic();
        assert!(re.is_palindromic());
        let seq = DNAsequence::from_sequence("GAATTC").unwrap();
        let sites = re.get_sites(&seq, None);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].offset, 0);
        assert!(sites[0].forward_strand);
    }

    #[test]
    fn test_restriction_enzyme_sites() {
        let mut re = RestrictionEnzyme {
            name: "EcoRI".to_string(),
            sequence: "GAATTC".to_string(),
            note: None,
            cut: 1,
            overlap: 1,
            is_palindromic: false,
        };
        re.check_palimdromic();
        assert!(re.is_palindromic());
        let seq = DNAsequence::from_sequence("GAATTCGAATTC").unwrap();
        let sites = re.get_sites(&seq, None);
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].offset, 0);
        assert!(sites[0].forward_strand);
        assert_eq!(sites[1].offset, 6);
        assert!(sites[1].forward_strand);
    }

    #[test]
    fn recessed_end_offsets_use_midpoint_for_blunt_cutters() {
        let re = RestrictionEnzyme {
            name: "SmaI".to_string(),
            sequence: "CCCGGG".to_string(),
            note: None,
            cut: 1,
            overlap: 0,
            is_palindromic: true,
        };
        assert_eq!(re.recessed_end_offsets(), (3, 3));
    }

    #[test]
    fn strand_cut_offsets_preserve_sticky_orientation() {
        let re = RestrictionEnzyme {
            name: "EcoRI".to_string(),
            sequence: "GAATTC".to_string(),
            note: None,
            cut: 1,
            overlap: 4,
            is_palindromic: true,
        };
        assert_eq!(re.strand_cut_offsets(), (1, 5));
        assert_eq!(
            re.end_geometry(),
            RestrictionEndGeometry::FivePrimeOverhang(4)
        );

        let re = RestrictionEnzyme {
            name: "KpnI".to_string(),
            sequence: "GGTACC".to_string(),
            note: None,
            cut: 5,
            overlap: -4,
            is_palindromic: true,
        };
        assert_eq!(re.strand_cut_offsets(), (5, 1));
        assert_eq!(
            re.end_geometry(),
            RestrictionEndGeometry::ThreePrimeOverhang(4)
        );
    }

    #[test]
    fn recessed_opening_window_tracks_sticky_end_recessed_termini() {
        let site = RestrictionEnzymeSite {
            offset: 10,
            enzyme: RestrictionEnzyme {
                name: "EcoRI".to_string(),
                sequence: "GAATTC".to_string(),
                note: None,
                cut: 1,
                overlap: 4,
                is_palindromic: true,
            },
            forward_strand: true,
        };
        assert_eq!(site.recognition_bounds_0based(100), Some((10, 16)));
        assert_eq!(site.strand_cut_positions_0based(100), Some((11, 15)));
        assert_eq!(site.recessed_opening_window_0based(100), Some((11, 15)));
    }

    #[test]
    fn recessed_opening_window_allows_zero_length_blunt_cutpoint() {
        let site = RestrictionEnzymeSite {
            offset: 10,
            enzyme: RestrictionEnzyme {
                name: "SmaI".to_string(),
                sequence: "CCCGGG".to_string(),
                note: None,
                cut: 1,
                overlap: 0,
                is_palindromic: true,
            },
            forward_strand: true,
        };
        assert_eq!(site.recognition_bounds_0based(100), Some((10, 16)));
        assert_eq!(site.strand_cut_positions_0based(100), Some((13, 13)));
        assert_eq!(site.recessed_opening_window_0based(100), Some((13, 13)));
    }
}
