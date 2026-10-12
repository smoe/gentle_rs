//! Restriction-enzyme site model and cut geometry utilities.

use crate::{
    dna_sequence::DNAsequence, engine::RestrictionEnzymeDisplayMode, iupac_code::IupacCode,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, fmt};

pub fn normalize_restriction_enzyme_name(name: &str) -> String {
    name.chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_ascii_uppercase())
        .collect()
}

pub fn normalize_preferred_restriction_enzyme_names(names: &[String]) -> Vec<String> {
    let mut out = Vec::new();
    let mut seen = BTreeSet::new();
    for raw in names {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            continue;
        }
        let normalized = normalize_restriction_enzyme_name(trimmed);
        if normalized.is_empty() || !seen.insert(normalized) {
            continue;
        }
        out.push(trimmed.to_string());
    }
    out
}

pub fn restriction_group_matches_display_mode(
    mode: RestrictionEnzymeDisplayMode,
    preferred_names: &[String],
    key: &RestrictionEnzymeKey,
    names: &[String],
) -> bool {
    let preferred = preferred_names
        .iter()
        .map(|name| normalize_restriction_enzyme_name(name))
        .collect::<BTreeSet<_>>();
    let is_preferred = names.iter().any(|name| {
        let normalized = normalize_restriction_enzyme_name(name);
        !normalized.is_empty() && preferred.contains(&normalized)
    });
    let is_unique = key.number_of_cuts() == 1;
    match mode {
        RestrictionEnzymeDisplayMode::PreferredOnly => is_preferred,
        RestrictionEnzymeDisplayMode::PreferredAndUnique => is_preferred || is_unique,
        RestrictionEnzymeDisplayMode::UniqueOnly => is_unique,
        RestrictionEnzymeDisplayMode::AllInView => true,
    }
}

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
        (
            self.cut_bounds(),
            self.pos,
            self.mate_pos,
            self.cut_size,
            self.number_of_cuts,
            self.from,
            self.to,
        )
            .cmp(&(
                other.cut_bounds(),
                other.pos,
                other.mate_pos,
                other.cut_size,
                other.number_of_cuts,
                other.from,
                other.to,
            ))
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
        self.is_palindromic = self.is_palindromic();
    }

    /// Derive identity from the current motif, including after deserialization or editing.
    pub fn is_palindromic(&self) -> bool {
        !self.sequence.is_empty()
            && self
                .sequence
                .bytes()
                .zip(self.sequence.bytes().rev())
                .all(|(a, b)| {
                    let a = IupacCode::from_letter(a);
                    !a.is_empty() && a == IupacCode::from_letter(b).complement()
                })
    }

    fn get_sequence_rc(&self) -> String {
        self.sequence
            .bytes()
            .rev()
            .map(|base| IupacCode::from_letter(base).complement().to_letter() as char)
            .collect()
    }

    /// Find definite recognition sites, ordered by start and then forward/reverse orientation.
    /// Every possible template base must be allowed by the motif's IUPAC code; an unknown
    /// template N does not establish a specific A/C/G/T match. Palindromes are emitted once.
    /// Circular sites may cross the origin, but the motif must fit within one molecule.
    /// `max_sites` retains the historical all-or-nothing filter, not a truncated hit list.
    pub fn get_sites(
        &self,
        seq: &DNAsequence,
        max_sites: Option<usize>,
    ) -> Vec<RestrictionEnzymeSite> {
        let mut ret = vec![];
        let recognition_len = self.sequence.len();
        if recognition_len == 0 || seq.len() < recognition_len {
            return ret;
        }
        let motif: Vec<_> = self.sequence.bytes().map(IupacCode::from_letter).collect();
        if motif.iter().any(IupacCode::is_empty) {
            return ret;
        }
        let reverse: Vec<_> = self
            .get_sequence_rc()
            .bytes()
            .map(IupacCode::from_letter)
            .collect();
        let palindrome = motif == reverse;
        let seq_len = if seq.is_circular() {
            seq.len()
        } else {
            seq.len() - recognition_len + 1
        };
        for start in 0..seq_len {
            for (forward_strand, pattern) in [(true, &motif), (false, &reverse)] {
                if !forward_strand && palindrome {
                    continue;
                }
                if pattern.iter().enumerate().all(|(index, allowed)| {
                    let base = IupacCode::from_letter(seq.get_base_or_n(start + index));
                    !base.is_empty() && base.subset(*allowed) == base
                }) {
                    ret.push(RestrictionEnzymeSite {
                        offset: start as isize,
                        enzyme: self.to_owned(),
                        forward_strand,
                    });
                    if max_sites.is_some_and(|max| ret.len() > max) {
                        return vec![];
                    }
                }
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
    /// Absolute top/bottom cut coordinates in the displayed reference orientation.
    /// Coordinates are unwrapped and may lie outside the recognition motif (Type IIS).
    pub fn strand_cut_positions_unwrapped(&self) -> Option<(isize, isize)> {
        let length = isize::try_from(self.enzyme.sequence.len()).ok()?;
        let (forward, reverse) = self.enzyme.strand_cut_offsets();
        let (forward, reverse) = if self.forward_strand {
            (forward, reverse)
        } else {
            (length.checked_sub(reverse)?, length.checked_sub(forward)?)
        };
        Some((
            self.offset.checked_add(forward)?,
            self.offset.checked_add(reverse)?,
        ))
    }

    /// Recognition alone does not imply that both cuts fit on a linear molecule.
    pub fn can_cleave(&self, seq_len: usize, circular: bool) -> bool {
        if self
            .recognition_bounds_for_topology(seq_len, circular)
            .is_none()
        {
            return false;
        }
        let Some((forward, reverse)) = self.strand_cut_positions_unwrapped() else {
            return false;
        };
        if circular {
            forward.abs_diff(reverse) < seq_len
        } else {
            forward >= 0
                && reverse >= 0
                && forward as usize <= seq_len
                && reverse as usize <= seq_len
                && (forward != reverse || (forward > 0 && (forward as usize) < seq_len))
        }
    }

    pub fn recognition_bounds_0based(&self, seq_len: usize) -> Option<(usize, usize)> {
        self.recognition_bounds_for_topology(seq_len, false)
    }

    /// Circular intervals are unrolled: an end above `seq_len` crosses the origin.
    pub fn recognition_bounds_for_topology(
        &self,
        seq_len: usize,
        circular: bool,
    ) -> Option<(usize, usize)> {
        let length = self.enzyme.sequence.len();
        let start = usize::try_from(self.offset).ok()?;
        let end = start.checked_add(length)?;
        (length > 0 && length <= seq_len && start < seq_len && (circular || end <= seq_len))
            .then_some((start, end))
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
        self.strand_cut_positions_for_topology(seq_len, false)
    }

    /// Circular opening coordinates are unrolled together, preserving overhang length/order.
    pub fn strand_cut_positions_for_topology(
        &self,
        seq_len: usize,
        circular: bool,
    ) -> Option<(usize, usize)> {
        self.recognition_bounds_for_topology(seq_len, circular)?;
        let (forward, reverse) = self.strand_cut_positions_unwrapped()?;
        if circular {
            let len = isize::try_from(seq_len).ok()?;
            let low = forward.min(reverse);
            let width = forward.abs_diff(reverse);
            if width >= seq_len {
                return None;
            }
            let start = low.rem_euclid(len) as usize;
            let end = start.checked_add(width)?;
            Some(if forward <= reverse {
                (start, end)
            } else {
                (end, start)
            })
        } else {
            let forward = usize::try_from(forward).ok()?;
            let reverse = usize::try_from(reverse).ok()?;
            (forward <= seq_len && reverse <= seq_len).then_some((forward, reverse))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence::DNAsequence;

    // Hand-crafted motifs/templates exercise strand and ambiguity rules, not vendor cut data.
    fn synthetic_enzyme(motif: &str, cut: isize, overlap: isize) -> RestrictionEnzyme {
        let mut enzyme = RestrictionEnzyme {
            name: "synthetic".into(),
            sequence: motif.into(),
            note: None,
            cut,
            overlap,
            is_palindromic: false,
        };
        enzyme.check_palimdromic();
        enzyme
    }

    #[test]
    fn restriction_nonpalindromic_reverse_only_regression() {
        let enzyme = synthetic_enzyme("GGTCTC", 1, 4);
        assert!(!enzyme.is_palindromic());
        assert_eq!(enzyme.get_sequence_rc(), "GAGACC");
        let sites = enzyme.get_sites(&DNAsequence::from_sequence("GAGACC").unwrap(), None);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].offset, 0);
        assert!(!sites[0].forward_strand);
    }

    #[test]
    fn restriction_both_strands_are_ordered_and_count_towards_the_limit() {
        let enzyme = synthetic_enzyme("ggtctc", 1, 4);
        let dna = DNAsequence::from_sequence("ttGGTCTCttGAGACCtt").unwrap();
        let sites = enzyme.get_sites(&dna, Some(2));
        assert_eq!(
            sites
                .iter()
                .map(|site| (site.offset, site.forward_strand))
                .collect::<Vec<_>>(),
            vec![(2, true), (10, false)]
        );
        assert!(enzyme.get_sites(&dna, Some(1)).is_empty());
        assert!(enzyme.get_sites(&dna, Some(0)).is_empty());
        let overlapping = synthetic_enzyme("AAGC", 1, 1)
            .get_sites(&DNAsequence::from_sequence("AAGCTT").unwrap(), None);
        assert_eq!(
            overlapping
                .iter()
                .map(|site| (site.offset, site.forward_strand))
                .collect::<Vec<_>>(),
            vec![(0, true), (2, false)]
        );
        let palindrome = synthetic_enzyme("gaatTc", 1, 4);
        assert!(palindrome.is_palindromic());
        assert_eq!(
            palindrome
                .get_sites(&DNAsequence::from_sequence("gaattc").unwrap(), Some(1))
                .len(),
            1
        );
    }

    #[test]
    fn restriction_palindrome_identity_does_not_trust_a_serialized_or_stale_flag() {
        let mut enzyme: RestrictionEnzyme = serde_json::from_value(serde_json::json!({
            "name":"synthetic", "sequence":"GGTCTC", "cut":1, "overlap":4,
            "is_palindromic":true
        }))
        .unwrap();
        assert!(!enzyme.is_palindromic());
        enzyme.sequence = "GAATTC".into();
        assert!(enzyme.is_palindromic());
        assert_eq!(
            enzyme
                .get_sites(&DNAsequence::from_sequence("GAATTC").unwrap(), None)
                .len(),
            1
        );
    }

    #[test]
    fn restriction_iupac_recognition_is_definite_not_merely_compatible() {
        for (code, allowed) in [
            ("A", "A"),
            ("C", "C"),
            ("G", "G"),
            ("T", "T"),
            ("R", "AG"),
            ("Y", "CT"),
            ("S", "CG"),
            ("W", "AT"),
            ("K", "GT"),
            ("M", "AC"),
            ("B", "CGT"),
            ("D", "AGT"),
            ("H", "ACT"),
            ("V", "ACG"),
            ("N", "ACGT"),
        ] {
            let enzyme = synthetic_enzyme(code, 0, 0);
            for base in ["A", "C", "G", "T"] {
                let sites = enzyme.get_sites(&DNAsequence::from_sequence(base).unwrap(), None);
                assert_eq!(
                    sites.iter().any(|s| s.forward_strand),
                    allowed.contains(base),
                    "motif {code}, base {base}"
                );
            }
        }
        let enzyme = synthetic_enzyme("GCNNGC", 1, 4);
        assert_eq!(
            enzyme
                .get_sites(&DNAsequence::from_sequence("GCNNGC").unwrap(), None)
                .len(),
            1
        );
        assert!(
            enzyme
                .get_sites(&DNAsequence::from_sequence("NCNNGC").unwrap(), None)
                .is_empty()
        );
        let mixed = synthetic_enzyme("ACGTWSMKRYBDHVN", 1, 1);
        assert_eq!(mixed.get_sequence_rc(), "NBDHVRYMKSWACGT");
        let sites = mixed.get_sites(
            &DNAsequence::from_sequence("NBDHVRYMKSWACGT").unwrap(),
            None,
        );
        assert_eq!(sites.len(), 1);
        assert!(!sites[0].forward_strand);
        let padded = synthetic_enzyme("GAAGACNNNNNN", 6, 4);
        assert!(
            !padded
                .get_sites(&DNAsequence::from_sequence("GAAGACACGTAC").unwrap(), None)
                .is_empty()
        );
        assert!(
            !padded
                .get_sites(&DNAsequence::from_sequence("GTACGTGTCTTC").unwrap(), None)
                .is_empty()
        );
    }

    #[test]
    fn restriction_circular_sites_cross_origin_once_on_either_strand() {
        let enzyme = synthetic_enzyme("GGTCTC", 1, 4);
        for (sequence, forward) in [("CTCAAGGT", true), ("ACCAAGAG", false)] {
            let mut dna = DNAsequence::from_sequence(sequence).unwrap();
            assert!(enzyme.get_sites(&dna, None).is_empty());
            dna.set_circular(true);
            let sites = enzyme.get_sites(&dna, None);
            assert_eq!(sites.len(), 1);
            assert_eq!((sites[0].offset, sites[0].forward_strand), (5, forward));
            assert_eq!(
                sites[0].recognition_bounds_for_topology(8, true),
                Some((5, 11))
            );
            assert_eq!(
                sites[0].strand_cut_positions_for_topology(8, true),
                Some((6, 10))
            );
        }
    }

    #[test]
    fn restriction_empty_short_or_invalid_input_cannot_match_or_panic() {
        for sequence in [b"".as_slice(), b"GGT", b"GG\xffCTC"] {
            for circular in [false, true] {
                let mut record = DNAsequence::from_sequence("").unwrap().clone_seq_record();
                record.seq = sequence.to_vec();
                record.len = Some(sequence.len());
                let mut dna = DNAsequence::from_genbank_seq(record);
                dna.set_circular(circular);
                assert!(
                    synthetic_enzyme("GGTCTC", 1, 4)
                        .get_sites(&dna, None)
                        .is_empty()
                );
            }
        }
        let dna = DNAsequence::from_sequence("GAATTC").unwrap();
        for motif in ["", "?", "G\u{00e4}T"] {
            let enzyme = synthetic_enzyme(motif, 0, 0);
            assert!(!enzyme.is_palindromic());
            assert!(enzyme.get_sites(&dna, None).is_empty());
        }
    }

    #[test]
    fn restriction_reverse_cut_geometry_swaps_strands_and_allows_external_offsets() {
        let enzyme = synthetic_enzyme("AAGC", 6, 2);
        let forward = RestrictionEnzymeSite {
            offset: 2,
            enzyme: enzyme.clone(),
            forward_strand: true,
        };
        assert_eq!(forward.strand_cut_positions_0based(12), Some((8, 10)));
        let reverse = RestrictionEnzymeSite {
            offset: 10,
            enzyme,
            forward_strand: false,
        };
        assert_eq!(reverse.strand_cut_positions_0based(18), Some((6, 8)));
        let mut edge = reverse;
        edge.offset = 0;
        assert_eq!(edge.strand_cut_positions_0based(18), None);
        assert!(!edge.can_cleave(18, false));
        assert_eq!(
            edge.strand_cut_positions_for_topology(18, true),
            Some((14, 16))
        );
        assert!(edge.can_cleave(18, true));
    }

    #[test]
    fn preferred_restriction_enzyme_names_are_normalized_and_deduplicated() {
        let names = vec![
            " EcoRI ".to_string(),
            "eco-ri".to_string(),
            "BamHI".to_string(),
            "".to_string(),
        ];

        assert_eq!(
            normalize_preferred_restriction_enzyme_names(&names),
            vec!["EcoRI".to_string(), "BamHI".to_string()]
        );
    }

    #[test]
    fn restriction_groups_follow_shared_display_mode_policy() {
        let preferred = vec!["EcoRI".to_string()];
        let unique = RestrictionEnzymeKey::new(41, 41, 0, 1, 38, 44);
        let repeated = RestrictionEnzymeKey::new(41, 41, 0, 2, 38, 44);
        let eco_ri = vec!["eco-ri".to_string()];
        let bam_hi = vec!["BamHI".to_string()];

        assert!(restriction_group_matches_display_mode(
            RestrictionEnzymeDisplayMode::PreferredOnly,
            &preferred,
            &repeated,
            &eco_ri,
        ));
        assert!(restriction_group_matches_display_mode(
            RestrictionEnzymeDisplayMode::PreferredAndUnique,
            &preferred,
            &unique,
            &bam_hi,
        ));
        assert!(!restriction_group_matches_display_mode(
            RestrictionEnzymeDisplayMode::UniqueOnly,
            &preferred,
            &repeated,
            &eco_ri,
        ));
        assert!(restriction_group_matches_display_mode(
            RestrictionEnzymeDisplayMode::AllInView,
            &preferred,
            &repeated,
            &bam_hi,
        ));
    }

    #[test]
    fn restriction_enzyme_key_order_distinguishes_shared_cut_coordinates() {
        let short_site = RestrictionEnzymeKey::new(41, 41, 0, 1, 38, 44);
        let long_site = RestrictionEnzymeKey::new(41, 41, 0, 1, 35, 47);

        assert_ne!(short_site, long_site);
        assert_ne!(short_site.cmp(&long_site), std::cmp::Ordering::Equal);
    }

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
