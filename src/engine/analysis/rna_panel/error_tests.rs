//! Synthetic A/C/G/T error regressions for the panel's shared nomination core.
//!
//! Origin/recreation: `tests::bases(1000 + locus, 2000)` makes two unrelated toy
//! loci; each has two 1200-base transcripts sharing exon [0, 600), with distinct
//! second exons [700, 1300) and [1400, 2000). No private or biological reads.
//! Literal scripts below use original, zero-based transcript coordinates:
//! substitution rotates A->C->G->T->A, insertion precedes the coordinate, and
//! deletion removes the indicated source bases. Indel events are 1--3 bases.
//! Error denominators count original bases, not edit events or alignment length:
//! 2% = (12 substitutions + 6 inserted + 6 deleted bases) / 1200 = 24/1200;
//! 8% = (48 substitutions + 24 inserted + 24 deleted bases) / 1200 = 96/1200.
//! Balanced mixed scripts leave reads 1200 bases long. Isolated indel cases use
//! 1, 2 or 3 edited bases / 1200, yielding 1197--1203 bases. The junction control
//! isolates one exact 15-base word; its flanked negative has one substitution
//! / 45 original bases and retains exact anchors outside that word.
//! These fixtures exercise `build_index`/`nominate`, not basecaller sensitivity,
//! alignment confirmation, exclusive isoform assignment or real-library rates.

use super::tests::{bases, starter};
use super::*;

const READ_BASES: usize = 1200;
const SHARED_EXON_BASES: usize = 600;

#[derive(Clone, Copy)]
enum Edit {
    Substitute(usize),
    Insert(usize, &'static [u8]),
    Delete(usize, usize),
}

const TWO_PERCENT_EDITS: &[Edit] = &[
    Edit::Substitute(35),
    Edit::Substitute(105),
    Edit::Substitute(175),
    Edit::Insert(200, b"A"),
    Edit::Substitute(275),
    Edit::Insert(350, b"CG"),
    Edit::Substitute(405),
    Edit::Insert(510, b"TGA"),
    Edit::Substitute(575),
    Edit::Substitute(635),
    Edit::Delete(700, 1),
    Edit::Substitute(735),
    Edit::Substitute(825),
    Edit::Delete(850, 2),
    Edit::Substitute(925),
    Edit::Delete(1000, 3),
    Edit::Substitute(1085),
    Edit::Substitute(1165),
];

const EIGHT_PERCENT_EDITS: &[Edit] = &[
    Edit::Substitute(11),
    Edit::Substitute(31),
    Edit::Insert(43, b"A"),
    Edit::Substitute(57),
    Edit::Delete(73, 1),
    Edit::Substitute(89),
    Edit::Insert(103, b"CG"),
    Edit::Substitute(119),
    Edit::Delete(133, 2),
    Edit::Substitute(149),
    Edit::Insert(163, b"TGA"),
    Edit::Substitute(179),
    Edit::Delete(193, 3),
    Edit::Substitute(209),
    Edit::Substitute(229),
    Edit::Substitute(249),
    Edit::Substitute(269),
    Edit::Substitute(289),
    Edit::Substitute(311),
    Edit::Substitute(331),
    Edit::Insert(343, b"A"),
    Edit::Substitute(357),
    Edit::Delete(373, 1),
    Edit::Substitute(389),
    Edit::Insert(403, b"CG"),
    Edit::Substitute(419),
    Edit::Delete(433, 2),
    Edit::Substitute(449),
    Edit::Insert(463, b"TGA"),
    Edit::Substitute(479),
    Edit::Delete(493, 3),
    Edit::Substitute(509),
    Edit::Substitute(529),
    Edit::Substitute(549),
    Edit::Substitute(569),
    Edit::Substitute(589),
    Edit::Substitute(611),
    Edit::Substitute(631),
    Edit::Insert(643, b"A"),
    Edit::Substitute(657),
    Edit::Delete(673, 1),
    Edit::Substitute(689),
    Edit::Insert(703, b"CG"),
    Edit::Substitute(719),
    Edit::Delete(733, 2),
    Edit::Substitute(749),
    Edit::Insert(763, b"TGA"),
    Edit::Substitute(779),
    Edit::Delete(793, 3),
    Edit::Substitute(809),
    Edit::Substitute(829),
    Edit::Substitute(849),
    Edit::Substitute(869),
    Edit::Substitute(889),
    Edit::Substitute(911),
    Edit::Substitute(931),
    Edit::Insert(943, b"A"),
    Edit::Substitute(957),
    Edit::Delete(973, 1),
    Edit::Substitute(989),
    Edit::Insert(1003, b"CG"),
    Edit::Substitute(1019),
    Edit::Delete(1033, 2),
    Edit::Substitute(1049),
    Edit::Insert(1063, b"TGA"),
    Edit::Substitute(1079),
    Edit::Delete(1093, 3),
    Edit::Substitute(1109),
    Edit::Substitute(1129),
    Edit::Substitute(1149),
    Edit::Substitute(1169),
    Edit::Substitute(1189),
];

fn edited_read(reference: &[u8], script: &[Edit], expected_counts: [usize; 3]) -> Vec<u8> {
    assert!(reference.iter().all(|base| b"ACGT".contains(base)));
    let mut read = Vec::new();
    let mut cursor = 0;
    let mut counts = [0; 3];
    for &edit in script {
        let at = match edit {
            Edit::Substitute(at) | Edit::Insert(at, _) | Edit::Delete(at, _) => at,
        };
        assert!(cursor <= at && at <= reference.len());
        read.extend_from_slice(&reference[cursor..at]);
        match edit {
            Edit::Substitute(_) => {
                read.push(match reference[at] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    b'T' => b'A',
                    _ => unreachable!(),
                });
                cursor = at + 1;
                counts[0] += 1;
            }
            Edit::Insert(_, inserted) => {
                assert!((1..=3).contains(&inserted.len()));
                read.extend_from_slice(inserted);
                cursor = at;
                counts[1] += inserted.len();
            }
            Edit::Delete(_, length) => {
                assert!((1..=3).contains(&length));
                assert!(at + length <= reference.len());
                cursor = at + length;
                counts[2] += length;
            }
        }
    }
    read.extend_from_slice(&reference[cursor..]);
    assert_eq!(counts, expected_counts);
    assert_eq!(read.len(), reference.len() + counts[1] - counts[2]);
    assert!(read.iter().all(|base| b"ACGT".contains(base)));
    read
}

fn long_read_panel(directory: &Path) -> RnaReadPanel {
    let (mut engine, build) = starter(2, directory);
    for locus in 0..2 {
        let mut dna = DNAsequence::from_sequence(&bases(1000 + locus as u64, 2000)).unwrap();
        for (id, second_start) in [("error_tx", 700), ("error_alt", 1400)] {
            dna.features_mut().push(gb_io::seq::Feature {
                kind: "mRNA".into(),
                location: gb_io::seq::Location::Join(vec![
                    gb_io::seq::Location::simple_range(0, SHARED_EXON_BASES as i64),
                    gb_io::seq::Location::simple_range(
                        second_start,
                        second_start + SHARED_EXON_BASES as i64,
                    ),
                ]),
                qualifiers: vec![
                    ("gene".into(), Some(format!("SYN_{locus}"))),
                    ("transcript_id".into(), Some(id.into())),
                ],
            });
        }
        engine
            .state_mut()
            .sequences
            .insert(format!("source_{locus}"), dna);
    }
    let panel = engine.build_rna_read_panel(&build).unwrap();
    assert_eq!(panel.loci.len(), 2);
    assert_eq!(panel.templates.len(), 4);
    assert!(panel.overlapping_source_pairs.is_empty());
    for pair in panel.templates.chunks_exact(2) {
        assert_eq!(pair[0].locus_id, pair[1].locus_id);
        assert_eq!(
            &pair[0].sequence[..SHARED_EXON_BASES],
            &pair[1].sequence[..SHARED_EXON_BASES]
        );
        assert_ne!(
            &pair[0].sequence[SHARED_EXON_BASES..],
            &pair[1].sequence[SHARED_EXON_BASES..]
        );
    }
    for template in &panel.templates {
        assert_eq!(template.sequence.len(), READ_BASES);
        assert_eq!(template.junction_offsets, vec![SHARED_EXON_BASES]);
    }
    panel
}

fn nominate_repeatably(
    panel: &RnaReadPanel,
    index: &Index,
    policy: &RnaReadPanelPolicy,
    read: &[u8],
) -> (String, Vec<RnaReadPanelCandidate>) {
    let first = nominate(panel, index, policy, read, &mut || false).unwrap();
    let second = nominate(panel, index, policy, read, &mut || false).unwrap();
    assert_eq!(encode(&first).unwrap(), encode(&second).unwrap());
    first
}

#[test]
fn rna_panel_two_percent_acgt_edits_retain_members_without_cross_locus_nominations() {
    let directory = tempfile::tempdir().unwrap();
    let panel = long_read_panel(directory.path());
    let policy = RnaReadPanelPolicy::default();
    let (index, _) = build_index(&panel, &policy).unwrap();
    for template in &panel.templates {
        let read = edited_read(template.sequence.as_bytes(), TWO_PERCENT_EDITS, [12, 6, 6]);
        assert_eq!(read.len(), READ_BASES);
        let (class, candidates) = nominate_repeatably(&panel, &index, &policy, &read);
        assert_eq!(class, "single_locus_candidate", "{}", template.member_id);
        assert!(candidates.iter().all(|c| c.locus_id == template.locus_id));
        // The true member survives; the shared exon also retains its alternative.
        for member in panel
            .templates
            .iter()
            .filter(|member| member.locus_id == template.locus_id)
        {
            assert!(
                candidates
                    .iter()
                    .any(|c| c.member_id == member.member_id && c.read_orientation == "as_read"),
                "read from {} lost {}: {candidates:?}",
                template.member_id,
                member.member_id
            );
        }
    }
}

#[test]
fn rna_panel_eight_percent_acgt_edits_characterize_without_sensitivity_requirement() {
    let directory = tempfile::tempdir().unwrap();
    let panel = long_read_panel(directory.path());
    let policy = RnaReadPanelPolicy::default();
    let (index, _) = build_index(&panel, &policy).unwrap();
    for template in &panel.templates {
        let read = edited_read(
            template.sequence.as_bytes(),
            EIGHT_PERCENT_EDITS,
            [48, 24, 24],
        );
        assert_eq!(read.len(), READ_BASES);
        let (class, candidates) = nominate_repeatably(&panel, &index, &policy, &read);
        eprintln!(
            "synthetic 8% RNA errors: source={}, edited_bases=96/1200, class={class}, candidates={candidates:?}",
            template.member_id
        );
        // Dropout is allowed: this records the core's outcome, not an 8% recall floor.
        assert!(matches!(
            class.as_str(),
            "single_locus_candidate" | "partial" | "unmatched"
        ));
        assert!(candidates.iter().all(|c| c.locus_id == template.locus_id));
    }
}

#[test]
fn rna_panel_one_locus_one_to_three_bp_indels_do_not_make_multi_segment_candidates() {
    let directory = tempfile::tempdir().unwrap();
    let panel = long_read_panel(directory.path());
    let policy = RnaReadPanelPolicy::default();
    let (index, _) = build_index(&panel, &policy).unwrap();
    for template in &panel.templates {
        for length in 1..=3 {
            for (edit, counts) in [
                (
                    Edit::Insert(SHARED_EXON_BASES, &b"ACG"[..length]),
                    [0, length, 0],
                ),
                (Edit::Delete(SHARED_EXON_BASES, length), [0, 0, length]),
            ] {
                let read = edited_read(template.sequence.as_bytes(), &[edit], counts);
                assert!(read.len() >= 1000);
                let (class, candidates) = nominate_repeatably(&panel, &index, &policy, &read);
                assert_eq!(class, "single_locus_candidate");
                assert_ne!(class, "multi_segment_candidate");
                assert!(
                    candidates
                        .iter()
                        .any(|c| c.member_id == template.member_id
                            && c.read_orientation == "as_read")
                );
                assert!(candidates.iter().all(|c| c.locus_id == template.locus_id));
            }
        }
    }
}

#[test]
fn rna_panel_mutated_junction_word_loses_support_but_exact_control_has_it() {
    let directory = tempfile::tempdir().unwrap();
    let panel = long_read_panel(directory.path());
    let template = &panel.templates[0];
    let junction = template.junction_offsets[0];
    let policy = RnaReadPanelPolicy {
        min_chain_bases: 15,
        ..Default::default()
    };
    let (index, _) = build_index(&panel, &policy).unwrap();
    // Isolate the middle 15-mer: k equals control length, so its sole anchor
    // must span the junction with 7/8-base flanks, regardless of LCSk++ ties.
    let exact = &template.sequence.as_bytes()[junction - 22..junction + 23];
    let exact_word = &exact[15..30];
    assert_eq!(exact_word.len(), policy.kmer_len);
    let (class, control) = nominate_repeatably(&panel, &index, &policy, exact_word);
    assert_eq!(class, "single_locus_candidate");
    let control = control
        .iter()
        .find(|c| c.member_id == template.member_id && c.read_orientation == "as_read")
        .expect("exact junction control must nominate the true member");
    assert_eq!(control.chain_bases, 15);
    assert_eq!(control.supported_junction_offsets, vec![junction]);

    let mutated = edited_read(exact, &[Edit::Substitute(22)], [1, 0, 0]);
    // No qualifying exact junction word may survive elsewhere in the edited read.
    for start in junction + policy.min_junction_flank_bp - policy.kmer_len
        ..=junction - policy.min_junction_flank_bp
    {
        let word = &template.sequence.as_bytes()[start..start + policy.kmer_len];
        assert!(
            !mutated
                .windows(policy.kmer_len)
                .any(|window| window == word)
        );
    }
    let (class, candidates) = nominate_repeatably(&panel, &index, &policy, &mutated);
    assert_eq!(class, "single_locus_candidate");
    let candidate = candidates
        .iter()
        .find(|c| c.member_id == template.member_id && c.read_orientation == "as_read")
        .expect("flanking exact anchors must still nominate the edited member");
    assert!(candidate.supported_junction_offsets.is_empty());
}
