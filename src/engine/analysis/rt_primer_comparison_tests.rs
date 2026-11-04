//! Hand-crafted synthetic transcript/CDS/oligo fixtures, reconstructed by the helpers below.
//! These tests use no real genes, private reads, biological references or physical oligos.
//! All temporary request files are recreated deterministically for shared-shell parity.

use super::*;
use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};

const TAIL: &str = "ACGTCAGTACGA";
const HANDLE: &str = "TTTTTT";

fn word() -> String {
    GentleEngine::reverse_complement(TAIL)
}

fn primer(id: &str, tail: &str) -> RnaLibraryPrimer {
    RnaLibraryPrimer {
        primer_id: id.into(),
        assay_label: "Synthetic RT comparison".into(),
        sequence: format!("{HANDLE}{tail}"),
        components: vec![
            RnaPrimerComponent {
                component_id: "handle".into(),
                role: "common_handle".into(),
                start: 0,
                end: HANDLE.len(),
            },
            RnaPrimerComponent {
                component_id: "tail".into(),
                role: "target_core".into(),
                start: HANDLE.len(),
                end: HANDLE.len() + tail.len(),
            },
        ],
        core_orientation: PrimerCoreOrientation::Antisense3prime,
        expected_targets: vec![],
        variant: "synthetic-v1".into(),
        source: "Literal invented regression oligo, no physical stock".into(),
        source_sha256: None,
    }
}

fn request() -> CompareRtPrimerCandidatesRequest {
    CompareRtPrimerCandidatesRequest {
        schema: RT_PRIMER_COMPARISON_REQUEST_SCHEMA.into(),
        assessment: AssessRtPrimerCandidatesRequest {
            schema: RT_PRIMER_ASSESSMENT_REQUEST_SCHEMA.into(),
            primers: vec![primer("pilot", TAIL)],
            references: vec![RtPrimerBackground {
                reference_id: "not_supplied".into(),
                scope: RtPrimerBackgroundScope::Transcriptome,
                source: "Explicitly unavailable synthetic background".into(),
                fasta: None,
                unavailable_reason: Some(
                    "No prospective background supplied by this fixture".into(),
                ),
                circular_record_ids: vec![],
            }],
            intended_sites: vec![],
            budgets: RtPrimerAssessmentBudgets::default(),
        },
        sources: vec![RtPrimerComparisonSource {
            source: TranscriptCaptureSource {
                seq_id: "locus".into(),
                source_feature_id: 0,
                coverage_universe: TranscriptAssayCoverageUniverse::default(),
                annotation_release: Some("synthetic-v1".into()),
            },
            expected_source_record_sha256: None,
        }],
        budgets: RtPrimerComparisonBudgets::default(),
    }
}

/// Each tuple is a transcript ID, biological 5′–3′ cDNA, and explicit local CDS ranges.
/// None means that CDS annotation was not supplied; it does not assert noncoding biology.
fn fixture(
    transcripts: &[(&str, &str, Option<Vec<(usize, usize)>>)],
    reverse: bool,
) -> GentleEngine {
    let mut sequence = String::new();
    let mut features = Vec::new();
    for (id, cdna, cds) in transcripts {
        let start = sequence.len();
        sequence.push_str(&if reverse {
            GentleEngine::reverse_complement(cdna)
        } else {
            cdna.to_string()
        });
        let end = sequence.len();
        let range = gb_io::seq::Location::simple_range(start as i64, end as i64);
        let mut qualifiers = vec![
            ("gene".into(), Some("synthetic_group".into())),
            ("transcript_id".into(), Some(id.to_string())),
        ];
        if let Some(cds) = cds {
            let ranges = cds
                .iter()
                .map(|&(a, z)| {
                    let (a, z) = if reverse {
                        (end - z, end - a)
                    } else {
                        (start + a, start + z)
                    };
                    format!("{}-{}", a + 1, z)
                })
                .collect::<Vec<_>>()
                .join(",");
            qualifiers.push(("cds_ranges_1based".into(), Some(ranges)));
        }
        features.push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: if reverse {
                gb_io::seq::Location::Complement(Box::new(range))
            } else {
                range
            },
            qualifiers,
        });
        sequence.push_str("NNNNNNNNNNNNNNNN");
    }
    let mut dna = DNAsequence::from_sequence(&sequence).unwrap();
    *dna.features_mut() = features;
    let mut state = ProjectState::default();
    state.sequences.insert("locus".into(), dna);
    GentleEngine::from_state(state)
}

fn run(engine: &GentleEngine, req: CompareRtPrimerCandidatesRequest) -> RtPrimerComparisonReport {
    engine
        .compare_rt_primer_candidates(req, &mut |_| true)
        .unwrap()
}

fn row<'a>(
    report: &'a RtPrimerComparisonReport,
    transcript: &str,
) -> &'a RtPrimerIsoformAssessment {
    let member = report
        .members
        .iter()
        .find(|m| m.transcript_id == transcript)
        .unwrap();
    report
        .isoforms
        .iter()
        .find(|r| r.member_id == member.member_id && r.primer_id == "pilot")
        .unwrap()
}

#[test]
fn rt_primer_comparison_per_isoform_cds_utr_mixed_unknown_on_both_strands() {
    for reverse in [false, true] {
        let cdna = format!("AAAA{}TTTT", word());
        let engine = fixture(
            &[
                ("five", &cdna, Some(vec![(16, 20)])),
                ("coding", &cdna, Some(vec![(4, 16)])),
                ("three", &cdna, Some(vec![(0, 4)])),
                ("boundary", &cdna, Some(vec![(8, 20)])),
                ("unknown", &cdna, None),
                ("gapped", &cdna, Some(vec![(4, 8), (12, 16)])),
            ],
            reverse,
        );
        let report = run(&engine, request());
        assert_eq!(report.members.len(), 6);
        for (tx, expected) in [
            ("five", RtPrimerBindingRegion::FivePrimeUtr),
            ("coding", RtPrimerBindingRegion::Cds),
            ("three", RtPrimerBindingRegion::ThreePrimeUtr),
            ("boundary", RtPrimerBindingRegion::Mixed),
            ("unknown", RtPrimerBindingRegion::Unknown),
            ("gapped", RtPrimerBindingRegion::Unknown),
        ] {
            let r = row(&report, tx);
            assert_eq!(r.presence, RtPrimerSitePresence::Present, "{tx}");
            assert_eq!(r.bindings.len(), 1, "{tx}");
            assert_eq!(r.bindings[0].region, expected, "{tx}");
            assert_eq!(
                (
                    r.bindings[0].transcript_start_0based,
                    r.bindings[0].transcript_end_0based_exclusive
                ),
                (4, 16)
            );
            let binding = &r.bindings[0];
            assert_eq!(
                binding.cds_bases
                    + binding.five_prime_utr_bases
                    + binding.three_prime_utr_bases
                    + binding.unknown_bases,
                12
            );
        }
        let boundary = &row(&report, "boundary").bindings[0];
        assert_eq!((boundary.five_prime_utr_bases, boundary.cds_bases), (4, 8));
        let gap = &row(&report, "gapped").bindings[0];
        assert_eq!((gap.cds_bases, gap.unknown_bases), (8, 4));
        assert_eq!(row(&report, "unknown").bindings[0].unknown_bases, 12);
        assert_eq!(report.candidates[0].members_with_exact_site, 6);
        for m in &report.members {
            assert_eq!(m.strand, if reverse { "-" } else { "+" });
        }
    }
}

#[test]
fn rt_primer_comparison_junction_spanning_bindings_use_shared_source_projection() {
    for reverse in [false, true] {
        let cdna = format!("AAAA{}TTTT", word());
        let genomic = format!("{}NNNNNNNN{}", &cdna[..10], &cdna[10..]);
        let mut dna = DNAsequence::from_sequence(&if reverse {
            GentleEngine::reverse_complement(&genomic)
        } else {
            genomic
        })
        .unwrap();
        let join = gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(0, 10),
            gb_io::seq::Location::simple_range(18, 28),
        ]);
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: if reverse {
                gb_io::seq::Location::Complement(Box::new(join))
            } else {
                join
            },
            qualifiers: vec![
                ("gene".into(), Some("synthetic_group".into())),
                ("transcript_id".into(), Some("split".into())),
                ("cds_ranges_1based".into(), Some("1-10,19-28".into())),
            ],
        });
        let mut state = ProjectState::default();
        state.sequences.insert("locus".into(), dna);
        let report = run(&GentleEngine::from_state(state), request());
        let binding = &row(&report, "split").bindings[0];
        assert_eq!(binding.region, RtPrimerBindingRegion::Cds);
        assert_eq!(binding.cds_bases, 12);
        assert_eq!(
            binding.source_ranges_0based,
            if reverse {
                vec![(18, 24), (4, 10)]
            } else {
                vec![(4, 10), (18, 24)]
            }
        );
    }
}

#[test]
fn rt_primer_comparison_absent_possible_and_exact_presence_are_distinct() {
    let exact = word();
    let mut possible = exact.clone();
    possible.replace_range(3..4, "N");
    let engine = fixture(
        &[
            ("exact", &exact, None),
            ("possible", &possible, None),
            ("mismatching_ambiguity", "AAAAAAAAAAAN", None),
            ("absent", "AAAAAAAAAAAA", None),
            ("short", "ACGT", None),
        ],
        false,
    );
    let report = run(&engine, request());
    assert_eq!(
        row(&report, "exact").presence,
        RtPrimerSitePresence::Present
    );
    let r = row(&report, "possible");
    assert_eq!(r.presence, RtPrimerSitePresence::Unresolved);
    assert_eq!(r.possible_ambiguous_sites, 1);
    assert!(r.bindings.is_empty());
    for tx in ["absent", "short", "mismatching_ambiguity"] {
        let r = row(&report, tx);
        assert_eq!(r.presence, RtPrimerSitePresence::Absent, "{tx}");
        assert!(r.bindings.is_empty());
        assert_eq!(r.possible_ambiguous_sites, 0);
    }
    let c = &report.candidates[0];
    assert_eq!(
        (
            c.total_members,
            c.members_with_exact_site,
            c.members_with_possible_only_site,
            c.members_without_site
        ),
        (5, 1, 1, 3)
    );
}

#[test]
fn rt_primer_comparison_repeated_sites_preserve_regions_and_count_member_once() {
    let cdna = format!("{}AAA{}", word(), word());
    let engine = fixture(&[("repeat", &cdna, Some(vec![(0, 12)]))], false);
    let report = run(&engine, request());
    let r = row(&report, "repeat");
    assert_eq!(r.bindings.len(), 2);
    assert_eq!(r.bindings[0].region, RtPrimerBindingRegion::Cds);
    assert_eq!(r.bindings[1].region, RtPrimerBindingRegion::ThreePrimeUtr);
    let c = &report.candidates[0];
    assert_eq!((c.total_members, c.members_with_exact_site), (1, 1));
    assert_eq!(c.region_member_counts[&RtPrimerBindingRegion::Cds], 1);
    assert_eq!(
        c.region_member_counts[&RtPrimerBindingRegion::ThreePrimeUtr],
        1
    );
}

#[test]
fn rt_primer_comparison_explicit_universe_preserves_requested_denominator_and_refuses_unknowns() {
    let cdna = word();
    let engine = fixture(
        &[
            ("present", &cdna, None),
            ("absent", "AAAAAAAAAAAA", None),
            ("omitted", &cdna, None),
        ],
        false,
    );
    let mut req = request();
    req.sources[0].source.coverage_universe = TranscriptAssayCoverageUniverse {
        kind: TranscriptAssayCoverageUniverseKind::ExplicitTranscripts,
        required_transcript_ids: vec!["present".into(), "absent".into()],
        ..Default::default()
    };
    let report = run(&engine, req.clone());
    assert_eq!(report.members.len(), 2);
    assert_eq!(
        (
            report.candidates[0].total_members,
            report.candidates[0].members_with_exact_site,
            report.candidates[0].members_without_site
        ),
        (2, 1, 1)
    );
    assert!(!report.members.iter().any(|m| m.transcript_id == "omitted"));
    req.sources[0]
        .source
        .coverage_universe
        .required_transcript_ids
        .push("unresolved".into());
    assert!(
        engine
            .compare_rt_primer_candidates(req, &mut |_| true)
            .is_err()
    );
}

#[test]
fn rt_primer_comparison_source_record_guard_detects_annotation_only_changes() {
    let cdna = word();
    let mut engine = fixture(&[("tx", &cdna, Some(vec![(0, 12)]))], false);
    let baseline = run(&engine, request());
    let mut req = request();
    req.sources[0].expected_source_record_sha256 =
        Some(baseline.sources[0].source_record_sha256.clone());
    assert_eq!(
        run(&engine, req.clone()).content_sha256,
        run(&engine, req.clone()).content_sha256
    );
    engine
        .state_mut()
        .sequences
        .get_mut("locus")
        .unwrap()
        .features_mut()[0]
        .qualifiers
        .push((
            "note".into(),
            Some("annotation changed without changing sequence letters".into()),
        ));
    assert!(
        engine
            .compare_rt_primer_candidates(req, &mut |_| true)
            .is_err()
    );
    let fresh = run(&engine, request());
    assert_ne!(
        baseline.sources[0].source_record_sha256,
        fresh.sources[0].source_record_sha256
    );
    assert_eq!(
        baseline.members[0].cdna_sha256,
        fresh.members[0].cdna_sha256
    );
}

#[test]
fn rt_primer_comparison_candidate_metrics_preserve_input_order_without_ranking() {
    let cdna = word();
    let engine = fixture(&[("tx", &cdna, Some(vec![(0, 12)]))], false);
    let mut req = request();
    req.assessment
        .primers
        .push(primer("uncovered", "AAAAAAAAAAAA"));
    let report = run(&engine, req);
    assert_eq!(
        report
            .candidates
            .iter()
            .map(|c| c.primer_id.as_str())
            .collect::<Vec<_>>(),
        ["pilot", "uncovered"]
    );
    assert_eq!(report.candidates[0].target_tail.length_bp, 12);
    assert_eq!(report.candidates[0].full_oligo.length_bp, 18);
    assert!((report.candidates[0].target_tail.gc_fraction - 0.5).abs() < 1e-12);
    assert!((report.candidates[0].full_oligo.gc_fraction - 1.0 / 3.0).abs() < 1e-12);
    assert!(report.candidates[0].target_tail.tm_c.is_finite());
    assert!(!report.tm_model.is_empty());
    assert_eq!(report.candidates[1].members_with_exact_site, 0);
}

#[test]
fn rt_primer_comparison_composition_preserves_assessment_bytes_and_determinism() {
    let cdna = word();
    let engine = fixture(&[("tx", &cdna, None)], false);
    let req = request();
    let old = engine
        .assess_rt_primer_candidates(req.assessment.clone(), &mut |_| true)
        .unwrap();
    let report = run(&engine, req.clone());
    assert_eq!(
        serde_json::to_value(&*report.assessment).unwrap(),
        serde_json::to_value(old).unwrap()
    );
    assert_eq!(report.content_sha256, run(&engine, req).content_sha256);
    assert_eq!(report.assessment.references[0].status, "not_assessed");
}

#[test]
fn rt_primer_comparison_limits_cancellation_schema_and_duplicate_sources_fail_closed() {
    let cdna = format!("{}AAA{}", word(), word());
    let engine = fixture(&[("tx", &cdna, None)], false);
    for which in 0..5 {
        let mut req = request();
        match which {
            0 => req.budgets.max_transcript_base_comparisons = 1,
            1 => req.budgets.max_bindings = 1,
            2 => req.sources.push(req.sources[0].clone()),
            3 => req.sources[0].source.source_feature_id = usize::MAX,
            _ => {
                req.sources[0].expected_source_record_sha256 =
                    Some(format!("sha256:{}", "0".repeat(64)))
            }
        }
        assert!(
            engine
                .compare_rt_primer_candidates(req, &mut |_| true)
                .is_err(),
            "case {which}"
        );
    }
    assert!(
        engine
            .compare_rt_primer_candidates(request(), &mut |_| false)
            .is_err()
    );
    let mut value = serde_json::to_value(request()).unwrap();
    value["unexpected_fallback"] = serde_json::json!(true);
    assert!(serde_json::from_value::<CompareRtPrimerCandidatesRequest>(value).is_err());
}

#[test]
fn rt_primer_comparison_shared_shell_and_typed_route_preserve_content_history_and_revision() {
    let cdna = word();
    let mut engine = fixture(&[("tx", &cdna, Some(vec![(0, 12)]))], true);
    let req = request();
    let before = serde_json::to_value(engine.snapshot()).unwrap();
    let revision = engine.mutation_revision();
    let undo = engine.undo_available();
    let result = engine
        .apply(Operation::CompareRtPrimerCandidates {
            request: Box::new(req.clone()),
        })
        .unwrap();
    assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), before);
    assert_eq!(engine.mutation_revision(), revision);
    assert_eq!(engine.undo_available(), undo);
    let journal = engine.journal_len();
    let expected = serde_json::to_value(result.rt_primer_comparison.unwrap()).unwrap();
    let command = parse_shell_line(&format!(
        "primers compare-rt-candidates {}",
        quote_shell_arg(&serde_json::to_string(&req).unwrap())
    ))
    .unwrap();
    assert!(!command.is_state_mutating());
    let actual = execute_shell_command(&mut engine, &command).unwrap();
    assert!(!actual.state_changed);
    assert_eq!(actual.output["result"]["rt_primer_comparison"], expected);
    assert_eq!(engine.journal_len(), journal);
    assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), before);
    let td = tempfile::tempdir().unwrap();
    let path = td.path().join("comparison request's input.json");
    std::fs::write(&path, serde_json::to_vec(&req).unwrap()).unwrap();
    let command = parse_shell_line(&format!(
        "primers compare-rt-candidates {}",
        quote_shell_arg(&format!("@{}", path.display()))
    ))
    .unwrap();
    assert_eq!(
        execute_shell_command(&mut engine, &command).unwrap().output["result"]["rt_primer_comparison"],
        expected
    );
}

#[test]
fn rt_primer_comparison_partial_malformed_cds_refuses_without_invented_utr() {
    let cdna = word();
    let mut engine = fixture(&[("tx", &cdna, None)], false);
    engine
        .state_mut()
        .sequences
        .get_mut("locus")
        .unwrap()
        .features_mut()[0]
        .qualifiers
        .push(("cds_ranges_1based".into(), Some("1-4,invalid".into())));
    assert!(
        engine
            .compare_rt_primer_candidates(request(), &mut |_| true)
            .is_err()
    );
}

#[test]
fn rt_primer_comparison_unparseable_or_empty_cds_remains_unknown() {
    let cdna = word();
    for value in ["invalid", "", " , "] {
        let mut engine = fixture(&[("tx", &cdna, None)], false);
        engine
            .state_mut()
            .sequences
            .get_mut("locus")
            .unwrap()
            .features_mut()[0]
            .qualifiers
            .push(("cds_ranges_1based".into(), Some(value.into())));
        let report = run(&engine, request());
        let r = row(&report, "tx");
        assert_eq!(r.presence, RtPrimerSitePresence::Present);
        assert_eq!(r.bindings[0].region, RtPrimerBindingRegion::Unknown);
        assert_eq!(r.bindings[0].unknown_bases, TAIL.len());
        assert!(report.members[0].cds_ranges_0based.is_empty());
    }
}

#[test]
fn rt_primer_comparison_conflicting_cds_models_and_clipped_exons_refuse() {
    let cdna = word();
    let mut engine = fixture(&[("tx", &cdna, None)], false);
    for (start, end, protein) in [(0, 8, "protein_a"), (4, 12, "protein_b")] {
        engine
            .state_mut()
            .sequences
            .get_mut("locus")
            .unwrap()
            .features_mut()
            .push(gb_io::seq::Feature {
                kind: "CDS".into(),
                location: gb_io::seq::Location::simple_range(start, end),
                qualifiers: vec![
                    ("gene".into(), Some("synthetic_group".into())),
                    ("transcript_id".into(), Some("tx".into())),
                    ("protein_id".into(), Some(protein.into())),
                ],
            });
    }
    assert!(
        engine
            .compare_rt_primer_candidates(request(), &mut |_| true)
            .is_err()
    );
    let mut clipped = fixture(&[("tx", &cdna, None)], false);
    // Fixture locus contains the 12-base transcript followed by 16 spacer Ns.
    // This annotated exon extends beyond all 28 loaded source bases.
    clipped
        .state_mut()
        .sequences
        .get_mut("locus")
        .unwrap()
        .features_mut()[0]
        .location = gb_io::seq::Location::simple_range(0, 40);
    assert!(
        clipped
            .compare_rt_primer_candidates(request(), &mut |_| true)
            .is_err()
    );
}

#[test]
fn rt_primer_comparison_annotation_complexity_is_refused_before_shared_splicing_view() {
    let cdna = word();
    let mut repeated = fixture(&[("tx", &cdna, None)], false);
    repeated
        .state_mut()
        .sequences
        .get_mut("locus")
        .unwrap()
        .features_mut()[0]
        .location = gb_io::seq::Location::Join(
        (0..17_000)
            .map(|_| gb_io::seq::Location::simple_range(0, 1))
            .collect(),
    );
    let error = repeated
        .compare_rt_primer_candidates(request(), &mut |_| true)
        .unwrap_err();
    assert!(
        error.message.contains("annotation_complexity"),
        "{}",
        error.message
    );

    let mut deep = fixture(&[("tx", &cdna, None)], false);
    let mut location = gb_io::seq::Location::simple_range(0, 12);
    for _ in 0..40 {
        location = gb_io::seq::Location::Complement(Box::new(location));
    }
    deep.state_mut()
        .sequences
        .get_mut("locus")
        .unwrap()
        .features_mut()[0]
        .location = location;
    let error = deep
        .compare_rt_primer_candidates(request(), &mut |_| true)
        .unwrap_err();
    assert!(
        error.message.contains("annotation_complexity"),
        "{}",
        error.message
    );
}
