//! Synthetic source-review counterexamples. No biological data or upstream impossible DAG paths.

use super::tests::{fixture, model, originals_fixture, path, preferred, run, run_originals};
use super::*;

fn input<'a>(
    sample: &'a str,
    analysis: &'a ObservedSamAnalysis,
    annotation: Option<&'a PrimerPayloadAnnotation>,
    cds: Option<&'a ModelDocument>,
    originals: Option<oriented::Input<'a>>,
) -> Input<'a> {
    Input {
        sample_id: sample,
        bindings: Bindings {
            report_content_sha256: format!("sha256:{}", "1".repeat(64)),
            data_sha256: pileup::hash(analysis).unwrap(),
            annotation_sha256: annotation.map(|n| pileup::hash(n).unwrap()),
            model_sha256: cds.map(|m| pileup::hash(m).unwrap()),
        },
        analysis,
        annotation,
        model: cds,
        original_bases: originals,
        diagnostics: vec![],
    }
}
fn build_input(region: &Region, input: Input<'_>, limits: Limits) -> Result<Report, EngineError> {
    build(
        "synthetic",
        region,
        &[input],
        oriented::Policy::default(),
        limits,
        &mut || Ok(()),
    )
}
fn candidate<'a>(report: &'a Report, id: &str) -> &'a Assignment {
    report.reads[0]
        .assignments
        .iter()
        .find(|a| a.model.transcript_id == id)
        .unwrap()
}
fn render_bound(report: &Report) -> (Page, gentle_render::read_end_pileup::Scene) {
    let p = page(report, &View::default()).unwrap();
    let s = scene(report, &p).unwrap();
    assert!(
        s.rows
            .iter()
            .flat_map(|r| &r.runs)
            .all(|r| (s.start..s.end).contains(&r.reference_anchor))
    );
    gentle_render::read_end_pileup::render(&s).unwrap();
    (p, s)
}

#[test]
fn read_end_locus_off_window_softclip_and_aligned_internal_runs_stay_in_json_not_glyphs() {
    let (a, n, mut r, original) = originals_fixture();
    r.flank_left = 0;
    let hash = pileup::hash(&original).unwrap();
    let report = run_originals(
        &a,
        Some(&n),
        &r,
        Some(oriented::Input {
            document: &original,
            file_sha256: &hash,
        }),
    );
    let (p, s) = render_bound(&report);
    let off = &s.rows[0].off_axis;
    assert!(
        off.contains("12A")
            && off.contains("genomic attachment 100")
            && off.contains("outside the locus window")
    );
    let serialized = serde_json::to_value(&p).unwrap();
    let run = &serialized["reads"][0]["genomic_basecalls"][0]["variants"][0]["runs"][0];
    assert_eq!(run["query_interval"], serde_json::json!([0, 12]));
    assert_eq!(run["reference_anchor"], 100);
    assert!(tsv(&p).unwrap().contains("\"reference_anchor\":100"));

    let (mut a, n, mut r, mut original) = originals_fixture();
    r.flank_left = 0;
    path(&mut a, &[((0, 150), (100, 250))], vec![]);
    a.units[0].rows[0].axis_length = Some(150);
    a.units[0].rows[0].sequence_length = Some(150);
    original.reads[0].sequence = format!("{}{}{}", "C".repeat(65), "A".repeat(12), "C".repeat(73));
    let hash = pileup::hash(&original).unwrap();
    let report = run_originals(
        &a,
        Some(&n),
        &r,
        Some(oriented::Input {
            document: &original,
            file_sha256: &hash,
        }),
    );
    let (p, s) = render_bound(&report);
    let run = &p.reads[0].genomic_basecalls[0].variants[0].runs[0];
    assert!(!run.terminal && run.aligned_to_reference);
    assert_eq!(run.reference_anchor, Some(165));
    assert!(s.rows[0].off_axis.contains("genomic attachment 165"));
    assert!(s.rows[0].runs.is_empty());
}

#[test]
fn read_end_locus_exon_containment_and_one_base_overhang_cannot_decide_molecule_strand() {
    for extra in [1, 3] {
        let (mut a, mut n, r) = fixture();
        n.transcripts[1].genomic_exons[0] = (200, 405);
        path(&mut a, &[((0, 150 + extra), (250, 400 + extra))], vec![]);
        let report = run(&a, Some(&n), &r);
        assert_eq!(preferred(&report), ["opposite"]);
        let winner = candidate(&report, "opposite");
        assert_eq!(winner.score.discriminating_exon_bp, extra);
        assert_eq!(winner.class, Class::AntisenseCompatible);
        assert_eq!(report.reads[0].molecule_strand.strand, None);
        assert_eq!(report.reads[0].molecule_strand.status, "unresolved");
        assert!(
            report.reads[0]
                .molecule_strand
                .reasons
                .contains(&"exon_containment_only_no_junction_support".into())
        );
        assert!(
            winner.paths[0]
                .evidence
                .contains(&"exon_containment_only_no_junction_support".into())
        );
    }
    let (mut a, mut n, r) = fixture();
    n.transcripts[1] = model("opposite", '-', vec![(160, 190)]);
    path(&mut a, &[((0, 20), (165, 185))], vec![]);
    let report = run(&a, Some(&n), &r);
    assert_eq!(preferred(&report), ["opposite"]);
    assert_eq!(report.reads[0].molecule_strand.strand, None);
    assert!(
        report.reads[0]
            .molecule_strand
            .reasons
            .contains(&"exon_containment_only_no_junction_support".into())
    );
}

#[test]
fn read_end_locus_mirrored_identical_chains_stay_alternatives_with_an_unrelated_third_model() {
    for third in [false, true] {
        let (mut a, mut n, r) = fixture();
        n.transcripts[1].genomic_exons = n.transcripts[0].genomic_exons.clone();
        if third {
            n.transcripts.push(model("unrelated", '+', vec![(70, 90)]));
        }
        let report = run(&a, Some(&n), &r);
        for id in ["sense", "opposite"] {
            assert_eq!(candidate(&report, id).class, Class::UninformativeOverlap);
            assert_eq!(candidate(&report, id).score.discriminating_exon_bp, 0);
        }
        path(
            &mut a,
            &[((0, 50), (100, 150)), ((50, 100), (200, 250))],
            vec![(150, 200)],
        );
        let spliced = run(&a, Some(&n), &r);
        assert_eq!(
            preferred(&spliced).into_iter().collect::<BTreeSet<_>>(),
            BTreeSet::from(["sense", "opposite"])
        );
        assert_eq!(spliced.reads[0].molecule_strand.strand, None);
    }
}

#[test]
fn read_end_locus_full_chain_ignores_deletion_gaps_and_short_ends_but_not_missing_junctions() {
    let (mut a, mut n, mut r) = fixture();
    n.transcripts = vec![model(
        "sense",
        '+',
        vec![(100, 150), (200, 250), (300, 400)],
    )];
    r.exon = (300, 400);
    path(
        &mut a,
        &[
            ((0, 15), (110, 125)),
            ((15, 39), (126, 150)),
            ((39, 89), (200, 250)),
            ((89, 119), (300, 330)),
        ],
        vec![(150, 200), (250, 300)],
    );
    let report = run(&a, Some(&n), &r);
    let m = candidate(&report, "sense");
    assert!(m.score.full_chain);
    assert_eq!(m.class, Class::ExactTranscriptMatch);
    assert_eq!(m.score.endpoint_matches, 0);
    assert_eq!(m.score.endpoint_distance, 80);
    path(
        &mut a,
        &[((0, 30), (120, 150)), ((30, 60), (200, 230))],
        vec![(150, 200)],
    );
    let partial = run(&a, Some(&n), &r);
    assert!(!candidate(&partial, "sense").score.full_chain);
    assert_eq!(
        candidate(&partial, "sense").class,
        Class::CompatibleTranscripts
    );
}

#[test]
fn read_end_locus_upstream_unavailable_units_are_not_mislabelled_as_outside_locus() {
    let (mut a, n, r) = fixture();
    a.observations.clear();
    a.units[0].hypotheses.clear();
    a.units[0].unclipped_axis_length = None;
    a.units[0].unavailable = vec![ObservedSamUnavailableReason::QueryAxisUnsupported];
    let report = run(&a, Some(&n), &r);
    assert!(report.reads.is_empty());
    assert_eq!(
        report.samples[0].excluded[0].1,
        "upstream_unit_unavailable:query_axis_unsupported"
    );
    assert_eq!(
        unavailable_names(&[ObservedSamUnavailableReason::NotProvided]).unwrap(),
        ["not_provided"]
    );
}

#[test]
fn read_end_locus_query_disjoint_collinear_and_noncollinear_segments_have_distinct_chimera_states()
{
    for (second_target, chimeric) in [((350, 400), false), ((200, 250), true)] {
        let (mut a, n, r) = fixture();
        path(&mut a, &[((0, 50), (250, 300))], vec![]);
        a.observations[0].query_axis_length = 100;
        a.observations[0].alignment_clipped = [false, true];
        a.units[0].unclipped_axis_length = Some(100);
        let mut second = a.observations[0].clone();
        second.observation_id = "query_disjoint_supplement".into();
        second.source_ordinal = 2;
        second.aligned_query_interval = (50, 100);
        second.observed_span = second_target;
        second.alignment_clipped = [true, false];
        second.blocks = vec![SupportedBlock {
            query: (50, 100),
            target: second_target,
        }];
        a.units[0].hypotheses[0].push(second.observation_id.clone());
        a.observations.push(second);
        let report = run(&a, Some(&n), &r);
        assert_eq!(report.reads[0].hypotheses[0].chimeric, chimeric);
        assert_eq!(
            candidate(&report, "sense").class == Class::Chimeric,
            chimeric
        );
    }
}

#[test]
fn read_end_locus_positive_novel_junction_and_disagreeing_hypotheses_are_preserved() {
    let (mut a, mut n, mut r) = fixture();
    n.transcripts = vec![model(
        "sense",
        '+',
        vec![(100, 150), (200, 250), (300, 400)],
    )];
    r.exon = (300, 400);
    path(
        &mut a,
        &[((0, 25), (125, 150)), ((25, 75), (300, 350))],
        vec![(150, 300)],
    );
    let novel = run(&a, Some(&n), &r);
    assert_eq!(
        candidate(&novel, "sense").class,
        Class::NovelIsoformCompatible
    );
    assert!(
        candidate(&novel, "sense")
            .conflicts
            .contains(&"junction_not_in_model_chain".into())
    );
    let (mut a, n, r) = fixture();
    path(&mut a, &[((0, 100), (250, 350))], vec![]);
    let mut other = a.observations[0].clone();
    other.observation_id = "alternative_locus_path".into();
    other.source_ordinal = 2;
    other.blocks[0].target = (450, 550);
    a.units[0]
        .hypotheses
        .push(vec![other.observation_id.clone()]);
    a.observations.push(other);
    let mixed = run(&a, Some(&n), &r);
    let m = candidate(&mixed, "sense");
    assert_eq!(m.paths.len(), 2);
    assert_eq!(m.class, Class::Unresolved);
    assert_ne!(m.paths[0].class, m.paths[1].class);
}

fn cds_fixture() -> (
    ObservedSamAnalysis,
    PrimerPayloadAnnotation,
    Region,
    ModelDocument,
) {
    let (mut a, n, r) = fixture();
    path(
        &mut a,
        &[((0, 40), (110, 150)), ((40, 110), (200, 270))],
        vec![(150, 200)],
    );
    let m = TranscriptModel {
        transcript_id: "sense".into(),
        gene_id: "gene_sense".into(),
        contig: r.contig.clone(),
        strand: '+',
        genomic_exons_0based_half_open: n.transcripts[0].genomic_exons.clone(),
        genomic_cds_0based_half_open: vec![(120, 150), (200, 250)],
        genomic_stop_codon_0based_half_open: vec![(250, 253)],
    };
    let doc = ModelDocument {
        schema: MODEL_SCHEMA.into(),
        source_id: "synthetic".into(),
        assembly: a.reference.assembly.clone(),
        reference_fasta_sha256: a.reference.fasta_sha256.clone(),
        annotation_sha256: pileup::hash(&n).unwrap(),
        coordinate_convention: "zero_based_half_open".into(),
        model_source: ModelSource {
            format: "gtf".into(),
            name: "Hand-crafted synthetic CDS and stop geometry".into(),
            sha256: format!("sha256:{}", "b".repeat(64)),
        },
        transcripts: vec![m],
    };
    (a, n, r, doc)
}
#[test]
fn read_end_locus_cds_metrics_and_conflicting_paths_keep_distinct_utr_reasons() {
    let (mut a, n, r, doc) = cds_fixture();
    let report = build_input(
        &r,
        input("sample", &a, Some(&n), Some(&doc), None),
        Limits::default(),
    )
    .unwrap();
    let m = candidate(&report, "sense");
    assert_eq!(m.paths[0].utr5, Metric::Observed { value_bp: 10 });
    assert_eq!(m.paths[0].utr3, Metric::Observed { value_bp: 17 });
    let mut other = a.observations[0].clone();
    other.observation_id = "second_utr_path".into();
    other.source_ordinal = 2;
    other.blocks = vec![
        SupportedBlock {
            query: (0, 45),
            target: (105, 150),
        },
        SupportedBlock {
            query: (45, 110),
            target: (200, 265),
        },
    ];
    other.observed_span = (105, 265);
    other.skips = vec![(150, 200)];
    a.units[0]
        .hypotheses
        .push(vec![other.observation_id.clone()]);
    a.observations.push(other);
    let report = build_input(
        &r,
        input("sample", &a, Some(&n), Some(&doc), None),
        Limits::default(),
    )
    .unwrap();
    let p = page(&report, &View::default()).unwrap();
    let key = &candidate(&report, "sense").model.key;
    assert_eq!(
        p.summaries[0].utr_by_model[key].five.get("alternatives"),
        Some(&1)
    );
    assert_eq!(
        p.summaries[0].utr_by_model[key].three.get("alternatives"),
        Some(&1)
    );
    path(
        &mut a,
        &[((0, 40), (110, 150)), ((40, 110), (450, 520))],
        vec![(150, 450)],
    );
    a.units[0].hypotheses.truncate(1);
    a.observations.truncate(1);
    let bad = build_input(
        &r,
        input("sample", &a, Some(&n), Some(&doc), None),
        Limits::default(),
    )
    .unwrap();
    assert_eq!(
        candidate(&bad, "sense").paths[0].utr5,
        Metric::Missing {
            reason: Reason::ModelProjectionIncompatible
        }
    );
}

#[test]
fn read_end_locus_basecall_variants_count_runs_once_and_keep_alignment_alternatives() {
    let (mut a, n, r, mut original) = originals_fixture();
    let mut variant = original.reads[0].clone();
    variant.sequence.replace_range(30..31, "G");
    original.reads.push(variant);
    let mut other = a.observations[0].clone();
    other.observation_id = "alternate_coverage".into();
    other.source_ordinal = 2;
    other.aligned_query_interval = (0, 100);
    other.alignment_clipped = [false, true];
    other.clipping = Some(super::super::super::ObservedSamClipping {
        hard: [0, 0],
        soft: [0, 30],
        padding_operations: 0,
    });
    other.blocks = vec![
        SupportedBlock {
            query: (0, 50),
            target: (100, 150),
        },
        SupportedBlock {
            query: (50, 100),
            target: (200, 250),
        },
    ];
    a.units[0]
        .hypotheses
        .push(vec![other.observation_id.clone()]);
    a.observations.push(other);
    let hash = pileup::hash(&original).unwrap();
    let report = run_originals(
        &a,
        Some(&n),
        &r,
        Some(oriented::Input {
            document: &original,
            file_sha256: &hash,
        }),
    );
    let p = page(&report, &View::default()).unwrap();
    let c = &p.summaries[0];
    assert_eq!(
        (c.possible_terminal_a_runs, c.possible_terminal_t_runs),
        (1, 1)
    );
    assert_eq!(c.units_with_base_alternatives, 1);
    assert_eq!(c.run_alignment.get("alternatives"), Some(&1));
}

#[test]
fn read_end_locus_budgets_refuse_named_run_assignment_and_retained_geometry_expansion() {
    let (a, n, r, original) = originals_fixture();
    let hash = pileup::hash(&original).unwrap();
    let err = build_input(
        &r,
        input(
            "sample",
            &a,
            Some(&n),
            None,
            Some(oriented::Input {
                document: &original,
                file_sha256: &hash,
            }),
        ),
        Limits {
            max_observations: 3,
            ..Limits::default()
        },
    )
    .unwrap_err();
    assert!(err.message.ends_with("locus_run_budget"));
    let (a, n, r) = fixture();
    let err = build_input(
        &r,
        input("sample", &a, Some(&n), None, None),
        Limits {
            max_observations: 3,
            ..Limits::default()
        },
    )
    .unwrap_err();
    assert!(err.message.ends_with("locus_retained_geometry_budget"));
    let (mut a, mut n, mut r) = fixture();
    r.anchor_transcript_id = None;
    n.transcripts = vec![
        model("one", '+', vec![(200, 400)]),
        model("two", '-', vec![(200, 400)]),
    ];
    let mut other = a.observations[0].clone();
    other.observation_id = "second_hypothesis".into();
    other.source_ordinal = 2;
    a.units[0]
        .hypotheses
        .push(vec![other.observation_id.clone()]);
    a.observations.push(other);
    let err = build_input(
        &r,
        input("sample", &a, Some(&n), None, None),
        Limits {
            max_observations: 2,
            ..Limits::default()
        },
    )
    .unwrap_err();
    assert!(err.message.ends_with("locus_assignment_budget"));
}

#[test]
fn read_end_locus_two_samples_deduplicate_tracks_and_filter_without_rewriting_denominators() {
    let (a, n, r) = fixture();
    let mut b = a.clone();
    b.units[0].unit.source_sha256 = format!("sha256:{}", "b".repeat(64));
    b.observations[0].unit = b.units[0].unit.clone();
    let report = build(
        "synthetic",
        &r,
        &[
            input("first", &a, Some(&n), None, None),
            input("second", &b, Some(&n), None, None),
        ],
        oriented::Policy::default(),
        Limits::default(),
        &mut || Ok(()),
    )
    .unwrap();
    let p = page(&report, &View::default()).unwrap();
    let s = scene(&report, &p).unwrap();
    assert_eq!(s.top_tracks.len(), 1);
    assert_eq!(s.bottom_tracks.len(), 1);
    assert_eq!(s.rows.len(), 2);
    let details: serde_json::Value = serde_json::from_str(&s.top_tracks[0].details).unwrap();
    assert_eq!(details["model_keys"].as_array().unwrap().len(), 2);
    let selected = page(
        &report,
        &View {
            sample_ids: vec!["second".into()],
            ..View::default()
        },
    )
    .unwrap();
    assert_eq!(selected.reads.len(), 1);
    assert_eq!(selected.total_input_units, 2);
    assert_eq!(selected.total_locus_units, 2);
    assert_eq!(selected.summaries[0].visible_units, 0);
    assert_eq!(selected.summaries[1].visible_units, 1);
    let minus = report.models.iter().find(|m| m.strand == '-').unwrap();
    let p = page(
        &report,
        &View {
            projection: ViewProjection::Transcript {
                model_key: minus.key.clone(),
            },
            ..View::default()
        },
    )
    .unwrap();
    let s = scene(&report, &p).unwrap();
    assert_eq!(s.top_tracks.len(), 1);
    assert!(s.bottom_tracks.is_empty());
}

#[test]
fn read_end_locus_nonadmitted_only_overlap_has_a_separate_membership_count() {
    let (mut a, n, r) = fixture();
    path(&mut a, &[((0, 50), (700, 750))], vec![]);
    let mut other = a.observations[0].clone();
    other.observation_id = "nonadmitted_at_locus".into();
    other.on_optimal_path = false;
    other.blocks[0].target = (250, 300);
    a.observations.push(other);
    let report = run(&a, Some(&n), &r);
    assert_eq!(
        report.reads[0].membership,
        Membership::NonAdmittedObservationOnly
    );
    assert_eq!(report.samples[0].non_admitted_observation_only_units, 1);
    let p = page(&report, &View::default()).unwrap();
    assert_eq!(p.summaries[0].non_admitted_observation_only_units, 1);
    assert!(
        serde_json::to_string(&p)
            .unwrap()
            .contains("non-admitted-observation-only")
    );
}

#[test]
fn read_end_locus_scale_30_models_12_exons_20_deletion_split_reads_records_work() {
    let (mut a, mut n, mut r) = fixture();
    r.contig_length = 2000;
    r.exon = (100, 140);
    r.flank_left = 100;
    r.flank_right = 1000;
    r.anchor_transcript_id = Some("model_00".into());
    let exons: Vec<_> = (0..12).map(|i| (100 + i * 60, 140 + i * 60)).collect();
    n.transcripts = (0..30)
        .map(|i| {
            model(
                &format!("model_{i:02}"),
                if i < 10 { '-' } else { '+' },
                exons.clone(),
            )
        })
        .collect();
    for m in &mut n.transcripts {
        m.target_length = 2000;
    }
    let mut blocks = vec![];
    let mut q = 0;
    for (i, &(start, _)) in exons.iter().enumerate() {
        let chunks = if i < 4 {
            vec![(0, 10), (11, 25), (26, 32), (33, 40)]
        } else {
            vec![(0, 10), (11, 25), (26, 40)]
        };
        for (l, z) in chunks {
            let len = z - l;
            blocks.push(SupportedBlock {
                query: (q, q + len),
                target: (start + l, start + z),
            });
            q += len;
        }
    }
    assert_eq!(blocks.len(), 40);
    let template = a.observations[0].clone();
    let unit = a.units[0].clone();
    a.observations.clear();
    a.units.clear();
    for i in 0..20 {
        let mut u = unit.clone();
        u.unit.unit_id = format!("read_{i:02}");
        u.unclipped_axis_length = Some(q);
        let mut o = template.clone();
        o.observation_id = format!("obs_{i:02}");
        o.unit = u.unit.clone();
        o.source_ordinal = i + 1;
        o.blocks = blocks.clone();
        o.observed_span = (exons.first().unwrap().0, exons.last().unwrap().1);
        o.aligned_query_interval = (0, q);
        o.alignment_clipped = [false, false];
        o.query_axis_length = q;
        o.skips = exons.windows(2).map(|p| (p[0].1, p[1].0)).collect();
        u.hypotheses = vec![vec![o.observation_id.clone()]];
        a.units.push(u);
        a.observations.push(o);
    }
    let report = build_input(
        &r,
        input("scale", &a, Some(&n), None, None),
        Limits {
            max_work: 100_000_000,
            ..Limits::default()
        },
    )
    .unwrap();
    assert_eq!(report.reads.len(), 20);
    assert_eq!(report.models.len(), 30);
    assert!(report.work > 0 && report.work <= 100_000_000);
    assert!(
        report
            .reads
            .iter()
            .flat_map(|r| &r.assignments)
            .all(|a| a.score.discriminating_exon_bp == 0)
    );
    eprintln!("read_end_locus_scale work={}", report.work);
}
