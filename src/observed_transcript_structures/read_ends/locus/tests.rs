//! Hand-crafted general sense/antisense loci, not a biological or calibrated strand oracle.

use super::*;
use crate::engine::PrimerPayloadTranscript;

pub(super) fn fixture() -> (ObservedSamAnalysis, PrimerPayloadAnnotation, Region) {
    let (document, mut analysis) = super::super::tests::fixture("plus");
    let mut annotation = super::super::tests::normalized(&document, &analysis);
    annotation.transcripts = vec![
        model("sense", '+', vec![(100, 150), (200, 400)]),
        model("opposite", '-', vec![(200, 400), (500, 550)]),
    ];
    analysis.observations.truncate(1);
    analysis.units.truncate(1);
    let o = &mut analysis.observations[0];
    o.blocks = vec![SupportedBlock {
        query: (0, 100),
        target: (250, 350),
    }];
    o.skips.clear();
    o.query_axis_length = 100;
    o.query_to_reference_strand = '+';
    o.observed_span = (250, 350);
    o.aligned_query_interval = (0, 100);
    o.alignment_clipped = [false, false];
    o.on_optimal_path = true;
    o.clipping = None;
    analysis.units[0].unclipped_axis_length = Some(100);
    analysis.units[0].hypotheses = vec![vec![o.observation_id.clone()]];
    analysis.units[0].rows.clear();
    let region = Region {
        contig: "synthetic_contig".into(),
        contig_length: 1000,
        exon: (200, 400),
        flank_left: 150,
        flank_right: 200,
        anchor_transcript_id: Some("sense".into()),
    };
    (analysis, annotation, region)
}
pub(super) fn model(id: &str, strand: char, exons: Vec<Interval>) -> PrimerPayloadTranscript {
    PrimerPayloadTranscript {
        transcript_id: id.into(),
        gene_ids: vec![format!("gene_{id}")],
        target_name: "synthetic_contig".into(),
        target_length: 1000,
        genomic_contig: Some("synthetic_contig".into()),
        genomic_strand: Some(strand),
        genomic_exons: exons,
    }
}
pub(super) fn run(
    analysis: &ObservedSamAnalysis,
    annotation: Option<&PrimerPayloadAnnotation>,
    region: &Region,
) -> Report {
    run_originals(analysis, annotation, region, None)
}
pub(super) fn run_originals(
    analysis: &ObservedSamAnalysis,
    annotation: Option<&PrimerPayloadAnnotation>,
    region: &Region,
    originals: Option<oriented::Input<'_>>,
) -> Report {
    let input = Input {
        sample_id: "sample",
        bindings: Bindings {
            report_content_sha256: format!("sha256:{}", "1".repeat(64)),
            data_sha256: format!("sha256:{}", "2".repeat(64)),
            annotation_sha256: annotation.map(|n| pileup::hash(n).unwrap()),
            model_sha256: None,
        },
        analysis,
        annotation,
        model: None,
        original_bases: originals,
        diagnostics: vec![],
    };
    build(
        "synthetic",
        region,
        &[input],
        oriented::Policy::default(),
        Limits::default(),
        &mut || Ok(()),
    )
    .unwrap()
}
pub(super) fn path(
    a: &mut ObservedSamAnalysis,
    blocks: &[(Interval, Interval)],
    skips: Vec<Interval>,
) {
    let o = &mut a.observations[0];
    o.blocks = blocks
        .iter()
        .map(|(query, target)| SupportedBlock {
            query: *query,
            target: *target,
        })
        .collect();
    o.skips = skips;
    let n = blocks.iter().map(|b| b.0.1).max().unwrap();
    o.observed_span = (
        blocks.iter().map(|b| b.1.0).min().unwrap(),
        blocks.iter().map(|b| b.1.1).max().unwrap(),
    );
    o.aligned_query_interval = (blocks.iter().map(|b| b.0.0).min().unwrap(), n);
    o.alignment_clipped = [o.aligned_query_interval.0 > 0, false];
    o.query_axis_length = n;
    a.units[0].unclipped_axis_length = Some(n);
}
pub(super) fn preferred(report: &Report) -> Vec<&str> {
    report.reads[0]
        .assignments
        .iter()
        .filter(|a| a.preferred)
        .map(|a| a.model.transcript_id.as_str())
        .collect()
}
#[test]
fn read_end_locus_shared_exon_is_uninformative_for_both_bam_directions() {
    let (mut a, n, r) = fixture();
    for strand in ['+', '-'] {
        a.observations[0].query_to_reference_strand = strand;
        let report = run(&a, Some(&n), &r);
        assert_eq!(preferred(&report).len(), 2);
        assert!(
            report.reads[0]
                .assignments
                .iter()
                .all(|m| m.class == Class::UninformativeOverlap)
        );
        assert_eq!(report.reads[0].molecule_strand.strand, None);
        assert!(!report.reads[0].search_declared_complete);
        assert!(
            report.reads[0]
                .molecule_strand
                .reasons
                .iter()
                .any(|r| r.contains("incomplete_other_loci_not_excluded"))
        );
    }
}
#[test]
fn read_end_locus_specific_sense_and_antisense_junctions_outweigh_endpoints() {
    for (blocks, skip, winner) in [
        (
            vec![((0, 50), (100, 150)), ((50, 150), (200, 300))],
            (150, 200),
            "sense",
        ),
        (
            vec![((0, 100), (300, 400)), ((100, 150), (500, 550))],
            (400, 500),
            "opposite",
        ),
    ] {
        let (mut a, n, r) = fixture();
        path(&mut a, &blocks, vec![skip]);
        let first = run(&a, Some(&n), &r);
        assert_eq!(preferred(&first), [winner]);
        assert_eq!(
            first.reads[0].assignments[0].score.discriminating_exon_bp,
            50
        );
        assert_eq!(first.reads[0].molecule_strand.status, "model_supported");
        assert_eq!(
            first.reads[0].molecule_strand.strand,
            Some(if winner == "sense" { '+' } else { '-' })
        );
        let len = a.units[0].unclipped_axis_length.unwrap();
        for o in &mut a.observations {
            o.query_to_reference_strand = '-';
            for b in &mut o.blocks {
                b.query = (len - b.query.1, len - b.query.0);
            }
        }
        let reverse = run(&a, Some(&n), &r);
        assert_eq!(preferred(&reverse), [winner]);
        assert_eq!(
            reverse.reads[0].molecule_strand.strand,
            first.reads[0].molecule_strand.strand
        );
        for (x, y) in first.reads[0]
            .assignments
            .iter()
            .zip(&reverse.reads[0].assignments)
        {
            assert_eq!(x.score, y.score);
            assert_eq!(x.class, y.class);
        }
        if winner == "opposite" {
            assert_eq!(
                first.reads[0].assignments[0].class,
                Class::AntisenseCompatible
            );
        }
    }
}
#[test]
fn read_end_locus_identical_models_remain_ordered_alternatives_without_duplicate_reads() {
    let (a, mut n, r) = fixture();
    let mut duplicate = n.transcripts[0].clone();
    duplicate.transcript_id = "sense_copy".into();
    n.transcripts.push(duplicate);
    let report = run(&a, Some(&n), &r);
    assert_eq!(report.reads.len(), 1);
    assert_eq!(report.samples[0].locus_units, 1);
    let names: Vec<_> = report.reads[0]
        .assignments
        .iter()
        .map(|a| a.model.transcript_id.as_str())
        .collect();
    assert_eq!(names.len(), 3);
    assert!(
        names.iter().position(|n| *n == "sense").unwrap()
            < names.iter().position(|n| *n == "sense_copy").unwrap()
    );
    assert_eq!(preferred(&report).len(), 3);
}
#[test]
fn read_end_locus_novel_junction_is_not_forced_to_nearest_model() {
    let (mut a, n, r) = fixture();
    path(
        &mut a,
        &[((0, 25), (100, 125)), ((25, 75), (225, 275))],
        vec![(125, 225)],
    );
    let report = run(&a, Some(&n), &r);
    assert!(preferred(&report).is_empty());
    assert!(
        report.reads[0]
            .assignments
            .iter()
            .all(|a| a.class == Class::Unresolved)
    );
}
#[test]
fn read_end_locus_projection_switch_preserves_geometry_ids_and_denominators() {
    let (a, n, r) = fixture();
    let report = run(&a, Some(&n), &r);
    let genome = page(&report, &View::default()).unwrap();
    for model in &report.models {
        let view = View {
            projection: ViewProjection::Transcript {
                model_key: model.key.clone(),
            },
            ..View::default()
        };
        let p = page(&report, &view).unwrap();
        assert_eq!(p.total_input_units, genome.total_input_units);
        assert_eq!(p.reads[0].read_key, genome.reads[0].read_key);
        assert_eq!(
            serde_json::to_value(&p.reads[0].observations).unwrap(),
            serde_json::to_value(&genome.reads[0].observations).unwrap()
        );
        let s = scene(&report, &p).unwrap();
        assert_eq!(s.rows.len(), 1);
        assert_eq!(s.axis_kind.as_deref(), Some("transcript"));
    }
    let filtered = page(
        &report,
        &View {
            class: Some(Class::Chimeric),
            ..View::default()
        },
    )
    .unwrap();
    assert_eq!(filtered.total_matching, 0);
    assert_eq!(filtered.summaries[0].input_units, 1);
    assert_eq!(filtered.summaries[0].locus_units, 1);
}
#[test]
fn read_end_locus_mirrored_models_preserve_scores_and_boundary_coordinates() {
    let (mut a, mut n, mut r) = fixture();
    path(
        &mut a,
        &[((0, 50), (100, 150)), ((50, 100), (200, 250))],
        vec![(150, 200)],
    );
    let first = run(&a, Some(&n), &r);
    for m in &mut n.transcripts {
        m.genomic_strand = Some(if m.genomic_strand == Some('+') {
            '-'
        } else {
            '+'
        });
        m.genomic_exons = m
            .genomic_exons
            .iter()
            .rev()
            .map(|&(a, z)| (1000 - z, 1000 - a))
            .collect();
    }
    for o in &mut a.observations {
        o.query_to_reference_strand = '-';
        for b in &mut o.blocks {
            b.target = (1000 - b.target.1, 1000 - b.target.0);
        }
        o.skips = o
            .skips
            .iter()
            .rev()
            .map(|&(a, z)| (1000 - z, 1000 - a))
            .collect();
    }
    r.exon = (600, 800);
    std::mem::swap(&mut r.flank_left, &mut r.flank_right);
    let mirror = run(&a, Some(&n), &r);
    for (x, y) in first.reads[0]
        .assignments
        .iter()
        .zip(&mirror.reads[0].assignments)
    {
        assert_eq!(x.score, y.score);
        assert_eq!(x.class, y.class);
    }
}
#[test]
fn read_end_locus_missing_changed_annotation_and_chimeras_are_explicit() {
    let (mut a, mut n, r) = fixture();
    let absent = run(&a, None, &r);
    assert_eq!(absent.reads.len(), 1);
    assert!(absent.models.is_empty());
    assert_eq!(
        scene(&absent, &page(&absent, &View::default()).unwrap())
            .unwrap()
            .rows
            .len(),
        1
    );
    let old = run(&a, Some(&n), &r);
    n.release = "new_snapshot".into();
    n.transcripts[0].genomic_exons[0] = (100, 140);
    let changed = run(&a, Some(&n), &r);
    assert_ne!(pileup::hash(&old).unwrap(), pileup::hash(&changed).unwrap());
    let mut second = a.observations[0].clone();
    second.observation_id = "chimera_second".into();
    second.source_ordinal = 2;
    second.contig = "other_contig".into();
    second.blocks[0].query = (100, 200);
    second.query_axis_length = 200;
    a.observations[0].query_axis_length = 200;
    a.units[0].unclipped_axis_length = Some(200);
    a.units[0].hypotheses[0].push(second.observation_id.clone());
    a.observations.push(second);
    let chimera = run(&a, Some(&n), &r);
    assert!(
        chimera.reads[0]
            .assignments
            .iter()
            .all(|a| a.class == Class::Chimeric)
    );
    assert!(preferred(&chimera).is_empty());
}
#[test]
fn read_end_locus_deterministic_exports_retain_every_candidate_rule_and_binding() {
    let (a, n, r) = fixture();
    let report = run(&a, Some(&n), &r);
    let p = page(&report, &View::default()).unwrap();
    let s = scene(&report, &p).unwrap();
    assert_eq!(s.top_tracks.len(), 1);
    assert_eq!(s.bottom_tracks.len(), 1);
    assert_eq!(s.rows.len(), 1);
    assert_eq!(
        tsv(&p).unwrap(),
        tsv(&page(&report, &View::default()).unwrap()).unwrap()
    );
    let svg = gentle_render::read_end_pileup::render(&s).unwrap();
    assert_eq!(svg, gentle_render::read_end_pileup::render(&s).unwrap());
    let table = tsv(&p).unwrap();
    assert!(
        table.contains(RULE_VERSION)
            && table.contains("structural_conflicts")
            && table.contains("annotation_sha256")
    );
}

#[test]
fn read_end_locus_no_legal_path_preserves_unassessed_geometry_without_assignments() {
    let (mut a, n, r) = fixture();
    a.units[0].hypotheses.clear();
    let report = run(&a, Some(&n), &r);
    let read = &report.reads[0];
    assert!(read.assignments.is_empty());
    assert_eq!(read.observations.len(), 1);
    assert_eq!(
        read.genomic_basecalls[0].availability,
        "no_legal_hypothesis_reference_coverage_unavailable"
    );
    let mut view = View::default();
    view.projection = ViewProjection::Transcript {
        model_key: report.models[0].key.clone(),
    };
    let p = page(&report, &view).unwrap();
    let scene = scene(&report, &p).unwrap();
    assert_eq!(
        scene.rows[0].orientation,
        Some(gentle_render::read_end_pileup::Orientation::Alternatives)
    );
    assert!(!scene.rows[0].sublanes[0].is_empty());
}

#[test]
fn read_end_locus_connections_outside_window_remain_marked_and_untrimmed() {
    let (mut a, n, mut r) = fixture();
    path(
        &mut a,
        &[((0, 50), (100, 150)), ((50, 100), (200, 250))],
        vec![(150, 200)],
    );
    r.flank_left = 0;
    r.flank_right = 0;
    let report = run(&a, Some(&n), &r);
    let p = page(&report, &View::default()).unwrap();
    let s = scene(&report, &p).unwrap();
    assert!(
        s.rows[0].sublanes[0]
            .iter()
            .any(|b| b.kind == "continuation")
    );
    assert_eq!(p.reads[0].observations[0].blocks[0].target, (100, 150));
    assert_eq!(p.reads[0].observations[0].skips, [(150, 200)]);
}

pub(super) fn originals_fixture() -> (
    ObservedSamAnalysis,
    PrimerPayloadAnnotation,
    Region,
    oriented::Document,
) {
    let (mut a, n, r) = fixture();
    path(
        &mut a,
        &[((12, 62), (100, 150)), ((62, 112), (200, 250))],
        vec![(150, 200)],
    );
    a.observations[0].query_axis_length = 130;
    a.observations[0].alignment_clipped = [true, true];
    a.units[0].unclipped_axis_length = Some(130);
    a.observations[0].clipping = Some(super::super::super::ObservedSamClipping {
        hard: [0, 0],
        soft: [12, 18],
        padding_operations: 0,
    });
    a.units[0].rows = vec![super::super::super::ObservedSamRowAssessment {
        ordinal: 1,
        qname: "synthetic_original".into(),
        read_group: None,
        flag: 0x800,
        target: Some("synthetic_contig".into()),
        position_1based: 101,
        reported_mapq: None,
        cigar: None,
        sequence_length: Some(130),
        quality_length: None,
        axis_length: Some(130),
        clipping: None,
        unavailable: None,
    }];
    let original = oriented::Document {
        schema: oriented::BASES_SCHEMA.into(),
        analysis_context_sha256: a.context_sha256.clone(),
        orientation: "original_sequencing_5_to_3".into(),
        provenance: "Hand-crafted synthetic original, not SAM reconstruction".into(),
        reads: vec![oriented::OriginalRead {
            unit: a.units[0].unit.clone(),
            read_id: "synthetic_original".into(),
            sequence: format!("{}{}{}", "A".repeat(12), "C".repeat(100), "T".repeat(18)),
        }],
    };
    (a, n, r, original)
}
#[test]
fn read_end_locus_original_softclipped_at_and_rc_cannot_overrule_specific_junctions() {
    let (mut a, n, r, mut original) = originals_fixture();
    let hash = pileup::hash(&original).unwrap();
    let first = run_originals(
        &a,
        Some(&n),
        &r,
        Some(oriented::Input {
            document: &original,
            file_sha256: &hash,
        }),
    );
    assert_eq!(preferred(&first), ["sense"]);
    assert!(first.reads[0].supplementary_present);
    assert!(!first.reads[0].hypotheses[0].chimeric);
    let runs = &first.reads[0].genomic_basecalls[0].variants[0].runs;
    assert_eq!(
        runs.iter().map(|r| (r.base, r.length)).collect::<Vec<_>>(),
        [('A', 12), ('T', 18)]
    );
    assert!(
        runs.iter()
            .all(|r| !r.aligned_to_reference && r.aligned_query_intervals.is_empty())
    );
    for b in &mut a.observations[0].blocks {
        b.query = (130 - b.query.1, 130 - b.query.0);
    }
    a.observations[0].query_to_reference_strand = '-';
    a.units[0].rows[0].flag |= 16;
    let flipped = run_originals(
        &a,
        Some(&n),
        &r,
        Some(oriented::Input {
            document: &original,
            file_sha256: &hash,
        }),
    );
    assert_eq!(preferred(&flipped), ["sense"]);
    assert_eq!(
        flipped.reads[0].genomic_basecalls[0].variants[0]
            .runs
            .iter()
            .map(|r| (r.base, r.length))
            .collect::<Vec<_>>(),
        [('A', 12), ('T', 18)]
    );
    original.reads[0].sequence = original.reads[0]
        .sequence
        .chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            _ => 'C',
        })
        .collect();
    let hash = pileup::hash(&original).unwrap();
    let rc = run_originals(
        &a,
        Some(&n),
        &r,
        Some(oriented::Input {
            document: &original,
            file_sha256: &hash,
        }),
    );
    assert_eq!(preferred(&rc), ["sense"]);
    for (x, y) in first.reads[0]
        .assignments
        .iter()
        .zip(&rc.reads[0].assignments)
    {
        assert_eq!(x.score, y.score);
        assert_eq!(x.class, y.class);
    }
    assert_eq!(
        rc.reads[0].genomic_basecalls[0].variants[0]
            .runs
            .iter()
            .map(|r| (r.base, r.length))
            .collect::<Vec<_>>(),
        [('A', 18), ('T', 12)]
    );
}

#[test]
fn read_end_locus_duplicate_units_and_work_exhaustion_refuse_without_truncation() {
    let (a, n, r) = fixture();
    let make = |sample: &'static str| Input {
        sample_id: sample,
        bindings: Bindings {
            report_content_sha256: format!("sha256:{}", "1".repeat(64)),
            data_sha256: format!("sha256:{}", "2".repeat(64)),
            annotation_sha256: Some(pileup::hash(&n).unwrap()),
            model_sha256: None,
        },
        analysis: &a,
        annotation: Some(&n),
        model: None,
        original_bases: None,
        diagnostics: vec![],
    };
    let error = build(
        "synthetic",
        &r,
        &[make("first"), make("second")],
        oriented::Policy::default(),
        Limits::default(),
        &mut || Ok(()),
    )
    .unwrap_err();
    assert!(error.message.contains("locus_duplicate_source_unit"));
    let error = build(
        "synthetic",
        &r,
        &[make("first")],
        oriented::Policy::default(),
        Limits {
            max_work: 0,
            ..Limits::default()
        },
        &mut || Ok(()),
    )
    .unwrap_err();
    assert!(error.message.contains("work_budget"));
    let cancelled = build(
        "synthetic",
        &r,
        &[make("first")],
        oriented::Policy::default(),
        Limits::default(),
        &mut || Err(fail("synthetic_cancel")),
    )
    .unwrap_err();
    assert_eq!(cancelled.message, "rna_read_end_synthetic_cancel");
}

#[test]
fn read_end_locus_nonadmitted_observations_remain_visible_and_consume_sublane_budget() {
    let (mut a, n, r) = fixture();
    let mut other = a.observations[0].clone();
    other.observation_id = "nonadmitted".into();
    other.on_optimal_path = false;
    a.observations.push(other);
    let report = run(&a, Some(&n), &r);
    let p = page(&report, &View::default()).unwrap();
    let scene = scene(&report, &p).unwrap();
    assert_eq!(scene.rows.len(), 1);
    assert_eq!(scene.rows[0].sublanes.len(), 2);
    assert_eq!(scene.rows[0].sublanes[1][0].kind, "unassessed");
    assert!(
        page(
            &report,
            &View {
                max_sublanes: 1,
                ..View::default()
            }
        )
        .is_err()
    );
}
