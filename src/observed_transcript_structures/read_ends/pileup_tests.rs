//! Synthetic saved-evidence controls; fixtures share the independently stated geometry corpus.

use super::super::tests::{fixture, normalized};
use super::*;
use evidence::{CompositionPolicy, Provenance};

fn report(case: &str, extra_hypothesis: bool) -> Report {
    let (document, mut analysis) = fixture(case);
    if extra_hypothesis {
        let h = analysis.units[0].hypotheses[0].clone();
        analysis.units[0].hypotheses.push(h);
    }
    let annotation = normalized(&document, &analysis);
    let projection =
        project_read_ends(&document, &analysis, &annotation, Limits::default()).unwrap();
    let hash = format!("sha256:{}", "a".repeat(64));
    let bindings = InputBindings {
        report_content_sha256: hash.clone(),
        data_sha256: hash.clone(),
        model_sha256: hash.clone(),
        projection_sha256: hash,
    };
    let composition = evidence::unavailable_composition(CompositionPolicy {
        orientation: Provenance::Undeclared,
        trimming: Provenance::Undeclared,
        adapter_sequences: vec![],
        primer_sequences: vec![],
        window_bp: 20,
        minimum_run_bp: 10,
        maximum_mismatches: 1,
    })
    .unwrap();
    let transcript = if case == "minus" {
        "tx_minus"
    } else {
        "tx_plus"
    };
    build(
        &document.source_id,
        "GENE_SYN_A",
        transcript,
        &[0.0, 0.5, 1.0],
        &[Input {
            sample: SampleInput {
                sample_id: "synthetic_sample",
                gene_id: "GENE_SYN_A",
                bindings: &bindings,
                analysis: &analysis,
                projection: &projection,
            },
            model: &document,
            annotation: &annotation,
            composition: &composition,
            tails: &[],
            tail_file_sha256: None,
            original_bases: None,
            terminal_policy: oriented::Policy::default(),
        }],
        Limits::default(),
    )
    .unwrap()
}

fn view() -> View {
    View {
        schema: VIEW_SCHEMA.into(),
        sort: SortKey::UnitKey,
        descending: false,
        filters: vec![],
        offset: 0,
        count: 10,
        max_sublanes: 10,
        width: 1400,
        scale: "genomic".into(),
        channels: Channels {
            utr5: true,
            utr3: true,
            read_length: true,
            clips: true,
            composition: true,
            tail: true,
        },
    }
}

#[test]
fn read_end_pileup_stacks_hypotheses_without_changing_denominators() {
    let report = report("plus", true);
    assert_eq!(report.lanes.len(), 1);
    assert_eq!(report.lanes[0].hypotheses.len(), 2);
    let p = page(&report, &view()).unwrap();
    let scene = scene(&report, &p).unwrap();
    assert_eq!(scene.rows.len(), 2);
    assert_eq!(p.coverage.samples[0].gene_assigned_units, 1);
    assert_eq!(
        p.orientation_summary.as_ref().unwrap().gene_assigned_units,
        1
    );
    assert_eq!(
        p.coverage.samples[0].selected_transcript_compatible_units,
        1
    );
    assert!(scene.rows[0].off_axis.contains("no_bound_basecalls"));
}

#[test]
fn read_end_oriented_legacy_report_hash_and_missing_evidence_remain_unchanged() {
    let mut report = report("plus", false);
    report.orientation_summary = None;
    assert!(page(&report, &view()).is_err());
    report.schema = LEGACY_SCHEMA.into();
    report.orientation_summary = None;
    for h in &mut report.lanes[0].hypotheses {
        h.oriented = None;
    }
    let bytes = serde_json::to_vec(&report).unwrap();
    let decoded: Report = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(bytes, serde_json::to_vec(&decoded).unwrap());
    assert_eq!(hash(&report).unwrap(), hash(&decoded).unwrap());
    let page = page(&decoded, &view()).unwrap();
    assert!(page.orientation_summary.is_none());
    assert!(
        scene(&decoded, &page)
            .unwrap()
            .rows
            .iter()
            .all(|r| r.runs.is_empty())
    );
}

#[test]
fn read_end_oriented_summary_distinguishes_unavailable_from_bound_nonterminal_runs() {
    let mut report = report("plus", true);
    assert_eq!(
        report
            .orientation_summary
            .as_ref()
            .unwrap()
            .units_with_bound_bases,
        0
    );
    for h in &mut report.lanes[0].hypotheses {
        let e = h.oriented.as_mut().unwrap();
        e.availability = "available".into();
        e.variants = vec![oriented::Variant {
            read_id: "synthetic".into(),
            sequence_sha256: format!("sha256:{}", "a".repeat(64)),
            query_length: 152,
            runs: vec![oriented::Run {
                base: 'A',
                length: 12,
                query_interval: (70, 82),
                terminal: false,
                read_ends: vec![],
                transcript_ends: None,
                aligned_to_reference: true,
                aligned_query_intervals: vec![(70, 82)],
                reference_intervals: vec![(100, 112)],
                reference_anchor: Some(106),
            }],
        }];
    }
    let summary = oriented::summarize(&report.lanes);
    assert_eq!(summary.units_with_bound_bases, 1);
    assert_eq!(summary.units_with_terminal_a, 0);
    assert_eq!(summary.gene_assigned_units, 1);
}

#[test]
fn read_end_oriented_scene_separates_directions_but_filtering_never_changes_denominators() {
    let mut report = report("plus", false);
    let mut opposite = report.lanes[0].clone();
    opposite.unit.unit.unit_id = "opposite".into();
    opposite.hypotheses[0].oriented.as_mut().unwrap().direction = oriented::Direction::Opposite;
    report.lanes.insert(0, opposite);
    report.orientation_summary = Some(oriented::summarize(&report.lanes));
    let page_all = page(&report, &view()).unwrap();
    let drawn = scene(&report, &page_all).unwrap();
    assert!(drawn.rows[0].orientation == Some(Orientation::Same));
    assert!(drawn.rows[1].orientation == Some(Orientation::Opposite));
    let mut filtered = view();
    filtered
        .filters
        .push(Filter::Composition { available: true });
    let page_filtered = page(&report, &filtered).unwrap();
    assert!(page_filtered.lanes.is_empty());
    assert_eq!(
        page_filtered.total_gene_assigned,
        page_all.total_gene_assigned
    );
    assert_eq!(
        json(&page_filtered.orientation_summary).unwrap(),
        json(&page_all.orientation_summary).unwrap()
    );
    assert_eq!(
        json(&page_filtered.coverage).unwrap(),
        json(&page_all.coverage).unwrap()
    );
}

#[test]
fn read_end_oriented_conflicting_hypotheses_stay_in_alternatives_without_double_counting() {
    let mut report = report("plus", true);
    report.lanes[0].hypotheses[1]
        .oriented
        .as_mut()
        .unwrap()
        .direction = oriented::Direction::Opposite;
    report.orientation_summary = Some(oriented::summarize(&report.lanes));
    let page = page(&report, &view()).unwrap();
    assert!(
        scene(&report, &page)
            .unwrap()
            .rows
            .iter()
            .all(|r| r.orientation == Some(Orientation::Alternatives))
    );
    assert_eq!(page.total_gene_assigned, 1);
    assert_eq!(
        page.orientation_summary
            .as_ref()
            .unwrap()
            .alternative_direction_units,
        1
    );
}

#[test]
fn read_end_pileup_filters_keep_coverage_and_use_explicit_ranges() {
    let report = report("plus", false);
    let mut v = view();
    v.filters.push(Filter::Composition { available: true });
    let p = page(&report, &v).unwrap();
    assert!(p.lanes.is_empty());
    assert_eq!(p.total_gene_assigned, 1);
    let filtered_scene = scene(&report, &p).unwrap();
    assert_eq!(
        (filtered_scene.start, filtered_scene.end),
        (
            report.display_interval.0 as u64,
            report.display_interval.1 as u64
        )
    );
    assert_eq!(p.coverage.samples[0].gene_assigned_units, 1);
    v.offset = 1;
    assert!(page(&report, &v).is_err());
    v = view();
    v.max_sublanes = 1;
    assert!(
        page(&self::report("plus", true), &v)
            .unwrap_err()
            .message
            .contains("required=2")
    );
    v.scale = "compressed".into();
    assert!(page(&report, &v).is_err());
}

#[test]
fn read_end_pileup_sorts_stably_and_never_drops_alternative_hypotheses() {
    let mut report = report("plus", true);
    let mut second = report.lanes[0].clone();
    second.unit.unit.unit_id = "a_synthetic".into();
    report.lanes.push(second);
    let mut v = view();
    v.sort = SortKey::ReadLength;
    v.count = 1;
    let p = page(&report, &v).unwrap();
    assert_eq!(p.lanes[0].unit.unit.unit_id, "a_synthetic");
    assert_eq!(p.lanes[0].hypotheses.len(), 2);
    assert_eq!(p.total_matching, 2);
    v.offset = 1;
    assert_ne!(
        p.lanes[0].unit.unit,
        page(&report, &v).unwrap().lanes[0].unit.unit
    );
}

#[test]
fn read_end_pileup_connectors_keep_n_deletion_and_segment_links_distinct() {
    let (_, mut analysis) = fixture("plus");
    let observation = &mut analysis.observations[0];
    observation.blocks = vec![
        SupportedBlock {
            query: (0, 10),
            target: (100, 110),
        },
        SupportedBlock {
            query: (10, 20),
            target: (130, 140),
        },
    ];
    observation.skips = vec![(115, 125)];
    let connectors = bands(observation);
    assert!(
        connectors
            .iter()
            .any(|c| c.kind == ConnectorKind::SpliceN && c.interval == (115, 125))
    );
    assert!(
        connectors
            .iter()
            .any(|c| c.kind == ConnectorKind::DeletionOrGap && c.interval == (110, 115))
    );
    assert!(
        connectors
            .iter()
            .any(|c| c.kind == ConnectorKind::DeletionOrGap && c.interval == (125, 130))
    );
    let mut report = report("plus", false);
    report.lanes[0].hypotheses[0].segment_links.push(Connector {
        interval: (110, 130),
        kind: ConnectorKind::SegmentLink,
    });
    let p = page(&report, &view()).unwrap();
    assert!(
        scene(&report, &p).unwrap().rows[0]
            .bands
            .iter()
            .any(|b| b.kind == "segment_link")
    );
}

#[test]
fn read_end_pileup_svg_and_tsv_are_deterministic_and_preserve_provenance() {
    for name in ["plus", "minus"] {
        let report = report(name, false);
        let p = page(&report, &view()).unwrap();
        assert_eq!(tsv(&p).unwrap(), tsv(&p).unwrap());
        let scene = scene(&report, &p).unwrap();
        let svg = gentle_render::read_end_pileup::render(&scene).unwrap();
        assert_eq!(svg, gentle_render::read_end_pileup::render(&scene).unwrap());
        assert!(tsv(&p).unwrap().contains("annotation_gene_ids"));
        assert!(tsv(&p).unwrap().contains("model_gene_ids"));
        assert!(svg.contains(&p.view_sha256));
        let mut tampered = p.clone();
        tampered.view.count += 1;
        assert!(self::scene(&report, &tampered).is_err());
    }
}

#[test]
fn read_end_pileup_selected_geometry_is_checked_even_without_compatible_reads() {
    let (document, analysis) = fixture("plus");
    let annotation = normalized(&document, &analysis);
    let mut meter = Meter {
        limits: Limits::default(),
        work: 0,
        control: None,
    };
    let index = projection::annotation_index(&document, &annotation, &mut meter).unwrap();
    let mut model = document.transcripts[0].clone();
    model.genomic_exons_0based_half_open[0].0 -= 1;
    let prepared = prepare(&model, &mut meter).unwrap();
    assert!(
        projection::checked_annotation(&model.transcript_id, Some(&prepared), &index, &mut meter)
            .is_err()
    );
}
