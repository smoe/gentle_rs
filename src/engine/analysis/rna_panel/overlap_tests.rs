//! Hand-crafted synthetic exon/anchor fixtures, recreated from tests::bases(778, 600).
//! No human genes or private reads. Exercises builder, screening, inspection and summary paths.

use super::tests::{bases, request, starter};
use super::*;
use serde_json::json;

#[derive(Clone)]
struct Source {
    start: usize,
    length: usize,
    extraction: Option<char>,
    genomic_strand: char,
    exons: Vec<(usize, usize)>,
}

fn plus_source() -> Source {
    Source {
        start: 1001,
        length: 150,
        extraction: Some('+'),
        genomic_strand: '+',
        exons: vec![(30, 60), (90, 120)],
    }
}

fn antisense_source(extraction: char) -> Source {
    Source {
        start: 1021,
        length: 180,
        extraction: Some(extraction),
        genomic_strand: '-',
        exons: vec![(10, 40), (70, 100)],
    }
}

fn sources(directory: &Path, specs: [Source; 2]) -> (GentleEngine, BuildRnaReadPanelRequest) {
    let (mut engine, build) = starter(2, directory);
    let reference = bases(778, 600);
    for (i, source) in specs.into_iter().enumerate() {
        let end = source.start - 1001 + source.length;
        let mut sequence = reference.as_bytes()[source.start - 1001..end].to_vec();
        let reverse = source.extraction == Some('-');
        if reverse {
            sequence = GentleEngine::reverse_complement_bytes(&sequence);
        }
        let mut dna = DNAsequence::from_sequence(&String::from_utf8(sequence).unwrap()).unwrap();
        let mut exons = source
            .exons
            .into_iter()
            .map(|(s, e)| {
                if reverse {
                    (source.length - e, source.length - s)
                } else {
                    (s, e)
                }
            })
            .collect::<Vec<_>>();
        exons.sort_unstable();
        let mut location = gb_io::seq::Location::Join(
            exons
                .into_iter()
                .map(|(s, e)| gb_io::seq::Location::simple_range(s as i64, e as i64))
                .collect(),
        );
        if (source.genomic_strand == '-') != reverse {
            location = gb_io::seq::Location::Complement(Box::new(location));
        }
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location,
            qualifiers: vec![
                ("gene".into(), Some(format!("SYN_{i}"))),
                ("transcript_id".into(), Some(format!("annotation_{i}"))),
            ],
        });
        let mut qualifiers = vec![
            (
                "genomic_start_1based".into(),
                Some(source.start.to_string()),
            ),
            (
                "genomic_end_1based".into(),
                Some((source.start + source.length - 1).to_string()),
            ),
            ("chromosome".into(), Some("synthetic_chr".into())),
            ("genome_id".into(), Some("synthetic".into())),
        ];
        if let Some(strand) = source.extraction {
            qualifiers.push(("strand".into(), Some(strand.to_string())));
        }
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "source".into(),
            location: gb_io::seq::Location::simple_range(0, source.length as i64),
            qualifiers,
        });
        engine
            .state_mut()
            .sequences
            .insert(format!("source_{i}"), dna);
    }
    (engine, build)
}

fn refresh(panel: &mut RnaReadPanel) {
    let (pairs, assessment) = overlap::assess(panel);
    panel.overlapping_source_pairs = pairs;
    panel.source_overlap_assessment = Some(assessment);
    seal_panel(panel).unwrap();
}

#[test]
fn rna_panel_opposite_record_orientations_cannot_hide_same_strand_duplicate_exons() {
    let td = tempfile::tempdir().unwrap();
    let mut second = plus_source();
    second.extraction = Some('-');
    second.exons[0].0 += 1; // Different annotation, not identical selected mature RNA.
    let (engine, build) = sources(td.path(), [plus_source(), second]);
    let panel = engine.build_rna_read_panel(&build).unwrap();
    assert_ne!(
        panel.loci[0].source_sequence_sha256,
        panel.loci[1].source_sequence_sha256
    );
    assert_ne!(
        panel.templates[0].sequence_sha256,
        panel.templates[1].sequence_sha256
    );
    assert!(
        panel
            .shared_transcript_sequences
            .as_ref()
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        panel.overlapping_source_pairs,
        vec![("locus_0".into(), "locus_1".into())]
    );
    let run = request(&panel, td.path(), &panel.templates[0].sequence);
    let unreadable_fastq = b"not FASTQ: overlapping sources must refuse before streaming\n";
    // Path admission requires an existing file; invalid contents prove streaming never starts.
    fs::write(&run.files[0].path, unreadable_fastq).unwrap();
    let mut progress_calls = 0;
    let error = engine
        .screen_rna_read_panel(run.clone(), &mut |_| {
            progress_calls += 1;
            true
        })
        .unwrap_err();
    assert!(
        error.message.starts_with("panel_overlapping_sources:"),
        "{error:?}"
    );
    assert_eq!(progress_calls, 0);
    assert_eq!(
        fs::read(&run.files[0].path).unwrap(),
        unreadable_fastq.to_vec()
    );
    assert!(!td.path().join("run").exists());
}

#[test]
fn rna_panel_antisense_assessment_and_nominations_are_extraction_orientation_invariant() {
    let td = tempfile::tempdir().unwrap();
    let mut observations = Vec::new();
    for (i, orientation) in ['+', '-'].into_iter().enumerate() {
        let dir = td.path().join(format!("case_{i}"));
        fs::create_dir(&dir).unwrap();
        let (engine, build) = sources(&dir, [plus_source(), antisense_source(orientation)]);
        let panel = engine.build_rna_read_panel(&build).unwrap();
        assert!(panel.overlapping_source_pairs.is_empty());
        let assessment = panel.source_overlap_assessment.clone().unwrap();
        assert!(assessment.conservative_source_pairs.is_empty());
        assert_eq!(
            assessment.antisense_overlapping_loci,
            vec![RnaReadPanelAntisenseOverlap {
                locus_ids: ("locus_0".into(), "locus_1".into()),
                assembly: "synthetic".into(),
                chromosome: "synthetic_chr".into(),
                genomic_strands: ('+', '-'),
                member_ids_a: vec!["locus_0:0".into()],
                member_ids_b: vec!["locus_1:0".into()],
                overlapping_intervals_1based: vec![(1031, 1060), (1091, 1120)],
            }]
        );
        let run = request(&panel, &dir, &panel.templates[0].sequence);
        let report = engine.screen_rna_read_panel(run, &mut |_| true).unwrap();
        let file = &report.files_completed[0];
        assert_eq!(file.classes.get("multi_locus_ambiguous"), Some(&1));
        assert_eq!(file.equivalence_classes.len(), 1);
        observations.push((
            assessment,
            file.classes.clone(),
            file.equivalence_classes.clone(),
        ));
    }
    assert_eq!(observations[0], observations[1]);
}

#[test]
fn rna_panel_adjacent_exons_are_not_duplicates_despite_overlapping_record_flanks() {
    let td = tempfile::tempdir().unwrap();
    let mut first = plus_source();
    first.exons = vec![(0, 20), (30, 50)];
    let mut second = antisense_source('+');
    second.genomic_strand = '+';
    second.exons = vec![(90, 110), (140, 160)];
    let (engine, build) = sources(td.path(), [first, second]);
    let panel = engine.build_rna_read_panel(&build).unwrap();
    assert_eq!(
        source_overlaps(&panel.loci).len(),
        1,
        "legacy whole-record rule would refuse"
    );
    assert!(panel.overlapping_source_pairs.is_empty());
    let assessment = panel.source_overlap_assessment.as_ref().unwrap();
    assert!(assessment.antisense_overlapping_loci.is_empty());
    assert!(assessment.conservative_source_pairs.is_empty());
    let run = request(&panel, td.path(), &panel.templates[0].sequence);
    engine.screen_rna_read_panel(run, &mut |_| true).unwrap();
}

#[test]
fn rna_panel_genomic_exon_endpoints_distinguish_adjacency_from_one_base_overlap() {
    let td = tempfile::tempdir().unwrap();
    for (i, (start, orientation)) in [(1031, '+'), (1031, '-'), (1030, '+'), (1030, '-')]
        .into_iter()
        .enumerate()
    {
        let dir = td.path().join(format!("boundary_{i}"));
        fs::create_dir(&dir).unwrap();
        let mut first = plus_source();
        first.exons = vec![(0, 30)];
        let mut second = antisense_source(orientation);
        second.genomic_strand = '+';
        second.start = start;
        second.exons = vec![(0, 30)];
        let (engine, build) = sources(&dir, [first, second]);
        let panel = engine.build_rna_read_panel(&build).unwrap();
        assert_eq!(
            panel.overlapping_source_pairs.len(),
            usize::from(start == 1030)
        );
    }
}

#[test]
fn rna_panel_identity_checks_still_block_antisense_despite_distinct_strands() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = sources(td.path(), [plus_source(), antisense_source('+')]);
    let original = engine.build_rna_read_panel(&build).unwrap();
    for mutation in 0..3 {
        let mut panel = original.clone();
        match mutation {
            0 => panel.loci[1].source.seq_id = panel.loci[0].source.seq_id.clone(),
            1 => panel.loci[1].source.gene_id = panel.loci[0].source.gene_id.clone(),
            _ => {
                panel.loci[1].source_sequence_sha256 = panel.loci[0].source_sequence_sha256.clone()
            }
        }
        refresh(&mut panel);
        if mutation == 1 {
            panel.gene_count = 1;
            seal_panel(&mut panel).unwrap();
        }
        validate_panel(&panel).unwrap();
        assert_eq!(
            panel.overlapping_source_pairs,
            vec![("locus_0".into(), "locus_1".into())]
        );
    }
}

#[test]
fn rna_panel_unknown_or_unmappable_anchor_uses_disclosed_conservative_overlap() {
    let td = tempfile::tempdir().unwrap();
    let mut second = antisense_source('+');
    second.extraction = None;
    let (engine, build) = sources(td.path(), [plus_source(), second]);
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let assessment = panel.source_overlap_assessment.as_ref().unwrap();
    assert_eq!(
        assessment.conservative_source_pairs[0].reasons,
        vec![RnaReadPanelOverlapFallbackReason::UnknownAnchorStrand]
    );
    assert!(assessment.antisense_overlapping_loci.is_empty());
    assert_eq!(panel.overlapping_source_pairs.len(), 1);
    let run = request(&panel, td.path(), &panel.templates[0].sequence);
    assert!(
        engine
            .screen_rna_read_panel(run, &mut |_| true)
            .unwrap_err()
            .message
            .contains("panel_overlapping_sources")
    );
    let mut panel = panel;
    let anchor = panel.loci[1].genome_anchor.as_mut().unwrap();
    anchor.strand = Some('-');
    anchor.end_1based = anchor.start_1based + 2;
    refresh(&mut panel);
    validate_panel(&panel).unwrap();
    assert_eq!(
        panel
            .source_overlap_assessment
            .as_ref()
            .unwrap()
            .conservative_source_pairs[0]
            .reasons,
        vec![RnaReadPanelOverlapFallbackReason::UnmappableExons]
    );
    panel.loci[1].genome_anchor.as_mut().unwrap().start_1based = 0;
    refresh(&mut panel);
    validate_panel(&panel).unwrap();
    assert_eq!(
        panel
            .source_overlap_assessment
            .as_ref()
            .unwrap()
            .conservative_source_pairs[0]
            .reasons,
        vec![RnaReadPanelOverlapFallbackReason::InvalidAnchorGeometry]
    );
}

#[test]
fn rna_panel_exon_overlap_requires_exact_assembly_and_chromosome_context() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = sources(td.path(), [plus_source(), antisense_source('+')]);
    let original = engine.build_rna_read_panel(&build).unwrap();
    for change_chromosome in [false, true] {
        let mut panel = original.clone();
        if change_chromosome {
            panel.loci[1].genome_anchor.as_mut().unwrap().chromosome = "other_chr".into();
        } else {
            panel.loci[1].source.assembly = "other_assembly".into();
        }
        refresh(&mut panel);
        validate_panel(&panel).unwrap();
        assert!(panel.overlapping_source_pairs.is_empty());
        assert!(
            panel
                .source_overlap_assessment
                .as_ref()
                .unwrap()
                .antisense_overlapping_loci
                .is_empty()
        );
    }
}

#[test]
fn rna_panel_rehashed_antisense_geometry_and_refusal_tampering_is_rejected() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = sources(td.path(), [plus_source(), antisense_source('+')]);
    let original = engine.build_rna_read_panel(&build).unwrap();
    for mutation in 0..5 {
        let mut panel = original.clone();
        match mutation {
            0 => {
                panel
                    .source_overlap_assessment
                    .as_mut()
                    .unwrap()
                    .antisense_overlapping_loci[0]
                    .overlapping_intervals_1based[0]
                    .0 += 1
            }
            1 => panel
                .source_overlap_assessment
                .as_mut()
                .unwrap()
                .antisense_overlapping_loci[0]
                .member_ids_a
                .clear(),
            2 => {
                panel
                    .source_overlap_assessment
                    .as_mut()
                    .unwrap()
                    .antisense_overlapping_loci[0]
                    .genomic_strands = ('+', '+')
            }
            3 => panel
                .overlapping_source_pairs
                .push(("locus_0".into(), "locus_1".into())),
            _ => panel.source_overlap_assessment = None, // No silent downgrade to legacy admission.
        }
        seal_panel(&mut panel).unwrap();
        assert_eq!(
            validate_panel(&panel).unwrap_err().message,
            "panel_source_overlap_binding"
        );
    }
    let mut json = serde_json::to_value(original).unwrap();
    json["source_overlap_assessment"]["rule"] = json!("guessed_strand");
    assert!(serde_json::from_value::<RnaReadPanel>(json).is_err());
}

#[test]
fn rna_panel_historical_overlap_rule_and_wire_digest_are_not_reinterpreted() {
    // Hand-authored historical wire shape, not a new builder result stripped of its rule.
    let mut value = json!({"schema":"gentle.rna_read_panel.v1", "content_sha256":"",
        "loci":[], "templates":[], "locus_count":2, "gene_count":2,
        "selected_transcript_count":2, "annotated_transcript_count":2, "overlapping_source_pairs":[]});
    for i in 0..2 {
        let sequence = if i == 0 {
            "ACGT".repeat(10)
        } else {
            "GGCA".repeat(10)
        };
        value["loci"].as_array_mut().unwrap().push(json!({
            "source":{"locus_id":format!("locus_{i}"), "seq_id":format!("old_{i}"), "source_feature_id":0,
                "source":"synthetic", "assembly":"synthetic", "annotation_release":"old",
                "gene_id":format!("GENE_{i}"), "gene_symbol":format!("GENE_{i}"), "selected_transcript_ids":[]},
            "source_record_sha256":format!("sha256:{}", if i == 0 { "a" } else { "b" }.repeat(64)),
            "genome_anchor":{"seq_id":format!("old_{i}"), "genome_id":"synthetic", "chromosome":"synthetic_chr",
                "start_1based":if i == 0 {1001} else {961}, "end_1based":if i == 0 {1100} else {1060},
                "strand":if i == 0 {"+"} else {"-"}},
            "annotated_transcript_count":1, "annotated_transcript_ids":[format!("tx_{i}")]}));
        value["templates"].as_array_mut().unwrap().push(json!({
            "member_id":format!("locus_{i}:0"), "locus_id":format!("locus_{i}"), "transcript_feature_id":0,
            "transcript_id":format!("tx_{i}"), "strand":if i == 0 {"+"} else {"-"},
            "sequence_sha256":sha256_prefixed_str(&sequence), "sequence":sequence,
            "junction_offsets":[], "source_exons":[[10,50]]}));
    }
    value["content_sha256"] = json!(sha256_prefixed_bytes(&encode(&value).unwrap()));
    let bytes = encode(&value).unwrap();
    let panel: RnaReadPanel = serde_json::from_slice(&bytes).unwrap();
    validate_panel(&panel).unwrap();
    assert!(panel.source_overlap_assessment.is_none());
    assert!(panel.overlapping_source_pairs.is_empty());
    assert_eq!(
        overlap::assess(&panel).0,
        vec![("locus_0".into(), "locus_1".into())]
    );
    assert_eq!(encode(&panel).unwrap(), bytes);
}

#[test]
fn rna_panel_antisense_evidence_survives_shared_inspection_summary_and_legacy_absence() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = sources(td.path(), [plus_source(), antisense_source('+')]);
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let shown = engine
        .inspect_rna_read_panel(&build.output_path, None, 0)
        .unwrap();
    assert_eq!(
        shown.0.unwrap().source_overlap_assessment,
        panel.source_overlap_assessment
    );
    let run = request(&panel, td.path(), &panel.templates[0].sequence);
    engine.screen_rna_read_panel(run, &mut |_| true).unwrap();
    let export = ExportRnaReadPanelSummaryRequest {
        schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
        report_paths: vec![td.path().join("run/report.json").display().to_string()],
        output_dir: td.path().join("summary").display().to_string(),
    };
    let summary = engine.export_rna_read_panel_summary(&export).unwrap();
    assert_eq!(
        summary.source_overlap_assessment,
        panel.source_overlap_assessment
    );
    let path = td.path().join("summary/report.json");
    let shown = engine
        .inspect_rna_read_panel(&path.to_string_lossy(), None, 0)
        .unwrap();
    assert_eq!(
        shown.3.unwrap().source_overlap_assessment,
        panel.source_overlap_assessment
    );
    engine
        .verify_rna_read_panel_summary_artifacts(&path.to_string_lossy(), &summary)
        .unwrap();
    let mut invalid = summary.clone();
    invalid
        .source_overlap_assessment
        .as_mut()
        .unwrap()
        .antisense_overlapping_loci[0]
        .genomic_strands = ('+', '+');
    invalid.content_sha256.clear();
    invalid.content_sha256 = sha256_prefixed_bytes(&encode(&invalid).unwrap());
    fs::write(&path, serde_json::to_vec_pretty(&invalid).unwrap()).unwrap();
    assert!(
        engine
            .inspect_rna_read_panel(&path.to_string_lossy(), None, 0)
            .unwrap_err()
            .message
            .contains("panel_summary_source_overlap_binding")
    );
    let mut legacy = summary;
    legacy.source_overlap_assessment = None;
    legacy.content_sha256.clear();
    legacy.content_sha256 = sha256_prefixed_bytes(&encode(&legacy).unwrap());
    let bytes = encode(&legacy).unwrap();
    let decoded: RnaReadPanelSummary = serde_json::from_slice(&bytes).unwrap();
    summary::validate_sequence_sharing(&decoded).unwrap();
    assert_eq!(encode(&decoded).unwrap(), bytes);
}

#[test]
fn rna_panel_antisense_summary_retains_source_geometry_without_verified_files() {
    let td = tempfile::tempdir().unwrap();
    let (engine, build) = sources(td.path(), [plus_source(), antisense_source('+')]);
    let panel = engine.build_rna_read_panel(&build).unwrap();
    let run = request(&panel, td.path(), &panel.templates[0].sequence);
    let report = engine.screen_rna_read_panel(run, &mut |_| false).unwrap();
    assert!(report.files_completed.is_empty());
    let export = ExportRnaReadPanelSummaryRequest {
        schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
        report_paths: vec![td.path().join("run/report.json").display().to_string()],
        output_dir: td.path().join("summary").display().to_string(),
    };
    let summary = engine.export_rna_read_panel_summary(&export).unwrap();
    assert!(summary.transcripts.is_empty());
    assert_eq!(
        summary.source_overlap_assessment,
        panel.source_overlap_assessment
    );
    let path = td.path().join("summary/report.json");
    engine
        .inspect_rna_read_panel(&path.to_string_lossy(), None, 0)
        .unwrap();
    engine
        .verify_rna_read_panel_summary_artifacts(&path.to_string_lossy(), &summary)
        .unwrap();
}
