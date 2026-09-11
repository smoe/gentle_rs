//! Synthetic context joins, including both genomic and loaded-sequence strands.
//!
//! All sequences, names, coordinates and scores are hand-crafted software tests,
//! not experimental data. `fixture` deterministically recreates temporary inputs;
//! no private loci, vendor tables or network resources are used.

use super::*;
use crate::engine::Engine;
use serde_json::json;

fn fixture(
    strand: TssStrand,
    loaded_reverse: bool,
) -> (
    TssRecord,
    TssReference,
    GeneLocusEvidenceDisplayReport,
    String,
    TataBoxScreenReport,
) {
    let genomic = "ACGTTGAACCTA".repeat(10);
    let sequence = if loaded_reverse {
        GentleEngine::reverse_complement(&genomic)
    } else {
        genomic.clone()
    };
    let window_sequence = if strand == TssStrand::Minus {
        GentleEngine::reverse_complement(&genomic[10..110])
    } else {
        genomic[10..110].to_string()
    };
    let local = |start: usize, end: usize| {
        if loaded_reverse {
            (220 - end, 220 - start)
        } else {
            (start - 99, end - 99)
        }
    };
    let local_reverse = loaded_reverse ^ (strand == TssStrand::Minus);
    let exon = |start, end| {
        let (start, end) = local(start, end);
        json!({"start_1based":start, "end_1based":end})
    };
    let seq_hash = sha256_hex_bytes(sequence.as_bytes());
    let record = TssRecord {
        promoter_id: "synthetic-detail".into(),
        gene_id: "synthetic-gene".into(),
        gene_symbol: "SYNTH".into(),
        geometry: TssGeometry {
            chromosome: "chrSynthetic".into(),
            strand,
            tss_1based: if strand == TssStrand::Plus { 125 } else { 194 },
            start_1based: 110,
            end_1based: 209,
            upstream_bp: 15,
            downstream_bp: 84,
        },
        transcripts: vec!["synthetic.tx".into()],
        sequence_sha256: sha256_hex_bytes(window_sequence.as_bytes()),
    };
    let reference = TssReference {
        genome_id: "synthetic-genome".into(),
        assembly: "synthetic-assembly".into(),
        annotation_release: Some("synthetic-release".into()),
    };
    let (occ_start, occ_end) = local(105, 125);
    let (start, stop) = if strand == TssStrand::Plus {
        (130, 188)
    } else {
        (190, 132)
    };
    let locus: GeneLocusEvidenceDisplayReport = serde_json::from_value(json!({
        "schema": gentle_protocol::GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA, "seq_id":"synthetic-locus", "gene_symbol":"SYNTH", "gene_strand":strand.as_str(),
        "locus_local_start_1based":1,"locus_local_end_1based":120,
        "sequence_binding":{"sequence_sha256":format!("sha256:{seq_hash}"), "sequence_length_bp":120,
            "genome_anchor":{"genome_id":reference.genome_id,"chromosome":"Synthetic","start_1based":100,"end_1based":219,"strand":if loaded_reverse {"-"} else {"+"}}},
        "isoform_evidence":{"assembly":reference.assembly,"chromosome":"Synthetic","gene_symbol":"SYNTH","annotation_release":reference.annotation_release,
            "splicing":{"seq_id":"synthetic-locus","target_feature_id":0,"group_label":"SYNTH","strand":if local_reverse {"-"} else {"+"},
                "region_start_1based":1,"region_end_1based":120,"transcript_count":1,"unique_exon_count":2,"instruction":"synthetic",
                "transcripts":[{"transcript_feature_id":0,"transcript_id":"synthetic.tx","label":"Synthetic transcript", "strand":if local_reverse {"-"} else {"+"},"exons":[exon(125,147),exon(183,194)],"introns":[],"has_target_feature":true}],
                "unique_exons":[],"matrix_rows":[],"boundaries":[],"junctions":[],"events":[]}},
        "transcript_metrics":[{"transcript_id":"synthetic.tx","coding_status":"complete_cds","cds_ranges_local_1based":[local(130,147),local(183,190)]}],
        "codon_markers":[{"transcript_id":"synthetic.tx","kind":"start","local_position_1based":local(start,start).0,"genomic_position_1based":start,"strand":if local_reverse {"-"} else {"+"},"genomic_strand":strand.as_str(),"basis":"annotated CDS boundary"},
            {"transcript_id":"synthetic.tx","kind":"stop","local_position_1based":local(stop,stop).0,"genomic_position_1based":stop,"strand":if local_reverse {"-"} else {"+"},"genomic_strand":strand.as_str(),"basis":"annotated CDS boundary"}],
        "occupancy_groups":[{"group_id":"sample","label":"Synthetic CUT&RUN","scale_mode":"shared_group","group_abs_max_score":10.0,"lanes":[
            {"state":"available","source_id":"synthetic-signal","assay":"CUT&RUN","factor":"synthetic-factor","role":"experimental","display_abs_max_score":10.0,
             "lane":{"lane_id":"signal","track_name":"signal","display_label":"Synthetic signal","source_kind":"bedgraph","intervals":[{"interval_id":"peak","local_start_1based":occ_start,"local_end_1based":occ_end,"genomic_start_1based":105,"genomic_end_1based":125,"score":4.0}]}},
            {"state":"not_prepared","source_id":"missing-control","role":"input_control","display_abs_max_score":10.0,"lane":{"lane_id":"control","display_label":"Unavailable control"}}]}]
    })).unwrap();
    let anchor = locus
        .sequence_binding
        .as_ref()
        .unwrap()
        .genome_anchor
        .as_ref()
        .unwrap();
    let (box_start, box_end) = local(115, 120);
    let mut tata: TataBoxScreenReport = serde_json::from_value(json!({
        "schema":TATA_BOX_SCREEN_SCHEMA,"report_id":"synthetic-tata","request":{"seq_id":"synthetic-locus"},"content_sha256":"",
        "sequence_sha256":format!("sha256:{seq_hash}"),"annotation_sha256":format!("sha256:{}","1".repeat(64)),"anchor_sha256":format!("sha256:{}","2".repeat(64)),
        "genome_id":anchor.genome_id,"chromosome":anchor.chromosome,"genomic_start_1based":100,"genomic_end_1based":219,"genomic_reverse":loaded_reverse,
        "motif_id":null,"matrix_sha256":null,"score_policy":"synthetic scores, not a biological prediction","epd_status":"not_requested","epd_bed_sha256":null,"epd_motifs_sha256":null,"tss":[],
        "rows":[{"row_id":"synthetic-box","evidence_kind":"source_annotation","label":"Synthetic TATA annotation","start_0based":box_start-1,"end_0based_exclusive":box_end,"reverse":local_reverse,"geometry_kind":"source_feature_interval","source_feature_id":null,"source_qualifiers":[],"tata_positive":null,"llr_bits":null,"sequence_5prime_to_3prime":null,"tss_associations":[]}],
        "scored_windows":0,"ambiguous_windows":0,"warnings":[],"non_claim":TATA_BOX_NON_CLAIM
    })).unwrap();
    tata.content_sha256 = format!(
        "sha256:{}",
        sha256_hex_bytes(&serde_json::to_vec(&tata).unwrap())
    );
    (record, reference, locus, sequence, tata)
}

fn project_fixture(strand: TssStrand, loaded_reverse: bool) -> (TssRecord, TssDetailContext) {
    let (record, reference, locus, sequence, tata) = fixture(strand, loaded_reverse);
    let bindings = vec![
        TssInputBinding {
            role: "tss_detail_locus_report".into(),
            name: "locus.json".into(),
            sha256: "a".repeat(64),
        },
        TssInputBinding {
            role: "tss_detail_tata_report".into(),
            name: "tata.json".into(),
            sha256: "b".repeat(64),
        },
    ];
    let context = project(
        &record,
        &reference,
        &locus,
        &sequence,
        &"a".repeat(64),
        bindings,
        Some((&tata, &"b".repeat(64))),
    )
    .unwrap();
    (record, context)
}

#[test]
fn tss_context_anchor_endpoint_arithmetic_stays_within_validated_bounds() {
    let (_, _, locus, _, _) = fixture(TssStrand::Plus, false);
    let mut anchor = locus.sequence_binding.unwrap().genome_anchor.unwrap();
    anchor.start_1based = usize::MAX - 9;
    anchor.end_1based = usize::MAX;
    for strand in ['+', '-'] {
        anchor.strand = Some(strand);
        assert_eq!(
            local_genomic(&anchor, 1, 10).unwrap(),
            ((usize::MAX - 9) as u64, usize::MAX as u64)
        );
        assert!(local_genomic(&anchor, 0, 10).is_err());
        assert!(local_genomic(&anchor, 1, 11).is_err());
    }
}

#[test]
fn tss_context_projects_transcript_cds_codons_signal_and_tata_in_all_orientations() {
    for strand in [TssStrand::Plus, TssStrand::Minus] {
        for loaded_reverse in [false, true] {
            let (_, context) = project_fixture(strand, loaded_reverse);
            let row = &context.transcripts[0];
            assert_eq!(row.exons.len(), 2);
            assert_eq!(row.exons[0].number_5prime_to_3prime, 1);
            assert_eq!(
                row.exons[0].span.genomic_start_1based,
                if strand == TssStrand::Plus { 125 } else { 183 }
            );
            assert_eq!(row.cds.len(), 2);
            assert_eq!(row.codons.len(), 2);
            assert_eq!(
                row.codons[0].position_0based,
                if strand == TssStrand::Plus { 20 } else { 19 }
            );
            assert_eq!(row.codons[1].kind, GeneLocusCodonKind::Stop);
            let peak = &context.occupancy[0].intervals[0];
            assert_eq!(
                peak.span.start_0based,
                if strand == TssStrand::Plus { 0 } else { 84 }
            );
            assert_eq!(
                peak.span.end_0based_exclusive,
                if strand == TssStrand::Plus { 16 } else { 100 }
            );
            assert!(peak.span.clipped);
            assert_eq!(peak.score, Some(4.0));
            assert_eq!(context.occupancy[0].display_abs_max_score, 10.0);
            assert_eq!(
                context.occupancy[1].state,
                GeneLocusOccupancyLaneState::NotPrepared
            );
            assert!(context.occupancy[1].intervals.is_empty());
            let tata = &context.tata.as_ref().unwrap().rows[0];
            assert_eq!(
                tata.span.start_0based,
                if strand == TssStrand::Plus { 5 } else { 89 }
            );
            assert_eq!(tata.genomic_strand, strand);
        }
    }
}

#[test]
fn tss_context_rejects_reference_sequence_geometry_and_tata_corruption() {
    for case in 0..9 {
        let (mut record, mut reference, mut locus, sequence, mut tata) =
            fixture(TssStrand::Plus, false);
        match case {
            0 => reference.assembly = "wrong".into(),
            1 => reference.annotation_release = Some("stale".into()),
            2 => record.sequence_sha256 = "f".repeat(64),
            3 => locus.gene_symbol = "wrong".into(),
            4 => {
                locus
                    .sequence_binding
                    .as_mut()
                    .unwrap()
                    .genome_anchor
                    .as_mut()
                    .unwrap()
                    .strand = None
            }
            5 => locus.occupancy_groups[0].lanes[0].lane.intervals[0].genomic_start_1based += 1,
            6 => locus.codon_markers[0].genomic_position_1based += 1,
            7 => tata.rows[0].start_0based += 1,
            _ => tata.sequence_sha256 = "f".repeat(64),
        }
        assert!(
            project(
                &record,
                &reference,
                &locus,
                &sequence,
                &"a".repeat(64),
                vec![],
                Some((&tata, &"b".repeat(64)))
            )
            .is_err(),
            "case {case}"
        );
    }
}

#[test]
fn tss_context_missing_evidence_is_unavailable_not_negative_and_does_not_invent_transcripts() {
    let (mut record, reference, mut locus, sequence, _) = fixture(TssStrand::Plus, false);
    record.transcripts.push("missing.tx".into());
    locus.occupancy_groups.clear();
    let c = project(
        &record,
        &reference,
        &locus,
        &sequence,
        &"a".repeat(64),
        vec![],
        None,
    )
    .unwrap();
    assert_eq!(c.transcripts.len(), 1);
    assert!(c.occupancy.is_empty());
    assert!(c.tata.is_none());
    assert!(c.warnings.iter().any(|s| s.contains("missing.tx")));
}

#[test]
fn tss_context_never_promotes_an_inferred_or_unclassified_orf_to_annotated_cds() {
    for status in ["inferred_orf", "", "partial_cds"] {
        let (record, reference, mut locus, sequence, _) = fixture(TssStrand::Plus, false);
        locus.transcript_metrics[0].coding_status = status.into();
        locus.codon_markers.clear();
        let c = project(
            &record,
            &reference,
            &locus,
            &sequence,
            &"a".repeat(64),
            vec![],
            None,
        )
        .unwrap();
        assert_eq!(c.transcripts[0].exons.len(), 2);
        assert!(c.transcripts[0].codons.is_empty());
        if status == "partial_cds" {
            assert_eq!(c.transcripts[0].cds.len(), 2);
        } else {
            assert!(c.transcripts[0].cds.is_empty());
            assert!(
                c.warnings
                    .iter()
                    .any(|w| w.contains("not displayed as annotated CDS"))
            );
        }
    }
}

#[test]
fn tss_context_bundled_example_uses_resolver_and_renders_all_tata_evidence_kinds() {
    let directory = Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/fixtures/tss_profiles");
    let manifest: TssBundleManifest =
        serde_json::from_slice(&std::fs::read(directory.join("manifest.json")).unwrap()).unwrap();
    let mut report = crate::tss_profile_export::tests::synthetic_report();
    report.reference = manifest.reference;
    report.windows.truncate(2);
    for (window, record) in report.windows.iter_mut().zip(manifest.records) {
        window.record = record;
    }
    let original = serde_json::to_vec(&report.windows[0].tracks).unwrap();
    attach(
        &mut report,
        &directory.join("context/manifest.json"),
        &mut || true,
    )
    .unwrap();
    assert!(report.windows[0].detail_context.is_none());
    assert_eq!(
        serde_json::to_vec(&report.windows[0].tracks).unwrap(),
        original
    );
    let context = report.windows[1].detail_context.as_ref().unwrap();
    assert_eq!(context.tata.as_ref().unwrap().rows.len(), 3);
    for window in &mut report.windows {
        for track in &mut window.tracks {
            let n = 11 - track.motif_length_bp + 1;
            track.forward_scores = vec![Some(0.0); n];
            track.reverse_scores = vec![Some(0.0); n];
            track.forward_peaks.clear();
            track.reverse_peaks.clear();
            track.forward_maximum = None;
            track.reverse_maximum = None;
        }
        window.comparisons = GentleEngine::tss_matrix_comparisons(window, &report.panel_resolution);
    }
    let pages = gentle_render::tss_profiles::render_tss_profile_pages(&report, &Default::default())
        .unwrap();
    let svg = pages
        .iter()
        .find(|p| p.gene_symbol == "SYNMINUS")
        .unwrap()
        .svg
        .as_str();
    for role in [
        "tata-annotation",
        "tata-prediction",
        "tata-epd-tss",
        "context-exon",
        "context-translation-start",
    ] {
        assert!(svg.contains(role));
    }
    assert!(svg.contains("TATA-positive: yes"));
    let mut bad = report.clone();
    bad.windows[1].detail_context.as_mut().unwrap().occupancy[0].intervals[0]
        .span
        .start_0based += 1;
    assert!(
        gentle_render::tss_profiles::render_tss_profile_pages(&bad, &Default::default()).is_err()
    );
    let mut canceled = crate::tss_profile_export::tests::synthetic_report();
    canceled.reference = report.reference.clone();
    canceled.windows = report.windows.clone();
    for w in &mut canceled.windows {
        w.detail_context = None;
    }
    assert!(
        attach(
            &mut canceled,
            &directory.join("context/manifest.json"),
            &mut || false
        )
        .is_err()
    );
    assert!(canceled.windows.iter().all(|w| w.detail_context.is_none()));
}

#[test]
fn tss_context_export_is_read_only_hash_bound_and_replays_without_source_files() {
    let (record, reference, locus, sequence, tata) = fixture(TssStrand::Minus, true);
    let mut report = crate::tss_profile_export::tests::synthetic_report();
    report.reference = reference.clone();
    report.windows.truncate(1);
    let window = &mut report.windows[0];
    window.record = record;
    window.selected = false;
    window.selection_evidence = None;
    window.comparisons.clear();
    for track in &mut window.tracks {
        let n = 100 - track.motif_length_bp + 1;
        track.forward_scores = vec![Some(0.0); n];
        track.reverse_scores = vec![Some(0.0); n];
        track.forward_maximum = None;
        track.reverse_maximum = None;
        track.forward_peaks.clear();
        track.reverse_peaks.clear();
    }
    window.comparisons = GentleEngine::tss_matrix_comparisons(window, &report.panel_resolution);
    let original_tracks = serde_json::to_value(&report.windows[0].tracks).unwrap();
    let retain = std::env::var_os("GENTLE_TEST_TSS_PROFILE_OUTPUT_DIR");
    let temp = match &retain {
        Some(path) => tempfile::Builder::new()
            .prefix("gentle-tss-context-")
            .tempdir_in(path),
        None => tempfile::tempdir(),
    }
    .unwrap();
    let root = std::fs::canonicalize(temp.path()).unwrap();
    let inputs = root.join("inputs");
    std::fs::create_dir(&inputs).unwrap();
    let write = |name: &str, bytes: Vec<u8>| {
        std::fs::write(inputs.join(name), &bytes).unwrap();
        TssContextFile {
            path: name.into(),
            sha256: sha256_hex_bytes(&bytes),
        }
    };
    let locus_file = write("locus.json", serde_json::to_vec(&locus).unwrap());
    let fasta_file = write("locus.fa", format!(">synthetic\n{sequence}\n").into_bytes());
    let tata_file = write("tata.json", serde_json::to_vec(&tata).unwrap());
    let manifest = TssContextManifest {
        schema: CONTEXT_INPUT_SCHEMA.into(),
        reference,
        genes: vec![TssContextSource {
            gene_id: "synthetic-gene".into(),
            locus_report: locus_file,
            locus_fasta: fasta_file,
            tata_report: Some(tata_file),
        }],
    };
    let manifest_path = inputs.join("context.json");
    std::fs::write(&manifest_path, serde_json::to_vec(&manifest).unwrap()).unwrap();
    let mut engine = GentleEngine::new();
    let before = serde_json::to_value(engine.state()).unwrap();
    let output = root.join("output");
    let request = ExportTssProfilesRequest {
        output_dir: output.to_str().unwrap().into(),
        context_manifest: Some(manifest_path.to_str().unwrap().into()),
        rendering: Default::default(),
        formats: vec![
            TssExportFormat::Svg,
            TssExportFormat::Png,
            TssExportFormat::Pdf,
        ],
    };
    let result = engine
        .apply(crate::engine::Operation::ExportTssTfbsProfiles {
            report: Box::new(report.clone()),
            request: request.clone(),
        })
        .unwrap();
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    let enriched: TssProfileReport =
        serde_json::from_slice(&std::fs::read(output.join("report.json")).unwrap()).unwrap();
    assert_eq!(
        serde_json::to_value(&enriched.windows[0].tracks).unwrap(),
        original_tracks
    );
    let receipt = result.tss_tfbs_profile_receipt.unwrap();
    assert!(
        receipt
            .inputs
            .iter()
            .any(|b| b.role == "tss_detail_tata_report")
    );
    let svg_name = receipt
        .outputs
        .keys()
        .find(|name| name.ends_with(".svg"))
        .unwrap();
    let svg = std::fs::read_to_string(output.join(svg_name)).unwrap();
    for marker in [
        "context-exon",
        "context-cds",
        "context-translation-start",
        "context-translation-stop",
        "tss-context-occupancy",
        "tata-annotation",
        "data-state=\"not_prepared\"",
    ] {
        assert!(svg.contains(marker), "{marker}");
    }
    assert!(!svg.contains("NaN"));
    let mut changed = request.clone();
    changed.output_dir = root.join("bad").to_str().unwrap().into();
    std::fs::write(inputs.join("locus.fa"), b">wrong\nAAAA\n").unwrap();
    assert!(
        engine
            .apply(crate::engine::Operation::ExportTssTfbsProfiles {
                report: Box::new(report),
                request: changed.clone()
            })
            .is_err()
    );
    assert!(!Path::new(&changed.output_dir).exists());
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    std::fs::remove_dir_all(inputs).unwrap();
    let replay_request = ExportTssProfilesRequest {
        output_dir: root.join("replay").to_str().unwrap().into(),
        context_manifest: None,
        ..request
    };
    let replay =
        crate::tss_profile_export::export_tss_profiles(&enriched, &replay_request).unwrap();
    assert_eq!(receipt.report_sha256, replay.report_sha256);
    assert_eq!(
        std::fs::read(output.join(svg_name)).unwrap(),
        std::fs::read(Path::new(&replay_request.output_dir).join(svg_name)).unwrap()
    );
    if retain.is_some() {
        eprintln!(
            "Retained synthetic TSS context artifacts: {}",
            root.display()
        );
        let _ = temp.keep();
    }
}
