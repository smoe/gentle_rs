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
fn tss_context_accepts_exact_assembly_anchor_but_rejects_catalog_tokens_and_near_matches() {
    let (record, reference, mut locus, sequence, _) = fixture(TssStrand::Plus, false);
    let anchor = locus
        .sequence_binding
        .as_mut()
        .unwrap()
        .genome_anchor
        .as_mut()
        .unwrap();
    anchor.genome_id = reference.assembly.clone();
    assert!(
        project(
            &record,
            &reference,
            &locus,
            &sequence,
            &"a".repeat(64),
            vec![],
            None,
        )
        .is_ok()
    );

    for incompatible in [
        "synthetic",
        "assembly",
        "synthetic-assembl",
        "Human",
        "Ensembl",
        "116",
    ] {
        locus
            .sequence_binding
            .as_mut()
            .unwrap()
            .genome_anchor
            .as_mut()
            .unwrap()
            .genome_id = incompatible.into();
        assert!(
            project(
                &record,
                &reference,
                &locus,
                &sequence,
                &"a".repeat(64),
                vec![],
                None,
            )
            .is_err(),
            "accepted incompatible anchor identifier {incompatible}"
        );
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
fn tss_context_distinct_and_overlapping_windows_keep_their_own_cutrun_intervals() {
    // One synthetic locus, four windows, two overlapping windows and an empty
    // window. Exercise the manifest join, both loaded/genomic strands and SVG.
    for strand in [TssStrand::Plus, TssStrand::Minus] {
        for loaded_reverse in [false, true] {
            let (record, reference, mut locus, sequence, _) = fixture(strand, loaded_reverse);
            let prototype = locus.occupancy_groups[0].lanes[0].lane.intervals[0].clone();
            locus.occupancy_groups[0].lanes[0].lane.intervals = [
                ("early", 105, 115, 2.0),
                ("shared", 125, 132, 7.0),
                ("late", 185, 191, 5.0),
            ]
            .into_iter()
            .map(|(id, start, end, score)| {
                let mut interval = prototype.clone();
                interval.interval_id = id.into();
                interval.genomic_start_1based = start;
                interval.genomic_end_1based = end;
                (interval.local_start_1based, interval.local_end_1based) = if loaded_reverse {
                    (220 - end, 220 - start)
                } else {
                    (start - 99, end - 99)
                };
                interval.score = Some(score);
                interval
            })
            .collect();
            let mut report = crate::tss_profile_export::tests::synthetic_report();
            report.reference = reference.clone();
            let template = report.windows[0].clone();
            report.windows = [110, 120, 180, 195]
                .into_iter()
                .enumerate()
                .map(|(index, start)| {
                    let mut window = template.clone();
                    window.record = record.clone();
                    window.record.promoter_id = format!("synthetic-window-{index}");
                    window.record.geometry = TssGeometry {
                        chromosome: record.geometry.chromosome.clone(),
                        strand,
                        start_1based: start,
                        end_1based: start + 24,
                        tss_1based: if strand == TssStrand::Plus {
                            start + 8
                        } else {
                            start + 16
                        },
                        upstream_bp: 8,
                        downstream_bp: 16,
                    };
                    let genomic = "ACGTTGAACCTA".repeat(10);
                    let bases = &genomic[(start - 100) as usize..(start - 75) as usize];
                    let bases = if strand == TssStrand::Minus {
                        GentleEngine::reverse_complement(bases)
                    } else {
                        bases.into()
                    };
                    window.record.sequence_sha256 = sha256_hex_bytes(bases.as_bytes());
                    window.selected = false;
                    window.selection_evidence = None;
                    for track in &mut window.tracks {
                        track.forward_scores = vec![Some(0.0); 26 - track.motif_length_bp];
                        track.reverse_scores = track.forward_scores.clone();
                        track.forward_maximum = None;
                        track.reverse_maximum = None;
                        track.forward_peaks.clear();
                        track.reverse_peaks.clear();
                    }
                    window.comparisons =
                        GentleEngine::tss_matrix_comparisons(&window, &report.panel_resolution);
                    window
                })
                .collect();
            let original = serde_json::to_value(&report).unwrap();
            let temp = tempfile::tempdir().unwrap();
            let write = |name: &str, bytes: Vec<u8>| {
                std::fs::write(temp.path().join(name), &bytes).unwrap();
                TssContextFile {
                    path: name.into(),
                    sha256: sha256_hex_bytes(&bytes),
                }
            };
            let manifest = TssContextManifest {
                schema: CONTEXT_INPUT_SCHEMA.into(),
                reference,
                genes: vec![TssContextSource {
                    transcript_annotation_sources: vec![],
                    gene_id: record.gene_id,
                    locus_report: write("locus.json", serde_json::to_vec(&locus).unwrap()),
                    locus_fasta: write(
                        "locus.fa",
                        format!(">synthetic\n{sequence}\n").into_bytes(),
                    ),
                    tata_report: None,
                }],
            };
            write("context.json", serde_json::to_vec(&manifest).unwrap());
            attach(&mut report, &temp.path().join("context.json"), &mut || true).unwrap();
            let mut expected = std::collections::BTreeMap::new();
            for (index, window) in report.windows.iter().enumerate() {
                let context = window.detail_context.as_ref().unwrap();
                let lane = &context.occupancy[0];
                let mut ids = lane
                    .intervals
                    .iter()
                    .map(|i| i.interval_id.as_str())
                    .collect::<Vec<_>>();
                ids.sort();
                assert_eq!(
                    ids,
                    [
                        vec!["early", "shared"],
                        vec!["shared"],
                        vec!["late"],
                        vec![]
                    ][index]
                );
                assert_eq!(lane.display_abs_max_score, 10.0);
                assert_eq!(
                    context.occupancy[1].state,
                    GeneLocusOccupancyLaneState::NotPrepared
                );
                assert!(context.occupancy[1].intervals.is_empty());
                for interval in &lane.intervals {
                    let start = interval
                        .span
                        .genomic_start_1based
                        .max(window.record.geometry.start_1based);
                    let end = interval
                        .span
                        .genomic_end_1based
                        .min(window.record.geometry.end_1based);
                    let (a, b) = if strand == TssStrand::Plus {
                        (
                            start - window.record.geometry.start_1based,
                            end - window.record.geometry.start_1based + 1,
                        )
                    } else {
                        (
                            window.record.geometry.end_1based - end,
                            window.record.geometry.end_1based - start + 1,
                        )
                    };
                    assert_eq!(
                        (
                            interval.span.start_0based,
                            interval.span.end_0based_exclusive
                        ),
                        (a as usize, b as usize)
                    );
                    let score = match interval.interval_id.as_str() {
                        "early" => 2.0,
                        "shared" => 7.0,
                        "late" => 5.0,
                        _ => unreachable!(),
                    };
                    assert_eq!(interval.score, Some(score));
                    expected.insert(
                        (
                            window.record.promoter_id.clone(),
                            interval.interval_id.clone(),
                        ),
                        (a, b, 22.0 * score / 10.0),
                    );
                }
            }
            for panels_per_page in [1, 2] {
                let pages = gentle_render::tss_profiles::render_tss_profile_pages(
                    &report,
                    &TssProfileRenderOptions {
                        panels_per_page,
                        ..Default::default()
                    },
                )
                .unwrap();
                let mut seen = BTreeSet::new();
                for page in pages {
                    let mut promoter = String::new();
                    for event in svg::read(&page.svg).unwrap() {
                        if let svg::parser::Event::Tag(_, _, attributes) = event {
                            let role = attributes
                                .get("data-role")
                                .map(ToString::to_string)
                                .unwrap_or_default();
                            if role == "tss-panel" {
                                promoter = attributes["data-promoter-id"].to_string();
                            }
                            if role == "context-occupancy-interval" {
                                let key =
                                    (promoter.clone(), attributes["data-interval-id"].to_string());
                                let (start, end, height) = expected[&key];
                                assert_eq!(
                                    attributes["data-start-0based"].to_string(),
                                    start.to_string()
                                );
                                assert_eq!(
                                    attributes["data-end-0based-exclusive"].to_string(),
                                    end.to_string()
                                );
                                assert!(
                                    (attributes["height"].to_string().parse::<f64>().unwrap()
                                        - height)
                                        .abs()
                                        < 1e-6
                                );
                                assert!(seen.insert(key));
                            }
                        }
                    }
                }
                assert_eq!(seen, expected.keys().cloned().collect());
            }
            let mut stripped = report.clone();
            for window in &mut stripped.windows {
                window.detail_context = None;
            }
            assert_eq!(
                serde_json::to_value(stripped).unwrap(),
                original,
                "context must not change TFBS scores or window identity"
            );
            let mut stale = report.clone();
            stale.windows[1].detail_context = stale.windows[0].detail_context.clone();
            assert!(
                gentle_render::tss_profiles::render_tss_profile_pages(&stale, &Default::default())
                    .is_err()
            );
        }
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
    let coverage = c.transcript_payload_coverage.as_ref().unwrap();
    assert_eq!(coverage.unassessed_transcript_ids, vec!["missing.tx"]);
    assert!(coverage.statement.contains("unassessed, not missing"));
    assert!(
        !c.warnings
            .iter()
            .any(|s| s.contains("missing.tx") || s.contains("transcript geometry unavailable"))
    );
}

#[test]
fn tss_context_source_join_runs_bound_gff3_adapters_and_preserves_signal_on_both_strands() {
    use gentle_protocol::transcript_presentation::*;
    for strand in [TssStrand::Plus, TssStrand::Minus] {
        for loaded_reverse in [false, true] {
            let (mut record, reference, locus, sequence, _) = fixture(strand, loaded_reverse);
            record
                .transcripts
                .extend(["other.1".into(), "omitted.1".into(), "omitted.2".into()]);
            let dir = tempfile::tempdir().unwrap();
            let locus_bytes = serde_json::to_vec(&locus).unwrap();
            std::fs::write(dir.path().join("locus.json"), &locus_bytes).unwrap();
            let fasta = format!(">synthetic\n{sequence}\n");
            std::fs::write(dir.path().join("locus.fa"), &fasta).unwrap();
            let mut annotations = Vec::new();
            for (provider, id) in [
                (TranscriptProvider::Ensembl, "synthetic.tx"),
                (TranscriptProvider::RefSeq, "other.1"),
            ] {
                let content = format!(
                    "##gff-version 3\n#!genome-build {}\nSynthetic\tannotation\tmRNA\t125\t194\t.\t{}\t.\tID=tx;Parent=gene;transcript_id={id}\nSynthetic\tannotation\texon\t125\t147\t.\t{}\t.\tID=e1;Parent=tx\nSynthetic\tannotation\texon\t183\t194\t.\t{}\t.\tID=e2;Parent=tx\n",
                    reference.assembly,
                    strand.as_str(),
                    strand.as_str(),
                    strand.as_str()
                );
                let name = format!("{id}.gff3");
                std::fs::write(dir.path().join(&name), &content).unwrap();
                annotations.push(TranscriptAnnotationSource {
                    path: name,
                    sha256: sha256_hex_bytes(content.as_bytes()),
                    provider,
                    format: TranscriptAnnotationFormat::Gff3,
                    assembly: reference.assembly.clone(),
                    release: format!("{provider:?}-release"),
                    accession: format!("{provider:?}-source"),
                    chromosome: "Synthetic".into(),
                    locus_sequence_sha256: sha256_hex_bytes(sequence.as_bytes()),
                    gene_ids: vec!["gene".into()],
                });
            }
            let manifest = TssContextManifest {
                schema: CONTEXT_INPUT_SCHEMA.into(),
                reference: reference.clone(),
                genes: vec![TssContextSource {
                    gene_id: record.gene_id.clone(),
                    locus_report: TssContextFile {
                        path: "locus.json".into(),
                        sha256: sha256_hex_bytes(&locus_bytes),
                    },
                    locus_fasta: TssContextFile {
                        path: "locus.fa".into(),
                        sha256: sha256_hex_bytes(fasta.as_bytes()),
                    },
                    tata_report: None,
                    transcript_annotation_sources: annotations,
                }],
            };
            let manifest_path = dir.path().join("context.json");
            std::fs::write(&manifest_path, serde_json::to_vec(&manifest).unwrap()).unwrap();
            let mut report = crate::tss_profile_export::tests::synthetic_report();
            report.reference = reference;
            report.windows.truncate(1);
            report.windows[0].record = record;
            let original = serde_json::to_vec(&report.windows[0].tracks).unwrap();
            attach(&mut report, &manifest_path, &mut || true).unwrap();
            assert_eq!(
                serde_json::to_vec(&report.windows[0].tracks).unwrap(),
                original
            );
            let c = report.windows[0].detail_context.as_ref().unwrap();
            let p = c.transcript_presentation.as_ref().unwrap();
            assert_eq!(p.physical_exons.len(), 2);
            assert_eq!(p.structure_groups.len(), 1);
            assert_eq!(p.records.len(), 2);
            assert_eq!(
                p.tss_ticks[0].genomic_position_1based,
                if strand == TssStrand::Plus { 125 } else { 194 }
            );
            assert_eq!(p.tss_deltas[0].transcript_oriented_delta_bp, 0);
            assert_eq!(
                c.transcript_payload_coverage
                    .as_ref()
                    .unwrap()
                    .unassessed_transcript_ids,
                vec!["omitted.1", "omitted.2"]
            );
            assert!(c.warnings.iter().all(|w| !w.contains("omitted")));
            assert!(!c.occupancy[0].intervals.is_empty());
            let untouched = std::fs::read(dir.path().join("locus.json")).unwrap();
            assert_eq!(untouched, locus_bytes);
        }
    }
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
            transcript_annotation_sources: vec![],
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
        genomic_motif_evidence: vec![],
        output_dir: output.to_str().unwrap().into(),
        context_manifest: Some(manifest_path.to_str().unwrap().into()),
        rendering: Default::default(),
        formats: vec![
            TssExportFormat::Svg,
            TssExportFormat::Png,
            TssExportFormat::Pdf,
            TssExportFormat::Genbank,
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
    let visible_text = svg::parser::Parser::new(&svg)
        .filter_map(|event| match event {
            svg::parser::Event::Text(text) => Some(text),
            _ => None,
        })
        .collect::<Vec<_>>()
        .join(" ");
    let visible_text = visible_text
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ");
    assert!(visible_text.contains("cropped"));
    assert!(visible_text.contains("not individual read ends"));
    let gb_name = receipt.outputs.keys().find(|n| n.ends_with(".gb")).unwrap();
    let gb = std::fs::read(output.join(gb_name)).unwrap();
    let records = gb_io::reader::SeqReader::new(gb.as_slice())
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(
        sha256_hex_bytes(
            &records[0]
                .seq
                .iter()
                .map(u8::to_ascii_uppercase)
                .collect::<Vec<_>>()
        ),
        enriched.windows[0].record.sequence_sha256
    );
    assert!(records[0].features.iter().any(|f| f.kind == "exon"));
    assert!(!records[0].features.iter().any(|f| f.kind == "CDS"));
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
        gb,
        std::fs::read(Path::new(&replay_request.output_dir).join(gb_name)).unwrap()
    );
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
