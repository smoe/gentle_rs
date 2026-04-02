//! Engine regression tests for the shared operation/state contract.
//!
//! These tests intentionally exercise the public engine surface across helper
//! module boundaries instead of re-testing each private extraction in isolation.
//!
//! Start here when you want deterministic examples of:
//! - complete `Operation` application flows
//! - report/export behavior that must stay adapter-neutral
//! - bug-regression coverage that spans multiple extracted engine submodules

use super::*;
use crate::genomes::BlastHit;
use crate::lineage_export::{LineageSvgNodeKind, build_lineage_svg_graph};
use bio::io::fasta;
use flate2::{Compression, write::GzEncoder};
use std::collections::HashMap;
use std::env;
use std::fs;
use std::io::Write;
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::path::Path;
use tempfile::tempdir;

const EXTERNAL_PRIMER_BINARY_TEST_ENV: &str = "GENTLE_TEST_EXTERNAL_BINARIES";

fn seq(s: &str) -> DNAsequence {
    DNAsequence::from_sequence(s).unwrap()
}

fn formula_test_sequence() -> DNAsequence {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(200)).expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "CDS".into(),
        location: gb_io::seq::Location::simple_range(20, 80),
        qualifiers: vec![("label".into(), Some("CDS_A".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(10, 40),
        qualifiers: vec![("label".into(), Some("TP73".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(80, 140),
        qualifiers: vec![("label".into(), Some("TP73".to_string()))],
    });
    dna
}

fn sequencing_confirmation_fixture_path(name: &str) -> String {
    format!(
        "{}/test_files/fixtures/sequencing_confirmation/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

fn sequencing_confirmation_junction_target(
    target_id: &str,
    start_0based: usize,
    end_0based_exclusive: usize,
    junction_left_end_0based: usize,
    label: &str,
) -> SequencingConfirmationTargetSpec {
    SequencingConfirmationTargetSpec {
        target_id: target_id.to_string(),
        label: label.to_string(),
        kind: SequencingConfirmationTargetKind::Junction,
        start_0based,
        end_0based_exclusive,
        junction_left_end_0based: Some(junction_left_end_0based),
        expected_bases: None,
        baseline_bases: None,
        required: true,
    }
}

fn sequencing_confirmation_trace_record(
    trace_id: &str,
    called_bases: &str,
    called_base_confidence_values: &[u8],
    peak_locations: &[u32],
) -> SequencingTraceRecord {
    SequencingTraceRecord {
        schema: SEQUENCING_TRACE_RECORD_SCHEMA.to_string(),
        trace_id: trace_id.to_string(),
        format: SequencingTraceFormat::Scf,
        source_path: format!("{trace_id}.scf"),
        imported_at_unix_ms: 1,
        seq_id: None,
        sample_name: Some(trace_id.to_string()),
        sample_well: None,
        run_name: None,
        machine_name: None,
        machine_model: None,
        called_bases: called_bases.to_string(),
        called_base_confidence_values: called_base_confidence_values.to_vec(),
        peak_locations: peak_locations.to_vec(),
        channel_data: vec![],
        channel_summaries: vec![],
        clip_start_base_index: None,
        clip_end_base_index_exclusive: None,
        comments_text: None,
    }
}

fn demo_blast_report() -> GenomeBlastReport {
    GenomeBlastReport {
        genome_id: "demo".to_string(),
        query_length: 12,
        max_hits: 25,
        task: "blastn-short".to_string(),
        blastn_executable: "blastn".to_string(),
        blast_db_prefix: "/tmp/demo".to_string(),
        command: vec![],
        hit_count: 2,
        hits: vec![
            BlastHit {
                subject_id: "chr1".to_string(),
                identity_percent: 99.0,
                alignment_length: 12,
                mismatches: 0,
                gap_opens: 0,
                query_start: 1,
                query_end: 12,
                subject_start: 100,
                subject_end: 111,
                evalue: 1e-9,
                bit_score: 42.0,
                query_coverage_percent: Some(100.0),
            },
            BlastHit {
                subject_id: "chr2".to_string(),
                identity_percent: 93.5,
                alignment_length: 10,
                mismatches: 1,
                gap_opens: 0,
                query_start: 2,
                query_end: 11,
                subject_start: 200,
                subject_end: 209,
                evalue: 1e-4,
                bit_score: 27.0,
                query_coverage_percent: Some(83.3),
            },
        ],
        warnings: vec![],
        stderr: String::new(),
        options_override_json: None,
        effective_options_json: None,
    }
}

#[test]
fn feature_coordinate_formula_helper_resolves_boundary_and_offset() {
    let dna = formula_test_sequence();

    let start = parse_required_usize_or_formula_text_on_sequence(
        &dna,
        "=CDS.start+10",
        "primer_design.roi_start_0based",
    )
    .expect("formula start");
    let end = parse_required_usize_or_formula_text_on_sequence(
        &dna,
        "=CDS.end-5",
        "primer_design.roi_end_0based",
    )
    .expect("formula end");

    assert_eq!(start, 30);
    assert_eq!(end, 75);
}

#[test]
fn feature_coordinate_formula_helper_supports_label_filter_and_occurrence() {
    let dna = formula_test_sequence();

    let start = parse_required_usize_or_formula_text_on_sequence(
        &dna,
        "=gene[label=TP73].start",
        "primer_design.roi_start_0based",
    )
    .expect("label formula");
    let second_end = parse_required_usize_or_formula_text_on_sequence(
        &dna,
        "=gene[2].end",
        "primer_design.roi_end_0based",
    )
    .expect("occurrence formula");

    assert_eq!(start, 10);
    assert_eq!(second_end, 140);
}

#[test]
fn feature_coordinate_formula_helper_resolves_selection_range_expression() {
    let dna = formula_test_sequence();

    let (start, end_exclusive) =
        resolve_selection_formula_range_0based_on_sequence(&dna, "=CDS.start+10 .. CDS.end-5")
            .expect("selection formula");

    assert_eq!((start, end_exclusive), (30, 75));
}

fn splicing_test_sequence() -> DNAsequence {
    let mut bases = vec![b'A'; 40];
    for (idx, base) in [
        (8usize, b'G'),
        (9, b'T'),
        (10, b'A'),
        (11, b'G'),
        (14, b'A'),
        (15, b'G'),
        (20, b'G'),
        (21, b'T'),
        (24, b'A'),
        (25, b'G'),
    ] {
        bases[idx] = base;
    }
    let sequence = String::from_utf8(bases).expect("valid ASCII DNA");
    let mut dna = DNAsequence::from_sequence(&sequence).expect("valid DNA sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(2, 8),
            gb_io::seq::Location::simple_range(12, 20),
            gb_io::seq::Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_1".to_string())),
            ("label".into(), Some("NM_TEST_1".to_string())),
        ],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(2, 8),
            gb_io::seq::Location::simple_range(16, 20),
            gb_io::seq::Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_TEST_2".to_string())),
            ("label".into(), Some("NM_TEST_2".to_string())),
        ],
    });
    dna
}

fn splicing_seed_feature_sequence() -> DNAsequence {
    let mut dna = splicing_test_sequence();
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(2, 34),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("label".into(), Some("GENE1".to_string())),
        ],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "CDS".into(),
        location: gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(2, 8),
            gb_io::seq::Location::simple_range(12, 20),
            gb_io::seq::Location::simple_range(26, 34),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("label".into(), Some("GENE1 CDS".to_string())),
            ("product".into(), Some("GENE1 protein".to_string())),
        ],
    });
    dna
}

fn splicing_misc_rna_sequence() -> DNAsequence {
    let mut dna = splicing_test_sequence();
    for feature in dna.features_mut().iter_mut() {
        feature.kind = "misc_RNA".into();
    }
    dna
}

fn splicing_multi_gene_test_sequence() -> DNAsequence {
    let mut bases = vec![b'A'; 120];
    let mut fill_range = |start_0: usize, end_0: usize, pattern: &[u8]| {
        for idx in start_0..end_0 {
            bases[idx] = pattern[(idx - start_0) % pattern.len()];
        }
    };
    fill_range(4, 20, b"ACGTCAGGTA");
    fill_range(30, 45, b"CCTAGGTACA");
    fill_range(35, 45, b"GGATCCGTTA");
    fill_range(60, 75, b"TTAACCGGCT");
    fill_range(85, 100, b"GGCCATTAGC");
    let sequence = String::from_utf8(bases).expect("valid DNA sequence");
    let mut dna = DNAsequence::from_sequence(&sequence).expect("valid DNA");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(5, 20),
            gb_io::seq::Location::simple_range(31, 45),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_GENE1_1".to_string())),
            ("label".into(), Some("NM_GENE1_1".to_string())),
        ],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(5, 20),
            gb_io::seq::Location::simple_range(36, 45),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("NM_GENE1_2".to_string())),
            ("label".into(), Some("NM_GENE1_2".to_string())),
        ],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(61, 75),
            gb_io::seq::Location::simple_range(86, 100),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE2".to_string())),
            ("transcript_id".into(), Some("NM_GENE2_1".to_string())),
            ("label".into(), Some("NM_GENE2_1".to_string())),
        ],
    });
    dna
}

fn build_rna_read_gene_support_test_engine() -> GentleEngine {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_multi_gene_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let features = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .clone();
    let seed_feature_id = features
        .iter()
        .enumerate()
        .find(|(idx, feature)| {
            feature.kind.to_string().eq_ignore_ascii_case("mRNA")
                && GentleEngine::splicing_group_label(feature, *idx).eq_ignore_ascii_case("GENE1")
        })
        .map(|(idx, _)| idx)
        .expect("GENE1 seed feature");
    let gene2_feature_id = features
        .iter()
        .enumerate()
        .find(|(idx, feature)| {
            feature.kind.to_string().eq_ignore_ascii_case("mRNA")
                && GentleEngine::splicing_group_label(feature, *idx).eq_ignore_ascii_case("GENE2")
        })
        .map(|(idx, _)| idx)
        .expect("GENE2 feature");
    let gene1_view = engine
        .build_splicing_expert_view(
            "seq_a",
            seed_feature_id,
            SplicingScopePreset::TargetGroupTargetStrand,
        )
        .expect("GENE1 splicing view");
    let gene2_view = engine
        .build_splicing_expert_view(
            "seq_a",
            gene2_feature_id,
            SplicingScopePreset::TargetGroupTargetStrand,
        )
        .expect("GENE2 splicing view");
    let gene1_tx_full = gene1_view
        .transcripts
        .iter()
        .find(|row| row.transcript_id == "NM_GENE1_1")
        .expect("GENE1 full transcript")
        .clone();
    let gene1_tx_skip = gene1_view
        .transcripts
        .iter()
        .find(|row| row.transcript_id == "NM_GENE1_2")
        .expect("GENE1 skip transcript")
        .clone();
    let gene2_tx = gene2_view
        .transcripts
        .iter()
        .find(|row| row.transcript_id == "NM_GENE2_1")
        .expect("GENE2 transcript")
        .clone();

    let transcript_length_bp = |lane: &SplicingTranscriptLane| -> usize {
        lane.exons
            .iter()
            .map(|exon| {
                exon.end_1based
                    .saturating_sub(exon.start_1based)
                    .saturating_add(1)
            })
            .sum::<usize>()
    };
    let gene1_tx_full_len = transcript_length_bp(&gene1_tx_full);
    let gene1_tx_skip_len = transcript_length_bp(&gene1_tx_skip);
    let gene2_tx_len = transcript_length_bp(&gene2_tx);
    let gene1_full_start = gene1_tx_full
        .exons
        .first()
        .expect("GENE1 full first exon")
        .start_1based;
    let gene1_full_end = gene1_tx_full
        .exons
        .last()
        .expect("GENE1 full last exon")
        .end_1based;
    let gene1_partial_end = gene1_tx_full
        .exons
        .first()
        .expect("GENE1 partial exon")
        .end_1based;
    let gene1_skip_start = gene1_tx_skip
        .exons
        .first()
        .expect("GENE1 skip first exon")
        .start_1based;
    let gene1_skip_end = gene1_tx_skip
        .exons
        .last()
        .expect("GENE1 skip last exon")
        .end_1based;
    let gene2_start = gene2_tx
        .exons
        .first()
        .expect("GENE2 first exon")
        .start_1based;
    let gene2_end = gene2_tx.exons.last().expect("GENE2 last exon").end_1based;

    engine
        .upsert_rna_read_report(RnaReadInterpretationReport {
            schema: "gentle.rna_read_report.v1".to_string(),
            report_id: "rna_reads_gene_support".to_string(),
            seq_id: "seq_a".to_string(),
            seed_feature_id,
            scope: SplicingScopePreset::TargetGroupTargetStrand,
            origin_mode: RnaReadOriginMode::MultiGeneSparse,
            target_gene_ids: vec!["GENE2".to_string()],
            align_config: RnaReadAlignConfig {
                min_identity_fraction: 0.60,
                ..RnaReadAlignConfig::default()
            },
            hits: vec![
                RnaReadInterpretationHit {
                    record_index: 0,
                    header_id: "gene1_full".to_string(),
                    sequence: "A".repeat(gene1_tx_full_len),
                    read_length_bp: gene1_tx_full_len,
                    best_mapping: Some(RnaReadMappingHit {
                        transcript_feature_id: gene1_tx_full.transcript_feature_id,
                        transcript_id: gene1_tx_full.transcript_id.clone(),
                        transcript_label: gene1_tx_full.label.clone(),
                        strand: gene1_tx_full.strand.clone(),
                        target_start_1based: gene1_full_start,
                        target_end_1based: gene1_full_end,
                        target_start_offset_0based: 0,
                        target_end_offset_0based_exclusive: gene1_tx_full_len,
                        score: 200,
                        identity_fraction: 0.99,
                        query_coverage_fraction: 1.0,
                        ..RnaReadMappingHit::default()
                    }),
                    ..RnaReadInterpretationHit::default()
                },
                RnaReadInterpretationHit {
                    record_index: 1,
                    header_id: "gene1_fragment".to_string(),
                    sequence: "A".repeat(
                        gene1_tx_full
                            .exons
                            .first()
                            .expect("GENE1 first exon")
                            .end_1based
                            .saturating_sub(
                                gene1_tx_full
                                    .exons
                                    .first()
                                    .expect("GENE1 first exon")
                                    .start_1based,
                            )
                            .saturating_add(1),
                    ),
                    read_length_bp: gene1_tx_full
                        .exons
                        .first()
                        .expect("GENE1 first exon")
                        .end_1based
                        .saturating_sub(
                            gene1_tx_full
                                .exons
                                .first()
                                .expect("GENE1 first exon")
                                .start_1based,
                        )
                        .saturating_add(1),
                    best_mapping: Some(RnaReadMappingHit {
                        transcript_feature_id: gene1_tx_full.transcript_feature_id,
                        transcript_id: gene1_tx_full.transcript_id.clone(),
                        transcript_label: gene1_tx_full.label.clone(),
                        strand: gene1_tx_full.strand.clone(),
                        target_start_1based: gene1_full_start,
                        target_end_1based: gene1_partial_end,
                        target_start_offset_0based: 0,
                        target_end_offset_0based_exclusive: gene1_tx_full
                            .exons
                            .first()
                            .expect("GENE1 first exon")
                            .end_1based
                            .saturating_sub(
                                gene1_tx_full
                                    .exons
                                    .first()
                                    .expect("GENE1 first exon")
                                    .start_1based,
                            )
                            .saturating_add(1),
                        score: 120,
                        identity_fraction: 0.97,
                        query_coverage_fraction: 0.60,
                        ..RnaReadMappingHit::default()
                    }),
                    ..RnaReadInterpretationHit::default()
                },
                RnaReadInterpretationHit {
                    record_index: 2,
                    header_id: "gene1_skip_complete".to_string(),
                    sequence: "A".repeat(gene1_tx_skip_len),
                    read_length_bp: gene1_tx_skip_len,
                    best_mapping: Some(RnaReadMappingHit {
                        transcript_feature_id: gene1_tx_skip.transcript_feature_id,
                        transcript_id: gene1_tx_skip.transcript_id.clone(),
                        transcript_label: gene1_tx_skip.label.clone(),
                        strand: gene1_tx_skip.strand.clone(),
                        target_start_1based: gene1_skip_start,
                        target_end_1based: gene1_skip_end,
                        target_start_offset_0based: 0,
                        target_end_offset_0based_exclusive: gene1_tx_skip_len,
                        score: 170,
                        identity_fraction: 0.50,
                        query_coverage_fraction: 1.0,
                        ..RnaReadMappingHit::default()
                    }),
                    ..RnaReadInterpretationHit::default()
                },
                RnaReadInterpretationHit {
                    record_index: 3,
                    header_id: "gene2_full".to_string(),
                    sequence: "A".repeat(gene2_tx_len),
                    read_length_bp: gene2_tx_len,
                    best_mapping: Some(RnaReadMappingHit {
                        transcript_feature_id: gene2_tx.transcript_feature_id,
                        transcript_id: gene2_tx.transcript_id.clone(),
                        transcript_label: gene2_tx.label.clone(),
                        strand: gene2_tx.strand.clone(),
                        target_start_1based: gene2_start,
                        target_end_1based: gene2_end,
                        target_start_offset_0based: 0,
                        target_end_offset_0based_exclusive: gene2_tx_len,
                        score: 190,
                        identity_fraction: 0.98,
                        query_coverage_fraction: 1.0,
                        ..RnaReadMappingHit::default()
                    }),
                    ..RnaReadInterpretationHit::default()
                },
            ],
            ..RnaReadInterpretationReport::default()
        })
        .expect("upsert gene-support report");
    engine
}

fn find_gene_support_exon_row<'a>(
    rows: &'a [RnaReadGeneExonSupportRow],
    gene_id: &str,
    exon_ordinal: usize,
) -> Option<&'a RnaReadGeneExonSupportRow> {
    rows.iter()
        .find(|row| row.gene_id.eq_ignore_ascii_case(gene_id) && row.exon_ordinal == exon_ordinal)
}

fn find_gene_support_pair_row<'a>(
    rows: &'a [RnaReadGeneExonPairSupportRow],
    gene_id: &str,
    from_exon_ordinal: usize,
    to_exon_ordinal: usize,
) -> Option<&'a RnaReadGeneExonPairSupportRow> {
    rows.iter().find(|row| {
        row.gene_id.eq_ignore_ascii_case(gene_id)
            && row.from_exon_ordinal == from_exon_ordinal
            && row.to_exon_ordinal == to_exon_ordinal
    })
}

fn synth_oligo(desc: &str, sequence: &[u8]) -> DNAsequence {
    let record = fasta::Record::with_attrs("synthetic", Some(desc), sequence);
    DNAsequence::from_fasta_record(&record)
}

fn assert_fasta_roundtrip(expected: DNAsequence) {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("synth".to_string(), expected.clone());
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("fa");
    let path_text = path.display().to_string();

    engine
        .apply(Operation::SaveFile {
            seq_id: "synth".to_string(),
            path: path_text.clone(),
            format: ExportFormat::Fasta,
        })
        .unwrap();
    engine
        .apply(Operation::LoadFile {
            path: path_text,
            as_id: Some("roundtrip".to_string()),
        })
        .unwrap();

    let actual = engine
        .state()
        .sequences
        .get("roundtrip")
        .expect("roundtrip sequence should exist");
    assert_eq!(actual.assert_sequence_equality(&expected), Ok(()));
}

fn write_gzip(path: &Path, text: &str) {
    let file = std::fs::File::create(path).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(text.as_bytes()).unwrap();
    encoder.finish().unwrap();
}

fn write_multi_member_gzip(path: &Path, members: &[&str]) {
    std::fs::File::create(path).unwrap();
    for member in members {
        let file = std::fs::OpenOptions::new().append(true).open(path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(member.as_bytes()).unwrap();
        encoder.finish().unwrap();
    }
}

#[cfg(unix)]
fn install_fake_rnapkin(path: &Path) -> String {
    let script_path = path.join("fake_rnapkin.sh");
    let script = r#"#!/bin/sh
if [ "$1" = "-v" ] && [ "$2" = "-p" ]; then
  seq="$3"
  echo "rnapkin textual report"
  echo "sequence_length=${#seq}"
  echo "points:"
  echo "0.0 0.0"
  echo "1.0 1.0"
  exit 0
fi

if [ "$#" -eq 2 ]; then
  seq="$1"
  out="$2"
  cat > "$out" <<EOF
<svg xmlns="http://www.w3.org/2000/svg" width="160" height="80"><text x="10" y="40">$seq</text></svg>
EOF
  echo "wrote $out" >&2
  exit 0
fi

echo "unexpected args: $@" >&2
exit 2
"#;
    std::fs::write(&script_path, script).expect("write fake rnapkin");
    let mut perms = std::fs::metadata(&script_path)
        .expect("metadata fake rnapkin")
        .permissions();
    perms.set_mode(0o755);
    std::fs::set_permissions(&script_path, perms).expect("chmod fake rnapkin");
    script_path.display().to_string()
}

#[cfg(unix)]
fn install_fake_bigwig_to_bedgraph(path: &Path, bedgraph_source: &Path) -> String {
    let script_path = path.join("fake_bigwig_to_bedgraph.sh");
    let script = format!(
        "#!/bin/sh\nif [ \"$#\" -ne 2 ]; then\n  echo \"expected INPUT.bw OUTPUT.bedGraph\" >&2\n  exit 2\nfi\ncp \"{}\" \"$2\"\n",
        bedgraph_source.display()
    );
    std::fs::write(&script_path, script).expect("write fake bigwig converter");
    let mut perms = std::fs::metadata(&script_path)
        .expect("metadata fake bigwig converter")
        .permissions();
    perms.set_mode(0o755);
    std::fs::set_permissions(&script_path, perms).expect("chmod fake bigwig converter");
    script_path.display().to_string()
}

#[cfg(unix)]
fn install_fake_primer3(path: &Path) -> String {
    let script_path = path.join("fake_primer3.sh");
    let script = r#"#!/bin/sh
if [ "$1" = "--version" ]; then
  echo "primer3_core synthetic-fixture 2.6.1"
  exit 0
fi
echo "unexpected args: $@" >&2
exit 2
"#;
    std::fs::write(&script_path, script).expect("write fake primer3");
    let mut perms = std::fs::metadata(&script_path)
        .expect("metadata fake primer3")
        .permissions();
    perms.set_mode(0o755);
    std::fs::set_permissions(&script_path, perms).expect("chmod fake primer3");
    script_path.display().to_string()
}

#[cfg(unix)]
fn install_fake_primer3_zero_pairs(path: &Path) -> (String, String) {
    let script_path = path.join("fake_primer3_zero_pairs.sh");
    let request_capture_path = path.join("fake_primer3_request_capture.txt");
    let script = format!(
        r#"#!/bin/sh
if [ "$1" = "--version" ]; then
  echo "primer3_core synthetic-fixture 2.6.1"
  exit 0
fi
cat > "{}"
echo "PRIMER_PAIR_NUM_RETURNED=0"
echo "PRIMER_LEFT_NUM_RETURNED=0"
echo "PRIMER_RIGHT_NUM_RETURNED=0"
echo "PRIMER_PAIR_EXPLAIN=considered 16, no acceptable primer pairs"
echo "="
"#,
        request_capture_path.display()
    );
    std::fs::write(&script_path, script).expect("write fake primer3 zero-pairs");
    let mut perms = std::fs::metadata(&script_path)
        .expect("metadata fake primer3 zero-pairs")
        .permissions();
    perms.set_mode(0o755);
    std::fs::set_permissions(&script_path, perms).expect("chmod fake primer3 zero-pairs");
    (
        script_path.display().to_string(),
        request_capture_path.display().to_string(),
    )
}

fn file_url(path: &Path) -> String {
    format!("file://{}", path.display())
}

struct EnvVarGuard {
    key: &'static str,
    previous: Option<String>,
}

fn candidate_store_env_lock() -> &'static std::sync::Mutex<()> {
    use std::sync::{Mutex, OnceLock};
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(()))
}

impl EnvVarGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let previous = env::var(key).ok();
        unsafe {
            env::set_var(key, value);
        }
        Self { key, previous }
    }
}

impl Drop for EnvVarGuard {
    fn drop(&mut self) {
        match &self.previous {
            Some(value) => unsafe {
                env::set_var(self.key, value);
            },
            None => unsafe {
                env::remove_var(self.key);
            },
        }
    }
}

#[test]
fn test_set_display_visibility() {
    let mut engine = GentleEngine::new();
    assert!(!engine.state().display.show_tfbs);
    assert!(engine.state().display.show_cds_features);
    let res = engine
        .apply(Operation::SetDisplayVisibility {
            target: DisplayTarget::Features,
            visible: false,
        })
        .unwrap();
    assert!(res.messages.iter().any(|m| m.contains("features")));
    assert!(!engine.state().display.show_features);
    let res_tfbs = engine
        .apply(Operation::SetDisplayVisibility {
            target: DisplayTarget::Tfbs,
            visible: true,
        })
        .unwrap();
    assert!(res_tfbs.messages.iter().any(|m| m.contains("tfbs")));
    assert!(engine.state().display.show_tfbs);

    let res_cds = engine
        .apply(Operation::SetDisplayVisibility {
            target: DisplayTarget::CdsFeatures,
            visible: false,
        })
        .unwrap();
    assert!(res_cds.messages.iter().any(|m| m.contains("cds_features")));
    assert!(!engine.state().display.show_cds_features);
}

#[test]
fn test_set_linear_viewport() {
    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::SetLinearViewport {
            start_bp: 123,
            span_bp: 456,
        })
        .unwrap();
    assert!(res.messages.iter().any(|m| m.contains("linear viewport")));
    assert_eq!(engine.state().display.linear_view_start_bp, 123);
    assert_eq!(engine.state().display.linear_view_span_bp, 456);
}

#[test]
fn test_extract_region() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("x".to_string(), seq(&"ATGC".repeat(100)));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::ExtractRegion {
            input: "x".to_string(),
            from: 2,
            to: 7,
            output_id: Some("part".to_string()),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids, vec!["part".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("part")
            .unwrap()
            .get_forward_string(),
        "GCATG"
    );
}

#[test]
fn test_extract_region_preserves_overlapping_features() {
    let mut state = ProjectState::default();
    let mut dna = seq("ATGCATGCATGC");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(2, 8),
        qualifiers: vec![("label".into(), Some("gene_a".to_string()))],
    });
    state.sequences.insert("x".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::ExtractRegion {
            input: "x".to_string(),
            from: 2,
            to: 10,
            output_id: Some("part".to_string()),
        })
        .unwrap();
    let part = engine.state().sequences.get("part").expect("part sequence");
    let gene = part
        .features()
        .iter()
        .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
        .expect("gene should be preserved in extracted region");
    assert_eq!(gene.location.find_bounds().expect("gene bounds"), (0, 6));
}

#[test]
fn test_extract_anchored_region_with_constraints_unique() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq("TTTGAATTCACGTACGTACGTTAAACCC"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::ExtractAnchoredRegion {
            input: "tpl".to_string(),
            anchor: AnchoredRegionAnchor::Position { zero_based: 21 },
            direction: AnchorDirection::Upstream,
            target_length_bp: 18,
            length_tolerance_bp: 2,
            required_re_sites: vec!["EcoRI".to_string()],
            required_tf_motifs: vec!["CACGTA".to_string()],
            forward_primer: Some("GAATTC".to_string()),
            reverse_primer: Some("ACGT".to_string()),
            output_prefix: Some("prom".to_string()),
            unique: Some(true),
            max_candidates: Some(5),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids, vec!["prom_1".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("prom_1")
            .unwrap()
            .get_forward_string(),
        "GAATTCACGTACGTACGT"
    );
}

#[test]
fn test_extract_anchored_region_from_feature_boundary() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::LoadFile {
            path: "test_files/pGEX-3X.gb".to_string(),
            as_id: Some("pgex".to_string()),
        })
        .unwrap();
    let res = engine
        .apply(Operation::ExtractAnchoredRegion {
            input: "pgex".to_string(),
            anchor: AnchoredRegionAnchor::FeatureBoundary {
                feature_kind: Some("CDS".to_string()),
                feature_label: None,
                boundary: AnchorBoundary::Start,
                occurrence: Some(0),
            },
            direction: AnchorDirection::Upstream,
            target_length_bp: 60,
            length_tolerance_bp: 0,
            required_re_sites: vec![],
            required_tf_motifs: vec![],
            forward_primer: None,
            reverse_primer: None,
            output_prefix: Some("cds_up".to_string()),
            unique: Some(true),
            max_candidates: Some(1),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids, vec!["cds_up_1".to_string()]);
    assert_eq!(engine.state().sequences.get("cds_up_1").unwrap().len(), 60);
}

#[test]
fn test_ligation_simple_concatenation() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("a".to_string(), seq(&"ATGC".repeat(40)));
    state
        .sequences
        .insert("b".to_string(), seq(&"TTAA".repeat(40)));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::Ligation {
            inputs: vec!["a".to_string(), "b".to_string()],
            circularize_if_possible: false,
            output_id: None,
            protocol: LigationProtocol::Blunt,
            output_prefix: Some("ab".to_string()),
            unique: None,
        })
        .unwrap();
    assert_eq!(res.created_seq_ids.len(), 2);
    let out = engine.state().sequences.get("ab_1").unwrap();
    assert_eq!(
        out.get_forward_string(),
        format!("{}{}", "ATGC".repeat(40), "TTAA".repeat(40))
    );
    assert!(!out.is_circular());
}

#[test]
fn test_merge_containers_creates_pool_copies() {
    let mut state = ProjectState::default();
    state.sequences.insert("a".to_string(), seq("ATGC"));
    state.sequences.insert("b".to_string(), seq("TTAA"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::MergeContainers {
            inputs: vec!["a".to_string(), "b".to_string()],
            output_prefix: Some("pool".to_string()),
        })
        .unwrap();
    assert_eq!(
        res.created_seq_ids,
        vec!["pool_1".to_string(), "pool_2".to_string()]
    );
    assert_eq!(
        engine
            .state()
            .sequences
            .get("pool_1")
            .unwrap()
            .get_forward_string(),
        "ATGC"
    );
    assert_eq!(
        engine
            .state()
            .sequences
            .get("pool_2")
            .unwrap()
            .get_forward_string(),
        "TTAA"
    );
}

#[test]
fn test_merge_containers_respects_max_fragments() {
    let mut state = ProjectState::default();
    state.parameters.max_fragments_per_container = 1;
    state.sequences.insert("a".to_string(), seq("ATGC"));
    state.sequences.insert("b".to_string(), seq("TTAA"));
    let mut engine = GentleEngine::from_state(state);
    let err = engine
        .apply(Operation::MergeContainers {
            inputs: vec!["a".to_string(), "b".to_string()],
            output_prefix: Some("pool".to_string()),
        })
        .unwrap_err();
    assert!(err.message.contains("max_fragments_per_container"));
}

#[test]
fn test_ligation_protocol_blunt_enumerates_ordered_pairs() {
    let mut state = ProjectState::default();
    state.sequences.insert("a".to_string(), seq("ATGC"));
    state.sequences.insert("b".to_string(), seq("TTAA"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::Ligation {
            inputs: vec!["a".to_string(), "b".to_string()],
            circularize_if_possible: false,
            output_id: None,
            protocol: LigationProtocol::Blunt,
            output_prefix: Some("lig".to_string()),
            unique: Some(false),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids.len(), 2);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("lig_1")
            .unwrap()
            .get_forward_string(),
        "ATGCTTAA"
    );
    assert_eq!(
        engine
            .state()
            .sequences
            .get("lig_2")
            .unwrap()
            .get_forward_string(),
        "TTAAATGC"
    );
}

#[test]
fn test_ligation_protocol_sticky_uses_overhang_compatibility() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("x".to_string(), seq("ATGGATCCGCATGGATCCGCATGGATCCGC"));
    let mut engine = GentleEngine::from_state(state);
    let digest_res = engine
        .apply(Operation::Digest {
            input: "x".to_string(),
            enzymes: vec!["BamHI".to_string()],
            output_prefix: Some("frag".to_string()),
        })
        .unwrap();
    assert!(digest_res.created_seq_ids.len() >= 2);
    let a = digest_res.created_seq_ids[0].clone();
    let b = digest_res.created_seq_ids[1].clone();

    let lig_res = engine
        .apply(Operation::Ligation {
            inputs: vec![a, b],
            circularize_if_possible: false,
            output_id: None,
            protocol: LigationProtocol::Sticky,
            output_prefix: Some("st".to_string()),
            unique: Some(false),
        })
        .unwrap();
    assert!(!lig_res.created_seq_ids.is_empty());
}

#[test]
fn test_workflow_digest_merge_ligation_is_deterministic() {
    let mut base = ProjectState::default();
    base.sequences
        .insert("x".to_string(), seq("ATGGATCCGCATGGATCCGCATGGATCCGC"));

    let run_once = |state: ProjectState| {
        let mut engine = GentleEngine::from_state(state);
        let digest = engine
            .apply(Operation::Digest {
                input: "x".to_string(),
                enzymes: vec!["BamHI".to_string(), "EcoRI".to_string()],
                output_prefix: Some("d".to_string()),
            })
            .unwrap();
        let merge = engine
            .apply(Operation::MergeContainers {
                inputs: digest.created_seq_ids.clone(),
                output_prefix: Some("m".to_string()),
            })
            .unwrap();
        let lig = engine
            .apply(Operation::Ligation {
                inputs: merge.created_seq_ids.clone(),
                circularize_if_possible: false,
                output_id: None,
                protocol: LigationProtocol::Sticky,
                output_prefix: Some("lig".to_string()),
                unique: Some(false),
            })
            .unwrap();
        lig.created_seq_ids
    };

    let a = run_once(base.clone());
    let b = run_once(base.clone());
    assert_eq!(a, b);
    assert!(!a.is_empty());
    assert_eq!(a.first().unwrap(), "lig_1");
}

#[test]
fn test_from_state_reseeds_operation_counter_from_existing_project_ops() {
    let mut state = ProjectState::default();
    state.sequences.insert("seed".to_string(), seq("ATGCATGC"));
    state.lineage.next_node_counter = 1;
    state.lineage.nodes.insert(
        "n-1".to_string(),
        LineageNode {
            node_id: "n-1".to_string(),
            seq_id: "seed".to_string(),
            created_by_op: Some("op-41".to_string()),
            origin: SequenceOrigin::ImportedSynthetic,
            created_at_unix_ms: 1,
        },
    );
    state
        .lineage
        .seq_to_node
        .insert("seed".to_string(), "n-1".to_string());
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("seed".to_string()),
            members: vec!["seed".to_string()],
            created_by_op: Some("op-57".to_string()),
            created_at_unix_ms: 1,
        },
    );
    state
        .container_state
        .seq_to_latest_container
        .insert("seed".to_string(), "container-1".to_string());
    state.container_state.next_container_counter = 1;

    let mut engine = GentleEngine::from_state(state);
    let result = engine
        .apply(Operation::SetParameter {
            name: "max_fragments_per_container".to_string(),
            value: json!(123),
        })
        .expect("set parameter");
    assert_eq!(result.op_id, "op-58");
}

#[test]
fn test_design_primer_pairs_persists_report() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        seq(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        ),
    );
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;
    let result = engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                non_annealing_5prime_tail: None,
                fixed_5prime: None,
                fixed_3prime: None,
                required_motifs: vec![],
                forbidden_motifs: vec![],
                locked_positions: vec![],
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(90),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                non_annealing_5prime_tail: None,
                fixed_5prime: None,
                fixed_3prime: None,
                required_motifs: vec![],
                forbidden_motifs: vec![],
                locked_positions: vec![],
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("tp73_roi".to_string()),
        })
        .expect("design primer pairs");
    assert!(
        result
            .messages
            .iter()
            .any(|line| line.contains("primer-design report"))
    );
    let report = engine
        .get_primer_design_report("tp73_roi")
        .expect("report by id");
    assert_eq!(report.report_id, "tp73_roi");
    assert_eq!(report.template, "tpl");
    assert!(!report.pairs.is_empty());
    assert_eq!(
        result.created_seq_ids.len(),
        report.pair_count.saturating_mul(3)
    );
    assert!(
        result
            .created_seq_ids
            .iter()
            .any(|seq_id| seq_id.ends_with("_amplicon"))
    );
    assert!(
        result
            .created_seq_ids
            .iter()
            .all(|seq_id| engine.state().sequences.contains_key(seq_id))
    );
    let lineage = &engine.state().lineage;
    let template_node = lineage
        .seq_to_node
        .get("tpl")
        .expect("template lineage node should exist")
        .clone();
    for seq_id in &result.created_seq_ids {
        let node_id = lineage
            .seq_to_node
            .get(seq_id)
            .expect("primer lineage node should exist")
            .clone();
        let node = lineage
            .nodes
            .get(&node_id)
            .expect("primer lineage node payload should exist");
        assert_eq!(node.created_by_op.as_deref(), Some(result.op_id.as_str()));
        assert!(
            lineage.edges.iter().any(|edge| {
                edge.from_node_id == template_node
                    && edge.to_node_id == node_id
                    && edge.op_id == result.op_id
            }),
            "missing lineage edge tpl -> {seq_id} for op {}",
            result.op_id
        );
        let container_id = engine
            .state()
            .container_state
            .seq_to_latest_container
            .get(seq_id)
            .expect("primer sequence should have a latest container");
        let container = engine
            .state()
            .container_state
            .containers
            .get(container_id)
            .expect("primer container should exist");
        assert_eq!(
            container.created_by_op.as_deref(),
            Some(result.op_id.as_str())
        );
        assert!(matches!(container.kind, ContainerKind::Pool));
        assert_eq!(container.members.len(), 3);
        assert!(
            container
                .name
                .as_deref()
                .unwrap_or("")
                .starts_with("Primer pair tp73_roi r")
        );
    }
    let created_containers = engine
        .state()
        .container_state
        .containers
        .values()
        .filter(|container| container.created_by_op.as_deref() == Some(result.op_id.as_str()))
        .count();
    assert_eq!(created_containers, report.pair_count);
    assert_eq!(report.backend.requested, "internal");
    assert_eq!(report.backend.used, "internal");
    let listed = engine.list_primer_design_reports();
    assert!(listed.iter().any(|row| row.report_id == "tp73_roi"));
}

#[test]
fn test_design_primer_pairs_auto_backend_falls_back_to_internal() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        seq(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        ),
    );
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Auto;
    engine.state_mut().parameters.primer3_executable =
        "/definitely/missing/primer3_core".to_string();
    let result = engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                non_annealing_5prime_tail: None,
                fixed_5prime: None,
                fixed_3prime: None,
                required_motifs: vec![],
                forbidden_motifs: vec![],
                locked_positions: vec![],
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(90),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                non_annealing_5prime_tail: None,
                fixed_5prime: None,
                fixed_3prime: None,
                required_motifs: vec![],
                forbidden_motifs: vec![],
                locked_positions: vec![],
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("tp73_roi_auto".to_string()),
        })
        .expect("design primer pairs");
    assert!(
        result
            .warnings
            .iter()
            .any(|line| line.contains("Primer3 backend unavailable"))
    );
    let report = engine
        .get_primer_design_report("tp73_roi_auto")
        .expect("report by id");
    assert_eq!(report.backend.requested, "auto");
    assert_eq!(report.backend.used, "internal");
    assert!(report.backend.fallback_reason.is_some());
    assert!(!report.pairs.is_empty());
}

#[test]
fn test_design_primer_pairs_internal_reports_pair_evaluation_budget_truncation() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq(&"ACGT".repeat(180)));
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;
    let side = PrimerDesignSideConstraint {
        min_length: 18,
        max_length: 22,
        location_0based: None,
        start_0based: None,
        end_0based: None,
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 10_000,
        non_annealing_5prime_tail: None,
        fixed_5prime: None,
        fixed_3prime: None,
        required_motifs: vec![],
        forbidden_motifs: vec![],
        locked_positions: vec![],
    };
    let result = engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 220,
            roi_end_0based: 320,
            forward: side.clone(),
            reverse: side,
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 80,
            max_amplicon_bp: 360,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(5),
            report_id: Some("budget_capped_pairs".to_string()),
        })
        .expect("internal primer-pair design should complete");

    let report = engine
        .get_primer_design_report("budget_capped_pairs")
        .expect("report by id");
    assert_eq!(report.backend.used, "internal");
    assert!(report.rejection_summary.pair_evaluation_limit_skipped > 0);
    assert!(report.pair_count <= 5);
    assert!(
        result
            .warnings
            .iter()
            .any(|line| line.contains("evaluation limit")),
        "expected warning about internal pair-evaluation limit, got {:?}",
        result.warnings
    );
}

#[test]
fn test_design_qpcr_assays_internal_reports_pair_evaluation_budget_truncation() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq(&"ACGT".repeat(180)));
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;
    let side = PrimerDesignSideConstraint {
        min_length: 18,
        max_length: 22,
        location_0based: None,
        start_0based: None,
        end_0based: None,
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 10_000,
        non_annealing_5prime_tail: None,
        fixed_5prime: None,
        fixed_3prime: None,
        required_motifs: vec![],
        forbidden_motifs: vec![],
        locked_positions: vec![],
    };
    let probe = PrimerDesignSideConstraint {
        location_0based: Some(260),
        ..side.clone()
    };
    let result = engine
        .apply(Operation::DesignQpcrAssays {
            template: "tpl".to_string(),
            roi_start_0based: 220,
            roi_end_0based: 320,
            forward: side.clone(),
            reverse: side,
            probe,
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 80,
            max_amplicon_bp: 360,
            max_tm_delta_c: Some(100.0),
            max_probe_tm_delta_c: Some(100.0),
            max_assays: Some(5),
            report_id: Some("budget_capped_qpcr".to_string()),
        })
        .expect("internal qPCR design should complete");

    let report = engine
        .get_qpcr_design_report("budget_capped_qpcr")
        .expect("qPCR report by id");
    assert_eq!(report.backend.used, "internal");
    assert!(
        report
            .rejection_summary
            .primer_pair
            .pair_evaluation_limit_skipped
            > 0
    );
    assert!(report.assay_count <= 5);
    assert!(
        result
            .warnings
            .iter()
            .any(|line| line.contains("evaluation limit")),
        "expected warning about internal pair-evaluation limit, got {:?}",
        result.warnings
    );
}

#[test]
fn test_design_primer_pairs_primer3_backend_requires_executable() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        seq("ATGCGTACGATCGTAGCTAGCTAGCTAGCATCGATCGATGCGTACGATCGTAGCTAGCTAGCTAGCATCGATCG"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Primer3;
    engine.state_mut().parameters.primer3_executable =
        "/definitely/missing/primer3_core".to_string();
    let err = engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 10,
            roi_end_0based: 30,
            forward: PrimerDesignSideConstraint::default(),
            reverse: PrimerDesignSideConstraint::default(),
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 20,
            max_amplicon_bp: 80,
            max_tm_delta_c: Some(10.0),
            max_pairs: Some(10),
            report_id: Some("tp73_roi_primer3".to_string()),
        })
        .expect_err("missing primer3 executable should fail");
    assert!(err.message.contains("Primer3 backend executable"));
}

#[cfg(unix)]
#[test]
fn test_design_primer_pairs_primer3_zero_pairs_persists_request_and_explain() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        seq("ATGCGTACGATCGTAGCTAGCTAGCTAGCATCGATCGATGCGTACGATCGTAGCTAGCTAGCTAGCATCGATCG"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Primer3;
    let tmp = tempdir().expect("tempdir");
    let (fake_primer3, request_capture) = install_fake_primer3_zero_pairs(tmp.path());
    engine.state_mut().parameters.primer3_executable = fake_primer3.clone();

    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 10,
            roi_end_0based: 30,
            forward: PrimerDesignSideConstraint::default(),
            reverse: PrimerDesignSideConstraint::default(),
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 20,
            max_amplicon_bp: 80,
            max_tm_delta_c: Some(10.0),
            max_pairs: Some(10),
            report_id: Some("tp73_roi_primer3_zero".to_string()),
        })
        .expect("fake primer3 zero-pair design should complete");

    let report = engine
        .get_primer_design_report("tp73_roi_primer3_zero")
        .expect("primer3 zero-pair report");
    assert_eq!(report.backend.used, "primer3");
    assert_eq!(report.pair_count, 0);
    assert!(
        report
            .backend
            .primer3_explain
            .as_deref()
            .unwrap_or_default()
            .contains("PRIMER_PAIR_EXPLAIN")
    );
    let request_payload = report
        .backend
        .primer3_request_boulder_io
        .as_deref()
        .expect("request payload should be captured");
    assert!(request_payload.contains("SEQUENCE_TEMPLATE="));
    assert!(request_payload.contains("PRIMER_EXPLAIN_FLAG=1"));

    let captured = std::fs::read_to_string(&request_capture).expect("captured stdin payload");
    assert_eq!(captured, request_payload);
}

#[cfg(unix)]
#[test]
fn test_primer3_preflight_report_success() {
    let mut engine = GentleEngine::new();
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Primer3;
    let tmp = tempdir().expect("tempdir");
    let fake_primer3 = install_fake_primer3(tmp.path());
    let report = engine.primer3_preflight_report(None, Some(fake_primer3.as_str()));
    assert_eq!(report.backend, "primer3");
    assert_eq!(report.executable, fake_primer3);
    assert!(report.reachable);
    assert!(report.version_probe_ok);
    assert_eq!(
        report.version.as_deref(),
        Some("primer3_core synthetic-fixture 2.6.1")
    );
    assert!(report.error.is_none());
}

#[test]
fn test_primer3_preflight_report_missing_executable() {
    let mut engine = GentleEngine::new();
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Primer3;
    let report = engine.primer3_preflight_report(None, Some("/definitely/missing/primer3_core"));
    assert_eq!(report.backend, "primer3");
    assert_eq!(report.executable, "/definitely/missing/primer3_core");
    assert!(!report.reachable);
    assert!(!report.version_probe_ok);
    assert!(report.error.is_some());
}

#[test]
fn test_real_primer3_preflight_is_opt_in() {
    if env::var_os(EXTERNAL_PRIMER_BINARY_TEST_ENV).is_none() {
        return;
    }
    let mut engine = GentleEngine::new();
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Primer3;
    let report = engine.primer3_preflight_report(None, None);
    assert!(
        report.reachable,
        "expected a reachable primer3 executable when {} is set, got {:?}",
        EXTERNAL_PRIMER_BINARY_TEST_ENV, report
    );
}

#[test]
fn test_real_primer3_design_primer_pairs_is_opt_in() {
    if env::var_os(EXTERNAL_PRIMER_BINARY_TEST_ENV).is_none() {
        return;
    }
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        seq(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        ),
    );
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Primer3;
    let result = engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 24,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                non_annealing_5prime_tail: None,
                fixed_5prime: None,
                fixed_3prime: None,
                required_motifs: vec![],
                forbidden_motifs: vec![],
                locked_positions: vec![],
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 24,
                location_0based: Some(90),
                start_0based: None,
                end_0based: None,
                min_tm_c: 0.0,
                max_tm_c: 100.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                non_annealing_5prime_tail: None,
                fixed_5prime: None,
                fixed_3prime: None,
                required_motifs: vec![],
                forbidden_motifs: vec![],
                locked_positions: vec![],
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(5),
            report_id: Some("real_primer3_opt_in".to_string()),
        })
        .expect("real primer3 design should complete when opt-in env is set");
    let report = engine
        .get_primer_design_report("real_primer3_opt_in")
        .expect("real primer3 report");
    assert_eq!(report.backend.requested, "primer3");
    assert_eq!(report.backend.used, "primer3");
    assert!(!report.pairs.is_empty());
    assert!(!result.created_seq_ids.is_empty());
}

#[test]
fn test_design_primer_pairs_rejects_invalid_roi() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        seq("ATGCGTACGATCGTAGCTAGCTAGCTAGCATCGATCG"),
    );
    let mut engine = GentleEngine::from_state(state);
    let err = engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 10,
            roi_end_0based: 200,
            forward: PrimerDesignSideConstraint::default(),
            reverse: PrimerDesignSideConstraint::default(),
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 120,
            max_tm_delta_c: Some(2.0),
            max_pairs: Some(100),
            report_id: Some("bad_roi".to_string()),
        })
        .expect_err("invalid roi");
    assert!(err.message.contains("ROI"));
}

#[test]
fn test_design_qpcr_assays_persists_report() {
    let template_seq = "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG";
    let mut state = ProjectState::default();
    state.sequences.insert("tpl".to_string(), seq(template_seq));
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;
    let result = engine
        .apply(Operation::DesignQpcrAssays {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                non_annealing_5prime_tail: None,
                fixed_5prime: None,
                fixed_3prime: None,
                required_motifs: vec![],
                forbidden_motifs: vec![],
                locked_positions: vec![],
            },
            reverse: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(90),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                non_annealing_5prime_tail: None,
                fixed_5prime: None,
                fixed_3prime: None,
                required_motifs: vec![],
                forbidden_motifs: vec![],
                locked_positions: vec![],
            },
            probe: PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(50),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 100,
                non_annealing_5prime_tail: None,
                fixed_5prime: None,
                fixed_3prime: None,
                required_motifs: vec![],
                forbidden_motifs: vec![],
                locked_positions: vec![],
            },
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_probe_tm_delta_c: Some(100.0),
            max_assays: Some(10),
            report_id: Some("tp73_qpcr".to_string()),
        })
        .expect("design qpcr assays");
    assert!(
        result
            .messages
            .iter()
            .any(|line| line.contains("qPCR-design report"))
    );
    let report = engine
        .get_qpcr_design_report("tp73_qpcr")
        .expect("qpcr report by id");
    assert_eq!(report.report_id, "tp73_qpcr");
    assert_eq!(report.template, "tpl");
    assert!(!report.assays.is_empty());
    assert_eq!(report.backend.requested, "internal");
    assert_eq!(report.backend.used, "internal");
    let listed = engine.list_qpcr_design_reports();
    assert!(listed.iter().any(|row| row.report_id == "tp73_qpcr"));
}

#[test]
fn test_design_primer_pairs_enforces_side_sequence_constraints() {
    let template_seq = "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG";
    let mut state = ProjectState::default();
    state.sequences.insert("tpl".to_string(), seq(template_seq));
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;

    let base_forward = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(5),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };
    let base_reverse = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(90),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };

    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: base_forward.clone(),
            reverse: base_reverse.clone(),
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("baseline_side_constraints".to_string()),
        })
        .expect("baseline primer design");
    let baseline = engine
        .get_primer_design_report("baseline_side_constraints")
        .expect("baseline report");
    assert_eq!(baseline.pair_count, 1);
    let pair = baseline.pairs[0].clone();

    let mut constrained_forward = base_forward.clone();
    constrained_forward.fixed_5prime = Some(pair.forward.sequence[..4].to_string());
    constrained_forward.fixed_3prime =
        Some(pair.forward.sequence[pair.forward.sequence.len() - 4..].to_string());
    constrained_forward.required_motifs = vec![pair.forward.sequence[5..9].to_string()];
    constrained_forward.locked_positions = vec![PrimerDesignBaseLock {
        offset_0based: 3,
        base: pair.forward.sequence[3..4].to_string(),
    }];
    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: constrained_forward.clone(),
            reverse: base_reverse.clone(),
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("matching_side_constraints".to_string()),
        })
        .expect("matching side constraints");
    let matching = engine
        .get_primer_design_report("matching_side_constraints")
        .expect("matching report");
    assert_eq!(matching.pair_count, 1);

    let mut failing_forward = constrained_forward;
    let first = pair.forward.sequence.as_bytes()[0] as char;
    let mismatch = ['A', 'C', 'G', 'T']
        .into_iter()
        .find(|base| *base != first)
        .expect("mismatch base");
    failing_forward.fixed_5prime = Some(mismatch.to_string());
    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: failing_forward,
            reverse: base_reverse,
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("failing_side_constraints".to_string()),
        })
        .expect("failing side constraints design");
    let failing = engine
        .get_primer_design_report("failing_side_constraints")
        .expect("failing report");
    assert_eq!(failing.pair_count, 0);
    assert!(failing.rejection_summary.primer_constraint_failure > 0);
}

#[test]
fn test_non_annealing_5prime_tail_is_included_but_tm_gc_use_anneal_segment() {
    let template_seq = "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG";
    let mut state = ProjectState::default();
    state.sequences.insert("tpl".to_string(), seq(template_seq));
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;

    let mut forward = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(5),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };
    let mut reverse = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(90),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };

    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: forward.clone(),
            reverse: reverse.clone(),
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("no_tail".to_string()),
        })
        .expect("no-tail primer design");
    let no_tail = engine
        .get_primer_design_report("no_tail")
        .expect("no-tail report");
    assert_eq!(no_tail.pair_count, 1);
    let no_tail_pair = &no_tail.pairs[0];

    forward.non_annealing_5prime_tail = Some("GAATTC".to_string());
    reverse.non_annealing_5prime_tail = Some("CTCGAG".to_string());
    let with_tail_result = engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward,
            reverse,
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("with_tail".to_string()),
        })
        .expect("tail primer design");
    let with_tail = engine
        .get_primer_design_report("with_tail")
        .expect("tail report");
    assert_eq!(with_tail.pair_count, 1);
    let with_tail_pair = &with_tail.pairs[0];

    assert!(with_tail_pair.forward.sequence.starts_with("GAATTC"));
    assert!(with_tail_pair.reverse.sequence.starts_with("CTCGAG"));
    assert_eq!(with_tail_pair.forward.anneal_length_bp, 20);
    assert_eq!(with_tail_pair.reverse.anneal_length_bp, 20);
    assert_eq!(with_tail_pair.forward.non_annealing_5prime_tail_bp, 6);
    assert_eq!(with_tail_pair.reverse.non_annealing_5prime_tail_bp, 6);
    assert_eq!(
        with_tail_pair.forward.start_0based,
        no_tail_pair.forward.start_0based
    );
    assert_eq!(
        with_tail_pair.reverse.start_0based,
        no_tail_pair.reverse.start_0based
    );
    assert!(
        (with_tail_pair.forward.tm_c - no_tail_pair.forward.tm_c).abs() < f64::EPSILON,
        "forward anneal Tm should stay unchanged when adding non-annealing tail"
    );
    assert!(
        (with_tail_pair.reverse.tm_c - no_tail_pair.reverse.tm_c).abs() < f64::EPSILON,
        "reverse anneal Tm should stay unchanged when adding non-annealing tail"
    );
    assert!(
        (with_tail_pair.forward.gc_fraction - no_tail_pair.forward.gc_fraction).abs()
            < f64::EPSILON,
        "forward anneal GC should stay unchanged when adding non-annealing tail"
    );
    assert!(
        (with_tail_pair.reverse.gc_fraction - no_tail_pair.reverse.gc_fraction).abs()
            < f64::EPSILON,
        "reverse anneal GC should stay unchanged when adding non-annealing tail"
    );
    let amplicon_seq_id = with_tail_result
        .created_seq_ids
        .iter()
        .find(|seq_id| seq_id.ends_with("_amplicon"))
        .expect("tail run should materialize one predicted amplicon")
        .clone();
    let amplicon = engine
        .state()
        .sequences
        .get(&amplicon_seq_id)
        .expect("predicted amplicon sequence should exist");
    let amplicon_seq = amplicon.get_forward_string().to_ascii_uppercase();
    assert!(amplicon_seq.starts_with("GAATTC"));
    let reverse_tail_rc = GentleEngine::reverse_complement("CTCGAG");
    assert!(amplicon_seq.ends_with(&reverse_tail_rc));
}

#[test]
fn test_design_insertion_primer_pairs_records_anchor_shift_compensation() {
    let template_seq = "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG";
    let mut state = ProjectState::default();
    state.sequences.insert("tpl".to_string(), seq(template_seq));
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;

    let base_forward = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(5),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };
    let base_reverse = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(90),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };

    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: base_forward.clone(),
            reverse: base_reverse.clone(),
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("insert_base".to_string()),
        })
        .expect("baseline primer design");
    let base_report = engine
        .get_primer_design_report("insert_base")
        .expect("baseline report");
    assert_eq!(base_report.pair_count, 1);
    let base_pair = &base_report.pairs[0];

    let requested_forward_end = base_pair.forward.end_0based_exclusive + 3;
    let requested_reverse_start = base_pair.reverse.start_0based.saturating_sub(4);
    let insertion = PrimerInsertionIntent {
        requested_forward_3prime_end_0based_exclusive: requested_forward_end,
        requested_reverse_3prime_start_0based: requested_reverse_start,
        forward_extension_5prime: "GAATTC".to_string(),
        reverse_extension_5prime: "CTCGAG".to_string(),
        forward_window_start_0based: 0,
        forward_window_end_0based_exclusive: 70,
        reverse_window_start_0based: 70,
        reverse_window_end_0based_exclusive: template_seq.len(),
        max_anchor_shift_bp: Some(12),
    };
    engine
        .apply(Operation::DesignInsertionPrimerPairs {
            template: "tpl".to_string(),
            insertion,
            forward: base_forward,
            reverse: base_reverse,
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("insert_mode".to_string()),
        })
        .expect("insertion primer design");

    let report = engine
        .get_primer_design_report("insert_mode")
        .expect("insertion report");
    let insertion_context = report
        .insertion_context
        .as_ref()
        .expect("insertion context should be present");
    assert_eq!(
        insertion_context.requested_forward_3prime_end_0based_exclusive,
        requested_forward_end
    );
    assert_eq!(
        insertion_context.requested_reverse_3prime_start_0based,
        requested_reverse_start
    );
    assert_eq!(insertion_context.max_anchor_shift_bp, 12);
    assert_eq!(insertion_context.uncompensable_pair_count, 0);
    assert_eq!(insertion_context.out_of_shift_budget_pair_count, 0);
    let row = &insertion_context.pairs[0];
    assert_eq!(row.forward_anchor_shift_bp, 3);
    assert_eq!(row.reverse_anchor_shift_bp, 4);
    assert!(row.within_shift_budget);
    assert!(row.compensable);
    assert_eq!(
        row.forward_compensation_5prime,
        template_seq[base_pair.forward.end_0based_exclusive..requested_forward_end].to_string()
    );
    assert_eq!(
        row.reverse_compensation_5prime,
        GentleEngine::reverse_complement(
            &template_seq[requested_reverse_start..base_pair.reverse.start_0based]
        )
    );
    assert!(row.compensated_forward_5prime_tail.starts_with("GAATTC"));
    assert!(row.compensated_reverse_5prime_tail.starts_with("CTCGAG"));
}

#[test]
fn test_pcr_overlap_extension_mutagenesis_insertion_materializes_staged_products() {
    let template_seq = "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG";
    let mut state = ProjectState::default();
    state.sequences.insert("tpl".to_string(), seq(template_seq));
    let mut engine = GentleEngine::from_state(state);

    let fixed_side = |location_0based: usize, length_bp: usize| PrimerDesignSideConstraint {
        min_length: length_bp,
        max_length: length_bp,
        location_0based: Some(location_0based),
        start_0based: None,
        end_0based: None,
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        non_annealing_5prime_tail: None,
        fixed_5prime: None,
        fixed_3prime: None,
        required_motifs: vec![],
        forbidden_motifs: vec![],
        locked_positions: vec![],
    };

    let result = engine
        .apply(Operation::PcrOverlapExtensionMutagenesis {
            template: "tpl".to_string(),
            edit_start_0based: 70,
            edit_end_0based_exclusive: 70,
            insert_sequence: "GGTACC".to_string(),
            constraints: OverlapExtensionMutagenesisConstraints {
                overlap_bp: 16,
                outer_forward: fixed_side(8, 20),
                outer_reverse: fixed_side(105, 20),
                inner_forward: fixed_side(88, 16),
                inner_reverse: fixed_side(45, 20),
            },
            output_prefix: Some("oe_insert".to_string()),
        })
        .expect("overlap-extension insertion mutagenesis");

    let find_created = |token: &str| {
        result
            .created_seq_ids
            .iter()
            .find(|seq_id| seq_id.contains(token))
            .cloned()
            .expect("expected created artifact")
    };
    let inner_forward_id = find_created("_inner_fwd");
    let stage1_left_id = find_created("_stage1_left");
    let stage1_right_id = find_created("_stage1_right");
    let mutant_id = find_created("_mutant");

    let inner_forward_seq = engine
        .state()
        .sequences
        .get(&inner_forward_id)
        .expect("inner forward sequence exists")
        .get_forward_string()
        .to_ascii_uppercase();
    let stage1_left_seq = engine
        .state()
        .sequences
        .get(&stage1_left_id)
        .expect("stage1 left sequence exists")
        .get_forward_string()
        .to_ascii_uppercase();
    let stage1_right_seq = engine
        .state()
        .sequences
        .get(&stage1_right_id)
        .expect("stage1 right sequence exists")
        .get_forward_string()
        .to_ascii_uppercase();
    let mutant_seq = engine
        .state()
        .sequences
        .get(&mutant_id)
        .expect("mutant sequence exists")
        .get_forward_string()
        .to_ascii_uppercase();

    let mutated_template = format!("{}{}{}", &template_seq[..70], "GGTACC", &template_seq[70..]);
    let expected_overlap = &mutated_template[65..94];
    let expected_mutant = &mutated_template[8..131];
    assert!(stage1_left_seq.ends_with(expected_overlap));
    assert!(stage1_right_seq.starts_with(expected_overlap));
    assert!(inner_forward_seq.contains("GGTACC"));
    assert_eq!(mutant_seq, expected_mutant.to_string());
    assert!(
        result
            .messages
            .iter()
            .any(|line| line.contains("Overlap-extension insertion mutagenesis"))
    );
    let preview = result
        .protocol_cartoon_preview
        .as_ref()
        .expect("insertion mutagenesis should include cartoon preview");
    assert_eq!(preview.protocol, "pcr.oe.substitution");
    assert_eq!(preview.insert_bp, 6);
    assert_eq!(preview.overlap_bp, expected_overlap.len());
    assert_eq!(preview.flank_bp, 64);
    assert_eq!(
        preview.bindings.template_id.as_deref(),
        Some("pcr.oe.substitution")
    );
    assert!(!preview.bindings.feature_overrides.is_empty());
}

#[test]
fn test_pcr_overlap_extension_mutagenesis_deletion_materializes_staged_products() {
    let template_seq = "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG";
    let mut state = ProjectState::default();
    state.sequences.insert("tpl".to_string(), seq(template_seq));
    let mut engine = GentleEngine::from_state(state);

    let fixed_side = |location_0based: usize, length_bp: usize| PrimerDesignSideConstraint {
        min_length: length_bp,
        max_length: length_bp,
        location_0based: Some(location_0based),
        start_0based: None,
        end_0based: None,
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        non_annealing_5prime_tail: None,
        fixed_5prime: None,
        fixed_3prime: None,
        required_motifs: vec![],
        forbidden_motifs: vec![],
        locked_positions: vec![],
    };

    let result = engine
        .apply(Operation::PcrOverlapExtensionMutagenesis {
            template: "tpl".to_string(),
            edit_start_0based: 70,
            edit_end_0based_exclusive: 80,
            insert_sequence: String::new(),
            constraints: OverlapExtensionMutagenesisConstraints {
                overlap_bp: 12,
                outer_forward: fixed_side(8, 20),
                outer_reverse: fixed_side(105, 20),
                inner_forward: fixed_side(88, 16),
                inner_reverse: fixed_side(45, 20),
            },
            output_prefix: Some("oe_delete".to_string()),
        })
        .expect("overlap-extension deletion mutagenesis");

    let find_created = |token: &str| {
        result
            .created_seq_ids
            .iter()
            .find(|seq_id| seq_id.contains(token))
            .cloned()
            .expect("expected created artifact")
    };
    let stage1_left_id = find_created("_stage1_left");
    let stage1_right_id = find_created("_stage1_right");
    let mutant_id = find_created("_mutant");

    let stage1_left_seq = engine
        .state()
        .sequences
        .get(&stage1_left_id)
        .expect("stage1 left sequence exists")
        .get_forward_string()
        .to_ascii_uppercase();
    let stage1_right_seq = engine
        .state()
        .sequences
        .get(&stage1_right_id)
        .expect("stage1 right sequence exists")
        .get_forward_string()
        .to_ascii_uppercase();
    let mutant_seq = engine
        .state()
        .sequences
        .get(&mutant_id)
        .expect("mutant sequence exists")
        .get_forward_string()
        .to_ascii_uppercase();

    let mutated_template = format!("{}{}", &template_seq[..70], &template_seq[80..]);
    let expected_overlap = &mutated_template[65..78];
    let expected_mutant = &mutated_template[8..115];
    assert!(stage1_left_seq.ends_with(expected_overlap));
    assert!(stage1_right_seq.starts_with(expected_overlap));
    assert_eq!(mutant_seq, expected_mutant.to_string());
    assert!(
        result
            .messages
            .iter()
            .any(|line| line.contains("Overlap-extension deletion mutagenesis"))
    );
    assert!(
        result.protocol_cartoon_preview.is_none(),
        "deletion mutagenesis should not attach substitution cartoon preview"
    );
}

#[test]
fn test_design_primer_pairs_enforces_pair_constraints() {
    let template_seq = "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG";
    let mut state = ProjectState::default();
    state.sequences.insert("tpl".to_string(), seq(template_seq));
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;

    let forward = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(5),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };
    let reverse = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(90),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };

    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: forward.clone(),
            reverse: reverse.clone(),
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("baseline_pair_constraints".to_string()),
        })
        .expect("baseline primer design");
    let baseline = engine
        .get_primer_design_report("baseline_pair_constraints")
        .expect("baseline report");
    assert_eq!(baseline.pair_count, 1);
    let pair = baseline.pairs[0].clone();
    let amplicon =
        &template_seq.as_bytes()[pair.amplicon_start_0based..pair.amplicon_end_0based_exclusive];
    let motif = String::from_utf8(amplicon[8..12].to_vec()).expect("motif text");

    let matching_constraints = PrimerDesignPairConstraint {
        require_roi_flanking: true,
        required_amplicon_motifs: vec![motif.clone()],
        forbidden_amplicon_motifs: vec![],
        fixed_amplicon_start_0based: Some(pair.amplicon_start_0based),
        fixed_amplicon_end_0based_exclusive: Some(pair.amplicon_end_0based_exclusive),
    };
    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: forward.clone(),
            reverse: reverse.clone(),
            pair_constraints: matching_constraints,
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("matching_pair_constraints".to_string()),
        })
        .expect("matching pair constraints design");
    let matching = engine
        .get_primer_design_report("matching_pair_constraints")
        .expect("matching report");
    assert_eq!(matching.pair_count, 1);

    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward,
            reverse,
            pair_constraints: PrimerDesignPairConstraint {
                require_roi_flanking: false,
                required_amplicon_motifs: vec![],
                forbidden_amplicon_motifs: vec![motif],
                fixed_amplicon_start_0based: None,
                fixed_amplicon_end_0based_exclusive: None,
            },
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("failing_pair_constraints".to_string()),
        })
        .expect("failing pair constraints design");
    let failing = engine
        .get_primer_design_report("failing_pair_constraints")
        .expect("failing report");
    assert_eq!(failing.pair_count, 0);
    assert!(failing.rejection_summary.pair_constraint_failure > 0);
}

#[test]
fn test_design_qpcr_assays_enforces_probe_sequence_constraints() {
    let template_seq = "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG";
    let mut state = ProjectState::default();
    state.sequences.insert("tpl".to_string(), seq(template_seq));
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;

    let forward = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(5),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };
    let reverse = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(90),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };
    let probe = PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(50),
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };

    engine
        .apply(Operation::DesignQpcrAssays {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: forward.clone(),
            reverse: reverse.clone(),
            probe: probe.clone(),
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_probe_tm_delta_c: Some(100.0),
            max_assays: Some(10),
            report_id: Some("baseline_qpcr_probe_constraints".to_string()),
        })
        .expect("baseline qpcr design");
    let baseline = engine
        .get_qpcr_design_report("baseline_qpcr_probe_constraints")
        .expect("baseline qpcr report");
    assert_eq!(baseline.assay_count, 1);
    let assay = baseline.assays[0].clone();

    let mut matching_probe = probe.clone();
    matching_probe.fixed_5prime = Some(assay.probe.sequence[..4].to_string());
    matching_probe.fixed_3prime =
        Some(assay.probe.sequence[assay.probe.sequence.len() - 4..].to_string());
    matching_probe.required_motifs = vec![assay.probe.sequence[5..9].to_string()];
    matching_probe.locked_positions = vec![PrimerDesignBaseLock {
        offset_0based: 3,
        base: assay.probe.sequence[3..4].to_string(),
    }];
    engine
        .apply(Operation::DesignQpcrAssays {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: forward.clone(),
            reverse: reverse.clone(),
            probe: matching_probe.clone(),
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_probe_tm_delta_c: Some(100.0),
            max_assays: Some(10),
            report_id: Some("matching_qpcr_probe_constraints".to_string()),
        })
        .expect("matching qpcr design");
    let matching = engine
        .get_qpcr_design_report("matching_qpcr_probe_constraints")
        .expect("matching qpcr report");
    assert_eq!(matching.assay_count, 1);

    let mut failing_probe = matching_probe;
    let first = assay.probe.sequence.as_bytes()[0] as char;
    let mismatch = ['A', 'C', 'G', 'T']
        .into_iter()
        .find(|base| *base != first)
        .expect("mismatch base");
    failing_probe.fixed_5prime = Some(mismatch.to_string());
    engine
        .apply(Operation::DesignQpcrAssays {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward,
            reverse,
            probe: failing_probe,
            pair_constraints: PrimerDesignPairConstraint::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_probe_tm_delta_c: Some(100.0),
            max_assays: Some(10),
            report_id: Some("failing_qpcr_probe_constraints".to_string()),
        })
        .expect("failing qpcr design");
    let failing = engine
        .get_qpcr_design_report("failing_qpcr_probe_constraints")
        .expect("failing qpcr report");
    assert_eq!(failing.assay_count, 0);
}

#[test]
fn test_primer_design_defaults_follow_20_to_30bp_window() {
    let side = PrimerDesignSideConstraint::default();
    assert_eq!(side.min_length, 20);
    assert_eq!(side.max_length, 30);
}

#[test]
fn test_primer_heuristics_detect_gc_clamp_and_homopolymer_runs() {
    let metrics = GentleEngine::compute_primer_heuristic_metrics(b"ATTTTTCG");
    assert_eq!(metrics.length_bp, 8);
    assert_eq!(metrics.three_prime_base, b'G');
    assert!(metrics.three_prime_gc_clamp);
    assert_eq!(metrics.longest_homopolymer_run_bp, 5);
}

#[test]
fn test_primer_pair_dimer_metrics_detect_3prime_complementarity() {
    let metrics = GentleEngine::compute_primer_pair_dimer_metrics(b"ACGTTTGGGG", b"CCCCAAAA");
    assert!(metrics.max_complementary_run_bp >= 4);
    assert!(metrics.max_3prime_complementary_run_bp >= 4);
}

#[test]
fn test_primer_pair_scoring_penalizes_dimer_prone_pairs() {
    let forward = PrimerDesignPrimerRecord {
        sequence: "ATGCGTACGCGTACGCGTAC".to_string(),
        start_0based: 10,
        end_0based_exclusive: 31,
        tm_c: 62.0,
        gc_fraction: 0.57,
        anneal_hits: 1,
        ..PrimerDesignPrimerRecord::default()
    };
    let reverse_good = PrimerDesignPrimerRecord {
        sequence: "TATATGCGATATATGCGATC".to_string(),
        start_0based: 70,
        end_0based_exclusive: 91,
        tm_c: 62.0,
        gc_fraction: 0.45,
        anneal_hits: 1,
        ..PrimerDesignPrimerRecord::default()
    };
    let reverse_bad = PrimerDesignPrimerRecord {
        sequence: GentleEngine::reverse_complement(&forward.sequence),
        start_0based: 70,
        end_0based_exclusive: 91,
        tm_c: 62.0,
        gc_fraction: 0.57,
        anneal_hits: 1,
        ..PrimerDesignPrimerRecord::default()
    };

    let good = GentleEngine::build_primer_design_pair_record(
        forward.clone(),
        reverse_good,
        40,
        60,
        80,
        220,
        3.0,
        120,
    )
    .expect("good pair");
    let bad = GentleEngine::build_primer_design_pair_record(
        forward,
        reverse_bad,
        40,
        60,
        80,
        220,
        3.0,
        120,
    )
    .expect("bad pair");
    assert!(good.score > bad.score);
    assert!(good.rule_flags.primer_pair_dimer_risk_low);
    assert!(!bad.rule_flags.primer_pair_dimer_risk_low);
}

#[test]
fn test_lineage_extract_creates_parent_child_edge() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("x".to_string(), seq(&"ATGC".repeat(100)));
    let mut engine = GentleEngine::from_state(state);
    let _ = engine
        .apply(Operation::ExtractRegion {
            input: "x".to_string(),
            from: 2,
            to: 7,
            output_id: Some("part".to_string()),
        })
        .unwrap();

    let lineage = &engine.state().lineage;
    let x_node = lineage.seq_to_node.get("x").unwrap();
    let part_node = lineage.seq_to_node.get("part").unwrap();
    assert!(
        lineage
            .edges
            .iter()
            .any(|e| e.from_node_id == *x_node && e.to_node_id == *part_node)
    );
}

#[test]
fn test_lineage_ligation_has_two_parents() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("a".to_string(), seq(&"ATGC".repeat(40)));
    state
        .sequences
        .insert("b".to_string(), seq(&"TTAA".repeat(40)));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::Ligation {
            inputs: vec!["a".to_string(), "b".to_string()],
            circularize_if_possible: false,
            output_id: None,
            protocol: LigationProtocol::Blunt,
            output_prefix: Some("ab".to_string()),
            unique: None,
        })
        .unwrap();

    let lineage = &engine.state().lineage;
    let a_node = lineage.seq_to_node.get("a").unwrap();
    let b_node = lineage.seq_to_node.get("b").unwrap();
    let ab_node = lineage.seq_to_node.get(&res.created_seq_ids[0]).unwrap();
    assert!(
        lineage
            .edges
            .iter()
            .any(|e| e.from_node_id == *a_node && e.to_node_id == *ab_node)
    );
    assert!(
        lineage
            .edges
            .iter()
            .any(|e| e.from_node_id == *b_node && e.to_node_id == *ab_node)
    );
}

#[test]
fn test_select_candidate_creates_in_silico_node() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("frag".to_string(), seq(&"ATGC".repeat(40)));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::SelectCandidate {
            input: "frag".to_string(),
            criterion: "band_size_range:150-170bp".to_string(),
            output_id: Some("picked".to_string()),
        })
        .unwrap();

    assert_eq!(res.created_seq_ids, vec!["picked".to_string()]);
    assert!(!res.warnings.is_empty());
    let lineage = &engine.state().lineage;
    let picked_node = lineage.seq_to_node.get("picked").unwrap();
    let picked_origin = &lineage.nodes.get(picked_node).unwrap().origin;
    assert!(matches!(picked_origin, SequenceOrigin::InSilicoSelection));
}

#[test]
fn test_reverse_complement_reverse_complement_and_branch() {
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), seq("ATGCCA"));
    let mut engine = GentleEngine::from_state(state);

    let res_rev = engine
        .apply(Operation::Reverse {
            input: "s".to_string(),
            output_id: Some("s_rev".to_string()),
        })
        .unwrap();
    assert_eq!(res_rev.created_seq_ids, vec!["s_rev".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("s_rev")
            .unwrap()
            .get_forward_string(),
        "ACCGTA"
    );

    let res_comp = engine
        .apply(Operation::Complement {
            input: "s".to_string(),
            output_id: Some("s_comp".to_string()),
        })
        .unwrap();
    assert_eq!(res_comp.created_seq_ids, vec!["s_comp".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("s_comp")
            .unwrap()
            .get_forward_string(),
        "TACGGT"
    );

    let res_rc = engine
        .apply(Operation::ReverseComplement {
            input: "s".to_string(),
            output_id: Some("s_rc".to_string()),
        })
        .unwrap();
    assert_eq!(res_rc.created_seq_ids, vec!["s_rc".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("s_rc")
            .unwrap()
            .get_forward_string(),
        "TGGCAT"
    );

    let res_split = engine
        .apply(Operation::Branch {
            input: "s".to_string(),
            output_id: Some("s_branch".to_string()),
        })
        .unwrap();
    assert_eq!(res_split.created_seq_ids, vec!["s_branch".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("s_branch")
            .unwrap()
            .get_forward_string(),
        "ATGCCA"
    );

    let lineage = &engine.state().lineage;
    let s_node = lineage.seq_to_node.get("s").unwrap();
    for derived in ["s_rev", "s_comp", "s_rc", "s_branch"] {
        let dnode = lineage.seq_to_node.get(derived).unwrap();
        assert!(
            lineage
                .edges
                .iter()
                .any(|e| e.from_node_id == *s_node && e.to_node_id == *dnode)
        );
    }
}

#[test]
fn test_default_max_fragments_per_container_is_80000() {
    let state = ProjectState::default();
    assert_eq!(state.parameters.max_fragments_per_container, 80_000);
}

#[test]
fn test_set_parameter_max_fragments_per_container() {
    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::SetParameter {
            name: "max_fragments_per_container".to_string(),
            value: serde_json::json!(1234),
        })
        .unwrap();
    assert!(
        res.messages
            .iter()
            .any(|m| m.contains("max_fragments_per_container"))
    );
    assert_eq!(engine.state().parameters.max_fragments_per_container, 1234);
}

#[test]
fn test_set_parameter_require_verified_genome_anchor_for_extension() {
    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::SetParameter {
            name: "require_verified_genome_anchor_for_extension".to_string(),
            value: serde_json::json!(true),
        })
        .unwrap();
    assert!(
        res.messages
            .iter()
            .any(|m| { m.contains("require_verified_genome_anchor_for_extension") })
    );
    assert!(
        engine
            .state()
            .parameters
            .require_verified_genome_anchor_for_extension
    );
}

#[test]
fn test_set_parameter_primer_backend_controls() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::SetParameter {
            name: "primer_design_backend".to_string(),
            value: serde_json::json!("primer3"),
        })
        .expect("set primer backend");
    assert_eq!(
        engine.state().parameters.primer_design_backend,
        PrimerDesignBackend::Primer3
    );
    engine
        .apply(Operation::SetParameter {
            name: "primer3_executable".to_string(),
            value: serde_json::json!("/opt/primer3/primer3_core"),
        })
        .expect("set primer3 executable");
    assert_eq!(
        engine.state().parameters.primer3_executable,
        "/opt/primer3/primer3_core"
    );
    engine
        .apply(Operation::SetParameter {
            name: "primer3_executable".to_string(),
            value: serde_json::Value::Null,
        })
        .expect("reset primer3 executable");
    assert_eq!(engine.state().parameters.primer3_executable, "primer3_core");
}

#[test]
fn test_set_parameter_feature_details_font_size() {
    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::SetParameter {
            name: "feature_details_font_size".to_string(),
            value: serde_json::json!(9.5),
        })
        .unwrap();
    assert!(
        res.messages
            .iter()
            .any(|m| m.contains("feature_details_font_size"))
    );
    assert!((engine.state().display.feature_details_font_size - 9.5).abs() < f32::EPSILON);
}

#[test]
fn test_set_parameter_linear_external_feature_label_style() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::SetParameter {
            name: "linear_external_feature_label_font_size".to_string(),
            value: serde_json::json!(12.5),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "linear_external_feature_label_background_opacity".to_string(),
            value: serde_json::json!(0.6),
        })
        .unwrap();
    assert!(
        (engine
            .state()
            .display
            .linear_external_feature_label_font_size
            - 12.5)
            .abs()
            < f32::EPSILON
    );
    assert!(
        (engine
            .state()
            .display
            .linear_external_feature_label_background_opacity
            - 0.6)
            .abs()
            < f32::EPSILON
    );
}

#[test]
fn test_set_parameter_regulatory_feature_max_view_span_bp() {
    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::SetParameter {
            name: "regulatory_feature_max_view_span_bp".to_string(),
            value: serde_json::json!(50000),
        })
        .unwrap();
    assert!(
        res.messages
            .iter()
            .any(|m| m.contains("regulatory_feature_max_view_span_bp"))
    );
    assert_eq!(
        engine.state().display.regulatory_feature_max_view_span_bp,
        50_000
    );
}

#[test]
fn test_set_parameter_gc_content_bin_size_bp() {
    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::SetParameter {
            name: "gc_content_bin_size_bp".to_string(),
            value: serde_json::json!(250),
        })
        .unwrap();
    assert!(
        res.messages
            .iter()
            .any(|m| m.contains("gc_content_bin_size_bp"))
    );
    assert_eq!(engine.state().display.gc_content_bin_size_bp, 250);
}

#[test]
fn test_linear_letter_defaults_are_auto_adaptive_and_compressed_enabled() {
    let engine = GentleEngine::new();
    assert_eq!(
        engine.state().display.linear_sequence_letter_layout_mode,
        LinearSequenceLetterLayoutMode::AutoAdaptive
    );
    assert!(
        engine
            .state()
            .display
            .linear_sequence_helical_letters_enabled
    );
}

#[test]
fn test_set_parameter_linear_sequence_base_text_max_view_span_bp() {
    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::SetParameter {
            name: "linear_sequence_base_text_max_view_span_bp".to_string(),
            value: serde_json::json!(1200),
        })
        .unwrap();
    assert!(res.messages.iter().any(|m| m.contains("deprecated no-op")));
    assert_eq!(
        engine
            .state()
            .display
            .linear_sequence_base_text_max_view_span_bp,
        500
    );
}

#[test]
fn test_set_parameter_sequence_panel_max_text_length_bp() {
    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::SetParameter {
            name: "sequence_panel_max_text_length_bp".to_string(),
            value: serde_json::json!(200_000),
        })
        .unwrap();
    assert!(
        res.messages
            .iter()
            .any(|m| m.contains("sequence_panel_max_text_length_bp"))
    );
    assert_eq!(
        engine.state().display.sequence_panel_max_text_length_bp,
        200_000
    );
    engine
        .apply(Operation::SetParameter {
            name: "sequence_panel_max_text_length_bp".to_string(),
            value: serde_json::json!(0),
        })
        .unwrap();
    assert_eq!(engine.state().display.sequence_panel_max_text_length_bp, 0);
}

#[test]
fn test_set_parameter_linear_helical_display_controls() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::SetParameter {
            name: "linear_sequence_helical_letters_enabled".to_string(),
            value: serde_json::json!(true),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "linear_sequence_helical_max_view_span_bp".to_string(),
            value: serde_json::json!(2500),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "linear_sequence_condensed_max_view_span_bp".to_string(),
            value: serde_json::json!(1600),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "linear_sequence_helical_phase_offset_bp".to_string(),
            value: serde_json::json!(4),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "linear_sequence_letter_layout_mode".to_string(),
            value: serde_json::json!("condensed_10_row"),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "linear_hide_backbone_when_sequence_bases_visible".to_string(),
            value: serde_json::json!(true),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "linear_helical_parallel_strands".to_string(),
            value: serde_json::json!(false),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "linear_show_reverse_strand_bases".to_string(),
            value: serde_json::json!(false),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "reverse_strand_visual_opacity".to_string(),
            value: serde_json::json!(0.45),
        })
        .unwrap();
    assert!(
        engine
            .state()
            .display
            .linear_sequence_helical_letters_enabled
    );
    assert_eq!(
        engine
            .state()
            .display
            .linear_sequence_helical_max_view_span_bp,
        2000
    );
    assert_eq!(
        engine
            .state()
            .display
            .linear_sequence_condensed_max_view_span_bp,
        1500
    );
    assert_eq!(
        engine
            .state()
            .display
            .linear_sequence_helical_phase_offset_bp,
        4
    );
    assert_eq!(
        engine.state().display.linear_sequence_letter_layout_mode,
        LinearSequenceLetterLayoutMode::Condensed10Row
    );
    assert!(
        engine
            .state()
            .display
            .linear_hide_backbone_when_sequence_bases_visible
    );
    assert!(!engine.state().display.linear_helical_parallel_strands);
    assert!(!engine.state().display.linear_show_double_strand_bases);
    assert!((engine.state().display.reverse_strand_visual_opacity - 0.45).abs() < f32::EPSILON);
}

#[test]
fn test_set_parameter_linear_layout_mode_accepts_new_aliases() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::SetParameter {
            name: "linear_sequence_letter_layout_mode".to_string(),
            value: serde_json::json!("auto"),
        })
        .unwrap();
    assert_eq!(
        engine.state().display.linear_sequence_letter_layout_mode,
        LinearSequenceLetterLayoutMode::AutoAdaptive
    );
    engine
        .apply(Operation::SetParameter {
            name: "linear_sequence_letter_layout_mode".to_string(),
            value: serde_json::json!("standard"),
        })
        .unwrap();
    assert_eq!(
        engine.state().display.linear_sequence_letter_layout_mode,
        LinearSequenceLetterLayoutMode::StandardLinear
    );
    engine
        .apply(Operation::SetParameter {
            name: "linear_sequence_letter_layout_mode".to_string(),
            value: serde_json::json!("helical"),
        })
        .unwrap();
    assert_eq!(
        engine.state().display.linear_sequence_letter_layout_mode,
        LinearSequenceLetterLayoutMode::ContinuousHelical
    );
}

#[test]
fn test_set_parameter_vcf_display_controls() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::SetParameter {
            name: "vcf_display_show_snp".to_string(),
            value: serde_json::json!(false),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "vcf_display_pass_only".to_string(),
            value: serde_json::json!(true),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "vcf_display_use_min_qual".to_string(),
            value: serde_json::json!(true),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "vcf_display_min_qual".to_string(),
            value: serde_json::json!(42.5),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "vcf_display_required_info_keys".to_string(),
            value: serde_json::json!("ac,ann"),
        })
        .unwrap();

    let display = &engine.state().display;
    assert!(!display.vcf_display_show_snp);
    assert!(display.vcf_display_pass_only);
    assert!(display.vcf_display_use_min_qual);
    assert!((display.vcf_display_min_qual - 42.5).abs() < f64::EPSILON);
    assert_eq!(
        display.vcf_display_required_info_keys,
        vec!["AC".to_string(), "ANN".to_string()]
    );
}

#[test]
fn test_digest_respects_max_fragments_per_container() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("x".to_string(), seq("ATGGATCCGCATGGATCCGCATGGATCCGC"));
    state.parameters.max_fragments_per_container = 2;
    let mut engine = GentleEngine::from_state(state);
    let err = engine
        .apply(Operation::Digest {
            input: "x".to_string(),
            enzymes: vec!["BamHI".to_string()],
            output_prefix: Some("frag".to_string()),
        })
        .unwrap_err();
    assert!(err.message.contains("max_fragments_per_container"));
}

#[test]
fn test_filter_by_molecular_weight_with_error_range() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("a".to_string(), seq(&"A".repeat(100)));
    state
        .sequences
        .insert("b".to_string(), seq(&"A".repeat(180)));
    state
        .sequences
        .insert("c".to_string(), seq(&"A".repeat(260)));
    let mut engine = GentleEngine::from_state(state);

    let res = engine
        .apply(Operation::FilterByMolecularWeight {
            inputs: vec!["a".to_string(), "b".to_string(), "c".to_string()],
            min_bp: 150,
            max_bp: 200,
            error: 0.10,
            unique: false,
            output_prefix: Some("mw".to_string()),
        })
        .unwrap();

    assert_eq!(res.created_seq_ids.len(), 1);
    let out = engine.state().sequences.get("mw_1").unwrap();
    assert_eq!(out.len(), 180);
}

#[test]
fn test_filter_by_molecular_weight_unique_fails_on_multiple_matches() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("a".to_string(), seq(&"A".repeat(100)));
    state
        .sequences
        .insert("b".to_string(), seq(&"A".repeat(105)));
    state
        .sequences
        .insert("c".to_string(), seq(&"A".repeat(200)));
    let mut engine = GentleEngine::from_state(state);

    let err = engine
        .apply(Operation::FilterByMolecularWeight {
            inputs: vec!["a".to_string(), "b".to_string(), "c".to_string()],
            min_bp: 95,
            max_bp: 105,
            error: 0.10,
            unique: true,
            output_prefix: Some("mw".to_string()),
        })
        .unwrap_err();

    assert!(err.message.contains("exactly one match"));
}

#[test]
fn test_filter_by_design_constraints_practical_filters() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("good".to_string(), seq("GACTGACTGACTGACTGACT"));
    state
        .sequences
        .insert("low_gc".to_string(), seq("ATATATATATATATATATAT"));
    state
        .sequences
        .insert("homopoly".to_string(), seq("GACAAAAAGACTGACTGACT"));
    state
        .sequences
        .insert("u6_t4".to_string(), seq("GACTTTTGACTGACTGACT"));
    state
        .sequences
        .insert("amb".to_string(), seq("GACTNNACTGACTGACTGAC"));
    let mut engine = GentleEngine::from_state(state);

    let res = engine
        .apply(Operation::FilterByDesignConstraints {
            inputs: vec![
                "good".to_string(),
                "low_gc".to_string(),
                "homopoly".to_string(),
                "u6_t4".to_string(),
                "amb".to_string(),
            ],
            gc_min: Some(0.30),
            gc_max: Some(0.70),
            max_homopolymer_run: Some(4),
            reject_ambiguous_bases: Some(true),
            avoid_u6_terminator_tttt: Some(true),
            forbidden_motifs: vec![],
            unique: false,
            output_prefix: Some("design".to_string()),
        })
        .unwrap();

    assert_eq!(res.created_seq_ids.len(), 1);
    assert_eq!(res.created_seq_ids[0], "design_1".to_string());
    assert!(engine.state().sequences.contains_key("design_1"));
}

#[test]
fn test_filter_by_design_constraints_unique_fails_on_multiple_matches() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("a".to_string(), seq("GACTGACTGACTGACTGACT"));
    state
        .sequences
        .insert("b".to_string(), seq("GACCGACTGACTGACTGACC"));
    let mut engine = GentleEngine::from_state(state);

    let err = engine
        .apply(Operation::FilterByDesignConstraints {
            inputs: vec!["a".to_string(), "b".to_string()],
            gc_min: Some(0.20),
            gc_max: Some(0.80),
            max_homopolymer_run: Some(6),
            reject_ambiguous_bases: Some(true),
            avoid_u6_terminator_tttt: Some(true),
            forbidden_motifs: vec![],
            unique: true,
            output_prefix: Some("design".to_string()),
        })
        .unwrap_err();

    assert!(err.message.contains("exactly one match"));
}

#[test]
fn test_filter_by_design_constraints_accepts_legacy_operation_name() {
    let json = r#"{
            "FilterBySequenceQuality": {
                "inputs": ["a"],
                "gc_min": 0.30,
                "gc_max": 0.70,
                "max_homopolymer_run": 4,
                "reject_ambiguous_bases": true,
                "avoid_u6_terminator_tttt": true,
                "forbidden_motifs": [],
                "unique": false,
                "output_prefix": "legacy"
            }
        }"#;
    let op: Operation = serde_json::from_str(json).expect("legacy op json parses");
    match op {
        Operation::FilterByDesignConstraints {
            inputs,
            gc_min,
            gc_max,
            max_homopolymer_run,
            unique,
            output_prefix,
            ..
        } => {
            assert_eq!(inputs, vec!["a".to_string()]);
            assert_eq!(gc_min, Some(0.30));
            assert_eq!(gc_max, Some(0.70));
            assert_eq!(max_homopolymer_run, Some(4));
            assert!(!unique);
            assert_eq!(output_prefix, Some("legacy".to_string()));
        }
        other => panic!("unexpected operation variant: {:?}", other),
    }
}

#[test]
fn test_pcr_single_amplicon() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::Pcr {
            template: "tpl".to_string(),
            forward_primer: "ATGAAA".to_string(),
            reverse_primer: "AAACCC".to_string(),
            output_id: Some("amp".to_string()),
            unique: Some(true),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids, vec!["amp".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("amp")
            .unwrap()
            .get_forward_string(),
        "ATGAAACCCGGGTTT"
    );
}

#[test]
fn test_pcr_unique_fails_on_multiple_amplicons() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq("AAAACCCCGGGGAAAACCCCGGGG"));
    let mut engine = GentleEngine::from_state(state);
    let err = engine
        .apply(Operation::Pcr {
            template: "tpl".to_string(),
            forward_primer: "AAAA".to_string(),
            reverse_primer: "CCCC".to_string(),
            output_id: None,
            unique: Some(true),
        })
        .unwrap_err();
    assert!(err.message.contains("unique=true"));
}

#[test]
fn test_pcr_multiple_amplicons_without_unique() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq("AAAACCCCGGGGAAAACCCCGGGG"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::Pcr {
            template: "tpl".to_string(),
            forward_primer: "AAAA".to_string(),
            reverse_primer: "CCCC".to_string(),
            output_id: None,
            unique: Some(false),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids.len(), 3);
}

#[test]
fn test_pcr_advanced_inserts_5prime_tail() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::PcrAdvanced {
            template: "tpl".to_string(),
            forward_primer: PcrPrimerSpec {
                sequence: "GGATCCATGAAA".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(0),
                require_3prime_exact_bases: Some(4),
                library_mode: None,
                max_variants: None,
                sample_seed: None,
            },
            reverse_primer: PcrPrimerSpec {
                sequence: "AAACCC".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(0),
                require_3prime_exact_bases: Some(4),
                library_mode: None,
                max_variants: None,
                sample_seed: None,
            },
            output_id: Some("amp_adv".to_string()),
            unique: Some(true),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids, vec!["amp_adv".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("amp_adv")
            .unwrap()
            .get_forward_string(),
        "GGATCCATGAAACCCGGGTTT"
    );
}

#[test]
fn test_pcr_advanced_allows_partial_match_and_introduces_mutation() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::PcrAdvanced {
            template: "tpl".to_string(),
            forward_primer: PcrPrimerSpec {
                sequence: "ATCAAA".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(1),
                require_3prime_exact_bases: Some(3),
                library_mode: None,
                max_variants: None,
                sample_seed: None,
            },
            reverse_primer: PcrPrimerSpec {
                sequence: "AAACCC".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(0),
                require_3prime_exact_bases: Some(4),
                library_mode: None,
                max_variants: None,
                sample_seed: None,
            },
            output_id: Some("amp_mut".to_string()),
            unique: Some(true),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids, vec!["amp_mut".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("amp_mut")
            .unwrap()
            .get_forward_string(),
        "ATCAAACCCGGGTTT"
    );
}

#[test]
fn test_pcr_mutagenesis_single_snp_success() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::PcrMutagenesis {
            template: "tpl".to_string(),
            forward_primer: PcrPrimerSpec {
                sequence: "ATCAAA".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(1),
                require_3prime_exact_bases: Some(3),
                library_mode: None,
                max_variants: None,
                sample_seed: None,
            },
            reverse_primer: PcrPrimerSpec {
                sequence: "AAACCC".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(0),
                require_3prime_exact_bases: Some(4),
                library_mode: None,
                max_variants: None,
                sample_seed: None,
            },
            mutations: vec![SnpMutationSpec {
                zero_based_position: 2,
                reference: "G".to_string(),
                alternate: "C".to_string(),
            }],
            output_id: Some("mut1".to_string()),
            unique: Some(true),
            require_all_mutations: Some(true),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids, vec!["mut1".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("mut1")
            .unwrap()
            .get_forward_string(),
        "ATCAAACCCGGGTTT"
    );
}

#[test]
fn test_pcr_mutagenesis_fails_when_requested_snp_not_introduced() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
    let mut engine = GentleEngine::from_state(state);
    let err = engine
        .apply(Operation::PcrMutagenesis {
            template: "tpl".to_string(),
            forward_primer: PcrPrimerSpec {
                sequence: "ATCAAA".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(1),
                require_3prime_exact_bases: Some(3),
                library_mode: None,
                max_variants: None,
                sample_seed: None,
            },
            reverse_primer: PcrPrimerSpec {
                sequence: "AAACCC".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(0),
                require_3prime_exact_bases: Some(4),
                library_mode: None,
                max_variants: None,
                sample_seed: None,
            },
            mutations: vec![SnpMutationSpec {
                zero_based_position: 2,
                reference: "G".to_string(),
                alternate: "T".to_string(),
            }],
            output_id: Some("mut_fail".to_string()),
            unique: Some(true),
            require_all_mutations: Some(true),
        })
        .unwrap_err();
    assert!(
        err.message
            .contains("No amplicon introduced all requested mutations")
    );
}

#[test]
fn test_pcr_advanced_degenerate_primer_enumerate() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::PcrAdvanced {
            template: "tpl".to_string(),
            forward_primer: PcrPrimerSpec {
                sequence: "ATNAAA".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(1),
                require_3prime_exact_bases: Some(3),
                library_mode: Some(PrimerLibraryMode::Enumerate),
                max_variants: Some(4),
                sample_seed: None,
            },
            reverse_primer: PcrPrimerSpec {
                sequence: "AAACCC".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(0),
                require_3prime_exact_bases: Some(4),
                library_mode: None,
                max_variants: None,
                sample_seed: None,
            },
            output_id: None,
            unique: Some(false),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids.len(), 4);
}

#[test]
fn test_pcr_advanced_degenerate_primer_sample_mode() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("tpl".to_string(), seq("ATGAAACCCGGGTTT"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::PcrAdvanced {
            template: "tpl".to_string(),
            forward_primer: PcrPrimerSpec {
                sequence: "ATNAAA".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(1),
                require_3prime_exact_bases: Some(3),
                library_mode: Some(PrimerLibraryMode::Sample),
                max_variants: Some(2),
                sample_seed: Some(42),
            },
            reverse_primer: PcrPrimerSpec {
                sequence: "AAACCC".to_string(),
                anneal_len: Some(6),
                max_mismatches: Some(0),
                require_3prime_exact_bases: Some(4),
                library_mode: None,
                max_variants: None,
                sample_seed: None,
            },
            output_id: None,
            unique: Some(false),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids.len(), 2);
}

#[test]
fn test_load_file_operation() {
    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::LoadFile {
            path: "test_files/pGEX_3X.fa".to_string(),
            as_id: Some("pgex".to_string()),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids, vec!["pgex".to_string()]);
    assert!(engine.state().sequences.contains_key("pgex"));
}

#[test]
fn test_fetch_genbank_accession_operation_loads_sequence_and_anchor() {
    let _guard = crate::genomes::genbank_env_lock().lock().unwrap();
    let td = tempdir().unwrap();
    let mock_dir = td.path().join("mock");
    fs::create_dir_all(&mock_dir).unwrap();
    fs::copy(
        "test_files/tp73.ncbi.gb",
        mock_dir.join("NC_000001.gbwithparts"),
    )
    .unwrap();
    let efetch_template = format!("file://{}/{{accession}}.{{rettype}}", mock_dir.display());
    let _efetch_env = EnvVarGuard::set("GENTLE_NCBI_EFETCH_URL", &efetch_template);

    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::FetchGenBankAccession {
            accession: "NC_000001".to_string(),
            as_id: Some("tp73_fetch".to_string()),
        })
        .expect("fetch genbank accession");
    assert_eq!(res.created_seq_ids, vec!["tp73_fetch".to_string()]);
    assert!(res.messages.iter().any(|m| {
        m.contains("Fetched GenBank accession 'NC_000001'") && m.contains("tp73_fetch")
    }));
    assert!(
        engine
            .list_sequences_with_genome_anchor()
            .iter()
            .any(|seq_id| seq_id == "tp73_fetch")
    );
    let provenance = engine
        .state()
        .metadata
        .get(PROVENANCE_METADATA_KEY)
        .and_then(|v| v.get(GENOME_EXTRACTIONS_METADATA_KEY))
        .and_then(|v| v.as_array())
        .and_then(|records| {
            records.iter().find(|entry| {
                entry
                    .get("seq_id")
                    .and_then(|v| v.as_str())
                    .map(|id| id == "tp73_fetch")
                    .unwrap_or(false)
            })
        })
        .cloned()
        .expect("provenance record");
    assert_eq!(
        provenance
            .get("sequence_source_type")
            .and_then(|v| v.as_str())
            .unwrap_or_default(),
        "genbank_accession"
    );
}

#[test]
fn test_fetch_dbsnp_region_operation_extracts_annotated_slice_and_provenance() {
    let _guard = crate::genomes::genbank_env_lock().lock().unwrap();
    let td = tempdir().unwrap();
    let fasta_path = td.path().join("toy.fa");
    let ann_path = td.path().join("toy.gtf");
    fs::write(
        &fasta_path,
        ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
    )
    .unwrap();
    fs::write(
        &ann_path,
        "chr1\tsrc\tgene\t1\t60\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"ONE\";\nchr1\tsrc\ttranscript\t1\t60\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"; gene_name \"ONE\";\nchr1\tsrc\texon\t1\t60\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"; exon_number \"1\";\nchr1\tsrc\tCDS\t10\t45\t.\t+\t0\tgene_id \"GENE1\"; transcript_id \"TX1\";\n",
    )
    .unwrap();
    let cache_dir = td.path().join("cache");
    let cache_dir_str = cache_dir.display().to_string();
    let catalog_path = td.path().join("catalog.json");
    fs::write(
        &catalog_path,
        format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy dbsnp genome",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            fasta_path.display(),
            ann_path.display(),
            cache_dir.display()
        ),
    )
    .unwrap();
    let catalog_path_str = catalog_path.display().to_string();
    let mock_dir = td.path().join("mock_dbsnp");
    fs::create_dir_all(&mock_dir).unwrap();
    fs::write(
        mock_dir.join("123.json"),
        r#"{
  "refsnp_id": "123",
  "primary_snapshot_data": {
    "placements_with_allele": [
      {
        "seq_id": "NC_000001.11",
        "is_ptlp": true,
        "placement_annot": {
          "seq_id_traits_by_assembly": [
            {
              "assembly_name": "ToyGenome.1",
              "is_top_level": true,
              "is_alt": false,
              "is_patch": false,
              "is_chromosome": true
            }
          ]
        },
        "alleles": [
          {
            "allele": {
              "spdi": {
                "seq_id": "NC_000001.11",
                "position": 29,
                "deleted_sequence": "A",
                "inserted_sequence": "G"
              }
            }
          }
        ]
      }
    ],
    "allele_annotations": [
      {
        "assembly_annotation": [
          {
            "genes": [
              {
                "locus": "tagA"
              }
            ]
          }
        ]
      }
    ]
  }
}
"#,
    )
    .unwrap();
    let refsnp_template = format!("file://{}/{{refsnp_id}}.json", mock_dir.display());
    let _dbsnp_env = EnvVarGuard::set("GENTLE_NCBI_DBSNP_REFSNP_URL", &refsnp_template);

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: Some(cache_dir_str.clone()),
            timeout_seconds: None,
        })
        .expect("prepare ToyGenome");

    let result = engine
        .apply(Operation::FetchDbSnpRegion {
            rs_id: "rs123".to_string(),
            genome_id: "ToyGenome".to_string(),
            flank_bp: Some(20),
            output_id: Some("rs123_local".to_string()),
            annotation_scope: Some(GenomeAnnotationScope::Full),
            max_annotation_features: Some(0),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: Some(cache_dir_str.clone()),
        })
        .expect("fetch dbsnp region");

    assert_eq!(result.created_seq_ids, vec!["rs123_local".to_string()]);
    assert!(
        result.messages.iter().any(|message| {
            message.contains("Resolved dbSNP 'rs123'") && message.contains("1 [NC_000001.11]:30")
        }),
        "messages were: {:?}",
        result.messages
    );
    let telemetry = result
        .genome_annotation_projection
        .as_ref()
        .expect("annotation telemetry");
    assert_eq!(telemetry.requested_scope, "full");
    assert!(telemetry.attached_feature_count > 0);
    let extracted = engine
        .state()
        .sequences
        .get("rs123_local")
        .expect("extracted sequence");
    assert_eq!(extracted.len(), 41);
    assert!(!extracted.features().is_empty());
    let marker = extracted
        .features()
        .iter()
        .find(|feature| {
            feature.kind.to_string().eq_ignore_ascii_case("variation")
                && feature
                    .qualifier_values("label")
                    .any(|value| value == "rs123")
        })
        .expect("dbSNP marker feature");
    assert_eq!(
        marker.location.find_bounds().expect("marker bounds"),
        (20, 21)
    );
    assert_eq!(
        marker.qualifier_values("db_xref").next(),
        Some("dbSNP:rs123")
    );
    assert_eq!(
        marker.qualifier_values("genomic_position_1based").next(),
        Some("30")
    );
    assert_eq!(
        marker.qualifier_values("refseq_accession").next(),
        Some("NC_000001.11")
    );
    assert_eq!(
        marker.qualifier_values("gentle_generated").next(),
        Some("dbsnp_variant_marker")
    );

    let provenance = engine
        .state()
        .metadata
        .get(PROVENANCE_METADATA_KEY)
        .and_then(|v| v.get(GENOME_EXTRACTIONS_METADATA_KEY))
        .and_then(|v| v.as_array())
        .and_then(|records| {
            records.iter().find(|entry| {
                entry
                    .get("seq_id")
                    .and_then(|v| v.as_str())
                    .map(|id| id == "rs123_local")
                    .unwrap_or(false)
            })
        })
        .cloned()
        .expect("dbsnp provenance record");
    assert_eq!(
        provenance
            .get("operation")
            .and_then(|v| v.as_str())
            .unwrap_or_default(),
        "FetchDbSnpRegion"
    );
    assert_eq!(
        provenance
            .get("chromosome")
            .and_then(|v| v.as_str())
            .unwrap_or_default(),
        "NC_000001.11"
    );
}

#[test]
fn test_fetch_dbsnp_region_operation_emits_staged_progress_updates() {
    let td = tempdir().expect("tempdir");
    let cache_dir = td.path().join("cache");
    fs::create_dir_all(&cache_dir).expect("cache dir");
    let cache_dir_str = cache_dir.to_string_lossy().to_string();
    let catalog_path = td.path().join("catalog.json");
    let catalog_path_str = catalog_path.to_string_lossy().to_string();
    fs::write(
        &catalog_path,
        format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy dbsnp genome",
    "sequence_local": "{root}/toy.fa",
    "annotations_local": "{root}/toy.gtf",
    "cache_dir": "{root}/cache"
  }}
}}"#,
            root = td.path().display()
        ),
    )
    .expect("catalog");
    fs::write(
        td.path().join("toy.fa"),
        ">1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
    )
    .expect("toy fasta");
    fs::write(
        td.path().join("toy.gtf"),
        "1\ttoy\tgene\t10\t60\t.\t+\t.\tgene_id \"tagA\"; gene_name \"tagA\";\n",
    )
    .expect("toy gtf");
    let mock_dir = td.path().join("mock_dbsnp");
    fs::create_dir_all(&mock_dir).expect("mock dir");
    fs::write(
        mock_dir.join("123.json"),
        r#"{
  "refsnp_id": "123",
  "primary_snapshot_data": {
    "placements_with_allele": [
      {
        "seq_id": "NC_000001.11",
        "is_ptlp": true,
        "placement_annot": {
          "seq_id_traits_by_assembly": [
            {
              "assembly_name": "ToyGenome.1",
              "is_top_level": true,
              "is_alt": false,
              "is_patch": false,
              "is_chromosome": true
            }
          ]
        },
        "alleles": [
          {
            "allele": {
              "spdi": {
                "seq_id": "NC_000001.11",
                "position": 29,
                "deleted_sequence": "A",
                "inserted_sequence": "G"
              }
            }
          }
        ]
      }
    ]
  }
}
"#,
    )
    .expect("mock dbsnp json");
    let refsnp_template = format!("file://{}/{{refsnp_id}}.json", mock_dir.display());
    let _dbsnp_env = EnvVarGuard::set("GENTLE_NCBI_DBSNP_REFSNP_URL", &refsnp_template);

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: Some(cache_dir_str.clone()),
            timeout_seconds: None,
        })
        .expect("prepare ToyGenome");

    let mut stages: Vec<String> = vec![];
    engine
        .apply_with_progress(
            Operation::FetchDbSnpRegion {
                rs_id: "rs123".to_string(),
                genome_id: "ToyGenome".to_string(),
                flank_bp: Some(20),
                output_id: Some("rs123_progress".to_string()),
                annotation_scope: Some(GenomeAnnotationScope::Full),
                max_annotation_features: Some(0),
                catalog_path: Some(catalog_path_str),
                cache_dir: Some(cache_dir_str),
            },
            |progress| {
                if let OperationProgress::DbSnpFetch(progress) = progress {
                    stages.push(progress.stage.as_str().to_string());
                }
                true
            },
        )
        .expect("fetch dbsnp region with progress");

    assert!(stages.contains(&"validate_input".to_string()));
    assert!(stages.contains(&"inspect_prepared_genome".to_string()));
    assert!(stages.contains(&"contact_server".to_string()));
    assert!(stages.contains(&"wait_response".to_string()));
    assert!(stages.contains(&"parse_response".to_string()));
    assert!(stages.contains(&"resolve_placement".to_string()));
    assert!(stages.contains(&"extract_region".to_string()));
    assert!(stages.contains(&"attach_variant_marker".to_string()));
}

#[test]
fn test_genome_chromosome_matches_accepts_refseq_accessions_for_sex_and_mito_contigs() {
    assert!(GentleEngine::genome_chromosome_matches(
        "chrX",
        "NC_000023.11"
    ));
    assert!(GentleEngine::genome_chromosome_matches(
        "chrY",
        "NC_000024.10"
    ));
    assert!(GentleEngine::genome_chromosome_matches(
        "chrM",
        "NC_012920.1"
    ));
    assert!(!GentleEngine::genome_chromosome_matches(
        "1",
        "NC_000016.10"
    ));
}

#[test]
fn test_load_file_operation_genbank_region_anchor_enables_bed_import() {
    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::LoadFile {
            path: "test_files/tp73.ncbi.gb".to_string(),
            as_id: Some("tp73".to_string()),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids, vec!["tp73".to_string()]);
    assert!(
        engine
            .list_sequences_with_genome_anchor()
            .iter()
            .any(|seq_id| seq_id == "tp73")
    );
    let provenance_catalog = engine
        .state()
        .metadata
        .get(PROVENANCE_METADATA_KEY)
        .and_then(|v| v.get(GENOME_EXTRACTIONS_METADATA_KEY))
        .and_then(|v| v.as_array())
        .and_then(|records| {
            records.iter().find(|entry| {
                entry
                    .get("seq_id")
                    .and_then(|v| v.as_str())
                    .map(|id| id == "tp73")
                    .unwrap_or(false)
            })
        })
        .and_then(|entry| entry.get("catalog_path"))
        .and_then(|v| v.as_str())
        .unwrap_or_default()
        .to_string();
    assert_eq!(provenance_catalog, DEFAULT_GENOME_CATALOG_PATH);

    let td = tempdir().unwrap();
    let bed_path = td.path().join("tp73_anchor_test.bed");
    std::fs::write(&bed_path, "chr1\t3652515\t3652525\tpeak1\t100\t+\n").unwrap();

    let import_res = engine
        .apply(Operation::ImportGenomeBedTrack {
            seq_id: "tp73".to_string(),
            path: bed_path.display().to_string(),
            track_name: Some("anchor-test".to_string()),
            min_score: None,
            max_score: None,
            clear_existing: Some(true),
        })
        .unwrap();
    assert!(
        import_res
            .changed_seq_ids
            .iter()
            .any(|seq_id| seq_id == "tp73")
    );

    let tp73 = engine
        .state()
        .sequences
        .get("tp73")
        .expect("tp73 should exist");
    assert!(
        tp73.features()
            .iter()
            .any(GentleEngine::is_generated_genome_bed_feature)
    );
}

#[test]
fn test_parse_genbank_accession_region_supports_complement() {
    let td = tempdir().unwrap();
    let gb_path = td.path().join("complement_header.gb");
    let gb_text = "\
LOCUS       TEST000001              11 bp    DNA     linear   CON 01-JAN-2000
DEFINITION  Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly.
ACCESSION   NC_000001 REGION: complement(3652516..3652526)
ORIGIN
        1 atgcgatgcga
//
";
    std::fs::write(&gb_path, gb_text).unwrap();
    let parsed = GentleEngine::parse_genbank_accession_region(gb_path.to_string_lossy().as_ref())
        .expect("region should parse");
    assert_eq!(parsed.0, "NC_000001");
    assert_eq!(parsed.1, 3652516);
    assert_eq!(parsed.2, 3652526);
    assert_eq!(parsed.4, '-');
}

#[test]
fn test_import_genome_bed_track_remaps_for_reverse_anchor() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("rev".to_string(), seq("ACGTACGTACG"));
    state.metadata.insert(
        PROVENANCE_METADATA_KEY.to_string(),
        serde_json::json!({
            GENOME_EXTRACTIONS_METADATA_KEY: [
                {
                    "seq_id": "rev",
                    "recorded_at_unix_ms": 1,
                    "operation": "LoadFileGenBankRegion",
                    "genome_id": "GRCh38.p14",
                    "catalog_path": "synthetic",
                    "cache_dir": null,
                    "chromosome": "1",
                    "start_1based": 100,
                    "end_1based": 110,
                    "gene_query": null,
                    "occurrence": null,
                    "gene_id": null,
                    "gene_name": null,
                    "strand": null,
                    "anchor_strand": "-",
                    "sequence_source_type": "genbank_file",
                    "annotation_source_type": "genbank_file",
                    "sequence_source": "synthetic",
                    "annotation_source": "synthetic",
                    "sequence_sha1": null,
                    "annotation_sha1": null
                }
            ]
        }),
    );
    let mut engine = GentleEngine::from_state(state);

    let td = tempdir().unwrap();
    let bed_path = td.path().join("reverse_anchor.bed");
    std::fs::write(&bed_path, "chr1\t99\t101\tpeak1\t100\t+\n").unwrap();

    let result = engine
        .apply(Operation::ImportGenomeBedTrack {
            seq_id: "rev".to_string(),
            path: bed_path.display().to_string(),
            track_name: Some("rev-track".to_string()),
            min_score: None,
            max_score: None,
            clear_existing: Some(true),
        })
        .unwrap();
    assert!(result.changed_seq_ids.iter().any(|id| id == "rev"));

    let dna = engine
        .state()
        .sequences
        .get("rev")
        .expect("rev sequence should exist");
    let feature = dna
        .features()
        .iter()
        .find(|f| GentleEngine::is_generated_genome_bed_feature(f))
        .expect("generated BED feature should exist");
    assert!(crate::feature_location::feature_is_reverse(feature));
    let ranges = crate::feature_location::feature_ranges_sorted_i64(feature);
    assert_eq!(ranges, vec![(9, 11)]);
    assert_eq!(feature.qualifier_values("bed_strand").next(), Some("+"));
    assert_eq!(feature.qualifier_values("strand").next(), Some("-"));
}

#[test]
fn test_load_file_operation_fasta_synthetic_oligo_metadata() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("oligo.fa");
    std::fs::write(&path, ">oligo1 molecule=ssdna\nATGCATGC\n").unwrap();

    let mut engine = GentleEngine::new();
    let res = engine
        .apply(Operation::LoadFile {
            path: path.display().to_string(),
            as_id: Some("oligo".to_string()),
        })
        .unwrap();
    assert_eq!(res.created_seq_ids, vec!["oligo".to_string()]);

    let dna = engine
        .state()
        .sequences
        .get("oligo")
        .expect("oligo sequence should exist");
    assert_eq!(dna.molecule_type(), Some("ssDNA"));
    assert!(dna.overhang().is_blunt());
}

#[test]
fn test_load_file_operation_fasta_dsdna_with_overhang_metadata() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("sticky.fa");
    std::fs::write(&path, ">sticky molecule=dsdna f5=GATC r5=CTAG\nATGCATGC\n").unwrap();

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::LoadFile {
            path: path.display().to_string(),
            as_id: Some("sticky".to_string()),
        })
        .unwrap();

    let dna = engine
        .state()
        .sequences
        .get("sticky")
        .expect("sticky sequence should exist");
    assert_eq!(dna.molecule_type(), Some("dsDNA"));
    assert_eq!(dna.overhang().forward_5, b"GATC".to_vec());
    assert_eq!(dna.overhang().reverse_5, b"CTAG".to_vec());
}

#[test]
fn test_import_and_project_uniprot_swiss_prot() {
    let dir = tempfile::tempdir().unwrap();
    let swiss_path = dir.path().join("toy_uniprot.txt");
    let swiss_text = r#"ID   TOY1_HUMAN              Reviewed;         30 AA.
AC   PTEST1;
DE   RecName: Full=Toy DNA-binding protein;
GN   Name=TOY1;
OS   Homo sapiens (Human).
DR   Ensembl; TX1; ENSPTOY1; ENSGTOY1.
FT   DOMAIN          2..8
FT                   /note="toy domain"
SQ   SEQUENCE   30 AA;  3333 MW;  0000000000000000 CRC64;
     MEEPQSDPSV EPPLSQETFSDLWKLLPEN
//
"#;
    std::fs::write(&swiss_path, swiss_text).unwrap();

    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(300)).expect("valid DNA");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::simple_range(99, 360),
        qualifiers: vec![
            ("gene".into(), Some("TOY1".to_string())),
            ("transcript_id".into(), Some("TX1".to_string())),
            ("label".into(), Some("TX1".to_string())),
            (
                "cds_ranges_1based".into(),
                Some("100-180,300-360".to_string()),
            ),
        ]
        .into_iter()
        .collect(),
    });
    state.sequences.insert("toy_seq".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let import = engine
        .apply(Operation::ImportUniprotSwissProt {
            path: swiss_path.display().to_string(),
            entry_id: None,
        })
        .expect("import uniprot swiss");
    assert!(
        import
            .messages
            .iter()
            .any(|message| message.contains("Imported UniProt SWISS-PROT entry"))
    );
    let entries = engine.list_uniprot_entries();
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].entry_id, "PTEST1");

    let map = engine
        .apply(Operation::ProjectUniprotToGenome {
            seq_id: "toy_seq".to_string(),
            entry_id: "PTEST1".to_string(),
            projection_id: None,
            transcript_id: None,
        })
        .expect("project uniprot");
    assert!(
        map.messages
            .iter()
            .any(|message| message.contains("Projected UniProt entry"))
    );

    let projection = engine
        .get_uniprot_genome_projection("PTEST1@toy_seq")
        .expect("projection should exist");
    assert_eq!(projection.entry_id, "PTEST1");
    assert_eq!(projection.seq_id, "toy_seq");
    assert!(!projection.transcript_projections.is_empty());
    assert!(
        projection
            .transcript_projections
            .iter()
            .any(|row| !row.feature_projections.is_empty())
    );
}

#[test]
fn test_import_uniprot_entry_sequence_is_currently_unsupported() {
    let dir = tempfile::tempdir().unwrap();
    let swiss_path = dir.path().join("toy_uniprot_use.txt");
    let swiss_text = r#"ID   TOY2_HUMAN              Reviewed;         12 AA.
AC   PUSE1;
DE   RecName: Full=Toy use protein;
GN   Name=TOY2;
OS   Homo sapiens (Human).
FT   DOMAIN          2..6
FT                   /note="toy segment"
SQ   SEQUENCE   12 AA;  1200 MW;  0000000000000000 CRC64;
     MEEPQSDPSVEP
//
"#;
    std::fs::write(&swiss_path, swiss_text).unwrap();

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::ImportUniprotSwissProt {
            path: swiss_path.display().to_string(),
            entry_id: Some("TOY_USE".to_string()),
        })
        .expect("import uniprot swiss");

    let import_sequence_err = engine
        .apply(Operation::ImportUniprotEntrySequence {
            entry_id: "TOY_USE".to_string(),
            output_id: None,
        })
        .expect_err("uniprot sequence import should be disabled");
    assert!(matches!(import_sequence_err.code, ErrorCode::Unsupported));
    assert!(
        import_sequence_err
            .message
            .contains("ImportUniprotEntrySequence is currently disabled"),
        "unexpected error message: {}",
        import_sequence_err.message
    );
    assert!(!engine.state().sequences.contains_key("TOY_USE"));
}

#[test]
fn test_fasta_roundtrip_synthetic_dsdna_blunt() {
    let expected = synth_oligo("molecule=dsdna topology=linear", b"ATGCATGC");
    assert_fasta_roundtrip(expected);
}

#[test]
fn test_fasta_roundtrip_synthetic_ssdna() {
    let expected = synth_oligo("molecule=ssdna topology=linear", b"ATGCATGC");
    assert_fasta_roundtrip(expected);
}

#[test]
fn test_fasta_roundtrip_synthetic_rna() {
    let expected = synth_oligo("molecule=rna topology=linear", b"AUGTT");
    assert_fasta_roundtrip(expected);
}

#[test]
fn test_fasta_roundtrip_synthetic_dsdna_with_overhangs() {
    let expected = synth_oligo(
        "molecule=dsdna f5=GATC r5=CTAG topology=linear",
        b"ATGCATGC",
    );
    assert_fasta_roundtrip(expected);
}

#[cfg(unix)]
#[test]
fn test_inspect_rna_structure_returns_rnapkin_textual_output() {
    let td = tempdir().unwrap();
    let fake_rnapkin = install_fake_rnapkin(td.path());
    let _bin_guard = EnvVarGuard::set("GENTLE_RNAPKIN_BIN", &fake_rnapkin);

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("rna".to_string(), synth_oligo("molecule=rna", b"AUGCAU"));
    let engine = GentleEngine::from_state(state);
    let report = engine.inspect_rna_structure("rna").unwrap();

    assert_eq!(report.tool, "rnapkin");
    assert!(report.stdout.contains("rnapkin textual report"));
    assert!(report.stdout.contains("points:"));
}

#[cfg(unix)]
#[test]
fn test_render_rna_structure_svg_operation() {
    let td = tempdir().unwrap();
    let fake_rnapkin = install_fake_rnapkin(td.path());
    let _bin_guard = EnvVarGuard::set("GENTLE_RNAPKIN_BIN", &fake_rnapkin);

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("rna".to_string(), synth_oligo("molecule=rna", b"AUGCAU"));
    let mut engine = GentleEngine::from_state(state);
    let output = td.path().join("rna.structure.svg");
    let output_text = output.display().to_string();

    let res = engine
        .apply(Operation::RenderRnaStructureSvg {
            seq_id: "rna".to_string(),
            path: output_text.clone(),
        })
        .unwrap();

    assert!(res.messages.iter().any(|m| m.contains("RNA structure SVG")));
    let svg = std::fs::read_to_string(output_text).unwrap();
    assert!(svg.contains("<svg"));
}

#[test]
fn test_render_rna_structure_svg_requires_rna_biotype() {
    let mut state = ProjectState::default();
    state.sequences.insert("dna".to_string(), seq("ATGCATGC"));
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let output = tmp.path().with_extension("svg");
    let err = engine
        .apply(Operation::RenderRnaStructureSvg {
            seq_id: "dna".to_string(),
            path: output.display().to_string(),
        })
        .unwrap_err();
    assert!(matches!(err.code, ErrorCode::InvalidInput));
    assert!(err.message.to_ascii_lowercase().contains("rna"));
}

#[test]
fn test_save_file_operation_genbank() {
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), seq("ATGCCA"));
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("gb");
    let path_text = path.display().to_string();
    let res = engine
        .apply(Operation::SaveFile {
            seq_id: "s".to_string(),
            path: path_text.clone(),
            format: ExportFormat::GenBank,
        })
        .unwrap();
    assert!(res.changed_seq_ids.contains(&"s".to_string()));
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("LOCUS"));
}

#[test]
fn test_set_topology_operation() {
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), seq("ATGCCA"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::SetTopology {
            seq_id: "s".to_string(),
            circular: true,
        })
        .unwrap();
    assert!(res.changed_seq_ids.contains(&"s".to_string()));
    assert!(engine.state().sequences.get("s").unwrap().is_circular());
}

#[test]
fn test_recompute_features_operation() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), seq("ATGAAACCCGGGTTT"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::RecomputeFeatures {
            seq_id: "s".to_string(),
        })
        .unwrap();
    assert!(res.changed_seq_ids.contains(&"s".to_string()));
}

#[test]
fn test_prepare_scoring_matrices_avoid_negative_infinity() {
    let matrix = vec![[10.0, 0.0, 0.0, 0.0], [5.0, 0.0, 0.0, 0.0]];
    let (llr, true_log_odds) = GentleEngine::prepare_scoring_matrices(&matrix);
    assert_eq!(llr.len(), 2);
    for col in llr {
        for v in col {
            assert!(v.is_finite());
        }
    }
    assert_eq!(true_log_odds.len(), 2);
    for col in true_log_odds {
        for v in col {
            assert!(v.is_finite());
        }
    }
}

#[test]
fn test_annotate_tfbs_adds_scored_features() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), seq("TTTACGTAAACGTGGG"));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::AnnotateTfbs {
            seq_id: "s".to_string(),
            motifs: vec!["ACGT".to_string()],
            min_llr_bits: Some(0.0),
            min_llr_quantile: Some(0.0),
            per_tf_thresholds: vec![],
            clear_existing: Some(true),
            max_hits: None,
        })
        .unwrap();
    assert!(res.changed_seq_ids.contains(&"s".to_string()));
    let dna = engine.state().sequences.get("s").unwrap();
    let tfbs_features: Vec<_> = dna
        .features()
        .iter()
        .filter(|f| {
            f.qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case("tfbs"))
        })
        .collect();
    assert!(!tfbs_features.is_empty());
    assert!(tfbs_features.iter().all(|f| {
        let motif_len = f
            .qualifier_values("motif_length_bp")
            .next()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(0);
        let span_matches = f
            .location
            .find_bounds()
            .ok()
            .map(|(from, to)| (to - from) as usize == motif_len)
            .unwrap_or(false);
        span_matches
            && motif_len == 4
            && f.qualifier_values("llr_bits").next().is_some()
            && f.qualifier_values("llr_quantile").next().is_some()
            && f.qualifier_values("true_log_odds_bits").next().is_some()
            && f.qualifier_values("true_log_odds_quantile")
                .next()
                .is_some()
            && f.qualifier_values("quantile_scope").next().is_some()
    }));
}

#[test]
fn test_annotate_tfbs_progress_reaches_completion() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), seq("TTTACGTAAACGTGGG"));
    let mut engine = GentleEngine::from_state(state);
    let mut progress_events: Vec<TfbsProgress> = vec![];
    let res = engine
        .apply_with_progress(
            Operation::AnnotateTfbs {
                seq_id: "s".to_string(),
                motifs: vec!["ACGT".to_string(), "TATAAA".to_string()],
                min_llr_bits: Some(0.0),
                min_llr_quantile: Some(0.0),
                per_tf_thresholds: vec![],
                clear_existing: Some(true),
                max_hits: None,
            },
            |progress| {
                if let OperationProgress::Tfbs(p) = progress {
                    progress_events.push(p);
                }
                true
            },
        )
        .unwrap();
    assert!(res.changed_seq_ids.contains(&"s".to_string()));
    assert!(!progress_events.is_empty());
    let last = progress_events.last().unwrap();
    assert_eq!(last.motif_index, last.motif_count);
    assert!((last.motif_percent - 100.0).abs() < f64::EPSILON);
    assert!((last.total_percent - 100.0).abs() < f64::EPSILON);
}

#[test]
fn test_annotate_tfbs_per_tf_override_changes_quantile_threshold() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), seq("TTTACGTAAACGTGGG"));
    let mut engine = GentleEngine::from_state(state);
    let strict = engine
        .apply(Operation::AnnotateTfbs {
            seq_id: "s".to_string(),
            motifs: vec!["ACGT".to_string()],
            min_llr_bits: Some(0.0),
            min_llr_quantile: Some(1.0),
            per_tf_thresholds: vec![],
            clear_existing: Some(true),
            max_hits: None,
        })
        .unwrap();
    let strict_count = engine
        .state()
        .sequences
        .get("s")
        .unwrap()
        .features()
        .iter()
        .filter(|f| {
            f.qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case("tfbs"))
        })
        .count();
    assert!(strict.changed_seq_ids.contains(&"s".to_string()));

    let relaxed = engine
        .apply(Operation::AnnotateTfbs {
            seq_id: "s".to_string(),
            motifs: vec!["ACGT".to_string()],
            min_llr_bits: Some(0.0),
            min_llr_quantile: Some(1.0),
            per_tf_thresholds: vec![TfThresholdOverride {
                tf: "ACGT".to_string(),
                min_llr_bits: None,
                min_llr_quantile: Some(0.0),
            }],
            clear_existing: Some(true),
            max_hits: None,
        })
        .unwrap();
    let relaxed_count = engine
        .state()
        .sequences
        .get("s")
        .unwrap()
        .features()
        .iter()
        .filter(|f| {
            f.qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case("tfbs"))
        })
        .count();
    assert!(relaxed.changed_seq_ids.contains(&"s".to_string()));
    assert!(relaxed_count >= strict_count);
}

#[test]
fn test_annotate_tfbs_max_hits_cap_and_unlimited() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), seq("ACGTACGTACGTACGTACGT"));
    let mut engine = GentleEngine::from_state(state);

    let capped = engine
        .apply(Operation::AnnotateTfbs {
            seq_id: "s".to_string(),
            motifs: vec!["ACGT".to_string()],
            min_llr_bits: Some(0.0),
            min_llr_quantile: Some(0.0),
            per_tf_thresholds: vec![],
            clear_existing: Some(true),
            max_hits: Some(2),
        })
        .unwrap();
    assert!(capped.changed_seq_ids.contains(&"s".to_string()));
    let capped_count = engine
        .state()
        .sequences
        .get("s")
        .unwrap()
        .features()
        .iter()
        .filter(|f| {
            f.qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case("tfbs"))
        })
        .count();
    assert_eq!(capped_count, 2);

    let unlimited = engine
        .apply(Operation::AnnotateTfbs {
            seq_id: "s".to_string(),
            motifs: vec!["ACGT".to_string()],
            min_llr_bits: Some(0.0),
            min_llr_quantile: Some(0.0),
            per_tf_thresholds: vec![],
            clear_existing: Some(true),
            max_hits: Some(0),
        })
        .unwrap();
    assert!(unlimited.changed_seq_ids.contains(&"s".to_string()));
    let unlimited_count = engine
        .state()
        .sequences
        .get("s")
        .unwrap()
        .features()
        .iter()
        .filter(|f| {
            f.qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case("tfbs"))
        })
        .count();
    assert!(unlimited_count > capped_count);
}

#[test]
fn test_inspect_tfbs_feature_expert_view() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), seq("TTTACGTAAACGTGGG"));
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::AnnotateTfbs {
            seq_id: "s".to_string(),
            motifs: vec!["ACGT".to_string()],
            min_llr_bits: Some(0.0),
            min_llr_quantile: Some(0.0),
            per_tf_thresholds: vec![],
            clear_existing: Some(true),
            max_hits: Some(1),
        })
        .unwrap();
    let feature_id = engine
        .state()
        .sequences
        .get("s")
        .unwrap()
        .features()
        .iter()
        .position(|feature| {
            feature
                .qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case("tfbs"))
        })
        .expect("a tfbs feature should exist");
    let view = engine
        .inspect_feature_expert("s", &FeatureExpertTarget::TfbsFeature { feature_id })
        .unwrap();
    match view {
        FeatureExpertView::Tfbs(tfbs) => {
            assert_eq!(tfbs.seq_id, "s");
            assert_eq!(tfbs.feature_id, feature_id);
            assert_eq!(tfbs.motif_length, 4);
            assert_eq!(tfbs.columns.len(), 4);
            assert_eq!(tfbs.matched_sequence.len(), 4);
            assert_eq!(tfbs.instruction, TFBS_EXPERT_INSTRUCTION);
            assert!(
                tfbs.columns
                    .iter()
                    .all(|column| column.information_content_bits.is_finite())
            );
        }
        other => panic!("expected tfbs expert view, got {other:?}"),
    }
}

#[test]
fn test_inspect_restriction_site_expert_view() {
    let mut dna = seq("AAGAATTCTT");
    *dna.restriction_enzymes_mut() = active_restriction_enzymes();
    dna.update_computed_features();
    let (key, names) = dna
        .restriction_enzyme_groups()
        .iter()
        .find(|(_, names)| names.iter().any(|name| name.eq_ignore_ascii_case("EcoRI")))
        .map(|(key, names)| (key.clone(), names.clone()))
        .expect("EcoRI site should exist in test sequence");
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), dna);
    let engine = GentleEngine::from_state(state);
    let view = engine
        .inspect_feature_expert(
            "s",
            &FeatureExpertTarget::RestrictionSite {
                cut_pos_1based: key.pos() as usize + 1,
                enzyme: Some("EcoRI".to_string()),
                recognition_start_1based: Some(key.from() as usize + 1),
                recognition_end_1based: Some(key.to() as usize),
            },
        )
        .unwrap();
    match view {
        FeatureExpertView::RestrictionSite(re) => {
            assert_eq!(re.seq_id, "s");
            assert_eq!(re.cut_pos_1based, key.pos() as usize + 1);
            assert_eq!(re.paired_cut_pos_1based, key.mate_pos() as usize + 1);
            assert_eq!(re.recognition_start_1based, key.from() as usize + 1);
            assert_eq!(re.recognition_end_1based, key.to() as usize);
            assert_eq!(re.cut_index_0based, 1);
            assert_eq!(re.paired_cut_index_0based, 5);
            assert_eq!(re.end_geometry, "5prime_overhang");
            assert_eq!(re.selected_enzyme.as_deref(), Some("EcoRI"));
            assert_eq!(re.recognition_iupac.as_deref(), Some("GAATTC"));
            assert_eq!(re.enzyme_cut_offset_0based, Some(1));
            assert_eq!(re.overlap_bp, Some(4));
            assert_eq!(re.enzyme_note, None);
            assert_eq!(
                re.rebase_url.as_deref(),
                Some("https://rebase.neb.com/rebase/enz/EcoRI.html")
            );
            assert!(
                re.enzyme_names
                    .iter()
                    .any(|name| name.eq_ignore_ascii_case("EcoRI"))
            );
            for name in names {
                assert!(
                    re.enzyme_names
                        .iter()
                        .any(|entry| entry.eq_ignore_ascii_case(&name))
                );
            }
            assert_eq!(re.instruction, RESTRICTION_EXPERT_INSTRUCTION);
        }
        other => panic!("expected restriction expert view, got {other:?}"),
    }
}

#[test]
fn test_inspect_splicing_feature_expert_view() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), splicing_test_sequence());
    let engine = GentleEngine::from_state(state);
    let view = engine
        .inspect_feature_expert(
            "s",
            &FeatureExpertTarget::SplicingFeature {
                feature_id: 0,
                scope: SplicingScopePreset::AllOverlappingBothStrands,
            },
        )
        .unwrap();
    match view {
        FeatureExpertView::Splicing(splicing) => {
            assert_eq!(splicing.seq_id, "s");
            assert_eq!(splicing.group_label, "GENE1");
            assert_eq!(splicing.transcript_count, 2);
            assert_eq!(splicing.transcripts.len(), 2);
            assert!(splicing.unique_exon_count >= 3);
            assert_eq!(splicing.matrix_rows.len(), 2);
            assert!(!splicing.unique_exons.is_empty());
            assert!(!splicing.boundaries.is_empty());
            assert!(!splicing.junctions.is_empty());
            assert!(
                splicing
                    .boundaries
                    .iter()
                    .any(|m| m.side.eq_ignore_ascii_case("donor") && m.canonical)
            );
            assert!(
                splicing
                    .boundaries
                    .iter()
                    .any(|m| m.side.eq_ignore_ascii_case("acceptor") && m.canonical)
            );
            assert_eq!(splicing.instruction, SPLICING_EXPERT_INSTRUCTION);
        }
        other => panic!("expected splicing expert view, got {other:?}"),
    }
}

#[test]
fn test_splicing_expert_view_reports_cds_flank_phase_forward() {
    let mut dna = DNAsequence::from_sequence(&"A".repeat(64)).expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(0, 5),
            gb_io::seq::Location::simple_range(10, 14),
            gb_io::seq::Location::simple_range(20, 26),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("TX_PHASE_PLUS".to_string())),
            ("label".into(), Some("TX_PHASE_PLUS".to_string())),
            (
                "cds_ranges_1based".into(),
                Some("1-5,11-14,21-26".to_string()),
            ),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), dna);
    let engine = GentleEngine::from_state(state);
    let view = engine
        .inspect_feature_expert(
            "s",
            &FeatureExpertTarget::SplicingFeature {
                feature_id: 0,
                scope: SplicingScopePreset::AllOverlappingBothStrands,
            },
        )
        .expect("inspect splicing");
    let FeatureExpertView::Splicing(splicing) = view else {
        panic!("expected splicing view");
    };
    let lane = &splicing.transcripts[0];
    assert_eq!(lane.exon_cds_phases.len(), 3);
    assert_eq!(lane.exon_cds_phases[0].left_cds_phase, Some(0));
    assert_eq!(lane.exon_cds_phases[0].right_cds_phase, Some(1));
    assert_eq!(lane.exon_cds_phases[1].left_cds_phase, Some(2));
    assert_eq!(lane.exon_cds_phases[1].right_cds_phase, Some(2));
    assert_eq!(lane.exon_cds_phases[2].left_cds_phase, Some(0));
    assert_eq!(lane.exon_cds_phases[2].right_cds_phase, Some(2));
}

#[test]
fn test_splicing_expert_view_reports_cds_flank_phase_reverse() {
    let mut dna = DNAsequence::from_sequence(&"A".repeat(64)).expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::Complement(Box::new(gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(0, 5),
            gb_io::seq::Location::simple_range(10, 14),
            gb_io::seq::Location::simple_range(20, 26),
        ]))),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("TX_PHASE_MINUS".to_string())),
            ("label".into(), Some("TX_PHASE_MINUS".to_string())),
            (
                "cds_ranges_1based".into(),
                Some("1-5,11-14,21-26".to_string()),
            ),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), dna);
    let engine = GentleEngine::from_state(state);
    let view = engine
        .inspect_feature_expert(
            "s",
            &FeatureExpertTarget::SplicingFeature {
                feature_id: 0,
                scope: SplicingScopePreset::AllOverlappingBothStrands,
            },
        )
        .expect("inspect splicing");
    let FeatureExpertView::Splicing(splicing) = view else {
        panic!("expected splicing view");
    };
    let lane = &splicing.transcripts[0];
    assert_eq!(lane.exon_cds_phases.len(), 3);
    assert_eq!(lane.exon_cds_phases[0].left_cds_phase, Some(2));
    assert_eq!(lane.exon_cds_phases[0].right_cds_phase, Some(1));
    assert_eq!(lane.exon_cds_phases[1].left_cds_phase, Some(0));
    assert_eq!(lane.exon_cds_phases[1].right_cds_phase, Some(0));
    assert_eq!(lane.exon_cds_phases[2].left_cds_phase, Some(2));
    assert_eq!(lane.exon_cds_phases[2].right_cds_phase, Some(0));
}

#[test]
fn test_is_mrna_feature_accepts_transcript_alias() {
    let feature = gb_io::seq::Feature {
        kind: "transcript".into(),
        location: gb_io::seq::Location::simple_range(1, 10),
        qualifiers: vec![],
    };
    assert!(GentleEngine::is_mrna_feature(&feature));
}

#[test]
fn test_derive_transcript_sequences_derives_all_mrna_features() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let source = engine
        .state()
        .sequences
        .get("s")
        .expect("source sequence")
        .get_forward_string();
    let result = engine
        .apply(Operation::DeriveTranscriptSequences {
            seq_id: "s".to_string(),
            feature_ids: vec![],
            scope: None,
            output_prefix: Some("tx".to_string()),
        })
        .expect("derive transcripts");
    assert_eq!(result.created_seq_ids.len(), 2);
    for derived_seq_id in result.created_seq_ids {
        let derived = engine
            .state()
            .sequences
            .get(&derived_seq_id)
            .expect("derived sequence");
        let features = derived.features();
        assert!(
            features
                .iter()
                .any(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        );
        let exon_count = features
            .iter()
            .filter(|feature| feature.kind.to_string().eq_ignore_ascii_case("exon"))
            .count();
        assert!(exon_count >= 2);
        for feature in features {
            if feature.kind.to_string().eq_ignore_ascii_case("mRNA") {
                let source_seq = feature
                    .qualifier_values("source_seq_id")
                    .next()
                    .unwrap_or_default()
                    .to_string();
                assert_eq!(source_seq, "s");
            }
        }
        assert!(
            source.len() >= derived.len(),
            "derived transcript should not exceed source length"
        );
    }
}

#[test]
fn test_derive_transcript_sequences_reverse_strand_uses_reverse_complement() {
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::Complement(Box::new(gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(2, 6),
            gb_io::seq::Location::simple_range(10, 14),
        ]))),
        qualifiers: vec![
            ("gene".into(), Some("GENE1".to_string())),
            ("transcript_id".into(), Some("TX_MINUS".to_string())),
            ("label".into(), Some("TX_MINUS".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);
    let source = engine
        .state()
        .sequences
        .get("s")
        .expect("source")
        .get_forward_string()
        .to_ascii_uppercase();
    let joined = format!("{}{}", &source[2..6], &source[10..14]);
    let expected = GentleEngine::reverse_complement(&joined);
    let result = engine
        .apply(Operation::DeriveTranscriptSequences {
            seq_id: "s".to_string(),
            feature_ids: vec![0],
            scope: None,
            output_prefix: Some("tx".to_string()),
        })
        .expect("derive reverse transcript");
    assert_eq!(result.created_seq_ids.len(), 1);
    let derived = engine
        .state()
        .sequences
        .get(&result.created_seq_ids[0])
        .expect("derived");
    assert_eq!(derived.get_forward_string(), expected);
}

#[test]
fn test_render_feature_expert_svg_operation() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), seq("TTTACGTAAACGTGGG"));
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::AnnotateTfbs {
            seq_id: "s".to_string(),
            motifs: vec!["ACGT".to_string()],
            min_llr_bits: Some(0.0),
            min_llr_quantile: Some(0.0),
            per_tf_thresholds: vec![],
            clear_existing: Some(true),
            max_hits: Some(1),
        })
        .unwrap();
    let feature_id = engine
        .state()
        .sequences
        .get("s")
        .unwrap()
        .features()
        .iter()
        .position(|feature| {
            feature
                .qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case("tfbs"))
        })
        .expect("tfbs feature");
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("feature.expert.svg");
    let path_text = path.display().to_string();
    let result = engine
        .apply(Operation::RenderFeatureExpertSvg {
            seq_id: "s".to_string(),
            target: FeatureExpertTarget::TfbsFeature { feature_id },
            path: path_text.clone(),
        })
        .unwrap();
    assert!(result.messages.iter().any(|m| m.contains("feature expert")));
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("<svg"));
    assert!(text.contains("TFBS expert"));
}

#[test]
fn test_render_splicing_feature_expert_svg_operation() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("splicing.feature.expert.svg");
    let path_text = path.display().to_string();
    let result = engine
        .apply(Operation::RenderFeatureExpertSvg {
            seq_id: "s".to_string(),
            target: FeatureExpertTarget::SplicingFeature {
                feature_id: 0,
                scope: SplicingScopePreset::AllOverlappingBothStrands,
            },
            path: path_text.clone(),
        })
        .unwrap();
    assert!(result.messages.iter().any(|m| m.contains("feature expert")));
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("<svg"));
    assert!(text.contains("Splicing expert"));
    assert!(text.contains("Junction transition support"));
    assert!(text.contains("cell color intensity encodes exon support frequency"));
    assert!(text.contains("exon transition matrix"));
    assert!(text.contains("length %3"));
}

#[test]
fn test_validate_isoform_panel_resource_reports_summary() {
    let report = GentleEngine::validate_isoform_panel_resource(
        "assets/panels/tp53_isoforms_v1.json",
        Some("tp53_isoforms_v1"),
    )
    .expect("validate panel");
    assert_eq!(report.schema, "gentle.isoform_panel_validation_report.v1");
    assert_eq!(report.panel_id, "tp53_isoforms_v1");
    assert_eq!(report.gene_symbol, "TP53");
    assert!(report.isoform_count >= 1);
    assert!(report.transcript_probe_count >= 1);
    assert!(matches!(report.status.as_str(), "ok" | "warning"));
}

#[test]
fn test_import_isoform_panel_allows_unmapped_when_strict_false_and_no_mrna_features() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), seq(&"ATGC".repeat(80)));
    let mut engine = GentleEngine::from_state(state);
    let import = engine
        .apply(Operation::ImportIsoformPanel {
            seq_id: "s".to_string(),
            panel_path: "assets/panels/tp53_isoforms_v1.json".to_string(),
            panel_id: Some("tp53_isoforms_v1".to_string()),
            strict: false,
        })
        .expect("strict=false import should succeed even without transcript mapping");
    assert!(
        import
            .messages
            .iter()
            .any(|m| m.contains("Imported isoform panel"))
    );
    let view = engine
        .inspect_feature_expert(
            "s",
            &FeatureExpertTarget::IsoformArchitecture {
                panel_id: "tp53_isoforms_v1".to_string(),
            },
        )
        .expect("inspect isoform architecture");
    match view {
        FeatureExpertView::IsoformArchitecture(isoform) => {
            assert!(!isoform.transcript_lanes.is_empty());
            assert!(
                isoform
                    .warnings
                    .iter()
                    .any(|w| w.contains("No mRNA/transcript features"))
            );
        }
        other => panic!("expected isoform architecture view, got {other:?}"),
    }
}

#[test]
fn test_isoform_panel_cds_geometry_mode_uses_cds_ranges_when_available() {
    let mut dna = seq(&"ATGC".repeat(120));
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(0, 30),
            gb_io::seq::Location::simple_range(40, 80),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("TP53".to_string())),
            ("transcript_id".into(), Some("TX1".to_string())),
            ("label".into(), Some("TX1".to_string())),
            ("cds_ranges_1based".into(), Some("6-20,51-68".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let tmp = tempdir().expect("tempdir");
    let panel_path = tmp.path().join("panel_cds.json");
    fs::write(
        &panel_path,
        r##"{
  "schema": "gentle.isoform_panel_resource.v1",
  "panel_id": "tp53_cds_demo",
  "gene_symbol": "TP53",
  "transcript_geometry_mode": "cds",
  "isoforms": [
    {
      "isoform_id": "alpha",
      "label": "TAp53α",
      "transcript_ids": ["TX1"],
      "domains": [{"name": "DBD", "start_aa": 10, "end_aa": 50, "color_hex": "#2563eb"}]
    }
  ]
}"##,
    )
    .expect("write panel");

    engine
        .apply(Operation::ImportIsoformPanel {
            seq_id: "s".to_string(),
            panel_path: panel_path.display().to_string(),
            panel_id: None,
            strict: true,
        })
        .expect("import panel");
    let view = engine
        .inspect_feature_expert(
            "s",
            &FeatureExpertTarget::IsoformArchitecture {
                panel_id: "tp53_cds_demo".to_string(),
            },
        )
        .expect("inspect isoform");
    match view {
        FeatureExpertView::IsoformArchitecture(isoform) => {
            assert_eq!(isoform.transcript_geometry_mode, "cds");
            assert_eq!(isoform.transcript_lanes.len(), 1);
            assert_eq!(isoform.transcript_lanes[0].exons[0].start_1based, 6);
            assert_eq!(isoform.transcript_lanes[0].exons[0].end_1based, 20);
            assert_eq!(isoform.transcript_lanes[0].exons[1].start_1based, 51);
            assert_eq!(isoform.transcript_lanes[0].exons[1].end_1based, 68);
            assert_eq!(isoform.transcript_lanes[0].cds_to_protein_segments.len(), 2);
            assert_eq!(
                isoform.transcript_lanes[0].cds_to_protein_segments[0].aa_start,
                1
            );
            assert_eq!(
                isoform.transcript_lanes[0].cds_to_protein_segments[0].aa_end,
                5
            );
            assert_eq!(
                isoform.transcript_lanes[0].cds_to_protein_segments[1].aa_start,
                6
            );
            assert_eq!(
                isoform.transcript_lanes[0].cds_to_protein_segments[1].aa_end,
                11
            );
        }
        other => panic!("expected isoform architecture view, got {other:?}"),
    }
}

#[test]
fn test_cds_ranges_to_reference_aa_segments_respects_reverse_transcript_order() {
    let segments = GentleEngine::cds_ranges_to_reference_aa_segments(
        vec![(100, 112), (200, 218)],
        true,
        Some(160),
    );
    assert_eq!(segments.len(), 2);
    assert_eq!(segments[0].genomic_start_1based, 201);
    assert_eq!(segments[0].genomic_end_1based, 218);
    assert_eq!(segments[0].aa_start, 160);
    assert_eq!(segments[0].aa_end, 165);
    assert_eq!(segments[1].genomic_start_1based, 101);
    assert_eq!(segments[1].genomic_end_1based, 112);
    assert_eq!(segments[1].aa_start, 166);
    assert_eq!(segments[1].aa_end, 169);
}

#[test]
fn test_transcript_feature_from_genome_record_includes_cds_ranges_qualifier() {
    let record = GenomeTranscriptRecord {
        chromosome: "chr17".to_string(),
        transcript_id: "TX1".to_string(),
        gene_id: Some("GENE1".to_string()),
        gene_name: Some("GENE1".to_string()),
        strand: Some('+'),
        transcript_start_1based: 100,
        transcript_end_1based: 180,
        exons_1based: vec![(100, 130), (151, 180)],
        cds_1based: vec![(105, 120), (160, 175)],
    };
    let feature =
        GentleEngine::transcript_feature_from_genome_record(&record, 100, 180).expect("feature");
    let cds = GentleEngine::feature_qualifier_text(&feature, "cds_ranges_1based")
        .expect("cds ranges qualifier");
    assert_eq!(cds, "6-21,61-76");
}

#[test]
fn test_validate_isoform_panel_resource_detects_curation_issues() {
    let tmp = tempdir().expect("temp dir");
    let panel_path = tmp.path().join("panel.json");
    let raw = r##"{
  "schema": "gentle.isoform_panel_resource.v1",
  "panel_id": "demo_panel",
  "gene_symbol": "DEMO1",
  "isoforms": [
    {
      "isoform_id": "alpha",
      "transcript_ids": ["ENST000001.1"],
      "expected_length_aa": 100,
      "domains": [
        {"name": "dbd", "start_aa": 50, "end_aa": 90, "color_hex": "#12GG00"},
        {"name": "ctd", "start_aa": 80, "end_aa": 120, "color_hex": "#112233"}
      ]
    },
    {
      "isoform_id": "beta",
      "transcript_ids": ["ENST000001.1"],
      "domains": []
    },
    {
      "isoform_id": "alpha",
      "transcript_ids": ["ENST000003.1"],
      "domains": []
    }
  ]
}"##;
    fs::write(&panel_path, raw).expect("write panel fixture");
    let report = GentleEngine::validate_isoform_panel_resource(panel_path.to_str().unwrap(), None)
        .expect("validate panel");

    assert_eq!(report.status, "warning");
    let mut codes = report
        .issues
        .iter()
        .map(|issue| issue.code.as_str())
        .collect::<Vec<_>>();
    codes.sort();
    assert!(codes.contains(&"invalid_color_hex"));
    assert!(codes.contains(&"overlapping_domains"));
    assert!(codes.contains(&"expected_length_below_domain_end"));
    assert!(codes.contains(&"duplicate_isoform_id"));
    assert!(codes.contains(&"shared_transcript_probe"));
    assert!(codes.contains(&"missing_domains"));
}

#[test]
fn test_render_sequence_svg_operation() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), seq(&"ATGC".repeat(40)));
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("svg");
    let path_text = path.display().to_string();
    let res = engine
        .apply(Operation::RenderSequenceSvg {
            seq_id: "s".to_string(),
            mode: RenderSvgMode::Linear,
            path: path_text.clone(),
        })
        .unwrap();
    assert!(res.messages.iter().any(|m| m.contains("SVG")));
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("<svg"));
}

#[test]
fn test_render_sequence_svg_operation_honors_linear_viewport() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("s".to_string(), seq(&"ATGC".repeat(40)));
    state.display.linear_view_start_bp = 40;
    state.display.linear_view_span_bp = 20;
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("svg");
    let path_text = path.display().to_string();
    engine
        .apply(Operation::RenderSequenceSvg {
            seq_id: "s".to_string(),
            mode: RenderSvgMode::Linear,
            path: path_text.clone(),
        })
        .unwrap();
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("41..60 (20 bp view of 160 bp)"));
}

#[test]
fn test_render_dotplot_svg_operation() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "query".to_string(),
        DNAsequence::from_sequence("ATGCATGCATGCATGC").expect("query sequence"),
    );
    state.sequences.insert(
        "ref".to_string(),
        DNAsequence::from_sequence("TTTATGCATGCATGCAAAA").expect("reference sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::ComputeDotplot {
            seq_id: "query".to_string(),
            reference_seq_id: Some("ref".to_string()),
            span_start_0based: Some(0),
            span_end_0based: Some(16),
            reference_span_start_0based: Some(3),
            reference_span_end_0based: Some(19),
            mode: DotplotMode::PairForward,
            word_size: 4,
            step_bp: 1,
            max_mismatches: 0,
            tile_bp: None,
            store_as: Some("pair_dotplot".to_string()),
        })
        .expect("compute dotplot");

    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("dotplot.svg");
    let path_text = path.display().to_string();
    let res = engine
        .apply(Operation::RenderDotplotSvg {
            seq_id: "query".to_string(),
            dotplot_id: "pair_dotplot".to_string(),
            path: path_text.clone(),
            flex_track_id: None,
            display_density_threshold: Some(0.0),
            display_intensity_gain: Some(1.0),
        })
        .expect("render dotplot svg");
    assert!(res.messages.iter().any(|m| m.contains("dotplot SVG")));
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("<svg"));
    assert!(text.contains("Dotplot workspace export"));
    assert!(text.contains("overlap by 3 bp"));
    assert!(text.contains("4 consecutive ordered windows"));
}

#[test]
fn test_render_lineage_svg_operation() {
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), seq("ATGCCA"));
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("svg");
    let path_text = path.display().to_string();
    let res = engine
        .apply(Operation::RenderLineageSvg {
            path: path_text.clone(),
        })
        .unwrap();
    assert!(res.messages.iter().any(|m| m.contains("lineage SVG")));
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("<svg"));
    let height_marker = "height=\"";
    let height_start = text
        .find(height_marker)
        .expect("svg height attribute should exist")
        + height_marker.len();
    let height_rest = &text[height_start..];
    let height_end = height_rest
        .find('"')
        .expect("svg height attribute should terminate");
    let height_value = height_rest[..height_end]
        .parse::<f32>()
        .expect("svg height should parse as number");
    assert!(
        height_value < 400.0,
        "simple lineage export should size canvas to content, got height={height_value}"
    );
}

#[test]
fn test_render_lineage_svg_includes_arrangement_nodes() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("a".to_string(), seq(&"ATGC".repeat(40)));
    state
        .sequences
        .insert("b".to_string(), seq(&"ATGC".repeat(55)));
    let mut engine = GentleEngine::from_state(state);

    let mut container_ids: Vec<String> = engine
        .state()
        .container_state
        .containers
        .keys()
        .cloned()
        .collect();
    container_ids.sort();
    engine
        .apply(Operation::CreateArrangementSerial {
            container_ids: container_ids.into_iter().take(2).collect(),
            arrangement_id: Some("arr-viz".to_string()),
            name: Some("Digest run".to_string()),
            ladders: Some(vec!["NEB 100bp DNA Ladder".to_string()]),
        })
        .unwrap();

    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("lineage.arr.svg");
    let path_text = path.display().to_string();
    engine
        .apply(Operation::RenderLineageSvg {
            path: path_text.clone(),
        })
        .unwrap();
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("arr-viz"));
    assert!(text.contains("lanes=2"));
}

#[test]
fn test_render_lineage_svg_includes_macro_instance_nodes() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("a".to_string(), seq(&"ATGC".repeat(40)));
    let mut engine = GentleEngine::from_state(state);
    let reverse = engine
        .apply(Operation::Reverse {
            input: "a".to_string(),
            output_id: Some("a_rev".to_string()),
        })
        .expect("reverse");
    let macro_instance_id = engine.record_lineage_macro_instance(LineageMacroInstance {
        macro_instance_id: String::new(),
        routine_id: Some("sequence.reverse".to_string()),
        routine_title: Some("Reverse helper".to_string()),
        template_name: Some("reverse_helper".to_string()),
        run_id: "macro".to_string(),
        created_at_unix_ms: 0,
        bound_inputs: vec![LineageMacroPortBinding {
            port_id: "seq_id".to_string(),
            kind: "sequence".to_string(),
            required: true,
            cardinality: "one".to_string(),
            values: vec!["a".to_string()],
            description: Some("Input sequence".to_string()),
        }],
        bound_outputs: vec![LineageMacroPortBinding {
            port_id: "out_id".to_string(),
            kind: "sequence".to_string(),
            required: false,
            cardinality: "one".to_string(),
            values: vec!["a_rev".to_string()],
            description: Some("Output sequence".to_string()),
        }],
        expanded_op_ids: vec![reverse.op_id],
        status: MacroInstanceStatus::Ok,
        status_message: None,
    });
    assert!(!macro_instance_id.trim().is_empty());

    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("lineage.macro.svg");
    let path_text = path.display().to_string();
    engine
        .apply(Operation::RenderLineageSvg {
            path: path_text.clone(),
        })
        .unwrap();
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("Reverse helper"));
    assert!(text.contains("in:seq_id"));
    assert!(text.contains("out:out_id"));
}

#[test]
fn test_render_lineage_svg_projects_single_insert_gibson_as_operation_hub() {
    let mut state = ProjectState::default();
    let mut destination =
        DNAsequence::from_sequence("AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT")
            .expect("destination sequence");
    destination.set_name("destination_vector".to_string());
    destination.set_circular(true);
    state
        .sequences
        .insert("destination_vector".to_string(), destination);

    let mut insert = DNAsequence::from_sequence("ATGCGTACGTTAGCGTACGATCGTACGTAGCTAGCTAGCATCGATCGA")
        .expect("insert sequence");
    insert.set_name("insert_x_amplicon".to_string());
    insert.set_circular(false);
    state
        .sequences
        .insert("insert_x_amplicon".to_string(), insert);

    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().container_state.containers.insert(
        "container-pool".to_string(),
        Container {
            container_id: "container-pool".to_string(),
            kind: ContainerKind::Pool,
            name: Some("Unrelated pooled inputs".to_string()),
            members: vec![
                "destination_vector".to_string(),
                "insert_x_amplicon".to_string(),
            ],
            created_by_op: None,
            created_at_unix_ms: 1,
        },
    );
    engine
        .state_mut()
        .container_state
        .seq_to_latest_container
        .insert(
            "destination_vector".to_string(),
            "container-pool".to_string(),
        );
    engine
        .state_mut()
        .container_state
        .seq_to_latest_container
        .insert(
            "insert_x_amplicon".to_string(),
            "container-pool".to_string(),
        );
    let plan_json = r#"{
  "schema": "gentle.gibson_assembly_plan.v1",
  "id": "lineage_cli_projection_test",
  "title": "CLI projection test",
  "summary": "single insert",
  "destination": {
    "seq_id": "destination_vector",
    "topology_before_opening": "circular",
    "opening": {
      "mode": "defined_site",
      "label": "selected window",
      "start_0based": 12,
      "end_0based_exclusive": 18,
      "left_end_id": "dest_left",
      "right_end_id": "dest_right",
      "uniqueness_requirement": "must_be_unambiguous"
    }
  },
  "product": {"topology": "circular", "output_id_hint": "lineage_out"},
  "fragments": [
    {
      "id": "insert_x",
      "seq_id": "insert_x_amplicon",
      "role": "insert",
      "orientation": "forward",
      "left_end_strategy": {"mode": "primer_added_overlap", "target_junction_id": "junction_left"},
      "right_end_strategy": {"mode": "primer_added_overlap", "target_junction_id": "junction_right"}
    }
  ],
  "assembly_order": [
    {"kind": "destination_end", "id": "dest_left"},
    {"kind": "fragment", "id": "insert_x"},
    {"kind": "destination_end", "id": "dest_right"}
  ],
  "junctions": [
    {
      "id": "junction_left",
      "left_member": {"kind": "destination_end", "id": "dest_left"},
      "right_member": {"kind": "fragment", "id": "insert_x"},
      "required_overlap_bp": 20,
      "overlap_partition": {"left_member_bp": 20, "right_member_bp": 0},
      "overlap_source": "derive_from_destination_left_flank",
      "distinct_from": ["junction_right"]
    },
    {
      "id": "junction_right",
      "left_member": {"kind": "fragment", "id": "insert_x"},
      "right_member": {"kind": "destination_end", "id": "dest_right"},
      "required_overlap_bp": 20,
      "overlap_partition": {"left_member_bp": 0, "right_member_bp": 20},
      "overlap_source": "derive_from_destination_right_flank",
      "distinct_from": ["junction_left"]
    }
  ],
  "validation_policy": {
    "require_unambiguous_destination_opening": true,
    "require_distinct_terminal_junctions": true,
    "adjacency_overlap_mismatch": "error",
    "design_targets": {
      "overlap_bp_min": 18,
      "overlap_bp_max": 24,
      "minimum_overlap_tm_celsius": 48.0,
      "priming_segment_tm_min_celsius": 48.0,
      "priming_segment_tm_max_celsius": 70.0,
      "priming_segment_min_length_bp": 18,
      "priming_segment_max_length_bp": 30,
      "max_anneal_hits": 4
    },
    "uniqueness_checks": {
      "destination_context": "warn",
      "participating_fragments": "warn",
      "reference_contexts": []
    }
  }
}"#;
    let apply_result = engine
        .apply(Operation::ApplyGibsonAssemblyPlan {
            plan_json: plan_json.to_string(),
        })
        .expect("apply Gibson operation");
    let gibson_arrangements = engine
        .state()
        .container_state
        .arrangements
        .values()
        .collect::<Vec<_>>();
    assert_eq!(
        gibson_arrangements.len(),
        1,
        "one Gibson arrangement should be created"
    );
    let arrangement = gibson_arrangements[0];
    assert_eq!(arrangement.lane_container_ids.len(), 3);
    assert!(
        !arrangement.ladders.is_empty(),
        "Gibson arrangement should carry recommended ladders"
    );
    let lane_members = arrangement
        .lane_container_ids
        .iter()
        .map(|container_id| {
            engine
                .state()
                .container_state
                .containers
                .get(container_id)
                .map(|container| container.members.clone())
                .unwrap_or_default()
        })
        .collect::<Vec<_>>();
    assert_eq!(
        lane_members,
        vec![
            vec!["destination_vector".to_string()],
            vec!["insert_x_amplicon".to_string()],
            vec!["lineage_out".to_string()]
        ],
        "Gibson arrangement should use exact per-sequence lane containers, not unrelated pools"
    );

    let (nodes, edges) = build_lineage_svg_graph(engine.state(), engine.operation_log());
    let hub_nodes: Vec<_> = nodes
        .iter()
        .filter(|node| {
            node.kind == LineageSvgNodeKind::OperationHub && node.title == "Gibson cloning"
        })
        .collect();
    assert_eq!(
        hub_nodes.len(),
        1,
        "projected lineage should contain one Gibson hub"
    );
    let hub_node_id = hub_nodes[0].node_id.clone();
    assert_eq!(
        edges
            .iter()
            .filter(|edge| edge.from_node_id == hub_node_id || edge.to_node_id == hub_node_id)
            .count(),
        5,
        "single-insert Gibson should render 2 input edges and 3 output edges via the hub"
    );
    assert!(
        !edges.iter().any(|edge| edge.label == apply_result.op_id),
        "projected Gibson edges should use the shared Gibson label instead of the raw op id"
    );

    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("lineage.gibson.svg");
    let path_text = path.display().to_string();
    engine
        .apply(Operation::RenderLineageSvg {
            path: path_text.clone(),
        })
        .expect("render lineage svg");
    let text = std::fs::read_to_string(path_text).unwrap();
    assert_eq!(
        text.matches("Gibson cloning").count(),
        1,
        "hero export should label the Gibson operation once on the hub node"
    );
    assert!(text.contains("lineage_out"));
    assert!(!text.contains("op=op-3"));

    let temp = tempdir().expect("tempdir");
    let state_path = temp.path().join("lineage.state.json");
    engine
        .state()
        .save_to_path(state_path.to_string_lossy().as_ref())
        .expect("save state");
    let reloaded_state =
        ProjectState::load_from_path(state_path.to_string_lossy().as_ref()).expect("reload state");
    let (reloaded_nodes, _) = build_lineage_svg_graph(&reloaded_state, &[]);
    assert_eq!(
        reloaded_nodes
            .iter()
            .filter(|node| {
                node.kind == LineageSvgNodeKind::OperationHub && node.title == "Gibson cloning"
            })
            .count(),
        1,
        "saved-state lineage export should still infer the Gibson hub after reload"
    );

    let mut reloaded_engine = GentleEngine::from_state(reloaded_state);
    let reloaded_svg_path = temp.path().join("lineage.reloaded.gibson.svg");
    reloaded_engine
        .apply(Operation::RenderLineageSvg {
            path: reloaded_svg_path.display().to_string(),
        })
        .expect("render lineage svg after reload");
    let reloaded_text = std::fs::read_to_string(reloaded_svg_path).unwrap();
    assert_eq!(
        reloaded_text.matches("Gibson cloning").count(),
        1,
        "reloaded hero export should still label the Gibson operation once"
    );
    assert!(reloaded_text.contains("lineage_out"));
    assert!(!reloaded_text.contains("op=op-3"));
}

#[test]
fn test_render_pool_gel_svg_operation() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("a".to_string(), seq(&"ATGC".repeat(80)));
    state
        .sequences
        .insert("b".to_string(), seq(&"ATGC".repeat(150)));
    state
        .sequences
        .insert("c".to_string(), seq(&"ATGC".repeat(260)));
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("pool.gel.svg");
    let path_text = path.display().to_string();
    let res = engine
        .apply(Operation::RenderPoolGelSvg {
            inputs: vec!["a".to_string(), "b".to_string(), "c".to_string()],
            path: path_text.clone(),
            ladders: None,
            container_ids: None,
            arrangement_id: None,
        })
        .unwrap();
    assert!(res.messages.iter().any(|m| m.contains("serial gel SVG")));
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("<svg"));
    assert!(text.contains("Serial Gel Preview"));
}

#[test]
fn test_create_arrangement_serial_operation() {
    let mut state = ProjectState::default();
    state.sequences.insert("a".to_string(), seq("ATGCATGC"));
    state.sequences.insert("b".to_string(), seq("ATGCATGCATGC"));
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Tube A".to_string()),
            members: vec!["a".to_string()],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.containers.insert(
        "container-2".to_string(),
        Container {
            container_id: "container-2".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Tube B".to_string()),
            members: vec!["b".to_string()],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);
    let result = engine
        .apply(Operation::CreateArrangementSerial {
            container_ids: vec!["container-1".to_string(), "container-2".to_string()],
            arrangement_id: Some("arr-test".to_string()),
            name: Some("Digest run".to_string()),
            ladders: Some(vec!["NEB 1kb DNA Ladder".to_string()]),
        })
        .unwrap();
    assert!(
        result
            .messages
            .iter()
            .any(|m| m.contains("Created serial arrangement 'arr-test'"))
    );
    let arrangement = engine
        .state()
        .container_state
        .arrangements
        .get("arr-test")
        .expect("arrangement was created");
    assert_eq!(arrangement.mode, ArrangementMode::Serial);
    assert_eq!(
        arrangement.lane_container_ids,
        vec!["container-1".to_string(), "container-2".to_string()]
    );
    assert_eq!(arrangement.ladders, vec!["NEB 1kb DNA Ladder".to_string()]);
}

#[test]
fn test_render_pool_gel_svg_operation_from_containers_and_arrangement() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("a".to_string(), seq(&"ATGC".repeat(70)));
    state
        .sequences
        .insert("b".to_string(), seq(&"ATGC".repeat(110)));
    state
        .sequences
        .insert("c".to_string(), seq(&"ATGC".repeat(150)));
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Tube A".to_string()),
            members: vec!["a".to_string()],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.containers.insert(
        "container-2".to_string(),
        Container {
            container_id: "container-2".to_string(),
            kind: ContainerKind::Pool,
            name: Some("Tube B".to_string()),
            members: vec!["b".to_string(), "c".to_string()],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-1".to_string(),
        Arrangement {
            arrangement_id: "arr-1".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Run 1".to_string()),
            lane_container_ids: vec!["container-1".to_string(), "container-2".to_string()],
            ladders: vec!["NEB 100bp DNA Ladder".to_string()],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);

    let tmp_container = tempfile::NamedTempFile::new().unwrap();
    let path_container = tmp_container.path().with_extension("container.gel.svg");
    let path_container_text = path_container.display().to_string();
    let res_container = engine
        .apply(Operation::RenderPoolGelSvg {
            inputs: vec![],
            path: path_container_text.clone(),
            ladders: None,
            container_ids: Some(vec!["container-2".to_string()]),
            arrangement_id: None,
        })
        .unwrap();
    assert!(
        res_container
            .messages
            .iter()
            .any(|m| m.contains("serial gel SVG"))
    );
    let svg_container = std::fs::read_to_string(path_container_text).unwrap();
    assert!(svg_container.contains("Serial Gel Preview"));
    assert!(svg_container.contains("Tube B"));

    let tmp_arrangement = tempfile::NamedTempFile::new().unwrap();
    let path_arrangement = tmp_arrangement.path().with_extension("arrangement.gel.svg");
    let path_arrangement_text = path_arrangement.display().to_string();
    let res_arrangement = engine
        .apply(Operation::RenderPoolGelSvg {
            inputs: vec![],
            path: path_arrangement_text.clone(),
            ladders: None,
            container_ids: None,
            arrangement_id: Some("arr-1".to_string()),
        })
        .unwrap();
    assert!(
        res_arrangement
            .messages
            .iter()
            .any(|m| m.contains("2 sample lane(s)"))
    );
    let svg_arrangement = std::fs::read_to_string(path_arrangement_text).unwrap();
    assert!(svg_arrangement.contains("Serial Gel Preview"));
    assert!(svg_arrangement.contains("Tube A"));
    assert!(svg_arrangement.contains("Tube B"));
}

#[test]
fn test_render_pool_gel_svg_operation_missing_input_fails() {
    let mut state = ProjectState::default();
    state.sequences.insert("a".to_string(), seq("ATGC"));
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("pool.gel.svg");
    let path_text = path.display().to_string();
    let err = engine
        .apply(Operation::RenderPoolGelSvg {
            inputs: vec!["missing".to_string()],
            path: path_text,
            ladders: None,
            container_ids: None,
            arrangement_id: None,
        })
        .unwrap_err();
    assert!(err.message.contains("not found"));
}

#[test]
fn test_export_pool_operation() {
    let mut state = ProjectState::default();
    state.sequences.insert("a".to_string(), seq("ATGC"));
    state.sequences.insert("b".to_string(), seq("TTAA"));
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("pool.gentle.json");
    let path_text = path.display().to_string();
    let res = engine
        .apply(Operation::ExportPool {
            inputs: vec!["a".to_string(), "b".to_string()],
            path: path_text.clone(),
            pool_id: Some("pool_1".to_string()),
            human_id: Some("test pool".to_string()),
        })
        .unwrap();
    assert!(res.messages.iter().any(|m| m.contains("Wrote pool export")));
    let text = std::fs::read_to_string(path_text).unwrap();
    let v: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(v["schema"], "gentle.pool.v1");
    assert_eq!(v["pool_id"], "pool_1");
    assert_eq!(v["human_id"], "test pool");
    assert_eq!(v["member_count"], 2);
}

#[test]
fn test_export_process_run_bundle_operation() {
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), seq("ATGCCA"));
    state.metadata.insert(
        ROUTINE_DECISION_TRACES_METADATA_KEY.to_string(),
        serde_json::to_value(RoutineDecisionTraceStore {
            schema: ROUTINE_DECISION_TRACE_STORE_SCHEMA.to_string(),
            traces: vec![RoutineDecisionTrace {
                schema: ROUTINE_DECISION_TRACE_SCHEMA.to_string(),
                trace_id: "trace_1".to_string(),
                source: "gui_routine_assistant".to_string(),
                status: "executed".to_string(),
                created_at_unix_ms: 10,
                updated_at_unix_ms: 20,
                goal_text: "Assemble reporter".to_string(),
                query_text: "golden gate".to_string(),
                candidate_routine_ids: vec!["golden_gate.type_iis_single_insert".to_string()],
                selected_routine_id: Some("golden_gate.type_iis_single_insert".to_string()),
                selected_routine_title: Some("Golden Gate Type IIS Single Insert".to_string()),
                selected_routine_family: Some("golden_gate".to_string()),
                alternatives_presented: vec!["gibson.two_fragment_overlap_preview".to_string()],
                comparisons: vec![RoutineDecisionTraceComparison {
                    left_routine_id: "golden_gate.type_iis_single_insert".to_string(),
                    right_routine_id: "gibson.two_fragment_overlap_preview".to_string(),
                }],
                disambiguation_questions_presented: vec![
                    RoutineDecisionTraceDisambiguationQuestion {
                        question_id: "termini_expectation".to_string(),
                        question_text:
                            "What molecule types and termini are expected for insert and vector?"
                                .to_string(),
                    },
                    RoutineDecisionTraceDisambiguationQuestion {
                        question_id: "directionality_constraints".to_string(),
                        question_text:
                            "Do you need directionality constraints or compatible overhang control?"
                                .to_string(),
                    },
                ],
                disambiguation_answers: vec![RoutineDecisionTraceDisambiguationAnswer {
                    question_id: "termini_expectation".to_string(),
                    answer_text: "vector and insert are type IIS compatible".to_string(),
                }],
                bindings_snapshot: std::collections::BTreeMap::from([(
                    "vector_seq_id".to_string(),
                    "s".to_string(),
                )]),
                preflight_history: vec![RoutineDecisionTracePreflightSnapshot {
                    can_execute: true,
                    warnings: vec!["no additional warnings".to_string()],
                    errors: vec![],
                    contract_source: Some("routine_catalog".to_string()),
                }],
                preflight_snapshot: Some(RoutineDecisionTracePreflightSnapshot {
                    can_execute: true,
                    warnings: vec![],
                    errors: vec![],
                    contract_source: Some("routine_catalog".to_string()),
                }),
                execution_attempted: true,
                execution_success: Some(true),
                transactional: Some(true),
                macro_instance_id: Some("macro_1".to_string()),
                emitted_operation_ids: vec!["op_1".to_string()],
                execution_error: None,
                export_events: vec![RoutineDecisionTraceExportEvent {
                    run_bundle_path: "run_bundle.json".to_string(),
                    exported_at_unix_ms: 25,
                }],
            }],
        })
        .expect("trace store json"),
    );
    let mut engine = GentleEngine::from_state(state);
    let reverse = engine
        .apply(Operation::Reverse {
            input: "s".to_string(),
            output_id: Some("s_rev".to_string()),
        })
        .expect("reverse");
    let _set_param = engine
        .apply(Operation::SetParameter {
            name: "max_fragments_per_container".to_string(),
            value: serde_json::json!(12345),
        })
        .expect("set parameter");
    let tmp = tempfile::NamedTempFile::new().expect("tmp");
    let path = tmp.path().with_extension("run.bundle.json");
    let path_text = path.display().to_string();
    let export = engine
        .apply(Operation::ExportProcessRunBundle {
            path: path_text.clone(),
            run_id: Some("interactive".to_string()),
        })
        .expect("export run bundle");
    assert!(
        export
            .messages
            .iter()
            .any(|line| line.contains("process run bundle"))
    );

    let text = std::fs::read_to_string(path_text).expect("read bundle");
    let bundle: ProcessRunBundleExport =
        serde_json::from_str(&text).expect("parse run bundle json");
    assert_eq!(bundle.schema, PROCESS_RUN_BUNDLE_SCHEMA);
    assert_eq!(bundle.run_id_filter.as_deref(), Some("interactive"));
    assert!(bundle.selected_record_count >= 2);
    assert!(
        bundle
            .operation_log
            .iter()
            .any(|record| record.result.op_id == reverse.op_id)
    );
    assert!(
        bundle
            .parameter_overrides
            .iter()
            .any(|row| row.name == "max_fragments_per_container"
                && row.value == serde_json::json!(12345))
    );
    assert!(
        bundle
            .outputs
            .created_seq_ids
            .iter()
            .any(|seq_id| seq_id == "s_rev")
    );
    assert!(
        bundle
            .inputs
            .root_sequence_ids
            .iter()
            .any(|seq_id| seq_id == "s")
    );
    assert_eq!(bundle.decision_traces.len(), 1);
    assert_eq!(bundle.decision_traces[0].trace_id, "trace_1");
    assert_eq!(
        bundle.decision_traces[0]
            .disambiguation_questions_presented
            .len(),
        2
    );
    assert_eq!(bundle.decision_traces[0].disambiguation_answers.len(), 1);
    assert_eq!(bundle.decision_traces[0].preflight_history.len(), 1);
    assert_eq!(
        bundle.decision_traces[0]
            .preflight_snapshot
            .as_ref()
            .map(|row| row.contract_source.as_deref()),
        Some(Some("routine_catalog"))
    );
}

#[test]
fn test_export_process_run_bundle_decision_trace_partial_statuses_and_ordering() {
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), seq("ATGCCA"));
    state.metadata.insert(
        ROUTINE_DECISION_TRACES_METADATA_KEY.to_string(),
        serde_json::to_value(RoutineDecisionTraceStore {
            schema: ROUTINE_DECISION_TRACE_STORE_SCHEMA.to_string(),
            traces: vec![
                RoutineDecisionTrace {
                    schema: ROUTINE_DECISION_TRACE_SCHEMA.to_string(),
                    trace_id: "trace_z".to_string(),
                    source: "gui_routine_assistant".to_string(),
                    status: "execution_failed".to_string(),
                    created_at_unix_ms: 30,
                    updated_at_unix_ms: 31,
                    disambiguation_answers: vec![
                        RoutineDecisionTraceDisambiguationAnswer {
                            question_id: "question_b".to_string(),
                            answer_text: "answer b".to_string(),
                        },
                        RoutineDecisionTraceDisambiguationAnswer {
                            question_id: "question_a".to_string(),
                            answer_text: "answer a".to_string(),
                        },
                    ],
                    ..RoutineDecisionTrace::default()
                },
                RoutineDecisionTrace {
                    schema: ROUTINE_DECISION_TRACE_SCHEMA.to_string(),
                    trace_id: "trace_a".to_string(),
                    source: "gui_routine_assistant".to_string(),
                    status: "preflight_failed".to_string(),
                    created_at_unix_ms: 10,
                    updated_at_unix_ms: 11,
                    preflight_history: vec![RoutineDecisionTracePreflightSnapshot {
                        can_execute: false,
                        warnings: vec![],
                        errors: vec!["missing sequence".to_string()],
                        contract_source: Some("routine_catalog".to_string()),
                    }],
                    preflight_snapshot: None,
                    ..RoutineDecisionTrace::default()
                },
                RoutineDecisionTrace {
                    schema: ROUTINE_DECISION_TRACE_SCHEMA.to_string(),
                    trace_id: "trace_m".to_string(),
                    source: "gui_routine_assistant".to_string(),
                    status: "aborted".to_string(),
                    created_at_unix_ms: 10,
                    updated_at_unix_ms: 12,
                    ..RoutineDecisionTrace::default()
                },
            ],
        })
        .expect("trace store json"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::Reverse {
            input: "s".to_string(),
            output_id: Some("s_rev".to_string()),
        })
        .expect("reverse");
    let tmp = tempfile::NamedTempFile::new().expect("tmp");
    let path = tmp.path().with_extension("run.bundle.partial_traces.json");
    let path_text = path.display().to_string();
    let export = engine
        .apply(Operation::ExportProcessRunBundle {
            path: path_text.clone(),
            run_id: Some("interactive".to_string()),
        })
        .expect("export run bundle");
    assert!(
        export
            .messages
            .iter()
            .any(|line| line.contains("process run bundle"))
    );

    let text = std::fs::read_to_string(path_text).expect("read bundle");
    let bundle: ProcessRunBundleExport =
        serde_json::from_str(&text).expect("parse run bundle json");
    assert_eq!(bundle.decision_traces.len(), 3);
    assert_eq!(bundle.decision_traces[0].trace_id, "trace_a");
    assert_eq!(bundle.decision_traces[0].status, "preflight_failed");
    assert_eq!(bundle.decision_traces[1].trace_id, "trace_m");
    assert_eq!(bundle.decision_traces[1].status, "aborted");
    assert_eq!(bundle.decision_traces[2].trace_id, "trace_z");
    assert_eq!(bundle.decision_traces[2].status, "execution_failed");
    assert_eq!(
        bundle.decision_traces[0]
            .preflight_snapshot
            .as_ref()
            .expect("preflight snapshot from history")
            .errors,
        vec!["missing sequence".to_string()]
    );
    assert_eq!(
        bundle.decision_traces[2]
            .disambiguation_answers
            .iter()
            .map(|row| row.question_id.as_str())
            .collect::<Vec<_>>(),
        vec!["question_a", "question_b"]
    );
}

#[test]
fn test_export_process_run_bundle_run_id_not_found_fails() {
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), seq("ATGCCA"));
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::Reverse {
            input: "s".to_string(),
            output_id: Some("s_rev".to_string()),
        })
        .expect("reverse");
    let tmp = tempfile::NamedTempFile::new().expect("tmp");
    let path = tmp.path().with_extension("run.bundle.none.json");
    let path_text = path.display().to_string();
    let err = engine
        .apply(Operation::ExportProcessRunBundle {
            path: path_text,
            run_id: Some("missing_run".to_string()),
        })
        .expect_err("missing run should fail");
    assert!(matches!(err.code, ErrorCode::NotFound));
    assert!(err.message.contains("No operation records found"));
}

#[test]
fn test_inspect_dna_ladders() {
    let catalog = GentleEngine::inspect_dna_ladders(None);
    assert_eq!(catalog.schema, "gentle.dna_ladders.v1");
    assert!(catalog.ladder_count > 0);
    assert_eq!(catalog.ladder_count, catalog.ladders.len());
    assert!(
        catalog
            .ladders
            .iter()
            .any(|ladder| ladder.name == "NEB 100bp DNA Ladder")
    );
}

#[test]
fn test_export_dna_ladders_operation() {
    let mut engine = GentleEngine::new();
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("dna.ladders.json");
    let path_text = path.display().to_string();
    let res = engine
        .apply(Operation::ExportDnaLadders {
            path: path_text.clone(),
            name_filter: Some("NEB".to_string()),
        })
        .unwrap();
    assert!(
        res.messages
            .iter()
            .any(|m| m.contains("DNA ladders catalog"))
    );
    let text = std::fs::read_to_string(path_text).unwrap();
    let catalog: DnaLadderCatalog = serde_json::from_str(&text).unwrap();
    assert_eq!(catalog.schema, "gentle.dna_ladders.v1");
    assert!(catalog.ladder_count > 0);
    assert!(
        catalog
            .ladders
            .iter()
            .all(|ladder| ladder.name.to_ascii_lowercase().contains("neb"))
    );
}

#[test]
fn test_inspect_rna_ladders() {
    let catalog = GentleEngine::inspect_rna_ladders(None);
    assert_eq!(catalog.schema, "gentle.rna_ladders.v1");
    assert!(catalog.ladder_count > 0);
    assert_eq!(catalog.ladder_count, catalog.ladders.len());
    assert!(
        catalog
            .ladders
            .iter()
            .any(|ladder| ladder.name.contains("RNA"))
    );
}

#[test]
fn test_export_rna_ladders_operation() {
    let mut engine = GentleEngine::new();
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("rna.ladders.json");
    let path_text = path.display().to_string();
    let res = engine
        .apply(Operation::ExportRnaLadders {
            path: path_text.clone(),
            name_filter: Some("NEB".to_string()),
        })
        .unwrap();
    assert!(
        res.messages
            .iter()
            .any(|m| m.contains("RNA ladders catalog"))
    );
    let text = std::fs::read_to_string(path_text).unwrap();
    let catalog: RnaLadderCatalog = serde_json::from_str(&text).unwrap();
    assert_eq!(catalog.schema, "gentle.rna_ladders.v1");
    assert!(catalog.ladder_count > 0);
    assert!(
        catalog
            .ladders
            .iter()
            .all(|ladder| ladder.name.to_ascii_lowercase().contains("neb"))
    );
}

#[test]
fn test_save_file_operation_fasta() {
    let mut state = ProjectState::default();
    state.sequences.insert("s".to_string(), seq("ATGCCA"));
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("fa");
    let path_text = path.display().to_string();
    let res = engine
        .apply(Operation::SaveFile {
            seq_id: "s".to_string(),
            path: path_text.clone(),
            format: ExportFormat::Fasta,
        })
        .unwrap();
    assert!(res.changed_seq_ids.contains(&"s".to_string()));
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.starts_with(">"));
    assert!(text.contains("ATGCCA"));
}

#[test]
fn test_save_file_operation_fasta_includes_synthetic_metadata() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "s".to_string(),
        synth_oligo("molecule=ssdna topology=linear", b"ATGCCA"),
    );
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("fa");
    let path_text = path.display().to_string();
    engine
        .apply(Operation::SaveFile {
            seq_id: "s".to_string(),
            path: path_text.clone(),
            format: ExportFormat::Fasta,
        })
        .unwrap();
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("molecule=ssdna"));
    assert!(text.contains("topology=linear"));
}

#[test]
fn test_render_sequence_svg_operation_circular() {
    let mut state = ProjectState::default();
    let mut s = seq(&"ATGC".repeat(40));
    s.set_circular(true);
    state.sequences.insert("s".to_string(), s);
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("svg");
    let path_text = path.display().to_string();
    let res = engine
        .apply(Operation::RenderSequenceSvg {
            seq_id: "s".to_string(),
            mode: RenderSvgMode::Circular,
            path: path_text.clone(),
        })
        .unwrap();
    assert!(res.messages.iter().any(|m| m.contains("SVG")));
    let text = std::fs::read_to_string(path_text).unwrap();
    assert!(text.contains("<svg"));
}

#[test]
fn test_export_pool_operation_defaults() {
    let mut state = ProjectState::default();
    state.sequences.insert("a".to_string(), seq("ATGC"));
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("pool.gentle.json");
    let path_text = path.display().to_string();
    engine
        .apply(Operation::ExportPool {
            inputs: vec!["a".to_string()],
            path: path_text.clone(),
            pool_id: None,
            human_id: None,
        })
        .unwrap();
    let text = std::fs::read_to_string(path_text).unwrap();
    let v: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(v["pool_id"], "pool_export");
    assert!(v["human_id"].as_str().unwrap().starts_with("Pool("));
}

#[test]
fn test_export_pool_operation_empty_inputs_fails() {
    let mut engine = GentleEngine::new();
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("pool.gentle.json");
    let path_text = path.display().to_string();
    let err = engine
        .apply(Operation::ExportPool {
            inputs: vec![],
            path: path_text,
            pool_id: None,
            human_id: None,
        })
        .unwrap_err();
    assert!(err.message.contains("at least one input"));
}

#[test]
fn test_export_pool_operation_missing_sequence_fails() {
    let mut engine = GentleEngine::new();
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().with_extension("pool.gentle.json");
    let path_text = path.display().to_string();
    let err = engine
        .apply(Operation::ExportPool {
            inputs: vec!["missing".to_string()],
            path: path_text,
            pool_id: None,
            human_id: None,
        })
        .unwrap_err();
    assert!(err.message.contains("not found"));
}

#[test]
fn test_set_parameter_unknown_name_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "unknown_param".to_string(),
            value: serde_json::json!(1),
        })
        .unwrap_err();
    assert!(err.message.contains("Unknown parameter"));
}

#[test]
fn test_set_parameter_invalid_type_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "max_fragments_per_container".to_string(),
            value: serde_json::json!("not-a-number"),
        })
        .unwrap_err();
    assert!(err.message.contains("requires a positive integer"));
}

#[test]
fn test_set_parameter_zero_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "max_fragments_per_container".to_string(),
            value: serde_json::json!(0),
        })
        .unwrap_err();
    assert!(err.message.contains("must be >= 1"));
}

#[test]
fn test_set_parameter_tfbs_display_filter_fields() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::SetParameter {
            name: "show_tfbs".to_string(),
            value: serde_json::json!(true),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "tfbs_display_use_llr_bits".to_string(),
            value: serde_json::json!(false),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "tfbs_display_min_llr_bits".to_string(),
            value: serde_json::json!(-2.5),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "tfbs_display_use_true_log_odds_quantile".to_string(),
            value: serde_json::json!(true),
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "tfbs_display_min_true_log_odds_quantile".to_string(),
            value: serde_json::json!(0.8),
        })
        .unwrap();
    assert!(engine.state().display.show_tfbs);
    assert!(!engine.state().display.tfbs_display_use_llr_bits);
    assert!((engine.state().display.tfbs_display_min_llr_bits + 2.5).abs() < f64::EPSILON);
    assert!(
        engine
            .state()
            .display
            .tfbs_display_use_true_log_odds_quantile
    );
    assert!(
        (engine
            .state()
            .display
            .tfbs_display_min_true_log_odds_quantile
            - 0.8)
            .abs()
            < f64::EPSILON
    );
}

#[test]
fn test_set_parameter_tfbs_display_quantile_out_of_range_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "tfbs_display_min_llr_quantile".to_string(),
            value: serde_json::json!(1.1),
        })
        .unwrap_err();
    assert!(err.message.contains("between 0.0 and 1.0"));
}

#[test]
fn test_set_parameter_feature_details_font_size_invalid_type_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "feature_details_font_size".to_string(),
            value: serde_json::json!("small"),
        })
        .unwrap_err();
    assert!(
        err.message
            .contains("feature_details_font_size requires a number")
    );
}

#[test]
fn test_set_parameter_feature_details_font_size_out_of_range_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "feature_details_font_size".to_string(),
            value: serde_json::json!(3.0),
        })
        .unwrap_err();
    assert!(
        err.message
            .contains("feature_details_font_size must be between 8.0 and 24.0")
    );
}

#[test]
fn test_set_parameter_linear_external_feature_label_font_size_out_of_range_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "linear_external_feature_label_font_size".to_string(),
            value: serde_json::json!(30.0),
        })
        .unwrap_err();
    assert!(
        err.message
            .contains("linear_external_feature_label_font_size must be between 8.0 and 24.0")
    );
}

#[test]
fn test_set_parameter_linear_external_feature_label_background_opacity_out_of_range_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "linear_external_feature_label_background_opacity".to_string(),
            value: serde_json::json!(1.5),
        })
        .unwrap_err();
    assert!(
        err.message.contains(
            "linear_external_feature_label_background_opacity must be between 0.0 and 1.0"
        )
    );
}

#[test]
fn test_set_parameter_reverse_strand_visual_opacity_out_of_range_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "reverse_strand_visual_opacity".to_string(),
            value: serde_json::json!(0.05),
        })
        .unwrap_err();
    assert!(
        err.message
            .contains("reverse_strand_visual_opacity must be between 0.2 and 1.0")
    );
}

#[test]
fn test_set_parameter_linear_sequence_helical_phase_offset_out_of_range_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "linear_sequence_helical_phase_offset_bp".to_string(),
            value: serde_json::json!(10),
        })
        .unwrap_err();
    assert!(
        err.message
            .contains("linear_sequence_helical_phase_offset_bp must be between 0 and 9")
    );
}

#[test]
fn test_set_parameter_regulatory_feature_max_view_span_bp_invalid_type_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "regulatory_feature_max_view_span_bp".to_string(),
            value: serde_json::json!("wide"),
        })
        .unwrap_err();
    assert!(
        err.message
            .contains("regulatory_feature_max_view_span_bp requires a non-negative integer")
    );
}

#[test]
fn test_set_parameter_gc_content_bin_size_bp_invalid_type_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "gc_content_bin_size_bp".to_string(),
            value: serde_json::json!("fine"),
        })
        .unwrap_err();
    assert!(
        err.message
            .contains("gc_content_bin_size_bp requires a positive integer")
    );
}

#[test]
fn test_set_parameter_gc_content_bin_size_bp_zero_fails() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::SetParameter {
            name: "gc_content_bin_size_bp".to_string(),
            value: serde_json::json!(0),
        })
        .unwrap_err();
    assert!(err.message.contains("gc_content_bin_size_bp must be >= 1"));
}

#[test]
fn test_containers_created_on_load_and_digest() {
    let mut engine = GentleEngine::new();
    let load = engine
        .apply(Operation::LoadFile {
            path: "test_files/pGEX_3X.fa".to_string(),
            as_id: Some("pgex".to_string()),
        })
        .unwrap();
    let seq_id = load.created_seq_ids.first().unwrap().clone();
    let latest = engine
        .state()
        .container_state
        .seq_to_latest_container
        .get(&seq_id)
        .cloned();
    assert!(latest.is_some());
    let digest = engine
        .apply(Operation::Digest {
            input: "pgex".to_string(),
            enzymes: vec!["BamHI".to_string(), "EcoRI".to_string()],
            output_prefix: Some("frag".to_string()),
        })
        .unwrap();
    assert!(digest.created_seq_ids.len() >= 2);
    let container_id = engine
        .state()
        .container_state
        .seq_to_latest_container
        .get(digest.created_seq_ids.first().unwrap())
        .unwrap();
    let container = engine
        .state()
        .container_state
        .containers
        .get(container_id)
        .unwrap();
    assert!(matches!(container.kind, ContainerKind::Pool));
    assert_eq!(container.members.len(), digest.created_seq_ids.len());
}

#[test]
fn test_select_candidate_creates_selection_container_kind() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("frag".to_string(), seq(&"ATGC".repeat(40)));
    let mut engine = GentleEngine::from_state(state);
    let res = engine
        .apply(Operation::SelectCandidate {
            input: "frag".to_string(),
            criterion: "manual".to_string(),
            output_id: Some("picked".to_string()),
        })
        .unwrap();
    let picked = res.created_seq_ids.first().unwrap();
    let cid = engine
        .state()
        .container_state
        .seq_to_latest_container
        .get(picked)
        .unwrap();
    let c = engine.state().container_state.containers.get(cid).unwrap();
    assert!(matches!(c.kind, ContainerKind::Selection));
}

#[test]
fn test_container_operations_map_to_core_ops() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("a".to_string(), seq(&"ATGC".repeat(40)));
    state
        .sequences
        .insert("b".to_string(), seq(&"TTAA".repeat(40)));
    let mut engine = GentleEngine::from_state(state);
    let merge = engine
        .apply(Operation::MergeContainersById {
            container_ids: vec![
                engine
                    .state()
                    .container_state
                    .seq_to_latest_container
                    .get("a")
                    .unwrap()
                    .clone(),
                engine
                    .state()
                    .container_state
                    .seq_to_latest_container
                    .get("b")
                    .unwrap()
                    .clone(),
            ],
            output_prefix: Some("pool".to_string()),
        })
        .unwrap();
    assert_eq!(merge.created_seq_ids.len(), 2);

    let pool_container = engine
        .state()
        .container_state
        .seq_to_latest_container
        .get(merge.created_seq_ids.first().unwrap())
        .unwrap()
        .clone();
    let lig = engine
        .apply(Operation::LigationContainer {
            container_id: pool_container.clone(),
            circularize_if_possible: false,
            output_id: None,
            protocol: LigationProtocol::Blunt,
            output_prefix: Some("lig".to_string()),
            unique: None,
        })
        .unwrap();
    assert!(!lig.created_seq_ids.is_empty());

    let filtered = engine
        .apply(Operation::FilterContainerByMolecularWeight {
            container_id: pool_container,
            min_bp: 120,
            max_bp: 400,
            error: 0.0,
            unique: false,
            output_prefix: Some("mw".to_string()),
        })
        .unwrap();
    assert!(!filtered.created_seq_ids.is_empty());
}

#[test]
fn test_digest_container_and_state_summary_include_containers() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::LoadFile {
            path: "test_files/pGEX_3X.fa".to_string(),
            as_id: Some("pgex".to_string()),
        })
        .unwrap();
    let cid = engine
        .state()
        .container_state
        .seq_to_latest_container
        .get("pgex")
        .unwrap()
        .clone();
    let res = engine
        .apply(Operation::DigestContainer {
            container_id: cid,
            enzymes: vec!["BamHI".to_string(), "EcoRI".to_string()],
            output_prefix: Some("dig".to_string()),
        })
        .unwrap();
    assert!(!res.created_seq_ids.is_empty());
    let summary = engine.summarize_state();
    assert!(summary.container_count > 0);
    assert!(!summary.containers.is_empty());
}

#[test]
fn test_prepare_genome_and_extract_region_operations() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta_gz = root.join("toy.fa.gz");
    let ann_gz = root.join("toy.gtf.gz");
    write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\nACGT\n");
    write_gzip(
        &ann_gz,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
    );
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "description": "toy genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        file_url(&fasta_gz),
        file_url(&ann_gz),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();

    let mut engine = GentleEngine::new();
    let prep = engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    assert!(
        prep.messages
            .iter()
            .any(|m| m.contains("Prepared genome 'ToyGenome'"))
    );
    assert!(
        prep.messages
            .iter()
            .any(|m| m.contains("cached_sequence=1 contigs")),
        "messages were: {:?}",
        prep.messages
    );
    let catalog_path_str = catalog_path.to_string_lossy().to_string();
    let catalog_names = GentleEngine::list_reference_genomes(Some(&catalog_path_str)).unwrap();
    assert!(catalog_names.contains(&"ToyGenome".to_string()));
    let prepared =
        GentleEngine::is_reference_genome_prepared(Some(&catalog_path_str), "ToyGenome", None)
            .unwrap();
    assert!(prepared);
    let listed_genes =
        GentleEngine::list_reference_genome_genes(Some(&catalog_path_str), "ToyGenome", None)
            .unwrap();
    assert_eq!(listed_genes.len(), 1);

    let extract = engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 3,
            end_1based: 10,
            output_id: Some("toy_slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
        })
        .unwrap();
    assert_eq!(extract.created_seq_ids, vec!["toy_slice".to_string()]);
    let loaded = engine.state().sequences.get("toy_slice").unwrap();
    assert_eq!(loaded.get_forward_string(), "GTACGTAC");

    let extract_gene = engine
        .apply(Operation::ExtractGenomeGene {
            genome_id: "ToyGenome".to_string(),
            gene_query: "MYGENE".to_string(),
            occurrence: None,
            output_id: Some("toy_gene".to_string()),
            extract_mode: None,
            promoter_upstream_bp: None,
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .unwrap();
    assert!(
        extract_gene
            .created_seq_ids
            .iter()
            .any(|id| id == "toy_gene")
    );
    let loaded_gene = engine.state().sequences.get("toy_gene").unwrap();
    assert_eq!(loaded_gene.get_forward_string(), "ACGTACGTACGT");
    let provenance = engine
        .state()
        .metadata
        .get("provenance")
        .and_then(|v| v.as_object())
        .expect("provenance metadata object");
    let extractions = provenance
        .get("genome_extractions")
        .and_then(|v| v.as_array())
        .expect("genome_extractions array");
    assert_eq!(extractions.len(), 2);
    assert!(extractions.iter().any(|entry| {
        entry
            .get("seq_id")
            .and_then(|v| v.as_str())
            .map(|v| v == "toy_slice")
            .unwrap_or(false)
            && entry
                .get("operation")
                .and_then(|v| v.as_str())
                .map(|v| v == "ExtractGenomeRegion")
                .unwrap_or(false)
    }));
    assert!(extractions.iter().any(|entry| {
        entry
            .get("seq_id")
            .and_then(|v| v.as_str())
            .map(|v| v == "toy_gene")
            .unwrap_or(false)
            && entry
                .get("operation")
                .and_then(|v| v.as_str())
                .map(|v| v == "ExtractGenomeGene")
                .unwrap_or(false)
            && entry
                .get("sequence_sha1")
                .and_then(|v| v.as_str())
                .map(|v| !v.is_empty())
                .unwrap_or(false)
    }));
}

#[test]
fn test_extract_genome_gene_reports_alias_guidance_for_contig_mismatch() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    fs::write(&fasta, ">NC_000017.11\nACGTACGTACGT\n").unwrap();
    fs::write(
        &ann,
        "17\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
    )
    .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();

    let err = engine
        .apply(Operation::ExtractGenomeGene {
            genome_id: "ToyGenome".to_string(),
            gene_query: "MYGENE".to_string(),
            occurrence: None,
            output_id: Some("toy_gene".to_string()),
            extract_mode: None,
            promoter_upstream_bp: None,
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
        })
        .expect_err("contig mismatch should fail");
    assert!(
        err.message
            .contains("Could not load gene extraction interval 17:1-12 from 'ToyGenome'")
    );
    assert!(err.message.contains("Tried aliases:"));
    assert!(err.message.contains("Available contigs"));
    assert!(err.message.contains("Suggested matching contigs"));
    assert!(err.message.contains("NC_000017.11"));
}

#[test]
fn test_extract_genome_gene_attaches_transcript_features_from_annotation() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGT\n").unwrap();
    fs::write(
            &ann,
            concat!(
                "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
                "chr1\tsrc\ttranscript\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\texon\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"1\";\n",
                "chr1\tsrc\tCDS\t2\t4\t.\t+\t0\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\texon\t9\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"2\";\n",
                "chr1\tsrc\tCDS\t9\t11\t.\t+\t2\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
            ),
        )
        .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    let extract_gene = engine
        .apply(Operation::ExtractGenomeGene {
            genome_id: "ToyGenome".to_string(),
            gene_query: "MYGENE".to_string(),
            occurrence: None,
            output_id: Some("toy_gene".to_string()),
            extract_mode: None,
            promoter_upstream_bp: None,
            annotation_scope: Some(GenomeAnnotationScope::Full),
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
        })
        .unwrap();
    assert!(
        extract_gene
            .created_seq_ids
            .iter()
            .any(|id| id == "toy_gene__exons")
    );
    assert!(extract_gene.messages.iter().any(|m| m.contains("Attached")));
    let loaded_gene = engine.state().sequences.get("toy_gene").unwrap();
    assert_eq!(loaded_gene.name().as_deref(), Some("toy_gene"));
    let gene_feature = loaded_gene
        .features()
        .iter()
        .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
        .expect("expected extracted gene feature");
    assert!(
        gene_feature
            .qualifier_values("gene_id")
            .any(|value| value == "GENE1")
    );
    assert!(
        gene_feature
            .qualifier_values("gentle_context_layer")
            .any(|value| value.eq_ignore_ascii_case("contextual_gene"))
    );
    assert!(
        gene_feature
            .qualifier_values("gentle_generated")
            .any(|value| value.eq_ignore_ascii_case("genome_annotation_projection"))
    );
    let tx_feature = loaded_gene
        .features()
        .iter()
        .find(|feature| {
            feature.kind.to_string().eq_ignore_ascii_case("mRNA")
                && feature
                    .qualifier_values("transcript_id")
                    .any(|value| value == "TX1")
        })
        .expect("expected extracted transcript feature");
    assert!(
        tx_feature
            .qualifier_values("gentle_context_layer")
            .any(|value| value.eq_ignore_ascii_case("contextual_transcript"))
    );
    assert!(
        tx_feature
            .qualifier_values("gentle_generated")
            .any(|value| value.eq_ignore_ascii_case("genome_annotation_projection"))
    );
    let mut exons = vec![];
    collect_location_ranges_usize(&tx_feature.location, &mut exons);
    exons.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    assert_eq!(exons, vec![(0, 4), (8, 12)]);
    let projected_exon = loaded_gene
        .features()
        .iter()
        .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("exon"))
        .expect("expected projected exon feature");
    assert!(
        projected_exon
            .qualifier_values("gentle_context_layer")
            .any(|value| value.eq_ignore_ascii_case("contextual_transcript"))
    );
    assert!(
        loaded_gene
            .features()
            .iter()
            .any(|feature| feature.kind.to_string().eq_ignore_ascii_case("exon"))
    );
    assert!(
        loaded_gene
            .features()
            .iter()
            .any(|feature| feature.kind.to_string().eq_ignore_ascii_case("cds"))
    );
    let exon_concat = engine
        .state()
        .sequences
        .get("toy_gene__exons")
        .expect("expected exon-concatenated synthetic sequence");
    assert_eq!(
        exon_concat.get_forward_string(),
        format!("ACGT{}ACGT", "N".repeat(24))
    );
    assert_eq!(
        exon_concat
            .features()
            .iter()
            .filter(|feature| feature.kind.to_string().eq_ignore_ascii_case("exon"))
            .count(),
        2
    );
}

#[test]
fn test_extract_genome_gene_exon_concat_respects_negative_strand_orientation() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nAAAACCCCGGGGTTTT\n").unwrap();
    fs::write(
        &ann,
        concat!(
            "chr1\tsrc\tgene\t1\t12\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\";\n",
            "chr1\tsrc\ttranscript\t1\t12\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\";\n",
            "chr1\tsrc\texon\t1\t4\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\"; exon_number \"1\";\n",
            "chr1\tsrc\texon\t9\t12\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\"; exon_number \"2\";\n",
        ),
    )
    .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    engine
        .apply(Operation::ExtractGenomeGene {
            genome_id: "ToyGenome".to_string(),
            gene_query: "NEG1".to_string(),
            occurrence: None,
            output_id: Some("toy_gene_neg".to_string()),
            extract_mode: None,
            promoter_upstream_bp: None,
            annotation_scope: Some(GenomeAnnotationScope::Core),
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
        })
        .unwrap();
    let loaded_gene = engine
        .state()
        .sequences
        .get("toy_gene_neg")
        .expect("expected extracted negative-strand gene sequence");
    let gene_feature = loaded_gene
        .features()
        .iter()
        .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
        .expect("expected projected gene feature");
    assert!(crate::feature_location::feature_is_reverse(gene_feature));
    assert!(matches!(
        gene_feature.location,
        gb_io::seq::Location::Complement(_)
    ));
    let tx_feature = loaded_gene
        .features()
        .iter()
        .find(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("expected projected transcript feature");
    assert!(crate::feature_location::feature_is_reverse(tx_feature));
    let exon_concat = engine
        .state()
        .sequences
        .get("toy_gene_neg__exons")
        .expect("expected negative-strand exon-concatenated synthetic sequence");
    assert_eq!(
        exon_concat.get_forward_string(),
        format!("CCCC{}TTTT", " ".repeat(24))
    );
    let exon_features: Vec<_> = exon_concat
        .features()
        .iter()
        .filter(|feature| feature.kind.to_string().eq_ignore_ascii_case("exon"))
        .collect();
    assert_eq!(exon_features.len(), 2);
    let first_genomic_start = exon_features[0]
        .qualifier_values("genomic_start_1based")
        .next()
        .unwrap_or_default();
    assert_eq!(first_genomic_start, "9");
}

#[test]
fn test_extract_genome_gene_coding_with_promoter_extends_plus_strand_from_first_cds_base() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGT\n").unwrap();
    fs::write(
        &ann,
        concat!(
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
            "chr1\tsrc\ttranscript\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
            "chr1\tsrc\texon\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"1\";\n",
            "chr1\tsrc\tCDS\t4\t6\t.\t+\t0\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
            "chr1\tsrc\texon\t9\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"2\";\n",
            "chr1\tsrc\tCDS\t9\t11\t.\t+\t2\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
        ),
    )
    .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    engine
        .apply(Operation::ExtractGenomeGene {
            genome_id: "ToyGenome".to_string(),
            gene_query: "MYGENE".to_string(),
            occurrence: None,
            output_id: Some("toy_gene_promoter".to_string()),
            extract_mode: Some(GenomeGeneExtractMode::CodingWithPromoter),
            promoter_upstream_bp: Some(2),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
        })
        .unwrap();

    let loaded_gene = engine.state().sequences.get("toy_gene_promoter").unwrap();
    assert_eq!(loaded_gene.get_forward_string(), "CGTACGTACG");
    let provenance = engine
        .state()
        .metadata
        .get("provenance")
        .and_then(|value| value.as_object())
        .and_then(|object| object.get("genome_extractions"))
        .and_then(|value| value.as_array())
        .expect("genome extraction provenance");
    let entry = provenance
        .iter()
        .find(|entry| {
            entry
                .get("seq_id")
                .and_then(|value| value.as_str())
                .map(|value| value == "toy_gene_promoter")
                .unwrap_or(false)
        })
        .expect("provenance entry for promoter extract");
    assert_eq!(
        entry
            .get("gene_extract_mode")
            .and_then(|value| value.as_str()),
        Some("coding_with_promoter")
    );
    assert_eq!(
        entry
            .get("promoter_upstream_bp")
            .and_then(|value| value.as_u64()),
        Some(2)
    );
}

#[test]
fn test_extract_genome_gene_coding_with_promoter_extends_minus_strand_on_five_prime_side() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nAAAACCCCGGGGTTTT\n").unwrap();
    fs::write(
        &ann,
        concat!(
            "chr1\tsrc\tgene\t1\t12\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\";\n",
            "chr1\tsrc\ttranscript\t1\t12\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\";\n",
            "chr1\tsrc\texon\t1\t4\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\"; exon_number \"1\";\n",
            "chr1\tsrc\tCDS\t2\t4\t.\t-\t0\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\";\n",
            "chr1\tsrc\texon\t9\t12\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\"; exon_number \"2\";\n",
            "chr1\tsrc\tCDS\t9\t11\t.\t-\t2\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\";\n",
        ),
    )
    .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    engine
        .apply(Operation::ExtractGenomeGene {
            genome_id: "ToyGenome".to_string(),
            gene_query: "NEG1".to_string(),
            occurrence: None,
            output_id: Some("toy_gene_neg_promoter".to_string()),
            extract_mode: Some(GenomeGeneExtractMode::CodingWithPromoter),
            promoter_upstream_bp: Some(2),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
        })
        .unwrap();

    let loaded_gene = engine
        .state()
        .sequences
        .get("toy_gene_neg_promoter")
        .unwrap();
    assert_eq!(loaded_gene.get_forward_string(), "AAACCCCGGGGT");
}

#[test]
fn test_extract_genome_gene_include_annotation_false_disables_projection() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGT\n").unwrap();
    fs::write(
            &ann,
            concat!(
                "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
                "chr1\tsrc\ttranscript\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\texon\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"1\";\n",
                "chr1\tsrc\tCDS\t2\t4\t.\t+\t0\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\texon\t9\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"2\";\n",
                "chr1\tsrc\tCDS\t9\t11\t.\t+\t2\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
            ),
        )
        .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    let extract_gene = engine
        .apply(Operation::ExtractGenomeGene {
            genome_id: "ToyGenome".to_string(),
            gene_query: "MYGENE".to_string(),
            occurrence: None,
            output_id: Some("toy_gene_no_annotation".to_string()),
            extract_mode: None,
            promoter_upstream_bp: None,
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: Some(false),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
        })
        .unwrap();

    let loaded_gene = engine
        .state()
        .sequences
        .get("toy_gene_no_annotation")
        .unwrap();
    assert!(
        loaded_gene.features().is_empty(),
        "annotation-disabled extraction should not attach features"
    );
    let projection = extract_gene
        .genome_annotation_projection
        .as_ref()
        .expect("genome_annotation_projection telemetry");
    assert_eq!(projection.requested_scope, "none");
    assert_eq!(projection.effective_scope, "none");
    assert_eq!(projection.attached_feature_count, 0);
    assert!(!projection.fallback_applied);
}

#[test]
fn test_extract_genome_gene_annotation_cap_falls_back_to_core() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGT\n").unwrap();
    fs::write(
            &ann,
            concat!(
                "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
                "chr1\tsrc\ttranscript\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\texon\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"1\";\n",
                "chr1\tsrc\tCDS\t2\t4\t.\t+\t0\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\texon\t9\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"2\";\n",
                "chr1\tsrc\tCDS\t9\t11\t.\t+\t2\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
            ),
        )
        .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    let extract_gene = engine
        .apply(Operation::ExtractGenomeGene {
            genome_id: "ToyGenome".to_string(),
            gene_query: "MYGENE".to_string(),
            occurrence: None,
            output_id: Some("toy_gene_cap_fallback".to_string()),
            extract_mode: None,
            promoter_upstream_bp: None,
            annotation_scope: Some(GenomeAnnotationScope::Full),
            max_annotation_features: Some(2),
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
        })
        .unwrap();

    let loaded_gene = engine
        .state()
        .sequences
        .get("toy_gene_cap_fallback")
        .unwrap();
    assert!(
        loaded_gene
            .features()
            .iter()
            .any(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
    );
    assert!(
        loaded_gene
            .features()
            .iter()
            .any(|feature| feature.kind.to_string().eq_ignore_ascii_case("mrna"))
    );
    assert!(
        !loaded_gene
            .features()
            .iter()
            .any(|feature| feature.kind.to_string().eq_ignore_ascii_case("exon"))
    );
    assert!(
        !loaded_gene
            .features()
            .iter()
            .any(|feature| feature.kind.to_string().eq_ignore_ascii_case("cds"))
    );
    let projection = extract_gene
        .genome_annotation_projection
        .as_ref()
        .expect("genome_annotation_projection telemetry");
    assert_eq!(projection.requested_scope, "full");
    assert_eq!(projection.effective_scope, "core");
    assert_eq!(projection.max_features_cap, Some(2));
    assert_eq!(projection.attached_feature_count, 2);
    assert!(projection.fallback_applied);
    assert!(
        projection
            .fallback_reason
            .as_deref()
            .unwrap_or_default()
            .contains("fell back to core projection")
    );
}

#[test]
fn test_extract_genome_region_include_annotation_attaches_features_and_sets_name() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGT\n").unwrap();
    fs::write(
            &ann,
            concat!(
                "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
                "chr1\tsrc\ttranscript\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\texon\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"1\";\n",
                "chr1\tsrc\texon\t9\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"2\";\n",
            ),
        )
        .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    let catalog_path_str = catalog_path.to_string_lossy().to_string();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();

    let with_annotation = engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 1,
            end_1based: 12,
            output_id: Some("toy_slice_ann".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: Some(true),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
        })
        .unwrap();
    assert_eq!(
        with_annotation.created_seq_ids,
        vec!["toy_slice_ann".to_string()]
    );
    assert!(
        with_annotation
            .messages
            .iter()
            .any(|m| m.contains("Attached"))
    );
    let telemetry = with_annotation
        .genome_annotation_projection
        .as_ref()
        .expect("annotation telemetry");
    assert_eq!(telemetry.requested_scope, "core");
    assert_eq!(telemetry.effective_scope, "core");
    assert!(telemetry.attached_feature_count >= 2);
    let annotated = engine.state().sequences.get("toy_slice_ann").unwrap();
    assert_eq!(annotated.name().as_deref(), Some("toy_slice_ann"));
    assert!(annotated.features().iter().any(|feature| {
        feature.kind.to_string().eq_ignore_ascii_case("gene")
            && feature
                .qualifier_values("gene_id")
                .any(|value| value == "GENE1")
    }));
    assert!(annotated.features().iter().any(|feature| {
        feature.kind.to_string().eq_ignore_ascii_case("mRNA")
            && feature
                .qualifier_values("transcript_id")
                .any(|value| value == "TX1")
    }));

    let default_result = engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 1,
            end_1based: 12,
            output_id: Some("toy_slice_default".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
        })
        .unwrap();
    let default_telemetry = default_result
        .genome_annotation_projection
        .as_ref()
        .expect("default annotation telemetry");
    assert_eq!(default_telemetry.requested_scope, "core");
    assert_eq!(default_telemetry.effective_scope, "core");
    let default_slice = engine.state().sequences.get("toy_slice_default").unwrap();
    assert_eq!(default_slice.name().as_deref(), Some("toy_slice_default"));
    assert!(default_slice.features().iter().any(|feature| {
        feature.kind.to_string().eq_ignore_ascii_case("gene")
            && feature
                .qualifier_values("gene_id")
                .any(|value| value == "GENE1")
    }));

    let plain_result = engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 1,
            end_1based: 12,
            output_id: Some("toy_slice_plain".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: Some(false),
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .unwrap();
    let plain_telemetry = plain_result
        .genome_annotation_projection
        .as_ref()
        .expect("plain annotation telemetry");
    assert_eq!(plain_telemetry.requested_scope, "none");
    assert_eq!(plain_telemetry.effective_scope, "none");
    assert_eq!(plain_telemetry.attached_feature_count, 0);
    let plain = engine.state().sequences.get("toy_slice_plain").unwrap();
    assert_eq!(plain.name().as_deref(), Some("toy_slice_plain"));
    assert!(plain.features().is_empty());
}

#[test]
fn test_extract_genome_region_full_scope_feature_cap_falls_back_to_core() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGT\n").unwrap();
    fs::write(
            &ann,
            concat!(
                "chr1\tsrc\tgene\t1\t20\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
                "chr1\tsrc\ttranscript\t1\t20\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\texon\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"1\";\n",
                "chr1\tsrc\texon\t9\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\"; exon_number \"2\";\n",
                "chr1\tsrc\tCDS\t2\t4\t.\t+\t0\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\tCDS\t9\t11\t.\t+\t0\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX1\";\n",
            ),
        )
        .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();

    let result = engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 1,
            end_1based: 20,
            output_id: Some("toy_slice_capped".to_string()),
            annotation_scope: Some(GenomeAnnotationScope::Full),
            max_annotation_features: Some(2),
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .unwrap();
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.contains("fell back to core"))
    );
    let telemetry = result
        .genome_annotation_projection
        .as_ref()
        .expect("annotation telemetry");
    assert_eq!(telemetry.requested_scope, "full");
    assert_eq!(telemetry.effective_scope, "core");
    assert!(telemetry.fallback_applied);
    assert!(telemetry.candidate_feature_count > telemetry.attached_feature_count);
    assert!(telemetry.exons_attached == 0);
    assert!(telemetry.cds_attached == 0);

    let seq = engine.state().sequences.get("toy_slice_capped").unwrap();
    assert!(
        seq.features()
            .iter()
            .any(|f| f.kind.to_string().eq_ignore_ascii_case("gene"))
    );
    assert!(
        seq.features()
            .iter()
            .any(|f| f.kind.to_string().eq_ignore_ascii_case("mRNA"))
    );
    assert!(
        !seq.features()
            .iter()
            .any(|f| f.kind.to_string().eq_ignore_ascii_case("exon"))
    );
    assert!(
        !seq.features()
            .iter()
            .any(|f| f.kind.to_string().eq_ignore_ascii_case("CDS"))
    );
}

#[test]
fn test_prepare_genome_operation_supports_timeout_seconds() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGT\nACGT\n").unwrap();
    fs::write(
        &ann,
        "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
    )
    .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();

    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: Some(0),
        })
        .unwrap_err();
    assert!(err.message.to_ascii_lowercase().contains("timed out"));
}

#[test]
fn test_extend_genome_anchor_plus_strand_adds_lineage_and_provenance() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta_gz = root.join("toy.fa.gz");
    let ann_gz = root.join("toy.gtf.gz");
    write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\nACGT\n");
    write_gzip(
        &ann_gz,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
    );
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "description": "toy genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        file_url(&fasta_gz),
        file_url(&ann_gz),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 3,
            end_1based: 10,
            output_id: Some("toy_slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
        })
        .unwrap();

    let extended = engine
        .apply(Operation::ExtendGenomeAnchor {
            seq_id: "toy_slice".to_string(),
            side: GenomeAnchorSide::FivePrime,
            length_bp: 2,
            output_id: Some("toy_slice_ext5".to_string()),
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
            prepared_genome_id: None,
        })
        .unwrap();
    assert_eq!(extended.created_seq_ids, vec!["toy_slice_ext5".to_string()]);
    let extended_seq = engine
        .state()
        .sequences
        .get("toy_slice_ext5")
        .expect("extended sequence should exist");
    assert_eq!(extended_seq.get_forward_string(), "ACGTACGTAC");
    assert_eq!(extended_seq.name().as_deref(), Some("toy_slice_ext5"));

    let lineage = &engine.state().lineage;
    let parent = lineage
        .seq_to_node
        .get("toy_slice")
        .expect("parent lineage node should exist");
    let child = lineage
        .seq_to_node
        .get("toy_slice_ext5")
        .expect("child lineage node should exist");
    assert!(
        lineage
            .edges
            .iter()
            .any(|edge| edge.from_node_id == *parent && edge.to_node_id == *child)
    );

    let provenance = engine
        .state()
        .metadata
        .get(PROVENANCE_METADATA_KEY)
        .and_then(|v| v.as_object())
        .expect("provenance metadata object");
    let extractions = provenance
        .get(GENOME_EXTRACTIONS_METADATA_KEY)
        .and_then(|v| v.as_array())
        .expect("genome_extractions array");
    let entry = extractions
        .iter()
        .find(|entry| {
            entry
                .get("seq_id")
                .and_then(|v| v.as_str())
                .map(|v| v == "toy_slice_ext5")
                .unwrap_or(false)
        })
        .expect("extended provenance entry");
    assert_eq!(
        entry.get("operation").and_then(|v| v.as_str()),
        Some("ExtendGenomeAnchor")
    );
    assert_eq!(entry.get("start_1based").and_then(|v| v.as_u64()), Some(1));
    assert_eq!(entry.get("end_1based").and_then(|v| v.as_u64()), Some(10));
    assert_eq!(
        entry.get("anchor_strand").and_then(|v| v.as_str()),
        Some("+")
    );
}

#[test]
fn test_extend_genome_anchor_reverse_strand_respects_5prime_and_3prime_physical_direction() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta_gz = root.join("toy.fa.gz");
    let ann_gz = root.join("toy.gtf.gz");
    write_gzip(&fasta_gz, ">chr1\nACGTTGCAATGCCGTA\n");
    write_gzip(
        &ann_gz,
        "chr1\tsrc\tgene\t1\t16\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
    );
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "description": "toy genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        file_url(&fasta_gz),
        file_url(&ann_gz),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("rev_anchor".to_string(), seq("ATTGCA"));
    state.metadata.insert(
        PROVENANCE_METADATA_KEY.to_string(),
        serde_json::json!({
            GENOME_EXTRACTIONS_METADATA_KEY: [
                {
                    "seq_id": "rev_anchor",
                    "recorded_at_unix_ms": 1,
                    "operation": "LoadFileGenBankRegion",
                    "genome_id": "ToyGenome",
                    "catalog_path": catalog_path_str,
                    "cache_dir": null,
                    "chromosome": "chr1",
                    "start_1based": 5,
                    "end_1based": 10,
                    "gene_query": null,
                    "occurrence": null,
                    "gene_id": null,
                    "gene_name": null,
                    "strand": null,
                    "anchor_strand": "-",
                    "sequence_source_type": "synthetic",
                    "annotation_source_type": "synthetic",
                    "sequence_source": "synthetic",
                    "annotation_source": "synthetic",
                    "sequence_sha1": null,
                    "annotation_sha1": null
                }
            ]
        }),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();

    let ext5 = engine
        .apply(Operation::ExtendGenomeAnchor {
            seq_id: "rev_anchor".to_string(),
            side: GenomeAnchorSide::FivePrime,
            length_bp: 3,
            output_id: Some("rev_ext5".to_string()),
            catalog_path: None,
            cache_dir: None,
            prepared_genome_id: None,
        })
        .unwrap();
    assert_eq!(ext5.created_seq_ids, vec!["rev_ext5".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("rev_ext5")
            .expect("rev_ext5 sequence")
            .get_forward_string(),
        "GGCATTGCA"
    );

    let ext3 = engine
        .apply(Operation::ExtendGenomeAnchor {
            seq_id: "rev_anchor".to_string(),
            side: GenomeAnchorSide::ThreePrime,
            length_bp: 2,
            output_id: Some("rev_ext3".to_string()),
            catalog_path: None,
            cache_dir: None,
            prepared_genome_id: None,
        })
        .unwrap();
    assert_eq!(ext3.created_seq_ids, vec!["rev_ext3".to_string()]);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("rev_ext3")
            .expect("rev_ext3 sequence")
            .get_forward_string(),
        "ATTGCAAC"
    );

    let provenance = engine
        .state()
        .metadata
        .get(PROVENANCE_METADATA_KEY)
        .and_then(|v| v.as_object())
        .expect("provenance metadata object");
    let extractions = provenance
        .get(GENOME_EXTRACTIONS_METADATA_KEY)
        .and_then(|v| v.as_array())
        .expect("genome_extractions array");

    let ext5_entry = extractions
        .iter()
        .find(|entry| {
            entry
                .get("seq_id")
                .and_then(|v| v.as_str())
                .map(|v| v == "rev_ext5")
                .unwrap_or(false)
        })
        .expect("rev_ext5 provenance");
    assert_eq!(
        ext5_entry.get("start_1based").and_then(|v| v.as_u64()),
        Some(5)
    );
    assert_eq!(
        ext5_entry.get("end_1based").and_then(|v| v.as_u64()),
        Some(13)
    );
    assert_eq!(
        ext5_entry.get("anchor_strand").and_then(|v| v.as_str()),
        Some("-")
    );

    let ext3_entry = extractions
        .iter()
        .find(|entry| {
            entry
                .get("seq_id")
                .and_then(|v| v.as_str())
                .map(|v| v == "rev_ext3")
                .unwrap_or(false)
        })
        .expect("rev_ext3 provenance");
    assert_eq!(
        ext3_entry.get("start_1based").and_then(|v| v.as_u64()),
        Some(3)
    );
    assert_eq!(
        ext3_entry.get("end_1based").and_then(|v| v.as_u64()),
        Some(10)
    );
    assert_eq!(
        ext3_entry.get("anchor_strand").and_then(|v| v.as_str()),
        Some("-")
    );
}

#[test]
fn test_extend_genome_anchor_uses_compatible_prepared_assembly_fallback() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("grch38.fa");
    let ann = root.join("grch38.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGT\n").unwrap();
    fs::write(
        &ann,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
    )
    .unwrap();
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "Human GRCh38 Ensembl 116": {{
    "ncbi_taxonomy_id": 9606,
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }},
  "Human GRCh38 NCBI RefSeq GCF_000001405.40": {{
    "ncbi_taxonomy_id": 9606,
    "ncbi_assembly_accession": "GCF_000001405.40",
    "ncbi_assembly_name": "GRCh38.p14",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        ann.display(),
        cache_dir.display(),
        fasta.display(),
        ann.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "Human GRCh38 Ensembl 116".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "GRCh38.p14".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 3,
            end_1based: 10,
            output_id: Some("alias_slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
        })
        .unwrap();

    let extended = engine
        .apply(Operation::ExtendGenomeAnchor {
            seq_id: "alias_slice".to_string(),
            side: GenomeAnchorSide::FivePrime,
            length_bp: 2,
            output_id: Some("alias_slice_ext5".to_string()),
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
            prepared_genome_id: None,
        })
        .unwrap();
    assert_eq!(
        extended.created_seq_ids,
        vec!["alias_slice_ext5".to_string()]
    );
    assert!(
        extended
            .warnings
            .iter()
            .any(|w| w.contains("compatible prepared genome"))
    );
    let extended_seq = engine
        .state()
        .sequences
        .get("alias_slice_ext5")
        .expect("extended sequence should exist");
    assert_eq!(extended_seq.get_forward_string(), "ACGTACGTAC");
}

#[test]
fn test_extend_genome_anchor_rejects_zero_length() {
    let mut engine = GentleEngine::new();
    let err = engine
        .apply(Operation::ExtendGenomeAnchor {
            seq_id: "missing".to_string(),
            side: GenomeAnchorSide::FivePrime,
            length_bp: 0,
            output_id: None,
            catalog_path: None,
            cache_dir: None,
            prepared_genome_id: None,
        })
        .unwrap_err();
    assert!(matches!(err.code, ErrorCode::InvalidInput));
    assert!(err.message.contains("length_bp >= 1"));
}

#[test]
fn test_extend_genome_anchor_requires_genome_anchor_provenance() {
    let mut state = ProjectState::default();
    state.sequences.insert("plain".to_string(), seq("ACGTACGT"));
    let mut engine = GentleEngine::from_state(state);
    let err = engine
        .apply(Operation::ExtendGenomeAnchor {
            seq_id: "plain".to_string(),
            side: GenomeAnchorSide::ThreePrime,
            length_bp: 10,
            output_id: None,
            catalog_path: None,
            cache_dir: None,
            prepared_genome_id: None,
        })
        .unwrap_err();
    assert!(matches!(err.code, ErrorCode::NotFound));
    assert!(err.message.contains("no genome anchor provenance"));
}

#[test]
fn test_extend_genome_anchor_strict_verification_rejects_unverified_anchor() {
    let mut state = ProjectState::default();
    state.sequences.insert("anch".to_string(), seq("ACGTACGT"));
    state.metadata.insert(
        PROVENANCE_METADATA_KEY.to_string(),
        serde_json::json!({
            GENOME_EXTRACTIONS_METADATA_KEY: [
                {
                    "seq_id": "anch",
                    "recorded_at_unix_ms": 1,
                    "operation": "LoadFileGenBankRegion",
                    "genome_id": "ToyGenome",
                    "catalog_path": "assets/genomes.json",
                    "cache_dir": null,
                    "chromosome": "chr1",
                    "start_1based": 1,
                    "end_1based": 8,
                    "gene_query": null,
                    "occurrence": null,
                    "gene_id": null,
                    "gene_name": null,
                    "strand": null,
                    "anchor_strand": "+",
                    "anchor_verified": false,
                    "sequence_source_type": "genbank_file",
                    "annotation_source_type": "genbank_file",
                    "sequence_source": "synthetic",
                    "annotation_source": "synthetic",
                    "sequence_sha1": null,
                    "annotation_sha1": null
                }
            ]
        }),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::SetParameter {
            name: "require_verified_genome_anchor_for_extension".to_string(),
            value: serde_json::json!(true),
        })
        .expect("enable strict anchor verification");

    let err = engine
        .apply(Operation::ExtendGenomeAnchor {
            seq_id: "anch".to_string(),
            side: GenomeAnchorSide::ThreePrime,
            length_bp: 2,
            output_id: Some("anch_ext".to_string()),
            catalog_path: None,
            cache_dir: None,
            prepared_genome_id: None,
        })
        .expect_err("strict verification should reject unverified anchor");
    assert!(matches!(err.code, ErrorCode::InvalidInput));
    assert!(err.message.contains("requires a verified genome anchor"));
}

#[test]
fn test_extend_genome_anchor_strict_verification_accepts_verified_anchor() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let gtf = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGT\n").unwrap();
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .unwrap();
    let catalog_path = root.join("catalog.json");
    let cache_dir = root.join("cache");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        gtf.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 2,
            end_1based: 7,
            output_id: Some("anch".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
        })
        .unwrap();
    engine
        .apply(Operation::SetParameter {
            name: "require_verified_genome_anchor_for_extension".to_string(),
            value: serde_json::json!(true),
        })
        .expect("enable strict anchor verification");

    let result = engine
        .apply(Operation::ExtendGenomeAnchor {
            seq_id: "anch".to_string(),
            side: GenomeAnchorSide::ThreePrime,
            length_bp: 3,
            output_id: Some("anch_ext".to_string()),
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
            prepared_genome_id: None,
        })
        .expect("verified anchor should be accepted under strict mode");
    assert_eq!(result.created_seq_ids, vec!["anch_ext".to_string()]);
}

#[test]
fn test_extend_genome_anchor_warns_when_clipped_at_chromosome_start() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let gtf = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGT\n").unwrap();
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .unwrap();
    let catalog_path = root.join("catalog.json");
    let cache_dir = root.join("cache");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        gtf.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut state = ProjectState::default();
    state.sequences.insert("anch".to_string(), seq("CGTAC"));
    state.metadata.insert(
        PROVENANCE_METADATA_KEY.to_string(),
        serde_json::json!({
            GENOME_EXTRACTIONS_METADATA_KEY: [
                {
                    "seq_id": "anch",
                    "recorded_at_unix_ms": 1,
                    "operation": "ExtractGenomeRegion",
                    "genome_id": "ToyGenome",
                    "catalog_path": catalog_path_str,
                    "cache_dir": null,
                    "chromosome": "chr1",
                    "start_1based": 2,
                    "end_1based": 6,
                    "gene_query": null,
                    "occurrence": null,
                    "gene_id": null,
                    "gene_name": null,
                    "strand": null,
                    "anchor_strand": "+",
                    "sequence_source_type": "local",
                    "annotation_source_type": "local",
                    "sequence_source": "local",
                    "annotation_source": "local",
                    "sequence_sha1": null,
                    "annotation_sha1": null
                }
            ]
        }),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();

    let result = engine
        .apply(Operation::ExtendGenomeAnchor {
            seq_id: "anch".to_string(),
            side: GenomeAnchorSide::FivePrime,
            length_bp: 10,
            output_id: Some("anch_ext".to_string()),
            catalog_path: None,
            cache_dir: None,
            prepared_genome_id: None,
        })
        .unwrap();
    assert_eq!(result.created_seq_ids, vec!["anch_ext".to_string()]);
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.contains("clipped at chromosome start position 1"))
    );
}

#[test]
fn test_verify_genome_anchor_records_unverified_status_in_provenance() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta = root.join("toy.fa");
    let gtf = root.join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGT\n").unwrap();
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .unwrap();
    let catalog_path = root.join("catalog.json");
    let cache_dir = root.join("cache");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta.display(),
        gtf.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 3,
            end_1based: 10,
            output_id: Some("anch".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
        })
        .unwrap();

    // Force mismatch before re-verification.
    engine
        .state_mut()
        .sequences
        .insert("anch".to_string(), seq("AAAAAAAA"));

    let result = engine
        .apply(Operation::VerifyGenomeAnchor {
            seq_id: "anch".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            prepared_genome_id: None,
        })
        .expect("verify anchor should succeed");
    assert_eq!(result.changed_seq_ids, vec!["anch".to_string()]);
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.contains("is unverified against")),
        "expected unverified warning in operation result"
    );
    let anchor_summary = engine
        .sequence_genome_anchor_summary("anch")
        .expect("anchor summary");
    assert_eq!(anchor_summary.anchor_verified, Some(false));

    let provenance = engine
        .state()
        .metadata
        .get(PROVENANCE_METADATA_KEY)
        .and_then(|v| v.as_object())
        .expect("provenance metadata object");
    let extractions = provenance
        .get(GENOME_EXTRACTIONS_METADATA_KEY)
        .and_then(|v| v.as_array())
        .expect("genome_extractions array");
    let verify_entry = extractions
        .iter()
        .rev()
        .find(|entry| {
            entry
                .get("seq_id")
                .and_then(|v| v.as_str())
                .map(|v| v == "anch")
                .unwrap_or(false)
                && entry
                    .get("operation")
                    .and_then(|v| v.as_str())
                    .map(|v| v == "VerifyGenomeAnchor")
                    .unwrap_or(false)
        })
        .expect("VerifyGenomeAnchor provenance entry");
    assert_eq!(
        verify_entry
            .get("anchor_verified")
            .and_then(|v| v.as_bool()),
        Some(false)
    );
    assert_eq!(
        verify_entry.get("catalog_path").and_then(|v| v.as_str()),
        Some(catalog_path_str.as_str())
    );
}

#[test]
fn test_import_genome_bed_track_supports_plain_and_gzip() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta_gz = root.join("toy.fa.gz");
    let ann_gz = root.join("toy.gtf.gz");
    write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\nACGT\n");
    write_gzip(
        &ann_gz,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
    );
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "description": "toy genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        file_url(&fasta_gz),
        file_url(&ann_gz),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 3,
            end_1based: 10,
            output_id: Some("toy_slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .unwrap();

    let bed_path = root.join("signals.bed");
    fs::write(
            &bed_path,
            "track name=toy\nchr1\t1\t4\tpeak_a\t42\t+\nchr1\t5\t12\tpeak_b\t7\t-\nchr2\t1\t4\twrong_chr\t50\t+\nchr1\tbad\t9\tbroken\n",
        )
        .unwrap();
    let plain = engine
        .apply(Operation::ImportGenomeBedTrack {
            seq_id: "toy_slice".to_string(),
            path: bed_path.to_string_lossy().to_string(),
            track_name: Some("chipseq_plain".to_string()),
            min_score: None,
            max_score: None,
            clear_existing: Some(true),
        })
        .unwrap();
    assert!(plain.changed_seq_ids.contains(&"toy_slice".to_string()));
    assert!(plain.warnings.iter().any(|w| w.contains("BED line")));

    let dna_plain = engine.state().sequences.get("toy_slice").unwrap();
    let generated_plain: Vec<_> = dna_plain
        .features()
        .iter()
        .filter(|f| {
            f.qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case(GENOME_BED_TRACK_GENERATED_TAG))
        })
        .collect();
    assert_eq!(generated_plain.len(), 2);
    assert!(generated_plain.iter().any(|f| {
        f.qualifier_values("label")
            .next()
            .map(|v| v.contains("peak_a"))
            .unwrap_or(false)
    }));

    let bed_gz = root.join("signals.bed.gz");
    write_gzip(
        &bed_gz,
        "chr1\t1\t4\tpeak_a\t42\t+\nchr1\t5\t12\tpeak_b\t.\t-\n",
    );
    let gz = engine
        .apply(Operation::ImportGenomeBedTrack {
            seq_id: "toy_slice".to_string(),
            path: bed_gz.to_string_lossy().to_string(),
            track_name: Some("chipseq_gz".to_string()),
            min_score: Some(10.0),
            max_score: None,
            clear_existing: Some(true),
        })
        .unwrap();
    assert!(gz.changed_seq_ids.contains(&"toy_slice".to_string()));
    assert!(gz.warnings.iter().any(|w| w.contains("score column")));

    let dna_gz = engine.state().sequences.get("toy_slice").unwrap();
    let generated_gz: Vec<_> = dna_gz
        .features()
        .iter()
        .filter(|f| {
            f.qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case(GENOME_BED_TRACK_GENERATED_TAG))
        })
        .collect();
    assert_eq!(generated_gz.len(), 1);
    assert!(
        generated_gz[0]
            .qualifier_values("label")
            .next()
            .map(|v| v.contains("peak_a"))
            .unwrap_or(false)
    );
}

#[test]
fn test_import_genome_bed_track_supports_concatenated_gzip_members() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta_gz = root.join("toy.fa.gz");
    let ann_gz = root.join("toy.gtf.gz");
    write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\nACGT\n");
    write_gzip(
        &ann_gz,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
    );
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "description": "toy genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        file_url(&fasta_gz),
        file_url(&ann_gz),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 1,
            end_1based: 12,
            output_id: Some("toy_slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .unwrap();

    let bed_gz = root.join("signals_concat.bed.gz");
    write_multi_member_gzip(
        &bed_gz,
        &[
            "chr1\t1\t4\tpeak_a\t42\t+\n",
            "chr1\t5\t8\tpeak_b\t25\t-\n",
            "chr1\t9\t12\tpeak_c\t18\t+\n",
        ],
    );
    let result = engine
        .apply(Operation::ImportGenomeBedTrack {
            seq_id: "toy_slice".to_string(),
            path: bed_gz.to_string_lossy().to_string(),
            track_name: Some("chipseq_concat_gz".to_string()),
            min_score: Some(20.0),
            max_score: None,
            clear_existing: Some(true),
        })
        .unwrap();
    assert!(result.changed_seq_ids.contains(&"toy_slice".to_string()));

    let dna = engine.state().sequences.get("toy_slice").unwrap();
    let generated: Vec<_> = dna
        .features()
        .iter()
        .filter(|f| {
            f.qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case(GENOME_BED_TRACK_GENERATED_TAG))
        })
        .collect();
    assert_eq!(generated.len(), 2);
    assert!(generated.iter().any(|f| {
        f.qualifier_values("label")
            .next()
            .map(|v| v.contains("peak_a"))
            .unwrap_or(false)
    }));
    assert!(generated.iter().any(|f| {
        f.qualifier_values("label")
            .next()
            .map(|v| v.contains("peak_b"))
            .unwrap_or(false)
    }));
}

#[cfg(unix)]
#[test]
fn test_import_genome_bigwig_track_uses_converter_and_filters_scores() {
    let td = tempdir().unwrap();
    let root = td.path();
    let fasta_gz = root.join("toy.fa.gz");
    let ann_gz = root.join("toy.gtf.gz");
    write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\nACGT\n");
    write_gzip(
        &ann_gz,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
    );
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "description": "toy genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        file_url(&fasta_gz),
        file_url(&ann_gz),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 3,
            end_1based: 10,
            output_id: Some("toy_slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .unwrap();

    let bed_path = root.join("baseline.bed");
    fs::write(&bed_path, "chr1\t1\t4\tpeak_a\t42\t+\n").unwrap();
    engine
        .apply(Operation::ImportGenomeBedTrack {
            seq_id: "toy_slice".to_string(),
            path: bed_path.to_string_lossy().to_string(),
            track_name: Some("baseline".to_string()),
            min_score: None,
            max_score: None,
            clear_existing: Some(true),
        })
        .unwrap();

    let converted_bedgraph = root.join("signals.from_bigwig.bedgraph");
    fs::write(
        &converted_bedgraph,
        "chr1\t1\t4\t0.5\nchr1\t5\t12\t2.0\nchr2\t1\t4\t7.0\nchr1\tbad\t9\tbroken\n",
    )
    .unwrap();
    let fake_converter = install_fake_bigwig_to_bedgraph(root, &converted_bedgraph);
    let _converter_guard = EnvVarGuard::set(BIGWIG_TO_BEDGRAPH_ENV_BIN, &fake_converter);

    let fake_bigwig = root.join("signals.bw");
    fs::write(&fake_bigwig, "placeholder").unwrap();
    let result = engine
        .apply(Operation::ImportGenomeBigWigTrack {
            seq_id: "toy_slice".to_string(),
            path: fake_bigwig.to_string_lossy().to_string(),
            track_name: Some("signal".to_string()),
            min_score: Some(1.0),
            max_score: Some(3.0),
            clear_existing: Some(true),
        })
        .unwrap();
    assert!(result.changed_seq_ids.contains(&"toy_slice".to_string()));
    assert!(result.warnings.iter().any(|w| w.contains("bedGraph line")));

    let dna = engine.state().sequences.get("toy_slice").unwrap();
    let bed_features: Vec<_> = dna
        .features()
        .iter()
        .filter(|f| {
            f.qualifier_values("gentle_generated")
                .any(|v| v.eq_ignore_ascii_case(GENOME_BED_TRACK_GENERATED_TAG))
        })
        .collect();
    assert_eq!(bed_features.len(), 0);
    let bigwig_features: Vec<_> = dna
        .features()
        .iter()
        .filter(|f| GentleEngine::is_generated_genome_bigwig_feature(f))
        .collect();
    assert_eq!(bigwig_features.len(), 1);
    assert_eq!(
        bigwig_features[0]
            .qualifier_values("score")
            .next()
            .unwrap_or_default(),
        "2.000000"
    );
}

#[test]
fn test_import_genome_vcf_track_supports_multiallelic_and_qual_filters() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("toy_slice".to_string(), seq("ACGTACGTACGTACGTACGT"));
    state.metadata.insert(
        PROVENANCE_METADATA_KEY.to_string(),
        serde_json::json!({
            GENOME_EXTRACTIONS_METADATA_KEY: [
                {
                    "seq_id": "toy_slice",
                    "recorded_at_unix_ms": 1,
                    "operation": "ExtractGenomeRegion",
                    "genome_id": "ToyGenome",
                    "catalog_path": "synthetic",
                    "cache_dir": null,
                    "chromosome": "chr1",
                    "start_1based": 100,
                    "end_1based": 119,
                    "gene_query": null,
                    "occurrence": null,
                    "gene_id": null,
                    "gene_name": null,
                    "strand": null,
                    "anchor_strand": "+",
                    "sequence_source_type": "synthetic",
                    "annotation_source_type": "synthetic",
                    "sequence_source": "synthetic",
                    "annotation_source": "synthetic",
                    "sequence_sha1": null,
                    "annotation_sha1": null
                }
            ]
        }),
    );
    let mut engine = GentleEngine::from_state(state);

    let td = tempdir().unwrap();
    let vcf_path = td.path().join("variants.vcf");
    fs::write(
            &vcf_path,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t101\trs1\tA\tG\t50\tPASS\tAC=1\nchr1\t105\t.\tA\tT,<DEL>\t.\tq10\tDP=5\nchr2\t101\trs2\tA\tC\t60\tPASS\tAC=1\nchr1\tbad\tbroken\n",
        )
        .unwrap();

    let plain = engine
        .apply(Operation::ImportGenomeVcfTrack {
            seq_id: "toy_slice".to_string(),
            path: vcf_path.to_string_lossy().to_string(),
            track_name: Some("variants".to_string()),
            min_score: None,
            max_score: None,
            clear_existing: Some(true),
        })
        .unwrap();
    assert!(plain.changed_seq_ids.contains(&"toy_slice".to_string()));
    assert!(plain.warnings.iter().any(|w| w.contains("VCF line")));
    assert!(
        plain
            .warnings
            .iter()
            .any(|w| w.contains("did not match anchor chromosome"))
    );
    assert!(plain.warnings.iter().any(|w| w.contains("chr2")));

    let dna_plain = engine.state().sequences.get("toy_slice").unwrap();
    let plain_features: Vec<_> = dna_plain
        .features()
        .iter()
        .filter(|f| GentleEngine::is_generated_genome_vcf_feature(f))
        .collect();
    assert_eq!(plain_features.len(), 3);
    assert!(
        plain_features
            .iter()
            .any(|f| { f.qualifier_values("vcf_alt").any(|v| v == "<DEL>") })
    );

    let vcf_gz = td.path().join("variants.vcf.gz");
    write_gzip(
        &vcf_gz,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t110\trs3\tC\tT\t25\tPASS\tAC=1\nchr1\t111\trs4\tC\tG\t.\tPASS\tAC=1\n",
    );
    let filtered = engine
        .apply(Operation::ImportGenomeVcfTrack {
            seq_id: "toy_slice".to_string(),
            path: vcf_gz.to_string_lossy().to_string(),
            track_name: Some("variants_gz".to_string()),
            min_score: Some(20.0),
            max_score: Some(30.0),
            clear_existing: Some(true),
        })
        .unwrap();
    assert!(filtered.changed_seq_ids.contains(&"toy_slice".to_string()));
    assert!(
        filtered
            .warnings
            .iter()
            .any(|w| w.contains("QUAL-based score filters"))
    );

    let dna_filtered = engine.state().sequences.get("toy_slice").unwrap();
    let filtered_features: Vec<_> = dna_filtered
        .features()
        .iter()
        .filter(|f| GentleEngine::is_generated_genome_vcf_feature(f))
        .collect();
    assert_eq!(filtered_features.len(), 1);
    assert_eq!(
        filtered_features[0]
            .qualifier_values("vcf_id")
            .next()
            .unwrap_or_default(),
        "rs3"
    );
}

#[test]
fn test_import_genome_vcf_track_chrom_alias_and_genotype_summary() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("toy_slice".to_string(), seq("ACGTACGTACGTACGTACGT"));
    state.metadata.insert(
        PROVENANCE_METADATA_KEY.to_string(),
        serde_json::json!({
            GENOME_EXTRACTIONS_METADATA_KEY: [
                {
                    "seq_id": "toy_slice",
                    "recorded_at_unix_ms": 1,
                    "operation": "ExtractGenomeRegion",
                    "genome_id": "ToyGenome",
                    "catalog_path": "synthetic",
                    "cache_dir": null,
                    "chromosome": "1",
                    "start_1based": 100,
                    "end_1based": 119,
                    "gene_query": null,
                    "occurrence": null,
                    "gene_id": null,
                    "gene_name": null,
                    "strand": null,
                    "anchor_strand": "+",
                    "sequence_source_type": "synthetic",
                    "annotation_source_type": "synthetic",
                    "sequence_source": "synthetic",
                    "annotation_source": "synthetic",
                    "sequence_sha1": null,
                    "annotation_sha1": null
                }
            ]
        }),
    );
    let mut engine = GentleEngine::from_state(state);
    let td = tempdir().unwrap();
    let vcf_path = td.path().join("genotype_variants.vcf");
    fs::write(
            &vcf_path,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_a\tsample_b\nchr01\t101\trsGT\tA\tG,TT\t51\tPASS\tAC=2;AN=4\tGT:DP\t0/1:12\t1|1:20\n",
        )
        .unwrap();
    let result = engine
        .apply(Operation::ImportGenomeVcfTrack {
            seq_id: "toy_slice".to_string(),
            path: vcf_path.to_string_lossy().to_string(),
            track_name: Some("gt_track".to_string()),
            min_score: None,
            max_score: None,
            clear_existing: Some(true),
        })
        .unwrap();
    assert!(result.changed_seq_ids.contains(&"toy_slice".to_string()));
    let dna = engine.state().sequences.get("toy_slice").unwrap();
    let features: Vec<_> = dna
        .features()
        .iter()
        .filter(|f| GentleEngine::is_generated_genome_vcf_feature(f))
        .collect();
    assert_eq!(features.len(), 2);
    let alt1 = features
        .iter()
        .find(|f| f.qualifier_values("vcf_alt").any(|v| v == "G"))
        .expect("ALT=G feature");
    assert_eq!(
        alt1.qualifier_values("vcf_variant_class")
            .next()
            .unwrap_or_default(),
        "SNP"
    );
    assert_eq!(
        alt1.qualifier_values("vcf_alt_allele_index")
            .next()
            .unwrap_or_default(),
        "1"
    );
    assert_eq!(
        alt1.qualifier_values("vcf_alt_carriers")
            .next()
            .unwrap_or_default(),
        "2"
    );
    assert_eq!(
        alt1.qualifier_values("vcf_gt_het")
            .next()
            .unwrap_or_default(),
        "1"
    );
    assert_eq!(
        alt1.qualifier_values("vcf_gt_hom_alt")
            .next()
            .unwrap_or_default(),
        "1"
    );
    assert_eq!(
        alt1.qualifier_values("vcf_phase")
            .next()
            .unwrap_or_default(),
        "mixed"
    );
    assert_eq!(
        alt1.qualifier_values("vcf_alt_carrier_samples")
            .next()
            .unwrap_or_default(),
        "sample_a,sample_b"
    );
}

#[test]
fn test_import_genome_vcf_track_can_cancel_via_progress_callback() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("toy_slice".to_string(), seq("ACGTACGTACGTACGTACGT"));
    state.metadata.insert(
        PROVENANCE_METADATA_KEY.to_string(),
        serde_json::json!({
            GENOME_EXTRACTIONS_METADATA_KEY: [
                {
                    "seq_id": "toy_slice",
                    "recorded_at_unix_ms": 1,
                    "operation": "ExtractGenomeRegion",
                    "genome_id": "ToyGenome",
                    "catalog_path": "synthetic",
                    "cache_dir": null,
                    "chromosome": "chr1",
                    "start_1based": 100,
                    "end_1based": 119,
                    "gene_query": null,
                    "occurrence": null,
                    "gene_id": null,
                    "gene_name": null,
                    "strand": null,
                    "anchor_strand": "+",
                    "sequence_source_type": "synthetic",
                    "annotation_source_type": "synthetic",
                    "sequence_source": "synthetic",
                    "annotation_source": "synthetic",
                    "sequence_sha1": null,
                    "annotation_sha1": null
                }
            ]
        }),
    );
    let mut engine = GentleEngine::from_state(state);
    let td = tempdir().unwrap();
    let vcf_path = td.path().join("many_variants.vcf");
    let mut payload =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for i in 0..600usize {
        payload.push_str(&format!(
            "chr1\t{}\trs{}\tA\tG\t50\tPASS\tAC=1\n",
            100 + (i % 20),
            i + 1
        ));
    }
    fs::write(&vcf_path, payload).unwrap();

    let mut cancelled = false;
    let result = engine
        .apply_with_progress(
            Operation::ImportGenomeVcfTrack {
                seq_id: "toy_slice".to_string(),
                path: vcf_path.to_string_lossy().to_string(),
                track_name: Some("many".to_string()),
                min_score: None,
                max_score: None,
                clear_existing: Some(true),
            },
            |progress| match progress {
                OperationProgress::GenomeTrackImport(p) => {
                    if !p.done && p.parsed_records >= 250 {
                        cancelled = true;
                        false
                    } else {
                        true
                    }
                }
                _ => true,
            },
        )
        .unwrap();
    assert!(cancelled);
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.contains("import cancelled"))
    );
    let dna = engine.state().sequences.get("toy_slice").unwrap();
    let features: Vec<_> = dna
        .features()
        .iter()
        .filter(|f| GentleEngine::is_generated_genome_vcf_feature(f))
        .collect();
    assert!(features.len() < 600);
}

#[test]
fn test_prepare_helper_genome_via_genbank_accession_and_extract() {
    let _guard = crate::genomes::genbank_env_lock().lock().unwrap();
    let td = tempdir().unwrap();
    let root = td.path();

    let mock_dir = root.join("mock");
    fs::create_dir_all(&mock_dir).unwrap();
    fs::copy("test_files/pGEX_3X.fa", mock_dir.join("L09137.fasta")).unwrap();
    fs::copy("test_files/pGEX-3X.gb", mock_dir.join("L09137.gbwithparts")).unwrap();

    let efetch_template = format!("file://{}/{{accession}}.{{rettype}}", mock_dir.display());
    let _efetch_env = EnvVarGuard::set("GENTLE_NCBI_EFETCH_URL", &efetch_template);

    let cache_dir = root.join("cache");
    let catalog_path = root.join("helper_catalog.json");
    let catalog_json = format!(
        r#"{{
  "Helper pUC19": {{
    "description": "helper vector from GenBank accession",
    "genbank_accession": "L09137",
    "cache_dir": "{}"
  }}
}}"#,
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    let prep = engine
        .apply(Operation::PrepareGenome {
            genome_id: "Helper pUC19".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    assert!(
        prep.messages
            .iter()
            .any(|m| m.contains("[genbank_accession]"))
    );

    let plan = GentleEngine::describe_reference_genome_sources(
        Some(&catalog_path_str),
        "Helper pUC19",
        None,
    )
    .unwrap();
    assert_eq!(plan.sequence_source_type, "genbank_accession");
    assert_eq!(plan.annotation_source_type, "genbank_accession");

    let genes =
        GentleEngine::list_reference_genome_genes(Some(&catalog_path_str), "Helper pUC19", None)
            .unwrap();
    assert!(!genes.is_empty());
    assert!(genes.iter().any(|g| {
        g.gene_name
            .as_ref()
            .map(|v| v.eq_ignore_ascii_case("bla"))
            .unwrap_or(false)
            || g.gene_id
                .as_ref()
                .map(|v| v.eq_ignore_ascii_case("bla"))
                .unwrap_or(false)
    }));

    let extract_gene = engine
        .apply(Operation::ExtractGenomeGene {
            genome_id: "Helper pUC19".to_string(),
            gene_query: "bla".to_string(),
            occurrence: Some(1),
            output_id: Some("helper_bla".to_string()),
            extract_mode: None,
            promoter_upstream_bp: None,
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
        })
        .unwrap();
    assert!(
        extract_gene
            .created_seq_ids
            .iter()
            .any(|id| id == "helper_bla")
    );
    let seq = engine.state().sequences.get("helper_bla").unwrap();
    assert!(seq.len() > 0);

    let extract_region = engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "Helper pUC19".to_string(),
            chromosome: "U13852.1".to_string(),
            start_1based: 1,
            end_1based: 40,
            output_id: Some("helper_head".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .unwrap();
    assert_eq!(
        extract_region.created_seq_ids,
        vec!["helper_head".to_string()]
    );
}

#[test]
fn test_extract_rebase_enzyme_names_from_text_handles_spaced_ecori() {
    let names = GentleEngine::extract_rebase_enzyme_names_from_text(
        "Multiple Cloning Site (MCS); contains BamHI, SmaI and EcoR I",
    );
    assert!(names.iter().any(|name| name == "BamHI"));
    assert!(names.iter().any(|name| name == "SmaI"));
    assert!(names.iter().any(|name| name == "EcoRI"));
}

#[test]
fn test_extract_helper_region_auto_annotates_puc_mcs() {
    for (genome_id, expected_preset, motif) in [
        ("Helper pUC19", "pUC19", PUC19_MCS_SEQUENCE),
        ("Helper pUC18", "pUC18", PUC18_MCS_SEQUENCE),
    ] {
        let td = tempdir().unwrap();
        let root = td.path();
        let cache_dir = root.join("cache");
        let seq_path = root.join("helper.fa");
        let ann_path = root.join("helper.gff3");
        let catalog_path = root.join("helper_catalog.json");

        let flank = "A".repeat(120);
        let tail = "T".repeat(120);
        let sequence = format!("{flank}{motif}{tail}");
        fs::write(&seq_path, format!(">chr_helper\n{}\n", sequence.as_str())).unwrap();
        fs::write(
            &ann_path,
            "##gff-version 3\nchr_helper\tGENtle\tgene\t5\t40\t.\t+\t.\tID=bla;Name=bla\n",
        )
        .unwrap();
        fs::write(
            &catalog_path,
            format!(
                r#"{{
  "{genome_id}": {{
    "description": "synthetic helper",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
                seq_path.display(),
                ann_path.display(),
                cache_dir.display()
            ),
        )
        .unwrap();
        let catalog_path_str = catalog_path.to_string_lossy().to_string();
        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::PrepareGenome {
                genome_id: genome_id.to_string(),
                catalog_path: Some(catalog_path_str.clone()),
                cache_dir: None,
                timeout_seconds: None,
            })
            .unwrap();
        let extract = engine
            .apply(Operation::ExtractGenomeRegion {
                genome_id: genome_id.to_string(),
                chromosome: "chr_helper".to_string(),
                start_1based: 1,
                end_1based: sequence.len(),
                output_id: Some("helper_slice".to_string()),
                annotation_scope: None,
                max_annotation_features: None,
                include_genomic_annotation: None,
                catalog_path: Some(catalog_path_str),
                cache_dir: None,
            })
            .unwrap();
        assert!(
            extract
                .messages
                .iter()
                .any(|line| line.contains("Annotated helper MCS")),
            "missing helper MCS message for {genome_id}: {:?}",
            extract.messages
        );
        let helper_slice = engine.state().sequences.get("helper_slice").unwrap();
        let mcs_features: Vec<_> = helper_slice
            .features()
            .iter()
            .filter(|feature| GentleEngine::is_generated_helper_mcs_feature(feature))
            .collect();
        assert_eq!(mcs_features.len(), 1, "expected one generated MCS feature");
        let mcs_feature = mcs_features[0];
        assert_eq!(
            GentleEngine::feature_qualifier_text(mcs_feature, "mcs_preset").as_deref(),
            Some(expected_preset)
        );
        assert!(
            GentleEngine::feature_qualifier_text(mcs_feature, "note")
                .as_deref()
                .map(|value| value.contains(PUC_MCS_EXPECTED_SITES))
                .unwrap_or(false)
        );
    }
}

#[test]
fn test_extract_helper_region_skips_mcs_annotation_when_motif_not_unique() {
    let td = tempdir().unwrap();
    let root = td.path();
    let cache_dir = root.join("cache");
    let seq_path = root.join("helper.fa");
    let ann_path = root.join("helper.gff3");
    let catalog_path = root.join("helper_catalog.json");
    let sequence = format!(
        "{}{}{}{}{}",
        "A".repeat(80),
        PUC19_MCS_SEQUENCE,
        "C".repeat(40),
        PUC19_MCS_SEQUENCE,
        "T".repeat(80)
    );
    fs::write(&seq_path, format!(">chr_helper\n{}\n", sequence.as_str())).unwrap();
    fs::write(
        &ann_path,
        "##gff-version 3\nchr_helper\tGENtle\tgene\t5\t40\t.\t+\t.\tID=bla;Name=bla\n",
    )
    .unwrap();
    fs::write(
        &catalog_path,
        format!(
            r#"{{
  "Helper pUC19": {{
    "description": "synthetic helper",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            seq_path.display(),
            ann_path.display(),
            cache_dir.display()
        ),
    )
    .unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "Helper pUC19".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    let extract = engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "Helper pUC19".to_string(),
            chromosome: "chr_helper".to_string(),
            start_1based: 1,
            end_1based: sequence.len(),
            output_id: Some("helper_slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .unwrap();
    assert!(
        extract
            .warnings
            .iter()
            .any(|line| line.contains("expected exactly one")),
        "missing uniqueness warning: {:?}",
        extract.warnings
    );
    let helper_slice = engine.state().sequences.get("helper_slice").unwrap();
    assert!(
        helper_slice
            .features()
            .iter()
            .all(|feature| !GentleEngine::is_generated_helper_mcs_feature(feature))
    );
}

#[test]
fn test_extract_helper_region_prefers_existing_mcs_annotation() {
    let td = tempdir().unwrap();
    let root = td.path();
    let cache_dir = root.join("cache");
    let seq_path = root.join("helper.fa");
    let ann_path = root.join("helper.gff3");
    let catalog_path = root.join("helper_catalog.json");
    let sequence = format!(
        "{}{}{}",
        "A".repeat(120),
        PUC19_MCS_SEQUENCE,
        "T".repeat(120)
    );
    fs::write(&seq_path, format!(">chr_helper\n{}\n", sequence.as_str())).unwrap();
    fs::write(
            &ann_path,
            "##gff-version 3\nchr_helper\tGENtle\tgene\t121\t180\t.\t+\t.\tID=mcs1;Name=Multiple Cloning Site (MCS) BamHI SmaI EcoR I\n",
        )
        .unwrap();
    fs::write(
        &catalog_path,
        format!(
            r#"{{
  "Helper pUC19": {{
    "description": "synthetic helper",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            seq_path.display(),
            ann_path.display(),
            cache_dir.display()
        ),
    )
    .unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "Helper pUC19".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    let extract = engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "Helper pUC19".to_string(),
            chromosome: "chr_helper".to_string(),
            start_1based: 1,
            end_1based: sequence.len(),
            output_id: Some("helper_slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .unwrap();
    assert!(
        extract
            .messages
            .iter()
            .any(|line| line.contains("Detected existing MCS annotation")),
        "missing source-annotation preference message: {:?}",
        extract.messages
    );
    let helper_slice = engine.state().sequences.get("helper_slice").unwrap();
    let source_mcs = helper_slice
        .features()
        .iter()
        .find(|feature| GentleEngine::feature_looks_like_mcs(feature))
        .expect("source MCS feature");
    let linked =
        GentleEngine::feature_qualifier_text(source_mcs, "mcs_expected_sites").unwrap_or_default();
    assert!(linked.contains("BamHI"), "missing BamHI in '{linked}'");
    assert!(linked.contains("EcoRI"), "missing EcoRI in '{linked}'");
    assert!(linked.contains("SmaI"), "missing SmaI in '{linked}'");
    assert!(
        helper_slice
            .features()
            .iter()
            .all(|feature| !GentleEngine::is_generated_helper_mcs_feature(feature))
    );
}

#[test]
fn test_extract_helper_region_does_not_apply_mcs_fallback_for_non_puc_ids() {
    let td = tempdir().unwrap();
    let root = td.path();
    let cache_dir = root.join("cache");
    let seq_path = root.join("helper.fa");
    let ann_path = root.join("helper.gff3");
    let catalog_path = root.join("helper_catalog.json");
    let sequence = format!(
        "{}{}{}",
        "A".repeat(120),
        PUC19_MCS_SEQUENCE,
        "T".repeat(120)
    );
    fs::write(&seq_path, format!(">chr_helper\n{}\n", sequence.as_str())).unwrap();
    fs::write(
        &ann_path,
        "##gff-version 3\nchr_helper\tGENtle\tgene\t5\t40\t.\t+\t.\tID=bla;Name=bla\n",
    )
    .unwrap();
    fs::write(
        &catalog_path,
        format!(
            r#"{{
  "Helper pGEX-like": {{
    "description": "synthetic helper",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            seq_path.display(),
            ann_path.display(),
            cache_dir.display()
        ),
    )
    .unwrap();
    let catalog_path_str = catalog_path.to_string_lossy().to_string();
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::PrepareGenome {
            genome_id: "Helper pGEX-like".to_string(),
            catalog_path: Some(catalog_path_str.clone()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .unwrap();
    let extract = engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: "Helper pGEX-like".to_string(),
            chromosome: "chr_helper".to_string(),
            start_1based: 1,
            end_1based: sequence.len(),
            output_id: Some("helper_slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path_str),
            cache_dir: None,
        })
        .unwrap();
    assert!(
        extract
            .messages
            .iter()
            .all(|line| !line.contains("helper MCS")),
        "unexpected helper MCS message: {:?}",
        extract.messages
    );
    let helper_slice = engine.state().sequences.get("helper_slice").unwrap();
    assert!(
        helper_slice
            .features()
            .iter()
            .all(|feature| !GentleEngine::is_generated_helper_mcs_feature(feature))
    );
}

#[test]
fn test_sync_tracked_genome_track_subscriptions_only_applies_to_new_anchors() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("anch1".to_string(), seq("ACGTACGTACGT"));
    state.metadata.insert(
        PROVENANCE_METADATA_KEY.to_string(),
        serde_json::json!({
            GENOME_EXTRACTIONS_METADATA_KEY: [
                {
                    "seq_id": "anch1",
                    "recorded_at_unix_ms": 1,
                    "operation": "ExtractGenomeRegion",
                    "genome_id": "ToyGenome",
                    "catalog_path": "synthetic",
                    "cache_dir": null,
                    "chromosome": "chr1",
                    "start_1based": 1,
                    "end_1based": 12,
                    "gene_query": null,
                    "occurrence": null,
                    "gene_id": null,
                    "gene_name": null,
                    "strand": null,
                    "anchor_strand": "+",
                    "sequence_source_type": null,
                    "annotation_source_type": null,
                    "sequence_source": null,
                    "annotation_source": null,
                    "sequence_sha1": null,
                    "annotation_sha1": null
                }
            ]
        }),
    );
    let mut engine = GentleEngine::from_state(state);

    let td = tempdir().unwrap();
    let bed_path = td.path().join("tracked.bed");
    std::fs::write(&bed_path, "chr1\t0\t4\tpeak1\t100\t+\n").unwrap();

    let inserted = engine
        .add_genome_track_subscription(GenomeTrackSubscription {
            source: GenomeTrackSource::Bed,
            path: bed_path.to_string_lossy().to_string(),
            track_name: Some("tracked".to_string()),
            min_score: None,
            max_score: None,
            clear_existing: true,
        })
        .unwrap();
    assert!(inserted);

    let first = engine
        .sync_tracked_genome_track_subscriptions(false)
        .unwrap();
    assert_eq!(first.target_sequences, 1);
    assert_eq!(first.applied_imports, 1);
    assert_eq!(first.failed_imports, 0);

    let anch1_features_after_first = engine
        .state()
        .sequences
        .get("anch1")
        .unwrap()
        .features()
        .iter()
        .filter(|f| GentleEngine::is_generated_genome_bed_feature(f))
        .count();
    assert_eq!(anch1_features_after_first, 1);

    engine
        .state_mut()
        .sequences
        .insert("anch2".to_string(), seq("ACGTACGTACGT"));
    let provenance = engine
        .state_mut()
        .metadata
        .get_mut(PROVENANCE_METADATA_KEY)
        .and_then(|v| v.as_object_mut())
        .unwrap();
    let records = provenance
        .get_mut(GENOME_EXTRACTIONS_METADATA_KEY)
        .and_then(|v| v.as_array_mut())
        .unwrap();
    records.push(serde_json::json!({
        "seq_id": "anch2",
        "recorded_at_unix_ms": 2,
        "operation": "ExtractGenomeRegion",
        "genome_id": "ToyGenome",
        "catalog_path": "synthetic",
        "cache_dir": null,
        "chromosome": "chr1",
        "start_1based": 1,
        "end_1based": 12,
        "gene_query": null,
        "occurrence": null,
        "gene_id": null,
        "gene_name": null,
        "strand": null,
        "anchor_strand": "+",
        "sequence_source_type": null,
        "annotation_source_type": null,
        "sequence_source": null,
        "annotation_source": null,
        "sequence_sha1": null,
        "annotation_sha1": null
    }));

    let second = engine
        .sync_tracked_genome_track_subscriptions(true)
        .unwrap();
    assert_eq!(second.target_sequences, 1);
    assert_eq!(second.applied_imports, 1);
    assert_eq!(second.failed_imports, 0);

    let anch1_features_after_second = engine
        .state()
        .sequences
        .get("anch1")
        .unwrap()
        .features()
        .iter()
        .filter(|f| GentleEngine::is_generated_genome_bed_feature(f))
        .count();
    assert_eq!(anch1_features_after_second, 1);
    let anch2_features = engine
        .state()
        .sequences
        .get("anch2")
        .unwrap()
        .features()
        .iter()
        .filter(|f| GentleEngine::is_generated_genome_bed_feature(f))
        .count();
    assert_eq!(anch2_features, 1);
}

#[test]
fn test_sync_tracked_genome_track_subscriptions_does_not_persist_empty_known_anchor_key() {
    let mut engine = GentleEngine::new();
    let report = engine
        .sync_tracked_genome_track_subscriptions(true)
        .expect("sync should succeed for empty state");
    assert_eq!(report.subscriptions_considered, 0);
    assert_eq!(report.target_sequences, 0);
    assert_eq!(report.applied_imports, 0);
    assert!(
        !engine
            .state()
            .metadata
            .contains_key(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY)
    );
}

#[test]
fn test_import_blast_hits_track_operation_adds_features_and_clears_previous() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("query".to_string(), seq("ACGTACGTACGTACGTACGTACGTACGT"));
    let mut engine = GentleEngine::from_state(state);

    let first = engine
        .apply(Operation::ImportBlastHitsTrack {
            seq_id: "query".to_string(),
            hits: vec![
                BlastHitFeatureInput {
                    subject_id: "chr1".to_string(),
                    query_start_1based: 1,
                    query_end_1based: 8,
                    subject_start_1based: 100,
                    subject_end_1based: 107,
                    identity_percent: 99.5,
                    bit_score: 42.0,
                    evalue: 1e-8,
                    query_coverage_percent: Some(100.0),
                },
                BlastHitFeatureInput {
                    subject_id: "chr2".to_string(),
                    query_start_1based: 20,
                    query_end_1based: 40,
                    subject_start_1based: 500,
                    subject_end_1based: 480,
                    identity_percent: 95.0,
                    bit_score: 33.0,
                    evalue: 1e-4,
                    query_coverage_percent: Some(75.0),
                },
            ],
            track_name: Some("blast_hits_demo".to_string()),
            clear_existing: Some(true),
            blast_provenance: Some(BlastInvocationProvenance {
                genome_id: "grch38".to_string(),
                query_label: "query".to_string(),
                query_length: 28,
                max_hits: 20,
                task: "blastn-short".to_string(),
                blastn_executable: "blastn".to_string(),
                blast_db_prefix: "/tmp/blastdb/genome".to_string(),
                command: vec![
                    "-db".to_string(),
                    "/tmp/blastdb/genome".to_string(),
                    "-query".to_string(),
                    "/tmp/query.fa".to_string(),
                ],
                command_line: "blastn -db /tmp/blastdb/genome -query /tmp/query.fa".to_string(),
                catalog_path: Some("assets/genomes.json".to_string()),
                cache_dir: Some("data/genomes".to_string()),
                options_override_json: None,
                effective_options_json: None,
            }),
        })
        .unwrap();
    assert!(first.changed_seq_ids.contains(&"query".to_string()));
    assert!(first.warnings.iter().any(|w| w.contains("was clamped")));
    assert!(
        first
            .messages
            .iter()
            .any(|m| m.contains("BLAST provenance") && m.contains("blastn -db"))
    );

    let dna = engine.state().sequences.get("query").unwrap();
    let blast_features: Vec<_> = dna
        .features()
        .iter()
        .filter(|f| GentleEngine::is_generated_blast_hit_feature(f))
        .collect();
    assert_eq!(blast_features.len(), 2);
    assert!(
        blast_features
            .iter()
            .any(|f| { f.qualifier_values("blast_subject_id").any(|v| v == "chr1") })
    );
    assert!(
        blast_features
            .iter()
            .any(|f| { f.qualifier_values("blast_subject_id").any(|v| v == "chr2") })
    );

    let second = engine
        .apply(Operation::ImportBlastHitsTrack {
            seq_id: "query".to_string(),
            hits: vec![BlastHitFeatureInput {
                subject_id: "chr3".to_string(),
                query_start_1based: 5,
                query_end_1based: 12,
                subject_start_1based: 1000,
                subject_end_1based: 1007,
                identity_percent: 98.0,
                bit_score: 50.0,
                evalue: 1e-12,
                query_coverage_percent: Some(90.0),
            }],
            track_name: Some("blast_hits_demo".to_string()),
            clear_existing: Some(true),
            blast_provenance: None,
        })
        .unwrap();
    assert!(second.changed_seq_ids.contains(&"query".to_string()));

    let dna = engine.state().sequences.get("query").unwrap();
    let blast_features_after_clear: Vec<_> = dna
        .features()
        .iter()
        .filter(|f| GentleEngine::is_generated_blast_hit_feature(f))
        .collect();
    assert_eq!(blast_features_after_clear.len(), 1);
    assert_eq!(
        blast_features_after_clear[0]
            .qualifier_values("blast_subject_id")
            .next()
            .unwrap_or_default(),
        "chr3"
    );
}

#[test]
fn test_resolve_blast_options_for_request_uses_project_legacy_and_request_layers() {
    let mut state = ProjectState::default();
    state.metadata.insert(
        BLAST_OPTIONS_OVERRIDE_METADATA_KEY.to_string(),
        serde_json::json!({
            "task": "blastn",
            "max_hits": 31,
            "thresholds": {
                "min_identity_percent": 91.0
            }
        }),
    );
    let engine = GentleEngine::from_state(state);
    let request = serde_json::json!({
        "max_hits": 7,
        "thresholds": {
            "min_identity_percent": 99.0
        }
    });
    let resolved = engine
        .resolve_blast_options_for_request(Some(&request), Some("blastn-short"), Some(5))
        .expect("resolve blast options");
    assert_eq!(resolved.task, "blastn-short");
    assert_eq!(resolved.max_hits, 7);
    assert_eq!(resolved.thresholds.min_identity_percent, Some(99.0));
}

#[test]
fn test_resolve_blast_options_rejects_unknown_keys() {
    let request = serde_json::json!({
        "unknown_key": 1
    });
    let err =
        GentleEngine::resolve_blast_options_with_layers(None, None, Some(&request), None, None)
            .expect_err("unknown key should be rejected");
    assert!(matches!(err.code, ErrorCode::InvalidInput));
    assert!(err.message.contains("Invalid BLAST options"));
}

#[test]
fn test_set_parameter_blast_options_metadata_roundtrip() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::SetParameter {
            name: "blast_options_override".to_string(),
            value: serde_json::json!({
                "max_hits": 13,
                "thresholds": {
                    "min_bit_score": 42.0
                }
            }),
        })
        .expect("set blast_options_override");
    assert_eq!(
        engine
            .state()
            .metadata
            .get(BLAST_OPTIONS_OVERRIDE_METADATA_KEY)
            .and_then(|v| v.get("max_hits"))
            .and_then(|v| v.as_u64()),
        Some(13)
    );

    engine
        .apply(Operation::SetParameter {
            name: "blast_options_defaults_path".to_string(),
            value: serde_json::json!("config/blast_defaults.json"),
        })
        .expect("set blast_options_defaults_path");
    assert_eq!(
        engine
            .state()
            .metadata
            .get(BLAST_OPTIONS_DEFAULTS_PATH_METADATA_KEY)
            .and_then(|v| v.as_str()),
        Some("config/blast_defaults.json")
    );

    engine
        .apply(Operation::SetParameter {
            name: "blast_options_override".to_string(),
            value: serde_json::Value::Null,
        })
        .expect("clear blast_options_override");
    assert!(
        !engine
            .state()
            .metadata
            .contains_key(BLAST_OPTIONS_OVERRIDE_METADATA_KEY)
    );
}

#[test]
fn test_apply_blast_thresholds_filters_hits_and_enforces_unique_best_hit() {
    let mut filtered = demo_blast_report();
    let thresholds = BlastThresholdOptions {
        min_identity_percent: Some(98.0),
        ..Default::default()
    };
    GentleEngine::apply_blast_thresholds_to_report(&mut filtered, &thresholds)
        .expect("threshold filtering");
    assert_eq!(filtered.hit_count, 1);
    assert_eq!(filtered.hits.len(), 1);
    assert_eq!(filtered.hits[0].subject_id, "chr1");
    assert!(
        filtered
            .warnings
            .iter()
            .any(|w| w.contains("thresholds removed 1 hit"))
    );

    let mut strict = demo_blast_report();
    let strict_thresholds = BlastThresholdOptions {
        min_identity_percent: Some(100.0),
        unique_best_hit: Some(true),
        ..Default::default()
    };
    let err = GentleEngine::apply_blast_thresholds_to_report(&mut strict, &strict_thresholds)
        .expect_err("unique_best_hit should fail when no hit remains");
    assert!(matches!(err.code, ErrorCode::InvalidInput));
    assert!(err.message.contains("unique_best_hit"));
}

#[test]
fn test_candidate_store_save_externalizes_and_load_hydrates() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seqA".to_string(), seq("ACGTACGTACGT"));
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "windows".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 4,
            step_bp: 4,
            feature_kinds: vec![],
            feature_label_regex: None,
            max_distance_bp: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(32),
        })
        .expect("generate candidates");

    let td = tempdir().expect("tempdir");
    let project_path = td.path().join("demo.gentle.json");
    engine
        .state()
        .save_to_path(project_path.to_string_lossy().as_ref())
        .expect("save project");

    let project_text = std::fs::read_to_string(&project_path).expect("read project");
    let project_json: serde_json::Value =
        serde_json::from_str(&project_text).expect("parse project");
    let candidate_meta = project_json
        .get("metadata")
        .and_then(|m| m.get(CANDIDATE_SETS_METADATA_KEY))
        .expect("candidate metadata present");
    assert_eq!(
        candidate_meta
            .get("schema")
            .and_then(|v| v.as_str())
            .unwrap_or_default(),
        CANDIDATE_SETS_REF_SCHEMA
    );

    let index_rel = candidate_meta
        .get("index_path")
        .and_then(|v| v.as_str())
        .expect("index path present");
    let index_abs = project_path
        .parent()
        .expect("project parent")
        .join(index_rel);
    assert!(index_abs.exists(), "candidate index file should exist");

    let loaded =
        ProjectState::load_from_path(project_path.to_string_lossy().as_ref()).expect("load");
    let loaded_meta = loaded
        .metadata
        .get(CANDIDATE_SETS_METADATA_KEY)
        .expect("hydrated candidate metadata");
    assert_eq!(
        loaded_meta
            .get("schema")
            .and_then(|v| v.as_str())
            .unwrap_or_default(),
        CANDIDATE_SETS_SCHEMA
    );
    let loaded_engine = GentleEngine::from_state(loaded);
    let summaries = loaded_engine.list_candidate_sets();
    assert_eq!(summaries.len(), 1);
    assert_eq!(summaries[0].name, "windows");
    assert_eq!(summaries[0].candidate_count, 3);
}

#[test]
fn test_project_state_load_degrades_when_candidate_sidecar_is_missing() {
    let _lock = candidate_store_env_lock().lock().unwrap();
    let _guard = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, "0");
    let td = tempdir().expect("tempdir");
    let project_path = td.path().join("broken.gentle.json");
    let project_json = serde_json::json!({
        "sequences": {},
        "metadata": {
            CANDIDATE_SETS_METADATA_KEY: {
                "schema": CANDIDATE_SETS_REF_SCHEMA,
                "storage": "jsonl_indexed",
                "index_path": "missing_sidecar/index.json",
                "set_count": 1,
                "updated_at_unix_ms": 0
            }
        }
    });
    std::fs::write(
        &project_path,
        serde_json::to_string_pretty(&project_json).expect("serialize project"),
    )
    .expect("write project");
    let loaded = ProjectState::load_from_path(project_path.to_string_lossy().as_ref())
        .expect("load in degraded mode");
    assert!(!loaded.metadata.contains_key(CANDIDATE_SETS_METADATA_KEY));
    let warning = loaded
        .metadata
        .get(CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY)
        .expect("degraded-load warning metadata");
    let warning_message = warning
        .get("message")
        .and_then(|v| v.as_str())
        .unwrap_or_default();
    assert!(warning_message.contains("candidate-store index"));
}

#[test]
fn test_project_state_load_strict_mode_errors_when_candidate_sidecar_is_missing() {
    let _lock = candidate_store_env_lock().lock().unwrap();
    let _guard = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, "1");
    let td = tempdir().expect("tempdir");
    let project_path = td.path().join("broken_strict.gentle.json");
    let project_json = serde_json::json!({
        "sequences": {},
        "metadata": {
            CANDIDATE_SETS_METADATA_KEY: {
                "schema": CANDIDATE_SETS_REF_SCHEMA,
                "storage": "jsonl_indexed",
                "index_path": "missing_sidecar/index.json",
                "set_count": 1,
                "updated_at_unix_ms": 0
            }
        }
    });
    std::fs::write(
        &project_path,
        serde_json::to_string_pretty(&project_json).expect("serialize project"),
    )
    .expect("write project");
    let err = ProjectState::load_from_path(project_path.to_string_lossy().as_ref()).unwrap_err();
    assert!(err.message.contains("candidate-store index"));
}

#[test]
fn test_candidate_store_save_replaces_sidecar_and_removes_stale_files() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seqA".to_string(), seq("ACGTACGTACGT"));
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "windows".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 4,
            step_bp: 2,
            feature_kinds: vec![],
            feature_label_regex: None,
            max_distance_bp: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(16),
        })
        .expect("generate candidates");

    let td = tempdir().expect("tempdir");
    let project_path = td.path().join("replace_sidecar.gentle.json");
    engine
        .state()
        .save_to_path(project_path.to_string_lossy().as_ref())
        .expect("initial save");

    let sidecar_dir = ProjectState::candidate_store_sidecar_dir(&project_path);
    let stale_path = sidecar_dir.join("stale.jsonl");
    std::fs::write(&stale_path, "{\"stale\":true}\n").expect("write stale sidecar file");
    assert!(
        stale_path.exists(),
        "stale file should exist before re-save"
    );

    engine
        .state()
        .save_to_path(project_path.to_string_lossy().as_ref())
        .expect("second save");

    assert!(
        !stale_path.exists(),
        "stale sidecar file should be removed by directory replacement"
    );
}

#[test]
fn test_candidate_generation_regex_anchor_and_filter_quantile_edges() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(0, 0),
        qualifiers: vec![("label".into(), Some("TP53".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(11, 11),
        qualifiers: vec![("label".into(), Some("TP53-AS1".to_string()))],
    });
    state.sequences.insert("seqA".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "gene_all".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 1,
            step_bp: 1,
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: None,
            max_distance_bp: Some(0),
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(64),
        })
        .expect("generate all gene-anchored candidates");
    let gene_all_count = engine
        .list_candidate_sets()
        .into_iter()
        .find(|s| s.name == "gene_all")
        .map(|s| s.candidate_count)
        .expect("gene_all set exists");

    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "tp53_only".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 1,
            step_bp: 1,
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: Some("^TP53$".to_string()),
            max_distance_bp: Some(0),
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(64),
        })
        .expect("generate regex-anchored candidates");
    let regex_count = engine
        .list_candidate_sets()
        .into_iter()
        .find(|s| s.name == "tp53_only")
        .map(|s| s.candidate_count)
        .expect("regex set exists");
    assert!(regex_count >= 1);
    assert!(regex_count < gene_all_count);

    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "windows".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 4,
            step_bp: 2,
            feature_kinds: vec![],
            feature_label_regex: None,
            max_distance_bp: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(64),
        })
        .expect("generate windows");
    let windows_count = engine
        .list_candidate_sets()
        .into_iter()
        .find(|s| s.name == "windows")
        .map(|s| s.candidate_count)
        .expect("windows set exists");

    engine
        .apply(Operation::FilterCandidateSet {
            input_set: "windows".to_string(),
            output_set: "windows_all_q".to_string(),
            metric: "gc_fraction".to_string(),
            min: None,
            max: None,
            min_quantile: Some(0.0),
            max_quantile: Some(1.0),
        })
        .expect("quantile full-range filter");
    let windows_all_q_count = engine
        .list_candidate_sets()
        .into_iter()
        .find(|s| s.name == "windows_all_q")
        .map(|s| s.candidate_count)
        .expect("windows_all_q exists");
    assert_eq!(windows_all_q_count, windows_count);

    engine
        .apply(Operation::FilterCandidateSet {
            input_set: "windows".to_string(),
            output_set: "windows_top_q".to_string(),
            metric: "gc_fraction".to_string(),
            min: None,
            max: None,
            min_quantile: Some(1.0),
            max_quantile: Some(1.0),
        })
        .expect("top quantile filter");
    let top_count = engine
        .list_candidate_sets()
        .into_iter()
        .find(|s| s.name == "windows_top_q")
        .map(|s| s.candidate_count)
        .expect("windows_top_q exists");
    assert!(top_count >= 1);
    assert!(top_count <= windows_count);

    let err = engine
        .apply(Operation::FilterCandidateSet {
            input_set: "windows".to_string(),
            output_set: "bad_q".to_string(),
            metric: "gc_fraction".to_string(),
            min: None,
            max: None,
            min_quantile: Some(1.1),
            max_quantile: None,
        })
        .unwrap_err();
    assert!(err.message.contains("between 0 and 1"));
}

#[test]
fn test_extract_anchored_region_supports_middle_feature_boundary() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(6, 14),
        qualifiers: vec![("label".into(), Some("MID_GENE".to_string()))],
    });
    state.sequences.insert("seqA".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let res = engine
        .apply(Operation::ExtractAnchoredRegion {
            input: "seqA".to_string(),
            anchor: SequenceAnchor::FeatureBoundary {
                feature_kind: Some("gene".to_string()),
                feature_label: Some("MID_GENE".to_string()),
                boundary: AnchorBoundary::Middle,
                occurrence: Some(0),
            },
            direction: AnchorDirection::Upstream,
            target_length_bp: 4,
            length_tolerance_bp: 0,
            required_re_sites: vec![],
            required_tf_motifs: vec![],
            forward_primer: None,
            reverse_primer: None,
            output_prefix: Some("mid".to_string()),
            unique: Some(true),
            max_candidates: Some(1),
        })
        .expect("extract using middle boundary");
    assert_eq!(res.created_seq_ids, vec!["mid_1".to_string()]);
    let seq = engine
        .state()
        .sequences
        .get("mid_1")
        .expect("created middle-boundary extract");
    assert_eq!(seq.len(), 4);
}

#[test]
fn test_generate_candidate_set_between_two_sequence_anchors() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGTACGT").unwrap(),
    );
    let mut engine = GentleEngine::from_state(state);

    let res = engine
        .apply(Operation::GenerateCandidateSetBetweenAnchors {
            set_name: "between".to_string(),
            seq_id: "seqA".to_string(),
            anchor_a: SequenceAnchor::Position { zero_based: 2 },
            anchor_b: SequenceAnchor::Position { zero_based: 10 },
            length_bp: 4,
            step_bp: 2,
            limit: Some(32),
        })
        .expect("generate between anchors");
    assert!(res.messages.iter().any(|m| m.contains("between anchors")));

    let summary = engine
        .list_candidate_sets()
        .into_iter()
        .find(|set| set.name == "between")
        .expect("between set summary");
    assert_eq!(summary.candidate_count, 3);

    let (page, total, _) = engine
        .inspect_candidate_set_page("between", 64, 0)
        .expect("inspect between candidates");
    assert_eq!(total, 3);
    assert_eq!(page.candidates[0].start_0based, 2);
    assert_eq!(page.candidates[0].end_0based, 6);
    assert_eq!(
        page.candidates[0]
            .metrics
            .get("distance_to_anchor_a_bp")
            .copied(),
        Some(0.0)
    );
    assert_eq!(
        page.candidates[0]
            .metrics
            .get("anchor_interval_span_bp")
            .copied(),
        Some(8.0)
    );
}

#[test]
fn test_candidate_generation_feature_parts_ignores_multipart_gaps() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "exon".into(),
        location: gb_io::seq::Location::Join(vec![
            gb_io::seq::Location::simple_range(2, 4),
            gb_io::seq::Location::simple_range(10, 12),
        ]),
        qualifiers: vec![("label".into(), Some("EXON_JOIN".to_string()))],
    });
    state.sequences.insert("seqA".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "span_mode".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 1,
            step_bp: 1,
            feature_kinds: vec!["exon".to_string()],
            feature_label_regex: Some("^EXON_JOIN$".to_string()),
            max_distance_bp: Some(0),
            feature_geometry_mode: Some(CandidateFeatureGeometryMode::FeatureSpan),
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(256),
        })
        .expect("generate span-mode candidates");

    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "parts_mode".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 1,
            step_bp: 1,
            feature_kinds: vec!["exon".to_string()],
            feature_label_regex: Some("^EXON_JOIN$".to_string()),
            max_distance_bp: Some(0),
            feature_geometry_mode: Some(CandidateFeatureGeometryMode::FeatureParts),
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(256),
        })
        .expect("generate parts-mode candidates");

    let span_count = engine
        .list_candidate_sets()
        .into_iter()
        .find(|s| s.name == "span_mode")
        .map(|s| s.candidate_count)
        .expect("span_mode exists");
    let parts_count = engine
        .list_candidate_sets()
        .into_iter()
        .find(|s| s.name == "parts_mode")
        .map(|s| s.candidate_count)
        .expect("parts_mode exists");
    assert!(
        parts_count < span_count,
        "feature_parts should not fill multipart feature gaps"
    );
}

#[test]
fn test_candidate_distance_feature_boundaries_respects_five_prime_and_three_prime() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(5, 8),
        qualifiers: vec![("label".into(), Some("PLUS_GENE".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::Complement(Box::new(gb_io::seq::Location::simple_range(
            12, 15,
        ))),
        qualifiers: vec![("label".into(), Some("MINUS_GENE".to_string()))],
    });
    state.sequences.insert("seqA".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "windows".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 1,
            step_bp: 1,
            feature_kinds: vec![],
            feature_label_regex: None,
            max_distance_bp: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(512),
        })
        .expect("generate windows");

    engine
        .apply(Operation::ScoreCandidateSetDistance {
            set_name: "windows".to_string(),
            metric: "dist_5p".to_string(),
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: None,
            feature_geometry_mode: Some(CandidateFeatureGeometryMode::FeatureBoundaries),
            feature_boundary_mode: Some(CandidateFeatureBoundaryMode::FivePrime),
            feature_strand_relation: None,
        })
        .expect("score five-prime distance");
    engine
        .apply(Operation::ScoreCandidateSetDistance {
            set_name: "windows".to_string(),
            metric: "dist_3p".to_string(),
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: None,
            feature_geometry_mode: Some(CandidateFeatureGeometryMode::FeatureBoundaries),
            feature_boundary_mode: Some(CandidateFeatureBoundaryMode::ThreePrime),
            feature_strand_relation: None,
        })
        .expect("score three-prime distance");

    let (page, _, _) = engine
        .inspect_candidate_set_page("windows", 2048, 0)
        .expect("inspect windows");
    let at_pos = |pos: usize, metric: &str| -> f64 {
        page.candidates
            .iter()
            .find(|candidate| candidate.start_0based == pos)
            .and_then(|candidate| candidate.metrics.get(metric).copied())
            .unwrap_or(f64::NAN)
    };
    let plus_start_dist_5p = at_pos(5, "dist_5p");
    let plus_end_dist_5p = at_pos(7, "dist_5p");
    let plus_start_dist_3p = at_pos(5, "dist_3p");
    let plus_end_dist_3p = at_pos(7, "dist_3p");
    assert!(plus_start_dist_5p < plus_end_dist_5p);
    assert!(plus_end_dist_3p < plus_start_dist_3p);

    let minus_end_dist_5p = at_pos(14, "dist_5p");
    let minus_start_dist_5p = at_pos(12, "dist_5p");
    let minus_end_dist_3p = at_pos(14, "dist_3p");
    let minus_start_dist_3p = at_pos(12, "dist_3p");
    assert!(minus_end_dist_5p < minus_start_dist_5p);
    assert!(minus_start_dist_3p < minus_end_dist_3p);
}

#[test]
fn test_candidate_generation_feature_strand_relation_filters_plus_and_minus() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(2, 5),
        qualifiers: vec![("label".into(), Some("PLUS_GENE".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::Complement(Box::new(gb_io::seq::Location::simple_range(
            12, 15,
        ))),
        qualifiers: vec![("label".into(), Some("MINUS_GENE".to_string()))],
    });
    state.sequences.insert("seqA".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    for (name, strand_relation) in [
        ("any", CandidateFeatureStrandRelation::Any),
        ("same", CandidateFeatureStrandRelation::Same),
        ("opposite", CandidateFeatureStrandRelation::Opposite),
    ] {
        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: format!("gene_{name}"),
                seq_id: "seqA".to_string(),
                length_bp: 1,
                step_bp: 1,
                feature_kinds: vec!["gene".to_string()],
                feature_label_regex: None,
                max_distance_bp: Some(0),
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: Some(strand_relation),
                limit: Some(256),
            })
            .expect("generate gene candidates with strand relation");
    }

    let count_for = |set_name: &str| -> usize {
        engine
            .list_candidate_sets()
            .into_iter()
            .find(|set| set.name == set_name)
            .map(|set| set.candidate_count)
            .expect("candidate set exists")
    };
    let any_count = count_for("gene_any");
    let same_count = count_for("gene_same");
    let opposite_count = count_for("gene_opposite");
    assert!(same_count > 0);
    assert!(opposite_count > 0);
    assert_eq!(same_count + opposite_count, any_count);
}

#[test]
fn test_candidate_distance_feature_strand_relation_prefers_matching_strand_features() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(2, 5),
        qualifiers: vec![("label".into(), Some("PLUS_GENE".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::Complement(Box::new(gb_io::seq::Location::simple_range(
            12, 15,
        ))),
        qualifiers: vec![("label".into(), Some("MINUS_GENE".to_string()))],
    });
    state.sequences.insert("seqA".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "windows".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 1,
            step_bp: 1,
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: None,
            max_distance_bp: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(256),
        })
        .expect("generate candidate windows");

    engine
        .apply(Operation::ScoreCandidateSetDistance {
            set_name: "windows".to_string(),
            metric: "dist_any".to_string(),
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: Some(CandidateFeatureStrandRelation::Any),
        })
        .expect("score distance any");
    engine
        .apply(Operation::ScoreCandidateSetDistance {
            set_name: "windows".to_string(),
            metric: "dist_same".to_string(),
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: Some(CandidateFeatureStrandRelation::Same),
        })
        .expect("score distance same");
    engine
        .apply(Operation::ScoreCandidateSetDistance {
            set_name: "windows".to_string(),
            metric: "dist_opposite".to_string(),
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: Some(CandidateFeatureStrandRelation::Opposite),
        })
        .expect("score distance opposite");

    let (page, _, _) = engine
        .inspect_candidate_set_page("windows", 4096, 0)
        .expect("inspect scored windows");
    let at_pos = |pos: usize, metric: &str| -> f64 {
        page.candidates
            .iter()
            .find(|candidate| candidate.start_0based == pos)
            .and_then(|candidate| candidate.metrics.get(metric).copied())
            .unwrap_or(f64::NAN)
    };

    let plus_any = at_pos(2, "dist_any");
    let plus_same = at_pos(2, "dist_same");
    let plus_opposite = at_pos(2, "dist_opposite");
    assert_eq!(plus_any, 0.0);
    assert_eq!(plus_same, 0.0);
    assert!(plus_opposite > 0.0);

    let minus_any = at_pos(14, "dist_any");
    let minus_same = at_pos(14, "dist_same");
    let minus_opposite = at_pos(14, "dist_opposite");
    assert_eq!(minus_any, 0.0);
    assert!(minus_same > 0.0);
    assert_eq!(minus_opposite, 0.0);
}

#[test]
fn test_candidate_optimizer_weighted_topk_and_pareto() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("GCATGAAA").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "cand".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 2,
            step_bp: 2,
            feature_kinds: vec![],
            feature_label_regex: None,
            max_distance_bp: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(64),
        })
        .expect("generate candidates");

    engine
        .apply(Operation::ScoreCandidateSetWeightedObjective {
            set_name: "cand".to_string(),
            metric: "objective".to_string(),
            objectives: vec![
                CandidateWeightedObjectiveTerm {
                    metric: "gc_fraction".to_string(),
                    weight: 0.7,
                    direction: CandidateObjectiveDirection::Maximize,
                },
                CandidateWeightedObjectiveTerm {
                    metric: "distance_to_seq_start_bp".to_string(),
                    weight: 0.3,
                    direction: CandidateObjectiveDirection::Minimize,
                },
            ],
            normalize_metrics: Some(true),
        })
        .expect("score weighted objective");

    engine
        .apply(Operation::TopKCandidateSet {
            input_set: "cand".to_string(),
            output_set: "cand_top1".to_string(),
            metric: "objective".to_string(),
            k: 1,
            direction: Some(CandidateObjectiveDirection::Maximize),
            tie_break: Some(CandidateTieBreakPolicy::SeqStartEnd),
        })
        .expect("top-k selection");
    let (top_page, _, _) = engine
        .inspect_candidate_set_page("cand_top1", 10, 0)
        .expect("inspect top-k result");
    assert_eq!(top_page.candidates.len(), 1);
    assert_eq!(top_page.candidates[0].start_0based, 0);

    engine
        .apply(Operation::ParetoFrontierCandidateSet {
            input_set: "cand".to_string(),
            output_set: "cand_pareto".to_string(),
            objectives: vec![
                CandidateObjectiveSpec {
                    metric: "gc_fraction".to_string(),
                    direction: CandidateObjectiveDirection::Maximize,
                },
                CandidateObjectiveSpec {
                    metric: "distance_to_seq_start_bp".to_string(),
                    direction: CandidateObjectiveDirection::Minimize,
                },
            ],
            max_candidates: None,
            tie_break: Some(CandidateTieBreakPolicy::SeqStartEnd),
        })
        .expect("pareto frontier selection");
    let pareto_summary = engine
        .list_candidate_sets()
        .into_iter()
        .find(|set| set.name == "cand_pareto")
        .expect("pareto set summary");
    assert_eq!(pareto_summary.candidate_count, 1);
}

#[test]
fn test_topk_tie_break_policy_is_deterministic() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::GenerateCandidateSet {
            set_name: "cand".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 2,
            step_bp: 1,
            feature_kinds: vec![],
            feature_label_regex: None,
            max_distance_bp: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: Some(64),
        })
        .expect("generate candidates");

    engine
        .apply(Operation::TopKCandidateSet {
            input_set: "cand".to_string(),
            output_set: "cand_top2".to_string(),
            metric: "length_bp".to_string(),
            k: 2,
            direction: Some(CandidateObjectiveDirection::Maximize),
            tie_break: Some(CandidateTieBreakPolicy::SeqStartEnd),
        })
        .expect("top-k tie-break");
    let (page, _, _) = engine
        .inspect_candidate_set_page("cand_top2", 10, 0)
        .expect("inspect top2");
    let starts = page
        .candidates
        .iter()
        .map(|candidate| candidate.start_0based)
        .collect::<Vec<_>>();
    assert_eq!(starts, vec![0, 1]);
}

#[test]
fn test_candidate_macro_template_store_and_render() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    engine
        .apply(Operation::UpsertCandidateMacroTemplate {
            name: "scan_tp53".to_string(),
            description: Some("demo template".to_string()),
            details_url: Some("https://example.org/candidates/scan-tp53".to_string()),
            parameters: vec![
                CandidateMacroTemplateParam {
                    name: "set_name".to_string(),
                    default_value: None,
                    required: true,
                },
                CandidateMacroTemplateParam {
                    name: "seq_id".to_string(),
                    default_value: Some("seqA".to_string()),
                    required: false,
                },
                CandidateMacroTemplateParam {
                    name: "len".to_string(),
                    default_value: Some("20".to_string()),
                    required: false,
                },
            ],
            script: "generate ${set_name} ${seq_id} --length ${len} --step 1".to_string(),
        })
        .expect("upsert template");

    let templates = engine.list_candidate_macro_templates();
    assert_eq!(templates.len(), 1);
    assert_eq!(templates[0].name, "scan_tp53");
    assert_eq!(
        templates[0].details_url.as_deref(),
        Some("https://example.org/candidates/scan-tp53")
    );
    let template = engine
        .get_candidate_macro_template("scan_tp53")
        .expect("get candidate template");
    assert_eq!(
        template.details_url.as_deref(),
        Some("https://example.org/candidates/scan-tp53")
    );

    let rendered = engine
        .render_candidate_macro_template_script(
            "scan_tp53",
            &HashMap::from([("set_name".to_string(), "tp53_hits".to_string())]),
        )
        .expect("render template");
    assert_eq!(rendered, "generate tp53_hits seqA --length 20 --step 1");

    engine
        .apply(Operation::DeleteCandidateMacroTemplate {
            name: "scan_tp53".to_string(),
        })
        .expect("delete template");
    assert!(engine.list_candidate_macro_templates().is_empty());
}

#[test]
fn test_candidate_macro_template_rejects_non_http_details_url() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let err = engine
        .apply(Operation::UpsertCandidateMacroTemplate {
            name: "bad_url".to_string(),
            description: Some("template".to_string()),
            details_url: Some("ftp://example.org/template".to_string()),
            parameters: vec![CandidateMacroTemplateParam {
                name: "set_name".to_string(),
                default_value: None,
                required: true,
            }],
            script: "generate ${set_name} seqA --length 20".to_string(),
        })
        .expect_err("expected invalid details_url error");
    assert!(
        err.message.contains("must start with http:// or https://"),
        "unexpected error: {}",
        err.message
    );
}

#[test]
fn test_workflow_macro_template_store_and_render() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    engine
        .apply(Operation::UpsertWorkflowMacroTemplate {
            name: "clone_step".to_string(),
            description: Some("demo workflow template".to_string()),
            details_url: Some("https://example.org/cloning/clone-step".to_string()),
            parameters: vec![
                WorkflowMacroTemplateParam {
                    name: "seq_id".to_string(),
                    default_value: None,
                    required: true,
                },
                WorkflowMacroTemplateParam {
                    name: "out_id".to_string(),
                    default_value: Some("seqA_rev".to_string()),
                    required: false,
                },
            ],
            input_ports: vec![],
            output_ports: vec![],
            script: "op {\"Reverse\":{\"input\":\"${seq_id}\",\"output_id\":\"${out_id}\"}}"
                .to_string(),
        })
        .expect("upsert workflow template");

    let templates = engine.list_workflow_macro_templates();
    assert_eq!(templates.len(), 1);
    assert_eq!(templates[0].name, "clone_step");
    let template = engine
        .get_workflow_macro_template("clone_step")
        .expect("get workflow template");
    assert_eq!(template.template_schema, CLONING_MACRO_TEMPLATE_SCHEMA);
    assert_eq!(
        template.details_url.as_deref(),
        Some("https://example.org/cloning/clone-step")
    );

    let rendered = engine
        .render_workflow_macro_template_script(
            "clone_step",
            &HashMap::from([("seq_id".to_string(), "seqX".to_string())]),
        )
        .expect("render workflow template");
    assert_eq!(
        rendered,
        "op {\"Reverse\":{\"input\":\"seqX\",\"output_id\":\"seqA_rev\"}}"
    );

    engine
        .apply(Operation::DeleteWorkflowMacroTemplate {
            name: "clone_step".to_string(),
        })
        .expect("delete workflow template");
    assert!(engine.list_workflow_macro_templates().is_empty());
}

#[test]
fn test_workflow_macro_template_rejects_non_http_details_url() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let err = engine
        .apply(Operation::UpsertWorkflowMacroTemplate {
            name: "bad_url".to_string(),
            description: Some("template".to_string()),
            details_url: Some("ftp://example.org/template".to_string()),
            parameters: vec![WorkflowMacroTemplateParam {
                name: "seq_id".to_string(),
                default_value: None,
                required: true,
            }],
            input_ports: vec![],
            output_ports: vec![],
            script: "op {\"Reverse\":{\"input\":\"${seq_id}\"}}".to_string(),
        })
        .expect_err("expected invalid details_url error");
    assert!(
        err.message.contains("must start with http:// or https://"),
        "unexpected error: {}",
        err.message
    );
}

#[test]
fn test_guide_design_filter_generate_and_export() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    engine
        .apply(Operation::UpsertGuideSet {
            guide_set_id: "tp73_guides".to_string(),
            guides: vec![
                GuideCandidate {
                    guide_id: "g1".to_string(),
                    seq_id: "tp73".to_string(),
                    start_0based: 100,
                    end_0based_exclusive: 120,
                    strand: "+".to_string(),
                    protospacer: "GACCTGTTGACGATGTTCCA".to_string(),
                    pam: "AGG".to_string(),
                    nuclease: "SpCas9".to_string(),
                    cut_offset_from_protospacer_start: 17,
                    rank: Some(1),
                },
                GuideCandidate {
                    guide_id: "g2".to_string(),
                    seq_id: "tp73".to_string(),
                    start_0based: 200,
                    end_0based_exclusive: 220,
                    strand: "+".to_string(),
                    protospacer: "ATATATATATATATATATAT".to_string(),
                    pam: "AGG".to_string(),
                    nuclease: "SpCas9".to_string(),
                    cut_offset_from_protospacer_start: 17,
                    rank: Some(2),
                },
                GuideCandidate {
                    guide_id: "g3".to_string(),
                    seq_id: "tp73".to_string(),
                    start_0based: 300,
                    end_0based_exclusive: 320,
                    strand: "+".to_string(),
                    protospacer: "GACTTTTGACTGACTGACT".to_string(),
                    pam: "AGG".to_string(),
                    nuclease: "SpCas9".to_string(),
                    cut_offset_from_protospacer_start: 17,
                    rank: Some(3),
                },
            ],
        })
        .expect("upsert guide set");

    engine
        .apply(Operation::FilterGuidesPractical {
            guide_set_id: "tp73_guides".to_string(),
            config: GuidePracticalFilterConfig {
                gc_min: Some(0.30),
                gc_max: Some(0.70),
                max_homopolymer_run: Some(4),
                max_homopolymer_run_per_base: HashMap::new(),
                reject_ambiguous_bases: true,
                avoid_u6_terminator_tttt: true,
                u6_terminator_window: GuideU6TerminatorWindow::SpacerPlusTail,
                max_dinucleotide_repeat_units: None,
                forbidden_motifs: vec![],
                required_5prime_base: None,
                allow_5prime_g_extension: true,
            },
            output_guide_set_id: Some("tp73_guides_pass".to_string()),
        })
        .expect("filter guides");
    let report = engine
        .get_guide_practical_filter_report("tp73_guides")
        .expect("guide report");
    assert_eq!(report.passed_count, 1);
    assert_eq!(report.rejected_count, 2);

    engine
        .apply(Operation::GenerateGuideOligos {
            guide_set_id: "tp73_guides".to_string(),
            template_id: "lenti_bsmbi_u6_default".to_string(),
            apply_5prime_g_extension: Some(true),
            output_oligo_set_id: Some("tp73_oligos".to_string()),
            passed_only: Some(true),
        })
        .expect("generate guide oligos");

    let dir = tempdir().expect("tempdir");
    let csv_path = dir.path().join("guides.csv").display().to_string();
    let fasta_path = dir.path().join("guides.fa").display().to_string();
    let protocol_path = dir.path().join("guides.protocol.txt").display().to_string();

    engine
        .apply(Operation::ExportGuideOligos {
            guide_set_id: "tp73_guides".to_string(),
            oligo_set_id: Some("tp73_oligos".to_string()),
            format: GuideOligoExportFormat::CsvTable,
            path: csv_path.clone(),
            plate_format: None,
        })
        .expect("export guide oligos csv");
    engine
        .apply(Operation::ExportGuideOligos {
            guide_set_id: "tp73_guides".to_string(),
            oligo_set_id: Some("tp73_oligos".to_string()),
            format: GuideOligoExportFormat::Fasta,
            path: fasta_path.clone(),
            plate_format: None,
        })
        .expect("export guide oligos fasta");
    engine
        .apply(Operation::ExportGuideProtocolText {
            guide_set_id: "tp73_guides".to_string(),
            oligo_set_id: Some("tp73_oligos".to_string()),
            path: protocol_path.clone(),
            include_qc_checklist: Some(true),
        })
        .expect("export guide protocol");

    let csv = fs::read_to_string(csv_path).expect("read csv export");
    assert!(csv.contains("guide_id,rank,forward_oligo,reverse_oligo,notes"));
    assert!(csv.contains("g1"));
    let fasta = fs::read_to_string(fasta_path).expect("read fasta export");
    assert!(fasta.contains(">g1|forward"));
    let protocol = fs::read_to_string(protocol_path).expect("read protocol export");
    assert!(protocol.contains("GENtle Guide Oligo Protocol"));
}

#[test]
fn test_guide_set_duplicate_ids_rejected() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let err = engine
        .apply(Operation::UpsertGuideSet {
            guide_set_id: "dup".to_string(),
            guides: vec![
                GuideCandidate {
                    guide_id: "dup_1".to_string(),
                    seq_id: "s1".to_string(),
                    start_0based: 0,
                    end_0based_exclusive: 20,
                    strand: "+".to_string(),
                    protospacer: "GACCTGTTGACGATGTTCCA".to_string(),
                    pam: "AGG".to_string(),
                    nuclease: "SpCas9".to_string(),
                    cut_offset_from_protospacer_start: 17,
                    rank: Some(1),
                },
                GuideCandidate {
                    guide_id: "dup_1".to_string(),
                    seq_id: "s1".to_string(),
                    start_0based: 20,
                    end_0based_exclusive: 40,
                    strand: "+".to_string(),
                    protospacer: "GACCTGTTGACGATGTTCCC".to_string(),
                    pam: "AGG".to_string(),
                    nuclease: "SpCas9".to_string(),
                    cut_offset_from_protospacer_start: 17,
                    rank: Some(2),
                },
            ],
        })
        .expect_err("duplicate guide ids should be rejected");
    assert!(err.message.contains("Duplicate guide_id"));
}

#[test]
fn test_compute_dotplot_stores_and_retrieves_view() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "s".to_string(),
        DNAsequence::from_sequence("ATGCATGCATGCATGC").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let result = engine
        .apply(Operation::ComputeDotplot {
            seq_id: "s".to_string(),
            reference_seq_id: None,
            span_start_0based: Some(0),
            span_end_0based: Some(16),
            reference_span_start_0based: None,
            reference_span_end_0based: None,
            mode: DotplotMode::SelfForward,
            word_size: 4,
            step_bp: 2,
            max_mismatches: 0,
            tile_bp: Some(128),
            store_as: Some("promoter_dotplot".to_string()),
        })
        .expect("compute dotplot");
    assert!(
        result
            .messages
            .iter()
            .any(|message| message.contains("promoter_dotplot")),
        "expected dotplot id in message, got: {:?}",
        result.messages
    );
    let rows = engine.list_dotplot_views(Some("s"));
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].dotplot_id, "promoter_dotplot");
    let view = engine
        .get_dotplot_view("promoter_dotplot")
        .expect("dotplot view");
    assert_eq!(view.schema, DOTPLOT_VIEW_SCHEMA);
    assert_eq!(view.seq_id, "s");
    assert!(view.reference_seq_id.is_none());
    assert_eq!(view.word_size, 4);
    assert!(!view.points.is_empty());
    assert_eq!(view.point_count, view.points.len());
    assert_eq!(view.boxplot_bin_count, view.boxplot_bins.len());
    assert!(view.boxplot_bin_count > 0);
    assert!(view.boxplot_bins.iter().any(|bin| bin.hit_count > 0));
}

#[test]
fn test_compute_pair_dotplot_uses_reference_sequence_span() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "query".to_string(),
        DNAsequence::from_sequence("ATGCATGCATGC").expect("query"),
    );
    state.sequences.insert(
        "ref".to_string(),
        DNAsequence::from_sequence("TTTATGCATGCAAAA").expect("ref"),
    );
    let mut engine = GentleEngine::from_state(state);
    let result = engine
        .apply(Operation::ComputeDotplot {
            seq_id: "query".to_string(),
            reference_seq_id: Some("ref".to_string()),
            span_start_0based: Some(0),
            span_end_0based: Some(12),
            reference_span_start_0based: Some(3),
            reference_span_end_0based: Some(15),
            mode: DotplotMode::PairForward,
            word_size: 4,
            step_bp: 2,
            max_mismatches: 0,
            tile_bp: None,
            store_as: Some("pair_dotplot".to_string()),
        })
        .expect("compute pair dotplot");
    assert!(
        result
            .messages
            .iter()
            .any(|message| message.contains("pair_dotplot"))
    );
    let view = engine
        .get_dotplot_view("pair_dotplot")
        .expect("pair dotplot view");
    assert_eq!(view.seq_id, "query");
    assert_eq!(view.reference_seq_id.as_deref(), Some("ref"));
    assert_eq!(view.span_start_0based, 0);
    assert_eq!(view.span_end_0based, 12);
    assert_eq!(view.reference_span_start_0based, 3);
    assert_eq!(view.reference_span_end_0based, 15);
    assert_eq!(view.mode, DotplotMode::PairForward);
    assert!(!view.points.is_empty());
    assert_eq!(view.boxplot_bin_count, view.boxplot_bins.len());
    assert!(view.boxplot_bins.iter().any(|bin| bin.hit_count > 0));
}

#[test]
fn test_preview_pair_dotplot_view_builds_transient_view() {
    let view = GentleEngine::preview_pair_dotplot_view(
        "read_1",
        "ACGTACGT",
        "roi_1",
        "TTTACGTACGTTTT",
        3,
        11,
        DotplotMode::PairForward,
        4,
        1,
        0,
        None,
    )
    .expect("preview pair dotplot view");
    assert_eq!(view.schema, DOTPLOT_VIEW_SCHEMA);
    assert_eq!(view.dotplot_id, "");
    assert_eq!(view.seq_id, "read_1");
    assert_eq!(view.reference_seq_id.as_deref(), Some("roi_1"));
    assert_eq!(view.span_start_0based, 0);
    assert_eq!(view.span_end_0based, 8);
    assert_eq!(view.reference_span_start_0based, 3);
    assert_eq!(view.reference_span_end_0based, 11);
    assert_eq!(view.word_size, 4);
    assert_eq!(view.step_bp, 1);
    assert!(!view.points.is_empty());
}

#[test]
fn test_compute_pair_reverse_complement_maps_expected_antidiagonal_hits() {
    let query = "ACGTTGCAAGTC";
    let reference = GentleEngine::reverse_complement(query);

    let mut state = ProjectState::default();
    state.sequences.insert(
        "query".to_string(),
        DNAsequence::from_sequence(query).expect("query"),
    );
    state.sequences.insert(
        "ref".to_string(),
        DNAsequence::from_sequence(&reference).expect("reference"),
    );
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::ComputeDotplot {
            seq_id: "query".to_string(),
            reference_seq_id: Some("ref".to_string()),
            span_start_0based: Some(0),
            span_end_0based: Some(query.len()),
            reference_span_start_0based: Some(0),
            reference_span_end_0based: Some(reference.len()),
            mode: DotplotMode::PairReverseComplement,
            word_size: 4,
            step_bp: 1,
            max_mismatches: 0,
            tile_bp: None,
            store_as: Some("pair_revcomp".to_string()),
        })
        .expect("compute pair reverse-complement dotplot");

    let view = engine
        .get_dotplot_view("pair_revcomp")
        .expect("pair reverse-complement view");
    assert_eq!(view.mode, DotplotMode::PairReverseComplement);
    let expected_windows = query.len().saturating_sub(4).saturating_add(1);
    assert_eq!(view.point_count, expected_windows);
    for x in 0..expected_windows {
        let expected_y = query.len().saturating_sub(4).saturating_sub(x);
        assert!(
            view.points
                .iter()
                .any(|point| point.x_0based == x && point.y_0based == expected_y),
            "missing expected anti-diagonal hit for x={x}, y={expected_y}"
        );
    }
}

#[test]
fn test_compute_dotplot_exact_seed_large_pair_grid_bypasses_pair_eval_guard() {
    let len = 12_500usize;
    let mut seed = 0x9e37_79b9_7f4a_7c15u64;
    let mut query = String::with_capacity(len);
    for _ in 0..len {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let base = match (seed & 0b11) as u8 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        };
        query.push(base);
    }

    let mut state = ProjectState::default();
    state.sequences.insert(
        "query".to_string(),
        DNAsequence::from_sequence(&query).expect("query"),
    );
    state.sequences.insert(
        "ref".to_string(),
        DNAsequence::from_sequence(&query).expect("reference"),
    );
    let mut engine = GentleEngine::from_state(state);

    let result = engine.apply(Operation::ComputeDotplot {
        seq_id: "query".to_string(),
        reference_seq_id: Some("ref".to_string()),
        span_start_0based: Some(0),
        span_end_0based: Some(len),
        reference_span_start_0based: Some(0),
        reference_span_end_0based: Some(len),
        mode: DotplotMode::PairForward,
        word_size: 10,
        step_bp: 1,
        max_mismatches: 0,
        tile_bp: None,
        store_as: Some("pair_exact_large".to_string()),
    });
    assert!(
        result.is_ok(),
        "exact-seed mode should bypass brute-force pair-evaluation guard"
    );
    let view = engine
        .get_dotplot_view("pair_exact_large")
        .expect("pair exact large view");
    assert!(
        view.point_count > 0,
        "expected deterministic exact-seed hits"
    );
}

#[test]
fn test_compute_flexibility_track_stores_and_retrieves_track() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "s".to_string(),
        DNAsequence::from_sequence("AAAATTTTCCCCGGGGAAAATTTT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let result = engine
        .apply(Operation::ComputeFlexibilityTrack {
            seq_id: "s".to_string(),
            span_start_0based: Some(0),
            span_end_0based: Some(24),
            model: FlexibilityModel::AtRichness,
            bin_bp: 4,
            smoothing_bp: Some(8),
            store_as: Some("promoter_flex".to_string()),
        })
        .expect("compute flexibility track");
    assert!(
        result
            .messages
            .iter()
            .any(|message| message.contains("promoter_flex")),
        "expected track id in message, got: {:?}",
        result.messages
    );
    let rows = engine.list_flexibility_tracks(Some("s"));
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].track_id, "promoter_flex");
    let track = engine
        .get_flexibility_track("promoter_flex")
        .expect("flexibility track");
    assert_eq!(track.schema, FLEXIBILITY_TRACK_SCHEMA);
    assert_eq!(track.model, FlexibilityModel::AtRichness);
    assert_eq!(track.bin_bp, 4);
    assert_eq!(track.bins.len(), 6);
    assert!(track.max_score >= track.min_score);
}

#[test]
fn test_derive_splicing_references_from_window_without_seed_feature() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);

    let result = engine
        .apply(Operation::DeriveSplicingReferences {
            seq_id: "seq_a".to_string(),
            span_start_0based: 1,
            span_end_0based: 34,
            seed_feature_id: None,
            scope: SplicingScopePreset::TargetGroupTargetStrand,
            output_prefix: Some("seq_a_splicing_refs".to_string()),
        })
        .expect("derive splicing references");

    assert_eq!(result.created_seq_ids.len(), 4);
    let dna_window_id = result
        .created_seq_ids
        .iter()
        .find(|seq_id| seq_id.ends_with("_dna"))
        .expect("dna window id");
    let dna_window = engine
        .state()
        .sequences
        .get(dna_window_id)
        .expect("dna window sequence");
    assert_eq!(dna_window.len(), 33);
    assert!(!dna_window.is_circular());

    let mrna_ids = result
        .created_seq_ids
        .iter()
        .filter(|seq_id| seq_id.contains("_mrna_"))
        .cloned()
        .collect::<Vec<_>>();
    assert_eq!(mrna_ids.len(), 2);
    let view = engine
        .build_splicing_expert_view("seq_a", 0, SplicingScopePreset::TargetGroupTargetStrand)
        .expect("splicing view");
    let source_dna = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("source sequence");
    let mut mrna_lengths = vec![];
    for mrna_id in &mrna_ids {
        let mrna = engine
            .state()
            .sequences
            .get(mrna_id)
            .expect("mRNA sequence present");
        assert!(!mrna.is_circular());
        mrna_lengths.push(mrna.len());
    }
    let mut expected_mrna_lengths = view
        .transcripts
        .iter()
        .map(|lane| {
            GentleEngine::make_transcript_template(source_dna, lane, 0)
                .sequence
                .len()
        })
        .collect::<Vec<_>>();
    mrna_lengths.sort_unstable();
    expected_mrna_lengths.sort_unstable();
    assert_eq!(mrna_lengths, expected_mrna_lengths);

    let exon_reference_id = result
        .created_seq_ids
        .iter()
        .find(|seq_id| seq_id.ends_with("_exon_reference"))
        .expect("exon reference id");
    let exon_reference = engine
        .state()
        .sequences
        .get(exon_reference_id)
        .expect("exon reference sequence");
    let expected_exon_ref_len = view
        .unique_exons
        .iter()
        .map(|exon| exon.end_1based + 1 - exon.start_1based)
        .sum::<usize>();
    assert_eq!(exon_reference.len(), expected_exon_ref_len);
    assert!(!exon_reference.get_forward_string().contains('U'));
}

#[test]
fn test_build_splicing_expert_view_accepts_gene_seed_feature() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_seed".to_string(), splicing_seed_feature_sequence());
    let engine = GentleEngine::from_state(state);

    let view = engine
        .build_splicing_expert_view("seq_seed", 2, SplicingScopePreset::TargetGroupTargetStrand)
        .expect("splicing view from gene seed");

    assert_eq!(view.target_feature_id, 2);
    assert_eq!(view.group_label, "GENE1");
    assert_eq!(view.transcript_count, 2);
    assert!(!view.transcripts.is_empty());
}

#[test]
fn test_build_splicing_expert_view_accepts_cds_seed_feature() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_seed".to_string(), splicing_seed_feature_sequence());
    let engine = GentleEngine::from_state(state);

    let view = engine
        .build_splicing_expert_view("seq_seed", 3, SplicingScopePreset::TargetGroupTargetStrand)
        .expect("splicing view from cds seed");

    assert_eq!(view.target_feature_id, 3);
    assert_eq!(view.group_label, "GENE1");
    assert_eq!(view.transcript_count, 2);
    assert!(!view.unique_exons.is_empty());
}

#[test]
fn test_build_splicing_expert_view_accepts_misc_rna_seed_feature() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_misc".to_string(), splicing_misc_rna_sequence());
    let engine = GentleEngine::from_state(state);

    let view = engine
        .build_splicing_expert_view("seq_misc", 0, SplicingScopePreset::TargetGroupTargetStrand)
        .expect("splicing view from misc_RNA seed");

    assert_eq!(view.target_feature_id, 0);
    assert_eq!(view.group_label, "GENE1");
    assert_eq!(view.transcript_count, 2);
    assert!(!view.transcripts.is_empty());
}

#[test]
fn test_align_sequences_global_sets_structured_result() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "query".to_string(),
        DNAsequence::from_sequence("ACGTTT").expect("query"),
    );
    state.sequences.insert(
        "target".to_string(),
        DNAsequence::from_sequence("ACGTAT").expect("target"),
    );
    let mut engine = GentleEngine::from_state(state);
    let result = engine
        .apply(Operation::AlignSequences {
            query_seq_id: "query".to_string(),
            target_seq_id: "target".to_string(),
            query_span_start_0based: None,
            query_span_end_0based: None,
            target_span_start_0based: None,
            target_span_end_0based: None,
            mode: PairwiseAlignmentMode::Global,
            match_score: 2,
            mismatch_score: -1,
            gap_open: -5,
            gap_extend: -1,
        })
        .expect("align sequences globally");

    assert!(result.created_seq_ids.is_empty());
    let report = result
        .sequence_alignment
        .expect("structured sequence alignment report");
    assert_eq!(report.schema, SEQUENCE_ALIGNMENT_REPORT_SCHEMA);
    assert_eq!(report.mode, PairwiseAlignmentMode::Global);
    assert_eq!(report.query_seq_id, "query");
    assert_eq!(report.target_seq_id, "target");
    assert_eq!(report.matches, 5);
    assert_eq!(report.mismatches, 1);
    assert_eq!(report.insertions, 0);
    assert_eq!(report.deletions, 0);
    assert_eq!(report.score, 9);
    assert_eq!(report.cigar, "4=1X1=");
    assert!((report.identity_fraction - (5.0 / 6.0)).abs() < 1e-9);
    assert!((report.query_coverage_fraction - 1.0).abs() < 1e-9);
    assert!((report.target_coverage_fraction - 1.0).abs() < 1e-9);
}

#[test]
fn test_align_sequences_local_reports_partial_coverage() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "query".to_string(),
        DNAsequence::from_sequence("TTTACGTAA").expect("query"),
    );
    state.sequences.insert(
        "target".to_string(),
        DNAsequence::from_sequence("GGGACGTCCC").expect("target"),
    );
    let mut engine = GentleEngine::from_state(state);
    let result = engine
        .apply(Operation::AlignSequences {
            query_seq_id: "query".to_string(),
            target_seq_id: "target".to_string(),
            query_span_start_0based: None,
            query_span_end_0based: None,
            target_span_start_0based: None,
            target_span_end_0based: None,
            mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
        })
        .expect("align sequences locally");

    let report = result
        .sequence_alignment
        .expect("structured local alignment report");
    assert_eq!(report.mode, PairwiseAlignmentMode::Local);
    assert!(report.matches >= 4);
    assert!(report.aligned_columns >= 4);
    assert!(report.query_coverage_fraction < 1.0);
    assert!(report.target_coverage_fraction < 1.0);
}

#[test]
fn test_confirm_construct_reads_confirms_junction_target_and_lists_report() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    state.sequences.insert(
        "read_junction".to_string(),
        DNAsequence::from_sequence("CCGTAACC").expect("read"),
    );
    let mut engine = GentleEngine::from_state(state);
    let result = engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_junction".to_string()],
            trace_ids: vec![],
            targets: vec![sequencing_confirmation_junction_target(
                "junction_1",
                4,
                12,
                8,
                "Insert junction",
            )],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("construct_confirm".to_string()),
        })
        .expect("confirm construct");

    let report = result
        .sequencing_confirmation_report
        .expect("sequencing confirmation report");
    assert_eq!(report.report_id, "construct_confirm");
    assert_eq!(
        report.overall_status,
        SequencingConfirmationStatus::Confirmed
    );
    assert_eq!(report.targets.len(), 1);
    assert_eq!(
        report.targets[0].status,
        SequencingConfirmationStatus::Confirmed
    );
    assert_eq!(report.targets[0].support_read_ids, vec!["read_junction"]);
    assert!(report.trace_ids.is_empty());
    assert_eq!(report.reads.len(), 1);
    assert_eq!(
        report.reads[0].evidence_kind,
        SequencingConfirmationEvidenceKind::Sequence
    );
    assert_eq!(report.reads[0].evidence_id, "read_junction");
    assert_eq!(
        report.reads[0].orientation,
        SequencingReadOrientation::Forward
    );
    assert_eq!(report.reads[0].confirmed_target_ids, vec!["junction_1"]);

    let fetched = engine
        .get_sequencing_confirmation_report("construct_confirm")
        .expect("stored report");
    assert_eq!(
        fetched.overall_status,
        SequencingConfirmationStatus::Confirmed
    );
    let summaries = engine.list_sequencing_confirmation_reports(Some("construct"));
    assert_eq!(summaries.len(), 1);
    assert_eq!(summaries[0].report_id, "construct_confirm");
    assert_eq!(
        summaries[0].overall_status,
        SequencingConfirmationStatus::Confirmed
    );
}

#[test]
fn test_confirm_construct_reads_supports_reverse_complement_reads() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    state.sequences.insert(
        "read_rc".to_string(),
        DNAsequence::from_sequence("GGTTACGG").expect("reverse-complement read"),
    );
    let mut engine = GentleEngine::from_state(state);
    let result = engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_rc".to_string()],
            trace_ids: vec![],
            targets: vec![sequencing_confirmation_junction_target(
                "junction_1",
                4,
                12,
                8,
                "Insert junction",
            )],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("construct_rc".to_string()),
        })
        .expect("confirm construct from reverse-complement read");

    let report = result
        .sequencing_confirmation_report
        .expect("sequencing confirmation report");
    assert_eq!(
        report.overall_status,
        SequencingConfirmationStatus::Confirmed
    );
    assert_eq!(
        report.reads[0].orientation,
        SequencingReadOrientation::ReverseComplement
    );
    assert_eq!(
        report.targets[0].status,
        SequencingConfirmationStatus::Confirmed
    );
}

#[test]
fn test_confirm_construct_reads_reports_insufficient_evidence_for_truncated_read() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    state.sequences.insert(
        "read_short".to_string(),
        DNAsequence::from_sequence("CCGT").expect("short read"),
    );
    let mut engine = GentleEngine::from_state(state);
    let result = engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_short".to_string()],
            trace_ids: vec![],
            targets: vec![sequencing_confirmation_junction_target(
                "junction_1",
                4,
                12,
                8,
                "Insert junction",
            )],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("construct_short".to_string()),
        })
        .expect("confirm construct from truncated read");

    let report = result
        .sequencing_confirmation_report
        .expect("sequencing confirmation report");
    assert_eq!(
        report.overall_status,
        SequencingConfirmationStatus::InsufficientEvidence
    );
    assert_eq!(
        report.targets[0].status,
        SequencingConfirmationStatus::InsufficientEvidence
    );
    assert!(report.targets[0].support_read_ids.is_empty());
    assert!(report.targets[0].contradicting_read_ids.is_empty());
    assert_eq!(report.reads[0].covered_target_ids, vec!["junction_1"]);
    assert!(report.targets[0].covered_bp < report.targets[0].target_length_bp);
}

#[test]
fn test_export_sequencing_confirmation_support_tsv_writes_target_rows() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    state.sequences.insert(
        "read_junction".to_string(),
        DNAsequence::from_sequence("CCGTAACC").expect("read"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_junction".to_string()],
            trace_ids: vec![],
            targets: vec![sequencing_confirmation_junction_target(
                "junction_1",
                4,
                12,
                8,
                "Insert junction",
            )],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("construct_tsv".to_string()),
        })
        .expect("confirm construct");

    let td = tempdir().expect("tempdir");
    let output = td.path().join("support.tsv");
    engine
        .export_sequencing_confirmation_support_tsv(
            "construct_tsv",
            output.to_str().expect("utf-8 path"),
        )
        .expect("export sequencing-confirmation TSV");
    let text = fs::read_to_string(output).expect("read TSV");
    assert!(text.contains("report_id\texpected_seq_id\toverall_status"));
    assert!(text.contains("construct_tsv\tconstruct\tconfirmed"));
    assert!(text.contains("junction_1"));
}

#[test]
fn test_import_sequencing_trace_stores_record_without_mutating_sequences() {
    let fixture = sequencing_confirmation_fixture_path("3100.ab1");
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    let original_sequence_count = state.sequences.len();
    let mut engine = GentleEngine::from_state(state);

    let result = engine
        .apply(Operation::ImportSequencingTrace {
            path: fixture.clone(),
            trace_id: Some("abi_trace".to_string()),
            seq_id: Some("construct".to_string()),
        })
        .expect("import ABI trace");

    assert_eq!(engine.state().sequences.len(), original_sequence_count);
    assert!(result.created_seq_ids.is_empty());
    assert!(result.changed_seq_ids.is_empty());
    let import_report = result
        .sequencing_trace_import_report
        .expect("trace import report");
    let record = result.sequencing_trace_record.expect("trace record");
    assert_eq!(import_report.trace_id, "abi_trace");
    assert_eq!(record.trace_id, "abi_trace");
    assert_eq!(record.format, SequencingTraceFormat::AbiAb1);
    assert_eq!(record.seq_id.as_deref(), Some("construct"));
    assert!(!record.called_bases.is_empty());
    assert_eq!(import_report.called_base_count, record.called_bases.len());
    assert_eq!(
        record.called_base_confidence_values.len(),
        record.called_bases.len()
    );
    assert_eq!(record.peak_locations.len(), record.called_bases.len());
    assert!(!record.channel_data.is_empty());
    assert_eq!(record.channel_data.len(), 4);
    assert!(record.channel_data.iter().all(|row| !row.points.is_empty()));
    assert_eq!(record.channel_summaries.len(), 4);
    assert!(import_report.has_curve_data);
    assert_eq!(import_report.channel_count, 4);
    assert!(
        record
            .sample_name
            .as_deref()
            .is_some_and(|value| !value.trim().is_empty())
    );

    let listed = engine.list_sequencing_traces(Some("construct"));
    assert_eq!(listed.len(), 1);
    assert_eq!(listed[0].trace_id, "abi_trace");
    assert!(listed[0].has_curve_data);
    assert_eq!(listed[0].channel_count, 4);
    let shown = engine
        .get_sequencing_trace("abi_trace")
        .expect("stored trace");
    assert_eq!(shown.called_bases, record.called_bases);
    assert_eq!(shown.channel_data.len(), 4);
}

#[test]
fn test_confirm_construct_reads_accepts_imported_trace_evidence() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .upsert_sequencing_trace(SequencingTraceRecord {
            schema: SEQUENCING_TRACE_RECORD_SCHEMA.to_string(),
            trace_id: "junction_trace".to_string(),
            format: SequencingTraceFormat::Scf,
            source_path: "synthetic.scf".to_string(),
            imported_at_unix_ms: 1,
            seq_id: None,
            sample_name: Some("junction_trace".to_string()),
            sample_well: None,
            run_name: None,
            machine_name: None,
            machine_model: None,
            called_bases: "CCGTAACC".to_string(),
            called_base_confidence_values: vec![80; 8],
            peak_locations: (1..=8).collect(),
            channel_data: vec![],
            channel_summaries: vec![],
            clip_start_base_index: None,
            clip_end_base_index_exclusive: None,
            comments_text: None,
        })
        .expect("store trace");

    let result = engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec![],
            trace_ids: vec!["junction_trace".to_string()],
            targets: vec![sequencing_confirmation_junction_target(
                "junction_1",
                4,
                12,
                8,
                "Insert junction",
            )],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("construct_trace".to_string()),
        })
        .expect("confirm construct from trace");

    let report = result
        .sequencing_confirmation_report
        .expect("sequencing confirmation report");
    assert_eq!(report.read_seq_ids, Vec::<String>::new());
    assert_eq!(report.trace_ids, vec!["junction_trace"]);
    assert_eq!(
        report.overall_status,
        SequencingConfirmationStatus::Confirmed
    );
    assert_eq!(report.targets[0].support_read_ids, vec!["junction_trace"]);
    assert_eq!(
        report.reads[0].evidence_kind,
        SequencingConfirmationEvidenceKind::Trace
    );
    assert_eq!(report.reads[0].evidence_id, "junction_trace");
    assert_eq!(report.reads[0].trace_id.as_deref(), Some("junction_trace"));
    assert_eq!(report.reads[0].read_seq_id, "junction_trace");
    assert_eq!(report.reads[0].confirmed_target_ids, vec!["junction_1"]);
}

#[test]
fn test_import_sequencing_trace_rejects_malformed_abi_fixture() {
    let fixture = sequencing_confirmation_fixture_path("fake.ab1");
    let mut engine = GentleEngine::default();
    let err = engine
        .apply(Operation::ImportSequencingTrace {
            path: fixture,
            trace_id: Some("bad_trace".to_string()),
            seq_id: None,
        })
        .expect_err("malformed ABI should fail");
    assert_eq!(err.code, ErrorCode::InvalidInput);
    assert!(
        err.message.contains("Unsupported sequencing trace format") || err.message.contains("ABIF"),
        "unexpected error: {}",
        err.message
    );
}

#[test]
fn test_list_and_show_sequencing_traces_use_persisted_store() {
    let fixture = sequencing_confirmation_fixture_path("3100.ab1");
    let mut engine = GentleEngine::default();
    engine
        .apply(Operation::ImportSequencingTrace {
            path: fixture.clone(),
            trace_id: Some("abi_trace".to_string()),
            seq_id: None,
        })
        .expect("import ABI trace");
    engine
        .apply(Operation::ImportSequencingTrace {
            path: fixture,
            trace_id: Some("abi_trace".to_string()),
            seq_id: None,
        })
        .expect("repeat import ABI trace");

    let list_result = engine
        .apply(Operation::ListSequencingTraces { seq_id: None })
        .expect("list sequencing traces");
    let rows = list_result
        .sequencing_trace_summaries
        .expect("trace summary rows");
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].trace_id, "abi_trace");

    let show_result = engine
        .apply(Operation::ShowSequencingTrace {
            trace_id: "abi_trace".to_string(),
        })
        .expect("show sequencing trace");
    let shown = show_result
        .sequencing_trace_record
        .expect("shown sequencing trace");
    assert_eq!(shown.trace_id, "abi_trace");
    assert_eq!(shown.format, SequencingTraceFormat::AbiAb1);
    assert!(!shown.called_bases.is_empty());
}

#[test]
fn test_old_v1_sequencing_trace_records_without_curve_data_remain_readable() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    state.metadata.insert(
        SEQUENCING_TRACES_METADATA_KEY.to_string(),
        serde_json::json!({
            "schema": SEQUENCING_TRACES_SCHEMA,
            "updated_at_unix_ms": 1u128,
            "traces": {
                "legacy_trace": {
                    "schema": "gentle.sequencing_trace_record.v1",
                    "trace_id": "legacy_trace",
                    "format": "scf",
                    "source_path": "legacy.scf",
                    "imported_at_unix_ms": 1u128,
                    "called_bases": "CCGTAACC",
                    "called_base_confidence_values": [80, 80, 80, 80, 80, 80, 80, 80],
                    "peak_locations": [10, 20, 30, 40, 50, 60, 70, 80],
                    "channel_summaries": []
                }
            }
        }),
    );
    let mut engine = GentleEngine::from_state(state);

    let legacy = engine
        .get_sequencing_trace("legacy_trace")
        .expect("legacy trace should deserialize");
    assert_eq!(legacy.schema, "gentle.sequencing_trace_record.v1");
    assert!(legacy.channel_data.is_empty());

    let result = engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec![],
            trace_ids: vec!["legacy_trace".to_string()],
            targets: vec![sequencing_confirmation_junction_target(
                "junction_1",
                4,
                12,
                8,
                "Insert junction",
            )],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("legacy_trace_confirm".to_string()),
        })
        .expect("legacy trace should still confirm construct");

    let report = result
        .sequencing_confirmation_report
        .expect("sequencing confirmation report");
    assert_eq!(
        report.overall_status,
        SequencingConfirmationStatus::Confirmed
    );
    assert_eq!(report.trace_ids, vec!["legacy_trace"]);
    assert!(report.variants.is_empty());
}

#[test]
fn test_confirm_construct_reads_baseline_infers_intended_edit_variant_for_trace() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "expected".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("expected"),
    );
    state.sequences.insert(
        "baseline".to_string(),
        DNAsequence::from_sequence("AAAACCGTGACCTTTT").expect("baseline"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .upsert_sequencing_trace(sequencing_confirmation_trace_record(
            "trace_expected_snp",
            "CCGTAACC",
            &[80, 80, 80, 80, 80, 80, 80, 80],
            &[10, 20, 30, 40, 50, 60, 70, 80],
        ))
        .expect("store trace");

    let result = engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "expected".to_string(),
            baseline_seq_id: Some("baseline".to_string()),
            read_seq_ids: vec![],
            trace_ids: vec!["trace_expected_snp".to_string()],
            targets: vec![],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("expected_edit_confirmed".to_string()),
        })
        .expect("confirm intended edit");

    let report = result
        .sequencing_confirmation_report
        .expect("sequencing confirmation report");
    assert_eq!(report.baseline_seq_id.as_deref(), Some("baseline"));
    assert_eq!(
        report.overall_status,
        SequencingConfirmationStatus::Confirmed
    );
    assert_eq!(report.targets.len(), 1);
    assert_eq!(
        report.targets[0].kind,
        SequencingConfirmationTargetKind::ExpectedEdit
    );
    assert_eq!(
        report.targets[0].status,
        SequencingConfirmationStatus::Confirmed
    );
    assert_eq!(report.variants.len(), 1);
    assert_eq!(
        report.variants[0].classification,
        SequencingConfirmationVariantClassification::IntendedEditConfirmed
    );
    assert_eq!(
        report.variants[0].status,
        SequencingConfirmationStatus::Confirmed
    );
    assert_eq!(
        report.variants[0].trace_id.as_deref(),
        Some("trace_expected_snp")
    );
}

#[test]
fn test_confirm_construct_reads_baseline_detects_reference_reversion_for_trace() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "expected".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("expected"),
    );
    state.sequences.insert(
        "baseline".to_string(),
        DNAsequence::from_sequence("AAAACCGTGACCTTTT").expect("baseline"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .upsert_sequencing_trace(sequencing_confirmation_trace_record(
            "trace_baseline_snp",
            "CCGTGACC",
            &[80, 80, 80, 80, 80, 80, 80, 80],
            &[10, 20, 30, 40, 50, 60, 70, 80],
        ))
        .expect("store trace");

    let report = engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "expected".to_string(),
            baseline_seq_id: Some("baseline".to_string()),
            read_seq_ids: vec![],
            trace_ids: vec!["trace_baseline_snp".to_string()],
            targets: vec![],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("reference_reversion".to_string()),
        })
        .expect("confirm reference reversion")
        .sequencing_confirmation_report
        .expect("sequencing confirmation report");

    assert_eq!(
        report.overall_status,
        SequencingConfirmationStatus::Contradicted
    );
    assert_eq!(report.targets.len(), 1);
    assert_eq!(
        report.targets[0].status,
        SequencingConfirmationStatus::Contradicted
    );
    assert_eq!(
        report.variants[0].classification,
        SequencingConfirmationVariantClassification::ReferenceReversion
    );
}

#[test]
fn test_confirm_construct_reads_baseline_detects_unexpected_difference_for_trace() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "expected".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("expected"),
    );
    state.sequences.insert(
        "baseline".to_string(),
        DNAsequence::from_sequence("AAAACCGTGACCTTTT").expect("baseline"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .upsert_sequencing_trace(sequencing_confirmation_trace_record(
            "trace_unexpected_snp",
            "CCGTTACC",
            &[80, 80, 80, 80, 80, 80, 80, 80],
            &[10, 20, 30, 40, 50, 60, 70, 80],
        ))
        .expect("store trace");

    let report = engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "expected".to_string(),
            baseline_seq_id: Some("baseline".to_string()),
            read_seq_ids: vec![],
            trace_ids: vec!["trace_unexpected_snp".to_string()],
            targets: vec![],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("unexpected_difference".to_string()),
        })
        .expect("confirm unexpected difference")
        .sequencing_confirmation_report
        .expect("sequencing confirmation report");

    assert_eq!(
        report.overall_status,
        SequencingConfirmationStatus::Contradicted
    );
    assert_eq!(
        report.variants[0].classification,
        SequencingConfirmationVariantClassification::UnexpectedDifference
    );
}

#[test]
fn test_confirm_construct_reads_baseline_insertion_counts_as_intended_edit_confirmation() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "expected".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("expected"),
    );
    state.sequences.insert(
        "baseline".to_string(),
        DNAsequence::from_sequence("AAAACCGTACCTTTT").expect("baseline"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .upsert_sequencing_trace(sequencing_confirmation_trace_record(
            "trace_expected_insertion",
            "AAAACCGTAACCTTTT",
            &[80; 16],
            &[
                10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160,
            ],
        ))
        .expect("store trace");

    let report = engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "expected".to_string(),
            baseline_seq_id: Some("baseline".to_string()),
            read_seq_ids: vec![],
            trace_ids: vec!["trace_expected_insertion".to_string()],
            targets: vec![],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("expected_insertion".to_string()),
        })
        .expect("confirm intended insertion")
        .sequencing_confirmation_report
        .expect("sequencing confirmation report");

    assert_eq!(
        report.overall_status,
        SequencingConfirmationStatus::Confirmed
    );
    assert_eq!(
        report.variants[0].classification,
        SequencingConfirmationVariantClassification::IntendedEditConfirmed
    );
    assert_eq!(
        report.targets[0].kind,
        SequencingConfirmationTargetKind::ExpectedEdit
    );
    assert_eq!(
        report.targets[0].status,
        SequencingConfirmationStatus::Confirmed
    );
}

#[test]
fn test_confirm_construct_reads_low_confidence_trace_softens_reversion_to_insufficient_evidence() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "expected".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("expected"),
    );
    state.sequences.insert(
        "baseline".to_string(),
        DNAsequence::from_sequence("AAAACCGTGACCTTTT").expect("baseline"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .upsert_sequencing_trace(sequencing_confirmation_trace_record(
            "trace_low_confidence_reversion",
            "CCGTGACC",
            &[10, 10, 10, 10, 10, 10, 10, 10],
            &[10, 20, 30, 40, 50, 60, 70, 80],
        ))
        .expect("store trace");

    let report = engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "expected".to_string(),
            baseline_seq_id: Some("baseline".to_string()),
            read_seq_ids: vec![],
            trace_ids: vec!["trace_low_confidence_reversion".to_string()],
            targets: vec![],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("low_confidence_reversion".to_string()),
        })
        .expect("confirm low-confidence reversion")
        .sequencing_confirmation_report
        .expect("sequencing confirmation report");

    assert_eq!(
        report.overall_status,
        SequencingConfirmationStatus::InsufficientEvidence
    );
    assert_eq!(
        report.targets[0].status,
        SequencingConfirmationStatus::InsufficientEvidence
    );
    assert_eq!(
        report.variants[0].classification,
        SequencingConfirmationVariantClassification::LowConfidenceOrAmbiguous
    );
    assert_eq!(
        report.variants[0].status,
        SequencingConfirmationStatus::InsufficientEvidence
    );
}

#[test]
fn test_suggest_sequencing_primers_reports_forward_and_reverse_overlay_hits() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    state.sequences.insert(
        "primer_fwd".to_string(),
        DNAsequence::from_sequence("TTTACCGT").expect("primer"),
    );
    state.sequences.insert(
        "primer_rev".to_string(),
        DNAsequence::from_sequence("TTTGGTT").expect("primer"),
    );
    let mut engine = GentleEngine::from_state(state);

    let report = engine
        .apply(Operation::SuggestSequencingPrimers {
            expected_seq_id: "construct".to_string(),
            primer_seq_ids: vec!["primer_fwd".to_string(), "primer_rev".to_string()],
            confirmation_report_id: None,
            min_3prime_anneal_bp: 4,
            predicted_read_length_bp: 10,
        })
        .expect("suggest primers")
        .sequencing_primer_overlay_report
        .expect("overlay report");

    assert_eq!(report.expected_seq_id, "construct");
    assert!(report.suggestion_count >= 2);
    assert!(report.suggestions.iter().any(|row| row.orientation
        == SequencingPrimerOrientation::ForwardRead
        && row.primer_seq_id == "primer_fwd"
        && row.predicted_read_span_start_0based == 4
        && row.predicted_read_span_end_0based_exclusive == 14));
    assert!(report.suggestions.iter().any(|row| row.orientation
        == SequencingPrimerOrientation::ReverseRead
        && row.primer_seq_id == "primer_rev"
        && row.predicted_read_span_start_0based == 2
        && row.predicted_read_span_end_0based_exclusive == 12));
}

#[test]
fn test_suggest_sequencing_primers_annotates_confirmation_report_coverage() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    state.sequences.insert(
        "read_junction".to_string(),
        DNAsequence::from_sequence("CCGTAACC").expect("read"),
    );
    state.sequences.insert(
        "primer_fwd".to_string(),
        DNAsequence::from_sequence("TTTACCGT").expect("primer"),
    );
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_junction".to_string()],
            trace_ids: vec![],
            targets: vec![sequencing_confirmation_junction_target(
                "junction_1",
                4,
                12,
                8,
                "Insert junction",
            )],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("construct_confirm".to_string()),
        })
        .expect("persist sequencing confirmation report");

    let report = engine
        .apply(Operation::SuggestSequencingPrimers {
            expected_seq_id: "construct".to_string(),
            primer_seq_ids: vec!["primer_fwd".to_string()],
            confirmation_report_id: Some("construct_confirm".to_string()),
            min_3prime_anneal_bp: 4,
            predicted_read_length_bp: 10,
        })
        .expect("suggest primers with report coverage")
        .sequencing_primer_overlay_report
        .expect("overlay report");

    assert_eq!(
        report.confirmation_report_id.as_deref(),
        Some("construct_confirm")
    );
    assert!(report.suggestion_count >= 1);
    assert!(report.suggestions.iter().any(|row| {
        row.primer_seq_id == "primer_fwd"
            && row.covered_target_ids == vec!["junction_1"]
            && row.covered_problem_target_ids.is_empty()
            && row.covered_variant_ids.is_empty()
            && row.covered_problem_variant_ids.is_empty()
    }));
}

#[test]
fn test_suggest_sequencing_primers_recommends_best_existing_primer_for_unresolved_target() {
    let mut state = ProjectState::default();
    let construct_text = format!(
        "{}ACCGTA{}TTGCAA{}",
        "G".repeat(20),
        "C".repeat(70),
        "A".repeat(120)
    );
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence(&construct_text).expect("construct"),
    );
    state.sequences.insert(
        "read_early".to_string(),
        DNAsequence::from_sequence("GGGGGGGGGGGG").expect("read"),
    );
    state.sequences.insert(
        "primer_good".to_string(),
        DNAsequence::from_sequence("TTTACCGTA").expect("primer"),
    );
    state.sequences.insert(
        "primer_near".to_string(),
        DNAsequence::from_sequence("TTTTTGCAA").expect("primer"),
    );
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_early".to_string()],
            trace_ids: vec![],
            targets: vec![SequencingConfirmationTargetSpec {
                target_id: "gap_target".to_string(),
                label: "Gap locus".to_string(),
                kind: SequencingConfirmationTargetKind::Junction,
                start_0based: 118,
                end_0based_exclusive: 122,
                junction_left_end_0based: Some(120),
                expected_bases: None,
                baseline_bases: None,
                required: true,
            }],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("unresolved_gap".to_string()),
        })
        .expect("persist unresolved report");

    let report = engine
        .apply(Operation::SuggestSequencingPrimers {
            expected_seq_id: "construct".to_string(),
            primer_seq_ids: vec!["primer_good".to_string(), "primer_near".to_string()],
            confirmation_report_id: Some("unresolved_gap".to_string()),
            min_3prime_anneal_bp: 6,
            predicted_read_length_bp: 150,
        })
        .expect("suggest primers for unresolved target")
        .sequencing_primer_overlay_report
        .expect("overlay report");

    assert_eq!(report.problem_guidance_count, 1);
    assert_eq!(report.problem_guidance[0].problem_id, "gap_target");
    assert_eq!(
        report.problem_guidance[0]
            .recommended_primer_seq_id
            .as_deref(),
        Some("primer_good")
    );
    assert_eq!(report.problem_guidance[0].candidate_count, 2);
    assert_eq!(
        report.problem_guidance[0].recommended_orientation,
        Some(SequencingPrimerOrientation::ForwardRead)
    );
    assert_eq!(
        report.problem_guidance[0].recommended_three_prime_distance_bp,
        Some(95)
    );
    assert_eq!(
        report.problem_guidance[0].recommended_in_read_direction,
        Some(true)
    );
}

#[test]
fn test_suggest_sequencing_primers_marks_unresolved_target_without_covering_primer() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    state.sequences.insert(
        "read_short".to_string(),
        DNAsequence::from_sequence("AAAA").expect("read"),
    );
    state.sequences.insert(
        "primer_other".to_string(),
        DNAsequence::from_sequence("TTTGGGGGG").expect("primer"),
    );
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_short".to_string()],
            trace_ids: vec![],
            targets: vec![sequencing_confirmation_junction_target(
                "junction_1",
                4,
                12,
                8,
                "Insert junction",
            )],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("uncovered_gap".to_string()),
        })
        .expect("persist unresolved report");

    let report = engine
        .apply(Operation::SuggestSequencingPrimers {
            expected_seq_id: "construct".to_string(),
            primer_seq_ids: vec!["primer_other".to_string()],
            confirmation_report_id: Some("uncovered_gap".to_string()),
            min_3prime_anneal_bp: 6,
            predicted_read_length_bp: 150,
        })
        .expect("suggest primers without coverage")
        .sequencing_primer_overlay_report
        .expect("overlay report");

    assert_eq!(report.problem_guidance_count, 1);
    assert_eq!(report.problem_guidance[0].problem_id, "junction_1");
    assert_eq!(report.problem_guidance[0].recommended_primer_seq_id, None);
    assert_eq!(report.problem_guidance[0].candidate_count, 0);
}

#[test]
fn test_suggest_sequencing_primers_report_only_mode_proposes_fresh_primer_for_unresolved_target() {
    let mut state = ProjectState::default();
    let construct_text = [
        "ACGTTGCAAGTCCTAGTGAC",
        "TTACCGGATGCTACGATCGA",
        "GCTTACAGGATCCGTTAGCA",
        "CGATTCGGAACCTGACTTGA",
        "TGCAGATCCGTACGTTACGA",
        "AGTCGATGGCATTCAGTGCA",
        "CAGTTCGACGGTATGCACTA",
        "TACGAGCTTGACCGTATGGA",
        "GATTCAGCGTACCTGATGCA",
        "CTAGTGACCGTTAGCATGGC",
    ]
    .concat();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence(&construct_text).expect("construct"),
    );
    state.sequences.insert(
        "read_early".to_string(),
        DNAsequence::from_sequence("ACGTTGCAAGTC").expect("read"),
    );
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_early".to_string()],
            trace_ids: vec![],
            targets: vec![SequencingConfirmationTargetSpec {
                target_id: "gap_target".to_string(),
                label: "Gap locus".to_string(),
                kind: SequencingConfirmationTargetKind::Junction,
                start_0based: 118,
                end_0based_exclusive: 122,
                junction_left_end_0based: Some(120),
                expected_bases: None,
                baseline_bases: None,
                required: true,
            }],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("report_only_gap".to_string()),
        })
        .expect("persist unresolved report");

    let report = engine
        .apply(Operation::SuggestSequencingPrimers {
            expected_seq_id: "construct".to_string(),
            primer_seq_ids: vec![],
            confirmation_report_id: Some("report_only_gap".to_string()),
            min_3prime_anneal_bp: 6,
            predicted_read_length_bp: 80,
        })
        .expect("suggest report-only fresh primers")
        .sequencing_primer_overlay_report
        .expect("overlay report");

    assert_eq!(report.suggestion_count, 0);
    assert_eq!(report.problem_guidance_count, 1);
    assert_eq!(report.problem_guidance[0].recommended_primer_seq_id, None);
    assert_eq!(report.proposal_count, 1);
    assert_eq!(report.proposals[0].problem_id, "gap_target");
    assert_eq!(
        report.proposals[0].problem_kind,
        SequencingPrimerProblemKind::Target
    );
    assert_eq!(report.proposals[0].anneal_hits, 1);
    assert!(
        report.proposals[0]
            .reason
            .contains("No existing primer covered")
    );
}

#[test]
fn test_suggest_sequencing_primers_adds_fresh_proposal_when_best_existing_hit_is_outside_window() {
    let mut state = ProjectState::default();
    let construct_text = [
        "ACGTTGCAAGTCCTAGTGAC",
        "TTACCGGATGCTACGATCGA",
        "GCTTACAGGATCCGTTAGCA",
        "CGATTCGGAACCTGACTTGA",
        "TGCAGATCCGTACGTTACGA",
        "AGTCGATGGCATTCAGTGCA",
        "CAGTTCGACGGTATGCACTA",
        "TACGAGCTTGACCGTATGGA",
        "GATTCAGCGTACCTGATGCA",
        "CTAGTGACCGTTAGCATGGC",
    ]
    .concat();
    let primer_near = construct_text[101..110].to_string();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence(&construct_text).expect("construct"),
    );
    state.sequences.insert(
        "read_early".to_string(),
        DNAsequence::from_sequence("ACGTTGCAAGTC").expect("read"),
    );
    state.sequences.insert(
        "primer_near".to_string(),
        DNAsequence::from_sequence(&primer_near).expect("primer"),
    );
    let mut engine = GentleEngine::from_state(state);

    engine
        .apply(Operation::ConfirmConstructReads {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_early".to_string()],
            trace_ids: vec![],
            targets: vec![SequencingConfirmationTargetSpec {
                target_id: "gap_target".to_string(),
                label: "Gap locus".to_string(),
                kind: SequencingConfirmationTargetKind::Junction,
                start_0based: 118,
                end_0based_exclusive: 122,
                junction_left_end_0based: Some(120),
                expected_bases: None,
                baseline_bases: None,
                required: true,
            }],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("window_gap".to_string()),
        })
        .expect("persist unresolved report");

    let report = engine
        .apply(Operation::SuggestSequencingPrimers {
            expected_seq_id: "construct".to_string(),
            primer_seq_ids: vec!["primer_near".to_string()],
            confirmation_report_id: Some("window_gap".to_string()),
            min_3prime_anneal_bp: 6,
            predicted_read_length_bp: 80,
        })
        .expect("suggest primers with fresh fallback")
        .sequencing_primer_overlay_report
        .expect("overlay report");

    assert_eq!(
        report.problem_guidance[0]
            .recommended_primer_seq_id
            .as_deref(),
        Some("primer_near")
    );
    assert_eq!(report.proposal_count, 1);
    assert_eq!(report.proposals[0].problem_id, "gap_target");
    assert!(
        report.proposals[0]
            .reason
            .contains("Existing primer 'primer_near'")
    );
}

#[test]
fn test_parse_fasta_records_with_offsets_supports_gzip_input() {
    let td = tempdir().expect("tempdir");
    let fasta_gz = td.path().join("reads.fa.gz");
    write_gzip(&fasta_gz, ">read_1 comment\nAUGT\n>read_2\nCC\nGG\n");

    let records =
        GentleEngine::parse_fasta_records_with_offsets(fasta_gz.to_str().expect("utf-8 path"))
            .expect("parse gzip FASTA");
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].header_id, "read_1");
    assert_eq!(records[0].sequence, b"ATGT".to_vec());
    assert_eq!(records[1].header_id, "read_2");
    assert_eq!(records[1].sequence, b"CCGG".to_vec());
    assert!(records[1].source_byte_offset > records[0].source_byte_offset);
}

#[test]
fn test_parse_fasta_records_with_offsets_supports_concatenated_gzip_input() {
    let td = tempdir().expect("tempdir");
    let fasta_gz = td.path().join("reads_concat.fa.gz");
    write_multi_member_gzip(
        &fasta_gz,
        &[
            ">read_1 comment\nAUGT\n",
            ">read_2\nCC\nGG\n",
            ">read_3\nUUAA\n",
        ],
    );

    let records =
        GentleEngine::parse_fasta_records_with_offsets(fasta_gz.to_str().expect("utf-8 path"))
            .expect("parse concatenated gzip FASTA");
    assert_eq!(records.len(), 3);
    assert_eq!(records[0].header_id, "read_1");
    assert_eq!(records[0].sequence, b"ATGT".to_vec());
    assert_eq!(records[1].header_id, "read_2");
    assert_eq!(records[1].sequence, b"CCGG".to_vec());
    assert_eq!(records[2].header_id, "read_3");
    assert_eq!(records[2].sequence, b"TTAA".to_vec());
    assert!(records[2].source_byte_offset > records[1].source_byte_offset);
}

#[test]
fn test_rna_read_seed_filter_default_kmer_len_is_10() {
    assert_eq!(RnaReadSeedFilterConfig::default().kmer_len, 10);
}

#[test]
fn test_summarize_rna_read_gene_support_filters_requested_gene_in_multi_gene_sparse_report() {
    let mut engine = build_rna_read_gene_support_test_engine();
    let result = engine
        .apply(Operation::SummarizeRnaReadGeneSupport {
            report_id: "rna_reads_gene_support".to_string(),
            gene_ids: vec!["GENE1".to_string()],
            selected_record_indices: vec![],
            complete_rule: RnaReadGeneSupportCompleteRule::Near,
            path: None,
        })
        .expect("summarize RNA-read gene support");
    let summary = result
        .rna_read_gene_support_summary
        .expect("gene-support summary payload");

    assert_eq!(
        summary.schema,
        "gentle.rna_read_gene_support_summary.v1".to_string()
    );
    assert_eq!(summary.report_id, "rna_reads_gene_support");
    assert_eq!(summary.seq_id, "seq_a");
    assert_eq!(summary.requested_gene_ids, vec!["GENE1".to_string()]);
    assert_eq!(summary.matched_gene_ids, vec!["GENE1".to_string()]);
    assert!(summary.missing_gene_ids.is_empty());
    assert_eq!(summary.aligned_base_count, 4);
    assert_eq!(summary.accepted_target_count, 3);
    assert_eq!(summary.all_target.read_count, 3);
    assert!(
        summary
            .all_target
            .exon_support
            .iter()
            .all(|row| row.gene_id == "GENE1")
    );
    assert!(
        result
            .messages
            .iter()
            .any(|message| message.contains("accepted_target_reads=3"))
    );
}

#[test]
fn test_summarize_rna_read_gene_support_splits_fragments_from_complete_reads() {
    let engine = build_rna_read_gene_support_test_engine();
    let summary = engine
        .summarize_rna_read_gene_support(
            "rna_reads_gene_support",
            &[String::from("GENE1")],
            &[],
            RnaReadGeneSupportCompleteRule::Near,
        )
        .expect("summarize RNA-read gene support");

    assert_eq!(summary.complete_rule, RnaReadGeneSupportCompleteRule::Near);
    assert_eq!(summary.accepted_target_count, 3);
    assert_eq!(summary.fragment_count, 1);
    assert_eq!(summary.complete_count, 2);
    assert_eq!(summary.complete_strict_count, 1);
    assert_eq!(summary.complete_exact_count, 2);
    assert_eq!(summary.fragments.read_count, 1);
    assert_eq!(summary.complete.read_count, 2);
}

#[test]
fn test_summarize_rna_read_gene_support_counts_skipped_exon_pairs_without_direct_transition() {
    let engine = build_rna_read_gene_support_test_engine();
    let summary = engine
        .summarize_rna_read_gene_support(
            "rna_reads_gene_support",
            &[String::from("GENE1")],
            &[2],
            RnaReadGeneSupportCompleteRule::Near,
        )
        .expect("summarize RNA-read gene support");

    let pair = find_gene_support_pair_row(&summary.all_target.exon_pair_support, "GENE1", 1, 3)
        .expect("skipped exon pair support");
    assert_eq!(pair.support_read_count, 1);
    assert!((pair.support_fraction - 1.0).abs() < f64::EPSILON);
    assert!(
        find_gene_support_pair_row(&summary.all_target.direct_transition_support, "GENE1", 1, 3)
            .is_none()
    );
}

#[test]
fn test_summarize_rna_read_gene_support_counts_adjacent_pairs_as_direct_transitions() {
    let engine = build_rna_read_gene_support_test_engine();
    let summary = engine
        .summarize_rna_read_gene_support(
            "rna_reads_gene_support",
            &[String::from("GENE1")],
            &[0],
            RnaReadGeneSupportCompleteRule::Near,
        )
        .expect("summarize RNA-read gene support");

    let exon1 = find_gene_support_exon_row(&summary.all_target.exon_support, "GENE1", 1)
        .expect("exon 1 support");
    let exon2 = find_gene_support_exon_row(&summary.all_target.exon_support, "GENE1", 2)
        .expect("exon 2 support");
    let adjacent_pair =
        find_gene_support_pair_row(&summary.all_target.exon_pair_support, "GENE1", 1, 2)
            .expect("adjacent exon pair support");
    let direct_transition =
        find_gene_support_pair_row(&summary.all_target.direct_transition_support, "GENE1", 1, 2)
            .expect("adjacent direct transition support");

    assert_eq!(exon1.support_read_count, 1);
    assert_eq!(exon2.support_read_count, 1);
    assert_eq!(adjacent_pair.support_read_count, 1);
    assert!((adjacent_pair.support_fraction - 1.0).abs() < f64::EPSILON);
    assert_eq!(direct_transition.support_read_count, 1);
    assert!((direct_transition.support_fraction - 1.0).abs() < f64::EPSILON);
}

#[test]
fn test_summarize_rna_read_gene_support_honors_selected_record_indices_before_gene_filtering() {
    let engine = build_rna_read_gene_support_test_engine();
    let summary = engine
        .summarize_rna_read_gene_support(
            "rna_reads_gene_support",
            &[String::from("GENE1")],
            &[3, 2, 2],
            RnaReadGeneSupportCompleteRule::Near,
        )
        .expect("summarize RNA-read gene support");

    assert_eq!(summary.selected_record_indices, vec![2, 3]);
    assert_eq!(summary.aligned_base_count, 2);
    assert_eq!(summary.accepted_target_count, 1);
    assert_eq!(summary.fragment_count, 0);
    assert_eq!(summary.complete_count, 1);
    assert_eq!(summary.all_target.read_count, 1);
    assert!(
        find_gene_support_pair_row(&summary.all_target.exon_pair_support, "GENE1", 1, 3).is_some()
    );
    assert!(
        summary
            .all_target
            .exon_support
            .iter()
            .all(|row| row.gene_id == "GENE1")
    );
}

#[test]
fn test_rna_read_top_hit_preview_propagates_alignment_summary() {
    let best_mapping = RnaReadMappingHit {
        alignment_mode: RnaReadAlignmentMode::Semiglobal,
        transcript_feature_id: 7,
        transcript_id: "NM_TEST_1".to_string(),
        transcript_label: "TP73-201".to_string(),
        strand: "+".to_string(),
        query_start_0based: 0,
        query_end_0based_exclusive: 120,
        query_reverse_complemented: false,
        target_start_1based: 101,
        target_end_1based: 220,
        target_start_offset_0based: 0,
        target_end_offset_0based_exclusive: 120,
        matches: 114,
        mismatches: 6,
        score: 342,
        identity_fraction: 0.95,
        query_coverage_fraction: 0.92,
    };
    let secondary_mapping = RnaReadMappingHit {
        alignment_mode: RnaReadAlignmentMode::Local,
        transcript_feature_id: 8,
        transcript_id: "NM_TEST_2".to_string(),
        transcript_label: "TP73-202".to_string(),
        strand: "+".to_string(),
        query_start_0based: 4,
        query_end_0based_exclusive: 104,
        query_reverse_complemented: false,
        target_start_1based: 140,
        target_end_1based: 239,
        target_start_offset_0based: 40,
        target_end_offset_0based_exclusive: 140,
        matches: 90,
        mismatches: 10,
        score: 210,
        identity_fraction: 0.90,
        query_coverage_fraction: 0.83,
    };
    let hit = RnaReadInterpretationHit {
        record_index: 41,
        header_id: "read_42".to_string(),
        sequence: "ACGT".repeat(40),
        read_length_bp: 160,
        best_mapping: Some(best_mapping),
        secondary_mappings: vec![secondary_mapping],
        ..Default::default()
    };

    let preview = GentleEngine::make_rna_read_top_hit_preview(&hit);
    assert!(preview.aligned);
    assert_eq!(preview.best_alignment_mode, "semiglobal");
    assert_eq!(preview.best_alignment_transcript_id, "NM_TEST_1");
    assert_eq!(preview.best_alignment_transcript_label, "TP73-201");
    assert_eq!(preview.best_alignment_strand, "+");
    assert_eq!(preview.best_alignment_target_start_1based, 101);
    assert_eq!(preview.best_alignment_target_end_1based, 220);
    assert!((preview.best_alignment_identity_fraction - 0.95).abs() < f64::EPSILON);
    assert!((preview.best_alignment_query_coverage_fraction - 0.92).abs() < f64::EPSILON);
    assert_eq!(preview.best_alignment_score, 342);
    assert_eq!(preview.secondary_mapping_count, 1);

    let preview_without_mapping =
        GentleEngine::make_rna_read_top_hit_preview(&RnaReadInterpretationHit::default());
    assert!(!preview_without_mapping.aligned);
    assert!(preview_without_mapping.best_alignment_mode.is_empty());
    assert!(
        preview_without_mapping
            .best_alignment_transcript_id
            .is_empty()
    );
    assert_eq!(
        preview_without_mapping.best_alignment_target_start_1based,
        0
    );
    assert_eq!(preview_without_mapping.best_alignment_target_end_1based, 0);
    assert_eq!(preview_without_mapping.secondary_mapping_count, 0);
}

#[test]
fn test_interpret_rna_reads_accepts_gzip_fasta_input() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };

    let td = tempdir().expect("tempdir");
    let input_gz = td.path().join("reads.fa.gz");
    write_gzip(&input_gz, &format!(">read_1\n{read_sequence}\n"));
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.min_seed_hit_fraction = 0.0;
    seed_filter.min_weighted_seed_hit_fraction = 0.0;
    seed_filter.min_unique_matched_kmers = 0;
    seed_filter.min_chain_consistency_fraction = 0.0;
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;

    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_gz.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_gz".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret gzip FASTA");

    let report = engine
        .get_rna_read_report("rna_reads_gz")
        .expect("stored RNA-read report");
    assert_eq!(report.read_count_total, 1);
    assert_eq!(report.hits[0].header_id, "read_1");
    assert!(report.hits[0].passed_seed_filter);
}

#[test]
fn test_align_rna_read_report_selected_record_indices_overrides_selection() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_align_selected.fa");
    fs::write(
        &input_path,
        format!(">read_0\n{read_sequence}\n>read_1\n{read_sequence}\n"),
    )
    .expect("write reads");
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.min_seed_hit_fraction = 0.0;
    seed_filter.min_weighted_seed_hit_fraction = 0.0;
    seed_filter.min_unique_matched_kmers = 0;
    seed_filter.min_chain_consistency_fraction = 0.0;
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;

    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_align_selected".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret FASTA");
    let before_align = engine
        .get_rna_read_report("rna_reads_align_selected")
        .expect("stored RNA-read report");
    assert_eq!(before_align.read_count_total, 2);
    assert_eq!(before_align.read_count_aligned, 0);
    assert!(
        before_align
            .hits
            .iter()
            .all(|hit| hit.best_mapping.is_none())
    );

    engine
        .apply(Operation::AlignRnaReadReport {
            report_id: "rna_reads_align_selected".to_string(),
            selection: RnaReadHitSelection::All,
            align_config_override: Some(RnaReadAlignConfig {
                band_width_bp: 24,
                min_identity_fraction: 0.60,
                max_secondary_mappings: 0,
            }),
            selected_record_indices: vec![0],
        })
        .expect("align selected record only");
    let aligned_report = engine
        .get_rna_read_report("rna_reads_align_selected")
        .expect("aligned report");
    assert_eq!(aligned_report.read_count_total, 2);
    assert_eq!(aligned_report.read_count_aligned, 1);
    let aligned_record_indices = aligned_report
        .hits
        .iter()
        .filter(|hit| hit.best_mapping.is_some())
        .map(|hit| hit.record_index)
        .collect::<Vec<_>>();
    assert_eq!(aligned_record_indices, vec![0]);
    assert!(
        aligned_report
            .warnings
            .iter()
            .any(|warning| warning.contains("explicit_record_indices=1"))
    );
}

#[test]
fn test_inspect_rna_read_alignment_detail_reconstructs_template_alignment() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_align_detail.fa");
    fs::write(&input_path, format!(">read_0\n{read_sequence}\n")).expect("write reads");
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.min_seed_hit_fraction = 0.0;
    seed_filter.min_weighted_seed_hit_fraction = 0.0;
    seed_filter.min_unique_matched_kmers = 0;
    seed_filter.min_chain_consistency_fraction = 0.0;
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;

    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_align_detail".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret FASTA");

    engine
        .apply(Operation::AlignRnaReadReport {
            report_id: "rna_reads_align_detail".to_string(),
            selection: RnaReadHitSelection::SeedPassed,
            align_config_override: None,
            selected_record_indices: vec![],
        })
        .expect("align retained row");

    let detail = engine
        .inspect_rna_read_alignment_detail("rna_reads_align_detail", 0)
        .expect("alignment detail");
    assert_eq!(detail.report_id, "rna_reads_align_detail");
    assert_eq!(detail.record_index, 0);
    assert_eq!(detail.header_id, "read_0");
    assert!(detail.aligned_columns > 0);
    assert_eq!(detail.matches, detail.aligned_columns);
    assert_eq!(detail.mismatches, 0);
    assert_eq!(detail.insertions, 0);
    assert_eq!(detail.deletions, 0);
    assert!((detail.identity_fraction - 1.0).abs() < 1e-9);
    assert!((detail.query_coverage_fraction - 1.0).abs() < 1e-9);
    assert_eq!(detail.target_length_bp, read_sequence.len());
    assert!((detail.target_coverage_fraction - 1.0).abs() < 1e-9);
    assert_eq!(detail.aligned_query, detail.aligned_target);
    assert!(detail.aligned_relation.chars().all(|ch| ch == '|'));
}

#[test]
fn test_transition_support_counts_only_seed_passed_reads() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let (read_sequence, exon_count) = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let transcript = splicing
            .transcripts
            .iter()
            .find(|tx| tx.exons.len() > 1)
            .expect("expected multi-exon transcript for transition test");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template = GentleEngine::make_transcript_template(dna, transcript, kmer_len);
        (
            String::from_utf8(template.sequence).expect("template sequence utf-8"),
            transcript.exons.len(),
        )
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_transition_gate.fa");
    fs::write(&input_path, format!(">read_1\n{read_sequence}\n")).expect("write reads");

    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.min_confirmed_exon_transitions = exon_count; // guaranteed above exon_count-1 transitions

    let mut final_progress: Option<RnaReadInterpretProgress> = None;
    let report = engine
        .compute_rna_read_report_with_progress_and_cancel(
            "seq_a",
            feature_id,
            RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path.to_str().expect("path"),
            RnaReadInputFormat::Fasta,
            SplicingScopePreset::AllOverlappingBothStrands,
            RnaReadOriginMode::SingleGene,
            &[],
            false,
            &seed_filter,
            &RnaReadAlignConfig::default(),
            Some("rna_reads_transition_gate"),
            &mut |progress| {
                if let OperationProgress::RnaReadInterpret(p) = progress {
                    if p.done {
                        final_progress = Some(p.clone());
                    }
                }
                true
            },
            &mut || true,
        )
        .expect("compute report");

    assert_eq!(report.read_count_total, 1);
    assert_eq!(report.read_count_seed_passed, 0);
    assert!(report.hits.iter().all(|hit| !hit.passed_seed_filter));
    assert!(
        report
            .hits
            .iter()
            .any(|hit| hit.exon_transitions_total > 0 && hit.exon_transitions_confirmed > 0),
        "expected transition-bearing read hit for regression guard"
    );
    assert!(
        report
            .transition_support_rows
            .iter()
            .all(|row| row.support_read_count == 0),
        "transition rows must not count reads that fail seed filter"
    );

    let progress = final_progress.expect("final progress event");
    assert_eq!(progress.seed_passed, 0);
    assert_eq!(progress.reads_with_transition_support, 0);
    assert_eq!(progress.transition_confirmations, 0);
    assert!(
        progress
            .transition_support_rows
            .iter()
            .all(|row| row.support_read_count == 0),
        "progress transition rows must not count reads that fail seed filter"
    );
}

fn seed_failed_but_alignable_rna_read_report() -> (GentleEngine, String) {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let report_id = "rna_reads_seed_failed_alignment".to_string();
    let strict_seed_filter = RnaReadSeedFilterConfig {
        min_confirmed_exon_transitions: 99,
        ..RnaReadSeedFilterConfig::default()
    };
    engine
        .upsert_rna_read_report(RnaReadInterpretationReport {
            schema: RNA_READ_REPORT_SCHEMA.to_string(),
            report_id: report_id.clone(),
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            seed_filter: strict_seed_filter,
            align_config: RnaReadAlignConfig::default(),
            hits: vec![RnaReadInterpretationHit {
                record_index: 0,
                header_id: "read_0".to_string(),
                sequence: read_sequence.clone(),
                read_length_bp: read_sequence.len(),
                tested_kmers: read_sequence
                    .len()
                    .saturating_sub(kmer_len)
                    .saturating_add(1),
                matched_kmers: read_sequence
                    .len()
                    .saturating_sub(kmer_len)
                    .saturating_add(1),
                seed_hit_fraction: 1.0,
                weighted_seed_hit_fraction: 1.0,
                weighted_matched_kmers: read_sequence.len() as f64,
                passed_seed_filter: false,
                ..RnaReadInterpretationHit::default()
            }],
            ..RnaReadInterpretationReport::default()
        })
        .expect("upsert seed-failed report");
    (engine, report_id)
}

#[test]
fn test_align_rna_read_report_aligns_selected_retained_rows_even_when_seed_gate_stays_false() {
    let (mut engine, report_id) = seed_failed_but_alignable_rna_read_report();

    engine
        .apply(Operation::AlignRnaReadReport {
            report_id: report_id.clone(),
            selection: RnaReadHitSelection::All,
            align_config_override: Some(RnaReadAlignConfig {
                band_width_bp: 24,
                min_identity_fraction: 0.60,
                max_secondary_mappings: 0,
            }),
            selected_record_indices: vec![],
        })
        .expect("align all retained rows");

    let report = engine
        .get_rna_read_report(&report_id)
        .expect("stored RNA-read report");
    assert_eq!(report.read_count_aligned, 1);
    assert_eq!(report.hits.len(), 1);
    assert!(!report.hits[0].passed_seed_filter);
    assert!(report.hits[0].best_mapping.is_some());
}

#[test]
fn test_align_rna_read_report_seed_passed_falls_back_to_raw_min_hit_rows() {
    let (mut engine, report_id) = seed_failed_but_alignable_rna_read_report();

    engine
        .apply(Operation::AlignRnaReadReport {
            report_id: report_id.clone(),
            selection: RnaReadHitSelection::SeedPassed,
            align_config_override: Some(RnaReadAlignConfig {
                band_width_bp: 24,
                min_identity_fraction: 0.60,
                max_secondary_mappings: 0,
            }),
            selected_record_indices: vec![],
        })
        .expect("align with seed_passed fallback");

    let report = engine
        .get_rna_read_report(&report_id)
        .expect("stored RNA-read report");
    assert_eq!(report.read_count_aligned, 1);
    assert!(report.hits[0].best_mapping.is_some());
    assert!(report.warnings.iter().any(|warning| {
        warning.contains("selection 'seed_passed' matched no retained hits")
            && warning.contains("fell back to 1 retained row(s) at or above raw min_hit")
    }));
}

#[test]
fn test_align_rna_read_report_progress_uses_stride_one() {
    let (mut engine, report_id) = seed_failed_but_alignable_rna_read_report();
    let mut progress_events = Vec::<RnaReadInterpretProgress>::new();

    engine
        .apply_with_progress(
            Operation::AlignRnaReadReport {
                report_id,
                selection: RnaReadHitSelection::All,
                align_config_override: Some(RnaReadAlignConfig {
                    band_width_bp: 24,
                    min_identity_fraction: 0.60,
                    max_secondary_mappings: 0,
                }),
                selected_record_indices: vec![],
            },
            |progress| {
                if let OperationProgress::RnaReadInterpret(p) = progress {
                    progress_events.push(p);
                }
                true
            },
        )
        .expect("align all retained rows with progress");

    assert!(!progress_events.is_empty());
    assert!(progress_events.iter().all(|p| p.update_every_reads == 1));
    let last = progress_events.last().expect("final progress");
    assert!(last.done);
    assert_eq!(last.reads_total, 1);
    assert_eq!(last.reads_processed, 1);
}

#[test]
fn test_seed_passed_only_report_mode_keeps_raw_min_hit_candidates() {
    let kept = GentleEngine::apply_rna_read_report_mode_to_hits(
        RnaReadReportMode::SeedPassedOnly,
        0.30,
        vec![
            RnaReadInterpretationHit {
                record_index: 1,
                header_id: "raw_candidate".to_string(),
                seed_hit_fraction: 0.31,
                passed_seed_filter: false,
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 2,
                header_id: "below_min_hit".to_string(),
                seed_hit_fraction: 0.12,
                passed_seed_filter: false,
                ..RnaReadInterpretationHit::default()
            },
        ],
    );

    assert_eq!(kept.len(), 1);
    assert_eq!(kept[0].header_id, "raw_candidate");
    assert!(!kept[0].passed_seed_filter);
}

#[test]
fn test_interpret_rna_reads_poly_t_cdna_flip_sets_rc_flag_and_sequence() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_poly_t.fa");
    fs::write(&input_path, ">read_1\nTTTACGTACGT\n").expect("write reads");

    let mut cdna_seed_filter = RnaReadSeedFilterConfig::default();
    cdna_seed_filter.poly_t_prefix_min_bp = 3;
    cdna_seed_filter.kmer_len = 3;
    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: Default::default(),
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter: cdna_seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_poly_t_cdna".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret cDNA-style reads");
    let cdna_report = engine
        .get_rna_read_report("rna_reads_poly_t_cdna")
        .expect("cdna report");
    assert_eq!(cdna_report.hits.len(), 1);
    assert!(cdna_report.hits[0].reverse_complement_applied);
    assert_eq!(cdna_report.hits[0].sequence, "ACGTACGTAAA");

    let mut direct_seed_filter = RnaReadSeedFilterConfig::default();
    direct_seed_filter.poly_t_prefix_min_bp = 3;
    direct_seed_filter.kmer_len = 3;
    direct_seed_filter.cdna_poly_t_flip_enabled = false;
    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: Default::default(),
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter: direct_seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_poly_t_direct".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret direct-RNA-style reads");
    let direct_report = engine
        .get_rna_read_report("rna_reads_poly_t_direct")
        .expect("direct report");
    assert_eq!(direct_report.hits.len(), 1);
    assert!(!direct_report.hits[0].reverse_complement_applied);
    assert_eq!(direct_report.hits[0].sequence, "TTTACGTACGT");
}

#[test]
fn test_align_read_to_template_supports_reverse_complement_query_orientation() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let engine = GentleEngine::from_state(state);
    let dna = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present");
    let feature_id = dna
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let splicing = engine
        .build_splicing_expert_view(
            "seq_a",
            feature_id,
            SplicingScopePreset::AllOverlappingBothStrands,
        )
        .expect("splicing view");
    let lane = splicing.transcripts.first().expect("first transcript lane");
    let template = GentleEngine::make_transcript_template(dna, lane, 3);
    let read = GentleEngine::reverse_complement_bytes(&template.sequence);

    let mapping = GentleEngine::align_read_to_template(
        &read,
        &template,
        &RnaReadAlignConfig {
            min_identity_fraction: 0.60,
            ..RnaReadAlignConfig::default()
        },
        3,
    )
    .expect("reverse-complement alignment");

    assert!(mapping.query_reverse_complemented);
    assert_eq!(mapping.query_start_0based, 0);
    assert_eq!(mapping.query_end_0based_exclusive, read.len());
    assert_eq!(mapping.transcript_feature_id, lane.transcript_feature_id);
    assert_eq!(mapping.alignment_mode, RnaReadAlignmentMode::Semiglobal);
    assert!((mapping.identity_fraction - 1.0).abs() < f64::EPSILON);
    assert!((mapping.query_coverage_fraction - 1.0).abs() < f64::EPSILON);
}

#[test]
fn test_poly_t_cdna_flip_accepts_disrupted_t_rich_heads() {
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.poly_t_prefix_min_bp = 18;
    seed_filter.cdna_poly_t_flip_enabled = true;

    let disrupted_examples = [
        "TTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTAAGGTGGCAGGCTTTTAATTTCC",
        "TTTTTTTTTTTTTTTTTATTTTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGGGATTCTGCCAAAAGGA",
        "TTTTTTTTTTTTTGTTTTTTTGTTTTTTGAGTGTGAAAAATAAACTATTTTATTTCA",
        "TTTTTTTTTTTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACCTCATTT",
        "TTGTCATATACTTATTTCTTTTTTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTAAATTGTATTTATAGATATCTGGACTCAGAACT",
    ];
    for sequence in disrupted_examples {
        let (_normalized, rc_applied) = GentleEngine::normalize_rna_read_sequence_for_scoring(
            sequence.as_bytes(),
            &seed_filter,
        );
        assert!(
            rc_applied,
            "expected disrupted T-rich cDNA head to be normalized: {}",
            sequence
        );
    }

    let (negative_normalized, negative_rc_applied) =
        GentleEngine::normalize_rna_read_sequence_for_scoring(
            b"ACGTCAGGACTGACCGTACCGATCGATCGATCGAC",
            &seed_filter,
        );
    assert!(
        !negative_rc_applied,
        "non poly-T heads must stay in original orientation"
    );
    assert_eq!(
        negative_normalized,
        b"ACGTCAGGACTGACCGTACCGATCGATCGATCGAC".to_vec()
    );
}

#[test]
fn test_interpret_rna_reads_multi_gene_sparse_persists_and_warns_for_missing_targets() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_sparse_mode.fa");
    fs::write(&input_path, format!(">read_1\n{read_sequence}\n")).expect("write reads");

    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: RnaReadOriginMode::MultiGeneSparse,
            target_gene_ids: vec!["TP73".to_string(), "TP53".to_string(), "TP73".to_string()],
            roi_seed_capture_enabled: true,
            seed_filter: RnaReadSeedFilterConfig::default(),
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_sparse_mode".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret reads with sparse-origin scaffolding");
    let report = engine
        .get_rna_read_report("rna_reads_sparse_mode")
        .expect("stored report");
    assert_eq!(report.origin_mode, RnaReadOriginMode::MultiGeneSparse);
    assert_eq!(
        report.target_gene_ids,
        vec!["TP53".to_string(), "TP73".to_string()]
    );
    assert!(report.roi_seed_capture_enabled);
    assert!(
        report
            .warnings
            .iter()
            .any(|warning| warning.contains("origin_mode=multi_gene_sparse active"))
    );
    assert!(
        report
            .warnings
            .iter()
            .any(|warning| warning.contains("target genes not found in local annotation"))
    );
    assert!(
        report
            .warnings
            .iter()
            .any(|warning| { warning.contains("roi_seed_capture_enabled=true requested") })
    );
    let summaries = engine.list_rna_read_reports(Some("seq_a"));
    assert_eq!(summaries.len(), 1);
    let summary = &summaries[0];
    assert_eq!(summary.report_id, "rna_reads_sparse_mode");
    assert_eq!(summary.origin_mode, RnaReadOriginMode::MultiGeneSparse);
    assert_eq!(summary.target_gene_count, 2);
    assert!(summary.roi_seed_capture_enabled);
}

#[test]
fn test_interpret_rna_reads_multi_gene_sparse_adds_target_gene_templates() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_multi_gene_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let features = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features();
    let seed_feature_id = features
        .iter()
        .enumerate()
        .find(|(idx, feature)| {
            feature.kind.to_string().eq_ignore_ascii_case("mRNA")
                && GentleEngine::splicing_group_label(feature, *idx).eq_ignore_ascii_case("GENE1")
        })
        .map(|(idx, _)| idx)
        .expect("GENE1 seed feature");
    let gene2_feature_id = features
        .iter()
        .enumerate()
        .find(|(idx, feature)| {
            feature.kind.to_string().eq_ignore_ascii_case("mRNA")
                && GentleEngine::splicing_group_label(feature, *idx).eq_ignore_ascii_case("GENE2")
        })
        .map(|(idx, _)| idx)
        .expect("GENE2 feature");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                gene2_feature_id,
                SplicingScopePreset::TargetGroupTargetStrand,
            )
            .expect("GENE2 splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_gene2.fa");
    fs::write(&input_path, format!(">read_gene2\n{read_sequence}\n")).expect("write reads");

    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::TargetGroupTargetStrand,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec!["GENE2".to_string()],
            roi_seed_capture_enabled: false,
            seed_filter: RnaReadSeedFilterConfig::default(),
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_single_gene2".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("single-gene baseline run");
    let baseline = engine
        .get_rna_read_report("rna_reads_single_gene2")
        .expect("baseline report");
    assert_eq!(baseline.read_count_total, 1);
    assert_eq!(baseline.read_count_seed_passed, 0);

    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::TargetGroupTargetStrand,
            origin_mode: RnaReadOriginMode::MultiGeneSparse,
            target_gene_ids: vec!["GENE2".to_string()],
            roi_seed_capture_enabled: false,
            seed_filter: RnaReadSeedFilterConfig::default(),
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_multi_gene2".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("multi-gene sparse run");
    let multi = engine
        .get_rna_read_report("rna_reads_multi_gene2")
        .expect("multi report");
    assert_eq!(multi.read_count_total, 1);
    assert_eq!(multi.read_count_seed_passed, 1);
    assert!(
        multi
            .warnings
            .iter()
            .any(|warning| warning.contains("added 1 transcript lane(s)"))
    );
}

#[test]
fn test_list_and_show_rna_read_reports_messages_include_origin_provenance() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_sparse_message.fa");
    fs::write(&input_path, format!(">read_1\n{read_sequence}\n")).expect("write reads");

    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: RnaReadOriginMode::MultiGeneSparse,
            target_gene_ids: vec!["TP73".to_string(), "TP53".to_string()],
            roi_seed_capture_enabled: true,
            seed_filter: RnaReadSeedFilterConfig::default(),
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_sparse_message".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret reads");

    let listed = engine
        .apply(Operation::ListRnaReadReports {
            seq_id: Some("seq_a".to_string()),
        })
        .expect("list rna-read reports");
    let listed_text = listed.messages.join("\n");
    assert!(listed_text.contains("origin=multi_gene_sparse"));
    assert!(listed_text.contains("targets=2"));
    assert!(listed_text.contains("roi_capture=true"));

    let shown = engine
        .apply(Operation::ShowRnaReadReport {
            report_id: "rna_reads_sparse_message".to_string(),
        })
        .expect("show rna-read report");
    let shown_text = shown.messages.join("\n");
    assert!(shown_text.contains("origin=multi_gene_sparse"));
    assert!(shown_text.contains("targets=2"));
    assert!(shown_text.contains("roi_capture=true"));
    assert!(shown_text.contains("target_genes=TP53,TP73"));
}

#[test]
fn test_interpret_rna_reads_accepts_reverse_strand_tp73_ncrna_seed() {
    let mut engine = GentleEngine::default();
    engine
        .apply(Operation::LoadFile {
            path: "test_files/tp73.ncbi.gb".to_string(),
            as_id: Some("tp73".to_string()),
        })
        .expect("load tp73 fixture");

    let ncrna_feature_id = engine
        .state()
        .sequences
        .get("tp73")
        .expect("tp73 sequence present")
        .features()
        .iter()
        .position(|feature| {
            feature.kind.to_string().eq_ignore_ascii_case("ncRNA")
                && feature
                    .qualifier_values("gene")
                    .any(|value| value.eq_ignore_ascii_case("TP73-AS3"))
        })
        .expect("TP73-AS3 ncRNA feature");

    let splicing = engine
        .build_splicing_expert_view(
            "tp73",
            ncrna_feature_id,
            SplicingScopePreset::TargetGroupTargetStrand,
        )
        .expect("build ncRNA splicing view");
    assert_eq!(splicing.strand, "-");
    assert!(!splicing.transcripts.is_empty());

    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let dna = engine
            .state()
            .sequences
            .get("tp73")
            .expect("tp73 sequence present");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence")
    };

    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("tp73_as3_reads.fa");
    fs::write(&input_path, format!(">read_1\n{read_sequence}\n")).expect("write reads");

    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "tp73".to_string(),
            seed_feature_id: ncrna_feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::TargetGroupTargetStrand,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter: RnaReadSeedFilterConfig::default(),
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("tp73_as3_reverse".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret ncRNA reads");

    let report = engine
        .get_rna_read_report("tp73_as3_reverse")
        .expect("stored report");
    assert_eq!(report.read_count_total, 1);
    assert_eq!(report.read_count_seed_passed, 1);
}

#[test]
fn test_interpret_rna_reads_accepts_misc_rna_seed() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_misc".to_string(), splicing_misc_rna_sequence());
    let mut engine = GentleEngine::from_state(state);

    let splicing = engine
        .build_splicing_expert_view("seq_misc", 0, SplicingScopePreset::TargetGroupTargetStrand)
        .expect("build misc_RNA splicing view");
    assert_eq!(splicing.transcript_count, 2);

    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let dna = engine
            .state()
            .sequences
            .get("seq_misc")
            .expect("misc_RNA sequence present");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence")
    };

    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("misc_rna_reads.fa");
    fs::write(&input_path, format!(">read_1\n{read_sequence}\n")).expect("write reads");
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.min_seed_hit_fraction = 0.0;
    seed_filter.min_weighted_seed_hit_fraction = 0.0;
    seed_filter.min_unique_matched_kmers = 0;
    seed_filter.min_chain_consistency_fraction = 0.0;
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;

    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_misc".to_string(),
            seed_feature_id: 0,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::TargetGroupTargetStrand,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("seq_misc_reads".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret misc_RNA reads");

    let report = engine
        .get_rna_read_report("seq_misc_reads")
        .expect("stored report");
    assert_eq!(report.read_count_total, 1);
    assert_eq!(report.read_count_seed_passed, 1);
    assert!(report.hits.iter().any(|hit| hit.passed_seed_filter));
}

#[test]
fn test_rna_read_progress_emit_policy_supports_read_and_timer_triggers() {
    assert!(GentleEngine::should_emit_rna_read_progress(
        1,
        Duration::from_millis(50),
        RNA_READ_PROGRESS_UPDATE_EVERY_READS,
    ));
    assert!(!GentleEngine::should_emit_rna_read_progress(
        RNA_READ_PROGRESS_UPDATE_EVERY_READS.saturating_sub(1),
        Duration::from_secs(1),
        RNA_READ_PROGRESS_UPDATE_EVERY_READS,
    ));
    assert!(GentleEngine::should_emit_rna_read_progress(
        RNA_READ_PROGRESS_UPDATE_EVERY_READS,
        Duration::from_millis(10),
        RNA_READ_PROGRESS_UPDATE_EVERY_READS,
    ));
    assert!(GentleEngine::should_emit_rna_read_progress(
        RNA_READ_PROGRESS_UPDATE_EVERY_READS.saturating_sub(1),
        Duration::from_secs(2),
        RNA_READ_PROGRESS_UPDATE_EVERY_READS,
    ));
}

#[test]
#[cfg(debug_assertions)]
fn test_rna_read_progress_update_stride_default_debug() {
    assert_eq!(RNA_READ_PROGRESS_UPDATE_EVERY_READS, 1_000);
}

#[test]
fn test_interpret_rna_reads_progress_reports_histogram_updates() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_many.fa");
    let mut fasta = String::new();
    for idx in 0..1001usize {
        fasta.push_str(&format!(">read_{}\n{}\n", idx + 1, read_sequence));
    }
    fs::write(&input_path, fasta).expect("write reads");

    let mut progress_events = Vec::<RnaReadInterpretProgress>::new();
    engine
        .apply_with_progress(
            Operation::InterpretRnaReads {
                seq_id: "seq_a".to_string(),
                seed_feature_id: feature_id,
                profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
                input_path: input_path.display().to_string(),
                input_format: RnaReadInputFormat::Fasta,
                scope: SplicingScopePreset::AllOverlappingBothStrands,
                origin_mode: Default::default(),
                target_gene_ids: vec![],
                roi_seed_capture_enabled: false,
                seed_filter: RnaReadSeedFilterConfig::default(),
                align_config: RnaReadAlignConfig::default(),
                report_id: Some("rna_reads_progress".to_string()),
                report_mode: RnaReadReportMode::Full,
                checkpoint_path: None,
                checkpoint_every_reads: 10_000,
                resume_from_checkpoint: false,
            },
            |progress| {
                if let OperationProgress::RnaReadInterpret(p) = progress {
                    progress_events.push(p);
                }
                true
            },
        )
        .expect("interpret reads with progress");

    assert!(!progress_events.is_empty());
    assert_eq!(progress_events.first().map(|p| p.reads_processed), Some(0));
    assert!(progress_events.iter().any(|p| p.reads_processed >= 1000));
    assert!(
        progress_events.iter().any(|p| {
            !p.done
                && p.reads_total == 0
                && p.reads_processed >= RNA_READ_PROGRESS_UPDATE_EVERY_READS
                && p.score_density_bins.iter().any(|count| *count > 0)
        }),
        "expected non-final streaming progress event with non-empty score density bins"
    );
    assert!(
        progress_events.iter().any(|p| {
            !p.done && p.reads_total == 0 && p.input_bytes_total > 0 && p.input_bytes_processed > 0
        }),
        "expected byte-level streaming progress events while total reads are still unknown"
    );
    let last = progress_events.last().expect("last progress");
    assert!(last.done);
    assert_eq!(last.reads_processed, 1001);
    assert_eq!(last.reads_total, 1001);
    assert!(last.input_bytes_total > 0);
    assert_eq!(last.input_bytes_processed, last.input_bytes_total);
    assert!(!last.bins.is_empty());
    assert!(
        last.bins
            .iter()
            .any(|bin| bin.confirmed_plus > 0 || bin.confirmed_minus > 0)
    );
    assert!(!last.score_density_bins.is_empty());
    assert!(last.score_density_bins.iter().any(|count| *count > 0));
}

#[test]
fn test_interpret_rna_reads_cancel_check_stops_quickly() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_cancel.fa");
    let mut fasta = String::new();
    for idx in 0..200usize {
        fasta.push_str(&format!(">read_{}\n{}\n", idx + 1, read_sequence));
    }
    fs::write(&input_path, fasta).expect("write reads");

    let mut progress_events = 0usize;
    let mut checks = 0usize;
    let result = engine.compute_rna_read_report_with_progress_and_cancel(
        "seq_a",
        feature_id,
        RnaReadInterpretationProfile::NanoporeCdnaV1,
        input_path.to_str().expect("path"),
        RnaReadInputFormat::Fasta,
        SplicingScopePreset::AllOverlappingBothStrands,
        RnaReadOriginMode::SingleGene,
        &[],
        false,
        &RnaReadSeedFilterConfig::default(),
        &RnaReadAlignConfig::default(),
        Some("rna_reads_cancel"),
        &mut |_progress| {
            progress_events = progress_events.saturating_add(1);
            true
        },
        &mut || {
            checks = checks.saturating_add(1);
            checks <= 2
        },
    );
    let err = result.expect_err("cancellation should abort run");
    assert!(err.message.to_ascii_lowercase().contains("cancel"));
    assert!(progress_events <= 3);
}

#[test]
fn test_rna_read_seed_hit_metrics_threshold_boundaries() {
    let (fraction, perfect, passed) = GentleEngine::seed_hit_metrics(10, 10, 0.30);
    assert!((fraction - 1.0).abs() < f64::EPSILON);
    assert!(perfect);
    assert!(passed);

    let (fraction, perfect, passed) = GentleEngine::seed_hit_metrics(10, 3, 0.30);
    assert!((fraction - 0.30).abs() < 1e-12);
    assert!(!perfect);
    assert!(passed);

    let (fraction, perfect, passed) = GentleEngine::seed_hit_metrics(10, 2, 0.30);
    assert!((fraction - 0.20).abs() < 1e-12);
    assert!(!perfect);
    assert!(!passed);

    let (fraction, perfect, passed) = GentleEngine::seed_hit_metrics(0, 0, 0.30);
    assert_eq!(fraction, 0.0);
    assert!(!perfect);
    assert!(!passed);
}

#[test]
fn test_rna_read_seed_filter_passes_uses_composite_gate() {
    let cfg = RnaReadSeedFilterConfig {
        min_seed_hit_fraction: 0.30,
        min_weighted_seed_hit_fraction: 0.05,
        min_unique_matched_kmers: 12,
        min_chain_consistency_fraction: 0.40,
        min_confirmed_exon_transitions: 1,
        min_transition_support_fraction: 0.20,
        ..RnaReadSeedFilterConfig::default()
    };
    assert!(GentleEngine::seed_filter_passes(
        0.30, 0.05, 40, 12, 0.90, 1.0, 8, 2, 4, &cfg
    ));
    assert!(!GentleEngine::seed_filter_passes(
        0.29, 0.90, 40, 99, 0.95, 1.0, 20, 5, 7, &cfg
    ));
    assert!(!GentleEngine::seed_filter_passes(
        0.30, 0.07, 40, 2, 0.80, 2.0, 12, 2, 4, &cfg
    ));
    assert!(!GentleEngine::seed_filter_passes(
        0.35, 0.01, 40, 13, 0.80, 3.5, 14, 2, 4, &cfg
    ));
    assert!(!GentleEngine::seed_filter_passes(
        0.30, 0.02, 40, 11, 0.80, 1.0, 10, 2, 4, &cfg
    ));
    assert!(!GentleEngine::seed_filter_passes(
        0.30, 0.09, 40, 9, 0.80, 5.0, 30, 2, 4, &cfg
    ));
    assert!(!GentleEngine::seed_filter_passes(
        0.30, 0.09, 40, 19, 0.80, 1.0, 30, 0, 8, &cfg
    ));
    assert!(!GentleEngine::seed_filter_passes(
        0.30, 0.09, 40, 19, 0.80, 1.0, 30, 1, 10, &cfg
    ));
    assert!(!GentleEngine::seed_filter_passes(
        0.90, 0.70, 40, 40, 0.10, 1.0, 39, 4, 4, &cfg
    ));
    assert!(GentleEngine::seed_filter_passes(
        0.95, 0.60, 8, 8, 1.0, 1.0, 7, 1, 1, &cfg
    ));
}

#[test]
fn test_classify_rna_read_origin_target_coherent_when_strict_gate_passes() {
    let cfg = RnaReadSeedFilterConfig::default();
    let spacing = SeedChainSpacingMetrics {
        support_kmers: 120,
        support_fraction: 0.82,
        transcript_id: "NM_005427.4".to_string(),
        ..SeedChainSpacingMetrics::default()
    };
    let path = ReadExonPathInference {
        transcript_id: "NM_005427.4".to_string(),
        strand: "+".to_string(),
        confirmed_transitions: 3,
        total_transitions: 4,
        strand_diagnostics: RnaReadStrandAssignmentDiagnostics {
            selected_strand: "+".to_string(),
            selected_exon_hits: 9,
            selected_transition_hits: 3,
            ..RnaReadStrandAssignmentDiagnostics::default()
        },
        ..ReadExonPathInference::default()
    };
    let classified =
        GentleEngine::classify_rna_read_origin(0.72, 0.66, true, &spacing, &path, "+", &cfg);
    assert_eq!(classified.origin_class, RnaReadOriginClass::TargetCoherent);
    assert!(classified.origin_confidence > 0.5);
    assert!((classified.strand_confidence - 1.0).abs() < f64::EPSILON);
}

#[test]
fn test_classify_rna_read_origin_marks_reverse_strand_local_block() {
    let cfg = RnaReadSeedFilterConfig::default();
    let spacing = SeedChainSpacingMetrics {
        support_kmers: 14,
        support_fraction: 0.21,
        transcript_id: "NM_reverse".to_string(),
        ..SeedChainSpacingMetrics::default()
    };
    let path = ReadExonPathInference {
        transcript_id: "NM_reverse".to_string(),
        strand: "-".to_string(),
        confirmed_transitions: 0,
        total_transitions: 0,
        strand_diagnostics: RnaReadStrandAssignmentDiagnostics {
            selected_strand: "-".to_string(),
            selected_exon_hits: 2,
            selected_transition_hits: 0,
            ..RnaReadStrandAssignmentDiagnostics::default()
        },
        ..ReadExonPathInference::default()
    };
    let classified =
        GentleEngine::classify_rna_read_origin(0.41, 0.18, false, &spacing, &path, "+", &cfg);
    assert_eq!(
        classified.origin_class,
        RnaReadOriginClass::RoiReverseStrandLocalBlock
    );
}

#[test]
fn test_build_rna_read_origin_candidates_collects_ranked_roles() {
    let path = ReadExonPathInference {
        transcript_id: "tx_selected".to_string(),
        strand: "+".to_string(),
        confirmed_transitions: 2,
        total_transitions: 4,
        strand_diagnostics: RnaReadStrandAssignmentDiagnostics {
            selected_exon_hits: 6,
            plus_best_transcript_id: "tx_selected".to_string(),
            plus_best_transition_hits: 2,
            plus_best_exon_hits: 6,
            minus_best_transcript_id: "tx_minus".to_string(),
            minus_best_transition_hits: 1,
            minus_best_exon_hits: 4,
            ..RnaReadStrandAssignmentDiagnostics::default()
        },
        ..ReadExonPathInference::default()
    };
    let spacing = SeedChainSpacingMetrics {
        transcript_id: "tx_chain".to_string(),
        support_kmers: 55,
        support_fraction: 0.64,
        ..SeedChainSpacingMetrics::default()
    };
    let transcript_models_by_id = HashMap::from([
        (
            "tx_chain".to_string(),
            TranscriptExonPathModel {
                transcript_id: "tx_chain".to_string(),
                strand: "-".to_string(),
                ..TranscriptExonPathModel::default()
            },
        ),
        (
            "tx_selected".to_string(),
            TranscriptExonPathModel {
                transcript_id: "tx_selected".to_string(),
                strand: "+".to_string(),
                ..TranscriptExonPathModel::default()
            },
        ),
        (
            "tx_minus".to_string(),
            TranscriptExonPathModel {
                transcript_id: "tx_minus".to_string(),
                strand: "-".to_string(),
                ..TranscriptExonPathModel::default()
            },
        ),
    ]);
    let candidates =
        GentleEngine::build_rna_read_origin_candidates(&path, &spacing, &transcript_models_by_id);
    let roles = candidates
        .iter()
        .map(|row| row.candidate_role.as_str())
        .collect::<HashSet<_>>();
    assert!(roles.contains("selected_path"));
    assert!(roles.contains("minus_best"));
    assert!(roles.contains("seed_chain_best"));
    assert_eq!(
        candidates
            .iter()
            .filter(|row| row.transcript_id == "tx_selected")
            .count(),
        1
    );
}

#[test]
fn test_compute_seed_chain_spacing_metrics_reports_median_transcript_gap() {
    let templates = vec![SplicingTranscriptTemplate {
        transcript_feature_id: 1,
        transcript_id: "tx1".to_string(),
        transcript_label: "tx1".to_string(),
        strand: "+".to_string(),
        ..SplicingTranscriptTemplate::default()
    }];
    let seed_template_positions = HashMap::from([
        (
            11_u32,
            vec![SeedTemplatePosition {
                template_idx: 0,
                template_pos: 100,
            }],
        ),
        (
            12_u32,
            vec![SeedTemplatePosition {
                template_idx: 0,
                template_pos: 101,
            }],
        ),
        (
            13_u32,
            vec![SeedTemplatePosition {
                template_idx: 0,
                template_pos: 102,
            }],
        ),
        (
            14_u32,
            vec![SeedTemplatePosition {
                template_idx: 0,
                template_pos: 103,
            }],
        ),
    ]);
    let contiguous_observations = vec![
        SeedMatchObservation {
            read_start: 0,
            bits: 11,
        },
        SeedMatchObservation {
            read_start: 1,
            bits: 12,
        },
        SeedMatchObservation {
            read_start: 2,
            bits: 13,
        },
        SeedMatchObservation {
            read_start: 3,
            bits: 14,
        },
    ];
    let contiguous = GentleEngine::compute_seed_chain_spacing_metrics(
        &contiguous_observations,
        &seed_template_positions,
        &templates,
    );
    assert_eq!(contiguous.transcript_id, "tx1");
    assert_eq!(contiguous.support_kmers, 4);
    assert_eq!(contiguous.transcript_gap_count, 3);
    assert!((contiguous.median_transcript_gap - 1.0).abs() < f64::EPSILON);

    let sparse_positions = HashMap::from([
        (
            11_u32,
            vec![SeedTemplatePosition {
                template_idx: 0,
                template_pos: 100,
            }],
        ),
        (
            12_u32,
            vec![SeedTemplatePosition {
                template_idx: 0,
                template_pos: 105,
            }],
        ),
        (
            13_u32,
            vec![SeedTemplatePosition {
                template_idx: 0,
                template_pos: 110,
            }],
        ),
        (
            14_u32,
            vec![SeedTemplatePosition {
                template_idx: 0,
                template_pos: 115,
            }],
        ),
    ]);
    let sparse_observations = vec![
        SeedMatchObservation {
            read_start: 0,
            bits: 11,
        },
        SeedMatchObservation {
            read_start: 5,
            bits: 12,
        },
        SeedMatchObservation {
            read_start: 10,
            bits: 13,
        },
        SeedMatchObservation {
            read_start: 15,
            bits: 14,
        },
    ];
    let sparse = GentleEngine::compute_seed_chain_spacing_metrics(
        &sparse_observations,
        &sparse_positions,
        &templates,
    );
    assert_eq!(sparse.support_kmers, 4);
    assert_eq!(sparse.transcript_gap_count, 3);
    assert!((sparse.median_transcript_gap - 5.0).abs() < f64::EPSILON);
}

fn deterministic_dna_sequence(mut state: u64, len: usize) -> Vec<u8> {
    let mut out = Vec::<u8>::with_capacity(len);
    for _ in 0..len {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        let base = match ((state >> 62) & 0x3) as u8 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        out.push(base);
    }
    out
}

fn delete_fraction_window(sequence: &[u8], start: usize, fraction: f64) -> Vec<u8> {
    if sequence.is_empty() {
        return vec![];
    }
    let delete_len = ((sequence.len() as f64) * fraction).round() as usize;
    let delete_len = delete_len.clamp(1, sequence.len().saturating_sub(1).max(1));
    let start = start.min(sequence.len().saturating_sub(delete_len));
    let mut out = Vec::<u8>::with_capacity(sequence.len().saturating_sub(delete_len));
    out.extend_from_slice(&sequence[..start]);
    out.extend_from_slice(&sequence[start + delete_len..]);
    out
}

#[test]
fn test_tp73_seed_filter_keeps_30pct_deleted_subsequences_and_rejects_random() {
    let mut engine = GentleEngine::default();
    engine
        .apply(Operation::LoadFile {
            path: "test_files/tp73.ncbi.gb".to_string(),
            as_id: Some("tp73".to_string()),
        })
        .expect("load tp73 fixture");
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;
    let (feature_id, dna_len) = {
        let dna = engine
            .state()
            .sequences
            .get("tp73")
            .expect("tp73 sequence present");
        let feature_id = dna
            .features()
            .iter()
            .position(GentleEngine::is_mrna_feature)
            .expect("tp73 mRNA feature");
        (feature_id, dna.len())
    };
    let splicing = engine
        .build_splicing_expert_view(
            "tp73",
            feature_id,
            SplicingScopePreset::TargetGroupTargetStrand,
        )
        .expect("build tp73 splicing view");
    assert!(
        !splicing.transcripts.is_empty(),
        "tp73 fixture should provide transcript templates"
    );
    let template = {
        let dna = engine
            .state()
            .sequences
            .get("tp73")
            .expect("tp73 sequence present");
        GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], seed_filter.kmer_len)
    };
    assert!(
        template.sequence.len() > 120,
        "expected sufficiently long tp73 transcript template"
    );
    let seed_catalog_rows = GentleEngine::collect_rna_seed_hash_catalog_rows(
        std::slice::from_ref(&template),
        seed_filter.kmer_len,
    );
    let seed_index = seed_catalog_rows
        .iter()
        .map(|row| row.seed_bits)
        .collect::<HashSet<_>>();
    let mut bins = GentleEngine::build_rna_read_seed_histogram_bins(dna_len);
    let histogram = GentleEngine::build_rna_read_seed_histogram_index(
        std::slice::from_ref(&template),
        dna_len,
        &bins,
    );
    let seed_template_positions =
        GentleEngine::build_seed_template_position_index(std::slice::from_ref(&template));
    let seed_occurrence_counts = HashMap::<u32, usize>::new();
    let mut score_read = |read: &[u8]| -> (bool, f64, f64, usize) {
        let (tested, matched, matched_bits, matched_observations) =
            GentleEngine::count_seed_hits_in_window_with_histogram(
                read,
                0,
                seed_filter.kmer_len,
                seed_filter.seed_stride_bp,
                &seed_index,
                &histogram,
                &mut bins,
                &seed_occurrence_counts,
            );
        let (fraction, _, _) =
            GentleEngine::seed_hit_metrics(tested, matched, seed_filter.min_seed_hit_fraction);
        let weighted = if tested == 0 {
            0.0
        } else {
            GentleEngine::weighted_seed_support_from_occurrences(
                &matched_bits,
                &seed_occurrence_counts,
            ) / tested as f64
        };
        let spacing = GentleEngine::compute_seed_chain_spacing_metrics(
            &matched_observations,
            &seed_template_positions,
            std::slice::from_ref(&template),
        );
        let pass = GentleEngine::seed_filter_passes(
            fraction,
            weighted,
            tested,
            matched_bits.len(),
            spacing.support_fraction,
            spacing.median_transcript_gap,
            spacing.transcript_gap_count,
            0,
            0,
            &seed_filter,
        );
        (pass, fraction, weighted, matched_bits.len())
    };

    let base_len = template.sequence.len().min(1500);
    let base_start = (template.sequence.len().saturating_sub(base_len)) / 2;
    let base_fragment = template.sequence[base_start..base_start + base_len].to_vec();
    let deletion_starts = [
        0usize,
        base_fragment.len() / 4,
        base_fragment.len() / 2,
        (base_fragment.len() * 3) / 4,
    ];
    for (idx, start) in deletion_starts.iter().copied().enumerate() {
        let subseq = delete_fraction_window(&base_fragment, start, 0.30);
        let (pass, fraction, weighted, unique_bits) = score_read(&subseq);
        assert!(
            pass,
            "tp73-derived subsequence #{idx} should pass seed filter: raw={fraction:.3} weighted={weighted:.3} unique={unique_bits}"
        );
    }

    let mut low_complexity = vec![b'T'; base_fragment.len()];
    if base_fragment.len() > 32 {
        low_complexity[base_fragment.len() / 2] = b'C';
    }
    let (pass, fraction, weighted, unique_bits) = score_read(&low_complexity);
    assert!(
        !pass,
        "low-complexity sequence unexpectedly passed: raw={fraction:.3} weighted={weighted:.3} unique={unique_bits}"
    );

    for idx in 0..12usize {
        let random = deterministic_dna_sequence(0xDEADBEEF_u64 + idx as u64, base_fragment.len());
        let (pass, fraction, weighted, unique_bits) = score_read(&random);
        assert!(
            !pass,
            "random sequence #{idx} unexpectedly passed: raw={fraction:.3} weighted={weighted:.3} unique={unique_bits}"
        );
    }
}

#[test]
fn test_tp73_seed_filter_cross_species_and_tp53_specificity_sets() {
    let mut engine = GentleEngine::default();
    engine
        .apply(Operation::LoadFile {
            path: "test_files/tp73.ncbi.gb".to_string(),
            as_id: Some("tp73".to_string()),
        })
        .expect("load tp73 fixture");
    let seed_filter = RnaReadSeedFilterConfig::default();
    let (feature_id, dna_len) = {
        let dna = engine
            .state()
            .sequences
            .get("tp73")
            .expect("tp73 sequence present");
        let feature_id = dna
            .features()
            .iter()
            .position(GentleEngine::is_mrna_feature)
            .expect("tp73 mRNA feature");
        (feature_id, dna.len())
    };
    let splicing = engine
        .build_splicing_expert_view(
            "tp73",
            feature_id,
            SplicingScopePreset::TargetGroupTargetStrand,
        )
        .expect("build tp73 splicing view");
    assert!(
        !splicing.transcripts.is_empty(),
        "tp73 fixture should provide transcript templates"
    );
    let templates = {
        let dna = engine
            .state()
            .sequences
            .get("tp73")
            .expect("tp73 sequence present");
        splicing
            .transcripts
            .iter()
            .map(|tx| GentleEngine::make_transcript_template(dna, tx, seed_filter.kmer_len))
            .collect::<Vec<_>>()
    };
    let seed_catalog_rows =
        GentleEngine::collect_rna_seed_hash_catalog_rows(&templates, seed_filter.kmer_len);
    assert!(
        !seed_catalog_rows.is_empty(),
        "tp73 seed catalog should not be empty"
    );
    let seed_index = seed_catalog_rows
        .iter()
        .map(|row| row.seed_bits)
        .collect::<HashSet<_>>();
    let mut seed_occurrence_counts = HashMap::<u32, usize>::new();
    for row in &seed_catalog_rows {
        *seed_occurrence_counts.entry(row.seed_bits).or_insert(0) += 1;
    }
    let seed_template_positions = GentleEngine::build_seed_template_position_index(&templates);
    let seed_support_exons =
        GentleEngine::collect_seed_support_exon_summaries(&splicing.transcripts);
    let (seed_to_exons, seed_to_transitions, transcript_models, _transition_rows) =
        GentleEngine::build_seed_support_indexes(
            &seed_support_exons,
            &templates,
            seed_filter.kmer_len,
        );
    let mut bins = GentleEngine::build_rna_read_seed_histogram_bins(dna_len);
    let histogram = GentleEngine::build_rna_read_seed_histogram_index(&templates, dna_len, &bins);

    let score_read = |sequence: &[u8],
                      bins: &mut [RnaReadSeedHistogramBin]|
     -> (bool, f64, f64, usize, usize, usize, usize) {
        let (normalized, _) =
            GentleEngine::normalize_rna_read_sequence_for_scoring(sequence, &seed_filter);
        let windows = GentleEngine::full_read_hash_windows(normalized.len());
        let mut tested_kmers = 0usize;
        let mut matched_kmers = 0usize;
        let mut matched_seed_bits = HashSet::<u32>::new();
        let mut matched_seed_observations = Vec::<SeedMatchObservation>::new();
        for (start, end) in windows {
            let (tested, matched, matched_bits, matched_observations) =
                GentleEngine::count_seed_hits_in_window_with_histogram(
                    &normalized[start..end],
                    start,
                    seed_filter.kmer_len,
                    seed_filter.seed_stride_bp,
                    &seed_index,
                    &histogram,
                    bins,
                    &seed_occurrence_counts,
                );
            tested_kmers = tested_kmers.saturating_add(tested);
            matched_kmers = matched_kmers.saturating_add(matched);
            matched_seed_bits.extend(matched_bits);
            matched_seed_observations.extend(matched_observations);
        }
        let (fraction, _, _) = GentleEngine::seed_hit_metrics(
            tested_kmers,
            matched_kmers,
            seed_filter.min_seed_hit_fraction,
        );
        let weighted = if tested_kmers == 0 {
            0.0
        } else {
            GentleEngine::weighted_seed_support_from_occurrences(
                &matched_seed_bits,
                &seed_occurrence_counts,
            ) / tested_kmers as f64
        };
        let spacing = GentleEngine::compute_seed_chain_spacing_metrics(
            &matched_seed_observations,
            &seed_template_positions,
            &templates,
        );
        let mut supported_exons = HashSet::<usize>::new();
        let mut supported_transitions = HashSet::<(usize, usize)>::new();
        for bits in &matched_seed_bits {
            if let Some(exons) = seed_to_exons.get(bits) {
                supported_exons.extend(exons.iter().copied());
            }
            if let Some(transitions) = seed_to_transitions.get(bits) {
                supported_transitions.extend(transitions.iter().copied());
            }
        }
        let path_inference = GentleEngine::infer_read_exon_path(
            &transcript_models,
            &supported_exons,
            &supported_transitions,
            &spacing.transcript_id,
        );
        let pass = GentleEngine::seed_filter_passes(
            fraction,
            weighted,
            tested_kmers,
            matched_seed_bits.len(),
            spacing.support_fraction,
            spacing.median_transcript_gap,
            spacing.transcript_gap_count,
            path_inference.confirmed_transitions,
            path_inference.total_transitions,
            &seed_filter,
        );
        (
            pass,
            fraction,
            weighted,
            matched_seed_bits.len(),
            spacing.transcript_gap_count,
            path_inference.confirmed_transitions,
            path_inference.total_transitions,
        )
    };

    let evaluate_set = |path: &str,
                        bins: &mut [RnaReadSeedHistogramBin]|
     -> (usize, usize, f64, f64, Vec<String>) {
        let records =
            GentleEngine::parse_fasta_records_with_offsets(path).expect("parse fasta set");
        assert!(!records.is_empty(), "fasta set should not be empty: {path}");
        let mut passed = 0usize;
        let mut best_raw = 0.0_f64;
        let mut best_weighted = 0.0_f64;
        let mut failed = Vec::<String>::new();
        for record in &records {
            let (pass, raw, weighted, unique, gap_n, confirmed, total) =
                score_read(&record.sequence, bins);
            if pass {
                passed = passed.saturating_add(1);
            } else {
                failed.push(format!(
                            "{} raw={raw:.3} weighted={weighted:.4} unique={unique} gap_n={gap_n} confirmed={confirmed}/{total}",
                            record.header_id
                        ));
            }
            best_raw = best_raw.max(raw);
            best_weighted = best_weighted.max(weighted);
        }
        (records.len(), passed, best_raw, best_weighted, failed)
    };

    let (chimp_total, chimp_passed, chimp_best_raw, _chimp_best_weighted, chimp_failed) =
        evaluate_set(
            "test_files/fixtures/mapping/ensembl_chimp_tp73_all.fasta",
            &mut bins,
        );
    let (human_total, human_passed, human_best_raw, _human_best_weighted, human_failed) =
        evaluate_set(
            "test_files/fixtures/mapping/ensembl_human_tp73_all.fasta",
            &mut bins,
        );
    let (tp53_total, tp53_passed, tp53_best_raw, tp53_best_weighted, tp53_failed) = evaluate_set(
        "test_files/fixtures/mapping/ensembl_human_tp53_all.fasta",
        &mut bins,
    );
    let (mouse_total, _mouse_passed, _mouse_best_raw, _mouse_best_weighted, _mouse_failed) =
        evaluate_set(
            "test_files/fixtures/mapping/ensembl_mouse_trp73_all.fasta",
            &mut bins,
        );
    assert!(mouse_total > 0, "mouse trp73 set should contain reads");

    assert_eq!(
        chimp_passed, chimp_total,
        "chimp TP73 sequences should all pass seed filter (passed={chimp_passed}/{chimp_total}, best_raw={chimp_best_raw:.3}, failed={chimp_failed:?})"
    );
    assert_eq!(
        human_passed, human_total,
        "human TP73 sequences should all pass seed filter (passed={human_passed}/{human_total}, best_raw={human_best_raw:.3}, failed={human_failed:?})"
    );
    assert_eq!(
        tp53_passed, 0,
        "human TP53 sequences should be rejected (passed={tp53_passed}/{tp53_total}, best_raw={tp53_best_raw:.3}, best_weighted={tp53_best_weighted:.4}, failed={tp53_failed:?})"
    );
}

#[test]
fn test_infer_read_exon_path_prefers_chain_strand_on_cross_strand_tie() {
    let transcript_models = vec![
        TranscriptExonPathModel {
            transcript_feature_id: 10,
            transcript_id: "tx_plus".to_string(),
            transcript_label: "TX+".to_string(),
            strand: "+".to_string(),
            exon_ordinals: vec![1, 2],
            transitions: vec![(1, 2)],
        },
        TranscriptExonPathModel {
            transcript_feature_id: 5,
            transcript_id: "tx_minus".to_string(),
            transcript_label: "TX-".to_string(),
            strand: "-".to_string(),
            exon_ordinals: vec![1, 2],
            transitions: vec![(1, 2)],
        },
    ];
    let supported_exons = [1usize, 2usize].into_iter().collect::<HashSet<_>>();
    let supported_transitions = [(1usize, 2usize)].into_iter().collect::<HashSet<_>>();

    let inferred = GentleEngine::infer_read_exon_path(
        &transcript_models,
        &supported_exons,
        &supported_transitions,
        "tx_plus",
    );
    assert_eq!(inferred.transcript_id, "tx_plus");
    assert_eq!(inferred.strand, "+");
    assert_eq!(inferred.confirmed_transitions, 1);
    assert_eq!(inferred.total_transitions, 1);
    assert_eq!(inferred.strand_diagnostics.chain_preferred_strand, "+");
    assert!(inferred.strand_diagnostics.competing_opposite_strand);
    assert!(inferred.strand_diagnostics.ambiguous_near_tie);
    assert_eq!(inferred.strand_diagnostics.selected_strand, "+");
}

#[test]
fn test_full_read_hash_windows_phase1_hashes_full_read() {
    assert_eq!(
        GentleEngine::full_read_hash_windows(0),
        Vec::<(usize, usize)>::new()
    );
    assert_eq!(GentleEngine::full_read_hash_windows(37), vec![(0, 37)]);
    assert_eq!(GentleEngine::full_read_hash_windows(420), vec![(0, 420)]);
    assert_eq!(GentleEngine::full_read_hash_windows(1200), vec![(0, 1200)]);
}

#[test]
fn test_rna_read_seed_filter_config_ignores_removed_legacy_window_fields_in_json() {
    let config: RnaReadSeedFilterConfig = serde_json::from_value(serde_json::json!({
        "kmer_len": 9,
        "seed_stride_bp": 1,
        "short_full_hash_max_bp": 420,
        "long_window_bp": 140,
        "long_window_count": 3,
        "min_seed_hit_fraction": 0.30,
        "min_weighted_seed_hit_fraction": 0.05,
        "min_unique_matched_kmers": 12,
        "max_median_transcript_gap": 4.0,
        "min_chain_consistency_fraction": 0.40,
        "min_confirmed_exon_transitions": 1,
        "min_transition_support_fraction": 0.05,
        "cdna_poly_t_flip_enabled": true,
        "poly_t_prefix_min_bp": 18
    }))
    .expect("legacy JSON with removed window fields should still deserialize");
    assert_eq!(config.kmer_len, 9);
    assert_eq!(config.seed_stride_bp, 1);
    assert!((config.min_seed_hit_fraction - 0.30).abs() < f64::EPSILON);
}

#[test]
fn test_summarize_read_lengths_reports_mean_median_and_p95() {
    let mut counts = vec![0u64; 12];
    for len in [4usize, 4, 5, 8, 8, 10] {
        GentleEngine::update_read_length_counts(&mut counts, len);
    }
    let total_reads = 6usize;
    let total_bases = (4 + 4 + 5 + 8 + 8 + 10) as u64;
    let (mean, median, p95) =
        GentleEngine::summarize_read_lengths(&counts, total_reads, total_bases);
    assert!((mean - (39.0 / 6.0)).abs() < 1e-12);
    assert_eq!(median, 5);
    assert_eq!(p95, 10);
}

#[test]
fn test_rna_read_sequence_scoring_from_seed_index_is_deterministic() {
    let read = b"ACGTTGCAACGT";
    let kmer_len = 3usize;
    let mut bins = GentleEngine::build_rna_read_seed_histogram_bins(12);
    let histogram_index: HashMap<u32, Vec<SeedHistogramWeight>> = HashMap::new();
    let seed_occurrence_counts = HashMap::<u32, usize>::new();

    let mut all_seed_bits = HashSet::new();
    for start in 0..=read.len() - kmer_len {
        let bits = GentleEngine::encode_kmer_bits(&read[start..start + kmer_len]).expect("bits");
        all_seed_bits.insert(bits);
    }
    assert!(!all_seed_bits.is_empty());

    let (tested_all, matched_all, matched_bits_all, matched_obs_all) =
        GentleEngine::count_seed_hits_in_window_with_histogram(
            read,
            0,
            kmer_len,
            1,
            &all_seed_bits,
            &histogram_index,
            &mut bins,
            &seed_occurrence_counts,
        );
    assert!(tested_all > 0);
    assert_eq!(matched_all, tested_all);
    assert!(!matched_bits_all.is_empty());
    assert!(!matched_obs_all.is_empty());
    let (fraction_all, perfect_all, passed_all) =
        GentleEngine::seed_hit_metrics(tested_all, matched_all, 0.30);
    assert!((fraction_all - 1.0).abs() < 1e-12);
    assert!(perfect_all);
    assert!(passed_all);

    let empty_seed_index = HashSet::new();
    let (tested_none, matched_none, matched_bits_none, matched_obs_none) =
        GentleEngine::count_seed_hits_in_window_with_histogram(
            read,
            0,
            kmer_len,
            1,
            &empty_seed_index,
            &histogram_index,
            &mut bins,
            &seed_occurrence_counts,
        );
    assert_eq!(tested_none, tested_all);
    assert_eq!(matched_none, 0);
    assert!(matched_bits_none.is_empty());
    assert!(matched_obs_none.is_empty());
    let (fraction_none, perfect_none, passed_none) =
        GentleEngine::seed_hit_metrics(tested_none, matched_none, 0.30);
    assert_eq!(fraction_none, 0.0);
    assert!(!perfect_none);
    assert!(!passed_none);
}

#[test]
fn test_rna_read_seed_stride_reduces_tested_kmers() {
    let read = b"ACGTACGT";
    let kmer_len = 2usize;
    let mut bins = GentleEngine::build_rna_read_seed_histogram_bins(read.len());
    let histogram_index: HashMap<u32, Vec<SeedHistogramWeight>> = HashMap::new();
    let seed_occurrence_counts = HashMap::<u32, usize>::new();
    let mut all_seed_bits = HashSet::new();
    for start in 0..=read.len() - kmer_len {
        let bits = GentleEngine::encode_kmer_bits(&read[start..start + kmer_len]).expect("bits");
        all_seed_bits.insert(bits);
    }

    let (tested_dense, matched_dense, _, matched_obs_dense) =
        GentleEngine::count_seed_hits_in_window_with_histogram(
            read,
            0,
            kmer_len,
            1,
            &all_seed_bits,
            &histogram_index,
            &mut bins,
            &seed_occurrence_counts,
        );
    let (tested_sparse, matched_sparse, _, matched_obs_sparse) =
        GentleEngine::count_seed_hits_in_window_with_histogram(
            read,
            0,
            kmer_len,
            2,
            &all_seed_bits,
            &histogram_index,
            &mut bins,
            &seed_occurrence_counts,
        );

    assert_eq!(tested_dense, 7);
    assert_eq!(matched_dense, 7);
    assert_eq!(tested_sparse, 4);
    assert_eq!(matched_sparse, 4);
    assert_eq!(
        matched_obs_sparse
            .iter()
            .map(|row| row.read_start)
            .collect::<Vec<_>>(),
        vec![0, 2, 4, 6]
    );
    assert_eq!(matched_obs_dense.len(), tested_dense);
}

#[test]
fn test_mapped_support_counts_follow_transcript_offsets_not_genomic_span_overlap() {
    let splicing = SplicingExpertView {
        seq_id: "seq".to_string(),
        target_feature_id: 1,
        group_label: "GENE1".to_string(),
        strand: "+".to_string(),
        region_start_1based: 1,
        region_end_1based: 25,
        transcript_count: 2,
        unique_exon_count: 3,
        instruction: String::new(),
        transcripts: vec![
            SplicingTranscriptLane {
                transcript_feature_id: 1,
                transcript_id: "tx_main".to_string(),
                label: "TX_MAIN".to_string(),
                strand: "+".to_string(),
                exons: vec![
                    SplicingRange {
                        start_1based: 1,
                        end_1based: 5,
                    },
                    SplicingRange {
                        start_1based: 21,
                        end_1based: 25,
                    },
                ],
                exon_cds_phases: vec![],
                introns: vec![SplicingRange {
                    start_1based: 6,
                    end_1based: 20,
                }],
                has_target_feature: true,
            },
            SplicingTranscriptLane {
                transcript_feature_id: 2,
                transcript_id: "tx_alt".to_string(),
                label: "TX_ALT".to_string(),
                strand: "+".to_string(),
                exons: vec![
                    SplicingRange {
                        start_1based: 1,
                        end_1based: 5,
                    },
                    SplicingRange {
                        start_1based: 11,
                        end_1based: 15,
                    },
                    SplicingRange {
                        start_1based: 21,
                        end_1based: 25,
                    },
                ],
                exon_cds_phases: vec![],
                introns: vec![
                    SplicingRange {
                        start_1based: 6,
                        end_1based: 10,
                    },
                    SplicingRange {
                        start_1based: 16,
                        end_1based: 20,
                    },
                ],
                has_target_feature: false,
            },
        ],
        unique_exons: vec![
            SplicingExonSummary {
                start_1based: 1,
                end_1based: 5,
                support_transcript_count: 2,
                constitutive: false,
            },
            SplicingExonSummary {
                start_1based: 11,
                end_1based: 15,
                support_transcript_count: 1,
                constitutive: false,
            },
            SplicingExonSummary {
                start_1based: 21,
                end_1based: 25,
                support_transcript_count: 2,
                constitutive: false,
            },
        ],
        matrix_rows: vec![],
        boundaries: vec![],
        junctions: vec![
            SplicingJunctionArc {
                donor_1based: 5,
                acceptor_1based: 11,
                support_transcript_count: 1,
                transcript_feature_ids: vec![2],
            },
            SplicingJunctionArc {
                donor_1based: 5,
                acceptor_1based: 21,
                support_transcript_count: 1,
                transcript_feature_ids: vec![1],
            },
            SplicingJunctionArc {
                donor_1based: 15,
                acceptor_1based: 21,
                support_transcript_count: 1,
                transcript_feature_ids: vec![2],
            },
        ],
        events: vec![],
    };
    let mapping = RnaReadMappingHit {
        transcript_feature_id: 1,
        transcript_id: "tx_main".to_string(),
        transcript_label: "TX_MAIN".to_string(),
        strand: "+".to_string(),
        target_start_1based: 1,
        target_end_1based: 25,
        target_start_offset_0based: 0,
        target_end_offset_0based_exclusive: 10,
        ..RnaReadMappingHit::default()
    };
    let mut exon_counts = vec![0usize; splicing.unique_exons.len()];
    let mut junction_counts = vec![0usize; splicing.junctions.len()];

    GentleEngine::accumulate_support_counts_for_mapping(
        &mapping,
        &splicing,
        &mut exon_counts,
        &mut junction_counts,
    );

    assert_eq!(exon_counts, vec![1, 0, 1]);
    assert_eq!(junction_counts, vec![0, 1, 0]);
}

#[test]
fn test_rna_read_retention_rank_prefers_aligned_hits_over_unaligned_seed_only_hits() {
    let aligned = RnaReadInterpretationHit {
        record_index: 2,
        header_id: "aligned_read".to_string(),
        sequence: "ACGTACGT".to_string(),
        read_length_bp: 8,
        tested_kmers: 6,
        matched_kmers: 2,
        seed_hit_fraction: 0.33,
        weighted_seed_hit_fraction: 0.20,
        weighted_matched_kmers: 1.2,
        passed_seed_filter: true,
        best_mapping: Some(RnaReadMappingHit {
            transcript_id: "tx_aligned".to_string(),
            transcript_label: "TX_ALIGNED".to_string(),
            strand: "+".to_string(),
            query_start_0based: 0,
            query_end_0based_exclusive: 8,
            target_start_1based: 101,
            target_end_1based: 108,
            score: 75,
            identity_fraction: 0.95,
            query_coverage_fraction: 0.95,
            ..RnaReadMappingHit::default()
        }),
        ..RnaReadInterpretationHit::default()
    };
    let unaligned_seed_heavy = RnaReadInterpretationHit {
        record_index: 1,
        header_id: "seed_heavy_unaligned".to_string(),
        sequence: "ACGTACGT".to_string(),
        read_length_bp: 8,
        tested_kmers: 6,
        matched_kmers: 6,
        seed_hit_fraction: 1.0,
        weighted_seed_hit_fraction: 1.0,
        weighted_matched_kmers: 6.0,
        passed_seed_filter: true,
        best_mapping: None,
        ..RnaReadInterpretationHit::default()
    };
    let mut hits = vec![unaligned_seed_heavy, aligned];
    GentleEngine::sort_rna_read_hits_by_retention_rank(&mut hits);
    assert_eq!(hits[0].header_id, "aligned_read");
    assert_eq!(hits[1].header_id, "seed_heavy_unaligned");
}

#[test]
fn test_inspect_and_export_rna_read_alignment_dotplot_follow_alignment_rank() {
    let mut engine = GentleEngine::default();
    let report = RnaReadInterpretationReport {
        schema: "gentle.rna_read_report.v1".to_string(),
        report_id: "rna_reads_alignment_rank".to_string(),
        seq_id: "seq_alignment".to_string(),
        align_config: RnaReadAlignConfig {
            min_identity_fraction: 0.6,
            max_secondary_mappings: 0,
            ..RnaReadAlignConfig::default()
        },
        hits: vec![
            RnaReadInterpretationHit {
                record_index: 0,
                header_id: "aligned_hi_id".to_string(),
                sequence: "ACGTACGTACGT".to_string(),
                read_length_bp: 12,
                tested_kmers: 9,
                matched_kmers: 8,
                seed_hit_fraction: 0.88,
                weighted_seed_hit_fraction: 0.82,
                weighted_matched_kmers: 7.4,
                seed_chain_transcript_id: "tx_hi_id".to_string(),
                seed_chain_support_fraction: 1.0,
                seed_median_transcript_gap: 0.0,
                seed_transcript_gap_count: 1,
                exon_path_transcript_id: "tx_hi_id".to_string(),
                exon_path: "1:2".to_string(),
                exon_transitions_confirmed: 1,
                exon_transitions_total: 1,
                reverse_complement_applied: true,
                passed_seed_filter: true,
                best_mapping: Some(RnaReadMappingHit {
                    transcript_id: "tx_hi_id".to_string(),
                    transcript_label: "TX_HI_ID".to_string(),
                    strand: "+".to_string(),
                    query_start_0based: 0,
                    query_end_0based_exclusive: 12,
                    target_start_1based: 201,
                    target_end_1based: 212,
                    score: 150,
                    identity_fraction: 0.98,
                    query_coverage_fraction: 0.72,
                    ..RnaReadMappingHit::default()
                }),
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 1,
                header_id: "aligned_hi_cov".to_string(),
                sequence: "ACGTACGTACGT".to_string(),
                read_length_bp: 12,
                tested_kmers: 9,
                matched_kmers: 3,
                seed_hit_fraction: 0.33,
                weighted_seed_hit_fraction: 0.27,
                weighted_matched_kmers: 2.5,
                seed_chain_transcript_id: "tx_seed_cov".to_string(),
                exon_path_transcript_id: "tx_seed_cov".to_string(),
                exon_path: "1-2".to_string(),
                exon_transitions_confirmed: 0,
                exon_transitions_total: 1,
                passed_seed_filter: true,
                best_mapping: Some(RnaReadMappingHit {
                    transcript_id: "tx_hi_cov".to_string(),
                    transcript_label: "TX_HI_COV".to_string(),
                    strand: "+".to_string(),
                    query_start_0based: 0,
                    query_end_0based_exclusive: 12,
                    target_start_1based: 301,
                    target_end_1based: 312,
                    score: 90,
                    identity_fraction: 0.84,
                    query_coverage_fraction: 0.96,
                    ..RnaReadMappingHit::default()
                }),
                ..RnaReadInterpretationHit::default()
            },
            RnaReadInterpretationHit {
                record_index: 2,
                header_id: "seed_only".to_string(),
                sequence: "ACGTACGTACGT".to_string(),
                read_length_bp: 12,
                tested_kmers: 9,
                matched_kmers: 9,
                seed_hit_fraction: 1.0,
                weighted_seed_hit_fraction: 1.0,
                weighted_matched_kmers: 9.0,
                passed_seed_filter: true,
                best_mapping: None,
                ..RnaReadInterpretationHit::default()
            },
        ],
        ..RnaReadInterpretationReport::default()
    };
    engine
        .upsert_rna_read_report(report)
        .expect("upsert synthetic RNA-read report");

    let inspection = engine
        .inspect_rna_read_alignments("rna_reads_alignment_rank", RnaReadHitSelection::All, 10)
        .expect("inspect alignment rows");
    assert_eq!(inspection.aligned_count, 2);
    assert_eq!(inspection.row_count, 2);
    assert_eq!(inspection.rows[0].rank, 1);
    assert_eq!(inspection.rows[0].header_id, "aligned_hi_cov");
    assert_eq!(
        inspection.rows[0].alignment_effect,
        RnaReadAlignmentEffect::ReassignedTranscript
    );
    assert_eq!(inspection.rows[1].rank, 2);
    assert_eq!(inspection.rows[1].header_id, "aligned_hi_id");
    assert_eq!(
        inspection.rows[1].alignment_effect,
        RnaReadAlignmentEffect::ConfirmedAssignment
    );
    assert_eq!(inspection.rows[1].phase1_primary_transcript_id, "tx_hi_id");
    assert!(inspection.rows[1].reverse_complement_applied);

    let filtered = engine
        .inspect_rna_read_alignments_with_subset(
            "rna_reads_alignment_rank",
            RnaReadHitSelection::All,
            10,
            Some(RnaReadAlignmentInspectionSubsetSpec {
                effect_filter: RnaReadAlignmentInspectionEffectFilter::DisagreementOnly,
                sort_key: RnaReadAlignmentInspectionSortKey::Score,
                search: "aligned_hi_cov".to_string(),
                selected_record_indices: vec![1, 1],
                score_density_variant: RnaReadScoreDensityVariant::AllScored,
                score_density_seed_filter_override: None,
                score_bin_index: None,
                score_bin_count: 0,
            }),
        )
        .expect("inspect filtered alignment rows");
    assert_eq!(filtered.aligned_count, 2);
    assert_eq!(filtered.subset_match_count, 1);
    assert_eq!(filtered.row_count, 1);
    assert_eq!(
        filtered.subset_spec.effect_filter,
        RnaReadAlignmentInspectionEffectFilter::DisagreementOnly
    );
    assert_eq!(
        filtered.subset_spec.sort_key,
        RnaReadAlignmentInspectionSortKey::Score
    );
    assert_eq!(filtered.subset_spec.search, "aligned_hi_cov");
    assert_eq!(filtered.subset_spec.selected_record_indices, vec![1]);
    assert_eq!(filtered.subset_spec.score_bin_index, None);
    assert_eq!(filtered.rows[0].header_id, "aligned_hi_cov");
    assert_eq!(filtered.rows[0].rank, 1);

    let sorted = engine
        .inspect_rna_read_alignments_with_subset(
            "rna_reads_alignment_rank",
            RnaReadHitSelection::All,
            10,
            Some(RnaReadAlignmentInspectionSubsetSpec {
                effect_filter: RnaReadAlignmentInspectionEffectFilter::AllAligned,
                sort_key: RnaReadAlignmentInspectionSortKey::Score,
                search: String::new(),
                selected_record_indices: vec![],
                score_density_variant: RnaReadScoreDensityVariant::AllScored,
                score_density_seed_filter_override: None,
                score_bin_index: None,
                score_bin_count: 0,
            }),
        )
        .expect("inspect score-sorted alignment rows");
    assert_eq!(sorted.subset_match_count, 2);
    assert_eq!(sorted.rows[0].header_id, "aligned_hi_id");
    assert_eq!(sorted.rows[0].rank, 2);
    assert_eq!(sorted.rows[1].header_id, "aligned_hi_cov");
    assert_eq!(sorted.rows[1].rank, 1);

    let target_score_bin_index =
        ((inspection.rows[0].seed_hit_fraction.clamp(0.0, 1.0) * 40.0).floor() as usize).min(39);
    let score_bin_filtered = engine
        .inspect_rna_read_alignments_with_subset(
            "rna_reads_alignment_rank",
            RnaReadHitSelection::All,
            10,
            Some(RnaReadAlignmentInspectionSubsetSpec {
                effect_filter: RnaReadAlignmentInspectionEffectFilter::AllAligned,
                sort_key: RnaReadAlignmentInspectionSortKey::Rank,
                search: String::new(),
                selected_record_indices: vec![],
                score_density_variant: RnaReadScoreDensityVariant::AllScored,
                score_density_seed_filter_override: None,
                score_bin_index: Some(target_score_bin_index),
                score_bin_count: 40,
            }),
        )
        .expect("inspect score-bin filtered alignment rows");
    assert_eq!(
        score_bin_filtered.subset_spec.score_bin_index,
        Some(target_score_bin_index)
    );
    assert_eq!(score_bin_filtered.subset_spec.score_bin_count, 40);
    assert!(score_bin_filtered.rows.iter().all(|row| {
        ((row.seed_hit_fraction.clamp(0.0, 1.0) * 40.0).floor() as usize).min(39)
            == target_score_bin_index
    }));

    let replay_target_score_bin_index = inspection
        .rows
        .iter()
        .find(|row| row.header_id == "aligned_hi_id")
        .map(|row| ((row.seed_hit_fraction.clamp(0.0, 1.0) * 40.0).floor() as usize).min(39))
        .expect("replay-filtered row");
    let replay_filter = RnaReadSeedFilterConfig {
        min_seed_hit_fraction: 0.75,
        min_weighted_seed_hit_fraction: 0.70,
        min_unique_matched_kmers: 8,
        min_chain_consistency_fraction: 0.90,
        max_median_transcript_gap: 4.0,
        min_confirmed_exon_transitions: 1,
        min_transition_support_fraction: 0.50,
        ..RnaReadSeedFilterConfig::default()
    };
    let replay_filtered = engine
        .inspect_rna_read_alignments_with_subset(
            "rna_reads_alignment_rank",
            RnaReadHitSelection::All,
            10,
            Some(RnaReadAlignmentInspectionSubsetSpec {
                effect_filter: RnaReadAlignmentInspectionEffectFilter::AllAligned,
                sort_key: RnaReadAlignmentInspectionSortKey::Rank,
                search: String::new(),
                selected_record_indices: vec![],
                score_density_variant: RnaReadScoreDensityVariant::RetainedReplayCurrentControls,
                score_density_seed_filter_override: Some(replay_filter),
                score_bin_index: Some(replay_target_score_bin_index),
                score_bin_count: 40,
            }),
        )
        .expect("inspect replay-filtered alignment rows");
    assert_eq!(replay_filtered.subset_match_count, 1);
    assert_eq!(replay_filtered.rows.len(), 1);
    assert_eq!(replay_filtered.rows[0].header_id, "aligned_hi_id");

    let td = tempdir().expect("tempdir");
    let tsv_path = td.path().join("alignment_rows.tsv");
    let tsv_export = engine
        .export_rna_read_alignments_tsv(
            "rna_reads_alignment_rank",
            tsv_path.to_str().expect("tsv path"),
            RnaReadHitSelection::All,
            Some(1),
            &[],
            None,
        )
        .expect("export alignment tsv");
    assert_eq!(tsv_export.row_count, 1);
    assert_eq!(tsv_export.aligned_count, 2);
    assert_eq!(tsv_export.limit, Some(1));
    let tsv_text = fs::read_to_string(&tsv_path).expect("read alignment tsv");
    assert!(tsv_text.contains("# report_id=rna_reads_alignment_rank"));
    assert!(tsv_text.contains(
        "selection=all selected_record_indices=none subset_spec=none limit=1 row_count=1 aligned_count=2"
    ));
    assert!(tsv_text.contains("# seed_filter: k=10 stride=1"));
    assert!(tsv_text.contains("adjacent windows overlap by 9 bp"));
    assert!(tsv_text.contains("# align_config: min_identity_fraction=0.60"));
    assert!(tsv_text.contains("alignment_effect"));
    assert!(tsv_text.contains("phase1_primary_transcript_id"));
    assert!(tsv_text.contains("mapped_exon_support"));
    assert!(tsv_text.contains("alignment_mode"));
    assert!(tsv_text.contains("aligned_hi_cov"));

    let selected_tsv_path = td.path().join("alignment_rows_selected.tsv");
    let selected_tsv_export = engine
        .export_rna_read_alignments_tsv(
            "rna_reads_alignment_rank",
            selected_tsv_path.to_str().expect("selected tsv path"),
            RnaReadHitSelection::All,
            None,
            &[0],
            Some("filter=selected only | sort=score | search=<none>"),
        )
        .expect("export selected alignment tsv");
    assert_eq!(selected_tsv_export.row_count, 1);
    assert_eq!(selected_tsv_export.aligned_count, 1);
    let selected_tsv_text =
        fs::read_to_string(&selected_tsv_path).expect("read selected alignment tsv");
    assert!(selected_tsv_text.contains("selected_record_indices=0"));
    assert!(
        selected_tsv_text.contains("subset_spec=filter=selected only | sort=score | search=<none>")
    );
    assert!(selected_tsv_text.contains("aligned_hi_id"));
    assert!(!selected_tsv_text.contains("aligned_hi_cov\t"));

    let svg_path = td.path().join("alignment_dotplot.svg");
    let export = engine
        .export_rna_read_alignment_dotplot_svg(
            "rna_reads_alignment_rank",
            svg_path.to_str().expect("svg path"),
            RnaReadHitSelection::All,
            1,
        )
        .expect("export alignment dotplot svg");
    assert_eq!(export.point_count, 2);
    assert_eq!(export.rendered_point_count, 1);
    let svg_text = fs::read_to_string(&svg_path).expect("read alignment dotplot svg");
    assert!(svg_text.contains("alignment dotplot"));
    assert!(svg_text.contains("rendered_points=1"));
    assert!(svg_text.contains("total_points=2"));
}

#[test]
fn test_export_rna_read_hits_fasta_selected_record_indices_override_selection() {
    let mut engine = GentleEngine::default();
    engine
        .upsert_rna_read_report(RnaReadInterpretationReport {
            schema: "gentle.rna_read_report.v1".to_string(),
            report_id: "rna_reads_selected_fasta".to_string(),
            seq_id: "seq_selected".to_string(),
            hits: vec![
                RnaReadInterpretationHit {
                    record_index: 0,
                    header_id: "read_alpha".to_string(),
                    sequence: "AACCGG".to_string(),
                    read_length_bp: 6,
                    best_mapping: Some(RnaReadMappingHit {
                        transcript_id: "tx_alpha".to_string(),
                        ..RnaReadMappingHit::default()
                    }),
                    ..RnaReadInterpretationHit::default()
                },
                RnaReadInterpretationHit {
                    record_index: 1,
                    header_id: "read_beta".to_string(),
                    sequence: "TTAACC".to_string(),
                    read_length_bp: 6,
                    ..RnaReadInterpretationHit::default()
                },
            ],
            ..RnaReadInterpretationReport::default()
        })
        .expect("upsert selected fasta report");

    let td = tempdir().expect("tempdir");
    let fasta_path = td.path().join("selected_reads.fa");
    let written = engine
        .export_rna_read_hits_fasta(
            "rna_reads_selected_fasta",
            fasta_path.to_str().expect("fasta path"),
            RnaReadHitSelection::Aligned,
            &[1],
            Some("filter=selected only | sort=score | search=<none>"),
        )
        .expect("export selected fasta");
    assert_eq!(written, 1);
    let fasta_text = fs::read_to_string(&fasta_path).expect("read selected fasta");
    assert!(fasta_text.contains(">read_beta "));
    assert!(fasta_text.contains("record_index=1"));
    assert!(fasta_text.contains("subset_spec=filter=selected only | sort=score | search=<none>"));
    assert!(!fasta_text.contains("read_alpha"));
}

#[test]
fn test_inspect_rna_read_alignments_reports_phase1_vs_phase2_effects_and_support_attribution() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let splicing = engine
        .build_splicing_expert_view(
            "seq_a",
            feature_id,
            SplicingScopePreset::AllOverlappingBothStrands,
        )
        .expect("splicing view");
    assert!(splicing.transcripts.len() >= 2);
    let tx_confirm = &splicing.transcripts[0];
    let tx_alt = &splicing.transcripts[1];
    let tx_confirm_len = tx_confirm
        .exons
        .iter()
        .map(|exon| {
            exon.end_1based
                .saturating_sub(exon.start_1based)
                .saturating_add(1)
        })
        .sum::<usize>();
    let tx_confirm_start = tx_confirm.exons.first().expect("first exon").start_1based;
    let tx_confirm_end = tx_confirm.exons.last().expect("last exon").end_1based;

    engine
        .upsert_rna_read_report(RnaReadInterpretationReport {
            schema: "gentle.rna_read_report.v1".to_string(),
            report_id: "rna_reads_effects".to_string(),
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            hits: vec![
                RnaReadInterpretationHit {
                    record_index: 0,
                    header_id: "confirmed".to_string(),
                    sequence: "ACGTACGTACGT".to_string(),
                    read_length_bp: 12,
                    passed_seed_filter: true,
                    seed_chain_transcript_id: tx_confirm.transcript_id.clone(),
                    exon_path_transcript_id: tx_confirm.transcript_id.clone(),
                    exon_path: "1:2:3".to_string(),
                    exon_transitions_confirmed: 2,
                    exon_transitions_total: 2,
                    best_mapping: Some(RnaReadMappingHit {
                        transcript_feature_id: tx_confirm.transcript_feature_id,
                        transcript_id: tx_confirm.transcript_id.clone(),
                        transcript_label: tx_confirm.label.clone(),
                        strand: tx_confirm.strand.clone(),
                        target_start_1based: tx_confirm_start,
                        target_end_1based: tx_confirm_end,
                        target_start_offset_0based: 0,
                        target_end_offset_0based_exclusive: tx_confirm_len,
                        identity_fraction: 0.97,
                        query_coverage_fraction: 1.0,
                        score: 220,
                        ..RnaReadMappingHit::default()
                    }),
                    ..RnaReadInterpretationHit::default()
                },
                RnaReadInterpretationHit {
                    record_index: 1,
                    header_id: "reassigned".to_string(),
                    sequence: "ACGTACGTACGT".to_string(),
                    read_length_bp: 12,
                    passed_seed_filter: true,
                    seed_chain_transcript_id: tx_alt.transcript_id.clone(),
                    exon_path_transcript_id: tx_alt.transcript_id.clone(),
                    exon_path: "1-2-3".to_string(),
                    exon_transitions_confirmed: 0,
                    exon_transitions_total: 2,
                    best_mapping: Some(RnaReadMappingHit {
                        transcript_feature_id: tx_confirm.transcript_feature_id,
                        transcript_id: tx_confirm.transcript_id.clone(),
                        transcript_label: tx_confirm.label.clone(),
                        strand: tx_confirm.strand.clone(),
                        target_start_1based: tx_confirm_start,
                        target_end_1based: tx_confirm_end,
                        target_start_offset_0based: 0,
                        target_end_offset_0based_exclusive: tx_confirm_len,
                        identity_fraction: 0.89,
                        query_coverage_fraction: 0.92,
                        score: 180,
                        ..RnaReadMappingHit::default()
                    }),
                    ..RnaReadInterpretationHit::default()
                },
                RnaReadInterpretationHit {
                    record_index: 2,
                    header_id: "no_phase1".to_string(),
                    sequence: "ACGTACGTACGT".to_string(),
                    read_length_bp: 12,
                    passed_seed_filter: true,
                    best_mapping: Some(RnaReadMappingHit {
                        transcript_feature_id: tx_confirm.transcript_feature_id,
                        transcript_id: tx_confirm.transcript_id.clone(),
                        transcript_label: tx_confirm.label.clone(),
                        strand: tx_confirm.strand.clone(),
                        target_start_1based: tx_confirm_start,
                        target_end_1based: tx_confirm_end,
                        target_start_offset_0based: 0,
                        target_end_offset_0based_exclusive: tx_confirm_len,
                        identity_fraction: 0.83,
                        query_coverage_fraction: 0.88,
                        score: 150,
                        ..RnaReadMappingHit::default()
                    }),
                    ..RnaReadInterpretationHit::default()
                },
            ],
            ..RnaReadInterpretationReport::default()
        })
        .expect("upsert RNA-read report");

    let inspection = engine
        .inspect_rna_read_alignments("rna_reads_effects", RnaReadHitSelection::All, 10)
        .expect("inspect read effects");
    assert_eq!(inspection.row_count, 3);
    let rows_by_header = inspection
        .rows
        .iter()
        .map(|row| (row.header_id.as_str(), row))
        .collect::<HashMap<_, _>>();

    let confirmed = rows_by_header.get("confirmed").expect("confirmed row");
    assert_eq!(
        confirmed.alignment_effect,
        RnaReadAlignmentEffect::ConfirmedAssignment
    );
    assert_eq!(
        confirmed.phase1_primary_transcript_id,
        tx_confirm.transcript_id
    );
    assert_eq!(confirmed.target_length_bp, tx_confirm_len);
    assert!((confirmed.target_coverage_fraction - 1.0).abs() < 1e-9);
    assert_eq!(confirmed.mapped_exon_support.len(), tx_confirm.exons.len());
    assert_eq!(
        confirmed
            .mapped_exon_support
            .iter()
            .map(|row| (row.start_1based, row.end_1based))
            .collect::<Vec<_>>(),
        tx_confirm
            .exons
            .iter()
            .map(|exon| (exon.start_1based, exon.end_1based))
            .collect::<Vec<_>>()
    );

    let reassigned = rows_by_header.get("reassigned").expect("reassigned row");
    assert_eq!(
        reassigned.alignment_effect,
        RnaReadAlignmentEffect::ReassignedTranscript
    );
    assert_eq!(
        reassigned.phase1_primary_transcript_id,
        tx_alt.transcript_id
    );
    assert_eq!(reassigned.target_length_bp, tx_confirm_len);
    assert!((reassigned.target_coverage_fraction - 1.0).abs() < 1e-9);

    let no_phase1 = rows_by_header.get("no_phase1").expect("no-phase1 row");
    assert_eq!(
        no_phase1.alignment_effect,
        RnaReadAlignmentEffect::AlignedWithoutPhase1Assignment
    );
    assert!(no_phase1.phase1_primary_transcript_id.is_empty());
    assert_eq!(no_phase1.target_length_bp, tx_confirm_len);
    assert!((no_phase1.target_coverage_fraction - 1.0).abs() < 1e-9);
}

#[test]
fn test_build_rna_read_alignment_display_reports_query_orientation_and_exact_midline() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let dna = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present");
    let feature_id = dna
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let splicing = engine
        .build_splicing_expert_view(
            "seq_a",
            feature_id,
            SplicingScopePreset::AllOverlappingBothStrands,
        )
        .expect("splicing view");
    let lane = splicing.transcripts.first().expect("first transcript lane");
    let template = GentleEngine::make_transcript_template(dna, lane, 3);
    let read = GentleEngine::reverse_complement_bytes(&template.sequence);
    let mapping = GentleEngine::align_read_to_template(
        &read,
        &template,
        &RnaReadAlignConfig {
            min_identity_fraction: 0.60,
            ..RnaReadAlignConfig::default()
        },
        3,
    )
    .expect("reverse-complement alignment");

    engine
        .upsert_rna_read_report(RnaReadInterpretationReport {
            schema: "gentle.rna_read_report.v1".to_string(),
            report_id: "rna_alignment_display".to_string(),
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            hits: vec![RnaReadInterpretationHit {
                record_index: 0,
                header_id: "rc_read".to_string(),
                sequence: String::from_utf8(read.clone()).expect("ASCII DNA"),
                read_length_bp: read.len(),
                passed_seed_filter: true,
                best_mapping: Some(mapping),
                ..RnaReadInterpretationHit::default()
            }],
            ..RnaReadInterpretationReport::default()
        })
        .expect("upsert RNA-read report");

    let display = engine
        .build_rna_read_alignment_display("rna_alignment_display", 0)
        .expect("alignment display");

    assert!(display.query_reverse_complemented);
    assert_eq!(display.alignment_mode, RnaReadAlignmentMode::Semiglobal);
    assert_eq!(display.matches, template.sequence.len());
    assert_eq!(display.mismatches, 0);
    assert_eq!(display.insertions, 0);
    assert_eq!(display.deletions, 0);
    assert_eq!(display.aligned_query.len(), display.aligned_midline.len());
    assert_eq!(display.aligned_query.len(), display.aligned_target.len());
    assert!(
        display.aligned_midline.chars().all(|ch| ch == '|'),
        "midline should show exact matches only: {}",
        display.aligned_midline
    );
}

#[test]
fn test_collect_mapped_isoform_support_rows_prefers_best_supported_mapping() {
    let hits = vec![
        RnaReadInterpretationHit {
            record_index: 0,
            header_id: "tx_a_1".to_string(),
            msa_eligible: true,
            secondary_mappings: vec![RnaReadMappingHit {
                transcript_id: "tx_b".to_string(),
                transcript_label: "TX_B".to_string(),
                strand: "-".to_string(),
                score: 80,
                identity_fraction: 0.88,
                query_coverage_fraction: 0.70,
                ..RnaReadMappingHit::default()
            }],
            best_mapping: Some(RnaReadMappingHit {
                transcript_feature_id: 11,
                transcript_id: "tx_a".to_string(),
                transcript_label: "TX_A".to_string(),
                strand: "+".to_string(),
                score: 120,
                identity_fraction: 0.92,
                query_coverage_fraction: 0.75,
                ..RnaReadMappingHit::default()
            }),
            ..RnaReadInterpretationHit::default()
        },
        RnaReadInterpretationHit {
            record_index: 1,
            header_id: "tx_a_2".to_string(),
            msa_eligible: false,
            best_mapping: Some(RnaReadMappingHit {
                transcript_feature_id: 11,
                transcript_id: "tx_a".to_string(),
                transcript_label: "TX_A".to_string(),
                strand: "+".to_string(),
                score: 105,
                identity_fraction: 0.86,
                query_coverage_fraction: 0.65,
                ..RnaReadMappingHit::default()
            }),
            ..RnaReadInterpretationHit::default()
        },
        RnaReadInterpretationHit {
            record_index: 2,
            header_id: "tx_b_1".to_string(),
            msa_eligible: true,
            best_mapping: Some(RnaReadMappingHit {
                transcript_feature_id: 12,
                transcript_id: "tx_b".to_string(),
                transcript_label: "TX_B".to_string(),
                strand: "-".to_string(),
                score: 140,
                identity_fraction: 0.98,
                query_coverage_fraction: 0.96,
                ..RnaReadMappingHit::default()
            }),
            ..RnaReadInterpretationHit::default()
        },
    ];

    let rows = GentleEngine::collect_mapped_isoform_support_rows(&hits);
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].transcript_id, "tx_a");
    assert_eq!(rows[0].aligned_read_count, 2);
    assert_eq!(rows[0].msa_eligible_read_count, 1);
    assert!((rows[0].mean_identity_fraction - 0.89).abs() < 1e-9);
    assert!((rows[0].mean_query_coverage_fraction - 0.70).abs() < 1e-9);
    assert_eq!(rows[0].best_alignment_score, 120);
    assert_eq!(rows[0].secondary_mapping_total, 1);

    assert_eq!(rows[1].transcript_id, "tx_b");
    assert_eq!(rows[1].aligned_read_count, 1);
    assert_eq!(rows[1].msa_eligible_read_count, 1);
    assert_eq!(rows[1].best_alignment_score, 140);
}

#[test]
fn test_interpret_rna_reads_populates_exon_and_junction_support_frequencies() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads.fa");
    fs::write(&input_path, format!(">read_1\n{read_sequence}\n")).expect("write reads");
    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: Default::default(),
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter: RnaReadSeedFilterConfig::default(),
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_support".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret reads");
    let report = engine
        .get_rna_read_report("rna_reads_support")
        .expect("report");
    assert!(!report.exon_support_frequencies.is_empty());
    assert!(
        report
            .exon_support_frequencies
            .iter()
            .all(|row| row.support_read_count == 0)
    );
    assert!(!report.junction_support_frequencies.is_empty());
    assert_eq!(report.read_count_aligned, 0);
    assert!(
        report
            .warnings
            .iter()
            .any(|warning| { warning.contains("phase-1 profile runs seed filtering only") })
    );
}

#[test]
fn test_export_rna_read_sample_sheet_writes_frequency_columns() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads.fa");
    fs::write(&input_path, format!(">read_1\n{read_sequence}\n")).expect("write reads");
    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: Default::default(),
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter: RnaReadSeedFilterConfig::default(),
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_sheet".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret reads");
    let sheet_path = td.path().join("sample_sheet.tsv");
    let export = engine
        .export_rna_read_sample_sheet(
            sheet_path.to_str().expect("sheet path"),
            Some("seq_a"),
            &[],
            false,
        )
        .expect("export sample sheet");
    assert_eq!(export.report_count, 1);
    let text = fs::read_to_string(&sheet_path).expect("read sample sheet");
    assert!(text.contains("exon_support_frequencies_json"));
    assert!(text.contains("junction_support_frequencies_json"));
    assert!(text.contains("origin_mode"));
    assert!(text.contains("target_gene_ids_json"));
    assert!(text.contains("origin_class_counts_json"));
    assert!(text.contains("rna_reads_sheet"));
}

#[test]
fn test_export_rna_read_exon_paths_and_abundance_tsv() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads.fa");
    fs::write(&input_path, format!(">read_1\n{read_sequence}\n")).expect("write reads");
    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: Default::default(),
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter: RnaReadSeedFilterConfig::default(),
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_paths".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret reads");

    let path_tsv = td.path().join("paths.tsv");
    let path_export = engine
        .export_rna_read_exon_paths_tsv(
            "rna_reads_paths",
            path_tsv.to_str().expect("path tsv"),
            RnaReadHitSelection::All,
            &[],
            None,
        )
        .expect("export paths");
    assert_eq!(path_export.row_count, 1);
    let path_text = fs::read_to_string(path_tsv).expect("read path tsv");
    assert!(path_text.contains("# report_id=rna_reads_paths"));
    assert!(path_text.contains("# seed_filter: k=10 stride=1"));
    assert!(path_text.contains("adjacent windows overlap by 9 bp"));
    assert!(path_text.contains("exon_path"));
    assert!(path_text.contains("reverse_complement_applied"));

    let abundance_tsv = td.path().join("abundance.tsv");
    let abundance_export = engine
        .export_rna_read_exon_abundance_tsv(
            "rna_reads_paths",
            abundance_tsv.to_str().expect("abundance tsv"),
            RnaReadHitSelection::All,
            &[],
            None,
        )
        .expect("export abundance");
    assert_eq!(abundance_export.selected_read_count, 1);
    let abundance_text = fs::read_to_string(abundance_tsv).expect("read abundance tsv");
    assert!(abundance_text.contains("# report_id=rna_reads_paths"));
    assert!(abundance_text.contains("# seed_filter: k=10 stride=1"));
    assert!(abundance_text.contains("adjacent windows overlap by 9 bp"));
    assert!(abundance_text.contains("row_kind"));
}

#[test]
fn test_export_rna_read_exon_paths_and_abundance_selected_record_indices_override_selection() {
    let mut engine = GentleEngine::default();
    engine
        .upsert_rna_read_report(RnaReadInterpretationReport {
            schema: "gentle.rna_read_report.v1".to_string(),
            report_id: "rna_reads_subset_tables".to_string(),
            seq_id: "seq_subset".to_string(),
            hits: vec![
                RnaReadInterpretationHit {
                    record_index: 0,
                    header_id: "read_alpha".to_string(),
                    sequence: "AACCGG".to_string(),
                    read_length_bp: 6,
                    exon_path_transcript_id: "tx_alpha".to_string(),
                    exon_path: "1:2".to_string(),
                    ..RnaReadInterpretationHit::default()
                },
                RnaReadInterpretationHit {
                    record_index: 1,
                    header_id: "read_beta".to_string(),
                    sequence: "TTAACC".to_string(),
                    read_length_bp: 6,
                    exon_path_transcript_id: "tx_beta".to_string(),
                    exon_path: "2:3".to_string(),
                    ..RnaReadInterpretationHit::default()
                },
            ],
            ..RnaReadInterpretationReport::default()
        })
        .expect("upsert subset tables report");

    let td = tempdir().expect("tempdir");
    let paths_tsv = td.path().join("subset_paths.tsv");
    let path_export = engine
        .export_rna_read_exon_paths_tsv(
            "rna_reads_subset_tables",
            paths_tsv.to_str().expect("paths tsv"),
            RnaReadHitSelection::All,
            &[1],
            Some("filter=disagreement only | sort=score | search=tp53"),
        )
        .expect("export selected paths");
    assert_eq!(path_export.row_count, 1);
    let paths_text = fs::read_to_string(&paths_tsv).expect("read selected paths");
    assert!(paths_text.contains("selected_record_indices=1"));
    assert!(paths_text.contains("subset_spec=filter=disagreement only | sort=score | search=tp53"));
    assert!(paths_text.contains("read_beta"));
    assert!(!paths_text.contains("read_alpha"));

    let abundance_tsv = td.path().join("subset_abundance.tsv");
    let abundance_export = engine
        .export_rna_read_exon_abundance_tsv(
            "rna_reads_subset_tables",
            abundance_tsv.to_str().expect("abundance tsv"),
            RnaReadHitSelection::All,
            &[1],
            Some("filter=disagreement only | sort=score | search=tp53"),
        )
        .expect("export selected abundance");
    assert_eq!(abundance_export.selected_read_count, 1);
    let abundance_text = fs::read_to_string(&abundance_tsv).expect("read selected abundance");
    assert!(abundance_text.contains("selected_record_indices=1"));
    assert!(
        abundance_text.contains("subset_spec=filter=disagreement only | sort=score | search=tp53")
    );
    assert!(abundance_text.contains("\ttransition\t\t2\t3\t"));
    assert!(!abundance_text.contains("\ttransition\t\t1\t2\t"));
}

#[test]
fn test_export_rna_read_score_density_svg() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads.fa");
    fs::write(&input_path, format!(">read_1\n{read_sequence}\n")).expect("write reads");
    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: Default::default(),
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter: RnaReadSeedFilterConfig::default(),
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_density".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret reads");

    let log_svg = td.path().join("score_density_log.svg");
    let log_export = engine
        .export_rna_read_score_density_svg(
            "rna_reads_density",
            log_svg.to_str().expect("log svg path"),
            RnaReadScoreDensityScale::Log,
            RnaReadScoreDensityVariant::AllScored,
            None,
        )
        .expect("export score-density log svg");
    assert_eq!(log_export.scale, RnaReadScoreDensityScale::Log);
    assert!(!log_export.derived_from_report_hits_only);
    let log_text = fs::read_to_string(&log_svg).expect("read log svg");
    assert!(log_text.contains("log(1+count)"));
    assert!(log_text.contains("report=rna_reads_density"));
    assert!(log_text.contains("profile=nanopore_cdna_v1"));
    assert!(log_text.contains("seed_filter: k=10 stride=1"));
    assert!(log_text.contains("adjacent windows overlap by 9 bp"));
    assert!(log_text.contains("score-density bins stored in report"));

    let mut report = engine
        .get_rna_read_report("rna_reads_density")
        .expect("report");
    report.score_density_bins.clear();
    engine
        .upsert_rna_read_report(report)
        .expect("upsert report without bins");

    let linear_svg = td.path().join("score_density_linear.svg");
    let linear_export = engine
        .export_rna_read_score_density_svg(
            "rna_reads_density",
            linear_svg.to_str().expect("linear svg path"),
            RnaReadScoreDensityScale::Linear,
            RnaReadScoreDensityVariant::AllScored,
            None,
        )
        .expect("export score-density linear svg");
    assert_eq!(linear_export.scale, RnaReadScoreDensityScale::Linear);
    assert!(linear_export.derived_from_report_hits_only);
    let linear_text = fs::read_to_string(&linear_svg).expect("read linear svg");
    assert!(linear_text.contains("linear count"));
    assert!(linear_text.contains("score-density bins derived from retained hits"));
}

#[test]
fn test_export_rna_read_score_density_svg_includes_compact_bar_labels() {
    let mut engine = GentleEngine::default();
    engine
        .upsert_rna_read_report(RnaReadInterpretationReport {
            schema: "gentle.rna_read_report.v1".to_string(),
            report_id: "rna_reads_density_labels".to_string(),
            seq_id: "seq_density".to_string(),
            score_density_bins: vec![0, 1_234, 15, 0],
            ..RnaReadInterpretationReport::default()
        })
        .expect("upsert density-label report");

    let td = tempdir().expect("tempdir");
    let svg_path = td.path().join("score_density_labels.svg");
    engine
        .export_rna_read_score_density_svg(
            "rna_reads_density_labels",
            svg_path.to_str().expect("svg path"),
            RnaReadScoreDensityScale::Log,
            RnaReadScoreDensityVariant::AllScored,
            None,
        )
        .expect("export score-density svg with labels");
    let svg_text = fs::read_to_string(&svg_path).expect("read svg");
    assert!(svg_text.contains("1.2k"));
    assert!(svg_text.contains(">15<"));
}

#[test]
fn test_export_rna_read_score_density_svg_composite_variant_uses_seed_pass_bins() {
    let mut engine = GentleEngine::default();
    engine
        .upsert_rna_read_report(RnaReadInterpretationReport {
            schema: "gentle.rna_read_report.v1".to_string(),
            report_id: "rna_reads_density_composite".to_string(),
            seq_id: "seq_density".to_string(),
            score_density_bins: vec![9, 0, 0, 0],
            seed_pass_score_density_bins: vec![2, 0, 0, 0],
            ..RnaReadInterpretationReport::default()
        })
        .expect("upsert composite density report");

    let td = tempdir().expect("tempdir");
    let svg_path = td.path().join("score_density_composite.svg");
    let export = engine
        .export_rna_read_score_density_svg(
            "rna_reads_density_composite",
            svg_path.to_str().expect("svg path"),
            RnaReadScoreDensityScale::Linear,
            RnaReadScoreDensityVariant::CompositeSeedGate,
            None,
        )
        .expect("export composite score-density svg");
    assert_eq!(
        export.variant,
        RnaReadScoreDensityVariant::CompositeSeedGate
    );
    assert_eq!(export.total_scored_reads, 2);
    let svg_text = fs::read_to_string(&svg_path).expect("read svg");
    assert!(svg_text.contains("composite seed-gate reads"));
    assert!(svg_text.contains("variant=composite_seed_gate"));
}

#[test]
fn test_export_rna_read_score_density_svg_retained_replay_variant_uses_override() {
    let mut engine = GentleEngine::default();
    engine
        .upsert_rna_read_report(RnaReadInterpretationReport {
            schema: "gentle.rna_read_report.v1".to_string(),
            report_id: "rna_reads_density_replay".to_string(),
            seq_id: "seq_density".to_string(),
            hits: vec![
                RnaReadInterpretationHit {
                    record_index: 0,
                    seed_hit_fraction: 0.84,
                    weighted_seed_hit_fraction: 0.84,
                    tested_kmers: 12,
                    matched_kmers: 12,
                    seed_chain_support_fraction: 1.0,
                    seed_median_transcript_gap: 1.0,
                    seed_transcript_gap_count: 1,
                    exon_transitions_confirmed: 1,
                    exon_transitions_total: 1,
                    ..RnaReadInterpretationHit::default()
                },
                RnaReadInterpretationHit {
                    record_index: 1,
                    seed_hit_fraction: 0.87,
                    weighted_seed_hit_fraction: 0.87,
                    tested_kmers: 12,
                    matched_kmers: 12,
                    seed_chain_support_fraction: 0.05,
                    seed_median_transcript_gap: 1.0,
                    seed_transcript_gap_count: 1,
                    exon_transitions_confirmed: 1,
                    exon_transitions_total: 1,
                    ..RnaReadInterpretationHit::default()
                },
            ],
            ..RnaReadInterpretationReport::default()
        })
        .expect("upsert replay density report");

    let replay_filter = RnaReadSeedFilterConfig {
        min_seed_hit_fraction: 0.80,
        min_weighted_seed_hit_fraction: 0.80,
        min_unique_matched_kmers: 10,
        min_chain_consistency_fraction: 0.50,
        max_median_transcript_gap: 4.0,
        min_confirmed_exon_transitions: 1,
        min_transition_support_fraction: 0.50,
        ..RnaReadSeedFilterConfig::default()
    };

    let td = tempdir().expect("tempdir");
    let svg_path = td.path().join("score_density_replay.svg");
    let export = engine
        .export_rna_read_score_density_svg(
            "rna_reads_density_replay",
            svg_path.to_str().expect("svg path"),
            RnaReadScoreDensityScale::Linear,
            RnaReadScoreDensityVariant::RetainedReplayCurrentControls,
            Some(&replay_filter),
        )
        .expect("export retained replay score-density svg");
    assert_eq!(
        export.variant,
        RnaReadScoreDensityVariant::RetainedReplayCurrentControls
    );
    assert_eq!(export.total_scored_reads, 1);
    let svg_text = fs::read_to_string(&svg_path).expect("read svg");
    assert!(svg_text.contains("retained replay under current controls"));
    assert!(svg_text.contains("variant=retained_replay_current_controls"));
    assert!(svg_text.contains("replay_seed_filter=seed_filter: k=10 stride=1"));
}

#[test]
fn test_materialize_rna_read_hit_sequences_creates_selected_sequences() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_origin".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("origin sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .upsert_rna_read_report(RnaReadInterpretationReport {
            schema: "gentle.rna_read_report.v1".to_string(),
            report_id: "rna_reads_materialize".to_string(),
            seq_id: "seq_origin".to_string(),
            hits: vec![
                RnaReadInterpretationHit {
                    record_index: 0,
                    header_id: "read_alpha".to_string(),
                    sequence: "AACCGG".to_string(),
                    read_length_bp: 6,
                    ..RnaReadInterpretationHit::default()
                },
                RnaReadInterpretationHit {
                    record_index: 1,
                    header_id: "read_beta".to_string(),
                    sequence: "TTAACC".to_string(),
                    read_length_bp: 6,
                    ..RnaReadInterpretationHit::default()
                },
            ],
            ..RnaReadInterpretationReport::default()
        })
        .expect("upsert materialize report");

    let result = engine
        .apply(Operation::MaterializeRnaReadHitSequences {
            report_id: "rna_reads_materialize".to_string(),
            selection: RnaReadHitSelection::All,
            selected_record_indices: vec![1],
            output_prefix: Some("tp73_outliers".to_string()),
        })
        .expect("materialize RNA-read hit sequence");
    assert_eq!(result.created_seq_ids.len(), 1);
    let created_seq_id = &result.created_seq_ids[0];
    assert!(created_seq_id.starts_with("tp73_outliers_r2_read_beta"));
    let created = engine
        .state()
        .sequences
        .get(created_seq_id)
        .expect("created read sequence present");
    assert_eq!(created.get_forward_string(), "TTAACC");
    assert!(result.messages.iter().any(|message| {
        message.contains("Materialized 1 RNA-read hit sequence(s)")
            && message.contains("rna_reads_materialize")
    }));
}

#[test]
fn test_interpret_rna_reads_retention_rescue_can_exceed_baseline_top_5000_in_memory() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let read_sequence = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_many.fa");
    let mut fasta = String::new();
    for idx in 0..5105usize {
        fasta.push_str(&format!(">read_{}\n{}\n", idx + 1, read_sequence));
    }
    fs::write(&input_path, fasta).expect("write reads");
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.min_seed_hit_fraction = 0.0;
    seed_filter.min_weighted_seed_hit_fraction = 0.0;
    seed_filter.min_unique_matched_kmers = 0;
    seed_filter.min_chain_consistency_fraction = 0.0;
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;

    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: Default::default(),
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_top5000".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret reads");
    let report = engine
        .get_rna_read_report("rna_reads_top5000")
        .expect("report");
    assert_eq!(report.read_count_total, 5105);
    assert!(report.hits.len() > 5000);
    assert_eq!(report.hits.len(), report.read_count_total);
    assert!(report.read_count_seed_passed <= report.read_count_total);
    assert_eq!(report.read_count_aligned, 0);
    assert!(
        report
            .warnings
            .iter()
            .any(|warning| { warning.contains("phase-1 profile runs seed filtering only") })
    );
    assert!(report.warnings.iter().any(|warning| {
        warning.contains("retention rescue kept") && warning.contains("baseline top-5000")
    }));
    assert_eq!(
        report.hits.iter().map(|hit| hit.record_index).max(),
        Some(report.read_count_total - 1)
    );
}

#[test]
fn test_interpret_rna_reads_report_mode_seed_passed_only_filters_hits() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let pass_read = {
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let fail_read = "NNNNNNNNNNNNNNNNNNNNNNNNNNNN".to_string();
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_mode.fa");
    fs::write(
        &input_path,
        format!(">read_pass\n{pass_read}\n>read_fail\n{fail_read}\n"),
    )
    .expect("write reads");

    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.min_seed_hit_fraction = 0.5;
    seed_filter.min_weighted_seed_hit_fraction = 0.1;
    seed_filter.min_unique_matched_kmers = 8;
    seed_filter.min_chain_consistency_fraction = 0.3;
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;
    engine
        .apply(Operation::InterpretRnaReads {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: input_path.display().to_string(),
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_seed_only".to_string()),
            report_mode: RnaReadReportMode::SeedPassedOnly,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        })
        .expect("interpret reads");
    let report = engine
        .get_rna_read_report("rna_reads_seed_only")
        .expect("report");
    assert_eq!(report.report_mode, RnaReadReportMode::SeedPassedOnly);
    assert_eq!(report.read_count_total, 2);
    assert_eq!(report.read_count_seed_passed, 1);
    assert_eq!(report.hits.len(), 1);
    assert!(report.hits.iter().all(|hit| hit.passed_seed_filter));
    assert_eq!(report.hits[0].header_id, "read_pass");
}

#[test]
fn test_interpret_rna_reads_checkpoint_resume_matches_uninterrupted_run() {
    let mut base_state = ProjectState::default();
    base_state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let feature_id = base_state
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let kmer_len = RnaReadSeedFilterConfig::default().kmer_len;
    let pass_read = {
        let engine = GentleEngine::from_state(base_state.clone());
        let splicing = engine
            .build_splicing_expert_view(
                "seq_a",
                feature_id,
                SplicingScopePreset::AllOverlappingBothStrands,
            )
            .expect("splicing view");
        let dna = engine
            .state()
            .sequences
            .get("seq_a")
            .expect("sequence for transcript template");
        let template =
            GentleEngine::make_transcript_template(dna, &splicing.transcripts[0], kmer_len);
        String::from_utf8(template.sequence).expect("template sequence utf-8")
    };
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("reads_checkpoint.fa");
    let mut fasta = String::new();
    for idx in 0..25usize {
        fasta.push_str(&format!(">read_{}\n{}\n", idx + 1, pass_read));
    }
    fs::write(&input_path, fasta).expect("write reads");
    let checkpoint_path = td.path().join("reads_checkpoint.json");

    let full_engine = GentleEngine::from_state(base_state.clone());
    let report_full = full_engine
        .compute_rna_read_report_with_progress(
            "seq_a",
            feature_id,
            RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path.to_str().expect("path"),
            RnaReadInputFormat::Fasta,
            SplicingScopePreset::AllOverlappingBothStrands,
            RnaReadOriginMode::SingleGene,
            &[],
            false,
            &RnaReadSeedFilterConfig::default(),
            &RnaReadAlignConfig::default(),
            Some("rna_reads_checkpoint"),
            &mut |_progress| true,
        )
        .expect("full report");

    let interrupted_engine = GentleEngine::from_state(base_state.clone());
    let mut continue_budget = 40usize;
    let interrupted = interrupted_engine
        .compute_rna_read_report_with_options_and_progress_and_cancel(
            "seq_a",
            feature_id,
            RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path.to_str().expect("path"),
            RnaReadInputFormat::Fasta,
            SplicingScopePreset::AllOverlappingBothStrands,
            RnaReadOriginMode::SingleGene,
            &[],
            false,
            &RnaReadSeedFilterConfig::default(),
            &RnaReadAlignConfig::default(),
            Some("rna_reads_checkpoint"),
            &RnaReadInterpretOptions {
                report_mode: RnaReadReportMode::Full,
                checkpoint_path: Some(checkpoint_path.display().to_string()),
                checkpoint_every_reads: 5,
                resume_from_checkpoint: false,
            },
            &mut |_progress| true,
            &mut || {
                if continue_budget == 0 {
                    false
                } else {
                    continue_budget = continue_budget.saturating_sub(1);
                    true
                }
            },
        );
    assert!(
        interrupted.is_err(),
        "interrupted run should return cancellation error"
    );
    assert!(
        checkpoint_path.exists(),
        "checkpoint file should be written"
    );

    let resumed_engine = GentleEngine::from_state(base_state);
    let report_resumed = resumed_engine
        .compute_rna_read_report_with_options_and_progress_and_cancel(
            "seq_a",
            feature_id,
            RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path.to_str().expect("path"),
            RnaReadInputFormat::Fasta,
            SplicingScopePreset::AllOverlappingBothStrands,
            RnaReadOriginMode::SingleGene,
            &[],
            false,
            &RnaReadSeedFilterConfig::default(),
            &RnaReadAlignConfig::default(),
            Some("rna_reads_checkpoint"),
            &RnaReadInterpretOptions {
                report_mode: RnaReadReportMode::Full,
                checkpoint_path: Some(checkpoint_path.display().to_string()),
                checkpoint_every_reads: 5,
                resume_from_checkpoint: true,
            },
            &mut |_progress| true,
            &mut || true,
        )
        .expect("resumed report");

    assert_eq!(
        report_resumed.read_count_total,
        report_full.read_count_total
    );
    assert_eq!(
        report_resumed.read_count_seed_passed,
        report_full.read_count_seed_passed
    );
    assert_eq!(
        report_resumed.read_count_aligned,
        report_full.read_count_aligned
    );
    assert_eq!(
        report_resumed.score_density_bins,
        report_full.score_density_bins
    );
    assert_eq!(
        report_resumed
            .hits
            .iter()
            .map(|hit| hit.record_index)
            .collect::<Vec<_>>(),
        report_full
            .hits
            .iter()
            .map(|hit| hit.record_index)
            .collect::<Vec<_>>()
    );
    assert!(report_resumed.resumed_from_checkpoint);
}

#[test]
fn test_export_rna_seed_hash_catalog_writes_hash_sequences_and_positions() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), splicing_test_sequence());
    let engine = GentleEngine::from_state(state);
    let feature_id = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present")
        .features()
        .iter()
        .position(|feature| feature.kind.to_string().eq_ignore_ascii_case("mRNA"))
        .expect("mRNA feature id");
    let td = tempdir().expect("tempdir");
    let output_path = td.path().join("seed_hash_catalog.tsv");
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.kmer_len = 9;
    let (rows, unique_hashes) = engine
        .export_rna_seed_hash_catalog(
            "seq_a",
            feature_id,
            SplicingScopePreset::AllOverlappingBothStrands,
            &seed_filter,
            output_path.to_str().expect("catalog path"),
        )
        .expect("export seed hash catalog");
    assert!(rows > 0);
    assert!(unique_hashes > 0);
    let text = fs::read_to_string(&output_path).expect("read seed hash catalog");
    assert!(text.contains("seed_bits"));
    assert!(text.contains("kmer_sequence"));
    assert!(text.contains("genomic_pos_1based"));
}

#[test]
fn test_collect_rna_seed_hash_catalog_rows_deduplicates_overlapping_templates() {
    let kmer_len = 3usize;
    let bits = GentleEngine::encode_kmer_bits(b"ACG").expect("seed bits");
    let template_a = SplicingTranscriptTemplate {
        transcript_feature_id: 10,
        transcript_id: "tx_a".to_string(),
        transcript_label: "tx_a".to_string(),
        strand: "+".to_string(),
        sequence: b"ACG".to_vec(),
        genomic_positions_1based: vec![100, 101, 102],
        kmer_positions: HashMap::from([(bits, vec![0])]),
    };
    let template_b = SplicingTranscriptTemplate {
        transcript_feature_id: 11,
        transcript_id: "tx_b".to_string(),
        transcript_label: "tx_b".to_string(),
        strand: "+".to_string(),
        sequence: b"ACG".to_vec(),
        genomic_positions_1based: vec![100, 101, 102],
        kmer_positions: HashMap::from([(bits, vec![0])]),
    };
    let rows =
        GentleEngine::collect_rna_seed_hash_catalog_rows(&[template_a, template_b], kmer_len);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].seed_bits, bits);
    assert_eq!(rows[0].kmer_sequence, "ACG");
    assert_eq!(rows[0].genomic_pos_1based, 100);
}

#[test]
fn test_build_rna_read_seed_histogram_index_deduplicates_shared_template_positions() {
    let bits = GentleEngine::encode_kmer_bits(b"ACG").expect("seed bits");
    let template_a = SplicingTranscriptTemplate {
        transcript_feature_id: 10,
        transcript_id: "tx_a".to_string(),
        transcript_label: "tx_a".to_string(),
        strand: "+".to_string(),
        sequence: b"ACG".to_vec(),
        genomic_positions_1based: vec![100, 101, 102],
        kmer_positions: HashMap::from([(bits, vec![0])]),
    };
    let template_b = SplicingTranscriptTemplate {
        transcript_feature_id: 11,
        transcript_id: "tx_b".to_string(),
        transcript_label: "tx_b".to_string(),
        strand: "+".to_string(),
        sequence: b"ACG".to_vec(),
        genomic_positions_1based: vec![100, 101, 102],
        kmer_positions: HashMap::from([(bits, vec![0])]),
    };
    let bins = GentleEngine::build_rna_read_seed_histogram_bins(300);
    let histogram =
        GentleEngine::build_rna_read_seed_histogram_index(&[template_a, template_b], 300, &bins);
    let weights = histogram.get(&bits).expect("weights for seed bit");
    assert_eq!(weights.len(), 1);
    assert!(!weights[0].strand_minus);
}

#[test]
fn test_candidate_metric_expression_fuzz_smoke_does_not_panic() {
    let metrics = HashMap::from([
        ("gc_fraction".to_string(), 0.5),
        ("at_fraction".to_string(), 0.5),
        ("length_bp".to_string(), 20.0),
        ("candidate_index".to_string(), 3.0),
        ("distance_to_feature_start_bp".to_string(), 7.0),
    ]);
    let atoms = [
        "gc_fraction",
        "at_fraction",
        "length_bp",
        "candidate_index",
        "distance_to_feature_start_bp",
        "1",
        "2.5",
    ];
    let ops = ["+", "-", "*", "/"];

    for i in 0..512usize {
        let left = atoms[i % atoms.len()];
        let mid = atoms[(i / 7) % atoms.len()];
        let right = atoms[(i / 29) % atoms.len()];
        let op_a = ops[i % ops.len()];
        let op_b = ops[(i / 11) % ops.len()];
        let expr = format!("({left} {op_a} ({mid} {op_b} {right}))");
        let run = std::panic::catch_unwind(|| {
            if let Ok(parsed) = GentleEngine::parse_metric_expression(&expr) {
                let _ = GentleEngine::evaluate_metric_expression(&parsed, &metrics);
            }
        });
        assert!(run.is_ok(), "expression caused panic: {expr}");
    }
}

#[test]
fn query_sequence_features_filters_by_kind_range_strand_and_label() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(150)).expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "SOURCE".into(),
        location: gb_io::seq::Location::simple_range(0, 600),
        qualifiers: vec![("label".into(), Some("source".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(20, 160),
        qualifiers: vec![("label".into(), Some("TP73".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "CDS".into(),
        location: gb_io::seq::Location::Complement(Box::new(gb_io::seq::Location::simple_range(
            40, 140,
        ))),
        qualifiers: vec![
            ("label".into(), Some("TP73_DBD".to_string())),
            (
                "product".into(),
                Some("TP73 DNA-binding region".to_string()),
            ),
        ],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "promoter".into(),
        location: gb_io::seq::Location::simple_range(180, 260),
        qualifiers: vec![("note".into(), Some("TP73-AS2 promoter".to_string()))],
    });

    let mut state = ProjectState::default();
    state.sequences.insert("tp73".to_string(), dna);
    let engine = GentleEngine::from_state(state);

    let query = SequenceFeatureQuery {
        seq_id: "tp73".to_string(),
        kind_in: vec!["CDS".to_string()],
        start_0based: Some(30),
        end_0based_exclusive: Some(180),
        range_relation: SequenceFeatureRangeRelation::Within,
        strand: SequenceFeatureStrandFilter::Reverse,
        label_contains: Some("tp73".to_string()),
        sort_by: SequenceFeatureSortBy::Start,
        ..SequenceFeatureQuery::default()
    };
    let result = engine
        .query_sequence_features(query)
        .expect("feature query should succeed");
    assert_eq!(result.schema, "gentle.sequence_feature_query_result.v1");
    assert_eq!(result.total_feature_count, 4);
    assert_eq!(result.matched_count, 1);
    assert_eq!(result.returned_count, 1);
    assert_eq!(result.rows[0].feature_id, 2);
    assert_eq!(result.rows[0].kind.to_ascii_uppercase(), "CDS");
    assert_eq!(result.rows[0].strand, "reverse");
}

#[test]
fn query_sequence_features_applies_qualifier_filters_and_pagination() {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(120)).expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(10, 70),
        qualifiers: vec![
            ("label".into(), Some("TP73".to_string())),
            ("gene".into(), Some("TP73".to_string())),
            ("note".into(), Some("tumor protein p73".to_string())),
        ],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "gene".into(),
        location: gb_io::seq::Location::simple_range(80, 130),
        qualifiers: vec![
            ("label".into(), Some("TP53".to_string())),
            ("gene".into(), Some("TP53".to_string())),
            ("note".into(), Some("tumor protein p53".to_string())),
        ],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "misc_feature".into(),
        location: gb_io::seq::Location::simple_range(140, 200),
        qualifiers: vec![
            ("label".into(), Some("TP73_promoter".to_string())),
            ("note".into(), Some("promoter-like region".to_string())),
        ],
    });

    let mut state = ProjectState::default();
    state.sequences.insert("tp73".to_string(), dna);
    let engine = GentleEngine::from_state(state);

    let query = SequenceFeatureQuery {
        seq_id: "tp73".to_string(),
        qualifier_filters: vec![
            SequenceFeatureQualifierFilter {
                key: "gene".to_string(),
                value_contains: Some("tp".to_string()),
                value_regex: Some("^TP7".to_string()),
                case_sensitive: false,
            },
            SequenceFeatureQualifierFilter {
                key: "note".to_string(),
                value_contains: Some("tumor".to_string()),
                value_regex: None,
                case_sensitive: false,
            },
        ],
        sort_by: SequenceFeatureSortBy::FeatureId,
        descending: false,
        offset: 0,
        limit: Some(1),
        include_qualifiers: true,
        ..SequenceFeatureQuery::default()
    };
    let result = engine
        .query_sequence_features(query)
        .expect("feature query should succeed");
    assert_eq!(result.matched_count, 1);
    assert_eq!(result.returned_count, 1);
    assert_eq!(result.rows[0].feature_id, 0);
    assert_eq!(
        result.rows[0]
            .qualifiers
            .get("gene")
            .expect("gene qualifier")
            .first()
            .expect("gene value"),
        "TP73"
    );
}
