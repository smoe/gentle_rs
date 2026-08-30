//! Dedicated Criterion coverage for deterministic primer-specificity interpretation.
//!
//! The generated rows are synthetic, complete-query `blastn-short` HSPs. Each
//! subject has one inward-facing forward/reverse pair, with the first product
//! declared intended and every remaining exact product classified off-target.

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use gentle::{
    about::GENTLE_SOURCE_REVISION,
    engine::{
        GentleEngine, PrimerSpecificityExpectedProduct, PrimerSpecificityFullAlignmentMode,
        PrimerSpecificityHandoff, PrimerSpecificityHandoffCommand, PrimerSpecificityInputPrimer,
        PrimerSpecificityIntendedTarget, PrimerSpecificityIntendedTargetModel,
        PrimerSpecificityPolicy, PrimerSpecificityPrimerRole, PrimerSpecificityReport,
        PrimerSpecificityReportDetailMode, PrimerSpecificitySubjectRange,
    },
};
use ring::digest::{SHA256, digest};
use std::{collections::BTreeMap, fs, hint::black_box};
use tempfile::TempDir;

const TOTAL_HSP_COUNTS: [usize; 3] = [100, 1_000, 6_600];
const FORWARD_SEQUENCE: &str = "ACGTACGTACGTACGTACGT";
const REVERSE_SEQUENCE: &str = "TGCATGCATGCATGCATGCA";
const TARGET_GENOME_ID: &str = "Synthetic specificity benchmark";

struct SpecificityFixture {
    total_hsp_count: usize,
    paired_subject_count: usize,
    fixture_sha256: String,
    handoff: PrimerSpecificityHandoff,
    outputs: BTreeMap<String, String>,
}

fn sha256_prefixed(bytes: &[u8]) -> String {
    let hex = digest(&SHA256, bytes)
        .as_ref()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    format!("sha256:{hex}")
}

fn benchmark_catalog() -> (TempDir, String) {
    let root = tempfile::tempdir().expect("create benchmark catalog directory");
    let fasta = root.path().join("synthetic.fa");
    let annotations = root.path().join("synthetic.gff3");
    let cache = root.path().join("cache");
    let catalog = root.path().join("genomes.json");
    fs::write(&fasta, ">bench_chr\nACGTACGTACGTACGTACGT\n")
        .expect("write synthetic benchmark FASTA");
    fs::write(&annotations, "##gff-version 3\n").expect("write synthetic benchmark annotations");
    fs::write(
        &catalog,
        serde_json::to_vec_pretty(&serde_json::json!({
            (TARGET_GENOME_ID): {
                "description": "Generated non-biological Criterion fixture",
                "sequence_local": fasta,
                "annotations_local": annotations,
                "cache_dir": cache,
            }
        }))
        .expect("serialize benchmark catalog"),
    )
    .expect("write synthetic benchmark catalog");
    (root, catalog.to_string_lossy().to_string())
}

fn primer(role: PrimerSpecificityPrimerRole, sequence: &str) -> PrimerSpecificityInputPrimer {
    PrimerSpecificityInputPrimer {
        role,
        full_sequence: sequence.to_string(),
        annealing_sequence: sequence.to_string(),
        annealing_length_bp: sequence.len(),
        non_annealing_5prime_tail: String::new(),
        non_annealing_5prime_tail_bp: 0,
    }
}

fn command(
    role: PrimerSpecificityPrimerRole,
    query_label: &str,
    output_name: &str,
    query_length_bp: usize,
    max_target_seqs: usize,
) -> PrimerSpecificityHandoffCommand {
    PrimerSpecificityHandoffCommand {
        command_id: format!("{}_benchmark_blast", role.as_str()),
        role,
        query_label: query_label.to_string(),
        query_length_bp,
        query_fasta_path: format!("benchmark://{query_label}.fa"),
        output_tsv_path: format!("benchmark://{output_name}"),
        program: "blastn".to_string(),
        args: vec![
            "-task".to_string(),
            "blastn-short".to_string(),
            "-max_target_seqs".to_string(),
            max_target_seqs.to_string(),
        ],
        command_line: "synthetic benchmark; no process is executed".to_string(),
        success_exit_codes: vec![0],
    }
}

fn aligned_hsp_rows(
    query_label: &str,
    sequence: &str,
    subject_count: usize,
    reverse: bool,
) -> String {
    let mut output = String::with_capacity(subject_count * 140);
    for subject_index in 0..subject_count {
        let subject_id = format!("bench_tx_{subject_index:06}");
        let (subject_start, subject_end) = if reverse { (240, 221) } else { (100, 119) };
        output.push_str(&format!(
            "{query_label}\t{subject_id}\t100\t20\t0\t0\t1\t20\t{subject_start}\t{subject_end}\t1e-20\t80\t100\t20\t{sequence}\t{sequence}\n"
        ));
    }
    output
}

fn build_fixture(total_hsp_count: usize, catalog_path: &str) -> SpecificityFixture {
    assert_eq!(total_hsp_count % 2, 0, "fixture needs paired HSP rows");
    let paired_subject_count = total_hsp_count / 2;
    let forward_command = command(
        PrimerSpecificityPrimerRole::Forward,
        "forward_annealing_segment",
        "forward.tsv",
        FORWARD_SEQUENCE.len(),
        paired_subject_count,
    );
    let reverse_command = command(
        PrimerSpecificityPrimerRole::Reverse,
        "reverse_annealing_segment",
        "reverse.tsv",
        REVERSE_SEQUENCE.len(),
        paired_subject_count,
    );
    let forward_output = aligned_hsp_rows(
        &forward_command.query_label,
        FORWARD_SEQUENCE,
        paired_subject_count,
        false,
    );
    let reverse_output = aligned_hsp_rows(
        &reverse_command.query_label,
        REVERSE_SEQUENCE,
        paired_subject_count,
        true,
    );
    let mut fixture_bytes = Vec::with_capacity(forward_output.len() + reverse_output.len() + 1);
    fixture_bytes.extend_from_slice(forward_output.as_bytes());
    fixture_bytes.push(0);
    fixture_bytes.extend_from_slice(reverse_output.as_bytes());
    let fixture_sha256 = sha256_prefixed(&fixture_bytes);

    let intended_range = PrimerSpecificitySubjectRange {
        start_1based: 100,
        end_1based: 240,
    };
    let mut policy = PrimerSpecificityPolicy {
        report_detail_mode: PrimerSpecificityReportDetailMode::Compact,
        max_hits_per_primer: paired_subject_count,
        ..PrimerSpecificityPolicy::default()
    };
    policy.full_alignment.mode = PrimerSpecificityFullAlignmentMode::BestEffort;
    let handoff = PrimerSpecificityHandoff {
        schema: "gentle.primer_specificity_handoff.v1".to_string(),
        handoff_id: format!("criterion_specificity_{total_hsp_count}"),
        requested_target_genome_id: TARGET_GENOME_ID.to_string(),
        resolved_target_genome_id: TARGET_GENOME_ID.to_string(),
        catalog_path: Some(catalog_path.to_string()),
        policy,
        primers: vec![
            primer(PrimerSpecificityPrimerRole::Forward, FORWARD_SEQUENCE),
            primer(PrimerSpecificityPrimerRole::Reverse, REVERSE_SEQUENCE),
        ],
        blast_db_prefix: "benchmark://not-prepared".to_string(),
        intended_target: PrimerSpecificityIntendedTarget {
            model: PrimerSpecificityIntendedTargetModel::GenomicInterval,
            subject_id: Some("bench_tx_000000".to_string()),
            forward_binding_ranges: vec![PrimerSpecificitySubjectRange {
                start_1based: 100,
                end_1based: 119,
            }],
            reverse_binding_ranges: vec![PrimerSpecificitySubjectRange {
                start_1based: 221,
                end_1based: 240,
            }],
            expected_product_range: Some(intended_range.clone()),
            genomic_target_geometry_known: true,
            contiguous_genomic_product_expected: true,
            expected_products: vec![PrimerSpecificityExpectedProduct {
                target_space: "genomic_dna".to_string(),
                subject_id: "bench_tx_000000".to_string(),
                expected_product_range: Some(intended_range),
                source_transcript_id: Some("bench_tx_000000".to_string()),
            }],
            source: "deterministic Criterion fixture".to_string(),
            warnings: vec![],
        },
        commands: vec![forward_command.clone(), reverse_command.clone()],
        completion_policy: "all_commands_success".to_string(),
        import_command: vec![],
        import_command_line: String::new(),
        ..PrimerSpecificityHandoff::default()
    };
    let outputs = BTreeMap::from([
        (forward_command.command_id, forward_output),
        (reverse_command.command_id, reverse_output),
    ]);
    SpecificityFixture {
        total_hsp_count,
        paired_subject_count,
        fixture_sha256,
        handoff,
        outputs,
    }
}

fn assert_fixture_report(fixture: &SpecificityFixture, report: &PrimerSpecificityReport) {
    assert_eq!(
        report.compaction.raw_forward_hit_count + report.compaction.raw_reverse_hit_count,
        fixture.total_hsp_count
    );
    assert_eq!(
        report.compaction.raw_amplicon_count,
        fixture.paired_subject_count
    );
    assert_eq!(
        report.compaction.pairing_candidate_comparison_count,
        fixture.paired_subject_count
    );
    assert_eq!(report.summary.intended_amplicon_count, 1);
    assert_eq!(
        report.summary.failing_unintended_amplicon_count,
        fixture.paired_subject_count - 1
    );
}

fn fixture_id(fixture: &SpecificityFixture) -> String {
    let digest_token = fixture
        .fixture_sha256
        .strip_prefix("sha256:")
        .unwrap_or(&fixture.fixture_sha256);
    format!(
        "{}_{}",
        fixture.total_hsp_count,
        &digest_token[..12.min(digest_token.len())]
    )
}

fn benchmark_specificity_finalization(c: &mut Criterion) {
    let (_catalog_root, catalog_path) = benchmark_catalog();
    let fixtures = TOTAL_HSP_COUNTS
        .into_iter()
        .map(|count| build_fixture(count, &catalog_path))
        .collect::<Vec<_>>();
    let engine = GentleEngine::default();
    eprintln!("GENtle source revision: {GENTLE_SOURCE_REVISION}");

    let mut parser_group = c.benchmark_group("specificity_hsp_parsing");
    for fixture in &fixtures {
        let parsed_count = fixture
            .outputs
            .values()
            .map(|output| GentleEngine::benchmark_parse_primer_specificity_blast_output(output).0)
            .sum::<usize>();
        assert_eq!(parsed_count, fixture.total_hsp_count);
        parser_group.throughput(Throughput::Elements(fixture.total_hsp_count as u64));
        parser_group.bench_with_input(
            BenchmarkId::from_parameter(fixture_id(fixture)),
            fixture,
            |b, fixture| {
                b.iter(|| {
                    let parsed = fixture
                        .outputs
                        .values()
                        .map(|output| {
                            GentleEngine::benchmark_parse_primer_specificity_blast_output(
                                black_box(output),
                            )
                        })
                        .collect::<Vec<_>>();
                    black_box(parsed)
                });
            },
        );
    }
    parser_group.finish();

    let mut interpretation_group = c.benchmark_group("specificity_finalization_in_memory");
    for fixture in &fixtures {
        let report = engine
            .benchmark_primer_specificity_report_from_handoff_outputs(
                &fixture.handoff,
                &fixture.outputs,
            )
            .expect("interpret deterministic specificity fixture");
        assert_fixture_report(fixture, &report);
        interpretation_group.throughput(Throughput::Elements(fixture.total_hsp_count as u64));
        interpretation_group.bench_with_input(
            BenchmarkId::from_parameter(fixture_id(fixture)),
            fixture,
            |b, fixture| {
                b.iter(|| {
                    let report = engine
                        .benchmark_primer_specificity_report_from_handoff_outputs(
                            black_box(&fixture.handoff),
                            black_box(&fixture.outputs),
                        )
                        .expect("interpret deterministic specificity fixture");
                    black_box(report)
                });
            },
        );
    }
    interpretation_group.finish();
}

criterion_group!(benches, benchmark_specificity_finalization);
criterion_main!(benches);
