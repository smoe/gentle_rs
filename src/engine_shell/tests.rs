//! Shared-shell regression tests for parser and executor parity.
//!
//! These tests exercise the public shell surface across parser helpers,
//! execution routing, and engine integration instead of treating each helper as
//! an isolated API.
//!
//! Start here when you need examples of:
//! - shell-token parsing expectations
//! - `execute_shell_command` behavior over real engine state
//! - CLI-like round trips that should stay in parity with GUI shell mode

use super::*;
use crate::dna_sequence::DNAsequence;
use crate::engine::{
    AdapterCaptureProtectionMode, AdapterCaptureStyle, AdapterRestrictionCapturePlan, Arrangement,
    ArrangementMode, AttractPwmMappingPolicy, AttractSplicingEvidenceSettings,
    BIGWIG_TO_BEDGRAPH_ENV_BIN, ConstructObjective, ConstructRole, Container, ContainerKind,
    CutRunAlignConfig, CutRunCoverageKind, CutRunInputFormat, CutRunReadLayout,
    CutRunSeedFilterConfig, EditableStatus, PrimerDesignProgress, PromoterTfbsGeneQuery,
    ProteinExternalOpinionSource, ProteinFeatureFilter, Rack, RackAuthoringTemplate,
    RackCarrierLabelPreset, RackFillDirection, RackLabelSheetPreset, RackOccupant,
    RackPhysicalTemplateKind, RackPlacementEntry, RackProfileKind, RackProfileSnapshot,
    RestrictionCloningPcrHandoffMode, RnaReadAlignConfig, RnaReadInterpretationHit,
    RnaReadInterpretationReport, RnaReadMappingHit, RnaReadOriginClass, SequenceScanTarget,
    TfThresholdOverride, TfbsScoreTrackCorrelationSignalSource, TfbsScoreTrackValueKind,
    TfbsTrackSimilarityRankingMetric,
};
use crate::ensembl_gene::{EnsemblGeneEntry, EnsemblGeneTranscriptSummary};
use crate::ensembl_protein::{
    EnsemblProteinEntry, EnsemblProteinFeature, EnsemblTranscriptExon, EnsemblTranscriptTranslation,
};
use crate::test_support::{
    decision_trace_fixture_state, decision_trace_with_construct_reasoning_fixture_state,
    write_demo_pool_json, write_demo_workflow_json, write_demo_workflow_with_shebang,
};
use gb_io::seq::{Feature, Location};
use std::env;
use std::fs;
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::path::Path;
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock};
use tempfile::tempdir;

static JASPAR_RELOAD_TEST_COUNTER: AtomicUsize = AtomicUsize::new(1);
static ATTRACT_RELOAD_TEST_COUNTER: AtomicUsize = AtomicUsize::new(1);

struct EnvVarGuard {
    key: String,
    original: Option<String>,
}

impl EnvVarGuard {
    fn set(key: &str, value: &str) -> Self {
        let original = env::var(key).ok();
        unsafe {
            env::set_var(key, value);
        }
        Self {
            key: key.to_string(),
            original,
        }
    }
}

impl Drop for EnvVarGuard {
    fn drop(&mut self) {
        match &self.original {
            Some(value) => unsafe {
                env::set_var(&self.key, value);
            },
            None => unsafe {
                env::remove_var(&self.key);
            },
        }
    }
}

fn cache_dir_env_lock() -> &'static Mutex<()> {
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(()))
}

fn cutrun_test_env_lock() -> &'static Mutex<()> {
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(()))
}

fn attract_test_lock() -> &'static Mutex<()> {
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(()))
}

fn jaspar_test_lock() -> &'static Mutex<()> {
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(()))
}

struct JasparReloadResetGuard;

impl Drop for JasparReloadResetGuard {
    fn drop(&mut self) {
        crate::tf_motifs::reload();
    }
}

fn lock_jaspar_tests() -> std::sync::MutexGuard<'static, ()> {
    jaspar_test_lock()
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

fn with_blast_async_test_overrides<R>(
    max_concurrent: usize,
    worker_delay_ms: u64,
    f: impl FnOnce() -> R,
) -> R {
    let _guard = BLAST_ASYNC_TEST_MUTEX
        .lock()
        .expect("blast async test mutex lock");
    let previous_max =
        BLAST_ASYNC_MAX_CONCURRENT_TEST_OVERRIDE.swap(max_concurrent, Ordering::SeqCst);
    let previous_delay =
        BLAST_ASYNC_WORKER_DELAY_MS_TEST_OVERRIDE.swap(worker_delay_ms, Ordering::SeqCst);
    clear_blast_async_jobs_for_test();
    let result = f();
    clear_blast_async_jobs_for_test();
    BLAST_ASYNC_WORKER_DELAY_MS_TEST_OVERRIDE.store(previous_delay, Ordering::SeqCst);
    BLAST_ASYNC_MAX_CONCURRENT_TEST_OVERRIDE.store(previous_max, Ordering::SeqCst);
    result
}

fn resource_fixture_path(name: &str) -> String {
    format!(
        "{}/test_files/fixtures/resources/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

fn write_synthetic_attract_zip(dir: &Path) -> std::path::PathBuf {
    let unique = ATTRACT_RELOAD_TEST_COUNTER.fetch_add(1, Ordering::Relaxed);
    let case_dir = dir.join(format!("attract_case_{unique}"));
    fs::create_dir_all(&case_dir).expect("create attract case dir");
    let db_path = case_dir.join("ATtRACT_db.txt");
    let pwm_path = case_dir.join("pwm.txt");
    fs::write(
        &db_path,
        "Gene_name\tGene_id\tOrganism\tMatrix_id\tMotif\tLen\tExperiment\tDomain\tPubmed\tQuality_score\nSRSF1\tENSG00000136450\tHomo sapiens\tM001\tGAAGAA\t6\tSELEX\tRRM\t12345\t4.2\nPTBP1\tENSG00000011304\tHomo sapiens\tM002\tUCUU\t4\tCLIP\tRRM\t23456\t3.1\n",
    )
    .expect("write attract db");
    fs::write(
        &pwm_path,
        "M001\nA\t5\t0\t0\t0\t5\t5\nC\t0\t5\t0\t0\t0\t0\nG\t0\t0\t5\t0\t0\t0\nT\t0\t0\t0\t5\t0\t0\n",
    )
    .expect("write pwm placeholder");
    let archive_path = dir.join(format!("attract_{unique}.zip"));
    let output = Command::new("zip")
        .arg("-jq")
        .arg(&archive_path)
        .arg(&db_path)
        .arg(&pwm_path)
        .output()
        .expect("run zip");
    assert!(
        output.status.success(),
        "zip command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    archive_path
}

fn primer3_fixture_path(name: &str) -> String {
    format!(
        "{}/test_files/fixtures/primer3/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

fn sequencing_confirmation_fixture_path(name: &str) -> String {
    format!(
        "{}/test_files/fixtures/sequencing_confirmation/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

fn shell_file_url(path: &Path) -> String {
    format!("file://{}", path.display())
}

fn write_cutrun_shell_reference_catalog(root: &Path, genome_id: &str) -> String {
    write_cutrun_shell_reference_catalog_with_sequence(root, genome_id, "ACGTACGTACGT")
}

fn write_cutrun_shell_reference_catalog_with_sequence(
    root: &Path,
    genome_id: &str,
    sequence: &str,
) -> String {
    let fasta_gz = root.join("toy.fa.gz");
    let ann_gz = root.join("toy.gtf.gz");
    let file = std::fs::File::create(&fasta_gz).expect("create fasta");
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    use std::io::Write as _;
    encoder
        .write_all(format!(">chr1\n{sequence}\n").as_bytes())
        .expect("write fasta");
    encoder.finish().expect("finish fasta");
    let file = std::fs::File::create(&ann_gz).expect("create gtf");
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    encoder
        .write_all(
            format!(
                "chr1\tsrc\tgene\t1\t{}\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
                sequence.len()
            )
            .as_bytes(),
        )
        .expect("write gtf");
    encoder.finish().expect("finish gtf");
    let cache_dir = root.join("reference_cache");
    let catalog_path = root.join("reference_catalog.json");
    fs::write(
        &catalog_path,
        format!(
            r#"{{
  "{genome_id}": {{
    "description": "toy genome",
    "sequence_remote": "{}",
    "annotations_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            shell_file_url(&fasta_gz),
            shell_file_url(&ann_gz),
            cache_dir.display()
        ),
    )
    .expect("write reference catalog");
    catalog_path.to_string_lossy().to_string()
}

fn prepare_cutrun_shell_anchor(
    engine: &mut GentleEngine,
    genome_id: &str,
    catalog_path: &str,
    seq_id: &str,
) {
    prepare_cutrun_shell_anchor_range(engine, genome_id, catalog_path, seq_id, 3, 10);
}

fn prepare_cutrun_shell_anchor_range(
    engine: &mut GentleEngine,
    genome_id: &str,
    catalog_path: &str,
    seq_id: &str,
    start_1based: usize,
    end_1based: usize,
) {
    engine
        .apply(Operation::PrepareGenome {
            genome_id: genome_id.to_string(),
            catalog_path: Some(catalog_path.to_string()),
            cache_dir: None,
            timeout_seconds: None,
        })
        .expect("prepare CUT&RUN shell genome");
    engine
        .apply(Operation::ExtractGenomeRegion {
            genome_id: genome_id.to_string(),
            chromosome: "chr1".to_string(),
            start_1based,
            end_1based,
            output_id: Some(seq_id.to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path.to_string()),
            cache_dir: None,
        })
        .expect("extract CUT&RUN shell anchor");
}

fn write_cutrun_shell_prepare_activity(
    install_dir: &Path,
    dataset_id: &str,
    lifecycle_status: &str,
    updated_at_unix_ms: u128,
    last_error: Option<&str>,
) {
    fs::create_dir_all(install_dir).expect("create CUT&RUN shell install dir for activity");
    let status_path = install_dir.join(".prepare_activity.json");
    let lock_path = install_dir.join(".prepare_activity.lock");
    let status = crate::engine::SharedAssetActivityStatus {
        resource_key: format!("cutrun_dataset:{dataset_id}"),
        display_name: dataset_id.to_string(),
        status_path: status_path.display().to_string(),
        lock_path: Some(lock_path.display().to_string()),
        lifecycle_status: lifecycle_status.to_string(),
        phase: Some("materializing_asset".to_string()),
        item: Some("toy_peaks.bed".to_string()),
        bytes_done: 64,
        bytes_total: Some(256),
        percent: Some(25.0),
        started_at_unix_ms: updated_at_unix_ms.saturating_sub(1000),
        updated_at_unix_ms,
        finished_at_unix_ms: None,
        last_error: last_error.map(str::to_string),
        owner_pid: Some(std::process::id()),
    };
    let text = serde_json::to_string_pretty(&status)
        .expect("serialize CUT&RUN shell prepare activity status");
    fs::write(&status_path, &text).expect("write CUT&RUN shell prepare activity status");
    if lifecycle_status == "running" {
        fs::write(&lock_path, text).expect("write CUT&RUN shell prepare activity lock");
    }
}

#[cfg(unix)]
fn install_fake_bigwig_to_bedgraph(path: &Path, bedgraph_source: &Path) -> String {
    let script_path = path.join("fake_bigwig_to_bedgraph.sh");
    let script = format!(
        "#!/bin/sh\nif [ \"$#\" -ne 2 ]; then\n  echo \"expected INPUT.bw OUTPUT.bedGraph\" >&2\n  exit 2\nfi\ncp \"{}\" \"$2\"\n",
        bedgraph_source.display()
    );
    fs::write(&script_path, script).expect("write fake bigwig converter");
    let mut perms = fs::metadata(&script_path)
        .expect("metadata fake bigwig converter")
        .permissions();
    perms.set_mode(0o755);
    fs::set_permissions(&script_path, perms).expect("chmod fake bigwig converter");
    script_path.display().to_string()
}

fn synthetic_scf_bytes(called_bases: &[u8]) -> Vec<u8> {
    let samples = 4u32;
    let sample_size = 1u32;
    let samples_offset = 128u32;
    let sample_section_len = samples * 4 * sample_size;
    let bases = called_bases.len() as u32;
    let bases_offset = samples_offset + sample_section_len;
    let base_section_len = bases * 12;
    let comments = b"NAME=shell_trace\nRUN=shell_demo\nMACH=staden_test\n\0";
    let comments_offset = bases_offset + base_section_len;
    let total_len = comments_offset as usize + comments.len();
    let mut bytes = vec![0u8; total_len];
    bytes[0..4].copy_from_slice(b".scf");
    bytes[4..8].copy_from_slice(&samples.to_be_bytes());
    bytes[8..12].copy_from_slice(&samples_offset.to_be_bytes());
    bytes[12..16].copy_from_slice(&bases.to_be_bytes());
    bytes[16..20].copy_from_slice(&0u32.to_be_bytes());
    bytes[20..24].copy_from_slice(&0u32.to_be_bytes());
    bytes[24..28].copy_from_slice(&bases_offset.to_be_bytes());
    bytes[28..32].copy_from_slice(&(comments.len() as u32).to_be_bytes());
    bytes[32..36].copy_from_slice(&comments_offset.to_be_bytes());
    bytes[36..40].copy_from_slice(b"3.00");
    bytes[40..44].copy_from_slice(&sample_size.to_be_bytes());
    bytes[44..48].copy_from_slice(&0u32.to_be_bytes());
    bytes[48..52].copy_from_slice(&0u32.to_be_bytes());
    bytes[52..56].copy_from_slice(&0u32.to_be_bytes());

    let base_start = bases_offset as usize;
    for idx in 0..bases as usize {
        let peak_index = (idx as u32) + 1;
        let offset = base_start + idx * 4;
        bytes[offset..offset + 4].copy_from_slice(&peak_index.to_be_bytes());
    }
    let prob_a_offset = base_start + bases as usize * 4;
    let prob_c_offset = prob_a_offset + bases as usize;
    let prob_g_offset = prob_c_offset + bases as usize;
    let prob_t_offset = prob_g_offset + bases as usize;
    let called_bases_offset = prob_t_offset + bases as usize;
    for (idx, base) in called_bases.iter().enumerate() {
        match *base {
            b'A' => bytes[prob_a_offset + idx] = 70,
            b'C' => bytes[prob_c_offset + idx] = 71,
            b'G' => bytes[prob_g_offset + idx] = 72,
            b'T' => bytes[prob_t_offset + idx] = 73,
            _ => {}
        }
    }
    bytes[called_bases_offset..called_bases_offset + called_bases.len()]
        .copy_from_slice(called_bases);
    bytes[comments_offset as usize..comments_offset as usize + comments.len()]
        .copy_from_slice(comments);
    bytes
}

fn write_shell_prepared_cache_install(root: &Path, genome_id: &str) -> std::path::PathBuf {
    let install_dir = root.join(genome_id.to_ascii_lowercase());
    std::fs::create_dir_all(install_dir.join("blastdb")).expect("create install dir");
    std::fs::write(install_dir.join("sequence.fa"), ">chr1\nACGTACGT\n").expect("sequence");
    std::fs::write(
        install_dir.join("annotation.gtf"),
        "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\";\n",
    )
    .expect("annotation");
    std::fs::write(install_dir.join("sequence.fa.fai"), "chr1\t8\t6\t8\t9\n").expect("fai");
    std::fs::write(install_dir.join("genes.json"), "[]").expect("genes");
    std::fs::write(install_dir.join("blastdb").join("genome.nhr"), "nhr").expect("nhr");
    std::fs::write(install_dir.join("blastdb").join("genome.nin"), "nin").expect("nin");
    std::fs::write(install_dir.join("blastdb").join("genome.nsq"), "nsq").expect("nsq");
    let manifest = serde_json::json!({
        "genome_id": genome_id,
        "sequence_source": install_dir.join("sequence.fa").display().to_string(),
        "annotation_source": install_dir.join("annotation.gtf").display().to_string(),
        "sequence_source_type": "local",
        "annotation_source_type": "local",
        "sequence_path": install_dir.join("sequence.fa").display().to_string(),
        "annotation_path": install_dir.join("annotation.gtf").display().to_string(),
        "fasta_index_path": install_dir.join("sequence.fa.fai").display().to_string(),
        "gene_index_path": install_dir.join("genes.json").display().to_string(),
        "blast_db_prefix": install_dir.join("blastdb").join("genome").display().to_string(),
        "blast_index_executable": "makeblastdb",
        "blast_indexed_at_unix_ms": 123,
        "installed_at_unix_ms": 456
    });
    std::fs::write(
        install_dir.join("manifest.json"),
        serde_json::to_string_pretty(&manifest).expect("serialize manifest"),
    )
    .expect("manifest");
    install_dir
}

#[cfg(unix)]
fn install_fake_primer3(path: &Path, fixture_path: &Path) -> String {
    let script_path = path.join("fake_primer3.sh");
    let script = format!(
        "#!/bin/sh\nif [ \"$1\" = \"--version\" ]; then\n  echo \"primer3_core synthetic-fixture 2.6.1\"\n  exit 0\nfi\ncat \"{}\"\n",
        fixture_path.display()
    );
    std::fs::write(&script_path, script).expect("write fake primer3");
    let mut perms = std::fs::metadata(&script_path)
        .expect("metadata fake primer3")
        .permissions();
    perms.set_mode(0o755);
    std::fs::set_permissions(&script_path, perms).expect("chmod fake primer3");
    script_path.display().to_string()
}

#[test]
fn parse_help_with_topic_and_options() {
    let cmd = parse_shell_line("help candidates generate --format json --interface cli-shell")
        .expect("parse help with topic/options");
    match cmd {
        ShellCommand::Help {
            topic,
            format,
            interface_filter,
        } => {
            assert_eq!(
                topic,
                vec!["candidates".to_string(), "generate".to_string()]
            );
            assert_eq!(format, HelpOutputFormat::Json);
            assert_eq!(interface_filter.as_deref(), Some("cli-shell"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_cache_inspect_command() {
    let cmd = parse_shell_line(
        "cache inspect --both --cache-dir data/genomes --cache-dir data/helper_genomes",
    )
    .expect("parse cache inspect");
    match cmd {
        ShellCommand::CacheInspect { scope, cache_dirs } => {
            assert_eq!(scope.label(), "both");
            assert_eq!(cache_dirs.len(), 2);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_cache_clear_command() {
    let cmd = parse_shell_line(
        "cache clear derived-indexes-only --helpers --cache-dir data/helper_genomes --prepared-id localproject",
    )
    .expect("parse cache clear");
    match cmd {
        ShellCommand::CacheClear {
            mode,
            scope,
            cache_dirs,
            prepared_ids,
            prepared_paths,
            include_orphans,
        } => {
            assert_eq!(
                mode,
                crate::genomes::PreparedCacheCleanupMode::DerivedIndexesOnly
            );
            assert_eq!(scope.label(), "helpers");
            assert_eq!(cache_dirs, vec!["data/helper_genomes".to_string()]);
            assert_eq!(prepared_ids, vec!["localproject".to_string()]);
            assert!(prepared_paths.is_empty());
            assert!(!include_orphans);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_cache_clear_command_with_prepared_path() {
    let cmd = parse_shell_line(
        "cache clear selected-prepared --references --cache-dir data/genomes --prepared-path data/genomes/localproject",
    )
    .expect("parse cache clear with path");
    match cmd {
        ShellCommand::CacheClear {
            mode,
            scope,
            cache_dirs,
            prepared_ids,
            prepared_paths,
            include_orphans,
        } => {
            assert_eq!(
                mode,
                crate::genomes::PreparedCacheCleanupMode::SelectedPreparedInstalls
            );
            assert_eq!(scope.label(), "references");
            assert_eq!(cache_dirs, vec!["data/genomes".to_string()]);
            assert!(prepared_ids.is_empty());
            assert_eq!(
                prepared_paths,
                vec!["data/genomes/localproject".to_string()]
            );
            assert!(!include_orphans);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_cache_inspect_and_clear_are_structured() {
    let td = tempdir().expect("tempdir");
    let root = td.path().join("cache");
    let install_dir = write_shell_prepared_cache_install(&root, "ToyGenome");
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let inspect = execute_shell_command(
        &mut engine,
        &ShellCommand::CacheInspect {
            scope: crate::engine_shell::CacheCleanupScope::References,
            cache_dirs: vec![root.display().to_string()],
        },
    )
    .expect("inspect cache");
    assert!(!inspect.state_changed);
    assert_eq!(
        inspect.output["schema"].as_str(),
        Some("gentle.prepared_cache_inspection.v1")
    );
    assert_eq!(inspect.output["entry_count"].as_u64(), Some(1));

    let clear = execute_shell_command(
        &mut engine,
        &ShellCommand::CacheClear {
            mode: crate::genomes::PreparedCacheCleanupMode::DerivedIndexesOnly,
            scope: crate::engine_shell::CacheCleanupScope::References,
            cache_dirs: vec![root.display().to_string()],
            prepared_ids: vec!["ToyGenome".to_string()],
            prepared_paths: vec![],
            include_orphans: false,
        },
    )
    .expect("clear cache");
    assert!(!clear.state_changed);
    assert_eq!(
        clear.output["schema"].as_str(),
        Some("gentle.prepared_cache_cleanup.v1")
    );
    assert!(install_dir.join("sequence.fa").exists());
    assert!(!install_dir.join("sequence.fa.fai").exists());
    assert!(!install_dir.join("genes.json").exists());
}

#[test]
fn execute_cache_clear_can_target_duplicate_ids_by_prepared_path() {
    let td = tempdir().expect("tempdir");
    let root_a = td.path().join("cache_a");
    let root_b = td.path().join("cache_b");
    let install_a = write_shell_prepared_cache_install(&root_a, "ToyGenome");
    let install_b = write_shell_prepared_cache_install(&root_b, "ToyGenome");
    let expected_selected_path = std::fs::canonicalize(&install_a)
        .unwrap_or_else(|_| install_a.clone())
        .to_string_lossy()
        .into_owned();
    let mut engine = GentleEngine::from_state(ProjectState::default());

    let clear = execute_shell_command(
        &mut engine,
        &ShellCommand::CacheClear {
            mode: crate::genomes::PreparedCacheCleanupMode::SelectedPreparedInstalls,
            scope: crate::engine_shell::CacheCleanupScope::References,
            cache_dirs: vec![root_a.display().to_string(), root_b.display().to_string()],
            prepared_ids: vec![],
            prepared_paths: vec![install_a.display().to_string()],
            include_orphans: false,
        },
    )
    .expect("clear cache by path");

    assert!(!clear.state_changed);
    let selected_path = clear.output["selected_prepared_paths"][0]
        .as_str()
        .expect("selected prepared path");
    assert_eq!(selected_path, expected_selected_path);
    assert!(!install_a.exists());
    assert!(install_b.exists());
}

#[test]
fn execute_help_topic_json() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::Help {
            topic: vec!["candidates".to_string(), "score-distance".to_string()],
            format: HelpOutputFormat::Json,
            interface_filter: None,
        },
    )
    .expect("execute help topic json");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["doc"]["path"].as_str(),
        Some("candidates score-distance")
    );
}

#[test]
fn parse_workflow_payload_keeps_whitespace() {
    let cmd = parse_shell_line("workflow '{ \"run_id\": \"x\", \"ops\": [] }'")
        .expect("workflow command parse");
    match cmd {
        ShellCommand::Workflow { payload } => {
            assert!(payload.contains("\"run_id\""));
            assert!(payload.contains("\"ops\""));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_seq_confirm_run_command() {
    let cmd = parse_shell_line(
        "seq-confirm run construct --baseline baseline_ref --reads read_a,read_b --trace-id abi_trace --trace-ids abi_trace_2,abi_trace_3 --junction 8 --junction-flank 4 --mode local --min-identity 0.90 --min-target-coverage 0.75 --no-reverse-complement --report-id construct_check",
    )
    .expect("parse seq-confirm run");
    match cmd {
        ShellCommand::SeqConfirmRun {
            expected_seq_id,
            baseline_seq_id,
            read_seq_ids,
            trace_ids,
            targets,
            alignment_mode,
            min_identity_fraction,
            min_target_coverage_fraction,
            allow_reverse_complement,
            report_id,
            ..
        } => {
            assert_eq!(expected_seq_id, "construct");
            assert_eq!(baseline_seq_id, Some("baseline_ref".to_string()));
            assert_eq!(read_seq_ids, vec!["read_a", "read_b"]);
            assert_eq!(trace_ids, vec!["abi_trace", "abi_trace_2", "abi_trace_3"]);
            assert_eq!(targets.len(), 1);
            assert_eq!(targets[0].kind, SequencingConfirmationTargetKind::Junction);
            assert_eq!(targets[0].start_0based, 4);
            assert_eq!(targets[0].end_0based_exclusive, 12);
            assert_eq!(targets[0].junction_left_end_0based, Some(8));
            assert_eq!(alignment_mode, PairwiseAlignmentMode::Local);
            assert!((min_identity_fraction - 0.90).abs() < 1e-9);
            assert!((min_target_coverage_fraction - 0.75).abs() < 1e-9);
            assert!(!allow_reverse_complement);
            assert_eq!(report_id.as_deref(), Some("construct_check"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_seq_trace_import_and_list_commands() {
    let import =
        parse_shell_line("seq-trace import /tmp/demo.ab1 --trace-id abi_trace --seq-id construct")
            .expect("parse seq-trace import");
    match import {
        ShellCommand::SeqTraceImport {
            path,
            trace_id,
            seq_id,
        } => {
            assert_eq!(path, "/tmp/demo.ab1");
            assert_eq!(trace_id.as_deref(), Some("abi_trace"));
            assert_eq!(seq_id.as_deref(), Some("construct"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let list = parse_shell_line("seq-trace list construct").expect("parse seq-trace list");
    match list {
        ShellCommand::SeqTraceList { seq_id } => {
            assert_eq!(seq_id.as_deref(), Some("construct"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let show = parse_shell_line("seq-trace show abi_trace").expect("parse seq-trace show");
    match show {
        ShellCommand::SeqTraceShow { trace_id } => {
            assert_eq!(trace_id, "abi_trace");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_seq_primer_suggest_command() {
    let cmd = parse_shell_line(
        "seq-primer suggest construct --primers primer_a,primer_b --confirmation-report construct_report --min-3prime-anneal-bp 20 --predicted-read-length-bp 650",
    )
    .expect("parse seq-primer suggest");
    match cmd {
        ShellCommand::SeqPrimerSuggest {
            expected_seq_id,
            primer_seq_ids,
            confirmation_report_id,
            min_3prime_anneal_bp,
            predicted_read_length_bp,
        } => {
            assert_eq!(expected_seq_id, "construct");
            assert_eq!(primer_seq_ids, vec!["primer_a", "primer_b"]);
            assert_eq!(confirmation_report_id.as_deref(), Some("construct_report"));
            assert_eq!(min_3prime_anneal_bp, 20);
            assert_eq!(predicted_read_length_bp, 650);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_reverse_translate_command_family() {
    let run = parse_shell_line(
        "reverse-translate run prot --output-id prot_coding --speed-profile ecoli --speed-mark slow --translation-table 11 --target-anneal-tm-c 58.0 --anneal-window-bp 9",
    )
    .expect("parse reverse-translate run");
    assert!(matches!(
        run,
        ShellCommand::ReverseTranslateRun {
            seq_id,
            output_id,
            speed_profile,
            speed_mark,
            translation_table,
            target_anneal_tm_c,
            anneal_window_bp,
        }
            if seq_id == "prot"
                && output_id.as_deref() == Some("prot_coding")
                && speed_profile == Some(TranslationSpeedProfile::Ecoli)
                && speed_mark == Some(TranslationSpeedMark::Slow)
                && translation_table == Some(11)
                && target_anneal_tm_c == Some(58.0)
                && anneal_window_bp == Some(9)
    ));

    let list = parse_shell_line("reverse-translate list-reports prot")
        .expect("parse reverse-translate list-reports");
    assert!(matches!(
        list,
        ShellCommand::ReverseTranslateListReports { protein_seq_id }
            if protein_seq_id.as_deref() == Some("prot")
    ));

    let show = parse_shell_line("reverse-translate show-report revtx_1")
        .expect("parse reverse-translate show-report");
    assert!(matches!(
        show,
        ShellCommand::ReverseTranslateShowReport { report_id } if report_id == "revtx_1"
    ));

    let export = parse_shell_line("reverse-translate export-report revtx_1 out.json")
        .expect("parse reverse-translate export-report");
    assert!(matches!(
        export,
        ShellCommand::ReverseTranslateExportReport { report_id, path }
            if report_id == "revtx_1" && path == "out.json"
    ));
}

#[test]
fn parse_construct_reasoning_protein_handoff_command_family() {
    let build = parse_shell_line(
        "construct-reasoning build-protein-dna-handoff seq_a prot_a --transcript TX1 --projection-id proj_1 --ensembl-entry ENSPTOY1 --feature-query junction --ranking-goal balanced_provenance --speed-profile ecoli --speed-mark slow --translation-table 11 --target-anneal-tm-c 58.0 --anneal-window-bp 9 --objective-id obj_1 --graph-id graph_1",
    )
    .expect("parse construct-reasoning build-protein-dna-handoff");
    assert!(matches!(
        build,
        ShellCommand::ConstructReasoningBuildProteinDnaHandoff {
            seq_id,
            protein_seq_id,
            transcript_filter,
            projection_id,
            ensembl_entry_id,
            feature_query,
            ranking_goal,
            speed_profile,
            speed_mark,
            translation_table,
            target_anneal_tm_c,
            anneal_window_bp,
            objective_id,
            graph_id,
        }
            if seq_id == "seq_a"
                && protein_seq_id == "prot_a"
                && transcript_filter.as_deref() == Some("TX1")
                && projection_id.as_deref() == Some("proj_1")
                && ensembl_entry_id.as_deref() == Some("ENSPTOY1")
                && feature_query.as_deref() == Some("junction")
                && ranking_goal == ProteinToDnaHandoffRankingGoal::BalancedProvenance
                && speed_profile == Some(TranslationSpeedProfile::Ecoli)
                && speed_mark == Some(TranslationSpeedMark::Slow)
                && translation_table == Some(11)
                && target_anneal_tm_c == Some(58.0)
                && anneal_window_bp == Some(9)
                && objective_id.as_deref() == Some("obj_1")
                && graph_id.as_deref() == Some("graph_1")
    ));

    let list = parse_shell_line("construct-reasoning list-graphs seq_a")
        .expect("parse construct-reasoning list-graphs");
    assert!(matches!(
        list,
        ShellCommand::ConstructReasoningListGraphs { seq_id }
            if seq_id.as_deref() == Some("seq_a")
    ));

    let show = parse_shell_line("construct-reasoning show-graph graph_1")
        .expect("parse construct-reasoning show-graph");
    assert!(matches!(
        show,
        ShellCommand::ConstructReasoningShowGraph { graph_id } if graph_id == "graph_1"
    ));

    let set_status = parse_shell_line(
        "construct-reasoning set-annotation-status graph_1 annotation_promoter accepted",
    )
    .expect("parse construct-reasoning set-annotation-status");
    assert!(matches!(
        set_status,
        ShellCommand::ConstructReasoningSetAnnotationStatus {
            graph_id,
            annotation_id,
            editable_status,
        } if graph_id == "graph_1"
            && annotation_id == "annotation_promoter"
            && editable_status == EditableStatus::Accepted
    ));

    let write_annotation =
        parse_shell_line("construct-reasoning write-annotation graph_1 annotation_promoter")
            .expect("parse construct-reasoning write-annotation");
    assert!(matches!(
        write_annotation,
        ShellCommand::ConstructReasoningWriteAnnotation {
            graph_id,
            annotation_id,
        } if graph_id == "graph_1" && annotation_id == "annotation_promoter"
    ));

    let export = parse_shell_line("construct-reasoning export-graph graph_1 out.json")
        .expect("parse construct-reasoning export-graph");
    assert!(matches!(
        export,
        ShellCommand::ConstructReasoningExportGraph { graph_id, path }
            if graph_id == "graph_1" && path == "out.json"
    ));
}

#[test]
fn parse_seq_primer_suggest_command_allows_report_only_mode() {
    let cmd = parse_shell_line(
        "seq-primer suggest construct --confirmation-report construct_report --min-3prime-anneal-bp 20 --predicted-read-length-bp 650",
    )
    .expect("parse report-only seq-primer suggest");
    match cmd {
        ShellCommand::SeqPrimerSuggest {
            expected_seq_id,
            primer_seq_ids,
            confirmation_report_id,
            min_3prime_anneal_bp,
            predicted_read_length_bp,
        } => {
            assert_eq!(expected_seq_id, "construct");
            assert!(primer_seq_ids.is_empty());
            assert_eq!(confirmation_report_id.as_deref(), Some("construct_report"));
            assert_eq!(min_3prime_anneal_bp, 20);
            assert_eq!(predicted_read_length_bp, 650);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_seq_trace_import_list_and_show_round_trip() {
    let fixture = sequencing_confirmation_fixture_path("3100.ab1");
    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    let mut engine = GentleEngine::from_state(state);

    let imported = execute_shell_command(
        &mut engine,
        &ShellCommand::SeqTraceImport {
            path: fixture,
            trace_id: Some("abi_trace".to_string()),
            seq_id: Some("construct".to_string()),
        },
    )
    .expect("execute seq-trace import");
    assert!(imported.state_changed);
    assert_eq!(
        imported.output["import_report"]["trace_id"].as_str(),
        Some("abi_trace")
    );
    assert_eq!(imported.output["trace"]["format"].as_str(), Some("abi_ab1"));

    let listed = execute_shell_command(
        &mut engine,
        &ShellCommand::SeqTraceList {
            seq_id: Some("construct".to_string()),
        },
    )
    .expect("execute seq-trace list");
    assert_eq!(listed.output["trace_count"].as_u64(), Some(1));
    assert_eq!(
        listed.output["traces"][0]["trace_id"].as_str(),
        Some("abi_trace")
    );

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::SeqTraceShow {
            trace_id: "abi_trace".to_string(),
        },
    )
    .expect("execute seq-trace show");
    assert_eq!(
        shown.output["trace"]["trace_id"].as_str(),
        Some("abi_trace")
    );
    assert!(
        shown.output["trace"]["called_bases"]
            .as_str()
            .is_some_and(|value| !value.is_empty())
    );
}

#[test]
fn execute_seq_confirm_run_returns_persisted_report() {
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
    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::SeqConfirmRun {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_junction".to_string()],
            trace_ids: vec![],
            targets: vec![SequencingConfirmationTargetSpec {
                target_id: "junction_1".to_string(),
                label: "Insert junction".to_string(),
                kind: SequencingConfirmationTargetKind::Junction,
                start_0based: 4,
                end_0based_exclusive: 12,
                junction_left_end_0based: Some(8),
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
            report_id: Some("construct_check".to_string()),
        },
    )
    .expect("execute seq-confirm run");
    assert!(run.state_changed);
    assert_eq!(
        run.output["report"]["report_id"].as_str(),
        Some("construct_check")
    );
    assert_eq!(
        run.output["report"]["overall_status"].as_str(),
        Some("confirmed")
    );
    assert_eq!(
        run.output["report"]["targets"][0]["status"].as_str(),
        Some("confirmed")
    );

    let listed = execute_shell_command(
        &mut engine,
        &ShellCommand::SeqConfirmListReports {
            expected_seq_id: Some("construct".to_string()),
        },
    )
    .expect("list sequencing-confirmation reports");
    assert_eq!(listed.output["report_count"].as_u64(), Some(1));
    assert_eq!(
        listed.output["reports"][0]["overall_status"].as_str(),
        Some("confirmed")
    );
}

#[test]
fn execute_seq_confirm_run_accepts_imported_trace_evidence() {
    let td = tempdir().expect("tempdir");
    let trace_path = td.path().join("junction_trace.scf");
    fs::write(&trace_path, synthetic_scf_bytes(b"CCGTAACC")).expect("write scf");

    let mut state = ProjectState::default();
    state.sequences.insert(
        "construct".to_string(),
        DNAsequence::from_sequence("AAAACCGTAACCTTTT").expect("construct"),
    );
    let mut engine = GentleEngine::from_state(state);

    execute_shell_command(
        &mut engine,
        &ShellCommand::SeqTraceImport {
            path: trace_path.to_str().expect("utf-8 path").to_string(),
            trace_id: Some("junction_trace".to_string()),
            seq_id: None,
        },
    )
    .expect("import trace");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::SeqConfirmRun {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec![],
            trace_ids: vec!["junction_trace".to_string()],
            targets: vec![SequencingConfirmationTargetSpec {
                target_id: "junction_1".to_string(),
                label: "Insert junction".to_string(),
                kind: SequencingConfirmationTargetKind::Junction,
                start_0based: 4,
                end_0based_exclusive: 12,
                junction_left_end_0based: Some(8),
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
            report_id: Some("construct_trace".to_string()),
        },
    )
    .expect("execute seq-confirm from trace");

    assert_eq!(
        run.output["report"]["overall_status"].as_str(),
        Some("confirmed")
    );
    assert_eq!(
        run.output["report"]["trace_ids"][0].as_str(),
        Some("junction_trace")
    );
    assert_eq!(
        run.output["report"]["reads"][0]["evidence_kind"].as_str(),
        Some("trace")
    );
    assert_eq!(
        run.output["report"]["reads"][0]["trace_id"].as_str(),
        Some("junction_trace")
    );
}

#[test]
fn execute_seq_confirm_run_with_baseline_exposes_variant_rows() {
    let td = tempdir().expect("tempdir");
    let trace_path = td.path().join("expected_edit_trace.scf");
    fs::write(&trace_path, synthetic_scf_bytes(b"CCGTAACC")).expect("write scf");

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

    execute_shell_command(
        &mut engine,
        &ShellCommand::SeqTraceImport {
            path: trace_path.to_str().expect("utf-8 path").to_string(),
            trace_id: Some("expected_edit_trace".to_string()),
            seq_id: None,
        },
    )
    .expect("import trace");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::SeqConfirmRun {
            expected_seq_id: "expected".to_string(),
            baseline_seq_id: Some("baseline".to_string()),
            read_seq_ids: vec![],
            trace_ids: vec!["expected_edit_trace".to_string()],
            targets: vec![],
            alignment_mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            min_identity_fraction: 0.80,
            min_target_coverage_fraction: 1.0,
            allow_reverse_complement: true,
            report_id: Some("expected_edit_report".to_string()),
        },
    )
    .expect("execute seq-confirm from baseline-aware trace");

    assert_eq!(
        run.output["report"]["baseline_seq_id"].as_str(),
        Some("baseline")
    );
    assert_eq!(
        run.output["report"]["variants"].as_array().map(Vec::len),
        Some(1)
    );
    assert_eq!(
        run.output["report"]["variants"][0]["classification"].as_str(),
        Some("intended_edit_confirmed")
    );
}

#[test]
fn execute_reverse_translate_commands_store_list_show_and_export_reports() {
    let td = tempdir().expect("tempdir");
    let export_path = td.path().join("reverse_translation_report.json");

    let mut protein = DNAsequence::from_sequence("MKP").expect("protein");
    protein.set_name("Toy protein");
    protein.set_molecule_type("protein");

    let mut state = ProjectState::default();
    state.sequences.insert("prot".to_string(), protein);
    let mut engine = GentleEngine::from_state(state);

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::ReverseTranslateRun {
            seq_id: "prot".to_string(),
            output_id: Some("prot_coding".to_string()),
            speed_profile: Some(TranslationSpeedProfile::Ecoli),
            speed_mark: Some(TranslationSpeedMark::Slow),
            translation_table: Some(11),
            target_anneal_tm_c: Some(58.0),
            anneal_window_bp: Some(9),
        },
    )
    .expect("execute reverse-translate run");
    assert!(run.state_changed);
    assert_eq!(
        run.output["report"]["protein_seq_id"].as_str(),
        Some("prot")
    );
    assert_eq!(
        run.output["report"]["coding_seq_id"].as_str(),
        Some("prot_coding")
    );
    let report_id = run.output["report"]["report_id"]
        .as_str()
        .expect("reverse-translation report id")
        .to_string();

    let list = execute_shell_command(
        &mut engine,
        &ShellCommand::ReverseTranslateListReports {
            protein_seq_id: Some("prot".to_string()),
        },
    )
    .expect("list reverse-translation reports");
    assert_eq!(list.output["report_count"].as_u64(), Some(1));
    assert!(
        list.output["summary_rows"][0]
            .as_str()
            .is_some_and(|value| value.contains("gc="))
    );

    let show = execute_shell_command(
        &mut engine,
        &ShellCommand::ReverseTranslateShowReport {
            report_id: report_id.clone(),
        },
    )
    .expect("show reverse-translation report");
    assert!(
        show.output["summary"]
            .as_str()
            .is_some_and(|value| value.contains("speed=ecoli:slow"))
    );

    let export = execute_shell_command(
        &mut engine,
        &ShellCommand::ReverseTranslateExportReport {
            report_id: report_id.clone(),
            path: export_path.display().to_string(),
        },
    )
    .expect("export reverse-translation report");
    assert_eq!(
        export.output["report_id"].as_str(),
        Some(report_id.as_str())
    );
    assert!(export_path.exists());
}

#[test]
fn execute_construct_reasoning_protein_handoff_commands_store_list_show_and_export_graphs() {
    let td = tempdir().expect("tempdir");
    let export_path = td.path().join("protein_handoff_graph.json");

    let mut protein = DNAsequence::from_sequence("MKP").expect("protein");
    protein.set_name("Toy protein");
    protein.set_molecule_type("protein");

    let mut state = ProjectState::default();
    state.sequences.insert(
        "target".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("dna"),
    );
    state.sequences.insert("prot".to_string(), protein);
    let mut engine = GentleEngine::from_state(state);

    let build = execute_shell_command(
        &mut engine,
        &ShellCommand::ConstructReasoningBuildProteinDnaHandoff {
            seq_id: "target".to_string(),
            protein_seq_id: "prot".to_string(),
            transcript_filter: None,
            projection_id: None,
            ensembl_entry_id: Some("ENSPTOY1".to_string()),
            feature_query: None,
            ranking_goal: ProteinToDnaHandoffRankingGoal::BalancedProvenance,
            speed_profile: Some(TranslationSpeedProfile::Ecoli),
            speed_mark: Some(TranslationSpeedMark::Slow),
            translation_table: Some(11),
            target_anneal_tm_c: Some(58.0),
            anneal_window_bp: Some(9),
            objective_id: None,
            graph_id: Some("protein_handoff_shell".to_string()),
        },
    )
    .expect("execute construct-reasoning build-protein-dna-handoff");
    assert!(build.state_changed);
    assert_eq!(
        build.output["graph"]["graph_id"].as_str(),
        Some("protein_handoff_shell")
    );
    assert!(
        build.output["graph"]["candidates"]
            .as_array()
            .is_some_and(|rows| !rows.is_empty())
    );

    let list = execute_shell_command(
        &mut engine,
        &ShellCommand::ConstructReasoningListGraphs {
            seq_id: Some("target".to_string()),
        },
    )
    .expect("list construct-reasoning graphs");
    assert_eq!(list.output["graph_count"].as_u64(), Some(1));
    assert_eq!(
        list.output["graphs"][0]["contains_protein_to_dna_handoff"].as_bool(),
        Some(true)
    );

    let show = execute_shell_command(
        &mut engine,
        &ShellCommand::ConstructReasoningShowGraph {
            graph_id: "protein_handoff_shell".to_string(),
        },
    )
    .expect("show construct-reasoning graph");
    assert_eq!(
        show.output["graph"]["graph_id"].as_str(),
        Some("protein_handoff_shell")
    );
    assert_eq!(
        show.output["graph"]["candidates"][0]["protein_to_dna_handoff"]["strategy"].as_str(),
        Some("reverse_translated_synthetic")
    );
    assert!(
        show.output["summary"]["summary_lines"]
            .as_array()
            .is_some_and(|rows| !rows.is_empty())
    );

    let export = execute_shell_command(
        &mut engine,
        &ShellCommand::ConstructReasoningExportGraph {
            graph_id: "protein_handoff_shell".to_string(),
            path: export_path.display().to_string(),
        },
    )
    .expect("export construct-reasoning graph");
    assert_eq!(
        export.output["graph_id"].as_str(),
        Some("protein_handoff_shell")
    );
    assert!(export_path.exists());
}

#[test]
fn execute_construct_reasoning_show_graph_includes_adapter_capture_summary() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "adapter_capture_cli".to_string(),
        DNAsequence::from_sequence("GAATTCTCTAGAGCGGCCGCTTT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let objective = engine
        .upsert_construct_objective(ConstructObjective {
            title: "CLI adapter capture".to_string(),
            goal: "Expose adapter capture review in shell summaries".to_string(),
            adapter_restriction_capture_plans: vec![AdapterRestrictionCapturePlan {
                capture_id: "mcs_capture".to_string(),
                restriction_enzyme_name: "EcoRI".to_string(),
                blunt_insert_required: true,
                adapter_style: AdapterCaptureStyle::McsLike,
                protection_mode: AdapterCaptureProtectionMode::InsertMethylation,
                extra_retrieval_enzyme_names: vec!["XbaI".to_string(), "NotI".to_string()],
                notes: vec![],
            }],
            ..ConstructObjective::default()
        })
        .expect("objective");
    let graph = engine
        .build_construct_reasoning_graph(
            "adapter_capture_cli",
            Some(&objective.objective_id),
            Some("adapter_capture_cli_graph"),
        )
        .expect("build graph");

    let show = execute_shell_command(
        &mut engine,
        &ShellCommand::ConstructReasoningShowGraph {
            graph_id: graph.graph_id.clone(),
        },
    )
    .expect("show construct-reasoning graph");

    assert!(
        show.output["summary"]["summary_lines"]
            .as_array()
            .map(|rows| rows.iter().any(|row| {
                row.as_str().is_some_and(|text| {
                    text.contains("Adapter/linker restriction capture requires review")
                        && text.contains("review_needed")
                })
            }))
            .unwrap_or(false)
    );
    assert!(
        show.output["summary"]["fact_summaries"]
            .as_array()
            .map(|rows| rows.iter().any(|row| {
                row["fact_type"].as_str() == Some("adapter_restriction_capture_context")
                    && row["detail_lines"]
                        .as_array()
                        .map(|detail_rows| {
                            detail_rows.iter().any(|entry| {
                                entry.as_str().is_some_and(|text| {
                                    text.contains(
                                        "aggregate: mcs_capture: all_named_motifs_present",
                                    )
                                })
                            })
                        })
                        .unwrap_or(false)
                    && row["warning_lines"]
                        .as_array()
                        .map(|warning_rows| {
                            warning_rows.iter().any(|entry| {
                                entry.as_str().is_some_and(|text| {
                                    text == "All adapter motifs already occur on the insert"
                                })
                            })
                        })
                        .unwrap_or(false)
            }))
            .unwrap_or(false)
    );
}

#[test]
fn execute_construct_reasoning_set_annotation_status_updates_graph_and_summary() {
    let mut dna = DNAsequence::from_sequence("ATGCGTATGCGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "exon".into(),
        location: gb_io::seq::Location::simple_range(2, 10),
        qualifiers: vec![
            ("label".into(), Some("Confirmed exon".to_string())),
            ("evidence".into(), Some("supported by cDNA".to_string())),
        ],
    });
    dna.update_computed_features();

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("construct_reasoning_cli_status".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);
    let graph = engine
        .build_construct_reasoning_graph(
            "construct_reasoning_cli_status",
            None,
            Some("graph_annotation_status_cli"),
        )
        .expect("build graph");
    let annotation_id = graph
        .annotation_candidates
        .iter()
        .find(|candidate| candidate.role == ConstructRole::Exon)
        .map(|candidate| candidate.annotation_id.clone())
        .expect("annotation candidate");

    let output = execute_shell_command(
        &mut engine,
        &ShellCommand::ConstructReasoningSetAnnotationStatus {
            graph_id: graph.graph_id.clone(),
            annotation_id: annotation_id.clone(),
            editable_status: EditableStatus::Accepted,
        },
    )
    .expect("set construct-reasoning annotation status");

    assert!(output.state_changed);
    assert_eq!(output.output["editable_status"].as_str(), Some("accepted"));
    assert_eq!(
        output.output["annotation_candidate"]["editable_status"].as_str(),
        Some("accepted")
    );
    assert!(
        output.output["summary"]["summary_lines"]
            .as_array()
            .map(|rows| rows.iter().any(|row| {
                row.as_str()
                    .is_some_and(|text| text.contains("Annotation candidates: 1 accepted"))
            }))
            .unwrap_or(false)
    );
}

#[test]
fn execute_construct_reasoning_write_annotation_materializes_generated_feature() {
    let mut dna = DNAsequence::from_sequence(&"A".repeat(6001)).expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::Complement(Box::new(gb_io::seq::Location::simple_range(
            2000, 2613,
        ))),
        qualifiers: vec![
            ("gene".into(), Some("VKORC1".to_string())),
            ("transcript_id".into(), Some("ENSTVKORC1".to_string())),
            ("label".into(), Some("VKORC1-201".to_string())),
        ],
    });
    dna.update_computed_features();

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("construct_reasoning_cli_writeback".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);
    let graph = engine
        .build_construct_reasoning_graph(
            "construct_reasoning_cli_writeback",
            None,
            Some("graph_annotation_writeback_cli"),
        )
        .expect("build graph");
    let annotation_id = graph
        .annotation_candidates
        .iter()
        .find(|candidate| candidate.role == ConstructRole::Promoter)
        .map(|candidate| candidate.annotation_id.clone())
        .expect("generated promoter annotation candidate");
    engine
        .set_construct_reasoning_annotation_candidate_status(
            &graph.graph_id,
            &annotation_id,
            EditableStatus::Accepted,
        )
        .expect("accept annotation candidate");

    let output = execute_shell_command(
        &mut engine,
        &ShellCommand::ConstructReasoningWriteAnnotation {
            graph_id: graph.graph_id.clone(),
            annotation_id: annotation_id.clone(),
        },
    )
    .expect("write back construct-reasoning annotation");

    assert!(output.state_changed);
    assert_eq!(output.output["writeback"]["created"].as_bool(), Some(true));
    assert_eq!(
        output.output["writeback"]["annotation_id"].as_str(),
        Some(annotation_id.as_str())
    );
    assert!(
        output.output["summary"]["summary_lines"]
            .as_array()
            .map(|rows| rows.iter().any(|row| {
                row.as_str()
                    .is_some_and(|text| text.contains("Annotation candidates: 1 accepted"))
            }))
            .unwrap_or(false)
    );

    let persisted_dna = engine
        .state()
        .sequences
        .get("construct_reasoning_cli_writeback")
        .expect("persisted sequence");
    assert!(persisted_dna.features().iter().any(|feature| {
        feature
            .qualifier_values("construct_reasoning_annotation_id")
            .any(|value| value == annotation_id.as_str())
    }));
}

#[test]
fn execute_seq_primer_suggest_returns_overlay_report() {
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
    state.sequences.insert(
        "primer_rev".to_string(),
        DNAsequence::from_sequence("TTTGGTT").expect("primer"),
    );
    let mut engine = GentleEngine::from_state(state);

    execute_shell_command(
        &mut engine,
        &ShellCommand::SeqConfirmRun {
            expected_seq_id: "construct".to_string(),
            baseline_seq_id: None,
            read_seq_ids: vec!["read_junction".to_string()],
            trace_ids: vec![],
            targets: vec![SequencingConfirmationTargetSpec {
                target_id: "junction_1".to_string(),
                label: "Insert junction".to_string(),
                kind: SequencingConfirmationTargetKind::Junction,
                start_0based: 4,
                end_0based_exclusive: 12,
                junction_left_end_0based: Some(8),
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
            report_id: Some("construct_check".to_string()),
        },
    )
    .expect("persist confirmation report");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::SeqPrimerSuggest {
            expected_seq_id: "construct".to_string(),
            primer_seq_ids: vec!["primer_fwd".to_string(), "primer_rev".to_string()],
            confirmation_report_id: Some("construct_check".to_string()),
            min_3prime_anneal_bp: 4,
            predicted_read_length_bp: 10,
        },
    )
    .expect("execute seq-primer suggest");

    assert!(!run.state_changed);
    assert!(
        run.output["suggestion_count"]
            .as_u64()
            .is_some_and(|value| value >= 2)
    );
    assert_eq!(
        run.output["report"]["confirmation_report_id"].as_str(),
        Some("construct_check")
    );
    assert!(
        run.output["report"]["suggestions"]
            .as_array()
            .is_some_and(|rows| rows.iter().any(|row| row["covered_target_ids"]
                .as_array()
                .is_some_and(|targets| targets
                    .iter()
                    .any(|value| value.as_str() == Some("junction_1")))))
    );
}

#[test]
fn execute_seq_primer_suggest_reports_guidance_for_unresolved_target() {
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

    execute_shell_command(
        &mut engine,
        &ShellCommand::SeqConfirmRun {
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
        },
    )
    .expect("persist unresolved report");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::SeqPrimerSuggest {
            expected_seq_id: "construct".to_string(),
            primer_seq_ids: vec!["primer_good".to_string(), "primer_near".to_string()],
            confirmation_report_id: Some("unresolved_gap".to_string()),
            min_3prime_anneal_bp: 6,
            predicted_read_length_bp: 150,
        },
    )
    .expect("execute seq-primer suggest with guidance");

    assert_eq!(
        run.output["report"]["problem_guidance"][0]["problem_id"].as_str(),
        Some("gap_target")
    );
    assert_eq!(
        run.output["report"]["problem_guidance"][0]["recommended_primer_seq_id"].as_str(),
        Some("primer_good")
    );
    assert_eq!(
        run.output["report"]["problem_guidance"][0]["candidate_count"].as_u64(),
        Some(2)
    );
}

#[test]
fn execute_seq_primer_suggest_report_only_mode_returns_fresh_proposals() {
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

    execute_shell_command(
        &mut engine,
        &ShellCommand::SeqConfirmRun {
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
        },
    )
    .expect("persist unresolved report");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::SeqPrimerSuggest {
            expected_seq_id: "construct".to_string(),
            primer_seq_ids: vec![],
            confirmation_report_id: Some("report_only_gap".to_string()),
            min_3prime_anneal_bp: 6,
            predicted_read_length_bp: 80,
        },
    )
    .expect("execute report-only seq-primer suggest");

    assert_eq!(run.output["suggestion_count"].as_u64(), Some(0));
    assert_eq!(run.output["proposal_count"].as_u64(), Some(1));
    assert_eq!(
        run.output["report"]["proposals"][0]["problem_id"].as_str(),
        Some("gap_target")
    );
    assert!(
        run.output["report"]["proposals"][0]["reason"]
            .as_str()
            .is_some_and(|value| value.contains("No existing primer covered"))
    );
}

#[test]
fn parse_workflow_json_payload_accepts_raw_workflow() {
    let payload = r#"{ "run_id": "raw", "ops": [] }"#;
    let workflow = parse_workflow_json_payload(payload).expect("parse raw workflow");
    assert_eq!(workflow.run_id, "raw");
    assert!(workflow.ops.is_empty());
}

#[test]
fn parse_workflow_json_payload_accepts_wrapped_example() {
    let payload = r#"{
            "schema": "gentle.workflow_example.v1",
            "id": "demo",
            "title": "Demo",
            "workflow": {
                "run_id": "wrapped",
                "ops": []
            }
        }"#;
    let workflow = parse_workflow_json_payload(payload).expect("parse wrapped workflow");
    assert_eq!(workflow.run_id, "wrapped");
    assert!(workflow.ops.is_empty());
}

#[test]
fn parse_json_payload_reads_existing_file_without_at_prefix() {
    let dir = tempdir().expect("temp dir");
    let path = write_demo_workflow_json(dir.path(), "workflow.json", "from_file");
    let loaded = parse_json_payload(path.to_str().expect("utf-8 path")).expect("parse payload");
    assert_eq!(loaded, r#"{"run_id":"from_file","ops":[]}"#);
}

#[test]
fn parse_json_payload_strips_shebang_from_file() {
    let dir = tempdir().expect("temp dir");
    let path = write_demo_workflow_with_shebang(dir.path(), "workflow.gsh", "from_shebang");
    let loaded = parse_json_payload(path.to_str().expect("utf-8 path")).expect("parse payload");
    assert_eq!(loaded.trim(), r#"{"run_id":"from_shebang","ops":[]}"#);
}

#[test]
fn load_macro_script_reads_existing_file_without_at_prefix() {
    let dir = tempdir().expect("temp dir");
    let path = dir.path().join("macro.gsh");
    fs::write(&path, "state-summary").expect("write macro");
    let loaded = load_macro_script(path.to_str().expect("utf-8 path"), "macros run")
        .expect("load macro script");
    assert_eq!(loaded, "state-summary");
}

#[test]
fn parse_render_pool_gel_with_ladders() {
    let cmd = parse_shell_line("render-pool-gel-svg a,b out.svg --ladders 1kb,100bp")
        .expect("parse command");
    match cmd {
        ShellCommand::RenderPoolGelSvg {
            inputs,
            output,
            ladders,
            container_ids,
            arrangement_id,
            conditions,
        } => {
            assert_eq!(inputs, vec!["a".to_string(), "b".to_string()]);
            assert_eq!(output, "out.svg".to_string());
            assert_eq!(ladders, Some(vec!["1kb".to_string(), "100bp".to_string()]));
            assert_eq!(container_ids, None);
            assert_eq!(arrangement_id, None);
            assert_eq!(conditions.agarose_percent, 1.0);
            assert!(conditions.topology_aware);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_render_pool_gel_from_arrangement() {
    let cmd = parse_shell_line("render-pool-gel-svg - out.svg --arrangement arrangement-2")
        .expect("parse command");
    match cmd {
        ShellCommand::RenderPoolGelSvg {
            inputs,
            output,
            ladders,
            container_ids,
            arrangement_id,
            conditions,
        } => {
            assert!(inputs.is_empty());
            assert_eq!(output, "out.svg".to_string());
            assert_eq!(ladders, None);
            assert_eq!(container_ids, None);
            assert_eq!(arrangement_id, Some("arrangement-2".to_string()));
            assert_eq!(conditions.buffer_model.as_str(), "tae");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_render_gel_svg_alias() {
    let cmd = parse_shell_line("render-gel-svg - out.svg --arrangement arrangement-2")
        .expect("parse command");
    match cmd {
        ShellCommand::RenderPoolGelSvg {
            inputs,
            output,
            ladders,
            container_ids,
            arrangement_id,
            conditions,
        } => {
            assert!(inputs.is_empty());
            assert_eq!(output, "out.svg".to_string());
            assert_eq!(ladders, None);
            assert_eq!(container_ids, None);
            assert_eq!(arrangement_id, Some("arrangement-2".to_string()));
            assert!(conditions.topology_aware);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_render_pool_gel_with_conditions() {
    let cmd = parse_shell_line(
        "render-pool-gel-svg - out.svg --arrangement arrangement-2 --agarose-pct 1.6 --buffer tbe --topology-aware false",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::RenderPoolGelSvg {
            arrangement_id,
            conditions,
            ..
        } => {
            assert_eq!(arrangement_id, Some("arrangement-2".to_string()));
            assert_eq!(conditions.agarose_percent, 1.6);
            assert_eq!(conditions.buffer_model.as_str(), "tbe");
            assert!(!conditions.topology_aware);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_arrange_serial_command() {
    let cmd = parse_shell_line(
        "arrange-serial container-1,container-2 --id arr-x --name test --ladders 100bp,1kb",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::CreateArrangementSerial {
            container_ids,
            arrangement_id,
            name,
            ladders,
        } => {
            assert_eq!(
                container_ids,
                vec!["container-1".to_string(), "container-2".to_string()]
            );
            assert_eq!(arrangement_id, Some("arr-x".to_string()));
            assert_eq!(name, Some("test".to_string()));
            assert_eq!(ladders, Some(vec!["100bp".to_string(), "1kb".to_string()]));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_arrange_set_ladders_command() {
    let cmd =
        parse_shell_line("arrange-set-ladders arr-x --ladders 100bp,1kb").expect("parse command");
    match cmd {
        ShellCommand::SetArrangementLadders {
            arrangement_id,
            ladders,
        } => {
            assert_eq!(arrangement_id, "arr-x".to_string());
            assert_eq!(ladders, Some(vec!["100bp".to_string(), "1kb".to_string()]));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_arrange_set_ladders_updates_existing_arrangement() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_a".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Lane A".to_string()),
            members: vec!["seq_a".to_string()],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-x".to_string(),
        Arrangement {
            arrangement_id: "arr-x".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Demo".to_string()),
            lane_container_ids: vec!["container-1".to_string()],
            ladders: vec!["Old ladder".to_string()],
            lane_role_labels: vec!["lane_1".to_string()],
            default_rack_id: None,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::SetArrangementLadders {
            arrangement_id: "arr-x".to_string(),
            ladders: Some(vec![
                "Plasmid Factory 1kb DNA Ladder".to_string(),
                "GeneRuler 100bp DNA Ladder Plus".to_string(),
            ]),
        },
    )
    .expect("execute arrange-set-ladders");

    assert!(out.state_changed);
    let arrangement = engine
        .state()
        .container_state
        .arrangements
        .get("arr-x")
        .expect("arrangement");
    assert_eq!(
        arrangement.ladders,
        vec![
            "Plasmid Factory 1kb DNA Ladder".to_string(),
            "GeneRuler 100bp DNA Ladder Plus".to_string()
        ]
    );
}

#[test]
fn parse_racks_create_from_arrangement_command() {
    let cmd = parse_shell_line(
        "racks create-from-arrangement arr-x --rack-id rack-x --name Bench --profile plate_96",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::RacksCreateFromArrangement {
            arrangement_id,
            rack_id,
            name,
            profile,
        } => {
            assert_eq!(arrangement_id, "arr-x".to_string());
            assert_eq!(rack_id, Some("rack-x".to_string()));
            assert_eq!(name, Some("Bench".to_string()));
            assert_eq!(profile, Some(RackProfileKind::Plate96));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_move_command() {
    let cmd =
        parse_shell_line("racks move rack-1 --from A1 --to B2 --block").expect("parse command");
    match cmd {
        ShellCommand::RacksMove {
            rack_id,
            from_coordinate,
            to_coordinate,
            move_block,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(from_coordinate, "A1".to_string());
            assert_eq!(to_coordinate, "B2".to_string());
            assert!(move_block);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_move_samples_command() {
    let cmd = parse_shell_line("racks move-samples rack-1 --from A1 --from A3 --to B2")
        .expect("parse command");
    match cmd {
        ShellCommand::RacksMoveSamples {
            rack_id,
            from_coordinates,
            to_coordinate,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(from_coordinates, vec!["A1".to_string(), "A3".to_string()]);
            assert_eq!(to_coordinate, "B2".to_string());
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_move_blocks_command() {
    let cmd = parse_shell_line(
        "racks move-blocks rack-1 --arrangement arr-a --arrangement arr-b --to B2",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::RacksMoveBlocks {
            rack_id,
            arrangement_ids,
            to_coordinate,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(
                arrangement_ids,
                vec!["arr-a".to_string(), "arr-b".to_string()]
            );
            assert_eq!(to_coordinate, "B2".to_string());
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_set_custom_profile_command() {
    let cmd = parse_shell_line("racks set-custom-profile rack-1 3 10").expect("parse command");
    match cmd {
        ShellCommand::RacksSetCustomProfile {
            rack_id,
            rows,
            columns,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(rows, 3);
            assert_eq!(columns, 10);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_set_fill_direction_command() {
    let cmd =
        parse_shell_line("racks set-fill-direction rack-1 column_major").expect("parse command");
    match cmd {
        ShellCommand::RacksSetFillDirection {
            rack_id,
            fill_direction,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(fill_direction, RackFillDirection::ColumnMajor);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_apply_template_command() {
    let cmd = parse_shell_line("racks apply-template rack-1 plate_edge_avoidance")
        .expect("parse command");
    match cmd {
        ShellCommand::RacksApplyTemplate { rack_id, template } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(template, RackAuthoringTemplate::PlateEdgeAvoidance);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_set_blocked_command() {
    let cmd = parse_shell_line("racks set-blocked rack-1 A1 B2,AA3").expect("parse command");
    match cmd {
        ShellCommand::RacksSetBlocked {
            rack_id,
            blocked_coordinates,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(
                blocked_coordinates,
                vec!["A1".to_string(), "B2".to_string(), "AA3".to_string()]
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_labels_svg_command_with_preset() {
    let cmd = parse_shell_line(
        "racks labels-svg rack-1 labels.svg --arrangement arr-x --preset print_a4",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::RacksLabelsSvg {
            rack_id,
            output,
            arrangement_id,
            preset,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(output, "labels.svg".to_string());
            assert_eq!(arrangement_id, Some("arr-x".to_string()));
            assert_eq!(preset, RackLabelSheetPreset::PrintA4);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_fabrication_svg_command() {
    let cmd = parse_shell_line(
        "racks fabrication-svg rack-1 rack.svg --template pipetting_pcr_tube_rack",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::RacksFabricationSvg {
            rack_id,
            output,
            template,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(output, "rack.svg".to_string());
            assert_eq!(template, RackPhysicalTemplateKind::PipettingPcrTubeRack);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_isometric_svg_command() {
    let cmd = parse_shell_line(
        "racks isometric-svg rack-1 rack.iso.svg --template storage_pcr_tube_rack",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::RacksIsometricSvg {
            rack_id,
            output,
            template,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(output, "rack.iso.svg".to_string());
            assert_eq!(template, RackPhysicalTemplateKind::StoragePcrTubeRack);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_openscad_command() {
    let cmd = parse_shell_line("racks openscad rack-1 rack.scad --template storage_pcr_tube_rack")
        .expect("parse command");
    match cmd {
        ShellCommand::RacksOpenScad {
            rack_id,
            output,
            template,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(output, "rack.scad".to_string());
            assert_eq!(template, RackPhysicalTemplateKind::StoragePcrTubeRack);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_carrier_labels_svg_command() {
    let cmd = parse_shell_line(
        "racks carrier-labels-svg rack-1 carrier.svg --arrangement arr-x --template pipetting_pcr_tube_rack --preset front_strip_only",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::RacksCarrierLabelsSvg {
            rack_id,
            output,
            arrangement_id,
            template,
            preset,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(output, "carrier.svg".to_string());
            assert_eq!(arrangement_id, Some("arr-x".to_string()));
            assert_eq!(template, RackPhysicalTemplateKind::PipettingPcrTubeRack);
            assert_eq!(preset, RackCarrierLabelPreset::FrontStripOnly);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_racks_simulation_json_command() {
    let cmd =
        parse_shell_line("racks simulation-json rack-1 rack.json --template storage_pcr_tube_rack")
            .expect("parse command");
    match cmd {
        ShellCommand::RacksSimulationJson {
            rack_id,
            output,
            template,
        } => {
            assert_eq!(rack_id, "rack-1".to_string());
            assert_eq!(output, "rack.json".to_string());
            assert_eq!(template, RackPhysicalTemplateKind::StoragePcrTubeRack);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_racks_move_samples_updates_snapshot() {
    let mut state = ProjectState::default();
    for (idx, seq_id) in ["seq_a", "seq_b", "seq_c", "seq_d"].iter().enumerate() {
        state.sequences.insert(
            (*seq_id).to_string(),
            DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
        );
        state.container_state.containers.insert(
            format!("container-{}", idx + 1),
            Container {
                container_id: format!("container-{}", idx + 1),
                kind: ContainerKind::Singleton,
                name: Some(format!("Lane {}", idx + 1)),
                members: vec![(*seq_id).to_string()],
                declared_contents_exclusive: true,
                created_by_op: None,
                created_at_unix_ms: 0,
            },
        );
    }
    state.container_state.arrangements.insert(
        "arr-a".to_string(),
        Arrangement {
            arrangement_id: "arr-a".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("A".to_string()),
            lane_container_ids: vec![
                "container-1".to_string(),
                "container-2".to_string(),
                "container-3".to_string(),
                "container-4".to_string(),
            ],
            ladders: vec![],
            lane_role_labels: vec![
                "lane_1".to_string(),
                "lane_2".to_string(),
                "lane_3".to_string(),
                "lane_4".to_string(),
            ],
            default_rack_id: Some("rack-1".to_string()),
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.racks.insert(
        "rack-1".to_string(),
        Rack {
            rack_id: "rack-1".to_string(),
            name: "Bench".to_string(),
            profile: RackProfileSnapshot::custom(1, 4),
            placements: vec![
                RackPlacementEntry {
                    coordinate: "A1".to_string(),
                    occupant: Some(RackOccupant::Container {
                        container_id: "container-1".to_string(),
                    }),
                    arrangement_id: "arr-a".to_string(),
                    order_index: 0,
                    role_label: "lane_1".to_string(),
                },
                RackPlacementEntry {
                    coordinate: "A2".to_string(),
                    occupant: Some(RackOccupant::Container {
                        container_id: "container-2".to_string(),
                    }),
                    arrangement_id: "arr-a".to_string(),
                    order_index: 1,
                    role_label: "lane_2".to_string(),
                },
                RackPlacementEntry {
                    coordinate: "A3".to_string(),
                    occupant: Some(RackOccupant::Container {
                        container_id: "container-3".to_string(),
                    }),
                    arrangement_id: "arr-a".to_string(),
                    order_index: 2,
                    role_label: "lane_3".to_string(),
                },
                RackPlacementEntry {
                    coordinate: "A4".to_string(),
                    occupant: Some(RackOccupant::Container {
                        container_id: "container-4".to_string(),
                    }),
                    arrangement_id: "arr-a".to_string(),
                    order_index: 3,
                    role_label: "lane_4".to_string(),
                },
            ],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RacksMoveSamples {
            rack_id: "rack-1".to_string(),
            from_coordinates: vec!["A1".to_string(), "A3".to_string()],
            to_coordinate: "A2".to_string(),
        },
    )
    .expect("move samples");

    assert!(out.state_changed);
    let rack = engine
        .state()
        .container_state
        .racks
        .get("rack-1")
        .expect("rack");
    assert_eq!(
        rack.placements
            .iter()
            .filter_map(|entry| match entry.occupant.as_ref() {
                Some(RackOccupant::Container { container_id }) => Some(container_id.as_str()),
                _ => None,
            })
            .collect::<Vec<_>>(),
        vec!["container-2", "container-1", "container-3", "container-4"]
    );
}

#[test]
fn execute_racks_set_fill_direction_updates_snapshot() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_a".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Lane A".to_string()),
            members: vec!["seq_a".to_string()],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-x".to_string(),
        Arrangement {
            arrangement_id: "arr-x".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Demo".to_string()),
            lane_container_ids: vec!["container-1".to_string()],
            ladders: vec![],
            lane_role_labels: vec!["vector".to_string()],
            default_rack_id: Some("rack-1".to_string()),
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.racks.insert(
        "rack-1".to_string(),
        Rack {
            rack_id: "rack-1".to_string(),
            name: "Bench".to_string(),
            profile: RackProfileSnapshot::from_kind(RackProfileKind::SmallTube4x6),
            placements: vec![RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: Some(RackOccupant::Container {
                    container_id: "container-1".to_string(),
                }),
                arrangement_id: "arr-x".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            }],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RacksSetFillDirection {
            rack_id: "rack-1".to_string(),
            fill_direction: RackFillDirection::ColumnMajor,
        },
    )
    .expect("set fill direction");
    assert!(out.state_changed);
    let rack = engine
        .state()
        .container_state
        .racks
        .get("rack-1")
        .expect("rack");
    assert_eq!(rack.profile.fill_direction, RackFillDirection::ColumnMajor);
}

#[test]
fn execute_racks_move_blocks_updates_snapshot() {
    let mut state = ProjectState::default();
    for (idx, seq_id) in ["seq_a", "seq_b", "seq_c", "seq_d"].iter().enumerate() {
        state.sequences.insert(
            (*seq_id).to_string(),
            DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
        );
        state.container_state.containers.insert(
            format!("container-{}", idx + 1),
            Container {
                container_id: format!("container-{}", idx + 1),
                kind: ContainerKind::Singleton,
                name: Some(format!("Lane {}", idx + 1)),
                members: vec![(*seq_id).to_string()],
                declared_contents_exclusive: true,
                created_by_op: None,
                created_at_unix_ms: 0,
            },
        );
    }
    state.container_state.arrangements.insert(
        "arr-a".to_string(),
        Arrangement {
            arrangement_id: "arr-a".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("A".to_string()),
            lane_container_ids: vec!["container-1".to_string(), "container-2".to_string()],
            ladders: vec![],
            lane_role_labels: vec!["lane_1".to_string(), "lane_2".to_string()],
            default_rack_id: Some("rack-1".to_string()),
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-b".to_string(),
        Arrangement {
            arrangement_id: "arr-b".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("B".to_string()),
            lane_container_ids: vec!["container-3".to_string()],
            ladders: vec![],
            lane_role_labels: vec!["lane_1".to_string()],
            default_rack_id: None,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-c".to_string(),
        Arrangement {
            arrangement_id: "arr-c".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("C".to_string()),
            lane_container_ids: vec!["container-4".to_string()],
            ladders: vec![],
            lane_role_labels: vec!["lane_1".to_string()],
            default_rack_id: None,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.racks.insert(
        "rack-1".to_string(),
        Rack {
            rack_id: "rack-1".to_string(),
            name: "Bench".to_string(),
            profile: RackProfileSnapshot::custom(1, 4),
            placements: vec![
                RackPlacementEntry {
                    coordinate: "A1".to_string(),
                    occupant: Some(RackOccupant::Container {
                        container_id: "container-1".to_string(),
                    }),
                    arrangement_id: "arr-a".to_string(),
                    order_index: 0,
                    role_label: "lane_1".to_string(),
                },
                RackPlacementEntry {
                    coordinate: "A2".to_string(),
                    occupant: Some(RackOccupant::Container {
                        container_id: "container-2".to_string(),
                    }),
                    arrangement_id: "arr-a".to_string(),
                    order_index: 1,
                    role_label: "lane_2".to_string(),
                },
                RackPlacementEntry {
                    coordinate: "A3".to_string(),
                    occupant: Some(RackOccupant::Container {
                        container_id: "container-3".to_string(),
                    }),
                    arrangement_id: "arr-b".to_string(),
                    order_index: 0,
                    role_label: "lane_1".to_string(),
                },
                RackPlacementEntry {
                    coordinate: "A4".to_string(),
                    occupant: Some(RackOccupant::Container {
                        container_id: "container-4".to_string(),
                    }),
                    arrangement_id: "arr-c".to_string(),
                    order_index: 0,
                    role_label: "lane_1".to_string(),
                },
            ],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RacksMoveBlocks {
            rack_id: "rack-1".to_string(),
            arrangement_ids: vec!["arr-c".to_string(), "arr-a".to_string()],
            to_coordinate: "A3".to_string(),
        },
    )
    .expect("move blocks");

    assert!(out.state_changed);
    let rack = engine
        .state()
        .container_state
        .racks
        .get("rack-1")
        .expect("rack");
    assert_eq!(
        rack.placements
            .iter()
            .map(|entry| entry.arrangement_id.as_str())
            .collect::<Vec<_>>(),
        vec!["arr-a", "arr-a", "arr-c", "arr-b"]
    );
}

#[test]
fn execute_racks_apply_template_updates_snapshot() {
    execute_racks_apply_template_updates_snapshot_impl();
}

#[inline(never)]
fn execute_racks_apply_template_updates_snapshot_impl() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_a".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Lane A".to_string()),
            members: vec!["seq_a".to_string()],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-x".to_string(),
        Arrangement {
            arrangement_id: "arr-x".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Demo".to_string()),
            lane_container_ids: vec!["container-1".to_string()],
            ladders: vec![],
            lane_role_labels: vec!["vector".to_string()],
            default_rack_id: Some("rack-1".to_string()),
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.racks.insert(
        "rack-1".to_string(),
        Rack {
            rack_id: "rack-1".to_string(),
            name: "Bench".to_string(),
            profile: RackProfileSnapshot::custom(4, 4),
            placements: vec![RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: Some(RackOccupant::Container {
                    container_id: "container-1".to_string(),
                }),
                arrangement_id: "arr-x".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            }],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RacksApplyTemplate {
            rack_id: "rack-1".to_string(),
            template: RackAuthoringTemplate::PlateEdgeAvoidance,
        },
    )
    .expect("apply rack template");
    assert!(out.state_changed);
    let rack = engine
        .state()
        .container_state
        .racks
        .get("rack-1")
        .expect("rack");
    assert_eq!(rack.profile.fill_direction, RackFillDirection::ColumnMajor);
    assert_eq!(
        rack.profile.blocked_coordinates.first().map(String::as_str),
        Some("A1")
    );
    assert!(rack.profile.blocked_coordinates.contains(&"D4".to_string()));
}

#[test]
fn execute_racks_set_custom_profile_updates_snapshot() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_a".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Lane A".to_string()),
            members: vec!["seq_a".to_string()],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-x".to_string(),
        Arrangement {
            arrangement_id: "arr-x".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Demo".to_string()),
            lane_container_ids: vec!["container-1".to_string()],
            ladders: vec![],
            lane_role_labels: vec!["vector".to_string()],
            default_rack_id: Some("rack-1".to_string()),
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.racks.insert(
        "rack-1".to_string(),
        Rack {
            rack_id: "rack-1".to_string(),
            name: "Bench".to_string(),
            profile: RackProfileSnapshot::from_kind(RackProfileKind::SmallTube4x6),
            placements: vec![RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: Some(RackOccupant::Container {
                    container_id: "container-1".to_string(),
                }),
                arrangement_id: "arr-x".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            }],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RacksSetCustomProfile {
            rack_id: "rack-1".to_string(),
            rows: 3,
            columns: 10,
        },
    )
    .expect("set custom profile");
    assert!(out.state_changed);
    let rack = engine
        .state()
        .container_state
        .racks
        .get("rack-1")
        .expect("rack");
    assert_eq!(rack.profile.kind, RackProfileKind::Custom);
    assert_eq!(rack.profile.rows, 3);
    assert_eq!(rack.profile.columns, 10);
}

#[test]
fn execute_racks_set_blocked_updates_snapshot() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_a".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Lane A".to_string()),
            members: vec!["seq_a".to_string()],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-x".to_string(),
        Arrangement {
            arrangement_id: "arr-x".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Demo".to_string()),
            lane_container_ids: vec!["container-1".to_string()],
            ladders: vec![],
            lane_role_labels: vec!["vector".to_string()],
            default_rack_id: Some("rack-1".to_string()),
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.racks.insert(
        "rack-1".to_string(),
        Rack {
            rack_id: "rack-1".to_string(),
            name: "Bench".to_string(),
            profile: RackProfileSnapshot::custom(28, 3),
            placements: vec![RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: Some(RackOccupant::Container {
                    container_id: "container-1".to_string(),
                }),
                arrangement_id: "arr-x".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            }],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RacksSetBlocked {
            rack_id: "rack-1".to_string(),
            blocked_coordinates: vec!["B2".to_string(), "AA3".to_string()],
        },
    )
    .expect("set blocked");
    assert!(out.state_changed);
    let rack = engine
        .state()
        .container_state
        .racks
        .get("rack-1")
        .expect("rack");
    assert_eq!(
        rack.profile.blocked_coordinates,
        vec!["B2".to_string(), "AA3".to_string()]
    );
}

#[test]
fn execute_racks_show_returns_structured_rack_state() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_a".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Lane A".to_string()),
            members: vec!["seq_a".to_string()],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-x".to_string(),
        Arrangement {
            arrangement_id: "arr-x".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Demo".to_string()),
            lane_container_ids: vec!["container-1".to_string()],
            ladders: vec!["NEB 1kb DNA Ladder".to_string()],
            lane_role_labels: vec!["vector".to_string()],
            default_rack_id: Some("rack-1".to_string()),
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.racks.insert(
        "rack-1".to_string(),
        Rack {
            rack_id: "rack-1".to_string(),
            name: "Bench".to_string(),
            profile: RackProfileSnapshot::from_kind(RackProfileKind::SmallTube4x6),
            placements: vec![RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: Some(RackOccupant::Container {
                    container_id: "container-1".to_string(),
                }),
                arrangement_id: "arr-x".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            }],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RacksShow {
            rack_id: "rack-1".to_string(),
        },
    )
    .expect("show rack");
    assert!(!out.state_changed);
    assert_eq!(out.output["schema"], "gentle.rack_state.v1");
    assert_eq!(out.output["rack"]["rack_id"], "rack-1");
    assert_eq!(out.output["placements"][0]["coordinate"], "A1");
    assert_eq!(out.output["placements"][0]["occupant"]["kind"], "container");
}

#[test]
fn execute_racks_labels_svg_with_preset_writes_marker() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_a".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Lane A".to_string()),
            members: vec!["seq_a".to_string()],
            declared_contents_exclusive: true,
            created_by_op: Some("op-1".to_string()),
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-x".to_string(),
        Arrangement {
            arrangement_id: "arr-x".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Demo".to_string()),
            lane_container_ids: vec!["container-1".to_string()],
            ladders: vec![],
            lane_role_labels: vec!["vector".to_string()],
            default_rack_id: Some("rack-1".to_string()),
            created_by_op: Some("op-2".to_string()),
            created_at_unix_ms: 0,
        },
    );
    state.container_state.racks.insert(
        "rack-1".to_string(),
        Rack {
            rack_id: "rack-1".to_string(),
            name: "Bench".to_string(),
            profile: RackProfileSnapshot::from_kind(RackProfileKind::SmallTube4x6),
            placements: vec![RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: Some(RackOccupant::Container {
                    container_id: "container-1".to_string(),
                }),
                arrangement_id: "arr-x".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            }],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempfile::NamedTempFile::new().expect("temp");
    let output = tmp.path().with_extension("labels.svg");
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RacksLabelsSvg {
            rack_id: "rack-1".to_string(),
            output: output.display().to_string(),
            arrangement_id: Some("arr-x".to_string()),
            preset: RackLabelSheetPreset::WideCards,
        },
    )
    .expect("export labels");
    assert!(!out.state_changed);
    let svg = fs::read_to_string(&output).expect("labels svg");
    assert!(svg.contains("data-label-preset=\"wide_cards\""));
    assert!(svg.contains("viewBox=\"0 0 556 156\""));
    assert!(svg.contains("role vector"));
    assert!(svg.contains("arrangement Demo"));
}

#[test]
fn parse_render_rna_svg() {
    let cmd = parse_shell_line("render-rna-svg rna_seq rna.svg").expect("parse command");
    match cmd {
        ShellCommand::RenderRnaSvg { seq_id, output } => {
            assert_eq!(seq_id, "rna_seq".to_string());
            assert_eq!(output, "rna.svg".to_string());
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_racks_physical_exports_write_markers() {
    std::thread::Builder::new()
        .name("rack-physical-export-shell-test".to_string())
        .stack_size(32 * 1024 * 1024)
        .spawn(|| {
            let mut state = ProjectState::default();
            state.sequences.insert(
                "seq_a".to_string(),
                DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
            );
            state.container_state.containers.insert(
                "container-1".to_string(),
                Container {
                    container_id: "container-1".to_string(),
                    kind: ContainerKind::Singleton,
                    name: Some("Vector".to_string()),
                    members: vec!["seq_a".to_string()],
                    declared_contents_exclusive: true,
                    created_by_op: None,
                    created_at_unix_ms: 0,
                },
            );
            state.container_state.arrangements.insert(
                "arr-x".to_string(),
                Arrangement {
                    arrangement_id: "arr-x".to_string(),
                    mode: ArrangementMode::Serial,
                    name: Some("Demo".to_string()),
                    lane_container_ids: vec!["container-1".to_string()],
                    ladders: vec!["1 kb Ladder".to_string()],
                    lane_role_labels: vec!["vector".to_string()],
                    default_rack_id: Some("rack-1".to_string()),
                    created_by_op: None,
                    created_at_unix_ms: 0,
                },
            );
            state.container_state.racks.insert(
                "rack-1".to_string(),
                Rack {
                    rack_id: "rack-1".to_string(),
                    name: "Bench".to_string(),
                    profile: RackProfileSnapshot::from_kind(RackProfileKind::SmallTube4x6),
                    placements: vec![
                        RackPlacementEntry {
                            coordinate: "A1".to_string(),
                            occupant: Some(RackOccupant::LadderReference {
                                ladder_name: "1 kb Ladder".to_string(),
                            }),
                            arrangement_id: "arr-x".to_string(),
                            order_index: 0,
                            role_label: "ladder_left".to_string(),
                        },
                        RackPlacementEntry {
                            coordinate: "A2".to_string(),
                            occupant: Some(RackOccupant::Container {
                                container_id: "container-1".to_string(),
                            }),
                            arrangement_id: "arr-x".to_string(),
                            order_index: 1,
                            role_label: "vector".to_string(),
                        },
                    ],
                    created_by_op: None,
                    created_at_unix_ms: 0,
                },
            );
            let mut engine = GentleEngine::from_state(state);
            let temp = tempdir().expect("tempdir");
            let svg_path = temp.path().join("rack.fabrication.svg");
            let isometric_path = temp.path().join("rack.isometric.svg");
            let scad_path = temp.path().join("rack.scad");
            let carrier_path = temp.path().join("rack.carrier.svg");
            let simulation_path = temp.path().join("rack.simulation.json");

            let svg_result = execute_shell_command(
                &mut engine,
                &ShellCommand::RacksFabricationSvg {
                    rack_id: "rack-1".to_string(),
                    output: svg_path.display().to_string(),
                    template: RackPhysicalTemplateKind::StoragePcrTubeRack,
                },
            )
            .expect("fabrication export");
            assert!(!svg_result.state_changed);
            let svg = fs::read_to_string(&svg_path).expect("fabrication svg");
            assert!(svg.contains("data-rack-physical-template=\"storage_pcr_tube_rack\""));

            let isometric_result = execute_shell_command(
                &mut engine,
                &ShellCommand::RacksIsometricSvg {
                    rack_id: "rack-1".to_string(),
                    output: isometric_path.display().to_string(),
                    template: RackPhysicalTemplateKind::StoragePcrTubeRack,
                },
            )
            .expect("isometric export");
            assert!(!isometric_result.state_changed);
            let isometric = fs::read_to_string(&isometric_path).expect("rack isometric svg");
            assert!(isometric.contains("data-rack-isometric-template=\"storage_pcr_tube_rack\""));
            assert!(isometric.contains("GENtle rack isometric sketch"));

            let scad_result = execute_shell_command(
                &mut engine,
                &ShellCommand::RacksOpenScad {
                    rack_id: "rack-1".to_string(),
                    output: scad_path.display().to_string(),
                    template: RackPhysicalTemplateKind::PipettingPcrTubeRack,
                },
            )
            .expect("openscad export");
            assert!(!scad_result.state_changed);
            let scad = fs::read_to_string(&scad_path).expect("rack scad");
            assert!(scad.contains("template=pipetting_pcr_tube_rack"));

            let carrier_result = execute_shell_command(
                &mut engine,
                &ShellCommand::RacksCarrierLabelsSvg {
                    rack_id: "rack-1".to_string(),
                    output: carrier_path.display().to_string(),
                    arrangement_id: Some("arr-x".to_string()),
                    template: RackPhysicalTemplateKind::StoragePcrTubeRack,
                    preset: RackCarrierLabelPreset::ModuleCardsOnly,
                },
            )
            .expect("carrier labels export");
            assert!(!carrier_result.state_changed);
            let carrier = fs::read_to_string(&carrier_path).expect("carrier svg");
            assert!(carrier.contains("data-rack-carrier-template=\"storage_pcr_tube_rack\""));
            assert!(carrier.contains("data-rack-carrier-preset=\"module_cards_only\""));

            let simulation_result = execute_shell_command(
                &mut engine,
                &ShellCommand::RacksSimulationJson {
                    rack_id: "rack-1".to_string(),
                    output: simulation_path.display().to_string(),
                    template: RackPhysicalTemplateKind::PipettingPcrTubeRack,
                },
            )
            .expect("simulation export");
            assert!(!simulation_result.state_changed);
            let simulation = fs::read_to_string(&simulation_path).expect("simulation json");
            assert!(simulation.contains("\"schema\": \"gentle.rack_simulation_export.v1\""));
            assert!(simulation.contains("\"kind\": \"pipetting_pcr_tube_rack\""));
        })
        .expect("spawn shell test thread")
        .join()
        .expect("join shell test thread");
}

#[test]
fn parse_protocol_cartoon_list() {
    let cmd = parse_shell_line("protocol-cartoon list").expect("parse command");
    match cmd {
        ShellCommand::ProtocolCartoonList => {}
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_protocol_cartoon_render_svg() {
    let cmd = parse_shell_line("protocol-cartoon render-svg gibson.two_fragment out.svg")
        .expect("parse command");
    match cmd {
        ShellCommand::RenderProtocolCartoonSvg { protocol, output } => {
            assert_eq!(protocol, ProtocolCartoonKind::GibsonTwoFragment);
            assert_eq!(output, "out.svg");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn reject_render_protocol_cartoon_svg_alias() {
    assert!(parse_shell_line("render-protocol-cartoon-svg gibson out.svg").is_err());
}

#[test]
fn parse_protocol_cartoon_render_template_svg() {
    let cmd = parse_shell_line("protocol-cartoon render-template-svg template.json out.svg")
        .expect("parse command");
    match cmd {
        ShellCommand::RenderProtocolCartoonTemplateSvg {
            template_path,
            output,
        } => {
            assert_eq!(template_path, "template.json");
            assert_eq!(output, "out.svg");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn reject_render_protocol_cartoon_template_svg_alias() {
    assert!(
        parse_shell_line("render-protocol-cartoon-template-svg template.json out.svg").is_err()
    );
}

#[test]
fn parse_protocol_cartoon_template_export() {
    let cmd = parse_shell_line("protocol-cartoon template-export gibson out.json")
        .expect("parse command");
    match cmd {
        ShellCommand::ExportProtocolCartoonTemplateJson { protocol, output } => {
            assert_eq!(protocol, ProtocolCartoonKind::GibsonTwoFragment);
            assert_eq!(output, "out.json");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_protocol_cartoon_template_validate() {
    let cmd = parse_shell_line("protocol-cartoon template-validate template.json")
        .expect("parse command");
    match cmd {
        ShellCommand::ValidateProtocolCartoonTemplate { template_path } => {
            assert_eq!(template_path, "template.json");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_protocol_cartoon_render_with_bindings() {
    let cmd = parse_shell_line(
        "protocol-cartoon render-with-bindings template.json bindings.json out.svg",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::RenderProtocolCartoonTemplateWithBindingsSvg {
            template_path,
            bindings_path,
            output,
        } => {
            assert_eq!(template_path, "template.json");
            assert_eq!(bindings_path, "bindings.json");
            assert_eq!(output, "out.svg");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_gibson_preview() {
    let cmd =
        parse_shell_line("gibson preview @plan.json --output preview.json").expect("parse command");
    match cmd {
        ShellCommand::GibsonPreview {
            request_json,
            output_path,
        } => {
            assert_eq!(request_json, "@plan.json");
            assert_eq!(output_path.as_deref(), Some("preview.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_gibson_apply() {
    let cmd = parse_shell_line("gibson apply @plan.json").expect("parse command");
    match cmd {
        ShellCommand::GibsonApply { request_json } => {
            assert_eq!(request_json, "@plan.json");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_rna_info() {
    let cmd = parse_shell_line("rna-info rna_seq").expect("parse command");
    match cmd {
        ShellCommand::RnaInfo { seq_id } => {
            assert_eq!(seq_id, "rna_seq".to_string());
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_protocol_cartoon_list() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(&mut engine, &ShellCommand::ProtocolCartoonList)
        .expect("execute protocol-cartoon list");
    assert!(!out.state_changed);
    let rows = out.output.as_array().expect("catalog rows array");
    assert!(!rows.is_empty());
    assert_eq!(rows[0]["id"].as_str(), Some("gibson.two_fragment"));
}

#[test]
fn execute_render_protocol_cartoon_svg() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let tmp = tempdir().expect("temp dir");
    let output = tmp.path().join("gibson.svg");
    let output_path = output.display().to_string();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RenderProtocolCartoonSvg {
            protocol: ProtocolCartoonKind::GibsonTwoFragment,
            output: output_path.clone(),
        },
    )
    .expect("execute render protocol cartoon");
    assert!(!out.state_changed);
    let svg = fs::read_to_string(output).expect("read protocol cartoon svg");
    assert!(svg.contains("<svg"));
    assert!(svg.contains("Chew-back"));
}

#[test]
fn execute_render_protocol_cartoon_template_svg() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let tmp = tempdir().expect("temp dir");
    let template = tmp.path().join("template.json");
    fs::write(
        &template,
        r#"{
            "id":"demo.protocol",
            "events":[
                {"id":"step_a","title":"Step A"}
            ]
        }"#,
    )
    .expect("write template json");
    let output = tmp.path().join("template.svg");
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RenderProtocolCartoonTemplateSvg {
            template_path: template.display().to_string(),
            output: output.display().to_string(),
        },
    )
    .expect("execute render protocol cartoon template");
    assert!(!out.state_changed);
    let svg = fs::read_to_string(output).expect("read protocol cartoon template svg");
    assert!(svg.contains("<svg"));
    assert!(svg.contains("demo.protocol"));
}

#[test]
fn execute_protocol_cartoon_template_validate() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let tmp = tempdir().expect("temp dir");
    let template = tmp.path().join("template.validate.json");
    fs::write(
        &template,
        r#"{
            "id":"demo.protocol",
            "events":[
                {"id":"step_a","title":"Step A"}
            ]
        }"#,
    )
    .expect("write template json");
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ValidateProtocolCartoonTemplate {
            template_path: template.display().to_string(),
        },
    )
    .expect("execute protocol cartoon template validate");
    assert!(!out.state_changed);
    let message = out.output["result"]["messages"][0]
        .as_str()
        .expect("validate message string");
    assert!(message.contains("Validated protocol cartoon template 'demo.protocol'"));
    assert!(message.contains("(events=1)"));
}

#[test]
fn execute_render_protocol_cartoon_template_with_bindings_svg() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let tmp = tempdir().expect("temp dir");
    let template = tmp.path().join("template.bindings.json");
    fs::write(
        &template,
        r#"{
            "id":"demo.protocol",
            "events":[
                {
                    "id":"step_a",
                    "title":"Step A",
                    "molecules":[
                        {
                            "id":"mol_a",
                            "features":[
                                {"id":"frag_a","label":"Fragment A","length_bp":42}
                            ]
                        }
                    ]
                }
            ]
        }"#,
    )
    .expect("write template json");
    let bindings = tmp.path().join("bindings.json");
    fs::write(
        &bindings,
        r#"{
            "template_id":"demo.protocol",
            "event_overrides":[
                {"event_id":"step_a","title":"Bound Step","caption":"Bound caption"}
            ]
        }"#,
    )
    .expect("write bindings json");
    let output = tmp.path().join("template.bound.svg");
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RenderProtocolCartoonTemplateWithBindingsSvg {
            template_path: template.display().to_string(),
            bindings_path: bindings.display().to_string(),
            output: output.display().to_string(),
        },
    )
    .expect("execute render protocol cartoon with bindings");
    assert!(!out.state_changed);
    let svg = fs::read_to_string(output).expect("read protocol cartoon template svg");
    assert!(svg.contains("<svg"));
    assert!(svg.contains("Bound Step"));
    assert!(svg.contains("Bound caption"));
}

#[test]
fn execute_protocol_cartoon_template_export() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let tmp = tempdir().expect("temp dir");
    let output = tmp.path().join("gibson.template.json");
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ExportProtocolCartoonTemplateJson {
            protocol: ProtocolCartoonKind::GibsonTwoFragment,
            output: output.display().to_string(),
        },
    )
    .expect("execute protocol cartoon template export");
    assert!(!out.state_changed);
    let json = fs::read_to_string(output).expect("read protocol template json");
    let template: serde_json::Value =
        serde_json::from_str(&json).expect("parse protocol template json");
    assert_eq!(
        template.get("schema").and_then(|v| v.as_str()),
        Some("gentle.protocol_cartoon_template.v1")
    );
    assert_eq!(
        template.get("id").and_then(|v| v.as_str()),
        Some("gibson.two_fragment")
    );
}

#[test]
fn execute_gibson_preview() {
    let mut state = ProjectState::default();
    let mut destination =
        DNAsequence::from_sequence("AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT")
            .expect("destination sequence");
    destination.set_circular(true);
    destination.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(2, 8),
        qualifiers: vec![("label".into(), Some("LEFT".to_string()))],
    });
    destination.features_mut().push(Feature {
        kind: "misc_feature".into(),
        location: Location::simple_range(10, 20),
        qualifiers: vec![
            (
                "label".into(),
                Some("Multiple Cloning Site (MCS)".to_string()),
            ),
            ("note".into(), Some("contains SmaI".to_string())),
        ],
    });
    state
        .sequences
        .insert("destination_vector".to_string(), destination);
    let mut insert = DNAsequence::from_sequence("ATGCGTACGTTAGCGTACGATCGTACGTAGCTAGCTAGCATCGATCGA")
        .expect("insert sequence");
    insert.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(4, 10),
        qualifiers: vec![("label".into(), Some("INSERT".to_string()))],
    });
    state
        .sequences
        .insert("insert_x_amplicon".to_string(), insert);
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempdir().expect("temp dir");
    let plan_path = tmp.path().join("gibson.plan.json");
    fs::write(
        &plan_path,
        r#"{
  "schema": "gentle.gibson_assembly_plan.v1",
  "id": "preview_test",
  "title": "Preview test",
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
  "product": {"topology": "circular", "output_id_hint": "out"},
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
}"#,
    )
    .expect("write plan json");
    let output_path = tmp.path().join("gibson.preview.json");
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::GibsonPreview {
            request_json: format!("@{}", plan_path.display()),
            output_path: Some(output_path.display().to_string()),
        },
    )
    .expect("execute gibson preview");
    assert!(!out.state_changed);
    assert_eq!(
        out.output.get("schema").and_then(|value| value.as_str()),
        Some("gentle.gibson_assembly_preview.v1")
    );
    assert_eq!(
        out.output
            .get("can_execute")
            .and_then(|value| value.as_bool()),
        Some(true)
    );
    assert_eq!(
        out.output
            .get("primer_suggestions")
            .and_then(|value| value.as_array())
            .map(Vec::len),
        Some(2)
    );
    let written: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(output_path).expect("read gibson preview output"))
            .expect("parse written preview json");
    assert_eq!(
        written.get("schema").and_then(|value| value.as_str()),
        Some("gentle.gibson_assembly_preview.v1")
    );
}

#[test]
fn execute_gibson_apply_creates_output_sequences() {
    let mut state = ProjectState::default();
    let mut destination =
        DNAsequence::from_sequence("AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT")
            .expect("destination sequence");
    destination.set_circular(true);
    destination.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(2, 8),
        qualifiers: vec![("label".into(), Some("LEFT".to_string()))],
    });
    destination.features_mut().push(Feature {
        kind: "misc_feature".into(),
        location: Location::simple_range(10, 20),
        qualifiers: vec![
            (
                "label".into(),
                Some("Multiple Cloning Site (MCS)".to_string()),
            ),
            ("note".into(), Some("contains SmaI".to_string())),
        ],
    });
    state
        .sequences
        .insert("destination_vector".to_string(), destination);
    let mut insert = DNAsequence::from_sequence("ATGCGTACGTTAGCGTACGATCGTACGTAGCTAGCTAGCATCGATCGA")
        .expect("insert sequence");
    insert.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(4, 10),
        qualifiers: vec![("label".into(), Some("INSERT".to_string()))],
    });
    state
        .sequences
        .insert("insert_x_amplicon".to_string(), insert);
    let mut engine = GentleEngine::from_state(state);
    let plan = r#"{
  "schema": "gentle.gibson_assembly_plan.v1",
  "id": "apply_test",
  "title": "Apply test",
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
  "product": {"topology": "circular", "output_id_hint": "out"},
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
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::GibsonApply {
            request_json: plan.to_string(),
        },
    )
    .expect("execute gibson apply");
    assert!(out.state_changed);
    let created = out.output["result"]["created_seq_ids"]
        .as_array()
        .expect("created ids array");
    assert_eq!(created.len(), 3);
    assert!(engine.state().sequences.contains_key("out"));
    let product = engine
        .state()
        .sequences
        .get("out")
        .expect("assembled product output");
    let labels = product
        .features()
        .iter()
        .filter_map(|feature| feature.qualifier_values("label").next().map(str::to_string))
        .collect::<Vec<_>>();
    let feature_debug = product
        .features()
        .iter()
        .map(|feature| {
            format!(
                "{} {:?}",
                feature.kind,
                feature.qualifier_values("label").next().map(str::to_string)
            )
        })
        .collect::<Vec<_>>();
    assert!(
        labels.iter().any(|label| label == "LEFT"),
        "product labels: {labels:?}; features: {feature_debug:?}"
    );
    assert!(
        labels.iter().any(|label| label == "INSERT"),
        "product labels: {labels:?}; features: {feature_debug:?}"
    );
    assert!(
        labels
            .iter()
            .any(|label| label == "Multiple Cloning Site (MCS)"),
        "product labels: {labels:?}; features: {feature_debug:?}"
    );
    let mcs = product
        .features()
        .iter()
        .find(|feature| {
            feature.qualifier_values("label").next() == Some("Multiple Cloning Site (MCS)")
        })
        .expect("assembled product MCS feature");
    assert_eq!(
        mcs.qualifier_values("mcs_crosscheck_status").next(),
        Some("validated_against_assembled_product")
    );
    assert!(
        engine
            .operation_log()
            .iter()
            .any(|record| matches!(record.op, Operation::ApplyGibsonAssemblyPlan { .. }))
    );
    let gibson_op_id = engine
        .operation_log()
        .last()
        .map(|record| record.result.op_id.clone())
        .expect("gibson op id");
    let created_by_gibson = engine
        .state()
        .container_state
        .containers
        .values()
        .filter(|container| container.created_by_op.as_deref() == Some(gibson_op_id.as_str()))
        .collect::<Vec<_>>();
    assert_eq!(created_by_gibson.len(), 3);
    assert!(created_by_gibson.iter().all(|container| {
        matches!(container.kind, ContainerKind::Singleton) && container.members.len() == 1
    }));
    let latest = &engine.state().container_state.seq_to_latest_container;
    for seq_id in [
        "insert_x_left_insert_primer",
        "insert_x_right_insert_primer",
        "out",
    ] {
        let container_id = latest
            .get(seq_id)
            .unwrap_or_else(|| panic!("missing latest container for {seq_id}"));
        let container = engine
            .state()
            .container_state
            .containers
            .get(container_id)
            .unwrap_or_else(|| panic!("missing container {container_id}"));
        assert!(matches!(container.kind, ContainerKind::Singleton));
        assert_eq!(container.members, vec![seq_id.to_string()]);
    }
    let arrangements = engine
        .state()
        .container_state
        .arrangements
        .values()
        .collect::<Vec<_>>();
    assert_eq!(arrangements.len(), 1);
    assert_eq!(arrangements[0].lane_container_ids.len(), 3);
    assert!(!arrangements[0].ladders.is_empty());
}

#[test]
fn parse_ladders_list_with_filter() {
    let cmd = parse_shell_line("ladders list --filter NEB").expect("parse command");
    match cmd {
        ShellCommand::LaddersList {
            molecule,
            name_filter,
        } => {
            assert_eq!(molecule, LadderMolecule::Dna);
            assert_eq!(name_filter, Some("NEB".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_ladders_export_with_filter() {
    let cmd =
        parse_shell_line("ladders export ladders.json --filter ruler").expect("parse command");
    match cmd {
        ShellCommand::LaddersExport {
            molecule,
            output,
            name_filter,
        } => {
            assert_eq!(molecule, LadderMolecule::Dna);
            assert_eq!(output, "ladders.json".to_string());
            assert_eq!(name_filter, Some("ruler".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_export_run_bundle_with_run_id() {
    let cmd = parse_shell_line("export-run-bundle run_bundle.json --run-id demo_run")
        .expect("parse command");
    match cmd {
        ShellCommand::ExportRunBundle { output, run_id } => {
            assert_eq!(output, "run_bundle.json");
            assert_eq!(run_id.as_deref(), Some("demo_run"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_ladders_list_rna() {
    let cmd = parse_shell_line("ladders list --molecule rna --filter ss").expect("parse");
    match cmd {
        ShellCommand::LaddersList {
            molecule,
            name_filter,
        } => {
            assert_eq!(molecule, LadderMolecule::Rna);
            assert_eq!(name_filter, Some("ss".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_proteases_list_with_filter_and_output() {
    let cmd = parse_shell_line("proteases list --filter trypsin --output proteases.json")
        .expect("parse proteases list");
    match cmd {
        ShellCommand::ProteasesList { filter, output } => {
            assert_eq!(filter.as_deref(), Some("trypsin"));
            assert_eq!(output.as_deref(), Some("proteases.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_proteases_show_with_output() {
    let cmd = parse_shell_line("proteases show Trypsin --output trypsin.json")
        .expect("parse proteases show");
    match cmd {
        ShellCommand::ProteasesShow { query, output } => {
            assert_eq!(query, "Trypsin");
            assert_eq!(output.as_deref(), Some("trypsin.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_reference_genes_with_regex_and_biotypes() {
    let cmd = parse_shell_line(
            "helpers genes Helper --filter '^bla$' --biotype promoter --biotype cds --limit 10 --offset 3",
        )
        .expect("parse command");
    match cmd {
        ShellCommand::ReferenceGenes {
            helper_mode,
            genome_id,
            filter,
            biotypes,
            limit,
            offset,
            ..
        } => {
            assert!(helper_mode);
            assert_eq!(genome_id, "Helper");
            assert_eq!(filter, "^bla$");
            assert_eq!(biotypes, vec!["promoter".to_string(), "cds".to_string()]);
            assert_eq!(limit, 10);
            assert_eq!(offset, 3);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_blast_with_options() {
    let cmd = parse_shell_line(
            "genomes blast ToyGenome ACGTACGT --max-hits 12 --task blastn --catalog c.json --cache-dir cache",
        )
        .expect("parse command");
    match cmd {
        ShellCommand::ReferenceBlast {
            helper_mode,
            genome_id,
            query_sequence,
            max_hits,
            max_hits_explicit,
            task,
            request_options_json,
            catalog_path,
            cache_dir,
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome");
            assert_eq!(query_sequence, "ACGTACGT");
            assert_eq!(max_hits, 12);
            assert!(max_hits_explicit);
            assert_eq!(task.as_deref(), Some("blastn"));
            assert!(request_options_json.is_none());
            assert_eq!(catalog_path.as_deref(), Some("c.json"));
            assert_eq!(cache_dir.as_deref(), Some("cache"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_helpers_blast_defaults() {
    let cmd = parse_shell_line("helpers blast pUC19 ACGTAG").expect("parse command");
    match cmd {
        ShellCommand::ReferenceBlast {
            helper_mode,
            genome_id,
            query_sequence,
            max_hits,
            max_hits_explicit,
            task,
            request_options_json,
            ..
        } => {
            assert!(helper_mode);
            assert_eq!(genome_id, "pUC19");
            assert_eq!(query_sequence, "ACGTAG");
            assert_eq!(max_hits, 25);
            assert!(!max_hits_explicit);
            assert!(task.is_none());
            assert!(request_options_json.is_none());
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_blast_start_with_options() {
    let cmd = parse_shell_line(
            "genomes blast-start ToyGenome ACGTACGT --max-hits 12 --task blastn --catalog c.json --cache-dir cache",
        )
        .expect("parse command");
    match cmd {
        ShellCommand::ReferenceBlastAsyncStart {
            helper_mode,
            genome_id,
            query_sequence,
            max_hits,
            max_hits_explicit,
            task,
            request_options_json,
            catalog_path,
            cache_dir,
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome");
            assert_eq!(query_sequence, "ACGTACGT");
            assert_eq!(max_hits, 12);
            assert!(max_hits_explicit);
            assert_eq!(task.as_deref(), Some("blastn"));
            assert!(request_options_json.is_none());
            assert_eq!(catalog_path.as_deref(), Some("c.json"));
            assert_eq!(cache_dir.as_deref(), Some("cache"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_blast_status_and_cancel() {
    let status = parse_shell_line("genomes blast-status blast-job-42 --with-report")
        .expect("parse status command");
    match status {
        ShellCommand::ReferenceBlastAsyncStatus {
            helper_mode,
            job_id,
            include_report,
        } => {
            assert!(!helper_mode);
            assert_eq!(job_id, "blast-job-42");
            assert!(include_report);
        }
        other => panic!("unexpected command: {other:?}"),
    }
    let cancel =
        parse_shell_line("helpers blast-cancel blast-job-42").expect("parse cancel command");
    match cancel {
        ShellCommand::ReferenceBlastAsyncCancel {
            helper_mode,
            job_id,
        } => {
            assert!(helper_mode);
            assert_eq!(job_id, "blast-job-42");
        }
        other => panic!("unexpected command: {other:?}"),
    }
    let listed = parse_shell_line("helpers blast-list").expect("parse list command");
    match listed {
        ShellCommand::ReferenceBlastAsyncList { helper_mode } => {
            assert!(helper_mode);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_blast_with_options_json_override() {
    let cmd = parse_shell_line(
            "genomes blast ToyGenome ACGTACGT --options-json '{\"max_hits\":7,\"thresholds\":{\"min_identity_percent\":97.5}}'",
        )
        .expect("parse command");
    match cmd {
        ShellCommand::ReferenceBlast {
            helper_mode,
            genome_id,
            max_hits,
            max_hits_explicit,
            task,
            request_options_json,
            ..
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome");
            assert_eq!(max_hits, 25);
            assert!(!max_hits_explicit);
            assert!(task.is_none());
            let options = request_options_json.expect("options json");
            assert_eq!(options.get("max_hits").and_then(|v| v.as_u64()), Some(7));
            let min_ident = options
                .get("thresholds")
                .and_then(|v| v.get("min_identity_percent"))
                .and_then(|v| v.as_f64());
            assert_eq!(min_ident, Some(97.5));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_helpers_blast_with_options_file_override() {
    let td = tempdir().expect("tempdir");
    let options_path = td.path().join("blast_options.json");
    fs::write(
        &options_path,
        r#"{"max_hits":9,"thresholds":{"min_bit_score":55.0}}"#,
    )
    .expect("write options file");
    let cmd = parse_shell_line(&format!(
        "helpers blast pUC19 ACGTAG --options-file {}",
        options_path.to_string_lossy()
    ))
    .expect("parse command");
    match cmd {
        ShellCommand::ReferenceBlast {
            helper_mode,
            request_options_json,
            max_hits_explicit,
            ..
        } => {
            assert!(helper_mode);
            assert!(!max_hits_explicit);
            let options = request_options_json.expect("options json");
            assert_eq!(options.get("max_hits").and_then(|v| v.as_u64()), Some(9));
            let min_bit = options
                .get("thresholds")
                .and_then(|v| v.get("min_bit_score"))
                .and_then(|v| v.as_f64());
            assert_eq!(min_bit, Some(55.0));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_blast_track_with_options() {
    let cmd = parse_shell_line(
            "genomes blast-track ToyGenome ACGTACGT query_seq --max-hits 8 --task blastn --track-name Hits --clear-existing --catalog c.json --cache-dir cache",
        )
        .expect("parse command");
    match cmd {
        ShellCommand::ReferenceBlastTrack {
            helper_mode,
            genome_id,
            query_sequence,
            target_seq_id,
            max_hits,
            max_hits_explicit,
            task,
            request_options_json,
            track_name,
            clear_existing,
            catalog_path,
            cache_dir,
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome");
            assert_eq!(query_sequence, "ACGTACGT");
            assert_eq!(target_seq_id, "query_seq");
            assert_eq!(max_hits, 8);
            assert!(max_hits_explicit);
            assert_eq!(task.as_deref(), Some("blastn"));
            assert!(request_options_json.is_none());
            assert_eq!(track_name.as_deref(), Some("Hits"));
            assert!(clear_existing);
            assert_eq!(catalog_path.as_deref(), Some("c.json"));
            assert_eq!(cache_dir.as_deref(), Some("cache"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_helpers_blast_track_defaults() {
    let cmd =
        parse_shell_line("helpers blast-track pUC19 ACGTAG target_seq").expect("parse command");
    match cmd {
        ShellCommand::ReferenceBlastTrack {
            helper_mode,
            genome_id,
            query_sequence,
            target_seq_id,
            max_hits,
            max_hits_explicit,
            task,
            request_options_json,
            track_name,
            clear_existing,
            ..
        } => {
            assert!(helper_mode);
            assert_eq!(genome_id, "pUC19");
            assert_eq!(query_sequence, "ACGTAG");
            assert_eq!(target_seq_id, "target_seq");
            assert_eq!(max_hits, 25);
            assert!(!max_hits_explicit);
            assert!(task.is_none());
            assert!(request_options_json.is_none());
            assert!(track_name.is_none());
            assert!(!clear_existing);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_primers_design_with_backend_overrides() {
    let cmd = parse_shell_line(
        "primers design @request.json --backend primer3 --primer3-exec /opt/primer3/primer3_core",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::PrimersDesign {
            request_json,
            backend,
            primer3_executable,
        } => {
            assert_eq!(request_json, "@request.json");
            assert_eq!(backend, Some(PrimerDesignBackend::Primer3));
            assert_eq!(
                primer3_executable.as_deref(),
                Some("/opt/primer3/primer3_core")
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_primers_design_qpcr_with_backend_overrides() {
    let cmd = parse_shell_line(
            "primers design-qpcr @request.json --backend primer3 --primer3-exec /opt/primer3/primer3_core",
        )
        .expect("parse command");
    match cmd {
        ShellCommand::PrimersDesignQpcr {
            request_json,
            backend,
            primer3_executable,
        } => {
            assert_eq!(request_json, "@request.json");
            assert_eq!(backend, Some(PrimerDesignBackend::Primer3));
            assert_eq!(
                primer3_executable.as_deref(),
                Some("/opt/primer3/primer3_core")
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_primers_prepare_restriction_cloning_request() {
    let cmd = parse_shell_line("primers prepare-restriction-cloning @handoff_request.json")
        .expect("parse restriction-cloning handoff");
    match cmd {
        ShellCommand::PrimersPrepareRestrictionCloning { request_json } => {
            assert_eq!(request_json, "@handoff_request.json");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_primers_restriction_cloning_vector_suggestions_request() {
    let cmd = parse_shell_line("primers restriction-cloning-vector-suggestions vec")
        .expect("parse restriction-cloning vector suggestions");
    match cmd {
        ShellCommand::PrimersRestrictionCloningVectorSuggestions { seq_id } => {
            assert_eq!(seq_id, "vec");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_primers_restriction_cloning_handoff_report_commands() {
    let seed = parse_shell_line(
        "primers seed-restriction-cloning-handoff primer_report_1 vec_a --pair-rank 3 --mode directed_pair --forward-enzyme EcoRI --reverse-enzyme HindIII --forward-leader GC --reverse-leader AT",
    )
    .expect("parse restriction-cloning handoff seed");
    assert!(matches!(
        seed,
        ShellCommand::PrimersSeedRestrictionCloningHandoff {
            primer_report_id,
            destination_vector_seq_id,
            pair_rank_1based,
            mode,
            forward_enzyme,
            reverse_enzyme,
            forward_leader_5prime,
            reverse_leader_5prime,
        }
            if primer_report_id == "primer_report_1"
                && destination_vector_seq_id == "vec_a"
                && pair_rank_1based == 3
                && mode == RestrictionCloningPcrHandoffMode::DirectedPair
                && forward_enzyme.as_deref() == Some("EcoRI")
                && reverse_enzyme.as_deref() == Some("HindIII")
                && forward_leader_5prime.as_deref() == Some("GC")
                && reverse_leader_5prime.as_deref() == Some("AT")
    ));

    let list = parse_shell_line("primers list-restriction-cloning-handoffs")
        .expect("parse restriction-cloning handoff list");
    assert!(matches!(
        list,
        ShellCommand::PrimersListRestrictionCloningHandoffs
    ));

    let show = parse_shell_line("primers show-restriction-cloning-handoff handoff_1")
        .expect("parse restriction-cloning handoff show");
    assert!(matches!(
        show,
        ShellCommand::PrimersShowRestrictionCloningHandoff { report_id } if report_id == "handoff_1"
    ));

    let export =
        parse_shell_line("primers export-restriction-cloning-handoff handoff_1 handoff_1.json")
            .expect("parse restriction-cloning handoff export");
    assert!(matches!(
        export,
        ShellCommand::PrimersExportRestrictionCloningHandoff { report_id, path }
            if report_id == "handoff_1" && path == "handoff_1.json"
    ));
}

#[test]
fn parse_primers_preflight_with_backend_overrides() {
    let cmd = parse_shell_line(
        "primers preflight --backend primer3 --primer3-exec /opt/primer3/primer3_core",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::PrimersPreflight {
            backend,
            primer3_executable,
        } => {
            assert_eq!(backend, Some(PrimerDesignBackend::Primer3));
            assert_eq!(
                primer3_executable.as_deref(),
                Some("/opt/primer3/primer3_core")
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_primers_seed_from_feature_and_splicing() {
    let feature =
        parse_shell_line("primers seed-from-feature seq_a 7").expect("parse seed-from-feature");
    assert!(matches!(
        feature,
        ShellCommand::PrimersSeedFromFeature {
            seq_id,
            feature_id
        } if seq_id == "seq_a" && feature_id == 7
    ));

    let splicing =
        parse_shell_line("primers seed-from-splicing seq_a 11").expect("parse seed-from-splicing");
    assert!(matches!(
        splicing,
        ShellCommand::PrimersSeedFromSplicing {
            seq_id,
            feature_id
        } if seq_id == "seq_a" && feature_id == 11
    ));

    let qpcr_feature = parse_shell_line("primers seed-qpcr-from-feature seq_a 13")
        .expect("parse seed-qpcr-from-feature");
    assert!(matches!(
        qpcr_feature,
        ShellCommand::PrimersSeedQpcrFromFeature {
            seq_id,
            feature_id
        } if seq_id == "seq_a" && feature_id == 13
    ));

    let qpcr_splicing = parse_shell_line("primers seed-qpcr-from-splicing seq_a 17")
        .expect("parse seed-qpcr-from-splicing");
    assert!(matches!(
        qpcr_splicing,
        ShellCommand::PrimersSeedQpcrFromSplicing {
            seq_id,
            feature_id
        } if seq_id == "seq_a" && feature_id == 17
    ));
}

#[test]
fn parse_features_query_with_filters() {
    let cmd = parse_shell_line(
        "features query seq_a --kind CDS --kind-not source --range 10..200 --within --strand reverse --label TP73 --label-regex ^TP73 --qual gene --qual-contains note=promoter --qual-regex product=TP73.* --min-len 30 --max-len 500 --limit 50 --offset 5 --sort start --desc --include-source --include-qualifiers",
    )
    .expect("parse features query");
    match cmd {
        ShellCommand::FeaturesQuery { query } => {
            assert_eq!(query.seq_id, "seq_a");
            assert_eq!(query.kind_in, vec!["CDS".to_string()]);
            assert_eq!(query.kind_not_in, vec!["source".to_string()]);
            assert_eq!(query.start_0based, Some(10));
            assert_eq!(query.end_0based_exclusive, Some(200));
            assert_eq!(query.range_relation, SequenceFeatureRangeRelation::Within);
            assert_eq!(query.strand, SequenceFeatureStrandFilter::Reverse);
            assert_eq!(query.label_contains.as_deref(), Some("TP73"));
            assert_eq!(query.label_regex.as_deref(), Some("^TP73"));
            assert_eq!(query.qualifier_filters.len(), 3);
            assert_eq!(query.min_len_bp, Some(30));
            assert_eq!(query.max_len_bp, Some(500));
            assert_eq!(query.limit, Some(50));
            assert_eq!(query.offset, 5);
            assert_eq!(query.sort_by, SequenceFeatureSortBy::Start);
            assert!(query.descending);
            assert!(query.include_source);
            assert!(query.include_qualifiers);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_features_export_bed_with_restriction_options() {
    let cmd = parse_shell_line(
        "features export-bed seq_a /tmp/seq_a.features.bed --coordinate-mode genomic --include-restriction-sites --restriction-enzyme EcoRI --kind TFBS --label SP1 --sort start",
    )
    .expect("parse features export-bed");
    match cmd {
        ShellCommand::FeaturesExportBed {
            query,
            output,
            coordinate_mode,
            include_restriction_sites,
            restriction_enzymes,
        } => {
            assert_eq!(query.seq_id, "seq_a");
            assert_eq!(output, "/tmp/seq_a.features.bed");
            assert_eq!(coordinate_mode, FeatureBedCoordinateMode::Genomic);
            assert!(include_restriction_sites);
            assert_eq!(restriction_enzymes, vec!["EcoRI".to_string()]);
            assert_eq!(query.kind_in, vec!["TFBS".to_string()]);
            assert_eq!(query.label_contains.as_deref(), Some("SP1"));
            assert_eq!(query.sort_by, SequenceFeatureSortBy::Start);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_features_tfbs_summary_with_context_filters() {
    let cmd = parse_shell_line(
        "features tfbs-summary seq_a --focus 2900..3100 --context 0..6001 --min-focus-count 2 --min-context-count 3 --limit 25",
    )
    .expect("parse features tfbs-summary");
    match cmd {
        ShellCommand::FeaturesTfbsSummary { request } => {
            assert_eq!(request.seq_id, "seq_a");
            assert_eq!(request.focus_start_0based, 2900);
            assert_eq!(request.focus_end_0based_exclusive, 3100);
            assert_eq!(request.context_start_0based, Some(0));
            assert_eq!(request.context_end_0based_exclusive, Some(6001));
            assert_eq!(request.min_focus_occurrences, 2);
            assert_eq!(request.min_context_occurrences, 3);
            assert_eq!(request.limit, Some(25));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_features_tfbs_score_tracks_svg_with_range_and_negative_scores() {
    let cmd = parse_shell_line(
        "features tfbs-score-tracks-svg seq_a /tmp/seq_a.tfbs.svg --motif SP1 --motif TP73 --range 2900..3100 --score-kind true_log_odds_quantile --allow-negative",
    )
    .expect("parse features tfbs-score-tracks-svg");
    match cmd {
        ShellCommand::FeaturesTfbsScoreTracksSvg {
            target,
            motifs,
            score_kind,
            clip_negative,
            output,
        } => {
            assert_eq!(
                target,
                SequenceScanTarget::SeqId {
                    seq_id: "seq_a".to_string(),
                    span_start_0based: Some(2900),
                    span_end_0based_exclusive: Some(3100),
                }
            );
            assert_eq!(motifs, vec!["SP1".to_string(), "TP73".to_string()]);
            assert_eq!(score_kind, TfbsScoreTrackValueKind::TrueLogOddsQuantile);
            assert!(!clip_negative);
            assert_eq!(output, "/tmp/seq_a.tfbs.svg");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_features_tfbs_score_tracks_svg_with_background_tail_score_kind() {
    let cmd = parse_shell_line(
        "features tfbs-score-tracks-svg seq_a /tmp/seq_a.tfbs.svg --motif SP1 --range 2900..3100 --score-kind llr_background_tail_log10",
    )
    .expect("parse features tfbs-score-tracks-svg background-tail");
    match cmd {
        ShellCommand::FeaturesTfbsScoreTracksSvg {
            score_kind,
            clip_negative,
            ..
        } => {
            assert_eq!(score_kind, TfbsScoreTrackValueKind::LlrBackgroundTailLog10);
            assert!(clip_negative);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_features_tfbs_score_tracks_svg_allows_full_span_without_explicit_range() {
    let cmd = parse_shell_line(
        "features tfbs-score-tracks-svg seq_a /tmp/seq_a.tfbs.svg --motif \"Yamanaka factors\" --motif SP1",
    )
    .expect("parse full-span tfbs-score-tracks-svg");
    match cmd {
        ShellCommand::FeaturesTfbsScoreTracksSvg {
            target,
            motifs,
            score_kind,
            clip_negative,
            output,
        } => {
            assert_eq!(
                target,
                SequenceScanTarget::SeqId {
                    seq_id: "seq_a".to_string(),
                    span_start_0based: None,
                    span_end_0based_exclusive: None,
                }
            );
            assert_eq!(
                motifs,
                vec!["Yamanaka factors".to_string(), "SP1".to_string()]
            );
            assert_eq!(score_kind, TfbsScoreTrackValueKind::LlrBits);
            assert!(clip_negative);
            assert_eq!(output, "/tmp/seq_a.tfbs.svg");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_features_tfbs_score_tracks_svg_for_inline_target() {
    let cmd = parse_shell_line(
        "features tfbs-score-tracks-svg --sequence-text ACGTACGTACGT --topology circular --id-hint inline_tp73 --output /tmp/inline.tfbs.svg --motif SP1 --range 2..10",
    )
    .expect("parse inline features tfbs-score-tracks-svg");
    match cmd {
        ShellCommand::FeaturesTfbsScoreTracksSvg {
            target,
            motifs,
            score_kind,
            clip_negative,
            output,
        } => {
            assert_eq!(
                target,
                SequenceScanTarget::InlineSequence {
                    sequence_text: "ACGTACGTACGT".to_string(),
                    topology: InlineSequenceTopology::Circular,
                    id_hint: Some("inline_tp73".to_string()),
                    span_start_0based: Some(2),
                    span_end_0based_exclusive: Some(10),
                }
            );
            assert_eq!(motifs, vec!["SP1".to_string()]);
            assert_eq!(score_kind, TfbsScoreTrackValueKind::LlrBits);
            assert!(clip_negative);
            assert_eq!(output, "/tmp/inline.tfbs.svg");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_features_tfbs_score_track_correlation_svg_with_range_and_negative_scores() {
    let cmd = parse_shell_line(
        "features tfbs-score-track-correlation-svg seq_a /tmp/seq_a.tfbs_corr.svg --motif SP1 --motif MYC --range 2900..3100 --score-kind llr_background_tail_log10 --correlation-metric spearman --signal-source reverse_only --allow-negative",
    )
    .expect("parse features tfbs-score-track-correlation-svg");
    match cmd {
        ShellCommand::FeaturesTfbsScoreTrackCorrelationSvg {
            seq_id,
            motifs,
            start_0based,
            end_0based_exclusive,
            score_kind,
            correlation_metric,
            correlation_signal_source,
            clip_negative,
            output,
        } => {
            assert_eq!(seq_id, "seq_a");
            assert_eq!(motifs, vec!["SP1".to_string(), "MYC".to_string()]);
            assert_eq!(start_0based, 2900);
            assert_eq!(end_0based_exclusive, 3100);
            assert_eq!(score_kind, TfbsScoreTrackValueKind::LlrBackgroundTailLog10);
            assert_eq!(
                correlation_metric,
                TfbsScoreTrackCorrelationMetric::Spearman
            );
            assert_eq!(
                correlation_signal_source,
                TfbsScoreTrackCorrelationSignalSource::ReverseOnly
            );
            assert!(!clip_negative);
            assert_eq!(output, "/tmp/seq_a.tfbs_corr.svg");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_features_tfbs_track_similarity_for_inline_targets_with_species_filter() {
    let cmd = parse_shell_line(
        "features tfbs-track-similarity --sequence-text ACGTACGTACGT --topology circular --id-hint inline_tp73 --anchor-motif SP1 --candidate-motif MYC --candidate-motif PATZ1 --range 2..10 --ranking-metric smoothed_spearman --score-kind llr_background_tail_log10 --species human --include-remote-metadata --limit 5 --path /tmp/inline.tfbs_similarity.json",
    )
    .expect("parse inline features tfbs-track-similarity");
    match cmd {
        ShellCommand::FeaturesTfbsTrackSimilarity {
            target,
            anchor_motif,
            candidate_motifs,
            ranking_metric,
            score_kind,
            clip_negative,
            species_filters,
            include_remote_metadata,
            limit,
            path,
        } => {
            assert_eq!(
                target,
                SequenceScanTarget::InlineSequence {
                    sequence_text: "ACGTACGTACGT".to_string(),
                    topology: InlineSequenceTopology::Circular,
                    id_hint: Some("inline_tp73".to_string()),
                    span_start_0based: Some(2),
                    span_end_0based_exclusive: Some(10),
                }
            );
            assert_eq!(anchor_motif, "SP1");
            assert_eq!(
                candidate_motifs,
                vec!["MYC".to_string(), "PATZ1".to_string()]
            );
            assert_eq!(
                ranking_metric,
                TfbsTrackSimilarityRankingMetric::SmoothedSpearman
            );
            assert_eq!(score_kind, TfbsScoreTrackValueKind::LlrBackgroundTailLog10);
            assert!(clip_negative);
            assert_eq!(species_filters, vec!["human".to_string()]);
            assert!(include_remote_metadata);
            assert_eq!(limit, Some(5));
            assert_eq!(path.as_deref(), Some("/tmp/inline.tfbs_similarity.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_features_restriction_scan_for_stored_and_inline_targets() {
    let cmd = parse_shell_line(
        "features restriction-scan seq_a --range 10..120 --enzyme EcoRI --enzyme SmaI --max-sites-per-enzyme 3 --no-cut-geometry --path /tmp/seq_a.restriction_scan.json",
    )
    .expect("parse features restriction-scan");
    match cmd {
        ShellCommand::FeaturesRestrictionScan {
            target,
            enzymes,
            max_sites_per_enzyme,
            include_cut_geometry,
            path,
        } => {
            match target {
                SequenceScanTarget::SeqId {
                    seq_id,
                    span_start_0based,
                    span_end_0based_exclusive,
                } => {
                    assert_eq!(seq_id, "seq_a");
                    assert_eq!(span_start_0based, Some(10));
                    assert_eq!(span_end_0based_exclusive, Some(120));
                }
                other => panic!("unexpected target: {other:?}"),
            }
            assert_eq!(enzymes, vec!["EcoRI".to_string(), "SmaI".to_string()]);
            assert_eq!(max_sites_per_enzyme, Some(3));
            assert!(!include_cut_geometry);
            assert_eq!(path.as_deref(), Some("/tmp/seq_a.restriction_scan.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let cmd = parse_shell_line(
        "features restriction-scan --sequence-text AAGAATTCCCGGGTTA --topology circular --id-hint toy_seq --range 1..12 --enzyme EcoRI",
    )
    .expect("parse inline features restriction-scan");
    match cmd {
        ShellCommand::FeaturesRestrictionScan { target, .. } => match target {
            SequenceScanTarget::InlineSequence {
                sequence_text,
                topology,
                id_hint,
                span_start_0based,
                span_end_0based_exclusive,
            } => {
                assert_eq!(sequence_text, "AAGAATTCCCGGGTTA");
                assert_eq!(topology, InlineSequenceTopology::Circular);
                assert_eq!(id_hint.as_deref(), Some("toy_seq"));
                assert_eq!(span_start_0based, Some(1));
                assert_eq!(span_end_0based_exclusive, Some(12));
            }
            other => panic!("unexpected target: {other:?}"),
        },
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_features_tfbs_scan_for_stored_and_inline_targets() {
    let cmd = parse_shell_line(
        "features tfbs-scan seq_a --range 10..120 --motif SP1 --motif TP73 --min-llr-bits -1.5 --min-llr-quantile 0.9 --per-tf-min-llr-bits SP1=0.5 --per-tf-min-llr-quantile TP73=0.95 --max-hits 25 --path /tmp/seq_a.tfbs_scan.json",
    )
    .expect("parse features tfbs-scan");
    match cmd {
        ShellCommand::FeaturesTfbsScan {
            target,
            motifs,
            min_llr_bits,
            min_llr_quantile,
            per_tf_thresholds,
            max_hits,
            path,
        } => {
            match target {
                SequenceScanTarget::SeqId {
                    seq_id,
                    span_start_0based,
                    span_end_0based_exclusive,
                } => {
                    assert_eq!(seq_id, "seq_a");
                    assert_eq!(span_start_0based, Some(10));
                    assert_eq!(span_end_0based_exclusive, Some(120));
                }
                other => panic!("unexpected target: {other:?}"),
            }
            assert_eq!(motifs, vec!["SP1".to_string(), "TP73".to_string()]);
            assert_eq!(min_llr_bits, Some(-1.5));
            assert_eq!(min_llr_quantile, Some(0.9));
            assert_eq!(
                per_tf_thresholds,
                vec![
                    TfThresholdOverride {
                        tf: "SP1".to_string(),
                        min_llr_bits: Some(0.5),
                        min_llr_quantile: None,
                    },
                    TfThresholdOverride {
                        tf: "TP73".to_string(),
                        min_llr_bits: None,
                        min_llr_quantile: Some(0.95),
                    },
                ]
            );
            assert_eq!(max_hits, Some(25));
            assert_eq!(path.as_deref(), Some("/tmp/seq_a.tfbs_scan.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let cmd = parse_shell_line(
        "features tfbs-scan --sequence-text TTTTATAAAGGGTATAAATTT --topology circular --id-hint toy_seq --range 1..12 --motif TATAAA --max-hits 5",
    )
    .expect("parse inline features tfbs-scan");
    match cmd {
        ShellCommand::FeaturesTfbsScan {
            target,
            motifs,
            max_hits,
            ..
        } => {
            match target {
                SequenceScanTarget::InlineSequence {
                    sequence_text,
                    topology,
                    id_hint,
                    span_start_0based,
                    span_end_0based_exclusive,
                } => {
                    assert_eq!(sequence_text, "TTTTATAAAGGGTATAAATTT");
                    assert_eq!(topology, InlineSequenceTopology::Circular);
                    assert_eq!(id_hint.as_deref(), Some("toy_seq"));
                    assert_eq!(span_start_0based, Some(1));
                    assert_eq!(span_end_0based_exclusive, Some(12));
                }
                other => panic!("unexpected target: {other:?}"),
            }
            assert_eq!(motifs, vec!["TATAAA".to_string()]);
            assert_eq!(max_hits, Some(5));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_variant_reporter_fragments_with_context_options() {
    let cmd = parse_shell_line(
        "variant reporter-fragments seq_a --variant rs9923231 --gene-label VKORC1 --retain-downstream-from-tss-bp 200 --retain-upstream-beyond-variant-bp 500 --max-candidates 3 --path out.json",
    )
    .expect("parse variant reporter-fragments");
    match cmd {
        ShellCommand::VariantReporterFragments {
            seq_id,
            variant_label_or_id,
            gene_label,
            transcript_id,
            retain_downstream_from_tss_bp,
            retain_upstream_beyond_variant_bp,
            max_candidates,
            path,
        } => {
            assert_eq!(seq_id, "seq_a");
            assert_eq!(variant_label_or_id.as_deref(), Some("rs9923231"));
            assert_eq!(gene_label.as_deref(), Some("VKORC1"));
            assert!(transcript_id.is_none());
            assert_eq!(retain_downstream_from_tss_bp, 200);
            assert_eq!(retain_upstream_beyond_variant_bp, 500);
            assert_eq!(max_candidates, 3);
            assert_eq!(path.as_deref(), Some("out.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_variant_materialize_allele_command() {
    let cmd = parse_shell_line(
        "variant materialize-allele seq_a --allele alternate --variant rs9923231 --output-id seq_alt",
    )
    .expect("parse variant materialize-allele");
    match cmd {
        ShellCommand::VariantMaterializeAllele {
            seq_id,
            variant_label_or_id,
            allele,
            output_id,
        } => {
            assert_eq!(seq_id, "seq_a");
            assert_eq!(variant_label_or_id.as_deref(), Some("rs9923231"));
            assert_eq!(allele, VariantAlleleChoice::Alternate);
            assert_eq!(output_id.as_deref(), Some("seq_alt"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_tracks_import_bed() {
    let cmd = parse_shell_line(
            "tracks import-bed toy_slice test_files/data/peaks.bed.gz --name ChIP --min-score 5 --max-score 50 --clear-existing",
        )
        .expect("parse command");
    match cmd {
        ShellCommand::TracksImportBed {
            seq_id,
            path,
            track_name,
            min_score,
            max_score,
            clear_existing,
        } => {
            assert_eq!(seq_id, "toy_slice".to_string());
            assert_eq!(path, "test_files/data/peaks.bed.gz".to_string());
            assert_eq!(track_name, Some("ChIP".to_string()));
            assert_eq!(min_score, Some(5.0));
            assert_eq!(max_score, Some(50.0));
            assert!(clear_existing);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_tracks_import_bigwig() {
    let cmd = parse_shell_line(
            "tracks import-bigwig toy_slice test_files/data/signal.bw --name RNA --min-score 0.5 --max-score 2.5 --clear-existing",
        )
        .expect("parse command");
    match cmd {
        ShellCommand::TracksImportBigWig {
            seq_id,
            path,
            track_name,
            min_score,
            max_score,
            clear_existing,
        } => {
            assert_eq!(seq_id, "toy_slice".to_string());
            assert_eq!(path, "test_files/data/signal.bw".to_string());
            assert_eq!(track_name, Some("RNA".to_string()));
            assert_eq!(min_score, Some(0.5));
            assert_eq!(max_score, Some(2.5));
            assert!(clear_existing);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_cutrun_list() {
    let cmd = parse_shell_line("cutrun list --catalog assets/cutrun.json --filter ctcf")
        .expect("parse CUT&RUN list");
    match cmd {
        ShellCommand::CutRunList {
            filter,
            catalog_path,
        } => {
            assert_eq!(filter.as_deref(), Some("ctcf"));
            assert_eq!(catalog_path.as_deref(), Some("assets/cutrun.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_cutrun_project() {
    let cmd = parse_shell_line(
        "cutrun project toy_slice toy_ctcf --no-signal --clear-existing --catalog assets/cutrun.json --cache-dir data/cutrun",
    )
    .expect("parse CUT&RUN project");
    match cmd {
        ShellCommand::CutRunProject {
            seq_id,
            dataset_id,
            include_peaks,
            include_signal,
            clear_existing,
            catalog_path,
            cache_dir,
        } => {
            assert_eq!(seq_id, "toy_slice");
            assert_eq!(dataset_id, "toy_ctcf");
            assert!(include_peaks);
            assert!(!include_signal);
            assert!(clear_existing);
            assert_eq!(catalog_path.as_deref(), Some("assets/cutrun.json"));
            assert_eq!(cache_dir.as_deref(), Some("data/cutrun"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_cutrun_interpret() {
    let cmd = parse_shell_line(
        "cutrun interpret toy_slice reads_r1.fastq.gz reads_r2.fastq.gz --format fastq --layout paired_end --flank-bp 0 --report-id toy_cutrun_reads --seed-kmer-len 3 --min-seed-matches 1 --max-mismatches 0 --min-identity 1.0 --max-fragment-span-bp 64 --no-deduplicate-fragments",
    )
    .expect("parse CUT&RUN interpret");
    match cmd {
        ShellCommand::CutRunInterpret {
            seq_id,
            input_r1_path,
            input_r2_path,
            dataset_id,
            catalog_path,
            cache_dir,
            input_format,
            read_layout,
            roi_flank_bp,
            seed_filter,
            align_config,
            deduplicate_fragments,
            report_id,
            checkpoint_path,
            checkpoint_every_reads,
        } => {
            assert_eq!(seq_id, "toy_slice");
            assert_eq!(input_r1_path.as_deref(), Some("reads_r1.fastq.gz"));
            assert_eq!(input_r2_path.as_deref(), Some("reads_r2.fastq.gz"));
            assert!(dataset_id.is_none());
            assert!(catalog_path.is_none());
            assert!(cache_dir.is_none());
            assert_eq!(input_format, CutRunInputFormat::Fastq);
            assert_eq!(read_layout, CutRunReadLayout::PairedEnd);
            assert_eq!(roi_flank_bp, 0);
            assert_eq!(seed_filter.kmer_len, 3);
            assert_eq!(seed_filter.min_seed_matches, 1);
            assert_eq!(align_config.max_mismatches, 0);
            assert!((align_config.min_identity_fraction - 1.0).abs() < f64::EPSILON);
            assert_eq!(align_config.max_fragment_span_bp, 64);
            assert!(!deduplicate_fragments);
            assert_eq!(report_id.as_deref(), Some("toy_cutrun_reads"));
            assert!(checkpoint_path.is_none());
            assert_eq!(checkpoint_every_reads, 500);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_cutrun_interpret_dataset_backed() {
    let cmd = parse_shell_line(
        "cutrun interpret toy_slice --dataset toy_ctcf_reads --catalog assets/cutrun.json --cache-dir data/cutrun --report-id toy_cutrun_reads",
    )
    .expect("parse dataset-backed CUT&RUN interpret");
    match cmd {
        ShellCommand::CutRunInterpret {
            seq_id,
            input_r1_path,
            input_r2_path,
            dataset_id,
            catalog_path,
            cache_dir,
            report_id,
            ..
        } => {
            assert_eq!(seq_id, "toy_slice");
            assert!(input_r1_path.is_none());
            assert!(input_r2_path.is_none());
            assert_eq!(dataset_id.as_deref(), Some("toy_ctcf_reads"));
            assert_eq!(catalog_path.as_deref(), Some("assets/cutrun.json"));
            assert_eq!(cache_dir.as_deref(), Some("data/cutrun"));
            assert_eq!(report_id.as_deref(), Some("toy_cutrun_reads"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_cutrun_list_show_and_export_read_reports() {
    let list_cmd = parse_shell_line("cutrun list-read-reports toy_slice")
        .expect("parse CUT&RUN list-read-reports");
    assert!(matches!(
        list_cmd,
        ShellCommand::CutRunListReadReports { seq_id } if seq_id.as_deref() == Some("toy_slice")
    ));

    let show_cmd = parse_shell_line("cutrun show-read-report toy_cutrun_reads")
        .expect("parse CUT&RUN show-read-report");
    assert!(matches!(
        show_cmd,
        ShellCommand::CutRunShowReadReport { report_id } if report_id == "toy_cutrun_reads"
    ));

    let export_cmd =
        parse_shell_line("cutrun export-coverage toy_cutrun_reads coverage.tsv --kind fragments")
            .expect("parse CUT&RUN export-coverage");
    assert!(matches!(
        export_cmd,
        ShellCommand::CutRunExportCoverage {
            report_id,
            path,
            kind
        } if report_id == "toy_cutrun_reads"
            && path == "coverage.tsv"
            && kind == CutRunCoverageKind::Fragments
    ));
}

#[test]
fn parse_cutrun_inspect_regulatory_support() {
    let cmd = parse_shell_line(
        "cutrun inspect-regulatory-support toy_cutrun_roi --dataset toy_ctcf --read-report toy_cutrun_reads --promoter-search-start 10 --promoter-search-end 180 --neighbor-window-bp 42 --species-filter human --species-filter mouse --path support.json",
    )
    .expect("parse CUT&RUN inspect-regulatory-support");
    assert!(matches!(
        cmd,
        ShellCommand::CutRunInspectRegulatorySupport {
            seq_id,
            dataset_ids,
            read_report_ids,
            promoter_search_start_0based,
            promoter_search_end_0based_exclusive,
            neighbor_window_bp,
            species_filters,
            output_path,
        } if seq_id == "toy_cutrun_roi"
            && dataset_ids == vec!["toy_ctcf".to_string()]
            && read_report_ids == vec!["toy_cutrun_reads".to_string()]
            && promoter_search_start_0based == Some(10)
            && promoter_search_end_0based_exclusive == Some(180)
            && neighbor_window_bp == 42
            && species_filters == vec!["human".to_string(), "mouse".to_string()]
            && output_path.as_deref() == Some("support.json")
    ));
}

#[test]
fn parse_tracks_import_vcf() {
    let cmd = parse_shell_line(
            "tracks import-vcf toy_slice test_files/data/variants.vcf.gz --name SNPs --min-score 10 --max-score 60 --clear-existing",
        )
        .expect("parse command");
    match cmd {
        ShellCommand::TracksImportVcf {
            seq_id,
            path,
            track_name,
            min_score,
            max_score,
            clear_existing,
        } => {
            assert_eq!(seq_id, "toy_slice".to_string());
            assert_eq!(path, "test_files/data/variants.vcf.gz".to_string());
            assert_eq!(track_name, Some("SNPs".to_string()));
            assert_eq!(min_score, Some(10.0));
            assert_eq!(max_score, Some(60.0));
            assert!(clear_existing);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_tracks_tracked_add_with_options() {
    let cmd = parse_shell_line(
            "tracks tracked add test_files/data/signal.bw --source bed --name ChIP --min-score 0.5 --max-score 2.5 --clear-existing",
        )
        .expect("parse command");
    match cmd {
        ShellCommand::TracksTrackedAdd { subscription } => {
            assert_eq!(subscription.source, GenomeTrackSource::Bed);
            assert_eq!(subscription.path, "test_files/data/signal.bw");
            assert_eq!(subscription.track_name.as_deref(), Some("ChIP"));
            assert_eq!(subscription.min_score, Some(0.5));
            assert_eq!(subscription.max_score, Some(2.5));
            assert!(subscription.clear_existing);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_tracks_tracked_add_with_vcf_source() {
    let cmd = parse_shell_line(
        "tracks tracked add test_files/data/variants.vcf.gz --source vcf --name Variants",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::TracksTrackedAdd { subscription } => {
            assert_eq!(subscription.source, GenomeTrackSource::Vcf);
            assert_eq!(subscription.path, "test_files/data/variants.vcf.gz");
            assert_eq!(subscription.track_name.as_deref(), Some("Variants"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_tracks_tracked_apply_with_index() {
    let cmd = parse_shell_line("tracks tracked apply --index 3 --only-new-anchors")
        .expect("parse command");
    match cmd {
        ShellCommand::TracksTrackedApply {
            index,
            only_new_anchors,
        } => {
            assert_eq!(index, Some(3));
            assert!(only_new_anchors);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_tracks_tracked_sync_alias() {
    let cmd = parse_shell_line("tracks tracked sync").expect("parse command");
    match cmd {
        ShellCommand::TracksTrackedApply {
            index,
            only_new_anchors,
        } => {
            assert_eq!(index, None);
            assert!(only_new_anchors);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_routines_list_with_filters() {
    let cmd = parse_shell_line(
            "routines list --catalog assets/cloning_routines.json --family crispr --status implemented --tag guide --query scan",
        )
        .expect("parse routines list");
    match cmd {
        ShellCommand::RoutinesList {
            catalog_path,
            family,
            status,
            tag,
            query,
            seq_id,
        } => {
            assert_eq!(
                catalog_path.as_deref(),
                Some("assets/cloning_routines.json")
            );
            assert_eq!(family.as_deref(), Some("crispr"));
            assert_eq!(status.as_deref(), Some("implemented"));
            assert_eq!(tag.as_deref(), Some("guide"));
            assert_eq!(query.as_deref(), Some("scan"));
            assert_eq!(seq_id, None);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_routines_explain_with_catalog() {
    let cmd = parse_shell_line(
        "routines explain golden_gate.type_iis_single_insert --catalog catalog.json",
    )
    .expect("parse routines explain");
    match cmd {
        ShellCommand::RoutinesExplain {
            catalog_path,
            routine_id,
            seq_id,
        } => {
            assert_eq!(catalog_path.as_deref(), Some("catalog.json"));
            assert_eq!(routine_id, "golden_gate.type_iis_single_insert".to_string());
            assert_eq!(seq_id, None);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_routines_explain_with_sequence_context() {
    let cmd = parse_shell_line(
        "routines explain golden_gate.type_iis_single_insert --seq-id seq_variant",
    )
    .expect("parse routines explain");
    match cmd {
        ShellCommand::RoutinesExplain { seq_id, .. } => {
            assert_eq!(seq_id.as_deref(), Some("seq_variant"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_routines_list_with_sequence_context() {
    let cmd = parse_shell_line("routines list --query reporter --seq-id seq_variant")
        .expect("parse routines list");
    match cmd {
        ShellCommand::RoutinesList { query, seq_id, .. } => {
            assert_eq!(query.as_deref(), Some("reporter"));
            assert_eq!(seq_id.as_deref(), Some("seq_variant"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_routines_compare_with_catalog() {
    let cmd = parse_shell_line(
            "routines compare golden_gate.type_iis_single_insert gibson.two_fragment_overlap_preview --catalog catalog.json",
        )
        .expect("parse routines compare");
    match cmd {
        ShellCommand::RoutinesCompare {
            catalog_path,
            left_routine_id,
            right_routine_id,
            seq_id,
        } => {
            assert_eq!(catalog_path.as_deref(), Some("catalog.json"));
            assert_eq!(
                left_routine_id,
                "golden_gate.type_iis_single_insert".to_string()
            );
            assert_eq!(
                right_routine_id,
                "gibson.two_fragment_overlap_preview".to_string()
            );
            assert_eq!(seq_id, None);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_routines_compare_with_sequence_context() {
    let cmd = parse_shell_line(
        "routines compare golden_gate.type_iis_single_insert gibson.two_fragment_overlap_preview --seq-id seq_variant",
    )
    .expect("parse routines compare");
    match cmd {
        ShellCommand::RoutinesCompare { seq_id, .. } => {
            assert_eq!(seq_id.as_deref(), Some("seq_variant"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_planning_commands() {
    let profile = parse_shell_line("planning profile set @profile.json --scope global")
        .expect("parse planning profile set");
    assert!(matches!(
        profile,
        ShellCommand::PlanningProfileSet { scope, payload_json }
            if scope == PlanningProfileScope::Global && payload_json == "@profile.json"
    ));

    let objective =
        parse_shell_line("planning objective show").expect("parse planning objective show");
    assert!(matches!(objective, ShellCommand::PlanningObjectiveShow));

    let suggestions = parse_shell_line("planning suggestions list --status pending")
        .expect("parse planning suggestions list");
    assert!(matches!(
        suggestions,
        ShellCommand::PlanningSuggestionsList { status }
            if status == Some(PlanningSuggestionStatus::Pending)
    ));

    let sync = parse_shell_line(
        "planning sync pull @sync.json --source lab_manager --confidence 0.75 --snapshot-id snap_1",
    )
    .expect("parse planning sync pull");
    assert!(matches!(
        sync,
        ShellCommand::PlanningSyncPull {
            payload_json,
            source,
            confidence,
            snapshot_id
        } if payload_json == "@sync.json"
            && source.as_deref() == Some("lab_manager")
            && confidence == Some(0.75)
            && snapshot_id.as_deref() == Some("snap_1")
    ));
}

#[test]
fn parse_candidates_generate_with_feature_filters() {
    let cmd = parse_shell_line(
            "candidates generate set1 seqA --length 20 --step 5 --feature-kind gene --feature-kind CDS --feature-label-regex '^TP53$' --max-distance 100 --limit 50",
        )
        .expect("parse candidates generate");
    match cmd {
        ShellCommand::CandidatesGenerate {
            set_name,
            seq_id,
            length_bp,
            step_bp,
            feature_kinds,
            feature_label_regex,
            max_distance_bp,
            feature_geometry_mode,
            feature_boundary_mode,
            feature_strand_relation,
            limit,
        } => {
            assert_eq!(set_name, "set1");
            assert_eq!(seq_id, "seqA");
            assert_eq!(length_bp, 20);
            assert_eq!(step_bp, 5);
            assert_eq!(feature_kinds, vec!["gene".to_string(), "CDS".to_string()]);
            assert_eq!(feature_label_regex.as_deref(), Some("^TP53$"));
            assert_eq!(max_distance_bp, Some(100));
            assert_eq!(feature_geometry_mode, None);
            assert_eq!(feature_boundary_mode, None);
            assert_eq!(feature_strand_relation, None);
            assert_eq!(limit, 50);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_generate_with_strand_relation() {
    let cmd = parse_shell_line("candidates generate set1 seqA --length 20 --strand-relation same")
        .expect("parse candidates generate with strand relation");
    match cmd {
        ShellCommand::CandidatesGenerate {
            set_name,
            seq_id,
            length_bp,
            step_bp,
            feature_kinds,
            feature_label_regex,
            max_distance_bp,
            feature_geometry_mode,
            feature_boundary_mode,
            feature_strand_relation,
            limit,
        } => {
            assert_eq!(set_name, "set1");
            assert_eq!(seq_id, "seqA");
            assert_eq!(length_bp, 20);
            assert_eq!(step_bp, 1);
            assert!(feature_kinds.is_empty());
            assert_eq!(feature_label_regex, None);
            assert_eq!(max_distance_bp, None);
            assert_eq!(feature_geometry_mode, None);
            assert_eq!(feature_boundary_mode, None);
            assert_eq!(
                feature_strand_relation,
                Some(CandidateFeatureStrandRelation::Same)
            );
            assert_eq!(limit, DEFAULT_CANDIDATE_SET_LIMIT);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_generate_between_anchors_with_positions() {
    let cmd = parse_shell_line(
            "candidates generate-between-anchors set1 seqA --length 20 --anchor-a-pos 5 --anchor-b-pos 105 --step 5 --limit 50",
        )
        .expect("parse candidates generate-between-anchors");
    match cmd {
        ShellCommand::CandidatesGenerateBetweenAnchors {
            set_name,
            seq_id,
            anchor_a,
            anchor_b,
            length_bp,
            step_bp,
            limit,
        } => {
            assert_eq!(set_name, "set1");
            assert_eq!(seq_id, "seqA");
            assert_eq!(length_bp, 20);
            assert_eq!(step_bp, 5);
            assert_eq!(limit, 50);
            assert_eq!(anchor_a, SequenceAnchor::Position { zero_based: 5 });
            assert_eq!(anchor_b, SequenceAnchor::Position { zero_based: 105 });
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_score_distance_with_geometry_and_boundary() {
    let cmd = parse_shell_line(
            "candidates score-distance set1 dist --feature-kind exon --feature-geometry feature_boundaries --feature-boundary five_prime",
        )
        .expect("parse candidates score-distance");
    match cmd {
        ShellCommand::CandidatesScoreDistance {
            set_name,
            metric,
            feature_kinds,
            feature_label_regex,
            feature_geometry_mode,
            feature_boundary_mode,
            feature_strand_relation,
        } => {
            assert_eq!(set_name, "set1");
            assert_eq!(metric, "dist");
            assert_eq!(feature_kinds, vec!["exon".to_string()]);
            assert_eq!(feature_label_regex, None);
            assert_eq!(
                feature_geometry_mode,
                Some(CandidateFeatureGeometryMode::FeatureBoundaries)
            );
            assert_eq!(
                feature_boundary_mode,
                Some(CandidateFeatureBoundaryMode::FivePrime)
            );
            assert_eq!(feature_strand_relation, None);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_score_distance_with_strand_relation() {
    let cmd = parse_shell_line(
        "candidates score-distance set1 dist --feature-kind gene --strand-relation opposite",
    )
    .expect("parse candidates score-distance with strand relation");
    match cmd {
        ShellCommand::CandidatesScoreDistance {
            set_name,
            metric,
            feature_kinds,
            feature_label_regex,
            feature_geometry_mode,
            feature_boundary_mode,
            feature_strand_relation,
        } => {
            assert_eq!(set_name, "set1");
            assert_eq!(metric, "dist");
            assert_eq!(feature_kinds, vec!["gene".to_string()]);
            assert_eq!(feature_label_regex, None);
            assert_eq!(feature_geometry_mode, None);
            assert_eq!(feature_boundary_mode, None);
            assert_eq!(
                feature_strand_relation,
                Some(CandidateFeatureStrandRelation::Opposite)
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_score_expression_preserves_formula() {
    let cmd =
        parse_shell_line("candidates score set1 custom_score '(gc_fraction*100)+(length_bp/2)'")
            .expect("parse candidates score");
    match cmd {
        ShellCommand::CandidatesScoreExpression {
            set_name,
            metric,
            expression,
        } => {
            assert_eq!(set_name, "set1");
            assert_eq!(metric, "custom_score");
            assert_eq!(expression, "(gc_fraction*100)+(length_bp/2)");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_score_weighted() {
    let cmd = parse_shell_line(
            "candidates score-weighted set1 objective --term gc_fraction:0.7:max --term distance_to_seq_start_bp:0.3:min --normalize",
        )
        .expect("parse candidates score-weighted");
    match cmd {
        ShellCommand::CandidatesScoreWeightedObjective {
            set_name,
            metric,
            objectives,
            normalize_metrics,
        } => {
            assert_eq!(set_name, "set1");
            assert_eq!(metric, "objective");
            assert_eq!(objectives.len(), 2);
            assert!(normalize_metrics);
            assert_eq!(objectives[0].metric, "gc_fraction");
            assert_eq!(
                objectives[1].direction,
                CandidateObjectiveDirection::Minimize
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_top_k() {
    let cmd = parse_shell_line(
            "candidates top-k in_set out_set --metric objective --k 5 --direction max --tie-break length_descending",
        )
        .expect("parse candidates top-k");
    match cmd {
        ShellCommand::CandidatesTopK {
            input_set,
            output_set,
            metric,
            k,
            direction,
            tie_break,
        } => {
            assert_eq!(input_set, "in_set");
            assert_eq!(output_set, "out_set");
            assert_eq!(metric, "objective");
            assert_eq!(k, 5);
            assert_eq!(direction, CandidateObjectiveDirection::Maximize);
            assert_eq!(tie_break, CandidateTieBreakPolicy::LengthDescending);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_pareto() {
    let cmd = parse_shell_line(
            "candidates pareto in_set out_set --objective gc_fraction:max --objective distance_to_seq_start_bp:min --max-candidates 10 --tie-break seq_start_end",
        )
        .expect("parse candidates pareto");
    match cmd {
        ShellCommand::CandidatesParetoFrontier {
            input_set,
            output_set,
            objectives,
            max_candidates,
            tie_break,
        } => {
            assert_eq!(input_set, "in_set");
            assert_eq!(output_set, "out_set");
            assert_eq!(objectives.len(), 2);
            assert_eq!(max_candidates, Some(10));
            assert_eq!(tie_break, CandidateTieBreakPolicy::SeqStartEnd);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_template_put_and_run() {
    let put = parse_shell_line(
            "candidates template-put scan --script 'generate ${set_name} ${seq_id} --length ${len}' --details-url https://example.org/candidates/scan --param set_name --param seq_id=seqA --param len=20",
        )
        .expect("parse template-put");
    match put {
        ShellCommand::CandidatesTemplateUpsert {
            name,
            details_url,
            parameters,
            script,
            ..
        } => {
            assert_eq!(name, "scan");
            assert_eq!(
                details_url.as_deref(),
                Some("https://example.org/candidates/scan")
            );
            assert_eq!(parameters.len(), 3);
            assert_eq!(script, "generate ${set_name} ${seq_id} --length ${len}");
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let run = parse_shell_line(
        "candidates template-run scan --bind set_name=hits --bind seq_id=seqB --transactional",
    )
    .expect("parse template-run");
    match run {
        ShellCommand::CandidatesTemplateRun {
            name,
            bindings,
            transactional,
        } => {
            assert_eq!(name, "scan");
            assert_eq!(bindings.get("set_name"), Some(&"hits".to_string()));
            assert_eq!(bindings.get("seq_id"), Some(&"seqB".to_string()));
            assert!(transactional);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_template_put_accepts_url_alias() {
    let put = parse_shell_line(
            "candidates template-put scan --script 'generate ${set_name} ${seq_id} --length 20' --url https://example.org/candidates/scan --param set_name --param seq_id",
        )
        .expect("parse candidates template-put with --url");
    match put {
        ShellCommand::CandidatesTemplateUpsert { details_url, .. } => {
            assert_eq!(
                details_url.as_deref(),
                Some("https://example.org/candidates/scan")
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_guides_filter_and_oligos_export_with_options() {
    let filter = parse_shell_line(
        "guides filter tp73 --config '{\"gc_min\":0.35,\"gc_max\":0.7}' --output-set tp73_pass",
    )
    .expect("parse guides filter");
    match filter {
        ShellCommand::GuidesFilter {
            guide_set_id,
            config_json,
            output_guide_set_id,
        } => {
            assert_eq!(guide_set_id, "tp73".to_string());
            assert!(config_json.is_some());
            assert_eq!(output_guide_set_id, Some("tp73_pass".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let export = parse_shell_line(
        "guides oligos-export tp73 out.csv --format plate_csv --plate 96 --oligo-set tp73_set",
    )
    .expect("parse guides oligos-export");
    match export {
        ShellCommand::GuidesOligosExport {
            guide_set_id,
            oligo_set_id,
            format,
            path,
            plate_format,
        } => {
            assert_eq!(guide_set_id, "tp73");
            assert_eq!(oligo_set_id.as_deref(), Some("tp73_set"));
            assert_eq!(format, GuideOligoExportFormat::PlateCsv);
            assert_eq!(path, "out.csv");
            assert_eq!(plate_format, Some(GuideOligoPlateFormat::Plate96));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_macros_run_and_template_commands() {
    let run = parse_shell_line("macros run --transactional --file test_files/workflow_plan.gsh")
        .expect("parse macros run");
    match run {
        ShellCommand::MacrosRun {
            script,
            transactional,
        } => {
            assert_eq!(script, "@test_files/workflow_plan.gsh");
            assert!(transactional);
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let put = parse_shell_line(
            r#"macros template-put clone --script 'op {"Reverse":{"input":"${seq_id}","output_id":"${out_id}"}}' --details-url https://example.org/clone --param seq_id --param out_id=seqA_rev"#,
        )
        .expect("parse macros template-put");
    match put {
        ShellCommand::MacrosTemplateUpsert {
            name,
            details_url,
            parameters,
            script,
            ..
        } => {
            assert_eq!(name, "clone");
            assert_eq!(details_url.as_deref(), Some("https://example.org/clone"));
            assert_eq!(parameters.len(), 2);
            assert_eq!(
                script,
                r#"op {"Reverse":{"input":"${seq_id}","output_id":"${out_id}"}}"#
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let run_template = parse_shell_line(
        "macros template-run clone --bind seq_id=seqB --bind out_id=seqB_rev --transactional",
    )
    .expect("parse macros template-run");
    match run_template {
        ShellCommand::MacrosTemplateRun {
            name,
            bindings,
            transactional,
            validate_only,
        } => {
            assert_eq!(name, "clone");
            assert_eq!(bindings.get("seq_id"), Some(&"seqB".to_string()));
            assert_eq!(bindings.get("out_id"), Some(&"seqB_rev".to_string()));
            assert!(transactional);
            assert!(!validate_only);
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let import = parse_shell_line("macros template-import assets/cloning_patterns.json")
        .expect("parse macros template-import");
    match import {
        ShellCommand::MacrosTemplateImport { path } => {
            assert_eq!(path, "assets/cloning_patterns.json");
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let instance_list =
        parse_shell_line("macros instance-list").expect("parse macros instance-list");
    match instance_list {
        ShellCommand::MacrosInstanceList => {}
        other => panic!("unexpected command: {other:?}"),
    }

    let instance_show =
        parse_shell_line("macros instance-show macro-1").expect("parse macros instance-show");
    match instance_show {
        ShellCommand::MacrosInstanceShow { macro_instance_id } => {
            assert_eq!(macro_instance_id, "macro-1");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_macros_template_put_accepts_url_alias() {
    let put = parse_shell_line(
            r#"macros template-put clone --script 'op {"Reverse":{"input":"${seq_id}"}}' --url https://example.org/clone --param seq_id"#,
        )
        .expect("parse macros template-put with --url");
    match put {
        ShellCommand::MacrosTemplateUpsert { details_url, .. } => {
            assert_eq!(details_url.as_deref(), Some("https://example.org/clone"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_macros_template_put_accepts_port_contracts() {
    let put = parse_shell_line(
            r#"macros template-put clone --script 'op {"Reverse":{"input":"${seq_id}","output_id":"${out_id}"}}' --param seq_id --param out_id=seqA_rev --input-port seq_id:sequence:one:required:Template --output-port out_id:sequence:one:optional:Derived"#,
        )
        .expect("parse macros template-put with ports");
    match put {
        ShellCommand::MacrosTemplateUpsert {
            input_ports,
            output_ports,
            ..
        } => {
            assert_eq!(input_ports.len(), 1);
            assert_eq!(input_ports[0].port_id, "seq_id");
            assert_eq!(input_ports[0].kind, "sequence");
            assert_eq!(input_ports[0].cardinality, "one");
            assert!(input_ports[0].required);
            assert_eq!(output_ports.len(), 1);
            assert_eq!(output_ports[0].port_id, "out_id");
            assert_eq!(output_ports[0].kind, "sequence");
            assert_eq!(output_ports[0].cardinality, "one");
            assert!(!output_ports[0].required);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_macros_template_run_validate_only_flag() {
    let cmd = parse_shell_line("macros template-run clone --validate-only --bind seq_id=seqA")
        .expect("parse macros template-run validate-only");
    match cmd {
        ShellCommand::MacrosTemplateRun {
            name,
            bindings,
            transactional,
            validate_only,
        } => {
            assert_eq!(name, "clone");
            assert_eq!(bindings.get("seq_id").map(String::as_str), Some("seqA"));
            assert!(!transactional);
            assert!(validate_only);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_macro_file_reference() {
    let cmd = parse_shell_line("candidates macro --file test_files/candidates_plan.gsh")
        .expect("parse candidates macro --file");
    match cmd {
        ShellCommand::CandidatesMacro {
            script,
            transactional,
        } => {
            assert_eq!(script, "@test_files/candidates_plan.gsh");
            assert!(!transactional);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_candidates_macro_transactional_flag() {
    let cmd =
        parse_shell_line("candidates macro --transactional --file test_files/candidates_plan.gsh")
            .expect("parse transactional candidates macro");
    match cmd {
        ShellCommand::CandidatesMacro {
            script,
            transactional,
        } => {
            assert_eq!(script, "@test_files/candidates_plan.gsh");
            assert!(transactional);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn split_candidates_macro_statements_handles_quotes_and_comments() {
    let script = r#"
# comment line
generate set1 seqA --length 4 --step 2;
score set1 score "(gc_fraction*100)+1";
filter set1 set2 --metric score --min 10
"#;
    let statements =
        split_candidates_macro_statements(script).expect("split candidates macro statements");
    assert_eq!(statements.len(), 3);
    assert_eq!(statements[0], "generate set1 seqA --length 4 --step 2");
    assert_eq!(statements[1], "score set1 score \"(gc_fraction*100)+1\"");
    assert_eq!(statements[2], "filter set1 set2 --metric score --min 10");
}

#[test]
fn execute_candidates_generate_score_distance_and_filter() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(5, 6),
        qualifiers: vec![("label".into(), Some("target".to_string()))],
    });
    state.sequences.insert("seqA".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let generated = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesGenerate {
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
            limit: 10,
        },
    )
    .expect("generate candidates");
    assert!(generated.state_changed);
    assert_eq!(generated.output["set_name"].as_str(), Some("windows"));
    assert!(
        generated.output["result"]["messages"]
            .as_array()
            .map(|messages| !messages.is_empty())
            .unwrap_or(false)
    );

    let score_distance = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesScoreDistance {
            set_name: "windows".to_string(),
            metric: "dist_gene".to_string(),
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
        },
    )
    .expect("score distance");
    assert!(score_distance.state_changed);

    let filtered = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesFilter {
            input_set: "windows".to_string(),
            output_set: "near_gene".to_string(),
            metric: "dist_gene".to_string(),
            min: None,
            max: Some(0.0),
            min_quantile: None,
            max_quantile: None,
        },
    )
    .expect("filter by distance");
    assert!(filtered.state_changed);
    assert_eq!(filtered.output["output_set"].as_str(), Some("near_gene"));

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesShow {
            set_name: "near_gene".to_string(),
            limit: 10,
            offset: 0,
        },
    )
    .expect("show near_gene");
    assert!(!shown.state_changed);
    assert_eq!(shown.output["candidate_count"].as_u64(), Some(1));
}

#[test]
fn execute_candidates_generate_between_anchors_creates_set() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesGenerateBetweenAnchors {
            set_name: "between".to_string(),
            seq_id: "seqA".to_string(),
            anchor_a: SequenceAnchor::Position { zero_based: 2 },
            anchor_b: SequenceAnchor::Position { zero_based: 10 },
            length_bp: 4,
            step_bp: 2,
            limit: 64,
        },
    )
    .expect("generate candidates between anchors");
    assert!(out.state_changed);
    assert_eq!(out.output["set_name"].as_str(), Some("between"));

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesShow {
            set_name: "between".to_string(),
            limit: 64,
            offset: 0,
        },
    )
    .expect("show between");
    assert_eq!(shown.output["candidate_count"].as_u64(), Some(3));
}

#[test]
fn execute_candidates_score_distance_honors_strand_relation() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(2, 5),
        qualifiers: vec![("label".into(), Some("PLUS_GENE".to_string()))],
    });
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::Complement(Box::new(Location::simple_range(12, 15))),
        qualifiers: vec![("label".into(), Some("MINUS_GENE".to_string()))],
    });
    state.sequences.insert("seqA".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesGenerate {
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
            limit: 128,
        },
    )
    .expect("generate candidate windows");

    execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesScoreDistance {
            set_name: "windows".to_string(),
            metric: "dist_same".to_string(),
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: Some(CandidateFeatureStrandRelation::Same),
        },
    )
    .expect("score distance same");
    execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesScoreDistance {
            set_name: "windows".to_string(),
            metric: "dist_opposite".to_string(),
            feature_kinds: vec!["gene".to_string()],
            feature_label_regex: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: Some(CandidateFeatureStrandRelation::Opposite),
        },
    )
    .expect("score distance opposite");

    let (page, _, _) = engine
        .inspect_candidate_set_page("windows", 1024, 0)
        .expect("inspect scored windows");
    let at_pos = |pos: usize, metric: &str| -> f64 {
        page.candidates
            .iter()
            .find(|candidate| candidate.start_0based == pos)
            .and_then(|candidate| candidate.metrics.get(metric).copied())
            .unwrap_or(f64::NAN)
    };

    let plus_same = at_pos(2, "dist_same");
    let plus_opposite = at_pos(2, "dist_opposite");
    assert_eq!(plus_same, 0.0);
    assert!(plus_opposite > 0.0);

    let minus_same = at_pos(14, "dist_same");
    let minus_opposite = at_pos(14, "dist_opposite");
    assert!(minus_same > 0.0);
    assert_eq!(minus_opposite, 0.0);
}

#[test]
fn execute_candidates_set_operations() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesGenerate {
            set_name: "left".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 4,
            step_bp: 4,
            feature_kinds: vec![],
            feature_label_regex: None,
            max_distance_bp: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: 10,
        },
    )
    .expect("generate left");
    execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesGenerate {
            set_name: "right".to_string(),
            seq_id: "seqA".to_string(),
            length_bp: 4,
            step_bp: 2,
            feature_kinds: vec![],
            feature_label_regex: None,
            max_distance_bp: None,
            feature_geometry_mode: None,
            feature_boundary_mode: None,
            feature_strand_relation: None,
            limit: 10,
        },
    )
    .expect("generate right");

    let intersect = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesSetOp {
            op: CandidateSetOperator::Intersect,
            left_set: "left".to_string(),
            right_set: "right".to_string(),
            output_set: "inter".to_string(),
        },
    )
    .expect("set intersect");
    assert!(intersect.state_changed);
    assert_eq!(intersect.output["operator"].as_str(), Some("intersect"));
    assert_eq!(intersect.output["output_set"].as_str(), Some("inter"));
}

#[test]
fn execute_candidates_optimizer_primitives() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("GCATGAAA").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesGenerate {
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
            limit: 64,
        },
    )
    .expect("generate candidate set");

    let weighted = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesScoreWeightedObjective {
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
            normalize_metrics: true,
        },
    )
    .expect("weighted objective");
    assert!(weighted.state_changed);

    let topk = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesTopK {
            input_set: "cand".to_string(),
            output_set: "top".to_string(),
            metric: "objective".to_string(),
            k: 1,
            direction: CandidateObjectiveDirection::Maximize,
            tie_break: CandidateTieBreakPolicy::SeqStartEnd,
        },
    )
    .expect("top-k");
    assert!(topk.state_changed);

    let pareto = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesParetoFrontier {
            input_set: "cand".to_string(),
            output_set: "front".to_string(),
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
            tie_break: CandidateTieBreakPolicy::SeqStartEnd,
        },
    )
    .expect("pareto");
    assert!(pareto.state_changed);

    let top_set = engine
        .inspect_candidate_set_page("top", 10, 0)
        .expect("inspect top")
        .0;
    assert_eq!(top_set.candidates.len(), 1);
    assert_eq!(top_set.candidates[0].start_0based, 0);
}

#[test]
fn execute_candidates_template_registry_and_run() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesTemplateUpsert {
            name: "scan".to_string(),
            description: Some("scan template".to_string()),
            details_url: Some("https://example.org/candidates/scan".to_string()),
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
            ],
            script: "generate ${set_name} ${seq_id} --length 4 --step 2".to_string(),
        },
    )
    .expect("upsert template");

    let listed = execute_shell_command(&mut engine, &ShellCommand::CandidatesTemplateList)
        .expect("list templates");
    assert_eq!(listed.output["template_count"].as_u64(), Some(1));
    let details_url = listed
        .output
        .get("templates")
        .and_then(|v| v.as_array())
        .and_then(|arr| arr.first())
        .and_then(|template| template.get("details_url"))
        .and_then(|v| v.as_str());
    assert_eq!(details_url, Some("https://example.org/candidates/scan"));

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesTemplateRun {
            name: "scan".to_string(),
            bindings: HashMap::from([("set_name".to_string(), "hits".to_string())]),
            transactional: false,
        },
    )
    .expect("run template");
    assert!(run.state_changed);
    assert_eq!(run.output["template_name"].as_str(), Some("scan"));
    assert!(
        engine
            .list_candidate_sets()
            .iter()
            .any(|set| set.name == "hits")
    );
}

#[test]
fn execute_routines_list_filters_and_searches() {
    let mut engine = GentleEngine::default();
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("cloning_routines.json");
    fs::write(
        &catalog_path,
        r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "restriction.digest_basic",
      "title": "Restriction Digest Basic",
      "family": "restriction",
      "status": "implemented",
      "vocabulary_tags": ["restriction", "digest", "sticky"],
      "template_name": "digest_ligate_extract_sticky",
      "input_ports": [
        { "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    },
    {
      "routine_id": "crispr.guides.scan_basic",
      "title": "Guide Candidate Scan",
      "family": "crispr",
      "status": "planned",
      "vocabulary_tags": ["crispr", "guide", "scan"],
      "template_name": "grna_candidate_priority_scan",
      "input_ports": [
        { "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "set_name", "kind": "candidate_set", "required": true, "cardinality": "one" }
      ]
    }
  ]
}"#,
    )
    .expect("write routines catalog");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::RoutinesList {
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            family: Some("crispr".to_string()),
            status: Some("planned".to_string()),
            tag: Some("guide".to_string()),
            query: Some("scan".to_string()),
            seq_id: None,
        },
    )
    .expect("list routines");
    assert!(!run.state_changed);
    assert_eq!(
        run.output["schema"].as_str(),
        Some(CLONING_ROUTINE_LIST_SCHEMA)
    );
    assert_eq!(run.output["routine_count"].as_u64(), Some(1));
    let rows = run
        .output
        .get("routines")
        .and_then(|value| value.as_array())
        .cloned()
        .unwrap_or_default();
    assert_eq!(rows.len(), 1);
    assert_eq!(
        rows[0].get("routine_id").and_then(|value| value.as_str()),
        Some("crispr.guides.scan_basic")
    );
}

#[test]
fn execute_planning_profile_merge_precedence_global_agent_project() {
    let mut engine = GentleEngine::default();

    execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningProfileSet {
                scope: PlanningProfileScope::Global,
                payload_json: r#"{
                  "schema":"gentle.planning_profile.v1",
                  "capabilities":["gel_imager"],
                  "inventory":{"enzymes":{"available":false,"procurement_business_days":3}},
                  "machine_availability":{"thermocycler":{"available":true,"queue_business_days":0.5}}
                }"#
                .to_string(),
            },
        )
        .expect("set global profile");

    execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningProfileSet {
            scope: PlanningProfileScope::ConfirmedAgentOverlay,
            payload_json: r#"{
                  "schema":"gentle.planning_profile.v1",
                  "capabilities":["thermocycler"],
                  "inventory":{"enzymes":{"available":true}}
                }"#
            .to_string(),
        },
    )
    .expect("set agent overlay");

    execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningProfileSet {
            scope: PlanningProfileScope::ProjectOverride,
            payload_json: r#"{
                  "schema":"gentle.planning_profile.v1",
                  "inventory":{"enzymes":{"available":false,"procurement_business_days":2}}
                }"#
            .to_string(),
        },
    )
    .expect("set project override");

    let effective = engine.planning_effective_profile();
    let enzymes = effective
        .inventory
        .get("enzymes")
        .expect("enzymes inventory");
    assert!(
        !enzymes.available,
        "project override should win over agent/global"
    );
    assert_eq!(enzymes.procurement_business_days, Some(2.0));
    assert!(
        effective.capabilities.iter().any(|cap| cap == "gel_imager"),
        "global capability should remain in merged profile"
    );
    assert!(
        effective
            .capabilities
            .iter()
            .any(|cap| cap == "thermocycler"),
        "agent overlay capability should be merged"
    );
    let machine = effective
        .machine_availability
        .get("thermocycler")
        .expect("thermocycler machine entry");
    assert!(
        (machine.queue_business_days - 0.5).abs() < f64::EPSILON,
        "global machine queue assumption should survive merge"
    );
}

#[test]
fn execute_planning_rejects_schema_mismatch() {
    let mut engine = GentleEngine::default();
    let err = execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningProfileSet {
            scope: PlanningProfileScope::Global,
            payload_json: r#"{"schema":"gentle.planning_profile.v0"}"#.to_string(),
        },
    )
    .expect_err("profile schema mismatch must fail");
    assert!(err.contains("Unsupported planning profile schema"));

    let err = execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningObjectiveSet {
            payload_json: r#"{"schema":"gentle.planning_objective.v0"}"#.to_string(),
        },
    )
    .expect_err("objective schema mismatch must fail");
    assert!(err.contains("Unsupported planning objective schema"));
}

#[test]
fn execute_planning_sync_payload_validation_rejects_unknown_fields() {
    let mut engine = GentleEngine::default();
    let err = execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningSyncPull {
            payload_json: r#"{"unknown_patch":{"x":1}}"#.to_string(),
            source: None,
            confidence: None,
            snapshot_id: None,
        },
    )
    .expect_err("unknown sync payload fields must fail");
    assert!(err.contains("Invalid planning sync pull JSON payload"));
}

#[test]
fn execute_routines_list_without_planning_preserves_legacy_order() {
    let mut engine = GentleEngine::default();
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("routines_legacy_order.json");
    fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "zeta.workflow",
      "title": "Zeta Routine",
      "family": "zeta",
      "status": "implemented",
      "vocabulary_tags": ["zeta"],
      "template_name": "zeta_template",
      "base_time_hours": 1.0,
      "required_material_classes": ["missing_item"],
      "input_ports": [{ "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    },
    {
      "routine_id": "alpha.workflow",
      "title": "Alpha Routine",
      "family": "alpha",
      "status": "implemented",
      "vocabulary_tags": ["alpha"],
      "template_name": "alpha_template",
      "base_time_hours": 10.0,
      "input_ports": [{ "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    }
  ]
}"#,
        )
        .expect("write legacy-order catalog");

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RoutinesList {
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            family: None,
            status: None,
            tag: None,
            query: None,
            seq_id: None,
        },
    )
    .expect("list routines");
    assert_eq!(out.output["planning"]["enabled"].as_bool(), Some(false));
    let rows = out.output["routines"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    assert_eq!(rows.len(), 2);
    assert_eq!(
        rows[0].get("routine_id").and_then(|v| v.as_str()),
        Some("alpha.workflow"),
        "legacy non-planning order should remain family/title based"
    );
}

#[test]
fn execute_routines_list_applies_missing_material_penalty_default_business_days() {
    let mut engine = GentleEngine::default();
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("routines_planning_penalty.json");
    fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "restriction.missing_mix",
      "title": "Restriction with Missing Mix",
      "family": "restriction",
      "status": "implemented",
      "vocabulary_tags": ["restriction"],
      "template_name": "restriction_missing_mix",
      "base_time_hours": 4.0,
      "base_cost": 8.0,
      "required_material_classes": ["custom_mix", "custom_mix"],
      "input_ports": [{ "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    }
  ]
}"#,
        )
        .expect("write planning catalog");

    execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningProfileSet {
            scope: PlanningProfileScope::ProjectOverride,
            payload_json: r#"{
                  "schema":"gentle.planning_profile.v1",
                  "procurement_business_days_default":10,
                  "inventory":{"custom_mix":{"available":false}}
                }"#
            .to_string(),
        },
    )
    .expect("set planning profile");

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RoutinesList {
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            family: None,
            status: None,
            tag: None,
            query: None,
            seq_id: None,
        },
    )
    .expect("list routines");
    let rows = out.output["routines"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    assert_eq!(rows.len(), 1);
    let estimate = rows[0]["planning_estimate"].clone();
    let estimated_time_hours = estimate["estimated_time_hours"]
        .as_f64()
        .unwrap_or(f64::NAN);
    assert!(
        (estimated_time_hours - 340.0).abs() < 1e-9,
        "expected base 4h + 10 business days (weekend-aware elapsed-time) penalty"
    );
    let procurement_days = estimate["explanation"]["procurement_delay_business_days"]
        .as_f64()
        .unwrap_or(f64::NAN);
    assert!(
        (procurement_days - 10.0).abs() < 1e-9,
        "expected exactly one deduplicated missing material penalty"
    );
    let missing = estimate["explanation"]["missing_material_classes"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    assert_eq!(missing.len(), 1);
    assert_eq!(missing[0].as_str(), Some("custom_mix"));
}

#[test]
fn execute_routines_list_prefers_helper_compatible_routine_family() {
    let mut engine = GentleEngine::default();
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("routines_helper_preference.json");
    fs::write(
        &catalog_path,
        r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "gibson.demo_overlap",
      "title": "Gibson Demo Overlap",
      "family": "gibson",
      "status": "implemented",
      "vocabulary_tags": ["gibson", "assembly"],
      "template_name": "gibson_demo_overlap",
      "base_time_hours": 4.0,
      "base_cost": 2.0,
      "input_ports": [{ "port_id": "left_seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    },
    {
      "routine_id": "restriction.demo_subcloning",
      "title": "Restriction Demo Subcloning",
      "family": "restriction",
      "status": "implemented",
      "vocabulary_tags": ["restriction", "ligation"],
      "template_name": "restriction_demo_subcloning",
      "base_time_hours": 4.0,
      "base_cost": 2.0,
      "input_ports": [{ "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    }
  ]
}"#,
    )
    .expect("write planning catalog");

    execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningObjectiveSet {
            payload_json: r#"{
              "schema":"gentle.planning_objective.v1",
              "helper_profile_id":"pUC19"
            }"#
            .to_string(),
        },
    )
    .expect("set planning objective");

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RoutinesList {
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            family: None,
            status: None,
            tag: None,
            query: None,
            seq_id: None,
        },
    )
    .expect("list routines");
    let rows = out.output["routines"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    assert_eq!(rows.len(), 2);
    assert_eq!(
        rows[0]["routine_id"].as_str(),
        Some("restriction.demo_subcloning"),
        "helper-derived pUC19 preferences should boost restriction-family routines"
    );
    assert_eq!(
        rows[0]["planning_estimate"]["explanation"]["routine_family_alignment_bonus"].as_f64(),
        Some(0.12)
    );
    assert!(
        rows[0]["planning_estimate"]["explanation"]["routine_preference_context"]
            ["helper_derived_preferred_routine_families"]
            .as_array()
            .map(|rows| rows.iter().any(|row| row.as_str() == Some("restriction")))
            .unwrap_or(false)
    );
    assert_eq!(
        rows[1]["planning_estimate"]["explanation"]["routine_family_alignment_bonus"].as_f64(),
        Some(0.0)
    );
}

#[test]
fn execute_routines_list_prefers_variant_derived_routine_family_for_sequence() {
    let mut engine = GentleEngine::default();
    engine
        .set_planning_objective(Some(PlanningObjective::default()))
        .expect("set planning objective");
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "promoter".into(),
        location: gb_io::seq::Location::simple_range(0, 8),
        qualifiers: vec![("label".into(), Some("VKORC1 promoter".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "TFBS".into(),
        location: gb_io::seq::Location::simple_range(4, 7),
        qualifiers: vec![("label".into(), Some("SP1 motif".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "variation".into(),
        location: gb_io::seq::Location::simple_range(5, 6),
        qualifiers: vec![
            ("label".into(), Some("rsVariant".to_string())),
            (
                "gentle_generated".into(),
                Some("dbsnp_variant_marker".to_string()),
            ),
        ],
    });
    dna.update_computed_features();
    engine
        .state_mut()
        .sequences
        .insert("variant_routine_seq".to_string(), dna);

    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("routines_variant_preference.json");
    fs::write(
        &catalog_path,
        r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "pcr.demo_screen",
      "title": "PCR Demo Screen",
      "family": "pcr",
      "status": "implemented",
      "vocabulary_tags": ["pcr", "screen"],
      "template_name": "pcr_demo_screen",
      "base_time_hours": 4.0,
      "base_cost": 2.0,
      "input_ports": [{ "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    },
    {
      "routine_id": "gibson.demo_overlap",
      "title": "Gibson Demo Overlap",
      "family": "gibson",
      "status": "implemented",
      "vocabulary_tags": ["gibson", "assembly"],
      "template_name": "gibson_demo_overlap",
      "base_time_hours": 4.0,
      "base_cost": 2.0,
      "input_ports": [{ "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    }
  ]
}"#,
    )
    .expect("write planning catalog");

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::RoutinesList {
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            family: None,
            status: None,
            tag: None,
            query: None,
            seq_id: Some("variant_routine_seq".to_string()),
        },
    )
    .expect("list routines");
    let rows = out.output["routines"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    assert_eq!(rows.len(), 2);
    assert_eq!(
        rows[0]["routine_id"].as_str(),
        Some("gibson.demo_overlap"),
        "variant-derived reporter/expression preferences should boost compatible assembly families"
    );
    assert_eq!(
        rows[0]["planning_estimate"]["explanation"]["routine_family_alignment_bonus"].as_f64(),
        Some(0.16)
    );
    assert!(
        rows[0]["planning_estimate"]["explanation"]["routine_family_alignment_sources"]
            .as_array()
            .map(|rows| rows
                .iter()
                .any(|row| row.as_str() == Some("variant_derived")))
            .unwrap_or(false)
    );
    assert!(
        rows[0]["planning_estimate"]["explanation"]["routine_preference_context"]
            ["variant_suggested_assay_ids"]
            .as_array()
            .map(|rows| rows
                .iter()
                .any(|row| { row.as_str() == Some("allele_paired_promoter_luciferase_reporter") }))
            .unwrap_or(false)
    );
}

#[test]
fn execute_planning_suggestion_lifecycle_pending_accept_reject() {
    let mut engine = GentleEngine::default();

    let pull = execute_shell_command(
            &mut engine,
            &ShellCommand::PlanningSyncPull {
                payload_json: r#"{
                  "profile_patch":{"schema":"gentle.planning_profile.v1","capabilities":["pcr_machine"]},
                  "message":"sync pull"
                }"#
                .to_string(),
                source: Some("lab_agent".to_string()),
                confidence: Some(0.8),
                snapshot_id: Some("snap_pull_1".to_string()),
            },
        )
        .expect("sync pull suggestion");
    let suggestion_id = pull.output["suggestion"]["suggestion_id"]
        .as_str()
        .expect("pull suggestion id")
        .to_string();

    let pending = execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningSuggestionsList {
            status: Some(PlanningSuggestionStatus::Pending),
        },
    )
    .expect("list pending suggestions");
    assert_eq!(pending.output["suggestion_count"].as_u64(), Some(1));

    let accepted = execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningSuggestionAccept {
            suggestion_id: suggestion_id.clone(),
        },
    )
    .expect("accept suggestion");
    assert_eq!(accepted.output["status"].as_str(), Some("accepted"));
    assert!(
        engine
            .planning_effective_profile()
            .capabilities
            .iter()
            .any(|cap| cap == "pcr_machine"),
        "accepted suggestion should update confirmed overlay profile"
    );

    let push = execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningSyncPush {
            payload_json: r#"{
                  "objective_patch":{"schema":"gentle.planning_objective.v1","weight_time":2.5},
                  "message":"sync push"
                }"#
            .to_string(),
            source: Some("lab_agent".to_string()),
            confidence: Some(0.6),
            snapshot_id: Some("snap_push_1".to_string()),
        },
    )
    .expect("sync push suggestion");
    let reject_id = push.output["suggestion"]["suggestion_id"]
        .as_str()
        .expect("push suggestion id")
        .to_string();

    let rejected = execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningSuggestionReject {
            suggestion_id: reject_id.clone(),
            reason: Some("manual_override".to_string()),
        },
    )
    .expect("reject suggestion");
    assert_eq!(rejected.output["status"].as_str(), Some("rejected"));
    assert_eq!(
        rejected.output["suggestion"]["rejection_reason"].as_str(),
        Some("manual_override")
    );

    let accepted_rows = execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningSuggestionsList {
            status: Some(PlanningSuggestionStatus::Accepted),
        },
    )
    .expect("list accepted suggestions");
    assert_eq!(accepted_rows.output["suggestion_count"].as_u64(), Some(1));
    let rejected_rows = execute_shell_command(
        &mut engine,
        &ShellCommand::PlanningSuggestionsList {
            status: Some(PlanningSuggestionStatus::Rejected),
        },
    )
    .expect("list rejected suggestions");
    assert_eq!(rejected_rows.output["suggestion_count"].as_u64(), Some(1));
}

#[test]
fn execute_routines_explain_with_explicit_alternative_metadata() {
    let mut engine = GentleEngine::default();
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("cloning_routines_explain.json");
    fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "golden_gate.type_iis_single_insert",
      "title": "Golden Gate Type IIS Single Insert",
      "family": "golden_gate",
      "status": "implemented",
      "vocabulary_tags": ["golden_gate", "type_iis", "ligation"],
      "summary": "Type IIS one-pot assembly for one insert.",
      "purpose": "Assemble one insert into one destination vector with directional overhangs.",
      "mechanism": "Digest + ligate cycles with Type IIS enzymes and programmed overhang junctions.",
      "requires": ["Type IIS-compatible enzymes", "Defined non-conflicting overhang plan"],
      "contraindications": ["Ambiguous or duplicated internal overhangs"],
      "confusing_alternatives": ["gibson.two_fragment_overlap_preview"],
      "difference_matrix": [
        { "axis": "junction_constraint", "value": "explicit overhang grammar" },
        { "axis": "assembly_mode", "value": "restriction-ligation cycling" }
      ],
      "disambiguation_questions": ["Do you require Type IIS junction tokens?"],
      "failure_modes": ["duplicate junction overhang token"],
      "template_name": "golden_gate_single_insert",
      "input_ports": [
        { "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    },
    {
      "routine_id": "gibson.two_fragment_overlap_preview",
      "title": "Gibson Two-Fragment Overlap Preview",
      "family": "gibson",
      "status": "implemented",
      "vocabulary_tags": ["gibson", "overlap", "assembly"],
      "summary": "Overlap assembly preflight and preview path.",
      "template_name": "gibson_two_fragment_overlap_preview",
      "input_ports": [
        { "port_id": "left_seq_id", "kind": "sequence", "required": true, "cardinality": "one" },
        { "port_id": "right_seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    }
  ]
}"#,
        )
        .expect("write routines catalog");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::RoutinesExplain {
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            routine_id: "golden_gate.type_iis_single_insert".to_string(),
            seq_id: None,
        },
    )
    .expect("routines explain");
    assert!(!run.state_changed);
    assert_eq!(
        run.output["schema"].as_str(),
        Some(CLONING_ROUTINE_EXPLAIN_SCHEMA)
    );
    assert_eq!(
        run.output["routine"]["routine_id"].as_str(),
        Some("golden_gate.type_iis_single_insert")
    );
    let alternatives = run
        .output
        .get("alternatives")
        .and_then(|value| value.as_array())
        .cloned()
        .unwrap_or_default();
    assert_eq!(alternatives.len(), 1);
    assert_eq!(
        alternatives[0]
            .get("routine_id")
            .and_then(|value| value.as_str()),
        Some("gibson.two_fragment_overlap_preview")
    );
    let explanation_requires = run
        .output
        .get("explanation")
        .and_then(|value| value.get("requires"))
        .and_then(|value| value.as_array())
        .cloned()
        .unwrap_or_default();
    assert_eq!(explanation_requires.len(), 2);
}

#[test]
fn execute_routines_compare_returns_difference_matrix() {
    let mut engine = GentleEngine::default();
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("cloning_routines_compare.json");
    fs::write(
        &catalog_path,
        r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "golden_gate.type_iis_single_insert",
      "title": "Golden Gate Type IIS Single Insert",
      "family": "golden_gate",
      "status": "implemented",
      "vocabulary_tags": ["golden_gate", "assembly"],
      "summary": "Type IIS one-pot assembly for one insert.",
      "difference_matrix": [
        { "axis": "assembly_mode", "value": "restriction-ligation cycling" },
        { "axis": "junction_constraint", "value": "explicit overhang grammar" }
      ],
      "confusing_alternatives": ["gibson.two_fragment_overlap_preview"],
      "template_name": "golden_gate_single_insert",
      "input_ports": [
        { "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    },
    {
      "routine_id": "gibson.two_fragment_overlap_preview",
      "title": "Gibson Two-Fragment Overlap Preview",
      "family": "gibson",
      "status": "implemented",
      "vocabulary_tags": ["gibson", "assembly"],
      "summary": "Overlap assembly preflight and preview path.",
      "difference_matrix": [
        { "axis": "assembly_mode", "value": "homology-overlap assembly" },
        { "axis": "junction_constraint", "value": "homology overlap sequence identity" }
      ],
      "confusing_alternatives": ["golden_gate.type_iis_single_insert"],
      "template_name": "gibson_two_fragment_overlap_preview",
      "input_ports": [
        { "port_id": "left_seq_id", "kind": "sequence", "required": true, "cardinality": "one" },
        { "port_id": "right_seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    }
  ]
}"#,
    )
    .expect("write routines catalog");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::RoutinesCompare {
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            left_routine_id: "golden_gate.type_iis_single_insert".to_string(),
            right_routine_id: "gibson.two_fragment_overlap_preview".to_string(),
            seq_id: None,
        },
    )
    .expect("routines compare");
    assert!(!run.state_changed);
    assert_eq!(
        run.output["schema"].as_str(),
        Some(CLONING_ROUTINE_COMPARE_SCHEMA)
    );
    assert_eq!(
        run.output["comparison"]["cross_referenced_as_alternatives"].as_bool(),
        Some(true)
    );
    let matrix = run
        .output
        .get("comparison")
        .and_then(|value| value.get("difference_matrix"))
        .and_then(|value| value.as_array())
        .cloned()
        .unwrap_or_default();
    assert!(
        matrix
            .iter()
            .any(|row| row.get("axis").and_then(|value| value.as_str()) == Some("assembly_mode"))
    );
}

#[test]
fn execute_construct_reasoning_show_graph_includes_similarity_predictor_summary() {
    let sequence = format!(
        "{}{}{}{}{}",
        "ACGT".repeat(12),
        "AAAAAAAAAAAAAA",
        "ATATATATATATATATATAT",
        "GATTACAGATTACCCGGGGATTACAGATTA",
        "GCGTACGCTATTTTTAGCGTACGC"
    );
    let mut state = ProjectState::default();
    state.sequences.insert(
        "similarity_predictor_cli".to_string(),
        DNAsequence::from_sequence(&sequence).expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let graph = engine
        .build_construct_reasoning_graph(
            "similarity_predictor_cli",
            None,
            Some("similarity_predictor_cli_graph"),
        )
        .expect("build graph");

    let show = execute_shell_command(
        &mut engine,
        &ShellCommand::ConstructReasoningShowGraph {
            graph_id: graph.graph_id.clone(),
        },
    )
    .expect("show construct-reasoning graph");

    assert!(
        show.output["summary"]["summary_lines"]
            .as_array()
            .map(|rows| rows.iter().any(|row| {
                row.as_str().is_some_and(|text| {
                    text.contains("PCR/amplification review suggested")
                        && text.contains("review_needed")
                })
            }))
            .unwrap_or(false)
    );
    assert!(
        show.output["summary"]["summary_lines"]
            .as_array()
            .map(|rows| rows.iter().any(|row| {
                row.as_str().is_some_and(|text| {
                    text.contains("Nanopore/direct-sequencing review suggested")
                })
            }))
            .unwrap_or(false)
    );
    assert!(
        show.output["summary"]["fact_summaries"]
            .as_array()
            .map(|rows| rows.iter().any(|row| {
                row["fact_type"].as_str() == Some("mapping_operational_risk_context")
                    && row["detail_lines"]
                        .as_array()
                        .map(|detail_rows| {
                            detail_rows.iter().any(|entry| {
                                entry.as_str().is_some_and(|text| {
                                    text.contains("mapping_review: repeat_ambiguity=")
                                })
                            })
                        })
                        .unwrap_or(false)
            }))
            .unwrap_or(false)
    );
    assert!(
        show.output["summary"]["fact_summaries"]
            .as_array()
            .map(|rows| rows.iter().any(|row| {
                row["fact_type"].as_str() == Some("cloning_stability_context")
                    && row["warning_lines"]
                        .as_array()
                        .map(|warning_rows| {
                            warning_rows.iter().any(|entry| {
                                entry.as_str().is_some_and(|text| {
                                    text.contains("inversion/recombination review")
                                })
                            })
                        })
                        .unwrap_or(false)
            }))
            .unwrap_or(false)
    );
}

#[test]
fn execute_routines_explain_includes_variant_planning_context_for_sequence() {
    let mut engine = GentleEngine::default();
    engine
        .set_planning_objective(Some(PlanningObjective::default()))
        .expect("set planning objective");
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "promoter".into(),
        location: gb_io::seq::Location::simple_range(0, 8),
        qualifiers: vec![("label".into(), Some("VKORC1 promoter".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "TFBS".into(),
        location: gb_io::seq::Location::simple_range(4, 7),
        qualifiers: vec![("label".into(), Some("SP1 motif".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "variation".into(),
        location: gb_io::seq::Location::simple_range(5, 6),
        qualifiers: vec![
            ("label".into(), Some("rsExplain".to_string())),
            (
                "gentle_generated".into(),
                Some("dbsnp_variant_marker".to_string()),
            ),
        ],
    });
    dna.update_computed_features();
    engine
        .state_mut()
        .sequences
        .insert("variant_explain_seq".to_string(), dna);

    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("cloning_routines_explain_variant.json");
    fs::write(
        &catalog_path,
        r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "gibson.demo_overlap",
      "title": "Gibson Demo Overlap",
      "family": "gibson",
      "status": "implemented",
      "vocabulary_tags": ["gibson", "assembly"],
      "summary": "Assembly-compatible demo routine.",
      "template_name": "gibson_demo_overlap",
      "base_time_hours": 4.0,
      "base_cost": 2.0,
      "input_ports": [{ "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }],
      "output_ports": []
    }
  ]
}"#,
    )
    .expect("write routines catalog");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::RoutinesExplain {
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            routine_id: "gibson.demo_overlap".to_string(),
            seq_id: Some("variant_explain_seq".to_string()),
        },
    )
    .expect("routines explain");
    assert_eq!(
        run.output["planning"]["seq_id"].as_str(),
        Some("variant_explain_seq")
    );
    assert!(
        run.output["planning"]["routine_preference_context"]["variant_effect_tags"]
            .as_array()
            .map(|rows| {
                rows.iter()
                    .any(|row| row.as_str() == Some("promoter_tfbs_regulatory_candidate"))
            })
            .unwrap_or(false)
    );
    assert!(
        run.output["planning"]["estimate"]["explanation"]["routine_family_alignment_sources"]
            .as_array()
            .map(|rows| {
                rows.iter()
                    .any(|row| row.as_str() == Some("variant_derived"))
            })
            .unwrap_or(false)
    );
}

#[test]
fn load_cloning_routine_catalog_rejects_unknown_confusing_alternative() {
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("cloning_routines_bad_alt.json");
    fs::write(
        &catalog_path,
        r#"{
  "schema": "gentle.cloning_routines.v1",
  "routines": [
    {
      "routine_id": "golden_gate.type_iis_single_insert",
      "title": "Golden Gate Type IIS Single Insert",
      "family": "golden_gate",
      "status": "implemented",
      "vocabulary_tags": ["golden_gate"],
      "summary": "baseline",
      "confusing_alternatives": ["missing.routine.id"],
      "template_name": "golden_gate_single_insert",
      "input_ports": [
        { "port_id": "seq_id", "kind": "sequence", "required": true, "cardinality": "one" }
      ],
      "output_ports": [
        { "port_id": "output_id", "kind": "sequence", "required": false, "cardinality": "one" }
      ]
    }
  ]
}"#,
    )
    .expect("write routines catalog");

    let err = load_cloning_routine_catalog(catalog_path.to_string_lossy().as_ref())
        .expect_err("unknown alternatives must be rejected");
    assert!(err.contains("unknown confusing_alternative"));
}

#[test]
fn execute_macros_template_registry_and_run() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateUpsert {
            name: "clone".to_string(),
            description: Some("reverse helper".to_string()),
            details_url: Some("https://example.org/macros/clone".to_string()),
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
            script: r#"op {"Reverse":{"input":"${seq_id}","output_id":"${out_id}"}}"#.to_string(),
        },
    )
    .expect("upsert macros template");

    let listed = execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateList)
        .expect("list macros templates");
    assert_eq!(listed.output["template_count"].as_u64(), Some(1));
    let details_url = listed
        .output
        .get("templates")
        .and_then(|v| v.as_array())
        .and_then(|arr| arr.first())
        .and_then(|template| template.get("details_url"))
        .and_then(|v| v.as_str());
    assert_eq!(details_url, Some("https://example.org/macros/clone"));

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "clone".to_string(),
            bindings: HashMap::from([("seq_id".to_string(), "seqA".to_string())]),
            transactional: false,
            validate_only: false,
        },
    )
    .expect("run macros template");
    assert!(run.state_changed);
    assert_eq!(run.output["template_name"].as_str(), Some("clone"));
    assert!(engine.state().sequences.contains_key("seqA_rev"));
}

#[test]
fn execute_macros_template_run_records_macro_instance() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateUpsert {
            name: "clone".to_string(),
            description: Some("reverse helper".to_string()),
            details_url: Some("https://example.org/macros/clone".to_string()),
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
            input_ports: vec![WorkflowMacroTemplatePort {
                port_id: "seq_id".to_string(),
                kind: "sequence".to_string(),
                required: true,
                cardinality: "one".to_string(),
                description: Some("Template input".to_string()),
            }],
            output_ports: vec![WorkflowMacroTemplatePort {
                port_id: "out_id".to_string(),
                kind: "sequence".to_string(),
                required: false,
                cardinality: "one".to_string(),
                description: Some("Derived sequence id".to_string()),
            }],
            script: r#"op {"Reverse":{"input":"${seq_id}","output_id":"${out_id}"}}"#.to_string(),
        },
    )
    .expect("upsert template");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "clone".to_string(),
            bindings: HashMap::from([
                ("seq_id".to_string(), "seqA".to_string()),
                ("out_id".to_string(), "seqA_rev".to_string()),
            ]),
            transactional: false,
            validate_only: false,
        },
    )
    .expect("run template");
    assert!(run.state_changed);
    assert!(run.output["macro_instance_id"].as_str().is_some());
    assert_eq!(engine.state().lineage.macro_instances.len(), 1);
    let instance = &engine.state().lineage.macro_instances[0];
    assert_eq!(instance.template_name.as_deref(), Some("clone"));
    assert!(!instance.expanded_op_ids.is_empty());
    assert_eq!(instance.bound_inputs.len(), 1);
    assert_eq!(instance.bound_inputs[0].port_id, "seq_id");
    assert_eq!(instance.bound_outputs.len(), 1);
    assert_eq!(instance.bound_outputs[0].port_id, "out_id");
}

#[test]
fn execute_macros_template_cloning_digest_ligation_extract_fixture() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "x".to_string(),
        DNAsequence::from_sequence("ATGGATCCGCATGGATCCGCATGGATCCGC").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let put = parse_shell_line(
            "macros template-put clone_slice --file test_files/cloning_digest_ligation_extract.gsh --param seq_id=x --param digest_prefix=d --param ligation_prefix=lig --param extract_from=0 --param extract_to=1 --param output_id=slice",
        )
        .expect("parse macros template-put clone fixture");
    let upsert = execute_shell_command(&mut engine, &put).expect("upsert clone fixture template");
    assert!(upsert.state_changed);

    let run_cmd = parse_shell_line("macros template-run clone_slice --transactional")
        .expect("parse macros template-run clone fixture");
    let run = execute_shell_command(&mut engine, &run_cmd).expect("run clone fixture template");
    assert!(run.state_changed);
    assert_eq!(run.output["template_name"].as_str(), Some("clone_slice"));
    assert!(engine.state().sequences.contains_key("d_1"));
    assert!(engine.state().sequences.contains_key("lig_1"));
    let slice = engine
        .state()
        .sequences
        .get("slice")
        .expect("expected extracted slice output");
    assert_eq!(slice.len(), 1);
}

#[test]
fn execute_macros_template_import_file_and_run() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let tmp = tempdir().expect("tempdir");
    let pattern_path = tmp.path().join("patterns.json");
    fs::write(
        &pattern_path,
        r#"{
  "schema": "gentle.cloning_patterns.v1",
  "templates": [
    {
      "name": "reverse_default",
      "description": "reverse one sequence",
      "details_url": "https://example.org/reverse-default",
      "parameters": [
        { "name": "seq_id", "default_value": "seqA", "required": false },
        { "name": "output_id", "default_value": "seqA_rev", "required": false }
      ],
      "script": "op {\"Reverse\":{\"input\":\"${seq_id}\",\"output_id\":\"${output_id}\"}}"
    }
  ]
}"#,
    )
    .expect("write patterns file");

    let import = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateImport {
            path: pattern_path.to_string_lossy().to_string(),
        },
    )
    .expect("import patterns");
    assert!(import.state_changed);
    assert_eq!(import.output["imported_count"].as_u64(), Some(1));
    let imported_template = engine
        .get_workflow_macro_template("reverse_default")
        .expect("get imported template");
    assert_eq!(
        imported_template.details_url.as_deref(),
        Some("https://example.org/reverse-default")
    );

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "reverse_default".to_string(),
            bindings: HashMap::new(),
            transactional: false,
            validate_only: false,
        },
    )
    .expect("run imported template");
    assert!(run.state_changed);
    assert!(engine.state().sequences.contains_key("seqA_rev"));
}

#[test]
fn execute_macros_template_import_single_template_schema_file_and_run() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let tmp = tempdir().expect("tempdir");
    let pattern_path = tmp.path().join("reverse_one.json");
    fs::write(
        &pattern_path,
        r#"{
  "schema": "gentle.cloning_pattern_template.v1",
  "name": "reverse_one",
  "description": "reverse one sequence",
  "details_url": "https://example.org/reverse-one",
  "parameters": [
    { "name": "seq_id", "default_value": "seqA", "required": false },
    { "name": "output_id", "default_value": "seqA_rev", "required": false }
  ],
  "script": "op {\"Reverse\":{\"input\":\"${seq_id}\",\"output_id\":\"${output_id}\"}}"
}"#,
    )
    .expect("write single-template file");

    let import = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateImport {
            path: pattern_path.to_string_lossy().to_string(),
        },
    )
    .expect("import single-template file");
    assert!(import.state_changed);
    assert_eq!(import.output["imported_count"].as_u64(), Some(1));
    assert_eq!(
        import
            .output
            .get("source_files")
            .and_then(|v| v.as_array())
            .map(|rows| rows.len()),
        Some(1)
    );

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "reverse_one".to_string(),
            bindings: HashMap::new(),
            transactional: false,
            validate_only: false,
        },
    )
    .expect("run imported single-template");
    assert!(run.state_changed);
    assert!(engine.state().sequences.contains_key("seqA_rev"));
}

#[test]
fn execute_macros_template_import_directory_recursive() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let tmp = tempdir().expect("tempdir");
    let root = tmp.path().join("catalog");
    fs::create_dir_all(root.join("a")).expect("mkdir a");
    fs::create_dir_all(root.join("b").join("nested")).expect("mkdir nested");
    fs::write(
        root.join("a").join("reverse_one.json"),
        r#"{
  "schema": "gentle.cloning_pattern_template.v1",
  "name": "reverse_one",
  "parameters": [
    { "name": "seq_id", "default_value": "seqA", "required": false },
    { "name": "output_id", "default_value": "seqA_rev", "required": false }
  ],
  "script": "op {\"Reverse\":{\"input\":\"${seq_id}\",\"output_id\":\"${output_id}\"}}"
}"#,
    )
    .expect("write reverse");
    fs::write(
        root.join("b").join("nested").join("complement_one.json"),
        r#"{
  "schema": "gentle.cloning_pattern_template.v1",
  "name": "complement_one",
  "parameters": [
    { "name": "seq_id", "default_value": "seqA", "required": false },
    { "name": "output_id", "default_value": "seqA_comp", "required": false }
  ],
  "script": "op {\"Complement\":{\"input\":\"${seq_id}\",\"output_id\":\"${output_id}\"}}"
}"#,
    )
    .expect("write complement");

    let import = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateImport {
            path: root.to_string_lossy().to_string(),
        },
    )
    .expect("import catalog directory");
    assert!(import.state_changed);
    assert_eq!(import.output["imported_count"].as_u64(), Some(2));
    assert_eq!(
        import
            .output
            .get("source_files")
            .and_then(|v| v.as_array())
            .map(|rows| rows.len()),
        Some(2)
    );
    let templates = import
        .output
        .get("templates")
        .and_then(|v| v.as_array())
        .cloned()
        .unwrap_or_default();
    assert!(templates.iter().any(|v| v.as_str() == Some("reverse_one")));
    assert!(
        templates
            .iter()
            .any(|v| v.as_str() == Some("complement_one"))
    );
}

#[test]
fn execute_macros_template_import_builtin_patterns_and_run_template() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let path = format!(
        "{}/assets/cloning_patterns.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let import = execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
        .expect("import built-in patterns");
    assert!(import.state_changed);
    assert!(import.output["imported_count"].as_u64().unwrap_or(0) >= 6);
    let templates = import
        .output
        .get("templates")
        .and_then(|v| v.as_array())
        .cloned()
        .unwrap_or_default();
    assert!(templates.iter().any(|v| {
        v.as_str()
            .map(|s| s == "grna_candidate_priority_scan")
            .unwrap_or(false)
    }));
    let branch_template = engine
        .get_workflow_macro_template("branch_reverse_complement")
        .expect("imported branch template should exist");
    assert_eq!(
        branch_template.details_url.as_deref(),
        Some("https://www.bioinformatics.org/sms/rev_comp.html")
    );

    let mut bindings = HashMap::new();
    bindings.insert("seq_id".to_string(), "seqA".to_string());
    bindings.insert("branch_copy_id".to_string(), "seqA_branch".to_string());
    bindings.insert(
        "reverse_complement_id".to_string(),
        "seqA_branch_rc".to_string(),
    );
    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "branch_reverse_complement".to_string(),
            bindings,
            transactional: false,
            validate_only: false,
        },
    )
    .expect("run imported branch template");
    assert!(run.state_changed);
    assert!(engine.state().sequences.contains_key("seqA_branch"));
    assert!(engine.state().sequences.contains_key("seqA_branch_rc"));
}

#[test]
fn execute_macros_template_run_validate_only_reports_preflight_without_mutation() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let path = format!(
        "{}/assets/cloning_patterns_catalog/sequence/transform/branch_reverse_complement.json",
        env!("CARGO_MANIFEST_DIR")
    );
    execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
        .expect("import branch template");
    assert!(!engine.state().sequences.contains_key("seqA_branch"));

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "branch_reverse_complement".to_string(),
            bindings: HashMap::new(),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("validate-only run");
    assert!(!run.state_changed);
    assert_eq!(run.output["validate_only"].as_bool(), Some(true));
    assert_eq!(run.output["can_execute"].as_bool(), Some(false));
    assert!(
        run.output
            .get("preflight")
            .and_then(|v| v.get("errors"))
            .and_then(|v| v.as_array())
            .map(|rows| !rows.is_empty())
            .unwrap_or(false)
    );
    assert!(!engine.state().sequences.contains_key("seqA_branch"));
}

#[test]
fn execute_macros_template_run_preflight_failure_records_macro_instance() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let path = format!(
        "{}/assets/cloning_patterns_catalog/sequence/transform/branch_reverse_complement.json",
        env!("CARGO_MANIFEST_DIR")
    );
    execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
        .expect("import branch template");
    let err = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "branch_reverse_complement".to_string(),
            bindings: HashMap::new(),
            transactional: false,
            validate_only: false,
        },
    )
    .expect_err("preflight should fail");
    assert!(err.contains("macro_instance_id="));
    assert_eq!(engine.state().lineage.macro_instances.len(), 1);
    let instance = &engine.state().lineage.macro_instances[0];
    assert_eq!(instance.status, MacroInstanceStatus::Failed);
    assert_eq!(
        instance.template_name.as_deref(),
        Some("branch_reverse_complement")
    );
    assert!(instance.status_message.is_some());
}

#[test]
fn execute_macros_template_validate_only_checks_anchor_semantics() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateUpsert {
            name: "anchor_semantics".to_string(),
            description: Some("anchor semantics test".to_string()),
            details_url: None,
            parameters: vec![
                WorkflowMacroTemplateParam {
                    name: "seq_id".to_string(),
                    default_value: None,
                    required: true,
                },
                WorkflowMacroTemplateParam {
                    name: "anchor_a".to_string(),
                    default_value: None,
                    required: true,
                },
            ],
            input_ports: vec![
                WorkflowMacroTemplatePort {
                    port_id: "seq_id".to_string(),
                    kind: "sequence".to_string(),
                    required: true,
                    cardinality: "one".to_string(),
                    description: Some("sequence input".to_string()),
                },
                WorkflowMacroTemplatePort {
                    port_id: "anchor_a".to_string(),
                    kind: "sequence_anchor".to_string(),
                    required: true,
                    cardinality: "one".to_string(),
                    description: Some("anchor input".to_string()),
                },
            ],
            output_ports: vec![],
            script: r#"op {"Reverse":{"input":"${seq_id}"}}"#.to_string(),
        },
    )
    .expect("upsert template");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "anchor_semantics".to_string(),
            bindings: HashMap::from([
                ("seq_id".to_string(), "seqA".to_string()),
                ("anchor_a".to_string(), "999".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("validate-only");
    assert!(!run.state_changed);
    assert_eq!(run.output["can_execute"].as_bool(), Some(false));
    let preflight_errors = run
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|row| row.as_str().map(|s| s.to_string()))
        .collect::<Vec<_>>();
    assert!(
        preflight_errors
            .iter()
            .any(|message| message.contains("anchor check failed")),
        "expected semantic anchor validation error, got: {:?}",
        preflight_errors
    );
}

#[test]
fn execute_macros_template_validate_only_rejects_output_alias_collisions() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateUpsert {
            name: "alias_collision".to_string(),
            description: Some("output alias collision".to_string()),
            details_url: None,
            parameters: vec![
                WorkflowMacroTemplateParam {
                    name: "seq_id".to_string(),
                    default_value: None,
                    required: true,
                },
                WorkflowMacroTemplateParam {
                    name: "out_a".to_string(),
                    default_value: None,
                    required: true,
                },
                WorkflowMacroTemplateParam {
                    name: "out_b".to_string(),
                    default_value: None,
                    required: true,
                },
            ],
            input_ports: vec![WorkflowMacroTemplatePort {
                port_id: "seq_id".to_string(),
                kind: "sequence".to_string(),
                required: true,
                cardinality: "one".to_string(),
                description: None,
            }],
            output_ports: vec![
                WorkflowMacroTemplatePort {
                    port_id: "out_a".to_string(),
                    kind: "sequence".to_string(),
                    required: true,
                    cardinality: "one".to_string(),
                    description: None,
                },
                WorkflowMacroTemplatePort {
                    port_id: "out_b".to_string(),
                    kind: "sequence".to_string(),
                    required: true,
                    cardinality: "one".to_string(),
                    description: None,
                },
            ],
            script: r#"op {"Reverse":{"input":"${seq_id}","output_id":"${out_a}"}}"#.to_string(),
        },
    )
    .expect("upsert template");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "alias_collision".to_string(),
            bindings: HashMap::from([
                ("seq_id".to_string(), "seqA".to_string()),
                ("out_a".to_string(), "same_out".to_string()),
                ("out_b".to_string(), "same_out".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("validate-only");
    assert_eq!(run.output["can_execute"].as_bool(), Some(false));
    let preflight_errors = run
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|row| row.as_str().map(|s| s.to_string()))
        .collect::<Vec<_>>();
    assert!(
        preflight_errors
            .iter()
            .any(|message| message.contains("Output alias conflict")),
        "expected output alias conflict error, got: {:?}",
        preflight_errors
    );
}

#[test]
fn execute_macros_template_validate_only_applies_gibson_family_overlap_checks() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "left".to_string(),
        DNAsequence::from_sequence("ATGCCGTTACCGGTTAAACCCGGGTTT").expect("sequence"),
    );
    state.sequences.insert(
        "right_ok".to_string(),
        DNAsequence::from_sequence("AACCCGGGTTTGGGAAATTTCCCGGG").expect("sequence"),
    );
    state.sequences.insert(
        "right_bad".to_string(),
        DNAsequence::from_sequence("TTTAAACCCGGGGGGAAATTTCCCGGG").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let path = format!(
        "{}/assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json",
        env!("CARGO_MANIFEST_DIR")
    );
    execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
        .expect("import gibson template");

    let ok_run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "gibson_two_fragment_overlap_preview".to_string(),
            bindings: HashMap::from([
                ("left_seq_id".to_string(), "left".to_string()),
                ("right_seq_id".to_string(), "right_ok".to_string()),
                ("overlap_bp".to_string(), "11".to_string()),
                ("assembly_prefix".to_string(), "gib_ok".to_string()),
                ("output_id".to_string(), "gib_ok_forward".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("gibson validate-only success");
    assert_eq!(ok_run.output["can_execute"].as_bool(), Some(true));
    let ok_errors = ok_run
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default();
    assert!(
        ok_errors.is_empty(),
        "did not expect Gibson preflight errors for matching overlap: {:?}",
        ok_errors
    );

    let bad_run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "gibson_two_fragment_overlap_preview".to_string(),
            bindings: HashMap::from([
                ("left_seq_id".to_string(), "left".to_string()),
                ("right_seq_id".to_string(), "right_bad".to_string()),
                ("overlap_bp".to_string(), "11".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("gibson validate-only mismatch");
    assert_eq!(bad_run.output["can_execute"].as_bool(), Some(false));
    let bad_errors = bad_run
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|row| row.as_str().map(|s| s.to_string()))
        .collect::<Vec<_>>();
    assert!(
        bad_errors
            .iter()
            .any(|message| message.contains("Gibson overlap mismatch")),
        "expected Gibson overlap mismatch error, got: {:?}",
        bad_errors
    );
}

#[test]
fn execute_macros_template_validate_only_applies_restriction_family_checks() {
    let enzymes = crate::enzymes::active_restriction_enzymes();
    assert!(
        enzymes.len() >= 2,
        "restriction enzyme catalog should provide at least two enzymes"
    );
    let enzyme_a = enzymes
        .iter()
        .find(|enzyme| !enzyme.sequence.is_empty())
        .expect("enzyme_a with sequence");
    let enzyme_b = enzymes
        .iter()
        .find(|enzyme| enzyme.name != enzyme_a.name && !enzyme.sequence.is_empty())
        .expect("enzyme_b with sequence");
    let sequence_text = format!("AAA{}TTT{}GGG", enzyme_a.sequence, enzyme_b.sequence);

    let mut state = ProjectState::default();
    state.sequences.insert(
        "vector".to_string(),
        DNAsequence::from_sequence(&sequence_text).expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let path = format!(
        "{}/assets/cloning_patterns_catalog/restriction/digest_ligation/digest_ligate_extract_sticky.json",
        env!("CARGO_MANIFEST_DIR")
    );
    execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
        .expect("import restriction template");

    let ok_run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "digest_ligate_extract_sticky".to_string(),
            bindings: HashMap::from([
                ("seq_id".to_string(), "vector".to_string()),
                ("enzyme_a".to_string(), enzyme_a.name.clone()),
                ("enzyme_b".to_string(), enzyme_b.name.clone()),
                ("left_fragment".to_string(), "1".to_string()),
                ("right_fragment".to_string(), "2".to_string()),
                ("extract_from".to_string(), "0".to_string()),
                ("extract_to".to_string(), "20".to_string()),
                ("digest_prefix".to_string(), "digest_ok".to_string()),
                ("ligation_prefix".to_string(), "lig_ok".to_string()),
                ("output_id".to_string(), "restriction_ok".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("restriction validate-only success");
    assert_eq!(ok_run.output["can_execute"].as_bool(), Some(true));
    let ok_errors = ok_run
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default();
    assert!(
        ok_errors.is_empty(),
        "did not expect restriction preflight errors for valid inputs: {:?}",
        ok_errors
    );

    let duplicate_run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "digest_ligate_extract_sticky".to_string(),
            bindings: HashMap::from([
                ("seq_id".to_string(), "vector".to_string()),
                ("enzyme_a".to_string(), enzyme_a.name.clone()),
                ("enzyme_b".to_string(), enzyme_a.name.clone()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("restriction validate-only duplicate enzyme");
    assert_eq!(duplicate_run.output["can_execute"].as_bool(), Some(false));
    let duplicate_errors = duplicate_run
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|row| row.as_str().map(|s| s.to_string()))
        .collect::<Vec<_>>();
    assert!(
        duplicate_errors
            .iter()
            .any(|message| message.contains("expects distinct enzyme inputs")),
        "expected duplicate-enzyme restriction error, got: {:?}",
        duplicate_errors
    );

    let unknown_run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "digest_ligate_extract_sticky".to_string(),
            bindings: HashMap::from([
                ("seq_id".to_string(), "vector".to_string()),
                ("enzyme_a".to_string(), "__missing_enzyme__".to_string()),
                ("enzyme_b".to_string(), enzyme_b.name.clone()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("restriction validate-only unknown enzyme");
    assert_eq!(unknown_run.output["can_execute"].as_bool(), Some(false));
    let unknown_errors = unknown_run
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|row| row.as_str().map(|s| s.to_string()))
        .collect::<Vec<_>>();
    assert!(
        unknown_errors
            .iter()
            .any(|message| message.contains("unknown enzyme name")),
        "expected unknown-enzyme restriction error, got: {:?}",
        unknown_errors
    );
}

#[test]
fn execute_macros_template_validate_only_applies_golden_gate_family_checks() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "gg_vector".to_string(),
        DNAsequence::from_sequence("AAAAGGTCTCTTTTGGTCTCAAAA").expect("sequence"),
    );
    state.sequences.insert(
        "gg_insert".to_string(),
        DNAsequence::from_sequence("CCCCGGTCTCGGGG").expect("sequence"),
    );
    state.sequences.insert(
        "plain_seq".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let base = env!("CARGO_MANIFEST_DIR");
    let single_path = format!(
        "{}/assets/cloning_patterns_catalog/golden_gate/type_iis/golden_gate_single_insert.json",
        base
    );
    let multi_path = format!(
        "{}/assets/cloning_patterns_catalog/golden_gate/type_iis/golden_gate_multi_insert.json",
        base
    );
    execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateImport { path: single_path },
    )
    .expect("import golden gate single");
    execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateImport { path: multi_path },
    )
    .expect("import golden gate multi");

    let ok_run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "golden_gate_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "gg_vector".to_string()),
                ("insert_seq_id".to_string(), "gg_insert".to_string()),
                ("type_iis_enzyme".to_string(), "Eco31".to_string()),
                ("junction_overhangs".to_string(), "AATG".to_string()),
                ("vector_fragment".to_string(), "1".to_string()),
                ("insert_fragment".to_string(), "1".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("golden gate validate-only success");
    assert_eq!(ok_run.output["can_execute"].as_bool(), Some(true));

    let unknown_or_non_type_iis = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "golden_gate_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "gg_vector".to_string()),
                ("insert_seq_id".to_string(), "gg_insert".to_string()),
                ("type_iis_enzyme".to_string(), "EcoRI".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("golden gate validate-only non-type-iis");
    assert_eq!(
        unknown_or_non_type_iis.output["can_execute"].as_bool(),
        Some(false)
    );
    let non_type_iis_errors = unknown_or_non_type_iis
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|row| row.as_str().map(|s| s.to_string()))
        .collect::<Vec<_>>();
    assert!(
        non_type_iis_errors
            .iter()
            .any(|message| message.contains("Type IIS-capable")),
        "expected Type IIS-capable error, got: {:?}",
        non_type_iis_errors
    );

    let missing_site = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "golden_gate_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "plain_seq".to_string()),
                ("insert_seq_id".to_string(), "gg_insert".to_string()),
                ("type_iis_enzyme".to_string(), "Eco31".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("golden gate validate-only missing site");
    assert_eq!(missing_site.output["can_execute"].as_bool(), Some(false));
    let missing_site_errors = missing_site
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|row| row.as_str().map(|s| s.to_string()))
        .collect::<Vec<_>>();
    assert!(
        missing_site_errors
            .iter()
            .any(|message| message.contains("no recognition site")),
        "expected missing-site error, got: {:?}",
        missing_site_errors
    );

    let invalid_junction_or_fragment = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "golden_gate_multi_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "gg_vector".to_string()),
                ("insert_a_seq_id".to_string(), "gg_insert".to_string()),
                ("insert_b_seq_id".to_string(), "gg_insert".to_string()),
                ("type_iis_enzyme".to_string(), "Eco31".to_string()),
                ("junction_overhangs".to_string(), "AATG,,GCTT".to_string()),
                ("vector_fragment".to_string(), "0".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("golden gate validate-only invalid overhang/fragment");
    assert_eq!(
        invalid_junction_or_fragment.output["can_execute"].as_bool(),
        Some(false)
    );
    let invalid_errors = invalid_junction_or_fragment
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|row| row.as_str().map(|s| s.to_string()))
        .collect::<Vec<_>>();
    assert!(
        invalid_errors.iter().any(|message| {
            message.contains("contains an empty token")
                || message.contains("expects vector_fragment >= 1")
        }),
        "expected invalid-junction/fragment errors, got: {:?}",
        invalid_errors
    );
}

#[test]
fn execute_macros_template_validate_only_applies_gateway_topo_and_ta_gc_checks() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "gw_donor".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGA").expect("sequence"),
    );
    state.sequences.insert(
        "gw_entry".to_string(),
        DNAsequence::from_sequence("TTTTACGTACGT").expect("sequence"),
    );
    state.sequences.insert(
        "topo_vector_t".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    state.sequences.insert(
        "topo_insert_a".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGA").expect("sequence"),
    );
    state.sequences.insert(
        "topo_insert_bad".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGG").expect("sequence"),
    );
    state.sequences.insert(
        "topo_insert_cacc".to_string(),
        DNAsequence::from_sequence("CACCACGTACGTACGA").expect("sequence"),
    );
    state.sequences.insert(
        "gc_vector_c".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGC").expect("sequence"),
    );
    state.sequences.insert(
        "gc_insert_g".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGG").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let base = env!("CARGO_MANIFEST_DIR");
    for path in [
        format!(
            "{}/assets/cloning_patterns_catalog/gateway/bp_lr/gateway_bp_single_insert.json",
            base
        ),
        format!(
            "{}/assets/cloning_patterns_catalog/topo/entry/topo_ta_single_insert.json",
            base
        ),
        format!(
            "{}/assets/cloning_patterns_catalog/topo/entry/topo_directional_cacc_single_insert.json",
            base
        ),
        format!(
            "{}/assets/cloning_patterns_catalog/ta_gc/entry/ta_clone_single_insert.json",
            base
        ),
        format!(
            "{}/assets/cloning_patterns_catalog/ta_gc/entry/gc_clone_single_insert.json",
            base
        ),
    ] {
        execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
            .expect("import family template");
    }

    let gateway_ok = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "gateway_bp_single_insert".to_string(),
            bindings: HashMap::from([
                ("donor_seq_id".to_string(), "gw_donor".to_string()),
                ("entry_vector_seq_id".to_string(), "gw_entry".to_string()),
                ("gateway_phase".to_string(), "bp".to_string()),
                ("att_tokens".to_string(), "attB,attP".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("gateway bp validate-only ok");
    assert_eq!(gateway_ok.output["can_execute"].as_bool(), Some(true));

    let gateway_bad = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "gateway_bp_single_insert".to_string(),
            bindings: HashMap::from([
                ("donor_seq_id".to_string(), "gw_donor".to_string()),
                ("entry_vector_seq_id".to_string(), "gw_entry".to_string()),
                ("gateway_phase".to_string(), "bp".to_string()),
                ("att_tokens".to_string(), "attL,attR".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("gateway bp validate-only bad");
    assert_eq!(gateway_bad.output["can_execute"].as_bool(), Some(false));

    let topo_ta_ok = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "topo_ta_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                ("insert_seq_id".to_string(), "topo_insert_a".to_string()),
                ("topo_mode".to_string(), "ta".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("topo ta validate-only ok");
    assert_eq!(topo_ta_ok.output["can_execute"].as_bool(), Some(true));

    let topo_dir_bad = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "topo_directional_cacc_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                ("insert_seq_id".to_string(), "topo_insert_bad".to_string()),
                ("topo_mode".to_string(), "directional_cacc".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("topo directional bad");
    assert_eq!(topo_dir_bad.output["can_execute"].as_bool(), Some(false));

    let ta_ok = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "ta_clone_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                ("insert_seq_id".to_string(), "topo_insert_a".to_string()),
                ("tail_mode".to_string(), "ta".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("ta clone validate-only ok");
    assert_eq!(ta_ok.output["can_execute"].as_bool(), Some(true));

    let gc_bad = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "gc_clone_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                ("insert_seq_id".to_string(), "topo_insert_a".to_string()),
                ("tail_mode".to_string(), "gc".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("gc clone bad");
    assert_eq!(gc_bad.output["can_execute"].as_bool(), Some(false));

    let gc_ok = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "gc_clone_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "gc_vector_c".to_string()),
                ("insert_seq_id".to_string(), "gc_insert_g".to_string()),
                ("tail_mode".to_string(), "gc".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("gc clone ok");
    assert_eq!(gc_ok.output["can_execute"].as_bool(), Some(true));

    let topo_dir_ok = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "topo_directional_cacc_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                ("insert_seq_id".to_string(), "topo_insert_cacc".to_string()),
                ("topo_mode".to_string(), "directional_cacc".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("topo directional ok");
    assert_eq!(topo_dir_ok.output["can_execute"].as_bool(), Some(true));
}

#[test]
fn execute_macros_template_validate_only_applies_infusion_and_nebuilder_checks() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "left".to_string(),
        DNAsequence::from_sequence("ACACACGGGGAAAA").expect("sequence"),
    );
    state.sequences.insert(
        "middle".to_string(),
        DNAsequence::from_sequence("GGGGAAAATTTTCCCC").expect("sequence"),
    );
    state.sequences.insert(
        "right".to_string(),
        DNAsequence::from_sequence("TTTTCCCCGAGAGA").expect("sequence"),
    );
    state.sequences.insert(
        "right_bad".to_string(),
        DNAsequence::from_sequence("CCCCAAAAGAGAGA").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let base = env!("CARGO_MANIFEST_DIR");
    for path in [
        format!(
            "{}/assets/cloning_patterns_catalog/infusion/overlap/infusion_two_fragment_overlap.json",
            base
        ),
        format!(
            "{}/assets/cloning_patterns_catalog/nebuilder_hifi/overlap/nebuilder_multi_fragment_overlap.json",
            base
        ),
    ] {
        execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
            .expect("import overlap template");
    }

    let infusion_ok = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "infusion_two_fragment_overlap".to_string(),
            bindings: HashMap::from([
                ("left_seq_id".to_string(), "left".to_string()),
                ("right_seq_id".to_string(), "middle".to_string()),
                ("overlap_bp".to_string(), "8".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("infusion validate-only ok");
    assert_eq!(infusion_ok.output["can_execute"].as_bool(), Some(true));

    let infusion_bad = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "infusion_two_fragment_overlap".to_string(),
            bindings: HashMap::from([
                ("left_seq_id".to_string(), "left".to_string()),
                ("right_seq_id".to_string(), "right_bad".to_string()),
                ("overlap_bp".to_string(), "8".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("infusion validate-only bad");
    assert_eq!(infusion_bad.output["can_execute"].as_bool(), Some(false));
    let infusion_bad_errors = infusion_bad
        .output
        .get("preflight")
        .and_then(|preflight| preflight.get("errors"))
        .and_then(|rows| rows.as_array())
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|row| row.as_str().map(|s| s.to_string()))
        .collect::<Vec<_>>();
    assert!(
        infusion_bad_errors
            .iter()
            .any(|message| message.contains("In-Fusion overlap mismatch")),
        "expected In-Fusion overlap mismatch error, got: {:?}",
        infusion_bad_errors
    );

    let nebuilder_ok = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "nebuilder_multi_fragment_overlap".to_string(),
            bindings: HashMap::from([
                ("left_seq_id".to_string(), "left".to_string()),
                ("middle_seq_id".to_string(), "middle".to_string()),
                ("right_seq_id".to_string(), "right".to_string()),
                ("overlap_bp".to_string(), "8".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("nebuilder validate-only ok");
    assert_eq!(nebuilder_ok.output["can_execute"].as_bool(), Some(true));

    let nebuilder_bad = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "nebuilder_multi_fragment_overlap".to_string(),
            bindings: HashMap::from([
                ("left_seq_id".to_string(), "left".to_string()),
                ("middle_seq_id".to_string(), "middle".to_string()),
                ("right_seq_id".to_string(), "right".to_string()),
                ("overlap_bp".to_string(), "40".to_string()),
            ]),
            transactional: false,
            validate_only: true,
        },
    )
    .expect("nebuilder validate-only bad");
    assert_eq!(nebuilder_bad.output["can_execute"].as_bool(), Some(false));
}

#[test]
fn execute_macros_template_run_new_family_packs_transactionally() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "gg_vector".to_string(),
        DNAsequence::from_sequence("AAAAGGTCTCTTTTGGTCTCAAAA").expect("sequence"),
    );
    state.sequences.insert(
        "gg_insert".to_string(),
        DNAsequence::from_sequence("CCCCGGTCTCGGGG").expect("sequence"),
    );
    state.sequences.insert(
        "gw_donor".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGA").expect("sequence"),
    );
    state.sequences.insert(
        "gw_entry".to_string(),
        DNAsequence::from_sequence("TTTTACGTACGT").expect("sequence"),
    );
    state.sequences.insert(
        "topo_vector_t".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    state.sequences.insert(
        "topo_insert_a".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGA").expect("sequence"),
    );
    state.sequences.insert(
        "inf_left".to_string(),
        DNAsequence::from_sequence("ACACACGGGGAAAA").expect("sequence"),
    );
    state.sequences.insert(
        "inf_right".to_string(),
        DNAsequence::from_sequence("GGGGAAAATTTTCCCC").expect("sequence"),
    );
    state.sequences.insert(
        "neb_left".to_string(),
        DNAsequence::from_sequence("TTTTCCCCAAAAGGGG").expect("sequence"),
    );
    state.sequences.insert(
        "neb_right".to_string(),
        DNAsequence::from_sequence("AAAAGGGGCCCTTTAA").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let base = env!("CARGO_MANIFEST_DIR");
    for path in [
        format!(
            "{}/assets/cloning_patterns_catalog/golden_gate/type_iis/golden_gate_single_insert.json",
            base
        ),
        format!(
            "{}/assets/cloning_patterns_catalog/gateway/bp_lr/gateway_bp_single_insert.json",
            base
        ),
        format!(
            "{}/assets/cloning_patterns_catalog/topo/entry/topo_ta_single_insert.json",
            base
        ),
        format!(
            "{}/assets/cloning_patterns_catalog/ta_gc/entry/ta_clone_single_insert.json",
            base
        ),
        format!(
            "{}/assets/cloning_patterns_catalog/infusion/overlap/infusion_two_fragment_overlap.json",
            base
        ),
        format!(
            "{}/assets/cloning_patterns_catalog/nebuilder_hifi/overlap/nebuilder_two_fragment_overlap.json",
            base
        ),
    ] {
        execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
            .expect("import template");
    }

    // Golden Gate sticky-end compatibility depends on fragment pair selection.
    // Probe deterministic digest fragment index combinations and keep the first
    // transactional run that yields at least one ligation product.
    let mut gg_output_id: Option<String> = None;
    let mut gg_last_error: Option<String> = None;
    'gg_search: for vector_fragment in 1..=3 {
        for insert_fragment in 1..=2 {
            let output_id = format!("gg_tx_out_v{}_i{}", vector_fragment, insert_fragment);
            let ligation_prefix = format!("gg_run_v{}_i{}", vector_fragment, insert_fragment);
            let run = execute_shell_command(
                &mut engine,
                &ShellCommand::MacrosTemplateRun {
                    name: "golden_gate_single_insert".to_string(),
                    bindings: HashMap::from([
                        ("vector_seq_id".to_string(), "gg_vector".to_string()),
                        ("insert_seq_id".to_string(), "gg_insert".to_string()),
                        ("type_iis_enzyme".to_string(), "Eco31".to_string()),
                        ("vector_fragment".to_string(), vector_fragment.to_string()),
                        ("insert_fragment".to_string(), insert_fragment.to_string()),
                        ("ligation_prefix".to_string(), ligation_prefix),
                        ("output_id".to_string(), output_id.clone()),
                    ]),
                    transactional: true,
                    validate_only: false,
                },
            );
            match run {
                Ok(out) if out.state_changed => {
                    gg_output_id = Some(output_id);
                    break 'gg_search;
                }
                Ok(_) => {
                    gg_last_error = Some(format!(
                        "golden_gate_single_insert v{} i{} returned state_changed=false",
                        vector_fragment, insert_fragment
                    ));
                }
                Err(err) => {
                    gg_last_error = Some(format!(
                        "golden_gate_single_insert v{} i{} failed: {}",
                        vector_fragment, insert_fragment, err
                    ));
                }
            }
        }
    }
    let gg_output_id = gg_output_id.unwrap_or_else(|| {
        panic!(
            "{}",
            gg_last_error.unwrap_or_else(|| {
                "No Golden Gate fragment pair produced a ligation product".to_string()
            })
        )
    });

    let runs = vec![
        ShellCommand::MacrosTemplateRun {
            name: "gateway_bp_single_insert".to_string(),
            bindings: HashMap::from([
                ("donor_seq_id".to_string(), "gw_donor".to_string()),
                ("entry_vector_seq_id".to_string(), "gw_entry".to_string()),
                ("gateway_phase".to_string(), "bp".to_string()),
                ("att_tokens".to_string(), "attB,attP".to_string()),
                ("assembly_prefix".to_string(), "gw_tx".to_string()),
                ("output_id".to_string(), "gw_tx_out".to_string()),
            ]),
            transactional: true,
            validate_only: false,
        },
        ShellCommand::MacrosTemplateRun {
            name: "topo_ta_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                ("insert_seq_id".to_string(), "topo_insert_a".to_string()),
                ("topo_mode".to_string(), "ta".to_string()),
                ("assembly_prefix".to_string(), "topo_tx".to_string()),
                ("output_id".to_string(), "topo_tx_out".to_string()),
            ]),
            transactional: true,
            validate_only: false,
        },
        ShellCommand::MacrosTemplateRun {
            name: "ta_clone_single_insert".to_string(),
            bindings: HashMap::from([
                ("vector_seq_id".to_string(), "topo_vector_t".to_string()),
                ("insert_seq_id".to_string(), "topo_insert_a".to_string()),
                ("tail_mode".to_string(), "ta".to_string()),
                ("assembly_prefix".to_string(), "ta_tx".to_string()),
                ("output_id".to_string(), "ta_tx_out".to_string()),
            ]),
            transactional: true,
            validate_only: false,
        },
        ShellCommand::MacrosTemplateRun {
            name: "infusion_two_fragment_overlap".to_string(),
            bindings: HashMap::from([
                ("left_seq_id".to_string(), "inf_left".to_string()),
                ("right_seq_id".to_string(), "inf_right".to_string()),
                ("overlap_bp".to_string(), "8".to_string()),
                ("assembly_prefix".to_string(), "inf_tx".to_string()),
                ("output_id".to_string(), "inf_tx_out".to_string()),
            ]),
            transactional: true,
            validate_only: false,
        },
        ShellCommand::MacrosTemplateRun {
            name: "nebuilder_two_fragment_overlap".to_string(),
            bindings: HashMap::from([
                ("left_seq_id".to_string(), "neb_left".to_string()),
                ("right_seq_id".to_string(), "neb_right".to_string()),
                ("overlap_bp".to_string(), "8".to_string()),
                ("assembly_prefix".to_string(), "neb_tx".to_string()),
                ("output_id".to_string(), "neb_tx_out".to_string()),
            ]),
            transactional: true,
            validate_only: false,
        },
    ];
    for run_cmd in runs {
        let out = execute_shell_command(&mut engine, &run_cmd).expect("transactional run");
        assert!(out.state_changed);
    }

    assert!(engine.state().sequences.contains_key(&gg_output_id));
    for seq_id in [
        "gw_tx_out",
        "topo_tx_out",
        "ta_tx_out",
        "inf_tx_out",
        "neb_tx_out",
    ] {
        assert!(engine.state().sequences.contains_key(seq_id));
    }
}

#[test]
fn execute_macros_template_run_gibson_preview_creates_deterministic_outputs() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "left".to_string(),
        DNAsequence::from_sequence("ATGCCGTTACCGGTTAAACCCGGGTTT").expect("sequence"),
    );
    state.sequences.insert(
        "right".to_string(),
        DNAsequence::from_sequence("AACCCGGGTTTGGGAAATTTCCCGGG").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let path = format!(
        "{}/assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json",
        env!("CARGO_MANIFEST_DIR")
    );
    execute_shell_command(&mut engine, &ShellCommand::MacrosTemplateImport { path })
        .expect("import gibson template");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosTemplateRun {
            name: "gibson_two_fragment_overlap_preview".to_string(),
            bindings: HashMap::from([
                ("left_seq_id".to_string(), "left".to_string()),
                ("right_seq_id".to_string(), "right".to_string()),
                ("overlap_bp".to_string(), "11".to_string()),
                ("assembly_prefix".to_string(), "gib_preview".to_string()),
                ("output_id".to_string(), "gib_forward".to_string()),
            ]),
            transactional: true,
            validate_only: false,
        },
    )
    .expect("run gibson template");
    assert!(run.state_changed);
    assert!(engine.state().sequences.contains_key("gib_preview_1"));
    assert!(engine.state().sequences.contains_key("gib_preview_2"));
    assert!(engine.state().sequences.contains_key("gib_forward"));
}

#[test]
fn execute_candidates_macro_runs_multiple_statements() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let out = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesMacro {
                script: "generate win seqA --length 4 --step 2; score win x '(gc_fraction+1)'; filter win win_top --metric x --min-quantile 1.0".to_string(),
                transactional: false,
            },
        )
        .expect("execute candidates macro");
    assert!(out.state_changed);
    assert_eq!(out.output["executed"].as_u64(), Some(3));
    assert_eq!(out.output["transactional"].as_bool(), Some(false));
    let sets = engine.list_candidate_sets();
    assert!(sets.iter().any(|s| s.name == "win"));
    assert!(sets.iter().any(|s| s.name == "win_top"));
}

#[test]
fn execute_candidates_macro_transactional_rolls_back_on_error() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let err = execute_shell_command(
            &mut engine,
            &ShellCommand::CandidatesMacro {
                script: "generate win seqA --length 4 --step 2; filter missing out --metric gc_fraction --min 0.1".to_string(),
                transactional: true,
            },
        )
        .unwrap_err();
    assert!(err.contains("rolled back"));
    let sets = engine.list_candidate_sets();
    assert!(
        !sets.iter().any(|s| s.name == "win"),
        "transactional macro should roll back earlier statements"
    );
}

#[test]
fn execute_candidates_macro_rejects_nested_macro() {
    let mut engine = GentleEngine::new();
    let err = execute_shell_command(
        &mut engine,
        &ShellCommand::CandidatesMacro {
            script: "macro x".to_string(),
            transactional: false,
        },
    )
    .unwrap_err();
    assert!(err.contains("Nested candidates macro"));
}

#[test]
fn execute_guides_commands_end_to_end() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let tmp = tempdir().expect("tempdir");
    let csv_path = tmp.path().join("guides.csv");
    let protocol_path = tmp.path().join("guides.protocol.txt");

    let guides_json = serde_json::to_string(&vec![
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
            start_0based: 220,
            end_0based_exclusive: 240,
            strand: "+".to_string(),
            protospacer: "TTTTGCCATGTTGACCTGAA".to_string(),
            pam: "TGG".to_string(),
            nuclease: "SpCas9".to_string(),
            cut_offset_from_protospacer_start: 17,
            rank: Some(2),
        },
    ])
    .expect("serialize guides");

    let put = execute_shell_command(
        &mut engine,
        &ShellCommand::GuidesPut {
            guide_set_id: "tp73_guides".to_string(),
            guides_json,
        },
    )
    .expect("guides put");
    assert!(put.state_changed);

    let filter = execute_shell_command(
        &mut engine,
        &ShellCommand::GuidesFilter {
            guide_set_id: "tp73_guides".to_string(),
            config_json: Some(
                "{\"gc_min\":0.3,\"gc_max\":0.7,\"avoid_u6_terminator_tttt\":true}".to_string(),
            ),
            output_guide_set_id: Some("tp73_pass".to_string()),
        },
    )
    .expect("guides filter");
    assert!(filter.state_changed);

    let report = execute_shell_command(
        &mut engine,
        &ShellCommand::GuidesFilterShow {
            guide_set_id: "tp73_guides".to_string(),
        },
    )
    .expect("guides filter-show");
    assert!(!report.state_changed);
    let report_rows = report.output["report"]["results"]
        .as_array()
        .map(|rows| rows.len())
        .unwrap_or_default();
    assert_eq!(report_rows, 2, "expected filter report rows for all guides");

    let generated = execute_shell_command(
        &mut engine,
        &ShellCommand::GuidesOligosGenerate {
            guide_set_id: "tp73_guides".to_string(),
            template_id: "lenti_bsmbi_u6_default".to_string(),
            apply_5prime_g_extension: true,
            output_oligo_set_id: Some("tp73_oligos".to_string()),
            passed_only: true,
        },
    )
    .expect("guides oligos-generate");
    assert!(generated.state_changed);

    let listed = execute_shell_command(
        &mut engine,
        &ShellCommand::GuidesOligosList {
            guide_set_id: Some("tp73_guides".to_string()),
        },
    )
    .expect("guides oligos-list");
    assert!(!listed.state_changed);
    assert_eq!(listed.output["oligo_set_count"].as_u64(), Some(1));

    let exported = execute_shell_command(
        &mut engine,
        &ShellCommand::GuidesOligosExport {
            guide_set_id: "tp73_guides".to_string(),
            oligo_set_id: Some("tp73_oligos".to_string()),
            format: GuideOligoExportFormat::CsvTable,
            path: csv_path.to_string_lossy().to_string(),
            plate_format: None,
        },
    )
    .expect("guides oligos-export");
    assert!(exported.state_changed);
    let csv = fs::read_to_string(&csv_path).expect("read csv export");
    assert!(csv.contains("guide_id,rank,forward_oligo,reverse_oligo,notes"));

    let protocol = execute_shell_command(
        &mut engine,
        &ShellCommand::GuidesProtocolExport {
            guide_set_id: "tp73_guides".to_string(),
            oligo_set_id: Some("tp73_oligos".to_string()),
            path: protocol_path.to_string_lossy().to_string(),
            include_qc_checklist: true,
        },
    )
    .expect("guides protocol-export");
    assert!(protocol.state_changed);
    let protocol_text = fs::read_to_string(&protocol_path).expect("read protocol export");
    assert!(protocol_text.contains("GENtle Guide Oligo Protocol"));
}

#[test]
fn execute_primers_design_list_show_export() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        )
        .expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempdir().expect("tempdir");
    let export_path = tmp.path().join("primer_report.json");
    let request = serde_json::to_string(&Operation::DesignPrimerPairs {
        template: "tpl".to_string(),
        roi_start_0based: 30,
        roi_end_0based: 70,
        forward: crate::engine::PrimerDesignSideConstraint {
            min_length: 20,
            max_length: 20,
            location_0based: Some(5),
            start_0based: None,
            end_0based: None,
            min_tm_c: 40.0,
            max_tm_c: 90.0,
            min_gc_fraction: 0.0,
            max_gc_fraction: 1.0,
            max_anneal_hits: 10,
            ..Default::default()
        },
        reverse: crate::engine::PrimerDesignSideConstraint {
            min_length: 20,
            max_length: 20,
            location_0based: Some(90),
            start_0based: None,
            end_0based: None,
            min_tm_c: 40.0,
            max_tm_c: 90.0,
            min_gc_fraction: 0.0,
            max_gc_fraction: 1.0,
            max_anneal_hits: 10,
            ..Default::default()
        },
        min_amplicon_bp: 40,
        max_amplicon_bp: 150,
        pair_constraints: crate::engine::PrimerDesignPairConstraint::default(),
        max_tm_delta_c: Some(100.0),
        max_pairs: Some(10),
        report_id: Some("tp73_roi".to_string()),
    })
    .expect("serialize request");

    let design = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersDesign {
            request_json: request,
            backend: Some(PrimerDesignBackend::Internal),
            primer3_executable: None,
        },
    )
    .expect("primers design");
    assert!(design.state_changed);
    let report_id = design.output["report"]["report_id"]
        .as_str()
        .unwrap_or_default()
        .to_string();
    assert_eq!(report_id, "tp73_roi");
    assert!(design.output["simple_pcr_pairs"].is_array());
    assert_eq!(
        design.output["simple_pcr_pairs"][0]["flanks_core_cleanly"].as_bool(),
        Some(true)
    );

    let listed = execute_shell_command(&mut engine, &ShellCommand::PrimersListReports)
        .expect("primers list-reports");
    assert!(!listed.state_changed);
    assert_eq!(
        listed.output["schema"].as_str(),
        Some("gentle.primer_design_report_list.v1")
    );
    assert_eq!(listed.output["report_count"].as_u64(), Some(1));

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersShowReport {
            report_id: report_id.clone(),
        },
    )
    .expect("primers show-report");
    assert!(!shown.state_changed);
    assert_eq!(
        shown.output["report"]["report_id"].as_str(),
        Some("tp73_roi")
    );
    assert!(shown.output["simple_pcr_pairs"].is_array());
    assert_eq!(
        shown.output["simple_pcr_pairs"][0]["left_to_core_label"]
            .as_str()
            .map(|value| value.ends_with("bp")),
        Some(true)
    );
    assert_eq!(
        shown.output["simple_pcr_pairs"][0]["right_to_core_label"]
            .as_str()
            .map(|value| value.ends_with("bp")),
        Some(true)
    );

    let exported = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersExportReport {
            report_id,
            path: export_path.to_string_lossy().to_string(),
        },
    )
    .expect("primers export-report");
    assert!(!exported.state_changed);
    assert_eq!(
        exported.output["schema"].as_str(),
        Some("gentle.primer_design_report_export.v1")
    );
    let text = fs::read_to_string(&export_path).expect("read export");
    assert!(text.contains("gentle.primer_design_report.v1"));
}

#[test]
fn execute_primers_prepare_restriction_cloning_returns_saved_report() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        )
        .expect("template"),
    );
    let mut vector = DNAsequence::from_sequence("AAAAGAATTCGGGGGAAGCTTTTTT").expect("vector");
    *vector.restriction_enzymes_mut() = crate::enzymes::active_restriction_enzymes();
    let vector_len_i64 = vector.len().try_into().unwrap();
    vector.features_mut().push(Feature {
        kind: "misc_feature".into(),
        location: Location::simple_range(0, vector_len_i64),
        qualifiers: vec![
            ("label".into(), Some("MCS".to_string())),
            (
                "mcs_expected_sites".into(),
                Some("EcoRI,HindIII".to_string()),
            ),
        ],
    });
    vector.update_computed_features();
    state.sequences.insert("vec".to_string(), vector);
    let mut engine = GentleEngine::from_state(state);

    let primer_request = serde_json::to_string(&Operation::DesignPrimerPairs {
        template: "tpl".to_string(),
        roi_start_0based: 40,
        roi_end_0based: 80,
        forward: crate::engine::PrimerDesignSideConstraint {
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
            ..Default::default()
        },
        reverse: crate::engine::PrimerDesignSideConstraint {
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
            ..Default::default()
        },
        min_amplicon_bp: 40,
        max_amplicon_bp: 150,
        pair_constraints: crate::engine::PrimerDesignPairConstraint::default(),
        max_tm_delta_c: Some(100.0),
        max_pairs: Some(10),
        report_id: Some("shell_handoff_pairs".to_string()),
    })
    .expect("serialize primer design request");
    execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersDesign {
            request_json: primer_request,
            backend: Some(PrimerDesignBackend::Internal),
            primer3_executable: None,
        },
    )
    .expect("design primer pairs");

    let handoff_request = serde_json::to_string(&Operation::PrepareRestrictionCloningPcrHandoff {
        template: "tpl".to_string(),
        primer_report_id: "shell_handoff_pairs".to_string(),
        pair_index: 0,
        destination_vector_seq_id: "vec".to_string(),
        mode: crate::engine::RestrictionCloningPcrHandoffMode::DirectedPair,
        forward_enzyme: "EcoRI".to_string(),
        reverse_enzyme: Some("HindIII".to_string()),
        forward_leader_5prime: Some("GC".to_string()),
        reverse_leader_5prime: Some("AT".to_string()),
    })
    .expect("serialize handoff request");
    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersPrepareRestrictionCloning {
            request_json: handoff_request,
        },
    )
    .expect("execute handoff shell command");

    assert!(run.state_changed);
    assert_eq!(
        run.output["report"]["schema"].as_str(),
        Some("gentle.restriction_cloning_pcr_handoff.v1")
    );
    assert_eq!(
        run.output["report"]["destination_vector_seq_id"].as_str(),
        Some("vec")
    );
    assert_eq!(run.output["report"]["mode"].as_str(), Some("directed_pair"));
    assert_eq!(
        run.output["report"]["compatibility"]["status"].as_str(),
        Some("compatible")
    );
    assert_eq!(
        run.output["report"]["workflow_hints"]["pcr_advanced_operation"]["PcrAdvanced"]["forward_primer"]["anneal_len"].as_u64(),
        Some(20)
    );
}

#[test]
fn execute_primers_restriction_cloning_vector_suggestions_prefers_mcs_then_unique() {
    let mut state = ProjectState::default();
    let mut vector =
        DNAsequence::from_sequence("AAAAGAATTCGGGGGAAGCTTTTTTGCGGCCGCTTTT").expect("vector");
    *vector.restriction_enzymes_mut() = crate::enzymes::active_restriction_enzymes();
    let vector_len_i64 = vector.len().try_into().unwrap();
    vector.features_mut().push(Feature {
        kind: "misc_feature".into(),
        location: Location::simple_range(0, vector_len_i64),
        qualifiers: vec![
            ("label".into(), Some("MCS".to_string())),
            (
                "mcs_expected_sites".into(),
                Some("EcoRI,HindIII,BamHI".to_string()),
            ),
        ],
    });
    vector.update_computed_features();
    state.sequences.insert("vec".to_string(), vector);
    let mut engine = GentleEngine::from_state(state);

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersRestrictionCloningVectorSuggestions {
            seq_id: "vec".to_string(),
        },
    )
    .expect("execute restriction-cloning vector suggestions");

    assert!(!run.state_changed);
    assert_eq!(
        run.output["schema"].as_str(),
        Some("gentle.restriction_cloning_vector_enzyme_suggestions.v1")
    );
    assert_eq!(run.output["suggestions"]["seq_id"].as_str(), Some("vec"));
    assert_eq!(
        run.output["suggestions"]["selected_mcs"][0].as_str(),
        Some("EcoRI")
    );
    assert_eq!(
        run.output["suggestions"]["selected_mcs"][1].as_str(),
        Some("HindIII")
    );
    let other_unique = run.output["suggestions"]["other_unique"]
        .as_array()
        .expect("other_unique array");
    assert!(
        other_unique
            .iter()
            .filter_map(|value| value.as_str())
            .any(|name| name == "NotI"),
        "expected NotI in other_unique suggestions, got {other_unique:?}"
    );
    assert_eq!(
        run.output["suggestions"]["missing_mcs"][0].as_str(),
        Some("BamHI")
    );
    assert_eq!(
        run.output["suggestions"]["recommended_single_site"][0]["enzyme"].as_str(),
        Some("EcoRI")
    );
    assert_eq!(
        run.output["suggestions"]["recommended_directed_pairs"][0]["forward_enzyme"].as_str(),
        Some("EcoRI")
    );
    assert_eq!(
        run.output["suggestions"]["recommended_directed_pairs"][0]["reverse_enzyme"].as_str(),
        Some("HindIII")
    );
}

#[test]
fn execute_primers_seed_restriction_cloning_handoff_uses_recommendations_or_validates_order() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        )
        .expect("template"),
    );
    let mut vector =
        DNAsequence::from_sequence("AAAAGAATTCGGGGGAAGCTTTTTTGCGGCCGCTTTT").expect("vector");
    *vector.restriction_enzymes_mut() = crate::enzymes::active_restriction_enzymes();
    let vector_len_i64 = vector.len().try_into().unwrap();
    vector.features_mut().push(Feature {
        kind: "misc_feature".into(),
        location: Location::simple_range(0, vector_len_i64),
        qualifiers: vec![
            ("label".into(), Some("MCS".to_string())),
            (
                "mcs_expected_sites".into(),
                Some("EcoRI,HindIII".to_string()),
            ),
        ],
    });
    vector.update_computed_features();
    state.sequences.insert("vec".to_string(), vector);
    let mut engine = GentleEngine::from_state(state);
    engine.state_mut().parameters.primer_design_backend = PrimerDesignBackend::Internal;

    engine
        .apply(Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: crate::engine::PrimerDesignSideConstraint {
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
                ..Default::default()
            },
            reverse: crate::engine::PrimerDesignSideConstraint {
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
                ..Default::default()
            },
            pair_constraints: Default::default(),
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("primer_seed_handoff".to_string()),
        })
        .expect("design primer pairs");

    let single = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersSeedRestrictionCloningHandoff {
            primer_report_id: "primer_seed_handoff".to_string(),
            destination_vector_seq_id: "vec".to_string(),
            pair_rank_1based: 1,
            mode: RestrictionCloningPcrHandoffMode::SingleSite,
            forward_enzyme: None,
            reverse_enzyme: None,
            forward_leader_5prime: None,
            reverse_leader_5prime: None,
        },
    )
    .expect("seed single-site handoff");
    assert!(!single.state_changed);
    assert_eq!(
        single.output["schema"].as_str(),
        Some("gentle.restriction_cloning_pcr_handoff_seed.v1")
    );
    assert_eq!(
        single.output["selection_source"].as_str(),
        Some("recommended_single_site")
    );
    assert_eq!(single.output["forward_enzyme"].as_str(), Some("EcoRI"));
    assert_eq!(single.output["reverse_enzyme"].as_str(), Some("EcoRI"));
    assert_eq!(
        single.output["operation"]["PrepareRestrictionCloningPcrHandoff"]["forward_enzyme"]
            .as_str(),
        Some("EcoRI")
    );

    let directed = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersSeedRestrictionCloningHandoff {
            primer_report_id: "primer_seed_handoff".to_string(),
            destination_vector_seq_id: "vec".to_string(),
            pair_rank_1based: 1,
            mode: RestrictionCloningPcrHandoffMode::DirectedPair,
            forward_enzyme: None,
            reverse_enzyme: None,
            forward_leader_5prime: Some("GC".to_string()),
            reverse_leader_5prime: Some("AT".to_string()),
        },
    )
    .expect("seed directed handoff");
    assert_eq!(
        directed.output["selection_source"].as_str(),
        Some("recommended_directed_pair")
    );
    assert_eq!(
        directed.output["suggestion_order_source"].as_str(),
        Some("mcs_expected_sites")
    );
    assert_eq!(directed.output["forward_enzyme"].as_str(), Some("EcoRI"));
    assert_eq!(directed.output["reverse_enzyme"].as_str(), Some("HindIII"));
    assert_eq!(
        directed.output["operation"]["PrepareRestrictionCloningPcrHandoff"]["reverse_enzyme"]
            .as_str(),
        Some("HindIII")
    );
    assert_eq!(
        directed.output["operation"]["PrepareRestrictionCloningPcrHandoff"]
            ["forward_leader_5prime"]
            .as_str(),
        Some("GC")
    );

    let reversed_err = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersSeedRestrictionCloningHandoff {
            primer_report_id: "primer_seed_handoff".to_string(),
            destination_vector_seq_id: "vec".to_string(),
            pair_rank_1based: 1,
            mode: RestrictionCloningPcrHandoffMode::DirectedPair,
            forward_enzyme: Some("HindIII".to_string()),
            reverse_enzyme: Some("EcoRI".to_string()),
            forward_leader_5prime: None,
            reverse_leader_5prime: None,
        },
    )
    .expect_err("reversed directed-pair order should fail");
    assert!(reversed_err.contains("does not match the valid order"));
    assert!(reversed_err.contains("HindIII -> EcoRI"));

    let unavailable_err = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersSeedRestrictionCloningHandoff {
            primer_report_id: "primer_seed_handoff".to_string(),
            destination_vector_seq_id: "vec".to_string(),
            pair_rank_1based: 1,
            mode: RestrictionCloningPcrHandoffMode::SingleSite,
            forward_enzyme: Some("BamHI".to_string()),
            reverse_enzyme: None,
            forward_leader_5prime: None,
            reverse_leader_5prime: None,
        },
    )
    .expect_err("unavailable cutter should fail");
    assert!(unavailable_err.contains("is not a unique usable cutter on vector 'vec'"));
}

#[test]
fn execute_primers_restriction_cloning_handoff_report_commands() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(
            "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
        )
        .expect("template"),
    );
    let mut vector = DNAsequence::from_sequence("AAAAGAATTCGGGGGAAGCTTTTTT").expect("vector");
    *vector.restriction_enzymes_mut() = crate::enzymes::active_restriction_enzymes();
    let vector_len_i64 = vector.len().try_into().unwrap();
    vector.features_mut().push(Feature {
        kind: "misc_feature".into(),
        location: Location::simple_range(0, vector_len_i64),
        qualifiers: vec![
            ("label".into(), Some("MCS".to_string())),
            (
                "mcs_expected_sites".into(),
                Some("EcoRI,HindIII".to_string()),
            ),
        ],
    });
    vector.update_computed_features();
    state.sequences.insert("vec".to_string(), vector);
    let mut engine = GentleEngine::from_state(state);

    let primer_request = serde_json::to_string(&Operation::DesignPrimerPairs {
        template: "tpl".to_string(),
        roi_start_0based: 40,
        roi_end_0based: 80,
        forward: crate::engine::PrimerDesignSideConstraint {
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
            ..Default::default()
        },
        reverse: crate::engine::PrimerDesignSideConstraint {
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
            ..Default::default()
        },
        min_amplicon_bp: 40,
        max_amplicon_bp: 150,
        pair_constraints: crate::engine::PrimerDesignPairConstraint::default(),
        max_tm_delta_c: Some(100.0),
        max_pairs: Some(10),
        report_id: Some("shell_handoff_list".to_string()),
    })
    .expect("serialize primer design request");
    execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersDesign {
            request_json: primer_request,
            backend: Some(PrimerDesignBackend::Internal),
            primer3_executable: None,
        },
    )
    .expect("design primer pairs");

    let handoff_request = serde_json::to_string(&Operation::PrepareRestrictionCloningPcrHandoff {
        template: "tpl".to_string(),
        primer_report_id: "shell_handoff_list".to_string(),
        pair_index: 0,
        destination_vector_seq_id: "vec".to_string(),
        mode: crate::engine::RestrictionCloningPcrHandoffMode::DirectedPair,
        forward_enzyme: "EcoRI".to_string(),
        reverse_enzyme: Some("HindIII".to_string()),
        forward_leader_5prime: Some("GC".to_string()),
        reverse_leader_5prime: Some("AT".to_string()),
    })
    .expect("serialize handoff request");
    let create = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersPrepareRestrictionCloning {
            request_json: handoff_request,
        },
    )
    .expect("create restriction-cloning handoff");
    let report_id = create.output["report"]["report_id"]
        .as_str()
        .expect("created handoff report id")
        .to_string();

    let list = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersListRestrictionCloningHandoffs,
    )
    .expect("list restriction-cloning handoffs");
    assert_eq!(
        list.output["schema"].as_str(),
        Some("gentle.restriction_cloning_pcr_handoff_report_list.v1")
    );
    assert!(list.output["report_count"].as_u64().unwrap_or(0) >= 1);

    let show = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersShowRestrictionCloningHandoff {
            report_id: report_id.clone(),
        },
    )
    .expect("show restriction-cloning handoff");
    assert_eq!(
        show.output["report"]["report_id"].as_str(),
        Some(report_id.as_str())
    );

    let td = tempdir().expect("tempdir");
    let export_path = td
        .path()
        .join("restriction_cloning_handoff_export.json")
        .to_string_lossy()
        .to_string();
    let export = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersExportRestrictionCloningHandoff {
            report_id: report_id.clone(),
            path: export_path.clone(),
        },
    )
    .expect("export restriction-cloning handoff");
    assert_eq!(
        export.output["schema"].as_str(),
        Some("gentle.restriction_cloning_pcr_handoff_report_export.v1")
    );
    let text = fs::read_to_string(&export_path).expect("read handoff export");
    assert!(text.contains("gentle.restriction_cloning_pcr_handoff.v1"));
}

#[test]
fn execute_features_query_returns_structured_rows() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(120)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "gene".into(),
        location: Location::simple_range(10, 90),
        qualifiers: vec![
            ("label".into(), Some("TP73".to_string())),
            ("gene".into(), Some("TP73".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "CDS".into(),
        location: Location::Complement(Box::new(Location::simple_range(20, 80))),
        qualifiers: vec![
            ("label".into(), Some("TP73_DBD".to_string())),
            ("product".into(), Some("TP73 DNA-binding".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "misc_feature".into(),
        location: Location::simple_range(120, 170),
        qualifiers: vec![("note".into(), Some("promoter".to_string()))],
    });
    state.sequences.insert("tp73".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let query = SequenceFeatureQuery {
        seq_id: "tp73".to_string(),
        kind_in: vec!["CDS".to_string()],
        strand: SequenceFeatureStrandFilter::Reverse,
        start_0based: Some(0),
        end_0based_exclusive: Some(120),
        range_relation: SequenceFeatureRangeRelation::Within,
        include_qualifiers: true,
        limit: Some(10),
        ..SequenceFeatureQuery::default()
    };
    let run = execute_shell_command(&mut engine, &ShellCommand::FeaturesQuery { query })
        .expect("features query");
    assert!(!run.state_changed);
    assert_eq!(
        run.output["schema"].as_str(),
        Some("gentle.sequence_feature_query_result.v1")
    );
    assert_eq!(run.output["matched_count"].as_u64(), Some(1));
    assert_eq!(run.output["returned_count"].as_u64(), Some(1));
    assert_eq!(run.output["rows"][0]["kind"].as_str(), Some("CDS"));
    assert_eq!(run.output["rows"][0]["strand"].as_str(), Some("reverse"));
}

#[test]
fn execute_features_export_bed_writes_tfbs_and_restriction_rows() {
    let td = tempdir().expect("tempdir");
    let output = td.path().join("features.bed");

    let mut dna = DNAsequence::from_sequence("AAGAATTCTTTACGTAA").expect("sequence");
    *dna.restriction_enzymes_mut() = crate::enzymes::active_restriction_enzymes();
    dna.update_computed_features();
    dna.features_mut().push(Feature {
        kind: "TFBS".into(),
        location: Location::simple_range(11, 15),
        qualifiers: vec![
            ("bound_moiety".into(), Some("SP1".to_string())),
            ("tf_id".into(), Some("MA0079.3".to_string())),
            ("llr_quantile".into(), Some("0.95".to_string())),
        ],
    });

    let mut state = ProjectState::default();
    state.sequences.insert("seq".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesExportBed {
            query: SequenceFeatureQuery {
                seq_id: "seq".to_string(),
                sort_by: SequenceFeatureSortBy::Start,
                ..SequenceFeatureQuery::default()
            },
            output: output.display().to_string(),
            coordinate_mode: FeatureBedCoordinateMode::Auto,
            include_restriction_sites: true,
            restriction_enzymes: vec!["EcoRI".to_string()],
        },
    )
    .expect("features export-bed");

    assert!(!run.state_changed);
    assert_eq!(
        run.output["schema"].as_str(),
        Some("gentle.sequence_feature_bed_export.v1")
    );
    assert_eq!(
        run.output["matched_sequence_feature_count"].as_u64(),
        Some(1)
    );
    assert_eq!(
        run.output["matched_restriction_site_count"].as_u64(),
        Some(1)
    );
    assert_eq!(run.output["exported_row_count"].as_u64(), Some(2));

    let bed = fs::read_to_string(&output).expect("read BED");
    assert!(bed.contains("\tTFBS\tfeature:0\t"));
    assert!(bed.contains("\trestriction_site\trestriction_site:ecori:2:8:0\t"));
    assert!(bed.contains("\t950\t+\tTFBS\t"));
}

#[test]
fn execute_features_tfbs_summary_returns_grouped_counts_and_densities() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(250)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "TFBS".into(),
        location: Location::simple_range(95, 105),
        qualifiers: vec![
            ("bound_moiety".into(), Some("SP1".to_string())),
            ("tf_id".into(), Some("MA0079.3".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "TFBS".into(),
        location: Location::simple_range(150, 160),
        qualifiers: vec![
            ("bound_moiety".into(), Some("SP1".to_string())),
            ("tf_id".into(), Some("MA0079.4".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "TFBS".into(),
        location: Location::simple_range(400, 410),
        qualifiers: vec![
            ("bound_moiety".into(), Some("SP1".to_string())),
            ("tf_id".into(), Some("MA0079.3".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "TFBS".into(),
        location: Location::simple_range(210, 220),
        qualifiers: vec![
            ("bound_moiety".into(), Some("CTCF".to_string())),
            ("tf_id".into(), Some("MA0139.1".to_string())),
        ],
    });
    state.sequences.insert("promoter".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesTfbsSummary {
            request: TfbsRegionSummaryRequest {
                seq_id: "promoter".to_string(),
                focus_start_0based: 90,
                focus_end_0based_exclusive: 200,
                context_start_0based: Some(0),
                context_end_0based_exclusive: Some(500),
                ..TfbsRegionSummaryRequest::default()
            },
        },
    )
    .expect("features tfbs-summary");
    assert!(!run.state_changed);
    assert_eq!(
        run.output["schema"].as_str(),
        Some("gentle.tfbs_region_summary.v1")
    );
    assert_eq!(run.output["matched_tf_count"].as_u64(), Some(1));
    assert_eq!(run.output["rows"][0]["tf_name"].as_str(), Some("SP1"));
    assert_eq!(run.output["rows"][0]["focus_occurrences"].as_u64(), Some(2));
    assert_eq!(
        run.output["rows"][0]["outside_focus_occurrences"].as_u64(),
        Some(1)
    );
}

#[test]
fn execute_features_tfbs_score_tracks_svg_shell_command_writes_svg_and_report() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tp73_context".to_string(),
        DNAsequence::from_sequence(&"ACGT".repeat(80)).expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let td = tempdir().expect("tempdir");
    let output = td.path().join("tp73_tracks.svg");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesTfbsScoreTracksSvg {
            target: SequenceScanTarget::SeqId {
                seq_id: "tp73_context".to_string(),
                span_start_0based: Some(0),
                span_end_0based_exclusive: Some(120),
            },
            motifs: vec!["TP73".to_string(), "SP1".to_string(), "BACH2".to_string()],
            score_kind: TfbsScoreTrackValueKind::LlrBits,
            clip_negative: true,
            output: output.display().to_string(),
        },
    )
    .expect("features tfbs-score-tracks-svg");

    assert!(!run.state_changed);
    assert_eq!(
        run.output["result"]["tfbs_score_tracks"]["schema"].as_str(),
        Some("gentle.tfbs_score_tracks.v1")
    );
    assert_eq!(
        run.output["result"]["tfbs_score_tracks"]["tracks"]
            .as_array()
            .map(Vec::len),
        Some(3)
    );
    let svg = fs::read_to_string(&output).expect("read svg");
    assert!(svg.contains("<svg"));
    assert!(svg.contains("Continuous TF motif score tracks"));
    assert!(svg.contains("tp73_context"));
    assert!(svg.contains("SP1"));
}

#[test]
fn execute_features_tfbs_score_tracks_svg_matches_inline_and_stored_targets() {
    let sequence_text = "ACGT".repeat(60);
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tp73_context".to_string(),
        DNAsequence::from_sequence(&sequence_text).expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let td = tempdir().expect("tempdir");
    let stored_output = td.path().join("stored_tracks.svg");
    let inline_output = td.path().join("inline_tracks.svg");

    let stored = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesTfbsScoreTracksSvg {
            target: SequenceScanTarget::SeqId {
                seq_id: "tp73_context".to_string(),
                span_start_0based: Some(10),
                span_end_0based_exclusive: Some(140),
            },
            motifs: vec!["SP1".to_string(), "TP73".to_string()],
            score_kind: TfbsScoreTrackValueKind::LlrQuantile,
            clip_negative: true,
            output: stored_output.display().to_string(),
        },
    )
    .expect("stored tfbs-score-tracks-svg");
    let inline = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesTfbsScoreTracksSvg {
            target: SequenceScanTarget::InlineSequence {
                sequence_text,
                topology: InlineSequenceTopology::Linear,
                id_hint: Some("tp73_context".to_string()),
                span_start_0based: Some(10),
                span_end_0based_exclusive: Some(140),
            },
            motifs: vec!["SP1".to_string(), "TP73".to_string()],
            score_kind: TfbsScoreTrackValueKind::LlrQuantile,
            clip_negative: true,
            output: inline_output.display().to_string(),
        },
    )
    .expect("inline tfbs-score-tracks-svg");

    assert_eq!(
        stored.output["result"]["tfbs_score_tracks"]["tracks"],
        inline.output["result"]["tfbs_score_tracks"]["tracks"]
    );
    assert_eq!(
        stored.output["result"]["tfbs_score_tracks"]["target_kind"].as_str(),
        Some("seq_id")
    );
    assert_eq!(
        inline.output["result"]["tfbs_score_tracks"]["target_kind"].as_str(),
        Some("inline_sequence")
    );
    assert!(
        fs::read_to_string(&stored_output)
            .expect("stored svg")
            .contains("tp73_context")
    );
    assert!(
        fs::read_to_string(&inline_output)
            .expect("inline svg")
            .contains("tp73_context")
    );
}

#[test]
fn execute_resources_resolve_tf_query_reports_builtin_groups_and_aliases() {
    let mut engine = GentleEngine::new();
    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::ResourcesResolveTfQuery {
            queries: vec!["Yamanaka factors".to_string(), "OCT4".to_string()],
            output: None,
        },
    )
    .expect("resources resolve-tf-query");

    assert!(!run.state_changed);
    assert_eq!(
        run.output["result"]["tf_query_resolution_report"]["schema"].as_str(),
        Some("gentle.tf_query_resolution.v1")
    );
    assert_eq!(
        run.output["result"]["tf_query_resolution_report"]["returned_query_count"].as_u64(),
        Some(2)
    );
    assert_eq!(
        run.output["result"]["tf_query_resolution_report"]["entries"][0]["resolution_kind"]
            .as_str(),
        Some("builtin_group")
    );
    assert_eq!(
        run.output["result"]["tf_query_resolution_report"]["entries"][1]["resolution_kind"]
            .as_str(),
        Some("exact_motif")
    );
}

#[test]
fn execute_features_tfbs_score_track_correlation_svg_shell_command_writes_svg_and_report() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tert_context".to_string(),
        DNAsequence::from_sequence(&"ACGT".repeat(90)).expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let td = tempdir().expect("tempdir");
    let output = td.path().join("tert_track_correlation.svg");

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesTfbsScoreTrackCorrelationSvg {
            seq_id: "tert_context".to_string(),
            motifs: vec!["SP1".to_string(), "MYC".to_string(), "PATZ1".to_string()],
            start_0based: 0,
            end_0based_exclusive: 160,
            score_kind: TfbsScoreTrackValueKind::LlrBackgroundTailLog10,
            correlation_metric: TfbsScoreTrackCorrelationMetric::Spearman,
            correlation_signal_source: TfbsScoreTrackCorrelationSignalSource::MaxStrands,
            clip_negative: true,
            output: output.display().to_string(),
        },
    )
    .expect("features tfbs-score-track-correlation-svg");

    assert!(!run.state_changed);
    assert_eq!(
        run.output["result"]["tfbs_score_tracks"]["schema"].as_str(),
        Some("gentle.tfbs_score_tracks.v1")
    );
    assert_eq!(
        run.output["result"]["tfbs_score_tracks"]["correlation_summary"]["pair_count"].as_u64(),
        Some(3)
    );
    let svg = fs::read_to_string(&output).expect("read svg");
    assert!(svg.contains("<svg"));
    assert!(svg.contains("TFBS track correlation"));
    assert!(svg.contains("Smoothed Spearman rho"));
    assert!(svg.contains("Raw Spearman rho"));
    assert!(svg.contains("SP1"));
    assert!(svg.contains("MYC"));
}

#[test]
fn execute_features_tfbs_track_similarity_shell_command_returns_report() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tert_context".to_string(),
        DNAsequence::from_sequence("GGGGCGGGGCACGTGGGGGCGGGGCACGTG".repeat(6).as_str())
            .expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesTfbsTrackSimilarity {
            target: SequenceScanTarget::SeqId {
                seq_id: "tert_context".to_string(),
                span_start_0based: Some(0),
                span_end_0based_exclusive: Some(120),
            },
            anchor_motif: "SP1".to_string(),
            candidate_motifs: vec!["MYC".to_string(), "PATZ1".to_string()],
            ranking_metric: TfbsTrackSimilarityRankingMetric::SmoothedSpearman,
            score_kind: TfbsScoreTrackValueKind::LlrBackgroundTailLog10,
            clip_negative: true,
            species_filters: vec![],
            include_remote_metadata: false,
            limit: Some(2),
            path: None,
        },
    )
    .expect("features tfbs-track-similarity");

    assert!(!run.state_changed);
    assert_eq!(
        run.output["result"]["tfbs_track_similarity"]["schema"].as_str(),
        Some("gentle.tfbs_track_similarity.v1")
    );
    assert_eq!(
        run.output["report"]["anchor_requested"].as_str(),
        Some("SP1")
    );
    assert_eq!(
        run.output["report"]["ranking_metric"].as_str(),
        Some("smoothed_spearman")
    );
    assert_eq!(
        run.output["report"]["returned_candidate_count"].as_u64(),
        Some(2)
    );
}

#[test]
fn execute_features_restriction_scan_matches_inline_and_stored_sequence_targets() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq".to_string(),
        DNAsequence::from_sequence("AAGAATTCCCGGGTTGAATTC").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let stored = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesRestrictionScan {
            target: SequenceScanTarget::SeqId {
                seq_id: "seq".to_string(),
                span_start_0based: None,
                span_end_0based_exclusive: None,
            },
            enzymes: vec!["EcoRI".to_string(), "SmaI".to_string()],
            max_sites_per_enzyme: None,
            include_cut_geometry: true,
            path: None,
        },
    )
    .expect("stored restriction scan");
    let inline = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesRestrictionScan {
            target: SequenceScanTarget::InlineSequence {
                sequence_text: "AAGAATTCCCGGGTTGAATTC".to_string(),
                topology: InlineSequenceTopology::Linear,
                id_hint: Some("seq".to_string()),
                span_start_0based: None,
                span_end_0based_exclusive: None,
            },
            enzymes: vec!["EcoRI".to_string(), "SmaI".to_string()],
            max_sites_per_enzyme: None,
            include_cut_geometry: true,
            path: None,
        },
    )
    .expect("inline restriction scan");

    assert!(!stored.state_changed);
    assert!(!inline.state_changed);
    assert_eq!(
        stored.output["result"]["restriction_site_scan"]["schema"].as_str(),
        Some("gentle.restriction_site_scan.v1")
    );
    assert_eq!(
        stored.output["report"]["matched_site_count"].as_u64(),
        Some(3)
    );
    assert_eq!(
        stored.output["report"]["rows"],
        inline.output["report"]["rows"]
    );
    assert_eq!(
        stored.output["report"]["skipped_enzyme_names_due_to_max_sites"],
        inline.output["report"]["skipped_enzyme_names_due_to_max_sites"]
    );
}

#[test]
fn execute_features_tfbs_scan_matches_inline_and_stored_sequence_targets() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq".to_string(),
        DNAsequence::from_sequence("TTTTATAAAGGGTATAAATTT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);

    let stored = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesTfbsScan {
            target: SequenceScanTarget::SeqId {
                seq_id: "seq".to_string(),
                span_start_0based: None,
                span_end_0based_exclusive: None,
            },
            motifs: vec!["TATAAA".to_string()],
            min_llr_bits: None,
            min_llr_quantile: Some(0.95),
            per_tf_thresholds: vec![],
            max_hits: Some(8),
            path: None,
        },
    )
    .expect("stored tfbs scan");
    let inline = execute_shell_command(
        &mut engine,
        &ShellCommand::FeaturesTfbsScan {
            target: SequenceScanTarget::InlineSequence {
                sequence_text: "TTTTATAAAGGGTATAAATTT".to_string(),
                topology: InlineSequenceTopology::Linear,
                id_hint: Some("seq".to_string()),
                span_start_0based: None,
                span_end_0based_exclusive: None,
            },
            motifs: vec!["TATAAA".to_string()],
            min_llr_bits: None,
            min_llr_quantile: Some(0.95),
            per_tf_thresholds: vec![],
            max_hits: Some(8),
            path: None,
        },
    )
    .expect("inline tfbs scan");

    assert!(!stored.state_changed);
    assert!(!inline.state_changed);
    assert_eq!(
        stored.output["result"]["tfbs_hit_scan"]["schema"].as_str(),
        Some("gentle.tfbs_hit_scan.v1")
    );
    assert!(
        stored.output["report"]["matched_hit_count"]
            .as_u64()
            .is_some_and(|value| value > 0)
    );
    assert_eq!(
        stored.output["report"]["rows"],
        inline.output["report"]["rows"]
    );
    assert_eq!(
        stored.output["report"]["motifs_scanned"],
        inline.output["report"]["motifs_scanned"]
    );
}

#[test]
fn execute_variant_promoter_context_shell_command_returns_report_payload() {
    let mut dna = DNAsequence::from_sequence(&"A".repeat(6001)).expect("sequence");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Complement(Box::new(Location::simple_range(2000, 2613))),
        qualifiers: vec![
            ("gene".into(), Some("VKORC1".to_string())),
            ("transcript_id".into(), Some("ENSTVKORC1".to_string())),
            ("label".into(), Some("VKORC1-201".to_string())),
        ],
    });
    dna.features_mut().push(Feature {
        kind: "variation".into(),
        location: Location::simple_range(3000, 3001),
        qualifiers: vec![
            ("label".into(), Some("rs9923231".to_string())),
            ("db_xref".into(), Some("dbSNP:rs9923231".to_string())),
            ("vcf_ref".into(), Some("C".to_string())),
            ("vcf_alt".into(), Some("T".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("vkorc1_demo".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::VariantPromoterContext {
            seq_id: "vkorc1_demo".to_string(),
            variant_label_or_id: Some("rs9923231".to_string()),
            gene_label: Some("VKORC1".to_string()),
            transcript_id: None,
            promoter_upstream_bp: 1000,
            promoter_downstream_bp: 200,
            tfbs_focus_half_window_bp: 100,
            path: None,
        },
    )
    .expect("variant promoter-context");

    assert!(!run.state_changed);
    assert_eq!(
        run.output["report"]["schema"].as_str(),
        Some("gentle.variant_promoter_context.v1")
    );
    assert_eq!(
        run.output["report"]["promoter_overlap"].as_bool(),
        Some(true)
    );
    assert_eq!(
        run.output["report"]["chosen_gene_label"].as_str(),
        Some("VKORC1")
    );
}

#[test]
fn execute_variant_materialize_allele_shell_command_creates_sequence() {
    let mut dna = DNAsequence::from_sequence("ACCGT").expect("sequence");
    dna.features_mut().push(Feature {
        kind: "variation".into(),
        location: Location::simple_range(2, 3),
        qualifiers: vec![
            ("label".into(), Some("rsDemo".to_string())),
            ("db_xref".into(), Some("dbSNP:rsDemo".to_string())),
            ("vcf_ref".into(), Some("C".to_string())),
            ("vcf_alt".into(), Some("A".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("demo".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::VariantMaterializeAllele {
            seq_id: "demo".to_string(),
            variant_label_or_id: Some("rsDemo".to_string()),
            allele: VariantAlleleChoice::Alternate,
            output_id: Some("demo_alt".to_string()),
        },
    )
    .expect("variant materialize-allele");

    assert!(run.state_changed);
    assert_eq!(
        engine
            .state()
            .sequences
            .get("demo_alt")
            .expect("materialized sequence")
            .get_forward_string(),
        "ACAGT"
    );
    assert_eq!(
        run.output["result"]["created_seq_ids"][0].as_str(),
        Some("demo_alt")
    );
}

#[cfg(unix)]
#[test]
fn execute_primers_design_internal_vs_primer3_fixture_normalization_parity() {
    let mut state = ProjectState::default();
    state.sequences.insert(
            "tpl".to_string(),
            DNAsequence::from_sequence(
                "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
            )
            .expect("sequence"),
        );
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempdir().expect("tempdir");
    let fixture_path = Path::new(&primer3_fixture_path("pairs.location_5_60.kv")).to_path_buf();
    let fake_primer3 = install_fake_primer3(tmp.path(), &fixture_path);

    let forward = crate::engine::PrimerDesignSideConstraint {
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
        ..Default::default()
    };
    let reverse = crate::engine::PrimerDesignSideConstraint {
        min_length: 20,
        max_length: 20,
        location_0based: Some(60),
        start_0based: None,
        end_0based: None,
        min_tm_c: 0.0,
        max_tm_c: 100.0,
        min_gc_fraction: 0.0,
        max_gc_fraction: 1.0,
        max_anneal_hits: 1000,
        ..Default::default()
    };

    let internal_request = serde_json::to_string(&Operation::DesignPrimerPairs {
        template: "tpl".to_string(),
        roi_start_0based: 40,
        roi_end_0based: 80,
        forward: forward.clone(),
        reverse: reverse.clone(),
        min_amplicon_bp: 40,
        max_amplicon_bp: 150,
        pair_constraints: crate::engine::PrimerDesignPairConstraint::default(),
        max_tm_delta_c: Some(100.0),
        max_pairs: Some(10),
        report_id: Some("internal_norm".to_string()),
    })
    .expect("serialize internal request");
    execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersDesign {
            request_json: internal_request,
            backend: Some(PrimerDesignBackend::Internal),
            primer3_executable: None,
        },
    )
    .expect("primers design internal");
    let internal = engine
        .get_primer_design_report("internal_norm")
        .expect("internal report");

    let primer3_request = serde_json::to_string(&Operation::DesignPrimerPairs {
        template: "tpl".to_string(),
        roi_start_0based: 40,
        roi_end_0based: 80,
        forward,
        reverse,
        min_amplicon_bp: 40,
        max_amplicon_bp: 150,
        pair_constraints: crate::engine::PrimerDesignPairConstraint::default(),
        max_tm_delta_c: Some(100.0),
        max_pairs: Some(10),
        report_id: Some("primer3_norm".to_string()),
    })
    .expect("serialize primer3 request");
    execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersDesign {
            request_json: primer3_request,
            backend: Some(PrimerDesignBackend::Primer3),
            primer3_executable: Some(fake_primer3.clone()),
        },
    )
    .expect("primers design primer3");
    let primer3 = engine
        .get_primer_design_report("primer3_norm")
        .expect("primer3 report");

    assert_eq!(internal.template, primer3.template);
    assert_eq!(internal.roi_start_0based, primer3.roi_start_0based);
    assert_eq!(internal.roi_end_0based, primer3.roi_end_0based);
    assert_eq!(internal.min_amplicon_bp, primer3.min_amplicon_bp);
    assert_eq!(internal.max_amplicon_bp, primer3.max_amplicon_bp);
    assert_eq!(internal.pair_count, 1);
    assert_eq!(primer3.pair_count, 1);
    let internal_pair = internal.pairs.first().expect("internal pair");
    let primer3_pair = primer3.pairs.first().expect("primer3 pair");
    assert_eq!(
        internal_pair.forward.start_0based,
        primer3_pair.forward.start_0based
    );
    assert_eq!(
        internal_pair.forward.end_0based_exclusive,
        primer3_pair.forward.end_0based_exclusive
    );
    assert_eq!(
        internal_pair.reverse.start_0based,
        primer3_pair.reverse.start_0based
    );
    assert_eq!(
        internal_pair.reverse.end_0based_exclusive,
        primer3_pair.reverse.end_0based_exclusive
    );
    assert_eq!(
        internal_pair.amplicon_length_bp,
        primer3_pair.amplicon_length_bp
    );
    assert_eq!(primer3.backend.requested, "primer3");
    assert_eq!(primer3.backend.used, "primer3");
    assert_eq!(
        primer3.backend.primer3_executable.as_deref(),
        Some(fake_primer3.as_str())
    );
    assert_eq!(
        primer3.backend.primer3_version.as_deref(),
        Some("primer3_core synthetic-fixture 2.6.1")
    );
}

#[cfg(unix)]
#[test]
fn execute_primers_preflight_reports_reachable_primer3() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let tmp = tempdir().expect("tempdir");
    let fixture_path = Path::new(&primer3_fixture_path("pairs.location_5_60.kv")).to_path_buf();
    let fake_primer3 = install_fake_primer3(tmp.path(), &fixture_path);
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersPreflight {
            backend: Some(PrimerDesignBackend::Primer3),
            primer3_executable: Some(fake_primer3.clone()),
        },
    )
    .expect("primers preflight");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["schema"].as_str(),
        Some("gentle.primer3_preflight.v1")
    );
    assert_eq!(out.output["preflight"]["reachable"].as_bool(), Some(true));
    assert_eq!(
        out.output["preflight"]["version_probe_ok"].as_bool(),
        Some(true)
    );
    assert_eq!(out.output["preflight"]["backend"].as_str(), Some("primer3"));
    assert_eq!(
        out.output["preflight"]["executable"].as_str(),
        Some(fake_primer3.as_str())
    );
}

#[test]
fn execute_primers_design_qpcr_list_show_export() {
    let mut state = ProjectState::default();
    state.sequences.insert(
            "tpl".to_string(),
            DNAsequence::from_sequence(
                "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            )
            .expect("sequence"),
        );
    let mut engine = GentleEngine::from_state(state);
    let tmp = tempdir().expect("tempdir");
    let export_path = tmp.path().join("qpcr_report.json");
    let request = serde_json::to_string(&Operation::DesignQpcrAssays {
        template: "tpl".to_string(),
        roi_start_0based: 30,
        roi_end_0based: 70,
        forward: crate::engine::PrimerDesignSideConstraint {
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
            ..Default::default()
        },
        reverse: crate::engine::PrimerDesignSideConstraint {
            min_length: 20,
            max_length: 20,
            location_0based: Some(60),
            start_0based: None,
            end_0based: None,
            min_tm_c: 40.0,
            max_tm_c: 90.0,
            min_gc_fraction: 0.0,
            max_gc_fraction: 1.0,
            max_anneal_hits: 100,
            ..Default::default()
        },
        probe: crate::engine::PrimerDesignSideConstraint {
            min_length: 20,
            max_length: 20,
            location_0based: Some(35),
            start_0based: None,
            end_0based: None,
            min_tm_c: 40.0,
            max_tm_c: 90.0,
            min_gc_fraction: 0.0,
            max_gc_fraction: 1.0,
            max_anneal_hits: 100,
            ..Default::default()
        },
        min_amplicon_bp: 40,
        max_amplicon_bp: 130,
        pair_constraints: crate::engine::PrimerDesignPairConstraint::default(),
        max_tm_delta_c: Some(50.0),
        max_probe_tm_delta_c: Some(50.0),
        max_assays: Some(10),
        report_id: Some("tp73_qpcr".to_string()),
    })
    .expect("serialize request");

    let design = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersDesignQpcr {
            request_json: request,
            backend: Some(PrimerDesignBackend::Internal),
            primer3_executable: None,
        },
    )
    .expect("primers design-qpcr");
    assert!(design.state_changed);
    let report_id = design.output["report"]["report_id"]
        .as_str()
        .unwrap_or_default()
        .to_string();
    assert_eq!(report_id, "tp73_qpcr");

    let listed = execute_shell_command(&mut engine, &ShellCommand::PrimersListQpcrReports)
        .expect("primers list-qpcr-reports");
    assert!(!listed.state_changed);
    assert_eq!(
        listed.output["schema"].as_str(),
        Some("gentle.qpcr_design_report_list.v1")
    );
    assert_eq!(listed.output["report_count"].as_u64(), Some(1));
    assert!(
        listed.output["reports"][0]["best_assay_summary"]
            .as_str()
            .is_some_and(|value| !value.trim().is_empty())
    );

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersShowQpcrReport {
            report_id: report_id.clone(),
        },
    )
    .expect("primers show-qpcr-report");
    assert!(!shown.state_changed);
    assert_eq!(
        shown.output["report"]["report_id"].as_str(),
        Some("tp73_qpcr")
    );
    assert!(
        shown.output["report"]["best_assay_probe_placement"]
            .as_str()
            .is_some_and(|value| !value.is_empty())
    );
    assert!(
        shown.output["report"]["best_assay_summary"]
            .as_str()
            .is_some_and(|value| !value.trim().is_empty())
    );

    let exported = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersExportQpcrReport {
            report_id,
            path: export_path.to_string_lossy().to_string(),
        },
    )
    .expect("primers export-qpcr-report");
    assert!(!exported.state_changed);
    assert_eq!(
        exported.output["schema"].as_str(),
        Some("gentle.qpcr_design_report_export.v1")
    );
    let text = fs::read_to_string(&export_path).expect("read export");
    assert!(text.contains("gentle.qpcr_design_report.v1"));
}

#[test]
fn execute_primers_seed_from_feature_and_splicing() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), tp53_isoform_test_sequence());
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

    let seeded_feature = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersSeedFromFeature {
            seq_id: "seq_a".to_string(),
            feature_id,
        },
    )
    .expect("seed from feature");
    assert!(!seeded_feature.state_changed);
    assert_eq!(
        seeded_feature.output["schema"].as_str(),
        Some("gentle.primer_seed_request.v1")
    );
    assert_eq!(
        seeded_feature.output["source"]["kind"].as_str(),
        Some("feature")
    );
    let feature_start = seeded_feature.output["roi_start_0based"]
        .as_u64()
        .expect("feature start roi") as usize;
    let feature_end = seeded_feature.output["roi_end_0based_exclusive"]
        .as_u64()
        .expect("feature end roi") as usize;
    assert!(feature_end > feature_start);
    assert_eq!(
        seeded_feature.output["operations"]["design_primer_pairs"]["DesignPrimerPairs"]["template"]
            .as_str(),
        Some("seq_a")
    );
    assert_eq!(
        seeded_feature.output["operations"]["design_primer_pairs"]["DesignPrimerPairs"]
            ["roi_start_0based"]
            .as_u64(),
        Some(feature_start as u64)
    );
    assert_eq!(
        seeded_feature.output["operations"]["design_primer_pairs"]["DesignPrimerPairs"]
            ["roi_end_0based"]
            .as_u64(),
        Some(feature_end as u64)
    );

    let seeded_splicing = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersSeedFromSplicing {
            seq_id: "seq_a".to_string(),
            feature_id,
        },
    )
    .expect("seed from splicing");
    assert!(!seeded_splicing.state_changed);
    assert_eq!(
        seeded_splicing.output["schema"].as_str(),
        Some("gentle.primer_seed_request.v1")
    );
    assert_eq!(
        seeded_splicing.output["source"]["kind"].as_str(),
        Some("splicing")
    );
    let splicing_start = seeded_splicing.output["roi_start_0based"]
        .as_u64()
        .expect("splicing start roi") as usize;
    let splicing_end = seeded_splicing.output["roi_end_0based_exclusive"]
        .as_u64()
        .expect("splicing end roi") as usize;
    assert!(splicing_end > splicing_start);
    assert_eq!(splicing_start, feature_start);
    assert_eq!(splicing_end, feature_end);
    assert_eq!(
        seeded_splicing.output["operations"]["design_qpcr_assays"]["DesignQpcrAssays"]["template"]
            .as_str(),
        Some("seq_a")
    );

    let seeded_qpcr_feature = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersSeedQpcrFromFeature {
            seq_id: "seq_a".to_string(),
            feature_id,
        },
    )
    .expect("seed qpcr from feature");
    assert!(!seeded_qpcr_feature.state_changed);
    assert_eq!(
        seeded_qpcr_feature.output["schema"].as_str(),
        Some("gentle.qpcr_seed_request.v1")
    );
    assert_eq!(
        seeded_qpcr_feature.output["operation"]["DesignQpcrAssays"]["template"].as_str(),
        Some("seq_a")
    );
    assert_eq!(
        seeded_qpcr_feature.output["protocol_cartoon"]["protocol"].as_str(),
        Some("pcr.assay.qpcr")
    );
    assert!(
        seeded_qpcr_feature.output["rationale"]["summary"]
            .as_str()
            .is_some_and(|value| value.contains("Feature span reused"))
    );
    assert_eq!(
        seeded_qpcr_feature.output["rationale"]["recommended_defaults"]["max_assays"].as_u64(),
        Some(200)
    );

    let seeded_qpcr_splicing = execute_shell_command(
        &mut engine,
        &ShellCommand::PrimersSeedQpcrFromSplicing {
            seq_id: "seq_a".to_string(),
            feature_id,
        },
    )
    .expect("seed qpcr from splicing");
    assert!(!seeded_qpcr_splicing.state_changed);
    assert_eq!(
        seeded_qpcr_splicing.output["schema"].as_str(),
        Some("gentle.qpcr_seed_request.v1")
    );
    assert_eq!(
        seeded_qpcr_splicing.output["source"]["kind"].as_str(),
        Some("splicing")
    );
    assert_eq!(
        seeded_qpcr_splicing.output["operation"]["DesignQpcrAssays"]["template"].as_str(),
        Some("seq_a")
    );
    assert!(
        seeded_qpcr_splicing.output["rationale"]["summary"]
            .as_str()
            .is_some_and(|value| value.contains("Splicing-group region reused"))
    );
    assert!(
        seeded_qpcr_splicing.output["rationale"]["why_this_roi"]
            .as_str()
            .is_some_and(|value| value.contains("group"))
    );
}

#[test]
fn execute_primers_design_with_options_emits_primer_progress() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "tpl".to_string(),
        DNAsequence::from_sequence(&"ACGT".repeat(180)).unwrap(),
    );
    let mut engine = GentleEngine::from_state(state);
    let request = serde_json::to_string(&Operation::DesignPrimerPairs {
        template: "tpl".to_string(),
        roi_start_0based: 220,
        roi_end_0based: 320,
        forward: PrimerDesignSideConstraint {
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
        },
        reverse: PrimerDesignSideConstraint {
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
        },
        pair_constraints: PrimerDesignPairConstraint::default(),
        min_amplicon_bp: 80,
        max_amplicon_bp: 360,
        max_tm_delta_c: Some(100.0),
        max_pairs: Some(5),
        report_id: Some("shell_progress_pairs".to_string()),
    })
    .expect("serialize request");
    let progress_events = Arc::new(Mutex::new(Vec::<PrimerDesignProgress>::new()));
    let callback_events = Arc::clone(&progress_events);
    let progress_callback: ShellProgressCallback =
        Arc::new(Mutex::new(Box::new(move |progress: OperationProgress| {
            if let OperationProgress::PrimerDesign(p) = progress {
                callback_events
                    .lock()
                    .expect("progress events lock")
                    .push(p);
            }
            true
        })));

    let out = execute_shell_command_with_options(
        &mut engine,
        &ShellCommand::PrimersDesign {
            request_json: request,
            backend: Some(PrimerDesignBackend::Internal),
            primer3_executable: None,
        },
        &ShellExecutionOptions {
            allow_screenshots: false,
            allow_agent_commands: true,
            progress_callback: Some(progress_callback),
        },
    )
    .expect("execute primers design with progress");

    assert!(out.state_changed);
    let progress_events = progress_events.lock().expect("progress events lock");
    assert!(!progress_events.is_empty());
    assert!(
        progress_events
            .iter()
            .any(|progress| progress.stage == "forward_candidates")
    );
    assert_eq!(
        progress_events
            .last()
            .map(|progress| progress.stage.as_str()),
        Some("pair_search_complete")
    );
}

#[test]
fn execute_async_blast_start_and_status_reports_failure_for_missing_genome() {
    let _guard = BLAST_ASYNC_TEST_MUTEX.lock().expect("blast mutex");
    clear_blast_async_jobs_for_test();
    let mut engine = GentleEngine::new();
    let start = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceBlastAsyncStart {
            helper_mode: true,
            genome_id: "missing_helper".to_string(),
            query_sequence: "ACGTACGT".to_string(),
            max_hits: 5,
            max_hits_explicit: true,
            task: Some("blastn-short".to_string()),
            request_options_json: None,
            catalog_path: None,
            cache_dir: None,
        },
    )
    .expect("start async blast");
    assert!(start.state_changed);
    assert_eq!(
        start.output["binary_preflight"]["schema"].as_str(),
        Some("gentle.blast_external_binary_preflight.v1")
    );
    let job_id = start
        .output
        .get("job")
        .and_then(|job| job.get("job_id"))
        .and_then(|value| value.as_str())
        .unwrap_or_default()
        .to_string();
    assert!(!job_id.is_empty());

    let mut terminal_state = String::new();
    for _ in 0..40 {
        let status = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceBlastAsyncStatus {
                helper_mode: true,
                job_id: job_id.clone(),
                include_report: true,
            },
        )
        .expect("status");
        terminal_state = status
            .output
            .get("job")
            .and_then(|job| job.get("state"))
            .and_then(|value| value.as_str())
            .unwrap_or_default()
            .to_string();
        if matches!(
            terminal_state.as_str(),
            "completed" | "failed" | "cancelled"
        ) {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(20));
    }
    assert!(
        matches!(terminal_state.as_str(), "failed" | "cancelled"),
        "unexpected async blast terminal state: {}",
        terminal_state
    );
    clear_blast_async_jobs_for_test();
}

#[test]
fn execute_export_run_bundle_writes_schema_json() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(crate::engine::Operation::Reverse {
            input: "seqA".to_string(),
            output_id: Some("seqA_rev".to_string()),
        })
        .expect("reverse");

    let tmp = tempfile::NamedTempFile::new().expect("tmp file");
    let path = tmp.path().with_extension("shell.run_bundle.json");
    let path_text = path.display().to_string();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ExportRunBundle {
            output: path_text.clone(),
            run_id: Some("interactive".to_string()),
        },
    )
    .expect("export run bundle");
    assert!(!out.state_changed);
    let text = fs::read_to_string(path_text).expect("read bundle");
    let value: serde_json::Value = serde_json::from_str(&text).expect("parse bundle");
    assert_eq!(
        value.get("schema").and_then(|v| v.as_str()),
        Some("gentle.process_run_bundle.v1")
    );
    assert_eq!(
        value
            .get("run_id_filter")
            .and_then(|v| v.as_str())
            .unwrap_or_default(),
        "interactive"
    );
}

#[test]
fn execute_export_run_bundle_matches_engine_decision_traces() {
    let td = tempdir().expect("tempdir");
    let shell_path = td.path().join("shell.run_bundle.json");
    let engine_path = td.path().join("engine.run_bundle.json");

    let mut shell_engine = GentleEngine::from_state(decision_trace_fixture_state());
    let shell_out = execute_shell_command(
        &mut shell_engine,
        &ShellCommand::ExportRunBundle {
            output: shell_path.to_string_lossy().to_string(),
            run_id: None,
        },
    )
    .expect("shell export run bundle");
    assert!(!shell_out.state_changed);

    let shell_bundle_text = fs::read_to_string(&shell_path).expect("read shell bundle output");
    let shell_bundle: crate::engine::ProcessRunBundleExport =
        serde_json::from_str(&shell_bundle_text).expect("parse shell bundle");

    let mut engine = GentleEngine::from_state(decision_trace_fixture_state());
    engine
        .apply(crate::engine::Operation::ExportProcessRunBundle {
            path: engine_path.to_string_lossy().to_string(),
            run_id: None,
        })
        .expect("engine export run bundle");
    let engine_bundle_text = fs::read_to_string(&engine_path).expect("read engine bundle");
    let engine_bundle: crate::engine::ProcessRunBundleExport =
        serde_json::from_str(&engine_bundle_text).expect("parse engine bundle");

    assert_eq!(
        serde_json::to_value(&shell_bundle.decision_traces).expect("serialize shell traces"),
        serde_json::to_value(&engine_bundle.decision_traces).expect("serialize engine traces")
    );
    assert_eq!(shell_bundle.decision_traces.len(), 1);
    assert_eq!(
        shell_bundle.decision_traces[0].status,
        "preflight_failed".to_string()
    );
    assert_eq!(
        shell_bundle.decision_traces[0]
            .preflight_snapshot
            .as_ref()
            .expect("snapshot derived from history")
            .errors,
        vec!["missing sequence".to_string()]
    );
}

#[test]
fn execute_export_run_bundle_matches_engine_construct_reasoning_block() {
    let td = tempdir().expect("tempdir");
    let shell_path = td.path().join("shell.run_bundle.reasoning.json");
    let engine_path = td.path().join("engine.run_bundle.reasoning.json");
    let fixture_state = decision_trace_with_construct_reasoning_fixture_state();

    let mut shell_engine = GentleEngine::from_state(fixture_state.clone());
    crate::engine::Engine::apply(
        &mut shell_engine,
        crate::engine::Operation::Reverse {
            input: "s".to_string(),
            output_id: Some("s_reasoning_rev".to_string()),
        },
    )
    .expect("shell fixture reverse");
    let shell_out = execute_shell_command(
        &mut shell_engine,
        &ShellCommand::ExportRunBundle {
            output: shell_path.to_string_lossy().to_string(),
            run_id: None,
        },
    )
    .expect("shell export run bundle");
    assert!(!shell_out.state_changed);

    let shell_bundle_text = fs::read_to_string(&shell_path).expect("read shell bundle output");
    let shell_bundle: crate::engine::ProcessRunBundleExport =
        serde_json::from_str(&shell_bundle_text).expect("parse shell bundle");

    let mut engine = GentleEngine::from_state(fixture_state);
    crate::engine::Engine::apply(
        &mut engine,
        crate::engine::Operation::Reverse {
            input: "s".to_string(),
            output_id: Some("s_reasoning_rev".to_string()),
        },
    )
    .expect("engine fixture reverse");
    engine
        .apply(crate::engine::Operation::ExportProcessRunBundle {
            path: engine_path.to_string_lossy().to_string(),
            run_id: None,
        })
        .expect("engine export run bundle");
    let engine_bundle_text = fs::read_to_string(&engine_path).expect("read engine bundle");
    let engine_bundle: crate::engine::ProcessRunBundleExport =
        serde_json::from_str(&engine_bundle_text).expect("parse engine bundle");

    assert_eq!(
        serde_json::to_value(&shell_bundle.construct_reasoning)
            .expect("serialize shell construct reasoning"),
        serde_json::to_value(&engine_bundle.construct_reasoning)
            .expect("serialize engine construct reasoning")
    );
    assert_eq!(shell_bundle.construct_reasoning.graphs.len(), 1);
    assert_eq!(shell_bundle.construct_reasoning.summaries.len(), 1);
    assert_eq!(
        shell_bundle.construct_reasoning.summaries[0]
            .fact_statuses
            .get("selection_context")
            .map(String::as_str),
        Some("supported")
    );
    assert!(
        shell_bundle.construct_reasoning.summaries[0]
            .supported_selection_rule_ids
            .iter()
            .any(|rule_id| rule_id == "ampicillin_vector_selection")
    );
}

#[test]
fn execute_async_blast_start_queues_when_capacity_is_reached() {
    with_blast_async_test_overrides(1, 2000, || {
        let mut engine = GentleEngine::new();
        let start_one = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceBlastAsyncStart {
                helper_mode: true,
                genome_id: "missing_helper".to_string(),
                query_sequence: "ACGTACGT".to_string(),
                max_hits: 5,
                max_hits_explicit: true,
                task: Some("blastn-short".to_string()),
                request_options_json: None,
                catalog_path: None,
                cache_dir: None,
            },
        )
        .expect("start first async blast");
        let job_one = start_one.output["job"]["job_id"]
            .as_str()
            .unwrap_or_default()
            .to_string();
        assert!(!job_one.is_empty());
        assert_eq!(
            start_one.output["job"]["max_concurrent_jobs"].as_u64(),
            Some(1)
        );
        let mut first_running = false;
        for _ in 0..40 {
            let status = execute_shell_command(
                &mut engine,
                &ShellCommand::ReferenceBlastAsyncStatus {
                    helper_mode: true,
                    job_id: job_one.clone(),
                    include_report: false,
                },
            )
            .expect("status first job before queue test");
            if status.output["job"]["state"].as_str() == Some("running") {
                first_running = true;
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        assert!(
            first_running,
            "expected first async BLAST job to occupy the only scheduler slot"
        );

        let start_two = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceBlastAsyncStart {
                helper_mode: true,
                genome_id: "missing_helper".to_string(),
                query_sequence: "TTTTAAAA".to_string(),
                max_hits: 5,
                max_hits_explicit: true,
                task: Some("blastn-short".to_string()),
                request_options_json: None,
                catalog_path: None,
                cache_dir: None,
            },
        )
        .expect("start second async blast");
        let job_two = start_two.output["job"]["job_id"]
            .as_str()
            .unwrap_or_default()
            .to_string();
        assert!(!job_two.is_empty());
        let start_two_state = start_two.output["job"]["state"].as_str();
        assert!(
            matches!(start_two_state, Some("queued") | Some("running")),
            "unexpected second-job state: {:?}",
            start_two_state
        );
        if start_two_state == Some("queued") {
            assert_eq!(start_two.output["job"]["queue_position"].as_u64(), Some(1));
        } else {
            assert_eq!(start_two.output["job"]["queue_position"].as_u64(), None);
        }

        let _cancel = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceBlastAsyncCancel {
                helper_mode: true,
                job_id: job_one,
            },
        )
        .expect("cancel first job");

        let mut observed_non_queued = start_two_state == Some("running");
        if !observed_non_queued {
            for _ in 0..80 {
                let status = execute_shell_command(
                    &mut engine,
                    &ShellCommand::ReferenceBlastAsyncStatus {
                        helper_mode: true,
                        job_id: job_two.clone(),
                        include_report: false,
                    },
                )
                .expect("status second job");
                let state = status.output["job"]["state"]
                    .as_str()
                    .unwrap_or_default()
                    .to_string();
                if state != "queued" {
                    observed_non_queued = true;
                    break;
                }
                std::thread::sleep(std::time::Duration::from_millis(15));
            }
        }
        assert!(
            observed_non_queued,
            "expected second job to leave the queue once scheduler capacity was available"
        );
    });
}

#[test]
fn execute_async_blast_list_reports_scheduler_metadata() {
    with_blast_async_test_overrides(2, 0, || {
        let mut engine = GentleEngine::new();
        let start = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceBlastAsyncStart {
                helper_mode: true,
                genome_id: "missing_helper".to_string(),
                query_sequence: "ACGTACGT".to_string(),
                max_hits: 5,
                max_hits_explicit: true,
                task: Some("blastn-short".to_string()),
                request_options_json: None,
                catalog_path: None,
                cache_dir: None,
            },
        )
        .expect("start async blast");
        assert!(start.state_changed);

        let listed = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceBlastAsyncList { helper_mode: true },
        )
        .expect("list async blast jobs");
        assert_eq!(listed.output["job_count"].as_u64(), Some(1));
        let jobs = listed.output["jobs"]
            .as_array()
            .expect("jobs should be array");
        assert_eq!(jobs.len(), 1);
        assert_eq!(jobs[0]["max_concurrent_jobs"].as_u64(), Some(2));
        assert!(jobs[0]["running_jobs"].as_u64().unwrap_or(0) <= 2);
        assert!(jobs[0]["queued_jobs"].as_u64().unwrap_or(0) <= 1);
    });
}

#[test]
fn execute_async_blast_restart_recovery_marks_orphaned_nonterminal_job_failed() {
    with_blast_async_test_overrides(1, 200, || {
        let mut engine = GentleEngine::new();
        let start = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceBlastAsyncStart {
                helper_mode: true,
                genome_id: "missing_helper".to_string(),
                query_sequence: "ACGTACGT".to_string(),
                max_hits: 5,
                max_hits_explicit: true,
                task: Some("blastn-short".to_string()),
                request_options_json: None,
                catalog_path: None,
                cache_dir: None,
            },
        )
        .expect("start async blast");
        assert!(start.state_changed);
        let job_id = start.output["job"]["job_id"]
            .as_str()
            .unwrap_or_default()
            .to_string();
        assert!(!job_id.is_empty());

        let persisted_state = engine.state().clone();
        {
            let mut jobs = BLAST_ASYNC_JOBS
                .lock()
                .expect("blast async registry lock for restart simulation");
            jobs.clear();
        }

        let mut restarted_engine = GentleEngine::from_state(persisted_state);
        let status = execute_shell_command(
            &mut restarted_engine,
            &ShellCommand::ReferenceBlastAsyncStatus {
                helper_mode: true,
                job_id: job_id.clone(),
                include_report: false,
            },
        )
        .expect("status after restart");
        assert!(status.state_changed);
        assert_eq!(
            status.output["job"]["job_id"].as_str(),
            Some(job_id.as_str())
        );
        assert_eq!(status.output["job"]["state"].as_str(), Some("failed"));
        assert_eq!(
            status.output["job"]["error"].as_str(),
            Some(BLAST_ASYNC_RESTART_INTERRUPTED_ERROR)
        );
        assert!(
            !status.output["job"]["result_available"]
                .as_bool()
                .unwrap_or(true)
        );
    });
}

#[test]
fn execute_async_blast_cancel_after_restart_is_deterministic() {
    with_blast_async_test_overrides(1, 200, || {
        let mut engine = GentleEngine::new();
        let start = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceBlastAsyncStart {
                helper_mode: true,
                genome_id: "missing_helper".to_string(),
                query_sequence: "ACGTACGT".to_string(),
                max_hits: 5,
                max_hits_explicit: true,
                task: Some("blastn-short".to_string()),
                request_options_json: None,
                catalog_path: None,
                cache_dir: None,
            },
        )
        .expect("start async blast");
        assert!(start.state_changed);
        let job_id = start.output["job"]["job_id"]
            .as_str()
            .unwrap_or_default()
            .to_string();
        assert!(!job_id.is_empty());

        let persisted_state = engine.state().clone();
        {
            let mut jobs = BLAST_ASYNC_JOBS
                .lock()
                .expect("blast async registry lock for restart simulation");
            jobs.clear();
        }

        let mut restarted_engine = GentleEngine::from_state(persisted_state);
        let cancel = execute_shell_command(
            &mut restarted_engine,
            &ShellCommand::ReferenceBlastAsyncCancel {
                helper_mode: true,
                job_id: job_id.clone(),
            },
        )
        .expect("cancel after restart");
        assert!(cancel.state_changed);
        assert_eq!(
            cancel.output["job"]["job_id"].as_str(),
            Some(job_id.as_str())
        );
        assert_eq!(cancel.output["job"]["state"].as_str(), Some("cancelled"));
        assert_eq!(
            cancel.output["job"]["error"].as_str(),
            Some("BLAST search cancelled by caller")
        );

        let status = execute_shell_command(
            &mut restarted_engine,
            &ShellCommand::ReferenceBlastAsyncStatus {
                helper_mode: true,
                job_id,
                include_report: false,
            },
        )
        .expect("status after restart-cancel");
        assert_eq!(status.output["job"]["state"].as_str(), Some("cancelled"));
    });
}

#[test]
fn execute_async_blast_restart_snapshot_race_still_normalizes_deterministically() {
    with_blast_async_test_overrides(1, 80, || {
        let mut engine = GentleEngine::new();
        let start = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceBlastAsyncStart {
                helper_mode: true,
                genome_id: "missing_helper".to_string(),
                query_sequence: "ACGTACGT".to_string(),
                max_hits: 5,
                max_hits_explicit: true,
                task: Some("blastn-short".to_string()),
                request_options_json: None,
                catalog_path: None,
                cache_dir: None,
            },
        )
        .expect("start async blast");
        let job_id = start.output["job"]["job_id"]
            .as_str()
            .unwrap_or_default()
            .to_string();
        assert!(!job_id.is_empty());
        let pre_completion_snapshot = engine.state().clone();

        let mut terminal_state = String::new();
        for _ in 0..60 {
            let status = execute_shell_command(
                &mut engine,
                &ShellCommand::ReferenceBlastAsyncStatus {
                    helper_mode: true,
                    job_id: job_id.clone(),
                    include_report: false,
                },
            )
            .expect("status before restart");
            terminal_state = status.output["job"]["state"]
                .as_str()
                .unwrap_or_default()
                .to_string();
            if matches!(
                terminal_state.as_str(),
                "completed" | "failed" | "cancelled"
            ) {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(15));
        }
        assert!(
            matches!(terminal_state.as_str(), "failed" | "cancelled"),
            "expected original process to reach terminal state, got {}",
            terminal_state
        );

        {
            let mut jobs = BLAST_ASYNC_JOBS
                .lock()
                .expect("blast async registry lock for restart simulation");
            jobs.clear();
        }

        let mut restarted_engine = GentleEngine::from_state(pre_completion_snapshot);
        let status = execute_shell_command(
            &mut restarted_engine,
            &ShellCommand::ReferenceBlastAsyncStatus {
                helper_mode: true,
                job_id,
                include_report: false,
            },
        )
        .expect("status from stale snapshot after restart");
        assert_eq!(status.output["job"]["state"].as_str(), Some("failed"));
        assert_eq!(
            status.output["job"]["error"].as_str(),
            Some(BLAST_ASYNC_RESTART_INTERRUPTED_ERROR)
        );
    });
}

#[test]
fn execute_workflow_macro_transactional_rolls_back_on_error() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let err = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosRun {
            script: r#"op {"Reverse":{"input":"seqA","output_id":"tmp_rev"}}
op {"Reverse":{"input":"missing","output_id":"bad"}}"#
                .to_string(),
            transactional: true,
        },
    )
    .unwrap_err();
    assert!(err.contains("rolled back"));
    assert!(!engine.state().sequences.contains_key("tmp_rev"));
}

#[test]
fn execute_macros_run_records_inline_macro_instance() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosRun {
            script: r#"op {"Reverse":{"input":"seqA","output_id":"seqA_rev"}}"#.to_string(),
            transactional: false,
        },
    )
    .expect("run inline macro");
    assert!(run.state_changed);
    assert!(run.output["macro_instance_id"].as_str().is_some());
    assert_eq!(engine.state().lineage.macro_instances.len(), 1);
    let instance = &engine.state().lineage.macro_instances[0];
    assert!(instance.template_name.is_none());
    assert!(!instance.expanded_op_ids.is_empty());
}

#[test]
fn execute_macros_run_failed_records_macro_instance() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let err = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosRun {
            script: r#"op {"Reverse":{"input":"missing","output_id":"bad"}}"#.to_string(),
            transactional: false,
        },
    )
    .expect_err("macro should fail");
    assert!(err.contains("macro_instance_id="));
    assert_eq!(engine.state().lineage.macro_instances.len(), 1);
    let instance = &engine.state().lineage.macro_instances[0];
    assert_eq!(instance.status, MacroInstanceStatus::Failed);
    assert!(instance.status_message.is_some());
    assert!(instance.expanded_op_ids.is_empty());
}

#[test]
fn execute_macros_instance_list_and_show() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seqA".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosRun {
            script: r#"op {"Reverse":{"input":"seqA","output_id":"seqA_rev"}}"#.to_string(),
            transactional: false,
        },
    )
    .expect("run macro");

    let listed = execute_shell_command(&mut engine, &ShellCommand::MacrosInstanceList)
        .expect("list macro instances");
    assert!(!listed.state_changed);
    assert_eq!(
        listed.output["schema"].as_str(),
        Some("gentle.lineage_macro_instances.v1")
    );
    let instance_id = listed
        .output
        .get("instances")
        .and_then(|rows| rows.as_array())
        .and_then(|rows| rows.first())
        .and_then(|row| row.get("macro_instance_id"))
        .and_then(|v| v.as_str())
        .unwrap_or_default()
        .to_string();
    assert!(!instance_id.is_empty());

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::MacrosInstanceShow {
            macro_instance_id: instance_id,
        },
    )
    .expect("show macro instance");
    assert!(!shown.state_changed);
    assert_eq!(
        shown.output["schema"].as_str(),
        Some("gentle.lineage_macro_instance.v1")
    );
    assert!(
        shown
            .output
            .get("instance")
            .and_then(|row| row.get("expanded_op_ids"))
            .and_then(|v| v.as_array())
            .map(|rows| !rows.is_empty())
            .unwrap_or(false)
    );
}

#[test]
fn execute_tracks_tracked_add_and_list() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let subscription = GenomeTrackSubscription {
        source: GenomeTrackSource::Bed,
        path: "test_files/data/peaks.bed.gz".to_string(),
        track_name: Some("ChIP".to_string()),
        min_score: Some(1.0),
        max_score: Some(10.0),
        clear_existing: true,
    };
    let add = execute_shell_command(
        &mut engine,
        &ShellCommand::TracksTrackedAdd {
            subscription: subscription.clone(),
        },
    )
    .expect("execute tracked add");
    assert!(add.state_changed);
    assert_eq!(add.output["inserted"].as_bool(), Some(true));
    assert_eq!(add.output["count"].as_u64(), Some(1));

    let list = execute_shell_command(&mut engine, &ShellCommand::TracksTrackedList)
        .expect("execute tracked list");
    assert!(!list.state_changed);
    assert_eq!(list.output["count"].as_u64(), Some(1));
}

#[test]
fn parse_agents_ask_command() {
    let cmd = parse_shell_line(
            "agents ask builtin_echo --prompt 'auto: state-summary' --catalog catalog.json --base-url http://localhost:11964 --model deepseek-r1:8b --timeout-secs 600 --connect-timeout-secs 20 --read-timeout-secs 900 --max-retries 4 --max-response-bytes 2097152 --allow-auto-exec --execute-index 2 --no-state-summary",
        )
        .expect("parse agents ask");
    match cmd {
        ShellCommand::AgentsAsk {
            system_id,
            prompt,
            catalog_path,
            base_url_override,
            model_override,
            timeout_seconds,
            connect_timeout_seconds,
            read_timeout_seconds,
            max_retries,
            max_response_bytes,
            include_state_summary,
            allow_auto_exec,
            execute_all,
            execute_indices,
        } => {
            assert_eq!(system_id, "builtin_echo");
            assert_eq!(prompt, "auto: state-summary");
            assert_eq!(catalog_path.as_deref(), Some("catalog.json"));
            assert_eq!(base_url_override.as_deref(), Some("http://localhost:11964"));
            assert_eq!(model_override.as_deref(), Some("deepseek-r1:8b"));
            assert_eq!(timeout_seconds, Some(600));
            assert_eq!(connect_timeout_seconds, Some(20));
            assert_eq!(read_timeout_seconds, Some(900));
            assert_eq!(max_retries, Some(4));
            assert_eq!(max_response_bytes, Some(2_097_152));
            assert!(!include_state_summary);
            assert!(allow_auto_exec);
            assert!(!execute_all);
            assert_eq!(execute_indices, vec![2]);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_agents_ask_runs_auto_suggestion_when_enabled() {
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("agents.json");
    let catalog_json = r#"{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    {
      "id": "builtin_echo",
      "label": "Builtin Echo",
      "transport": "builtin_echo"
    }
  ]
}"#;
    fs::write(&catalog_path, catalog_json).expect("write catalog");

    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::AgentsAsk {
            system_id: "builtin_echo".to_string(),
            prompt: "auto: state-summary\nask: capabilities".to_string(),
            catalog_path: Some(catalog_path.display().to_string()),
            base_url_override: None,
            model_override: None,
            timeout_seconds: None,
            connect_timeout_seconds: None,
            read_timeout_seconds: None,
            max_retries: None,
            max_response_bytes: None,
            include_state_summary: true,
            allow_auto_exec: true,
            execute_all: false,
            execute_indices: vec![],
        },
    )
    .expect("execute agents ask");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["summary"]["suggested_command_count"].as_u64(),
        Some(2)
    );
    assert_eq!(out.output["summary"]["executed_count"].as_u64(), Some(1));
}

#[test]
fn execute_agents_ask_rejected_when_context_disallows_agent_commands() {
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("agents.json");
    let catalog_json = r#"{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    {
      "id": "builtin_echo",
      "label": "Builtin Echo",
      "transport": "builtin_echo"
    }
  ]
}"#;
    fs::write(&catalog_path, catalog_json).expect("write catalog");
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let err = execute_shell_command_with_options(
        &mut engine,
        &ShellCommand::AgentsAsk {
            system_id: "builtin_echo".to_string(),
            prompt: "ask: state-summary".to_string(),
            catalog_path: Some(catalog_path.display().to_string()),
            base_url_override: None,
            model_override: None,
            timeout_seconds: None,
            connect_timeout_seconds: None,
            read_timeout_seconds: None,
            max_retries: None,
            max_response_bytes: None,
            include_state_summary: true,
            allow_auto_exec: false,
            execute_all: false,
            execute_indices: vec![],
        },
        &ShellExecutionOptions {
            allow_screenshots: false,
            allow_agent_commands: false,
            progress_callback: None,
        },
    )
    .expect_err("agents ask should be blocked in this execution context");
    assert!(
        err.contains("agent-to-agent recursion guardrail"),
        "unexpected error: {err}"
    );
}

#[test]
fn execute_agents_ask_blocks_nested_agent_call_inside_macro_suggestion() {
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("agents.json");
    let catalog_json = r#"{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    {
      "id": "builtin_echo",
      "label": "Builtin Echo",
      "transport": "builtin_echo"
    }
  ]
}"#;
    fs::write(&catalog_path, catalog_json).expect("write catalog");

    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::AgentsAsk {
            system_id: "builtin_echo".to_string(),
            prompt: "auto: macros run 'agents ask builtin_echo --prompt \"ask: capabilities\"'"
                .to_string(),
            catalog_path: Some(catalog_path.display().to_string()),
            base_url_override: None,
            model_override: None,
            timeout_seconds: None,
            connect_timeout_seconds: None,
            read_timeout_seconds: None,
            max_retries: None,
            max_response_bytes: None,
            include_state_summary: true,
            allow_auto_exec: true,
            execute_all: false,
            execute_indices: vec![],
        },
    )
    .expect("execute agents ask");
    assert!(!out.state_changed);
    assert_eq!(out.output["summary"]["executed_count"].as_u64(), Some(1));
    assert_eq!(
        out.output["summary"]["executed_error_count"].as_u64(),
        Some(1)
    );
    assert_eq!(out.output["executions"][0]["ok"].as_bool(), Some(false));
    let error = out.output["executions"][0]["error"]
        .as_str()
        .expect("error string");
    assert!(
        error.contains("agent-to-agent recursion guardrail"),
        "unexpected error: {error}"
    );
}

#[test]
fn execute_macros_run_blocks_nested_agents_ask_when_agent_commands_are_disabled() {
    let tmp = tempdir().expect("tempdir");
    let catalog_path = tmp.path().join("agents.json");
    let catalog_json = r#"{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    {
      "id": "builtin_echo",
      "label": "Builtin Echo",
      "transport": "builtin_echo"
    }
  ]
}"#;
    fs::write(&catalog_path, catalog_json).expect("write catalog");

    let mut engine = GentleEngine::from_state(ProjectState::default());
    let err = execute_shell_command_with_options(
        &mut engine,
        &ShellCommand::MacrosRun {
            script: format!(
                "agents ask builtin_echo --catalog {} --prompt \"ask: capabilities\"",
                catalog_path.display()
            ),
            transactional: false,
        },
        &ShellExecutionOptions {
            allow_agent_commands: false,
            ..ShellExecutionOptions::default()
        },
    )
    .expect_err("nested agents ask should be blocked inside macros run");
    assert!(
        err.contains("agent-to-agent recursion guardrail"),
        "unexpected error: {err}"
    );
}

#[test]
fn execute_agent_suggestions_allows_blast_shell_route() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let suggestions = vec![crate::agent_bridge::AgentSuggestedCommand {
        title: Some("BLAST".to_string()),
        rationale: Some("Check specificity".to_string()),
        command: "genomes blast __missing_genome__ ACGT".to_string(),
        execution: AgentExecutionIntent::Auto,
    }];
    let execute_indices = std::collections::BTreeSet::new();
    let (changed, reports) = execute_agent_suggested_commands(
        &mut engine,
        &suggestions,
        false,
        &execute_indices,
        true,
        &ShellExecutionOptions::default(),
    );
    assert!(!changed);
    assert_eq!(reports.len(), 1);
    let row = &reports[0];
    assert!(row.executed);
    assert!(!row.ok);
    let error = row
        .error
        .as_deref()
        .unwrap_or_default()
        .to_ascii_lowercase();
    assert!(
        !error.contains("agent-to-agent"),
        "BLAST route should not be blocked by recursion guardrail: {error}"
    );
    assert!(
        !error.contains("blocked in this context"),
        "BLAST route should be executable in suggested-command context: {error}"
    );
}

#[test]
fn execute_state_summary_returns_json() {
    let mut state = ProjectState::default();
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Measured sample".to_string()),
            members: vec![],
            declared_contents_exclusive: false,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);
    let out = execute_shell_command(&mut engine, &ShellCommand::StateSummary)
        .expect("execute state summary");
    assert!(!out.state_changed);
    assert!(out.output.get("sequence_count").is_some());
    assert_eq!(
        out.output["containers"][0]["declared_contents_exclusive"],
        serde_json::json!(false)
    );
}

#[test]
fn parse_containers_set_exclusive_command() {
    let cmd = parse_shell_line("containers set-exclusive container-1 false")
        .expect("parse containers set-exclusive");
    match cmd {
        ShellCommand::ContainersSetExclusive {
            container_id,
            exclusive,
        } => {
            assert_eq!(container_id, "container-1".to_string());
            assert!(!exclusive);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_containers_set_exclusive_updates_state() {
    let mut state = ProjectState::default();
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Measured sample".to_string()),
            members: vec![],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut engine = GentleEngine::from_state(state);
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ContainersSetExclusive {
            container_id: "container-1".to_string(),
            exclusive: false,
        },
    )
    .expect("execute containers set-exclusive");
    assert!(out.state_changed);
    assert_eq!(
        engine
            .state()
            .container_state
            .containers
            .get("container-1")
            .expect("container")
            .declared_contents_exclusive,
        false
    );
}

#[test]
fn parse_resources_sync_rebase_with_flag() {
    let cmd = parse_shell_line(
        "resources sync-rebase https://example.org/rebase.withrefm out.json --commercial-only",
    )
    .expect("parse resources sync-rebase");
    match cmd {
        ShellCommand::ResourcesSyncRebase {
            input,
            output,
            commercial_only,
        } => {
            assert_eq!(input, "https://example.org/rebase.withrefm".to_string());
            assert_eq!(output, Some("out.json".to_string()));
            assert!(commercial_only);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_resources_sync_jaspar_with_output() {
    let cmd = parse_shell_line("resources sync-jaspar motifs.pfm out.motifs.json")
        .expect("parse resources sync-jaspar");
    match cmd {
        ShellCommand::ResourcesSyncJaspar { input, output } => {
            assert_eq!(input, "motifs.pfm".to_string());
            assert_eq!(output, Some("out.motifs.json".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_resources_summarize_jaspar_with_options() {
    let cmd = parse_shell_line(
        "resources summarize-jaspar --motif SP1 --motifs REST,PATZ1 --random-length 512 --seed 99 --output jaspar.json",
    )
    .expect("parse resources summarize-jaspar");
    match cmd {
        ShellCommand::ResourcesSummarizeJaspar {
            motifs,
            random_sequence_length_bp,
            random_seed,
            output,
        } => {
            assert_eq!(
                motifs,
                vec!["SP1".to_string(), "REST".to_string(), "PATZ1".to_string()]
            );
            assert_eq!(random_sequence_length_bp, 512);
            assert_eq!(random_seed, 99);
            assert_eq!(output.as_deref(), Some("jaspar.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_resources_benchmark_jaspar_with_options() {
    let cmd = parse_shell_line(
        "resources benchmark-jaspar --random-length 4096 --seed 11 --output jaspar.benchmark.json",
    )
    .expect("parse resources benchmark-jaspar");
    match cmd {
        ShellCommand::ResourcesBenchmarkJaspar {
            random_sequence_length_bp,
            random_seed,
            output,
        } => {
            assert_eq!(random_sequence_length_bp, 4096);
            assert_eq!(random_seed, 11);
            assert_eq!(output.as_deref(), Some("jaspar.benchmark.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_resources_sync_attract_with_output() {
    let cmd = parse_shell_line("resources sync-attract ATtRACT.zip out.attract.json")
        .expect("parse resources sync-attract");
    match cmd {
        ShellCommand::ResourcesSyncAttract { input, output } => {
            assert_eq!(input, "ATtRACT.zip".to_string());
            assert_eq!(output, Some("out.attract.json".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_resources_inspect_jaspar_with_options() {
    let cmd = parse_shell_line(
        "resources inspect-jaspar SP1 --random-length 2048 --seed 5 --fetch-remote --output jaspar.expert.json",
    )
    .expect("parse resources inspect-jaspar");
    match cmd {
        ShellCommand::ResourcesInspectJaspar {
            motif,
            random_sequence_length_bp,
            random_seed,
            fetch_remote,
            output,
        } => {
            assert_eq!(motif, "SP1".to_string());
            assert_eq!(random_sequence_length_bp, 2048);
            assert_eq!(random_seed, 5);
            assert!(fetch_remote);
            assert_eq!(output.as_deref(), Some("jaspar.expert.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_resources_list_jaspar_with_options() {
    let cmd = parse_shell_line(
        "resources list-jaspar --filter TP --limit 25 --fetch-remote --output jaspar.catalog.json",
    )
    .expect("parse resources list-jaspar");
    match cmd {
        ShellCommand::ResourcesListJaspar {
            filter,
            limit,
            fetch_remote,
            output,
        } => {
            assert_eq!(filter.as_deref(), Some("TP"));
            assert_eq!(limit, Some(25));
            assert!(fetch_remote);
            assert_eq!(output.as_deref(), Some("jaspar.catalog.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_resources_sync_jaspar_remote_metadata_with_options() {
    let cmd = parse_shell_line(
        "resources sync-jaspar-remote-metadata --motif SP1 --motif REST --limit 20 --output jaspar.remote.json",
    )
    .expect("parse resources sync-jaspar-remote-metadata");
    match cmd {
        ShellCommand::ResourcesSyncJasparRemoteMetadata {
            motifs,
            filter,
            limit,
            output,
        } => {
            assert_eq!(motifs, vec!["SP1".to_string(), "REST".to_string()]);
            assert_eq!(filter, None);
            assert_eq!(limit, Some(20));
            assert_eq!(output.as_deref(), Some("jaspar.remote.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_resources_status() {
    let cmd = parse_shell_line("resources status").expect("parse resources status");
    match cmd {
        ShellCommand::ResourcesStatus => {}
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_services_status() {
    let cmd = parse_shell_line("services status").expect("parse services status");
    match cmd {
        ShellCommand::ServicesStatus => {}
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_services_handoff_with_scope_and_output() {
    let cmd = parse_shell_line(
        "services handoff --scope clawbio --output artifacts/service_handoff.json",
    )
    .expect("parse services handoff");
    match cmd {
        ShellCommand::ServicesHandoff { scope, output } => {
            assert_eq!(scope.as_deref(), Some("clawbio"));
            assert_eq!(output.as_deref(), Some("artifacts/service_handoff.json"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_attract_inspect_splicing() {
    let cmd = parse_shell_line(
        "attract inspect-splicing seq_1 5 --scope target_group_target_strand --organism \"Homo sapiens\" --flank-bp 12 --min-score 3.0 --min-match-quantile 0.975 --pwm-mapping windowed_submatrix --compare-policies --all-transcripts --no-fallback",
    )
    .expect("parse attract inspect-splicing");
    match cmd {
        ShellCommand::InspectSplicingAttract {
            seq_id,
            feature_id,
            settings,
        } => {
            assert_eq!(seq_id, "seq_1".to_string());
            assert_eq!(feature_id, 5);
            assert_eq!(settings.scope, SplicingScopePreset::TargetGroupTargetStrand);
            assert_eq!(settings.requested_organism.as_deref(), Some("Homo sapiens"));
            assert_eq!(settings.boundary_flank_bp, 12);
            assert_eq!(settings.minimum_quality_score, 3.0);
            assert_eq!(settings.minimum_match_quantile, 0.975);
            assert_eq!(
                settings.pwm_mapping_policy,
                AttractPwmMappingPolicy::WindowedSubmatrix
            );
            assert!(settings.compare_alternate_policy);
            assert!(!settings.transcript_strand_only);
            assert!(!settings.allow_species_fallback);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_resources_sync_rebase_with_local_fixture() {
    let td = tempdir().expect("tempdir");
    let output_path = td.path().join("rebase.sync.json");
    let mut engine = GentleEngine::from_state(ProjectState::default());

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ResourcesSyncRebase {
            input: resource_fixture_path("rebase.edge.withrefm"),
            output: Some(output_path.to_string_lossy().to_string()),
            commercial_only: true,
        },
    )
    .expect("execute resources sync-rebase");

    assert!(!out.state_changed);
    assert_eq!(
        out.output["report"]["resource"].as_str(),
        Some("rebase-commercial")
    );
    assert_eq!(out.output["report"]["item_count"].as_u64(), Some(2));
    assert_eq!(
        out.output["report"]["output"].as_str(),
        Some(output_path.to_string_lossy().as_ref())
    );

    let written = fs::read_to_string(&output_path).expect("read synced REBASE output");
    let enzymes = serde_json::from_str::<serde_json::Value>(&written).expect("parse JSON");
    let names = enzymes
        .as_array()
        .expect("enzyme array")
        .iter()
        .filter_map(|entry| entry.get("name").and_then(|v| v.as_str()))
        .collect::<Vec<_>>();
    assert_eq!(names, vec!["BsaI", "EcoRI"]);
}

#[test]
fn execute_resources_sync_jaspar_with_local_fixture() {
    let _serial = lock_jaspar_tests();
    crate::tf_motifs::reload();
    let _guard = JasparReloadResetGuard;
    let td = tempdir().expect("tempdir");
    let output_path = td.path().join("jaspar.sync.json");
    let mut engine = GentleEngine::from_state(ProjectState::default());

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ResourcesSyncJaspar {
            input: resource_fixture_path("jaspar.edge.pfm"),
            output: Some(output_path.to_string_lossy().to_string()),
        },
    )
    .expect("execute resources sync-jaspar");

    assert!(!out.state_changed);
    assert_eq!(out.output["report"]["resource"].as_str(), Some("jaspar"));
    assert_eq!(out.output["report"]["item_count"].as_u64(), Some(2));
    assert_eq!(
        out.output["report"]["output"].as_str(),
        Some(output_path.to_string_lossy().as_ref())
    );

    let written = fs::read_to_string(&output_path).expect("read synced JASPAR output");
    let snapshot = serde_json::from_str::<serde_json::Value>(&written).expect("parse JSON");
    assert_eq!(
        snapshot.get("schema").and_then(|v| v.as_str()),
        Some("gentle.tf_motifs.v2")
    );
    assert_eq!(
        snapshot.get("motif_count").and_then(|v| v.as_u64()),
        Some(2)
    );
}

#[test]
fn execute_resources_summarize_jaspar_writes_report_and_returns_payload() {
    let _serial = lock_jaspar_tests();
    crate::tf_motifs::reload();
    let td = tempdir().expect("tempdir");
    let output_path = td.path().join("jaspar.summary.json");
    let mut engine = GentleEngine::from_state(ProjectState::default());

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ResourcesSummarizeJaspar {
            motifs: vec!["SP1".to_string()],
            random_sequence_length_bp: 512,
            random_seed: 42,
            output: Some(output_path.to_string_lossy().to_string()),
        },
    )
    .expect("execute resources summarize-jaspar");

    assert!(!out.state_changed);
    assert_eq!(
        out.output["result"]["jaspar_entry_presentation"]["schema"].as_str(),
        Some("gentle.jaspar_entry_presentation.v1")
    );
    assert_eq!(
        out.output["result"]["jaspar_entry_presentation"]["resolved_entry_count"].as_u64(),
        Some(1)
    );
    assert_eq!(
        out.output["result"]["jaspar_entry_presentation"]["rows"][0]["maximizing_sequence"]
            .as_str(),
        Some("GGGGCGGGG")
    );
    let written = fs::read_to_string(&output_path).expect("read JASPAR summary JSON");
    let json: serde_json::Value = serde_json::from_str(&written).expect("parse JSON");
    assert_eq!(
        json.get("random_sequence_length_bp")
            .and_then(|value| value.as_u64()),
        Some(512)
    );
}

#[test]
fn execute_resources_benchmark_jaspar_writes_report_and_returns_payload() {
    let _serial = lock_jaspar_tests();
    crate::tf_motifs::reload();
    let td = tempdir().expect("tempdir");
    let output_path = td.path().join("jaspar.benchmark.json");
    let mut engine = GentleEngine::from_state(ProjectState::default());

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ResourcesBenchmarkJaspar {
            random_sequence_length_bp: 512,
            random_seed: 42,
            output: Some(output_path.to_string_lossy().to_string()),
        },
    )
    .expect("execute resources benchmark-jaspar");

    assert!(!out.state_changed);
    assert_eq!(
        out.output["result"]["jaspar_registry_benchmark"]["schema"].as_str(),
        Some("gentle.jaspar_registry_benchmark.v1")
    );
    assert_eq!(
        out.output["result"]["jaspar_registry_benchmark"]["score_family_summaries"]
            .as_array()
            .map(|rows| rows.len()),
        Some(2)
    );
    let written = fs::read_to_string(&output_path).expect("read JASPAR benchmark JSON");
    let json: serde_json::Value = serde_json::from_str(&written).expect("parse JSON");
    assert_eq!(
        json.get("random_sequence_length_bp")
            .and_then(|value| value.as_u64()),
        Some(512)
    );
}

#[test]
fn execute_resources_sync_attract_with_local_fixture() {
    let _serial = attract_test_lock().lock().expect("attract test lock");
    struct ReloadResetGuard;
    impl Drop for ReloadResetGuard {
        fn drop(&mut self) {
            crate::attract_motifs::reload();
        }
    }
    let _guard = ReloadResetGuard;

    let td = tempdir().expect("tempdir");
    let archive_path = write_synthetic_attract_zip(td.path());
    let output_path = td.path().join("attract.sync.json");
    let mut engine = GentleEngine::from_state(ProjectState::default());

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ResourcesSyncAttract {
            input: archive_path.to_string_lossy().to_string(),
            output: Some(output_path.to_string_lossy().to_string()),
        },
    )
    .expect("execute resources sync-attract");

    assert!(!out.state_changed);
    assert_eq!(out.output["report"]["resource"].as_str(), Some("attract"));
    assert_eq!(out.output["report"]["item_count"].as_u64(), Some(2));
    assert_eq!(
        out.output["report"]["output"].as_str(),
        Some(output_path.to_string_lossy().as_ref())
    );

    let written = fs::read_to_string(&output_path).expect("read synced ATtRACT output");
    let snapshot = serde_json::from_str::<serde_json::Value>(&written).expect("parse JSON");
    assert_eq!(
        snapshot.get("schema").and_then(|v| v.as_str()),
        Some("gentle.attract_motifs.v1")
    );
    assert_eq!(
        snapshot.get("motif_count").and_then(|v| v.as_u64()),
        Some(2)
    );
    assert_eq!(crate::attract_motifs::list_motif_summaries().len(), 2);
}

#[test]
fn execute_resources_inspect_jaspar_writes_report_and_returns_payload() {
    let _serial = lock_jaspar_tests();
    crate::tf_motifs::reload();
    let td = tempdir().expect("tempdir");
    let output_path = td.path().join("jaspar.expert.json");
    let mut engine = GentleEngine::from_state(ProjectState::default());

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ResourcesInspectJaspar {
            motif: "SP1".to_string(),
            random_sequence_length_bp: 512,
            random_seed: 42,
            fetch_remote: false,
            output: Some(output_path.to_string_lossy().to_string()),
        },
    )
    .expect("execute resources inspect-jaspar");

    assert!(!out.state_changed);
    assert_eq!(
        out.output["result"]["jaspar_entry_expert_view"]["schema"].as_str(),
        Some("gentle.jaspar_entry_expert.v1")
    );
    assert_eq!(
        out.output["result"]["jaspar_entry_expert_view"]["motif_name"].as_str(),
        Some("SP1")
    );
    let written = fs::read_to_string(&output_path).expect("read JASPAR expert JSON");
    let json: serde_json::Value = serde_json::from_str(&written).expect("parse JSON");
    assert_eq!(
        json.get("score_panels")
            .and_then(|value| value.as_array())
            .map(|rows| rows.len()),
        Some(2)
    );
}

#[test]
fn execute_resources_list_jaspar_writes_report_and_returns_payload() {
    let _serial = lock_jaspar_tests();
    crate::tf_motifs::reload();
    let td = tempdir().expect("tempdir");
    let output_path = td.path().join("jaspar.catalog.json");
    let mut engine = GentleEngine::from_state(ProjectState::default());

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ResourcesListJaspar {
            filter: Some("SP1".to_string()),
            limit: Some(10),
            fetch_remote: false,
            output: Some(output_path.to_string_lossy().to_string()),
        },
    )
    .expect("execute resources list-jaspar");

    assert!(!out.state_changed);
    assert_eq!(
        out.output["result"]["jaspar_catalog_report"]["schema"].as_str(),
        Some("gentle.jaspar_catalog.v1")
    );
    let returned = out.output["result"]["jaspar_catalog_report"]["returned_entry_count"]
        .as_u64()
        .expect("returned count");
    assert!(returned >= 1);
    assert!(returned <= 10);
    let written = fs::read_to_string(&output_path).expect("read JASPAR catalog JSON");
    let json: serde_json::Value = serde_json::from_str(&written).expect("parse JSON");
    assert_eq!(
        json.get("filter").and_then(|value| value.as_str()),
        Some("SP1")
    );
}

#[test]
fn execute_attract_inspect_splicing_uses_shared_engine_view() {
    let _serial = attract_test_lock().lock().expect("attract test lock");
    struct ReloadResetGuard;
    impl Drop for ReloadResetGuard {
        fn drop(&mut self) {
            crate::attract_motifs::reload();
        }
    }
    let _guard = ReloadResetGuard;
    let td = tempdir().expect("tempdir");
    let archive_path = write_synthetic_attract_zip(td.path());
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence("TTTTGAAGAACCCTCTTGGGAAAAAAAAAA").expect("dna");
    dna.features_mut().push(Feature {
        kind: "SOURCE".into(),
        location: Location::simple_range(0, 30),
        qualifiers: vec![("organism".into(), Some("Homo sapiens".to_string()))],
    });
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::Join(vec![
            Location::simple_range(0, 10),
            Location::simple_range(20, 30),
        ]),
        qualifiers: vec![
            ("gene".into(), Some("TP73".to_string())),
            ("transcript_id".into(), Some("NM_TP73_1".to_string())),
            ("label".into(), Some("NM_TP73_1".to_string())),
        ],
    });
    state.sequences.insert("tp73".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);
    execute_shell_command(
        &mut engine,
        &ShellCommand::ResourcesSyncAttract {
            input: archive_path.to_string_lossy().to_string(),
            output: Some(
                td.path()
                    .join("attract.sync.json")
                    .to_string_lossy()
                    .to_string(),
            ),
        },
    )
    .expect("sync attract");
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::InspectSplicingAttract {
            seq_id: "tp73".to_string(),
            feature_id: 1,
            settings: AttractSplicingEvidenceSettings {
                scope: SplicingScopePreset::TargetGroupTargetStrand,
                transcript_strand_only: true,
                boundary_flank_bp: 2,
                requested_organism: None,
                allow_species_fallback: true,
                minimum_quality_score: 0.0,
                minimum_match_quantile: 0.99,
                pwm_mapping_policy: AttractPwmMappingPolicy::StrictSameLength,
                compare_alternate_policy: true,
            },
        },
    )
    .expect("inspect attract splicing");
    assert_eq!(
        out.output["schema"].as_str(),
        Some("gentle.attract_splicing_evidence.v1")
    );
    assert!(
        out.output["active_resource_fingerprint"]
            .as_str()
            .map(|value| value.starts_with("sha1:"))
            .unwrap_or(false)
    );
    assert_eq!(
        out.output["alternate_policy_summary"]["pwm_mapping_policy"].as_str(),
        Some("windowed_submatrix")
    );
    assert_eq!(out.output["hit_count"].as_u64(), Some(2));
    assert_eq!(
        out.output["summary_rows"][0]["gene_name"].as_str(),
        Some("PTBP1")
    );
    let hit_rows = out.output["hit_rows"].as_array().expect("hit rows array");
    assert!(hit_rows.iter().any(|row| {
        row.get("match_score_kind").and_then(|value| value.as_str()) == Some("llr_bits")
    }));
    assert!(hit_rows.iter().any(|row| {
        row.get("match_score_kind").and_then(|value| value.as_str()) == Some("exact_match_bp")
    }));
}

#[test]
fn execute_resources_status_reports_builtin_or_runtime_sources() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(&mut engine, &ShellCommand::ResourcesStatus)
        .expect("execute resources status");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["schema"].as_str(),
        Some("gentle.resource_status.v1")
    );
    assert_eq!(out.output["rebase"]["resource_id"].as_str(), Some("rebase"));
    assert_eq!(out.output["jaspar"]["resource_id"].as_str(), Some("jaspar"));
    assert_eq!(
        out.output["rebase"]["active_source"].as_str(),
        Some("builtin")
    );
    assert_eq!(
        out.output["jaspar"]["active_source"].as_str(),
        Some("builtin")
    );
    assert_eq!(
        out.output["jaspar"]["remote_metadata_snapshot"]["path"].as_str(),
        Some("data/resources/jaspar.remote_metadata.json")
    );
    let jaspar_remote_metadata_exists =
        std::path::Path::new("data/resources/jaspar.remote_metadata.json").exists();
    assert_eq!(
        out.output["jaspar"]["remote_metadata_snapshot"]["exists"].as_bool(),
        Some(jaspar_remote_metadata_exists)
    );
    assert_eq!(
        out.output["jaspar"]["remote_metadata_snapshot"]["valid"].as_bool(),
        Some(jaspar_remote_metadata_exists)
    );
    assert_eq!(
        out.output["attract"]["support_status"].as_str(),
        Some(
            if std::path::Path::new(crate::attract_motifs::DEFAULT_ATTRACT_RESOURCE_PATH).exists() {
                "runtime_snapshot"
            } else {
                "known_external_only"
            },
        )
    );
    assert_eq!(
        out.output["attract"]["display_name"].as_str(),
        Some("ATtRACT")
    );
    if out.output["attract"]["support_status"].as_str() == Some("runtime_snapshot") {
        assert!(
            out.output["attract"]["active_pwm_row_count"]
                .as_u64()
                .unwrap_or(0)
                > 0
                || out.output["attract"]["active_consensus_only_row_count"]
                    .as_u64()
                    .unwrap_or(0)
                    > 0
        );
        assert_eq!(
            out.output["attract"]["active_source"].as_str(),
            Some("runtime")
        );
    } else {
        assert_eq!(
            out.output["attract"]["active_pwm_row_count"].as_u64(),
            Some(0)
        );
        assert_eq!(
            out.output["attract"]["active_source"].as_str(),
            Some("unavailable")
        );
    }
    assert_eq!(
        out.output["attract"]["download_url"].as_str(),
        Some("https://attract.cnic.es/attract/static/ATtRACT.zip")
    );
}

#[test]
fn execute_proteases_list_writes_catalog_json() {
    let td = tempdir().expect("tempdir");
    let output_path = td.path().join("proteases.catalog.json");
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ProteasesList {
            filter: Some("tag_removal".to_string()),
            output: Some(output_path.to_string_lossy().to_string()),
        },
    )
    .expect("execute proteases list");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["schema"].as_str(),
        Some("gentle.protease_catalog.v1")
    );
    assert!(out.output["proteases"].as_array().is_some_and(|rows| {
        rows.iter()
            .any(|row| row["name"].as_str() == Some("TEV protease"))
    }));
    let written = fs::read_to_string(&output_path).expect("read protease catalog");
    let json: serde_json::Value = serde_json::from_str(&written).expect("parse JSON");
    assert_eq!(
        json.get("filter").and_then(|value| value.as_str()),
        Some("tag_removal")
    );
}

#[test]
fn execute_proteases_show_returns_trypsin_metadata() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ProteasesShow {
            query: "Trypsin".to_string(),
            output: None,
        },
    )
    .expect("execute proteases show");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["schema"].as_str(),
        Some("gentle.protease_catalog_entry.v1")
    );
    assert_eq!(out.output["name"].as_str(), Some("Trypsin"));
    assert_eq!(
        out.output["specificity"].as_str(),
        Some("Cleaves C-terminal to Lys or Arg unless the next residue is Pro.")
    );
}

#[test]
fn execute_services_status_reports_combined_readiness() {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(&mut engine, &ShellCommand::ServicesStatus)
        .expect("execute services status");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["schema"].as_str(),
        Some("gentle.service_readiness.v1")
    );
    assert_eq!(
        out.output["references"][0]["genome_id"].as_str(),
        Some("Human GRCh38 Ensembl 116")
    );
    assert!(
        out.output["references"][0]["resource_key"]
            .as_str()
            .is_some()
    );
    assert!(
        out.output["references"][0]["lifecycle_status"]
            .as_str()
            .is_some()
    );
    assert_eq!(
        out.output["helpers"][0]["genome_id"].as_str(),
        Some("Plasmid pUC19 (online)")
    );
    assert_eq!(
        out.output["resources"]["attract"]["display_name"].as_str(),
        Some("ATtRACT")
    );
    assert!(
        out.output["summary_lines"]
            .as_array()
            .map(|rows| !rows.is_empty())
            .unwrap_or(false)
    );
}

#[test]
fn execute_services_handoff_reports_actions_and_writes_json() {
    let td = tempdir().expect("tempdir");
    let output_path = td.path().join("service_handoff.json");
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ServicesHandoff {
            scope: Some("clawbio".to_string()),
            output: Some(output_path.to_string_lossy().to_string()),
        },
    )
    .expect("execute services handoff");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["schema"].as_str(),
        Some("gentle.service_handoff.v1")
    );
    assert_eq!(out.output["scope"].as_str(), Some("clawbio"));
    assert_eq!(
        out.output["service_readiness"]["schema"].as_str(),
        Some("gentle.service_readiness.v1")
    );
    assert!(
        out.output["readiness"]
            .as_array()
            .map(|rows| rows.iter().any(|row| row["resource_key"].as_str()
                == Some("reference_genome:Human GRCh38 Ensembl 116")))
            .unwrap_or(false)
    );
    assert!(out.output["suggested_actions"].as_array().is_some());
    assert!(out.output["blocked_actions"].as_array().is_some());
    assert!(
        out.output["preferred_demo_actions"]
            .as_array()
            .map(|rows| rows
                .iter()
                .any(|row| row["kind"].as_str() == Some("demo_graphic")))
            .unwrap_or(false)
    );
    assert!(output_path.exists());
    let persisted: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(output_path).expect("read handoff report"))
            .expect("parse persisted handoff report");
    assert_eq!(
        persisted["preferred_artifacts"][0]["recommended_use"].as_str(),
        Some("installation_doctor_payload")
    );
}

#[test]
fn execute_resources_sync_jaspar_reloads_motif_registry_from_synced_output() {
    let _serial = lock_jaspar_tests();
    crate::tf_motifs::reload();
    let _guard = JasparReloadResetGuard;

    let td = tempdir().expect("tempdir");
    let unique = JASPAR_RELOAD_TEST_COUNTER.fetch_add(1, Ordering::Relaxed);
    let motif_id = format!("MTEST{unique}.1");
    let motif_name = format!("CODEX_TEST_MOTIF_{unique}");
    let input_path = td.path().join("custom_reload.pfm");
    let output_path = td.path().join("custom_reload.motifs.json");
    let motif_text =
        format!(">{motif_id} {motif_name}\nA [1 0 0 0]\nC [0 1 0 0]\nG [0 0 1 0]\nT [0 0 0 1]\n");
    fs::write(&input_path, motif_text).expect("write custom JASPAR input");

    assert_eq!(crate::tf_motifs::resolve_motif(&motif_name), None);

    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ResourcesSyncJaspar {
            input: input_path.to_string_lossy().to_string(),
            output: Some(output_path.to_string_lossy().to_string()),
        },
    )
    .expect("execute resources sync-jaspar");

    assert!(!out.state_changed);
    assert_eq!(out.output["report"]["item_count"].as_u64(), Some(1));
    assert_eq!(
        crate::tf_motifs::resolve_motif(&motif_name).as_deref(),
        Some("ACGT")
    );
}

#[test]
fn parse_import_pool_with_prefix() {
    let cmd = parse_shell_line("import-pool test_files/demo.pool.gentle.json imported")
        .expect("parse import-pool");
    match cmd {
        ShellCommand::ImportPool { input, prefix } => {
            assert_eq!(input, "test_files/demo.pool.gentle.json".to_string());
            assert_eq!(prefix, "imported".to_string());
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_import_pool_loads_members_into_engine_state() {
    let td = tempdir().expect("tempdir");
    let pool_path = write_demo_pool_json(td.path());

    let mut state = ProjectState::default();
    state.sequences.insert(
        "imported_1".to_string(),
        DNAsequence::from_sequence("AAAA").expect("seed sequence"),
    );
    let mut engine = GentleEngine::from_state(state);
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ImportPool {
            input: pool_path.to_string_lossy().to_string(),
            prefix: "imported".to_string(),
        },
    )
    .expect("execute import-pool");

    assert!(out.state_changed);
    assert_eq!(out.output["pool_id"].as_str(), Some("demo_pool"));
    assert_eq!(out.output["member_count"].as_u64(), Some(1));
    let imported_ids = out
        .output
        .get("imported_ids")
        .and_then(|v| v.as_array())
        .expect("imported_ids array");
    assert_eq!(imported_ids.len(), 1);
    let imported_id = imported_ids[0].as_str().expect("imported id string");
    assert!(
        imported_id.starts_with("imported_1"),
        "import id should be a collision-resolved variant of imported_1, got {imported_id}"
    );
    assert_ne!(imported_id, "imported_1");
    assert!(engine.state().sequences.contains_key(imported_id));
}

#[test]
fn execute_import_pool_preserves_supercoiled_topology_hint() {
    let td = tempdir().expect("tempdir");
    let pool_path = td.path().join("supercoiled.pool.gentle.json");
    let pool_json = serde_json::json!({
        "schema": "gentle.pool.v1",
        "pool_id": "supercoiled_pool",
        "human_id": "demo",
        "member_count": 1,
        "members": [
            {
                "seq_id": "member_1",
                "human_id": "member_1",
                "name": "Vector supercoiled",
                "sequence": "ATGCATGC",
                "length_bp": 8,
                "topology": "supercoiled",
                "ends": {
                    "end_type": "blunt",
                    "forward_5": "",
                    "forward_3": "",
                    "reverse_5": "",
                    "reverse_3": ""
                }
            }
        ]
    });
    fs::write(
        &pool_path,
        serde_json::to_string_pretty(&pool_json).unwrap(),
    )
    .unwrap();

    let mut engine = GentleEngine::from_state(ProjectState::default());
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ImportPool {
            input: pool_path.to_string_lossy().to_string(),
            prefix: "imported".to_string(),
        },
    )
    .expect("execute import-pool");

    let imported_id = out.output["imported_ids"][0].as_str().expect("imported id");
    let imported = engine
        .state()
        .sequences
        .get(imported_id)
        .expect("imported sequence");
    assert!(imported.is_circular());
    assert!(
        imported
            .description()
            .iter()
            .any(|line| line.contains("gel_topology=supercoiled"))
    );
}

#[test]
fn parse_helpers_extract_gene() {
    let cmd = parse_shell_line(
        "helpers extract-gene pUC19 bla --occurrence 2 --output-id out --cache-dir cache",
    )
    .expect("parse helpers extract-gene");
    match cmd {
        ShellCommand::ReferenceExtractGene {
            helper_mode,
            genome_id,
            gene_query,
            occurrence,
            output_id,
            cache_dir,
            ..
        } => {
            assert!(helper_mode);
            assert_eq!(genome_id, "pUC19".to_string());
            assert_eq!(gene_query, "bla".to_string());
            assert_eq!(occurrence, Some(2));
            assert_eq!(output_id, Some("out".to_string()));
            assert_eq!(cache_dir, Some("cache".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_extract_gene_with_scope_and_cap() {
    let cmd = parse_shell_line(
            "genomes extract-gene ToyGenome MYGENE --annotation-scope full --max-annotation-features 120 --output-id out",
        )
        .expect("parse genomes extract-gene with scope and cap");
    match cmd {
        ShellCommand::ReferenceExtractGene {
            helper_mode,
            genome_id,
            gene_query,
            output_id,
            annotation_scope,
            max_annotation_features,
            include_genomic_annotation,
            ..
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome".to_string());
            assert_eq!(gene_query, "MYGENE".to_string());
            assert_eq!(output_id, Some("out".to_string()));
            assert_eq!(annotation_scope, Some(GenomeAnnotationScope::Full));
            assert_eq!(max_annotation_features, Some(120));
            assert_eq!(include_genomic_annotation, None);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_extract_gene_with_coding_promoter_mode() {
    let cmd = parse_shell_line(
        "genomes extract-gene ToyGenome MYGENE --extract-mode coding_with_promoter --promoter-upstream-bp 2000 --output-id out",
    )
    .expect("parse genomes extract-gene with coding promoter mode");
    match cmd {
        ShellCommand::ReferenceExtractGene {
            helper_mode,
            genome_id,
            gene_query,
            output_id,
            extract_mode,
            promoter_upstream_bp,
            ..
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome".to_string());
            assert_eq!(gene_query, "MYGENE".to_string());
            assert_eq!(output_id, Some("out".to_string()));
            assert_eq!(
                extract_mode,
                Some(crate::engine::GenomeGeneExtractMode::CodingWithPromoter)
            );
            assert_eq!(promoter_upstream_bp, Some(2000));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_extract_promoter_with_scope_and_transcript() {
    let cmd = parse_shell_line(
        "genomes extract-promoter ToyGenome MYGENE --occurrence 2 --transcript-id TX2 --output-id promoter_out --upstream-bp 750 --downstream-bp 80 --annotation-scope full --max-annotation-features 25 --catalog c.json --cache-dir cache",
    )
    .expect("parse genomes extract-promoter");
    match cmd {
        ShellCommand::ReferenceExtractPromoter {
            helper_mode,
            genome_id,
            gene_query,
            occurrence,
            transcript_id,
            output_id,
            upstream_bp,
            downstream_bp,
            annotation_scope,
            max_annotation_features,
            include_genomic_annotation,
            catalog_path,
            cache_dir,
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome".to_string());
            assert_eq!(gene_query, "MYGENE".to_string());
            assert_eq!(occurrence, Some(2));
            assert_eq!(transcript_id, Some("TX2".to_string()));
            assert_eq!(output_id, Some("promoter_out".to_string()));
            assert_eq!(upstream_bp, Some(750));
            assert_eq!(downstream_bp, Some(80));
            assert_eq!(annotation_scope, Some(GenomeAnnotationScope::Full));
            assert_eq!(max_annotation_features, Some(25));
            assert_eq!(include_genomic_annotation, None);
            assert_eq!(catalog_path, Some("c.json".to_string()));
            assert_eq!(cache_dir, Some("cache".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_promoter_tfbs_summary_with_multiple_gene_specs() {
    let cmd = parse_shell_line(
        "genomes promoter-tfbs-summary ToyGenome --gene POS1 --gene NEG1::2@TX_NEG#NEG1_panel --motif stemness --motifs SP1,TP73 --upstream-bp 750 --downstream-bp 80 --score-kind llr_background_tail_log10 --allow-negative --catalog c.json --cache-dir cache --path out.json",
    )
    .expect("parse genomes promoter-tfbs-summary");
    match cmd {
        ShellCommand::ReferencePromoterTfbsSummary {
            helper_mode,
            genome_id,
            genes,
            motifs,
            upstream_bp,
            downstream_bp,
            score_kind,
            clip_negative,
            catalog_path,
            cache_dir,
            path,
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome".to_string());
            assert_eq!(
                genes,
                vec![
                    PromoterTfbsGeneQuery {
                        gene_query: "POS1".to_string(),
                        occurrence: None,
                        transcript_id: None,
                        display_label: None,
                    },
                    PromoterTfbsGeneQuery {
                        gene_query: "NEG1".to_string(),
                        occurrence: Some(2),
                        transcript_id: Some("TX_NEG".to_string()),
                        display_label: Some("NEG1_panel".to_string()),
                    },
                ]
            );
            assert_eq!(
                motifs,
                vec![
                    "stemness".to_string(),
                    "SP1".to_string(),
                    "TP73".to_string()
                ]
            );
            assert_eq!(upstream_bp, 750);
            assert_eq!(downstream_bp, 80);
            assert_eq!(score_kind, TfbsScoreTrackValueKind::LlrBackgroundTailLog10);
            assert!(!clip_negative);
            assert_eq!(catalog_path, Some("c.json".to_string()));
            assert_eq!(cache_dir, Some("cache".to_string()));
            assert_eq!(path, Some("out.json".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_extract_gene_with_annotation_flag() {
    let cmd = parse_shell_line(
        "genomes extract-gene ToyGenome MYGENE --include-genomic-annotation --output-id out",
    )
    .expect("parse genomes extract-gene with annotation flag");
    match cmd {
        ShellCommand::ReferenceExtractGene {
            helper_mode,
            genome_id,
            gene_query,
            output_id,
            annotation_scope,
            include_genomic_annotation,
            ..
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome".to_string());
            assert_eq!(gene_query, "MYGENE".to_string());
            assert_eq!(output_id, Some("out".to_string()));
            assert_eq!(annotation_scope, Some(GenomeAnnotationScope::Core));
            assert_eq!(include_genomic_annotation, Some(true));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_extract_region_with_annotation_flag() {
    let cmd = parse_shell_line(
            "genomes extract-region ToyGenome chr1 100 250 --include-genomic-annotation --output-id out",
        )
        .expect("parse genomes extract-region with annotation flag");
    match cmd {
        ShellCommand::ReferenceExtractRegion {
            helper_mode,
            genome_id,
            chromosome,
            start_1based,
            end_1based,
            output_id,
            annotation_scope,
            include_genomic_annotation,
            ..
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome".to_string());
            assert_eq!(chromosome, "chr1".to_string());
            assert_eq!(start_1based, 100);
            assert_eq!(end_1based, 250);
            assert_eq!(output_id, Some("out".to_string()));
            assert_eq!(annotation_scope, Some(GenomeAnnotationScope::Core));
            assert_eq!(include_genomic_annotation, Some(true));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_extract_region_with_scope_and_cap() {
    let cmd = parse_shell_line(
            "genomes extract-region ToyGenome chr1 100 250 --annotation-scope full --max-annotation-features 1200 --output-id out",
        )
        .expect("parse genomes extract-region with scope and cap");
    match cmd {
        ShellCommand::ReferenceExtractRegion {
            helper_mode,
            genome_id,
            chromosome,
            start_1based,
            end_1based,
            output_id,
            annotation_scope,
            max_annotation_features,
            include_genomic_annotation,
            ..
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome".to_string());
            assert_eq!(chromosome, "chr1".to_string());
            assert_eq!(start_1based, 100);
            assert_eq!(end_1based, 250);
            assert_eq!(output_id, Some("out".to_string()));
            assert_eq!(annotation_scope, Some(GenomeAnnotationScope::Full));
            assert_eq!(max_annotation_features, Some(1200));
            assert_eq!(include_genomic_annotation, None);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_prepare_with_timeout() {
    let cmd = parse_shell_line(
        "genomes prepare ToyGenome --catalog c.json --cache-dir cache --timeout-secs 90",
    )
    .expect("parse genomes prepare");
    match cmd {
        ShellCommand::ReferencePrepare {
            helper_mode,
            genome_id,
            catalog_path,
            cache_dir,
            timeout_seconds,
        } => {
            assert!(!helper_mode);
            assert_eq!(genome_id, "ToyGenome".to_string());
            assert_eq!(catalog_path, Some("c.json".to_string()));
            assert_eq!(cache_dir, Some("cache".to_string()));
            assert_eq!(timeout_seconds, Some(90));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_extend_anchor_five_prime() {
    let cmd = parse_shell_line(
        "genomes extend-anchor tp73 5p 150 --output-id tp73_ext --catalog c.json --cache-dir cache",
    )
    .expect("parse genomes extend-anchor");
    match cmd {
        ShellCommand::ReferenceExtendAnchor {
            helper_mode,
            seq_id,
            side,
            length_bp,
            output_id,
            catalog_path,
            cache_dir,
            prepared_genome_id,
        } => {
            assert!(!helper_mode);
            assert_eq!(seq_id, "tp73".to_string());
            assert_eq!(side, GenomeAnchorSide::FivePrime);
            assert_eq!(length_bp, 150);
            assert_eq!(output_id, Some("tp73_ext".to_string()));
            assert_eq!(catalog_path, Some("c.json".to_string()));
            assert_eq!(cache_dir, Some("cache".to_string()));
            assert_eq!(prepared_genome_id, None);
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_extend_anchor_with_prepared_genome() {
    let cmd = parse_shell_line(
        "genomes extend-anchor tp73 3p 200 --prepared-genome Human_GRCh38 --catalog c.json",
    )
    .expect("parse genomes extend-anchor with prepared genome");
    match cmd {
        ShellCommand::ReferenceExtendAnchor {
            helper_mode,
            seq_id,
            side,
            length_bp,
            output_id,
            catalog_path,
            cache_dir,
            prepared_genome_id,
        } => {
            assert!(!helper_mode);
            assert_eq!(seq_id, "tp73".to_string());
            assert_eq!(side, GenomeAnchorSide::ThreePrime);
            assert_eq!(length_bp, 200);
            assert_eq!(output_id, None);
            assert_eq!(catalog_path, Some("c.json".to_string()));
            assert_eq!(cache_dir, None);
            assert_eq!(prepared_genome_id, Some("Human_GRCh38".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_verify_anchor_with_prepared_genome() {
    let cmd = parse_shell_line(
            "genomes verify-anchor tp73 --catalog c.json --cache-dir cache --prepared-genome Human_GRCh38",
        )
        .expect("parse genomes verify-anchor");
    match cmd {
        ShellCommand::ReferenceVerifyAnchor {
            helper_mode,
            seq_id,
            catalog_path,
            cache_dir,
            prepared_genome_id,
        } => {
            assert!(!helper_mode);
            assert_eq!(seq_id, "tp73".to_string());
            assert_eq!(catalog_path, Some("c.json".to_string()));
            assert_eq!(cache_dir, Some("cache".to_string()));
            assert_eq!(prepared_genome_id, Some("Human_GRCh38".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_genomes_extend_anchor_creates_sequence() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
    let cache_dir = td.path().join("cache");
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
    fs::write(&catalog, catalog_json).expect("write catalog");
    let catalog_path = catalog.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    let prepare = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferencePrepare {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.clone()),
            cache_dir: None,
            timeout_seconds: None,
        },
    )
    .expect("prepare genome");
    assert_eq!(
        prepare.output["binary_preflight"]["schema"].as_str(),
        Some("gentle.blast_external_binary_preflight.v1")
    );

    execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceExtractRegion {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 3,
            end_1based: 10,
            output_id: Some("slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path.clone()),
            cache_dir: None,
        },
    )
    .expect("extract region");

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceExtendAnchor {
            helper_mode: false,
            seq_id: "slice".to_string(),
            side: GenomeAnchorSide::FivePrime,
            length_bp: 2,
            output_id: Some("slice_ext5".to_string()),
            catalog_path: Some(catalog_path),
            cache_dir: None,
            prepared_genome_id: None,
        },
    )
    .expect("execute extend-anchor");

    assert!(out.state_changed);
    let created = out
        .output
        .get("result")
        .and_then(|v| v.get("created_seq_ids"))
        .and_then(|v| v.as_array())
        .expect("created_seq_ids");
    assert!(
        created
            .iter()
            .any(|v| v.as_str().map(|id| id == "slice_ext5").unwrap_or(false))
    );
    let seq = engine
        .state()
        .sequences
        .get("slice_ext5")
        .expect("extended sequence in state");
    assert_eq!(seq.get_forward_string(), "ACGTACGTAC");
}

#[test]
fn execute_genomes_extract_region_default_scope_core_with_telemetry() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGTACGT\n").expect("write fasta");
    fs::write(
            &gtf,
            concat!(
                "chr1\tsrc\tgene\t1\t16\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
                "chr1\tsrc\ttranscript\t1\t16\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\"; transcript_id \"TX1\";\n",
                "chr1\tsrc\texon\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\"; transcript_id \"TX1\"; exon_number \"1\";\n",
            ),
        )
        .expect("write gtf");
    let catalog = td.path().join("catalog.json");
    let cache_dir = td.path().join("cache");
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
    fs::write(&catalog, catalog_json).expect("write catalog");
    let catalog_path = catalog.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    execute_shell_command(
        &mut engine,
        &ShellCommand::ReferencePrepare {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.clone()),
            cache_dir: None,
            timeout_seconds: None,
        },
    )
    .expect("prepare genome");

    let command = parse_shell_line(&format!(
        "genomes extract-region ToyGenome chr1 1 16 --output-id shell_slice_default --catalog {}",
        catalog_path
    ))
    .expect("parse extract-region");
    let out = execute_shell_command(&mut engine, &command).expect("execute extract-region");
    assert!(out.state_changed);
    let telemetry = out
        .output
        .get("result")
        .and_then(|v| v.get("genome_annotation_projection"))
        .and_then(|v| v.as_object())
        .expect("annotation projection telemetry");
    assert_eq!(
        telemetry.get("requested_scope").and_then(|v| v.as_str()),
        Some("core")
    );
    assert_eq!(
        telemetry.get("effective_scope").and_then(|v| v.as_str()),
        Some("core")
    );
    assert!(
        telemetry
            .get("attached_feature_count")
            .and_then(|v| v.as_u64())
            .unwrap_or(0)
            >= 2
    );
    let seq = engine
        .state()
        .sequences
        .get("shell_slice_default")
        .expect("sequence created by extract-region");
    assert!(
        seq.features()
            .iter()
            .any(|f| f.kind.to_string().eq_ignore_ascii_case("gene"))
    );
}

#[test]
fn execute_genomes_extract_promoter_uses_transcript_tss_on_reverse_strand() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    let sequence: String = (0..3000)
        .map(|idx| match ((idx * 37) + (idx / 7)) % 4 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        })
        .collect();
    fs::write(&fasta, format!(">chr1\n{sequence}\n")).expect("write fasta");
    fs::write(
        &gtf,
        concat!(
            "chr1\tsrc\tgene\t1001\t2000\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\";\n",
            "chr1\tsrc\ttranscript\t1201\t1950\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX1\";\n",
            "chr1\tsrc\texon\t1201\t1400\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX1\"; exon_number \"1\";\n",
            "chr1\tsrc\texon\t1801\t1950\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX1\"; exon_number \"2\";\n",
        ),
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
    let cache_dir = td.path().join("cache");
    fs::write(
        &catalog,
        format!(
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
        ),
    )
    .expect("write catalog");
    let catalog_path = catalog.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    execute_shell_command(
        &mut engine,
        &ShellCommand::ReferencePrepare {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.clone()),
            cache_dir: None,
            timeout_seconds: None,
        },
    )
    .expect("prepare genome");

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceExtractPromoter {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            gene_query: "NEG1".to_string(),
            occurrence: None,
            transcript_id: Some("TX1".to_string()),
            output_id: Some("neg1_promoter_shell".to_string()),
            upstream_bp: Some(100),
            downstream_bp: Some(20),
            annotation_scope: Some(GenomeAnnotationScope::None),
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path),
            cache_dir: None,
        },
    )
    .expect("extract promoter");

    assert!(out.state_changed);
    let created = out
        .output
        .get("result")
        .and_then(|value| value.get("created_seq_ids"))
        .and_then(|value| value.as_array())
        .expect("created_seq_ids");
    assert!(
        created
            .iter()
            .any(|value| value.as_str() == Some("neg1_promoter_shell"))
    );

    let seq = engine
        .state()
        .sequences
        .get("neg1_promoter_shell")
        .expect("sequence created by extract-promoter");
    let expected = sequence[(1930 - 1)..2050].to_string();
    assert_eq!(seq.get_forward_string(), expected);
}

#[test]
fn execute_genomes_promoter_tfbs_summary_and_svg_return_multi_gene_payloads() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    let summary_path = td.path().join("promoter_tfbs_summary.json");
    let svg_path = td.path().join("promoter_tfbs.svg");
    let sequence: String = (0..5000)
        .map(|idx| match ((idx * 17) + (idx / 5)) % 4 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        })
        .collect();
    fs::write(&fasta, format!(">chr1\n{sequence}\n")).expect("write fasta");
    fs::write(
        &gtf,
        concat!(
            "chr1\tsrc\tgene\t501\t950\t.\t+\t.\tgene_id \"GENE_POS\"; gene_name \"POS1\";\n",
            "chr1\tsrc\ttranscript\t601\t920\t.\t+\t.\tgene_id \"GENE_POS\"; gene_name \"POS1\"; transcript_id \"TX_POS\";\n",
            "chr1\tsrc\texon\t601\t700\t.\t+\t.\tgene_id \"GENE_POS\"; gene_name \"POS1\"; transcript_id \"TX_POS\"; exon_number \"1\";\n",
            "chr1\tsrc\texon\t801\t920\t.\t+\t.\tgene_id \"GENE_POS\"; gene_name \"POS1\"; transcript_id \"TX_POS\"; exon_number \"2\";\n",
            "chr1\tsrc\tgene\t2001\t2600\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\";\n",
            "chr1\tsrc\ttranscript\t2051\t2500\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\";\n",
            "chr1\tsrc\texon\t2051\t2200\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\"; exon_number \"1\";\n",
            "chr1\tsrc\texon\t2351\t2500\t.\t-\t.\tgene_id \"GENE_NEG\"; gene_name \"NEG1\"; transcript_id \"TX_NEG\"; exon_number \"2\";\n",
        ),
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
    let cache_dir = td.path().join("cache");
    fs::write(
        &catalog,
        format!(
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
        ),
    )
    .expect("write catalog");
    let catalog_path = catalog.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    let _guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    execute_shell_command(
        &mut engine,
        &ShellCommand::ReferencePrepare {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.clone()),
            cache_dir: None,
            timeout_seconds: None,
        },
    )
    .expect("prepare genome");

    let summary_cmd = parse_shell_line(&format!(
        "genomes promoter-tfbs-summary ToyGenome --gene POS1@TX_POS --gene NEG1@TX_NEG#NEG1_panel --motif SP1 --motif \"Yamanaka factors\" --upstream-bp 100 --downstream-bp 20 --score-kind llr_background_tail_log10 --catalog {} --path {}",
        catalog_path,
        summary_path.to_string_lossy(),
    ))
    .expect("parse promoter-tfbs-summary");
    let summary = execute_shell_command(&mut engine, &summary_cmd).expect("run promoter summary");
    assert!(!summary.state_changed);
    assert_eq!(
        summary.output["result"]["multi_gene_promoter_tfbs"]["schema"].as_str(),
        Some("gentle.multi_gene_promoter_tfbs.v1")
    );
    assert_eq!(
        summary.output["result"]["multi_gene_promoter_tfbs"]["returned_gene_count"].as_u64(),
        Some(2)
    );
    assert_eq!(
        summary.output["result"]["multi_gene_promoter_tfbs"]["genes"][1]["sequence_orientation"]
            .as_str(),
        Some("transcription_aligned")
    );
    let summary_json = fs::read_to_string(&summary_path).expect("read summary json");
    assert!(summary_json.contains("\"returned_gene_count\": 2"));

    let svg_cmd = parse_shell_line(&format!(
        "genomes promoter-tfbs-svg ToyGenome --gene POS1@TX_POS --gene NEG1@TX_NEG#NEG1_panel --motif SP1 --motif stemness --upstream-bp 100 --downstream-bp 20 --score-kind llr_background_tail_log10 --catalog {} {}",
        catalog_path,
        svg_path.to_string_lossy(),
    ))
    .expect("parse promoter-tfbs-svg");
    let svg = execute_shell_command(&mut engine, &svg_cmd).expect("run promoter svg");
    assert!(!svg.state_changed);
    assert_eq!(
        svg.output["result"]["multi_gene_promoter_tfbs"]["returned_gene_count"].as_u64(),
        Some(2)
    );
    let svg_text = fs::read_to_string(&svg_path).expect("read promoter svg");
    assert!(svg_text.contains("<svg"));
    assert!(svg_text.contains("POS1"));
    assert!(svg_text.contains("NEG1_panel"));
    assert!(svg_text.contains("TSS"));
}

#[test]
fn execute_genomes_verify_anchor_updates_verification_status() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    fs::write(&fasta, ">chr1\nACGTACGTACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
    let cache_dir = td.path().join("cache");
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
    fs::write(&catalog, catalog_json).expect("write catalog");
    let catalog_path = catalog.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    execute_shell_command(
        &mut engine,
        &ShellCommand::ReferencePrepare {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.clone()),
            cache_dir: None,
            timeout_seconds: None,
        },
    )
    .expect("prepare genome");
    execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceExtractRegion {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            chromosome: "chr1".to_string(),
            start_1based: 3,
            end_1based: 10,
            output_id: Some("slice".to_string()),
            annotation_scope: None,
            max_annotation_features: None,
            include_genomic_annotation: None,
            catalog_path: Some(catalog_path.clone()),
            cache_dir: None,
        },
    )
    .expect("extract region");

    // Mutate sequence to force a mismatch and validate unverified status recording.
    engine.state_mut().sequences.insert(
        "slice".to_string(),
        DNAsequence::from_sequence("AAAAAAAA").expect("mutated sequence"),
    );

    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceVerifyAnchor {
            helper_mode: false,
            seq_id: "slice".to_string(),
            catalog_path: Some(catalog_path),
            cache_dir: None,
            prepared_genome_id: None,
        },
    )
    .expect("verify-anchor");

    assert!(out.state_changed);
    let changed = out
        .output
        .get("result")
        .and_then(|v| v.get("changed_seq_ids"))
        .and_then(|v| v.as_array())
        .expect("changed_seq_ids");
    assert!(
        changed
            .iter()
            .any(|v| v.as_str().map(|id| id == "slice").unwrap_or(false))
    );
    let anchor_status = engine
        .describe_sequence_genome_anchor("slice")
        .expect("anchor status");
    assert!(
        anchor_status.contains("unverified"),
        "expected unverified status, got: {}",
        anchor_status
    );
}

#[test]
fn parse_genomes_validate_catalog() {
    let cmd = parse_shell_line("genomes validate-catalog --catalog assets/genomes.json")
        .expect("parse command");
    match cmd {
        ShellCommand::ReferenceValidateCatalog {
            helper_mode,
            catalog_path,
        } => {
            assert!(!helper_mode);
            assert_eq!(catalog_path, Some("assets/genomes.json".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_helpers_list_with_filter() {
    let cmd = parse_shell_line("helpers list --catalog assets/helper_genomes.json --filter gst")
        .expect("parse command");
    match cmd {
        ShellCommand::ReferenceList {
            helper_mode,
            catalog_path,
            filter,
        } => {
            assert!(helper_mode);
            assert_eq!(catalog_path, Some("assets/helper_genomes.json".to_string()));
            assert_eq!(filter, Some("gst".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_hosts_list_with_filter() {
    let cmd = parse_shell_line("hosts list --catalog assets/host_profiles.json --filter deoR")
        .expect("parse command");
    match cmd {
        ShellCommand::HostsList {
            catalog_path,
            filter,
        } => {
            assert_eq!(catalog_path, Some("assets/host_profiles.json".to_string()));
            assert_eq!(filter, Some("deoR".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_ensembl_available_with_filter() {
    let cmd = parse_shell_line("genomes ensembl-available --collection vertebrates --filter human")
        .expect("parse command");
    match cmd {
        ShellCommand::ReferenceEnsemblAvailable {
            helper_mode,
            collection,
            filter,
        } => {
            assert!(!helper_mode);
            assert_eq!(collection, Some("vertebrates".to_string()));
            assert_eq!(filter, Some("human".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_helpers_install_ensembl_with_overrides() {
    let cmd = parse_shell_line(
        "helpers install-ensembl drosophila_melanogaster --collection metazoa --catalog assets/helper_genomes.json --output-catalog data/helper_overlay.json --genome-id 'Drosophila helper' --cache-dir data/helper_genomes --timeout-secs 90",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::ReferenceInstallEnsembl {
            helper_mode,
            species_dir,
            collection,
            catalog_path,
            output_catalog_path,
            genome_id,
            cache_dir,
            timeout_seconds,
        } => {
            assert!(helper_mode);
            assert_eq!(species_dir, "drosophila_melanogaster");
            assert_eq!(collection, Some("metazoa".to_string()));
            assert_eq!(catalog_path, Some("assets/helper_genomes.json".to_string()));
            assert_eq!(
                output_catalog_path,
                Some("data/helper_overlay.json".to_string())
            );
            assert_eq!(genome_id, Some("Drosophila helper".to_string()));
            assert_eq!(cache_dir, Some("data/helper_genomes".to_string()));
            assert_eq!(timeout_seconds, Some(90));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn helper_operation_catalog_path_defaults_to_discovery_token() {
    assert_eq!(
        operation_catalog_path(&None, true),
        Some(crate::genomes::default_catalog_discovery_token(true).to_string())
    );
    assert_eq!(operation_catalog_path(&None, false), None);
}

#[test]
fn parse_genomes_update_ensembl_specs() {
    let cmd = parse_shell_line(
        "genomes update-ensembl-specs --catalog assets/genomes.json --output-catalog exports/genomes.updated.json",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::ReferenceUpdateEnsemblSpecs {
            helper_mode,
            catalog_path,
            output_catalog_path,
        } => {
            assert!(!helper_mode);
            assert_eq!(catalog_path, Some("assets/genomes.json".to_string()));
            assert_eq!(
                output_catalog_path,
                Some("exports/genomes.updated.json".to_string())
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_genomes_preview_ensembl_specs() {
    let cmd = parse_shell_line("genomes preview-ensembl-specs --catalog assets/genomes.json")
        .expect("parse command");
    match cmd {
        ShellCommand::ReferencePreviewEnsemblSpecs {
            helper_mode,
            catalog_path,
        } => {
            assert!(!helper_mode);
            assert_eq!(catalog_path, Some("assets/genomes.json".to_string()));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_helpers_remove_catalog_entry() {
    let cmd = parse_shell_line(
        "helpers remove-catalog-entry HelperGenome --catalog helpers.json --output-catalog helpers.updated.json",
    )
    .expect("parse command");
    match cmd {
        ShellCommand::ReferenceRemoveCatalogEntry {
            helper_mode,
            genome_id,
            catalog_path,
            output_catalog_path,
        } => {
            assert!(helper_mode);
            assert_eq!(genome_id, "HelperGenome".to_string());
            assert_eq!(catalog_path, Some("helpers.json".to_string()));
            assert_eq!(
                output_catalog_path,
                Some("helpers.updated.json".to_string())
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn execute_genomes_validate_catalog_reports_valid() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    let cache = td.path().join("cache");
    fs::write(&fasta, ">chr1\nACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
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
        cache.display()
    );
    fs::write(&catalog, catalog_json).expect("write catalog");
    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceValidateCatalog {
            helper_mode: false,
            catalog_path: Some(catalog.to_string_lossy().to_string()),
        },
    )
    .expect("execute validate-catalog");
    assert!(!out.state_changed);
    assert_eq!(out.output["valid"].as_bool(), Some(true));
    assert_eq!(out.output["genome_count"].as_u64(), Some(1));
}

#[test]
fn execute_cutrun_list_reads_catalog_rows() {
    let _serial = cutrun_test_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let td = tempdir().expect("tempdir");
    let catalog = td.path().join("cutrun.catalog.json");
    fs::write(
        &catalog,
        r#"{
  "toy_ctcf": {
    "summary": "Toy CTCF CUT&RUN",
    "target_factor": "CTCF",
    "tissue_or_cell_type": "K562"
  },
  "toy_gata3": {
    "summary": "Toy GATA3 CUT&RUN",
    "target_factor": "GATA3",
    "tissue_or_cell_type": "Jurkat"
  }
}"#,
    )
    .expect("write CUT&RUN catalog");
    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunList {
            filter: Some("k562".to_string()),
            catalog_path: Some(catalog.to_string_lossy().to_string()),
        },
    )
    .expect("execute CUT&RUN list");
    assert!(!out.state_changed);
    assert_eq!(out.output["dataset_count"].as_u64(), Some(1));
    assert_eq!(
        out.output["datasets"][0]["dataset_id"].as_str(),
        Some("toy_ctcf")
    );
}

#[cfg(unix)]
#[test]
fn execute_cutrun_prepare_and_project_return_shared_payloads() {
    let _serial = cutrun_test_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let reference_catalog_path = write_cutrun_shell_reference_catalog(root, "ToyGenome");
    let peaks_path = root.join("toy_peaks.bed");
    let signal_path = root.join("toy_signal.bw");
    let converted_bedgraph = root.join("toy_signal.bedgraph");
    fs::write(&peaks_path, "chr1\t1\t4\tpeak_a\t42\t+\n").expect("write peaks");
    fs::write(&converted_bedgraph, "chr1\t1\t4\t0.5\nchr1\t5\t12\t2.0\n").expect("write bedgraph");
    fs::write(&signal_path, "placeholder").expect("write signal");
    let cutrun_catalog = root.join("cutrun.catalog.json");
    fs::write(
        &cutrun_catalog,
        format!(
            r#"{{
  "toy_ctcf": {{
    "summary": "Toy CUT&RUN",
    "supported_reference_genome_ids": ["ToyGenome"],
    "peaks_local": "{}",
    "signal_local": "{}"
  }}
}}"#,
            peaks_path.display(),
            signal_path.display()
        ),
    )
    .expect("write CUT&RUN catalog");
    let cutrun_cache = root.join("cutrun_cache");
    let fake_converter = install_fake_bigwig_to_bedgraph(root, &converted_bedgraph);
    let _converter_guard = EnvVarGuard::set(BIGWIG_TO_BEDGRAPH_ENV_BIN, &fake_converter);
    let _makeblastdb_guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );

    let mut engine = GentleEngine::new();
    prepare_cutrun_shell_anchor(
        &mut engine,
        "ToyGenome",
        &reference_catalog_path,
        "toy_slice",
    );

    let prepare = execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunPrepare {
            dataset_id: "toy_ctcf".to_string(),
            catalog_path: Some(cutrun_catalog.to_string_lossy().to_string()),
            cache_dir: Some(cutrun_cache.to_string_lossy().to_string()),
        },
    )
    .expect("execute CUT&RUN prepare");
    assert!(!prepare.state_changed);
    assert_eq!(prepare.output["dataset_id"].as_str(), Some("toy_ctcf"));
    assert_eq!(prepare.output["prepared"].as_bool(), Some(true));
    assert_eq!(prepare.output["lifecycle_status"].as_str(), Some("ready"));
    assert_eq!(
        prepare.output["resource_key"].as_str(),
        Some("cutrun_dataset:toy_ctcf")
    );

    let project = execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunProject {
            seq_id: "toy_slice".to_string(),
            dataset_id: "toy_ctcf".to_string(),
            include_peaks: true,
            include_signal: true,
            clear_existing: true,
            catalog_path: Some(cutrun_catalog.to_string_lossy().to_string()),
            cache_dir: Some(cutrun_cache.to_string_lossy().to_string()),
        },
    )
    .expect("execute CUT&RUN project");
    assert!(project.state_changed);
    assert_eq!(project.output["seq_id"].as_str(), Some("toy_slice"));
    assert_eq!(project.output["dataset_id"].as_str(), Some("toy_ctcf"));
    assert!(
        project.output["projected_peak_features"]
            .as_u64()
            .unwrap_or(0)
            > 0
    );
    assert!(
        project.output["projected_signal_features"]
            .as_u64()
            .unwrap_or(0)
            > 0
    );
}

#[test]
fn execute_cutrun_status_reports_running_lifecycle_from_active_prepare_marker() {
    let _serial = cutrun_test_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let peaks_path = root.join("toy_peaks.bed");
    fs::write(&peaks_path, "chr1\t1\t4\tpeak_a\t42\t+\n").expect("write peaks");
    let cutrun_catalog = root.join("cutrun.catalog.json");
    fs::write(
        &cutrun_catalog,
        format!(
            r#"{{
  "toy_ctcf": {{
    "summary": "Toy CUT&RUN",
    "target_factor": "CTCF",
    "peaks_local": "{}"
  }}
}}"#,
            peaks_path.display()
        ),
    )
    .expect("write CUT&RUN catalog");
    let cutrun_cache = root.join("cutrun_cache");
    let install_dir = cutrun_cache.join("toy_ctcf");
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    write_cutrun_shell_prepare_activity(&install_dir, "toy_ctcf", "running", now, None);

    let mut engine = GentleEngine::new();
    let status = execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunStatus {
            dataset_id: "toy_ctcf".to_string(),
            catalog_path: Some(cutrun_catalog.to_string_lossy().to_string()),
            cache_dir: Some(cutrun_cache.to_string_lossy().to_string()),
        },
    )
    .expect("execute CUT&RUN status");
    assert!(!status.state_changed);
    assert_eq!(status.output["dataset_id"].as_str(), Some("toy_ctcf"));
    assert_eq!(status.output["prepared"].as_bool(), Some(false));
    assert_eq!(status.output["lifecycle_status"].as_str(), Some("running"));
    assert_eq!(
        status.output["resource_key"].as_str(),
        Some("cutrun_dataset:toy_ctcf")
    );
    assert_eq!(
        status.output["current_activity"]["phase"].as_str(),
        Some("materializing_asset")
    );
}

#[test]
fn execute_cutrun_interpret_list_show_and_export_coverage() {
    let _serial = cutrun_test_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let reference_catalog_path = write_cutrun_shell_reference_catalog_with_sequence(
        root,
        "ToyGenome",
        "AAAGCCGTAGCTTACGGAACCTTT",
    );
    let input_r1 = root.join("reads_r1.fastq");
    let input_r2 = root.join("reads_r2.fastq");
    fs::write(
        &input_r1,
        "@pair1/1\nGCCG\n+\nIIII\n@orphan_r1/1\nTAGC\n+\nIIII\n",
    )
    .expect("write CUT&RUN shell R1 FASTQ");
    fs::write(
        &input_r2,
        "@pair1/2\nGTAA\n+\nIIII\n@orphan_r2/2\nACGG\n+\nIIII\n",
    )
    .expect("write CUT&RUN shell R2 FASTQ");
    let _makeblastdb_guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );

    let mut engine = GentleEngine::new();
    prepare_cutrun_shell_anchor_range(
        &mut engine,
        "ToyGenome",
        &reference_catalog_path,
        "toy_cutrun_roi",
        4,
        15,
    );

    let interpret = execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunInterpret {
            seq_id: "toy_cutrun_roi".to_string(),
            input_r1_path: Some(input_r1.to_string_lossy().to_string()),
            input_r2_path: Some(input_r2.to_string_lossy().to_string()),
            dataset_id: None,
            catalog_path: None,
            cache_dir: None,
            input_format: CutRunInputFormat::Fastq,
            read_layout: CutRunReadLayout::PairedEnd,
            roi_flank_bp: 0,
            seed_filter: CutRunSeedFilterConfig {
                kmer_len: 2,
                min_seed_matches: 1,
            },
            align_config: CutRunAlignConfig {
                max_mismatches: 0,
                min_identity_fraction: 1.0,
                max_fragment_span_bp: 32,
            },
            deduplicate_fragments: false,
            report_id: Some("toy_cutrun_reads".to_string()),
            checkpoint_path: None,
            checkpoint_every_reads: 10,
        },
    )
    .expect("execute CUT&RUN interpret");
    assert!(interpret.state_changed);
    assert_eq!(
        interpret.output["report"]["report_id"].as_str(),
        Some("toy_cutrun_reads")
    );
    assert_eq!(
        interpret.output["report"]["fragment_count"].as_u64(),
        Some(3)
    );
    assert_eq!(
        interpret.output["report"]["concordant_pair_count"].as_u64(),
        Some(1)
    );

    let listed = execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunListReadReports {
            seq_id: Some("toy_cutrun_roi".to_string()),
        },
    )
    .expect("execute CUT&RUN list-read-reports");
    assert!(!listed.state_changed);
    assert_eq!(listed.output["report_count"].as_u64(), Some(1));
    assert_eq!(
        listed.output["reports"][0]["report_id"].as_str(),
        Some("toy_cutrun_reads")
    );

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunShowReadReport {
            report_id: "toy_cutrun_reads".to_string(),
        },
    )
    .expect("execute CUT&RUN show-read-report");
    assert!(!shown.state_changed);
    assert_eq!(
        shown.output["report"]["seq_id"].as_str(),
        Some("toy_cutrun_roi")
    );
    assert_eq!(shown.output["report"]["mapped_units"].as_u64(), Some(3));
    assert_eq!(shown.output["report"]["orphan_r1_count"].as_u64(), Some(1));
    assert_eq!(shown.output["report"]["orphan_r2_count"].as_u64(), Some(1));

    let export_path = root.join("coverage.tsv");
    let exported = execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunExportCoverage {
            report_id: "toy_cutrun_reads".to_string(),
            path: export_path.to_string_lossy().to_string(),
            kind: CutRunCoverageKind::Coverage,
        },
    )
    .expect("execute CUT&RUN export-coverage");
    assert!(!exported.state_changed);
    assert_eq!(
        exported.output["report_id"].as_str(),
        Some("toy_cutrun_reads")
    );
    assert_eq!(exported.output["row_count"].as_u64(), Some(12));
    let coverage_text = fs::read_to_string(&export_path).expect("read coverage export");
    assert_eq!(
        coverage_text.lines().next(),
        Some("local_pos_1based\tgenomic_pos_1based\tcoverage")
    );
    assert!(coverage_text.contains("1\t4\t1"));
    assert!(coverage_text.contains("5\t8\t3"));
    assert!(coverage_text.contains("12\t15\t1"));
}

#[test]
fn execute_cutrun_interpret_from_prepared_dataset_reads() {
    let _serial = cutrun_test_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let reference_catalog_path = write_cutrun_shell_reference_catalog_with_sequence(
        root,
        "ToyGenome",
        "AAAGCCGTAGCTTACGGAACCTTT",
    );
    let reads_r1_path = root.join("toy_reads_r1.fastq");
    let reads_r2_path = root.join("toy_reads_r2.fastq");
    fs::write(
        &reads_r1_path,
        "@pair1/1\nGCCG\n+\nIIII\n@orphan_r1/1\nTAGC\n+\nIIII\n",
    )
    .expect("write shell dataset reads r1");
    fs::write(
        &reads_r2_path,
        "@pair1/2\nGTAA\n+\nIIII\n@orphan_r2/2\nACGG\n+\nIIII\n",
    )
    .expect("write shell dataset reads r2");
    let cutrun_catalog = root.join("cutrun.catalog.json");
    fs::write(
        &cutrun_catalog,
        format!(
            r#"{{
  "toy_ctcf_reads": {{
    "summary": "Toy CUT&RUN raw reads",
    "supported_reference_genome_ids": ["ToyGenome"],
    "reads_r1_local": "{}",
    "reads_r2_local": "{}",
    "read_layout": "paired_end"
  }}
}}"#,
            reads_r1_path.display(),
            reads_r2_path.display()
        ),
    )
    .expect("write shell CUT&RUN raw-read catalog");
    let cutrun_cache = root.join("cutrun_cache");
    let _makeblastdb_guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );

    let mut engine = GentleEngine::new();
    prepare_cutrun_shell_anchor_range(
        &mut engine,
        "ToyGenome",
        &reference_catalog_path,
        "toy_cutrun_roi",
        4,
        15,
    );
    execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunPrepare {
            dataset_id: "toy_ctcf_reads".to_string(),
            catalog_path: Some(cutrun_catalog.to_string_lossy().to_string()),
            cache_dir: Some(cutrun_cache.to_string_lossy().to_string()),
        },
    )
    .expect("prepare shell CUT&RUN raw-read dataset");

    let interpret = execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunInterpret {
            seq_id: "toy_cutrun_roi".to_string(),
            input_r1_path: None,
            input_r2_path: None,
            dataset_id: Some("toy_ctcf_reads".to_string()),
            catalog_path: Some(cutrun_catalog.to_string_lossy().to_string()),
            cache_dir: Some(cutrun_cache.to_string_lossy().to_string()),
            input_format: CutRunInputFormat::Fasta,
            read_layout: CutRunReadLayout::SingleEnd,
            roi_flank_bp: 0,
            seed_filter: CutRunSeedFilterConfig {
                kmer_len: 2,
                min_seed_matches: 1,
            },
            align_config: CutRunAlignConfig {
                max_mismatches: 0,
                min_identity_fraction: 1.0,
                max_fragment_span_bp: 32,
            },
            deduplicate_fragments: false,
            report_id: Some("toy_cutrun_reads_from_dataset".to_string()),
            checkpoint_path: None,
            checkpoint_every_reads: 10,
        },
    )
    .expect("execute dataset-backed CUT&RUN interpret");
    assert!(interpret.state_changed);
    assert_eq!(
        interpret.output["report"]["report_id"].as_str(),
        Some("toy_cutrun_reads_from_dataset")
    );
    assert_eq!(
        interpret.output["report"]["input_format"].as_str(),
        Some("fastq")
    );
    assert_eq!(
        interpret.output["report"]["read_layout"].as_str(),
        Some("paired_end")
    );
    assert!(
        interpret.output["report"]["warnings"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(|value| value.as_str())
            .any(|value| value.contains("resolved CUT&RUN raw reads from prepared dataset"))
    );
}

#[test]
fn execute_cutrun_inspect_regulatory_support_writes_json_payload() {
    let _serial = cutrun_test_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let td = tempdir().expect("tempdir");
    let root = td.path();
    let reference_catalog_path = write_cutrun_shell_reference_catalog_with_sequence(
        root,
        "ToyGenome",
        "AAAGCCGTAGCTTACGGAACCTTT",
    );
    let input_r1 = root.join("reads.fasta");
    fs::write(&input_r1, ">read_a\nGCCGTAGC\n>read_b\nGCCGTAGC\n")
        .expect("write CUT&RUN shell FASTA");
    let _makeblastdb_guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );

    let mut engine = GentleEngine::new();
    prepare_cutrun_shell_anchor_range(
        &mut engine,
        "ToyGenome",
        &reference_catalog_path,
        "toy_cutrun_roi",
        4,
        15,
    );
    execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunInterpret {
            seq_id: "toy_cutrun_roi".to_string(),
            input_r1_path: Some(input_r1.to_string_lossy().to_string()),
            input_r2_path: None,
            dataset_id: None,
            catalog_path: None,
            cache_dir: None,
            input_format: CutRunInputFormat::Fasta,
            read_layout: CutRunReadLayout::SingleEnd,
            roi_flank_bp: 0,
            seed_filter: CutRunSeedFilterConfig {
                kmer_len: 3,
                min_seed_matches: 1,
            },
            align_config: CutRunAlignConfig {
                max_mismatches: 0,
                min_identity_fraction: 1.0,
                max_fragment_span_bp: 32,
            },
            deduplicate_fragments: false,
            report_id: Some("toy_cutrun_reads".to_string()),
            checkpoint_path: None,
            checkpoint_every_reads: 10,
        },
    )
    .expect("execute shell CUT&RUN interpret");

    let output_path = root.join("cutrun_support.json");
    let inspected = execute_shell_command(
        &mut engine,
        &ShellCommand::CutRunInspectRegulatorySupport {
            seq_id: "toy_cutrun_roi".to_string(),
            dataset_ids: vec![],
            read_report_ids: vec!["toy_cutrun_reads".to_string()],
            promoter_search_start_0based: None,
            promoter_search_end_0based_exclusive: None,
            neighbor_window_bp: 150,
            species_filters: vec![],
            output_path: Some(output_path.to_string_lossy().to_string()),
        },
    )
    .expect("execute CUT&RUN inspect-regulatory-support");
    assert!(!inspected.state_changed);
    assert_eq!(
        inspected.output["schema"].as_str(),
        Some("gentle.cutrun_regulatory_support.v1")
    );
    assert_eq!(
        inspected.output["evidence_sources"][0]["report_id"].as_str(),
        Some("toy_cutrun_reads")
    );
    assert_eq!(
        inspected.output["support_windows"].as_array().map(Vec::len),
        Some(1)
    );
    let written = fs::read_to_string(&output_path).expect("read written regulatory-support JSON");
    assert!(written.contains("\"schema\": \"gentle.cutrun_regulatory_support.v1\""));
    assert!(written.contains("\"report_id\": \"toy_cutrun_reads\""));
}

#[test]
fn execute_genomes_update_ensembl_specs_with_no_templates_reports_no_updates() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    let cache = td.path().join("cache");
    fs::write(&fasta, ">chr1\nACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
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
        cache.display()
    );
    fs::write(&catalog, catalog_json).expect("write catalog");

    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceUpdateEnsemblSpecs {
            helper_mode: false,
            catalog_path: Some(catalog.to_string_lossy().to_string()),
            output_catalog_path: None,
        },
    )
    .expect("execute update-ensembl-specs");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["report"]["updates"].as_array().map(Vec::len),
        Some(0)
    );
    assert_eq!(out.output["report"]["wrote_catalog"].as_bool(), Some(true));
    assert_eq!(
        out.output["report"]["output_catalog_path"].as_str(),
        Some(catalog.to_string_lossy().as_ref())
    );
}

#[test]
fn execute_genomes_preview_ensembl_specs_with_no_templates_reports_no_updates() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    let cache = td.path().join("cache");
    fs::write(&fasta, ">chr1\nACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
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
        cache.display()
    );
    fs::write(&catalog, catalog_json).expect("write catalog");

    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferencePreviewEnsemblSpecs {
            helper_mode: false,
            catalog_path: Some(catalog.to_string_lossy().to_string()),
        },
    )
    .expect("execute preview-ensembl-specs");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["preview"]["updates"].as_array().map(Vec::len),
        Some(0)
    );
    assert_eq!(
        out.output["preview"]["writable_catalog"].as_bool(),
        Some(true)
    );
}

#[test]
fn execute_helpers_list_filters_semantic_metadata() {
    let td = tempdir().expect("tempdir");
    let catalog = td.path().join("helpers.json");
    fs::write(
        &catalog,
        r#"{
  "pGEX_like_vector": {
    "sequence_remote": "https://example.invalid/pgex.fa.gz",
    "annotations_remote": "https://example.invalid/pgex.gb.gz",
    "summary": "GST fusion vector with Factor Xa cleavage site",
    "aliases": ["pGEX"],
    "helper_kind": "plasmid_vector",
    "host_system": "Escherichia coli",
    "search_terms": ["factor xa", "affinity purification"],
    "semantics": {
      "schema": "gentle.helper_semantics.v1",
      "affordances": ["bacterial_expression", "protease_tag_removal"],
      "components": [
        {
          "id": "mcs",
          "kind": "cloning_site",
          "label": "MCS"
        }
      ]
    }
  },
  "neutral_backbone": {
    "sequence_remote": "https://example.invalid/neutral.fa.gz",
    "annotations_remote": "https://example.invalid/neutral.gb.gz",
    "summary": "Simple backbone"
  }
}"#,
    )
    .expect("write catalog");
    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceList {
            helper_mode: true,
            catalog_path: Some(catalog.to_string_lossy().to_string()),
            filter: Some("factor xa".to_string()),
        },
    )
    .expect("execute helpers list");
    assert!(!out.state_changed);
    assert_eq!(out.output["genome_count"].as_u64(), Some(1));
    assert_eq!(out.output["genomes"][0].as_str(), Some("pGEX_like_vector"));
    assert_eq!(
        out.output["entries"][0]["helper_kind"].as_str(),
        Some("plasmid_vector")
    );
    let offered_functions = out.output["entries"][0]["interpretation"]["offered_functions"]
        .as_array()
        .expect("interpretation functions");
    assert_eq!(
        offered_functions
            .iter()
            .filter_map(|value| value.as_str())
            .collect::<Vec<_>>(),
        vec![
            "bacterial_expression",
            "insert_cloning",
            "protease_tag_removal"
        ]
    );
}

#[test]
fn execute_helpers_status_includes_normalized_interpretation() {
    let td = tempdir().expect("tempdir");
    let catalog = td.path().join("helpers.json");
    fs::write(
        &catalog,
        r#"{
  "Plasmid pUC19 (online)": {
    "sequence_remote": "https://example.invalid/puc19.fa.gz",
    "annotations_remote": "https://example.invalid/puc19.gb.gz",
    "summary": "High-copy cloning vector",
    "helper_kind": "plasmid_vector",
    "host_system": "Escherichia coli",
    "semantics": {
      "schema": "gentle.helper_semantics.v1",
      "affordances": ["insert_cloning", "ampicillin_selection"],
      "components": [
        {
          "id": "bla",
          "kind": "selectable_marker",
          "label": "AmpR",
          "attributes": {"agent": "ampicillin"}
        }
      ]
    }
  }
}"#,
    )
    .expect("write catalog");

    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceStatus {
            helper_mode: true,
            genome_id: "Plasmid pUC19 (online)".to_string(),
            catalog_path: Some(catalog.to_string_lossy().to_string()),
            cache_dir: None,
        },
    )
    .expect("execute helpers status");

    assert!(!out.state_changed);
    assert_eq!(
        out.output["interpretation"]["helper_kinds"][0].as_str(),
        Some("plasmid_vector")
    );
    assert_eq!(
        out.output["interpretation"]["host_systems"][0].as_str(),
        Some("Escherichia coli")
    );
    assert_eq!(
        out.output["interpretation"]["offered_functions"]
            .as_array()
            .expect("offered functions")
            .iter()
            .filter_map(|value| value.as_str())
            .collect::<Vec<_>>(),
        vec!["ampicillin_selection", "insert_cloning", "selection"]
    );
}

#[test]
fn execute_hosts_list_filters_trait_metadata() {
    let td = tempdir().expect("tempdir");
    let catalog = td.path().join("host_profiles.json");
    fs::write(
        &catalog,
        r#"{
  "schema": "gentle.host_profile_catalog.v1",
  "profiles": [
    {
      "profile_id": "ecoli_dh5alpha",
      "species": "Escherichia coli",
      "strain": "DH5alpha",
      "aliases": ["DH5α"],
      "genotype_tags": ["hsdR17", "deoR", "relA1"],
      "phenotype_tags": ["cloning_host", "large_insert_friendly"],
      "notes": ["Routine cloning host with deoR-associated large-insert friendliness."],
      "source_notes": ["Synthetic test fixture"]
    },
    {
      "profile_id": "ecoli_k12_restriction_positive",
      "species": "Escherichia coli",
      "strain": "K-12 restriction-positive background",
      "aliases": ["E. coli K-12"],
      "genotype_tags": ["hsdR+", "hsdM+"],
      "phenotype_tags": ["type_i_restriction_barrier"],
      "notes": ["Restriction-positive route review host."],
      "source_notes": ["Synthetic test fixture"]
    }
  ]
}"#,
    )
    .expect("write host catalog");

    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::HostsList {
            catalog_path: Some(catalog.to_string_lossy().to_string()),
            filter: Some("deoR".to_string()),
        },
    )
    .expect("execute hosts list");

    assert!(!out.state_changed);
    assert_eq!(out.output["profile_count"].as_u64(), Some(1));
    assert_eq!(
        out.output["profile_ids"][0].as_str(),
        Some("ecoli_dh5alpha")
    );
    assert_eq!(
        out.output["entries"][0]["phenotype_tags"][0].as_str(),
        Some("cloning_host")
    );
}

#[test]
fn execute_genomes_remove_prepared_reports_removed_and_status_clears() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    let cache = td.path().join("cache");
    fs::write(&fasta, ">chr1\nACGTACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
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
        cache.display()
    );
    fs::write(&catalog, catalog_json).expect("write catalog");
    let catalog_path = catalog.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    execute_shell_command(
        &mut engine,
        &ShellCommand::ReferencePrepare {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.clone()),
            cache_dir: None,
            timeout_seconds: None,
        },
    )
    .expect("prepare genome");

    let removed = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceRemovePrepared {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.clone()),
            cache_dir: None,
        },
    )
    .expect("remove prepared genome");
    assert!(!removed.state_changed);
    assert_eq!(removed.output["report"]["removed"].as_bool(), Some(true));

    let status = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceStatus {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path),
            cache_dir: None,
        },
    )
    .expect("status after remove");
    assert_eq!(status.output["prepared"].as_bool(), Some(false));
}

#[test]
fn execute_genomes_remove_catalog_entry_updates_catalog_file() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    let cache = td.path().join("cache");
    fs::write(&fasta, ">chr1\nACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
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
        cache.display()
    );
    fs::write(&catalog, catalog_json).expect("write catalog");
    let catalog_path = catalog.to_string_lossy().to_string();

    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceRemoveCatalogEntry {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.clone()),
            output_catalog_path: None,
        },
    )
    .expect("remove catalog entry");
    assert!(!out.state_changed);
    assert_eq!(out.output["report"]["removed"].as_bool(), Some(true));

    let genomes =
        GentleEngine::list_reference_genomes(Some(&catalog_path)).expect("reload edited catalog");
    assert!(genomes.is_empty());
}

#[test]
fn execute_genomes_status_reports_length_and_mass_metadata() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    let cache = td.path().join("cache");
    fs::write(&fasta, ">chr1\nACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}",
    "nucleotide_length_bp": 4
  }}
}}"#,
        fasta.display(),
        gtf.display(),
        cache.display()
    );
    fs::write(&catalog, catalog_json).expect("write catalog");
    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceStatus {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog.to_string_lossy().to_string()),
            cache_dir: None,
        },
    )
    .expect("execute status");
    assert!(!out.state_changed);
    assert_eq!(out.output["nucleotide_length_bp"].as_u64(), Some(4));
    let expected_mass = 4.0 * 617.96 + 36.04;
    let observed_mass = out.output["molecular_mass_da"]
        .as_f64()
        .expect("molecular_mass_da should be present");
    assert!((observed_mass - expected_mass).abs() < 1e-6);
    assert_eq!(
        out.output["molecular_mass_source"].as_str(),
        Some("estimated_from_nucleotide_length")
    );
}

#[test]
fn execute_genomes_status_reports_effective_cache_dir_and_prepare_hint_when_unprepared() {
    let _lock = cache_dir_env_lock().lock().expect("cache dir env lock");
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    let shared_cache = td.path().join("shared_ref_cache");
    fs::write(&fasta, ">chr1\nACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
    fs::write(
        &catalog,
        format!(
            r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "data/genomes"
  }}
}}"#,
            fasta.display(),
            gtf.display()
        ),
    )
    .expect("write catalog");
    let _guard = EnvVarGuard::set(
        crate::genomes::REFERENCE_GENOME_CACHE_DIR_ENV,
        shared_cache.to_string_lossy().as_ref(),
    );

    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceStatus {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog.to_string_lossy().to_string()),
            cache_dir: None,
        },
    )
    .expect("execute status");

    assert_eq!(out.output["prepared"].as_bool(), Some(false));
    assert_eq!(
        out.output["effective_cache_dir"].as_str(),
        Some(shared_cache.to_string_lossy().as_ref())
    );
    assert_eq!(
        out.output["requested_catalog_key"].as_str(),
        Some("ToyGenome")
    );
    assert_eq!(
        out.output["compatible_prepared_options"].as_array(),
        Some(&vec![])
    );
    let status_message = out.output["status_message"]
        .as_str()
        .expect("status message");
    assert!(status_message.contains("not prepared"));
    assert!(status_message.contains(shared_cache.to_string_lossy().as_ref()));
    let prepare_command = out.output["prepare_command"]
        .as_str()
        .expect("prepare command");
    assert!(prepare_command.contains("genomes prepare"));
    assert!(prepare_command.contains("ToyGenome"));
    assert!(prepare_command.contains(shared_cache.to_string_lossy().as_ref()));
    assert_eq!(out.output["lifecycle_status"].as_str(), Some("missing"));
}

#[test]
fn execute_genomes_status_reports_running_lifecycle_and_suppresses_prepare_hint() {
    let td = tempdir().expect("tempdir");
    let fasta = td.path().join("toy.fa");
    let gtf = td.path().join("toy.gtf");
    let cache = td.path().join("cache");
    fs::write(&fasta, ">chr1\nACGT\n").expect("write fasta");
    fs::write(
        &gtf,
        "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
    )
    .expect("write gtf");
    let catalog = td.path().join("catalog.json");
    fs::write(
        &catalog,
        format!(
            r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            fasta.display(),
            gtf.display(),
            cache.display()
        ),
    )
    .expect("write catalog");
    let install_dir = cache.join("toygenome");
    fs::create_dir_all(&install_dir).expect("install dir");
    let status_path = install_dir.join(".prepare_activity.json");
    let lock_path = install_dir.join(".prepare_activity.lock");
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("now")
        .as_millis();
    let activity = serde_json::json!({
        "genome_id": "ToyGenome",
        "status_path": status_path.display().to_string(),
        "lock_path": lock_path.display().to_string(),
        "lifecycle_status": "running",
        "prepare_mode": "prepare_or_reuse",
        "phase": "download_sequence",
        "item": "toy.fa",
        "bytes_done": 10,
        "bytes_total": 100,
        "percent": 10.0,
        "step_id": null,
        "step_label": null,
        "started_at_unix_ms": now.saturating_sub(10),
        "updated_at_unix_ms": now,
        "finished_at_unix_ms": null,
        "last_error": null,
        "owner_pid": 12345
    });
    fs::write(
        &status_path,
        serde_json::to_string_pretty(&activity).expect("serialize activity"),
    )
    .expect("status");
    fs::write(
        &lock_path,
        serde_json::to_string_pretty(&activity).expect("serialize activity lock"),
    )
    .expect("lock");

    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::ReferenceStatus {
            helper_mode: false,
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog.to_string_lossy().to_string()),
            cache_dir: None,
        },
    )
    .expect("execute status");

    assert_eq!(out.output["lifecycle_status"].as_str(), Some("running"));
    assert_eq!(
        out.output["current_activity"]["phase"].as_str(),
        Some("download_sequence")
    );
    assert!(out.output["prepare_command"].is_null());
    assert!(
        out.output["status_message"]
            .as_str()
            .is_some_and(|value| value.contains("already being prepared"))
    );
}

#[test]
fn parse_ui_open_and_prepared_commands() {
    let open = parse_shell_line("ui open prepared-references --genome-id \"Human GRCh38\"")
        .expect("parse ui open");
    match open {
        ShellCommand::UiIntent {
            action,
            target,
            genome_id,
            helper_mode,
            catalog_path,
            cache_dir,
            filter,
            species,
            latest,
        } => {
            assert_eq!(action, UiIntentAction::Open);
            assert_eq!(target, UiIntentTarget::PreparedReferences);
            assert_eq!(genome_id.as_deref(), Some("Human GRCh38"));
            assert!(!helper_mode);
            assert!(catalog_path.is_none());
            assert!(cache_dir.is_none());
            assert!(filter.is_none());
            assert!(species.is_none());
            assert!(!latest);
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let open_latest = parse_shell_line("ui open prepared-references --species human --latest")
        .expect("parse ui open latest");
    match open_latest {
        ShellCommand::UiIntent {
            action,
            target,
            genome_id,
            helper_mode,
            species,
            latest,
            ..
        } => {
            assert_eq!(action, UiIntentAction::Open);
            assert_eq!(target, UiIntentTarget::PreparedReferences);
            assert!(genome_id.is_none());
            assert!(!helper_mode);
            assert_eq!(species.as_deref(), Some("human"));
            assert!(latest);
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let open_pcr = parse_shell_line("ui open pcr-design").expect("parse ui open pcr-design");
    match open_pcr {
        ShellCommand::UiIntent { action, target, .. } => {
            assert_eq!(action, UiIntentAction::Open);
            assert_eq!(target, UiIntentTarget::PcrDesign);
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let focus_pcr = parse_shell_line("ui focus pcr-design").expect("parse ui focus pcr-design");
    match focus_pcr {
        ShellCommand::UiIntent { action, target, .. } => {
            assert_eq!(action, UiIntentAction::Focus);
            assert_eq!(target, UiIntentTarget::PcrDesign);
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let open_seq_confirm = parse_shell_line("ui open sequencing-confirmation")
        .expect("parse ui open sequencing-confirmation");
    match open_seq_confirm {
        ShellCommand::UiIntent { action, target, .. } => {
            assert_eq!(action, UiIntentAction::Open);
            assert_eq!(target, UiIntentTarget::SequencingConfirmation);
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let focus_seq_confirm =
        parse_shell_line("ui focus seq-confirm").expect("parse ui focus seq-confirm");
    match focus_seq_confirm {
        ShellCommand::UiIntent { action, target, .. } => {
            assert_eq!(action, UiIntentAction::Focus);
            assert_eq!(target, UiIntentTarget::SequencingConfirmation);
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let prepared = parse_shell_line("ui prepared-genomes --species human --latest")
        .expect("parse ui prepared-genomes");
    match prepared {
        ShellCommand::UiPreparedGenomes {
            helper_mode,
            species,
            latest,
            ..
        } => {
            assert!(!helper_mode);
            assert_eq!(species.as_deref(), Some("human"));
            assert!(latest);
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let err = parse_shell_line("ui open agent-assistant --latest")
        .expect_err("ui open agent-assistant --latest should fail");
    assert!(
            err.contains("only supports --helpers/--catalog/--cache-dir/--filter/--species/--latest when TARGET is prepared-references"),
            "unexpected parse error: {err}"
        );
}

#[test]
fn execute_ui_prepared_and_latest_prepared_queries() {
    let td = tempdir().expect("tempdir");
    let fasta_113 = td.path().join("h113.fa");
    let ann_113 = td.path().join("h113.gtf");
    let fasta_116 = td.path().join("h116.fa");
    let ann_116 = td.path().join("h116.gtf");
    let cache_dir = td.path().join("cache");
    fs::write(&fasta_113, ">chr1\nACGT\n").expect("write fasta 113");
    fs::write(
        &ann_113,
        "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"G113\"; gene_name \"G113\";\n",
    )
    .expect("write ann 113");
    fs::write(&fasta_116, ">chr1\nACGTACGT\n").expect("write fasta 116");
    fs::write(
        &ann_116,
        "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"G116\"; gene_name \"G116\";\n",
    )
    .expect("write ann 116");
    let catalog_path = td.path().join("catalog.json");
    let catalog_json = format!(
        r#"{{
  "Human GRCh38 Ensembl 113": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }},
  "Human GRCh38 Ensembl 116": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
        fasta_113.display(),
        ann_113.display(),
        cache_dir.display(),
        fasta_116.display(),
        ann_116.display(),
        cache_dir.display()
    );
    fs::write(&catalog_path, catalog_json).expect("write catalog");

    let catalog = GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref()).unwrap();
    catalog
        .prepare_genome_once("Human GRCh38 Ensembl 113")
        .expect("prepare 113");
    std::thread::sleep(std::time::Duration::from_millis(2));
    catalog
        .prepare_genome_once("Human GRCh38 Ensembl 116")
        .expect("prepare 116");

    let mut engine = GentleEngine::new();
    let prepared = execute_shell_command(
        &mut engine,
        &ShellCommand::UiPreparedGenomes {
            helper_mode: false,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            filter: None,
            species: Some("human".to_string()),
            latest: false,
        },
    )
    .expect("execute ui prepared-genomes");
    assert!(!prepared.state_changed);
    assert_eq!(prepared.output["prepared_count"].as_u64(), Some(2));
    assert_eq!(
        prepared.output["genomes"][0]["genome_id"].as_str(),
        Some("Human GRCh38 Ensembl 116")
    );

    let latest = execute_shell_command(
        &mut engine,
        &ShellCommand::UiLatestPrepared {
            helper_mode: false,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            species: "human".to_string(),
        },
    )
    .expect("execute ui latest-prepared");
    assert!(!latest.state_changed);
    assert_eq!(
        latest.output["selected_genome_id"].as_str(),
        Some("Human GRCh38 Ensembl 116")
    );

    let intent_latest = execute_shell_command(
        &mut engine,
        &ShellCommand::UiIntent {
            action: UiIntentAction::Open,
            target: UiIntentTarget::PreparedReferences,
            genome_id: None,
            helper_mode: false,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            filter: None,
            species: Some("human".to_string()),
            latest: true,
        },
    )
    .expect("execute ui intent latest");
    assert!(!intent_latest.state_changed);
    assert_eq!(
        intent_latest.output["selected_genome_id"].as_str(),
        Some("Human GRCh38 Ensembl 116")
    );

    let intent_explicit = execute_shell_command(
        &mut engine,
        &ShellCommand::UiIntent {
            action: UiIntentAction::Open,
            target: UiIntentTarget::PreparedReferences,
            genome_id: Some("Human GRCh38 Ensembl 113".to_string()),
            helper_mode: false,
            catalog_path: Some(catalog_path.to_string_lossy().to_string()),
            cache_dir: None,
            filter: None,
            species: Some("human".to_string()),
            latest: true,
        },
    )
    .expect("execute ui intent explicit genome");
    assert!(!intent_explicit.state_changed);
    assert_eq!(
        intent_explicit.output["selected_genome_id"].as_str(),
        Some("Human GRCh38 Ensembl 113")
    );
    assert!(
        intent_explicit.output["prepared_query"].is_null(),
        "explicit genome_id should bypass prepared query resolution"
    );
}

#[test]
fn parse_set_param_command() {
    let cmd = parse_shell_line("set-param vcf_display_pass_only true").expect("parse set-param");
    match cmd {
        ShellCommand::SetParameter { name, value_json } => {
            assert_eq!(name, "vcf_display_pass_only");
            assert_eq!(value_json, "true");
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_feature_expert_commands() {
    let inspect =
        parse_shell_line("inspect-feature-expert s tfbs 7").expect("parse inspect-feature-expert");
    match inspect {
        ShellCommand::InspectFeatureExpert { seq_id, target } => {
            assert_eq!(seq_id, "s");
            assert_eq!(target, FeatureExpertTarget::TfbsFeature { feature_id: 7 });
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let splicing =
        parse_shell_line("inspect-feature-expert s splicing 11").expect("parse splicing");
    match splicing {
        ShellCommand::InspectFeatureExpert { seq_id, target } => {
            assert_eq!(seq_id, "s");
            assert_eq!(
                target,
                FeatureExpertTarget::SplicingFeature {
                    feature_id: 11,
                    scope: SplicingScopePreset::AllOverlappingBothStrands,
                }
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let isoform = parse_shell_line("inspect-feature-expert s isoform tp53_isoforms_v1")
        .expect("parse isoform feature target");
    match isoform {
        ShellCommand::InspectFeatureExpert { seq_id, target } => {
            assert_eq!(seq_id, "s");
            assert_eq!(
                target,
                FeatureExpertTarget::IsoformArchitecture {
                    panel_id: "tp53_isoforms_v1".to_string()
                }
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let protein_compare = parse_shell_line("inspect-feature-expert s protein-comparison")
        .expect("parse transcript protein comparison target");
    match protein_compare {
        ShellCommand::InspectFeatureExpert { seq_id, target } => {
            assert_eq!(seq_id, "s");
            assert_eq!(
                target,
                FeatureExpertTarget::ProteinComparison {
                    transcript_id_filter: None,
                    protein_feature_filter: Default::default(),
                    external_source: None,
                    external_entry_id: None,
                }
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let filtered_protein_compare =
        parse_shell_line("inspect-feature-expert s protein-comparison --transcript TX_TPROT")
            .expect("parse filtered transcript protein comparison target");
    match filtered_protein_compare {
        ShellCommand::InspectFeatureExpert { seq_id, target } => {
            assert_eq!(seq_id, "s");
            assert_eq!(
                target,
                FeatureExpertTarget::ProteinComparison {
                    transcript_id_filter: Some("TX_TPROT".to_string()),
                    protein_feature_filter: Default::default(),
                    external_source: None,
                    external_entry_id: None,
                }
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let ensembl_protein_compare = parse_shell_line(
        "inspect-feature-expert s protein-comparison --transcript TX_TPROT --ensembl-entry ENSPTOY1 --feature-key PF02196 --feature-key-not SEG",
    )
    .expect("parse Ensembl-backed transcript protein comparison target");
    match ensembl_protein_compare {
        ShellCommand::InspectFeatureExpert { seq_id, target } => {
            assert_eq!(seq_id, "s");
            assert_eq!(
                target,
                FeatureExpertTarget::ProteinComparison {
                    transcript_id_filter: Some("TX_TPROT".to_string()),
                    protein_feature_filter: ProteinFeatureFilter {
                        include_feature_keys: vec!["PF02196".to_string()],
                        exclude_feature_keys: vec!["SEG".to_string()],
                    },
                    external_source: Some(ProteinExternalOpinionSource::Ensembl),
                    external_entry_id: Some("ENSPTOY1".to_string()),
                }
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let uniprot = parse_shell_line("inspect-feature-expert s uniprot-projection PTEST1@seq_u")
        .expect("parse uniprot projection feature target");
    match uniprot {
        ShellCommand::InspectFeatureExpert { seq_id, target } => {
            assert_eq!(seq_id, "s");
            assert_eq!(
                target,
                FeatureExpertTarget::UniprotProjection {
                    projection_id: "PTEST1@seq_u".to_string(),
                    protein_feature_filter: Default::default(),
                }
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let filtered_uniprot = parse_shell_line(
        "inspect-feature-expert s uniprot-projection PTEST1@seq_u --feature-key DOMAIN --feature-key REGION --feature-key-not VARIANT",
    )
    .expect("parse filtered uniprot projection feature target");
    match filtered_uniprot {
        ShellCommand::InspectFeatureExpert { seq_id, target } => {
            assert_eq!(seq_id, "s");
            assert_eq!(
                target,
                FeatureExpertTarget::UniprotProjection {
                    projection_id: "PTEST1@seq_u".to_string(),
                    protein_feature_filter: ProteinFeatureFilter {
                        include_feature_keys: vec!["DOMAIN".to_string(), "REGION".to_string()],
                        exclude_feature_keys: vec!["VARIANT".to_string()],
                    },
                }
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let render = parse_shell_line(
        "render-feature-expert-svg s restriction 123 --enzyme EcoRI --start 100 --end 106 out.svg",
    )
    .expect("parse render-feature-expert-svg");
    match render {
        ShellCommand::RenderFeatureExpertSvg {
            seq_id,
            target,
            output,
        } => {
            assert_eq!(seq_id, "s");
            assert_eq!(output, "out.svg");
            assert_eq!(
                target,
                FeatureExpertTarget::RestrictionSite {
                    cut_pos_1based: 123,
                    enzyme: Some("EcoRI".to_string()),
                    recognition_start_1based: Some(100),
                    recognition_end_1based: Some(106),
                }
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let render_protein_compare = parse_shell_line(
        "render-feature-expert-svg s protein-comparison --transcript TX_TPROT out.svg",
    )
    .expect("parse render-feature-expert-svg transcript protein comparison");
    match render_protein_compare {
        ShellCommand::RenderFeatureExpertSvg {
            seq_id,
            target,
            output,
        } => {
            assert_eq!(seq_id, "s");
            assert_eq!(output, "out.svg");
            assert_eq!(
                target,
                FeatureExpertTarget::ProteinComparison {
                    transcript_id_filter: Some("TX_TPROT".to_string()),
                    protein_feature_filter: Default::default(),
                    external_source: None,
                    external_entry_id: None,
                }
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_panels_commands() {
    let import = parse_shell_line(
            "panels import-isoform seq_a assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1 --strict",
        )
        .expect("parse panels import-isoform");
    match import {
        ShellCommand::PanelsImportIsoform {
            seq_id,
            panel_path,
            panel_id,
            strict,
        } => {
            assert_eq!(seq_id, "seq_a");
            assert_eq!(panel_path, "assets/panels/tp53_isoforms_v1.json");
            assert_eq!(panel_id.as_deref(), Some("tp53_isoforms_v1"));
            assert!(strict);
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let inspect = parse_shell_line("panels inspect-isoform seq_a tp53_isoforms_v1")
        .expect("parse panels inspect-isoform");
    assert!(matches!(inspect, ShellCommand::PanelsInspectIsoform { .. }));

    let render =
        parse_shell_line("panels render-isoform-svg seq_a tp53_isoforms_v1 tp53.panel.svg")
            .expect("parse panels render-isoform-svg");
    assert!(matches!(
        render,
        ShellCommand::PanelsRenderIsoformSvg { .. }
    ));

    let validate = parse_shell_line(
        "panels validate-isoform assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1",
    )
    .expect("parse panels validate-isoform");
    assert!(matches!(
        validate,
        ShellCommand::PanelsValidateIsoform { .. }
    ));
}

#[test]
fn parse_uniprot_commands() {
    let genbank = parse_shell_line("genbank fetch NC_000001 --as-id tp73_fetch")
        .expect("parse genbank fetch");
    match genbank {
        ShellCommand::GenbankFetch { accession, as_id } => {
            assert_eq!(accession, "NC_000001");
            assert_eq!(as_id.as_deref(), Some("tp73_fetch"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let fetch = parse_shell_line("uniprot fetch P04637 --entry-id TP53_UNIPROT")
        .expect("parse uniprot fetch");
    match fetch {
        ShellCommand::UniprotFetch { query, entry_id } => {
            assert_eq!(query, "P04637");
            assert_eq!(entry_id.as_deref(), Some("TP53_UNIPROT"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let import = parse_shell_line("uniprot import-swissprot tp53.txt --entry-id TP53_FILE")
        .expect("parse uniprot import-swissprot");
    match import {
        ShellCommand::UniprotImportSwissProt { path, entry_id } => {
            assert_eq!(path, "tp53.txt");
            assert_eq!(entry_id.as_deref(), Some("TP53_FILE"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let map = parse_shell_line(
        "uniprot map P04637 seq_a --projection-id tp53_map --transcript ENST00000269305.9",
    )
    .expect("parse uniprot map");
    match map {
        ShellCommand::UniprotMap {
            entry_id,
            seq_id,
            projection_id,
            transcript_id,
        } => {
            assert_eq!(entry_id, "P04637");
            assert_eq!(seq_id, "seq_a");
            assert_eq!(projection_id.as_deref(), Some("tp53_map"));
            assert_eq!(transcript_id.as_deref(), Some("ENST00000269305.9"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let projections =
        parse_shell_line("uniprot projection-list --seq seq_a").expect("parse projection-list");
    assert!(matches!(
        projections,
        ShellCommand::UniprotProjectionList { .. }
    ));

    let feature = parse_shell_line(
        "uniprot feature-coding-dna tp53_map DNA-binding --transcript ENST00000269305.9 --mode both --speed-profile human",
    )
    .expect("parse feature-coding-dna");
    match feature {
        ShellCommand::UniprotFeatureCodingDna {
            projection_id,
            feature_query,
            transcript_id,
            mode,
            translation_speed_profile,
        } => {
            assert_eq!(projection_id, "tp53_map");
            assert_eq!(feature_query, "DNA-binding");
            assert_eq!(transcript_id.as_deref(), Some("ENST00000269305.9"));
            assert_eq!(mode, UniprotFeatureCodingDnaQueryMode::Both);
            assert_eq!(
                translation_speed_profile,
                Some(TranslationSpeedProfile::Human)
            );
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let resolve =
        parse_shell_line("uniprot resolve-ensembl-links tp53_map --transcript ENST00000269305.9")
            .expect("parse resolve-ensembl-links");
    assert!(matches!(
        resolve,
        ShellCommand::UniprotResolveEnsemblLinks { .. }
    ));

    let audit = parse_shell_line(
        "uniprot audit-projection tp53_map --transcript ENST00000269305.9 --ensembl-entry ENSP00000269305 --report-id tp53_audit",
    )
    .expect("parse audit-projection");
    match audit {
        ShellCommand::UniprotAuditProjection {
            projection_id,
            transcript_id,
            ensembl_entry_id,
            report_id,
        } => {
            assert_eq!(projection_id, "tp53_map");
            assert_eq!(transcript_id.as_deref(), Some("ENST00000269305.9"));
            assert_eq!(ensembl_entry_id.as_deref(), Some("ENSP00000269305"));
            assert_eq!(report_id.as_deref(), Some("tp53_audit"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let parity_export = parse_shell_line(
        "uniprot audit-parity-export tp53_audit_parity /tmp/tp53_audit_parity.json",
    )
    .expect("parse audit-parity-export");
    assert!(matches!(
        parity_export,
        ShellCommand::UniprotAuditParityExport { .. }
    ));
}

#[test]
fn parse_ensembl_protein_commands() {
    let fetch = parse_shell_line("ensembl-protein fetch ENSP00000288602 --entry-id BRAF_ENS")
        .expect("parse ensembl-protein fetch");
    match fetch {
        ShellCommand::EnsemblProteinFetch { query, entry_id } => {
            assert_eq!(query, "ENSP00000288602");
            assert_eq!(entry_id.as_deref(), Some("BRAF_ENS"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let show = parse_shell_line("ensembl-protein show ENSP00000288602")
        .expect("parse ensembl-protein show");
    assert!(matches!(show, ShellCommand::EnsemblProteinShow { .. }));

    let list = parse_shell_line("ensembl-protein list").expect("parse ensembl-protein list");
    assert!(matches!(list, ShellCommand::EnsemblProteinList));

    let import = parse_shell_line(
        "ensembl-protein import-sequence ENSP00000288602 --output-id braf_ens_protein",
    )
    .expect("parse ensembl-protein import-sequence");
    match import {
        ShellCommand::EnsemblProteinImportSequence {
            entry_id,
            output_id,
        } => {
            assert_eq!(entry_id, "ENSP00000288602");
            assert_eq!(output_id.as_deref(), Some("braf_ens_protein"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

#[test]
fn parse_ensembl_gene_commands() {
    let fetch =
        parse_shell_line("ensembl-gene fetch TP53 --species homo_sapiens --entry-id tp53_gene")
            .expect("parse ensembl-gene fetch");
    match fetch {
        ShellCommand::EnsemblGeneFetch {
            query,
            species,
            entry_id,
        } => {
            assert_eq!(query, "TP53");
            assert_eq!(species.as_deref(), Some("homo_sapiens"));
            assert_eq!(entry_id.as_deref(), Some("tp53_gene"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    let show = parse_shell_line("ensembl-gene show TP53_GENE").expect("parse ensembl-gene show");
    assert!(matches!(show, ShellCommand::EnsemblGeneShow { .. }));

    let list = parse_shell_line("ensembl-gene list").expect("parse ensembl-gene list");
    assert!(matches!(list, ShellCommand::EnsemblGeneList));

    let import =
        parse_shell_line("ensembl-gene import-sequence TP53_GENE --output-id tp53_gene_seq")
            .expect("parse ensembl-gene import-sequence");
    match import {
        ShellCommand::EnsemblGeneImportSequence {
            entry_id,
            output_id,
        } => {
            assert_eq!(entry_id, "TP53_GENE");
            assert_eq!(output_id.as_deref(), Some("tp53_gene_seq"));
        }
        other => panic!("unexpected command: {other:?}"),
    }
}

fn tp53_isoform_test_sequence() -> DNAsequence {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(800)).expect("valid dna");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::simple_range(99, 560),
        qualifiers: vec![
            ("gene".into(), Some("TP53".to_string())),
            (
                "transcript_id".into(),
                Some("ENST00000269305.9".to_string()),
            ),
            ("label".into(), Some("TP53-201".to_string())),
        ]
        .into_iter()
        .collect(),
    });
    dna
}

fn uniprot_projection_test_sequence() -> DNAsequence {
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(300)).expect("valid dna");
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::simple_range(99, 360),
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
    dna
}

fn synthetic_ensembl_entry() -> EnsemblProteinEntry {
    EnsemblProteinEntry {
        schema: "gentle.ensembl_protein_entry.v1".to_string(),
        entry_id: "ENSPTOY1".to_string(),
        protein_id: "ENSPTOY1".to_string(),
        protein_version: Some(1),
        transcript_id: "TX_TPROT".to_string(),
        transcript_version: Some(1),
        gene_id: Some("ENSGTOY1".to_string()),
        gene_symbol: Some("TPROT".to_string()),
        transcript_display_name: Some("TPROT-201".to_string()),
        species: Some("homo_sapiens".to_string()),
        transcript_exons: vec![EnsemblTranscriptExon {
            exon_id: "TX_TPROT_exon_1".to_string(),
            start_1based: 100,
            end_1based: 360,
            strand: 1,
            seq_region_name: Some("1".to_string()),
            version: Some(1),
        }],
        transcript_translation: Some(EnsemblTranscriptTranslation {
            protein_id: "ENSPTOY1".to_string(),
            protein_version: Some(1),
            length_aa: Some(60),
            genomic_start_1based: Some(100),
            genomic_end_1based: Some(279),
            species: Some("homo_sapiens".to_string()),
        }),
        sequence: "M".repeat(60),
        sequence_length: 60,
        features: vec![EnsemblProteinFeature {
            feature_key: "PF02196".to_string(),
            feature_type: "Pfam".to_string(),
            description: Some("synthetic Ensembl domain".to_string()),
            start_aa: Some(10),
            end_aa: Some(20),
            interpro_id: Some("IPR003116".to_string()),
            qualifiers: HashMap::from([
                ("feature_id".to_string(), "PF02196".to_string()),
                ("hit_name".to_string(), "PF02196".to_string()),
            ])
            .into_iter()
            .collect(),
        }],
        aliases: vec![],
        source: "ensembl_rest".to_string(),
        source_query: Some("ENSPTOY1".to_string()),
        imported_at_unix_ms: 0,
        transcript_lookup_source_url: "https://rest.ensembl.org/lookup/id/TX_TPROT?content-type=application/json;expand=1".to_string(),
        protein_lookup_source_url: Some(
            "https://rest.ensembl.org/lookup/id/ENSPTOY1?content-type=application/json"
                .to_string(),
        ),
        sequence_source_url:
            "https://rest.ensembl.org/sequence/id/ENSPTOY1?type=protein;content-type=application/json"
                .to_string(),
        feature_source_url:
            "https://rest.ensembl.org/overlap/translation/ENSPTOY1?feature=protein_feature;content-type=application/json"
                .to_string(),
        raw_transcript_lookup_json: "{}".to_string(),
        raw_protein_lookup_json: Some("{}".to_string()),
        raw_sequence_json: "{}".to_string(),
        raw_feature_json: "[]".to_string(),
    }
}

fn synthetic_ensembl_gene_entry() -> EnsemblGeneEntry {
    EnsemblGeneEntry {
        schema: "gentle.ensembl_gene_entry.v1".to_string(),
        entry_id: "TP53_GENE".to_string(),
        gene_id: "ENSG00000141510".to_string(),
        gene_version: Some(19),
        gene_symbol: Some("TP53".to_string()),
        gene_display_name: Some("tumor protein p53".to_string()),
        species: Some("homo_sapiens".to_string()),
        assembly_name: Some("GRCh38".to_string()),
        biotype: Some("protein_coding".to_string()),
        strand: Some(-1),
        seq_region_name: Some("17".to_string()),
        genomic_start_1based: Some(7668402),
        genomic_end_1based: Some(7668513),
        sequence: "ACGT".repeat(28),
        sequence_length: 112,
        transcripts: vec![EnsemblGeneTranscriptSummary {
            transcript_id: "ENST00000269305".to_string(),
            transcript_version: Some(9),
            display_name: Some("TP53-201".to_string()),
            biotype: Some("protein_coding".to_string()),
            start_1based: Some(7668402),
            end_1based: Some(7668513),
            strand: Some(-1),
        }],
        aliases: vec![],
        source: "ensembl_rest".to_string(),
        source_query: Some("TP53".to_string()),
        imported_at_unix_ms: 0,
        lookup_source_url:
            "https://rest.ensembl.org/lookup/symbol/homo_sapiens/TP53?content-type=application/json;expand=1"
                .to_string(),
        sequence_source_url:
            "https://rest.ensembl.org/sequence/id/ENSG00000141510?content-type=application/json;type=genomic"
                .to_string(),
        raw_lookup_json: "{\"id\":\"ENSG00000141510\"}".to_string(),
        raw_sequence_json: "{\"id\":\"ENSG00000141510\",\"seq\":\"ACGT\"}".to_string(),
    }
}

fn synthetic_audit_ensembl_entry() -> EnsemblProteinEntry {
    EnsemblProteinEntry {
        entry_id: "ENSPTOY1".to_string(),
        protein_id: "ENSPTOY1".to_string(),
        transcript_id: "TX1".to_string(),
        transcript_display_name: Some("TX1".to_string()),
        gene_id: Some("ENSGTOY1".to_string()),
        gene_symbol: Some("TOY1".to_string()),
        transcript_exons: vec![
            EnsemblTranscriptExon {
                exon_id: "TX1_exon_1".to_string(),
                start_1based: 100,
                end_1based: 180,
                strand: 1,
                seq_region_name: Some("1".to_string()),
                version: Some(1),
            },
            EnsemblTranscriptExon {
                exon_id: "TX1_exon_2".to_string(),
                start_1based: 300,
                end_1based: 360,
                strand: 1,
                seq_region_name: Some("1".to_string()),
                version: Some(1),
            },
        ],
        transcript_translation: Some(EnsemblTranscriptTranslation {
            protein_id: "ENSPTOY1".to_string(),
            protein_version: Some(1),
            length_aa: Some(47),
            genomic_start_1based: Some(100),
            genomic_end_1based: Some(360),
            species: Some("homo_sapiens".to_string()),
        }),
        sequence: "M".repeat(47),
        sequence_length: 47,
        ..synthetic_ensembl_entry()
    }
}

#[test]
fn execute_inspect_feature_expert_returns_view_json() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "s".to_string(),
        DNAsequence::from_sequence("TTTACGTAAACGTGGG").expect("valid dna"),
    );
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
        .expect("annotate tfbs");
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

    let output = execute_shell_command(
        &mut engine,
        &ShellCommand::InspectFeatureExpert {
            seq_id: "s".to_string(),
            target: FeatureExpertTarget::TfbsFeature { feature_id },
        },
    )
    .expect("execute inspect-feature-expert");
    assert!(!output.state_changed);
    assert_eq!(output.output["kind"].as_str(), Some("tfbs"));
    assert_eq!(
        output.output["data"]["feature_id"].as_u64(),
        Some(feature_id as u64)
    );
}

#[test]
fn execute_panels_import_and_inspect_isoform_architecture() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), tp53_isoform_test_sequence());
    let mut engine = GentleEngine::from_state(state);

    let import = execute_shell_command(
        &mut engine,
        &ShellCommand::PanelsImportIsoform {
            seq_id: "seq_a".to_string(),
            panel_path: "assets/panels/tp53_isoforms_v1.json".to_string(),
            panel_id: Some("tp53_isoforms_v1".to_string()),
            strict: false,
        },
    )
    .expect("execute panels import");
    assert!(import.state_changed);

    let inspect = execute_shell_command(
        &mut engine,
        &ShellCommand::PanelsInspectIsoform {
            seq_id: "seq_a".to_string(),
            panel_id: "tp53_isoforms_v1".to_string(),
        },
    )
    .expect("execute panels inspect");
    assert!(!inspect.state_changed);
    assert_eq!(
        inspect.output["kind"].as_str(),
        Some("isoform_architecture")
    );
    assert_eq!(
        inspect.output["data"]["panel_id"].as_str(),
        Some("tp53_isoforms_v1")
    );
}

#[test]
fn execute_uniprot_import_list_and_projection() {
    let td = tempdir().expect("tempdir");
    let swiss_path = td.path().join("toy_uniprot.txt");
    let swiss_text = r#"ID   TOY1_HUMAN              Reviewed;         30 AA.
AC   PTEST1;
DE   RecName: Full=Toy DNA-binding protein;
GN   Name=TOY1;
OS   Homo sapiens (Human).
DR   Ensembl; TX1; ENSPTOY1; ENSGTOY1.
FT   DOMAIN          2..8
FT                   /note="toy domain"
FT   REGION          25..30
FT                   /note="junction feature"
SQ   SEQUENCE   30 AA;  3333 MW;  0000000000000000 CRC64;
     MEEPQSDPSV EPPLSQETFSDLWKLLPENA
//
"#;
    fs::write(&swiss_path, swiss_text).expect("write swiss file");

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_u".to_string(), uniprot_projection_test_sequence());
    let mut engine = GentleEngine::from_state(state);

    let import = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotImportSwissProt {
            path: swiss_path.to_string_lossy().to_string(),
            entry_id: None,
        },
    )
    .expect("execute uniprot import");
    assert!(import.state_changed);

    let listed = execute_shell_command(&mut engine, &ShellCommand::UniprotList)
        .expect("execute uniprot list");
    assert!(!listed.state_changed);
    let rows = listed
        .output
        .as_array()
        .expect("uniprot list output must be array");
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0]["entry_id"].as_str(), Some("PTEST1"));

    let mapped = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotMap {
            entry_id: "PTEST1".to_string(),
            seq_id: "seq_u".to_string(),
            projection_id: None,
            transcript_id: None,
        },
    )
    .expect("execute uniprot map");
    assert!(mapped.state_changed);

    let projection_id = "PTEST1@seq_u".to_string();
    let projection = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotProjectionShow {
            projection_id: projection_id.clone(),
        },
    )
    .expect("execute uniprot projection-show");
    assert!(!projection.state_changed);
    assert_eq!(
        projection.output["projection_id"].as_str(),
        Some(projection_id.as_str())
    );
    assert_eq!(projection.output["entry_id"].as_str(), Some("PTEST1"));
    assert_eq!(projection.output["seq_id"].as_str(), Some("seq_u"));
    assert!(
        projection.output["transcript_projections"]
            .as_array()
            .map(|v| !v.is_empty())
            .unwrap_or(false)
    );

    let inspect = execute_shell_command(
        &mut engine,
        &ShellCommand::InspectFeatureExpert {
            seq_id: "seq_u".to_string(),
            target: FeatureExpertTarget::UniprotProjection {
                projection_id: projection_id.clone(),
                protein_feature_filter: Default::default(),
            },
        },
    )
    .expect("inspect uniprot projection expert");
    assert!(!inspect.state_changed);
    assert_eq!(
        inspect.output["kind"].as_str(),
        Some("isoform_architecture")
    );
    assert_eq!(
        inspect.output["data"]["panel_id"].as_str(),
        Some(projection_id.as_str())
    );
    assert_eq!(inspect.output["data"]["gene_symbol"].as_str(), Some("TOY1"));
    assert_eq!(
        inspect.output["data"]["panel_source"].as_str(),
        Some("Transcript-native proteins with optional UniProt opinion PTEST1 (PTEST1)")
    );
    assert!(
        inspect.output["data"]["protein_lanes"]
            .as_array()
            .map(|rows| !rows.is_empty())
            .unwrap_or(false)
    );

    let feature_query = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotFeatureCodingDna {
            projection_id: projection_id.clone(),
            feature_query: "junction".to_string(),
            transcript_id: Some("TX1".to_string()),
            mode: UniprotFeatureCodingDnaQueryMode::Both,
            translation_speed_profile: Some(TranslationSpeedProfile::Human),
        },
    )
    .expect("query feature coding dna");
    assert!(!feature_query.state_changed);
    assert_eq!(
        feature_query.output["projection_id"].as_str(),
        Some(projection_id.as_str())
    );
    assert_eq!(feature_query.output["match_count"].as_u64(), Some(1));
    assert_eq!(
        feature_query.output["matches"][0]["primary_exon_pair"]["from_exon_ordinal"].as_u64(),
        Some(1)
    );
    assert_eq!(
        feature_query.output["matches"][0]["primary_exon_pair"]["to_exon_ordinal"].as_u64(),
        Some(2)
    );
    assert_eq!(
        feature_query.output["matches"][0]["translation_speed_optimized_dna"]
            .as_str()
            .map(str::len),
        Some(18)
    );
}

#[test]
fn execute_uniprot_projection_audit_commands() {
    let td = tempdir().expect("tempdir");
    let swiss_path = td.path().join("toy_uniprot_audit.txt");
    let export_path = td.path().join("toy_uniprot_audit.json");
    let parity_export_path = td.path().join("toy_uniprot_audit_parity.json");
    let swiss_text = r#"ID   TOY1_HUMAN              Reviewed;         30 AA.
AC   PTEST1;
DE   RecName: Full=Toy DNA-binding protein;
GN   Name=TOY1;
OS   Homo sapiens (Human).
DR   Ensembl; TX1; ENSPTOY1; ENSGTOY1.
FT   VARIANT         2
FT                   /note="synthetic variant"
FT   CONFLICT        10
FT                   /note="synthetic conflict"
SQ   SEQUENCE   30 AA;  3333 MW;  0000000000000000 CRC64;
     MEEPQSDPSV EPPLSQETFSDLWKLLPENA
//
"#;
    fs::write(&swiss_path, swiss_text).expect("write swiss audit file");

    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_u".to_string(), uniprot_projection_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::ImportUniprotSwissProt {
            path: swiss_path.to_string_lossy().to_string(),
            entry_id: None,
        })
        .expect("import swiss");
    engine
        .apply(Operation::ProjectUniprotToGenome {
            seq_id: "seq_u".to_string(),
            entry_id: "PTEST1".to_string(),
            projection_id: None,
            transcript_id: None,
        })
        .expect("map projection");
    engine.state_mut().metadata.insert(
        "ensembl_protein_entries".to_string(),
        serde_json::json!({
            "schema": "gentle.ensembl_protein_entries.v1",
            "updated_at_unix_ms": 1u128,
            "entries": {
                "ENSPTOY1": synthetic_audit_ensembl_entry()
            }
        }),
    );

    let resolved = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotResolveEnsemblLinks {
            projection_id: "PTEST1@seq_u".to_string(),
            transcript_id: Some("TX1".to_string()),
        },
    )
    .expect("resolve links");
    assert!(!resolved.state_changed);
    assert_eq!(
        resolved.output["schema"].as_str(),
        Some("gentle.uniprot_ensembl_links.v1")
    );
    assert_eq!(resolved.output["rows"].as_array().map(Vec::len), Some(1));

    let accounting = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotTranscriptAccounting {
            projection_id: "PTEST1@seq_u".to_string(),
            transcript_id: Some("TX1".to_string()),
        },
    )
    .expect("transcript accounting");
    assert_eq!(
        accounting.output["schema"].as_str(),
        Some("gentle.uniprot_projection_transcript_accounting.v1")
    );

    let exon_compare = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotCompareEnsemblExons {
            projection_id: "PTEST1@seq_u".to_string(),
            transcript_id: Some("TX1".to_string()),
            ensembl_entry_id: Some("ENSPTOY1".to_string()),
        },
    )
    .expect("compare exons");
    assert_eq!(
        exon_compare.output["schema"].as_str(),
        Some("gentle.uniprot_ensembl_exon_compare.v1")
    );

    let peptide_compare = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotCompareEnsemblPeptide {
            projection_id: "PTEST1@seq_u".to_string(),
            transcript_id: Some("TX1".to_string()),
            ensembl_entry_id: Some("ENSPTOY1".to_string()),
        },
    )
    .expect("compare peptide");
    assert_eq!(
        peptide_compare.output["schema"].as_str(),
        Some("gentle.uniprot_ensembl_peptide_compare.v1")
    );

    let audit = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotAuditProjection {
            projection_id: "PTEST1@seq_u".to_string(),
            transcript_id: Some("TX1".to_string()),
            ensembl_entry_id: Some("ENSPTOY1".to_string()),
            report_id: Some("toy_audit".to_string()),
        },
    )
    .expect("audit projection");
    assert!(audit.state_changed);
    assert_eq!(
        audit.output["report"]["schema"].as_str(),
        Some("gentle.uniprot_projection_audit.v1")
    );
    assert!(
        audit.output["report"]["maintainer_email_draft"]["body"]
            .as_str()
            .map(|body| body.contains("presumed inconsistency"))
            .unwrap_or(false)
    );

    let listed = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotAuditList {
            seq_id: Some("seq_u".to_string()),
        },
    )
    .expect("list audits");
    assert_eq!(listed.output.as_array().map(Vec::len), Some(1));

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotAuditShow {
            report_id: "toy_audit".to_string(),
        },
    )
    .expect("show audit");
    assert_eq!(shown.output["report_id"].as_str(), Some("toy_audit"));

    let exported = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotAuditExport {
            report_id: "toy_audit".to_string(),
            path: export_path.display().to_string(),
        },
    )
    .expect("export audit");
    assert_eq!(
        exported.output["schema"].as_str(),
        Some("gentle.uniprot_projection_audit_export.v1")
    );
    assert!(export_path.exists());

    let parity = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotAuditParity {
            projection_id: "PTEST1@seq_u".to_string(),
            transcript_id: Some("TX1".to_string()),
            ensembl_entry_id: Some("ENSPTOY1".to_string()),
            report_id: Some("toy_audit_parity".to_string()),
        },
    )
    .expect("audit parity");
    assert!(parity.state_changed);
    assert_eq!(
        parity.output["report"]["schema"].as_str(),
        Some("gentle.uniprot_projection_audit_parity.v1")
    );

    let parity_list = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotAuditParityList {
            seq_id: Some("seq_u".to_string()),
        },
    )
    .expect("list audit parity");
    assert_eq!(parity_list.output.as_array().map(Vec::len), Some(1));

    let parity_show = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotAuditParityShow {
            report_id: "toy_audit_parity".to_string(),
        },
    )
    .expect("show audit parity");
    assert_eq!(
        parity_show.output["report_id"].as_str(),
        Some("toy_audit_parity")
    );

    let parity_export = execute_shell_command(
        &mut engine,
        &ShellCommand::UniprotAuditParityExport {
            report_id: "toy_audit_parity".to_string(),
            path: parity_export_path.display().to_string(),
        },
    )
    .expect("export audit parity");
    assert_eq!(
        parity_export.output["schema"].as_str(),
        Some("gentle.uniprot_projection_audit_parity_export.v1")
    );
    assert!(parity_export_path.exists());
}

#[test]
fn execute_transcript_protein_comparison_expert_without_uniprot() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence(&"ATGAAAACC".repeat(20)).expect("valid DNA");
    let dna_len_i64 = i64::try_from(dna.len()).unwrap();
    dna.features_mut().push(Feature {
        kind: "source".into(),
        location: Location::simple_range(0, dna_len_i64),
        qualifiers: vec![("organism".into(), Some("Homo sapiens".to_string()))]
            .into_iter()
            .collect(),
    });
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::simple_range(0, 180),
        qualifiers: vec![
            ("gene".into(), Some("TPROT".to_string())),
            ("transcript_id".into(), Some("TX_TPROT".to_string())),
            ("label".into(), Some("TX_TPROT".to_string())),
            ("cds_ranges_1based".into(), Some("1-180".to_string())),
        ]
        .into_iter()
        .collect(),
    });
    state.sequences.insert("toy_tx_only".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let inspect = execute_shell_command(
        &mut engine,
        &ShellCommand::InspectFeatureExpert {
            seq_id: "toy_tx_only".to_string(),
            target: FeatureExpertTarget::ProteinComparison {
                transcript_id_filter: Some("TX_TPROT".to_string()),
                protein_feature_filter: Default::default(),
                external_source: None,
                external_entry_id: None,
            },
        },
    )
    .expect("inspect transcript-native protein comparison expert");
    assert!(!inspect.state_changed);
    assert_eq!(
        inspect.output["kind"].as_str(),
        Some("isoform_architecture")
    );
    assert_eq!(
        inspect.output["data"]["panel_id"].as_str(),
        Some("protein_compare:TX_TPROT@toy_tx_only")
    );
    assert_eq!(
        inspect.output["data"]["panel_source"].as_str(),
        Some("Transcript-native protein derivation (external protein opinions optional)")
    );
    let protein_lanes = inspect.output["data"]["protein_lanes"]
        .as_array()
        .expect("protein lanes array");
    assert_eq!(protein_lanes.len(), 1);
    assert_eq!(
        protein_lanes[0]["comparison"]["status"].as_str(),
        Some("derived_only")
    );
    assert!(protein_lanes[0]["comparison"]["external_opinion"].is_null());
}

#[test]
fn execute_ensembl_protein_list_show_import_and_compare() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence(&"ATGAAAACC".repeat(20)).expect("valid DNA");
    let dna_len_i64 = i64::try_from(dna.len()).unwrap();
    dna.features_mut().push(Feature {
        kind: "source".into(),
        location: Location::simple_range(0, dna_len_i64),
        qualifiers: vec![("organism".into(), Some("Homo sapiens".to_string()))]
            .into_iter()
            .collect(),
    });
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::simple_range(0, 180),
        qualifiers: vec![
            ("gene".into(), Some("TPROT".to_string())),
            ("transcript_id".into(), Some("TX_TPROT".to_string())),
            ("label".into(), Some("TX_TPROT".to_string())),
            ("cds_ranges_1based".into(), Some("1-180".to_string())),
        ]
        .into_iter()
        .collect(),
    });
    state.sequences.insert("toy_tx_only".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);
    engine
        .upsert_ensembl_protein_entry(synthetic_ensembl_entry())
        .expect("upsert Ensembl protein entry");

    let listed = execute_shell_command(&mut engine, &ShellCommand::EnsemblProteinList)
        .expect("execute ensembl-protein list");
    assert!(!listed.state_changed);
    assert_eq!(listed.output[0]["entry_id"].as_str(), Some("ENSPTOY1"));

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::EnsemblProteinShow {
            entry_id: "ENSPTOY1".to_string(),
        },
    )
    .expect("execute ensembl-protein show");
    assert!(!shown.state_changed);
    assert_eq!(shown.output["transcript_id"].as_str(), Some("TX_TPROT"));

    let imported = execute_shell_command(
        &mut engine,
        &ShellCommand::EnsemblProteinImportSequence {
            entry_id: "ENSPTOY1".to_string(),
            output_id: Some("toy_ens_import".to_string()),
        },
    )
    .expect("execute ensembl-protein import-sequence");
    assert!(imported.state_changed);
    assert!(
        engine.state().sequences.contains_key("toy_ens_import"),
        "imported Ensembl protein sequence should exist in state"
    );

    let inspect = execute_shell_command(
        &mut engine,
        &ShellCommand::InspectFeatureExpert {
            seq_id: "toy_tx_only".to_string(),
            target: FeatureExpertTarget::ProteinComparison {
                transcript_id_filter: Some("TX_TPROT".to_string()),
                protein_feature_filter: ProteinFeatureFilter {
                    include_feature_keys: vec!["PF02196".to_string()],
                    exclude_feature_keys: vec![],
                },
                external_source: Some(ProteinExternalOpinionSource::Ensembl),
                external_entry_id: Some("ENSPTOY1".to_string()),
            },
        },
    )
    .expect("inspect Ensembl-backed transcript protein comparison expert");
    assert!(!inspect.state_changed);
    assert_eq!(
        inspect.output["data"]["panel_id"].as_str(),
        Some("ENSPTOY1")
    );
    assert_eq!(
        inspect.output["data"]["protein_lanes"][0]["comparison"]["external_opinion"]["source"]
            .as_str(),
        Some("ensembl")
    );
    assert!(
        inspect.output["data"]["protein_lanes"][0]["domains"]
            .as_array()
            .map(|rows| rows.iter().any(|row| row["name"]
                .as_str()
                .map(|name| name.contains("synthetic Ensembl domain"))
                .unwrap_or(false)))
            .unwrap_or(false)
    );
}

#[test]
fn execute_ensembl_gene_list_show_and_import() {
    let mut engine = GentleEngine::default();
    engine
        .upsert_ensembl_gene_entry(synthetic_ensembl_gene_entry())
        .expect("upsert Ensembl gene entry");

    let listed = execute_shell_command(&mut engine, &ShellCommand::EnsemblGeneList)
        .expect("execute ensembl-gene list");
    assert!(!listed.state_changed);
    assert_eq!(listed.output[0]["entry_id"].as_str(), Some("tp53_gene"));

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::EnsemblGeneShow {
            entry_id: "TP53_GENE".to_string(),
        },
    )
    .expect("execute ensembl-gene show");
    assert!(!shown.state_changed);
    assert_eq!(shown.output["gene_id"].as_str(), Some("ENSG00000141510"));

    let imported = execute_shell_command(
        &mut engine,
        &ShellCommand::EnsemblGeneImportSequence {
            entry_id: "TP53_GENE".to_string(),
            output_id: Some("tp53_gene_seq".to_string()),
        },
    )
    .expect("execute ensembl-gene import-sequence");
    assert!(imported.state_changed);
    assert!(
        engine.state().sequences.contains_key("tp53_gene_seq"),
        "imported Ensembl gene sequence should exist in state"
    );
}

#[test]
fn execute_transcript_protein_comparison_expert_svg_shell_and_operation_routes_match() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence(&"ATGAAAACC".repeat(20)).expect("valid DNA");
    let dna_len_i64 = i64::try_from(dna.len()).unwrap();
    dna.features_mut().push(Feature {
        kind: "source".into(),
        location: Location::simple_range(0, dna_len_i64),
        qualifiers: vec![("organism".into(), Some("Homo sapiens".to_string()))]
            .into_iter()
            .collect(),
    });
    dna.features_mut().push(Feature {
        kind: "mRNA".into(),
        location: Location::simple_range(0, 180),
        qualifiers: vec![
            ("gene".into(), Some("TPROT".to_string())),
            ("transcript_id".into(), Some("TX_TPROT".to_string())),
            ("label".into(), Some("TX_TPROT".to_string())),
            ("cds_ranges_1based".into(), Some("1-180".to_string())),
        ]
        .into_iter()
        .collect(),
    });
    state.sequences.insert("toy_tx_only".to_string(), dna);
    let mut engine = GentleEngine::from_state(state);

    let tmp = tempdir().expect("temp dir");
    let op_svg = tmp.path().join("protein_compare.op.svg");
    let shell_svg = tmp.path().join("protein_compare.shell.svg");
    let op_path = op_svg.display().to_string();
    let shell_path = shell_svg.display().to_string();

    engine
        .apply(Operation::RenderFeatureExpertSvg {
            seq_id: "toy_tx_only".to_string(),
            target: FeatureExpertTarget::ProteinComparison {
                transcript_id_filter: Some("TX_TPROT".to_string()),
                protein_feature_filter: Default::default(),
                external_source: None,
                external_entry_id: None,
            },
            path: op_path.clone(),
        })
        .expect("render transcript protein comparison op route");

    execute_shell_command(
        &mut engine,
        &ShellCommand::RenderFeatureExpertSvg {
            seq_id: "toy_tx_only".to_string(),
            target: FeatureExpertTarget::ProteinComparison {
                transcript_id_filter: Some("TX_TPROT".to_string()),
                protein_feature_filter: Default::default(),
                external_source: None,
                external_entry_id: None,
            },
            output: shell_path.clone(),
        },
    )
    .expect("render transcript protein comparison shell route");

    let op_text = fs::read_to_string(op_path).expect("read op svg");
    let shell_text = fs::read_to_string(shell_path).expect("read shell svg");
    assert_eq!(
        op_text, shell_text,
        "RenderFeatureExpertSvg operation and shell routes must match for transcript protein comparison"
    );
    assert!(op_text.contains("protein-rail"));
}

#[test]
fn execute_panels_validate_isoform_returns_report() {
    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::PanelsValidateIsoform {
            panel_path: "assets/panels/tp53_isoforms_v1.json".to_string(),
            panel_id: Some("tp53_isoforms_v1".to_string()),
        },
    )
    .expect("validate panel");
    assert!(!out.state_changed);
    assert_eq!(
        out.output["schema"].as_str(),
        Some("gentle.isoform_panel_validation_report.v1")
    );
    assert_eq!(out.output["panel_id"].as_str(), Some("tp53_isoforms_v1"));
    assert!(matches!(
        out.output["status"].as_str(),
        Some("ok") | Some("warning")
    ));
    assert!(out.output["isoform_count"].as_u64().unwrap_or(0) >= 1);
}

#[test]
fn execute_isoform_svg_routes_are_byte_identical() {
    execute_isoform_svg_routes_are_byte_identical_impl();
}

#[inline(never)]
fn execute_isoform_svg_routes_are_byte_identical_impl() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), tp53_isoform_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::ImportIsoformPanel {
            seq_id: "seq_a".to_string(),
            panel_path: "assets/panels/tp53_isoforms_v1.json".to_string(),
            panel_id: Some("tp53_isoforms_v1".to_string()),
            strict: false,
        })
        .expect("import panel");

    let tmp = tempdir().expect("temp dir");
    let op_svg = tmp.path().join("isoform.op.svg");
    let shell_svg = tmp.path().join("isoform.shell.svg");
    let expert_svg = tmp.path().join("isoform.expert.svg");
    let op_path = op_svg.display().to_string();
    let shell_path = shell_svg.display().to_string();
    let expert_path = expert_svg.display().to_string();

    engine
        .apply(Operation::RenderIsoformArchitectureSvg {
            seq_id: "seq_a".to_string(),
            panel_id: "tp53_isoforms_v1".to_string(),
            path: op_path.clone(),
        })
        .expect("render op route");

    execute_shell_command(
        &mut engine,
        &ShellCommand::PanelsRenderIsoformSvg {
            seq_id: "seq_a".to_string(),
            panel_id: "tp53_isoforms_v1".to_string(),
            output: shell_path.clone(),
        },
    )
    .expect("render shell route");

    execute_shell_command(
        &mut engine,
        &ShellCommand::RenderFeatureExpertSvg {
            seq_id: "seq_a".to_string(),
            target: FeatureExpertTarget::IsoformArchitecture {
                panel_id: "tp53_isoforms_v1".to_string(),
            },
            output: expert_path.clone(),
        },
    )
    .expect("render expert route");

    let op_text = fs::read_to_string(op_path).expect("read op svg");
    let shell_text = fs::read_to_string(shell_path).expect("read shell svg");
    let expert_text = fs::read_to_string(expert_path).expect("read expert svg");
    assert_eq!(
        op_text, shell_text,
        "RenderIsoformArchitectureSvg and panels render-isoform-svg outputs must match"
    );
    assert_eq!(
        op_text, expert_text,
        "RenderIsoformArchitectureSvg and render-feature-expert-svg isoform outputs must match"
    );
}

#[test]
fn execute_splicing_feature_expert_svg_shell_and_operation_routes_match() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), tp53_isoform_test_sequence());
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

    let tmp = tempdir().expect("temp dir");
    let op_svg = tmp.path().join("splicing.op.svg");
    let shell_svg = tmp.path().join("splicing.shell.svg");
    let op_path = op_svg.display().to_string();
    let shell_path = shell_svg.display().to_string();

    engine
        .apply(Operation::RenderFeatureExpertSvg {
            seq_id: "seq_a".to_string(),
            target: FeatureExpertTarget::SplicingFeature {
                feature_id,
                scope: SplicingScopePreset::AllOverlappingBothStrands,
            },
            path: op_path.clone(),
        })
        .expect("render splicing op route");

    execute_shell_command(
        &mut engine,
        &ShellCommand::RenderFeatureExpertSvg {
            seq_id: "seq_a".to_string(),
            target: FeatureExpertTarget::SplicingFeature {
                feature_id,
                scope: SplicingScopePreset::AllOverlappingBothStrands,
            },
            output: shell_path.clone(),
        },
    )
    .expect("render splicing shell route");

    let op_text = fs::read_to_string(op_path).expect("read op svg");
    let shell_text = fs::read_to_string(shell_path).expect("read shell svg");
    assert_eq!(
        op_text, shell_text,
        "RenderFeatureExpertSvg operation and shell routes must match for splicing"
    );
    assert!(op_text.contains("cell color intensity encodes exon support frequency"));
}

#[test]
fn parse_dotplot_and_flex_commands() {
    let dotplot = parse_shell_line(
            "dotplot compute seq_a --start 10 --end 110 --mode self_reverse_complement --word-size 8 --step 4 --max-mismatches 1 --tile-bp 256 --id promoter_dp",
        )
        .expect("parse dotplot compute");
    match dotplot {
        ShellCommand::DotplotCompute {
            seq_id,
            reference_seq_id,
            span_start_0based,
            span_end_0based,
            reference_span_start_0based,
            reference_span_end_0based,
            mode,
            word_size,
            step_bp,
            max_mismatches,
            tile_bp,
            dotplot_id,
        } => {
            assert_eq!(seq_id, "seq_a");
            assert!(reference_seq_id.is_none());
            assert_eq!(span_start_0based, Some(10));
            assert_eq!(span_end_0based, Some(110));
            assert!(reference_span_start_0based.is_none());
            assert!(reference_span_end_0based.is_none());
            assert_eq!(mode, DotplotMode::SelfReverseComplement);
            assert_eq!(word_size, 8);
            assert_eq!(step_bp, 4);
            assert_eq!(max_mismatches, 1);
            assert_eq!(tile_bp, Some(256));
            assert_eq!(dotplot_id.as_deref(), Some("promoter_dp"));
        }
        other => panic!("expected DotplotCompute, got {other:?}"),
    }

    let pair_dotplot = parse_shell_line(
        "dotplot compute seq_a --reference-seq seq_ref --start 10 --end 110 --ref-start 200 --ref-end 500 --mode pair_forward --word-size 10 --step 5 --id pair_dp",
    )
    .expect("parse pair dotplot compute");
    match pair_dotplot {
        ShellCommand::DotplotCompute {
            seq_id,
            reference_seq_id,
            span_start_0based,
            span_end_0based,
            reference_span_start_0based,
            reference_span_end_0based,
            mode,
            word_size,
            step_bp,
            max_mismatches,
            tile_bp,
            dotplot_id,
        } => {
            assert_eq!(seq_id, "seq_a");
            assert_eq!(reference_seq_id.as_deref(), Some("seq_ref"));
            assert_eq!(span_start_0based, Some(10));
            assert_eq!(span_end_0based, Some(110));
            assert_eq!(reference_span_start_0based, Some(200));
            assert_eq!(reference_span_end_0based, Some(500));
            assert_eq!(mode, DotplotMode::PairForward);
            assert_eq!(word_size, 10);
            assert_eq!(step_bp, 5);
            assert_eq!(max_mismatches, 0);
            assert_eq!(tile_bp, None);
            assert_eq!(dotplot_id.as_deref(), Some("pair_dp"));
        }
        other => panic!("expected DotplotCompute, got {other:?}"),
    }

    let overlay_dotplot = parse_shell_line(
        "dotplot overlay-compute locus_a --reference-seq genome_ref --query-spec '{\"seq_id\":\"iso_1\",\"label\":\"Isoform A\",\"mode\":\"pair_forward\",\"color_rgb\":[29,78,216]}' --query-spec '{\"seq_id\":\"iso_2\",\"label\":\"Isoform B\",\"span_start_0based\":5,\"span_end_0based\":17,\"mode\":\"pair_reverse_complement\"}' --ref-start 200 --ref-end 500 --word-size 10 --step 5 --id overlay_dp",
    )
    .expect("parse overlay dotplot compute");
    match overlay_dotplot {
        ShellCommand::DotplotOverlayCompute {
            owner_seq_id,
            reference_seq_id,
            reference_span_start_0based,
            reference_span_end_0based,
            queries,
            word_size,
            step_bp,
            max_mismatches,
            tile_bp,
            dotplot_id,
        } => {
            assert_eq!(owner_seq_id, "locus_a");
            assert_eq!(reference_seq_id.as_deref(), Some("genome_ref"));
            assert_eq!(reference_span_start_0based, Some(200));
            assert_eq!(reference_span_end_0based, Some(500));
            assert_eq!(queries.len(), 2);
            assert_eq!(queries[0].seq_id, "iso_1");
            assert_eq!(queries[0].label, "Isoform A");
            assert_eq!(queries[0].mode, DotplotMode::PairForward);
            assert_eq!(queries[0].color_rgb, Some([29, 78, 216]));
            assert_eq!(queries[1].seq_id, "iso_2");
            assert_eq!(queries[1].label, "Isoform B");
            assert_eq!(queries[1].span_start_0based, Some(5));
            assert_eq!(queries[1].span_end_0based, Some(17));
            assert_eq!(queries[1].mode, DotplotMode::PairReverseComplement);
            assert_eq!(word_size, 10);
            assert_eq!(step_bp, 5);
            assert_eq!(max_mismatches, 0);
            assert_eq!(tile_bp, None);
            assert_eq!(dotplot_id.as_deref(), Some("overlay_dp"));
        }
        other => panic!("expected DotplotOverlayCompute, got {other:?}"),
    }

    let render_dotplot = parse_shell_line(
        "dotplot render-svg seq_a pair_dp /tmp/pair_dp.svg --flex-track flex_1 --display-threshold 0.2 --intensity-gain 1.7 --overlay-x-axis shared_exon_anchor --overlay-anchor-exon 27..34",
    )
    .expect("parse dotplot render-svg");
    match render_dotplot {
        ShellCommand::RenderDotplotSvg {
            seq_id,
            dotplot_id,
            output,
            flex_track_id,
            display_density_threshold,
            display_intensity_gain,
            overlay_x_axis_mode,
            overlay_anchor_exon,
        } => {
            assert_eq!(seq_id, "seq_a");
            assert_eq!(dotplot_id, "pair_dp");
            assert_eq!(output, "/tmp/pair_dp.svg");
            assert_eq!(flex_track_id.as_deref(), Some("flex_1"));
            assert_eq!(display_density_threshold, Some(0.2));
            assert_eq!(display_intensity_gain, Some(1.7));
            assert_eq!(
                overlay_x_axis_mode,
                DotplotOverlayXAxisMode::SharedExonAnchor
            );
            assert_eq!(
                overlay_anchor_exon,
                Some(DotplotOverlayAnchorExonRef {
                    start_1based: 27,
                    end_1based: 34,
                })
            );
        }
        other => panic!("expected RenderDotplotSvg, got {other:?}"),
    }

    let render_dotplot_shared_anchor = parse_shell_line(
        "dotplot render-svg seq_a pair_dp /tmp/pair_dp.anchor.svg --overlay-x-axis shared_exon_anchor",
    )
    .expect("parse dotplot render-svg shared exon anchor");
    match render_dotplot_shared_anchor {
        ShellCommand::RenderDotplotSvg {
            overlay_x_axis_mode,
            ..
        } => {
            assert_eq!(
                overlay_x_axis_mode,
                DotplotOverlayXAxisMode::SharedExonAnchor
            );
        }
        other => panic!("expected RenderDotplotSvg, got {other:?}"),
    }

    let render_dotplot_query_anchor = parse_shell_line(
        "dotplot render-svg seq_a pair_dp /tmp/pair_dp.query-anchor.svg --overlay-x-axis query_anchor_bp",
    )
    .expect("parse dotplot render-svg query anchor");
    match render_dotplot_query_anchor {
        ShellCommand::RenderDotplotSvg {
            overlay_x_axis_mode,
            ..
        } => {
            assert_eq!(overlay_x_axis_mode, DotplotOverlayXAxisMode::QueryAnchorBp);
        }
        other => panic!("expected RenderDotplotSvg, got {other:?}"),
    }

    let flex = parse_shell_line(
            "flex compute seq_a --start 25 --end 325 --model at_skew --bin-bp 20 --smoothing-bp 60 --id promoter_flex",
        )
        .expect("parse flex compute");
    match flex {
        ShellCommand::FlexCompute {
            seq_id,
            span_start_0based,
            span_end_0based,
            model,
            bin_bp,
            smoothing_bp,
            track_id,
        } => {
            assert_eq!(seq_id, "seq_a");
            assert_eq!(span_start_0based, Some(25));
            assert_eq!(span_end_0based, Some(325));
            assert_eq!(model, FlexibilityModel::AtSkew);
            assert_eq!(bin_bp, 20);
            assert_eq!(smoothing_bp, Some(60));
            assert_eq!(track_id.as_deref(), Some("promoter_flex"));
        }
        other => panic!("expected FlexCompute, got {other:?}"),
    }

    let splicing_refs = parse_shell_line(
        "splicing-refs derive seq_a 10 220 --seed-feature-id 7 --scope target_group_any_strand --output-prefix tp53_refs",
    )
    .expect("parse splicing-refs derive");
    match splicing_refs {
        ShellCommand::SplicingRefsDerive {
            seq_id,
            span_start_0based,
            span_end_0based,
            seed_feature_id,
            scope,
            output_prefix,
        } => {
            assert_eq!(seq_id, "seq_a");
            assert_eq!(span_start_0based, 10);
            assert_eq!(span_end_0based, 220);
            assert_eq!(seed_feature_id, Some(7));
            assert_eq!(scope, SplicingScopePreset::TargetGroupAnyStrand);
            assert_eq!(output_prefix.as_deref(), Some("tp53_refs"));
        }
        other => panic!("expected SplicingRefsDerive, got {other:?}"),
    }

    let align = parse_shell_line(
        "align compute query target --query-start 5 --query-end 105 --target-start 10 --target-end 120 --mode local --match 3 --mismatch -4 --gap-open -6 --gap-extend -2",
    )
    .expect("parse align compute");
    match align {
        ShellCommand::AlignCompute {
            query_seq_id,
            target_seq_id,
            query_span_start_0based,
            query_span_end_0based,
            target_span_start_0based,
            target_span_end_0based,
            mode,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        } => {
            assert_eq!(query_seq_id, "query");
            assert_eq!(target_seq_id, "target");
            assert_eq!(query_span_start_0based, Some(5));
            assert_eq!(query_span_end_0based, Some(105));
            assert_eq!(target_span_start_0based, Some(10));
            assert_eq!(target_span_end_0based, Some(120));
            assert_eq!(mode, PairwiseAlignmentMode::Local);
            assert_eq!(match_score, 3);
            assert_eq!(mismatch_score, -4);
            assert_eq!(gap_open, -6);
            assert_eq!(gap_extend, -2);
        }
        other => panic!("expected AlignCompute, got {other:?}"),
    }
}

#[test]
fn parse_transcripts_derive_command() {
    let parsed = parse_shell_line(
        "transcripts derive seq_a --feature-id 11 --scope all_overlapping_both_strands --output-prefix seq_a__mrna",
    )
    .expect("parse transcripts derive");
    match parsed {
        ShellCommand::TranscriptsDerive {
            seq_id,
            feature_ids,
            scope,
            output_prefix,
        } => {
            assert_eq!(seq_id, "seq_a");
            assert_eq!(feature_ids, vec![11]);
            assert_eq!(scope, Some(SplicingScopePreset::AllOverlappingBothStrands));
            assert_eq!(output_prefix.as_deref(), Some("seq_a__mrna"));
        }
        other => panic!("expected TranscriptsDerive, got {other:?}"),
    }
}

#[test]
fn parse_rna_reads_commands() {
    let interpret = parse_shell_line(
            "rna-reads interpret seq_a 7 reads.fa --report-id tp73_reads --scope target_group_any_strand --kmer-len 9 --seed-stride-bp 1 --min-seed-hit-fraction 0.30 --min-weighted-seed-hit-fraction 0.05 --min-unique-matched-kmers 12 --min-chain-consistency-fraction 0.55 --max-median-transcript-gap 4.0 --min-confirmed-transitions 1 --min-transition-support-fraction 0.20 --no-cdna-poly-t-flip --poly-t-prefix-min-bp 20 --align-band-bp 24 --align-min-identity 0.60 --max-secondary-mappings 2",
        )
        .expect("parse rna-reads interpret");
    match interpret {
        ShellCommand::RnaReadsInterpret {
            seq_id,
            seed_feature_id,
            input_path,
            scope,
            origin_mode,
            target_gene_ids,
            roi_seed_capture_enabled,
            seed_filter,
            align_config,
            report_id,
            report_mode,
            checkpoint_path,
            checkpoint_every_reads,
            resume_from_checkpoint,
            ..
        } => {
            assert_eq!(seq_id, "seq_a");
            assert_eq!(seed_feature_id, 7);
            assert_eq!(input_path, "reads.fa");
            assert_eq!(scope, SplicingScopePreset::TargetGroupAnyStrand);
            assert_eq!(origin_mode, RnaReadOriginMode::SingleGene);
            assert!(target_gene_ids.is_empty());
            assert!(!roi_seed_capture_enabled);
            assert_eq!(seed_filter.kmer_len, 9);
            assert_eq!(seed_filter.seed_stride_bp, 1);
            assert!((seed_filter.min_seed_hit_fraction - 0.30).abs() < f64::EPSILON);
            assert!((seed_filter.min_weighted_seed_hit_fraction - 0.05).abs() < f64::EPSILON);
            assert_eq!(seed_filter.min_unique_matched_kmers, 12);
            assert!((seed_filter.min_chain_consistency_fraction - 0.55).abs() < f64::EPSILON);
            assert!((seed_filter.max_median_transcript_gap - 4.0).abs() < f64::EPSILON);
            assert_eq!(seed_filter.min_confirmed_exon_transitions, 1);
            assert!((seed_filter.min_transition_support_fraction - 0.20).abs() < f64::EPSILON);
            assert!(!seed_filter.cdna_poly_t_flip_enabled);
            assert_eq!(seed_filter.poly_t_prefix_min_bp, 20);
            assert_eq!(align_config.band_width_bp, 24);
            assert!((align_config.min_identity_fraction - 0.60).abs() < f64::EPSILON);
            assert_eq!(align_config.max_secondary_mappings, 2);
            assert_eq!(report_id.as_deref(), Some("tp73_reads"));
            assert_eq!(report_mode, RnaReadReportMode::Full);
            assert!(checkpoint_path.is_none());
            assert_eq!(checkpoint_every_reads, 10_000);
            assert!(!resume_from_checkpoint);
        }
        other => panic!("expected RnaReadsInterpret, got {other:?}"),
    }

    let interpret_multi = parse_shell_line(
            "rna-reads interpret seq_a 7 reads.fa --origin-mode multi_gene_sparse --target-gene TP73 --target-gene TP53 --roi-seed-capture",
        )
        .expect("parse rna-reads interpret multi-gene scaffold");
    match interpret_multi {
        ShellCommand::RnaReadsInterpret {
            origin_mode,
            target_gene_ids,
            roi_seed_capture_enabled,
            ..
        } => {
            assert_eq!(origin_mode, RnaReadOriginMode::MultiGeneSparse);
            assert_eq!(
                target_gene_ids,
                vec!["TP73".to_string(), "TP53".to_string()]
            );
            assert!(roi_seed_capture_enabled);
        }
        other => panic!("expected RnaReadsInterpret, got {other:?}"),
    }

    let interpret_checkpoint = parse_shell_line(
            "rna-reads interpret seq_a 7 reads.fa --report-mode seed_passed_only --checkpoint-path /tmp/tp53.chk.json --checkpoint-every-reads 2500 --resume-from-checkpoint",
        )
        .expect("parse rna-reads interpret checkpoint options");
    match interpret_checkpoint {
        ShellCommand::RnaReadsInterpret {
            report_mode,
            checkpoint_path,
            checkpoint_every_reads,
            resume_from_checkpoint,
            ..
        } => {
            assert_eq!(report_mode, RnaReadReportMode::SeedPassedOnly);
            assert_eq!(checkpoint_path.as_deref(), Some("/tmp/tp53.chk.json"));
            assert_eq!(checkpoint_every_reads, 2500);
            assert!(resume_from_checkpoint);
        }
        other => panic!("expected RnaReadsInterpret, got {other:?}"),
    }

    let list =
        parse_shell_line("rna-reads list-reports seq_a").expect("parse rna-reads list-reports");
    assert!(matches!(
        list,
        ShellCommand::RnaReadsListReports { seq_id } if seq_id.as_deref() == Some("seq_a")
    ));

    let show = parse_shell_line("rna-reads show-report tp73_reads").expect("parse rna-reads show");
    assert!(matches!(
        show,
        ShellCommand::RnaReadsShowReport { report_id } if report_id == "tp73_reads"
    ));

    let summarize = parse_shell_line(
        "rna-reads summarize-gene-support tp73_reads --gene TP73 --gene TP53 --record-indices 6,8 --complete-rule strict --output tp53_support.json",
    )
    .expect("parse rna-reads summarize-gene-support");
    assert!(matches!(
        summarize,
        ShellCommand::RnaReadsSummarizeGeneSupport {
            report_id,
            gene_ids,
            selected_record_indices,
            complete_rule,
            output_path,
        }
            if report_id == "tp73_reads"
                && gene_ids == vec!["TP73".to_string(), "TP53".to_string()]
                && selected_record_indices == vec![6, 8]
                && complete_rule == RnaReadGeneSupportCompleteRule::Strict
                && output_path.as_deref() == Some("tp53_support.json")
    ));

    let inspect_gene_support = parse_shell_line(
        "rna-reads inspect-gene-support tp73_reads --gene TP73 --gene TP53 --record-indices 6,8 --complete-rule exact --cohort rejected --output tp53_support_rows.json",
    )
    .expect("parse rna-reads inspect-gene-support");
    assert!(matches!(
        inspect_gene_support,
        ShellCommand::RnaReadsInspectGeneSupport {
            report_id,
            gene_ids,
            selected_record_indices,
            complete_rule,
            cohort_filter,
            output_path,
        }
            if report_id == "tp73_reads"
                && gene_ids == vec!["TP73".to_string(), "TP53".to_string()]
                && selected_record_indices == vec![6, 8]
                && complete_rule == RnaReadGeneSupportCompleteRule::Exact
                && cohort_filter == RnaReadGeneSupportAuditCohortFilter::Rejected
                && output_path.as_deref() == Some("tp53_support_rows.json")
    ));

    let inspect = parse_shell_line(
        "rna-reads inspect-alignments tp73_reads --selection aligned --limit 25 --effect-filter disagreement_only --sort score --search tp53 --record-indices 2,7,9 --score-bin-variant composite_seed_gate --score-bin-index 39 --score-bin-count 40",
    )
    .expect("parse rna-reads inspect-alignments");
    assert!(matches!(
        inspect,
        ShellCommand::RnaReadsInspectAlignments { report_id, selection, limit, effect_filter, sort_key, search, selected_record_indices, score_density_variant, score_bin_index, score_bin_count }
            if report_id == "tp73_reads"
                && selection == RnaReadHitSelection::Aligned
                && limit == 25
                && effect_filter == RnaReadAlignmentInspectionEffectFilter::DisagreementOnly
                && sort_key == RnaReadAlignmentInspectionSortKey::Score
                && search == "tp53"
                && selected_record_indices == vec![2, 7, 9]
                && score_density_variant == RnaReadScoreDensityVariant::CompositeSeedGate
                && score_bin_index == Some(39)
                && score_bin_count == 40
    ));

    let inspect_concatemers = parse_shell_line(
        "rna-reads inspect-concatemers tp73_reads --selection seed_passed --limit 15 --record-indices 2,7,9 --internal-homopolymer-min-bp 21 --end-margin-bp 40 --max-primary-query-cov 0.80 --min-secondary-identity 0.90 --max-secondary-query-overlap 0.10 --adapter-fasta data/resources/nanopore_direct_cdna_kit14_adapters.fasta --adapter-min-match-bp 18 --fragment-min-bp 75 --fragment-max-parts 3 --fragment-min-identity 0.88 --fragment-min-query-cov 0.45 --transcript-fasta data/transcriptome_cdna.fa.gz --transcript-fasta data/transcriptome_ncrna.fa.gz --transcript-index data/transcriptome_joined.index.json",
    )
    .expect("parse rna-reads inspect-concatemers");
    assert!(matches!(
        inspect_concatemers,
        ShellCommand::RnaReadsInspectConcatemers { report_id, selection, limit, selected_record_indices, settings }
            if report_id == "tp73_reads"
                && selection == RnaReadHitSelection::SeedPassed
                && limit == 15
                && selected_record_indices == vec![2, 7, 9]
                && settings.internal_homopolymer_min_bp == 21
                && settings.end_margin_bp == 40
                && (settings.max_primary_query_coverage_fraction - 0.80).abs() < 1e-9
                && (settings.min_secondary_identity_fraction - 0.90).abs() < 1e-9
                && (settings.max_secondary_query_overlap_fraction - 0.10).abs() < 1e-9
                && settings.adapter_fasta_path.as_deref() == Some("data/resources/nanopore_direct_cdna_kit14_adapters.fasta")
                && settings.adapter_min_match_bp == 18
                && settings.fragment_min_bp == 75
                && settings.fragment_max_parts == 3
                && (settings.fragment_min_identity_fraction - 0.88).abs() < 1e-9
                && (settings.fragment_min_query_coverage_fraction - 0.45).abs() < 1e-9
                && settings.transcript_fasta_paths == vec![
                    "data/transcriptome_cdna.fa.gz".to_string(),
                    "data/transcriptome_ncrna.fa.gz".to_string(),
                ]
                && settings.transcript_index_paths == vec![
                    "data/transcriptome_joined.index.json".to_string(),
                ]
    ));

    let build_transcript_index = parse_shell_line(
        "rna-reads build-transcript-index data/transcriptome.index.json --kmer-len 11 --transcript-fasta data/transcriptome_cdna.fa.gz --transcript-fasta data/transcriptome_ncrna.fa.gz",
    )
    .expect("parse rna-reads build-transcript-index");
    assert!(matches!(
        build_transcript_index,
        ShellCommand::RnaReadsBuildTranscriptIndex { path, seed_kmer_len, transcript_fasta_paths }
            if path == "data/transcriptome.index.json"
                && seed_kmer_len == 11
                && transcript_fasta_paths == vec![
                    "data/transcriptome_cdna.fa.gz".to_string(),
                    "data/transcriptome_ncrna.fa.gz".to_string(),
                ]
    ));

    let export = parse_shell_line("rna-reads export-report tp73_reads out.json")
        .expect("parse rna-reads export-report");
    assert!(matches!(
        export,
        ShellCommand::RnaReadsExportReport { report_id, path } if report_id == "tp73_reads" && path == "out.json"
    ));

    let export_hits = parse_shell_line(
        "rna-reads export-hits-fasta tp73_reads hits.fa --selection aligned --record-indices 4,9 --subset-spec filtered_tp53",
    )
    .expect("parse rna-reads export-hits-fasta");
    assert!(matches!(
        export_hits,
        ShellCommand::RnaReadsExportHitsFasta { report_id, path, selection, selected_record_indices, subset_spec }
            if report_id == "tp73_reads"
                && path == "hits.fa"
                && selection == RnaReadHitSelection::Aligned
                && selected_record_indices == vec![4, 9]
                && subset_spec.as_deref() == Some("filtered_tp53")
    ));

    let export_sheet = parse_shell_line(
        "rna-reads export-sample-sheet samples.tsv --seq-id seq_a --report-id tp73_reads --gene TP53 --complete-rule strict --append",
    )
    .expect("parse rna-reads export-sample-sheet");
    assert!(matches!(
        export_sheet,
        ShellCommand::RnaReadsExportSampleSheet { path, seq_id, report_ids, gene_ids, complete_rule, append }
            if path == "samples.tsv"
                && seq_id.as_deref() == Some("seq_a")
                && report_ids == vec!["tp73_reads".to_string()]
                && gene_ids == vec!["TP53".to_string()]
                && complete_rule == RnaReadGeneSupportCompleteRule::Strict
                && append
    ));

    let export_target_quality = parse_shell_line(
        "rna-reads export-target-quality tp73_reads target_quality.json --gene TP73 --gene TP53 --complete-rule exact",
    )
    .expect("parse rna-reads export-target-quality");
    assert!(matches!(
        export_target_quality,
        ShellCommand::RnaReadsExportTargetQuality {
            report_id,
            path,
            gene_ids,
            complete_rule
        } if report_id == "tp73_reads"
            && path == "target_quality.json"
            && gene_ids == vec!["TP73".to_string(), "TP53".to_string()]
            && complete_rule == RnaReadGeneSupportCompleteRule::Exact
    ));

    let export_paths =
        parse_shell_line(
            "rna-reads export-paths-tsv tp73_reads paths.tsv --selection seed_passed --record-indices 3,5 --subset-spec filtered_tp53",
        )
            .expect("parse rna-reads export-paths-tsv");
    assert!(matches!(
        export_paths,
        ShellCommand::RnaReadsExportExonPathsTsv { report_id, path, selection, selected_record_indices, subset_spec }
            if report_id == "tp73_reads"
                && path == "paths.tsv"
                && selection == RnaReadHitSelection::SeedPassed
                && selected_record_indices == vec![3, 5]
                && subset_spec.as_deref() == Some("filtered_tp53")
    ));

    let export_abundance = parse_shell_line(
        "rna-reads export-abundance-tsv tp73_reads abundance.tsv --selection aligned --record-indices 6,8 --subset-spec filtered_tp53",
    )
    .expect("parse rna-reads export-abundance-tsv");
    assert!(matches!(
        export_abundance,
        ShellCommand::RnaReadsExportExonAbundanceTsv { report_id, path, selection, selected_record_indices, subset_spec }
            if report_id == "tp73_reads"
                && path == "abundance.tsv"
                && selection == RnaReadHitSelection::Aligned
                && selected_record_indices == vec![6, 8]
                && subset_spec.as_deref() == Some("filtered_tp53")
    ));

    let export_score_density = parse_shell_line(
        "rna-reads export-score-density-svg tp73_reads score_density.svg --scale linear --variant composite_seed_gate",
    )
    .expect("parse rna-reads export-score-density-svg");
    assert!(matches!(
        export_score_density,
        ShellCommand::RnaReadsExportScoreDensitySvg { report_id, path, scale, variant }
            if report_id == "tp73_reads"
                && path == "score_density.svg"
                && scale == RnaReadScoreDensityScale::Linear
                && variant == RnaReadScoreDensityVariant::CompositeSeedGate
    ));

    let export_alignments_tsv = parse_shell_line(
        "rna-reads export-alignments-tsv tp73_reads alignments.tsv --selection aligned --limit 200 --record-indices 7,8 --subset-spec filtered_tp53",
    )
    .expect("parse rna-reads export-alignments-tsv");
    assert!(matches!(
        export_alignments_tsv,
        ShellCommand::RnaReadsExportAlignmentsTsv {
            report_id,
            path,
            selection,
            limit,
            selected_record_indices,
            subset_spec
        } if report_id == "tp73_reads"
            && path == "alignments.tsv"
            && selection == RnaReadHitSelection::Aligned
            && limit == Some(200)
            && selected_record_indices == vec![7, 8]
            && subset_spec.as_deref() == Some("filtered_tp53")
    ));

    let export_alignment_dotplot = parse_shell_line(
        "rna-reads export-alignment-dotplot-svg tp73_reads alignment_dotplot.svg --selection aligned --max-points 500",
    )
    .expect("parse rna-reads export-alignment-dotplot-svg");
    assert!(matches!(
        export_alignment_dotplot,
        ShellCommand::RnaReadsExportAlignmentDotplotSvg {
            report_id,
            path,
            selection,
            max_points
        } if report_id == "tp73_reads"
            && path == "alignment_dotplot.svg"
            && selection == RnaReadHitSelection::Aligned
            && max_points == 500
    ));

    let align = parse_shell_line(
        "rna-reads align-report tp73_reads --selection seed_passed --align-band-bp 32 --align-min-identity 0.70 --max-secondary-mappings 1",
    )
    .expect("parse rna-reads align-report");
    match align {
        ShellCommand::RnaReadsAlignReport {
            report_id,
            selection,
            align_config_override,
            selected_record_indices,
        } => {
            assert_eq!(report_id, "tp73_reads");
            assert_eq!(selection, RnaReadHitSelection::SeedPassed);
            assert!(selected_record_indices.is_empty());
            let cfg = align_config_override.expect("align config override");
            assert_eq!(cfg.band_width_bp, 32);
            assert!((cfg.min_identity_fraction - 0.70).abs() < f64::EPSILON);
            assert_eq!(cfg.max_secondary_mappings, 1);
        }
        other => panic!("expected RnaReadsAlignReport, got {other:?}"),
    }

    let align_selected =
        parse_shell_line("rna-reads align-report tp73_reads --record-indices 2,7,2,9")
            .expect("parse rna-reads align-report explicit record indices");
    match align_selected {
        ShellCommand::RnaReadsAlignReport {
            report_id,
            selection,
            align_config_override,
            selected_record_indices,
        } => {
            assert_eq!(report_id, "tp73_reads");
            assert_eq!(selection, RnaReadHitSelection::SeedPassed);
            assert!(align_config_override.is_none());
            assert_eq!(selected_record_indices, vec![2, 7, 9]);
        }
        other => panic!("expected RnaReadsAlignReport, got {other:?}"),
    }

    let align_selected_invalid =
        parse_shell_line("rna-reads align-report tp73_reads --record-indices ,,")
            .expect_err("expected invalid empty --record-indices to fail");
    assert!(
        align_selected_invalid.contains("Invalid --record-indices value"),
        "unexpected parse error: {align_selected_invalid}"
    );
}

#[test]
fn parse_rna_reads_interpret_defaults_to_engine_kmer_len() {
    let cmd = parse_shell_line(
        "rna-reads interpret seq_a 7 reads.fa --scope all_overlapping_both_strands",
    )
    .expect("parse rna-reads interpret defaults");
    match cmd {
        ShellCommand::RnaReadsInterpret {
            origin_mode,
            target_gene_ids,
            roi_seed_capture_enabled,
            seed_filter,
            report_mode,
            checkpoint_path,
            checkpoint_every_reads,
            resume_from_checkpoint,
            ..
        } => {
            assert_eq!(origin_mode, RnaReadOriginMode::SingleGene);
            assert!(target_gene_ids.is_empty());
            assert!(!roi_seed_capture_enabled);
            assert_eq!(
                seed_filter.kmer_len,
                RnaReadSeedFilterConfig::default().kmer_len
            );
            assert_eq!(seed_filter.kmer_len, 10);
            assert_eq!(report_mode, RnaReadReportMode::Full);
            assert!(checkpoint_path.is_none());
            assert_eq!(checkpoint_every_reads, 10_000);
            assert!(!resume_from_checkpoint);
        }
        other => panic!("expected RnaReadsInterpret, got {other:?}"),
    }
}

#[test]
fn execute_dotplot_and_flex_commands_store_payloads() {
    std::thread::Builder::new()
        .name("execute_dotplot_and_flex_commands_store_payloads".to_string())
        .stack_size(16 * 1024 * 1024)
        .spawn(|| {
            let mut state = ProjectState::default();
            state.sequences.insert(
                "seq_a".to_string(),
                DNAsequence::from_sequence("ATGCATGCATGCATGCATGCATGC").expect("sequence"),
            );
            state.sequences.insert(
                "ref".to_string(),
                DNAsequence::from_sequence("ATGCATGCATGCATGCATGCATGC").expect("reference"),
            );
            state.sequences.insert(
                "iso_a".to_string(),
                DNAsequence::from_sequence("ATGCATGCATGC").expect("isoform a"),
            );
            state.sequences.insert(
                "iso_b".to_string(),
                DNAsequence::from_sequence("GCATGCATGCAT").expect("isoform b"),
            );
            let mut engine = GentleEngine::from_state(state);

            let dotplot = execute_shell_command(
                &mut engine,
                &ShellCommand::DotplotCompute {
                    seq_id: "seq_a".to_string(),
                    reference_seq_id: None,
                    span_start_0based: Some(0),
                    span_end_0based: Some(24),
                    reference_span_start_0based: None,
                    reference_span_end_0based: None,
                    mode: DotplotMode::SelfForward,
                    word_size: 4,
                    step_bp: 2,
                    max_mismatches: 0,
                    tile_bp: None,
                    dotplot_id: Some("dp_1".to_string()),
                },
            )
            .expect("execute dotplot compute");
            assert!(dotplot.state_changed);
            assert_eq!(
                dotplot.output["dotplot"]["dotplot_id"].as_str(),
                Some("dp_1")
            );

            let overlay = execute_shell_command(
                &mut engine,
                &ShellCommand::DotplotOverlayCompute {
                    owner_seq_id: "ref".to_string(),
                    reference_seq_id: None,
                    reference_span_start_0based: Some(0),
                    reference_span_end_0based: Some(24),
                    queries: vec![
                        DotplotOverlayQuerySpec {
                            seq_id: "iso_a".to_string(),
                            label: "Isoform A".to_string(),
                            transcript_feature_id: None,
                            query_anchor_0based: None,
                            query_anchor_label: None,
                            span_start_0based: None,
                            span_end_0based: None,
                            mode: DotplotMode::PairForward,
                            color_rgb: Some([29, 78, 216]),
                        },
                        DotplotOverlayQuerySpec {
                            seq_id: "iso_b".to_string(),
                            label: "Isoform B".to_string(),
                            transcript_feature_id: None,
                            query_anchor_0based: None,
                            query_anchor_label: None,
                            span_start_0based: None,
                            span_end_0based: None,
                            mode: DotplotMode::PairForward,
                            color_rgb: Some([220, 38, 38]),
                        },
                    ],
                    word_size: 4,
                    step_bp: 2,
                    max_mismatches: 0,
                    tile_bp: None,
                    dotplot_id: Some("overlay_1".to_string()),
                },
            )
            .expect("execute overlay dotplot compute");
            assert!(overlay.state_changed);
            assert_eq!(
                overlay.output["dotplot"]["dotplot_id"].as_str(),
                Some("overlay_1")
            );
            assert_eq!(overlay.output["dotplot"]["series_count"].as_u64(), Some(2));

            let overlay_show = execute_shell_command(
                &mut engine,
                &ShellCommand::DotplotShow {
                    dotplot_id: "overlay_1".to_string(),
                },
            )
            .expect("show overlay dotplot");
            assert_eq!(
                overlay_show.output["dotplot"]["owner_seq_id"].as_str(),
                Some("ref")
            );
            assert_eq!(
                overlay_show.output["dotplot"]["query_series"]
                    .as_array()
                    .map(Vec::len),
                Some(2)
            );

            let overlay_svg_path = tempfile::NamedTempFile::new()
                .expect("tmp")
                .path()
                .with_extension("overlay.anchor.svg");
            let render_overlay = execute_shell_command(
                &mut engine,
                &ShellCommand::RenderDotplotSvg {
                    seq_id: "ref".to_string(),
                    dotplot_id: "overlay_1".to_string(),
                    output: overlay_svg_path.display().to_string(),
                    flex_track_id: None,
                    display_density_threshold: Some(0.0),
                    display_intensity_gain: Some(1.0),
                    overlay_x_axis_mode: DotplotOverlayXAxisMode::SharedExonAnchor,
                    overlay_anchor_exon: Some(DotplotOverlayAnchorExonRef {
                        start_1based: 10,
                        end_1based: 20,
                    }),
                },
            )
            .expect("render overlay dotplot shell command");
            assert!(!render_overlay.state_changed);

            let flex = execute_shell_command(
                &mut engine,
                &ShellCommand::FlexCompute {
                    seq_id: "seq_a".to_string(),
                    span_start_0based: Some(0),
                    span_end_0based: Some(24),
                    model: FlexibilityModel::AtRichness,
                    bin_bp: 6,
                    smoothing_bp: Some(12),
                    track_id: Some("flex_1".to_string()),
                },
            )
            .expect("execute flex compute");
            assert!(flex.state_changed);
            assert_eq!(flex.output["track"]["track_id"].as_str(), Some("flex_1"));

            let listed = execute_shell_command(
                &mut engine,
                &ShellCommand::DotplotList {
                    seq_id: Some("seq_a".to_string()),
                },
            )
            .expect("list dotplots");
            assert_eq!(listed.output["dotplot_count"].as_u64(), Some(1));
            let overlay_list = execute_shell_command(
                &mut engine,
                &ShellCommand::DotplotList {
                    seq_id: Some("ref".to_string()),
                },
            )
            .expect("list overlay dotplots");
            assert_eq!(overlay_list.output["dotplot_count"].as_u64(), Some(1));
            let flex_list = execute_shell_command(
                &mut engine,
                &ShellCommand::FlexList {
                    seq_id: Some("seq_a".to_string()),
                },
            )
            .expect("list flex tracks");
            assert_eq!(flex_list.output["track_count"].as_u64(), Some(1));
        })
        .expect("spawn dotplot shell test")
        .join()
        .expect("join dotplot shell test");
}

#[test]
fn execute_transcripts_derive_creates_sequences() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), tp53_isoform_test_sequence());
    let mut engine = GentleEngine::from_state(state);
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::TranscriptsDerive {
            seq_id: "seq_a".to_string(),
            feature_ids: vec![],
            scope: None,
            output_prefix: Some("seq_a__mrna".to_string()),
        },
    )
    .expect("execute transcripts derive");
    assert!(out.state_changed);
    let transcript_count = out.output["transcript_count"]
        .as_u64()
        .expect("transcript_count");
    assert!(transcript_count >= 1);
    let rows = out.output["transcripts"]
        .as_array()
        .expect("transcripts array");
    assert_eq!(rows.len() as u64, transcript_count);
}

#[test]
fn execute_splicing_refs_and_align_commands() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), tp53_isoform_test_sequence());
    state.sequences.insert(
        "query".to_string(),
        DNAsequence::from_sequence("ATGCGTAA").expect("query sequence"),
    );
    state.sequences.insert(
        "target".to_string(),
        DNAsequence::from_sequence("TTATGCGTAACCG").expect("target sequence"),
    );
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

    let derived = execute_shell_command(
        &mut engine,
        &ShellCommand::SplicingRefsDerive {
            seq_id: "seq_a".to_string(),
            span_start_0based: 0,
            span_end_0based: 400,
            seed_feature_id: Some(feature_id),
            scope: SplicingScopePreset::TargetGroupAnyStrand,
            output_prefix: Some("tp53_refs".to_string()),
        },
    )
    .expect("execute splicing-refs derive");
    assert!(derived.state_changed);
    let derived_ids = derived
        .output
        .get("derived_sequence_ids")
        .and_then(|value| value.as_array())
        .expect("derived sequence id array");
    assert!(
        derived_ids
            .iter()
            .filter_map(|value| value.as_str())
            .any(|id| id.starts_with("tp53_refs_dna")),
        "expected derived DNA-window sequence id in payload"
    );

    let align = execute_shell_command(
        &mut engine,
        &ShellCommand::AlignCompute {
            query_seq_id: "query".to_string(),
            target_seq_id: "target".to_string(),
            query_span_start_0based: Some(0),
            query_span_end_0based: Some(8),
            target_span_start_0based: None,
            target_span_end_0based: None,
            mode: PairwiseAlignmentMode::Local,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
        },
    )
    .expect("execute align compute");
    assert!(!align.state_changed);
    assert_eq!(align.output["alignment"]["mode"].as_str(), Some("local"));
    assert_eq!(
        align.output["alignment"]["query_seq_id"].as_str(),
        Some("query")
    );
    assert_eq!(
        align.output["alignment"]["target_seq_id"].as_str(),
        Some("target")
    );
    assert_eq!(
        align.output["result"]["created_seq_ids"]
            .as_array()
            .map(|v| v.len()),
        Some(0)
    );
}

#[test]
fn execute_rna_reads_commands_store_and_export_reports() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), tp53_isoform_test_sequence());
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
    let fasta_dir = tempdir().expect("tempdir");
    let input_path = fasta_dir.path().join("reads.fa");
    fs::write(
            &input_path,
            ">read_1\nATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTAATGGGCCCGGATTCCTTTTCTCTGTGAACCTTCCCGATGATGATGGAGGTGGAATGGAGGAGCCGCAGTCA\n",
        )
        .expect("write input fasta");
    let report_id = "rna_reads_test".to_string();
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.kmer_len = 3;
    seed_filter.min_seed_hit_fraction = 0.0;
    seed_filter.min_weighted_seed_hit_fraction = 0.0;
    seed_filter.min_unique_matched_kmers = 0;
    seed_filter.min_chain_consistency_fraction = 0.0;
    seed_filter.max_median_transcript_gap = 10_000.0;
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;
    let run = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsInterpret {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            input_path: input_path.display().to_string(),
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some(report_id.clone()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        },
    )
    .expect("execute rna-reads interpret");
    assert!(run.state_changed);
    assert_eq!(
        run.output["report"]["report_id"].as_str(),
        Some(report_id.as_str())
    );
    assert_eq!(run.output["report"]["read_count_total"].as_u64(), Some(1));
    assert_eq!(run.output["report"]["read_count_aligned"].as_u64(), Some(0));

    let align_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsAlignReport {
            report_id: report_id.clone(),
            selection: RnaReadHitSelection::SeedPassed,
            align_config_override: Some(RnaReadAlignConfig {
                band_width_bp: 24,
                min_identity_fraction: 0.60,
                max_secondary_mappings: 0,
            }),
            selected_record_indices: vec![],
        },
    )
    .expect("execute rna-reads align-report");
    assert!(align_result.state_changed);
    assert_eq!(
        align_result.output["report"]["report_id"].as_str(),
        Some(report_id.as_str())
    );
    assert!(
        align_result.output["report"]["read_count_aligned"]
            .as_u64()
            .is_some()
    );
    assert!(
        align_result.output["report"]["retained_count_msa_eligible"]
            .as_u64()
            .is_some()
    );

    let listed = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsListReports {
            seq_id: Some("seq_a".to_string()),
        },
    )
    .expect("list rna-read reports");
    assert_eq!(listed.output["report_count"].as_u64(), Some(1));
    assert_eq!(
        listed.output["reports"][0]["target_gene_count"].as_u64(),
        Some(0)
    );
    assert_eq!(
        listed.output["reports"][0]["roi_seed_capture_enabled"].as_bool(),
        Some(false)
    );
    assert!(
        listed.output["summary_rows"][0]
            .as_str()
            .is_some_and(|line| line.contains("origin=single_gene"))
    );
    assert!(
        listed.output["summary_rows"][0]
            .as_str()
            .is_some_and(|line| line.contains("msa_eligible(retained)="))
    );

    let shown = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsShowReport {
            report_id: report_id.clone(),
        },
    )
    .expect("show rna-read report");
    assert_eq!(
        shown.output["report"]["report_id"].as_str(),
        Some(report_id.as_str())
    );
    assert!(
        shown.output["summary"]
            .as_str()
            .is_some_and(|line| line.contains("mode=full") && line.contains("origin=single_gene"))
    );

    let inspected = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsInspectAlignments {
            report_id: report_id.clone(),
            selection: RnaReadHitSelection::Aligned,
            limit: 10,
            effect_filter: RnaReadAlignmentInspectionEffectFilter::SelectedOnly,
            sort_key: RnaReadAlignmentInspectionSortKey::Score,
            search: "read_1".to_string(),
            selected_record_indices: vec![0],
            score_density_variant: RnaReadScoreDensityVariant::CompositeSeedGate,
            score_bin_index: Some(39),
            score_bin_count: 40,
        },
    )
    .expect("inspect rna-read alignments");
    assert_eq!(
        inspected.output["inspection"]["report_id"].as_str(),
        Some(report_id.as_str())
    );
    assert_eq!(
        inspected.output["inspection"]["selection"].as_str(),
        Some("aligned")
    );
    assert_eq!(
        inspected.output["inspection"]["subset_spec"]["effect_filter"].as_str(),
        Some("selected_only")
    );
    assert_eq!(
        inspected.output["inspection"]["subset_spec"]["sort_key"].as_str(),
        Some("score")
    );
    assert_eq!(
        inspected.output["inspection"]["subset_spec"]["search"].as_str(),
        Some("read_1")
    );
    assert_eq!(
        inspected.output["inspection"]["subset_spec"]["selected_record_indices"]
            .as_array()
            .map(|values| values.len()),
        Some(1)
    );
    assert_eq!(
        inspected.output["inspection"]["subset_spec"]["score_density_variant"].as_str(),
        Some("composite_seed_gate")
    );
    assert_eq!(
        inspected.output["inspection"]["subset_spec"]["score_bin_index"].as_u64(),
        Some(39)
    );
    assert_eq!(
        inspected.output["inspection"]["subset_spec"]["score_bin_count"].as_u64(),
        Some(40)
    );
    assert_eq!(
        inspected.output["inspection"]["subset_match_count"].as_u64(),
        Some(0)
    );

    let concatemer_report = RnaReadInterpretationReport {
        schema: "gentle.rna_read_report.v1".to_string(),
        report_id: "rna_reads_concatemers_shell".to_string(),
        seq_id: "seq_concat".to_string(),
        align_config: RnaReadAlignConfig {
            max_secondary_mappings: 2,
            ..RnaReadAlignConfig::default()
        },
        hits: vec![RnaReadInterpretationHit {
            record_index: 0,
            header_id: "concatemer_candidate".to_string(),
            sequence: format!("{}{}{}", "C".repeat(35), "A".repeat(20), "G".repeat(35)),
            read_length_bp: 90,
            origin_class: RnaReadOriginClass::TargetPartialLocalBlock,
            best_mapping: Some(RnaReadMappingHit {
                transcript_id: "tx_primary".to_string(),
                transcript_label: "TX_PRIMARY".to_string(),
                strand: "+".to_string(),
                query_start_0based: 0,
                query_end_0based_exclusive: 38,
                query_coverage_fraction: 0.42,
                identity_fraction: 0.96,
                score: 180,
                ..RnaReadMappingHit::default()
            }),
            secondary_mappings: vec![RnaReadMappingHit {
                transcript_id: "tx_secondary".to_string(),
                transcript_label: "TX_SECONDARY".to_string(),
                strand: "+".to_string(),
                query_start_0based: 46,
                query_end_0based_exclusive: 86,
                query_coverage_fraction: 0.44,
                identity_fraction: 0.92,
                score: 150,
                ..RnaReadMappingHit::default()
            }],
            ..RnaReadInterpretationHit::default()
        }],
        ..RnaReadInterpretationReport::default()
    };
    let mut report_store_value = engine
        .state()
        .metadata
        .get("rna_read_reports")
        .cloned()
        .unwrap_or_else(|| {
            serde_json::json!({
                "schema": "gentle.rna_read_reports.v1",
                "reports": {}
            })
        });
    report_store_value["reports"]["rna_reads_concatemers_shell"] =
        serde_json::to_value(concatemer_report).expect("serialize concatemer report");
    engine
        .state_mut()
        .metadata
        .insert("rna_read_reports".to_string(), report_store_value);

    let concatemer_inspection = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsInspectConcatemers {
            report_id: "rna_reads_concatemers_shell".to_string(),
            selection: RnaReadHitSelection::All,
            limit: 10,
            selected_record_indices: vec![],
            settings: RnaReadConcatemerInspectionSettings::default(),
        },
    )
    .expect("inspect rna-read concatemers");
    assert_eq!(
        concatemer_inspection.output["inspection"]["report_id"].as_str(),
        Some("rna_reads_concatemers_shell")
    );
    assert_eq!(
        concatemer_inspection.output["inspection"]["suspicious_count"].as_u64(),
        Some(1)
    );
    assert_eq!(
        concatemer_inspection.output["inspection"]["strong_count"].as_u64(),
        Some(1)
    );
    assert_eq!(
        concatemer_inspection.output["inspection"]["partner_gene_summaries"]
            .as_array()
            .map(|values| values.len()),
        Some(0)
    );
    assert_eq!(
        concatemer_inspection.output["inspection"]["selected_record_indices"]
            .as_array()
            .map(|values| values.len()),
        Some(0)
    );
    assert_eq!(
        concatemer_inspection.output["inspection"]["rows"][0]["suspicion_level"].as_str(),
        Some("strong")
    );
    assert_eq!(
        concatemer_inspection.output["inspection"]["rows"][0]["top_disjoint_secondary_transcript_id"]
            .as_str(),
        Some("tx_secondary")
    );

    let transcript_file_gene3 = fasta_dir.path().join("transcriptome_gene3.fa");
    fs::write(
        &transcript_file_gene3,
        ">NM_GENE3_1 gene=GENE3 transcript=NM_GENE3_1\nATGCGTACGTTAGGCCATGCGTAC\n",
    )
    .expect("write gene3 transcript fasta");
    let transcript_file_gene4 = fasta_dir.path().join("transcriptome_gene4.fa");
    fs::write(
        &transcript_file_gene4,
        ">NM_GENE4_1 gene=GENE4 transcript=NM_GENE4_1\nCGTTAACCGGTTACCGGTTAACCG\n",
    )
    .expect("write gene4 transcript fasta");
    let transcript_index_output = fasta_dir.path().join("transcriptome.index.json");
    let build_index_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsBuildTranscriptIndex {
            path: transcript_index_output.display().to_string(),
            seed_kmer_len: 10,
            transcript_fasta_paths: vec![
                transcript_file_gene3.display().to_string(),
                transcript_file_gene4.display().to_string(),
            ],
        },
    )
    .expect("build transcript index");
    assert_eq!(
        build_index_result.output["path"].as_str(),
        Some(transcript_index_output.to_string_lossy().as_ref())
    );
    assert_eq!(
        build_index_result.output["index"]["schema"].as_str(),
        Some("gentle.rna_read_transcript_catalog_index.v1")
    );
    assert_eq!(
        build_index_result.output["index"]["transcript_count"].as_u64(),
        Some(2)
    );
    assert_eq!(
        build_index_result.output["index"]["gene_count"].as_u64(),
        Some(2)
    );
    assert!(transcript_index_output.exists());

    let exported_report = fasta_dir.path().join("report.json");
    let export_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsExportReport {
            report_id: report_id.clone(),
            path: exported_report.display().to_string(),
        },
    )
    .expect("export rna-read report");
    assert_eq!(
        export_result.output["report_id"].as_str(),
        Some(report_id.as_str())
    );
    assert!(exported_report.exists());

    let exported_hits = fasta_dir.path().join("hits.fa");
    let export_hits_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsExportHitsFasta {
            report_id,
            path: exported_hits.display().to_string(),
            selection: RnaReadHitSelection::All,
            selected_record_indices: vec![0],
            subset_spec: Some("filter=all aligned | sort=rank | search=<none>".to_string()),
        },
    )
    .expect("export rna-read hits");
    assert_eq!(
        export_hits_result.output["written_records"].as_u64(),
        Some(1)
    );
    assert_eq!(
        export_hits_result.output["subset_spec"].as_str(),
        Some("filter=all aligned | sort=rank | search=<none>")
    );
    let fasta_text = fs::read_to_string(exported_hits).expect("read exported hits");
    assert!(fasta_text.contains(">read_1"));
    assert!(fasta_text.contains("subset_spec=filter=all aligned | sort=rank | search=<none>"));
    assert!(
        fasta_text.contains(" mode=") || fasta_text.contains(" best=none"),
        "FASTA headers should include alignment mode when a best mapping exists"
    );

    let exported_sheet = fasta_dir.path().join("samples.tsv");
    let export_sheet_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsExportSampleSheet {
            path: exported_sheet.display().to_string(),
            seq_id: Some("seq_a".to_string()),
            report_ids: vec![],
            gene_ids: vec![],
            complete_rule: RnaReadGeneSupportCompleteRule::Near,
            append: false,
        },
    )
    .expect("export rna-read sample sheet");
    assert_eq!(export_sheet_result.output["report_count"].as_u64(), Some(1));
    assert_eq!(
        export_sheet_result.output["complete_rule"].as_str(),
        Some("near")
    );
    let sheet_text = fs::read_to_string(exported_sheet).expect("read sample sheet");
    assert!(sheet_text.contains("sample_id"));
    assert!(sheet_text.contains("mean_read_length_bp"));
    assert!(sheet_text.contains("exon_support_frequencies_json"));
    assert!(sheet_text.contains("origin_mode"));
    assert!(sheet_text.contains("target_gene_ids_json"));

    let exported_target_quality = fasta_dir.path().join("target_quality.json");
    let export_target_quality_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsExportTargetQuality {
            report_id: "rna_reads_test".to_string(),
            path: exported_target_quality.display().to_string(),
            gene_ids: vec!["TP53".to_string()],
            complete_rule: RnaReadGeneSupportCompleteRule::Near,
        },
    )
    .expect("export rna-read target quality");
    assert_eq!(
        export_target_quality_result.output["format"].as_str(),
        Some("json")
    );
    assert_eq!(
        export_target_quality_result.output["gene_ids"]
            .as_array()
            .map(|values| values.len()),
        Some(1)
    );
    assert_eq!(
        export_target_quality_result.output["entry_count"].as_u64(),
        Some(1)
    );
    let target_quality_text =
        fs::read_to_string(exported_target_quality).expect("read target quality export");
    assert!(target_quality_text.contains("gentle.rna_read_target_quality_comparison_bundle.v1"));
    assert!(target_quality_text.contains("\"requested_gene_ids\""));
    assert!(target_quality_text.contains("\"gentle_version\""));

    let exported_paths = fasta_dir.path().join("paths.tsv");
    let export_paths_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsExportExonPathsTsv {
            report_id: "rna_reads_test".to_string(),
            path: exported_paths.display().to_string(),
            selection: RnaReadHitSelection::All,
            selected_record_indices: vec![0],
            subset_spec: Some("filter=selected only | sort=score | search=tp53".to_string()),
        },
    )
    .expect("export rna-read exon paths");
    assert_eq!(export_paths_result.output["row_count"].as_u64(), Some(1));
    assert_eq!(
        export_paths_result.output["selected_record_indices"]
            .as_array()
            .map(|values| values.len()),
        Some(1)
    );
    assert_eq!(
        export_paths_result.output["subset_spec"].as_str(),
        Some("filter=selected only | sort=score | search=tp53")
    );
    let paths_text = fs::read_to_string(exported_paths).expect("read path sheet");
    assert!(paths_text.contains("selected_record_indices=0"));
    assert!(paths_text.contains("subset_spec=filter=selected only | sort=score | search=tp53"));
    assert!(paths_text.contains("exon_path"));
    assert!(paths_text.contains("reverse_complement_applied"));
    assert!(paths_text.contains("best_alignment_mode"));

    let exported_abundance = fasta_dir.path().join("abundance.tsv");
    let export_abundance_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsExportExonAbundanceTsv {
            report_id: "rna_reads_test".to_string(),
            path: exported_abundance.display().to_string(),
            selection: RnaReadHitSelection::All,
            selected_record_indices: vec![0],
            subset_spec: Some("filter=selected only | sort=score | search=tp53".to_string()),
        },
    )
    .expect("export rna-read abundance");
    assert_eq!(
        export_abundance_result.output["selected_read_count"].as_u64(),
        Some(1)
    );
    assert_eq!(
        export_abundance_result.output["selected_record_indices"]
            .as_array()
            .map(|values| values.len()),
        Some(1)
    );
    assert_eq!(
        export_abundance_result.output["subset_spec"].as_str(),
        Some("filter=selected only | sort=score | search=tp53")
    );
    let abundance_text = fs::read_to_string(exported_abundance).expect("read abundance sheet");
    assert!(abundance_text.contains("selected_record_indices=0"));
    assert!(abundance_text.contains("subset_spec=filter=selected only | sort=score | search=tp53"));
    assert!(abundance_text.contains("row_kind"));

    let exported_density_svg = fasta_dir.path().join("score_density.svg");
    let export_density_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsExportScoreDensitySvg {
            report_id: "rna_reads_test".to_string(),
            path: exported_density_svg.display().to_string(),
            scale: RnaReadScoreDensityScale::Log,
            variant: RnaReadScoreDensityVariant::CompositeSeedGate,
        },
    )
    .expect("export rna-read score density svg");
    assert_eq!(export_density_result.output["scale"].as_str(), Some("log"));
    let density_text = fs::read_to_string(exported_density_svg).expect("read density svg");
    assert!(density_text.contains("<svg"));
    assert!(density_text.contains("seed-hit score density"));

    let exported_alignments_tsv = fasta_dir.path().join("alignments.tsv");
    let export_alignments_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsExportAlignmentsTsv {
            report_id: "rna_reads_test".to_string(),
            path: exported_alignments_tsv.display().to_string(),
            selection: RnaReadHitSelection::Aligned,
            limit: Some(50),
            selected_record_indices: vec![0],
            subset_spec: Some("filter=selected only | sort=score | search=tp53".to_string()),
        },
    )
    .expect("export rna-read alignments tsv");
    assert_eq!(
        export_alignments_result.output["selection"].as_str(),
        Some("aligned")
    );
    assert_eq!(
        export_alignments_result.output["selected_record_indices"]
            .as_array()
            .map(|values| values.len()),
        Some(1)
    );
    assert_eq!(
        export_alignments_result.output["subset_spec"].as_str(),
        Some("filter=selected only | sort=score | search=tp53")
    );
    let alignments_text = fs::read_to_string(exported_alignments_tsv).expect("read alignments tsv");
    assert!(alignments_text.contains("alignment_mode"));
    assert!(alignments_text.contains("identity_fraction"));
    assert!(alignments_text.contains("selected_record_indices=0"));
    assert!(
        alignments_text.contains("subset_spec=filter=selected only | sort=score | search=tp53")
    );

    let exported_alignment_dotplot_svg = fasta_dir.path().join("alignment_dotplot.svg");
    let export_alignment_dotplot_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsExportAlignmentDotplotSvg {
            report_id: "rna_reads_test".to_string(),
            path: exported_alignment_dotplot_svg.display().to_string(),
            selection: RnaReadHitSelection::Aligned,
            max_points: 1_000,
        },
    )
    .expect("export rna-read alignment dotplot svg");
    assert_eq!(
        export_alignment_dotplot_result.output["selection"].as_str(),
        Some("aligned")
    );
    let alignment_dotplot_text =
        fs::read_to_string(exported_alignment_dotplot_svg).expect("read alignment dotplot svg");
    assert!(alignment_dotplot_text.contains("<svg"));
    assert!(alignment_dotplot_text.contains("alignment dotplot"));
}

#[test]
fn execute_rna_reads_export_sample_sheet_supports_target_gene_columns() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), tp53_isoform_test_sequence());
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
    let fasta_dir = tempdir().expect("tempdir");
    let input_path = fasta_dir.path().join("reads.fa");
    fs::write(
        &input_path,
        ">read_1\nATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTAATGGGCCCGGATTCCTTTTCTCTGTGAACCTTCCCGATGATGATGGAGGTGGAATGGAGGAGCCGCAGTCA\n",
    )
    .expect("write input fasta");
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.kmer_len = 3;
    seed_filter.min_seed_hit_fraction = 0.0;
    seed_filter.min_weighted_seed_hit_fraction = 0.0;
    seed_filter.min_unique_matched_kmers = 0;
    seed_filter.min_chain_consistency_fraction = 0.0;
    seed_filter.max_median_transcript_gap = 10_000.0;
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;
    execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsInterpret {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            input_path: input_path.display().to_string(),
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("rna_reads_gene_sheet".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        },
    )
    .expect("interpret RNA reads");
    let exported_sheet = fasta_dir.path().join("gene_samples.tsv");
    let export_sheet_result = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsExportSampleSheet {
            path: exported_sheet.display().to_string(),
            seq_id: Some("seq_a".to_string()),
            report_ids: vec!["rna_reads_gene_sheet".to_string()],
            gene_ids: vec!["TP53".to_string()],
            complete_rule: RnaReadGeneSupportCompleteRule::Near,
            append: false,
        },
    )
    .expect("export target-gene RNA-read sample sheet");
    assert_eq!(export_sheet_result.output["report_count"].as_u64(), Some(1));
    assert_eq!(
        export_sheet_result.output["gene_ids"]
            .as_array()
            .map(|values| values.len()),
        Some(1)
    );
    assert_eq!(
        export_sheet_result.output["complete_rule"].as_str(),
        Some("near")
    );
    let sheet_text = fs::read_to_string(exported_sheet).expect("read sample sheet");
    assert!(sheet_text.contains("gene_support_accepted_target_count"));
    assert!(sheet_text.contains("gene_support_exon_pair_support_json"));
    assert!(sheet_text.contains("gene_support_mean_assigned_read_length_bp"));
    assert!(sheet_text.contains("[\"TP53\"]"));
}

#[test]
fn execute_rna_reads_summarize_gene_support_returns_summary_json_and_writes_file() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), tp53_isoform_test_sequence());
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
    let dna = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present");
    let read_sequence =
        String::from_utf8(dna.forward_bytes()[98..560].to_vec()).expect("TP53 transcript slice");
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("tp53_gene_support_reads.fa");
    fs::write(&input_path, format!(">tp53_full\n{read_sequence}\n")).expect("write input fasta");
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.kmer_len = 3;
    seed_filter.min_seed_hit_fraction = 0.0;
    seed_filter.min_weighted_seed_hit_fraction = 0.0;
    seed_filter.min_unique_matched_kmers = 0;
    seed_filter.min_chain_consistency_fraction = 0.0;
    seed_filter.max_median_transcript_gap = 10_000.0;
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;

    execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsInterpret {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            input_path: input_path.display().to_string(),
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("tp53_gene_support".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        },
    )
    .expect("interpret TP53 RNA reads");
    execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsAlignReport {
            report_id: "tp53_gene_support".to_string(),
            selection: RnaReadHitSelection::SeedPassed,
            align_config_override: Some(RnaReadAlignConfig {
                band_width_bp: 24,
                min_identity_fraction: 0.60,
                max_secondary_mappings: 0,
            }),
            selected_record_indices: vec![],
        },
    )
    .expect("align TP53 RNA-read report");
    let source_report = engine
        .get_rna_read_report("tp53_gene_support")
        .expect("read persisted source report");

    let output_path = td.path().join("tp53_gene_support.json");
    let output = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsSummarizeGeneSupport {
            report_id: "tp53_gene_support".to_string(),
            gene_ids: vec!["TP53".to_string(), "TP73".to_string()],
            selected_record_indices: vec![0],
            complete_rule: RnaReadGeneSupportCompleteRule::Near,
            output_path: Some(output_path.display().to_string()),
        },
    )
    .expect("execute gene-support summary");
    assert!(!output.state_changed);
    assert_eq!(
        output.output["schema"].as_str(),
        Some("gentle.rna_read_gene_support_summary.v1")
    );
    assert_eq!(
        output.output["report_id"].as_str(),
        Some("tp53_gene_support")
    );
    assert!(output.output["generated_at_unix_ms"].as_u64().unwrap_or(0) > 0);
    assert!(output.output["op_id"].as_str().is_some());
    assert!(output.output["run_id"].as_str().is_some());
    assert_eq!(
        output.output["source_report_generated_at_unix_ms"].as_u64(),
        Some(source_report.generated_at_unix_ms as u64)
    );
    assert_eq!(
        output.output["source_report_op_id"].as_str(),
        source_report.op_id.as_deref()
    );
    assert_eq!(
        output.output["source_report_run_id"].as_str(),
        source_report.run_id.as_deref()
    );
    assert_eq!(
        output.output["requested_gene_ids"]
            .as_array()
            .map(|values| values.len()),
        Some(2)
    );
    assert_eq!(output.output["matched_gene_ids"][0].as_str(), Some("TP53"));
    assert_eq!(output.output["missing_gene_ids"][0].as_str(), Some("TP73"));
    assert_eq!(
        output.output["selected_record_indices"][0].as_u64(),
        Some(0)
    );
    assert_eq!(output.output["complete_rule"].as_str(), Some("near"));
    assert_eq!(output.output["accepted_target_count"].as_u64(), Some(1));
    assert_eq!(output.output["complete_count"].as_u64(), Some(1));
    assert_eq!(output.output["all_target"]["read_count"].as_u64(), Some(1));
    assert!(output_path.exists());
    let written = fs::read_to_string(&output_path).expect("read written summary");
    let written_json: serde_json::Value =
        serde_json::from_str(&written).expect("written summary json");
    assert_eq!(written_json, output.output);
}

#[test]
fn execute_rna_reads_inspect_gene_support_returns_audit_json_and_writes_file() {
    let mut state = ProjectState::default();
    state
        .sequences
        .insert("seq_a".to_string(), tp53_isoform_test_sequence());
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
    let dna = engine
        .state()
        .sequences
        .get("seq_a")
        .expect("sequence present");
    let read_sequence =
        String::from_utf8(dna.forward_bytes()[98..560].to_vec()).expect("TP53 transcript slice");
    let td = tempdir().expect("tempdir");
    let input_path = td.path().join("tp53_gene_support_reads.fa");
    fs::write(&input_path, format!(">tp53_full\n{read_sequence}\n")).expect("write input fasta");
    let mut seed_filter = RnaReadSeedFilterConfig::default();
    seed_filter.kmer_len = 3;
    seed_filter.min_seed_hit_fraction = 0.0;
    seed_filter.min_weighted_seed_hit_fraction = 0.0;
    seed_filter.min_unique_matched_kmers = 0;
    seed_filter.min_chain_consistency_fraction = 0.0;
    seed_filter.max_median_transcript_gap = 10_000.0;
    seed_filter.min_confirmed_exon_transitions = 0;
    seed_filter.min_transition_support_fraction = 0.0;

    execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsInterpret {
            seq_id: "seq_a".to_string(),
            seed_feature_id: feature_id,
            input_path: input_path.display().to_string(),
            profile: RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_format: RnaReadInputFormat::Fasta,
            scope: SplicingScopePreset::AllOverlappingBothStrands,
            origin_mode: RnaReadOriginMode::SingleGene,
            target_gene_ids: vec![],
            roi_seed_capture_enabled: false,
            seed_filter,
            align_config: RnaReadAlignConfig::default(),
            report_id: Some("tp53_gene_support".to_string()),
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: 10_000,
            resume_from_checkpoint: false,
        },
    )
    .expect("interpret TP53 RNA reads");
    execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsAlignReport {
            report_id: "tp53_gene_support".to_string(),
            selection: RnaReadHitSelection::SeedPassed,
            align_config_override: Some(RnaReadAlignConfig {
                band_width_bp: 24,
                min_identity_fraction: 0.60,
                max_secondary_mappings: 0,
            }),
            selected_record_indices: vec![],
        },
    )
    .expect("align TP53 RNA-read report");
    let source_report = engine
        .get_rna_read_report("tp53_gene_support")
        .expect("read persisted source report");

    let output_path = td.path().join("tp53_gene_support_audit.json");
    let output = execute_shell_command(
        &mut engine,
        &ShellCommand::RnaReadsInspectGeneSupport {
            report_id: "tp53_gene_support".to_string(),
            gene_ids: vec!["TP53".to_string(), "TP73".to_string()],
            selected_record_indices: vec![0],
            complete_rule: RnaReadGeneSupportCompleteRule::Near,
            cohort_filter: RnaReadGeneSupportAuditCohortFilter::All,
            output_path: Some(output_path.display().to_string()),
        },
    )
    .expect("execute gene-support audit");
    assert!(!output.state_changed);
    assert_eq!(
        output.output["schema"].as_str(),
        Some("gentle.rna_read_gene_support_audit.v1")
    );
    assert_eq!(
        output.output["report_id"].as_str(),
        Some("tp53_gene_support")
    );
    assert!(output.output["generated_at_unix_ms"].as_u64().unwrap_or(0) > 0);
    assert!(output.output["op_id"].as_str().is_some());
    assert!(output.output["run_id"].as_str().is_some());
    assert_eq!(
        output.output["source_report_generated_at_unix_ms"].as_u64(),
        Some(source_report.generated_at_unix_ms as u64)
    );
    assert_eq!(
        output.output["source_report_op_id"].as_str(),
        source_report.op_id.as_deref()
    );
    assert_eq!(
        output.output["source_report_run_id"].as_str(),
        source_report.run_id.as_deref()
    );
    assert_eq!(
        output.output["requested_gene_ids"]
            .as_array()
            .map(|values| values.len()),
        Some(2)
    );
    assert_eq!(output.output["matched_gene_ids"][0].as_str(), Some("TP53"));
    assert_eq!(output.output["missing_gene_ids"][0].as_str(), Some("TP73"));
    assert_eq!(
        output.output["selected_record_indices"][0].as_u64(),
        Some(0)
    );
    assert_eq!(output.output["complete_rule"].as_str(), Some("near"));
    assert_eq!(output.output["cohort_filter"].as_str(), Some("all"));
    assert_eq!(output.output["evaluated_row_count"].as_u64(), Some(1));
    assert_eq!(output.output["row_count"].as_u64(), Some(1));
    assert_eq!(
        output.output["accepted_target_record_indices"][0].as_u64(),
        Some(0)
    );
    assert_eq!(
        output.output["complete_record_indices"][0].as_u64(),
        Some(0)
    );
    assert_eq!(
        output.output["rows"][0]["status"].as_str(),
        Some("accepted_complete")
    );
    assert_eq!(
        output.output["rows"][0]["status_reason"].as_str(),
        Some("requested_gene_meets_complete_rule")
    );
    assert!(output_path.exists());
    let written = fs::read_to_string(&output_path).expect("read written audit");
    let mut written_json: serde_json::Value =
        serde_json::from_str(&written).expect("written audit json");
    let mut output_json = output.output.clone();
    let written_identity = written_json["rows"][0]["identity_fraction"]
        .as_f64()
        .expect("written identity fraction");
    let output_identity = output_json["rows"][0]["identity_fraction"]
        .as_f64()
        .expect("output identity fraction");
    assert!(
        (written_identity - output_identity).abs() <= 1e-12,
        "written/output identity_fraction drifted more than expected: {written_identity} vs {output_identity}"
    );
    written_json["rows"][0]["identity_fraction"] = serde_json::Value::Null;
    output_json["rows"][0]["identity_fraction"] = serde_json::Value::Null;
    assert_eq!(written_json, output_json);
}

#[test]
fn execute_set_param_updates_display_state() {
    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::SetParameter {
            name: "vcf_display_min_qual".to_string(),
            value_json: "33.5".to_string(),
        },
    )
    .expect("execute set-param");
    assert!(out.state_changed);
    assert!(
        (engine.state().display.vcf_display_min_qual - 33.5).abs() < f64::EPSILON,
        "vcf_display_min_qual should be updated by set-param"
    );
}

#[test]
fn execute_set_param_updates_tfbs_display_state() {
    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::SetParameter {
            name: "tfbs_display_min_llr_quantile".to_string(),
            value_json: "0.85".to_string(),
        },
    )
    .expect("execute set-param");
    assert!(out.state_changed);
    assert!(
        (engine.state().display.tfbs_display_min_llr_quantile - 0.85).abs() < f64::EPSILON,
        "tfbs_display_min_llr_quantile should be updated by set-param"
    );
}

#[test]
fn execute_set_param_updates_restriction_display_state() {
    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::SetParameter {
            name: "preferred_restriction_enzymes".to_string(),
            value_json: r#"[" BamHI ","EcoRI","BamHI"]"#.to_string(),
        },
    )
    .expect("execute restriction set-param");
    assert!(out.state_changed);
    assert_eq!(
        engine.state().display.preferred_restriction_enzymes,
        vec!["BamHI".to_string(), "EcoRI".to_string()],
        "preferred_restriction_enzymes should be updated by set-param"
    );
}

#[test]
fn execute_op_set_display_visibility_marks_state_changed() {
    let mut engine = GentleEngine::new();
    let out = execute_shell_command(
        &mut engine,
        &ShellCommand::Op {
            payload: r#"{"SetDisplayVisibility":{"target":"Tfbs","visible":true}}"#.to_string(),
        },
    )
    .expect("execute op");
    assert!(out.state_changed);
    assert!(engine.state().display.show_tfbs);
}
