//! Genome and helper-genome infrastructure.
//!
//! This module owns catalog-driven source resolution, preparation/indexing,
//! annotation normalization, extraction, and BLAST helpers used across GUI, CLI,
//! JavaScript, and Lua adapters.
//!
//! Catalog and source-resolution model:
//! - Catalog entries describe sequence/annotation origins using local paths,
//!   remote URLs, and optional NCBI assembly/GenBank accessions.
//! - Preparation resolves those descriptors into concrete, reproducible local
//!   artifacts under the configured cache/install layout.
//! - Helper-genome catalog flows share the same preparation contracts so
//!   retrieval/extraction behavior remains adapter-equivalent.
//!
//! Preparation/indexing lifecycle and side effects:
//! - Fetch/copy sequence + annotation inputs into a prepared install directory.
//! - Build FASTA index sidecars (`.fai`) and gene-index sidecars for fast lookup.
//! - Optionally build BLAST databases and record executable/index metadata.
//! - Persist install manifest fields (sources, checksums, paths, timestamps,
//!   source types) so provenance and reuse decisions are auditable.
//!
//! Annotation parsing behavior:
//! - GenBank remains the canonical annotation import path.
//! - XML (`GBSet/GBSeq`) and tabular GFF/GTF parsing are additive and normalized into shared
//!   `GenomeGeneRecord` structures consumed by engine operations.
//! - Malformed annotation lines are skipped with bounded warning summaries
//!   (including capped file/line context) rather than causing silent drift.
//!
//! Extraction and BLAST contracts:
//! - Region/gene extraction is coordinate-driven and deterministic by prepared
//!   source + chromosome + range.
//! - BLAST wrappers expose typed report records and preserve tool invocation
//!   context needed for diagnostics.
//! - Expected failure modes include missing catalogs/sources, missing external
//!   binaries, malformed records, network/download issues, and caller-requested
//!   cancellation/timebox exits.

use crate::feature_location::feature_is_reverse;
use crate::ncbi_genbank_xml::parse_gbseq_xml_file;
use flate2::read::MultiGzDecoder;
use serde::{Deserialize, Serialize};
use sha1::{Digest, Sha1};
use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Output, Stdio};
use std::thread;
use std::time::Duration;
use tempfile::NamedTempFile;

pub const DEFAULT_GENOME_CATALOG_PATH: &str = "assets/genomes.json";
pub const DEFAULT_HELPER_GENOME_CATALOG_PATH: &str = "assets/helper_genomes.json";
pub const DEFAULT_GENOME_CACHE_DIR: &str = "data/genomes";
pub const DEFAULT_MAKEBLASTDB_BIN: &str = "makeblastdb";
pub const DEFAULT_BLASTN_BIN: &str = "blastn";
pub const MAKEBLASTDB_ENV_BIN: &str = "GENTLE_MAKEBLASTDB_BIN";
pub const BLASTN_ENV_BIN: &str = "GENTLE_BLASTN_BIN";
pub const PREPARE_GENOME_TIMEOUT_SECS_ENV: &str = "GENTLE_PREPARE_GENOME_TIMEOUT_SECS";
const DEFAULT_NCBI_EFETCH_ENDPOINT: &str =
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
const NCBI_EFETCH_ENV_VAR: &str = "GENTLE_NCBI_EFETCH_URL";
const HTTP_RETRY_ATTEMPTS: usize = 4;
const HTTP_RETRY_BASE_BACKOFF_MS: u64 = 1000;
const HTTP_CONNECT_TIMEOUT_SECS: u64 = 20;
const HTTP_READ_TIMEOUT_SECS: u64 = 120;
const BLASTN_OUTFMT_FIELDS: &str =
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs";
const PREPARE_CANCELLED_BY_CALLER: &str = "Genome preparation cancelled by caller";
const BLAST_CANCELLED_BY_CALLER: &str = "BLAST search cancelled by caller";
const ANNOTATION_PARSE_ISSUE_LIMIT: usize = 12;
const ANNOTATION_PARSE_CONTEXT_CHARS: usize = 140;

#[cfg(test)]
pub(crate) fn genbank_env_lock() -> &'static std::sync::Mutex<()> {
    use std::sync::{Mutex, OnceLock};
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(()))
}

/// Catalog entry describing where to fetch one genome assembly and annotation.
#[derive(Default, Deserialize, Serialize, Debug, Clone)]
pub struct GenomeCatalogEntry {
    pub ncbi_taxonomy_id: Option<u32>,
    pub description: Option<String>,
    pub ncbi_assembly_accession: Option<String>,
    pub ncbi_assembly_name: Option<String>,
    pub genbank_accession: Option<String>,
    pub genbank_reference: Option<String>,
    pub local_variant_unpublished: Option<bool>,
    pub sequence_remote: Option<String>,
    pub annotations_remote: Option<String>,
    pub sequence_local: Option<String>,
    pub annotations_local: Option<String>,
    #[serde(default = "default_cache_dir")]
    pub cache_dir: Option<String>,
    #[serde(default)]
    pub nucleotide_length_bp: Option<usize>,
    #[serde(default)]
    pub molecular_mass_da: Option<f64>,
}

fn default_cache_dir() -> Option<String> {
    Some(DEFAULT_GENOME_CACHE_DIR.to_string())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GenomeInstallManifest {
    genome_id: String,
    sequence_source: String,
    annotation_source: String,
    #[serde(default)]
    sequence_source_type: Option<String>,
    #[serde(default)]
    annotation_source_type: Option<String>,
    #[serde(default)]
    sequence_sha1: Option<String>,
    #[serde(default)]
    annotation_sha1: Option<String>,
    sequence_path: String,
    annotation_path: String,
    fasta_index_path: String,
    gene_index_path: Option<String>,
    #[serde(default)]
    blast_db_prefix: Option<String>,
    #[serde(default)]
    blast_index_executable: Option<String>,
    #[serde(default)]
    blast_indexed_at_unix_ms: Option<u128>,
    installed_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrepareGenomeReport {
    pub genome_id: String,
    pub reused_existing: bool,
    pub sequence_path: String,
    pub annotation_path: String,
    pub fasta_index_path: String,
    pub sequence_source: Option<String>,
    pub annotation_source: Option<String>,
    pub sequence_source_type: Option<String>,
    pub annotation_source_type: Option<String>,
    pub blast_db_prefix: Option<String>,
    pub blast_index_ready: bool,
    pub blast_index_executable: Option<String>,
    #[serde(default)]
    pub cached_contig_count: usize,
    #[serde(default)]
    pub cached_total_span_bp: u64,
    #[serde(default)]
    pub cached_longest_contig: Option<String>,
    #[serde(default)]
    pub cached_longest_contig_bp: Option<u64>,
    #[serde(default)]
    pub cached_contig_preview: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub annotation_parse_report: Option<AnnotationParseReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalBinaryPreflightProbe {
    pub tool: String,
    pub env_var: String,
    pub executable: String,
    pub resolved_path: Option<String>,
    pub found: bool,
    pub version_probe_ok: bool,
    pub status_code: Option<i32>,
    pub version: Option<String>,
    pub detail: Option<String>,
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct BlastExternalBinaryPreflightReport {
    pub schema: String,
    pub blastn: ExternalBinaryPreflightProbe,
    pub makeblastdb: ExternalBinaryPreflightProbe,
}

pub fn blast_external_binary_preflight_report() -> BlastExternalBinaryPreflightReport {
    blast_external_binary_preflight_report_with_executables(
        resolve_tool_executable(BLASTN_ENV_BIN, DEFAULT_BLASTN_BIN),
        resolve_tool_executable(MAKEBLASTDB_ENV_BIN, DEFAULT_MAKEBLASTDB_BIN),
    )
}

fn blast_external_binary_preflight_report_with_executables(
    blastn_executable: String,
    makeblastdb_executable: String,
) -> BlastExternalBinaryPreflightReport {
    BlastExternalBinaryPreflightReport {
        schema: "gentle.blast_external_binary_preflight.v1".to_string(),
        blastn: probe_external_binary(BLASTN_ENV_BIN, "blastn", blastn_executable, &["-version"]),
        makeblastdb: probe_external_binary(
            MAKEBLASTDB_ENV_BIN,
            "makeblastdb",
            makeblastdb_executable,
            &["-version"],
        ),
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PreparedGenomeInspection {
    pub genome_id: String,
    pub install_dir: String,
    pub manifest_path: String,
    pub sequence_source_type: String,
    pub annotation_source_type: String,
    pub sequence_source: String,
    pub annotation_source: String,
    pub sequence_path: String,
    pub annotation_path: String,
    pub fasta_index_path: String,
    pub gene_index_path: String,
    pub blast_db_prefix: Option<String>,
    pub blast_index_ready: bool,
    pub sequence_sha1: Option<String>,
    pub annotation_sha1: Option<String>,
    pub sequence_present: bool,
    pub annotation_present: bool,
    pub fasta_index_ready: bool,
    pub gene_index_ready: bool,
    pub total_size_bytes: u64,
    pub installed_at_unix_ms: u128,
    #[serde(default)]
    pub cached_contig_count: usize,
    #[serde(default)]
    pub cached_total_span_bp: u64,
    #[serde(default)]
    pub cached_longest_contig: Option<String>,
    #[serde(default)]
    pub cached_longest_contig_bp: Option<u64>,
    #[serde(default)]
    pub cached_contig_preview: Vec<String>,
}

#[derive(Debug, Clone, Default)]
struct PreparedSequenceCacheSummary {
    contig_count: usize,
    total_span_bp: u64,
    longest_contig: Option<String>,
    longest_contig_bp: Option<u64>,
    contig_preview: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PrepareGenomeMode {
    PrepareOrReuse,
    ReindexCachedFiles,
    RefreshFromSources,
}

/// Result of resolving a requested genome id to a prepared install.
///
/// When the exact requested id is not prepared, a compatible prepared assembly
/// can be selected deterministically if and only if there is a unique match.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreparedGenomeResolution {
    pub requested_genome_id: String,
    pub resolved_genome_id: String,
    pub fallback_warning: Option<String>,
}

/// Policy controlling how prepared-genome compatibility fallback is handled.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PreparedGenomeFallbackPolicy {
    /// Require exact prepared genome id (no compatibility fallback).
    Off,
    /// Auto-fallback only when exactly one compatible prepared entry exists.
    SingleCompatible,
    /// Never auto-fallback; require explicit user/tool choice from options.
    AlwaysExplicit,
}

impl Default for PreparedGenomeFallbackPolicy {
    fn default() -> Self {
        Self::SingleCompatible
    }
}

/// Inspection payload used by GUI/CLI preflight for prepared-genome selection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreparedGenomeCompatibilityInspection {
    pub requested_genome_id: String,
    pub requested_catalog_key: String,
    pub requested_family: Option<String>,
    pub exact_prepared: bool,
    pub compatible_prepared_options: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlastHit {
    pub subject_id: String,
    pub identity_percent: f64,
    pub alignment_length: usize,
    pub mismatches: usize,
    pub gap_opens: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub subject_start: usize,
    pub subject_end: usize,
    pub evalue: f64,
    pub bit_score: f64,
    pub query_coverage_percent: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeBlastReport {
    pub genome_id: String,
    pub query_length: usize,
    pub max_hits: usize,
    pub task: String,
    pub blastn_executable: String,
    pub blast_db_prefix: String,
    pub command: Vec<String>,
    pub hit_count: usize,
    pub hits: Vec<BlastHit>,
    #[serde(default)]
    pub warnings: Vec<String>,
    pub stderr: String,
    #[serde(default)]
    pub options_override_json: Option<serde_json::Value>,
    #[serde(default)]
    pub effective_options_json: Option<serde_json::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrepareGenomeProgress {
    pub genome_id: String,
    pub phase: String,
    pub item: String,
    pub bytes_done: u64,
    pub bytes_total: Option<u64>,
    pub percent: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AnnotationParseIssue {
    pub line: usize,
    pub reason: String,
    pub context: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AnnotationParseReport {
    pub path: String,
    pub source_format: String,
    pub total_lines: usize,
    pub parsed_gene_records: usize,
    pub skipped_lines: usize,
    pub malformed_lines: usize,
    #[serde(default)]
    pub issue_examples: Vec<AnnotationParseIssue>,
    pub truncated_issue_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeGeneRecord {
    pub chromosome: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: Option<char>,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub biotype: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeTranscriptRecord {
    pub chromosome: String,
    pub transcript_id: String,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub strand: Option<char>,
    pub transcript_start_1based: usize,
    pub transcript_end_1based: usize,
    pub exons_1based: Vec<(usize, usize)>,
    #[serde(default)]
    pub cds_1based: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeChromosomeRecord {
    pub chromosome: String,
    pub length_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeSourcePlan {
    pub sequence_source: String,
    pub annotation_source: String,
    pub sequence_source_type: String,
    pub annotation_source_type: String,
    #[serde(default)]
    pub nucleotide_length_bp: Option<usize>,
    #[serde(default)]
    pub molecular_mass_da: Option<f64>,
    #[serde(default)]
    pub molecular_mass_source: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SourceType {
    Local,
    NcbiAssembly,
    GenbankAccession,
    RemoteHttp,
}

#[derive(Debug, Clone)]
struct SourceResolution {
    source: String,
    source_type: SourceType,
}

#[derive(Debug, Clone)]
struct FastaIndexEntry {
    length: u64,
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

#[derive(Debug, Clone, Default)]
struct BlastIndexOutcome {
    ready: bool,
    executable: Option<String>,
    warnings: Vec<String>,
}

#[derive(Debug, Clone, Default)]
pub struct GenomeCatalog {
    entries: HashMap<String, GenomeCatalogEntry>,
    catalog_base_dir: PathBuf,
}

impl GenomeCatalog {
    /// Load and validate a genome catalog JSON file.
    ///
    /// Validation enforces source consistency rules (local/remote/assembly/
    /// accession combinations) so downstream preparation logic can stay strict.
    pub fn from_json_file(path: &str) -> Result<Self, String> {
        let text = fs::read_to_string(path)
            .map_err(|e| format!("Could not read genome catalog '{path}': {e}"))?;
        let entries: HashMap<String, GenomeCatalogEntry> = serde_json::from_str(&text)
            .map_err(|e| format!("Could not parse genome catalog '{path}': {e}"))?;
        validate_catalog_entries(path, &entries)?;
        let base = Path::new(path)
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        Ok(Self {
            entries,
            catalog_base_dir: base,
        })
    }

    pub fn list_genomes(&self) -> Vec<String> {
        let mut names: Vec<String> = self.entries.keys().cloned().collect();
        names.sort_unstable();
        names
    }

    /// Resolve concrete sequence/annotation sources for one genome entry.
    ///
    /// The returned plan captures normalized source paths/URLs, inferred source
    /// types, and optional physical metadata (length/mass) including estimated
    /// molecular mass when only nucleotide length is available.
    pub fn source_plan(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<GenomeSourcePlan, String> {
        let entry = self.entry(genome_id)?;
        let sequence = self.resolve_source_with_type(
            genome_id,
            "sequence",
            entry.sequence_local.as_ref(),
            entry.sequence_remote.as_ref(),
            entry,
        )?;
        let annotation = self.resolve_source_with_type(
            genome_id,
            "annotation",
            entry.annotations_local.as_ref(),
            entry.annotations_remote.as_ref(),
            entry,
        )?;
        let nucleotide_length_bp = entry
            .nucleotide_length_bp
            .or_else(|| self.prepared_sequence_length_bp(genome_id, entry, cache_dir_override));
        let (molecular_mass_da, molecular_mass_source) =
            if let Some(mass_da) = entry.molecular_mass_da {
                (Some(mass_da), Some("catalog".to_string()))
            } else if let Some(length_bp) = nucleotide_length_bp {
                (
                    estimate_double_stranded_dna_mass_da(length_bp),
                    Some("estimated_from_nucleotide_length".to_string()),
                )
            } else {
                (None, None)
            };
        Ok(GenomeSourcePlan {
            sequence_source_type: source_type_label(sequence.source_type).to_string(),
            annotation_source_type: source_type_label(annotation.source_type).to_string(),
            sequence_source: sequence.source,
            annotation_source: annotation.source,
            nucleotide_length_bp,
            molecular_mass_da,
            molecular_mass_source,
        })
    }

    fn prepared_sequence_length_bp(
        &self,
        genome_id: &str,
        entry: &GenomeCatalogEntry,
        cache_dir_override: Option<&str>,
    ) -> Option<usize> {
        let install_dir = self.install_dir(genome_id, entry, cache_dir_override);
        let manifest_path = install_dir.join("manifest.json");
        if !manifest_path.exists() {
            return None;
        }
        let manifest = Self::load_manifest(&manifest_path).ok()?;
        let index = load_fasta_index(Path::new(&manifest.fasta_index_path)).ok()?;
        let total_bp = index
            .values()
            .fold(0u128, |acc, entry| acc.saturating_add(entry.length as u128));
        usize::try_from(total_bp).ok()
    }

    /// Return whether a prepared install appears complete for this genome id.
    pub fn is_prepared(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<bool, String> {
        let entry = self.entry(genome_id)?;
        let install_dir = self.install_dir(genome_id, entry, cache_dir_override);
        let manifest_path = install_dir.join("manifest.json");
        if !manifest_path.exists() {
            return Ok(false);
        }
        let manifest = Self::load_manifest(&manifest_path)?;
        if Self::validate_manifest_files(&manifest).is_err() {
            return Ok(false);
        }
        let gene_index_path = manifest
            .gene_index_path
            .as_ref()
            .map(PathBuf::from)
            .unwrap_or_else(|| install_dir.join("genes.json"));
        Ok(gene_index_path.exists())
    }

    /// Inspect prepared-install metadata and filesystem state.
    ///
    /// Returns `Ok(None)` when no install manifest exists yet.
    pub fn inspect_prepared_genome(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<Option<PreparedGenomeInspection>, String> {
        let entry = self.entry(genome_id)?;
        let install_dir = self.install_dir(genome_id, entry, cache_dir_override);
        let manifest_path = install_dir.join("manifest.json");
        if !manifest_path.exists() {
            return Ok(None);
        }
        let manifest = Self::load_manifest(&manifest_path)?;
        let gene_index_path = manifest
            .gene_index_path
            .as_ref()
            .map(PathBuf::from)
            .unwrap_or_else(|| install_dir.join("genes.json"));
        let sequence_path = PathBuf::from(&manifest.sequence_path);
        let annotation_path = PathBuf::from(&manifest.annotation_path);
        let fasta_index_path = PathBuf::from(&manifest.fasta_index_path);
        let blast_db_prefix_path = manifest
            .blast_db_prefix
            .as_deref()
            .map(PathBuf::from)
            .or_else(|| Some(default_blast_db_prefix(&install_dir)));
        let blast_db_prefix = manifest.blast_db_prefix.clone().or_else(|| {
            blast_db_prefix_path
                .as_ref()
                .map(|p| canonical_or_display(p))
        });
        let blast_index_files = blast_db_prefix_path
            .as_ref()
            .map(|p| collect_blast_index_files(p))
            .unwrap_or_default();

        let sequence_present = sequence_path.exists();
        let annotation_present = annotation_path.exists();
        let fasta_index_ready = fasta_index_path.exists();
        let gene_index_ready = gene_index_path.exists();
        let blast_index_ready = is_blast_index_ready(&blast_index_files);
        let cache_summary = if fasta_index_ready {
            load_fasta_index(&fasta_index_path)
                .map(|index| summarize_fasta_index(&index))
                .unwrap_or_default()
        } else {
            PreparedSequenceCacheSummary::default()
        };
        let total_size_bytes = [
            &sequence_path,
            &annotation_path,
            &fasta_index_path,
            &gene_index_path,
        ]
        .into_iter()
        .filter_map(|path| fs::metadata(path).ok())
        .map(|meta| meta.len())
        .sum::<u64>()
            + blast_index_files
                .iter()
                .filter_map(|path| fs::metadata(path).ok())
                .map(|meta| meta.len())
                .sum::<u64>();

        Ok(Some(PreparedGenomeInspection {
            genome_id: genome_id.to_string(),
            install_dir: canonical_or_display(&install_dir),
            manifest_path: canonical_or_display(&manifest_path),
            sequence_source_type: manifest.sequence_source_type.clone().unwrap_or_else(|| {
                classify_source_type_label(&manifest.sequence_source).to_string()
            }),
            annotation_source_type: manifest.annotation_source_type.clone().unwrap_or_else(|| {
                classify_source_type_label(&manifest.annotation_source).to_string()
            }),
            sequence_source: manifest.sequence_source.clone(),
            annotation_source: manifest.annotation_source.clone(),
            sequence_path: manifest.sequence_path.clone(),
            annotation_path: manifest.annotation_path.clone(),
            fasta_index_path: manifest.fasta_index_path.clone(),
            gene_index_path: canonical_or_display(&gene_index_path),
            blast_db_prefix,
            blast_index_ready,
            sequence_sha1: manifest.sequence_sha1.clone(),
            annotation_sha1: manifest.annotation_sha1.clone(),
            sequence_present,
            annotation_present,
            fasta_index_ready,
            gene_index_ready,
            total_size_bytes,
            installed_at_unix_ms: manifest.installed_at_unix_ms,
            cached_contig_count: cache_summary.contig_count,
            cached_total_span_bp: cache_summary.total_span_bp,
            cached_longest_contig: cache_summary.longest_contig,
            cached_longest_contig_bp: cache_summary.longest_contig_bp,
            cached_contig_preview: cache_summary.contig_preview,
        }))
    }

    /// Prepare sequence+annotation+indexes once using catalog/default cache.
    pub fn prepare_genome_once(&self, genome_id: &str) -> Result<PrepareGenomeReport, String> {
        self.prepare_genome_once_with_cache(genome_id, None)
    }

    /// Prepare sequence+annotation+indexes once with an explicit cache override.
    pub fn prepare_genome_once_with_cache(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<PrepareGenomeReport, String> {
        let mut noop = |_p: PrepareGenomeProgress| true;
        self.prepare_genome_once_with_progress(genome_id, cache_dir_override, &mut noop)
    }

    /// Rebuild indexes from already cached local sequence/annotation files
    /// without re-downloading sources.
    pub fn reindex_genome_once(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<PrepareGenomeReport, String> {
        let mut noop = |_p: PrepareGenomeProgress| true;
        self.reindex_genome_once_with_progress(genome_id, cache_dir_override, &mut noop)
    }

    /// Reinstall sequence+annotation+indexes from sources, replacing any stale
    /// prepared files already present in the install dir.
    pub fn reinstall_genome_once(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<PrepareGenomeReport, String> {
        let mut noop = |_p: PrepareGenomeProgress| true;
        self.reinstall_genome_once_with_progress(genome_id, cache_dir_override, &mut noop)
    }

    /// Prepare sequence+annotation+indexes while reporting progress.
    ///
    /// The callback is cooperative: returning `false` requests cancellation.
    pub fn prepare_genome_once_with_progress(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, String> {
        self.prepare_genome_once_with_progress_options(
            genome_id,
            cache_dir_override,
            PrepareGenomeMode::PrepareOrReuse,
            on_progress,
        )
    }

    /// Rebuild indexes from cached local sequence/annotation files while
    /// reporting progress.
    pub fn reindex_genome_once_with_progress(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, String> {
        self.prepare_genome_once_with_progress_options(
            genome_id,
            cache_dir_override,
            PrepareGenomeMode::ReindexCachedFiles,
            on_progress,
        )
    }

    /// Reinstall sequence+annotation+indexes while reporting progress.
    ///
    /// This forces a fresh source materialization and index rebuild instead of
    /// reusing prepared files already present in the install dir.
    pub fn reinstall_genome_once_with_progress(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, String> {
        self.prepare_genome_once_with_progress_options(
            genome_id,
            cache_dir_override,
            PrepareGenomeMode::RefreshFromSources,
            on_progress,
        )
    }

    fn prepare_genome_once_with_progress_options(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
        mode: PrepareGenomeMode,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, String> {
        let entry = self.entry(genome_id)?;
        let install_dir = self.install_dir(genome_id, entry, cache_dir_override);
        let manifest_path = install_dir.join("manifest.json");
        let sequence_resolution = self.resolve_source_with_type(
            genome_id,
            "sequence",
            entry.sequence_local.as_ref(),
            entry.sequence_remote.as_ref(),
            entry,
        )?;
        let annotation_resolution = self.resolve_source_with_type(
            genome_id,
            "annotation",
            entry.annotations_local.as_ref(),
            entry.annotations_remote.as_ref(),
            entry,
        )?;
        let sequence_source = sequence_resolution.source.clone();
        let annotation_source = annotation_resolution.source.clone();
        let reindex_from_cached_files = matches!(mode, PrepareGenomeMode::ReindexCachedFiles);
        let mut force_refresh_from_sources = matches!(mode, PrepareGenomeMode::RefreshFromSources);
        let mut cached_source_drift_warning: Option<String> = None;
        let existing_manifest = if manifest_path.exists() {
            Some(Self::load_manifest(&manifest_path)?)
        } else {
            None
        };

        if let Some(manifest) = existing_manifest.as_ref() {
            let sources_drifted = !sources_equivalent(&manifest.sequence_source, &sequence_source)
                || !sources_equivalent(&manifest.annotation_source, &annotation_source);
            if sources_drifted {
                if reindex_from_cached_files {
                    cached_source_drift_warning = Some(format!(
                        "Prepared cache sources differ from the current catalog entry for '{}'; reindex used the cached local files without deleting or re-downloading them.",
                        genome_id
                    ));
                } else {
                    force_refresh_from_sources = true;
                }
            }
            if !force_refresh_from_sources {
                Self::validate_manifest_files(&manifest)?;
            }
        } else if reindex_from_cached_files {
            return Err(format!(
                "Genome '{}' is not prepared locally; prepare it once before reindexing cached files",
                genome_id
            ));
        }

        if force_refresh_from_sources {
            reset_genome_install_dir(&install_dir)?;
        }
        fs::create_dir_all(&install_dir).map_err(|e| {
            format!(
                "Could not create genome cache dir '{}': {e}",
                install_dir.display()
            )
        })?;

        if let Some(mut manifest) = existing_manifest.clone() {
            if reindex_from_cached_files || force_refresh_from_sources {
                // Fall through to the rebuild path below.
            } else {
            let checksum_changed = ensure_manifest_checksums(&mut manifest)?;
            let mut warnings: Vec<String> = vec![];
            if let Some(warning) = cached_source_drift_warning.clone() {
                warnings.push(warning);
            }
            let mut annotation_parse_report: Option<AnnotationParseReport> = None;
            let gene_index_path = manifest
                .gene_index_path
                .as_ref()
                .map(PathBuf::from)
                .unwrap_or_else(|| install_dir.join("genes.json"));
            if !gene_index_path.exists() {
                let report = build_gene_index_file(
                    Path::new(&manifest.annotation_path),
                    &gene_index_path,
                    |done, total| {
                        on_progress(PrepareGenomeProgress {
                            genome_id: genome_id.to_string(),
                            phase: "index_genes".to_string(),
                            item: canonical_or_display(Path::new(&manifest.annotation_path)),
                            bytes_done: done,
                            bytes_total: total,
                            percent: total.and_then(|t| {
                                if t == 0 {
                                    None
                                } else {
                                    Some((done as f64 / t as f64) * 100.0)
                                }
                            }),
                        })
                    },
                )?;
                warnings.extend(summarize_annotation_parse_warnings(&report));
                annotation_parse_report = Some(report);
                manifest.gene_index_path = Some(canonical_or_display(&gene_index_path));
            }
            let fasta_index = load_fasta_index(Path::new(&manifest.fasta_index_path))?;
            let gene_records = load_gene_index_file(&gene_index_path)?;
            validate_annotation_gene_contigs_against_fasta_index(
                genome_id,
                &fasta_index,
                &gene_records,
                Path::new(&manifest.sequence_path),
                Path::new(&manifest.fasta_index_path),
                &gene_index_path,
            )?;
            let blast_prefix_path = manifest
                .blast_db_prefix
                .as_deref()
                .map(PathBuf::from)
                .unwrap_or_else(|| default_blast_db_prefix(&install_dir));
            forward_prepare_progress(
                on_progress,
                PrepareGenomeProgress {
                    genome_id: genome_id.to_string(),
                    phase: "index_blast".to_string(),
                    item: canonical_or_display(&blast_prefix_path),
                    bytes_done: 0,
                    bytes_total: None,
                    percent: None,
                },
            )?;
            let blast_outcome =
                ensure_blast_index(Path::new(&manifest.sequence_path), &blast_prefix_path);
            if !blast_outcome.warnings.is_empty() {
                warnings.extend(blast_outcome.warnings.clone());
            }
            forward_prepare_progress(
                on_progress,
                PrepareGenomeProgress {
                    genome_id: genome_id.to_string(),
                    phase: "index_blast".to_string(),
                    item: canonical_or_display(&blast_prefix_path),
                    bytes_done: if blast_outcome.ready { 1 } else { 0 },
                    bytes_total: Some(1),
                    percent: Some(if blast_outcome.ready { 100.0 } else { 0.0 }),
                },
            )?;
            let mut manifest_changed = checksum_changed;
            let blast_prefix = canonical_or_display(&blast_prefix_path);
            if manifest.blast_db_prefix.as_deref() != Some(blast_prefix.as_str()) {
                manifest.blast_db_prefix = Some(blast_prefix.clone());
                manifest_changed = true;
            }
            if manifest.blast_index_executable != blast_outcome.executable {
                manifest.blast_index_executable = blast_outcome.executable.clone();
                manifest_changed = true;
            }
            if blast_outcome.ready {
                if manifest.blast_indexed_at_unix_ms.is_none() {
                    manifest.blast_indexed_at_unix_ms = Some(now_unix_ms());
                    manifest_changed = true;
                }
            } else if manifest.blast_indexed_at_unix_ms.is_some() {
                manifest.blast_indexed_at_unix_ms = None;
                manifest_changed = true;
            }
            if manifest_changed {
                Self::write_manifest(&manifest_path, &manifest)?;
            }
            forward_prepare_progress(
                on_progress,
                PrepareGenomeProgress {
                    genome_id: genome_id.to_string(),
                    phase: "ready".to_string(),
                    item: "manifest".to_string(),
                    bytes_done: 0,
                    bytes_total: None,
                    percent: Some(100.0),
                },
            )?;
            let sequence_source_type = manifest.sequence_source_type.clone().unwrap_or_else(|| {
                classify_source_type_label(&manifest.sequence_source).to_string()
            });
            let annotation_source_type =
                manifest.annotation_source_type.clone().unwrap_or_else(|| {
                    classify_source_type_label(&manifest.annotation_source).to_string()
                });
            let cache_summary = summarize_fasta_index(&fasta_index);
            return Ok(PrepareGenomeReport {
                genome_id: genome_id.to_string(),
                reused_existing: true,
                sequence_path: manifest.sequence_path,
                annotation_path: manifest.annotation_path,
                fasta_index_path: manifest.fasta_index_path,
                sequence_source: Some(manifest.sequence_source.clone()),
                annotation_source: Some(manifest.annotation_source.clone()),
                sequence_source_type: Some(sequence_source_type),
                annotation_source_type: Some(annotation_source_type),
                blast_db_prefix: manifest.blast_db_prefix.clone(),
                blast_index_ready: blast_outcome.ready,
                blast_index_executable: blast_outcome.executable,
                cached_contig_count: cache_summary.contig_count,
                cached_total_span_bp: cache_summary.total_span_bp,
                cached_longest_contig: cache_summary.longest_contig,
                cached_longest_contig_bp: cache_summary.longest_contig_bp,
                cached_contig_preview: cache_summary.contig_preview,
                warnings,
                annotation_parse_report,
            });
            }
        }

        let existing_paths = if force_refresh_from_sources {
            None
        } else {
            existing_manifest.as_ref()
        };
        let manifest_sequence_source = if reindex_from_cached_files {
            existing_manifest
                .as_ref()
                .map(|manifest| manifest.sequence_source.clone())
                .unwrap_or_else(|| sequence_source.clone())
        } else {
            sequence_source.clone()
        };
        let manifest_annotation_source = if reindex_from_cached_files {
            existing_manifest
                .as_ref()
                .map(|manifest| manifest.annotation_source.clone())
                .unwrap_or_else(|| annotation_source.clone())
        } else {
            annotation_source.clone()
        };
        let manifest_sequence_source_type = if reindex_from_cached_files {
            existing_manifest
                .as_ref()
                .and_then(|manifest| manifest.sequence_source_type.clone())
                .unwrap_or_else(|| {
                    source_type_label(sequence_resolution.source_type).to_string()
                })
        } else {
            source_type_label(sequence_resolution.source_type).to_string()
        };
        let manifest_annotation_source_type = if reindex_from_cached_files {
            existing_manifest
                .as_ref()
                .and_then(|manifest| manifest.annotation_source_type.clone())
                .unwrap_or_else(|| {
                    source_type_label(annotation_resolution.source_type).to_string()
                })
        } else {
            source_type_label(annotation_resolution.source_type).to_string()
        };
        let sequence_path = existing_paths
            .map(|manifest| PathBuf::from(&manifest.sequence_path))
            .unwrap_or_else(|| install_dir.join("sequence.fa"));
        let annotation_ext = infer_annotation_extension(&manifest_annotation_source);
        let annotation_path = existing_paths
            .map(|manifest| PathBuf::from(&manifest.annotation_path))
            .unwrap_or_else(|| install_dir.join(format!("annotation.{annotation_ext}")));
        let fasta_index_path = existing_paths
            .map(|manifest| PathBuf::from(&manifest.fasta_index_path))
            .unwrap_or_else(|| install_dir.join("sequence.fa.fai"));
        let gene_index_path = existing_paths
            .and_then(|manifest| manifest.gene_index_path.as_ref().map(PathBuf::from))
            .unwrap_or_else(|| install_dir.join("genes.json"));
        let blast_prefix_path = existing_paths
            .and_then(|manifest| manifest.blast_db_prefix.as_deref().map(PathBuf::from))
            .unwrap_or_else(|| default_blast_db_prefix(&install_dir));

        if reindex_from_cached_files {
            forward_prepare_progress(
                on_progress,
                PrepareGenomeProgress {
                    genome_id: genome_id.to_string(),
                    phase: "reset_indexes".to_string(),
                    item: install_dir.display().to_string(),
                    bytes_done: 0,
                    bytes_total: None,
                    percent: Some(0.0),
                },
            )?;
            reset_prepared_genome_index_artifacts(
                &fasta_index_path,
                &gene_index_path,
                &blast_prefix_path,
            )?;
        }

        if reindex_from_cached_files {
            if !non_empty_regular_file_exists(&sequence_path) {
                return Err(format!(
                    "Cannot reindex genome '{}' from cached files because prepared sequence '{}' is missing. Use source refresh/re-download instead.",
                    genome_id,
                    sequence_path.display()
                ));
            }
            let bytes = fs::metadata(&sequence_path)
                .map(|meta| meta.len())
                .unwrap_or(0);
            forward_prepare_progress(
                on_progress,
                PrepareGenomeProgress {
                    genome_id: genome_id.to_string(),
                    phase: "reuse_sequence".to_string(),
                    item: canonical_or_display(&sequence_path),
                    bytes_done: bytes,
                    bytes_total: Some(bytes),
                    percent: Some(100.0),
                },
            )?;
        } else if !force_refresh_from_sources && non_empty_regular_file_exists(&sequence_path) {
            let bytes = fs::metadata(&sequence_path)
                .map(|meta| meta.len())
                .unwrap_or(0);
            forward_prepare_progress(
                on_progress,
                PrepareGenomeProgress {
                    genome_id: genome_id.to_string(),
                    phase: "reuse_sequence".to_string(),
                    item: canonical_or_display(&sequence_path),
                    bytes_done: bytes,
                    bytes_total: Some(bytes),
                    percent: Some(100.0),
                },
            )?;
        } else {
            materialize_source_with_progress(&sequence_source, &sequence_path, |done, total| {
                on_progress(PrepareGenomeProgress {
                    genome_id: genome_id.to_string(),
                    phase: "download_sequence".to_string(),
                    item: sequence_source.clone(),
                    bytes_done: done,
                    bytes_total: total,
                    percent: total.and_then(|t| {
                        if t == 0 {
                            None
                        } else {
                            Some((done as f64 / t as f64) * 100.0)
                        }
                    }),
                })
            })?;
        }
        if reindex_from_cached_files {
            if !non_empty_regular_file_exists(&annotation_path) {
                return Err(format!(
                    "Cannot reindex genome '{}' from cached files because prepared annotation '{}' is missing. Use source refresh/re-download instead.",
                    genome_id,
                    annotation_path.display()
                ));
            }
            let bytes = fs::metadata(&annotation_path)
                .map(|meta| meta.len())
                .unwrap_or(0);
            forward_prepare_progress(
                on_progress,
                PrepareGenomeProgress {
                    genome_id: genome_id.to_string(),
                    phase: "reuse_annotation".to_string(),
                    item: canonical_or_display(&annotation_path),
                    bytes_done: bytes,
                    bytes_total: Some(bytes),
                    percent: Some(100.0),
                },
            )?;
        } else if !force_refresh_from_sources && non_empty_regular_file_exists(&annotation_path) {
            let bytes = fs::metadata(&annotation_path)
                .map(|meta| meta.len())
                .unwrap_or(0);
            forward_prepare_progress(
                on_progress,
                PrepareGenomeProgress {
                    genome_id: genome_id.to_string(),
                    phase: "reuse_annotation".to_string(),
                    item: canonical_or_display(&annotation_path),
                    bytes_done: bytes,
                    bytes_total: Some(bytes),
                    percent: Some(100.0),
                },
            )?;
        } else {
            materialize_source_with_progress(
                &annotation_source,
                &annotation_path,
                |done, total| {
                    on_progress(PrepareGenomeProgress {
                        genome_id: genome_id.to_string(),
                        phase: "download_annotation".to_string(),
                        item: annotation_source.clone(),
                        bytes_done: done,
                        bytes_total: total,
                        percent: total.and_then(|t| {
                            if t == 0 {
                                None
                            } else {
                                Some((done as f64 / t as f64) * 100.0)
                            }
                        }),
                    })
                },
            )?;
        }
        build_fasta_index_with_progress(&sequence_path, &fasta_index_path, |done, total| {
            on_progress(PrepareGenomeProgress {
                genome_id: genome_id.to_string(),
                phase: "index_fasta".to_string(),
                item: canonical_or_display(&sequence_path),
                bytes_done: done,
                bytes_total: total,
                percent: total.and_then(|t| {
                    if t == 0 {
                        None
                    } else {
                        Some((done as f64 / t as f64) * 100.0)
                    }
                }),
            })
        })?;
        let annotation_parse_report =
            if !reindex_from_cached_files
                && !force_refresh_from_sources
                && non_empty_regular_file_exists(&gene_index_path)
            {
                let bytes = fs::metadata(&gene_index_path)
                    .map(|meta| meta.len())
                    .unwrap_or(0);
                forward_prepare_progress(
                    on_progress,
                    PrepareGenomeProgress {
                        genome_id: genome_id.to_string(),
                        phase: "reuse_gene_index".to_string(),
                        item: canonical_or_display(&gene_index_path),
                        bytes_done: bytes,
                        bytes_total: Some(bytes),
                        percent: Some(100.0),
                    },
                )?;
                None
            } else {
                Some(build_gene_index_file(
                    &annotation_path,
                    &gene_index_path,
                    |done, total| {
                        on_progress(PrepareGenomeProgress {
                            genome_id: genome_id.to_string(),
                            phase: "index_genes".to_string(),
                            item: canonical_or_display(&annotation_path),
                            bytes_done: done,
                            bytes_total: total,
                            percent: total.and_then(|t| {
                                if t == 0 {
                                    None
                                } else {
                                    Some((done as f64 / t as f64) * 100.0)
                                }
                            }),
                        })
                    },
                )?)
            };
        let fasta_index = load_fasta_index(&fasta_index_path)?;
        let gene_records = load_gene_index_file(&gene_index_path)?;
        validate_annotation_gene_contigs_against_fasta_index(
            genome_id,
            &fasta_index,
            &gene_records,
            &sequence_path,
            &fasta_index_path,
            &gene_index_path,
        )?;
        forward_prepare_progress(
            on_progress,
            PrepareGenomeProgress {
                genome_id: genome_id.to_string(),
                phase: "index_blast".to_string(),
                item: canonical_or_display(&blast_prefix_path),
                bytes_done: 0,
                bytes_total: None,
                percent: None,
            },
        )?;
        let blast_outcome = ensure_blast_index(&sequence_path, &blast_prefix_path);
        forward_prepare_progress(
            on_progress,
            PrepareGenomeProgress {
                genome_id: genome_id.to_string(),
                phase: "index_blast".to_string(),
                item: canonical_or_display(&blast_prefix_path),
                bytes_done: if blast_outcome.ready { 1 } else { 0 },
                bytes_total: Some(1),
                percent: Some(if blast_outcome.ready { 100.0 } else { 0.0 }),
            },
        )?;

        let manifest = GenomeInstallManifest {
            genome_id: genome_id.to_string(),
            sequence_source: manifest_sequence_source,
            annotation_source: manifest_annotation_source,
            sequence_source_type: Some(manifest_sequence_source_type),
            annotation_source_type: Some(manifest_annotation_source_type),
            sequence_sha1: Some(compute_file_sha1(&sequence_path)?),
            annotation_sha1: Some(compute_file_sha1(&annotation_path)?),
            sequence_path: canonical_or_display(&sequence_path),
            annotation_path: canonical_or_display(&annotation_path),
            fasta_index_path: canonical_or_display(&fasta_index_path),
            gene_index_path: Some(canonical_or_display(&gene_index_path)),
            blast_db_prefix: Some(canonical_or_display(&blast_prefix_path)),
            blast_index_executable: blast_outcome.executable.clone(),
            blast_indexed_at_unix_ms: blast_outcome.ready.then_some(now_unix_ms()),
            installed_at_unix_ms: now_unix_ms(),
        };
        Self::write_manifest(&manifest_path, &manifest)?;

        forward_prepare_progress(
            on_progress,
            PrepareGenomeProgress {
                genome_id: genome_id.to_string(),
                phase: "ready".to_string(),
                item: "manifest".to_string(),
                bytes_done: 0,
                bytes_total: None,
                percent: Some(100.0),
            },
        )?;

        let mut warnings = blast_outcome.warnings;
        if let Some(warning) = cached_source_drift_warning {
            warnings.push(warning);
        }
        if let Some(report) = annotation_parse_report.as_ref() {
            warnings.extend(summarize_annotation_parse_warnings(report));
        }
        let cache_summary = summarize_fasta_index(&fasta_index);

        Ok(PrepareGenomeReport {
            genome_id: genome_id.to_string(),
            reused_existing: reindex_from_cached_files,
            sequence_path: manifest.sequence_path,
            annotation_path: manifest.annotation_path,
            fasta_index_path: manifest.fasta_index_path,
            sequence_source: Some(manifest.sequence_source.clone()),
            annotation_source: Some(manifest.annotation_source.clone()),
            sequence_source_type: manifest.sequence_source_type.clone(),
            annotation_source_type: manifest.annotation_source_type.clone(),
            blast_db_prefix: manifest.blast_db_prefix.clone(),
            blast_index_ready: blast_outcome.ready,
            blast_index_executable: blast_outcome.executable,
            cached_contig_count: cache_summary.contig_count,
            cached_total_span_bp: cache_summary.total_span_bp,
            cached_longest_contig: cache_summary.longest_contig,
            cached_longest_contig_bp: cache_summary.longest_contig_bp,
            cached_contig_preview: cache_summary.contig_preview,
            warnings,
            annotation_parse_report,
        })
    }

    /// Extract a 1-based inclusive genomic interval from prepared sequence.
    pub fn get_sequence_region(
        &self,
        genome_id: &str,
        chromosome: &str,
        start_1based: usize,
        end_1based: usize,
    ) -> Result<String, String> {
        self.get_sequence_region_with_cache(genome_id, chromosome, start_1based, end_1based, None)
    }

    /// Extract a 1-based inclusive genomic interval with explicit cache
    /// override.
    ///
    /// Chromosome aliases with/without `chr` prefix are tried automatically.
    pub fn get_sequence_region_with_cache(
        &self,
        genome_id: &str,
        chromosome: &str,
        start_1based: usize,
        end_1based: usize,
        cache_dir_override: Option<&str>,
    ) -> Result<String, String> {
        if start_1based == 0 {
            return Err("Coordinates must be 1-based (start >= 1)".to_string());
        }
        if end_1based < start_1based {
            return Err(format!(
                "Invalid interval: start ({start_1based}) is greater than end ({end_1based})"
            ));
        }
        let prepared = self.resolve_prepared_genome_id(genome_id, cache_dir_override)?;
        let resolved_genome_id = prepared.resolved_genome_id;
        let entry = self.entry(&resolved_genome_id)?;
        let manifest_path = self
            .install_dir(&resolved_genome_id, entry, cache_dir_override)
            .join("manifest.json");

        let manifest = Self::load_manifest(&manifest_path)?;
        Self::validate_manifest_files(&manifest)?;

        let index = load_fasta_index(Path::new(&manifest.fasta_index_path))?;
        let mut chr_candidates = Vec::new();
        chr_candidates.push(chromosome.to_string());
        if let Some(trimmed) = chromosome.strip_prefix("chr") {
            chr_candidates.push(trimmed.to_string());
        } else {
            chr_candidates.push(format!("chr{chromosome}"));
        }

        let index_entry = chr_candidates
            .iter()
            .find_map(|name| index.get(name))
            .ok_or_else(|| {
                let mut tried_aliases = chr_candidates.clone();
                tried_aliases.sort_unstable();
                tried_aliases.dedup();

                let mut available: Vec<String> = index.keys().cloned().collect();
                available.sort_unstable();
                let query_lower = chromosome.trim().to_ascii_lowercase();
                let query_normalized = normalize_chromosome_token(chromosome);
                let mut suggestions = available
                    .iter()
                    .filter_map(|candidate| {
                        let candidate_lower = candidate.to_ascii_lowercase();
                        let candidate_normalized = normalize_chromosome_token(candidate);
                        let query_numeric = !query_normalized.is_empty()
                            && query_normalized.chars().all(|ch| ch.is_ascii_digit());
                        let candidate_numeric = !candidate_normalized.is_empty()
                            && candidate_normalized.chars().all(|ch| ch.is_ascii_digit());
                        let allow_shorter_partial = !(query_numeric && candidate_numeric);
                        let score = if candidate_normalized == query_normalized {
                            0u8
                        } else if !query_normalized.is_empty()
                            && (candidate_normalized.starts_with(&query_normalized)
                                || (allow_shorter_partial
                                    && query_normalized.starts_with(&candidate_normalized)))
                        {
                            1u8
                        } else if !query_lower.is_empty()
                            && (candidate_lower.contains(&query_lower)
                                || (allow_shorter_partial
                                    && query_lower.contains(&candidate_lower)))
                        {
                            2u8
                        } else {
                            return None;
                        };
                        Some((score, candidate.clone()))
                    })
                    .collect::<Vec<_>>();
                suggestions.sort_by(|left, right| {
                    left.0
                        .cmp(&right.0)
                        .then(left.1.len().cmp(&right.1.len()))
                        .then(left.1.cmp(&right.1))
                });
                let mut suggested_names = Vec::new();
                for (_, name) in suggestions {
                    if suggested_names.iter().any(|v| v == &name) {
                        continue;
                    }
                    suggested_names.push(name);
                    if suggested_names.len() >= 8 {
                        break;
                    }
                }
                let preview_limit = 12usize;
                let preview_items = available
                    .iter()
                    .take(preview_limit)
                    .cloned()
                    .collect::<Vec<_>>();
                let mut preview = preview_items.join(", ");
                let hidden = available.len().saturating_sub(preview_items.len());
                if hidden > 0 {
                    preview.push_str(&format!(", ... (+{hidden} more)"));
                }
                if preview.is_empty() {
                    preview = "(none)".to_string();
                }
                let case_hint = available
                    .iter()
                    .find(|name| name.eq_ignore_ascii_case(chromosome))
                    .map(|name| format!(" Case-sensitive match exists as '{name}'."))
                    .unwrap_or_default();
                let suggestions_hint = if suggested_names.is_empty() {
                    String::new()
                } else {
                    format!(
                        " Suggested matching contigs: {}.",
                        suggested_names.join(", ")
                    )
                };
                format!(
                    "Chromosome/contig '{}' not found in genome '{}'. Tried aliases: {}. \
Available contigs ({}): {}. This can happen when prepared sequence/annotation naming differs \
or cache contents are stale; re-run PrepareGenome for this genome/cache. Prepared sequence='{}'; \
FASTA index='{}'.{}{}",
                    chromosome,
                    resolved_genome_id,
                    tried_aliases.join(", "),
                    available.len(),
                    preview,
                    manifest.sequence_path,
                    manifest.fasta_index_path,
                    case_hint,
                    suggestions_hint
                )
            })?;

        let end_u64 = end_1based as u64;
        if end_u64 > index_entry.length {
            return Err(format!(
                "Requested end {} exceeds chromosome length {}",
                end_1based, index_entry.length
            ));
        }

        read_fasta_slice(
            Path::new(&manifest.sequence_path),
            index_entry,
            start_1based as u64,
            end_1based as u64,
        )
    }

    /// List prepared chromosome/contig lengths, sorted by descending size.
    pub fn list_chromosome_lengths(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<Vec<GenomeChromosomeRecord>, String> {
        let prepared = self.resolve_prepared_genome_id(genome_id, cache_dir_override)?;
        let resolved_genome_id = prepared.resolved_genome_id;
        let entry = self.entry(&resolved_genome_id)?;
        let manifest_path = self
            .install_dir(&resolved_genome_id, entry, cache_dir_override)
            .join("manifest.json");
        let manifest = Self::load_manifest(&manifest_path)?;
        Self::validate_manifest_files(&manifest)?;
        let index = load_fasta_index(Path::new(&manifest.fasta_index_path))?;
        let mut chromosomes = index
            .into_iter()
            .map(|(chromosome, entry)| {
                let length_bp = usize::try_from(entry.length).map_err(|_| {
                    format!(
                        "Chromosome/contig '{}' length exceeds supported platform range",
                        chromosome
                    )
                })?;
                Ok(GenomeChromosomeRecord {
                    chromosome,
                    length_bp,
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        chromosomes.sort_by(|a, b| {
            b.length_bp
                .cmp(&a.length_bp)
                .then(a.chromosome.cmp(&b.chromosome))
        });
        Ok(chromosomes)
    }

    /// List gene-like annotation records from the prepared gene index.
    ///
    /// If the index is missing but manifest files are valid, it is rebuilt from
    /// the prepared annotation source.
    pub fn list_gene_regions(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<Vec<GenomeGeneRecord>, String> {
        let prepared = self.resolve_prepared_genome_id(genome_id, cache_dir_override)?;
        let resolved_genome_id = prepared.resolved_genome_id;
        let entry = self.entry(&resolved_genome_id)?;
        let install_dir = self.install_dir(&resolved_genome_id, entry, cache_dir_override);
        let manifest_path = install_dir.join("manifest.json");

        let mut manifest = Self::load_manifest(&manifest_path)?;
        Self::validate_manifest_files(&manifest)?;
        let gene_index_path = manifest
            .gene_index_path
            .as_ref()
            .map(PathBuf::from)
            .unwrap_or_else(|| install_dir.join("genes.json"));
        if !gene_index_path.exists() {
            let _ = build_gene_index_file(
                Path::new(&manifest.annotation_path),
                &gene_index_path,
                |_done, _total| true,
            )?;
            manifest.gene_index_path = Some(canonical_or_display(&gene_index_path));
            Self::write_manifest(&manifest_path, &manifest)?;
        }
        load_gene_index_file(&gene_index_path)
    }

    /// List transcript/exon records overlapping one extracted gene interval.
    ///
    /// This is used to project annotation-backed transcript structures into
    /// extracted gene slices (for example for isoform architecture views).
    /// For non-tabular annotation sources (GenBank/XML), this currently returns
    /// an empty list.
    pub fn list_gene_transcript_records(
        &self,
        genome_id: &str,
        chromosome: &str,
        start_1based: usize,
        end_1based: usize,
        gene_id_hint: Option<&str>,
        gene_name_hint: Option<&str>,
        cache_dir_override: Option<&str>,
    ) -> Result<Vec<GenomeTranscriptRecord>, String> {
        if start_1based == 0 || end_1based < start_1based {
            return Err(format!(
                "Invalid gene interval {}:{}-{}",
                chromosome, start_1based, end_1based
            ));
        }

        let prepared = self.resolve_prepared_genome_id(genome_id, cache_dir_override)?;
        let resolved_genome_id = prepared.resolved_genome_id;
        let entry = self.entry(&resolved_genome_id)?;
        let install_dir = self.install_dir(&resolved_genome_id, entry, cache_dir_override);
        let manifest_path = install_dir.join("manifest.json");

        let manifest = Self::load_manifest(&manifest_path)?;
        Self::validate_manifest_files(&manifest)?;
        let annotation_path = PathBuf::from(&manifest.annotation_path);
        if !annotation_path.exists() {
            return Err(format!(
                "Prepared annotation file '{}' is missing",
                annotation_path.display()
            ));
        }
        if is_genbank_annotation_path(&annotation_path) || is_xml_annotation_path(&annotation_path)
        {
            return Ok(vec![]);
        }

        #[derive(Default)]
        struct TranscriptAccum {
            chromosome: String,
            transcript_id: String,
            gene_id: Option<String>,
            gene_name: Option<String>,
            strand: Option<char>,
            transcript_start_1based: Option<usize>,
            transcript_end_1based: Option<usize>,
            exons_1based: Vec<(usize, usize)>,
            cds_1based: Vec<(usize, usize)>,
        }

        let file = File::open(&annotation_path).map_err(|e| {
            format!(
                "Could not open annotation file '{}': {e}",
                annotation_path.display()
            )
        })?;
        let reader = BufReader::new(file);
        let gene_id_hint_norm = gene_id_hint
            .map(|v| normalize_gene_match_token(v))
            .filter(|v| !v.is_empty());
        let gene_name_hint_norm = gene_name_hint
            .map(|v| normalize_gene_match_token(v))
            .filter(|v| !v.is_empty());

        let mut transcripts: HashMap<String, TranscriptAccum> = HashMap::new();
        for line_result in reader.lines() {
            let line = line_result.map_err(|e| {
                format!(
                    "Could not read annotation file '{}': {e}",
                    annotation_path.display()
                )
            })?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = trimmed.split('\t').collect();
            if cols.len() < 9 {
                continue;
            }
            if !chromosome_names_match(chromosome, cols[0]) {
                continue;
            }
            let feature_kind = cols[2].trim().to_ascii_lowercase();
            if !matches!(
                feature_kind.as_str(),
                "transcript" | "mrna" | "exon" | "cds"
            ) {
                continue;
            }
            let Some(row_start) = parse_annotation_coordinate(cols[3]) else {
                continue;
            };
            let Some(row_end) = parse_annotation_coordinate(cols[4]) else {
                continue;
            };
            if row_end < row_start {
                continue;
            }
            if row_end < start_1based || row_start > end_1based {
                continue;
            }
            let attrs = parse_annotation_attributes(cols[8]);
            let row_gene_id = pick_gene_id(&attrs);
            let row_gene_name = pick_gene_name(&attrs);
            let matches_gene_id = gene_id_hint_norm
                .as_ref()
                .map(|hint| {
                    row_gene_id
                        .as_deref()
                        .map(|value| normalize_gene_match_token(value) == *hint)
                        .unwrap_or(false)
                })
                .unwrap_or(true);
            let matches_gene_name = gene_name_hint_norm
                .as_ref()
                .map(|hint| {
                    row_gene_name
                        .as_deref()
                        .map(|value| normalize_gene_match_token(value) == *hint)
                        .unwrap_or(false)
                })
                .unwrap_or(true);
            if !matches_gene_id && !matches_gene_name {
                continue;
            }

            let transcript_id =
                pick_annotation_attribute(&attrs, &["transcript_id", "transcript", "id", "name"])
                    .map(|v| normalize_transcript_id(&v))
                    .filter(|v| !v.is_empty())
                    .or_else(|| {
                        if feature_kind == "transcript" || feature_kind == "mrna" {
                            Some(format!(
                                "{}:{}-{}",
                                row_gene_id
                                    .as_deref()
                                    .or(row_gene_name.as_deref())
                                    .unwrap_or("transcript"),
                                row_start,
                                row_end
                            ))
                        } else {
                            None
                        }
                    });
            let Some(transcript_id) = transcript_id else {
                continue;
            };
            let strand = cols[6].chars().next().and_then(|c| match c {
                '+' | '-' => Some(c),
                _ => None,
            });
            let entry =
                transcripts
                    .entry(transcript_id.clone())
                    .or_insert_with(|| TranscriptAccum {
                        chromosome: cols[0].trim().to_string(),
                        transcript_id: transcript_id.clone(),
                        gene_id: row_gene_id.clone(),
                        gene_name: row_gene_name.clone(),
                        strand,
                        transcript_start_1based: None,
                        transcript_end_1based: None,
                        exons_1based: vec![],
                        cds_1based: vec![],
                    });
            if entry.gene_id.is_none() {
                entry.gene_id = row_gene_id.clone();
            }
            if entry.gene_name.is_none() {
                entry.gene_name = row_gene_name.clone();
            }
            if entry.strand.is_none() {
                entry.strand = strand;
            }
            if feature_kind == "transcript" || feature_kind == "mrna" {
                entry.transcript_start_1based = Some(row_start);
                entry.transcript_end_1based = Some(row_end);
            } else if feature_kind == "exon" {
                entry.exons_1based.push((row_start, row_end));
            } else {
                entry.cds_1based.push((row_start, row_end));
            }
        }

        let mut out: Vec<GenomeTranscriptRecord> = transcripts
            .into_values()
            .filter_map(|mut tx| {
                tx.exons_1based
                    .sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
                tx.exons_1based.dedup();
                let mut clipped_exons: Vec<(usize, usize)> = tx
                    .exons_1based
                    .into_iter()
                    .filter_map(|(start, end)| {
                        let clipped_start = start.max(start_1based);
                        let clipped_end = end.min(end_1based);
                        (clipped_end >= clipped_start).then_some((clipped_start, clipped_end))
                    })
                    .collect();
                clipped_exons.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
                clipped_exons.dedup();
                let mut clipped_cds: Vec<(usize, usize)> = tx
                    .cds_1based
                    .into_iter()
                    .filter_map(|(start, end)| {
                        let clipped_start = start.max(start_1based);
                        let clipped_end = end.min(end_1based);
                        (clipped_end >= clipped_start).then_some((clipped_start, clipped_end))
                    })
                    .collect();
                clipped_cds.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
                clipped_cds.dedup();

                let transcript_start = tx
                    .transcript_start_1based
                    .or_else(|| clipped_exons.first().map(|(start, _)| *start))?;
                let transcript_end = tx
                    .transcript_end_1based
                    .or_else(|| clipped_exons.last().map(|(_, end)| *end))?;
                let clipped_start = transcript_start.max(start_1based);
                let clipped_end = transcript_end.min(end_1based);
                if clipped_end < clipped_start {
                    return None;
                }
                if clipped_exons.is_empty() {
                    clipped_exons.push((clipped_start, clipped_end));
                }
                Some(GenomeTranscriptRecord {
                    chromosome: tx.chromosome,
                    transcript_id: tx.transcript_id,
                    gene_id: tx.gene_id,
                    gene_name: tx.gene_name,
                    strand: tx.strand,
                    transcript_start_1based: clipped_start,
                    transcript_end_1based: clipped_end,
                    exons_1based: clipped_exons,
                    cds_1based: clipped_cds,
                })
            })
            .collect();
        out.sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
        Ok(out)
    }

    /// Run BLASTN of query sequence against the prepared genome index.
    pub fn blast_sequence(
        &self,
        genome_id: &str,
        query_sequence: &str,
        max_hits: usize,
        task: Option<&str>,
    ) -> Result<GenomeBlastReport, String> {
        self.blast_sequence_with_cache(genome_id, query_sequence, max_hits, task, None)
    }

    pub fn blast_sequence_with_cache(
        &self,
        genome_id: &str,
        query_sequence: &str,
        max_hits: usize,
        task: Option<&str>,
        cache_dir_override: Option<&str>,
    ) -> Result<GenomeBlastReport, String> {
        let mut never_cancel = || false;
        self.blast_sequence_with_cache_and_cancel(
            genome_id,
            query_sequence,
            max_hits,
            task,
            cache_dir_override,
            &mut never_cancel,
        )
    }

    pub fn blast_sequence_with_cache_and_cancel(
        &self,
        genome_id: &str,
        query_sequence: &str,
        max_hits: usize,
        task: Option<&str>,
        cache_dir_override: Option<&str>,
        should_cancel: &mut dyn FnMut() -> bool,
    ) -> Result<GenomeBlastReport, String> {
        if should_cancel() {
            return Err(blast_cancelled_error("before blast start"));
        }
        if max_hits == 0 {
            return Err("BLAST search requires max_hits >= 1".to_string());
        }
        let task = normalize_blast_task(task)?;
        let query = normalize_blast_query_sequence(query_sequence)?;

        let prepared = self.resolve_prepared_genome_id(genome_id, cache_dir_override)?;
        let resolved_genome_id = prepared.resolved_genome_id;
        let mut additional_warnings: Vec<String> = vec![];
        if let Some(warning) = prepared.fallback_warning {
            additional_warnings.push(warning);
        }
        let entry = self.entry(&resolved_genome_id)?;
        let install_dir = self.install_dir(&resolved_genome_id, entry, cache_dir_override);
        let manifest_path = install_dir.join("manifest.json");
        let mut manifest = Self::load_manifest(&manifest_path)?;
        Self::validate_manifest_files(&manifest)?;

        let blast_prefix_path = manifest
            .blast_db_prefix
            .as_deref()
            .map(PathBuf::from)
            .unwrap_or_else(|| default_blast_db_prefix(&install_dir));
        let mut blast_outcome =
            ensure_blast_index(Path::new(&manifest.sequence_path), &blast_prefix_path);
        let mut manifest_changed = false;
        let blast_prefix = canonical_or_display(&blast_prefix_path);
        if manifest.blast_db_prefix.as_deref() != Some(blast_prefix.as_str()) {
            manifest.blast_db_prefix = Some(blast_prefix.clone());
            manifest_changed = true;
        }
        if manifest.blast_index_executable != blast_outcome.executable {
            manifest.blast_index_executable = blast_outcome.executable.clone();
            manifest_changed = true;
        }
        if blast_outcome.ready {
            if manifest.blast_indexed_at_unix_ms.is_none() {
                manifest.blast_indexed_at_unix_ms = Some(now_unix_ms());
                manifest_changed = true;
            }
        } else if manifest.blast_indexed_at_unix_ms.is_some() {
            manifest.blast_indexed_at_unix_ms = None;
            manifest_changed = true;
        }
        if manifest_changed {
            Self::write_manifest(&manifest_path, &manifest)?;
        }
        if !blast_outcome.ready {
            let detail = if blast_outcome.warnings.is_empty() {
                "no index files found and makeblastdb could not build one".to_string()
            } else {
                blast_outcome.warnings.join(" | ")
            };
            return Err(format!(
                "BLAST index for genome '{}' is not ready (db='{}'): {}",
                resolved_genome_id, blast_prefix, detail
            ));
        }
        if should_cancel() {
            return Err(blast_cancelled_error("before blast launch"));
        }

        let blastn_executable = resolve_tool_executable(BLASTN_ENV_BIN, DEFAULT_BLASTN_BIN);
        let mut query_file = NamedTempFile::new_in(&install_dir).map_err(|e| {
            format!(
                "Could not create temporary BLAST query file in '{}': {e}",
                install_dir.display()
            )
        })?;
        write!(query_file, ">query\n{}\n", query)
            .map_err(|e| format!("Could not write temporary BLAST query file: {e}"))?;
        query_file
            .flush()
            .map_err(|e| format!("Could not flush temporary BLAST query file: {e}"))?;
        let query_path = canonical_or_display(query_file.path());

        let args = vec![
            "-db".to_string(),
            blast_prefix.clone(),
            "-query".to_string(),
            query_path,
            "-task".to_string(),
            task.clone(),
            "-outfmt".to_string(),
            BLASTN_OUTFMT_FIELDS.to_string(),
            "-max_target_seqs".to_string(),
            max_hits.to_string(),
        ];
        let mut child = Command::new(&blastn_executable)
            .args(&args)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| {
                if e.kind() == std::io::ErrorKind::NotFound {
                    format!(
                        "Could not find blastn executable '{}'. Install BLAST+ or set {}",
                        blastn_executable, BLASTN_ENV_BIN
                    )
                } else {
                    format!(
                        "Could not run blastn executable '{}' with args [{}]: {}",
                        blastn_executable,
                        args.join(" "),
                        e
                    )
                }
            })?;
        let output = wait_for_child_output_with_cancel(
            &mut child,
            &blastn_executable,
            &args,
            should_cancel,
        )?;
        let stdout = String::from_utf8_lossy(&output.stdout).to_string();
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
        if !output.status.success() {
            return Err(format!(
                "blastn search failed: {} {} (status={:?}, stderr='{}')",
                blastn_executable,
                args.join(" "),
                output.status.code(),
                stderr.trim()
            ));
        }

        let (hits, parse_warnings) = parse_blastn_tabular_hits(&stdout);
        if !parse_warnings.is_empty() {
            blast_outcome.warnings.extend(parse_warnings);
        }
        if !additional_warnings.is_empty() {
            blast_outcome.warnings.extend(additional_warnings);
        }
        Ok(GenomeBlastReport {
            genome_id: resolved_genome_id,
            query_length: query.len(),
            max_hits,
            task,
            blastn_executable,
            blast_db_prefix: blast_prefix,
            command: args,
            hit_count: hits.len(),
            hits,
            warnings: blast_outcome.warnings,
            stderr,
            options_override_json: None,
            effective_options_json: None,
        })
    }

    fn normalize_genome_lookup_token(raw: &str) -> String {
        raw.trim().to_ascii_lowercase()
    }

    fn infer_assembly_family_token(raw: &str) -> Option<String> {
        let token = raw.trim();
        if token.is_empty() {
            return None;
        }
        for chunk in token.split(|c: char| !c.is_ascii_alphanumeric() && c != '.') {
            let lower = chunk.to_ascii_lowercase();
            if lower.is_empty() {
                continue;
            }
            let core = lower
                .split_once(".p")
                .map(|(prefix, _)| prefix)
                .unwrap_or(lower.as_str())
                .trim();
            if core.is_empty() {
                continue;
            }
            let has_alpha = core.chars().any(|c| c.is_ascii_alphabetic());
            let has_digit = core.chars().any(|c| c.is_ascii_digit());
            if has_alpha && has_digit {
                return Some(core.to_string());
            }
        }
        None
    }

    fn entry_assembly_family(genome_key: &str, entry: &GenomeCatalogEntry) -> Option<String> {
        entry
            .ncbi_assembly_name
            .as_deref()
            .and_then(Self::infer_assembly_family_token)
            .or_else(|| Self::infer_assembly_family_token(genome_key))
            .or_else(|| {
                entry
                    .description
                    .as_deref()
                    .and_then(Self::infer_assembly_family_token)
            })
    }

    fn find_compatible_prepared_entries(
        &self,
        requested_key: &str,
        requested_entry: &GenomeCatalogEntry,
        requested_family: Option<&str>,
        cache_dir_override: Option<&str>,
    ) -> Result<Vec<String>, String> {
        let mut compatible: Vec<String> = vec![];
        for (candidate_key, candidate_entry) in &self.entries {
            if candidate_key == requested_key {
                continue;
            }
            if let (Some(req_tax), Some(candidate_tax)) = (
                requested_entry.ncbi_taxonomy_id,
                candidate_entry.ncbi_taxonomy_id,
            ) {
                if req_tax != candidate_tax {
                    continue;
                }
            }
            if let Some(family) = requested_family {
                let candidate_family =
                    Self::entry_assembly_family(candidate_key, candidate_entry).unwrap_or_default();
                if candidate_family != family {
                    continue;
                }
            }
            if self.has_sequence_ready_install(candidate_key, cache_dir_override)? {
                compatible.push(candidate_key.clone());
            }
        }
        compatible.sort_unstable();
        compatible.dedup();
        Ok(compatible)
    }

    /// Inspect prepared compatibility options for a requested genome id.
    pub fn inspect_prepared_genome_compatibility(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<PreparedGenomeCompatibilityInspection, String> {
        let requested_key = self.resolve_entry_key(genome_id)?;
        let requested_entry = self
            .entries
            .get(&requested_key)
            .ok_or_else(|| format!("Genome '{genome_id}' is not present in the catalog"))?;
        let requested_family = Self::infer_assembly_family_token(genome_id)
            .or_else(|| Self::entry_assembly_family(&requested_key, requested_entry));
        let exact_prepared = self.has_sequence_ready_install(&requested_key, cache_dir_override)?;
        let mut compatible_prepared = self.find_compatible_prepared_entries(
            &requested_key,
            requested_entry,
            requested_family.as_deref(),
            cache_dir_override,
        )?;
        if exact_prepared {
            compatible_prepared.push(requested_key.clone());
            compatible_prepared.sort_unstable();
            compatible_prepared.dedup();
        }
        Ok(PreparedGenomeCompatibilityInspection {
            requested_genome_id: genome_id.to_string(),
            requested_catalog_key: requested_key,
            requested_family,
            exact_prepared,
            compatible_prepared_options: compatible_prepared,
        })
    }

    /// Resolve a requested genome id to a prepared local install.
    ///
    /// If the exact id is not prepared, this can fall back to a unique,
    /// assembly-compatible prepared entry (for example `GRCh38.p14` ->
    /// `Human GRCh38 Ensembl 116`). If multiple compatible prepared entries
    /// exist, this returns an error that lists explicit options.
    pub fn resolve_prepared_genome_id(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<PreparedGenomeResolution, String> {
        self.resolve_prepared_genome_id_with_policy(
            genome_id,
            cache_dir_override,
            PreparedGenomeFallbackPolicy::SingleCompatible,
        )
    }

    /// Resolve prepared genome id with explicit compatibility policy.
    pub fn resolve_prepared_genome_id_with_policy(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
        fallback_policy: PreparedGenomeFallbackPolicy,
    ) -> Result<PreparedGenomeResolution, String> {
        let inspection =
            self.inspect_prepared_genome_compatibility(genome_id, cache_dir_override)?;
        if inspection.exact_prepared {
            return Ok(PreparedGenomeResolution {
                requested_genome_id: inspection.requested_genome_id,
                resolved_genome_id: inspection.requested_catalog_key,
                fallback_warning: None,
            });
        }
        let requested_key = inspection.requested_catalog_key;
        let requested_family = inspection.requested_family;
        let compatible_prepared = inspection.compatible_prepared_options;
        if compatible_prepared.is_empty() {
            return Err(format!(
                "Genome '{genome_id}' is not prepared locally. Run PrepareGenome first."
            ));
        }
        let family_label = requested_family.unwrap_or_else(|| "same-family".to_string());
        match fallback_policy {
            PreparedGenomeFallbackPolicy::Off => Err(format!(
                "Genome '{genome_id}' is not prepared locally and compatibility fallback policy is 'off'. Compatible prepared options for assembly family '{}' are: {}. Choose one explicitly or run PrepareGenome for '{}'.",
                family_label,
                compatible_prepared.join(", "),
                requested_key
            )),
            PreparedGenomeFallbackPolicy::AlwaysExplicit => Err(format!(
                "Genome '{genome_id}' is not prepared locally. Compatibility fallback policy requires explicit selection. Compatible prepared options for assembly family '{}' are: {}. Choose one explicitly or run PrepareGenome for '{}'.",
                family_label,
                compatible_prepared.join(", "),
                requested_key
            )),
            PreparedGenomeFallbackPolicy::SingleCompatible => {
                if compatible_prepared.len() > 1 {
                    return Err(format!(
                        "Genome '{genome_id}' is not prepared locally. Compatible prepared options for assembly family '{}' are: {}. Choose one explicitly or run PrepareGenome for '{}'.",
                        family_label,
                        compatible_prepared.join(", "),
                        requested_key
                    ));
                }
                let resolved_genome_id = compatible_prepared[0].clone();
                Ok(PreparedGenomeResolution {
                    requested_genome_id: genome_id.to_string(),
                    resolved_genome_id: resolved_genome_id.clone(),
                    fallback_warning: Some(format!(
                        "Genome '{}' is not prepared; using compatible prepared genome '{}' (assembly family '{}').",
                        genome_id, resolved_genome_id, family_label
                    )),
                })
            }
        }
    }

    fn has_sequence_ready_install(
        &self,
        genome_id: &str,
        cache_dir_override: Option<&str>,
    ) -> Result<bool, String> {
        let entry = self.entry(genome_id)?;
        let install_dir = self.install_dir(genome_id, entry, cache_dir_override);
        let manifest_path = install_dir.join("manifest.json");
        if !manifest_path.exists() {
            return Ok(false);
        }
        let manifest = Self::load_manifest(&manifest_path)?;
        Ok(Self::validate_manifest_files(&manifest).is_ok())
    }

    fn resolve_entry_key(&self, genome_id: &str) -> Result<String, String> {
        if self.entries.contains_key(genome_id) {
            return Ok(genome_id.to_string());
        }
        let needle = Self::normalize_genome_lookup_token(genome_id);
        if needle.is_empty() {
            return Err("Genome id must not be empty".to_string());
        }
        let mut matches: Vec<String> = vec![];
        for (key, entry) in &self.entries {
            let key_match = Self::normalize_genome_lookup_token(key) == needle;
            let alias_match = entry
                .ncbi_assembly_name
                .as_deref()
                .map(Self::normalize_genome_lookup_token)
                .map(|value| value == needle)
                .unwrap_or(false)
                || entry
                    .ncbi_assembly_accession
                    .as_deref()
                    .map(Self::normalize_genome_lookup_token)
                    .map(|value| value == needle)
                    .unwrap_or(false)
                || entry
                    .genbank_accession
                    .as_deref()
                    .map(Self::normalize_genome_lookup_token)
                    .map(|value| value == needle)
                    .unwrap_or(false);
            if key_match || alias_match {
                matches.push(key.clone());
            }
        }
        matches.sort_unstable();
        matches.dedup();
        if matches.len() == 1 {
            return Ok(matches.remove(0));
        }
        if matches.len() > 1 {
            return Err(format!(
                "Genome alias '{}' matches multiple catalog entries: {}",
                genome_id,
                matches.join(", ")
            ));
        }
        Err(format!(
            "Genome '{genome_id}' is not present in the catalog"
        ))
    }

    fn entry(&self, genome_id: &str) -> Result<&GenomeCatalogEntry, String> {
        let key = self.resolve_entry_key(genome_id)?;
        self.entries
            .get(&key)
            .ok_or_else(|| format!("Genome '{genome_id}' is not present in the catalog"))
    }

    fn install_dir(
        &self,
        genome_id: &str,
        entry: &GenomeCatalogEntry,
        cache_dir_override: Option<&str>,
    ) -> PathBuf {
        let base = cache_dir_override
            .map(|raw| self.resolve_local_path(raw))
            .or_else(|| {
                entry
                    .cache_dir
                    .as_ref()
                    .map(|raw| self.resolve_local_path(raw))
            })
            .unwrap_or_else(|| PathBuf::from(DEFAULT_GENOME_CACHE_DIR));
        base.join(sanitize_for_path(genome_id))
    }

    #[cfg(test)]
    fn resolve_source(
        &self,
        genome_id: &str,
        kind: &str,
        local: Option<&String>,
        remote: Option<&String>,
        entry: &GenomeCatalogEntry,
    ) -> Result<String, String> {
        self.resolve_source_with_type(genome_id, kind, local, remote, entry)
            .map(|r| r.source)
    }

    fn resolve_source_with_type(
        &self,
        genome_id: &str,
        kind: &str,
        local: Option<&String>,
        remote: Option<&String>,
        entry: &GenomeCatalogEntry,
    ) -> Result<SourceResolution, String> {
        let mut missing_local_path: Option<PathBuf> = None;
        if let Some(local_raw) = local {
            let local_path = self.resolve_local_path(local_raw);
            if local_path.exists() {
                return Ok(SourceResolution {
                    source: canonical_or_display(&local_path),
                    source_type: SourceType::Local,
                });
            }
            missing_local_path = Some(local_path);
        }
        if let Some(remote_src) = remote {
            return Ok(SourceResolution {
                source: remote_src.clone(),
                source_type: classify_source_type(remote_src),
            });
        }
        if let Some(ncbi_source) = resolve_ncbi_assembly_source(kind, entry)? {
            return Ok(SourceResolution {
                source: ncbi_source,
                source_type: SourceType::NcbiAssembly,
            });
        }
        if let Some(genbank_source) = resolve_genbank_accession_source(kind, entry)? {
            return Ok(SourceResolution {
                source: genbank_source,
                source_type: SourceType::GenbankAccession,
            });
        }
        if let Some(local_path) = missing_local_path {
            return Err(format!(
                "Genome '{genome_id}' has a {kind}_local path '{}', but that file does not exist and no remote/NCBI fallback is configured",
                local_path.display()
            ));
        }
        Err(format!(
            "Genome '{genome_id}' does not provide {kind}_local/{kind}_remote and has no resolvable NCBI assembly/GenBank accession source"
        ))
    }

    fn resolve_local_path(&self, raw: &str) -> PathBuf {
        if let Some(stripped) = raw.strip_prefix("file://") {
            return PathBuf::from(stripped);
        }
        let p = Path::new(raw);
        if p.is_absolute() {
            p.to_path_buf()
        } else {
            self.catalog_base_dir.join(p)
        }
    }

    fn load_manifest(path: &Path) -> Result<GenomeInstallManifest, String> {
        let text = fs::read_to_string(path)
            .map_err(|e| format!("Could not read genome manifest '{}': {e}", path.display()))?;
        serde_json::from_str(&text)
            .map_err(|e| format!("Could not parse genome manifest '{}': {e}", path.display()))
    }

    fn write_manifest(path: &Path, manifest: &GenomeInstallManifest) -> Result<(), String> {
        let text = serde_json::to_string_pretty(manifest)
            .map_err(|e| format!("Could not serialize genome manifest: {e}"))?;
        fs::write(path, text)
            .map_err(|e| format!("Could not write genome manifest '{}': {e}", path.display()))
    }

    fn validate_manifest_files(manifest: &GenomeInstallManifest) -> Result<(), String> {
        for path in [
            &manifest.sequence_path,
            &manifest.annotation_path,
            &manifest.fasta_index_path,
        ] {
            if !Path::new(path).exists() {
                return Err(format!(
                    "Genome installation is incomplete; missing file '{}'",
                    path
                ));
            }
        }
        Ok(())
    }
}

fn now_unix_ms() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

fn canonical_or_display(path: &Path) -> String {
    fs::canonicalize(path)
        .unwrap_or_else(|_| path.to_path_buf())
        .to_string_lossy()
        .into_owned()
}

pub fn is_prepare_cancelled_error(message: &str) -> bool {
    message.contains(PREPARE_CANCELLED_BY_CALLER)
}

fn prepare_cancelled_error(context: &str) -> String {
    format!("{PREPARE_CANCELLED_BY_CALLER}: {context}")
}

pub fn is_blast_cancelled_error(message: &str) -> bool {
    message.contains(BLAST_CANCELLED_BY_CALLER)
}

fn blast_cancelled_error(context: &str) -> String {
    format!("{BLAST_CANCELLED_BY_CALLER}: {context}")
}

fn forward_prepare_progress(
    on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    progress: PrepareGenomeProgress,
) -> Result<(), String> {
    if on_progress(progress) {
        Ok(())
    } else {
        Err(prepare_cancelled_error("progress callback"))
    }
}

fn summarize_annotation_parse_warnings(report: &AnnotationParseReport) -> Vec<String> {
    let mut warnings: Vec<String> = vec![];
    if report.malformed_lines == 0 {
        return warnings;
    }
    warnings.push(format!(
        "Annotation parse for '{}' skipped {} malformed line(s) while extracting genes (source={}, total_lines={})",
        report.path,
        report.malformed_lines,
        report.source_format,
        report.total_lines
    ));
    for issue in &report.issue_examples {
        warnings.push(format!(
            "Annotation parse issue at line {}: {} | {}",
            issue.line, issue.reason, issue.context
        ));
    }
    if report.truncated_issue_count > 0 {
        warnings.push(format!(
            "{} additional malformed annotation line(s) omitted from warning output",
            report.truncated_issue_count
        ));
    }
    warnings
}

fn compute_file_sha1(path: &Path) -> Result<String, String> {
    let mut file = File::open(path).map_err(|e| {
        format!(
            "Could not open '{}' to compute checksum: {e}",
            path.display()
        )
    })?;
    let mut hasher = Sha1::new();
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = file.read(&mut buf).map_err(|e| {
            format!(
                "Could not read '{}' while computing checksum: {e}",
                path.display()
            )
        })?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn ensure_manifest_file_checksum(
    slot: &mut Option<String>,
    path: &Path,
    label: &str,
) -> Result<bool, String> {
    let actual = compute_file_sha1(path)?;
    match slot {
        Some(expected) => {
            let expected_norm = expected.trim().to_ascii_lowercase();
            if expected_norm != actual {
                return Err(format!(
                    "Integrity check failed for {} file '{}': expected sha1 {}, got {}",
                    label,
                    path.display(),
                    expected,
                    actual
                ));
            }
            Ok(false)
        }
        None => {
            *slot = Some(actual);
            Ok(true)
        }
    }
}

fn ensure_manifest_checksums(manifest: &mut GenomeInstallManifest) -> Result<bool, String> {
    let mut changed = false;
    changed |= ensure_manifest_file_checksum(
        &mut manifest.sequence_sha1,
        Path::new(&manifest.sequence_path),
        "sequence",
    )?;
    changed |= ensure_manifest_file_checksum(
        &mut manifest.annotation_sha1,
        Path::new(&manifest.annotation_path),
        "annotation",
    )?;
    Ok(changed)
}

fn sanitize_for_path(s: &str) -> String {
    let mut out = String::new();
    for c in s.chars() {
        if c.is_ascii_alphanumeric() {
            out.push(c.to_ascii_lowercase());
        } else if matches!(c, ' ' | '-' | '_' | '.') {
            if !out.ends_with('_') {
                out.push('_');
            }
        }
    }
    let trimmed = out.trim_matches('_');
    if trimmed.is_empty() {
        "genome".to_string()
    } else {
        trimmed.to_string()
    }
}

fn resolve_ncbi_assembly_source(
    kind: &str,
    entry: &GenomeCatalogEntry,
) -> Result<Option<String>, String> {
    let accession = entry
        .ncbi_assembly_accession
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty());
    let assembly_name = entry
        .ncbi_assembly_name
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty());
    match (accession, assembly_name) {
        (None, None) => Ok(None),
        (Some(_), None) | (None, Some(_)) => Err(
            "NCBI assembly source requires both 'ncbi_assembly_accession' and 'ncbi_assembly_name'"
                .to_string(),
        ),
        (Some(accession), Some(assembly_name)) => {
            let normalized_accession = accession.trim().to_ascii_uppercase();
            let normalized_assembly_name = normalize_ncbi_assembly_name(assembly_name);
            let base =
                build_ncbi_assembly_ftp_base(&normalized_accession, &normalized_assembly_name)?;
            let suffix = match kind {
                "sequence" => "_genomic.fna.gz",
                "annotation" => "_genomic.gff.gz",
                _ => return Err(format!("Unsupported NCBI assembly source kind '{kind}'")),
            };
            Ok(Some(format!(
                "{base}/{}_{}{}",
                normalized_accession, normalized_assembly_name, suffix
            )))
        }
    }
}

fn resolve_genbank_accession_source(
    kind: &str,
    entry: &GenomeCatalogEntry,
) -> Result<Option<String>, String> {
    let accession = entry
        .genbank_accession
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty());
    let Some(raw_accession) = accession else {
        return Ok(None);
    };
    let accession = validate_genbank_accession(raw_accession)?;
    let rettype = match kind {
        "sequence" => "fasta",
        // Keep complete features/parts for vector-oriented annotation searches.
        "annotation" => "gbwithparts",
        _ => {
            return Err(format!(
                "Unsupported GenBank accession source kind '{kind}'"
            ));
        }
    };
    Ok(Some(build_genbank_efetch_url(&accession, rettype)))
}

pub(crate) fn validate_genbank_accession(raw: &str) -> Result<String, String> {
    let value = raw.trim();
    if value.is_empty() {
        return Err("GenBank accession is empty".to_string());
    }
    if is_unpublished_genbank_placeholder(value) {
        return Err(format!(
            "GenBank accession '{value}' is marked as unpublished/local-only and cannot be fetched remotely"
        ));
    }
    let first = value.chars().next().unwrap_or_default();
    if !first.is_ascii_alphabetic() {
        return Err(format!(
            "Invalid GenBank accession '{value}' (must start with an ASCII letter)"
        ));
    }
    if !value
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-'))
    {
        return Err(format!(
            "Invalid GenBank accession '{value}' (allowed characters: A-Z, 0-9, '_', '.', '-')"
        ));
    }
    Ok(value.to_string())
}

fn is_unpublished_genbank_placeholder(raw: &str) -> bool {
    let upper = raw.trim().to_ascii_uppercase();
    matches!(
        upper.as_str(),
        "LOCAL_UNPUBLISHED" | "NOT_UPLOADED_TO_GENBANK" | "NOT_UPLOADED"
    )
}

pub(crate) fn build_genbank_efetch_url(accession: &str, rettype: &str) -> String {
    let override_url = std::env::var(NCBI_EFETCH_ENV_VAR)
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty());
    if let Some(template) = override_url {
        if template.contains("{accession}")
            || template.contains("{rettype}")
            || template.contains("{db}")
        {
            return template
                .replace("{accession}", accession)
                .replace("{rettype}", rettype)
                .replace("{db}", "nuccore");
        }
        let sep = if template.contains('?') { "&" } else { "?" };
        return format!("{template}{sep}db=nuccore&id={accession}&rettype={rettype}&retmode=text");
    }
    format!(
        "{DEFAULT_NCBI_EFETCH_ENDPOINT}?db=nuccore&id={accession}&rettype={rettype}&retmode=text"
    )
}

fn build_ncbi_assembly_ftp_base(accession: &str, assembly_name: &str) -> Result<String, String> {
    let normalized_accession = accession.trim().to_ascii_uppercase();
    let Some((prefix, rest)) = normalized_accession.split_once('_') else {
        return Err(format!(
            "Invalid NCBI assembly accession '{accession}' (expected GCA_/GCF_ prefix)"
        ));
    };
    if prefix != "GCA" && prefix != "GCF" {
        return Err(format!(
            "Invalid NCBI assembly accession '{accession}' (expected GCA_/GCF_ prefix)"
        ));
    }
    let Some((digits, version)) = rest.split_once('.') else {
        return Err(format!(
            "Invalid NCBI assembly accession '{accession}' (expected version suffix like '.2')"
        ));
    };
    if digits.is_empty()
        || !digits.chars().all(|c| c.is_ascii_digit())
        || version.is_empty()
        || !version.chars().all(|c| c.is_ascii_digit())
    {
        return Err(format!(
            "Invalid NCBI assembly accession '{accession}' (non-numeric accession/version segment)"
        ));
    }
    let mut padded_digits = digits.to_string();
    let rem = padded_digits.len() % 3;
    if rem != 0 {
        padded_digits = format!("{}{}", "0".repeat(3 - rem), padded_digits);
    }
    let segments: Vec<String> = padded_digits
        .as_bytes()
        .chunks(3)
        .map(|chunk| String::from_utf8_lossy(chunk).to_string())
        .collect();
    let assembly_name = normalize_ncbi_assembly_name(assembly_name);
    if assembly_name.is_empty() || assembly_name.contains('/') || assembly_name.contains('\\') {
        return Err(format!(
            "Invalid NCBI assembly name '{assembly_name}' for accession '{accession}'"
        ));
    }
    Ok(format!(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}_{}",
        prefix,
        segments.join("/"),
        normalized_accession,
        assembly_name
    ))
}

fn normalize_ncbi_assembly_name(assembly_name: &str) -> String {
    assembly_name
        .trim()
        .chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .collect::<String>()
}

fn classify_source_type_label(source: &str) -> &'static str {
    source_type_label(classify_source_type(source))
}

fn source_type_label(source_type: SourceType) -> &'static str {
    match source_type {
        SourceType::Local => "local",
        SourceType::NcbiAssembly => "ncbi_assembly",
        SourceType::GenbankAccession => "genbank_accession",
        SourceType::RemoteHttp => "remote_http",
    }
}

fn classify_source_type(source: &str) -> SourceType {
    if !is_http_source(source) {
        return SourceType::Local;
    }
    let lower = source.to_ascii_lowercase();
    if lower.contains("entrez/eutils/efetch.fcgi")
        && lower.contains("db=nuccore")
        && lower.contains("rettype=")
    {
        return SourceType::GenbankAccession;
    }
    if lower.contains("ftp.ncbi.nlm.nih.gov/genomes/all/") && lower.contains("_genomic.") {
        return SourceType::NcbiAssembly;
    }
    SourceType::RemoteHttp
}

fn summarize_fasta_index(index: &HashMap<String, FastaIndexEntry>) -> PreparedSequenceCacheSummary {
    let mut contigs = index
        .iter()
        .map(|(name, entry)| (name.clone(), entry.length))
        .collect::<Vec<_>>();
    contigs.sort_by(|left, right| right.1.cmp(&left.1).then(left.0.cmp(&right.0)));
    let total_span_bp = contigs
        .iter()
        .fold(0u64, |acc, (_, length)| acc.saturating_add(*length));
    let (longest_contig, longest_contig_bp) = contigs
        .first()
        .map(|(name, length)| (Some(name.clone()), Some(*length)))
        .unwrap_or((None, None));
    PreparedSequenceCacheSummary {
        contig_count: contigs.len(),
        total_span_bp,
        longest_contig,
        longest_contig_bp,
        contig_preview: contigs
            .iter()
            .take(8)
            .map(|(name, _)| name.clone())
            .collect(),
    }
}

fn validate_annotation_gene_contigs_against_fasta_index(
    genome_id: &str,
    fasta_index: &HashMap<String, FastaIndexEntry>,
    gene_records: &[GenomeGeneRecord],
    sequence_path: &Path,
    fasta_index_path: &Path,
    gene_index_path: &Path,
) -> Result<(), String> {
    let mut referenced_contigs = gene_records
        .iter()
        .map(|gene| gene.chromosome.trim())
        .filter(|chromosome| !chromosome.is_empty())
        .map(|chromosome| chromosome.to_string())
        .collect::<Vec<_>>();
    if referenced_contigs.is_empty() {
        return Ok(());
    }
    referenced_contigs.sort_unstable();
    referenced_contigs.dedup();

    let mut missing = referenced_contigs
        .iter()
        .filter(|chromosome| {
            !fasta_index
                .keys()
                .any(|candidate| chromosome_names_match(candidate, chromosome))
        })
        .cloned()
        .collect::<Vec<_>>();
    if missing.is_empty() {
        return Ok(());
    }

    missing.sort_unstable();
    let preview_items = missing.iter().take(12).cloned().collect::<Vec<_>>();
    let hidden = missing.len().saturating_sub(preview_items.len());
    let mut preview = preview_items.join(", ");
    if hidden > 0 {
        preview.push_str(&format!(", ... (+{hidden} more)"));
    }
    Err(format!(
        "Prepared genome '{}' is inconsistent: annotation gene index '{}' references contigs \
missing from prepared sequence '{}' (missing {} of {} gene-bearing contigs: {}). FASTA \
index='{}'. This often indicates truncated gzip decode, mismatched sequence/annotation \
sources, or stale cache; reinstall or verify the configured sources.",
        genome_id,
        gene_index_path.display(),
        sequence_path.display(),
        missing.len(),
        referenced_contigs.len(),
        preview,
        fasta_index_path.display()
    ))
}

fn has_non_empty(value: &Option<String>) -> bool {
    value
        .as_ref()
        .map(|v| !v.trim().is_empty())
        .unwrap_or(false)
}

const ESTIMATED_DOUBLE_STRANDED_DNA_DA_PER_BP: f64 = 617.96;
const ESTIMATED_DOUBLE_STRANDED_DNA_DA_TERMINAL: f64 = 36.04;

fn estimate_double_stranded_dna_mass_da(length_bp: usize) -> Option<f64> {
    if length_bp == 0 {
        None
    } else {
        Some(
            (length_bp as f64) * ESTIMATED_DOUBLE_STRANDED_DNA_DA_PER_BP
                + ESTIMATED_DOUBLE_STRANDED_DNA_DA_TERMINAL,
        )
    }
}

fn validate_catalog_entries(
    catalog_path: &str,
    entries: &HashMap<String, GenomeCatalogEntry>,
) -> Result<(), String> {
    let mut errors: Vec<String> = vec![];
    for (genome_id, entry) in entries {
        let mut entry_errors: Vec<String> = vec![];
        let has_assembly_accession = has_non_empty(&entry.ncbi_assembly_accession);
        let has_assembly_name = has_non_empty(&entry.ncbi_assembly_name);
        let has_assembly = has_assembly_accession && has_assembly_name;
        if has_assembly_accession != has_assembly_name {
            entry_errors.push(
                "NCBI assembly source requires both 'ncbi_assembly_accession' and 'ncbi_assembly_name'"
                    .to_string(),
            );
        }

        let genbank_accession = entry
            .genbank_accession
            .as_ref()
            .map(|v| v.trim())
            .filter(|v| !v.is_empty());
        if has_assembly && genbank_accession.is_some() {
            entry_errors.push(
                "Entry cannot declare both NCBI assembly fields and 'genbank_accession'"
                    .to_string(),
            );
        }

        let local_unpublished = entry.local_variant_unpublished.unwrap_or(false);
        if let Some(accession) = genbank_accession {
            if is_unpublished_genbank_placeholder(accession) {
                if !local_unpublished {
                    entry_errors.push(
                        "Placeholder 'genbank_accession' requires 'local_variant_unpublished: true'"
                            .to_string(),
                    );
                }
                if has_non_empty(&entry.sequence_remote)
                    || has_non_empty(&entry.annotations_remote)
                    || has_assembly
                {
                    entry_errors.push(
                        "Placeholder 'genbank_accession' cannot be combined with remote/assembly sources"
                            .to_string(),
                    );
                }
                if !has_non_empty(&entry.sequence_local) && !has_non_empty(&entry.annotations_local)
                {
                    entry_errors.push(
                        "Placeholder 'genbank_accession' requires local sequence/annotation paths"
                            .to_string(),
                    );
                }
            } else if let Err(e) = validate_genbank_accession(accession) {
                entry_errors.push(e);
            } else if local_unpublished {
                entry_errors.push(
                    "'local_variant_unpublished: true' requires placeholder 'genbank_accession'"
                        .to_string(),
                );
            }
        } else if local_unpublished {
            entry_errors.push(
                "'local_variant_unpublished: true' requires placeholder 'genbank_accession'"
                    .to_string(),
            );
        }

        let has_sequence_source = has_non_empty(&entry.sequence_local)
            || has_non_empty(&entry.sequence_remote)
            || has_assembly
            || genbank_accession.is_some();
        let has_annotation_source = has_non_empty(&entry.annotations_local)
            || has_non_empty(&entry.annotations_remote)
            || has_assembly
            || genbank_accession.is_some();

        if !has_sequence_source {
            entry_errors.push(
                "Missing sequence source: provide sequence_local/sequence_remote or NCBI assembly/GenBank accession fields"
                    .to_string(),
            );
        }
        if !has_annotation_source {
            entry_errors.push(
                "Missing annotation source: provide annotations_local/annotations_remote or NCBI assembly/GenBank accession fields"
                    .to_string(),
            );
        }
        if let Some(length_bp) = entry.nucleotide_length_bp {
            if length_bp == 0 {
                entry_errors.push(
                    "'nucleotide_length_bp' must be a positive integer when provided".to_string(),
                );
            }
        }
        if let Some(mass_da) = entry.molecular_mass_da {
            if !mass_da.is_finite() || mass_da <= 0.0 {
                entry_errors.push(
                    "'molecular_mass_da' must be a finite positive number when provided"
                        .to_string(),
                );
            }
        }

        for err in entry_errors {
            errors.push(format!("{genome_id}: {err}"));
        }
    }
    if errors.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "Could not validate genome catalog '{catalog_path}':\n- {}",
            errors.join("\n- ")
        ))
    }
}

fn infer_annotation_extension(source: &str) -> &'static str {
    let lower = source.to_ascii_lowercase();
    if lower.contains(".gbff")
        || lower.contains(".gbk")
        || lower.contains(".genbank")
        || lower.contains(".gb")
        || lower.contains("rettype=gbwithparts")
        || lower.contains("rettype=gb")
        || lower.contains("rettype=genbank")
    {
        "gbff"
    } else if lower.contains(".gff3") {
        "gff3"
    } else if lower.contains(".gff") {
        "gff"
    } else {
        "gtf"
    }
}

fn is_http_source(source: &str) -> bool {
    source.starts_with("http://") || source.starts_with("https://")
}

fn source_compare_key(source: &str) -> String {
    if is_http_source(source) {
        return source.trim().to_string();
    }
    let path = if let Some(stripped) = source.strip_prefix("file://") {
        PathBuf::from(stripped)
    } else {
        PathBuf::from(source)
    };
    fs::canonicalize(&path)
        .unwrap_or(path)
        .to_string_lossy()
        .to_string()
}

fn sources_equivalent(left: &str, right: &str) -> bool {
    source_compare_key(left) == source_compare_key(right)
}

fn is_gzip_source(source: &str) -> bool {
    source.to_ascii_lowercase().ends_with(".gz")
}

struct SourceReader {
    reader: Box<dyn Read>,
    total_bytes: Option<u64>,
}

fn open_source_reader(source: &str) -> Result<SourceReader, String> {
    if is_http_source(source) {
        let response = fetch_http_source_with_retry(source)?;
        let total_bytes = response.content_length();
        return Ok(SourceReader {
            reader: Box::new(response),
            total_bytes,
        });
    }
    let path = if let Some(stripped) = source.strip_prefix("file://") {
        PathBuf::from(stripped)
    } else {
        PathBuf::from(source)
    };
    let file = File::open(&path)
        .map_err(|e| format!("Could not open source file '{}': {e}", path.display()))?;
    let total_bytes = file.metadata().ok().map(|m| m.len());
    Ok(SourceReader {
        reader: Box::new(file),
        total_bytes,
    })
}

fn build_http_client() -> Result<reqwest::blocking::Client, String> {
    reqwest::blocking::Client::builder()
        .connect_timeout(Duration::from_secs(HTTP_CONNECT_TIMEOUT_SECS))
        .timeout(Duration::from_secs(HTTP_READ_TIMEOUT_SECS))
        .user_agent("gentle_rs/genome_prepare")
        .build()
        .map_err(|e| format!("Could not construct HTTP client: {e}"))
}

fn is_retryable_http_status(status: reqwest::StatusCode) -> bool {
    status.is_server_error() || status.as_u16() == 429 || status.as_u16() == 408
}

fn is_retryable_http_error(err: &reqwest::Error) -> bool {
    err.is_connect() || err.is_timeout() || err.is_request() || err.is_body()
}

fn fetch_http_source_with_retry(source: &str) -> Result<reqwest::blocking::Response, String> {
    let client = build_http_client()?;
    let source_hint = if source.to_ascii_lowercase().contains("ncbi.nlm.nih.gov") {
        " (NCBI source)"
    } else {
        ""
    };
    let mut last_error: Option<String> = None;
    for attempt in 1..=HTTP_RETRY_ATTEMPTS {
        match client.get(source).send() {
            Ok(response) => {
                let status = response.status();
                if status.is_success() {
                    return Ok(response);
                }
                let retryable = is_retryable_http_status(status);
                let msg = format!(
                    "HTTP {} for '{}' (attempt {}/{})",
                    status.as_u16(),
                    source,
                    attempt,
                    HTTP_RETRY_ATTEMPTS
                );
                last_error = Some(msg.clone());
                if retryable && attempt < HTTP_RETRY_ATTEMPTS {
                    let delay_ms = HTTP_RETRY_BASE_BACKOFF_MS.saturating_mul(1 << (attempt - 1));
                    thread::sleep(Duration::from_millis(delay_ms));
                    continue;
                }
                return Err(format!(
                    "Could not fetch '{}'{}: {}{}",
                    source,
                    source_hint,
                    msg,
                    if retryable {
                        " (retries exhausted)"
                    } else {
                        ""
                    }
                ));
            }
            Err(err) => {
                let retryable = is_retryable_http_error(&err);
                let msg = format!("{} (attempt {}/{})", err, attempt, HTTP_RETRY_ATTEMPTS);
                last_error = Some(msg.clone());
                if retryable && attempt < HTTP_RETRY_ATTEMPTS {
                    let delay_ms = HTTP_RETRY_BASE_BACKOFF_MS.saturating_mul(1 << (attempt - 1));
                    thread::sleep(Duration::from_millis(delay_ms));
                    continue;
                }
                return Err(format!(
                    "Could not fetch '{}'{}: {}{}",
                    source,
                    source_hint,
                    msg,
                    if retryable {
                        " (retries exhausted)"
                    } else {
                        ""
                    }
                ));
            }
        }
    }
    Err(format!(
        "Could not fetch '{}'{}: {}",
        source,
        source_hint,
        last_error.unwrap_or_else(|| "unknown error".to_string())
    ))
}

fn parse_content_range_header(value: &str) -> Option<(u64, Option<u64>)> {
    let rest = value.trim().strip_prefix("bytes ")?;
    let (range, total) = rest.split_once('/')?;
    let (start, _end) = range.split_once('-')?;
    let start = start.trim().parse::<u64>().ok()?;
    let total = {
        let trimmed = total.trim();
        if trimmed == "*" {
            None
        } else {
            Some(trimmed.parse::<u64>().ok()?)
        }
    };
    Some((start, total))
}

fn append_path_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut out: OsString = path.as_os_str().to_os_string();
    out.push(suffix);
    PathBuf::from(out)
}

fn download_http_source_with_resume<F>(
    source: &str,
    download_path: &Path,
    on_progress: &mut F,
) -> Result<Option<u64>, String>
where
    F: FnMut(u64, Option<u64>) -> bool,
{
    if let Some(parent) = download_path.parent() {
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Could not create download directory '{}': {e}",
                parent.display()
            )
        })?;
    }
    let client = build_http_client()?;
    let source_hint = if source.to_ascii_lowercase().contains("ncbi.nlm.nih.gov") {
        " (NCBI source)"
    } else {
        ""
    };
    let mut resume_from = fs::metadata(download_path)
        .ok()
        .map(|m| m.len())
        .unwrap_or(0);
    let mut last_error: Option<String> = None;

    for attempt in 1..=HTTP_RETRY_ATTEMPTS {
        let mut request = client.get(source);
        if resume_from > 0 {
            request = request.header(reqwest::header::RANGE, format!("bytes={resume_from}-"));
        }
        let mut response = match request.send() {
            Ok(response) => response,
            Err(err) => {
                let retryable = is_retryable_http_error(&err);
                let msg = format!(
                    "request failed at offset {} (attempt {}/{}): {}",
                    resume_from, attempt, HTTP_RETRY_ATTEMPTS, err
                );
                last_error = Some(msg.clone());
                if retryable && attempt < HTTP_RETRY_ATTEMPTS {
                    let delay_ms = HTTP_RETRY_BASE_BACKOFF_MS.saturating_mul(1 << (attempt - 1));
                    thread::sleep(Duration::from_millis(delay_ms));
                    continue;
                }
                return Err(format!(
                    "Could not download '{}'{}: {}{}",
                    source,
                    source_hint,
                    msg,
                    if retryable {
                        " (retries exhausted)"
                    } else {
                        ""
                    }
                ));
            }
        };

        let status = response.status();
        if !(status.is_success() || status == reqwest::StatusCode::PARTIAL_CONTENT) {
            let retryable = is_retryable_http_status(status);
            let msg = format!(
                "HTTP {} at offset {} (attempt {}/{})",
                status.as_u16(),
                resume_from,
                attempt,
                HTTP_RETRY_ATTEMPTS
            );
            last_error = Some(msg.clone());
            if retryable && attempt < HTTP_RETRY_ATTEMPTS {
                let delay_ms = HTTP_RETRY_BASE_BACKOFF_MS.saturating_mul(1 << (attempt - 1));
                thread::sleep(Duration::from_millis(delay_ms));
                continue;
            }
            return Err(format!(
                "Could not download '{}'{}: {}{}",
                source,
                source_hint,
                msg,
                if retryable {
                    " (retries exhausted)"
                } else {
                    ""
                }
            ));
        }

        let mut append = resume_from > 0 && status == reqwest::StatusCode::PARTIAL_CONTENT;
        let mut base_done = if append { resume_from } else { 0 };
        let overall_total: Option<u64> = if append {
            let parsed_range = response
                .headers()
                .get(reqwest::header::CONTENT_RANGE)
                .and_then(|v| v.to_str().ok())
                .and_then(parse_content_range_header);
            if let Some((range_start, range_total)) = parsed_range {
                if range_start != resume_from {
                    // Server resumed from unexpected offset. Fall back to restart.
                    append = false;
                    base_done = 0;
                    response.content_length()
                } else {
                    range_total.or_else(|| {
                        response
                            .content_length()
                            .map(|v| range_start.saturating_add(v))
                    })
                }
            } else {
                response
                    .content_length()
                    .map(|v| base_done.saturating_add(v))
            }
        } else {
            response.content_length()
        };

        let mut writer = if append {
            OpenOptions::new()
                .create(true)
                .append(true)
                .open(download_path)
                .map_err(|e| {
                    format!(
                        "Could not open partial download '{}' for append: {e}",
                        download_path.display()
                    )
                })?
        } else {
            OpenOptions::new()
                .create(true)
                .truncate(true)
                .write(true)
                .open(download_path)
                .map_err(|e| {
                    format!(
                        "Could not create download file '{}': {e}",
                        download_path.display()
                    )
                })?
        };

        if !on_progress(base_done, overall_total) {
            return Err(prepare_cancelled_error("HTTP download progress"));
        }
        let mut transferred: u64 = 0;
        let mut copy_error: Option<String> = None;
        let mut buf = [0u8; 64 * 1024];
        loop {
            match response.read(&mut buf) {
                Ok(0) => break,
                Ok(n) => {
                    if let Err(e) = writer.write_all(&buf[..n]) {
                        copy_error = Some(format!(
                            "Could not write partial download '{}': {e}",
                            download_path.display()
                        ));
                        break;
                    }
                    transferred = transferred.saturating_add(n as u64);
                    if !on_progress(base_done.saturating_add(transferred), overall_total) {
                        return Err(prepare_cancelled_error("HTTP download progress"));
                    }
                }
                Err(e) => {
                    copy_error = Some(format!("Could not read HTTP response body: {e}"));
                    break;
                }
            }
        }
        if let Err(e) = writer.flush() {
            copy_error = Some(format!(
                "Could not flush partial download '{}': {e}",
                download_path.display()
            ));
        }
        let final_done = fs::metadata(download_path)
            .ok()
            .map(|m| m.len())
            .unwrap_or_else(|| base_done.saturating_add(transferred));

        if let Some(err) = copy_error {
            let msg = format!(
                "{} (offset {} attempt {}/{})",
                err, final_done, attempt, HTTP_RETRY_ATTEMPTS
            );
            last_error = Some(msg.clone());
            if attempt < HTTP_RETRY_ATTEMPTS {
                resume_from = final_done;
                let delay_ms = HTTP_RETRY_BASE_BACKOFF_MS.saturating_mul(1 << (attempt - 1));
                thread::sleep(Duration::from_millis(delay_ms));
                continue;
            }
            return Err(format!(
                "Could not download '{}'{}: {} (retries exhausted)",
                source, source_hint, msg
            ));
        }

        if let Some(total) = overall_total {
            if final_done < total {
                let msg = format!(
                    "incomplete download for '{}' (got {} of {} bytes at attempt {}/{})",
                    source, final_done, total, attempt, HTTP_RETRY_ATTEMPTS
                );
                last_error = Some(msg.clone());
                if attempt < HTTP_RETRY_ATTEMPTS {
                    resume_from = final_done;
                    let delay_ms = HTTP_RETRY_BASE_BACKOFF_MS.saturating_mul(1 << (attempt - 1));
                    thread::sleep(Duration::from_millis(delay_ms));
                    continue;
                }
                return Err(format!(
                    "Could not download '{}'{}: {} (retries exhausted)",
                    source, source_hint, msg
                ));
            }
            if final_done > total {
                let msg = format!(
                    "download size mismatch for '{}' (got {} bytes, expected {})",
                    source, final_done, total
                );
                last_error = Some(msg.clone());
                if attempt < HTTP_RETRY_ATTEMPTS {
                    resume_from = 0;
                    let delay_ms = HTTP_RETRY_BASE_BACKOFF_MS.saturating_mul(1 << (attempt - 1));
                    thread::sleep(Duration::from_millis(delay_ms));
                    continue;
                }
                return Err(format!(
                    "Could not download '{}'{}: {}",
                    source, source_hint, msg
                ));
            }
            if !on_progress(total, Some(total)) {
                return Err(prepare_cancelled_error("HTTP download completion"));
            }
            return Ok(Some(total));
        }

        if !on_progress(final_done, Some(final_done)) {
            return Err(prepare_cancelled_error("HTTP download completion"));
        }
        return Ok(Some(final_done));
    }

    Err(format!(
        "Could not download '{}'{}: {}",
        source,
        source_hint,
        last_error.unwrap_or_else(|| "unknown error".to_string())
    ))
}

struct ProgressReader<R, F> {
    inner: R,
    callback: F,
    bytes_done: u64,
}

impl<R, F> ProgressReader<R, F> {
    fn new(inner: R, callback: F) -> Self {
        Self {
            inner,
            callback,
            bytes_done: 0,
        }
    }
}

impl<R: Read, F: FnMut(u64) -> bool> Read for ProgressReader<R, F> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.bytes_done += n as u64;
        if !(self.callback)(self.bytes_done) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Interrupted,
                PREPARE_CANCELLED_BY_CALLER,
            ));
        }
        Ok(n)
    }
}

fn materialize_source_with_progress<F>(
    source: &str,
    destination: &Path,
    mut on_progress: F,
) -> Result<(), String>
where
    F: FnMut(u64, Option<u64>) -> bool,
{
    if let Some(parent) = destination.parent() {
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Could not create destination directory '{}': {e}",
                parent.display()
            )
        })?;
    }
    let tmp_path = append_path_suffix(destination, ".part");

    if is_http_source(source) {
        if is_gzip_source(source) {
            let compressed_path = append_path_suffix(destination, ".download.part");
            let compressed_total =
                download_http_source_with_resume(source, &compressed_path, &mut on_progress)?;
            let compressed_size = compressed_total
                .or_else(|| fs::metadata(&compressed_path).ok().map(|m| m.len()))
                .unwrap_or(0);
            if !on_progress(compressed_size, compressed_total.or(Some(compressed_size))) {
                return Err(prepare_cancelled_error("annotation/sequence gzip download"));
            }
            let mut writer = BufWriter::new(
                File::create(&tmp_path)
                    .map_err(|e| format!("Could not create '{}': {e}", tmp_path.display()))?,
            );
            let mut decoder = MultiGzDecoder::new(BufReader::new(
                File::open(&compressed_path).map_err(|e| {
                    format!(
                        "Could not open downloaded gzip '{}' for decode: {e}",
                        compressed_path.display()
                    )
                })?,
            ));
            if let Err(e) = std::io::copy(&mut decoder, &mut writer)
                .map_err(|e| format!("Could not decompress '{source}': {e}"))
            {
                let _ = fs::remove_file(&tmp_path);
                return Err(e);
            }
            writer
                .flush()
                .map_err(|e| format!("Could not flush '{}': {e}", tmp_path.display()))?;
            fs::rename(&tmp_path, destination).map_err(|e| {
                format!(
                    "Could not finalize destination '{}': {e}",
                    destination.display()
                )
            })?;
            let _ = fs::remove_file(&compressed_path);
            return Ok(());
        }

        let downloaded_total =
            download_http_source_with_resume(source, &tmp_path, &mut on_progress)?;
        let done = downloaded_total
            .or_else(|| fs::metadata(&tmp_path).ok().map(|m| m.len()))
            .unwrap_or(0);
        if !on_progress(done, downloaded_total.or(Some(done))) {
            return Err(prepare_cancelled_error("annotation/sequence download"));
        }
        fs::rename(&tmp_path, destination).map_err(|e| {
            format!(
                "Could not finalize destination '{}': {e}",
                destination.display()
            )
        })?;
        return Ok(());
    }

    let SourceReader {
        reader,
        total_bytes,
    } = open_source_reader(source)?;
    if !on_progress(0, total_bytes) {
        return Err(prepare_cancelled_error("local source copy"));
    }
    let mut writer = BufWriter::new(
        File::create(&tmp_path)
            .map_err(|e| format!("Could not create '{}': {e}", tmp_path.display()))?,
    );

    let copy_result = if is_gzip_source(source) {
        let progress_reader = ProgressReader::new(reader, |done| on_progress(done, total_bytes));
        let mut decoder = MultiGzDecoder::new(progress_reader);
        std::io::copy(&mut decoder, &mut writer)
            .map_err(|e| format!("Could not decompress '{source}': {e}"))
    } else {
        let mut progress_reader =
            ProgressReader::new(reader, |done| on_progress(done, total_bytes));
        std::io::copy(&mut progress_reader, &mut writer)
            .map_err(|e| format!("Could not copy '{source}': {e}"))
    };

    if let Err(e) = copy_result {
        let _ = fs::remove_file(&tmp_path);
        return Err(e);
    }

    writer
        .flush()
        .map_err(|e| format!("Could not flush '{}': {e}", tmp_path.display()))?;
    let final_done = fs::metadata(&tmp_path).ok().map(|m| m.len()).unwrap_or(0);
    if !on_progress(total_bytes.unwrap_or(final_done), total_bytes) {
        return Err(prepare_cancelled_error("local source copy completion"));
    }
    fs::rename(&tmp_path, destination).map_err(|e| {
        format!(
            "Could not finalize destination '{}': {e}",
            destination.display()
        )
    })?;
    Ok(())
}

fn non_empty_regular_file_exists(path: &Path) -> bool {
    fs::metadata(path)
        .map(|meta| meta.is_file() && meta.len() > 0)
        .unwrap_or(false)
}

fn resolve_tool_executable(env_var: &str, default_bin: &str) -> String {
    crate::tool_overrides::resolve_tool_executable(env_var, default_bin)
}

fn resolve_executable_path(executable: &str) -> Option<String> {
    let trimmed = executable.trim();
    if trimmed.is_empty() {
        return None;
    }
    let candidate = Path::new(trimmed);
    if candidate.is_absolute() || candidate.components().count() > 1 {
        if candidate.exists() {
            return Some(canonical_or_display(candidate));
        }
        return None;
    }
    let path_env = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path_env) {
        let joined = dir.join(trimmed);
        if joined.exists() {
            return Some(canonical_or_display(&joined));
        }
    }
    None
}

fn probe_external_binary(
    env_var: &str,
    tool_name: &str,
    executable: String,
    version_args: &[&str],
) -> ExternalBinaryPreflightProbe {
    let mut report = ExternalBinaryPreflightProbe {
        tool: tool_name.to_string(),
        env_var: env_var.to_string(),
        executable: executable.clone(),
        resolved_path: resolve_executable_path(&executable),
        ..ExternalBinaryPreflightProbe::default()
    };
    match Command::new(&executable).args(version_args).output() {
        Ok(output) => {
            report.found = true;
            report.status_code = output.status.code();
            report.version_probe_ok = output.status.success();
            let stdout = String::from_utf8_lossy(&output.stdout).to_string();
            let stderr = String::from_utf8_lossy(&output.stderr).to_string();
            let first_line = first_non_empty_output_line(&stdout, &stderr);
            if first_line != "no output" {
                report.version = Some(first_line.clone());
            }
            if !report.version_probe_ok {
                report.detail = Some(first_line);
            }
        }
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {
            report.found = false;
            report.error = Some(format!(
                "Executable '{}' was not found in PATH or configured path",
                executable
            ));
        }
        Err(err) => {
            report.found = report.resolved_path.is_some();
            report.error = Some(err.to_string());
        }
    }
    report
}

fn default_blast_db_prefix(install_dir: &Path) -> PathBuf {
    install_dir.join("blastdb").join("genome")
}

fn reset_genome_install_dir(install_dir: &Path) -> Result<(), String> {
    if !install_dir.exists() {
        return Ok(());
    }
    fs::remove_dir_all(install_dir).map_err(|e| {
        format!(
            "Could not remove stale genome install dir '{}': {e}",
            install_dir.display()
        )
    })
}

fn remove_optional_file(path: &Path, label: &str) -> Result<(), String> {
    if !path.exists() {
        return Ok(());
    }
    fs::remove_file(path)
        .map_err(|e| format!("Could not remove stale {label} '{}': {e}", path.display()))
}

fn reset_prepared_genome_index_artifacts(
    fasta_index_path: &Path,
    gene_index_path: &Path,
    blast_prefix_path: &Path,
) -> Result<(), String> {
    remove_optional_file(fasta_index_path, "FASTA index")?;
    remove_optional_file(gene_index_path, "gene index")?;
    for blast_index_path in collect_blast_index_files(blast_prefix_path) {
        fs::remove_file(&blast_index_path).map_err(|e| {
            format!(
                "Could not remove stale BLAST index '{}': {e}",
                blast_index_path.display()
            )
        })?;
    }
    Ok(())
}

fn is_blast_index_suffix(suffix: &str) -> bool {
    matches!(
        suffix,
        "nhr" | "nin" | "nsq" | "ndb" | "not" | "ntf" | "nto" | "nog" | "nos" | "nsd" | "nsi"
    )
}

fn collect_blast_index_files(db_prefix: &Path) -> Vec<PathBuf> {
    let Some(parent) = db_prefix.parent() else {
        return vec![];
    };
    let Some(base_name) = db_prefix.file_name().and_then(|v| v.to_str()) else {
        return vec![];
    };
    let prefix = format!("{base_name}.");
    let mut files: Vec<PathBuf> = fs::read_dir(parent)
        .ok()
        .into_iter()
        .flat_map(|iter| iter.flatten())
        .filter_map(|entry| {
            let name = entry.file_name();
            let name = name.to_str()?;
            if !name.starts_with(&prefix) {
                return None;
            }
            let suffix = name.strip_prefix(&prefix)?;
            if !is_blast_index_suffix(suffix) {
                return None;
            }
            Some(entry.path())
        })
        .collect();
    files.sort();
    files
}

fn is_blast_index_ready(index_files: &[PathBuf]) -> bool {
    let mut has_nhr = false;
    let mut has_nin = false;
    let mut has_nsq = false;
    let mut has_ndb = false;
    for path in index_files {
        let suffix = path
            .extension()
            .and_then(|v| v.to_str())
            .unwrap_or_default();
        match suffix {
            "nhr" => has_nhr = true,
            "nin" => has_nin = true,
            "nsq" => has_nsq = true,
            "ndb" => has_ndb = true,
            _ => {}
        }
    }
    (has_nhr && has_nin && has_nsq) || has_ndb
}

fn first_non_empty_output_line(stdout: &str, stderr: &str) -> String {
    stdout
        .lines()
        .chain(stderr.lines())
        .find(|line| !line.trim().is_empty())
        .map(|line| line.trim().to_string())
        .unwrap_or_else(|| "no output".to_string())
}

fn wait_for_child_output_with_cancel(
    child: &mut Child,
    executable: &str,
    args: &[String],
    should_cancel: &mut dyn FnMut() -> bool,
) -> Result<Output, String> {
    let Some(stdout_pipe) = child.stdout.take() else {
        return Err(format!(
            "Could not capture blastn stdout for '{}' with args [{}]",
            executable,
            args.join(" ")
        ));
    };
    let Some(stderr_pipe) = child.stderr.take() else {
        return Err(format!(
            "Could not capture blastn stderr for '{}' with args [{}]",
            executable,
            args.join(" ")
        ));
    };

    let stdout_reader = thread::spawn(move || {
        let mut buf = Vec::new();
        let mut reader = BufReader::new(stdout_pipe);
        let _ = reader.read_to_end(&mut buf);
        buf
    });
    let stderr_reader = thread::spawn(move || {
        let mut buf = Vec::new();
        let mut reader = BufReader::new(stderr_pipe);
        let _ = reader.read_to_end(&mut buf);
        buf
    });

    loop {
        if should_cancel() {
            let _ = child.kill();
            let _ = child.wait();
            let _ = stdout_reader.join();
            let _ = stderr_reader.join();
            return Err(blast_cancelled_error("running blastn process"));
        }

        match child.try_wait() {
            Ok(Some(status)) => {
                let stdout = stdout_reader.join().unwrap_or_else(|_| vec![]);
                let stderr = stderr_reader.join().unwrap_or_else(|_| vec![]);
                return Ok(Output {
                    status,
                    stdout,
                    stderr,
                });
            }
            Ok(None) => {
                thread::sleep(Duration::from_millis(50));
            }
            Err(e) => {
                let _ = child.kill();
                let _ = child.wait();
                let _ = stdout_reader.join();
                let _ = stderr_reader.join();
                return Err(format!(
                    "Could not poll blastn process '{}' with args [{}]: {}",
                    executable,
                    args.join(" "),
                    e
                ));
            }
        }
    }
}

fn ensure_blast_index(sequence_path: &Path, db_prefix: &Path) -> BlastIndexOutcome {
    let db_prefix_string = canonical_or_display(db_prefix);
    let mut outcome = BlastIndexOutcome {
        ready: false,
        executable: None,
        warnings: vec![],
    };

    let existing = collect_blast_index_files(db_prefix);
    if is_blast_index_ready(&existing) {
        outcome.ready = true;
        return outcome;
    }

    let executable = resolve_tool_executable(MAKEBLASTDB_ENV_BIN, DEFAULT_MAKEBLASTDB_BIN);
    outcome.executable = Some(executable.clone());

    if let Some(parent) = db_prefix.parent() {
        if let Err(e) = fs::create_dir_all(parent) {
            outcome.warnings.push(format!(
                "BLAST indexing skipped for '{}': could not create directory '{}': {}",
                db_prefix_string,
                parent.display(),
                e
            ));
            return outcome;
        }
    }

    let args = vec![
        "-in".to_string(),
        canonical_or_display(sequence_path),
        "-dbtype".to_string(),
        "nucl".to_string(),
        "-out".to_string(),
        db_prefix_string.clone(),
        "-parse_seqids".to_string(),
    ];
    match Command::new(&executable).args(&args).output() {
        Ok(output) => {
            let stdout = String::from_utf8_lossy(&output.stdout).to_string();
            let stderr = String::from_utf8_lossy(&output.stderr).to_string();
            if !output.status.success() {
                outcome.warnings.push(format!(
                    "BLAST indexing skipped for '{}': makeblastdb failed (status {:?}): {}",
                    db_prefix_string,
                    output.status.code(),
                    first_non_empty_output_line(&stdout, &stderr)
                ));
                return outcome;
            }
            let produced = collect_blast_index_files(db_prefix);
            if is_blast_index_ready(&produced) {
                outcome.ready = true;
            } else {
                outcome.warnings.push(format!(
                    "BLAST indexing skipped for '{}': makeblastdb reported success but index files were not produced",
                    db_prefix_string
                ));
            }
        }
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            outcome.warnings.push(format!(
                "BLAST indexing skipped for '{}': makeblastdb executable '{}' not found (set {} to override)",
                db_prefix_string, executable, MAKEBLASTDB_ENV_BIN
            ));
        }
        Err(e) => {
            outcome.warnings.push(format!(
                "BLAST indexing skipped for '{}': could not run makeblastdb executable '{}': {}",
                db_prefix_string, executable, e
            ));
        }
    }
    outcome
}

fn normalize_blast_task(task: Option<&str>) -> Result<String, String> {
    let normalized = task
        .map(str::trim)
        .filter(|v| !v.is_empty())
        .unwrap_or("blastn-short")
        .to_ascii_lowercase();
    match normalized.as_str() {
        "blastn-short" | "blastn" => Ok(normalized),
        other => Err(format!(
            "Unsupported BLAST task '{}'. Expected one of: blastn-short, blastn",
            other
        )),
    }
}

fn normalize_blast_query_sequence(raw: &str) -> Result<String, String> {
    let normalized: String = raw
        .chars()
        .filter(|c| !c.is_ascii_whitespace())
        .map(|c| match c.to_ascii_uppercase() {
            'U' => 'T',
            other => other,
        })
        .collect();
    if normalized.is_empty() {
        return Err("BLAST query sequence is empty".to_string());
    }
    if normalized
        .chars()
        .any(|c| !c.is_ascii_alphabetic() && c != '*')
    {
        return Err(
            "BLAST query contains invalid characters; only nucleotide/IUPAC letters are supported"
                .to_string(),
        );
    }
    Ok(normalized)
}

fn parse_blastn_tabular_hits(stdout: &str) -> (Vec<BlastHit>, Vec<String>) {
    let mut hits: Vec<BlastHit> = vec![];
    let mut warnings: Vec<String> = vec![];
    for (idx, line) in stdout.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 13 {
            if warnings.len() < 20 {
                warnings.push(format!(
                    "BLAST output line {} skipped: expected 13 tab-separated fields, got {}",
                    idx + 1,
                    cols.len()
                ));
            }
            continue;
        }
        let parse_usize = |raw: &str, name: &str| -> Result<usize, String> {
            raw.parse::<usize>()
                .map_err(|e| format!("could not parse {name}='{raw}': {e}"))
        };
        let parse_f64 = |raw: &str, name: &str| -> Result<f64, String> {
            raw.parse::<f64>()
                .map_err(|e| format!("could not parse {name}='{raw}': {e}"))
        };
        let parsed: Result<BlastHit, String> = (|| {
            Ok(BlastHit {
                subject_id: cols[1].to_string(),
                identity_percent: parse_f64(cols[2], "pident")?,
                alignment_length: parse_usize(cols[3], "length")?,
                mismatches: parse_usize(cols[4], "mismatch")?,
                gap_opens: parse_usize(cols[5], "gapopen")?,
                query_start: parse_usize(cols[6], "qstart")?,
                query_end: parse_usize(cols[7], "qend")?,
                subject_start: parse_usize(cols[8], "sstart")?,
                subject_end: parse_usize(cols[9], "send")?,
                evalue: parse_f64(cols[10], "evalue")?,
                bit_score: parse_f64(cols[11], "bitscore")?,
                query_coverage_percent: cols[12].parse::<f64>().ok(),
            })
        })();
        match parsed {
            Ok(hit) => hits.push(hit),
            Err(e) => {
                if warnings.len() < 20 {
                    warnings.push(format!("BLAST output line {} skipped: {}", idx + 1, e));
                }
            }
        }
    }
    (hits, warnings)
}

fn build_fasta_index_with_progress<F>(
    fasta_path: &Path,
    index_path: &Path,
    mut on_progress: F,
) -> Result<(), String>
where
    F: FnMut(u64, Option<u64>) -> bool,
{
    let total_bytes = fs::metadata(fasta_path).ok().map(|m| m.len());
    let file = File::open(fasta_path)
        .map_err(|e| format!("Could not open FASTA '{}': {e}", fasta_path.display()))?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut byte_offset: u64 = 0;
    let mut entries: Vec<(String, FastaIndexEntry, bool)> = Vec::new();
    // tuple fields: (name, entry, saw_short_line)
    let mut active: Option<(String, FastaIndexEntry, bool)> = None;
    if !on_progress(0, total_bytes) {
        return Err(prepare_cancelled_error("FASTA indexing start"));
    }

    loop {
        line.clear();
        let bytes_read = reader
            .read_line(&mut line)
            .map_err(|e| format!("Could not read FASTA '{}': {e}", fasta_path.display()))?;
        if bytes_read == 0 {
            break;
        }
        let stripped = line.trim_end_matches(&['\n', '\r'][..]);

        if stripped.starts_with('>') {
            if let Some((name, entry, _)) = active.take() {
                if entry.length == 0 {
                    return Err(format!(
                        "FASTA '{}' has empty sequence record '{}'",
                        fasta_path.display(),
                        name
                    ));
                }
                entries.push((name, entry, false));
            }
            let name = stripped
                .trim_start_matches('>')
                .split_whitespace()
                .next()
                .ok_or_else(|| {
                    format!("FASTA '{}' has malformed header line", fasta_path.display())
                })?
                .to_string();
            active = Some((
                name,
                FastaIndexEntry {
                    length: 0,
                    offset: 0,
                    line_bases: 0,
                    line_bytes: 0,
                },
                false,
            ));
        } else if !stripped.is_empty() {
            let seq_bases = stripped.len() as u64;
            if let Some((_, ref mut entry, ref mut saw_short_line)) = active {
                if entry.line_bases == 0 {
                    entry.offset = byte_offset;
                    entry.line_bases = seq_bases;
                    entry.line_bytes = bytes_read as u64;
                } else {
                    if seq_bases > entry.line_bases {
                        return Err(format!(
                            "FASTA '{}' has inconsistent line length; line longer than first line",
                            fasta_path.display()
                        ));
                    }
                    if *saw_short_line && seq_bases == entry.line_bases {
                        return Err(format!(
                            "FASTA '{}' has inconsistent line length after a short line",
                            fasta_path.display()
                        ));
                    }
                    if seq_bases < entry.line_bases {
                        *saw_short_line = true;
                    }
                }
                entry.length += seq_bases;
            } else {
                return Err(format!(
                    "FASTA '{}' contains sequence data before first header",
                    fasta_path.display()
                ));
            }
        }
        byte_offset += bytes_read as u64;
        if !on_progress(byte_offset, total_bytes) {
            return Err(prepare_cancelled_error("FASTA indexing progress"));
        }
    }

    if let Some((name, entry, _)) = active {
        if entry.length == 0 {
            return Err(format!(
                "FASTA '{}' has empty sequence record '{}'",
                fasta_path.display(),
                name
            ));
        }
        entries.push((name, entry, false));
    }
    if entries.is_empty() {
        return Err(format!(
            "FASTA '{}' does not contain any sequence records",
            fasta_path.display()
        ));
    }

    if let Some(parent) = index_path.parent() {
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Could not create FASTA index directory '{}': {e}",
                parent.display()
            )
        })?;
    }
    let mut writer = BufWriter::new(File::create(index_path).map_err(|e| {
        format!(
            "Could not create FASTA index '{}': {e}",
            index_path.display()
        )
    })?);
    for (name, entry, _) in entries {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            name, entry.length, entry.offset, entry.line_bases, entry.line_bytes
        )
        .map_err(|e| {
            format!(
                "Could not write FASTA index '{}': {e}",
                index_path.display()
            )
        })?;
    }
    writer.flush().map_err(|e| {
        format!(
            "Could not flush FASTA index '{}': {e}",
            index_path.display()
        )
    })?;
    if !on_progress(total_bytes.unwrap_or(byte_offset), total_bytes) {
        return Err(prepare_cancelled_error("FASTA indexing completion"));
    }
    Ok(())
}

fn load_fasta_index(path: &Path) -> Result<HashMap<String, FastaIndexEntry>, String> {
    let file = File::open(path)
        .map_err(|e| format!("Could not open FASTA index '{}': {e}", path.display()))?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for (i, line) in reader.lines().enumerate() {
        let line =
            line.map_err(|e| format!("Could not read FASTA index '{}': {e}", path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 5 {
            return Err(format!(
                "Invalid FASTA index line {} in '{}': expected 5 tab-separated fields",
                i + 1,
                path.display()
            ));
        }
        let name = cols[0].to_string();
        let length = cols[1].parse::<u64>().map_err(|e| {
            format!(
                "Invalid FASTA index length '{}' at line {} in '{}': {e}",
                cols[1],
                i + 1,
                path.display()
            )
        })?;
        let offset = cols[2].parse::<u64>().map_err(|e| {
            format!(
                "Invalid FASTA index offset '{}' at line {} in '{}': {e}",
                cols[2],
                i + 1,
                path.display()
            )
        })?;
        let line_bases = cols[3].parse::<u64>().map_err(|e| {
            format!(
                "Invalid FASTA index line_bases '{}' at line {} in '{}': {e}",
                cols[3],
                i + 1,
                path.display()
            )
        })?;
        let line_bytes = cols[4].parse::<u64>().map_err(|e| {
            format!(
                "Invalid FASTA index line_bytes '{}' at line {} in '{}': {e}",
                cols[4],
                i + 1,
                path.display()
            )
        })?;
        if line_bases == 0 || line_bytes == 0 {
            return Err(format!(
                "Invalid FASTA index line {} in '{}': line_bases/line_bytes must be > 0",
                i + 1,
                path.display()
            ));
        }
        map.insert(
            name,
            FastaIndexEntry {
                length,
                offset,
                line_bases,
                line_bytes,
            },
        );
    }
    if map.is_empty() {
        return Err(format!("FASTA index '{}' is empty", path.display()));
    }
    Ok(map)
}

fn read_fasta_slice(
    fasta_path: &Path,
    entry: &FastaIndexEntry,
    start_1based: u64,
    end_1based: u64,
) -> Result<String, String> {
    let start0 = start_1based - 1;
    let target_len = (end_1based - start_1based + 1) as usize;
    let row = start0 / entry.line_bases;
    let col = start0 % entry.line_bases;
    let seek_pos = entry
        .offset
        .checked_add(row.saturating_mul(entry.line_bytes))
        .and_then(|v| v.checked_add(col))
        .ok_or_else(|| "Computed FASTA seek position overflowed".to_string())?;

    let mut file = File::open(fasta_path)
        .map_err(|e| format!("Could not open FASTA '{}': {e}", fasta_path.display()))?;
    file.seek(SeekFrom::Start(seek_pos))
        .map_err(|e| format!("Could not seek FASTA '{}': {e}", fasta_path.display()))?;
    let mut reader = BufReader::new(file);

    let mut out: Vec<u8> = Vec::with_capacity(target_len);
    let mut chunk = [0u8; 8192];
    while out.len() < target_len {
        let n = reader
            .read(&mut chunk)
            .map_err(|e| format!("Could not read FASTA '{}': {e}", fasta_path.display()))?;
        if n == 0 {
            break;
        }
        for b in &chunk[..n] {
            if *b == b'\n' || *b == b'\r' {
                continue;
            }
            out.push(b.to_ascii_uppercase());
            if out.len() >= target_len {
                break;
            }
        }
    }
    if out.len() != target_len {
        return Err(format!(
            "Could not read requested interval from FASTA '{}'; expected {} bases, got {}",
            fasta_path.display(),
            target_len,
            out.len()
        ));
    }
    String::from_utf8(out).map_err(|e| format!("Extracted sequence is not valid UTF-8: {e}"))
}

fn build_gene_index_file<F>(
    annotation_path: &Path,
    gene_index_path: &Path,
    on_progress: F,
) -> Result<AnnotationParseReport, String>
where
    F: FnMut(u64, Option<u64>) -> bool,
{
    let (genes, report) =
        parse_annotation_gene_records_with_progress(annotation_path, on_progress)?;
    if let Some(parent) = gene_index_path.parent() {
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Could not create gene-index directory '{}': {e}",
                parent.display()
            )
        })?;
    }
    let file = File::create(gene_index_path).map_err(|e| {
        format!(
            "Could not create gene-index file '{}': {e}",
            gene_index_path.display()
        )
    })?;
    serde_json::to_writer(BufWriter::new(file), &genes).map_err(|e| {
        format!(
            "Could not write gene index '{}': {e}",
            gene_index_path.display()
        )
    })?;
    Ok(report)
}

fn load_gene_index_file(gene_index_path: &Path) -> Result<Vec<GenomeGeneRecord>, String> {
    let file = File::open(gene_index_path).map_err(|e| {
        format!(
            "Could not open gene index '{}': {e}",
            gene_index_path.display()
        )
    })?;
    serde_json::from_reader(BufReader::new(file)).map_err(|e| {
        format!(
            "Could not parse gene index '{}': {e}",
            gene_index_path.display()
        )
    })
}

fn parse_annotation_gene_records_with_progress<F>(
    path: &Path,
    on_progress: F,
) -> Result<(Vec<GenomeGeneRecord>, AnnotationParseReport), String>
where
    F: FnMut(u64, Option<u64>) -> bool,
{
    if is_genbank_annotation_path(path) {
        parse_genbank_gene_records_with_progress(path, on_progress)
    } else if is_xml_annotation_path(path) {
        parse_xml_gbseq_gene_records_with_progress(path, on_progress)
    } else {
        parse_tabular_annotation_gene_records_with_progress(path, on_progress)
    }
}

fn is_genbank_annotation_path(path: &Path) -> bool {
    let lower = path
        .file_name()
        .and_then(|v| v.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    lower.ends_with(".gb")
        || lower.ends_with(".gbk")
        || lower.ends_with(".gbff")
        || lower.ends_with(".genbank")
}

fn is_xml_annotation_path(path: &Path) -> bool {
    path.file_name()
        .and_then(|v| v.to_str())
        .map(|name| name.to_ascii_lowercase().ends_with(".xml"))
        .unwrap_or(false)
}

fn parse_tabular_annotation_gene_records_with_progress<F>(
    path: &Path,
    mut on_progress: F,
) -> Result<(Vec<GenomeGeneRecord>, AnnotationParseReport), String>
where
    F: FnMut(u64, Option<u64>) -> bool,
{
    let total_bytes = fs::metadata(path).ok().map(|m| m.len());
    let file = File::open(path)
        .map_err(|e| format!("Could not open annotation file '{}': {e}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut genes: Vec<GenomeGeneRecord> = vec![];
    let mut total_lines: usize = 0;
    let mut skipped_lines: usize = 0;
    let mut malformed_lines: usize = 0;
    let mut issue_examples: Vec<AnnotationParseIssue> = vec![];
    let mut bytes_read_total: u64 = 0;
    let mut raw_line: Vec<u8> = vec![];
    if !on_progress(0, total_bytes) {
        return Err(prepare_cancelled_error("annotation parse start"));
    }
    loop {
        raw_line.clear();
        let line_bytes = reader
            .read_until(b'\n', &mut raw_line)
            .map_err(|e| format!("Could not read annotation file '{}': {e}", path.display()))?;
        if line_bytes == 0 {
            break;
        }
        total_lines = total_lines.saturating_add(1);
        bytes_read_total = bytes_read_total.saturating_add(line_bytes as u64);
        if !on_progress(bytes_read_total, total_bytes) {
            return Err(prepare_cancelled_error("annotation parse progress"));
        }
        if raw_line.ends_with(b"\n") {
            raw_line.pop();
        }
        if raw_line.ends_with(b"\r") {
            raw_line.pop();
        }
        let line = match std::str::from_utf8(&raw_line) {
            Ok(v) => v.to_string(),
            Err(_) => {
                malformed_lines = malformed_lines.saturating_add(1);
                skipped_lines = skipped_lines.saturating_add(1);
                push_annotation_parse_issue(
                    &mut issue_examples,
                    total_lines,
                    "line is not valid UTF-8",
                    &String::from_utf8_lossy(&raw_line),
                );
                continue;
            }
        };
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            skipped_lines = skipped_lines.saturating_add(1);
            continue;
        }

        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 9 {
            malformed_lines = malformed_lines.saturating_add(1);
            skipped_lines = skipped_lines.saturating_add(1);
            push_annotation_parse_issue(
                &mut issue_examples,
                total_lines,
                "expected at least 9 tab-separated columns",
                trimmed,
            );
            continue;
        }

        let chromosome = cols[0].trim();
        if chromosome.is_empty() {
            malformed_lines = malformed_lines.saturating_add(1);
            skipped_lines = skipped_lines.saturating_add(1);
            push_annotation_parse_issue(
                &mut issue_examples,
                total_lines,
                "chromosome/contig column is empty",
                trimmed,
            );
            continue;
        }

        let feature_kind = cols[2];
        let start_1based = match parse_annotation_coordinate(cols[3]) {
            Some(v) => v,
            None => {
                malformed_lines = malformed_lines.saturating_add(1);
                skipped_lines = skipped_lines.saturating_add(1);
                push_annotation_parse_issue(
                    &mut issue_examples,
                    total_lines,
                    "invalid start coordinate",
                    trimmed,
                );
                continue;
            }
        };
        let end_1based = match parse_annotation_coordinate(cols[4]) {
            Some(v) => v,
            None => {
                malformed_lines = malformed_lines.saturating_add(1);
                skipped_lines = skipped_lines.saturating_add(1);
                push_annotation_parse_issue(
                    &mut issue_examples,
                    total_lines,
                    "invalid end coordinate",
                    trimmed,
                );
                continue;
            }
        };
        let strand_raw = cols[6];
        let attrs_raw = cols[8];

        if !is_gene_feature_kind(feature_kind) {
            skipped_lines = skipped_lines.saturating_add(1);
            continue;
        }
        if start_1based == 0 || end_1based < start_1based {
            malformed_lines = malformed_lines.saturating_add(1);
            skipped_lines = skipped_lines.saturating_add(1);
            push_annotation_parse_issue(
                &mut issue_examples,
                total_lines,
                "gene interval has start/end inconsistency",
                trimmed,
            );
            continue;
        }
        let attrs = parse_annotation_attributes(attrs_raw);
        let gene_name = pick_gene_name(&attrs);
        let mut gene_id = pick_gene_id(&attrs);
        let biotype = pick_gene_biotype(&attrs);
        if gene_name.is_none() && gene_id.is_none() {
            gene_id = Some(format!("{}_{}_{}", chromosome, start_1based, end_1based));
        }
        let strand = strand_raw.chars().next().and_then(|c| match c {
            '+' | '-' => Some(c),
            _ => None,
        });
        genes.push(GenomeGeneRecord {
            chromosome: chromosome.to_string(),
            start_1based,
            end_1based,
            strand,
            gene_id,
            gene_name,
            biotype,
        });
    }
    genes.sort_by(|a, b| {
        let a_name = a
            .gene_name
            .as_deref()
            .or(a.gene_id.as_deref())
            .unwrap_or(&a.chromosome);
        let b_name = b
            .gene_name
            .as_deref()
            .or(b.gene_id.as_deref())
            .unwrap_or(&b.chromosome);
        a_name
            .cmp(b_name)
            .then(a.chromosome.cmp(&b.chromosome))
            .then(a.start_1based.cmp(&b.start_1based))
            .then(a.end_1based.cmp(&b.end_1based))
    });
    if !on_progress(total_bytes.unwrap_or(bytes_read_total), total_bytes) {
        return Err(prepare_cancelled_error("annotation parse completion"));
    }
    let issue_examples_len = issue_examples.len();
    let report = AnnotationParseReport {
        path: canonical_or_display(path),
        source_format: "tabular".to_string(),
        total_lines,
        parsed_gene_records: genes.len(),
        skipped_lines,
        malformed_lines,
        issue_examples,
        truncated_issue_count: malformed_lines.saturating_sub(issue_examples_len),
    };
    Ok((genes, report))
}

fn parse_genbank_gene_records_with_progress<F>(
    path: &Path,
    mut on_progress: F,
) -> Result<(Vec<GenomeGeneRecord>, AnnotationParseReport), String>
where
    F: FnMut(u64, Option<u64>) -> bool,
{
    let total_bytes = fs::metadata(path).ok().map(|m| m.len());
    if !on_progress(0, total_bytes) {
        return Err(prepare_cancelled_error("GenBank annotation parse start"));
    }
    let parsed = gb_io::reader::parse_file(path.to_string_lossy().as_ref()).map_err(|e| {
        format!(
            "Could not parse GenBank annotation file '{}': {e}",
            path.display()
        )
    })?;
    let records = collect_genbank_like_gene_records(&parsed);
    let total_done = total_bytes.unwrap_or(0);
    if !on_progress(total_done, total_bytes) {
        return Err(prepare_cancelled_error(
            "GenBank annotation parse completion",
        ));
    }
    let report = AnnotationParseReport {
        path: canonical_or_display(path),
        source_format: "genbank".to_string(),
        total_lines: 0,
        parsed_gene_records: records.len(),
        skipped_lines: 0,
        malformed_lines: 0,
        issue_examples: vec![],
        truncated_issue_count: 0,
    };
    Ok((records, report))
}

fn parse_xml_gbseq_gene_records_with_progress<F>(
    path: &Path,
    mut on_progress: F,
) -> Result<(Vec<GenomeGeneRecord>, AnnotationParseReport), String>
where
    F: FnMut(u64, Option<u64>) -> bool,
{
    let total_bytes = fs::metadata(path).ok().map(|m| m.len());
    if !on_progress(0, total_bytes) {
        return Err(prepare_cancelled_error("XML annotation parse start"));
    }
    let parsed = parse_gbseq_xml_file(path.to_string_lossy().as_ref()).map_err(|e| {
        format!(
            "Could not parse XML annotation file '{}': {e}",
            path.display()
        )
    })?;
    let records = collect_genbank_like_gene_records(&parsed);
    let total_done = total_bytes.unwrap_or(0);
    if !on_progress(total_done, total_bytes) {
        return Err(prepare_cancelled_error("XML annotation parse completion"));
    }
    let report = AnnotationParseReport {
        path: canonical_or_display(path),
        source_format: "xml_gbseq".to_string(),
        total_lines: 0,
        parsed_gene_records: records.len(),
        skipped_lines: 0,
        malformed_lines: 0,
        issue_examples: vec![],
        truncated_issue_count: 0,
    };
    Ok((records, report))
}

fn collect_genbank_like_gene_records(parsed: &[gb_io::seq::Seq]) -> Vec<GenomeGeneRecord> {
    let mut records: Vec<GenomeGeneRecord> = vec![];
    for (record_idx, seq) in parsed.iter().enumerate() {
        let chromosome = genbank_record_chromosome(seq, record_idx);
        for feature in &seq.features {
            let kind = feature.kind.to_string();
            if !is_genbank_feature_kind_indexed(&kind) {
                continue;
            }
            let Ok((raw_from, raw_to)) = feature.location.find_bounds() else {
                continue;
            };
            if raw_from < 0 || raw_to < 0 {
                continue;
            }
            let mut start_0based = raw_from as usize;
            let mut end_0based = raw_to as usize;
            if end_0based < start_0based {
                std::mem::swap(&mut start_0based, &mut end_0based);
            }
            let attrs = genbank_feature_attributes(feature);
            let biotype = normalize_genbank_feature_biotype(&kind, &attrs);
            let gene_name = pick_gene_name(&attrs).or_else(|| {
                pick_annotation_attribute(
                    &attrs,
                    &[
                        "label",
                        "product",
                        "note",
                        "standard_name",
                        "bound_moiety",
                        "regulatory_class",
                    ],
                )
            });
            let mut gene_id = pick_gene_id(&attrs);
            if gene_id.is_none() {
                gene_id = Some(format!(
                    "{}:{}:{}-{}",
                    kind.to_ascii_lowercase(),
                    chromosome,
                    start_0based + 1,
                    end_0based + 1
                ));
            }
            let strand = Some(if feature_is_reverse(feature) {
                '-'
            } else {
                '+'
            });
            records.push(GenomeGeneRecord {
                chromosome: chromosome.clone(),
                start_1based: start_0based + 1,
                end_1based: end_0based + 1,
                strand,
                gene_id,
                gene_name,
                biotype,
            });
        }
    }
    records.sort_by(|a, b| {
        let a_name = a
            .gene_name
            .as_deref()
            .or(a.gene_id.as_deref())
            .unwrap_or(&a.chromosome);
        let b_name = b
            .gene_name
            .as_deref()
            .or(b.gene_id.as_deref())
            .unwrap_or(&b.chromosome);
        a_name
            .cmp(b_name)
            .then(a.chromosome.cmp(&b.chromosome))
            .then(a.start_1based.cmp(&b.start_1based))
            .then(a.end_1based.cmp(&b.end_1based))
    });
    records
}

fn push_annotation_parse_issue(
    issue_examples: &mut Vec<AnnotationParseIssue>,
    line: usize,
    reason: &str,
    context: &str,
) {
    if issue_examples.len() >= ANNOTATION_PARSE_ISSUE_LIMIT {
        return;
    }
    let mut normalized = context.replace('\t', " ");
    normalized = normalized.trim().to_string();
    let context = if normalized.chars().count() > ANNOTATION_PARSE_CONTEXT_CHARS {
        let clipped: String = normalized
            .chars()
            .take(ANNOTATION_PARSE_CONTEXT_CHARS)
            .collect();
        format!("{clipped}...")
    } else {
        normalized
    };
    issue_examples.push(AnnotationParseIssue {
        line,
        reason: reason.to_string(),
        context,
    });
}

fn genbank_record_chromosome(seq: &gb_io::seq::Seq, record_idx: usize) -> String {
    let from_version = seq
        .version
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty())
        .map(str::to_string);
    let from_accession = seq
        .accession
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty())
        .map(str::to_string);
    let from_name = seq
        .name
        .as_ref()
        .map(|v| v.trim())
        .filter(|v| !v.is_empty())
        .map(str::to_string);
    from_version
        .or(from_accession)
        .or(from_name)
        .unwrap_or_else(|| format!("record_{}", record_idx + 1))
}

fn is_genbank_feature_kind_indexed(raw: &str) -> bool {
    let lower = raw.trim().to_ascii_lowercase();
    if lower.is_empty() {
        return false;
    }
    if lower == "source" {
        return false;
    }
    true
}

fn genbank_feature_attributes(feature: &gb_io::seq::Feature) -> HashMap<String, String> {
    let mut attrs = HashMap::new();
    for (key, value) in &feature.qualifiers {
        let key = key.trim().to_ascii_lowercase();
        if key.is_empty() {
            continue;
        }
        let Some(raw) = value.as_ref() else {
            continue;
        };
        if let Some(normalized) = normalize_annotation_attribute_value(raw) {
            attrs.entry(key).or_insert(normalized);
        }
    }
    attrs
}

fn normalize_genbank_feature_biotype(
    feature_kind: &str,
    attrs: &HashMap<String, String>,
) -> Option<String> {
    let kind = feature_kind.trim().to_ascii_lowercase();
    if kind.is_empty() {
        return None;
    }
    if kind == "regulatory" {
        if let Some(class) = attrs.get("regulatory_class").map(|v| v.trim()) {
            if !class.is_empty() {
                return Some(class.to_ascii_lowercase());
            }
        }
    }
    if kind == "rep_origin" {
        return Some("origin_of_replication".to_string());
    }
    if kind == "misc_feature" {
        for key in ["label", "note"] {
            if let Some(value) = attrs.get(key) {
                let lower = value.to_ascii_lowercase();
                if lower.contains("ltr") {
                    return Some("ltr".to_string());
                }
                if lower.contains("itr") {
                    return Some("itr".to_string());
                }
                if lower.contains("origin") || lower.contains(" ori") || lower.starts_with("ori") {
                    return Some("origin_of_replication".to_string());
                }
                if lower.contains("promoter") {
                    return Some("promoter".to_string());
                }
                if lower.contains("marker") || lower.contains("selection") {
                    return Some("selectable_marker".to_string());
                }
            }
        }
    }
    Some(kind)
}

fn is_gene_feature_kind(raw: &str) -> bool {
    let lower = raw.trim().to_ascii_lowercase();
    lower == "gene" || lower == "pseudogene" || lower.ends_with("_gene")
}

fn parse_annotation_coordinate(raw: &str) -> Option<usize> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let stripped = trimmed.trim_start_matches(['<', '>', '=']);
    if stripped.is_empty() {
        return None;
    }
    let no_commas = stripped.replace(',', "");
    if let Ok(v) = no_commas.parse::<usize>() {
        return (v > 0).then_some(v);
    }
    let prefix_digits: String = no_commas
        .chars()
        .take_while(|ch| ch.is_ascii_digit())
        .collect();
    if prefix_digits.is_empty() {
        None
    } else {
        prefix_digits.parse::<usize>().ok().filter(|v| *v > 0)
    }
}

fn normalize_transcript_id(raw: &str) -> String {
    raw.trim().to_string()
}

fn normalize_gene_match_token(raw: &str) -> String {
    let trimmed = raw.trim();
    let core = trimmed
        .split_once('.')
        .map(|(left, _)| left)
        .unwrap_or(trimmed);
    core.to_ascii_uppercase()
}

fn normalize_chromosome_token(raw: &str) -> String {
    let trimmed = raw.trim().to_ascii_lowercase();
    let without_chr = trimmed.strip_prefix("chr").unwrap_or(&trimmed);
    let stripped_zeros = without_chr.trim_start_matches('0');
    if stripped_zeros.is_empty() {
        "0".to_string()
    } else {
        stripped_zeros.to_string()
    }
}

fn chromosome_names_match(expected: &str, actual: &str) -> bool {
    if expected.eq_ignore_ascii_case(actual) {
        return true;
    }
    normalize_chromosome_token(expected) == normalize_chromosome_token(actual)
}

fn pick_gene_name(attrs: &HashMap<String, String>) -> Option<String> {
    pick_annotation_attribute(
        attrs,
        &[
            "gene_name",
            "name",
            "gene_symbol",
            "symbol",
            "gene",
            "locus_tag",
            "standard_name",
        ],
    )
}

fn pick_gene_id(attrs: &HashMap<String, String>) -> Option<String> {
    let direct = pick_annotation_attribute(
        attrs,
        &[
            "gene_id",
            "id",
            "gene",
            "locus_tag",
            "geneid",
            "transcript_id",
        ],
    );
    direct
        .or_else(|| pick_gene_id_from_dbxref(attrs))
        .map(normalize_gene_id)
}

fn pick_gene_biotype(attrs: &HashMap<String, String>) -> Option<String> {
    pick_annotation_attribute(
        attrs,
        &[
            "gene_biotype",
            "gene_type",
            "biotype",
            "transcript_biotype",
            "transcript_type",
        ],
    )
}

fn pick_gene_id_from_dbxref(attrs: &HashMap<String, String>) -> Option<String> {
    let raw = attrs
        .get("dbxref")
        .or_else(|| attrs.get("db_xref"))
        .or_else(|| attrs.get("xref"))?;
    for token in raw.split(',').map(str::trim).filter(|v| !v.is_empty()) {
        if let Some((key, value)) = token.split_once(':') {
            let k = key.trim().to_ascii_lowercase();
            let v = value.trim();
            if v.is_empty() {
                continue;
            }
            if matches!(
                k.as_str(),
                "geneid" | "gene" | "ensembl_gene_id" | "ensemblgene" | "hgnc" | "locus_tag"
            ) {
                return Some(v.to_string());
            }
        } else {
            return Some(token.to_string());
        }
    }
    None
}

fn normalize_gene_id(raw: String) -> String {
    let trimmed = raw.trim();
    if let Some((prefix, rest)) = trimmed.split_once(':') {
        let p = prefix.trim().to_ascii_lowercase();
        if matches!(p.as_str(), "gene" | "geneid" | "id") {
            let rest = rest.trim();
            if !rest.is_empty() {
                return rest.to_string();
            }
        }
    }
    trimmed.to_string()
}

fn pick_annotation_attribute(attrs: &HashMap<String, String>, keys: &[&str]) -> Option<String> {
    keys.iter()
        .filter_map(|key| attrs.get(*key))
        .find(|value| !value.trim().is_empty())
        .cloned()
}

fn split_annotation_attributes(raw: &str) -> Vec<String> {
    let mut fields: Vec<String> = vec![];
    let mut current = String::new();
    let mut in_quotes = false;
    let mut escape = false;
    for ch in raw.chars() {
        if escape {
            current.push(ch);
            escape = false;
            continue;
        }
        match ch {
            '\\' if in_quotes => {
                current.push(ch);
                escape = true;
            }
            '"' => {
                in_quotes = !in_quotes;
                current.push(ch);
            }
            ';' if !in_quotes => {
                let field = current.trim();
                if !field.is_empty() {
                    fields.push(field.to_string());
                }
                current.clear();
            }
            _ => current.push(ch),
        }
    }
    let tail = current.trim();
    if !tail.is_empty() {
        fields.push(tail.to_string());
    }
    fields
}

fn normalize_annotation_attribute_value(raw: &str) -> Option<String> {
    let mut cleaned = raw
        .trim()
        .trim_matches('"')
        .trim()
        .trim_end_matches(',')
        .trim()
        .to_string();
    cleaned = cleaned
        .replace("\\\"", "\"")
        .replace("\\;", ";")
        .replace("\\,", ",");
    cleaned = cleaned.trim().trim_matches('"').trim().to_string();
    if cleaned.is_empty() {
        return None;
    }
    if cleaned.is_empty() {
        None
    } else {
        Some(cleaned)
    }
}

fn parse_annotation_attributes(raw: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for part in split_annotation_attributes(raw) {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((key, value)) = part.split_once('=') {
            let key = key.trim().to_ascii_lowercase();
            if !key.is_empty() {
                if let Some(value) = normalize_annotation_attribute_value(value) {
                    map.entry(key).or_insert(value);
                }
            }
            continue;
        }
        let mut pieces = part.splitn(2, char::is_whitespace);
        let key = pieces
            .next()
            .unwrap_or_default()
            .trim()
            .to_ascii_lowercase();
        let value = pieces.next().unwrap_or_default();
        if !key.is_empty() {
            if let Some(value) = normalize_annotation_attribute_value(value) {
                map.entry(key).or_insert(value);
            }
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{Compression, write::GzEncoder};
    use std::env;
    #[cfg(unix)]
    use std::os::unix::fs::PermissionsExt;
    use tempfile::tempdir;

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

    fn write_gzip(path: &Path, text: &str) {
        let file = File::create(path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(text.as_bytes()).unwrap();
        encoder.finish().unwrap();
    }

    fn write_multi_member_gzip(path: &Path, members: &[&str]) {
        File::create(path).unwrap();
        for member in members {
            let file = OpenOptions::new().append(true).open(path).unwrap();
            let mut encoder = GzEncoder::new(file, Compression::default());
            encoder.write_all(member.as_bytes()).unwrap();
            encoder.finish().unwrap();
        }
    }

    fn file_url(path: &Path) -> String {
        format!("file://{}", path.display())
    }

    #[cfg(unix)]
    fn write_executable_script(path: &Path, body: &str) {
        fs::write(path, body).expect("write executable script");
        let mut perms = fs::metadata(path).expect("script metadata").permissions();
        perms.set_mode(0o755);
        fs::set_permissions(path, perms).expect("set executable permissions");
    }

    #[test]
    fn test_blast_external_binary_preflight_reports_missing_tools() {
        let report = blast_external_binary_preflight_report_with_executables(
            "__gentle_missing_blastn__".to_string(),
            "__gentle_missing_makeblastdb__".to_string(),
        );
        assert_eq!(report.schema, "gentle.blast_external_binary_preflight.v1");
        assert_eq!(report.blastn.tool, "blastn");
        assert_eq!(report.makeblastdb.tool, "makeblastdb");
        assert!(!report.blastn.found);
        assert!(!report.makeblastdb.found);
        assert_eq!(report.blastn.executable, "__gentle_missing_blastn__");
        assert_eq!(
            report.makeblastdb.executable,
            "__gentle_missing_makeblastdb__"
        );
    }

    #[cfg(unix)]
    #[test]
    fn test_blast_external_binary_preflight_reports_version_and_path_for_configured_tools() {
        let td = tempdir().expect("tempdir");
        let fake_blastn = td.path().join("fake_blastn");
        let fake_makeblastdb = td.path().join("fake_makeblastdb");
        write_executable_script(
            &fake_blastn,
            "#!/bin/sh\necho \"blastn: fake 1.2.3\"; exit 0\n",
        );
        write_executable_script(
            &fake_makeblastdb,
            "#!/bin/sh\necho \"makeblastdb: fake 4.5.6\"; exit 0\n",
        );

        let report = blast_external_binary_preflight_report_with_executables(
            fake_blastn.to_string_lossy().to_string(),
            fake_makeblastdb.to_string_lossy().to_string(),
        );
        assert!(report.blastn.found);
        assert!(report.makeblastdb.found);
        assert!(report.blastn.version_probe_ok);
        assert!(report.makeblastdb.version_probe_ok);
        assert!(
            report
                .blastn
                .version
                .as_deref()
                .unwrap_or_default()
                .contains("fake 1.2.3")
        );
        assert!(
            report
                .makeblastdb
                .version
                .as_deref()
                .unwrap_or_default()
                .contains("fake 4.5.6")
        );
        assert!(report.blastn.resolved_path.is_some());
        assert!(report.makeblastdb.resolved_path.is_some());
    }

    #[test]
    fn test_prepare_genome_once_and_extract_region() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta_gz = root.join("toy.fa.gz");
        let ann_gz = root.join("toy.gtf.gz");
        write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\nACGT\n>II\nTTTT\nGGGG\n");
        write_gzip(
            &ann_gz,
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; gene_biotype \"protein_coding\";\n",
        );

        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
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

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let first = catalog.prepare_genome_once("ToyGenome").unwrap();
        assert!(!first.reused_existing);
        assert!(Path::new(&first.sequence_path).exists());
        assert!(Path::new(&first.annotation_path).exists());
        assert!(Path::new(&first.fasta_index_path).exists());
        assert_eq!(first.cached_contig_count, 2);
        assert_eq!(first.cached_total_span_bp, 20);
        assert_eq!(first.cached_longest_contig.as_deref(), Some("chr1"));
        assert_eq!(first.cached_longest_contig_bp, Some(12));
        let manifest =
            GenomeCatalog::load_manifest(&cache_dir.join("toygenome").join("manifest.json"))
                .unwrap();
        assert!(
            manifest
                .sequence_sha1
                .as_ref()
                .map(|v| !v.is_empty())
                .unwrap_or(false)
        );
        assert!(
            manifest
                .annotation_sha1
                .as_ref()
                .map(|v| !v.is_empty())
                .unwrap_or(false)
        );
        let gene_index_path = manifest
            .gene_index_path
            .expect("gene index path should be set");
        assert!(Path::new(&gene_index_path).exists());

        let inspection = catalog
            .inspect_prepared_genome("ToyGenome", None)
            .unwrap()
            .expect("inspection should exist");
        assert!(inspection.sequence_present);
        assert!(inspection.annotation_present);
        assert!(inspection.fasta_index_ready);
        assert!(inspection.gene_index_ready);
        assert_eq!(inspection.cached_contig_count, 2);
        assert_eq!(inspection.cached_total_span_bp, 20);
        assert_eq!(inspection.cached_longest_contig.as_deref(), Some("chr1"));
        assert_eq!(inspection.cached_longest_contig_bp, Some(12));
        assert_eq!(
            inspection.cached_contig_preview,
            vec!["chr1".to_string(), "II".to_string()]
        );
        assert!(
            inspection
                .sequence_sha1
                .as_ref()
                .map(|v| !v.is_empty())
                .unwrap_or(false)
        );
        assert!(
            inspection
                .annotation_sha1
                .as_ref()
                .map(|v| !v.is_empty())
                .unwrap_or(false)
        );

        let second = catalog.prepare_genome_once("ToyGenome").unwrap();
        assert!(second.reused_existing);
        assert_eq!(second.cached_contig_count, 2);
        assert_eq!(second.cached_total_span_bp, 20);
        assert_eq!(second.cached_longest_contig.as_deref(), Some("chr1"));
        assert_eq!(second.cached_longest_contig_bp, Some(12));

        let seq = catalog
            .get_sequence_region("ToyGenome", "1", 3, 10)
            .unwrap();
        assert_eq!(seq, "GTACGTAC");

        let genes = catalog.list_gene_regions("ToyGenome", None).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].chromosome, "chr1");
        assert_eq!(genes[0].start_1based, 1);
        assert_eq!(genes[0].end_1based, 12);
        assert_eq!(genes[0].strand, Some('+'));
        assert_eq!(genes[0].gene_id.as_deref(), Some("GENE1"));
        assert_eq!(genes[0].gene_name.as_deref(), Some("MYGENE"));
        assert_eq!(genes[0].biotype.as_deref(), Some("protein_coding"));
    }

    #[test]
    fn test_prepare_genome_reports_cache_summary_for_worm_like_lettered_contigs() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta_gz = root.join("celegans.fa.gz");
        let ann_gz = root.join("celegans.gtf.gz");
        write_gzip(
            &fasta_gz,
            ">I\nACGTACGTACGT\n>II\nTTTTGGGGCCCCAAAATTTT\n>X\nGATTACAGATTACAGATTACA\n>MtDNA\nATGCATGC\n",
        );
        write_gzip(
            &ann_gz,
            concat!(
                "X\tsrc\tgene\t3\t18\t.\t+\t.\tgene_id \"WBGene00000001\"; gene_name \"lin-4\"; gene_biotype \"ncRNA\";\n",
                "II\tsrc\tgene\t5\t16\t.\t-\t.\tgene_id \"WBGene00000002\"; gene_name \"unc-54\"; gene_biotype \"protein_coding\";\n"
            ),
        );

        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let genome_id = "Caenorhabditis elegans WBcel235 Ensembl 115";
        let catalog_json = format!(
            r#"{{
  "{genome_id}": {{
    "description": "worm test genome",
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

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        let report = catalog.prepare_genome_once(genome_id).unwrap();
        assert_eq!(report.cached_contig_count, 4);
        assert_eq!(report.cached_total_span_bp, 61);
        assert_eq!(report.cached_longest_contig.as_deref(), Some("X"));
        assert_eq!(report.cached_longest_contig_bp, Some(21));
        assert_eq!(
            report.cached_contig_preview,
            vec![
                "X".to_string(),
                "II".to_string(),
                "I".to_string(),
                "MtDNA".to_string()
            ]
        );

        let seq = catalog.get_sequence_region(genome_id, "X", 5, 10).unwrap();
        assert_eq!(seq, "ACAGAT");
        let inspection = catalog
            .inspect_prepared_genome(genome_id, None)
            .unwrap()
            .expect("inspection should exist");
        assert_eq!(inspection.cached_contig_count, 4);
        assert_eq!(inspection.cached_total_span_bp, 61);
        assert_eq!(inspection.cached_longest_contig.as_deref(), Some("X"));
        assert_eq!(inspection.cached_longest_contig_bp, Some(21));
        let genes = catalog.list_gene_regions(genome_id, None).unwrap();
        assert!(genes.iter().any(|gene| gene.chromosome == "X"));
    }

    #[test]
    fn test_prepare_genome_reads_all_members_from_concatenated_gzip_sources() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta_gz = root.join("celegans_concat.fa.gz");
        let ann_gz = root.join("celegans_concat.gtf.gz");
        write_multi_member_gzip(
            &fasta_gz,
            &[
                ">I\nACGTACGTACGT\n",
                ">III\nTTTTCCCCAAAA\n",
                ">X\nGATTACA\n",
            ],
        );
        write_multi_member_gzip(
            &ann_gz,
            &[
                "III\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"WBGene00010001\"; gene_name \"abc-1\"; gene_biotype \"protein_coding\";\n",
            ],
        );

        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let genome_id = "Caenorhabditis elegans concat test";
        let catalog_json = format!(
            r#"{{
  "{genome_id}": {{
    "description": "worm concatenated gzip test genome",
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

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        let report = catalog.prepare_genome_once(genome_id).unwrap();
        assert_eq!(report.cached_contig_count, 3);
        assert_eq!(report.cached_total_span_bp, 31);
        assert_eq!(
            report.cached_contig_preview,
            vec!["I".to_string(), "III".to_string(), "X".to_string()]
        );

        let seq = catalog.get_sequence_region(genome_id, "III", 3, 8).unwrap();
        assert_eq!(seq, "TTCCCC");
        let genes = catalog.list_gene_regions(genome_id, None).unwrap();
        assert!(genes.iter().any(|gene| gene.chromosome == "III"));
    }

    #[test]
    fn test_prepare_genome_rejects_annotation_contigs_missing_from_sequence_cache() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta = root.join("toy.fa");
        let ann = root.join("toy.gtf");
        fs::write(&fasta, ">I\nACGTACGTACGT\n").unwrap();
        fs::write(
            &ann,
            "III\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"WBGene00020001\"; gene_name \"def-1\"; gene_biotype \"protein_coding\";\n",
        )
        .unwrap();

        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let genome_id = "Synthetic mismatch genome";
        let catalog_json = format!(
            r#"{{
  "{genome_id}": {{
    "description": "synthetic mismatch genome",
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

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        let err = catalog
            .prepare_genome_once(genome_id)
            .expect_err("mismatched annotation contigs should fail prepare");
        assert!(err.contains("annotation gene index"));
        assert!(err.contains("III"));
        assert!(err.contains("sequence.fa"));
        assert!(err.contains("sequence.fa.fai"));
    }

    #[test]
    fn test_list_gene_transcript_records_from_gtf() {
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
                "chr1\tsrc\ttranscript\t2\t10\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX2\";\n",
                "chr1\tsrc\texon\t2\t10\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX2\"; exon_number \"1\";\n",
                "chr1\tsrc\tCDS\t3\t9\t.\t+\t0\tgene_id \"GENE1\"; gene_name \"MYGENE\"; transcript_id \"TX2\";\n",
                "chr1\tsrc\tgene\t14\t18\t.\t+\t.\tgene_id \"GENE2\"; gene_name \"OTHER\";\n",
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

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog.prepare_genome_once("ToyGenome").unwrap();

        let records = catalog
            .list_gene_transcript_records(
                "ToyGenome",
                "chr1",
                1,
                12,
                Some("GENE1"),
                Some("MYGENE"),
                None,
            )
            .unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].transcript_id, "TX1");
        assert_eq!(records[0].exons_1based, vec![(1, 4), (9, 12)]);
        assert_eq!(records[0].cds_1based, vec![(2, 4), (9, 11)]);
        assert_eq!(records[1].transcript_id, "TX2");
        assert_eq!(records[1].exons_1based, vec![(2, 10)]);
        assert_eq!(records[1].cds_1based, vec![(3, 9)]);
    }

    #[test]
    fn test_get_sequence_region_missing_contig_reports_aliases_and_available_contigs() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta = root.join("toy.fa");
        let ann = root.join("toy.gtf");
        fs::write(&fasta, ">chr1\nACGTACGT\n>chrM\nTTTT\n").unwrap();
        fs::write(
            &ann,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
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

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog.prepare_genome_once("ToyGenome").unwrap();

        let err = catalog
            .get_sequence_region("ToyGenome", "17", 1, 4)
            .expect_err("missing contig should fail");
        assert!(err.contains("Chromosome/contig '17' not found in genome 'ToyGenome'"));
        assert!(err.contains("Tried aliases:"));
        assert!(err.contains("chr17"));
        assert!(err.contains("Available contigs"));
        assert!(err.contains("chr1"));
        assert!(err.contains("PrepareGenome"));
        assert!(err.contains("Prepared sequence='"));
        assert!(err.contains("FASTA index='"));
    }

    #[test]
    fn test_get_sequence_region_missing_contig_reports_suggestions_for_accession_style_names() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta = root.join("toy.fa");
        let ann = root.join("toy.gtf");
        fs::write(&fasta, ">NC_000017.11\nACGTACGT\n").unwrap();
        fs::write(
            &ann,
            "17\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
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

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog.prepare_genome_once("ToyGenome").unwrap();

        let err = catalog
            .get_sequence_region("ToyGenome", "17", 1, 4)
            .expect_err("missing contig should fail");
        assert!(err.contains("Suggested matching contigs"));
        assert!(err.contains("NC_000017.11"));
    }

    #[test]
    fn test_get_sequence_region_missing_contig_does_not_suggest_shorter_numeric_prefix() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta = root.join("toy.fa");
        let ann = root.join("toy.gtf");
        fs::write(&fasta, ">1\nACGTACGT\n").unwrap();
        fs::write(
            &ann,
            "1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
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

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog.prepare_genome_once("ToyGenome").unwrap();

        let err = catalog
            .get_sequence_region("ToyGenome", "17", 1, 4)
            .expect_err("missing contig should fail");
        assert!(err.contains("Available contigs (1): 1"));
        assert!(!err.contains("Suggested matching contigs"));
    }

    #[test]
    fn test_list_chromosome_lengths_reads_prepared_fasta_index() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta = root.join("toy.fa");
        let ann = root.join("toy.gtf");
        fs::write(&fasta, ">chr1\nACGTACGT\n>chrM\nTTTT\n").unwrap();
        fs::write(
            &ann,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
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

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog.prepare_genome_once("ToyGenome").unwrap();

        let chromosomes = catalog.list_chromosome_lengths("ToyGenome", None).unwrap();
        assert_eq!(chromosomes.len(), 2);
        assert_eq!(chromosomes[0].chromosome, "chr1");
        assert_eq!(chromosomes[0].length_bp, 8);
        assert_eq!(chromosomes[1].chromosome, "chrM");
        assert_eq!(chromosomes[1].length_bp, 4);
    }

    #[test]
    fn test_prepare_fails_when_annotation_source_missing() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta_gz = root.join("toy.fa.gz");
        write_gzip(&fasta_gz, ">chr1\nACGT\n");

        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_remote": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            file_url(&fasta_gz),
            root.join("cache").display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();

        let err = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap_err();
        assert!(err.contains("annotation"));
    }

    #[test]
    fn test_prepare_fails_fast_when_annotation_path_is_invalid_before_sequence_download() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta_gz = root.join("toy.fa.gz");
        write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\n");

        let cache_dir = root.join("cache");
        let missing_annotation = root.join("missing_annotation.gtf");
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_remote": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            file_url(&fasta_gz),
            missing_annotation.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();

        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();

        let mut first_phases: Vec<String> = vec![];
        let first_err = catalog
            .prepare_genome_once_with_progress("ToyGenome", None, &mut |progress| {
                first_phases.push(progress.phase);
                true
            })
            .unwrap_err();
        assert!(first_err.contains("missing_annotation.gtf"));
        assert!(
            !first_phases
                .iter()
                .any(|phase| phase == "download_sequence")
        );
        assert!(first_phases.is_empty());

        let sequence_path = cache_dir.join("toygenome").join("sequence.fa");
        assert!(!sequence_path.exists());

        let mut second_phases: Vec<String> = vec![];
        let second_err = catalog
            .prepare_genome_once_with_progress("ToyGenome", None, &mut |progress| {
                second_phases.push(progress.phase);
                true
            })
            .unwrap_err();
        assert!(second_err.contains("missing_annotation.gtf"));
        assert!(
            !second_phases
                .iter()
                .any(|phase| phase == "download_sequence")
        );
        assert!(second_phases.is_empty());
    }

    #[test]
    fn test_prepare_second_run_reuses_sequence_and_annotation_without_redownload() {
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
    "description": "toy test genome",
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
        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );

        let mut first_phases: Vec<String> = vec![];
        let first = catalog
            .prepare_genome_once_with_progress("ToyGenome", None, &mut |progress| {
                first_phases.push(progress.phase);
                true
            })
            .unwrap();
        assert!(!first.reused_existing);
        assert!(
            first_phases
                .iter()
                .any(|phase| phase == "download_sequence")
        );
        assert!(
            first_phases
                .iter()
                .any(|phase| phase == "download_annotation")
        );

        let mut second_phases: Vec<String> = vec![];
        let second = catalog
            .prepare_genome_once_with_progress("ToyGenome", None, &mut |progress| {
                second_phases.push(progress.phase);
                true
            })
            .unwrap();
        assert!(second.reused_existing);
        assert!(
            !second_phases
                .iter()
                .any(|phase| phase == "download_sequence")
        );
        assert!(
            !second_phases
                .iter()
                .any(|phase| phase == "download_annotation")
        );
        assert!(second_phases.iter().any(|phase| phase == "index_blast"));
        assert!(second_phases.iter().any(|phase| phase == "ready"));
    }

    #[test]
    fn test_prepare_rebuilds_cache_when_catalog_sources_change() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta_a = root.join("toy_a.fa");
        let ann_a = root.join("toy_a.gtf");
        fs::write(&fasta_a, ">chr1\nACGTACGT\n").unwrap();
        fs::write(
            &ann_a,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
        )
        .unwrap();

        let fasta_b = root.join("toy_b.fa");
        let ann_b = root.join("toy_b.gtf");
        fs::write(&fasta_b, ">chr1\nTTTTGGGG\n").unwrap();
        fs::write(
            &ann_b,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE2\"; gene_name \"GENE2\";\n",
        )
        .unwrap();

        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let catalog_a_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            fasta_a.display(),
            ann_a.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_a_json).unwrap();
        let catalog_a = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        let first = catalog_a.prepare_genome_once("ToyGenome").unwrap();
        assert!(!first.reused_existing);

        let catalog_b_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            fasta_b.display(),
            ann_b.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_b_json).unwrap();
        let catalog_b = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();

        let mut second_phases: Vec<String> = vec![];
        let second = catalog_b
            .prepare_genome_once_with_progress("ToyGenome", None, &mut |progress| {
                second_phases.push(progress.phase);
                true
            })
            .unwrap();
        assert!(!second.reused_existing);
        assert!(
            second_phases
                .iter()
                .any(|phase| phase == "download_sequence")
        );
        assert!(
            second_phases
                .iter()
                .any(|phase| phase == "download_annotation")
        );

        let installed_sequence =
            fs::read_to_string(cache_dir.join("toygenome").join("sequence.fa"))
                .expect("installed sequence");
        assert!(installed_sequence.contains("TTTTGGGG"));
        let installed_annotation =
            fs::read_to_string(cache_dir.join("toygenome").join("annotation.gtf"))
                .expect("installed annotation");
        assert!(installed_annotation.contains("GENE2"));
    }

    #[test]
    fn test_prepare_first_run_reuses_existing_annotation_and_gene_index_in_install_dir() {
        let td = tempdir().unwrap();
        let root = td.path();

        let source_fasta = root.join("source.fa");
        let source_ann = root.join("source.gtf");
        fs::write(&source_fasta, ">chr1\nACGT\nACGT\n").unwrap();
        fs::write(
            &source_ann,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
        )
        .unwrap();

        let cache_dir = root.join("cache");
        let install_dir = cache_dir.join("toygenome");
        fs::create_dir_all(&install_dir).unwrap();
        fs::copy(&source_fasta, install_dir.join("sequence.fa")).unwrap();
        fs::copy(&source_ann, install_dir.join("annotation.gtf")).unwrap();
        fs::write(install_dir.join("genes.json"), "[]").unwrap();

        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            source_fasta.display(),
            source_ann.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );

        let mut phases: Vec<String> = vec![];
        let report = catalog
            .prepare_genome_once_with_progress("ToyGenome", None, &mut |progress| {
                phases.push(progress.phase);
                true
            })
            .unwrap();
        assert!(!report.reused_existing);
        assert!(phases.iter().any(|phase| phase == "reuse_sequence"));
        assert!(phases.iter().any(|phase| phase == "reuse_annotation"));
        assert!(phases.iter().any(|phase| phase == "reuse_gene_index"));
        assert!(!phases.iter().any(|phase| phase == "download_sequence"));
        assert!(!phases.iter().any(|phase| phase == "download_annotation"));
        assert!(!phases.iter().any(|phase| phase == "index_genes"));
        assert!(
            report.annotation_parse_report.is_none(),
            "existing gene index should avoid reparsing annotation"
        );
    }

    #[test]
    fn test_reinstall_replaces_stale_install_dir_contents_from_sources() {
        let td = tempdir().unwrap();
        let root = td.path();

        let source_fasta = root.join("source.fa");
        let source_ann = root.join("source.gtf");
        fs::write(&source_fasta, ">chr1\nACGTACGT\n>chr17\nTTTTCCCC\n").unwrap();
        fs::write(
            &source_ann,
            concat!(
                "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"ONE\";\n",
                "chr17\tsrc\tgene\t1\t8\t.\t-\t.\tgene_id \"GENE17\"; gene_name \"SEVENTEEN\";\n"
            ),
        )
        .unwrap();

        let cache_dir = root.join("cache");
        let install_dir = cache_dir.join("toygenome");
        let blast_dir = install_dir.join("blastdb");
        fs::create_dir_all(&blast_dir).unwrap();
        fs::write(install_dir.join("sequence.fa"), ">chr1\nAAAA\n").unwrap();
        fs::write(
            install_dir.join("annotation.gtf"),
            "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"STALE\"; gene_name \"STALE\";\n",
        )
        .unwrap();
        fs::write(install_dir.join("genes.json"), "[]").unwrap();
        fs::write(blast_dir.join("genome.nhr"), "stale").unwrap();
        fs::write(blast_dir.join("genome.nin"), "stale").unwrap();
        fs::write(blast_dir.join("genome.nsq"), "stale").unwrap();

        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            source_fasta.display(),
            source_ann.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );

        let mut phases: Vec<String> = vec![];
        let report = catalog
            .reinstall_genome_once_with_progress("ToyGenome", None, &mut |progress| {
                phases.push(progress.phase);
                true
            })
            .unwrap();

        assert!(!report.reused_existing);
        assert!(phases.iter().any(|phase| phase == "download_sequence"));
        assert!(phases.iter().any(|phase| phase == "download_annotation"));
        assert!(phases.iter().any(|phase| phase == "index_genes"));
        assert!(!phases.iter().any(|phase| phase == "reuse_sequence"));
        assert!(!phases.iter().any(|phase| phase == "reuse_annotation"));
        assert!(!phases.iter().any(|phase| phase == "reuse_gene_index"));

        let chr17 = catalog
            .get_sequence_region_with_cache("ToyGenome", "17", 1, 4, None)
            .expect("chr17 should be available after reinstall");
        assert_eq!(chr17, "TTTT");

        assert!(!blast_dir.join("genome.nhr").exists());
        let installed_annotation =
            fs::read_to_string(cache_dir.join("toygenome").join("annotation.gtf"))
                .expect("installed annotation");
        assert!(installed_annotation.contains("GENE17"));
    }

    #[test]
    fn test_reindex_reuses_cached_files_without_redownloading_sources() {
        let td = tempdir().unwrap();
        let root = td.path();

        let source_fasta = root.join("source.fa");
        let source_ann = root.join("source.gtf");
        fs::write(&source_fasta, ">chr1\nACGTACGT\n").unwrap();
        fs::write(
            &source_ann,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"ONE\";\n",
        )
        .unwrap();

        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            source_fasta.display(),
            source_ann.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );

        let first = catalog.prepare_genome_once("ToyGenome").unwrap();
        assert!(!first.reused_existing);
        assert_eq!(
            catalog.get_sequence_region("ToyGenome", "chr1", 1, 4).unwrap(),
            "ACGT"
        );

        fs::write(&source_fasta, ">chr1\nTTTTTTTT\n>chr17\nCCCCAAAA\n").unwrap();
        fs::write(
            &source_ann,
            concat!(
                "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1B\"; gene_name \"ONE_B\";\n",
                "chr17\tsrc\tgene\t1\t8\t.\t-\t.\tgene_id \"GENE17\"; gene_name \"SEVENTEEN\";\n"
            ),
        )
        .unwrap();

        let mut phases: Vec<String> = vec![];
        let report = catalog
            .reindex_genome_once_with_progress("ToyGenome", None, &mut |progress| {
                phases.push(progress.phase);
                true
            })
            .unwrap();
        assert!(report.reused_existing);
        assert!(phases.iter().any(|phase| phase == "reset_indexes"));
        assert!(phases.iter().any(|phase| phase == "reuse_sequence"));
        assert!(phases.iter().any(|phase| phase == "reuse_annotation"));
        assert!(phases.iter().any(|phase| phase == "index_fasta"));
        assert!(phases.iter().any(|phase| phase == "index_genes"));
        assert!(!phases.iter().any(|phase| phase == "download_sequence"));
        assert!(!phases.iter().any(|phase| phase == "download_annotation"));
        assert!(!phases.iter().any(|phase| phase == "reuse_gene_index"));

        assert_eq!(
            catalog.get_sequence_region("ToyGenome", "chr1", 1, 4).unwrap(),
            "ACGT"
        );
        assert!(
            catalog.get_sequence_region("ToyGenome", "chr17", 1, 4).is_err(),
            "reindex should not have pulled in the changed source FASTA"
        );
        let genes = catalog.list_gene_regions("ToyGenome", None).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].gene_id.as_deref(), Some("GENE1"));
    }

    #[test]
    fn test_reindex_keeps_cached_files_even_when_catalog_sources_drift() {
        let td = tempdir().unwrap();
        let root = td.path();

        let source_fasta_v1 = root.join("source_v1.fa");
        let source_ann_v1 = root.join("source_v1.gtf");
        fs::write(&source_fasta_v1, ">chr1\nACGTACGT\n").unwrap();
        fs::write(
            &source_ann_v1,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"ONE\";\n",
        )
        .unwrap();

        let cache_dir = root.join("cache");
        let catalog_path_v1 = root.join("catalog_v1.json");
        let catalog_json_v1 = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            source_fasta_v1.display(),
            source_ann_v1.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path_v1, catalog_json_v1).unwrap();
        let catalog_v1 = GenomeCatalog::from_json_file(&catalog_path_v1.to_string_lossy()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );

        let first = catalog_v1.prepare_genome_once("ToyGenome").unwrap();
        assert!(!first.reused_existing);
        assert_eq!(
            catalog_v1.get_sequence_region("ToyGenome", "chr1", 1, 4).unwrap(),
            "ACGT"
        );

        let source_fasta_v2 = root.join("source_v2.fa");
        let source_ann_v2 = root.join("source_v2.gtf");
        fs::write(&source_fasta_v2, ">chr1\nTTTTTTTT\n>chr17\nCCCCAAAA\n").unwrap();
        fs::write(
            &source_ann_v2,
            concat!(
                "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1B\"; gene_name \"ONE_B\";\n",
                "chr17\tsrc\tgene\t1\t8\t.\t-\t.\tgene_id \"GENE17\"; gene_name \"SEVENTEEN\";\n"
            ),
        )
        .unwrap();

        let catalog_path_v2 = root.join("catalog_v2.json");
        let catalog_json_v2 = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            source_fasta_v2.display(),
            source_ann_v2.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path_v2, catalog_json_v2).unwrap();
        let catalog_v2 = GenomeCatalog::from_json_file(&catalog_path_v2.to_string_lossy()).unwrap();

        let mut phases: Vec<String> = vec![];
        let report = catalog_v2
            .reindex_genome_once_with_progress("ToyGenome", None, &mut |progress| {
                phases.push(progress.phase);
                true
            })
            .unwrap();
        assert!(report.reused_existing);
        assert!(phases.iter().any(|phase| phase == "reset_indexes"));
        assert!(phases.iter().any(|phase| phase == "reuse_sequence"));
        assert!(phases.iter().any(|phase| phase == "reuse_annotation"));
        assert!(!phases.iter().any(|phase| phase == "download_sequence"));
        assert!(!phases.iter().any(|phase| phase == "download_annotation"));
        assert!(
            report
                .sequence_source
                .as_deref()
                .is_some_and(|source| source.ends_with("source_v1.fa")),
            "report sequence source should preserve the cached-source provenance: {:?}",
            report.sequence_source
        );
        assert!(
            report
                .annotation_source
                .as_deref()
                .is_some_and(|source| source.ends_with("source_v1.gtf")),
            "report annotation source should preserve the cached-source provenance: {:?}",
            report.annotation_source
        );
        assert!(report.warnings.iter().any(|warning| {
            warning.contains("reindex used the cached local files without deleting or re-downloading them")
        }));

        assert_eq!(
            catalog_v2.get_sequence_region("ToyGenome", "chr1", 1, 4).unwrap(),
            "ACGT"
        );
        assert!(
            catalog_v2.get_sequence_region("ToyGenome", "chr17", 1, 4).is_err(),
            "reindex should keep using the cached local files when the catalog sources drift"
        );
        let genes = catalog_v2.list_gene_regions("ToyGenome", None).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].gene_id.as_deref(), Some("GENE1"));
    }

    #[test]
    fn test_prepare_genome_progress_callback_can_cancel() {
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
    "description": "toy test genome",
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
        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();

        let err = catalog
            .prepare_genome_once_with_progress("ToyGenome", None, &mut |_progress| false)
            .unwrap_err();
        assert!(is_prepare_cancelled_error(&err));
    }

    #[test]
    fn test_prepare_warns_but_succeeds_when_makeblastdb_missing() {
        let td = tempdir().unwrap();
        let root = td.path();

        let fasta_gz = root.join("toy.fa.gz");
        let ann_gz = root.join("toy.gtf.gz");
        write_gzip(&fasta_gz, ">chr1\nACGT\nACGT\n");
        write_gzip(
            &ann_gz,
            "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"MYGENE\";\n",
        );

        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
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
        let catalog = GenomeCatalog::from_json_file(&catalog_path.to_string_lossy()).unwrap();

        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        let report = catalog.prepare_genome_once("ToyGenome").unwrap();
        assert!(!report.blast_index_ready);
        assert!(!report.warnings.is_empty());
        assert!(
            report
                .warnings
                .iter()
                .any(|w| w.to_ascii_lowercase().contains("makeblastdb"))
        );

        let inspection = catalog
            .inspect_prepared_genome("ToyGenome", None)
            .unwrap()
            .expect("inspection should exist");
        assert!(!inspection.blast_index_ready);
    }

    #[test]
    fn test_parse_annotation_attributes_handles_quoted_semicolons() {
        let attrs = parse_annotation_attributes(
            r#"gene_id "GENE1"; gene_name "A;B"; Dbxref=GeneID:12345,HGNC:HGNC:5;"#,
        );
        assert_eq!(attrs.get("gene_id").map(String::as_str), Some("GENE1"));
        assert_eq!(attrs.get("gene_name").map(String::as_str), Some("A;B"));
        assert_eq!(
            attrs.get("dbxref").map(String::as_str),
            Some("GeneID:12345,HGNC:HGNC:5")
        );
    }

    #[test]
    fn test_parse_annotation_gene_records_tolerates_unexpected_lines() {
        let td = tempdir().unwrap();
        let path = td.path().join("mixed.gff");
        let mut bytes = vec![];
        bytes.extend_from_slice(b"##gff-version 3\n");
        bytes.extend_from_slice(
            b"chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene:GENE1;Name=Gene One;biotype=protein_coding\n",
        );
        bytes.extend_from_slice(
            b"chr1\tsrc\tprotein_coding_gene\t200\t300\t.\t-\t.\tgene_id \"GENE2\"; gene_name \"A;B\"; gene_biotype \"lncRNA\";\n",
        );
        bytes.extend_from_slice(b"chr1\tsrc\tgene\tbad\t400\t.\t+\t.\tID=gene:BAD\n");
        bytes.extend_from_slice(b"chr1\tsrc\tgene\t500\t400\t.\t+\t.\tID=gene:INV\n");
        bytes.extend_from_slice(
            b"chr1\tsrc\tpseudogene\t800\t900\t.\t+\t.\tDbxref=GeneID:12345,HGNC:HGNC:5\n",
        );
        bytes.extend_from_slice(b"chr1\tsrc\tgene\t<1001\t>1050\t.\t+\t.\tID=gene:GENE4\n");
        bytes.extend_from_slice(b"chr1\tsrc\tgene\t1200\t1300\t.\t+\t.\t.\n");
        bytes.extend_from_slice(&[0xff, 0xfe, 0x00, b'\n']);
        fs::write(&path, bytes).unwrap();

        let (genes, report) =
            parse_annotation_gene_records_with_progress(&path, |_done, _total| true).unwrap();
        assert_eq!(genes.len(), 5);
        assert_eq!(report.parsed_gene_records, 5);
        assert_eq!(report.malformed_lines, 3);
        assert!(!report.issue_examples.is_empty());
        assert!(genes.iter().any(|g| g.gene_id.as_deref() == Some("GENE1")
            && g.gene_name.as_deref() == Some("Gene One")
            && g.biotype.as_deref() == Some("protein_coding")));
        assert!(genes.iter().any(|g| g.gene_id.as_deref() == Some("GENE2")
            && g.gene_name.as_deref() == Some("A;B")
            && g.biotype.as_deref() == Some("lncRNA")
            && g.strand == Some('-')));
        assert!(
            genes
                .iter()
                .any(|g| g.gene_id.as_deref() == Some("12345") && g.gene_name.is_none())
        );
        assert!(genes.iter().any(|g| g.gene_id.as_deref() == Some("GENE4")
            && g.start_1based == 1001
            && g.end_1based == 1050));
        assert!(
            genes
                .iter()
                .any(|g| g.gene_id.as_deref() == Some("chr1_1200_1300"))
        );
    }

    #[test]
    fn test_parse_annotation_gene_records_handles_malformed_gtf_fields() {
        let td = tempdir().unwrap();
        let path = td.path().join("malformed.gtf");
        let mut bytes = vec![];
        bytes.extend_from_slice(b"chr1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id \"G1\"; gene_name \"Good One\"; gene_biotype \"protein_coding\";\n");
        bytes.extend_from_slice(b"chr1 src gene 11 20 . + . gene_id \"BROKEN_SPACES\"\n");
        bytes.extend_from_slice(b"chr1\tsrc\tgene\t21\tthirty\t.\t+\t.\tgene_id \"BROKEN_END\";\n");
        bytes.extend_from_slice(
            b"chr1\tsrc\tgene\t184467440737095516151\t190\t.\t+\t.\tgene_id \"TOO_BIG\";\n",
        );
        bytes.extend_from_slice(b"chr1\tsrc\tgene\t31\t40\t.\t?\t.\tgene_name \"No strand\";\n");
        bytes.extend_from_slice(
            b"chr1\tsrc\tgene\t41\t50\t.\t-\t.\tID=gene:G2;Name=Good Two;gene_biotype=lncRNA\n",
        );
        bytes.extend_from_slice(b"chr1\tsrc\texon\t41\t50\t.\t-\t.\tID=exon:E1\n");
        fs::write(&path, bytes).unwrap();

        let (genes, report) =
            parse_annotation_gene_records_with_progress(&path, |_done, _total| true).unwrap();
        assert_eq!(genes.len(), 3);
        assert_eq!(report.malformed_lines, 3);
        assert!(!report.issue_examples.is_empty());
        assert!(genes.iter().any(|g| {
            g.gene_id.as_deref() == Some("G1")
                && g.gene_name.as_deref() == Some("Good One")
                && g.biotype.as_deref() == Some("protein_coding")
                && g.strand == Some('+')
        }));
        assert!(genes.iter().any(|g| {
            g.gene_id.as_deref() == Some("G2")
                && g.gene_name.as_deref() == Some("Good Two")
                && g.biotype.as_deref() == Some("lncRNA")
                && g.strand == Some('-')
        }));
        assert!(genes.iter().any(|g| {
            g.gene_name.as_deref() == Some("No strand")
                && g.strand.is_none()
                && g.start_1based == 31
                && g.end_1based == 40
        }));
    }

    #[test]
    fn test_build_ncbi_assembly_ftp_urls() {
        let entry = GenomeCatalogEntry {
            ncbi_assembly_accession: Some("gcf_000005845.2".to_string()),
            ncbi_assembly_name: Some("ASM584v2".to_string()),
            ..Default::default()
        };
        let seq = resolve_ncbi_assembly_source("sequence", &entry)
            .unwrap()
            .expect("sequence URL");
        let ann = resolve_ncbi_assembly_source("annotation", &entry)
            .unwrap()
            .expect("annotation URL");
        assert_eq!(
            seq,
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
        );
        assert_eq!(
            ann,
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz"
        );
    }

    #[test]
    fn test_ncbi_assembly_source_requires_accession_and_name() {
        let accession_only = GenomeCatalogEntry {
            ncbi_assembly_accession: Some("GCF_000005845.2".to_string()),
            ..Default::default()
        };
        let err = resolve_ncbi_assembly_source("sequence", &accession_only).unwrap_err();
        assert!(err.contains("requires both"));
    }

    #[test]
    fn test_build_genbank_efetch_urls() {
        let _lock = genbank_env_lock().lock().unwrap();
        let _env = EnvVarGuard::set(NCBI_EFETCH_ENV_VAR, DEFAULT_NCBI_EFETCH_ENDPOINT);
        let entry = GenomeCatalogEntry {
            genbank_accession: Some("L09137".to_string()),
            ..Default::default()
        };
        let seq = resolve_genbank_accession_source("sequence", &entry)
            .unwrap()
            .expect("sequence URL");
        let ann = resolve_genbank_accession_source("annotation", &entry)
            .unwrap()
            .expect("annotation URL");
        assert_eq!(
            seq,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137&rettype=fasta&retmode=text"
        );
        assert_eq!(
            ann,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137&rettype=gbwithparts&retmode=text"
        );
    }

    #[test]
    fn test_source_plan_keeps_genbank_source_type_with_file_override() {
        let _lock = genbank_env_lock().lock().unwrap();
        let td = tempdir().unwrap();
        let root = td.path();
        let mock_dir = root.join("mock");
        fs::create_dir_all(&mock_dir).unwrap();
        let efetch_template = format!("file://{}/{{accession}}.{{rettype}}", mock_dir.display());
        let _env = EnvVarGuard::set(NCBI_EFETCH_ENV_VAR, &efetch_template);

        let catalog_path = root.join("catalog.json");
        let catalog = r#"{
  "Helper": {
    "genbank_accession": "L09137"
  }
}"#;
        fs::write(&catalog_path, catalog).unwrap();
        let catalog =
            GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref()).unwrap();
        let plan = catalog.source_plan("Helper", None).unwrap();
        assert_eq!(plan.sequence_source_type, "genbank_accession");
        assert_eq!(plan.annotation_source_type, "genbank_accession");
        assert!(plan.sequence_source.contains("L09137.fasta"));
        assert!(plan.annotation_source.contains("L09137.gbwithparts"));
    }

    #[test]
    fn test_genbank_accession_rejects_unpublished_placeholders() {
        let entry = GenomeCatalogEntry {
            genbank_accession: Some("LOCAL_UNPUBLISHED".to_string()),
            ..Default::default()
        };
        let err = resolve_genbank_accession_source("sequence", &entry).unwrap_err();
        assert!(err.contains("unpublished"));
    }

    #[test]
    fn test_catalog_validation_rejects_inconsistent_helper_fields() {
        let td = tempdir().unwrap();
        let catalog_path = td.path().join("bad_helper_catalog.json");
        let catalog = r#"{
  "BadHelper1": {
    "genbank_accession": "LOCAL_UNPUBLISHED",
    "sequence_local": "local.fa.gz"
  },
  "BadHelper2": {
    "genbank_accession": "L09137",
    "ncbi_assembly_accession": "GCF_000005845.2",
    "ncbi_assembly_name": "ASM584v2"
  }
}"#;
        fs::write(&catalog_path, catalog).unwrap();
        let err = GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref())
            .expect_err("catalog validation should fail");
        assert!(err.contains("BadHelper1"));
        assert!(err.contains("local_variant_unpublished"));
        assert!(err.contains("BadHelper2"));
        assert!(err.contains("cannot declare both"));
    }

    #[test]
    fn test_catalog_validation_accepts_unpublished_local_variant() {
        let td = tempdir().unwrap();
        let root = td.path();
        let seq = root.join("vector.fa");
        let ann = root.join("vector.gb");
        fs::write(&seq, ">v\nACGT\n").unwrap();
        fs::write(
            &ann,
            "LOCUS       VEC\nFEATURES             Location/Qualifiers\n",
        )
        .unwrap();
        let catalog_path = root.join("helper_catalog.json");
        let catalog = format!(
            r#"{{
  "Local Helper": {{
    "genbank_accession": "LOCAL_UNPUBLISHED",
    "local_variant_unpublished": true,
    "sequence_local": "{}",
    "annotations_local": "{}"
  }}
}}"#,
            seq.display(),
            ann.display()
        );
        fs::write(&catalog_path, catalog).unwrap();
        GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref())
            .expect("catalog validation should accept unpublished local variant");
    }

    #[test]
    fn test_catalog_validation_rejects_unpublished_flag_without_placeholder() {
        let td = tempdir().unwrap();
        let catalog_path = td.path().join("bad_helper_catalog.json");
        let catalog = r#"{
  "BadHelper": {
    "genbank_accession": "L09137",
    "local_variant_unpublished": true
  }
}"#;
        fs::write(&catalog_path, catalog).unwrap();
        let err = GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref())
            .expect_err("catalog validation should fail");
        assert!(err.contains("local_variant_unpublished"));
        assert!(err.contains("placeholder"));
    }

    #[test]
    fn test_local_missing_falls_back_to_genbank_accession_source() {
        let _lock = genbank_env_lock().lock().unwrap();
        let _env = EnvVarGuard::set(NCBI_EFETCH_ENV_VAR, DEFAULT_NCBI_EFETCH_ENDPOINT);
        let entry = GenomeCatalogEntry {
            sequence_local: Some("missing/local.fa.gz".to_string()),
            genbank_accession: Some("L09137".to_string()),
            ..Default::default()
        };
        let mut entries = HashMap::new();
        entries.insert("Vec".to_string(), entry.clone());
        let catalog = GenomeCatalog {
            entries,
            catalog_base_dir: PathBuf::from("."),
        };
        let source = catalog
            .resolve_source(
                "Vec",
                "sequence",
                entry.sequence_local.as_ref(),
                entry.sequence_remote.as_ref(),
                &entry,
            )
            .unwrap();
        assert!(source.contains("db=nuccore&id=L09137"));
        assert!(source.contains("rettype=fasta"));
    }

    #[test]
    fn test_entry_resolves_ncbi_assembly_alias_and_accession() {
        let mut entries = HashMap::new();
        entries.insert(
            "Human GRCh38 NCBI RefSeq GCF_000001405.40".to_string(),
            GenomeCatalogEntry {
                ncbi_assembly_accession: Some("GCF_000001405.40".to_string()),
                ncbi_assembly_name: Some("GRCh38.p14".to_string()),
                sequence_remote: Some("https://example.invalid/grch38.fa.gz".to_string()),
                annotations_remote: Some("https://example.invalid/grch38.gtf.gz".to_string()),
                ..Default::default()
            },
        );
        let catalog = GenomeCatalog {
            entries,
            catalog_base_dir: PathBuf::from("."),
        };

        assert!(catalog.source_plan("GRCh38.p14", None).is_ok());
        assert!(catalog.source_plan("gcf_000001405.40", None).is_ok());
    }

    #[test]
    fn test_entry_alias_rejects_ambiguous_matches() {
        let mut entries = HashMap::new();
        entries.insert(
            "GenomeA".to_string(),
            GenomeCatalogEntry {
                ncbi_assembly_name: Some("GRCh38.p14".to_string()),
                sequence_remote: Some("https://example.invalid/a.fa.gz".to_string()),
                annotations_remote: Some("https://example.invalid/a.gtf.gz".to_string()),
                ..Default::default()
            },
        );
        entries.insert(
            "GenomeB".to_string(),
            GenomeCatalogEntry {
                ncbi_assembly_name: Some("GRCh38.p14".to_string()),
                sequence_remote: Some("https://example.invalid/b.fa.gz".to_string()),
                annotations_remote: Some("https://example.invalid/b.gtf.gz".to_string()),
                ..Default::default()
            },
        );
        let catalog = GenomeCatalog {
            entries,
            catalog_base_dir: PathBuf::from("."),
        };

        let err = catalog.source_plan("GRCh38.p14", None).unwrap_err();
        assert!(err.contains("matches multiple catalog entries"));
    }

    #[test]
    fn test_resolve_prepared_genome_id_falls_back_to_unique_compatible_assembly() {
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
        let catalog = GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref())
            .expect("load catalog");
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog
            .prepare_genome_once("Human GRCh38 Ensembl 116")
            .expect("prepare compatible assembly");

        let resolved = catalog
            .resolve_prepared_genome_id("GRCh38.p14", None)
            .expect("resolve compatible prepared fallback");
        assert_eq!(resolved.resolved_genome_id, "Human GRCh38 Ensembl 116");
        assert!(
            resolved
                .fallback_warning
                .as_deref()
                .map(|msg| msg.contains("compatible prepared genome"))
                .unwrap_or(false)
        );

        let seq = catalog
            .get_sequence_region("GRCh38.p14", "chr1", 1, 4)
            .expect("extract sequence via fallback");
        assert_eq!(seq, "ACGT");
    }

    #[test]
    fn test_resolve_prepared_genome_id_lists_options_for_ambiguous_compatible_entries() {
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
  "Human GRCh38 Ensembl 113": {{
    "ncbi_taxonomy_id": 9606,
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }},
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
            cache_dir.display(),
            fasta.display(),
            ann.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog = GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref())
            .expect("load catalog");
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog
            .prepare_genome_once("Human GRCh38 Ensembl 113")
            .expect("prepare first compatible assembly");
        catalog
            .prepare_genome_once("Human GRCh38 Ensembl 116")
            .expect("prepare second compatible assembly");

        let err = catalog
            .resolve_prepared_genome_id("GRCh38.p14", None)
            .expect_err("ambiguous compatible options should require explicit selection");
        assert!(err.contains("Compatible prepared options"));
        assert!(err.contains("Human GRCh38 Ensembl 113"));
        assert!(err.contains("Human GRCh38 Ensembl 116"));
    }

    #[test]
    fn test_resolve_prepared_genome_id_policy_off_requires_explicit_selection() {
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
        let catalog = GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref())
            .expect("load catalog");
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog
            .prepare_genome_once("Human GRCh38 Ensembl 116")
            .expect("prepare compatible assembly");

        let err = catalog
            .resolve_prepared_genome_id_with_policy(
                "GRCh38.p14",
                None,
                PreparedGenomeFallbackPolicy::Off,
            )
            .expect_err("policy off should reject compatibility fallback");
        assert!(err.contains("policy is 'off'"));
        assert!(err.contains("Human GRCh38 Ensembl 116"));
    }

    #[test]
    fn test_resolve_prepared_genome_id_policy_always_explicit_requires_selection_even_when_unique()
    {
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
        let catalog = GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref())
            .expect("load catalog");
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog
            .prepare_genome_once("Human GRCh38 Ensembl 116")
            .expect("prepare compatible assembly");

        let err = catalog
            .resolve_prepared_genome_id_with_policy(
                "GRCh38.p14",
                None,
                PreparedGenomeFallbackPolicy::AlwaysExplicit,
            )
            .expect_err("always_explicit should require explicit selection");
        assert!(err.contains("requires explicit selection"));
        assert!(err.contains("Human GRCh38 Ensembl 116"));
    }

    #[test]
    fn test_inspect_prepared_genome_compatibility_reports_exact_and_options() {
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
        let catalog = GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref())
            .expect("load catalog");
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog
            .prepare_genome_once("Human GRCh38 Ensembl 116")
            .expect("prepare compatible assembly");

        let inspection = catalog
            .inspect_prepared_genome_compatibility("GRCh38.p14", None)
            .expect("inspect compatibility");
        assert!(!inspection.exact_prepared);
        assert_eq!(
            inspection.requested_catalog_key,
            "Human GRCh38 NCBI RefSeq GCF_000001405.40"
        );
        assert_eq!(
            inspection.compatible_prepared_options,
            vec!["Human GRCh38 Ensembl 116".to_string()]
        );
    }

    #[test]
    fn test_parse_genbank_annotation_records_pgex() {
        let (genes, report) = parse_annotation_gene_records_with_progress(
            Path::new("test_files/pGEX-3X.gb"),
            |_done, _total| true,
        )
        .unwrap();
        assert!(!genes.is_empty());
        assert_eq!(report.source_format, "genbank");
        assert!(
            genes
                .iter()
                .any(|g| g.biotype.as_deref() == Some("promoter"))
        );
        assert!(genes.iter().any(|g| {
            g.biotype.as_deref() == Some("origin_of_replication")
                && g.start_1based <= 2978
                && g.end_1based >= 2978
        }));
        assert!(genes.iter().any(|g| {
            let is_bla = g
                .gene_name
                .as_ref()
                .map(|name| name.eq_ignore_ascii_case("bla"))
                .unwrap_or(false)
                || g.gene_id
                    .as_ref()
                    .map(|id| id.eq_ignore_ascii_case("bla"))
                    .unwrap_or(false);
            is_bla
                && (1289..=1291).contains(&g.start_1based)
                && (2219..=2221).contains(&g.end_1based)
        }));
    }

    #[test]
    fn test_parse_xml_annotation_records_toy_small() {
        let (xml_genes, report) = parse_annotation_gene_records_with_progress(
            Path::new("test_files/fixtures/import_parity/toy.small.gbseq.xml"),
            |_done, _total| true,
        )
        .expect("parse XML GBSeq annotations");
        let (genbank_genes, _) = parse_annotation_gene_records_with_progress(
            Path::new("test_files/fixtures/import_parity/toy.small.gb"),
            |_done, _total| true,
        )
        .expect("parse GenBank toy annotations");
        assert_eq!(report.source_format, "xml_gbseq");
        for gene_name in ["toyA", "toyB"] {
            let from_xml = xml_genes
                .iter()
                .find(|g| g.gene_name.as_deref() == Some(gene_name))
                .expect("gene should exist in XML parse");
            let from_genbank = genbank_genes
                .iter()
                .find(|g| g.gene_name.as_deref() == Some(gene_name))
                .expect("gene should exist in GenBank parse");
            assert_eq!(from_xml.start_1based, from_genbank.start_1based);
            assert_eq!(from_xml.end_1based, from_genbank.end_1based);
            assert_eq!(from_xml.strand, from_genbank.strand);
        }
    }

    #[test]
    fn test_parse_xml_annotation_rejects_insd_dialect() {
        let td = tempdir().unwrap();
        let xml_path = td.path().join("insd.annotation.xml");
        fs::write(&xml_path, "<INSDSet><INSDSeq/></INSDSet>").expect("write INSD XML");
        let err =
            parse_annotation_gene_records_with_progress(xml_path.as_path(), |_done, _total| true)
                .expect_err("INSDSet should be rejected");
        assert!(
            err.contains("Unsupported XML dialect"),
            "expected unsupported-dialect error, got: {err}"
        );
    }

    #[test]
    fn test_normalize_genbank_biotypes_for_helper_vectors() {
        let mut attrs = HashMap::new();
        attrs.insert("label".to_string(), "5' LTR".to_string());
        assert_eq!(
            normalize_genbank_feature_biotype("misc_feature", &attrs).as_deref(),
            Some("ltr")
        );

        attrs.insert("label".to_string(), "AAV ITR left".to_string());
        assert_eq!(
            normalize_genbank_feature_biotype("misc_feature", &attrs).as_deref(),
            Some("itr")
        );

        attrs.insert("label".to_string(), "pUC ori".to_string());
        assert_eq!(
            normalize_genbank_feature_biotype("misc_feature", &attrs).as_deref(),
            Some("origin_of_replication")
        );

        attrs.insert("label".to_string(), "CMV promoter".to_string());
        assert_eq!(
            normalize_genbank_feature_biotype("misc_feature", &attrs).as_deref(),
            Some("promoter")
        );
    }

    #[test]
    fn test_source_plan_reports_local_ncbi_and_genbank_types() {
        let td = tempdir().unwrap();
        let root = td.path();
        let local_seq = root.join("local.fa");
        let local_ann = root.join("local.gff3");
        fs::write(&local_seq, ">chr1\nACGT\n").unwrap();
        fs::write(&local_ann, "##gff-version 3\n").unwrap();
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "LocalGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}"
  }},
  "AssemblyGenome": {{
    "ncbi_assembly_accession": "GCF_000005845.2",
    "ncbi_assembly_name": "ASM584v2"
  }},
  "GenbankGenome": {{
    "genbank_accession": "L09137"
  }}
}}"#,
            local_seq.display(),
            local_ann.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog =
            GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref()).unwrap();
        let local = catalog.source_plan("LocalGenome", None).unwrap();
        assert_eq!(local.sequence_source_type, "local");
        assert_eq!(local.annotation_source_type, "local");
        let assembly = catalog.source_plan("AssemblyGenome", None).unwrap();
        assert_eq!(assembly.sequence_source_type, "ncbi_assembly");
        assert_eq!(assembly.annotation_source_type, "ncbi_assembly");
        let gb = catalog.source_plan("GenbankGenome", None).unwrap();
        assert_eq!(gb.sequence_source_type, "genbank_accession");
        assert_eq!(gb.annotation_source_type, "genbank_accession");
    }

    #[test]
    fn test_source_plan_uses_catalog_molecular_mass_when_provided() {
        let td = tempdir().unwrap();
        let root = td.path();
        let local_seq = root.join("helper.fa");
        let local_ann = root.join("helper.gff3");
        fs::write(&local_seq, ">chrH\nACGT\n").unwrap();
        fs::write(&local_ann, "##gff-version 3\n").unwrap();
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "HelperMassKnown": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "nucleotide_length_bp": 5000,
    "molecular_mass_da": 3089836.04
  }}
}}"#,
            local_seq.display(),
            local_ann.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog =
            GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref()).unwrap();
        let plan = catalog.source_plan("HelperMassKnown", None).unwrap();
        assert_eq!(plan.nucleotide_length_bp, Some(5000));
        assert_eq!(plan.molecular_mass_da, Some(3_089_836.04));
        assert_eq!(plan.molecular_mass_source.as_deref(), Some("catalog"));
    }

    #[test]
    fn test_source_plan_estimates_mass_from_catalog_nucleotide_length() {
        let td = tempdir().unwrap();
        let root = td.path();
        let local_seq = root.join("helper.fa");
        let local_ann = root.join("helper.gff3");
        fs::write(&local_seq, ">chrH\nACGT\n").unwrap();
        fs::write(&local_ann, "##gff-version 3\n").unwrap();
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "HelperMassEstimated": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "nucleotide_length_bp": 2686
  }}
}}"#,
            local_seq.display(),
            local_ann.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog =
            GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref()).unwrap();
        let plan = catalog.source_plan("HelperMassEstimated", None).unwrap();
        assert_eq!(plan.nucleotide_length_bp, Some(2686));
        let expected = estimate_double_stranded_dna_mass_da(2686).unwrap();
        let observed = plan.molecular_mass_da.unwrap();
        assert!((observed - expected).abs() < 1e-6);
        assert_eq!(
            plan.molecular_mass_source.as_deref(),
            Some("estimated_from_nucleotide_length")
        );
    }

    #[test]
    fn test_source_plan_estimates_mass_from_prepared_sequence_length() {
        let td = tempdir().unwrap();
        let root = td.path();
        let local_seq = root.join("toy.fa");
        let local_ann = root.join("toy.gtf");
        fs::write(&local_seq, ">chr1\nACGTACGTACGT\n").unwrap();
        fs::write(
            &local_ann,
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tgene_id \"G1\"; gene_name \"G1\";\n",
        )
        .unwrap();
        let cache_dir = root.join("cache");
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "PreparedHelper": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            local_seq.display(),
            local_ann.display(),
            cache_dir.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let catalog =
            GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref()).unwrap();
        let _guard = EnvVarGuard::set(
            MAKEBLASTDB_ENV_BIN,
            "__gentle_makeblastdb_missing_for_test__",
        );
        catalog
            .prepare_genome_once("PreparedHelper")
            .expect("prepare should succeed");
        let plan = catalog.source_plan("PreparedHelper", None).unwrap();
        assert_eq!(plan.nucleotide_length_bp, Some(12));
        let expected = estimate_double_stranded_dna_mass_da(12).unwrap();
        let observed = plan.molecular_mass_da.unwrap();
        assert!((observed - expected).abs() < 1e-6);
        assert_eq!(
            plan.molecular_mass_source.as_deref(),
            Some("estimated_from_nucleotide_length")
        );
    }

    #[test]
    fn test_catalog_validation_rejects_invalid_length_and_mass_fields() {
        let td = tempdir().unwrap();
        let root = td.path();
        let local_seq = root.join("local.fa");
        let local_ann = root.join("local.gff3");
        fs::write(&local_seq, ">chr1\nACGT\n").unwrap();
        fs::write(&local_ann, "##gff-version 3\n").unwrap();
        let catalog_path = root.join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "BadLength": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "nucleotide_length_bp": 0
  }},
  "BadMass": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "molecular_mass_da": -3.0
  }}
}}"#,
            local_seq.display(),
            local_ann.display(),
            local_seq.display(),
            local_ann.display()
        );
        fs::write(&catalog_path, catalog_json).unwrap();
        let err = GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref())
            .expect_err("catalog validation should fail");
        assert!(err.contains("BadLength"));
        assert!(err.contains("'nucleotide_length_bp' must be a positive integer"));
        assert!(err.contains("BadMass"));
        assert!(err.contains("'molecular_mass_da' must be a finite positive number"));
    }

    #[test]
    fn test_repo_assets_genome_catalog_is_valid() {
        let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let catalog_path = root.join("assets").join("genomes.json");
        let catalog =
            GenomeCatalog::from_json_file(catalog_path.to_string_lossy().as_ref()).unwrap();
        let genomes = catalog.list_genomes();
        for expected in [
            "Human GRCh38 Ensembl 113",
            "Human GRCh38 Ensembl 116",
            "Human GRCh38 NCBI RefSeq GCF_000001405.40",
            "Mouse GRCm39 Ensembl 116",
            "Rat GRCr8 Ensembl 116",
            "Caenorhabditis elegans WBcel235 Ensembl 115",
            "Saccharomyces cerevisiae S288c Ensembl 113",
            "Saccharomyces cerevisiae S288c Ensembl 115",
            "LocalProject",
        ] {
            assert!(
                genomes.iter().any(|id| id == expected),
                "Expected catalog entry '{expected}' missing in assets/genomes.json"
            );
        }

        let local = catalog.source_plan("LocalProject", None).unwrap();
        assert_eq!(local.sequence_source_type, "local");
        assert_eq!(local.annotation_source_type, "local");
    }

    #[test]
    fn test_normalize_blast_task_accepts_and_rejects_expected_values() {
        assert_eq!(
            normalize_blast_task(None).expect("default task"),
            "blastn-short"
        );
        assert_eq!(
            normalize_blast_task(Some(" blastn ")).expect("blastn task"),
            "blastn"
        );
        assert!(normalize_blast_task(Some("megablast")).is_err());
    }

    #[test]
    fn test_is_blast_cancelled_error_matches_prefix() {
        let msg = blast_cancelled_error("test");
        assert!(is_blast_cancelled_error(&msg));
        assert!(!is_blast_cancelled_error("regular blast error"));
    }

    #[test]
    fn test_normalize_blast_query_sequence_validates_and_normalizes() {
        assert_eq!(
            normalize_blast_query_sequence("acgu n").expect("query should normalize"),
            "ACGTN"
        );
        assert!(normalize_blast_query_sequence("").is_err());
        assert!(normalize_blast_query_sequence("ACGT-").is_err());
    }

    #[test]
    fn test_parse_blastn_tabular_hits_reads_expected_columns() {
        let stdout = "query\tsubject1\t99.1\t20\t0\t0\t1\t20\t100\t119\t1e-10\t50.2\t100\n";
        let (hits, warnings) = parse_blastn_tabular_hits(stdout);
        assert!(warnings.is_empty());
        assert_eq!(hits.len(), 1);
        let hit = &hits[0];
        assert_eq!(hit.subject_id, "subject1");
        assert_eq!(hit.identity_percent, 99.1);
        assert_eq!(hit.alignment_length, 20);
        assert_eq!(hit.mismatches, 0);
        assert_eq!(hit.gap_opens, 0);
        assert_eq!(hit.query_start, 1);
        assert_eq!(hit.query_end, 20);
        assert_eq!(hit.subject_start, 100);
        assert_eq!(hit.subject_end, 119);
        assert_eq!(hit.evalue, 1e-10);
        assert_eq!(hit.bit_score, 50.2);
        assert_eq!(hit.query_coverage_percent, Some(100.0));
    }
}
