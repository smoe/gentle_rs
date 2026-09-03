//! Optional read-only queries over finalized `jaspar-mapping` sparse packages.
//!
//! GENtle intentionally shells out to DuckDB only for an explicit query. The
//! package and executable are not startup dependencies, and every discovery or
//! execution failure is returned as typed unavailable evidence rather than
//! disabling local motif scoring.

use crate::digest_utils::{sha256_prefixed_bytes, short_sha256_id};
use gentle_protocol::{
    GENOMIC_MOTIF_EVIDENCE_SCHEMA, GenomicMotifEvidenceAvailability,
    GenomicMotifEvidenceCompatibilityStatus, GenomicMotifEvidenceCoverageStatus,
    GenomicMotifEvidenceHit, GenomicMotifEvidenceMotifCoverage,
    GenomicMotifEvidencePayloadProvenance, GenomicMotifEvidenceProviderProvenance,
    GenomicMotifEvidenceReport, GenomicMotifEvidenceRequest, GenomicMotifEvidenceResolvedRegion,
    MAX_GENOMIC_MOTIF_EVIDENCE_QUERY_MOTIFS,
};
use serde_json::Value;
use std::{
    collections::{BTreeMap, BTreeSet},
    env, fs,
    io::{BufReader, Read},
    path::{Path, PathBuf},
    process::{Command, Output, Stdio},
    thread,
    time::{Duration, Instant},
};

pub const JASPAR_GENOME_SCAN_PACKAGE_ENV: &str = "GENTLE_JASPAR_GENOME_SCAN_PACKAGE";
pub const DUCKDB_BIN_ENV: &str = "GENTLE_DUCKDB_BIN";

const MAX_REGIONS: usize = 256;
const MAX_ROWS_HARD: usize = 1_000_000;
const MAX_PAYLOAD_FILES_HARD: usize = 4_096;
const MAX_TIMEOUT_SECONDS: u64 = 300;
const MAX_MANIFEST_BYTES: u64 = 4 * 1024 * 1024;
const MAX_DUCKDB_OUTPUT_BYTES: usize = 64 * 1024 * 1024;

#[derive(Debug, Clone)]
pub(crate) struct GenomicMotifQueryRegion {
    pub interval_id: String,
    pub label: Option<String>,
    pub chromosome: String,
    pub start_0based: u64,
    pub end_0based_exclusive: u64,
    pub source_seq_id: Option<String>,
    pub source_start_0based: Option<usize>,
    pub source_end_0based_exclusive: Option<usize>,
    pub source_sequence_length_bp: Option<usize>,
    pub source_anchor_start_1based: Option<usize>,
    pub source_anchor_end_1based: Option<usize>,
    pub source_anchor_reverse: bool,
}

#[derive(Debug)]
struct BoundedOutput {
    output: Output,
    stdout_truncated: bool,
    stderr_truncated: bool,
}

#[derive(Debug, Clone)]
struct PackagePaths {
    root: PathBuf,
    manifest_path: PathBuf,
    database_path: PathBuf,
    manifest: Value,
    manifest_bytes: Vec<u8>,
}

#[derive(Debug, Clone)]
struct InventoryRow {
    task_id: String,
    output_relative_path: String,
    sha256: String,
    emitted_hits: u64,
}

#[derive(Debug, Clone)]
struct MotifThresholdRow {
    motif_id: String,
    motif_name: Option<String>,
    threshold_set_id: Option<String>,
    informative_threshold: Option<f64>,
    final_minimum_score: Option<f64>,
    density_limited: Option<bool>,
}

fn declared_content_fingerprint(paths: &PackagePaths, inventory_rows: &[InventoryRow]) -> String {
    let mut identity = sha256_prefixed_bytes(&paths.manifest_bytes);
    for row in inventory_rows {
        identity.push('\n');
        identity.push_str(&row.task_id);
        identity.push('\t');
        identity.push_str(&row.output_relative_path);
        identity.push('\t');
        identity.push_str(&row.sha256);
        identity.push('\t');
        identity.push_str(&row.emitted_hits.to_string());
    }
    sha256_prefixed_bytes(identity.as_bytes())
}

fn report_id(request: &GenomicMotifEvidenceRequest, manifest_sha256: Option<&str>) -> String {
    let request_json = serde_json::to_string(request).unwrap_or_default();
    let identity = match manifest_sha256 {
        Some(manifest_sha256) => format!("{request_json}\n{manifest_sha256}"),
        None => request_json,
    };
    short_sha256_id("genomic_motif_evidence", &identity)
}

fn unavailable_report(
    request: &GenomicMotifEvidenceRequest,
    regions: &[GenomicMotifQueryRegion],
    availability: GenomicMotifEvidenceAvailability,
    warning: impl Into<String>,
) -> GenomicMotifEvidenceReport {
    GenomicMotifEvidenceReport {
        schema: GENOMIC_MOTIF_EVIDENCE_SCHEMA.to_string(),
        report_id: report_id(request, None),
        availability,
        request: request.clone(),
        regions: regions
            .iter()
            .map(|region| GenomicMotifEvidenceResolvedRegion {
                interval_id: region.interval_id.clone(),
                label: region.label.clone(),
                requested_chromosome: region.chromosome.clone(),
                resolved_chromosome: None,
                start_0based: region.start_0based,
                end_0based_exclusive: region.end_0based_exclusive,
                source_seq_id: region.source_seq_id.clone(),
                source_start_0based: region.source_start_0based,
                source_end_0based_exclusive: region.source_end_0based_exclusive,
                compatibility_status: GenomicMotifEvidenceCompatibilityStatus::NotAssessed,
                package_contig_length_bp: None,
            })
            .collect(),
        warnings: vec![warning.into()],
        ..GenomicMotifEvidenceReport::default()
    }
}

pub(crate) fn validate_genomic_motif_request(
    request: &GenomicMotifEvidenceRequest,
    regions: &[GenomicMotifQueryRegion],
) -> Result<(), String> {
    if request.motif_ids.is_empty() {
        return Err("genomic motif evidence requires at least one motif_id".to_string());
    }
    if request.motif_ids.len() > MAX_GENOMIC_MOTIF_EVIDENCE_QUERY_MOTIFS {
        return Err(format!(
            "genomic motif evidence accepts at most {MAX_GENOMIC_MOTIF_EVIDENCE_QUERY_MOTIFS} motif IDs per row query"
        ));
    }
    let mut unique_motifs = BTreeSet::new();
    for raw_motif in &request.motif_ids {
        let motif = raw_motif.trim();
        if motif.is_empty() {
            return Err("genomic motif evidence motif IDs must not be blank".to_string());
        }
        if motif != raw_motif {
            return Err(format!(
                "genomic motif evidence motif ID '{raw_motif}' must not contain surrounding whitespace"
            ));
        }
        if motif.eq_ignore_ascii_case("all") {
            return Err(
                "genomic motif evidence row queries require an explicit bounded motif set; ALL is not accepted"
                    .to_string(),
            );
        }
        if !motif
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"._:-".contains(&byte))
        {
            return Err(format!(
                "genomic motif evidence motif ID '{motif}' contains unsupported characters"
            ));
        }
        if !unique_motifs.insert(motif) {
            return Err(format!(
                "genomic motif evidence motif ID '{motif}' is duplicated"
            ));
        }
    }
    if regions.is_empty() {
        return Err("genomic motif evidence requires at least one genomic interval".to_string());
    }
    if regions.len() > MAX_REGIONS {
        return Err(format!(
            "genomic motif evidence accepts at most {MAX_REGIONS} intervals per query"
        ));
    }
    let mut ids = BTreeSet::new();
    for region in regions {
        if region.interval_id.trim().is_empty() || !ids.insert(region.interval_id.clone()) {
            return Err(
                "genomic motif evidence interval IDs must be non-empty and unique".to_string(),
            );
        }
        if region.chromosome.trim().is_empty()
            || !region
                .chromosome
                .bytes()
                .all(|byte| byte.is_ascii_alphanumeric() || b"._:-".contains(&byte))
        {
            return Err(format!(
                "genomic motif evidence chromosome '{}' is blank or contains unsupported characters",
                region.chromosome
            ));
        }
        if region.end_0based_exclusive <= region.start_0based {
            return Err(format!(
                "genomic motif evidence interval '{}' must have start < end",
                region.interval_id
            ));
        }
    }
    if request.max_rows == 0 || request.max_rows > MAX_ROWS_HARD {
        return Err(format!(
            "genomic motif evidence max_rows must be within 1..={MAX_ROWS_HARD}"
        ));
    }
    if request.max_payload_files == 0 || request.max_payload_files > MAX_PAYLOAD_FILES_HARD {
        return Err(format!(
            "genomic motif evidence max_payload_files must be within 1..={MAX_PAYLOAD_FILES_HARD}"
        ));
    }
    if request.timeout_seconds == 0 || request.timeout_seconds > MAX_TIMEOUT_SECONDS {
        return Err(format!(
            "genomic motif evidence timeout_seconds must be within 1..={MAX_TIMEOUT_SECONDS}"
        ));
    }
    for (name, value) in [
        ("minimum_score", request.minimum_score),
        (
            "minimum_pwm_relative_score",
            request.minimum_pwm_relative_score,
        ),
    ] {
        if value.is_some_and(|number| !number.is_finite()) {
            return Err(format!("genomic motif evidence {name} must be finite"));
        }
    }
    Ok(())
}

fn sql_string(value: &str) -> String {
    format!("'{}'", value.replace('\'', "''"))
}

fn sql_string_list<'a>(values: impl IntoIterator<Item = &'a str>) -> String {
    format!(
        "[{}]",
        values
            .into_iter()
            .map(sql_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn read_capped<R: Read>(reader: R, cap: usize) -> (Vec<u8>, bool) {
    let mut reader = BufReader::new(reader);
    let mut kept = Vec::with_capacity(cap.min(64 * 1024));
    let mut buffer = [0u8; 16 * 1024];
    let mut truncated = false;
    loop {
        match reader.read(&mut buffer) {
            Ok(0) | Err(_) => break,
            Ok(count) => {
                let remaining = cap.saturating_sub(kept.len());
                let take = remaining.min(count);
                kept.extend_from_slice(&buffer[..take]);
                truncated |= take < count;
            }
        }
    }
    (kept, truncated)
}

fn bounded_command_output(
    command: &mut Command,
    timeout: Duration,
    max_output_bytes: usize,
) -> Result<BoundedOutput, String> {
    command.stdout(Stdio::piped()).stderr(Stdio::piped());
    let mut child = command
        .spawn()
        .map_err(|error| format!("could not start DuckDB: {error}"))?;
    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| "could not capture DuckDB stdout".to_string())?;
    let stderr = child
        .stderr
        .take()
        .ok_or_else(|| "could not capture DuckDB stderr".to_string())?;
    let stdout_reader = thread::spawn(move || read_capped(stdout, max_output_bytes));
    let stderr_reader = thread::spawn(move || read_capped(stderr, max_output_bytes));
    let started = Instant::now();
    loop {
        match child.try_wait() {
            Ok(Some(status)) => {
                let (stdout, stdout_truncated) = stdout_reader.join().unwrap_or_default();
                let (stderr, stderr_truncated) = stderr_reader.join().unwrap_or_default();
                return Ok(BoundedOutput {
                    output: Output {
                        status,
                        stdout,
                        stderr,
                    },
                    stdout_truncated,
                    stderr_truncated,
                });
            }
            Ok(None) if started.elapsed() < timeout => thread::sleep(Duration::from_millis(25)),
            Ok(None) => {
                let _ = child.kill();
                let _ = child.wait();
                let _ = stdout_reader.join();
                let _ = stderr_reader.join();
                return Err(format!(
                    "DuckDB query exceeded the {} second timeout",
                    timeout.as_secs()
                ));
            }
            Err(error) => {
                let _ = child.kill();
                let _ = child.wait();
                let _ = stdout_reader.join();
                let _ = stderr_reader.join();
                return Err(format!("could not poll DuckDB: {error}"));
            }
        }
    }
}

fn run_duckdb(
    executable: &str,
    database: Option<&Path>,
    query: Option<&str>,
    timeout: Duration,
) -> Result<BoundedOutput, String> {
    let mut command = Command::new(executable);
    if let Some(database) = database {
        command
            .arg("-readonly")
            .arg("-json")
            .arg(database)
            .arg("-c")
            .arg(query.unwrap_or_default());
    } else {
        command.arg("--version");
    }
    bounded_command_output(&mut command, timeout, MAX_DUCKDB_OUTPUT_BYTES)
}

fn query_json(
    executable: &str,
    database: &Path,
    query: &str,
    timeout: Duration,
) -> Result<Vec<Value>, String> {
    let bounded = run_duckdb(executable, Some(database), Some(query), timeout)?;
    if bounded.stdout_truncated || bounded.stderr_truncated {
        return Err(format!(
            "DuckDB output exceeded the {} MiB safety cap",
            MAX_DUCKDB_OUTPUT_BYTES / 1024 / 1024
        ));
    }
    if !bounded.output.status.success() {
        let detail = String::from_utf8_lossy(&bounded.output.stderr)
            .trim()
            .to_string();
        return Err(format!(
            "DuckDB query failed with status {:?}: {}",
            bounded.output.status.code(),
            if detail.is_empty() {
                "no diagnostic output"
            } else {
                &detail
            }
        ));
    }
    serde_json::from_slice::<Vec<Value>>(&bounded.output.stdout)
        .map_err(|error| format!("DuckDB returned invalid JSON: {error}"))
}

fn path_is_within(root: &Path, candidate: &Path) -> bool {
    candidate.starts_with(root)
}

fn resolve_package_paths(request: &GenomicMotifEvidenceRequest) -> Result<PackagePaths, String> {
    let raw_root = request
        .package_root
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .or_else(|| {
            env::var(JASPAR_GENOME_SCAN_PACKAGE_ENV)
                .ok()
                .filter(|value| !value.trim().is_empty())
        })
        .ok_or_else(|| "package_not_configured".to_string())?;
    let root = PathBuf::from(raw_root);
    if !root.is_dir() {
        return Err(format!("package_missing:{}", root.display()));
    }
    let root = root
        .canonicalize()
        .map_err(|error| format!("invalid_package:could not resolve package root: {error}"))?;
    let manifest_path = root.join("manifest.json");
    let metadata = fs::metadata(&manifest_path).map_err(|error| {
        format!(
            "invalid_package:could not inspect '{}': {error}",
            manifest_path.display()
        )
    })?;
    if metadata.len() > MAX_MANIFEST_BYTES {
        return Err("invalid_package:manifest exceeds the 4 MiB safety cap".to_string());
    }
    let manifest_bytes = fs::read(&manifest_path).map_err(|error| {
        format!(
            "invalid_package:could not read '{}': {error}",
            manifest_path.display()
        )
    })?;
    let manifest: Value = serde_json::from_slice(&manifest_bytes)
        .map_err(|error| format!("invalid_package:manifest JSON is invalid: {error}"))?;
    if manifest.get("schema_version").and_then(Value::as_u64) != Some(2) {
        return Err(
            "invalid_package:finalized genome-scan manifest schema_version must be 2".to_string(),
        );
    }
    if manifest.get("state").and_then(Value::as_str) != Some("complete") {
        return Err("invalid_package:genome-scan manifest must declare state=complete".to_string());
    }
    let database_value = request
        .database_path
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .or_else(|| manifest.get("database").and_then(Value::as_str))
        .ok_or_else(|| "invalid_package:manifest does not name a DuckDB catalog".to_string())?;
    let database_candidate = PathBuf::from(database_value);
    let database_candidate = if database_candidate.is_absolute() {
        database_candidate
    } else {
        root.join(database_candidate)
    };
    let database_path = database_candidate.canonicalize().map_err(|error| {
        format!(
            "invalid_package:could not resolve DuckDB catalog '{}': {error}",
            database_candidate.display()
        )
    })?;
    if !path_is_within(&root, &database_path) {
        return Err("invalid_package:DuckDB catalog resolves outside the package root".to_string());
    }
    if !database_path.is_file() {
        return Err(format!(
            "invalid_package:DuckDB catalog is not a file: {}",
            database_path.display()
        ));
    }
    Ok(PackagePaths {
        root,
        manifest_path,
        database_path,
        manifest,
        manifest_bytes,
    })
}

fn optional_string(row: &Value, key: &str) -> Option<String> {
    row.get(key)
        .and_then(Value::as_str)
        .map(str::to_string)
        .filter(|value| !value.trim().is_empty())
}

fn required_string(row: &Value, key: &str) -> Result<String, String> {
    optional_string(row, key).ok_or_else(|| format!("DuckDB row is missing string field '{key}'"))
}

fn optional_u64(row: &Value, key: &str) -> Option<u64> {
    row.get(key)
        .and_then(|value| value.as_u64().or_else(|| value.as_i64()?.try_into().ok()))
}

fn optional_f64(row: &Value, key: &str) -> Option<f64> {
    row.get(key).and_then(Value::as_f64)
}

fn optional_bool(row: &Value, key: &str) -> Option<bool> {
    row.get(key).and_then(Value::as_bool)
}

fn motif_coverage_rows(
    request: &GenomicMotifEvidenceRequest,
    thresholds: &BTreeMap<String, MotifThresholdRow>,
    source_minimum_pwm_relative_score: Option<f64>,
    hits: &[GenomicMotifEvidenceHit],
    truncated: bool,
) -> Vec<GenomicMotifEvidenceMotifCoverage> {
    let hit_counts = hits.iter().fold(BTreeMap::new(), |mut counts, hit| {
        *counts.entry(hit.motif_id.as_str()).or_insert(0usize) += 1;
        counts
    });
    request
        .motif_ids
        .iter()
        .map(|motif_id| {
            let threshold = thresholds.get(motif_id);
            let status = match threshold {
                None => GenomicMotifEvidenceCoverageStatus::MotifNotInPackage,
                Some(_) if truncated => GenomicMotifEvidenceCoverageStatus::TruncatedAtMaxRows,
                Some(threshold)
                    if request
                        .minimum_score
                        .zip(threshold.final_minimum_score)
                        .is_some_and(|(requested, retained)| requested < retained)
                        || request
                            .minimum_pwm_relative_score
                            .zip(source_minimum_pwm_relative_score)
                            .is_some_and(|(requested, retained)| requested < retained) =>
                {
                    GenomicMotifEvidenceCoverageStatus::IncompleteBelowStorageFloor
                }
                Some(_)
                    if request.minimum_score.is_some()
                        || request.minimum_pwm_relative_score.is_some() =>
                {
                    GenomicMotifEvidenceCoverageStatus::CompleteForRequestedThreshold
                }
                Some(threshold) if threshold.density_limited == Some(true) => {
                    GenomicMotifEvidenceCoverageStatus::DensityLimitedAtSource
                }
                Some(_) => GenomicMotifEvidenceCoverageStatus::CompleteForPackageRetention,
            };
            GenomicMotifEvidenceMotifCoverage {
                motif_id: motif_id.clone(),
                motif_name: threshold.and_then(|row| row.motif_name.clone()),
                threshold_set_id: threshold.and_then(|row| row.threshold_set_id.clone()),
                informative_threshold: threshold.and_then(|row| row.informative_threshold),
                source_minimum_score: threshold.and_then(|row| row.final_minimum_score),
                source_minimum_pwm_relative_score,
                requested_minimum_score: request.minimum_score,
                requested_minimum_pwm_relative_score: request.minimum_pwm_relative_score,
                density_limited: threshold.and_then(|row| row.density_limited),
                status,
                returned_hit_count: hit_counts.get(motif_id.as_str()).copied().unwrap_or(0),
            }
        })
        .collect()
}

fn motif_coverage_warnings(coverage_rows: &[GenomicMotifEvidenceMotifCoverage]) -> Vec<String> {
    coverage_rows
        .iter()
        .filter_map(|coverage| match coverage.status {
            GenomicMotifEvidenceCoverageStatus::IncompleteBelowStorageFloor => Some(format!(
                "Motif '{}' was requested below the package retention floor; missing lower-scoring rows cannot be recovered from this sparse package",
                coverage.motif_id
            )),
            GenomicMotifEvidenceCoverageStatus::DensityLimitedAtSource => Some(format!(
                "Motif '{}' uses a density-limited package retention floor",
                coverage.motif_id
            )),
            GenomicMotifEvidenceCoverageStatus::MotifNotInPackage => Some(format!(
                "Motif '{}' is not represented in the package threshold registry",
                coverage.motif_id
            )),
            _ => None,
        })
        .collect()
}

fn motif_coverage_is_query_complete(coverage_rows: &[GenomicMotifEvidenceMotifCoverage]) -> bool {
    coverage_rows.iter().all(|coverage| {
        matches!(
            coverage.status,
            GenomicMotifEvidenceCoverageStatus::CompleteForPackageRetention
                | GenomicMotifEvidenceCoverageStatus::CompleteForRequestedThreshold
                | GenomicMotifEvidenceCoverageStatus::DensityLimitedAtSource
        )
    })
}

fn chromosome_aliases(chromosome: &str) -> Vec<String> {
    let chromosome = chromosome.trim();
    let mut aliases = vec![chromosome.to_string()];
    if let Some(unprefixed) = chromosome.strip_prefix("chr") {
        aliases.push(unprefixed.to_string());
    } else {
        aliases.push(format!("chr{chromosome}"));
    }
    if chromosome.eq_ignore_ascii_case("MT") || chromosome.eq_ignore_ascii_case("chrM") {
        aliases.extend(["MT".to_string(), "M".to_string(), "chrM".to_string()]);
    }
    aliases.sort();
    aliases.dedup();
    aliases
}

fn resolve_payload_path(root: &Path, row: &InventoryRow) -> Result<PathBuf, String> {
    let relative = Path::new(&row.output_relative_path);
    if relative.is_absolute() {
        return Err(format!(
            "inventory payload path is absolute: {}",
            row.output_relative_path
        ));
    }
    let candidate = root
        .join("task_data")
        .join(format!("task_id={}", row.task_id))
        .join(relative);
    if !candidate.is_file() {
        return Err(format!(
            "inventory payload is missing for task '{}': {}",
            row.task_id,
            candidate.display()
        ));
    }
    let canonical = candidate.canonicalize().map_err(|error| {
        format!(
            "could not resolve payload '{}': {error}",
            candidate.display()
        )
    })?;
    if !path_is_within(root, &canonical) {
        return Err(format!(
            "inventory payload resolves outside package root: {}",
            candidate.display()
        ));
    }
    Ok(canonical)
}

fn source_coordinates(
    region: &GenomicMotifQueryRegion,
    hit_start: u64,
    hit_end: u64,
    hit_strand: &str,
) -> (Option<usize>, Option<usize>, Option<bool>) {
    let (Some(source_start), Some(source_end)) = (
        region.source_start_0based,
        region.source_end_0based_exclusive,
    ) else {
        return (None, None, None);
    };
    if hit_start < region.start_0based || hit_end > region.end_0based_exclusive {
        return (None, None, None);
    }
    let offset_start = usize::try_from(hit_start - region.start_0based).ok();
    let offset_end = usize::try_from(hit_end - region.start_0based).ok();
    let (Some(offset_start), Some(offset_end)) = (offset_start, offset_end) else {
        return (None, None, None);
    };
    let coordinates = if region.source_anchor_reverse {
        (
            source_end.checked_sub(offset_end),
            source_end.checked_sub(offset_start),
        )
    } else {
        (
            source_start.checked_add(offset_start),
            source_start.checked_add(offset_end),
        )
    };
    let local_forward = match hit_strand {
        "+" => Some(!region.source_anchor_reverse),
        "-" => Some(region.source_anchor_reverse),
        _ => None,
    };
    (coordinates.0, coordinates.1, local_forward)
}

pub(crate) fn query_genomic_motif_evidence(
    request: &GenomicMotifEvidenceRequest,
    regions: &[GenomicMotifQueryRegion],
) -> Result<GenomicMotifEvidenceReport, String> {
    validate_genomic_motif_request(request, regions)?;
    let paths = match resolve_package_paths(request) {
        Ok(paths) => paths,
        Err(error) if error == "package_not_configured" => {
            return Ok(unavailable_report(
                request,
                regions,
                GenomicMotifEvidenceAvailability::PackageNotConfigured,
                format!(
                    "No genome motif package was requested. Set request.package_root or {JASPAR_GENOME_SCAN_PACKAGE_ENV}; local GENtle TFBS scoring remains available."
                ),
            ));
        }
        Err(error) if error.starts_with("package_missing:") => {
            return Ok(unavailable_report(
                request,
                regions,
                GenomicMotifEvidenceAvailability::PackageMissing,
                format!(
                    "Genome motif package directory is missing: {}; local GENtle TFBS scoring remains available",
                    error.trim_start_matches("package_missing:")
                ),
            ));
        }
        Err(error) => {
            return Ok(unavailable_report(
                request,
                regions,
                GenomicMotifEvidenceAvailability::InvalidPackage,
                format!(
                    "{}; local GENtle TFBS scoring remains available",
                    error.trim_start_matches("invalid_package:")
                ),
            ));
        }
    };
    let executable = request
        .duckdb_executable
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .or_else(|| {
            env::var(DUCKDB_BIN_ENV)
                .ok()
                .filter(|value| !value.trim().is_empty())
        })
        .unwrap_or_else(|| "duckdb".to_string());
    let timeout = Duration::from_secs(request.timeout_seconds);
    let version_output = match run_duckdb(&executable, None, None, timeout) {
        Ok(output) if output.output.status.success() && !output.stdout_truncated => {
            String::from_utf8_lossy(&output.output.stdout)
                .lines()
                .find(|line| !line.trim().is_empty())
                .map(|line| line.trim().to_string())
        }
        Ok(output) => {
            let diagnostic = String::from_utf8_lossy(&output.output.stderr);
            return Ok(unavailable_report(
                request,
                regions,
                GenomicMotifEvidenceAvailability::DuckdbUnavailable,
                format!(
                    "DuckDB version probe failed with status {:?}: {}; local GENtle TFBS scoring remains available",
                    output.output.status.code(),
                    diagnostic.trim()
                ),
            ));
        }
        Err(error) => {
            return Ok(unavailable_report(
                request,
                regions,
                GenomicMotifEvidenceAvailability::DuckdbUnavailable,
                format!("{error}; local GENtle TFBS scoring remains available"),
            ));
        }
    };

    let metadata_query = "SELECT r.run_id::VARCHAR AS run_id, r.genome_id::VARCHAR AS genome_id, r.motif_set_id::VARCHAR AS motif_set_id, r.score_mode::VARCHAR AS score_mode, CAST(r.pseudocount AS DOUBLE) AS pseudocount, r.pseudocount_scheme::VARCHAR AS pseudocount_scheme, r.background_model_id::VARCHAR AS background_model_id, CAST(r.minimum_pwm_relative_score AS DOUBLE) AS minimum_pwm_relative_score, CAST(r.maximum_pwm_relative_score AS DOUBLE) AS maximum_pwm_relative_score, r.coordinate_mode::VARCHAR AS coordinate_mode, r.n_policy::VARCHAR AS n_policy, r.matched_sequence_policy::VARCHAR AS matched_sequence_policy, g.assembly_name::VARCHAR AS assembly_name, g.assembly_accession::VARCHAR AS assembly_accession, CAST(g.ensembl_release AS VARCHAR) AS ensembl_release, m.jaspar_version::VARCHAR AS jaspar_version FROM scan_run r JOIN genome g USING (genome_id) JOIN motif_set m USING (motif_set_id) ORDER BY r.run_id LIMIT 2;";
    let metadata_rows = match query_json(&executable, &paths.database_path, metadata_query, timeout)
    {
        Ok(rows) if rows.len() == 1 => rows,
        Ok(rows) => {
            return Ok(unavailable_report(
                request,
                regions,
                GenomicMotifEvidenceAvailability::InvalidPackage,
                format!(
                    "Expected exactly one scan_run metadata row, found {}; local GENtle TFBS scoring remains available",
                    rows.len()
                ),
            ));
        }
        Err(error) => {
            return Ok(unavailable_report(
                request,
                regions,
                GenomicMotifEvidenceAvailability::QueryFailed,
                format!("{error}; local GENtle TFBS scoring remains available"),
            ));
        }
    };
    let metadata = &metadata_rows[0];
    let package_genome_id = required_string(metadata, "genome_id")?;
    if let Some(expected) = request
        .expected_genome_id
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        && expected != package_genome_id
    {
        return Ok(unavailable_report(
            request,
            regions,
            GenomicMotifEvidenceAvailability::IncompatiblePackage,
            format!(
                "Requested genome_id '{expected}' does not match package genome_id '{package_genome_id}'; local GENtle TFBS scoring remains available"
            ),
        ));
    }

    let threshold_query = format!(
        "SELECT t.motif_id::VARCHAR AS motif_id, m.motif_name::VARCHAR AS motif_name, t.threshold_set_id::VARCHAR AS threshold_set_id, CAST(t.informative_threshold AS DOUBLE) AS informative_threshold, CAST(t.final_minimum_score AS DOUBLE) AS final_minimum_score, t.density_limited::BOOLEAN AS density_limited FROM scan_motif_threshold t LEFT JOIN motif_metadata m USING (motif_set_id, motif_id) WHERE t.motif_id IN ({}) ORDER BY t.motif_id;",
        request
            .motif_ids
            .iter()
            .map(|value| sql_string(value.trim()))
            .collect::<Vec<_>>()
            .join(",")
    );
    let threshold_values =
        query_json(&executable, &paths.database_path, &threshold_query, timeout)?;
    let mut thresholds = BTreeMap::new();
    for row in &threshold_values {
        let threshold = MotifThresholdRow {
            motif_id: required_string(row, "motif_id")?,
            motif_name: optional_string(row, "motif_name"),
            threshold_set_id: optional_string(row, "threshold_set_id"),
            informative_threshold: optional_f64(row, "informative_threshold"),
            final_minimum_score: optional_f64(row, "final_minimum_score"),
            density_limited: optional_bool(row, "density_limited"),
        };
        if thresholds
            .insert(threshold.motif_id.clone(), threshold)
            .is_some()
        {
            return Err("DuckDB returned duplicate motif threshold rows".to_string());
        }
    }
    let source_minimum_pwm_relative_score = optional_f64(metadata, "minimum_pwm_relative_score");

    let mut all_aliases = BTreeSet::new();
    for region in regions {
        all_aliases.extend(chromosome_aliases(&region.chromosome));
    }
    let contig_query = format!(
        "SELECT CAST(chrom AS VARCHAR) AS chrom, CAST(length AS UBIGINT) AS length FROM sequence_region WHERE CAST(chrom AS VARCHAR) IN ({} ) ORDER BY chrom;",
        all_aliases
            .iter()
            .map(|value| sql_string(value))
            .collect::<Vec<_>>()
            .join(",")
    );
    let contig_rows = query_json(&executable, &paths.database_path, &contig_query, timeout)?;
    let contigs = contig_rows
        .iter()
        .filter_map(|row| {
            Some((
                required_string(row, "chrom").ok()?,
                optional_u64(row, "length")?,
            ))
        })
        .collect::<BTreeMap<_, _>>();
    let mut resolved_regions = vec![];
    let mut query_regions = vec![];
    let mut warnings = vec![];
    for region in regions {
        let aliases = chromosome_aliases(&region.chromosome);
        let resolved = aliases.iter().find_map(|alias| {
            contigs
                .get(alias)
                .copied()
                .map(|length| (alias.clone(), length))
        });
        let (resolved_chromosome, contig_length, status) = match resolved {
            Some((chromosome, length)) => {
                let anchor_geometry_ok = region.end_0based_exclusive <= length
                    && match (
                        region.source_sequence_length_bp,
                        region.source_anchor_start_1based,
                        region.source_anchor_end_1based,
                    ) {
                        (Some(sequence_length), Some(anchor_start), Some(anchor_end)) => {
                            anchor_end >= anchor_start
                                && anchor_end - anchor_start + 1 == sequence_length
                        }
                        _ => true,
                    };
                (
                    Some(chromosome),
                    Some(length),
                    if anchor_geometry_ok {
                        GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMatchedOnly
                    } else {
                        GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMismatch
                    },
                )
            }
            None => (
                None,
                None,
                GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMismatch,
            ),
        };
        if status == GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMismatch {
            warnings.push(format!(
                "Interval '{}' could not be matched safely to package contig '{}' geometry",
                region.interval_id, region.chromosome
            ));
        } else if let Some(chromosome) = &resolved_chromosome {
            query_regions.push((region.clone(), chromosome.clone()));
        }
        resolved_regions.push(GenomicMotifEvidenceResolvedRegion {
            interval_id: region.interval_id.clone(),
            label: region.label.clone(),
            requested_chromosome: region.chromosome.clone(),
            resolved_chromosome,
            start_0based: region.start_0based,
            end_0based_exclusive: region.end_0based_exclusive,
            source_seq_id: region.source_seq_id.clone(),
            source_start_0based: region.source_start_0based,
            source_end_0based_exclusive: region.source_end_0based_exclusive,
            compatibility_status: status,
            package_contig_length_bp: contig_length,
        });
    }
    if query_regions.is_empty() {
        let mut report = unavailable_report(
            request,
            regions,
            GenomicMotifEvidenceAvailability::IncompatiblePackage,
            "No requested interval matched package contig geometry; local GENtle TFBS scoring remains available",
        );
        report.regions = resolved_regions;
        report.warnings.extend(warnings);
        return Ok(report);
    }
    let regions_complete = resolved_regions.iter().all(|region| {
        region.compatibility_status
            == GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMatchedOnly
    });

    let selected_chromosomes = query_regions
        .iter()
        .map(|(_, chromosome)| chromosome.as_str())
        .collect::<BTreeSet<_>>();
    let inventory_query = format!(
        "SELECT task_id::VARCHAR AS task_id, output_relative_path::VARCHAR AS output_relative_path, sha256::VARCHAR AS sha256, CAST(emitted_hits AS UBIGINT) AS emitted_hits FROM scan_file_inventory WHERE motif_id IN ({}) AND CAST(chrom AS VARCHAR) IN ({}) ORDER BY motif_id, chrom, strand, task_id, output_relative_path;",
        request
            .motif_ids
            .iter()
            .map(|value| sql_string(value.trim()))
            .collect::<Vec<_>>()
            .join(","),
        selected_chromosomes
            .iter()
            .map(|value| sql_string(value))
            .collect::<Vec<_>>()
            .join(",")
    );
    let inventory_values =
        query_json(&executable, &paths.database_path, &inventory_query, timeout)?;
    let inventory_rows = inventory_values
        .iter()
        .map(|row| {
            Ok(InventoryRow {
                task_id: required_string(row, "task_id")?,
                output_relative_path: required_string(row, "output_relative_path")?,
                sha256: required_string(row, "sha256")?,
                emitted_hits: optional_u64(row, "emitted_hits").unwrap_or(0),
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    if inventory_rows.len() > request.max_payload_files {
        let mut report = unavailable_report(
            request,
            regions,
            GenomicMotifEvidenceAvailability::QueryFailed,
            format!(
                "Query selected {} payload files, exceeding max_payload_files={}; reduce motifs or intervals",
                inventory_rows.len(),
                request.max_payload_files
            ),
        );
        report.regions = resolved_regions;
        return Ok(report);
    }
    let selected_payload_emitted_hit_count = inventory_rows
        .iter()
        .map(|row| row.emitted_hits)
        .sum::<u64>();
    if inventory_rows.is_empty() {
        let content_fingerprint = declared_content_fingerprint(&paths, &inventory_rows);
        let motif_coverage = motif_coverage_rows(
            request,
            &thresholds,
            source_minimum_pwm_relative_score,
            &[],
            false,
        );
        let mut no_payload_warnings = motif_coverage_warnings(&motif_coverage);
        no_payload_warnings.push(
            "No inventory payloads matched the requested motif/contig set; this is not evidence that the motifs are biologically absent outside the package's retained-score policy"
                .to_string(),
        );
        return Ok(GenomicMotifEvidenceReport {
            schema: GENOMIC_MOTIF_EVIDENCE_SCHEMA.to_string(),
            report_id: report_id(request, Some(&content_fingerprint)),
            availability: GenomicMotifEvidenceAvailability::Available,
            request: request.clone(),
            provider: Some(provider_provenance(
                &paths,
                &executable,
                version_output,
                metadata,
                &inventory_rows,
            )?),
            regions: resolved_regions,
            query_complete: regions_complete && motif_coverage_is_query_complete(&motif_coverage),
            motif_coverage,
            warnings: no_payload_warnings,
            ..GenomicMotifEvidenceReport::default()
        });
    }
    let payload_paths = inventory_rows
        .iter()
        .map(|row| resolve_payload_path(&paths.root, row))
        .collect::<Result<Vec<_>, _>>()?;
    let region_values = query_regions
        .iter()
        .map(|(region, chromosome)| {
            format!(
                "({}, {}, {}::UBIGINT, {}::UBIGINT)",
                sql_string(&region.interval_id),
                sql_string(chromosome),
                region.start_0based,
                region.end_0based_exclusive
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let mut filters = vec!["TRUE".to_string()];
    if let Some(minimum) = request.minimum_score {
        filters.push(format!("h.score >= {minimum:.17}"));
    }
    if let Some(minimum) = request.minimum_pwm_relative_score {
        filters.push(format!("h.pwm_relative_score >= {minimum:.17}"));
    }
    let payload_path_strings = payload_paths
        .iter()
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>();
    let hit_query = format!(
        "WITH requested_regions(interval_id, chrom, region_start, region_end) AS (VALUES {region_values}) SELECT r.interval_id, h.chrom, h.start, h.\"end\", h.motif_id, h.motif_name, h.strand, h.score, h.pwm_relative_score, h.score_mode, h.minimum_score, h.matched_seq FROM motif_hit_files({paths}) h JOIN requested_regions r ON h.chrom = r.chrom AND h.\"end\" > r.region_start AND h.start < r.region_end WHERE {filters} ORDER BY r.interval_id, h.chrom, h.start, h.\"end\", h.motif_id, h.strand LIMIT {limit};",
        paths = sql_string_list(payload_path_strings.iter().map(String::as_str)),
        filters = filters.join(" AND "),
        limit = request.max_rows.saturating_add(1),
    );
    let mut hit_values = query_json(&executable, &paths.database_path, &hit_query, timeout)?;
    let truncated = hit_values.len() > request.max_rows;
    if truncated {
        hit_values.truncate(request.max_rows);
    }
    let regions_by_id = query_regions
        .iter()
        .map(|(region, _)| (region.interval_id.as_str(), region))
        .collect::<BTreeMap<_, _>>();
    let hits = hit_values
        .iter()
        .map(|row| {
            let interval_id = required_string(row, "interval_id")?;
            let region = regions_by_id
                .get(interval_id.as_str())
                .ok_or_else(|| format!("DuckDB returned unknown interval_id '{interval_id}'"))?;
            let start = optional_u64(row, "start")
                .ok_or_else(|| "DuckDB hit is missing start".to_string())?;
            let end =
                optional_u64(row, "end").ok_or_else(|| "DuckDB hit is missing end".to_string())?;
            let strand = required_string(row, "strand")?;
            let (source_start, source_end, source_forward) =
                source_coordinates(region, start, end, &strand);
            Ok(GenomicMotifEvidenceHit {
                interval_id,
                chromosome: required_string(row, "chrom")?,
                start_0based: start,
                end_0based_exclusive: end,
                motif_id: required_string(row, "motif_id")?,
                motif_name: optional_string(row, "motif_name"),
                strand,
                score: optional_f64(row, "score")
                    .ok_or_else(|| "DuckDB hit is missing numeric score".to_string())?,
                pwm_relative_score: optional_f64(row, "pwm_relative_score"),
                score_mode: required_string(row, "score_mode")?,
                minimum_score: optional_f64(row, "minimum_score"),
                matched_sequence: optional_string(row, "matched_seq"),
                source_seq_id: region.source_seq_id.clone(),
                source_start_0based: source_start,
                source_end_0based_exclusive: source_end,
                source_forward_strand: source_forward,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    if truncated {
        warnings.push(format!(
            "Hit rows were truncated at max_rows={}; narrow the regions, motifs, or score threshold for a complete row report",
            request.max_rows
        ));
    }
    let motif_coverage = motif_coverage_rows(
        request,
        &thresholds,
        source_minimum_pwm_relative_score,
        &hits,
        truncated,
    );
    let query_complete = regions_complete && motif_coverage_is_query_complete(&motif_coverage);
    warnings.extend(motif_coverage_warnings(&motif_coverage));
    warnings.push(
        "Package retention floors are storage/calibration policy, not biological binding thresholds; package scores are not placed on the same axis as GENtle's local scorer"
            .to_string(),
    );
    if resolved_regions.iter().any(|region| {
        region.compatibility_status
            == GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMatchedOnly
    }) {
        warnings.push(
            "Compatibility is limited to contig alias, interval containment, and anchor geometry; GENtle has not verified the package per-contig sequence SHA-256"
                .to_string(),
        );
    }
    let content_fingerprint = declared_content_fingerprint(&paths, &inventory_rows);
    Ok(GenomicMotifEvidenceReport {
        schema: GENOMIC_MOTIF_EVIDENCE_SCHEMA.to_string(),
        report_id: report_id(request, Some(&content_fingerprint)),
        availability: GenomicMotifEvidenceAvailability::Available,
        request: request.clone(),
        provider: Some(provider_provenance(
            &paths,
            &executable,
            version_output,
            metadata,
            &inventory_rows,
        )?),
        regions: resolved_regions,
        motif_coverage,
        selected_payload_file_count: inventory_rows.len(),
        selected_payload_emitted_hit_count,
        matched_hit_count: hits.len(),
        returned_hit_count: hits.len(),
        truncated,
        query_complete,
        hits,
        warnings,
        ..GenomicMotifEvidenceReport::default()
    })
}

fn provider_provenance(
    paths: &PackagePaths,
    executable: &str,
    version: Option<String>,
    metadata: &Value,
    inventory_rows: &[InventoryRow],
) -> Result<GenomicMotifEvidenceProviderProvenance, String> {
    Ok(GenomicMotifEvidenceProviderProvenance {
        provider_kind: "precomputed_sparse_genome_scan_package".to_string(),
        package_root: paths.root.display().to_string(),
        manifest_path: paths.manifest_path.display().to_string(),
        manifest_sha256: sha256_prefixed_bytes(&paths.manifest_bytes),
        declared_content_fingerprint_sha256: declared_content_fingerprint(paths, inventory_rows),
        manifest_schema_version: paths
            .manifest
            .get("schema_version")
            .and_then(Value::as_u64)
            .unwrap_or_default(),
        database_path: paths.database_path.display().to_string(),
        duckdb_executable: executable.to_string(),
        duckdb_version: version,
        run_id: required_string(metadata, "run_id")?,
        genome_id: required_string(metadata, "genome_id")?,
        assembly_name: optional_string(metadata, "assembly_name"),
        assembly_accession: optional_string(metadata, "assembly_accession"),
        ensembl_release: optional_string(metadata, "ensembl_release"),
        motif_set_id: required_string(metadata, "motif_set_id")?,
        jaspar_version: optional_string(metadata, "jaspar_version"),
        score_mode: required_string(metadata, "score_mode")?,
        pseudocount: optional_f64(metadata, "pseudocount").unwrap_or_default(),
        pseudocount_scheme: required_string(metadata, "pseudocount_scheme")?,
        background_model_id: required_string(metadata, "background_model_id")?,
        minimum_pwm_relative_score: optional_f64(metadata, "minimum_pwm_relative_score"),
        maximum_pwm_relative_score: optional_f64(metadata, "maximum_pwm_relative_score"),
        coordinate_mode: required_string(metadata, "coordinate_mode")?,
        n_policy: required_string(metadata, "n_policy")?,
        matched_sequence_policy: required_string(metadata, "matched_sequence_policy")?,
        selected_payloads: inventory_rows
            .iter()
            .map(|row| GenomicMotifEvidencePayloadProvenance {
                task_id: row.task_id.clone(),
                output_relative_path: row.output_relative_path.clone(),
                declared_sha256: row.sha256.clone(),
                emitted_hits: row.emitted_hits,
            })
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(unix)]
    fn synthetic_package_with_fake_duckdb() -> (tempfile::TempDir, PathBuf) {
        use std::os::unix::fs::PermissionsExt;

        // Hand-authored synthetic package used only to exercise the provider contract.
        let package = tempfile::tempdir().expect("temporary package");
        fs::write(
            package.path().join("manifest.json"),
            r#"{"schema_version":2,"state":"complete","database":"jaspar_genome_scan.duckdb"}"#,
        )
        .expect("manifest");
        fs::write(
            package.path().join("jaspar_genome_scan.duckdb"),
            b"synthetic",
        )
        .expect("database placeholder");
        let payload_path = package
            .path()
            .join("task_data")
            .join("task_id=task-1")
            .join("hits.parquet");
        fs::create_dir_all(payload_path.parent().expect("payload parent"))
            .expect("payload directory");
        fs::write(&payload_path, b"synthetic parquet placeholder").expect("payload placeholder");
        let executable = package.path().join("fake_duckdb.sh");
        fs::write(
            &executable,
            r#"#!/bin/sh
if [ "$1" = "--version" ]; then
  printf '%s\n' 'DuckDB synthetic 1.0'
  exit 0
fi
case "$5" in
  *"FROM scan_run r"*)
    printf '%s\n' '[{"run_id":"synthetic_run","genome_id":"synthetic_grch38","motif_set_id":"synthetic_jaspar","score_mode":"log2_relative_risk","pseudocount":1.0,"pseudocount_scheme":"additive_per_base","background_model_id":"uniform_acgt_v1","minimum_pwm_relative_score":null,"maximum_pwm_relative_score":null,"coordinate_mode":"bed_0based_half_open","n_policy":"skip","matched_sequence_policy":"forward_reference","assembly_name":"GRCh38","assembly_accession":"GCA_000001405.15","ensembl_release":"116","jaspar_version":"2026"}]'
    ;;
  *"FROM scan_motif_threshold"*)
    printf '%s\n' '[{"motif_id":"MA0525.2","motif_name":"TP63","threshold_set_id":"synthetic_thresholds","informative_threshold":2.0,"final_minimum_score":-1.0,"density_limited":false}]'
    ;;
  *"FROM sequence_region"*)
    printf '%s\n' '[{"chrom":"1","length":1000}]'
    ;;
  *"FROM scan_file_inventory"*)
    printf '%s\n' '[{"task_id":"task-1","output_relative_path":"hits.parquet","sha256":"sha256:synthetic-payload","emitted_hits":2}]'
    ;;
  *"motif_hit_files"*)
    printf '%s\n' '[{"interval_id":"region-1","chrom":"1","start":12,"end":20,"motif_id":"MA0525.2","motif_name":"TP63","strand":"+","score":3.25,"pwm_relative_score":0.91,"score_mode":"log2_relative_risk","minimum_score":0.0,"matched_seq":"ACGTACGT"}]'
    ;;
  *)
    printf '%s\n' '[]'
    ;;
esac
"#,
        )
        .expect("fake DuckDB");
        let mut permissions = fs::metadata(&executable)
            .expect("fake metadata")
            .permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&executable, permissions).expect("fake executable permissions");
        (package, executable)
    }

    #[test]
    fn request_rejects_unbounded_all_motif_row_queries() {
        let request = GenomicMotifEvidenceRequest {
            motif_ids: vec!["ALL".to_string()],
            ..GenomicMotifEvidenceRequest::default()
        };
        let regions = vec![GenomicMotifQueryRegion {
            interval_id: "region-1".to_string(),
            label: None,
            chromosome: "1".to_string(),
            start_0based: 10,
            end_0based_exclusive: 20,
            source_seq_id: None,
            source_start_0based: None,
            source_end_0based_exclusive: None,
            source_sequence_length_bp: None,
            source_anchor_start_1based: None,
            source_anchor_end_1based: None,
            source_anchor_reverse: false,
        }];
        let error = validate_genomic_motif_request(&request, &regions).expect_err("ALL fails");
        assert!(error.contains("explicit bounded motif set"));
    }

    #[test]
    fn reverse_anchor_maps_forward_reference_hit_to_local_coordinates() {
        let region = GenomicMotifQueryRegion {
            interval_id: "reverse".to_string(),
            label: None,
            chromosome: "22".to_string(),
            start_0based: 1_000,
            end_0based_exclusive: 1_100,
            source_seq_id: Some("seq".to_string()),
            source_start_0based: Some(0),
            source_end_0based_exclusive: Some(100),
            source_sequence_length_bp: Some(100),
            source_anchor_start_1based: Some(1_001),
            source_anchor_end_1based: Some(1_100),
            source_anchor_reverse: true,
        };
        assert_eq!(
            source_coordinates(&region, 1_010, 1_020, "+"),
            (Some(80), Some(90), Some(false))
        );
        assert_eq!(
            source_coordinates(&region, 1_010, 1_020, "-"),
            (Some(80), Some(90), Some(true))
        );
    }

    #[test]
    fn motif_coverage_distinguishes_storage_floor_density_and_missing_motif() {
        let thresholds = BTreeMap::from([
            (
                "MA0001.1".to_string(),
                MotifThresholdRow {
                    motif_id: "MA0001.1".to_string(),
                    motif_name: Some("TF1".to_string()),
                    threshold_set_id: Some("thresholds-v1".to_string()),
                    informative_threshold: Some(2.0),
                    final_minimum_score: Some(-1.0),
                    density_limited: Some(false),
                },
            ),
            (
                "MA0002.1".to_string(),
                MotifThresholdRow {
                    motif_id: "MA0002.1".to_string(),
                    motif_name: Some("TF2".to_string()),
                    threshold_set_id: Some("thresholds-v1".to_string()),
                    informative_threshold: Some(-4.0),
                    final_minimum_score: Some(0.5),
                    density_limited: Some(true),
                },
            ),
        ]);
        let request = GenomicMotifEvidenceRequest {
            motif_ids: vec![
                "MA0001.1".to_string(),
                "MA0002.1".to_string(),
                "MA9999.1".to_string(),
            ],
            minimum_score: Some(-2.0),
            ..GenomicMotifEvidenceRequest::default()
        };
        let rows = motif_coverage_rows(&request, &thresholds, None, &[], false);
        assert_eq!(
            rows[0].status,
            GenomicMotifEvidenceCoverageStatus::IncompleteBelowStorageFloor
        );
        assert_eq!(
            rows[1].status,
            GenomicMotifEvidenceCoverageStatus::IncompleteBelowStorageFloor
        );
        assert_eq!(
            rows[2].status,
            GenomicMotifEvidenceCoverageStatus::MotifNotInPackage
        );

        let density_rows = motif_coverage_rows(
            &GenomicMotifEvidenceRequest {
                motif_ids: vec!["MA0002.1".to_string()],
                ..GenomicMotifEvidenceRequest::default()
            },
            &thresholds,
            None,
            &[],
            false,
        );
        assert_eq!(
            density_rows[0].status,
            GenomicMotifEvidenceCoverageStatus::DensityLimitedAtSource
        );
    }

    #[test]
    fn missing_package_is_an_optional_unavailable_report() {
        let request = GenomicMotifEvidenceRequest {
            package_root: Some("/definitely/absent/gentle-genome-scan-package".to_string()),
            motif_ids: vec!["MA0525.2".to_string()],
            ..GenomicMotifEvidenceRequest::default()
        };
        let regions = vec![GenomicMotifQueryRegion {
            interval_id: "region-1".to_string(),
            label: None,
            chromosome: "1".to_string(),
            start_0based: 10,
            end_0based_exclusive: 20,
            source_seq_id: None,
            source_start_0based: None,
            source_end_0based_exclusive: None,
            source_sequence_length_bp: None,
            source_anchor_start_1based: None,
            source_anchor_end_1based: None,
            source_anchor_reverse: false,
        }];
        let report = query_genomic_motif_evidence(&request, &regions).expect("typed report");
        assert_eq!(
            report.availability,
            GenomicMotifEvidenceAvailability::PackageMissing
        );
        assert!(report.hits.is_empty());
        assert!(report.warnings[0].contains("local GENtle TFBS scoring remains available"));
    }

    #[cfg(unix)]
    #[test]
    fn exact_inventory_query_is_deterministic_and_content_bound() {
        let (package, executable) = synthetic_package_with_fake_duckdb();
        let request = GenomicMotifEvidenceRequest {
            package_root: Some(package.path().display().to_string()),
            duckdb_executable: Some(executable.display().to_string()),
            expected_genome_id: Some("synthetic_grch38".to_string()),
            motif_ids: vec!["MA0525.2".to_string()],
            ..GenomicMotifEvidenceRequest::default()
        };
        let regions = vec![GenomicMotifQueryRegion {
            interval_id: "region-1".to_string(),
            label: Some("synthetic target".to_string()),
            chromosome: "chr1".to_string(),
            start_0based: 10,
            end_0based_exclusive: 30,
            source_seq_id: None,
            source_start_0based: None,
            source_end_0based_exclusive: None,
            source_sequence_length_bp: None,
            source_anchor_start_1based: None,
            source_anchor_end_1based: None,
            source_anchor_reverse: false,
        }];

        let first = query_genomic_motif_evidence(&request, &regions).expect("first query");
        let second = query_genomic_motif_evidence(&request, &regions).expect("second query");

        assert_eq!(first, second);
        assert_eq!(
            first.availability,
            GenomicMotifEvidenceAvailability::Available
        );
        assert_eq!(first.selected_payload_file_count, 1);
        assert_eq!(first.returned_hit_count, 1);
        assert_eq!(first.motif_coverage.len(), 1);
        assert_eq!(
            first.motif_coverage[0].status,
            GenomicMotifEvidenceCoverageStatus::CompleteForPackageRetention
        );
        assert_eq!(first.motif_coverage[0].source_minimum_score, Some(-1.0));
        assert_eq!(first.hits[0].start_0based, 12);
        assert_eq!(first.hits[0].motif_id, "MA0525.2");
        let provider = first.provider.expect("provider provenance");
        assert_eq!(provider.selected_payloads.len(), 1);
        assert_eq!(
            provider.selected_payloads[0].declared_sha256,
            "sha256:synthetic-payload"
        );
        assert!(
            provider
                .declared_content_fingerprint_sha256
                .starts_with("sha256:")
        );
    }

    #[cfg(unix)]
    #[test]
    fn missing_duckdb_is_an_optional_unavailable_report() {
        let (package, _) = synthetic_package_with_fake_duckdb();
        let request = GenomicMotifEvidenceRequest {
            package_root: Some(package.path().display().to_string()),
            duckdb_executable: Some(
                package
                    .path()
                    .join("not-installed-duckdb")
                    .display()
                    .to_string(),
            ),
            motif_ids: vec!["MA0525.2".to_string()],
            ..GenomicMotifEvidenceRequest::default()
        };
        let regions = vec![GenomicMotifQueryRegion {
            interval_id: "region-1".to_string(),
            label: None,
            chromosome: "1".to_string(),
            start_0based: 10,
            end_0based_exclusive: 20,
            source_seq_id: None,
            source_start_0based: None,
            source_end_0based_exclusive: None,
            source_sequence_length_bp: None,
            source_anchor_start_1based: None,
            source_anchor_end_1based: None,
            source_anchor_reverse: false,
        }];

        let report = query_genomic_motif_evidence(&request, &regions).expect("typed report");
        assert_eq!(
            report.availability,
            GenomicMotifEvidenceAvailability::DuckdbUnavailable
        );
        assert!(report.hits.is_empty());
        assert!(report.warnings[0].contains("local GENtle TFBS scoring remains available"));
    }
}
