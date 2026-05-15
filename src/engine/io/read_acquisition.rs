//! Shared SRA-backed read acquisition for RNA and CUT&RUN workflows.
//!
//! This module owns the deterministic setup boundary around external SRA
//! Toolkit commands. Downstream analysis code consumes prepared FASTA/FASTQ
//! paths and reports from here rather than each workflow learning how to
//! download or dump SRA runs independently.

use super::*;
use sha1::{Digest, Sha1};
use std::fs;
#[cfg(unix)]
use std::os::unix::process::CommandExt;

const READ_ACQUIRE_ACTIVITY_STATUS_FILE: &str = ".read_acquisition_activity.json";
const READ_ACQUIRE_ACTIVITY_LOCK_FILE: &str = ".read_acquisition_activity.lock";
const READ_ACQUIRE_CANCEL_FILE: &str = ".read_acquisition_cancel";
const READ_ACQUIRE_ACTIVITY_STALE_AFTER_MS: u128 = 6 * 60 * 60 * 1000;
const READ_ACQUIRE_MONITOR_INTERVAL_MS: u64 = 500;

#[derive(Debug, Clone)]
struct ReadAcquisitionEffectiveRow {
    manifest: ReadAcquisitionManifestRow,
    read_layout: ReadAcquisitionReadLayout,
    analysis_format: ReadAcquisitionAnalysisFormat,
}

enum ReadAcquisitionActivityStart {
    Acquired(ReadAcquisitionActivityTracker),
    Running,
}

struct ReadAcquisitionActivityTracker {
    status_path: PathBuf,
    lock_path: PathBuf,
    cancel_path: PathBuf,
    status: SharedAssetActivityStatus,
}

impl ReadAcquisitionActivityTracker {
    fn start(run_dir: &Path, sra_accession: &str) -> Result<ReadAcquisitionActivityStart, String> {
        fs::create_dir_all(run_dir).map_err(|e| {
            format!(
                "Could not create read-acquisition run directory '{}': {e}",
                run_dir.display()
            )
        })?;
        let status_path = read_acquire_activity_status_path(run_dir);
        let lock_path = read_acquire_activity_lock_path(run_dir);
        let cancel_path = read_acquire_cancel_path(run_dir);
        if let Some(existing) =
            inspect_read_acquire_activity_status_paths(&status_path, &lock_path)?
        {
            if existing.lifecycle_status == "running" {
                return Ok(ReadAcquisitionActivityStart::Running);
            }
        }

        let now = GentleEngine::now_unix_ms();
        let status = SharedAssetActivityStatus {
            resource_key: read_acquisition_resource_key(sra_accession),
            display_name: sra_accession.to_string(),
            status_path: read_acquisition_path_string(&status_path),
            lock_path: Some(read_acquisition_path_string(&lock_path)),
            cancel_path: Some(read_acquisition_path_string(&cancel_path)),
            lifecycle_status: "running".to_string(),
            phase: Some("prefetch".to_string()),
            item: Some(sra_accession.to_string()),
            started_at_unix_ms: now,
            updated_at_unix_ms: now,
            owner_pid: Some(std::process::id()),
            ..Default::default()
        };
        let tracker = Self {
            status_path,
            lock_path,
            cancel_path,
            status,
        };
        if !create_read_acquire_activity_lock(&tracker.lock_path, &tracker.status)? {
            if let Some(existing) = inspect_read_acquire_activity_status_paths(
                &tracker.status_path,
                &tracker.lock_path,
            )? {
                if existing.lifecycle_status == "running" {
                    return Ok(ReadAcquisitionActivityStart::Running);
                }
            }
            return Err(format!(
                "Could not acquire read-acquisition lock '{}'",
                tracker.lock_path.display()
            ));
        }
        remove_read_acquire_cancel_marker(&tracker.cancel_path);
        tracker.write_status();
        Ok(ReadAcquisitionActivityStart::Acquired(tracker))
    }

    fn update_phase(&mut self, phase: &str) {
        self.status.lifecycle_status = "running".to_string();
        self.status.phase = Some(phase.to_string());
        self.status.cancel_path = Some(read_acquisition_path_string(&self.cancel_path));
        self.status.updated_at_unix_ms = GentleEngine::now_unix_ms();
        self.write_status();
    }

    fn update_monitor(
        &mut self,
        phase: &str,
        produced_bytes: u64,
        monitored_free_bytes: Option<u64>,
        minimum_free_bytes: Option<u64>,
    ) {
        self.status.lifecycle_status = "running".to_string();
        self.status.phase = Some(phase.to_string());
        self.status.bytes_done = produced_bytes;
        self.status.bytes_total = None;
        self.status.percent = None;
        self.status.monitored_free_bytes = monitored_free_bytes;
        self.status.minimum_free_bytes = minimum_free_bytes;
        self.status.cancel_path = Some(read_acquisition_path_string(&self.cancel_path));
        self.status.updated_at_unix_ms = GentleEngine::now_unix_ms();
        self.write_status();
    }

    fn finish_failure(&mut self, message: &str) {
        self.status.lifecycle_status = "failed".to_string();
        self.status.last_error = Some(message.to_string());
        self.status.cancel_path = Some(read_acquisition_path_string(&self.cancel_path));
        self.status.updated_at_unix_ms = GentleEngine::now_unix_ms();
        self.status.finished_at_unix_ms = Some(GentleEngine::now_unix_ms());
        self.write_status();
        remove_read_acquire_activity_lock(&self.lock_path);
    }

    fn finish_cancelled(&mut self, message: &str) {
        self.status.lifecycle_status = "cancelled".to_string();
        self.status.last_error = Some(message.to_string());
        self.status.cancel_path = Some(read_acquisition_path_string(&self.cancel_path));
        self.status.updated_at_unix_ms = GentleEngine::now_unix_ms();
        self.status.finished_at_unix_ms = Some(GentleEngine::now_unix_ms());
        self.write_status();
        remove_read_acquire_activity_lock(&self.lock_path);
        remove_read_acquire_cancel_marker(&self.cancel_path);
    }

    fn finish_success(&self) {
        remove_read_acquire_activity_status(&self.status_path);
        remove_read_acquire_activity_lock(&self.lock_path);
        remove_read_acquire_cancel_marker(&self.cancel_path);
    }

    fn write_status(&self) {
        let _ = write_read_acquire_activity_status(&self.status_path, &self.status);
        let _ = write_read_acquire_activity_status(&self.lock_path, &self.status);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn parse_read_acquisition_manifest_reads_optional_fields() {
        let td = tempdir().expect("tempdir");
        let manifest = td.path().join("reads.tsv");
        fs::write(
            &manifest,
            "sample_id\tsample_name\tsra_accession\tassay_kind\tread_layout\tanalysis_format\tnote\nsample_a\tSample A\tSRRFAKE001\trna\tpaired_end\tfastq\tprimary\n",
        )
        .expect("write read acquisition manifest");

        let rows = parse_read_acquisition_manifest(&manifest.display().to_string())
            .expect("parse manifest");

        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].sample_id, "sample_a");
        assert_eq!(rows[0].sample_name.as_deref(), Some("Sample A"));
        assert_eq!(rows[0].sra_accession, "SRRFAKE001");
        assert_eq!(rows[0].assay_kind.as_deref(), Some("rna"));
        assert_eq!(
            rows[0].read_layout,
            Some(ReadAcquisitionReadLayout::PairedEnd)
        );
        assert_eq!(
            rows[0].analysis_format,
            Some(ReadAcquisitionAnalysisFormat::Fastq)
        );
        assert_eq!(rows[0].note.as_deref(), Some("primary"));
    }

    #[test]
    fn parse_read_acquisition_manifest_rejects_missing_required_columns() {
        let td = tempdir().expect("tempdir");
        let manifest = td.path().join("reads.tsv");
        fs::write(&manifest, "sample_id\tnote\nsample_a\tmissing accession\n")
            .expect("write invalid manifest");

        let error = parse_read_acquisition_manifest(&manifest.display().to_string())
            .expect_err("missing sra_accession should fail");

        assert_eq!(error.code, ErrorCode::InvalidInput);
        assert!(error.message.contains("sra_accession"));
    }

    #[test]
    fn read_acquisition_status_reports_missing_without_external_tools() {
        let td = tempdir().expect("tempdir");
        let manifest = td.path().join("reads.tsv");
        fs::write(
            &manifest,
            "sample_id\tsra_accession\nsample_a\tSRRFAKE001\n",
        )
        .expect("write manifest");
        let engine = GentleEngine::default();

        let report = engine
            .read_acquisition_status(
                &manifest.display().to_string(),
                &td.path().join("cache").display().to_string(),
                &td.path().join("work").display().to_string(),
            )
            .expect("status report");

        assert_eq!(report.lifecycle_status, "missing");
        assert_eq!(report.sample_count, 1);
        assert_eq!(report.missing_count, 1);
        assert_eq!(report.rows[0].resource_key, "read_acquisition:SRRFAKE001");
    }

    #[test]
    fn read_acquisition_cancel_writes_marker_for_running_activity() {
        let td = tempdir().expect("tempdir");
        let cache_dir = td.path().join("cache");
        let work_dir = td.path().join("work");
        let run_dir = read_acquire_run_dir(&work_dir, "SRRFAKE001");
        fs::create_dir_all(&run_dir).expect("create run dir");
        let status_path = read_acquire_activity_status_path(&run_dir);
        let lock_path = read_acquire_activity_lock_path(&run_dir);
        let cancel_path = read_acquire_cancel_path(&run_dir);
        let status = SharedAssetActivityStatus {
            resource_key: read_acquisition_resource_key("SRRFAKE001"),
            display_name: "SRRFAKE001".to_string(),
            status_path: read_acquisition_path_string(&status_path),
            lock_path: Some(read_acquisition_path_string(&lock_path)),
            cancel_path: Some(read_acquisition_path_string(&cancel_path)),
            lifecycle_status: "running".to_string(),
            phase: Some("prefetch".to_string()),
            item: Some("SRRFAKE001".to_string()),
            started_at_unix_ms: GentleEngine::now_unix_ms(),
            updated_at_unix_ms: GentleEngine::now_unix_ms(),
            owner_pid: Some(std::process::id()),
            ..Default::default()
        };
        write_read_acquire_activity_status(&status_path, &status).expect("write status");
        write_read_acquire_activity_status(&lock_path, &status).expect("write lock");
        let engine = GentleEngine::default();

        let report = engine
            .read_acquisition_cancel(
                "SRRFAKE001",
                &cache_dir.display().to_string(),
                &work_dir.display().to_string(),
            )
            .expect("cancel running acquisition");

        assert!(cancel_path.exists(), "cancel marker should be written");
        assert_eq!(report.lifecycle_status, "running");
        assert!(
            report
                .warnings
                .iter()
                .any(|warning| warning.contains("Cancellation requested")),
            "cancel report should tell the user cancellation was requested"
        );
    }

    #[test]
    fn read_acquisition_monitor_rejects_disk_space_drop() {
        let td = tempdir().expect("tempdir");
        if available_disk_bytes(td.path()).is_none() {
            return;
        }
        let cache_dir = td.path().join("cache");
        let work_dir = td.path().join("work");
        let run_dir = read_acquire_run_dir(&work_dir, "SRRFAKE001");
        fs::create_dir_all(&run_dir).expect("create run dir");
        let mut tracker = ReadAcquisitionActivityTracker {
            status_path: read_acquire_activity_status_path(&run_dir),
            lock_path: read_acquire_activity_lock_path(&run_dir),
            cancel_path: read_acquire_cancel_path(&run_dir),
            status: SharedAssetActivityStatus {
                resource_key: read_acquisition_resource_key("SRRFAKE001"),
                display_name: "SRRFAKE001".to_string(),
                status_path: read_acquisition_path_string(&read_acquire_activity_status_path(
                    &run_dir,
                )),
                lock_path: Some(read_acquisition_path_string(
                    &read_acquire_activity_lock_path(&run_dir),
                )),
                cancel_path: Some(read_acquisition_path_string(&read_acquire_cancel_path(
                    &run_dir,
                ))),
                lifecycle_status: "running".to_string(),
                phase: Some("prefetch".to_string()),
                item: Some("SRRFAKE001".to_string()),
                started_at_unix_ms: GentleEngine::now_unix_ms(),
                updated_at_unix_ms: GentleEngine::now_unix_ms(),
                owner_pid: Some(std::process::id()),
                ..Default::default()
            },
        };
        let mut progress_called = false;

        let error = update_read_acquisition_command_monitor(
            &mut tracker,
            "prefetch",
            &cache_dir,
            &work_dir,
            std::slice::from_ref(&run_dir),
            Some(u64::MAX),
            &mut |_progress| {
                progress_called = true;
                true
            },
        )
        .expect_err("huge min-free-gb should trip monitored disk threshold");

        assert!(!progress_called, "disk failure should stop before progress");
        assert!(error.message.contains("free disk dropped"));
        assert!(tracker.status.monitored_free_bytes.is_some());
        assert_eq!(tracker.status.minimum_free_bytes, Some(u64::MAX));
    }

    #[test]
    fn read_acquisition_inspect_infers_paired_fastq_outputs() {
        let td = tempdir().expect("tempdir");
        let work_dir = td.path().join("work");
        let fastq_dir = work_dir.join("SRRFAKE001/fastq");
        fs::create_dir_all(&fastq_dir).expect("create FASTQ dir");
        fs::write(
            fastq_dir.join("SRRFAKE001_1.fastq"),
            "@read/1\nACGT\n+\n!!!!\n",
        )
        .expect("write R1");
        fs::write(
            fastq_dir.join("SRRFAKE001_2.fastq"),
            "@read/2\nTGCA\n+\n!!!!\n",
        )
        .expect("write R2");
        let engine = GentleEngine::default();

        let report = engine
            .read_acquisition_inspect(
                "SRRFAKE001",
                &td.path().join("cache").display().to_string(),
                &work_dir.display().to_string(),
            )
            .expect("inspect report");

        assert_eq!(report.lifecycle_status, "ready");
        assert_eq!(report.ready_count, 1);
        assert_eq!(
            report.rows[0].analysis_format,
            ReadAcquisitionAnalysisFormat::Fastq
        );
        assert_eq!(
            report.rows[0].read_layout,
            ReadAcquisitionReadLayout::PairedEnd
        );
        assert_eq!(report.rows[0].output_paths.len(), 2);
    }
}

fn read_acquisition_resource_key(sra_accession: &str) -> String {
    format!("read_acquisition:{}", sra_accession.trim())
}

fn read_acquisition_path_string(path: &Path) -> String {
    fs::canonicalize(path)
        .unwrap_or_else(|_| path.to_path_buf())
        .display()
        .to_string()
}

fn sanitize_read_acquisition_token(value: &str) -> String {
    let token = value
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || matches!(ch, '-' | '_' | '.') {
                ch
            } else {
                '_'
            }
        })
        .collect::<String>();
    let trimmed = token.trim_matches('_');
    if trimmed.is_empty() {
        "read_run".to_string()
    } else {
        trimmed.to_string()
    }
}

fn read_acquire_run_dir(work_dir: &Path, sra_accession: &str) -> PathBuf {
    work_dir.join(sanitize_read_acquisition_token(sra_accession))
}

fn read_acquire_sra_dir(cache_dir: &Path, sra_accession: &str) -> PathBuf {
    cache_dir.join(sanitize_read_acquisition_token(sra_accession))
}

fn read_acquire_sra_path(cache_dir: &Path, sra_accession: &str) -> PathBuf {
    read_acquire_sra_dir(cache_dir, sra_accession).join(format!("{sra_accession}.sra"))
}

fn read_acquire_activity_status_path(run_dir: &Path) -> PathBuf {
    run_dir.join(READ_ACQUIRE_ACTIVITY_STATUS_FILE)
}

fn read_acquire_activity_lock_path(run_dir: &Path) -> PathBuf {
    run_dir.join(READ_ACQUIRE_ACTIVITY_LOCK_FILE)
}

fn read_acquire_cancel_path(run_dir: &Path) -> PathBuf {
    run_dir.join(READ_ACQUIRE_CANCEL_FILE)
}

fn load_read_acquire_activity_status(path: &Path) -> Result<SharedAssetActivityStatus, String> {
    let text = fs::read_to_string(path).map_err(|e| {
        format!(
            "Could not read read-acquisition activity status '{}': {e}",
            path.display()
        )
    })?;
    serde_json::from_str(&text).map_err(|e| {
        format!(
            "Could not parse read-acquisition activity status '{}': {e}",
            path.display()
        )
    })
}

fn write_read_acquire_activity_status(
    path: &Path,
    status: &SharedAssetActivityStatus,
) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Could not create read-acquisition activity directory '{}': {e}",
                parent.display()
            )
        })?;
    }
    let text = serde_json::to_string_pretty(status)
        .map_err(|e| format!("Could not serialize read-acquisition activity status: {e}"))?;
    fs::write(path, text).map_err(|e| {
        format!(
            "Could not write read-acquisition activity status '{}': {e}",
            path.display()
        )
    })
}

fn create_read_acquire_activity_lock(
    path: &Path,
    status: &SharedAssetActivityStatus,
) -> Result<bool, String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| {
            format!(
                "Could not create read-acquisition lock directory '{}': {e}",
                parent.display()
            )
        })?;
    }
    match OpenOptions::new().write(true).create_new(true).open(path) {
        Ok(mut file) => {
            let text = serde_json::to_string_pretty(status)
                .map_err(|e| format!("Could not serialize read-acquisition lock status: {e}"))?;
            file.write_all(text.as_bytes()).map_err(|e| {
                format!(
                    "Could not write read-acquisition lock '{}': {e}",
                    path.display()
                )
            })?;
            Ok(true)
        }
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => Ok(false),
        Err(e) => Err(format!(
            "Could not create read-acquisition lock '{}': {e}",
            path.display()
        )),
    }
}

fn remove_read_acquire_activity_status(path: &Path) {
    if path.exists() {
        let _ = fs::remove_file(path);
    }
}

fn remove_read_acquire_activity_lock(path: &Path) {
    if path.exists() {
        let _ = fs::remove_file(path);
    }
}

fn remove_read_acquire_cancel_marker(path: &Path) {
    if path.exists() {
        let _ = fs::remove_file(path);
    }
}

fn read_acquire_activity_is_stale(status: &SharedAssetActivityStatus, now: u128) -> bool {
    status.lifecycle_status == "running"
        && now.saturating_sub(status.updated_at_unix_ms) > READ_ACQUIRE_ACTIVITY_STALE_AFTER_MS
}

#[cfg(unix)]
fn read_acquire_owner_pid_is_running(pid: u32) -> Option<bool> {
    if pid == 0 {
        return None;
    }
    Command::new("ps")
        .arg("-p")
        .arg(pid.to_string())
        .arg("-o")
        .arg("pid=")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .ok()
        .map(|status| status.success())
}

#[cfg(not(unix))]
fn read_acquire_owner_pid_is_running(_pid: u32) -> Option<bool> {
    None
}

fn read_acquire_owner_pid_stale_reason(status: &SharedAssetActivityStatus) -> Option<String> {
    if status.lifecycle_status != "running" {
        return None;
    }
    let pid = status.owner_pid?;
    match read_acquire_owner_pid_is_running(pid) {
        Some(false) => Some(format!(
            "Read-acquisition owner_pid {pid} is no longer running and is treated as stale"
        )),
        Some(true) | None => None,
    }
}

fn mark_read_acquire_activity_stale(
    status_path: &Path,
    lock_path: &Path,
    mut status: SharedAssetActivityStatus,
    reason: &str,
) -> SharedAssetActivityStatus {
    status.lifecycle_status = "stale".to_string();
    status.last_error = Some(reason.to_string());
    status.updated_at_unix_ms = GentleEngine::now_unix_ms();
    let _ = write_read_acquire_activity_status(status_path, &status);
    remove_read_acquire_activity_lock(lock_path);
    status
}

fn inspect_read_acquire_activity_status_paths(
    status_path: &Path,
    lock_path: &Path,
) -> Result<Option<SharedAssetActivityStatus>, String> {
    let status = if status_path.exists() {
        load_read_acquire_activity_status(status_path)?
    } else if lock_path.exists() {
        load_read_acquire_activity_status(lock_path)?
    } else {
        return Ok(None);
    };
    if status.lifecycle_status == "running" {
        if !lock_path.exists() {
            return Ok(Some(mark_read_acquire_activity_stale(
                status_path,
                lock_path,
                status,
                "Read-acquisition activity lost its active lock and is treated as stale",
            )));
        }
        if let Some(reason) = read_acquire_owner_pid_stale_reason(&status) {
            return Ok(Some(mark_read_acquire_activity_stale(
                status_path,
                lock_path,
                status,
                &reason,
            )));
        }
        if read_acquire_activity_is_stale(&status, GentleEngine::now_unix_ms()) {
            return Ok(Some(mark_read_acquire_activity_stale(
                status_path,
                lock_path,
                status,
                "Read-acquisition heartbeat is stale and can be retried safely",
            )));
        }
    } else {
        remove_read_acquire_activity_lock(lock_path);
    }
    Ok(Some(status))
}

fn read_acquisition_io_error(path: &Path, action: &str, error: std::io::Error) -> EngineError {
    EngineError {
        code: ErrorCode::Io,
        message: format!("Could not {action} '{}': {error}", path.display()),
        cause_chain: vec![],
    }
}

fn read_acquisition_invalid_input(message: impl Into<String>) -> EngineError {
    EngineError {
        code: ErrorCode::InvalidInput,
        message: message.into(),
        cause_chain: vec![],
    }
}

fn parse_read_acquisition_analysis_format(
    raw: &str,
) -> Result<ReadAcquisitionAnalysisFormat, EngineError> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "fasta" | "fa" => Ok(ReadAcquisitionAnalysisFormat::Fasta),
        "fastq" | "fq" => Ok(ReadAcquisitionAnalysisFormat::Fastq),
        other => Err(read_acquisition_invalid_input(format!(
            "Unknown read acquisition analysis_format '{other}' (expected fasta|fastq)"
        ))),
    }
}

fn parse_read_acquisition_layout(raw: &str) -> Result<ReadAcquisitionReadLayout, EngineError> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "single_end" | "single-end" | "single" | "se" => Ok(ReadAcquisitionReadLayout::SingleEnd),
        "paired_end" | "paired-end" | "paired" | "pe" => Ok(ReadAcquisitionReadLayout::PairedEnd),
        "split_spot" | "split-spot" | "splitspot" => Ok(ReadAcquisitionReadLayout::SplitSpot),
        other => Err(read_acquisition_invalid_input(format!(
            "Unknown read acquisition read_layout '{other}' (expected single_end|paired_end|split_spot)"
        ))),
    }
}

fn parse_read_acquisition_manifest(
    manifest_path: &str,
) -> Result<Vec<ReadAcquisitionManifestRow>, EngineError> {
    let manifest_path = manifest_path.trim();
    if manifest_path.is_empty() {
        return Err(read_acquisition_invalid_input(
            "reads acquire requires a non-empty manifest path",
        ));
    }
    let text = fs::read_to_string(manifest_path).map_err(|e| EngineError {
        code: ErrorCode::Io,
        message: format!("Could not read read-acquisition manifest '{manifest_path}': {e}"),

        cause_chain: vec![],
    })?;
    let mut meaningful_lines = text.lines().enumerate().filter(|(_, line)| {
        let trimmed = line.trim();
        !trimmed.is_empty() && !trimmed.starts_with('#')
    });
    let Some((header_index, header_line)) = meaningful_lines.next() else {
        return Err(read_acquisition_invalid_input(format!(
            "Read-acquisition manifest '{manifest_path}' is empty"
        )));
    };
    let headers = header_line
        .trim_start_matches('\u{feff}')
        .split('\t')
        .map(|value| value.trim().to_ascii_lowercase())
        .collect::<Vec<_>>();
    let header_index_of = |name: &str| headers.iter().position(|header| header == name);
    let sample_id_idx = header_index_of("sample_id").ok_or_else(|| {
        read_acquisition_invalid_input("Read-acquisition manifest requires 'sample_id'")
    })?;
    let sra_accession_idx = header_index_of("sra_accession").ok_or_else(|| {
        read_acquisition_invalid_input("Read-acquisition manifest requires 'sra_accession'")
    })?;
    let sample_name_idx = header_index_of("sample_name");
    let assay_kind_idx = header_index_of("assay_kind");
    let read_layout_idx = header_index_of("read_layout");
    let analysis_format_idx = header_index_of("analysis_format");
    let note_idx = header_index_of("note");
    let mut rows = Vec::<ReadAcquisitionManifestRow>::new();
    for (line_index, line) in meaningful_lines {
        let columns = line.split('\t').collect::<Vec<_>>();
        let field = |idx: Option<usize>| -> Option<String> {
            idx.and_then(|idx| columns.get(idx))
                .map(|value| value.trim())
                .filter(|value| !value.is_empty())
                .map(ToString::to_string)
        };
        let sample_id = columns
            .get(sample_id_idx)
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .ok_or_else(|| {
                read_acquisition_invalid_input(format!(
                    "Read-acquisition manifest line {} has empty sample_id",
                    line_index + 1
                ))
            })?
            .to_string();
        let sra_accession = columns
            .get(sra_accession_idx)
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .ok_or_else(|| {
                read_acquisition_invalid_input(format!(
                    "Read-acquisition manifest line {} sample '{}' has empty sra_accession",
                    line_index + 1,
                    sample_id
                ))
            })?
            .to_string();
        let read_layout = field(read_layout_idx)
            .as_deref()
            .map(parse_read_acquisition_layout)
            .transpose()?;
        let analysis_format = field(analysis_format_idx)
            .as_deref()
            .map(parse_read_acquisition_analysis_format)
            .transpose()?;
        rows.push(ReadAcquisitionManifestRow {
            row_number: line_index.saturating_add(1).max(header_index + 2),
            sample_id,
            sra_accession,
            sample_name: field(sample_name_idx),
            assay_kind: field(assay_kind_idx),
            read_layout,
            analysis_format,
            note: field(note_idx),
        });
    }
    if rows.is_empty() {
        return Err(read_acquisition_invalid_input(format!(
            "Read-acquisition manifest '{manifest_path}' contains no sample rows"
        )));
    }
    Ok(rows)
}

fn effective_read_acquisition_rows(
    rows: Vec<ReadAcquisitionManifestRow>,
    default_format: ReadAcquisitionAnalysisFormat,
    default_layout: ReadAcquisitionReadLayout,
) -> Vec<ReadAcquisitionEffectiveRow> {
    rows.into_iter()
        .map(|row| ReadAcquisitionEffectiveRow {
            read_layout: row.read_layout.unwrap_or(default_layout),
            analysis_format: row.analysis_format.unwrap_or(default_format),
            manifest: row,
        })
        .collect()
}

fn expected_fastq_paths(
    dump_dir: &Path,
    accession: &str,
    layout: ReadAcquisitionReadLayout,
) -> Vec<(String, PathBuf)> {
    match layout {
        ReadAcquisitionReadLayout::PairedEnd => vec![
            (
                "r1".to_string(),
                dump_dir.join(format!("{accession}_1.fastq")),
            ),
            (
                "r2".to_string(),
                dump_dir.join(format!("{accession}_2.fastq")),
            ),
        ],
        ReadAcquisitionReadLayout::SingleEnd | ReadAcquisitionReadLayout::SplitSpot => {
            vec![(
                "single".to_string(),
                dump_dir.join(format!("{accession}.fastq")),
            )]
        }
    }
}

fn expected_fasta_paths(
    fasta_dir: &Path,
    accession: &str,
    layout: ReadAcquisitionReadLayout,
) -> Vec<(String, PathBuf)> {
    match layout {
        ReadAcquisitionReadLayout::PairedEnd => vec![
            (
                "r1".to_string(),
                fasta_dir.join(format!("{accession}_1.fasta")),
            ),
            (
                "r2".to_string(),
                fasta_dir.join(format!("{accession}_2.fasta")),
            ),
        ],
        ReadAcquisitionReadLayout::SingleEnd | ReadAcquisitionReadLayout::SplitSpot => {
            vec![(
                "single".to_string(),
                fasta_dir.join(format!("{accession}.fasta")),
            )]
        }
    }
}

fn shell_join(program: &str, args: &[String]) -> String {
    std::iter::once(program.to_string())
        .chain(args.iter().cloned())
        .map(|part| {
            if part
                .chars()
                .all(|ch| ch.is_ascii_alphanumeric() || "-_./:=+".contains(ch))
            {
                part
            } else {
                format!("'{}'", part.replace('\'', "'\\''"))
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn directory_size_bytes(path: &Path) -> u64 {
    let Ok(metadata) = fs::metadata(path) else {
        return 0;
    };
    if metadata.is_file() {
        return metadata.len();
    }
    let Ok(read_dir) = fs::read_dir(path) else {
        return 0;
    };
    read_dir
        .flatten()
        .map(|entry| directory_size_bytes(&entry.path()))
        .sum()
}

fn min_free_bytes(min_free_gb: Option<u64>) -> Option<u64> {
    min_free_gb.and_then(|value| {
        if value == 0 {
            None
        } else {
            Some(value.saturating_mul(1024 * 1024 * 1024))
        }
    })
}

fn monitored_free_bytes(cache_dir: &Path, work_dir: &Path) -> Option<u64> {
    [
        available_disk_bytes(cache_dir),
        available_disk_bytes(work_dir),
    ]
    .into_iter()
    .flatten()
    .min()
}

fn produced_progress_bytes(paths: &[PathBuf]) -> u64 {
    paths.iter().map(|path| directory_size_bytes(path)).sum()
}

fn spawn_read_acquisition_log_copy<R>(
    mut reader: R,
    path: PathBuf,
) -> std::thread::JoinHandle<Result<(), String>>
where
    R: Read + Send + 'static,
{
    std::thread::spawn(move || {
        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .map_err(|e| {
                format!(
                    "Could not open read-acquisition stream log '{}': {e}",
                    path.display()
                )
            })?;
        std::io::copy(&mut reader, &mut file).map_err(|e| {
            format!(
                "Could not write read-acquisition stream log '{}': {e}",
                path.display()
            )
        })?;
        Ok(())
    })
}

fn join_read_acquisition_log_copy(
    handle: std::thread::JoinHandle<Result<(), String>>,
) -> Result<(), EngineError> {
    handle
        .join()
        .map_err(|_| EngineError {
            code: ErrorCode::Io,
            message: "Read-acquisition log-copy worker panicked".to_string(),

            cause_chain: vec![],
        })?
        .map_err(|message| EngineError {
            code: ErrorCode::Io,
            message,

            cause_chain: vec![],
        })
}

fn terminate_read_acquisition_child(child: &mut std::process::Child) {
    #[cfg(unix)]
    {
        let process_group = format!("-{}", child.id());
        let _ = Command::new("kill")
            .arg("-TERM")
            .arg(&process_group)
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status();
        std::thread::sleep(std::time::Duration::from_millis(250));
        if child.try_wait().ok().flatten().is_none() {
            let _ = Command::new("kill")
                .arg("-KILL")
                .arg(&process_group)
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status();
        }
    }
    let _ = child.kill();
    let _ = child.wait();
}

fn update_read_acquisition_command_monitor(
    tracker: &mut ReadAcquisitionActivityTracker,
    phase: &str,
    cache_dir: &Path,
    work_dir: &Path,
    progress_paths: &[PathBuf],
    min_free_gb: Option<u64>,
    on_progress: &mut dyn FnMut(OperationProgress) -> bool,
) -> Result<bool, EngineError> {
    let produced_bytes = produced_progress_bytes(progress_paths);
    let monitored_free = monitored_free_bytes(cache_dir, work_dir);
    let minimum_free = min_free_bytes(min_free_gb);
    tracker.update_monitor(phase, produced_bytes, monitored_free, minimum_free);
    if let (Some(available), Some(required)) = (monitored_free, minimum_free) {
        if available < required {
            return Err(read_acquisition_invalid_input(format!(
                "Read acquisition free disk dropped below the configured minimum during phase '{phase}': required {:.2} GB, available {:.2} GB",
                required as f64 / 1024.0 / 1024.0 / 1024.0,
                available as f64 / 1024.0 / 1024.0 / 1024.0
            )));
        }
    }
    Ok(on_progress(OperationProgress::ReadAcquisition(
        tracker.status.clone(),
    )))
}

fn run_read_acquisition_command(
    log_dir: &Path,
    phase: &str,
    program: &str,
    args: &[String],
    tracker: &mut ReadAcquisitionActivityTracker,
    cache_dir: &Path,
    work_dir: &Path,
    progress_paths: &[PathBuf],
    min_free_gb: Option<u64>,
    on_progress: &mut dyn FnMut(OperationProgress) -> bool,
) -> Result<ReadAcquisitionCommandProvenance, EngineError> {
    fs::create_dir_all(log_dir)
        .map_err(|e| read_acquisition_io_error(log_dir, "create read-acquisition log dir", e))?;
    let stdout_log_path = log_dir.join(format!("{phase}.stdout.log"));
    let stderr_log_path = log_dir.join(format!("{phase}.stderr.log"));
    let command_log_path = log_dir.join(format!("{phase}.command.txt"));
    let command = shell_join(program, args);
    fs::write(&command_log_path, format!("{command}\n")).map_err(|e| {
        read_acquisition_io_error(&command_log_path, "write read-acquisition command log", e)
    })?;
    let mut command_builder = Command::new(program);
    command_builder
        .args(args)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    #[cfg(unix)]
    {
        command_builder.process_group(0);
    }
    let mut child = command_builder.spawn().map_err(|e| EngineError {
        code: ErrorCode::Io,
        message: format!("Could not run SRA Toolkit command for phase '{phase}': {command}: {e}"),

        cause_chain: vec![],
    })?;
    let stdout_handle = child
        .stdout
        .take()
        .map(|stdout| spawn_read_acquisition_log_copy(stdout, stdout_log_path.clone()));
    let stderr_handle = child
        .stderr
        .take()
        .map(|stderr| spawn_read_acquisition_log_copy(stderr, stderr_log_path.clone()));
    let mut last_monitor = std::time::Instant::now()
        .checked_sub(std::time::Duration::from_millis(
            READ_ACQUIRE_MONITOR_INTERVAL_MS,
        ))
        .unwrap_or_else(std::time::Instant::now);
    let status = loop {
        if tracker.cancel_path.exists() {
            terminate_read_acquisition_child(&mut child);
            break Err(EngineError {
                code: ErrorCode::Io,
                message: format!("Read acquisition cancelled during phase '{phase}'"),

                cause_chain: vec![],
            });
        }
        if last_monitor.elapsed()
            >= std::time::Duration::from_millis(READ_ACQUIRE_MONITOR_INTERVAL_MS)
        {
            last_monitor = std::time::Instant::now();
            match update_read_acquisition_command_monitor(
                tracker,
                phase,
                cache_dir,
                work_dir,
                progress_paths,
                min_free_gb,
                on_progress,
            ) {
                Ok(true) => {}
                Ok(false) => {
                    terminate_read_acquisition_child(&mut child);
                    break Err(EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Read acquisition cancelled by progress callback during phase '{phase}'"
                        ),

                        cause_chain: vec![],
                    });
                }
                Err(error) => {
                    terminate_read_acquisition_child(&mut child);
                    break Err(error);
                }
            }
        }
        match child.try_wait().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not poll SRA Toolkit command '{command}': {e}"),

            cause_chain: vec![],
        })? {
            Some(status) => break Ok(status),
            None => std::thread::sleep(std::time::Duration::from_millis(100)),
        }
    };
    if let Some(handle) = stdout_handle {
        join_read_acquisition_log_copy(handle)?;
    }
    if let Some(handle) = stderr_handle {
        join_read_acquisition_log_copy(handle)?;
    }
    let status = status?;
    let provenance = ReadAcquisitionCommandProvenance {
        phase: phase.to_string(),
        command,
        exit_code: status.code(),
        stdout_log_path: read_acquisition_path_string(&stdout_log_path),
        stderr_log_path: read_acquisition_path_string(&stderr_log_path),
    };
    if !status.success() {
        return Err(EngineError {
            code: ErrorCode::Io,
            message: format!(
                "SRA Toolkit command failed in phase '{}' with exit_code={:?}; see '{}'",
                phase,
                status.code(),
                stderr_log_path.display()
            ),

            cause_chain: vec![],
        });
    }
    Ok(provenance)
}

fn find_sra_file(root: &Path, accession: &str) -> Option<PathBuf> {
    let mut stack = vec![root.to_path_buf()];
    while let Some(path) = stack.pop() {
        let Ok(read_dir) = fs::read_dir(&path) else {
            continue;
        };
        for entry in read_dir.flatten() {
            let entry_path = entry.path();
            if entry_path.is_dir() {
                stack.push(entry_path);
            } else if entry_path
                .file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name == format!("{accession}.sra"))
            {
                return Some(entry_path);
            }
        }
    }
    None
}

fn ensure_sra_at_expected_path(sra_dir: &Path, accession: &str) -> Result<PathBuf, EngineError> {
    let expected = sra_dir.join(format!("{accession}.sra"));
    if expected.exists() {
        return Ok(expected);
    }
    let Some(found) = find_sra_file(sra_dir, accession) else {
        return Err(EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "prefetch did not create an SRA file for '{}' under '{}'",
                accession,
                sra_dir.display()
            ),

            cause_chain: vec![],
        });
    };
    fs::copy(&found, &expected).map_err(|e| {
        read_acquisition_io_error(
            &expected,
            "copy prefetched SRA into canonical cache path",
            e,
        )
    })?;
    Ok(expected)
}

fn sha1_file(path: &Path) -> Result<String, EngineError> {
    let mut file =
        File::open(path).map_err(|e| read_acquisition_io_error(path, "open output for SHA1", e))?;
    let mut hasher = Sha1::new();
    let mut buffer = [0u8; 8192];
    loop {
        let bytes = file
            .read(&mut buffer)
            .map_err(|e| read_acquisition_io_error(path, "read output for SHA1", e))?;
        if bytes == 0 {
            break;
        }
        hasher.update(&buffer[..bytes]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn stats_from_lengths(lengths: &[usize]) -> ReadAcquisitionReadLengthStats {
    let total_bases = lengths.iter().sum::<usize>();
    let read_count = lengths.len();
    ReadAcquisitionReadLengthStats {
        read_count,
        total_bases,
        min_length_bp: lengths.iter().copied().min().unwrap_or(0),
        max_length_bp: lengths.iter().copied().max().unwrap_or(0),
        mean_length_bp: if read_count == 0 {
            0.0
        } else {
            total_bases as f64 / read_count as f64
        },
    }
}

fn read_fastq_length_stats(path: &Path) -> Result<ReadAcquisitionReadLengthStats, EngineError> {
    let file =
        File::open(path).map_err(|e| read_acquisition_io_error(path, "open FASTQ for stats", e))?;
    let mut lines = BufReader::new(file).lines();
    let mut lengths = Vec::<usize>::new();
    loop {
        let Some(header) = lines.next() else {
            break;
        };
        let _header =
            header.map_err(|e| read_acquisition_io_error(path, "read FASTQ header", e))?;
        let Some(seq) = lines.next() else {
            break;
        };
        let seq = seq.map_err(|e| read_acquisition_io_error(path, "read FASTQ sequence", e))?;
        let _plus = lines
            .next()
            .transpose()
            .map_err(|e| read_acquisition_io_error(path, "read FASTQ plus line", e))?;
        let _qual = lines
            .next()
            .transpose()
            .map_err(|e| read_acquisition_io_error(path, "read FASTQ quality line", e))?;
        lengths.push(seq.trim().len());
    }
    Ok(stats_from_lengths(&lengths))
}

fn read_fasta_length_stats(path: &Path) -> Result<ReadAcquisitionReadLengthStats, EngineError> {
    let file =
        File::open(path).map_err(|e| read_acquisition_io_error(path, "open FASTA for stats", e))?;
    let reader = BufReader::new(file);
    let mut lengths = Vec::<usize>::new();
    let mut current = 0usize;
    let mut saw_record = false;
    for line in reader.lines() {
        let line = line.map_err(|e| read_acquisition_io_error(path, "read FASTA", e))?;
        if line.starts_with('>') {
            if saw_record {
                lengths.push(current);
            }
            saw_record = true;
            current = 0;
        } else {
            current = current.saturating_add(line.trim().len());
        }
    }
    if saw_record {
        lengths.push(current);
    }
    Ok(stats_from_lengths(&lengths))
}

fn convert_fastq_to_fasta(fastq_path: &Path, fasta_path: &Path) -> Result<(), EngineError> {
    if let Some(parent) = fasta_path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| read_acquisition_io_error(parent, "create FASTA output dir", e))?;
    }
    let input = File::open(fastq_path)
        .map_err(|e| read_acquisition_io_error(fastq_path, "open FASTQ for FASTA conversion", e))?;
    let mut lines = BufReader::new(input).lines();
    let mut output =
        BufWriter::new(File::create(fasta_path).map_err(|e| {
            read_acquisition_io_error(fasta_path, "create FASTA conversion output", e)
        })?);
    loop {
        let Some(header) = lines.next() else {
            break;
        };
        let header =
            header.map_err(|e| read_acquisition_io_error(fastq_path, "read FASTQ header", e))?;
        let Some(seq) = lines.next() else {
            break;
        };
        let seq =
            seq.map_err(|e| read_acquisition_io_error(fastq_path, "read FASTQ sequence", e))?;
        let _plus = lines
            .next()
            .transpose()
            .map_err(|e| read_acquisition_io_error(fastq_path, "read FASTQ plus line", e))?;
        let _qual = lines
            .next()
            .transpose()
            .map_err(|e| read_acquisition_io_error(fastq_path, "read FASTQ quality line", e))?;
        let id = header
            .trim_start_matches('@')
            .split_whitespace()
            .next()
            .unwrap_or("read");
        writeln!(output, ">{id}")
            .map_err(|e| read_acquisition_io_error(fasta_path, "write FASTA header", e))?;
        writeln!(output, "{}", seq.trim())
            .map_err(|e| read_acquisition_io_error(fasta_path, "write FASTA sequence", e))?;
    }
    output
        .flush()
        .map_err(|e| read_acquisition_io_error(fasta_path, "flush FASTA output", e))
}

fn output_path_report(
    role: &str,
    path: &Path,
    analysis_format: ReadAcquisitionAnalysisFormat,
) -> Result<ReadAcquisitionOutputPath, EngineError> {
    let metadata = fs::metadata(path)
        .map_err(|e| read_acquisition_io_error(path, "stat prepared read output", e))?;
    let read_length_stats = match analysis_format {
        ReadAcquisitionAnalysisFormat::Fasta => Some(read_fasta_length_stats(path)?),
        ReadAcquisitionAnalysisFormat::Fastq => Some(read_fastq_length_stats(path)?),
    };
    Ok(ReadAcquisitionOutputPath {
        role: role.to_string(),
        path: read_acquisition_path_string(path),
        analysis_format,
        file_size_bytes: metadata.len(),
        checksum_sha1: sha1_file(path)?,
        read_length_stats,
    })
}

fn available_disk_bytes(path: &Path) -> Option<u64> {
    fs::create_dir_all(path).ok()?;
    let output = Command::new("df").arg("-Pk").arg(path).output().ok()?;
    if !output.status.success() {
        return None;
    }
    let text = String::from_utf8_lossy(&output.stdout);
    let line = text.lines().nth(1)?;
    let fields = line.split_whitespace().collect::<Vec<_>>();
    let available_kb = fields.get(3)?.parse::<u64>().ok()?;
    Some(available_kb.saturating_mul(1024))
}

fn ensure_min_free_gb(
    path: &Path,
    min_free_gb: Option<u64>,
) -> Result<Option<String>, EngineError> {
    let Some(min_free_gb) = min_free_gb else {
        return Ok(None);
    };
    if min_free_gb == 0 {
        return Ok(None);
    }
    let required = min_free_gb.saturating_mul(1024 * 1024 * 1024);
    match available_disk_bytes(path) {
        Some(available) if available < required => Err(read_acquisition_invalid_input(format!(
            "Read acquisition requires at least {min_free_gb} GB free at '{}', but only {:.2} GB are available",
            path.display(),
            available as f64 / 1024.0 / 1024.0 / 1024.0
        ))),
        Some(_) => Ok(None),
        None => Ok(Some(format!(
            "Could not determine free disk space at '{}'; continuing without min-free-gb enforcement",
            path.display()
        ))),
    }
}

fn read_acquisition_report_status(rows: &[ReadAcquisitionRunReport]) -> String {
    if rows.iter().any(|row| row.lifecycle_status == "running") {
        "running".to_string()
    } else if rows.iter().any(|row| row.lifecycle_status == "failed") {
        "failed".to_string()
    } else if rows.iter().any(|row| row.lifecycle_status == "cancelled") {
        "cancelled".to_string()
    } else if rows.iter().any(|row| row.lifecycle_status == "stale") {
        "stale".to_string()
    } else if rows.iter().all(|row| row.lifecycle_status == "ready") && !rows.is_empty() {
        "ready".to_string()
    } else {
        "missing".to_string()
    }
}

fn build_read_acquisition_report(
    manifest_path: Option<String>,
    cache_dir: &Path,
    work_dir: &Path,
    rows: Vec<ReadAcquisitionRunReport>,
    warnings: Vec<String>,
) -> ReadAcquisitionReport {
    let ready_count = rows
        .iter()
        .filter(|row| row.lifecycle_status == "ready")
        .count();
    let running_count = rows
        .iter()
        .filter(|row| row.lifecycle_status == "running")
        .count();
    let failed_count = rows
        .iter()
        .filter(|row| row.lifecycle_status == "failed")
        .count();
    let cancelled_count = rows
        .iter()
        .filter(|row| row.lifecycle_status == "cancelled")
        .count();
    let stale_count = rows
        .iter()
        .filter(|row| row.lifecycle_status == "stale")
        .count();
    let missing_count = rows
        .iter()
        .filter(|row| row.lifecycle_status == "missing")
        .count();
    let lifecycle_status = read_acquisition_report_status(&rows);
    ReadAcquisitionReport {
        schema: READ_ACQUISITION_REPORT_SCHEMA.to_string(),
        manifest_path,
        cache_dir: read_acquisition_path_string(cache_dir),
        work_dir: read_acquisition_path_string(work_dir),
        generated_at_unix_ms: GentleEngine::now_unix_ms(),
        lifecycle_status,
        sample_count: rows.len(),
        ready_count,
        running_count,
        failed_count,
        cancelled_count,
        stale_count,
        missing_count,
        rows,
        warnings,
    }
}

impl GentleEngine {
    fn read_acquisition_status_for_row(
        row: &ReadAcquisitionEffectiveRow,
        cache_dir: &Path,
        work_dir: &Path,
    ) -> Result<ReadAcquisitionRunReport, EngineError> {
        let accession = row.manifest.sra_accession.trim();
        let run_dir = read_acquire_run_dir(work_dir, accession);
        let sra_path = read_acquire_sra_path(cache_dir, accession);
        let current_activity = inspect_read_acquire_activity_status_paths(
            &read_acquire_activity_status_path(&run_dir),
            &read_acquire_activity_lock_path(&run_dir),
        )
        .map_err(|message| EngineError {
            code: ErrorCode::Io,
            message,

            cause_chain: vec![],
        })?;
        let final_paths = match row.analysis_format {
            ReadAcquisitionAnalysisFormat::Fastq => {
                expected_fastq_paths(&run_dir.join("fastq"), accession, row.read_layout)
            }
            ReadAcquisitionAnalysisFormat::Fasta => {
                expected_fasta_paths(&run_dir.join("fasta"), accession, row.read_layout)
            }
        };
        let mut warnings = Vec::<String>::new();
        let mut outputs = Vec::<ReadAcquisitionOutputPath>::new();
        let mut all_outputs_ready = true;
        for (role, path) in final_paths {
            if path.exists() {
                outputs.push(output_path_report(&role, &path, row.analysis_format)?);
            } else {
                all_outputs_ready = false;
            }
        }
        if !sra_path.exists() && outputs.is_empty() {
            warnings.push(format!(
                "No cached SRA file was found at '{}'",
                sra_path.display()
            ));
        }
        let lifecycle_status = if let Some(activity) = current_activity.as_ref() {
            activity.lifecycle_status.clone()
        } else if all_outputs_ready && !outputs.is_empty() {
            "ready".to_string()
        } else {
            "missing".to_string()
        };
        Ok(ReadAcquisitionRunReport {
            sample_id: row.manifest.sample_id.clone(),
            sra_accession: accession.to_string(),
            sample_name: row.manifest.sample_name.clone(),
            assay_kind: row.manifest.assay_kind.clone(),
            note: row.manifest.note.clone(),
            resource_key: read_acquisition_resource_key(accession),
            lifecycle_status,
            read_layout: row.read_layout,
            analysis_format: row.analysis_format,
            sra_path: read_acquisition_path_string(&sra_path),
            run_dir: read_acquisition_path_string(&run_dir),
            output_paths: outputs,
            current_activity,
            warnings,
            ..Default::default()
        })
    }

    fn prepare_read_acquisition_row(
        &self,
        row: &ReadAcquisitionEffectiveRow,
        cache_dir: &Path,
        work_dir: &Path,
        threads: Option<usize>,
        max_size: Option<&str>,
        min_free_gb: Option<u64>,
        drop_intermediate_fastq: bool,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<ReadAcquisitionRunReport, EngineError> {
        let accession = row.manifest.sra_accession.trim();
        let run_dir = read_acquire_run_dir(work_dir, accession);
        let sra_dir = read_acquire_sra_dir(cache_dir, accession);
        let log_dir = run_dir.join("logs");
        fs::create_dir_all(&sra_dir)
            .map_err(|e| read_acquisition_io_error(&sra_dir, "create SRA cache dir", e))?;
        fs::create_dir_all(&run_dir)
            .map_err(|e| read_acquisition_io_error(&run_dir, "create read run dir", e))?;
        let mut warnings = Vec::<String>::new();
        if let Some(warning) = ensure_min_free_gb(cache_dir, min_free_gb)? {
            warnings.push(warning);
        }
        if let Some(warning) = ensure_min_free_gb(work_dir, min_free_gb)? {
            warnings.push(warning);
        }
        let mut tracker =
            match ReadAcquisitionActivityTracker::start(&run_dir, accession).map_err(|message| {
                EngineError {
                    code: ErrorCode::Io,
                    message,

                    cause_chain: vec![],
                }
            })? {
                ReadAcquisitionActivityStart::Acquired(tracker) => tracker,
                ReadAcquisitionActivityStart::Running => {
                    return Self::read_acquisition_status_for_row(row, cache_dir, work_dir);
                }
            };
        let already_ready = Self::read_acquisition_status_for_row(row, cache_dir, work_dir)?;
        if already_ready.lifecycle_status == "ready" {
            tracker.finish_success();
            return Ok(already_ready);
        }

        let mut command_provenance = Vec::<ReadAcquisitionCommandProvenance>::new();
        let prepared = (|| -> Result<ReadAcquisitionRunReport, EngineError> {
            tracker.update_phase("prefetch");
            let mut prefetch_args = vec![
                accession.to_string(),
                "-O".to_string(),
                read_acquisition_path_string(&sra_dir),
            ];
            if let Some(max_size) = max_size.filter(|value| !value.trim().is_empty()) {
                prefetch_args.push("--max-size".to_string());
                prefetch_args.push(max_size.to_string());
            }
            command_provenance.push(run_read_acquisition_command(
                &log_dir,
                "prefetch",
                "prefetch",
                &prefetch_args,
                &mut tracker,
                cache_dir,
                work_dir,
                &[sra_dir.clone(), run_dir.clone()],
                min_free_gb,
                on_progress,
            )?);
            let sra_path = ensure_sra_at_expected_path(&sra_dir, accession)?;

            tracker.update_phase("validate_sra");
            command_provenance.push(run_read_acquisition_command(
                &log_dir,
                "validate_sra",
                "vdb-validate",
                &[read_acquisition_path_string(&sra_path)],
                &mut tracker,
                cache_dir,
                work_dir,
                &[sra_dir.clone(), run_dir.clone()],
                min_free_gb,
                on_progress,
            )?);

            tracker.update_phase("dump_fastq");
            let dump_dir = run_dir.join("fastq");
            fs::create_dir_all(&dump_dir)
                .map_err(|e| read_acquisition_io_error(&dump_dir, "create FASTQ output dir", e))?;
            let mut dump_args = vec![
                read_acquisition_path_string(&sra_path),
                "--outdir".to_string(),
                read_acquisition_path_string(&dump_dir),
            ];
            if let Some(threads) = threads.filter(|value| *value > 0) {
                dump_args.push("--threads".to_string());
                dump_args.push(threads.to_string());
            }
            match row.read_layout {
                ReadAcquisitionReadLayout::PairedEnd => {
                    dump_args.push("--split-files".to_string());
                }
                ReadAcquisitionReadLayout::SplitSpot => {
                    dump_args.push("--split-spot".to_string());
                }
                ReadAcquisitionReadLayout::SingleEnd => {}
            }
            command_provenance.push(run_read_acquisition_command(
                &log_dir,
                "dump_fastq",
                "fasterq-dump",
                &dump_args,
                &mut tracker,
                cache_dir,
                work_dir,
                &[sra_dir.clone(), run_dir.clone()],
                min_free_gb,
                on_progress,
            )?);
            let fastq_paths = expected_fastq_paths(&dump_dir, accession, row.read_layout);
            for (_, path) in &fastq_paths {
                if !path.exists() {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "fasterq-dump did not create expected FASTQ output '{}'",
                            path.display()
                        ),

                        cause_chain: vec![],
                    });
                }
            }

            let output_paths = match row.analysis_format {
                ReadAcquisitionAnalysisFormat::Fastq => {
                    tracker.update_phase("verify_output");
                    fastq_paths
                        .iter()
                        .map(|(role, path)| {
                            output_path_report(role, path, ReadAcquisitionAnalysisFormat::Fastq)
                        })
                        .collect::<Result<Vec<_>, _>>()?
                }
                ReadAcquisitionAnalysisFormat::Fasta => {
                    tracker.update_phase("convert_fasta");
                    let fasta_dir = run_dir.join("fasta");
                    let fasta_paths = expected_fasta_paths(&fasta_dir, accession, row.read_layout);
                    for ((_, fastq_path), (_, fasta_path)) in
                        fastq_paths.iter().zip(fasta_paths.iter())
                    {
                        convert_fastq_to_fasta(fastq_path, fasta_path)?;
                    }
                    tracker.update_phase("verify_output");
                    let outputs = fasta_paths
                        .iter()
                        .map(|(role, path)| {
                            output_path_report(role, path, ReadAcquisitionAnalysisFormat::Fasta)
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    if drop_intermediate_fastq {
                        for (_, fastq_path) in &fastq_paths {
                            if let Err(error) = fs::remove_file(fastq_path) {
                                warnings.push(format!(
                                    "Could not remove intermediate FASTQ '{}': {error}",
                                    fastq_path.display()
                                ));
                            }
                        }
                    }
                    outputs
                }
            };
            tracker.update_phase("complete");
            Ok(ReadAcquisitionRunReport {
                sample_id: row.manifest.sample_id.clone(),
                sra_accession: accession.to_string(),
                sample_name: row.manifest.sample_name.clone(),
                assay_kind: row.manifest.assay_kind.clone(),
                note: row.manifest.note.clone(),
                resource_key: read_acquisition_resource_key(accession),
                lifecycle_status: "ready".to_string(),
                read_layout: row.read_layout,
                analysis_format: row.analysis_format,
                sra_path: read_acquisition_path_string(&sra_path),
                run_dir: read_acquisition_path_string(&run_dir),
                output_paths,
                command_provenance,
                warnings,
                ..Default::default()
            })
        })();

        match prepared {
            Ok(report) => {
                tracker.finish_success();
                Ok(report)
            }
            Err(error) => {
                let cancelled = error.message.to_ascii_lowercase().contains("cancel");
                if cancelled {
                    tracker.finish_cancelled(&error.message);
                } else {
                    tracker.finish_failure(&error.message);
                }
                let mut failed = Self::read_acquisition_status_for_row(row, cache_dir, work_dir)
                    .unwrap_or_else(|_| ReadAcquisitionRunReport {
                        sample_id: row.manifest.sample_id.clone(),
                        sra_accession: accession.to_string(),
                        resource_key: read_acquisition_resource_key(accession),
                        lifecycle_status: if cancelled {
                            "cancelled".to_string()
                        } else {
                            "failed".to_string()
                        },
                        read_layout: row.read_layout,
                        analysis_format: row.analysis_format,
                        ..Default::default()
                    });
                failed.lifecycle_status = if cancelled {
                    "cancelled".to_string()
                } else {
                    "failed".to_string()
                };
                failed.error = Some(error.message.clone());
                Err(EngineError {
                    code: error.code,
                    message: failed.error.clone().unwrap_or(error.message),

                    cause_chain: vec![],
                })
            }
        }
    }

    pub fn read_acquisition_status(
        &self,
        manifest_path: &str,
        cache_dir: &str,
        work_dir: &str,
    ) -> Result<ReadAcquisitionReport, EngineError> {
        let cache_dir = PathBuf::from(cache_dir);
        let work_dir = PathBuf::from(work_dir);
        let rows = effective_read_acquisition_rows(
            parse_read_acquisition_manifest(manifest_path)?,
            ReadAcquisitionAnalysisFormat::Fasta,
            ReadAcquisitionReadLayout::SingleEnd,
        );
        let mut reports = Vec::<ReadAcquisitionRunReport>::new();
        for row in rows {
            reports.push(Self::read_acquisition_status_for_row(
                &row, &cache_dir, &work_dir,
            )?);
        }
        Ok(build_read_acquisition_report(
            Some(manifest_path.to_string()),
            &cache_dir,
            &work_dir,
            reports,
            vec![],
        ))
    }

    pub fn read_acquisition_inspect(
        &self,
        sra_accession: &str,
        cache_dir: &str,
        work_dir: &str,
    ) -> Result<ReadAcquisitionReport, EngineError> {
        let accession = sra_accession.trim();
        if accession.is_empty() {
            return Err(read_acquisition_invalid_input(
                "reads acquire inspect requires a non-empty RUN_ACCESSION",
            ));
        }
        let cache_dir = PathBuf::from(cache_dir);
        let work_dir = PathBuf::from(work_dir);
        let manifest = ReadAcquisitionManifestRow {
            row_number: 1,
            sample_id: accession.to_string(),
            sra_accession: accession.to_string(),
            ..Default::default()
        };
        let candidates = [
            (
                ReadAcquisitionAnalysisFormat::Fasta,
                ReadAcquisitionReadLayout::SingleEnd,
            ),
            (
                ReadAcquisitionAnalysisFormat::Fasta,
                ReadAcquisitionReadLayout::PairedEnd,
            ),
            (
                ReadAcquisitionAnalysisFormat::Fastq,
                ReadAcquisitionReadLayout::SingleEnd,
            ),
            (
                ReadAcquisitionAnalysisFormat::Fastq,
                ReadAcquisitionReadLayout::PairedEnd,
            ),
            (
                ReadAcquisitionAnalysisFormat::Fasta,
                ReadAcquisitionReadLayout::SplitSpot,
            ),
            (
                ReadAcquisitionAnalysisFormat::Fastq,
                ReadAcquisitionReadLayout::SplitSpot,
            ),
        ];
        let mut reports = Vec::<ReadAcquisitionRunReport>::new();
        for (analysis_format, read_layout) in candidates {
            let row = ReadAcquisitionEffectiveRow {
                manifest: manifest.clone(),
                read_layout,
                analysis_format,
            };
            reports.push(Self::read_acquisition_status_for_row(
                &row, &cache_dir, &work_dir,
            )?);
        }
        let report = reports
            .iter()
            .find(|report| report.lifecycle_status == "ready")
            .or_else(|| {
                reports.iter().find(|report| {
                    matches!(
                        report.lifecycle_status.as_str(),
                        "running" | "failed" | "cancelled" | "stale"
                    )
                })
            })
            .cloned()
            .unwrap_or_else(|| reports.remove(0));
        Ok(build_read_acquisition_report(
            None,
            &cache_dir,
            &work_dir,
            vec![report],
            vec![],
        ))
    }

    pub fn read_acquisition_cancel(
        &self,
        sra_accession: &str,
        cache_dir: &str,
        work_dir: &str,
    ) -> Result<ReadAcquisitionReport, EngineError> {
        let accession = sra_accession.trim();
        if accession.is_empty() {
            return Err(read_acquisition_invalid_input(
                "reads acquire cancel requires a non-empty RUN_ACCESSION",
            ));
        }
        let run_dir = read_acquire_run_dir(&PathBuf::from(work_dir), accession);
        let status = inspect_read_acquire_activity_status_paths(
            &read_acquire_activity_status_path(&run_dir),
            &read_acquire_activity_lock_path(&run_dir),
        )
        .map_err(|message| EngineError {
            code: ErrorCode::Io,
            message,

            cause_chain: vec![],
        })?;
        let mut report = self.read_acquisition_inspect(accession, cache_dir, work_dir)?;
        if status
            .as_ref()
            .is_some_and(|status| status.lifecycle_status == "running")
        {
            fs::create_dir_all(&run_dir).map_err(|e| {
                read_acquisition_io_error(&run_dir, "create read-acquisition run dir", e)
            })?;
            let cancel_path = read_acquire_cancel_path(&run_dir);
            fs::write(
                &cancel_path,
                format!(
                    "cancel requested at {}\n",
                    read_acquisition_path_string(&cancel_path)
                ),
            )
            .map_err(|e| {
                read_acquisition_io_error(&cancel_path, "write read-acquisition cancel marker", e)
            })?;
            report.warnings.push(format!(
                "Cancellation requested for SRA run '{}' via '{}'",
                accession,
                cancel_path.display()
            ));
        } else {
            report.warnings.push(format!(
                "No running read acquisition was found for SRA run '{}'",
                accession
            ));
        }
        Ok(report)
    }

    pub fn read_acquisition_prepare(
        &self,
        manifest_path: &str,
        cache_dir: &str,
        work_dir: &str,
        analysis_format: ReadAcquisitionAnalysisFormat,
        read_layout: ReadAcquisitionReadLayout,
        threads: Option<usize>,
        max_size: Option<&str>,
        min_free_gb: Option<u64>,
        drop_intermediate_fastq: bool,
        continue_on_error: bool,
    ) -> Result<ReadAcquisitionReport, EngineError> {
        let mut on_progress = |_progress: OperationProgress| true;
        self.read_acquisition_prepare_with_progress(
            manifest_path,
            cache_dir,
            work_dir,
            analysis_format,
            read_layout,
            threads,
            max_size,
            min_free_gb,
            drop_intermediate_fastq,
            continue_on_error,
            &mut on_progress,
        )
    }

    pub fn read_acquisition_prepare_with_progress(
        &self,
        manifest_path: &str,
        cache_dir: &str,
        work_dir: &str,
        analysis_format: ReadAcquisitionAnalysisFormat,
        read_layout: ReadAcquisitionReadLayout,
        threads: Option<usize>,
        max_size: Option<&str>,
        min_free_gb: Option<u64>,
        drop_intermediate_fastq: bool,
        continue_on_error: bool,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<ReadAcquisitionReport, EngineError> {
        let cache_dir = PathBuf::from(cache_dir);
        let work_dir = PathBuf::from(work_dir);
        fs::create_dir_all(&cache_dir)
            .map_err(|e| read_acquisition_io_error(&cache_dir, "create read cache dir", e))?;
        fs::create_dir_all(&work_dir)
            .map_err(|e| read_acquisition_io_error(&work_dir, "create read work dir", e))?;
        let rows = effective_read_acquisition_rows(
            parse_read_acquisition_manifest(manifest_path)?,
            analysis_format,
            read_layout,
        );
        let mut reports_by_accession = HashMap::<
            (
                String,
                ReadAcquisitionAnalysisFormat,
                ReadAcquisitionReadLayout,
            ),
            ReadAcquisitionRunReport,
        >::new();
        let mut reports = Vec::<ReadAcquisitionRunReport>::new();
        let mut warnings = Vec::<String>::new();
        for row in rows {
            let key = (
                row.manifest.sra_accession.to_ascii_uppercase(),
                row.analysis_format,
                row.read_layout,
            );
            let mut report = if let Some(existing) = reports_by_accession.get(&key) {
                let mut duplicate = existing.clone();
                duplicate.sample_id = row.manifest.sample_id.clone();
                duplicate.sample_name = row.manifest.sample_name.clone();
                duplicate.assay_kind = row.manifest.assay_kind.clone();
                duplicate.note = row.manifest.note.clone();
                duplicate.warnings.push(format!(
                    "Reused prepared output from duplicate accession '{}'",
                    row.manifest.sra_accession
                ));
                duplicate
            } else {
                match self.prepare_read_acquisition_row(
                    &row,
                    &cache_dir,
                    &work_dir,
                    threads,
                    max_size,
                    min_free_gb,
                    drop_intermediate_fastq,
                    on_progress,
                ) {
                    Ok(report) => {
                        reports_by_accession.insert(key, report.clone());
                        report
                    }
                    Err(error) => {
                        if !continue_on_error {
                            return Err(error);
                        }
                        let cancelled = error.message.to_ascii_lowercase().contains("cancel");
                        let lifecycle_status = if cancelled { "cancelled" } else { "failed" };
                        warnings.push(format!(
                            "sample '{}' accession '{}' {} but continue-on-error is enabled: {}",
                            row.manifest.sample_id,
                            row.manifest.sra_accession,
                            if cancelled { "was cancelled" } else { "failed" },
                            error.message
                        ));
                        ReadAcquisitionRunReport {
                            sample_id: row.manifest.sample_id.clone(),
                            sra_accession: row.manifest.sra_accession.clone(),
                            sample_name: row.manifest.sample_name.clone(),
                            assay_kind: row.manifest.assay_kind.clone(),
                            note: row.manifest.note.clone(),
                            resource_key: read_acquisition_resource_key(
                                &row.manifest.sra_accession,
                            ),
                            lifecycle_status: lifecycle_status.to_string(),
                            read_layout: row.read_layout,
                            analysis_format: row.analysis_format,
                            sra_path: read_acquisition_path_string(&read_acquire_sra_path(
                                &cache_dir,
                                &row.manifest.sra_accession,
                            )),
                            run_dir: read_acquisition_path_string(&read_acquire_run_dir(
                                &work_dir,
                                &row.manifest.sra_accession,
                            )),
                            error: Some(error.message),
                            ..Default::default()
                        }
                    }
                }
            };
            if report.lifecycle_status == "ready" && report.output_paths.is_empty() {
                report.lifecycle_status = "failed".to_string();
                report.error = Some("Prepared read acquisition returned no output paths".into());
            }
            reports.push(report);
        }
        Ok(build_read_acquisition_report(
            Some(manifest_path.to_string()),
            &cache_dir,
            &work_dir,
            reports,
            warnings,
        ))
    }

    pub(crate) fn read_acquisition_prepare_single_accession(
        &self,
        sample_id: &str,
        sra_accession: &str,
        cache_dir: &Path,
        work_dir: &Path,
        analysis_format: ReadAcquisitionAnalysisFormat,
        read_layout: ReadAcquisitionReadLayout,
        drop_intermediate_fastq: bool,
    ) -> Result<ReadAcquisitionRunReport, EngineError> {
        let mut on_progress = |_progress: OperationProgress| true;
        self.read_acquisition_prepare_single_accession_with_progress(
            sample_id,
            sra_accession,
            cache_dir,
            work_dir,
            analysis_format,
            read_layout,
            drop_intermediate_fastq,
            &mut on_progress,
        )
    }

    pub(crate) fn read_acquisition_prepare_single_accession_with_progress(
        &self,
        sample_id: &str,
        sra_accession: &str,
        cache_dir: &Path,
        work_dir: &Path,
        analysis_format: ReadAcquisitionAnalysisFormat,
        read_layout: ReadAcquisitionReadLayout,
        drop_intermediate_fastq: bool,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<ReadAcquisitionRunReport, EngineError> {
        let row = ReadAcquisitionEffectiveRow {
            manifest: ReadAcquisitionManifestRow {
                row_number: 1,
                sample_id: sample_id.to_string(),
                sra_accession: sra_accession.to_string(),
                ..Default::default()
            },
            read_layout,
            analysis_format,
        };
        self.prepare_read_acquisition_row(
            &row,
            cache_dir,
            work_dir,
            None,
            None,
            None,
            drop_intermediate_fastq,
            on_progress,
        )
    }
}
