//! Separate payload evidence, explicitly bound to Phase 1 bytes and caller-supplied references.

use super::*;
use ring::digest::{Context, SHA256};
use std::collections::BTreeSet;
use std::io::{BufRead, BufReader};

#[cfg(test)]
#[path = "primer_payload/alignment_tests.rs"]
mod alignment_tests;
#[path = "primer_payload/evidence.rs"]
mod evidence;
#[path = "primer_payload/inspect.rs"]
mod inspect;
#[path = "primer_payload/observed.rs"]
mod observed;
#[path = "primer_payload/junction_assessment.rs"]
mod junction_assessment;
#[cfg(test)]
#[path = "primer_payload/observed_tests.rs"]
mod observed_tests;
#[path = "primer_payload/paf.rs"]
mod paf;
#[path = "primer_payload/reference.rs"]
mod reference;
#[path = "primer_payload/summarize.rs"]
mod summarize;

#[path = "primer_payload/binding.rs"]
mod binding;
#[cfg(test)]
#[path = "primer_payload/binding_tests.rs"]
mod binding_tests;
#[path = "primer_payload/prepare.rs"]
mod prepare;
#[cfg(test)]
#[path = "primer_payload/tests.rs"]
mod tests;

const INTERRUPTED: &str = "primer_payload_interrupted";

#[path = "primer_payload/annotation_prepare.rs"]
mod annotation_prepare;

fn hash<T: Serialize>(value: &T) -> Result<String, EngineError> {
    Ok(sha256_prefixed_bytes(&encode(value)?))
}

fn is_hash(value: &str) -> bool {
    value
        .strip_prefix("sha256:")
        .is_some_and(|s| valid_digest(s, 64))
}

fn budgets(b: &PrimerPayloadBudgets) -> Result<(), EngineError> {
    if !(1..=1_000_000_000_000).contains(&b.max_input_bytes)
        || !(1..=4_000_000_000).contains(&b.max_output_bytes)
        || !(1..=1_000_000).contains(&b.max_records)
        || !(1..=1_000_000).contains(&b.max_queries)
        || !(1..=2_000_000).contains(&b.max_query_bases)
        || !(128..=8_000_000).contains(&b.max_line_bytes)
        || !(1..=8192).contains(&b.max_alignments_per_query)
        || !(1..=1_000_000).contains(&b.max_cigar_operations)
        || !(1..=4_000_000).contains(&b.max_annotation_rows)
        || !(1..=128_000_000).contains(&b.max_reference_context_bases)
        || !(1..=100_000_000).contains(&b.max_chain_comparisons)
    {
        return Err(bad("primer_payload_budget_bounds"));
    }
    Ok(())
}

fn regular(path: &Path) -> Result<fs::Metadata, EngineError> {
    let meta = fs::symlink_metadata(path).map_err(|e| bad(format!("primer_payload_file: {e}")))?;
    if !meta.is_file() || meta.file_type().is_symlink() {
        return Err(bad("primer_payload_regular_file_required"));
    }
    Ok(meta)
}

/// Hash exactly the bytes consumed by a bounded parser; EOF and path stability are mandatory.
struct VerifiedReader {
    file: fs::File,
    path: PathBuf,
    before: fs::Metadata,
    expected: PrimerPayloadFile,
    hasher: Context,
    bytes: u64,
}

impl VerifiedReader {
    fn open(expected: &PrimerPayloadFile, limit: u64) -> Result<Self, EngineError> {
        if !is_hash(&expected.sha256) || expected.bytes > limit {
            return Err(bad("primer_payload_file_binding_or_budget"));
        }
        let path = PathBuf::from(&expected.path);
        let before = regular(&path)?;
        if before.len() != expected.bytes {
            return Err(bad("primer_payload_file_size"));
        }
        let file = fs::File::open(&path).map_err(|e| bad(e.to_string()))?;
        let opened = file.metadata().map_err(|e| bad(e.to_string()))?;
        if opened.len() != before.len() || opened.modified().ok() != before.modified().ok() {
            return Err(bad("primer_payload_file_changed"));
        }
        Ok(Self {
            file,
            path,
            before,
            expected: expected.clone(),
            hasher: Context::new(&SHA256),
            bytes: 0,
        })
    }

    fn finish(self) -> Result<(), EngineError> {
        let after = regular(&self.path)?;
        let opened = self.file.metadata().map_err(|e| bad(e.to_string()))?;
        if self.bytes != self.expected.bytes
            || after.len() != self.before.len()
            || after.modified().ok() != self.before.modified().ok()
            || opened.len() != after.len()
            || opened.modified().ok() != after.modified().ok()
            || format!("sha256:{}", input::hex(self.hasher.finish().as_ref()))
                != self.expected.sha256
        {
            return Err(bad("primer_payload_file_hash_or_change"));
        }
        Ok(())
    }
}

impl Read for VerifiedReader {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        let n = self.file.read(bytes)?;
        self.bytes = self
            .bytes
            .checked_add(n as u64)
            .ok_or_else(|| std::io::Error::other("primer_payload_byte_overflow"))?;
        if self.bytes > self.expected.bytes {
            return Err(std::io::Error::other("primer_payload_file_grew"));
        }
        self.hasher.update(&bytes[..n]);
        Ok(n)
    }
}

fn read_json<T: serde::de::DeserializeOwned>(
    file: &PrimerPayloadFile,
    b: &PrimerPayloadBudgets,
) -> Result<T, EngineError> {
    let mut reader = BufReader::new(VerifiedReader::open(
        file,
        b.max_input_bytes.min(128_000_000),
    )?);
    let value = serde_json::from_reader(&mut reader)
        .map_err(|e| bad(format!("primer_payload_json: {e}")))?;
    reader.into_inner().finish()?;
    Ok(value)
}

fn artifact_file(root: &Path, a: &RnaReadPanelArtifact) -> Result<PrimerPayloadFile, EngineError> {
    let meta = fs::symlink_metadata(root).map_err(|e| bad(e.to_string()))?;
    if !meta.is_dir() || meta.file_type().is_symlink() {
        return Err(bad("primer_payload_artifact_directory"));
    }
    let mut current = root.to_path_buf();
    for part in Path::new(&a.path).components() {
        let std::path::Component::Normal(name) = part else {
            return Err(bad("primer_payload_artifact_path"));
        };
        current.push(name);
        let meta = fs::symlink_metadata(&current).map_err(|e| bad(e.to_string()))?;
        if meta.file_type().is_symlink() {
            return Err(bad("primer_payload_artifact_symlink"));
        }
    }
    Ok(PrimerPayloadFile {
        path: current.display().to_string(),
        sha256: a.sha256.clone(),
        bytes: a.bytes,
    })
}

fn verified_bytes(
    file: &PrimerPayloadFile,
    b: &PrimerPayloadBudgets,
) -> Result<Vec<u8>, EngineError> {
    let mut reader = VerifiedReader::open(file, b.max_input_bytes.min(128_000_000))?;
    let mut bytes = Vec::new();
    reader
        .read_to_end(&mut bytes)
        .map_err(|e| bad(e.to_string()))?;
    reader.finish()?;
    Ok(bytes)
}

fn line(reader: &mut impl BufRead, limit: usize) -> Result<Option<String>, EngineError> {
    let mut bytes = Vec::new();
    if reader
        .take(limit as u64 + 1)
        .read_until(b'\n', &mut bytes)
        .map_err(|e| bad(e.to_string()))?
        == 0
    {
        return Ok(None);
    }
    if bytes.len() > limit {
        return Err(bad("primer_payload_line_budget"));
    }
    if bytes.last() == Some(&b'\n') {
        bytes.pop();
    }
    if bytes.last() == Some(&b'\r') {
        bytes.pop();
    }
    String::from_utf8(bytes)
        .map(Some)
        .map_err(|_| bad("primer_payload_text_utf8"))
}

fn union(mut spans: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    spans.sort_unstable();
    let mut result: Vec<(usize, usize)> = Vec::new();
    for (start, end) in spans {
        if start >= end {
            continue;
        }
        if let Some(last) = result.last_mut()
            && start <= last.1
        {
            last.1 = last.1.max(end);
        } else {
            result.push((start, end));
        }
    }
    result
}

fn span_bases(spans: &[(usize, usize)]) -> usize {
    spans.iter().map(|(a, b)| b - a).sum()
}

struct Control<'a> {
    marker: Option<String>,
    identity: String,
    progress: &'a mut dyn FnMut(OperationProgress) -> bool,
    gate: control::PollGate,
    completed: usize,
    total: usize,
    stopped: bool,
    error: Option<String>,
}

impl<'a> Control<'a> {
    fn new(
        marker: Option<String>,
        identity: String,
        total: usize,
        progress: &'a mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<Self, EngineError> {
        let stopped =
            control::marker_requested(marker.as_deref(), std::slice::from_ref(&identity))?;
        Ok(Self {
            marker,
            identity,
            progress,
            gate: control::PollGate::default(),
            completed: 0,
            total,
            stopped,
            error: None,
        })
    }
    fn stop(&mut self, force: bool) -> bool {
        if self.stopped {
            return true;
        }
        if !self.gate.due(std::time::Instant::now(), force) {
            return false;
        }
        match control::marker_requested(
            self.marker.as_deref(),
            std::slice::from_ref(&self.identity),
        ) {
            Ok(true) => self.stopped = true,
            Ok(false) => {
                self.stopped = !(self.progress)(OperationProgress::Workflow {
                    completed: self.completed,
                    total: self.total,
                });
                // Progress/control observers may change a marker, including at finalization.
                match control::marker_requested(
                    self.marker.as_deref(),
                    std::slice::from_ref(&self.identity),
                ) {
                    Ok(requested) => self.stopped |= requested,
                    Err(error) => {
                        self.error = Some(error.message);
                        self.stopped = true;
                    }
                }
            }
            Err(error) => {
                self.error = Some(error.message);
                self.stopped = true;
            }
        }
        self.stopped
    }
    fn check(&mut self) -> Result<(), EngineError> {
        if self.stop(false) {
            Err(bad(INTERRUPTED))
        } else {
            Ok(())
        }
    }
    fn status(&self) -> &str {
        if self.error.is_some() {
            "control_error"
        } else if self.stopped {
            "cancelled"
        } else {
            "completed"
        }
    }
}

fn staging(output: &str) -> Result<tempfile::TempDir, EngineError> {
    if output.trim().is_empty()
        || Path::new(output).exists()
        || fs::symlink_metadata(output).is_ok()
    {
        return Err(bad("primer_payload_output_exists_or_empty"));
    }
    let parent = Path::new(output)
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    tempfile::Builder::new()
        .prefix(".gentle-primer-payload-")
        .tempdir_in(parent)
        .map_err(|e| bad(e.to_string()))
}

fn publish<T: Serialize>(
    stage: tempfile::TempDir,
    output: &str,
    report: &T,
    b: &PrimerPayloadBudgets,
) -> Result<(), EngineError> {
    let bytes = serde_json::to_vec_pretty(report).map_err(|e| bad(e.to_string()))?;
    fn size(path: &Path) -> Result<u64, EngineError> {
        let mut total = 0u64;
        for entry in fs::read_dir(path).map_err(|e| bad(e.to_string()))? {
            let entry = entry.map_err(|e| bad(e.to_string()))?;
            let meta = fs::symlink_metadata(entry.path()).map_err(|e| bad(e.to_string()))?;
            let n = if meta.is_dir() {
                size(&entry.path())?
            } else if meta.is_file() {
                meta.len()
            } else {
                return Err(bad("primer_payload_staging_file"));
            };
            total = total
                .checked_add(n)
                .ok_or_else(|| bad("primer_payload_output_budget"))?;
        }
        Ok(total)
    }
    if bytes.len() > 128_000_000
        || size(stage.path())?
            .checked_add(bytes.len() as u64)
            .is_none_or(|n| n > b.max_output_bytes)
    {
        return Err(bad("primer_payload_report_budget"));
    }
    let out = Path::new(output);
    fs::create_dir(out).map_err(|e| bad(format!("primer_payload_output_reserve: {e}")))?;
    if let Err(e) = fs::rename(stage.path(), out.join("artifacts")) {
        let _ = fs::remove_dir(out);
        return Err(bad(format!("primer_payload_publish: {e}")));
    }
    let mut file = tempfile::NamedTempFile::new_in(out).map_err(|e| bad(e.to_string()))?;
    file.write_all(&bytes).map_err(|e| bad(e.to_string()))?;
    file.as_file().sync_all().map_err(|e| bad(e.to_string()))?;
    file.persist_noclobber(out.join("report.json"))
        .map_err(|e| bad(e.to_string()))?;
    Ok(())
}

fn nonclaims() -> Vec<String> {
    vec![
        "Content/receipt consistency is not authenticated producer provenance or an independent alignment search.".into(),
        "Read support is not original molecule count, expression, exclusive isoform identity or whole-reference uniqueness.".into(),
        "Primer-payload edges are not splice junctions; foreign partners do not establish fusion, trans-splicing or an experimental mechanism.".into(),
        "Primer masking affects new query copies and support accounting only, never original FASTQ or Phase 1 evidence.".into(),
    ]
}
