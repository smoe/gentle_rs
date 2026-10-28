//! Preparation-only path admission and owned, no-clobber/report-last publication.

use super::*;

struct Location {
    resolved: PathBuf,
    identity: Option<same_file::Handle>,
}

fn io(e: std::io::Error) -> EngineError {
    bad(format!("observed_annotation_io: {e}"))
}

/// Validate the original spelling before components/canonicalize can erase it.
fn raw_path(raw: &str) -> Result<&Path, EngineError> {
    if !text(raw, 32768) {
        return Err(bad("observed_annotation_raw_path"));
    }
    let p = Path::new(raw);
    let prefix = match p.components().next() {
        Some(std::path::Component::Prefix(v)) => v.as_os_str().len(),
        _ => 0,
    };
    let sep = |c| c == '/' || (cfg!(windows) && c == '\\');
    let mut body = &raw[prefix..];
    if body.starts_with(sep) {
        body = &body[1..];
    }
    if body.is_empty()
        || body.split(sep).any(|s| {
            s.is_empty()
                || s == "."
                || s == ".."
                || (cfg!(windows)
                    && (s.contains(':')
                        || s.ends_with(['.', ' '])
                        || matches!(
                            s.split('.')
                                .next()
                                .unwrap_or("")
                                .to_ascii_uppercase()
                                .as_str(),
                            "CON"
                                | "PRN"
                                | "AUX"
                                | "NUL"
                                | "COM1"
                                | "COM2"
                                | "COM3"
                                | "COM4"
                                | "COM5"
                                | "COM6"
                                | "COM7"
                                | "COM8"
                                | "COM9"
                                | "LPT1"
                                | "LPT2"
                                | "LPT3"
                                | "LPT4"
                                | "LPT5"
                                | "LPT6"
                                | "LPT7"
                                | "LPT8"
                                | "LPT9"
                        )))
        })
    {
        return Err(bad("observed_annotation_raw_path"));
    }
    Ok(p)
}

fn location(raw: &str, input: bool) -> Result<Location, EngineError> {
    location_kind(raw, input, false)
}

fn location_kind(raw: &str, input: bool, directory: bool) -> Result<Location, EngineError> {
    let path = raw_path(raw)?;
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let parent = fs::canonicalize(parent).map_err(io)?;
    if !parent.is_dir() {
        return Err(bad("observed_annotation_parent_directory"));
    }
    let leaf = path
        .file_name()
        .ok_or_else(|| bad("observed_annotation_path_leaf"))?;
    let resolved = parent.join(leaf);
    let identity = match fs::symlink_metadata(&resolved) {
        Ok(meta) => {
            if !(meta.is_file() || (directory && meta.is_dir())) || meta.file_type().is_symlink() {
                return Err(bad("observed_annotation_regular_file"));
            }
            Some(same_file::Handle::from_path(&resolved).map_err(io)?)
        }
        Err(e) if !input && e.kind() == std::io::ErrorKind::NotFound => None,
        Err(e) => return Err(io(e)),
    };
    Ok(Location { resolved, identity })
}

fn equivalent(a: &Path, b: &Path) -> bool {
    a == b
        || ((cfg!(windows) || cfg!(target_os = "macos"))
            && a.to_string_lossy().to_lowercase() == b.to_string_lossy().to_lowercase())
}

fn aliases(a: &Location, b: &Location) -> bool {
    equivalent(&a.resolved, &b.resolved)
        || a.identity
            .as_ref()
            .zip(b.identity.as_ref())
            .is_some_and(|(a, b)| a == b)
}

pub(super) struct Paths {
    inputs: [(String, Location); 2],
    output_raw: String,
    output: Location,
    marker: Option<(String, Location)>,
}

impl Paths {
    pub(super) fn new(r: &PrepareObservedTranscriptAnnotationRequest) -> Result<Self, EngineError> {
        Self::admit(r, true)
    }

    /// Control validates aliases without re-opening source bytes or requiring them to remain.
    pub(super) fn for_cancel(
        r: &PrepareObservedTranscriptAnnotationRequest,
    ) -> Result<Self, EngineError> {
        Self::admit(r, false)
    }

    fn admit(
        r: &PrepareObservedTranscriptAnnotationRequest,
        prepare: bool,
    ) -> Result<Self, EngineError> {
        let a = location(&r.gtf.path, prepare)?;
        let b = location(&r.dictionary.path, prepare)?;
        let output = location_kind(&r.output_dir, false, !prepare)?;
        if (prepare && output.identity.is_some())
            || aliases(&a, &b)
            || aliases(&output, &a)
            || aliases(&output, &b)
        {
            return Err(bad("observed_annotation_path_collision_or_existing_output"));
        }
        let marker = r
            .cancel_path
            .as_ref()
            .map(|p| -> Result<_, EngineError> { Ok((p.clone(), location(p, false)?)) })
            .transpose()?;
        if let Some((_, m)) = &marker {
            // Case-fold only as a conservative additional refusal on common
            // case-insensitive platforms; existing aliases use file identities.
            let contains = m.resolved.starts_with(&output.resolved)
                || ((cfg!(windows) || cfg!(target_os = "macos"))
                    && m.resolved.to_string_lossy().to_lowercase().starts_with(
                        &(output.resolved.to_string_lossy().to_lowercase()
                            + std::path::MAIN_SEPARATOR_STR),
                    ));
            if aliases(m, &a) || aliases(m, &b) || aliases(m, &output) || contains {
                return Err(bad("observed_annotation_control_collision"));
            }
        }
        Ok(Self {
            inputs: [(r.gtf.path.clone(), a), (r.dictionary.path.clone(), b)],
            output_raw: r.output_dir.clone(),
            output,
            marker,
        })
    }

    pub(super) fn recheck_cancel(&self) -> Result<(), EngineError> {
        let marker = self
            .marker
            .as_ref()
            .ok_or_else(|| bad("observed_annotation_cancel_path_required"))?;
        let output = location_kind(&self.output_raw, false, true)?;
        let m = location(&marker.0, false)?;
        if output.resolved != self.output.resolved
            || m.resolved != marker.1.resolved
            || aliases(&m, &output)
        {
            return Err(bad("observed_annotation_control_path_changed"));
        }
        for (raw, _) in &self.inputs {
            if aliases(&m, &location(raw, false)?) {
                return Err(bad("observed_annotation_control_collision"));
            }
        }
        Ok(())
    }

    pub(super) fn unchanged(&self, input: &str) -> Result<(), EngineError> {
        let original = self
            .inputs
            .iter()
            .find(|(raw, _)| raw == input)
            .ok_or_else(|| bad("observed_annotation_input_path"))?;
        let current = location(input, true)?;
        if !equivalent(&original.1.resolved, &current.resolved)
            || original.1.identity != current.identity
        {
            return Err(bad("observed_annotation_input_changed"));
        }
        Ok(())
    }

    pub(super) fn reader(
        &self,
        f: &PrimerPayloadFile,
        limit: u64,
    ) -> Result<VerifiedReader, EngineError> {
        self.unchanged(&f.path)?;
        let reader = VerifiedReader::open(f, limit)?;
        let handle =
            same_file::Handle::from_file(reader.file.try_clone().map_err(io)?).map_err(io)?;
        let bound = &self
            .inputs
            .iter()
            .find(|(raw, _)| raw == &f.path)
            .unwrap()
            .1;
        if bound.identity.as_ref() != Some(&handle) {
            return Err(bad("observed_annotation_input_changed"));
        }
        Ok(reader)
    }

    fn recheck_output(&self) -> Result<(), EngineError> {
        let now = location(&self.output_raw, false)?;
        if now.identity.is_some() || now.resolved != self.output.resolved {
            return Err(bad("observed_annotation_output_changed"));
        }
        if let Some((raw, original)) = &self.marker {
            // Malformed/non-regular late markers are execution control failures.
            // Still refuse a changed parent or a new alias to an input/output.
            let p = raw_path(raw)?;
            let parent = p
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            let resolved = fs::canonicalize(parent)
                .map_err(io)?
                .join(p.file_name().unwrap());
            if resolved != original.resolved {
                return Err(bad("observed_annotation_control_path_changed"));
            }
            let identity = same_file::Handle::from_path(&resolved).ok();
            let m = Location { resolved, identity };
            if aliases(&m, &now) || self.inputs.iter().any(|(_, p)| aliases(&m, p)) {
                return Err(bad("observed_annotation_control_collision"));
            }
        }
        Ok(())
    }
}

struct BoundedWriter {
    writer: HashedWriter,
    left: u64,
}
impl Write for BoundedWriter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() as u64 > self.left {
            return Err(std::io::Error::other("observed_annotation_output_budget"));
        }
        let n = self.writer.write(b)?;
        self.left -= n as u64;
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }
}

pub(super) fn write_json<T: Serialize>(
    path: &Path,
    name: &str,
    value: &T,
    limit: u64,
) -> Result<RnaReadPanelArtifact, EngineError> {
    let mut writer = BoundedWriter {
        writer: HashedWriter::create(path)?,
        left: limit,
    };
    serde_json::to_writer_pretty(&mut writer, value).map_err(|e| bad(e.to_string()))?;
    writer.writer.receipt(name.into())
}

#[derive(Clone, Copy, Debug)]
pub(super) enum Phase {
    BeforeReserve,
    Reserved,
    ArtifactsVisible,
    ReportReady,
}

/// Finalization hooks are internal fault/visibility probes, never adapter callbacks.
pub(super) fn publish(
    stage: tempfile::TempDir,
    paths: &Paths,
    report: &ObservedTranscriptAnnotationReport,
    hook: &mut dyn FnMut(Phase) -> Result<(), EngineError>,
) -> Result<PrimerPayloadFile, EngineError> {
    let artifact_bytes = report.annotation.as_ref().map_or(0, |a| a.bytes);
    let remaining = report
        .request
        .budgets
        .max_output_bytes
        .checked_sub(artifact_bytes)
        .ok_or_else(|| bad("observed_annotation_output_budget"))?;
    // Write/sync before reserving the public destination; the final report
    // rename is the only completion marker, not directory iteration order.
    let report_stage = tempfile::Builder::new()
        .prefix(".gentle-annotation-report-")
        .tempdir_in(stage.path().parent().unwrap())
        .map_err(io)?;
    let receipt = write_json(
        &report_stage.path().join("report.json"),
        "report.json",
        report,
        remaining,
    )?;
    hook(Phase::BeforeReserve)?;
    paths.recheck_output()?;
    let out = &paths.output.resolved;
    let stage_owner = same_file::Handle::from_path(stage.path()).map_err(io)?;
    let annotation_owner = report
        .annotation
        .as_ref()
        .map(|_| same_file::Handle::from_path(stage.path().join("annotation.json")).map_err(io))
        .transpose()?;
    fs::create_dir(out).map_err(io)?;
    let owner = match same_file::Handle::from_path(out) {
        Ok(h) => h,
        Err(e) => {
            let _ = fs::remove_dir(out);
            return Err(bad(format!(
                "observed_annotation_reservation_identity: {e}; output reservation may remain at {}",
                out.display()
            )));
        }
    };
    let result: Result<PrimerPayloadFile, EngineError> = (|| {
        hook(Phase::Reserved)?;
        fs::rename(stage.path(), out.join("artifacts")).map_err(io)?;
        // The old temporary name no longer belongs to us. Disarm its automatic
        // recursive cleanup before any observer could recreate that name.
        let _ = stage.keep();
        hook(Phase::ArtifactsVisible)?;
        hook(Phase::ReportReady)?;
        // No fallible work follows this no-clobber completion marker.
        // Hard linking is same-filesystem and refuses an existing report.
        fs::hard_link(
            report_stage.path().join("report.json"),
            out.join("report.json"),
        )
        .map_err(io)?;
        Ok(PrimerPayloadFile {
            path: out.join("report.json").display().to_string(),
            sha256: receipt.sha256,
            bytes: receipt.bytes,
        })
    })();
    if result.is_err() {
        // Never recursively delete a directory that another actor could have
        // populated/replaced. Remove only known, identity-matching owned files.
        let same = |p: &Path, h: &same_file::Handle| {
            fs::symlink_metadata(p).is_ok_and(|m| !m.file_type().is_symlink())
                && same_file::Handle::from_path(p).is_ok_and(|now| &now == h)
        };
        if same(out, &owner) {
            let artifacts = out.join("artifacts");
            if same(&artifacts, &stage_owner) {
                if let Some(h) = &annotation_owner {
                    let p = artifacts.join("annotation.json");
                    if same(&p, h) {
                        let _ = fs::remove_file(p);
                    }
                }
                let _ = fs::remove_dir(artifacts);
            }
            let _ = fs::remove_dir(out);
        }
        if fs::symlink_metadata(out).is_ok() {
            let e = result.unwrap_err();
            return Err(bad(format!(
                "{}; observed_annotation_orphaned_output: {}",
                e.message,
                out.display()
            )));
        }
    }
    result
}
