//! J1 bounded artifact publication and saved-only inspection; never reinterprets S1 evidence.

use super::*;

pub(super) fn raw_path(raw: &str) -> Result<&Path, EngineError> {
    if raw.trim().is_empty() || raw.len() > 32768 || raw.chars().any(char::is_control) {
        return Err(bad("junction_assessment_raw_path"));
    }
    let p = Path::new(raw);
    let prefix = match p.components().next() {
        Some(std::path::Component::Prefix(p)) => p.as_os_str().len(),
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
        return Err(bad("junction_assessment_raw_path"));
    }
    Ok(p)
}
pub(super) fn resolved(raw: &str) -> Result<PathBuf, EngineError> {
    let p = raw_path(raw)?;
    let parent = p
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let parent = fs::canonicalize(parent).map_err(|e| bad(e.to_string()))?;
    Ok(parent.join(
        p.file_name()
            .ok_or_else(|| bad("junction_assessment_raw_path"))?,
    ))
}
pub(super) fn execution_paths(r: &AssessObservedSpliceJunctionsRequest) -> Result<(), EngineError> {
    let Some(marker) = &r.cancel_path else {
        return Ok(());
    };
    let marker = resolved(marker)?;
    let out = resolved(&r.output_dir)?;
    let fold = |p: &Path| p.to_string_lossy().to_lowercase();
    if marker == out
        || marker.starts_with(&out)
        || ((cfg!(windows) || cfg!(target_os = "macos"))
            && (fold(&marker) == fold(&out)
                || fold(&marker).starts_with(&(fold(&out) + std::path::MAIN_SEPARATOR_STR))))
    {
        return Err(bad("junction_assessment_control_collision"));
    }
    let ObservedJunctionSource::S1Report { report, .. } = &r.source;
    for file in std::iter::once(report)
        .chain(r.reference_sequence.iter().map(|f| &f.file))
        .chain(r.sample_manifest.iter())
    {
        if marker == resolved(&file.path)?
            || same_file::is_same_file(&marker, &file.path).unwrap_or(false)
        {
            return Err(bad("junction_assessment_control_collision"));
        }
    }
    Ok(())
}

struct Output {
    writer: HashedWriter,
    left: u64,
}
impl Write for Output {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() as u64 > self.left {
            return Err(std::io::Error::other("junction_assessment_output_budget"));
        }
        let n = self.writer.write(b)?;
        self.left -= n as u64;
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }
}
fn json_file<T: Serialize>(
    root: &Path,
    name: &str,
    value: &T,
    left: &mut u64,
) -> Result<RnaReadPanelArtifact, EngineError> {
    let mut w = Output {
        writer: HashedWriter::create(&root.join(name))?,
        left: *left,
    };
    serde_json::to_writer(&mut w, value).map_err(|e| bad(e.to_string()))?;
    *left = w.left;
    w.writer.receipt(name.into())
}
fn table<T: Serialize>(
    root: &Path,
    name: &str,
    rows: &[T],
    left: &mut u64,
) -> Result<RnaReadPanelArtifact, EngineError> {
    let mut w = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Output {
            writer: HashedWriter::create(&root.join(name))?,
            left: *left,
        });
    w.write_record(["row_json"])
        .map_err(|e| bad(e.to_string()))?;
    for r in rows {
        let json = serde_json::to_string(r).map_err(|e| bad(e.to_string()))?;
        if json.len() > 8_000_000 {
            return Err(bad("junction_assessment_output_row_budget"));
        }
        w.write_record([json]).map_err(|e| bad(e.to_string()))?;
    }
    w.flush().map_err(|e| bad(e.to_string()))?;
    let w = w.into_inner().map_err(|e| bad(e.to_string()))?;
    *left = w.left;
    w.writer.receipt(name.into())
}
pub(super) fn artifacts(
    root: &Path,
    data: &ObservedJunctionAssessmentData,
    r: &ObservedSpliceJunctionAssessmentReport,
) -> Result<Vec<RnaReadPanelArtifact>, EngineError> {
    let mut left = r.request.budgets.max_output_bytes;
    let mut artifacts = vec![json_file(root, "assessment.json", data, &mut left)?];
    artifacts.push(table(root, TABLES[0], &data.junctions, &mut left)?);
    artifacts.push(table(root, TABLES[1], &data.chains, &mut left)?);
    let samples = data
        .libraries
        .iter()
        .map(|v| json!({"kind":"library_denominator","value":v}))
        .chain(data.junctions.iter().flat_map(|j| {
            j.grouped_support.iter().map(
                move |v| json!({"kind":"junction_support","junction_id":j.junction_id,"value":v}),
            )
        }))
        .chain(
            data.unused_manifest_rows
                .iter()
                .map(|v| json!({"kind":"unused_manifest_row","value":v})),
        )
        .collect::<Vec<_>>();
    artifacts.push(table(root, TABLES[2], &samples, &mut left)?);
    let mut unavailable = data
        .unassessed_contexts
        .iter()
        .map(|v| json!({"kind":"unassessed_context","value":v}))
        .chain(
            data.source_unavailable
                .iter()
                .map(|v| json!({"kind":"s1_unavailable","value":v})),
        )
        .collect::<Vec<_>>();
    if r.counts.is_none() {
        unavailable.push(json!({"kind":"assessment_unavailable","status":r.status,"source_status":r.source_status}));
    }
    artifacts.push(table(root, TABLES[3], &unavailable, &mut left)?);
    Ok(artifacts)
}

#[derive(Clone, Copy, Debug)]
pub(super) enum Phase {
    BeforeReserve,
    ArtifactsVisible,
    ReportReady,
}
pub(super) fn publish(
    stage: tempfile::TempDir,
    r: &ObservedSpliceJunctionAssessmentReport,
    hook: &mut dyn FnMut(Phase) -> Result<(), EngineError>,
) -> Result<(), EngineError> {
    let remaining = r
        .artifacts
        .iter()
        .try_fold(r.request.budgets.max_output_bytes, |left, a| {
            left.checked_sub(a.bytes)
        })
        .ok_or_else(|| bad("junction_assessment_output_budget"))?;
    let bytes = serde_json::to_vec_pretty(r).map_err(|e| bad(e.to_string()))?;
    if bytes.len() as u64 > remaining {
        return Err(bad("junction_assessment_output_budget"));
    }
    let parent = stage
        .path()
        .parent()
        .ok_or_else(|| bad("junction_assessment_output_parent"))?;
    let report_stage = tempfile::Builder::new()
        .prefix(".gentle-junction-report-")
        .tempdir_in(parent)
        .map_err(|e| bad(e.to_string()))?;
    {
        let mut f = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(report_stage.path().join("report.json"))
            .map_err(|e| bad(e.to_string()))?;
        f.write_all(&bytes).map_err(|e| bad(e.to_string()))?;
        f.sync_all().map_err(|e| bad(e.to_string()))?;
    }
    hook(Phase::BeforeReserve)?;
    let out = Path::new(&r.request.output_dir);
    fs::create_dir(out).map_err(|e| bad(e.to_string()))?;
    let owner = match same_file::Handle::from_path(out) {
        Ok(h) => h,
        Err(e) => {
            let _ = fs::remove_dir(out);
            return Err(bad(format!(
                "junction_assessment_reservation_identity: {e}"
            )));
        }
    };
    let mut artifacts_owner = None;
    let result: Result<(), EngineError> = (|| {
        fs::rename(stage.path(), out.join("artifacts")).map_err(|e| bad(e.to_string()))?;
        let _ = stage.keep();
        let dir =
            same_file::Handle::from_path(out.join("artifacts")).map_err(|e| bad(e.to_string()))?;
        let files = r
            .artifacts
            .iter()
            .map(|a| {
                Ok((
                    a.path.clone(),
                    same_file::Handle::from_path(out.join("artifacts").join(&a.path))
                        .map_err(|e| bad(e.to_string()))?,
                ))
            })
            .collect::<Result<Vec<_>, EngineError>>()?;
        artifacts_owner = Some((dir, files));
        hook(Phase::ArtifactsVisible)?;
        hook(Phase::ReportReady)?;
        fs::hard_link(
            report_stage.path().join("report.json"),
            out.join("report.json"),
        )
        .map_err(|e| bad(e.to_string()))?;
        Ok(())
    })();
    if let Err(e) = result {
        let same = |p: &Path, h: &same_file::Handle| {
            fs::symlink_metadata(p).is_ok_and(|m| !m.file_type().is_symlink())
                && same_file::Handle::from_path(p).is_ok_and(|now| &now == h)
        };
        if same(out, &owner) {
            if let Some((dir, files)) = &artifacts_owner {
                let a = out.join("artifacts");
                if same(&a, dir) {
                    for (name, h) in files {
                        let p = a.join(name);
                        if same(&p, h) {
                            let _ = fs::remove_file(p);
                        }
                    }
                    let _ = fs::remove_dir(a);
                }
            }
            let _ = fs::remove_dir(out);
        }
        return if fs::symlink_metadata(out).is_ok() {
            Err(bad(format!(
                "{}; junction_assessment_orphaned_output: {}",
                e.message,
                out.display()
            )))
        } else {
            Err(e)
        };
    }
    Ok(())
}

impl GentleEngine {
    /// Saved J1 inspection does not reopen source reports, annotations, samples or reference FASTA.
    pub fn inspect_observed_splice_junctions(
        path: &str,
        expected: Option<&ObservedSpliceJunctionResultBinding>,
        table: &str,
        rows: usize,
        verify: bool,
    ) -> Result<ObservedSpliceJunctionInspection, EngineError> {
        raw_path(path)?;
        regular(Path::new(path))?;
        if !TABLES.contains(&table) || rows > 1000 {
            return Err(bad("junction_assessment_inspection_bounds"));
        }
        let report: ObservedSpliceJunctionAssessmentReport = load(path)?;
        let binding = validate_report(&report)?;
        if fs::canonicalize(path).map_err(|e| bad(e.to_string()))?
            != fs::canonicalize(&binding.report_path).map_err(|e| bad(e.to_string()))?
            || expected.is_some_and(|e| encode(e).ok() != encode(&binding).ok())
        {
            return Err(bad("junction_assessment_inspection_binding"));
        }
        let root = Path::new(path)
            .parent()
            .unwrap_or(Path::new("."))
            .join("artifacts");
        let mut selected = Vec::new();
        let (mut total, mut size) = (0u64, 0usize);
        for a in &report.artifacts {
            if !verify && !(a.path == table && rows > 0) {
                continue;
            }
            total = total
                .checked_add(a.bytes)
                .filter(|n| *n <= report.request.budgets.max_output_bytes)
                .ok_or_else(|| bad("junction_assessment_inspection_budget"))?;
            let reader = VerifiedReader::open(
                &artifact_file(&root, a)?,
                report.request.budgets.max_output_bytes,
            )?;
            if a.path == table && rows > 0 {
                let mut csv = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .from_reader(reader);
                if csv
                    .headers()
                    .map_err(|e| bad(e.to_string()))?
                    .iter()
                    .collect::<Vec<_>>()
                    != ["row_json"]
                {
                    return Err(bad("junction_assessment_table_columns"));
                }
                for row in csv.records() {
                    let row = row.map_err(|e| bad(e.to_string()))?;
                    if row.len() != 1 || row[0].len() > 8_000_000 {
                        return Err(bad("junction_assessment_row_budget"));
                    }
                    if selected.len() < rows {
                        size = size
                            .checked_add(row[0].len())
                            .filter(|n| *n <= 8_000_000)
                            .ok_or_else(|| bad("junction_assessment_response_budget"))?;
                        selected
                            .push(serde_json::from_str(&row[0]).map_err(|e| bad(e.to_string()))?);
                    }
                }
                csv.into_inner().finish()?;
            } else {
                let mut reader = reader;
                std::io::copy(&mut reader, &mut std::io::sink()).map_err(|e| bad(e.to_string()))?;
                reader.finish()?;
            }
        }
        Ok(ObservedSpliceJunctionInspection {
            report,
            artifacts_verified: verify,
            table: table.into(),
            rows: selected,
        })
    }

    pub fn cancel_observed_splice_junctions(
        r: &AssessObservedSpliceJunctionsRequest,
    ) -> Result<String, EngineError> {
        admit(r)?;
        execution_paths(r)?;
        let identity = hash(r)?;
        let path = r
            .cancel_path
            .as_ref()
            .ok_or_else(|| bad("junction_assessment_cancel_path_required"))?;
        if control::marker_requested(Some(path), std::slice::from_ref(&identity))? {
            return Ok("Cancellation already requested; termination not confirmed".into());
        }
        let parent = Path::new(path)
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let mut f = tempfile::NamedTempFile::new_in(parent).map_err(|e| bad(e.to_string()))?;
        f.write_all(identity.as_bytes())
            .map_err(|e| bad(e.to_string()))?;
        f.as_file().sync_all().map_err(|e| bad(e.to_string()))?;
        execution_paths(r)?;
        f.persist_noclobber(path).map_err(|e| bad(e.to_string()))?;
        Ok("Cancellation requested for the exact J1 request; termination not confirmed".into())
    }
}
