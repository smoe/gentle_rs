//! SAM-only bounded artifact writing, owned report-last publication and saved-only inspection.

use super::*;

struct Output {
    writer: HashedWriter,
    left: u64,
}
impl Write for Output {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() as u64 > self.left {
            return Err(std::io::Error::other("observed_sam_output_budget"));
        }
        let n = self.writer.write(bytes)?;
        self.left -= n as u64;
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }
}

fn table<T: Serialize>(
    root: &Path,
    name: &str,
    rows: impl IntoIterator<Item = T>,
    r: &Request,
    left: &mut u64,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<RnaReadPanelArtifact, EngineError> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(Output {
            writer: HashedWriter::create(&root.join(name))?,
            left: *left,
        });
    writer
        .write_record(["row_json"])
        .map_err(|e| bad(e.to_string()))?;
    for row in rows {
        check()?;
        size(
            &row,
            (r.budgets.max_output_row_bytes as u64)
                .min((r.budgets.max_retained_bytes as u64 - 65536) / 64),
        )?;
        let text = serde_json::to_string(&row).map_err(|e| bad(e.to_string()))?;
        writer
            .write_record([text])
            .map_err(|e| bad(e.to_string()))?;
    }
    writer.flush().map_err(|e| bad(e.to_string()))?;
    let writer = writer.into_inner().map_err(|e| bad(e.to_string()))?;
    *left = writer.left;
    writer.writer.receipt(name.into())
}

pub(super) fn artifacts(
    root: &Path,
    data: &ObservedSamData,
    status: &str,
    r: &Request,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Vec<RnaReadPanelArtifact>, EngineError> {
    check()?;
    let mut writer = Output {
        writer: HashedWriter::create(&root.join(DATA))?,
        left: r.budgets.max_output_bytes,
    };
    serde_json::to_writer(&mut writer, data).map_err(|e| bad(e.to_string()))?;
    let mut left = writer.left;
    let mut artifacts = vec![writer.writer.receipt(DATA.into())?];
    artifacts.push(table(
        root,
        TABLES[0],
        data.analysis.iter().flat_map(|a| &a.units),
        r,
        &mut left,
        check,
    )?);
    artifacts.push(table(
        root,
        TABLES[1],
        data.analysis.iter().flat_map(|a| &a.observations),
        r,
        &mut left,
        check,
    )?);
    artifacts.push(table(
        root,
        TABLES[2],
        data.analysis.iter().flat_map(|a| &a.chain_families),
        r,
        &mut left,
        check,
    )?);
    let compatibility = data.analysis.iter().flat_map(|a| a.gene_resolutions.iter().map(|c| json!({"kind":"source","comparison":c}))
        .chain(a.observations.iter().flat_map(|o| o.comparisons.iter().map(move |c| json!({"kind":"observation","observation_id":o.observation_id,"comparison":c})))));
    artifacts.push(table(root, TABLES[3], compatibility, r, &mut left, check)?);
    artifacts.push(table(
        root,
        TABLES[4],
        data.analysis.iter().flat_map(|a| &a.junctions),
        r,
        &mut left,
        check,
    )?);
    let unavailable = data
        .analysis
        .iter()
        .flat_map(|a| {
            a.units
                .iter()
                .filter(|u| !u.all_query_views_assessed)
                .map(|u| json!({"kind":"unit","unit":u}))
        })
        .chain(
            data.intake
                .iter()
                .flat_map(|i| &i.missing_sa_companions)
                .map(|c| json!({"kind":"missing_sa_companion","value":c})),
        )
        .chain(
            data.annotations
                .iter()
                .filter(|a| a.availability != ObservedSamAnnotationAvailability::Available)
                .map(|a| json!({"kind":"annotation","value":a})),
        )
        .chain(
            data.analysis
                .is_none()
                .then(|| json!({"kind":"analysis_unavailable","execution_status":status})),
        );
    artifacts.push(table(root, TABLES[5], unavailable, r, &mut left, check)?);
    Ok(artifacts)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Phase {
    BeforeReserve,
    ReservationCreated,
    ArtifactsVisible,
    ReportReady,
}

fn same(path: &Path, owner: &same_file::Handle) -> bool {
    fs::symlink_metadata(path).is_ok_and(|m| !m.file_type().is_symlink())
        && same_file::Handle::from_path(path).is_ok_and(|now| &now == owner)
}

pub(super) fn publish(
    stage: tempfile::TempDir,
    report: &ObservedSamReport,
    hook: &mut dyn FnMut(Phase) -> Result<(), EngineError>,
) -> Result<(), EngineError> {
    let left = report
        .artifacts
        .iter()
        .try_fold(report.request.budgets.max_output_bytes, |n, a| {
            n.checked_sub(a.bytes)
        })
        .ok_or_else(|| bad("observed_sam_output_budget"))?;
    let raw = serde_json::to_vec_pretty(report).map_err(|e| bad(e.to_string()))?;
    if raw.len() as u64 > left
        || raw
            .len()
            .checked_mul(64)
            .and_then(|n| n.checked_add(65536))
            .is_none_or(|n| n > report.request.budgets.max_retained_bytes)
    {
        return Err(bad("observed_sam_output_budget"));
    }
    let parent = stage
        .path()
        .parent()
        .ok_or_else(|| bad("observed_sam_output_parent"))?;
    let temporary = tempfile::Builder::new()
        .prefix(".gentle-sam-report-")
        .tempdir_in(parent)
        .map_err(|e| bad(e.to_string()))?;
    {
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(temporary.path().join("report.json"))
            .map_err(|e| bad(e.to_string()))?;
        file.write_all(&raw).map_err(|e| bad(e.to_string()))?;
        file.sync_all().map_err(|e| bad(e.to_string()))?;
    }
    hook(Phase::BeforeReserve)?;
    let out = Path::new(&report.resolved_output_dir);
    fs::create_dir(out).map_err(|e| bad(e.to_string()))?;
    let owner = same_file::Handle::from_path(out).map_err(|e| {
        bad(format!(
            "observed_sam_reservation_identity: {e}; observed_sam_orphaned_output: {}",
            out.display()
        ))
    })?;
    let mut owned = None;
    let result = (|| {
        hook(Phase::ReservationCreated)?;
        if !same(out, &owner) {
            return Err(bad("observed_sam_output_identity"));
        }
        // No staging-directory or artifact handle is open across this rename.
        fs::rename(stage.path(), out.join("artifacts")).map_err(|e| bad(e.to_string()))?;
        let _ = stage.keep();
        let directory =
            same_file::Handle::from_path(out.join("artifacts")).map_err(|e| bad(e.to_string()))?;
        let files = report
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
        owned = Some((directory, files));
        hook(Phase::ArtifactsVisible)?;
        hook(Phase::ReportReady)?;
        if !same(out, &owner) {
            return Err(bad("observed_sam_output_identity"));
        }
        let (directory, files) = owned
            .as_ref()
            .ok_or_else(|| bad("observed_sam_artifact_identity"))?;
        let root = out.join("artifacts");
        if !same(&root, directory) {
            return Err(bad("observed_sam_artifact_identity"));
        }
        for (name, handle) in files {
            if !same(&root.join(name), handle) {
                return Err(bad("observed_sam_artifact_identity"));
            }
        }
        for artifact in &report.artifacts {
            let mut reader = VerifiedReader::open(
                &artifact_file(&root, artifact)?,
                report.request.budgets.max_output_bytes,
            )?;
            std::io::copy(&mut reader, &mut std::io::sink()).map_err(|e| bad(e.to_string()))?;
            reader.finish()?;
        }
        fs::hard_link(
            temporary.path().join("report.json"),
            out.join("report.json"),
        )
        .map_err(|e| bad(e.to_string()))?;
        Ok(())
    })();
    if let Err(error) = result {
        if same(out, &owner) {
            if let Some((directory, files)) = owned.take() {
                let root = out.join("artifacts");
                if same(&root, &directory) {
                    for (name, handle) in files {
                        let path = root.join(name);
                        let remove = same(&path, &handle);
                        drop(handle);
                        if remove {
                            let _ = fs::remove_file(path);
                        }
                    }
                    drop(directory);
                    let _ = fs::remove_dir(root);
                }
            }
            drop(owner);
            let _ = fs::remove_dir(out);
        }
        return if fs::symlink_metadata(out).is_ok() {
            Err(bad(format!(
                "{}; observed_sam_orphaned_output: {}",
                error.message,
                out.display()
            )))
        } else {
            Err(error)
        };
    }
    Ok(())
}

impl GentleEngine {
    /// Inspect only saved report/artifact bytes; SAM, GTF and preparation sources may be offline.
    pub fn inspect_observed_transcript_structures_from_sam(
        path: &str,
        expected: Option<&ObservedSamResultBinding>,
        table: &str,
        rows: usize,
        verify: bool,
    ) -> Result<ObservedSamInspection, EngineError> {
        kernel::validate_sam_raw_path(path)?;
        if !TABLES.contains(&table) || rows > 1000 {
            return Err(bad("observed_sam_inspection_bounds"));
        }
        let metadata = regular(Path::new(path))?;
        let limit = (CAP - 65536) / 64;
        if metadata.len() > limit {
            return Err(bad("observed_sam_report_memory_budget"));
        }
        let mut raw = Vec::new();
        fs::File::open(path)
            .map_err(|e| bad(e.to_string()))?
            .take(limit + 1)
            .read_to_end(&mut raw)
            .map_err(|e| bad(e.to_string()))?;
        if raw.len() as u64 > limit {
            return Err(bad("observed_sam_report_memory_budget"));
        }
        let file = PrimerPayloadFile {
            path: path.into(),
            bytes: raw.len() as u64,
            sha256: sha256_prefixed_bytes(&raw),
        };
        drop(raw);
        let report: ObservedSamReport =
            kernel::read_sam_bound_json(&file, CAP as usize, &mut || Ok(()))?;
        let binding = validate_report(&report)?;
        if fs::canonicalize(path).map_err(|e| bad(e.to_string()))?
            != fs::canonicalize(&binding.report_path).map_err(|e| bad(e.to_string()))?
            || expected.is_some_and(|e| encode(e).ok() != encode(&binding).ok())
        {
            return Err(bad("observed_sam_inspection_binding"));
        }
        let root = Path::new(path)
            .parent()
            .unwrap_or(Path::new("."))
            .join("artifacts");
        let mut selected = Vec::new();
        let mut memory = 65536usize;
        let mut total = 0u64;
        for artifact in &report.artifacts {
            if !verify && !(artifact.path == table && rows > 0) {
                continue;
            }
            total = total
                .checked_add(artifact.bytes)
                .filter(|n| *n <= report.request.budgets.max_output_bytes)
                .ok_or_else(|| bad("observed_sam_inspection_budget"))?;
            let reader = VerifiedReader::open(
                &artifact_file(&root, artifact)?,
                report.request.budgets.max_output_bytes,
            )?;
            if artifact.path == table && rows > 0 {
                let mut reader = BufReader::new(reader);
                let mut header = true;
                loop {
                    let mut line = Vec::new();
                    let encoded_limit = 2 * report.request.budgets.max_output_row_bytes as u64 + 8;
                    let n = reader
                        .by_ref()
                        .take(encoded_limit + 1)
                        .read_until(b'\n', &mut line)
                        .map_err(|e| bad(e.to_string()))?;
                    if n == 0 {
                        break;
                    }
                    if line.len() as u64 > encoded_limit {
                        return Err(bad("observed_sam_row_budget"));
                    }
                    let mut csv = csv::ReaderBuilder::new()
                        .delimiter(b'\t')
                        .has_headers(false)
                        .from_reader(line.as_slice());
                    let mut records = csv.records();
                    let row = records
                        .next()
                        .ok_or_else(|| bad("observed_sam_table_columns"))?
                        .map_err(|e| bad(e.to_string()))?;
                    if records.next().is_some() || row.len() != 1 {
                        return Err(bad("observed_sam_table_columns"));
                    }
                    if header {
                        if &row[0] != "row_json" {
                            return Err(bad("observed_sam_table_columns"));
                        }
                        header = false;
                        continue;
                    }
                    if row[0].len() > report.request.budgets.max_output_row_bytes {
                        return Err(bad("observed_sam_row_budget"));
                    }
                    if selected.len() < rows {
                        memory = row[0]
                            .len()
                            .checked_mul(64)
                            .and_then(|n| n.checked_add(memory))
                            .filter(|n| *n <= report.request.budgets.max_retained_bytes)
                            .ok_or_else(|| bad("observed_sam_response_budget"))?;
                        selected
                            .push(serde_json::from_str(&row[0]).map_err(|e| bad(e.to_string()))?);
                    }
                }
                if header {
                    return Err(bad("observed_sam_table_columns"));
                }
                reader.into_inner().finish()?;
            } else {
                let mut reader = reader;
                std::io::copy(&mut reader, &mut std::io::sink()).map_err(|e| bad(e.to_string()))?;
                reader.finish()?;
            }
        }
        Ok(ObservedSamInspection {
            report,
            artifacts_verified: verify,
            table: table.into(),
            rows: selected,
        })
    }

    pub fn cancel_observed_transcript_structures_from_sam(
        r: &Request,
    ) -> Result<String, EngineError> {
        admit(r)?;
        paths(r)?;
        let path = r
            .cancel_path
            .as_ref()
            .ok_or_else(|| bad("observed_sam_cancel_path_required"))?;
        let identity = digest(r)?;
        if control::marker_requested(Some(path), std::slice::from_ref(&identity))? {
            return Ok("Cancellation already requested; termination not confirmed".into());
        }
        let parent = Path::new(path)
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let mut file = tempfile::NamedTempFile::new_in(parent).map_err(|e| bad(e.to_string()))?;
        file.write_all(identity.as_bytes())
            .map_err(|e| bad(e.to_string()))?;
        file.as_file().sync_all().map_err(|e| bad(e.to_string()))?;
        paths(r)?;
        file.persist_noclobber(path)
            .map_err(|e| bad(e.to_string()))?;
        Ok("Cancellation requested for the exact SAM request; termination not confirmed".into())
    }
}
