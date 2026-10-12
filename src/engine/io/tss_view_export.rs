//! Shared native TSS export with validated collection membership and staged pages.
//! Presentation never changes loaded sequences, approval or scientific verdicts.

use super::*;
use crate::{
    digest_utils::sha256_hex_bytes,
    tss_sequence_view::{
        TSS_SVG_HOVER_LIMIT_PER_STRAND, TssLocalScoreRequest, TssSequenceView, TssViewSvgOptions,
        read_tss_profile_file, write_tss_view_svg,
    },
};
use gentle_protocol::tss_workspace::{TssViewSvgExportPage, TssViewSvgExportReceipt};
use std::{fs, io::Write, path::Path};

const MAX_MEMBERS: usize = 32;
const MAX_BUNDLE_BYTES: u64 = 32 * 1024 * 1024;
const SCALE_POLICY: &str = "per_window_per_lane";

pub(super) struct NativeTssExportOutcome {
    pub receipt: TssViewSvgExportReceipt,
    pub messages: Vec<String>,
    pub warnings: Vec<String>,
}

fn digest(value: &impl Serialize) -> Result<String, EngineError> {
    serde_json::to_vec(value)
        .map(|bytes| sha256_hex_bytes(&bytes))
        .map_err(|e| EngineError::internal(format!("TSS SVG receipt serialization: {e}")))
}

fn io_error(error: impl std::fmt::Display) -> EngineError {
    EngineError::new(ErrorCode::Io, format!("Native TSS export: {error}"))
}

fn checkpoint(
    on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    seq_id: &str,
    completed: usize,
    total: usize,
) -> Result<(), EngineError> {
    let percent = 100.0 * completed as f64 / total as f64;
    if !on_progress(OperationProgress::Tfbs(TfbsProgress {
        seq_id: seq_id.into(),
        motif_id: String::new(),
        motif_index: 0,
        motif_count: 0,
        scanned_steps: completed,
        total_steps: total,
        motif_percent: 0.0,
        total_percent: percent,
        task_kind: Some("tss_svg_export".into()),
        stage_label: Some("staging_pages".into()),
        detail: Some(format!(
            "{completed}/{total} pages staged; not yet published"
        )),
        stage_percent: Some(percent),
    })) {
        return Err(EngineError::invalid_input(
            "Native TSS SVG export cancelled before publication; no success receipt published",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_root() -> tempfile::TempDir {
        // Avoid macOS's system /var alias; deliberate symlinks are tested below.
        tempfile::tempdir_in(fs::canonicalize(std::env::temp_dir()).unwrap()).unwrap()
    }

    // Reuse entirely synthetic annotated loci; no private study/production data.
    fn collection_engine(reverse: bool) -> GentleEngine {
        let mut engine = synthetic_tss_engine(reverse);
        let request = synthetic_tss_approval(&engine);
        engine
            .apply(Operation::MaterializeTssWindows { request })
            .unwrap();
        engine
    }

    fn operation(path: &Path) -> Operation {
        Operation::ExportTssViewSvg {
            seq_id: String::new(),
            collection_id: Some("toy_tss".into()),
            path: path.to_string_lossy().into_owned(),
            report: None,
            local_motifs: vec![],
            score_kind: TfbsScoreTrackValueKind::LlrBits,
            clip_negative: true,
            start_0based: None,
            end_0based_exclusive: None,
            width_px: None,
        }
    }

    #[test]
    fn native_tss_collection_export_preserves_members_scales_hashes_and_state() {
        for reverse in [false, true] {
            let mut engine = collection_engine(reverse);
            let dir = temp_root();
            let output = dir.path().join("ordered-pages");
            let before = serde_json::to_value(engine.state()).unwrap();
            let op = operation(&output);
            let request_sha256 = digest(&op).unwrap();
            let result = engine.apply(op).unwrap();
            let receipt = result.tss_view_svg_export.unwrap();
            assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
            assert_eq!(receipt.request_sha256, request_sha256);
            assert_eq!(receipt.scale_policy, SCALE_POLICY);
            let collection = engine.get_tss_collection("toy_tss").unwrap();
            assert_eq!(
                receipt.collection_sha256.as_deref(),
                Some(digest(&collection).unwrap().as_str())
            );
            assert_eq!(receipt.pages.len(), collection.members.len());
            for (i, (page, member)) in receipt.pages.iter().zip(&collection.members).enumerate() {
                assert_eq!(page.seq_id, member.tss.output_seq_id);
                assert_eq!(page.filename, format!("{:03}.svg", i + 1));
                let bytes = fs::read(output.join(&page.filename)).unwrap();
                assert_eq!(page.svg_sha256, sha256_hex_bytes(&bytes));
                let view =
                    TssSequenceView::from_dna(&engine.state.sequences[&page.seq_id]).unwrap();
                assert_eq!(page.sequence_sha256, view.sequence_sha256);
                assert_eq!(page.view_sha256, digest(&view).unwrap());
                assert_eq!(page.geometry, view.geometry);
                let svg = String::from_utf8(bytes).unwrap();
                assert!(svg.contains("preserve each window&apos;s supplied lane scales"));
                assert!(svg.contains(&page.geometry.genomic_at(0).unwrap().to_string()));
            }
            let index = fs::read(output.join("index.html")).unwrap();
            assert_eq!(receipt.index_html_sha256, Some(sha256_hex_bytes(&index)));
            let disk: TssViewSvgExportReceipt =
                serde_json::from_slice(&fs::read(output.join("receipt.json")).unwrap()).unwrap();
            assert_eq!(
                serde_json::to_value(disk).unwrap(),
                serde_json::to_value(&receipt).unwrap()
            );
            assert!(
                engine.apply(operation(&output)).is_err(),
                "never overwrite an existing bundle"
            );
        }
    }

    #[test]
    fn native_tss_collection_export_rejects_stale_and_mixed_targets_without_artifacts() {
        let mut engine = collection_engine(false);
        let dir = temp_root();
        let output = dir.path().join("refused");
        let mut mixed = operation(&output);
        if let Operation::ExportTssViewSvg { seq_id, .. } = &mut mixed {
            *seq_id = "locus".into();
        }
        assert!(
            engine
                .apply(mixed)
                .unwrap_err()
                .message
                .contains("exactly one")
        );
        assert!(!output.exists());
        let member = engine.get_tss_collection("toy_tss").unwrap().members[0]
            .tss
            .output_seq_id
            .clone();
        let mut record = engine.state.sequences[&member].clone_seq_record();
        record.seq[0] = if record.seq[0] == b'A' { b'T' } else { b'A' };
        engine
            .state
            .sequences
            .insert(member, DNAsequence::from_genbank_seq(record));
        assert!(
            engine
                .apply(operation(&output))
                .unwrap_err()
                .message
                .contains("edited")
        );
        assert!(!output.exists());
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn native_tss_collection_export_rejects_invalid_shared_span_before_staging() {
        let mut engine = collection_engine(false);
        let dir = temp_root();
        let output = dir.path().join("invalid-span");
        let mut op = operation(&output);
        if let Operation::ExportTssViewSvg {
            end_0based_exclusive,
            ..
        } = &mut op
        {
            *end_0based_exclusive = Some(1001);
        }
        assert!(
            engine
                .apply(op)
                .unwrap_err()
                .message
                .contains("Invalid TSS SVG span")
        );
        assert!(!output.exists());
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn native_tss_collection_export_cancellation_cleans_staged_pages() {
        let engine = collection_engine(false);
        let dir = temp_root();
        let output = dir.path().join("cancelled");
        let mut saw_staged_page = false;
        let result = engine.export_native_tss_svg(&operation(&output), &mut |progress| {
            if let OperationProgress::Tfbs(p) = progress {
                if p.task_kind.as_deref() == Some("tss_svg_export") && p.scanned_steps > 0 {
                    saw_staged_page = true;
                    return false;
                }
            }
            true
        });
        assert!(saw_staged_page);
        assert!(
            result
                .err()
                .unwrap()
                .message
                .contains("cancelled before publication")
        );
        assert!(!output.exists());
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn native_tss_collection_export_refuses_changed_destination_and_cleans_staging() {
        let engine = collection_engine(false);
        let dir = temp_root();
        let output = dir.path().join("appeared-during-export");
        let before = serde_json::to_value(engine.state()).unwrap();
        let mut created = false;
        let result = engine.export_native_tss_svg(&operation(&output), &mut |progress| {
            if let OperationProgress::Tfbs(p) = progress {
                if p.task_kind.as_deref() == Some("tss_svg_export")
                    && p.scanned_steps > 0
                    && !created
                {
                    fs::create_dir(&output).unwrap();
                    fs::write(output.join("sentinel"), b"do not overwrite").unwrap();
                    created = true;
                }
            }
            true
        });
        assert!(created);
        assert!(result.is_err());
        assert_eq!(
            fs::read(output.join("sentinel")).unwrap(),
            b"do not overwrite"
        );
        assert_eq!(fs::read_dir(&output).unwrap().count(), 1);
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    }

    #[test]
    fn native_tss_collection_export_checks_raw_paths_through_the_shared_shell() {
        use crate::engine_shell::{execute_shell_command, parse_shell_line, shell_quote};
        let mut engine = collection_engine(false);
        let dir = temp_root();
        fs::create_dir(dir.path().join("child")).unwrap();
        // Build the raw string, not a Windows verbatim PathBuf::join that can
        // normalize away the parent step before the validator receives it.
        let raw = format!(
            "{}{sep}child{sep}..{sep}refused",
            dir.path().display(),
            sep = std::path::MAIN_SEPARATOR
        );
        assert!(Path::new(&raw).components().any(|component| {
            matches!(component, std::path::Component::ParentDir)
                || matches!(component, std::path::Component::Normal(value) if value == std::ffi::OsStr::new(".."))
        }));
        let command = parse_shell_line(&format!(
            "promoters tss-view-svg --collection toy_tss {}",
            shell_quote(&raw)
        ))
        .unwrap();
        let error = execute_shell_command(&mut engine, &command).unwrap_err();
        assert!(error.contains("parent traversal"), "{error}");
        assert!(!dir.path().join("refused").exists());

        let output = dir.path().join("existing");
        fs::create_dir(&output).unwrap();
        assert!(engine.apply(operation(&output)).is_err());
        assert_eq!(fs::read_dir(&output).unwrap().count(), 0);
    }

    #[cfg(unix)]
    #[test]
    fn native_tss_collection_export_refuses_symlink_ancestors_before_staging() {
        let mut engine = collection_engine(false);
        let dir = temp_root();
        let target = dir.path().join("target");
        fs::create_dir(&target).unwrap();
        let alias = dir.path().join("alias");
        std::os::unix::fs::symlink(&target, &alias).unwrap();
        let error = engine.apply(operation(&alias.join("pages"))).unwrap_err();
        assert!(error.message.contains("symlinks"));
        assert_eq!(fs::read_dir(&target).unwrap().count(), 0);
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 2);
    }

    #[test]
    fn native_tss_collection_export_keeps_the_32_member_limit() {
        let mut engine = synthetic_tss_engine(false);
        for i in 0..33 {
            engine
                .state
                .sequences
                .get_mut("locus")
                .unwrap()
                .features_mut()
                .push(gb_io::seq::Feature {
                    kind: "mRNA".into(),
                    location: gb_io::seq::Location::simple_range(100 + i, 600),
                    qualifiers: vec![
                        ("gene".into(), Some("TOY".into())),
                        ("gene_id".into(), Some("gene_TOY".into())),
                        ("transcript_id".into(), Some(format!("additional_{i}"))),
                        ("source".into(), Some("synthetic".into())),
                    ],
                });
        }
        let request = synthetic_tss_approval(&engine);
        assert!(request.selected_tss_ids.len() > MAX_MEMBERS);
        engine
            .apply(Operation::MaterializeTssWindows { request })
            .unwrap();
        let dir = temp_root();
        let output = dir.path().join("too-many");
        assert!(
            engine
                .apply(operation(&output))
                .unwrap_err()
                .message
                .contains("at most 32")
        );
        assert!(!output.exists());
    }

    #[test]
    fn native_tss_export_index_escapes_labels_but_retains_svg_links() {
        let engine = collection_engine(false);
        let dir = temp_root();
        let mut receipt = engine
            .export_native_tss_svg(&operation(&dir.path().join("pages")), &mut |_| true)
            .unwrap()
            .receipt;
        receipt.pages[0].title = "<script>not executable</script> & DNA".into();
        let index = index_html(&receipt);
        assert!(index.contains("&lt;script&gt;not executable&lt;/script&gt; &amp; DNA"));
        assert!(!index.contains("<script>"));
        assert!(index.contains("href=\"001.svg\""));
    }
}

fn write_file(path: &Path, bytes: &[u8]) -> Result<(), EngineError> {
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(io_error)?;
    file.write_all(bytes).map_err(io_error)?;
    file.sync_all().map_err(io_error)
}

fn html_text(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

fn index_html(receipt: &TssViewSvgExportReceipt) -> String {
    let mut html = String::from(
        "<!doctype html>\n<html lang=\"en\"><meta charset=\"utf-8\"><title>Native TSS collection</title><body><h1>Native TSS collection</h1><p>Ordered SVG pages retain hover text. Scale policy: per window and per lane, preserving supplied report scales. No new cross-window calibration or cross-matrix ranking. Scores do not establish binding or promoter activity.</p><p><a href=\"receipt.json\">Content-bound receipt</a></p><ol>",
    );
    for page in &receipt.pages {
        html.push_str(&format!(
            "<li><a href=\"{}\">{}</a> | {} | {}:{}..{} ({}) | local {}..{}<br>SVG SHA-256 <code>{}</code></li>",
            page.filename, html_text(&page.title), html_text(&page.seq_id),
            html_text(&page.geometry.chromosome), page.geometry.start_1based,
            page.geometry.end_1based, page.geometry.strand.as_str(),
            page.start_0based + 1, page.end_0based_exclusive, page.svg_sha256
        ));
    }
    html.push_str("</ol></body></html>\n");
    html
}

impl GentleEngine {
    pub(super) fn export_native_tss_svg(
        &self,
        operation: &Operation,
        on_progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<NativeTssExportOutcome, EngineError> {
        let Operation::ExportTssViewSvg {
            seq_id,
            collection_id,
            path,
            report,
            local_motifs,
            score_kind,
            clip_negative,
            start_0based,
            end_0based_exclusive,
            width_px,
        } = operation
        else {
            return Err(EngineError::internal("Expected ExportTssViewSvg operation"));
        };
        if seq_id.is_empty() == collection_id.is_none() {
            return Err(EngineError::invalid_input(
                "TSS SVG export requires exactly one seq_id or collection_id",
            ));
        }
        let collection = collection_id
            .as_deref()
            .map(|id| self.get_tss_collection(id))
            .transpose()?;
        let ids = if let Some(collection) = &collection {
            if collection.members.len() > MAX_MEMBERS {
                return Err(EngineError::invalid_input(
                    "Native TSS collection export supports at most 32 validated windows",
                ));
            }
            collection
                .members
                .iter()
                .map(|member| member.tss.output_seq_id.clone())
                .collect::<Vec<_>>()
        } else {
            vec![seq_id.clone()]
        };
        let output = if collection.is_some() {
            // Reuse the raw-path, ancestor and Windows-verbatim traversal checks.
            let output = crate::tss_profile_export::destination(path)?;
            if fs::symlink_metadata(&output).is_ok() {
                return Err(EngineError::invalid_input(
                    "Native TSS collection export requires a new output directory",
                ));
            }
            Some(output)
        } else {
            None
        };
        let width = width_px.unwrap_or(1600);
        if !(1000..=5000).contains(&width) {
            return Err(EngineError::invalid_input(
                "TSS SVG width must be between 1000 and 5000 pixels",
            ));
        }
        let profile = report
            .as_deref()
            .map(|path| read_tss_profile_file(Path::new(path)))
            .transpose()
            .map_err(EngineError::invalid_input)?;
        let local_request = (!local_motifs.is_empty()).then(|| TssLocalScoreRequest {
            matrix_ids: local_motifs.clone(),
            score_kind: *score_kind,
            clip_negative: *clip_negative,
        });
        // Decode and validate every member and attachment before scoring any member.
        let mut views = Vec::with_capacity(ids.len());
        for id in &ids {
            let dna = self.state.sequences.get(id).ok_or_else(|| {
                EngineError::new(ErrorCode::NotFound, format!("Sequence '{id}' not found"))
            })?;
            let mut view = TssSequenceView::from_dna(dna).map_err(EngineError::invalid_input)?;
            let length = view
                .geometry
                .length()
                .ok_or_else(|| EngineError::invalid_input("Invalid TSS geometry"))?;
            let span = start_0based.unwrap_or(0)..end_0based_exclusive.unwrap_or(length);
            if span.start >= span.end || span.end > length {
                return Err(EngineError::invalid_input(format!(
                    "Invalid TSS SVG span within annotated window '{id}'"
                )));
            }
            if let Some((report, hash)) = &profile {
                view = view
                    .with_profile(report)
                    .map_err(EngineError::invalid_input)?;
                view.profile.as_mut().unwrap().file_sha256.clone_from(hash);
            }
            if let Some(request) = &local_request {
                request
                    .validate_budget(length)
                    .map_err(EngineError::invalid_input)?;
            }
            views.push((id.clone(), view, span));
        }
        let staging = output
            .as_ref()
            .map(|output| {
                tempfile::Builder::new()
                    .prefix(".gentle-tss-view-")
                    .tempdir_in(output.parent().unwrap())
                    .map_err(io_error)
            })
            .transpose()?;
        let mut outcome = NativeTssExportOutcome {
            receipt: TssViewSvgExportReceipt {
                schema: "gentle.tss_view_svg_export.v1".into(),
                request: serde_json::to_value(operation).map_err(io_error)?,
                request_sha256: digest(operation)?,
                exporter_revision: option_env!("GENTLE_SOURCE_REVISION")
                    .unwrap_or("unknown")
                    .into(),
                collection_sha256: collection.as_ref().map(digest).transpose()?,
                collection,
                scale_policy: SCALE_POLICY.into(),
                hover_limit_per_lane_strand: TSS_SVG_HOVER_LIMIT_PER_STRAND,
                pages: Vec::with_capacity(ids.len()),
                index_html_sha256: None,
            },
            messages: Vec::new(),
            warnings: Vec::new(),
        };
        let mut staged_bytes = 0u64;
        for (ordinal, (id, mut view, span)) in views.into_iter().enumerate() {
            checkpoint(on_progress, &id, ordinal, ids.len())?;
            if let Some(request) = &local_request {
                let sequence = self.state.sequences[&id].get_forward_string();
                let tracks = view
                    .compute_local_scores_in_span(&sequence, request, span.clone(), on_progress)
                    .map_err(EngineError::invalid_input)?;
                view = view
                    .with_local_scores_in_span(&sequence, request, &tracks, span.clone())
                    .map_err(EngineError::invalid_input)?;
                outcome.messages.push(format!(
                    "Computed local {} curves for {} exact matrix accession(s) in local {}..{}; only complete footprints within the span, excluding boundary-crossing windows; admission budget uses the full {}-bp annotated window ('{id}')",
                    score_kind.as_str(), local_motifs.len(), span.start + 1, span.end, view.geometry.length().unwrap()
                ));
            }
            let options = TssViewSvgOptions {
                start_0based: span.start,
                end_0based_exclusive: span.end,
                lane_indices: (0..view.lanes.len()).collect(),
                width_px: width,
                print_size_mm: None,
            };
            let view_sha256 = digest(&view)?;
            let filename = if staging.is_some() {
                format!("{:03}.svg", ordinal + 1)
            } else {
                path.clone()
            };
            let page_path = staging
                .as_ref()
                .map(|dir| dir.path().join(&filename))
                .unwrap_or_else(|| Path::new(path).to_path_buf());
            let svg_sha256 = write_tss_view_svg(&view, &options, &page_path)
                .map_err(EngineError::invalid_input)?;
            staged_bytes += fs::metadata(&page_path).map_err(io_error)?.len();
            if staging.is_some() && staged_bytes > MAX_BUNDLE_BYTES {
                return Err(EngineError::invalid_input(
                    "Native TSS collection export exceeds 32 MiB; use a smaller span or collection",
                ));
            }
            outcome.warnings.extend(view.warnings.iter().cloned());
            outcome.receipt.pages.push(TssViewSvgExportPage {
                seq_id: id.clone(),
                title: view.title.clone(),
                filename,
                svg_sha256: svg_sha256.clone(),
                view_sha256,
                sequence_sha256: view.sequence_sha256.clone(),
                geometry: view.geometry.clone(),
                start_0based: span.start,
                end_0based_exclusive: span.end,
                width_px: width,
                lane_count: view.lanes.len(),
                profile_file_sha256: view.profile.as_ref().map(|p| p.file_sha256.clone()),
                local_score_report_sha256: view
                    .local_scoring
                    .as_ref()
                    .map(|p| p.report_sha256.clone()),
            });
            outcome.messages.push(format!(
                "Wrote the native TSS view for '{id}' ({} lanes, local {}..{}) {} (SVG SHA-256 {svg_sha256})",
                view.lanes.len(), span.start + 1, span.end,
                if staging.is_some() { "as a staged collection page".into() } else { format!("to '{path}'") }
            ));
        }
        if let (Some(staging), Some(output)) = (staging, output) {
            let index = index_html(&outcome.receipt);
            outcome.receipt.index_html_sha256 = Some(sha256_hex_bytes(index.as_bytes()));
            let receipt_bytes = serde_json::to_vec_pretty(&outcome.receipt).map_err(io_error)?;
            if staged_bytes + index.len() as u64 + receipt_bytes.len() as u64 > MAX_BUNDLE_BYTES {
                return Err(EngineError::invalid_input(
                    "Native TSS collection export exceeds the 32 MiB bundle budget",
                ));
            }
            write_file(&staging.path().join("index.html"), index.as_bytes())?;
            write_file(&staging.path().join("receipt.json"), &receipt_bytes)?;
            checkpoint(
                on_progress,
                collection_id.as_deref().unwrap(),
                ids.len(),
                ids.len(),
            )?;
            // Check ancestors and freshness again; never remove an existing destination.
            let current_output = crate::tss_profile_export::destination(path)?;
            if current_output != output || fs::symlink_metadata(&output).is_ok() {
                return Err(EngineError::invalid_input(
                    "Native TSS collection destination changed before publication",
                ));
            }
            fs::rename(staging.path(), &output).map_err(io_error)?;
            let _ = staging.keep();
            outcome.messages.push(format!(
                "Published {} ordered TSS SVG pages, index.html and receipt.json to '{path}'; scale policy {SCALE_POLICY}; receipt SHA-256 {}",
                ids.len(), sha256_hex_bytes(&receipt_bytes)
            ));
        }
        if let Some(report) = report {
            outcome.messages.push(format!(
                "Attached validated TSS profile report '{report}'; no rescoring or database query"
            ));
        }
        outcome.messages.push(format!("SVG hover policy: at most {TSS_SVG_HOVER_LIMIT_PER_STRAND} titles per lane/strand; omissions disclosed without removing curve points. Highest raw positions are ranked only within each matrix lane."));
        Ok(outcome)
    }
}
