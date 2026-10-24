//! Bounded inspection, explicit cancellation and a synthetic rehearsal for every adapter.

use super::*;
use md5::{Digest, Md5};
use serde_json::Value;
use std::io::Read;

impl GentleEngine {
    /// Inspect or atomically save draft JSON without opening a reference, manifest or reads.
    pub fn rna_read_panel_screen_draft(
        &self,
        request: &ScreenRnaReadPanelRequest,
        output_path: Option<&str>,
    ) -> Result<RnaReadPanelScreenDraft, EngineError> {
        if request.schema != RNA_PANEL_SCREEN_REQUEST_SCHEMA {
            return Err(bad("panel_draft_schema"));
        }
        let artifact = if let Some(path) = output_path {
            let destination = Path::new(path);
            if path.trim().is_empty()
                || destination.exists()
                || destination.file_name().is_none()
                || destination.components().any(|c| c.as_os_str() == "..")
            {
                return Err(bad(
                    "panel_draft_output: choose a new file without parent traversal",
                ));
            }
            let parent = destination
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            let bytes = serde_json::to_vec_pretty(request).map_err(|e| bad(e.to_string()))?;
            let mut file =
                tempfile::NamedTempFile::new_in(parent).map_err(|e| bad(e.to_string()))?;
            file.write_all(&bytes).map_err(|e| bad(e.to_string()))?;
            file.as_file().sync_all().map_err(|e| bad(e.to_string()))?;
            file.persist_noclobber(destination)
                .map_err(|e| bad(format!("panel_draft_output: {e}")))?;
            Some(RnaReadPanelArtifact {
                path: path.into(),
                sha256: sha256_prefixed_bytes(&bytes),
                bytes: bytes.len() as u64,
                locus_id: None,
            })
        } else {
            None
        };
        Ok(RnaReadPanelScreenDraft {
            schema: "gentle.rna_read_panel_screen_draft.v1".into(),
            request: request.clone(),
            artifact,
        })
    }

    /// Create only synthetic files in a fresh directory; never replace the live project.
    pub fn prepare_rna_read_panel_demo(
        &self,
        output_dir: &str,
    ) -> Result<Vec<String>, EngineError> {
        let output = Path::new(output_dir);
        if output_dir.trim().is_empty() || output.exists() {
            return Err(bad("panel_demo_output_exists_or_empty"));
        }
        let parent = output
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let output = fs::canonicalize(parent)
            .map_err(|e| bad(e.to_string()))?
            .join(
                output
                    .file_name()
                    .ok_or_else(|| bad("panel_demo_output_name"))?,
            );
        let directory = output.join("artifacts");
        let stage = tempfile::tempdir_in(parent).map_err(|e| bad(e.to_string()))?;
        let mut state = ProjectState::default();
        let mut sources = vec![];
        let mut sequences = vec![];
        // Hand-crafted fixed-LCG DNA, not a human gene or a biological control.
        for i in 0..5 {
            let mut value = 100 + i as u64;
            let sequence: String = (0..150)
                .map(|_| {
                    value = value.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[((value >> 32) % 4) as usize] as char
                })
                .collect();
            let seq_id = format!("syn_{i}");
            let mut dna = DNAsequence::from_sequence(&sequence).map_err(|e| bad(e.to_string()))?;
            for (id, end) in [("reused_tx", 100), ("reused_alt", 150)] {
                dna.features_mut().push(gb_io::seq::Feature {
                    kind: "mRNA".into(),
                    location: gb_io::seq::Location::Join(vec![
                        gb_io::seq::Location::simple_range(0, 40),
                        gb_io::seq::Location::simple_range(end - 40, end),
                    ]),
                    qualifiers: vec![
                        ("gene".into(), Some(format!("SYN_{i}"))),
                        ("transcript_id".into(), Some(id.into())),
                    ],
                });
            }
            state.sequences.insert(seq_id.clone(), dna);
            sources.push(RnaReadPanelSource {
                locus_id: format!("locus_{i}"),
                seq_id,
                source_feature_id: 0,
                source: "hand_crafted".into(),
                assembly: "synthetic".into(),
                annotation_release: "v1".into(),
                gene_id: format!("SYN_{i}"),
                gene_symbol: format!("SYN_{i}"),
                selected_transcript_ids: vec![],
                expected_source_record_sha256: None,
            });
            sequences.push(sequence);
        }
        GentleEngine::from_state(state)
            .state()
            .save_to_path(&stage.path().join("starter.json").display().to_string())?;
        let build = BuildRnaReadPanelRequest {
            gene_set: None,
            register_in_project: None,
            base_panel_path: None,
            schema: RNA_PANEL_BUILD_SCHEMA.into(),
            sources,
            output_path: directory.join("panel.json").display().to_string(),
        };
        write_json(&stage.path().join("build.json"), &build)?;
        let mut files = vec![];
        let mut manifest = String::new();
        for number in 1..=2 {
            let name = format!("demo-{number}.fastq.gz");
            let reads = [
                sequences[0][..40].to_string(),
                format!("{}{}", &sequences[0][..40], &sequences[0][60..100]),
                format!("{}{}", &sequences[1][..40], &sequences[1][110..150]),
                "N".repeat(80),
            ];
            let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            for (i, sequence) in reads.iter().enumerate() {
                writeln!(
                    gzip,
                    "@demo_{i} synthetic_library={number}\n{sequence}\n+\n{}",
                    ">".repeat(sequence.len())
                )
                .map_err(|e| bad(e.to_string()))?;
            }
            let bytes = gzip.finish().map_err(|e| bad(e.to_string()))?;
            fs::write(stage.path().join(&name), &bytes).map_err(|e| bad(e.to_string()))?;
            manifest.push_str(&format!(
                "{}  /synthetic/facility/{name}\n",
                input::hex(&Md5::digest(&bytes))
            ));
            files.push(RnaReadPanelInput {
                sample_id: format!("DEMO_{number}"),
                barcode: format!("demo_barcode_{number}"),
                priming: "synthetic_custom".into(),
                path: directory.join(name).display().to_string(),
                declared_complete: true,
                expected_md5: None,
                expected_sha256: None,
                checksum_manifest_key: None,
            });
        }
        fs::write(stage.path().join("md5_checksums.txt"), manifest)
            .map_err(|e| bad(e.to_string()))?;
        let screen = ScreenRnaReadPanelRequest {
            expected_panel_sha256: None,
            schema: RNA_PANEL_SCREEN_REQUEST_SCHEMA.into(),
            panel_path: build.output_path.clone(),
            files,
            md5_manifest_path: Some(directory.join("md5_checksums.txt").display().to_string()),
            output_dir: output.join("screen").display().to_string(),
            policy: RnaReadPanelPolicy {
                kmer_len: 12,
                min_chain_bases: 24,
                min_junction_flank_bp: 4,
                write_locus_fastas: true,
                ..Default::default()
            },
            cancel_path: Some(output.join("CANCEL").display().to_string()),
        };
        write_json(&stage.path().join("screen.json"), &screen)?;
        write_json(
            &stage.path().join("summary.json"),
            &ExportRnaReadPanelSummaryRequest {
                schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
                report_paths: vec![output.join("screen/report.json").display().to_string()],
                output_dir: output.join("summary").display().to_string(),
            },
        )?;
        fs::create_dir(&output).map_err(|e| bad(e.to_string()))?;
        if let Err(error) = fs::rename(stage.path(), &directory) {
            let _ = fs::remove_dir(&output);
            return Err(bad(error.to_string()));
        }
        Ok([
            "starter.json",
            "build.json",
            "screen.json",
            "summary.json",
            "md5_checksums.txt",
        ]
        .iter()
        .map(|name| directory.join(name).display().to_string())
        .collect())
    }

    /// Request cancellation without altering inputs, existing markers or project state.
    pub fn cancel_rna_read_panel(
        &self,
        request: &ScreenRnaReadPanelRequest,
    ) -> Result<String, EngineError> {
        if request.schema != RNA_PANEL_SCREEN_REQUEST_SCHEMA {
            return Err(bad("panel_cancel_schema"));
        }
        let path = request
            .cancel_path
            .as_deref()
            .filter(|p| !p.trim().is_empty())
            .ok_or_else(|| bad("panel_cancel_path_required"))?;
        let identity = control::cancel_identity(request)?;
        if control::marker_requested(Some(path), std::slice::from_ref(&identity))? {
            return Ok(
                "Cancellation marker already exists for this exact request; no bytes changed"
                    .into(),
            );
        }
        let parent = Path::new(path)
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        // Publish the whole identity atomically; a running reader must never see an empty marker.
        let mut marker = tempfile::NamedTempFile::new_in(parent)
            .map_err(|e| bad(format!("panel_cancel_marker: {e}")))?;
        marker
            .write_all(identity.as_bytes())
            .map_err(|e| bad(e.to_string()))?;
        marker
            .as_file()
            .sync_all()
            .map_err(|e| bad(e.to_string()))?;
        match marker.persist_noclobber(path) {
            Ok(_) => {}
            Err(error) if error.error.kind() == std::io::ErrorKind::AlreadyExists => {
                if !control::marker_requested(Some(path), std::slice::from_ref(&identity))? {
                    return Err(bad(
                        "panel_cancel_marker_changed: retry the explicit request",
                    ));
                }
                return Ok(
                    "Cancellation marker already exists for this exact request; no bytes changed"
                        .into(),
                );
            }
            Err(error) => return Err(bad(format!("panel_cancel_marker: {error}"))),
        }
        Ok(
            "Cancellation requested, not confirmed terminated; inspect the eventual report status"
                .into(),
        )
    }

    /// Inspect content-bound metadata and optionally a bounded, explicitly unverified TSV prefix.
    pub fn inspect_rna_read_panel(
        &self,
        path: &str,
        file_id: Option<&str>,
        rows: usize,
    ) -> Result<
        (
            Option<RnaReadPanel>,
            Option<RnaReadPanelScreenReport>,
            Option<RnaReadPanelInspection>,
            Option<RnaReadPanelSummary>,
        ),
        EngineError,
    > {
        if rows > 1000 {
            return Err(bad("panel_inspection_limit: at most 1000 rows"));
        }
        if rows == 0 && file_id.is_some() {
            return Err(bad(
                "panel_inspection_file_requires_rows: select a positive --rows count",
            ));
        }
        let value: Value = load(path)?;
        if value.get("schema").and_then(Value::as_str) == Some(RNA_PANEL_SCHEMA) {
            if file_id.is_some() || rows != 0 {
                return Err(bad("panel_inspection_rows_require_screen_report"));
            }
            let panel: RnaReadPanel =
                serde_json::from_value(value).map_err(|e| bad(e.to_string()))?;
            validate_panel(&panel)?;
            return Ok((Some(panel), None, None, None));
        }
        if matches!(
            value.get("schema").and_then(Value::as_str),
            Some(RNA_PANEL_SUMMARY_SCHEMA | RNA_PANEL_SUMMARY_LEGACY_SCHEMA)
        ) {
            if file_id.is_some() || rows != 0 {
                return Err(bad("panel_inspection_rows_require_screen_report"));
            }
            let summary: RnaReadPanelSummary =
                serde_json::from_value(value).map_err(|e| bad(e.to_string()))?;
            let mut basis = summary.clone();
            basis.content_sha256.clear();
            if summary.content_sha256 != sha256_prefixed_bytes(&encode(&basis)?)
                || !matches!(summary.completeness.as_str(), "complete" | "partial")
                || summary.original_inputs_reverified
                || summary.nomination_sidecars_reverified
                || (summary.schema == RNA_PANEL_SUMMARY_SCHEMA
                    && (summary.ambiguous_sets.is_some()
                        || summary.co_nomination_sets.is_none()
                        || summary.sample_coverage.is_none()
                        || summary.class_member_counts_complete.is_none()
                        || summary.samples.iter().any(|s| {
                            s.co_nominated_reads.is_none() || s.ambiguous_reads.is_some()
                        })
                        || summary.transcripts.iter().any(|s| {
                            s.co_nominated_reads.is_none() || s.ambiguous_compatible_reads.is_some()
                        })))
            {
                return Err(bad("panel_summary_binding"));
            }
            return Ok((None, None, None, Some(summary)));
        }
        let report: RnaReadPanelScreenReport =
            serde_json::from_value(value).map_err(|e| bad(e.to_string()))?;
        validate_screen_report_binding(&report)?;
        if rows == 0 {
            return Ok((None, Some(report), None, None));
        }
        let receipt = file_id
            .and_then(|id| report.files_completed.iter().find(|r| r.file_id == id))
            .or_else(|| {
                if file_id.is_none() {
                    report.files_completed.first()
                } else {
                    None
                }
            })
            .ok_or_else(|| bad("panel_inspection_file_missing"))?;
        let ordinal = report
            .files_completed
            .iter()
            .position(|r| r.file_id == receipt.file_id)
            .unwrap();
        if receipt.file_id != format!("file_{ordinal:04}") {
            return Err(bad("panel_inspection_file_identity"));
        }
        let relative = format!("artifacts/{}/nominations.tsv", receipt.file_id);
        let artifact = receipt
            .artifacts
            .iter()
            .find(|a| a.path == relative)
            .ok_or_else(|| bad("panel_inspection_artifact_missing"))?;
        let root = Path::new(path).parent().unwrap_or(Path::new("."));
        let file = fs::File::open(root.join(relative)).map_err(|e| bad(e.to_string()))?;
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(file.take(8_000_000));
        let mut inspected = vec![];
        for row in reader.records().take(rows) {
            let row = row.map_err(|e| bad(format!("panel_inspection_tsv: {e}")))?;
            if row.len() != 9 {
                return Err(bad("panel_inspection_tsv_columns"));
            }
            inspected.push(RnaReadPanelNomination {
                sample_id: row[0].into(),
                barcode: row[1].into(),
                file_id: row[2].into(),
                expected_input_digest: row[3].into(),
                record_index: row[4]
                    .parse()
                    .map_err(|_| bad("panel_inspection_record_index"))?,
                read_id: row[5].into(),
                decompressed_record_offset: row[6]
                    .parse()
                    .map_err(|_| bad("panel_inspection_offset"))?,
                class: row[7].into(),
                candidates: serde_json::from_str(&row[8])
                    .map_err(|e| bad(format!("panel_inspection_candidates: {e}")))?,
            });
        }
        let inspection = RnaReadPanelInspection {
            schema: "gentle.rna_read_panel_inspection.v1".into(),
            report_sha256: report.content_sha256.clone(),
            file_id: receipt.file_id.clone(),
            expected_artifact_sha256: artifact.sha256.clone(),
            artifact_content_verified: false,
            rows: inspected,
        };
        Ok((None, Some(report), Some(inspection), None))
    }
}
