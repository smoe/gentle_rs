//! Bounded, report-only read-support aggregation with ambiguity and input identity preserved.

use super::*;

const MAX_METADATA_BYTES: u64 = 256_000_000;
const MAX_TRANSCRIPT_ROWS: usize = 250_000;
const MAX_AMBIGUOUS_SETS: usize = 65_536;
const MAX_MEMBER_VISITS: usize = 2_000_000;

fn add(total: &mut u64, value: u64) -> Result<(), EngineError> {
    *total = total
        .checked_add(value)
        .ok_or_else(|| bad("panel_summary_count_overflow"))?;
    Ok(())
}

fn bounded_metadata<T: serde::de::DeserializeOwned>(
    path: &Path,
    remaining: &mut u64,
) -> Result<T, EngineError> {
    let bytes = fs::metadata(path)
        .map_err(|e| bad(format!("panel_summary_metadata: {e}")))?
        .len();
    *remaining = remaining.checked_sub(bytes).ok_or_else(|| {
        bad("panel_summary_metadata_budget: at most 256 MB of report/panel metadata")
    })?;
    load(&path.to_string_lossy())
}

fn validate_run(report: &RnaReadPanelScreenReport) -> Result<(), EngineError> {
    validate_screen_report_binding(report)?;
    let p = &report.request.policy;
    if report.request.schema != RNA_PANEL_SCREEN_REQUEST_SCHEMA
        || report.request.files.is_empty()
        || report.request.files.len() > 4096
        || report.files_completed.len() > report.request.files.len()
        || !(8..=16).contains(&p.kmer_len)
        || p.min_chain_bases < p.kmer_len as u32
        || p.min_junction_flank_bp == 0
        || p.min_junction_flank_bp > p.kmer_len / 2
        || !(1..=4_000_000).contains(&p.max_index_postings)
        || !(1..=250_000).contains(&p.max_anchors_per_read)
        || !(1..=16_384).contains(&p.max_equivalence_classes)
        || report.indexed_postings == 0
        || report.indexed_postings > p.max_index_postings
        || (report.status == "completed"
            && report.files_completed.len() != report.request.files.len())
        || report.files_not_completed
            != report.request.files[report.files_completed.len()..]
                .iter()
                .map(|input| input.path.clone())
                .collect::<Vec<_>>()
    {
        return Err(bad("panel_summary_run_geometry_or_policy"));
    }
    for (ordinal, file) in report.files_completed.iter().enumerate() {
        if file.file_id != format!("file_{ordinal:04}")
            || encode(&file.input)? != encode(&report.request.files[ordinal])?
            || !safe_id(&file.input.sample_id)
            || file.input.barcode.trim().is_empty()
            || file.input.priming.trim().is_empty()
            || !file.input.declared_complete
            || !valid_digest(&file.observed_md5, 32)
            || !file.observed_sha256.starts_with("sha256:")
            || !valid_digest(&file.observed_sha256[7..], 64)
            || (file.input.expected_md5.is_none() && file.input.expected_sha256.is_none())
            || file
                .input
                .expected_md5
                .as_ref()
                .is_some_and(|expected| expected != &file.observed_md5)
            || file.input.expected_sha256.as_ref().is_some_and(|expected| {
                expected.trim_start_matches("sha256:").to_ascii_lowercase()
                    != file.observed_sha256[7..]
            })
            || file.equivalence_classes.len() > p.max_equivalence_classes
        {
            return Err(bad("panel_summary_file_binding"));
        }
        let mut class_reads = 0;
        for (class, count) in &file.classes {
            if !matches!(
                class.as_str(),
                "single_locus_candidate"
                    | "multi_locus_ambiguous"
                    | "multi_segment_candidate"
                    | "partial"
                    | "unmatched"
            ) || *count == 0
            {
                return Err(bad("panel_summary_unknown_or_empty_class"));
            }
            add(&mut class_reads, *count)?;
        }
        if class_reads != file.reads {
            return Err(bad("panel_summary_class_total"));
        }
    }
    Ok(())
}

fn write_tsv<T: Serialize>(
    directory: &Path,
    name: &str,
    header: &[&str],
    rows: impl IntoIterator<Item = T>,
) -> Result<RnaReadPanelArtifact, EngineError> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .terminator(csv::Terminator::Any(b'\n'))
        .from_writer(HashedWriter::create(&directory.join(name))?);
    writer
        .write_record(header)
        .map_err(|e| bad(e.to_string()))?;
    for row in rows {
        writer
            .serialize(row)
            .map_err(|e| bad(format!("panel_summary_tsv: {e}")))?;
    }
    writer
        .into_inner()
        .map_err(|e| bad(e.to_string()))?
        .receipt(format!("artifacts/{name}"))
}

impl GentleEngine {
    /// Export only stored read support. Repeated compressed-input hashes are refused, never pooled twice.
    pub fn export_rna_read_panel_summary(
        &self,
        request: &ExportRnaReadPanelSummaryRequest,
    ) -> Result<RnaReadPanelSummary, EngineError> {
        if request.schema != RNA_PANEL_SUMMARY_REQUEST_SCHEMA
            || !(1..=128).contains(&request.report_paths.len())
            || request.report_paths.iter().any(|p| p.trim().is_empty())
        {
            return Err(bad(
                "panel_summary_request: expected summary v1 and 1..=128 explicit reports",
            ));
        }
        let output = Path::new(&request.output_dir);
        if request.output_dir.trim().is_empty()
            || output.exists()
            || output.file_name().is_none()
            || output
                .components()
                .any(|c| c == std::path::Component::ParentDir || c.as_os_str() == "..")
        {
            return Err(bad(
                "panel_summary_output: choose a fresh directory without parent traversal",
            ));
        }
        let mut remaining = MAX_METADATA_BYTES;
        let mut reports = vec![];
        let mut identities = BTreeSet::new();
        for path in &request.report_paths {
            let report: RnaReadPanelScreenReport =
                bounded_metadata(Path::new(path), &mut remaining)?;
            validate_run(&report)?;
            if !identities.insert(report.content_sha256.clone()) {
                return Err(bad("panel_summary_duplicate_report"));
            }
            reports.push((path.clone(), report));
        }
        reports.sort_by(|a, b| a.1.content_sha256.cmp(&b.1.content_sha256));
        let mut panel: Option<RnaReadPanel> = None;
        let policy = reports[0].1.request.policy.clone();
        let mut runs = vec![];
        let mut files = vec![];
        let mut input_hashes = BTreeSet::new();
        let mut samples: BTreeMap<String, RnaReadPanelSummarySample> = BTreeMap::new();
        let mut transcripts: BTreeMap<(String, String), RnaReadPanelSummaryTranscript> =
            BTreeMap::new();
        let mut sets: BTreeMap<(String, Vec<String>), u64> = BTreeMap::new();
        let mut member_visits = 0usize;
        let mut completeness = "complete";
        for (path, report) in reports {
            let retained_panel: RnaReadPanel = bounded_metadata(
                &Path::new(&path)
                    .parent()
                    .unwrap_or(Path::new("."))
                    .join("artifacts/panel.json"),
                &mut remaining,
            )?;
            validate_panel(&retained_panel)?;
            if retained_panel.content_sha256 != report.panel_sha256
                || panel
                    .as_ref()
                    .is_some_and(|p| p.content_sha256 != report.panel_sha256)
                || encode(&policy)? != encode(&report.request.policy)?
                || !retained_panel.overlapping_source_pairs.is_empty()
            {
                return Err(bad("panel_summary_incompatible_panel_or_policy"));
            }
            if panel.is_none() {
                panel = Some(retained_panel);
            }
            let reference = panel.as_ref().unwrap();
            let members: BTreeMap<_, _> = reference
                .templates
                .iter()
                .map(|t| (t.member_id.as_str(), t))
                .collect();
            if report.status != "completed" || !report.files_not_completed.is_empty() {
                completeness = "partial";
            }
            runs.push(RnaReadPanelSummaryRun {
                report_path: path,
                report_sha256: report.content_sha256.clone(),
                status: report.status,
                completed_files: report.files_completed.len(),
                files_not_completed: report.files_not_completed,
            });
            for file in report.files_completed {
                if files.len() == 4096 {
                    return Err(bad(
                        "panel_summary_file_budget: at most 4096 verified files",
                    ));
                }
                if !input_hashes.insert(file.observed_sha256.to_ascii_lowercase()) {
                    return Err(bad(
                        "panel_summary_duplicate_input: select each compressed-input SHA-256 once, even across samples/runs",
                    ));
                }
                let sample_id = file.input.sample_id.clone();
                if !samples.contains_key(&sample_id) {
                    if samples.len() == 128
                        || (samples.len() + 1) * reference.templates.len() > MAX_TRANSCRIPT_ROWS
                    {
                        return Err(bad("panel_summary_sample_or_row_budget"));
                    }
                    samples.insert(
                        sample_id.clone(),
                        RnaReadPanelSummarySample {
                            sample_id: sample_id.clone(),
                            barcode: file.input.barcode.clone(),
                            priming: file.input.priming.clone(),
                            verified_files: 0,
                            reads: 0,
                            single_candidate_reads: 0,
                            ambiguous_reads: 0,
                            no_candidate_reads: 0,
                            classes: BTreeMap::new(),
                        },
                    );
                    for template in &reference.templates {
                        let source = &reference
                            .loci
                            .iter()
                            .find(|l| l.source.locus_id == template.locus_id)
                            .unwrap()
                            .source;
                        transcripts.insert(
                            (sample_id.clone(), template.member_id.clone()),
                            RnaReadPanelSummaryTranscript {
                                sample_id: sample_id.clone(),
                                gene_id: source.gene_id.clone(),
                                gene_symbol: source.gene_symbol.clone(),
                                locus_id: template.locus_id.clone(),
                                member_id: template.member_id.clone(),
                                transcript_id: template.transcript_id.clone(),
                                single_candidate_reads: 0,
                                ambiguous_compatible_reads: 0,
                            },
                        );
                    }
                }
                let sample = samples.get_mut(&sample_id).unwrap();
                if sample.barcode != file.input.barcode || sample.priming != file.input.priming {
                    return Err(bad("panel_summary_sample_binding_conflict"));
                }
                sample.verified_files += 1;
                add(&mut sample.reads, file.reads)?;
                for (class, count) in &file.classes {
                    add(sample.classes.entry(class.clone()).or_default(), *count)?;
                }
                let mut counted = 0;
                let mut empty = 0;
                for (key, count) in file.equivalence_classes {
                    let member_ids: Vec<String> = serde_json::from_str(&key)
                        .map_err(|e| bad(format!("panel_summary_member_set: {e}")))?;
                    if count == 0
                        || member_ids.windows(2).any(|w| w[0] >= w[1])
                        || serde_json::to_string(&member_ids).map_err(|e| bad(e.to_string()))?
                            != key
                        || member_ids
                            .iter()
                            .any(|id| !members.contains_key(id.as_str()))
                    {
                        return Err(bad("panel_summary_unknown_unsorted_or_duplicate_member"));
                    }
                    member_visits += member_ids.len();
                    if member_visits > MAX_MEMBER_VISITS {
                        return Err(bad("panel_summary_member_budget"));
                    }
                    add(&mut counted, count)?;
                    match member_ids.len() {
                        0 => {
                            add(&mut sample.no_candidate_reads, count)?;
                            add(&mut empty, count)?;
                        }
                        1 => {
                            add(&mut sample.single_candidate_reads, count)?;
                        }
                        _ => {
                            add(&mut sample.ambiguous_reads, count)?;
                            let key = (sample_id.clone(), member_ids.clone());
                            if !sets.contains_key(&key) && sets.len() == MAX_AMBIGUOUS_SETS {
                                return Err(bad("panel_summary_set_budget"));
                            }
                            add(sets.entry(key).or_default(), count)?;
                        }
                    }
                    for id in &member_ids {
                        let row = transcripts
                            .get_mut(&(sample_id.clone(), id.clone()))
                            .unwrap();
                        if member_ids.len() == 1 {
                            add(&mut row.single_candidate_reads, count)?;
                        } else {
                            add(&mut row.ambiguous_compatible_reads, count)?;
                        }
                    }
                }
                let mut expected_empty = *file.classes.get("partial").unwrap_or(&0);
                add(
                    &mut expected_empty,
                    *file.classes.get("unmatched").unwrap_or(&0),
                )?;
                if counted != file.reads || empty != expected_empty {
                    return Err(bad("panel_summary_equivalence_total"));
                }
                files.push(RnaReadPanelSummaryFile {
                    report_sha256: report.content_sha256.clone(),
                    file_id: file.file_id,
                    sample_id,
                    input_path: file.input.path,
                    observed_sha256: file.observed_sha256,
                    observed_md5: file.observed_md5,
                    reads: file.reads,
                });
            }
        }
        let mut summary = RnaReadPanelSummary {
            schema: RNA_PANEL_SUMMARY_SCHEMA.into(), content_sha256: String::new(), panel_sha256: panel.unwrap().content_sha256,
            policy, completeness: completeness.into(), source_runs: runs, source_files: files,
            samples: samples.into_values().collect(), transcripts: transcripts.into_values().collect(),
            ambiguous_sets: sets.into_iter().map(|((sample_id, member_ids), reads)|
                RnaReadPanelSummaryAmbiguousSet { sample_id, member_ids, reads }).collect(),
            original_inputs_reverified: false, nomination_sidecars_reverified: false, artifacts: vec![], nonclaims: vec![
                "Counts are supporting reads, not original molecule copy numbers, expression or exclusive isoform identification; no alignment confirmation is added.".into(),
                "Single-candidate means one compatible panel member under the stored nomination policy, not transcriptome-wide uniqueness.".into(),
                "Ambiguous sets count each read once; per-transcript ambiguous_compatible_reads overlap and must not be summed as read totals. No fractional redistribution or winner is inferred.".into(),
                "Zero denotes no stored support among verified files, not biological absence. Completeness covers only the explicit report selection; cancelled reports keep this summary partial.".into(),
                "Only stored report/panel content bindings and count consistency are checked. Original FASTQ and nomination sidecars are not reopened or reverified; content hashes are not authenticated provenance.".into(),
                "Duplicate compressed-input hashes are refused. Distinct files can still contain duplicate reads; no read-level or PCR/UMI deduplication is performed.".into(),
            ],
        };
        let parent = output
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let staging = tempfile::Builder::new()
            .prefix(".gentle-rna-summary-")
            .tempdir_in(parent)
            .map_err(|e| bad(e.to_string()))?;
        summary.artifacts.push(write_tsv(
            staging.path(),
            "samples.tsv",
            &[
                "sample_id",
                "barcode",
                "priming",
                "verified_files",
                "reads",
                "single_candidate_reads",
                "ambiguous_reads",
                "no_candidate_reads",
                "classes_json",
            ],
            summary
                .samples
                .iter()
                .map(|s| {
                    Ok::<_, EngineError>((
                        &s.sample_id,
                        &s.barcode,
                        &s.priming,
                        s.verified_files,
                        s.reads,
                        s.single_candidate_reads,
                        s.ambiguous_reads,
                        s.no_candidate_reads,
                        serde_json::to_string(&s.classes).map_err(|e| bad(e.to_string()))?,
                    ))
                })
                .collect::<Result<Vec<_>, _>>()?,
        )?);
        summary.artifacts.push(write_tsv(
            staging.path(),
            "transcripts.tsv",
            &[
                "sample_id",
                "gene_id",
                "gene_symbol",
                "locus_id",
                "member_id",
                "transcript_id",
                "single_candidate_reads",
                "ambiguous_compatible_reads",
            ],
            summary.transcripts.iter().map(|r| {
                (
                    &r.sample_id,
                    &r.gene_id,
                    &r.gene_symbol,
                    &r.locus_id,
                    &r.member_id,
                    &r.transcript_id,
                    r.single_candidate_reads,
                    r.ambiguous_compatible_reads,
                )
            }),
        )?);
        summary.artifacts.push(write_tsv(
            staging.path(),
            "ambiguous_sets.tsv",
            &["sample_id", "member_ids_json", "reads"],
            summary
                .ambiguous_sets
                .iter()
                .map(|r| {
                    Ok::<_, EngineError>((
                        &r.sample_id,
                        serde_json::to_string(&r.member_ids).map_err(|e| bad(e.to_string()))?,
                        r.reads,
                    ))
                })
                .collect::<Result<Vec<_>, _>>()?,
        )?);
        let basis = encode(&summary)?;
        if basis.len() > 128_000_000 {
            return Err(bad("panel_summary_json_budget: at most 128 MB"));
        }
        summary.content_sha256 = sha256_prefixed_bytes(&basis);
        fs::create_dir(output).map_err(|e| bad(format!("panel_summary_output_reserve: {e}")))?;
        if let Err(error) = fs::rename(staging.path(), output.join("artifacts")) {
            let _ = fs::remove_dir(output);
            return Err(bad(format!("panel_summary_publish: {error}")));
        }
        let mut final_file =
            tempfile::NamedTempFile::new_in(output).map_err(|e| bad(e.to_string()))?;
        serde_json::to_writer_pretty(final_file.as_file_mut(), &summary)
            .map_err(|e| bad(e.to_string()))?;
        final_file
            .as_file()
            .sync_all()
            .map_err(|e| bad(e.to_string()))?;
        final_file
            .persist_noclobber(output.join("report.json"))
            .map_err(|e| bad(e.to_string()))?;
        Ok(summary)
    }
}
