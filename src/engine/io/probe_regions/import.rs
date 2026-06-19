use super::*;

impl GentleEngine {
    pub fn import_apt_probe_region_output(
        &self,
        apt_summary_path: &str,
        annotation_path: &str,
        output_dir: &str,
        metadata_path: Option<&str>,
        condition_column: Option<&str>,
        sample_column: Option<&str>,
        probe_intensity_path: Option<&str>,
        probe_id_column: Option<&str>,
        platform: Option<&str>,
        normalization: Option<&str>,
        coordinate_system: Option<&str>,
        genome_build: Option<&str>,
    ) -> Result<ProbeRegionAptImportReport, EngineError> {
        let apt_summary_path = apt_summary_path.trim();
        let annotation_path = annotation_path.trim();
        let output_dir = output_dir.trim();
        if apt_summary_path.is_empty() || annotation_path.is_empty() || output_dir.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "APT probe-region import requires SUMMARY.tsv ANNOTATION.csv OUTPUT_DIR"
                    .to_string(),
                cause_chain: vec![],
            });
        }

        let platform = platform
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("Affymetrix APT")
            .to_string();
        let normalization = normalization
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("rma-sketch")
            .to_string();
        let coordinate_system = coordinate_system
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("not_declared")
            .to_string();
        let genome_build = genome_build
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("not_declared")
            .to_string();

        let annotation = Self::probe_region_apt_annotation_rows(annotation_path)?;
        let annotation_row_count = annotation.regions.len();
        let mut warnings = Vec::new();
        let mut summary_reader = csv::ReaderBuilder::new()
            .delimiter(Self::probe_region_metadata_delimiter(apt_summary_path))
            .trim(csv::Trim::All)
            .flexible(true)
            .from_path(apt_summary_path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not open APT summary '{apt_summary_path}': {e}"),
                cause_chain: vec![],
            })?;
        let summary_headers = summary_reader
            .headers()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not read APT summary header '{apt_summary_path}': {e}"),
                cause_chain: vec![],
            })?
            .iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        let feature_idx = Self::probe_region_metadata_column_index(
            &summary_headers,
            None,
            &[
                "probeset_id",
                "probeset",
                "probeset id",
                "feature_id",
                "probeset_or_region_id",
                "id",
            ],
        )
        .unwrap_or(0);
        let sample_indices = summary_headers
            .iter()
            .enumerate()
            .filter(|(idx, header)| *idx != feature_idx && !header.trim().is_empty())
            .map(|(idx, _)| idx)
            .collect::<Vec<_>>();
        let sample_columns = sample_indices
            .iter()
            .filter_map(|idx| summary_headers.get(*idx))
            .map(|header| header.trim().to_string())
            .collect::<Vec<_>>();
        if sample_columns.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("APT summary '{apt_summary_path}' has no sample columns"),
                cause_chain: vec![],
            });
        }
        let metadata_summary = if let Some(path) = metadata_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            let summary = Self::probe_region_apt_metadata_summary(
                path,
                condition_column,
                sample_column,
                &sample_columns,
            )?;
            warnings.extend(summary.warnings.iter().cloned());
            Some(summary)
        } else {
            None
        };
        let probe_intensity = if let Some(path) = probe_intensity_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            let table = Self::probe_region_apt_probe_intensity_rows(
                path,
                probe_id_column,
                metadata_path,
                condition_column,
                sample_column,
            )?;
            warnings.extend(table.warnings.iter().cloned());
            Some(table)
        } else {
            None
        };

        let mut output_rows: Vec<Vec<String>> = Vec::new();
        let mut output_value_rows: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut summary_row_count = 0usize;
        let mut missing_annotation_count = 0usize;
        let mut skipped_invalid_count = 0usize;
        for record in summary_reader.records() {
            let record = record.map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not read APT summary row: {e}"),
                cause_chain: vec![],
            })?;
            summary_row_count += 1;
            let feature_id = record.get(feature_idx).unwrap_or("").trim();
            if feature_id.is_empty() {
                skipped_invalid_count += 1;
                continue;
            }
            let Some(annotation_row) = annotation.regions.get(feature_id) else {
                missing_annotation_count += 1;
                continue;
            };
            let mut row = vec![
                annotation_row.chromosome.clone(),
                annotation_row.start_1based.to_string(),
                annotation_row.end_1based.to_string(),
                annotation_row.strand.clone(),
                annotation_row.feature_id.clone(),
                annotation_row.transcript_cluster_id.clone(),
                annotation_row.number_of_probes.clone(),
                annotation_row.gene_symbol.clone(),
            ];
            let sample_values = sample_indices
                .iter()
                .map(|idx| record.get(*idx).unwrap_or("").trim().to_string())
                .collect::<Vec<_>>();
            let mut value_row = sample_values.clone();
            for value in &sample_values {
                row.push(value.clone());
            }
            if let Some(metadata) = metadata_summary.as_ref() {
                let stats = metadata
                    .conditions
                    .iter()
                    .map(|condition| {
                        Self::probe_region_apt_condition_stats(
                            &sample_values,
                            &condition.sample_indices,
                        )
                    })
                    .collect::<Vec<_>>();
                for stat in &stats {
                    let mean = stat
                        .map(|(mean, _)| Self::probe_region_format_number(mean))
                        .unwrap_or_default();
                    let sd = stat
                        .map(|(_, sd)| Self::probe_region_format_number(sd))
                        .unwrap_or_default();
                    row.push(mean.clone());
                    row.push(sd.clone());
                    value_row.push(mean);
                    value_row.push(sd);
                }
                for contrast in &metadata.contrasts {
                    let numerator = stats
                        .get(contrast.numerator_index)
                        .and_then(|value| value.map(|(mean, _)| mean));
                    let denominator = stats
                        .get(contrast.denominator_index)
                        .and_then(|value| value.map(|(mean, _)| mean));
                    let logfc = numerator
                        .zip(denominator)
                        .map(|(num, den)| Self::probe_region_format_number(num - den))
                        .unwrap_or_default();
                    row.push(logfc.clone());
                    value_row.push(logfc);
                }
            }
            output_value_rows.insert(feature_id.to_string(), value_row);
            output_rows.push(row);
        }
        output_rows.sort_by(|a, b| {
            a[0].cmp(&b[0])
                .then_with(|| {
                    a[1].parse::<usize>()
                        .unwrap_or(usize::MAX)
                        .cmp(&b[1].parse::<usize>().unwrap_or(usize::MAX))
                })
                .then_with(|| a[4].cmp(&b[4]))
        });

        if missing_annotation_count > 0 {
            warnings.push(format!(
                "{} APT summary row(s) lacked matching annotation rows",
                missing_annotation_count
            ));
        }
        if skipped_invalid_count > 0 {
            warnings.push(format!(
                "{} APT summary row(s) were skipped because the probeset/region id was empty",
                skipped_invalid_count
            ));
        }

        std::fs::create_dir_all(output_dir).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not create APT probe-region output dir '{output_dir}': {e}"),
            cause_chain: vec![],
        })?;
        let region_table_path = Path::new(output_dir).join(PROBE_REGION_TABLE_FILE);
        let mut writer = csv::Writer::from_path(&region_table_path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create probe-region table '{}': {e}",
                region_table_path.to_string_lossy()
            ),
            cause_chain: vec![],
        })?;
        let mut value_header = sample_columns.clone();
        let mut header = vec![
            "chromosome".to_string(),
            "start".to_string(),
            "stop".to_string(),
            "strand".to_string(),
            "probeset_or_region_id".to_string(),
            "transcript_cluster_id".to_string(),
            "number_of_probes".to_string(),
            "gene_symbol".to_string(),
        ];
        header.extend(sample_columns.iter().cloned());
        if let Some(metadata) = metadata_summary.as_ref() {
            for condition in &metadata.conditions {
                let mean_header = format!("mean_log2_{}", condition.column_label);
                let sd_header = format!("sd_log2_{}", condition.column_label);
                header.push(mean_header.clone());
                header.push(sd_header.clone());
                value_header.push(mean_header);
                value_header.push(sd_header);
            }
            for contrast in &metadata.contrasts {
                let logfc_header = format!("log2FC_{}", contrast.label);
                header.push(logfc_header.clone());
                value_header.push(logfc_header);
            }
        }
        writer.write_record(&header).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write probe-region table header: {e}"),
            cause_chain: vec![],
        })?;
        for row in &output_rows {
            writer.write_record(row).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write probe-region row: {e}"),
                cause_chain: vec![],
            })?;
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not flush probe-region table: {e}"),
            cause_chain: vec![],
        })?;

        let mut probe_row_count = 0usize;
        let mut missing_probe_intensity_count = 0usize;
        let mut probe_intensity_source: Option<&str> = None;
        if !annotation.probes.is_empty() {
            let using_probe_intensity = probe_intensity.is_some();
            let probe_value_header = probe_intensity
                .as_ref()
                .map(|table| table.value_header.clone())
                .unwrap_or_else(|| value_header.clone());
            let row_intensity_source = if using_probe_intensity {
                "probe_level_input"
            } else {
                "parent_probeset_summary"
            };
            let probe_table_path = Path::new(output_dir).join(PROBE_REGION_PROBE_TABLE_FILE);
            let mut probe_writer =
                csv::Writer::from_path(&probe_table_path).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not create probe-coordinate table '{}': {e}",
                        probe_table_path.to_string_lossy()
                    ),
                    cause_chain: vec![],
                })?;
            let mut probe_header = vec![
                "chromosome".to_string(),
                "start".to_string(),
                "stop".to_string(),
                "strand".to_string(),
                "probe_id".to_string(),
                "x".to_string(),
                "y".to_string(),
                "parent_probeset_or_region_id".to_string(),
                "transcript_cluster_id".to_string(),
                "gene_symbol".to_string(),
                "intensity_source".to_string(),
            ];
            probe_header.extend(probe_value_header.iter().cloned());
            probe_writer
                .write_record(&probe_header)
                .map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write probe-coordinate table header: {e}"),
                    cause_chain: vec![],
                })?;
            let mut probes = annotation.probes.iter().collect::<Vec<_>>();
            probes.sort_by(|a, b| {
                a.chromosome
                    .cmp(&b.chromosome)
                    .then(a.start_1based.cmp(&b.start_1based))
                    .then(a.end_1based.cmp(&b.end_1based))
                    .then(a.probe_id.cmp(&b.probe_id))
            });
            let mut skipped_probe_rows = 0usize;
            for probe in probes {
                if !output_value_rows.contains_key(&probe.parent_feature_id) {
                    skipped_probe_rows += 1;
                    continue;
                };
                let Some(values) = probe_intensity
                    .as_ref()
                    .map(|table| table.rows.get(&probe.probe_id))
                    .unwrap_or_else(|| output_value_rows.get(&probe.parent_feature_id))
                else {
                    missing_probe_intensity_count += 1;
                    continue;
                };
                let mut row = vec![
                    probe.chromosome.clone(),
                    probe.start_1based.to_string(),
                    probe.end_1based.to_string(),
                    probe.strand.clone(),
                    probe.probe_id.clone(),
                    probe.x.clone(),
                    probe.y.clone(),
                    probe.parent_feature_id.clone(),
                    probe.transcript_cluster_id.clone(),
                    probe.gene_symbol.clone(),
                    row_intensity_source.to_string(),
                ];
                row.extend(values.iter().cloned());
                probe_writer.write_record(&row).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not write probe-coordinate row: {e}"),
                    cause_chain: vec![],
                })?;
                probe_row_count += 1;
            }
            if probe_row_count > 0 {
                probe_intensity_source = Some(row_intensity_source);
            }
            probe_writer.flush().map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not flush probe-coordinate table: {e}"),
                cause_chain: vec![],
            })?;
            if !using_probe_intensity && probe_row_count > 0 {
                warnings.push(
                    "probe_intensity_chrom_order.csv uses parent probeset-summary intensities; true PM probe intensities require probe-level intensity input"
                        .to_string(),
                );
            }
            if missing_probe_intensity_count > 0 {
                warnings.push(format!(
                    "{} probe-coordinate annotation row(s) were skipped because their probe_id was absent from the supplied probe-intensity table",
                    missing_probe_intensity_count
                ));
            }
            if skipped_probe_rows > 0 {
                warnings.push(format!(
                    "{} probe-coordinate annotation row(s) were skipped because their parent probeset was absent from the imported APT summary",
                    skipped_probe_rows
                ));
            }
        } else if probe_intensity.is_some() {
            warnings.push(
                "Probe-intensity table was supplied, but the annotation table did not provide probe_id plus probe_start/probe_stop coordinates"
                    .to_string(),
            );
        }

        let manifest_path = Path::new(output_dir).join(PROBE_REGION_MATRIX_MANIFEST_FILE);
        let provenance_path = Path::new(output_dir).join(PROBE_REGION_PROVENANCE_FILE);
        let mut artifacts = vec![PROBE_REGION_TABLE_FILE.to_string()];
        let mut targets = vec!["probeset".to_string()];
        if probe_row_count > 0 {
            artifacts.push(PROBE_REGION_PROBE_TABLE_FILE.to_string());
            targets.push("pm_probe".to_string());
        }
        let manifest = json!({
            "schema": PROBE_REGION_MATRIX_MANIFEST_SCHEMA,
            "platform": platform,
            "platform_package": null,
            "coordinate_system": coordinate_system,
            "genome_build": genome_build,
            "normalization": normalization,
            "metadata_path": metadata_summary.as_ref().map(|metadata| metadata.path.as_str()),
            "condition_column": metadata_summary
                .as_ref()
                .map(|metadata| metadata.condition_column.as_str()),
            "sample_column": metadata_summary
                .as_ref()
                .map(|metadata| metadata.sample_column.as_str()),
            "probe_intensity_path": probe_intensity.as_ref().map(|table| table.path.as_str()),
            "probe_intensity_probe_id_column": probe_intensity
                .as_ref()
                .map(|table| table.probe_id_column.as_str()),
            "probe_intensity_source": probe_intensity_source,
            "targets": targets,
            "artifacts": artifacts.clone()
        });
        let tool_version_checks = Self::probe_region_backend_tool_version_checks();
        let mut input_fingerprints = vec![
            Self::probe_region_input_fingerprint(apt_summary_path, "apt_summary"),
            Self::probe_region_input_fingerprint(annotation_path, "annotation"),
        ];
        if let Some(metadata) = metadata_summary.as_ref() {
            input_fingerprints.push(Self::probe_region_input_fingerprint(
                &metadata.path,
                "metadata",
            ));
        }
        if let Some(table) = probe_intensity.as_ref() {
            input_fingerprints.push(Self::probe_region_input_fingerprint(
                &table.path,
                "probe_intensity",
            ));
        }
        let rendered_backend_command = Self::probe_region_apt_import_replay_command(
            apt_summary_path,
            annotation_path,
            output_dir,
            metadata_summary
                .as_ref()
                .map(|metadata| metadata.path.as_str()),
            metadata_summary
                .as_ref()
                .map(|metadata| metadata.condition_column.as_str()),
            metadata_summary
                .as_ref()
                .map(|metadata| metadata.sample_column.as_str()),
            probe_intensity.as_ref().map(|table| table.path.as_str()),
            probe_id_column,
            Some(&platform),
            Some(&normalization),
            Some(&coordinate_system),
            Some(&genome_build),
        );
        let provenance = json!({
            "schema": PROBE_REGION_BACKEND_PROVENANCE_SCHEMA,
            "backend": "affymetrix_power_tools",
            "declared_backend": "affymetrix_power_tools",
            "selected_backend": "affymetrix_power_tools_imported_tables",
            "backend_execution_policy": "external_explicit_only_not_run_by_gentle",
            "backend_command_source": "GENtle import replay command; upstream APT/R summarization command is planned separately and not reconstructed from imported tables",
            "rendered_backend_command": rendered_backend_command,
            "tool_version_checks": tool_version_checks,
            "input_fingerprints": input_fingerprints,
            "apt_summary_path": apt_summary_path,
            "annotation_path": annotation_path,
            "metadata_path": metadata_summary.as_ref().map(|metadata| metadata.path.as_str()),
            "condition_column": metadata_summary
                .as_ref()
                .map(|metadata| metadata.condition_column.as_str()),
            "sample_column": metadata_summary
                .as_ref()
                .map(|metadata| metadata.sample_column.as_str()),
            "coordinate_system": coordinate_system,
            "genome_build": genome_build,
            "normalization": normalization,
            "probe_intensity_path": probe_intensity.as_ref().map(|table| table.path.as_str()),
            "probe_intensity_probe_id_column": probe_intensity
                .as_ref()
                .map(|table| table.probe_id_column.as_str()),
            "probe_intensity_source": probe_intensity_source,
            "artifacts": artifacts
        });
        std::fs::write(
            &manifest_path,
            serde_json::to_string_pretty(&manifest).unwrap_or_else(|_| "{}".to_string()),
        )
        .map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write normalized matrix manifest '{}': {e}",
                manifest_path.to_string_lossy()
            ),
            cause_chain: vec![],
        })?;
        std::fs::write(
            &provenance_path,
            serde_json::to_string_pretty(&provenance).unwrap_or_else(|_| "{}".to_string()),
        )
        .map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write probe-region provenance '{}': {e}",
                provenance_path.to_string_lossy()
            ),
            cause_chain: vec![],
        })?;

        let inspection = self.inspect_probe_region_output(output_dir)?;
        Ok(ProbeRegionAptImportReport {
            schema: PROBE_REGION_APT_IMPORT_REPORT_SCHEMA.to_string(),
            apt_summary_path: apt_summary_path.to_string(),
            annotation_path: annotation_path.to_string(),
            output_dir: output_dir.to_string(),
            platform,
            normalization,
            coordinate_system,
            genome_build,
            summary_row_count,
            annotation_row_count,
            written_row_count: output_rows.len(),
            missing_annotation_count,
            skipped_invalid_count,
            probe_row_count,
            sample_columns,
            probe_intensity_path: probe_intensity.as_ref().map(|table| table.path.clone()),
            probe_intensity_source: probe_intensity_source.map(str::to_string),
            probe_intensity_sample_columns: probe_intensity
                .as_ref()
                .map(|table| table.sample_columns.clone())
                .unwrap_or_default(),
            missing_probe_intensity_count,
            metadata_path: metadata_summary
                .as_ref()
                .map(|metadata| metadata.path.clone()),
            condition_column: metadata_summary
                .as_ref()
                .map(|metadata| metadata.condition_column.clone()),
            sample_column: metadata_summary
                .as_ref()
                .map(|metadata| metadata.sample_column.clone()),
            condition_columns: metadata_summary
                .as_ref()
                .map(|metadata| {
                    metadata
                        .conditions
                        .iter()
                        .map(|condition| condition.label.clone())
                        .collect()
                })
                .unwrap_or_default(),
            logfc_columns: metadata_summary
                .as_ref()
                .map(|metadata| {
                    metadata
                        .contrasts
                        .iter()
                        .map(|contrast| contrast.label.clone())
                        .collect()
                })
                .unwrap_or_default(),
            warnings,
            inspection,
        })
    }

    /// Inspect outputs written by the explicit R/oligo probe-region helper.

    pub(super) fn probe_region_apt_annotation_rows(
        annotation_path: &str,
    ) -> Result<ProbeRegionAptAnnotationTable, EngineError> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(Self::probe_region_metadata_delimiter(annotation_path))
            .trim(csv::Trim::All)
            .flexible(true)
            .from_path(annotation_path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not open APT annotation table '{annotation_path}': {e}"),
                cause_chain: vec![],
            })?;
        let headers = reader
            .headers()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read APT annotation table header '{annotation_path}': {e}"
                ),
                cause_chain: vec![],
            })?
            .iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        let feature_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &[
                "probeset_or_region_id",
                "probeset_id",
                "probeset",
                "feature_id",
                "id",
            ],
        )
        .ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: "APT annotation table is missing probeset/feature id column".to_string(),
            cause_chain: vec![],
        })?;
        let chromosome_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["chromosome", "chrom"])
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "APT annotation table is missing chromosome column".to_string(),
                    cause_chain: vec![],
                })?;
        let start_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["start", "start_1based", "genomic_start_1based"],
        )
        .ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: "APT annotation table is missing start column".to_string(),
            cause_chain: vec![],
        })?;
        let stop_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["stop", "end", "end_1based", "genomic_end_1based"],
        )
        .ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: "APT annotation table is missing stop/end column".to_string(),
            cause_chain: vec![],
        })?;
        let strand_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["strand", "orientation"]);
        let transcript_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["transcript_cluster_id", "transcript_cluster"],
        );
        let number_of_probes_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["number_of_probes", "num_probes", "probe_count"],
        );
        let gene_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["gene_symbol", "gene"]);
        let probe_id_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["probe_id", "probeid", "probe", "pm_probe_id"],
        );
        let probe_chromosome_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &[
                "probe_chromosome",
                "probe_chrom",
                "pm_chromosome",
                "pm_chrom",
            ],
        );
        let probe_start_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &[
                "probe_start",
                "probe_start_1based",
                "probe_genomic_start_1based",
                "pm_start",
                "pm_start_1based",
            ],
        );
        let probe_stop_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &[
                "probe_stop",
                "probe_end",
                "probe_stop_1based",
                "probe_end_1based",
                "probe_genomic_end_1based",
                "pm_stop",
                "pm_end",
                "pm_stop_1based",
                "pm_end_1based",
            ],
        );
        let probe_x_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["x", "probe_x", "x_coord", "x_coordinate", "feature_x"],
        );
        let probe_y_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["y", "probe_y", "y_coord", "y_coordinate", "feature_y"],
        );

        let mut rows = HashMap::new();
        let mut probes = Vec::new();
        for (line_offset, record) in reader.records().enumerate() {
            let record = record.map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not read APT annotation row: {e}"),
                cause_chain: vec![],
            })?;
            let line_no = line_offset + 2;
            let feature_id = record.get(feature_idx).unwrap_or("").trim().to_string();
            if feature_id.is_empty() {
                continue;
            }
            let start_1based = record
                .get(start_idx)
                .unwrap_or("")
                .trim()
                .parse::<usize>()
                .map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Invalid APT annotation start coordinate on line {line_no}: {e}"
                    ),
                    cause_chain: vec![],
                })?;
            let end_1based = record
                .get(stop_idx)
                .unwrap_or("")
                .trim()
                .parse::<usize>()
                .map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Invalid APT annotation stop coordinate on line {line_no}: {e}"
                    ),
                    cause_chain: vec![],
                })?;
            if start_1based == 0 || end_1based < start_1based {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Invalid APT annotation interval on line {line_no}: start={start_1based}, end={end_1based}"
                    ),
                    cause_chain: vec![],
                });
            }
            let chromosome = record.get(chromosome_idx).unwrap_or("").trim().to_string();
            let strand = strand_idx
                .and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| matches!(*value, "+" | "-"))
                .unwrap_or("+")
                .to_string();
            let transcript_cluster_id = transcript_idx
                .and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or("")
                .to_string();
            let number_of_probes = number_of_probes_idx
                .and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or("")
                .to_string();
            let gene_symbol = gene_idx
                .and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or("")
                .to_string();
            rows.entry(feature_id.clone())
                .or_insert_with(|| ProbeRegionAptAnnotationRow {
                    chromosome: chromosome.clone(),
                    start_1based,
                    end_1based,
                    strand: strand.clone(),
                    feature_id: feature_id.clone(),
                    transcript_cluster_id: transcript_cluster_id.clone(),
                    number_of_probes: number_of_probes.clone(),
                    gene_symbol: gene_symbol.clone(),
                });
            if let (Some(probe_id_idx), Some(probe_start_idx), Some(probe_stop_idx)) =
                (probe_id_idx, probe_start_idx, probe_stop_idx)
            {
                let probe_id = record.get(probe_id_idx).unwrap_or("").trim();
                if !probe_id.is_empty() {
                    let probe_start_1based = record
                        .get(probe_start_idx)
                        .unwrap_or("")
                        .trim()
                        .parse::<usize>()
                        .map_err(|e| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Invalid APT probe start coordinate on line {line_no}: {e}"
                            ),
                            cause_chain: vec![],
                        })?;
                    let probe_end_1based = record
                        .get(probe_stop_idx)
                        .unwrap_or("")
                        .trim()
                        .parse::<usize>()
                        .map_err(|e| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Invalid APT probe stop coordinate on line {line_no}: {e}"
                            ),
                            cause_chain: vec![],
                        })?;
                    if probe_start_1based == 0 || probe_end_1based < probe_start_1based {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Invalid APT probe interval on line {line_no}: start={probe_start_1based}, end={probe_end_1based}"
                            ),
                            cause_chain: vec![],
                        });
                    }
                    probes.push(ProbeRegionAptProbeAnnotationRow {
                        chromosome: probe_chromosome_idx
                            .and_then(|idx| record.get(idx))
                            .map(str::trim)
                            .filter(|value| !value.is_empty())
                            .unwrap_or(&chromosome)
                            .to_string(),
                        start_1based: probe_start_1based,
                        end_1based: probe_end_1based,
                        strand: strand.clone(),
                        parent_feature_id: feature_id.clone(),
                        transcript_cluster_id: transcript_cluster_id.clone(),
                        gene_symbol: gene_symbol.clone(),
                        probe_id: probe_id.to_string(),
                        x: probe_x_idx
                            .and_then(|idx| record.get(idx))
                            .map(str::trim)
                            .unwrap_or("")
                            .to_string(),
                        y: probe_y_idx
                            .and_then(|idx| record.get(idx))
                            .map(str::trim)
                            .unwrap_or("")
                            .to_string(),
                    });
                }
            }
        }
        if rows.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("APT annotation table '{annotation_path}' had no usable rows"),
                cause_chain: vec![],
            });
        }
        Ok(ProbeRegionAptAnnotationTable {
            regions: rows,
            probes,
        })
    }

    pub(super) fn probe_region_apt_probe_intensity_rows(
        probe_intensity_path: &str,
        requested_probe_id_column: Option<&str>,
        metadata_path: Option<&str>,
        condition_column: Option<&str>,
        sample_column: Option<&str>,
    ) -> Result<ProbeRegionAptProbeIntensityTable, EngineError> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(Self::probe_region_metadata_delimiter(probe_intensity_path))
            .trim(csv::Trim::All)
            .flexible(true)
            .from_path(probe_intensity_path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not open APT probe-intensity table '{probe_intensity_path}': {e}"
                ),
                cause_chain: vec![],
            })?;
        let headers = reader
            .headers()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read APT probe-intensity table header '{probe_intensity_path}': {e}"
                ),
                cause_chain: vec![],
            })?
            .iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        let probe_id_idx = Self::probe_region_metadata_column_index(
            &headers,
            requested_probe_id_column,
            &["probe_id", "probe id", "probeid", "probe", "id"],
        )
        .ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: "APT probe-intensity table is missing a probe id column".to_string(),
            cause_chain: vec![],
        })?;
        let sample_indices = headers
            .iter()
            .enumerate()
            .filter(|(idx, header)| *idx != probe_id_idx && !header.trim().is_empty())
            .map(|(idx, _)| idx)
            .collect::<Vec<_>>();
        let sample_columns = sample_indices
            .iter()
            .filter_map(|idx| headers.get(*idx))
            .map(|header| header.trim().to_string())
            .collect::<Vec<_>>();
        if sample_columns.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "APT probe-intensity table '{probe_intensity_path}' has no sample intensity columns"
                ),
                cause_chain: vec![],
            });
        }
        let metadata_summary = if let Some(path) = metadata_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            Some(Self::probe_region_apt_metadata_summary(
                path,
                condition_column,
                sample_column,
                &sample_columns,
            )?)
        } else {
            None
        };
        let mut warnings = metadata_summary
            .as_ref()
            .map(|summary| summary.warnings.clone())
            .unwrap_or_default();
        let mut value_header = sample_columns.clone();
        if let Some(metadata) = metadata_summary.as_ref() {
            for condition in &metadata.conditions {
                value_header.push(format!("mean_log2_{}", condition.column_label));
                value_header.push(format!("sd_log2_{}", condition.column_label));
            }
            for contrast in &metadata.contrasts {
                value_header.push(format!("log2FC_{}", contrast.label));
            }
        }

        let mut rows = BTreeMap::new();
        let mut skipped_invalid_count = 0usize;
        let mut duplicate_probe_count = 0usize;
        for record in reader.records() {
            let record = record.map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not read APT probe-intensity row: {e}"),
                cause_chain: vec![],
            })?;
            let probe_id = record.get(probe_id_idx).unwrap_or("").trim();
            if probe_id.is_empty() {
                skipped_invalid_count += 1;
                continue;
            }
            let sample_values = sample_indices
                .iter()
                .map(|idx| record.get(*idx).unwrap_or("").trim().to_string())
                .collect::<Vec<_>>();
            let mut value_row = sample_values.clone();
            if let Some(metadata) = metadata_summary.as_ref() {
                let stats = metadata
                    .conditions
                    .iter()
                    .map(|condition| {
                        Self::probe_region_apt_condition_stats(
                            &sample_values,
                            &condition.sample_indices,
                        )
                    })
                    .collect::<Vec<_>>();
                for stat in &stats {
                    value_row.push(
                        stat.map(|(mean, _)| Self::probe_region_format_number(mean))
                            .unwrap_or_default(),
                    );
                    value_row.push(
                        stat.map(|(_, sd)| Self::probe_region_format_number(sd))
                            .unwrap_or_default(),
                    );
                }
                for contrast in &metadata.contrasts {
                    let numerator = stats
                        .get(contrast.numerator_index)
                        .and_then(|value| value.map(|(mean, _)| mean));
                    let denominator = stats
                        .get(contrast.denominator_index)
                        .and_then(|value| value.map(|(mean, _)| mean));
                    value_row.push(
                        numerator
                            .zip(denominator)
                            .map(|(num, den)| Self::probe_region_format_number(num - den))
                            .unwrap_or_default(),
                    );
                }
            }
            if rows.insert(probe_id.to_string(), value_row).is_some() {
                duplicate_probe_count += 1;
            }
        }
        if rows.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "APT probe-intensity table '{probe_intensity_path}' had no usable probe rows"
                ),
                cause_chain: vec![],
            });
        }
        if skipped_invalid_count > 0 {
            warnings.push(format!(
                "{} APT probe-intensity row(s) were skipped because probe id was empty",
                skipped_invalid_count
            ));
        }
        if duplicate_probe_count > 0 {
            warnings.push(format!(
                "{} duplicate APT probe-intensity row(s) reused the last value for that probe id",
                duplicate_probe_count
            ));
        }

        Ok(ProbeRegionAptProbeIntensityTable {
            path: probe_intensity_path.to_string(),
            probe_id_column: headers.get(probe_id_idx).cloned().unwrap_or_default(),
            sample_columns,
            value_header,
            rows,
            warnings,
        })
    }

    pub(super) fn probe_region_apt_metadata_summary(
        metadata_path: &str,
        requested_condition_column: Option<&str>,
        requested_sample_column: Option<&str>,
        sample_columns: &[String],
    ) -> Result<ProbeRegionAptMetadataSummary, EngineError> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(Self::probe_region_metadata_delimiter(metadata_path))
            .trim(csv::Trim::All)
            .flexible(true)
            .from_path(metadata_path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not open APT sample metadata '{metadata_path}': {e}"),
                cause_chain: vec![],
            })?;
        let headers = reader
            .headers()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read APT sample metadata header '{metadata_path}': {e}"
                ),
                cause_chain: vec![],
            })?
            .iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        let sample_idx = Self::probe_region_metadata_column_index(
            &headers,
            requested_sample_column,
            &[
                "file",
                "cel",
                "cel_file",
                "cel file",
                "array data file",
                "sample",
                "sample_name",
                "source name",
            ],
        )
        .ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: "APT sample metadata is missing a sample/CEL column".to_string(),
            cause_chain: vec![],
        })?;
        let condition_idx = Self::probe_region_metadata_column_index(
            &headers,
            requested_condition_column,
            &[
                "condition",
                "group",
                "treatment",
                "sample group",
                "factor value condition",
                "factor value treatment",
                "characteristics condition",
            ],
        )
        .ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: "APT sample metadata is missing a condition/group column".to_string(),
            cause_chain: vec![],
        })?;
        let (sample_lookup, mut warnings) = Self::probe_region_apt_sample_lookup(sample_columns);
        let mut condition_samples: BTreeMap<String, BTreeSet<usize>> = BTreeMap::new();
        let mut sample_conditions: BTreeMap<usize, String> = BTreeMap::new();
        let mut unmatched_samples = Vec::new();
        let mut row_count = 0usize;
        for record in reader.records() {
            let record = record.map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not read APT sample metadata row: {e}"),
                cause_chain: vec![],
            })?;
            row_count += 1;
            let sample = record.get(sample_idx).unwrap_or("").trim();
            let condition = record.get(condition_idx).unwrap_or("").trim();
            if sample.is_empty() || condition.is_empty() {
                continue;
            }
            let matched_idx = Self::probe_region_apt_sample_match_keys(sample)
                .iter()
                .find_map(|key| sample_lookup.get(key).copied());
            let Some(sample_column_idx) = matched_idx else {
                if unmatched_samples.len() < 8 {
                    unmatched_samples.push(sample.to_string());
                }
                continue;
            };
            if let Some(previous) = sample_conditions.get(&sample_column_idx) {
                if previous != condition {
                    warnings.push(format!(
                        "Sample '{}' appears with multiple conditions ('{}' and '{}'); keeping '{}'",
                        sample, previous, condition, previous
                    ));
                    continue;
                }
            } else {
                sample_conditions.insert(sample_column_idx, condition.to_string());
            }
            condition_samples
                .entry(condition.to_string())
                .or_default()
                .insert(sample_column_idx);
        }
        if row_count == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("APT sample metadata '{metadata_path}' had no rows"),
                cause_chain: vec![],
            });
        }
        if !unmatched_samples.is_empty() {
            warnings.push(format!(
                "APT sample metadata rows did not match summary columns for: {}{}",
                unmatched_samples.join(", "),
                if unmatched_samples.len() >= 8 {
                    " ..."
                } else {
                    ""
                }
            ));
        }
        if condition_samples.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "APT sample metadata '{metadata_path}' had no rows matching summary sample columns"
                ),
                cause_chain: vec![],
            });
        }

        let mut used_labels: BTreeMap<String, usize> = BTreeMap::new();
        let conditions = condition_samples
            .into_iter()
            .map(|(label, sample_set)| {
                let base_label = Self::probe_region_apt_output_label(&label);
                let count = used_labels.entry(base_label.clone()).or_insert(0);
                *count += 1;
                let column_label = if *count == 1 {
                    base_label
                } else {
                    format!("{base_label}_{count}")
                };
                ProbeRegionAptConditionGroup {
                    label,
                    column_label,
                    sample_indices: sample_set.into_iter().collect(),
                }
            })
            .collect::<Vec<_>>();
        let contrasts = if conditions.len() < 2 {
            Vec::new()
        } else {
            let denominator = &conditions[0];
            conditions
                .iter()
                .enumerate()
                .skip(1)
                .map(|(idx, numerator)| ProbeRegionAptContrast {
                    label: format!("{}-{}", numerator.column_label, denominator.column_label),
                    numerator_index: idx,
                    denominator_index: 0,
                })
                .collect()
        };

        Ok(ProbeRegionAptMetadataSummary {
            path: metadata_path.to_string(),
            sample_column: headers.get(sample_idx).cloned().unwrap_or_default(),
            condition_column: headers.get(condition_idx).cloned().unwrap_or_default(),
            conditions,
            contrasts,
            warnings,
        })
    }

    pub(super) fn probe_region_apt_sample_lookup(
        sample_columns: &[String],
    ) -> (BTreeMap<String, usize>, Vec<String>) {
        let mut lookup = BTreeMap::new();
        let mut warnings = Vec::new();
        for (idx, sample) in sample_columns.iter().enumerate() {
            for key in Self::probe_region_apt_sample_match_keys(sample) {
                if let Some(previous) = lookup.insert(key.clone(), idx)
                    && previous != idx
                {
                    warnings.push(format!(
                        "APT summary sample columns '{}' and '{}' share matching key '{}'; metadata matching will use the latter",
                        sample_columns
                            .get(previous)
                            .map(String::as_str)
                            .unwrap_or("<unknown>"),
                        sample,
                        key
                    ));
                }
            }
        }
        (lookup, warnings)
    }

    pub(super) fn probe_region_apt_sample_match_keys(value: &str) -> Vec<String> {
        let trimmed = value.trim().trim_matches('"').trim_matches('\'');
        if trimmed.is_empty() {
            return Vec::new();
        }
        let basename = trimmed.rsplit(['/', '\\']).next().unwrap_or(trimmed).trim();
        let mut keys = Vec::new();
        for candidate in [trimmed, basename] {
            let lower = candidate.to_ascii_lowercase();
            if !lower.is_empty() && !keys.contains(&lower) {
                keys.push(lower.clone());
            }
            if let Some(stripped) = lower.strip_suffix(".cel")
                && !stripped.is_empty()
                && !keys.iter().any(|key| key == stripped)
            {
                keys.push(stripped.to_string());
            }
        }
        keys
    }

    pub(super) fn probe_region_apt_output_label(value: &str) -> String {
        let mut out = String::new();
        let mut last_was_separator = false;
        for ch in value.trim().chars() {
            if ch.is_ascii_alphanumeric() || ch == '-' || ch == '.' {
                out.push(ch);
                last_was_separator = false;
            } else if !last_was_separator {
                out.push('_');
                last_was_separator = true;
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "condition".to_string()
        } else {
            trimmed.to_string()
        }
    }

    pub(super) fn probe_region_apt_condition_stats(
        sample_values: &[String],
        sample_indices: &[usize],
    ) -> Option<(f64, f64)> {
        let values = sample_indices
            .iter()
            .filter_map(|idx| sample_values.get(*idx))
            .filter_map(|value| value.trim().parse::<f64>().ok())
            .filter(|value| value.is_finite())
            .collect::<Vec<_>>();
        if values.is_empty() {
            return None;
        }
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let sd = if values.len() > 1 {
            let variance = values
                .iter()
                .map(|value| {
                    let delta = value - mean;
                    delta * delta
                })
                .sum::<f64>()
                / (values.len() - 1) as f64;
            variance.sqrt()
        } else {
            0.0
        };
        Some((mean, sd))
    }
}
