use super::*;

impl GentleEngine {
    pub fn inspect_probe_region_output(
        &self,
        output_dir: &str,
    ) -> Result<ProbeRegionOutputInspection, EngineError> {
        let output_dir = output_dir.trim();
        if output_dir.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Probe-region output directory must not be empty".to_string(),
                cause_chain: vec![],
            });
        }

        let output_path = Path::new(output_dir);
        let region_table_path = output_path.join(PROBE_REGION_TABLE_FILE);
        let probe_table_path = output_path.join(PROBE_REGION_PROBE_TABLE_FILE);
        let sample_table_path = output_path.join(PROBE_REGION_SAMPLE_TABLE_FILE);
        let manifest_path = output_path.join(PROBE_REGION_MATRIX_MANIFEST_FILE);
        let provenance_path = output_path.join(PROBE_REGION_PROVENANCE_FILE);

        let region_table =
            Self::probe_region_file_status(&region_table_path.to_string_lossy(), "region_table");
        let probe_table = Self::probe_region_optional_file_status(&probe_table_path, "probe_table");
        let sample_table =
            Self::probe_region_optional_file_status(&sample_table_path, "sample_table");
        let normalized_matrix_manifest =
            Self::probe_region_optional_file_status(&manifest_path, "normalized_matrix_manifest");
        let provenance = Self::probe_region_optional_file_status(&provenance_path, "provenance");

        let mut warnings = Vec::new();
        let mut errors = Vec::new();
        if !output_path.exists() {
            errors.push(format!(
                "Probe-region output directory '{}' does not exist",
                output_dir
            ));
        } else if !output_path.is_dir() {
            errors.push(format!(
                "Probe-region output path '{}' is not a directory",
                output_dir
            ));
        }
        if !region_table.exists {
            errors.push(format!(
                "Missing required probe-region table '{}'",
                region_table.path
            ));
        } else if !region_table.is_file {
            errors.push(format!(
                "Probe-region table path '{}' is not a file",
                region_table.path
            ));
        }
        if sample_table.is_none() {
            warnings.push(format!(
                "Optional helper output '{}' is missing; sample columns will be inferred from the region table",
                sample_table_path.to_string_lossy()
            ));
        }
        if normalized_matrix_manifest.is_none() {
            warnings.push(format!(
                "Optional helper output '{}' is missing; normalization provenance is incomplete",
                manifest_path.to_string_lossy()
            ));
        }
        if provenance.is_none() {
            warnings.push(format!(
                "Optional helper output '{}' is missing; backend provenance is incomplete",
                provenance_path.to_string_lossy()
            ));
        }

        let sample_columns_from_table = if sample_table
            .as_ref()
            .is_some_and(|status| status.exists && status.is_file)
        {
            match Self::probe_region_sample_ids_from_sample_table(&sample_table_path) {
                Ok(columns) => columns,
                Err(err) => {
                    warnings.push(format!("Could not parse sample table: {err}"));
                    Vec::new()
                }
            }
        } else {
            Vec::new()
        };

        let mut backend = None;
        let mut platform = None;
        let mut platform_package = None;
        let mut normalization = None;
        let mut coordinate_system = None;
        let mut genome_build = None;
        let mut coordinate_projections = Vec::new();
        let mut target_levels = Vec::new();
        let mut artifact_paths = Vec::new();
        if normalized_matrix_manifest
            .as_ref()
            .is_some_and(|status| status.exists && status.is_file)
        {
            match Self::probe_region_json_file(&manifest_path) {
                Ok(value) => {
                    Self::probe_region_warn_unexpected_json_schema(
                        &value,
                        PROBE_REGION_MATRIX_MANIFEST_SCHEMA,
                        &mut warnings,
                        "normalized matrix manifest",
                    );
                    platform = Self::probe_region_json_string(&value, "platform");
                    platform_package = Self::probe_region_json_string(&value, "platform_package");
                    normalization = Self::probe_region_json_string(&value, "normalization");
                    coordinate_system = Self::probe_region_json_string(&value, "coordinate_system");
                    genome_build = Self::probe_region_json_string(&value, "genome_build")
                        .or_else(|| Self::probe_region_json_string(&value, "reference_genome_id"));
                    coordinate_projections = Self::probe_region_json_coordinate_projections(&value);
                    target_levels = Self::probe_region_json_string_array(&value, "targets");
                    artifact_paths = Self::probe_region_json_string_array(&value, "artifacts");
                }
                Err(err) => {
                    errors.push(format!("Could not parse normalized matrix manifest: {err}"))
                }
            }
        }
        if provenance
            .as_ref()
            .is_some_and(|status| status.exists && status.is_file)
        {
            match Self::probe_region_json_file(&provenance_path) {
                Ok(value) => {
                    Self::probe_region_warn_unexpected_json_schema(
                        &value,
                        PROBE_REGION_BACKEND_PROVENANCE_SCHEMA,
                        &mut warnings,
                        "backend provenance",
                    );
                    backend = Self::probe_region_json_string(&value, "backend");
                    if platform_package.is_none() {
                        platform_package =
                            Self::probe_region_json_string(&value, "platform_package");
                    }
                    if normalization.is_none() {
                        normalization = Self::probe_region_json_string(&value, "normalization");
                    }
                    if coordinate_system.is_none() {
                        coordinate_system =
                            Self::probe_region_json_string(&value, "coordinate_system");
                    }
                    if genome_build.is_none() {
                        genome_build = Self::probe_region_json_string(&value, "genome_build")
                            .or_else(|| {
                                Self::probe_region_json_string(&value, "reference_genome_id")
                            });
                    }
                    Self::probe_region_extend_unique_coordinate_projections(
                        &mut coordinate_projections,
                        Self::probe_region_json_coordinate_projections(&value),
                    );
                    Self::probe_region_extend_unique(
                        &mut artifact_paths,
                        Self::probe_region_json_string_array(&value, "artifacts"),
                    );
                }
                Err(err) => errors.push(format!("Could not parse backend provenance: {err}")),
            }
        }

        let table_summary = if region_table.exists && region_table.is_file {
            match Self::inspect_probe_region_table(&region_table_path, &sample_columns_from_table) {
                Ok(summary) => summary,
                Err(err) => {
                    errors.push(format!("Could not parse probe-region table: {err}"));
                    ProbeRegionTableSummary::default()
                }
            }
        } else {
            ProbeRegionTableSummary::default()
        };
        warnings.extend(table_summary.warnings.iter().cloned());
        errors.extend(
            table_summary
                .required_columns_missing
                .iter()
                .map(|column| format!("Probe-region table is missing required column '{column}'")),
        );
        let probe_table_summary = if probe_table
            .as_ref()
            .is_some_and(|status| status.exists && status.is_file)
        {
            match Self::inspect_probe_region_probe_table(&probe_table_path) {
                Ok(summary) => summary,
                Err(err) => {
                    errors.push(format!("Could not parse probe-coordinate table: {err}"));
                    ProbeRegionProbeTableSummary::default()
                }
            }
        } else {
            ProbeRegionProbeTableSummary::default()
        };
        errors.extend(
            probe_table_summary
                .required_columns_missing
                .iter()
                .map(|column| {
                    format!("Probe-coordinate table is missing required column '{column}'")
                }),
        );

        let mut projection_blockers = Vec::new();
        if coordinate_system.is_none() {
            projection_blockers.push(
                "Missing coordinate_system; genome-anchored projection must be refused until the helper output declares its coordinate basis"
                    .to_string(),
            );
        }
        if genome_build.is_none() {
            projection_blockers.push(
                "Missing genome_build; genome-anchored projection must be refused until the helper output declares its reference build"
                    .to_string(),
            );
        }
        if !table_summary.required_columns_missing.is_empty() {
            projection_blockers.push(format!(
                "Missing required region table column(s): {}",
                table_summary.required_columns_missing.join(", ")
            ));
        }
        let projection_ready = projection_blockers.is_empty();

        Ok(ProbeRegionOutputInspection {
            schema: PROBE_REGION_OUTPUT_INSPECTION_SCHEMA.to_string(),
            output_dir: output_dir.to_string(),
            usable: errors.is_empty(),
            region_table,
            probe_table,
            sample_table,
            normalized_matrix_manifest,
            provenance,
            backend,
            platform,
            platform_package,
            normalization,
            coordinate_system,
            genome_build,
            coordinate_projections,
            projection_ready,
            projection_blockers,
            target_levels,
            artifact_paths,
            row_count: table_summary.row_count,
            column_count: table_summary.column_count,
            probe_row_count: probe_table_summary.row_count,
            probe_parent_feature_count: probe_table_summary.parent_feature_count,
            feature_count: table_summary.feature_count,
            transcript_cluster_count: table_summary.transcript_cluster_count,
            chromosome_count: table_summary.chromosome_count,
            chromosomes: table_summary.chromosomes,
            gene_symbols: table_summary.gene_symbols,
            sample_columns: table_summary.sample_columns,
            condition_summary_columns: table_summary.condition_summary_columns,
            logfc_columns: table_summary.logfc_columns,
            preview_rows: table_summary.preview_rows,
            required_columns_missing: table_summary.required_columns_missing,
            warnings,
            errors,
        })
    }

    /// Render inspected probe-region output as a deterministic SVG plot.

    pub(super) fn inspect_probe_region_probe_table(
        path: &Path,
    ) -> Result<ProbeRegionProbeTableSummary, String> {
        let mut reader = csv::ReaderBuilder::new()
            .trim(csv::Trim::All)
            .from_path(path)
            .map_err(|e| format!("could not open '{}': {e}", path.to_string_lossy()))?;
        let headers = reader
            .headers()
            .map_err(|e| format!("could not read header: {e}"))?
            .iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        let chromosome_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["chromosome", "chrom"]);
        let start_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["start", "start_1based"]);
        let stop_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["stop", "end", "end_1based"],
        );
        let probe_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["probe_id", "probe"]);
        let parent_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &[
                "parent_probeset_or_region_id",
                "parent_probeset_id",
                "probeset_id",
            ],
        );
        let mut required_columns_missing = Vec::new();
        if chromosome_idx.is_none() {
            required_columns_missing.push("chromosome".to_string());
        }
        if start_idx.is_none() {
            required_columns_missing.push("start".to_string());
        }
        if stop_idx.is_none() {
            required_columns_missing.push("stop".to_string());
        }
        if probe_idx.is_none() {
            required_columns_missing.push("probe_id".to_string());
        }
        if parent_idx.is_none() {
            required_columns_missing.push("parent_probeset_or_region_id".to_string());
        }
        if !required_columns_missing.is_empty() {
            return Ok(ProbeRegionProbeTableSummary {
                required_columns_missing,
                ..Default::default()
            });
        }

        let chromosome_idx = chromosome_idx.unwrap_or_default();
        let start_idx = start_idx.unwrap_or_default();
        let stop_idx = stop_idx.unwrap_or_default();
        let parent_idx = parent_idx.unwrap_or_default();
        let mut row_count = 0usize;
        let mut parents = BTreeSet::new();
        for record in reader.records() {
            let record = record.map_err(|e| format!("could not read probe row: {e}"))?;
            row_count += 1;
            let chromosome = record.get(chromosome_idx).unwrap_or("").trim();
            let start = record
                .get(start_idx)
                .unwrap_or("")
                .trim()
                .parse::<usize>()
                .map_err(|e| format!("invalid probe start on row {}: {e}", row_count + 1))?;
            let stop = record
                .get(stop_idx)
                .unwrap_or("")
                .trim()
                .parse::<usize>()
                .map_err(|e| format!("invalid probe stop on row {}: {e}", row_count + 1))?;
            if chromosome.is_empty() || start == 0 || stop < start {
                return Err(format!(
                    "invalid probe interval on row {}: {chromosome}:{start}-{stop}",
                    row_count + 1
                ));
            }
            if let Some(parent) = record
                .get(parent_idx)
                .map(str::trim)
                .filter(|v| !v.is_empty())
            {
                parents.insert(parent.to_string());
            }
        }
        Ok(ProbeRegionProbeTableSummary {
            row_count,
            parent_feature_count: parents.len(),
            required_columns_missing,
        })
    }

    pub(super) fn inspect_probe_region_table(
        path: &Path,
        sample_ids_from_table: &[String],
    ) -> Result<ProbeRegionTableSummary, String> {
        let mut reader = csv::ReaderBuilder::new()
            .trim(csv::Trim::All)
            .from_path(path)
            .map_err(|e| format!("could not open '{}': {e}", path.to_string_lossy()))?;
        let headers = reader
            .headers()
            .map_err(|e| format!("could not read header: {e}"))?
            .iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        let chromosome_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["chromosome", "chrom"]);
        let start_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["start", "start_1based"]);
        let stop_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["stop", "end", "end_1based"],
        );
        let feature_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &[
                "probeset_or_region_id",
                "feature_id",
                "probeset_id",
                "psr_id",
            ],
        );
        let transcript_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["transcript_cluster_id", "transcript_cluster"],
        );
        let gene_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["gene_symbol", "gene"]);

        let mut required_columns_missing = Vec::new();
        if chromosome_idx.is_none() {
            required_columns_missing.push("chromosome".to_string());
        }
        if start_idx.is_none() {
            required_columns_missing.push("start".to_string());
        }
        if stop_idx.is_none() {
            required_columns_missing.push("stop".to_string());
        }
        if feature_idx.is_none() {
            required_columns_missing.push("probeset_or_region_id".to_string());
        }

        let condition_summary_columns = headers
            .iter()
            .filter(|header| header.starts_with("mean_log2_") || header.starts_with("sd_log2_"))
            .cloned()
            .collect::<Vec<_>>();
        let logfc_columns = headers
            .iter()
            .filter(|header| header.starts_with("log2FC_"))
            .cloned()
            .collect::<Vec<_>>();
        let sample_columns =
            Self::probe_region_region_sample_columns(&headers, sample_ids_from_table);
        let mut warnings = Vec::new();
        if !sample_ids_from_table.is_empty() {
            let header_keys = headers
                .iter()
                .map(|header| header.to_ascii_lowercase())
                .collect::<BTreeSet<_>>();
            let missing_samples = sample_ids_from_table
                .iter()
                .filter(|sample_id| !header_keys.contains(&sample_id.to_ascii_lowercase()))
                .cloned()
                .collect::<Vec<_>>();
            if !missing_samples.is_empty() {
                warnings.push(format!(
                    "sample_table.tsv sample(s) absent from region table columns: {}",
                    missing_samples.join(", ")
                ));
            }
        }

        let mut chromosomes = BTreeSet::new();
        let mut genes = BTreeSet::new();
        let mut features = BTreeSet::new();
        let mut transcripts = BTreeSet::new();
        let mut row_count = 0usize;
        let mut preview_rows = Vec::new();
        for record in reader.records() {
            let record = record.map_err(|e| format!("could not read region row: {e}"))?;
            row_count += 1;
            if preview_rows.len() < 12 {
                let value_at = |idx: Option<usize>| -> String {
                    idx.and_then(|idx| record.get(idx))
                        .unwrap_or("")
                        .trim()
                        .to_string()
                };
                let parse_position =
                    |idx: Option<usize>| -> Option<usize> { value_at(idx).parse().ok() };
                preview_rows.push(ProbeRegionOutputPreviewRow {
                    chromosome: value_at(chromosome_idx),
                    start_1based: parse_position(start_idx),
                    stop_1based: parse_position(stop_idx),
                    probeset_or_region_id: value_at(feature_idx),
                    transcript_cluster_id: value_at(transcript_idx),
                    gene_symbol: value_at(gene_idx),
                });
            }
            if let Some(idx) = chromosome_idx {
                let value = record.get(idx).unwrap_or("").trim();
                if !value.is_empty() {
                    chromosomes.insert(value.to_string());
                }
            }
            if let Some(idx) = gene_idx {
                let value = record.get(idx).unwrap_or("").trim();
                if !value.is_empty() {
                    genes.insert(value.to_string());
                }
            }
            if let Some(idx) = feature_idx {
                let value = record.get(idx).unwrap_or("").trim();
                if !value.is_empty() {
                    features.insert(value.to_string());
                }
            }
            if let Some(idx) = transcript_idx {
                let value = record.get(idx).unwrap_or("").trim();
                if !value.is_empty() {
                    transcripts.insert(value.to_string());
                }
            }
        }

        Ok(ProbeRegionTableSummary {
            row_count,
            column_count: headers.len(),
            feature_count: features.len(),
            transcript_cluster_count: transcripts.len(),
            chromosome_count: chromosomes.len(),
            chromosomes: Self::probe_region_preview_values(chromosomes),
            gene_symbols: Self::probe_region_preview_values(genes),
            sample_columns,
            condition_summary_columns,
            logfc_columns,
            preview_rows,
            required_columns_missing,
            warnings,
        })
    }
}
