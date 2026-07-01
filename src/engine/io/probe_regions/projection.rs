use super::*;

impl GentleEngine {
    pub(in crate::engine) fn project_probe_region_output(
        dna: &mut DNAsequence,
        anchor: &GenomeSequenceAnchor,
        seq_id: &str,
        output_dir: &str,
        contrasts: &[String],
        level: Option<&str>,
        min_abs_logfc: Option<f64>,
        max_features: Option<usize>,
        clear_existing: bool,
    ) -> Result<MicroarrayProjectionReport, EngineError> {
        let output_dir = output_dir.trim();
        if output_dir.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ProjectProbeRegionOutput requires a non-empty output_dir".to_string(),
                cause_chain: vec![],
            });
        }
        let inspection = GentleEngine::default().inspect_probe_region_output(output_dir)?;
        let coordinate_system =
            inspection
                .coordinate_system
                .clone()
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Probe-region output '{}' does not declare coordinate_system",
                        output_dir
                    ),
                    cause_chain: vec![],
                })?;
        let genome_build = inspection.genome_build.clone().ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Probe-region output '{}' does not declare genome_build",
                output_dir
            ),
            cause_chain: vec![],
        })?;
        let anchor_matches =
            Self::probe_region_output_supports_anchor(&coordinate_system, &genome_build, anchor);
        let projection_spec = if anchor_matches {
            None
        } else {
            Self::resolve_probe_region_coordinate_projection_spec(&inspection, anchor)
        };
        let projection_path = projection_spec
            .map(|spec| Self::resolve_probe_region_output_path(output_dir, &spec.path));
        let projection_blocks =
            if let (Some(spec), Some(path)) = (projection_spec, &projection_path) {
                Some(Self::load_genome_coordinate_projection_blocks(
                    path,
                    &spec.method,
                    &spec.source_genome_id,
                    &spec.target_genome_id,
                )?)
            } else {
                None
            };
        if !anchor_matches && projection_blocks.is_none() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region output '{}' coordinate_system '{}' / genome_build '{}' is not compatible with sequence anchor genome_id '{}' (projection maps: {})",
                    output_dir,
                    coordinate_system,
                    genome_build,
                    anchor.genome_id,
                    if inspection.coordinate_projections.is_empty() {
                        "none".to_string()
                    } else {
                        inspection
                            .coordinate_projections
                            .iter()
                            .map(|projection| {
                                format!(
                                    "{}->{}",
                                    projection.source_genome_id, projection.target_genome_id
                                )
                            })
                            .collect::<Vec<_>>()
                            .join(", ")
                    }
                ),
                cause_chain: vec![],
            });
        }
        if let Some(min_abs_logfc) = min_abs_logfc
            && (!min_abs_logfc.is_finite() || min_abs_logfc < 0.0)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ProjectProbeRegionOutput min_abs_logfc must be >= 0".to_string(),
                cause_chain: vec![],
            });
        }
        let projection_level = Self::probe_region_projection_level(level)?;
        let feature_limit = max_features.unwrap_or(MAX_IMPORTED_SIGNAL_FEATURES);
        if feature_limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ProjectProbeRegionOutput max_features must be greater than zero"
                    .to_string(),
                cause_chain: vec![],
            });
        }

        let region_table_path = Path::new(output_dir).join(match projection_level.as_str() {
            "pm_probe" => PROBE_REGION_PROBE_TABLE_FILE,
            _ => PROBE_REGION_TABLE_FILE,
        });
        let (rows, available_contrasts, parse_warnings) =
            Self::probe_region_projection_rows_from_table(&region_table_path, &projection_level)?;
        if rows.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: if projection_level == "pm_probe" {
                    format!(
                        "Probe-region output '{}' has no projectable PM probe rows marked as probe_level_input",
                        output_dir
                    )
                } else {
                    format!(
                        "Probe-region output '{}' has no projectable rows",
                        output_dir
                    )
                },
                cause_chain: vec![],
            });
        }
        let requested_contrasts = contrasts
            .iter()
            .map(|contrast| contrast.trim())
            .filter(|contrast| !contrast.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>();
        let selected_contrasts = if requested_contrasts.is_empty() {
            available_contrasts.clone()
        } else {
            let mut selected = Vec::new();
            let mut missing = Vec::new();
            for requested in &requested_contrasts {
                if let Some(found) = available_contrasts
                    .iter()
                    .find(|contrast| contrast.eq_ignore_ascii_case(requested))
                {
                    selected.push(found.clone());
                } else {
                    missing.push(requested.clone());
                }
            }
            if !missing.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Probe-region output lacks requested log2FC contrast(s): {}",
                        missing.join(", ")
                    ),
                    cause_chain: vec![],
                });
            }
            selected
        };
        if selected_contrasts.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region output '{}' has no log2FC_* columns to project",
                    output_dir
                ),
                cause_chain: vec![],
            });
        }

        let mut report = MicroarrayProjectionReport {
            schema: MICROARRAY_PROJECTION_REPORT_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            manifest_path: output_dir.to_string(),
            dataset: "probe_region_output".to_string(),
            platform: inspection
                .platform
                .clone()
                .unwrap_or_else(|| "probe-region helper output".to_string()),
            normalization: inspection
                .normalization
                .clone()
                .unwrap_or_else(|| "not declared".to_string()),
            coordinate_system: coordinate_system.clone(),
            coordinate_projection_used: projection_blocks.is_some(),
            coordinate_projection_method: projection_spec.map(|spec| spec.method.clone()),
            coordinate_projection_path: projection_path.clone(),
            anchor_genome_id: anchor.genome_id.clone(),
            anchor_chromosome: anchor.chromosome.clone(),
            anchor_start_1based: anchor.start_1based,
            anchor_end_1based: anchor.end_1based,
            anchor_strand: anchor.strand.unwrap_or('+').to_string(),
            requested_contrasts,
            projected_contrasts: selected_contrasts.clone(),
            level: projection_level.clone(),
            parsed_rows: rows.len(),
            warnings: parse_warnings,
            ..Default::default()
        };

        let mut pending_features = Vec::new();
        let mut grouped_values: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut mismatch_counts: HashMap<String, usize> = HashMap::new();
        'rows: for mut row in rows {
            let native_row = row.clone();
            let projected_interval = if let Some(blocks) = projection_blocks.as_deref() {
                match Self::project_genome_interval_with_blocks(
                    blocks,
                    &native_row.chromosome,
                    native_row.start_1based,
                    native_row.end_1based,
                    native_row.strand,
                ) {
                    Some(projected) => {
                        row.chromosome = projected.target_chrom.clone();
                        row.start_1based = projected.target_start_1based;
                        row.end_1based = projected.target_end_1based;
                        row.strand = projected.target_strand;
                        Some(projected)
                    }
                    None => {
                        report.skipped_rows += selected_contrasts.len();
                        report.skipped_projection_unmapped += selected_contrasts.len();
                        continue;
                    }
                }
            } else {
                None
            };
            if !Self::chromosomes_match(&row.chromosome, &anchor.chromosome) {
                report.skipped_rows += selected_contrasts.len();
                report.skipped_wrong_chromosome += selected_contrasts.len();
                *mismatch_counts.entry(row.chromosome.clone()).or_insert(0) += 1;
                continue;
            }
            if row.end_1based < anchor.start_1based || row.start_1based > anchor.end_1based {
                report.skipped_rows += selected_contrasts.len();
                report.skipped_non_overlap += selected_contrasts.len();
                continue;
            }
            let overlap_start_1based = row.start_1based.max(anchor.start_1based);
            let overlap_end_1based = row.end_1based.min(anchor.end_1based);
            if overlap_end_1based < overlap_start_1based {
                report.skipped_rows += selected_contrasts.len();
                report.skipped_non_overlap += selected_contrasts.len();
                continue;
            }
            let (local_start_0based, local_end_0based_exclusive) = if anchor.strand == Some('-') {
                (
                    anchor.end_1based.saturating_sub(overlap_end_1based),
                    anchor.end_1based.saturating_sub(overlap_start_1based) + 1,
                )
            } else {
                (
                    overlap_start_1based - anchor.start_1based,
                    overlap_end_1based - anchor.start_1based + 1,
                )
            };
            let local_strand = match (row.strand, anchor.strand) {
                (Some('+'), Some('-')) => Some('-'),
                (Some('-'), Some('-')) => Some('+'),
                (Some(strand), _) => Some(strand),
                _ => None,
            };
            for contrast in &selected_contrasts {
                let Some(logfc) = row.logfc_values.get(contrast).copied() else {
                    report.skipped_rows += 1;
                    report.skipped_invalid += 1;
                    continue;
                };
                if min_abs_logfc
                    .map(|threshold| logfc.abs() < threshold)
                    .unwrap_or(false)
                {
                    report.skipped_rows += 1;
                    report.skipped_filter += 1;
                    continue;
                }
                if pending_features.len() >= feature_limit {
                    report.truncated_at_limit = true;
                    break 'rows;
                }
                grouped_values
                    .entry(row.feature_id.clone())
                    .or_default()
                    .push(format!(
                        "{} logFC={}",
                        contrast,
                        Self::probe_region_format_number(logfc)
                    ));
                pending_features.push(Self::build_probe_region_output_feature(
                    output_dir,
                    &inspection,
                    anchor,
                    contrast,
                    &native_row,
                    &row,
                    &projection_level,
                    projected_interval.as_ref(),
                    logfc,
                    local_start_0based,
                    local_end_0based_exclusive,
                    local_strand,
                ));
            }
        }

        if !mismatch_counts.is_empty() {
            let mut sorted = mismatch_counts
                .iter()
                .map(|(chromosome, count)| (chromosome.clone(), *count))
                .collect::<Vec<_>>();
            sorted.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
            let seen = sorted
                .iter()
                .take(3)
                .map(|(chromosome, count)| format!("{chromosome} ({count})"))
                .collect::<Vec<_>>()
                .join(", ");
            report.warnings.push(format!(
                "{} probe-region projection row(s) did not match anchor chromosome '{}' (examples: {})",
                report.skipped_wrong_chromosome, anchor.chromosome, seen
            ));
        }

        if clear_existing {
            Self::remove_generated_microarray_track_features(dna.features_mut());
        }
        report.imported_features = pending_features.len();
        for mut feature in pending_features {
            let feature_id = feature
                .qualifier_values("gentle_array_feature_id")
                .next()
                .map(str::to_string);
            if let Some(values) = feature_id.as_ref().and_then(|id| grouped_values.get(id)) {
                feature
                    .qualifiers
                    .push(("gentle_array_value_summary".into(), Some(values.join("; "))));
            }
            dna.features_mut().push(feature);
        }
        if report.truncated_at_limit {
            report.warnings.push(format!(
                "Probe-region output projection was truncated after {} features (limit={})",
                report.imported_features, feature_limit
            ));
        }

        Ok(report)
    }

    pub(super) fn probe_region_output_supports_anchor(
        coordinate_system: &str,
        genome_build: &str,
        anchor: &GenomeSequenceAnchor,
    ) -> bool {
        let anchor_id = anchor.genome_id.trim();
        !anchor_id.is_empty()
            && (coordinate_system.trim().eq_ignore_ascii_case(anchor_id)
                || genome_build.trim().eq_ignore_ascii_case(anchor_id))
    }

    pub(super) fn resolve_probe_region_coordinate_projection_spec<'a>(
        inspection: &'a ProbeRegionOutputInspection,
        anchor: &GenomeSequenceAnchor,
    ) -> Option<&'a GenomeCoordinateProjectionSpec> {
        let coordinate_system = inspection.coordinate_system.as_deref().unwrap_or("");
        let genome_build = inspection.genome_build.as_deref().unwrap_or("");
        inspection.coordinate_projections.iter().find(|projection| {
            (Self::genome_build_tokens_match(&projection.source_genome_id, coordinate_system)
                || Self::genome_build_tokens_match(&projection.source_genome_id, genome_build))
                && Self::genome_build_tokens_match(&projection.target_genome_id, &anchor.genome_id)
        })
    }

    pub(super) fn resolve_probe_region_output_path(output_dir: &str, path: &str) -> String {
        let path = Path::new(path);
        if path.is_absolute() {
            return path.to_string_lossy().to_string();
        }
        Path::new(output_dir)
            .join(path)
            .to_string_lossy()
            .to_string()
    }

    pub(super) fn probe_region_projection_rows_from_table(
        path: &Path,
        projection_level: &str,
    ) -> Result<(Vec<ProbeRegionProjectionRow>, Vec<String>, Vec<String>), EngineError> {
        let table_label = if projection_level == "pm_probe" {
            "Probe-intensity table"
        } else {
            "Probe-region table"
        };
        let mut reader = csv::ReaderBuilder::new()
            .trim(csv::Trim::All)
            .from_path(path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not open {} '{}': {e}",
                    table_label.to_ascii_lowercase(),
                    path.to_string_lossy()
                ),
                cause_chain: vec![],
            })?;
        let headers = reader
            .headers()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read {} header: {e}",
                    table_label.to_ascii_lowercase()
                ),
                cause_chain: vec![],
            })?
            .iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        let chromosome_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["chromosome", "chrom"])
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("{table_label} is missing chromosome column"),
                    cause_chain: vec![],
                })?;
        let start_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["start", "start_1based"])
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("{table_label} is missing start column"),
                    cause_chain: vec![],
                })?;
        let stop_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["stop", "end", "end_1based"],
        )
        .ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("{table_label} is missing stop column"),
            cause_chain: vec![],
        })?;
        let feature_idx = if projection_level == "pm_probe" {
            Self::probe_region_metadata_column_index(&headers, None, &["probe_id", "probe"])
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "Probe-intensity table is missing probe_id column".to_string(),
                    cause_chain: vec![],
                })?
        } else {
            Self::probe_region_metadata_column_index(
                &headers,
                None,
                &[
                    "probeset_or_region_id",
                    "feature_id",
                    "probeset_id",
                    "psr_id",
                ],
            )
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: "Probe-region table is missing probeset_or_region_id column".to_string(),
                cause_chain: vec![],
            })?
        };
        let parent_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &[
                "parent_probeset_or_region_id",
                "parent_probeset_id",
                "probeset_id",
            ],
        );
        let intensity_source_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["intensity_source"]);
        let strand_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["strand", "orientation"]);
        let transcript_idx = Self::probe_region_metadata_column_index(
            &headers,
            None,
            &["transcript_cluster_id", "transcript_cluster"],
        );
        let gene_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["gene_symbol", "gene"]);
        let logfc_indices = headers
            .iter()
            .enumerate()
            .filter(|(_, header)| header.starts_with("log2FC_"))
            .map(|(idx, header)| (idx, Self::probe_region_plot_track_label(header)))
            .collect::<Vec<_>>();
        if logfc_indices.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("{table_label} has no log2FC_* columns to project"),
                cause_chain: vec![],
            });
        }
        let available_contrasts = logfc_indices
            .iter()
            .map(|(_, contrast)| contrast.clone())
            .collect::<Vec<_>>();
        let mut rows = Vec::new();
        let mut warnings = Vec::new();
        let mut skipped_non_probe_input = 0usize;
        for (line_offset, record) in reader.records().enumerate() {
            let record = record.map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read {} row: {e}",
                    table_label.to_ascii_lowercase()
                ),
                cause_chain: vec![],
            })?;
            let line_no = line_offset + 2;
            let intensity_source = intensity_source_idx
                .and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string);
            if projection_level == "pm_probe"
                && intensity_source.as_deref() != Some("probe_level_input")
            {
                skipped_non_probe_input += 1;
                continue;
            }
            let chromosome = record.get(chromosome_idx).unwrap_or("").trim().to_string();
            let start_1based = record
                .get(start_idx)
                .unwrap_or("")
                .trim()
                .parse::<usize>()
                .map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Invalid start coordinate on probe-region line {line_no}: {e}"
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
                    message: format!("Invalid stop coordinate on probe-region line {line_no}: {e}"),
                    cause_chain: vec![],
                })?;
            if end_1based < start_1based {
                warnings.push(format!(
                    "Probe-region line {line_no} has stop before start; row skipped"
                ));
                continue;
            }
            let strand = strand_idx
                .and_then(|idx| record.get(idx))
                .and_then(|value| value.trim().chars().next())
                .filter(|value| matches!(value, '+' | '-'));
            let feature_id = record.get(feature_idx).unwrap_or("").trim().to_string();
            let parent_feature_id = parent_idx
                .and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string);
            let transcript_cluster_id = transcript_idx
                .and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string);
            let gene_symbol = gene_idx
                .and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string);
            let mut logfc_values = BTreeMap::new();
            for (idx, contrast) in &logfc_indices {
                if let Some(value) = record
                    .get(*idx)
                    .and_then(Self::probe_region_parse_plot_value)
                {
                    logfc_values.insert(contrast.clone(), value);
                }
            }
            rows.push(ProbeRegionProjectionRow {
                chromosome,
                start_1based,
                end_1based,
                strand,
                feature_id,
                parent_feature_id,
                intensity_source,
                transcript_cluster_id,
                gene_symbol,
                logfc_values,
            });
        }
        if skipped_non_probe_input > 0 {
            warnings.push(format!(
                "{} PM probe table row(s) were skipped because intensity_source was not probe_level_input",
                skipped_non_probe_input
            ));
        }
        Ok((rows, available_contrasts, warnings))
    }

    #[allow(clippy::too_many_arguments)]

    pub(super) fn build_probe_region_output_feature(
        output_dir: &str,
        inspection: &ProbeRegionOutputInspection,
        anchor: &GenomeSequenceAnchor,
        contrast: &str,
        native_row: &ProbeRegionProjectionRow,
        row: &ProbeRegionProjectionRow,
        projection_level: &str,
        projection: Option<&ProjectedGenomeInterval>,
        logfc: f64,
        local_start_0based: usize,
        local_end_0based_exclusive: usize,
        local_strand: Option<char>,
    ) -> gb_io::seq::Feature {
        let platform = inspection
            .platform
            .as_deref()
            .unwrap_or("probe-region helper output");
        let normalization = inspection
            .normalization
            .as_deref()
            .unwrap_or("not declared");
        let coordinate_system = inspection
            .coordinate_system
            .as_deref()
            .unwrap_or("not declared");
        let genome_build = inspection.genome_build.as_deref().unwrap_or("not declared");
        let assembly_check = if projection.is_some() {
            "projected_from_native_coordinate_system"
        } else {
            "helper_output_coordinate_system_matches_anchor"
        };
        let projection_status = projection
            .map(|projection| projection.status.as_str())
            .unwrap_or("direct_helper_output_coordinate_match");
        let gene_label = row.gene_symbol.as_deref().unwrap_or("");
        let label_suffix = if gene_label.is_empty() {
            row.feature_id.as_str()
        } else {
            gene_label
        };
        let mut qualifiers = vec![
            (
                "label".into(),
                Some(format!("{contrast} {}", label_suffix).trim().to_string()),
            ),
            (
                "note".into(),
                Some(format!(
                    "Probe-region helper output '{}' [{}:{}-{}]",
                    output_dir, row.chromosome, row.start_1based, row.end_1based
                )),
            ),
            ("gentle_track_source".into(), Some("Array".to_string())),
            (
                "gentle_generated".into(),
                Some(MICROARRAY_TRACK_GENERATED_TAG.to_string()),
            ),
            (
                "gentle_track_name".into(),
                Some(format!("{platform} {contrast}")),
            ),
            ("gentle_track_file".into(), Some(output_dir.to_string())),
            (
                "gentle_array_dataset".into(),
                Some("probe_region_output".to_string()),
            ),
            ("gentle_array_platform".into(), Some(platform.to_string())),
            (
                "gentle_array_normalization".into(),
                Some(normalization.to_string()),
            ),
            (
                "gentle_array_coordinate_system".into(),
                Some(coordinate_system.to_string()),
            ),
            (
                "gentle_array_supported_genome_ids".into(),
                Some(format!("{coordinate_system},{genome_build}")),
            ),
            (
                "gentle_array_anchor_genome_id".into(),
                Some(anchor.genome_id.clone()),
            ),
            (
                "gentle_array_anchor_chromosome".into(),
                Some(anchor.chromosome.clone()),
            ),
            (
                "gentle_array_anchor_start_1based".into(),
                Some(anchor.start_1based.to_string()),
            ),
            (
                "gentle_array_anchor_end_1based".into(),
                Some(anchor.end_1based.to_string()),
            ),
            (
                "gentle_array_anchor_strand".into(),
                Some(anchor.strand.unwrap_or('+').to_string()),
            ),
            (
                "gentle_array_assembly_check".into(),
                Some(assembly_check.to_string()),
            ),
            (
                "gentle_array_projection_status".into(),
                Some(projection_status.to_string()),
            ),
            ("gentle_array_contrast".into(), Some(contrast.to_string())),
            (
                "gentle_array_level".into(),
                Some(projection_level.to_string()),
            ),
            ("chromosome".into(), Some(row.chromosome.clone())),
            ("start_1based".into(), Some(row.start_1based.to_string())),
            ("end_1based".into(), Some(row.end_1based.to_string())),
            (
                "genomic_start_1based".into(),
                Some(row.start_1based.to_string()),
            ),
            (
                "genomic_end_1based".into(),
                Some(row.end_1based.to_string()),
            ),
            (
                "gentle_array_native_chromosome".into(),
                Some(native_row.chromosome.clone()),
            ),
            (
                "gentle_array_native_start_1based".into(),
                Some(native_row.start_1based.to_string()),
            ),
            (
                "gentle_array_native_end_1based".into(),
                Some(native_row.end_1based.to_string()),
            ),
        ];
        if let Some(projection) = projection {
            qualifiers.push((
                "gentle_array_projection_method".into(),
                Some(projection.method.clone()),
            ));
        }
        if !row.feature_id.trim().is_empty() {
            qualifiers.push(("feature_id".into(), Some(row.feature_id.clone())));
            qualifiers.push((
                "gentle_array_feature_id".into(),
                Some(row.feature_id.clone()),
            ));
        }
        if let Some(value) = &row.parent_feature_id {
            qualifiers.push(("gentle_array_parent_feature_id".into(), Some(value.clone())));
        }
        if let Some(value) = &row.intensity_source {
            qualifiers.push(("gentle_array_intensity_source".into(), Some(value.clone())));
        }
        if let Some(value) = &row.transcript_cluster_id {
            qualifiers.push(("transcript_cluster_id".into(), Some(value.clone())));
        }
        if let Some(value) = &row.gene_symbol {
            qualifiers.push(("gene".into(), Some(value.clone())));
            qualifiers.push(("gene_symbol".into(), Some(value.clone())));
        }
        if let Some(strand) = native_row.strand {
            qualifiers.push(("array_strand".into(), Some(strand.to_string())));
        }
        if let Some(strand) = local_strand {
            qualifiers.push(("strand".into(), Some(strand.to_string())));
        }
        let formatted = Self::probe_region_format_number(logfc);
        qualifiers.push(("logFC".into(), Some(formatted.clone())));
        qualifiers.push(("score".into(), Some(formatted)));

        let base_location = gb_io::seq::Location::simple_range(
            local_start_0based as i64,
            local_end_0based_exclusive as i64,
        );
        let location = if local_strand == Some('-') {
            gb_io::seq::Location::Complement(Box::new(base_location))
        } else {
            base_location
        };
        gb_io::seq::Feature {
            kind: "track".into(),
            location,
            qualifiers,
        }
    }

    pub(super) fn probe_region_format_number(value: f64) -> String {
        if value.abs() >= 100.0 || (value != 0.0 && value.abs() < 0.001) {
            format!("{value:.3e}")
        } else {
            format!("{value:.3}")
        }
    }

    pub(super) fn probe_region_projection_level(
        level: Option<&str>,
    ) -> Result<String, EngineError> {
        let normalized = level
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("probe_region")
            .replace('-', "_")
            .to_ascii_lowercase();
        match normalized.as_str() {
            "probe_region" | "region" | "probeset" | "probeset_region" | "psr" => {
                Ok("probe_region".to_string())
            }
            "pm_probe" | "probe" | "probe_level" | "pm" => Ok("pm_probe".to_string()),
            other => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "ProjectProbeRegionOutput level '{other}' is not supported; use probe_region or pm_probe"
                ),
                cause_chain: vec![],
            }),
        }
    }
}
