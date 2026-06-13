//! Affymetrix probe/probeset-region planning helpers.
//!
//! This first slice is intentionally plan-only: it checks user-supplied CEL,
//! metadata, annotation/library, platform, and local tool availability without
//! running R, APT, or any summarization backend.

use super::*;

const PROBE_REGION_STAGE: &str = "plan_only";
const PROBE_REGION_IMPLEMENTATION_STATUS: &str =
    "stage_1_preflight_only_no_cel_summarization_backend";
const PROBE_REGION_OLIGO_HELPER: &str = "scripts/probe_regions_oligo.R";

impl GentleEngine {
    /// Build a reusable preflight plan for chromosome-ordered probe inspection.
    pub fn plan_probe_regions(&self, request: ProbeRegionRequest) -> ProbeRegionPlan {
        let mut request = Self::normalize_probe_region_request(request);
        if request.normalization.is_empty() {
            request.normalization = "rma".to_string();
        }

        let mut warnings = Vec::new();
        let mut errors = Vec::new();
        let input_mode = Self::probe_region_input_mode(&request);
        if request.cel_paths.is_empty() && request.dataset.is_none() {
            errors.push(
                "arrays probe-regions requires at least one --cel path or a --dataset id"
                    .to_string(),
            );
        }
        if !request.cel_paths.is_empty() && request.dataset.is_some() {
            warnings.push(
                "Both explicit CEL files and --dataset were supplied; the current stage records both but does not discover dataset-owned CEL files"
                    .to_string(),
            );
        }
        if request.dataset.is_some() && request.cel_paths.is_empty() {
            warnings.push(
                "Publication-resource dataset discovery is planned but not executed in this stage"
                    .to_string(),
            );
        }
        if !Self::probe_region_has_selector(&request) {
            errors.push(
                "arrays probe-regions requires at least one --gene, --locus, --transcript-cluster-id, or --probeset-id selector"
                    .to_string(),
            );
        }

        let normalization = request.normalization.to_ascii_lowercase();
        match normalization.as_str() {
            "rma" | "none" => {}
            "quantile-feature" => warnings.push(
                "normalization=quantile-feature is accepted as a planned backend mode, but no execution backend is wired in this stage"
                    .to_string(),
            ),
            other => errors.push(format!(
                "Unsupported --normalization '{other}' (expected rma, quantile-feature, or none)"
            )),
        }

        let cel_files: Vec<_> = request
            .cel_paths
            .iter()
            .map(|path| Self::probe_region_file_status(path, "cel"))
            .collect();
        for cel in &cel_files {
            if !cel.exists {
                errors.push(format!("CEL file '{}' does not exist", cel.path));
            } else if !cel.is_file {
                errors.push(format!("CEL path '{}' is not a file", cel.path));
            }
        }

        let metadata = request
            .metadata_path
            .as_deref()
            .map(|path| Self::probe_region_file_status(path, "metadata"));
        if let Some(metadata) = &metadata {
            if !metadata.exists {
                errors.push(format!("Metadata file '{}' does not exist", metadata.path));
            } else if !metadata.is_file {
                errors.push(format!("Metadata path '{}' is not a file", metadata.path));
            }
        }
        let metadata_plan = Self::probe_region_metadata_plan(&request, metadata.as_ref());
        if let Some(plan) = &metadata_plan {
            errors.extend(
                plan.errors
                    .iter()
                    .map(|err| format!("Metadata preflight: {err}")),
            );
            warnings.extend(
                plan.warnings
                    .iter()
                    .map(|warning| format!("Metadata preflight: {warning}")),
            );
        }

        let platform = Self::probe_region_platform_plan(request.platform.as_deref());
        if platform.requested.is_none() {
            warnings.push(
                "No --platform was supplied; execution will need to infer chip type from CEL headers or a dataset manifest"
                    .to_string(),
            );
        } else if platform.confidence == "unknown" {
            warnings.push(format!(
                "Platform '{}' is not in GENtle's current Affymetrix backend map",
                platform.requested.clone().unwrap_or_default()
            ));
        }

        let mut dependencies = vec![
            Self::probe_region_command_dependency(
                "Rscript",
                "command",
                normalization != "none",
                &["--version"],
            ),
            Self::probe_region_command_dependency(
                "apt-probeset-summarize",
                "command",
                false,
                &["--version"],
            ),
        ];
        let rscript_available = dependencies
            .iter()
            .any(|row| row.name == "Rscript" && row.status == "present");
        dependencies.push(Self::probe_region_r_package_dependency(
            "oligo",
            normalization != "none",
            rscript_available,
        ));
        dependencies.push(Self::probe_region_r_package_dependency(
            "limma",
            normalization != "none",
            rscript_available,
        ));
        if let Some(package) = platform.bioconductor_package.as_deref() {
            dependencies.push(Self::probe_region_r_package_dependency(
                package,
                request.annotation_library_path.is_none(),
                rscript_available,
            ));
        }

        let annotation_source =
            Self::probe_region_annotation_source(&request, &platform, &dependencies);
        if let Some(path) = &annotation_source.path {
            if !path.exists {
                errors.push(format!(
                    "Annotation/library path '{}' does not exist",
                    path.path
                ));
            }
        } else if annotation_source.required_r_package.is_none() {
            errors.push(
                "No --annotation-library was supplied and no known Bioconductor platform package could be inferred; provide NetAffx/APT library files or --platform for a supported chip"
                    .to_string(),
            );
        }
        if !annotation_source.usable {
            warnings.push(
                annotation_source
                    .detail
                    .clone()
                    .unwrap_or_else(|| "No usable annotation source was confirmed".to_string()),
            );
        }

        if normalization != "none"
            && !dependencies
                .iter()
                .any(|row| row.name == "Rscript" && row.status == "present")
            && !dependencies
                .iter()
                .any(|row| row.name == "apt-probeset-summarize" && row.status == "present")
        {
            errors.push(
                "No local Affymetrix summarization backend was detected: install/provide Rscript with oligo or Affymetrix Power Tools, then rerun the preflight"
                    .to_string(),
            );
        }

        let output_status = request
            .output_dir
            .as_deref()
            .map(|path| Self::probe_region_file_status(path, "output_dir"));
        let cache_status = request
            .cache_dir
            .as_deref()
            .map(|path| Self::probe_region_file_status(path, "cache_dir"));
        if let Some(status) = &output_status {
            if status.exists && !status.is_dir {
                errors.push(format!(
                    "Output path '{}' exists but is not a directory",
                    status.path
                ));
            }
        }
        if let Some(status) = &cache_status {
            if status.exists && !status.is_dir {
                errors.push(format!(
                    "Cache path '{}' exists but is not a directory",
                    status.path
                ));
            }
        }
        let backend_candidates = Self::probe_region_backend_candidates(
            &request,
            &platform,
            &annotation_source,
            &dependencies,
            normalization.as_str(),
        );
        let contrasts = Self::probe_region_contrast_plan(metadata_plan.as_ref());

        let planned_outputs = vec![
            "region_intensity_chrom_order.csv".to_string(),
            "probe_intensity_chrom_order.csv when PM probe coordinates are available".to_string(),
            "plot.png and plot.svg when --plot is set".to_string(),
            "provenance.json".to_string(),
        ];
        let mut cache_compatibility_keys = vec![
            "CEL file path".to_string(),
            "CEL file size_bytes".to_string(),
            "CEL file modified_unix_seconds".to_string(),
            "normalization".to_string(),
            "annotation source path or Bioconductor package".to_string(),
            "platform".to_string(),
        ];
        if request.metadata_path.is_some() {
            cache_compatibility_keys.push("metadata path/size/mtime".to_string());
        }

        if !request.dry_run {
            warnings.push(
                "This command is currently plan-only; rerunning with --dry-run records the same preflight without implying execution"
                    .to_string(),
            );
        }

        let selectors = ProbeRegionSelectorPlan {
            genes: request.genes.clone(),
            loci: request.loci.clone(),
            transcript_cluster_ids: request.transcript_cluster_ids.clone(),
            probeset_ids: request.probeset_ids.clone(),
        };

        let preflight_ok = errors.is_empty();
        ProbeRegionPlan {
            schema: PROBE_REGION_PLAN_SCHEMA.to_string(),
            stage: PROBE_REGION_STAGE.to_string(),
            implementation_status: PROBE_REGION_IMPLEMENTATION_STATUS.to_string(),
            input_mode,
            request,
            selectors,
            cel_files,
            metadata,
            metadata_plan,
            annotation_source,
            platform,
            dependencies,
            backend_candidates,
            contrasts,
            output_dir_status: output_status,
            cache_dir_status: cache_status,
            planned_outputs,
            cache_compatibility_keys,
            warnings,
            errors,
            preflight_ok,
        }
    }

    fn normalize_probe_region_request(mut request: ProbeRegionRequest) -> ProbeRegionRequest {
        request.cel_paths = Self::probe_region_dedupe_nonempty(&request.cel_paths);
        request.dataset = Self::probe_region_trim_option(request.dataset);
        request.metadata_path = Self::probe_region_trim_option(request.metadata_path);
        request.genes = Self::probe_region_dedupe_nonempty(&request.genes);
        request.loci = Self::probe_region_dedupe_nonempty(&request.loci);
        request.transcript_cluster_ids =
            Self::probe_region_dedupe_nonempty(&request.transcript_cluster_ids);
        request.probeset_ids = Self::probe_region_dedupe_nonempty(&request.probeset_ids);
        request.platform = Self::probe_region_trim_option(request.platform);
        request.annotation_library_path =
            Self::probe_region_trim_option(request.annotation_library_path);
        request.condition_column = Self::probe_region_trim_option(request.condition_column);
        request.sample_column = Self::probe_region_trim_option(request.sample_column);
        request.block_column = Self::probe_region_trim_option(request.block_column);
        request.normalization = request.normalization.trim().to_string();
        request.output_dir = Self::probe_region_trim_option(request.output_dir);
        request.cache_dir = Self::probe_region_trim_option(request.cache_dir);
        request
    }

    fn probe_region_has_selector(request: &ProbeRegionRequest) -> bool {
        !request.genes.is_empty()
            || !request.loci.is_empty()
            || !request.transcript_cluster_ids.is_empty()
            || !request.probeset_ids.is_empty()
    }

    fn probe_region_input_mode(request: &ProbeRegionRequest) -> String {
        match (!request.cel_paths.is_empty(), request.dataset.is_some()) {
            (true, true) => "explicit_cel_plus_dataset".to_string(),
            (true, false) => "explicit_cel".to_string(),
            (false, true) => "publication_resource_dataset".to_string(),
            (false, false) => "missing".to_string(),
        }
    }

    fn probe_region_dedupe_nonempty(values: &[String]) -> Vec<String> {
        let mut out = Vec::new();
        let mut seen = BTreeSet::new();
        for value in values {
            let trimmed = value.trim();
            if trimmed.is_empty() {
                continue;
            }
            if seen.insert(trimmed.to_ascii_lowercase()) {
                out.push(trimmed.to_string());
            }
        }
        out
    }

    fn probe_region_trim_option(value: Option<String>) -> Option<String> {
        value
            .map(|text| text.trim().to_string())
            .filter(|text| !text.is_empty())
    }

    fn probe_region_file_status(path: &str, role: &str) -> ProbeRegionFileStatus {
        let trimmed = path.trim().to_string();
        match std::fs::metadata(&trimmed) {
            Ok(metadata) => ProbeRegionFileStatus {
                path: trimmed,
                role: role.to_string(),
                exists: true,
                is_file: metadata.is_file(),
                is_dir: metadata.is_dir(),
                size_bytes: Some(metadata.len()),
                modified_unix_seconds: metadata.modified().ok().and_then(|time| {
                    time.duration_since(std::time::UNIX_EPOCH)
                        .ok()
                        .map(|duration| duration.as_secs())
                }),
                detail: None,
            },
            Err(e) => ProbeRegionFileStatus {
                path: trimmed,
                role: role.to_string(),
                exists: false,
                is_file: false,
                is_dir: false,
                size_bytes: None,
                modified_unix_seconds: None,
                detail: Some(e.to_string()),
            },
        }
    }

    fn probe_region_platform_plan(platform: Option<&str>) -> ProbeRegionPlatformPlan {
        let requested = platform
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string);
        let normalized_key = requested
            .as_deref()
            .map(|value| {
                value
                    .chars()
                    .filter(|c| c.is_ascii_alphanumeric())
                    .flat_map(char::to_lowercase)
                    .collect::<String>()
            })
            .unwrap_or_default();
        match normalized_key.as_str() {
            "clariomdhuman" => ProbeRegionPlatformPlan {
                requested,
                normalized: "Clariom_D_Human".to_string(),
                backend_hint: "R/oligo".to_string(),
                bioconductor_package: Some("pd.clariom.d.human".to_string()),
                confidence: "known".to_string(),
            },
            "" => ProbeRegionPlatformPlan {
                requested,
                normalized: "unknown".to_string(),
                backend_hint: "infer_from_cel_or_dataset".to_string(),
                bioconductor_package: None,
                confidence: "unknown".to_string(),
            },
            _ => ProbeRegionPlatformPlan {
                normalized: requested.clone().unwrap_or_else(|| "unknown".to_string()),
                requested,
                backend_hint: "unknown".to_string(),
                bioconductor_package: None,
                confidence: "unknown".to_string(),
            },
        }
    }

    fn probe_region_annotation_source(
        request: &ProbeRegionRequest,
        platform: &ProbeRegionPlatformPlan,
        dependencies: &[ProbeRegionDependencyCheck],
    ) -> ProbeRegionAnnotationSourcePlan {
        let path = request
            .annotation_library_path
            .as_deref()
            .map(|path| Self::probe_region_file_status(path, "annotation_library"));
        let required_r_package = platform.bioconductor_package.clone();
        if let Some(path_status) = path {
            let source_kind = if path_status.is_dir {
                "annotation_directory"
            } else {
                Self::probe_region_annotation_file_kind(&path_status.path)
            };
            return ProbeRegionAnnotationSourcePlan {
                usable: path_status.exists,
                detail: if path_status.exists {
                    Some(format!(
                        "Using user-supplied {} '{}'",
                        source_kind, path_status.path
                    ))
                } else {
                    Some(format!(
                        "User-supplied annotation/library path '{}' could not be read",
                        path_status.path
                    ))
                },
                path: Some(path_status),
                source_kind: source_kind.to_string(),
                required_r_package,
            };
        }

        let package_usable = required_r_package.as_deref().is_some_and(|package| {
            dependencies
                .iter()
                .any(|row| row.name == package && row.status == "present")
        });
        ProbeRegionAnnotationSourcePlan {
            path: None,
            source_kind: if required_r_package.is_some() {
                "bioconductor_package".to_string()
            } else {
                "unspecified".to_string()
            },
            usable: package_usable,
            detail: required_r_package.as_ref().map(|package| {
                if package_usable {
                    format!("Using detected Bioconductor platform package '{package}'")
                } else {
                    format!(
                        "Bioconductor platform package '{package}' was not confirmed; provide --annotation-library or install the package before execution"
                    )
                }
            }),
            required_r_package,
        }
    }

    fn probe_region_annotation_file_kind(path: &str) -> &'static str {
        let lower = Path::new(path)
            .extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        match lower.as_str() {
            "pgf" | "clf" | "mps" => "apt_library_file",
            "zip" => "netaffx_zip",
            "csv" | "tsv" => "annotation_table",
            "sqlite" | "db" => "sqlite_annotation",
            _ => "annotation_file",
        }
    }

    fn probe_region_metadata_plan(
        request: &ProbeRegionRequest,
        metadata: Option<&ProbeRegionFileStatus>,
    ) -> Option<ProbeRegionMetadataPlan> {
        let metadata = metadata?;
        let delimiter_byte = Self::probe_region_metadata_delimiter(&metadata.path);
        let delimiter = if delimiter_byte == b',' {
            "comma"
        } else {
            "tab"
        }
        .to_string();
        if !metadata.exists || !metadata.is_file {
            return Some(ProbeRegionMetadataPlan {
                status: "unavailable".to_string(),
                delimiter,
                errors: vec![format!(
                    "Metadata path '{}' is not a readable file",
                    metadata.path
                )],
                ..Default::default()
            });
        }

        let mut reader = match csv::ReaderBuilder::new()
            .delimiter(delimiter_byte)
            .flexible(true)
            .trim(csv::Trim::All)
            .from_path(&metadata.path)
        {
            Ok(reader) => reader,
            Err(e) => {
                return Some(ProbeRegionMetadataPlan {
                    status: "parse_error".to_string(),
                    delimiter,
                    errors: vec![format!(
                        "Could not open metadata table '{}': {e}",
                        metadata.path
                    )],
                    ..Default::default()
                });
            }
        };

        let headers = match reader.headers() {
            Ok(headers) => headers.clone(),
            Err(e) => {
                return Some(ProbeRegionMetadataPlan {
                    status: "parse_error".to_string(),
                    delimiter,
                    errors: vec![format!(
                        "Could not read metadata header '{}': {e}",
                        metadata.path
                    )],
                    ..Default::default()
                });
            }
        };
        let columns: Vec<String> = headers.iter().map(str::to_string).collect();
        let sample_idx = Self::probe_region_metadata_column_index(
            &columns,
            request.sample_column.as_deref(),
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
        );
        let condition_idx = Self::probe_region_metadata_column_index(
            &columns,
            request.condition_column.as_deref(),
            &[
                "condition",
                "group",
                "treatment",
                "sample group",
                "factor value condition",
                "factor value treatment",
                "characteristics condition",
            ],
        );
        let block_idx = Self::probe_region_metadata_column_index(
            &columns,
            request.block_column.as_deref(),
            &["block", "batch", "replicate", "biological replicate"],
        );

        let mut warnings = Vec::new();
        let mut errors = Vec::new();
        if request.sample_column.is_some() && sample_idx.is_none() {
            errors.push(format!(
                "Requested sample column '{}' was not found",
                request.sample_column.as_deref().unwrap_or_default()
            ));
        } else if sample_idx.is_none() {
            warnings.push(
                "No sample/CEL column was inferred; sample-to-CEL validation was skipped"
                    .to_string(),
            );
        }
        if request.condition_column.is_some() && condition_idx.is_none() {
            errors.push(format!(
                "Requested condition column '{}' was not found",
                request.condition_column.as_deref().unwrap_or_default()
            ));
        } else if condition_idx.is_none() {
            warnings.push(
                "No condition column was inferred; group summaries and default contrasts will be unavailable"
                    .to_string(),
            );
        }
        if request.block_column.is_some() && block_idx.is_none() {
            errors.push(format!(
                "Requested block column '{}' was not found",
                request.block_column.as_deref().unwrap_or_default()
            ));
        }

        let mut row_count = 0usize;
        let mut sample_count = 0usize;
        let mut condition_counts: BTreeMap<String, usize> = BTreeMap::new();
        for record in reader.records() {
            let record = match record {
                Ok(record) => record,
                Err(e) => {
                    errors.push(format!(
                        "Could not parse metadata row {}: {e}",
                        row_count + 2
                    ));
                    continue;
                }
            };
            row_count += 1;
            if let Some(idx) = sample_idx
                && record
                    .get(idx)
                    .map(str::trim)
                    .is_some_and(|v| !v.is_empty())
            {
                sample_count += 1;
            }
            if let Some(idx) = condition_idx
                && let Some(condition) = record.get(idx).map(str::trim).filter(|v| !v.is_empty())
            {
                *condition_counts.entry(condition.to_string()).or_insert(0) += 1;
            }
        }
        if sample_idx.is_some() && sample_count == 0 && row_count > 0 {
            warnings
                .push("Sample column was present but no non-empty samples were found".to_string());
        }
        if condition_idx.is_some() && condition_counts.is_empty() && row_count > 0 {
            warnings.push(
                "Condition column was present but no non-empty conditions were found".to_string(),
            );
        }

        Some(ProbeRegionMetadataPlan {
            status: if errors.is_empty() {
                "parsed".to_string()
            } else {
                "parse_error".to_string()
            },
            delimiter,
            columns,
            row_count,
            sample_column: sample_idx.map(|idx| headers.get(idx).unwrap_or_default().to_string()),
            condition_column: condition_idx
                .map(|idx| headers.get(idx).unwrap_or_default().to_string()),
            block_column: block_idx.map(|idx| headers.get(idx).unwrap_or_default().to_string()),
            sample_count,
            conditions: condition_counts
                .into_iter()
                .map(|(condition, sample_count)| ProbeRegionConditionSummary {
                    condition,
                    sample_count,
                })
                .collect(),
            warnings,
            errors,
        })
    }

    fn probe_region_metadata_delimiter(path: &str) -> u8 {
        match Path::new(path)
            .extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("")
            .to_ascii_lowercase()
            .as_str()
        {
            "csv" => b',',
            _ => b'\t',
        }
    }

    fn probe_region_metadata_column_index(
        columns: &[String],
        requested: Option<&str>,
        aliases: &[&str],
    ) -> Option<usize> {
        if let Some(requested) = requested {
            let requested_key = Self::probe_region_header_key(requested);
            return columns
                .iter()
                .position(|column| Self::probe_region_header_key(column) == requested_key);
        }
        let alias_keys: BTreeSet<String> = aliases
            .iter()
            .map(|alias| Self::probe_region_header_key(alias))
            .collect();
        columns
            .iter()
            .position(|column| alias_keys.contains(&Self::probe_region_header_key(column)))
    }

    fn probe_region_header_key(value: &str) -> String {
        value
            .chars()
            .filter(|ch| ch.is_ascii_alphanumeric())
            .flat_map(char::to_lowercase)
            .collect()
    }

    fn probe_region_contrast_plan(
        metadata_plan: Option<&ProbeRegionMetadataPlan>,
    ) -> Vec<ProbeRegionContrastPlan> {
        let Some(metadata_plan) = metadata_plan else {
            return Vec::new();
        };
        if metadata_plan.conditions.len() < 2 {
            return Vec::new();
        }
        let baseline = &metadata_plan.conditions[0];
        metadata_plan
            .conditions
            .iter()
            .skip(1)
            .map(|condition| ProbeRegionContrastPlan {
                contrast: format!("{}-{}", condition.condition, baseline.condition),
                numerator_condition: condition.condition.clone(),
                denominator_condition: baseline.condition.clone(),
                numerator_sample_count: condition.sample_count,
                denominator_sample_count: baseline.sample_count,
                status: if condition.sample_count > 0 && baseline.sample_count > 0 {
                    "planned".to_string()
                } else {
                    "insufficient_samples".to_string()
                },
            })
            .collect()
    }

    fn probe_region_backend_candidates(
        request: &ProbeRegionRequest,
        platform: &ProbeRegionPlatformPlan,
        annotation_source: &ProbeRegionAnnotationSourcePlan,
        dependencies: &[ProbeRegionDependencyCheck],
        normalization: &str,
    ) -> Vec<ProbeRegionBackendCandidate> {
        let rscript_present = Self::probe_region_dependency_present(dependencies, "Rscript");
        let oligo_present = Self::probe_region_dependency_present(dependencies, "oligo");
        let limma_present = Self::probe_region_dependency_present(dependencies, "limma");
        let apt_present =
            Self::probe_region_dependency_present(dependencies, "apt-probeset-summarize");
        let platform_package_present = platform
            .bioconductor_package
            .as_deref()
            .is_some_and(|package| Self::probe_region_dependency_present(dependencies, package));
        let annotation_path_usable = annotation_source
            .path
            .as_ref()
            .is_some_and(|path| path.exists);
        let r_oligo_normalization_supported = normalization == "rma";

        let mut r_missing = Vec::new();
        if !r_oligo_normalization_supported {
            r_missing.push(format!("normalization=rma (requested {normalization})"));
        }
        if !rscript_present {
            r_missing.push("Rscript command".to_string());
        }
        if !oligo_present && normalization != "none" {
            r_missing.push("R package oligo".to_string());
        }
        if !limma_present && normalization != "none" {
            r_missing.push("R package limma".to_string());
        }
        if platform.bioconductor_package.is_none() {
            r_missing.push("known Bioconductor platform package mapping".to_string());
        } else if !platform_package_present {
            r_missing.push(
                platform
                    .bioconductor_package
                    .clone()
                    .unwrap_or_else(|| "Bioconductor platform package".to_string()),
            );
        }

        let mut apt_missing = Vec::new();
        if !apt_present {
            apt_missing.push("apt-probeset-summarize command".to_string());
        }
        if !annotation_path_usable {
            apt_missing.push("user-supplied APT/NetAffx annotation-library path".to_string());
        }
        if request.cel_paths.is_empty() {
            apt_missing.push("explicit CEL file paths".to_string());
        }

        vec![
            ProbeRegionBackendCandidate {
                backend: "r_oligo".to_string(),
                status: if r_missing.is_empty() {
                    "ready".to_string()
                } else {
                    "missing_inputs".to_string()
                },
                required_inputs: vec![
                    "Rscript command".to_string(),
                    "R package oligo".to_string(),
                    "R package limma".to_string(),
                    "Bioconductor platform package".to_string(),
                    "normalization=rma".to_string(),
                ],
                missing: r_missing,
                helper_script: Some(PROBE_REGION_OLIGO_HELPER.to_string()),
                suggested_command: Self::probe_region_oligo_suggested_command(
                    request,
                    platform,
                    normalization,
                ),
                detail: Some(
                    "Preferred first execution backend for Clariom/whole-transcript arrays"
                        .to_string(),
                ),
            },
            ProbeRegionBackendCandidate {
                backend: "affymetrix_power_tools".to_string(),
                status: if apt_missing.is_empty() {
                    "ready".to_string()
                } else {
                    "missing_inputs".to_string()
                },
                required_inputs: vec![
                    "apt-probeset-summarize command".to_string(),
                    "APT/NetAffx annotation-library path".to_string(),
                    "explicit CEL file paths".to_string(),
                ],
                missing: apt_missing,
                helper_script: None,
                suggested_command: None,
                detail: Some(
                    "Useful when PGF/CLF/MPS or compatible vendor libraries are supplied"
                        .to_string(),
                ),
            },
            ProbeRegionBackendCandidate {
                backend: "plan_only".to_string(),
                status: "ready".to_string(),
                required_inputs: vec![],
                missing: vec![],
                helper_script: None,
                suggested_command: None,
                detail: Some(
                    "Current implemented mode: preflight without CEL summarization".to_string(),
                ),
            },
        ]
    }

    fn probe_region_oligo_suggested_command(
        request: &ProbeRegionRequest,
        platform: &ProbeRegionPlatformPlan,
        normalization: &str,
    ) -> Option<String> {
        if request.cel_paths.is_empty() || normalization != "rma" {
            return None;
        }
        let mut args = vec!["Rscript".to_string(), PROBE_REGION_OLIGO_HELPER.to_string()];
        for cel in &request.cel_paths {
            args.push("--cel".to_string());
            args.push(cel.clone());
        }
        if let Some(metadata) = &request.metadata_path {
            args.push("--metadata".to_string());
            args.push(metadata.clone());
        }
        if let Some(sample_column) = &request.sample_column {
            args.push("--sample-column".to_string());
            args.push(sample_column.clone());
        }
        if let Some(condition_column) = &request.condition_column {
            args.push("--condition-column".to_string());
            args.push(condition_column.clone());
        }
        if let Some(block_column) = &request.block_column {
            args.push("--block-column".to_string());
            args.push(block_column.clone());
        }
        if let Some(package) = &platform.bioconductor_package {
            args.push("--platform-package".to_string());
            args.push(package.clone());
        }
        args.push("--normalization".to_string());
        args.push(normalization.to_string());
        args.push("--output".to_string());
        args.push(
            request
                .output_dir
                .clone()
                .unwrap_or_else(|| "analysis/probe_regions".to_string()),
        );
        if let Some(cache_dir) = &request.cache_dir {
            args.push("--cache-dir".to_string());
            args.push(cache_dir.clone());
        }
        for gene in &request.genes {
            args.push("--gene".to_string());
            args.push(gene.clone());
        }
        for locus in &request.loci {
            args.push("--locus".to_string());
            args.push(locus.clone());
        }
        for transcript_cluster_id in &request.transcript_cluster_ids {
            args.push("--transcript-cluster-id".to_string());
            args.push(transcript_cluster_id.clone());
        }
        for probeset_id in &request.probeset_ids {
            args.push("--probeset-id".to_string());
            args.push(probeset_id.clone());
        }
        Some(
            args.iter()
                .map(|arg| Self::probe_region_shell_quote(arg))
                .collect::<Vec<_>>()
                .join(" "),
        )
    }

    fn probe_region_shell_quote(value: &str) -> String {
        if value
            .chars()
            .all(|ch| ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.' | '/' | ':' | '='))
        {
            return value.to_string();
        }
        format!("'{}'", value.replace('\'', "'\"'\"'"))
    }

    fn probe_region_dependency_present(
        dependencies: &[ProbeRegionDependencyCheck],
        name: &str,
    ) -> bool {
        dependencies
            .iter()
            .any(|row| row.name == name && row.status == "present")
    }

    fn probe_region_command_dependency(
        name: &str,
        kind: &str,
        required: bool,
        version_args: &[&str],
    ) -> ProbeRegionDependencyCheck {
        match Command::new(name).args(version_args).output() {
            Ok(output) => {
                let text = if output.stderr.is_empty() {
                    String::from_utf8_lossy(&output.stdout).to_string()
                } else {
                    String::from_utf8_lossy(&output.stderr).to_string()
                };
                let version = text
                    .lines()
                    .map(str::trim)
                    .find(|line| !line.is_empty())
                    .map(str::to_string);
                ProbeRegionDependencyCheck {
                    name: name.to_string(),
                    kind: kind.to_string(),
                    required,
                    status: if output.status.success() {
                        "present".to_string()
                    } else {
                        "missing".to_string()
                    },
                    version,
                    detail: if output.status.success() {
                        None
                    } else {
                        Some(format!("{name} exited with status {}", output.status))
                    },
                }
            }
            Err(e) => ProbeRegionDependencyCheck {
                name: name.to_string(),
                kind: kind.to_string(),
                required,
                status: "missing".to_string(),
                version: None,
                detail: Some(e.to_string()),
            },
        }
    }

    fn probe_region_r_package_dependency(
        package: &str,
        required: bool,
        rscript_available: bool,
    ) -> ProbeRegionDependencyCheck {
        if !rscript_available {
            return ProbeRegionDependencyCheck {
                name: package.to_string(),
                kind: "r_package".to_string(),
                required,
                status: "unchecked".to_string(),
                version: None,
                detail: Some(
                    "Rscript was not available, so the R package check was skipped".to_string(),
                ),
            };
        }
        let expression = format!(
            "quit(status = if (requireNamespace({:?}, quietly = TRUE)) 0 else 1)",
            package
        );
        match Command::new("Rscript").args(["-e", &expression]).output() {
            Ok(output) if output.status.success() => ProbeRegionDependencyCheck {
                name: package.to_string(),
                kind: "r_package".to_string(),
                required,
                status: "present".to_string(),
                version: None,
                detail: None,
            },
            Ok(output) => ProbeRegionDependencyCheck {
                name: package.to_string(),
                kind: "r_package".to_string(),
                required,
                status: "missing".to_string(),
                version: None,
                detail: Some(format!(
                    "R package check exited with status {}",
                    output.status
                )),
            },
            Err(e) => ProbeRegionDependencyCheck {
                name: package.to_string(),
                kind: "r_package".to_string(),
                required,
                status: "missing".to_string(),
                version: None,
                detail: Some(e.to_string()),
            },
        }
    }
}
