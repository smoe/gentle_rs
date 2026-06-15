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
const PROBE_REGION_TABLE_FILE: &str = "region_intensity_chrom_order.csv";
const PROBE_REGION_SAMPLE_TABLE_FILE: &str = "sample_table.tsv";
const PROBE_REGION_MATRIX_MANIFEST_FILE: &str = "normalized_feature_matrix_manifest.json";
const PROBE_REGION_PROVENANCE_FILE: &str = "provenance.json";
const PROBE_REGION_MATRIX_MANIFEST_SCHEMA: &str =
    "gentle.probe_region_normalized_matrix_manifest.v1";
const PROBE_REGION_BACKEND_PROVENANCE_SCHEMA: &str = "gentle.probe_region_backend_provenance.v1";
const CLARIOM_D_HUMAN_VENDOR_SUPPORT_DIR: &str =
    "data/resources/affymetrix/clariom_d_human_na36_hg38";
const CLARIOM_D_HUMAN_PROBESET_ZIP: &str = "Clariom_D_Human-na36-hg38-probeset-csv.zip";
const CLARIOM_D_HUMAN_TRANSCRIPT_ZIP: &str = "Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip";

#[derive(Default)]
struct ProbeRegionTableSummary {
    row_count: usize,
    column_count: usize,
    feature_count: usize,
    transcript_cluster_count: usize,
    chromosome_count: usize,
    chromosomes: Vec<String>,
    gene_symbols: Vec<String>,
    sample_columns: Vec<String>,
    condition_summary_columns: Vec<String>,
    logfc_columns: Vec<String>,
    preview_rows: Vec<ProbeRegionOutputPreviewRow>,
    required_columns_missing: Vec<String>,
    warnings: Vec<String>,
}

#[derive(Clone)]
struct ProbeRegionPlotTrack {
    label: String,
    values: Vec<Option<f64>>,
}

#[derive(Clone)]
struct ProbeRegionPlotRow {
    chromosome: String,
    start_1based: usize,
    stop_1based: Option<usize>,
    probeset_or_region_id: String,
    gene_symbol: String,
}

struct ProbeRegionPlotData {
    rows: Vec<ProbeRegionPlotRow>,
    intensity_tracks: Vec<ProbeRegionPlotTrack>,
    logfc_tracks: Vec<ProbeRegionPlotTrack>,
    chromosome_count: usize,
    warnings: Vec<String>,
}

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

    /// Inspect outputs written by the explicit R/oligo probe-region helper.
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
        let sample_table_path = output_path.join(PROBE_REGION_SAMPLE_TABLE_FILE);
        let manifest_path = output_path.join(PROBE_REGION_MATRIX_MANIFEST_FILE);
        let provenance_path = output_path.join(PROBE_REGION_PROVENANCE_FILE);

        let region_table =
            Self::probe_region_file_status(&region_table_path.to_string_lossy(), "region_table");
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
            sample_table,
            normalized_matrix_manifest,
            provenance,
            backend,
            platform,
            platform_package,
            normalization,
            coordinate_system,
            genome_build,
            projection_ready,
            projection_blockers,
            target_levels,
            artifact_paths,
            row_count: table_summary.row_count,
            column_count: table_summary.column_count,
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
    pub fn export_probe_region_output_svg(
        &self,
        output_dir: &str,
        svg_path: &str,
    ) -> Result<ProbeRegionOutputSvgExport, EngineError> {
        let output_dir = output_dir.trim();
        let svg_path = svg_path.trim();
        if output_dir.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Probe-region output directory must not be empty".to_string(),
                cause_chain: vec![],
            });
        }
        if svg_path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Probe-region SVG output path must not be empty".to_string(),
                cause_chain: vec![],
            });
        }

        let inspection = self.inspect_probe_region_output(output_dir)?;
        if !inspection.usable {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region output '{}' is not usable: {}",
                    output_dir,
                    inspection.errors.join("; ")
                ),
                cause_chain: vec![],
            });
        }

        let region_table_path = Path::new(output_dir).join(PROBE_REGION_TABLE_FILE);
        let plot_data = Self::probe_region_plot_data_from_table(&region_table_path)?;
        let svg = Self::render_probe_region_output_svg_text(&inspection, &plot_data);
        if let Some(parent) = Path::new(svg_path).parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not create probe-region SVG output directory '{}': {e}",
                        parent.to_string_lossy()
                    ),
                    cause_chain: vec![],
                })?;
            }
        }
        std::fs::write(svg_path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write probe-region SVG '{}': {e}", svg_path),
            cause_chain: vec![],
        })?;

        let mut warnings = inspection.warnings.clone();
        warnings.extend(plot_data.warnings.clone());
        warnings.extend(inspection.projection_blockers.iter().map(|blocker| {
            format!("Projection remains blocked for genome-anchored feature projection: {blocker}")
        }));

        Ok(ProbeRegionOutputSvgExport {
            schema: PROBE_REGION_OUTPUT_SVG_EXPORT_SCHEMA.to_string(),
            output_dir: output_dir.to_string(),
            svg_path: svg_path.to_string(),
            row_count: plot_data.rows.len(),
            intensity_track_count: plot_data.intensity_tracks.len(),
            logfc_track_count: plot_data.logfc_tracks.len(),
            chromosome_count: plot_data.chromosome_count,
            projection_ready: inspection.projection_ready,
            warnings,
        })
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

    fn probe_region_optional_file_status(path: &Path, role: &str) -> Option<ProbeRegionFileStatus> {
        let status = Self::probe_region_file_status(&path.to_string_lossy(), role);
        status.exists.then_some(status)
    }

    fn probe_region_json_file(path: &Path) -> Result<serde_json::Value, String> {
        let path_text = path.to_string_lossy();
        let file = File::open(path).map_err(|e| format!("could not open '{path_text}': {e}"))?;
        serde_json::from_reader(BufReader::new(file))
            .map_err(|e| format!("could not parse '{path_text}' as JSON: {e}"))
    }

    fn probe_region_json_string(value: &serde_json::Value, key: &str) -> Option<String> {
        value
            .get(key)
            .and_then(|value| value.as_str())
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
    }

    fn probe_region_json_string_array(value: &serde_json::Value, key: &str) -> Vec<String> {
        value
            .get(key)
            .and_then(|value| value.as_array())
            .map(|values| {
                values
                    .iter()
                    .filter_map(|value| value.as_str())
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(str::to_string)
                    .collect()
            })
            .unwrap_or_default()
    }

    fn probe_region_extend_unique(target: &mut Vec<String>, values: Vec<String>) {
        let mut seen = target
            .iter()
            .map(|value| value.to_ascii_lowercase())
            .collect::<BTreeSet<_>>();
        for value in values {
            if seen.insert(value.to_ascii_lowercase()) {
                target.push(value);
            }
        }
    }

    fn probe_region_warn_unexpected_json_schema(
        value: &serde_json::Value,
        expected_schema: &str,
        warnings: &mut Vec<String>,
        label: &str,
    ) {
        match Self::probe_region_json_string(value, "schema") {
            Some(schema) if schema == expected_schema => {}
            Some(schema) => warnings.push(format!(
                "{} uses schema '{}' but '{}' was expected",
                label, schema, expected_schema
            )),
            None => warnings.push(format!("{label} does not declare a schema")),
        }
    }

    fn probe_region_sample_ids_from_sample_table(path: &Path) -> Result<Vec<String>, String> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .trim(csv::Trim::All)
            .from_path(path)
            .map_err(|e| format!("could not open '{}': {e}", path.to_string_lossy()))?;
        let headers = reader
            .headers()
            .map_err(|e| format!("could not read header: {e}"))?
            .iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        let sample_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["sample_id", "sample"])
                .ok_or_else(|| {
                    "sample_table.tsv does not include a sample_id column".to_string()
                })?;
        let mut sample_ids = Vec::new();
        let mut seen = BTreeSet::new();
        for record in reader.records() {
            let record = record.map_err(|e| format!("could not read sample row: {e}"))?;
            let sample_id = record.get(sample_idx).unwrap_or("").trim();
            if !sample_id.is_empty() && seen.insert(sample_id.to_ascii_lowercase()) {
                sample_ids.push(sample_id.to_string());
            }
        }
        Ok(sample_ids)
    }

    fn inspect_probe_region_table(
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

    fn probe_region_plot_data_from_table(path: &Path) -> Result<ProbeRegionPlotData, EngineError> {
        let mut reader = csv::ReaderBuilder::new()
            .trim(csv::Trim::All)
            .from_path(path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not open probe-region table '{}': {e}",
                    path.to_string_lossy()
                ),
                cause_chain: vec![],
            })?;
        let headers = reader
            .headers()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not read probe-region table header: {e}"),
                cause_chain: vec![],
            })?
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
        let gene_idx =
            Self::probe_region_metadata_column_index(&headers, None, &["gene_symbol", "gene"]);
        let intensity_indices = headers
            .iter()
            .enumerate()
            .filter(|(_, header)| header.starts_with("mean_log2_"))
            .map(|(idx, header)| (idx, Self::probe_region_plot_track_label(header)))
            .collect::<Vec<_>>();
        let logfc_indices = headers
            .iter()
            .enumerate()
            .filter(|(_, header)| header.starts_with("log2FC_"))
            .map(|(idx, header)| (idx, Self::probe_region_plot_track_label(header)))
            .collect::<Vec<_>>();

        let mut rows = Vec::new();
        let mut intensity_tracks = intensity_indices
            .iter()
            .map(|(_, label)| ProbeRegionPlotTrack {
                label: label.clone(),
                values: Vec::new(),
            })
            .collect::<Vec<_>>();
        let mut logfc_tracks = logfc_indices
            .iter()
            .map(|(_, label)| ProbeRegionPlotTrack {
                label: label.clone(),
                values: Vec::new(),
            })
            .collect::<Vec<_>>();
        let mut chromosomes = BTreeSet::new();
        let mut warnings = Vec::new();

        for record in reader.records() {
            let record = record.map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not read probe-region row: {e}"),
                cause_chain: vec![],
            })?;
            let row_index = rows.len();
            let chromosome = chromosome_idx
                .and_then(|idx| record.get(idx))
                .unwrap_or("")
                .trim()
                .to_string();
            if !chromosome.is_empty() {
                chromosomes.insert(chromosome.clone());
            }
            let start_1based = start_idx
                .and_then(|idx| record.get(idx))
                .and_then(|value| value.trim().parse::<usize>().ok())
                .unwrap_or_else(|| row_index + 1);
            if start_idx.is_some() && start_1based == row_index + 1 {
                let raw = start_idx
                    .and_then(|idx| record.get(idx))
                    .unwrap_or("")
                    .trim();
                if !raw.is_empty() && raw.parse::<usize>().is_err() {
                    warnings.push(format!(
                        "Could not parse start coordinate '{}' on row {}; plotted by row order",
                        raw,
                        row_index + 1
                    ));
                }
            }
            let stop_1based = stop_idx
                .and_then(|idx| record.get(idx))
                .and_then(|value| value.trim().parse::<usize>().ok());
            let probeset_or_region_id = feature_idx
                .and_then(|idx| record.get(idx))
                .unwrap_or("")
                .trim()
                .to_string();
            let gene_symbol = gene_idx
                .and_then(|idx| record.get(idx))
                .unwrap_or("")
                .trim()
                .to_string();
            for (track, (idx, _)) in intensity_tracks.iter_mut().zip(intensity_indices.iter()) {
                track.values.push(
                    record
                        .get(*idx)
                        .and_then(Self::probe_region_parse_plot_value),
                );
            }
            for (track, (idx, _)) in logfc_tracks.iter_mut().zip(logfc_indices.iter()) {
                track.values.push(
                    record
                        .get(*idx)
                        .and_then(Self::probe_region_parse_plot_value),
                );
            }
            rows.push(ProbeRegionPlotRow {
                chromosome,
                start_1based,
                stop_1based,
                probeset_or_region_id,
                gene_symbol,
            });
        }

        if rows.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region table '{}' has no data rows",
                    path.to_string_lossy()
                ),
                cause_chain: vec![],
            });
        }
        if intensity_tracks.is_empty() {
            warnings.push(
                "No mean_log2_* columns found; upper intensity panel will be annotated as unavailable"
                    .to_string(),
            );
        }
        if logfc_tracks.is_empty() {
            warnings.push(
                "No log2FC_* columns found; lower fold-change panel will be annotated as unavailable"
                    .to_string(),
            );
        }
        if chromosomes.len() > 1 {
            warnings.push(
                "Multiple chromosomes are present; x-axis uses chromosome-ordered row positions"
                    .to_string(),
            );
        }

        Ok(ProbeRegionPlotData {
            rows,
            intensity_tracks,
            logfc_tracks,
            chromosome_count: chromosomes.len(),
            warnings,
        })
    }

    fn probe_region_plot_track_label(header: &str) -> String {
        header
            .strip_prefix("mean_log2_")
            .or_else(|| header.strip_prefix("log2FC_"))
            .unwrap_or(header)
            .replace('_', " ")
    }

    fn probe_region_parse_plot_value(value: &str) -> Option<f64> {
        let parsed = value.trim().parse::<f64>().ok()?;
        parsed.is_finite().then_some(parsed)
    }

    fn probe_region_plot_range(
        tracks: &[ProbeRegionPlotTrack],
        include_zero: bool,
    ) -> Option<(f64, f64)> {
        let mut min_value = if include_zero { Some(0.0) } else { None };
        let mut max_value = if include_zero { Some(0.0) } else { None };
        for value in tracks
            .iter()
            .flat_map(|track| track.values.iter().filter_map(|value| *value))
        {
            min_value = Some(min_value.map_or(value, |current: f64| current.min(value)));
            max_value = Some(max_value.map_or(value, |current: f64| current.max(value)));
        }
        let (mut min_value, mut max_value) = (min_value?, max_value?);
        if (max_value - min_value).abs() < f64::EPSILON {
            min_value -= 1.0;
            max_value += 1.0;
        } else {
            let pad = (max_value - min_value) * 0.08;
            min_value -= pad;
            max_value += pad;
        }
        Some((min_value, max_value))
    }

    fn probe_region_plot_x(
        row_idx: usize,
        row: &ProbeRegionPlotRow,
        row_count: usize,
        use_genomic_x: bool,
        min_start: usize,
        max_start: usize,
        left: f64,
        width: f64,
    ) -> f64 {
        if row_count <= 1 {
            return left + width / 2.0;
        }
        let frac = if use_genomic_x && max_start > min_start {
            (row.start_1based.saturating_sub(min_start)) as f64 / (max_start - min_start) as f64
        } else {
            row_idx as f64 / (row_count - 1) as f64
        };
        left + width * frac.clamp(0.0, 1.0)
    }

    fn probe_region_plot_y(
        value: f64,
        min_value: f64,
        max_value: f64,
        top: f64,
        height: f64,
    ) -> f64 {
        let frac = ((value - min_value) / (max_value - min_value)).clamp(0.0, 1.0);
        top + height - (height * frac)
    }

    fn probe_region_svg_escape(raw: &str) -> String {
        raw.replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
            .replace('"', "&quot;")
    }

    fn render_probe_region_output_svg_text(
        inspection: &ProbeRegionOutputInspection,
        plot: &ProbeRegionPlotData,
    ) -> String {
        use std::fmt::Write as _;

        let width = 1120.0;
        let height = 680.0;
        let left = 90.0;
        let right = 40.0;
        let plot_width = width - left - right;
        let top_panel_top = 125.0;
        let panel_height = 210.0;
        let bottom_panel_top = 405.0;
        let palette = [
            "#2563eb", "#dc2626", "#16a34a", "#9333ea", "#d97706", "#0891b2", "#be123c", "#4f46e5",
        ];
        let min_start = plot
            .rows
            .iter()
            .map(|row| row.start_1based)
            .min()
            .unwrap_or(0);
        let max_start = plot
            .rows
            .iter()
            .map(|row| row.start_1based)
            .max()
            .unwrap_or(min_start);
        let use_genomic_x = plot.chromosome_count <= 1 && max_start > min_start;
        let intensity_range = Self::probe_region_plot_range(&plot.intensity_tracks, false);
        let logfc_range = Self::probe_region_plot_range(&plot.logfc_tracks, true);
        let target_label = if inspection.gene_symbols.is_empty() {
            "probe-region output".to_string()
        } else {
            inspection.gene_symbols.join(", ")
        };
        let provenance = format!(
            "backend={} platform={} normalization={} coordinate_system={} genome_build={}",
            inspection.backend.as_deref().unwrap_or("not declared"),
            inspection.platform.as_deref().unwrap_or("not declared"),
            inspection
                .normalization
                .as_deref()
                .unwrap_or("not declared"),
            inspection
                .coordinate_system
                .as_deref()
                .unwrap_or("not declared"),
            inspection.genome_build.as_deref().unwrap_or("not declared")
        );

        let mut svg = String::new();
        let _ = writeln!(
            svg,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width:.0}\" height=\"{height:.0}\" viewBox=\"0 0 {width:.0} {height:.0}\" data-gentle-schema=\"{}\">",
            Self::probe_region_svg_escape(PROBE_REGION_OUTPUT_SVG_EXPORT_SCHEMA)
        );
        svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n");
        let _ = writeln!(
            svg,
            "<text x=\"36\" y=\"42\" font-family=\"Inter,Arial,sans-serif\" font-size=\"24\" font-weight=\"700\" fill=\"#111827\">Probe-region intensity evidence</text>"
        );
        let _ = writeln!(
            svg,
            "<text x=\"36\" y=\"70\" font-family=\"Inter,Arial,sans-serif\" font-size=\"15\" fill=\"#374151\">{}</text>",
            Self::probe_region_svg_escape(&target_label)
        );
        let _ = writeln!(
            svg,
            "<text x=\"36\" y=\"94\" font-family=\"Inter,Arial,sans-serif\" font-size=\"12\" fill=\"#6b7280\">{}</text>",
            Self::probe_region_svg_escape(&provenance)
        );

        Self::render_probe_region_svg_panel(
            &mut svg,
            "Mean log2 intensity by condition",
            &plot.intensity_tracks,
            &plot.rows,
            intensity_range,
            use_genomic_x,
            min_start,
            max_start,
            left,
            top_panel_top,
            plot_width,
            panel_height,
            &palette,
        );
        Self::render_probe_region_svg_panel(
            &mut svg,
            "Log2 fold-change tracks",
            &plot.logfc_tracks,
            &plot.rows,
            logfc_range,
            use_genomic_x,
            min_start,
            max_start,
            left,
            bottom_panel_top,
            plot_width,
            panel_height,
            &palette,
        );

        let axis_label = if use_genomic_x {
            let chromosome = plot
                .rows
                .first()
                .map(|row| row.chromosome.as_str())
                .unwrap_or("");
            format!("{chromosome}:{min_start}-{max_start}")
        } else {
            "chromosome-ordered helper rows".to_string()
        };
        let _ = writeln!(
            svg,
            "<text x=\"{:.1}\" y=\"650\" text-anchor=\"middle\" font-family=\"Inter,Arial,sans-serif\" font-size=\"13\" fill=\"#374151\">x-axis: {}</text>",
            left + plot_width / 2.0,
            Self::probe_region_svg_escape(&axis_label)
        );
        if let (Some(first), Some(last)) = (plot.rows.first(), plot.rows.last()) {
            let _ = writeln!(
                svg,
                "<text x=\"{left:.1}\" y=\"632\" font-family=\"Inter,Arial,sans-serif\" font-size=\"11\" fill=\"#6b7280\">{} {} {}</text>",
                Self::probe_region_svg_escape(&first.chromosome),
                first.start_1based,
                Self::probe_region_svg_escape(&first.gene_symbol)
            );
            let _ = writeln!(
                svg,
                "<text x=\"{:.1}\" y=\"632\" text-anchor=\"end\" font-family=\"Inter,Arial,sans-serif\" font-size=\"11\" fill=\"#6b7280\">{} {} {}</text>",
                left + plot_width,
                Self::probe_region_svg_escape(&last.chromosome),
                last.stop_1based.unwrap_or(last.start_1based),
                Self::probe_region_svg_escape(&last.gene_symbol)
            );
        }
        let _ = writeln!(
            svg,
            "<text x=\"36\" y=\"670\" font-family=\"Inter,Arial,sans-serif\" font-size=\"11\" fill=\"#6b7280\">rows={} intensity_tracks={} logfc_tracks={} projection_ready={}</text>",
            plot.rows.len(),
            plot.intensity_tracks.len(),
            plot.logfc_tracks.len(),
            inspection.projection_ready
        );
        svg.push_str("</svg>\n");
        svg
    }

    #[allow(clippy::too_many_arguments)]
    fn render_probe_region_svg_panel(
        svg: &mut String,
        title: &str,
        tracks: &[ProbeRegionPlotTrack],
        rows: &[ProbeRegionPlotRow],
        range: Option<(f64, f64)>,
        use_genomic_x: bool,
        min_start: usize,
        max_start: usize,
        left: f64,
        top: f64,
        width: f64,
        height: f64,
        palette: &[&str],
    ) {
        use std::fmt::Write as _;

        let _ = writeln!(
            svg,
            "<text x=\"36\" y=\"{:.1}\" font-family=\"Inter,Arial,sans-serif\" font-size=\"16\" font-weight=\"700\" fill=\"#111827\">{}</text>",
            top - 18.0,
            Self::probe_region_svg_escape(title)
        );
        let _ = writeln!(
            svg,
            "<rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{width:.1}\" height=\"{height:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
        );
        for frac in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let y = top + height * frac;
            let _ = writeln!(
                svg,
                "<line x1=\"{left:.1}\" y1=\"{y:.1}\" x2=\"{:.1}\" y2=\"{y:.1}\" stroke=\"#e5e7eb\"/>",
                left + width
            );
        }
        let Some((min_value, max_value)) = range else {
            let _ = writeln!(
                svg,
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"Inter,Arial,sans-serif\" font-size=\"13\" fill=\"#6b7280\">No numeric track columns available for this panel</text>",
                left + width / 2.0,
                top + height / 2.0
            );
            return;
        };
        for (value, frac) in [
            (max_value, 0.0),
            ((min_value + max_value) / 2.0, 0.5),
            (min_value, 1.0),
        ] {
            let y = top + height * frac;
            let _ = writeln!(
                svg,
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"Inter,Arial,sans-serif\" font-size=\"11\" fill=\"#6b7280\">{value:.2}</text>",
                left - 8.0,
                y + 4.0
            );
        }
        for (track_idx, track) in tracks.iter().enumerate() {
            let color = palette[track_idx % palette.len()];
            let mut points = String::new();
            for (row_idx, row) in rows.iter().enumerate() {
                let Some(value) = track.values.get(row_idx).and_then(|value| *value) else {
                    continue;
                };
                let x = Self::probe_region_plot_x(
                    row_idx,
                    row,
                    rows.len(),
                    use_genomic_x,
                    min_start,
                    max_start,
                    left,
                    width,
                );
                let y = Self::probe_region_plot_y(value, min_value, max_value, top, height);
                let _ = write!(points, "{x:.1},{y:.1} ");
            }
            if !points.trim().is_empty() {
                let _ = writeln!(
                    svg,
                    "<polyline points=\"{}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"2\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>",
                    points.trim()
                );
            }
            for (row_idx, row) in rows.iter().enumerate() {
                let Some(value) = track.values.get(row_idx).and_then(|value| *value) else {
                    continue;
                };
                let x = Self::probe_region_plot_x(
                    row_idx,
                    row,
                    rows.len(),
                    use_genomic_x,
                    min_start,
                    max_start,
                    left,
                    width,
                );
                let y = Self::probe_region_plot_y(value, min_value, max_value, top, height);
                let tooltip = format!(
                    "{} {}:{} {} {}={value:.3}",
                    row.probeset_or_region_id,
                    row.chromosome,
                    row.start_1based,
                    row.gene_symbol,
                    track.label
                );
                let _ = writeln!(
                    svg,
                    "<circle cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"2.6\" fill=\"{color}\"><title>{}</title></circle>",
                    Self::probe_region_svg_escape(&tooltip)
                );
            }
            let legend_x = left + 12.0 + (track_idx % 4) as f64 * 230.0;
            let legend_y = top + height + 28.0 + (track_idx / 4) as f64 * 18.0;
            let _ = writeln!(
                svg,
                "<line x1=\"{legend_x:.1}\" y1=\"{legend_y:.1}\" x2=\"{:.1}\" y2=\"{legend_y:.1}\" stroke=\"{color}\" stroke-width=\"3\"/>",
                legend_x + 22.0
            );
            let _ = writeln!(
                svg,
                "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"Inter,Arial,sans-serif\" font-size=\"11\" fill=\"#374151\">{}</text>",
                legend_x + 30.0,
                legend_y + 4.0,
                Self::probe_region_svg_escape(&track.label)
            );
        }
    }

    fn probe_region_region_sample_columns(
        headers: &[String],
        sample_ids_from_table: &[String],
    ) -> Vec<String> {
        if !sample_ids_from_table.is_empty() {
            let sample_keys = sample_ids_from_table
                .iter()
                .map(|sample_id| sample_id.to_ascii_lowercase())
                .collect::<BTreeSet<_>>();
            return headers
                .iter()
                .filter(|header| sample_keys.contains(&header.to_ascii_lowercase()))
                .cloned()
                .collect();
        }

        let fixed_keys = [
            "chromosome",
            "chrom",
            "start",
            "stop",
            "end",
            "strand",
            "probeset_or_region_id",
            "feature_id",
            "probeset_id",
            "psr_id",
            "transcript_cluster_id",
            "transcript_cluster",
            "exon_id",
            "number_of_probes",
            "gene_symbol",
            "gene",
        ]
        .iter()
        .map(|value| Self::probe_region_header_key(value))
        .collect::<BTreeSet<_>>();
        headers
            .iter()
            .filter(|header| {
                let key = Self::probe_region_header_key(header);
                !fixed_keys.contains(&key)
                    && !header.starts_with("mean_log2_")
                    && !header.starts_with("sd_log2_")
                    && !header.starts_with("log2FC_")
            })
            .cloned()
            .collect()
    }

    fn probe_region_preview_values(values: BTreeSet<String>) -> Vec<String> {
        values.into_iter().take(64).collect()
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
        let vendor_support_files = Self::probe_region_vendor_support_files(platform);
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
                vendor_support_files,
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
                let vendor_hint = Self::probe_region_vendor_support_detail(&vendor_support_files);
                let package_detail = if package_usable {
                    format!("Using detected Bioconductor platform package '{package}'")
                } else {
                    format!(
                        "Bioconductor platform package '{package}' was not confirmed; provide --annotation-library or install the package before execution"
                    )
                };
                if let Some(vendor_hint) = vendor_hint {
                    format!("{package_detail}. {vendor_hint}")
                } else {
                    package_detail
                }
            }),
            vendor_support_files,
            required_r_package,
        }
    }

    fn probe_region_vendor_support_files(
        platform: &ProbeRegionPlatformPlan,
    ) -> Vec<ProbeRegionFileStatus> {
        if platform.normalized != "Clariom_D_Human" {
            return Vec::new();
        }
        [
            (
                CLARIOM_D_HUMAN_PROBESET_ZIP,
                "thermofisher_clariom_d_human_hg38_probeset_zip",
            ),
            (
                CLARIOM_D_HUMAN_TRANSCRIPT_ZIP,
                "thermofisher_clariom_d_human_hg38_transcript_zip",
            ),
        ]
        .iter()
        .map(|(file_name, role)| {
            let path = Path::new(CLARIOM_D_HUMAN_VENDOR_SUPPORT_DIR).join(file_name);
            let mut status = Self::probe_region_file_status(&path.to_string_lossy(), role);
            if !status.exists {
                status.detail = Some(
                    "Manual Thermo Fisher login download required; GENtle never auto-downloads this support file".to_string(),
                );
            }
            status
        })
        .collect()
    }

    fn probe_region_vendor_support_detail(
        vendor_support_files: &[ProbeRegionFileStatus],
    ) -> Option<String> {
        if vendor_support_files.is_empty() {
            return None;
        }
        let present = vendor_support_files
            .iter()
            .filter(|file| file.exists)
            .count();
        let expected = vendor_support_files.len();
        if present == expected {
            Some(format!(
                "Thermo Fisher Clariom D hg38 support ZIPs are present under {CLARIOM_D_HUMAN_VENDOR_SUPPORT_DIR}"
            ))
        } else {
            Some(format!(
                "Thermo Fisher Clariom D hg38 support ZIPs are expected under {CLARIOM_D_HUMAN_VENDOR_SUPPORT_DIR} ({present}/{expected} present); these login-walled files must be placed manually"
            ))
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
