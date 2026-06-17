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
const PROBE_REGION_PROBE_TABLE_FILE: &str = "probe_intensity_chrom_order.csv";
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

#[derive(Default)]
struct ProbeRegionProbeTableSummary {
    row_count: usize,
    parent_feature_count: usize,
    required_columns_missing: Vec<String>,
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

#[derive(Clone)]
struct ProbeRegionProjectionRow {
    chromosome: String,
    start_1based: usize,
    end_1based: usize,
    strand: Option<char>,
    feature_id: String,
    parent_feature_id: Option<String>,
    intensity_source: Option<String>,
    transcript_cluster_id: Option<String>,
    gene_symbol: Option<String>,
    logfc_values: BTreeMap<String, f64>,
}

#[derive(Clone)]
struct ProbeRegionProjectedEvidence {
    evidence_id: String,
    level: String,
    feature_id: String,
    parent_feature_id: Option<String>,
    intensity_source: Option<String>,
    chromosome: Option<String>,
    start_1based: Option<usize>,
    end_1based: Option<usize>,
    strand: Option<String>,
    logfc: Option<f64>,
    assembly_check: Option<String>,
    ranges_0based: Vec<(usize, usize)>,
}

struct ProbeRegionTranscriptModel {
    transcript_id: String,
    gene: Option<String>,
    label: Option<String>,
    strand: Option<String>,
    exons: Vec<ProbeRegionTranscriptExon>,
    exon_ranges_0based: Vec<(usize, usize)>,
    span_0based: Option<(usize, usize)>,
}

#[derive(Clone)]
struct ProbeRegionTranscriptExon {
    ordinal: usize,
    range_0based: (usize, usize),
}

#[derive(Default)]
struct ProbeRegionTranscriptEvidenceCounts {
    compatible: usize,
    constraining: usize,
    shared: usize,
    unique: usize,
    compatible_score: f64,
    constraining_score: f64,
    shared_score: f64,
    unique_score: f64,
}

#[derive(Default)]
struct ProbeRegionAptLibraryPlan {
    pgf_path: Option<String>,
    clf_path: Option<String>,
    mps_path: Option<String>,
    source_detail: Option<String>,
}

struct ProbeRegionAptAnnotationRow {
    chromosome: String,
    start_1based: usize,
    end_1based: usize,
    strand: String,
    feature_id: String,
    transcript_cluster_id: String,
    number_of_probes: String,
    gene_symbol: String,
}

#[derive(Default)]
struct ProbeRegionAptAnnotationTable {
    regions: HashMap<String, ProbeRegionAptAnnotationRow>,
    probes: Vec<ProbeRegionAptProbeAnnotationRow>,
}

struct ProbeRegionAptProbeAnnotationRow {
    chromosome: String,
    start_1based: usize,
    end_1based: usize,
    strand: String,
    parent_feature_id: String,
    transcript_cluster_id: String,
    gene_symbol: String,
    probe_id: String,
    x: String,
    y: String,
}

#[derive(Default)]
struct ProbeRegionAptMetadataSummary {
    path: String,
    sample_column: String,
    condition_column: String,
    conditions: Vec<ProbeRegionAptConditionGroup>,
    contrasts: Vec<ProbeRegionAptContrast>,
    warnings: Vec<String>,
}

struct ProbeRegionAptConditionGroup {
    label: String,
    column_label: String,
    sample_indices: Vec<usize>,
}

struct ProbeRegionAptContrast {
    label: String,
    numerator_index: usize,
    denominator_index: usize,
}

struct ProbeRegionAptProbeIntensityTable {
    path: String,
    probe_id_column: String,
    sample_columns: Vec<String>,
    value_header: Vec<String>,
    rows: BTreeMap<String, Vec<String>>,
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
        let provenance = json!({
            "schema": PROBE_REGION_BACKEND_PROVENANCE_SCHEMA,
            "backend": "affymetrix_power_tools",
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

    pub(super) fn project_probe_region_output(
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

    pub(super) fn interpret_probe_region_evidence(
        dna: &DNAsequence,
        seq_id: &str,
        gene_label: Option<&str>,
        level: Option<&str>,
        min_abs_logfc: Option<f64>,
    ) -> Result<ProbeRegionEvidenceInterpretationReport, EngineError> {
        if let Some(threshold) = min_abs_logfc
            && (!threshold.is_finite() || threshold < 0.0)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "InterpretProbeRegionEvidence min_abs_logfc must be >= 0".to_string(),
                cause_chain: vec![],
            });
        }
        let level = Self::probe_region_interpretation_level(level)?;
        let gene_label = gene_label
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string);
        let mut warnings = Vec::new();
        let evidence = Self::projected_probe_region_evidence(dna, &level, min_abs_logfc);
        let transcripts = Self::probe_region_transcript_models(dna, gene_label.as_deref());

        if evidence.is_empty() {
            warnings.push(
                "No projected probe-region array features matched the requested interpretation filter"
                    .to_string(),
            );
        }
        if transcripts.is_empty() {
            warnings.push(
                "No transcript/exon models matched the requested gene filter; evidence is reported without transcript compatibility calls"
                    .to_string(),
            );
        }

        let mut transcript_counts = transcripts
            .iter()
            .map(|tx| {
                (
                    tx.transcript_id.clone(),
                    ProbeRegionTranscriptEvidenceCounts::default(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let mut rows = Vec::new();
        for item in &evidence {
            let mut overlapping_transcripts = Vec::new();
            let mut overlapping_exon_count = 0usize;
            let mut span_overlaps = Vec::new();
            let mut transcript_mappings = Vec::new();
            for tx in &transcripts {
                if let Some(mapping) = Self::probe_region_evidence_transcript_mapping(item, tx) {
                    if !mapping.exon_ordinals.is_empty() {
                        overlapping_exon_count += mapping.exon_ordinals.len();
                        overlapping_transcripts.push(tx.transcript_id.clone());
                    } else {
                        span_overlaps.push(tx.transcript_id.clone());
                    }
                    transcript_mappings.push(mapping);
                }
            }
            overlapping_transcripts.sort();
            overlapping_transcripts.dedup();
            span_overlaps.sort();
            span_overlaps.dedup();
            transcript_mappings.sort_by(|left, right| {
                left.transcript_id
                    .cmp(&right.transcript_id)
                    .then(left.mapping_kind.cmp(&right.mapping_kind))
            });

            let mapping_status = if transcripts.is_empty() {
                "no_transcript_models"
            } else if overlapping_transcripts.len() > 1 {
                "shared_exon_overlap"
            } else if overlapping_transcripts.len() == 1 && overlapping_exon_count > 1 {
                "unique_multi_exon_overlap"
            } else if overlapping_transcripts.len() == 1 {
                "unique_exon_overlap"
            } else if !span_overlaps.is_empty() {
                "no_exon_overlap_inside_transcript_span"
            } else {
                "no_transcript_overlap"
            }
            .to_string();
            let relationship = if transcripts.is_empty() {
                "no_transcript_models"
            } else if !overlapping_transcripts.is_empty() {
                "compatible_with_exon_geometry"
            } else if !span_overlaps.is_empty() {
                "constrains_by_non_exonic_overlap"
            } else {
                "unmapped_to_transcript_models"
            }
            .to_string();

            let mut ambiguity_tags = BTreeSet::from([
                "multi_hit_not_assessed".to_string(),
                "probe_sequence_alignment_not_assessed".to_string(),
                "isoform_support_not_inferred".to_string(),
            ]);
            if overlapping_transcripts.len() > 1 {
                ambiguity_tags.insert("shared_transcript_overlap".to_string());
            }
            if item.parent_feature_id.is_some() {
                ambiguity_tags.insert("parent_probeset_context".to_string());
            }
            if item.intensity_source.as_deref() == Some("probe_level_input") {
                ambiguity_tags.insert("pm_probe_input".to_string());
            }
            if item.assembly_check.as_deref() == Some("projected_from_native_coordinate_system") {
                ambiguity_tags.insert("coordinate_projection_used".to_string());
            }
            if transcripts.is_empty() {
                ambiguity_tags.insert("no_transcript_models".to_string());
            }
            if !span_overlaps.is_empty() && overlapping_transcripts.is_empty() {
                ambiguity_tags.insert("non_exonic_transcript_span_overlap".to_string());
            }

            for mapping in &transcript_mappings {
                if let Some(counts) = transcript_counts.get_mut(&mapping.transcript_id) {
                    if mapping.exon_ordinals.is_empty() {
                        counts.constraining += 1;
                        counts.constraining_score += mapping.geometry_score;
                    } else {
                        counts.compatible += 1;
                        counts.compatible_score += mapping.geometry_score;
                        if overlapping_transcripts.len() > 1 {
                            counts.shared += 1;
                            counts.shared_score += mapping.geometry_score;
                        } else {
                            counts.unique += 1;
                            counts.unique_score += mapping.geometry_score;
                        }
                    }
                }
            }

            rows.push(ProbeRegionEvidenceMappingRow {
                evidence_id: item.evidence_id.clone(),
                level: item.level.clone(),
                feature_id: item.feature_id.clone(),
                parent_feature_id: item.parent_feature_id.clone(),
                intensity_source: item.intensity_source.clone(),
                chromosome: item.chromosome.clone(),
                start_1based: item.start_1based,
                end_1based: item.end_1based,
                strand: item.strand.clone(),
                logfc: item.logfc,
                overlapping_transcript_ids: overlapping_transcripts,
                overlapping_exon_count,
                transcript_mappings,
                mapping_status,
                ambiguity_tags: ambiguity_tags.into_iter().collect(),
                relationship,
            });
        }

        let total_evidence = evidence.len();
        let transcript_rows = transcripts
            .iter()
            .map(|tx| {
                let counts = transcript_counts
                    .remove(&tx.transcript_id)
                    .unwrap_or_default();
                let unmapped =
                    total_evidence.saturating_sub(counts.compatible + counts.constraining);
                let relationship_summary = if counts.unique > 0 {
                    "has_unique_compatible_evidence"
                } else if counts.shared > 0 {
                    "only_shared_compatible_evidence"
                } else if counts.constraining > 0 {
                    "has_non_exonic_constraining_evidence"
                } else {
                    "no_mapped_evidence"
                }
                .to_string();
                let review_status = if counts.unique > 0 {
                    "unique_geometry_for_review"
                } else if counts.shared > 0 {
                    "shared_geometry_for_review"
                } else if counts.constraining > 0 {
                    "non_exonic_constraint_for_review"
                } else {
                    "no_probe_region_evidence"
                }
                .to_string();
                ProbeRegionEvidenceTranscriptRow {
                    transcript_id: tx.transcript_id.clone(),
                    gene: tx.gene.clone(),
                    label: tx.label.clone(),
                    strand: tx.strand.clone(),
                    exon_count: tx.exon_ranges_0based.len(),
                    compatible_evidence_count: counts.compatible,
                    constraining_evidence_count: counts.constraining,
                    shared_evidence_count: counts.shared,
                    unique_evidence_count: counts.unique,
                    unmapped_evidence_count: unmapped,
                    compatible_geometry_score: counts.compatible_score,
                    shared_geometry_score: counts.shared_score,
                    unique_geometry_score: counts.unique_score,
                    constraining_geometry_score: counts.constraining_score,
                    review_status,
                    relationship_summary,
                }
            })
            .collect::<Vec<_>>();

        Ok(ProbeRegionEvidenceInterpretationReport {
            schema: PROBE_REGION_EVIDENCE_INTERPRETATION_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            gene_label,
            level,
            min_abs_logfc,
            array_feature_count: evidence.len(),
            transcript_count: transcripts.len(),
            evidence_rows: rows,
            transcript_rows,
            warnings,
        })
    }

    fn probe_region_interpretation_level(level: Option<&str>) -> Result<String, EngineError> {
        let normalized = level
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("all")
            .replace('-', "_")
            .to_ascii_lowercase();
        match normalized.as_str() {
            "all" | "*" => Ok("all".to_string()),
            "probe_region" | "region" | "probeset" | "probeset_region" | "psr" => {
                Ok("probe_region".to_string())
            }
            "pm_probe" | "probe" | "probe_level" | "pm" => Ok("pm_probe".to_string()),
            other => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "InterpretProbeRegionEvidence level '{other}' is not supported; use all, probe_region, or pm_probe"
                ),
                cause_chain: vec![],
            }),
        }
    }

    fn projected_probe_region_evidence(
        dna: &DNAsequence,
        level: &str,
        min_abs_logfc: Option<f64>,
    ) -> Vec<ProbeRegionProjectedEvidence> {
        let mut out = Vec::new();
        for (idx, feature) in dna.features().iter().enumerate() {
            if Self::probe_region_first_qualifier(feature, &["gentle_track_source"]).as_deref()
                != Some("Array")
                || Self::probe_region_first_qualifier(feature, &["gentle_array_dataset"]).as_deref()
                    != Some("probe_region_output")
            {
                continue;
            }
            let feature_level =
                Self::probe_region_first_qualifier(feature, &["gentle_array_level"])
                    .unwrap_or_else(|| "probe_region".to_string());
            if level != "all" && feature_level != level {
                continue;
            }
            let logfc = Self::probe_region_first_qualifier(feature, &["logFC", "score"]).and_then(
                |value| {
                    value
                        .trim()
                        .parse::<f64>()
                        .ok()
                        .filter(|value| value.is_finite())
                },
            );
            if let Some(threshold) = min_abs_logfc
                && logfc.map(|value| value.abs() < threshold).unwrap_or(true)
            {
                continue;
            }
            let mut ranges = Vec::new();
            collect_location_ranges_usize(&feature.location, &mut ranges);
            if ranges.is_empty() {
                continue;
            }
            let feature_id = Self::probe_region_first_qualifier(
                feature,
                &["gentle_array_feature_id", "feature_id", "label"],
            )
            .unwrap_or_else(|| format!("array_feature_{idx}"));
            let contrast = Self::probe_region_first_qualifier(feature, &["gentle_array_contrast"]);
            let evidence_id = contrast
                .as_deref()
                .map(|contrast| format!("{feature_id}:{contrast}"))
                .unwrap_or_else(|| feature_id.clone());
            let (local_start, local_end) = Self::probe_region_range_span(&ranges).unwrap_or((0, 0));
            out.push(ProbeRegionProjectedEvidence {
                evidence_id,
                level: feature_level,
                feature_id,
                parent_feature_id: Self::probe_region_first_qualifier(
                    feature,
                    &["gentle_array_parent_feature_id"],
                ),
                intensity_source: Self::probe_region_first_qualifier(
                    feature,
                    &["gentle_array_intensity_source"],
                ),
                chromosome: Self::probe_region_first_qualifier(feature, &["chromosome"]),
                start_1based: Self::probe_region_first_qualifier(
                    feature,
                    &["genomic_start_1based", "start_1based"],
                )
                .and_then(|value| value.parse::<usize>().ok())
                .or_else(|| Some(local_start + 1)),
                end_1based: Self::probe_region_first_qualifier(
                    feature,
                    &["genomic_end_1based", "end_1based"],
                )
                .and_then(|value| value.parse::<usize>().ok())
                .or_else(|| Some(local_end)),
                strand: Self::probe_region_first_qualifier(feature, &["strand", "array_strand"]),
                logfc,
                assembly_check: Self::probe_region_first_qualifier(
                    feature,
                    &["gentle_array_assembly_check"],
                ),
                ranges_0based: ranges,
            });
        }
        out
    }

    fn probe_region_transcript_models(
        dna: &DNAsequence,
        gene_label: Option<&str>,
    ) -> Vec<ProbeRegionTranscriptModel> {
        let mut exon_ranges_by_transcript: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
        for feature in dna.features() {
            if !feature.kind.eq_ignore_ascii_case("exon") {
                continue;
            }
            let mut ranges = Vec::new();
            collect_location_ranges_usize(&feature.location, &mut ranges);
            if ranges.is_empty() {
                continue;
            }
            for transcript_id in Self::probe_region_feature_text_values(
                feature,
                &["transcript_id", "Parent", "parent", "transcript"],
            ) {
                exon_ranges_by_transcript
                    .entry(transcript_id)
                    .or_default()
                    .extend(ranges.iter().copied());
            }
        }

        let mut out = Vec::new();
        let mut seen = BTreeSet::new();
        for (idx, feature) in dna.features().iter().enumerate() {
            if !Self::probe_region_is_transcript_feature(feature) {
                continue;
            }
            if !Self::probe_region_feature_matches_gene(feature, gene_label) {
                continue;
            }
            let transcript_id = Self::probe_region_first_qualifier(
                feature,
                &[
                    "transcript_id",
                    "transcript",
                    "transcript_cluster_id",
                    "ID",
                    "id",
                    "Name",
                    "label",
                ],
            )
            .unwrap_or_else(|| format!("transcript_{idx}"));
            if !seen.insert(transcript_id.clone()) {
                continue;
            }
            let mut ranges = Vec::new();
            collect_location_ranges_usize(&feature.location, &mut ranges);
            if ranges.len() <= 1
                && let Some(exon_ranges) = exon_ranges_by_transcript.get(&transcript_id)
            {
                ranges = exon_ranges.clone();
            }
            if ranges.is_empty() {
                continue;
            }
            ranges.sort();
            ranges.dedup();
            let strand = Self::probe_region_first_qualifier(feature, &["strand"]);
            let is_reverse = strand.as_deref() == Some("-");
            let exon_count = ranges.len();
            let exons = ranges
                .iter()
                .enumerate()
                .map(|(idx, range)| ProbeRegionTranscriptExon {
                    ordinal: if is_reverse {
                        exon_count.saturating_sub(idx)
                    } else {
                        idx.saturating_add(1)
                    },
                    range_0based: *range,
                })
                .collect::<Vec<_>>();
            out.push(ProbeRegionTranscriptModel {
                transcript_id,
                gene: Self::probe_region_first_qualifier(
                    feature,
                    &["gene", "gene_symbol", "gene_name", "gene_id"],
                ),
                label: Self::probe_region_first_qualifier(feature, &["label", "Name"]),
                strand,
                exons,
                span_0based: Self::probe_region_range_span(&ranges),
                exon_ranges_0based: ranges,
            });
        }
        out
    }

    fn probe_region_is_transcript_feature(feature: &gb_io::seq::Feature) -> bool {
        matches!(
            feature.kind.to_ascii_lowercase().as_str(),
            "mrna" | "transcript" | "ncrna" | "misc_rna" | "rrna" | "trna"
        )
    }

    fn probe_region_feature_matches_gene(
        feature: &gb_io::seq::Feature,
        gene_label: Option<&str>,
    ) -> bool {
        let Some(requested) = gene_label.map(str::trim).filter(|value| !value.is_empty()) else {
            return true;
        };
        let requested_lower = requested.to_ascii_lowercase();
        let values = Self::probe_region_feature_text_values(
            feature,
            &[
                "gene",
                "gene_symbol",
                "gene_name",
                "gene_id",
                "label",
                "Name",
                "transcript_id",
                "transcript",
            ],
        );
        values.iter().any(|value| {
            let lower = value.to_ascii_lowercase();
            lower == requested_lower || lower.contains(&requested_lower)
        })
    }

    fn probe_region_feature_text_values(
        feature: &gb_io::seq::Feature,
        keys: &[&str],
    ) -> Vec<String> {
        let mut out = Vec::new();
        let mut seen = BTreeSet::new();
        for key in keys {
            for value in feature.qualifier_values(key) {
                for part in value
                    .split([',', ';'])
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    if seen.insert(part.to_ascii_lowercase()) {
                        out.push(part.to_string());
                    }
                }
            }
        }
        out
    }

    fn probe_region_first_qualifier(
        feature: &gb_io::seq::Feature,
        keys: &[&str],
    ) -> Option<String> {
        keys.iter().find_map(|key| {
            feature
                .qualifier_values(key)
                .find(|value| !value.trim().is_empty())
                .map(|value| value.trim().to_string())
        })
    }

    fn probe_region_ranges_overlap(left: &(usize, usize), right: &(usize, usize)) -> bool {
        left.0 < right.1 && right.0 < left.1
    }

    fn probe_region_range_overlap_bp(left: &(usize, usize), right: &(usize, usize)) -> usize {
        let start = left.0.max(right.0);
        let end = left.1.min(right.1);
        end.saturating_sub(start)
    }

    fn probe_region_range_span(ranges: &[(usize, usize)]) -> Option<(usize, usize)> {
        let start = ranges.iter().map(|(start, _)| *start).min()?;
        let end = ranges.iter().map(|(_, end)| *end).max()?;
        Some((start, end))
    }

    fn probe_region_evidence_transcript_mapping(
        item: &ProbeRegionProjectedEvidence,
        tx: &ProbeRegionTranscriptModel,
    ) -> Option<ProbeRegionEvidenceTranscriptMapping> {
        let mut exon_hits = tx
            .exons
            .iter()
            .filter_map(|exon| {
                let overlap_bp = item
                    .ranges_0based
                    .iter()
                    .map(|evidence_range| {
                        Self::probe_region_range_overlap_bp(&exon.range_0based, evidence_range)
                    })
                    .sum::<usize>();
                (overlap_bp > 0).then_some((exon, overlap_bp))
            })
            .collect::<Vec<_>>();
        exon_hits.sort_by(|left, right| left.0.ordinal.cmp(&right.0.ordinal));

        if !exon_hits.is_empty() {
            let exon_ordinals = exon_hits
                .iter()
                .map(|(exon, _)| exon.ordinal)
                .collect::<Vec<_>>();
            let exon_ranges_1based = exon_hits
                .iter()
                .map(|(exon, _)| {
                    format!(
                        "{}..{}",
                        exon.range_0based.0.saturating_add(1),
                        exon.range_0based.1
                    )
                })
                .collect::<Vec<_>>();
            let overlap_bp = exon_hits.iter().map(|(_, overlap)| *overlap).sum();
            let junction_spans = Self::probe_region_evidence_junction_spans(item, tx, &exon_hits);
            let mapping_kind = if !junction_spans.is_empty() {
                "junction_spanning_exon_overlap"
            } else if exon_ordinals.len() > 1 {
                "multi_exon_overlap"
            } else {
                "exon_overlap"
            }
            .to_string();
            let (geometry_score, geometry_score_class, score_basis) =
                Self::probe_region_evidence_geometry_score(
                    &mapping_kind,
                    overlap_bp,
                    junction_spans.len(),
                );
            return Some(ProbeRegionEvidenceTranscriptMapping {
                transcript_id: tx.transcript_id.clone(),
                mapping_kind,
                geometry_score,
                geometry_score_class,
                score_basis,
                exon_ordinals,
                exon_ranges_1based,
                junction_spans,
                overlap_bp,
            });
        }

        let span = tx.span_0based?;
        let overlap_bp = item
            .ranges_0based
            .iter()
            .map(|range| Self::probe_region_range_overlap_bp(&span, range))
            .sum::<usize>();
        (overlap_bp > 0).then(|| ProbeRegionEvidenceTranscriptMapping {
            transcript_id: tx.transcript_id.clone(),
            mapping_kind: "transcript_span_non_exonic".to_string(),
            geometry_score: -0.25,
            geometry_score_class: "constraining_non_exonic_geometry".to_string(),
            score_basis: vec![
                "mapping_kind=transcript_span_non_exonic".to_string(),
                format!("overlap_bp={overlap_bp}"),
                "probe_sequence_alignment_not_assessed".to_string(),
                "isoform_support_not_inferred".to_string(),
            ],
            exon_ordinals: Vec::new(),
            exon_ranges_1based: Vec::new(),
            junction_spans: Vec::new(),
            overlap_bp,
        })
    }

    fn probe_region_evidence_geometry_score(
        mapping_kind: &str,
        overlap_bp: usize,
        junction_span_count: usize,
    ) -> (f64, String, Vec<String>) {
        let (score, class) = match mapping_kind {
            "junction_spanning_exon_overlap" => (0.75, "junction_spanning_geometry"),
            "multi_exon_overlap" => (0.60, "multi_exon_geometry"),
            "exon_overlap" => (0.50, "exon_geometry"),
            _ => (0.0, "unscored_geometry"),
        };
        (
            score,
            class.to_string(),
            vec![
                format!("mapping_kind={mapping_kind}"),
                format!("overlap_bp={overlap_bp}"),
                format!("junction_spans={junction_span_count}"),
                "probe_sequence_alignment_not_assessed".to_string(),
                "isoform_support_not_inferred".to_string(),
            ],
        )
    }

    fn probe_region_evidence_junction_spans(
        item: &ProbeRegionProjectedEvidence,
        tx: &ProbeRegionTranscriptModel,
        exon_hits: &[(&ProbeRegionTranscriptExon, usize)],
    ) -> Vec<ProbeRegionEvidenceJunctionSpan> {
        let hit_ordinals = exon_hits
            .iter()
            .map(|(exon, _)| exon.ordinal)
            .collect::<BTreeSet<_>>();
        let mut out = Vec::new();
        let mut genomic_exons = tx.exons.iter().collect::<Vec<_>>();
        genomic_exons.sort_by(|left, right| {
            left.range_0based
                .0
                .cmp(&right.range_0based.0)
                .then(left.range_0based.1.cmp(&right.range_0based.1))
        });
        for pair in genomic_exons.windows(2) {
            let left = pair[0];
            let right = pair[1];
            if !hit_ordinals.contains(&left.ordinal) || !hit_ordinals.contains(&right.ordinal) {
                continue;
            }
            let gap = (left.range_0based.1, right.range_0based.0);
            let one_range_spans_boundary = item
                .ranges_0based
                .iter()
                .any(|range| range.0 < left.range_0based.1 && range.1 > right.range_0based.0);
            let joined_evidence_hits_both_sides = item.ranges_0based.len() > 1
                && item
                    .ranges_0based
                    .iter()
                    .any(|range| Self::probe_region_ranges_overlap(&left.range_0based, range))
                && item
                    .ranges_0based
                    .iter()
                    .any(|range| Self::probe_region_ranges_overlap(&right.range_0based, range));
            if !one_range_spans_boundary && !joined_evidence_hits_both_sides {
                continue;
            }
            let is_reverse = tx.strand.as_deref() == Some("-");
            let (from_exon_ordinal, to_exon_ordinal) = if is_reverse {
                (right.ordinal, left.ordinal)
            } else {
                (left.ordinal, right.ordinal)
            };
            out.push(ProbeRegionEvidenceJunctionSpan {
                from_exon_ordinal,
                to_exon_ordinal,
                genomic_start_1based: gap.0.saturating_add(1),
                genomic_end_1based: gap.1,
            });
        }
        out.sort_by(|left, right| {
            left.from_exon_ordinal
                .cmp(&right.from_exon_ordinal)
                .then(left.to_exon_ordinal.cmp(&right.to_exon_ordinal))
        });
        out
    }

    fn probe_region_output_supports_anchor(
        coordinate_system: &str,
        genome_build: &str,
        anchor: &GenomeSequenceAnchor,
    ) -> bool {
        let anchor_id = anchor.genome_id.trim();
        !anchor_id.is_empty()
            && (coordinate_system.trim().eq_ignore_ascii_case(anchor_id)
                || genome_build.trim().eq_ignore_ascii_case(anchor_id))
    }

    fn resolve_probe_region_coordinate_projection_spec<'a>(
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

    fn resolve_probe_region_output_path(output_dir: &str, path: &str) -> String {
        let path = Path::new(path);
        if path.is_absolute() {
            return path.to_string_lossy().to_string();
        }
        Path::new(output_dir)
            .join(path)
            .to_string_lossy()
            .to_string()
    }

    fn probe_region_apt_annotation_rows(
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

    fn probe_region_apt_probe_intensity_rows(
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

    fn probe_region_apt_metadata_summary(
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

    fn probe_region_apt_sample_lookup(
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

    fn probe_region_apt_sample_match_keys(value: &str) -> Vec<String> {
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

    fn probe_region_apt_output_label(value: &str) -> String {
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

    fn probe_region_apt_condition_stats(
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

    fn probe_region_json_alias_string(value: &serde_json::Value, keys: &[&str]) -> Option<String> {
        keys.iter()
            .find_map(|key| Self::probe_region_json_string(value, key))
    }

    fn probe_region_json_coordinate_projections(
        value: &serde_json::Value,
    ) -> Vec<GenomeCoordinateProjectionSpec> {
        let mut out = Vec::new();
        let mut seen = BTreeSet::new();
        for key in ["coordinate_projections", "projection_maps"] {
            let Some(items) = value.get(key).and_then(|items| items.as_array()) else {
                continue;
            };
            for item in items {
                let source_genome_id = Self::probe_region_json_alias_string(
                    item,
                    &["source_genome_id", "source_genome", "from_genome_id"],
                )
                .unwrap_or_default();
                let target_genome_id = Self::probe_region_json_alias_string(
                    item,
                    &["target_genome_id", "target_genome", "to_genome_id"],
                )
                .unwrap_or_default();
                let method = Self::probe_region_json_string(item, "method")
                    .unwrap_or_else(|| "interval_map".to_string());
                let path = Self::probe_region_json_alias_string(
                    item,
                    &["path", "projection_path", "map_path"],
                )
                .unwrap_or_default();
                if source_genome_id.trim().is_empty()
                    || target_genome_id.trim().is_empty()
                    || path.trim().is_empty()
                {
                    continue;
                }
                let spec = GenomeCoordinateProjectionSpec {
                    source_genome_id: source_genome_id.trim().to_string(),
                    target_genome_id: target_genome_id.trim().to_string(),
                    method: method.trim().to_string(),
                    path: path.trim().to_string(),
                };
                let dedupe_key = format!(
                    "{}\t{}\t{}",
                    spec.source_genome_id.to_ascii_lowercase(),
                    spec.target_genome_id.to_ascii_lowercase(),
                    spec.path
                );
                if seen.insert(dedupe_key) {
                    out.push(spec);
                }
            }
        }
        out
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

    fn probe_region_extend_unique_coordinate_projections(
        target: &mut Vec<GenomeCoordinateProjectionSpec>,
        values: Vec<GenomeCoordinateProjectionSpec>,
    ) {
        let mut seen = target
            .iter()
            .map(|spec| {
                format!(
                    "{}\t{}\t{}",
                    spec.source_genome_id.to_ascii_lowercase(),
                    spec.target_genome_id.to_ascii_lowercase(),
                    spec.path
                )
            })
            .collect::<BTreeSet<_>>();
        for spec in values {
            let key = format!(
                "{}\t{}\t{}",
                spec.source_genome_id.to_ascii_lowercase(),
                spec.target_genome_id.to_ascii_lowercase(),
                spec.path
            );
            if seen.insert(key) {
                target.push(spec);
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

    fn inspect_probe_region_probe_table(
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

    fn probe_region_projection_rows_from_table(
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
    fn build_probe_region_output_feature(
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

    fn probe_region_format_number(value: f64) -> String {
        if value.abs() >= 100.0 || (value != 0.0 && value.abs() < 0.001) {
            format!("{value:.3e}")
        } else {
            format!("{value:.3}")
        }
    }

    fn probe_region_plot_track_label(header: &str) -> String {
        header
            .strip_prefix("mean_log2_")
            .or_else(|| header.strip_prefix("log2FC_"))
            .unwrap_or(header)
            .replace('_', " ")
    }

    fn probe_region_projection_level(level: Option<&str>) -> Result<String, EngineError> {
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
        let apt_normalization_supported = normalization == "rma";
        let apt_library = Self::probe_region_apt_library_plan(annotation_source);

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
        if !apt_normalization_supported {
            apt_missing.push(format!("normalization=rma (requested {normalization})"));
        }
        if !apt_present {
            apt_missing.push("apt-probeset-summarize command".to_string());
        }
        if !annotation_path_usable {
            apt_missing.push("user-supplied APT annotation-library path".to_string());
        }
        if apt_library.pgf_path.is_none() {
            apt_missing.push("APT PGF library file".to_string());
        }
        if apt_library.clf_path.is_none() {
            apt_missing.push("APT CLF library file".to_string());
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
                    "normalization=rma".to_string(),
                    "APT PGF library file".to_string(),
                    "APT CLF library file".to_string(),
                    "optional APT MPS/meta-probesets file".to_string(),
                    "explicit CEL file paths".to_string(),
                ],
                missing: apt_missing,
                helper_script: None,
                suggested_command: Self::probe_region_apt_suggested_command(
                    request,
                    &apt_library,
                    normalization,
                ),
                detail: Some(Self::probe_region_apt_backend_detail(&apt_library)),
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

    fn probe_region_apt_library_plan(
        annotation_source: &ProbeRegionAnnotationSourcePlan,
    ) -> ProbeRegionAptLibraryPlan {
        let Some(path_status) = annotation_source.path.as_ref() else {
            return ProbeRegionAptLibraryPlan::default();
        };
        if !path_status.exists {
            return ProbeRegionAptLibraryPlan::default();
        }
        let path = Path::new(&path_status.path);
        let mut plan = ProbeRegionAptLibraryPlan {
            source_detail: Some(format!("APT library source: {}", path_status.path)),
            ..Default::default()
        };
        if path.is_file() {
            Self::probe_region_note_apt_library_file(&mut plan, path);
            return plan;
        }
        if path.is_dir()
            && let Ok(entries) = std::fs::read_dir(path)
        {
            let mut candidates = entries
                .filter_map(Result::ok)
                .map(|entry| entry.path())
                .filter(|entry_path| entry_path.is_file())
                .collect::<Vec<_>>();
            candidates.sort();
            for candidate in candidates {
                Self::probe_region_note_apt_library_file(&mut plan, &candidate);
            }
        }
        plan
    }

    fn probe_region_note_apt_library_file(plan: &mut ProbeRegionAptLibraryPlan, path: &Path) {
        let ext = path
            .extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        let path_text = path.to_string_lossy().to_string();
        match ext.as_str() {
            "pgf" if plan.pgf_path.is_none() => plan.pgf_path = Some(path_text),
            "clf" if plan.clf_path.is_none() => plan.clf_path = Some(path_text),
            "mps" if plan.mps_path.is_none() => plan.mps_path = Some(path_text),
            _ => {}
        }
    }

    fn probe_region_apt_backend_detail(library: &ProbeRegionAptLibraryPlan) -> String {
        let mut parts = vec![
            "Explicit Affymetrix Power Tools preflight; GENtle reports the command but does not run APT implicitly".to_string(),
        ];
        if let Some(detail) = &library.source_detail {
            parts.push(detail.clone());
        }
        if let Some(pgf) = &library.pgf_path {
            parts.push(format!("PGF={pgf}"));
        }
        if let Some(clf) = &library.clf_path {
            parts.push(format!("CLF={clf}"));
        }
        if let Some(mps) = &library.mps_path {
            parts.push(format!("MPS={mps}"));
        }
        parts.join("; ")
    }

    fn probe_region_apt_suggested_command(
        request: &ProbeRegionRequest,
        library: &ProbeRegionAptLibraryPlan,
        normalization: &str,
    ) -> Option<String> {
        if normalization != "rma" || request.cel_paths.is_empty() {
            return None;
        }
        let pgf = library.pgf_path.as_ref()?;
        let clf = library.clf_path.as_ref()?;
        let output_dir = request
            .output_dir
            .clone()
            .unwrap_or_else(|| "analysis/probe_regions/apt".to_string());
        let mut args = vec![
            "apt-probeset-summarize".to_string(),
            "-a".to_string(),
            "rma-sketch".to_string(),
            "-p".to_string(),
            pgf.clone(),
            "-c".to_string(),
            clf.clone(),
        ];
        if let Some(mps) = &library.mps_path {
            args.push("-m".to_string());
            args.push(mps.clone());
        }
        args.push("-o".to_string());
        args.push(output_dir);
        args.extend(request.cel_paths.iter().cloned());
        Some(
            args.iter()
                .map(|arg| Self::probe_region_shell_quote(arg))
                .collect::<Vec<_>>()
                .join(" "),
        )
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
