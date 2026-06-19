use super::*;

impl GentleEngine {
    pub fn plan_probe_regions(&self, request: ProbeRegionRequest) -> ProbeRegionPlan {
        let mut request = Self::normalize_probe_region_request(request);
        if request.normalization.is_empty() {
            request.normalization = "rma".to_string();
        }

        let mut warnings = Vec::new();
        let mut errors = Vec::new();
        let user_supplied_cel_paths = !request.cel_paths.is_empty();
        let input_mode = Self::probe_region_input_mode(&request);
        Self::probe_region_apply_publication_dataset_defaults(
            &mut request,
            &mut warnings,
            &mut errors,
        );
        if request.cel_paths.is_empty() && request.dataset.is_none() {
            errors.push(
                "arrays probe-regions requires at least one --cel path or a --dataset id"
                    .to_string(),
            );
        }
        if user_supplied_cel_paths && request.dataset.is_some() {
            warnings.push(
                "Both explicit CEL files and --dataset were supplied; the current stage records the explicit CEL files and uses the dataset only for resource context"
                    .to_string(),
            );
        }
        if request.dataset.is_some() && request.cel_paths.is_empty() {
            warnings.push(
                "Publication-resource dataset did not declare any CEL files for probe-region planning"
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

        let mut planned_outputs = vec![
            "region_intensity_chrom_order.csv".to_string(),
            "probe_intensity_chrom_order.csv when PM probe coordinates are available".to_string(),
            "plot.png and plot.svg when --plot is set".to_string(),
            "provenance.json".to_string(),
        ];
        if request.output_dir.is_some() {
            planned_outputs.insert(0, PROBE_REGION_PLAN_FILE.to_string());
        }
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

    /// Persist a stage-one probe-region plan beside later backend outputs.

    pub fn write_probe_region_plan(
        &self,
        plan: &ProbeRegionPlan,
        output_dir: &str,
    ) -> Result<String, EngineError> {
        let output_dir = output_dir.trim();
        if output_dir.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Probe-region plan output directory must not be empty".to_string(),
                cause_chain: vec![],
            });
        }
        let output_dir_path = Path::new(output_dir);
        std::fs::create_dir_all(output_dir_path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create probe-region plan output directory '{}': {e}",
                output_dir_path.to_string_lossy()
            ),
            cause_chain: vec![],
        })?;
        let plan_path = output_dir_path.join(PROBE_REGION_PLAN_FILE);
        std::fs::write(
            &plan_path,
            serde_json::to_string_pretty(plan).unwrap_or_else(|_| "{}".to_string()),
        )
        .map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write probe-region plan '{}': {e}",
                plan_path.to_string_lossy()
            ),
            cause_chain: vec![],
        })?;
        Ok(plan_path.to_string_lossy().to_string())
    }

    /// Run an explicitly approved external R/oligo or APT backend from plan.json.

    pub fn run_probe_region_backend(
        &self,
        plan_path: &str,
        backend: Option<&str>,
        allow_external_execution: bool,
    ) -> Result<ProbeRegionBackendRunReport, EngineError> {
        let plan_path = plan_path.trim();
        if plan_path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "arrays run-probe-region-backend requires PLAN.json".to_string(),
                cause_chain: vec![],
            });
        }
        let plan_text = std::fs::read_to_string(plan_path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read probe-region plan '{plan_path}': {e}"),
            cause_chain: vec![],
        })?;
        let plan: ProbeRegionPlan = serde_json::from_str(&plan_text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse probe-region plan '{plan_path}': {e}"),
            cause_chain: vec![],
        })?;
        if plan.schema != PROBE_REGION_PLAN_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region plan '{plan_path}' has schema '{}' but expected '{}'",
                    plan.schema, PROBE_REGION_PLAN_SCHEMA
                ),
                cause_chain: vec![],
            });
        }
        let candidate = Self::probe_region_select_backend_candidate(&plan, backend)?;
        let command = candidate
            .suggested_command
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region backend '{}' has no suggested command in plan '{}'",
                    candidate.backend, plan_path
                ),
                cause_chain: vec![],
            })?
            .to_string();
        let output_dir = Self::probe_region_backend_output_dir(&plan, candidate);
        let mut report = ProbeRegionBackendRunReport {
            schema: PROBE_REGION_BACKEND_RUN_SCHEMA.to_string(),
            plan_path: plan_path.to_string(),
            backend: candidate.backend.clone(),
            command: command.clone(),
            output_dir: Some(output_dir.clone()),
            allow_external_execution,
            ..Default::default()
        };
        if !allow_external_execution {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "arrays run-probe-region-backend refuses to run external R/APT without --allow-external-execution".to_string(),
                cause_chain: vec![],
            });
        }
        if !plan.preflight_ok {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region plan preflight is not OK; fix before running backend: {}",
                    plan.errors.join("; ")
                ),
                cause_chain: vec![],
            });
        }
        let readiness_errors = Self::probe_region_backend_readiness_errors(candidate);
        if !readiness_errors.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region backend '{}' is not ready: {}",
                    candidate.backend,
                    readiness_errors.join("; ")
                ),
                cause_chain: vec![],
            });
        }
        let argv = Self::probe_region_split_backend_command(&command)?;
        let Some((program, args)) = argv.split_first() else {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region backend '{}' rendered an empty command",
                    candidate.backend
                ),
                cause_chain: vec![],
            });
        };
        let output = Command::new(program)
            .args(args)
            .output()
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not execute probe-region backend '{}' command '{}': {e}",
                    candidate.backend, command
                ),
                cause_chain: vec![],
            })?;
        report.executed = true;
        report.exit_code = output.status.code();
        report.stdout = String::from_utf8_lossy(&output.stdout).to_string();
        report.stderr = String::from_utf8_lossy(&output.stderr).to_string();
        if !output.status.success() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region backend '{}' failed with status {:?}: {}{}",
                    candidate.backend,
                    output.status.code(),
                    report.stderr.trim(),
                    if report.stdout.trim().is_empty() {
                        ""
                    } else {
                        "; see stdout in backend run report"
                    }
                ),
                cause_chain: vec![],
            });
        }

        let inspection = self.inspect_probe_region_output(&output_dir)?;
        if !inspection.usable {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region backend '{}' completed but output '{}' failed the helper-output contract: {}",
                    candidate.backend,
                    output_dir,
                    inspection.errors.join("; ")
                ),
                cause_chain: vec![],
            });
        }
        Self::write_probe_region_backend_run_provenance(&plan, candidate, &report, &inspection)?;
        report.inspection = Some(inspection);
        Ok(report)
    }

    pub(super) fn probe_region_select_backend_candidate<'a>(
        plan: &'a ProbeRegionPlan,
        requested_backend: Option<&str>,
    ) -> Result<&'a ProbeRegionBackendCandidate, EngineError> {
        if let Some(requested) = requested_backend
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            return plan
                .backend_candidates
                .iter()
                .find(|candidate| candidate.backend.eq_ignore_ascii_case(requested))
                .ok_or_else(|| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Probe-region plan has no backend candidate '{requested}' (available: {})",
                        plan.backend_candidates
                            .iter()
                            .map(|candidate| candidate.backend.as_str())
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                    cause_chain: vec![],
                });
        }
        plan.backend_candidates
            .iter()
            .find(|candidate| {
                candidate.backend != "plan_only"
                    && candidate.status == "ready"
                    && candidate.suggested_command.is_some()
            })
            .or_else(|| {
                plan.backend_candidates.iter().find(|candidate| {
                    candidate.backend != "plan_only" && candidate.suggested_command.is_some()
                })
            })
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Probe-region plan has no executable backend candidate with a suggested command"
                        .to_string(),
                cause_chain: vec![],
            })
    }

    pub(super) fn probe_region_backend_readiness_errors(
        candidate: &ProbeRegionBackendCandidate,
    ) -> Vec<String> {
        let mut errors = Vec::new();
        if candidate.backend == "plan_only" {
            errors.push("plan_only is not an executable backend".to_string());
        }
        if candidate.status != "ready" {
            errors.push(format!("status={}", candidate.status));
        }
        errors.extend(
            candidate
                .missing
                .iter()
                .map(|missing| format!("missing {missing}")),
        );
        if candidate
            .suggested_command
            .as_deref()
            .map(str::trim)
            .unwrap_or_default()
            .is_empty()
        {
            errors.push("missing suggested_command".to_string());
        }
        errors
    }

    pub(super) fn probe_region_backend_output_dir(
        plan: &ProbeRegionPlan,
        candidate: &ProbeRegionBackendCandidate,
    ) -> String {
        plan.request.output_dir.clone().unwrap_or_else(|| {
            if candidate.backend == "affymetrix_power_tools" {
                "analysis/probe_regions/apt".to_string()
            } else {
                "analysis/probe_regions".to_string()
            }
        })
    }

    pub(super) fn probe_region_split_backend_command(
        command: &str,
    ) -> Result<Vec<String>, EngineError> {
        let mut out = Vec::new();
        let mut current = String::new();
        let mut chars = command.chars().peekable();
        let mut in_single = false;
        let mut in_double = false;
        while let Some(ch) = chars.next() {
            match ch {
                '\'' if !in_double => in_single = !in_single,
                '"' if !in_single => in_double = !in_double,
                '\\' if !in_single => {
                    if let Some(next) = chars.next() {
                        current.push(next);
                    }
                }
                ch if ch.is_whitespace() && !in_single && !in_double => {
                    if !current.is_empty() {
                        out.push(std::mem::take(&mut current));
                    }
                }
                _ => current.push(ch),
            }
        }
        if in_single || in_double {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Unterminated quote in probe-region backend command: {command}"),
                cause_chain: vec![],
            });
        }
        if !current.is_empty() {
            out.push(current);
        }
        Ok(out)
    }

    pub(super) fn write_probe_region_backend_run_provenance(
        plan: &ProbeRegionPlan,
        candidate: &ProbeRegionBackendCandidate,
        report: &ProbeRegionBackendRunReport,
        inspection: &ProbeRegionOutputInspection,
    ) -> Result<(), EngineError> {
        let output_dir = report.output_dir.as_deref().unwrap_or_default();
        let output_dir_path = Path::new(output_dir);
        let provenance_path = output_dir_path.join(PROBE_REGION_PROVENANCE_FILE);
        let mut input_fingerprints = Vec::new();
        for cel in &plan.request.cel_paths {
            input_fingerprints.push(Self::probe_region_input_fingerprint(cel, "cel"));
        }
        if let Some(metadata) = &plan.request.metadata_path {
            input_fingerprints.push(Self::probe_region_input_fingerprint(metadata, "metadata"));
        }
        if let Some(path) = plan.annotation_source.path.as_ref() {
            input_fingerprints.push(Self::probe_region_input_fingerprint(
                &path.path,
                "annotation_library",
            ));
        }
        let output_fingerprints = [
            PROBE_REGION_TABLE_FILE,
            PROBE_REGION_PROBE_TABLE_FILE,
            PROBE_REGION_MATRIX_MANIFEST_FILE,
        ]
        .iter()
        .map(|file_name| {
            Self::probe_region_input_fingerprint(
                &output_dir_path.join(file_name).to_string_lossy(),
                "backend_output",
            )
        })
        .collect::<Vec<_>>();
        let provenance = json!({
            "schema": PROBE_REGION_BACKEND_PROVENANCE_SCHEMA,
            "backend": candidate.backend,
            "declared_backend": candidate.backend,
            "selected_backend": candidate.backend,
            "backend_execution_policy": "explicit_user_invoked_by_arrays_run_probe_region_backend",
            "backend_command_source": "persisted ProbeRegionPlan backend_candidates[].suggested_command",
            "rendered_backend_command": report.command,
            "actual_backend_command": report.command,
            "exit_code": report.exit_code,
            "tool_version_checks": plan.dependencies,
            "input_fingerprints": input_fingerprints,
            "output_fingerprints": output_fingerprints,
            "output_dir": output_dir,
            "platform": plan.platform.normalized,
            "normalization": plan.request.normalization,
            "coordinate_system": inspection.coordinate_system,
            "genome_build": inspection.genome_build,
            "artifacts": inspection.artifact_paths,
            "stdout": report.stdout,
            "stderr": report.stderr,
            "warnings": report.warnings,
            "errors": report.errors
        });
        std::fs::write(
            &provenance_path,
            serde_json::to_string_pretty(&provenance).unwrap_or_else(|_| "{}".to_string()),
        )
        .map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write probe-region backend provenance '{}': {e}",
                provenance_path.to_string_lossy()
            ),
            cause_chain: vec![],
        })
    }

    pub(super) fn normalize_probe_region_request(
        mut request: ProbeRegionRequest,
    ) -> ProbeRegionRequest {
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

    pub(super) fn probe_region_has_selector(request: &ProbeRegionRequest) -> bool {
        !request.genes.is_empty()
            || !request.loci.is_empty()
            || !request.transcript_cluster_ids.is_empty()
            || !request.probeset_ids.is_empty()
    }

    pub(super) fn probe_region_input_mode(request: &ProbeRegionRequest) -> String {
        match (!request.cel_paths.is_empty(), request.dataset.is_some()) {
            (true, true) => "explicit_cel_plus_dataset".to_string(),
            (true, false) => "explicit_cel".to_string(),
            (false, true) => "publication_resource_dataset".to_string(),
            (false, false) => "missing".to_string(),
        }
    }

    pub(super) fn probe_region_apply_publication_dataset_defaults(
        request: &mut ProbeRegionRequest,
        warnings: &mut Vec<String>,
        errors: &mut Vec<String>,
    ) {
        let Some(dataset_id) = request
            .dataset
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
        else {
            return;
        };

        let status =
            match publication_resources::publication_dataset_status(&dataset_id, None, None) {
                Ok(status) => status,
                Err(err) => {
                    errors.push(format!(
                        "Publication-resource dataset '{dataset_id}' could not be resolved: {err}"
                    ));
                    return;
                }
            };

        let declared_cel_paths: Vec<String> = status
            .files
            .iter()
            .filter(|file| {
                file.category.eq_ignore_ascii_case("raw_microarray")
                    && file.file_name.to_ascii_lowercase().ends_with(".cel")
            })
            .map(|file| file.local_path.clone())
            .collect();
        let missing_cel_count = status
            .files
            .iter()
            .filter(|file| {
                file.category.eq_ignore_ascii_case("raw_microarray")
                    && file.file_name.to_ascii_lowercase().ends_with(".cel")
                    && !file.exists
            })
            .count();

        let using_declared_cel_paths = request.cel_paths.is_empty();
        if using_declared_cel_paths && !declared_cel_paths.is_empty() {
            let declared_count = declared_cel_paths.len();
            request.cel_paths = declared_cel_paths;
            warnings.push(format!(
                "Publication-resource dataset '{}' resolved {declared_count} declared CEL file(s) under '{}'",
                status.accession, status.install_dir
            ));
        }
        if using_declared_cel_paths && missing_cel_count > 0 {
            warnings.push(format!(
                "Publication-resource dataset '{}' is not locally ready: {missing_cel_count} declared CEL file(s) are missing; run `resources prepare-publication-dataset {} --categories raw_microarray --download-files` explicitly before execution",
                status.accession, status.accession
            ));
        }

        if request.metadata_path.is_none() {
            if let Some(sdrf) = status.files.iter().find(|file| {
                file.category.eq_ignore_ascii_case("metadata")
                    && file.file_name.to_ascii_lowercase().contains("sdrf")
            }) {
                if sdrf.exists {
                    request.metadata_path = Some(sdrf.local_path.clone());
                    if request.sample_column.is_none() {
                        request.sample_column = Some("Array Data File".to_string());
                    }
                    if request.condition_column.is_none() {
                        request.condition_column =
                            Some("Characteristics[genetic modification]".to_string());
                    }
                    warnings.push(format!(
                        "Using publication-resource SDRF metadata '{}' with sample column 'Array Data File' and condition column 'Characteristics[genetic modification]'",
                        sdrf.local_path
                    ));
                } else {
                    warnings.push(format!(
                        "Publication-resource dataset '{}' declares SDRF metadata '{}', but it is not present locally; metadata parsing is skipped",
                        status.accession, sdrf.local_path
                    ));
                }
            }
        }
    }

    pub(super) fn probe_region_dedupe_nonempty(values: &[String]) -> Vec<String> {
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

    pub(super) fn probe_region_trim_option(value: Option<String>) -> Option<String> {
        value
            .map(|text| text.trim().to_string())
            .filter(|text| !text.is_empty())
    }

    pub(super) fn probe_region_file_status(path: &str, role: &str) -> ProbeRegionFileStatus {
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

    pub(super) fn probe_region_input_fingerprint(path: &str, role: &str) -> serde_json::Value {
        let status = Self::probe_region_file_status(path, role);
        let sha1 = if status.is_file
            && status
                .size_bytes
                .is_some_and(|size| size <= PROBE_REGION_FINGERPRINT_SHA1_MAX_BYTES)
        {
            std::fs::read(&status.path)
                .ok()
                .map(|bytes| format!("{:x}", Sha1::digest(&bytes)))
        } else {
            None
        };
        json!({
            "path": status.path,
            "role": status.role,
            "exists": status.exists,
            "is_file": status.is_file,
            "size_bytes": status.size_bytes,
            "modified_unix_seconds": status.modified_unix_seconds,
            "sha1": sha1
        })
    }

    pub(super) fn probe_region_backend_tool_version_checks() -> Vec<ProbeRegionDependencyCheck> {
        vec![
            Self::probe_region_command_dependency("Rscript", "command", false, &["--version"]),
            Self::probe_region_command_dependency(
                "apt-probeset-summarize",
                "command",
                false,
                &["--version"],
            ),
        ]
    }

    #[allow(clippy::too_many_arguments)]

    pub(super) fn probe_region_apt_import_replay_command(
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
    ) -> String {
        let mut args = vec![
            "arrays".to_string(),
            "import-apt-probe-region-output".to_string(),
            apt_summary_path.to_string(),
            annotation_path.to_string(),
            output_dir.to_string(),
        ];
        if let Some(value) = metadata_path {
            args.push("--metadata".to_string());
            args.push(value.to_string());
        }
        if let Some(value) = condition_column {
            args.push("--condition-column".to_string());
            args.push(value.to_string());
        }
        if let Some(value) = sample_column {
            args.push("--sample-column".to_string());
            args.push(value.to_string());
        }
        if let Some(value) = probe_intensity_path {
            args.push("--probe-intensity".to_string());
            args.push(value.to_string());
        }
        if let Some(value) = probe_id_column {
            args.push("--probe-id-column".to_string());
            args.push(value.to_string());
        }
        if let Some(value) = platform {
            args.push("--platform".to_string());
            args.push(value.to_string());
        }
        if let Some(value) = normalization {
            args.push("--normalization".to_string());
            args.push(value.to_string());
        }
        if let Some(value) = coordinate_system {
            args.push("--coordinate-system".to_string());
            args.push(value.to_string());
        }
        if let Some(value) = genome_build {
            args.push("--genome-build".to_string());
            args.push(value.to_string());
        }
        args.iter()
            .map(|arg| Self::probe_region_shell_quote(arg))
            .collect::<Vec<_>>()
            .join(" ")
    }

    pub(super) fn probe_region_optional_file_status(
        path: &Path,
        role: &str,
    ) -> Option<ProbeRegionFileStatus> {
        let status = Self::probe_region_file_status(&path.to_string_lossy(), role);
        status.exists.then_some(status)
    }

    pub(super) fn probe_region_json_file(path: &Path) -> Result<serde_json::Value, String> {
        let path_text = path.to_string_lossy();
        let file = File::open(path).map_err(|e| format!("could not open '{path_text}': {e}"))?;
        serde_json::from_reader(BufReader::new(file))
            .map_err(|e| format!("could not parse '{path_text}' as JSON: {e}"))
    }

    pub(super) fn probe_region_json_string(value: &serde_json::Value, key: &str) -> Option<String> {
        value
            .get(key)
            .and_then(|value| value.as_str())
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
    }

    pub(super) fn probe_region_json_string_array(
        value: &serde_json::Value,
        key: &str,
    ) -> Vec<String> {
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

    pub(super) fn probe_region_json_alias_string(
        value: &serde_json::Value,
        keys: &[&str],
    ) -> Option<String> {
        keys.iter()
            .find_map(|key| Self::probe_region_json_string(value, key))
    }

    pub(super) fn probe_region_json_coordinate_projections(
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

    pub(super) fn probe_region_extend_unique(target: &mut Vec<String>, values: Vec<String>) {
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

    pub(super) fn probe_region_extend_unique_coordinate_projections(
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

    pub(super) fn probe_region_warn_unexpected_json_schema(
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

    pub(super) fn probe_region_sample_ids_from_sample_table(
        path: &Path,
    ) -> Result<Vec<String>, String> {
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

    pub(super) fn probe_region_region_sample_columns(
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

    pub(super) fn probe_region_preview_values(values: BTreeSet<String>) -> Vec<String> {
        values.into_iter().take(64).collect()
    }

    pub(super) fn probe_region_platform_plan(platform: Option<&str>) -> ProbeRegionPlatformPlan {
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

    pub(super) fn probe_region_annotation_source(
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

    pub(super) fn probe_region_vendor_support_files(
        platform: &ProbeRegionPlatformPlan,
    ) -> Vec<ProbeRegionFileStatus> {
        if platform.normalized != "Clariom_D_Human" {
            return Vec::new();
        }
        [
            (
                CLARIOM_D_HUMAN_PROBESET_ZIP,
                CLARIOM_D_HUMAN_PROBESET_ZIP_ALIASES,
                "thermofisher_clariom_d_human_hg38_probeset_zip",
            ),
            (
                CLARIOM_D_HUMAN_TRANSCRIPT_ZIP,
                CLARIOM_D_HUMAN_TRANSCRIPT_ZIP_ALIASES,
                "thermofisher_clariom_d_human_hg38_transcript_zip",
            ),
        ]
        .iter()
        .map(|(file_name, aliases, role)| {
            Self::probe_region_vendor_support_file_status(
                Path::new(CLARIOM_D_HUMAN_VENDOR_SUPPORT_DIR),
                file_name,
                aliases,
                role,
            )
        })
        .collect()
    }

    pub(super) fn probe_region_vendor_support_file_status(
        base_dir: &Path,
        canonical_file_name: &str,
        aliases: &[&str],
        role: &str,
    ) -> ProbeRegionFileStatus {
        let canonical_path = base_dir.join(canonical_file_name);
        let mut canonical_status =
            Self::probe_region_file_status(&canonical_path.to_string_lossy(), role);
        if canonical_status.exists {
            return canonical_status;
        }

        for alias in aliases {
            let alias_path = base_dir.join(alias);
            let mut alias_status =
                Self::probe_region_file_status(&alias_path.to_string_lossy(), role);
            if alias_status.exists {
                alias_status.detail = Some(format!(
                    "Using Thermo Fisher download filename '{alias}'; canonical local filename is '{canonical_file_name}'"
                ));
                return alias_status;
            }
        }

        let accepted_names = std::iter::once(canonical_file_name)
            .chain(aliases.iter().copied())
            .collect::<Vec<_>>()
            .join(", ");
        canonical_status.detail = Some(format!(
            "Manual Thermo Fisher login download required; GENtle never auto-downloads this support file. Accepted local filenames: {accepted_names}"
        ));
        canonical_status
    }

    pub(super) fn probe_region_vendor_support_detail(
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
                "Thermo Fisher Clariom D hg38 support ZIPs are expected under {CLARIOM_D_HUMAN_VENDOR_SUPPORT_DIR} ({present}/{expected} present); these login-walled files must be placed manually, using either the canonical names or the browser-preserved Thermo download names"
            ))
        }
    }

    pub(super) fn probe_region_annotation_file_kind(path: &str) -> &'static str {
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

    pub(super) fn probe_region_metadata_plan(
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

    pub(super) fn probe_region_metadata_delimiter(path: &str) -> u8 {
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

    pub(super) fn probe_region_metadata_column_index(
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

    pub(super) fn probe_region_header_key(value: &str) -> String {
        value
            .chars()
            .filter(|ch| ch.is_ascii_alphanumeric())
            .flat_map(char::to_lowercase)
            .collect()
    }

    pub(super) fn probe_region_contrast_plan(
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

    pub(super) fn probe_region_backend_candidates(
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

    pub(super) fn probe_region_apt_library_plan(
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

    pub(super) fn probe_region_note_apt_library_file(
        plan: &mut ProbeRegionAptLibraryPlan,
        path: &Path,
    ) {
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

    pub(super) fn probe_region_apt_backend_detail(library: &ProbeRegionAptLibraryPlan) -> String {
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

    pub(super) fn probe_region_apt_suggested_command(
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

    pub(super) fn probe_region_oligo_suggested_command(
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

    pub(super) fn probe_region_shell_quote(value: &str) -> String {
        if value
            .chars()
            .all(|ch| ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.' | '/' | ':' | '='))
        {
            return value.to_string();
        }
        format!("'{}'", value.replace('\'', "'\"'\"'"))
    }

    pub(super) fn probe_region_dependency_present(
        dependencies: &[ProbeRegionDependencyCheck],
        name: &str,
    ) -> bool {
        dependencies
            .iter()
            .any(|row| row.name == name && row.status == "present")
    }

    pub(super) fn probe_region_command_dependency(
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

    pub(super) fn probe_region_r_package_dependency(
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

#[cfg(test)]
mod tests {
    #[test]
    fn probe_region_vendor_support_file_status_accepts_browser_download_filename() {
        let dir = tempfile::tempdir().expect("tempdir");
        std::fs::write(dir.path().join("browser_name.zip"), b"zip").expect("write alias");

        let status = crate::engine::GentleEngine::probe_region_vendor_support_file_status(
            dir.path(),
            "canonical_name.zip",
            &["browser_name.zip"],
            "vendor_zip",
        );

        assert!(status.exists);
        assert!(status.is_file);
        assert!(status.path.ends_with("browser_name.zip"));
        assert_eq!(status.role, "vendor_zip");
        assert!(
            status
                .detail
                .as_deref()
                .is_some_and(|detail| detail.contains("canonical_name.zip"))
        );
    }

    #[test]
    fn probe_region_vendor_support_file_status_reports_canonical_when_absent() {
        let dir = tempfile::tempdir().expect("tempdir");
        let status = crate::engine::GentleEngine::probe_region_vendor_support_file_status(
            dir.path(),
            "canonical_name.zip",
            &["browser_name.zip"],
            "vendor_zip",
        );

        assert!(!status.exists);
        assert!(status.path.ends_with("canonical_name.zip"));
        assert!(
            status
                .detail
                .as_deref()
                .is_some_and(|detail| detail.contains("browser_name.zip"))
        );
    }
}
