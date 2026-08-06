//! Resolve one portable gene-set manifest into mutually consistent web and print artifacts.

use gentle_protocol::{
    GENE_ISOFORM_ASSAY_PUBLICATION_PROJECTION_SCHEMA,
    GENE_ISOFORM_ASSAY_PUBLICATION_REQUEST_SCHEMA, GENE_ISOFORM_ASSAY_PUBLICATION_SCHEMA,
    GENE_SET_PUBLICATION_BUNDLE_MANIFEST_SCHEMA, GENE_SET_PUBLICATION_GENERATION_SCHEMA,
    GENE_SET_PUBLICATION_REPORT_SCHEMA, GENE_SET_PUBLICATION_REQUEST_SCHEMA,
    GeneIsoformAssayPublicationArtifact, GeneIsoformAssayPublicationBlock,
    GeneIsoformAssayPublicationBoundReport, GeneIsoformAssayPublicationGene,
    GeneIsoformAssayPublicationGeneRequest, GeneIsoformAssayPublicationGeneStatus,
    GeneIsoformAssayPublicationParameter, GeneIsoformAssayPublicationParameterOverride,
    GeneIsoformAssayPublicationProfile, GeneIsoformAssayPublicationProjectionReport,
    GeneIsoformAssayPublicationReport, GeneIsoformAssayPublicationReportRef,
    GeneIsoformAssayPublicationRequest, GeneSetPublicationBundleArtifact,
    GeneSetPublicationBundleManifest, GeneSetPublicationDownload, GeneSetPublicationFigure,
    GeneSetPublicationGene, GeneSetPublicationGenerationReport, GeneSetPublicationPrimerColumnMap,
    GeneSetPublicationPrimerRow, GeneSetPublicationReport, GeneSetPublicationRequest,
};
use gentle_render::{
    render_gene_isoform_assay_publication_gene, render_gene_isoform_assay_publication_index,
    render_gene_isoform_assay_publication_print, render_gene_set_publication_html,
    render_gene_set_publication_markdown,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::ErrorKind,
    path::{Component, Path, PathBuf},
    process::Command,
};

use crate::{
    digest_utils::sha256_prefixed_bytes,
    engine::protocol::{
        ExperimentalAssayHandoffReport, ExperimentalAssayReadinessState, OligoOrderLineProvenance,
    },
    svg_pdf::render_svg_file_to_pdf,
    svg_png::SvgPngRenderOptions,
    tool_overrides::{configured_or_env, resolve_tool_executable},
};

#[cfg(test)]
use crate::engine::protocol::ExperimentalAssayReadinessPolicy;

fn source_path(base: &Path, raw: &str) -> PathBuf {
    let path = Path::new(raw);
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base.join(path)
    }
}

fn safe_filename(path: &Path, what: &str) -> Result<String, String> {
    path.file_name()
        .and_then(|value| value.to_str())
        .filter(|value| !value.is_empty())
        .map(ToOwned::to_owned)
        .ok_or_else(|| format!("{what} must name a file: '{}'", path.display()))
}

fn bundle_target_key(path: &Path) -> String {
    path.to_string_lossy().replace('\\', "/").to_lowercase()
}

fn reserve_bundle_target(
    targets: &mut BTreeMap<String, String>,
    target: &Path,
    owner: String,
) -> Result<bool, String> {
    let key = bundle_target_key(target);
    if let Some(existing_owner) = targets.get(&key) {
        if existing_owner == &owner {
            return Ok(false);
        }
        return Err(format!(
            "Bundle path collision at '{}': {existing_owner} conflicts with {owner}",
            target.display()
        ));
    }
    targets.insert(key, owner);
    Ok(true)
}

fn copy_asset(
    source: &Path,
    target_dir: &Path,
    what: &str,
    targets: &mut BTreeMap<String, String>,
) -> Result<(PathBuf, String, bool), String> {
    if !source.is_file() {
        return Err(format!("{what} does not exist: '{}'", source.display()));
    }
    let filename = safe_filename(source, what)?;
    let target = target_dir.join(&filename);
    let canonical_source = fs::canonicalize(source)
        .map_err(|error| format!("Could not resolve {what} '{}': {error}", source.display()))?;
    let should_copy = reserve_bundle_target(
        targets,
        &target,
        format!("{what} source '{}'", canonical_source.display()),
    )?;
    if should_copy {
        fs::copy(source, &target).map_err(|error| {
            format!(
                "Could not copy {what} '{}' to '{}': {error}",
                source.display(),
                target.display()
            )
        })?;
    }
    Ok((target, filename, should_copy))
}

fn delimiter(raw: &str) -> Result<u8, String> {
    match raw {
        "\\t" | "\t" => Ok(b'\t'),
        value if value.as_bytes().len() == 1 => Ok(value.as_bytes()[0]),
        _ => Err("primer_table.delimiter must be one byte or \\t".to_string()),
    }
}

fn cell(
    row: &csv::StringRecord,
    headers: &csv::StringRecord,
    name: &str,
    required: bool,
) -> Result<String, String> {
    if name.is_empty() {
        return Ok(String::new());
    }
    let Some(index) = headers.iter().position(|header| header == name) else {
        return if required {
            Err(format!("Primer table is missing required column '{name}'"))
        } else {
            Ok(String::new())
        };
    };
    Ok(row.get(index).unwrap_or_default().trim().to_string())
}

fn read_primers(
    path: &Path,
    delimiter: u8,
    columns: &GeneSetPublicationPrimerColumnMap,
) -> Result<Vec<GeneSetPublicationPrimerRow>, String> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .from_path(path)
        .map_err(|error| format!("Could not read primer table '{}': {error}", path.display()))?;
    let headers = reader
        .headers()
        .map_err(|error| format!("Could not read primer headers: {error}"))?
        .clone();
    let mut primers = Vec::new();
    for record in reader.records() {
        let row = record.map_err(|error| {
            format!("Could not parse primer table '{}': {error}", path.display())
        })?;
        primers.push(GeneSetPublicationPrimerRow {
            gene: cell(&row, &headers, &columns.gene, true)?,
            pair_id: cell(&row, &headers, &columns.pair_id, true)?,
            role: cell(&row, &headers, &columns.role, true)?,
            status: cell(&row, &headers, &columns.status, true)?,
            forward_sequence_5_to_3: cell(&row, &headers, &columns.forward_sequence_5_to_3, true)?,
            reverse_sequence_5_to_3: cell(&row, &headers, &columns.reverse_sequence_5_to_3, true)?,
            forward_tm_c: cell(&row, &headers, &columns.forward_tm_c, false)?,
            reverse_tm_c: cell(&row, &headers, &columns.reverse_tm_c, false)?,
            cdna_specificity: cell(&row, &headers, &columns.cdna_specificity, false)?,
            genome_assessment: cell(&row, &headers, &columns.genome_assessment, false)?,
            note: cell(&row, &headers, &columns.note, false)?,
        });
    }
    if primers.is_empty() {
        return Err("Primer table contains no data rows".to_string());
    }
    Ok(primers)
}

fn run(command: &mut Command, description: &str, override_env: &str) -> Result<(), String> {
    let executable = command.get_program().to_string_lossy().into_owned();
    let output = command
        .output()
        .map_err(|error| {
            format!(
                "Could not run {description} with executable '{executable}': {error}. Install the tool or set {override_env}"
            )
        })?;
    if output.status.success() {
        return Ok(());
    }
    let stderr = String::from_utf8_lossy(&output.stderr).trim().to_string();
    let stdout = String::from_utf8_lossy(&output.stdout).trim().to_string();
    let detail = if !stderr.is_empty() {
        stderr
    } else if !stdout.is_empty() {
        stdout
    } else {
        "no diagnostic output".to_string()
    };
    Err(format!(
        "{description} failed with executable '{executable}' ({status}): {detail}. Override it with {override_env}",
        status = output.status
    ))
}

fn require_generated_file(path: &Path, description: &str) -> Result<(), String> {
    let metadata = fs::metadata(path).map_err(|error| {
        format!(
            "{description} reported success but did not create '{}': {error}",
            path.display()
        )
    })?;
    if !metadata.is_file() || metadata.len() == 0 {
        return Err(format!(
            "{description} reported success but produced an empty or non-file output at '{}'",
            path.display()
        ));
    }
    Ok(())
}

fn remove_generated_file_if_present(path: &Path) -> Result<(), String> {
    match fs::remove_file(path) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == ErrorKind::NotFound => Ok(()),
        Err(error) => Err(format!(
            "Could not remove stale generated file '{}': {error}",
            path.display()
        )),
    }
}

fn cleanup_generated_file(path: &Path, warnings: &mut Vec<String>) {
    if let Err(error) = fs::remove_file(path)
        && error.kind() != ErrorKind::NotFound
    {
        warnings.push(format!(
            "Could not remove temporary file '{}': {error}",
            path.display()
        ));
    }
}

fn validate_leaf_filename(value: &str, label: &str) -> Result<(), String> {
    let mut components = Path::new(value).components();
    let is_single_normal_component =
        matches!(components.next(), Some(Component::Normal(_))) && components.next().is_none();
    if value.is_empty()
        || value.trim() != value
        || value == "."
        || value == ".."
        || value.contains('/')
        || value.contains('\\')
        || Path::new(value).is_absolute()
        || !is_single_normal_component
    {
        return Err(format!(
            "{label} must be one filename without directory components: '{value}'"
        ));
    }
    Ok(())
}

fn validate_output_filenames(request: &GeneSetPublicationRequest) -> Result<(), String> {
    const RESERVED: &[&str] = &[
        "assembled-uncompressed.pdf",
        "bundle-manifest.json",
        "data",
        "favicon.ico",
        "figures",
        "generation-report.json",
        "narrative.pdf",
        "resolved-report.json",
        ".gentle-publication-optimized.pdf",
    ];
    let names = [
        ("html_filename", request.html_filename.as_str()),
        ("markdown_filename", request.markdown_filename.as_str()),
        ("pdf_filename", request.pdf_filename.as_str()),
    ];
    let mut seen = BTreeMap::new();
    for (label, value) in names {
        validate_leaf_filename(value, label)?;
        let key = value.to_lowercase();
        if RESERVED
            .iter()
            .any(|reserved| reserved.eq_ignore_ascii_case(value))
        {
            return Err(format!(
                "{label} uses reserved bundle path '{value}'; choose another filename"
            ));
        }
        if let Some(existing_label) = seen.insert(key, label) {
            return Err(format!(
                "{existing_label} and {label} must use distinct filenames: '{value}'"
            ));
        }
    }
    Ok(())
}

fn validate_request(request: &GeneSetPublicationRequest) -> Result<(), String> {
    if request.schema != GENE_SET_PUBLICATION_REQUEST_SCHEMA {
        return Err(format!("Unsupported schema '{}'", request.schema));
    }
    if request.report_id.trim().is_empty() || request.title.trim().is_empty() {
        return Err("report_id and title are required".to_string());
    }
    if request.genes.is_empty() {
        return Err("At least one gene is required".to_string());
    }
    validate_output_filenames(request)?;
    let mut genes = BTreeSet::new();
    let mut figure_ids = BTreeSet::new();
    for gene in &request.genes {
        if gene.gene_symbol.trim().is_empty() || !genes.insert(gene.gene_symbol.clone()) {
            return Err(format!(
                "Gene symbols must be non-empty and unique: '{}'",
                gene.gene_symbol
            ));
        }
        for figure in &gene.figures {
            validate_leaf_filename(&figure.figure_id, "figure_id")?;
            if !figure_ids.insert(figure.figure_id.clone()) {
                return Err(format!(
                    "Figure IDs must be non-empty and unique: '{}'",
                    figure.figure_id
                ));
            }
        }
    }
    Ok(())
}

/// Generate a self-contained responsive HTML report and optional companion PDF.
pub fn generate_gene_set_publication(
    request_path: &Path,
    output_directory: &Path,
    generate_pdf: bool,
) -> Result<GeneSetPublicationGenerationReport, String> {
    let request_text = fs::read_to_string(request_path).map_err(|error| {
        format!(
            "Could not read request '{}': {error}",
            request_path.display()
        )
    })?;
    let request_sha256 = sha256_prefixed_bytes(request_text.as_bytes());
    let request: GeneSetPublicationRequest =
        serde_json::from_str(&request_text).map_err(|error| {
            format!(
                "Could not parse request '{}': {error}",
                request_path.display()
            )
        })?;
    validate_request(&request)?;
    let base = request_path.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_directory)
        .map_err(|error| format!("Could not create '{}': {error}", output_directory.display()))?;
    let figures_dir = output_directory.join("figures");
    let data_dir = output_directory.join("data");
    fs::create_dir_all(&figures_dir).map_err(|error| error.to_string())?;
    fs::create_dir_all(&data_dir).map_err(|error| error.to_string())?;
    let mut bundle_targets = BTreeMap::new();

    let primer_source = source_path(base, &request.primer_table.path);
    let primers = read_primers(
        &primer_source,
        delimiter(&request.primer_table.delimiter)?,
        &request.primer_table.columns,
    )?;
    let (_, primer_name, _) = copy_asset(
        &primer_source,
        &data_dir,
        "primer table",
        &mut bundle_targets,
    )?;
    let mut copied_files = vec![format!("data/{primer_name}")];
    let mut pdf_pages = Vec::new();
    let mut genes = Vec::new();
    for gene in &request.genes {
        let mut figures = Vec::new();
        for figure in &gene.figures {
            let source = source_path(base, &figure.source_path);
            let (copied, name, was_copied) =
                copy_asset(&source, &figures_dir, "figure", &mut bundle_targets)?;
            if was_copied {
                copied_files.push(format!("figures/{name}"));
            }
            let pdf_path = if generate_pdf && figure.include_in_pdf {
                let pdf_target = if let Some(pdf_source) = &figure.pdf_source_path {
                    let source = source_path(base, pdf_source);
                    let (target, name, was_copied) =
                        copy_asset(&source, &figures_dir, "figure PDF", &mut bundle_targets)?;
                    if was_copied {
                        copied_files.push(format!("figures/{name}"));
                    }
                    target
                } else {
                    let target = figures_dir.join(format!("{}.pdf", figure.figure_id));
                    reserve_bundle_target(
                        &mut bundle_targets,
                        &target,
                        format!("generated PDF for figure '{}'", figure.figure_id),
                    )?;
                    render_svg_file_to_pdf(
                        &copied,
                        &target,
                        SvgPngRenderOptions {
                            scale: 2.0,
                            drop_dotplot_metadata: false,
                        },
                    )?;
                    copied_files.push(format!("figures/{}.pdf", figure.figure_id));
                    target
                };
                pdf_pages.push(pdf_target.clone());
                Some(format!(
                    "figures/{}",
                    safe_filename(&pdf_target, "figure PDF")?
                ))
            } else {
                None
            };
            figures.push(GeneSetPublicationFigure {
                figure_id: figure.figure_id.clone(),
                kind: figure.kind.clone(),
                label: figure.label.clone(),
                web_path: format!("figures/{name}"),
                pdf_path,
                caption: figure.caption.clone(),
                alt_text: figure.alt_text.clone(),
                include_in_pdf: figure.include_in_pdf,
            });
        }
        genes.push(GeneSetPublicationGene {
            gene_symbol: gene.gene_symbol.clone(),
            heading: gene.heading.clone(),
            summary: gene.summary.clone(),
            metrics: gene.metrics.clone(),
            figures,
        });
    }
    let mut downloads = Vec::new();
    for download in &request.downloads {
        let source = source_path(base, &download.source_path);
        let (_, name, was_copied) =
            copy_asset(&source, &data_dir, "download", &mut bundle_targets)?;
        if was_copied {
            copied_files.push(format!("data/{name}"));
        }
        downloads.push(GeneSetPublicationDownload {
            label: download.label.clone(),
            web_path: format!("data/{name}"),
            description: download.description.clone(),
        });
    }
    let favicon_path = if let Some(raw) = &request.favicon_source_path {
        let source = source_path(base, raw);
        fs::copy(&source, output_directory.join("favicon.ico"))
            .map_err(|error| format!("Could not copy favicon '{}': {error}", source.display()))?;
        copied_files.push("favicon.ico".to_string());
        Some("favicon.ico".to_string())
    } else {
        None
    };
    let mut report = GeneSetPublicationReport {
        schema: GENE_SET_PUBLICATION_REPORT_SCHEMA.to_string(),
        report_id: request.report_id.clone(),
        title: request.title.clone(),
        subtitle: request.subtitle.clone(),
        eyebrow: request.eyebrow.clone(),
        generated_date: request.generated_date.clone(),
        summary: request.summary.clone(),
        conditions: request.conditions.clone(),
        metrics: request.metrics.clone(),
        sections: request.sections.clone(),
        primers,
        genes,
        downloads,
        provenance: request.provenance.clone(),
        footer: request.footer.clone(),
        favicon_path,
        pdf_path: None,
    };
    let html_path = output_directory.join(&request.html_filename);
    let markdown_path = output_directory.join(&request.markdown_filename);
    let resolved_path = output_directory.join("resolved-report.json");
    fs::write(&html_path, render_gene_set_publication_html(&report))
        .map_err(|error| error.to_string())?;
    fs::write(
        &markdown_path,
        render_gene_set_publication_markdown(&report),
    )
    .map_err(|error| error.to_string())?;
    fs::write(
        &resolved_path,
        serde_json::to_string_pretty(&report).map_err(|error| error.to_string())?,
    )
    .map_err(|error| error.to_string())?;
    let (mut bundle_manifest_path, mut bundle_manifest_sha256, mut artifacts) =
        write_gene_set_bundle_manifest(
            output_directory,
            &request.report_id,
            &request_sha256,
            &copied_files,
            &request.html_filename,
            &request.markdown_filename,
            None,
        )?;
    let mut warnings = Vec::new();
    let pdf_path = if generate_pdf {
        let narrative = output_directory.join("narrative.pdf");
        let assembled_pdf = output_directory.join("assembled-uncompressed.pdf");
        let optimized_pdf = output_directory.join(".gentle-publication-optimized.pdf");
        for path in [&narrative, &assembled_pdf, &optimized_pdf] {
            remove_generated_file_if_present(path)?;
        }
        let pandoc = resolve_tool_executable("GENTLE_PANDOC_BIN", "pandoc");
        run(
            Command::new(pandoc)
                .current_dir(output_directory)
                .arg(&request.markdown_filename)
                .args([
                    "--pdf-engine=lualatex",
                    "--toc",
                    "--number-sections",
                    "-o",
                    "narrative.pdf",
                ]),
            "pandoc PDF generation",
            "GENTLE_PANDOC_BIN",
        )?;
        require_generated_file(&narrative, "pandoc PDF generation")?;
        let final_pdf = output_directory.join(&request.pdf_filename);
        let unite = resolve_tool_executable("GENTLE_PDFUNITE_BIN", "pdfunite");
        let mut command = Command::new(unite);
        command.arg(&narrative);
        for page in &pdf_pages {
            command.arg(page);
        }
        command.arg(&assembled_pdf);
        run(&mut command, "PDF page assembly", "GENTLE_PDFUNITE_BIN")?;
        require_generated_file(&assembled_pdf, "PDF page assembly")?;
        let ghostscript = resolve_tool_executable("GENTLE_GHOSTSCRIPT_BIN", "gs");
        run(
            Command::new(ghostscript)
                .args([
                    "-q",
                    "-dNOPAUSE",
                    "-dBATCH",
                    "-dSAFER",
                    "-sDEVICE=pdfwrite",
                    "-dCompatibilityLevel=1.7",
                    "-dPDFSETTINGS=/prepress",
                ])
                .arg(format!("-sOutputFile={}", optimized_pdf.display()))
                .arg(&assembled_pdf),
            "Ghostscript PDF optimization",
            "GENTLE_GHOSTSCRIPT_BIN",
        )?;
        require_generated_file(&optimized_pdf, "Ghostscript PDF optimization")?;
        fs::copy(&optimized_pdf, &final_pdf).map_err(|error| {
            format!(
                "Could not install optimized PDF '{}' as '{}': {error}",
                optimized_pdf.display(),
                final_pdf.display()
            )
        })?;
        require_generated_file(&final_pdf, "PDF bundle installation")?;
        report.pdf_path = Some(request.pdf_filename.clone());
        fs::write(&html_path, render_gene_set_publication_html(&report))
            .map_err(|error| error.to_string())?;
        fs::write(
            &resolved_path,
            serde_json::to_string_pretty(&report).map_err(|error| error.to_string())?,
        )
        .map_err(|error| error.to_string())?;
        let refreshed_manifest = write_gene_set_bundle_manifest(
            output_directory,
            &request.report_id,
            &request_sha256,
            &copied_files,
            &request.html_filename,
            &request.markdown_filename,
            report.pdf_path.as_deref(),
        )?;
        bundle_manifest_path = refreshed_manifest.0;
        bundle_manifest_sha256 = refreshed_manifest.1;
        artifacts = refreshed_manifest.2;
        cleanup_generated_file(&narrative, &mut warnings);
        cleanup_generated_file(&assembled_pdf, &mut warnings);
        cleanup_generated_file(&optimized_pdf, &mut warnings);
        Some(final_pdf.to_string_lossy().into_owned())
    } else {
        report.pdf_path = None;
        None
    };
    let generation = GeneSetPublicationGenerationReport {
        schema: GENE_SET_PUBLICATION_GENERATION_SCHEMA.to_string(),
        request_path: request_path.to_string_lossy().into_owned(),
        output_directory: output_directory.to_string_lossy().into_owned(),
        html_path: html_path.to_string_lossy().into_owned(),
        markdown_path: markdown_path.to_string_lossy().into_owned(),
        resolved_report_path: resolved_path.to_string_lossy().into_owned(),
        pdf_path,
        copied_files,
        bundle_manifest_path: bundle_manifest_path.to_string_lossy().into_owned(),
        bundle_manifest_sha256,
        artifacts,
        warnings,
    };
    fs::write(
        output_directory.join("generation-report.json"),
        serde_json::to_string_pretty(&generation).map_err(|error| error.to_string())?,
    )
    .map_err(|error| error.to_string())?;
    Ok(generation)
}

fn load_publication_bound_report(
    base: &Path,
    reference: &GeneIsoformAssayPublicationReportRef,
    expected_schema: &str,
    id_fields: &[&str],
) -> Result<GeneIsoformAssayPublicationBoundReport, String> {
    if reference.path.trim().is_empty() || reference.expected_sha256.trim().is_empty() {
        return Err("Publication report references require path and expected_sha256".to_string());
    }
    let path = source_path(base, &reference.path);
    let bytes = fs::read(&path)
        .map_err(|error| format!("Could not read report '{}': {error}", path.display()))?;
    let sha256 = sha256_prefixed_bytes(&bytes);
    if sha256 != reference.expected_sha256 {
        return Err(format!(
            "Report '{}' has digest '{}' but '{}' was required",
            path.display(),
            sha256,
            reference.expected_sha256
        ));
    }
    let value: serde_json::Value = serde_json::from_slice(&bytes)
        .map_err(|error| format!("Could not parse report '{}': {error}", path.display()))?;
    let schema = value
        .get("schema")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if schema != expected_schema {
        return Err(format!(
            "Report '{}' uses schema '{}'; expected '{}'",
            path.display(),
            schema,
            expected_schema
        ));
    }
    let report_id = id_fields
        .iter()
        .find_map(|field| value.get(*field).and_then(serde_json::Value::as_str))
        .filter(|value| !value.trim().is_empty())
        .ok_or_else(|| {
            format!(
                "Report '{}' has none of the required identity fields {:?}",
                path.display(),
                id_fields
            )
        })?;
    Ok(GeneIsoformAssayPublicationBoundReport {
        source_path: reference.path.clone(),
        sha256,
        schema: schema.to_string(),
        report_id: report_id.to_string(),
        value,
    })
}

fn publication_gene_slug(symbol: &str) -> Result<String, String> {
    let slug = symbol
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
                ch
            } else {
                '-'
            }
        })
        .collect::<String>();
    if slug.is_empty() || slug == "." || slug == ".." {
        return Err(format!(
            "Gene symbol '{symbol}' cannot form a safe page name"
        ));
    }
    Ok(slug)
}

fn publication_parameter_entries(
    genes: &[GeneIsoformAssayPublicationGene],
) -> (
    Vec<GeneIsoformAssayPublicationParameter>,
    Vec<GeneIsoformAssayPublicationParameterOverride>,
) {
    let mut by_name: BTreeMap<String, Vec<(usize, serde_json::Value)>> = BTreeMap::new();
    for (gene_index, gene) in genes.iter().enumerate() {
        if let Some(policy) = gene
            .study_plan
            .value
            .pointer("/normalized_request/policy")
            .and_then(serde_json::Value::as_object)
        {
            for (name, value) in policy {
                by_name
                    .entry(name.clone())
                    .or_default()
                    .push((gene_index, value.clone()));
            }
        }
    }
    let mut common = vec![];
    let mut overrides = vec![];
    for (name, values) in by_name {
        let all_genes_present = values.len() == genes.len();
        let first = values.first().map(|(_, value)| value);
        let all_equal = first.is_some_and(|first| values.iter().all(|(_, value)| value == first));
        if all_genes_present && all_equal {
            common.push(GeneIsoformAssayPublicationParameter {
                name: name.clone(),
                value: first.cloned().unwrap_or(serde_json::Value::Null),
                source_pointers: values
                    .iter()
                    .map(|(gene_index, _)| {
                        format!(
                            "/genes/{gene_index}/study_plan/value/normalized_request/policy/{name}"
                        )
                    })
                    .collect(),
            });
        } else {
            for (gene_index, value) in values {
                let gene = &genes[gene_index];
                let reason = gene
                    .study_plan
                    .value
                    .pointer("/profile_override/reason")
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or(
                        "GENtle retained this gene-specific effective policy value from the approved normalized request.",
                    )
                    .to_string();
                overrides.push(GeneIsoformAssayPublicationParameterOverride {
                    gene_symbol: gene.gene_symbol.clone(),
                    name: name.clone(),
                    common_value: None,
                    effective_value: value,
                    reason,
                    source_pointer: format!(
                        "/genes/{gene_index}/study_plan/value/normalized_request/policy/{name}"
                    ),
                });
            }
        }
    }
    (common, overrides)
}

fn validate_order_form_handoff_binding(
    gene: &GeneIsoformAssayPublicationGene,
) -> Result<(), String> {
    let handoffs = gene
        .handoffs
        .iter()
        .map(|handoff| {
            let report: ExperimentalAssayHandoffReport =
                serde_json::from_value(handoff.value.clone()).map_err(|error| {
                    format!(
                        "Could not validate experimental handoff '{}': {error}",
                        handoff.report_id
                    )
                })?;
            let policy_text = serde_json::to_string(&report.policy).map_err(|error| {
                format!(
                    "Could not serialize readiness policy for handoff '{}': {error}",
                    handoff.report_id
                )
            })?;
            let policy_sha256 = sha256_prefixed_bytes(policy_text.as_bytes());
            let expected_policy_id = format!(
                "readiness_policy_sha256_{}",
                policy_sha256.trim_start_matches("sha256:")
            );
            if report.policy_id != expected_policy_id {
                return Err(format!(
                    "Experimental handoff '{}' policy_id '{}' does not match embedded policy '{}'",
                    handoff.report_id, report.policy_id, expected_policy_id
                ));
            }
            Ok((handoff, report, policy_sha256))
        })
        .collect::<Result<Vec<_>, String>>()?;
    for order_form in &gene.order_forms {
        let lines = order_form
            .value
            .get("line_items")
            .and_then(serde_json::Value::as_array)
            .ok_or_else(|| {
                format!(
                    "Order form '{}' has no line_items array",
                    order_form.report_id
                )
            })?;
        for line in lines {
            let provenance_value = line.get("provenance").ok_or_else(|| {
                format!(
                    "Order form '{}' line lacks readiness-bound provenance",
                    order_form.report_id
                )
            })?;
            let provenance: OligoOrderLineProvenance =
                serde_json::from_value(provenance_value.clone()).map_err(|error| {
                    format!(
                        "Order form '{}' line has invalid readiness provenance: {error}",
                        order_form.report_id
                    )
                })?;
            let report_sha256 = provenance.report_sha256.as_deref().unwrap_or_default();
            let policy_id = provenance
                .readiness_policy_id
                .as_deref()
                .unwrap_or_default();
            let policy_sha256 = provenance
                .readiness_policy_sha256
                .as_deref()
                .unwrap_or_default();
            let matching_handoff = handoffs.iter().find(|(bound, report, digest)| {
                bound.report_id == provenance.report_id
                    && bound.sha256 == report_sha256
                    && report.policy_id == policy_id
                    && digest == policy_sha256
            });
            let Some((_, handoff, _)) = matching_handoff else {
                return Err(format!(
                    "Order form '{}' line '{}' is not bound to an included order-ready experimental handoff and policy",
                    order_form.report_id,
                    line.get("line_id")
                        .and_then(serde_json::Value::as_str)
                        .unwrap_or("unknown")
                ));
            };
            if provenance.source_kind != "experimental_assay_handoff"
                || provenance.report_schema != "gentle.experimental_assay_handoff.v1"
                || provenance.readiness_state.as_deref() != Some("order_ready")
            {
                return Err(format!(
                    "Order form '{}' line '{}' does not identify an order-ready experimental-assay handoff source",
                    order_form.report_id,
                    line.get("line_id")
                        .and_then(serde_json::Value::as_str)
                        .unwrap_or("unknown")
                ));
            }
            let readiness_row = provenance.readiness_row.as_ref().ok_or_else(|| {
                format!(
                    "Order form '{}' line '{}' does not retain its readiness row",
                    order_form.report_id,
                    line.get("line_id")
                        .and_then(serde_json::Value::as_str)
                        .unwrap_or("unknown")
                )
            })?;
            let card_id = provenance.readiness_card_id.as_deref().unwrap_or_default();
            let assay_id = provenance.assay_id.as_deref().unwrap_or_default();
            let pair_id = provenance.pair_id.as_deref().unwrap_or_default();
            let oligo_id = provenance.oligo_id.as_deref().unwrap_or_default();
            let sequence = line
                .get("sequence_5_to_3")
                .and_then(serde_json::Value::as_str)
                .unwrap_or_default();
            let row_matches_provenance = readiness_row.order_ready
                && readiness_row.readiness_state == ExperimentalAssayReadinessState::OrderReady
                && readiness_row.card_id == card_id
                && readiness_row.assay_id == assay_id
                && readiness_row.pair_id == pair_id
                && readiness_row
                    .oligo_ids
                    .iter()
                    .any(|value| value == oligo_id)
                && readiness_row
                    .sequences_5_to_3
                    .iter()
                    .any(|value| value == sequence);
            let row_matches_handoff = handoff.order_readiness_table.iter().any(|row| {
                row.order_ready
                    && row.readiness_state == ExperimentalAssayReadinessState::OrderReady
                    && row.card_id == readiness_row.card_id
                    && row.assay_id == readiness_row.assay_id
                    && row.pair_id == readiness_row.pair_id
                    && row.oligo_ids == readiness_row.oligo_ids
                    && row.sequences_5_to_3 == readiness_row.sequences_5_to_3
            });
            if !row_matches_provenance || !row_matches_handoff {
                return Err(format!(
                    "Order form '{}' line '{}' readiness row does not match its provenance or approved handoff",
                    order_form.report_id,
                    line.get("line_id")
                        .and_then(serde_json::Value::as_str)
                        .unwrap_or("unknown")
                ));
            }
        }
    }
    Ok(())
}

/// Decide whether one dossier gene shows results, is still expected, or could
/// not be addressed automatically.
///
/// An omitted status is derived from the presence of assay handoffs so that
/// dossiers written before this field existed report honestly instead of
/// presenting a gene without results as finished. A declared status is checked
/// against the gene's actual contents rather than trusted, because the whole
/// point of the label is that a reader can rely on it.
fn resolve_publication_gene_status(
    request: &GeneIsoformAssayPublicationGeneRequest,
    handoffs: &[GeneIsoformAssayPublicationBoundReport],
) -> Result<(GeneIsoformAssayPublicationGeneStatus, Option<String>), String> {
    let symbol = request.gene_symbol.trim();
    let reason = request
        .status_reason
        .as_deref()
        .map(str::trim)
        .filter(|reason| !reason.is_empty())
        .map(str::to_string);
    let status = match request.status {
        Some(declared) => declared,
        None if handoffs.is_empty() => GeneIsoformAssayPublicationGeneStatus::Pending,
        None => GeneIsoformAssayPublicationGeneStatus::Resolved,
    };
    match status {
        GeneIsoformAssayPublicationGeneStatus::Resolved if handoffs.is_empty() => {
            return Err(format!(
                "Gene '{symbol}' is declared 'resolved' but carries no experimental assay handoff; declare 'pending' or 'unresolved' with a status_reason instead of publishing it as complete"
            ));
        }
        GeneIsoformAssayPublicationGeneStatus::Resolved if reason.is_some() => {
            return Err(format!(
                "Gene '{symbol}' is 'resolved' and its status_reason would never be shown; remove the reason or declare 'pending' or 'unresolved'"
            ));
        }
        GeneIsoformAssayPublicationGeneStatus::Unresolved if reason.is_none() => {
            return Err(format!(
                "Gene '{symbol}' is declared 'unresolved' and requires a status_reason stating why the established automatism could not address it"
            ));
        }
        GeneIsoformAssayPublicationGeneStatus::Pending
        | GeneIsoformAssayPublicationGeneStatus::Unresolved
            if !handoffs.is_empty() =>
        {
            return Err(format!(
                "Gene '{symbol}' is declared '{}' but carries {} experimental assay handoff(s); publish it as 'resolved' or remove the handoffs",
                status.as_str(),
                handoffs.len()
            ));
        }
        _ => {}
    }
    Ok((status, reason))
}

fn validate_handoff_plan_binding(gene: &GeneIsoformAssayPublicationGene) -> Result<(), String> {
    let mut planned_panels = BTreeMap::new();
    for step in gene
        .study_plan
        .value
        .get("planned_operations")
        .and_then(serde_json::Value::as_array)
        .into_iter()
        .flatten()
    {
        let Some(payload) = step
            .get("operation")
            .and_then(|operation| operation.get("DesignTranscriptAssayPanel"))
        else {
            continue;
        };
        let report_id = payload
            .get("report_id")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        let seq_id = payload
            .get("seq_id")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        let source_feature_id = payload
            .get("source_feature_id")
            .and_then(serde_json::Value::as_u64)
            .and_then(|value| usize::try_from(value).ok());
        if report_id.is_empty() || seq_id.is_empty() || source_feature_id.is_none() {
            return Err(format!(
                "Study plan '{}' has an incomplete DesignTranscriptAssayPanel payload",
                gene.study_plan.report_id
            ));
        }
        if planned_panels
            .insert(
                report_id.to_string(),
                (seq_id.to_string(), source_feature_id.unwrap_or_default()),
            )
            .is_some()
        {
            return Err(format!(
                "Study plan '{}' repeats panel report id '{}'",
                gene.study_plan.report_id, report_id
            ));
        }
    }
    for bound in &gene.handoffs {
        let handoff: ExperimentalAssayHandoffReport = serde_json::from_value(bound.value.clone())
            .map_err(|error| {
            format!(
                "Could not validate experimental handoff '{}': {error}",
                bound.report_id
            )
        })?;
        let Some((seq_id, source_feature_id)) = planned_panels.get(&handoff.source_panel_report_id)
        else {
            return Err(format!(
                "Experimental handoff '{}' source panel '{}' is absent from approved study plan '{}'",
                bound.report_id, handoff.source_panel_report_id, gene.study_plan.report_id
            ));
        };
        if handoff.source_seq_id != *seq_id || handoff.source_feature_id != *source_feature_id {
            return Err(format!(
                "Experimental handoff '{}' source sequence/feature does not match approved panel '{}'",
                bound.report_id, handoff.source_panel_report_id
            ));
        }
    }
    Ok(())
}

fn write_publication_order_sheet(
    path: &Path,
    gene: &GeneIsoformAssayPublicationGene,
) -> Result<bool, String> {
    if gene.order_forms.is_empty() {
        return Ok(false);
    }
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|error| format!("Could not create order sheet '{}': {error}", path.display()))?;
    writer
        .write_record([
            "form_id",
            "line_id",
            "line_no",
            "name",
            "role",
            "sequence_5_to_3",
            "length_nt",
            "modifications",
            "scale",
            "purification",
            "handoff_sha256",
            "readiness_policy_id",
            "readiness_policy_sha256",
            "readiness_card_id",
            "assay_id",
            "pair_id",
            "oligo_id",
        ])
        .map_err(|error| error.to_string())?;
    for form in &gene.order_forms {
        for line in form
            .value
            .get("line_items")
            .and_then(serde_json::Value::as_array)
            .into_iter()
            .flatten()
        {
            let provenance = line.get("provenance").unwrap_or(&serde_json::Value::Null);
            writer
                .write_record([
                    form.report_id.clone(),
                    json_string(line.get("line_id")),
                    json_string(line.get("line_no")),
                    json_string(line.get("name")),
                    json_string(line.get("role")),
                    json_string(line.get("sequence_5_to_3")),
                    json_string(line.get("length_nt")),
                    json_string(line.get("modifications")),
                    json_string(line.get("scale")),
                    json_string(line.get("purification")),
                    json_string(provenance.get("report_sha256")),
                    json_string(provenance.get("readiness_policy_id")),
                    json_string(provenance.get("readiness_policy_sha256")),
                    json_string(provenance.get("readiness_card_id")),
                    json_string(provenance.get("assay_id")),
                    json_string(provenance.get("pair_id")),
                    json_string(provenance.get("oligo_id")),
                ])
                .map_err(|error| error.to_string())?;
        }
    }
    writer.flush().map_err(|error| error.to_string())?;
    Ok(true)
}

fn json_string(value: Option<&serde_json::Value>) -> String {
    match value {
        Some(serde_json::Value::String(value)) => value.clone(),
        Some(serde_json::Value::Null) | None => String::new(),
        Some(value) => serde_json::to_string(value).unwrap_or_default(),
    }
}

fn gene_set_bundle_artifact(
    output_directory: &Path,
    relative_path: &str,
) -> Result<GeneSetPublicationBundleArtifact, String> {
    let path = output_directory.join(relative_path);
    let bytes = fs::read(&path).map_err(|error| {
        format!(
            "Could not hash gene-set publication artifact '{}': {error}",
            path.display()
        )
    })?;
    Ok(GeneSetPublicationBundleArtifact {
        path: relative_path.replace('\\', "/"),
        sha256: sha256_prefixed_bytes(&bytes),
        media_type: publication_media_type(&path).to_string(),
    })
}

fn write_gene_set_bundle_manifest(
    output_directory: &Path,
    report_id: &str,
    request_sha256: &str,
    copied_files: &[String],
    html_filename: &str,
    markdown_filename: &str,
    pdf_filename: Option<&str>,
) -> Result<(PathBuf, String, Vec<GeneSetPublicationBundleArtifact>), String> {
    let mut paths = copied_files.iter().cloned().collect::<BTreeSet<_>>();
    paths.insert(html_filename.to_string());
    paths.insert(markdown_filename.to_string());
    paths.insert("resolved-report.json".to_string());
    if let Some(pdf_filename) = pdf_filename {
        paths.insert(pdf_filename.to_string());
    }
    let artifacts = paths
        .iter()
        .map(|path| gene_set_bundle_artifact(output_directory, path))
        .collect::<Result<Vec<_>, _>>()?;
    let manifest = GeneSetPublicationBundleManifest {
        schema: GENE_SET_PUBLICATION_BUNDLE_MANIFEST_SCHEMA.to_string(),
        report_id: report_id.to_string(),
        request_sha256: request_sha256.to_string(),
        artifacts: artifacts.clone(),
    };
    let bytes = serde_json::to_vec_pretty(&manifest).map_err(|error| error.to_string())?;
    let path = output_directory.join("bundle-manifest.json");
    fs::write(&path, &bytes).map_err(|error| {
        format!(
            "Could not write gene-set publication bundle manifest '{}': {error}",
            path.display()
        )
    })?;
    Ok((path, sha256_prefixed_bytes(&bytes), artifacts))
}

fn publication_artifact(
    path: &Path,
    media_type: &str,
) -> Result<GeneIsoformAssayPublicationArtifact, String> {
    let bytes = fs::read(path).map_err(|error| {
        format!(
            "Could not hash publication artifact '{}': {error}",
            path.display()
        )
    })?;
    Ok(GeneIsoformAssayPublicationArtifact {
        path: path.to_string_lossy().into_owned(),
        sha256: sha256_prefixed_bytes(&bytes),
        media_type: media_type.to_string(),
    })
}

fn publication_media_type(path: &Path) -> &'static str {
    match path
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase()
        .as_str()
    {
        "svg" => "image/svg+xml",
        "png" => "image/png",
        "jpg" | "jpeg" => "image/jpeg",
        "pdf" => "application/pdf",
        "html" | "htm" => "text/html",
        "md" | "markdown" => "text/markdown",
        "json" => "application/json",
        "tsv" => "text/tab-separated-values",
        "csv" => "text/csv",
        "ico" => "image/x-icon",
        _ => "application/octet-stream",
    }
}

fn choose_publication_browser(
    configured: Option<&str>,
    installed_candidates: &[PathBuf],
    default_bin: &str,
) -> String {
    if let Some(configured) = configured.map(str::trim).filter(|value| !value.is_empty()) {
        return configured.to_string();
    }
    installed_candidates
        .iter()
        .find(|path| path.is_file())
        .map(|path| path.to_string_lossy().into_owned())
        .unwrap_or_else(|| default_bin.to_string())
}

fn publication_browser_candidates() -> Vec<PathBuf> {
    let mut candidates = Vec::new();
    #[cfg(target_os = "macos")]
    candidates.extend([
        PathBuf::from("/Applications/Chromium.app/Contents/MacOS/Chromium"),
        PathBuf::from("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"),
        PathBuf::from("/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge"),
    ]);
    #[cfg(target_os = "linux")]
    candidates.extend([
        PathBuf::from("/usr/bin/chromium"),
        PathBuf::from("/usr/bin/chromium-browser"),
        PathBuf::from("/usr/bin/google-chrome"),
        PathBuf::from("/usr/bin/google-chrome-stable"),
        PathBuf::from("/usr/bin/microsoft-edge"),
        PathBuf::from("/usr/bin/microsoft-edge-stable"),
        PathBuf::from("/snap/bin/chromium"),
    ]);
    #[cfg(target_os = "windows")]
    for root in ["PROGRAMFILES", "PROGRAMFILES(X86)", "LOCALAPPDATA"] {
        if let Some(root) = std::env::var_os(root) {
            let root = PathBuf::from(root);
            candidates.extend([
                root.join("Google/Chrome/Application/chrome.exe"),
                root.join("Chromium/Application/chrome.exe"),
                root.join("Microsoft/Edge/Application/msedge.exe"),
            ]);
        }
    }

    #[cfg(target_os = "windows")]
    const PATH_NAMES: &[&str] = &["chrome.exe", "chromium.exe", "msedge.exe"];
    #[cfg(not(target_os = "windows"))]
    const PATH_NAMES: &[&str] = &[
        "chromium",
        "chromium-browser",
        "google-chrome",
        "google-chrome-stable",
        "microsoft-edge",
        "microsoft-edge-stable",
    ];
    if let Some(path) = std::env::var_os("PATH") {
        for directory in std::env::split_paths(&path) {
            candidates.extend(PATH_NAMES.iter().map(|name| directory.join(name)));
        }
    }
    candidates
}

fn resolve_publication_browser_executable() -> String {
    let configured = configured_or_env("GENTLE_BROWSER_BIN");
    let installed_candidates = publication_browser_candidates();
    choose_publication_browser(
        (!configured.is_empty()).then_some(configured.as_str()),
        &installed_candidates,
        "chromium",
    )
}

fn print_html_to_pdf(html_path: &Path, pdf_path: &Path) -> Result<(), String> {
    let browser = resolve_publication_browser_executable();
    let canonical_html = fs::canonicalize(html_path).map_err(|error| {
        format!(
            "Could not resolve printable HTML '{}': {error}",
            html_path.display()
        )
    })?;
    let file_url = format!("file://{}", canonical_html.to_string_lossy());
    run(
        Command::new(browser)
            .arg("--headless")
            .arg("--disable-gpu")
            .arg("--no-pdf-header-footer")
            .arg(format!("--print-to-pdf={}", pdf_path.display()))
            .arg(file_url),
        "browser HTML-to-PDF rendering",
        "GENTLE_BROWSER_BIN",
    )?;
    require_generated_file(pdf_path, "browser HTML-to-PDF rendering")
}

/// Generate a canonical multi-gene assay dossier, one HTML page per gene, an
/// index/meta page, a print HTML projection, and an optional browser-rendered
/// PDF. Explicit block selection is presentation-only and cannot add content.
pub fn generate_gene_isoform_assay_publication(
    request_path: &Path,
    output_directory: &Path,
    selected_profile: Option<&str>,
    selected_block_ids: &[String],
    generate_pdf: bool,
) -> Result<GeneIsoformAssayPublicationProjectionReport, String> {
    let request_bytes = fs::read(request_path).map_err(|error| {
        format!(
            "Could not read isoform-assay publication request '{}': {error}",
            request_path.display()
        )
    })?;
    let request: GeneIsoformAssayPublicationRequest = serde_json::from_slice(&request_bytes)
        .map_err(|error| {
            format!(
                "Could not parse isoform-assay publication request '{}': {error}",
                request_path.display()
            )
        })?;
    if request.schema != GENE_ISOFORM_ASSAY_PUBLICATION_REQUEST_SCHEMA {
        return Err(format!("Unsupported schema '{}'", request.schema));
    }
    if request.report_id.trim().is_empty()
        || request.title.trim().is_empty()
        || request.genes.is_empty()
    {
        return Err("report_id, title, and at least one gene are required".to_string());
    }
    for (label, filename) in [
        ("index_filename", request.index_filename.as_str()),
        ("print_filename", request.print_filename.as_str()),
        ("pdf_filename", request.pdf_filename.as_str()),
    ] {
        validate_leaf_filename(filename, label)?;
    }
    if request
        .index_filename
        .eq_ignore_ascii_case(&request.print_filename)
        || request
            .index_filename
            .eq_ignore_ascii_case(&request.pdf_filename)
        || request
            .print_filename
            .eq_ignore_ascii_case(&request.pdf_filename)
    {
        return Err("index, print HTML, and PDF filenames must be distinct".to_string());
    }

    let base = request_path.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_directory)
        .map_err(|error| format!("Could not create '{}': {error}", output_directory.display()))?;
    let figures_dir = output_directory.join("figures");
    let data_dir = output_directory.join("data");
    fs::create_dir_all(&figures_dir).map_err(|error| error.to_string())?;
    fs::create_dir_all(&data_dir).map_err(|error| error.to_string())?;
    let mut targets = BTreeMap::new();
    let mut symbols = BTreeSet::new();
    let mut genes = vec![];
    let mut artifacts = vec![];
    for gene_request in &request.genes {
        let symbol = gene_request.gene_symbol.trim();
        if symbol.is_empty() || !symbols.insert(symbol.to_ascii_uppercase()) {
            return Err(format!(
                "Gene symbols must be non-empty and unique: '{}'",
                gene_request.gene_symbol
            ));
        }
        let study_plan = load_publication_bound_report(
            base,
            &gene_request.study_plan,
            "gentle.gene_isoform_assay_study_plan.v1",
            &["plan_id"],
        )?;
        let plan_gene = study_plan
            .value
            .get("gene_symbol")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        if !plan_gene.eq_ignore_ascii_case(symbol) {
            return Err(format!(
                "Study plan '{}' belongs to gene '{}' rather than '{}'",
                study_plan.report_id, plan_gene, symbol
            ));
        }
        let handoffs = gene_request
            .handoffs
            .iter()
            .map(|reference| {
                load_publication_bound_report(
                    base,
                    reference,
                    "gentle.experimental_assay_handoff.v1",
                    &["package_id"],
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let order_forms = gene_request
            .order_forms
            .iter()
            .map(|reference| {
                load_publication_bound_report(
                    base,
                    reference,
                    "gentle.oligo_order_form.v1",
                    &["form_id"],
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut figures = vec![];
        for figure in &gene_request.figures {
            validate_leaf_filename(&figure.figure_id, "figure_id")?;
            let source = source_path(base, &figure.source_path);
            let (target, name, copied) = copy_asset(&source, &figures_dir, "figure", &mut targets)?;
            if copied {
                artifacts.push(publication_artifact(
                    &target,
                    publication_media_type(&target),
                )?);
            }
            figures.push(GeneSetPublicationFigure {
                figure_id: figure.figure_id.clone(),
                kind: figure.kind.clone(),
                label: figure.label.clone(),
                web_path: format!("figures/{name}"),
                pdf_path: None,
                caption: figure.caption.clone(),
                alt_text: figure.alt_text.clone(),
                include_in_pdf: figure.include_in_pdf,
            });
        }
        let slug = publication_gene_slug(symbol)?;
        let order_sheet_name = format!("{slug}-oligo-order.tsv");
        let order_sheet_path =
            (!order_forms.is_empty()).then(|| format!("data/{order_sheet_name}"));
        let (status, status_reason) = resolve_publication_gene_status(gene_request, &handoffs)?;
        let gene = GeneIsoformAssayPublicationGene {
            gene_symbol: symbol.to_string(),
            page_path: format!("gene-{slug}.html"),
            study_plan,
            handoffs,
            order_forms,
            figures,
            order_sheet_path,
            warnings: vec![],
            status,
            status_reason,
        };
        validate_handoff_plan_binding(&gene)?;
        validate_order_form_handoff_binding(&gene)?;
        if let Some(relative) = &gene.order_sheet_path {
            let path = output_directory.join(relative);
            if write_publication_order_sheet(&path, &gene)? {
                artifacts.push(publication_artifact(&path, "text/tab-separated-values")?);
            }
        }
        genes.push(gene);
    }
    let favicon_path = if let Some(raw) = &request.favicon_source_path {
        let source = source_path(base, raw);
        let target = output_directory.join("favicon.ico");
        fs::copy(&source, &target)
            .map_err(|error| format!("Could not copy favicon '{}': {error}", source.display()))?;
        artifacts.push(publication_artifact(&target, "image/x-icon")?);
        Some("favicon.ico".to_string())
    } else {
        None
    };
    let (common_parameters, parameter_overrides) = publication_parameter_entries(&genes);
    let mut content_blocks = vec![GeneIsoformAssayPublicationBlock {
        block_id: "run.parameters".to_string(),
        label: "Run parameters".to_string(),
        projection: "run_parameters".to_string(),
        source_pointer: "/common_parameters".to_string(),
        gene_symbol: None,
    }];
    let mut full = vec!["run.parameters".to_string()];
    let mut summary = vec!["run.parameters".to_string()];
    let mut ordering = vec!["run.parameters".to_string()];
    for (gene_index, gene) in genes.iter().enumerate() {
        let prefix = gene.gene_symbol.to_ascii_lowercase();
        let definitions = [
            (
                "overview",
                "Study decision",
                "study_overview",
                "study_plan/value",
                true,
                true,
                false,
            ),
            (
                "evidence",
                "Evidence summary",
                "evidence_summary",
                "study_plan/value/evidence_summary",
                true,
                true,
                false,
            ),
            (
                "decision",
                "Decision factors",
                "decision_factors",
                "study_plan/value/decision_factors",
                true,
                false,
                false,
            ),
            (
                "operations",
                "Planned operations",
                "planned_operations",
                "study_plan/value/planned_operations",
                true,
                false,
                false,
            ),
            (
                "assays",
                "Assay readiness",
                "assay_handoffs",
                "handoffs",
                true,
                true,
                true,
            ),
            (
                "quality",
                "Quality assurance",
                "quality_assurance",
                "handoffs",
                true,
                true,
                true,
            ),
            (
                "order",
                "Oligo order sheet",
                "order_sheet",
                "order_forms",
                true,
                false,
                true,
            ),
            (
                "figures", "Figures", "figures", "figures", true, true, false,
            ),
            (
                "provenance",
                "Provenance",
                "provenance",
                "study_plan",
                true,
                true,
                true,
            ),
        ];
        for (suffix, label, projection, pointer, in_full, in_summary, in_ordering) in definitions {
            let block_id = format!("gene.{prefix}.{suffix}");
            content_blocks.push(GeneIsoformAssayPublicationBlock {
                block_id: block_id.clone(),
                label: label.to_string(),
                projection: projection.to_string(),
                source_pointer: format!("/genes/{gene_index}/{pointer}"),
                gene_symbol: Some(gene.gene_symbol.clone()),
            });
            if in_full {
                full.push(block_id.clone());
            }
            if in_summary {
                summary.push(block_id.clone());
            }
            if in_ordering {
                ordering.push(block_id);
            }
        }
    }
    let profiles = vec![
        GeneIsoformAssayPublicationProfile {
            profile_id: "full".to_string(),
            label: "Complete scientific and procurement dossier".to_string(),
            block_ids: full,
        },
        GeneIsoformAssayPublicationProfile {
            profile_id: "review".to_string(),
            label: "Compact biological and assay review".to_string(),
            block_ids: summary,
        },
        GeneIsoformAssayPublicationProfile {
            profile_id: "ordering".to_string(),
            label: "Readiness and procurement review".to_string(),
            block_ids: ordering,
        },
    ];
    let gene_count = genes.len();
    let resolved_gene_count = genes
        .iter()
        .filter(|gene| gene.status.is_resolved())
        .count();
    let pending_gene_count = genes
        .iter()
        .filter(|gene| gene.status == GeneIsoformAssayPublicationGeneStatus::Pending)
        .count();
    let unresolved_gene_count = genes
        .iter()
        .filter(|gene| gene.status == GeneIsoformAssayPublicationGeneStatus::Unresolved)
        .count();
    let report = GeneIsoformAssayPublicationReport {
        schema: GENE_ISOFORM_ASSAY_PUBLICATION_SCHEMA.to_string(),
        report_id: request.report_id.clone(),
        title: request.title.clone(),
        subtitle: request.subtitle.clone(),
        generated_date: request.generated_date.clone(),
        default_profile: request.default_profile.clone(),
        index_path: request.index_filename.clone(),
        print_path: request.print_filename.clone(),
        common_parameters,
        parameter_overrides,
        genes,
        content_blocks,
        profiles,
        footer: request.footer.clone(),
        favicon_path,
        warnings: {
            let mut warnings = vec![
                "GENtle report values remain sequence/evidence/design records; approval authorizes an exact operation payload but does not validate its biological interpretation."
                    .to_string(),
            ];
            if pending_gene_count > 0 || unresolved_gene_count > 0 {
                warnings.push(format!(
                    "This dossier is partial: {resolved_gene_count} of {gene_count} genes carry designed assay results, {pending_gene_count} are pending, and {unresolved_gene_count} could not be addressed with the established automatism. Regenerate it once the outstanding results are in."
                ));
            }
            warnings
        },
        resolved_gene_count,
        pending_gene_count,
        unresolved_gene_count,
        complete: pending_gene_count == 0 && unresolved_gene_count == 0,
    };
    let selected_profile = selected_profile
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .unwrap_or(&report.default_profile)
        .to_string();
    let canonical_path = output_directory.join("canonical-report.json");
    let canonical_bytes = serde_json::to_vec_pretty(&report).map_err(|error| error.to_string())?;
    fs::write(&canonical_path, &canonical_bytes).map_err(|error| error.to_string())?;
    let canonical_report_sha256 = sha256_prefixed_bytes(&canonical_bytes);
    artifacts.push(publication_artifact(&canonical_path, "application/json")?);

    let index_path = output_directory.join(&request.index_filename);
    fs::write(
        &index_path,
        render_gene_isoform_assay_publication_index(&report),
    )
    .map_err(|error| error.to_string())?;
    artifacts.push(publication_artifact(&index_path, "text/html")?);
    for (gene_index, gene) in report.genes.iter().enumerate() {
        let path = output_directory.join(&gene.page_path);
        fs::write(
            &path,
            render_gene_isoform_assay_publication_gene(
                &report,
                gene_index,
                &selected_profile,
                selected_block_ids,
            )?,
        )
        .map_err(|error| error.to_string())?;
        artifacts.push(publication_artifact(&path, "text/html")?);
    }
    let print_path = output_directory.join(&request.print_filename);
    fs::write(
        &print_path,
        render_gene_isoform_assay_publication_print(
            &report,
            &selected_profile,
            selected_block_ids,
        )?,
    )
    .map_err(|error| error.to_string())?;
    artifacts.push(publication_artifact(&print_path, "text/html")?);
    if generate_pdf {
        let pdf_path = output_directory.join(&request.pdf_filename);
        remove_generated_file_if_present(&pdf_path)?;
        print_html_to_pdf(&print_path, &pdf_path)?;
        artifacts.push(publication_artifact(&pdf_path, "application/pdf")?);
    }
    artifacts.sort_by(|left, right| left.path.cmp(&right.path));
    let selected_profile_block_ids = report
        .profiles
        .iter()
        .find(|profile| profile.profile_id == selected_profile)
        .map(|profile| profile.block_ids.clone())
        .ok_or_else(|| format!("Unknown publication profile '{selected_profile}'"))?;
    let receipt = GeneIsoformAssayPublicationProjectionReport {
        schema: GENE_ISOFORM_ASSAY_PUBLICATION_PROJECTION_SCHEMA.to_string(),
        canonical_report_path: canonical_path.to_string_lossy().into_owned(),
        canonical_report_sha256,
        selected_profile,
        selected_block_ids: if selected_block_ids.is_empty() {
            selected_profile_block_ids
        } else {
            selected_block_ids.to_vec()
        },
        artifacts,
        warnings: vec![
            "Block/profile selection is a presentation projection of the immutable canonical report and does not alter scientific content."
                .to_string(),
        ],
    };
    fs::write(
        output_directory.join("projection-report.json"),
        serde_json::to_vec_pretty(&receipt).map_err(|error| error.to_string())?,
    )
    .map_err(|error| error.to_string())?;
    Ok(receipt)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_minimal_primer_table(directory: &Path) {
        fs::write(
            directory.join("primers.tsv"),
            "gene\tpair_alias\tselection_role\tselection_status\tforward_sequence_5_to_3\treverse_sequence_5_to_3\nGENE1\tG1_P1\tprimary\tselected\tACGT\tTGCA\n",
        )
        .expect("primer table");
    }

    fn write_minimal_isoform_publication_request(directory: &Path) -> PathBuf {
        let plan = serde_json::json!({
            "schema": "gentle.gene_isoform_assay_study_plan.v1",
            "plan_id": "plan_gene1",
            "gene_symbol": "GENE1",
            "recommended_profile": "comprehensive_isoform_dossier",
            "selected_profile": "isoform_discrimination",
            "profile_override": {
                "selected_profile": "isoform_discrimination",
                "reason": "Review the compact synthetic panel first."
            },
            "iteration": 2,
            "prior_plan_id": "plan_gene1_initial",
            "observations": [{
                "observation_id": "observation_1",
                "statement": "One synthetic band was observed.",
                "source": "synthetic test notebook",
                "validation_status": "user_supplied_unvalidated"
            }],
            "operation_batch_sha256": "sha256:operations",
            "normalized_request": {
                "policy": {
                    "schema": "gentle.gene_isoform_assay_study_policy.v1",
                    "policy_version": "1",
                    "short_min_amplicon_bp": 70,
                    "short_max_amplicon_bp": 250
                }
            },
            "evidence_summary": {
                "transcript_count": 3,
                "exact_cdna_equivalence_group_count": 2,
                "responsive_region_count": 1
            },
            "decision_factors": [{
                "rule_id": "multiple_cdna_classes",
                "triggered": true,
                "summary": "Two exact mature-cDNA classes require discrimination."
            }],
            "planned_operations": [{
                "step_index": 0,
                "step_id": "short_discrimination",
                "purpose": "Distinguish exact mature-cDNA classes",
                "operation_sha256": "sha256:operation",
                "operation": {
                    "DesignTranscriptAssayPanel": {
                        "seq_id": "gene1_seq",
                        "source_feature_id": 0,
                        "report_id": "panel_gene1"
                    }
                }
            }]
        });
        let plan_bytes = serde_json::to_vec_pretty(&plan).unwrap();
        fs::write(directory.join("plan-minimal.json"), &plan_bytes).unwrap();
        let request = serde_json::json!({
            "schema": GENE_ISOFORM_ASSAY_PUBLICATION_REQUEST_SCHEMA,
            "report_id": "dossier_minimal",
            "title": "Synthetic isoform assay dossier",
            "genes": [{
                "gene_symbol": "GENE1",
                "study_plan": {
                    "path": "plan-minimal.json",
                    "expected_sha256": sha256_prefixed_bytes(&plan_bytes)
                }
            }]
        });
        let request_path = directory.join("request-minimal.json");
        fs::write(&request_path, serde_json::to_vec_pretty(&request).unwrap()).unwrap();
        request_path
    }

    #[test]
    fn publication_browser_resolution_prefers_override_then_installed_candidate() {
        let temp = tempfile::tempdir().expect("tempdir");
        let installed = temp.path().join("Chromium");
        fs::write(&installed, b"synthetic browser marker").expect("browser marker");
        let candidates = vec![temp.path().join("missing"), installed.clone()];

        assert_eq!(
            choose_publication_browser(Some("/configured/browser"), &candidates, "chromium"),
            "/configured/browser"
        );
        assert_eq!(
            choose_publication_browser(None, &candidates, "chromium"),
            installed.to_string_lossy()
        );
        assert_eq!(
            choose_publication_browser(None, &[temp.path().join("missing")], "chromium"),
            "chromium"
        );
    }

    #[test]
    fn one_manifest_drives_html_and_markdown() {
        let temp = tempfile::tempdir().expect("tempdir");
        write_minimal_primer_table(temp.path());
        fs::write(
            temp.path().join("figure.svg"),
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"10\" height=\"10\"></svg>",
        )
        .expect("figure");
        let request = serde_json::json!({
            "schema": GENE_SET_PUBLICATION_REQUEST_SCHEMA, "report_id": "synthetic", "title": "Synthetic", "primer_table": {"path": "primers.tsv"},
            "genes": [{"gene_symbol": "GENE1", "figures": [{"figure_id": "gene1-map", "source_path": "figure.svg", "label": "Map"}]}]
        });
        let request_bytes = serde_json::to_vec_pretty(&request).unwrap();
        fs::write(temp.path().join("request.json"), &request_bytes).expect("request");
        let output = temp.path().join("out");
        let generation =
            generate_gene_set_publication(&temp.path().join("request.json"), &output, false)
                .expect("generate");
        let html = fs::read_to_string(output.join("index.html")).unwrap();
        let markdown = fs::read_to_string(output.join("report.md")).unwrap();
        assert!(html.contains("G1_P1") && markdown.contains("G1_P1"));
        assert!(html.contains("figures/figure.svg") && markdown.contains("figures/figure.svg"));
        assert!(html.contains("1 primary pairs"));
        let manifest_bytes = fs::read(output.join("bundle-manifest.json")).unwrap();
        let manifest: GeneSetPublicationBundleManifest =
            serde_json::from_slice(&manifest_bytes).unwrap();
        assert_eq!(manifest.schema, GENE_SET_PUBLICATION_BUNDLE_MANIFEST_SCHEMA);
        assert_eq!(
            manifest.request_sha256,
            sha256_prefixed_bytes(&request_bytes)
        );
        assert_eq!(
            generation.bundle_manifest_sha256,
            sha256_prefixed_bytes(&manifest_bytes)
        );
        assert_eq!(generation.artifacts, manifest.artifacts);
        for expected in [
            "data/primers.tsv",
            "figures/figure.svg",
            "index.html",
            "report.md",
            "resolved-report.json",
        ] {
            let artifact = manifest
                .artifacts
                .iter()
                .find(|artifact| artifact.path == expected)
                .unwrap_or_else(|| panic!("missing content-bound artifact {expected}"));
            assert_eq!(
                artifact.sha256,
                sha256_prefixed_bytes(&fs::read(output.join(expected)).unwrap())
            );
        }
        assert!(
            !manifest
                .artifacts
                .iter()
                .any(|artifact| artifact.path == "bundle-manifest.json"
                    || artifact.path == "generation-report.json")
        );
    }

    #[test]
    fn manifest_output_filenames_cannot_escape_or_replace_bundle_internals() {
        let mut request = GeneSetPublicationRequest {
            report_id: "synthetic".to_string(),
            title: "Synthetic".to_string(),
            genes: vec![gentle_protocol::GeneSetPublicationGeneRequest {
                gene_symbol: "GENE1".to_string(),
                ..Default::default()
            }],
            ..Default::default()
        };

        request.html_filename = "../outside.html".to_string();
        assert!(
            validate_request(&request)
                .expect_err("parent traversal must be rejected")
                .contains("without directory components")
        );

        request.html_filename = "resolved-report.json".to_string();
        assert!(
            validate_request(&request)
                .expect_err("internal report path must be reserved")
                .contains("reserved bundle path")
        );

        request.html_filename = "report.md".to_string();
        assert!(
            validate_request(&request)
                .expect_err("root outputs must not collide")
                .contains("distinct filenames")
        );
    }

    #[test]
    fn duplicate_asset_basenames_are_rejected_instead_of_overwritten() {
        let temp = tempfile::tempdir().expect("tempdir");
        write_minimal_primer_table(temp.path());
        for (directory, fill) in [("one", "red"), ("two", "blue")] {
            let directory = temp.path().join(directory);
            fs::create_dir(&directory).expect("asset directory");
            fs::write(
                directory.join("shared.svg"),
                format!("<svg xmlns=\"http://www.w3.org/2000/svg\"><rect fill=\"{fill}\"/></svg>"),
            )
            .expect("figure");
        }
        let request = serde_json::json!({
            "schema": GENE_SET_PUBLICATION_REQUEST_SCHEMA,
            "report_id": "synthetic",
            "title": "Synthetic",
            "primer_table": {"path": "primers.tsv"},
            "genes": [{
                "gene_symbol": "GENE1",
                "figures": [
                    {"figure_id": "one", "source_path": "one/shared.svg"},
                    {"figure_id": "two", "source_path": "two/shared.svg"}
                ]
            }]
        });
        let request_path = temp.path().join("request.json");
        fs::write(&request_path, serde_json::to_vec_pretty(&request).unwrap()).expect("request");

        let error = generate_gene_set_publication(&request_path, &temp.path().join("out"), false)
            .expect_err("colliding bundle targets must fail");
        assert!(error.contains("Bundle path collision"), "{error}");
        assert!(error.contains("shared.svg"), "{error}");
    }

    #[test]
    fn missing_pdf_tool_reports_override_and_leaves_web_bundle_without_pdf_link() {
        let temp = tempfile::tempdir().expect("tempdir");
        write_minimal_primer_table(temp.path());
        let request = serde_json::json!({
            "schema": GENE_SET_PUBLICATION_REQUEST_SCHEMA,
            "report_id": "synthetic",
            "title": "Synthetic",
            "primer_table": {"path": "primers.tsv"},
            "genes": [{"gene_symbol": "GENE1"}]
        });
        let request_path = temp.path().join("request.json");
        fs::write(&request_path, serde_json::to_vec_pretty(&request).unwrap()).expect("request");
        let missing_pandoc = temp.path().join("missing-pandoc");
        let _override = crate::tool_overrides::ScopedToolOverrideGuard::set(
            "GENTLE_PANDOC_BIN",
            &missing_pandoc.to_string_lossy(),
        );
        let output = temp.path().join("out");

        let error = generate_gene_set_publication(&request_path, &output, true)
            .expect_err("missing pandoc must fail clearly");
        assert!(error.contains("pandoc PDF generation"), "{error}");
        assert!(error.contains("GENTLE_PANDOC_BIN"), "{error}");
        assert!(
            error.contains(&missing_pandoc.to_string_lossy().to_string()),
            "{error}"
        );
        let html = fs::read_to_string(output.join("index.html")).expect("fallback web report");
        assert!(!html.contains("Printable companion"));
        let resolved: GeneSetPublicationReport = serde_json::from_slice(
            &fs::read(output.join("resolved-report.json")).expect("resolved report"),
        )
        .expect("parse resolved report");
        assert_eq!(resolved.pdf_path, None);
        let manifest: GeneSetPublicationBundleManifest = serde_json::from_slice(
            &fs::read(output.join("bundle-manifest.json")).expect("web bundle manifest"),
        )
        .expect("parse web bundle manifest");
        assert!(
            manifest
                .artifacts
                .iter()
                .any(|artifact| artifact.path == "index.html")
        );
        assert!(
            !manifest
                .artifacts
                .iter()
                .any(|artifact| artifact.path.ends_with(".pdf"))
        );
    }

    #[test]
    fn isoform_assay_missing_browser_fails_closed_after_writing_html() {
        let temp = tempfile::tempdir().expect("tempdir");
        let request_path = write_minimal_isoform_publication_request(temp.path());
        let output = temp.path().join("out");
        let missing_browser = temp.path().join("missing-browser");
        let _override = crate::tool_overrides::ScopedToolOverrideGuard::set(
            "GENTLE_BROWSER_BIN",
            &missing_browser.to_string_lossy(),
        );

        let error = generate_gene_isoform_assay_publication(
            &request_path,
            &output,
            Some("full"),
            &[],
            true,
        )
        .expect_err("missing browser must fail clearly");
        assert!(error.contains("browser HTML-to-PDF rendering"), "{error}");
        assert!(error.contains("GENTLE_BROWSER_BIN"), "{error}");
        assert!(
            error.contains(&missing_browser.to_string_lossy().to_string()),
            "{error}"
        );
        for expected in [
            "canonical-report.json",
            "index.html",
            "gene-GENE1.html",
            "print.html",
        ] {
            assert!(output.join(expected).is_file(), "missing {expected}");
        }
        assert!(!output.join("report.pdf").exists());
        assert!(!output.join("projection-report.json").exists());
    }

    /// Rewrite the single gene of the minimal request so a test can declare a
    /// status without rebuilding the whole fixture.
    fn set_minimal_request_gene_status(
        request_path: &Path,
        status: &str,
        status_reason: Option<&str>,
    ) {
        let mut request: serde_json::Value =
            serde_json::from_slice(&fs::read(request_path).expect("read request"))
                .expect("parse request");
        request["genes"][0]["status"] = serde_json::Value::String(status.to_string());
        if let Some(reason) = status_reason {
            request["genes"][0]["status_reason"] = serde_json::Value::String(reason.to_string());
        }
        fs::write(
            request_path,
            serde_json::to_vec_pretty(&request).expect("serialize request"),
        )
        .expect("write request");
    }

    #[test]
    fn isoform_assay_gene_without_results_is_published_as_pending() {
        let temp = tempfile::tempdir().expect("tempdir");
        let request_path = write_minimal_isoform_publication_request(temp.path());
        let output = temp.path().join("out");
        generate_gene_isoform_assay_publication(&request_path, &output, Some("full"), &[], false)
            .expect("publish a dossier whose only gene has no results yet");

        let canonical: serde_json::Value =
            serde_json::from_slice(&fs::read(output.join("canonical-report.json")).unwrap())
                .expect("parse canonical report");
        assert_eq!(canonical["genes"][0]["status"], "pending");
        assert_eq!(canonical["complete"], false);
        assert_eq!(canonical["resolved_gene_count"], 0);
        assert_eq!(canonical["pending_gene_count"], 1);
        assert_eq!(canonical["unresolved_gene_count"], 0);

        let index = fs::read_to_string(output.join("index.html")).expect("index page");
        assert!(index.contains("Partial dossier"), "{index}");
        assert!(index.contains("Pending"), "{index}");
        assert!(index.contains("0 of 1 genes carry designed assay results"));
        let gene_page = fs::read_to_string(output.join("gene-GENE1.html")).expect("gene page");
        assert!(gene_page.contains("Pending"), "{gene_page}");
        assert!(
            gene_page.contains("has not produced results yet"),
            "{gene_page}"
        );
    }

    #[test]
    fn isoform_assay_unresolved_gene_states_why_the_automatism_could_not_address_it() {
        let temp = tempfile::tempdir().expect("tempdir");
        let request_path = write_minimal_isoform_publication_request(temp.path());
        set_minimal_request_gene_status(&request_path, "unresolved", None);
        let output = temp.path().join("out");
        let error = generate_gene_isoform_assay_publication(
            &request_path,
            &output,
            Some("full"),
            &[],
            false,
        )
        .expect_err("an unresolved gene without a reason must be rejected");
        assert!(error.contains("requires a status_reason"), "{error}");

        set_minimal_request_gene_status(
            &request_path,
            "unresolved",
            Some("Primer3 found no admissible bounded region for the junction target."),
        );
        let output = temp.path().join("out-with-reason");
        generate_gene_isoform_assay_publication(&request_path, &output, Some("full"), &[], false)
            .expect("publish an unresolved gene that states its reason");

        let canonical: serde_json::Value =
            serde_json::from_slice(&fs::read(output.join("canonical-report.json")).unwrap())
                .expect("parse canonical report");
        assert_eq!(canonical["genes"][0]["status"], "unresolved");
        assert_eq!(canonical["unresolved_gene_count"], 1);
        assert_eq!(canonical["complete"], false);

        for page in ["index.html", "gene-GENE1.html", "print.html"] {
            let html = fs::read_to_string(output.join(page)).expect(page);
            assert!(
                html.contains("no admissible bounded region"),
                "{page} must state the reason: {html}"
            );
        }
        let gene_page = fs::read_to_string(output.join("gene-GENE1.html")).unwrap();
        assert!(
            gene_page.contains("Not addressed automatically"),
            "{gene_page}"
        );
    }

    #[test]
    fn isoform_assay_gene_declared_resolved_without_handoffs_is_rejected() {
        let temp = tempfile::tempdir().expect("tempdir");
        let request_path = write_minimal_isoform_publication_request(temp.path());
        set_minimal_request_gene_status(&request_path, "resolved", None);
        let output = temp.path().join("out");
        let error = generate_gene_isoform_assay_publication(
            &request_path,
            &output,
            Some("full"),
            &[],
            false,
        )
        .expect_err("a gene without results must not claim to be resolved");
        assert!(
            error.contains("carries no experimental assay handoff"),
            "{error}"
        );
        assert!(error.contains("publishing it as complete"), "{error}");
    }

    #[test]
    fn isoform_assay_profiles_and_explicit_blocks_filter_only_presentation() {
        let temp = tempfile::tempdir().expect("tempdir");
        let request_path = write_minimal_isoform_publication_request(temp.path());
        let review_output = temp.path().join("review");
        let review = generate_gene_isoform_assay_publication(
            &request_path,
            &review_output,
            Some("review"),
            &[],
            false,
        )
        .expect("render review profile");
        let review_html =
            fs::read_to_string(review_output.join("gene-GENE1.html")).expect("review page");
        assert!(review_html.contains("Study decision"));
        assert!(review_html.contains("Evidence summary"));
        assert!(!review_html.contains("Approved operation candidates"));
        assert!(review_html.contains("Provenance"));
        assert!(review_html.contains("GENtle study plan"));
        assert!(
            review
                .selected_block_ids
                .contains(&"gene.gene1.overview".to_string())
        );
        assert!(
            !review
                .selected_block_ids
                .contains(&"gene.gene1.operations".to_string())
        );
        assert!(
            review
                .selected_block_ids
                .contains(&"gene.gene1.provenance".to_string())
        );

        let operations_output = temp.path().join("operations-only");
        let requested_blocks = vec!["gene.gene1.operations".to_string()];
        let operations = generate_gene_isoform_assay_publication(
            &request_path,
            &operations_output,
            Some("full"),
            &requested_blocks,
            false,
        )
        .expect("render declared operation block");
        let operations_html =
            fs::read_to_string(operations_output.join("gene-GENE1.html")).expect("operations page");
        assert!(operations_html.contains("Approved operation candidates"));
        assert!(operations_html.contains("Distinguish exact mature-cDNA classes"));
        assert!(!operations_html.contains("Study decision"));
        assert!(!operations_html.contains("Evidence summary"));
        assert_eq!(operations.selected_block_ids, requested_blocks);
        assert_eq!(
            operations.canonical_report_sha256, review.canonical_report_sha256,
            "presentation selection must not change canonical scientific content"
        );
    }

    #[test]
    fn isoform_assay_publication_uses_bound_reports_and_declared_blocks() {
        let temp = tempfile::tempdir().expect("tempdir");
        let plan = serde_json::json!({
            "schema": "gentle.gene_isoform_assay_study_plan.v1",
            "plan_id": "plan_gene1",
            "gene_symbol": "GENE1",
            "recommended_profile": "isoform_discrimination",
            "selected_profile": "isoform_discrimination",
            "iteration": 1,
            "operation_batch_sha256": "sha256:operations",
            "normalized_request": {
                "policy": {
                    "schema": "gentle.gene_isoform_assay_study_policy.v1",
                    "policy_version": "1",
                    "short_min_amplicon_bp": 70,
                    "short_max_amplicon_bp": 250
                }
            },
            "evidence_summary": {
                "transcript_count": 3,
                "exact_cdna_equivalence_group_count": 2,
                "responsive_region_count": 1
            },
            "decision_factors": [{
                "rule_id": "multiple_cdna_classes",
                "triggered": true,
                "summary": "Two exact mature-cDNA classes require discrimination."
            }],
            "planned_operations": [{
                "step_index": 0,
                "step_id": "short_discrimination",
                "purpose": "Distinguish exact mature-cDNA classes",
                "operation_sha256": "sha256:operation",
                "operation": {
                    "DesignTranscriptAssayPanel": {
                        "seq_id": "gene1_seq",
                        "source_feature_id": 0,
                        "report_id": "panel_gene1"
                    }
                }
            }]
        });
        let plan_bytes = serde_json::to_vec_pretty(&plan).unwrap();
        fs::write(temp.path().join("plan.json"), &plan_bytes).unwrap();
        let readiness_policy = ExperimentalAssayReadinessPolicy::default();
        let readiness_policy_text = serde_json::to_string(&readiness_policy).unwrap();
        let readiness_policy_sha256 = sha256_prefixed_bytes(readiness_policy_text.as_bytes());
        let readiness_policy_id = format!(
            "readiness_policy_sha256_{}",
            readiness_policy_sha256.trim_start_matches("sha256:")
        );
        let readiness_row = serde_json::json!({
            "card_id": "card_1",
            "assay_id": "assay_1",
            "pair_id": "pair_1",
            "pair_rank": 1,
            "display_label": "GENE1 pair 1",
            "readiness_state": "order_ready",
            "order_ready": true,
            "oligo_ids": ["oligo_1"],
            "sequences_5_to_3": ["ACGTACGT"]
        });
        let handoff = serde_json::json!({
            "schema": "gentle.experimental_assay_handoff.v1",
            "package_id": "handoff_gene1",
            "source_panel_report_id": "panel_gene1",
            "source_panel_schema": "gentle.transcript_assay_panel.v2",
            "source_panel_sha256": "sha256:panel_gene1",
            "source_seq_id": "gene1_seq",
            "source_feature_id": 0,
            "policy": readiness_policy,
            "policy_id": readiness_policy_id,
            "order_readiness_table": [readiness_row.clone()],
            "cards": [{
                "schema": "gentle.experimental_assay_card.v1",
                "card_id": "card_1",
                "panel_report_id": "panel_gene1",
                "assay_id": "assay_1",
                "pair_id": "pair_1",
                "pair_rank": 1,
                "display_label": "GENE1 pair 1",
                "readiness_state": "order_ready",
                "policy_schema": "gentle.experimental_assay_readiness_policy.v1",
                "policy_version": "1",
                "policy_id": readiness_policy_id,
                "gate_outcomes": [
                    {
                        "gate": "critical_oligo_qc",
                        "required": true,
                        "status": "pass",
                        "summary": "GENtle critical oligo QC passed.",
                        "evidence_ids": ["panel_gene1"]
                    },
                    {
                        "gate": "genomic_carryover",
                        "required": true,
                        "status": "pass",
                        "summary": "The prepared genomic-DNA search passed.",
                        "evidence_ids": ["primer_specificity_genomic_1"]
                    },
                    {
                        "gate": "transcriptome_specificity",
                        "required": true,
                        "status": "pass",
                        "summary": "The prepared whole-transcriptome cDNA search passed.",
                        "evidence_ids": ["primer_specificity_cdna_1"]
                    },
                    {
                        "gate": "variant_evidence",
                        "required": false,
                        "status": "not_evaluated",
                        "summary": "No provenance-bound variant evidence was supplied.",
                        "evidence_ids": []
                    }
                ],
                "blockers": [],
                "variant_evidence_status": "not_evaluated",
                "warnings": ["Synthetic QA fixture warning."]
            }]
        });
        let handoff_bytes = serde_json::to_vec_pretty(&handoff).unwrap();
        let handoff_sha256 = sha256_prefixed_bytes(&handoff_bytes);
        fs::write(temp.path().join("handoff.json"), &handoff_bytes).unwrap();
        let order = serde_json::json!({
            "schema": "gentle.oligo_order_form.v1",
            "form_id": "order_gene1",
            "line_items": [{
                "line_id": "line_1",
                "line_no": 1,
                "name": "GENE1_F",
                "role": "forward",
                "sequence_5_to_3": "ACGTACGT",
                "length_nt": 8,
                "modifications": [],
                "scale": "25 nmol",
                "purification": "desalt",
                "provenance": {
                    "source_kind": "experimental_assay_handoff",
                    "report_id": "handoff_gene1",
                    "report_schema": "gentle.experimental_assay_handoff.v1",
                    "report_sha256": handoff_sha256,
                    "readiness_policy_id": readiness_policy_id,
                    "readiness_policy_sha256": readiness_policy_sha256,
                    "readiness_card_id": "card_1",
                    "readiness_state": "order_ready",
                    "assay_id": "assay_1",
                    "pair_id": "pair_1",
                    "oligo_id": "oligo_1",
                    "readiness_row": readiness_row
                }
            }]
        });
        let order_bytes = serde_json::to_vec_pretty(&order).unwrap();
        fs::write(temp.path().join("order.json"), &order_bytes).unwrap();
        let request = serde_json::json!({
            "schema": GENE_ISOFORM_ASSAY_PUBLICATION_REQUEST_SCHEMA,
            "report_id": "dossier_1",
            "title": "Synthetic isoform assay dossier",
            "genes": [{
                "gene_symbol": "GENE1",
                "study_plan": {
                    "path": "plan.json",
                    "expected_sha256": sha256_prefixed_bytes(&plan_bytes)
                },
                "handoffs": [{
                    "path": "handoff.json",
                    "expected_sha256": sha256_prefixed_bytes(&handoff_bytes)
                }],
                "order_forms": [{
                    "path": "order.json",
                    "expected_sha256": sha256_prefixed_bytes(&order_bytes)
                }]
            }]
        });
        let request_path = temp.path().join("request.json");
        fs::write(&request_path, serde_json::to_vec_pretty(&request).unwrap()).unwrap();
        let output = temp.path().join("out");
        let receipt = generate_gene_isoform_assay_publication(
            &request_path,
            &output,
            Some("review"),
            &[],
            false,
        )
        .expect("generate canonical dossier");
        assert_eq!(
            receipt.schema,
            GENE_ISOFORM_ASSAY_PUBLICATION_PROJECTION_SCHEMA
        );
        assert!(output.join("index.html").is_file());
        assert!(output.join("gene-GENE1.html").is_file());
        assert!(output.join("print.html").is_file());
        assert!(output.join("data/GENE1-oligo-order.tsv").is_file());
        let canonical: GeneIsoformAssayPublicationReport =
            serde_json::from_slice(&fs::read(output.join("canonical-report.json")).unwrap())
                .unwrap();
        assert_eq!(
            canonical.genes[0].study_plan.sha256,
            sha256_prefixed_bytes(&plan_bytes)
        );
        assert!(
            canonical
                .content_blocks
                .iter()
                .any(|block| block.block_id == "gene.gene1.order")
        );
        assert!(
            canonical
                .content_blocks
                .iter()
                .any(|block| block.block_id == "gene.gene1.quality"
                    && block.projection == "quality_assurance")
        );
        assert!(canonical.profiles.iter().all(|profile| {
            profile
                .block_ids
                .contains(&"gene.gene1.quality".to_string())
                && profile
                    .block_ids
                    .contains(&"gene.gene1.provenance".to_string())
        }));
        let gene_html = fs::read_to_string(output.join("gene-GENE1.html")).unwrap();
        assert!(gene_html.contains("isoform_discrimination"));
        assert!(gene_html.contains("assay_1"));
        assert!(gene_html.contains("Quality assurance"));
        assert!(gene_html.contains("data-gentle-qa-gate=\"transcriptome_specificity\""));
        assert!(gene_html.contains("primer_specificity_cdna_1"));
        assert!(gene_html.contains("not-evaluated gate is not a pass"));
        assert!(gene_html.contains("data-gentle-provenance-kind=\"study_plan\""));
        assert!(gene_html.contains("data-gentle-provenance-kind=\"experimental_handoff\""));
        assert!(gene_html.contains("data-gentle-provenance-kind=\"oligo_order_form\""));
        assert!(gene_html.contains("Selects study depth from declared evidence"));
        assert!(gene_html.contains("Evaluates the selected assay panel"));
        assert!(gene_html.contains("Projects readiness-qualified assay oligos"));
        assert!(gene_html.contains(&handoff_sha256));
        assert!(gene_html.contains(&readiness_policy_sha256));
        assert!(!gene_html.contains("Sequence 5′→3′"));
        let unknown = render_gene_isoform_assay_publication_gene(
            &canonical,
            0,
            "review",
            &["gene.gene1.fabricated".to_string()],
        )
        .expect_err("undeclared presentation block must fail");
        assert!(unknown.contains("Unknown publication block"), "{unknown}");

        let mut tampered_order = order.clone();
        tampered_order["line_items"][0]["provenance"]["readiness_row"]["order_ready"] =
            serde_json::Value::Bool(false);
        let tampered_order_bytes = serde_json::to_vec_pretty(&tampered_order).unwrap();
        fs::write(
            temp.path().join("order-tampered.json"),
            &tampered_order_bytes,
        )
        .unwrap();
        let mut tampered_request = request.clone();
        tampered_request["genes"][0]["order_forms"][0]["path"] =
            serde_json::Value::String("order-tampered.json".to_string());
        tampered_request["genes"][0]["order_forms"][0]["expected_sha256"] =
            serde_json::Value::String(sha256_prefixed_bytes(&tampered_order_bytes));
        let tampered_request_path = temp.path().join("request-tampered.json");
        fs::write(
            &tampered_request_path,
            serde_json::to_vec_pretty(&tampered_request).unwrap(),
        )
        .unwrap();
        let error = generate_gene_isoform_assay_publication(
            &tampered_request_path,
            &temp.path().join("out-tampered"),
            None,
            &[],
            false,
        )
        .expect_err("tampered readiness row must fail");
        assert!(error.contains("readiness row does not match"), "{error}");

        let mut unrelated_handoff = handoff;
        unrelated_handoff["source_panel_report_id"] =
            serde_json::Value::String("unapproved_panel".to_string());
        let unrelated_handoff_bytes = serde_json::to_vec_pretty(&unrelated_handoff).unwrap();
        fs::write(
            temp.path().join("handoff-unapproved.json"),
            &unrelated_handoff_bytes,
        )
        .unwrap();
        let mut unrelated_request = request;
        unrelated_request["genes"][0]["handoffs"][0]["path"] =
            serde_json::Value::String("handoff-unapproved.json".to_string());
        unrelated_request["genes"][0]["handoffs"][0]["expected_sha256"] =
            serde_json::Value::String(sha256_prefixed_bytes(&unrelated_handoff_bytes));
        let unrelated_request_path = temp.path().join("request-unapproved.json");
        fs::write(
            &unrelated_request_path,
            serde_json::to_vec_pretty(&unrelated_request).unwrap(),
        )
        .unwrap();
        let error = generate_gene_isoform_assay_publication(
            &unrelated_request_path,
            &temp.path().join("out-unapproved"),
            None,
            &[],
            false,
        )
        .expect_err("handoff outside approved operation batch must fail");
        assert!(error.contains("absent from approved study plan"), "{error}");
    }
}
