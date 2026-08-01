//! Resolve one portable gene-set manifest into mutually consistent web and print artifacts.

use gentle_protocol::{
    GENE_SET_PUBLICATION_GENERATION_SCHEMA, GENE_SET_PUBLICATION_REPORT_SCHEMA,
    GENE_SET_PUBLICATION_REQUEST_SCHEMA, GeneSetPublicationDownload, GeneSetPublicationFigure,
    GeneSetPublicationGene, GeneSetPublicationGenerationReport, GeneSetPublicationPrimerColumnMap,
    GeneSetPublicationPrimerRow, GeneSetPublicationReport, GeneSetPublicationRequest,
};
use gentle_render::{render_gene_set_publication_html, render_gene_set_publication_markdown};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::ErrorKind,
    path::{Component, Path, PathBuf},
    process::Command,
};

use crate::{
    svg_pdf::render_svg_file_to_pdf, svg_png::SvgPngRenderOptions,
    tool_overrides::resolve_tool_executable,
};

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
        warnings,
    };
    fs::write(
        output_directory.join("generation-report.json"),
        serde_json::to_string_pretty(&generation).map_err(|error| error.to_string())?,
    )
    .map_err(|error| error.to_string())?;
    Ok(generation)
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
        fs::write(
            temp.path().join("request.json"),
            serde_json::to_vec_pretty(&request).unwrap(),
        )
        .expect("request");
        let output = temp.path().join("out");
        generate_gene_set_publication(&temp.path().join("request.json"), &output, false)
            .expect("generate");
        let html = fs::read_to_string(output.join("index.html")).unwrap();
        let markdown = fs::read_to_string(output.join("report.md")).unwrap();
        assert!(html.contains("G1_P1") && markdown.contains("G1_P1"));
        assert!(html.contains("figures/figure.svg") && markdown.contains("figures/figure.svg"));
        assert!(html.contains("1 primary pairs"));
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
    }
}
