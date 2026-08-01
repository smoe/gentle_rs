//! Resolve one portable gene-set manifest into mutually consistent web and print artifacts.

use gentle_protocol::{
    GENE_SET_PUBLICATION_GENERATION_SCHEMA, GENE_SET_PUBLICATION_REPORT_SCHEMA,
    GENE_SET_PUBLICATION_REQUEST_SCHEMA, GeneSetPublicationDownload, GeneSetPublicationFigure,
    GeneSetPublicationGene, GeneSetPublicationGenerationReport, GeneSetPublicationPrimerColumnMap,
    GeneSetPublicationPrimerRow, GeneSetPublicationReport, GeneSetPublicationRequest,
};
use gentle_render::{render_gene_set_publication_html, render_gene_set_publication_markdown};
use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
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

fn copy_asset(source: &Path, target_dir: &Path, what: &str) -> Result<(PathBuf, String), String> {
    if !source.is_file() {
        return Err(format!("{what} does not exist: '{}'", source.display()));
    }
    let filename = safe_filename(source, what)?;
    let target = target_dir.join(&filename);
    fs::copy(source, &target).map_err(|error| {
        format!(
            "Could not copy {what} '{}' to '{}': {error}",
            source.display(),
            target.display()
        )
    })?;
    Ok((target, filename))
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

fn run(command: &mut Command, description: &str) -> Result<(), String> {
    let output = command
        .output()
        .map_err(|error| format!("Could not run {description}: {error}"))?;
    if output.status.success() {
        return Ok(());
    }
    Err(format!(
        "{description} failed: {}",
        String::from_utf8_lossy(&output.stderr).trim()
    ))
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
    let mut genes = BTreeSet::new();
    for gene in &request.genes {
        if gene.gene_symbol.trim().is_empty() || !genes.insert(gene.gene_symbol.clone()) {
            return Err(format!(
                "Gene symbols must be non-empty and unique: '{}'",
                gene.gene_symbol
            ));
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

    let primer_source = source_path(base, &request.primer_table.path);
    let primers = read_primers(
        &primer_source,
        delimiter(&request.primer_table.delimiter)?,
        &request.primer_table.columns,
    )?;
    let (_, primer_name) = copy_asset(&primer_source, &data_dir, "primer table")?;
    let mut copied_files = vec![format!("data/{primer_name}")];
    let mut pdf_pages = Vec::new();
    let mut figure_ids = BTreeSet::new();
    let mut genes = Vec::new();
    for gene in &request.genes {
        let mut figures = Vec::new();
        for figure in &gene.figures {
            if figure.figure_id.trim().is_empty() || !figure_ids.insert(figure.figure_id.clone()) {
                return Err(format!(
                    "Figure IDs must be non-empty and unique: '{}'",
                    figure.figure_id
                ));
            }
            let source = source_path(base, &figure.source_path);
            let (copied, name) = copy_asset(&source, &figures_dir, "figure")?;
            copied_files.push(format!("figures/{name}"));
            let pdf_path = if generate_pdf && figure.include_in_pdf {
                let pdf_target = if let Some(pdf_source) = &figure.pdf_source_path {
                    let source = source_path(base, pdf_source);
                    let (target, name) = copy_asset(&source, &figures_dir, "figure PDF")?;
                    copied_files.push(format!("figures/{name}"));
                    target
                } else {
                    let target = figures_dir.join(format!("{}.pdf", figure.figure_id));
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
        let (_, name) = copy_asset(&source, &data_dir, "download")?;
        copied_files.push(format!("data/{name}"));
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
        pdf_path: generate_pdf.then(|| request.pdf_filename.clone()),
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
    let pdf_path = if generate_pdf {
        let narrative = output_directory.join("narrative.pdf");
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
        )?;
        let assembled_pdf = output_directory.join("assembled-uncompressed.pdf");
        let final_pdf = output_directory.join(&request.pdf_filename);
        let unite = resolve_tool_executable("GENTLE_PDFUNITE_BIN", "pdfunite");
        let mut command = Command::new(unite);
        command.arg(&narrative);
        for page in &pdf_pages {
            command.arg(page);
        }
        command.arg(&assembled_pdf);
        run(&mut command, "PDF page assembly")?;
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
                .arg(format!("-sOutputFile={}", final_pdf.display()))
                .arg(&assembled_pdf),
            "Ghostscript PDF optimization",
        )?;
        fs::remove_file(narrative).map_err(|error| error.to_string())?;
        fs::remove_file(assembled_pdf).map_err(|error| error.to_string())?;
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
        warnings: vec![],
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

    #[test]
    fn one_manifest_drives_html_and_markdown() {
        let temp = tempfile::tempdir().expect("tempdir");
        fs::write(temp.path().join("primers.tsv"), "gene\tpair_alias\tselection_role\tselection_status\tforward_sequence_5_to_3\treverse_sequence_5_to_3\nGENE1\tG1_P1\tprimary\tselected\tACGT\tTGCA\n").expect("primer table");
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
    }
}
