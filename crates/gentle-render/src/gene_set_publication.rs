//! Responsive HTML and printable Markdown for multi-gene publication reports.
//!
//! Both formats consume the same resolved protocol record. This prevents a web
//! landing page, primer table, and archival PDF narrative from drifting into
//! contradictory analysis summaries.

use gentle_protocol::{
    GeneSetPublicationGene, GeneSetPublicationPrimerRow, GeneSetPublicationReport,
};

fn html_escape(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

fn markdown_cell(value: &str) -> String {
    value.replace('|', "\\|").replace('\n', " ")
}

fn status_class(status: &str) -> &'static str {
    let normalized = status.to_ascii_lowercase();
    if normalized.contains("reject") || normalized.contains("fail") {
        "danger"
    } else if normalized.contains("review")
        || normalized.contains("gate")
        || normalized.contains("warning")
    {
        "warning"
    } else {
        "ready"
    }
}

fn primer_rows_for_gene<'a>(
    report: &'a GeneSetPublicationReport,
    gene: &GeneSetPublicationGene,
) -> Vec<&'a GeneSetPublicationPrimerRow> {
    report
        .primers
        .iter()
        .filter(|row| row.gene.eq_ignore_ascii_case(&gene.gene_symbol))
        .collect()
}

fn render_html_primer_table(rows: &[&GeneSetPublicationPrimerRow]) -> String {
    if rows.is_empty() {
        return "<p class=\"muted\">No primer rows were supplied for this gene.</p>".to_string();
    }
    let mut html = String::from(
        "<div class=\"table-scroll\"><table><thead><tr><th>Pair</th><th>Role / status</th><th>Forward primer 5′→3′</th><th>Reverse primer 5′→3′</th><th>Tm, °C</th><th>Specificity / genome</th></tr></thead><tbody>",
    );
    for row in rows {
        let note = if row.note.is_empty() {
            String::new()
        } else {
            format!("<br><small>{}</small>", html_escape(&row.note))
        };
        html.push_str(&format!(
            "<tr data-gentle-primer-pair=\"{}\"><td><b>{}</b></td><td>{}<br><span class=\"pill {}\">{}</span></td><td><code>{}</code></td><td><code>{}</code></td><td>{} / {}</td><td>{} / {}{}</td></tr>",
            html_escape(&row.pair_id),
            html_escape(&row.pair_id),
            html_escape(&row.role),
            status_class(&row.status),
            html_escape(&row.status),
            html_escape(&row.forward_sequence_5_to_3),
            html_escape(&row.reverse_sequence_5_to_3),
            html_escape(&row.forward_tm_c),
            html_escape(&row.reverse_tm_c),
            html_escape(&row.cdna_specificity),
            html_escape(&row.genome_assessment),
            note,
        ));
    }
    html.push_str("</tbody></table></div>");
    html
}

/// Render one responsive, print-aware web page from the resolved report.
pub fn render_gene_set_publication_html(report: &GeneSetPublicationReport) -> String {
    let favicon = report
        .favicon_path
        .as_ref()
        .map(|path| {
            format!(
                "<link rel=\"icon\" href=\"{}\" sizes=\"any\">",
                html_escape(path)
            )
        })
        .unwrap_or_default();
    let conditions = if report.conditions.is_empty() {
        String::new()
    } else {
        format!(
            "<p class=\"conditions\"><b>Conditions:</b> {}</p>",
            report
                .conditions
                .iter()
                .map(|value| format!("<code>{}</code>", html_escape(value)))
                .collect::<Vec<_>>()
                .join(" · ")
        )
    };
    let metrics = report
        .metrics
        .iter()
        .map(|metric| {
            format!(
                "<div class=\"metric\"><b>{}</b><strong>{}</strong><small>{}</small></div>",
                html_escape(&metric.value),
                html_escape(&metric.label),
                html_escape(&metric.detail)
            )
        })
        .collect::<String>();
    let figure_count = report
        .genes
        .iter()
        .map(|gene| gene.figures.len())
        .sum::<usize>();
    let primary_count = report
        .primers
        .iter()
        .filter(|primer| {
            primer.role.eq_ignore_ascii_case("primary")
                || primer.status.eq_ignore_ascii_case("primary")
        })
        .count();
    let overview = format!(
        "<section class=\"overview\" aria-label=\"Analysis overview\"><div class=\"overview-card\"><small>Evidence</small><b>{} labelled conditions</b><span>Transcript models, abundance, differential activity, and specificity inputs</span></div><div class=\"overview-arrow\">→</div><div class=\"overview-card\"><small>Assay panel</small><b>{} pairs across {} genes</b><span>{} primary pairs; review gates and rejection states remain visible</span></div><div class=\"overview-arrow\">→</div><div class=\"overview-card\"><small>Interpretation</small><b>{} linked figures</b><span>Locus evidence, primer/exon maps, and expected virtual-gel products</span></div></section>",
        report.conditions.len(),
        report.primers.len(),
        report.genes.len(),
        primary_count,
        figure_count,
    );

    let mut sections = String::new();
    for section in &report.sections {
        sections.push_str(&format!(
            "<section id=\"{}\"><h2>{}</h2>",
            html_escape(&section.section_id),
            html_escape(&section.heading)
        ));
        for paragraph in &section.paragraphs {
            sections.push_str(&format!("<p>{}</p>", html_escape(paragraph)));
        }
        if !section.bullets.is_empty() {
            sections.push_str("<ul>");
            for bullet in &section.bullets {
                sections.push_str(&format!("<li>{}</li>", html_escape(bullet)));
            }
            sections.push_str("</ul>");
        }
        sections.push_str("</section>");
    }

    let all_primer_refs = report.primers.iter().collect::<Vec<_>>();
    let all_primers = render_html_primer_table(&all_primer_refs);

    let mut genes = String::new();
    for gene in &report.genes {
        let heading = if gene.heading.is_empty() {
            &gene.gene_symbol
        } else {
            &gene.heading
        };
        let gene_metrics = gene
            .metrics
            .iter()
            .map(|metric| {
                format!(
                    "<span><b>{}</b> {} <small>{}</small></span>",
                    html_escape(&metric.value),
                    html_escape(&metric.label),
                    html_escape(&metric.detail)
                )
            })
            .collect::<String>();
        let figures = gene
            .figures
            .iter()
            .map(|figure| {
                format!(
                    "<figure data-gentle-figure-kind=\"{}\"><a href=\"{}\"><img src=\"{}\" alt=\"{}\"></a><figcaption><b>{}</b><br>{}</figcaption></figure>",
                    html_escape(&figure.kind),
                    html_escape(&figure.web_path),
                    html_escape(&figure.web_path),
                    html_escape(&figure.alt_text),
                    html_escape(&figure.label),
                    html_escape(&figure.caption),
                )
            })
            .collect::<String>();
        let primer_rows = primer_rows_for_gene(report, gene);
        genes.push_str(&format!(
            "<article class=\"gene\" id=\"gene-{}\"><header><h2>{}</h2><p>{}</p><div class=\"gene-metrics\">{}</div></header><div class=\"figure-grid\">{}</div><h3>Primer pairs</h3>{}</article>",
            html_escape(&gene.gene_symbol),
            html_escape(heading),
            html_escape(&gene.summary),
            gene_metrics,
            figures,
            render_html_primer_table(&primer_rows),
        ));
    }

    let downloads = report
        .downloads
        .iter()
        .map(|download| {
            format!(
                "<a class=\"button\" href=\"{}\"><b>{}</b><small>{}</small></a>",
                html_escape(&download.web_path),
                html_escape(&download.label),
                html_escape(&download.description),
            )
        })
        .collect::<String>();
    let pdf_button = report
        .pdf_path
        .as_ref()
        .map(|path| {
            format!(
                "<a class=\"button\" href=\"{}\"><b>Printable companion</b><small>PDF generated from this same resolved record</small></a>",
                html_escape(path)
            )
        })
        .unwrap_or_default();
    let provenance = if report.provenance.is_empty() {
        String::new()
    } else {
        format!(
            "<section><h2>Provenance</h2><ul>{}</ul></section>",
            report
                .provenance
                .iter()
                .map(|row| format!("<li><code>{}</code></li>", html_escape(row)))
                .collect::<String>()
        )
    };

    format!(
        r##"<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>{favicon}
<style>
:root{{--ink:#172231;--muted:#5d6978;--line:#d9e0e8;--paper:#fff;--wash:#f4f7fa;--blue:#285fb8;--green:#16856b;--amber:#d38b17;--red:#b42318}}
*{{box-sizing:border-box}} body{{margin:0;background:#eaf0f5;color:var(--ink);font:15px/1.48 system-ui,-apple-system,Segoe UI,sans-serif}}
main{{max-width:1240px;margin:28px auto;background:var(--paper);padding:34px 38px 46px;border-radius:16px;box-shadow:0 12px 35px #18334c20}}
h1{{font-size:32px;margin:0 0 6px}} h2{{margin-top:34px;border-bottom:2px solid var(--line);padding-bottom:7px}} h3{{margin-top:22px}}
.eyebrow{{color:var(--blue);font-weight:800;letter-spacing:.06em;text-transform:uppercase}} .subtitle,.lede{{color:var(--muted);max-width:980px}} .lede{{font-size:17px}}
.metrics{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:24px 0}} .metric{{background:var(--wash);border:1px solid var(--line);border-radius:11px;padding:14px}}
.metric>b{{font-size:24px;display:block;color:var(--blue)}} .metric strong,.metric small,.button small{{display:block}} small,.muted{{color:var(--muted)}}
.overview{{display:grid;grid-template-columns:1fr auto 1fr auto 1fr;gap:12px;align-items:center;margin:26px 0;padding:18px;background:linear-gradient(135deg,#eef6ff,#f3faf7);border:1px solid var(--line);border-radius:13px}} .overview-card{{min-height:124px;padding:16px;background:#ffffffd9;border:1px solid var(--line);border-radius:10px}} .overview-card small,.overview-card b,.overview-card span{{display:block}} .overview-card small{{color:var(--blue);font-weight:800;text-transform:uppercase;letter-spacing:.05em}} .overview-card b{{font-size:19px;margin:5px 0}} .overview-card span{{color:var(--muted)}} .overview-arrow{{color:var(--blue);font-size:28px;font-weight:800}}
.conditions code{{margin-right:5px}} .actions{{display:flex;flex-wrap:wrap;gap:10px;margin:20px 0}} .button{{text-decoration:none;color:var(--ink);border:1px solid var(--line);border-radius:9px;padding:10px 13px;background:var(--wash)}} .button.primary{{background:var(--blue);color:white;border-color:var(--blue)}} .button.primary small{{color:#eaf2ff}}
.gene{{margin:34px 0 48px;padding-top:2px;break-before:page}} .gene-metrics{{display:flex;gap:10px;flex-wrap:wrap}} .gene-metrics span{{background:var(--wash);border-radius:999px;padding:5px 10px}}
.figure-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(290px,1fr));gap:14px;align-items:start;margin:18px 0}} figure{{margin:0;border:1px solid var(--line);border-radius:10px;background:#fbfcfd;overflow:hidden}} figure img{{display:block;width:100%;height:auto;max-height:620px;object-fit:contain;background:white}} figcaption{{padding:10px 12px;color:var(--muted)}}
.table-scroll{{overflow:auto;border:1px solid var(--line);border-radius:9px}} table{{width:100%;border-collapse:collapse;font-size:12.5px}} th,td{{text-align:left;vertical-align:top;border-bottom:1px solid var(--line);padding:8px 7px}} th{{background:var(--wash);position:sticky;top:0}} code{{font-size:11.5px;overflow-wrap:anywhere}}
.pill{{display:inline-block;border-radius:999px;padding:2px 8px;color:white;font-size:10px}} .pill.ready{{background:var(--green)}} .pill.warning{{background:var(--amber)}} .pill.danger{{background:var(--red)}} footer{{margin-top:38px;color:var(--muted);border-top:1px solid var(--line);padding-top:16px}}
@media(max-width:760px){{main{{margin:0;padding:22px;border-radius:0}}.figure-grid,.overview{{grid-template-columns:1fr}}.overview-arrow{{transform:rotate(90deg);text-align:center}}}}
@media print{{body{{background:white;font-size:10pt}}main{{max-width:none;margin:0;padding:0;box-shadow:none}}.actions{{display:none}}.figure-grid{{grid-template-columns:1fr}}figure{{break-inside:avoid}}th{{position:static}}}}
</style></head><body><main data-gentle-schema="{schema}">
<header><div class="eyebrow">{eyebrow}</div><h1>{title}</h1><p class="subtitle">{subtitle}</p><p class="lede">{summary}</p>{conditions}<div class="metrics">{metrics}</div>{overview}<div class="actions">{pdf_button}{downloads}</div></header>
{sections}
<section id="primer-pairs"><h2>Primer pairs</h2>{all_primers}</section>
<section id="gene-reports"><h2>Gene reports and figures</h2>{genes}</section>
{provenance}<footer>{footer}</footer>
</main></body></html>"##,
        schema = html_escape(&report.schema),
        eyebrow = html_escape(&report.eyebrow),
        title = html_escape(&report.title),
        subtitle = html_escape(&report.subtitle),
        summary = html_escape(&report.summary),
        conditions = conditions,
        metrics = metrics,
        overview = overview,
        pdf_button = pdf_button,
        downloads = downloads,
        sections = sections,
        all_primers = all_primers,
        genes = genes,
        provenance = provenance,
        footer = html_escape(&report.footer),
        favicon = favicon,
    )
}

fn render_markdown_primer_list(rows: &[&GeneSetPublicationPrimerRow]) -> String {
    if rows.is_empty() {
        return "_No primer rows were supplied for this gene._\n".to_string();
    }
    let mut markdown = String::new();
    let mut previous_gene = "";
    for row in rows {
        if row.gene != previous_gene {
            markdown.push_str(&format!("## {}\n\n", markdown_cell(&row.gene)));
            previous_gene = &row.gene;
        }
        markdown.push_str(&format!(
            "**{}** — {} (`{}`)\n\n- Forward 5′→3′: `{}`; Tm {} °C\n- Reverse 5′→3′: `{}`; Tm {} °C\n- Specificity: cDNA `{}`; genome `{}`{}\n\n",
            markdown_cell(&row.pair_id),
            markdown_cell(&row.role),
            markdown_cell(&row.status),
            markdown_cell(&row.forward_sequence_5_to_3),
            markdown_cell(&row.forward_tm_c),
            markdown_cell(&row.reverse_sequence_5_to_3),
            markdown_cell(&row.reverse_tm_c),
            markdown_cell(&row.cdna_specificity),
            markdown_cell(&row.genome_assessment),
            if row.note.is_empty() { String::new() } else { format!(" — {}", markdown_cell(&row.note)) }
        ));
    }
    markdown
}

/// Render the print narrative and full primer list from the resolved report.
///
/// Figure assets are listed in their append order. The bundle generator joins
/// their PDF pages after this narrative; explicitly supplied PDF sources can
/// retain vector detail, while SVG-only sources use GENtle's rasterized PDF
/// fallback.
pub fn render_gene_set_publication_markdown(report: &GeneSetPublicationReport) -> String {
    let mut markdown = format!(
        "---\ntitle: \"{}\"\nsubtitle: \"{}\"\ndate: \"{}\"\nlang: en-GB\ngeometry: margin=20mm\nfontsize: 9pt\npapersize: a4\ncolorlinks: true\nheader-includes:\n  - |\n    ```{{=latex}}\n    \\usepackage{{booktabs}}\n    \\usepackage{{longtable}}\n    \\usepackage{{array}}\n    \\usepackage{{microtype}}\n    \\setlength{{\\emergencystretch}}{{3em}}\n    ```\n---\n\n# Purpose and summary\n\n{}\n\n",
        report.title.replace('"', "\\\""),
        report.subtitle.replace('"', "\\\""),
        report.generated_date.replace('"', "\\\""),
        report.summary
    );
    if !report.conditions.is_empty() {
        markdown.push_str(&format!(
            "**Conditions:** {}\n\n",
            report
                .conditions
                .iter()
                .map(|condition| format!("`{condition}`"))
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    if !report.metrics.is_empty() {
        markdown.push_str("| Value | Metric | Detail |\n|---:|---|---|\n");
        for metric in &report.metrics {
            markdown.push_str(&format!(
                "| {} | {} | {} |\n",
                markdown_cell(&metric.value),
                markdown_cell(&metric.label),
                markdown_cell(&metric.detail)
            ));
        }
        markdown.push('\n');
    }
    for section in &report.sections {
        markdown.push_str(&format!("# {}\n\n", section.heading));
        for paragraph in &section.paragraphs {
            markdown.push_str(paragraph);
            markdown.push_str("\n\n");
        }
        for bullet in &section.bullets {
            markdown.push_str(&format!("- {bullet}\n"));
        }
        markdown.push('\n');
    }

    markdown.push_str("# Complete primer-pair list\n\n");
    let all_rows = report.primers.iter().collect::<Vec<_>>();
    markdown.push_str(&render_markdown_primer_list(&all_rows));
    markdown.push('\n');

    markdown.push_str("# Gene-by-gene guide to appended figures\n\n");
    for gene in &report.genes {
        let heading = if gene.heading.is_empty() {
            &gene.gene_symbol
        } else {
            &gene.heading
        };
        markdown.push_str(&format!("## {heading}\n\n{}\n\n", gene.summary));
        for figure in &gene.figures {
            markdown.push_str(&format!(
                "- **{} ({})** — {} ([web asset]({}))\n",
                figure.label, figure.kind, figure.caption, figure.web_path
            ));
        }
        markdown.push('\n');
    }
    if !report.provenance.is_empty() {
        markdown.push_str("# Provenance\n\n");
        for row in &report.provenance {
            markdown.push_str(&format!("- `{row}`\n"));
        }
        markdown.push('\n');
    }
    markdown.push_str(
        "# Figure appendix\n\nThe following pages contain the complete figures in the same gene and figure order used by the HTML report. Vector source files remain linked from the web version.\n",
    );
    markdown
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle_protocol::{
        GeneSetPublicationFigure, GeneSetPublicationGene, GeneSetPublicationMetric,
        GeneSetPublicationPrimerRow,
    };

    fn fixture_report() -> GeneSetPublicationReport {
        GeneSetPublicationReport {
            schema: "gentle.gene_set_publication_report.v1".to_string(),
            title: "Synthetic two-gene panel".to_string(),
            subtitle: "One source, two renderers".to_string(),
            summary: "A deterministic report fixture.".to_string(),
            metrics: vec![GeneSetPublicationMetric {
                value: "2".to_string(),
                label: "genes".to_string(),
                detail: "synthetic".to_string(),
            }],
            primers: vec![GeneSetPublicationPrimerRow {
                gene: "GENE1".to_string(),
                pair_id: "G1_P1".to_string(),
                role: "common abundance".to_string(),
                status: "primary".to_string(),
                forward_sequence_5_to_3: "ACGT".to_string(),
                reverse_sequence_5_to_3: "TGCA".to_string(),
                ..Default::default()
            }],
            genes: vec![GeneSetPublicationGene {
                gene_symbol: "GENE1".to_string(),
                summary: "First synthetic gene.".to_string(),
                figures: vec![GeneSetPublicationFigure {
                    figure_id: "g1_locus".to_string(),
                    kind: "locus_evidence".to_string(),
                    label: "Locus evidence".to_string(),
                    web_path: "figures/GENE1_locus.svg".to_string(),
                    caption: "Synthetic locus.".to_string(),
                    alt_text: "GENE1 locus".to_string(),
                    include_in_pdf: true,
                    ..Default::default()
                }],
                ..Default::default()
            }],
            ..Default::default()
        }
    }

    #[test]
    fn html_and_markdown_share_primer_and_figure_semantics() {
        let report = fixture_report();
        let html = render_gene_set_publication_html(&report);
        let markdown = render_gene_set_publication_markdown(&report);
        for marker in ["GENE1", "G1_P1", "ACGT", "Synthetic locus"] {
            assert!(html.contains(marker), "HTML missing {marker}");
            assert!(markdown.contains(marker), "Markdown missing {marker}");
        }
        assert!(html.contains("data-gentle-figure-kind=\"locus_evidence\""));
        assert!(html.contains("aria-label=\"Analysis overview\""));
        assert!(html.contains("1 primary pairs"));
        assert!(markdown.contains("Complete primer-pair list"));
    }
}
