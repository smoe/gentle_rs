//! Shared SVG rendering for transcript-aware promoter-reporter comparisons.
//!
//! The renderer consumes the portable engine report directly. It does not
//! rescore motifs, reinterpret transcript geometry, or infer promoter usage.

use crate::engine::{
    PromoterReporterArchitectureComparisonReport, PromoterReporterArchitectureKind,
    PromoterReporterAtgDisposition, PromoterReporterCutRunLaneState,
};
use gentle_render::{
    GeneLocusEvidenceOverlay, GeneLocusEvidenceOverlayRow, GeneLocusEvidenceOverlaySchematicTail,
    GeneLocusEvidenceOverlaySegment, render_gene_locus_evidence_with_overlay_svg,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

const SVG_WIDTH: f64 = 1680.0;
const LABEL_WIDTH: f64 = 500.0;
const RIGHT_MARGIN: f64 = 36.0;
const HEADER_HEIGHT: f64 = 132.0;
const AXIS_HEIGHT: f64 = 42.0;
const ROW_HEIGHT: f64 = 38.0;
const SECTION_GAP: f64 = 32.0;
const SECTION_SPACER: f64 = 10.0;
const FOOTER_HEIGHT: f64 = 92.0;

fn escape_svg_text(raw: &str) -> String {
    raw.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn architecture_label(kind: PromoterReporterArchitectureKind) -> &'static str {
    match kind {
        PromoterReporterArchitectureKind::TssProximalTranscriptional => "TSS-proximal",
        PromoterReporterArchitectureKind::Spliced5utrLucAtgReplacement => "spliced 5' UTR",
        PromoterReporterArchitectureKind::Genomic5utrLucAtgReplacement => "genomic 5' UTR",
        PromoterReporterArchitectureKind::EndogenousAtgRetainedFusion => "ATG-retained fusion",
    }
}

fn lane_state_label(state: PromoterReporterCutRunLaneState) -> &'static str {
    match state {
        PromoterReporterCutRunLaneState::Unevaluated => "unevaluated",
        PromoterReporterCutRunLaneState::NotPrepared => "not prepared",
        PromoterReporterCutRunLaneState::PreparedNoCompatibleEvidence => {
            "prepared; no compatible evidence"
        }
        PromoterReporterCutRunLaneState::Evaluated => "evaluated",
    }
}

fn genomic_coordinate_label(
    report: &PromoterReporterArchitectureComparisonReport,
    local_0based: usize,
) -> String {
    let Some(anchor) = report.source_provenance.genome_anchor.as_ref() else {
        return format!("{} bp", local_0based.saturating_add(1));
    };
    let genomic = if anchor.strand == Some('-') {
        anchor.end_1based.saturating_sub(local_0based)
    } else {
        anchor.start_1based.saturating_add(local_0based)
    };
    format!("{}:{}", anchor.chromosome, genomic)
}

fn visible_bounds(report: &PromoterReporterArchitectureComparisonReport) -> (usize, usize) {
    let mut ranges = Vec::new();
    for transcript in &report.transcripts {
        ranges.extend(transcript.transcript_exon_ranges_0based.iter().copied());
        ranges.extend(transcript.cds_ranges_0based.iter().copied());
        ranges.push((
            transcript.tss_local_0based,
            transcript.tss_local_0based.saturating_add(1),
        ));
    }
    for architecture in &report.architectures {
        ranges.extend(architecture.segments.iter().filter_map(|segment| {
            Some((
                segment.source_start_0based?,
                segment.source_end_0based_exclusive?,
            ))
        }));
    }
    if let Some(support) = report.cutrun_evidence.regulatory_support.as_ref() {
        ranges.extend(
            support
                .support_windows
                .iter()
                .map(|window| (window.local_start_0based, window.local_end_0based_exclusive)),
        );
    }
    if let Some(motifs) = report.theoretical_motif_hits.as_ref() {
        ranges.extend(motifs.rows.iter().map(|hit| {
            (
                hit.source_match_start_0based,
                hit.source_match_end_0based_exclusive,
            )
        }));
    }
    if let Some(cds) = report.common_cds_start_local_0based {
        ranges.push((cds, cds.saturating_add(1)));
    }
    ranges.retain(|(start, end)| end > start);
    let start = ranges.iter().map(|(start, _)| *start).min().unwrap_or(0);
    let end = ranges
        .iter()
        .map(|(_, end)| *end)
        .max()
        .unwrap_or_else(|| start.saturating_add(1));
    (start, end.max(start.saturating_add(1)))
}

fn position_to_x(position: usize, start: usize, end: usize) -> f64 {
    let plot_width = SVG_WIDTH - LABEL_WIDTH - RIGHT_MARGIN;
    let span = end.saturating_sub(start).max(1) as f64;
    let fraction = position.saturating_sub(start) as f64 / span;
    LABEL_WIDTH + fraction.clamp(0.0, 1.0) * plot_width
}

fn range_rect(
    svg: &mut String,
    start: usize,
    end: usize,
    view_start: usize,
    view_end: usize,
    y: f64,
    height: f64,
    fill: &str,
    class_name: &str,
) {
    if end <= start {
        return;
    }
    let x1 = position_to_x(start, view_start, view_end);
    let x2 = position_to_x(end, view_start, view_end);
    let _ = write!(
        svg,
        "<rect class=\"{}\" x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" rx=\"2\" fill=\"{}\"/>",
        class_name,
        x1,
        y,
        (x2 - x1).max(1.5),
        height,
        fill
    );
}

fn render_section_title(svg: &mut String, title: &str, y: f64) {
    let _ = write!(
        svg,
        "<text x=\"28\" y=\"{:.2}\" class=\"section-title\">{}</text><line x1=\"28\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" class=\"section-line\"/>",
        y,
        escape_svg_text(title),
        y + 7.0,
        SVG_WIDTH - 28.0,
        y + 7.0
    );
}

fn normalized_locus_overlay(
    report: &PromoterReporterArchitectureComparisonReport,
) -> GeneLocusEvidenceOverlay {
    let rows = report
        .architectures
        .iter()
        .map(|architecture| {
            let terminal_anchor_local_1based =
                architecture.segments.iter().rev().find_map(|segment| {
                    let start = segment.source_start_0based?;
                    let end = segment.source_end_0based_exclusive?;
                    (end > start).then_some(if architecture.strand == "-" {
                        start.saturating_add(1)
                    } else {
                        end
                    })
                });
            let segments = architecture
                .segments
                .iter()
                .filter_map(|segment| {
                    let start = segment.source_start_0based?.saturating_add(1);
                    let end = segment.source_end_0based_exclusive?;
                    (end >= start).then(|| GeneLocusEvidenceOverlaySegment {
                        segment_id: segment.segment_id.clone(),
                        material: segment.material.clone(),
                        local_start_1based: start,
                        local_end_1based: end,
                        fill: match segment.material.as_str() {
                            "spliced_cdna" | "spliced_genomic_exon" => "#d58a3a",
                            "reporter_sequence" => "#b54b43",
                            _ => "#27847c",
                        }
                        .to_string(),
                    })
                })
                .collect::<Vec<_>>();
            let atg_position = architecture
                .cds_start_codon_source_ranges_0based
                .first()
                .map(|(start, end)| {
                    if architecture.strand == "-" {
                        end.saturating_sub(1)
                    } else {
                        *start
                    }
                })
                .map(|position| position.saturating_add(1));
            GeneLocusEvidenceOverlayRow {
                row_id: architecture.architecture_id.clone(),
                label: format!(
                    "{} · {}",
                    architecture.transcript_id,
                    architecture_label(architecture.architecture_kind)
                ),
                detail: format!(
                    "{} bp | {} | TSS class {}",
                    architecture.modeled_insert_length_bp,
                    architecture.junction.endogenous_atg_disposition.as_str(),
                    architecture.tss_class_id
                ),
                segments,
                marker_local_1based: atg_position,
                marker_label: Some("ATG boundary".to_string()),
                schematic_tail: terminal_anchor_local_1based.map(|anchor_local_1based| {
                    GeneLocusEvidenceOverlaySchematicTail {
                        segment_id: format!("{}_luciferase", architecture.architecture_id),
                        label: "LUC".to_string(),
                        detail: format!(
                            "{}; luciferase coding body is schematic and not to genomic scale",
                            architecture.junction.reporter_atg_source
                        ),
                        fill: "#b54b43".to_string(),
                        anchor_local_1based,
                    }
                }),
            }
        })
        .collect();
    GeneLocusEvidenceOverlay {
        overlay_id: "promoter_reporter_architecture_comparison".to_string(),
        title: "Proposed reporter architectures".to_string(),
        document_title: format!(
            "{} promoter-reporter architectures",
            report
                .gene_label
                .as_deref()
                .filter(|label| !label.trim().is_empty())
                .unwrap_or(&report.seq_id)
        ),
        summary: format!(
            "Canonical reporter comparison: {} TSS class(es), {} transcript(s), {} explicit architecture row(s); normalized evidence remains observational/predictive context.",
            report.tss_classes.len(),
            report.transcripts.len(),
            report.architectures.len()
        ),
        rows,
        non_claims: report.non_claims.clone(),
    }
}

/// Render one portable comparison report on a shared local/genomic axis.
pub fn render_promoter_reporter_architecture_svg(
    report: &PromoterReporterArchitectureComparisonReport,
) -> String {
    if let Some(locus_evidence) = report.locus_evidence.as_ref() {
        let overlay = normalized_locus_overlay(report);
        return render_gene_locus_evidence_with_overlay_svg(locus_evidence, Some(&overlay));
    }
    let motif_lane_count = report
        .theoretical_motif_hits
        .as_ref()
        .map(|motifs| {
            motifs
                .rows
                .iter()
                .map(|hit| hit.tf_id.clone())
                .collect::<BTreeSet<_>>()
                .len()
        })
        .unwrap_or_default();
    let cutrun_rows = usize::from(report.cutrun_evidence.regulatory_support.is_some())
        + report.cutrun_evidence.lanes.len().max(1);
    let total_rows = report.transcripts.len()
        + report.architectures.len()
        + cutrun_rows
        + motif_lane_count.max(1);
    let section_count = 4.0;
    let height = HEADER_HEIGHT
        + AXIS_HEIGHT
        + total_rows as f64 * ROW_HEIGHT
        + section_count * SECTION_GAP
        + (section_count - 1.0) * SECTION_SPACER
        + FOOTER_HEIGHT;
    let (view_start, view_end) = visible_bounds(report);
    let title_gene = report
        .gene_label
        .as_deref()
        .filter(|label| !label.trim().is_empty())
        .unwrap_or(&report.seq_id);
    let tss_labels = report
        .tss_classes
        .iter()
        .enumerate()
        .map(|(index, class)| (class.class_id.clone(), format!("TSS {}", index + 1)))
        .collect::<BTreeMap<_, _>>();
    let mut svg = String::new();
    let _ = write!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{:.0}\" viewBox=\"0 0 {} {:.0}\" role=\"img\" aria-labelledby=\"title desc\"><title id=\"title\">Promoter-reporter architecture comparison for {}</title><desc id=\"desc\">Transcript annotations, proposed reporter architectures, CUT&amp;RUN occupancy support, and theoretical motif hits aligned on one coordinate axis. Occupancy does not establish promoter use.</desc><style>text{{font-family:'Avenir Next','Gill Sans',sans-serif;fill:#243238}}.title{{font-size:25px;font-weight:700}}.subtitle{{font-size:12px;fill:#526268}}.section-title{{font-size:15px;font-weight:700;fill:#173f43}}.section-line{{stroke:#b8c8c8;stroke-width:1}}.row-label{{font-size:11px;fill:#33464b}}.row-meta{{font-size:9px;fill:#69797d}}.axis{{stroke:#60767a;stroke-width:1}}.axis-label{{font-size:10px;fill:#596d72}}.guide{{stroke:#d9e2df;stroke-width:1;stroke-dasharray:3 4}}.tss{{stroke:#db6b3f;stroke-width:2}}.cds-start{{stroke:#d29d21;stroke-width:2}}.transcript-line{{stroke:#7a8d91;stroke-width:1.2}}.synthetic-link{{stroke:#ba6b32;stroke-width:1.3;stroke-dasharray:4 3;fill:none}}.nonclaim{{font-size:10px;fill:#6c4d43}}</style><rect width=\"100%\" height=\"100%\" fill=\"#fbfaf4\"/><rect x=\"0\" y=\"0\" width=\"100%\" height=\"102\" fill=\"#eaf1ed\"/>",
        SVG_WIDTH as usize,
        height,
        SVG_WIDTH as usize,
        height,
        escape_svg_text(title_gene)
    );
    let _ = write!(
        svg,
        "<text x=\"28\" y=\"38\" class=\"title\">{} promoter-reporter architectures</text><text x=\"28\" y=\"64\" class=\"subtitle\">sequence {} | {} TSS class(es) | {} transcript(s) | {} explicit architecture row(s)</text><text x=\"28\" y=\"86\" class=\"subtitle\">orange TSS | gold common CDS start | teal genomic DNA | amber spliced cDNA | blue CDS | violet occupancy support</text>",
        escape_svg_text(title_gene),
        escape_svg_text(&report.seq_id),
        report.tss_classes.len(),
        report.transcripts.len(),
        report.architectures.len()
    );

    let axis_y = HEADER_HEIGHT;
    let _ = write!(
        svg,
        "<line x1=\"{}\" y1=\"{:.2}\" x2=\"{}\" y2=\"{:.2}\" class=\"axis\"/>",
        LABEL_WIDTH,
        axis_y,
        SVG_WIDTH - RIGHT_MARGIN,
        axis_y
    );
    for tick in 0..=4 {
        let position = view_start
            + (view_end.saturating_sub(view_start) as f64 * tick as f64 / 4.0).round() as usize;
        let x = position_to_x(position, view_start, view_end);
        let _ = write!(
            svg,
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" class=\"guide\"/><text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"middle\" class=\"axis-label\">{}</text>",
            x,
            axis_y,
            x,
            height - FOOTER_HEIGHT,
            x,
            axis_y - 8.0,
            escape_svg_text(&genomic_coordinate_label(report, position))
        );
    }
    for (class_index, class) in report.tss_classes.iter().enumerate() {
        let x = position_to_x(class.representative_tss_local_0based, view_start, view_end);
        let _ = write!(
            svg,
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" class=\"tss\" opacity=\"0.35\"/><text x=\"{:.2}\" y=\"{:.2}\" class=\"axis-label\">TSS {}</text>",
            x,
            axis_y,
            x,
            height - FOOTER_HEIGHT,
            x + 3.0,
            axis_y + 20.0,
            class_index + 1
        );
    }
    if let Some(cds) = report.common_cds_start_local_0based {
        let x = position_to_x(cds, view_start, view_end);
        let _ = write!(
            svg,
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" class=\"cds-start\" opacity=\"0.45\"/><text x=\"{:.2}\" y=\"{:.2}\" class=\"axis-label\">common CDS start</text>",
            x,
            axis_y,
            x,
            height - FOOTER_HEIGHT,
            x + 5.0,
            axis_y + 34.0
        );
    }

    let mut y = HEADER_HEIGHT + AXIS_HEIGHT;
    render_section_title(&mut svg, "Annotated transcript models", y);
    y += SECTION_GAP;
    for transcript in &report.transcripts {
        let center = y + ROW_HEIGHT * 0.5;
        let tss_label = tss_labels
            .get(&transcript.tss_class_id)
            .map(String::as_str)
            .unwrap_or("TSS ?");
        let _ = write!(
            svg,
            "<g data-transcript-id=\"{}\"><text x=\"28\" y=\"{:.2}\" class=\"row-label\">{}</text><text x=\"28\" y=\"{:.2}\" class=\"row-meta\">{} | spliced UTR {} bp | introns {} bp</text><line x1=\"{}\" y1=\"{:.2}\" x2=\"{}\" y2=\"{:.2}\" class=\"transcript-line\"/>",
            escape_svg_text(&transcript.transcript_id),
            y + 13.0,
            escape_svg_text(&transcript.transcript_id),
            y + 28.0,
            tss_label,
            transcript.five_prime_utr_spliced_length_bp,
            transcript.five_prime_utr_total_intron_bp,
            LABEL_WIDTH,
            center,
            SVG_WIDTH - RIGHT_MARGIN,
            center
        );
        for (start, end) in &transcript.transcript_exon_ranges_0based {
            range_rect(
                &mut svg,
                *start,
                *end,
                view_start,
                view_end,
                center - 4.0,
                8.0,
                "#78a9a2",
                "transcript-exon",
            );
        }
        for (start, end) in &transcript.cds_ranges_0based {
            range_rect(
                &mut svg,
                *start,
                *end,
                view_start,
                view_end,
                center - 6.0,
                12.0,
                "#225f73",
                "transcript-cds",
            );
        }
        let tss_x = position_to_x(transcript.tss_local_0based, view_start, view_end);
        let _ = write!(
            svg,
            "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" class=\"tss\"/></g>",
            tss_x,
            center - 8.0,
            tss_x,
            center + 8.0
        );
        y += ROW_HEIGHT;
    }

    y += SECTION_SPACER;
    render_section_title(&mut svg, "Proposed reporter architectures", y);
    y += SECTION_GAP;
    for architecture in &report.architectures {
        let center = y + ROW_HEIGHT * 0.5;
        let atg_label = match architecture.junction.endogenous_atg_disposition {
            PromoterReporterAtgDisposition::EndogenousExcludedVectorLuciferaseRetained => {
                "endogenous ATG excluded"
            }
            PromoterReporterAtgDisposition::EndogenousReplacedByLuciferase => {
                "endogenous ATG replaced"
            }
            PromoterReporterAtgDisposition::EndogenousRetainedVectorLuciferaseRemoved => {
                "endogenous ATG retained"
            }
            PromoterReporterAtgDisposition::EndogenousRetainedExtraUpstreamOfVectorLuciferase => {
                "two ATGs retained"
            }
        };
        let _ = write!(
            svg,
            "<g data-architecture-id=\"{}\"><text x=\"28\" y=\"{:.2}\" class=\"row-label\">{} · {}</text><text x=\"28\" y=\"{:.2}\" class=\"row-meta\">{} bp | {}</text>",
            escape_svg_text(&architecture.architecture_id),
            y + 13.0,
            escape_svg_text(&architecture.transcript_id),
            architecture_label(architecture.architecture_kind),
            y + 28.0,
            architecture.modeled_insert_length_bp,
            atg_label
        );
        let mut source_segments = architecture
            .segments
            .iter()
            .filter_map(|segment| {
                Some((
                    segment,
                    segment.source_start_0based?,
                    segment.source_end_0based_exclusive?,
                ))
            })
            .collect::<Vec<_>>();
        source_segments.sort_by_key(|(_, start, end)| (*start, *end));
        for (segment, start, end) in &source_segments {
            let fill = match segment.material.as_str() {
                "spliced_cdna" | "spliced_genomic_exon" => "#d58a3a",
                "reporter_sequence" => "#b54b43",
                _ => "#27847c",
            };
            range_rect(
                &mut svg,
                *start,
                *end,
                view_start,
                view_end,
                center - 5.0,
                10.0,
                fill,
                "architecture-segment",
            );
        }
        for pair in source_segments.windows(2) {
            let x1 = position_to_x(pair[0].2, view_start, view_end);
            let x2 = position_to_x(pair[1].1, view_start, view_end);
            if x2 > x1 + 1.0 {
                let _ = write!(
                    svg,
                    "<path d=\"M {:.2} {:.2} Q {:.2} {:.2} {:.2} {:.2}\" class=\"synthetic-link\"/>",
                    x1,
                    center,
                    (x1 + x2) * 0.5,
                    center - 8.0,
                    x2,
                    center
                );
            }
        }
        let atg_position = architecture
            .cds_start_codon_source_ranges_0based
            .first()
            .map(|(start, end)| {
                if architecture.strand == "-" {
                    end.saturating_sub(1)
                } else {
                    *start
                }
            })
            .unwrap_or(architecture.tss_local_0based);
        let atg_x = position_to_x(atg_position, view_start, view_end);
        let _ = write!(
            svg,
            "<path d=\"M {:.2} {:.2} l 5 5 l -5 5 l -5 -5 z\" fill=\"#d29d21\"/></g>",
            atg_x,
            center - 5.0
        );
        y += ROW_HEIGHT;
    }

    y += SECTION_SPACER;
    render_section_title(&mut svg, "CUT&RUN occupancy evidence", y);
    y += SECTION_GAP;
    if let Some(support) = report.cutrun_evidence.regulatory_support.as_ref() {
        let center = y + ROW_HEIGHT * 0.5;
        let _ = write!(
            svg,
            "<text x=\"28\" y=\"{:.2}\" class=\"row-label\">aggregate support windows</text><text x=\"28\" y=\"{:.2}\" class=\"row-meta\">qualitative; {} window(s)</text>",
            y + 13.0,
            y + 28.0,
            support.support_windows.len()
        );
        for window in &support.support_windows {
            range_rect(
                &mut svg,
                window.local_start_0based,
                window.local_end_0based_exclusive,
                view_start,
                view_end,
                center - 6.0,
                12.0,
                "#8063a6",
                "cutrun-support-window",
            );
        }
        y += ROW_HEIGHT;
    }
    if report.cutrun_evidence.lanes.is_empty() {
        let _ = write!(
            svg,
            "<text x=\"28\" y=\"{:.2}\" class=\"row-label\">No CUT&amp;RUN lane selected</text><text x=\"{}\" y=\"{:.2}\" class=\"row-meta\">occupancy unevaluated; no source acquired automatically</text>",
            y + 15.5,
            LABEL_WIDTH,
            y + 15.5
        );
        y += ROW_HEIGHT;
    } else {
        for lane in &report.cutrun_evidence.lanes {
            let details = [
                lane.tissue_or_cell_type.as_deref(),
                lane.condition.as_deref(),
                lane.sample_label.as_deref(),
            ]
            .into_iter()
            .flatten()
            .filter(|value| !value.trim().is_empty())
            .collect::<Vec<_>>()
            .join(" · ");
            let _ = write!(
                svg,
                "<g data-cutrun-source-id=\"{}\"><text x=\"28\" y=\"{:.2}\" class=\"row-label\">{}</text><text x=\"28\" y=\"{:.2}\" class=\"row-meta\">{} | {}{}</text></g>",
                escape_svg_text(&lane.source_id),
                y + 13.0,
                escape_svg_text(&lane.source_id),
                y + 28.0,
                escape_svg_text(&lane.role),
                lane_state_label(lane.state),
                if details.is_empty() {
                    String::new()
                } else {
                    format!(" | {}", escape_svg_text(&details))
                }
            );
            y += ROW_HEIGHT;
        }
    }

    y += SECTION_SPACER;
    render_section_title(&mut svg, "Theoretical motif layer", y);
    y += SECTION_GAP;
    if let Some(motifs) = report.theoretical_motif_hits.as_ref() {
        let mut grouped = BTreeMap::<String, Vec<_>>::new();
        for hit in &motifs.rows {
            grouped.entry(hit.tf_id.clone()).or_default().push(hit);
        }
        if grouped.is_empty() {
            let _ = write!(
                svg,
                "<text x=\"28\" y=\"{:.2}\" class=\"row-label\">No motif hit passed the requested threshold</text>",
                y + 15.5
            );
        } else {
            for (tf_id, hits) in grouped {
                let center = y + ROW_HEIGHT * 0.5;
                let _ = write!(
                    svg,
                    "<text x=\"28\" y=\"{:.2}\" class=\"row-label\">{}</text><text x=\"28\" y=\"{:.2}\" class=\"row-meta\">{} theoretical hit(s)</text>",
                    y + 13.0,
                    escape_svg_text(&tf_id),
                    y + 28.0,
                    hits.len()
                );
                for hit in hits {
                    range_rect(
                        &mut svg,
                        hit.source_match_start_0based,
                        hit.source_match_end_0based_exclusive,
                        view_start,
                        view_end,
                        center - 4.0,
                        8.0,
                        if hit.forward_strand {
                            "#2c8194"
                        } else {
                            "#c86c3c"
                        },
                        "theoretical-motif-hit",
                    );
                }
                y += ROW_HEIGHT;
            }
        }
    } else {
        let _ = write!(
            svg,
            "<text x=\"28\" y=\"{:.2}\" class=\"row-label\">No theoretical motif scan requested</text>",
            y + 15.5
        );
    }

    let footer_y = height - FOOTER_HEIGHT + 14.0;
    let _ = write!(
        svg,
        "<rect x=\"0\" y=\"{:.2}\" width=\"100%\" height=\"{}\" fill=\"#f3ece4\"/><text x=\"28\" y=\"{:.2}\" class=\"nonclaim\">Interpretation boundary: {}</text><text x=\"28\" y=\"{:.2}\" class=\"nonclaim\">{}</text><text x=\"28\" y=\"{:.2}\" class=\"nonclaim\">CUT&amp;RUN status: {} | quantitative comparison: {}</text></svg>",
        height - FOOTER_HEIGHT,
        FOOTER_HEIGHT as usize,
        footer_y,
        escape_svg_text(
            report
                .non_claims
                .first()
                .map(String::as_str)
                .unwrap_or("Annotation does not establish promoter usage.")
        ),
        footer_y + 19.0,
        escape_svg_text(
            report
                .non_claims
                .get(1)
                .map(String::as_str)
                .unwrap_or("Occupancy does not establish reporter activity.")
        ),
        footer_y + 38.0,
        escape_svg_text(&report.cutrun_evidence.evaluation_state),
        escape_svg_text(&report.cutrun_evidence.quantitative_comparison_status)
    );
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{
        PromoterReporterArchitectureJunctionAudit, PromoterReporterArchitectureRow,
        PromoterReporterArchitectureSegment, PromoterReporterCutRunEvidence,
        PromoterReporterTranscriptArchitectureAudit,
    };

    #[test]
    fn normalized_overlay_anchors_schematic_luciferase_after_each_transcript_oriented_insert() {
        let source_segment =
            |id: &str, start: usize, end: usize| PromoterReporterArchitectureSegment {
                segment_id: id.to_string(),
                material: "genomic_dna".to_string(),
                source_start_0based: Some(start),
                source_end_0based_exclusive: Some(end),
                length_bp: end - start,
                ..Default::default()
            };
        let report = PromoterReporterArchitectureComparisonReport {
            architectures: vec![
                PromoterReporterArchitectureRow {
                    architecture_id: "plus_arch".to_string(),
                    transcript_id: "PLUS-201".to_string(),
                    strand: "+".to_string(),
                    segments: vec![source_segment("plus_insert", 100, 120)],
                    junction: PromoterReporterArchitectureJunctionAudit {
                        reporter_atg_source: "vector_luciferase_atg".to_string(),
                        ..Default::default()
                    },
                    ..Default::default()
                },
                PromoterReporterArchitectureRow {
                    architecture_id: "minus_arch".to_string(),
                    transcript_id: "MINUS-201".to_string(),
                    strand: "-".to_string(),
                    segments: vec![
                        source_segment("minus_promoter", 300, 320),
                        source_segment("minus_leader", 250, 270),
                    ],
                    junction: PromoterReporterArchitectureJunctionAudit {
                        reporter_atg_source: "luciferase_cds_atg".to_string(),
                        ..Default::default()
                    },
                    ..Default::default()
                },
            ],
            ..Default::default()
        };

        let overlay = normalized_locus_overlay(&report);
        assert_eq!(
            overlay.rows[0]
                .schematic_tail
                .as_ref()
                .expect("plus LUC tail")
                .anchor_local_1based,
            120
        );
        assert_eq!(
            overlay.rows[1]
                .schematic_tail
                .as_ref()
                .expect("minus LUC tail")
                .anchor_local_1based,
            251
        );
    }

    #[test]
    fn svg_keeps_architecture_and_non_claims_visible() {
        let report = PromoterReporterArchitectureComparisonReport {
            seq_id: "synthetic_promoter".to_string(),
            gene_label: Some("SYN1".to_string()),
            transcripts: vec![PromoterReporterTranscriptArchitectureAudit {
                transcript_id: "SYN1-201".to_string(),
                tss_class_id: "main_tss".to_string(),
                tss_local_0based: 100,
                transcript_exon_ranges_0based: vec![(100, 180), (250, 400)],
                cds_ranges_0based: vec![(300, 400)],
                ..Default::default()
            }],
            architectures: vec![PromoterReporterArchitectureRow {
                architecture_id: "arch_main_spliced".to_string(),
                architecture_kind: PromoterReporterArchitectureKind::Spliced5utrLucAtgReplacement,
                transcript_id: "SYN1-201".to_string(),
                cds_start_codon_source_ranges_0based: vec![(300, 303)],
                junction: PromoterReporterArchitectureJunctionAudit {
                    endogenous_atg_disposition:
                        PromoterReporterAtgDisposition::EndogenousReplacedByLuciferase,
                    ..Default::default()
                },
                ..Default::default()
            }],
            cutrun_evidence: PromoterReporterCutRunEvidence {
                evaluation_state: "unevaluated_no_sources_selected".to_string(),
                quantitative_comparison_status: "not_requested".to_string(),
                ..Default::default()
            },
            non_claims: vec![
                "Annotation does not establish promoter usage.".to_string(),
                "Occupancy does not establish reporter activity.".to_string(),
            ],
            ..Default::default()
        };

        let svg = render_promoter_reporter_architecture_svg(&report);
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("arch_main_spliced"));
        assert!(svg.contains("endogenous ATG replaced"));
        assert!(svg.contains("No CUT&amp;RUN lane selected"));
        assert!(svg.contains("Annotation does not establish promoter usage."));
    }
}
