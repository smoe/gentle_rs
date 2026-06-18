//! Deterministic SVG rendering for interpreted probe-region evidence reports.

use super::*;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::Path;

struct EvidenceSvgProjector {
    min_coord: usize,
    max_coord: usize,
    local_offset: Option<i128>,
}

impl EvidenceSvgProjector {
    fn project_report_coord(&self, coord: usize) -> usize {
        let Some(offset) = self.local_offset else {
            return coord;
        };
        (coord as i128 + offset).max(1) as usize
    }

    fn uses_local_alignment(&self) -> bool {
        self.local_offset.is_some()
    }
}

impl GentleEngine {
    pub fn export_probe_region_evidence_svg(
        &self,
        report_path: &str,
        svg_path: &str,
    ) -> Result<ProbeRegionEvidenceSvgExport, EngineError> {
        let report_path = report_path.trim();
        if report_path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Probe-region evidence SVG report path must not be empty".to_string(),
                cause_chain: vec![],
            });
        }
        let svg_path = svg_path.trim();
        if svg_path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Probe-region evidence SVG output path must not be empty".to_string(),
                cause_chain: vec![],
            });
        }
        let report_text = std::fs::read_to_string(report_path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not read probe-region evidence report '{}': {e}",
                report_path
            ),
            cause_chain: vec![],
        })?;
        let report: ProbeRegionEvidenceInterpretationReport = serde_json::from_str(&report_text)
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not parse probe-region evidence report '{}': {e}",
                    report_path
                ),
                cause_chain: vec![],
            })?;
        if report.schema != PROBE_REGION_EVIDENCE_INTERPRETATION_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Probe-region evidence report '{}' has schema '{}' but expected '{}'",
                    report_path, report.schema, PROBE_REGION_EVIDENCE_INTERPRETATION_SCHEMA
                ),
                cause_chain: vec![],
            });
        }

        let svg = Self::render_probe_region_evidence_svg_text(&report);
        if let Some(parent) = Path::new(svg_path).parent()
            && !parent.as_os_str().is_empty()
        {
            std::fs::create_dir_all(parent).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not create probe-region evidence SVG output directory '{}': {e}",
                    parent.to_string_lossy()
                ),
                cause_chain: vec![],
            })?;
        }
        std::fs::write(svg_path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write probe-region evidence SVG '{}': {e}",
                svg_path
            ),
            cause_chain: vec![],
        })?;

        Ok(Self::probe_region_evidence_svg_export_summary(
            &report,
            report_path,
            svg_path,
        ))
    }

    pub(crate) fn render_probe_region_evidence_svg_text(
        report: &ProbeRegionEvidenceInterpretationReport,
    ) -> String {
        let projector = Self::probe_region_evidence_svg_projector(report);
        let warnings = Self::probe_region_evidence_svg_warnings(report, &projector);
        let transcript_rows = Self::probe_region_evidence_svg_transcript_rows(report);
        let parent_groups = Self::probe_region_evidence_svg_parent_groups(report);
        let transcript_lane_height = 44.0;
        let parent_lane_height = 32.0;
        let transcript_top = 170.0;
        let probe_top =
            transcript_top + transcript_rows.len().max(1) as f64 * transcript_lane_height + 46.0;
        let legend_top = probe_top + parent_groups.len().max(1) as f64 * parent_lane_height + 72.0;
        let width = 1240.0;
        let height = (legend_top + 130.0).max(520.0);
        let left = 180.0;
        let right = 52.0;
        let plot_width = width - left - right;
        let mut svg = String::new();

        let _ = writeln!(
            svg,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width:.0}\" height=\"{height:.0}\" viewBox=\"0 0 {width:.0} {height:.0}\" data-gentle-schema=\"{}\" data-report-schema=\"{}\">",
            Self::probe_region_svg_escape(PROBE_REGION_EVIDENCE_SVG_EXPORT_SCHEMA),
            Self::probe_region_svg_escape(&report.schema)
        );
        svg.push_str("<rect x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n");
        svg.push_str("<style><![CDATA[\n");
        svg.push_str(
            ".axis{stroke:#475569;stroke-width:1}.grid{stroke:#e5e7eb;stroke-width:1}.transcript-baseline{stroke:#cbd5e1;stroke-width:2}.exon-segment{fill:#0f766e;stroke:#115e59;stroke-width:1}.junction-span{fill:none;stroke:#7c3aed;stroke-width:2;stroke-dasharray:4 3}.parent-probeset{fill:#fef3c7;stroke:#92400e;stroke-width:1}.pm-probe{stroke:#111827;stroke-width:.8}.region-summary{fill:#94a3b8;stroke:#475569;stroke-width:.8}.label{font-family:Inter,Arial,sans-serif;fill:#111827}.muted{font-family:Inter,Arial,sans-serif;fill:#64748b}.legend-box{fill:#f8fafc;stroke:#cbd5e1;stroke-width:1}\n",
        );
        svg.push_str("]]></style>\n");
        let title_gene = report
            .gene_label
            .as_deref()
            .or_else(|| {
                report
                    .transcript_rows
                    .iter()
                    .find_map(|row| row.gene.as_deref())
            })
            .unwrap_or("probe-region evidence");
        let _ = writeln!(
            svg,
            "<text class=\"label\" x=\"36\" y=\"42\" font-size=\"24\" font-weight=\"700\">Probe-region evidence geometry constraints</text>"
        );
        let _ = writeln!(
            svg,
            "<text class=\"label\" x=\"36\" y=\"70\" font-size=\"14\">{} level={} rows={} transcripts={}</text>",
            Self::probe_region_svg_escape(title_gene),
            Self::probe_region_svg_escape(&report.level),
            report.evidence_rows.len(),
            report.transcript_rows.len()
        );
        let _ = writeln!(
            svg,
            "<text class=\"muted\" x=\"36\" y=\"94\" font-size=\"12\">Review-only transcript/exon geometry constraints; no isoform support, probe specificity, or multi-hit verdict is inferred.</text>"
        );
        if projector.uses_local_alignment() {
            let _ = writeln!(
                svg,
                "<text class=\"muted\" x=\"36\" y=\"114\" font-size=\"11\">Overlapped exon/junction ranges are report-local and are aligned to the evidence coordinate axis for display; full gene models are not reconstructed.</text>"
            );
        }

        Self::render_probe_region_evidence_axis(
            &mut svg, &projector, left, 140.0, plot_width, title_gene,
        );
        Self::render_probe_region_evidence_transcripts(
            &mut svg,
            report,
            &projector,
            &transcript_rows,
            left,
            transcript_top,
            plot_width,
            transcript_lane_height,
        );
        Self::render_probe_region_evidence_parents(
            &mut svg,
            &projector,
            &parent_groups,
            left,
            probe_top,
            plot_width,
            parent_lane_height,
        );
        Self::render_probe_region_evidence_legend(
            &mut svg, report, &warnings, left, legend_top, plot_width,
        );
        svg.push_str("</svg>\n");
        svg
    }

    fn probe_region_evidence_svg_export_summary(
        report: &ProbeRegionEvidenceInterpretationReport,
        report_path: &str,
        svg_path: &str,
    ) -> ProbeRegionEvidenceSvgExport {
        let projector = Self::probe_region_evidence_svg_projector(report);
        ProbeRegionEvidenceSvgExport {
            schema: PROBE_REGION_EVIDENCE_SVG_EXPORT_SCHEMA.to_string(),
            report_path: report_path.to_string(),
            svg_path: svg_path.to_string(),
            evidence_row_count: report.evidence_rows.len(),
            transcript_count: report.transcript_rows.len(),
            parent_feature_count: Self::probe_region_evidence_svg_parent_groups(report).len(),
            junction_span_count: Self::probe_region_evidence_svg_junction_span_count(report),
            ambiguity_tags: Self::probe_region_evidence_svg_ambiguity_tags(report),
            warnings: Self::probe_region_evidence_svg_warnings(report, &projector),
        }
    }

    fn render_probe_region_evidence_axis(
        svg: &mut String,
        projector: &EvidenceSvgProjector,
        left: f64,
        y: f64,
        width: f64,
        label: &str,
    ) {
        let _ = writeln!(
            svg,
            "<g id=\"probe-region-evidence-axis\" class=\"coordinate-axis\" data-min=\"{}\" data-max=\"{}\">",
            projector.min_coord, projector.max_coord
        );
        let _ = writeln!(
            svg,
            "<line class=\"axis\" x1=\"{left:.1}\" y1=\"{y:.1}\" x2=\"{:.1}\" y2=\"{y:.1}\"/>",
            left + width
        );
        for frac in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let x = left + width * frac;
            let coord = projector.min_coord as f64
                + (projector.max_coord.saturating_sub(projector.min_coord) as f64 * frac);
            let _ = writeln!(
                svg,
                "<line class=\"grid\" x1=\"{x:.1}\" y1=\"{:.1}\" x2=\"{x:.1}\" y2=\"{:.1}\"/>",
                y + 6.0,
                y + 360.0
            );
            let _ = writeln!(
                svg,
                "<text class=\"muted\" x=\"{x:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"10\">{coord:.0}</text>",
                y - 8.0
            );
        }
        let _ = writeln!(
            svg,
            "<text class=\"muted\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"12\">{} coordinate span</text>",
            left + width / 2.0,
            y - 26.0,
            Self::probe_region_svg_escape(label)
        );
        svg.push_str("</g>\n");
    }

    #[allow(clippy::too_many_arguments)]
    fn render_probe_region_evidence_transcripts(
        svg: &mut String,
        report: &ProbeRegionEvidenceInterpretationReport,
        projector: &EvidenceSvgProjector,
        transcript_rows: &[ProbeRegionEvidenceTranscriptRow],
        left: f64,
        top: f64,
        width: f64,
        lane_height: f64,
    ) {
        let _ = writeln!(
            svg,
            "<g id=\"probe-region-evidence-transcripts\" class=\"transcript-lanes\">"
        );
        let _ = writeln!(
            svg,
            "<text class=\"label\" x=\"36\" y=\"{:.1}\" font-size=\"15\" font-weight=\"700\">Transcript overlap lanes</text>",
            top - 18.0
        );
        for (idx, tx) in transcript_rows.iter().enumerate() {
            let y = top + idx as f64 * lane_height;
            let lane_mid = y + 16.0;
            let safe_id = Self::probe_region_evidence_svg_id(&tx.transcript_id);
            let review = Self::probe_region_svg_escape(&tx.review_status);
            let label = tx.label.as_deref().unwrap_or(&tx.transcript_id);
            let _ = writeln!(
                svg,
                "<g id=\"transcript-{safe_id}\" class=\"transcript\" data-transcript-id=\"{}\" data-review-status=\"{review}\">",
                Self::probe_region_svg_escape(&tx.transcript_id)
            );
            let _ = writeln!(
                svg,
                "<text class=\"label\" x=\"36\" y=\"{:.1}\" font-size=\"11\">{} ({}) exons={} {}</text>",
                lane_mid + 4.0,
                Self::probe_region_svg_escape(label),
                Self::probe_region_svg_escape(tx.strand.as_deref().unwrap_or("?")),
                tx.exon_count,
                review
            );
            let _ = writeln!(
                svg,
                "<line class=\"transcript-baseline\" x1=\"{left:.1}\" y1=\"{lane_mid:.1}\" x2=\"{:.1}\" y2=\"{lane_mid:.1}\"/>",
                left + width
            );
            for exon in Self::probe_region_evidence_svg_exons_for_transcript(
                report,
                &tx.transcript_id,
                projector,
            ) {
                let x1 =
                    Self::probe_region_evidence_svg_x(exon.start_1based, projector, left, width);
                let x2 = Self::probe_region_evidence_svg_x(exon.end_1based, projector, left, width);
                let rect_x = x1.min(x2);
                let rect_w = (x1.max(x2) - rect_x).max(4.0);
                let _ = writeln!(
                    svg,
                    "<rect class=\"exon-segment\" x=\"{rect_x:.1}\" y=\"{:.1}\" width=\"{rect_w:.1}\" height=\"14\" rx=\"2\" data-exon-ordinal=\"{}\" data-original-range=\"{}\"><title>{} exon {} range {}</title></rect>",
                    lane_mid - 7.0,
                    exon.ordinal,
                    Self::probe_region_svg_escape(&exon.original_range),
                    Self::probe_region_svg_escape(&tx.transcript_id),
                    exon.ordinal,
                    Self::probe_region_svg_escape(&exon.original_range)
                );
            }
            for span in Self::probe_region_evidence_svg_junctions_for_transcript(
                report,
                &tx.transcript_id,
                projector,
            ) {
                let x1 =
                    Self::probe_region_evidence_svg_x(span.start_1based, projector, left, width);
                let x2 = Self::probe_region_evidence_svg_x(span.end_1based, projector, left, width);
                let c1 = x1.min(x2);
                let c2 = x1.max(x2);
                let peak = lane_mid - 18.0;
                let _ = writeln!(
                    svg,
                    "<path class=\"junction-span\" d=\"M {c1:.1} {lane_mid:.1} Q {:.1} {peak:.1} {c2:.1} {lane_mid:.1}\" data-from-exon=\"{}\" data-to-exon=\"{}\"><title>{} junction {} to {} {}</title></path>",
                    (c1 + c2) / 2.0,
                    span.from_exon,
                    span.to_exon,
                    Self::probe_region_svg_escape(&tx.transcript_id),
                    span.from_exon,
                    span.to_exon,
                    Self::probe_region_svg_escape(&span.original_range)
                );
            }
            svg.push_str("</g>\n");
        }
        svg.push_str("</g>\n");
    }

    #[allow(clippy::too_many_arguments)]
    fn render_probe_region_evidence_parents(
        svg: &mut String,
        projector: &EvidenceSvgProjector,
        parent_groups: &[EvidenceSvgParentGroup],
        left: f64,
        top: f64,
        width: f64,
        lane_height: f64,
    ) {
        let _ = writeln!(
            svg,
            "<g id=\"probe-region-evidence-probes\" class=\"probe-region-evidence-tracks\">"
        );
        let _ = writeln!(
            svg,
            "<text class=\"label\" x=\"36\" y=\"{:.1}\" font-size=\"15\" font-weight=\"700\">Probe and parent probeset intervals</text>",
            top - 20.0
        );
        for (idx, group) in parent_groups.iter().enumerate() {
            let y = top + idx as f64 * lane_height;
            let lane_mid = y + 13.0;
            let safe_id = Self::probe_region_evidence_svg_id(&group.parent_id);
            let parent_x1 =
                Self::probe_region_evidence_svg_x(group.start_1based, projector, left, width);
            let parent_x2 =
                Self::probe_region_evidence_svg_x(group.end_1based, projector, left, width);
            let parent_rect_x = parent_x1.min(parent_x2);
            let parent_rect_w = (parent_x1.max(parent_x2) - parent_rect_x).max(5.0);
            let _ = writeln!(
                svg,
                "<g id=\"parent-probeset-{safe_id}\" class=\"parent-probeset-group\" data-parent-feature-id=\"{}\">",
                Self::probe_region_svg_escape(&group.parent_id)
            );
            let _ = writeln!(
                svg,
                "<text class=\"label\" x=\"36\" y=\"{:.1}\" font-size=\"10\">{}</text>",
                lane_mid + 4.0,
                Self::probe_region_svg_escape(&group.parent_id)
            );
            let _ = writeln!(
                svg,
                "<rect class=\"parent-probeset\" x=\"{parent_rect_x:.1}\" y=\"{:.1}\" width=\"{parent_rect_w:.1}\" height=\"10\" rx=\"2\"><title>{} parent span {}..{}</title></rect>",
                lane_mid - 5.0,
                Self::probe_region_svg_escape(&group.parent_id),
                group.start_1based,
                group.end_1based
            );
            for evidence in &group.rows {
                let x1 = Self::probe_region_evidence_svg_x(
                    evidence.start_1based,
                    projector,
                    left,
                    width,
                );
                let x2 =
                    Self::probe_region_evidence_svg_x(evidence.end_1based, projector, left, width);
                let rect_x = x1.min(x2);
                let rect_w = (x1.max(x2) - rect_x).max(3.0);
                let class_name =
                    if evidence.intensity_source.as_deref() == Some("probe_level_input") {
                        "pm-probe"
                    } else {
                        "region-summary"
                    };
                let fill = Self::probe_region_evidence_svg_logfc_fill(evidence.logfc);
                let _ = writeln!(
                    svg,
                    "<rect class=\"{class_name}\" x=\"{rect_x:.1}\" y=\"{:.1}\" width=\"{rect_w:.1}\" height=\"18\" rx=\"2\" fill=\"{fill}\" data-feature-id=\"{}\" data-parent-feature-id=\"{}\"><title>{} {}..{} logFC={}</title></rect>",
                    lane_mid + 7.0,
                    Self::probe_region_svg_escape(&evidence.feature_id),
                    Self::probe_region_svg_escape(&group.parent_id),
                    Self::probe_region_svg_escape(&evidence.feature_id),
                    evidence.start_1based,
                    evidence.end_1based,
                    evidence
                        .logfc
                        .map(|value| format!("{value:.3}"))
                        .unwrap_or_else(|| "not available".to_string())
                );
            }
            svg.push_str("</g>\n");
        }
        svg.push_str("</g>\n");
    }

    fn render_probe_region_evidence_legend(
        svg: &mut String,
        report: &ProbeRegionEvidenceInterpretationReport,
        warnings: &[String],
        left: f64,
        top: f64,
        width: f64,
    ) {
        let tags = Self::probe_region_evidence_svg_ambiguity_tags(report);
        let review_statuses = report
            .transcript_rows
            .iter()
            .map(|row| row.review_status.clone())
            .collect::<BTreeSet<_>>();
        let _ = writeln!(
            svg,
            "<g id=\"probe-region-evidence-legend\" class=\"legend\" transform=\"translate({left:.1},{top:.1})\">"
        );
        let _ = writeln!(
            svg,
            "<rect class=\"legend-box\" x=\"0\" y=\"0\" width=\"{width:.1}\" height=\"92\" rx=\"4\"/>"
        );
        let _ = writeln!(
            svg,
            "<text class=\"label\" x=\"14\" y=\"24\" font-size=\"12\" font-weight=\"700\">Legend and guardrails</text>"
        );
        let _ = writeln!(
            svg,
            "<rect class=\"pm-probe\" x=\"16\" y=\"38\" width=\"20\" height=\"12\" fill=\"#dc2626\"/><text class=\"muted\" x=\"44\" y=\"49\" font-size=\"11\">true PM probe (probe_level_input); red/blue encode logFC sign</text>"
        );
        let _ = writeln!(
            svg,
            "<rect class=\"parent-probeset\" x=\"360\" y=\"38\" width=\"28\" height=\"10\"/><text class=\"muted\" x=\"396\" y=\"49\" font-size=\"11\">parent probeset span</text>"
        );
        let _ = writeln!(
            svg,
            "<path class=\"junction-span\" d=\"M 584 48 Q 606 30 628 48\"/><text class=\"muted\" x=\"638\" y=\"49\" font-size=\"11\">junction-spanning evidence from report</text>"
        );
        let _ = writeln!(
            svg,
            "<text class=\"muted\" x=\"14\" y=\"70\" font-size=\"11\">ambiguity labels: {}</text>",
            Self::probe_region_svg_escape(&tags.join(", "))
        );
        let _ = writeln!(
            svg,
            "<text class=\"muted\" x=\"14\" y=\"86\" font-size=\"11\">transcript review_status values: {}</text>",
            Self::probe_region_svg_escape(
                &review_statuses.into_iter().collect::<Vec<_>>().join(", ")
            )
        );
        for (idx, warning) in warnings.iter().enumerate() {
            let _ = writeln!(
                svg,
                "<text class=\"muted\" x=\"720\" y=\"{:.1}\" font-size=\"10\">{}</text>",
                24.0 + idx as f64 * 14.0,
                Self::probe_region_svg_escape(warning)
            );
        }
        svg.push_str("</g>\n");
    }

    fn probe_region_evidence_svg_projector(
        report: &ProbeRegionEvidenceInterpretationReport,
    ) -> EvidenceSvgProjector {
        let mut evidence_coords = Vec::new();
        for row in &report.evidence_rows {
            if let Some(start) = row.start_1based {
                evidence_coords.push(start);
            }
            if let Some(end) = row.end_1based {
                evidence_coords.push(end);
            }
        }
        let local_bounds = Self::probe_region_evidence_svg_local_geometry_bounds(report);
        let evidence_min = evidence_coords
            .iter()
            .copied()
            .min()
            .unwrap_or_else(|| local_bounds.map(|bounds| bounds.0).unwrap_or(1).max(1));
        let evidence_max = evidence_coords
            .iter()
            .copied()
            .max()
            .unwrap_or_else(|| local_bounds.map(|bounds| bounds.1).unwrap_or(evidence_min));
        let local_offset = local_bounds.and_then(|(local_min, local_max)| {
            (evidence_min > local_max.saturating_add(1000))
                .then_some(evidence_min as i128 - local_min as i128)
        });
        let mut coords = evidence_coords;
        if let Some((local_min, local_max)) = local_bounds {
            let projector = EvidenceSvgProjector {
                min_coord: evidence_min,
                max_coord: evidence_max,
                local_offset,
            };
            coords.push(projector.project_report_coord(local_min));
            coords.push(projector.project_report_coord(local_max));
        }
        for row in &report.evidence_rows {
            for mapping in &row.transcript_mappings {
                for span in &mapping.junction_spans {
                    let projector = EvidenceSvgProjector {
                        min_coord: evidence_min,
                        max_coord: evidence_max,
                        local_offset,
                    };
                    coords.push(projector.project_report_coord(span.genomic_start_1based));
                    coords.push(projector.project_report_coord(span.genomic_end_1based));
                }
            }
        }
        let mut min_coord = coords.iter().copied().min().unwrap_or(evidence_min).max(1);
        let mut max_coord = coords
            .iter()
            .copied()
            .max()
            .unwrap_or(evidence_max)
            .max(min_coord);
        if max_coord <= min_coord {
            min_coord = min_coord.saturating_sub(1).max(1);
            max_coord = max_coord.saturating_add(1);
        }
        EvidenceSvgProjector {
            min_coord,
            max_coord,
            local_offset,
        }
    }

    fn probe_region_evidence_svg_warnings(
        report: &ProbeRegionEvidenceInterpretationReport,
        projector: &EvidenceSvgProjector,
    ) -> Vec<String> {
        let mut warnings = report.warnings.clone();
        if projector.uses_local_alignment() {
            warnings.push(
                "report_local_geometry_aligned_to_evidence_axis_without_full_gene_model"
                    .to_string(),
            );
        }
        if report.evidence_rows.is_empty() {
            warnings.push("no_evidence_rows_to_render".to_string());
        }
        warnings
    }

    fn probe_region_evidence_svg_transcript_rows(
        report: &ProbeRegionEvidenceInterpretationReport,
    ) -> Vec<ProbeRegionEvidenceTranscriptRow> {
        let mut rows = report.transcript_rows.clone();
        rows.sort_by(|left, right| {
            left.transcript_id
                .cmp(&right.transcript_id)
                .then(left.label.cmp(&right.label))
        });
        rows
    }

    fn probe_region_evidence_svg_parent_groups(
        report: &ProbeRegionEvidenceInterpretationReport,
    ) -> Vec<EvidenceSvgParentGroup> {
        let mut grouped: BTreeMap<String, Vec<EvidenceSvgEvidenceRow>> = BTreeMap::new();
        for row in &report.evidence_rows {
            let (Some(start), Some(end)) = (row.start_1based, row.end_1based) else {
                continue;
            };
            let parent_id = row
                .parent_feature_id
                .clone()
                .unwrap_or_else(|| row.feature_id.clone());
            grouped
                .entry(parent_id)
                .or_default()
                .push(EvidenceSvgEvidenceRow {
                    feature_id: row.feature_id.clone(),
                    intensity_source: row.intensity_source.clone(),
                    start_1based: start.min(end),
                    end_1based: start.max(end),
                    logfc: row.logfc,
                });
        }
        grouped
            .into_iter()
            .map(|(parent_id, mut rows)| {
                rows.sort_by(|left, right| {
                    left.start_1based
                        .cmp(&right.start_1based)
                        .then(left.end_1based.cmp(&right.end_1based))
                        .then(left.feature_id.cmp(&right.feature_id))
                });
                let start_1based = rows.iter().map(|row| row.start_1based).min().unwrap_or(1);
                let end_1based = rows
                    .iter()
                    .map(|row| row.end_1based)
                    .max()
                    .unwrap_or(start_1based);
                EvidenceSvgParentGroup {
                    parent_id,
                    start_1based,
                    end_1based,
                    rows,
                }
            })
            .collect()
    }

    fn probe_region_evidence_svg_exons_for_transcript(
        report: &ProbeRegionEvidenceInterpretationReport,
        transcript_id: &str,
        projector: &EvidenceSvgProjector,
    ) -> Vec<EvidenceSvgExonSegment> {
        let mut seen = BTreeSet::new();
        let mut out = Vec::new();
        for row in &report.evidence_rows {
            for mapping in &row.transcript_mappings {
                if mapping.transcript_id != transcript_id {
                    continue;
                }
                for (idx, raw_range) in mapping.exon_ranges_1based.iter().enumerate() {
                    let Some((start, end)) = Self::probe_region_evidence_svg_parse_range(raw_range)
                    else {
                        continue;
                    };
                    let ordinal = mapping.exon_ordinals.get(idx).copied().unwrap_or(idx + 1);
                    if seen.insert((ordinal, start, end)) {
                        out.push(EvidenceSvgExonSegment {
                            ordinal,
                            start_1based: projector.project_report_coord(start.min(end)),
                            end_1based: projector.project_report_coord(start.max(end)),
                            original_range: raw_range.clone(),
                        });
                    }
                }
            }
        }
        out.sort_by(|left, right| {
            left.start_1based
                .cmp(&right.start_1based)
                .then(left.ordinal.cmp(&right.ordinal))
        });
        out
    }

    fn probe_region_evidence_svg_junctions_for_transcript(
        report: &ProbeRegionEvidenceInterpretationReport,
        transcript_id: &str,
        projector: &EvidenceSvgProjector,
    ) -> Vec<EvidenceSvgJunctionSegment> {
        let mut seen = BTreeSet::new();
        let mut out = Vec::new();
        for row in &report.evidence_rows {
            for mapping in &row.transcript_mappings {
                if mapping.transcript_id != transcript_id {
                    continue;
                }
                for span in &mapping.junction_spans {
                    if seen.insert((
                        span.from_exon_ordinal,
                        span.to_exon_ordinal,
                        span.genomic_start_1based,
                        span.genomic_end_1based,
                    )) {
                        out.push(EvidenceSvgJunctionSegment {
                            from_exon: span.from_exon_ordinal,
                            to_exon: span.to_exon_ordinal,
                            start_1based: projector.project_report_coord(span.genomic_start_1based),
                            end_1based: projector.project_report_coord(span.genomic_end_1based),
                            original_range: format!(
                                "{}..{}",
                                span.genomic_start_1based, span.genomic_end_1based
                            ),
                        });
                    }
                }
            }
        }
        out.sort_by(|left, right| {
            left.start_1based
                .cmp(&right.start_1based)
                .then(left.from_exon.cmp(&right.from_exon))
                .then(left.to_exon.cmp(&right.to_exon))
        });
        out
    }

    fn probe_region_evidence_svg_local_geometry_bounds(
        report: &ProbeRegionEvidenceInterpretationReport,
    ) -> Option<(usize, usize)> {
        let mut coords = Vec::new();
        for row in &report.evidence_rows {
            for mapping in &row.transcript_mappings {
                for raw_range in &mapping.exon_ranges_1based {
                    if let Some((start, end)) =
                        Self::probe_region_evidence_svg_parse_range(raw_range)
                    {
                        coords.push(start);
                        coords.push(end);
                    }
                }
                for span in &mapping.junction_spans {
                    coords.push(span.genomic_start_1based);
                    coords.push(span.genomic_end_1based);
                }
            }
        }
        Some((coords.iter().copied().min()?, coords.iter().copied().max()?))
    }

    fn probe_region_evidence_svg_junction_span_count(
        report: &ProbeRegionEvidenceInterpretationReport,
    ) -> usize {
        report
            .evidence_rows
            .iter()
            .flat_map(|row| &row.transcript_mappings)
            .map(|mapping| mapping.junction_spans.len())
            .sum()
    }

    fn probe_region_evidence_svg_ambiguity_tags(
        report: &ProbeRegionEvidenceInterpretationReport,
    ) -> Vec<String> {
        report
            .evidence_rows
            .iter()
            .flat_map(|row| row.ambiguity_tags.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    fn probe_region_evidence_svg_parse_range(raw: &str) -> Option<(usize, usize)> {
        let (left, right) = raw.split_once("..")?;
        let start = left.trim().parse::<usize>().ok()?;
        let end = right.trim().parse::<usize>().ok()?;
        Some((start.min(end), start.max(end)))
    }

    fn probe_region_evidence_svg_x(
        coord: usize,
        projector: &EvidenceSvgProjector,
        left: f64,
        width: f64,
    ) -> f64 {
        let span = projector
            .max_coord
            .saturating_sub(projector.min_coord)
            .max(1);
        let offset = coord.saturating_sub(projector.min_coord).min(span);
        left + width * (offset as f64 / span as f64)
    }

    fn probe_region_evidence_svg_logfc_fill(logfc: Option<f64>) -> &'static str {
        match logfc {
            Some(value) if value > 0.0 => "#dc2626",
            Some(value) if value < 0.0 => "#2563eb",
            Some(_) => "#64748b",
            None => "#94a3b8",
        }
    }

    fn probe_region_evidence_svg_id(raw: &str) -> String {
        let mut out = raw
            .chars()
            .map(|ch| {
                if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
                    ch
                } else {
                    '_'
                }
            })
            .collect::<String>();
        if out.is_empty() {
            out.push_str("unnamed");
        }
        out
    }
}

struct EvidenceSvgParentGroup {
    parent_id: String,
    start_1based: usize,
    end_1based: usize,
    rows: Vec<EvidenceSvgEvidenceRow>,
}

struct EvidenceSvgEvidenceRow {
    feature_id: String,
    intensity_source: Option<String>,
    start_1based: usize,
    end_1based: usize,
    logfc: Option<f64>,
}

struct EvidenceSvgExonSegment {
    ordinal: usize,
    start_1based: usize,
    end_1based: usize,
    original_range: String,
}

struct EvidenceSvgJunctionSegment {
    from_exon: usize,
    to_exon: usize,
    start_1based: usize,
    end_1based: usize,
    original_range: String,
}
