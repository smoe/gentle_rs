use super::*;

impl GentleEngine {
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

    pub(super) fn probe_region_plot_data_from_table(
        path: &Path,
    ) -> Result<ProbeRegionPlotData, EngineError> {
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

    pub(super) fn probe_region_plot_track_label(header: &str) -> String {
        header
            .strip_prefix("mean_log2_")
            .or_else(|| header.strip_prefix("log2FC_"))
            .unwrap_or(header)
            .replace('_', " ")
    }

    pub(super) fn probe_region_parse_plot_value(value: &str) -> Option<f64> {
        let parsed = value.trim().parse::<f64>().ok()?;
        parsed.is_finite().then_some(parsed)
    }

    pub(super) fn probe_region_plot_range(
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

    pub(super) fn probe_region_plot_x(
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

    pub(super) fn probe_region_plot_y(
        value: f64,
        min_value: f64,
        max_value: f64,
        top: f64,
        height: f64,
    ) -> f64 {
        let frac = ((value - min_value) / (max_value - min_value)).clamp(0.0, 1.0);
        top + height - (height * frac)
    }

    pub(crate) fn probe_region_svg_escape(raw: &str) -> String {
        raw.replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
            .replace('"', "&quot;")
    }

    pub(super) fn render_probe_region_output_svg_text(
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

    pub(super) fn render_probe_region_svg_panel(
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
}
