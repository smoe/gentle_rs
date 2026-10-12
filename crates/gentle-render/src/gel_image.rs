//! Portable measured-gel SVG presentation. No band sizing is performed here.

use gentle_protocol::gel_image::{GelImageAnalysisReport, GelMigrationDirection};
use svg::node::element::{Circle, Image, Line, Rectangle, Text};

fn text(x: f64, y: f64, value: impl Into<String>, size: u32) -> Text {
    Text::new(value.into())
        .set("x", x)
        .set("y", y)
        .set("font-size", size)
}

fn size_label(value: f64) -> String {
    if value >= 100.0 {
        format!("{value:.0}")
    } else {
        format!("{value:.1}")
    }
}

fn wrap(value: &str, width: usize) -> Vec<String> {
    let mut lines = vec![];
    let mut current = String::new();
    for word in value.split_whitespace() {
        if !current.is_empty() && current.len() + word.len() + 1 > width {
            lines.push(std::mem::take(&mut current));
        }
        if !current.is_empty() {
            current.push(' ');
        }
        current.push_str(word);
    }
    if !current.is_empty() {
        lines.push(current);
    }
    lines
}

/// Render a validated report over its bounded PNG preview. The original image
/// aspect ratio and pixel orientation are preserved; all labels are escaped.
pub fn render_gel_image_analysis_svg(
    report: &GelImageAnalysisReport,
    preview_png_base64: &str,
) -> String {
    let scale = (900.0 / f64::from(report.image.width)).min(600.0 / f64::from(report.image.height));
    let left = 120.0;
    let top = 170.0;
    let width = f64::from(report.image.width) * scale;
    let height = f64::from(report.image.height) * scale;
    let mut document = svg::Document::new()
        .set("font-family", "sans-serif")
        .set("fill", "#23313d");
    document = document.add(text(28.0, 32.0, "GENtle | Measured gel band sizing", 22));
    document = document.add(text(
        28.0,
        60.0,
        format!(
            "{} | {} | {}",
            report.request.report_id,
            report.image.source_name,
            report.request.size_kind.unit()
        ),
        14,
    ));
    document = document.add(text(
        28.0,
        84.0,
        format!(
            "{} | original {} x {} px",
            report.algorithm, report.image.width, report.image.height
        ),
        13,
    ));
    document = document.add(text(28.0, 107.0, "Orange: confirmed ladder | Teal: measured sample | Dashed: user-defined lane | Numbers: table rows", 13));
    document = document.add(text(28.0, 130.0, "Display preview only. No image warping, intensity quantification or automatic band identification.", 12));
    document = document.add(
        Image::new()
            .set("x", left)
            .set("y", top)
            .set("width", width)
            .set("height", height)
            .set(
                "href",
                format!("data:image/png;base64,{preview_png_base64}"),
            ),
    );
    for lane in &report.request.lanes {
        document = document.add(
            Rectangle::new()
                .set("x", left + lane.min.x * scale)
                .set("y", top + lane.min.y * scale)
                .set("width", (lane.max.x - lane.min.x) * scale)
                .set("height", (lane.max.y - lane.min.y) * scale)
                .set("fill", "none")
                .set("stroke", "#d7dfdf")
                .set("stroke-dasharray", "4,4"),
        );
    }
    let vertical = matches!(
        report.request.migration,
        GelMigrationDirection::Down | GelMigrationDirection::Up
    );
    let mut bands = report.calibration_bands.iter().collect::<Vec<_>>();
    bands.sort_by(|a, b| {
        let coordinate = |band: &&gentle_protocol::gel_image::GelLadderBand| {
            if vertical {
                band.center.y
            } else {
                band.center.x
            }
        };
        coordinate(a).total_cmp(&coordinate(b))
    });
    let mut last_label = f64::NEG_INFINITY;
    for band in bands {
        let x = left + band.center.x * scale;
        let y = top + band.center.y * scale;
        document = document.add(
            Circle::new()
                .set("cx", x)
                .set("cy", y)
                .set("r", 4)
                .set("fill", "#df861b"),
        );
        let coord = if vertical { y } else { x };
        if coord - last_label >= if vertical { 18.0 } else { 72.0 } {
            let label = format!(
                "{} {}",
                size_label(band.size),
                report.request.size_kind.unit()
            );
            if vertical {
                document = document.add(
                    Line::new()
                        .set("x1", left - 8.0)
                        .set("x2", left)
                        .set("y1", y)
                        .set("y2", y)
                        .set("stroke", "#23313d"),
                );
                document =
                    document.add(text(left - 14.0, y + 4.0, label, 12).set("text-anchor", "end"));
            } else {
                document = document.add(
                    Line::new()
                        .set("x1", x)
                        .set("x2", x)
                        .set("y1", top - 8.0)
                        .set("y2", top)
                        .set("stroke", "#23313d"),
                );
                document =
                    document.add(text(x, top - 15.0, label, 12).set("text-anchor", "middle"));
            }
            last_label = coord;
        }
    }
    for (index, band) in report.estimates.iter().enumerate() {
        let x = left + band.center.x * scale;
        let y = top + band.center.y * scale;
        document = document.add(
            Circle::new()
                .set("cx", x)
                .set("cy", y)
                .set("r", 4)
                .set("fill", "#16a19c"),
        );
        document = document
            .add(text(x + 7.0, y - 5.0, (index + 1).to_string(), 12).set("fill", "#16a19c"));
    }
    let mut y = top + height + 32.0;
    for line in wrap(
        &format!(
            "Ladder: {} | {} | system: {}",
            report.request.ladder.label,
            report.request.ladder.source,
            report
                .request
                .ladder
                .gel_system
                .as_deref()
                .unwrap_or("not specified")
        ),
        125,
    ) {
        document = document.add(text(28.0, y, line, 13));
        y += 18.0;
    }
    document = document.add(text(
        28.0,
        y,
        "Estimates (rounded for display; exact values and calibration assignments in JSON/TSV)",
        14,
    ));
    y += 26.0;
    for (index, band) in report.estimates.iter().enumerate() {
        let estimate = band
            .estimated_size
            .map(|size| format!("~{} {}", size_label(size), report.request.size_kind.unit()))
            .unwrap_or_else(|| "outside calibrated range (not extrapolated)".into());
        let row = format!(
            "{}. {} | {} | {} | {}",
            index + 1,
            band.lane_id,
            band.label,
            estimate,
            band.warnings.join("; ")
        );
        for line in wrap(&row, 120) {
            document = document.add(text(28.0, y, line, 13));
            y += 18.0;
        }
    }
    y += 10.0;
    document = document.add(text(28.0, y, "Calibration limits and interpretation", 14));
    y += 22.0;
    for warning in &report.warnings {
        for line in wrap(warning, 130) {
            document = document.add(text(28.0, y, line, 12));
            y += 17.0;
        }
    }
    y += 10.0;
    document = document.add(text(
        28.0,
        y,
        format!("Original SHA-256: {}", report.image.sha256),
        11,
    ));
    document
        .set("width", 1100)
        .set("height", y + 26.0)
        .set("viewBox", (0, 0, 1100, y + 26.0))
        .to_string()
}
