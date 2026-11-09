//! Headless transcript-direction pileup. All membership, statuses and paging are engine-owned.

use crate::transcript_presentation::Frame;
use serde::{Deserialize, Serialize};
use svg::{
    Node,
    node::element::{Group, Line, Rectangle, Text, Title},
};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Band {
    pub interval: (u64, u64),
    pub kind: String,
    pub label: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Row {
    pub label: String,
    pub details: String,
    pub off_axis: String,
    pub bands: Vec<Band>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Scene {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub width: u32,
    pub header: Vec<String>,
    pub model: Vec<Band>,
    pub rows: Vec<Row>,
}

fn text(x: f64, y: f64, value: &str) -> Text {
    Text::new(value)
        .set("x", x)
        .set("y", y)
        .set("font-size", 12)
        .set("font-family", "monospace")
        .set("textLength", value.chars().count() as f64 * 8.0)
        .set("lengthAdjust", "spacingAndGlyphs")
        .set("fill", "#172b34")
}

fn wrapped(value: &str, width: u32) -> Vec<String> {
    let columns = (width / 8).max(1) as usize;
    value
        .split('\n')
        .flat_map(|line| {
            let chars: Vec<_> = line.chars().collect();
            if chars.is_empty() {
                vec![String::new()]
            } else {
                chars
                    .chunks(columns)
                    .map(|chunk| chunk.iter().collect())
                    .collect()
            }
        })
        .collect()
}

fn block(x: f64, y: f64, lines: &[String]) -> Group {
    let mut group = Group::new();
    for (index, line) in lines.iter().enumerate() {
        if !line.is_empty() {
            group.append(text(x, y + index as f64 * 18.0, line));
        }
    }
    group
}

fn after_lines(top: u64, lines: usize, padding: u64) -> Result<u64, String> {
    u64::try_from(lines)
        .ok()
        .and_then(|count| count.checked_mul(18))
        .and_then(|height| top.checked_add(height))
        .and_then(|height| height.checked_add(padding))
        .ok_or_else(|| "read_end_pileup_render_height".into())
}

fn row_height(row: &Row, width: u32) -> Result<u32, String> {
    let label = wrapped(&row.label, 185).len();
    let status = wrapped(&row.off_axis, width - 230).len();
    let height = after_lines(31, label - 1, 0)?.max(after_lines(54, status - 1, 0)?);
    u32::try_from(height).map_err(|_| "read_end_pileup_render_height".into())
}

fn validate(scene: &Scene) -> Result<(), String> {
    if scene.start >= scene.end
        || !matches!(scene.strand, -1 | 1)
        || !(900..=8000).contains(&scene.width)
        || scene.rows.len() > 20_000
    {
        return Err("read_end_pileup_render_limits".into());
    }
    Ok(())
}

/// Cache these page-local heights once to keep native virtualization font-independent.
pub fn lane_heights(scene: &Scene) -> Result<Vec<u32>, String> {
    validate(scene)?;
    scene
        .rows
        .iter()
        .map(|row| row_height(row, scene.width))
        .collect()
}

fn row_group(frame: Frame, row: &Row, width: u32, top: f64) -> Group {
    let mut group = Group::new().add(Title::new(&row.details));
    group.append(block(10.0, top + 17.0, &wrapped(&row.label, 185)));
    draw(&mut group, frame, &row.bands, top + 13.0);
    group.append(block(
        205.0,
        top + 40.0,
        &wrapped(&row.off_axis, width - 230),
    ));
    group
}

struct PageLayout {
    heights: Vec<u32>,
    header_lines: Vec<String>,
    title: Vec<String>,
    legend: Vec<String>,
    title_y: u64,
    model_y: u64,
    legend_y: u64,
    rows_y: u64,
    height: u64,
}

fn page_layout(scene: &Scene) -> Result<PageLayout, String> {
    let heights = lane_heights(scene)?;
    let header_lines: Vec<_> = scene
        .header
        .iter()
        .flat_map(|line| wrapped(line, scene.width - 20))
        .collect();
    let title = wrapped(
        &format!(
            "{} {}..{} ({:+} transcript direction)",
            scene.contig, scene.start, scene.end, scene.strand
        ),
        scene.width - 20,
    );
    let legend = wrapped(
        "N = solid; other gap = dashed; segment link = dotted. Clips/evidence below axis, not aligned bases.",
        scene.width - 20,
    );
    let title_y = after_lines(40, header_lines.len(), 0)?;
    let model_y = after_lines(title_y, title.len(), 8)?;
    let legend_y = after_lines(model_y, 0, 32)?;
    let rows_y = after_lines(legend_y, legend.len(), 8)?;
    let height = heights
        .iter()
        .try_fold(after_lines(rows_y, 0, 8)?, |sum, height| {
            sum.checked_add(u64::from(*height))
                .ok_or_else(|| "read_end_pileup_render_height".to_owned())
        })?;
    Ok(PageLayout {
        heights,
        header_lines,
        title,
        legend,
        title_y,
        model_y,
        legend_y,
        rows_y,
        height,
    })
}

/// Exact intrinsic SVG pixel dimensions, using the same checked layout as `render`.
/// PNG/PDF callers must check their pixel budget before invoking a converter.
pub fn dimensions(scene: &Scene) -> Result<(u32, u64), String> {
    Ok((scene.width, page_layout(scene)?.height))
}

/// Deterministic SVG; dimensions depend only on the explicit page, never the host fonts.
pub fn render(scene: &Scene) -> Result<String, String> {
    let layout = page_layout(scene)?;
    let frame = Frame {
        start: scene.start + 1,
        end: scene.end,
        strand: scene.strand,
        left: 205.0,
        right: scene.width as f64 - 25.0,
        label_left: 10.0,
    };
    let mut document = svg::Document::new()
        .set("viewBox", (0, 0, scene.width, layout.height))
        .set("width", scene.width)
        .set("height", layout.height)
        .add(
            Rectangle::new()
                .set("width", "100%")
                .set("height", "100%")
                .set("fill", "#ffffff"),
        );
    document.append(block(10.0, 22.0, &layout.header_lines));
    document.append(block(10.0, layout.title_y as f64, &layout.title));
    draw(&mut document, frame, &scene.model, layout.model_y as f64);
    document.append(block(10.0, layout.legend_y as f64, &layout.legend));
    let mut top = layout.rows_y as f64;
    for (row, row_height) in scene.rows.iter().zip(layout.heights) {
        document.append(row_group(frame, row, scene.width, top));
        top += f64::from(row_height);
    }
    Ok(document.to_string())
}

fn draw(doc: &mut impl Node, frame: Frame, bands: &[Band], y: f64) {
    for band in bands {
        if let Some((x, width)) = frame.half_open(band.interval.0, band.interval.1) {
            match band.kind.as_str() {
                "splice_n" | "deletion_or_gap" | "segment_link" => {
                    let dash = match band.kind.as_str() {
                        "splice_n" => "none",
                        "segment_link" => "2,4",
                        _ => "8,4",
                    };
                    doc.append(
                        Line::new()
                            .set("x1", x)
                            .set("x2", x + width)
                            .set("y1", y)
                            .set("y2", y)
                            .set("stroke", "#5d6871")
                            .set("stroke-dasharray", dash)
                            .add(Title::new(format!("{}: {}", band.kind, band.label))),
                    );
                }
                kind => {
                    let (color, h) = match kind {
                        "cds" => ("#16485a", 16),
                        "stop" => ("#a54512", 20),
                        "exon" => ("#b7d5dc", 10),
                        _ => ("#2e7186", 10),
                    };
                    doc.append(
                        Rectangle::new()
                            .set("x", x)
                            .set("y", y - h as f64 / 2.0)
                            .set("width", width.max(1.0))
                            .set("height", h)
                            .set("fill", color)
                            .add(Title::new(format!("{}: {}", kind, band.label))),
                    );
                }
            }
        }
    }
}

/// Compact model row for the same horizontal frame as virtualized lanes.
fn model_scene(scene: &Scene) -> Scene {
    Scene {
        contig: scene.contig.clone(),
        start: scene.start,
        end: scene.end,
        strand: scene.strand,
        width: scene.width,
        header: vec![],
        model: vec![],
        rows: vec![Row {
            label: "Selected transcript".into(),
            details: scene.header.join("\n"),
            off_axis: format!(
                "Genomic {}:{}..{}; transcript direction {:+}. N solid; other gaps dashed; segment links dotted.",
                scene.contig, scene.start, scene.end, scene.strand
            ),
            bands: scene.model.clone(),
        }],
    }
}

/// Model row with the same wrapping and genomic frame as the virtualized lanes.
pub fn render_model(scene: &Scene) -> Result<String, String> {
    render_lane(&model_scene(scene), 0)
}

/// Height of the model row when a long contig description needs multiple lines.
pub fn model_height(scene: &Scene) -> Result<u32, String> {
    Ok(lane_heights(&model_scene(scene))?[0])
}

/// One visible sublane, using exactly the same frame and drawing code as the complete SVG.
pub fn render_lane(scene: &Scene, index: usize) -> Result<String, String> {
    validate(scene)?;
    let row = scene.rows.get(index).ok_or("read_end_pileup_lane_range")?;
    let height = row_height(row, scene.width)?;
    let frame = Frame {
        start: scene.start + 1,
        end: scene.end,
        strand: scene.strand,
        left: 205.0,
        right: scene.width as f64 - 25.0,
        label_left: 10.0,
    };
    let mut document = svg::Document::new()
        .set("viewBox", (0, 0, scene.width, height))
        .set("width", scene.width)
        .set("height", height)
        .add(
            Rectangle::new()
                .set("width", "100%")
                .set("height", "100%")
                .set("fill", "#ffffff"),
        );
    document.append(row_group(frame, row, scene.width, 0.0));
    Ok(document.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_end_page_height_arithmetic_refuses_overflow() {
        assert_eq!(after_lines(40, 2, 8).unwrap(), 84);
        assert!(after_lines(u64::MAX, 1, 0).is_err());
        assert!(after_lines(u64::MAX, 0, 1).is_err());
    }

    #[test]
    fn read_end_long_labels_wrap_without_loss_or_crossing_page_or_lane_bounds() {
        let long = "status observed 0 bp / right_censored / missing / \u{3b1} ".repeat(20);
        for width in [900, 1400, 8000] {
            let scene = Scene {
                contig: "synthetic_contig".repeat(20),
                start: 100,
                end: 200,
                strand: 1,
                width,
                header: vec![long.clone()],
                model: vec![],
                rows: vec![Row {
                    label: long.clone(),
                    details: "complete synthetic provenance".into(),
                    off_axis: long.clone(),
                    bands: vec![],
                }],
            };
            assert_eq!(wrapped(&long, width - 230).concat(), long);
            let heights = lane_heights(&scene).unwrap();
            assert!(heights[0] > 54);
            let (page_width, page_height) = dimensions(&scene).unwrap();
            let full = render(&scene).unwrap();
            let root = full.split('>').next().unwrap();
            assert!(root.contains(&format!(" width=\"{page_width}\"")));
            assert!(root.contains(&format!(" height=\"{page_height}\"")));
            for svg in [
                full,
                render_lane(&scene, 0).unwrap(),
                render_model(&scene).unwrap(),
            ] {
                let numeric = |tag: &str, name: &str| -> f64 {
                    // Leading space: `y="` must not match the tail of `font-family="`.
                    let pattern = format!(" {name}=\"");
                    tag.split(pattern.as_str())
                        .nth(1)
                        .unwrap()
                        .split('"')
                        .next()
                        .unwrap()
                        .parse()
                        .unwrap()
                };
                let probe = r#"<text data-x="not-a-number" data-y="not-a-number" font-family="monospace" x="10" y="17"/>"#;
                assert_eq!(numeric(probe, "x"), 10.0);
                assert_eq!(numeric(probe, "y"), 17.0);
                let root = svg.split('>').next().unwrap();
                let height = numeric(root, "height");
                for text in svg.split("<text ").skip(1) {
                    let tag = text.split('>').next().unwrap();
                    assert!(
                        numeric(tag, "x") + numeric(tag, "textLength") <= f64::from(width) - 10.0
                    );
                    assert!(numeric(tag, "y") <= height - 8.0);
                }
                assert_eq!(svg.matches("<svg").count(), 1);
            }
            let lane = render_lane(&scene, 0).unwrap();
            let shown: String = lane
                .split("<text ")
                .skip(1)
                .map(|part| {
                    part.split_once('>')
                        .unwrap()
                        .1
                        .split("</text>")
                        .next()
                        .unwrap()
                        .trim_matches('\n')
                })
                .collect();
            assert_eq!(shown, format!("{long}{long}"));
            assert!(lane.contains("complete synthetic provenance"));
            assert_eq!(lane, render_lane(&scene, 0).unwrap());
        }
    }

    #[test]
    fn read_end_virtual_lane_uses_the_same_block_frame_without_rendering_hidden_lanes() {
        let row = Row {
            label: "synthetic_unit".into(),
            details: "synthetic provenance".into(),
            off_axis: "observed: 0 bp; not missing".into(),
            bands: vec![Band {
                interval: (110, 130),
                kind: "aligned".into(),
                label: "saved block".into(),
            }],
        };
        for strand in [1, -1] {
            let scene = Scene {
                contig: "synthetic".into(),
                start: 100,
                end: 200,
                strand,
                width: 1000,
                header: vec![],
                model: vec![],
                rows: vec![row.clone()],
            };
            let full = render(&scene).unwrap();
            let lane = render_lane(&scene, 0).unwrap();
            assert_eq!(lane_heights(&scene).unwrap(), vec![54]);
            let x = |text: &str| {
                text.split("<rect ")
                    .skip(1)
                    .filter_map(|r| {
                        r.split(" x=\"")
                            .nth(1)
                            .map(|v| v.split('"').next().unwrap().to_string())
                    })
                    .collect::<Vec<_>>()
            };
            assert_eq!(x(&full), x(&lane));
            assert!(!x(&lane).is_empty());
            assert!(lane.contains("observed: 0 bp"));
            assert!(render_lane(&scene, 1).is_err());
        }
    }
    #[test]
    fn read_end_pileup_svg_is_deterministic_and_distinguishes_connectors() {
        let mut scene = Scene {
            contig: "synthetic".into(),
            start: 100,
            end: 200,
            strand: -1,
            width: 1000,
            header: vec!["0 observed is not missing".into()],
            model: vec![],
            rows: vec![Row {
                label: "synthetic_unit".into(),
                details: "no_bound_basecalls".into(),
                off_axis: "hard clip 5 bp".into(),
                bands: ["splice_n", "deletion_or_gap", "segment_link"]
                    .into_iter()
                    .map(|kind| Band {
                        interval: (110, 150),
                        kind: kind.into(),
                        label: kind.into(),
                    })
                    .collect(),
            }],
        };
        let a = render(&scene).unwrap();
        assert_eq!(a, render(&scene).unwrap());
        for needle in [
            "splice_n",
            "deletion_or_gap",
            "segment_link",
            "no_bound_basecalls",
            "hard clip",
        ] {
            assert!(a.contains(needle));
        }
        scene.strand = 1;
        assert_ne!(a, render(&scene).unwrap());
        scene.width = 1;
        assert!(render(&scene).is_err());
    }
}
