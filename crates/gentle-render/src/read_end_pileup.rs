//! Headless transcript-direction pileup. All membership, statuses and paging are engine-owned.

use crate::transcript_presentation::Frame;
use serde::{Deserialize, Serialize};
use svg::{
    Node,
    node::element::{Line, Rectangle, Text, Title},
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
        .set("fill", "#172b34")
}

/// Deterministic SVG; dimensions depend only on the explicit page, never the host fonts.
pub fn render(scene: &Scene) -> Result<String, String> {
    if scene.start >= scene.end
        || !matches!(scene.strand, -1 | 1)
        || !(900..=8000).contains(&scene.width)
        || scene.rows.len() > 20_000
    {
        return Err("read_end_pileup_render_limits".into());
    }
    let header = 30.0 + 19.0 * scene.header.len() as f64;
    let height = header + 95.0 + scene.rows.len() as f64 * 54.0;
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
    for (i, line) in scene.header.iter().enumerate() {
        document.append(text(10.0, 22.0 + i as f64 * 19.0, line));
    }
    document.append(text(
        10.0,
        header + 15.0,
        &format!(
            "{} {}..{} ({:+} transcript direction)",
            scene.contig, scene.start, scene.end, scene.strand
        ),
    ));
    draw(&mut document, frame, &scene.model, header + 36.0);
    document.append(text(10.0, header + 65.0, "N = solid; other gap = dashed; segment link = dotted. Clips/evidence below axis, not aligned bases."));
    for (i, row) in scene.rows.iter().enumerate() {
        let y = header + 93.0 + i as f64 * 54.0;
        document.append(text(10.0, y, &row.label).add(Title::new(&row.details)));
        draw(&mut document, frame, &row.bands, y - 4.0);
        document.append(text(205.0, y + 19.0, &row.off_axis).add(Title::new(&row.details)));
    }
    Ok(document.to_string())
}

fn draw(doc: &mut svg::Document, frame: Frame, bands: &[Band], y: f64) {
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
pub fn render_model(scene: &Scene) -> Result<String, String> {
    let model = Scene {
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
    };
    render_lane(&model, 0)
}

/// One visible sublane, using exactly the same frame and drawing code as the complete SVG.
pub fn render_lane(scene: &Scene, index: usize) -> Result<String, String> {
    if scene.start >= scene.end
        || !matches!(scene.strand, -1 | 1)
        || !(900..=8000).contains(&scene.width)
    {
        return Err("read_end_pileup_render_limits".into());
    }
    let row = scene.rows.get(index).ok_or("read_end_pileup_lane_range")?;
    let frame = Frame {
        start: scene.start + 1,
        end: scene.end,
        strand: scene.strand,
        left: 205.0,
        right: scene.width as f64 - 25.0,
        label_left: 10.0,
    };
    let mut document = svg::Document::new()
        .set("viewBox", (0, 0, scene.width, 54))
        .set("width", scene.width)
        .set("height", 54)
        .add(
            Rectangle::new()
                .set("width", "100%")
                .set("height", "100%")
                .set("fill", "#ffffff"),
        );
    document.append(text(10.0, 17.0, &row.label).add(Title::new(&row.details)));
    draw(&mut document, frame, &row.bands, 13.0);
    document.append(text(205.0, 40.0, &row.off_axis).add(Title::new(&row.details)));
    Ok(document.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
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
