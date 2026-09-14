//! Source-coherent drawing only: identities and agreement are engine-owned.

use gentle_protocol::transcript_presentation::*;
use std::collections::BTreeSet;
use svg::{
    Node,
    node::element::{Circle, Group, Line, Rectangle, Text, Title},
};

const PITCH: f64 = 26.0;

/// A genomic frame, always oriented in the transcript direction.
#[derive(Clone, Copy)]
pub struct Frame {
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub left: f64,
    pub right: f64,
    pub label_left: f64,
}
impl Frame {
    fn x(self, g: u64) -> f64 {
        let offset = if self.strand == -1 {
            self.end.saturating_sub(g)
        } else {
            g.saturating_sub(self.start)
        };
        if self.start == self.end {
            return (self.left + self.right) * 0.5;
        }
        self.left + offset as f64 / (self.end - self.start) as f64 * (self.right - self.left)
    }
    fn interval(self, i: TranscriptInterval) -> Option<(f64, f64)> {
        let a = i.start_1based.max(self.start);
        let b = i.end_1based.min(self.end);
        if a > b {
            return None;
        }
        if self.start == self.end {
            return Some((self.left, self.right - self.left));
        }
        let (x, y) = (self.x(a), self.x(b));
        let half_base = (self.right - self.left) / (self.end - self.start) as f64 * 0.5;
        let left = (x.min(y) - half_base).max(self.left);
        let right = (x.max(y) + half_base).min(self.right);
        Some((left, (right - left).max(0.0)))
    }
}

fn lanes(report: &TranscriptStructurePresentation, frame: Frame) -> (Vec<(usize, usize)>, usize) {
    let mut ends: Vec<u64> = Vec::new();
    let mut rows = Vec::new();
    for (index, e) in report
        .physical_exons
        .iter()
        .enumerate()
        .filter(|(_, e)| e.strand == frame.strand && frame.interval(e.interval).is_some())
    {
        let lane = ends
            .iter()
            .position(|end| *end < e.interval.start_1based)
            .unwrap_or(ends.len());
        if lane == ends.len() {
            ends.push(e.interval.end_1based);
        } else {
            ends[lane] = e.interval.end_1based;
        }
        rows.push((index, lane));
    }
    (rows, ends.len().max(1))
}

fn cds_groups(report: &TranscriptStructurePresentation) -> Vec<&TranscriptStructureGroup> {
    let mut ids = BTreeSet::new();
    report
        .structure_groups
        .iter()
        .filter(|g| g.cds.as_ref().is_some_and(|c| !c.is_empty()) && ids.insert(&g.cds_geometry_id))
        .collect()
}

pub fn height(report: &TranscriptStructurePresentation, frame: Frame) -> f64 {
    (lanes(report, frame).1
        + report.structure_groups.len()
        + cds_groups(report).len()
        + 5
        + delta_lines(report, frame).len()) as f64
        * PITCH
}

fn delta_lines(report: &TranscriptStructurePresentation, frame: Frame) -> Vec<String> {
    let deltas: Vec<_> = report
        .tss_deltas
        .iter()
        .filter(|d| {
            d.strand == frame.strand
                && ((frame.start..=frame.end).contains(&d.ensembl_position_1based)
                    || (frame.start..=frame.end).contains(&d.refseq_position_1based))
        })
        .collect();
    let mut lines = vec![format!(
        "{} cross-source start comparisons touching this view; {} exact. All pairings, not a nearest-start assignment.",
        deltas.len(),
        deltas
            .iter()
            .filter(|d| d.transcript_oriented_delta_bp == 0)
            .count()
    )];
    for d in deltas.iter().take(16) {
        lines.push(format!(
            "Ensembl {} -> RefSeq {}: {:+} bp (transcript-oriented)",
            d.ensembl_position_1based, d.refseq_position_1based, d.transcript_oriented_delta_bp
        ));
    }
    if deltas.len() > 16 {
        lines.push(
            "First 16 comparisons printed; complete deltas remain in JSON and SVG hover.".into(),
        );
    }
    lines
}

fn details(report: &TranscriptStructurePresentation, ids: &[String]) -> String {
    ids.iter()
        .filter_map(|id| report.records.iter().find(|r| &r.record_id == id))
        .map(|r| {
            let s = report
                .sources
                .iter()
                .find(|s| s.source_id == r.structure.source_id);
            format!(
                "{} | {} | designations: {} | record SHA-256 {} | {} | source exons: {} | {}",
                r.structure.transcript_id,
                r.structure.label,
                r.structure
                    .designations
                    .iter()
                    .map(|d| format!("{} ({}={})", d.label, d.field, d.value))
                    .collect::<Vec<_>>()
                    .join("; "),
                r.content_sha256,
                s.map(|s| format!(
                    "{:?} {} / {} / {} / reference {} / annotation SHA-256 {}",
                    s.provider,
                    s.assembly,
                    s.release,
                    s.accession,
                    s.chromosome,
                    s.annotation_sha256
                ))
                .unwrap_or_default(),
                r.structure
                    .exons
                    .iter()
                    .map(|e| format!(
                        "{} ({}..{})",
                        e.source_exon_id.as_deref().unwrap_or("ID not supplied"),
                        e.interval.start_1based,
                        e.interval.end_1based
                    ))
                    .collect::<Vec<_>>()
                    .join("; "),
                r.structure.notes.join("; ")
            )
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn label(text: impl Into<String>, x: f64, y: f64) -> Text {
    Text::new(text.into())
        .set("x", x)
        .set("y", y + 4.0)
        .set("font-family", "sans-serif")
        .set("font-size", 11)
        .set("fill", "#29323a")
}

/// Use the same genomic transform in both overview and detail views. No clipped
/// interval is ever used to decide structure equivalence.
pub fn render(report: &TranscriptStructurePresentation, frame: Frame, top: f64) -> Group {
    let mut group = Group::new()
        .set("data-role", "source-coherent-transcripts")
        .set("data-content-sha256", report.content_sha256.as_str());
    if let Err(error) = gentle_engine::transcript_presentation::validate(report) {
        group.append(label(
            format!("Invalid transcript source binding: {error}"),
            frame.label_left,
            top,
        ));
        return group;
    }
    let mut y = top;
    for (provider, name, color) in [
        (TranscriptProvider::Ensembl, "Ensembl TSSs", "#176b9a"),
        (TranscriptProvider::RefSeq, "RefSeq TSSs", "#b13b43"),
    ] {
        let supplied = report.sources.iter().any(|s| s.provider == provider);
        group.append(label(
            format!("{name}{}", if supplied { "" } else { " (unassessed)" }),
            frame.label_left,
            y,
        ));
        for tick in report.tss_ticks.iter().filter(|t| {
            t.provider == provider
                && t.strand == frame.strand
                && (frame.start..=frame.end).contains(&t.genomic_position_1based)
        }) {
            let x = frame.x(tick.genomic_position_1based);
            let deltas = report
                .tss_deltas
                .iter()
                .filter(|d| {
                    d.strand == tick.strand
                        && match provider {
                            TranscriptProvider::Ensembl => {
                                d.ensembl_position_1based == tick.genomic_position_1based
                            }
                            TranscriptProvider::RefSeq => {
                                d.refseq_position_1based == tick.genomic_position_1based
                            }
                        }
                })
                .map(|d| {
                    format!(
                        "Ensembl {} -> RefSeq {}: {:+} bp (transcript-oriented)",
                        d.ensembl_position_1based,
                        d.refseq_position_1based,
                        d.transcript_oriented_delta_bp
                    )
                })
                .collect::<Vec<_>>()
                .join("\n");
            let agreement = if tick.exact_cross_source_agreement {
                "exact agreement"
            } else if report.sources.iter().any(|s| s.provider != provider) {
                "source-only start in supplied annotations"
            } else {
                "other source unassessed"
            };
            let mut marker = Group::new()
                .set("data-role", "source-tss")
                .set("data-provider", format!("{provider:?}"))
                .set("data-genomic-position", tick.genomic_position_1based)
                .set("data-agreement", agreement);
            marker.append(Title::new(format!(
                "{name}: g {} | {agreement}\n{deltas}\n{}",
                tick.genomic_position_1based,
                details(report, &tick.member_record_ids)
            )));
            marker.append(
                Line::new()
                    .set("x1", x)
                    .set("x2", x)
                    .set("y1", y - 9.0)
                    .set("y2", y + 9.0)
                    .set("stroke", color)
                    .set("stroke-width", 2),
            );
            marker.append(
                Text::new(if tick.exact_cross_source_agreement {
                    "="
                } else if report.sources.iter().any(|s| s.provider != provider) {
                    "s"
                } else {
                    "?"
                })
                .set("x", x + 3.0)
                .set("y", y - 4.0)
                .set("font-size", 10)
                .set("fill", color),
            );
            group.append(marker);
        }
        y += PITCH;
    }
    group.append(label("= exact start; s source-only in supplied annotations; ? other source unassessed. No preferred TSS inferred.",frame.label_left,y));
    y += PITCH;
    for line in delta_lines(report, frame) {
        group.append(label(line, frame.label_left, y));
        y += PITCH;
    }
    let (exon_rows, nlanes) = lanes(report, frame);
    group.append(label("Physical exons (shared)", frame.label_left, y));
    for (index, lane) in exon_rows {
        let e = &report.physical_exons[index];
        let (x, w) = frame.interval(e.interval).expect("visible exon");
        let mut cell = Group::new()
            .set("data-role", "physical-exon")
            .set("data-exon-id", e.exon_id.as_str());
        cell.append(Title::new(format!(
            "E{}: {}..{} strand {}\n{}",
            index + 1,
            e.interval.start_1based,
            e.interval.end_1based,
            e.strand,
            details(report, &e.member_record_ids)
        )));
        cell.append(
            Rectangle::new()
                .set("x", x)
                .set("y", y + lane as f64 * PITCH - 7.0)
                .set("width", w)
                .set("height", 14)
                .set("fill", "#d9e9df")
                .set("stroke", "#28724c"),
        );
        if w > 20.0 {
            cell.append(label(
                format!("E{}", index + 1),
                x + 2.0,
                y + lane as f64 * PITCH,
            ));
        }
        group.append(cell);
    }
    y += nlanes as f64 * PITCH;
    let cds = cds_groups(report);
    for (i, g) in cds.iter().enumerate() {
        group.append(label(
            format!("CDS geometry C{}", i + 1),
            frame.label_left,
            y,
        ));
        let members: Vec<_> = report
            .structure_groups
            .iter()
            .filter(|s| s.cds_geometry_id == g.cds_geometry_id)
            .flat_map(|s| s.member_record_ids.clone())
            .collect();
        let mut row = Group::new()
            .set("data-role", "shared-cds")
            .set("data-cds-id", g.cds_geometry_id.as_str());
        row.append(Title::new(details(report, &members)));
        for c in g.cds.iter().flatten() {
            if let Some((x, w)) = frame.interval(c.interval) {
                row.append(
                    Rectangle::new()
                        .set("x", x)
                        .set("width", w)
                        .set("y", y - 5.0)
                        .set("height", 10)
                        .set("fill", "#353f78"),
                );
            }
        }
        group.append(row);
        y += PITCH;
    }
    for (i, g) in report.structure_groups.iter().enumerate() {
        let mut row = Group::new()
            .set("data-role", "transcript-structure")
            .set("data-structure-id", g.structure_id.as_str());
        let exon_indices: Vec<_> = g
            .exon_ids
            .iter()
            .filter_map(|id| report.physical_exons.iter().position(|e| &e.exon_id == id))
            .collect();
        let chain = exon_indices
            .iter()
            .map(|i| format!("E{}", i + 1))
            .collect::<Vec<_>>()
            .join(" > ");
        let cds_label = cds
            .iter()
            .position(|c| c.cds_geometry_id == g.cds_geometry_id)
            .map(|n| format!("C{}", n + 1))
            .unwrap_or_else(|| {
                if g.cds.is_none() {
                    "CDS unassessed".into()
                } else {
                    "noncoding".into()
                }
            });
        row.append(Title::new(format!(
            "{chain}; {cds_label}\n{}",
            details(report, &g.member_record_ids)
        )));
        row.append(label(
            format!(
                "S{}: {} records / {}",
                i + 1,
                g.member_record_ids.len(),
                cds_label
            ),
            frame.label_left,
            y,
        ));
        let positions: Vec<_> = exon_indices
            .iter()
            .filter_map(|i| {
                frame
                    .interval(report.physical_exons[*i].interval)
                    .map(|(x, w)| (*i, x + w / 2.0))
            })
            .collect();
        if let (Some(a), Some(b)) = (positions.first(), positions.last()) {
            row.append(
                Line::new()
                    .set("x1", a.1)
                    .set("x2", b.1)
                    .set("y1", y)
                    .set("y2", y)
                    .set("stroke", "#87948e"),
            );
        }
        for (idx, x) in positions {
            row.append(
                Circle::new()
                    .set("cx", x)
                    .set("cy", y)
                    .set("r", 3)
                    .set("fill", "#28724c"),
            );
            row.append(
                Text::new(format!("E{}", idx + 1))
                    .set("x", x + 4.0)
                    .set("y", y - 4.0)
                    .set("font-size", 9),
            );
        }
        group.append(row);
        y += PITCH;
    }
    group.append(label("Exon boxes and CDS geometries drawn once; structure rows reference full chains, including off-window exons (hover).",frame.label_left,y));
    group
}

#[cfg(test)]
mod tests {
    use super::*;
    // Publication-safe hand-crafted intervals. This fixture is regenerated in
    // memory, not copied from a gene or a public/private annotation dataset.
    fn fixture(strand: i8) -> TranscriptStructurePresentation {
        let sources = vec![
            TranscriptSourceBinding {
                source_id: "e".into(),
                provider: TranscriptProvider::Ensembl,
                assembly: "test".into(),
                release: "e-release".into(),
                accession: "e-accession".into(),
                chromosome: "1".into(),
                annotation_sha256: "a".repeat(64),
                locus_sequence_sha256: "b".repeat(64),
            },
            TranscriptSourceBinding {
                source_id: "r".into(),
                provider: TranscriptProvider::RefSeq,
                assembly: "test".into(),
                release: "r-release".into(),
                accession: "r-accession".into(),
                chromosome: "1".into(),
                annotation_sha256: "c".repeat(64),
                locus_sequence_sha256: "b".repeat(64),
            },
        ];
        let records = sources
            .iter()
            .map(|s| {
                let mut exons = vec![
                    TranscriptExon {
                        interval: TranscriptInterval {
                            start_1based: 200,
                            end_1based: 240,
                        },
                        source_exon_id: None,
                    },
                    TranscriptExon {
                        interval: TranscriptInterval {
                            start_1based: 500,
                            end_1based: 700,
                        },
                        source_exon_id: None,
                    },
                ];
                if strand == -1 {
                    exons.reverse();
                }
                SourceTranscriptStructure {
                    source_id: s.source_id.clone(),
                    transcript_id: format!("{}.1", s.source_id),
                    label: String::new(),
                    strand,
                    exons,
                    cds: Some(vec![TranscriptCds {
                        interval: TranscriptInterval {
                            start_1based: 520,
                            end_1based: 680,
                        },
                        phase: None,
                    }]),
                    designations: vec![TranscriptDesignation {
                        label: format!("{} source designation", s.source_id),
                        field: "tag".into(),
                        value: "provider-value".into(),
                    }],
                    notes: vec![],
                }
            })
            .collect();
        gentle_engine::transcript_presentation::build(
            "test",
            "1",
            &"b".repeat(64),
            sources,
            records,
        )
        .unwrap()
    }
    #[test]
    fn shared_geometry_once_with_full_hover_provenance_on_plus_and_minus() {
        for strand in [1, -1] {
            let p = fixture(strand);
            let frame = Frame {
                start: 100,
                end: 800,
                strand,
                left: 255.0,
                right: 1050.0,
                label_left: 34.0,
            };
            let svg = render(&p, frame, 10.0).to_string();
            assert_eq!(svg.matches("data-role=\"physical-exon\"").count(), 2);
            assert_eq!(svg.matches("data-role=\"transcript-structure\"").count(), 1);
            assert_eq!(svg.matches("data-role=\"shared-cds\"").count(), 1);
            assert_eq!(svg.matches("data-role=\"source-tss\"").count(), 2);
            for text in [
                "e.1",
                "r.1",
                "e-release",
                "r-release",
                "source designation",
                "exact agreement",
                "<title>",
            ] {
                assert!(svg.contains(text), "{text}");
            }
            assert!(svg.contains(&p.sources[0].annotation_sha256));
            assert!(svg.contains(&p.content_sha256));
        }
    }
    #[test]
    fn shared_renderer_uses_existing_tss_base_center_axis_including_minus() {
        let frame = Frame {
            start: 100,
            end: 800,
            strand: 1,
            left: 255.0,
            right: 1050.0,
            label_left: 34.0,
        };
        assert_eq!(frame.x(100), 255.0);
        assert_eq!(frame.x(800), 1050.0);
        assert_eq!(frame.x(450), 652.5);
        let reverse = Frame {
            strand: -1,
            ..frame
        };
        assert_eq!(reverse.x(800), 255.0);
        assert_eq!(reverse.x(100), 1050.0);
        let (x, w) = frame
            .interval(TranscriptInterval {
                start_1based: 200,
                end_1based: 240,
            })
            .unwrap();
        assert!((x - (255.0 + 99.5 / 700.0 * 795.0)).abs() < 1e-9);
        assert!((w - 41.0 / 700.0 * 795.0).abs() < 1e-9);
    }
    #[test]
    fn cropping_does_not_merge_distinct_full_chains_and_tampered_content_is_not_drawn() {
        let p = fixture(1);
        let mut records: Vec<_> = p.records.iter().map(|r| r.structure.clone()).collect();
        records[1].exons[1].interval.end_1based = 750;
        let p = gentle_engine::transcript_presentation::build(
            "test",
            "1",
            &p.locus_sequence_sha256,
            p.sources,
            records,
        )
        .unwrap();
        let frame = Frame {
            start: 180,
            end: 260,
            strand: 1,
            left: 255.0,
            right: 1050.0,
            label_left: 34.0,
        };
        let svg = render(&p, frame, 0.0).to_string();
        assert_eq!(svg.matches("data-role=\"physical-exon\"").count(), 1);
        assert_eq!(svg.matches("data-role=\"transcript-structure\"").count(), 2);
        let mut altered = p;
        altered.physical_exons.clear();
        let svg = render(&altered, frame, 0.0).to_string();
        assert!(svg.contains("Invalid transcript source binding"));
        assert!(!svg.contains("data-role=\"physical-exon\""));
    }

    #[test]
    fn source_coherent_preview_svgs() {
        for strand in [1, -1] {
            let p = fixture(strand);
            let mut records: Vec<_> = p.records.iter().map(|r| r.structure.clone()).collect();
            if strand == 1 {
                records[1].exons[0].interval.start_1based += 10;
            } else {
                records[1].exons[0].interval.end_1based += 10;
            }
            let mut alternative = records[1].clone();
            alternative.transcript_id = "r.2".into();
            alternative.cds.as_mut().unwrap()[0].interval.end_1based -= 20;
            records.push(alternative);
            let p = gentle_engine::transcript_presentation::build(
                "test",
                "1",
                &p.locus_sequence_sha256,
                p.sources,
                records,
            )
            .unwrap();
            let frame = Frame {
                start: 100,
                end: 800,
                strand,
                left: 255.0,
                right: 1050.0,
                label_left: 34.0,
            };
            let h = height(&p, frame) + 60.0;
            let svg = svg::Document::new()
                .set("viewBox", (0, 0, 1400, h))
                .set("width", 1400)
                .set("height", h)
                .add(
                    Rectangle::new()
                        .set("width", 1400)
                        .set("height", h)
                        .set("fill", "white"),
                )
                .add(render(&p, frame, 24.0))
                .to_string();
            assert_eq!(svg.matches("data-role=\"physical-exon\"").count(), 3);
            assert_eq!(svg.matches("data-role=\"transcript-structure\"").count(), 3);
            assert!(svg.contains(if strand == 1 { "+10 bp" } else { "-10 bp" }));
            if let Some(dir) = std::env::var_os("GENTLE_TEST_TRANSCRIPT_SVG_DIR") {
                let dir = std::path::Path::new(&dir);
                std::fs::create_dir_all(dir).unwrap();
                std::fs::write(
                    dir.join(if strand == 1 { "plus.svg" } else { "minus.svg" }),
                    svg,
                )
                .unwrap();
            }
        }
    }
}
