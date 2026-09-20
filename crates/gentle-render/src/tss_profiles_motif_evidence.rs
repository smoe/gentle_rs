//! Separate sparse-package lanes: triangle footprint = covered bases, not a peak
//! width; orientation = displayed strand; height = labelled imported-score range.

use super::*;
use gentle_protocol::tss_motif_evidence::{self as evidence, ProjectedMotifHit};
use gentle_protocol::tss_profiles::TssImportedMotifEvidence;
use svg::node::element::{ClipPath, Polygon};

const LABEL_LIMIT: usize = 12;
const HALF_HEIGHT: f64 = 44.0;

pub(super) struct Lanes<'a> {
    lanes: Vec<Lane<'a>>,
    missing: bool,
}

struct Lane<'a> {
    source: &'a TssImportedMotifEvidence,
    hits: Vec<ProjectedMotifHit<'a>>,
    summary: TextBlock,
    labels: Vec<TextBlock>,
    min: f64,
    max: f64,
    height: f64,
}

fn hit_label(projected: &ProjectedMotifHit<'_>, index: usize) -> String {
    let hit = projected.hit;
    format!(
        "{}: {}:{}..{} ({} bp), local {} / genomic {}, score {}{}",
        index + 1,
        hit.chromosome,
        hit.start_0based + 1,
        hit.end_0based_exclusive,
        hit.end_0based_exclusive - hit.start_0based,
        projected.local_strand.as_str(),
        hit.strand,
        number(hit.score),
        if projected.clipped {
            "; footprint clipped at window edge"
        } else {
            ""
        }
    )
}

impl<'a> Lanes<'a> {
    pub(super) fn new(
        report: &'a TssProfileReport,
        window: &TssProfileWindow,
        accession: &str,
    ) -> Self {
        let g = &window.record.geometry;
        let mut lanes = vec![];
        for source in &report.imported_motif_evidence {
            let r = &source.report;
            if !r.request.motif_ids.iter().any(|m| m == accession)
                || !r.regions.iter().any(|q| {
                    q.resolved_chromosome.as_deref() == Some(g.chromosome.as_str())
                        && q.end_0based_exclusive >= g.start_1based
                        && q.start_0based < g.end_1based
                })
            {
                continue;
            }
            let hits = evidence::project(r, g, accession);
            let (min, max) = evidence::score_range(r, accession);
            let coverage = r
                .motif_coverage
                .iter()
                .find(|c| c.motif_id == accession)
                .expect("validated coverage");
            let provider = r.provider.as_ref().expect("validated provider");
            let summary = TextBlock::new(
                &format!(
                    "jaspar-mapping retained hits | {} | {}\n{}; source floor {}; query complete {}; truncated {}; queried span {}.\nUp = local +; down = local - (relative to displayed DNA). Base = TFBS span; height = raw score {}..{} {}; 2 px minimum marker. Scale fixed per source/report/matrix across TSSs, separate from local PWM.\n{} overlapping hit rows; {} printed labels; every hit in SVG hover and imported-motif-hits.tsv. Missing hits are not evidence of absence.{}",
                    accession,
                    r.report_id,
                    coverage.status.as_str(),
                    coverage
                        .source_minimum_score
                        .map(number)
                        .unwrap_or("not recorded".into()),
                    r.query_complete,
                    r.truncated,
                    if evidence::covers_window(r, g) {
                        "covers this window"
                    } else {
                        "PARTIAL window coverage"
                    },
                    number(min),
                    number(max),
                    provider.score_mode,
                    hits.len(),
                    hits.len().min(LABEL_LIMIT),
                    format!(
                        "{}{}",
                        r.regulatory_subset
                            .as_ref()
                            .map(|s| format!("\n{}", s.summary()))
                            .unwrap_or_default(),
                        if r.warnings.is_empty() {
                            String::new()
                        } else {
                            format!("\nSource warnings: {}", r.warnings.join("; "))
                        }
                    )
                ),
                TEXT_WIDTH - 16.0,
                12.0,
            );
            let labels = hits
                .iter()
                .take(LABEL_LIMIT)
                .enumerate()
                .map(|(i, h)| TextBlock::new(&hit_label(h, i), PLOT_WIDTH, 12.0))
                .collect::<Vec<_>>();
            let height = summary.height()
                + 2.0 * HALF_HEIGHT
                + 40.0
                + labels.iter().map(|l| l.height() + 3.0).sum::<f64>()
                + 20.0;
            lanes.push(Lane {
                source,
                hits,
                summary,
                labels,
                min,
                max,
                height,
            });
        }
        Self {
            missing: lanes.is_empty() && !report.imported_motif_evidence.is_empty(),
            lanes,
        }
    }

    pub(super) fn height(&self) -> f64 {
        if self.missing {
            44.0
        } else {
            self.lanes.iter().map(|l| l.height).sum()
        }
    }

    pub(super) fn draw(&self, parent: &mut Group, axis: LocalAxis, mut y: f64) {
        if self.missing {
            parent.append(text_node(
                MARGIN,
                y + 20.0,
                "Imported hits: no query for this exact matrix/window. Not evidence of absence.",
                12.0,
            ));
        }
        for lane in &self.lanes {
            let mut group = Group::new()
                .set("data-role", "imported-motif-lane")
                .set("data-source-report", lane.source.report.report_id.as_str())
                .set("data-source-sha256", lane.source.source.sha256.as_str())
                .set("data-score-min", lane.min)
                .set("data-score-max", lane.max)
                .set(
                    "data-height-policy",
                    "source_report_matrix_raw_range_2px_floor",
                )
                .set("data-hit-count", lane.hits.len());
            lane.summary
                .draw(&mut group, MARGIN, y + 8.0, "imported-motif-policy");
            let baseline = y + lane.summary.height() + 18.0 + HALF_HEIGHT;
            let clip_id = format!("imported-motif-clip-{baseline}");
            group.append(
                ClipPath::new().set("id", clip_id.as_str()).add(
                    Rectangle::new()
                        .set("x", PLOT_LEFT)
                        .set("y", baseline - HALF_HEIGHT - 2.0)
                        .set("width", PLOT_WIDTH)
                        .set("height", 2.0 * HALF_HEIGHT + 4.0),
                ),
            );
            group.append(
                Line::new()
                    .set("x1", PLOT_LEFT)
                    .set("x2", PLOT_RIGHT)
                    .set("y1", baseline)
                    .set("y2", baseline)
                    .set("stroke", "#a6b2bb"),
            );
            group.append(text_node(MARGIN, baseline - 12.0, "local + (up)", 12.0));
            group.append(text_node(MARGIN, baseline + 20.0, "local - (down)", 12.0));
            for sign in [-1.0, 1.0] {
                for fraction in [0.5, 1.0] {
                    let tick_y = baseline + sign * HALF_HEIGHT * fraction;
                    let value = lane.min * (1.0 - fraction) + lane.max * fraction;
                    group.append(
                        Line::new()
                            .set("x1", PLOT_RIGHT)
                            .set("x2", PLOT_RIGHT + 4.0)
                            .set("y1", tick_y)
                            .set("y2", tick_y)
                            .set("stroke", "#8296a0"),
                    );
                    group.append(text_node(
                        PLOT_RIGHT + 8.0,
                        tick_y + 4.0,
                        &number(value),
                        11.0,
                    ));
                }
            }
            group.append(text_node(
                PLOT_RIGHT + 8.0,
                baseline + 4.0,
                &number(lane.min),
                11.0,
            ));
            for (i, hit) in lane.hits.iter().enumerate() {
                let (x, width) = if axis.length == 1 {
                    (
                        PLOT_LEFT + hit.full_start as f64 * PLOT_WIDTH,
                        (hit.full_end - hit.full_start) as f64 * PLOT_WIDTH,
                    )
                } else {
                    let x = axis.x(hit.full_start as f64 - 0.5);
                    (x, axis.x(hit.full_end as f64 - 0.5) - x)
                };
                let magnitude = lane.min.abs().max(lane.max.abs());
                let fraction = ((hit.hit.score / magnitude - lane.min / magnitude)
                    / (lane.max / magnitude - lane.min / magnitude))
                    .clamp(0.0, 1.0);
                let height = (HALF_HEIGHT * fraction).max(2.0);
                let up = hit.local_strand == TssStrand::Plus;
                let peak_y = baseline + if up { -height } else { height };
                let color = if up { "#9b4c13" } else { "#196a91" };
                group.append(Polygon::new()
                    .set("points", format!("{x},{baseline} {},{peak_y} {},{baseline}", x + width * 0.5, x + width))
                    .set("fill", color).set("fill-opacity", 0.22).set("stroke", color)
                    .set("clip-path", format!("url(#{clip_id})"))
                    .set("stroke-width", 0.8)
                    .set("stroke-dasharray", if hit.clipped { "2,2" } else { "none" })
                    .set("data-role", "imported-motif-hit")
                    .set("data-local-start", hit.start).set("data-local-end", hit.end)
                    .set("data-genomic-start-1based", hit.hit.start_0based + 1)
                    .set("data-genomic-end-1based", hit.hit.end_0based_exclusive)
                    .set("data-local-strand", hit.local_strand.as_str())
                    .set("data-genomic-strand", hit.hit.strand.as_str())
                    .set("data-score", hit.hit.score).set("data-clipped", hit.clipped.to_string())
                    .add(Title::new(format!("{}; score mode {}; query interval {}; pwm-relative score {:?}; original score {}",
                        hit_label(hit, i), hit.hit.score_mode, hit.hit.interval_id, hit.hit.pwm_relative_score, hit.hit.score))));
                if i < LABEL_LIMIT {
                    group.append(
                        text_node(
                            (x + width * 0.5).clamp(PLOT_LEFT, PLOT_RIGHT),
                            peak_y + if up { -3.0 } else { 12.0 },
                            &(i + 1).to_string(),
                            10.0,
                        )
                        .set("text-anchor", "middle"),
                    );
                }
            }
            let mut label_y = baseline + HALF_HEIGHT + 22.0;
            for label in &lane.labels {
                label.draw(&mut group, PLOT_LEFT, label_y, "imported-motif-label");
                label_y += label.height() + 3.0;
            }
            parent.append(group);
            y += lane.height;
        }
    }
}
