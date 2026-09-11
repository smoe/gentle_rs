//! Presentation of engine-projected TSS context; no joins or scientific rescoring.

use super::*;
use gentle_protocol::{isoform_evidence::*, tata_boxes::TataBoxEvidenceKind, tss_profiles::*};

pub(super) fn validate(window: &TssProfileWindow) -> Result<(), String> {
    let Some(context) = &window.detail_context else {
        return Ok(());
    };
    let g = &window.record.geometry;
    let length = g.length().ok_or("Invalid context geometry")?;
    let hash_ok = |hash: &str| {
        hash.len() == 64
            && hash
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    };
    if context.schema != CONTEXT_SCHEMA
        || context.promoter_id != window.record.promoter_id
        || context.geometry != *g
        || context.window_sequence_sha256 != window.record.sequence_sha256
        || !hash_ok(&context.locus_sequence_sha256)
        || !hash_ok(&context.locus_report_sha256)
        || context.non_claims != CONTEXT_NON_CLAIMS
        || context.transcripts.len() > 128
        || context.occupancy.len() > 128
        || context.bindings.len() > 16
        || context.bindings.iter().any(|b| !hash_ok(&b.sha256))
        || !context
            .bindings
            .iter()
            .any(|b| b.role == "tss_detail_locus_report" && b.sha256 == context.locus_report_sha256)
    {
        return Err("Invalid or mismatched TSS detail context binding/size".into());
    }
    let check_span = |span: &TssContextSpan| -> Result<(), String> {
        let start = span.genomic_start_1based.max(g.start_1based);
        let end = span.genomic_end_1based.min(g.end_1based);
        if span.genomic_start_1based == 0
            || end < start
            || span.genomic_end_1based < span.genomic_start_1based
        {
            return Err("Invalid genomic context span".into());
        }
        let expected = match g.strand {
            TssStrand::Plus => (start - g.start_1based, end - g.start_1based + 1),
            TssStrand::Minus => (g.end_1based - end, g.end_1based - start + 1),
        };
        if span.start_0based as u64 != expected.0
            || span.end_0based_exclusive as u64 != expected.1
            || span.end_0based_exclusive > length
            || span.clipped
                != (start != span.genomic_start_1based || end != span.genomic_end_1based)
        {
            return Err("Context span is not an exact clipped TSS projection".into());
        }
        Ok(())
    };
    let mut seen = BTreeSet::new();
    let mut total = 0usize;
    for transcript in &context.transcripts {
        if !window
            .record
            .transcripts
            .contains(&transcript.transcript_id)
            || !seen.insert(&transcript.transcript_id)
            || transcript.genomic_strand != g.strand
        {
            return Err("Unexpected, duplicated or opposite-strand context transcript".into());
        }
        for exon in &transcript.exons {
            if exon.number_5prime_to_3prime == 0 {
                return Err("Invalid exon ordinal".into());
            }
            check_span(&exon.span)?;
        }
        for cds in &transcript.cds {
            check_span(cds)?;
        }
        for codon in &transcript.codons {
            if g.genomic_at(codon.position_0based) != Some(codon.genomic_position_1based)
                || codon.basis.trim().is_empty()
            {
                return Err("Invalid annotation-backed translation boundary".into());
            }
        }
        total += transcript.exons.len() + transcript.cds.len() + transcript.codons.len();
    }
    let mut seen = BTreeSet::new();
    for lane in &context.occupancy {
        if !seen.insert((&lane.group_id, &lane.lane_id))
            || !lane.display_abs_max_score.is_finite()
            || lane.display_abs_max_score < 0.0
            || (lane.state != GeneLocusOccupancyLaneState::Available && !lane.intervals.is_empty())
        {
            return Err("Invalid occupancy identity, availability or scale".into());
        }
        for interval in &lane.intervals {
            check_span(&interval.span)?;
            if interval.score.is_some_and(|score| {
                !score.is_finite() || (score != 0.0 && lane.display_abs_max_score == 0.0)
            }) {
                return Err("Occupancy score outside its inherited locus scale".into());
            }
        }
        total += lane.intervals.len();
    }
    if let Some(tata) = &context.tata {
        if !hash_ok(&tata.report_sha256)
            || tata.rows.len() > 256
            || !context
                .bindings
                .iter()
                .any(|b| b.role == "tss_detail_tata_report" && b.sha256 == tata.report_sha256)
        {
            return Err("Invalid TATA evidence binding/size".into());
        }
        for row in &tata.rows {
            check_span(&row.span)?;
            if row.evidence.llr_bits.is_some_and(|s| !s.is_finite()) {
                return Err("Non-finite TATA score".into());
            }
        }
        total += tata.rows.len();
    }
    if total > 100_000 {
        return Err("TSS context exceeds 100,000 projected intervals".into());
    }
    Ok(())
}

enum Content<'a> {
    Transcript(&'a TssContextTranscript),
    Occupancy(&'a TssContextOccupancyLane),
    Tata(&'a TssContextTataRow),
    Unavailable(&'a str),
}

struct ContextRow<'a> {
    label: TextBlock,
    details: TextBlock,
    content: Content<'a>,
}

impl<'a> ContextRow<'a> {
    fn new(label: &str, details: &str, content: Content<'a>) -> Self {
        Self {
            label: TextBlock::new(label, LABEL_WIDTH, 12.0),
            details: TextBlock::new(details, PAGE_WIDTH - PLOT_RIGHT - MARGIN - 14.0, 11.0),
            content,
        }
    }
    fn height(&self) -> f64 {
        58.0_f64.max(self.label.height()).max(self.details.height()) + 12.0
    }
    fn draw(&self, group: &mut Group, top: f64, axis: LocalAxis) {
        self.label.draw(group, MARGIN, top, "context-row-label");
        self.details
            .draw(group, PLOT_RIGHT + 14.0, top, "context-row-details");
        let center = top + 26.0;
        let rect = |span: &TssContextSpan, y: f64, height: f64| {
            let (x, width) = axis.interval(span.start_0based, span.end_0based_exclusive);
            Rectangle::new()
                .set("x", x)
                .set("y", y)
                .set("width", width.max(1.0))
                .set("height", height)
                .set("data-start-0based", span.start_0based)
                .set("data-end-0based-exclusive", span.end_0based_exclusive)
                .set("data-genomic-start", span.genomic_start_1based)
                .set("data-genomic-end", span.genomic_end_1based)
        };
        match &self.content {
            Content::Transcript(t) => {
                let mut row = Group::new()
                    .set("data-role", "tss-context-transcript")
                    .set("data-transcript-id", t.transcript_id.as_str());
                if let (Some(first), Some(last)) = (t.exons.first(), t.exons.last()) {
                    row.append(
                        Line::new()
                            .set("x1", axis.x(first.span.start_0based as f64))
                            .set(
                                "x2",
                                axis.x(last.span.end_0based_exclusive.saturating_sub(1) as f64),
                            )
                            .set("y1", center)
                            .set("y2", center)
                            .set("stroke", "#65727a"),
                    );
                }
                for exon in &t.exons {
                    row.append(
                        rect(&exon.span, center - 5.0, 10.0)
                            .set("fill", "#bfdce5")
                            .set("stroke", "#39758c")
                            .set("data-role", "context-exon")
                            .set("data-exon-number", exon.number_5prime_to_3prime)
                            .add(Title::new(format!(
                                "Exon {}: genomic {}..{}",
                                exon.number_5prime_to_3prime,
                                exon.span.genomic_start_1based,
                                exon.span.genomic_end_1based
                            ))),
                    );
                    let (x, width) =
                        axis.interval(exon.span.start_0based, exon.span.end_0based_exclusive);
                    if width > 38.0 {
                        row.append(text_node(
                            x + 3.0,
                            center - 10.0,
                            &format!("E{}", exon.number_5prime_to_3prime),
                            10.0,
                        ));
                    }
                }
                for cds in &t.cds {
                    row.append(
                        rect(cds, center - 8.0, 16.0)
                            .set("fill", "#39758c")
                            .set("data-role", "context-cds"),
                    );
                }
                for codon in &t.codons {
                    let (kind, color) = match codon.kind {
                        GeneLocusCodonKind::Start => ("start", "#237c47"),
                        GeneLocusCodonKind::Stop => ("stop", "#bf3847"),
                    };
                    row.append(
                        Line::new()
                            .set("x1", axis.x(codon.position_0based as f64))
                            .set("x2", axis.x(codon.position_0based as f64))
                            .set("y1", center - 15.0)
                            .set("y2", center + 17.0)
                            .set("stroke", color)
                            .set("stroke-width", 2)
                            .set("data-role", format!("context-translation-{kind}"))
                            .set("data-genomic-position", codon.genomic_position_1based)
                            .add(Title::new(format!(
                                "Translation {kind}: {} ({})",
                                codon.genomic_position_1based, codon.basis
                            ))),
                    );
                }
                if t.exons.is_empty() {
                    row.append(text_node(
                        PLOT_LEFT + 8.0,
                        center,
                        "No exon geometry in this window",
                        12.0,
                    ));
                }
                group.append(row);
            }
            Content::Occupancy(lane) => {
                let mut row = Group::new()
                    .set("data-role", "tss-context-occupancy")
                    .set("data-source-id", lane.source_id.as_str())
                    .set("data-group-id", lane.group_id.as_str())
                    .set("data-lane-id", lane.lane_id.as_str())
                    .set("data-state", lane.state.as_str());
                if lane.intervals.is_empty() {
                    let text = if lane.state == GeneLocusOccupancyLaneState::Available {
                        "No supplied intervals here; not a measured zero"
                    } else {
                        lane.state.as_str()
                    };
                    row.append(text_node(PLOT_LEFT + 8.0, center, text, 12.0));
                } else {
                    row.append(
                        Line::new()
                            .set("x1", PLOT_LEFT)
                            .set("x2", PLOT_RIGHT)
                            .set("y1", center)
                            .set("y2", center)
                            .set("stroke", "#d5dcdf"),
                    );
                    for interval in &lane.intervals {
                        let (y, height, fill) = match interval.score {
                            Some(score) => {
                                let height = if lane.display_abs_max_score > 0.0 {
                                    (score.abs() / lane.display_abs_max_score).min(1.0) * 22.0
                                } else {
                                    0.0
                                };
                                (
                                    if score >= 0.0 {
                                        center - height
                                    } else {
                                        center
                                    },
                                    height.max(0.7),
                                    if lane.role == GeneLocusOccupancyLaneRole::Experimental {
                                        "#9e405c"
                                    } else {
                                        "#6b7580"
                                    },
                                )
                            }
                            None => (center - 4.0, 8.0, "none"),
                        };
                        row.append(
                            rect(&interval.span, y, height)
                                .set("fill", fill)
                                .set("stroke", if fill == "none" { "#677580" } else { fill })
                                .set("data-role", "context-occupancy-interval")
                                .set("data-interval-id", interval.interval_id.as_str())
                                .set(
                                    "data-clipped-to-locus-scale",
                                    interval
                                        .score
                                        .is_some_and(|s| s.abs() > lane.display_abs_max_score)
                                        .to_string(),
                                )
                                .add(Title::new(format!(
                                    "{}; raw score {:?}; inherited |max| {}; {}",
                                    interval.interval_id,
                                    interval.score,
                                    lane.display_abs_max_score,
                                    interval.label.as_deref().unwrap_or("")
                                ))),
                        );
                    }
                }
                group.append(row);
            }
            Content::Tata(row) => {
                let (role, color) = match row.evidence.evidence_kind {
                    TataBoxEvidenceKind::SourceAnnotation => ("tata-annotation", "#79579a"),
                    TataBoxEvidenceKind::MotifPrediction => ("tata-prediction", "#b58325"),
                    TataBoxEvidenceKind::EpdClassification => ("tata-epd-tss", "#79579a"),
                };
                group.append(
                    rect(&row.span, center - 7.0, 14.0)
                        .set(
                            "fill",
                            if role == "tata-prediction" {
                                "none"
                            } else {
                                color
                            },
                        )
                        .set("stroke", color)
                        .set("stroke-width", 2)
                        .set("data-role", role)
                        .set("data-evidence-id", row.evidence.row_id.as_str())
                        .set("data-genomic-strand", row.genomic_strand.as_str())
                        .add(Title::new(format!(
                            "{}; {}; genomic strand {}",
                            row.evidence.label,
                            row.evidence.geometry_kind,
                            row.genomic_strand.as_str()
                        ))),
                );
            }
            Content::Unavailable(message) => group.append(
                text_node(PLOT_LEFT + 8.0, center, message, 12.0)
                    .set("data-role", "context-unavailable"),
            ),
        }
    }
}

pub(super) struct ContextLayout<'a> {
    context: &'a TssDetailContext,
    intro: TextBlock,
    rows: Vec<ContextRow<'a>>,
    footer: TextBlock,
}

impl<'a> ContextLayout<'a> {
    pub(super) fn new(c: &'a TssDetailContext) -> Self {
        let mut rows = Vec::new();
        for t in &c.transcripts {
            rows.push(ContextRow::new(
                &format!("{}\n{}", t.transcript_id, t.label),
                &format!(
                    "Exons / CDS\nGenomic strand {}\nTranslation: green start, red stop",
                    t.genomic_strand.as_str()
                ),
                Content::Transcript(t),
            ));
        }
        if c.transcripts.is_empty() {
            rows.push(ContextRow::new(
                "Gene structure",
                "",
                Content::Unavailable("Transcript geometry unavailable"),
            ));
        }
        for lane in &c.occupancy {
            let role = serde_json::to_value(lane.role)
                .ok()
                .and_then(|v| v.as_str().map(str::to_string))
                .unwrap_or_default();
            let metadata = [
                Some(lane.source_kind.as_str()),
                lane.assay.as_deref(),
                lane.factor.as_deref(),
                lane.condition.as_deref(),
                lane.cell_line.as_deref(),
            ]
            .into_iter()
            .flatten()
            .collect::<Vec<_>>()
            .join(" | ");
            rows.push(ContextRow::new(
                &format!("{}\n{}", lane.group_label, lane.label),
                &format!(
                    "{metadata}\n{role}; locus |max| {}\nSource: {}",
                    lane.display_abs_max_score, lane.source_id
                ),
                Content::Occupancy(lane),
            ));
        }
        if c.occupancy.is_empty() {
            rows.push(ContextRow::new(
                "CUT&RUN / chromatin",
                "",
                Content::Unavailable("No occupancy lanes supplied"),
            ));
        }
        if let Some(tata) = &c.tata {
            for row in &tata.rows {
                let kind = match row.evidence.evidence_kind {
                    TataBoxEvidenceKind::SourceAnnotation => "TATA annotation",
                    TataBoxEvidenceKind::MotifPrediction => "TATA prediction",
                    TataBoxEvidenceKind::EpdClassification => "EPD TSS classification",
                };
                let detail = if row.evidence.evidence_kind == TataBoxEvidenceKind::EpdClassification
                {
                    format!(
                        "TATA-positive: {}\nTSS marker, not box coordinates",
                        match row.evidence.tata_positive {
                            Some(true) => "yes",
                            Some(false) => "no",
                            None => "unavailable",
                        }
                    )
                } else {
                    format!(
                        "Genomic strand {}\nLLR bits: {}",
                        row.genomic_strand.as_str(),
                        row.evidence
                            .llr_bits
                            .map(|s| format!("{s:.2}"))
                            .unwrap_or_else(|| "not scored".into())
                    )
                };
                rows.push(ContextRow::new(
                    &format!("{kind}\n{}", row.evidence.label),
                    &detail,
                    Content::Tata(row),
                ));
            }
            if tata.rows.is_empty() {
                rows.push(ContextRow::new(
                    "TATA evidence",
                    "",
                    Content::Unavailable(
                        "No reported TATA evidence in this window; not proof of absence",
                    ),
                ));
            }
        } else {
            rows.push(ContextRow::new(
                "TATA evidence",
                "",
                Content::Unavailable("TATA screen not supplied"),
            ));
        }
        let mut footer = format!(
            "Context: {} | locus report SHA-256 {}\n{}",
            c.locus_seq_id, c.locus_report_sha256, c.non_claims
        );
        if let Some(tata) = &c.tata {
            let _ = write!(
                footer,
                "\nTATA {} | SHA-256 {} | motif {} | {} | EPD {}",
                tata.report_id,
                tata.report_sha256,
                tata.motif_id.as_deref().unwrap_or("not used"),
                tata.score_policy,
                tata.epd_status
            );
            for warning in &tata.warnings {
                let _ = write!(footer, "\nTATA warning: {warning}");
            }
        }
        for warning in &c.warnings {
            let _ = write!(footer, "\nContext warning: {warning}");
        }
        Self {
            context: c,
            intro: TextBlock::new(
                "LOCUS CONTEXT ON THE SAME TSS AXIS\nThin boxes: exons; thick boxes: CDS. Translation markers are annotation-backed. Signal heights retain the source-locus scale; gaps are not interpolated. TATA predictions remain separate from annotations and EPD TSS classifications.",
                TEXT_WIDTH,
                13.0,
            ),
            rows,
            footer: TextBlock::new(&footer, TEXT_WIDTH, 11.0),
        }
    }
    pub(super) fn height(&self) -> f64 {
        self.intro.height()
            + 14.0
            + self.rows.iter().map(ContextRow::height).sum::<f64>()
            + self.footer.height()
            + 24.0
    }
    pub(super) fn draw(&self, top: f64, axis: LocalAxis, geometry: &TssGeometry) -> Group {
        let mut group = Group::new()
            .set("data-role", "tss-detail-context")
            .set(
                "data-locus-report-sha256",
                self.context.locus_report_sha256.as_str(),
            )
            .set("data-y", top)
            .set("data-height", self.height());
        self.intro.draw(&mut group, MARGIN, top, "context-legend");
        let mut y = top + self.intro.height() + 14.0;
        for row in &self.rows {
            group.append(
                Line::new()
                    .set("x1", axis.x(geometry.upstream_bp as f64))
                    .set("x2", axis.x(geometry.upstream_bp as f64))
                    .set("y1", y)
                    .set("y2", y + row.height())
                    .set("stroke", "#c7cdd1")
                    .set("stroke-dasharray", "2,3")
                    .set("data-role", "context-tss-guide"),
            );
            row.draw(&mut group, y, axis);
            y += row.height();
        }
        self.footer
            .draw(&mut group, MARGIN, y, "context-provenance");
        group
    }
}
