//! Engine-owned per-unit pileup and explicit views over saved projections, never source alignments.

use super::*;
use coverage::{Compatibility, CoverageReport, InputBindings, SampleInput};
use evidence::{CompositionEvidence, TailEvidence};
use gentle_render::read_end_pileup::{Band, Row, Scene};

pub const SCHEMA: &str = "gentle.rna_read_end_pileup.v1";
pub const VIEW_SCHEMA: &str = "gentle.rna_read_end_pileup_view.v1";

pub struct Input<'a> {
    pub sample: SampleInput<'a>,
    pub model: &'a ModelDocument,
    pub annotation: &'a crate::engine::PrimerPayloadAnnotation,
    pub composition: &'a CompositionEvidence,
    pub tails: &'a [TailEvidence],
    pub tail_file_sha256: Option<&'a str>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ConnectorKind {
    SpliceN,
    DeletionOrGap,
    SegmentLink,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Connector {
    pub interval: Interval,
    pub kind: ConnectorKind,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Observation {
    pub observation_id: String,
    pub row_ordinal: u64,
    pub blocks: Vec<SupportedBlock>,
    pub connectors: Vec<Connector>,
    pub query_to_reference_strand: char,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Hypothesis {
    pub projection: HypothesisProjection,
    pub observations: Vec<Observation>,
    pub segment_links: Vec<Connector>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Lane {
    pub sample_id: String,
    pub bindings: InputBindings,
    pub unit: UnitProjection,
    pub rows: Vec<ObservedSamRowAssessment>,
    pub compatibility: Option<Compatibility>,
    pub hypotheses: Vec<Hypothesis>,
    pub composition: CompositionEvidence,
    pub tails: Vec<TailEvidence>,
    pub tail_file_sha256: Option<String>,
    pub composition_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub schema: String,
    pub identity_scope: String,
    pub source_id: String,
    pub gene_id: String,
    pub transcript_id: String,
    pub model: TranscriptModel,
    pub boundary: ModelBoundary,
    /// Stable bounds include every gene-assigned saved block, not just the visible page.
    pub display_interval: Interval,
    pub coverage: CoverageReport,
    pub lanes: Vec<Lane>,
    pub nonclaims: Vec<String>,
}

fn bands(observation: &ObservedSamObservation) -> Vec<Connector> {
    let mut out: Vec<_> = observation
        .skips
        .iter()
        .map(|&interval| Connector {
            interval,
            kind: ConnectorKind::SpliceN,
        })
        .collect();
    let mut blocks: Vec<_> = observation.blocks.iter().map(|b| b.target).collect();
    blocks.sort_unstable();
    for pair in blocks.windows(2) {
        let gap = (pair[0].1, pair[1].0);
        if gap.0 >= gap.1 {
            continue;
        }
        // Remove actual recorded N intervals; any remainder is not an invented splice.
        let mut position = gap.0;
        let mut skips = observation.skips.clone();
        skips.sort_unstable();
        for (a, z) in skips {
            if z <= position || a >= gap.1 {
                continue;
            }
            if position < a {
                out.push(Connector {
                    interval: (position, a.min(gap.1)),
                    kind: ConnectorKind::DeletionOrGap,
                });
            }
            position = position.max(z).min(gap.1);
        }
        if position < gap.1 {
            out.push(Connector {
                interval: (position, gap.1),
                kind: ConnectorKind::DeletionOrGap,
            });
        }
    }
    out
}

/// All references are already authenticated by the service. This function does no I/O.
pub fn build(
    source: &str,
    gene: &str,
    transcript: &str,
    probabilities: &[f64],
    inputs: &[Input<'_>],
    limits: Limits,
) -> Result<Report, EngineError> {
    build_with_control(
        source,
        gene,
        transcript,
        probabilities,
        inputs,
        limits,
        &mut || Ok(()),
    )
}

pub fn build_with_control(
    source: &str,
    gene: &str,
    transcript: &str,
    probabilities: &[f64],
    inputs: &[Input<'_>],
    limits: Limits,
    control: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Report, EngineError> {
    control()?;
    if !token(transcript) || inputs.is_empty() {
        return Err(fail("pileup_scope"));
    }
    let coverage_inputs: Vec<_> = inputs
        .iter()
        .map(|i| SampleInput {
            sample_id: i.sample.sample_id,
            gene_id: i.sample.gene_id,
            bindings: i.sample.bindings,
            analysis: i.sample.analysis,
            projection: i.sample.projection,
        })
        .collect();
    let coverage = coverage::summarize_with_control(
        source,
        gene,
        Some(transcript),
        probabilities,
        &coverage_inputs,
        limits,
        control,
    )?;
    let mut meter = Meter {
        limits,
        work: 0,
        control: Some(control),
    };
    let (mut chosen_model, mut boundary, mut lanes) = (None, None, vec![]);
    let mut total_rows = 0usize;
    for input in inputs {
        let models: Vec<_> = input
            .model
            .transcripts
            .iter()
            .filter(|m| m.transcript_id == transcript)
            .collect();
        if models.len() != 1 || input.model.source_id != source || models[0].gene_id != gene {
            return Err(fail("pileup_selected_model"));
        }
        let model = models[0];
        let prepared = prepare(model, &mut meter)?;
        let annotation = projection::annotation_index(input.model, input.annotation, &mut meter)?;
        projection::checked_annotation(transcript, Some(&prepared), &annotation, &mut meter)?;
        let model_json = serde_json::to_value(model).map_err(|e| fail(&e.to_string()))?;
        if chosen_model
            .as_ref()
            .is_some_and(|m| serde_json::to_value(m).ok() != Some(model_json.clone()))
        {
            return Err(fail("pileup_model_scope"));
        }
        chosen_model = Some(model.clone());
        boundary = Some(prepared.boundary.clone());
        let upstream: BTreeMap<_, _> = input
            .sample
            .analysis
            .units
            .iter()
            .map(|u| (&u.unit, u))
            .collect();
        let observations: BTreeMap<_, _> = input
            .sample
            .analysis
            .observations
            .iter()
            .map(|o| (o.observation_id.as_str(), o))
            .collect();
        if observations.len() != input.sample.analysis.observations.len() {
            return Err(fail("pileup_observation_binding"));
        }
        for unit in &input.sample.projection.units {
            meter.charge(1)?;
            let saved = upstream
                .get(&unit.unit)
                .ok_or_else(|| fail("pileup_unit_binding"))?;
            if coverage::gene_assignment(unit, !saved.unavailable.is_empty(), gene).is_err() {
                continue;
            }
            if lanes.len() >= limits.max_units {
                return Err(fail("pileup_lane_budget"));
            }
            let mut hypotheses = vec![];
            for hypothesis in &unit.hypotheses {
                let mut rows = vec![];
                for id in &hypothesis.observation_ids {
                    meter.charge(1)?;
                    total_rows = total_rows
                        .checked_add(1)
                        .filter(|n| *n <= limits.max_observations)
                        .ok_or_else(|| fail("pileup_observation_budget"))?;
                    let o = observations
                        .get(id.as_str())
                        .ok_or_else(|| fail("pileup_observation_binding"))?;
                    if o.unit != unit.unit || !o.on_optimal_path || o.contig != model.contig {
                        return Err(fail("pileup_observation_binding"));
                    }
                    meter.charge(o.blocks.len() + o.skips.len())?;
                    rows.push(Observation {
                        observation_id: id.clone(),
                        row_ordinal: o.source_ordinal,
                        blocks: o.blocks.clone(),
                        connectors: bands(o),
                        query_to_reference_strand: o.query_to_reference_strand,
                    });
                }
                // Query order is provenance, not a new splice or transcript-orientation assertion.
                rows.sort_by_key(|o| {
                    (
                        o.blocks
                            .iter()
                            .map(|b| b.query.0)
                            .min()
                            .unwrap_or(usize::MAX),
                        o.observation_id.clone(),
                    )
                });
                let links = rows
                    .windows(2)
                    .filter_map(|p| {
                        let a = if p[0].query_to_reference_strand == '+' {
                            p[0].blocks.iter().map(|b| b.target.1).max()?
                        } else {
                            p[0].blocks.iter().map(|b| b.target.0).min()?
                        };
                        let z = if p[1].query_to_reference_strand == '+' {
                            p[1].blocks.iter().map(|b| b.target.0).min()?
                        } else {
                            p[1].blocks.iter().map(|b| b.target.1).max()?
                        };
                        (a != z).then_some(Connector {
                            interval: (a.min(z), a.max(z)),
                            kind: ConnectorKind::SegmentLink,
                        })
                    })
                    .collect();
                hypotheses.push(Hypothesis {
                    projection: hypothesis.clone(),
                    observations: rows,
                    segment_links: links,
                });
            }
            lanes.push(Lane {
                sample_id: input.sample.sample_id.into(),
                bindings: input.sample.bindings.clone(),
                unit: unit.clone(),
                rows: saved.rows.clone(),
                compatibility: coverage::compatibility(unit, transcript),
                hypotheses,
                composition: input.composition.clone(),
                composition_sha256: hash(input.composition)?,
                tail_file_sha256: input.tail_file_sha256.map(str::to_owned),
                tails: input
                    .tails
                    .iter()
                    .filter(|t| t.row.unit == unit.unit)
                    .cloned()
                    .collect(),
            });
        }
    }
    let model = chosen_model.ok_or_else(|| fail("pileup_model"))?;
    let mut display_interval = (
        model.genomic_exons_0based_half_open[0].0,
        model
            .genomic_exons_0based_half_open
            .last()
            .ok_or_else(|| fail("pileup_model"))?
            .1,
    );
    for block in lanes
        .iter()
        .flat_map(|l| &l.hypotheses)
        .flat_map(|h| &h.observations)
        .flat_map(|o| &o.blocks)
    {
        display_interval.0 = display_interval.0.min(block.target.0);
        display_interval.1 = display_interval.1.max(block.target.1);
    }
    Ok(Report {schema: SCHEMA.into(), identity_scope: "sam_declared_query_unit".into(), source_id: source.into(), gene_id: gene.into(),
        transcript_id: transcript.into(), model, display_interval, boundary: boundary.ok_or_else(|| fail("pileup_model"))?, coverage, lanes,
        nonclaims: vec!["Incomplete aligner/region views are not an exhaustive locus search.".into(),
            "Conditional coverage is compatibility support, not resolved isoforms or molecule counts.".into(),
            "Clip descriptors do not reinterpret censored ends; composition and signal estimates are separate evidence.".into()]})
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SortKey {
    Sample,
    End5,
    End3,
    Status5,
    Status3,
    Value5,
    Value3,
    Compatibility,
    ReadLength,
    Clip5,
    Clip3,
    UnitKey,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Status {
    Observed,
    RightCensored,
    Missing,
}

fn status(metric: &Metric) -> Status {
    match metric {
        Metric::Observed { .. } => Status::Observed,
        Metric::RightCensored { .. } => Status::RightCensored,
        Metric::Missing { .. } => Status::Missing,
    }
}
fn value(metric: &Metric) -> Option<usize> {
    match metric {
        Metric::Observed { value_bp } => Some(*value_bp),
        Metric::RightCensored { lower_bound_bp } => Some(*lower_bound_bp),
        _ => None,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Filter {
    Status5 { value: Status },
    Status3 { value: Status },
    Compatibility { value: Compatibility },
    Clip5 { value: ClipKind },
    Clip3 { value: ClipKind },
    Composition { available: bool },
    Tail { available: bool },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct View {
    pub schema: String,
    pub sort: SortKey,
    pub descending: bool,
    pub filters: Vec<Filter>,
    pub offset: usize,
    pub count: usize,
    pub max_sublanes: usize,
    pub width: u32,
    /// Genomic scale only; no implicit or unlabelled intron compression.
    pub scale: String,
    pub channels: Channels,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Channels {
    pub utr5: bool,
    pub utr3: bool,
    pub read_length: bool,
    pub clips: bool,
    pub composition: bool,
    pub tail: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Page {
    pub schema: String,
    pub view: View,
    pub view_sha256: String,
    pub report_sha256: String,
    pub total_gene_assigned: usize,
    pub total_selected_transcript_compatible: usize,
    pub total_matching: usize,
    pub shown_start: usize,
    pub shown_end: usize,
    pub coverage: CoverageReport,
    pub lanes: Vec<Lane>,
}

fn selected<'a>(lane: &'a Lane, id: &str) -> Vec<&'a CandidateProjection> {
    lane.hypotheses
        .iter()
        .filter_map(|h| {
            h.projection
                .candidates
                .iter()
                .find(|c| c.transcript_id == id)
        })
        .collect()
}
fn agreed<T: Eq + Copy>(values: impl Iterator<Item = T>) -> Option<T> {
    let mut values = values;
    let first = values.next()?;
    values.all(|v| v == first).then_some(first)
}
fn clip(lane: &Lane, id: &str, five: bool) -> Option<usize> {
    let candidates = selected(lane, id);
    if candidates.len() != lane.hypotheses.len() {
        return None;
    }
    agreed(candidates.iter().map(|c| {
        if five {
            c.clip5.unaligned_query_bases
        } else {
            c.clip3.unaligned_query_bases
        }
    }))
    .flatten()
}
fn sort_value(lane: &Lane, id: &str, key: SortKey) -> Option<usize> {
    match key {
        SortKey::Value5 => value(&lane.unit.utr5),
        SortKey::Value3 => value(&lane.unit.utr3),
        SortKey::ReadLength => value(&lane.unit.read_length),
        SortKey::Clip5 => clip(lane, id, true),
        SortKey::Clip3 => clip(lane, id, false),
        SortKey::Status5 => Some(status(&lane.unit.utr5) as usize),
        SortKey::Status3 => Some(status(&lane.unit.utr3) as usize),
        SortKey::Compatibility => lane.compatibility.map(|c| c as usize),
        SortKey::End5 | SortKey::End3 => {
            let candidates = selected(lane, id);
            if candidates.len() != lane.hypotheses.len() {
                return None;
            }
            agreed(candidates.iter().map(|c| {
                if key == SortKey::End5 {
                    c.transcript_intervals.first().map(|i| i.0)
                } else {
                    c.transcript_intervals.last().map(|i| i.1)
                }
            }))
            .flatten()
        }
        _ => None,
    }
}
fn matches(lane: &Lane, transcript: &str, filters: &[Filter]) -> bool {
    filters.iter().all(|filter| match filter {
        Filter::Status5 { value } => status(&lane.unit.utr5) == *value,
        Filter::Status3 { value } => status(&lane.unit.utr3) == *value,
        Filter::Compatibility { value } => lane.compatibility == Some(*value),
        Filter::Clip5 { value } | Filter::Clip3 { value } => {
            let candidates = selected(lane, transcript);
            !candidates.is_empty()
                && candidates.len() == lane.hypotheses.len()
                && candidates.iter().all(|c| {
                    if matches!(filter, Filter::Clip5 { .. }) {
                        c.clip5.kind == *value
                    } else {
                        c.clip3.kind == *value
                    }
                })
        }
        Filter::Composition { available } => {
            (lane.composition.availability == "available") == *available
        }
        Filter::Tail { available } => {
            lane.tails.iter().any(|t| t.accepted_estimate_nt.is_some()) == *available
        }
    })
}

pub fn hash<T: Serialize>(record: &T) -> Result<String, EngineError> {
    crate::digest_utils::canonical_json_string(record)
        .map(|s| crate::digest_utils::sha256_prefixed_str(&s))
        .map_err(|e| fail(&e.to_string()))
}

/// Sorting/filtering selects whole units, never one favoured hypothesis. Denominators are copied unchanged.
pub fn page(report: &Report, view: &View) -> Result<Page, EngineError> {
    if report.schema != SCHEMA
        || view.schema != VIEW_SCHEMA
        || view.scale != "genomic"
        || view.count == 0
        || view.count > 20_000
        || view.max_sublanes == 0
        || view.max_sublanes > 20_000
        || !(900..=8000).contains(&view.width)
        || view.filters.len() > 32
    {
        return Err(fail("pileup_view"));
    }
    let mut lanes: Vec<_> = report
        .lanes
        .iter()
        .filter(|l| matches(l, &report.transcript_id, &view.filters))
        .collect();
    lanes.sort_by(|a, b| {
        let primary = match view.sort {
            SortKey::Sample => a.sample_id.cmp(&b.sample_id),
            SortKey::UnitKey => a.unit.unit.cmp(&b.unit.unit),
            key => sort_value(a, &report.transcript_id, key).cmp(&sort_value(
                b,
                &report.transcript_id,
                key,
            )),
        };
        let primary = if view.descending {
            primary.reverse()
        } else {
            primary
        };
        primary.then_with(|| (&a.sample_id, &a.unit.unit).cmp(&(&b.sample_id, &b.unit.unit)))
    });
    if view.offset > lanes.len() {
        return Err(fail("pileup_page_range"));
    }
    let end = view
        .offset
        .checked_add(view.count)
        .ok_or_else(|| fail("pileup_page_range"))?
        .min(lanes.len());
    let rows: usize = lanes[view.offset..end]
        .iter()
        .map(|l| l.hypotheses.len())
        .sum();
    if rows > view.max_sublanes {
        return Err(fail(&format!(
            "pileup_sublane_budget: required={rows}, allowed={}",
            view.max_sublanes
        )));
    }
    Ok(Page {
        schema: "gentle.rna_read_end_pileup_page.v1".into(),
        view: view.clone(),
        view_sha256: hash(view)?,
        report_sha256: hash(report)?,
        total_gene_assigned: report.lanes.len(),
        total_selected_transcript_compatible: report
            .coverage
            .samples
            .iter()
            .try_fold(0usize, |n, s| {
                n.checked_add(s.selected_transcript_compatible_units)
            })
            .ok_or_else(|| fail("pileup_count_overflow"))?,
        total_matching: lanes.len(),
        shown_start: view.offset,
        shown_end: end,
        coverage: report.coverage.clone(),
        lanes: lanes[view.offset..end]
            .iter()
            .map(|l| (*l).clone())
            .collect(),
    })
}

fn band(interval: Interval, kind: &str, label: String) -> Band {
    Band {
        interval: (interval.0 as u64, interval.1 as u64),
        kind: kind.into(),
        label,
    }
}
fn json<T: Serialize>(value: &T) -> Result<String, EngineError> {
    serde_json::to_string(value).map_err(|e| fail(&e.to_string()))
}

pub fn scene(report: &Report, page: &Page) -> Result<Scene, EngineError> {
    if page.report_sha256 != hash(report)? || page.view_sha256 != hash(&page.view)? {
        return Err(fail("pileup_page_binding"));
    }
    let mut model = vec![];
    for (kind, ranges) in [
        ("exon", &report.model.genomic_exons_0based_half_open),
        ("cds", &report.model.genomic_cds_0based_half_open),
        ("stop", &report.model.genomic_stop_codon_0based_half_open),
    ] {
        model.extend(ranges.iter().map(|&i| {
            band(
                i,
                kind,
                format!(
                    "{}; transcript CDS {:?}; UTR3 boundary {:?}",
                    report.transcript_id, report.boundary.cds, report.boundary.utr3_start
                ),
            )
        }));
    }
    let mut rows = vec![];
    for lane in &page.lanes {
        for h in &lane.hypotheses {
            let mut bands = vec![];
            for o in &h.observations {
                bands.extend(
                    o.blocks
                        .iter()
                        .map(|b| band(b.target, "aligned", o.observation_id.clone())),
                );
                bands.extend(o.connectors.iter().map(|c| {
                    band(
                        c.interval,
                        match c.kind {
                            ConnectorKind::SpliceN => "splice_n",
                            ConnectorKind::DeletionOrGap => "deletion_or_gap",
                            ConnectorKind::SegmentLink => "segment_link",
                        },
                        o.observation_id.clone(),
                    )
                }));
            }
            bands.extend(h.segment_links.iter().map(|c| {
                band(
                    c.interval,
                    "segment_link",
                    "Between saved observations; not an intron".into(),
                )
            }));
            let candidate = h
                .projection
                .candidates
                .iter()
                .find(|c| c.transcript_id == report.transcript_id);
            let description = if let Some(c) = candidate {
                let channels = &page.view.channels;
                let mut fields = vec![format!(
                    "sense {:?} (query to annotation, not RNA strand)",
                    c.sense_to_annotation
                )];
                if channels.utr5 {
                    fields.push(format!(
                        "5' {} [unit {}]",
                        json(&c.utr5)?,
                        json(&lane.unit.utr5)?
                    ));
                }
                if channels.utr3 {
                    fields.push(format!(
                        "3' {} [unit {}]",
                        json(&c.utr3)?,
                        json(&lane.unit.utr3)?
                    ));
                }
                if channels.clips {
                    fields.push(format!("clips {} / {}", json(&c.clip5)?, json(&c.clip3)?));
                }
                fields.join("; ")
            } else {
                "Selected transcript incompatible under this hypothesis".into()
            };
            let mut off_axis = vec![description];
            if page.view.channels.read_length {
                off_axis.push(format!("length {}", json(&lane.unit.read_length)?));
            }
            if page.view.channels.composition {
                off_axis.push(format!("A/T {}", lane.composition.reason));
            }
            if page.view.channels.tail {
                off_axis.push(if lane.tails.is_empty() {
                    "tail unavailable: no_imported_tail_evidence".into()
                } else {
                    format!("tail {}", json(&lane.tails)?)
                });
            }
            rows.push(Row {
                label: format!(
                    "{} / {} / H{}",
                    lane.sample_id, lane.unit.unit.unit_id, h.projection.hypothesis_index
                ),
                details: json(lane)?,
                off_axis: off_axis.join("; "),
                bands,
            });
        }
    }
    let mut header = vec![
        format!(
            "{} / {} / {}: genomic scale, declared SAM units",
            report.source_id, report.gene_id, report.transcript_id
        ),
        format!(
            "Lanes {}..{} of {} matching / {} gene-assigned; view {}",
            page.shown_start,
            page.shown_end,
            page.total_matching,
            page.total_gene_assigned,
            page.view_sha256
        ),
    ];
    for sample in &report.coverage.samples {
        header.push(format!(
            "{}: gene-assigned {} / selected-transcript-compatible {} (declared {}, assessed {})",
            sample.sample_id,
            sample.gene_assigned_units,
            sample.selected_transcript_compatible_units,
            sample.declared_units,
            sample.assessed_units
        ));
    }
    Ok(Scene {
        contig: report.model.contig.clone(),
        start: report.display_interval.0 as u64,
        end: report.display_interval.1 as u64,
        strand: if report.model.strand == '+' { 1 } else { -1 },
        width: page.view.width,
        header,
        model,
        rows,
    })
}

/// JSON-encoded cells prevent tabs/newlines in declared names or provenance from changing TSV structure.
pub fn tsv(page: &Page) -> Result<String, EngineError> {
    let mut out = format!(
        "# view={} report={} gene_assigned={} selected_transcript_compatible={} matching={} shown={}..{}\nrecord_type\tsample_json\tunit_json\tevidence_json\n",
        page.view_sha256,
        page.report_sha256,
        page.total_gene_assigned,
        page.total_selected_transcript_compatible,
        page.total_matching,
        page.shown_start,
        page.shown_end
    );
    for lane in &page.lanes {
        out.push_str(&format!(
            "lane\t{}\t{}\t{}\n",
            json(&lane.sample_id)?,
            json(&lane.unit.unit)?,
            json(lane)?
        ));
    }
    for sample in &page.coverage.samples {
        out.push_str(&format!(
            "denominators\t{}\tnull\t{}\n",
            json(&sample.sample_id)?,
            json(sample)?
        ));
        for unit in &sample.excluded_units {
            out.push_str(&format!(
                "excluded\t{}\t{}\t{}\n",
                json(&sample.sample_id)?,
                json(&unit.unit)?,
                json(unit)?
            ));
        }
    }
    Ok(out)
}

#[cfg(test)]
#[path = "pileup_tests.rs"]
mod tests;
