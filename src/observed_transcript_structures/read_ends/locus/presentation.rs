//! One-read pages and linked scene projections; filters never recompute input denominators.

use super::*;
use gentle_render::read_end_pileup::{Band, Orientation, Row, RunGlyph, Scene};

pub const VIEW_SCHEMA: &str = "gentle.rna_read_end_locus_view.v1";
pub const PAGE_SCHEMA: &str = "gentle.rna_read_end_locus_page.v1";
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Projection {
    Genome,
    Transcript { model_key: String },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct View {
    pub schema: String,
    pub projection: Projection,
    pub sample_ids: Vec<String>,
    pub class: Option<Class>,
    pub offset: usize,
    pub count: usize,
    pub max_sublanes: usize,
    pub width: u32,
}
impl Default for View {
    fn default() -> Self {
        Self {
            schema: VIEW_SCHEMA.into(),
            projection: Projection::Genome,
            sample_ids: vec![],
            class: None,
            offset: 0,
            count: 100,
            max_sublanes: 2000,
            width: 1400,
        }
    }
}
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UtrCounts {
    pub five: BTreeMap<String, usize>,
    pub three: BTreeMap<String, usize>,
}
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Counters {
    pub sample_id: String,
    pub input_units: usize,
    pub locus_units: usize,
    pub non_admitted_observation_only_units: usize,
    pub visible_units: usize,
    pub same_direction: usize,
    pub opposite_direction: usize,
    pub alternative_direction: usize,
    pub units_with_bases: usize,
    pub units_with_base_alternatives: usize,
    pub possible_terminal_a_runs: usize,
    pub possible_terminal_t_runs: usize,
    pub run_alignment: BTreeMap<String, usize>,
    pub utr_by_model: BTreeMap<String, UtrCounts>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Page {
    pub schema: String,
    pub view: View,
    pub report_sha256: String,
    pub view_sha256: String,
    pub scoring_version: String,
    pub scoring_rules: Vec<String>,
    pub effective_limits: Limits,
    pub work: usize,
    pub models: Vec<Model>,
    pub nonclaims: Vec<String>,
    pub total_input_units: usize,
    pub total_locus_units: usize,
    pub total_matching: usize,
    pub shown_start: usize,
    pub shown_end: usize,
    pub summaries: Vec<Counters>,
    pub source_samples: Vec<SampleSummary>,
    pub reads: Vec<Read>,
}
fn selected<'a>(report: &'a Report, view: &View) -> Result<Option<&'a Model>, EngineError> {
    match &view.projection {
        Projection::Genome => Ok(None),
        Projection::Transcript { model_key } => report
            .models
            .iter()
            .find(|m| m.key == *model_key)
            .map(Some)
            .ok_or_else(|| fail("locus_view_model")),
    }
}
fn assignment<'a>(read: &'a Read, model: &Model) -> Option<&'a Assignment> {
    read.assignments.iter().find(|a| {
        a.model.transcript_id == model.transcript_id
            && a.model.annotation_sha256 == model.annotation_sha256
    })
}
fn direction(read: &Read, model: Option<&Model>) -> Option<bool> {
    let model = model?;
    if read.hypotheses.is_empty()
        || read.hypotheses.iter().any(|h| h.chimeric)
        || read.alignment_directions.len() != 1
    {
        return None;
    }
    Some(read.alignment_directions[0] == model.strand)
}
fn status(metric: &Metric) -> &'static str {
    match metric {
        Metric::Observed { .. } => "observed",
        Metric::RightCensored { .. } => "right_censored",
        Metric::Missing { .. } => "missing",
    }
}
fn metric_status(a: &Assignment, five: bool) -> &'static str {
    let mut values = a.paths.iter().map(|p| if five { &p.utr5 } else { &p.utr3 });
    let Some(first) = values.next() else {
        return "missing";
    };
    if values.all(|m| m == first) {
        status(first)
    } else {
        "alternatives"
    }
}
fn unassigned(read: &Read) -> impl Iterator<Item = &Observation> {
    let admitted: BTreeSet<_> = read
        .hypotheses
        .iter()
        .flat_map(|h| h.observation_ids.iter().map(String::as_str))
        .collect();
    read.observations
        .iter()
        .filter(move |o| !admitted.contains(o.observation_id.as_str()))
}
fn counters(report: &Report, view: &View, visible: &[Read]) -> Result<Vec<Counters>, EngineError> {
    let chosen = selected(report, view)?;
    let mut result = vec![];
    for sample in &report.samples {
        let mut c = Counters {
            sample_id: sample.sample_id.clone(),
            input_units: sample.input_units,
            locus_units: sample.locus_units,
            non_admitted_observation_only_units: sample.non_admitted_observation_only_units,
            visible_units: visible
                .iter()
                .filter(|r| r.sample_id == sample.sample_id)
                .count(),
            ..Default::default()
        };
        for read in report
            .reads
            .iter()
            .filter(|r| r.sample_id == sample.sample_id)
        {
            let reference = chosen.or_else(|| {
                report.models.iter().find(|m| {
                    m.sample_id == read.sample_id
                        && Some(&m.transcript_id) == report.region.anchor_transcript_id.as_ref()
                })
            });
            match direction(read, reference) {
                Some(true) => c.same_direction += 1,
                Some(false) => c.opposite_direction += 1,
                None => c.alternative_direction += 1,
            }
            c.units_with_bases += usize::from(
                read.genomic_basecalls
                    .iter()
                    .any(|e| !e.variants.is_empty()),
            );
            c.units_with_base_alternatives +=
                usize::from(read.genomic_basecalls.iter().any(|e| e.variants.len() > 1));
            let mut runs: BTreeMap<_, BTreeSet<_>> = BTreeMap::new();
            for e in &read.genomic_basecalls {
                for v in &e.variants {
                    for r in &v.runs {
                        if !r.terminal {
                            continue;
                        }
                        runs.entry((r.query_interval, r.base)).or_default().insert(
                            if read.hypotheses.is_empty() {
                                "unavailable"
                            } else if r.aligned_to_reference {
                                "full"
                            } else if r.aligned_query_intervals.is_empty() {
                                "none"
                            } else {
                                "partial"
                            },
                        );
                    }
                }
            }
            for ((_, base), states) in runs {
                if base == 'A' {
                    c.possible_terminal_a_runs += 1;
                } else {
                    c.possible_terminal_t_runs += 1;
                }
                let state = if states.len() == 1 {
                    *states.first().unwrap()
                } else {
                    "alternatives"
                };
                *c.run_alignment.entry(state.into()).or_default() += 1;
            }
            for a in &read.assignments {
                let entry = c.utr_by_model.entry(a.model.key.clone()).or_default();
                *entry.five.entry(metric_status(a, true).into()).or_default() += 1;
                *entry
                    .three
                    .entry(metric_status(a, false).into())
                    .or_default() += 1;
            }
        }
        result.push(c);
    }
    Ok(result)
}
pub fn page(report: &Report, view: &View) -> Result<Page, EngineError> {
    if report.schema != SCHEMA
        || report.scoring_version != RULE_VERSION
        || view.schema != VIEW_SCHEMA
        || !(900..=8000).contains(&view.width)
        || view.count == 0
        || view.count > 2000
        || view.max_sublanes == 0
        || view.max_sublanes > 20_000
        || view.sample_ids.len() > 64
    {
        return Err(fail("locus_view"));
    }
    selected(report, view)?;
    if view
        .sample_ids
        .iter()
        .any(|id| !report.samples.iter().any(|s| s.sample_id == *id))
    {
        return Err(fail("locus_view_sample"));
    }
    let matching: Vec<_> = report
        .reads
        .iter()
        .filter(|r| {
            (view.sample_ids.is_empty() || view.sample_ids.contains(&r.sample_id))
                && view
                    .class
                    .is_none_or(|c| r.assignments.iter().any(|a| a.class == c))
        })
        .collect();
    if view.offset > matching.len() {
        return Err(fail("locus_view_range"));
    }
    let end = view
        .offset
        .checked_add(view.count)
        .ok_or_else(|| fail("locus_view_range"))?
        .min(matching.len());
    let count = matching[view.offset..end]
        .iter()
        .try_fold(0usize, |n, r| {
            let lanes = r.hypotheses.len().checked_add(unassigned(r).count())?;
            n.checked_add(lanes.max(1))
        })
        .ok_or_else(|| fail("locus_sublane_budget"))?;
    if count > view.max_sublanes {
        return Err(fail("locus_sublane_budget"));
    }
    let reads: Vec<_> = matching[view.offset..end]
        .iter()
        .map(|r| (*r).clone())
        .collect();
    Ok(Page {
        schema: PAGE_SCHEMA.into(),
        view: view.clone(),
        report_sha256: pileup::hash(report)?,
        view_sha256: pileup::hash(view)?,
        scoring_version: report.scoring_version.clone(),
        scoring_rules: report.scoring_rules.clone(),
        effective_limits: report.effective_limits,
        work: report.work,
        models: report.models.clone(),
        nonclaims: report.nonclaims.clone(),
        total_input_units: report.samples.iter().map(|s| s.input_units).sum(),
        total_locus_units: report.reads.len(),
        total_matching: matching.len(),
        shown_start: view.offset,
        shown_end: end,
        summaries: counters(report, view, &reads)?,
        source_samples: report.samples.clone(),
        reads,
    })
}
fn json<T: Serialize>(v: &T) -> Result<String, EngineError> {
    serde_json::to_string(v).map_err(|e| fail(&e.to_string()))
}
fn band(interval: Interval, kind: &str, label: String) -> Band {
    Band {
        interval: (interval.0 as u64, interval.1 as u64),
        kind: kind.into(),
        label,
    }
}
fn mapped(range: Interval, model: Option<&Model>) -> Result<Vec<Interval>, EngineError> {
    if let Some(model) = model {
        let map = exon_map(&model.exons, model.strand)?;
        Ok(map
            .iter()
            .filter_map(|e| geometry::project(e, range, model.strand))
            .collect())
    } else {
        Ok(vec![range])
    }
}
pub fn scene(report: &Report, page: &Page) -> Result<Scene, EngineError> {
    if page.report_sha256 != pileup::hash(report)? || page.view_sha256 != pileup::hash(&page.view)?
    {
        return Err(fail("locus_page_binding"));
    }
    let chosen = selected(report, &page.view)?;
    let (start, end, axis) = if let Some(model) = chosen {
        (0, model.exons.iter().map(|i| i.1 - i.0).sum(), "transcript")
    } else {
        (report.window.0, report.window.1, "genomic")
    };
    let mut top_tracks = vec![];
    let mut bottom_tracks = vec![];
    let mut tracks: BTreeMap<_, Vec<_>> = BTreeMap::new();
    for m in &report.models {
        tracks
            .entry((&m.annotation_sha256, &m.transcript_id))
            .or_default()
            .push(m);
    }
    for group in tracks.values() {
        let model = group[0];
        if chosen.is_some_and(|c| {
            c.annotation_sha256 != model.annotation_sha256 || c.transcript_id != model.transcript_id
        }) {
            continue;
        }
        let model_keys: Vec<_> = group.iter().map(|m| &m.key).collect();
        let mut bands = vec![];
        for &exon in &model.exons {
            for i in mapped(exon, chosen)? {
                bands.push(band(
                    i,
                    "exon",
                    format!("{} genomic {:?}", model.transcript_id, exon),
                ));
            }
        }
        let row = Row {
            label: format!(
                "{} {}",
                if chosen.is_some() || model.strand == '+' {
                    "->"
                } else {
                    "<-"
                },
                model.transcript_id
            ),
            details: json(&serde_json::json!({"model_keys":model_keys,"sample_bindings":group}))?,
            off_axis: format!(
                "genes {:?}; genomic strand {}; annotation {} {}; model {:?}; model keys {:?}",
                model.gene_ids,
                model.strand,
                model.annotation_release,
                model.annotation_sha256,
                model.model_sha256,
                model_keys
            ),
            bands,
            orientation: None,
            runs: vec![],
            sublanes: vec![],
        };
        if chosen.is_some() || model.strand == '+' {
            top_tracks.push(row);
        } else {
            bottom_tracks.push(row);
        }
    }
    let mut rows = vec![];
    for read in &page.reads {
        let mut sublanes = vec![];
        let mut off = vec![format!(
            "read {}; molecule strand {} {:?}; supplementary={}; raw directions {:?}",
            read.read_key,
            read.molecule_strand.status,
            read.molecule_strand.strand,
            read.supplementary_present,
            read.alignment_directions
        )];
        off.push(format!("membership: {}", json(&read.membership)?));
        off.push(format!("conditional model support {:?}; search_complete={} views_complete={} all_views_assessed={}; not physical RNA-origin proof",
            read.molecule_strand.reasons,read.search_declared_complete,read.views_declared_complete,read.all_query_views_assessed));
        for h in &read.hypotheses {
            let mut bands = vec![];
            for o in read
                .observations
                .iter()
                .filter(|o| h.observation_ids.contains(&o.observation_id))
            {
                if o.contig != report.region.contig {
                    off.push(format!("H{} off-contig {}", h.index, o.contig));
                    continue;
                }
                if let Some(clipping) = &o.clipping {
                    off.push(format!("H{} saved clipping {}", h.index, json(clipping)?));
                }
                for b in &o.blocks {
                    for i in mapped(b.target, chosen)? {
                        bands.push(band(
                            i,
                            "aligned",
                            format!("H{} {} query {:?}", h.index, o.observation_id, b.query),
                        ));
                    }
                }
                for &skip in &o.skips {
                    if chosen.is_none() {
                        bands.push(band(
                            skip,
                            "splice_n",
                            format!("H{} saved N {:?}", h.index, skip),
                        ));
                    } else {
                        off.push(format!(
                            "H{} genomic N {:?} retained; transcript projection omits introns",
                            h.index, skip
                        ));
                    }
                }
                if o.blocks
                    .iter()
                    .any(|b| b.target.0 < report.window.0 || b.target.1 > report.window.1)
                {
                    off.push(format!(
                        "H{} alignment continues outside locus; full coordinates retained",
                        h.index
                    ));
                    if chosen.is_none() {
                        if o.blocks.iter().any(|b| b.target.0 < start) {
                            bands.push(band(
                                (start, start + 1),
                                "continuation",
                                "<- alignment continues outside window".into(),
                            ));
                        }
                        if o.blocks.iter().any(|b| b.target.1 > end) {
                            bands.push(band(
                                (end - 1, end),
                                "continuation",
                                "alignment continues outside window ->".into(),
                            ));
                        }
                    }
                }
            }
            if chosen.is_none() {
                let mut segments: Vec<_> = read
                    .observations
                    .iter()
                    .filter(|o| {
                        h.observation_ids.contains(&o.observation_id)
                            && o.contig == report.region.contig
                    })
                    .filter_map(|o| {
                        Some((
                            o.blocks.iter().map(|b| b.target.0).min()?,
                            o.blocks.iter().map(|b| b.target.1).max()?,
                        ))
                    })
                    .collect();
                segments.sort_unstable();
                for pair in segments.windows(2) {
                    if pair[0].1 < pair[1].0 {
                        bands.push(band(
                            (pair[0].1, pair[1].0),
                            "segment_link",
                            "saved separate segments; not an inferred splice".into(),
                        ));
                    }
                }
            }
            sublanes.push(bands);
        }
        if read.hypotheses.is_empty() {
            off.push(
                "No legal path; observed geometry remains inspectable, not a fabricated alignment"
                    .into(),
            );
        }
        for o in unassigned(read) {
            off.push(format!(
                "Saved observation {} not admitted to a legal hypothesis; contig {}; clipping {:?}",
                o.observation_id, o.contig, o.clipping
            ));
            if o.contig != report.region.contig {
                continue;
            }
            let mut raw = vec![];
            for b in &o.blocks {
                for i in mapped(b.target, chosen)? {
                    raw.push(band(
                        i,
                        "unassessed",
                        format!(
                            "unassessed raw observation {} query {:?}",
                            o.observation_id, b.query
                        ),
                    ));
                }
            }
            if chosen.is_none() {
                for &skip in &o.skips {
                    raw.push(band(
                        skip,
                        "deletion_or_gap",
                        format!("saved N {:?}, unassessed hypothesis", skip),
                    ));
                }
            }
            sublanes.push(raw);
        }
        for a in &read.assignments {
            off.push(format!(
                "{} {:?}; rank {}; preferred={}; reasons {:?}; conflicts {:?}",
                a.model.transcript_id,
                a.class,
                json(&a.score)?,
                a.preferred,
                a.evidence,
                a.conflicts
            ));
        }
        let evidence = chosen
            .and_then(|m| assignment(read, m))
            .map(|a| a.paths.iter().map(|p| &p.oriented).collect::<Vec<_>>())
            .unwrap_or_else(|| read.genomic_basecalls.iter().collect());
        let mut glyphs = vec![];
        let mut seen = BTreeSet::new();
        for e in evidence {
            for v in &e.variants {
                for r in &v.runs {
                    if !seen.insert((
                        v.sequence_sha256.clone(),
                        r.query_interval,
                        r.base,
                        r.aligned_to_reference,
                        r.reference_anchor,
                    )) {
                        continue;
                    }
                    if let Some(anchor) = r.reference_anchor {
                        if chosen.is_none() && !(start..end).contains(&anchor) {
                            off.push(format!("Original {}{} query {:?}; genomic attachment {}; outside the locus window; retained in JSON/TSV",r.length,r.base,r.query_interval,anchor));
                            continue;
                        }
                        let anchors = mapped((anchor, anchor + 1), chosen)?;
                        let anchor = anchors.first().map(|i| i.0);
                        if let Some(anchor) = anchor {
                            let label = format!(
                                "{}{} query {:?}; terminal={}; reference overlap={} / {}; read ends {:?}; transcript ends {:?}; attachment only",
                                r.length,
                                r.base,
                                r.query_interval,
                                r.terminal,
                                r.aligned_query_intervals
                                    .iter()
                                    .map(|i| i.1 - i.0)
                                    .sum::<usize>(),
                                r.length,
                                r.read_ends,
                                r.transcript_ends
                            );
                            glyphs.push(RunGlyph {
                                base: r.base,
                                length: r.length,
                                reference_anchor: anchor as u64,
                                aligned_to_reference: r.aligned_to_reference,
                                label,
                            });
                        } else {
                            off.push(format!(
                                "Unprojectable original run {}{} query {:?}; retained in JSON/TSV",
                                r.length, r.base, r.query_interval
                            ));
                        }
                    } else {
                        off.push(format!("Original {}{} query {:?}; no reference anchor or admitted coverage; retained in JSON/TSV",r.length,r.base,r.query_interval));
                    }
                }
            }
        }
        let orientation = chosen.map(|m| match direction(read, Some(m)) {
            Some(true) => Orientation::Same,
            Some(false) => Orientation::Opposite,
            None => Orientation::Alternatives,
        });
        rows.push(Row {
            label: format!("{} / {}", read.sample_id, read.unit.unit_id),
            details: json(read)?,
            off_axis: off.join("; "),
            bands: vec![],
            orientation,
            runs: glyphs,
            sublanes,
        });
    }
    if chosen.is_some() {
        rows.sort_by_key(|r| r.orientation);
    }
    Ok(Scene {
        contig: chosen.map_or_else(
            || report.region.contig.clone(),
            |m| format!("transcript {}", m.transcript_id),
        ),
        start: start as u64,
        end: end as u64,
        strand: 1,
        width: page.view.width,
        header: vec![
            format!(
                "Exon locus {}; reads {}..{} / {} matching / {} locus / {} declared; scores {}",
                report.source_id,
                page.shown_start,
                page.shown_end,
                page.total_matching,
                page.total_locus_units,
                page.total_input_units,
                report.scoring_version
            ),
            format!(
                "Counters are conditional per model, not additive molecules: {}",
                json(&page.summaries)?
            ),
        ],
        model: vec![],
        rows,
        top_tracks,
        bottom_tracks,
        axis_kind: Some(axis.into()),
    })
}
pub fn tsv(page: &Page) -> Result<String, EngineError> {
    let mut out = format!(
        "# report={} view={} input={} locus={} matching={} shown={}..{}\nrecord_type\tdata_json\n",
        page.report_sha256,
        page.view_sha256,
        page.total_input_units,
        page.total_locus_units,
        page.total_matching,
        page.shown_start,
        page.shown_end
    );
    out.push_str(&format!(
        "contract\t{}\n",
        json(&(
            page.scoring_version.clone(),
            page.effective_limits,
            page.work,
            &page.scoring_rules,
            &page.models,
            &page.nonclaims
        ))?
    ));
    for r in &page.reads {
        out.push_str(&format!("read\t{}\n", json(r)?));
    }
    for s in &page.summaries {
        out.push_str(&format!("sample_summary\t{}\n", json(s)?));
    }
    for s in &page.source_samples {
        out.push_str(&format!("source_sample\t{}\n", json(s)?));
    }
    Ok(out)
}
