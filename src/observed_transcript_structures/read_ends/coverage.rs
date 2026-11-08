//! Denominators and conditional transcript coverage from saved core projections only.
//! Alternative hypotheses are interval bounds, never fractional read assignments.

use super::*;

#[cfg(test)]
#[path = "coverage_tests.rs"]
mod tests;

pub const COVERAGE_SCHEMA: &str = "gentle.rna_read_end_coverage.v1";

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InputBindings {
    pub report_content_sha256: String,
    pub data_sha256: String,
    pub model_sha256: String,
    pub projection_sha256: String,
}

pub struct SampleInput<'a> {
    pub sample_id: &'a str,
    pub gene_id: &'a str,
    pub bindings: &'a InputBindings,
    pub analysis: &'a ObservedSamAnalysis,
    pub projection: &'a Projection,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ExclusionBucket {
    MultiGene,
    NoCompatibleTranscript,
    NoOptimalHypothesis,
    UnavailableUnit,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Compatibility {
    Unique,
    Shared,
    HypothesisDependent,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExcludedUnit {
    pub unit: ObservedSamQueryUnitKey,
    pub bucket: ExclusionBucket,
    pub reason: String,
}

/// One shared lane-membership rule; incomplete search/views alone do not exclude a unit.
pub fn gene_assignment(
    unit: &UnitProjection,
    unavailable: bool,
    gene: &str,
) -> Result<(), (ExclusionBucket, &'static str)> {
    if unavailable {
        return Err((
            ExclusionBucket::UnavailableUnit,
            "upstream_unit_unavailable",
        ));
    }
    if unit.hypotheses.is_empty() {
        return Err((
            ExclusionBucket::NoOptimalHypothesis,
            "no_optimal_hypothesis",
        ));
    }
    if unit
        .hypotheses
        .iter()
        .any(|h| h.unavailable.is_some() || h.candidates.is_empty())
    {
        return Err((
            ExclusionBucket::NoCompatibleTranscript,
            "incomplete_compatible_candidate_sets",
        ));
    }
    if unit.annotation_gene_ids.len() > 1 {
        return Err((ExclusionBucket::MultiGene, "multiple_annotation_genes"));
    }
    if unit.annotation_gene_ids != [gene] {
        return Err((
            ExclusionBucket::UnavailableUnit,
            "gene_scope_unavailable_or_mismatched",
        ));
    }
    Ok(())
}

pub fn compatibility(unit: &UnitProjection, selected: &str) -> Option<Compatibility> {
    let present = unit
        .hypotheses
        .iter()
        .filter(|h| h.candidates.iter().any(|c| c.transcript_id == selected))
        .count();
    if present == 0 {
        None
    } else if present < unit.hypotheses.len() {
        Some(Compatibility::HypothesisDependent)
    } else if unit.hypotheses.iter().all(|h| h.candidates.len() == 1) {
        Some(Compatibility::Unique)
    } else {
        Some(Compatibility::Shared)
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Fraction {
    pub numerator: usize,
    pub denominator: usize,
    pub value: Option<f64>,
}

impl Fraction {
    fn new(numerator: usize, denominator: usize) -> Self {
        Self {
            numerator,
            denominator,
            value: (denominator != 0).then(|| numerator as f64 / denominator as f64),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MetricSummary {
    pub observed: usize,
    pub right_censored: usize,
    pub missing: usize,
    pub missing_by_reason: BTreeMap<Reason, usize>,
    pub censored_among_observed_and_censored: Fraction,
    pub censored_among_all: Fraction,
    pub exact_quantiles: Vec<(f64, Option<f64>)>,
}

fn summary(values: Vec<&Metric>, probabilities: &[f64]) -> MetricSummary {
    let (mut observed, mut censored, mut missing) = (vec![], 0, BTreeMap::new());
    for value in values {
        match value {
            Metric::Observed { value_bp } => observed.push(*value_bp),
            Metric::RightCensored { .. } => censored += 1,
            Metric::Missing { reason } => *missing.entry(*reason).or_insert(0) += 1,
        }
    }
    observed.sort_unstable();
    let missing_count = missing.values().sum();
    MetricSummary {
        observed: observed.len(),
        right_censored: censored,
        missing: missing_count,
        missing_by_reason: missing,
        censored_among_observed_and_censored: Fraction::new(censored, observed.len() + censored),
        censored_among_all: Fraction::new(censored, observed.len() + censored + missing_count),
        exact_quantiles: probabilities
            .iter()
            .map(|&p| {
                let value = if observed.is_empty() {
                    None
                } else {
                    let index = (observed.len() - 1) as f64 * p;
                    let (a, b) = (index.floor() as usize, index.ceil() as usize);
                    Some(
                        observed[a] as f64 * (1.0 - index.fract())
                            + observed[b] as f64 * index.fract(),
                    )
                };
                (p, value)
            })
            .collect(),
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CoverageSegment {
    pub transcript_interval: Interval,
    pub minimum_conditional_depth: usize,
    pub maximum_conditional_depth: usize,
    pub minimum_fraction: Fraction,
    pub maximum_fraction: Fraction,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SampleCoverage {
    pub sample_id: String,
    pub gene_id: String,
    pub bindings: InputBindings,
    pub declared_units: usize,
    pub assessed_units: usize,
    pub gene_assigned_units: usize,
    pub selected_transcript_compatible_units: usize,
    pub compatibility_counts: BTreeMap<Compatibility, usize>,
    pub read_length: MetricSummary,
    pub utr5: MetricSummary,
    pub utr3: MetricSummary,
    pub excluded_units: Vec<ExcludedUnit>,
    pub coverage: Vec<CoverageSegment>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CoverageReport {
    pub schema: String,
    pub identity_scope: String,
    pub source_id: String,
    pub gene_id: String,
    pub selected_transcript: Option<String>,
    pub quantile_probabilities: Vec<f64>,
    pub samples: Vec<SampleCoverage>,
    pub nonclaims: Vec<String>,
}

fn intersection(
    a: &[Interval],
    b: &[Interval],
    meter: &mut Meter,
) -> Result<Vec<Interval>, EngineError> {
    let (mut i, mut j, mut out) = (0, 0, vec![]);
    while i < a.len() && j < b.len() {
        meter.charge(1)?;
        let (start, end) = (a[i].0.max(b[j].0), a[i].1.min(b[j].1));
        if start < end {
            out.push((start, end));
        }
        if a[i].1 < b[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    Ok(out)
}

fn curves(
    units: &[&UnitProjection],
    selected: &str,
    length: usize,
    meter: &mut Meter,
) -> Result<Vec<CoverageSegment>, EngineError> {
    let mut events: BTreeMap<usize, (isize, isize)> = BTreeMap::new();
    events.insert(0, (0, 0));
    events.insert(length, (0, 0));
    for unit in units {
        let mut lower = None;
        let mut upper = vec![];
        for h in &unit.hypotheses {
            meter.charge(1)?;
            let ranges = h
                .candidates
                .iter()
                .find(|c| c.transcript_id == selected && c.projection_available)
                .map(|c| c.transcript_intervals.clone())
                .unwrap_or_default();
            if ranges.iter().any(|&(a, z)| a >= z || z > length) {
                return Err(fail("coverage_interval"));
            }
            if !ranges.windows(2).all(|w| w[0].1 < w[1].0) {
                return Err(fail("coverage_interval_order"));
            }
            lower = Some(match lower {
                None => ranges.clone(),
                Some(previous) => intersection(&previous, &ranges, meter)?,
            });
            if upper
                .len()
                .checked_add(ranges.len())
                .is_none_or(|n| n > meter.limits.max_intervals)
            {
                return Err(fail("coverage_interval_budget"));
            }
            upper.extend(ranges);
        }
        let upper = meter.merge(upper)?;
        for (kind, ranges) in [lower.unwrap_or_default(), upper].into_iter().enumerate() {
            for (a, z) in ranges {
                meter.charge(1)?;
                let e = events.entry(a).or_default();
                if kind == 0 {
                    e.0 += 1;
                } else {
                    e.1 += 1;
                }
                let e = events.entry(z).or_default();
                if kind == 0 {
                    e.0 -= 1;
                } else {
                    e.1 -= 1;
                }
            }
        }
    }
    let (mut low, mut high, mut previous, mut out) = (0isize, 0isize, 0usize, vec![]);
    for (position, (a, b)) in events {
        meter.charge(1)?;
        if previous < position {
            if low < 0 || high < low {
                return Err(fail("coverage_counts"));
            }
            out.push(CoverageSegment {
                transcript_interval: (previous, position),
                minimum_conditional_depth: low as usize,
                maximum_conditional_depth: high as usize,
                minimum_fraction: Fraction::new(low as usize, units.len()),
                maximum_fraction: Fraction::new(high as usize, units.len()),
            });
        }
        low += a;
        high += b;
        previous = position;
    }
    Ok(out)
}

/// Saved-input composition only. Callers authenticate files; no table/source is reopened here.
pub fn summarize(
    source_id: &str,
    gene_id: &str,
    selected: Option<&str>,
    probabilities: &[f64],
    inputs: &[SampleInput<'_>],
    limits: Limits,
) -> Result<CoverageReport, EngineError> {
    summarize_with_control(
        source_id,
        gene_id,
        selected,
        probabilities,
        inputs,
        limits,
        &mut || Ok(()),
    )
}

pub fn summarize_with_control(
    source_id: &str,
    gene_id: &str,
    selected: Option<&str>,
    probabilities: &[f64],
    inputs: &[SampleInput<'_>],
    limits: Limits,
    control: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<CoverageReport, EngineError> {
    control()?;
    if !token(source_id)
        || !token(gene_id)
        || inputs.is_empty()
        || probabilities.len() > 1000
        || probabilities
            .iter()
            .any(|p| !p.is_finite() || !(0.0..=1.0).contains(p))
        || !probabilities.windows(2).all(|w| w[0] < w[1])
    {
        return Err(fail("coverage_parameters"));
    }
    let mut meter = Meter {
        limits,
        work: 0,
        control: Some(control),
    };
    let mut seen = BTreeSet::new();
    let mut samples = vec![];
    for input in inputs {
        meter.charge(1)?;
        if !token(input.sample_id)
            || input.gene_id != gene_id
            || input.projection.source_id != source_id
            || !seen.insert((input.sample_id, input.gene_id))
        {
            return Err(fail("coverage_sample_gene_analysis"));
        }
        for hash in [
            &input.bindings.report_content_sha256,
            &input.bindings.data_sha256,
            &input.bindings.model_sha256,
            &input.bindings.projection_sha256,
        ] {
            if !sha(hash) {
                return Err(fail("coverage_binding"));
            }
        }
        let declared_gene = if input.analysis.reference_gene.source_id == source_id {
            Some(input.analysis.reference_gene.gene_id.as_str())
        } else {
            input
                .analysis
                .gene_resolutions
                .iter()
                .find(|s| s.source_id == source_id)
                .and_then(|s| s.comparison.as_ref())
                .filter(|c| c.resolution == GeneResolution::Resolved)
                .and_then(|c| c.gene_id.as_deref())
        };
        if declared_gene != Some(gene_id) {
            return Err(fail("coverage_gene_binding"));
        }
        if input.analysis.units.len() > limits.max_units
            || input.projection.units.len() != input.analysis.units.len()
        {
            return Err(fail("coverage_unit_binding"));
        }
        let upstream: BTreeMap<_, _> = input.analysis.units.iter().map(|u| (&u.unit, u)).collect();
        if upstream.len() != input.analysis.units.len() {
            return Err(fail("coverage_unit_binding"));
        }
        let (mut assigned, mut excluded, mut compatible, mut counts) =
            (vec![], vec![], vec![], BTreeMap::new());
        if selected.is_some() {
            for class in [
                Compatibility::Unique,
                Compatibility::Shared,
                Compatibility::HypothesisDependent,
            ] {
                counts.insert(class, 0);
            }
        }
        let mut unit_ids = BTreeSet::new();
        for unit in &input.projection.units {
            meter.charge(1)?;
            if !unit_ids.insert(&unit.unit) {
                return Err(fail("coverage_unit_binding"));
            }
            let source = upstream
                .get(&unit.unit)
                .ok_or_else(|| fail("coverage_unit_binding"))?;
            if let Err((bucket, reason)) =
                gene_assignment(unit, !source.unavailable.is_empty(), gene_id)
            {
                excluded.push(ExcludedUnit {
                    unit: unit.unit.clone(),
                    bucket,
                    reason: reason.into(),
                });
                continue;
            }
            assigned.push(unit);
            if let Some(class) = selected.and_then(|id| compatibility(unit, id)) {
                *counts.entry(class).or_insert(0) += 1;
                compatible.push(unit);
            }
        }
        let coverage = if let Some(id) = selected {
            let model = input
                .projection
                .models
                .iter()
                .find(|m| m.transcript_id == id)
                .ok_or_else(|| fail("coverage_selected_model"))?;
            curves(&compatible, id, model.transcript_length, &mut meter)?
        } else {
            vec![]
        };
        samples.push(SampleCoverage {
            sample_id: input.sample_id.into(),
            gene_id: gene_id.into(),
            bindings: input.bindings.clone(),
            declared_units: input.analysis.units.len(),
            assessed_units: input
                .analysis
                .units
                .iter()
                .filter(|u| u.all_query_views_assessed)
                .count(),
            gene_assigned_units: assigned.len(),
            selected_transcript_compatible_units: compatible.len(),
            compatibility_counts: counts,
            read_length: summary(
                assigned.iter().map(|u| &u.read_length).collect(),
                probabilities,
            ),
            utr5: summary(assigned.iter().map(|u| &u.utr5).collect(), probabilities),
            utr3: summary(assigned.iter().map(|u| &u.utr3).collect(), probabilities),
            excluded_units: excluded,
            coverage,
        });
    }
    Ok(CoverageReport {schema:COVERAGE_SCHEMA.into(),identity_scope:"sam_declared_query_unit".into(),source_id:source_id.into(),
        gene_id:gene_id.into(),selected_transcript:selected.map(str::to_owned),quantile_probabilities:probabilities.to_vec(),samples,
        nonclaims:vec!["Coverage is conditional on model compatibility, not isoform allocation or molecule abundance.".into(),
            "Minimum/maximum depths retain alternative-hypothesis uncertainty; no fractional assignment or smoothing.".into(),
            "Gene and selected-transcript denominators differ; missing values never become zero.".into()]})
}
