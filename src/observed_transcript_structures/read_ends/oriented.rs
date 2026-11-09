//! Original-sequencing-orientation basecalls, literal terminal runs and mRNA-relative lanes.
//! A supplied basecall document declares the original orientation; it is not a tail caller.

use super::*;

pub const BASES_SCHEMA: &str = "gentle.rna_read_end_original_bases.v1";
pub(super) const MAX_RUN_RECORDS: usize = 100_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Policy {
    pub window_bp: usize,
    pub minimum_run_bp: usize,
}
impl Default for Policy {
    fn default() -> Self {
        Self {
            window_bp: 60,
            minimum_run_bp: 10,
        }
    }
}

/// Explicit unit association, never a join inferred from a QNAME or a current experiment.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OriginalRead {
    pub unit: ObservedSamQueryUnitKey,
    pub read_id: String,
    /// Full original-query DNA including unavailable-to-SAM hard clips; never a reconstructed SEQ.
    pub sequence: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Document {
    pub schema: String,
    pub analysis_context_sha256: String,
    pub orientation: String,
    pub provenance: String,
    pub reads: Vec<OriginalRead>,
}

#[derive(Clone, Copy)]
pub struct Input<'a> {
    pub document: &'a Document,
    pub file_sha256: &'a str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Direction {
    Same,
    Opposite,
    Alternatives,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum End {
    FivePrime,
    ThreePrime,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Run {
    pub base: char,
    pub length: usize,
    /// Zero-based half-open coordinates before any alignment reverse complementation.
    pub query_interval: Interval,
    pub terminal: bool,
    pub read_ends: Vec<End>,
    /// None means the hypothesis cannot select one mRNA-relative orientation.
    pub transcript_ends: Option<Vec<End>>,
    pub aligned_to_reference: bool,
    pub aligned_query_intervals: Vec<Interval>,
    pub reference_intervals: Vec<Interval>,
    /// An attachment point, never an extrapolated alignment for clipped bases.
    pub reference_anchor: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Variant {
    pub read_id: String,
    pub sequence_sha256: String,
    pub query_length: usize,
    pub runs: Vec<Run>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Evidence {
    pub read_key: String,
    pub direction: Direction,
    pub availability: String,
    pub original_bases_sha256: Option<String>,
    pub provenance: Option<String>,
    pub policy: Policy,
    pub variants: Vec<Variant>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Summary {
    pub gene_assigned_units: usize,
    pub same_direction_units: usize,
    pub opposite_direction_units: usize,
    pub alternative_direction_units: usize,
    pub units_with_bound_bases: usize,
    pub units_with_basecall_alternatives: usize,
    pub units_with_terminal_a: usize,
    pub units_with_terminal_t: usize,
    pub units_with_unaligned_runs: usize,
}

fn literals(sequence: &str, policy: Policy, max_runs: usize) -> Result<Vec<Run>, EngineError> {
    let bytes = sequence.as_bytes();
    let mut runs = vec![];
    let mut start = 0;
    while start < bytes.len() {
        let base = bytes[start].to_ascii_uppercase();
        let mut end = start + 1;
        while end < bytes.len() && bytes[end].to_ascii_uppercase() == base {
            end += 1;
        }
        if matches!(base, b'A' | b'T') && end - start >= policy.minimum_run_bp {
            let mut read_ends = vec![];
            if start < policy.window_bp {
                read_ends.push(End::FivePrime);
            }
            if end > bytes.len().saturating_sub(policy.window_bp) {
                read_ends.push(End::ThreePrime);
            }
            if runs.len() >= max_runs {
                return Err(fail("terminal_run_record_budget"));
            }
            runs.push(Run {
                base: base as char,
                length: end - start,
                query_interval: (start, end),
                terminal: !read_ends.is_empty(),
                read_ends,
                transcript_ends: None,
                aligned_to_reference: false,
                aligned_query_intervals: vec![],
                reference_intervals: vec![],
                reference_anchor: None,
            });
        }
        start = end;
    }
    Ok(runs)
}

/// Validate the authenticated document once; do not choose among conflicting sequence variants.
pub(super) fn prepare(
    input: Option<Input<'_>>,
    analysis: &ObservedSamAnalysis,
    policy: Policy,
    meter: &mut Meter<'_>,
) -> Result<BTreeMap<ObservedSamQueryUnitKey, Vec<Variant>>, EngineError> {
    if policy.window_bp == 0
        || policy.window_bp > 1_000_000
        || policy.minimum_run_bp == 0
        || policy.minimum_run_bp > policy.window_bp
    {
        return Err(fail("terminal_run_policy"));
    }
    let Some(input) = input else {
        return Ok(BTreeMap::new());
    };
    let document = input.document;
    if document.schema != BASES_SCHEMA
        || document.analysis_context_sha256 != analysis.context_sha256
        || document.orientation != "original_sequencing_5_to_3"
        || !sha(input.file_sha256)
        || !identity(&document.provenance)
        || document.provenance.trim().is_empty()
        || document.reads.len() > meter.limits.max_observations
    {
        return Err(fail("original_bases_binding"));
    }
    let units: BTreeMap<_, _> = analysis.units.iter().map(|u| (&u.unit, u)).collect();
    let mut seen = BTreeSet::new();
    let mut out: BTreeMap<_, Vec<Variant>> = BTreeMap::new();
    let mut remaining_runs = meter.limits.max_observations.min(MAX_RUN_RECORDS);
    for read in &document.reads {
        meter.charge(
            read.sequence
                .len()
                .checked_add(1)
                .ok_or_else(|| fail("original_bases_budget"))?,
        )?;
        let unit = units
            .get(&read.unit)
            .ok_or_else(|| fail("original_bases_unit"))?;
        if unit.unclipped_axis_length != Some(read.sequence.len())
            || read.sequence.is_empty()
            || unit.rows.is_empty()
            || unit.rows.iter().any(|r| r.qname != read.read_id)
            || !read
                .sequence
                .bytes()
                .all(|b| b"ACGTNRYKMSWBDHVacgtnrykmswbdhv".contains(&b))
        {
            return Err(fail("original_bases_axis_or_read"));
        }
        let sequence_sha256 = crate::digest_utils::sha256_prefixed_bytes(read.sequence.as_bytes());
        if !seen.insert((read.unit.clone(), sequence_sha256.clone())) {
            return Err(fail("original_bases_duplicate"));
        }
        let runs = literals(&read.sequence, policy, remaining_runs)?;
        remaining_runs -= runs.len();
        out.entry(read.unit.clone()).or_default().push(Variant {
            read_id: read.read_id.clone(),
            sequence_sha256,
            query_length: read.sequence.len(),
            runs,
        });
    }
    for variants in out.values_mut() {
        variants.sort_by(|a, b| a.sequence_sha256.cmp(&b.sequence_sha256));
    }
    if out.values().any(|v| v.len() > 64) {
        return Err(fail("original_bases_variant_budget"));
    }
    Ok(out)
}

/// Relative query/model geometry only, never pore travel or physical RNA strand.
pub fn direction(strand: char, observations: &[pileup::Observation]) -> Direction {
    let values: BTreeSet<_> = observations
        .iter()
        .map(|o| o.query_to_reference_strand == strand)
        .collect();
    match (values.len(), values.first()) {
        (1, Some(true)) => Direction::Same,
        (1, Some(false)) => Direction::Opposite,
        _ => Direction::Alternatives,
    }
}

pub(super) fn project(
    sample: &str,
    unit: &ObservedSamQueryUnitKey,
    strand: char,
    observations: &[pileup::Observation],
    variants: &[Variant],
    input: Option<Input<'_>>,
    policy: Policy,
    meter: &mut Meter<'_>,
) -> Result<Evidence, EngineError> {
    let direction = direction(strand, observations);
    let mut variants = variants.to_vec();
    for variant in &mut variants {
        for run in &mut variant.runs {
            run.transcript_ends =
                (run.terminal && direction != Direction::Alternatives).then(|| {
                    run.read_ends
                        .iter()
                        .map(|end| match (direction, end) {
                            (Direction::Opposite, End::FivePrime) => End::ThreePrime,
                            (Direction::Opposite, End::ThreePrime) => End::FivePrime,
                            (_, end) => *end,
                        })
                        .collect()
                });
            let mut anchors = vec![];
            for observation in observations {
                for block in &observation.blocks {
                    meter.charge(1)?;
                    let (a, z) = block.query;
                    if a >= z || block.target.1.checked_sub(block.target.0) != Some(z - a) {
                        return Err(fail("original_bases_block"));
                    }
                    let (l, r) = (a.max(run.query_interval.0), z.min(run.query_interval.1));
                    if l < r {
                        run.aligned_query_intervals.push((l, r));
                        run.reference_intervals.push(
                            if observation.query_to_reference_strand == '+' {
                                (block.target.0 + l - a, block.target.0 + r - a)
                            } else {
                                (block.target.1 - (r - a), block.target.1 - (l - a))
                            },
                        );
                    }
                    let before = run.query_interval.1 <= a;
                    let distance = if before {
                        a - run.query_interval.1
                    } else {
                        run.query_interval.0.saturating_sub(z)
                    };
                    let anchor = if before == (observation.query_to_reference_strand == '+') {
                        block.target.0
                    } else {
                        block.target.1 - 1
                    };
                    anchors.push((distance, anchor));
                }
            }
            run.aligned_query_intervals.sort_unstable();
            let mut merged: Vec<Interval> = vec![];
            for &(a, z) in &run.aligned_query_intervals {
                if let Some(last) = merged.last_mut()
                    && a <= last.1
                {
                    last.1 = last.1.max(z);
                } else {
                    merged.push((a, z));
                }
            }
            run.aligned_to_reference = merged.len() == 1 && merged[0] == run.query_interval;
            run.aligned_query_intervals = merged;
            run.reference_intervals.sort_unstable();
            run.reference_intervals.dedup();
            run.reference_anchor = run
                .reference_intervals
                .first()
                .map(|&(a, z)| a + (z - a) / 2)
                .or_else(|| anchors.into_iter().min().map(|(_, anchor)| anchor));
        }
    }
    Ok(Evidence {
        read_key: pileup::hash(&(sample, unit))?,
        direction,
        availability: match variants.len() {
            0 => "no_bound_basecalls",
            1 => "available",
            _ => "alternatives",
        }
        .into(),
        original_bases_sha256: input.map(|i| i.file_sha256.to_owned()),
        provenance: input.map(|i| i.document.provenance.clone()),
        policy,
        variants,
    })
}

/// Counts are per declared unit, not per hypothesis/variant; A/T counters can overlap.
pub fn summarize(lanes: &[pileup::Lane]) -> Summary {
    let mut out = Summary {
        gene_assigned_units: lanes.len(),
        ..Summary::default()
    };
    for lane in lanes {
        let evidence: Vec<_> = lane
            .hypotheses
            .iter()
            .filter_map(|h| h.oriented.as_ref())
            .collect();
        let directions: BTreeSet<_> = evidence.iter().map(|e| e.direction).collect();
        match (directions.len(), directions.first()) {
            (1, Some(Direction::Same)) => out.same_direction_units += 1,
            (1, Some(Direction::Opposite)) => out.opposite_direction_units += 1,
            _ => out.alternative_direction_units += 1,
        }
        out.units_with_bound_bases += usize::from(evidence.iter().any(|e| !e.variants.is_empty()));
        out.units_with_basecall_alternatives +=
            usize::from(evidence.iter().any(|e| e.variants.len() > 1));
        for (base, count) in [
            ('A', &mut out.units_with_terminal_a),
            ('T', &mut out.units_with_terminal_t),
        ] {
            *count += usize::from(
                evidence
                    .iter()
                    .flat_map(|e| &e.variants)
                    .flat_map(|v| &v.runs)
                    .any(|r| r.terminal && r.base == base),
            );
        }
        out.units_with_unaligned_runs += usize::from(
            evidence
                .iter()
                .flat_map(|e| &e.variants)
                .flat_map(|v| &v.runs)
                .any(|r| !r.aligned_to_reference),
        );
    }
    out
}

#[cfg(test)]
#[path = "oriented_tests.rs"]
mod tests;
