//! Pure, bounded synonymous motif/GC design, not a general DNA Chisel implementation.
//!
//! The caller supplies an explicit complete codon mapping and expected protein.
//! DNA is uppercase, linear, with one forward CDS and frozen ATG/terminal stop.
//! No files, state, network, serialization, random numbers or external tools are used.

use std::collections::BTreeMap;

mod gc;
pub use gc::{GcAssessment, GcBounds, GcMeasurement, assess_gc};
use gc::{GcViolation, PreparedGc};
mod gc_window;
use gc_window::PreparedWindowGc;
pub use gc_window::{WindowGcAssessment, WindowGcBounds, WindowGcMeasurement, assess_window_gc};

pub const ALGORITHM: &str = "synonymous_full_enumeration_v1";
pub const CONFLICT_ALGORITHM: &str = "synonymous_conflict_search_v1";
pub const GC_ALGORITHM: &str = "synonymous_full_enumeration_gc_v1";
pub const GC_CONFLICT_ALGORITHM: &str = "synonymous_conflict_search_gc_v1";
pub const WINDOW_GC_ALGORITHM: &str = "synonymous_full_enumeration_window_gc_v1";
pub const WINDOW_GC_CONFLICT_ALGORITHM: &str = "synonymous_conflict_search_window_gc_v1";
pub const MAX_SEQUENCE_BP: usize = 12_000;
pub const MAX_MOTIFS: usize = 16;
pub const MAX_MOTIF_BP: usize = 32;
pub const MAX_PROTECTED_INTERVALS: usize = 256;
pub const MAX_EVALUATIONS: u64 = 100_000;
pub const MAX_TOTAL_MOTIF_WORK: u64 = 50_000_000;
pub const MAX_REPORTED_MATCHES: usize = 4_096;

/// Enumeration remains the compatibility default; conflict search is explicit.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SearchStrategy {
    #[default]
    FullEnumeration,
    ConflictDirected,
}

impl SearchStrategy {
    pub fn algorithm(self) -> &'static str {
        match self {
            Self::FullEnumeration => ALGORITHM,
            Self::ConflictDirected => CONFLICT_ALGORITHM,
        }
    }

    pub fn algorithm_with_gc(self, gc_requested: bool) -> &'static str {
        match (self, gc_requested) {
            (_, false) => self.algorithm(),
            (Self::FullEnumeration, true) => GC_ALGORITHM,
            (Self::ConflictDirected, true) => GC_CONFLICT_ALGORITHM,
        }
    }

    pub fn algorithm_with_constraints(
        self,
        gc_requested: bool,
        window_requested: bool,
    ) -> &'static str {
        match (self, window_requested) {
            (_, false) => self.algorithm_with_gc(gc_requested),
            (Self::FullEnumeration, true) => WINDOW_GC_ALGORITHM,
            (Self::ConflictDirected, true) => WINDOW_GC_CONFLICT_ALGORITHM,
        }
    }
}

/// All coordinates are sequence-local, zero-based, half-open.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Interval {
    pub start: usize,
    pub end: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Both,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Motif {
    pub pattern: String,
    pub strand: Strand,
}

/// Elongation mapping; '*' denotes stop. Initiation is not inferred from it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneticCode {
    pub id: usize,
    pub codons: BTreeMap<[u8; 3], u8>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Input {
    pub sequence: String,
    pub cds: Interval,
    /// Includes initial M, excludes terminal stop.
    pub protein: String,
    pub code: GeneticCode,
    pub protected: Vec<Interval>,
    pub motifs: Vec<Motif>,
    pub gc_content: Option<GcBounds>,
    pub gc_window: Option<WindowGcBounds>,
    pub max_evaluations: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Invalid,
    Unsupported,
    Feasible,
    SearchExhausted,
    ProvenInfeasible,
    Cancelled,
}

/// A palindromic match requested on both strands is recorded once, as Both.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Match {
    pub motif_index: usize,
    pub interval: Interval,
    pub strand: Strand,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Edit {
    pub position: usize,
    pub before: u8,
    pub after: u8,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Outcome {
    pub status: Status,
    pub reason: String,
    pub sequence: Option<String>,
    pub edits: Vec<Edit>,
    pub initial_matches: Vec<Match>,
    pub evaluated_candidates: u64,
    pub effective_evaluation_budget: u64,
    /// None means overflow, not an infeasible or infinite search space.
    pub search_space: Option<u64>,
    pub optimization_complete: bool,
    pub minimum_edits_proven: bool,
}

impl Outcome {
    fn empty(status: Status, reason: impl Into<String>) -> Self {
        Self {
            status,
            reason: reason.into(),
            sequence: None,
            edits: vec![],
            initial_matches: vec![],
            evaluated_candidates: 0,
            effective_evaluation_budget: 0,
            search_space: None,
            optimization_complete: false,
            minimum_edits_proven: false,
        }
    }
}

struct Prepared<'a> {
    input: &'a Input,
    frozen: Vec<bool>,
    patterns: Vec<(Vec<u8>, Vec<u8>, Strand)>,
    choices: Vec<(usize, Vec<[u8; 3]>)>,
    budget: u64,
    space: Option<u64>,
    gc: Option<PreparedGc>,
    gc_window: Option<PreparedWindowGc>,
}

struct Violations {
    matches: Vec<Match>,
    gc: Option<GcViolation>,
    gc_window: Option<(Interval, GcViolation)>,
}

#[derive(Clone, Copy)]
enum Conflict {
    Motif(Interval),
    Gc(GcViolation),
    WindowGc(Interval, GcViolation),
}

impl Violations {
    fn satisfied(&self) -> bool {
        self.matches.is_empty() && self.gc.is_none() && self.gc_window.is_none()
    }

    fn first_conflict(&self) -> Option<Conflict> {
        self.matches
            .first()
            .map(|hit| Conflict::Motif(hit.interval))
            .or_else(|| self.gc.map(Conflict::Gc))
            .or_else(|| {
                self.gc_window
                    .map(|(interval, direction)| Conflict::WindowGc(interval, direction))
            })
    }
}

fn cancelled(mut outcome: Outcome) -> Outcome {
    outcome.status = Status::Cancelled;
    outcome.reason = "externally_cancelled".into();
    outcome.sequence = None;
    outcome.edits.clear();
    outcome.optimization_complete = false;
    outcome.minimum_edits_proven = false;
    outcome
}

fn mask(base: u8) -> Option<u8> {
    Some(match base {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => return None,
    })
}

fn complement(m: u8) -> u8 {
    ((m & 1) << 3) | ((m & 2) << 1) | ((m & 4) >> 1) | ((m & 8) >> 3)
}

fn translate(sequence: &[u8], cds: Interval, code: &GeneticCode) -> Option<Vec<u8>> {
    sequence
        .get(cds.start..cds.end)?
        .chunks_exact(3)
        .map(|c| code.codons.get(&[c[0], c[1], c[2]]).copied())
        .collect()
}

fn prepare(input: &Input) -> Result<Prepared<'_>, Outcome> {
    let fail = |reason: &str| Outcome::empty(Status::Invalid, reason);
    let n = input.sequence.len();
    if n == 0
        || n > MAX_SEQUENCE_BP
        || !input
            .sequence
            .bytes()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
    {
        return Err(fail(
            "sequence_requires_1_to_12000_uppercase_unambiguous_bases",
        ));
    }
    if input.code.id != 1 {
        return Err(Outcome::empty(
            Status::Unsupported,
            "only_standard_genetic_code_1_is_supported",
        ));
    }
    if input.code.codons.len() != 64
        || input.code.codons.iter().any(|(c, aa)| {
            !c.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
                || !b"ACDEFGHIKLMNPQRSTVWY*".contains(aa)
        })
    {
        return Err(fail("explicit_complete_64_codon_mapping_required"));
    }
    let cds = input.cds;
    if cds.start >= cds.end
        || cds.end > n
        || cds.end - cds.start < 6
        || (cds.end - cds.start) % 3 != 0
    {
        return Err(fail("one_complete_forward_cds_required"));
    }
    if let Some(bounds) = input.gc_content {
        bounds.validate().map_err(fail)?;
    }
    if let Some(bounds) = input.gc_window {
        bounds.validate(cds.end - cds.start).map_err(fail)?;
    }
    if &input.sequence.as_bytes()[cds.start..cds.start + 3] != b"ATG" {
        return Err(Outcome::empty(
            Status::Unsupported,
            "literal_ATG_start_required",
        ));
    }
    let protein = translate(input.sequence.as_bytes(), cds, &input.code)
        .ok_or_else(|| fail("unmapped_codon"))?;
    if protein.last() != Some(&b'*')
        || protein[..protein.len() - 1].contains(&b'*')
        || protein[0] != b'M'
    {
        return Err(fail("terminal_stop_without_internal_stops_required"));
    }
    if protein[..protein.len() - 1] != *input.protein.as_bytes() {
        return Err(fail("supplied_protein_translation_mismatch"));
    }
    if input.protected.len() > MAX_PROTECTED_INTERVALS
        || input.max_evaluations == 0
        || input.max_evaluations > MAX_EVALUATIONS
    {
        return Err(fail("protection_or_candidate_budget_limit_exceeded"));
    }
    let mut frozen = vec![true; n];
    frozen[cds.start + 3..cds.end - 3].fill(false);
    for interval in &input.protected {
        if interval.start >= interval.end || interval.end > n {
            return Err(fail("protected_interval_out_of_bounds"));
        }
        frozen[interval.start..interval.end].fill(true);
    }
    if input.motifs.is_empty() && input.gc_content.is_none() && input.gc_window.is_none() {
        return Err(fail("at_least_one_motif_or_explicit_gc_bounds_required"));
    }
    if input.motifs.len() > MAX_MOTIFS {
        return Err(fail("one_to_16_motifs_required"));
    }
    let mut patterns = vec![];
    // Also bound full-DNA copying, frozen-base checks and two translations.
    let mut work = n as u64 * 4;
    if input.gc_content.is_some() {
        work += (cds.end - cds.start) as u64;
    }
    if let Some(bounds) = input.gc_window {
        work +=
            (cds.end - cds.start + 1) as u64 + (cds.end - cds.start - bounds.window_bp + 1) as u64;
    }
    for motif in &input.motifs {
        if motif.pattern.is_empty() || motif.pattern.len() > MAX_MOTIF_BP {
            return Err(fail("motif_length_requires_1_to_32_bases"));
        }
        let forward = motif
            .pattern
            .bytes()
            .map(mask)
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| fail("uppercase_finite_IUPAC_pattern_required"))?;
        let reverse = forward
            .iter()
            .rev()
            .map(|m| complement(*m))
            .collect::<Vec<_>>();
        let strands = if motif.strand == Strand::Both && forward != reverse {
            2
        } else {
            1
        };
        work += n.saturating_sub(forward.len() - 1) as u64 * forward.len() as u64 * strands;
        if patterns
            .iter()
            .any(|p| p == &(forward.clone(), reverse.clone(), motif.strand))
        {
            return Err(fail("duplicate_motif_specification"));
        }
        patterns.push((forward, reverse, motif.strand));
    }
    let budget = input
        .max_evaluations
        .min(MAX_TOTAL_MOTIF_WORK / work.max(1));
    if budget == 0 {
        return Err(fail("motif_work_admission_limit_exceeded"));
    }
    let mut choices = vec![];
    let mut space = Some(1u64);
    for start in (cds.start + 3..cds.end - 3).step_by(3) {
        let original: [u8; 3] = input.sequence.as_bytes()[start..start + 3]
            .try_into()
            .unwrap();
        let aa = input.code.codons[&original];
        let mut alternatives = input
            .code
            .codons
            .iter()
            .filter_map(|(codon, residue)| {
                (*residue == aa && (0..3).all(|k| !frozen[start + k] || codon[k] == original[k]))
                    .then_some(*codon)
            })
            .collect::<Vec<_>>();
        alternatives.sort_by_key(|c| ((0..3).filter(|k| c[*k] != original[*k]).count(), *c));
        space = space.and_then(|s| s.checked_mul(alternatives.len() as u64));
        if alternatives.len() > 1 {
            choices.push((start, alternatives));
        }
    }
    let gc = input
        .gc_content
        .map(|_| PreparedGc::new(input, &choices))
        .transpose()
        .map_err(|reason| fail(&reason))?;
    Ok(Prepared {
        input,
        frozen,
        patterns,
        choices,
        budget,
        space,
        gc,
        gc_window: input.gc_window.map(|_| PreparedWindowGc::new(input)),
    })
}

impl Prepared<'_> {
    fn matches(&self, sequence: &[u8]) -> Result<Vec<Match>, String> {
        let mut hits = vec![];
        for (motif_index, (forward, reverse, strand)) in self.patterns.iter().enumerate() {
            if forward.len() > sequence.len() {
                continue;
            }
            for start in 0..=sequence.len() - forward.len() {
                let window = &sequence[start..start + forward.len()];
                let fits = |pattern: &[u8]| {
                    window
                        .iter()
                        .zip(pattern)
                        .all(|(b, m)| mask(*b).is_some_and(|b| b & m != 0))
                };
                let f = *strand != Strand::Reverse && fits(forward);
                let r = *strand != Strand::Forward && fits(reverse);
                let directions = if *strand == Strand::Both && forward == reverse && f {
                    vec![Strand::Both]
                } else {
                    [f.then_some(Strand::Forward), r.then_some(Strand::Reverse)]
                        .into_iter()
                        .flatten()
                        .collect()
                };
                for direction in directions {
                    if hits.len() == MAX_REPORTED_MATCHES {
                        return Err("match_diagnostic_limit_exceeded".into());
                    }
                    hits.push(Match {
                        motif_index,
                        interval: Interval {
                            start,
                            end: start + forward.len(),
                        },
                        strand: direction,
                    });
                }
            }
        }
        Ok(hits)
    }

    fn validate(&self, sequence: &[u8]) -> Result<Violations, String> {
        let original = self.input.sequence.as_bytes();
        if sequence.len() != original.len()
            || !sequence
                .iter()
                .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        {
            return Err("output_length_or_alphabet_changed".into());
        }
        if original
            .iter()
            .zip(sequence)
            .enumerate()
            .any(|(i, (a, b))| self.frozen[i] && a != b)
        {
            return Err("frozen_base_changed".into());
        }
        if translate(sequence, self.input.cds, &self.input.code)
            != translate(original, self.input.cds, &self.input.code)
        {
            return Err("translation_changed".into());
        }
        Ok(Violations {
            matches: self.matches(sequence)?,
            gc: self
                .gc
                .as_ref()
                .and_then(|gc| gc.violation(sequence, self.input.cds)),
            gc_window: self
                .gc_window
                .as_ref()
                .and_then(|gc| gc.violation(sequence)),
        })
    }

    fn initial_infeasibility(&self, matches: &[Match]) -> Option<&'static str> {
        if matches.iter().any(|hit| {
            self.frozen[hit.interval.start..hit.interval.end]
                .iter()
                .all(|b| *b)
        }) {
            Some("forbidden_match_entirely_in_frozen_bases")
        } else {
            self.gc
                .as_ref()
                .and_then(PreparedGc::infeasibility)
                .or_else(|| {
                    self.gc_window
                        .as_ref()
                        .and_then(|gc| gc.infeasibility(&self.frozen))
                })
        }
    }
}

/// Fresh full-sequence validation, suitable for approval without rerunning search.
pub fn validate_output(input: &Input, output: &str) -> Result<Vec<Edit>, String> {
    let prepared = prepare(input).map_err(|o| o.reason)?;
    let violations = prepared.validate(output.as_bytes())?;
    if !violations.satisfied() {
        return Err(if !violations.matches.is_empty() {
            "forbidden_motifs_remain"
        } else if violations.gc.is_some() {
            "gc_bounds_not_satisfied"
        } else {
            "gc_window_bounds_not_satisfied"
        }
        .into());
    }
    Ok(edits(input.sequence.as_bytes(), output.as_bytes()))
}

fn edits(original: &[u8], output: &[u8]) -> Vec<Edit> {
    original
        .iter()
        .zip(output)
        .enumerate()
        .filter_map(|(position, (before, after))| {
            (before != after).then_some(Edit {
                position,
                before: *before,
                after: *after,
            })
        })
        .collect()
}

/// Deterministic full-candidate enumeration with a bounded work budget.
///
/// Initially violating codons are visited first, then remaining codons in local
/// order. Original codons precede alternatives sorted by nucleotide edit count
/// and lexicographic DNA. Enumeration is complete only when all variants were
/// actually evaluated; an interrupted/budget-limited run proves no infeasibility.
/// `cancel` is polled before every candidate and once before publishing a result.
pub fn solve(input: &Input, mut cancel: impl FnMut(u64) -> bool) -> Outcome {
    let mut prepared = match prepare(input) {
        Ok(p) => p,
        Err(o) => return o,
    };
    let original = input.sequence.as_bytes();
    let initial = match prepared.matches(original) {
        Ok(hits) => hits,
        Err(reason) => return Outcome::empty(Status::Invalid, reason),
    };
    let mut outcome = Outcome::empty(
        Status::SearchExhausted,
        "candidate_or_motif_work_budget_exhausted",
    );
    outcome.initial_matches = initial.clone();
    outcome.effective_evaluation_budget = prepared.budget;
    outcome.search_space = prepared.space;
    if cancel(0) {
        return cancelled(outcome);
    }
    if let Some(reason) = prepared.initial_infeasibility(&initial) {
        outcome.status = Status::ProvenInfeasible;
        outcome.reason = reason.into();
        return outcome;
    }
    prepared.choices.sort_by_key(|(start, _)| {
        (
            !initial
                .iter()
                .any(|hit| hit.interval.start < start + 3 && hit.interval.end > *start),
            *start,
        )
    });
    let mut indices = vec![0usize; prepared.choices.len()];
    let mut best: Option<(usize, Vec<u8>)> = None;
    let mut complete = false;
    loop {
        if cancel(outcome.evaluated_candidates) {
            return cancelled(outcome);
        }
        if outcome.evaluated_candidates == prepared.budget {
            break;
        }
        let mut candidate = original.to_vec();
        for ((start, choices), index) in prepared.choices.iter().zip(&indices) {
            candidate[*start..*start + 3].copy_from_slice(&choices[*index]);
        }
        outcome.evaluated_candidates += 1;
        let hits = match prepared.validate(&candidate) {
            Ok(hits) => hits,
            Err(reason) => return Outcome::empty(Status::Invalid, reason),
        };
        if hits.satisfied() {
            let count = edits(original, &candidate).len();
            if best
                .as_ref()
                .is_none_or(|(old_count, old)| (count, &candidate) < (*old_count, old))
            {
                best = Some((count, candidate));
            }
            if count == 0 {
                complete = true;
                break;
            }
        }
        let mut position = 0;
        while position < indices.len() {
            indices[position] += 1;
            if indices[position] < prepared.choices[position].1.len() {
                break;
            }
            indices[position] = 0;
            position += 1;
        }
        if position == indices.len() {
            complete = true;
            break;
        }
    }
    if cancel(outcome.evaluated_candidates) {
        return cancelled(outcome);
    }
    outcome.optimization_complete = complete;
    if let Some((_, sequence)) = best {
        // Do not publish the search's cached evaluation as the final validator.
        if let Err(reason) = prepared.validate(&sequence).and_then(|h| {
            if h.satisfied() {
                Ok(h)
            } else {
                Err("final_motif_validation_failed".into())
            }
        }) {
            return Outcome::empty(Status::Invalid, reason);
        }
        outcome.edits = edits(original, &sequence);
        if cancel(outcome.evaluated_candidates) {
            return cancelled(outcome);
        }
        outcome.sequence = Some(String::from_utf8(sequence).unwrap());
        outcome.status = Status::Feasible;
        outcome.minimum_edits_proven = complete;
        outcome.reason = if complete {
            "complete_enumeration_or_unchanged_feasible_input"
        } else {
            "validated_best_candidate_with_incomplete_optimization"
        }
        .into();
    } else if complete {
        outcome.status = Status::ProvenInfeasible;
        outcome.reason = "complete_enumeration_found_no_feasible_variant".into();
    }
    outcome
}

/// Choose a versioned bounded solver without changing admission or validation.
pub fn solve_with_strategy(
    input: &Input,
    strategy: SearchStrategy,
    cancel: impl FnMut(u64) -> bool,
) -> Outcome {
    match strategy {
        SearchStrategy::FullEnumeration => solve(input, cancel),
        SearchStrategy::ConflictDirected => conflict_search::solve(input, cancel),
    }
}

mod conflict_search;

#[cfg(test)]
mod tests;
