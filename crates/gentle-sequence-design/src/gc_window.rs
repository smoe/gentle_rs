//! Exact GC bounds on every full, one-base-step window inside the declared CDS.

use super::*;

/// No partial edge windows, flanks, inferred defaults or configurable stride.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WindowGcBounds {
    pub window_bp: usize,
    pub min_basis_points: u16,
    pub max_basis_points: u16,
}

impl WindowGcBounds {
    pub(super) fn validate(self, cds_length: usize) -> Result<(), &'static str> {
        GcBounds {
            min_basis_points: self.min_basis_points,
            max_basis_points: self.max_basis_points,
        }
        .validate()?;
        if self.window_bp == 0 || self.window_bp > cds_length {
            Err("gc_window_length_requires_1_to_declared_cds_length")
        } else {
            Ok(())
        }
    }
}

/// Exact sequence-local window and count, including any frozen bases it contains.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WindowGcMeasurement {
    pub interval: Interval,
    pub gc_bases: usize,
    pub satisfies_bounds: bool,
}

/// All windows are retained in ascending sequence-local order; no sampled summary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WindowGcAssessment {
    pub window_bp: usize,
    pub window_count: usize,
    pub minimum_gc_bases: usize,
    pub maximum_gc_bases: usize,
    pub input: Vec<WindowGcMeasurement>,
    pub output: Option<Vec<WindowGcMeasurement>>,
}

pub(super) struct PreparedWindowGc {
    pub assessment: WindowGcAssessment,
    cds: Interval,
}

impl PreparedWindowGc {
    pub fn new(input: &Input) -> Self {
        // Admission already validated length/bounds before any window allocation.
        let bounds = input.gc_window.unwrap();
        let mut prepared = Self {
            assessment: WindowGcAssessment {
                window_bp: bounds.window_bp,
                window_count: input.cds.end - input.cds.start - bounds.window_bp + 1,
                minimum_gc_bases: (usize::from(bounds.min_basis_points) * bounds.window_bp)
                    .div_ceil(10_000),
                maximum_gc_bases: usize::from(bounds.max_basis_points) * bounds.window_bp / 10_000,
                input: vec![],
                output: None,
            },
            cds: input.cds,
        };
        prepared.assessment.input = prepared.measurements(input.sequence.as_bytes());
        prepared
    }

    fn counts(&self, sequence: &[u8]) -> impl Iterator<Item = (Interval, usize)> {
        // Prefix counts bound each complete evaluation to O(CDS length + windows),
        // even when windows overlap almost entirely. This is not a local evaluator.
        let mut prefix = Vec::with_capacity(self.cds.end - self.cds.start + 1);
        prefix.push(0usize);
        for base in &sequence[self.cds.start..self.cds.end] {
            prefix.push(prefix.last().unwrap() + usize::from(matches!(base, b'G' | b'C')));
        }
        let window_bp = self.assessment.window_bp;
        let start = self.cds.start;
        (0..self.assessment.window_count).map(move |offset| {
            (
                Interval {
                    start: start + offset,
                    end: start + offset + window_bp,
                },
                prefix[offset + window_bp] - prefix[offset],
            )
        })
    }

    pub fn measurements(&self, sequence: &[u8]) -> Vec<WindowGcMeasurement> {
        self.counts(sequence)
            .map(|(interval, gc_bases)| WindowGcMeasurement {
                interval,
                gc_bases,
                satisfies_bounds: self.direction(gc_bases).is_none(),
            })
            .collect()
    }

    fn direction(&self, count: usize) -> Option<GcViolation> {
        if count < self.assessment.minimum_gc_bases {
            Some(GcViolation::BelowMinimum)
        } else if count > self.assessment.maximum_gc_bases {
            Some(GcViolation::AboveMaximum)
        } else {
            None
        }
    }

    pub fn violation(&self, sequence: &[u8]) -> Option<(Interval, GcViolation)> {
        self.counts(sequence).find_map(|(interval, count)| {
            self.direction(count).map(|direction| (interval, direction))
        })
    }

    pub fn infeasibility(&self, frozen: &[bool]) -> Option<&'static str> {
        if self.assessment.minimum_gc_bases > self.assessment.maximum_gc_bases {
            Some("gc_window_integer_count_bounds_empty")
        } else if self.assessment.input.iter().any(|window| {
            !window.satisfies_bounds
                && frozen[window.interval.start..window.interval.end]
                    .iter()
                    .all(|b| *b)
        }) {
            Some("gc_window_violation_entirely_in_frozen_bases")
        } else {
            None
        }
    }
}

pub(super) fn direction_capable(
    start: usize,
    original: &[u8],
    alternatives: &[[u8; 3]],
    window: Interval,
    direction: GcViolation,
) -> bool {
    let left = window.start.saturating_sub(start).min(3);
    let right = window.end.saturating_sub(start).min(3);
    if left >= right {
        return false;
    }
    let original_count = gc::count(&original[left..right]);
    alternatives.iter().any(|codon| match direction {
        GcViolation::BelowMinimum => gc::count(&codon[left..right]) > original_count,
        GcViolation::AboveMaximum => gc::count(&codon[left..right]) < original_count,
    })
}

/// Absent policy or output means absent facts, never zero-valued observations.
pub fn assess_window_gc(
    input: &Input,
    output: Option<&str>,
) -> Result<Option<WindowGcAssessment>, String> {
    let prepared = prepare(input).map_err(|o| o.reason)?;
    let Some(gc) = &prepared.gc_window else {
        return Ok(None);
    };
    let mut assessment = gc.assessment.clone();
    if let Some(output) = output {
        if !prepared.validate(output.as_bytes())?.satisfied() {
            return Err("output_constraints_not_satisfied".into());
        }
        assessment.output = Some(gc.measurements(output.as_bytes()));
    }
    Ok(Some(assessment))
}
