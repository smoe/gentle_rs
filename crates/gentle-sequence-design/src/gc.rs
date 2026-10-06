//! Exact CDS-wide GC count bounds and conservative synonymous reachability.

use super::*;

/// Inclusive hundredths-of-a-percent bounds over the complete declared CDS.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GcBounds {
    pub min_basis_points: u16,
    pub max_basis_points: u16,
}

impl GcBounds {
    pub(super) fn validate(self) -> Result<(), &'static str> {
        if self.min_basis_points > self.max_basis_points || self.max_basis_points > 10_000 {
            Err("gc_bounds_require_0_to_10000_inclusive_and_minimum_not_above_maximum")
        } else {
            Ok(())
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GcMeasurement {
    pub gc_bases: usize,
    pub satisfies_bounds: bool,
}

/// Reachability extrema are a relaxation: intersection does not prove feasibility.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GcAssessment {
    pub denominator_bases: usize,
    pub minimum_gc_bases: usize,
    pub maximum_gc_bases: usize,
    pub reachable_minimum_gc_bases: usize,
    pub reachable_maximum_gc_bases: usize,
    pub input: GcMeasurement,
    pub output: Option<GcMeasurement>,
}

#[derive(Clone, Copy)]
pub(super) enum GcViolation {
    BelowMinimum,
    AboveMaximum,
}

pub(super) struct PreparedGc {
    pub assessment: GcAssessment,
}

pub(super) fn count(sequence: &[u8]) -> usize {
    sequence.iter().filter(|b| matches!(b, b'G' | b'C')).count()
}

impl PreparedGc {
    pub fn new(input: &Input, choices: &[(usize, Vec<[u8; 3]>)]) -> Result<Self, String> {
        let bounds = input.gc_content.unwrap();
        bounds.validate().map_err(String::from)?;
        let length = input.cds.end - input.cds.start;
        let minimum = (usize::from(bounds.min_basis_points) * length).div_ceil(10_000);
        let maximum = usize::from(bounds.max_basis_points) * length / 10_000;
        let original = &input.sequence.as_bytes()[input.cds.start..input.cds.end];
        let initial = count(original);
        let mut reachable_minimum = initial;
        let mut reachable_maximum = initial;
        for (start, alternatives) in choices {
            let original_count = count(&input.sequence.as_bytes()[*start..*start + 3]);
            let counts = alternatives.iter().map(|codon| count(codon));
            reachable_minimum = reachable_minimum - original_count + counts.clone().min().unwrap();
            reachable_maximum = reachable_maximum - original_count + counts.max().unwrap();
        }
        Ok(Self {
            assessment: GcAssessment {
                denominator_bases: length,
                minimum_gc_bases: minimum,
                maximum_gc_bases: maximum,
                reachable_minimum_gc_bases: reachable_minimum,
                reachable_maximum_gc_bases: reachable_maximum,
                input: GcMeasurement {
                    gc_bases: initial,
                    satisfies_bounds: minimum <= initial && initial <= maximum,
                },
                output: None,
            },
        })
    }

    pub fn measurement(&self, sequence: &[u8], cds: Interval) -> GcMeasurement {
        let gc_bases = count(&sequence[cds.start..cds.end]);
        GcMeasurement {
            gc_bases,
            satisfies_bounds: self.assessment.minimum_gc_bases <= gc_bases
                && gc_bases <= self.assessment.maximum_gc_bases,
        }
    }

    pub fn violation(&self, sequence: &[u8], cds: Interval) -> Option<GcViolation> {
        let count = self.measurement(sequence, cds).gc_bases;
        if count < self.assessment.minimum_gc_bases {
            Some(GcViolation::BelowMinimum)
        } else if count > self.assessment.maximum_gc_bases {
            Some(GcViolation::AboveMaximum)
        } else {
            None
        }
    }

    pub fn infeasibility(&self) -> Option<&'static str> {
        let a = &self.assessment;
        if a.minimum_gc_bases > a.maximum_gc_bases {
            Some("gc_count_window_empty_for_declared_cds_length")
        } else if a.reachable_maximum_gc_bases < a.minimum_gc_bases
            || a.reachable_minimum_gc_bases > a.maximum_gc_bases
        {
            Some("gc_bounds_unreachable_with_frozen_bases_and_synonymous_choices")
        } else {
            None
        }
    }
}

pub(super) fn direction_capable(
    original: &[u8],
    alternatives: &[[u8; 3]],
    violation: GcViolation,
) -> bool {
    let original_count = count(original);
    alternatives.iter().any(|codon| match violation {
        GcViolation::BelowMinimum => count(codon) > original_count,
        GcViolation::AboveMaximum => count(codon) < original_count,
    })
}

/// Derive portable GC facts from admitted input and, if present, validated output.
/// No requested constraint means no assessment; no output means no output facts.
pub fn assess_gc(input: &Input, output: Option<&str>) -> Result<Option<GcAssessment>, String> {
    let prepared = prepare(input).map_err(|o| o.reason)?;
    let Some(gc) = &prepared.gc else {
        return Ok(None);
    };
    let mut assessment = gc.assessment.clone();
    if let Some(output) = output {
        if !prepared.validate(output.as_bytes())?.satisfied() {
            return Err("output_constraints_not_satisfied".into());
        }
        assessment.output = Some(gc.measurement(output.as_bytes(), input.cds));
    }
    Ok(Some(assessment))
}
