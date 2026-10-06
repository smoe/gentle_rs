//! Bounded conflict-directed codon assignment with full-candidate evaluation.
//!
//! A feasible sequence must change an unassigned codon implicated by a current
//! motif or whole-CDS GC violation. Branch over its synonymous assignments, fixing one
//! per level. Complete traversal therefore covers a feasible/minimum-edit path
//! if one exists; budget interruption does not establish either conclusion.
//! The explicit depth-first stack avoids recursion and a sequence-sized frontier.

use super::*;

struct Frame {
    assigned_choice: Option<usize>,
    conflict: Option<Conflict>,
    next_codon: usize,
    next_synonym: usize,
}

impl Frame {
    fn next_branch(
        &mut self,
        prepared: &Prepared<'_>,
        assigned: &[bool],
    ) -> Option<(usize, usize)> {
        let conflict = self.conflict?;
        while let Some((start, alternatives)) = prepared.choices.get(self.next_codon) {
            if !assigned[self.next_codon] {
                let implicated = match conflict {
                    Conflict::Motif(interval) => {
                        *start < interval.end && start + 3 > interval.start
                    }
                    // Capability filters codons, never individual assignments:
                    // neutral/worsening synonyms remain branches.
                    Conflict::Gc(violation) => gc::direction_capable(
                        &prepared.input.sequence.as_bytes()[*start..*start + 3],
                        alternatives,
                        violation,
                    ),
                };
                if implicated && self.next_synonym < alternatives.len() {
                    let branch = (self.next_codon, self.next_synonym);
                    self.next_synonym += 1;
                    return Some(branch);
                }
            }
            self.next_codon += 1;
            self.next_synonym = 1;
        }
        None
    }
}

pub(super) fn solve(input: &Input, mut cancel: impl FnMut(u64) -> bool) -> Outcome {
    let prepared = match prepare(input) {
        Ok(prepared) => prepared,
        Err(outcome) => return outcome,
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

    let mut candidate = original.to_vec();
    let mut assigned = vec![false; prepared.choices.len()];
    let mut frames: Vec<Frame> = vec![];
    let mut assigned_choice = None;
    let mut best: Option<(usize, Vec<u8>)> = None;
    let mut complete = false;
    'search: loop {
        if cancel(outcome.evaluated_candidates) {
            return cancelled(outcome);
        }
        if outcome.evaluated_candidates == prepared.budget {
            break;
        }
        outcome.evaluated_candidates += 1;
        let hits = match prepared.validate(&candidate) {
            Ok(hits) => hits,
            Err(reason) => return Outcome::empty(Status::Invalid, reason),
        };
        let count = original
            .iter()
            .zip(&candidate)
            .filter(|(before, after)| before != after)
            .count();
        if hits.satisfied()
            && best
                .as_ref()
                .is_none_or(|(old_count, old)| (count, &candidate) < (*old_count, old))
        {
            best = Some((count, candidate.clone()));
            if count == 0 {
                complete = true;
                break;
            }
        }

        // Assigned codons cannot revert in a descendant, so edit cost increases.
        // Prune only infeasible nodes already as costly as a validated best.
        let conflict = if best
            .as_ref()
            .is_none_or(|(best_count, _)| count < *best_count)
        {
            hits.first_conflict()
        } else {
            None
        };
        frames.push(Frame {
            assigned_choice,
            conflict,
            next_codon: 0,
            next_synonym: 1,
        });
        loop {
            let Some(frame) = frames.last_mut() else {
                complete = true;
                break 'search;
            };
            if let Some((index, choice)) = frame.next_branch(&prepared, &assigned) {
                let (start, alternatives) = &prepared.choices[index];
                candidate[*start..*start + 3].copy_from_slice(&alternatives[choice]);
                assigned[index] = true;
                assigned_choice = Some(index);
                break;
            }
            if let Some(index) = frames.pop().unwrap().assigned_choice {
                let start = prepared.choices[index].0;
                candidate[start..start + 3].copy_from_slice(&original[start..start + 3]);
                assigned[index] = false;
            }
        }
    }
    if cancel(outcome.evaluated_candidates) {
        return cancelled(outcome);
    }
    outcome.optimization_complete = complete;
    if let Some((_, sequence)) = best {
        if let Err(reason) = prepared.validate(&sequence).and_then(|hits| {
            hits.satisfied()
                .then_some(())
                .ok_or_else(|| "final_motif_validation_failed".into())
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
            "complete_conflict_search_or_unchanged_feasible_input"
        } else {
            "validated_best_candidate_with_incomplete_optimization"
        }
        .into();
    } else if complete {
        outcome.status = Status::ProvenInfeasible;
        outcome.reason = "complete_conflict_search_found_no_feasible_variant".into();
    }
    outcome
}
