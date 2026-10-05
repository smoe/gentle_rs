//! Bounded conflict-directed codon assignment with full-candidate evaluation.
//!
//! A feasible sequence must change an unassigned codon overlapping the first
//! current violation. Branch over all such synonymous assignments, fixing one
//! per level. Complete traversal therefore covers a feasible/minimum-edit path
//! if one exists; budget interruption does not establish either conclusion.
//! The explicit depth-first stack avoids recursion and a sequence-sized frontier.

use super::*;

struct Frame {
    assigned_choice: Option<usize>,
    branches: Vec<(usize, usize)>,
    next: usize,
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
    if initial.iter().any(|hit| {
        prepared.frozen[hit.interval.start..hit.interval.end]
            .iter()
            .all(|base| *base)
    }) {
        outcome.status = Status::ProvenInfeasible;
        outcome.reason = "forbidden_match_entirely_in_frozen_bases".into();
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
        if hits.is_empty()
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

        let mut branches = vec![];
        // Assigned codons cannot revert in a descendant, so edit cost increases.
        // Prune only infeasible nodes already as costly as a validated best.
        if best
            .as_ref()
            .is_none_or(|(best_count, _)| count < *best_count)
        {
            if let Some(hit) = hits.first() {
                for (index, (start, alternatives)) in prepared.choices.iter().enumerate() {
                    if !assigned[index]
                        && *start < hit.interval.end
                        && start + 3 > hit.interval.start
                    {
                        // Unassigned codons remain original; alternative zero is original.
                        branches.extend((1..alternatives.len()).map(|choice| (index, choice)));
                    }
                }
            }
        }
        frames.push(Frame {
            assigned_choice,
            branches,
            next: 0,
        });
        loop {
            let Some(frame) = frames.last_mut() else {
                complete = true;
                break 'search;
            };
            if let Some(&(index, choice)) = frame.branches.get(frame.next) {
                frame.next += 1;
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
            hits.is_empty()
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
