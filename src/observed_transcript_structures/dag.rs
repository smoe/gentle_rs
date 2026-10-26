//! Best-equivalent query-interval DAG shared with Phase 2; paths are not node unions.

use crate::engine::EngineError;

pub(crate) struct Dag {
    pub weight: Vec<usize>,
    pub prefix: Vec<usize>,
    pub suffix: Vec<usize>,
    pub best: usize,
    pub optimal: Vec<usize>,
    intervals: Vec<(usize, usize)>,
}

pub(crate) fn charge(work: &mut usize, limit: usize) -> Result<(), EngineError> {
    charge_many(work, limit, 1)
}

pub(crate) fn charge_many(
    work: &mut usize,
    limit: usize,
    amount: usize,
) -> Result<(), EngineError> {
    *work = work
        .checked_add(amount)
        .filter(|n| *n <= limit)
        .ok_or_else(|| EngineError::invalid_input("primer_payload_chain_comparison_budget"))?;
    Ok(())
}

/// Callers supply query-sorted intervals and masked supported-base weights.
pub(crate) fn best_equivalent(
    intervals: Vec<(usize, usize)>,
    weight: Vec<usize>,
    work: &mut usize,
    limit: usize,
) -> Result<Dag, EngineError> {
    if intervals.len() != weight.len()
        || !intervals.windows(2).all(|w| w[0] <= w[1])
        || intervals
            .iter()
            .zip(&weight)
            .any(|(&(a, z), &n)| a >= z || n > z - a)
    {
        return Err(EngineError::invalid_input("observed_structure_dag_input"));
    }
    let n = intervals.len();
    let mut prefix = weight.clone();
    let mut suffix = weight.clone();
    for i in 0..n {
        for j in 0..i {
            charge(work, limit)?;
            if intervals[j].1 <= intervals[i].0 {
                prefix[i] = prefix[i].max(prefix[j] + weight[i]);
            }
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            charge(work, limit)?;
            if intervals[i].1 <= intervals[j].0 {
                suffix[i] = suffix[i].max(weight[i] + suffix[j]);
            }
        }
    }
    let best = prefix.iter().copied().max().unwrap_or(0);
    let optimal = (0..n)
        .filter(|i| weight[*i] > 0 && prefix[*i] + suffix[*i] - weight[*i] == best)
        .collect();
    Ok(Dag {
        weight,
        prefix,
        suffix,
        best,
        optimal,
        intervals,
    })
}

impl Dag {
    /// Enumerate only legal positive-weight best paths; resource exhaustion is an error.
    pub fn paths(
        &self,
        cap: usize,
        work: &mut usize,
        limit: usize,
    ) -> Result<Vec<Vec<usize>>, EngineError> {
        let mut pending = self
            .optimal
            .iter()
            .filter(|&&i| self.prefix[i] == self.weight[i])
            .map(|&i| vec![i])
            .collect::<Vec<_>>();
        let mut result = Vec::new();
        if pending.len() > cap {
            return Err(EngineError::invalid_input(
                "observed_structure_alternative_budget",
            ));
        }
        while let Some(path) = pending.pop() {
            charge(work, limit)?;
            let i = *path.last().expect("nonempty DAG path");
            if self.suffix[i] == self.weight[i] {
                result.push(path);
            } else {
                for &j in &self.optimal {
                    charge(work, limit)?;
                    if j > i
                        && self.intervals[i].1 <= self.intervals[j].0
                        && self.prefix[i] + self.suffix[j] == self.best
                    {
                        charge_many(work, limit, path.len())?;
                        let mut next = path.clone();
                        next.push(j);
                        pending.push(next);
                    }
                }
            }
            if result.len() + pending.len() > cap {
                return Err(EngineError::invalid_input(
                    "observed_structure_alternative_budget",
                ));
            }
        }
        result.sort();
        Ok(result)
    }
}
