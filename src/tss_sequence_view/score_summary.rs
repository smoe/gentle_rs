//! Deterministic within-matrix summaries and bounded SVG hover annotation.

use super::TssViewTrace;
use serde::Serialize;
use std::{collections::BTreeSet, ops::Range};

/// Raw-score positions within one matrix lane, never a ranking across factors.
pub const TSS_SCORE_SUMMARY_LIMIT: usize = 3;
/// Limits SVG annotations only; curves and validity bands retain every position.
pub const TSS_SVG_HOVER_LIMIT_PER_STRAND: usize = 128;

/// One evaluated motif-window start on the displayed, transcript-oriented DNA.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct TssScoreWindow {
    pub start_0based: usize,
    pub reverse: bool,
    pub raw_score: f64,
}

impl TssViewTrace {
    /// Highest finite, evaluated scores among visible window starts. Ties use
    /// the local start then forward before reverse; clipping never changes rank.
    pub fn strongest_windows(&self, span: Range<usize>) -> Vec<TssScoreWindow> {
        self.rank_windows(span, None)
    }

    fn rank_windows(&self, span: Range<usize>, strand: Option<bool>) -> Vec<TssScoreWindow> {
        let mut ranked = Vec::with_capacity(TSS_SCORE_SUMMARY_LIMIT + 1);
        for (reverse, scores) in [(false, &self.forward), (true, &self.reverse)] {
            if strand.is_some_and(|value| value != reverse) {
                continue;
            }
            let start = span
                .start
                .saturating_sub(self.start_0based)
                .min(scores.len());
            let end = span.end.saturating_sub(self.start_0based).min(scores.len());
            for (index, score) in scores[start..end.max(start)].iter().enumerate() {
                let Some(raw_score) = score.filter(|value| value.is_finite()) else {
                    continue;
                };
                ranked.push(TssScoreWindow {
                    start_0based: self.start_0based + start + index,
                    reverse,
                    raw_score,
                });
                ranked.sort_by(|a, b| {
                    b.raw_score
                        .total_cmp(&a.raw_score)
                        .then_with(|| a.start_0based.cmp(&b.start_0based))
                        .then_with(|| a.reverse.cmp(&b.reverse))
                });
                ranked.truncate(TSS_SCORE_SUMMARY_LIMIT);
            }
        }
        ranked
    }

    /// Retain strongest starts plus evenly spaced valid starts, filling duplicate
    /// sample slots in local order. Missing values never consume a hover slot.
    pub(crate) fn hover_starts(
        &self,
        span: Range<usize>,
        reverse: bool,
    ) -> (BTreeSet<usize>, usize) {
        let scores = if reverse {
            &self.reverse
        } else {
            &self.forward
        };
        let positions: Vec<_> = scores
            .iter()
            .enumerate()
            .filter_map(|(i, value)| {
                let start = self.start_0based + i;
                (span.contains(&start) && value.is_some_and(|value| value.is_finite()))
                    .then_some(start)
            })
            .collect();
        let total = positions.len();
        if total <= TSS_SVG_HOVER_LIMIT_PER_STRAND {
            return (positions.into_iter().collect(), total);
        }
        let mut selected: BTreeSet<_> = self
            .rank_windows(span, Some(reverse))
            .into_iter()
            .map(|window| window.start_0based)
            .collect();
        let slots = TSS_SVG_HOVER_LIMIT_PER_STRAND - selected.len();
        for i in 0..slots {
            selected.insert(positions[i * (total - 1) / (slots - 1)]);
        }
        for position in positions {
            if selected.len() == TSS_SVG_HOVER_LIMIT_PER_STRAND {
                break;
            }
            selected.insert(position);
        }
        (selected, total)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strongest_windows_preserve_raw_scores_zero_gaps_scope_and_tie_order() {
        let trace = TssViewTrace {
            start_0based: 100,
            end_0based_exclusive: 107,
            motif_length_bp: 3,
            forward: vec![Some(-1.0), None, Some(0.0), Some(3.0), Some(3.0)],
            reverse: vec![Some(99.0), Some(f64::NAN), None, Some(3.0), Some(-2.0)],
            clip_negative: true,
            range_is_fallback: false,
        };
        let rows = trace.strongest_windows(101..105);
        assert_eq!(
            rows.iter()
                .map(|w| (w.start_0based, w.reverse, w.raw_score))
                .collect::<Vec<_>>(),
            [(103, false, 3.0), (103, true, 3.0), (104, false, 3.0)]
        );
        assert_eq!(trace.strongest_windows(102..103)[0].raw_score, 0.0);
        assert_eq!(trace.strongest_windows(104..105)[1].raw_score, -2.0);
        assert!(trace.strongest_windows(0..50).is_empty());
        assert!(trace.strongest_windows(106..110).is_empty());
    }

    #[test]
    fn hover_sampling_retains_strongest_starts_and_has_a_fixed_bound() {
        for length in [130, 700, 10_000] {
            let mut trace = TssViewTrace {
                start_0based: 0,
                end_0based_exclusive: length + 2,
                motif_length_bp: 3,
                forward: vec![Some(0.0); length],
                reverse: vec![None; length],
                clip_negative: true,
                range_is_fallback: false,
            };
            trace.forward[length / 2] = Some(10.0);
            let (starts, total) = trace.hover_starts(0..length, false);
            assert_eq!(total, length);
            assert_eq!(starts.len(), TSS_SVG_HOVER_LIMIT_PER_STRAND);
            assert!(starts.contains(&(length / 2)));
            assert_eq!(
                (starts.clone(), total),
                trace.hover_starts(0..length, false)
            );
            assert_eq!(trace.hover_starts(0..length, true), (BTreeSet::new(), 0));
        }
    }
}
