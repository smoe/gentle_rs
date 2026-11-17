//! Unwired mismatch-rule primitive for a later opt-in supplied-primer profile.
//!
//! Inputs describe complete, unambiguous, ungapped alignments in each oligo's
//! own 5'-to-3' orientation. This module neither establishes that evidence nor
//! searches, pairs sites, identifies intended targets or declares specificity.

const TERMINAL_WINDOW_BP: usize = 5;
const MIN_TOTAL_MISMATCHES: usize = 2;
const MIN_TERMINAL_MISMATCHES: usize = 2;
const IGNORE_TOTAL_MISMATCHES: usize = 6;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct CompleteUngappedMismatchCounts {
    total: usize,
    terminal: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum InvalidMismatchCounts {
    EmptyPrimer,
    InconsistentCounts,
}

impl CompleteUngappedMismatchCounts {
    /// Validate count geometry, not the provenance/completeness of an alignment.
    pub(super) fn new(
        primer_length: usize,
        total: usize,
        terminal: usize,
    ) -> Result<Self, InvalidMismatchCounts> {
        if primer_length == 0 {
            return Err(InvalidMismatchCounts::EmptyPrimer);
        }
        let window = primer_length.min(TERMINAL_WINDOW_BP);
        if total > primer_length
            || terminal > window
            || terminal > total
            || total - terminal > primer_length - window
        {
            return Err(InvalidMismatchCounts::InconsistentCounts);
        }
        Ok(Self { total, terminal })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SiteMismatchClassification {
    pub(super) ignored_by_total_threshold: bool,
    pub(super) excluded_by_stringency: bool,
}

impl SiteMismatchClassification {
    pub(super) fn relevant(self) -> bool {
        !self.ignored_by_total_threshold && !self.excluded_by_stringency
    }
}

pub(super) fn classify_site(counts: CompleteUngappedMismatchCounts) -> SiteMismatchClassification {
    SiteMismatchClassification {
        ignored_by_total_threshold: counts.total >= IGNORE_TOTAL_MISMATCHES,
        excluded_by_stringency: counts.total >= MIN_TOTAL_MISMATCHES
            && counts.terminal >= MIN_TERMINAL_MISMATCHES,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ProductMismatchDisposition {
    RelevantOffTarget,
    ExcludedByStringency,
    IgnoredByTotalThreshold,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct ProductMismatchClassification {
    pub(super) first: SiteMismatchClassification,
    pub(super) second: SiteMismatchClassification,
    pub(super) disposition: ProductMismatchDisposition,
}

/// Classify an already paired, unintended product; preserve both sites' reasons.
pub(super) fn classify_product(
    first: CompleteUngappedMismatchCounts,
    second: CompleteUngappedMismatchCounts,
) -> ProductMismatchClassification {
    let first = classify_site(first);
    let second = classify_site(second);
    let disposition = if first.ignored_by_total_threshold || second.ignored_by_total_threshold {
        ProductMismatchDisposition::IgnoredByTotalThreshold
    } else if first.excluded_by_stringency || second.excluded_by_stringency {
        ProductMismatchDisposition::ExcludedByStringency
    } else {
        ProductMismatchDisposition::RelevantOffTarget
    };
    ProductMismatchClassification {
        first,
        second,
        disposition,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(total: usize, terminal: usize) -> CompleteUngappedMismatchCounts {
        CompleteUngappedMismatchCounts::new(20, total, terminal).unwrap()
    }

    #[test]
    fn primer_blast_mismatch_product_truth_table() {
        use ProductMismatchDisposition::*;
        for (first, second, expected) in [
            ((0, 0), (0, 0), RelevantOffTarget),
            ((1, 1), (0, 0), RelevantOffTarget),
            ((2, 0), (0, 0), RelevantOffTarget),
            ((1, 1), (1, 1), RelevantOffTarget),
            ((2, 1), (0, 0), RelevantOffTarget),
            ((2, 2), (0, 0), ExcludedByStringency),
            ((5, 1), (0, 0), RelevantOffTarget),
            ((5, 2), (0, 0), ExcludedByStringency),
            ((6, 0), (0, 0), IgnoredByTotalThreshold),
            ((0, 0), (6, 5), IgnoredByTotalThreshold),
            ((2, 2), (6, 0), IgnoredByTotalThreshold),
        ] {
            let first = counts(first.0, first.1);
            let second = counts(second.0, second.1);
            assert_eq!(classify_product(first, second).disposition, expected);
            assert_eq!(classify_product(second, first).disposition, expected);
        }
    }

    #[test]
    fn primer_blast_mismatch_counts_reject_impossible_geometry() {
        assert_eq!(
            CompleteUngappedMismatchCounts::new(0, 0, 0),
            Err(InvalidMismatchCounts::EmptyPrimer)
        );
        for (length, total, terminal) in
            [(20, 21, 0), (20, 1, 2), (20, 6, 6), (20, 20, 4), (5, 2, 1)]
        {
            assert_eq!(
                CompleteUngappedMismatchCounts::new(length, total, terminal),
                Err(InvalidMismatchCounts::InconsistentCounts)
            );
        }
        assert!(CompleteUngappedMismatchCounts::new(20, 20, 5).is_ok());
        assert!(CompleteUngappedMismatchCounts::new(3, 2, 2).is_ok());
        assert!(CompleteUngappedMismatchCounts::new(usize::MAX, usize::MAX, 5).is_ok());
    }

    #[test]
    fn primer_blast_mismatch_reasons_remain_separate() {
        let result = classify_product(counts(6, 2), counts(2, 2));
        assert!(result.first.ignored_by_total_threshold);
        assert!(result.first.excluded_by_stringency);
        assert!(!result.second.ignored_by_total_threshold);
        assert!(result.second.excluded_by_stringency);
        assert!(!result.first.relevant());
    }

    #[test]
    fn primer_blast_mismatch_product_decomposes_over_every_valid_twentymer_count() {
        let valid: Vec<_> = (0..=20)
            .flat_map(|total| {
                (0..=5).filter_map(move |terminal| {
                    CompleteUngappedMismatchCounts::new(20, total, terminal).ok()
                })
            })
            .collect();
        assert_eq!(valid.len(), 96);
        for first in &valid {
            for second in &valid {
                let result = classify_product(*first, *second);
                let ignored = first.total >= 6 || second.total >= 6;
                let excluded = (first.total >= 2 && first.terminal >= 2)
                    || (second.total >= 2 && second.terminal >= 2);
                let relevant = !ignored && !excluded;
                assert_eq!(
                    result.disposition == ProductMismatchDisposition::RelevantOffTarget,
                    relevant
                );
                assert_eq!(
                    result.first.relevant() && result.second.relevant(),
                    relevant
                );
                let swapped = classify_product(*second, *first);
                assert_eq!(swapped.disposition, result.disposition);
                assert_eq!(swapped.first, result.second);
                assert_eq!(swapped.second, result.first);
            }
        }
    }
}
