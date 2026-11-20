//! Bounded exact canonical-DNA reference recognition, shared by annotation consumers.
//!
//! Ambiguous query bases never match a canonical reference. Coordinates remain
//! on the supplied forward axis and a reverse hit refers to the reference's
//! reverse complement. This helper makes no biological interpretation.

#[derive(Debug, PartialEq, Eq)]
pub(super) struct ExactReferenceHit {
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
}

#[derive(Debug, Default, PartialEq, Eq)]
pub(super) struct ExactReferenceScan {
    pub hits: Vec<ExactReferenceHit>,
    pub scanned_windows: usize,
    pub ambiguous_windows: usize,
}

fn charge(work: &mut usize, max_work: usize) -> Result<(), String> {
    *work = work.checked_add(1).ok_or("Exact reference work overflow")?;
    if *work > max_work {
        return Err("budget_exceeded: exact reference comparison/window budget exceeded".into());
    }
    Ok(())
}

/// Enumerate exact matches with a shared, deterministic comparison budget.
/// Exceeding the budget is an error, never a successful negative search.
pub(super) fn exact_reference_matches(
    sequence: &[u8],
    reference: &[u8],
    work: &mut usize,
    max_work: usize,
    max_hits: usize,
    strand: Option<bool>,
) -> Result<ExactReferenceScan, String> {
    if reference.is_empty()
        || !reference
            .iter()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
    {
        return Err("Exact reference must contain nonempty canonical uppercase DNA".into());
    }
    let reverse = reference
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            _ => b'C',
        })
        .collect::<Vec<_>>();
    let mut scan = ExactReferenceScan::default();
    if sequence.len() < reference.len() {
        return Ok(scan);
    }
    for (is_reverse, pattern) in [(false, reference), (true, reverse.as_slice())] {
        if strand.is_some_and(|v| v != is_reverse) {
            continue;
        }
        let ambiguous_base = |b: u8| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T');
        let mut ambiguous = sequence[..reference.len()]
            .iter()
            .filter(|b| ambiguous_base(**b))
            .count();
        for start in 0..=sequence.len() - reference.len() {
            charge(work, max_work)?;
            scan.scanned_windows += 1;
            if start > 0 {
                ambiguous -= usize::from(ambiguous_base(sequence[start - 1]));
                ambiguous += usize::from(ambiguous_base(sequence[start + reference.len() - 1]));
            }
            if ambiguous > 0 {
                scan.ambiguous_windows += 1;
                continue;
            }
            let mut matched = true;
            for (a, b) in sequence[start..start + reference.len()].iter().zip(pattern) {
                charge(work, max_work)?;
                if a.to_ascii_uppercase() != *b {
                    matched = false;
                    break;
                }
            }
            if matched {
                scan.hits.push(ExactReferenceHit {
                    start,
                    end: start + reference.len(),
                    reverse: is_reverse,
                });
                if scan.hits.len() >= max_hits {
                    return Ok(scan);
                }
            }
        }
    }
    Ok(scan)
}

#[cfg(test)]
mod tests {
    use super::*;
    // Synthetic strings authored here, recreated by the tests; no biological references.
    #[test]
    fn exact_references_preserve_both_strands_and_ambiguity() {
        let hits =
            exact_reference_matches(b"NACGTTAACGN", b"AACG", &mut 0, 1000, 100, None).unwrap();
        assert_eq!(
            hits.hits,
            vec![
                ExactReferenceHit {
                    start: 6,
                    end: 10,
                    reverse: false
                },
                ExactReferenceHit {
                    start: 2,
                    end: 6,
                    reverse: true
                }
            ]
        );
        assert!(
            exact_reference_matches(b"ANCG", b"AACG", &mut 0, 1000, 100, None)
                .unwrap()
                .hits
                .is_empty()
        );
        assert!(exact_reference_matches(b"AAA", b"AN", &mut 0, 1000, 100, None).is_err());
    }
    #[test]
    fn budget_refusal_is_not_a_negative_search() {
        assert!(exact_reference_matches(b"AAAAAAAAA", b"AAAA", &mut 0, 1, 100, None).is_err());
    }
}
