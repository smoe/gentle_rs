//! Pure, bounded, uncalibrated termination candidate heuristics.
//! Input is a non-wrapping DNA window. Output uses window-local ascending source
//! coordinates with explicit strand; no external RNA folding program is called.
//! Bacterial candidates test one local stem/loop model, not co-transcriptional
//! folding, attenuation, Rho activity or measured termination. Mammalian candidates
//! require an authored hexamer/context rule; rare signals, chromatin, competition
//! and RNA-end usage are unmodelled. Cleavage windows enclose passing context
//! hypotheses and may contain intervening positions that did not individually
//! pass. Neither mode estimates a probability or defines a mature transcript.

use gentle_protocol::transcription_termination::*;
use std::collections::BTreeMap;

pub(super) struct PredictionScan {
    pub rows: Vec<TerminationEvidenceRow>,
    pub scanned_windows: usize,
    pub ambiguous_windows: usize,
    pub missing_context_windows: usize,
    pub work: usize,
    pub truncated: bool,
    pub warnings: Vec<String>,
}

fn canonical(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}
fn mapped_interval(start: usize, end: usize, len: usize, reverse: bool) -> (usize, usize) {
    if reverse {
        (len - end, len - start)
    } else {
        (start, end)
    }
}
fn row(
    sequence: &[u8],
    start: usize,
    end: usize,
    reverse: bool,
    kind: TerminationEvidenceKind,
    geometry: TerminationGeometryKind,
    label: &str,
    components: BTreeMap<String, f64>,
    related: Vec<TerminationInterval>,
) -> TerminationEvidenceRow {
    let (source_start, source_end) = mapped_interval(start, end, sequence.len(), reverse);
    TerminationEvidenceRow {
        row_id: String::new(),
        evidence_kind: kind,
        geometry_kind: geometry,
        label: label.into(),
        start_0based: source_start,
        end_0based_exclusive: source_end,
        reverse,
        element_id: None,
        source_feature_id: None,
        subject_feature_id: None,
        sequence_5prime_to_3prime: Some(String::from_utf8_lossy(&sequence[start..end]).into()),
        components,
        related_intervals: related
            .into_iter()
            .map(|mut interval| {
                (interval.start_0based, interval.end_0based_exclusive) = if reverse
                    && interval.coordinate_basis == TerminationCoordinateBasis::InterbaseBoundarySet
                {
                    (
                        sequence.len() - (interval.end_0based_exclusive - 1),
                        sequence.len() - interval.start_0based + 1,
                    )
                } else {
                    mapped_interval(
                        interval.start_0based,
                        interval.end_0based_exclusive,
                        sequence.len(),
                        reverse,
                    )
                };
                interval
            })
            .collect(),
    }
}
fn bounded_tick(scan: &mut PredictionScan, request: &TerminationScreenRequest) -> bool {
    if scan.work >= request.max_work || scan.rows.len() >= request.max_rows {
        scan.truncated = true;
        false
    } else {
        scan.work += 1;
        true
    }
}

/// Validate configurable authored policies before allocating or scanning.
pub(super) fn validate_prediction_policy(request: &TerminationScreenRequest) -> Result<(), String> {
    if request.max_rows == 0 || request.max_work == 0 {
        return Err("Prediction row/work limits must be positive".into());
    }
    if request.bacterial_intrinsic {
        let p = &request.bacterial_policy;
        if p.version != "bacterial_intrinsic_heuristic.v1"
            || p.min_stem_bp < 3
            || p.min_stem_bp > p.max_stem_bp
            || p.max_stem_bp > 32
            || p.min_loop_bp < 1
            || p.min_loop_bp > p.max_loop_bp
            || p.max_loop_bp > 64
            || p.max_mismatches >= p.min_stem_bp
            || p.min_gc_pairs > p.min_stem_bp
            || p.max_tract_gap_bp > 16
            || p.tract_window_bp < 3
            || p.tract_window_bp > 32
            || p.min_tract_t == 0
            || p.min_tract_t > p.tract_window_bp
        {
            return Err("Invalid bacterial intrinsic heuristic policy or version".into());
        }
    }
    if request.mammalian_pol_ii {
        let p = &request.mammalian_policy;
        if p.version != "mammalian_pol_ii_heuristic.v1"
            || p.signals.is_empty()
            || p.signals.len() > 32
            || p.signals
                .iter()
                .any(|s| s.len() != 6 || !s.bytes().all(canonical))
            || p.signals
                .iter()
                .collect::<std::collections::BTreeSet<_>>()
                .len()
                != p.signals.len()
            || p.min_cleavage_offset_bp == 0
            || p.min_cleavage_offset_bp > p.max_cleavage_offset_bp
            || p.max_cleavage_offset_bp > 100
            || p.downstream_window_bp < 4
            || p.downstream_window_bp > 100
            || p.min_downstream_gt == 0
            || p.min_downstream_gt > p.downstream_window_bp
        {
            return Err("Invalid mammalian Pol II heuristic policy or version".into());
        }
    }
    Ok(())
}

/// Score both strands with bounded work and no installation-dependent inputs.
/// Rows contain no stable identity yet: the engine binds them to source inputs.
pub(super) fn predict_termination(
    sequence: &[u8],
    request: &TerminationScreenRequest,
) -> Result<PredictionScan, String> {
    validate_prediction_policy(request)?;
    let forward: Vec<u8> = sequence.iter().map(u8::to_ascii_uppercase).collect();
    let reverse: Vec<u8> = forward.iter().rev().map(|b| complement(*b)).collect();
    let mut scan = PredictionScan {
        rows: vec![],
        scanned_windows: 0,
        ambiguous_windows: 0,
        missing_context_windows: 0,
        work: 0,
        truncated: false,
        warnings: vec![],
    };
    for (oriented, reverse_strand) in [(&forward, false), (&reverse, true)] {
        if request.bacterial_intrinsic {
            bacterial(oriented, reverse_strand, request, &mut scan);
        }
        if !scan.truncated && request.mammalian_pol_ii {
            mammalian(oriented, reverse_strand, request, &mut scan);
        }
        if scan.truncated {
            break;
        }
    }
    if scan.truncated {
        scan.warnings.push(
            "Prediction work or row limit reached; remaining windows are unassessed, not negative."
                .into(),
        );
    }
    Ok(scan)
}

fn bacterial(
    sequence: &[u8],
    reverse: bool,
    request: &TerminationScreenRequest,
    scan: &mut PredictionScan,
) {
    let p = &request.bacterial_policy;
    for start in 0..sequence.len() {
        let mut best: Option<(usize, usize, usize, usize, usize, usize, usize)> = None;
        for stem in p.min_stem_bp..=p.max_stem_bp {
            for loop_bp in p.min_loop_bp..=p.max_loop_bp {
                if !bounded_tick(scan, request) {
                    return;
                }
                let hairpin_end = start + 2 * stem + loop_bp;
                let Some(tract_end) = hairpin_end.checked_add(p.tract_window_bp) else {
                    return;
                };
                if tract_end > sequence.len() {
                    continue;
                }
                scan.scanned_windows += 1;
                if !sequence[start..hairpin_end].iter().copied().all(canonical) {
                    scan.ambiguous_windows += 1;
                    continue;
                }
                let mut mismatches = 0;
                let mut gc = 0;
                let mut wobble = 0;
                for j in 0..stem {
                    let a = sequence[start + j];
                    let b = sequence[hairpin_end - 1 - j];
                    if complement(a) == b {
                        if a == b'G' || a == b'C' {
                            gc += 1;
                        }
                    } else if p.allow_gu_wobble && matches!((a, b), (b'G', b'T') | (b'T', b'G')) {
                        wobble += 1;
                    } else {
                        mismatches += 1;
                    }
                }
                if mismatches > p.max_mismatches || gc < p.min_gc_pairs {
                    continue;
                }
                for gap in 0..=p.max_tract_gap_bp {
                    if !bounded_tick(scan, request) {
                        return;
                    }
                    let end = hairpin_end + gap + p.tract_window_bp;
                    if end > sequence.len() {
                        continue;
                    }
                    let tract = &sequence[hairpin_end + gap..end];
                    if !tract.iter().copied().all(canonical) {
                        scan.ambiguous_windows += 1;
                        continue;
                    }
                    let ts = tract.iter().filter(|b| **b == b'T').count();
                    if ts < p.min_tract_t {
                        continue;
                    }
                    // Stable preference: more exact pairs and T bases, then fewer mismatches;
                    // geometry enumeration supplies deterministic ties.
                    let quality = stem * 4 + gc * 2 + ts * 2 - mismatches * 3 - wobble;
                    if best.as_ref().is_none_or(|v| quality > v.0) {
                        best = Some((quality, end, stem, loop_bp, mismatches, gc, wobble));
                    }
                }
            }
        }
        if let Some((_, end, stem, loop_bp, mismatches, gc, wobble)) = best {
            let components = BTreeMap::from([
                ("stem_bp".into(), stem as f64),
                ("loop_bp".into(), loop_bp as f64),
                ("mismatches".into(), mismatches as f64),
                ("gc_pairs".into(), gc as f64),
                ("wobble_pairs".into(), wobble as f64),
                (
                    "tract_t_count".into(),
                    sequence[end - p.tract_window_bp..end]
                        .iter()
                        .filter(|b| **b == b'T')
                        .count() as f64,
                ),
            ]);
            scan.rows.push(row(
                sequence,
                start,
                end,
                reverse,
                TerminationEvidenceKind::BacterialIntrinsicCandidate,
                TerminationGeometryKind::IntrinsicTerminatorCandidate,
                "Bacterial intrinsic candidate (uncalibrated)",
                components,
                vec![],
            ));
        }
    }
}

fn mammalian(
    sequence: &[u8],
    reverse: bool,
    request: &TerminationScreenRequest,
    scan: &mut PredictionScan,
) {
    let p = &request.mammalian_policy;
    let mut insufficient = false;
    for start in 0..sequence.len().saturating_sub(5) {
        if !bounded_tick(scan, request) {
            return;
        }
        scan.scanned_windows += 1;
        if !sequence[start..start + 6].iter().copied().all(canonical) {
            scan.ambiguous_windows += 1;
            continue;
        }
        if !p
            .signals
            .iter()
            .any(|s| s.as_bytes() == &sequence[start..start + 6])
        {
            continue;
        }
        let signal_end = start + 6;
        let mut first = None;
        let mut last = 0;
        let mut maximum_gt = 0;
        let mut possible = 0;
        for offset in p.min_cleavage_offset_bp..=p.max_cleavage_offset_bp {
            if !bounded_tick(scan, request) {
                return;
            }
            let boundary = signal_end + offset;
            let Some(end) = boundary.checked_add(p.downstream_window_bp) else {
                return;
            };
            if end > sequence.len() {
                insufficient = true;
                scan.missing_context_windows += 1;
                continue;
            }
            let context = &sequence[boundary..end];
            if !context.iter().copied().all(canonical) {
                scan.ambiguous_windows += 1;
                continue;
            }
            let gt = context.iter().filter(|b| matches!(b, b'G' | b'T')).count();
            if gt >= p.min_downstream_gt {
                first.get_or_insert(boundary);
                last = boundary;
                maximum_gt = maximum_gt.max(gt);
                possible += 1;
            }
        }
        if let Some(first) = first {
            let components = BTreeMap::from([
                ("maximum_downstream_gt".into(), maximum_gt as f64),
                ("passing_context_windows".into(), possible as f64),
            ]);
            scan.rows.push(row(
                sequence,
                start,
                signal_end,
                reverse,
                TerminationEvidenceKind::MammalianPolIiCandidate,
                TerminationGeometryKind::ProcessingSignal,
                "Mammalian Pol II signal/context candidate (uncalibrated)",
                components,
                vec![TerminationInterval {
                    geometry_kind: TerminationGeometryKind::CleavageWindow,
                    coordinate_basis: TerminationCoordinateBasis::InterbaseBoundarySet,
                    start_0based: first,
                    end_0based_exclusive: last + 1,
                }],
            ));
        }
    }
    if insufficient
        && !scan
            .warnings
            .iter()
            .any(|v| v.starts_with("Some mammalian"))
    {
        scan.warnings.push("Some mammalian signal candidates lack complete downstream context; missing windows are unassessed. Cleavage windows are positional hypotheses, not observed endpoints.".into());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // All sequences below are hand-crafted, non-biological fixtures, recreated by
    // the literals in these tests and used only for detector geometry/bounds.
    fn request() -> TerminationScreenRequest {
        TerminationScreenRequest {
            bacterial_intrinsic: true,
            ..Default::default()
        }
    }
    #[test]
    fn intrinsic_positive_reverse_negative_and_unknown() {
        let dna = b"GCGCGCAAAAGCGCGCTTTTTTTT";
        let a = predict_termination(dna, &request()).unwrap();
        assert!(a.rows.iter().any(|r| !r.reverse && r.start_0based == 0));
        let rc: Vec<_> = dna.iter().rev().map(|b| complement(*b)).collect();
        let b = predict_termination(&rc, &request()).unwrap();
        assert!(
            b.rows
                .iter()
                .any(|r| r.reverse && r.end_0based_exclusive == dna.len())
        );
        assert!(
            predict_termination(b"ACACACAAAACACACATTTTTTTT", &request())
                .unwrap()
                .rows
                .is_empty()
        );
        assert!(
            predict_termination(b"GCGNGCAAAAGCGCGCTTTTTTTT", &request())
                .unwrap()
                .rows
                .is_empty()
        );
    }
    #[test]
    fn mammalian_requires_context_keeps_window_separate() {
        let request = TerminationScreenRequest {
            mammalian_pol_ii: true,
            ..Default::default()
        };
        let dna = [b"AATAAA".as_slice(), &[b'C'; 10], &[b'T'; 50]].concat();
        let scan = predict_termination(&dna, &request).unwrap();
        let row = scan.rows.iter().find(|r| !r.reverse).unwrap();
        assert_eq!((row.start_0based, row.end_0based_exclusive), (0, 6));
        assert_eq!(
            row.related_intervals[0].geometry_kind,
            TerminationGeometryKind::CleavageWindow
        );
        assert_eq!(
            (
                row.related_intervals[0].start_0based,
                row.related_intervals[0].end_0based_exclusive
            ),
            (16, 37)
        );
        assert!(
            predict_termination(&[b"AATAAA".as_slice(), &[b'C'; 60]].concat(), &request)
                .unwrap()
                .rows
                .is_empty()
        );
        assert!(
            predict_termination(b"AATAAA", &request)
                .unwrap()
                .warnings
                .iter()
                .any(|v| v.contains("lack complete"))
        );
        let rc: Vec<_> = dna.iter().rev().map(|b| complement(*b)).collect();
        let reverse = predict_termination(&rc, &request).unwrap();
        let r = reverse.rows.iter().find(|r| r.reverse).unwrap();
        assert_eq!(
            (r.start_0based, r.end_0based_exclusive),
            (dna.len() - 6, dna.len())
        );
        assert_eq!(
            r.related_intervals[0].coordinate_basis,
            TerminationCoordinateBasis::InterbaseBoundarySet
        );
        assert_eq!(r.related_intervals[0].start_0based, dna.len() - 36);
        assert_eq!(
            r.related_intervals[0].end_0based_exclusive,
            dna.len() - 16 + 1
        );
        let forward_cuts = (row.related_intervals[0].start_0based
            ..row.related_intervals[0].end_0based_exclusive)
            .collect::<std::collections::BTreeSet<_>>();
        let reverse_cuts_on_forward = (r.related_intervals[0].start_0based
            ..r.related_intervals[0].end_0based_exclusive)
            .map(|p| dna.len() - p)
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(forward_cuts, reverse_cuts_on_forward);
    }
    #[test]
    fn budgets_validation_and_determinism() {
        let mut r = request();
        r.max_work = 1;
        let scan = predict_termination(&[b'G'; 100], &r).unwrap();
        assert!(scan.truncated);
        assert_eq!(scan.work, 1);
        r.bacterial_policy.max_stem_bp = usize::MAX;
        assert!(predict_termination(b"ACGT", &r).is_err());
        let a = predict_termination(b"GCGCGCAAAAGCGCGCTTTTTTTT", &request()).unwrap();
        let b = predict_termination(b"GCGCGCAAAAGCGCGCTTTTTTTT", &request()).unwrap();
        assert_eq!(a.rows, b.rows);
        assert_eq!(a.work, b.work);
    }
    #[test]
    fn multiple_signals_and_ambiguous_context_remain_explicit() {
        let request = TerminationScreenRequest {
            mammalian_pol_ii: true,
            ..Default::default()
        };
        let unit = [b"AATAAA".as_slice(), &[b'C'; 10], &[b'T'; 50]].concat();
        let dna = [unit.as_slice(), unit.as_slice()].concat();
        let scan = predict_termination(&dna, &request).unwrap();
        let starts = scan
            .rows
            .iter()
            .filter(|r| !r.reverse)
            .map(|r| r.start_0based)
            .collect::<Vec<_>>();
        assert_eq!(starts, vec![0, unit.len()]);
        let ambiguous = [b"AATAAA".as_slice(), &[b'N'; 60]].concat();
        let scan = predict_termination(&ambiguous, &request).unwrap();
        assert!(scan.rows.is_empty());
        assert!(scan.ambiguous_windows > 0);
        let short = predict_termination(b"AATAAA", &request).unwrap();
        assert!(short.missing_context_windows > 0);
    }

    #[test]
    fn empty_short_inputs_and_excessive_policy_values_are_bounded() {
        for dna in [b"".as_slice(), b"A", b"AATAA", b"N"] {
            let request = TerminationScreenRequest {
                mammalian_pol_ii: true,
                bacterial_intrinsic: true,
                ..Default::default()
            };
            assert!(predict_termination(dna, &request).unwrap().rows.is_empty());
        }
        let mut request = TerminationScreenRequest {
            mammalian_pol_ii: true,
            ..Default::default()
        };
        request.mammalian_policy.max_cleavage_offset_bp = usize::MAX;
        assert!(predict_termination(b"AATAAA", &request).is_err());
        request.mammalian_policy = MammalianPolyAPolicy::default();
        request.mammalian_policy.signals.push("AATAAA".into());
        assert!(predict_termination(b"AATAAA", &request).is_err());
        request.mammalian_policy = MammalianPolyAPolicy::default();
        request.mammalian_policy.min_downstream_gt = 0;
        assert!(predict_termination(b"AATAAA", &request).is_err());
        let mut bacterial_request = super::tests::request();
        bacterial_request.bacterial_policy.tract_window_bp = usize::MAX;
        assert!(predict_termination(b"ACGT", &bacterial_request).is_err());
    }
}
