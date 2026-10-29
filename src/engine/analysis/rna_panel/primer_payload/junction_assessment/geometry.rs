//! Exact physical geometry and explicit strand-conditioned motifs, not RNA/protein inference.

use super::*;

pub(super) fn assigned_orientation(rows: &[ObservedJunctionSourceAssessment]) -> Option<char> {
    let strand = rows.first()?.model_strand?;
    rows.iter()
        .all(|r| r.model_strand == Some(strand))
        .then_some(strand)
}

pub(super) fn junction_id(
    reference: &s1::ReferenceBinding,
    contig: &str,
    skip: s1::Interval,
) -> Result<String, EngineError> {
    hash(&(
        "gentle.observed_splice_junction.v1",
        reference,
        contig,
        skip,
    ))
}

pub(super) fn near_hits_for_models(
    skip: s1::Interval,
    models: &[PrimerPayloadTranscript],
    window: usize,
) -> Vec<ObservedAnnotatedJunctionHit> {
    let mut hits = BTreeMap::<s1::Interval, BTreeSet<String>>::new();
    for m in models {
        for exons in m.genomic_exons.windows(2) {
            let pair = (exons[0].1, exons[1].0);
            if pair.0 < pair.1
                && pair.0.abs_diff(skip.0) <= window
                && pair.1.abs_diff(skip.1) <= window
            {
                hits.entry(pair)
                    .or_default()
                    .insert(m.transcript_id.clone());
            }
        }
    }
    hits.into_iter()
        .map(|(p, ids)| ObservedAnnotatedJunctionHit {
            skipped_interval: p,
            start_offset: p.0 as i64 - skip.0 as i64,
            end_offset: p.1 as i64 - skip.1 as i64,
            model_ids: ids.into_iter().collect(),
        })
        .collect()
}

pub(super) fn near_relation(hits: &[ObservedAnnotatedJunctionHit]) -> &'static str {
    if hits
        .iter()
        .any(|h| h.start_offset == 0 && h.end_offset == 0)
    {
        return "exact";
    }
    let start = hits.iter().any(|h| h.start_offset == 0);
    let end = hits.iter().any(|h| h.end_offset == 0);
    match (start, end, hits.is_empty()) {
        (true, false, _) => "shared_start_only",
        (false, true, _) => "shared_end_only",
        (_, _, true) => "none_within_window",
        _ => "near",
    }
}

pub(super) fn geometry_events(skip: s1::Interval, m: &PrimerPayloadTranscript) -> Vec<String> {
    let pairs = m
        .genomic_exons
        .windows(2)
        .filter_map(|e| (e[0].1 < e[1].0).then_some((e[0].1, e[1].0)))
        .collect::<Vec<_>>();
    let start = pairs.iter().any(|p| p.0 == skip.0);
    let end = pairs.iter().any(|p| p.1 == skip.1);
    let mut labels = BTreeSet::new();
    let plus = m.genomic_strand == Some('+');
    if start && !end {
        labels.insert(if plus {
            "alternative_acceptor"
        } else {
            "alternative_donor"
        });
    }
    if end && !start {
        labels.insert(if plus {
            "alternative_donor"
        } else {
            "alternative_acceptor"
        });
    }
    if m.genomic_exons.iter().enumerate().any(|(i, e)| {
        e.1 == skip.0
            && m.genomic_exons
                .iter()
                .enumerate()
                .any(|(j, t)| j >= i + 2 && t.0 == skip.1)
    }) {
        labels.insert("exon_skip");
    }
    let first = if plus { pairs.first() } else { pairs.last() };
    let last = if plus { pairs.last() } else { pairs.first() };
    if first.is_some_and(|p| {
        if plus {
            p.1 == skip.1 && p.0 != skip.0
        } else {
            p.0 == skip.0 && p.1 != skip.1
        }
    }) {
        labels.insert("alternative_first_junction");
    }
    if last.is_some_and(|p| {
        if plus {
            p.0 == skip.0 && p.1 != skip.1
        } else {
            p.1 == skip.1 && p.0 != skip.0
        }
    }) {
        labels.insert("alternative_terminal_acceptor");
    }
    if !start && !end {
        labels.insert("unannotated_both_ends");
    }
    labels.into_iter().map(str::to_owned).collect()
}

fn classify(first: &str, last: &str) -> ObservedMotifClass {
    match (first, last) {
        ("GT", "AG") => ObservedMotifClass::GtAg,
        ("GC", "AG") => ObservedMotifClass::GcAg,
        ("AT", "AC") => ObservedMotifClass::AtAc,
        _ => ObservedMotifClass::Other,
    }
}

fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

pub(super) fn motif_from_dinucleotides(
    first: Option<&str>,
    last: Option<&str>,
    orientation: Option<char>,
    authority: &str,
) -> ObservedJunctionMotif {
    let mut out = ObservedJunctionMotif {
        reference_authority: authority.into(),
        reference_forward_first: first.map(str::to_owned),
        reference_forward_last: last.map(str::to_owned),
        reference_forward_motif_class: ObservedMotifClass::NotAssessable,
        reference_reverse_motif_class: ObservedMotifClass::NotAssessable,
        motif_class: ObservedMotifClass::NotAssessable,
        opposite_orientation_motif_class: ObservedMotifClass::NotAssessable,
        reason: "ambiguous_or_unavailable_bases".into(),
    };
    let (Some(first), Some(last)) = (first, last) else {
        return out;
    };
    let first = first.to_ascii_uppercase();
    let last = last.to_ascii_uppercase();
    if first.len() != 2
        || last.len() != 2
        || !first
            .bytes()
            .chain(last.bytes())
            .all(|b| b"ACGT".contains(&b))
    {
        return out;
    }
    out.reference_forward_motif_class = classify(&first, &last);
    out.reference_reverse_motif_class = classify(&rc(&last), &rc(&first));
    match orientation {
        Some('+') => {
            out.motif_class = out.reference_forward_motif_class;
            out.opposite_orientation_motif_class = out.reference_reverse_motif_class;
            out.reason = "assigned_model_orientation".into();
        }
        Some('-') => {
            out.motif_class = out.reference_reverse_motif_class;
            out.opposite_orientation_motif_class = out.reference_forward_motif_class;
            out.reason = "assigned_model_orientation".into();
        }
        _ => out.reason = "orientation_unassigned".into(),
    }
    out
}
