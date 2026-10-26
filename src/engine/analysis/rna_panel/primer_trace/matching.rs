//! Bounded primer alignments and local template concordance in original-read coordinates.

use super::*;
use bio::alignment::{Alignment, AlignmentOperation};
use bio::pattern_matching::myers::Myers;

pub(super) fn rc(sequence: &[u8]) -> Vec<u8> {
    bio::alphabets::dna::revcomp(sequence)
}
fn canonical(sequence: &[u8]) -> bool {
    sequence.iter().all(|b| b"ACGT".contains(b))
}

pub(super) struct Work {
    used: usize,
    limit: usize,
}
impl Work {
    pub(super) fn new(limit: usize) -> Self {
        Self { used: 0, limit }
    }
    fn add(&mut self, n: usize) -> Result<(), EngineError> {
        self.used = self
            .used
            .checked_add(n)
            .ok_or_else(|| bad("primer_trace_work_budget"))?;
        if self.used > self.limit {
            return Err(bad("primer_trace_work_budget: no hits silently discarded"));
        }
        Ok(())
    }
}

#[derive(Clone)]
struct Match {
    start: usize,
    end: usize,
    edits: usize,
    core: Option<(usize, usize, usize)>,
}

fn aligned_core(aln: &Alignment, core: Option<(usize, usize)>) -> Option<(usize, usize, usize)> {
    let (a, b) = core?;
    let (mut x, mut y, mut first, mut end, mut edits) = (0, aln.ystart, None, 0, 0);
    for op in &aln.operations {
        let inside = a <= x && x < b;
        match op {
            AlignmentOperation::Match | AlignmentOperation::Subst => {
                if inside {
                    first.get_or_insert(y);
                    end = y + 1;
                    edits += usize::from(*op == AlignmentOperation::Subst);
                }
                x += 1;
                y += 1;
            }
            // Rust-Bio Ins consumes x (the pattern); Del consumes y (the read).
            AlignmentOperation::Ins => {
                if inside {
                    first.get_or_insert(y);
                    end = y;
                    edits += 1;
                }
                x += 1;
            }
            AlignmentOperation::Del => {
                if inside {
                    first.get_or_insert(y);
                    end = y + 1;
                    edits += 1;
                }
                y += 1;
            }
            _ => return None,
        }
    }
    first.map(|start| (start, end, edits))
}

fn windows(n: usize, width: usize, overlap: usize) -> Vec<(usize, usize)> {
    let mut result = vec![];
    let mut start = 0;
    while start < n {
        let end = (start + width).min(n);
        result.push((start, end));
        if end == n {
            break;
        }
        start = end - overlap;
    }
    result
}

/// Exact matches cover the whole subject. Approximate matches use explicit bounded windows.
fn find(
    pattern: &[u8],
    text: &[u8],
    edits: u8,
    core: Option<(usize, usize)>,
    ranges: &[(usize, usize)],
    cap: usize,
    collapse_read_occurrences: bool,
    work: &mut Work,
    cancel: &mut dyn FnMut() -> bool,
) -> Result<Vec<Match>, EngineError> {
    let mut hits = BTreeMap::new();
    work.add(text.len().saturating_mul(pattern.len()))?;
    if text.len() >= pattern.len() {
        for (start, word) in text.windows(pattern.len()).enumerate() {
            if start % 4096 == 0 && cancel() {
                return Err(bad("primer_trace_cancelled"));
            }
            if word == pattern {
                let end = start + pattern.len();
                hits.insert(
                    (start, end),
                    Match {
                        start,
                        end,
                        edits: 0,
                        core: core.map(|(a, b)| (start + a, start + b, 0)),
                    },
                );
                if hits.len() > cap {
                    return Err(bad("primer_trace_hit_budget"));
                }
            }
        }
    }
    if edits > 0 {
        let mut myers = Myers::<u64>::new(pattern);
        for &(a, b) in ranges {
            if cancel() {
                return Err(bad("primer_trace_cancelled"));
            }
            work.add((b - a).saturating_mul(pattern.len()))?;
            let mut matches = myers.find_all_lazy(&text[a..b], edits);
            while let Some((end, distance)) = matches.next() {
                if distance == 0 {
                    continue;
                }
                let mut aln = Alignment::default();
                if !matches.alignment_at(end, &mut aln) {
                    return Err(bad("primer_trace_alignment_unavailable"));
                }
                let start = a + aln.ystart;
                let end = a + aln.yend;
                if start == end {
                    continue;
                }
                let core = aligned_core(&aln, core).map(|(s, e, d)| (a + s, a + e, d));
                hits.entry((start, end)).or_insert(Match {
                    start,
                    end,
                    edits: distance as usize,
                    core,
                });
                if hits.len() > cap {
                    return Err(bad("primer_trace_hit_budget"));
                }
            }
        }
    }
    let mut hits = hits.into_values().collect::<Vec<_>>();
    hits.sort_by_key(|h| (h.edits, h.start, h.end));
    let mut kept: Vec<Match> = vec![];
    for hit in hits {
        if collapse_read_occurrences
            && kept.iter().any(|old| {
                (hit.edits > 0 || old.edits > 0) && hit.start < old.end && old.start < hit.end
            })
        {
            continue;
        }
        kept.push(hit);
    }
    Ok(kept)
}

pub(super) fn core(primer: &RnaLibraryPrimer) -> &RnaPrimerComponent {
    primer
        .components
        .iter()
        .find(|c| c.role == "target_core")
        .expect("validated core")
}

fn rank(
    hit: &RnaPrimerHit,
) -> (
    bool,
    usize,
    usize,
    usize,
    Option<&str>,
    PrimerReadStrand,
    usize,
) {
    (
        !hit.terminal,
        if hit.anchor_class.starts_with("full_") {
            0
        } else if hit.anchor_class.starts_with("core_") {
            1
        } else {
            2
        },
        hit.edits,
        hit.start,
        hit.primer_id.as_deref(),
        hit.strand,
        hit.end,
    )
}

pub(super) fn anchors(
    panel: &RnaPrimerPanel,
    read: &[u8],
    budgets: &RnaPrimerTraceBudgets,
    work: &mut Work,
    cancel: &mut dyn FnMut() -> bool,
) -> Result<Vec<RnaPrimerHit>, EngineError> {
    let window = panel.policy.end_window_nt.min(read.len());
    let ranges = if window == read.len() {
        vec![(0, read.len())]
    } else {
        vec![(0, window), (read.len() - window, read.len())]
    };
    let mut result = vec![];
    let mut handles = BTreeSet::new();
    for primer in &panel.primers {
        let c = core(primer);
        for (strand, reverse) in [
            (PrimerReadStrand::Forward, false),
            (PrimerReadStrand::ReverseComplement, true),
        ] {
            for full in [true, false] {
                let original = if full {
                    primer.sequence.as_bytes()
                } else {
                    &primer.sequence.as_bytes()[c.start..c.end]
                };
                let pattern = if reverse {
                    rc(original)
                } else {
                    original.to_vec()
                };
                let component = if full {
                    if reverse {
                        (pattern.len() - c.end, pattern.len() - c.start)
                    } else {
                        (c.start, c.end)
                    }
                } else {
                    (0, pattern.len())
                };
                let limit = if full {
                    panel.policy.max_full_edits
                } else {
                    panel.policy.max_core_edits_for_attribution
                };
                for h in find(
                    &pattern,
                    read,
                    limit,
                    Some(component),
                    &ranges,
                    budgets.max_hits_per_read,
                    true,
                    work,
                    cancel,
                )? {
                    let (cs, ce, cd) = h.core.ok_or_else(|| bad("primer_trace_core_alignment"))?;
                    result.push(RnaPrimerHit {
                        primer_id: Some(primer.primer_id.clone()),
                        anchor_class: format!(
                            "{}_{}",
                            if full { "full" } else { "core" },
                            if h.edits == 0 { "exact" } else { "approximate" }
                        ),
                        strand,
                        start: h.start,
                        end: h.end,
                        edits: h.edits,
                        core_start: Some(cs),
                        core_end: Some(ce),
                        core_edits: Some(cd),
                        terminal: if reverse {
                            read.len() - h.end <= panel.policy.end_window_nt
                        } else {
                            h.start <= panel.policy.end_window_nt
                        },
                    });
                }
            }
        }
        for h in primer
            .components
            .iter()
            .filter(|c| c.role == "common_handle")
        {
            handles.insert(primer.sequence[h.start..h.end].to_string());
        }
    }
    let full = result
        .iter()
        .filter(|h| h.anchor_class.starts_with("full_"))
        .cloned()
        .collect::<Vec<_>>();
    result.retain(|h| {
        h.anchor_class.starts_with("full_")
            || !full.iter().any(|f| {
                f.primer_id == h.primer_id
                    && f.strand == h.strand
                    && f.start <= h.start
                    && h.end <= f.end
            })
    });
    for handle in handles {
        for (strand, reverse) in [
            (PrimerReadStrand::Forward, false),
            (PrimerReadStrand::ReverseComplement, true),
        ] {
            let pattern = if reverse {
                rc(handle.as_bytes())
            } else {
                handle.as_bytes().to_vec()
            };
            for h in find(
                &pattern,
                read,
                0,
                None,
                &[],
                budgets.max_hits_per_read,
                true,
                work,
                cancel,
            )? {
                if full
                    .iter()
                    .any(|f| f.strand == strand && f.start <= h.start && h.end <= f.end)
                {
                    continue;
                }
                result.push(RnaPrimerHit {
                    primer_id: None,
                    anchor_class: "handle_only".into(),
                    strand,
                    start: h.start,
                    end: h.end,
                    edits: 0,
                    core_start: None,
                    core_end: None,
                    core_edits: None,
                    terminal: if reverse {
                        read.len() - h.end <= panel.policy.end_window_nt
                    } else {
                        h.start <= panel.policy.end_window_nt
                    },
                });
            }
        }
    }
    result.sort_by(|a, b| rank(a).cmp(&rank(b)));
    result.dedup();
    if result.len() > budgets.max_hits_per_read {
        return Err(bad("primer_trace_hit_budget"));
    }
    Ok(result)
}

pub(super) struct Site {
    primer: String,
    member: usize,
    start: usize,
    end: usize,
}

#[derive(Default)]
pub(super) struct ReferencePreparation {
    sites: Vec<Site>,
    index: Index,
    seed_nt: usize,
}
fn expected(primer: &RnaLibraryPrimer, template: &RnaReadPanelTemplate) -> bool {
    primer.expected_targets.iter().any(|t| {
        t.locus_id == template.locus_id
            && t.transcript_feature_id
                .is_none_or(|id| id == template.transcript_feature_id)
    })
}

pub(super) fn sites(
    primers: &RnaPrimerPanel,
    reference: &RnaReadPanel,
    budgets: &RnaPrimerTraceBudgets,
    cancel: &mut dyn FnMut() -> bool,
) -> Result<(ReferencePreparation, Vec<String>), EngineError> {
    let mut sites = vec![];
    let mut warnings = vec![];
    let mut work = Work::new(budgets.max_scan_cells);
    for primer in &primers.primers {
        let c = core(primer);
        let original = &primer.sequence.as_bytes()[c.start..c.end];
        let pattern = if primer.core_orientation == PrimerCoreOrientation::Antisense3prime {
            rc(original)
        } else {
            original.to_vec()
        };
        for (member, t) in reference
            .templates
            .iter()
            .enumerate()
            .filter(|(_, t)| expected(primer, t))
        {
            let ranges = windows(t.sequence.len(), 4096, pattern.len() + 3);
            for h in find(
                &pattern,
                t.sequence.as_bytes(),
                primers.policy.max_core_edits_for_attribution,
                None,
                &ranges,
                budgets.max_reference_sites,
                false,
                &mut work,
                cancel,
            )? {
                if !canonical(&t.sequence.as_bytes()[h.start..h.end]) {
                    continue;
                }
                sites.push(Site {
                    primer: primer.primer_id.clone(),
                    member,
                    start: h.start,
                    end: h.end,
                });
                if sites.len() > budgets.max_reference_sites {
                    return Err(bad("primer_trace_reference_site_budget"));
                }
            }
        }
        if !sites.iter().any(|s| s.primer == primer.primer_id) {
            warnings.push(format!("primer_core_site_unavailable:{}", primer.primer_id));
        }
        for target in &primer.expected_targets {
            if !reference.templates.iter().any(|t| {
                t.locus_id == target.locus_id
                    && target
                        .transcript_feature_id
                        .is_none_or(|id| id == t.transcript_feature_id)
            }) {
                return Err(bad("primer_trace_unknown_expected_target"));
            }
        }
    }
    if cancel() {
        return Err(bad("primer_trace_cancelled"));
    }
    let policy = RnaReadPanelPolicy {
        kmer_len: primers.policy.concordance_seed_nt,
        max_index_postings: budgets.max_index_postings,
        ..RnaReadPanelPolicy::default()
    };
    work.add(
        reference
            .templates
            .iter()
            .map(|t| t.sequence.len())
            .sum::<usize>()
            .saturating_mul(policy.kmer_len),
    )?;
    let (index, postings) = build_index(reference, &policy)?;
    work.add(postings)?;
    if cancel() {
        return Err(bad("primer_trace_cancelled"));
    }
    Ok((
        ReferencePreparation {
            sites,
            index,
            seed_nt: policy.kmer_len,
        },
        warnings,
    ))
}

/// Global, unit-cost Levenshtein distance; ambiguous bases cannot establish concordance.
fn distance(a: &[u8], b: &[u8], work: &mut Work) -> Result<usize, EngineError> {
    work.add(a.len().saturating_mul(b.len()))?;
    if !canonical(a) || !canonical(b) {
        return Ok(usize::MAX);
    }
    let mut prev = (0..=b.len()).collect::<Vec<_>>();
    for (i, x) in a.iter().enumerate() {
        let mut next = vec![i + 1; b.len() + 1];
        for (j, y) in b.iter().enumerate() {
            next[j + 1] = (prev[j] + usize::from(x != y))
                .min(prev[j + 1] + 1)
                .min(next[j] + 1);
        }
        prev = next;
    }
    Ok(prev[b.len()])
}

fn prefix(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}
fn oriented(read: &[u8], a: usize, b: usize, strand: PrimerReadStrand) -> Vec<u8> {
    if strand == PrimerReadStrand::ReverseComplement {
        rc(&read[a..b])
    } else {
        read[a..b].to_vec()
    }
}

pub(super) fn row(
    panel: &RnaPrimerPanel,
    reference: Option<&RnaReadPanel>,
    prepared: &ReferencePreparation,
    read: &input::FastqRecord,
    index: u64,
    budgets: &RnaPrimerTraceBudgets,
    cancel: &mut dyn FnMut() -> bool,
) -> Result<RnaPrimerTraceRow, EngineError> {
    let mut work = Work::new(budgets.max_scan_cells);
    if read.sequence.len() > budgets.max_read_nt {
        return Err(bad("primer_trace_read_budget"));
    }
    let hits = anchors(panel, &read.sequence, budgets, &mut work, cancel)?;
    let mut row = RnaPrimerTraceRow {
        record_index: index,
        read_id: read.read_id.clone(),
        read_length: read.sequence.len(),
        anchor_class: "no_primer".into(),
        anchor_position: None,
        primer_multiplicity: None,
        locus_relation: if reference.is_some() {
            "unresolved"
        } else {
            "not_evaluated"
        }
        .into(),
        evidence_tier: "no_anchor".into(),
        anchor: hits.first().cloned(),
        hits,
        leader: None,
        payload: None,
        leader_length: None,
        payload_length: None,
        payload_transcript_sense: None,
        payload_mean_phred: None,
        payload_fraction_below_q10: None,
        template_matches: vec![],
        artifact_flags: vec![],
        reasons: vec![],
    };
    let Some(hit) = row.anchor.clone() else {
        row.reasons.push("no_primer_anchor".into());
        return Ok(row);
    };
    row.anchor_class = hit.anchor_class.clone();
    row.anchor_position = Some(if hit.terminal { "terminal" } else { "internal" }.into());
    if row
        .hits
        .iter()
        .any(|h| h.primer_id.is_some() && !h.terminal)
    {
        row.artifact_flags.push("internal_primer".into());
    }
    let ids = row
        .hits
        .iter()
        .filter_map(|h| h.primer_id.as_ref())
        .collect::<BTreeSet<_>>();
    let strands = row
        .hits
        .iter()
        .filter(|h| h.primer_id.is_some())
        .map(|h| h.strand)
        .collect::<BTreeSet<_>>();
    if !ids.is_empty() {
        row.primer_multiplicity = Some(
            if strands.len() > 1 {
                "conflicting_strands"
            } else if ids.len() > 1 {
                "multiple_primers"
            } else if row.hits.iter().filter(|h| h.primer_id.is_some()).count() > 1 {
                "repeated_same_primer"
            } else {
                "single"
            }
            .into(),
        );
    }
    if ids.len() > 1 {
        row.artifact_flags.push("multiple_primers".into());
        row.reasons
            .push("multiple_primers_primary_anchor_only".into());
    }
    if panel.policy.check_handle_followed_by_poly_t {
        let n = panel.policy.poly_t_min_run;
        for h in &row.hits {
            if h.anchor_class != "handle_only" {
                continue;
            }
            let tail = if h.strand == PrimerReadStrand::Forward {
                read.sequence.get(h.end..h.end + n).map(|s| s.to_vec())
            } else {
                h.start
                    .checked_sub(n)
                    .and_then(|a| read.sequence.get(a..h.start))
                    .map(rc)
            };
            if tail.is_some_and(|t| {
                canonical(&t)
                    && t.iter().filter(|&&b| b != b'T').count()
                        <= panel.policy.poly_t_max_mismatches
            }) {
                row.artifact_flags.push("handle_followed_by_poly_t".into());
                break;
            }
        }
    }
    let Some(id) = &hit.primer_id else {
        row.evidence_tier = "handle_observation".into();
        row.reasons.push("handle_has_no_target_identity".into());
        return Ok(row);
    };
    let primer = panel
        .primers
        .iter()
        .find(|p| &p.primer_id == id)
        .expect("validated hit");
    row.evidence_tier = "primer_observation".into();
    let (a, b, leader) = if hit.strand == PrimerReadStrand::Forward {
        (hit.end, read.sequence.len(), (0, hit.start))
    } else {
        (0, hit.start, (hit.end, read.sequence.len()))
    };
    row.leader = Some(leader);
    row.payload = Some((a, b));
    row.leader_length = Some(leader.1 - leader.0);
    row.payload_length = Some(b - a);
    row.payload_transcript_sense =
        Some(primer.core_orientation == PrimerCoreOrientation::Sense5prime);
    if a < b {
        let q = &read.quality[a..b];
        row.payload_mean_phred =
            Some(q.iter().map(|q| u64::from(q - 33)).sum::<u64>() as f64 / q.len() as f64);
        row.payload_fraction_below_q10 =
            Some(q.iter().filter(|&&q| q < 43).count() as f64 / q.len() as f64);
    }
    if b - a < panel.policy.min_payload_nt {
        row.artifact_flags.push("short_payload".into());
        row.reasons.push("insufficient_payload".into());
        return Ok(row);
    }
    // A second observed anchor bounds interpretation. Phase 1 does not segment molecules.
    if row
        .hits
        .iter()
        .any(|h| h.primer_id.is_some() && a <= h.start && h.end <= b)
    {
        row.reasons
            .push("multiple_anchor_payload_not_segmented".into());
        return Ok(row);
    }
    let Some(reference) = reference else {
        return Ok(row);
    };
    if hit
        .core_edits
        .is_none_or(|e| e > panel.policy.max_core_edits_for_attribution as usize)
    {
        row.reasons.push("core_edit_limit".into());
        return Ok(row);
    }
    let payload = oriented(&read.sequence, a, b, hit.strand);
    let probe = &payload[..panel.policy.concordance_probe_nt];
    if !canonical(probe) {
        row.reasons.push("ambiguous_payload_probe".into());
        return Ok(row);
    }
    let observed_core = oriented(
        &read.sequence,
        hit.core_start.unwrap(),
        hit.core_end.unwrap(),
        hit.strand,
    );
    let mut assessed = false;
    let mut divergent = true;
    for site in prepared.sites.iter().filter(|s| &s.primer == id) {
        if cancel() {
            return Err(bad("primer_trace_cancelled"));
        }
        let template = &reference.templates[site.member];
        let seq = template.sequence.as_bytes();
        let anti = primer.core_orientation == PrimerCoreOrientation::Antisense3prime;
        let template_core = if anti {
            rc(&seq[site.start..site.end])
        } else {
            seq[site.start..site.end].to_vec()
        };
        if distance(&observed_core, &template_core, &mut work)?
            > panel.policy.max_core_edits_for_attribution as usize
        {
            continue;
        }
        let nominal = if anti {
            site.start
                .checked_sub(probe.len())
                .map(|a| &seq[a..site.start])
        } else {
            seq.get(site.end..site.end + probe.len())
        };
        if !nominal.is_some_and(canonical) {
            continue;
        }
        let mut had_flank = false;
        for width in probe.len().saturating_sub(panel.policy.max_probe_edits)
            ..=probe.len() + panel.policy.max_probe_edits
        {
            let flank = if anti {
                site.start
                    .checked_sub(width)
                    .map(|s| rc(&seq[s..site.start]))
            } else {
                seq.get(site.end..site.end + width).map(|s| s.to_vec())
            };
            let Some(flank) = flank else {
                continue;
            };
            if !canonical(&flank) {
                continue;
            }
            had_flank = true;
            let edits = distance(probe, &flank, &mut work)?;
            if edits <= panel.policy.max_probe_edits {
                row.template_matches.push(RnaPrimerTemplateMatch {
                    member_id: template.member_id.clone(),
                    template_start: if anti { site.start - width } else { site.start },
                    template_end: if anti { site.end } else { site.end + width },
                    expected_target: true,
                    core_and_probe_contiguous: true,
                    payload_direction_matches_declared: true,
                    probe_edits: edits,
                    concordant_nt: prefix(probe, &flank),
                    core_3prime_mismatches: {
                        let c = core(primer);
                        let oligo = &primer.sequence.as_bytes()[c.start..c.end];
                        (oligo.len() == template_core.len()).then(|| {
                            oligo[oligo.len() - 10..]
                                .iter()
                                .zip(&template_core[template_core.len() - 10..])
                                .filter(|(a, b)| a != b)
                                .count()
                        })
                    },
                });
                if row.template_matches.len() > budgets.max_reference_sites {
                    return Err(bad("primer_trace_reference_site_budget"));
                }
            }
            // Three discordant canonical bases exceed the default two-edit allowance;
            // a single basecall error never becomes a mechanism assertion.
            if probe
                .iter()
                .take(3)
                .zip(&flank)
                .filter(|(x, y)| x != y)
                .count()
                <= panel.policy.max_probe_edits
            {
                divergent = false;
            }
        }
        assessed |= had_flank;
    }
    if !row.template_matches.is_empty() {
        row.locus_relation = "expected_target_contiguous".into();
        row.evidence_tier = "local_template_concordance".into();
    } else {
        // Seeds nominate bounded candidate intervals only. Global verification below
        // consumes the entire probe; no seed match alone sets a locus relation.
        let mut candidates = BTreeSet::new();
        let edits = panel.policy.max_probe_edits;
        for anti in [false, true] {
            let query = if anti { rc(probe) } else { probe.to_vec() };
            for (offset, word) in query.windows(prepared.seed_nt).enumerate() {
                let Some(key) = GentleEngine::encode_kmer_bits(word) else {
                    continue;
                };
                let Some(postings) = prepared.index.get(&key) else {
                    continue;
                };
                work.add(postings.len())?;
                for &(member, position) in postings {
                    if cancel() {
                        return Err(bad("primer_trace_cancelled"));
                    }
                    let nominal = position as isize - offset as isize;
                    for drift in -(edits as isize)..=edits as isize {
                        let start = nominal + drift;
                        if start >= 0 {
                            candidates.insert((member as usize, start as usize, anti));
                            if candidates.len() > budgets.max_reference_sites {
                                return Err(bad("primer_trace_reference_candidate_budget"));
                            }
                        }
                    }
                }
            }
        }
        for (member, start, anti) in candidates {
            if cancel() {
                return Err(bad("primer_trace_cancelled"));
            }
            let template = &reference.templates[member];
            let seq = template.sequence.as_bytes();
            for width in probe.len() - edits..=probe.len() + edits {
                let end = start + width;
                let Some(raw) = seq.get(start..end) else {
                    continue;
                };
                if !canonical(raw) {
                    continue;
                }
                let flank = if anti { rc(raw) } else { raw.to_vec() };
                let found_edits = distance(probe, &flank, &mut work)?;
                if found_edits > edits {
                    continue;
                }
                let adjacent = if anti {
                    seq.get(end..end + 10).map(rc)
                } else {
                    start.checked_sub(10).map(|a| seq[a..start].to_vec())
                };
                let c = core(primer);
                let oligo = &primer.sequence.as_bytes()[c.end - 10..c.end];
                let mismatch = adjacent
                    .filter(|s| canonical(s))
                    .map(|s| oligo.iter().zip(s).filter(|(a, b)| **a != *b).count());
                let declared_anti =
                    primer.core_orientation == PrimerCoreOrientation::Antisense3prime;
                row.template_matches.push(RnaPrimerTemplateMatch {
                    member_id: template.member_id.clone(),
                    template_start: start,
                    template_end: end,
                    expected_target: expected(primer, template),
                    core_and_probe_contiguous: false,
                    payload_direction_matches_declared: anti == declared_anti,
                    probe_edits: found_edits,
                    concordant_nt: prefix(probe, &flank),
                    core_3prime_mismatches: mismatch,
                });
                if row.template_matches.len() > budgets.max_reference_sites {
                    return Err(bad("primer_trace_reference_site_budget"));
                }
            }
        }
        if row.template_matches.iter().any(|m| !m.expected_target) {
            row.locus_relation = "other_panel_template".into();
        } else if assessed && divergent {
            row.locus_relation = "boundary_divergent".into();
        }
        if assessed && divergent {
            row.artifact_flags.push("primer_boundary_divergence".into());
        }
        row.reasons.push(
            if !assessed {
                "expected_site_or_flank_unavailable"
            } else if row.template_matches.is_empty() {
                "no_panel_template"
            } else {
                "payload_match_without_core_concordance"
            }
            .into(),
        );
    }
    row.template_matches.sort_by(|a, b| {
        (
            &a.member_id,
            a.template_start,
            a.template_end,
            a.probe_edits,
        )
            .cmp(&(
                &b.member_id,
                b.template_start,
                b.template_end,
                b.probe_edits,
            ))
    });
    row.artifact_flags.sort();
    row.artifact_flags.dedup();
    Ok(row)
}
