//! Bounded canonical RT-tail generation composed with shared isoform comparison and seed assessment.

use super::rt_primer_comparison::{classify, metrics};
use super::*;

fn invalid(message: &str) -> EngineError {
    EngineError::invalid_input(message)
}
fn digest<T: Serialize>(value: &T) -> Result<String, EngineError> {
    let value = serde_json::to_value(value).map_err(|e| EngineError::internal(e.to_string()))?;
    Ok(sha256_prefixed_bytes(value.to_string().as_bytes()))
}
fn poll(
    progress: &mut dyn FnMut(OperationProgress) -> bool,
    completed: usize,
    total: usize,
) -> Result<(), EngineError> {
    if progress(OperationProgress::Workflow { completed, total }) {
        Ok(())
    } else {
        Err(invalid("rt_redesign_cancelled"))
    }
}
fn charge(work: &mut u64, amount: u64, max: u64) -> Result<(), EngineError> {
    *work = work
        .checked_add(amount)
        .filter(|n| *n <= max)
        .ok_or_else(|| invalid("rt_redesign_generation_work_budget"))?;
    Ok(())
}
fn member_key(m: &RtPrimerComparisonMember) -> RtPrimerRedesignMember {
    RtPrimerRedesignMember {
        seq_id: m.seq_id.clone(),
        transcript_id: m.transcript_id.clone(),
    }
}
fn primer(request: &RedesignRtPrimersRequest, tail: &str) -> RnaLibraryPrimer {
    let handle = &request.common_handle;
    let sequence = format!("{handle}{tail}");
    RnaLibraryPrimer {
        primer_id: short_sha256_id("rt_design", &sequence),
        assay_label: format!("RT redesign {}", request.group_id),
        sequence,
        components: vec![
            RnaPrimerComponent {
                component_id: "handle".into(),
                role: "common_handle".into(),
                start: 0,
                end: handle.len(),
            },
            RnaPrimerComponent {
                component_id: "tail".into(),
                role: "target_core".into(),
                start: handle.len(),
                end: handle.len() + tail.len(),
            },
        ],
        core_orientation: PrimerCoreOrientation::Antisense3prime,
        expected_targets: vec![],
        variant: "generated-v1".into(),
        source: format!(
            "Constrained mature-transcript enumeration for group {}",
            request.group_id
        ),
        source_sha256: None,
    }
}
fn admission(r: &RedesignRtPrimersRequest) -> Result<(), EngineError> {
    if !(1..=16).contains(&r.sources.len())
        || !(1..=16).contains(&r.references.len())
        || !(1..=256).contains(&r.required_members.len())
    {
        return Err(invalid("rt_redesign_schema_group_constraints_or_budget"));
    }
    let c = &r.constraints;
    let members = r.required_members.iter().collect::<BTreeSet<_>>();
    if r.schema != RT_PRIMER_REDESIGN_REQUEST_SCHEMA
        || r.group_id.is_empty()
        || r.group_id.len() > 128
        || !r
            .group_id
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_.-".contains(&b))
        || !(1..=256).contains(&r.required_members.len())
        || members.len() != r.required_members.len()
        || !members.contains(&r.anchor.member)
        || r.required_members.iter().any(|m| {
            m.seq_id.is_empty()
                || m.seq_id.len() > 512
                || m.transcript_id.is_empty()
                || m.transcript_id.len() > 512
        })
        || r.sources.iter().any(|s| {
            s.source.coverage_universe.kind
                != TranscriptAssayCoverageUniverseKind::ExplicitTranscripts
                || !(1..=256).contains(&s.source.coverage_universe.required_transcript_ids.len())
        })
        || r.common_handle.is_empty()
        || r.common_handle.len() > 54
        || !r.common_handle.bytes().all(|b| b"ACGT".contains(&b))
        || c.min_tail_length_bp < 10
        || c.max_tail_length_bp < c.min_tail_length_bp
        || c.max_tail_length_bp > 64 - r.common_handle.len()
        || !c.min_tail_gc_fraction.is_finite()
        || !c.max_tail_gc_fraction.is_finite()
        || !(0.0..=1.0).contains(&c.min_tail_gc_fraction)
        || !(c.min_tail_gc_fraction..=1.0).contains(&c.max_tail_gc_fraction)
        || !c.min_tail_tm_c.is_finite()
        || !c.max_tail_tm_c.is_finite()
        || c.min_tail_tm_c > c.max_tail_tm_c
        || c.max_homopolymer_run_bp > 64
        || c.max_self_complementary_run_bp > 64
        || c.max_self_3prime_complementary_run_bp > 64
        || !(1..=100_000).contains(&r.budgets.max_enumerated_windows)
        || !(1..=1_000_000_000).contains(&r.budgets.max_generation_work)
        || !(1..=16).contains(&r.budgets.max_candidates)
        || r.seed_limits.len() > 64
        || serde_json::to_vec(r)
            .map_err(|e| EngineError::internal(e.to_string()))?
            .len()
            > 1_000_000
    {
        return Err(invalid("rt_redesign_schema_group_constraints_or_budget"));
    }
    let mut limits = BTreeSet::new();
    for limit in &r.seed_limits {
        if ![6, 8, 10, 12].contains(&limit.seed_length)
            || limit.max_mismatches > limit.seed_length
            || !r
                .references
                .iter()
                .any(|v| v.reference_id == limit.reference_id)
            || !limits.insert((&limit.reference_id, limit.seed_length, limit.max_mismatches))
        {
            return Err(invalid("rt_redesign_seed_limit"));
        }
    }
    Ok(())
}
fn metric_pass(
    tail: &RtPrimerComparisonMetrics,
    full: &RtPrimerComparisonMetrics,
    c: &RtPrimerRedesignConstraints,
) -> bool {
    (c.min_tail_gc_fraction..=c.max_tail_gc_fraction).contains(&tail.gc_fraction)
        && (c.min_tail_tm_c..=c.max_tail_tm_c).contains(&tail.tm_c)
        && [tail, full].iter().all(|m| {
            m.longest_homopolymer_run_bp <= c.max_homopolymer_run_bp
                && m.self_complementary_run_bp <= c.max_self_complementary_run_bp
                && m.self_3prime_complementary_run_bp <= c.max_self_3prime_complementary_run_bp
        })
}

impl GentleEngine {
    /// Enumerate a declared interval and group; compare a deterministic capped shortlist without mutation.
    pub fn redesign_rt_primers(
        &self,
        request: RedesignRtPrimersRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<RtPrimerRedesignReport, EngineError> {
        admission(&request)?;
        poll(progress, 0, 1)?;
        let mut comparison_request = CompareRtPrimerCandidatesRequest {
            schema: RT_PRIMER_COMPARISON_REQUEST_SCHEMA.into(),
            assessment: AssessRtPrimerCandidatesRequest {
                schema: RT_PRIMER_ASSESSMENT_REQUEST_SCHEMA.into(),
                primers: vec![primer(
                    &request,
                    &"A".repeat(request.constraints.min_tail_length_bp),
                )],
                references: request.references.clone(),
                intended_sites: vec![],
                budgets: request.assessment_budgets.clone(),
            },
            sources: request.sources.clone(),
            budgets: request.comparison_budgets.clone(),
        };
        Self::validate_rt_primer_assessment_request(&comparison_request.assessment)?;
        let (templates, members, sources) =
            self.prepare_rt_comparison(&comparison_request, progress)?;
        let required = request
            .required_members
            .iter()
            .map(|m| {
                members
                    .iter()
                    .position(|v| &member_key(v) == m)
                    .ok_or_else(|| invalid("rt_redesign_required_member_not_resolved"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let anchor_index = members
            .iter()
            .position(|m| member_key(m) == request.anchor.member)
            .ok_or_else(|| invalid("rt_redesign_anchor_not_resolved"))?;
        let anchor = &templates[anchor_index];
        let start = request.anchor.start_0based;
        let end = request.anchor.end_0based_exclusive;
        if start >= end || end > anchor.sequence.len() {
            return Err(invalid("rt_redesign_anchor_interval"));
        }
        let mut report = RtPrimerRedesignReport {
            schema: RT_PRIMER_REDESIGN_SCHEMA.into(), content_sha256:String::new(), request_sha256:digest(&request)?, request:request.clone(),
            sources, required_member_ids:required.iter().map(|&i|members[i].member_id.clone()).collect(), members,
            enumerated_windows:0, skipped_ambiguous_windows:0, duplicate_tail_windows:0, metric_rejected_tails:0, group_rejected_tails:0,
            qualifying_distinct_tails:0, generation_work:0, shortlist_truncated:false, candidates:vec![], comparison:None,
            warnings:vec![
                "Required members form an explicit group, not pan-gene coverage or an expression denominator. Other selected members are comparison-only.".into(),
                "Every required member needs a canonical exact full-tail site. Three-prime-UTR mode additionally requires every exact or ambiguity-compatible group-member site to lie wholly after its resolved CDS; unknown annotation is never UTR.".into(),
                "Candidates are deduplicated tails in anchor start-coordinate then length order, not ranked for optimality. Additional qualifying tails outside the capped shortlist remain background-unassessed; failed seed checks do not refill the shortlist.".into(),
                "Seed limits use the all-sites marginal mismatch histogram, including intended sites. No genomic/background identity or intended-site exclusion is inferred from transcript geometry; ambiguous or unavailable evidence cannot pass a declared limit.".into(),
                "No declared seed limits is a vacuous constraint pass, not a specificity claim. Missing backgrounds, hairpin energies, abundance, retrospective mechanisms, experimental validation and order readiness remain independent.".into(),
                "An empty generation result has no comparison and does not open background files. Resource exhaustion is an error, never an empty biological result.".into(),
            ],
        };
        // Admit the entire canonical intake before collecting or scoring windows.
        for len in request.constraints.min_tail_length_bp..=request.constraints.max_tail_length_bp {
            let n = (end - start).checked_sub(len).map_or(0, |n| n + 1);
            report.enumerated_windows = report
                .enumerated_windows
                .checked_add(n)
                .filter(|n| *n <= request.budgets.max_enumerated_windows)
                .ok_or_else(|| invalid("rt_redesign_enumeration_budget"))?;
            charge(
                &mut report.generation_work,
                (n as u64) * (len as u64),
                request.budgets.max_generation_work,
            )?;
        }
        let mut windows = Vec::new();
        for len in request.constraints.min_tail_length_bp..=request.constraints.max_tail_length_bp {
            poll(progress, windows.len(), report.enumerated_windows)?;
            if len > end - start {
                continue;
            }
            let (mut canonical, ambiguous) = Self::enumerate_canonical_transcript_interval_windows(
                anchor,
                start,
                end,
                end - start,
                len,
            )?;
            report.skipped_ambiguous_windows += ambiguous;
            windows.append(&mut canonical);
        }
        windows.sort_by_key(|w| (w.local_start_0based, w.local_end_0based_exclusive));
        let mut seen = BTreeSet::new();
        let mut generated = Vec::new();
        for (index, window) in windows.iter().enumerate() {
            poll(progress, index, windows.len())?;
            let tail = Self::reverse_complement(&window.sequence_5_to_3);
            if !seen.insert(tail.clone()) {
                report.duplicate_tail_windows += 1;
                continue;
            }
            let oligo = primer(&request, &tail);
            // Suffix-in-target searches can be cubic; admit both tail/full heuristic upper bounds.
            charge(
                &mut report.generation_work,
                2 * (oligo.sequence.len() as u64).pow(3) + 8 * (oligo.sequence.len() as u64).pow(2),
                request.budgets.max_generation_work,
            )?;
            if !metric_pass(
                &metrics(tail.as_bytes()),
                &metrics(oligo.sequence.as_bytes()),
                &request.constraints,
            ) {
                report.metric_rejected_tails += 1;
                continue;
            }
            let word = window.sequence_5_to_3.as_bytes();
            let mut compatible = true;
            for &member_index in &required {
                let sequence = templates[member_index].sequence.as_bytes();
                let n = sequence.len().checked_sub(word.len()).map_or(0, |n| n + 1);
                // Exact equality and IUPAC compatibility each consume at most len comparisons.
                charge(
                    &mut report.generation_work,
                    2 * (n as u64) * (word.len() as u64),
                    request.budgets.max_generation_work,
                )?;
                let mut exact = false;
                for pos in 0..n {
                    if pos % 1024 == 0 {
                        poll(progress, index, windows.len())?;
                    }
                    let observed = &sequence[pos..pos + word.len()];
                    let equal = observed == word;
                    exact |= equal;
                    if request.mode == RtPrimerRedesignMode::ThreePrimeUtr
                        && (equal
                            || observed.iter().zip(word).all(|(&a, &b)| {
                                !IupacCode::from_letter(a)
                                    .subset(IupacCode::from_letter(b))
                                    .is_empty()
                            }))
                        && classify(
                            pos,
                            pos + word.len(),
                            &report.members[member_index].cds_ranges_0based,
                        )
                        .region
                            != RtPrimerBindingRegion::ThreePrimeUtr
                    {
                        compatible = false;
                        break;
                    }
                }
                if !exact {
                    compatible = false;
                }
                if !compatible {
                    break;
                }
            }
            if !compatible {
                report.group_rejected_tails += 1;
                continue;
            }
            report.qualifying_distinct_tails += 1;
            if generated.len() < request.budgets.max_candidates {
                report.candidates.push(RtPrimerRedesignCandidate {
                    primer_id: oligo.primer_id.clone(),
                    target_tail: tail,
                    anchor_start_0based: window.local_start_0based,
                    anchor_end_0based_exclusive: window.local_end_0based_exclusive,
                    seed_constraint_status: RtPrimerRedesignConstraintStatus::Passed,
                    seed_checks: vec![],
                });
                generated.push(oligo);
            }
        }
        report.shortlist_truncated = report.qualifying_distinct_tails > generated.len();
        if !generated.is_empty() {
            comparison_request.assessment.primers = generated;
            let comparison = self.compare_rt_primer_candidates(comparison_request, progress)?;
            for candidate in &mut report.candidates {
                for limit in &request.seed_limits {
                    let cell = comparison.assessment.cells.iter().find(|c| {
                        c.primer_id == candidate.primer_id
                            && c.reference_id == limit.reference_id
                            && c.seed_length == limit.seed_length
                    });
                    let counts = cell
                        .filter(|c| c.status == "completed")
                        .map(|c| &c.all_sites);
                    // A verified all-ambiguous reference retains its observed ambiguity count
                    // even though the existing assessment correctly leaves the cell unassessed.
                    let ambiguous_windows = cell
                        .filter(|_| {
                            candidate.target_tail.len() >= limit.seed_length
                                && comparison.assessment.references.iter().any(|r| {
                                    r.reference_id == limit.reference_id && r.status == "verified"
                                })
                        })
                        .map(|c| c.all_sites.ambiguous_windows);
                    let count = counts.map(|c| {
                        c.mismatch_histogram
                            .iter()
                            .take(limit.max_mismatches + 1)
                            .sum::<u64>()
                    });
                    let status = if count.is_some_and(|n| n > limit.max_compatible_windows) {
                        RtPrimerRedesignConstraintStatus::Failed
                    } else if counts.is_none_or(|c| c.ambiguous_windows > 0) {
                        RtPrimerRedesignConstraintStatus::Unresolved
                    } else {
                        RtPrimerRedesignConstraintStatus::Passed
                    };
                    candidate.seed_checks.push(RtPrimerRedesignSeedCheck {
                        limit: limit.clone(),
                        compatible_windows: count,
                        ambiguous_windows,
                        status,
                    });
                }
                candidate.seed_constraint_status = if candidate
                    .seed_checks
                    .iter()
                    .any(|c| c.status == RtPrimerRedesignConstraintStatus::Failed)
                {
                    RtPrimerRedesignConstraintStatus::Failed
                } else if candidate
                    .seed_checks
                    .iter()
                    .any(|c| c.status == RtPrimerRedesignConstraintStatus::Unresolved)
                {
                    RtPrimerRedesignConstraintStatus::Unresolved
                } else {
                    RtPrimerRedesignConstraintStatus::Passed
                };
            }
            report.comparison = Some(Box::new(comparison));
        }
        poll(
            progress,
            report.enumerated_windows,
            report.enumerated_windows,
        )?;
        report.content_sha256 = digest(&report)?;
        if serde_json::to_vec(&report)
            .map_err(|e| EngineError::internal(e.to_string()))?
            .len()
            > 8_000_000
        {
            return Err(invalid("rt_redesign_response_budget"));
        }
        Ok(report)
    }
}

#[cfg(test)]
#[path = "rt_primer_redesign_tests.rs"]
mod tests;
