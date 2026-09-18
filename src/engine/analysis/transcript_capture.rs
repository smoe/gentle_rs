//! Bounded exact-site capture discovery over shared annotated cDNA templates.
//!
//! Coverage means an exact eligible binding site, not observed RT reach, PCR
//! yield, whole-genome specificity or an approved oligo order.

use super::*;

struct CaptureTemplate {
    target_index: usize,
    member: TranscriptCaptureMember,
    template: TranscriptQpcrDesignTemplate,
}

fn canonical(sequence: &[u8]) -> bool {
    sequence
        .iter()
        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
}

fn normalize_stages(stages: &mut Vec<String>) -> Result<(), EngineError> {
    if stages.is_empty() || stages.len() > 16 || stages.iter().any(|s| s.trim().is_empty()) {
        return Err(EngineError::invalid_input(
            "Capture oligos require 1..=16 nonempty stage_ids",
        ));
    }
    for stage in stages.iter_mut() {
        *stage = stage.trim().to_string();
    }
    stages.sort();
    stages.dedup();
    Ok(())
}

fn shared_stages(left: &[String], right: &[String]) -> Vec<String> {
    left.iter()
        .filter(|stage| right.contains(stage))
        .cloned()
        .collect()
}

/// Canonical oligos use the existing heuristic verbatim. Degenerate protocol
/// oligos use possible IUPAC matches, without exponential variant expansion.
fn dimer(left: &[u8], right: &[u8]) -> PrimerPairDimerMetrics {
    if canonical(left) && canonical(right) {
        return GentleEngine::compute_primer_pair_dimer_metrics(left, right);
    }
    let compatible = |a: u8, b: u8| {
        !IupacCode::from_letter(a)
            .subset(IupacCode::from_letter(b))
            .is_empty()
    };
    let rc = GentleEngine::reverse_complement_bytes(right);
    let mut max_run = 0;
    for shift in -(rc.len() as isize)..left.len() as isize {
        let mut run = 0;
        for (i, a) in left.iter().enumerate() {
            let j = i as isize - shift;
            if j >= 0 && (j as usize) < rc.len() && compatible(*a, rc[j as usize]) {
                run += 1;
                max_run = max_run.max(run);
            } else {
                run = 0;
            }
        }
    }
    let suffix = |a: &[u8], b: &[u8]| {
        (1..=a.len().min(b.len()))
            .rev()
            .find(|&len| {
                b.windows(len).any(|window| {
                    a[a.len() - len..]
                        .iter()
                        .zip(window)
                        .all(|(x, y)| compatible(*x, *y))
                })
            })
            .unwrap_or(0)
    };
    PrimerPairDimerMetrics {
        max_complementary_run_bp: max_run,
        max_3prime_complementary_run_bp: suffix(left, &rc)
            .max(suffix(right, &GentleEngine::reverse_complement_bytes(left))),
    }
}

fn validate_request(request: &mut TranscriptCapturePoolRequest) -> Result<(), EngineError> {
    if request.schema != TRANSCRIPT_CAPTURE_REQUEST_SCHEMA {
        return Err(EngineError::invalid_input(
            "Unsupported transcript capture request schema",
        ));
    }
    request.report_id = GentleEngine::normalize_primer_design_report_id(&request.report_id)?;
    let policy = &request.search;
    if request.targets.is_empty()
        || request.targets.len() > 16
        || !(8..=60).contains(&policy.min_length_bp)
        || !(policy.min_length_bp..=60).contains(&policy.max_length_bp)
        || !(1..=32).contains(&policy.max_candidates_per_target)
        || !(1..=512).contains(&policy.beam_width)
        || !(4..=100).contains(&policy.internal_a_run_min_bp)
        || request.fixed_oligos.len() > 64
    {
        return Err(EngineError::invalid_input(
            "Capture limits: 1..=16 targets, oligo lengths 8..=60, 1..=32 retained candidates/target, beam 1..=512, A-run threshold 4..=100, at most 64 fixed oligos",
        ));
    }
    if let Some(tm) = &policy.tm_range
        && (!tm.min_c.is_finite() || !tm.max_c.is_finite() || tm.min_c > tm.max_c)
    {
        return Err(EngineError::invalid_input(
            "Capture tm_range must be finite with min_c <= max_c",
        ));
    }
    let mut ids = BTreeSet::new();
    let mut budget = 0;
    for target in &mut request.targets {
        if target.target_id.trim().is_empty()
            || !ids.insert(target.target_id.clone())
            || target.sources.is_empty()
            || target.sources.len() > 16
            || !(1..=4).contains(&target.max_primers)
        {
            return Err(EngineError::invalid_input(
                "Capture targets need unique nonempty IDs, 1..=16 sources and max_primers 1..=4",
            ));
        }
        budget += target.max_primers;
        normalize_stages(&mut target.stage_ids)?;
        target.tail_5prime = GentleEngine::normalize_iupac_text(&target.tail_5prime)?;
        if !canonical(target.tail_5prime.as_bytes()) || target.tail_5prime.len() > 60 {
            return Err(EngineError::invalid_input(
                "Capture tails must be canonical DNA, at most 60 bases",
            ));
        }
        if target
            .sharing_group
            .as_ref()
            .is_some_and(|key| key.trim().is_empty())
        {
            return Err(EngineError::invalid_input(
                "An explicit sharing_group must not be empty",
            ));
        }
        match target.window {
            TranscriptCaptureWindow::TranscriptRange {
                start_0based,
                end_0based_exclusive,
            } if end_0based_exclusive <= start_0based
                || end_0based_exclusive - start_0based > 5_000 =>
            {
                return Err(EngineError::invalid_input(
                    "Capture transcript ranges must be nonempty and at most 5000 bp",
                ));
            }
            TranscriptCaptureWindow::TerminalExonStart { search_window_bp }
                if !(1..=5_000).contains(&search_window_bp) =>
            {
                return Err(EngineError::invalid_input(
                    "Capture terminal-exon search window must be 1..=5000 bp",
                ));
            }
            _ => {}
        }
    }
    if budget > 16 {
        return Err(EngineError::invalid_input(
            "The total capture primer budget may not exceed 16",
        ));
    }
    let mut fixed_ids = BTreeSet::new();
    for fixed in &mut request.fixed_oligos {
        if fixed.oligo_id.trim().is_empty()
            || !fixed_ids.insert(fixed.oligo_id.clone())
            || fixed.provenance.trim().is_empty()
        {
            return Err(EngineError::invalid_input(
                "Fixed oligos need unique nonempty IDs and explicit provenance",
            ));
        }
        fixed.full_oligo_5_to_3 = GentleEngine::normalize_iupac_text(&fixed.full_oligo_5_to_3)?;
        if fixed.full_oligo_5_to_3.is_empty() || fixed.full_oligo_5_to_3.len() > 120 {
            return Err(EngineError::invalid_input(
                "Fixed complete oligos must have 1..=120 bases",
            ));
        }
        normalize_stages(&mut fixed.stage_ids)?;
    }
    Ok(())
}

impl GentleEngine {
    fn capture_templates(
        &self,
        request: &TranscriptCapturePoolRequest,
    ) -> Result<(Vec<CaptureTemplate>, Vec<TranscriptCaptureSourceBinding>), EngineError> {
        let mut out = Vec::new();
        let mut sources = Vec::new();
        let mut total_bases = 0;
        let mut seen = BTreeSet::new();
        for (target_index, target) in request.targets.iter().enumerate() {
            for source in &target.sources {
                let dna = self.state.sequences.get(&source.seq_id).ok_or_else(|| {
                    EngineError::invalid_input(format!(
                        "Capture source '{}' is not loaded",
                        source.seq_id
                    ))
                })?;
                if dna.is_protein_sequence()
                    || dna
                        .molecule_type()
                        .is_some_and(|kind| kind.to_ascii_lowercase().contains("rna"))
                {
                    return Err(EngineError::invalid_input(
                        "Capture sources must be annotated DNA loci",
                    ));
                }
                let splicing = self.build_splicing_expert_view(
                    &source.seq_id,
                    source.source_feature_id,
                    SplicingScopePreset::TargetGroupTargetStrand,
                )?;
                let all = Self::build_qpcr_transcript_design_templates(dna, &splicing)?;
                let anchor = self.transcript_qpcr_panel_source_anchor(&source.seq_id, dna);
                let assembly = Self::transcript_qpcr_panel_expected_assembly(dna, anchor.as_ref());
                // UniProt inventory-to-local digest joins need separate capture
                // accounting; do not silently downgrade them to local transcripts.
                if source.coverage_universe.kind
                    == TranscriptAssayCoverageUniverseKind::UniprotSupportedIsoforms
                {
                    return Err(EngineError::invalid_input(
                        "Capture discovery currently accepts all_annotated_cdna_classes or explicit_transcripts; resolve a UniProt inventory explicitly first",
                    ));
                }
                let (templates, coverage) = Self::resolve_transcript_assay_coverage_universe(
                    source.coverage_universe.clone(),
                    &source.seq_id,
                    &all,
                    true,
                    assembly.as_deref(),
                    source.annotation_release.as_deref(),
                )?;
                if templates.is_empty()
                    || !coverage.unresolved_target_ids.is_empty()
                    || !coverage.ambiguous_target_ids.is_empty()
                {
                    return Err(EngineError::invalid_input(format!(
                        "Capture universe is empty, unresolved or ambiguous: {}",
                        serde_json::to_string(&coverage).unwrap_or_default()
                    )));
                }
                sources.push(TranscriptCaptureSourceBinding {
                    target_id: target.target_id.clone(),
                    source: source.clone(),
                    genome_anchor: anchor,
                    source_record_sha256: sha256_prefixed_bytes(
                        serde_json::to_value(dna)
                            .map_err(|e| EngineError::internal(e.to_string()))?
                            .to_string()
                            .as_bytes(),
                    ),
                    coverage,
                });
                for mut template in templates {
                    if !seen.insert((
                        target_index,
                        source.seq_id.clone(),
                        template.transcript_id.clone(),
                    )) {
                        return Err(EngineError::invalid_input(
                            "A transcript occurs twice in the same capture target",
                        ));
                    }
                    template.sequence = Self::normalize_iupac_text(&template.sequence)?;
                    total_bases += template.sequence.len();
                    if out.len() >= 256 || total_bases > 1_000_000 {
                        return Err(EngineError::invalid_input(
                            "Capture discovery is bounded to 256 transcripts and 1,000,000 total cDNA bases",
                        ));
                    }
                    let feature = dna
                        .features()
                        .get(template.transcript_feature_id)
                        .ok_or_else(|| {
                            EngineError::internal("Capture transcript feature is absent")
                        })?;
                    let cds = Self::resolve_transcript_source_cds_ranges_0based(
                        feature,
                        dna.features(),
                        &template.exon_chain,
                    );
                    let mut exons = template.exon_chain.clone();
                    exons.sort_unstable();
                    let mut cursor = 0;
                    let forward = exons
                        .iter()
                        .map(|&(start, end)| {
                            let span = (start, end, cursor, cursor + end - start);
                            cursor += end - start;
                            span
                        })
                        .collect::<Vec<_>>();
                    if cursor != template.sequence.len() {
                        return Err(EngineError::invalid_input(
                            "Capture source has clipped/inconsistent exon geometry; retrieve the complete annotated transcript locus",
                        ));
                    }
                    let local_cds = Self::map_source_ranges_to_transcript_local_ranges_0based(
                        &cds,
                        &forward,
                        template.strand == "-",
                        cursor,
                    );
                    let utr_length = local_cds.iter().map(|range| range.0).min();
                    let mut notes = Vec::new();
                    if utr_length.is_none() {
                        notes.push("No resolvable annotated CDS boundary; 5-prime UTR is unknown, not assumed absent.".into());
                    }
                    let range = match target.window {
                        TranscriptCaptureWindow::FivePrimeUtr => utr_length.map(|end| (0, end)),
                        TranscriptCaptureWindow::TranscriptRange {
                            start_0based,
                            end_0based_exclusive,
                        } => {
                            if end_0based_exclusive > cursor {
                                notes.push("Requested transcript range extends beyond this member; no clipped fallback used.".into());
                                None
                            } else {
                                Some((start_0based, end_0based_exclusive))
                            }
                        }
                        TranscriptCaptureWindow::TerminalExonStart { search_window_bp } => {
                            template.local_exon_segments.last().map(|exon| {
                                (
                                    exon.local_start_0based,
                                    exon.local_end_0based_exclusive
                                        .min(exon.local_start_0based + search_window_bp),
                                )
                            })
                        }
                    };
                    if range.is_some_and(|(start, end)| end - start > 5_000) {
                        return Err(EngineError::invalid_input(
                            "A capture UTR exceeds 5000 bp; provide an explicit bounded transcript_range",
                        ));
                    }
                    if range.is_none_or(|(start, end)| end - start < request.search.min_length_bp) {
                        notes.push("No permitted window long enough for a candidate; member remains in the coverage denominator.".into());
                    }
                    let member = TranscriptCaptureMember {
                        member_id: short_sha256_id(
                            "capture_member",
                            &format!(
                                "{}|{}|{}",
                                target.target_id, source.seq_id, template.transcript_id
                            ),
                        ),
                        target_id: target.target_id.clone(),
                        seq_id: source.seq_id.clone(),
                        transcript_id: template.transcript_id.clone(),
                        transcript_feature_id: template.transcript_feature_id,
                        source_strand: template.strand.clone(),
                        cdna_length_bp: cursor,
                        cdna_sha256: sha256_prefixed_bytes(template.sequence.as_bytes()),
                        five_prime_utr_length_bp: utr_length,
                        search_range_0based: range,
                        notes,
                    };
                    out.push(CaptureTemplate {
                        target_index,
                        member,
                        template,
                    });
                }
            }
        }
        out.sort_by(|a, b| a.member.member_id.cmp(&b.member.member_id));
        Ok((out, sources))
    }

    pub(super) fn design_transcript_capture_pool(
        &mut self,
        mut request: TranscriptCapturePoolRequest,
        op_id: &str,
        run_id: &str,
    ) -> Result<TranscriptCapturePoolReport, EngineError> {
        validate_request(&mut request)?;
        let previous = self.read_primer_design_store();
        if previous.reports.contains_key(&request.report_id)
            || previous.qpcr_reports.contains_key(&request.report_id)
            || previous
                .terminal_exon_rt_primer_pools
                .contains_key(&request.report_id)
            || previous
                .primer_specificity_reports
                .contains_key(&request.report_id)
            || previous
                .transcript_assay_panels
                .contains_key(&request.report_id)
        {
            return Err(EngineError::invalid_input(
                "Capture report_id is already used by another primer-report family",
            ));
        }
        let (templates, sources) = self.capture_templates(&request)?;
        let (mut candidates, ambiguous_windows_skipped) = discover(&request, &templates)?;
        let distinct_candidates_evaluated = candidates.len();
        let mut retained_ids = BTreeSet::new();
        for target in &request.targets {
            let target_members = templates
                .iter()
                .filter(|template| template.member.target_id == target.target_id)
                .map(|template| &template.member.member_id)
                .collect::<BTreeSet<_>>();
            let mut eligible = candidates
                .iter()
                .filter(|candidate| {
                    candidate.permitted_target_ids.contains(&target.target_id)
                        && candidate
                            .covered_member_ids
                            .iter()
                            .any(|member_id| target_members.contains(member_id))
                })
                .collect::<Vec<_>>();
            eligible.sort_by_key(|candidate| candidate_rank(candidate));
            retained_ids.extend(
                eligible
                    .into_iter()
                    .take(request.search.max_candidates_per_target)
                    .map(|candidate| candidate.candidate_id.clone()),
            );
        }
        let candidate_retention_truncated = retained_ids.len() < candidates.len();
        candidates.retain(|candidate| retained_ids.contains(&candidate.candidate_id));
        if candidates.len() > 128 {
            return Err(EngineError::invalid_input(
                "Capture pool has more than 128 retained candidates; reduce max_candidates_per_target",
            ));
        }
        // Only retained candidates need the full transcript occurrence audit.
        // Do not limit it to eligible windows: repeats outside them also matter.
        let mut binding_count = 0;
        let mut retained_bases_audited = 0;
        for candidate in &mut candidates {
            candidate.bindings = candidate_bindings(
                candidate,
                &templates,
                request.search.internal_a_run_min_bp,
                &mut binding_count,
                &mut retained_bases_audited,
            )?;
        }
        let (indices, pool_states_evaluated, beam_truncated) =
            select_pool(&request, &templates, &candidates);
        let covered = indices
            .iter()
            .flat_map(|&i| candidates[i].covered_member_ids.iter().cloned())
            .collect::<BTreeSet<_>>();
        let uncovered_member_ids = templates
            .iter()
            .filter(|t| !covered.contains(&t.member.member_id))
            .map(|t| t.member.member_id.clone())
            .collect::<Vec<_>>();
        let proposed_candidate_ids = indices
            .iter()
            .map(|&i| candidates[i].candidate_id.clone())
            .collect();
        let mut interactions = Vec::new();
        let mut oligos = indices
            .iter()
            .map(|&i| {
                (
                    format!("candidate:{}", candidates[i].candidate_id),
                    &candidates[i].full_oligo_5_to_3,
                    &candidates[i].stage_ids,
                )
            })
            .collect::<Vec<_>>();
        oligos.extend(request.fixed_oligos.iter().map(|o| {
            (
                format!("fixed:{}", o.oligo_id),
                &o.full_oligo_5_to_3,
                &o.stage_ids,
            )
        }));
        for left in 0..oligos.len() {
            for right in left + 1..oligos.len() {
                let stages = shared_stages(oligos[left].2, oligos[right].2);
                if stages.is_empty() {
                    continue;
                }
                let metrics = dimer(oligos[left].1.as_bytes(), oligos[right].1.as_bytes());
                interactions.push(TranscriptCaptureInteraction {
                    left_oligo_id: oligos[left].0.clone(),
                    right_oligo_id: oligos[right].0.clone(),
                    shared_stage_ids: stages,
                    max_complementary_run_bp: metrics.max_complementary_run_bp,
                    max_3prime_complementary_run_bp: metrics.max_3prime_complementary_run_bp,
                });
            }
        }
        let mut groups = BTreeMap::<(String, usize), Vec<String>>::new();
        let mut lengths = Vec::new();
        for &index in &indices {
            for binding in &candidates[index].bindings {
                if !binding.within_requested_window || !binding.permitted_target {
                    continue;
                }
                lengths.push(binding.retained_length_bp);
                if binding.retained_sequence_canonical {
                    groups
                        .entry((
                            binding.retained_sequence_sha256.clone(),
                            binding.retained_length_bp,
                        ))
                        .or_default()
                        .push(format!(
                            "{}:{}:{}",
                            candidates[index].candidate_id,
                            binding.member_id,
                            binding.transcript_start_0based
                        ));
                }
            }
        }
        let mut warnings = vec![
            "Candidate discovery only: whole-genome and whole-transcriptome specificity, actual RT reach, amplification yield, library adapters and order readiness are unassessed.".into(),
            "Exact coverage is not full-length recovery. Primer-supplied bases are not independent observations of the original template; Nanopore read orientation is not inferred.".into(),
            "Retained intervals use annotated transcript ends, not a measured poly(A) tail or a guaranteed oligo(dT) endpoint. Equivalence excludes synthetic tails and is reported only for canonical sequences.".into(),
            "Internal exact A-runs are possible priming sites, not measured truncated products. Other A-rich patterns, RT falloff and length-dependent amplification bias are not predicted.".into(),
            "Fixed oligos affect joint stage-specific interaction scoring and retain reorder provenance; their target coverage is not inferred. IUPAC complementarity is a conservative possible-match heuristic, not a binding probability.".into(),
            "Bounded candidate retention/beam search reports the best found pool, not proof of minimum primer count or biological feasibility.".into(),
        ];
        if !uncovered_member_ids.is_empty() {
            warnings.push(format!("{} requested transcript members have no eligible site in the proposed pool; do not interpret this as a complete {} design.", uncovered_member_ids.len(), request.coverage_policy.as_str()));
        }
        if request.fixed_oligos.is_empty() {
            warnings.push("No retained/protocol oligos were supplied; interaction assessment is incomplete for the actual reaction.".into());
        }
        if request.cdna_synthesis == TranscriptAssayCdnaSynthesis::Unspecified {
            warnings.push("cDNA synthesis is unspecified; capture geometry alone does not establish a usable protocol.".into());
        }
        let report = TranscriptCapturePoolReport {
            schema: TRANSCRIPT_CAPTURE_REPORT_SCHEMA.into(),
            report_id: request.report_id.clone(),
            op_id: op_id.into(),
            run_id: run_id.into(),
            request_sha256: sha256_prefixed_bytes(
                &serde_json::to_vec(&request).map_err(|e| EngineError::internal(e.to_string()))?,
            ),
            request,
            sources,
            members: templates.into_iter().map(|t| t.member).collect(),
            candidates,
            proposed_candidate_ids,
            coverage_satisfied: uncovered_member_ids.is_empty(),
            uncovered_member_ids,
            distinct_candidates_evaluated,
            ambiguous_windows_skipped,
            candidate_retention_truncated,
            beam_truncated,
            pool_states_evaluated,
            interactions,
            captured_equivalence_groups: groups
                .into_iter()
                .map(|((hash, len), instances)| TranscriptCaptureEquivalence {
                    retained_sequence_sha256: hash,
                    retained_length_bp: len,
                    binding_instances: instances,
                })
                .collect(),
            retained_length_range_bp: lengths
                .iter()
                .min()
                .zip(lengths.iter().max())
                .map(|(&min, &max)| (min, max)),
            tm_model: Self::primer_tm_model_description(),
            specificity_status: "unassessed".into(),
            warnings,
        };
        let mut store = self.read_primer_design_store();
        store
            .transcript_capture_pools
            .insert(report.report_id.clone(), report.clone());
        self.write_primer_design_store(store)?;
        Ok(report)
    }

    /// Retrieve a saved discovery without recomputing or upgrading readiness.
    pub fn get_transcript_capture_pool_report(
        &self,
        report_id: &str,
    ) -> Result<TranscriptCapturePoolReport, EngineError> {
        self.read_primer_design_store()
            .transcript_capture_pools
            .get(report_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Transcript capture report '{report_id}' not found"),
                cause_chain: vec![],
            })
    }

    /// Stable report IDs for the shared primer-report listing.
    pub fn list_transcript_capture_pool_report_ids(&self) -> Vec<String> {
        let mut ids = self
            .read_primer_design_store()
            .transcript_capture_pools
            .into_keys()
            .collect::<Vec<_>>();
        ids.sort();
        ids
    }

    /// Export the stored report verbatim; this is not a new scientific screen.
    pub fn export_transcript_capture_pool_report(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<TranscriptCapturePoolReport, EngineError> {
        let report = self.get_transcript_capture_pool_report(report_id)?;
        let bytes =
            serde_json::to_vec_pretty(&report).map_err(|e| EngineError::internal(e.to_string()))?;
        std::fs::write(path, bytes).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not export capture report to '{path}': {e}"),
            cause_chain: vec![],
        })?;
        Ok(report)
    }
}

fn candidate_rank(
    candidate: &TranscriptCaptureCandidate,
) -> (std::cmp::Reverse<usize>, usize, usize, usize, usize, String) {
    (
        std::cmp::Reverse(candidate.covered_member_ids.len()),
        candidate
            .self_3prime_run_bp
            .max(candidate.max_fixed_3prime_run_bp),
        candidate
            .self_complementary_run_bp
            .max(candidate.max_fixed_complementary_run_bp),
        candidate.homopolymer_run_bp,
        candidate.annealing_5_to_3.len(),
        candidate.candidate_id.clone(),
    )
}

fn share_allowed(left: &TranscriptCaptureTarget, right: &TranscriptCaptureTarget) -> bool {
    left.target_id == right.target_id
        || (left.sharing_group.is_some()
            && left.sharing_group == right.sharing_group
            && left.role == right.role
            && left.tail_5prime == right.tail_5prime
            && left.stage_ids == right.stage_ids)
}

fn discover(
    request: &TranscriptCapturePoolRequest,
    templates: &[CaptureTemplate],
) -> Result<(Vec<TranscriptCaptureCandidate>, usize), EngineError> {
    let mut candidates = BTreeMap::<String, TranscriptCaptureCandidate>::new();
    let mut skipped = 0;
    for source in templates {
        let Some((start, end)) = source.member.search_range_0based else {
            continue;
        };
        let target = &request.targets[source.target_index];
        let permitted = request
            .targets
            .iter()
            .filter(|other| share_allowed(target, other))
            .map(|t| t.target_id.clone())
            .collect::<Vec<_>>();
        for length in request.search.min_length_bp..=request.search.max_length_bp {
            if end - start < length {
                continue;
            }
            for position in start..=end - length {
                let word = &source.template.sequence.as_bytes()[position..position + length];
                if !canonical(word) {
                    skipped += 1;
                    continue;
                }
                let annealing = match target.role {
                    TranscriptCaptureRole::SenseForward => {
                        String::from_utf8(word.to_vec()).unwrap()
                    }
                    TranscriptCaptureRole::AntisenseReverse => {
                        GentleEngine::reverse_complement(std::str::from_utf8(word).unwrap())
                    }
                };
                let full = format!("{}{}", target.tail_5prime, annealing);
                let key =
                    serde_json::to_string(&(&permitted, target.role, &target.stage_ids, &full))
                        .map_err(|e| EngineError::internal(e.to_string()))?;
                if !candidates.contains_key(&key) {
                    let tm = GentleEngine::estimate_primer_tm_c(annealing.as_bytes());
                    if request
                        .search
                        .tm_range
                        .as_ref()
                        .is_some_and(|range| tm < range.min_c || tm > range.max_c)
                    {
                        continue;
                    }
                    let metrics =
                        GentleEngine::compute_primer_heuristic_metrics(annealing.as_bytes());
                    let self_metrics = dimer(full.as_bytes(), full.as_bytes());
                    let fixed = request
                        .fixed_oligos
                        .iter()
                        .filter(|o| !shared_stages(&target.stage_ids, &o.stage_ids).is_empty())
                        .map(|o| dimer(full.as_bytes(), o.full_oligo_5_to_3.as_bytes()))
                        .collect::<Vec<_>>();
                    candidates.insert(
                        key.clone(),
                        TranscriptCaptureCandidate {
                            candidate_id: short_sha256_id("capture_oligo", &key),
                            role: target.role,
                            permitted_target_ids: permitted.clone(),
                            stage_ids: target.stage_ids.clone(),
                            annealing_5_to_3: annealing,
                            full_oligo_5_to_3: full,
                            tm_c: tm,
                            self_3prime_run_bp: self_metrics.max_3prime_complementary_run_bp,
                            self_complementary_run_bp: self_metrics.max_complementary_run_bp,
                            homopolymer_run_bp: metrics.longest_homopolymer_run_bp,
                            max_fixed_3prime_run_bp: fixed
                                .iter()
                                .map(|m| m.max_3prime_complementary_run_bp)
                                .max()
                                .unwrap_or(0),
                            max_fixed_complementary_run_bp: fixed
                                .iter()
                                .map(|m| m.max_complementary_run_bp)
                                .max()
                                .unwrap_or(0),
                            covered_member_ids: Vec::new(),
                            bindings: Vec::new(),
                        },
                    );
                    if candidates.len() > 20_000 {
                        return Err(EngineError::invalid_input(
                            "Capture search exceeded 20000 distinct candidates; narrow the windows or length range",
                        ));
                    }
                }
                let candidate = candidates.get_mut(&key).unwrap();
                if !candidate
                    .covered_member_ids
                    .contains(&source.member.member_id)
                {
                    candidate
                        .covered_member_ids
                        .push(source.member.member_id.clone());
                }
            }
        }
    }
    for candidate in candidates.values_mut() {
        candidate.covered_member_ids.sort();
    }
    Ok((candidates.into_values().collect(), skipped))
}

fn a_runs(sequence: &[u8], minimum: usize, offset: usize) -> Vec<(usize, usize)> {
    let mut ranges = Vec::new();
    let mut start = 0;
    while start < sequence.len() {
        if sequence[start] != b'A' {
            start += 1;
            continue;
        }
        let mut end = start + 1;
        while end < sequence.len() && sequence[end] == b'A' {
            end += 1;
        }
        if end - start >= minimum {
            ranges.push((start + offset, end + offset));
        }
        start = end;
    }
    ranges
}

fn candidate_bindings(
    candidate: &TranscriptCaptureCandidate,
    templates: &[CaptureTemplate],
    a_min: usize,
    binding_count: &mut usize,
    retained_bases_audited: &mut usize,
) -> Result<Vec<TranscriptCaptureBinding>, EngineError> {
    let word = match candidate.role {
        TranscriptCaptureRole::SenseForward => candidate.annealing_5_to_3.clone(),
        TranscriptCaptureRole::AntisenseReverse => {
            GentleEngine::reverse_complement(&candidate.annealing_5_to_3)
        }
    };
    let mut bindings = Vec::new();
    for source in templates {
        let sequence = source.template.sequence.as_bytes();
        for (start, window) in sequence.windows(word.len()).enumerate() {
            if window != word.as_bytes() {
                continue;
            }
            let end = start + word.len();
            let (retained_start, retained_end) = match candidate.role {
                TranscriptCaptureRole::SenseForward => (start, sequence.len()),
                TranscriptCaptureRole::AntisenseReverse => (0, end),
            };
            *binding_count += 1;
            *retained_bases_audited += retained_end - retained_start;
            if *binding_count > 50_000 || *retained_bases_audited > 100_000_000 {
                return Err(EngineError::invalid_input(
                    "Capture occurrence audit exceeded 50000 bindings or 100 million retained bases; narrow the request",
                ));
            }
            let mapped = GentleEngine::map_transcript_local_interval(
                &source.template.local_exon_segments,
                source.template.strand == "-",
                start,
                end,
            );
            bindings.push(TranscriptCaptureBinding {
                member_id: source.member.member_id.clone(),
                transcript_start_0based: start,
                transcript_end_0based_exclusive: end,
                source_ranges_0based: mapped
                    .source_ranges_0based
                    .into_iter()
                    .map(|r| (r.start_0based, r.end_0based_exclusive))
                    .collect(),
                within_requested_window: source
                    .member
                    .search_range_0based
                    .is_some_and(|(a, b)| a <= start && end <= b),
                permitted_target: candidate
                    .permitted_target_ids
                    .contains(&source.member.target_id),
                upstream_bases_omitted: retained_start,
                downstream_bases_omitted: sequence.len() - retained_end,
                retained_length_bp: retained_end - retained_start,
                retained_sequence_sha256: sha256_prefixed_bytes(
                    &sequence[retained_start..retained_end],
                ),
                retained_sequence_canonical: canonical(&sequence[retained_start..retained_end]),
                internal_a_runs_0based: if candidate.role == TranscriptCaptureRole::SenseForward {
                    a_runs(&sequence[end..], a_min, end)
                } else {
                    vec![]
                },
            });
        }
    }
    Ok(bindings)
}

#[derive(Clone)]
struct PoolState {
    indices: Vec<usize>,
    covered: BTreeSet<String>,
    target_counts: Vec<usize>,
    worst_3prime: usize,
    worst_general: usize,
    worst_homopolymer: usize,
    omitted_bases: usize,
}

fn pool_rank(
    state: &PoolState,
    count: usize,
) -> (usize, usize, usize, usize, usize, usize, Vec<usize>) {
    (
        count - state.covered.len(),
        state.indices.len(),
        state.worst_3prime,
        state.worst_general,
        state.worst_homopolymer,
        state.omitted_bases,
        state.indices.clone(),
    )
}

fn select_pool(
    request: &TranscriptCapturePoolRequest,
    templates: &[CaptureTemplate],
    candidates: &[TranscriptCaptureCandidate],
) -> (Vec<usize>, usize, bool) {
    let mut pair_metrics = BTreeMap::new();
    for (i, left) in candidates.iter().enumerate() {
        for (j, right) in candidates.iter().enumerate().skip(i + 1) {
            if !shared_stages(&left.stage_ids, &right.stage_ids).is_empty() {
                pair_metrics.insert(
                    (i, j),
                    dimer(
                        left.full_oligo_5_to_3.as_bytes(),
                        right.full_oligo_5_to_3.as_bytes(),
                    ),
                );
            }
        }
    }
    let memberships = candidates
        .iter()
        .map(|c| {
            request
                .targets
                .iter()
                .map(|t| {
                    templates.iter().any(|member| {
                        member.member.target_id == t.target_id
                            && c.covered_member_ids.contains(&member.member.member_id)
                    })
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let initial = PoolState {
        indices: vec![],
        covered: BTreeSet::new(),
        target_counts: vec![0; request.targets.len()],
        worst_3prime: 0,
        worst_general: 0,
        worst_homopolymer: 0,
        omitted_bases: 0,
    };
    let mut best = initial.clone();
    let mut states = vec![initial];
    let mut evaluated = 0;
    let mut truncated = false;
    for _ in 0..request.targets.iter().map(|t| t.max_primers).sum::<usize>() {
        let mut expanded = BTreeMap::<Vec<usize>, PoolState>::new();
        for state in &states {
            for (index, candidate) in candidates.iter().enumerate() {
                if state.indices.contains(&index)
                    || candidate
                        .covered_member_ids
                        .iter()
                        .all(|id| state.covered.contains(id))
                {
                    continue;
                }
                if memberships[index].iter().enumerate().any(|(i, &included)| {
                    included && state.target_counts[i] >= request.targets[i].max_primers
                }) {
                    continue;
                }
                let mut next = state.clone();
                for (i, &included) in memberships[index].iter().enumerate() {
                    if included {
                        next.target_counts[i] += 1;
                    }
                }
                next.indices.push(index);
                next.indices.sort_unstable();
                next.covered
                    .extend(candidate.covered_member_ids.iter().cloned());
                next.worst_3prime = next
                    .worst_3prime
                    .max(candidate.self_3prime_run_bp)
                    .max(candidate.max_fixed_3prime_run_bp);
                next.worst_general = next
                    .worst_general
                    .max(candidate.self_complementary_run_bp)
                    .max(candidate.max_fixed_complementary_run_bp);
                next.worst_homopolymer = next.worst_homopolymer.max(candidate.homopolymer_run_bp);
                next.omitted_bases += candidate
                    .bindings
                    .iter()
                    .filter(|b| b.permitted_target && b.within_requested_window)
                    .map(|b| b.upstream_bases_omitted + b.downstream_bases_omitted)
                    .sum::<usize>();
                for &other in &state.indices {
                    if let Some(metrics) = pair_metrics.get(&(index.min(other), index.max(other))) {
                        next.worst_3prime = next
                            .worst_3prime
                            .max(metrics.max_3prime_complementary_run_bp);
                        next.worst_general =
                            next.worst_general.max(metrics.max_complementary_run_bp);
                    }
                }
                evaluated += 1;
                expanded.entry(next.indices.clone()).or_insert(next);
            }
        }
        if expanded.is_empty() {
            break;
        }
        let mut ranked = expanded.into_values().collect::<Vec<_>>();
        ranked.sort_by_key(|state| pool_rank(state, templates.len()));
        if pool_rank(&ranked[0], templates.len()) < pool_rank(&best, templates.len()) {
            best = ranked[0].clone();
        }
        if best.covered.len() == templates.len() {
            break;
        }
        truncated |= ranked.len() > request.search.beam_width;
        ranked.truncate(request.search.beam_width);
        states = ranked;
    }
    (best.indices, evaluated, truncated)
}

#[cfg(test)]
#[path = "transcript_capture_tests.rs"]
mod tests;
