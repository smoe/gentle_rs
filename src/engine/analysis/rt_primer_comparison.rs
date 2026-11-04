//! Shared cDNA/CDS geometry and prospective evidence composed into an advisory RT comparison.

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
        Err(invalid("rt_comparison_cancelled"))
    }
}

/// Bound annotation expansion and the shared view's exon-pair/CDS projection work before invoking it.
fn admit_annotation(
    dna: &DNAsequence,
    nodes: &mut usize,
    text_bytes: &mut usize,
    geometry: &mut u64,
) -> Result<(), EngineError> {
    use gb_io::seq::Location;
    let mut intervals = BTreeSet::new();
    let mut seeds = 0u64;
    for feature in dna.features() {
        let seed = GentleEngine::is_splicing_seed_feature(feature);
        seeds += u64::from(seed);
        let mut stack = vec![(&feature.location, 0usize)];
        while let Some((location, depth)) = stack.pop() {
            *nodes += 1;
            if *nodes > 16_384 || depth > 32 {
                return Err(invalid("rt_comparison_annotation_complexity_budget"));
            }
            match location {
                Location::Range((a, _), (z, _)) => {
                    if seed {
                        intervals.insert((*a, *z));
                    }
                }
                Location::Complement(child) => stack.push((child, depth + 1)),
                Location::Join(children)
                | Location::Order(children)
                | Location::Bond(children)
                | Location::OneOf(children) => {
                    if children.len() > 16_384 {
                        return Err(invalid("rt_comparison_annotation_complexity_budget"));
                    }
                    stack.extend(children.iter().map(|c| (c, depth + 1)));
                }
                Location::External(id, child) => {
                    *text_bytes = text_bytes
                        .checked_add(id.len())
                        .filter(|n| *n <= 16_000_000)
                        .ok_or_else(|| invalid("rt_comparison_annotation_text_budget"))?;
                    if let Some(child) = child {
                        stack.push((child, depth + 1));
                    }
                }
                _ => {}
            }
        }
        for (key, value) in &feature.qualifiers {
            *text_bytes = text_bytes
                .checked_add(key.len())
                .and_then(|n| n.checked_add(value.as_ref().map_or(0, |v| v.len())))
                .filter(|n| *n <= 16_000_000)
                .ok_or_else(|| invalid("rt_comparison_annotation_text_budget"))?;
            if seed && matches!(key.as_ref(), "exon_ranges_1based" | "cds_ranges_1based") {
                if let Some(value) = value {
                    let tokens = value.split(',').filter(|s| !s.trim().is_empty()).count();
                    *nodes = nodes
                        .checked_add(tokens)
                        .filter(|n| *n <= 16_384)
                        .ok_or_else(|| invalid("rt_comparison_annotation_complexity_budget"))?;
                    for (a, z) in GentleEngine::parse_ranges_1based(value) {
                        intervals.insert(((a - 1) as i64, z as i64));
                    }
                }
            }
        }
    }
    let n = intervals.len() as u64;
    let work = n
        .checked_mul(n)
        .and_then(|n| n.checked_mul(seeds.max(1)))
        .ok_or_else(|| invalid("rt_comparison_geometry_work_budget"))?;
    *geometry = geometry
        .checked_add(work)
        .filter(|n| *n <= 25_000_000)
        .ok_or_else(|| invalid("rt_comparison_geometry_work_budget"))?;
    Ok(())
}

pub(super) fn metrics(s: &[u8]) -> RtPrimerComparisonMetrics {
    let m = GentleEngine::compute_primer_heuristic_metrics(s);
    RtPrimerComparisonMetrics {
        length_bp: s.len(),
        gc_fraction: GentleEngine::sequence_gc_fraction(s).unwrap_or(0.0),
        tm_c: GentleEngine::estimate_primer_tm_c(s),
        self_complementary_run_bp: m.self_complementary_run_bp,
        self_3prime_complementary_run_bp:
            GentleEngine::compute_primer_self_3prime_complementary_run(s),
        longest_homopolymer_run_bp: m.longest_homopolymer_run_bp,
    }
}

pub(super) fn classify(
    start: usize,
    end: usize,
    cds: &[(usize, usize)],
) -> RtPrimerClassifiedBinding {
    let overlap = |a: usize, z: usize| end.min(z).saturating_sub(start.max(a));
    let first = cds.partition_point(|r| r.1 <= start);
    let coding = cds[first..]
        .iter()
        .take_while(|r| r.0 < end)
        .map(|&(a, z)| overlap(a, z))
        .sum::<usize>();
    let five = cds.first().map_or(0, |r| overlap(0, r.0));
    let three = cds.last().map_or(0, |r| end.saturating_sub(start.max(r.1)));
    let unknown = end - start - coding - five - three;
    let region = if unknown > 0 {
        RtPrimerBindingRegion::Unknown
    } else if coding == end - start {
        RtPrimerBindingRegion::Cds
    } else if five == end - start {
        RtPrimerBindingRegion::FivePrimeUtr
    } else if three == end - start {
        RtPrimerBindingRegion::ThreePrimeUtr
    } else {
        RtPrimerBindingRegion::Mixed
    };
    RtPrimerClassifiedBinding {
        transcript_start_0based: start,
        transcript_end_0based_exclusive: end,
        source_ranges_0based: Vec::new(),
        region,
        cds_bases: coding,
        five_prime_utr_bases: five,
        three_prime_utr_bases: three,
        unknown_bases: unknown,
    }
}

impl GentleEngine {
    fn rt_comparison_templates(
        &self,
        request: &CompareRtPrimerCandidatesRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<
        (
            Vec<TranscriptQpcrDesignTemplate>,
            Vec<RtPrimerComparisonMember>,
            Vec<TranscriptCaptureSourceBinding>,
        ),
        EngineError,
    > {
        let mut templates = Vec::new();
        let mut members = Vec::new();
        let mut sources = Vec::new();
        let mut seen = BTreeSet::new();
        let (mut materialized_members, mut materialized_bases) = (0usize, 0usize);
        for (index, binding) in request.sources.iter().enumerate() {
            poll(progress, index, request.sources.len())?;
            let source = &binding.source;
            let dna = self
                .state
                .sequences
                .get(&source.seq_id)
                .ok_or_else(|| invalid("rt_comparison_source_not_loaded"))?;
            let record_hash = digest(dna)?;
            if binding
                .expected_source_record_sha256
                .as_ref()
                .is_some_and(|h| h != &record_hash)
            {
                return Err(invalid("rt_comparison_source_record_changed"));
            }
            let view = self.build_splicing_expert_view(
                &source.seq_id,
                source.source_feature_id,
                SplicingScopePreset::TargetGroupTargetStrand,
            )?;
            materialized_members = materialized_members
                .checked_add(view.transcripts.len())
                .filter(|v| *v <= 256)
                .ok_or_else(|| invalid("rt_comparison_transcript_count_budget"))?;
            for lane in &view.transcripts {
                if !matches!(lane.strand.as_str(), "+" | "-") {
                    return Err(invalid("rt_comparison_transcript_strand_unknown"));
                }
                for exon in &lane.exons {
                    let n = exon
                        .end_1based
                        .checked_sub(exon.start_1based)
                        .and_then(|n| n.checked_add(1))
                        .ok_or_else(|| invalid("rt_comparison_exon_geometry"))?;
                    materialized_bases = materialized_bases
                        .checked_add(n)
                        .filter(|v| *v <= 1_000_000)
                        .ok_or_else(|| invalid("rt_comparison_cdna_base_budget"))?;
                }
            }
            let all = Self::build_qpcr_transcript_design_templates(dna, &view)?;
            let anchor = self.transcript_qpcr_panel_source_anchor(&source.seq_id, dna);
            let assembly = Self::transcript_qpcr_panel_expected_assembly(dna, anchor.as_ref());
            if source.coverage_universe.kind
                == TranscriptAssayCoverageUniverseKind::UniprotSupportedIsoforms
            {
                return Err(invalid(
                    "rt_comparison_resolve_uniprot_inventory_explicitly_first",
                ));
            }
            let (selected, coverage) = Self::resolve_transcript_assay_coverage_universe(
                source.coverage_universe.clone(),
                &source.seq_id,
                &all,
                true,
                assembly.as_deref(),
                source.annotation_release.as_deref(),
            )?;
            if selected.is_empty()
                || !coverage.unresolved_target_ids.is_empty()
                || !coverage.ambiguous_target_ids.is_empty()
            {
                return Err(invalid(
                    "rt_comparison_unresolved_or_empty_transcript_universe",
                ));
            }
            sources.push(TranscriptCaptureSourceBinding {
                target_id: "rt_comparison".into(),
                source: source.clone(),
                source_record_sha256: record_hash,
                genome_anchor: anchor,
                coverage,
            });
            for mut template in selected {
                if !seen.insert((source.seq_id.clone(), template.transcript_id.clone())) {
                    return Err(invalid("rt_comparison_duplicate_transcript_member"));
                }
                template.sequence = Self::normalize_iupac_text(&template.sequence)?;
                if template.sequence.is_empty() {
                    return Err(invalid("rt_comparison_empty_transcript_geometry"));
                }
                let feature = &dna.features()[template.transcript_feature_id];
                let explicit = feature
                    .qualifier_values("cds_ranges_1based")
                    .collect::<Vec<_>>();
                if explicit.len() > 1 {
                    return Err(invalid("rt_comparison_multiple_explicit_cds_qualifiers"));
                }
                for raw in explicit {
                    let mut valid = 0usize;
                    let mut invalid_tokens = 0usize;
                    for token in raw.split(',').filter(|s| !s.trim().is_empty()) {
                        if Self::parse_ranges_1based(token).is_empty() {
                            invalid_tokens += 1;
                        } else {
                            valid += 1;
                        }
                        if valid + invalid_tokens > 4096 {
                            return Err(invalid("rt_comparison_cds_token_budget"));
                        }
                    }
                    if valid > 0 && invalid_tokens > 0 {
                        return Err(invalid("rt_comparison_partially_malformed_cds"));
                    }
                }
                let cds =
                    Self::transcript_design_local_cds_ranges(dna, &template, "RT comparison")?;
                let mut notes = Vec::new();
                if cds.is_empty() {
                    notes.push("No resolvable CDS annotation; noncoding and unknown CDS are not an all-UTR call.".into());
                }
                if cds.windows(2).any(|w| w[0].1 < w[1].0) {
                    notes.push("Gaps inside the annotated CDS union remain unknown.".into());
                }
                members.push(RtPrimerComparisonMember {
                    member_id: short_sha256_id(
                        "rt_member",
                        &serde_json::to_string(&(&source.seq_id, &template.transcript_id))
                            .map_err(|e| EngineError::internal(e.to_string()))?,
                    ),
                    source_index: index,
                    seq_id: source.seq_id.clone(),
                    transcript_id: template.transcript_id.clone(),
                    transcript_feature_id: template.transcript_feature_id,
                    strand: template.strand.clone(),
                    cdna_length_bp: template.sequence.len(),
                    cdna_sha256: sha256_prefixed_bytes(template.sequence.as_bytes()),
                    cds_ranges_0based: cds,
                    notes,
                });
                templates.push(template);
            }
        }
        Ok((templates, members, sources))
    }

    /// Shared bounded input preparation for comparison and constrained redesign.
    pub(super) fn prepare_rt_comparison(
        &self,
        request: &CompareRtPrimerCandidatesRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<
        (
            Vec<TranscriptQpcrDesignTemplate>,
            Vec<RtPrimerComparisonMember>,
            Vec<TranscriptCaptureSourceBinding>,
        ),
        EngineError,
    > {
        if request.schema != RT_PRIMER_COMPARISON_REQUEST_SCHEMA
            || !(1..=16).contains(&request.sources.len())
            || !(1..=1_000_000_000).contains(&request.budgets.max_transcript_base_comparisons)
            || !(1..=50_000).contains(&request.budgets.max_bindings)
        {
            return Err(invalid("rt_comparison_schema_or_budget"));
        }
        poll(progress, 0, request.sources.len())?;
        let (mut source_bases, mut source_features) = (0usize, 0usize);
        let (mut annotation_nodes, mut annotation_text, mut geometry_work) = (0usize, 0usize, 0u64);
        for binding in &request.sources {
            let dna = self
                .state
                .sequences
                .get(&binding.source.seq_id)
                .ok_or_else(|| invalid("rt_comparison_source_not_loaded"))?;
            if dna.is_protein_sequence()
                || dna
                    .molecule_type()
                    .is_some_and(|kind| kind.to_ascii_lowercase().contains("rna"))
            {
                return Err(invalid("rt_comparison_sources_require_annotated_dna"));
            }
            admit_annotation(
                dna,
                &mut annotation_nodes,
                &mut annotation_text,
                &mut geometry_work,
            )?;
            source_bases = source_bases
                .checked_add(dna.len())
                .filter(|v| *v <= 16_000_000)
                .ok_or_else(|| invalid("rt_comparison_source_base_budget"))?;
            source_features = source_features
                .checked_add(dna.features().len())
                .filter(|v| *v <= 100_000)
                .ok_or_else(|| invalid("rt_comparison_source_feature_budget"))?;
        }
        self.rt_comparison_templates(request, progress)
    }

    /// Compare supplied candidates in one explicit isoform universe; no mutation, ranking or order approval.
    pub fn compare_rt_primer_candidates(
        &self,
        request: CompareRtPrimerCandidatesRequest,
        progress: &mut dyn FnMut(OperationProgress) -> bool,
    ) -> Result<RtPrimerComparisonReport, EngineError> {
        let (templates, members, sources) = self.prepare_rt_comparison(&request, progress)?;
        // This establishes strict oligo/partition validation and exactly one common reference universe.
        let assessment = self.assess_rt_primer_candidates(request.assessment.clone(), progress)?;
        let mut report = RtPrimerComparisonReport { schema:RT_PRIMER_COMPARISON_SCHEMA.into(), content_sha256:String::new(),
            request_sha256:digest(&request)?, request:request.clone(), assessment:Box::new(assessment), sources, members,
            isoforms:Vec::new(), candidates:Vec::new(), transcript_base_comparisons:0, tm_model:Self::primer_tm_model_description(),
            secondary_structure_status:"not_assessed".into(), warnings:vec![
                "Each candidate is compared in the same explicit transcript universe and background request. Transcript membership is not RNA expression or molecular abundance.".into(),
                "All exact complementary full-tail sites are retained. Unknown CDS and internal annotation gaps are not UTR. Incomplete annotation is not an experimental coding-boundary validation.".into(),
                "Present sites can coexist with additional compatible ambiguous windows. Absent means no exact or ambiguity-compatible site in the selected mature transcript, not global specificity.".into(),
                "Region member counts overlap when an isoform has sites in several regions; repeated sites never multiply covered members. No single composite score or automatic ranking is assigned.".into(),
                "Tm/GC and complementary-run diagnostics use the existing engine heuristics; full-oligo Tm assumes a hypothetical full complement, not RT initiation at the target tail. Hairpin energies and calibrated secondary structure are not assessed.".into(),
                "Abundance, retrospective junction evidence, new reference searches beyond the supplied backgrounds, order readiness, GUI lanes and automatic redesign remain unassessed.".into(),
            ] };
        let mut binding_count = 0usize;
        let total = templates.len() * request.assessment.primers.len();
        let mut completed = 0;
        for primer in &request.assessment.primers {
            let core = primer
                .components
                .iter()
                .find(|c| c.role == "target_core")
                .expect("assessment validated");
            let tail = &primer.sequence.as_bytes()[core.start..core.end];
            let word = Self::reverse_complement_bytes(tail);
            let mut summary = RtPrimerCandidateComparison {
                primer_id: primer.primer_id.clone(),
                target_tail: metrics(tail),
                full_oligo: metrics(primer.sequence.as_bytes()),
                total_members: templates.len(),
                members_with_exact_site: 0,
                members_with_possible_only_site: 0,
                members_without_site: 0,
                region_member_counts: [
                    RtPrimerBindingRegion::Cds,
                    RtPrimerBindingRegion::FivePrimeUtr,
                    RtPrimerBindingRegion::ThreePrimeUtr,
                    RtPrimerBindingRegion::Mixed,
                    RtPrimerBindingRegion::Unknown,
                ]
                .into_iter()
                .map(|r| (r, 0))
                .collect(),
            };
            for (template, member) in templates.iter().zip(&report.members) {
                poll(progress, completed, total)?;
                let mut row = RtPrimerIsoformAssessment {
                    primer_id: primer.primer_id.clone(),
                    member_id: member.member_id.clone(),
                    presence: RtPrimerSitePresence::Absent,
                    possible_ambiguous_sites: 0,
                    bindings: Vec::new(),
                };
                let sequence = template.sequence.as_bytes();
                let windows = if sequence.len() < word.len() {
                    0
                } else {
                    sequence.len() - word.len() + 1
                };
                let cost = (windows as u64)
                    .checked_mul(word.len() as u64)
                    .ok_or_else(|| invalid("rt_comparison_work_budget"))?;
                report.transcript_base_comparisons = report
                    .transcript_base_comparisons
                    .checked_add(cost)
                    .filter(|v| *v <= request.budgets.max_transcript_base_comparisons)
                    .ok_or_else(|| invalid("rt_comparison_work_budget"))?;
                for start in 0..windows {
                    if start % 1024 == 0 {
                        poll(progress, completed, total)?;
                    }
                    let window = &sequence[start..start + word.len()];
                    if window == word {
                        binding_count += 1;
                        if binding_count > request.budgets.max_bindings {
                            return Err(invalid("rt_comparison_binding_budget"));
                        }
                        let end = start + word.len();
                        let mut binding = classify(start, end, &member.cds_ranges_0based);
                        let first = template
                            .local_exon_segments
                            .partition_point(|s| s.local_end_0based_exclusive <= start);
                        let last = template
                            .local_exon_segments
                            .partition_point(|s| s.local_start_0based < end);
                        let mapping = Self::map_transcript_local_interval(
                            &template.local_exon_segments[first..last],
                            template.strand == "-",
                            start,
                            end,
                        );
                        binding.source_ranges_0based = mapping
                            .source_ranges_0based
                            .into_iter()
                            .map(|r| (r.start_0based, r.end_0based_exclusive))
                            .collect();
                        row.bindings.push(binding);
                    } else if window.iter().zip(&word).all(|(&a, &b)| {
                        !IupacCode::from_letter(a)
                            .subset(IupacCode::from_letter(b))
                            .is_empty()
                    }) {
                        row.possible_ambiguous_sites += 1;
                    }
                }
                if !row.bindings.is_empty() {
                    row.presence = RtPrimerSitePresence::Present;
                    summary.members_with_exact_site += 1;
                    for region in row
                        .bindings
                        .iter()
                        .map(|b| b.region)
                        .collect::<BTreeSet<_>>()
                    {
                        *summary
                            .region_member_counts
                            .get_mut(&region)
                            .expect("all regions") += 1;
                    }
                } else if row.possible_ambiguous_sites > 0 {
                    row.presence = RtPrimerSitePresence::Unresolved;
                    summary.members_with_possible_only_site += 1;
                } else {
                    summary.members_without_site += 1;
                }
                report.isoforms.push(row);
                completed += 1;
            }
            report.candidates.push(summary);
        }
        poll(progress, total, total)?;
        report.content_sha256 = digest(&report)?;
        Ok(report)
    }
}

#[cfg(test)]
#[path = "rt_primer_comparison_tests.rs"]
mod tests;
