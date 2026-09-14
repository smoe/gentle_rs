//! Evidence-guided candidate ROIs upstream of exact reporter contrast planning.
//! Scores are inherited; only declared fixed-envelope signal comparisons are computed.

use super::*;
use gentle_protocol as gp;
use std::cmp::Reverse;
use std::io::Read;

fn invalid(message: impl Into<String>) -> EngineError {
    EngineError::invalid_input(message)
}

fn overlaps(a: &gp::GenomicRegionInterval, b: &gp::GenomicRegionInterval) -> bool {
    a.start_0based < b.end_0based_exclusive && b.start_0based < a.end_0based_exclusive
}

fn contains(a: &gp::GenomicRegionInterval, b: &gp::GenomicRegionInterval) -> bool {
    a.start_0based <= b.start_0based && a.end_0based_exclusive >= b.end_0based_exclusive
}

fn span(mut basis: gp::GenomicRegionInterval, start: u64, end: u64) -> gp::GenomicRegionInterval {
    basis.start_0based = start;
    basis.end_0based_exclusive = end;
    basis
}

fn hull(a: &gp::GenomicRegionInterval, b: &gp::GenomicRegionInterval) -> gp::GenomicRegionInterval {
    span(
        a.clone(),
        a.start_0based.min(b.start_0based),
        a.end_0based_exclusive.max(b.end_0based_exclusive),
    )
}

fn envelope(
    tss: &gp::GenomicRegionInterval,
    upstream: usize,
    downstream: usize,
) -> Result<gp::GenomicRegionInterval, EngineError> {
    let (left, right) = if tss.strand == gp::GenomicRegionStrand::Minus {
        (downstream, upstream)
    } else {
        (upstream, downstream)
    };
    let start = tss
        .start_0based
        .checked_sub(left as u64)
        .ok_or_else(|| invalid("Search envelope crosses the beginning of the contig"))?;
    let end = tss
        .end_0based_exclusive
        .checked_add(right as u64)
        .ok_or_else(|| invalid("Search envelope overflows genomic coordinates"))?;
    Ok(span(tss.clone(), start, end))
}

fn digest<T: Serialize>(value: &T) -> Result<String, EngineError> {
    serde_json::to_vec(value)
        .map(|v| crate::digest_utils::sha256_prefixed_bytes(&v))
        .map_err(|e| invalid(e.to_string()))
}

fn text(value: &str, label: &str) -> Result<(), EngineError> {
    if value.trim().is_empty() || value.len() > 4096 {
        return Err(invalid(format!("Missing or oversized {label}")));
    }
    Ok(())
}

fn normalize(
    mut request: FragmentSelectionRequest,
) -> Result<FragmentSelectionRequest, EngineError> {
    if request.schema != FRAGMENT_SELECTION_REQUEST_SCHEMA {
        return Err(invalid("Unsupported reporter-fragment selection schema"));
    }
    text(&request.locus.annotation_release, "annotation release")?;
    let p = &request.policy;
    if request.anchors.is_empty()
        || request.anchors.len() > 64
        || request.adjustments.len() > 128
        || request.signal_comparisons.len() > 64
        || p.upstream_bp > 100_000
        || p.downstream_bp > 100_000
        || p.maximum_extension_bp > 10_000
        || p.maximum_boundary_shift_bp > 10_000
        || p.maximum_candidates == 0
        || p.maximum_candidates > 512
        || p.maximum_evidence_intervals == 0
        || p.maximum_evidence_intervals > 100_000
        || p.preferred_length_bp == 0
        || p.maximum_length_bp < p.preferred_length_bp
        || p.maximum_length_bp > 100_000
        || p.flank_bp > 10_000
        || p.fallback_promoter_upstream_bp > p.upstream_bp
        || p.fallback_promoter_downstream_bp > p.downstream_bp
    {
        return Err(invalid("Invalid or excessive fragment-selection bounds"));
    }
    request
        .anchors
        .sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
    for anchor in &mut request.anchors {
        text(&anchor.transcript_id, "transcript id")?;
        anchor.seed_evidence_ids.sort();
        anchor.seed_evidence_ids.dedup();
        anchor.required_evidence_ids.sort();
        anchor.required_evidence_ids.dedup();
        if anchor.seed_region.is_some() && !anchor.seed_evidence_ids.is_empty() {
            return Err(invalid("Choose seed_region or seed_evidence_ids, not both"));
        }
    }
    if request
        .anchors
        .windows(2)
        .any(|w| w[0].transcript_id == w[1].transcript_id)
    {
        return Err(invalid("Duplicate anchor transcript"));
    }
    request
        .adjustments
        .sort_by(|a, b| a.adjustment_id.cmp(&b.adjustment_id));
    for adjustment in &mut request.adjustments {
        text(&adjustment.adjustment_id, "adjustment id")?;
        text(&adjustment.explanation, "boundary explanation")?;
        if adjustment.upstream_delta_bp.unsigned_abs() > p.maximum_boundary_shift_bp as u64
            || adjustment.downstream_delta_bp.unsigned_abs() > p.maximum_boundary_shift_bp as u64
        {
            return Err(invalid(
                "Boundary adjustment exceeds the declared shift limit",
            ));
        }
        adjustment.evidence_ids.sort();
        adjustment.evidence_ids.dedup();
    }
    if request
        .adjustments
        .windows(2)
        .any(|w| w[0].adjustment_id == w[1].adjustment_id)
    {
        return Err(invalid("Duplicate boundary adjustment id"));
    }
    request
        .signal_comparisons
        .sort_by(|a, b| a.comparison_id.cmp(&b.comparison_id));
    for c in &request.signal_comparisons {
        for (label, value) in [
            ("comparison id", &c.comparison_id),
            ("sample replicate", &c.sample_replicate_id),
            ("control replicate", &c.control_replicate_id),
            ("cell line", &c.cell_line),
            ("units", &c.units),
        ] {
            text(value, label)?;
        }
        if c.sample_lane_id == c.control_lane_id || !c.minimum_mean_difference.is_finite() {
            return Err(invalid("Invalid sample/control comparison"));
        }
    }
    if request
        .signal_comparisons
        .windows(2)
        .any(|w| w[0].comparison_id == w[1].comparison_id)
    {
        return Err(invalid("Duplicate signal comparison id"));
    }
    Ok(request)
}

impl GentleEngine {
    /// Propose immutable candidate geometry. Does not save ROIs, change existing plans,
    /// or authorize construct materialization.
    pub fn plan_reporter_fragment_selection(
        &self,
        request: FragmentSelectionRequest,
    ) -> Result<FragmentSelectionReport, EngineError> {
        let request = normalize(request)?;
        let document = self.fragment_selection_locus(&request.locus)?;
        let locus = document.locus();
        let p = &request.policy;
        let dna = self
            .state
            .sequences
            .get(&locus.seq_id)
            .ok_or_else(|| invalid("Source sequence is not loaded"))?;
        if dna.len() > 2_000_000 {
            return Err(invalid(
                "Load a locus of at most 2 Mb before fragment selection",
            ));
        }
        let evidence =
            self.fragment_selection_evidence(locus, &request.locus, p.maximum_evidence_intervals)?;
        let transcripts = self.fragment_selection_transcripts(locus)?;
        let vector = request
            .vector
            .as_ref()
            .map(|v| self.fragment_selection_vector(v, &locus.seq_id))
            .transpose()?;
        let mut anchors = Vec::new();
        let mut fixed_comparisons = Vec::new();
        let mut candidates = Vec::new();
        let mut findings = Vec::new();
        let mut used_adjustments = BTreeSet::new();
        // One isolated scratch engine for all candidate scans; never one project clone per row.
        let mut scratch = self.fork_detached_execution();
        for selected in &request.anchors {
            let tss = transcripts.get(&selected.transcript_id).ok_or_else(|| {
                invalid(format!(
                    "Unknown anchor transcript '{}'",
                    selected.transcript_id
                ))
            })?;
            let search = envelope(tss, p.upstream_bp, p.downstream_bp)?;
            self.local_projection_for_interval(&locus.seq_id, &search)?;
            let mut anchor = FragmentSelectionAnchorResult {
                transcript_id: selected.transcript_id.clone(),
                tss: tss.clone(),
                search_envelope: search.clone(),
                opposite_strand_transcript_ids: transcripts
                    .iter()
                    .filter(|(_, point)| point.strand != tss.strand && overlaps(&search, point))
                    .map(|(id, _)| id.clone())
                    .collect(),
                gene_structure: None,
                findings: Vec::new(),
            };
            if !anchor.opposite_strand_transcript_ids.is_empty() {
                anchor.findings.push("Opposite-strand TSSs overlap this envelope; their promoter identities remain separate.".into());
            }
            if p.purpose == FragmentSelectionPurpose::EndogenousPromoter {
                match self.fragment_gene_structure(locus, tss, &selected.transcript_id) {
                    Ok(audit) => anchor.gene_structure = Some(audit),
                    Err(e) => anchor
                        .findings
                        .push(format!("gene_structure_not_evaluated: {}", e.message)),
                }
            }
            for comparison in &request.signal_comparisons {
                fixed_comparisons.push(self.fragment_fixed_comparison(locus, &anchor, comparison)?);
            }
            let enriched = fixed_comparisons.iter().any(|c| {
                c.transcript_id == selected.transcript_id && c.passes_descriptive_rule == Some(true)
            });
            let mut required = Vec::<gp::GenomicRegionInterval>::new();
            if p.purpose == FragmentSelectionPurpose::EndogenousPromoter {
                let promoter_features = evidence
                    .iter()
                    .filter(|e| {
                        e.available
                            && e.kind == FragmentEvidenceKind::Annotation
                            && e.label.eq_ignore_ascii_case("promoter")
                            && contains(&e.interval, tss)
                    })
                    .collect::<Vec<_>>();
                if promoter_features.is_empty() {
                    required.push(envelope(
                        tss,
                        p.fallback_promoter_upstream_bp,
                        p.fallback_promoter_downstream_bp,
                    )?);
                    anchor.findings.push("promoter_context_provisional: no covering promoter annotation; configured TSS context retained, not a proven core promoter.".into());
                } else {
                    required.extend(promoter_features.iter().map(|e| e.interval.clone()));
                }
            }
            for id in &selected.required_evidence_ids {
                let row = evidence
                    .iter()
                    .find(|e| &e.evidence_id == id)
                    .ok_or_else(|| invalid(format!("Unknown required evidence '{id}'")))?;
                if !row.available {
                    return Err(invalid(format!("Required evidence '{id}' is unavailable")));
                }
                required.push(row.interval.clone());
            }
            let mut seeds = evidence
                .iter()
                .filter(|e| {
                    e.may_seed_boundary
                        && overlaps(&e.interval, &search)
                        && (selected.seed_evidence_ids.is_empty()
                            || selected.seed_evidence_ids.contains(&e.evidence_id))
                })
                .map(|e| (e.evidence_id.clone(), e.interval.clone(), true))
                .collect::<Vec<_>>();
            for id in &selected.seed_evidence_ids {
                if !seeds.iter().any(|(seed, _, _)| seed == id) {
                    return Err(invalid(format!(
                        "Seed '{id}' is absent, unavailable, raw-coverage-only, or outside the search envelope"
                    )));
                }
            }
            if let Some(region) = &selected.seed_region {
                let mut verified = region.clone();
                super::genomic_regions::recompute_region_digests(&mut verified)?;
                if verified.identity_sha256 != region.identity_sha256
                    || verified.content_sha256 != region.content_sha256
                {
                    return Err(invalid("Seed ROI content digest mismatch"));
                }
                let projection =
                    self.local_projection_for_interval(&locus.seq_id, &region.interval)?;
                if region
                    .local_projection
                    .as_ref()
                    .is_some_and(|old| old != &projection)
                {
                    return Err(invalid("Seed ROI projection is stale"));
                }
                if !overlaps(&region.interval, &search) {
                    return Err(invalid("Seed ROI is outside its declared TSS envelope"));
                }
                seeds = vec![(region.region_id.clone(), region.interval.clone(), false)];
            }
            if seeds.is_empty() {
                anchor.findings.push("No annotation/model seed within the envelope; raw coverage alone does not define a regulatory boundary.".into());
            }
            for (seed_id, seed, protect_seed) in seeds {
                let mut requirements = required.clone();
                if protect_seed {
                    requirements.push(seed.clone());
                }
                let core = requirements.iter().fold(seed.clone(), |a, b| hull(&a, b));
                let padded = span(
                    core.clone(),
                    core.start_0based.saturating_sub(p.flank_bp as u64),
                    core.end_0based_exclusive.saturating_add(p.flank_bp as u64),
                );
                let mut base_interval = if selected.seed_region.is_some() {
                    seed.clone()
                } else {
                    padded
                };
                base_interval.strand = tss.strand;
                let base_id = format!(
                    "fragment_{}",
                    &digest(&(
                        &selected.transcript_id,
                        &seed_id,
                        &base_interval,
                        "evidence_hull"
                    ))?[7..23]
                );
                let mut variants = vec![("evidence_hull".to_string(),base_interval.clone(),None,vec!["Evidence geometry plus declared flanks; not selected by maximizing signal.".into()])];
                if base_interval.end_0based_exclusive - base_interval.start_0based
                    > p.preferred_length_bp as u64
                {
                    if core.end_0based_exclusive - core.start_0based <= p.preferred_length_bp as u64
                    {
                        variants.push(("compact_without_optional_flanks".into(),core.clone(),Some(base_id.clone()),vec!["Optional padding removed; the selected evidence/context hull is unchanged.".into()]));
                    } else {
                        anchor.findings.push(format!("compact_unavailable:{seed_id}: required evidence/context spans {} bp; longer context retained explicitly",core.end_0based_exclusive-core.start_0based));
                    }
                }
                for adjustment in request
                    .adjustments
                    .iter()
                    .filter(|a| a.transcript_id == selected.transcript_id && a.seed_id == seed_id)
                {
                    used_adjustments.insert(adjustment.adjustment_id.clone());
                    for id in &adjustment.evidence_ids {
                        if !evidence.iter().any(|e| &e.evidence_id == id) {
                            return Err(invalid(format!(
                                "Adjustment cites unknown evidence '{id}'"
                            )));
                        }
                    }
                    let (left, right) = if tss.strand == gp::GenomicRegionStrand::Minus {
                        (adjustment.downstream_delta_bp, adjustment.upstream_delta_bp)
                    } else {
                        (adjustment.upstream_delta_bp, adjustment.downstream_delta_bp)
                    };
                    let from = base_interval
                        .start_0based
                        .checked_add_signed(-left)
                        .ok_or_else(|| invalid("Boundary adjustment crosses coordinate zero"))?;
                    let to = base_interval
                        .end_0based_exclusive
                        .checked_add_signed(right)
                        .ok_or_else(|| invalid("Boundary adjustment overflows coordinates"))?;
                    if from >= to {
                        return Err(invalid("Boundary adjustment would erase/invert the insert"));
                    }
                    variants.push((
                        format!("human:{}", adjustment.adjustment_id),
                        span(base_interval.clone(), from, to),
                        Some(base_id.clone()),
                        vec![format!(
                            "{:?}: {} (evidence: {})",
                            adjustment.reason,
                            adjustment.explanation,
                            adjustment.evidence_ids.join(", ")
                        )],
                    ));
                }
                if request
                    .vector
                    .as_ref()
                    .is_some_and(|v| v.suggest_restriction_adjustments)
                    && let Some(v) = &vector
                {
                    for site in &v.source_sites {
                        if !overlaps(&site.interval, &base_interval) {
                            continue;
                        }
                        for (name, from, to) in [
                            (
                                "trim_genomic_left",
                                site.interval.end_0based_exclusive,
                                base_interval.end_0based_exclusive,
                            ),
                            (
                                "trim_genomic_right",
                                base_interval.start_0based,
                                site.interval.start_0based,
                            ),
                        ] {
                            if from >= to
                                || from.abs_diff(base_interval.start_0based)
                                    > p.maximum_boundary_shift_bp as u64
                                || to.abs_diff(base_interval.end_0based_exclusive)
                                    > p.maximum_boundary_shift_bp as u64
                            {
                                continue;
                            }
                            let proposed = span(base_interval.clone(), from, to);
                            if requirements.iter().any(|r| !contains(&proposed, r)) {
                                continue;
                            }
                            variants.push((format!("restriction:{}:{name}:{}",site.enzyme,site.interval.start_0based),proposed,Some(base_id.clone()),vec![format!("Exclude internal {} recognition footprint without removing required evidence/context. PCR-added MCS sites remain separate from genomic insert boundaries.",site.enzyme)]));
                        }
                    }
                }
                for (variant, mut interval, parent, reasons) in variants {
                    interval.strand = tss.strand;
                    if candidates.len() >= p.maximum_candidates {
                        return Err(invalid(
                            "Candidate budget exceeded; select explicit seed_evidence_ids or increase the declared budget",
                        ));
                    }
                    let candidate = scratch.engine_mut().fragment_selection_candidate(
                        &request,
                        locus,
                        &anchor,
                        &evidence,
                        &requirements,
                        &seed_id,
                        &variant,
                        interval,
                        parent,
                        reasons,
                        enriched,
                        &transcripts,
                        vector.as_ref(),
                    )?;
                    candidates.push(candidate);
                }
            }
            anchors.push(anchor);
        }
        if request
            .adjustments
            .iter()
            .any(|a| !used_adjustments.contains(&a.adjustment_id))
        {
            return Err(invalid(
                "A boundary adjustment did not resolve to an anchor/seed; no adjustment was silently ignored",
            ));
        }
        candidates.sort_by_key(|c| {
            (
                Reverse(c.ranking.preserves_required_context),
                Reverse(c.ranking.descriptive_enrichment_supported),
                Reverse(c.ranking.retains_annotation),
                Reverse(c.ranking.retains_model_site),
                c.ranking.bisected_feature_count,
                Reverse(c.ranking.fits_preferred_length),
                c.ranking.length_bp,
                c.candidate_id.clone(),
            )
        });
        findings.push("Ranking compares explicit preserved-context, descriptive-support (any declared comparison passing), annotation, model, biological-feature bisection, compactness, length and stable-ID keys; raw signal-bin boundaries do not affect ranking. No independent-source count or composite biological score.".into());
        if request.signal_comparisons.is_empty() {
            findings.push("CUT&RUN enrichment not_evaluated: no explicit matched-control comparison supplied. Raw coverage remains inspectable.".into());
        }
        let request_sha256 = digest(&request)?;
        let mut region_set = gp::GenomicRegionSet {
            schema: gp::GENOMIC_REGION_SET_SCHEMA.into(),
            set_id: format!("fragment_selection_{}", &request_sha256[7..23]),
            regions: candidates
                .iter()
                .filter(|c| c.blockers.is_empty())
                .map(|c| c.region.clone())
                .collect(),
            ..Default::default()
        };
        region_set
            .regions
            .sort_by(|a, b| a.region_id.cmp(&b.region_id));
        region_set
            .regions
            .dedup_by(|a, b| a.region_id == b.region_id);
        super::genomic_regions::recompute_set_digest(&mut region_set)?;
        let mut report = FragmentSelectionReport { schema: FRAGMENT_SELECTION_REPORT_SCHEMA.into(), request, request_sha256, proposal_sha256: String::new(), source_sequence_sha256: crate::digest_utils::sha256_prefixed_bytes(dna.forward_bytes()), evidence, anchors, fixed_comparisons, candidates, proposed_region_set: region_set, vector, findings, non_claims: vec![
            "Candidate geometry does not prove autonomous promoter function, direct binding, cofactor cooperation or a change in luciferase activity.".into(),
            "Fixed-envelope enrichment may prioritize a neighbourhood; it does not localize the supporting bases or independently validate a selected compact insert.".into(),
            "Raw occupancy, model scores and provider annotations remain separate. A derived motif/occupancy join is not an additional independent experiment.".into(),
            "Only read-only proposals are produced. Saving a region and approving downstream exact constructs are separate actions; existing selections and approvals are unchanged.".into(),
        ] };
        report.proposal_sha256 = digest(&report)?;
        Ok(report)
    }

    #[allow(clippy::too_many_arguments)]
    fn fragment_selection_candidate(
        &mut self,
        request: &FragmentSelectionRequest,
        locus: &gp::GeneLocusEvidenceDisplayReport,
        anchor: &FragmentSelectionAnchorResult,
        evidence: &[FragmentSelectionEvidence],
        required: &[gp::GenomicRegionInterval],
        seed_id: &str,
        variant: &str,
        interval: gp::GenomicRegionInterval,
        parent: Option<String>,
        mut reasons: Vec<String>,
        enriched: bool,
        transcripts: &BTreeMap<String, gp::GenomicRegionInterval>,
        vector: Option<&FragmentSelectionVectorResult>,
    ) -> Result<ReporterFragmentCandidate, EngineError> {
        let p = &request.policy;
        let mut blockers = Vec::new();
        let outer = span(
            anchor.search_envelope.clone(),
            anchor
                .search_envelope
                .start_0based
                .saturating_sub(p.maximum_extension_bp as u64),
            anchor
                .search_envelope
                .end_0based_exclusive
                .saturating_add(p.maximum_extension_bp as u64),
        );
        if !contains(&outer, &interval) {
            blockers.push("outside_declared_extension_limit".into());
        }
        let length = interval.end_0based_exclusive - interval.start_0based;
        if length > p.maximum_length_bp as u64 {
            blockers.push("maximum_insert_length_exceeded".into());
        }
        if required.iter().any(|r| !contains(&interval, r)) {
            blockers.push("required_evidence_or_promoter_context_removed".into());
        }
        let projection = self.local_projection_for_interval(&locus.seq_id, &interval)?;
        let dna = self
            .state
            .sequences
            .get(&locus.seq_id)
            .ok_or_else(|| invalid("Source sequence missing"))?;
        let raw = &dna.get_forward_string()[projection.local_start_0based as usize
            ..projection.local_end_0based_exclusive as usize];
        let sequence = if projection.local_strand == gp::GenomicRegionStrand::Minus {
            Self::reverse_complement(raw)
        } else {
            raw.to_string()
        };
        let mut retained = Vec::new();
        let mut excluded = Vec::new();
        let mut bisected = Vec::new();
        let mut references = Vec::new();
        for e in evidence.iter().filter(|e| overlaps(&e.interval, &outer)) {
            if contains(&interval, &e.interval) {
                retained.push(e.evidence_id.clone());
            } else {
                excluded.push(e.evidence_id.clone());
                if overlaps(&interval, &e.interval) {
                    bisected.push(e.evidence_id.clone());
                }
            }
            if overlaps(&interval, &e.interval) {
                references.push(gp::GenomicRegionEvidenceReference {
                    evidence_id: e.evidence_id.clone(),
                    source_kind: format!("{:?}", e.kind),
                    source_id: e.source_id.clone(),
                    source_sha256: Some(e.source_sha256.clone()),
                    report_id: Some(request.locus.panel_id.clone()),
                    source_release: Some(request.locus.annotation_release.clone()),
                    availability: if e.available {
                        gp::GenomicRegionEvidenceAvailability::Available
                    } else {
                        gp::GenomicRegionEvidenceAvailability::Unavailable
                    },
                    evidence_statement: e.statement.clone(),
                    ..Default::default()
                });
            }
        }
        if !contains(&anchor.search_envelope, &interval) {
            reasons.push("Evidence/context crosses the initial search envelope; extension is explicit and within the separately checked limit.".into());
        }
        if !bisected.is_empty() {
            reasons.push("Some contextual features are bisected; their complete source bounds remain in the evidence inventory.".into());
        }
        if length > p.preferred_length_bp as u64 {
            reasons.push("Longer-context alternative exceeds the preferred compact length; no selected context was silently removed to meet that preference.".into());
        }
        if request.policy.purpose == FragmentSelectionPurpose::EndogenousPromoter {
            if let Some(audit) = &anchor.gene_structure {
                for warning in &audit.warnings {
                    reasons.push(format!(
                        "Transcript context, not an assertion for every insert: {}",
                        warning.detail
                    ));
                }
                if audit
                    .cds_start_codon_source_ranges_0based
                    .iter()
                    .any(|(start, end)| {
                        (*start as u64) < projection.local_end_0based_exclusive
                            && (*end as u64) > projection.local_start_0based
                    })
                {
                    reasons.push("Candidate includes all or part of the annotated endogenous start codon; inspect the reporter translation context before cloning.".into());
                }
            } else {
                reasons.push(
                    "Gene-structure audit unavailable; inspect CDS/UTR context before cloning."
                        .into(),
                );
            }
        }
        let candidate_id = format!(
            "fragment_{}",
            &digest(&(&anchor.transcript_id, seed_id, &interval, variant))?[7..23]
        );
        references.push(gp::GenomicRegionEvidenceReference { evidence_id: "selection_policy".into(), source_kind: "reporter_fragment_selection_request".into(), source_id: digest(request)?, source_sha256: Some(digest(request)?), evidence_statement: format!("Derived boundary proposal {variant}; exact policy, manual adjustments and any parent ROI are bound in the selection request. Not an unmodified set-theoretic hull."), ..Default::default() });
        let mut region = gp::GenomicRegionOfInterest {
            schema: gp::GENOMIC_REGION_OF_INTEREST_SCHEMA.into(),
            region_id: candidate_id.clone(),
            label: Some(format!("{} / {}", anchor.transcript_id, variant)),
            interval: interval.clone(),
            local_projection: Some(projection),
            purpose: gp::GenomicRegionPurpose::ReporterCandidate,
            selection_method: gp::GenomicRegionSelectionMethod::Derived,
            evidence: references,
            notes: reasons.clone(),
            ..Default::default()
        };
        super::genomic_regions::recompute_region_digests(&mut region)?;
        let cloning = match (&request.vector, vector) {
            (Some(v), Some(binding)) => {
                Some(self.fragment_candidate_cloning(&sequence, v, binding)?)
            }
            _ => None,
        };
        let retains = |kind| {
            evidence
                .iter()
                .any(|e| e.available && e.kind == kind && retained.contains(&e.evidence_id))
        };
        let ranking = FragmentCandidateRanking {
            preserves_required_context: required.iter().all(|r| contains(&interval, r)),
            descriptive_enrichment_supported: enriched,
            retains_annotation: retains(FragmentEvidenceKind::Annotation),
            retains_model_site: retains(FragmentEvidenceKind::ModelSite),
            bisected_feature_count: evidence
                .iter()
                .filter(|e| {
                    e.available
                        && e.kind != FragmentEvidenceKind::RawCoverage
                        && bisected.contains(&e.evidence_id)
                })
                .count(),
            fits_preferred_length: length <= p.preferred_length_bp as u64,
            length_bp: length as usize,
        };
        Ok(ReporterFragmentCandidate {
            candidate_id,
            transcript_id: anchor.transcript_id.clone(),
            member_transcript_ids: transcripts
                .iter()
                .filter(|(_, t)| t.strand == interval.strand && contains(&interval, t))
                .map(|(id, _)| id.clone())
                .collect(),
            seed_id: seed_id.into(),
            variant: variant.into(),
            parent_candidate_id: parent,
            region,
            sequence_sha256: crate::digest_utils::sha256_prefixed_bytes(sequence.as_bytes()),
            sequence_5prime_to_3prime: sequence,
            retained_evidence_ids: retained,
            excluded_evidence_ids: excluded,
            bisected_evidence_ids: bisected,
            ranking,
            blockers,
            reasons,
            cloning,
        })
    }

    fn fragment_fixed_comparison(
        &self,
        locus: &gp::GeneLocusEvidenceDisplayReport,
        anchor: &FragmentSelectionAnchorResult,
        comparison: &FragmentSignalComparison,
    ) -> Result<FragmentFixedComparison, EngineError> {
        let lanes = locus
            .occupancy_groups
            .iter()
            .flat_map(|g| &g.lanes)
            .collect::<Vec<_>>();
        let lane = |id: &str| -> Result<&gp::GeneLocusOccupancyLane, EngineError> {
            let matches = lanes
                .iter()
                .filter(|l| l.lane.lane_id == id)
                .copied()
                .collect::<Vec<_>>();
            if matches.len() != 1 {
                return Err(invalid(format!(
                    "Signal lane '{id}' must resolve exactly once"
                )));
            }
            Ok(matches[0])
        };
        let sample = lane(&comparison.sample_lane_id)?;
        let control = lane(&comparison.control_lane_id)?;
        if sample.source_sha256.is_some() && sample.source_sha256 == control.source_sha256 {
            return Err(invalid("Sample and control are the same source content"));
        }
        let projection =
            self.local_projection_for_interval(&locus.seq_id, &anchor.search_envelope)?;
        let start = projection.local_start_0based as usize;
        let end = projection.local_end_0based_exclusive as usize;
        let mean = |lane: &gp::GeneLocusOccupancyLane| -> Result<Option<f64>, EngineError> {
            if lane.state != gp::GeneLocusOccupancyLaneState::Available
                || lane.cell_line_label.as_deref() != Some(&comparison.cell_line)
                || lane.lane.interval_count != lane.lane.intervals.len()
            {
                return Ok(None);
            }
            let mut values = vec![None; end - start];
            for row in &lane.lane.intervals {
                let from = row.local_start_1based.saturating_sub(1).max(start);
                let to = row.local_end_1based.min(end);
                if from >= to {
                    continue;
                }
                let Some(score) = row.score else {
                    return Ok(None);
                };
                for value in &mut values[from - start..to - start] {
                    if value.is_some() {
                        return Err(invalid(
                            "Overlapping signal bins cannot be treated as per-base means",
                        ));
                    }
                    *value = Some(score);
                }
            }
            if comparison.missing_policy == FragmentSignalMissingPolicy::RequireComplete
                && values.iter().any(Option::is_none)
            {
                return Ok(None);
            }
            let mean =
                values.into_iter().map(|v| v.unwrap_or(0.0)).sum::<f64>() / (end - start) as f64;
            if !mean.is_finite() {
                return Err(invalid("Non-finite fixed-envelope signal mean"));
            }
            Ok(Some(mean))
        };
        let sample_mean = mean(sample)?;
        let control_mean = mean(control)?;
        let difference = sample_mean.zip(control_mean).map(|(a, b)| a - b);
        if difference.is_some_and(|d| !d.is_finite()) {
            return Err(invalid("Signal difference overflow"));
        }
        Ok(FragmentFixedComparison { comparison: comparison.clone(), transcript_id: anchor.transcript_id.clone(), measurement_window: anchor.search_envelope.clone(), sample_source_sha256: sample.source_sha256.clone(), control_source_sha256: control.source_sha256.clone(), sample_mean, control_mean, mean_difference: difference, passes_descriptive_rule: difference.map(|d| d > comparison.minimum_mean_difference), status: if difference.is_some() { "descriptive_enrichment; caller-declared replicate identities and compatible units; no significance or localization claim" } else { "not_evaluated: missing/incomplete signal or cell-line metadata" }.into() })
    }

    fn fragment_selection_vector(
        &self,
        request: &FragmentSelectionVector,
        seq_id: &str,
    ) -> Result<FragmentSelectionVectorResult, EngineError> {
        let resolved = self.resolve_reporter_backbone(
            &request.seq_id,
            None,
            true,
            Some(&request.catalog_id),
            request.helper_catalog_path.as_deref(),
        )?;
        let validation = resolved
            .validation
            .ok_or_else(|| invalid("Exact vector validation is required"))?;
        if validation.status != ReporterVectorValidationStatus::Verified {
            return Err(invalid(
                "Vector does not match the catalog's exact sequence/MCS expectations",
            ));
        }
        let vector = self
            .state
            .sequences
            .get(&request.seq_id)
            .ok_or_else(|| invalid("Vector is not loaded"))?;
        let (catalog, path) =
            Self::open_helper_genome_catalog(request.helper_catalog_path.as_deref())?;
        let catalog_bytes = std::fs::read(path).map_err(|e| invalid(e.to_string()))?;
        let expectation = catalog
            .helper_vector_sequence_expectation(&request.catalog_id)
            .map_err(invalid)?
            .ok_or_else(|| invalid("No exact vector expectation"))?;
        let mcs = expectation
            .required_features
            .iter()
            .find(|f| f.id.eq_ignore_ascii_case("multiple_cloning_region"))
            .ok_or_else(|| invalid("Vector catalog lacks the MCS"))?;
        let (from, to, _) = Self::find_reporter_vector_feature(vector, mcs)
            .ok_or_else(|| invalid("Verified MCS is absent"))?;
        let mcs_start = from.saturating_sub(1);
        let suggestions = self.restriction_cloning_vector_enzyme_suggestions(&request.seq_id)?;
        let enzymes = suggestions
            .recommended_directed_pairs
            .iter()
            .filter(|p| {
                (mcs_start..to).contains(&p.forward_cut_position_0based)
                    && (mcs_start..to).contains(&p.reverse_cut_position_0based)
            })
            .flat_map(|p| [p.forward_enzyme.clone(), p.reverse_enzyme.clone()])
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        if enzymes.len() > 64 {
            return Err(invalid("MCS enzyme budget exceeded"));
        }
        let mut source_sites = Vec::new();
        if !enzymes.is_empty() {
            let scan = self.find_restriction_sites(
                SequenceScanTarget::SeqId {
                    seq_id: seq_id.into(),
                    span_start_0based: None,
                    span_end_0based_exclusive: None,
                },
                &enzymes,
                None,
                false,
                None,
                None,
            )?;
            for row in scan.rows {
                let local_strand = if row.forward_strand {
                    gp::GenomicRegionStrand::Plus
                } else {
                    gp::GenomicRegionStrand::Minus
                };
                let (interval, _) = self.interval_and_projection_from_local(
                    seq_id,
                    row.source_recognition_start_0based as u64,
                    row.source_recognition_end_0based_exclusive as u64,
                    local_strand,
                    None,
                )?;
                source_sites.push(FragmentRestrictionSite {
                    enzyme: row.enzyme_name,
                    interval,
                    local_top_cut_0based: row.source_forward_cut_0based,
                    local_bottom_cut_0based: row.source_reverse_cut_0based,
                });
            }
        }
        source_sites.sort_by(|a, b| {
            (
                &a.enzyme,
                a.interval.start_0based,
                a.interval.end_0based_exclusive,
            )
                .cmp(&(
                    &b.enzyme,
                    b.interval.start_0based,
                    b.interval.end_0based_exclusive,
                ))
        });
        Ok(FragmentSelectionVectorResult {
            validation,
            sequence_sha256: crate::digest_utils::sha256_prefixed_bytes(vector.forward_bytes()),
            catalog_sha256: crate::digest_utils::sha256_prefixed_bytes(&catalog_bytes),
            mcs_start_0based: mcs_start,
            mcs_end_0based_exclusive: to,
            source_sites,
        })
    }

    fn fragment_candidate_cloning(
        &mut self,
        sequence: &str,
        request: &FragmentSelectionVector,
        vector: &FragmentSelectionVectorResult,
    ) -> Result<PromoterReporterPanelCloningStrategyReport, EngineError> {
        let id = format!(
            "fragment_selection_{}",
            crate::digest_utils::sha256_prefixed_bytes(sequence.as_bytes())
                .trim_start_matches("sha256:")
        );
        if self.state.sequences.contains_key(&id) {
            return Err(invalid(
                "Temporary candidate sequence ID collides with project state",
            ));
        }
        self.state.sequences.insert(
            id.clone(),
            DNAsequence::from_sequence(sequence).map_err(|e| invalid(e.to_string()))?,
        );
        let result =
            self.restriction_cloning_panel_strategy(&request.seq_id, std::slice::from_ref(&id));
        self.state.sequences.remove(&id);
        let mut report = result?;
        report.pair_evaluations.retain(|p| {
            (vector.mcs_start_0based..vector.mcs_end_0based_exclusive)
                .contains(&p.pair.forward_cut_position_0based)
                && (vector.mcs_start_0based..vector.mcs_end_0based_exclusive)
                    .contains(&p.pair.reverse_cut_position_0based)
        });
        report.selected_pair = report
            .pair_evaluations
            .iter()
            .find(|p| p.compatible)
            .map(|p| p.pair.clone());
        report.strategy = if report.selected_pair.is_some() {
            PromoterReporterPanelCloningStrategy::DirectionalRestriction
        } else {
            PromoterReporterPanelCloningStrategy::Gibson
        };
        report.fallback_reason = report.selected_pair.is_none().then(|| "No internal-site-free pair inside the validated MCS; consider another boundary or an explicitly chosen alternative assembly method.".into());
        report.warnings.push("This route assumes restriction sites added by PCR primers. Native genomic sites are context only; no native-fragment ligation, primer-tail design, end-digestion efficiency or reaction simulation is claimed.".into());
        Ok(report)
    }

    fn fragment_gene_structure(
        &self,
        locus: &gp::GeneLocusEvidenceDisplayReport,
        tss: &gp::GenomicRegionInterval,
        transcript_id: &str,
    ) -> Result<PromoterReporterPanelExtendedBoundaryAudit, EngineError> {
        let projection = self.local_projection_for_interval(&locus.seq_id, tss)?;
        let dna = self
            .state
            .sequences
            .get(&locus.seq_id)
            .ok_or_else(|| invalid("Missing source sequence"))?;
        let candidate = PromoterReporterFragmentCandidate {
            transcript_id: transcript_id.into(),
            strand: projection.local_strand.bed_value().into(),
            tss_local_0based: projection.local_start_0based as usize,
            start_0based: projection.local_start_0based as usize,
            end_0based_exclusive: projection.local_end_0based_exclusive as usize,
            variant_start_0based: projection.local_start_0based as usize,
            variant_end_0based_exclusive: projection.local_end_0based_exclusive as usize,
            ..Default::default()
        };
        Self::promoter_reporter_panel_extended_geometry(
            dna,
            &candidate,
            &PromoterReporterPanelExtendedBoundaryPolicy {
                kind: PromoterReporterPanelExtendedBoundaryKind::CanonicalCdsStartExclusive,
                transcript_id: transcript_id.into(),
            },
        )
        .map(|(_, _, audit)| audit)
    }

    fn fragment_selection_locus(
        &self,
        binding: &FragmentSelectionLocus,
    ) -> Result<crate::locus_report::LocusDocument, EngineError> {
        let mut bytes = Vec::new();
        std::fs::File::open(&binding.path)
            .map_err(|e| invalid(e.to_string()))?
            .take(128 * 1024 * 1024 + 1)
            .read_to_end(&mut bytes)
            .map_err(|e| invalid(e.to_string()))?;
        if bytes.len() > 128 * 1024 * 1024
            || crate::digest_utils::sha256_prefixed_bytes(&bytes) != binding.sha256
        {
            return Err(invalid("Locus file hash mismatch or size limit exceeded"));
        }
        let document = crate::locus_report::LocusDocument::from_json(&bytes).map_err(invalid)?;
        let locus = document.locus();
        if locus.panel_id != binding.panel_id
            || locus.isoform_evidence.annotation_release.as_deref()
                != Some(&binding.annotation_release)
        {
            return Err(invalid("Locus panel or annotation-release mismatch"));
        }
        let dna = self
            .state
            .sequences
            .get(&locus.seq_id)
            .ok_or_else(|| invalid("Locus source sequence is not loaded"))?;
        if dna.is_circular() {
            return Err(invalid(
                "Genomic reporter selection requires a linear source locus",
            ));
        }
        let anchor = self.sequence_genome_anchor_summary(&locus.seq_id)?;
        crate::locus_report::verify_live_binding(
            locus,
            &locus.seq_id,
            &crate::locus_report::sequence_binding(dna, Some(&anchor)),
        )
        .map_err(invalid)?;
        Ok(document)
    }

    fn fragment_selection_transcripts(
        &self,
        locus: &gp::GeneLocusEvidenceDisplayReport,
    ) -> Result<BTreeMap<String, gp::GenomicRegionInterval>, EngineError> {
        let mut result = BTreeMap::new();
        for row in &locus.isoform_evidence.transcripts {
            let first_id = row
                .exon_family_ids_5_to_3
                .first()
                .ok_or_else(|| invalid("Transcript lacks first-exon geometry"))?;
            let exon = locus
                .isoform_evidence
                .exon_families
                .iter()
                .find(|e| &e.exon_family_id == first_id)
                .ok_or_else(|| invalid("Missing transcript first-exon family"))?;
            let genomic_strand = match row.strand.as_str() {
                "+" => gp::GenomicRegionStrand::Plus,
                "-" => gp::GenomicRegionStrand::Minus,
                _ => return Err(invalid("Transcript strand is unavailable")),
            };
            if exon.local_start_1based == 0 {
                return Err(invalid("Invalid first-exon coordinates"));
            }
            let (mut interval, _) = self.interval_and_projection_from_local(
                &locus.seq_id,
                (exon.local_start_1based - 1) as u64,
                exon.local_end_1based as u64,
                gp::GenomicRegionStrand::Unstranded,
                None,
            )?;
            interval.strand = genomic_strand;
            if genomic_strand == gp::GenomicRegionStrand::Minus {
                interval.start_0based = interval.end_0based_exclusive - 1;
            } else {
                interval.end_0based_exclusive = interval.start_0based + 1;
            }
            let projection = self.local_projection_for_interval(&locus.seq_id, &interval)?;
            let dna = self
                .state
                .sequences
                .get(&locus.seq_id)
                .ok_or_else(|| invalid("Missing source DNA"))?;
            let live = dna
                .features()
                .iter()
                .enumerate()
                .filter(|(index, f)| {
                    Self::is_splicing_transcript_feature(f)
                        && Self::feature_transcript_id(f, *index) == row.transcript_id
                })
                .collect::<Vec<_>>();
            if live.len() != 1 {
                return Err(invalid(
                    "Transcript must resolve uniquely in current annotations; refresh the locus report",
                ));
            }
            let (_, feature) = live[0];
            let mut live_ranges = Vec::new();
            collect_location_ranges_usize(&feature.location, &mut live_ranges);
            live_ranges.sort_unstable();
            let live_tss = if feature_is_reverse(feature) {
                live_ranges.last().and_then(|(_, end)| end.checked_sub(1))
            } else {
                live_ranges.first().map(|(start, _)| *start)
            };
            if live_tss != Some(projection.local_start_0based as usize) {
                return Err(invalid(
                    "Reported first exon disagrees with the current transcript TSS",
                ));
            }
            let mut reported_ranges = row
                .exon_family_ids_5_to_3
                .iter()
                .map(|id| {
                    let e = locus
                        .isoform_evidence
                        .exon_families
                        .iter()
                        .find(|e| &e.exon_family_id == id)
                        .ok_or_else(|| invalid("Missing exon family"))?;
                    if e.local_start_1based == 0 {
                        return Err(invalid("Invalid exon geometry"));
                    }
                    Ok((e.local_start_1based - 1, e.local_end_1based))
                })
                .collect::<Result<Vec<_>, EngineError>>()?;
            reported_ranges.sort_unstable();
            if live_ranges != reported_ranges
                || feature_is_reverse(feature)
                    != (projection.local_strand == gp::GenomicRegionStrand::Minus)
            {
                return Err(invalid(
                    "Transcript geometry/strand changed; refresh the locus report",
                ));
            }
            if result.insert(row.transcript_id.clone(), interval).is_some() {
                return Err(invalid("Duplicate transcript identity in locus"));
            }
        }
        Ok(result)
    }

    fn fragment_selection_evidence(
        &self,
        locus: &gp::GeneLocusEvidenceDisplayReport,
        source: &FragmentSelectionLocus,
        limit: usize,
    ) -> Result<Vec<FragmentSelectionEvidence>, EngineError> {
        let mut rows = Vec::new();
        if let Some(regulation) = &locus.ensembl_regulation {
            let available =
                regulation.availability == gp::GeneLocusEnsemblRegulationAvailability::Available;
            if available
                && !regulation.source_binding.as_ref().is_some_and(|b| {
                    b.content_identity_verified
                        && b.source.assembly_name == locus.isoform_evidence.assembly
                })
            {
                return Err(invalid("Unverified Ensembl annotation source"));
            }
            for row in &regulation.rows {
                if row.assembly_name != locus.isoform_evidence.assembly
                    || row.core_genomic_start_1based == 0
                    || row.core_genomic_end_1based < row.core_genomic_start_1based
                {
                    return Err(invalid(
                        "Invalid regulatory annotation geometry or assembly",
                    ));
                }
                let (mut interval, _) = self.interval_and_projection_from_local(
                    &locus.seq_id,
                    (row.displayed_local_start_1based.saturating_sub(1)) as u64,
                    row.displayed_local_end_1based as u64,
                    gp::GenomicRegionStrand::Unstranded,
                    None,
                )?;
                if interval.start_0based + 1 != row.displayed_genomic_start_1based as u64
                    || interval.end_0based_exclusive != row.displayed_genomic_end_1based as u64
                {
                    return Err(invalid("Regulatory annotation local/genomic mismatch"));
                }
                interval.start_0based = row.core_genomic_start_1based as u64 - 1;
                interval.end_0based_exclusive = row.core_genomic_end_1based as u64;
                rows.push(FragmentSelectionEvidence { evidence_id: format!("annotation:{}:{}", row.source_id, row.feature_id), kind: FragmentEvidenceKind::Annotation, label: row.feature_type.clone(), interval, source_sha256: source.sha256.clone(), source_id: row.source_id.clone(), available, may_seed_boundary: available, score: None, statement: "Provider annotation, not demonstrated reporter activity; complete core bounds retained.".into() });
            }
        }
        for track in &locus.regulatory_score_tracks {
            let available = track.state == gp::GeneLocusRegulatoryScoreState::Available;
            let binding = locus
                .sequence_binding
                .as_ref()
                .ok_or_else(|| invalid("Missing sequence binding"))?;
            if available
                && (track.input_sequence_id != locus.seq_id
                    || track.input_sequence_sha256.trim_start_matches("sha256:")
                        != binding.sequence_sha256.trim_start_matches("sha256:")
                    || track.assembly != locus.isoform_evidence.assembly)
            {
                return Err(invalid("TFBS track sequence/assembly mismatch"));
            }
            if available
                && binding.genome_anchor.as_ref().is_none_or(|anchor| {
                    !Self::chromosomes_match(&anchor.chromosome, &track.chromosome)
                        || anchor.start_1based != track.anchor_start_1based
                        || anchor.end_1based != track.anchor_end_1based
                })
            {
                return Err(invalid("TFBS track genome-anchor mismatch"));
            }
            for site in &track.sites {
                if !site.score.is_finite() {
                    return Err(invalid("Non-finite motif score"));
                }
                let (interval, _) = self.interval_and_projection_from_local(
                    &locus.seq_id,
                    site.local_start_0based as u64,
                    site.local_end_0based_exclusive as u64,
                    gp::GenomicRegionStrand::Unstranded,
                    None,
                )?;
                if interval.start_0based + 1 != site.genomic_start_1based as u64
                    || interval.end_0based_exclusive != site.genomic_end_1based as u64
                {
                    return Err(invalid("TFBS local/genomic coordinates disagree"));
                }
                rows.push(FragmentSelectionEvidence { evidence_id: format!("motif:{}:{}",track.track_id,site.site_id), kind: FragmentEvidenceKind::ModelSite, label: site.label.clone().unwrap_or_else(|| track.label.clone()), interval, source_sha256: source.sha256.clone(), source_id: track.track_id.clone(), available, may_seed_boundary: available, score: Some(site.score), statement: "Stored model-site prediction; no rescoring, occupancy confirmation or activity inference.".into() });
            }
        }
        for lane in locus.occupancy_groups.iter().flat_map(|g| &g.lanes) {
            let available = lane.state == gp::GeneLocusOccupancyLaneState::Available;
            if available
                && (lane
                    .source_sha256
                    .as_deref()
                    .is_none_or(|s| !Self::regulatory_fragment_valid_sha256(s))
                    || lane.source_assembly.as_deref()
                        != Some(locus.isoform_evidence.assembly.as_str()))
            {
                return Err(invalid(
                    "Occupancy lane lacks source hash or matching assembly",
                ));
            }
            for interval in &lane.lane.intervals {
                if interval.score.is_some_and(|score| !score.is_finite()) {
                    return Err(invalid("Non-finite occupancy signal"));
                }
                let (geometry, _) = self.interval_and_projection_from_local(
                    &locus.seq_id,
                    interval.local_start_1based.saturating_sub(1) as u64,
                    interval.local_end_1based as u64,
                    gp::GenomicRegionStrand::Unstranded,
                    None,
                )?;
                if geometry.start_0based + 1 != interval.genomic_start_1based as u64
                    || geometry.end_0based_exclusive != interval.genomic_end_1based as u64
                {
                    return Err(invalid("Occupancy local/genomic coordinates disagree"));
                }
                rows.push(FragmentSelectionEvidence { evidence_id: format!("coverage:{}:{}:{}",lane.lane.lane_id,geometry.start_0based,geometry.end_0based_exclusive), kind: FragmentEvidenceKind::RawCoverage, label: lane.lane.display_label.clone(), interval: geometry, source_sha256: lane.source_sha256.clone().unwrap_or_else(|| source.sha256.clone()), source_id: lane.lane.lane_id.clone(), available, may_seed_boundary: false, score: interval.score, statement: "Raw coverage only; not a peak, enrichment test or independent motif confirmation.".into() });
            }
        }
        if rows.len() > limit {
            return Err(invalid(
                "Evidence interval budget exceeded; supply a bounded locus report",
            ));
        }
        rows.sort_by(|a, b| a.evidence_id.cmp(&b.evidence_id));
        if rows
            .windows(2)
            .any(|w| w[0].evidence_id == w[1].evidence_id)
        {
            return Err(invalid("Duplicate evidence interval identity"));
        }
        Ok(rows)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gb_io::seq::{Feature, Location};
    use std::fs;
    use tempfile::{TempDir, tempdir};

    // Hand-crafted, deterministic locus/lanes, not experimental or provider data.
    struct Fixture {
        temp: TempDir,
        engine: GentleEngine,
        locus: gp::GeneLocusEvidenceDisplayReport,
        request: FragmentSelectionRequest,
    }

    impl Fixture {
        fn bind(&mut self) {
            let bytes = serde_json::to_vec(&self.locus).unwrap();
            let path = self.temp.path().join("locus.json");
            fs::write(&path, &bytes).unwrap();
            self.request.locus.path = path.to_string_lossy().into();
            self.request.locus.sha256 = crate::digest_utils::sha256_prefixed_bytes(&bytes);
        }
        fn plan(&self) -> FragmentSelectionReport {
            self.engine
                .plan_reporter_fragment_selection(self.request.clone())
                .unwrap()
        }
    }

    fn fixture(reverse_anchor: bool, reverse_transcript: bool) -> Fixture {
        let mut engine = GentleEngine::default();
        let tss = if reverse_transcript { 1399 } else { 1100 };
        let mut sequence = "ACGT".repeat(750);
        sequence.replace_range(tss - 115..tss - 109, "AAGCTT");
        if reverse_transcript {
            sequence.replace_range(1297..1300, "CAT");
        } else {
            sequence.replace_range(1200..1203, "ATG");
        }
        let mut dna = DNAsequence::from_sequence(&sequence).unwrap();
        for (kind, start, end) in [("mRNA", 1100, 1400), ("CDS", 1200, 1300)] {
            let loc = Location::simple_range(start, end);
            dna.features_mut().push(Feature {
                kind: kind.into(),
                location: if reverse_transcript {
                    Location::Complement(Box::new(loc))
                } else {
                    loc
                },
                qualifiers: vec![("transcript_id".into(), Some("TX_A.1".into()))],
            });
        }
        GentleEngine::prepare_sequence(&mut dna);
        engine.state.sequences.insert("source".into(), dna);
        engine.state.metadata.insert(PROVENANCE_METADATA_KEY.into(),serde_json::json!({"genome_extractions":[{
            "seq_id":"source","genome_id":"GRCh38","chromosome":"chr7","start_1based":10001,"end_1based":13000,
            "anchor_strand":if reverse_anchor {"-"} else {"+"},"anchor_verified":true,"recorded_at_unix_ms":0
        }]}));
        let anchor = engine.sequence_genome_anchor_summary("source").unwrap();
        let binding =
            crate::locus_report::sequence_binding(&engine.state.sequences["source"], Some(&anchor));
        let geometry = |start: usize, end: usize| {
            engine
                .interval_and_projection_from_local(
                    "source",
                    start as u64,
                    end as u64,
                    gp::GenomicRegionStrand::Unstranded,
                    None,
                )
                .unwrap()
                .0
        };
        let promoter = geometry(tss - 100, tss + 100);
        let exons = geometry(1100, 1400);
        let genomic_strand = if reverse_anchor ^ reverse_transcript {
            "-"
        } else {
            "+"
        };
        let mut locus = gp::GeneLocusEvidenceDisplayReport {
            schema: gp::GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA.into(),
            seq_id: "source".into(),
            panel_id: "synthetic_locus".into(),
            sequence_binding: Some(binding.clone()),
            locus_local_start_1based: 1,
            locus_local_end_1based: 3000,
            ..Default::default()
        };
        locus.isoform_evidence = gp::GeneIsoformEvidenceReport {
            assembly: "GRCh38".into(),
            annotation_release: Some("synthetic-1".into()),
            transcripts: vec![gp::GeneIsoformTranscriptRow {
                transcript_id: "TX_A.1".into(),
                strand: genomic_strand.into(),
                exon_family_ids_5_to_3: vec!["exonA".into()],
                ..Default::default()
            }],
            exon_families: vec![gp::GeneIsoformExonFamilyRow {
                exon_family_id: "exonA".into(),
                local_start_1based: 1101,
                local_end_1based: 1400,
                start_1based: exons.start_0based as usize + 1,
                end_1based: exons.end_0based_exclusive as usize,
                ..Default::default()
            }],
            ..Default::default()
        };
        locus.ensembl_regulation = Some(serde_json::from_value(serde_json::json!({
            "availability":"available","source_binding":{"content_identity_verified":true,"source":{"assembly_name":"GRCh38"}},
            "rows":[{"source_id":"reg","feature_id":"promoterA","feature_type":"Promoter","assembly_name":"GRCh38",
                "core_genomic_start_1based":promoter.start_0based+1,"core_genomic_end_1based":promoter.end_0based_exclusive,
                "displayed_genomic_start_1based":promoter.start_0based+1,"displayed_genomic_end_1based":promoter.end_0based_exclusive,
                "displayed_local_start_1based":tss-99,"displayed_local_end_1based":tss+100}]
        })).unwrap());
        let site = geometry(tss - 70, tss - 50);
        locus
            .regulatory_score_tracks
            .push(gp::GeneLocusRegulatoryScoreTrack {
                track_id: "model".into(),
                label: "synthetic PWM".into(),
                state: gp::GeneLocusRegulatoryScoreState::Available,
                input_sequence_id: "source".into(),
                input_sequence_sha256: binding.sequence_sha256,
                assembly: "GRCh38".into(),
                chromosome: "chr7".into(),
                anchor_start_1based: 10001,
                anchor_end_1based: 13000,
                sites: vec![gp::GeneLocusRegulatoryScoreSite {
                    site_id: "siteA".into(),
                    local_start_0based: (tss - 70) as usize,
                    local_end_0based_exclusive: (tss - 50) as usize,
                    genomic_start_1based: site.start_0based as usize + 1,
                    genomic_end_1based: site.end_0based_exclusive as usize,
                    score: 7.5,
                    ..Default::default()
                }],
                ..Default::default()
            });
        locus.occupancy_groups.push(gp::GeneLocusOccupancyGroup {
            lanes: [("sample", 5.0), ("control", 1.0)]
                .into_iter()
                .map(|(id, score)| gp::GeneLocusOccupancyLane {
                    state: gp::GeneLocusOccupancyLaneState::Available,
                    source_sha256: Some(crate::digest_utils::sha256_prefixed_bytes(id.as_bytes())),
                    source_assembly: Some("GRCh38".into()),
                    cell_line_label: Some("synthetic cells".into()),
                    lane: gp::GeneIsoformOccupancyLane {
                        lane_id: id.into(),
                        interval_count: 1,
                        intervals: vec![gp::GeneIsoformOccupancyInterval {
                            local_start_1based: 1,
                            local_end_1based: 3000,
                            genomic_start_1based: 10001,
                            genomic_end_1based: 13000,
                            score: Some(score),
                            ..Default::default()
                        }],
                        ..Default::default()
                    },
                    ..Default::default()
                })
                .collect(),
            ..Default::default()
        });
        let request = FragmentSelectionRequest {
            schema: FRAGMENT_SELECTION_REQUEST_SCHEMA.into(),
            locus: FragmentSelectionLocus {
                path: String::new(),
                sha256: String::new(),
                panel_id: locus.panel_id.clone(),
                annotation_release: "synthetic-1".into(),
            },
            anchors: vec![FragmentSelectionAnchor {
                transcript_id: "TX_A.1".into(),
                seed_evidence_ids: vec!["annotation:reg:promoterA".into()],
                ..Default::default()
            }],
            policy: Default::default(),
            adjustments: vec![],
            signal_comparisons: vec![],
            vector: None,
        };
        let mut f = Fixture {
            temp: tempdir().unwrap(),
            engine,
            locus,
            request,
        };
        f.engine = GentleEngine::from_state(f.engine.snapshot().clone());
        f.bind();
        f
    }

    fn comparison() -> FragmentSignalComparison {
        FragmentSignalComparison {
            comparison_id: "sample-control".into(),
            sample_lane_id: "sample".into(),
            control_lane_id: "control".into(),
            cell_line: "synthetic cells".into(),
            sample_replicate_id: "s1".into(),
            control_replicate_id: "c1".into(),
            units: "synthetic per-base signal".into(),
            minimum_mean_difference: 2.0,
            missing_policy: FragmentSignalMissingPolicy::RequireComplete,
        }
    }

    fn adjustment(upstream: i64, downstream: i64) -> FragmentBoundaryAdjustment {
        FragmentBoundaryAdjustment {
            adjustment_id: "inspect-flank".into(),
            transcript_id: "TX_A.1".into(),
            seed_id: "annotation:reg:promoterA".into(),
            upstream_delta_bp: upstream,
            downstream_delta_bp: downstream,
            reason: FragmentBoundaryReason::TfbsContext,
            explanation: "Inspect additional model context".into(),
            evidence_ids: vec!["motif:model:siteA".into()],
        }
    }

    #[test]
    fn reporter_fragment_selection_is_bound_read_only_deterministic_and_importable() {
        let mut f = fixture(false, false);
        let before = digest(&f.engine.snapshot()).unwrap();
        let report = f.plan();
        assert_eq!(report.schema, FRAGMENT_SELECTION_REPORT_SCHEMA);
        assert_eq!(
            report.anchors[0].search_envelope.end_0based_exclusive
                - report.anchors[0].search_envelope.start_0based,
            1001
        );
        assert_eq!(report.candidates[0].ranking.length_bp, 240);
        assert_eq!(
            report.candidates[0].ranking.bisected_feature_count, 0,
            "signal-bin edges cannot penalize geometry"
        );
        assert!(report.anchors[0].gene_structure.is_some());
        assert_eq!(
            serde_json::to_vec(&report).unwrap(),
            serde_json::to_vec(&f.plan()).unwrap()
        );
        assert_eq!(before, digest(&f.engine.snapshot()).unwrap());
        let command = crate::engine_shell::parse_shell_tokens(&[
            "promoters".into(),
            "fragment-candidates".into(),
            serde_json::to_string(&f.request).unwrap(),
        ])
        .unwrap();
        let shell = crate::engine_shell::execute_shell_command(&mut f.engine, &command).unwrap();
        assert!(!shell.state_changed);
        assert_eq!(
            shell.output["result"],
            serde_json::to_value(&report).unwrap()
        );
        assert_eq!(before, digest(&f.engine.snapshot()).unwrap());
        let path = f.temp.path().join("regions.json");
        let output = f.temp.path().join("proposal.json");
        let result = f
            .engine
            .apply(Operation::PlanEvidenceGuidedFragmentCandidates {
                request: Box::new(f.request.clone()),
                path: Some(output.to_string_lossy().into()),
            })
            .unwrap();
        let written: serde_json::Value =
            serde_json::from_slice(&fs::read(output).unwrap()).unwrap();
        assert_eq!(
            written,
            serde_json::to_value(result.reporter_fragment_selection.as_ref().unwrap()).unwrap()
        );
        let input_bytes = fs::read(&f.request.locus.path).unwrap();
        assert!(
            f.engine
                .apply(Operation::PlanEvidenceGuidedFragmentCandidates {
                    request: Box::new(f.request.clone()),
                    path: Some(f.request.locus.path.clone())
                })
                .is_err()
        );
        assert_eq!(input_bytes, fs::read(&f.request.locus.path).unwrap());
        assert_eq!(before, digest(&f.engine.snapshot()).unwrap());
        fs::write(
            &path,
            serde_json::to_vec(&report.proposed_region_set).unwrap(),
        )
        .unwrap();
        f.engine
            .apply(Operation::ImportGenomicRegionSet {
                request: gp::GenomicRegionImportRequest {
                    path: path.to_string_lossy().into(),
                    ..Default::default()
                },
            })
            .unwrap();
    }

    #[test]
    fn reporter_fragment_selection_adjusts_both_ends_in_transcript_orientation() {
        for reverse_anchor in [false, true] {
            for reverse_transcript in [false, true] {
                let mut f = fixture(reverse_anchor, reverse_transcript);
                f.request.adjustments.push(adjustment(50, -10));
                let report = f.plan();
                let base = report
                    .candidates
                    .iter()
                    .find(|c| c.parent_candidate_id.is_none())
                    .unwrap();
                let child = report
                    .candidates
                    .iter()
                    .find(|c| c.parent_candidate_id.is_some())
                    .unwrap();
                assert_eq!(
                    child.parent_candidate_id.as_deref(),
                    Some(base.candidate_id.as_str())
                );
                assert!(child.blockers.is_empty());
                assert_eq!(child.ranking.length_bp, 280);
                let a = &base.region.interval;
                let b = &child.region.interval;
                if reverse_anchor ^ reverse_transcript {
                    assert_eq!(b.start_0based, a.start_0based + 10);
                    assert_eq!(b.end_0based_exclusive, a.end_0based_exclusive + 50);
                } else {
                    assert_eq!(b.start_0based + 50, a.start_0based);
                    assert_eq!(b.end_0based_exclusive + 10, a.end_0based_exclusive);
                }
                let projection = child.region.local_projection.as_ref().unwrap();
                let raw = &f.engine.state.sequences["source"].get_forward_string()[projection
                    .local_start_0based
                    as usize
                    ..projection.local_end_0based_exclusive as usize];
                assert_eq!(
                    child.sequence_5prime_to_3prime,
                    if reverse_transcript {
                        GentleEngine::reverse_complement(raw)
                    } else {
                        raw.into()
                    }
                );
            }
        }
    }

    #[test]
    fn reporter_fragment_selection_keeps_long_context_and_blocks_destructive_shortening() {
        let mut f = fixture(false, false);
        f.request.policy.preferred_length_bp = 100;
        f.request.adjustments.push(adjustment(-60, 0));
        let r = f.plan();
        assert!(
            r.candidates
                .iter()
                .any(|c| c.blockers.is_empty() && c.ranking.length_bp > 100)
        );
        let bad = r
            .candidates
            .iter()
            .find(|c| c.variant.starts_with("human:"))
            .unwrap();
        assert!(!bad.ranking.preserves_required_context);
        assert!(
            bad.blockers
                .contains(&"required_evidence_or_promoter_context_removed".into())
        );
        assert!(
            !r.proposed_region_set
                .regions
                .iter()
                .any(|c| c.region_id == bad.candidate_id)
        );
        assert!(
            bad.bisected_evidence_ids
                .contains(&"annotation:reg:promoterA".into())
        );
    }

    #[test]
    fn reporter_fragment_selection_compacts_only_optional_padding() {
        let mut f = fixture(false, false);
        f.request.policy.preferred_length_bp = 220;
        let r = f.plan();
        let compact = r
            .candidates
            .iter()
            .find(|c| c.variant == "compact_without_optional_flanks")
            .unwrap();
        let original = r
            .candidates
            .iter()
            .find(|c| c.parent_candidate_id.is_none())
            .unwrap();
        assert_eq!(compact.ranking.length_bp, 200);
        assert_eq!(original.ranking.length_bp, 240);
        assert_eq!(
            compact.parent_candidate_id.as_deref(),
            Some(original.candidate_id.as_str())
        );
        assert!(compact.blockers.is_empty());
        assert!(original.blockers.is_empty());
        for id in ["annotation:reg:promoterA", "motif:model:siteA"] {
            assert!(compact.retained_evidence_ids.iter().any(|e| e == id));
            assert!(!compact.bisected_evidence_ids.iter().any(|e| e == id));
        }
    }

    #[test]
    fn reporter_fragment_selection_preserves_tss_memberships_and_request_order_independence() {
        let mut f = fixture(false, false);
        for (id, reverse) in [("TX_B.1", false), ("TX_C.1", true)] {
            let location = Location::simple_range(1100, 1400);
            f.engine
                .state
                .sequences
                .get_mut("source")
                .unwrap()
                .features_mut()
                .push(Feature {
                    kind: "mRNA".into(),
                    location: if reverse {
                        Location::Complement(Box::new(location))
                    } else {
                        location
                    },
                    qualifiers: vec![("transcript_id".into(), Some(id.into()))],
                });
            let mut row = f.locus.isoform_evidence.transcripts[0].clone();
            row.transcript_id = id.into();
            row.strand = if reverse { "-" } else { "+" }.into();
            f.locus.isoform_evidence.transcripts.push(row);
        }
        let mut second_anchor = f.request.anchors[0].clone();
        second_anchor.transcript_id = "TX_B.1".into();
        f.request.anchors.push(second_anchor);
        f.request.adjustments.push(adjustment(40, 0));
        let mut second_adjustment = adjustment(0, 40);
        second_adjustment.adjustment_id = "inspect-other-flank".into();
        f.request.adjustments.push(second_adjustment);
        f.bind();
        let r = f.plan();
        for anchor in &r.anchors {
            assert_eq!(anchor.opposite_strand_transcript_ids, ["TX_C.1"]);
            assert!(
                anchor
                    .findings
                    .iter()
                    .any(|s| s.contains("Opposite-strand"))
            );
        }
        for candidate in &r.candidates {
            assert_eq!(candidate.member_transcript_ids, ["TX_A.1", "TX_B.1"]);
        }
        f.request.anchors.reverse();
        f.request.adjustments.reverse();
        assert_eq!(
            serde_json::to_vec(&r).unwrap(),
            serde_json::to_vec(&f.plan()).unwrap()
        );
    }

    #[test]
    fn reporter_fragment_selection_raw_coverage_never_seeds_or_claims_enrichment() {
        let mut f = fixture(false, false);
        f.locus.regulatory_score_tracks.clear();
        f.locus.ensembl_regulation = None;
        f.request.anchors[0].seed_evidence_ids.clear();
        f.bind();
        let r = f.plan();
        assert!(r.candidates.is_empty());
        assert!(r.fixed_comparisons.is_empty());
        assert!(
            r.evidence
                .iter()
                .all(|e| e.kind == FragmentEvidenceKind::RawCoverage && !e.may_seed_boundary)
        );
    }

    #[test]
    fn reporter_fragment_selection_comparisons_use_fixed_envelopes_and_explicit_missingness() {
        let mut f = fixture(false, false);
        f.request.signal_comparisons.push(comparison());
        f.request.adjustments.push(adjustment(50, -10));
        let r = f.plan();
        assert_eq!(r.fixed_comparisons.len(), 1);
        assert_eq!(r.fixed_comparisons[0].mean_difference, Some(4.0));
        assert_eq!(r.fixed_comparisons[0].passes_descriptive_rule, Some(true));
        assert_ne!(
            r.fixed_comparisons[0].measurement_window,
            r.candidates[0].region.interval
        );
        f.locus.occupancy_groups[0].lanes[1].lane.interval_count = 2;
        f.bind();
        let r = f.plan();
        assert_eq!(r.fixed_comparisons[0].passes_descriptive_rule, None);
        f.locus.occupancy_groups[0].lanes[1].lane.interval_count = 1;
        f.locus.occupancy_groups[0].lanes[1].source_sha256 =
            f.locus.occupancy_groups[0].lanes[0].source_sha256.clone();
        f.bind();
        assert!(
            f.engine
                .plan_reporter_fragment_selection(f.request.clone())
                .unwrap_err()
                .message
                .contains("same source")
        );
    }

    #[test]
    fn reporter_fragment_selection_rejects_stale_report_annotation_and_shifted_tracks() {
        let mut f = fixture(false, false);
        fs::write(&f.request.locus.path, b"{}").unwrap();
        assert!(
            f.engine
                .plan_reporter_fragment_selection(f.request.clone())
                .is_err()
        );
        f.bind();
        f.request.locus.annotation_release = "other".into();
        assert!(
            f.engine
                .plan_reporter_fragment_selection(f.request.clone())
                .is_err()
        );
        f.request.locus.annotation_release = "synthetic-1".into();
        f.locus.regulatory_score_tracks[0].sites[0].genomic_start_1based += 1;
        f.bind();
        assert!(
            f.engine
                .plan_reporter_fragment_selection(f.request.clone())
                .unwrap_err()
                .message
                .contains("TFBS local/genomic")
        );
        f.locus.regulatory_score_tracks[0].sites[0].genomic_start_1based -= 1;
        f.locus.occupancy_groups[0].lanes[0].lane.intervals[0].genomic_end_1based -= 1;
        f.bind();
        assert!(
            f.engine
                .plan_reporter_fragment_selection(f.request.clone())
                .unwrap_err()
                .message
                .contains("Occupancy local/genomic")
        );
        f.locus.occupancy_groups[0].lanes[0].lane.intervals[0].genomic_end_1based += 1;
        f.bind();
        f.engine
            .state
            .sequences
            .get_mut("source")
            .unwrap()
            .features_mut()[0]
            .location = Location::simple_range(1100, 1399);
        assert!(
            f.engine
                .plan_reporter_fragment_selection(f.request.clone())
                .unwrap_err()
                .message
                .contains("geometry/strand changed")
        );
    }

    #[test]
    fn reporter_fragment_selection_bounds_fail_explicitly() {
        let mut f = fixture(false, false);
        f.request.policy.maximum_candidates = 1;
        f.request.adjustments.push(adjustment(50, 0));
        assert!(
            f.engine
                .plan_reporter_fragment_selection(f.request.clone())
                .unwrap_err()
                .message
                .contains("Candidate budget")
        );
        f.request.policy.maximum_candidates = 128;
        f.request.adjustments[0].upstream_delta_bp = 201;
        assert!(
            f.engine
                .plan_reporter_fragment_selection(f.request.clone())
                .unwrap_err()
                .message
                .contains("shift limit")
        );
    }

    fn add_vector(f: &mut Fixture) {
        let fixture_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("test_files/fixtures/reporter_vectors/synthetic_mcs_backbone.gb");
        let gb = gb_io::reader::parse_file(&fixture_path).unwrap().remove(0);
        let mut vector = DNAsequence::from_genbank_seq(gb);
        GentleEngine::prepare_sequence(&mut vector);
        f.engine.state.sequences.insert("vector".into(), vector);
        let catalog = serde_json::json!({"Synthetic panel vector":{
            "description":"Hand-crafted layout-only validation fixture","sequence_local":fixture_path,"annotations_local":fixture_path,
            "usable_as_empty_backbone":true,"helper_kind":"plasmid_vector","sequence_expectation":{
                "schema":crate::genomes::HELPER_VECTOR_SEQUENCE_EXPECTATION_SCHEMA,"provider":"GENtle tests","product_name":"synthetic MCS backbone",
                "catalog_number":"SYNTH-MCS-1","accession_version":"GENTLE_SYNTHETIC_MCS.1","expected_length_bp":240,"expected_topology":"circular",
                "required_features":[
                    {"id":"multiple_cloning_region","feature_kinds":["misc_feature"],"qualifier_terms":["multiple cloning site region"],"expected_start_1based":1,"expected_end_1based":70},
                    {"id":"luc2","feature_kinds":["CDS"],"qualifier_terms":["luciferase luc2 marker"],"expected_start_1based":100,"expected_end_1based":180}
                ],"restriction_site_equivalences":[],"provenance":[{"source_id":"synthetic-test-fixture","source_url":"test_files/fixtures/reporter_vectors/synthetic_mcs_backbone.gb","asserted_on":"2026-09-14","note":"Repository-owned deterministic layout fixture, not a commercial vector."}]
            }
        }});
        let path = f.temp.path().join("vectors.json");
        fs::write(&path, serde_json::to_vec(&catalog).unwrap()).unwrap();
        f.request.vector = Some(FragmentSelectionVector {
            seq_id: "vector".into(),
            catalog_id: "Synthetic panel vector".into(),
            helper_catalog_path: Some(path.to_string_lossy().into()),
            suggest_restriction_adjustments: true,
        });
    }

    #[test]
    fn reporter_fragment_selection_verified_mcs_drives_safe_padding_trims() {
        let mut f = fixture(false, false);
        add_vector(&mut f);
        let before = digest(&f.engine.snapshot()).unwrap();
        let r = f.plan();
        let vector = r.vector.as_ref().unwrap();
        assert_eq!(
            serde_json::to_vec(&r).unwrap(),
            serde_json::to_vec(&f.plan()).unwrap()
        );
        assert_eq!(
            vector.validation.status,
            ReporterVectorValidationStatus::Verified
        );
        assert_eq!(
            (vector.mcs_start_0based, vector.mcs_end_0based_exclusive),
            (0, 70)
        );
        assert!(vector.source_sites.iter().any(|s| s.enzyme == "HindIII"));
        let trim = r
            .candidates
            .iter()
            .find(|c| c.variant.starts_with("restriction:HindIII:"))
            .expect("padding site should be avoidable");
        assert!(trim.blockers.is_empty());
        assert!(trim.ranking.preserves_required_context);
        assert!(trim.cloning.as_ref().unwrap().selected_pair.is_some());
        for c in &r.candidates {
            let cloning = c.cloning.as_ref().unwrap();
            assert!(
                cloning
                    .pair_evaluations
                    .iter()
                    .all(|p| p.pair.forward_cut_position_0based < 70
                        && p.pair.reverse_cut_position_0based < 70)
            );
            assert!(cloning.warnings.iter().any(|w| w.contains("added by PCR")));
        }
        assert_eq!(before, digest(&f.engine.snapshot()).unwrap());
        f.engine
            .state
            .sequences
            .get_mut("vector")
            .unwrap()
            .features_mut()
            .retain(|feature| {
                !feature
                    .kind
                    .to_string()
                    .eq_ignore_ascii_case("misc_feature")
            });
        assert!(
            f.engine
                .plan_reporter_fragment_selection(f.request.clone())
                .is_err(),
            "missing MCS must not quietly use vector-wide sites"
        );
    }

    #[test]
    fn reporter_fragment_selection_human_roi_can_be_revisited_without_changing_original() {
        let mut f = fixture(false, false);
        let old = f.plan().candidates.remove(0).region;
        let original = serde_json::to_vec(&old).unwrap();
        f.request.anchors[0].seed_evidence_ids.clear();
        f.request.anchors[0].seed_region = Some(old.clone());
        let mut edit = adjustment(70, -10);
        edit.seed_id = old.region_id.clone();
        edit.reason = FragmentBoundaryReason::CutrunContext;
        f.request.adjustments.push(edit);
        let r = f.plan();
        assert_eq!(r.candidates.len(), 2);
        assert_eq!(
            original,
            serde_json::to_vec(r.request.anchors[0].seed_region.as_ref().unwrap()).unwrap()
        );
        assert!(
            r.candidates
                .iter()
                .any(|c| c.ranking.length_bp == 300 && c.blockers.is_empty())
        );
    }
}
