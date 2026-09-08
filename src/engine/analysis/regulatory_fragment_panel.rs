//! Deterministic planning for bounded regulatory-fragment reporter contrasts.
//!
//! This slice consumes exact persisted genomic ROIs and an exact reporter
//! vector. It plans inserts and comparisons without writing files or mutating
//! live project state. Typed locus evidence stays source-bound; exact product
//! materialization is a separate, digest-approved atomic transition.

use super::*;
use gentle_protocol as gp;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs;

#[path = "regulatory_fragment_panel/external_evidence.rs"]
mod external_evidence;
#[path = "regulatory_fragment_panel/materialization.rs"]
mod materialization;

#[derive(Debug, Clone)]
struct ResolvedRegulatoryFragment {
    binding: RegulatoryFragmentBinding,
    projection: gp::GenomicRegionLocalProjection,
    canonical_sequence: String,
    warnings: Vec<RegulatoryFragmentFinding>,
}

#[derive(Debug, Clone)]
struct CandidateRegulatoryConstruct {
    member: RegulatoryFragmentPanelMember,
    sort_key: (usize, usize, String),
    mandatory: bool,
    eligible: bool,
    questions: BTreeSet<RegulatoryFragmentQuestion>,
}

#[derive(Debug, Clone)]
struct CandidateContrast {
    contrast: RegulatoryFragmentContrast,
}

impl GentleEngine {
    /// Build a content-addressed, read-only panel plan from exact persisted ROIs.
    pub fn plan_regulatory_fragment_panel(
        &self,
        request: RegulatoryFragmentPanelRequest,
    ) -> Result<RegulatoryFragmentPanelPlan, EngineError> {
        let request = Self::normalize_regulatory_fragment_panel_request(request)?;
        let (resolved, mut validation_warnings) =
            self.resolve_regulatory_fragment_bindings(&request)?;
        let (vector_validation, vector_context) =
            self.resolve_regulatory_fragment_vector_context(&request)?;
        let request_sha256 =
            Self::regulatory_fragment_value_sha256(&request, "regulatory-fragment panel request")?;
        let source_state_sha256 = Self::regulatory_fragment_value_sha256(
            &serde_json::json!({
                "regions": resolved.iter().map(|fragment| serde_json::json!({
                    "fragment_id": fragment.binding.fragment_id,
                    "region_set_id": fragment.binding.region_set_id,
                    "region_set_content_sha256": fragment.binding.region_set_content_sha256,
                    "region_identity_sha256": fragment.binding.region.identity_sha256,
                    "region_content_sha256": fragment.binding.region.content_sha256,
                    "source_seq_id": fragment.projection.seq_id,
                    "source_sequence_sha256": fragment.projection.sequence_sha256,
                })).collect::<Vec<_>>(),
                "vector_context": vector_context,
                "evidence_bindings": request.evidence_bindings,
            }),
            "regulatory-fragment source state",
        )?;

        let by_id = resolved
            .iter()
            .map(|fragment| (fragment.binding.fragment_id.as_str(), fragment))
            .collect::<HashMap<_, _>>();
        let similarity = self.regulatory_fragment_candidate_partner_similarity(&request, &by_id)?;
        let highly_similar = similarity
            .as_ref()
            .is_some_and(|audit| audit.highly_similar);

        let (mut candidates, mut omitted) = self.regulatory_fragment_construct_candidates(
            &request,
            &by_id,
            &vector_context,
            highly_similar,
        )?;
        let contrast_candidates =
            Self::regulatory_fragment_contrast_candidates(&request, &candidates);
        let selected_ids =
            Self::select_regulatory_fragment_members(&request, &candidates, &contrast_candidates);
        let selected_id_set = selected_ids.iter().cloned().collect::<BTreeSet<_>>();
        let covered_questions = Self::covered_regulatory_fragment_questions(
            &request.questions,
            &contrast_candidates,
            &selected_id_set,
        );
        let uncovered_questions = request
            .questions
            .iter()
            .copied()
            .filter(|question| !covered_questions.contains(question))
            .collect::<Vec<_>>();

        let mut contrasts = contrast_candidates
            .iter()
            .filter(|candidate| {
                selected_id_set.contains(&candidate.contrast.left_member_id)
                    && selected_id_set.contains(&candidate.contrast.right_member_id)
            })
            .map(|candidate| candidate.contrast.clone())
            .collect::<Vec<_>>();
        contrasts.sort_by(|left, right| left.contrast_id.cmp(&right.contrast_id));

        let planning_label = if !uncovered_questions.is_empty() {
            RegulatoryFragmentPlanningLabel::UnresolvedExperimentalComparisonRequired
        } else if highly_similar {
            RegulatoryFragmentPlanningLabel::ContextOrGeometryConfounded
        } else if request.questions.iter().any(|question| {
            matches!(
                question,
                RegulatoryFragmentQuestion::PartnerDependence
                    | RegulatoryFragmentQuestion::OrderDependence
                    | RegulatoryFragmentQuestion::OrientationDependence
                    | RegulatoryFragmentQuestion::SpacingDependence
            )
        }) {
            RegulatoryFragmentPlanningLabel::PartnerDependenceTestableHypothesis
        } else {
            RegulatoryFragmentPlanningLabel::StandaloneTestableCandidate
        };

        let contrast_by_member = Self::regulatory_fragment_contrasts_by_member(&contrasts);
        let mut members = vec![];
        for candidate in &mut candidates {
            if selected_id_set.contains(&candidate.member.member_id) {
                candidate.member.planning_label = planning_label;
                candidate.member.contrast_ids = contrast_by_member
                    .get(&candidate.member.member_id)
                    .cloned()
                    .unwrap_or_default();
                Self::attach_regulatory_fragment_inclusion_reasons(candidate, &contrasts, &request);
                members.push(candidate.member.clone());
            } else {
                omitted.push(RegulatoryFragmentOmittedVariant {
                    member_id: candidate.member.member_id.clone(),
                    construct_kind: candidate.member.construct_kind,
                    reason: if candidate.eligible && uncovered_questions.is_empty() {
                        RegulatoryFragmentOmissionReasonKind::RedundantForQuestionCoverage
                    } else if candidate.eligible {
                        RegulatoryFragmentOmissionReasonKind::PanelMemberBound
                    } else {
                        RegulatoryFragmentOmissionReasonKind::ConstructLengthBound
                    },
                    affected_questions: candidate.questions.iter().copied().collect(),
                    detail: if candidate.eligible && uncovered_questions.is_empty() {
                        "The deterministic minimal cover answered every requested question without this construct."
                            .to_string()
                    } else if candidate.eligible {
                        format!(
                            "The max_panel_members={} bound prevented retaining this construct in the best deterministic partial cover.",
                            request.max_panel_members
                        )
                    } else {
                        format!(
                            "The {} bp insert exceeds max_construct_length_bp={}.",
                            candidate.member.insert_length_bp, request.max_construct_length_bp
                        )
                    },
                });
            }
        }
        omitted.sort_by(|left, right| left.member_id.cmp(&right.member_id));

        let cloning_strategy =
            self.populate_regulatory_fragment_cloning_feasibility(&request, &mut members)?;
        let evidence_dimensions = self.populate_regulatory_fragment_evidence_dimensions(
            &request,
            &resolved,
            &members,
            cloning_strategy.as_ref(),
        )?;
        let mut blockers = vec![];
        if !uncovered_questions.is_empty() {
            blockers.push(RegulatoryFragmentFinding {
                code: "panel_member_bound_leaves_questions_uncovered".to_string(),
                subject_ids: uncovered_questions
                    .iter()
                    .map(|question| Self::regulatory_fragment_question_token(*question).to_string())
                    .collect(),
                detail: format!(
                    "The bounded deterministic cover cannot represent every requested comparison within {} panel members.",
                    request.max_panel_members
                ),
            });
        }
        if highly_similar {
            validation_warnings.push(RegulatoryFragmentFinding {
                code: "candidate_partner_high_similarity".to_string(),
                subject_ids: similarity
                    .iter()
                    .flat_map(|audit| {
                        [audit.left_fragment_id.clone(), audit.right_fragment_id.clone()]
                    })
                    .collect(),
                detail: "Candidate and partner sequences are highly similar under the declared global-alignment threshold; sequence context or geometry may confound interpretation."
                    .to_string(),
            });
        }
        let mut plan = RegulatoryFragmentPanelPlan {
            schema: REGULATORY_FRAGMENT_PANEL_PLAN_SCHEMA.to_string(),
            plan_id: request.plan_id.clone(),
            request_sha256,
            source_state_sha256,
            request,
            vector_validation,
            vector_context,
            planning_label,
            members,
            contrasts,
            omitted_variants: omitted,
            evidence_dimensions,
            candidate_partner_similarity: similarity,
            cloning_strategy,
            uncovered_questions,
            approval_required: true,
            materialization_supported: true,
            blockers,
            warnings: validation_warnings,
            nonclaims: vec![
                "This plan defines testable reporter comparisons; it does not establish regulatory sufficiency, enhancer or silencer activity, or partner dependence."
                    .to_string(),
                "Sequence similarity, annotations, occupancy evidence, and cloning feasibility are not proof of reporter activity or causal regulation."
                    .to_string(),
                "Exact ordered design constructs require a separate regulatory-fragment materialization proposal and approval; they are not simulated cloning products or evidence of experimental success."
                    .to_string(),
            ],
            ..RegulatoryFragmentPanelPlan::default()
        };
        plan.nonclaims
            .extend(plan.request.scientific_caveats.iter().cloned());
        plan.proposal_digest = Self::regulatory_fragment_panel_proposal_digest(&plan)?;
        Ok(plan)
    }

    /// Validate exact approval and current source state without materializing anything.
    pub fn validate_regulatory_fragment_panel_approval(
        &self,
        plan: &RegulatoryFragmentPanelPlan,
        approval_digest: &str,
    ) -> Result<(), EngineError> {
        Self::validate_regulatory_fragment_panel_document(plan)?;
        if approval_digest.trim() != plan.proposal_digest {
            return Err(Self::regulatory_fragment_error(
                "approval_digest_mismatch",
                [&plan.plan_id],
                "Approval must equal the exact current proposal_digest.",
            ));
        }
        let fresh = self.plan_regulatory_fragment_panel(plan.request.clone())?;
        if fresh.proposal_digest != plan.proposal_digest {
            return Err(Self::regulatory_fragment_error(
                "proposal_source_state_stale",
                [&plan.plan_id],
                "Current ROI, source-sequence, vector, evidence, or panel ordering state produces a different proposal digest.",
            ));
        }
        Ok(())
    }

    pub(crate) fn validate_regulatory_fragment_panel_document(
        plan: &RegulatoryFragmentPanelPlan,
    ) -> Result<(), EngineError> {
        if plan.schema != REGULATORY_FRAGMENT_PANEL_PLAN_SCHEMA {
            return Err(Self::regulatory_fragment_error(
                "unsupported_plan_schema",
                [&plan.plan_id],
                format!(
                    "Expected schema '{}', observed '{}'.",
                    REGULATORY_FRAGMENT_PANEL_PLAN_SCHEMA, plan.schema
                ),
            ));
        }
        let embedded = Self::regulatory_fragment_panel_proposal_digest(plan)?;
        if embedded != plan.proposal_digest {
            return Err(Self::regulatory_fragment_error(
                "proposal_content_digest_mismatch",
                [&plan.plan_id],
                "The supplied plan content no longer matches its embedded proposal_digest.",
            ));
        }
        Ok(())
    }

    fn normalize_regulatory_fragment_panel_request(
        mut request: RegulatoryFragmentPanelRequest,
    ) -> Result<RegulatoryFragmentPanelRequest, EngineError> {
        if request.schema.trim().is_empty() {
            request.schema = REGULATORY_FRAGMENT_PANEL_REQUEST_SCHEMA.to_string();
        }
        if request.schema != REGULATORY_FRAGMENT_PANEL_REQUEST_SCHEMA {
            return Err(Self::regulatory_fragment_error(
                "unsupported_request_schema",
                [&request.plan_id],
                format!(
                    "Expected schema '{}', observed '{}'.",
                    REGULATORY_FRAGMENT_PANEL_REQUEST_SCHEMA, request.schema
                ),
            ));
        }
        request.plan_id = Self::regulatory_fragment_slug(&request.plan_id);
        request.vector_seq_id = request.vector_seq_id.trim().to_string();
        request.vector_catalog_id = request.vector_catalog_id.trim().to_string();
        request.insertion_context_feature_id = request
            .insertion_context_feature_id
            .trim()
            .to_ascii_lowercase();
        if request.plan_id.is_empty()
            || request.vector_seq_id.is_empty()
            || request.vector_catalog_id.is_empty()
        {
            return Err(Self::regulatory_fragment_error(
                "missing_required_identity",
                [&request.plan_id],
                "plan_id, vector_seq_id, and vector_catalog_id must be non-empty.",
            ));
        }
        if request.insertion_context_feature_id != "multiple_cloning_region" {
            return Err(Self::regulatory_fragment_error(
                "unsupported_insertion_context",
                [&request.insertion_context_feature_id],
                "V1 supports only the catalog-validated multiple_cloning_region insertion context.",
            ));
        }
        if let Some(path) = request.helper_catalog_path.as_deref() {
            request.helper_catalog_path = Some(
                fs::canonicalize(path.trim())
                    .map_err(|error| {
                        Self::regulatory_fragment_error(
                            "helper_catalog_unresolvable",
                            [path],
                            format!("Could not resolve helper-vector catalog: {error}"),
                        )
                    })?
                    .to_string_lossy()
                    .to_string(),
            );
        }
        if request.max_panel_members == 0 || request.max_construct_length_bp == 0 {
            return Err(Self::regulatory_fragment_error(
                "invalid_panel_bound",
                [&request.plan_id],
                "max_panel_members and max_construct_length_bp must both be at least 1.",
            ));
        }
        if request.policy.max_candidate_constructs == 0
            || request.policy.max_candidate_constructs > 20
        {
            return Err(Self::regulatory_fragment_error(
                "invalid_candidate_search_bound",
                [&request.plan_id],
                "policy.max_candidate_constructs must be within 1..=20 for exhaustive deterministic selection.",
            ));
        }
        if !(1..=64).contains(&request.policy.sequence_word_size_bp)
            || request.policy.near_exact_max_mismatches > 8
            || !(1..=512).contains(&request.policy.junction_flank_bp)
            || !(1..=10_000).contains(&request.policy.max_evidence_observations_per_dimension)
        {
            return Err(Self::regulatory_fragment_error(
                "invalid_evidence_policy",
                [&request.plan_id],
                "Evidence policy requires sequence_word_size_bp within 1..=64, near_exact_max_mismatches <= 8, junction_flank_bp within 1..=512, and max_evidence_observations_per_dimension within 1..=10000.",
            ));
        }
        for (name, value) in [
            (
                "high_similarity_identity_fraction",
                request.policy.high_similarity_identity_fraction,
            ),
            (
                "high_similarity_coverage_fraction",
                request.policy.high_similarity_coverage_fraction,
            ),
        ] {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(Self::regulatory_fragment_error(
                    "invalid_similarity_threshold",
                    [name],
                    format!("{name} must be finite and within [0, 1]."),
                ));
            }
        }
        if request.fragments.is_empty() || request.questions.is_empty() {
            return Err(Self::regulatory_fragment_error(
                "empty_panel_request",
                [&request.plan_id],
                "At least one fragment and one requested experimental question are required.",
            ));
        }

        let mut fragment_ids = HashSet::new();
        for fragment in &mut request.fragments {
            fragment.fragment_id = fragment.fragment_id.trim().to_string();
            fragment.region_set_id = fragment.region_set_id.trim().to_string();
            fragment.region_set_content_sha256 =
                fragment.region_set_content_sha256.trim().to_string();
            fragment.reference_release = fragment.reference_release.trim().to_string();
            if fragment.fragment_id.is_empty()
                || fragment.region_set_id.is_empty()
                || fragment.reference_release.is_empty()
            {
                return Err(Self::regulatory_fragment_error(
                    "incomplete_fragment_binding",
                    [&fragment.fragment_id],
                    "Each fragment requires fragment_id, region_set_id, and reference_release.",
                ));
            }
            if !fragment_ids.insert(fragment.fragment_id.clone()) {
                return Err(Self::regulatory_fragment_error(
                    "duplicate_fragment_id",
                    [&fragment.fragment_id],
                    "Fragment ids must be unique.",
                ));
            }
        }
        request.fragments.sort_by(|left, right| {
            left.declared_order
                .cmp(&right.declared_order)
                .then(left.role.cmp(&right.role))
                .then(left.fragment_id.cmp(&right.fragment_id))
        });
        let count_role = |role| {
            request
                .fragments
                .iter()
                .filter(|fragment| fragment.role == role)
                .count()
        };
        if count_role(RegulatoryFragmentRole::Candidate) != 1
            || count_role(RegulatoryFragmentRole::Partner) > 1
            || count_role(RegulatoryFragmentRole::MinimalPromoter) > 1
            || count_role(RegulatoryFragmentRole::ReferenceControl) > 1
        {
            return Err(Self::regulatory_fragment_error(
                "unsupported_role_cardinality",
                [&request.plan_id],
                "V1 requires exactly one candidate and at most one partner, minimal_promoter, and reference_control.",
            ));
        }
        if request.policy.include_minimal_promoter_control
            && count_role(RegulatoryFragmentRole::MinimalPromoter) != 1
        {
            return Err(Self::regulatory_fragment_error(
                "minimal_promoter_binding_required",
                [&request.plan_id],
                "The effective policy requires exactly one persisted minimal_promoter ROI.",
            ));
        }

        request.questions.sort();
        request.questions.dedup();
        let has_partner = count_role(RegulatoryFragmentRole::Partner) == 1;
        if !has_partner
            && request.questions.iter().any(|question| {
                matches!(
                    question,
                    RegulatoryFragmentQuestion::PartnerDependence
                        | RegulatoryFragmentQuestion::OrderDependence
                        | RegulatoryFragmentQuestion::OrientationDependence
                        | RegulatoryFragmentQuestion::SpacingDependence
                )
            })
        {
            return Err(Self::regulatory_fragment_error(
                "partner_binding_required",
                [&request.plan_id],
                "Partner/order/orientation/spacing questions require an explicitly bound partner ROI.",
            ));
        }

        let known_ids = request
            .fragments
            .iter()
            .map(|fragment| fragment.fragment_id.clone())
            .collect::<HashSet<_>>();
        let mut geometry_ids = HashSet::new();
        if let Some(reference) = request.reference_combination.as_mut() {
            if reference.kind != RegulatoryFragmentGeometryKind::ReferenceCombination {
                return Err(Self::regulatory_fragment_error(
                    "invalid_reference_geometry_kind",
                    [&reference.variant_id],
                    "reference_combination must use kind=reference_combination.",
                ));
            }
            Self::normalize_regulatory_fragment_geometry(reference, &known_ids)?;
            geometry_ids.insert(reference.variant_id.clone());
        }
        for variant in &mut request.requested_variants {
            if variant.kind == RegulatoryFragmentGeometryKind::ReferenceCombination {
                return Err(Self::regulatory_fragment_error(
                    "invalid_requested_variant_kind",
                    [&variant.variant_id],
                    "requested_variants may not use kind=reference_combination.",
                ));
            }
            Self::normalize_regulatory_fragment_geometry(variant, &known_ids)?;
            if !geometry_ids.insert(variant.variant_id.clone()) {
                return Err(Self::regulatory_fragment_error(
                    "duplicate_geometry_id",
                    [&variant.variant_id],
                    "Geometry ids must be unique across reference and requested variants.",
                ));
            }
            if matches!(
                variant.kind,
                RegulatoryFragmentGeometryKind::BoundaryShift
                    | RegulatoryFragmentGeometryKind::Tiling
            ) && variant.extended_boundary.is_none()
            {
                return Err(Self::regulatory_fragment_error(
                    "boundary_variant_missing_basis",
                    [&variant.variant_id],
                    "Boundary-shift and tiling variants must cite the existing extended_boundary transcript policy.",
                ));
            }
        }
        request.requested_variants.sort_by(|left, right| {
            left.declared_order
                .cmp(&right.declared_order)
                .then(left.kind.cmp(&right.kind))
                .then(left.variant_id.cmp(&right.variant_id))
        });
        let needs_reference = request.questions.iter().any(|question| {
            matches!(
                question,
                RegulatoryFragmentQuestion::PartnerDependence
                    | RegulatoryFragmentQuestion::OrderDependence
                    | RegulatoryFragmentQuestion::OrientationDependence
                    | RegulatoryFragmentQuestion::SpacingDependence
            )
        });
        if needs_reference && request.reference_combination.is_none() {
            return Err(Self::regulatory_fragment_error(
                "reference_combination_required",
                [&request.plan_id],
                "The requested partner/geometry comparisons require one exact reference_combination.",
            ));
        }
        for (question, kind) in [
            (
                RegulatoryFragmentQuestion::OrderDependence,
                RegulatoryFragmentGeometryKind::ReversedOrder,
            ),
            (
                RegulatoryFragmentQuestion::OrientationDependence,
                RegulatoryFragmentGeometryKind::ReversedOrientation,
            ),
            (
                RegulatoryFragmentQuestion::SpacingDependence,
                RegulatoryFragmentGeometryKind::ControlledSpacing,
            ),
            (
                RegulatoryFragmentQuestion::BoundaryUncertainty,
                RegulatoryFragmentGeometryKind::BoundaryShift,
            ),
        ] {
            if request.questions.contains(&question)
                && !request.requested_variants.iter().any(|row| {
                    row.kind == kind
                        || (question == RegulatoryFragmentQuestion::BoundaryUncertainty
                            && row.kind == RegulatoryFragmentGeometryKind::Tiling)
                })
            {
                return Err(Self::regulatory_fragment_error(
                    "requested_question_has_no_explicit_variant",
                    [Self::regulatory_fragment_question_token(question)],
                    "Geometry and boundary questions require at least one explicit matching variant; the planner will not invent one.",
                ));
            }
        }

        if request
            .questions
            .contains(&RegulatoryFragmentQuestion::MotifDisruption)
        {
            if request.mutation_policy == PromoterReporterPanelMutationPolicy::Unspecified {
                return Err(Self::regulatory_fragment_error(
                    "mutation_policy_unspecified",
                    [&request.plan_id],
                    "Motif disruption requires the existing non-Unspecified promoter-panel mutation policy.",
                ));
            }
            if request.mutation_policy
                != PromoterReporterPanelMutationPolicy::P53FamilyCoreDisruptionV1
            {
                return Err(Self::regulatory_fragment_error(
                    "mutation_policy_does_not_define_edit",
                    [&request.plan_id],
                    "native_only_v1 does not define a motif-disruption edit.",
                ));
            }
            let target = request.motif_disruption.as_mut().ok_or_else(|| {
                Self::regulatory_fragment_error(
                    "motif_disruption_target_required",
                    [&request.plan_id],
                    "Motif disruption requires an exact fragment-local target interval.",
                )
            })?;
            target.control_id = Self::regulatory_fragment_slug(&target.control_id);
            target.fragment_id = target.fragment_id.trim().to_string();
            if target.control_id.is_empty()
                || !known_ids.contains(&target.fragment_id)
                || target.motif_start_in_fragment_0based
                    >= target.motif_end_in_fragment_0based_exclusive
            {
                return Err(Self::regulatory_fragment_error(
                    "invalid_motif_disruption_target",
                    [&target.fragment_id],
                    "The motif control needs a unique id, a bound fragment, and a non-empty 0-based half-open interval.",
                ));
            }
        } else if request.motif_disruption.is_some() {
            return Err(Self::regulatory_fragment_error(
                "unrequested_motif_disruption",
                [&request.plan_id],
                "motif_disruption input is valid only when motif_disruption is a requested question.",
            ));
        }

        for binding in &mut request.evidence_bindings {
            binding.report_id = binding.report_id.trim().to_string();
            binding.report_sha256 = binding.report_sha256.trim().to_ascii_lowercase();
            if binding
                .report_path
                .as_ref()
                .is_some_and(|path| path.trim().is_empty())
            {
                return Err(Self::regulatory_fragment_error(
                    "empty_evidence_path",
                    [&binding.report_id],
                    "Use a non-empty path or omit report_path for a citation-only binding.",
                ));
            }
            binding.row_id = binding
                .row_id
                .take()
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty());
            if binding.report_id.is_empty()
                || !Self::regulatory_fragment_valid_sha256(&binding.report_sha256)
            {
                return Err(Self::regulatory_fragment_error(
                    "invalid_evidence_binding",
                    [&binding.report_id],
                    "Evidence bindings require a report_id and exact sha256: digest.",
                ));
            }
        }
        request.evidence_bindings.sort_by(|left, right| {
            left.dimension
                .cmp(&right.dimension)
                .then(left.report_id.cmp(&right.report_id))
                .then(left.row_id.cmp(&right.row_id))
                .then(left.report_path.cmp(&right.report_path))
                .then(left.report_sha256.cmp(&right.report_sha256))
        });
        request.evidence_bindings.dedup();
        request.scientific_caveats = request
            .scientific_caveats
            .into_iter()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        Ok(request)
    }

    fn normalize_regulatory_fragment_geometry(
        geometry: &mut RegulatoryFragmentGeometryRequest,
        known_ids: &HashSet<String>,
    ) -> Result<(), EngineError> {
        geometry.variant_id = Self::regulatory_fragment_slug(&geometry.variant_id);
        if geometry.variant_id.is_empty() || geometry.instances.is_empty() {
            return Err(Self::regulatory_fragment_error(
                "empty_geometry",
                [&geometry.variant_id],
                "Every geometry requires a stable id and at least one ordered instance.",
            ));
        }
        let mut seen = HashSet::new();
        for (index, instance) in geometry.instances.iter_mut().enumerate() {
            instance.fragment_id = instance.fragment_id.trim().to_string();
            instance.spacer_before = Self::normalize_regulatory_fragment_dna(
                &instance.spacer_before,
                "geometry spacer",
            )?;
            if !known_ids.contains(&instance.fragment_id) {
                return Err(Self::regulatory_fragment_error(
                    "geometry_references_unknown_fragment",
                    [&geometry.variant_id, &instance.fragment_id],
                    "Every geometry instance must reference an exact bound fragment.",
                ));
            }
            if !seen.insert(instance.fragment_id.clone()) {
                return Err(Self::regulatory_fragment_error(
                    "geometry_repeats_fragment",
                    [&geometry.variant_id, &instance.fragment_id],
                    "V1 does not infer the meaning of repeated instances of one ROI.",
                ));
            }
            if index == 0 && !instance.spacer_before.is_empty() {
                return Err(Self::regulatory_fragment_error(
                    "leading_spacer_not_supported",
                    [&geometry.variant_id],
                    "The first geometry instance must have an empty spacer_before.",
                ));
            }
        }
        Ok(())
    }

    fn resolve_regulatory_fragment_bindings(
        &self,
        request: &RegulatoryFragmentPanelRequest,
    ) -> Result<
        (
            Vec<ResolvedRegulatoryFragment>,
            Vec<RegulatoryFragmentFinding>,
        ),
        EngineError,
    > {
        let store = self.genomic_region_store_snapshot()?;
        let mut resolved = vec![];
        let mut warnings = vec![];
        let mut assemblies = BTreeSet::new();
        let mut releases = BTreeSet::new();
        for binding in &request.fragments {
            let set = store
                .sets
                .iter()
                .find(|set| set.set_id == binding.region_set_id)
                .ok_or_else(|| {
                    Self::regulatory_fragment_error(
                        "region_set_not_found",
                        [&binding.fragment_id, &binding.region_set_id],
                        "The bound persisted genomic region set is unavailable.",
                    )
                })?;
            if set.content_sha256 != binding.region_set_content_sha256 {
                return Err(Self::regulatory_fragment_error(
                    "region_set_content_digest_mismatch",
                    [&binding.fragment_id, &binding.region_set_id],
                    format!(
                        "Bound digest '{}' differs from current digest '{}'.",
                        binding.region_set_content_sha256, set.content_sha256
                    ),
                ));
            }
            let persisted = set
                .regions
                .iter()
                .find(|region| region.region_id == binding.region.region_id)
                .ok_or_else(|| {
                    Self::regulatory_fragment_error(
                        "region_not_found",
                        [&binding.fragment_id, &binding.region.region_id],
                        "The bound persisted genomic ROI is unavailable.",
                    )
                })?;
            if persisted.identity_sha256 != binding.region.identity_sha256
                || persisted.content_sha256 != binding.region.content_sha256
                || persisted != &binding.region
            {
                return Err(Self::regulatory_fragment_error(
                    "region_content_digest_mismatch",
                    [&binding.fragment_id, &binding.region.region_id],
                    "The exact bound ROI snapshot differs from the current persisted ROI content.",
                ));
            }
            let projection = persisted.local_projection.as_ref().ok_or_else(|| {
                Self::regulatory_fragment_error(
                    "region_projection_absent",
                    [&binding.fragment_id, &binding.region.region_id],
                    "A bound regulatory fragment requires a local_projection.",
                )
            })?;
            if projection.status != gp::GenomicRegionLocalProjectionStatus::Current {
                return Err(Self::regulatory_fragment_error(
                    "region_projection_not_current",
                    [&binding.fragment_id, &binding.region.region_id],
                    format!(
                        "Stored local_projection status is '{}'.",
                        Self::regulatory_fragment_projection_status_token(projection.status)
                    ),
                ));
            }
            let mut effective = persisted.clone();
            self.refresh_projection_status(&mut effective)?;
            let effective_projection = effective.local_projection.as_ref().ok_or_else(|| {
                Self::regulatory_fragment_error(
                    "region_projection_absent",
                    [&binding.fragment_id, &binding.region.region_id],
                    "A bound regulatory fragment requires a local_projection.",
                )
            })?;
            if effective_projection.status != gp::GenomicRegionLocalProjectionStatus::Current {
                return Err(Self::regulatory_fragment_error(
                    "region_projection_not_current",
                    [&binding.fragment_id, &binding.region.region_id],
                    format!(
                        "Current local_projection status is '{}'.",
                        Self::regulatory_fragment_projection_status_token(
                            effective_projection.status
                        )
                    ),
                ));
            }
            let source = self
                .snapshot()
                .sequences
                .get(&effective_projection.seq_id)
                .ok_or_else(|| {
                    Self::regulatory_fragment_error(
                        "source_sequence_missing",
                        [&binding.fragment_id, &effective_projection.seq_id],
                        "The exact source sequence is not loaded.",
                    )
                })?;
            let observed_source_digest = sha256_prefixed_bytes(source.forward_bytes());
            if observed_source_digest != effective_projection.sequence_sha256 {
                return Err(Self::regulatory_fragment_error(
                    "source_sequence_digest_mismatch",
                    [&binding.fragment_id, &effective_projection.seq_id],
                    "The loaded source sequence differs from the projection-bound digest.",
                ));
            }
            let start = usize::try_from(effective_projection.local_start_0based).map_err(|_| {
                Self::regulatory_fragment_error(
                    "source_interval_overflow",
                    [&binding.fragment_id],
                    "The source interval does not fit this platform.",
                )
            })?;
            let end =
                usize::try_from(effective_projection.local_end_0based_exclusive).map_err(|_| {
                    Self::regulatory_fragment_error(
                        "source_interval_overflow",
                        [&binding.fragment_id],
                        "The source interval does not fit this platform.",
                    )
                })?;
            let source_sequence = source.get_forward_string();
            let slice = source_sequence.get(start..end).ok_or_else(|| {
                Self::regulatory_fragment_error(
                    "region_projection_not_current",
                    [&binding.fragment_id, &binding.region.region_id],
                    "Current projection resolves outside the loaded source sequence.",
                )
            })?;
            let canonical_sequence =
                if effective_projection.local_strand == gp::GenomicRegionStrand::Minus {
                    Self::reverse_complement(slice)
                } else {
                    slice.to_string()
                };
            let fragment_warnings = Self::regulatory_fragment_role_purpose_warnings(binding);
            warnings.extend(fragment_warnings.iter().cloned());
            assemblies.insert((
                binding.region.interval.reference.assembly_name.clone(),
                binding.region.interval.reference.assembly_accession.clone(),
            ));
            releases.insert(binding.reference_release.clone());
            resolved.push(ResolvedRegulatoryFragment {
                binding: binding.clone(),
                projection: effective_projection.clone(),
                canonical_sequence,
                warnings: fragment_warnings,
            });
        }
        if assemblies.len() > 1 {
            return Err(Self::regulatory_fragment_error(
                "mixed_reference_assemblies",
                request
                    .fragments
                    .iter()
                    .map(|fragment| fragment.fragment_id.as_str()),
                "All bound regulatory fragments must use one exact assembly; no liftover is performed.",
            ));
        }
        if releases.len() > 1 {
            return Err(Self::regulatory_fragment_error(
                "mixed_reference_releases",
                request
                    .fragments
                    .iter()
                    .map(|fragment| fragment.fragment_id.as_str()),
                "All bound regulatory fragments must use one exact annotation/reference release.",
            ));
        }
        Ok((resolved, warnings))
    }

    fn resolve_regulatory_fragment_vector_context(
        &self,
        request: &RegulatoryFragmentPanelRequest,
    ) -> Result<
        (
            ReporterVectorValidationReport,
            RegulatoryFragmentVectorContext,
        ),
        EngineError,
    > {
        let backbone = self.resolve_reporter_backbone(
            &request.vector_seq_id,
            None,
            true,
            Some(&request.vector_catalog_id),
            request.helper_catalog_path.as_deref(),
        )?;
        let validation = backbone.validation.ok_or_else(|| {
            Self::regulatory_fragment_error(
                "vector_identity_unavailable",
                [&request.vector_seq_id, &request.vector_catalog_id],
                "Catalog-owned exact vector validation did not produce a report.",
            )
        })?;
        if validation.status != ReporterVectorValidationStatus::Verified {
            return Err(Self::regulatory_fragment_error(
                "vector_identity_not_verified",
                [&request.vector_seq_id, &request.vector_catalog_id],
                format!(
                    "Exact vector validation status is '{:?}'.",
                    validation.status
                ),
            ));
        }
        let vector = self
            .snapshot()
            .sequences
            .get(&request.vector_seq_id)
            .ok_or_else(|| {
                Self::regulatory_fragment_error(
                    "vector_sequence_missing",
                    [&request.vector_seq_id],
                    "The exact catalog-validated reporter vector is not loaded.",
                )
            })?;
        let observed = validation
            .observed_multiple_cloning_region
            .as_deref()
            .ok_or_else(|| {
                Self::regulatory_fragment_error(
                    "vector_insertion_context_unresolved",
                    [
                        &request.vector_seq_id,
                        &request.insertion_context_feature_id,
                    ],
                    "Exact vector validation did not resolve the requested insertion feature.",
                )
            })?;
        let interval = observed.split_whitespace().next().unwrap_or(observed);
        let (start, end) = interval.split_once("..").ok_or_else(|| {
            Self::regulatory_fragment_error(
                "vector_insertion_context_unresolved",
                [&request.vector_seq_id],
                format!("Could not parse validated insertion interval '{observed}'."),
            )
        })?;
        let start_1based = start.parse::<usize>().map_err(|_| {
            Self::regulatory_fragment_error(
                "vector_insertion_context_unresolved",
                [&request.vector_seq_id],
                format!("Could not parse validated insertion interval '{observed}'."),
            )
        })?;
        let end_1based = end.parse::<usize>().map_err(|_| {
            Self::regulatory_fragment_error(
                "vector_insertion_context_unresolved",
                [&request.vector_seq_id],
                format!("Could not parse validated insertion interval '{observed}'."),
            )
        })?;
        if start_1based == 0 || start_1based > end_1based || end_1based > vector.len() {
            return Err(Self::regulatory_fragment_error(
                "vector_insertion_context_unresolved",
                [&request.vector_seq_id],
                format!("Validated insertion interval '{observed}' is outside the vector."),
            ));
        }
        let start_0based = start_1based - 1;
        let end_0based_exclusive = end_1based;
        let vector_sequence = vector.get_forward_string();
        let flank = request.policy.vector_context_flank_bp;
        let left = vector_sequence[start_0based.saturating_sub(flank)..start_0based].to_string();
        let right = vector_sequence[end_0based_exclusive
            ..end_0based_exclusive
                .saturating_add(flank)
                .min(vector_sequence.len())]
            .to_string();
        let vector_sequence_sha256 = sha256_prefixed_str(&vector_sequence);
        let context_sha256 = Self::regulatory_fragment_value_sha256(
            &serde_json::json!({
                "vector_seq_id": request.vector_seq_id,
                "vector_catalog_id": request.vector_catalog_id,
                "vector_sequence_sha256": vector_sequence_sha256,
                "insertion_context_feature_id": request.insertion_context_feature_id,
                "start_0based": start_0based,
                "end_0based_exclusive": end_0based_exclusive,
                "left_flank": left,
                "right_flank": right,
            }),
            "regulatory-fragment vector context",
        )?;
        Ok((
            validation,
            RegulatoryFragmentVectorContext {
                vector_seq_id: request.vector_seq_id.clone(),
                vector_catalog_id: request.vector_catalog_id.clone(),
                vector_sequence_sha256,
                insertion_context_feature_id: request.insertion_context_feature_id.clone(),
                insertion_start_0based: start_0based,
                insertion_end_0based_exclusive: end_0based_exclusive,
                left_flank_5prime_to_3prime: left,
                right_flank_5prime_to_3prime: right,
                context_sha256,
            },
        ))
    }

    fn regulatory_fragment_candidate_partner_similarity(
        &self,
        request: &RegulatoryFragmentPanelRequest,
        fragments: &HashMap<&str, &ResolvedRegulatoryFragment>,
    ) -> Result<Option<RegulatoryFragmentSimilarityAudit>, EngineError> {
        let candidate = request
            .fragments
            .iter()
            .find(|fragment| fragment.role == RegulatoryFragmentRole::Candidate)
            .and_then(|fragment| fragments.get(fragment.fragment_id.as_str()).copied());
        let partner = request
            .fragments
            .iter()
            .find(|fragment| fragment.role == RegulatoryFragmentRole::Partner)
            .and_then(|fragment| fragments.get(fragment.fragment_id.as_str()).copied());
        let (Some(candidate), Some(partner)) = (candidate, partner) else {
            return Ok(None);
        };
        let computed = Self::compute_pairwise_alignment_report(
            &candidate.binding.fragment_id,
            &candidate.canonical_sequence,
            None,
            None,
            &partner.binding.fragment_id,
            &partner.canonical_sequence,
            None,
            None,
            PairwiseAlignmentMode::Global,
            2,
            -3,
            -5,
            -1,
        )?;
        let highly_similar = computed.report.identity_fraction
            >= request.policy.high_similarity_identity_fraction
            && computed.report.query_coverage_fraction
                >= request.policy.high_similarity_coverage_fraction
            && computed.report.target_coverage_fraction
                >= request.policy.high_similarity_coverage_fraction;
        Ok(Some(RegulatoryFragmentSimilarityAudit {
            left_fragment_id: candidate.binding.fragment_id.clone(),
            right_fragment_id: partner.binding.fragment_id.clone(),
            highly_similar,
            alignment: computed.report,
        }))
    }

    fn regulatory_fragment_construct_candidates(
        &self,
        request: &RegulatoryFragmentPanelRequest,
        fragments: &HashMap<&str, &ResolvedRegulatoryFragment>,
        vector_context: &RegulatoryFragmentVectorContext,
        highly_similar: bool,
    ) -> Result<
        (
            Vec<CandidateRegulatoryConstruct>,
            Vec<RegulatoryFragmentOmittedVariant>,
        ),
        EngineError,
    > {
        let fragment_for_role = |role| {
            request
                .fragments
                .iter()
                .find(|fragment| fragment.role == role)
        };
        let candidate = fragment_for_role(RegulatoryFragmentRole::Candidate)
            .expect("normalization requires candidate");
        let partner = fragment_for_role(RegulatoryFragmentRole::Partner);
        let minimal = fragment_for_role(RegulatoryFragmentRole::MinimalPromoter);
        let reference_control = fragment_for_role(RegulatoryFragmentRole::ReferenceControl);
        let instance = |fragment: &RegulatoryFragmentBinding| RegulatoryFragmentInstanceRequest {
            fragment_id: fragment.fragment_id.clone(),
            orientation: RegulatoryFragmentOrientation::Forward,
            spacer_before: String::new(),
        };
        let with_minimal = |fragment: &RegulatoryFragmentBinding| {
            let mut instances = vec![instance(fragment)];
            if let Some(minimal) = minimal {
                instances.push(instance(minimal));
            }
            instances
        };
        let mut specs = vec![];
        if request.policy.include_promoterless_control {
            specs.push((
                "control_promoterless".to_string(),
                RegulatoryFragmentConstructKind::PromoterlessControl,
                vec![],
                true,
                None,
                0,
            ));
        }
        if request.policy.include_minimal_promoter_control {
            let minimal = minimal.expect("normalization requires minimal promoter");
            specs.push((
                "control_minimal_promoter".to_string(),
                RegulatoryFragmentConstructKind::MinimalPromoterControl,
                vec![instance(minimal)],
                true,
                None,
                minimal.declared_order,
            ));
        }
        specs.push((
            "candidate_alone".to_string(),
            RegulatoryFragmentConstructKind::CandidateAlone,
            with_minimal(candidate),
            false,
            None,
            candidate.declared_order,
        ));
        if let Some(partner) = partner {
            specs.push((
                "partner_alone".to_string(),
                RegulatoryFragmentConstructKind::PartnerAlone,
                with_minimal(partner),
                false,
                None,
                partner.declared_order,
            ));
        }
        if let Some(reference) = request.reference_combination.as_ref() {
            specs.push((
                "reference_combination".to_string(),
                RegulatoryFragmentConstructKind::ReferenceCombination,
                reference.instances.clone(),
                false,
                None,
                reference.declared_order,
            ));
        }
        if request.policy.include_reference_control_when_bound
            && let Some(reference_control) = reference_control
        {
            specs.push((
                "control_reference".to_string(),
                RegulatoryFragmentConstructKind::ReferenceControl,
                with_minimal(reference_control),
                true,
                None,
                reference_control.declared_order,
            ));
        }
        for variant in &request.requested_variants {
            specs.push((
                format!("variant_{}", variant.variant_id),
                RegulatoryFragmentConstructKind::RequestedGeometryVariant,
                variant.instances.clone(),
                false,
                None,
                variant.declared_order,
            ));
        }
        if let Some(motif) = request.motif_disruption.as_ref() {
            let base_instances = request
                .reference_combination
                .as_ref()
                .filter(|geometry| {
                    geometry
                        .instances
                        .iter()
                        .any(|instance| instance.fragment_id == motif.fragment_id)
                })
                .map(|geometry| geometry.instances.clone())
                .unwrap_or_else(|| {
                    let fragment = request
                        .fragments
                        .iter()
                        .find(|fragment| fragment.fragment_id == motif.fragment_id)
                        .expect("normalized motif fragment");
                    with_minimal(fragment)
                });
            specs.push((
                format!("motif_{}", motif.control_id),
                RegulatoryFragmentConstructKind::MotifDisruptionControl,
                base_instances,
                false,
                Some(motif.clone()),
                request
                    .fragments
                    .iter()
                    .find(|fragment| fragment.fragment_id == motif.fragment_id)
                    .map(|fragment| fragment.declared_order)
                    .unwrap_or_default(),
            ));
        }
        if specs.len() > request.policy.max_candidate_constructs {
            return Err(Self::regulatory_fragment_error(
                "candidate_construct_bound_exceeded",
                [&request.plan_id],
                format!(
                    "The explicit request yields {} construct candidates, exceeding policy.max_candidate_constructs={}; reduce explicit variants.",
                    specs.len(),
                    request.policy.max_candidate_constructs
                ),
            ));
        }

        let mut out = vec![];
        let mut omitted = vec![];
        let mut seen_geometry = BTreeMap::<String, String>::new();
        for (member_token, kind, instances, mandatory, motif, declared_order) in specs {
            let member_id = format!("{}_{}", request.plan_id, member_token);
            let mut member = self.resolve_regulatory_fragment_construct(
                request,
                fragments,
                vector_context,
                &member_id,
                kind,
                &instances,
                motif.as_ref(),
                highly_similar,
            )?;
            let geometry_key = Self::regulatory_fragment_value_sha256(
                &serde_json::json!({
                    "instances": member.instances,
                    "insert_sequence_sha256": member.insert_sequence_sha256,
                }),
                "regulatory-fragment construct geometry",
            )?;
            let duplicate_of = seen_geometry.get(&geometry_key).cloned();
            if duplicate_of.is_none() {
                seen_geometry.insert(geometry_key, member_id.clone());
            }
            let questions = Self::regulatory_fragment_questions_for_construct(
                request,
                kind,
                member_token.strip_prefix("variant_").unwrap_or_default(),
            );
            let eligible = member.insert_length_bp <= request.max_construct_length_bp
                && duplicate_of.is_none();
            if let Some(duplicate_of) = duplicate_of {
                omitted.push(RegulatoryFragmentOmittedVariant {
                    member_id: member_id.clone(),
                    construct_kind: kind,
                    reason: RegulatoryFragmentOmissionReasonKind::DuplicateGeometry,
                    affected_questions: questions.iter().copied().collect(),
                    detail: format!(
                        "This exact ordered insert duplicates '{}'; the first declared geometry is retained.",
                        duplicate_of
                    ),
                });
                continue;
            }
            member.planning_label = if highly_similar {
                RegulatoryFragmentPlanningLabel::ContextOrGeometryConfounded
            } else {
                RegulatoryFragmentPlanningLabel::StandaloneTestableCandidate
            };
            out.push(CandidateRegulatoryConstruct {
                member,
                sort_key: (
                    declared_order,
                    Self::regulatory_fragment_construct_role_rank(kind),
                    member_id,
                ),
                mandatory,
                eligible,
                questions,
            });
        }
        out.sort_by(|left, right| left.sort_key.cmp(&right.sort_key));
        Ok((out, omitted))
    }

    #[allow(clippy::too_many_arguments)]
    fn resolve_regulatory_fragment_construct(
        &self,
        request: &RegulatoryFragmentPanelRequest,
        fragments: &HashMap<&str, &ResolvedRegulatoryFragment>,
        vector_context: &RegulatoryFragmentVectorContext,
        member_id: &str,
        kind: RegulatoryFragmentConstructKind,
        instances: &[RegulatoryFragmentInstanceRequest],
        motif: Option<&RegulatoryFragmentMotifDisruptionRequest>,
        highly_similar: bool,
    ) -> Result<RegulatoryFragmentPanelMember, EngineError> {
        let mut assembled = String::new();
        let mut resolved_instances = vec![];
        let mut warnings = vec![];
        for instance in instances {
            let fragment = fragments
                .get(instance.fragment_id.as_str())
                .copied()
                .ok_or_else(|| {
                    Self::regulatory_fragment_error(
                        "construct_references_unknown_fragment",
                        [member_id, &instance.fragment_id],
                        "Construct resolution requires an exact bound fragment.",
                    )
                })?;
            assembled.push_str(&instance.spacer_before);
            let assembled_start = assembled.len();
            let mut sequence = fragment.canonical_sequence.clone();
            if motif.is_some_and(|target| target.fragment_id == instance.fragment_id) {
                let target = motif.expect("checked motif target");
                if target.motif_end_in_fragment_0based_exclusive > sequence.len() {
                    return Err(Self::regulatory_fragment_error(
                        "motif_disruption_target_outside_fragment",
                        [member_id, &target.fragment_id],
                        "The exact motif-disruption interval exceeds the resolved fragment sequence.",
                    ));
                }
                sequence = self
                    .design_promoter_reporter_panel_mutation(
                        request.mutation_policy,
                        &sequence,
                        target.motif_start_in_fragment_0based,
                        target.motif_end_in_fragment_0based_exclusive,
                        target.motif_forward_strand,
                        &[],
                    )?
                    .mutant_sequence;
            }
            if instance.orientation == RegulatoryFragmentOrientation::ReverseComplement {
                sequence = Self::reverse_complement(&sequence);
            }
            assembled.push_str(&sequence);
            let assembled_end = assembled.len();
            warnings.extend(fragment.warnings.iter().cloned());
            resolved_instances.push(RegulatoryFragmentResolvedInstance {
                fragment_id: fragment.binding.fragment_id.clone(),
                role: fragment.binding.role,
                region_set_id: fragment.binding.region_set_id.clone(),
                region_id: fragment.binding.region.region_id.clone(),
                region_identity_sha256: fragment.binding.region.identity_sha256.clone(),
                region_content_sha256: fragment.binding.region.content_sha256.clone(),
                assembly: fragment
                    .binding
                    .region
                    .interval
                    .reference
                    .assembly_name
                    .clone(),
                reference_release: fragment.binding.reference_release.clone(),
                contig: fragment
                    .binding
                    .region
                    .interval
                    .reference
                    .contig_name
                    .clone(),
                genomic_start_0based: fragment.binding.region.interval.start_0based,
                genomic_end_0based_exclusive: fragment.binding.region.interval.end_0based_exclusive,
                genomic_strand: fragment.binding.region.interval.strand,
                source_seq_id: fragment.projection.seq_id.clone(),
                source_sequence_sha256: fragment.projection.sequence_sha256.clone(),
                source_start_0based: fragment.projection.local_start_0based,
                source_end_0based_exclusive: fragment.projection.local_end_0based_exclusive,
                source_strand: fragment.projection.local_strand,
                orientation: instance.orientation,
                spacer_before: instance.spacer_before.clone(),
                assembled_start_0based: assembled_start,
                assembled_end_0based_exclusive: assembled_end,
            });
        }
        let output_id = format!("regulatory_insert_{member_id}");
        let planned_operations = if assembled.is_empty() {
            vec![]
        } else {
            vec![
                serde_json::to_value(Operation::CreateSequenceFromText {
                    sequence_text: assembled.clone(),
                    output_id: Some(output_id),
                    name: Some(member_id.to_string()),
                    circular: false,
                })
                .map_err(|error| {
                    Self::regulatory_fragment_error(
                        "planned_operation_serialization_failed",
                        [member_id],
                        format!("Could not serialize exact insert operation: {error}"),
                    )
                })?,
            ]
        };
        let mut nonclaims = vec![
            "This construct is a planned experimental comparison, not evidence of regulatory sufficiency or partner dependence."
                .to_string(),
        ];
        if highly_similar {
            nonclaims.push(
                "Candidate/partner sequence similarity can confound attribution and is not evidence of shared function."
                    .to_string(),
            );
        }
        Ok(RegulatoryFragmentPanelMember {
            member_id: member_id.to_string(),
            construct_kind: kind,
            instances: resolved_instances,
            insert_sequence_sha256: sha256_prefixed_str(&assembled),
            insert_length_bp: assembled.len(),
            insert_sequence_5prime_to_3prime: assembled,
            vector_context: vector_context.clone(),
            cloning_feasibility: if instances.is_empty() {
                RegulatoryFragmentCloningFeasibilityState::NotApplicableEmptyInsert
            } else {
                RegulatoryFragmentCloningFeasibilityState::NotEvaluated
            },
            final_product_audit_state: RegulatoryFragmentFinalProductAuditState::NotEvaluated,
            planned_operations,
            warnings,
            nonclaims,
            ..RegulatoryFragmentPanelMember::default()
        })
    }

    fn regulatory_fragment_contrast_candidates(
        request: &RegulatoryFragmentPanelRequest,
        candidates: &[CandidateRegulatoryConstruct],
    ) -> Vec<CandidateContrast> {
        let id = |suffix: &str| format!("{}_{}", request.plan_id, suffix);
        let existing = candidates
            .iter()
            .filter(|candidate| candidate.eligible)
            .map(|candidate| candidate.member.member_id.clone())
            .collect::<BTreeSet<_>>();
        let mut contrasts = vec![];
        let mut push = |question, left: String, right: String, interpretation: &str| {
            if existing.contains(&left) && existing.contains(&right) {
                contrasts.push(CandidateContrast {
                    contrast: RegulatoryFragmentContrast {
                        contrast_id: format!(
                            "{}_{}_{}_vs_{}",
                            request.plan_id,
                            Self::regulatory_fragment_question_token(question),
                            left.strip_prefix(&format!("{}_", request.plan_id))
                                .unwrap_or(&left),
                            right
                                .strip_prefix(&format!("{}_", request.plan_id))
                                .unwrap_or(&right)
                        ),
                        question,
                        left_member_id: left,
                        right_member_id: right,
                        interpretation: interpretation.to_string(),
                    },
                });
            }
        };
        for question in &request.questions {
            match question {
                RegulatoryFragmentQuestion::StandaloneCandidate => push(
                    *question,
                    id("candidate_alone"),
                    id("control_minimal_promoter"),
                    "Compare the candidate-containing insert with the matched minimal-promoter control.",
                ),
                RegulatoryFragmentQuestion::PartnerDependence => {
                    push(
                        *question,
                        id("reference_combination"),
                        id("candidate_alone"),
                        "Compare the declared A+B geometry with candidate A alone.",
                    );
                    push(
                        *question,
                        id("reference_combination"),
                        id("partner_alone"),
                        "Compare the declared A+B geometry with partner B alone.",
                    );
                }
                RegulatoryFragmentQuestion::OrderDependence
                | RegulatoryFragmentQuestion::OrientationDependence
                | RegulatoryFragmentQuestion::SpacingDependence => {
                    let target_kind = match question {
                        RegulatoryFragmentQuestion::OrderDependence => {
                            RegulatoryFragmentGeometryKind::ReversedOrder
                        }
                        RegulatoryFragmentQuestion::OrientationDependence => {
                            RegulatoryFragmentGeometryKind::ReversedOrientation
                        }
                        _ => RegulatoryFragmentGeometryKind::ControlledSpacing,
                    };
                    for variant in request
                        .requested_variants
                        .iter()
                        .filter(|variant| variant.kind == target_kind)
                    {
                        push(
                            *question,
                            id(&format!("variant_{}", variant.variant_id)),
                            id("reference_combination"),
                            "Compare only the explicitly requested geometry with the declared reference combination.",
                        );
                    }
                }
                RegulatoryFragmentQuestion::MotifDisruption => {
                    if let Some(motif) = request.motif_disruption.as_ref() {
                        let left = id(&format!("motif_{}", motif.control_id));
                        let right =
                            if request
                                .reference_combination
                                .as_ref()
                                .is_some_and(|geometry| {
                                    geometry
                                        .instances
                                        .iter()
                                        .any(|instance| instance.fragment_id == motif.fragment_id)
                                })
                            {
                                id("reference_combination")
                            } else if request.fragments.iter().any(|fragment| {
                                fragment.fragment_id == motif.fragment_id
                                    && fragment.role == RegulatoryFragmentRole::Partner
                            }) {
                                id("partner_alone")
                            } else {
                                id("candidate_alone")
                            };
                        push(
                            *question,
                            left,
                            right,
                            "Compare the stated-rule motif edit with its exact wild-type geometry.",
                        );
                    }
                }
                RegulatoryFragmentQuestion::BoundaryUncertainty => {
                    for variant in request.requested_variants.iter().filter(|variant| {
                        matches!(
                            variant.kind,
                            RegulatoryFragmentGeometryKind::BoundaryShift
                                | RegulatoryFragmentGeometryKind::Tiling
                        )
                    }) {
                        push(
                            *question,
                            id(&format!("variant_{}", variant.variant_id)),
                            id("candidate_alone"),
                            "Compare the explicitly bounded alternate ROI geometry with the candidate-alone geometry.",
                        );
                    }
                }
            }
        }
        contrasts.sort_by(|left, right| left.contrast.contrast_id.cmp(&right.contrast.contrast_id));
        contrasts
    }

    fn select_regulatory_fragment_members(
        request: &RegulatoryFragmentPanelRequest,
        candidates: &[CandidateRegulatoryConstruct],
        contrasts: &[CandidateContrast],
    ) -> Vec<String> {
        let eligible = candidates
            .iter()
            .filter(|candidate| candidate.eligible)
            .collect::<Vec<_>>();
        let mandatory = eligible
            .iter()
            .filter(|candidate| candidate.mandatory)
            .map(|candidate| candidate.member.member_id.clone())
            .collect::<BTreeSet<_>>();
        let optional = eligible
            .iter()
            .filter(|candidate| !candidate.mandatory)
            .collect::<Vec<_>>();
        let mut best: Option<(usize, usize, Vec<usize>, Vec<String>)> = None;
        let mask_limit = 1u64 << optional.len();
        for mask in 0..mask_limit {
            let mut selected = mandatory.clone();
            for (index, candidate) in optional.iter().enumerate() {
                if mask & (1u64 << index) != 0 {
                    selected.insert(candidate.member.member_id.clone());
                }
            }
            if selected.len() > request.max_panel_members {
                continue;
            }
            let covered = Self::covered_regulatory_fragment_questions(
                &request.questions,
                contrasts,
                &selected,
            );
            let ordered_indices = eligible
                .iter()
                .enumerate()
                .filter(|(_, candidate)| selected.contains(&candidate.member.member_id))
                .map(|(index, _)| index)
                .collect::<Vec<_>>();
            let ordered_ids = ordered_indices
                .iter()
                .map(|index| eligible[*index].member.member_id.clone())
                .collect::<Vec<_>>();
            let score = (covered.len(), selected.len(), ordered_indices, ordered_ids);
            let replace = match &best {
                None => true,
                Some((best_covered, best_len, best_indices, _)) => {
                    score.0 > *best_covered
                        || (score.0 == *best_covered && score.1 < *best_len)
                        || (score.0 == *best_covered
                            && score.1 == *best_len
                            && score.2 < *best_indices)
                }
            };
            if replace {
                best = Some(score);
            }
        }
        best.map(|(_, _, _, ids)| ids).unwrap_or_else(|| {
            eligible
                .iter()
                .filter(|candidate| candidate.mandatory)
                .take(request.max_panel_members)
                .map(|candidate| candidate.member.member_id.clone())
                .collect()
        })
    }

    fn covered_regulatory_fragment_questions(
        requested: &[RegulatoryFragmentQuestion],
        contrasts: &[CandidateContrast],
        selected: &BTreeSet<String>,
    ) -> BTreeSet<RegulatoryFragmentQuestion> {
        requested
            .iter()
            .copied()
            .filter(|question| {
                let rows = contrasts
                    .iter()
                    .filter(|candidate| candidate.contrast.question == *question)
                    .collect::<Vec<_>>();
                !rows.is_empty()
                    && rows.iter().all(|candidate| {
                        selected.contains(&candidate.contrast.left_member_id)
                            && selected.contains(&candidate.contrast.right_member_id)
                    })
            })
            .collect()
    }

    fn attach_regulatory_fragment_inclusion_reasons(
        candidate: &mut CandidateRegulatoryConstruct,
        contrasts: &[RegulatoryFragmentContrast],
        request: &RegulatoryFragmentPanelRequest,
    ) {
        if candidate.mandatory {
            candidate
                .member
                .inclusion_reasons
                .push(RegulatoryFragmentInclusionReason {
                    kind: if candidate.member.construct_kind
                        == RegulatoryFragmentConstructKind::ReferenceControl
                    {
                        RegulatoryFragmentInclusionReasonKind::ExplicitReferenceControl
                    } else {
                        RegulatoryFragmentInclusionReasonKind::RequiredControl
                    },
                    detail: "The complete effective policy requires this explicit control."
                        .to_string(),
                    ..RegulatoryFragmentInclusionReason::default()
                });
        }
        for contrast in contrasts.iter().filter(|contrast| {
            contrast.left_member_id == candidate.member.member_id
                || contrast.right_member_id == candidate.member.member_id
        }) {
            candidate
                .member
                .inclusion_reasons
                .push(RegulatoryFragmentInclusionReason {
                    kind: RegulatoryFragmentInclusionReasonKind::ExperimentalComparisonEndpoint,
                    question: Some(contrast.question),
                    contrast_id: Some(contrast.contrast_id.clone()),
                    detail: "This member is an endpoint of an induced pair required by the requested question."
                        .to_string(),
                });
        }
        if candidate.member.inclusion_reasons.is_empty() {
            candidate
                .member
                .inclusion_reasons
                .push(RegulatoryFragmentInclusionReason {
                    kind: RegulatoryFragmentInclusionReasonKind::DeterministicCoverageChoice,
                    detail: format!(
                        "Retained by deterministic tie-breaking under max_panel_members={}.",
                        request.max_panel_members
                    ),
                    ..RegulatoryFragmentInclusionReason::default()
                });
        }
    }

    fn populate_regulatory_fragment_cloning_feasibility(
        &self,
        request: &RegulatoryFragmentPanelRequest,
        members: &mut [RegulatoryFragmentPanelMember],
    ) -> Result<Option<PromoterReporterPanelCloningStrategyReport>, EngineError> {
        let mut detached = self.fork_detached_execution();
        let mut insert_ids = vec![];
        for member in members.iter() {
            if member.insert_sequence_5prime_to_3prime.is_empty() {
                continue;
            }
            let insert_id = format!("regulatory_insert_{}", member.member_id);
            detached
                .engine_mut()
                .apply(Operation::CreateSequenceFromText {
                    sequence_text: member.insert_sequence_5prime_to_3prime.clone(),
                    output_id: Some(insert_id.clone()),
                    name: Some(member.member_id.clone()),
                    circular: false,
                })?;
            insert_ids.push(insert_id);
        }
        if insert_ids.is_empty() {
            return Ok(None);
        }
        let strategy = detached
            .engine()
            .restriction_cloning_panel_strategy(&request.vector_seq_id, &insert_ids)?;
        let state = match strategy.strategy {
            PromoterReporterPanelCloningStrategy::DirectionalRestriction => {
                RegulatoryFragmentCloningFeasibilityState::DirectionalRestrictionCandidate
            }
            PromoterReporterPanelCloningStrategy::Gibson => {
                RegulatoryFragmentCloningFeasibilityState::GibsonFallbackCandidate
            }
        };
        for member in members.iter_mut() {
            if !member.insert_sequence_5prime_to_3prime.is_empty() {
                member.cloning_feasibility = state;
            }
        }
        Ok(Some(strategy))
    }

    fn populate_regulatory_fragment_evidence_dimensions(
        &self,
        request: &RegulatoryFragmentPanelRequest,
        fragments: &[ResolvedRegulatoryFragment],
        members: &[RegulatoryFragmentPanelMember],
        cloning_strategy: Option<&PromoterReporterPanelCloningStrategyReport>,
    ) -> Result<Vec<RegulatoryFragmentEvidenceDimension>, EngineError> {
        let max_observations = request.policy.max_evidence_observations_per_dimension;
        let source_sequences = fragments
            .iter()
            .map(|fragment| fragment.projection.seq_id.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .filter_map(|seq_id| {
                self.state
                    .sequences
                    .get(&seq_id)
                    .map(|sequence| (seq_id, sequence.get_forward_string()))
            })
            .collect::<Vec<_>>();

        let (reference_observations, reference_truncated) =
            Self::regulatory_fragment_reference_uniqueness_observations(
                request,
                fragments,
                &source_sequences,
                max_observations,
            );
        let (similarity_observations, similarity_truncated) = self
            .regulatory_fragment_panel_similarity_observations(
                request,
                fragments,
                max_observations,
            )?;
        let (repeat_observations, repeat_truncated) =
            Self::regulatory_fragment_repeat_observations(members, max_observations)?;
        let (junction_observations, junction_truncated) =
            Self::regulatory_fragment_junction_observations(
                request,
                members,
                &source_sequences,
                max_observations,
            );
        let (cloning_observations, cloning_truncated) =
            Self::regulatory_fragment_cloning_observations(
                members,
                cloning_strategy,
                max_observations,
            );

        let mut dimensions = vec![
            Self::regulatory_fragment_evidence_dimension(
                request,
                RegulatoryFragmentEvidenceDimensionKind::ReferenceGenomicUniqueness,
                RegulatoryFragmentEvidenceState::Evaluated,
                vec!["compute_dotplot_shared_point_engine".to_string()],
                reference_observations,
                reference_truncated,
                "Exact full-fragment and mismatch-tolerant full-window matches were assessed against the loaded source sequences bound by the ROIs. This is reference-context evidence, not a whole-genome uniqueness claim.",
            )?,
            Self::regulatory_fragment_evidence_dimension(
                request,
                RegulatoryFragmentEvidenceDimensionKind::PanelSequenceSimilarity,
                RegulatoryFragmentEvidenceState::Evaluated,
                vec![
                    "align_sequences_global".to_string(),
                    "compute_dotplot_shared_point_engine".to_string(),
                ],
                similarity_observations,
                similarity_truncated,
                "Global alignment and exact forward/inverted word matches compare each bound fragment with the other fragments and the validated reporter vector. Similarity is a geometry/context warning, not an activity verdict.",
            )?,
            Self::regulatory_fragment_evidence_dimension(
                request,
                RegulatoryFragmentEvidenceDimensionKind::RepeatsAndLowComplexity,
                RegulatoryFragmentEvidenceState::Evaluated,
                vec!["construct_reasoning_sequence_similarity".to_string()],
                repeat_observations,
                repeat_truncated,
                "The existing construct-reasoning repeat and low-complexity scanner was applied to each selected non-empty insert. An empty observation list means no configured pattern was detected; it is not a general experimental pass.",
            )?,
            Self::regulatory_fragment_evidence_dimension(
                request,
                RegulatoryFragmentEvidenceDimensionKind::PairSpecificJunctionUniqueness,
                RegulatoryFragmentEvidenceState::Evaluated,
                vec!["exact_junction_context_scan".to_string()],
                junction_observations,
                junction_truncated,
                "Each selected multi-fragment junction was searched exactly, in both orientations, against bound source sequences and other selected inserts. This bounded context is not whole-genome uniqueness.",
            )?,
            Self::regulatory_fragment_evidence_dimension(
                request,
                RegulatoryFragmentEvidenceDimensionKind::RestrictionAndCloningRisk,
                if cloning_strategy.is_some() {
                    RegulatoryFragmentEvidenceState::Evaluated
                } else {
                    RegulatoryFragmentEvidenceState::NotEvaluated
                },
                if cloning_strategy.is_some() {
                    vec!["promoter_reporter_panel_cloning_strategy".to_string()]
                } else {
                    vec![]
                },
                cloning_observations,
                cloning_truncated,
                if cloning_strategy.is_some() {
                    "Existing shared directional-restriction/Gibson feasibility was retained per insert. Feasibility is an operational planning result, not evidence of regulatory behavior."
                } else {
                    "No non-empty selected insert was available for cloning assessment; not_evaluated is never a pass."
                },
            )?,
        ];
        for kind in [
            RegulatoryFragmentEvidenceDimensionKind::EnsemblRegulatoryOverlap,
            RegulatoryFragmentEvidenceDimensionKind::TfbsModelScoreContext,
            RegulatoryFragmentEvidenceDimensionKind::CutrunAndChromatinContext,
        ] {
            dimensions.push(self.regulatory_fragment_external_dimension(request, fragments, kind)?);
        }
        Ok(dimensions)
    }

    fn regulatory_fragment_evidence_dimension(
        request: &RegulatoryFragmentPanelRequest,
        kind: RegulatoryFragmentEvidenceDimensionKind,
        state: RegulatoryFragmentEvidenceState,
        mut method_ids: Vec<String>,
        observations: Vec<RegulatoryFragmentEvidenceObservation>,
        truncated: bool,
        detail: &str,
    ) -> Result<RegulatoryFragmentEvidenceDimension, EngineError> {
        method_ids.sort();
        method_ids.dedup();
        let bindings = request
            .evidence_bindings
            .iter()
            .filter(|binding| binding.dimension == kind)
            .cloned()
            .collect::<Vec<_>>();
        let assessment_id = format!(
            "{}:{}",
            request.plan_id,
            Self::regulatory_fragment_evidence_dimension_token(kind)
        );
        let (blockers, mut warnings) =
            Self::regulatory_fragment_evidence_findings(kind, &observations);
        if truncated {
            warnings.push(RegulatoryFragmentFinding {
                code: "evidence_observation_limit_reached".to_string(),
                subject_ids: vec![assessment_id.clone()],
                detail: format!(
                    "The deterministic evidence output stopped at max_evidence_observations_per_dimension={}.",
                    request.policy.max_evidence_observations_per_dimension
                ),
            });
        }
        let assessment_sha256 = Self::regulatory_fragment_value_sha256(
            &serde_json::json!({
                "assessment_id": assessment_id,
                "kind": kind,
                "state": state,
                "method_ids": method_ids,
                "bindings": bindings,
                "observations": observations,
                "blockers": blockers,
                "warnings": warnings,
                "truncated": truncated,
                "detail": detail,
            }),
            "regulatory-fragment evidence assessment",
        )?;
        Ok(RegulatoryFragmentEvidenceDimension {
            kind,
            state,
            assessment_id,
            assessment_sha256,
            method_ids,
            bindings,
            observations,
            blockers,
            warnings,
            truncated,
            detail: detail.to_string(),
        })
    }

    fn regulatory_fragment_evidence_findings(
        kind: RegulatoryFragmentEvidenceDimensionKind,
        observations: &[RegulatoryFragmentEvidenceObservation],
    ) -> (
        Vec<RegulatoryFragmentFinding>,
        Vec<RegulatoryFragmentFinding>,
    ) {
        let blockers = vec![];
        let mut warnings = vec![];
        for observation in observations {
            match observation {
                RegulatoryFragmentEvidenceObservation::PanelSequenceSimilarity {
                    left_subject_id,
                    right_subject_id,
                    comparison_scope,
                    exact_forward_word_match_count,
                    exact_inverted_word_match_count,
                    word_match_counts_truncated,
                    ..
                } if kind == RegulatoryFragmentEvidenceDimensionKind::PanelSequenceSimilarity => {
                    if comparison_scope == "fragment_to_reporter_vector"
                        && (*exact_forward_word_match_count > 0
                            || *exact_inverted_word_match_count > 0)
                    {
                        warnings.push(RegulatoryFragmentFinding {
                            code: "fragment_reporter_vector_word_similarity".to_string(),
                            subject_ids: vec![left_subject_id.clone(), right_subject_id.clone()],
                            detail: "Exact local word reuse with the reporter vector was detected and should be reviewed as construct-context evidence."
                                .to_string(),
                        });
                    }
                    if *exact_inverted_word_match_count > 0 {
                        warnings.push(RegulatoryFragmentFinding {
                            code: "inverted_sequence_similarity_detected".to_string(),
                            subject_ids: vec![left_subject_id.clone(), right_subject_id.clone()],
                            detail: "Exact reverse-complement word reuse was detected; this remains separate from direct similarity and repeat evidence."
                                .to_string(),
                        });
                    }
                    if *word_match_counts_truncated {
                        warnings.push(RegulatoryFragmentFinding {
                            code: "sequence_similarity_word_matches_truncated".to_string(),
                            subject_ids: vec![left_subject_id.clone(), right_subject_id.clone()],
                            detail: "The shared dotplot point limit was reached, so word-match counts are lower bounds."
                                .to_string(),
                        });
                    }
                }
                RegulatoryFragmentEvidenceObservation::RepeatOrLowComplexity {
                    subject_id,
                    evidence,
                    ..
                } if kind == RegulatoryFragmentEvidenceDimensionKind::RepeatsAndLowComplexity => {
                    warnings.push(RegulatoryFragmentFinding {
                        code: "repeat_or_low_complexity_context_detected".to_string(),
                        subject_ids: vec![subject_id.clone(), evidence.evidence_id.clone()],
                        detail: format!("{}: {}", evidence.label, evidence.rationale),
                    });
                }
                RegulatoryFragmentEvidenceObservation::PairSpecificJunctionUniqueness {
                    member_id,
                    observation_id,
                    exact_unique_in_assessed_context: false,
                    ..
                } if kind
                    == RegulatoryFragmentEvidenceDimensionKind::PairSpecificJunctionUniqueness =>
                {
                    warnings.push(RegulatoryFragmentFinding {
                        code: "junction_not_unique_in_assessed_context".to_string(),
                        subject_ids: vec![member_id.clone(), observation_id.clone()],
                        detail: "The exact junction window also occurs in a bound source sequence or another selected insert."
                            .to_string(),
                    });
                }
                RegulatoryFragmentEvidenceObservation::RestrictionAndCloningRisk {
                    member_id,
                    blocker_count,
                    ..
                } if kind == RegulatoryFragmentEvidenceDimensionKind::RestrictionAndCloningRisk
                    && *blocker_count > 0 =>
                {
                    warnings.push(RegulatoryFragmentFinding {
                        code: "directional_restriction_pair_conflict".to_string(),
                        subject_ids: vec![member_id.clone()],
                        detail: format!(
                            "Existing cloning assessment found {blocker_count} insert/enzyme conflict(s); the panel-wide strategy may use another pair or Gibson fallback."
                        ),
                    });
                }
                _ => {}
            }
        }
        warnings.sort_by(|left, right| {
            left.code
                .cmp(&right.code)
                .then(left.subject_ids.cmp(&right.subject_ids))
        });
        warnings.dedup_by(|left, right| {
            left.code == right.code && left.subject_ids == right.subject_ids
        });
        (blockers, warnings)
    }

    fn regulatory_fragment_reference_uniqueness_observations(
        request: &RegulatoryFragmentPanelRequest,
        fragments: &[ResolvedRegulatoryFragment],
        source_sequences: &[(String, String)],
        max_observations: usize,
    ) -> (Vec<RegulatoryFragmentEvidenceObservation>, bool) {
        let mut observations = vec![];
        let mut truncated = false;
        'fragments: for fragment in fragments {
            for (source_seq_id, source_sequence) in source_sequences {
                if observations.len() >= max_observations {
                    truncated = true;
                    break 'fragments;
                }
                let exact_forward = Self::regulatory_fragment_exact_match_count(
                    source_sequence,
                    &fragment.canonical_sequence,
                );
                let reverse = Self::reverse_complement(&fragment.canonical_sequence);
                let exact_reverse =
                    Self::regulatory_fragment_exact_match_count(source_sequence, &reverse);
                let near_forward = Self::regulatory_fragment_full_window_match_count(
                    &fragment.binding.fragment_id,
                    &fragment.canonical_sequence,
                    source_seq_id,
                    source_sequence,
                    DotplotMode::PairForward,
                    request.policy.near_exact_max_mismatches,
                );
                let near_reverse = Self::regulatory_fragment_full_window_match_count(
                    &fragment.binding.fragment_id,
                    &fragment.canonical_sequence,
                    source_seq_id,
                    source_sequence,
                    DotplotMode::PairReverseComplement,
                    request.policy.near_exact_max_mismatches,
                );
                let near_state = if near_forward.is_some() && near_reverse.is_some() {
                    RegulatoryFragmentEvidenceState::Evaluated
                } else {
                    RegulatoryFragmentEvidenceState::Unavailable
                };
                observations.push(
                    RegulatoryFragmentEvidenceObservation::ReferenceGenomicUniqueness {
                        observation_id: format!(
                            "reference_uniqueness:{}:{}",
                            fragment.binding.fragment_id, source_seq_id
                        ),
                        fragment_id: fragment.binding.fragment_id.clone(),
                        source_seq_id: source_seq_id.clone(),
                        fragment_length_bp: fragment.canonical_sequence.len(),
                        exact_forward_match_count: exact_forward,
                        exact_reverse_complement_match_count: exact_reverse,
                        near_exact_max_mismatches: request.policy.near_exact_max_mismatches,
                        near_exact_forward_match_count: near_forward,
                        near_exact_reverse_complement_match_count: near_reverse,
                        near_exact_state: near_state,
                        detail: if near_state == RegulatoryFragmentEvidenceState::Evaluated {
                            "Counts use full-fragment windows over this loaded ROI-bound source sequence."
                                .to_string()
                        } else {
                            "Exact counts were evaluated, but the bounded shared mismatch-aware dotplot engine declined the near-exact scan at this source size."
                                .to_string()
                        },
                    },
                );
            }
        }
        (observations, truncated)
    }

    fn regulatory_fragment_panel_similarity_observations(
        &self,
        request: &RegulatoryFragmentPanelRequest,
        fragments: &[ResolvedRegulatoryFragment],
        max_observations: usize,
    ) -> Result<(Vec<RegulatoryFragmentEvidenceObservation>, bool), EngineError> {
        let vector_sequence = self
            .state
            .sequences
            .get(&request.vector_seq_id)
            .map(DNAsequence::get_forward_string)
            .ok_or_else(|| {
                Self::regulatory_fragment_error(
                    "vector_sequence_missing",
                    [&request.vector_seq_id],
                    "The validated reporter vector is no longer loaded during evidence assessment.",
                )
            })?;
        let mut subjects = fragments
            .iter()
            .map(|fragment| {
                (
                    fragment.binding.fragment_id.clone(),
                    fragment.canonical_sequence.clone(),
                    "fragment",
                )
            })
            .collect::<Vec<_>>();
        subjects.push((
            request.vector_seq_id.clone(),
            vector_sequence,
            "reporter_vector",
        ));
        let mut observations = vec![];
        let mut truncated = false;
        'outer: for left_index in 0..subjects.len() {
            for right_index in left_index + 1..subjects.len() {
                if observations.len() >= max_observations {
                    truncated = true;
                    break 'outer;
                }
                let (left_id, left_sequence, left_kind) = &subjects[left_index];
                let (right_id, right_sequence, right_kind) = &subjects[right_index];
                let alignment = Self::compute_pairwise_alignment_report(
                    left_id,
                    left_sequence,
                    None,
                    None,
                    right_id,
                    right_sequence,
                    None,
                    None,
                    PairwiseAlignmentMode::Global,
                    2,
                    -3,
                    -5,
                    -1,
                )?
                .report;
                let word_size = request
                    .policy
                    .sequence_word_size_bp
                    .min(left_sequence.len())
                    .min(right_sequence.len());
                let (forward_words, inverted_words, word_match_counts_truncated) = if word_size == 0
                {
                    (0, 0, false)
                } else {
                    let (forward, forward_truncated) =
                        Self::regulatory_fragment_dotplot_match_count(
                            left_id,
                            left_sequence,
                            right_id,
                            right_sequence,
                            DotplotMode::PairForward,
                            word_size,
                            0,
                        )?;
                    let (inverted, inverted_truncated) =
                        Self::regulatory_fragment_dotplot_match_count(
                            left_id,
                            left_sequence,
                            right_id,
                            right_sequence,
                            DotplotMode::PairReverseComplement,
                            word_size,
                            0,
                        )?;
                    (forward, inverted, forward_truncated || inverted_truncated)
                };
                observations.push(
                    RegulatoryFragmentEvidenceObservation::PanelSequenceSimilarity {
                        observation_id: format!("panel_similarity:{left_id}:{right_id}"),
                        left_subject_id: left_id.clone(),
                        right_subject_id: right_id.clone(),
                        comparison_scope: if *left_kind == "reporter_vector"
                            || *right_kind == "reporter_vector"
                        {
                            "fragment_to_reporter_vector".to_string()
                        } else {
                            "fragment_to_fragment".to_string()
                        },
                        word_size_bp: word_size,
                        exact_forward_word_match_count: forward_words,
                        exact_inverted_word_match_count: inverted_words,
                        word_match_counts_truncated,
                        alignment,
                        detail: "Global similarity and local exact-word reuse are reported separately; neither predicts regulatory activity."
                            .to_string(),
                    },
                );
            }
        }
        Ok((observations, truncated))
    }

    fn regulatory_fragment_repeat_observations(
        members: &[RegulatoryFragmentPanelMember],
        max_observations: usize,
    ) -> Result<(Vec<RegulatoryFragmentEvidenceObservation>, bool), EngineError> {
        let mut observations = vec![];
        let mut truncated = false;
        'members: for member in members {
            if member.insert_sequence_5prime_to_3prime.is_empty() {
                continue;
            }
            let dna = DNAsequence::from_sequence(&member.insert_sequence_5prime_to_3prime)
                .map_err(|error| {
                    Self::regulatory_fragment_error(
                        "evidence_sequence_invalid",
                        [&member.member_id],
                        format!("Could not inspect the planned insert sequence: {error}"),
                    )
                })?;
            let mut rows = Self::build_construct_reasoning_sequence_similarity_evidence(
                &member.member_id,
                &dna,
            );
            rows.sort_by(|left, right| left.evidence_id.cmp(&right.evidence_id));
            for evidence in rows {
                if observations.len() >= max_observations {
                    truncated = true;
                    break 'members;
                }
                observations.push(
                    RegulatoryFragmentEvidenceObservation::RepeatOrLowComplexity {
                        observation_id: format!(
                            "repeat_context:{}:{}",
                            member.member_id, evidence.evidence_id
                        ),
                        subject_id: member.member_id.clone(),
                        evidence,
                    },
                );
            }
        }
        Ok((observations, truncated))
    }

    fn regulatory_fragment_junction_observations(
        request: &RegulatoryFragmentPanelRequest,
        members: &[RegulatoryFragmentPanelMember],
        source_sequences: &[(String, String)],
        max_observations: usize,
    ) -> (Vec<RegulatoryFragmentEvidenceObservation>, bool) {
        let mut observations = vec![];
        let mut truncated = false;
        'members: for member in members {
            for (index, pair) in member.instances.windows(2).enumerate() {
                if observations.len() >= max_observations {
                    truncated = true;
                    break 'members;
                }
                let left = &pair[0];
                let right = &pair[1];
                let start = left
                    .assembled_end_0based_exclusive
                    .saturating_sub(request.policy.junction_flank_bp);
                let end = right
                    .assembled_start_0based
                    .saturating_add(request.policy.junction_flank_bp)
                    .min(member.insert_sequence_5prime_to_3prime.len());
                let junction = &member.insert_sequence_5prime_to_3prime[start..end];
                let reverse = Self::reverse_complement(junction);
                let reference_forward_match_count = source_sequences
                    .iter()
                    .map(|(_, sequence)| {
                        Self::regulatory_fragment_exact_match_count(sequence, junction)
                    })
                    .sum();
                let reference_reverse_complement_match_count = source_sequences
                    .iter()
                    .map(|(_, sequence)| {
                        Self::regulatory_fragment_exact_match_count(sequence, &reverse)
                    })
                    .sum();
                let panel_other_member_match_count = members
                    .iter()
                    .filter(|other| other.member_id != member.member_id)
                    .map(|other| {
                        Self::regulatory_fragment_exact_match_count(
                            &other.insert_sequence_5prime_to_3prime,
                            junction,
                        ) + Self::regulatory_fragment_exact_match_count(
                            &other.insert_sequence_5prime_to_3prime,
                            &reverse,
                        )
                    })
                    .sum();
                let exact_unique_in_assessed_context = reference_forward_match_count == 0
                    && reference_reverse_complement_match_count == 0
                    && panel_other_member_match_count == 0;
                observations.push(
                    RegulatoryFragmentEvidenceObservation::PairSpecificJunctionUniqueness {
                        observation_id: format!("junction:{}:{}", member.member_id, index + 1),
                        member_id: member.member_id.clone(),
                        left_fragment_id: left.fragment_id.clone(),
                        right_fragment_id: right.fragment_id.clone(),
                        assembled_start_0based: start,
                        assembled_end_0based_exclusive: end,
                        junction_sequence_sha256: sha256_prefixed_str(junction),
                        reference_forward_match_count,
                        reference_reverse_complement_match_count,
                        panel_other_member_match_count,
                        exact_unique_in_assessed_context,
                        detail: "The assayed window contains both fragment flanks and any declared spacer; the current construct itself is excluded from off-target counts."
                            .to_string(),
                    },
                );
            }
        }
        (observations, truncated)
    }

    fn regulatory_fragment_cloning_observations(
        members: &[RegulatoryFragmentPanelMember],
        cloning_strategy: Option<&PromoterReporterPanelCloningStrategyReport>,
        max_observations: usize,
    ) -> (Vec<RegulatoryFragmentEvidenceObservation>, bool) {
        let Some(strategy) = cloning_strategy else {
            return (vec![], false);
        };
        let summary_by_id = strategy
            .insert_site_summaries
            .iter()
            .map(|summary| (summary.insert_seq_id.as_str(), summary))
            .collect::<HashMap<_, _>>();
        let mut observations = vec![];
        let mut truncated = false;
        for member in members {
            if member.insert_sequence_5prime_to_3prime.is_empty() {
                continue;
            }
            if observations.len() >= max_observations {
                truncated = true;
                break;
            }
            let insert_id = format!("regulatory_insert_{}", member.member_id);
            let site_count_by_enzyme = summary_by_id
                .get(insert_id.as_str())
                .map(|summary| summary.site_count_by_enzyme.clone())
                .unwrap_or_default();
            let blocker_count = strategy
                .pair_evaluations
                .iter()
                .flat_map(|evaluation| &evaluation.blockers)
                .filter(|blocker| blocker.insert_seq_id == insert_id)
                .count();
            observations.push(
                RegulatoryFragmentEvidenceObservation::RestrictionAndCloningRisk {
                    observation_id: format!("cloning_risk:{}", member.member_id),
                    member_id: member.member_id.clone(),
                    cloning_feasibility: member.cloning_feasibility,
                    site_count_by_enzyme,
                    blocker_count,
                    selected_strategy: strategy.strategy,
                    detail: "Counts and strategy come from the existing panel-wide restriction-cloning feasibility helper."
                        .to_string(),
                },
            );
        }
        (observations, truncated)
    }

    fn regulatory_fragment_full_window_match_count(
        query_id: &str,
        query: &str,
        reference_id: &str,
        reference: &str,
        mode: DotplotMode,
        max_mismatches: usize,
    ) -> Option<usize> {
        if query.is_empty() || reference.len() < query.len() {
            return Some(0);
        }
        Self::regulatory_fragment_dotplot_match_count(
            query_id,
            query,
            reference_id,
            reference,
            mode,
            query.len(),
            max_mismatches,
        )
        .ok()
        .and_then(|(count, truncated)| (!truncated).then_some(count))
    }

    fn regulatory_fragment_dotplot_match_count(
        _query_id: &str,
        query: &str,
        _reference_id: &str,
        reference: &str,
        mode: DotplotMode,
        word_size: usize,
        max_mismatches: usize,
    ) -> Result<(usize, bool), EngineError> {
        Self::compute_dotplot_points(
            query.as_bytes(),
            reference.as_bytes(),
            0,
            0,
            mode,
            word_size,
            1,
            max_mismatches,
            MAX_DOTPLOT_POINTS,
        )
        .map(|(points, truncated)| (points.len(), truncated))
    }

    fn regulatory_fragment_exact_match_count(haystack: &str, needle: &str) -> usize {
        if needle.is_empty() || haystack.len() < needle.len() {
            return 0;
        }
        let haystack = haystack.as_bytes();
        let needle = needle.as_bytes();
        haystack
            .windows(needle.len())
            .filter(|window| {
                window
                    .iter()
                    .zip(needle)
                    .all(|(left, right)| left.eq_ignore_ascii_case(right))
            })
            .count()
    }

    fn regulatory_fragment_questions_for_construct(
        request: &RegulatoryFragmentPanelRequest,
        kind: RegulatoryFragmentConstructKind,
        variant_id: &str,
    ) -> BTreeSet<RegulatoryFragmentQuestion> {
        request
            .questions
            .iter()
            .copied()
            .filter(|question| match (kind, question) {
                (
                    RegulatoryFragmentConstructKind::CandidateAlone
                    | RegulatoryFragmentConstructKind::MinimalPromoterControl,
                    RegulatoryFragmentQuestion::StandaloneCandidate,
                ) => true,
                (
                    RegulatoryFragmentConstructKind::CandidateAlone
                    | RegulatoryFragmentConstructKind::PartnerAlone
                    | RegulatoryFragmentConstructKind::ReferenceCombination,
                    RegulatoryFragmentQuestion::PartnerDependence,
                ) => true,
                (
                    RegulatoryFragmentConstructKind::MotifDisruptionControl,
                    RegulatoryFragmentQuestion::MotifDisruption,
                ) => true,
                (RegulatoryFragmentConstructKind::RequestedGeometryVariant, question) => {
                    request.requested_variants.iter().any(|variant| {
                        variant.variant_id == variant_id
                            && matches!(
                                (variant.kind, question),
                                (
                                    RegulatoryFragmentGeometryKind::ReversedOrder,
                                    RegulatoryFragmentQuestion::OrderDependence
                                ) | (
                                    RegulatoryFragmentGeometryKind::ReversedOrientation,
                                    RegulatoryFragmentQuestion::OrientationDependence
                                ) | (
                                    RegulatoryFragmentGeometryKind::ControlledSpacing,
                                    RegulatoryFragmentQuestion::SpacingDependence
                                ) | (
                                    RegulatoryFragmentGeometryKind::BoundaryShift
                                        | RegulatoryFragmentGeometryKind::Tiling,
                                    RegulatoryFragmentQuestion::BoundaryUncertainty
                                )
                            )
                    })
                }
                _ => false,
            })
            .collect()
    }

    fn regulatory_fragment_contrasts_by_member(
        contrasts: &[RegulatoryFragmentContrast],
    ) -> BTreeMap<String, Vec<String>> {
        let mut out = BTreeMap::<String, Vec<String>>::new();
        for contrast in contrasts {
            for member_id in [&contrast.left_member_id, &contrast.right_member_id] {
                out.entry(member_id.clone())
                    .or_default()
                    .push(contrast.contrast_id.clone());
            }
        }
        for ids in out.values_mut() {
            ids.sort();
            ids.dedup();
        }
        out
    }

    fn regulatory_fragment_role_purpose_warnings(
        binding: &RegulatoryFragmentBinding,
    ) -> Vec<RegulatoryFragmentFinding> {
        let compatible = match binding.role {
            RegulatoryFragmentRole::Candidate | RegulatoryFragmentRole::Partner => matches!(
                binding.region.purpose,
                gp::GenomicRegionPurpose::CandidateCisRegulatoryRegion
                    | gp::GenomicRegionPurpose::ReporterCandidate
            ),
            RegulatoryFragmentRole::MinimalPromoter => matches!(
                binding.region.purpose,
                gp::GenomicRegionPurpose::PromoterRegion
                    | gp::GenomicRegionPurpose::ReporterCandidate
            ),
            RegulatoryFragmentRole::ReferenceControl => matches!(
                binding.region.purpose,
                gp::GenomicRegionPurpose::PromoterRegion
                    | gp::GenomicRegionPurpose::ReporterCandidate
                    | gp::GenomicRegionPurpose::CandidateCisRegulatoryRegion
            ),
        };
        if compatible {
            vec![]
        } else {
            vec![RegulatoryFragmentFinding {
                code: "declared_role_region_purpose_mismatch".to_string(),
                subject_ids: vec![
                    binding.fragment_id.clone(),
                    binding.region.region_id.clone(),
                ],
                detail: format!(
                    "Declared role '{:?}' is unusual for persisted ROI purpose '{}'; retained as a non-blocking review warning.",
                    binding.role,
                    binding.region.purpose.as_str()
                ),
            }]
        }
    }

    fn regulatory_fragment_construct_role_rank(kind: RegulatoryFragmentConstructKind) -> usize {
        match kind {
            RegulatoryFragmentConstructKind::PromoterlessControl => 0,
            RegulatoryFragmentConstructKind::MinimalPromoterControl => 1,
            RegulatoryFragmentConstructKind::ReferenceControl => 2,
            RegulatoryFragmentConstructKind::CandidateAlone => 3,
            RegulatoryFragmentConstructKind::PartnerAlone => 4,
            RegulatoryFragmentConstructKind::ReferenceCombination => 5,
            RegulatoryFragmentConstructKind::RequestedGeometryVariant => 6,
            RegulatoryFragmentConstructKind::MotifDisruptionControl => 7,
        }
    }

    fn regulatory_fragment_question_token(question: RegulatoryFragmentQuestion) -> &'static str {
        match question {
            RegulatoryFragmentQuestion::StandaloneCandidate => "standalone_candidate",
            RegulatoryFragmentQuestion::PartnerDependence => "partner_dependence",
            RegulatoryFragmentQuestion::OrderDependence => "order_dependence",
            RegulatoryFragmentQuestion::OrientationDependence => "orientation_dependence",
            RegulatoryFragmentQuestion::SpacingDependence => "spacing_dependence",
            RegulatoryFragmentQuestion::MotifDisruption => "motif_disruption",
            RegulatoryFragmentQuestion::BoundaryUncertainty => "boundary_uncertainty",
        }
    }

    fn regulatory_fragment_evidence_dimension_token(
        kind: RegulatoryFragmentEvidenceDimensionKind,
    ) -> &'static str {
        match kind {
            RegulatoryFragmentEvidenceDimensionKind::ReferenceGenomicUniqueness => {
                "reference_genomic_uniqueness"
            }
            RegulatoryFragmentEvidenceDimensionKind::PanelSequenceSimilarity => {
                "panel_sequence_similarity"
            }
            RegulatoryFragmentEvidenceDimensionKind::RepeatsAndLowComplexity => {
                "repeats_and_low_complexity"
            }
            RegulatoryFragmentEvidenceDimensionKind::PairSpecificJunctionUniqueness => {
                "pair_specific_junction_uniqueness"
            }
            RegulatoryFragmentEvidenceDimensionKind::RestrictionAndCloningRisk => {
                "restriction_and_cloning_risk"
            }
            RegulatoryFragmentEvidenceDimensionKind::EnsemblRegulatoryOverlap => {
                "ensembl_regulatory_overlap"
            }
            RegulatoryFragmentEvidenceDimensionKind::TfbsModelScoreContext => {
                "tfbs_model_score_context"
            }
            RegulatoryFragmentEvidenceDimensionKind::CutrunAndChromatinContext => {
                "cutrun_and_chromatin_context"
            }
        }
    }

    fn regulatory_fragment_projection_status_token(
        status: gp::GenomicRegionLocalProjectionStatus,
    ) -> &'static str {
        match status {
            gp::GenomicRegionLocalProjectionStatus::Current => "current",
            gp::GenomicRegionLocalProjectionStatus::SequenceUnavailable => "sequence_unavailable",
            gp::GenomicRegionLocalProjectionStatus::SequenceDigestMismatch => {
                "sequence_digest_mismatch"
            }
            gp::GenomicRegionLocalProjectionStatus::AnchorMismatch => "anchor_mismatch",
            gp::GenomicRegionLocalProjectionStatus::OutsideSequence => "outside_sequence",
        }
    }

    fn normalize_regulatory_fragment_dna(raw: &str, label: &str) -> Result<String, EngineError> {
        let normalized = raw
            .bytes()
            .filter(|byte| !byte.is_ascii_whitespace())
            .map(|byte| match byte.to_ascii_uppercase() {
                b'U' => b'T',
                other => other,
            })
            .collect::<Vec<_>>();
        if !normalized
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
        {
            return Err(Self::regulatory_fragment_error(
                "invalid_geometry_spacer",
                [label],
                "Spacer DNA may contain only A, C, G, T/U, N, and whitespace.",
            ));
        }
        String::from_utf8(normalized).map_err(|error| {
            Self::regulatory_fragment_error(
                "invalid_geometry_spacer",
                [label],
                format!("Spacer is not valid UTF-8: {error}"),
            )
        })
    }

    fn regulatory_fragment_slug(raw: &str) -> String {
        raw.trim()
            .to_ascii_lowercase()
            .chars()
            .map(|character| {
                if character.is_ascii_alphanumeric() {
                    character
                } else {
                    '_'
                }
            })
            .collect::<String>()
            .split('_')
            .filter(|part| !part.is_empty())
            .collect::<Vec<_>>()
            .join("_")
    }

    fn regulatory_fragment_valid_sha256(value: &str) -> bool {
        value
            .strip_prefix("sha256:")
            .is_some_and(|hex| hex.len() == 64 && hex.bytes().all(|byte| byte.is_ascii_hexdigit()))
    }

    fn regulatory_fragment_value_sha256<T: Serialize>(
        value: &T,
        label: &str,
    ) -> Result<String, EngineError> {
        let canonical = Self::construct_reasoning_canonical_json(value, label)?;
        let portable =
            serde_json::from_str::<serde_json::Value>(&canonical).map_err(|error| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Could not normalize regulatory-fragment {label} through portable JSON: {error}"
                ),
                cause_chain: vec![],
            })?;
        serde_json::to_string(&portable)
            .map(|portable| sha256_prefixed_str(&portable))
            .map_err(|error| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Could not serialize normalized regulatory-fragment {label}: {error}"
                ),
                cause_chain: vec![],
            })
    }

    fn regulatory_fragment_panel_proposal_digest(
        plan: &RegulatoryFragmentPanelPlan,
    ) -> Result<String, EngineError> {
        let mut basis = plan.clone();
        basis.proposal_digest.clear();
        Self::regulatory_fragment_value_sha256(&basis, "regulatory-fragment panel proposal")
    }

    fn regulatory_fragment_error<I, S>(
        code: &str,
        subject_ids: I,
        detail: impl Into<String>,
    ) -> EngineError
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let ids = subject_ids
            .into_iter()
            .map(|id| id.as_ref().trim().to_string())
            .filter(|id| !id.is_empty())
            .collect::<Vec<_>>();
        EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("{code}: subject_ids=[{}]: {}", ids.join(","), detail.into()),
            cause_chain: vec![],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence::DNAsequence;
    use std::path::{Path, PathBuf};
    use tempfile::{TempDir, tempdir};

    const CANDIDATE_SEQUENCE: &str = "ACGTTGCAAGTCCTGATCGATGCTAGCATCGTACGATTCGGAATCCGTCAGTACCG";
    const PARTNER_SEQUENCE: &str = "TTAGGCCATACCTTGGACTACGGCATGATCCGTAGGCTTAACCGTTCGATGACCAT";
    const MINIMAL_PROMOTER_SEQUENCE: &str = "TATAAACTGCGCGTCCGATCAGTA";

    struct PlannerFixture {
        _temp: TempDir,
        engine: GentleEngine,
        request: RegulatoryFragmentPanelRequest,
    }

    fn test_reference(assembly_name: &str) -> gp::GenomicRegionReference {
        gp::GenomicRegionReference {
            species_scientific_name: Some("Homo sapiens".to_string()),
            taxon_id: Some(9606),
            assembly_name: assembly_name.to_string(),
            assembly_accession: Some("GCA_000001405.15".to_string()),
            contig_name: "chr7".to_string(),
            contig_accession: Some("NC_000007.14".to_string()),
            contig_aliases: vec!["7".to_string()],
        }
    }

    fn write_helper_catalog(temp: &TempDir) -> PathBuf {
        let fixture_path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("test_files/fixtures/reporter_vectors/synthetic_mcs_backbone.gb");
        let path = temp.path().join("helper_vectors.json");
        let catalog = serde_json::json!({
            "Synthetic panel vector": {
                "description": "Repository-owned synthetic promoter-reporter test vector",
                "sequence_local": fixture_path.to_string_lossy(),
                "annotations_local": fixture_path.to_string_lossy(),
                "usable_as_empty_backbone": true,
                "helper_kind": "plasmid_vector",
                "sequence_expectation": {
                    "schema": crate::genomes::HELPER_VECTOR_SEQUENCE_EXPECTATION_SCHEMA,
                    "provider": "GENtle tests",
                    "product_name": "synthetic MCS backbone",
                    "catalog_number": "SYNTH-MCS-1",
                    "accession_version": "GENTLE_SYNTHETIC_MCS.1",
                    "expected_length_bp": 240,
                    "expected_topology": "circular",
                    "required_features": [
                        {
                            "id": "multiple_cloning_region",
                            "feature_kinds": ["misc_feature"],
                            "qualifier_terms": ["multiple cloning site region"],
                            "expected_start_1based": 1,
                            "expected_end_1based": 70
                        },
                        {
                            "id": "luc2",
                            "feature_kinds": ["CDS"],
                            "qualifier_terms": ["luciferase luc2 marker"],
                            "expected_start_1based": 100,
                            "expected_end_1based": 180
                        }
                    ],
                    "restriction_site_equivalences": [],
                    "provenance": [{
                        "source_id": "synthetic-test-fixture",
                        "source_url": "test_files/fixtures/reporter_vectors/synthetic_mcs_backbone.gb",
                        "asserted_on": "2026-08-28",
                        "note": "Repository-owned deterministic fixture."
                    }]
                }
            }
        });
        fs::write(
            &path,
            serde_json::to_vec_pretty(&catalog).expect("helper catalog JSON"),
        )
        .expect("write helper catalog");
        path
    }

    fn add_anchored_sequence(
        engine: &mut GentleEngine,
        seq_id: &str,
        sequence: &str,
        start_1based: u64,
    ) {
        let mut dna = DNAsequence::from_sequence(sequence).expect("synthetic DNA");
        GentleEngine::prepare_sequence(&mut dna);
        engine.state_mut().sequences.insert(seq_id.to_string(), dna);
        let provenance = engine
            .state_mut()
            .metadata
            .entry(PROVENANCE_METADATA_KEY.to_string())
            .or_insert_with(|| serde_json::json!({"genome_extractions": []}));
        provenance["genome_extractions"]
            .as_array_mut()
            .expect("genome extraction array")
            .push(serde_json::json!({
                "seq_id": seq_id,
                "genome_id": "GRCh38",
                "chromosome": "chr7",
                "start_1based": start_1based,
                "end_1based": start_1based + sequence.len() as u64 - 1,
                "anchor_strand": "+",
                "anchor_verified": true,
                "recorded_at_unix_ms": 1_706_000_000_000_u128
            }));
    }

    fn capture_fragment(
        engine: &mut GentleEngine,
        set_id: &str,
        seq_id: &str,
        region_id: &str,
        purpose: gp::GenomicRegionPurpose,
        strand: gp::GenomicRegionStrand,
    ) -> gp::GenomicRegionOfInterest {
        let sequence_len = engine
            .state()
            .sequences
            .get(seq_id)
            .expect("source sequence")
            .len();
        engine
            .apply(Operation::CaptureGenomicRegion {
                request: gp::GenomicRegionCaptureRequest {
                    set_id: set_id.to_string(),
                    set_label: Some("Synthetic regulatory panel regions".to_string()),
                    region_id: Some(region_id.to_string()),
                    label: Some(region_id.to_string()),
                    purpose,
                    source: gp::GenomicRegionCaptureSource::SequenceSelection {
                        seq_id: seq_id.to_string(),
                        local_start_0based: 0,
                        local_end_0based_exclusive: sequence_len as u64,
                        strand,
                        reference_override: Some(test_reference("GRCh38")),
                    },
                    created_at_unix_ms: Some(1_706_000_000_000),
                    collision_policy: gp::GenomicRegionCollisionPolicy::Reject,
                    ..gp::GenomicRegionCaptureRequest::default()
                },
            })
            .expect("capture synthetic ROI")
            .genomic_region_operation
            .expect("genomic-region report")
            .region
            .expect("captured ROI")
    }

    fn planner_fixture(
        include_partner: bool,
        candidate_strand: gp::GenomicRegionStrand,
        highly_similar_partner: bool,
    ) -> PlannerFixture {
        let temp = tempdir().expect("regulatory panel test directory");
        let helper_catalog_path = write_helper_catalog(&temp);
        let fixture_path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("test_files/fixtures/reporter_vectors/synthetic_mcs_backbone.gb");
        let mut vector = crate::dna_sequence::load_from_file(&fixture_path.to_string_lossy())
            .expect("load synthetic panel vector");
        GentleEngine::prepare_sequence(&mut vector);

        let mut engine = GentleEngine::default();
        engine
            .state_mut()
            .sequences
            .insert("synthetic_panel_vector".to_string(), vector);
        add_anchored_sequence(&mut engine, "candidate_source", CANDIDATE_SEQUENCE, 1_001);
        let partner_sequence = if highly_similar_partner {
            "ACGTTGCAAGTCCTGATCGATGCTAGCATCGTACGATTCGGAATCCGTCAGTTCCT"
        } else {
            PARTNER_SEQUENCE
        };
        if include_partner {
            add_anchored_sequence(&mut engine, "partner_source", partner_sequence, 2_001);
        }
        add_anchored_sequence(
            &mut engine,
            "minimal_source",
            MINIMAL_PROMOTER_SEQUENCE,
            3_001,
        );

        let candidate = capture_fragment(
            &mut engine,
            "candidate_regions",
            "candidate_source",
            "candidate_roi",
            gp::GenomicRegionPurpose::ReporterCandidate,
            candidate_strand,
        );
        let partner = include_partner.then(|| {
            capture_fragment(
                &mut engine,
                "partner_regions",
                "partner_source",
                "partner_roi",
                gp::GenomicRegionPurpose::CandidateCisRegulatoryRegion,
                gp::GenomicRegionStrand::Plus,
            )
        });
        let minimal = capture_fragment(
            &mut engine,
            "minimal_promoter_regions",
            "minimal_source",
            "minimal_promoter_roi",
            gp::GenomicRegionPurpose::PromoterRegion,
            gp::GenomicRegionStrand::Plus,
        );
        let store = engine
            .genomic_region_store_snapshot()
            .expect("region store");
        let binding = |fragment_id: &str,
                       declared_order: usize,
                       role: RegulatoryFragmentRole,
                       region: gp::GenomicRegionOfInterest| {
            let set = store
                .sets
                .iter()
                .find(|set| {
                    set.regions
                        .iter()
                        .any(|stored| stored.region_id == region.region_id)
                })
                .expect("region set containing fixture ROI");
            RegulatoryFragmentBinding {
                fragment_id: fragment_id.to_string(),
                declared_order,
                role,
                region_set_id: set.set_id.clone(),
                region_set_content_sha256: set.content_sha256.clone(),
                reference_release: "Ensembl 116".to_string(),
                region,
            }
        };
        let mut fragments = vec![binding(
            "candidate_a",
            0,
            RegulatoryFragmentRole::Candidate,
            candidate,
        )];
        if let Some(partner) = partner {
            fragments.push(binding(
                "partner_b",
                1,
                RegulatoryFragmentRole::Partner,
                partner,
            ));
        }
        fragments.push(binding(
            "minimal_promoter",
            2,
            RegulatoryFragmentRole::MinimalPromoter,
            minimal,
        ));
        let reference_combination = include_partner.then(|| RegulatoryFragmentGeometryRequest {
            variant_id: "reference_ab".to_string(),
            declared_order: 3,
            kind: RegulatoryFragmentGeometryKind::ReferenceCombination,
            instances: vec![
                RegulatoryFragmentInstanceRequest {
                    fragment_id: "candidate_a".to_string(),
                    orientation: RegulatoryFragmentOrientation::Forward,
                    spacer_before: String::new(),
                },
                RegulatoryFragmentInstanceRequest {
                    fragment_id: "partner_b".to_string(),
                    orientation: RegulatoryFragmentOrientation::Forward,
                    spacer_before: "GG".to_string(),
                },
                RegulatoryFragmentInstanceRequest {
                    fragment_id: "minimal_promoter".to_string(),
                    orientation: RegulatoryFragmentOrientation::Forward,
                    spacer_before: String::new(),
                },
            ],
            ..RegulatoryFragmentGeometryRequest::default()
        });
        let mut questions = vec![RegulatoryFragmentQuestion::StandaloneCandidate];
        if include_partner {
            questions.push(RegulatoryFragmentQuestion::PartnerDependence);
        }
        PlannerFixture {
            _temp: temp,
            engine,
            request: RegulatoryFragmentPanelRequest {
                plan_id: "synthetic_regulatory_panel".to_string(),
                vector_seq_id: "synthetic_panel_vector".to_string(),
                vector_catalog_id: "Synthetic panel vector".to_string(),
                helper_catalog_path: Some(helper_catalog_path.to_string_lossy().to_string()),
                fragments,
                questions,
                reference_combination,
                mutation_policy: PromoterReporterPanelMutationPolicy::NativeOnlyV1,
                ..RegulatoryFragmentPanelRequest::default()
            },
        }
    }

    fn member_kinds(
        plan: &RegulatoryFragmentPanelPlan,
    ) -> BTreeSet<RegulatoryFragmentConstructKind> {
        plan.members
            .iter()
            .map(|member| member.construct_kind)
            .collect()
    }

    fn mutate_bound_region(
        fixture: &mut PlannerFixture,
        fragment_id: &str,
        mutate: impl FnOnce(&mut gp::GenomicRegionOfInterest),
    ) {
        let mut store = fixture
            .engine
            .genomic_region_store_snapshot()
            .expect("region store");
        let binding_index = fixture
            .request
            .fragments
            .iter()
            .position(|fragment| fragment.fragment_id == fragment_id)
            .expect("fragment binding");
        let region_id = fixture.request.fragments[binding_index]
            .region
            .region_id
            .clone();
        let set = store
            .sets
            .iter_mut()
            .find(|set| set.set_id == fixture.request.fragments[binding_index].region_set_id)
            .expect("region set");
        let region = set
            .regions
            .iter_mut()
            .find(|region| region.region_id == region_id)
            .expect("persisted region");
        mutate(region);
        crate::engine::genomic_regions::recompute_region_digests(region)
            .expect("recompute region digests");
        let rebound_region = region.clone();
        crate::engine::genomic_regions::recompute_set_digest(set).expect("recompute set digest");
        fixture.request.fragments[binding_index].region = rebound_region;
        let set_digest = set.content_sha256.clone();
        for binding in &mut fixture.request.fragments {
            if binding.region_set_id == set.set_id {
                binding.region_set_content_sha256 = set_digest.clone();
            }
        }
        fixture.engine.state_mut().metadata.insert(
            GENOMIC_REGION_SETS_METADATA_KEY.to_string(),
            serde_json::to_value(store).expect("serialize region store"),
        );
    }

    #[test]
    fn regulatory_fragment_a_only_yields_candidate_and_explicit_controls() {
        let mut fixture = planner_fixture(false, gp::GenomicRegionStrand::Plus, false);
        let result = fixture
            .engine
            .apply(Operation::PlanRegulatoryFragmentPanel {
                request: Box::new(fixture.request.clone()),
                path: None,
            })
            .expect("plan through shared engine operation");
        let plan = result
            .regulatory_fragment_panel_plan
            .expect("regulatory-fragment plan");
        assert_eq!(plan.schema, REGULATORY_FRAGMENT_PANEL_PLAN_SCHEMA);
        assert_eq!(plan.members.len(), 3);
        assert_eq!(
            member_kinds(&plan),
            BTreeSet::from([
                RegulatoryFragmentConstructKind::PromoterlessControl,
                RegulatoryFragmentConstructKind::MinimalPromoterControl,
                RegulatoryFragmentConstructKind::CandidateAlone,
            ])
        );
        assert!(plan.uncovered_questions.is_empty());
        assert_eq!(
            plan.planning_label,
            RegulatoryFragmentPlanningLabel::StandaloneTestableCandidate
        );
        assert!(
            plan.members
                .iter()
                .all(|member| !member.inclusion_reasons.is_empty())
        );
        assert!(plan.materialization_supported);
    }

    #[test]
    fn regulatory_fragment_plan_operation_writes_the_exact_json_report() {
        let mut fixture = planner_fixture(false, gp::GenomicRegionStrand::Plus, false);
        let path = fixture._temp.path().join("regulatory_panel.json");
        let result = fixture
            .engine
            .apply(Operation::PlanRegulatoryFragmentPanel {
                request: Box::new(fixture.request.clone()),
                path: Some(path.to_string_lossy().to_string()),
            })
            .expect("plan and export through shared engine operation");
        let planned = result
            .regulatory_fragment_panel_plan
            .expect("regulatory-fragment plan");
        let exported: RegulatoryFragmentPanelPlan =
            serde_json::from_slice(&fs::read(&path).expect("exported regulatory-fragment plan"))
                .expect("parse exported regulatory-fragment plan");
        assert_eq!(exported.proposal_digest, planned.proposal_digest);
        assert_eq!(exported.request_sha256, planned.request_sha256);
    }

    #[test]
    fn regulatory_fragment_shell_routes_plan_and_render_exact_artifacts() {
        let mut fixture = planner_fixture(true, gp::GenomicRegionStrand::Plus, false);
        let request_path = fixture._temp.path().join("request.json");
        let plan_path = fixture._temp.path().join("plan.json");
        let svg_path = fixture._temp.path().join("plan.svg");
        fs::write(
            &request_path,
            serde_json::to_vec_pretty(&fixture.request).expect("request JSON"),
        )
        .expect("write request JSON");
        let plan_command = crate::engine_shell::parse_shell_line(&format!(
            "promoters regulatory-panel-plan @{} --path {}",
            request_path.display(),
            plan_path.display()
        ))
        .expect("parse regulatory panel shell plan");
        let out = crate::engine_shell::execute_shell_command(&mut fixture.engine, &plan_command)
            .expect("execute regulatory panel shell plan");
        assert!(!out.state_changed);
        let result_digest = out.output["result"]["proposal_digest"]
            .as_str()
            .expect("result proposal digest")
            .to_string();
        let plan: RegulatoryFragmentPanelPlan =
            serde_json::from_slice(&fs::read(&plan_path).expect("shell plan artifact"))
                .expect("parse shell plan artifact");
        assert_eq!(plan.proposal_digest, result_digest);
        GentleEngine::validate_regulatory_fragment_panel_document(&plan)
            .expect("exported plan retains a valid portable digest basis");

        let render_command = crate::engine_shell::parse_shell_line(&format!(
            "promoters regulatory-panel-render @{} --path {}",
            plan_path.display(),
            svg_path.display()
        ))
        .expect("parse regulatory panel shell render");
        let out = crate::engine_shell::execute_shell_command(&mut fixture.engine, &render_command)
            .expect("execute regulatory panel shell render");
        assert!(!out.state_changed);
        assert_eq!(
            out.output["result"]["proposal_digest"].as_str(),
            Some(result_digest.as_str())
        );
        assert!(
            fs::read_to_string(svg_path)
                .expect("shell SVG artifact")
                .contains("Independent evidence lanes")
        );
    }

    #[test]
    fn regulatory_fragment_ab_uses_minimal_identifiable_contrast_panel() {
        let fixture = planner_fixture(true, gp::GenomicRegionStrand::Plus, false);
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("plan A+B panel");
        assert_eq!(plan.members.len(), 5);
        assert_eq!(plan.contrasts.len(), 3);
        assert!(plan.uncovered_questions.is_empty());
        assert!(
            member_kinds(&plan).contains(&RegulatoryFragmentConstructKind::ReferenceCombination)
        );
        assert!(
            !member_kinds(&plan)
                .contains(&RegulatoryFragmentConstructKind::RequestedGeometryVariant)
        );
    }

    #[test]
    fn regulatory_fragment_svg_operation_renders_exact_plan_and_rejects_tampering() {
        let mut fixture = planner_fixture(true, gp::GenomicRegionStrand::Plus, false);
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request.clone())
            .expect("regulatory-fragment plan");
        let svg_path = fixture._temp.path().join("panel.svg");
        fixture
            .engine
            .apply(Operation::RenderRegulatoryFragmentPanelSvg {
                plan: Box::new(plan.clone()),
                path: svg_path.to_string_lossy().to_string(),
            })
            .expect("render exact plan");
        let svg = fs::read_to_string(&svg_path).expect("rendered SVG");
        assert!(svg.contains("Genome-anchored source fragments"));
        assert!(svg.contains("Independent evidence lanes"));

        let mut tampered = plan;
        tampered.members[0].member_id.push_str("_changed");
        let rejected_path = fixture._temp.path().join("tampered.svg");
        let error = fixture
            .engine
            .apply(Operation::RenderRegulatoryFragmentPanelSvg {
                plan: Box::new(tampered),
                path: rejected_path.to_string_lossy().to_string(),
            })
            .expect_err("tampered plan must not render");
        assert!(error.message.contains("proposal_content_digest_mismatch"));
        assert!(!rejected_path.exists());
    }

    #[test]
    fn regulatory_fragment_geometry_variants_are_strictly_opt_in() {
        let mut fixture = planner_fixture(true, gp::GenomicRegionStrand::Plus, false);
        fixture
            .request
            .questions
            .push(RegulatoryFragmentQuestion::OrderDependence);
        fixture
            .request
            .requested_variants
            .push(RegulatoryFragmentGeometryRequest {
                variant_id: "reverse_order_only".to_string(),
                declared_order: 4,
                kind: RegulatoryFragmentGeometryKind::ReversedOrder,
                instances: vec![
                    RegulatoryFragmentInstanceRequest {
                        fragment_id: "partner_b".to_string(),
                        ..RegulatoryFragmentInstanceRequest::default()
                    },
                    RegulatoryFragmentInstanceRequest {
                        fragment_id: "candidate_a".to_string(),
                        spacer_before: "GG".to_string(),
                        ..RegulatoryFragmentInstanceRequest::default()
                    },
                    RegulatoryFragmentInstanceRequest {
                        fragment_id: "minimal_promoter".to_string(),
                        ..RegulatoryFragmentInstanceRequest::default()
                    },
                ],
                ..RegulatoryFragmentGeometryRequest::default()
            });
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("plan requested order variant");
        let variants = plan
            .members
            .iter()
            .filter(|member| {
                member.construct_kind == RegulatoryFragmentConstructKind::RequestedGeometryVariant
            })
            .collect::<Vec<_>>();
        assert_eq!(variants.len(), 1);
        assert!(
            variants[0]
                .member_id
                .ends_with("variant_reverse_order_only")
        );
        assert!(
            plan.contrasts.iter().any(|contrast| {
                contrast.question == RegulatoryFragmentQuestion::OrderDependence
            })
        );
        assert!(!plan.contrasts.iter().any(|contrast| {
            matches!(
                contrast.question,
                RegulatoryFragmentQuestion::OrientationDependence
                    | RegulatoryFragmentQuestion::SpacingDependence
            )
        }));
    }

    #[test]
    fn regulatory_fragment_high_similarity_is_a_geometry_confound_not_a_verdict() {
        let fixture = planner_fixture(true, gp::GenomicRegionStrand::Plus, true);
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("plan similar A+B panel");
        assert_eq!(
            plan.planning_label,
            RegulatoryFragmentPlanningLabel::ContextOrGeometryConfounded
        );
        assert!(
            plan.candidate_partner_similarity
                .as_ref()
                .is_some_and(|audit| audit.highly_similar)
        );
        assert!(
            plan.warnings
                .iter()
                .any(|warning| { warning.code == "candidate_partner_high_similarity" })
        );
        assert!(plan.nonclaims.iter().all(|text| !text.contains("proves")));
    }

    #[test]
    fn regulatory_fragment_materialization_is_exact_approved_and_atomic() {
        let mut fixture = planner_fixture(true, gp::GenomicRegionStrand::Minus, false);
        fixture
            .request
            .reference_combination
            .as_mut()
            .expect("combination")
            .instances[1]
            .orientation = RegulatoryFragmentOrientation::ReverseComplement;
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("plan");
        let before = serde_json::to_value(fixture.engine.snapshot()).expect("baseline");
        let proposal = fixture
            .engine
            .plan_regulatory_fragment_materialization(plan.clone(), "design".into())
            .expect("proposal");
        assert_eq!(
            serde_json::to_value(fixture.engine.snapshot()).expect("state"),
            before
        );
        assert_eq!(
            serde_json::to_value(&proposal).expect("proposal"),
            serde_json::to_value(
                fixture
                    .engine
                    .plan_regulatory_fragment_materialization(plan.clone(), "design".into())
                    .expect("repeat")
            )
            .expect("proposal")
        );
        assert!(proposal.products.iter().any(|p| p.instances.len() >= 3));
        let plan_json = serde_json::to_string(&plan).expect("plan JSON");
        let command = crate::engine_shell::parse_shell_line(&format!(
            "promoters regulatory-products-plan '{plan_json}' --output-prefix design"
        ))
        .expect("shell parse");
        let result = crate::engine_shell::execute_shell_command(&mut fixture.engine, &command)
            .expect("shared shell");
        assert!(!result.state_changed);
        let actual = &result.output["result"]["regulatory_fragment_materialization_proposal"];
        // Match the existing portable-JSON digest convention, including its
        // float round-trip. Approval below consumes the transported object.
        let expected: serde_json::Value =
            serde_json::from_str(&serde_json::to_string(&proposal).expect("proposal JSON"))
                .expect("portable proposal");
        assert!(
            actual == &expected,
            "shared shell changed the portable proposal"
        );
        let proposal: RegulatoryFragmentMaterializationProposal =
            serde_json::from_str(&serde_json::to_string(actual).expect("returned JSON"))
                .expect("transported proposal");
        let vector_id = plan.vector_context.vector_seq_id.clone();
        let original_features = fixture.engine.state.sequences[&vector_id]
            .features()
            .clone();
        fixture
            .engine
            .state
            .sequences
            .get_mut(&vector_id)
            .expect("vector")
            .features_mut()[0]
            .qualifiers
            .push(("note".into(), Some("changed after review".into())));
        let modified_state = serde_json::to_value(fixture.engine.snapshot()).expect("state");
        assert!(
            fixture
                .engine
                .materialize_regulatory_fragment_panel(proposal.clone(), &proposal.proposal_digest)
                .is_err()
        );
        assert_eq!(
            serde_json::to_value(fixture.engine.snapshot()).expect("state"),
            modified_state
        );
        *fixture
            .engine
            .state
            .sequences
            .get_mut(&vector_id)
            .expect("vector")
            .features_mut() = original_features;
        for product in &proposal.products {
            let member = plan
                .members
                .iter()
                .find(|m| m.member_id == product.member_id)
                .expect("member");
            let start = plan.vector_context.insertion_start_0based;
            assert_eq!(
                &product.sequence_5prime_to_3prime[start..start + member.insert_length_bp],
                member.insert_sequence_5prime_to_3prime
            );
            assert_eq!(product.instances, member.instances);
            assert!(
                product
                    .features
                    .iter()
                    .filter(|f| f
                        .qualifiers
                        .iter()
                        .any(|(k, _)| k == "gentle_fragment_instance"))
                    .count()
                    == member.instances.len()
            );
        }
        let mut tampered = proposal.clone();
        tampered.products.reverse();
        assert!(
            fixture
                .engine
                .materialize_regulatory_fragment_panel(tampered, &proposal.proposal_digest)
                .is_err()
        );
        assert!(
            fixture
                .engine
                .materialize_regulatory_fragment_panel(proposal.clone(), "wrong")
                .is_err()
        );
        assert_eq!(
            serde_json::to_value(fixture.engine.snapshot()).expect("state"),
            before
        );
        // A collision in the last member must not leak earlier products.
        let collision = proposal
            .products
            .last()
            .expect("last")
            .output_seq_id
            .clone();
        fixture.engine.state.sequences.insert(
            collision.clone(),
            DNAsequence::from_sequence("ACGT").expect("DNA"),
        );
        let collision_state = serde_json::to_value(fixture.engine.snapshot()).expect("state");
        assert!(
            fixture
                .engine
                .materialize_regulatory_fragment_panel(proposal.clone(), &proposal.proposal_digest)
                .is_err()
        );
        assert_eq!(
            serde_json::to_value(fixture.engine.snapshot()).expect("state"),
            collision_state
        );
        fixture.engine.state.sequences.remove(&collision);
        let result = fixture
            .engine
            .apply(Operation::MaterializeRegulatoryFragmentPanel {
                approval_digest: proposal.proposal_digest.clone(),
                proposal: Box::new(proposal.clone()),
            })
            .expect("approved atomic operation");
        let receipt = result
            .regulatory_fragment_materialization_receipt
            .expect("receipt");
        assert_eq!(receipt.created_seq_ids.len(), plan.members.len());
        assert_eq!(
            receipt.final_product_audit_state,
            RegulatoryFragmentFinalProductAuditState::NotEvaluated
        );
        for product in &proposal.products {
            let actual = &fixture.engine.state.sequences[&product.output_seq_id];
            assert_eq!(
                actual.get_forward_string(),
                product.sequence_5prime_to_3prime
            );
            assert_eq!(actual.features(), &product.features);
        }
        fixture.engine.undo_last_operation().expect("one undo");
        assert_eq!(
            serde_json::to_value(fixture.engine.snapshot()).expect("state"),
            before
        );
    }

    #[test]
    fn regulatory_fragment_products_mcp_preserves_approval_and_results() {
        let fixture = planner_fixture(false, gp::GenomicRegionStrand::Plus, false);
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("plan");
        let state_path = fixture._temp.path().join("mcp_state.json");
        fixture
            .engine
            .snapshot()
            .save_to_path(&state_path.to_string_lossy())
            .expect("save state");
        let before = fs::read(&state_path).expect("state bytes");
        let call = |operation: &Operation, confirm: bool| {
            crate::mcp_server::mcp_tool_call_for_capability_surface_tests(
                &state_path.to_string_lossy(),
                "op",
                serde_json::json!({"operation": operation, "confirm": confirm}),
            )
        };
        let planned = call(
            &Operation::PlanRegulatoryFragmentMaterialization {
                plan: Box::new(plan),
                output_prefix: "mcp_design".into(),
            },
            true,
        );
        assert_eq!(planned["isError"], false, "{planned:#}");
        let proposal: RegulatoryFragmentMaterializationProposal = serde_json::from_value(
            planned["structuredContent"]["result"]["regulatory_fragment_materialization_proposal"]
                .clone(),
        )
        .expect("typed MCP proposal");
        // Generic MCP op saves after success. Reloading rebuilds this unordered
        // restriction cache; compare its content without mistaking row order for a mutation.
        let normalized = |bytes: &[u8]| {
            let mut state: serde_json::Value = serde_json::from_slice(bytes).expect("state JSON");
            for dna in state["sequences"]
                .as_object_mut()
                .expect("sequences")
                .values_mut()
            {
                dna["restriction_enzyme_groups"]
                    .as_array_mut()
                    .expect("cache groups")
                    .sort_by_cached_key(|row| serde_json::to_string(row).expect("cache row"));
            }
            state
        };
        let after_plan = fs::read(&state_path).expect("state bytes");
        assert!(
            normalized(&after_plan) == normalized(&before),
            "planning changed project content"
        );
        let before = after_plan;
        let operation = Operation::MaterializeRegulatoryFragmentPanel {
            proposal: Box::new(proposal.clone()),
            approval_digest: proposal.proposal_digest.clone(),
        };
        let denied = call(&operation, false);
        assert_eq!(denied["isError"], true, "{denied:#}");
        assert!(
            fs::read(&state_path).expect("state bytes") == before,
            "unconfirmed materialization wrote state"
        );
        let wrong_digest = call(
            &Operation::MaterializeRegulatoryFragmentPanel {
                proposal: Box::new(proposal.clone()),
                approval_digest: "wrong".into(),
            },
            true,
        );
        assert_eq!(wrong_digest["isError"], true, "{wrong_digest:#}");
        assert!(
            fs::read(&state_path).expect("state bytes") == before,
            "rejected digest wrote state"
        );
        let materialized = call(&operation, true);
        assert_eq!(materialized["isError"], false, "{materialized:#}");
        let receipt: RegulatoryFragmentMaterializationReceipt = serde_json::from_value(
            materialized["structuredContent"]["result"]["regulatory_fragment_materialization_receipt"].clone()
        ).expect("typed MCP receipt");
        assert_eq!(receipt.approved_proposal_digest, proposal.proposal_digest);
        assert_eq!(receipt.created_seq_ids.len(), proposal.products.len());
        assert_eq!(
            receipt.final_product_audit_state,
            RegulatoryFragmentFinalProductAuditState::NotEvaluated
        );
        let saved =
            ProjectState::load_from_path(&state_path.to_string_lossy()).expect("saved products");
        for product in proposal.products {
            let dna = &saved.sequences[&product.output_seq_id];
            assert_eq!(dna.get_forward_string(), product.sequence_5prime_to_3prime);
            assert_eq!(dna.features(), &product.features);
        }
    }

    #[test]
    fn regulatory_fragment_external_reports_are_bound_and_keep_unavailable_distinct() {
        let mut fixture = planner_fixture(false, gp::GenomicRegionStrand::Minus, false);
        let fragment = &fixture.request.fragments[0];
        let projection = fragment
            .region
            .local_projection
            .as_ref()
            .expect("projection");
        let seq_id = projection.seq_id.clone();
        let dna = &fixture.engine.state.sequences[&seq_id];
        let anchor = fixture
            .engine
            .sequence_genome_anchor_summary(&seq_id)
            .expect("anchor");
        let mut locus = gp::GeneLocusEvidenceDisplayReport {
            schema: gp::GENE_LOCUS_EVIDENCE_DISPLAY_SCHEMA.into(),
            panel_id: "external_panel".into(),
            seq_id: seq_id.clone(),
            sequence_binding: Some(crate::locus_report::sequence_binding(dna, Some(&anchor))),
            locus_local_start_1based: 1,
            locus_local_end_1based: dna.len(),
            isoform_evidence: gp::GeneIsoformEvidenceReport {
                assembly: fragment.region.interval.reference.assembly_name.clone(),
                annotation_release: Some(fragment.reference_release.clone()),
                ..Default::default()
            },
            regulatory_score_tracks: vec![gp::GeneLocusRegulatoryScoreTrack {
                track_id: "model".into(),
                state: gp::GeneLocusRegulatoryScoreState::NotAssessable,
                ..Default::default()
            }],
            ..Default::default()
        };
        let path = fixture._temp.path().join("locus.json");
        let write = |locus: &gp::GeneLocusEvidenceDisplayReport| {
            let bytes = serde_json::to_vec(locus).expect("synthetic typed source");
            fs::write(&path, &bytes).expect("write source");
            crate::digest_utils::sha256_prefixed_bytes(&bytes)
        };
        fixture.request.evidence_bindings = vec![RegulatoryFragmentEvidenceBinding {
            dimension: RegulatoryFragmentEvidenceDimensionKind::TfbsModelScoreContext,
            report_id: locus.panel_id.clone(),
            report_sha256: write(&locus),
            report_path: Some(path.to_string_lossy().into_owned()),
            row_id: Some("model".into()),
        }];
        let before = serde_json::to_value(fixture.engine.snapshot()).expect("state");
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request.clone())
            .expect("unavailable report");
        let dimension = &plan.evidence_dimensions[6];
        assert_eq!(
            dimension.state,
            RegulatoryFragmentEvidenceState::NotEvaluated
        );
        assert!(matches!(
            &dimension.observations[0],
            RegulatoryFragmentEvidenceObservation::ExternalLocusContext {
                state: RegulatoryFragmentEvidenceState::Unavailable,
                ..
            }
        ));
        let track = &mut locus.regulatory_score_tracks[0];
        track.state = gp::GeneLocusRegulatoryScoreState::Available;
        track.input_sequence_id = seq_id;
        track.input_sequence_sha256 = locus
            .sequence_binding
            .as_ref()
            .expect("binding")
            .sequence_sha256
            .clone();
        track.assembly = locus.isoform_evidence.assembly.clone();
        track.chromosome = anchor.chromosome.clone();
        track.anchor_start_1based = anchor.start_1based;
        track.anchor_end_1based = anchor.end_1based;
        track.window_length_bp = 4;
        track.stride_bp = 1;
        track.forward_scores = vec![1.5, 2.5];
        fixture.request.evidence_bindings[0].report_sha256 = write(&locus);
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request.clone())
            .expect("available report");
        assert_eq!(
            plan.evidence_dimensions[6].state,
            RegulatoryFragmentEvidenceState::Evaluated
        );
        assert_eq!(
            serde_json::to_value(fixture.engine.snapshot()).expect("state"),
            before
        );
        // The same envelope retains distinct annotation and experimental-signal payloads.
        let source = gp::EnsemblRegulationSourceDescriptor {
            source_id: "synthetic_regulation".into(),
            assembly_name: locus.isoform_evidence.assembly.clone(),
            annotation_release: "synthetic regulatory release".into(),
            ..Default::default()
        };
        locus.ensembl_regulation = Some(gp::GeneLocusEnsemblRegulationEvidence {
            availability: gp::GeneLocusEnsemblRegulationAvailability::Available,
            requested_source_id: source.source_id.clone(),
            source: Some(source.clone()),
            source_binding: Some(gp::GeneLocusEnsemblRegulationSourceBinding {
                source,
                content_identity_verified: true,
                index_sha256: format!("sha256:{}", "1".repeat(64)),
                intervals_sha256: format!("sha256:{}", "2".repeat(64)),
                ..Default::default()
            }),
            rows: vec![gp::GeneLocusEnsemblRegulationFeatureRow {
                feature_id: "regulatory_feature".into(),
                assembly_name: locus.isoform_evidence.assembly.clone(),
                displayed_local_start_1based: 3,
                displayed_local_end_1based: 12,
                ..Default::default()
            }],
            ..Default::default()
        });
        locus.occupancy_groups = vec![gp::GeneLocusOccupancyGroup {
            group_id: "synthetic_cells".into(),
            lanes: vec![gp::GeneLocusOccupancyLane {
                lane: gp::GeneIsoformOccupancyLane {
                    lane_id: "synthetic_cutrun".into(),
                    interval_count: 1,
                    intervals: vec![gp::GeneIsoformOccupancyInterval {
                        interval_id: "peak".into(),
                        local_start_1based: 5,
                        local_end_1based: 15,
                        score: Some(4.25),
                        ..Default::default()
                    }],
                    ..Default::default()
                },
                source_sha256: Some(format!("sha256:{}", "3".repeat(64))),
                source_assembly: Some(locus.isoform_evidence.assembly.clone()),
                ..Default::default()
            }],
            ..Default::default()
        }];
        for (kind, row) in [
            (
                RegulatoryFragmentEvidenceDimensionKind::EnsemblRegulatoryOverlap,
                "regulatory_feature",
            ),
            (
                RegulatoryFragmentEvidenceDimensionKind::CutrunAndChromatinContext,
                "synthetic_cutrun",
            ),
        ] {
            let mut binding = fixture.request.evidence_bindings[0].clone();
            binding.dimension = kind;
            binding.row_id = Some(row.into());
            fixture.request.evidence_bindings.push(binding);
        }
        let digest = write(&locus);
        for binding in &mut fixture.request.evidence_bindings {
            binding.report_sha256 = digest.clone();
        }
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request.clone())
            .expect("three typed lanes");
        for kind in [
            RegulatoryFragmentEvidenceDimensionKind::EnsemblRegulatoryOverlap,
            RegulatoryFragmentEvidenceDimensionKind::TfbsModelScoreContext,
            RegulatoryFragmentEvidenceDimensionKind::CutrunAndChromatinContext,
        ] {
            let lane = plan
                .evidence_dimensions
                .iter()
                .find(|lane| lane.kind == kind)
                .expect("lane");
            assert_eq!(lane.state, RegulatoryFragmentEvidenceState::Evaluated);
            assert_eq!(lane.observations.len(), 1);
        }
        // Exact bytes, source sequence, assembly, release and row selection fail independently.
        for field in [
            "bytes",
            "sequence",
            "assembly",
            "release",
            "row",
            "score_bounds",
            "ensembl_bounds",
            "occupancy_hash",
        ] {
            let mut changed = locus.clone();
            let mut request = fixture.request.clone();
            match field {
                "bytes" => changed.warnings.push("edited".into()),
                "sequence" => {
                    changed
                        .sequence_binding
                        .as_mut()
                        .expect("binding")
                        .sequence_sha256 = format!("sha256:{}", "0".repeat(64))
                }
                "assembly" => changed.isoform_evidence.assembly = "other".into(),
                "release" => changed.isoform_evidence.annotation_release = Some("other".into()),
                "row" => request.evidence_bindings[0].row_id = Some("missing".into()),
                "score_bounds" => {
                    changed.regulatory_score_tracks[0].track_start_0based = usize::MAX
                }
                "ensembl_bounds" => {
                    changed
                        .ensembl_regulation
                        .as_mut()
                        .expect("annotation")
                        .rows[0]
                        .displayed_local_start_1based = 0
                }
                "occupancy_hash" => changed.occupancy_groups[0].lanes[0].source_sha256 = None,
                _ => unreachable!(),
            }
            let digest = write(&changed);
            if field != "bytes" {
                for binding in &mut request.evidence_bindings {
                    binding.report_sha256 = digest.clone();
                }
            }
            assert!(
                fixture
                    .engine
                    .plan_regulatory_fragment_panel(request)
                    .is_err(),
                "{field}"
            );
        }
    }

    #[test]
    fn regulatory_fragment_evidence_lanes_are_populated_independently() {
        let mut fixture = planner_fixture(false, gp::GenomicRegionStrand::Plus, false);
        let digest = format!("sha256:{}", "0".repeat(64));
        fixture.request.evidence_bindings = vec![
            RegulatoryFragmentEvidenceBinding {
                dimension: RegulatoryFragmentEvidenceDimensionKind::PanelSequenceSimilarity,
                report_id: "vector_similarity".to_string(),
                report_sha256: digest.clone(),
                ..RegulatoryFragmentEvidenceBinding::default()
            },
            RegulatoryFragmentEvidenceBinding {
                dimension: RegulatoryFragmentEvidenceDimensionKind::RepeatsAndLowComplexity,
                report_id: "inverted_repeats".to_string(),
                report_sha256: digest.clone(),
                ..RegulatoryFragmentEvidenceBinding::default()
            },
            RegulatoryFragmentEvidenceBinding {
                dimension: RegulatoryFragmentEvidenceDimensionKind::PairSpecificJunctionUniqueness,
                report_id: "junction_uniqueness".to_string(),
                report_sha256: digest,
                ..RegulatoryFragmentEvidenceBinding::default()
            },
        ];
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("plan with evidence citations");
        assert_eq!(plan.evidence_dimensions.len(), 8);
        for evaluated in [
            RegulatoryFragmentEvidenceDimensionKind::ReferenceGenomicUniqueness,
            RegulatoryFragmentEvidenceDimensionKind::PanelSequenceSimilarity,
            RegulatoryFragmentEvidenceDimensionKind::RepeatsAndLowComplexity,
            RegulatoryFragmentEvidenceDimensionKind::PairSpecificJunctionUniqueness,
            RegulatoryFragmentEvidenceDimensionKind::RestrictionAndCloningRisk,
        ] {
            let lane = plan
                .evidence_dimensions
                .iter()
                .find(|lane| lane.kind == evaluated)
                .expect("evaluated evidence lane");
            assert_eq!(lane.state, RegulatoryFragmentEvidenceState::Evaluated);
            assert!(lane.assessment_sha256.starts_with("sha256:"));
        }
        for unevaluated in [
            RegulatoryFragmentEvidenceDimensionKind::EnsemblRegulatoryOverlap,
            RegulatoryFragmentEvidenceDimensionKind::TfbsModelScoreContext,
            RegulatoryFragmentEvidenceDimensionKind::CutrunAndChromatinContext,
        ] {
            let lane = plan
                .evidence_dimensions
                .iter()
                .find(|lane| lane.kind == unevaluated)
                .expect("unevaluated evidence lane");
            assert_eq!(lane.state, RegulatoryFragmentEvidenceState::NotEvaluated);
            assert!(lane.detail.contains("never a pass"));
            assert!(lane.observations.is_empty());
        }
        let kinds = plan
            .evidence_dimensions
            .iter()
            .map(|lane| lane.kind)
            .collect::<BTreeSet<_>>();
        assert_eq!(kinds.len(), 8);
        for expected in [
            RegulatoryFragmentEvidenceDimensionKind::PanelSequenceSimilarity,
            RegulatoryFragmentEvidenceDimensionKind::RepeatsAndLowComplexity,
            RegulatoryFragmentEvidenceDimensionKind::PairSpecificJunctionUniqueness,
        ] {
            assert_eq!(
                plan.evidence_dimensions
                    .iter()
                    .find(|lane| lane.kind == expected)
                    .expect("evidence lane")
                    .bindings
                    .len(),
                1
            );
        }
    }

    #[test]
    fn regulatory_fragment_similarity_repeats_and_junction_risks_stay_separate() {
        let mut fixture = planner_fixture(true, gp::GenomicRegionStrand::Plus, false);
        fixture
            .request
            .questions
            .push(RegulatoryFragmentQuestion::SpacingDependence);
        fixture.request.policy.sequence_word_size_bp = 3;
        fixture
            .request
            .requested_variants
            .push(RegulatoryFragmentGeometryRequest {
                variant_id: "long_a_spacer".to_string(),
                declared_order: 4,
                kind: RegulatoryFragmentGeometryKind::ControlledSpacing,
                instances: vec![
                    RegulatoryFragmentInstanceRequest {
                        fragment_id: "candidate_a".to_string(),
                        ..RegulatoryFragmentInstanceRequest::default()
                    },
                    RegulatoryFragmentInstanceRequest {
                        fragment_id: "partner_b".to_string(),
                        spacer_before: "A".repeat(40),
                        ..RegulatoryFragmentInstanceRequest::default()
                    },
                    RegulatoryFragmentInstanceRequest {
                        fragment_id: "minimal_promoter".to_string(),
                        ..RegulatoryFragmentInstanceRequest::default()
                    },
                ],
                ..RegulatoryFragmentGeometryRequest::default()
            });
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("plan with independent sequence evidence");
        let lane = |kind| {
            plan.evidence_dimensions
                .iter()
                .find(|lane| lane.kind == kind)
                .expect("evidence lane")
        };
        let similarity = lane(RegulatoryFragmentEvidenceDimensionKind::PanelSequenceSimilarity);
        let repeats = lane(RegulatoryFragmentEvidenceDimensionKind::RepeatsAndLowComplexity);
        let junctions =
            lane(RegulatoryFragmentEvidenceDimensionKind::PairSpecificJunctionUniqueness);
        let cloning = lane(RegulatoryFragmentEvidenceDimensionKind::RestrictionAndCloningRisk);

        assert!(similarity.observations.iter().any(|row| matches!(
            row,
            RegulatoryFragmentEvidenceObservation::PanelSequenceSimilarity {
                comparison_scope,
                ..
            } if comparison_scope == "fragment_to_reporter_vector"
        )));
        assert!(
            similarity
                .warnings
                .iter()
                .any(|warning| warning.code == "fragment_reporter_vector_word_similarity")
        );
        assert!(repeats.observations.iter().any(|row| matches!(
            row,
            RegulatoryFragmentEvidenceObservation::RepeatOrLowComplexity { evidence, .. }
                if evidence.context_tags.iter().any(|tag| tag == "low_complexity")
        )));
        assert!(
            repeats
                .warnings
                .iter()
                .all(|warning| warning.code == "repeat_or_low_complexity_context_detected")
        );
        assert!(junctions.observations.iter().any(|row| matches!(
            row,
            RegulatoryFragmentEvidenceObservation::PairSpecificJunctionUniqueness {
                panel_other_member_match_count,
                exact_unique_in_assessed_context: false,
                ..
            } if *panel_other_member_match_count > 0
        )));
        assert!(
            junctions
                .warnings
                .iter()
                .any(|warning| warning.code == "junction_not_unique_in_assessed_context")
        );
        assert!(cloning.observations.iter().all(|row| matches!(
            row,
            RegulatoryFragmentEvidenceObservation::RestrictionAndCloningRisk { .. }
        )));
        assert!(!similarity.warnings.iter().any(|warning| {
            warning.code == "repeat_or_low_complexity_context_detected"
                || warning.code == "junction_not_unique_in_assessed_context"
        }));
    }

    #[test]
    fn regulatory_fragment_mixed_assemblies_and_releases_fail_closed() {
        let mut assembly_fixture = planner_fixture(true, gp::GenomicRegionStrand::Plus, false);
        mutate_bound_region(&mut assembly_fixture, "partner_b", |region| {
            region.interval.reference.assembly_name = "T2T-CHM13v2.0".to_string();
            region.interval.reference.assembly_accession = Some("GCA_009914755.4".to_string());
            region
                .local_projection
                .as_mut()
                .expect("partner projection")
                .source_genome_id = "T2T-CHM13v2.0".to_string();
        });
        assembly_fixture
            .engine
            .state_mut()
            .metadata
            .get_mut(PROVENANCE_METADATA_KEY)
            .expect("provenance metadata")["genome_extractions"]
            .as_array_mut()
            .expect("genome extractions")
            .iter_mut()
            .find(|row| row["seq_id"] == "partner_source")
            .expect("partner provenance")["genome_id"] =
            serde_json::Value::String("T2T-CHM13v2.0".to_string());
        let error = assembly_fixture
            .engine
            .plan_regulatory_fragment_panel(assembly_fixture.request)
            .expect_err("mixed assemblies must fail");
        assert!(
            error.message.contains("mixed_reference_assemblies"),
            "unexpected mixed-assembly error: {}",
            error.message
        );

        let mut release_fixture = planner_fixture(true, gp::GenomicRegionStrand::Plus, false);
        release_fixture
            .request
            .fragments
            .iter_mut()
            .find(|fragment| fragment.fragment_id == "partner_b")
            .expect("partner")
            .reference_release = "Ensembl 115".to_string();
        let error = release_fixture
            .engine
            .plan_regulatory_fragment_panel(release_fixture.request)
            .expect_err("mixed releases must fail");
        assert!(error.message.contains("mixed_reference_releases"));
    }

    #[test]
    fn regulatory_fragment_each_non_current_projection_status_fails_closed() {
        for status in [
            gp::GenomicRegionLocalProjectionStatus::SequenceUnavailable,
            gp::GenomicRegionLocalProjectionStatus::SequenceDigestMismatch,
            gp::GenomicRegionLocalProjectionStatus::AnchorMismatch,
            gp::GenomicRegionLocalProjectionStatus::OutsideSequence,
        ] {
            let mut fixture = planner_fixture(false, gp::GenomicRegionStrand::Plus, false);
            mutate_bound_region(&mut fixture, "candidate_a", |region| {
                region.local_projection.as_mut().expect("projection").status = status;
            });
            let error = fixture
                .engine
                .plan_regulatory_fragment_panel(fixture.request)
                .expect_err("non-current projection must fail");
            assert!(error.message.contains("region_projection_not_current"));
            assert!(error.message.contains(
                GentleEngine::regulatory_fragment_projection_status_token(status)
            ));
        }
    }

    #[test]
    fn regulatory_fragment_missing_source_sequence_fails_before_planning() {
        let mut fixture = planner_fixture(false, gp::GenomicRegionStrand::Plus, false);
        fixture
            .engine
            .state_mut()
            .sequences
            .remove("candidate_source");
        let error = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect_err("missing source must fail");
        assert!(error.message.contains("sequence_unavailable"));
        assert!(!error.message.contains("cloning"));
    }

    #[test]
    fn regulatory_fragment_normalizes_input_order_and_is_byte_stable() {
        let fixture = planner_fixture(true, gp::GenomicRegionStrand::Plus, false);
        let mut reordered = fixture.request.clone();
        reordered.fragments.reverse();
        reordered.questions.reverse();
        let first = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request.clone())
            .expect("first plan");
        let second = fixture
            .engine
            .plan_regulatory_fragment_panel(reordered)
            .expect("reordered plan");
        let third = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("repeated plan");
        assert_eq!(
            serde_json::to_vec(&first).expect("serialize first"),
            serde_json::to_vec(&second).expect("serialize second")
        );
        assert_eq!(
            serde_json::to_vec(&first).expect("serialize first"),
            serde_json::to_vec(&third).expect("serialize third")
        );
    }

    #[test]
    fn regulatory_fragment_reordered_plan_members_cannot_pass_old_approval() {
        let fixture = planner_fixture(true, gp::GenomicRegionStrand::Plus, false);
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("plan");
        fixture
            .engine
            .validate_regulatory_fragment_panel_approval(&plan, &plan.proposal_digest)
            .expect("exact plan approval");
        let mut reordered = plan.clone();
        reordered.members.reverse();
        let error = fixture
            .engine
            .validate_regulatory_fragment_panel_approval(&reordered, &plan.proposal_digest)
            .expect_err("reordered proposal must fail");
        assert!(error.message.contains("proposal_content_digest_mismatch"));
    }

    #[test]
    fn regulatory_fragment_planning_does_not_mutate_state_or_write_files() {
        let fixture = planner_fixture(false, gp::GenomicRegionStrand::Plus, false);
        let state_before = serde_json::to_vec(fixture.engine.state()).expect("state before");
        let mut files_before = fs::read_dir(fixture._temp.path())
            .expect("test directory")
            .map(|entry| entry.expect("directory entry").file_name())
            .collect::<Vec<_>>();
        files_before.sort();
        fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("read-only plan");
        let state_after = serde_json::to_vec(fixture.engine.state()).expect("state after");
        let mut files_after = fs::read_dir(fixture._temp.path())
            .expect("test directory")
            .map(|entry| entry.expect("directory entry").file_name())
            .collect::<Vec<_>>();
        files_after.sort();
        assert_eq!(state_before, state_after);
        assert_eq!(files_before, files_after);
    }

    #[test]
    fn regulatory_fragment_plus_and_minus_rois_preserve_intended_insert_orientation() {
        let plus = planner_fixture(false, gp::GenomicRegionStrand::Plus, false);
        let plus_plan = plus
            .engine
            .plan_regulatory_fragment_panel(plus.request)
            .expect("plus plan");
        let minus = planner_fixture(false, gp::GenomicRegionStrand::Minus, false);
        let minus_plan = minus
            .engine
            .plan_regulatory_fragment_panel(minus.request)
            .expect("minus plan");
        let candidate_insert = |plan: &RegulatoryFragmentPanelPlan| {
            plan.members
                .iter()
                .find(|member| {
                    member.construct_kind == RegulatoryFragmentConstructKind::CandidateAlone
                })
                .expect("candidate construct")
                .insert_sequence_5prime_to_3prime
                .clone()
        };
        assert!(candidate_insert(&plus_plan).starts_with(CANDIDATE_SEQUENCE));
        assert!(
            candidate_insert(&minus_plan)
                .starts_with(&GentleEngine::reverse_complement(CANDIDATE_SEQUENCE))
        );
        assert_eq!(
            minus_plan
                .members
                .iter()
                .find(|member| {
                    member.construct_kind == RegulatoryFragmentConstructKind::CandidateAlone
                })
                .expect("minus candidate")
                .instances[0]
                .source_strand,
            gp::GenomicRegionStrand::Minus
        );
    }

    #[test]
    fn regulatory_fragment_member_bound_names_uncovered_questions() {
        let mut fixture = planner_fixture(false, gp::GenomicRegionStrand::Plus, false);
        fixture.request.max_panel_members = 2;
        let plan = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect("bounded partial plan");
        assert_eq!(plan.members.len(), 2);
        assert_eq!(
            plan.planning_label,
            RegulatoryFragmentPlanningLabel::UnresolvedExperimentalComparisonRequired
        );
        assert_eq!(
            plan.uncovered_questions,
            vec![RegulatoryFragmentQuestion::StandaloneCandidate]
        );
        assert!(plan.blockers.iter().any(|blocker| {
            blocker.code == "panel_member_bound_leaves_questions_uncovered"
                && blocker.subject_ids == ["standalone_candidate"]
        }));
    }

    #[test]
    fn regulatory_fragment_motif_control_requires_existing_explicit_policy() {
        let mut fixture = planner_fixture(false, gp::GenomicRegionStrand::Plus, false);
        fixture
            .request
            .questions
            .push(RegulatoryFragmentQuestion::MotifDisruption);
        fixture.request.mutation_policy = PromoterReporterPanelMutationPolicy::Unspecified;
        fixture.request.motif_disruption = Some(RegulatoryFragmentMotifDisruptionRequest {
            control_id: "candidate_motif_edit".to_string(),
            fragment_id: "candidate_a".to_string(),
            motif_start_in_fragment_0based: 4,
            motif_end_in_fragment_0based_exclusive: 12,
            motif_forward_strand: true,
        });
        let error = fixture
            .engine
            .plan_regulatory_fragment_panel(fixture.request)
            .expect_err("unspecified mutation policy must fail");
        assert!(error.message.contains("mutation_policy_unspecified"));
    }

    #[test]
    fn regulatory_fragment_default_policy_round_trips_to_explicit_bytes() {
        let implicit: RegulatoryFragmentPanelRequest =
            serde_json::from_value(serde_json::json!({})).expect("implicit defaults");
        let explicit: RegulatoryFragmentPanelRequest = serde_json::from_value(serde_json::json!({
            "schema": REGULATORY_FRAGMENT_PANEL_REQUEST_SCHEMA,
            "insertion_context_feature_id": "multiple_cloning_region",
            "mutation_policy": "unspecified",
            "max_panel_members": 8,
            "max_construct_length_bp": 5000,
            "policy": {
                "include_promoterless_control": true,
                "include_minimal_promoter_control": true,
                "include_reference_control_when_bound": true,
                "high_similarity_identity_fraction": 0.95,
                "high_similarity_coverage_fraction": 0.90,
                "vector_context_flank_bp": 24,
                "max_candidate_constructs": 20,
                "sequence_word_size_bp": 12,
                "near_exact_max_mismatches": 1,
                "junction_flank_bp": 16,
                "max_evidence_observations_per_dimension": 512
            }
        }))
        .expect("explicit defaults");
        assert_eq!(
            serde_json::to_vec(&implicit).expect("implicit bytes"),
            serde_json::to_vec(&explicit).expect("explicit bytes")
        );
    }

    #[test]
    fn legacy_promoter_reporter_proposal_serialization_snapshot_is_unchanged() {
        let proposal = PromoterReporterPanelProposal {
            schema: PROMOTER_REPORTER_PANEL_PROPOSAL_SCHEMA.to_string(),
            proposal_id: "legacy_vkorc1_serpine1".to_string(),
            proposal_digest: "sha256:legacy-proposal".to_string(),
            request_sha256: "sha256:legacy-request".to_string(),
            baseline_state_sha256: "sha256:legacy-state".to_string(),
            approval_required: true,
            nonclaims: vec!["legacy serialization sentinel".to_string()],
            ..PromoterReporterPanelProposal::default()
        };
        let bytes = serde_json::to_vec(&proposal).expect("legacy proposal JSON");
        assert_eq!(
            sha256_prefixed_bytes(&bytes),
            "sha256:fbc8540ce603c18d612d04d534450be964adee743e775af78537dd93868bac94"
        );
    }
}
