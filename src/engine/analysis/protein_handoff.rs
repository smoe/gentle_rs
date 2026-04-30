//! Protein-to-DNA construct-reasoning helpers.
//!
//! This module keeps the first protein-to-DNA handoff assistant engine-owned
//! and construct-graph-backed instead of duplicating selection logic in the
//! GUI specialist.

use super::*;
use serde_json::json;
use std::cmp::Ordering;

#[derive(Debug, Clone)]
struct ProteinHandoffDraftCandidate {
    candidate: ConstructCandidate,
    evidence_rows: Vec<DesignEvidence>,
    provenance_score: f64,
    coverage_fraction: f64,
    warning_count: usize,
    strategy_rank: usize,
    rationale: String,
}

#[derive(Debug, Clone)]
struct TranscriptNativeSource {
    report_id: String,
    transcript_id: String,
    transcript_feature_id: usize,
    derivation: TranscriptProteinDerivation,
}

impl GentleEngine {
    pub(crate) fn construct_reasoning_graph_has_protein_to_dna_handoff(
        graph: &ConstructReasoningGraph,
    ) -> bool {
        graph
            .candidates
            .iter()
            .any(|candidate| candidate.protein_to_dna_handoff.is_some())
    }

    pub(crate) fn construct_reasoning_graph_protein_to_dna_handoff_summary(
        graph: &ConstructReasoningGraph,
    ) -> (usize, Vec<String>) {
        let mut source_protein_seq_ids = graph
            .candidates
            .iter()
            .filter_map(|candidate| candidate.protein_to_dna_handoff.as_ref())
            .map(|detail| detail.source_protein_seq_id.clone())
            .collect::<Vec<_>>();
        source_protein_seq_ids.sort();
        source_protein_seq_ids.dedup();
        (
            graph
                .candidates
                .iter()
                .filter(|candidate| candidate.protein_to_dna_handoff.is_some())
                .count(),
            source_protein_seq_ids,
        )
    }

    pub(crate) fn build_protein_to_dna_handoff_reasoning_graph(
        &mut self,
        seq_id: &str,
        protein_seq_id: &str,
        transcript_filter: Option<&str>,
        projection_id: Option<&str>,
        ensembl_entry_id: Option<&str>,
        feature_query: Option<&str>,
        ranking_goal: ProteinToDnaHandoffRankingGoal,
        speed_profile: Option<TranslationSpeedProfile>,
        speed_mark: Option<TranslationSpeedMark>,
        translation_table: Option<usize>,
        target_anneal_tm_c: Option<f64>,
        anneal_window_bp: Option<usize>,
        objective_id: Option<&str>,
        graph_id: Option<&str>,
        op_id: Option<&str>,
        run_id: Option<&str>,
    ) -> Result<ConstructReasoningGraph, EngineError> {
        let seq_id = seq_id.trim();
        let protein_seq_id = protein_seq_id.trim();
        if seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "BuildProteinToDnaHandoffReasoning requires a non-empty seq_id"
                    .to_string(),
            });
        }
        if protein_seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "BuildProteinToDnaHandoffReasoning requires a non-empty protein_seq_id"
                    .to_string(),
            });
        }

        let protein = self
            .state
            .sequences
            .get(protein_seq_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Protein sequence '{}' not found", protein_seq_id),
            })?;
        if !protein.is_protein_sequence() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "BuildProteinToDnaHandoffReasoning requires a protein sequence; '{}' has molecule_type={:?}",
                    protein_seq_id,
                    protein.molecule_type()
                ),
            });
        }
        if !self.state.sequences.contains_key(seq_id) {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            });
        }

        let transcript_filter = transcript_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let projection_id = projection_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let ensembl_entry_id = ensembl_entry_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let feature_query = feature_query
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());

        let protein_sequence = protein.get_forward_string().to_ascii_uppercase();
        let requested_amino_acids = protein_sequence.len();
        let mut graph = self.assemble_construct_reasoning_graph(seq_id, objective_id, graph_id)?;
        graph.op_id = op_id.map(|value| value.to_string());
        graph.run_id = run_id.map(|value| value.to_string());
        if graph.objective.goal.trim().is_empty() {
            graph.objective.goal = format!(
                "Protein-to-DNA handoff for '{}' using '{}'",
                seq_id, protein_seq_id
            );
        }
        if graph.objective.title.trim().is_empty() {
            graph.objective.title = format!(
                "Protein-to-DNA handoff: {}",
                protein.name().as_deref().unwrap_or(protein_seq_id)
            );
        }
        if !graph
            .objective
            .required_roles
            .iter()
            .any(|role| *role == ConstructRole::Cds)
        {
            graph.objective.required_roles.push(ConstructRole::Cds);
        }

        let protein_source_evidence_id = format!(
            "protein_handoff_source_{}",
            Self::normalize_id_token(protein_seq_id)
        );
        graph.evidence.push(DesignEvidence {
            evidence_id: protein_source_evidence_id.clone(),
            seq_id: protein_seq_id.to_string(),
            scope: EvidenceScope::WholeConstruct,
            start_0based: 0,
            end_0based_exclusive: requested_amino_acids,
            role: ConstructRole::Other,
            evidence_class: EvidenceClass::ReliableAnnotation,
            label: format!("Protein source '{}'", protein_seq_id),
            rationale: "Selected protein sequence for ranked DNA handoff reasoning".to_string(),
            confidence: Some(1.0),
            context_tags: vec![
                "protein_to_dna_handoff".to_string(),
                format!("ranking_goal={}", ranking_goal.as_str()),
            ],
            provenance_kind: "protein_seq".to_string(),
            provenance_refs: vec![protein_seq_id.to_string()],
            editable_status: EditableStatus::Draft,
            ..DesignEvidence::default()
        });

        let mut facts = vec![Self::construct_reasoning_build_fact(
            "protein_to_dna_handoff_context",
            "protein_to_dna_handoff_context",
            "Protein-to-DNA handoff context".to_string(),
            format!(
                "Selected protein '{}' for DNA handoff on '{}' with ranking goal '{}'.",
                protein_seq_id,
                seq_id,
                ranking_goal.as_str()
            ),
            vec![protein_source_evidence_id.clone()],
            json!({
                "seq_id": seq_id,
                "protein_seq_id": protein_seq_id,
                "transcript_filter": transcript_filter.clone(),
                "projection_id": projection_id.clone(),
                "ensembl_entry_id": ensembl_entry_id.clone(),
                "feature_query": feature_query.clone(),
                "ranking_goal": ranking_goal.as_str(),
                "requested_speed_profile": speed_profile.map(TranslationSpeedProfile::as_str),
                "requested_speed_mark": speed_mark.map(TranslationSpeedMark::as_str),
                "requested_translation_table": translation_table,
                "target_anneal_tm_c": target_anneal_tm_c,
                "anneal_window_bp": anneal_window_bp,
            }),
        )];
        let mut decisions = vec![];
        let mut candidates: Vec<ProteinHandoffDraftCandidate> = vec![];

        let transcript_native = self.build_transcript_native_handoff_candidate(
            seq_id,
            protein_seq_id,
            transcript_filter.as_deref(),
            requested_amino_acids,
            ranking_goal,
        )?;
        Self::push_protein_handoff_family_fact_and_decision(
            &mut facts,
            &mut decisions,
            "transcript_native_reuse",
            "Transcript-native CDS reuse",
            transcript_native.as_ref(),
            vec![protein_source_evidence_id.clone()],
        );
        if let Some(candidate) = transcript_native {
            graph
                .evidence
                .extend(candidate.evidence_rows.iter().cloned());
            candidates.push(candidate);
        }

        let feature_coding = self.build_feature_coding_handoff_candidate(
            seq_id,
            protein_seq_id,
            projection_id.as_deref(),
            transcript_filter.as_deref(),
            feature_query.as_deref(),
            requested_amino_acids,
            ranking_goal,
            speed_profile,
            speed_mark,
        )?;
        Self::push_protein_handoff_family_fact_and_decision(
            &mut facts,
            &mut decisions,
            "feature_coding_dna",
            "Feature-coding DNA handoff",
            feature_coding.as_ref(),
            vec![protein_source_evidence_id.clone()],
        );
        if let Some(candidate) = feature_coding {
            graph
                .evidence
                .extend(candidate.evidence_rows.iter().cloned());
            candidates.push(candidate);
        }

        let reverse_translated = self.build_reverse_translated_handoff_candidate(
            seq_id,
            protein_seq_id,
            requested_amino_acids,
            ranking_goal,
            speed_profile,
            speed_mark,
            translation_table,
            target_anneal_tm_c,
            anneal_window_bp,
            ensembl_entry_id.as_deref(),
        )?;
        Self::push_protein_handoff_family_fact_and_decision(
            &mut facts,
            &mut decisions,
            "reverse_translated_synthetic",
            "Reverse-translated synthetic fallback",
            reverse_translated.as_ref(),
            vec![protein_source_evidence_id.clone()],
        );
        if let Some(candidate) = reverse_translated {
            graph
                .evidence
                .extend(candidate.evidence_rows.iter().cloned());
            candidates.push(candidate);
        }

        if candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Could not derive any protein-to-DNA handoff candidates for '{}' on '{}'",
                    protein_seq_id, seq_id
                ),
            });
        }

        candidates.sort_by(|left, right| {
            Self::protein_to_dna_candidate_ordering(ranking_goal, left, right)
        });
        let ranked_candidate_ids = candidates
            .iter()
            .map(|row| row.candidate.candidate_id.clone())
            .collect::<Vec<_>>();
        let ranking_rationale = candidates
            .iter()
            .enumerate()
            .map(|(idx, row)| format!("{}. {}: {}", idx + 1, row.candidate.title, row.rationale))
            .collect::<Vec<_>>();
        decisions.push(Self::construct_reasoning_build_decision(
            "protein_to_dna_handoff_rank",
            "protein_to_dna_handoff_rank",
            "Protein-to-DNA candidate ranking".to_string(),
            format!(
                "Ranked {} handoff candidates using '{}' preference rules.",
                candidates.len(),
                ranking_goal.as_str()
            ),
            graph
                .evidence
                .iter()
                .filter(|row| {
                    row.context_tags
                        .iter()
                        .any(|tag| tag == "protein_to_dna_handoff")
                })
                .map(|row| row.evidence_id.clone())
                .collect(),
            vec!["protein_to_dna_handoff_context".to_string()],
            json!({
                "ranking_goal": ranking_goal.as_str(),
                "candidate_ids": ranked_candidate_ids,
                "rationale_rows": ranking_rationale,
            }),
        ));

        graph.facts.extend(facts);
        graph.decisions.extend(decisions);
        graph.candidates = candidates.into_iter().map(|row| row.candidate).collect();
        graph.notes.push(
            "v1 protein-to-DNA handoff reasoning ranks transcript-native CDS reuse, mapped feature-coding DNA, and reverse-translated synthetic fallbacks without materializing DNA sequences."
                .to_string(),
        );
        self.upsert_construct_reasoning_graph(graph.clone())?;
        Ok(graph)
    }

    fn push_protein_handoff_family_fact_and_decision(
        facts: &mut Vec<DesignFact>,
        decisions: &mut Vec<DesignDecisionNode>,
        family_id: &str,
        family_title: &str,
        draft: Option<&ProteinHandoffDraftCandidate>,
        input_evidence_ids: Vec<String>,
    ) {
        let (status, rationale, candidate_ids, warnings, parameters_json) =
            if let Some(draft) = draft {
                (
                    "admitted",
                    draft.rationale.clone(),
                    vec![draft.candidate.candidate_id.clone()],
                    draft.candidate.warnings.clone(),
                    json!({
                        "status": "admitted",
                        "strategy": draft
                            .candidate
                            .protein_to_dna_handoff
                            .as_ref()
                            .map(|detail| detail.strategy.as_str()),
                        "candidate_id": draft.candidate.candidate_id.clone(),
                        "warnings": warnings_to_json(&draft.candidate.warnings),
                    }),
                )
            } else {
                (
                    "rejected",
                    format!(
                        "{family_title} could not be admitted from the currently selected evidence."
                    ),
                    vec![],
                    Vec::<String>::new(),
                    json!({
                        "status": "rejected",
                    }),
                )
            };
        let fact_id = format!("protein_to_dna_handoff_family_{family_id}");
        facts.push(Self::construct_reasoning_build_fact(
            &fact_id,
            "protein_to_dna_handoff_family",
            format!("{family_title} ({status})"),
            rationale.clone(),
            input_evidence_ids.clone(),
            json!({
                "family_id": family_id,
                "status": status,
                "candidate_ids": candidate_ids,
                "warnings": warnings_to_json(&warnings),
            }),
        ));
        decisions.push(Self::construct_reasoning_build_decision(
            &format!("protein_to_dna_handoff_decision_{family_id}"),
            "protein_to_dna_handoff_family_decision",
            family_title.to_string(),
            rationale,
            input_evidence_ids,
            vec![fact_id],
            parameters_json,
        ));
    }

    fn build_transcript_native_handoff_candidate(
        &self,
        seq_id: &str,
        protein_seq_id: &str,
        transcript_filter: Option<&str>,
        requested_amino_acids: usize,
        ranking_goal: ProteinToDnaHandoffRankingGoal,
    ) -> Result<Option<ProteinHandoffDraftCandidate>, EngineError> {
        let Some(source) =
            self.find_transcript_native_protein_source(seq_id, protein_seq_id, transcript_filter)?
        else {
            return Ok(None);
        };
        let source_seq = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        let source_sequence_upper = source_seq
            .get_forward_string()
            .to_ascii_uppercase()
            .into_bytes();
        let source_features = source_seq.features().to_vec();
        let transcript_feature = source_features
            .get(source.transcript_feature_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Transcript feature '{}' was not found in sequence '{}'",
                    source.transcript_feature_id + 1,
                    seq_id
                ),
            })?;
        let (derived_transcript, _, _, _, _, _) = Self::derive_transcript_sequence_from_feature(
            &source_sequence_upper,
            transcript_feature,
            &source_features,
            source.transcript_feature_id,
            seq_id,
        )?;
        let coding_sequence = Self::extract_ranges_1based_from_sequence(
            &derived_transcript.get_forward_string(),
            &source.derivation.cds_ranges_1based,
        );
        if coding_sequence.is_empty() {
            return Ok(None);
        }
        let preserved_constraints = vec![
            format!("reuse transcript '{}' CDS", source.transcript_id),
            "keep native coding provenance".to_string(),
            "avoid synthetic codon rewriting".to_string(),
        ];
        let detail = ProteinToDnaHandoffCandidate {
            strategy: ProteinToDnaHandoffStrategy::TranscriptNativeReuse,
            ranking_goal,
            source_protein_seq_id: protein_seq_id.to_string(),
            source_artifact_refs: vec![format!("protein_derivation_report:{}", source.report_id)],
            transcript_id: Some(source.transcript_id.clone()),
            projection_id: None,
            ensembl_entry_id: None,
            feature_query: None,
            matched_feature_label: None,
            coverage: ProteinToDnaHandoffCoverage {
                amino_acid_start_1based: Some(1),
                amino_acid_end_1based: Some(source.derivation.protein_length_aa),
                covered_amino_acids: source.derivation.protein_length_aa,
                requested_amino_acids,
                preserved_feature_labels: vec!["full_transcript_cds".to_string()],
                relaxed_feature_labels: vec![],
            },
            preserved_constraints: preserved_constraints.clone(),
            relaxed_constraints: vec![],
            translation_table: Some(source.derivation.translation_table),
            speed_profile: source
                .derivation
                .translation_speed_profile_hint
                .as_deref()
                .and_then(|raw| Self::infer_translation_speed_profile_enum(Some(raw))),
            speed_profile_source: source.derivation.translation_speed_profile_source,
            speed_mark: None,
            provenance_score: Some(1.0),
            codon_policy_summary: "Reuse existing CDS as encoded".to_string(),
            next_step_recommendations: vec![
                format!(
                    "Use transcript '{}' CDS as the primary cloning handoff.",
                    source.transcript_id
                ),
                "Treat synthetic reverse translation as fallback only if transcript-native CDS cannot be reused.".to_string(),
            ],
            future_materialization_seq_id: None,
        };
        let evidence_id = format!(
            "protein_handoff_native_{}",
            Self::normalize_id_token(protein_seq_id)
        );
        let candidate = ConstructCandidate {
            candidate_id: format!(
                "protein_handoff_native_{}",
                Self::normalize_id_token(protein_seq_id)
            ),
            objective_id: format!("protein_handoff_{}", Self::normalize_id_token(seq_id)),
            title: "Transcript-native CDS reuse".to_string(),
            component_ids: vec![seq_id.to_string(), protein_seq_id.to_string()],
            compactness_score: Some(1.0),
            confidence_score: Some(1.0),
            cloning_complexity_score: Some(0.15),
            host_fit_score: Some(0.8),
            protein_to_dna_handoff: Some(detail),
            notes: vec![
                format!("Transcript '{}'", source.transcript_id),
                format!("Coding length {} bp", coding_sequence.len()),
            ],
            ..ConstructCandidate::default()
        };
        Ok(Some(ProteinHandoffDraftCandidate {
            evidence_rows: vec![DesignEvidence {
                evidence_id,
                seq_id: seq_id.to_string(),
                scope: EvidenceScope::WholeConstruct,
                start_0based: 0,
                end_0based_exclusive: source_seq.len(),
                role: ConstructRole::Cds,
                evidence_class: EvidenceClass::HardFact,
                label: format!("Transcript-native CDS '{}'", source.transcript_id),
                rationale: format!(
                    "Protein '{}' was derived from transcript '{}' in '{}', so the existing CDS can be handed off directly.",
                    protein_seq_id, source.transcript_id, seq_id
                ),
                confidence: Some(1.0),
                context_tags: vec![
                    "protein_to_dna_handoff".to_string(),
                    "transcript_native_reuse".to_string(),
                ],
                provenance_kind: "protein_derivation_report".to_string(),
                provenance_refs: vec![source.report_id],
                editable_status: EditableStatus::Draft,
                ..DesignEvidence::default()
            }],
            provenance_score: 1.0,
            coverage_fraction: 1.0,
            warning_count: 0,
            strategy_rank: 0,
            rationale: preserved_constraints.join("; "),
            candidate,
        }))
    }

    fn build_feature_coding_handoff_candidate(
        &self,
        seq_id: &str,
        protein_seq_id: &str,
        projection_id: Option<&str>,
        transcript_filter: Option<&str>,
        feature_query: Option<&str>,
        requested_amino_acids: usize,
        ranking_goal: ProteinToDnaHandoffRankingGoal,
        speed_profile: Option<TranslationSpeedProfile>,
        speed_mark: Option<TranslationSpeedMark>,
    ) -> Result<Option<ProteinHandoffDraftCandidate>, EngineError> {
        let Some(projection_id) = projection_id else {
            return Ok(None);
        };
        let Some(feature_query) = feature_query else {
            return Ok(None);
        };
        let report = self.query_uniprot_feature_coding_dna(
            projection_id,
            feature_query,
            transcript_filter,
            UniprotFeatureCodingDnaQueryMode::Both,
            speed_profile,
        )?;
        let Some(best) = report
            .matches
            .iter()
            .max_by(|left, right| {
                let left_len = left.aa_end.saturating_sub(left.aa_start).saturating_add(1);
                let right_len = right
                    .aa_end
                    .saturating_sub(right.aa_start)
                    .saturating_add(1);
                left_len
                    .cmp(&right_len)
                    .then(left.transcript_id.cmp(&right.transcript_id))
            })
            .cloned()
        else {
            return Ok(None);
        };
        let matched_feature_label = Self::protein_handoff_feature_label(&best);
        let preserved_constraints = vec![
            format!("preserve mapped feature '{}'", matched_feature_label),
            "reuse genomic coding DNA provenance".to_string(),
        ];
        let relaxed_constraints = vec![
            "full-length CDS remains a separate follow-on if the handoff only targets one mapped feature".to_string(),
        ];
        let coverage_len = best.aa_end.saturating_sub(best.aa_start).saturating_add(1);
        let warning_count = best.warnings.len().saturating_add(report.warnings.len());
        let detail = ProteinToDnaHandoffCandidate {
            strategy: ProteinToDnaHandoffStrategy::FeatureCodingDna,
            ranking_goal,
            source_protein_seq_id: protein_seq_id.to_string(),
            source_artifact_refs: vec![format!("uniprot_projection:{}", report.projection_id)],
            transcript_id: Some(best.transcript_id.clone()),
            projection_id: Some(report.projection_id.clone()),
            ensembl_entry_id: None,
            feature_query: Some(report.feature_query.clone()),
            matched_feature_label: Some(matched_feature_label.clone()),
            coverage: ProteinToDnaHandoffCoverage {
                amino_acid_start_1based: Some(best.aa_start),
                amino_acid_end_1based: Some(best.aa_end),
                covered_amino_acids: coverage_len,
                requested_amino_acids,
                preserved_feature_labels: vec![matched_feature_label.clone()],
                relaxed_feature_labels: vec![],
            },
            preserved_constraints: preserved_constraints.clone(),
            relaxed_constraints: relaxed_constraints.clone(),
            translation_table: None,
            speed_profile: report.resolved_translation_speed_profile,
            speed_profile_source: report.resolved_translation_speed_profile_source,
            speed_mark,
            provenance_score: Some(0.8),
            codon_policy_summary: if best.translation_speed_optimized_dna.is_some() {
                "Prefer mapped coding DNA and keep a translation-speed-optimized alternative in reserve".to_string()
            } else {
                "Reuse mapped genomic coding DNA as encoded".to_string()
            },
            next_step_recommendations: vec![
                format!(
                    "Use '{}' on transcript '{}' if the cloning handoff only needs this mapped protein feature.",
                    matched_feature_label, best.transcript_id
                ),
                "Keep a full-length CDS handoff candidate available separately when downstream cloning needs the complete product.".to_string(),
            ],
            future_materialization_seq_id: None,
        };
        let candidate = ConstructCandidate {
            candidate_id: format!(
                "protein_handoff_feature_{}",
                Self::normalize_id_token(&report.projection_id)
            ),
            objective_id: format!("protein_handoff_{}", Self::normalize_id_token(seq_id)),
            title: format!("Feature-coding DNA: {}", matched_feature_label),
            component_ids: vec![seq_id.to_string(), protein_seq_id.to_string()],
            compactness_score: Some(coverage_len as f64 / requested_amino_acids.max(1) as f64),
            confidence_score: Some(0.85),
            cloning_complexity_score: Some(0.3),
            host_fit_score: Some(if best.translation_speed_optimized_dna.is_some() {
                0.75
            } else {
                0.55
            }),
            warnings: best
                .warnings
                .iter()
                .chain(report.warnings.iter())
                .cloned()
                .collect(),
            notes: vec![
                format!("Transcript '{}'", best.transcript_id),
                format!("AA span {}..{}", best.aa_start, best.aa_end),
                format!("Genomic coding length {} bp", best.genomic_coding_dna.len()),
            ],
            protein_to_dna_handoff: Some(detail),
            ..ConstructCandidate::default()
        };
        Ok(Some(ProteinHandoffDraftCandidate {
            evidence_rows: vec![DesignEvidence {
                evidence_id: format!(
                    "protein_handoff_feature_{}",
                    Self::normalize_id_token(&report.projection_id)
                ),
                seq_id: protein_seq_id.to_string(),
                scope: EvidenceScope::SequenceSpan,
                start_0based: best.aa_start.saturating_sub(1),
                end_0based_exclusive: best.aa_end,
                role: ConstructRole::Other,
                evidence_class: EvidenceClass::ReliableAnnotation,
                label: format!("Mapped feature '{}'", matched_feature_label),
                rationale: format!(
                    "Stored UniProt projection '{}' links feature '{}' back to coding DNA on transcript '{}'.",
                    report.projection_id, matched_feature_label, best.transcript_id
                ),
                confidence: Some(0.85),
                context_tags: vec![
                    "protein_to_dna_handoff".to_string(),
                    "feature_coding_dna".to_string(),
                ],
                provenance_kind: "uniprot_projection".to_string(),
                provenance_refs: vec![report.projection_id],
                warnings: report.warnings.clone(),
                editable_status: EditableStatus::Draft,
                ..DesignEvidence::default()
            }],
            provenance_score: 0.8,
            coverage_fraction: coverage_len as f64 / requested_amino_acids.max(1) as f64,
            warning_count,
            strategy_rank: 1,
            rationale: format!(
                "Mapped feature '{}' preserves aa {}..{} with coding-DNA-backed provenance.",
                matched_feature_label, best.aa_start, best.aa_end
            ),
            candidate,
        }))
    }

    fn build_reverse_translated_handoff_candidate(
        &self,
        seq_id: &str,
        protein_seq_id: &str,
        requested_amino_acids: usize,
        ranking_goal: ProteinToDnaHandoffRankingGoal,
        speed_profile: Option<TranslationSpeedProfile>,
        speed_mark: Option<TranslationSpeedMark>,
        translation_table: Option<usize>,
        target_anneal_tm_c: Option<f64>,
        anneal_window_bp: Option<usize>,
        ensembl_entry_id: Option<&str>,
    ) -> Result<Option<ProteinHandoffDraftCandidate>, EngineError> {
        let protein = self
            .state
            .sequences
            .get(protein_seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Protein sequence '{}' not found", protein_seq_id),
            })?;
        let speed_profile_resolution = speed_profile
            .map(|profile| {
                Self::translation_speed_profile_resolution_from_profile(
                    profile,
                    TranslationSpeedProfileSource::ExplicitRequest,
                    None,
                )
            })
            .or_else(|| Self::sequence_feature_translation_speed_resolution(protein));
        let preferred_species_label = speed_profile_resolution
            .as_ref()
            .map(|resolution| resolution.reference_species.as_str());
        let (
            effective_translation_table,
            effective_translation_table_label,
            _effective_translation_table_source,
            _translation_context_organism,
            _translation_context_organelle,
            translation_table_warnings,
        ) = Self::resolve_translation_table_for_reverse_translation(
            protein,
            translation_table.filter(|value| *value > 0),
        );
        let anneal_window_bp = anneal_window_bp.unwrap_or(20).max(6);
        let protein_sequence = protein.get_forward_string().to_ascii_uppercase();
        let (coding_sequence, reverse_translation_warnings) =
            Self::build_reverse_translated_coding_sequence(
                &protein_sequence,
                effective_translation_table,
                preferred_species_label,
                speed_mark,
                target_anneal_tm_c,
                anneal_window_bp,
            );
        let mut warnings = translation_table_warnings;
        warnings.extend(
            speed_profile_resolution
                .as_ref()
                .map(|resolution| resolution.warnings.clone())
                .unwrap_or_default(),
        );
        warnings.extend(reverse_translation_warnings);
        let (
            preferred_synonymous_choice_count,
            alternative_synonymous_choice_count,
            fallback_unknown_codon_count,
            _gc_fraction,
            realized_anneal_tm_c,
        ) = Self::reverse_translation_choice_diagnostics(
            &protein_sequence,
            &coding_sequence,
            effective_translation_table,
            preferred_species_label,
            target_anneal_tm_c,
            anneal_window_bp,
        );
        let provenance_score = 0.35;
        let detail = ProteinToDnaHandoffCandidate {
            strategy: ProteinToDnaHandoffStrategy::ReverseTranslatedSynthetic,
            ranking_goal,
            source_protein_seq_id: protein_seq_id.to_string(),
            source_artifact_refs: ensembl_entry_id
                .map(|value| vec![format!("ensembl_protein_entry:{value}")])
                .unwrap_or_default(),
            transcript_id: None,
            projection_id: None,
            ensembl_entry_id: ensembl_entry_id.map(|value| value.to_string()),
            feature_query: None,
            matched_feature_label: None,
            coverage: ProteinToDnaHandoffCoverage {
                amino_acid_start_1based: Some(1),
                amino_acid_end_1based: Some(protein_sequence.len()),
                covered_amino_acids: protein_sequence.len(),
                requested_amino_acids,
                preserved_feature_labels: vec!["full_protein_sequence".to_string()],
                relaxed_feature_labels: if fallback_unknown_codon_count > 0 {
                    vec!["explicit_codon_choice_missing".to_string()]
                } else {
                    vec![]
                },
            },
            preserved_constraints: vec!["retain full amino-acid sequence".to_string()],
            relaxed_constraints: vec![
                "native nucleotide provenance is replaced by synthetic codon choice".to_string(),
            ],
            translation_table: Some(effective_translation_table),
            speed_profile: speed_profile_resolution.as_ref().map(|row| row.profile),
            speed_profile_source: speed_profile_resolution.as_ref().map(|row| row.source),
            speed_mark,
            provenance_score: Some(provenance_score),
            codon_policy_summary: format!(
                "Synthetic reverse translation using table {} ('{}'), preferred={}, alternative={}, fallback={}, tm={}",
                effective_translation_table,
                effective_translation_table_label,
                preferred_synonymous_choice_count,
                alternative_synonymous_choice_count,
                fallback_unknown_codon_count,
                realized_anneal_tm_c
                    .map(|value| format!("{value:.1}C"))
                    .unwrap_or_else(|| "-".to_string())
            ),
            next_step_recommendations: vec![
                "Use this candidate when no native or mapped coding DNA is available, or when a portable synthetic CDS is explicitly desired.".to_string(),
                "Treat cloning and primer planning as a separate follow-on once one candidate is explicitly materialized.".to_string(),
            ],
            future_materialization_seq_id: None,
        };
        let candidate = ConstructCandidate {
            candidate_id: format!(
                "protein_handoff_reverse_{}",
                Self::normalize_id_token(protein_seq_id)
            ),
            objective_id: format!("protein_handoff_{}", Self::normalize_id_token(seq_id)),
            title: "Reverse-translated synthetic CDS".to_string(),
            component_ids: vec![protein_seq_id.to_string()],
            compactness_score: Some(1.0),
            confidence_score: Some(0.7),
            cloning_complexity_score: Some(0.55),
            host_fit_score: Some(if speed_profile_resolution.is_some() {
                0.7
            } else {
                0.45
            }),
            warnings: warnings.clone(),
            notes: vec![
                format!("Synthetic coding length {} bp", coding_sequence.len()),
                format!("Translation table {}", effective_translation_table),
            ],
            protein_to_dna_handoff: Some(detail),
            ..ConstructCandidate::default()
        };
        Ok(Some(ProteinHandoffDraftCandidate {
            evidence_rows: vec![DesignEvidence {
                evidence_id: format!(
                    "protein_handoff_reverse_{}",
                    Self::normalize_id_token(protein_seq_id)
                ),
                seq_id: protein_seq_id.to_string(),
                scope: EvidenceScope::WholeConstruct,
                start_0based: 0,
                end_0based_exclusive: protein_sequence.len(),
                role: ConstructRole::Cds,
                evidence_class: EvidenceClass::SoftHypothesis,
                label: "Synthetic reverse-translation fallback".to_string(),
                rationale: format!(
                    "Reverse translation can provide a full-length CDS handoff for '{}' when native or mapped coding DNA is unavailable or intentionally not reused.",
                    protein_seq_id
                ),
                confidence: Some(0.7),
                context_tags: vec![
                    "protein_to_dna_handoff".to_string(),
                    "reverse_translated_synthetic".to_string(),
                ],
                provenance_kind: "reverse_translation_policy".to_string(),
                provenance_refs: vec![protein_seq_id.to_string()],
                warnings: warnings.clone(),
                editable_status: EditableStatus::Draft,
                ..DesignEvidence::default()
            }],
            provenance_score,
            coverage_fraction: 1.0,
            warning_count: warnings.len(),
            strategy_rank: 2,
            rationale: "Provides a portable full-length fallback when native or mapped DNA evidence is unavailable.".to_string(),
            candidate,
        }))
    }

    fn find_transcript_native_protein_source(
        &self,
        seq_id: &str,
        protein_seq_id: &str,
        transcript_filter: Option<&str>,
    ) -> Result<Option<TranscriptNativeSource>, EngineError> {
        let normalized_filter = transcript_filter
            .map(Self::normalize_transcript_probe)
            .filter(|value| !value.is_empty());
        let store = self.read_protein_derivation_report_store();
        for report in store.reports.values() {
            if !report.seq_id.eq_ignore_ascii_case(seq_id) {
                continue;
            }
            for row in &report.rows {
                if !row.protein_seq_id.eq_ignore_ascii_case(protein_seq_id) {
                    continue;
                }
                if let Some(filter) = normalized_filter.as_deref()
                    && Self::normalize_transcript_probe(&row.derivation.transcript_id) != filter
                {
                    continue;
                }
                return Ok(Some(TranscriptNativeSource {
                    report_id: report.report_id.clone(),
                    transcript_id: row.derivation.transcript_id.clone(),
                    transcript_feature_id: row.transcript_feature_id,
                    derivation: row.derivation.clone(),
                }));
            }
        }
        Ok(None)
    }

    fn protein_handoff_feature_label(best: &UniprotFeatureCodingDnaMatch) -> String {
        let mut parts = vec![best.feature_key.trim().to_string()];
        if let Some(note) = best
            .feature_note
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            parts.push(note.to_string());
        }
        parts.retain(|value| !value.is_empty());
        if parts.is_empty() {
            format!("aa {}..{}", best.aa_start, best.aa_end)
        } else {
            parts.join(": ")
        }
    }

    fn extract_ranges_1based_from_sequence(sequence: &str, ranges: &[(usize, usize)]) -> String {
        let bytes = sequence.as_bytes();
        let mut out = String::new();
        for (start_1based, end_1based) in ranges {
            let start = start_1based.saturating_sub(1);
            let end = (*end_1based).min(bytes.len());
            if start >= end || start >= bytes.len() {
                continue;
            }
            out.push_str(&String::from_utf8_lossy(&bytes[start..end]));
        }
        out
    }

    fn protein_to_dna_candidate_ordering(
        ranking_goal: ProteinToDnaHandoffRankingGoal,
        left: &ProteinHandoffDraftCandidate,
        right: &ProteinHandoffDraftCandidate,
    ) -> Ordering {
        let cmp_desc_f64 =
            |left: f64, right: f64| right.partial_cmp(&left).unwrap_or(Ordering::Equal);
        let cmp_asc_usize = |left: usize, right: usize| left.cmp(&right);
        match ranking_goal {
            ProteinToDnaHandoffRankingGoal::BalancedProvenance => {
                cmp_desc_f64(left.provenance_score, right.provenance_score)
                    .then_with(|| cmp_desc_f64(left.coverage_fraction, right.coverage_fraction))
                    .then_with(|| cmp_asc_usize(left.warning_count, right.warning_count))
                    .then_with(|| cmp_asc_usize(left.strategy_rank, right.strategy_rank))
            }
            ProteinToDnaHandoffRankingGoal::NativeFidelity => {
                cmp_desc_f64(left.coverage_fraction, right.coverage_fraction)
                    .then_with(|| cmp_desc_f64(left.provenance_score, right.provenance_score))
                    .then_with(|| cmp_asc_usize(left.warning_count, right.warning_count))
                    .then_with(|| cmp_asc_usize(left.strategy_rank, right.strategy_rank))
            }
            ProteinToDnaHandoffRankingGoal::ExpressionOptimized => cmp_desc_f64(
                left.candidate.host_fit_score.unwrap_or(0.0),
                right.candidate.host_fit_score.unwrap_or(0.0),
            )
            .then_with(|| cmp_desc_f64(left.coverage_fraction, right.coverage_fraction))
            .then_with(|| cmp_asc_usize(left.warning_count, right.warning_count))
            .then_with(|| cmp_asc_usize(left.strategy_rank, right.strategy_rank)),
        }
    }
}

fn warnings_to_json(lines: &[String]) -> Vec<String> {
    lines
        .iter()
        .map(|line| line.trim())
        .filter(|line| !line.is_empty())
        .map(|line| line.to_string())
        .collect()
}
