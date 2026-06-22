//! Sequence/digestion and primer/overhang utility routines used by core
//! operations.
//!
//! Keep low-level sequence manipulation helpers here so restriction/ligation,
//! primer design, and export paths reuse one implementation.
//!
//! Look here for:
//! - restriction digest and fragment/overhang handling
//! - primer/overhang utility routines reused by cloning workflows
//! - shared sequence export/pool materialization helpers

use super::*;

pub(crate) struct PrimerDesignProgressContext<'a> {
    pub(crate) seq_id: &'a str,
    pub(crate) design_kind: &'static str,
    pub(crate) backend_requested: &'a str,
    pub(crate) backend_used: &'a str,
    pub(crate) roi_start_0based: usize,
    pub(crate) roi_end_0based_exclusive: usize,
    pub(crate) max_output: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ProjectFactWorld {
    Closed,
    Open,
}

struct ProjectFactTypeSpec {
    world: ProjectFactWorld,
    requires_basis: bool,
}

fn project_fact_type_spec(name: &str) -> Option<ProjectFactTypeSpec> {
    let (world, requires_basis) = match name {
        "sequence.exists" | "sequence.kind" | "sequence.length" | "sequence.circular" => {
            (ProjectFactWorld::Closed, false)
        }
        "ui.host_available" => (ProjectFactWorld::Closed, false),
        "report.exists" | "restriction_site.present" | "restriction_site.absent" => {
            (ProjectFactWorld::Open, true)
        }
        _ => return None,
    };
    Some(ProjectFactTypeSpec {
        world,
        requires_basis,
    })
}

fn fact_subject(kind: FactSubjectKind, id: impl Into<String>) -> FactSubject {
    FactSubject {
        kind,
        id: id.into(),
    }
}

fn fact_subject_sort_key(subject: &FactSubject) -> (&'static str, &str) {
    let kind = match subject.kind {
        FactSubjectKind::Sequence => "sequence",
        FactSubjectKind::Report => "report",
        FactSubjectKind::Container => "container",
        FactSubjectKind::Ui => "ui",
        FactSubjectKind::Other => "other",
    };
    (kind, subject.id.as_str())
}

fn fact_value_key(value: &Option<serde_json::Value>) -> String {
    value
        .as_ref()
        .map(|value| serde_json::to_string(value).unwrap_or_default())
        .unwrap_or_default()
}

fn fact_range_key(range: &Option<FactRange>) -> String {
    range
        .as_ref()
        .map(|range| serde_json::to_string(range).unwrap_or_default())
        .unwrap_or_default()
}

fn sort_project_facts(facts: &mut Vec<ProjectFact>) {
    facts.sort_by(|a, b| {
        a.fact
            .cmp(&b.fact)
            .then(fact_subject_sort_key(&a.subject).cmp(&fact_subject_sort_key(&b.subject)))
            .then(a.enzyme.cmp(&b.enzyme))
            .then(fact_range_key(&a.range).cmp(&fact_range_key(&b.range)))
            .then(fact_value_key(&a.value).cmp(&fact_value_key(&b.value)))
            .then(
                a.basis
                    .as_ref()
                    .map(|basis| basis.report_id.as_str())
                    .unwrap_or("")
                    .cmp(
                        b.basis
                            .as_ref()
                            .map(|basis| basis.report_id.as_str())
                            .unwrap_or(""),
                    ),
            )
    });
}

fn fact_range_contains(proven: Option<&FactRange>, required: Option<&FactRange>) -> bool {
    match (proven, required) {
        (_, None) => true,
        (None, Some(_)) => false,
        (Some(FactRange::WholeSequence), Some(_)) => true,
        (Some(FactRange::Span { .. }), Some(FactRange::WholeSequence)) => false,
        (
            Some(FactRange::Span {
                start_0based: proven_start,
                end_0based_exclusive: proven_end,
                topology: proven_topology,
            }),
            Some(FactRange::Span {
                start_0based: required_start,
                end_0based_exclusive: required_end,
                topology: required_topology,
            }),
        ) => {
            proven_topology == required_topology
                && proven_start <= required_start
                && proven_end >= required_end
        }
    }
}

fn fact_range_satisfies_atom(
    fact_name: &str,
    proven: Option<&FactRange>,
    required: Option<&FactRange>,
) -> bool {
    if fact_name == "restriction_site.present" {
        return required
            .map(|required| fact_range_contains(Some(required), proven))
            .unwrap_or(true);
    }
    fact_range_contains(proven, required)
}

fn normalize_fact_atom_subject(atom: &FactAtom) -> Option<FactSubject> {
    atom.subject.clone().or_else(|| {
        atom.id
            .as_ref()
            .map(|id| fact_subject(FactSubjectKind::Sequence, id.trim().to_string()))
    })
}

fn fact_atom_subject_matches(atom: &FactAtom, fact: &ProjectFact) -> bool {
    normalize_fact_atom_subject(atom)
        .map(|subject| subject == fact.subject)
        .unwrap_or(true)
}

fn fact_atom_scalar_matches(atom: &FactAtom, fact: &ProjectFact) -> bool {
    if let Some(expected) = &atom.equals
        && fact.value.as_ref() != Some(expected)
    {
        return false;
    }
    if let Some(compare) = &atom.compare {
        let Some(actual) = fact.value.as_ref().and_then(serde_json::Value::as_f64) else {
            return false;
        };
        let Some(expected) = compare.value.as_f64() else {
            return false;
        };
        let ok = match compare.op.as_str() {
            "eq" | "=" | "==" => (actual - expected).abs() < f64::EPSILON,
            "ne" | "!=" => (actual - expected).abs() >= f64::EPSILON,
            "gt" | ">" => actual > expected,
            "gte" | ">=" => actual >= expected,
            "lt" | "<" => actual < expected,
            "lte" | "<=" => actual <= expected,
            _ => false,
        };
        if !ok {
            return false;
        }
    }
    true
}

fn fact_atom_matches_fact(atom: &FactAtom, fact: &ProjectFact) -> bool {
    atom.fact == fact.fact
        && fact_atom_subject_matches(atom, fact)
        && atom
            .enzyme
            .as_ref()
            .map(|enzyme| {
                fact.enzyme
                    .as_ref()
                    .map(|value| value.eq_ignore_ascii_case(enzyme))
                    .unwrap_or(false)
            })
            .unwrap_or(true)
        && fact_range_satisfies_atom(&fact.fact, fact.range.as_ref(), atom.range.as_ref())
        && fact_atom_scalar_matches(atom, fact)
}

fn contradiction_fact_name(name: &str) -> Option<&'static str> {
    match name {
        "restriction_site.present" => Some("restriction_site.absent"),
        "restriction_site.absent" => Some("restriction_site.present"),
        _ => None,
    }
}

impl GentleEngine {
    pub(crate) fn summarize_process_run_bundle_construct_reasoning_graph(
        graph: &ConstructReasoningGraph,
    ) -> ProcessRunBundleConstructReasoningSummary {
        let fact_types = graph
            .facts
            .iter()
            .map(|fact| fact.fact_type.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let decision_types = graph
            .decisions
            .iter()
            .map(|node| node.decision_type.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let mut fact_statuses: BTreeMap<String, String> = BTreeMap::new();
        for fact in &graph.facts {
            if let Some(status) = fact
                .value_json
                .get("status")
                .and_then(serde_json::Value::as_str)
            {
                fact_statuses.insert(fact.fact_type.clone(), status.to_string());
            }
        }

        let mut host_profile_ids = graph
            .objective
            .propagation_host_profile_id
            .iter()
            .chain(graph.objective.expression_host_profile_id.iter())
            .cloned()
            .collect::<Vec<_>>();
        host_profile_ids.extend(
            graph
                .objective
                .host_route
                .iter()
                .map(|step| step.host_profile_id.clone()),
        );
        host_profile_ids.sort();
        host_profile_ids.dedup();

        let growth_condition_signals = graph
            .facts
            .iter()
            .find(|fact| fact.fact_type == "growth_condition_context")
            .and_then(|fact| {
                fact.value_json
                    .get("condition_signals")
                    .and_then(serde_json::Value::as_array)
            })
            .map(|rows| {
                rows.iter()
                    .filter_map(|row| {
                        row.get("label")
                            .or_else(|| row.get("token"))
                            .and_then(serde_json::Value::as_str)
                            .map(|value| value.to_string())
                    })
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let supported_selection_rule_ids = graph
            .facts
            .iter()
            .find(|fact| fact.fact_type == "selection_context")
            .and_then(|fact| {
                fact.value_json
                    .get("selection_rule_matches")
                    .and_then(serde_json::Value::as_array)
            })
            .map(|rows| {
                rows.iter()
                    .filter(|row| {
                        row.get("status").and_then(serde_json::Value::as_str) == Some("supported")
                    })
                    .filter_map(|row| {
                        row.get("rule_id")
                            .and_then(serde_json::Value::as_str)
                            .map(|value| value.to_string())
                    })
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let variant_effect_tags = graph
            .facts
            .iter()
            .find(|fact| fact.fact_type == "variant_effect_context")
            .map(|fact| {
                fact.value_json
                    .get("effect_tags")
                    .and_then(serde_json::Value::as_array)
                    .map(|rows| {
                        rows.iter()
                            .filter_map(serde_json::Value::as_str)
                            .map(|value| value.to_string())
                            .collect::<BTreeSet<_>>()
                            .into_iter()
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default()
            })
            .unwrap_or_default();
        let suggested_variant_assay_ids = graph
            .facts
            .iter()
            .find(|fact| fact.fact_type == "variant_assay_context")
            .map(|fact| {
                fact.value_json
                    .get("suggested_assay_ids")
                    .and_then(serde_json::Value::as_array)
                    .map(|rows| {
                        rows.iter()
                            .filter_map(serde_json::Value::as_str)
                            .map(|value| value.to_string())
                            .collect::<BTreeSet<_>>()
                            .into_iter()
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default()
            })
            .unwrap_or_default();
        let annotation_candidate_status_counts = graph.annotation_candidates.iter().fold(
            BTreeMap::<&'static str, usize>::new(),
            |mut counts, candidate| {
                let key = match candidate.editable_status {
                    EditableStatus::Draft => "draft",
                    EditableStatus::Accepted => "accepted",
                    EditableStatus::Rejected => "rejected",
                    EditableStatus::Locked => "locked",
                };
                *counts.entry(key).or_insert(0) += 1;
                counts
            },
        );

        let mut summary_lines = vec![];
        if !graph.objective.goal.trim().is_empty() {
            summary_lines.push(format!("Goal: {}", graph.objective.goal.trim()));
        }
        if !host_profile_ids.is_empty() {
            summary_lines.push(format!("Host profiles: {}", host_profile_ids.join(", ")));
        }
        if let Some(helper_profile_id) = graph.objective.helper_profile_id.as_ref() {
            summary_lines.push(format!("Helper profile: {helper_profile_id}"));
        }
        if !graph.objective.medium_conditions.is_empty() {
            summary_lines.push(format!(
                "Medium/growth conditions: {}",
                graph.objective.medium_conditions.join(", ")
            ));
        }
        for fact_type in [
            "host_transition_context",
            "host_restriction_methylation_context",
            "adapter_restriction_capture_context",
            "growth_condition_context",
            "helper_context",
            "selection_context",
            "variant_effect_context",
            "variant_assay_context",
            "pcr_operational_risk_context",
            "nanopore_operational_risk_context",
            "mapping_operational_risk_context",
            "cloning_stability_context",
        ] {
            if let Some(fact) = graph.facts.iter().find(|fact| fact.fact_type == fact_type) {
                if let Some(status) = fact_statuses.get(fact_type) {
                    summary_lines.push(format!("{} ({status})", fact.label.trim()));
                } else if !fact.label.trim().is_empty() {
                    summary_lines.push(fact.label.trim().to_string());
                }
            }
        }
        if !growth_condition_signals.is_empty() {
            summary_lines.push(format!(
                "Condition signals: {}",
                growth_condition_signals.join(", ")
            ));
        }
        if !supported_selection_rule_ids.is_empty() {
            summary_lines.push(format!(
                "Supported selection rules: {}",
                supported_selection_rule_ids.join(", ")
            ));
        }
        if !variant_effect_tags.is_empty() {
            summary_lines.push(format!(
                "Variant effect tags: {}",
                variant_effect_tags.join(", ")
            ));
        }
        if !suggested_variant_assay_ids.is_empty() {
            summary_lines.push(format!(
                "Suggested variant assays: {}",
                suggested_variant_assay_ids.join(", ")
            ));
        }
        if !graph.annotation_candidates.is_empty() {
            let mut status_parts = vec![];
            for key in ["accepted", "rejected", "draft", "locked"] {
                if let Some(count) = annotation_candidate_status_counts.get(key) {
                    status_parts.push(format!("{count} {key}"));
                }
            }
            if !status_parts.is_empty() {
                summary_lines.push(format!(
                    "Annotation candidates: {}",
                    status_parts.join(", ")
                ));
            }
        }
        let (protein_handoff_candidate_count, protein_handoff_source_protein_seq_ids) =
            Self::construct_reasoning_graph_protein_to_dna_handoff_summary(graph);
        if protein_handoff_candidate_count > 0 {
            summary_lines.push(format!(
                "Protein-to-DNA handoff candidates: {}",
                protein_handoff_candidate_count
            ));
        }
        if !protein_handoff_source_protein_seq_ids.is_empty() {
            summary_lines.push(format!(
                "Protein handoff sources: {}",
                protein_handoff_source_protein_seq_ids.join(", ")
            ));
        }
        summary_lines.retain(|line| !line.trim().is_empty());
        summary_lines.dedup();

        let mut warning_lines = graph
            .evidence
            .iter()
            .flat_map(|row| row.warnings.iter())
            .chain(
                graph
                    .facts
                    .iter()
                    .filter(|fact| {
                        fact.value_json
                            .get("status")
                            .and_then(serde_json::Value::as_str)
                            .is_some_and(|value| value.starts_with("review_needed"))
                    })
                    .map(|fact| &fact.label),
            )
            .filter_map(|line| {
                let trimmed = line.trim();
                (!trimmed.is_empty()).then(|| trimmed.to_string())
            })
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        warning_lines.sort();

        ProcessRunBundleConstructReasoningSummary {
            seq_id: graph.seq_id.clone(),
            graph_id: graph.graph_id.clone(),
            objective_id: graph.objective.objective_id.clone(),
            objective_title: graph.objective.title.clone(),
            objective_goal: graph.objective.goal.clone(),
            fact_types,
            decision_types,
            fact_statuses,
            host_profile_ids,
            helper_profile_id: graph.objective.helper_profile_id.clone(),
            medium_conditions: graph.objective.medium_conditions.clone(),
            growth_condition_signals,
            supported_selection_rule_ids,
            variant_effect_tags,
            suggested_variant_assay_ids,
            summary_lines,
            warning_lines,
        }
    }

    fn collect_process_run_bundle_construct_reasoning(
        &self,
        seq_ids: &BTreeSet<String>,
    ) -> ProcessRunBundleConstructReasoningExport {
        let seq_ids_considered = seq_ids.iter().cloned().collect::<Vec<_>>();
        let mut graphs = seq_ids_considered
            .iter()
            .filter_map(|seq_id| self.construct_reasoning_graph_for_seq_id(seq_id).ok())
            .collect::<Vec<_>>();
        graphs.sort_by(|left, right| {
            left.seq_id
                .to_ascii_lowercase()
                .cmp(&right.seq_id.to_ascii_lowercase())
                .then_with(|| left.graph_id.cmp(&right.graph_id))
        });
        let summaries = graphs
            .iter()
            .map(Self::summarize_process_run_bundle_construct_reasoning_graph)
            .collect::<Vec<_>>();
        ProcessRunBundleConstructReasoningExport {
            seq_ids_considered,
            summaries,
            graphs,
        }
    }

    pub(crate) fn infer_gel_topology_form_from_dna(
        dna: &DNAsequence,
    ) -> gentle_protocol::GelTopologyForm {
        if !dna.is_circular() {
            return gentle_protocol::GelTopologyForm::Linear;
        }

        let mut text = String::new();
        if let Some(name) = dna.name() {
            text.push_str(name);
            text.push(' ');
        }
        if let Some(definition) = dna.definition() {
            text.push_str(definition);
            text.push(' ');
        }
        for comment in dna.description() {
            text.push_str(comment);
            text.push(' ');
        }
        let lowered = text.to_ascii_lowercase();

        for token in lowered
            .split(|c: char| c.is_whitespace() || matches!(c, ',' | ';' | '(' | ')' | '[' | ']'))
        {
            let trimmed = token.trim_matches(|c: char| {
                !c.is_ascii_alphanumeric() && c != '_' && c != '=' && c != '-'
            });
            if let Some(value) = trimmed.strip_prefix("gel_topology=")
                && let Some(form) = gentle_protocol::GelTopologyForm::from_hint(value)
            {
                return form;
            }
            if let Some(form) = gentle_protocol::GelTopologyForm::from_hint(trimmed)
                && form.is_circular()
            {
                return form;
            }
        }
        if lowered.contains("open circular")
            || lowered.contains("open-circle")
            || lowered.contains("nicked circular")
        {
            return gentle_protocol::GelTopologyForm::NickedCircular;
        }
        gentle_protocol::GelTopologyForm::Circular
    }

    pub(super) fn digest_with_guard(
        dna: &DNAsequence,
        enzymes: Vec<RestrictionEnzyme>,
        max_fragments: usize,
    ) -> Result<Vec<DNAsequence>, EngineError> {
        let mut fragments = vec![dna.clone()];
        for enzyme in &enzymes {
            println!("Digesting with enzyme: {}", enzyme.name);
            let mut seen_states: HashSet<u64> = HashSet::new();
            let mut rounds: usize = 0;
            let mut last_fragment_count = fragments.len();
            let enzyme_started = Instant::now();
            // Conservative guard against non-converging digest loops.
            let max_rounds = max_fragments.min(1_024).max(64);
            loop {
                rounds += 1;
                if enzyme_started.elapsed().as_millis() > 750 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest timed out for enzyme '{}' (>{} ms)",
                            enzyme.name, 750
                        ),

                        cause_chain: vec![],
                    });
                }
                if rounds > max_rounds {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest did not converge for enzyme '{}' within {} rounds",
                            enzyme.name, max_rounds
                        ),

                        cause_chain: vec![],
                    });
                }

                // Detect cyclic/non-converging fragment states.
                let mut state_hasher = DefaultHasher::new();
                fragments.len().hash(&mut state_hasher);
                for seq in &fragments {
                    seq.get_forward_string().hash(&mut state_hasher);
                    seq.overhang().forward_3.hash(&mut state_hasher);
                    seq.overhang().forward_5.hash(&mut state_hasher);
                    seq.overhang().reverse_3.hash(&mut state_hasher);
                    seq.overhang().reverse_5.hash(&mut state_hasher);
                }
                let state_sig = state_hasher.finish();
                if !seen_states.insert(state_sig) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest entered a repeated state for enzyme '{}'",
                            enzyme.name
                        ),

                        cause_chain: vec![],
                    });
                }

                let mut found_one = false;
                let mut new_fragments: Vec<DNAsequence> = vec![];
                for seq in fragments.drain(..) {
                    if let Some(site) = enzyme.get_sites(&seq, None).first() {
                        let split = seq.split_at_restriction_enzyme_site(site);
                        found_one = true;
                        new_fragments.extend(split);
                    } else {
                        new_fragments.push(seq);
                    }

                    if new_fragments.len() > max_fragments {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Digest produced more than max_fragments_per_container={}",
                                max_fragments
                            ),

                            cause_chain: vec![],
                        });
                    }
                }
                fragments = new_fragments;
                let current_count = fragments.len();

                if !found_one {
                    println!(
                        "Digest enzyme '{}' completed in {} round(s), fragments: {}",
                        enzyme.name, rounds, current_count
                    );
                    break;
                }
                // For linear-digest progression, total fragment count should increase.
                // If not, we are likely cutting in a pathological cycle.
                if current_count <= last_fragment_count && rounds > 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Digest stalled for enzyme '{}' (no fragment-count progress)",
                            enzyme.name
                        ),

                        cause_chain: vec![],
                    });
                }
                last_fragment_count = current_count;
            }
        }
        Ok(fragments)
    }

    pub(super) fn max_fragments_per_container(&self) -> usize {
        self.state.parameters.max_fragments_per_container
    }

    pub(super) fn container_members(&self, container_id: &str) -> Result<Vec<SeqId>, EngineError> {
        let container = self
            .state
            .container_state
            .containers
            .get(container_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Container '{container_id}' not found"),

                cause_chain: vec![],
            })?;
        if container.members.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Container '{container_id}' has no members"),

                cause_chain: vec![],
            });
        }
        for seq_id in &container.members {
            if !self.state.sequences.contains_key(seq_id) {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Container '{container_id}' references unknown sequence '{seq_id}'"
                    ),

                    cause_chain: vec![],
                });
            }
        }
        Ok(container.members.clone())
    }

    pub(super) fn gel_samples_from_container_ids(
        &self,
        container_ids: &[ContainerId],
    ) -> Result<Vec<GelSampleInput>, EngineError> {
        if container_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "At least one container id is required for gel rendering".to_string(),

                cause_chain: vec![],
            });
        }
        let mut samples: Vec<GelSampleInput> = vec![];
        for container_id in container_ids {
            let container = self
                .state
                .container_state
                .containers
                .get(container_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Container '{container_id}' not found"),

                    cause_chain: vec![],
                })?;
            if container.members.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Container '{container_id}' has no members"),

                    cause_chain: vec![],
                });
            }
            let mut members: Vec<crate::pool_gel::GelSampleMember> =
                Vec::with_capacity(container.members.len());
            for seq_id in &container.members {
                let dna = self
                    .state
                    .sequences
                    .get(seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!(
                            "Container '{container_id}' references unknown sequence '{seq_id}'"
                        ),

                        cause_chain: vec![],
                    })?;
                members.push(crate::pool_gel::GelSampleMember {
                    seq_id: seq_id.clone(),
                    bp: dna.len(),
                    topology_form: Self::infer_gel_topology_form_from_dna(dna),
                });
            }
            let lane_name = container
                .name
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(ToString::to_string)
                .unwrap_or_else(|| container.container_id.clone());
            samples.push(GelSampleInput {
                name: lane_name,
                role_label: None,
                members,
            });
        }
        Ok(samples)
    }

    pub(super) fn gel_samples_from_arrangement(
        &self,
        arrangement_id: &str,
    ) -> Result<(Vec<GelSampleInput>, Vec<String>), EngineError> {
        let arrangement = self
            .state
            .container_state
            .arrangements
            .get(arrangement_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Arrangement '{arrangement_id}' not found"),

                cause_chain: vec![],
            })?;
        if arrangement.mode != ArrangementMode::Serial {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Arrangement '{}' is mode '{:?}', only serial arrangements can render gels",
                    arrangement_id, arrangement.mode
                ),

                cause_chain: vec![],
            });
        }
        let mut samples = self.gel_samples_from_container_ids(&arrangement.lane_container_ids)?;
        for (idx, sample) in samples.iter_mut().enumerate() {
            sample.role_label = arrangement
                .lane_role_labels
                .get(idx)
                .map(|label| label.trim())
                .filter(|label| !label.is_empty())
                .map(ToString::to_string);
        }
        Ok((samples, arrangement.ladders.clone()))
    }

    /// Resolve one serial-gel layout from direct inputs, containers, or one
    /// stored arrangement, optionally overriding ladder selection.
    pub fn build_serial_gel_layout_for_render(
        &self,
        inputs: &[SeqId],
        container_ids: Option<&[ContainerId]>,
        arrangement_id: Option<&str>,
        ladder_override: Option<&[String]>,
        conditions: Option<&gentle_protocol::GelRunConditions>,
    ) -> Result<crate::pool_gel::PoolGelLayout, EngineError> {
        let mut ladder_names = ladder_override
            .map(Self::normalize_serial_gel_ladders_slice)
            .unwrap_or_default();
        let samples: Vec<GelSampleInput> =
            if let Some(arrangement_id) = arrangement_id.map(str::trim) {
                if arrangement_id.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "arrangement_id cannot be empty".to_string(),

                        cause_chain: vec![],
                    });
                }
                let (arrangement_samples, arrangement_ladders) =
                    self.gel_samples_from_arrangement(arrangement_id)?;
                if ladder_override.is_none() {
                    ladder_names = arrangement_ladders;
                }
                arrangement_samples
            } else if let Some(container_ids) = container_ids {
                if container_ids.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "container_ids was provided but empty".to_string(),

                        cause_chain: vec![],
                    });
                }
                self.gel_samples_from_container_ids(container_ids)?
            } else {
                if inputs.is_empty() {
                    return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message:
                        "RenderPoolGelSvg requires either inputs, container_ids, or arrangement_id"
                            .to_string(),
                
                    cause_chain: vec![],});
                }
                let mut members: Vec<crate::pool_gel::GelSampleMember> =
                    Vec::with_capacity(inputs.len());
                for seq_id in inputs {
                    let dna = self
                        .state
                        .sequences
                        .get(seq_id)
                        .ok_or_else(|| EngineError {
                            code: ErrorCode::NotFound,
                            message: format!("Sequence '{seq_id}' not found"),

                            cause_chain: vec![],
                        })?;
                    members.push(crate::pool_gel::GelSampleMember {
                        seq_id: seq_id.clone(),
                        bp: dna.len(),
                        topology_form: Self::infer_gel_topology_form_from_dna(dna),
                    });
                }
                vec![GelSampleInput {
                    name: format!("Input tube (n={})", members.len()),
                    role_label: None,
                    members,
                }]
            };
        crate::pool_gel::build_serial_gel_layout(&samples, &ladder_names, conditions).map_err(|e| {
            EngineError {
                code: ErrorCode::InvalidInput,
                message: e,

                cause_chain: vec![],
            }
        })
    }

    pub(super) fn flatten_container_members(
        &self,
        container_ids: &[ContainerId],
    ) -> Result<Vec<SeqId>, EngineError> {
        if container_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "At least one container id is required".to_string(),

                cause_chain: vec![],
            });
        }
        let mut members = Vec::new();
        for container_id in container_ids {
            members.extend(self.container_members(container_id)?);
            if members.len() > self.max_fragments_per_container() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Container merge input count exceeds max_fragments_per_container={}",
                        self.max_fragments_per_container()
                    ),

                    cause_chain: vec![],
                });
            }
        }
        if members.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "No sequences found in selected containers".to_string(),

                cause_chain: vec![],
            });
        }
        Ok(members)
    }

    pub(super) fn resolve_enzymes(
        &self,
        enzymes: &[String],
    ) -> Result<(Vec<RestrictionEnzyme>, Vec<String>), EngineError> {
        if enzymes.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Digest requires at least one enzyme".to_string(),

                cause_chain: vec![],
            });
        }

        let mut by_name: HashMap<String, RestrictionEnzyme> = HashMap::new();
        for enzyme in active_restriction_enzymes() {
            by_name.entry(enzyme.name.clone()).or_insert(enzyme);
        }
        for dna in self.state.sequences.values() {
            for enzyme in dna.restriction_enzymes() {
                by_name
                    .entry(enzyme.name.clone())
                    .or_insert_with(|| enzyme.clone());
            }
        }

        let mut found = Vec::new();
        let mut missing = Vec::new();
        let mut seen_names: HashSet<String> = HashSet::new();
        for name in enzymes {
            if let Some(enzyme) = by_name.get(name) {
                if seen_names.insert(name.clone()) {
                    found.push(enzyme.clone());
                }
            } else if !missing.contains(name) {
                missing.push(name.clone());
            }
        }

        if found.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "None of the requested enzymes are known: {}",
                    enzymes.join(",")
                ),

                cause_chain: vec![],
            });
        }
        Ok((found, missing))
    }

    pub fn summarize_state(&self) -> EngineStateSummary {
        let mut sequences: Vec<EngineSequenceSummary> = self
            .state
            .sequences
            .iter()
            .map(|(id, dna)| EngineSequenceSummary {
                id: id.to_string(),
                name: dna.name().clone(),
                length: dna.len(),
                circular: dna.is_circular(),
            })
            .collect();
        sequences.sort_by(|a, b| a.id.cmp(&b.id));

        let mut containers: Vec<EngineContainerSummary> = self
            .state
            .container_state
            .containers
            .iter()
            .map(|(id, c)| EngineContainerSummary {
                id: id.to_string(),
                kind: format!("{:?}", c.kind),
                declared_contents_exclusive: c.declared_contents_exclusive,
                member_count: c.members.len(),
                members: c.members.clone(),
            })
            .collect();
        containers.sort_by(|a, b| a.id.cmp(&b.id));

        let mut arrangements: Vec<EngineArrangementSummary> = self
            .state
            .container_state
            .arrangements
            .iter()
            .map(|(id, arrangement)| EngineArrangementSummary {
                id: id.to_string(),
                mode: format!("{:?}", arrangement.mode),
                lane_count: arrangement.lane_container_ids.len(),
                lane_container_ids: arrangement.lane_container_ids.clone(),
                ladders: arrangement.ladders.clone(),
            })
            .collect();
        arrangements.sort_by(|a, b| a.id.cmp(&b.id));

        EngineStateSummary {
            sequence_count: sequences.len(),
            sequences,
            container_count: containers.len(),
            containers,
            arrangement_count: arrangements.len(),
            arrangements,
            display: self.state.display.clone(),
        }
    }

    pub fn project_fact_graph(&self) -> ProjectFactGraph {
        self.project_fact_graph_with_restriction_evidence(&[])
    }

    pub fn project_fact_graph_with_restriction_evidence(
        &self,
        restriction_reports: &[RestrictionSiteScanReport],
    ) -> ProjectFactGraph {
        let mut facts: Vec<ProjectFact> = vec![];

        for (id, dna) in &self.state.sequences {
            let subject = fact_subject(FactSubjectKind::Sequence, id.clone());
            facts.push(ProjectFact {
                fact: "sequence.exists".to_string(),
                subject: subject.clone(),
                ..ProjectFact::default()
            });
            facts.push(ProjectFact {
                fact: "sequence.kind".to_string(),
                subject: subject.clone(),
                value: Some(serde_json::json!("dna")),
                ..ProjectFact::default()
            });
            facts.push(ProjectFact {
                fact: "sequence.length".to_string(),
                subject: subject.clone(),
                value: Some(serde_json::json!(dna.len())),
                ..ProjectFact::default()
            });
            facts.push(ProjectFact {
                fact: "sequence.circular".to_string(),
                subject,
                value: Some(serde_json::json!(dna.is_circular())),
                ..ProjectFact::default()
            });
        }

        for report in restriction_reports {
            Self::push_restriction_scan_facts(&mut facts, report);
        }

        sort_project_facts(&mut facts);
        ProjectFactGraph {
            schema: PROJECT_FACT_GRAPH_SCHEMA.to_string(),
            facts,
        }
    }

    fn push_restriction_scan_facts(
        facts: &mut Vec<ProjectFact>,
        report: &RestrictionSiteScanReport,
    ) {
        let report_id = if report.report_id.trim().is_empty() {
            Self::restriction_site_scan_report_id(report)
        } else {
            report.report_id.trim().to_string()
        };
        let report_subject = fact_subject(FactSubjectKind::Report, report_id.clone());
        let basis = FactBasis {
            report_id: report_id.clone(),
            report_kind: "restriction_scan".to_string(),
            evidence_class: EvidenceClass::HardFact,
            op_id: report.op_id.clone(),
            run_id: report.run_id.clone(),
        };
        facts.push(ProjectFact {
            fact: "report.exists".to_string(),
            subject: report_subject,
            value: Some(serde_json::json!("restriction_scan")),
            basis: Some(basis.clone()),
            ..ProjectFact::default()
        });

        if report.target_kind != "seq_id" {
            return;
        }
        let subject = fact_subject(FactSubjectKind::Sequence, report.target_label.clone());
        let range = Self::restriction_report_fact_range(report);
        let skipped = report
            .skipped_enzyme_names_due_to_max_sites
            .iter()
            .map(|name| name.to_ascii_lowercase())
            .collect::<BTreeSet<_>>();

        for enzyme in &report.enzymes_scanned {
            if skipped.contains(&enzyme.to_ascii_lowercase()) {
                continue;
            }
            let matching_hits = report
                .rows
                .iter()
                .filter(|row| row.enzyme_name.eq_ignore_ascii_case(enzyme))
                .collect::<Vec<_>>();
            if matching_hits.is_empty() {
                facts.push(ProjectFact {
                    fact: "restriction_site.absent".to_string(),
                    subject: subject.clone(),
                    enzyme: Some(enzyme.clone()),
                    range: Some(range.clone()),
                    value: Some(serde_json::json!(0)),
                    basis: Some(basis.clone()),
                });
            } else {
                for hit in matching_hits {
                    facts.push(ProjectFact {
                        fact: "restriction_site.present".to_string(),
                        subject: subject.clone(),
                        enzyme: Some(hit.enzyme_name.clone()),
                        range: Some(FactRange::Span {
                            start_0based: hit.source_recognition_start_0based,
                            end_0based_exclusive: hit.source_recognition_end_0based_exclusive,
                            topology: report.scan_topology,
                        }),
                        value: Some(serde_json::json!(true)),
                        basis: Some(basis.clone()),
                    });
                }
            }
        }
    }

    fn restriction_report_fact_range(report: &RestrictionSiteScanReport) -> FactRange {
        if report.scan_start_0based == 0
            && report.scan_end_0based_exclusive == report.source_sequence_length_bp
        {
            FactRange::WholeSequence
        } else {
            FactRange::Span {
                start_0based: report.scan_start_0based,
                end_0based_exclusive: report.scan_end_0based_exclusive,
                topology: report.scan_topology,
            }
        }
    }

    pub fn restriction_site_scan_report_id(report: &RestrictionSiteScanReport) -> String {
        let mut enzymes = report.enzymes_scanned.clone();
        enzymes.sort();
        let raw = format!(
            "{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}",
            report.op_id.as_deref().unwrap_or(""),
            report.run_id.as_deref().unwrap_or(""),
            report.target_kind,
            report.target_label,
            report.scan_start_0based,
            report.scan_end_0based_exclusive,
            report.scan_topology.as_str(),
            enzymes.join(",")
        );
        Self::short_sha1_id("rsr", &raw)
    }

    pub fn evaluate_fact_expression_against_graph(
        expression: &FactExpression,
        graph: &ProjectFactGraph,
    ) -> FactEvaluationResult {
        let mut unmet_atoms = vec![];
        let mut unknown_atoms = vec![];
        let truth = Self::evaluate_fact_expression_inner(
            expression,
            graph,
            &mut unmet_atoms,
            &mut unknown_atoms,
        );
        FactEvaluationResult {
            schema: FACT_EVALUATION_SCHEMA.to_string(),
            truth,
            unmet_atoms,
            unknown_atoms,
        }
    }

    pub fn evaluate_fact_expression(
        &self,
        expression: &FactExpression,
        restriction_reports: &[RestrictionSiteScanReport],
    ) -> FactEvaluationResult {
        let graph = self.project_fact_graph_with_restriction_evidence(restriction_reports);
        Self::evaluate_fact_expression_against_graph(expression, &graph)
    }

    fn evaluate_fact_expression_inner(
        expression: &FactExpression,
        graph: &ProjectFactGraph,
        unmet_atoms: &mut Vec<FactAtom>,
        unknown_atoms: &mut Vec<FactAtom>,
    ) -> FactTruth {
        match expression {
            FactExpression::Atom(atom) => {
                let truth = Self::evaluate_fact_atom(atom, graph);
                match truth {
                    FactTruth::Unsatisfied => unmet_atoms.push(atom.clone()),
                    FactTruth::Unknown => unknown_atoms.push(atom.clone()),
                    FactTruth::Satisfied => {}
                }
                truth
            }
            FactExpression::All { all } => {
                let mut saw_unknown = false;
                for child in all {
                    match Self::evaluate_fact_expression_inner(
                        child,
                        graph,
                        unmet_atoms,
                        unknown_atoms,
                    ) {
                        FactTruth::Unsatisfied => return FactTruth::Unsatisfied,
                        FactTruth::Unknown => saw_unknown = true,
                        FactTruth::Satisfied => {}
                    }
                }
                if saw_unknown {
                    FactTruth::Unknown
                } else {
                    FactTruth::Satisfied
                }
            }
            FactExpression::Any { any } => {
                let mut saw_unknown = false;
                for child in any {
                    match Self::evaluate_fact_expression_inner(
                        child,
                        graph,
                        unmet_atoms,
                        unknown_atoms,
                    ) {
                        FactTruth::Satisfied => return FactTruth::Satisfied,
                        FactTruth::Unknown => saw_unknown = true,
                        FactTruth::Unsatisfied => {}
                    }
                }
                if saw_unknown {
                    FactTruth::Unknown
                } else {
                    FactTruth::Unsatisfied
                }
            }
            FactExpression::Not { not } => {
                let before_unmet = unmet_atoms.len();
                match Self::evaluate_fact_expression_inner(not, graph, unmet_atoms, unknown_atoms) {
                    FactTruth::Satisfied => FactTruth::Unsatisfied,
                    FactTruth::Unsatisfied => {
                        unmet_atoms.truncate(before_unmet);
                        FactTruth::Satisfied
                    }
                    FactTruth::Unknown => {
                        unmet_atoms.truncate(before_unmet);
                        FactTruth::Unknown
                    }
                }
            }
        }
    }

    fn evaluate_fact_atom(atom: &FactAtom, graph: &ProjectFactGraph) -> FactTruth {
        let Some(spec) = project_fact_type_spec(atom.fact.as_str()) else {
            return FactTruth::Unknown;
        };
        let has_match = graph.facts.iter().any(|fact| {
            fact_atom_matches_fact(atom, fact) && (!spec.requires_basis || fact.basis.is_some())
        });
        if has_match {
            return FactTruth::Satisfied;
        }
        if spec.world == ProjectFactWorld::Closed {
            return FactTruth::Unsatisfied;
        }
        if let Some(contradiction) = contradiction_fact_name(atom.fact.as_str()) {
            let mut contradiction_atom = atom.clone();
            contradiction_atom.fact = contradiction.to_string();
            let contradicted = graph.facts.iter().any(|fact| {
                fact_atom_matches_fact(&contradiction_atom, fact)
                    && (!spec.requires_basis || fact.basis.is_some())
            });
            if contradicted {
                return FactTruth::Unsatisfied;
            }
        }
        FactTruth::Unknown
    }

    pub(super) fn canonical_fasta_molecule(raw: Option<&str>) -> &'static str {
        let normalized = raw
            .unwrap_or("dsdna")
            .trim()
            .to_ascii_lowercase()
            .replace(['_', '-'], "")
            .replace(' ', "");
        match normalized.as_str() {
            "ssdna" | "singlestrandeddna" | "single" | "ssdnaoligo" | "ss" => "ssdna",
            "rna" | "ssrna" | "singlestrandedrna" | "transcript" | "mrna" | "cdna" => "rna",
            "protein" | "peptide" | "aa" | "aminoacid" | "aminoacids" => "protein",
            _ => "dsdna",
        }
    }

    pub(super) fn fasta_metadata_tokens(dna: &DNAsequence) -> Vec<String> {
        let molecule = Self::canonical_fasta_molecule(dna.molecule_type());
        let mut tokens = vec![
            format!("molecule={molecule}"),
            format!(
                "topology={}",
                if dna.is_circular() {
                    "circular"
                } else {
                    "linear"
                }
            ),
        ];

        if molecule == "dsdna" {
            let overhang = dna.overhang();
            if !overhang.forward_5.is_empty() {
                tokens.push(format!("f5={}", Self::overhang_text(&overhang.forward_5)));
            }
            if !overhang.forward_3.is_empty() {
                tokens.push(format!("f3={}", Self::overhang_text(&overhang.forward_3)));
            }
            if !overhang.reverse_5.is_empty() {
                tokens.push(format!("r5={}", Self::overhang_text(&overhang.reverse_5)));
            }
            if !overhang.reverse_3.is_empty() {
                tokens.push(format!("r3={}", Self::overhang_text(&overhang.reverse_3)));
            }
        }

        tokens
    }

    pub(super) fn save_as_fasta(
        seq_id: &str,
        dna: &DNAsequence,
        path: &str,
    ) -> Result<(), EngineError> {
        let mut file = File::create(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not create FASTA file '{path}': {e}"),

            cause_chain: vec![],
        })?;

        let header = dna
            .name()
            .clone()
            .unwrap_or_else(|| seq_id.to_string())
            .replace(' ', "_");
        let seq = dna.get_forward_string();

        let mut header_line = format!(">{header}");
        let metadata = Self::fasta_metadata_tokens(dna);
        if !metadata.is_empty() {
            header_line.push(' ');
            header_line.push_str(&metadata.join(" "));
        }

        writeln!(file, "{header_line}").map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write FASTA header to '{path}': {e}"),

            cause_chain: vec![],
        })?;

        for chunk in seq.as_bytes().chunks(80) {
            file.write_all(chunk).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write FASTA sequence to '{path}': {e}"),

                cause_chain: vec![],
            })?;
            file.write_all(b"\n").map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write FASTA newline to '{path}': {e}"),

                cause_chain: vec![],
            })?;
        }

        Ok(())
    }

    pub(super) fn overhang_text(v: &[u8]) -> String {
        String::from_utf8_lossy(v).to_string()
    }

    pub(super) fn infer_pool_end_type(end: &PoolEnd) -> String {
        let any = !end.forward_5.is_empty()
            || !end.forward_3.is_empty()
            || !end.reverse_5.is_empty()
            || !end.reverse_3.is_empty();
        if !any {
            return "blunt".to_string();
        }
        if !end.forward_5.is_empty() || !end.reverse_5.is_empty() {
            return "sticky_5p".to_string();
        }
        if !end.forward_3.is_empty() || !end.reverse_3.is_empty() {
            return "sticky_3p".to_string();
        }
        "unknown".to_string()
    }

    pub(super) fn default_pool_human_id(inputs: &[SeqId]) -> String {
        if inputs.len() <= 4 {
            format!("Pool({})", inputs.join(", "))
        } else {
            format!("Pool({} members, first: {})", inputs.len(), inputs[0])
        }
    }

    pub(super) fn export_pool_file(
        &self,
        inputs: &[SeqId],
        path: &str,
        pool_id: Option<String>,
        human_id: Option<String>,
    ) -> Result<(String, String, usize), EngineError> {
        if inputs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ExportPool requires at least one input sequence id".to_string(),

                cause_chain: vec![],
            });
        }
        let mut members: Vec<PoolMember> = Vec::with_capacity(inputs.len());
        for seq_id in inputs {
            let dna = self
                .state
                .sequences
                .get(seq_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Sequence '{seq_id}' not found"),

                    cause_chain: vec![],
                })?;
            let mut end = PoolEnd {
                end_type: String::new(),
                forward_5: Self::overhang_text(&dna.overhang().forward_5),
                forward_3: Self::overhang_text(&dna.overhang().forward_3),
                reverse_5: Self::overhang_text(&dna.overhang().reverse_5),
                reverse_3: Self::overhang_text(&dna.overhang().reverse_3),
            };
            end.end_type = Self::infer_pool_end_type(&end);
            members.push(PoolMember {
                seq_id: seq_id.clone(),
                human_id: dna.name().clone().unwrap_or_else(|| seq_id.clone()),
                name: dna.name().clone(),
                sequence: dna.get_forward_string(),
                length_bp: dna.len(),
                topology: Self::infer_gel_topology_form_from_dna(dna)
                    .as_str()
                    .to_string(),
                ends: end,
            });
        }

        let pool_id = pool_id.unwrap_or_else(|| "pool_export".to_string());
        let human_id = human_id.unwrap_or_else(|| Self::default_pool_human_id(inputs));
        let export = PoolExport {
            schema: "gentle.pool.v1".to_string(),
            pool_id: pool_id.clone(),
            human_id: human_id.clone(),
            member_count: members.len(),
            members,
        };
        let text = serde_json::to_string_pretty(&export).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize pool JSON: {e}"),

            cause_chain: vec![],
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write pool file '{path}': {e}"),

            cause_chain: vec![],
        })?;
        Ok((pool_id, human_id, export.member_count))
    }

    pub(super) fn push_unique_token(values: &mut Vec<String>, token: &str) {
        let trimmed = token.trim();
        if trimmed.is_empty() || values.iter().any(|existing| existing == trimmed) {
            return;
        }
        values.push(trimmed.to_string());
    }

    pub(super) fn operation_variant_name(op: &Operation) -> String {
        serde_json::to_value(op)
            .ok()
            .and_then(|value| value.as_object().and_then(|obj| obj.keys().next().cloned()))
            .unwrap_or_else(|| "UnknownOperation".to_string())
    }

    pub(super) fn summarize_run_bundle_operation_inputs(
        op: &Operation,
        op_id: &str,
        run_id: &str,
        record_index: usize,
    ) -> ProcessRunBundleOperationInputSummary {
        let mut summary = ProcessRunBundleOperationInputSummary {
            op_id: op_id.to_string(),
            run_id: run_id.to_string(),
            operation: Self::operation_variant_name(op),
            record_index,
            sequence_ids: vec![],
            container_ids: vec![],
            arrangement_ids: vec![],
            candidate_set_ids: vec![],
            guide_set_ids: vec![],
            genome_ids: vec![],
            file_paths: vec![],
        };
        match op {
            Operation::LoadFile { path, .. } => {
                Self::push_unique_token(&mut summary.file_paths, path);
            }
            Operation::CreateSequenceFromText { output_id, .. } => {
                if let Some(output_id) = output_id {
                    Self::push_unique_token(&mut summary.sequence_ids, output_id);
                }
            }
            Operation::SaveFile { seq_id, .. }
            | Operation::RenderSequenceSvg { seq_id, .. }
            | Operation::ExportSequenceContextBundle { seq_id, .. }
            | Operation::RenderDotplotSvg { seq_id, .. }
            | Operation::RenderFeatureExpertSvg { seq_id, .. }
            | Operation::RenderIsoformArchitectureSvg { seq_id, .. }
            | Operation::RenderRnaStructureSvg { seq_id, .. }
            | Operation::ExtendGenomeAnchor { seq_id, .. }
            | Operation::VerifyGenomeAnchor { seq_id, .. }
            | Operation::ImportGenomeBedTrack { seq_id, .. }
            | Operation::ImportGenomeBigWigTrack { seq_id, .. }
            | Operation::ImportGenomeVcfTrack { seq_id, .. }
            | Operation::ImportIsoformPanel { seq_id, .. }
            | Operation::ProjectUniprotToGenome { seq_id, .. }
            | Operation::QueryProteinResidueGenomicCoordinates { seq_id, .. }
            | Operation::ImportBlastHitsTrack { seq_id, .. }
            | Operation::GenerateCandidateSet { seq_id, .. }
            | Operation::GenerateCandidateSetBetweenAnchors { seq_id, .. }
            | Operation::DeriveTranscriptSequences { seq_id, .. }
            | Operation::ComputeDotplot { seq_id, .. }
            | Operation::ComputeFlexibilityTrack { seq_id, .. }
            | Operation::DeriveSplicingReferences { seq_id, .. }
            | Operation::InterpretRnaReads { seq_id, .. }
            | Operation::PreflightRnaReadIsoforms { seq_id, .. }
            | Operation::SetTopology { seq_id, .. }
            | Operation::RecomputeFeatures { seq_id, .. }
            | Operation::AnnotateTfbs { seq_id, .. } => {
                Self::push_unique_token(&mut summary.sequence_ids, seq_id);
            }
            Operation::FindRestrictionSites { target, .. }
            | Operation::ScanTfbsHits { target, .. }
            | Operation::SummarizeTfbsScoreTracks { target, .. }
            | Operation::RenderTfbsScoreTracksSvg { target, .. } => match target {
                SequenceScanTarget::SeqId { seq_id, .. } => {
                    Self::push_unique_token(&mut summary.sequence_ids, seq_id);
                }
                SequenceScanTarget::InlineSequence { .. } => {}
            },
            Operation::QueryRepeatAnnotations { genome_id, .. }
            | Operation::BuildRepeatEnvironmentCohort { genome_id, .. }
            | Operation::SummarizeWindowCohortTfbs {
                cohort: RepeatEnvironmentCohortReport { genome_id, .. },
                ..
            } => {
                Self::push_unique_token(&mut summary.sequence_ids, genome_id);
            }
            Operation::ComputeDotplotOverlay {
                owner_seq_id,
                reference_seq_id,
                queries,
                ..
            } => {
                Self::push_unique_token(&mut summary.sequence_ids, owner_seq_id);
                Self::push_unique_token(&mut summary.sequence_ids, reference_seq_id);
                for query in queries {
                    Self::push_unique_token(&mut summary.sequence_ids, &query.seq_id);
                }
            }
            Operation::AlignSequences {
                query,
                target,
                query_seq_id,
                target_seq_id,
                ..
            } => {
                if let Some(SequenceScanTarget::SeqId { seq_id, .. }) = query {
                    Self::push_unique_token(&mut summary.sequence_ids, seq_id);
                } else if let Some(seq_id) = query_seq_id.as_deref() {
                    Self::push_unique_token(&mut summary.sequence_ids, seq_id);
                }
                if let Some(SequenceScanTarget::SeqId { seq_id, .. }) = target {
                    Self::push_unique_token(&mut summary.sequence_ids, seq_id);
                } else if let Some(seq_id) = target_seq_id.as_deref() {
                    Self::push_unique_token(&mut summary.sequence_ids, seq_id);
                }
            }
            Operation::ListRnaReadReports { seq_id } => {
                if let Some(seq_id) = seq_id.as_deref() {
                    Self::push_unique_token(&mut summary.sequence_ids, seq_id);
                }
            }
            Operation::RenderPoolGelSvg {
                inputs,
                container_ids,
                arrangement_id,
                ..
            } => {
                for seq_id in inputs {
                    Self::push_unique_token(&mut summary.sequence_ids, seq_id);
                }
                if let Some(container_ids) = container_ids {
                    for container_id in container_ids {
                        Self::push_unique_token(&mut summary.container_ids, container_id);
                    }
                }
                if let Some(arrangement_id) = arrangement_id.as_deref() {
                    Self::push_unique_token(&mut summary.arrangement_ids, arrangement_id);
                }
            }
            Operation::CreateArrangementSerial { container_ids, .. }
            | Operation::MergeContainersById { container_ids, .. } => {
                for container_id in container_ids {
                    Self::push_unique_token(&mut summary.container_ids, container_id);
                }
            }
            Operation::SetContainerDeclaredContentsExclusive { container_id, .. } => {
                Self::push_unique_token(&mut summary.container_ids, container_id);
            }
            Operation::SetArrangementLadders { arrangement_id, .. } => {
                Self::push_unique_token(&mut summary.arrangement_ids, arrangement_id);
            }
            Operation::ExportPool { inputs, .. }
            | Operation::MergeContainers { inputs, .. }
            | Operation::Ligation { inputs, .. }
            | Operation::FilterByMolecularWeight { inputs, .. }
            | Operation::FilterByDesignConstraints { inputs, .. } => {
                for seq_id in inputs {
                    Self::push_unique_token(&mut summary.sequence_ids, seq_id);
                }
            }
            Operation::PrepareGenome { genome_id, .. }
            | Operation::ExtractGenomeRegion { genome_id, .. }
            | Operation::ExtractGenomeGene { genome_id, .. } => {
                Self::push_unique_token(&mut summary.genome_ids, genome_id);
            }
            Operation::AssessPrimerPairSpecificity {
                target_genome_id, ..
            } => {
                Self::push_unique_token(&mut summary.genome_ids, target_genome_id);
            }
            Operation::ImportUniprotSwissProt { path, .. } => {
                Self::push_unique_token(&mut summary.file_paths, path);
            }
            Operation::ListReporterCatalog {
                catalog_path, path, ..
            }
            | Operation::RecommendReporters {
                catalog_path, path, ..
            } => {
                if let Some(path) = catalog_path.as_deref() {
                    Self::push_unique_token(&mut summary.file_paths, path);
                }
                if let Some(path) = path.as_deref() {
                    Self::push_unique_token(&mut summary.file_paths, path);
                }
            }
            Operation::ExportReporterCorpus {
                catalog_path, path, ..
            } => {
                if let Some(path) = catalog_path.as_deref() {
                    Self::push_unique_token(&mut summary.file_paths, path);
                }
                Self::push_unique_token(&mut summary.file_paths, path);
            }
            Operation::PlanReporterConstructHandoff {
                candidate_set_path,
                catalog_path,
                reporter_backbone_load_path,
                path,
                ..
            } => {
                Self::push_unique_token(&mut summary.file_paths, candidate_set_path);
                if let Some(path) = catalog_path.as_deref() {
                    Self::push_unique_token(&mut summary.file_paths, path);
                }
                if let Some(path) = reporter_backbone_load_path.as_deref() {
                    Self::push_unique_token(&mut summary.file_paths, path);
                }
                if let Some(path) = path.as_deref() {
                    Self::push_unique_token(&mut summary.file_paths, path);
                }
            }
            Operation::ValidateProtocolCartoonTemplate { template_path } => {
                Self::push_unique_token(&mut summary.file_paths, template_path);
            }
            Operation::RenderProtocolCartoonTemplateWithBindingsSvg {
                template_path,
                bindings_path,
                ..
            } => {
                Self::push_unique_token(&mut summary.file_paths, template_path);
                Self::push_unique_token(&mut summary.file_paths, bindings_path);
            }
            Operation::FetchGenBankAccession { .. }
            | Operation::FetchDbSnpRegion { .. }
            | Operation::FetchUniprotLinkedGenBank { .. } => {}
            Operation::DigestContainer { container_id, .. }
            | Operation::LigationContainer { container_id, .. }
            | Operation::FilterContainerByMolecularWeight { container_id, .. } => {
                Self::push_unique_token(&mut summary.container_ids, container_id);
            }
            Operation::Digest { input, .. }
            | Operation::Pcr {
                template: input, ..
            }
            | Operation::PcrAdvanced {
                template: input, ..
            }
            | Operation::PcrMutagenesis {
                template: input, ..
            }
            | Operation::DesignPrimerPairs {
                template: input, ..
            }
            | Operation::DesignInsertionPrimerPairs {
                template: input, ..
            }
            | Operation::PcrOverlapExtensionMutagenesis {
                template: input, ..
            }
            | Operation::DesignQpcrAssays {
                template: input, ..
            }
            | Operation::TestCdnaPcr { seq_id: input, .. }
            | Operation::TestCdnaQpcr { seq_id: input, .. }
            | Operation::ExtractRegion { input, .. }
            | Operation::ExtractAnchoredRegion { input, .. }
            | Operation::SelectCandidate { input, .. }
            | Operation::Reverse { input, .. }
            | Operation::Complement { input, .. }
            | Operation::ReverseComplement { input, .. }
            | Operation::Branch { input, .. } => {
                Self::push_unique_token(&mut summary.sequence_ids, input);
            }
            Operation::DeleteCandidateSet { set_name }
            | Operation::ScoreCandidateSetExpression { set_name, .. }
            | Operation::ScoreCandidateSetDistance { set_name, .. }
            | Operation::ScoreCandidateSetWeightedObjective { set_name, .. } => {
                Self::push_unique_token(&mut summary.candidate_set_ids, set_name);
            }
            Operation::FilterCandidateSet {
                input_set,
                output_set,
                ..
            }
            | Operation::CandidateSetOp {
                left_set: input_set,
                right_set: output_set,
                ..
            }
            | Operation::TopKCandidateSet {
                input_set,
                output_set,
                ..
            }
            | Operation::ParetoFrontierCandidateSet {
                input_set,
                output_set,
                ..
            } => {
                Self::push_unique_token(&mut summary.candidate_set_ids, input_set);
                Self::push_unique_token(&mut summary.candidate_set_ids, output_set);
            }
            Operation::UpsertGuideSet { guide_set_id, .. }
            | Operation::DeleteGuideSet { guide_set_id }
            | Operation::FilterGuidesPractical { guide_set_id, .. }
            | Operation::GenerateGuideOligos { guide_set_id, .. }
            | Operation::ExportGuideOligos { guide_set_id, .. }
            | Operation::ExportGuideProtocolText { guide_set_id, .. } => {
                Self::push_unique_token(&mut summary.guide_set_ids, guide_set_id);
            }
            _ => {}
        }
        if let Operation::RenderIsoformArchitectureSvg {
            expression_tsv_path: Some(path),
            ..
        } = op
        {
            Self::push_unique_token(&mut summary.file_paths, path);
        }
        if let Operation::ImportGenomeBedTrack { path, .. }
        | Operation::ImportGenomeBigWigTrack { path, .. }
        | Operation::ImportGenomeVcfTrack { path, .. }
        | Operation::ImportIsoformPanel {
            panel_path: path, ..
        }
        | Operation::InterpretRnaReads {
            input_path: path, ..
        }
        | Operation::ExportRnaReadReport { path, .. }
        | Operation::ExportRnaReadHitsFasta { path, .. }
        | Operation::ExportRnaReadSampleSheet { path, .. }
        | Operation::ExportRnaReadExonPathsTsv { path, .. }
        | Operation::ExportRnaReadExonAbundanceTsv { path, .. }
        | Operation::ExportRnaReadScoreDensitySvg { path, .. }
        | Operation::ExportRnaReadAlignmentsTsv { path, .. }
        | Operation::ExportRnaReadIsoformTriageTsv { path, .. }
        | Operation::ExportRnaReadAlignmentDotplotSvg { path, .. }
        | Operation::RenderTfbsScoreTracksSvg { path, .. }
        | Operation::RenderProtocolCartoonSvg { path, .. }
        | Operation::RenderProtocolCartoonTemplateSvg { path, .. }
        | Operation::RenderProtocolCartoonTemplateWithBindingsSvg { path, .. }
        | Operation::ExportProtocolCartoonTemplateJson { path, .. } = op
        {
            Self::push_unique_token(&mut summary.file_paths, path);
        }
        if let Operation::PreflightRnaReadIsoforms {
            positive_transcript_fasta_paths,
            control_transcript_fasta_paths,
            ..
        } = op
        {
            for path in positive_transcript_fasta_paths {
                Self::push_unique_token(&mut summary.file_paths, path);
            }
            for path in control_transcript_fasta_paths {
                Self::push_unique_token(&mut summary.file_paths, path);
            }
        }
        if let Operation::SummarizeRnaReadGeneSupport {
            path: Some(path), ..
        } = op
        {
            Self::push_unique_token(&mut summary.file_paths, path);
        }
        if let Operation::InspectRnaReadGeneSupport {
            path: Some(path), ..
        } = op
        {
            Self::push_unique_token(&mut summary.file_paths, path);
        }
        if let Operation::FindRestrictionSites {
            path: Some(path), ..
        }
        | Operation::ScanTfbsHits {
            path: Some(path), ..
        }
        | Operation::QueryRepeatAnnotations {
            path: Some(path), ..
        }
        | Operation::BuildRepeatEnvironmentCohort {
            path: Some(path), ..
        }
        | Operation::SummarizeWindowCohortTfbs {
            path: Some(path), ..
        } = op
        {
            Self::push_unique_token(&mut summary.file_paths, path);
        }
        if let Operation::SummarizeJasparEntries {
            path: Some(path), ..
        } = op
        {
            Self::push_unique_token(&mut summary.file_paths, path);
        }
        summary
    }

    pub(super) fn collect_run_bundle_export_paths(op: &Operation) -> Vec<String> {
        let mut paths: Vec<String> = vec![];
        let mut push = |path: &str| Self::push_unique_token(&mut paths, path);
        match op {
            Operation::SaveFile { path, .. }
            | Operation::RenderSequenceSvg { path, .. }
            | Operation::RenderDotplotSvg { path, .. }
            | Operation::RenderTfbsScoreTracksSvg { path, .. }
            | Operation::RenderFeatureExpertSvg { path, .. }
            | Operation::RenderIsoformArchitectureSvg { path, .. }
            | Operation::RenderRnaStructureSvg { path, .. }
            | Operation::RenderLineageSvg { path }
            | Operation::RenderPoolGelSvg { path, .. }
            | Operation::RenderProteinGelSvg { path, .. }
            | Operation::RenderProteinGelReportsSvg { path, .. }
            | Operation::RenderProteaseDigestGelSvg { path, .. }
            | Operation::RenderProtein2dGelSvg { path, .. }
            | Operation::RenderProtocolCartoonSvg { path, .. }
            | Operation::RenderProtocolCartoonTemplateSvg { path, .. }
            | Operation::RenderProtocolCartoonTemplateWithBindingsSvg { path, .. }
            | Operation::ExportProtocolCartoonTemplateJson { path, .. }
            | Operation::ExportDnaLadders { path, .. }
            | Operation::ExportRnaLadders { path, .. }
            | Operation::ExportPool { path, .. }
            | Operation::ExportProcessRunBundle { path, .. }
            | Operation::ExportLabAssistantInstructions { path, .. }
            | Operation::ExportGuideOligos { path, .. }
            | Operation::ExportGuideProtocolText { path, .. }
            | Operation::ExportPrimerDesignReport { path, .. }
            | Operation::AssessPrimerPairSpecificity {
                path: Some(path), ..
            }
            | Operation::ExportRnaReadReport { path, .. }
            | Operation::ExportRnaReadHitsFasta { path, .. }
            | Operation::ExportRnaReadSampleSheet { path, .. }
            | Operation::ExportRnaReadExonPathsTsv { path, .. }
            | Operation::ExportRnaReadExonAbundanceTsv { path, .. }
            | Operation::ExportRnaReadScoreDensitySvg { path, .. }
            | Operation::ExportRnaReadAlignmentsTsv { path, .. }
            | Operation::ExportRnaReadIsoformTriageTsv { path, .. }
            | Operation::ExportRnaReadAlignmentDotplotSvg { path, .. } => {
                push(path);
            }
            Operation::ExportSequenceContextBundle { output_dir, .. } => push(output_dir),
            Operation::SummarizeRnaReadGeneSupport {
                path: Some(path), ..
            } => push(path),
            Operation::InspectRnaReadGeneSupport {
                path: Some(path), ..
            } => push(path),
            Operation::FindRestrictionSites {
                path: Some(path), ..
            } => push(path),
            Operation::ScanTfbsHits {
                path: Some(path), ..
            } => push(path),
            Operation::QueryRepeatAnnotations {
                path: Some(path), ..
            } => push(path),
            Operation::BuildRepeatEnvironmentCohort {
                path: Some(path), ..
            } => push(path),
            Operation::SummarizeWindowCohortTfbs {
                path: Some(path), ..
            } => push(path),
            Operation::SummarizeJasparEntries {
                path: Some(path), ..
            } => push(path),
            Operation::ListReporterCatalog {
                path: Some(path), ..
            }
            | Operation::RecommendReporters {
                path: Some(path), ..
            }
            | Operation::ExportReporterCorpus { path, .. }
            | Operation::PlanReporterConstructHandoff {
                path: Some(path), ..
            } => push(path),
            _ => {}
        }
        paths
    }

    fn normalize_routine_decision_trace_preflight_snapshot(
        snapshot: &mut RoutineDecisionTracePreflightSnapshot,
    ) {
        let mut warnings = vec![];
        for warning in std::mem::take(&mut snapshot.warnings) {
            let warning = warning.trim();
            if warning.is_empty() {
                continue;
            }
            Self::push_unique_token(&mut warnings, warning);
        }
        snapshot.warnings = warnings;

        let mut errors = vec![];
        for error in std::mem::take(&mut snapshot.errors) {
            let error = error.trim();
            if error.is_empty() {
                continue;
            }
            Self::push_unique_token(&mut errors, error);
        }
        snapshot.errors = errors;

        snapshot.contract_source = snapshot
            .contract_source
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
    }

    fn normalize_routine_decision_disambiguation_question_id_from_text(text: &str) -> String {
        let mut out = String::new();
        let mut last_was_sep = false;
        for ch in text.trim().chars().flat_map(|c| c.to_lowercase()) {
            if ch.is_ascii_alphanumeric() {
                out.push(ch);
                last_was_sep = false;
            } else if !last_was_sep {
                out.push('_');
                last_was_sep = true;
            }
            if out.len() >= 48 {
                break;
            }
        }
        let compact = out.trim_matches('_').to_string();
        if compact.is_empty() {
            "question".to_string()
        } else {
            compact
        }
    }

    fn normalize_routine_preference_context_record(
        mut context: RoutinePreferenceContextRecord,
    ) -> RoutinePreferenceContextRecord {
        context.helper_profile_id = context
            .helper_profile_id
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        context.helper_resolution_status = context.helper_resolution_status.trim().to_string();
        if context.helper_resolution_status.is_empty() {
            context.helper_resolution_status = "not_requested".to_string();
        }
        let normalize_vec = |values: &mut Vec<String>| {
            let mut normalized = vec![];
            for value in std::mem::take(values) {
                Self::push_unique_token(&mut normalized, &value);
            }
            *values = normalized;
        };
        normalize_vec(&mut context.explicit_preferred_routine_families);
        normalize_vec(&mut context.helper_derived_preferred_routine_families);
        normalize_vec(&mut context.effective_preferred_routine_families);
        normalize_vec(&mut context.helper_offered_functions);
        normalize_vec(&mut context.helper_component_labels);
        normalize_vec(&mut context.rationale);
        context
    }

    fn normalize_routine_decision_trace_candidate_score(
        mut score: RoutineDecisionTraceCandidateScore,
    ) -> Option<RoutineDecisionTraceCandidateScore> {
        score.routine_id = score.routine_id.trim().to_string();
        if score.routine_id.is_empty() {
            return None;
        }
        score.routine_title = score
            .routine_title
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        score.routine_family = score.routine_family.trim().to_string();
        score.estimated_time_hours = score
            .estimated_time_hours
            .filter(|value| value.is_finite() && *value >= 0.0);
        score.estimated_cost = score
            .estimated_cost
            .filter(|value| value.is_finite() && *value >= 0.0);
        score.local_fit_score = score
            .local_fit_score
            .filter(|value| value.is_finite() && *value >= 0.0 && *value <= 1.0);
        score.composite_meta_score = score.composite_meta_score.filter(|value| value.is_finite());
        score.routine_family_alignment_bonus = score
            .routine_family_alignment_bonus
            .filter(|value| value.is_finite());
        let mut sources = vec![];
        for source in std::mem::take(&mut score.routine_family_alignment_sources) {
            Self::push_unique_token(&mut sources, &source);
        }
        score.routine_family_alignment_sources = sources;
        Some(score)
    }

    fn normalize_macro_template_suggestion(
        mut suggestion: MacroTemplateSuggestion,
    ) -> Option<MacroTemplateSuggestion> {
        suggestion.macro_kind = suggestion.macro_kind.trim().to_string();
        if suggestion.macro_kind.is_empty() {
            return None;
        }
        suggestion.template_name = suggestion.template_name.trim().to_string();
        if suggestion.template_name.is_empty() {
            return None;
        }
        suggestion.description = suggestion
            .description
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        suggestion.details_url = suggestion
            .details_url
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        if !suggestion.score.is_finite() || suggestion.score < 0.0 {
            suggestion.score = 0.0;
        }
        let normalize_vec = |values: &mut Vec<String>| {
            let mut normalized = vec![];
            for value in std::mem::take(values) {
                Self::push_unique_token(&mut normalized, &value);
            }
            *values = normalized;
        };
        normalize_vec(&mut suggestion.matched_routine_families);
        normalize_vec(&mut suggestion.matched_terms);
        normalize_vec(&mut suggestion.rationale);
        Some(suggestion)
    }

    fn normalize_routine_decision_trace(
        mut trace: RoutineDecisionTrace,
    ) -> Option<RoutineDecisionTrace> {
        let schema = trace.schema.trim();
        if schema.is_empty() {
            trace.schema = ROUTINE_DECISION_TRACE_SCHEMA.to_string();
        } else if !schema.eq_ignore_ascii_case(ROUTINE_DECISION_TRACE_SCHEMA) {
            return None;
        } else {
            trace.schema = ROUTINE_DECISION_TRACE_SCHEMA.to_string();
        }

        trace.trace_id = trace.trace_id.trim().to_string();
        if trace.trace_id.is_empty() {
            return None;
        }
        trace.source = trace.source.trim().to_string();
        if trace.source.is_empty() {
            trace.source = "unknown".to_string();
        }
        trace.status = trace.status.trim().to_string();
        if trace.status.is_empty() {
            trace.status = "draft".to_string();
        }
        trace.goal_text = trace.goal_text.trim().to_string();
        trace.query_text = trace.query_text.trim().to_string();
        if trace.created_at_unix_ms == 0 {
            trace.created_at_unix_ms = trace.updated_at_unix_ms;
        }
        if trace.updated_at_unix_ms == 0 {
            trace.updated_at_unix_ms = trace.created_at_unix_ms;
        }

        let normalize_opt = |value: &mut Option<String>| {
            *value = value
                .take()
                .map(|v| v.trim().to_string())
                .filter(|v| !v.is_empty());
        };
        normalize_opt(&mut trace.selected_routine_id);
        normalize_opt(&mut trace.selected_routine_title);
        normalize_opt(&mut trace.selected_routine_family);
        normalize_opt(&mut trace.macro_instance_id);
        normalize_opt(&mut trace.execution_error);
        trace.routine_preference_context = trace
            .routine_preference_context
            .take()
            .map(Self::normalize_routine_preference_context_record);

        let mut dedup_candidates = vec![];
        for token in std::mem::take(&mut trace.candidate_routine_ids) {
            Self::push_unique_token(&mut dedup_candidates, &token);
        }
        trace.candidate_routine_ids = dedup_candidates;

        let mut candidate_scores: Vec<RoutineDecisionTraceCandidateScore> = vec![];
        let mut seen_candidate_ids: HashSet<String> = HashSet::new();
        for row in std::mem::take(&mut trace.candidate_planning_scores) {
            let Some(row) = Self::normalize_routine_decision_trace_candidate_score(row) else {
                continue;
            };
            if !seen_candidate_ids.insert(row.routine_id.to_ascii_lowercase()) {
                continue;
            }
            candidate_scores.push(row);
        }
        candidate_scores.sort_by(|left, right| {
            right
                .passes_guardrails
                .cmp(&left.passes_guardrails)
                .then_with(|| {
                    right
                        .composite_meta_score
                        .unwrap_or(f64::NEG_INFINITY)
                        .total_cmp(&left.composite_meta_score.unwrap_or(f64::NEG_INFINITY))
                })
                .then_with(|| left.routine_family.cmp(&right.routine_family))
                .then_with(|| left.routine_id.cmp(&right.routine_id))
        });
        trace.candidate_planning_scores = candidate_scores;

        let mut dedup_alternatives = vec![];
        for token in std::mem::take(&mut trace.alternatives_presented) {
            Self::push_unique_token(&mut dedup_alternatives, &token);
        }
        trace.alternatives_presented = dedup_alternatives;

        let mut macro_suggestions: Vec<MacroTemplateSuggestion> = vec![];
        let mut seen_macro_keys: HashSet<String> = HashSet::new();
        for row in std::mem::take(&mut trace.macro_suggestions) {
            let Some(row) = Self::normalize_macro_template_suggestion(row) else {
                continue;
            };
            let key = format!(
                "{}\u{1f}{}",
                row.macro_kind.to_ascii_lowercase(),
                row.template_name.to_ascii_lowercase()
            );
            if !seen_macro_keys.insert(key) {
                continue;
            }
            macro_suggestions.push(row);
        }
        macro_suggestions.sort_by(|left, right| {
            right
                .score
                .total_cmp(&left.score)
                .then_with(|| left.macro_kind.cmp(&right.macro_kind))
                .then_with(|| left.template_name.cmp(&right.template_name))
        });
        trace.macro_suggestions = macro_suggestions;

        let mut normalized_questions: Vec<RoutineDecisionTraceDisambiguationQuestion> = vec![];
        let mut used_question_ids: HashMap<String, usize> = HashMap::new();
        for mut row in std::mem::take(&mut trace.disambiguation_questions_presented) {
            row.question_id = row.question_id.trim().to_string();
            row.question_text = row.question_text.trim().to_string();
            if row.question_text.is_empty() {
                continue;
            }
            let base_id = if row.question_id.is_empty() {
                Self::normalize_routine_decision_disambiguation_question_id_from_text(
                    &row.question_text,
                )
            } else {
                row.question_id.clone()
            };
            let count = used_question_ids.entry(base_id.clone()).or_insert(0);
            *count += 1;
            row.question_id = if *count == 1 {
                base_id
            } else {
                format!("{}_{}", base_id, *count)
            };
            if normalized_questions.iter().any(|existing| {
                existing.question_id.eq_ignore_ascii_case(&row.question_id)
                    || existing
                        .question_text
                        .eq_ignore_ascii_case(&row.question_text)
            }) {
                continue;
            }
            normalized_questions.push(row);
        }
        trace.disambiguation_questions_presented = normalized_questions;

        let mut normalized_answers_by_question: BTreeMap<String, String> = BTreeMap::new();
        for mut row in std::mem::take(&mut trace.disambiguation_answers) {
            row.question_id = row.question_id.trim().to_string();
            row.answer_text = row.answer_text.trim().to_string();
            if row.question_id.is_empty() || row.answer_text.is_empty() {
                continue;
            }
            normalized_answers_by_question.insert(row.question_id, row.answer_text);
        }
        trace.disambiguation_answers = normalized_answers_by_question
            .into_iter()
            .map(
                |(question_id, answer_text)| RoutineDecisionTraceDisambiguationAnswer {
                    question_id,
                    answer_text,
                },
            )
            .collect();

        let mut dedup_emitted_op_ids = vec![];
        for token in std::mem::take(&mut trace.emitted_operation_ids) {
            Self::push_unique_token(&mut dedup_emitted_op_ids, &token);
        }
        trace.emitted_operation_ids = dedup_emitted_op_ids;

        let mut comparisons_out: Vec<RoutineDecisionTraceComparison> = vec![];
        let mut seen_comparisons: HashSet<String> = HashSet::new();
        for mut row in std::mem::take(&mut trace.comparisons) {
            row.left_routine_id = row.left_routine_id.trim().to_string();
            row.right_routine_id = row.right_routine_id.trim().to_string();
            if row.left_routine_id.is_empty() || row.right_routine_id.is_empty() {
                continue;
            }
            let key = format!("{}\u{1f}{}", row.left_routine_id, row.right_routine_id);
            if !seen_comparisons.insert(key) {
                continue;
            }
            comparisons_out.push(row);
        }
        trace.comparisons = comparisons_out;

        let mut normalized_bindings: BTreeMap<String, String> = BTreeMap::new();
        for (key, value) in std::mem::take(&mut trace.bindings_snapshot) {
            let key = key.trim().to_string();
            let value = value.trim().to_string();
            if key.is_empty() || value.is_empty() {
                continue;
            }
            normalized_bindings.insert(key, value);
        }
        trace.bindings_snapshot = normalized_bindings;

        let mut preflight_history: Vec<RoutineDecisionTracePreflightSnapshot> = vec![];
        for mut snapshot in std::mem::take(&mut trace.preflight_history) {
            Self::normalize_routine_decision_trace_preflight_snapshot(&mut snapshot);
            preflight_history.push(snapshot);
        }
        trace.preflight_history = preflight_history;

        if let Some(snapshot) = trace.preflight_snapshot.as_mut() {
            Self::normalize_routine_decision_trace_preflight_snapshot(snapshot);
        }
        if trace.preflight_history.is_empty()
            && let Some(snapshot) = trace.preflight_snapshot.clone()
        {
            trace.preflight_history.push(snapshot);
        }
        trace.preflight_snapshot = trace.preflight_history.last().cloned();

        let mut export_events_out: Vec<RoutineDecisionTraceExportEvent> = vec![];
        let mut seen_export_events: HashSet<String> = HashSet::new();
        for mut event in std::mem::take(&mut trace.export_events) {
            event.run_bundle_path = event.run_bundle_path.trim().to_string();
            if event.run_bundle_path.is_empty() {
                continue;
            }
            let key = format!(
                "{}\u{1f}{}",
                event.exported_at_unix_ms, event.run_bundle_path
            );
            if !seen_export_events.insert(key) {
                continue;
            }
            export_events_out.push(event);
        }
        export_events_out.sort_by(|left, right| {
            left.exported_at_unix_ms
                .cmp(&right.exported_at_unix_ms)
                .then_with(|| left.run_bundle_path.cmp(&right.run_bundle_path))
        });
        trace.export_events = export_events_out;
        Some(trace)
    }

    fn normalize_routine_decision_traces(
        traces: Vec<RoutineDecisionTrace>,
    ) -> Vec<RoutineDecisionTrace> {
        let mut by_trace_id: HashMap<String, RoutineDecisionTrace> = HashMap::new();
        for trace in traces {
            let Some(normalized) = Self::normalize_routine_decision_trace(trace) else {
                continue;
            };
            let replace_existing = by_trace_id
                .get(&normalized.trace_id)
                .map(|existing| {
                    (
                        normalized.updated_at_unix_ms,
                        normalized.created_at_unix_ms,
                        normalized.trace_id.as_str(),
                    ) > (
                        existing.updated_at_unix_ms,
                        existing.created_at_unix_ms,
                        existing.trace_id.as_str(),
                    )
                })
                .unwrap_or(true);
            if replace_existing {
                by_trace_id.insert(normalized.trace_id.clone(), normalized);
            }
        }
        let mut out = by_trace_id.into_values().collect::<Vec<_>>();
        out.sort_by(|left, right| {
            left.created_at_unix_ms
                .cmp(&right.created_at_unix_ms)
                .then_with(|| left.trace_id.cmp(&right.trace_id))
        });
        out
    }

    fn read_routine_decision_traces_from_metadata(&self) -> Vec<RoutineDecisionTrace> {
        let Some(raw) = self
            .state
            .metadata
            .get(ROUTINE_DECISION_TRACES_METADATA_KEY)
            .cloned()
        else {
            return vec![];
        };
        if let Ok(mut store) = serde_json::from_value::<RoutineDecisionTraceStore>(raw.clone()) {
            if store.schema.trim().is_empty() {
                store.schema = ROUTINE_DECISION_TRACE_STORE_SCHEMA.to_string();
            }
            if !store
                .schema
                .trim()
                .eq_ignore_ascii_case(ROUTINE_DECISION_TRACE_STORE_SCHEMA)
            {
                return vec![];
            }
            return Self::normalize_routine_decision_traces(store.traces);
        }
        serde_json::from_value::<Vec<RoutineDecisionTrace>>(raw)
            .map(Self::normalize_routine_decision_traces)
            .unwrap_or_default()
    }

    pub(super) fn export_process_run_bundle_file(
        &self,
        path: &str,
        run_id_filter: Option<&str>,
    ) -> Result<ProcessRunBundleExport, EngineError> {
        let path = path.trim();
        if path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ExportProcessRunBundle requires non-empty path".to_string(),

                cause_chain: vec![],
            });
        }
        let normalized_run_id = run_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty());
        let selected_records: Vec<OperationRecord> = self
            .journal
            .iter()
            .filter(|record| {
                normalized_run_id
                    .map(|run_id| record.run_id == run_id)
                    .unwrap_or(true)
            })
            .cloned()
            .collect();
        if normalized_run_id.is_some() && selected_records.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No operation records found for run_id '{}'",
                    normalized_run_id.unwrap_or_default()
                ),

                cause_chain: vec![],
            });
        }

        let mut created_seq_ids: BTreeSet<String> = BTreeSet::new();
        let mut changed_seq_ids: BTreeSet<String> = BTreeSet::new();
        let mut op_ids: BTreeSet<String> = BTreeSet::new();
        for record in &selected_records {
            op_ids.insert(record.result.op_id.clone());
            for seq_id in &record.result.created_seq_ids {
                created_seq_ids.insert(seq_id.clone());
            }
            for seq_id in &record.result.changed_seq_ids {
                changed_seq_ids.insert(seq_id.clone());
            }
        }

        let operation_inputs: Vec<ProcessRunBundleOperationInputSummary> = selected_records
            .iter()
            .enumerate()
            .map(|(index, record)| {
                Self::summarize_run_bundle_operation_inputs(
                    &record.op,
                    &record.result.op_id,
                    &record.run_id,
                    index + 1,
                )
            })
            .collect();
        let mut referenced_sequence_ids: BTreeSet<String> = BTreeSet::new();
        let mut referenced_container_ids: BTreeSet<String> = BTreeSet::new();
        let mut referenced_arrangement_ids: BTreeSet<String> = BTreeSet::new();
        let mut referenced_candidate_set_ids: BTreeSet<String> = BTreeSet::new();
        let mut referenced_guide_set_ids: BTreeSet<String> = BTreeSet::new();
        let mut referenced_genome_ids: BTreeSet<String> = BTreeSet::new();
        let mut file_inputs: BTreeSet<String> = BTreeSet::new();
        for summary in &operation_inputs {
            referenced_sequence_ids.extend(summary.sequence_ids.iter().cloned());
            referenced_container_ids.extend(summary.container_ids.iter().cloned());
            referenced_arrangement_ids.extend(summary.arrangement_ids.iter().cloned());
            referenced_candidate_set_ids.extend(summary.candidate_set_ids.iter().cloned());
            referenced_guide_set_ids.extend(summary.guide_set_ids.iter().cloned());
            referenced_genome_ids.extend(summary.genome_ids.iter().cloned());
            file_inputs.extend(summary.file_paths.iter().cloned());
        }
        let root_sequence_ids = referenced_sequence_ids
            .iter()
            .filter(|seq_id| !created_seq_ids.contains(*seq_id))
            .cloned()
            .collect::<Vec<_>>();

        let mut parameter_overrides: Vec<ProcessRunBundleParameterOverride> = vec![];
        for (index, record) in selected_records.iter().enumerate() {
            if let Operation::SetParameter { name, value } = &record.op {
                parameter_overrides.push(ProcessRunBundleParameterOverride {
                    op_id: record.result.op_id.clone(),
                    run_id: record.run_id.clone(),
                    record_index: index + 1,
                    name: name.clone(),
                    value: value.clone(),
                });
            }
        }

        let mut created_container_ids: Vec<String> = self
            .state
            .container_state
            .containers
            .iter()
            .filter_map(|(container_id, container)| {
                container
                    .created_by_op
                    .as_deref()
                    .filter(|op_id| op_ids.contains(*op_id))
                    .map(|_| container_id.clone())
            })
            .collect();
        created_container_ids.sort_by_key(|value| value.to_ascii_lowercase());

        let mut created_arrangement_ids: Vec<String> = self
            .state
            .container_state
            .arrangements
            .iter()
            .filter_map(|(arrangement_id, arrangement)| {
                arrangement
                    .created_by_op
                    .as_deref()
                    .filter(|op_id| op_ids.contains(*op_id))
                    .map(|_| arrangement_id.clone())
            })
            .collect();
        created_arrangement_ids.sort_by_key(|value| value.to_ascii_lowercase());

        let mut exported_paths: BTreeSet<String> = BTreeSet::new();
        for record in &selected_records {
            for path in Self::collect_run_bundle_export_paths(&record.op) {
                exported_paths.insert(path);
            }
        }

        let final_sequence_ids: BTreeSet<String> =
            created_seq_ids.union(&changed_seq_ids).cloned().collect();
        let construct_reasoning_seq_ids: BTreeSet<String> = referenced_sequence_ids
            .iter()
            .cloned()
            .chain(final_sequence_ids.iter().cloned())
            .collect();
        let mut final_sequences: Vec<EngineSequenceSummary> = vec![];
        for seq_id in final_sequence_ids {
            let Some(dna) = self.state.sequences.get(&seq_id) else {
                continue;
            };
            final_sequences.push(EngineSequenceSummary {
                id: seq_id,
                name: dna.name().clone(),
                length: dna.len(),
                circular: dna.is_circular(),
            });
        }
        final_sequences.sort_by(|left, right| {
            left.id
                .to_ascii_lowercase()
                .cmp(&right.id.to_ascii_lowercase())
                .then_with(|| left.id.cmp(&right.id))
        });

        let decision_traces = self.read_routine_decision_traces_from_metadata();
        let construct_reasoning =
            self.collect_process_run_bundle_construct_reasoning(&construct_reasoning_seq_ids);

        let bundle = ProcessRunBundleExport {
            schema: PROCESS_RUN_BUNDLE_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            run_id_filter: normalized_run_id.map(|value| value.to_string()),
            selected_record_count: selected_records.len(),
            inputs: ProcessRunBundleInputs {
                root_sequence_ids,
                referenced_sequence_ids: referenced_sequence_ids.into_iter().collect(),
                referenced_container_ids: referenced_container_ids.into_iter().collect(),
                referenced_arrangement_ids: referenced_arrangement_ids.into_iter().collect(),
                referenced_candidate_set_ids: referenced_candidate_set_ids.into_iter().collect(),
                referenced_guide_set_ids: referenced_guide_set_ids.into_iter().collect(),
                referenced_genome_ids: referenced_genome_ids.into_iter().collect(),
                file_inputs: file_inputs.into_iter().collect(),
                operation_inputs,
            },
            parameter_overrides,
            decision_traces,
            operation_log: selected_records,
            outputs: ProcessRunBundleOutputs {
                created_seq_ids: created_seq_ids.into_iter().collect(),
                changed_seq_ids: changed_seq_ids.into_iter().collect(),
                final_sequences,
                created_container_ids,
                created_arrangement_ids,
                exported_paths: exported_paths.into_iter().collect(),
            },
            parameter_snapshot: serde_json::to_value(&self.state.parameters).map_err(|e| {
                EngineError {
                    code: ErrorCode::Internal,
                    message: format!("Could not serialize engine parameter snapshot: {e}"),

                    cause_chain: vec![],
                }
            })?,
            construct_reasoning,
        };

        let text = serde_json::to_string_pretty(&bundle).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize run bundle JSON: {e}"),

            cause_chain: vec![],
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write run bundle file '{path}': {e}"),

            cause_chain: vec![],
        })?;
        Ok(bundle)
    }

    pub(super) fn export_lab_assistant_instructions_file(
        &self,
        path: &str,
        run_id_filter: Option<&str>,
        title: Option<&str>,
        audience: Option<&str>,
        format: Option<LabAssistantInstructionsFormat>,
    ) -> Result<LabAssistantInstructionsExport, EngineError> {
        let path = path.trim();
        if path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ExportLabAssistantInstructions requires non-empty path".to_string(),

                cause_chain: vec![],
            });
        }
        let output_format =
            format.unwrap_or_else(|| LabAssistantInstructionsFormat::infer_from_path(path));
        let normalized_run_id = run_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty());
        let selected_records: Vec<&OperationRecord> = self
            .journal
            .iter()
            .filter(|record| {
                normalized_run_id
                    .map(|run_id| record.run_id == run_id)
                    .unwrap_or(true)
            })
            .collect();
        if normalized_run_id.is_some() && selected_records.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No operation records found for run_id '{}'",
                    normalized_run_id.unwrap_or_default()
                ),

                cause_chain: vec![],
            });
        }

        let op_ids = selected_records
            .iter()
            .map(|record| record.result.op_id.clone())
            .collect::<BTreeSet<_>>();
        let title = title
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .or_else(|| Self::infer_lab_assistant_title(&selected_records))
            .unwrap_or_else(|| "GENtle cloning handoff".to_string());
        let audience = audience
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("Lab assistant")
            .to_string();

        let material_rows = self.lab_assistant_material_rows(&selected_records, &op_ids);
        let mut warning_lines = vec![];
        let embedded_images = if matches!(
            output_format,
            LabAssistantInstructionsFormat::Odt | LabAssistantInstructionsFormat::Docx
        ) {
            match self.lab_assistant_default_embedded_images() {
                Ok(images) => images,
                Err(err) => {
                    warning_lines.push(format!(
                        "Could not embed lineage overview graphic: {}",
                        err.message
                    ));
                    Vec::new()
                }
            }
        } else {
            Vec::new()
        };
        let embedded_visuals = embedded_images
            .iter()
            .map(|image| image.visual_row("operation lineage overview"))
            .collect::<Vec<_>>();
        if selected_records.is_empty() {
            warning_lines.push(
                "No recorded operations were available; this export is a general handoff scaffold."
                    .to_string(),
            );
        }
        let step_sections = Self::lab_assistant_instruction_sections(&selected_records);
        let no_cloning_specific_steps = step_sections
            .iter()
            .find(|section| section.heading == "Design-derived bench sequence")
            .is_some_and(|section| {
                section
                    .steps
                    .iter()
                    .any(|step| step.contains("No cloning-specific operations"))
            });
        if no_cloning_specific_steps {
            warning_lines.push(
                "No cloning-specific operation was recognized in the selected run; review the run bundle or project history before bench work."
                    .to_string(),
            );
        }

        let checkpoint_lines = vec![
            "Confirm every physical tube label matches the sequence/container IDs listed below."
                .to_string(),
            "Verify expected product lengths against the GENtle output and any gel/band summaries before proceeding."
                .to_string(),
            "If sequencing, restriction digest, or colony PCR confirmation is planned, keep sample IDs synchronized with this handoff."
                .to_string(),
        ];
        let safety_lines = vec![
            "Use institution-approved SOPs, kit manuals, and supervisor-approved conditions for all wet-lab execution."
                .to_string(),
            "GENtle exports design-derived intent and checkpoints; it does not replace local biosafety, chemical safety, or strain/vector approval."
                .to_string(),
            "Do not proceed when material labels, expected lengths, antibiotic/host context, or approval status are ambiguous."
                .to_string(),
        ];
        let record_keeping_lines = vec![
            "Record operator, date, reagent lot numbers, template/vector/insert tube IDs, and plate/tube coordinates."
                .to_string(),
            "Attach photos or instrument exports for gels, plates, colony PCR, digest checks, and sequencing chromatograms."
                .to_string(),
            "Store deviations from this handoff together with the GENtle project state or run bundle."
                .to_string(),
        ];
        let summary_lines = vec![
            format!(
                "Selected {} operation record(s){}.",
                selected_records.len(),
                normalized_run_id
                    .map(|run_id| format!(" for run_id `{run_id}`"))
                    .unwrap_or_default()
            ),
            format!(
                "Prepared {} material row(s) and {} instruction section(s).",
                material_rows.len(),
                step_sections.len()
            ),
            format!(
                "Embedded {} graphical overview(s) in {} output.",
                embedded_visuals.len(),
                output_format.as_str()
            ),
        ];

        let export = LabAssistantInstructionsExport {
            schema: LAB_ASSISTANT_INSTRUCTIONS_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            title,
            audience,
            output_path: path.to_string(),
            output_format,
            run_id_filter: normalized_run_id.map(str::to_string),
            selected_record_count: selected_records.len(),
            material_rows,
            step_sections,
            embedded_visuals,
            checkpoint_lines,
            safety_lines,
            record_keeping_lines,
            warning_lines,
            summary_lines,
        };

        let parent = Path::new(path)
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        std::fs::create_dir_all(parent).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create parent directory for lab assistant instructions '{path}': {e}"
            ),

            cause_chain: vec![],
        })?;
        let bytes = match output_format {
            LabAssistantInstructionsFormat::Markdown => {
                crate::engine::lab_assistant_export::render_markdown(&export).into_bytes()
            }
            LabAssistantInstructionsFormat::Odt => {
                crate::engine::lab_assistant_export::render_odt(&export, &embedded_images).map_err(
                    |e| EngineError {
                        code: ErrorCode::Internal,
                        message: format!("Could not render ODT lab assistant report: {e}"),

                        cause_chain: vec![],
                    },
                )?
            }
            LabAssistantInstructionsFormat::Docx => {
                crate::engine::lab_assistant_export::render_docx(&export, &embedded_images)
                    .map_err(|e| EngineError {
                        code: ErrorCode::Internal,
                        message: format!("Could not render DOCX lab assistant report: {e}"),

                        cause_chain: vec![],
                    })?
            }
        };
        std::fs::write(path, bytes).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write lab assistant instructions file '{path}': {e}"),

            cause_chain: vec![],
        })?;
        Ok(export)
    }

    fn lab_assistant_default_embedded_images(
        &self,
    ) -> Result<Vec<crate::engine::lab_assistant_export::LabAssistantEmbeddedImage>, EngineError>
    {
        let svg = crate::lineage_export::export_lineage_svg(&self.state, &self.journal);
        let rendered = crate::svg_png::render_svg_to_png_bytes(
            &svg,
            crate::svg_png::SvgPngRenderOptions {
                scale: 1.5,
                drop_dotplot_metadata: false,
            },
        )
        .map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not rasterize lineage overview SVG: {e}"),

            cause_chain: vec![],
        })?;
        Ok(vec![
            crate::engine::lab_assistant_export::LabAssistantEmbeddedImage {
                visual_id: "lineage_overview".to_string(),
                label: "GENtle operation lineage overview".to_string(),
                file_name: "lineage_overview.png".to_string(),
                bytes: rendered.bytes,
                width_px: rendered.width as usize,
                height_px: rendered.height as usize,
            },
        ])
    }

    fn infer_lab_assistant_title(records: &[&OperationRecord]) -> Option<String> {
        records.iter().find_map(|record| match &record.op {
            Operation::ApplyGibsonAssemblyPlan { plan_json } => {
                let plan: crate::gibson_planning::GibsonAssemblyPlan =
                    serde_json::from_str(plan_json).ok()?;
                let title = plan.title.trim();
                (!title.is_empty()).then(|| format!("Lab handoff: {title}"))
            }
            Operation::PrepareRestrictionCloningPcrHandoff { .. } => {
                Some("Lab handoff: restriction-cloning PCR".to_string())
            }
            Operation::Pcr { .. }
            | Operation::PcrAdvanced { .. }
            | Operation::PcrMutagenesis { .. } => {
                Some("Lab handoff: PCR product preparation".to_string())
            }
            _ => None,
        })
    }

    fn lab_assistant_material_rows(
        &self,
        records: &[&OperationRecord],
        op_ids: &BTreeSet<String>,
    ) -> Vec<LabAssistantMaterialRow> {
        let mut referenced_sequence_ids = BTreeSet::new();
        let mut created_sequence_ids = BTreeSet::new();
        let mut referenced_container_ids = BTreeSet::new();
        let mut referenced_arrangement_ids = BTreeSet::new();
        for (index, record) in records.iter().enumerate() {
            let summary = Self::summarize_run_bundle_operation_inputs(
                &record.op,
                &record.result.op_id,
                &record.run_id,
                index + 1,
            );
            referenced_sequence_ids.extend(summary.sequence_ids);
            referenced_container_ids.extend(summary.container_ids);
            referenced_arrangement_ids.extend(summary.arrangement_ids);
            created_sequence_ids.extend(record.result.created_seq_ids.iter().cloned());
        }

        let root_sequence_ids = referenced_sequence_ids
            .difference(&created_sequence_ids)
            .cloned()
            .collect::<Vec<_>>();
        let mut rows = vec![];
        for seq_id in root_sequence_ids {
            if let Some(row) = self.lab_assistant_sequence_material_row(&seq_id, "input_sequence") {
                rows.push(row);
            }
        }
        for seq_id in created_sequence_ids {
            if let Some(row) = self.lab_assistant_sequence_material_row(&seq_id, "designed_output")
            {
                rows.push(row);
            }
        }

        let created_container_ids =
            self.state
                .container_state
                .containers
                .iter()
                .filter_map(|(container_id, container)| {
                    container
                        .created_by_op
                        .as_deref()
                        .filter(|op_id| op_ids.contains(*op_id))
                        .map(|_| container_id.clone())
                });
        referenced_container_ids.extend(created_container_ids);
        for container_id in referenced_container_ids {
            let Some(container) = self.state.container_state.containers.get(&container_id) else {
                continue;
            };
            rows.push(LabAssistantMaterialRow {
                material_id: container_id.clone(),
                display_name: container
                    .name
                    .clone()
                    .unwrap_or_else(|| container_id.clone()),
                kind: format!("container_{:?}", container.kind).to_ascii_lowercase(),
                source: "container_state".to_string(),
                length_bp: None,
                topology: None,
                members: container.members.clone(),
                notes: vec![format!(
                    "Declared contents are {}.",
                    if container.declared_contents_exclusive {
                        "exclusive"
                    } else {
                        "not exclusive"
                    }
                )],
            });
        }

        let created_arrangement_ids = self.state.container_state.arrangements.iter().filter_map(
            |(arrangement_id, arrangement)| {
                arrangement
                    .created_by_op
                    .as_deref()
                    .filter(|op_id| op_ids.contains(*op_id))
                    .map(|_| arrangement_id.clone())
            },
        );
        referenced_arrangement_ids.extend(created_arrangement_ids);
        for arrangement_id in referenced_arrangement_ids {
            let Some(arrangement) = self.state.container_state.arrangements.get(&arrangement_id)
            else {
                continue;
            };
            rows.push(LabAssistantMaterialRow {
                material_id: arrangement_id.clone(),
                display_name: arrangement
                    .name
                    .clone()
                    .unwrap_or_else(|| arrangement_id.clone()),
                kind: "arrangement".to_string(),
                source: "container_state".to_string(),
                length_bp: None,
                topology: Some(format!("{:?}", arrangement.mode).to_ascii_lowercase()),
                members: arrangement.lane_container_ids.clone(),
                notes: arrangement
                    .ladders
                    .iter()
                    .map(|ladder| format!("Gel ladder: {ladder}"))
                    .collect(),
            });
        }

        rows.sort_by(|left, right| {
            left.kind
                .cmp(&right.kind)
                .then(left.material_id.cmp(&right.material_id))
        });
        rows
    }

    fn lab_assistant_sequence_material_row(
        &self,
        seq_id: &str,
        source: &str,
    ) -> Option<LabAssistantMaterialRow> {
        let dna = self.state.sequences.get(seq_id)?;
        Some(LabAssistantMaterialRow {
            material_id: seq_id.to_string(),
            display_name: dna.name().clone().unwrap_or_else(|| seq_id.to_string()),
            kind: "sequence".to_string(),
            source: source.to_string(),
            length_bp: Some(dna.len()),
            topology: Some(if dna.is_circular() {
                "circular".to_string()
            } else {
                "linear".to_string()
            }),
            members: vec![],
            notes: vec![],
        })
    }

    fn lab_assistant_instruction_sections(
        records: &[&OperationRecord],
    ) -> Vec<LabAssistantInstructionSection> {
        let mut prepare_steps = vec![
            "Read this handoff completely before touching samples; resolve any unclear material IDs with the designing scientist."
                .to_string(),
            "Gather the listed DNA samples, primers/oligos, vectors, enzymes/kits, competent cells, media, antibiotics, and confirmation reagents required by the local SOP."
                .to_string(),
            "Create tube/plate labels that exactly match the material IDs and planned output IDs.".to_string(),
        ];
        if records.is_empty() {
            prepare_steps.push(
                "No recorded design operations were selected; use this export only as a checklist scaffold."
                    .to_string(),
            );
        }

        let mut bench_steps = vec![];
        for record in records {
            bench_steps.extend(Self::lab_assistant_steps_for_record(record));
        }
        if bench_steps.is_empty() {
            bench_steps.push(
                "No cloning-specific operations were recognized in the selected run; ask the designing scientist for the intended bench route before proceeding."
                    .to_string(),
            );
        }

        vec![
            LabAssistantInstructionSection {
                heading: "Before starting".to_string(),
                steps: prepare_steps,
            },
            LabAssistantInstructionSection {
                heading: "Design-derived bench sequence".to_string(),
                steps: bench_steps,
            },
            LabAssistantInstructionSection {
                heading: "After the bench work".to_string(),
                steps: vec![
                    "Compare observed bands, colonies, digests, or sequencing reads with the expected GENtle outputs before declaring success."
                        .to_string(),
                    "Keep unsuccessful or ambiguous outcomes; they are useful for redesign and troubleshooting."
                        .to_string(),
                ],
            },
        ]
    }

    fn lab_assistant_steps_for_record(record: &OperationRecord) -> Vec<String> {
        match &record.op {
            Operation::LoadFile { path, as_id } => vec![format!(
                "Use design input `{}` from `{path}`.",
                as_id.as_deref().unwrap_or("derived sequence ID")
            )],
            Operation::CreateSequenceFromText {
                output_id,
                name,
                circular,
                ..
            } => vec![format!(
                "Use inline {} sequence `{}`{} as a project input.",
                if *circular { "circular" } else { "linear" },
                output_id
                    .as_deref()
                    .or(name.as_deref())
                    .unwrap_or("derived sequence ID"),
                name.as_deref()
                    .map(|value| format!(" ({value})"))
                    .unwrap_or_default()
            )],
            Operation::Digest {
                input,
                enzymes,
                output_prefix,
            } => vec![format!(
                "Digest `{input}` with {} and keep fragments under prefix `{}`.",
                Self::lab_join_or_dash(enzymes),
                output_prefix.as_deref().unwrap_or("GENtle-generated")
            )],
            Operation::DigestContainer {
                container_id,
                enzymes,
                output_prefix,
            } => vec![format!(
                "Digest all contents of container `{container_id}` with {} and keep products under prefix `{}`.",
                Self::lab_join_or_dash(enzymes),
                output_prefix.as_deref().unwrap_or("GENtle-generated")
            )],
            Operation::Pcr {
                template,
                forward_primer,
                reverse_primer,
                output_id,
                ..
            } => vec![format!(
                "Run PCR from template `{template}` with forward primer `{}` and reverse primer `{}`; expected product ID: `{}`.",
                Self::lab_oligo_display(forward_primer),
                Self::lab_oligo_display(reverse_primer),
                Self::lab_expected_output(output_id.as_deref(), &record.result.created_seq_ids)
            )],
            Operation::PcrAdvanced {
                template,
                forward_primer,
                reverse_primer,
                output_id,
                ..
            } => vec![format!(
                "Run constrained PCR from template `{template}` with forward primer `{}` and reverse primer `{}`; expected product ID: `{}`.",
                Self::lab_oligo_display(&forward_primer.sequence),
                Self::lab_oligo_display(&reverse_primer.sequence),
                Self::lab_expected_output(output_id.as_deref(), &record.result.created_seq_ids)
            )],
            Operation::PcrMutagenesis {
                template,
                forward_primer,
                reverse_primer,
                mutations,
                output_id,
                ..
            } => vec![format!(
                "Run mutagenic PCR from template `{template}` with {} requested edit(s), forward primer `{}`, reverse primer `{}`, and expected product ID `{}`.",
                mutations.len(),
                Self::lab_oligo_display(&forward_primer.sequence),
                Self::lab_oligo_display(&reverse_primer.sequence),
                Self::lab_expected_output(output_id.as_deref(), &record.result.created_seq_ids)
            )],
            Operation::DesignPrimerPairs {
                template,
                report_id,
                min_amplicon_bp,
                max_amplicon_bp,
                ..
            } => vec![format!(
                "Use primer-design report `{}` for template `{template}`; expected amplicon window: {min_amplicon_bp}..{max_amplicon_bp} bp. Order or pick primers only after review.",
                report_id.as_deref().unwrap_or("GENtle-generated")
            )],
            Operation::PrepareRestrictionCloningPcrHandoff {
                template,
                primer_report_id,
                pair_index,
                destination_vector_seq_id,
                forward_enzyme,
                reverse_enzyme,
                ..
            } => vec![format!(
                "Prepare restriction-cloning PCR from template `{template}` using primer report `{primer_report_id}` pair {pair_index}; destination vector `{destination_vector_seq_id}`; enzyme setup: {}{}.",
                forward_enzyme,
                reverse_enzyme
                    .as_deref()
                    .map(|enzyme| format!(" / {enzyme}"))
                    .unwrap_or_default()
            )],
            Operation::ApplyGibsonAssemblyPlan { plan_json } => {
                Self::lab_steps_for_gibson_plan(plan_json, &record.result.created_seq_ids)
            }
            Operation::Ligation {
                inputs,
                circularize_if_possible,
                output_id,
                protocol,
                ..
            } => vec![format!(
                "Ligate {} using {:?} ligation{}; expected product ID `{}`.",
                Self::lab_join_or_dash(inputs),
                protocol,
                if *circularize_if_possible {
                    " with circularization allowed"
                } else {
                    ""
                },
                Self::lab_expected_output(output_id.as_deref(), &record.result.created_seq_ids)
            )],
            Operation::LigationContainer {
                container_id,
                circularize_if_possible,
                output_id,
                protocol,
                ..
            } => vec![format!(
                "Ligate products in container `{container_id}` using {:?} ligation{}; expected product ID `{}`.",
                protocol,
                if *circularize_if_possible {
                    " with circularization allowed"
                } else {
                    ""
                },
                Self::lab_expected_output(output_id.as_deref(), &record.result.created_seq_ids)
            )],
            Operation::MergeContainers {
                inputs,
                output_prefix,
            } => vec![format!(
                "Pool sequence products {} into a shared container using output prefix `{}`.",
                Self::lab_join_or_dash(inputs),
                output_prefix.as_deref().unwrap_or("GENtle-generated")
            )],
            Operation::MergeContainersById {
                container_ids,
                output_prefix,
            } => vec![format!(
                "Pool containers {} using output prefix `{}`.",
                Self::lab_join_or_dash(container_ids),
                output_prefix.as_deref().unwrap_or("GENtle-generated")
            )],
            Operation::CreateArrangementSerial {
                container_ids,
                arrangement_id,
                ..
            } => vec![format!(
                "Arrange containers {} in the listed lane/order sequence; arrangement ID `{}`.",
                Self::lab_join_or_dash(container_ids),
                arrangement_id.as_deref().unwrap_or("GENtle-generated")
            )],
            Operation::RenderPoolGelSvg {
                inputs,
                container_ids,
                arrangement_id,
                ..
            } => vec![format!(
                "Use the rendered gel plan to verify expected band positions for inputs `{}`{}{}.",
                Self::lab_join_or_dash(inputs),
                container_ids
                    .as_ref()
                    .map(|ids| format!(", containers `{}`", Self::lab_join_or_dash(ids)))
                    .unwrap_or_default(),
                arrangement_id
                    .as_deref()
                    .map(|id| format!(", arrangement `{id}`"))
                    .unwrap_or_default()
            )],
            Operation::CreateRackFromArrangement {
                arrangement_id,
                rack_id,
                ..
            } => vec![format!(
                "Prepare the physical rack/plate layout from arrangement `{arrangement_id}`; rack ID `{}`.",
                rack_id.as_deref().unwrap_or("GENtle-generated")
            )],
            Operation::ExportRackLabelsSvg {
                rack_id,
                arrangement_id,
                ..
            } => vec![format!(
                "Print or copy rack labels for rack `{rack_id}`{} and apply them before sample transfer.",
                arrangement_id
                    .as_deref()
                    .map(|id| format!(" / arrangement `{id}`"))
                    .unwrap_or_default()
            )],
            Operation::SaveFile { seq_id, path, .. } => vec![format!(
                "Use exported sequence file `{path}` for designed sequence `{seq_id}` as the digital reference."
            )],
            _ => vec![],
        }
    }

    fn lab_steps_for_gibson_plan(plan_json: &str, created_seq_ids: &[String]) -> Vec<String> {
        let Ok(plan) =
            serde_json::from_str::<crate::gibson_planning::GibsonAssemblyPlan>(plan_json)
        else {
            return vec![
                "Run the recorded Gibson assembly plan; the plan JSON could not be summarized, so review the run bundle before bench work."
                    .to_string(),
            ];
        };
        let fragments = plan
            .fragments
            .iter()
            .map(|fragment| {
                format!(
                    "`{}` ({}; {} orientation)",
                    fragment.seq_id, fragment.role, fragment.orientation
                )
            })
            .collect::<Vec<_>>()
            .join(", ");
        let order = plan
            .assembly_order
            .iter()
            .map(|member| format!("{} `{}`", member.kind, member.id))
            .collect::<Vec<_>>()
            .join(" -> ");
        let overlaps = plan
            .junctions
            .iter()
            .filter_map(|junction| {
                junction
                    .required_overlap_bp
                    .map(|bp| format!("{}: {bp} bp", junction.id))
            })
            .collect::<Vec<_>>();
        let mut steps = vec![
            format!(
                "Prepare Gibson destination `{}` opened at `{}` and fragment(s): {}.",
                plan.destination.seq_id,
                plan.destination.opening.label,
                if fragments.is_empty() {
                    "none listed".to_string()
                } else {
                    fragments
                }
            ),
            format!(
                "Assemble in the designed order: {}.",
                if order.is_empty() {
                    "see GENtle plan".to_string()
                } else {
                    order
                }
            ),
        ];
        if !overlaps.is_empty() {
            steps.push(format!(
                "Check that ordered fragments provide the designed Gibson overlaps: {}.",
                overlaps.join(", ")
            ));
        }
        steps.push(format!(
            "Expected assembled product ID: `{}` (requested topology: {}).",
            Self::lab_expected_output(
                (!plan.product.output_id_hint.trim().is_empty())
                    .then_some(plan.product.output_id_hint.as_str()),
                created_seq_ids,
            ),
            if plan.product.topology.trim().is_empty() {
                "not specified"
            } else {
                plan.product.topology.trim()
            }
        ));
        steps
    }

    fn lab_expected_output(explicit: Option<&str>, created_seq_ids: &[String]) -> String {
        explicit
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .or_else(|| created_seq_ids.first().cloned())
            .unwrap_or_else(|| "not recorded".to_string())
    }

    fn lab_join_or_dash(values: &[String]) -> String {
        if values.is_empty() {
            "-".to_string()
        } else {
            values.join(", ")
        }
    }

    fn lab_oligo_display(sequence: &str) -> String {
        let trimmed = sequence.trim();
        if trimmed.len() <= 80 {
            trimmed.to_string()
        } else {
            format!("{}... ({} nt)", &trimmed[..80], trimmed.len())
        }
    }

    pub(crate) fn reverse_complement(seq: &str) -> String {
        seq.as_bytes()
            .iter()
            .rev()
            .map(|c| IupacCode::letter_complement(*c))
            .map(char::from)
            .collect()
    }

    pub(super) fn reverse_complement_bytes(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|c| IupacCode::letter_complement(*c))
            .collect()
    }

    pub(super) fn complement_bytes(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .map(|c| IupacCode::letter_complement(*c))
            .collect()
    }

    pub(super) fn normalize_dna_text(seq: &str) -> String {
        let cleaned = DNAsequence::validate_dna_sequence(seq.as_bytes());
        String::from_utf8_lossy(&cleaned).to_string()
    }

    pub(super) fn normalize_iupac_text(seq: &str) -> Result<String, EngineError> {
        let upper = seq.trim().to_ascii_uppercase();
        if upper.is_empty() {
            return Ok(upper);
        }
        for b in upper.as_bytes() {
            if !IupacCode::is_valid_letter(*b) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Invalid IUPAC nucleotide '{}'", *b as char),

                    cause_chain: vec![],
                });
            }
        }
        Ok(upper)
    }

    pub(super) fn resolve_tf_motif_or_iupac(token: &str) -> Result<String, EngineError> {
        let normalized = token.trim().to_ascii_uppercase();
        if normalized.is_empty() {
            return Ok(normalized);
        }
        if normalized
            .as_bytes()
            .iter()
            .all(|b| IupacCode::is_valid_letter(*b))
        {
            return Self::normalize_iupac_text(&normalized);
        }
        if let Some(motif) = tf_motifs::resolve_motif_definition(&normalized) {
            return Self::normalize_iupac_text(&motif.consensus_iupac);
        }
        Err(EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "TF motif '{}' was neither valid IUPAC text nor found in local motif registry",
                token
            ),

            cause_chain: vec![],
        })
    }

    pub(super) fn iupac_letter_counts(letter: u8) -> [f64; 4] {
        match letter.to_ascii_uppercase() {
            b'A' => [1.0, 0.0, 0.0, 0.0],
            b'C' => [0.0, 1.0, 0.0, 0.0],
            b'G' => [0.0, 0.0, 1.0, 0.0],
            b'T' | b'U' => [0.0, 0.0, 0.0, 1.0],
            b'M' => [1.0, 1.0, 0.0, 0.0], // A/C
            b'R' => [1.0, 0.0, 1.0, 0.0], // A/G
            b'W' => [1.0, 0.0, 0.0, 1.0], // A/T
            b'S' => [0.0, 1.0, 1.0, 0.0], // C/G
            b'Y' => [0.0, 1.0, 0.0, 1.0], // C/T
            b'K' => [0.0, 0.0, 1.0, 1.0], // G/T
            b'V' => [1.0, 1.0, 1.0, 0.0], // A/C/G
            b'H' => [1.0, 1.0, 0.0, 1.0], // A/C/T
            b'D' => [1.0, 0.0, 1.0, 1.0], // A/G/T
            b'B' => [0.0, 1.0, 1.0, 1.0], // C/G/T
            _ => [1.0, 1.0, 1.0, 1.0],    // N/unknown
        }
    }

    pub(super) fn matrix_from_iupac(consensus: &str) -> Vec<[f64; 4]> {
        consensus
            .as_bytes()
            .iter()
            .map(|b| Self::iupac_letter_counts(*b))
            .collect()
    }

    pub(super) fn resolve_tf_motif_for_scoring(
        token: &str,
    ) -> Result<(String, Option<String>, String, Vec<[f64; 4]>), EngineError> {
        let normalized = token.trim().to_ascii_uppercase();
        if normalized.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Empty TF motif token".to_string(),

                cause_chain: vec![],
            });
        }
        if let Some(motif) = tf_motifs::resolve_motif_definition(&normalized) {
            let consensus = Self::normalize_iupac_text(&motif.consensus_iupac)?;
            return Ok((motif.id, motif.name, consensus, motif.matrix_counts));
        }
        if normalized
            .as_bytes()
            .iter()
            .all(|b| IupacCode::is_valid_letter(*b))
        {
            let consensus = Self::normalize_iupac_text(&normalized)?;
            let matrix = Self::matrix_from_iupac(&consensus);
            return Ok((consensus.clone(), None, consensus, matrix));
        }
        Err(EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "TF motif '{}' was neither valid IUPAC text nor found in local motif registry",
                token
            ),

            cause_chain: vec![],
        })
    }

    pub(super) fn expand_tf_query_tokens(tokens: &[String]) -> Result<Vec<String>, EngineError> {
        let mut expanded = vec![];
        let mut seen = BTreeSet::new();
        for token in tokens {
            let trimmed = token.trim();
            if trimmed.is_empty() {
                continue;
            }
            if matches!(trimmed.to_ascii_uppercase().as_str(), "ALL" | "*") {
                for summary in crate::tf_motifs::list_motif_summaries() {
                    if seen.insert(summary.id.clone()) {
                        expanded.push(summary.name.unwrap_or(summary.id));
                    }
                }
                continue;
            }
            if tf_motifs::resolve_motif_definition(trimmed).is_some() {
                let (tf_id, tf_name, _consensus, _matrix_counts) =
                    Self::resolve_tf_motif_for_scoring(trimmed)?;
                if seen.insert(tf_id.clone()) {
                    expanded.push(tf_name.unwrap_or(tf_id));
                }
                continue;
            }
            if trimmed
                .as_bytes()
                .iter()
                .all(|b| IupacCode::is_valid_letter(*b))
            {
                if seen.insert(trimmed.to_ascii_uppercase()) {
                    expanded.push(trimmed.to_ascii_uppercase());
                }
                continue;
            }
            let resolution = tf_motifs::resolve_tf_query(trimmed);
            if resolution.matches.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "TF query '{}' did not match a local motif, a built-in TF group, or a family-like motif-name prefix",
                        trimmed
                    ),

                    cause_chain: vec![],
                });
            }
            for matched in resolution.matches {
                let candidate_token =
                    if tf_motifs::resolve_motif_definition(&matched.motif_id).is_some() {
                        matched.motif_id
                    } else if let Some(name) = matched
                        .motif_name
                        .as_deref()
                        .filter(|name| tf_motifs::resolve_motif_definition(name).is_some())
                    {
                        name.to_string()
                    } else {
                        matched.motif_id
                    };
                let (tf_id, tf_name, _consensus, _matrix_counts) =
                    Self::resolve_tf_motif_for_scoring(&candidate_token)?;
                if seen.insert(tf_id.clone()) {
                    expanded.push(tf_name.unwrap_or(tf_id));
                }
            }
        }
        if expanded.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "At least one TF motif or TF query is required".to_string(),

                cause_chain: vec![],
            });
        }
        Ok(expanded)
    }

    pub(super) fn resolve_tf_queries_report(
        queries: &[String],
    ) -> Result<TfQueryResolutionReport, EngineError> {
        let mut normalized_queries = queries
            .iter()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .collect::<Vec<_>>();
        if normalized_queries.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ResolveTfQueries requires at least one non-empty query".to_string(),

                cause_chain: vec![],
            });
        }
        normalized_queries.dedup();
        let entries = normalized_queries
            .iter()
            .map(|query| {
                let resolution = tf_motifs::resolve_tf_query(query);
                TfQueryResolutionEntry {
                    query: resolution.query,
                    normalized_query: resolution.normalized_query,
                    resolution_kind: resolution.resolution_kind,
                    label: resolution.label,
                    description: resolution.description,
                    aliases: resolution.aliases,
                    matches: resolution
                        .matches
                        .into_iter()
                        .map(|row| TfQueryResolutionMatch {
                            motif_id: row.motif_id,
                            motif_name: row.motif_name,
                            consensus_iupac: row.consensus_iupac,
                            motif_length_bp: row.motif_length_bp,
                        })
                        .collect(),
                    unresolved_reason: resolution.unresolved_reason,
                }
            })
            .collect::<Vec<_>>();
        let matched_motif_count = entries.iter().map(|entry| entry.matches.len()).sum();
        Ok(TfQueryResolutionReport {
            schema: TF_QUERY_RESOLUTION_REPORT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            queries: normalized_queries,
            returned_query_count: entries.len(),
            matched_motif_count,
            entries,
            warnings: vec![],
        })
    }

    pub(super) fn validate_tf_thresholds(min_llr_quantile: f64) -> Result<(), EngineError> {
        if !(0.0..=1.0).contains(&min_llr_quantile) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "TFBS min_llr_quantile must be between 0.0 and 1.0, got {}",
                    min_llr_quantile
                ),

                cause_chain: vec![],
            });
        }
        Ok(())
    }

    pub(super) fn format_tf_threshold_summary(min_llr_bits: f64, min_llr_quantile: f64) -> String {
        let mut parts: Vec<String> = vec![];
        if min_llr_bits.is_finite() {
            parts.push(format!("min_llr_bits={min_llr_bits}"));
        }
        if min_llr_quantile > 0.0 {
            parts.push(format!("min_llr_quantile={min_llr_quantile}"));
        }
        if parts.is_empty() {
            String::new()
        } else {
            format!(" with {}", parts.join(", "))
        }
    }

    pub(super) fn is_generated_tfbs_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated")
            .any(|v| v.eq_ignore_ascii_case("tfbs"))
    }

    pub(super) fn remove_generated_tfbs_features(features: &mut Vec<gb_io::seq::Feature>) {
        features.retain(|f| !Self::is_generated_tfbs_feature(f));
    }

    pub(super) fn build_tfbs_feature(
        start: usize,
        end: usize,
        reverse: bool,
        motif_len: usize,
        tf_id: &str,
        tf_name: Option<&str>,
        llr_bits: f64,
        llr_quantile: f64,
        true_log_odds_bits: f64,
        true_log_odds_quantile: f64,
    ) -> gb_io::seq::Feature {
        let base_location = gb_io::seq::Location::simple_range(start as i64, end as i64);
        let location = if reverse {
            gb_io::seq::Location::Complement(Box::new(base_location))
        } else {
            base_location
        };
        let mut qualifiers = vec![
            ("label".into(), Some(format!("TFBS {tf_id}"))),
            ("tf_id".into(), Some(tf_id.to_string())),
            ("motif_length_bp".into(), Some(motif_len.to_string())),
            ("llr_bits".into(), Some(format!("{llr_bits:.6}"))),
            ("llr_quantile".into(), Some(format!("{llr_quantile:.6}"))),
            (
                "true_log_odds_bits".into(),
                Some(format!("{true_log_odds_bits:.6}")),
            ),
            (
                "true_log_odds_quantile".into(),
                Some(format!("{true_log_odds_quantile:.6}")),
            ),
            (
                "log_odds_ratio_bits".into(),
                Some(format!("{true_log_odds_bits:.6}")),
            ),
            (
                "log_odds_ratio_quantile".into(),
                Some(format!("{true_log_odds_quantile:.6}")),
            ),
            (
                "quantile_scope".into(),
                Some("per_motif_windows_both_strands".to_string()),
            ),
            (
                "note".into(),
                Some(format!(
                    "tf_id={tf_id}; motif_length_bp={motif_len}; llr_bits={llr_bits:.4}; llr_quantile={llr_quantile:.4}; true_log_odds_bits={true_log_odds_bits:.4}; true_log_odds_quantile={true_log_odds_quantile:.4}"
                )),
            ),
            ("gentle_generated".into(), Some("tfbs".to_string())),
        ];
        if let Some(name) = tf_name
            && !name.trim().is_empty()
        {
            qualifiers.push(("bound_moiety".into(), Some(name.trim().to_string())));
        }
        gb_io::seq::Feature {
            kind: "TFBS".into(),
            location,
            qualifiers,
        }
    }

    pub(super) fn iupac_letter_complement(letter: u8) -> Option<u8> {
        match letter.to_ascii_uppercase() {
            b'A' => Some(b'T'),
            b'C' => Some(b'G'),
            b'G' => Some(b'C'),
            b'T' | b'U' => Some(b'A'),
            b'W' => Some(b'W'),
            b'S' => Some(b'S'),
            b'M' => Some(b'K'),
            b'K' => Some(b'M'),
            b'R' => Some(b'Y'),
            b'Y' => Some(b'R'),
            b'B' => Some(b'V'),
            b'D' => Some(b'H'),
            b'H' => Some(b'D'),
            b'V' => Some(b'B'),
            b'N' => Some(b'N'),
            _ => None,
        }
    }

    pub(super) fn reverse_complement_iupac(seq: &str) -> Result<String, EngineError> {
        let seq = Self::normalize_iupac_text(seq)?;
        let mut out = String::with_capacity(seq.len());
        for b in seq.as_bytes().iter().rev() {
            let c = Self::iupac_letter_complement(*b).ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid IUPAC nucleotide '{}'", *b as char),

                cause_chain: vec![],
            })?;
            out.push(c as char);
        }
        Ok(out)
    }

    pub(super) fn iupac_match_at(sequence: &[u8], pattern: &[u8], start: usize) -> bool {
        if pattern.is_empty() {
            return true;
        }
        if start > sequence.len() || start + pattern.len() > sequence.len() {
            return false;
        }
        for i in 0..pattern.len() {
            let s = IupacCode::from_letter(sequence[start + i]);
            let p = IupacCode::from_letter(pattern[i]);
            if s.is_empty() || p.is_empty() || s.subset(p).is_empty() {
                return false;
            }
        }
        true
    }

    pub(super) fn contains_iupac_pattern(sequence: &[u8], pattern: &[u8]) -> bool {
        if pattern.is_empty() {
            return true;
        }
        if sequence.len() < pattern.len() {
            return false;
        }
        for start in 0..=(sequence.len() - pattern.len()) {
            if Self::iupac_match_at(sequence, pattern, start) {
                return true;
            }
        }
        false
    }

    pub(super) fn contains_motif_any_strand(
        sequence: &[u8],
        motif: &str,
    ) -> Result<bool, EngineError> {
        let motif = Self::normalize_iupac_text(motif)?;
        if motif.is_empty() {
            return Ok(true);
        }
        let motif_bytes = motif.as_bytes();
        if Self::contains_iupac_pattern(sequence, motif_bytes) {
            return Ok(true);
        }
        let motif_rc = Self::reverse_complement_iupac(&motif)?;
        Ok(Self::contains_iupac_pattern(sequence, motif_rc.as_bytes()))
    }

    pub(super) fn normalized_sequence_for_quality(dna: &DNAsequence) -> Vec<u8> {
        dna.get_forward_string()
            .as_bytes()
            .iter()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|b| match b.to_ascii_uppercase() {
                b'U' => b'T',
                other => other,
            })
            .collect()
    }

    pub(crate) fn sequence_gc_fraction(sequence: &[u8]) -> Option<f64> {
        let mut canonical = 0usize;
        let mut gc = 0usize;
        for b in sequence {
            match b.to_ascii_uppercase() {
                b'G' | b'C' => {
                    canonical += 1;
                    gc += 1;
                }
                b'A' | b'T' => {
                    canonical += 1;
                }
                _ => {}
            }
        }
        if canonical == 0 {
            None
        } else {
            Some(gc as f64 / canonical as f64)
        }
    }

    pub(super) fn max_homopolymer_run(sequence: &[u8]) -> usize {
        let mut best = 0usize;
        let mut current = 0usize;
        let mut prev = 0u8;
        for b in sequence {
            let base = b.to_ascii_uppercase();
            if !matches!(base, b'A' | b'C' | b'G' | b'T') {
                current = 0;
                prev = 0;
                continue;
            }
            if base == prev {
                current += 1;
            } else {
                current = 1;
                prev = base;
            }
            if current > best {
                best = current;
            }
        }
        best
    }

    pub(super) fn contains_u6_terminator_t4(sequence: &[u8]) -> bool {
        sequence.windows(4).any(|w| w == b"TTTT")
    }

    pub(super) fn has_ambiguous_bases(sequence: &[u8]) -> bool {
        sequence
            .iter()
            .any(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    }

    pub(super) fn feature_labels(feature: &gb_io::seq::Feature) -> Vec<String> {
        let mut labels = Vec::new();
        for key in [
            "label",
            "gene",
            "locus_tag",
            "product",
            "standard_name",
            "note",
        ] {
            for value in feature.qualifier_values(key) {
                let v = value.trim();
                if !v.is_empty() {
                    labels.push(v.to_string());
                }
            }
        }
        labels
    }

    pub(super) fn resolve_sequence_anchor_position(
        dna: &DNAsequence,
        anchor: &SequenceAnchor,
        anchor_name: &str,
    ) -> Result<usize, EngineError> {
        match anchor {
            SequenceAnchor::Position { zero_based } => {
                if *zero_based > dna.len() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Sequence anchor '{anchor_name}' position {} is out of bounds for sequence length {}",
                            zero_based,
                            dna.len()
                        ),

                        cause_chain: vec![],
                    });
                }
                Ok(*zero_based)
            }
            SequenceAnchor::FeatureBoundary {
                feature_kind,
                feature_label,
                boundary,
                occurrence,
            } => {
                let kind_filter = feature_kind.as_ref().map(|s| s.to_ascii_uppercase());
                let label_filter = feature_label.as_ref().map(|s| s.to_ascii_uppercase());
                let mut matches = Vec::new();
                for feature in dna.features() {
                    if feature.kind.to_string().eq_ignore_ascii_case("SOURCE") {
                        continue;
                    }
                    if let Some(expected_kind) = &kind_filter
                        && feature.kind.to_string().to_ascii_uppercase() != *expected_kind
                    {
                        continue;
                    }
                    if let Some(expected_label) = &label_filter {
                        let labels = Self::feature_labels(feature);
                        let found_label = labels.iter().any(|label| {
                            let upper = label.to_ascii_uppercase();
                            upper == *expected_label || upper.contains(expected_label)
                        });
                        if !found_label {
                            continue;
                        }
                    }
                    let Ok((from, to)) = feature.location.find_bounds() else {
                        continue;
                    };
                    if from < 0 || to < 0 {
                        continue;
                    }
                    let start = from as usize;
                    let end = to as usize;
                    let pos = match boundary {
                        AnchorBoundary::Start => start,
                        AnchorBoundary::End => end,
                        AnchorBoundary::Middle => start + (end.saturating_sub(start) / 2),
                    };
                    if pos <= dna.len() {
                        matches.push(pos);
                    }
                }
                if matches.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("No feature matched sequence anchor '{}'", anchor_name),

                        cause_chain: vec![],
                    });
                }
                matches.sort_unstable();
                let idx = occurrence.unwrap_or(0);
                matches.get(idx).cloned().ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Sequence anchor '{}' occurrence {} was requested, but only {} match(es) found",
                        anchor_name, idx, matches.len()
                    ),
                
                    cause_chain: vec![],})
            }
        }
    }

    pub(super) fn anchored_range(
        anchor_pos: usize,
        len: usize,
        direction: &AnchorDirection,
        seq_len: usize,
        circular: bool,
    ) -> Option<(usize, usize)> {
        if len == 0 || seq_len == 0 {
            return None;
        }
        match direction {
            AnchorDirection::Upstream => {
                if !circular {
                    if anchor_pos > seq_len || anchor_pos < len {
                        return None;
                    }
                    Some((anchor_pos - len, anchor_pos))
                } else {
                    if anchor_pos > seq_len || len > seq_len {
                        return None;
                    }
                    let start = if anchor_pos >= len {
                        anchor_pos - len
                    } else {
                        seq_len - (len - anchor_pos)
                    };
                    Some((start, anchor_pos))
                }
            }
            AnchorDirection::Downstream => {
                if !circular {
                    if anchor_pos > seq_len || anchor_pos + len > seq_len {
                        return None;
                    }
                    Some((anchor_pos, anchor_pos + len))
                } else {
                    if anchor_pos > seq_len || len > seq_len {
                        return None;
                    }
                    Some((anchor_pos, anchor_pos + len))
                }
            }
        }
    }

    pub(super) fn primer_options(seq: &str) -> Result<Vec<Vec<u8>>, EngineError> {
        let mut ret = Vec::with_capacity(seq.len());
        for b in seq.as_bytes() {
            let opts = IupacCode::from_letter(*b).to_vec();
            if opts.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Invalid primer base '{}'", *b as char),

                    cause_chain: vec![],
                });
            }
            ret.push(opts);
        }
        Ok(ret)
    }

    pub(super) fn total_primer_variants(options: &[Vec<u8>]) -> usize {
        options
            .iter()
            .fold(1usize, |acc, v| acc.saturating_mul(v.len().max(1)))
    }

    pub(super) fn primer_variant_by_index(options: &[Vec<u8>], mut idx: usize) -> String {
        if options.is_empty() {
            return String::new();
        }
        let mut out: Vec<u8> = vec![b'A'; options.len()];
        for pos in (0..options.len()).rev() {
            let radix = options[pos].len();
            let choice = idx % radix;
            out[pos] = options[pos][choice];
            idx /= radix;
        }
        String::from_utf8(out).unwrap_or_default()
    }

    pub(super) fn expand_primer_variants(
        spec: &PcrPrimerSpec,
        cap: usize,
    ) -> Result<Vec<String>, EngineError> {
        let normalized = Self::normalize_iupac_text(&spec.sequence)?;
        if normalized.is_empty() {
            return Ok(vec![]);
        }
        let options = Self::primer_options(&normalized)?;
        let total = Self::total_primer_variants(&options);
        if total == 0 {
            return Ok(vec![]);
        }

        let max_variants = spec.max_variants.unwrap_or(total).min(cap).max(1);
        let mode = spec
            .library_mode
            .clone()
            .unwrap_or(PrimerLibraryMode::Enumerate);
        match mode {
            PrimerLibraryMode::Enumerate => {
                if total > max_variants {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Primer variant space ({total}) exceeds max_variants ({max_variants}); use Sample mode or raise limits"
                        ),

                        cause_chain: vec![],
                    });
                }
                let mut ret = Vec::with_capacity(total);
                for idx in 0..total {
                    ret.push(Self::primer_variant_by_index(&options, idx));
                }
                Ok(ret)
            }
            PrimerLibraryMode::Sample => {
                let target = max_variants.min(total);
                let mut chosen: Vec<usize> = Vec::with_capacity(target);
                let mut seen: HashSet<usize> = HashSet::with_capacity(target * 2);
                let mut state = spec.sample_seed.unwrap_or(0x9E3779B97F4A7C15);

                while chosen.len() < target {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    let idx = (state as usize) % total;
                    if seen.insert(idx) {
                        chosen.push(idx);
                    }
                    if seen.len() == total {
                        break;
                    }
                }
                chosen.sort_unstable();
                let mut ret = Vec::with_capacity(chosen.len());
                for idx in chosen {
                    ret.push(Self::primer_variant_by_index(&options, idx));
                }
                Ok(ret)
            }
        }
    }

    pub(super) fn find_subsequence(haystack: &[u8], needle: &[u8], start: usize) -> Option<usize> {
        if needle.is_empty() || haystack.len() < needle.len() || start >= haystack.len() {
            return None;
        }
        let end = haystack.len() - needle.len();
        (start..=end).find(|idx| &haystack[*idx..*idx + needle.len()] == needle)
    }

    pub(crate) fn find_all_subsequences(haystack: &[u8], needle: &[u8]) -> Vec<usize> {
        let mut ret = vec![];
        if needle.is_empty() || haystack.len() < needle.len() {
            return ret;
        }
        let mut start = 0usize;
        while let Some(pos) = Self::find_subsequence(haystack, needle, start) {
            ret.push(pos);
            start = pos + 1;
            if start >= haystack.len() {
                break;
            }
        }
        ret
    }

    pub(crate) fn primer_tm_model_description() -> String {
        format!(
            "Displayed Tm values use GENtle's shared Thermo Fisher-style modified Allawi/SantaLucia nearest-neighbor estimate with fixed assumptions: exact complement, {:.0} mM effective monovalent salt, {:.0} nM primer concentration, empirical high-fidelity adjustment, and no mismatch/dangling-end/Mg correction. Very short or ambiguous sequences fall back to the simple 2/4 estimate.",
            Self::primer_tm_effective_monovalent_salt_molar() * 1_000.0,
            Self::primer_tm_primer_concentration_molar() * 1_000_000_000.0
        )
    }

    pub(crate) fn estimate_primer_tm_c(primer: &[u8]) -> f64 {
        if primer.is_empty() {
            return 0.0;
        }
        let Some(canonical) = Self::canonical_dna_bases(primer) else {
            return Self::estimate_primer_tm_wallace_c(primer);
        };
        if canonical.len() < 8 {
            return Self::estimate_primer_tm_wallace_c(&canonical);
        }
        Self::estimate_primer_tm_nearest_neighbor_c(&canonical)
            .unwrap_or_else(|| Self::estimate_primer_tm_wallace_c(&canonical))
    }

    fn canonical_dna_bases(primer: &[u8]) -> Option<Vec<u8>> {
        let mut out = Vec::with_capacity(primer.len());
        for base in primer {
            let upper = base.to_ascii_uppercase();
            match upper {
                b'A' | b'C' | b'G' | b'T' => out.push(upper),
                _ => return None,
            }
        }
        Some(out)
    }

    fn estimate_primer_tm_wallace_c(primer: &[u8]) -> f64 {
        let mut at = 0usize;
        let mut gc = 0usize;
        for base in primer {
            match base.to_ascii_uppercase() {
                b'A' | b'T' => at += 1,
                b'C' | b'G' => gc += 1,
                _ => {}
            }
        }
        (2 * at + 4 * gc) as f64
    }

    fn estimate_primer_tm_nearest_neighbor_c(primer: &[u8]) -> Option<f64> {
        const THERMOFISHER_TM_MIN_C: f64 = 0.0;
        const THERMOFISHER_TM_MAX_C: f64 = 95.0;

        let mut delta_h_kcal_per_mol = 0.0;
        let mut delta_s_cal_per_mol_k = 0.0;

        for terminal_base in [primer.first().copied()?, primer.last().copied()?] {
            if matches!(terminal_base, b'A' | b'T') {
                delta_h_kcal_per_mol += 2.3;
                delta_s_cal_per_mol_k += 4.1;
            } else {
                delta_h_kcal_per_mol += 0.1;
                delta_s_cal_per_mol_k -= 2.8;
            }
        }

        for pair in primer.windows(2) {
            let (pair_h, pair_s) = Self::primer_tm_nearest_neighbor_parameters(pair[0], pair[1])?;
            delta_h_kcal_per_mol += pair_h;
            delta_s_cal_per_mol_k += pair_s;
        }

        let concentration = Self::primer_tm_primer_concentration_molar();
        let concentration_term =
            (concentration / 4.0).ln() * Self::primer_tm_gas_constant_cal_per_mol_k();
        let denominator = delta_s_cal_per_mol_k + concentration_term;
        if !denominator.is_finite() || denominator >= 0.0 {
            return None;
        }

        let salt_adjustment_c = 16.6 * Self::primer_tm_effective_monovalent_salt_molar().log10();
        let raw_tm_c = (delta_h_kcal_per_mol * 1_000.0) / denominator + salt_adjustment_c - 273.15;
        if !raw_tm_c.is_finite() {
            return None;
        }
        let adjusted_tm_c = (raw_tm_c + 3.0) * Self::primer_tm_thermofisher_adjustment_slope()
            + Self::primer_tm_thermofisher_adjustment_intercept();
        if !adjusted_tm_c.is_finite() {
            return None;
        }
        Some(adjusted_tm_c.clamp(THERMOFISHER_TM_MIN_C, THERMOFISHER_TM_MAX_C))
    }

    fn primer_tm_nearest_neighbor_parameters(left: u8, right: u8) -> Option<(f64, f64)> {
        match (left, right) {
            (b'A', b'A') | (b'T', b'T') => Some((-7.9, -22.2)),
            (b'A', b'T') => Some((-7.2, -20.4)),
            (b'T', b'A') => Some((-7.2, -21.3)),
            (b'C', b'A') | (b'T', b'G') => Some((-8.5, -22.7)),
            (b'G', b'T') | (b'A', b'C') => Some((-8.4, -22.4)),
            (b'C', b'T') | (b'A', b'G') => Some((-7.8, -21.0)),
            (b'G', b'A') | (b'T', b'C') => Some((-8.2, -22.2)),
            (b'C', b'G') => Some((-10.6, -27.2)),
            (b'G', b'C') => Some((-9.8, -24.4)),
            (b'G', b'G') | (b'C', b'C') => Some((-8.0, -19.9)),
            _ => None,
        }
    }

    fn primer_tm_effective_monovalent_salt_molar() -> f64 {
        // Fixed effective salt term used by Thermo Fisher's v4 high-fidelity All97 path.
        0.215273974689348
    }

    fn primer_tm_primer_concentration_molar() -> f64 {
        500e-9
    }

    fn primer_tm_thermofisher_adjustment_slope() -> f64 {
        // Empirical high-fidelity adjustment from Thermo Fisher's v4 All97 path.
        0.9376798568
    }

    fn primer_tm_thermofisher_adjustment_intercept() -> f64 {
        4.5185404499
    }

    fn primer_tm_gas_constant_cal_per_mol_k() -> f64 {
        1.9872
    }

    pub(super) fn max_contiguous_match_run_with_shift(left: &[u8], right: &[u8]) -> usize {
        if left.is_empty() || right.is_empty() {
            return 0;
        }
        let mut best = 0usize;
        let min_shift = -(right.len() as isize) + 1;
        let max_shift = left.len() as isize - 1;
        for shift in min_shift..=max_shift {
            let mut run = 0usize;
            for (i, left_base) in left.iter().enumerate() {
                let right_idx = i as isize - shift;
                if right_idx < 0 || right_idx >= right.len() as isize {
                    run = 0;
                    continue;
                }
                let right_base = right[right_idx as usize];
                if left_base.eq_ignore_ascii_case(&right_base) {
                    run += 1;
                    best = best.max(run);
                } else {
                    run = 0;
                }
            }
        }
        best
    }

    pub(super) fn longest_suffix_match_in_target(source: &[u8], target: &[u8]) -> usize {
        if source.is_empty() || target.is_empty() {
            return 0;
        }
        for len in (1..=source.len()).rev() {
            let suffix = &source[source.len() - len..];
            if Self::find_subsequence(target, suffix, 0).is_some() {
                return len;
            }
        }
        0
    }

    pub(super) fn compute_primer_heuristic_metrics(sequence: &[u8]) -> PrimerHeuristicMetrics {
        let three_prime_base = sequence
            .last()
            .copied()
            .unwrap_or(b'N')
            .to_ascii_uppercase();
        let three_prime_gc_clamp = matches!(three_prime_base, b'G' | b'C');
        let self_complementary_run_bp = if sequence.is_empty() {
            0
        } else {
            let rc = Self::reverse_complement_bytes(sequence);
            Self::max_contiguous_match_run_with_shift(sequence, &rc)
        };
        PrimerHeuristicMetrics {
            length_bp: sequence.len(),
            three_prime_gc_clamp,
            three_prime_base,
            longest_homopolymer_run_bp: Self::max_homopolymer_run(sequence),
            self_complementary_run_bp,
        }
    }

    pub(super) fn compute_primer_pair_dimer_metrics(
        forward_sequence: &[u8],
        reverse_sequence: &[u8],
    ) -> PrimerPairDimerMetrics {
        if forward_sequence.is_empty() || reverse_sequence.is_empty() {
            return PrimerPairDimerMetrics::default();
        }
        let reverse_rc = Self::reverse_complement_bytes(reverse_sequence);
        let forward_rc = Self::reverse_complement_bytes(forward_sequence);
        let max_complementary_run_bp =
            Self::max_contiguous_match_run_with_shift(forward_sequence, &reverse_rc);
        let max_3prime_complementary_run_bp =
            Self::longest_suffix_match_in_target(forward_sequence, &reverse_rc).max(
                Self::longest_suffix_match_in_target(reverse_sequence, &forward_rc),
            );
        PrimerPairDimerMetrics {
            max_complementary_run_bp,
            max_3prime_complementary_run_bp,
        }
    }

    pub(super) fn compute_primer_self_3prime_complementary_run(sequence: &[u8]) -> usize {
        if sequence.is_empty() {
            return 0;
        }
        let rc = Self::reverse_complement_bytes(sequence);
        Self::longest_suffix_match_in_target(sequence, &rc)
    }

    fn oligo_qc_is_primer_role(role: &str) -> bool {
        matches!(role, "forward_primer" | "reverse_primer" | "primer")
    }

    fn oligo_qc_status_rank(status: &str) -> usize {
        match status {
            "fail" => 3,
            "warning" => 2,
            "pass" => 1,
            _ => 0,
        }
    }

    fn oligo_qc_join_warnings(warnings: &[String], pass_summary: String) -> String {
        if warnings.is_empty() {
            pass_summary
        } else {
            warnings.join(" ")
        }
    }

    pub(super) fn build_oligo_qc_report(
        assay_kind: &str,
        oligos: &[(&str, &str, &str)],
    ) -> OligoQcReport {
        // Conceptually follows Primer3's public distinction between broad
        // oligo self/pair complementarity and 3'-end anchored complementarity
        // (`SELF_ANY`, `SELF_END`, `COMPL_ANY`, `COMPL_END`, plus the `_TH`
        // thermodynamic variants). This is an independent GENtle exact-run
        // implementation; no Primer3 source code is vendored or translated.
        let mut oligo_records = Vec::with_capacity(oligos.len());
        let mut report_warnings = vec![];
        for (label, role, sequence) in oligos {
            let sequence_bytes = sequence.as_bytes();
            let metrics = Self::compute_primer_heuristic_metrics(sequence_bytes);
            let self_3prime_complementary_run_bp =
                Self::compute_primer_self_3prime_complementary_run(sequence_bytes);
            let mut warnings = vec![];
            if metrics.longest_homopolymer_run_bp > PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP {
                warnings.push(format!(
                    "{label} has a {} bp homopolymer run (recommended <= {} bp).",
                    metrics.longest_homopolymer_run_bp, PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP
                ));
            }
            if metrics.self_complementary_run_bp > PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP
            {
                warnings.push(format!(
                    "{label} has a {} bp self-complementary reverse-complement run (recommended <= {} bp).",
                    metrics.self_complementary_run_bp,
                    PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP
                ));
            }
            if self_3prime_complementary_run_bp > PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP {
                warnings.push(format!(
                    "{label} has a {} bp 3'-anchored self-complementary run (recommended <= {} bp).",
                    self_3prime_complementary_run_bp,
                    PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP
                ));
            }
            let status = if warnings.is_empty() {
                "pass"
            } else {
                "warning"
            }
            .to_string();
            let pass_summary = format!(
                "{label} is {} bp; self-complementarity max={} bp, 3'-self max={} bp.",
                sequence.len(),
                metrics.self_complementary_run_bp,
                self_3prime_complementary_run_bp
            );
            let summary = Self::oligo_qc_join_warnings(&warnings, pass_summary);
            report_warnings.extend(warnings.iter().cloned());
            oligo_records.push(OligoQcOligoRecord {
                label: (*label).to_string(),
                role: (*role).to_string(),
                sequence: (*sequence).to_string(),
                length_bp: sequence.len(),
                gc_fraction: Self::sequence_gc_fraction(sequence_bytes).unwrap_or(0.0),
                three_prime_base: char::from(metrics.three_prime_base).to_string(),
                three_prime_gc_clamp: metrics.three_prime_gc_clamp,
                longest_homopolymer_run_bp: metrics.longest_homopolymer_run_bp,
                self_complementary_run_bp: metrics.self_complementary_run_bp,
                self_3prime_complementary_run_bp,
                status,
                summary,
                warnings,
            });
        }

        let mut interaction_records = vec![];
        for left_idx in 0..oligos.len() {
            for right_idx in left_idx + 1..oligos.len() {
                let (left_label, left_role, left_sequence) = oligos[left_idx];
                let (right_label, right_role, right_sequence) = oligos[right_idx];
                let left_sequence_bytes = left_sequence.as_bytes();
                let right_sequence_bytes = right_sequence.as_bytes();
                let dimer_metrics = Self::compute_primer_pair_dimer_metrics(
                    left_sequence_bytes,
                    right_sequence_bytes,
                );
                let right_rc = Self::reverse_complement_bytes(right_sequence_bytes);
                let left_rc = Self::reverse_complement_bytes(left_sequence_bytes);
                let left_3prime_complementary_run_bp =
                    Self::longest_suffix_match_in_target(left_sequence_bytes, &right_rc);
                let right_3prime_complementary_run_bp =
                    Self::longest_suffix_match_in_target(right_sequence_bytes, &left_rc);
                let left_is_primer = Self::oligo_qc_is_primer_role(left_role);
                let right_is_primer = Self::oligo_qc_is_primer_role(right_role);
                let primer_3prime_extension_risk = (left_is_primer
                    && left_3prime_complementary_run_bp
                        > PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP)
                    || (right_is_primer
                        && right_3prime_complementary_run_bp
                            > PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP);
                let mut warnings = vec![];
                if dimer_metrics.max_complementary_run_bp > PRIMER_RECOMMENDED_MAX_PAIR_DIMER_RUN_BP
                {
                    warnings.push(format!(
                        "{left_label}/{right_label} have a {} bp reverse-complement interaction run (recommended <= {} bp).",
                        dimer_metrics.max_complementary_run_bp,
                        PRIMER_RECOMMENDED_MAX_PAIR_DIMER_RUN_BP
                    ));
                }
                if left_is_primer
                    && left_3prime_complementary_run_bp
                        > PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP
                {
                    warnings.push(format!(
                        "{left_label} 3' end has {} bp complementarity to {right_label} (recommended <= {} bp).",
                        left_3prime_complementary_run_bp,
                        PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP
                    ));
                }
                if right_is_primer
                    && right_3prime_complementary_run_bp
                        > PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP
                {
                    warnings.push(format!(
                        "{right_label} 3' end has {} bp complementarity to {left_label} (recommended <= {} bp).",
                        right_3prime_complementary_run_bp,
                        PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP
                    ));
                }
                let status = if left_is_primer && right_is_primer && primer_3prime_extension_risk {
                    "fail"
                } else if warnings.is_empty() {
                    "pass"
                } else {
                    "warning"
                }
                .to_string();
                let pass_summary = format!(
                    "{left_label}/{right_label} max interaction={} bp, max 3' primer-side interaction={} bp.",
                    dimer_metrics.max_complementary_run_bp,
                    dimer_metrics.max_3prime_complementary_run_bp
                );
                let summary = Self::oligo_qc_join_warnings(&warnings, pass_summary);
                report_warnings.extend(warnings.iter().cloned());
                interaction_records.push(OligoQcInteractionRecord {
                    left_label: left_label.to_string(),
                    right_label: right_label.to_string(),
                    left_role: left_role.to_string(),
                    right_role: right_role.to_string(),
                    max_complementary_run_bp: dimer_metrics.max_complementary_run_bp,
                    left_3prime_complementary_run_bp,
                    right_3prime_complementary_run_bp,
                    max_3prime_complementary_run_bp: dimer_metrics.max_3prime_complementary_run_bp,
                    primer_3prime_extension_risk,
                    status,
                    summary,
                    warnings,
                });
            }
        }

        let status = oligo_records
            .iter()
            .map(|record| record.status.as_str())
            .chain(
                interaction_records
                    .iter()
                    .map(|record| record.status.as_str()),
            )
            .max_by_key(|status| Self::oligo_qc_status_rank(status))
            .unwrap_or("pass")
            .to_string();
        let summary = if report_warnings.is_empty() {
            format!(
                "Oligo QC passed exact reverse-complement run screening for {} oligo(s) and {} interaction(s).",
                oligo_records.len(),
                interaction_records.len()
            )
        } else {
            format!(
                "Oligo QC found {} warning/failure signal(s) across {} oligo(s) and {} interaction(s).",
                report_warnings.len(),
                oligo_records.len(),
                interaction_records.len()
            )
        };
        OligoQcReport {
            schema: OLIGO_QC_REPORT_SCHEMA.to_string(),
            assay_kind: assay_kind.to_string(),
            method: "exact_reverse_complement_run_screen".to_string(),
            method_reference:
                "Primer3-style oligo QC vocabulary: SELF_ANY/SELF_END and COMPL_ANY/COMPL_END; independent GENtle exact-run implementation"
                    .to_string(),
            status,
            summary,
            oligo_count: oligo_records.len(),
            interaction_count: interaction_records.len(),
            oligos: oligo_records,
            interactions: interaction_records,
            warnings: report_warnings,
        }
    }

    pub(super) fn preferred_primer_length_penalty(length_bp: usize) -> f64 {
        if length_bp < PRIMER_PREFERRED_MIN_LENGTH_BP {
            let delta = (PRIMER_PREFERRED_MIN_LENGTH_BP - length_bp) as f64;
            return delta * 3.0;
        }
        if length_bp > PRIMER_PREFERRED_MAX_LENGTH_BP {
            let delta = (length_bp - PRIMER_PREFERRED_MAX_LENGTH_BP) as f64;
            return delta * 1.5;
        }
        0.0
    }

    pub(super) fn primer_secondary_structure_penalty(metrics: PrimerHeuristicMetrics) -> f64 {
        let homopolymer_penalty = if metrics.longest_homopolymer_run_bp
            > PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP
        {
            (metrics.longest_homopolymer_run_bp - PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP) as f64
                * 10.0
        } else {
            0.0
        };
        let self_complementary_penalty = if metrics.self_complementary_run_bp
            > PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP
        {
            (metrics.self_complementary_run_bp - PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP)
                as f64
                * 12.0
        } else {
            0.0
        };
        homopolymer_penalty + self_complementary_penalty
    }

    pub(super) fn primer_pair_dimer_penalty(metrics: PrimerPairDimerMetrics) -> f64 {
        let broad = if metrics.max_complementary_run_bp > PRIMER_RECOMMENDED_MAX_PAIR_DIMER_RUN_BP {
            (metrics.max_complementary_run_bp - PRIMER_RECOMMENDED_MAX_PAIR_DIMER_RUN_BP) as f64
                * 20.0
        } else {
            0.0
        };
        let three_prime = if metrics.max_3prime_complementary_run_bp
            > PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP
        {
            (metrics.max_3prime_complementary_run_bp
                - PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP) as f64
                * 35.0
        } else {
            0.0
        };
        broad + three_prime
    }

    pub(super) fn validate_primer_design_side_constraints(
        label: &str,
        side: &PrimerDesignSideConstraint,
    ) -> Result<(), EngineError> {
        if side.min_length == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("{label}.min_length must be >= 1"),

                cause_chain: vec![],
            });
        }
        if side.min_length > side.max_length {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "{label}.min_length ({}) must be <= {label}.max_length ({})",
                    side.min_length, side.max_length
                ),

                cause_chain: vec![],
            });
        }
        if !(0.0..=1.0).contains(&side.min_gc_fraction) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "{label}.min_gc_fraction ({}) must be between 0.0 and 1.0",
                    side.min_gc_fraction
                ),

                cause_chain: vec![],
            });
        }
        if !(0.0..=1.0).contains(&side.max_gc_fraction) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "{label}.max_gc_fraction ({}) must be between 0.0 and 1.0",
                    side.max_gc_fraction
                ),

                cause_chain: vec![],
            });
        }
        if side.min_gc_fraction > side.max_gc_fraction {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "{label}.min_gc_fraction ({}) must be <= {label}.max_gc_fraction ({})",
                    side.min_gc_fraction, side.max_gc_fraction
                ),

                cause_chain: vec![],
            });
        }
        if side.min_tm_c > side.max_tm_c {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "{label}.min_tm_c ({}) must be <= {label}.max_tm_c ({})",
                    side.min_tm_c, side.max_tm_c
                ),

                cause_chain: vec![],
            });
        }
        if side.max_anneal_hits == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("{label}.max_anneal_hits must be >= 1"),

                cause_chain: vec![],
            });
        }
        let tail_5prime = side
            .non_annealing_5prime_tail
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(Self::normalize_iupac_text)
            .transpose()?;
        let tail_5prime_len = tail_5prime.as_ref().map(|v| v.len()).unwrap_or(0);
        let max_full_length = side.max_length.saturating_add(tail_5prime_len);
        if let Some(raw) = &side.fixed_5prime {
            let norm = Self::normalize_iupac_text(raw)?;
            if norm.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("{label}.fixed_5prime must not be empty"),

                    cause_chain: vec![],
                });
            }
            if norm.len() > max_full_length {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "{label}.fixed_5prime length ({}) must be <= max full primer length ({})",
                        norm.len(),
                        max_full_length
                    ),

                    cause_chain: vec![],
                });
            }
        }
        if let Some(raw) = &side.fixed_3prime {
            let norm = Self::normalize_iupac_text(raw)?;
            if norm.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("{label}.fixed_3prime must not be empty"),

                    cause_chain: vec![],
                });
            }
            if norm.len() > max_full_length {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "{label}.fixed_3prime length ({}) must be <= max full primer length ({})",
                        norm.len(),
                        max_full_length
                    ),

                    cause_chain: vec![],
                });
            }
        }
        for motif in &side.required_motifs {
            let norm = Self::normalize_iupac_text(motif)?;
            if norm.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("{label}.required_motifs contains an empty motif"),

                    cause_chain: vec![],
                });
            }
        }
        for motif in &side.forbidden_motifs {
            let norm = Self::normalize_iupac_text(motif)?;
            if norm.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("{label}.forbidden_motifs contains an empty motif"),

                    cause_chain: vec![],
                });
            }
        }
        let mut seen_locks: HashSet<usize> = HashSet::new();
        for lock in &side.locked_positions {
            if lock.offset_0based >= max_full_length {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "{label}.locked_positions offset {} must be < max full primer length ({})",
                        lock.offset_0based, max_full_length
                    ),

                    cause_chain: vec![],
                });
            }
            let norm = Self::normalize_iupac_text(&lock.base)?;
            if norm.len() != 1 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "{label}.locked_positions base '{}' must be a single IUPAC letter",
                        lock.base
                    ),

                    cause_chain: vec![],
                });
            }
            if !seen_locks.insert(lock.offset_0based) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "{label}.locked_positions has duplicate offset {}",
                        lock.offset_0based
                    ),

                    cause_chain: vec![],
                });
            }
        }
        if let (Some(start), Some(end)) = (side.start_0based, side.end_0based)
            && start >= end
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "{label}.start_0based ({start}) must be < {label}.end_0based ({end})"
                ),

                cause_chain: vec![],
            });
        }
        Ok(())
    }

    pub(super) fn normalize_primer_side_sequence_constraints(
        side: &PrimerDesignSideConstraint,
    ) -> Result<NormalizedPrimerSideSequenceConstraints, EngineError> {
        let non_annealing_5prime_tail = side
            .non_annealing_5prime_tail
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(Self::normalize_iupac_text)
            .transpose()?;
        let fixed_5prime = side
            .fixed_5prime
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(Self::normalize_iupac_text)
            .transpose()?;
        let fixed_3prime = side
            .fixed_3prime
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(Self::normalize_iupac_text)
            .transpose()?;
        let mut required_motifs = Vec::with_capacity(side.required_motifs.len());
        for motif in &side.required_motifs {
            let norm = Self::normalize_iupac_text(motif)?;
            if !norm.is_empty() {
                required_motifs.push(norm);
            }
        }
        let mut forbidden_motifs = Vec::with_capacity(side.forbidden_motifs.len());
        for motif in &side.forbidden_motifs {
            let norm = Self::normalize_iupac_text(motif)?;
            if !norm.is_empty() {
                forbidden_motifs.push(norm);
            }
        }
        let mut locked_positions = side
            .locked_positions
            .iter()
            .map(|lock| {
                let base = Self::normalize_iupac_text(&lock.base)?;
                Ok((lock.offset_0based, base))
            })
            .collect::<Result<Vec<_>, EngineError>>()?;
        locked_positions.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        Ok(NormalizedPrimerSideSequenceConstraints {
            non_annealing_5prime_tail,
            fixed_5prime,
            fixed_3prime,
            required_motifs,
            forbidden_motifs,
            locked_positions,
        })
    }

    pub(super) fn normalize_primer_pair_constraints(
        constraints: &PrimerDesignPairConstraint,
    ) -> Result<NormalizedPrimerPairConstraints, EngineError> {
        if let (Some(start), Some(end)) = (
            constraints.fixed_amplicon_start_0based,
            constraints.fixed_amplicon_end_0based_exclusive,
        ) && start >= end
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "pair_constraints.fixed_amplicon_start_0based ({start}) must be < pair_constraints.fixed_amplicon_end_0based_exclusive ({end})"
                ),

                cause_chain: vec![],
            });
        }
        let mut required_amplicon_motifs =
            Vec::with_capacity(constraints.required_amplicon_motifs.len());
        for motif in &constraints.required_amplicon_motifs {
            let norm = Self::normalize_iupac_text(motif)?;
            if !norm.is_empty() {
                required_amplicon_motifs.push(norm);
            }
        }
        let mut forbidden_amplicon_motifs =
            Vec::with_capacity(constraints.forbidden_amplicon_motifs.len());
        for motif in &constraints.forbidden_amplicon_motifs {
            let norm = Self::normalize_iupac_text(motif)?;
            if !norm.is_empty() {
                forbidden_amplicon_motifs.push(norm);
            }
        }
        Ok(NormalizedPrimerPairConstraints {
            require_roi_flanking: constraints.require_roi_flanking,
            required_amplicon_motifs,
            forbidden_amplicon_motifs,
            fixed_amplicon_start_0based: constraints.fixed_amplicon_start_0based,
            fixed_amplicon_end_0based_exclusive: constraints.fixed_amplicon_end_0based_exclusive,
        })
    }

    pub(super) fn generate_primer_side_candidates(
        template: &[u8],
        side: &PrimerDesignSideConstraint,
        sequence_constraints: &NormalizedPrimerSideSequenceConstraints,
        reverse_orientation: bool,
        rejections: &mut PrimerDesignRejectionSummary,
    ) -> Vec<PrimerDesignCandidate> {
        let mut ret = vec![];
        if template.is_empty() || side.min_length > template.len() {
            return ret;
        }
        let max_start = template.len().saturating_sub(side.min_length);
        for start in 0..=max_start {
            if let Some(location) = side.location_0based
                && start != location
            {
                continue;
            }
            if let Some(min_start) = side.start_0based
                && start < min_start
            {
                continue;
            }
            for length in side.min_length..=side.max_length {
                let Some(end) = start.checked_add(length) else {
                    continue;
                };
                if end > template.len() {
                    break;
                }
                if let Some(max_end) = side.end_0based
                    && end > max_end
                {
                    rejections.out_of_window = rejections.out_of_window.saturating_add(1);
                    continue;
                }
                let binding_window = &template[start..end];
                let anneal_bytes = if reverse_orientation {
                    Self::reverse_complement_bytes(binding_window)
                } else {
                    binding_window.to_vec()
                };
                let mut primer_bytes: Vec<u8> = Vec::new();
                if let Some(tail) = &sequence_constraints.non_annealing_5prime_tail {
                    primer_bytes.extend_from_slice(tail.as_bytes());
                }
                primer_bytes.extend_from_slice(&anneal_bytes);
                let gc_fraction = Self::sequence_gc_fraction(&anneal_bytes).unwrap_or(0.0);
                let tm_c = Self::estimate_primer_tm_c(&anneal_bytes);
                if gc_fraction < side.min_gc_fraction
                    || gc_fraction > side.max_gc_fraction
                    || tm_c < side.min_tm_c
                    || tm_c > side.max_tm_c
                {
                    rejections.gc_or_tm_out_of_bounds =
                        rejections.gc_or_tm_out_of_bounds.saturating_add(1);
                    continue;
                }
                let anneal_hits = Self::find_all_subsequences(template, binding_window).len();
                if anneal_hits == 0 || anneal_hits > side.max_anneal_hits {
                    rejections.non_unique_anneal = rejections.non_unique_anneal.saturating_add(1);
                    continue;
                }
                if !Self::primer_sequence_matches_side_constraints(
                    &primer_bytes,
                    sequence_constraints,
                ) {
                    rejections.primer_constraint_failure =
                        rejections.primer_constraint_failure.saturating_add(1);
                    continue;
                }
                ret.push(PrimerDesignCandidate {
                    sequence: String::from_utf8(primer_bytes).unwrap_or_default(),
                    start_0based: start,
                    end_0based_exclusive: end,
                    tm_c,
                    gc_fraction,
                    anneal_hits,
                });
            }
        }
        ret.sort_by(|a, b| {
            a.start_0based
                .cmp(&b.start_0based)
                .then(a.end_0based_exclusive.cmp(&b.end_0based_exclusive))
                .then(a.sequence.cmp(&b.sequence))
        });
        ret
    }

    pub(super) fn primer_sequence_matches_side_constraints(
        primer_sequence: &[u8],
        constraints: &NormalizedPrimerSideSequenceConstraints,
    ) -> bool {
        if let Some(prefix) = &constraints.fixed_5prime
            && (primer_sequence.len() < prefix.len()
                || !Self::iupac_match_at(primer_sequence, prefix.as_bytes(), 0))
        {
            return false;
        }
        if let Some(suffix) = &constraints.fixed_3prime
            && (primer_sequence.len() < suffix.len()
                || !Self::iupac_match_at(
                    primer_sequence,
                    suffix.as_bytes(),
                    primer_sequence.len().saturating_sub(suffix.len()),
                ))
        {
            return false;
        }
        for motif in &constraints.required_motifs {
            if !Self::contains_iupac_pattern(primer_sequence, motif.as_bytes()) {
                return false;
            }
        }
        for motif in &constraints.forbidden_motifs {
            if Self::contains_iupac_pattern(primer_sequence, motif.as_bytes()) {
                return false;
            }
        }
        for (offset, base) in &constraints.locked_positions {
            if *offset >= primer_sequence.len()
                || !Self::iupac_match_at(primer_sequence, base.as_bytes(), *offset)
            {
                return false;
            }
        }
        true
    }

    pub(super) fn render_primer_design_report_id(raw: Option<String>, template: &str) -> String {
        let candidate = raw
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .unwrap_or_else(|| {
                let stem = template
                    .chars()
                    .map(|ch| {
                        if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                            ch
                        } else {
                            '_'
                        }
                    })
                    .collect::<String>();
                format!("primer_report_{}_{}", stem, Self::now_unix_ms())
            });
        Self::normalize_primer_design_report_id(&candidate).unwrap_or_else(|_| {
            format!(
                "primer_report_{}_{}",
                template.replace(' ', "_"),
                Self::now_unix_ms()
            )
        })
    }

    pub(super) fn render_restriction_cloning_pcr_handoff_report_id(
        primer_report_id: &str,
        pair_rank: usize,
        destination_vector_seq_id: &str,
        mode: RestrictionCloningPcrHandoffMode,
        forward_enzyme: &str,
        reverse_enzyme: &str,
        forward_leader_5prime: &str,
        reverse_leader_5prime: &str,
    ) -> String {
        let sanitize = |text: &str| {
            text.chars()
                .map(|ch| {
                    if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                        ch
                    } else {
                        '_'
                    }
                })
                .collect::<String>()
        };
        let mut hasher = DefaultHasher::new();
        primer_report_id.hash(&mut hasher);
        pair_rank.hash(&mut hasher);
        destination_vector_seq_id.hash(&mut hasher);
        mode.as_str().hash(&mut hasher);
        forward_enzyme.hash(&mut hasher);
        reverse_enzyme.hash(&mut hasher);
        forward_leader_5prime.hash(&mut hasher);
        reverse_leader_5prime.hash(&mut hasher);
        let digest = hasher.finish();
        let base = format!(
            "{}_r{:02}_{}_{}_{}_{}",
            sanitize(primer_report_id),
            pair_rank.max(1),
            sanitize(destination_vector_seq_id),
            sanitize(forward_enzyme),
            sanitize(reverse_enzyme),
            mode.as_str()
        );
        let candidate = format!("{base}_{digest:08x}");
        Self::normalize_primer_design_report_id(&candidate)
            .unwrap_or_else(|_| format!("restriction_clone_handoff_{digest:08x}"))
    }

    pub(super) fn annotate_primer_record_heuristics(
        mut record: PrimerDesignPrimerRecord,
        anneal_length_bp: usize,
    ) -> PrimerDesignPrimerRecord {
        let metrics = Self::compute_primer_heuristic_metrics(record.sequence.as_bytes());
        record.length_bp = metrics.length_bp;
        record.anneal_length_bp = anneal_length_bp;
        record.non_annealing_5prime_tail_bp = metrics.length_bp.saturating_sub(anneal_length_bp);
        record.three_prime_base = if metrics.length_bp == 0 {
            "N".to_string()
        } else {
            (metrics.three_prime_base as char).to_string()
        };
        record.three_prime_gc_clamp = metrics.three_prime_gc_clamp;
        record.longest_homopolymer_run_bp = metrics.longest_homopolymer_run_bp;
        record.self_complementary_run_bp = metrics.self_complementary_run_bp;
        record
    }

    pub(super) fn sort_and_rank_primer_design_pairs(
        pairs: &mut Vec<PrimerDesignPairRecord>,
        max_pairs: usize,
    ) {
        pairs.sort_by(|a, b| {
            b.score
                .total_cmp(&a.score)
                .then(a.forward.start_0based.cmp(&b.forward.start_0based))
                .then(a.reverse.start_0based.cmp(&b.reverse.start_0based))
                .then(a.amplicon_length_bp.cmp(&b.amplicon_length_bp))
                .then(a.forward.sequence.cmp(&b.forward.sequence))
                .then(a.reverse.sequence.cmp(&b.reverse.sequence))
        });
        if pairs.len() > max_pairs {
            pairs.truncate(max_pairs);
        }
        for (idx, pair) in pairs.iter_mut().enumerate() {
            pair.rank = idx + 1;
        }
    }

    pub(super) fn build_primer_design_pair_record(
        forward: PrimerDesignPrimerRecord,
        reverse: PrimerDesignPrimerRecord,
        roi_start_0based: usize,
        roi_end_0based: usize,
        min_amplicon_bp: usize,
        max_amplicon_bp: usize,
        max_tm_delta_c: f64,
        target_amplicon_bp: usize,
    ) -> Option<PrimerDesignPairRecord> {
        if reverse.end_0based_exclusive <= forward.start_0based {
            return None;
        }
        let forward_anneal_len = forward
            .end_0based_exclusive
            .saturating_sub(forward.start_0based);
        let reverse_anneal_len = reverse
            .end_0based_exclusive
            .saturating_sub(reverse.start_0based);
        let forward = Self::annotate_primer_record_heuristics(forward, forward_anneal_len);
        let reverse = Self::annotate_primer_record_heuristics(reverse, reverse_anneal_len);
        let forward_metrics = Self::compute_primer_heuristic_metrics(forward.sequence.as_bytes());
        let reverse_metrics = Self::compute_primer_heuristic_metrics(reverse.sequence.as_bytes());
        let dimer_metrics = Self::compute_primer_pair_dimer_metrics(
            forward.sequence.as_bytes(),
            reverse.sequence.as_bytes(),
        );
        let amplicon_start = forward.start_0based;
        let amplicon_end = reverse.end_0based_exclusive;
        let amplicon_length_bp = amplicon_end.saturating_sub(amplicon_start);
        let roi_covered = amplicon_start <= roi_start_0based && amplicon_end >= roi_end_0based;
        let amplicon_size_ok =
            amplicon_length_bp >= min_amplicon_bp && amplicon_length_bp <= max_amplicon_bp;
        let tm_delta_c = (forward.tm_c - reverse.tm_c).abs();
        let tm_delta_ok = tm_delta_c <= max_tm_delta_c;
        let length_penalty = amplicon_length_bp.abs_diff(target_amplicon_bp) as f64;
        let primer_length_penalty = Self::preferred_primer_length_penalty(forward.anneal_length_bp)
            + Self::preferred_primer_length_penalty(reverse.anneal_length_bp);
        let hit_penalty =
            (forward.anneal_hits.saturating_sub(1) + reverse.anneal_hits.saturating_sub(1)) as f64;
        let secondary_penalty = Self::primer_secondary_structure_penalty(forward_metrics)
            + Self::primer_secondary_structure_penalty(reverse_metrics);
        let dimer_penalty = Self::primer_pair_dimer_penalty(dimer_metrics);
        let gc_clamp_bonus = if forward_metrics.three_prime_gc_clamp {
            8.0
        } else {
            -8.0
        } + if reverse_metrics.three_prime_gc_clamp {
            8.0
        } else {
            -8.0
        };
        let score = 1000.0
            - (tm_delta_c * 20.0)
            - (length_penalty * 0.1)
            - (hit_penalty * 10.0)
            - (primer_length_penalty * 6.0)
            - secondary_penalty
            - dimer_penalty
            + gc_clamp_bonus;
        let forward_secondary_ok = forward_metrics.longest_homopolymer_run_bp
            <= PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP
            && forward_metrics.self_complementary_run_bp
                <= PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP;
        let reverse_secondary_ok = reverse_metrics.longest_homopolymer_run_bp
            <= PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP
            && reverse_metrics.self_complementary_run_bp
                <= PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP;
        let primer_pair_dimer_risk_low = dimer_metrics.max_complementary_run_bp
            <= PRIMER_RECOMMENDED_MAX_PAIR_DIMER_RUN_BP
            && dimer_metrics.max_3prime_complementary_run_bp
                <= PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP;
        Some(PrimerDesignPairRecord {
            rank: 0,
            score,
            forward,
            reverse,
            amplicon_start_0based: amplicon_start,
            amplicon_end_0based_exclusive: amplicon_end,
            amplicon_length_bp,
            tm_delta_c,
            primer_pair_complementary_run_bp: dimer_metrics.max_complementary_run_bp,
            primer_pair_3prime_complementary_run_bp: dimer_metrics.max_3prime_complementary_run_bp,
            rule_flags: PrimerDesignPairRuleFlags {
                roi_covered,
                amplicon_size_in_range: amplicon_size_ok,
                tm_delta_in_range: tm_delta_ok,
                forward_secondary_structure_ok: forward_secondary_ok,
                reverse_secondary_structure_ok: reverse_secondary_ok,
                primer_pair_dimer_risk_low,
                forward_three_prime_gc_clamp: forward_metrics.three_prime_gc_clamp,
                reverse_three_prime_gc_clamp: reverse_metrics.three_prime_gc_clamp,
            },
        })
    }

    pub(super) fn build_tailed_amplicon_sequence_from_primer_pair(
        template_seq: &str,
        pair: &PrimerDesignPairRecord,
    ) -> Result<String, EngineError> {
        let template_len = template_seq.len();
        if pair.forward.anneal_length_bp == 0 || pair.reverse.anneal_length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Primer pair records must carry non-zero anneal lengths".to_string(),

                cause_chain: vec![],
            });
        }
        if pair.forward.start_0based >= template_len
            || pair.reverse.start_0based > template_len
            || pair.forward.end_0based_exclusive > template_len
            || pair.reverse.end_0based_exclusive > template_len
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Primer pair coordinates are outside template bounds (template_len={template_len})"
                ),

                cause_chain: vec![],
            });
        }
        let forward_anneal_end = pair
            .forward
            .start_0based
            .checked_add(pair.forward.anneal_length_bp)
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: "Forward anneal geometry overflows template coordinates".to_string(),

                cause_chain: vec![],
            })?;
        if forward_anneal_end > pair.reverse.start_0based {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Primer pair anneal geometry is inconsistent: forward_anneal_end={} exceeds reverse_start={}",
                    forward_anneal_end, pair.reverse.start_0based
                ),

                cause_chain: vec![],
            });
        }
        let interior = &template_seq[forward_anneal_end..pair.reverse.start_0based];
        let reverse_full_rc = Self::reverse_complement(&pair.reverse.sequence);
        Ok(format!(
            "{}{}{}",
            pair.forward.sequence, interior, reverse_full_rc
        ))
    }

    pub(super) fn primer_pair_matches_constraints(
        template: &[u8],
        pair: &PrimerDesignPairRecord,
        roi_start_0based: usize,
        roi_end_0based: usize,
        constraints: &NormalizedPrimerPairConstraints,
    ) -> bool {
        if constraints.require_roi_flanking {
            let flanking = pair.forward.end_0based_exclusive <= roi_start_0based
                && pair.reverse.start_0based >= roi_end_0based;
            if !flanking {
                return false;
            }
        }
        if let Some(expected_start) = constraints.fixed_amplicon_start_0based
            && pair.amplicon_start_0based != expected_start
        {
            return false;
        }
        if let Some(expected_end) = constraints.fixed_amplicon_end_0based_exclusive
            && pair.amplicon_end_0based_exclusive != expected_end
        {
            return false;
        }
        if pair.amplicon_end_0based_exclusive > template.len()
            || pair.amplicon_start_0based >= pair.amplicon_end_0based_exclusive
        {
            return false;
        }
        let amplicon = &template[pair.amplicon_start_0based..pair.amplicon_end_0based_exclusive];
        for motif in &constraints.required_amplicon_motifs {
            if !Self::contains_iupac_pattern(amplicon, motif.as_bytes()) {
                return false;
            }
        }
        for motif in &constraints.forbidden_amplicon_motifs {
            if Self::contains_iupac_pattern(amplicon, motif.as_bytes()) {
                return false;
            }
        }
        true
    }

    pub(super) fn primer_pair_heuristic_advisories(pair: &PrimerDesignPairRecord) -> Vec<String> {
        let mut issues: Vec<String> = Vec::new();
        if !pair.rule_flags.forward_three_prime_gc_clamp {
            issues.push("forward primer lacks 3' GC clamp".to_string());
        }
        if !pair.rule_flags.reverse_three_prime_gc_clamp {
            issues.push("reverse primer lacks 3' GC clamp".to_string());
        }
        if !pair.rule_flags.forward_secondary_structure_ok {
            issues.push(format!(
                "forward primer may form secondary structure (homopolymer_run={}, self_complement_run={})",
                pair.forward.longest_homopolymer_run_bp, pair.forward.self_complementary_run_bp
            ));
        }
        if !pair.rule_flags.reverse_secondary_structure_ok {
            issues.push(format!(
                "reverse primer may form secondary structure (homopolymer_run={}, self_complement_run={})",
                pair.reverse.longest_homopolymer_run_bp, pair.reverse.self_complementary_run_bp
            ));
        }
        if !pair.rule_flags.primer_pair_dimer_risk_low {
            issues.push(format!(
                "primer-pair dimer risk elevated (max_complement_run={}, max_3prime_complement_run={})",
                pair.primer_pair_complementary_run_bp, pair.primer_pair_3prime_complementary_run_bp
            ));
        }
        issues
    }

    fn primer_design_cancelled_error(context: &str) -> EngineError {
        EngineError {
            code: ErrorCode::Internal,
            message: format!("Primer design cancelled during {context}"),
            cause_chain: vec![],
        }
    }

    fn emit_primer_design_progress(
        progress_context: Option<&PrimerDesignProgressContext<'_>>,
        on_progress: &mut dyn FnMut(PrimerDesignProgress) -> bool,
        stage: &str,
        detail: impl Into<String>,
        done: bool,
        forward_candidate_count: Option<usize>,
        reverse_candidate_count: Option<usize>,
        probe_candidate_count: Option<usize>,
        pair_candidate_combinations: Option<usize>,
        pair_evaluated: Option<usize>,
        pair_evaluation_limit: Option<usize>,
        pair_evaluation_limited: Option<bool>,
        accepted_pair_count: Option<usize>,
        assay_candidate_combinations: Option<usize>,
        assays_evaluated: Option<usize>,
        accepted_assay_count: Option<usize>,
    ) -> Result<(), EngineError> {
        let Some(progress_context) = progress_context else {
            return Ok(());
        };
        let progress = PrimerDesignProgress {
            seq_id: progress_context.seq_id.to_string(),
            design_kind: progress_context.design_kind.to_string(),
            backend_requested: progress_context.backend_requested.to_string(),
            backend_used: progress_context.backend_used.to_string(),
            stage: stage.to_string(),
            detail: detail.into(),
            roi_start_0based: progress_context.roi_start_0based,
            roi_end_0based_exclusive: progress_context.roi_end_0based_exclusive,
            forward_candidate_count,
            reverse_candidate_count,
            probe_candidate_count,
            pair_candidate_combinations,
            pair_evaluated,
            pair_evaluation_limit,
            pair_evaluation_limited,
            accepted_pair_count,
            assay_candidate_combinations,
            assays_evaluated,
            accepted_assay_count,
            max_output: progress_context.max_output,
            done,
        };
        if on_progress(progress) {
            Ok(())
        } else {
            Err(Self::primer_design_cancelled_error(stage))
        }
    }

    pub(super) fn design_primer_pairs_internal_core(
        template_bytes: &[u8],
        roi_start_0based: usize,
        roi_end_0based: usize,
        forward: &PrimerDesignSideConstraint,
        forward_sequence_constraints: &NormalizedPrimerSideSequenceConstraints,
        reverse: &PrimerDesignSideConstraint,
        reverse_sequence_constraints: &NormalizedPrimerSideSequenceConstraints,
        pair_constraints: &NormalizedPrimerPairConstraints,
        min_amplicon_bp: usize,
        max_amplicon_bp: usize,
        max_tm_delta_c: f64,
        max_pairs: usize,
        progress_context: Option<&PrimerDesignProgressContext<'_>>,
        on_progress: &mut dyn FnMut(PrimerDesignProgress) -> bool,
    ) -> Result<(Vec<PrimerDesignPairRecord>, PrimerDesignRejectionSummary), EngineError> {
        let mut rejection_summary = PrimerDesignRejectionSummary::default();
        let forward_candidates = Self::generate_primer_side_candidates(
            template_bytes,
            forward,
            forward_sequence_constraints,
            false,
            &mut rejection_summary,
        );
        Self::emit_primer_design_progress(
            progress_context,
            on_progress,
            "forward_candidates",
            format!(
                "Generated {} forward candidate primer(s)",
                forward_candidates.len()
            ),
            false,
            Some(forward_candidates.len()),
            None,
            None,
            None,
            Some(0),
            None,
            None,
            Some(0),
            None,
            None,
            None,
        )?;
        let reverse_candidates = Self::generate_primer_side_candidates(
            template_bytes,
            reverse,
            reverse_sequence_constraints,
            true,
            &mut rejection_summary,
        );
        let total_candidate_combinations = forward_candidates
            .len()
            .saturating_mul(reverse_candidates.len());
        let pair_evaluation_limit = max_pairs
            .saturating_mul(1_000)
            .clamp(max_pairs, PRIMER_INTERNAL_MAX_PAIR_EVALUATIONS);
        Self::emit_primer_design_progress(
            progress_context,
            on_progress,
            "reverse_candidates",
            format!(
                "Generated {} reverse candidate primer(s)",
                reverse_candidates.len()
            ),
            false,
            Some(forward_candidates.len()),
            Some(reverse_candidates.len()),
            None,
            Some(total_candidate_combinations),
            Some(0),
            Some(pair_evaluation_limit),
            Some(false),
            Some(0),
            None,
            None,
            None,
        )?;
        let mut pairs: Vec<PrimerDesignPairRecord> = vec![];
        let target_amplicon_bp = (min_amplicon_bp + max_amplicon_bp) / 2;
        let pair_progress_stride = (pair_evaluation_limit / 100).max(1);
        let mut pair_evaluations = 0usize;
        let mut pair_evaluation_limited = false;
        Self::emit_primer_design_progress(
            progress_context,
            on_progress,
            "pair_search",
            format!(
                "Scanning up to {} candidate primer-pair combination(s)",
                pair_evaluation_limit
            ),
            false,
            Some(forward_candidates.len()),
            Some(reverse_candidates.len()),
            None,
            Some(total_candidate_combinations),
            Some(0),
            Some(pair_evaluation_limit),
            Some(false),
            Some(0),
            None,
            None,
            None,
        )?;
        'pair_search: for fwd in &forward_candidates {
            for rev in &reverse_candidates {
                if pair_evaluations >= pair_evaluation_limit {
                    pair_evaluation_limited = true;
                    break 'pair_search;
                }
                pair_evaluations = pair_evaluations.saturating_add(1);
                let Some(pair) = Self::build_primer_design_pair_record(
                    PrimerDesignPrimerRecord {
                        sequence: fwd.sequence.clone(),
                        start_0based: fwd.start_0based,
                        end_0based_exclusive: fwd.end_0based_exclusive,
                        tm_c: fwd.tm_c,
                        gc_fraction: fwd.gc_fraction,
                        anneal_hits: fwd.anneal_hits,
                        ..PrimerDesignPrimerRecord::default()
                    },
                    PrimerDesignPrimerRecord {
                        sequence: rev.sequence.clone(),
                        start_0based: rev.start_0based,
                        end_0based_exclusive: rev.end_0based_exclusive,
                        tm_c: rev.tm_c,
                        gc_fraction: rev.gc_fraction,
                        anneal_hits: rev.anneal_hits,
                        ..PrimerDesignPrimerRecord::default()
                    },
                    roi_start_0based,
                    roi_end_0based,
                    min_amplicon_bp,
                    max_amplicon_bp,
                    max_tm_delta_c,
                    target_amplicon_bp,
                ) else {
                    rejection_summary.amplicon_or_roi_failure =
                        rejection_summary.amplicon_or_roi_failure.saturating_add(1);
                    if pair_evaluations.is_multiple_of(pair_progress_stride) {
                        Self::emit_primer_design_progress(
                            progress_context,
                            on_progress,
                            "pair_search",
                            format!(
                                "Evaluated {} pair candidate combination(s)",
                                pair_evaluations
                            ),
                            false,
                            Some(forward_candidates.len()),
                            Some(reverse_candidates.len()),
                            None,
                            Some(total_candidate_combinations),
                            Some(pair_evaluations),
                            Some(pair_evaluation_limit),
                            Some(false),
                            Some(pairs.len()),
                            None,
                            None,
                            None,
                        )?;
                    }
                    continue;
                };
                if !(pair.rule_flags.amplicon_size_in_range
                    && pair.rule_flags.roi_covered
                    && pair.rule_flags.tm_delta_in_range)
                {
                    rejection_summary.amplicon_or_roi_failure =
                        rejection_summary.amplicon_or_roi_failure.saturating_add(1);
                    if pair_evaluations.is_multiple_of(pair_progress_stride) {
                        Self::emit_primer_design_progress(
                            progress_context,
                            on_progress,
                            "pair_search",
                            format!(
                                "Evaluated {} pair candidate combination(s)",
                                pair_evaluations
                            ),
                            false,
                            Some(forward_candidates.len()),
                            Some(reverse_candidates.len()),
                            None,
                            Some(total_candidate_combinations),
                            Some(pair_evaluations),
                            Some(pair_evaluation_limit),
                            Some(false),
                            Some(pairs.len()),
                            None,
                            None,
                            None,
                        )?;
                    }
                    continue;
                }
                if !Self::primer_pair_matches_constraints(
                    template_bytes,
                    &pair,
                    roi_start_0based,
                    roi_end_0based,
                    pair_constraints,
                ) {
                    rejection_summary.pair_constraint_failure =
                        rejection_summary.pair_constraint_failure.saturating_add(1);
                    if pair_evaluations.is_multiple_of(pair_progress_stride) {
                        Self::emit_primer_design_progress(
                            progress_context,
                            on_progress,
                            "pair_search",
                            format!(
                                "Evaluated {} pair candidate combination(s)",
                                pair_evaluations
                            ),
                            false,
                            Some(forward_candidates.len()),
                            Some(reverse_candidates.len()),
                            None,
                            Some(total_candidate_combinations),
                            Some(pair_evaluations),
                            Some(pair_evaluation_limit),
                            Some(false),
                            Some(pairs.len()),
                            None,
                            None,
                            None,
                        )?;
                    }
                    continue;
                }
                pairs.push(pair);
                if pair_evaluations.is_multiple_of(pair_progress_stride) {
                    Self::emit_primer_design_progress(
                        progress_context,
                        on_progress,
                        "pair_search",
                        format!(
                            "Evaluated {} pair candidate combination(s)",
                            pair_evaluations
                        ),
                        false,
                        Some(forward_candidates.len()),
                        Some(reverse_candidates.len()),
                        None,
                        Some(total_candidate_combinations),
                        Some(pair_evaluations),
                        Some(pair_evaluation_limit),
                        Some(false),
                        Some(pairs.len()),
                        None,
                        None,
                        None,
                    )?;
                }
            }
        }
        if pair_evaluation_limited {
            rejection_summary.pair_evaluation_limit_skipped =
                total_candidate_combinations.saturating_sub(pair_evaluations);
        }
        Self::sort_and_rank_primer_design_pairs(&mut pairs, max_pairs);
        let operation_done = progress_context
            .map(|context| context.design_kind == "primer_pairs")
            .unwrap_or(true);
        Self::emit_primer_design_progress(
            progress_context,
            on_progress,
            "pair_search_complete",
            if pair_evaluation_limited {
                format!(
                    "Ranked {} primer pair(s) after reaching the evaluation limit of {}",
                    pairs.len(),
                    pair_evaluation_limit
                )
            } else {
                format!("Ranked {} primer pair(s)", pairs.len())
            },
            operation_done,
            Some(forward_candidates.len()),
            Some(reverse_candidates.len()),
            None,
            Some(total_candidate_combinations),
            Some(pair_evaluations),
            Some(pair_evaluation_limit),
            Some(pair_evaluation_limited),
            Some(pairs.len()),
            None,
            None,
            None,
        )?;
        Ok((pairs, rejection_summary))
    }

    #[allow(dead_code)]
    pub(super) fn design_primer_pairs_internal(
        template_bytes: &[u8],
        roi_start_0based: usize,
        roi_end_0based: usize,
        forward: &PrimerDesignSideConstraint,
        forward_sequence_constraints: &NormalizedPrimerSideSequenceConstraints,
        reverse: &PrimerDesignSideConstraint,
        reverse_sequence_constraints: &NormalizedPrimerSideSequenceConstraints,
        pair_constraints: &NormalizedPrimerPairConstraints,
        min_amplicon_bp: usize,
        max_amplicon_bp: usize,
        max_tm_delta_c: f64,
        max_pairs: usize,
    ) -> (Vec<PrimerDesignPairRecord>, PrimerDesignRejectionSummary) {
        let mut noop = |_progress: PrimerDesignProgress| true;
        Self::design_primer_pairs_internal_core(
            template_bytes,
            roi_start_0based,
            roi_end_0based,
            forward,
            forward_sequence_constraints,
            reverse,
            reverse_sequence_constraints,
            pair_constraints,
            min_amplicon_bp,
            max_amplicon_bp,
            max_tm_delta_c,
            max_pairs,
            None,
            &mut noop,
        )
        .expect("noop primer progress callback cannot cancel")
    }

    pub(super) fn sort_and_rank_qpcr_assays(assays: &mut Vec<QpcrAssayRecord>, max_assays: usize) {
        assays.sort_by(|a, b| {
            b.score
                .total_cmp(&a.score)
                .then(a.amplicon_start_0based.cmp(&b.amplicon_start_0based))
                .then(
                    a.amplicon_end_0based_exclusive
                        .cmp(&b.amplicon_end_0based_exclusive),
                )
                .then(a.probe.start_0based.cmp(&b.probe.start_0based))
                .then(a.forward.sequence.cmp(&b.forward.sequence))
                .then(a.reverse.sequence.cmp(&b.reverse.sequence))
                .then(a.probe.sequence.cmp(&b.probe.sequence))
        });
        if assays.len() > max_assays {
            assays.truncate(max_assays);
        }
        for (idx, assay) in assays.iter_mut().enumerate() {
            assay.rank = idx + 1;
        }
    }

    #[allow(dead_code)]
    pub(super) fn design_qpcr_assays_from_pairs(
        template_bytes: &[u8],
        roi_start_0based: usize,
        roi_end_0based: usize,
        probe: &PrimerDesignSideConstraint,
        probe_sequence_constraints: &NormalizedPrimerSideSequenceConstraints,
        max_probe_tm_delta_c: f64,
        max_assays: usize,
        pairs: Vec<PrimerDesignPairRecord>,
        pair_rejections: PrimerDesignRejectionSummary,
    ) -> (Vec<QpcrAssayRecord>, QpcrDesignRejectionSummary) {
        let mut noop = |_progress: PrimerDesignProgress| true;
        Self::design_qpcr_assays_from_pairs_core(
            template_bytes,
            roi_start_0based,
            roi_end_0based,
            probe,
            probe_sequence_constraints,
            max_probe_tm_delta_c,
            max_assays,
            pairs,
            pair_rejections,
            None,
            &mut noop,
        )
        .expect("noop primer progress callback cannot cancel")
    }

    pub(super) fn design_qpcr_assays_from_pairs_core(
        template_bytes: &[u8],
        roi_start_0based: usize,
        roi_end_0based: usize,
        probe: &PrimerDesignSideConstraint,
        probe_sequence_constraints: &NormalizedPrimerSideSequenceConstraints,
        max_probe_tm_delta_c: f64,
        max_assays: usize,
        pairs: Vec<PrimerDesignPairRecord>,
        pair_rejections: PrimerDesignRejectionSummary,
        progress_context: Option<&PrimerDesignProgressContext<'_>>,
        on_progress: &mut dyn FnMut(PrimerDesignProgress) -> bool,
    ) -> Result<(Vec<QpcrAssayRecord>, QpcrDesignRejectionSummary), EngineError> {
        let mut rejection = QpcrDesignRejectionSummary {
            primer_pair: pair_rejections,
            ..QpcrDesignRejectionSummary::default()
        };
        let mut probe_candidate_rejections = PrimerDesignRejectionSummary::default();
        let probe_candidates = Self::generate_primer_side_candidates(
            template_bytes,
            probe,
            probe_sequence_constraints,
            false,
            &mut probe_candidate_rejections,
        );
        rejection.probe_out_of_window = probe_candidate_rejections.out_of_window;
        rejection.probe_gc_or_tm_out_of_bounds = probe_candidate_rejections.gc_or_tm_out_of_bounds;
        rejection.probe_non_unique_anneal = probe_candidate_rejections.non_unique_anneal;
        let assay_candidate_combinations = pairs.len().saturating_mul(probe_candidates.len());
        Self::emit_primer_design_progress(
            progress_context,
            on_progress,
            "probe_candidates",
            format!(
                "Generated {} probe candidate oligo(s)",
                probe_candidates.len()
            ),
            false,
            None,
            None,
            Some(probe_candidates.len()),
            None,
            None,
            None,
            None,
            None,
            Some(assay_candidate_combinations),
            Some(0),
            Some(0),
        )?;

        let mut assays: Vec<QpcrAssayRecord> = vec![];
        let assay_progress_stride = (assay_candidate_combinations / 100).max(1);
        let mut assays_evaluated = 0usize;
        Self::emit_primer_design_progress(
            progress_context,
            on_progress,
            "assay_search",
            format!(
                "Scanning {} pair/probe assay combination(s)",
                assay_candidate_combinations
            ),
            false,
            None,
            None,
            Some(probe_candidates.len()),
            None,
            None,
            None,
            None,
            None,
            Some(assay_candidate_combinations),
            Some(0),
            Some(0),
        )?;
        for pair in pairs {
            let amplicon_mid =
                (pair.amplicon_start_0based + pair.amplicon_end_0based_exclusive) / 2;
            for probe_candidate in &probe_candidates {
                assays_evaluated = assays_evaluated.saturating_add(1);
                let probe_start = probe_candidate.start_0based;
                let probe_end = probe_candidate.end_0based_exclusive;
                let probe_inside_amplicon = probe_start >= pair.forward.end_0based_exclusive
                    && probe_end <= pair.reverse.start_0based;
                if !probe_inside_amplicon {
                    rejection.probe_or_assay_failure =
                        rejection.probe_or_assay_failure.saturating_add(1);
                    if assays_evaluated.is_multiple_of(assay_progress_stride) {
                        Self::emit_primer_design_progress(
                            progress_context,
                            on_progress,
                            "assay_search",
                            format!("Evaluated {} assay combination(s)", assays_evaluated),
                            false,
                            None,
                            None,
                            Some(probe_candidates.len()),
                            None,
                            None,
                            None,
                            None,
                            None,
                            Some(assay_candidate_combinations),
                            Some(assays_evaluated),
                            Some(assays.len()),
                        )?;
                    }
                    continue;
                }
                let primer_mean_tm_c = (pair.forward.tm_c + pair.reverse.tm_c) / 2.0;
                let probe_tm_offset_c = probe_candidate.tm_c - primer_mean_tm_c;
                let probe_tm_delta_c = probe_tm_offset_c.abs();
                let probe_tm_ok = probe_tm_delta_c <= max_probe_tm_delta_c;
                if !probe_tm_ok {
                    rejection.probe_or_assay_failure =
                        rejection.probe_or_assay_failure.saturating_add(1);
                    if assays_evaluated.is_multiple_of(assay_progress_stride) {
                        Self::emit_primer_design_progress(
                            progress_context,
                            on_progress,
                            "assay_search",
                            format!("Evaluated {} assay combination(s)", assays_evaluated),
                            false,
                            None,
                            None,
                            Some(probe_candidates.len()),
                            None,
                            None,
                            None,
                            None,
                            None,
                            Some(assay_candidate_combinations),
                            Some(assays_evaluated),
                            Some(assays.len()),
                        )?;
                    }
                    continue;
                }
                let probe_mid = (probe_start + probe_end) / 2;
                let probe_mid_penalty = probe_mid.abs_diff(amplicon_mid) as f64;
                let probe_offset_penalty_c =
                    (probe_tm_offset_c - QPCR_PREFERRED_PROBE_TM_OFFSET_C).abs();
                let score =
                    pair.score - (probe_offset_penalty_c * 10.0) - (probe_mid_penalty * 0.05);
                let probe_record = Self::annotate_primer_record_heuristics(
                    PrimerDesignPrimerRecord {
                        sequence: probe_candidate.sequence.clone(),
                        start_0based: probe_start,
                        end_0based_exclusive: probe_end,
                        tm_c: probe_candidate.tm_c,
                        gc_fraction: probe_candidate.gc_fraction,
                        anneal_hits: probe_candidate.anneal_hits,
                        ..PrimerDesignPrimerRecord::default()
                    },
                    probe_end.saturating_sub(probe_start),
                );
                assays.push(QpcrAssayRecord {
                    rank: 0,
                    score,
                    forward: pair.forward.clone(),
                    reverse: pair.reverse.clone(),
                    probe: probe_record,
                    amplicon_start_0based: pair.amplicon_start_0based,
                    amplicon_end_0based_exclusive: pair.amplicon_end_0based_exclusive,
                    amplicon_length_bp: pair.amplicon_length_bp,
                    primer_tm_delta_c: pair.tm_delta_c,
                    probe_tm_delta_c,
                    transcript_context: None,
                    rule_flags: QpcrAssayRuleFlags {
                        roi_covered: pair.rule_flags.roi_covered,
                        amplicon_size_in_range: pair.rule_flags.amplicon_size_in_range,
                        primer_tm_delta_in_range: pair.rule_flags.tm_delta_in_range,
                        probe_inside_amplicon,
                        probe_tm_delta_in_range: probe_tm_ok,
                    },
                });
                if assays_evaluated.is_multiple_of(assay_progress_stride) {
                    Self::emit_primer_design_progress(
                        progress_context,
                        on_progress,
                        "assay_search",
                        format!("Evaluated {} assay combination(s)", assays_evaluated),
                        false,
                        None,
                        None,
                        Some(probe_candidates.len()),
                        None,
                        None,
                        None,
                        None,
                        None,
                        Some(assay_candidate_combinations),
                        Some(assays_evaluated),
                        Some(assays.len()),
                    )?;
                }
            }
        }
        let _ = (roi_start_0based, roi_end_0based); // reserved for future qPCR ROI/probe rules
        Self::sort_and_rank_qpcr_assays(&mut assays, max_assays);
        Self::emit_primer_design_progress(
            progress_context,
            on_progress,
            "assay_search_complete",
            format!("Ranked {} qPCR assay candidate(s)", assays.len()),
            true,
            None,
            None,
            Some(probe_candidates.len()),
            None,
            None,
            None,
            None,
            None,
            Some(assay_candidate_combinations),
            Some(assays_evaluated),
            Some(assays.len()),
        )?;
        Ok((assays, rejection))
    }

    pub(super) fn parse_primer3_coord_pair(
        raw: &str,
        key: &str,
    ) -> Result<(usize, usize), EngineError> {
        let (left, right) = raw.split_once(',').ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Primer3 output key '{key}' is not in start,len form"),

            cause_chain: vec![],
        })?;
        let start = left.trim().parse::<usize>().map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Primer3 output key '{key}' has invalid start '{}': {e}",
                left
            ),

            cause_chain: vec![],
        })?;
        let len = right.trim().parse::<usize>().map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Primer3 output key '{key}' has invalid len '{}': {e}",
                right
            ),

            cause_chain: vec![],
        })?;
        Ok((start, len))
    }

    pub(super) fn parse_primer3_kv_output(text: &str) -> HashMap<String, String> {
        let mut map: HashMap<String, String> = HashMap::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if line == "=" {
                break;
            }
            if let Some((key, value)) = line.split_once('=') {
                map.insert(key.trim().to_string(), value.trim().to_string());
            }
        }
        map
    }

    pub(super) fn primer3_explain_summary(map: &HashMap<String, String>) -> Option<String> {
        let mut parts: Vec<String> = vec![];
        for key in [
            "PRIMER_PAIR_EXPLAIN",
            "PRIMER_LEFT_EXPLAIN",
            "PRIMER_RIGHT_EXPLAIN",
        ] {
            if let Some(value) = map
                .get(key)
                .map(|raw| raw.trim())
                .filter(|raw| !raw.is_empty())
            {
                parts.push(format!("{key}={value}"));
            }
        }
        if parts.is_empty() {
            None
        } else {
            Some(parts.join(" | "))
        }
    }

    pub(super) fn first_nonempty_utf8_line(bytes: &[u8]) -> Option<String> {
        String::from_utf8_lossy(bytes)
            .lines()
            .map(str::trim)
            .find(|line| !line.is_empty())
            .map(|line| line.to_string())
    }

    fn primer3_preflight_display_path(path: &Path) -> String {
        std::fs::canonicalize(path)
            .unwrap_or_else(|_| path.to_path_buf())
            .display()
            .to_string()
    }

    fn primer3_preflight_resolved_path(executable: &str) -> Option<String> {
        let trimmed = executable.trim();
        if trimmed.is_empty() {
            return None;
        }
        let candidate = Path::new(trimmed);
        if candidate.is_absolute() || candidate.components().count() > 1 {
            return candidate
                .exists()
                .then(|| Self::primer3_preflight_display_path(candidate));
        }
        let path_env = std::env::var_os("PATH")?;
        for dir in std::env::split_paths(&path_env) {
            let joined = dir.join(trimmed);
            if joined.exists() {
                return Some(Self::primer3_preflight_display_path(&joined));
            }
        }
        None
    }

    pub(super) fn probe_primer3_executable_status(executable: &str) -> Primer3PreflightReport {
        let mut report = Primer3PreflightReport {
            executable: executable.to_string(),
            resolved_path: Self::primer3_preflight_resolved_path(executable),
            ..Primer3PreflightReport::default()
        };
        let started = Instant::now();
        match Command::new(executable).arg("--version").output() {
            Ok(output) => {
                report.reachable = true;
                report.status_code = output.status.code();
                report.version_probe_ok = output.status.success();
                let stdout_line = Self::first_nonempty_utf8_line(&output.stdout);
                let stderr_line = Self::first_nonempty_utf8_line(&output.stderr);
                report.version = stdout_line.or_else(|| stderr_line.clone());
                if !report.version_probe_ok {
                    report.detail = stderr_line.or_else(|| report.version.clone());
                }
            }
            Err(err) => {
                report.error = Some(err.to_string());
            }
        }
        report.probe_time_ms = started.elapsed().as_millis();
        report
    }

    pub fn primer3_preflight_report(
        &self,
        backend_override: Option<PrimerDesignBackend>,
        primer3_executable_override: Option<&str>,
    ) -> Primer3PreflightReport {
        let backend = backend_override.unwrap_or(self.state.parameters.primer_design_backend);
        let configured_executable = primer3_executable_override
            .map(str::trim)
            .unwrap_or_else(|| self.state.parameters.primer3_executable.trim());
        let configured_executable =
            (!configured_executable.is_empty()).then(|| configured_executable.to_string());
        let executable = configured_executable.as_deref().unwrap_or("primer3_core");
        let mut report = Self::probe_primer3_executable_status(executable);
        report.backend = backend.as_str().to_string();
        report.configured_executable = configured_executable;
        report.used_default_executable = report.configured_executable.is_none();
        report.working_directory = std::env::current_dir()
            .ok()
            .map(|path| path.display().to_string());
        report
    }

    pub(super) fn probe_primer3_version(executable: &str) -> Option<String> {
        Self::probe_primer3_executable_status(executable).version
    }

    pub(super) fn design_primer_pairs_primer3(
        template_seq: &str,
        roi_start_0based: usize,
        roi_end_0based: usize,
        forward: &PrimerDesignSideConstraint,
        forward_sequence_constraints: &NormalizedPrimerSideSequenceConstraints,
        reverse: &PrimerDesignSideConstraint,
        reverse_sequence_constraints: &NormalizedPrimerSideSequenceConstraints,
        pair_constraints: &NormalizedPrimerPairConstraints,
        min_amplicon_bp: usize,
        max_amplicon_bp: usize,
        max_tm_delta_c: f64,
        max_pairs: usize,
        primer3_executable: &str,
    ) -> Result<
        (
            Vec<PrimerDesignPairRecord>,
            PrimerDesignRejectionSummary,
            Option<String>,
            Option<String>,
            String,
        ),
        EngineError,
    > {
        let template_bytes = template_seq.as_bytes();
        let template_len = template_bytes.len();
        let target_len = roi_end_0based.saturating_sub(roi_start_0based);
        let min_size = forward.min_length.min(reverse.min_length);
        let max_size = forward.max_length.max(reverse.max_length);
        let min_tm = forward.min_tm_c.min(reverse.min_tm_c);
        let max_tm = forward.max_tm_c.max(reverse.max_tm_c);
        let min_gc_percent = forward.min_gc_fraction.min(reverse.min_gc_fraction) * 100.0;
        let max_gc_percent = forward.max_gc_fraction.max(reverse.max_gc_fraction) * 100.0;
        let num_return = max_pairs.saturating_mul(5).clamp(50, 1000);

        let input = format!(
            "SEQUENCE_ID=gentle_primer_design\nSEQUENCE_TEMPLATE={template_seq}\nSEQUENCE_TARGET={roi_start_0based},{target_len}\nPRIMER_TASK=generic\nPRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_RIGHT_PRIMER=1\nPRIMER_PICK_INTERNAL_OLIGO=0\nPRIMER_MIN_SIZE={min_size}\nPRIMER_MAX_SIZE={max_size}\nPRIMER_MIN_TM={min_tm:.3}\nPRIMER_MAX_TM={max_tm:.3}\nPRIMER_MIN_GC={min_gc_percent:.3}\nPRIMER_MAX_GC={max_gc_percent:.3}\nPRIMER_PRODUCT_SIZE_RANGE={min_amplicon_bp}-{max_amplicon_bp}\nPRIMER_NUM_RETURN={num_return}\nPRIMER_EXPLAIN_FLAG=1\n=\n"
        );
        let mut child = Command::new(primer3_executable)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| EngineError {
                code: if e.kind() == std::io::ErrorKind::NotFound {
                    ErrorCode::Unsupported
                } else {
                    ErrorCode::Io
                },
                message: format!(
                    "Primer3 backend executable '{}' is not available: {e}",
                    primer3_executable
                ),

                cause_chain: vec![],
            })?;
        if let Some(stdin) = child.stdin.as_mut() {
            stdin.write_all(input.as_bytes()).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write Primer3 request to stdin: {e}"),

                cause_chain: vec![],
            })?;
        }
        let output = child.wait_with_output().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not read Primer3 response from '{}': {e}",
                primer3_executable
            ),

            cause_chain: vec![],
        })?;
        let stdout_text = String::from_utf8_lossy(&output.stdout).to_string();
        let stderr_text = String::from_utf8_lossy(&output.stderr).to_string();
        let map = Self::parse_primer3_kv_output(&stdout_text);
        let primer3_explain = Self::primer3_explain_summary(&map);

        if let Some(primer_error) = map.get("PRIMER_ERROR") {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Primer3 backend rejected request: {primer_error}"),

                cause_chain: vec![],
            });
        }
        if !output.status.success() {
            let detail = stderr_text
                .lines()
                .map(str::trim)
                .find(|line| !line.is_empty())
                .unwrap_or("no stderr detail");
            return Err(EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Primer3 backend '{}' failed (status={}): {}",
                    primer3_executable,
                    output
                        .status
                        .code()
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "signal".to_string()),
                    detail
                ),

                cause_chain: vec![],
            });
        }

        let pair_count = map
            .get("PRIMER_PAIR_NUM_RETURNED")
            .and_then(|raw| raw.parse::<usize>().ok())
            .unwrap_or(0);
        let mut rejection_summary = PrimerDesignRejectionSummary::default();
        let mut pairs: Vec<PrimerDesignPairRecord> = vec![];
        let target_amplicon_bp = (min_amplicon_bp + max_amplicon_bp) / 2;

        for idx in 0..pair_count {
            let left_key = format!("PRIMER_LEFT_{idx}");
            let right_key = format!("PRIMER_RIGHT_{idx}");
            let left_seq_key = format!("PRIMER_LEFT_{idx}_SEQUENCE");
            let right_seq_key = format!("PRIMER_RIGHT_{idx}_SEQUENCE");
            let left_tm_key = format!("PRIMER_LEFT_{idx}_TM");
            let right_tm_key = format!("PRIMER_RIGHT_{idx}_TM");
            let left_gc_key = format!("PRIMER_LEFT_{idx}_GC_PERCENT");
            let right_gc_key = format!("PRIMER_RIGHT_{idx}_GC_PERCENT");

            let (left_start, left_len) = match map.get(&left_key) {
                Some(raw) => Self::parse_primer3_coord_pair(raw, &left_key)?,
                None => continue,
            };
            let (right_3prime_0based, right_len) = match map.get(&right_key) {
                Some(raw) => Self::parse_primer3_coord_pair(raw, &right_key)?,
                None => continue,
            };
            if left_len == 0 || right_len == 0 {
                rejection_summary.out_of_window = rejection_summary.out_of_window.saturating_add(1);
                continue;
            }
            let Some(left_end) = left_start.checked_add(left_len) else {
                rejection_summary.out_of_window = rejection_summary.out_of_window.saturating_add(1);
                continue;
            };
            if left_end > template_len
                || right_3prime_0based >= template_len
                || right_3prime_0based + 1 < right_len
            {
                rejection_summary.out_of_window = rejection_summary.out_of_window.saturating_add(1);
                continue;
            }
            let right_start = right_3prime_0based + 1 - right_len;
            let right_end_exclusive = right_3prime_0based + 1;
            if right_end_exclusive <= left_start || right_end_exclusive > template_len {
                rejection_summary.amplicon_or_roi_failure =
                    rejection_summary.amplicon_or_roi_failure.saturating_add(1);
                continue;
            }

            let forward_anneal_sequence = map
                .get(&left_seq_key)
                .cloned()
                .unwrap_or_else(|| template_seq[left_start..left_end].to_string())
                .to_ascii_uppercase();
            let reverse_anneal_sequence = map
                .get(&right_seq_key)
                .cloned()
                .unwrap_or_else(|| {
                    Self::reverse_complement(&template_seq[right_start..right_end_exclusive])
                })
                .to_ascii_uppercase();
            let forward_sequence =
                if let Some(tail) = &forward_sequence_constraints.non_annealing_5prime_tail {
                    format!("{tail}{forward_anneal_sequence}")
                } else {
                    forward_anneal_sequence.clone()
                };
            let reverse_sequence =
                if let Some(tail) = &reverse_sequence_constraints.non_annealing_5prime_tail {
                    format!("{tail}{reverse_anneal_sequence}")
                } else {
                    reverse_anneal_sequence.clone()
                };
            let forward_tm = map
                .get(&left_tm_key)
                .and_then(|raw| raw.parse::<f64>().ok())
                .unwrap_or_else(|| Self::estimate_primer_tm_c(forward_anneal_sequence.as_bytes()));
            let reverse_tm = map
                .get(&right_tm_key)
                .and_then(|raw| raw.parse::<f64>().ok())
                .unwrap_or_else(|| Self::estimate_primer_tm_c(reverse_anneal_sequence.as_bytes()));
            let forward_gc = map
                .get(&left_gc_key)
                .and_then(|raw| raw.parse::<f64>().ok())
                .map(|value| (value / 100.0).clamp(0.0, 1.0))
                .unwrap_or_else(|| {
                    Self::sequence_gc_fraction(forward_anneal_sequence.as_bytes()).unwrap_or(0.0)
                });
            let reverse_gc = map
                .get(&right_gc_key)
                .and_then(|raw| raw.parse::<f64>().ok())
                .map(|value| (value / 100.0).clamp(0.0, 1.0))
                .unwrap_or_else(|| {
                    Self::sequence_gc_fraction(reverse_anneal_sequence.as_bytes()).unwrap_or(0.0)
                });

            let left_window = &template_bytes[left_start..left_end];
            let right_window = &template_bytes[right_start..right_end_exclusive];
            let forward_anneal_hits =
                Self::find_all_subsequences(template_bytes, left_window).len();
            let reverse_anneal_hits =
                Self::find_all_subsequences(template_bytes, right_window).len();

            let forward_window_ok = forward
                .location_0based
                .is_none_or(|value| value == left_start)
                && forward.start_0based.is_none_or(|value| left_start >= value)
                && forward.end_0based.is_none_or(|value| left_end <= value)
                && left_len >= forward.min_length
                && left_len <= forward.max_length;
            let reverse_window_ok = reverse
                .location_0based
                .is_none_or(|value| value == right_start)
                && reverse
                    .start_0based
                    .is_none_or(|value| right_start >= value)
                && reverse
                    .end_0based
                    .is_none_or(|value| right_end_exclusive <= value)
                && right_len >= reverse.min_length
                && right_len <= reverse.max_length;
            if !(forward_window_ok && reverse_window_ok) {
                rejection_summary.out_of_window = rejection_summary.out_of_window.saturating_add(1);
                continue;
            }
            let forward_gc_tm_ok = forward_tm >= forward.min_tm_c
                && forward_tm <= forward.max_tm_c
                && forward_gc >= forward.min_gc_fraction
                && forward_gc <= forward.max_gc_fraction;
            let reverse_gc_tm_ok = reverse_tm >= reverse.min_tm_c
                && reverse_tm <= reverse.max_tm_c
                && reverse_gc >= reverse.min_gc_fraction
                && reverse_gc <= reverse.max_gc_fraction;
            if !(forward_gc_tm_ok && reverse_gc_tm_ok) {
                rejection_summary.gc_or_tm_out_of_bounds =
                    rejection_summary.gc_or_tm_out_of_bounds.saturating_add(1);
                continue;
            }
            let forward_anneal_ok =
                forward_anneal_hits >= 1 && forward_anneal_hits <= forward.max_anneal_hits;
            let reverse_anneal_ok =
                reverse_anneal_hits >= 1 && reverse_anneal_hits <= reverse.max_anneal_hits;
            if !(forward_anneal_ok && reverse_anneal_ok) {
                rejection_summary.non_unique_anneal =
                    rejection_summary.non_unique_anneal.saturating_add(1);
                continue;
            }
            if !Self::primer_sequence_matches_side_constraints(
                forward_sequence.as_bytes(),
                forward_sequence_constraints,
            ) || !Self::primer_sequence_matches_side_constraints(
                reverse_sequence.as_bytes(),
                reverse_sequence_constraints,
            ) {
                rejection_summary.primer_constraint_failure = rejection_summary
                    .primer_constraint_failure
                    .saturating_add(1);
                continue;
            }

            let Some(pair) = Self::build_primer_design_pair_record(
                PrimerDesignPrimerRecord {
                    sequence: forward_sequence,
                    start_0based: left_start,
                    end_0based_exclusive: left_end,
                    tm_c: forward_tm,
                    gc_fraction: forward_gc,
                    anneal_hits: forward_anneal_hits,
                    ..PrimerDesignPrimerRecord::default()
                },
                PrimerDesignPrimerRecord {
                    sequence: reverse_sequence,
                    start_0based: right_start,
                    end_0based_exclusive: right_end_exclusive,
                    tm_c: reverse_tm,
                    gc_fraction: reverse_gc,
                    anneal_hits: reverse_anneal_hits,
                    ..PrimerDesignPrimerRecord::default()
                },
                roi_start_0based,
                roi_end_0based,
                min_amplicon_bp,
                max_amplicon_bp,
                max_tm_delta_c,
                target_amplicon_bp,
            ) else {
                rejection_summary.amplicon_or_roi_failure =
                    rejection_summary.amplicon_or_roi_failure.saturating_add(1);
                continue;
            };
            if !(pair.rule_flags.amplicon_size_in_range
                && pair.rule_flags.roi_covered
                && pair.rule_flags.tm_delta_in_range)
            {
                rejection_summary.amplicon_or_roi_failure =
                    rejection_summary.amplicon_or_roi_failure.saturating_add(1);
                continue;
            }
            if !Self::primer_pair_matches_constraints(
                template_bytes,
                &pair,
                roi_start_0based,
                roi_end_0based,
                pair_constraints,
            ) {
                rejection_summary.pair_constraint_failure =
                    rejection_summary.pair_constraint_failure.saturating_add(1);
                continue;
            }
            pairs.push(pair);
        }
        Self::sort_and_rank_primer_design_pairs(&mut pairs, max_pairs);
        Ok((
            pairs,
            rejection_summary,
            Self::probe_primer3_version(primer3_executable),
            primer3_explain,
            input,
        ))
    }

    pub(super) fn right_end_overhangs(dna: &DNAsequence) -> Vec<Vec<u8>> {
        let mut ret = Vec::new();
        if !dna.overhang().forward_3.is_empty() {
            ret.push(dna.overhang().forward_3.clone());
        }
        if !dna.overhang().reverse_5.is_empty() {
            ret.push(dna.overhang().reverse_5.clone());
        }
        ret
    }

    pub(super) fn left_end_overhangs(dna: &DNAsequence) -> Vec<Vec<u8>> {
        let mut ret = Vec::new();
        if !dna.overhang().forward_5.is_empty() {
            ret.push(dna.overhang().forward_5.clone());
        }
        if !dna.overhang().reverse_3.is_empty() {
            ret.push(dna.overhang().reverse_3.clone());
        }
        ret
    }

    pub(super) fn right_end_is_blunt(dna: &DNAsequence) -> bool {
        dna.overhang().forward_3.is_empty() && dna.overhang().reverse_5.is_empty()
    }

    pub(super) fn left_end_is_blunt(dna: &DNAsequence) -> bool {
        dna.overhang().forward_5.is_empty() && dna.overhang().reverse_3.is_empty()
    }

    pub(super) fn sticky_compatible(left: &DNAsequence, right: &DNAsequence) -> bool {
        let right_opts = Self::right_end_overhangs(left);
        let left_opts = Self::left_end_overhangs(right);
        if right_opts.is_empty() || left_opts.is_empty() {
            return false;
        }
        for r in &right_opts {
            let rc_r = Self::reverse_complement_bytes(r);
            let c_r = Self::complement_bytes(r);
            if left_opts.iter().any(|l| *l == rc_r || *l == c_r) {
                return true;
            }
        }
        false
    }

    pub(super) fn find_anneal_sites(
        template: &[u8],
        anneal: &[u8],
        max_mismatches: usize,
        require_3prime_exact_bases: usize,
        three_prime_is_window_end: bool,
    ) -> Vec<usize> {
        let mut ret = vec![];
        if anneal.is_empty() || template.len() < anneal.len() {
            return ret;
        }
        let window_len = anneal.len();
        for start in 0..=(template.len() - window_len) {
            let window = &template[start..start + window_len];
            let mismatches = window
                .iter()
                .zip(anneal.iter())
                .filter(|(a, b)| a != b)
                .count();
            if mismatches > max_mismatches {
                continue;
            }
            if require_3prime_exact_bases > 0 {
                if require_3prime_exact_bases > window_len {
                    continue;
                }
                let exact_ok = if three_prime_is_window_end {
                    let from = window_len - require_3prime_exact_bases;
                    window[from..] == anneal[from..]
                } else {
                    window[..require_3prime_exact_bases] == anneal[..require_3prime_exact_bases]
                };
                if !exact_ok {
                    continue;
                }
            }
            ret.push(start);
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn estimate_primer_tm_uses_nearest_neighbor_for_canonical_oligos() {
        let primer = b"AAAATCGATCGATCGATCGATCGATCGATC";
        let tm = GentleEngine::estimate_primer_tm_c(primer);
        let wallace = GentleEngine::estimate_primer_tm_wallace_c(primer);
        assert!(
            (60.0..75.0).contains(&tm),
            "nearest-neighbor Tm should stay in a realistic range, got {tm}"
        );
        assert!(
            (tm - wallace).abs() > 5.0,
            "nearest-neighbor estimate should materially differ from the Wallace estimate ({tm} vs {wallace})"
        );
    }

    #[test]
    fn estimate_primer_tm_matches_thermofisher_allawi_santalucia_examples() {
        let cases = [
            (b"CAGTCAGTCGATC".as_slice(), 46.85),
            (b"ATGCGTACGTAGCTAGCTA".as_slice(), 60.80),
            (b"GCGCGCGCGCGCGCGCGC".as_slice(), 82.38),
            (b"AAAATCGATCGATCGATCGATCGATCGATC".as_slice(), 68.85),
        ];
        for (primer, expected_tm) in cases {
            let observed = GentleEngine::estimate_primer_tm_c(primer);
            assert!(
                (observed - expected_tm).abs() < 0.01,
                "Tm for {} should match Thermo Fisher-style Allawi/SantaLucia calculation ({expected_tm:.2}), got {observed:.2}",
                String::from_utf8_lossy(primer)
            );
        }
    }

    #[test]
    fn estimate_primer_tm_falls_back_for_too_short_oligos() {
        let primer = b"ATGC";
        assert_eq!(
            GentleEngine::estimate_primer_tm_c(primer),
            GentleEngine::estimate_primer_tm_wallace_c(primer)
        );
    }

    #[test]
    fn estimate_primer_tm_falls_back_for_ambiguous_oligos() {
        let primer = b"ACGTNN";
        assert_eq!(
            GentleEngine::estimate_primer_tm_c(primer),
            GentleEngine::estimate_primer_tm_wallace_c(primer)
        );
    }

    #[test]
    fn primer_tm_model_description_mentions_shared_assumptions() {
        let description = GentleEngine::primer_tm_model_description();
        assert!(description.contains("Thermo Fisher"));
        assert!(description.contains("SantaLucia"));
        assert!(description.contains("Displayed Tm values"));
        assert!(description.contains("215"));
        assert!(description.contains("500"));
        assert!(description.contains("high-fidelity"));
        assert!(description.contains("fall back"));
    }
}
