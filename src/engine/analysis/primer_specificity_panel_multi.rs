//! Bounded external specificity work for a saved panel. Children are ordinary
//! single-reference panel handoffs; only the parent import persists diagnostics.

use super::primer_specificity_multi::*;
use super::*;
use std::fs;

const LIMITS: PanelSpecificityLimits = PanelSpecificityLimits {
    max_references: MULTI_MAX_REFERENCES,
    max_prepared_references: MULTI_MAX_REFERENCES,
    max_assays: 24,
    max_commands: 384,
    max_bytes_per_output: MULTI_MAX_OUTPUT_BYTES,
    max_total_output_bytes: MULTI_MAX_OUTPUT_BYTES,
    max_parsed_rows: 25_000,
};

fn invalid(message: &str) -> EngineError {
    EngineError::invalid_input(format!("panel_specificity_multi: {message}"))
}

fn digest(value: &impl Serialize) -> Result<String, EngineError> {
    serde_json::to_vec(value)
        .map(|b| sha256_prefixed_bytes(&b))
        .map_err(|e| EngineError::internal(e.to_string()))
}

fn nonclaims() -> Vec<String> {
    let panel = [
        "This standalone summary never attaches evidence to the panel or changes readiness.",
        "Admission ceilings are for synthetic evaluation only. Single-pair verification does not establish real-reference panel capacity or a memory/time bound.",
    ];
    MULTI_NONCLAIMS
        .iter()
        .copied()
        .chain(panel)
        .map(str::to_owned)
        .collect()
}

impl PanelSpecificityLimits {
    fn validate(self) -> Result<(), EngineError> {
        if !(1..=LIMITS.max_references).contains(&self.max_references)
            || !(1..=LIMITS.max_prepared_references).contains(&self.max_prepared_references)
            || !(1..=LIMITS.max_assays).contains(&self.max_assays)
            || !(1..=LIMITS.max_commands).contains(&self.max_commands)
            || self.max_bytes_per_output > LIMITS.max_bytes_per_output
            || self.max_total_output_bytes > LIMITS.max_total_output_bytes
            || self.max_parsed_rows > LIMITS.max_parsed_rows
        {
            return Err(invalid("unsupported effective limits"));
        }
        Ok(())
    }

    fn commands(self, assays: usize, prepared: usize) -> Result<usize, EngineError> {
        if assays == 0
            || assays > self.max_assays
            || prepared == 0
            || prepared > self.max_prepared_references
        {
            return Err(invalid("assay or prepared-reference limit exceeded"));
        }
        assays
            .checked_mul(prepared)
            .and_then(|n| n.checked_mul(2))
            .filter(|n| *n <= self.max_commands)
            .ok_or_else(|| invalid("command limit exceeded"))
    }
}

fn request_refs(
    request: &TranscriptAssayPanelMultiReferenceRequest,
    limits: PanelSpecificityLimits,
) -> Result<Vec<PrimerSpecificityMultiReference>, EngineError> {
    limits.validate()?;
    if request.schema != PANEL_MULTI_REQUEST_SCHEMA
        || request.panel_report_id.trim().is_empty()
        || request.references.len() > limits.max_references
    {
        return Err(invalid(
            "invalid request schema, panel ID or reference count",
        ));
    }
    let refs: Vec<_> = request
        .references
        .iter()
        .map(|r| PrimerSpecificityMultiReference {
            target_genome_id: r.target_genome_id.clone(),
            expected_index_kind: r.expected_index_kind,
            required: r.required,
            intended_target: None,
        })
        .collect();
    GentleEngine::validate_specificity_references(&request.policy, &refs)?;
    Ok(refs)
}

fn ref_summary(r: &PanelSpecificityChild) -> PanelSpecificityReferenceSummary {
    PanelSpecificityReferenceSummary {
        requested_genome_id: r.requested_genome_id.clone(),
        resolved_genome_id: r.resolved_genome_id.clone(),
        expected_index_kind: r.expected_index_kind,
        required: r.required,
        reference: r.reference.clone(),
        availability: r.availability.clone(),
        diagnostic: r.diagnostic.clone(),
    }
}

fn scientific_id(h: &TranscriptAssayPanelMultiReferenceHandoff) -> Result<String, EngineError> {
    let references: Vec<_> = h
        .references
        .iter()
        .map(|r| {
            let assays = r.child.as_ref().map(|child| {
                child
                    .assays
                    .iter()
                    .map(|a| {
                        (
                            &a.assay_id,
                            a.assay_rank,
                            &a.primer_pair_digest,
                            &a.handoff.handoff_id,
                            &a.handoff.effective_blast_options,
                        )
                    })
                    .collect::<Vec<_>>()
            });
            json!([
                r.requested_genome_id,
                r.resolved_genome_id,
                r.required,
                r.expected_index_kind,
                r.availability,
                r.reference,
                assays
            ])
        })
        .collect();
    digest(&json!([
        h.schema,
        h.panel_report_id,
        h.panel_digest,
        h.source_snapshot_sha256,
        h.request.policy,
        references,
        h.effective_limits,
        h.expected_command_count
    ]))
}

/// Runners may use either the parent's inventory or the existing child handoff files.
fn check_child_file<T: Serialize>(
    expected: &Path,
    declared: &str,
    record: &T,
) -> Result<(), EngineError> {
    if Path::new(declared) != expected {
        return Err(invalid("child file outside declared bundle"));
    }
    let expected_json =
        serde_json::to_value(record).map_err(|e| EngineError::internal(e.to_string()))?;
    let actual: serde_json::Value = specificity_multi_read_json(declared)?;
    if actual != expected_json {
        return Err(invalid("published child file differs from bound record"));
    }
    Ok(())
}

fn bundle_digest(h: &TranscriptAssayPanelMultiReferenceHandoff) -> Result<String, EngineError> {
    let mut bound = h.clone();
    bound.bundle_sha256.clear();
    digest(&bound)
}

fn commands(
    h: &TranscriptAssayPanelMultiReferenceHandoff,
) -> impl Iterator<Item = &PrimerSpecificityHandoffCommand> {
    h.references
        .iter()
        .filter_map(|r| r.child.as_ref())
        .flat_map(|c| &c.assays)
        .flat_map(|a| &a.handoff.commands)
}

fn panel_manifest(
    h: &TranscriptAssayPanelMultiReferenceHandoff,
) -> Result<PanelSpecificityExecutionManifest, EngineError> {
    h.effective_limits.validate()?;
    if h.expected_command_count > h.effective_limits.max_commands
        || h.references.len() > h.effective_limits.max_references
    {
        return Err(invalid("manifest inventory exceeds admitted limits"));
    }
    let mut ids = HashSet::new();
    let commands = commands(h)
        .map(|c| {
            if ids.len() >= h.expected_command_count || !ids.insert(&c.command_id) {
                return Err(invalid("extra or duplicate child command ID"));
            }
            Ok(PrimerSpecificityMultiExecutionCommand {
                command_id: c.command_id.clone(),
                command_sha256: GentleEngine::specificity_multi_command_digest(c)?,
                output_path: c.output_tsv_path.clone(),
                state: PrimerSpecificityMultiExecutionState::Pending,
                exit_code: None,
                output_size_bytes: None,
                output_sha256: None,
            })
        })
        .collect::<Result<Vec<_>, EngineError>>()?;
    if commands.len() != h.expected_command_count {
        return Err(invalid("expected command count mismatch"));
    }
    Ok(PanelSpecificityExecutionManifest {
        schema: PANEL_MULTI_MANIFEST_SCHEMA.into(),
        handoff_id: h.handoff_id.clone(),
        handoff_content_sha256: h.content_sha256.clone(),
        bundle_sha256: h.bundle_sha256.clone(),
        panel_digest: h.panel_digest.clone(),
        commands,
    })
}

fn manifest_rows<'a>(
    h: &TranscriptAssayPanelMultiReferenceHandoff,
    m: &'a PanelSpecificityExecutionManifest,
) -> Result<
    (
        BTreeMap<String, &'a PrimerSpecificityMultiExecutionCommand>,
        u64,
    ),
    EngineError,
> {
    let expected = panel_manifest(h)?;
    if m.schema != expected.schema
        || m.handoff_id != expected.handoff_id
        || m.handoff_content_sha256 != expected.handoff_content_sha256
        || m.bundle_sha256 != expected.bundle_sha256
        || m.panel_digest != expected.panel_digest
        || m.commands.len() != expected.commands.len()
    {
        return Err(invalid("manifest binding or inventory length mismatch"));
    }
    let by_id: BTreeMap<_, _> = expected
        .commands
        .iter()
        .map(|c| (c.command_id.as_str(), c))
        .collect();
    let mut rows = BTreeMap::new();
    let mut total = 0u64;
    for c in &m.commands {
        let expected = by_id
            .get(c.command_id.as_str())
            .ok_or_else(|| invalid("unknown manifest command"))?;
        if expected.command_sha256 != c.command_sha256
            || expected.output_path != c.output_path
            || rows.insert(c.command_id.clone(), c).is_some()
        {
            return Err(invalid("changed or duplicate manifest command"));
        }
        if c.state == PrimerSpecificityMultiExecutionState::Completed {
            let n = c
                .output_size_bytes
                .ok_or_else(|| invalid("completed output size unknown"))?;
            if c.exit_code != Some(0)
                || !c
                    .output_sha256
                    .as_deref()
                    .and_then(|s| s.strip_prefix("sha256:"))
                    .is_some_and(|s| s.len() == 64 && s.bytes().all(|c| c.is_ascii_hexdigit()))
            {
                return Err(invalid(
                    "completed output requires a successful exit and SHA-256; use failed for unsuccessful commands",
                ));
            }
            if n > h.effective_limits.max_bytes_per_output {
                return Err(invalid("per-output byte limit exceeded"));
            }
            total = total
                .checked_add(n)
                .filter(|v| *v <= h.effective_limits.max_total_output_bytes)
                .ok_or_else(|| invalid("total declared output byte limit exceeded"))?;
        } else if c.output_size_bytes.is_some() || c.output_sha256.is_some() {
            return Err(invalid(
                "only completed successful commands can claim retained output",
            ));
        }
    }
    Ok((rows, total))
}

fn assay_order(
    panel: &TranscriptAssayPanelReport,
) -> Result<Vec<&TranscriptAssayPanelAssay>, EngineError> {
    let mut ordered: Vec<_> = panel.selected_assays.iter().collect();
    let mut ids = HashSet::new();
    if ordered
        .iter()
        .any(|a| a.assay_id.trim().is_empty() || !ids.insert(a.assay_id.as_str()))
    {
        return Err(invalid("empty or duplicate assay ID"));
    }
    ordered.sort_by(|a, b| {
        a.rank
            .cmp(&b.rank)
            .then_with(|| a.assay_id.cmp(&b.assay_id))
    });
    Ok(ordered)
}

impl GentleEngine {
    fn panel_multi_snapshot(
        &self,
        panel: &TranscriptAssayPanelReport,
    ) -> Result<String, EngineError> {
        let store = self.read_primer_design_store();
        let binding = store.primer_specificity_panel_sources.get(&panel.report_id);
        self.specificity_source_snapshot(&panel.source_seq_id, binding.map(|s| &s.source_reference))
    }

    /// Preflight all declared references, then build only inside a freshly reserved directory.
    pub fn prepare_panel_multi_specificity(
        &self,
        mut request: TranscriptAssayPanelMultiReferenceRequest,
        output_dir: &str,
    ) -> Result<TranscriptAssayPanelMultiReferenceHandoff, EngineError> {
        let mut requested = request_refs(&request, LIMITS)?;
        let panel = self.get_transcript_assay_panel_report(&request.panel_report_id)?;
        let assays = assay_order(&panel)?;
        LIMITS.commands(assays.len(), 1)?;
        let panel_digest = Self::transcript_assay_panel_specificity_digest(&panel)?;
        if request
            .expected_panel_digest
            .as_ref()
            .is_some_and(|s| s != &panel_digest)
        {
            return Err(invalid("expected panel digest mismatch"));
        }
        let source_snapshot_sha256 = self.panel_multi_snapshot(&panel)?;
        let (catalog, mut references): (_, Vec<PanelSpecificityChild>) =
            Self::specificity_preflight_references(
                &mut requested,
                request.catalog_path.as_deref(),
                request.cache_dir.as_deref(),
            )?;
        let expected_command_count = LIMITS.commands(
            assays.len(),
            references.iter().filter(|r| r.reference.is_some()).count(),
        )?;
        request.catalog_path = Some(catalog);
        request.panel_report_id = panel.report_id.clone();
        for (r, resolved) in request.references.iter_mut().zip(requested) {
            r.target_genome_id = resolved.target_genome_id;
        }
        with_fresh_specificity_bundle(output_dir, |root| {
            for (i, r) in references.iter_mut().enumerate() {
                if r.reference.is_none() {
                    continue;
                }
                let dir = root.join(format!("reference_{i:02}"));
                let child = self.prepare_transcript_assay_panel_specificity_handoff(
                    &panel.report_id,
                    &r.resolved_genome_id,
                    request.policy.clone(),
                    request.catalog_path.as_deref(),
                    request.cache_dir.as_deref(),
                    &dir.to_string_lossy(),
                )?;
                if child.panel_digest != panel_digest
                    || child.resolved_target_genome_id != r.resolved_genome_id
                    || child.assays.iter().any(|a| {
                        a.handoff
                            .blast_database
                            .as_ref()
                            .and_then(Self::primer_specificity_reference_identity)
                            != r.reference
                    })
                {
                    return Err(invalid(
                        "panel or exact reference changed during preparation",
                    ));
                }
                r.child = Some(child);
            }
            let mut h = TranscriptAssayPanelMultiReferenceHandoff {
                schema: PANEL_MULTI_HANDOFF_SCHEMA.into(),
                handoff_id: String::new(),
                content_sha256: String::new(),
                bundle_sha256: String::new(),
                panel_report_id: panel.report_id.clone(),
                panel_digest,
                source_snapshot_sha256,
                request,
                effective_limits: LIMITS,
                expected_command_count,
                references,
                manifest_path: root
                    .join("execution_manifest.json")
                    .to_string_lossy()
                    .into_owned(),
                nonclaims: nonclaims(),
            };
            h.content_sha256 = scientific_id(&h)?;
            h.handoff_id = format!(
                "panel_multi_{}",
                h.content_sha256.trim_start_matches("sha256:")
            );
            h.bundle_sha256 = bundle_digest(&h)?;
            let manifest = panel_manifest(&h)?;
            specificity_write_panel_bundle_json(Path::new(&h.manifest_path), &manifest)?;
            specificity_write_panel_bundle_json(&root.join("handoff.json"), &h)?;
            Ok(h)
        })
    }

    fn validate_panel_multi(
        &self,
        h: &TranscriptAssayPanelMultiReferenceHandoff,
        root: &Path,
    ) -> Result<TranscriptAssayPanelReport, EngineError> {
        request_refs(&h.request, h.effective_limits)?;
        let panel = self.get_transcript_assay_panel_report(&h.request.panel_report_id)?;
        let order = assay_order(&panel)?;
        let expected = h.effective_limits.commands(
            order.len(),
            h.references.iter().filter(|r| r.child.is_some()).count(),
        )?;
        if h.schema != PANEL_MULTI_HANDOFF_SCHEMA
            || h.content_sha256 != scientific_id(h)?
            || h.bundle_sha256 != bundle_digest(h)?
            || h.handoff_id
                != format!(
                    "panel_multi_{}",
                    h.content_sha256.trim_start_matches("sha256:")
                )
            || h.nonclaims != nonclaims()
            || h.panel_report_id != panel.report_id
            || h.request.panel_report_id != panel.report_id
            || h.panel_digest != Self::transcript_assay_panel_specificity_digest(&panel)?
            || h.request
                .expected_panel_digest
                .as_ref()
                .is_some_and(|d| d != &h.panel_digest)
            || h.source_snapshot_sha256 != self.panel_multi_snapshot(&panel)?
            || h.expected_command_count != expected
            || h.references.len() != h.request.references.len()
            || Path::new(&h.manifest_path) != root.join("execution_manifest.json")
        {
            return Err(invalid("parent, panel or source binding changed"));
        }
        for (i, (r, requested)) in h.references.iter().zip(&h.request.references).enumerate() {
            if r.resolved_genome_id != requested.target_genome_id
                || r.expected_index_kind != requested.expected_index_kind
                || r.required != requested.required
            {
                return Err(invalid("reference declaration mismatch"));
            }
            let Some(child) = &r.child else {
                if r.required
                    || r.reference.is_some()
                    || r.availability != PrimerSpecificityMultiAvailability::Unavailable
                    || r.diagnostic.is_none()
                {
                    return Err(invalid("invalid unavailable reference"));
                }
                continue;
            };
            let dir = root.join(format!("reference_{i:02}"));
            let (_, policy) = Self::panel_specificity_policy(
                &panel,
                &r.resolved_genome_id,
                h.request.policy.clone(),
            )?;
            if child.schema != TRANSCRIPT_ASSAY_PANEL_SPECIFICITY_HANDOFF_SCHEMA
                || child.panel_report_id != panel.report_id
                || child.panel_digest != h.panel_digest
                || child.source_seq_id != panel.source_seq_id
                || child.source_feature_id != panel.source_feature_id
                || child.selected_assay_count != order.len()
                || child.assays.len() != order.len()
                || child.requested_target_genome_id != r.resolved_genome_id
                || child.resolved_target_genome_id != r.resolved_genome_id
                || child.policy_schema != PRIMER_SPECIFICITY_POLICY_SCHEMA
                || !multi_equal(&child.policy, &policy)
                || Path::new(&child.bundle_dir) != dir
                || child.completion_policy != "all_assays_all_commands_success_and_specificity_pass"
                || child.execution_manifest_schema
                    != TRANSCRIPT_ASSAY_PANEL_SPECIFICITY_EXECUTION_MANIFEST_SCHEMA
                || child.handoff_id
                    != Self::transcript_assay_panel_handoff_id(
                        &panel.report_id,
                        &h.panel_digest,
                        &child.requested_target_genome_id,
                        &child.resolved_target_genome_id,
                        &child.policy_schema,
                        &policy,
                        &child.completion_policy,
                        &child.execution_manifest_schema,
                        &child.assays,
                    )?
            {
                return Err(invalid("child handoff, policy or membership changed"));
            }
            check_child_file(
                &dir.join(format!("{}.json", child.handoff_id)),
                &child.handoff_path,
                child,
            )?;
            for (a, saved) in child.assays.iter().zip(&order) {
                let f = Self::primer_specificity_input_from_record(
                    PrimerSpecificityPrimerRole::Forward,
                    &saved.primer_pair.forward,
                )?;
                let rev = Self::primer_specificity_input_from_record(
                    PrimerSpecificityPrimerRole::Reverse,
                    &saved.primer_pair.reverse,
                )?;
                if a.assay_id != saved.assay_id
                    || a.assay_rank != saved.rank
                    || a.primer_pair_digest
                        != Self::transcript_assay_panel_primer_pair_digest(
                            &saved.assay_id,
                            saved.rank,
                            &f,
                            &rev,
                        )?
                    || a.forward_annealing_sequence != f.annealing_sequence
                    || a.reverse_annealing_sequence != rev.annealing_sequence
                    || a.handoff.primary_seq_id.as_ref() != Some(&panel.source_seq_id)
                    || a.handoff.pair_rank != Some(saved.rank)
                    || a.handoff.expected_amplicon_length_bp
                        != Some(saved.primer_pair.amplicon_length_bp)
                {
                    return Err(invalid("assay/oligo binding changed"));
                }
                let target =
                    self.primer_specificity_intended_target_from_transcript_assay(&panel, saved)?;
                let target = Self::primer_specificity_bind_intended_target(
                    &target,
                    a.handoff.blast_database.as_ref(),
                );
                let bundle = dir.join(format!(
                    "assay_{:03}_{}",
                    saved.rank,
                    short_sha256_id("assay", &saved.assay_id)
                ));
                self.specificity_validate_bound_child(
                    &a.handoff,
                    r,
                    &[f, rev],
                    &target,
                    &policy,
                    h.request.catalog_path.as_deref(),
                    h.request.cache_dir.as_deref(),
                    &bundle,
                )?;
                check_child_file(
                    &bundle.join(format!("{}.json", a.handoff.handoff_id)),
                    &a.handoff.handoff_path,
                    &a.handoff,
                )?;
            }
        }
        Ok(panel)
    }

    /// Check exact manifest coverage and all output budgets before interpreting a cell.
    /// Late refusal cannot leave a subset of the reports or any selection in the store.
    pub fn import_panel_multi_specificity(
        &mut self,
        handoff_path: &str,
        manifest_path: &str,
    ) -> Result<TranscriptAssayPanelMultiReferenceSummary, EngineError> {
        let h: TranscriptAssayPanelMultiReferenceHandoff =
            specificity_multi_read_json(handoff_path)?;
        let mut manifest: PanelSpecificityExecutionManifest =
            specificity_multi_read_json(manifest_path)?;
        let path = fs::canonicalize(handoff_path)
            .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?;
        let root = path
            .parent()
            .ok_or_else(|| invalid("missing bundle directory"))?;
        manifest
            .commands
            .sort_by(|a, b| a.command_id.cmp(&b.command_id));
        let (executions, declared_output_bytes) = manifest_rows(&h, &manifest)?;
        let panel = self.validate_panel_multi(&h, root)?;
        let mut budget = SpecificityOutputBudget {
            bytes_left: h.effective_limits.max_total_output_bytes,
            per_output: h.effective_limits.max_bytes_per_output,
            rows_left: h.effective_limits.max_parsed_rows,
        };
        let mut retained = BTreeMap::new();
        let mut retained_command_ids = BTreeSet::new();
        for r in &h.references {
            let Some(child) = &r.child else {
                continue;
            };
            for a in &child.assays {
                let mut complete = a.handoff.commands.iter().all(|c| {
                    executions.get(&c.command_id).is_some_and(|e| {
                        e.state == PrimerSpecificityMultiExecutionState::Completed
                            && e.exit_code == Some(0)
                    })
                });
                let mut diagnostics = Vec::new();
                let outputs = Self::specificity_multi_read_outputs(
                    &a.handoff,
                    &a.handoff.primers,
                    &executions,
                    &mut diagnostics,
                    &mut complete,
                    Some(&mut budget),
                )?;
                retained_command_ids.extend(outputs.keys().cloned());
                retained.insert(
                    (&r.resolved_genome_id, &a.assay_id),
                    (outputs, diagnostics, complete),
                );
            }
        }
        let mut assays = Vec::new();
        for (idx, saved) in assay_order(&panel)?.into_iter().enumerate() {
            let mut cells = Vec::new();
            let mut execution_complete = true;
            let mut pair_digest = None;
            for r in &h.references {
                let mut row = PrimerSpecificityMultiSummaryReference {
                    genome_id: r.resolved_genome_id.clone(),
                    required: r.required,
                    index_kind: r.expected_index_kind,
                    reference: r.reference.clone(),
                    verdict: PrimerSpecificityMultiVerdict::Incomplete,
                    applicability: "unavailable_at_preparation".into(),
                    diagnostics: Vec::new(),
                    report: None,
                };
                if let Some(child) = &r.child {
                    let a = &child.assays[idx];
                    pair_digest = Some(a.primer_pair_digest.clone());
                    let (outputs, diagnostics, complete) = retained
                        .remove(&(&r.resolved_genome_id, &a.assay_id))
                        .ok_or_else(|| invalid("missing retained cell"))?;
                    execution_complete &= complete;
                    row.diagnostics = diagnostics;
                    match self.validate_primer_specificity_handoff_database(&a.handoff) {
                        Ok(_) => {
                            row.applicability = "current_at_import".into();
                            self.specificity_multi_interpret_outputs(
                                &a.handoff, &outputs, &mut row,
                            )?;
                        }
                        Err(error) => {
                            row.applicability = "stale_or_unavailable_reference".into();
                            row.diagnostics.push(error.message);
                        }
                    }
                } else {
                    row.diagnostics.extend(r.diagnostic.iter().cloned());
                }
                cells.push(row);
            }
            assays.push(PanelSpecificityAssaySummary {
                assay_id: saved.assay_id.clone(),
                assay_rank: saved.rank,
                primer_pair_digest: pair_digest.ok_or_else(|| invalid("no prepared reference"))?,
                genomic: Self::specificity_multi_dimension(
                    &cells,
                    BlastDatabaseIndexKind::GenomicDna,
                    execution_complete,
                ),
                transcriptome: Self::specificity_multi_dimension(
                    &cells,
                    BlastDatabaseIndexKind::TranscriptomeCdna,
                    execution_complete,
                ),
                cells,
                execution_complete,
            });
        }
        let mut summary = TranscriptAssayPanelMultiReferenceSummary {
            schema: PANEL_MULTI_SUMMARY_SCHEMA.into(),
            summary_id: String::new(),
            content_sha256: String::new(),
            handoff_id: h.handoff_id,
            panel_report_id: panel.report_id,
            panel_digest: h.panel_digest,
            source_snapshot_sha256: h.source_snapshot_sha256,
            request: h.request,
            effective_limits: h.effective_limits,
            execution_manifest: manifest,
            execution_complete: assays.iter().all(|a| a.execution_complete),
            declared_output_bytes,
            retained_command_ids: retained_command_ids.into_iter().collect(),
            retained_output_bytes: h.effective_limits.max_total_output_bytes - budget.bytes_left,
            parsed_rows: h.effective_limits.max_parsed_rows - budget.rows_left,
            references: h.references.iter().map(ref_summary).collect(),
            assays,
            nonclaims: h.nonclaims,
        };
        summary.content_sha256 = Self::panel_multi_summary_digest(&summary)?;
        summary.summary_id = format!(
            "panel_multi_summary_{}",
            summary.content_sha256.trim_start_matches("sha256:")
        );
        let mut store = self.read_primer_design_store();
        if let Some(existing) = store
            .transcript_assay_panel_multi_reference_summaries
            .get(&summary.summary_id)
        {
            if existing.summary_id != summary.summary_id
                || Self::panel_multi_summary_digest(existing)? != summary.content_sha256
                || existing.content_sha256 != summary.content_sha256
            {
                return Err(invalid("immutable summary content collision"));
            }
            return Ok(existing.clone());
        }
        for report in summary
            .assays
            .iter()
            .flat_map(|a| &a.cells)
            .filter_map(|r| r.report.as_ref())
        {
            if let Some(existing) = store.primer_specificity_reports.get(&report.report_id)
                && Self::specificity_multi_report_basis(existing)
                    != Self::specificity_multi_report_basis(report)
            {
                return Err(invalid("immutable child report content collision"));
            }
            store
                .primer_specificity_reports
                .entry(report.report_id.clone())
                .or_insert_with(|| report.clone());
        }
        store
            .transcript_assay_panel_multi_reference_summaries
            .insert(summary.summary_id.clone(), summary.clone());
        self.write_primer_design_store_snapshot(store)?;
        Ok(summary)
    }

    fn panel_multi_summary_digest(
        s: &TranscriptAssayPanelMultiReferenceSummary,
    ) -> Result<String, EngineError> {
        let assays = s
            .assays
            .iter()
            .map(|a| {
                json!([
                    a.assay_id,
                    a.assay_rank,
                    a.primer_pair_digest,
                    a.genomic,
                    a.transcriptome,
                    a.execution_complete,
                    a.cells
                        .iter()
                        .map(|r| json!([
                            r.genome_id,
                            r.index_kind,
                            r.required,
                            r.reference,
                            r.applicability,
                            r.verdict,
                            r.report.as_ref().map(Self::specificity_multi_report_basis)
                        ]))
                        .collect::<Vec<_>>()
                ])
            })
            .collect::<Vec<_>>();
        let executions = s
            .execution_manifest
            .commands
            .iter()
            .map(|c| {
                (
                    &c.command_id,
                    &c.state,
                    c.exit_code,
                    c.output_size_bytes,
                    &c.output_sha256,
                )
            })
            .collect::<Vec<_>>();
        let references = s
            .references
            .iter()
            .map(|r| {
                (
                    &r.requested_genome_id,
                    &r.resolved_genome_id,
                    r.expected_index_kind,
                    r.required,
                    &r.reference,
                    &r.availability,
                )
            })
            .collect::<Vec<_>>();
        digest(&json!([
            s.schema,
            s.handoff_id,
            s.panel_report_id,
            s.panel_digest,
            s.source_snapshot_sha256,
            s.effective_limits,
            s.execution_complete,
            s.request.schema,
            s.request.panel_report_id,
            s.request.expected_panel_digest,
            s.request.policy,
            s.request.references,
            references,
            s.declared_output_bytes,
            s.retained_output_bytes,
            s.retained_command_ids,
            s.parsed_rows,
            assays,
            executions,
            s.nonclaims
        ]))
    }

    /// Read a retained summary without resolving its historical reference paths.
    pub fn get_panel_multi_specificity(
        &self,
        id: &str,
    ) -> Result<TranscriptAssayPanelMultiReferenceSummary, EngineError> {
        self.read_primer_design_store()
            .transcript_assay_panel_multi_reference_summaries
            .remove(id)
            .ok_or_else(|| {
                EngineError::new(
                    ErrorCode::NotFound,
                    format!("Unknown panel multi-reference summary: {id}"),
                )
            })
    }

    /// List immutable summaries in their content-ID order.
    pub fn list_panel_multi_specificity(&self) -> Vec<TranscriptAssayPanelMultiReferenceSummary> {
        self.read_primer_design_store()
            .transcript_assay_panel_multi_reference_summaries
            .into_values()
            .collect()
    }
}

#[cfg(test)]
#[path = "primer_specificity_panel_multi_tests.rs"]
mod tests;
