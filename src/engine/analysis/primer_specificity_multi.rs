//! Standalone, non-executing multi-reference primer specificity and receipt validation.
//! Biology stays in the existing single-reference interpreter; references never
//! become interchangeable, and this family never attaches panel readiness.

use super::operation_handlers::PrimerSpecificityResolvedInput;
use super::*;
use std::io::Read;
use std::path::Component;

const MULTI_NONCLAIMS: &[&str] = &[
    "Specificity applies only to the explicitly bound prepared references and policy, not every genome or transcript.",
    "A hash binds retained content; it does not authenticate external process execution.",
    "This standalone summary does not establish isoform discrimination, oligo QC, laboratory validation or order readiness.",
    "Current applicability is checked at import; historical display performs no database probes.",
];

const MULTI_MAX_OUTPUT_BYTES: u64 = 128 * 1024 * 1024;

fn multi_equal<T: Serialize>(left: &T, right: &T) -> bool {
    match (serde_json::to_value(left), serde_json::to_value(right)) {
        (Ok(left), Ok(right)) => left == right,
        _ => false,
    }
}

impl GentleEngine {
    fn specificity_multi_source_snapshot(
        &self,
        input: &PrimerSpecificityResolvedInput,
    ) -> Result<Option<String>, EngineError> {
        input.primary_seq_id.as_ref().map(|id| {
            let dna = self.state.sequences.get(id).ok_or_else(|| EngineError::invalid_input("Saved primer-pair template is no longer loaded"))?;
            let anchor = self.sequence_genome_anchor_summary(id).ok();
            serde_json::to_vec(&json!({"record":dna.clone_seq_record(), "overhang":dna.overhang(), "anchor":anchor, "source_reference":input.intended_target.source_reference}))
                .map(|bytes| sha256_prefixed_bytes(&bytes)).map_err(|e| EngineError::internal(e.to_string()))
        }).transpose()
    }
    fn validate_specificity_multi_request(
        request: &PrimerSpecificityMultiRequest,
    ) -> Result<(), EngineError> {
        if request.schema != PRIMER_SPECIFICITY_MULTI_REQUEST_SCHEMA
            || !(1..=8).contains(&request.references.len())
            || !request
                .references
                .iter()
                .any(|reference| reference.required)
            || request.policy.specificity_target_genome_id.is_some()
        {
            return Err(EngineError::invalid_input(
                "specificity_multi_request_invalid: require v1, 1..=8 explicit references, at least one required reference and no common policy target selector",
            ));
        }
        let mut seen = HashSet::new();
        for reference in &request.references {
            let id = reference.target_genome_id.trim();
            if id.is_empty()
                || id.contains(['*', '?'])
                || id.chars().any(char::is_control)
                || !seen.insert(id.to_ascii_lowercase())
            {
                return Err(EngineError::invalid_input(
                    "specificity_multi_reference_invalid: empty, wildcard or duplicate reference",
                ));
            }
            if let Some(mapping) = &reference.intended_target {
                if mapping.model == PrimerSpecificityIntendedTargetModel::Unknown
                    || (mapping.expected_products.is_empty()
                        && mapping
                            .subject_id
                            .as_deref()
                            .is_none_or(|s| s.trim().is_empty()))
                    || mapping
                        .forward_binding_ranges
                        .iter()
                        .chain(&mapping.reverse_binding_ranges)
                        .any(|r| r.start_1based == 0 || r.end_1based < r.start_1based)
                    || mapping
                        .expected_product_range
                        .as_ref()
                        .is_some_and(|r| r.start_1based == 0 || r.end_1based < r.start_1based)
                    || mapping.expected_products.iter().any(|p| {
                        p.subject_id.trim().is_empty()
                            || !["genomic_dna", "transcriptome_cdna"]
                                .contains(&p.target_space.as_str())
                            || p.expected_product_range.as_ref().is_some_and(|r| {
                                r.start_1based == 0 || r.end_1based < r.start_1based
                            })
                    })
                {
                    return Err(EngineError::invalid_input(
                        "specificity_multi_mapping_invalid: explicit mappings need known geometry and valid subject/ranges",
                    ));
                }
            }
        }
        Ok(())
    }

    fn specificity_multi_input(
        &self,
        pair: &PrimerSpecificityMultiPair,
    ) -> Result<PrimerSpecificityResolvedInput, EngineError> {
        match pair {
            PrimerSpecificityMultiPair::SavedPair {
                primer_report_id,
                pair_rank,
                pair_index,
            } => {
                if pair_rank.is_some() && pair_index.is_some() {
                    return Err(EngineError::invalid_input(
                        "Select a saved pair by rank or index, not both",
                    ));
                }
                self.resolve_primer_specificity_input(
                    Some(primer_report_id),
                    *pair_rank,
                    *pair_index,
                    None,
                    None,
                )
            }
            PrimerSpecificityMultiPair::ExplicitPair { forward, reverse } => {
                let canonical = |primer: &PrimerSpecificityInputPrimer, role| {
                    let record = PrimerDesignPrimerRecord {
                        sequence: primer.full_sequence.clone(),
                        non_annealing_5prime_tail_bp: primer.non_annealing_5prime_tail_bp,
                        ..Default::default()
                    };
                    let value = Self::primer_specificity_input_from_record(role, &record)?;
                    if serde_json::to_value(&value).ok() != serde_json::to_value(primer).ok() {
                        return Err(EngineError::invalid_input(
                            "Explicit primers must give consistent normalized full/annealing sequences, roles, lengths and 5-prime tail boundaries",
                        ));
                    }
                    Ok(value)
                };
                let forward = canonical(forward, PrimerSpecificityPrimerRole::Forward)?;
                let reverse = canonical(reverse, PrimerSpecificityPrimerRole::Reverse)?;
                let mut input = self.resolve_primer_specificity_input(
                    None,
                    None,
                    None,
                    Some(&forward.full_sequence),
                    Some(&reverse.full_sequence),
                )?;
                input.forward = forward;
                input.reverse = reverse;
                Ok(input)
            }
        }
    }

    fn specificity_multi_handoff_digest(
        handoff: &PrimerSpecificityMultiHandoff,
    ) -> Result<String, EngineError> {
        let mut basis = handoff.clone();
        basis.handoff_id.clear();
        basis.content_sha256.clear();
        serde_json::to_vec(&basis)
            .map(|b| sha256_prefixed_bytes(&b))
            .map_err(|e| {
                EngineError::new(
                    ErrorCode::Internal,
                    format!("Could not bind multi-reference handoff: {e}"),
                )
            })
    }

    fn specificity_multi_command_digest(
        command: &PrimerSpecificityHandoffCommand,
    ) -> Result<String, EngineError> {
        serde_json::to_vec(command)
            .map(|b| sha256_prefixed_bytes(&b))
            .map_err(|e| {
                EngineError::new(
                    ErrorCode::Internal,
                    format!("Could not bind BLAST command: {e}"),
                )
            })
    }

    /// Inspect every exact reference before writing any file; never prepares resources.
    pub fn prepare_primer_pair_multi_reference_specificity_handoff(
        &self,
        mut request: PrimerSpecificityMultiRequest,
        output_dir: &str,
    ) -> Result<PrimerSpecificityMultiHandoff, EngineError> {
        Self::validate_specificity_multi_request(&request)?;
        let input = self.specificity_multi_input(&request.pair)?;
        let source_snapshot_sha256 = self.specificity_multi_source_snapshot(&input)?;
        let (catalog, catalog_path) =
            Self::open_reference_genome_catalog(request.catalog_path.as_deref())?;
        let cache_dir = request.cache_dir.clone();
        let cache = cache_dir.as_deref();
        let mut resolved = HashSet::new();
        let mut references = Vec::new();
        for reference in &mut request.references {
            let requested = reference.target_genome_id.trim().to_string();
            let (id, entry) = catalog
                .exact_catalog_entry(&requested)
                .map_err(EngineError::invalid_input)?;
            if entry
                .blast_index_kind
                .is_some_and(|kind| kind != reference.expected_index_kind)
            {
                return Err(EngineError::invalid_input(format!(
                    "specificity_multi_index_kind_mismatch: {id}"
                )));
            }
            if !resolved.insert(id.clone()) {
                return Err(EngineError::invalid_input(
                    "specificity_multi_duplicate_resolved_reference",
                ));
            }
            reference.target_genome_id = id.clone();
            let inspection = catalog.inspect_prepared_genome(&id, cache);
            let (identity, diagnostic) = match inspection {
                Ok(Some(inspection)) if inspection.blast_index_ready => {
                    match inspection.blast_database.as_ref() {
                        Some(database) if database.index_kind != reference.expected_index_kind => {
                            return Err(EngineError::invalid_input(format!(
                                "specificity_multi_index_kind_mismatch: {id}"
                            )));
                        }
                        Some(database)
                            if database
                                .sequence_count
                                .is_some_and(|n| n > 0 && n.checked_add(1).is_some()) =>
                        {
                            (
                                Self::primer_specificity_reference_identity(database),
                                Some(
                                    "BLAST reference validation or content fingerprint unavailable"
                                        .to_string(),
                                ),
                            )
                        }
                        _ => (
                            None,
                            Some("Validated positive BLAST sequence count unavailable".to_string()),
                        ),
                    }
                }
                Ok(_) => (
                    None,
                    Some("Prepared BLAST index unavailable; no resource was prepared".to_string()),
                ),
                Err(error) => (
                    None,
                    Some(format!("Prepared reference inspection failed: {error}")),
                ),
            };
            if identity.is_none() && reference.required {
                return Err(EngineError::invalid_input(format!(
                    "specificity_multi_required_reference_unavailable: {id}: {}",
                    diagnostic.as_deref().unwrap_or("unavailable")
                )));
            }
            if reference
                .intended_target
                .as_ref()
                .and_then(|mapping| mapping.reference_binding.as_ref())
                .is_some_and(|binding| identity.as_ref() != Some(binding))
            {
                return Err(EngineError::invalid_input(
                    "Caller-provided mapping has a stale or incompatible reference binding",
                ));
            }
            references.push(PrimerSpecificityMultiHandoffReference {
                requested_genome_id: requested,
                resolved_genome_id: id,
                required: reference.required,
                expected_index_kind: reference.expected_index_kind,
                availability: if identity.is_some() {
                    PrimerSpecificityMultiAvailability::Prepared
                } else {
                    PrimerSpecificityMultiAvailability::Unavailable
                },
                diagnostic: if identity.is_none() { diagnostic } else { None },
                reference: identity,
                child: None,
            });
        }
        references.sort_by(|a, b| a.resolved_genome_id.cmp(&b.resolved_genome_id));
        request
            .references
            .sort_by(|a, b| a.target_genome_id.cmp(&b.target_genome_id));
        request.catalog_path = Some(catalog_path);
        request.policy = Self::normalize_primer_specificity_policy(
            &references[0].resolved_genome_id,
            request.policy,
        )?
        .1;
        request.policy.specificity_target_genome_id = None;
        let destination = specificity_multi_destination(output_dir)?;
        fs::create_dir(&destination).map_err(|e| {
            EngineError::new(
                ErrorCode::Io,
                format!("Could not create multi-reference bundle: {e}"),
            )
        })?;
        let build = (|| {
            for (ordinal, row) in references.iter_mut().enumerate() {
                let Some(identity) = &row.reference else {
                    continue;
                };
                let mut child_input = input.clone();
                let reference = &request.references[ordinal];
                if let Some(mapping) = &reference.intended_target {
                    if mapping
                        .reference_binding
                        .as_ref()
                        .is_some_and(|binding| binding != identity)
                    {
                        return Err(EngineError::invalid_input(
                            "Caller-provided target binding disagrees with the inspected reference",
                        ));
                    }
                    child_input.intended_target = mapping.clone();
                    child_input.intended_target.reference_binding = Some(identity.clone());
                    child_input.intended_target.source =
                        "caller_provided_per_reference_mapping".to_string();
                    child_input.intended_target.warnings.push("Caller-provided geometry, not independently inferred orthology or liftover.".to_string());
                }
                let child_dir = destination.join(format!("reference_{ordinal:02}"));
                let child = self.prepare_primer_specificity_handoff_resolved(
                    child_input,
                    &row.resolved_genome_id,
                    request.policy.clone(),
                    request.catalog_path.as_deref(),
                    cache,
                    &child_dir.to_string_lossy(),
                )?;
                if child.resolved_target_genome_id != row.resolved_genome_id
                    || child
                        .blast_database
                        .as_ref()
                        .and_then(Self::primer_specificity_reference_identity)
                        .as_ref()
                        != Some(identity)
                {
                    return Err(EngineError::invalid_input(
                        "Reference changed or fallback resolution occurred while preparing the bundle",
                    ));
                }
                row.child = Some(child);
            }
            let mut handoff = PrimerSpecificityMultiHandoff {
                schema: PRIMER_SPECIFICITY_MULTI_HANDOFF_SCHEMA.into(),
                handoff_id: String::new(),
                content_sha256: String::new(),
                pair_binding_sha256: Self::primer_specificity_pair_binding(
                    &input.forward,
                    &input.reverse,
                )?,
                source_snapshot_sha256,
                primers: vec![input.forward, input.reverse],
                request,
                references,
                manifest_path: destination
                    .join("execution_manifest.json")
                    .to_string_lossy()
                    .to_string(),
                nonclaims: MULTI_NONCLAIMS.iter().map(|s| s.to_string()).collect(),
            };
            handoff.content_sha256 = Self::specificity_multi_handoff_digest(&handoff)?;
            handoff.handoff_id = format!(
                "primer_multi_{}",
                handoff.content_sha256.trim_start_matches("sha256:")
            );
            let manifest = PrimerSpecificityMultiExecutionManifest {
                schema: PRIMER_SPECIFICITY_MULTI_MANIFEST_SCHEMA.into(),
                handoff_id: handoff.handoff_id.clone(),
                handoff_content_sha256: handoff.content_sha256.clone(),
                pair_binding_sha256: handoff.pair_binding_sha256.clone(),
                commands: handoff
                    .references
                    .iter()
                    .filter_map(|r| r.child.as_ref())
                    .flat_map(|c| &c.commands)
                    .map(|command| {
                        Ok(PrimerSpecificityMultiExecutionCommand {
                            command_id: command.command_id.clone(),
                            command_sha256: Self::specificity_multi_command_digest(command)?,
                            output_path: command.output_tsv_path.clone(),
                            state: PrimerSpecificityMultiExecutionState::Pending,
                            exit_code: None,
                            output_size_bytes: None,
                            output_sha256: None,
                        })
                    })
                    .collect::<Result<Vec<_>, EngineError>>()?,
            };
            specificity_multi_write_json(&destination.join("handoff.json"), &handoff)?;
            specificity_multi_write_json(Path::new(&handoff.manifest_path), &manifest)?;
            Ok(handoff)
        })();
        if build.is_err() {
            let _ = fs::remove_dir_all(&destination);
        }
        build
    }
}

impl GentleEngine {
    fn specificity_multi_validate_handoff(
        &self,
        handoff: &PrimerSpecificityMultiHandoff,
    ) -> Result<(), EngineError> {
        Self::validate_specificity_multi_request(&handoff.request)?;
        let digest = Self::specificity_multi_handoff_digest(handoff)?;
        if handoff.schema != PRIMER_SPECIFICITY_MULTI_HANDOFF_SCHEMA
            || handoff.content_sha256 != digest
            || handoff.handoff_id
                != format!("primer_multi_{}", digest.trim_start_matches("sha256:"))
            || handoff.nonclaims
                != MULTI_NONCLAIMS
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>()
            || handoff.references.len() != handoff.request.references.len()
        {
            return Err(EngineError::invalid_input(
                "specificity_multi_handoff_binding_invalid",
            ));
        }
        let input = self.specificity_multi_input(&handoff.request.pair)?;
        if handoff.source_snapshot_sha256 != self.specificity_multi_source_snapshot(&input)?
            || handoff.pair_binding_sha256
                != Self::primer_specificity_pair_binding(&input.forward, &input.reverse)?
            || !multi_equal(
                &handoff.primers,
                &vec![input.forward.clone(), input.reverse.clone()],
            )
        {
            return Err(EngineError::invalid_input(
                "specificity_multi_pair_or_source_changed",
            ));
        }
        let mut previous: Option<&str> = None;
        for (ordinal, (row, requested)) in handoff
            .references
            .iter()
            .zip(&handoff.request.references)
            .enumerate()
        {
            if row.resolved_genome_id != requested.target_genome_id
                || row.required != requested.required
                || row.expected_index_kind != requested.expected_index_kind
                || previous.is_some_and(|id| id >= row.resolved_genome_id.as_str())
            {
                return Err(EngineError::invalid_input(
                    "specificity_multi_effective_reference_list_invalid",
                ));
            }
            previous = Some(&row.resolved_genome_id);
            let Some(child) = &row.child else {
                if row.required
                    || row.reference.is_some()
                    || row.availability != PrimerSpecificityMultiAvailability::Unavailable
                {
                    return Err(EngineError::invalid_input(
                        "specificity_multi_unavailable_reference_invalid",
                    ));
                }
                continue;
            };
            let reference = row
                .reference
                .as_ref()
                .ok_or_else(|| EngineError::invalid_input("Missing reference content identity"))?;
            let database = child.blast_database.as_ref().ok_or_else(|| {
                EngineError::invalid_input("Legacy child database identity unavailable")
            })?;
            let expected_policy = Self::normalize_primer_specificity_policy(
                &row.resolved_genome_id,
                handoff.request.policy.clone(),
            )?
            .1;
            let expected_options = self.resolve_blast_options_for_request(
                None,
                Some("blastn-short"),
                Some(expected_policy.max_hits_per_primer),
            )?;
            let mut expected_target = requested
                .intended_target
                .clone()
                .unwrap_or_else(|| input.intended_target.clone());
            if requested.intended_target.is_some() {
                if expected_target
                    .reference_binding
                    .as_ref()
                    .is_some_and(|b| b != reference)
                {
                    return Err(EngineError::invalid_input("Caller mapping is stale"));
                }
                expected_target.reference_binding = Some(reference.clone());
                expected_target.source = "caller_provided_per_reference_mapping".into();
                expected_target.warnings.push(
                    "Caller-provided geometry, not independently inferred orthology or liftover."
                        .into(),
                );
            }
            expected_target =
                Self::primer_specificity_bind_intended_target(&expected_target, Some(database));
            let bundle = Path::new(&handoff.manifest_path)
                .parent()
                .ok_or_else(|| EngineError::invalid_input("Invalid manifest path"))?
                .join(format!("reference_{ordinal:02}"));
            if child.schema != PRIMER_SPECIFICITY_HANDOFF_SCHEMA
                || child.handoff_id != Self::primer_specificity_handoff_id_from_record(child)?
                || row.availability != PrimerSpecificityMultiAvailability::Prepared
                || reference.index_kind != row.expected_index_kind
                || reference.genome_id != row.resolved_genome_id
                || Some(reference.clone()) != Self::primer_specificity_reference_identity(database)
                || child.resolved_target_genome_id != row.resolved_genome_id
                || child.requested_target_genome_id != row.resolved_genome_id
                || child.catalog_path != handoff.request.catalog_path
                || child.cache_dir != handoff.request.cache_dir
                || !multi_equal(&child.policy, &expected_policy)
                || !multi_equal(
                    &child.effective_blast_options,
                    &Some(expected_options.clone()),
                )
                || !multi_equal(&child.primers, &handoff.primers)
                || !multi_equal(&child.intended_target, &expected_target)
                || child.commands.len() != 2
                || Path::new(&child.bundle_dir) != bundle
                || child.blast_db_prefix != database.prefix
                || child.completion_policy != "all_commands_success"
            {
                return Err(EngineError::invalid_input(
                    "specificity_multi_child_binding_invalid",
                ));
            }
            let subject_limit = database
                .sequence_count
                .filter(|n| *n > 0)
                .and_then(|n| n.checked_add(1))
                .ok_or_else(|| {
                    EngineError::invalid_input("Missing complete database sequence count")
                })?;
            for (command, primer) in child.commands.iter().zip(&handoff.primers) {
                let role = primer.role.as_str();
                let query = bundle
                    .join(format!("{}.{role}.fa", child.handoff_id))
                    .to_string_lossy()
                    .to_string();
                let output = bundle
                    .join(format!("{}.{role}.blast.tsv", child.handoff_id))
                    .to_string_lossy()
                    .to_string();
                let expected_args = vec![
                    "-db".into(),
                    child.blast_db_prefix.clone(),
                    "-query".into(),
                    query.clone(),
                    "-task".into(),
                    expected_options.task.clone(),
                    "-outfmt".into(),
                    super::operation_handlers::PRIMER_SPECIFICITY_BLASTN_OUTFMT_FIELDS.into(),
                    "-evalue".into(),
                    "1000".into(),
                    "-dust".into(),
                    "no".into(),
                    "-soft_masking".into(),
                    "false".into(),
                    "-max_target_seqs".into(),
                    subject_limit.to_string(),
                    "-out".into(),
                    output.clone(),
                ];
                if command.role != primer.role
                    || command.command_id != format!("{}:{role}", child.handoff_id)
                    || command.query_label != format!("{role}_annealing_segment")
                    || command.query_length_bp != primer.annealing_length_bp
                    || command.query_fasta_path != query
                    || command.output_tsv_path != output
                    || command.args != expected_args
                    || command.success_exit_codes != vec![0]
                    || command.program.trim().is_empty()
                    || (!child.blast_preflight.blastn.executable.trim().is_empty()
                        && command.program != child.blast_preflight.blastn.executable)
                {
                    return Err(EngineError::invalid_input(
                        "specificity_multi_expected_command_changed",
                    ));
                }
            }
            let actual = Self::primer_specificity_search_completeness_for_commands(
                Some(database),
                &child
                    .commands
                    .iter()
                    .map(|c| c.args.clone())
                    .collect::<Vec<_>>(),
            );
            if !actual.complete || !multi_equal(&actual, &child.search_completeness) {
                return Err(EngineError::invalid_input(
                    "specificity_multi_search_completeness_invalid",
                ));
            }
        }
        Ok(())
    }

    fn specificity_multi_validate_manifest<'a>(
        handoff: &PrimerSpecificityMultiHandoff,
        manifest: &'a PrimerSpecificityMultiExecutionManifest,
    ) -> Result<BTreeMap<String, &'a PrimerSpecificityMultiExecutionCommand>, EngineError> {
        if manifest.schema != PRIMER_SPECIFICITY_MULTI_MANIFEST_SCHEMA
            || manifest.handoff_id != handoff.handoff_id
            || manifest.handoff_content_sha256 != handoff.content_sha256
            || manifest.pair_binding_sha256 != handoff.pair_binding_sha256
        {
            return Err(EngineError::invalid_input(
                "specificity_multi_manifest_binding_invalid",
            ));
        }
        let expected = handoff
            .references
            .iter()
            .filter_map(|r| r.child.as_ref())
            .flat_map(|child| &child.commands)
            .map(|c| (c.command_id.as_str(), c))
            .collect::<BTreeMap<_, _>>();
        let mut executions = BTreeMap::new();
        for execution in &manifest.commands {
            let command = expected
                .get(execution.command_id.as_str())
                .ok_or_else(|| EngineError::invalid_input("Unexpected manifest command"))?;
            if execution.command_sha256 != Self::specificity_multi_command_digest(command)?
                || execution.output_path != command.output_tsv_path
                || executions
                    .insert(execution.command_id.clone(), execution)
                    .is_some()
            {
                return Err(EngineError::invalid_input(
                    "specificity_multi_manifest_command_changed_or_duplicated",
                ));
            }
        }
        Ok(executions)
    }

    fn specificity_multi_report_basis(report: &PrimerSpecificityReport) -> Value {
        json!({"schema": report.schema, "pair": report.primers, "reference":report.blast_database.as_ref().and_then(Self::primer_specificity_reference_identity),
            "intended_target":report.intended_target, "policy":report.policy, "forward_hits":report.forward_hits, "reverse_hits":report.reverse_hits, "amplicons":report.amplicons,
            "completeness":report.search_completeness, "compaction":report.compaction, "summary":report.summary,
            "genomic":report.genomic_specificity, "transcriptome":report.transcriptome_specificity,
            "blast_options":report.blast_runs.iter().map(|r| &r.effective_options_json).collect::<Vec<_>>(),
            "raw_outputs":report.raw_detail_artifacts.iter().map(|a| json!({"kind":a.source_kind, "source_id":a.source_id, "checksum":a.checksum, "checksum_algorithm":a.checksum_algorithm})).collect::<Vec<_>>()})
    }

    fn specificity_multi_summary_digest(
        summary: &PrimerSpecificityMultiSummary,
    ) -> Result<String, EngineError> {
        let mut request = serde_json::to_value(&summary.request)
            .map_err(|e| EngineError::internal(e.to_string()))?;
        request.as_object_mut().unwrap().remove("catalog_path");
        request.as_object_mut().unwrap().remove("cache_dir");
        let basis = json!({"schema":summary.schema, "request":request, "pair_binding":summary.pair_binding_sha256,
            "source_snapshot":summary.source_snapshot_sha256,
            "execution_complete":summary.execution_complete, "genomic":summary.genomic, "transcriptome":summary.transcriptome, "nonclaims":summary.nonclaims,
            "executions": summary.execution_manifest.commands.iter().map(|c| json!({"command_id":c.command_id, "state":c.state, "exit_code":c.exit_code, "output_size_bytes":c.output_size_bytes, "output_sha256":c.output_sha256})).collect::<Vec<_>>(),
            "references":summary.references.iter().map(|r| json!({"genome_id":r.genome_id, "required":r.required,"kind":r.index_kind,"reference":r.reference,"verdict":r.verdict,"applicability":r.applicability,
                "report":r.report.as_ref().map(Self::specificity_multi_report_basis)})).collect::<Vec<_>>()});
        serde_json::to_vec(&basis)
            .map(|bytes| sha256_prefixed_bytes(&bytes))
            .map_err(|e| EngineError::internal(e.to_string()))
    }

    fn specificity_multi_dimension(
        rows: &[PrimerSpecificityMultiSummaryReference],
        kind: BlastDatabaseIndexKind,
        execution_complete: bool,
    ) -> PrimerSpecificityMultiVerdict {
        let rows = rows
            .iter()
            .filter(|r| r.index_kind == kind)
            .collect::<Vec<_>>();
        if rows.is_empty() {
            return PrimerSpecificityMultiVerdict::NotRequested;
        }
        let required = rows.iter().filter(|r| r.required).collect::<Vec<_>>();
        if required.is_empty() {
            return PrimerSpecificityMultiVerdict::NotRequired;
        }
        if required
            .iter()
            .any(|r| r.verdict == PrimerSpecificityMultiVerdict::Fail)
        {
            return PrimerSpecificityMultiVerdict::Fail;
        }
        if execution_complete
            && required
                .iter()
                .all(|r| r.verdict == PrimerSpecificityMultiVerdict::Pass)
        {
            PrimerSpecificityMultiVerdict::Pass
        } else {
            PrimerSpecificityMultiVerdict::Incomplete
        }
    }

    /// Validate retained outputs without persisting anything until every row is classified.
    pub fn import_primer_pair_multi_reference_specificity(
        &mut self,
        handoff_path: &str,
        manifest_path: &str,
    ) -> Result<PrimerSpecificityMultiSummary, EngineError> {
        let handoff: PrimerSpecificityMultiHandoff = specificity_multi_read_json(handoff_path)?;
        let mut manifest: PrimerSpecificityMultiExecutionManifest =
            specificity_multi_read_json(manifest_path)?;
        self.specificity_multi_validate_handoff(&handoff)?;
        manifest
            .commands
            .sort_by(|a, b| a.command_id.cmp(&b.command_id));
        let executions = Self::specificity_multi_validate_manifest(&handoff, &manifest)?;
        let mut execution_complete = handoff
            .references
            .iter()
            .filter_map(|r| r.child.as_ref())
            .flat_map(|child| &child.commands)
            .all(|command| {
                executions.get(&command.command_id).is_some_and(|e| {
                    e.state == PrimerSpecificityMultiExecutionState::Completed
                        && e.exit_code == Some(0)
                })
            });
        let mut rows = Vec::new();
        for reference in &handoff.references {
            let mut row = PrimerSpecificityMultiSummaryReference {
                genome_id: reference.resolved_genome_id.clone(),
                required: reference.required,
                index_kind: reference.expected_index_kind,
                reference: reference.reference.clone(),
                verdict: PrimerSpecificityMultiVerdict::Incomplete,
                applicability: "unavailable_at_preparation".into(),
                diagnostics: Vec::new(),
                report: None,
            };
            let Some(child) = &reference.child else {
                row.diagnostics
                    .push(reference.diagnostic.clone().unwrap_or_else(|| {
                        "Optional reference unavailable; no search was requested".into()
                    }));
                rows.push(row);
                continue;
            };
            if let Err(error) = self.validate_primer_specificity_handoff_database(child) {
                row.applicability = "stale_or_unavailable_reference".into();
                row.diagnostics.push(error.message);
                rows.push(row);
                continue;
            }
            row.applicability = "current_at_import".into();
            let mut outputs = BTreeMap::new();
            for (command, primer) in child.commands.iter().zip(&handoff.primers) {
                let Some(execution) = executions.get(&command.command_id) else {
                    row.diagnostics
                        .push(format!("Missing execution: {}", command.command_id));
                    continue;
                };
                if execution.state != PrimerSpecificityMultiExecutionState::Completed
                    || execution.exit_code != Some(0)
                {
                    row.diagnostics.push(format!(
                        "Execution is {:?}, exit {:?}: {}",
                        execution.state, execution.exit_code, command.command_id
                    ));
                    continue;
                }
                let bytes = match specificity_multi_read_bytes(
                    &command.output_tsv_path,
                    MULTI_MAX_OUTPUT_BYTES,
                ) {
                    Ok(bytes) => bytes,
                    Err(error) if error.code == ErrorCode::Io => {
                        execution_complete = false;
                        row.diagnostics.push(format!(
                            "Retained output unavailable for {}: {}",
                            command.command_id, error.message
                        ));
                        continue;
                    }
                    Err(error) => return Err(error),
                };
                if execution.output_size_bytes != Some(bytes.len() as u64)
                    || execution.output_sha256.as_deref()
                        != Some(sha256_prefixed_bytes(&bytes).as_str())
                {
                    return Err(EngineError::invalid_input(
                        "specificity_multi_output_size_or_hash_mismatch",
                    ));
                }
                let query = specificity_multi_read_bytes(&command.query_fasta_path, 1024 * 1024)?;
                if query
                    != format!(">{}\n{}\n", command.query_label, primer.annealing_sequence)
                        .as_bytes()
                {
                    return Err(EngineError::invalid_input(
                        "specificity_multi_query_sequence_changed",
                    ));
                }
                let text = String::from_utf8(bytes)
                    .map_err(|_| EngineError::invalid_input("BLAST output must be UTF-8"))?;
                if !crate::genomes::parse_blastn_tabular_hits(&text)
                    .1
                    .is_empty()
                {
                    return Err(EngineError::invalid_input(
                        "Malformed BLAST rows cannot establish complete specificity",
                    ));
                }
                outputs.insert(command.command_id.clone(), text);
            }
            if outputs.len() == child.commands.len() {
                match self.primer_specificity_report_from_handoff_outputs(child, Some(&outputs)) {
                    Ok(mut report) => {
                        // Recheck after subject-window/annotation lookups, not just before parsing.
                        if let Err(error) = self.validate_primer_specificity_handoff_database(child)
                        {
                            row.applicability = "stale_or_unavailable_reference".into();
                            row.diagnostics.push(error.message);
                        } else {
                            row.verdict = if !report.search_completeness.complete
                                || report.intended_target.reference_binding.is_none()
                            {
                                PrimerSpecificityMultiVerdict::Incomplete
                            } else if report.summary.status == "pass"
                                && report.summary.specificity_pass
                            {
                                PrimerSpecificityMultiVerdict::Pass
                            } else if report.summary.status == "fail" {
                                PrimerSpecificityMultiVerdict::Fail
                            } else {
                                PrimerSpecificityMultiVerdict::Incomplete
                            };
                            let digest =
                                serde_json::to_vec(&Self::specificity_multi_report_basis(&report))
                                    .map(|b| sha256_prefixed_bytes(&b))
                                    .map_err(|e| EngineError::internal(e.to_string()))?;
                            report.report_id = format!(
                                "primer_multi_child_{}",
                                digest.trim_start_matches("sha256:")
                            );
                            row.diagnostics.push(report.summary.summary.clone());
                            row.report = Some(report);
                        }
                    }
                    Err(error) => row.diagnostics.push(format!(
                        "Shared interpretation incomplete: {}",
                        error.message
                    )),
                }
            }
            rows.push(row);
        }
        let mut summary = PrimerSpecificityMultiSummary {
            schema: PRIMER_SPECIFICITY_MULTI_SUMMARY_SCHEMA.into(),
            summary_id: String::new(),
            content_sha256: String::new(),
            handoff_id: handoff.handoff_id,
            pair_binding_sha256: handoff.pair_binding_sha256,
            source_snapshot_sha256: handoff.source_snapshot_sha256,
            request: handoff.request,
            genomic: Self::specificity_multi_dimension(
                &rows,
                BlastDatabaseIndexKind::GenomicDna,
                execution_complete,
            ),
            transcriptome: Self::specificity_multi_dimension(
                &rows,
                BlastDatabaseIndexKind::TranscriptomeCdna,
                execution_complete,
            ),
            references: rows,
            execution_manifest: manifest,
            execution_complete,
            nonclaims: handoff.nonclaims,
        };
        summary.content_sha256 = Self::specificity_multi_summary_digest(&summary)?;
        summary.summary_id = format!(
            "primer_multi_summary_{}",
            summary.content_sha256.trim_start_matches("sha256:")
        );
        let mut store = self.read_primer_design_store();
        if let Some(existing) = store
            .primer_specificity_multi_summaries
            .get(&summary.summary_id)
        {
            if Self::specificity_multi_summary_digest(existing)? != summary.content_sha256 {
                return Err(EngineError::invalid_input(
                    "Immutable summary content collision",
                ));
            }
            return Ok(existing.clone());
        }
        for report in summary.references.iter().filter_map(|r| r.report.as_ref()) {
            store
                .primer_specificity_reports
                .entry(report.report_id.clone())
                .or_insert_with(|| report.clone());
        }
        store
            .primer_specificity_multi_summaries
            .insert(summary.summary_id.clone(), summary.clone());
        self.write_primer_design_store(store)?;
        Ok(summary)
    }

    pub fn get_primer_pair_multi_reference_specificity_summary(
        &self,
        id: &str,
    ) -> Result<PrimerSpecificityMultiSummary, EngineError> {
        self.read_primer_design_store()
            .primer_specificity_multi_summaries
            .get(id)
            .cloned()
            .ok_or_else(|| {
                EngineError::new(
                    ErrorCode::NotFound,
                    format!("Multi-reference summary '{id}' not found"),
                )
            })
    }

    pub fn list_primer_pair_multi_reference_specificity_summaries(
        &self,
    ) -> Vec<PrimerSpecificityMultiSummary> {
        self.read_primer_design_store()
            .primer_specificity_multi_summaries
            .into_values()
            .collect()
    }
}

fn specificity_multi_read_bytes(path: &str, limit: u64) -> Result<Vec<u8>, EngineError> {
    let path = Path::new(path);
    if path
        .components()
        .any(|c| c == Component::ParentDir || c.as_os_str() == "..")
    {
        return Err(EngineError::invalid_input(
            "Parent traversal is not permitted in receipt paths",
        ));
    }
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?
            .join(path)
    };
    specificity_multi_ordinary_parent(&absolute)?;
    let metadata = fs::symlink_metadata(path).map_err(|e| {
        EngineError::new(
            ErrorCode::Io,
            format!("Could not inspect '{}': {e}", path.display()),
        )
    })?;
    if !metadata.is_file() || metadata.file_type().is_symlink() || metadata.len() > limit {
        return Err(EngineError::invalid_input(
            "Receipt/output must be an ordinary file within the disclosed byte budget",
        ));
    }
    let mut bytes = Vec::new();
    File::open(path)
        .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?
        .take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?;
    if bytes.len() as u64 > limit {
        return Err(EngineError::invalid_input(
            "Receipt/output byte budget exceeded",
        ));
    }
    Ok(bytes)
}

fn specificity_multi_read_json<T: serde::de::DeserializeOwned>(
    path: &str,
) -> Result<T, EngineError> {
    serde_json::from_slice(&specificity_multi_read_bytes(path, 16 * 1024 * 1024)?)
        .map_err(|e| EngineError::invalid_input(format!("Invalid multi-reference JSON: {e}")))
}

fn specificity_multi_write_json<T: Serialize>(path: &Path, record: &T) -> Result<(), EngineError> {
    let bytes = serde_json::to_vec_pretty(record)
        .map_err(|e| EngineError::new(ErrorCode::Internal, e.to_string()))?;
    fs::write(path, bytes).map_err(|e| {
        EngineError::new(
            ErrorCode::Io,
            format!("Could not write '{}': {e}", path.display()),
        )
    })
}

// Validate raw components before normalization, including Windows verbatim paths.
fn specificity_multi_destination(raw: &str) -> Result<PathBuf, EngineError> {
    let path = Path::new(raw);
    let parent_component = |c: Component<'_>| c == Component::ParentDir || c.as_os_str() == "..";
    if raw.trim().is_empty()
        || raw.chars().any(char::is_control)
        || path.components().any(parent_component)
    {
        return Err(EngineError::invalid_input(
            "Unsafe multi-reference output path",
        ));
    }
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?
            .join(path)
    };
    let parent = specificity_multi_ordinary_parent(&absolute)?;
    let absolute = fs::canonicalize(parent)
        .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?
        .join(absolute.file_name().ok_or_else(|| {
            EngineError::invalid_input("Output must name a new bundle directory")
        })?);
    match fs::symlink_metadata(&absolute) {
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(absolute),
        Ok(_) => Err(EngineError::invalid_input(
            "Output directory already exists; evidence is never overwritten",
        )),
        Err(e) => Err(EngineError::new(ErrorCode::Io, e.to_string())),
    }
}

fn specificity_multi_ordinary_parent(path: &Path) -> Result<&Path, EngineError> {
    let parent = path
        .parent()
        .ok_or_else(|| EngineError::invalid_input("Path must name a file or bundle directory"))?;
    let mut ancestor = PathBuf::new();
    for component in parent.components() {
        if component == Component::ParentDir || component.as_os_str() == ".." {
            return Err(EngineError::invalid_input(
                "Output parent traversal is not permitted",
            ));
        }
        ancestor.push(component.as_os_str());
        if matches!(component, Component::Prefix(_) | Component::CurDir) {
            continue;
        }
        let metadata = fs::symlink_metadata(&ancestor)
            .map_err(|e| EngineError::new(ErrorCode::Io, e.to_string()))?;
        if !metadata.is_dir() || metadata.file_type().is_symlink() {
            return Err(EngineError::invalid_input(
                "Output ancestors must be ordinary existing directories",
            ));
        }
    }
    Ok(parent)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Hand-crafted synthetic oligos/references; no real gene or external resource.
    fn request() -> PrimerSpecificityMultiRequest {
        let primer = |role| PrimerSpecificityInputPrimer {
            role,
            full_sequence: "ACGTACGTACGTACGTACGT".into(),
            annealing_sequence: "ACGTACGTACGTACGTACGT".into(),
            annealing_length_bp: 20,
            ..Default::default()
        };
        PrimerSpecificityMultiRequest {
            schema: PRIMER_SPECIFICITY_MULTI_REQUEST_SCHEMA.into(),
            pair: PrimerSpecificityMultiPair::ExplicitPair {
                forward: primer(PrimerSpecificityPrimerRole::Forward),
                reverse: primer(PrimerSpecificityPrimerRole::Reverse),
            },
            policy: Default::default(),
            references: vec![PrimerSpecificityMultiReference {
                target_genome_id: "synthetic-a".into(),
                expected_index_kind: BlastDatabaseIndexKind::GenomicDna,
                required: true,
                intended_target: None,
            }],
            catalog_path: None,
            cache_dir: None,
        }
    }

    #[test]
    fn specificity_multi_admission_is_explicit_and_bounded() {
        let valid = request();
        GentleEngine::validate_specificity_multi_request(&valid).unwrap();
        for count in [0, 9] {
            let mut invalid = valid.clone();
            invalid.references = vec![valid.references[0].clone(); count];
            assert!(GentleEngine::validate_specificity_multi_request(&invalid).is_err());
        }
        let mut invalid = valid.clone();
        invalid.references[0].required = false;
        assert!(GentleEngine::validate_specificity_multi_request(&invalid).is_err());
        invalid = valid.clone();
        invalid.policy.specificity_target_genome_id = Some("another-reference".into());
        assert!(GentleEngine::validate_specificity_multi_request(&invalid).is_err());
        for id in ["", "all*", "?"] {
            invalid = valid.clone();
            invalid.references[0].target_genome_id = id.into();
            assert!(GentleEngine::validate_specificity_multi_request(&invalid).is_err());
        }
        let engine = GentleEngine::new();
        if let PrimerSpecificityMultiPair::ExplicitPair { forward, .. } = &mut invalid.pair {
            forward.non_annealing_5prime_tail_bp = 5;
        }
        assert!(engine.specificity_multi_input(&invalid.pair).is_err());
        let mut bytes = serde_json::to_value(&valid).unwrap();
        bytes["references"][0]["invented_mapping"] = json!(true);
        assert!(serde_json::from_value::<PrimerSpecificityMultiRequest>(bytes).is_err());
    }

    #[test]
    fn specificity_multi_required_preflight_leaves_no_files_or_state() {
        let temp = tempfile::tempdir().unwrap();
        let root = temp.path().canonicalize().unwrap();
        let catalog = root.join("catalog.json");
        fs::write(&catalog, serde_json::to_vec(&json!({"synthetic-a": {
            "description":"hand-crafted unavailable reference", "sequence_local":temp.path().join("not-installed.fa"),
            "annotations_local":temp.path().join("not-installed.gtf"), "cache_dir":temp.path().join("cache")
        }})).unwrap()).unwrap();
        let mut request = request();
        request.catalog_path = Some(catalog.to_string_lossy().into());
        let engine = GentleEngine::new();
        let before = serde_json::to_value(engine.snapshot()).unwrap();
        let output = root.join("new bundle");
        let error = engine
            .prepare_primer_pair_multi_reference_specificity_handoff(
                request,
                &output.to_string_lossy(),
            )
            .unwrap_err();
        assert!(error.message.contains("required_reference_unavailable"));
        assert!(!output.exists());
        assert_eq!(before, serde_json::to_value(engine.snapshot()).unwrap());
    }

    #[test]
    fn specificity_multi_output_paths_refuse_traversal_and_existing_evidence() {
        let temp = tempfile::tempdir().unwrap();
        let root = temp.path().canonicalize().unwrap();
        assert!(
            specificity_multi_destination(&temp.path().join("../escaped").to_string_lossy())
                .is_err()
        );
        assert!(specificity_multi_destination(&temp.path().to_string_lossy()).is_err());
        assert!(specificity_multi_destination(&root.join("new bundle").to_string_lossy()).is_ok());
        #[cfg(windows)]
        assert!(specificity_multi_destination(r"\\?\C:\existing\..\escaped").is_err());
        #[cfg(unix)]
        {
            let alias = root.join("alias");
            std::os::unix::fs::symlink(&root, &alias).unwrap();
            assert!(specificity_multi_destination(&alias.join("new").to_string_lossy()).is_err());
        }
    }

    #[test]
    fn specificity_multi_dimensions_never_vacuously_pass() {
        use PrimerSpecificityMultiVerdict::*;
        let row = |kind, required, verdict| PrimerSpecificityMultiSummaryReference {
            genome_id: "synthetic".into(),
            required,
            index_kind: kind,
            reference: None,
            verdict,
            applicability: "synthetic_test".into(),
            diagnostics: vec![],
            report: None,
        };
        for first in [Pass, Fail, Incomplete] {
            for second in [Pass, Fail, Incomplete] {
                let rows = vec![
                    row(BlastDatabaseIndexKind::GenomicDna, true, first),
                    row(BlastDatabaseIndexKind::GenomicDna, true, second),
                ];
                let expected = if [first, second].contains(&Fail) {
                    Fail
                } else if first == Pass && second == Pass {
                    Pass
                } else {
                    Incomplete
                };
                assert_eq!(
                    GentleEngine::specificity_multi_dimension(
                        &rows,
                        BlastDatabaseIndexKind::GenomicDna,
                        true
                    ),
                    expected
                );
                assert_eq!(
                    GentleEngine::specificity_multi_dimension(
                        &rows,
                        BlastDatabaseIndexKind::GenomicDna,
                        false
                    ),
                    if expected == Fail { Fail } else { Incomplete }
                );
            }
        }
        let rows = vec![
            row(BlastDatabaseIndexKind::GenomicDna, true, Pass),
            row(BlastDatabaseIndexKind::GenomicDna, false, Fail),
            row(BlastDatabaseIndexKind::TranscriptomeCdna, false, Pass),
        ];
        assert_eq!(
            GentleEngine::specificity_multi_dimension(
                &rows,
                BlastDatabaseIndexKind::GenomicDna,
                true
            ),
            Pass
        );
        assert_eq!(
            GentleEngine::specificity_multi_dimension(
                &rows,
                BlastDatabaseIndexKind::TranscriptomeCdna,
                true
            ),
            NotRequired
        );
        assert_eq!(
            GentleEngine::specificity_multi_dimension(
                &[],
                BlastDatabaseIndexKind::TranscriptomeCdna,
                true
            ),
            NotRequested
        );
    }

    #[test]
    fn specificity_multi_scientific_identity_excludes_paths_and_execution_ids() {
        let mut summary = PrimerSpecificityMultiSummary {
            schema: PRIMER_SPECIFICITY_MULTI_SUMMARY_SCHEMA.into(),
            summary_id: String::new(),
            content_sha256: String::new(),
            handoff_id: "synthetic-handoff".into(),
            pair_binding_sha256: "sha256:synthetic-pair".into(),
            source_snapshot_sha256: None,
            request: request(),
            execution_manifest: PrimerSpecificityMultiExecutionManifest {
                schema: PRIMER_SPECIFICITY_MULTI_MANIFEST_SCHEMA.into(),
                handoff_id: "synthetic-handoff".into(),
                handoff_content_sha256: "sha256:handoff".into(),
                pair_binding_sha256: "sha256:synthetic-pair".into(),
                commands: vec![],
            },
            execution_complete: false,
            references: vec![PrimerSpecificityMultiSummaryReference {
                genome_id: "synthetic-a".into(),
                required: true,
                index_kind: BlastDatabaseIndexKind::GenomicDna,
                reference: None,
                verdict: PrimerSpecificityMultiVerdict::Incomplete,
                applicability: "synthetic".into(),
                diagnostics: vec![],
                report: Some(PrimerSpecificityReport::default()),
            }],
            genomic: PrimerSpecificityMultiVerdict::Incomplete,
            transcriptome: PrimerSpecificityMultiVerdict::NotRequested,
            nonclaims: MULTI_NONCLAIMS.iter().map(|s| s.to_string()).collect(),
        };
        let digest = GentleEngine::specificity_multi_summary_digest(&summary).unwrap();
        summary.request.catalog_path = Some("another checkout/catalog.json".into());
        summary.request.cache_dir = Some("another cache".into());
        summary.handoff_id = "another-bundle".into();
        summary.execution_manifest.handoff_content_sha256 = "sha256:another-bundle".into();
        let report = summary.references[0].report.as_mut().unwrap();
        report.generated_at_unix_ms = 1234;
        report.op_id = Some("different-op".into());
        report.run_id = Some("different-run".into());
        report.catalog_path = Some("another checkout/catalog.json".into());
        report.cache_dir = Some("another cache".into());
        assert_eq!(
            digest,
            GentleEngine::specificity_multi_summary_digest(&summary).unwrap()
        );
        summary.source_snapshot_sha256 = Some("sha256:changed-template".into());
        assert_ne!(
            digest,
            GentleEngine::specificity_multi_summary_digest(&summary).unwrap()
        );
        summary.source_snapshot_sha256 = None;
        summary.request.policy.max_3prime_mismatches += 1;
        assert_ne!(
            digest,
            GentleEngine::specificity_multi_summary_digest(&summary).unwrap()
        );
    }

    #[test]
    fn specificity_multi_manifest_rejects_changed_and_duplicate_commands() {
        let command = PrimerSpecificityHandoffCommand {
            command_id: "synthetic:forward".into(),
            output_tsv_path: "output with spaces.tsv".into(),
            args: vec!["literal".into()],
            ..Default::default()
        };
        let handoff = PrimerSpecificityMultiHandoff {
            schema: PRIMER_SPECIFICITY_MULTI_HANDOFF_SCHEMA.into(),
            handoff_id: "bound".into(),
            content_sha256: "sha256:bound".into(),
            request: request(),
            pair_binding_sha256: "sha256:pair".into(),
            source_snapshot_sha256: None,
            primers: vec![],
            manifest_path: "manifest.json".into(),
            nonclaims: vec![],
            references: vec![PrimerSpecificityMultiHandoffReference {
                requested_genome_id: "synthetic-a".into(),
                resolved_genome_id: "synthetic-a".into(),
                required: true,
                expected_index_kind: BlastDatabaseIndexKind::GenomicDna,
                availability: PrimerSpecificityMultiAvailability::Prepared,
                diagnostic: None,
                reference: None,
                child: Some(PrimerSpecificityHandoff {
                    commands: vec![command.clone()],
                    ..Default::default()
                }),
            }],
        };
        let valid = PrimerSpecificityMultiExecutionManifest {
            schema: PRIMER_SPECIFICITY_MULTI_MANIFEST_SCHEMA.into(),
            handoff_id: handoff.handoff_id.clone(),
            handoff_content_sha256: handoff.content_sha256.clone(),
            pair_binding_sha256: handoff.pair_binding_sha256.clone(),
            commands: vec![PrimerSpecificityMultiExecutionCommand {
                command_id: command.command_id.clone(),
                command_sha256: GentleEngine::specificity_multi_command_digest(&command).unwrap(),
                output_path: command.output_tsv_path.clone(),
                state: PrimerSpecificityMultiExecutionState::Pending,
                exit_code: None,
                output_size_bytes: None,
                output_sha256: None,
            }],
        };
        GentleEngine::specificity_multi_validate_manifest(&handoff, &valid).unwrap();
        for field in ["handoff_content_sha256", "pair_binding_sha256"] {
            let mut value = serde_json::to_value(&valid).unwrap();
            value[field] = json!("sha256:tampered");
            assert!(
                GentleEngine::specificity_multi_validate_manifest(
                    &handoff,
                    &serde_json::from_value(value).unwrap()
                )
                .is_err()
            );
        }
        for field in ["command_id", "command_sha256", "output_path"] {
            let mut value = serde_json::to_value(&valid).unwrap();
            value["commands"][0][field] = json!("tampered");
            assert!(
                GentleEngine::specificity_multi_validate_manifest(
                    &handoff,
                    &serde_json::from_value(value).unwrap()
                )
                .is_err()
            );
        }
        let mut duplicated = valid.clone();
        duplicated.commands.push(valid.commands[0].clone());
        assert!(GentleEngine::specificity_multi_validate_manifest(&handoff, &duplicated).is_err());
    }

    #[test]
    fn specificity_multi_reads_retained_bytes_and_refuses_symlinks() {
        let temp = tempfile::tempdir().unwrap();
        let root = temp.path().canonicalize().unwrap();
        let path = root.join("rows with spaces.tsv");
        for newline in ["\n", "\r\n"] {
            let bytes = format!("q\tchr1\t100\t20\t0\t0\t1\t20\t100\t119\t1e-20\t80\t100{newline}")
                .into_bytes();
            fs::write(&path, &bytes).unwrap();
            let retained = specificity_multi_read_bytes(path.to_str().unwrap(), 1000).unwrap();
            assert_eq!(retained, bytes);
            let (hits, warnings) =
                crate::genomes::parse_blastn_tabular_hits(std::str::from_utf8(&retained).unwrap());
            assert_eq!(hits.len(), 1);
            assert!(warnings.is_empty());
            assert!(specificity_multi_read_bytes(path.to_str().unwrap(), 1).is_err());
        }
        #[cfg(unix)]
        {
            let alias = root.join("alias");
            std::os::unix::fs::symlink(&root, &alias).unwrap();
            assert!(
                specificity_multi_read_bytes(
                    alias.join("rows with spaces.tsv").to_str().unwrap(),
                    1000
                )
                .is_err()
            );
        }
    }

    #[cfg(unix)]
    #[test]
    fn specificity_multi_external_import_is_atomic_reference_bound_and_replayable() {
        use std::os::unix::fs::PermissionsExt;
        let _lock = crate::genomes::genbank_env_lock()
            .lock()
            .unwrap_or_else(|e| e.into_inner());
        let temp = tempfile::tempdir().unwrap();
        let root = temp.path().canonicalize().unwrap();
        // Entire fixture is hand-crafted: placeholder index files, fake tool
        // probes, invented HSPs. It checks contracts, never actual BLAST biology.
        let fasta = root.join("toy.fa");
        let annotation = root.join("toy.gtf");
        fs::write(&fasta, format!(">chr1\n{}\n", "ACGT".repeat(250))).unwrap();
        fs::write(
            &annotation,
            "chr1\tsynthetic\tgene\t1\t1000\t.\t+\t.\tgene_id \"TOY\";\n",
        )
        .unwrap();
        let catalog = root.join("catalog.json");
        let mut entries = serde_json::Map::new();
        for (id, kind) in [
            ("synthetic-a", "genomic_dna"),
            ("synthetic-b", "transcriptome_cdna"),
        ] {
            entries.insert(id.into(), json!({"description":"synthetic reference", "reference_name":"toy-assembly", "reference_release":"toy-release-1", "sequence_local":fasta, "annotations_local":annotation, "cache_dir":root.join(id), "blast_index_kind":kind}));
        }
        entries.insert("synthetic-optional".into(), json!({"description":"unavailable optional", "sequence_local":root.join("absent.fa"), "annotations_local":root.join("absent.gtf"), "cache_dir":root.join("optional")}));
        fs::write(&catalog, serde_json::to_vec(&entries).unwrap()).unwrap();
        let make = root.join("makeblastdb.sh");
        let blast = root.join("blastn.sh");
        let dbcmd = root.join("blastdbcmd.sh");
        fs::write(&make, "#!/bin/sh\nif [ \"$1\" = '-version' ]; then echo 'makeblastdb: synthetic'; exit 0; fi\nout=''\nwhile [ $# -gt 0 ]; do if [ \"$1\" = '-out' ]; then out=\"$2\"; shift 2; else shift; fi; done\nprintf nhr > \"${out}.nhr\"\nprintf nin > \"${out}.nin\"\nprintf nsq > \"${out}.nsq\"\n").unwrap();
        fs::write(&blast, "#!/bin/sh\nif [ \"$1\" = '-version' ]; then echo 'blastn: synthetic'; exit 0; fi\nexit 99\n").unwrap();
        fs::write(&dbcmd, "#!/bin/sh\nif [ \"$1\" = '-version' ]; then echo 'blastdbcmd: synthetic'; exit 0; fi\nif [ \"$3\" = '-info' ]; then printf 'Database: synthetic\\nBLASTDB Version: 5\\n 1 sequences; 1000 total letters\\n'; exit 0; fi\nexit 2\n").unwrap();
        for path in [&make, &blast, &dbcmd] {
            fs::set_permissions(path, fs::Permissions::from_mode(0o755)).unwrap();
        }
        let _make = crate::tool_overrides::ScopedToolOverrideGuard::set(
            crate::genomes::MAKEBLASTDB_ENV_BIN,
            make.to_str().unwrap(),
        );
        let _blast = crate::tool_overrides::ScopedToolOverrideGuard::set(
            crate::genomes::BLASTN_ENV_BIN,
            blast.to_str().unwrap(),
        );
        let _dbcmd = crate::tool_overrides::ScopedToolOverrideGuard::set(
            crate::genomes::BLASTDBCMD_ENV_BIN,
            dbcmd.to_str().unwrap(),
        );
        let mut engine = GentleEngine::new();
        for id in ["synthetic-a", "synthetic-b"] {
            engine
                .apply(Operation::PrepareGenome {
                    genome_id: id.into(),
                    catalog_path: Some(catalog.to_string_lossy().into()),
                    cache_dir: None,
                    timeout_seconds: None,
                })
                .unwrap();
        }
        let mut req = request();
        req.catalog_path = Some(catalog.to_string_lossy().into());
        req.policy.full_alignment.mode = PrimerSpecificityFullAlignmentMode::Disabled;
        req.references = [
            ("synthetic-a", BlastDatabaseIndexKind::GenomicDna, true),
            (
                "synthetic-b",
                BlastDatabaseIndexKind::TranscriptomeCdna,
                true,
            ),
            (
                "synthetic-optional",
                BlastDatabaseIndexKind::GenomicDna,
                false,
            ),
        ]
        .into_iter()
        .map(|(id, kind, required)| PrimerSpecificityMultiReference {
            target_genome_id: id.into(),
            expected_index_kind: kind,
            required,
            intended_target: Some(PrimerSpecificityIntendedTarget {
                model: if kind == BlastDatabaseIndexKind::GenomicDna {
                    PrimerSpecificityIntendedTargetModel::GenomicInterval
                } else {
                    PrimerSpecificityIntendedTargetModel::TranscriptSet
                },
                subject_id: Some("chr1".into()),
                expected_product_range: Some(PrimerSpecificitySubjectRange {
                    start_1based: 100,
                    end_1based: 219,
                }),
                forward_binding_ranges: vec![PrimerSpecificitySubjectRange {
                    start_1based: 100,
                    end_1based: 119,
                }],
                reverse_binding_ranges: vec![PrimerSpecificitySubjectRange {
                    start_1based: 200,
                    end_1based: 219,
                }],
                genomic_target_geometry_known: true,
                contiguous_genomic_product_expected: true,
                expected_products: vec![PrimerSpecificityExpectedProduct {
                    target_space: kind.as_str().into(),
                    subject_id: "chr1".into(),
                    expected_product_range: Some(PrimerSpecificitySubjectRange {
                        start_1based: 100,
                        end_1based: 219,
                    }),
                    source_transcript_id: None,
                }],
                ..Default::default()
            }),
        })
        .collect();
        let before = serde_json::to_value(engine.snapshot()).unwrap();
        let output = root.join("bundle with spaces");
        let handoff = engine
            .prepare_primer_pair_multi_reference_specificity_handoff(
                req.clone(),
                output.to_str().unwrap(),
            )
            .unwrap();
        engine.specificity_multi_validate_handoff(&handoff).unwrap();
        assert_eq!(before, serde_json::to_value(engine.snapshot()).unwrap());
        assert!(handoff.references.last().unwrap().child.is_none());
        let path = output.join("handoff.json");
        let manifest_path = Path::new(&handoff.manifest_path);
        let mut manifest: PrimerSpecificityMultiExecutionManifest =
            specificity_multi_read_json(&handoff.manifest_path).unwrap();
        let pending = engine
            .import_primer_pair_multi_reference_specificity(
                path.to_str().unwrap(),
                manifest_path.to_str().unwrap(),
            )
            .unwrap();
        assert_eq!(pending.genomic, PrimerSpecificityMultiVerdict::Incomplete);
        for child in handoff.references.iter().filter_map(|r| r.child.as_ref()) {
            for command in &child.commands {
                let (start, end) = if command.role == PrimerSpecificityPrimerRole::Forward {
                    (100, 119)
                } else {
                    (219, 200)
                };
                let bytes = format!(
                    "{}\tchr1\t100\t20\t0\t0\t1\t20\t{start}\t{end}\t1e-20\t80\t100\n",
                    command.query_label
                )
                .into_bytes();
                fs::write(&command.output_tsv_path, &bytes).unwrap();
                let execution = manifest
                    .commands
                    .iter_mut()
                    .find(|c| c.command_id == command.command_id)
                    .unwrap();
                execution.state = PrimerSpecificityMultiExecutionState::Completed;
                execution.exit_code = Some(0);
                execution.output_size_bytes = Some(bytes.len() as u64);
                execution.output_sha256 = Some(sha256_prefixed_bytes(&bytes));
            }
        }
        specificity_multi_write_json(manifest_path, &manifest).unwrap();
        let approved = engine
            .apply(Operation::ImportPrimerPairMultiReferenceSpecificity {
                handoff_path: path.to_string_lossy().into(),
                manifest_path: manifest_path.to_string_lossy().into(),
            })
            .unwrap()
            .primer_specificity_multi_summary
            .unwrap();
        assert_eq!(approved.genomic, PrimerSpecificityMultiVerdict::Pass);
        assert_eq!(approved.transcriptome, PrimerSpecificityMultiVerdict::Pass);
        assert_eq!(
            approved.references.last().unwrap().verdict,
            PrimerSpecificityMultiVerdict::Incomplete
        );
        let show = crate::engine_shell::parse_shell_tokens(&[
            "primers".into(),
            "specificity-multi-show".into(),
            approved.summary_id.clone(),
        ])
        .unwrap();
        let shown = crate::engine_shell::execute_shell_command(&mut engine, &show).unwrap();
        assert!(!shown.state_changed);
        assert_eq!(
            shown.output["result"]["primer_specificity_multi_summary"],
            serde_json::to_value(&approved).unwrap()
        );
        let raw = crate::engine_shell::ShellCommand::Op {
            payload: serde_json::to_string(
                &Operation::GetPrimerPairMultiReferenceSpecificitySummary {
                    summary_id: approved.summary_id.clone(),
                },
            )
            .unwrap(),
        };
        let raw = crate::engine_shell::execute_shell_command(&mut engine, &raw).unwrap();
        assert_eq!(
            raw.output["result"]["primer_specificity_multi_summary"],
            shown.output["result"]["primer_specificity_multi_summary"]
        );
        engine.undo_last_operation().unwrap();
        assert!(
            engine
                .get_primer_pair_multi_reference_specificity_summary(&approved.summary_id)
                .is_err()
        );
        engine.redo_last_operation().unwrap();
        let reloaded = GentleEngine::from_state(
            serde_json::from_value(serde_json::to_value(engine.snapshot()).unwrap()).unwrap(),
        );
        assert_eq!(
            reloaded
                .get_primer_pair_multi_reference_specificity_summary(&approved.summary_id)
                .unwrap()
                .content_sha256,
            approved.content_sha256
        );
        assert_eq!(
            engine
                .import_primer_pair_multi_reference_specificity(
                    path.to_str().unwrap(),
                    manifest_path.to_str().unwrap()
                )
                .unwrap()
                .summary_id,
            approved.summary_id
        );
        assert!(
            engine
                .read_primer_design_store()
                .active_primer_specificity_reference_selections
                .is_empty()
        );
        let missing_output =
            &handoff.references[0].child.as_ref().unwrap().commands[0].output_tsv_path;
        let retained = fs::read(missing_output).unwrap();
        fs::remove_file(missing_output).unwrap();
        let missing = engine
            .import_primer_pair_multi_reference_specificity(
                path.to_str().unwrap(),
                manifest_path.to_str().unwrap(),
            )
            .unwrap();
        assert!(!missing.execution_complete);
        assert_eq!(missing.genomic, PrimerSpecificityMultiVerdict::Incomplete);
        assert_eq!(
            missing.transcriptome,
            PrimerSpecificityMultiVerdict::Incomplete
        );
        assert!(missing.references[1].report.is_some());
        fs::write(missing_output, retained).unwrap();
        let atomic_before = serde_json::to_value(engine.snapshot()).unwrap();
        let mut tampered = manifest.clone();
        tampered.commands[0].output_sha256 = Some("sha256:wrong".into());
        specificity_multi_write_json(manifest_path, &tampered).unwrap();
        assert!(
            engine
                .import_primer_pair_multi_reference_specificity(
                    path.to_str().unwrap(),
                    manifest_path.to_str().unwrap()
                )
                .is_err()
        );
        assert_eq!(
            atomic_before,
            serde_json::to_value(engine.snapshot()).unwrap()
        );
        specificity_multi_write_json(manifest_path, &manifest).unwrap();
        let mut cancelled = manifest.clone();
        cancelled.commands[0].state = PrimerSpecificityMultiExecutionState::Cancelled;
        cancelled.commands[0].exit_code = None;
        specificity_multi_write_json(manifest_path, &cancelled).unwrap();
        let partial = engine
            .import_primer_pair_multi_reference_specificity(
                path.to_str().unwrap(),
                manifest_path.to_str().unwrap(),
            )
            .unwrap();
        assert!(!partial.execution_complete);
        assert_eq!(partial.genomic, PrimerSpecificityMultiVerdict::Incomplete);
        assert_eq!(
            partial.transcriptome,
            PrimerSpecificityMultiVerdict::Incomplete
        );
        assert!(
            partial.references[1].report.is_some(),
            "valid child evidence survives a partial run"
        );
        specificity_multi_write_json(manifest_path, &manifest).unwrap();
        for mutate in ["command", "pair", "fingerprint", "policy", "options"] {
            let mut altered = handoff.clone();
            match mutate {
                "command" => {
                    altered.references[0].child.as_mut().unwrap().commands[0].args[5] =
                        "blastn".into()
                }
                "pair" => altered.primers[0].non_annealing_5prime_tail_bp = 1,
                "fingerprint" => {
                    altered.references[0]
                        .child
                        .as_mut()
                        .unwrap()
                        .blast_database
                        .as_mut()
                        .unwrap()
                        .content_fingerprint = None
                }
                "options" => {
                    altered.references[0]
                        .child
                        .as_mut()
                        .unwrap()
                        .effective_blast_options
                        .as_mut()
                        .unwrap()
                        .max_hits += 1
                }
                _ => altered.request.policy.max_3prime_mismatches += 1,
            }
            altered.content_sha256 =
                GentleEngine::specificity_multi_handoff_digest(&altered).unwrap();
            altered.handoff_id = format!(
                "primer_multi_{}",
                altered.content_sha256.trim_start_matches("sha256:")
            );
            assert!(
                engine.specificity_multi_validate_handoff(&altered).is_err(),
                "{mutate}"
            );
        }
        if let Some(annotation_path) = handoff.references[0]
            .child
            .as_ref()
            .unwrap()
            .blast_database
            .as_ref()
            .unwrap()
            .subject_annotation_index_path
            .as_ref()
        {
            let old = fs::read(annotation_path).unwrap();
            let mut annotation: Value = serde_json::from_slice(&old).unwrap();
            annotation["synthetic_tamper"] = json!(true);
            fs::write(annotation_path, serde_json::to_vec(&annotation).unwrap()).unwrap();
            let changed = engine
                .import_primer_pair_multi_reference_specificity(
                    path.to_str().unwrap(),
                    manifest_path.to_str().unwrap(),
                )
                .unwrap();
            assert_eq!(changed.genomic, PrimerSpecificityMultiVerdict::Incomplete);
            fs::write(annotation_path, old).unwrap();
        }
        // Same prefix, different database bytes must be stale; display/list
        // remain historical and must not execute a probe to erase that history.
        let prefix = &handoff.references[0]
            .child
            .as_ref()
            .unwrap()
            .blast_db_prefix;
        fs::write(format!("{prefix}.nsq"), b"replaced-index").unwrap();
        let stale = engine
            .import_primer_pair_multi_reference_specificity(
                path.to_str().unwrap(),
                manifest_path.to_str().unwrap(),
            )
            .unwrap();
        assert_eq!(stale.genomic, PrimerSpecificityMultiVerdict::Incomplete);
        assert_eq!(
            stale.references[0].applicability,
            "stale_or_unavailable_reference"
        );
        assert_eq!(
            engine
                .get_primer_pair_multi_reference_specificity_summary(&approved.summary_id)
                .unwrap()
                .genomic,
            PrimerSpecificityMultiVerdict::Pass
        );
        assert!(
            engine
                .list_primer_pair_multi_reference_specificity_summaries()
                .len()
                >= 3
        );
    }

    #[test]
    fn specificity_multi_source_snapshot_detects_source_and_topology_changes() {
        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::CreateSequenceFromText {
                sequence_text: "ACGT".repeat(40),
                name: Some("synthetic".into()),
                output_id: Some("synthetic".into()),
                circular: false,
            })
            .unwrap();
        let mut input = engine.specificity_multi_input(&request().pair).unwrap();
        input.primary_seq_id = Some("synthetic".into());
        let original = engine.specificity_multi_source_snapshot(&input).unwrap();
        input.intended_target.source_reference = Some(PrimerSpecificitySourceReference {
            genome_id: "synthetic-a".into(),
            assembly: Some("toy".into()),
            release: Some("changed".into()),
            ..Default::default()
        });
        assert_ne!(
            original,
            engine.specificity_multi_source_snapshot(&input).unwrap()
        );
        input.intended_target.source_reference = None;
        engine
            .state
            .sequences
            .get_mut("synthetic")
            .unwrap()
            .set_circular(true);
        assert_ne!(
            original,
            engine.specificity_multi_source_snapshot(&input).unwrap()
        );
        engine.state.sequences.remove("synthetic");
        assert!(engine.specificity_multi_source_snapshot(&input).is_err());
    }
}
