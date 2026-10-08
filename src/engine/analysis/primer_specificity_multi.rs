//! Standalone, non-executing multi-reference primer specificity and receipt validation.
//! Biology stays in the existing single-reference interpreter; references never
//! become interchangeable, and this family never attaches panel readiness.

use super::operation_handlers::PrimerSpecificityResolvedInput;
use super::*;
use std::path::Component;

const MULTI_NONCLAIMS: &[&str] = &[
    "Specificity applies only to the explicitly bound prepared references and policy, not every genome or transcript.",
    "A hash binds retained content; it does not authenticate external process execution.",
    "This standalone summary does not establish isoform discrimination, oligo QC, laboratory validation or order readiness.",
    "Current applicability is checked at import; historical display performs no database probes.",
];

impl GentleEngine {
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
    let parent = absolute
        .parent()
        .ok_or_else(|| EngineError::invalid_input("Output must name a new bundle directory"))?;
    let mut ancestor = PathBuf::new();
    for component in parent.components() {
        if parent_component(component) {
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
    match fs::symlink_metadata(&absolute) {
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(absolute),
        Ok(_) => Err(EngineError::invalid_input(
            "Output directory already exists; evidence is never overwritten",
        )),
        Err(e) => Err(EngineError::new(ErrorCode::Io, e.to_string())),
    }
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
        let catalog = temp.path().join("catalog.json");
        fs::write(&catalog, serde_json::to_vec(&json!({"synthetic-a": {
            "description":"hand-crafted unavailable reference", "sequence_local":temp.path().join("not-installed.fa"),
            "annotations_local":temp.path().join("not-installed.gtf"), "cache_dir":temp.path().join("cache")
        }})).unwrap()).unwrap();
        let mut request = request();
        request.catalog_path = Some(catalog.to_string_lossy().into());
        let engine = GentleEngine::new();
        let before = serde_json::to_value(engine.snapshot()).unwrap();
        let output = temp.path().join("new bundle");
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
        assert!(
            specificity_multi_destination(&temp.path().join("../escaped").to_string_lossy())
                .is_err()
        );
        assert!(specificity_multi_destination(&temp.path().to_string_lossy()).is_err());
        assert!(
            specificity_multi_destination(&temp.path().join("new bundle").to_string_lossy())
                .is_ok()
        );
        #[cfg(windows)]
        assert!(specificity_multi_destination(r"\\?\C:\existing\..\escaped").is_err());
        #[cfg(unix)]
        {
            let alias = temp.path().join("alias");
            std::os::unix::fs::symlink(temp.path(), &alias).unwrap();
            assert!(specificity_multi_destination(&alias.join("new").to_string_lossy()).is_err());
        }
    }
}
