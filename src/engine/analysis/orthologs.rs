//! Ortholog promoter cohort resolution and comparison.
//!
//! This slice is intentionally offline-first: local ortholog resources map
//! genes across species, prepared genome catalogs resolve promoter windows, and
//! evidence summaries remain separated by evidence type.

use super::*;
use gentle_protocol::{
    BiologicalContext, BiologicalContextRegistry, GeneSetCohortRelationship,
    GeneSetCohortRelationshipFlag, ORTHOLOG_PROMOTER_COHORT_SCHEMA,
    ORTHOLOG_PROMOTER_COMPARISON_SCHEMA, ORTHOLOG_RESOURCE_SCHEMA, OrthologAmbiguityCandidate,
    OrthologConfidence, OrthologCutRunSupportRow, OrthologCutRunSupportStatus,
    OrthologExpressionAssignment, OrthologMappingRow, OrthologPairwiseTfbsSimilarity,
    OrthologPromoterCohortReport, OrthologPromoterCohortRequest, OrthologPromoterComparisonReport,
    OrthologPromoterRole, OrthologPromoterRow, OrthologResource, OrthologSequenceSimilarityRow,
    OrthologTfbsPeakSummary, OrthologTfbsSummaryRow, OrthologUnresolvedRow, OrthologyType,
};

#[derive(Debug, Clone)]
struct OrientedOrthologMapping {
    source_context_id: Option<String>,
    target_species: String,
    target_context_id: Option<String>,
    target_gene_id: Option<String>,
    target_gene_symbol: Option<String>,
    orthology_type: Option<OrthologyType>,
    confidence: Option<OrthologConfidence>,
    source: Option<String>,
    evidence: Vec<String>,
}

impl GentleEngine {
    pub(crate) fn load_ortholog_resource(path: &str) -> Result<OrthologResource, EngineError> {
        let raw = std::fs::read_to_string(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read ortholog resource '{}': {}", path, e),
            cause_chain: vec![],
        })?;
        let mut resource: OrthologResource =
            serde_json::from_str(&raw).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not parse ortholog resource '{}': {}", path, e),
                cause_chain: vec![],
            })?;
        if resource.schema.trim().is_empty() {
            resource.schema = ORTHOLOG_RESOURCE_SCHEMA.to_string();
        }
        resource
            .validate_context_references()
            .map_err(|err| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Ortholog resource '{}' has invalid biological-context references: {}",
                    path, err
                ),
                cause_chain: vec![],
            })?;
        Ok(resource)
    }

    fn ortholog_species_key(raw: &str) -> String {
        raw.chars()
            .filter(|ch| ch.is_ascii_alphanumeric())
            .flat_map(char::to_lowercase)
            .collect()
    }

    fn ortholog_gene_key(raw: &str) -> String {
        Self::normalize_id_token(raw)
    }

    fn ortholog_species_alias_map(resource: &OrthologResource) -> BTreeMap<String, String> {
        let mut map = BTreeMap::new();
        for alias in &resource.species_aliases {
            let species = alias.species.trim();
            if species.is_empty() {
                continue;
            }
            map.insert(Self::ortholog_species_key(species), species.to_string());
            for raw_alias in &alias.aliases {
                let raw_alias = raw_alias.trim();
                if !raw_alias.is_empty() {
                    map.insert(Self::ortholog_species_key(raw_alias), species.to_string());
                }
            }
        }
        for row in &resource.rows {
            for species in [&row.source_species, &row.target_species] {
                let species = species.trim();
                if !species.is_empty() {
                    map.entry(Self::ortholog_species_key(species))
                        .or_insert_with(|| species.to_string());
                }
            }
        }
        map
    }

    fn ortholog_canonical_species(raw: &str, aliases: &BTreeMap<String, String>) -> String {
        aliases
            .get(&Self::ortholog_species_key(raw))
            .cloned()
            .unwrap_or_else(|| raw.trim().to_string())
    }

    fn validate_ortholog_resource_context_species(
        resource: &OrthologResource,
        aliases: &BTreeMap<String, String>,
    ) -> Result<(), EngineError> {
        for (row_index, row) in resource.rows.iter().enumerate() {
            for (side, species, context_id) in [
                (
                    "source",
                    row.source_species.as_str(),
                    row.source_context_id.as_deref(),
                ),
                (
                    "target",
                    row.target_species.as_str(),
                    row.target_context_id.as_deref(),
                ),
            ] {
                let Some(context_id) = context_id else {
                    continue;
                };
                let context = resource
                    .biological_contexts
                    .context(context_id)
                    .map_err(|err| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Ortholog row {} has invalid {} context '{}': {}",
                            row_index + 1,
                            side,
                            context_id,
                            err
                        ),
                        cause_chain: vec![],
                    })?;
                let Some(context_organism) = context.organism.as_deref() else {
                    continue;
                };
                let row_species = Self::ortholog_canonical_species(species, aliases);
                let context_species = Self::ortholog_canonical_species(context_organism, aliases);
                if row_species != context_species {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Ortholog row {} {} species '{}' conflicts with biological context '{}' organism '{}'",
                            row_index + 1,
                            side,
                            species,
                            context_id,
                            context_organism
                        ),
                        cause_chain: vec![],
                    });
                }
            }
        }
        Ok(())
    }

    fn validate_ortholog_context_genome(
        resource: &OrthologResource,
        context_id: Option<&str>,
        expected_genome_id: &str,
        role: &str,
    ) -> Result<(), EngineError> {
        let Some(context_id) = context_id else {
            return Ok(());
        };
        let context = resource
            .biological_contexts
            .context(context_id)
            .map_err(|err| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Invalid {role} context '{context_id}': {err}"),
                cause_chain: vec![],
            })?;
        if let Some(context_genome_id) = context.genome_id.as_deref()
            && context_genome_id.trim() != expected_genome_id.trim()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Ortholog {role} context '{}' declares genome '{}' but the request uses '{}'",
                    context_id, context_genome_id, expected_genome_id
                ),
                cause_chain: vec![],
            });
        }
        Ok(())
    }

    fn ortholog_runtime_context_id(species: &str, genome_id: &str) -> String {
        format!(
            "ortholog_{}_{}",
            Self::ortholog_species_key(species),
            Self::ortholog_gene_key(genome_id)
        )
    }

    fn register_ortholog_report_context(
        report_contexts: &mut BiologicalContextRegistry,
        resource: &OrthologResource,
        explicit_context_id: Option<&str>,
        species: &str,
        genome_id: &str,
    ) -> Result<String, EngineError> {
        Self::validate_ortholog_context_genome(
            resource,
            explicit_context_id,
            genome_id,
            "mapping",
        )?;
        let context_id = explicit_context_id
            .map(str::to_string)
            .unwrap_or_else(|| Self::ortholog_runtime_context_id(species, genome_id));
        let mut context = if let Some(explicit_context_id) = explicit_context_id {
            resource
                .biological_contexts
                .context(explicit_context_id)
                .map_err(|err| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Invalid ortholog biological context '{}': {}",
                        explicit_context_id, err
                    ),
                    cause_chain: vec![],
                })?
                .clone()
        } else {
            BiologicalContext {
                context_id: context_id.clone(),
                ..BiologicalContext::default()
            }
        };
        if context.organism.is_none() {
            context.organism = Some(species.to_string());
        }
        if context.genome_id.is_none() {
            context.genome_id = Some(genome_id.to_string());
        }

        if let Some(existing) = report_contexts
            .contexts
            .iter()
            .find(|existing| existing.context_id == context_id)
        {
            if !existing.semantically_matches(&context) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Ortholog context id '{}' resolves to conflicting biological contexts",
                        context_id
                    ),
                    cause_chain: vec![],
                });
            }
        } else {
            report_contexts.contexts.push(context);
            report_contexts
                .contexts
                .sort_by(|left, right| left.context_id.cmp(&right.context_id));
        }
        Ok(context_id)
    }

    fn ortholog_gene_matches(
        query_key: &str,
        gene_id: Option<&str>,
        gene_symbol: Option<&str>,
    ) -> bool {
        gene_id
            .map(Self::ortholog_gene_key)
            .filter(|key| !key.is_empty() && key == query_key)
            .is_some()
            || gene_symbol
                .map(Self::ortholog_gene_key)
                .filter(|key| !key.is_empty() && key == query_key)
                .is_some()
    }

    fn orient_ortholog_mapping(
        row: &OrthologMappingRow,
        anchor_species: &str,
        anchor_gene_query: &str,
        target_species: &str,
        aliases: &BTreeMap<String, String>,
    ) -> Option<OrientedOrthologMapping> {
        let anchor_species = Self::ortholog_canonical_species(anchor_species, aliases);
        let target_species = Self::ortholog_canonical_species(target_species, aliases);
        let source_species = Self::ortholog_canonical_species(&row.source_species, aliases);
        let target_row_species = Self::ortholog_canonical_species(&row.target_species, aliases);
        let anchor_key = Self::ortholog_gene_key(anchor_gene_query);
        if anchor_key.is_empty() {
            return None;
        }
        if source_species == anchor_species
            && target_row_species == target_species
            && Self::ortholog_gene_matches(
                &anchor_key,
                row.source_gene_id.as_deref(),
                row.source_gene_symbol.as_deref(),
            )
        {
            return Some(OrientedOrthologMapping {
                source_context_id: row.source_context_id.clone(),
                target_species: target_row_species,
                target_context_id: row.target_context_id.clone(),
                target_gene_id: row.target_gene_id.clone(),
                target_gene_symbol: row.target_gene_symbol.clone(),
                orthology_type: row.orthology_type.clone(),
                confidence: row.confidence.clone(),
                source: row.source.clone(),
                evidence: row.evidence.clone(),
            });
        }
        if target_row_species == anchor_species
            && source_species == target_species
            && Self::ortholog_gene_matches(
                &anchor_key,
                row.target_gene_id.as_deref(),
                row.target_gene_symbol.as_deref(),
            )
        {
            return Some(OrientedOrthologMapping {
                source_context_id: row.target_context_id.clone(),
                target_species: source_species,
                target_context_id: row.source_context_id.clone(),
                target_gene_id: row.source_gene_id.clone(),
                target_gene_symbol: row.source_gene_symbol.clone(),
                orthology_type: row.orthology_type.clone(),
                confidence: row.confidence.clone(),
                source: row.source.clone(),
                evidence: row.evidence.clone(),
            });
        }
        None
    }

    fn ortholog_candidate_label(candidate: &OrientedOrthologMapping) -> String {
        let id = candidate
            .target_gene_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("-");
        let symbol = candidate
            .target_gene_symbol
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("-");
        format!(
            "{}:{}:{}:{}",
            candidate.target_species,
            symbol,
            id,
            candidate.orthology_type.as_deref().unwrap_or("-")
        )
    }

    fn ortholog_candidate_target_genome_id(
        resource: &OrthologResource,
        candidate: &OrientedOrthologMapping,
        requested_target_genome_id: Option<&str>,
    ) -> (String, bool) {
        if let Some(genome_id) = requested_target_genome_id {
            return (genome_id.to_string(), false);
        }
        if let Some(genome_id) = candidate
            .target_context_id
            .as_deref()
            .and_then(|context_id| resource.biological_contexts.context(context_id).ok())
            .and_then(|context| context.genome_id.clone())
        {
            return (genome_id, false);
        }
        (candidate.target_species.clone(), true)
    }

    fn lookup_ortholog_species_map_value(
        map: &BTreeMap<String, String>,
        species: &str,
        aliases: &BTreeMap<String, String>,
    ) -> Option<String> {
        let target_key = Self::ortholog_species_key(species);
        map.iter().find_map(|(key, value)| {
            let canonical_key =
                Self::ortholog_species_key(&Self::ortholog_canonical_species(key, aliases));
            (Self::ortholog_species_key(key) == target_key || canonical_key == target_key)
                .then(|| value.clone())
        })
    }

    fn lookup_ortholog_transcript_id(
        map: &BTreeMap<String, String>,
        species: &str,
        aliases: &BTreeMap<String, String>,
    ) -> Option<String> {
        Self::lookup_ortholog_species_map_value(map, species, aliases)
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
    }

    fn ortholog_target_gene_query(candidate: &OrientedOrthologMapping) -> Option<String> {
        candidate
            .target_gene_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .or_else(|| {
                candidate
                    .target_gene_symbol
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(str::to_string)
            })
    }

    fn resolve_ortholog_promoter_row(
        catalog: &GenomeCatalog,
        species: &str,
        genome_id: &str,
        context_id: Option<&str>,
        role: OrthologPromoterRole,
        gene_query: &str,
        transcript_id: Option<&str>,
        upstream_bp: usize,
        downstream_bp: usize,
        cache_dir: Option<&str>,
        orthology: Option<&OrientedOrthologMapping>,
    ) -> Result<OrthologPromoterRow, EngineError> {
        let resolved = Self::resolve_genome_promoter_slice_request(
            catalog,
            genome_id,
            gene_query,
            None,
            transcript_id,
            upstream_bp,
            downstream_bp,
            cache_dir,
        )?;
        let strand = resolved
            .selected_gene
            .strand
            .or(resolved.selected_transcript.strand)
            .unwrap_or('+');
        let promoter_sequence = catalog
            .get_sequence_region_with_cache(
                genome_id,
                &resolved.selected_transcript.chromosome,
                resolved.extract_start_1based,
                resolved.extract_end_1based,
                cache_dir,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Could not load promoter slice {}:{}-{} from '{}': {}",
                    resolved.selected_transcript.chromosome,
                    resolved.extract_start_1based,
                    resolved.extract_end_1based,
                    genome_id,
                    e
                ),
                cause_chain: vec![],
            })?;
        let oriented_sequence = Self::promoter_aligned_sequence(&promoter_sequence, Some(strand));
        let promoter_length_bp = oriented_sequence.len();
        let tss_position_0based = Self::promoter_oriented_tss_position_0based(
            promoter_length_bp,
            resolved.extract_start_1based,
            resolved.tss_1based,
            Some(strand),
        )
        .ok_or_else(|| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not place TSS for gene '{}' into promoter-oriented coordinates",
                Self::genome_gene_display_label(&resolved.selected_gene)
            ),
            cause_chain: vec![],
        })?;

        let mut warnings = resolved.warnings.clone();
        let desired_promoter_length = upstream_bp.saturating_add(downstream_bp).saturating_add(1);
        if promoter_length_bp != desired_promoter_length {
            warnings.push(format!(
                "Promoter slice for '{}' in species '{}' is {} bp instead of requested {} bp, likely due to contig-boundary clipping.",
                Self::genome_gene_display_label(&resolved.selected_gene),
                species,
                promoter_length_bp,
                desired_promoter_length
            ));
        }
        let gene_symbol = resolved.selected_gene.gene_name.clone();
        let display_label = gene_symbol
            .clone()
            .unwrap_or_else(|| resolved.query.clone());
        Ok(OrthologPromoterRow {
            species: species.to_string(),
            genome_id: genome_id.to_string(),
            context_id: context_id.map(str::to_string),
            role,
            gene_query: resolved.query.clone(),
            gene_symbol,
            gene_id: resolved.selected_gene.gene_id.clone(),
            transcript_id_requested: transcript_id.map(str::to_string),
            transcript_id: resolved.selected_transcript.transcript_id.clone(),
            display_label,
            chromosome: resolved.selected_transcript.chromosome.clone(),
            strand: strand.to_string(),
            promoter_start_1based: resolved.extract_start_1based,
            promoter_end_1based: resolved.extract_end_1based,
            promoter_length_bp,
            tss_1based: resolved.tss_1based,
            tss_position_0based,
            sequence_orientation: "transcription_aligned".to_string(),
            promoter_sequence: Some(oriented_sequence),
            orthology_type: orthology.and_then(|row| row.orthology_type.clone()),
            confidence: orthology.and_then(|row| row.confidence.clone()),
            orthology_source_context_id: orthology.and_then(|row| row.source_context_id.clone()),
            orthology_target_context_id: orthology.and_then(|row| row.target_context_id.clone()),
            orthology_source: orthology.and_then(|row| row.source.clone()),
            orthology_evidence: orthology
                .map(|row| row.evidence.clone())
                .unwrap_or_default(),
            warnings,
        })
    }

    pub(crate) fn resolve_ortholog_promoter_cohort(
        &self,
        anchor_species: &str,
        anchor_genome_id: &str,
        anchor_gene_query: &str,
        target_species: &[String],
        target_genome_ids: &BTreeMap<String, String>,
        transcript_ids: &BTreeMap<String, String>,
        ortholog_resource_path: &str,
        upstream_bp: usize,
        downstream_bp: usize,
        ambiguity_policy: OrthologAmbiguityPolicy,
        relationship: GeneSetCohortRelationship,
        genome_catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<OrthologPromoterCohortReport, EngineError> {
        let trimmed_anchor_species = anchor_species.trim();
        let trimmed_anchor_genome_id = anchor_genome_id.trim();
        let trimmed_anchor_gene_query = anchor_gene_query.trim();
        if trimmed_anchor_species.is_empty()
            || trimmed_anchor_genome_id.is_empty()
            || trimmed_anchor_gene_query.is_empty()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "ResolveOrthologPromoterCohort requires anchor species, genome id, and gene query"
                        .to_string(),
                cause_chain: vec![],
            });
        }
        if target_species.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ResolveOrthologPromoterCohort requires at least one target species"
                    .to_string(),
                cause_chain: vec![],
            });
        }

        let resource = Self::load_ortholog_resource(ortholog_resource_path)?;
        let aliases = Self::ortholog_species_alias_map(&resource);
        Self::validate_ortholog_resource_context_species(&resource, &aliases)?;
        let anchor_species = Self::ortholog_canonical_species(trimmed_anchor_species, &aliases);
        let mut anchor_context_ids = BTreeSet::new();
        for raw_target_species in target_species {
            let canonical_target_species =
                Self::ortholog_canonical_species(raw_target_species, &aliases);
            let requested_target_genome_id = Self::lookup_ortholog_species_map_value(
                target_genome_ids,
                raw_target_species,
                &aliases,
            );
            for candidate in resource.rows.iter().filter_map(|row| {
                Self::orient_ortholog_mapping(
                    row,
                    &anchor_species,
                    trimmed_anchor_gene_query,
                    &canonical_target_species,
                    &aliases,
                )
            }) {
                if let Some(context_id) = candidate.source_context_id.as_deref() {
                    anchor_context_ids.insert(context_id.to_string());
                }
                Self::validate_ortholog_context_genome(
                    &resource,
                    candidate.source_context_id.as_deref(),
                    trimmed_anchor_genome_id,
                    "source",
                )?;
                if let Some(target_genome_id) = requested_target_genome_id.as_deref() {
                    Self::validate_ortholog_context_genome(
                        &resource,
                        candidate.target_context_id.as_deref(),
                        target_genome_id,
                        "target",
                    )?;
                }
            }
        }
        let mut biological_contexts = BiologicalContextRegistry::default();
        let anchor_explicit_context_id = if anchor_context_ids.len() == 1 {
            anchor_context_ids.first().map(String::as_str)
        } else {
            None
        };
        let anchor_context_id = Self::register_ortholog_report_context(
            &mut biological_contexts,
            &resource,
            anchor_explicit_context_id,
            &anchor_species,
            trimmed_anchor_genome_id,
        )?;
        let effective_catalog_path =
            genome_catalog_path.unwrap_or(crate::genomes::default_catalog_discovery_token(false));
        let (catalog, _) = Self::open_reference_genome_catalog(Some(effective_catalog_path))?;
        let mut warnings = resource.warnings.clone();
        if resource.schema != ORTHOLOG_RESOURCE_SCHEMA {
            warnings.push(format!(
                "Ortholog resource '{}' declares schema '{}', expected '{}'",
                ortholog_resource_path, resource.schema, ORTHOLOG_RESOURCE_SCHEMA
            ));
        }

        let anchor_transcript =
            Self::lookup_ortholog_transcript_id(transcript_ids, &anchor_species, &aliases);
        let mut rows = vec![Self::resolve_ortholog_promoter_row(
            &catalog,
            &anchor_species,
            trimmed_anchor_genome_id,
            Some(&anchor_context_id),
            OrthologPromoterRole::Anchor,
            trimmed_anchor_gene_query,
            anchor_transcript.as_deref(),
            upstream_bp,
            downstream_bp,
            cache_dir,
            None,
        )?];
        warnings.extend(rows[0].warnings.clone());

        let mut unresolved_rows = vec![];
        for raw_target_species in target_species {
            let canonical_target_species =
                Self::ortholog_canonical_species(raw_target_species, &aliases);
            let requested_target_genome_id = Self::lookup_ortholog_species_map_value(
                target_genome_ids,
                raw_target_species,
                &aliases,
            );
            let mut candidates = resource
                .rows
                .iter()
                .filter_map(|row| {
                    Self::orient_ortholog_mapping(
                        row,
                        &anchor_species,
                        trimmed_anchor_gene_query,
                        &canonical_target_species,
                        &aliases,
                    )
                })
                .collect::<Vec<_>>();
            candidates.sort_by(|left, right| {
                Self::ortholog_candidate_label(left).cmp(&Self::ortholog_candidate_label(right))
            });
            if candidates.len() > 1 && ambiguity_policy == OrthologAmbiguityPolicy::Preserve {
                let labels = candidates
                    .iter()
                    .map(Self::ortholog_candidate_label)
                    .collect::<Vec<_>>();
                let mut candidate_mappings = Vec::with_capacity(candidates.len());
                let mut candidate_context_ids = BTreeSet::new();
                let mut candidate_genome_ids = BTreeSet::new();
                let mut used_species_as_genome_id = false;
                for (candidate_index, candidate) in candidates.iter().enumerate() {
                    let (target_genome_id, used_species_fallback) =
                        Self::ortholog_candidate_target_genome_id(
                            &resource,
                            candidate,
                            requested_target_genome_id.as_deref(),
                        );
                    used_species_as_genome_id |= used_species_fallback;
                    let source_context_id = Self::register_ortholog_report_context(
                        &mut biological_contexts,
                        &resource,
                        candidate.source_context_id.as_deref(),
                        &anchor_species,
                        trimmed_anchor_genome_id,
                    )?;
                    let target_context_id = Self::register_ortholog_report_context(
                        &mut biological_contexts,
                        &resource,
                        candidate.target_context_id.as_deref(),
                        &candidate.target_species,
                        &target_genome_id,
                    )?;
                    candidate_context_ids.insert(target_context_id.clone());
                    candidate_genome_ids.insert(target_genome_id.clone());
                    candidate_mappings.push(OrthologAmbiguityCandidate {
                        candidate_rank: candidate_index + 1,
                        candidate_label: labels[candidate_index].clone(),
                        target_species: candidate.target_species.clone(),
                        target_genome_id: Some(target_genome_id),
                        source_context_id: Some(source_context_id),
                        target_context_id: Some(target_context_id),
                        target_gene_id: candidate.target_gene_id.clone(),
                        target_gene_symbol: candidate.target_gene_symbol.clone(),
                        orthology_type: candidate.orthology_type.clone(),
                        confidence: candidate.confidence.clone(),
                        source: candidate.source.clone(),
                        evidence: candidate.evidence.clone(),
                    });
                }
                if used_species_as_genome_id {
                    warnings.push(format!(
                        "No target genome id or context genome was provided for ambiguous species '{}'; preserved candidate mappings use the species label as genome id",
                        canonical_target_species
                    ));
                }
                let reason = format!(
                    "Ambiguous local ortholog mapping from '{}' '{}' to '{}' ({} candidate(s)); preserving all candidates without selecting a promoter because ambiguity_policy=preserve",
                    anchor_species,
                    trimmed_anchor_gene_query,
                    canonical_target_species,
                    labels.len()
                );
                warnings.push(reason.clone());
                let context_id = if candidate_context_ids.len() == 1 {
                    candidate_context_ids.first().cloned()
                } else {
                    None
                };
                let genome_id = if candidate_genome_ids.len() == 1 {
                    candidate_genome_ids.first().cloned()
                } else {
                    None
                };
                unresolved_rows.push(OrthologUnresolvedRow {
                    species: canonical_target_species,
                    context_id,
                    genome_id,
                    gene_query: None,
                    reason,
                    candidates: labels,
                    candidate_mappings,
                });
                continue;
            }
            let mut candidate = match candidates.len() {
                0 => {
                    let reason = format!(
                        "No local ortholog mapping from '{}' '{}' to '{}'",
                        anchor_species, trimmed_anchor_gene_query, canonical_target_species
                    );
                    warnings.push(reason.clone());
                    unresolved_rows.push(OrthologUnresolvedRow {
                        species: canonical_target_species,
                        context_id: None,
                        genome_id: requested_target_genome_id,
                        gene_query: None,
                        reason,
                        candidates: vec![],
                        candidate_mappings: vec![],
                    });
                    continue;
                }
                1 => candidates.remove(0),
                _ if ambiguity_policy == OrthologAmbiguityPolicy::First => {
                    let labels = candidates
                        .iter()
                        .map(Self::ortholog_candidate_label)
                        .collect::<Vec<_>>();
                    warnings.push(format!(
                        "Ortholog mapping from '{}' '{}' to '{}' is ambiguous ({} candidate(s)); using first because ambiguity_policy=first: {}",
                        anchor_species,
                        trimmed_anchor_gene_query,
                        canonical_target_species,
                        labels.len(),
                        labels.join(", ")
                    ));
                    candidates.remove(0)
                }
                _ => {
                    let labels = candidates
                        .iter()
                        .map(Self::ortholog_candidate_label)
                        .collect::<Vec<_>>();
                    let reason = format!(
                        "Ambiguous local ortholog mapping from '{}' '{}' to '{}' ({} candidate(s)); set ambiguity_policy=first to choose deterministically",
                        anchor_species,
                        trimmed_anchor_gene_query,
                        canonical_target_species,
                        labels.len()
                    );
                    warnings.push(reason.clone());
                    unresolved_rows.push(OrthologUnresolvedRow {
                        species: canonical_target_species,
                        context_id: None,
                        genome_id: requested_target_genome_id,
                        gene_query: None,
                        reason,
                        candidates: labels,
                        candidate_mappings: vec![],
                    });
                    continue;
                }
            };
            let (target_genome_id, used_species_as_genome_id) =
                Self::ortholog_candidate_target_genome_id(
                    &resource,
                    &candidate,
                    requested_target_genome_id.as_deref(),
                );
            if used_species_as_genome_id {
                warnings.push(format!(
                    "No target genome id provided for species '{}'; using species label as genome id",
                    candidate.target_species
                ));
            }
            let source_relation_context_id = Self::register_ortholog_report_context(
                &mut biological_contexts,
                &resource,
                candidate.source_context_id.as_deref(),
                &anchor_species,
                trimmed_anchor_genome_id,
            )?;
            let target_context_id = Self::register_ortholog_report_context(
                &mut biological_contexts,
                &resource,
                candidate.target_context_id.as_deref(),
                &candidate.target_species,
                &target_genome_id,
            )?;
            candidate.source_context_id = Some(source_relation_context_id);
            candidate.target_context_id = Some(target_context_id.clone());
            let Some(gene_query) = Self::ortholog_target_gene_query(&candidate) else {
                let reason = format!(
                    "Ortholog mapping to '{}' has neither target_gene_id nor target_gene_symbol",
                    candidate.target_species
                );
                warnings.push(reason.clone());
                unresolved_rows.push(OrthologUnresolvedRow {
                    species: candidate.target_species.clone(),
                    context_id: Some(target_context_id),
                    genome_id: Some(target_genome_id),
                    gene_query: None,
                    reason,
                    candidates: vec![Self::ortholog_candidate_label(&candidate)],
                    candidate_mappings: vec![],
                });
                continue;
            };
            let target_transcript = Self::lookup_ortholog_transcript_id(
                transcript_ids,
                &candidate.target_species,
                &aliases,
            );
            match Self::resolve_ortholog_promoter_row(
                &catalog,
                &candidate.target_species,
                &target_genome_id,
                Some(&target_context_id),
                OrthologPromoterRole::Target,
                &gene_query,
                target_transcript.as_deref(),
                upstream_bp,
                downstream_bp,
                cache_dir,
                Some(&candidate),
            ) {
                Ok(row) => {
                    warnings.extend(row.warnings.clone());
                    rows.push(row);
                }
                Err(err) => {
                    let reason = format!(
                        "Could not resolve promoter for '{}' ortholog '{}': {}",
                        candidate.target_species, gene_query, err.message
                    );
                    warnings.push(reason.clone());
                    unresolved_rows.push(OrthologUnresolvedRow {
                        species: candidate.target_species.clone(),
                        context_id: Some(target_context_id),
                        genome_id: Some(target_genome_id),
                        gene_query: Some(gene_query),
                        reason,
                        candidates: vec![Self::ortholog_candidate_label(&candidate)],
                        candidate_mappings: vec![],
                    });
                }
            }
        }

        rows.sort_by(|left, right| {
            (left.role != OrthologPromoterRole::Anchor)
                .cmp(&(right.role != OrthologPromoterRole::Anchor))
                .then(left.species.cmp(&right.species))
                .then(left.display_label.cmp(&right.display_label))
        });
        let request = OrthologPromoterCohortRequest {
            anchor_species,
            anchor_genome_id: trimmed_anchor_genome_id.to_string(),
            anchor_gene_query: trimmed_anchor_gene_query.to_string(),
            target_species: target_species
                .iter()
                .map(|value| Self::ortholog_canonical_species(value, &aliases))
                .collect(),
            target_genome_ids: target_genome_ids.clone(),
            transcript_ids: transcript_ids.clone(),
            ortholog_resource_path: Some(ortholog_resource_path.to_string()),
            upstream_bp,
            downstream_bp,
            ambiguity_policy,
            relationship,
        };
        Ok(OrthologPromoterCohortReport {
            schema: ORTHOLOG_PROMOTER_COHORT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            request,
            ortholog_resource_label: resource.label.or(resource.id),
            biological_contexts,
            resolved_promoter_count: rows.len(),
            unresolved_count: unresolved_rows.len(),
            rows,
            unresolved_rows,
            relationship,
            relationship_flags: vec![],
            warnings,
            ..OrthologPromoterCohortReport::default()
        })
    }

    fn ortholog_tfbs_summary_row(
        promoter: &OrthologPromoterRow,
        track: &TfbsScoreTrackRow,
    ) -> OrthologTfbsSummaryRow {
        let peak_position_0based = track.max_position_0based;
        let strand = promoter.strand.chars().next();
        OrthologTfbsSummaryRow {
            species: promoter.species.clone(),
            gene_label: promoter.display_label.clone(),
            transcript_id: promoter.transcript_id.clone(),
            tf_id: track.tf_id.clone(),
            tf_name: track.tf_name.clone(),
            max_score: track.max_score,
            peak_position_0based,
            peak_position_promoter_relative_bp: peak_position_0based
                .map(|value| value as i64 - promoter.tss_position_0based as i64),
            peak_genomic_position_1based: peak_position_0based.and_then(|value| {
                Self::promoter_local_position_to_genomic_1based(
                    strand,
                    promoter.promoter_start_1based,
                    promoter.promoter_end_1based,
                    promoter.promoter_length_bp,
                    value,
                )
            }),
            positive_fraction: Self::tfbs_track_positive_support_fraction(track),
        }
    }

    fn ortholog_promoter_gene_report(
        promoter: &OrthologPromoterRow,
        tfbs_score_tracks: TfbsScoreTrackReport,
    ) -> MultiGenePromoterTfbsGeneReport {
        MultiGenePromoterTfbsGeneReport {
            gene_query: promoter.gene_query.clone(),
            occurrence: 1,
            transcript_id_requested: promoter.transcript_id_requested.clone(),
            display_label: promoter.display_label.clone(),
            gene_id: promoter.gene_id.clone(),
            gene_name: promoter.gene_symbol.clone(),
            transcript_id: promoter.transcript_id.clone(),
            chromosome: promoter.chromosome.clone(),
            strand: promoter.strand.clone(),
            promoter_start_1based: promoter.promoter_start_1based,
            promoter_end_1based: promoter.promoter_end_1based,
            promoter_length_bp: promoter.promoter_length_bp,
            tss_1based: promoter.tss_1based,
            sequence_orientation: promoter.sequence_orientation.clone(),
            used_fuzzy_gene_match: false,
            tfbs_score_tracks,
        }
    }

    fn ortholog_sequence_similarity(
        left: &OrthologPromoterRow,
        right: &OrthologPromoterRow,
    ) -> OrthologSequenceSimilarityRow {
        let left_sequence = left.promoter_sequence.as_deref().unwrap_or("");
        let right_sequence = right.promoter_sequence.as_deref().unwrap_or("");
        let compared_length_bp = left_sequence.len().min(right_sequence.len());
        let identical_bp = left_sequence
            .as_bytes()
            .iter()
            .zip(right_sequence.as_bytes())
            .take(compared_length_bp)
            .filter(|(left, right)| left.eq_ignore_ascii_case(right))
            .count();
        let identity_fraction = if compared_length_bp == 0 {
            0.0
        } else {
            ((identical_bp as f64 / compared_length_bp as f64) * 1_000_000.0).round() / 1_000_000.0
        };
        OrthologSequenceSimilarityRow {
            left_species: left.species.clone(),
            right_species: right.species.clone(),
            left_gene_label: left.display_label.clone(),
            right_gene_label: right.display_label.clone(),
            alignment_mode: "direct_promoter_aligned_prefix".to_string(),
            compared_length_bp,
            identical_bp,
            identity_fraction,
        }
    }

    fn ortholog_peak_summary(
        peak: PromoterCohortTfbsPeakSummary,
        species_by_label: &BTreeMap<String, String>,
    ) -> OrthologTfbsPeakSummary {
        let mut species = peak
            .gene_labels
            .iter()
            .filter_map(|label| species_by_label.get(label).cloned())
            .collect::<Vec<_>>();
        species.sort();
        species.dedup();
        OrthologTfbsPeakSummary {
            tf_id: peak.tf_id,
            tf_name: peak.tf_name,
            promoter_count: peak.promoter_count,
            species,
            gene_labels: peak.gene_labels,
            max_score: peak.max_score,
            peak_positions_promoter_relative_bp: peak.peak_positions_promoter_relative_bp,
        }
    }

    fn ortholog_expression_assignments(
        cohort: &OrthologPromoterCohortReport,
        expression_rows: &[PromoterExpressionEvidenceInput],
        expression_source_label: Option<&str>,
        warnings: &mut Vec<String>,
    ) -> Vec<OrthologExpressionAssignment> {
        let source_label = expression_source_label
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("expression_input");
        let mut assignments = vec![];
        for row in expression_rows {
            let matched = cohort.rows.iter().find(|promoter| {
                row.transcript_id
                    .as_deref()
                    .map(|transcript_id| transcript_id == promoter.transcript_id)
                    .unwrap_or(false)
                    || row
                        .gene_label
                        .as_deref()
                        .map(|label| {
                            label == promoter.display_label
                                || promoter
                                    .gene_symbol
                                    .as_deref()
                                    .map(|symbol| label == symbol)
                                    .unwrap_or(false)
                        })
                        .unwrap_or(false)
            });
            if let Some(promoter) = matched {
                assignments.push(OrthologExpressionAssignment {
                    species: promoter.species.clone(),
                    gene_label: promoter.display_label.clone(),
                    sample_id: row.sample_id.clone(),
                    condition: row.condition.clone(),
                    value: row.value,
                    unit: row.unit.clone(),
                    source: row
                        .source
                        .clone()
                        .unwrap_or_else(|| source_label.to_string()),
                    assignment_note:
                        "Supplied expression metadata is assigned by gene label or transcript id; no causal inference is made."
                            .to_string(),
                });
            } else {
                warnings.push(format!(
                    "Expression row for gene_label='{}' transcript_id='{}' did not match an ortholog promoter row",
                    row.gene_label.as_deref().unwrap_or("-"),
                    row.transcript_id.as_deref().unwrap_or("-")
                ));
            }
        }
        assignments
    }

    fn ortholog_member_symbol(species: &str, gene_label: &str) -> String {
        format!("{species}: {gene_label}")
    }

    fn ortholog_member_key(species: &str, gene_label: &str) -> String {
        format!(
            "ortholog:{}:{}",
            Self::normalize_id_token(species),
            Self::normalize_id_token(gene_label)
        )
    }

    pub(crate) fn ortholog_tfbs_relationship_flags(
        relationship: GeneSetCohortRelationship,
        pairwise_similarity: &[OrthologPairwiseTfbsSimilarity],
    ) -> Vec<GeneSetCohortRelationshipFlag> {
        match relationship {
            GeneSetCohortRelationship::Unspecified | GeneSetCohortRelationship::Manual => vec![],
            GeneSetCohortRelationship::CoRegulated => pairwise_similarity
                .iter()
                .filter(|row| {
                    row.mean_smoothed_spearman
                        < Self::PROMOTER_COHORT_DIVERGENCE_SIMILARITY_THRESHOLD
                })
                .map(|row| GeneSetCohortRelationshipFlag {
                    flag_kind: "unexpected_divergence".to_string(),
                    evidence_kind: "tfbs_score_track_similarity".to_string(),
                    member_symbols: vec![
                        Self::ortholog_member_symbol(&row.left_species, &row.left_gene_label),
                        Self::ortholog_member_symbol(&row.right_species, &row.right_gene_label),
                    ],
                    member_dedup_keys: vec![
                        Self::ortholog_member_key(&row.left_species, &row.left_gene_label),
                        Self::ortholog_member_key(&row.right_species, &row.right_gene_label),
                    ],
                    detail: format!(
                        "Declared co-regulated ortholog promoter cohort, but '{}' ({}) and '{}' ({}) have low mean smoothed TFBS-track Spearman similarity ({:+.3}; threshold < {:+.2}). This is an evidence-triage flag, not a regulatory verdict.",
                        row.left_gene_label,
                        row.left_species,
                        row.right_gene_label,
                        row.right_species,
                        row.mean_smoothed_spearman,
                        Self::PROMOTER_COHORT_DIVERGENCE_SIMILARITY_THRESHOLD,
                    ),
                })
                .collect(),
            GeneSetCohortRelationship::AntiCoRegulated => pairwise_similarity
                .iter()
                .filter(|row| {
                    row.mean_smoothed_spearman
                        >= Self::PROMOTER_COHORT_CONCORDANCE_SIMILARITY_THRESHOLD
                })
                .map(|row| GeneSetCohortRelationshipFlag {
                    flag_kind: "unexpected_concordance".to_string(),
                    evidence_kind: "tfbs_score_track_similarity".to_string(),
                    member_symbols: vec![
                        Self::ortholog_member_symbol(&row.left_species, &row.left_gene_label),
                        Self::ortholog_member_symbol(&row.right_species, &row.right_gene_label),
                    ],
                    member_dedup_keys: vec![
                        Self::ortholog_member_key(&row.left_species, &row.left_gene_label),
                        Self::ortholog_member_key(&row.right_species, &row.right_gene_label),
                    ],
                    detail: format!(
                        "Declared anti-co-regulated ortholog promoter cohort, but '{}' ({}) and '{}' ({}) have high mean smoothed TFBS-track Spearman similarity ({:+.3}; threshold >= {:+.2}). This is an evidence-triage flag, not a regulatory verdict.",
                        row.left_gene_label,
                        row.left_species,
                        row.right_gene_label,
                        row.right_species,
                        row.mean_smoothed_spearman,
                        Self::PROMOTER_COHORT_CONCORDANCE_SIMILARITY_THRESHOLD,
                    ),
                })
                .collect(),
        }
    }

    fn ortholog_cutrun_relationship_evidence_class(
        row: &OrthologCutRunSupportRow,
    ) -> Option<&'static str> {
        match row.status {
            OrthologCutRunSupportStatus::NotComparable => None,
            OrthologCutRunSupportStatus::Confirmed | OrthologCutRunSupportStatus::Nearby => {
                Some("motif_supported")
            }
            OrthologCutRunSupportStatus::OccupancyOnly => Some("occupancy_only"),
            OrthologCutRunSupportStatus::MotifOnly => Some("motif_only"),
            OrthologCutRunSupportStatus::NoData => Some("no_support_detected"),
        }
    }

    pub(crate) fn ortholog_cutrun_relationship_flags(
        relationship: GeneSetCohortRelationship,
        cutrun_support: &[OrthologCutRunSupportRow],
    ) -> Vec<GeneSetCohortRelationshipFlag> {
        if matches!(
            relationship,
            GeneSetCohortRelationship::Unspecified | GeneSetCohortRelationship::Manual
        ) {
            return vec![];
        }

        let mut groups = Vec::<(&'static str, Vec<&OrthologCutRunSupportRow>)>::new();
        for row in cutrun_support {
            let Some(class) = Self::ortholog_cutrun_relationship_evidence_class(row) else {
                continue;
            };
            if let Some((_, rows)) = groups.iter_mut().find(|(candidate, _)| *candidate == class) {
                rows.push(row);
            } else {
                groups.push((class, vec![row]));
            }
        }
        if groups.is_empty() {
            return vec![];
        }

        fn relationship_flag(
            flag_kind: &str,
            evidence_class: &str,
            rows: &[&OrthologCutRunSupportRow],
            detail: String,
        ) -> GeneSetCohortRelationshipFlag {
            GeneSetCohortRelationshipFlag {
                flag_kind: flag_kind.to_string(),
                evidence_kind: format!("ortholog_cutrun_{evidence_class}"),
                member_symbols: rows
                    .iter()
                    .map(|row| GentleEngine::ortholog_member_symbol(&row.species, &row.gene_label))
                    .collect(),
                member_dedup_keys: rows
                    .iter()
                    .map(|row| GentleEngine::ortholog_member_key(&row.species, &row.gene_label))
                    .collect(),
                detail,
            }
        }

        match relationship {
            GeneSetCohortRelationship::CoRegulated => {
                if groups.len() < 2 {
                    return vec![];
                }
                let max_group_size = groups
                    .iter()
                    .map(|(_, rows)| rows.len())
                    .max()
                    .unwrap_or_default();
                let max_group_count = groups
                    .iter()
                    .filter(|(_, rows)| rows.len() == max_group_size)
                    .count();
                groups
                    .into_iter()
                    .filter(|(_, rows)| max_group_count > 1 || rows.len() < max_group_size)
                    .map(|(class, rows)| {
                        let symbols = rows
                            .iter()
                            .map(|row| Self::ortholog_member_symbol(&row.species, &row.gene_label))
                            .collect::<Vec<_>>()
                            .join(", ");
                        relationship_flag(
                            "unexpected_divergence",
                            class,
                            &rows,
                            format!(
                                "Declared co_regulated ortholog expectation found divergent CUT&RUN support class '{class}' for: {symbols}"
                            ),
                        )
                    })
                    .collect()
            }
            GeneSetCohortRelationship::AntiCoRegulated => groups
                .into_iter()
                .filter(|(_, rows)| rows.len() > 1)
                .map(|(class, rows)| {
                    let symbols = rows
                        .iter()
                        .map(|row| Self::ortholog_member_symbol(&row.species, &row.gene_label))
                        .collect::<Vec<_>>()
                        .join(", ");
                    relationship_flag(
                        "unexpected_concordance",
                        class,
                        &rows,
                        format!(
                            "Declared anti_co_regulated ortholog expectation found concordant CUT&RUN support class '{class}' for: {symbols}"
                        ),
                    )
                })
                .collect(),
            GeneSetCohortRelationship::Unspecified | GeneSetCohortRelationship::Manual => vec![],
        }
    }

    pub(crate) fn summarize_ortholog_promoter_comparison(
        &self,
        mut cohort: OrthologPromoterCohortReport,
        motifs: &[String],
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        relationship: GeneSetCohortRelationship,
        expression_rows: &[PromoterExpressionEvidenceInput],
        expression_source_label: Option<&str>,
        cutrun_dataset_ids: &[String],
        cutrun_read_report_ids: &[String],
    ) -> Result<OrthologPromoterComparisonReport, EngineError> {
        if cohort.rows.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeOrthologPromoterComparison requires at least one resolved promoter row".to_string(),
                cause_chain: vec![],
            });
        }
        if motifs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "SummarizeOrthologPromoterComparison requires at least one TF motif query"
                    .to_string(),
                cause_chain: vec![],
            });
        }

        let effective_relationship = if relationship == GeneSetCohortRelationship::Unspecified {
            cohort.relationship
        } else {
            relationship
        };
        cohort.relationship = effective_relationship;
        cohort.request.relationship = effective_relationship;

        let mut warnings = cohort.warnings.clone();
        let mut gene_reports = vec![];
        let mut promoter_summaries = vec![];
        for promoter in &cohort.rows {
            let Some(sequence) = promoter.promoter_sequence.clone() else {
                warnings.push(format!(
                    "No promoter sequence was available for '{}' in '{}'; skipping TFBS scoring",
                    promoter.display_label, promoter.species
                ));
                continue;
            };
            let mut report = self.summarize_tfbs_score_tracks(
                SequenceScanTarget::InlineSequence {
                    sequence_text: sequence,
                    topology: InlineSequenceTopology::Linear,
                    id_hint: Some(format!("{} {}", promoter.species, promoter.display_label)),
                    span_start_0based: None,
                    span_end_0based_exclusive: None,
                },
                motifs,
                score_kind,
                clip_negative,
            )?;
            report.tss_markers = vec![TfbsScoreTrackTssMarker {
                feature_id: usize::MAX,
                feature_kind: "ortholog_promoter_slice".to_string(),
                label: promoter.transcript_id.clone(),
                position_0based: promoter.tss_position_0based,
                is_reverse: false,
            }];
            for track in &report.tracks {
                promoter_summaries.push(Self::ortholog_tfbs_summary_row(promoter, track));
            }
            gene_reports.push(Self::ortholog_promoter_gene_report(promoter, report));
        }

        let mut pairwise_tfbs_similarity = vec![];
        for left_idx in 0..gene_reports.len() {
            for right_idx in left_idx.saturating_add(1)..gene_reports.len() {
                if let Some(row) = Self::summarize_promoter_cohort_pairwise_similarity(
                    &gene_reports[left_idx],
                    &gene_reports[right_idx],
                ) {
                    let left_species = cohort
                        .rows
                        .iter()
                        .find(|promoter| promoter.display_label == row.left_gene_label)
                        .map(|promoter| promoter.species.clone())
                        .unwrap_or_default();
                    let right_species = cohort
                        .rows
                        .iter()
                        .find(|promoter| promoter.display_label == row.right_gene_label)
                        .map(|promoter| promoter.species.clone())
                        .unwrap_or_default();
                    pairwise_tfbs_similarity.push(OrthologPairwiseTfbsSimilarity {
                        left_species,
                        right_species,
                        left_gene_label: row.left_gene_label,
                        right_gene_label: row.right_gene_label,
                        shared_motif_count: row.shared_motif_count,
                        mean_raw_pearson: row.mean_raw_pearson,
                        mean_smoothed_spearman: row.mean_smoothed_spearman,
                        motif_ids: row.motif_ids,
                    });
                }
            }
        }
        pairwise_tfbs_similarity.sort_by(|left, right| {
            right
                .mean_smoothed_spearman
                .total_cmp(&left.mean_smoothed_spearman)
                .then(left.left_species.cmp(&right.left_species))
                .then(left.right_species.cmp(&right.right_species))
        });

        let summary_rows_for_peaks = promoter_summaries
            .iter()
            .map(|row| MultiGenePromoterTfbsSummaryRow {
                gene_label: row.gene_label.clone(),
                gene_query: row.gene_label.clone(),
                transcript_id: row.transcript_id.clone(),
                tf_id: row.tf_id.clone(),
                tf_name: row.tf_name.clone(),
                max_score: row.max_score,
                peak_position_0based: row.peak_position_0based,
                peak_position_promoter_relative_bp: row.peak_position_promoter_relative_bp,
                peak_genomic_position_1based: row.peak_genomic_position_1based,
                positive_fraction: row.positive_fraction,
            })
            .collect::<Vec<_>>();
        let species_by_label = cohort
            .rows
            .iter()
            .map(|row| (row.display_label.clone(), row.species.clone()))
            .collect::<BTreeMap<_, _>>();
        let (shared, specific) =
            Self::summarize_promoter_cohort_peak_sets(&summary_rows_for_peaks, gene_reports.len());
        let conserved_tfbs_peaks = shared
            .into_iter()
            .map(|peak| Self::ortholog_peak_summary(peak, &species_by_label))
            .collect();
        let species_specific_tfbs_peaks = specific
            .into_iter()
            .map(|peak| Self::ortholog_peak_summary(peak, &species_by_label))
            .collect();

        let mut sequence_similarity = vec![];
        for left_idx in 0..cohort.rows.len() {
            for right_idx in left_idx.saturating_add(1)..cohort.rows.len() {
                sequence_similarity.push(Self::ortholog_sequence_similarity(
                    &cohort.rows[left_idx],
                    &cohort.rows[right_idx],
                ));
            }
        }

        let expression_assignments = Self::ortholog_expression_assignments(
            &cohort,
            expression_rows,
            expression_source_label,
            &mut warnings,
        );
        let cutrun_support = self.inspect_ortholog_cutrun_support(
            &cohort,
            &promoter_summaries,
            cutrun_dataset_ids,
            cutrun_read_report_ids,
            &mut warnings,
        )?;
        let mut relationship_flags = Self::ortholog_tfbs_relationship_flags(
            effective_relationship,
            &pairwise_tfbs_similarity,
        );
        relationship_flags.extend(Self::ortholog_cutrun_relationship_flags(
            effective_relationship,
            &cutrun_support,
        ));
        warnings.extend(relationship_flags.iter().map(|flag| flag.detail.clone()));
        cohort.relationship_flags = relationship_flags.clone();

        Ok(OrthologPromoterComparisonReport {
            schema: ORTHOLOG_PROMOTER_COMPARISON_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            cohort,
            motifs_requested: motifs.to_vec(),
            score_kind: score_kind.as_str().to_string(),
            clip_negative,
            promoter_summaries,
            pairwise_tfbs_similarity,
            conserved_tfbs_peaks,
            species_specific_tfbs_peaks,
            sequence_similarity,
            cutrun_support,
            expression_assignments,
            relationship: effective_relationship,
            relationship_flags,
            warnings,
            ..OrthologPromoterComparisonReport::default()
        })
    }
}
