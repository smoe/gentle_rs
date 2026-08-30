//! Ortholog promoter cohort resolution and comparison.
//!
//! This slice is intentionally offline-first: local ortholog resources map
//! genes across species, prepared genome catalogs resolve promoter windows, and
//! evidence summaries remain separated by evidence type.

use super::*;
use gentle_protocol::{
    BiologicalContext, BiologicalContextRegistry, GeneSetCohortRelationship,
    GeneSetCohortRelationshipFlag, ORTHOLOG_PROMOTER_COHORT_SCHEMA,
    ORTHOLOG_PROMOTER_COMPARISON_SCHEMA, ORTHOLOG_PROMOTER_CONSERVATION_SCHEMA,
    ORTHOLOG_RESOURCE_SCHEMA, OrthologAmbiguityCandidate, OrthologConfidence,
    OrthologConservationPairwiseAlignment, OrthologConservedInterval,
    OrthologCutRunNormalizationInput, OrthologCutRunNormalizedAssignment,
    OrthologCutRunNormalizedValueInput, OrthologCutRunPairwiseQuantitativeComparison,
    OrthologCutRunQuantitativeComparison, OrthologCutRunQuantitativeComparisonStatus,
    OrthologCutRunSupportRow, OrthologCutRunSupportStatus, OrthologExpressionAssignment,
    OrthologMappingRow, OrthologPairwiseTfbsSimilarity, OrthologPromoterCohortReport,
    OrthologPromoterCohortRequest, OrthologPromoterComparisonReport,
    OrthologPromoterConservationReport, OrthologPromoterRole, OrthologPromoterRow,
    OrthologResource, OrthologSequenceSimilarityRow, OrthologTfbsPeakSummary,
    OrthologTfbsSummaryRow, OrthologUnresolvedRow, OrthologyType,
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

    fn ortholog_runtime_context_id(species: &str, genome_id: Option<&str>) -> String {
        let genome_key = genome_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(Self::ortholog_gene_key)
            .unwrap_or_else(|| "unspecified_genome".to_string());
        format!(
            "ortholog_{}_{}",
            Self::ortholog_species_key(species),
            genome_key
        )
    }

    fn register_ortholog_report_context(
        report_contexts: &mut BiologicalContextRegistry,
        resource: &OrthologResource,
        explicit_context_id: Option<&str>,
        species: &str,
        genome_id: Option<&str>,
    ) -> Result<String, EngineError> {
        if let Some(genome_id) = genome_id {
            Self::validate_ortholog_context_genome(
                resource,
                explicit_context_id,
                genome_id,
                "mapping",
            )?;
        }
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
            context.genome_id = genome_id.map(str::to_string);
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

    fn ortholog_candidate_sort_key(candidate: &OrientedOrthologMapping) -> String {
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

    fn ortholog_candidate_label(candidate: &OrientedOrthologMapping) -> String {
        let base = Self::ortholog_candidate_sort_key(candidate);
        candidate
            .source
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|source| format!("{base}:source={source}"))
            .unwrap_or(base)
    }

    fn ortholog_candidate_target_genome_id(
        resource: &OrthologResource,
        candidate: &OrientedOrthologMapping,
        requested_target_genome_id: Option<&str>,
    ) -> Option<String> {
        if let Some(genome_id) = requested_target_genome_id {
            return Some(genome_id.to_string());
        }
        candidate
            .target_context_id
            .as_deref()
            .and_then(|context_id| resource.biological_contexts.context(context_id).ok())
            .and_then(|context| context.genome_id.clone())
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
            Some(trimmed_anchor_genome_id),
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
                Self::ortholog_candidate_sort_key(left)
                    .cmp(&Self::ortholog_candidate_sort_key(right))
            });
            if candidates.len() > 1 && ambiguity_policy == OrthologAmbiguityPolicy::Preserve {
                let labels = candidates
                    .iter()
                    .map(Self::ortholog_candidate_label)
                    .collect::<Vec<_>>();
                let mut candidate_mappings = Vec::with_capacity(candidates.len());
                let mut candidate_context_ids = BTreeSet::new();
                let mut candidate_genome_ids = BTreeSet::new();
                let mut missing_target_genome_id = false;
                for (candidate_index, candidate) in candidates.iter().enumerate() {
                    let target_genome_id = Self::ortholog_candidate_target_genome_id(
                        &resource,
                        candidate,
                        requested_target_genome_id.as_deref(),
                    );
                    missing_target_genome_id |= target_genome_id.is_none();
                    let source_context_id = Self::register_ortholog_report_context(
                        &mut biological_contexts,
                        &resource,
                        candidate.source_context_id.as_deref(),
                        &anchor_species,
                        Some(trimmed_anchor_genome_id),
                    )?;
                    let target_context_id = Self::register_ortholog_report_context(
                        &mut biological_contexts,
                        &resource,
                        candidate.target_context_id.as_deref(),
                        &candidate.target_species,
                        target_genome_id.as_deref(),
                    )?;
                    candidate_context_ids.insert(target_context_id.clone());
                    if let Some(target_genome_id) = target_genome_id.as_ref() {
                        candidate_genome_ids.insert(target_genome_id.clone());
                    }
                    candidate_mappings.push(OrthologAmbiguityCandidate {
                        candidate_rank: candidate_index + 1,
                        candidate_label: labels[candidate_index].clone(),
                        target_species: candidate.target_species.clone(),
                        target_genome_id,
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
                if missing_target_genome_id {
                    warnings.push(format!(
                        "No target genome id or context genome was provided for ambiguous species '{}'; preserved candidate mappings do not claim a genome identity",
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
            let declared_target_genome_id = Self::ortholog_candidate_target_genome_id(
                &resource,
                &candidate,
                requested_target_genome_id.as_deref(),
            );
            let used_species_as_genome_id = declared_target_genome_id.is_none();
            let target_genome_id = declared_target_genome_id
                .clone()
                .unwrap_or_else(|| candidate.target_species.clone());
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
                Some(trimmed_anchor_genome_id),
            )?;
            let target_context_id = Self::register_ortholog_report_context(
                &mut biological_contexts,
                &resource,
                candidate.target_context_id.as_deref(),
                &candidate.target_species,
                declared_target_genome_id.as_deref(),
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
    ) -> Result<OrthologSequenceSimilarityRow, EngineError> {
        let left_sequence = left.promoter_sequence.as_deref().unwrap_or("");
        let right_sequence = right.promoter_sequence.as_deref().unwrap_or("");
        let alignment = Self::compute_pairwise_alignment_report(
            &format!("{} {}", left.species, left.display_label),
            left_sequence,
            None,
            None,
            &format!("{} {}", right.species, right.display_label),
            right_sequence,
            None,
            None,
            PairwiseAlignmentMode::Global,
            2,
            -3,
            -5,
            -1,
        )?;
        Ok(OrthologSequenceSimilarityRow {
            left_species: left.species.clone(),
            right_species: right.species.clone(),
            left_gene_label: left.display_label.clone(),
            right_gene_label: right.display_label.clone(),
            alignment_mode: "global_pairwise".to_string(),
            compared_length_bp: alignment.report.aligned_columns,
            identical_bp: alignment.report.matches,
            identity_fraction: alignment.report.identity_fraction,
        })
    }

    pub(crate) fn compare_ortholog_promoter_conservation(
        &self,
        cohort: OrthologPromoterCohortReport,
        min_conserved_bp: usize,
    ) -> Result<OrthologPromoterConservationReport, EngineError> {
        if min_conserved_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Ortholog promoter conservation requires min_conserved_bp >= 1"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        let anchor = cohort
            .rows
            .iter()
            .find(|row| row.role == OrthologPromoterRole::Anchor)
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: "Ortholog promoter conservation requires one resolved anchor row"
                    .to_string(),
                cause_chain: vec![],
            })?;
        let anchor_sequence = anchor
            .promoter_sequence
            .as_deref()
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: "Ortholog promoter conservation anchor has no sequence".to_string(),
                cause_chain: vec![],
            })?;
        let targets = cohort
            .rows
            .iter()
            .filter(|row| row.role == OrthologPromoterRole::Target)
            .collect::<Vec<_>>();
        if targets.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Ortholog promoter conservation requires at least one resolved target row"
                    .to_string(),
                cause_chain: vec![],
            });
        }

        let mut common_mask = vec![true; anchor_sequence.len()];
        let mut pairwise_alignments = Vec::with_capacity(targets.len());
        let mut supporting_species = Vec::with_capacity(targets.len() + 1);
        supporting_species.push(anchor.species.clone());
        for target in targets {
            let target_sequence =
                target
                    .promoter_sequence
                    .as_deref()
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Ortholog promoter conservation target '{}' has no sequence",
                            target.species
                        ),
                        cause_chain: vec![],
                    })?;
            let computed = Self::compute_pairwise_alignment_report(
                &format!("{} {}", anchor.species, anchor.display_label),
                anchor_sequence,
                None,
                None,
                &format!("{} {}", target.species, target.display_label),
                target_sequence,
                None,
                None,
                PairwiseAlignmentMode::Global,
                2,
                -3,
                -5,
                -1,
            )?;
            let mut target_mask = vec![false; anchor_sequence.len()];
            let mut anchor_index = 0usize;
            let mut target_index = 0usize;
            for operation in &computed.operations {
                match operation {
                    bio::alignment::AlignmentOperation::Match => {
                        if anchor_index < target_mask.len()
                            && target_index < target_sequence.len()
                            && anchor_sequence.as_bytes()[anchor_index]
                                .eq_ignore_ascii_case(&target_sequence.as_bytes()[target_index])
                        {
                            target_mask[anchor_index] = true;
                        }
                        anchor_index += 1;
                        target_index += 1;
                    }
                    bio::alignment::AlignmentOperation::Subst => {
                        anchor_index += 1;
                        target_index += 1;
                    }
                    bio::alignment::AlignmentOperation::Ins => anchor_index += 1,
                    bio::alignment::AlignmentOperation::Del => target_index += 1,
                    bio::alignment::AlignmentOperation::Xclip(length) => anchor_index += *length,
                    bio::alignment::AlignmentOperation::Yclip(length) => target_index += *length,
                }
            }
            for (common, target_match) in common_mask.iter_mut().zip(target_mask.iter()) {
                *common &= *target_match;
            }
            let identical_anchor_bp = target_mask.iter().filter(|matched| **matched).count();
            pairwise_alignments.push(OrthologConservationPairwiseAlignment {
                target_species: target.species.clone(),
                target_gene_label: target.display_label.clone(),
                target_transcript_id: target.transcript_id.clone(),
                anchor_coverage_bp: computed.report.matches + computed.report.mismatches,
                identical_anchor_bp,
                alignment: computed.report,
            });
            supporting_species.push(target.species.clone());
        }

        let mut conserved_intervals = vec![];
        let mut interval_start = None;
        for (position, conserved) in common_mask
            .iter()
            .copied()
            .chain(std::iter::once(false))
            .enumerate()
        {
            if conserved {
                interval_start.get_or_insert(position);
                continue;
            }
            let Some(start) = interval_start.take() else {
                continue;
            };
            if position.saturating_sub(start) < min_conserved_bp {
                continue;
            }
            let strand = anchor.strand.chars().next();
            let genomic_a = Self::promoter_local_position_to_genomic_1based(
                strand,
                anchor.promoter_start_1based,
                anchor.promoter_end_1based,
                anchor.promoter_length_bp,
                start,
            );
            let genomic_b = Self::promoter_local_position_to_genomic_1based(
                strand,
                anchor.promoter_start_1based,
                anchor.promoter_end_1based,
                anchor.promoter_length_bp,
                position.saturating_sub(1),
            );
            conserved_intervals.push(OrthologConservedInterval {
                anchor_start_0based: start,
                anchor_end_0based_exclusive: position,
                length_bp: position - start,
                promoter_relative_start_bp: start as i64 - anchor.tss_position_0based as i64,
                promoter_relative_end_bp_exclusive: position as i64
                    - anchor.tss_position_0based as i64,
                genomic_start_1based: genomic_a.zip(genomic_b).map(|(a, b)| a.min(b)),
                genomic_end_1based: genomic_a.zip(genomic_b).map(|(a, b)| a.max(b)),
                supporting_species: supporting_species.clone(),
            });
        }
        Ok(OrthologPromoterConservationReport {
            schema: ORTHOLOG_PROMOTER_CONSERVATION_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            anchor_species: anchor.species.clone(),
            anchor_gene_label: anchor.display_label.clone(),
            anchor_transcript_id: anchor.transcript_id.clone(),
            min_conserved_bp,
            pairwise_alignments,
            conserved_intervals,
            warnings: cohort.warnings.clone(),
            cohort,
            ..OrthologPromoterConservationReport::default()
        })
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

    fn ortholog_cutrun_normalized_value_matches(
        promoter: &OrthologPromoterRow,
        value: &OrthologCutRunNormalizedValueInput,
    ) -> bool {
        if Self::ortholog_species_key(&promoter.species)
            != Self::ortholog_species_key(&value.species)
        {
            return false;
        }
        if let Some(gene_label) = value
            .gene_label
            .as_deref()
            .map(str::trim)
            .filter(|label| !label.is_empty())
        {
            let requested = Self::ortholog_gene_key(gene_label);
            let gene_matches = [
                Some(promoter.display_label.as_str()),
                promoter.gene_symbol.as_deref(),
                Some(promoter.gene_query.as_str()),
            ]
            .into_iter()
            .flatten()
            .any(|candidate| Self::ortholog_gene_key(candidate) == requested);
            if !gene_matches {
                return false;
            }
        }
        if let Some(transcript_id) = value
            .transcript_id
            .as_deref()
            .map(str::trim)
            .filter(|transcript_id| !transcript_id.is_empty())
            && transcript_id != promoter.transcript_id
        {
            return false;
        }
        true
    }

    fn summarize_ortholog_cutrun_quantitative_comparison(
        cohort: &OrthologPromoterCohortReport,
        cutrun_support: &[OrthologCutRunSupportRow],
        cutrun_dataset_ids: &[String],
        cutrun_read_report_ids: &[String],
        normalization: Option<&OrthologCutRunNormalizationInput>,
        warnings: &mut Vec<String>,
    ) -> Option<OrthologCutRunQuantitativeComparison> {
        let has_selected_sources =
            !cutrun_dataset_ids.is_empty() || !cutrun_read_report_ids.is_empty();
        let Some(normalization) = normalization else {
            if !has_selected_sources {
                return None;
            }
            let detail = "Selected CUT&RUN evidence was interpreted qualitatively only because no explicit normalization metadata and normalized values were supplied. Raw peak scores, signal heights, fragment counts, and read depth were not compared across species.".to_string();
            warnings.push(detail.clone());
            return Some(OrthologCutRunQuantitativeComparison {
                status: OrthologCutRunQuantitativeComparisonStatus::NotComparable,
                detail,
                ..OrthologCutRunQuantitativeComparison::default()
            });
        };

        let selected_dataset_ids = cutrun_dataset_ids
            .iter()
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>();
        let selected_read_report_ids = cutrun_read_report_ids
            .iter()
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>();
        let mut reasons = Vec::<String>::new();
        for (label, value) in [
            (
                "normalization_method",
                normalization.normalization_method.as_str(),
            ),
            ("unit", normalization.unit.as_str()),
            (
                "comparison_reference",
                normalization.comparison_reference.as_str(),
            ),
            ("provenance", normalization.provenance.as_str()),
        ] {
            if value.trim().is_empty() {
                reasons.push(format!("normalization metadata field '{label}' is empty"));
            }
        }
        if normalization.values.is_empty() {
            reasons.push("no normalized promoter values were supplied".to_string());
        }
        if cohort.rows.len() < 2 {
            reasons.push(
                "at least two resolved promoter rows are required for quantitative comparison"
                    .to_string(),
            );
        }

        for (index, value) in normalization.values.iter().enumerate() {
            if value.species.trim().is_empty() {
                reasons.push(format!(
                    "normalized value {} has an empty species",
                    index + 1
                ));
            }
            if !value.normalized_value.is_finite() {
                reasons.push(format!(
                    "normalized value {} for '{}' is not finite",
                    index + 1,
                    value.species
                ));
            }
            if value.contributing_dataset_ids.is_empty()
                && value.contributing_read_report_ids.is_empty()
            {
                reasons.push(format!(
                    "normalized value {} for '{}' has no contributing dataset or read-report id",
                    index + 1,
                    value.species
                ));
            }
            for dataset_id in &value.contributing_dataset_ids {
                if !selected_dataset_ids.contains(dataset_id.trim()) {
                    reasons.push(format!(
                        "normalized value {} references CUT&RUN dataset '{}' that was not selected",
                        index + 1,
                        dataset_id
                    ));
                }
            }
            for report_id in &value.contributing_read_report_ids {
                if !selected_read_report_ids.contains(report_id.trim()) {
                    reasons.push(format!(
                        "normalized value {} references CUT&RUN read report '{}' that was not selected",
                        index + 1,
                        report_id
                    ));
                }
            }
            if !cohort
                .rows
                .iter()
                .any(|promoter| Self::ortholog_cutrun_normalized_value_matches(promoter, value))
            {
                reasons.push(format!(
                    "normalized value {} for species '{}' did not match a resolved promoter row",
                    index + 1,
                    value.species
                ));
            }
        }

        let mut assignments = Vec::<OrthologCutRunNormalizedAssignment>::new();
        for promoter in &cohort.rows {
            let matches = normalization
                .values
                .iter()
                .filter(|value| Self::ortholog_cutrun_normalized_value_matches(promoter, value))
                .collect::<Vec<_>>();
            if matches.is_empty() {
                reasons.push(format!(
                    "no normalized value matched '{}' ({})",
                    promoter.display_label, promoter.species
                ));
                continue;
            }
            if matches.len() > 1 {
                reasons.push(format!(
                    "{} normalized values matched '{}' ({}); one unambiguous value is required",
                    matches.len(),
                    promoter.display_label,
                    promoter.species
                ));
                continue;
            }
            let value = matches[0];
            let support = cutrun_support.iter().find(|support| {
                Self::ortholog_species_key(&support.species)
                    == Self::ortholog_species_key(&promoter.species)
                    && Self::ortholog_gene_key(&support.gene_label)
                        == Self::ortholog_gene_key(&promoter.display_label)
            });
            let Some(support) = support else {
                reasons.push(format!(
                    "no qualitative CUT&RUN support row was available for '{}' ({})",
                    promoter.display_label, promoter.species
                ));
                continue;
            };
            if support.status == OrthologCutRunSupportStatus::NotComparable {
                reasons.push(format!(
                    "CUT&RUN provenance did not match '{}' ({})",
                    promoter.display_label, promoter.species
                ));
            }
            for dataset_id in &value.contributing_dataset_ids {
                if !support
                    .contributing_dataset_ids
                    .iter()
                    .any(|candidate| candidate.trim() == dataset_id.trim())
                {
                    reasons.push(format!(
                        "dataset '{}' was not evaluated for '{}' ({})",
                        dataset_id, promoter.display_label, promoter.species
                    ));
                }
            }
            for report_id in &value.contributing_read_report_ids {
                if !support
                    .contributing_read_report_ids
                    .iter()
                    .any(|candidate| candidate.trim() == report_id.trim())
                {
                    reasons.push(format!(
                        "read report '{}' was not evaluated for '{}' ({})",
                        report_id, promoter.display_label, promoter.species
                    ));
                }
            }
            assignments.push(OrthologCutRunNormalizedAssignment {
                species: promoter.species.clone(),
                gene_label: promoter.display_label.clone(),
                transcript_id: promoter.transcript_id.clone(),
                normalized_value: value.normalized_value,
                contributing_dataset_ids: value.contributing_dataset_ids.clone(),
                contributing_read_report_ids: value.contributing_read_report_ids.clone(),
                provenance: value.provenance.clone(),
            });
        }

        reasons.sort();
        reasons.dedup();
        if !reasons.is_empty() {
            let detail = format!(
                "Quantitative CUT&RUN comparison was not performed: {}. Qualitative support states remain available; raw cross-species intensity was not compared.",
                reasons.join("; ")
            );
            warnings.push(detail.clone());
            return Some(OrthologCutRunQuantitativeComparison {
                status: OrthologCutRunQuantitativeComparisonStatus::NotComparable,
                normalization: normalization.clone(),
                assignments,
                pairwise_comparisons: vec![],
                detail,
            });
        }

        let mut pairwise_comparisons = Vec::new();
        for left_index in 0..assignments.len() {
            for right_index in left_index.saturating_add(1)..assignments.len() {
                let left = &assignments[left_index];
                let right = &assignments[right_index];
                let delta = right.normalized_value - left.normalized_value;
                pairwise_comparisons.push(OrthologCutRunPairwiseQuantitativeComparison {
                    left_species: left.species.clone(),
                    right_species: right.species.clone(),
                    left_gene_label: left.gene_label.clone(),
                    right_gene_label: right.gene_label.clone(),
                    left_normalized_value: left.normalized_value,
                    right_normalized_value: right.normalized_value,
                    delta_right_minus_left: delta,
                    absolute_delta: delta.abs(),
                });
            }
        }
        Some(OrthologCutRunQuantitativeComparison {
            status: OrthologCutRunQuantitativeComparisonStatus::Comparable,
            normalization: normalization.clone(),
            assignments,
            pairwise_comparisons,
            detail: format!(
                "Compared explicitly normalized CUT&RUN values using '{}' in '{}' against shared reference '{}'. Raw source intensities were not compared directly.",
                normalization.normalization_method,
                normalization.unit,
                normalization.comparison_reference
            ),
        })
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
        cutrun_normalization: Option<&OrthologCutRunNormalizationInput>,
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
                )?);
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
        let cutrun_quantitative_comparison =
            Self::summarize_ortholog_cutrun_quantitative_comparison(
                &cohort,
                &cutrun_support,
                cutrun_dataset_ids,
                cutrun_read_report_ids,
                cutrun_normalization,
                &mut warnings,
            );
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
            cutrun_quantitative_comparison,
            expression_assignments,
            relationship: effective_relationship,
            relationship_flags,
            warnings,
            ..OrthologPromoterComparisonReport::default()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quantitative_test_cohort() -> OrthologPromoterCohortReport {
        OrthologPromoterCohortReport {
            resolved_promoter_count: 2,
            rows: vec![
                OrthologPromoterRow {
                    species: "Homo sapiens".to_string(),
                    display_label: "TP73".to_string(),
                    gene_query: "TP73".to_string(),
                    transcript_id: "ENST_TP73".to_string(),
                    ..OrthologPromoterRow::default()
                },
                OrthologPromoterRow {
                    species: "Mus musculus".to_string(),
                    display_label: "Trp73".to_string(),
                    gene_query: "Trp73".to_string(),
                    transcript_id: "ENSMUST_Trp73".to_string(),
                    ..OrthologPromoterRow::default()
                },
            ],
            ..OrthologPromoterCohortReport::default()
        }
    }

    #[test]
    fn normalized_cutrun_comparison_requires_complete_provenance_and_matched_sources() {
        let cohort = quantitative_test_cohort();
        let support = vec![
            OrthologCutRunSupportRow {
                species: "Homo sapiens".to_string(),
                gene_label: "TP73".to_string(),
                status: OrthologCutRunSupportStatus::Confirmed,
                contributing_dataset_ids: vec!["human_cutrun".to_string()],
                ..OrthologCutRunSupportRow::default()
            },
            OrthologCutRunSupportRow {
                species: "Mus musculus".to_string(),
                gene_label: "Trp73".to_string(),
                status: OrthologCutRunSupportStatus::Confirmed,
                contributing_dataset_ids: vec!["mouse_cutrun".to_string()],
                ..OrthologCutRunSupportRow::default()
            },
        ];
        let normalization = OrthologCutRunNormalizationInput {
            normalization_method: "spike_in_scaled_cpm".to_string(),
            unit: "normalized_fragments_per_million".to_string(),
            comparison_reference: "shared_batch_1".to_string(),
            provenance: "synthetic reviewed normalization table".to_string(),
            values: vec![
                OrthologCutRunNormalizedValueInput {
                    species: "Homo sapiens".to_string(),
                    gene_label: Some("TP73".to_string()),
                    normalized_value: 4.0,
                    contributing_dataset_ids: vec!["human_cutrun".to_string()],
                    ..OrthologCutRunNormalizedValueInput::default()
                },
                OrthologCutRunNormalizedValueInput {
                    species: "Mus musculus".to_string(),
                    transcript_id: Some("ENSMUST_Trp73".to_string()),
                    normalized_value: 7.5,
                    contributing_dataset_ids: vec!["mouse_cutrun".to_string()],
                    ..OrthologCutRunNormalizedValueInput::default()
                },
            ],
        };
        let mut warnings = vec![];
        let comparison = GentleEngine::summarize_ortholog_cutrun_quantitative_comparison(
            &cohort,
            &support,
            &["human_cutrun".to_string(), "mouse_cutrun".to_string()],
            &[],
            Some(&normalization),
            &mut warnings,
        )
        .expect("quantitative comparison");
        assert_eq!(
            comparison.status,
            OrthologCutRunQuantitativeComparisonStatus::Comparable
        );
        assert_eq!(comparison.assignments.len(), 2);
        assert_eq!(comparison.pairwise_comparisons.len(), 1);
        assert_eq!(
            comparison.pairwise_comparisons[0].delta_right_minus_left,
            3.5
        );
        assert!(warnings.is_empty());

        let mut incomplete = normalization;
        incomplete.comparison_reference.clear();
        let rejected = GentleEngine::summarize_ortholog_cutrun_quantitative_comparison(
            &cohort,
            &support,
            &["human_cutrun".to_string(), "mouse_cutrun".to_string()],
            &[],
            Some(&incomplete),
            &mut warnings,
        )
        .expect("fail-closed quantitative comparison");
        assert_eq!(
            rejected.status,
            OrthologCutRunQuantitativeComparisonStatus::NotComparable
        );
        assert!(rejected.pairwise_comparisons.is_empty());
        assert!(rejected.detail.contains("comparison_reference"));
    }

    #[test]
    fn selected_cutrun_without_normalization_stays_qualitative() {
        let comparison = GentleEngine::summarize_ortholog_cutrun_quantitative_comparison(
            &quantitative_test_cohort(),
            &[],
            &["human_cutrun".to_string()],
            &[],
            None,
            &mut vec![],
        )
        .expect("qualitative-only status");
        assert_eq!(
            comparison.status,
            OrthologCutRunQuantitativeComparisonStatus::NotComparable
        );
        assert!(comparison.detail.contains("qualitatively only"));
    }
}
