//! Gene-set resolution and promoter-cohort construction.
//!
//! This module turns curated gene-group records, explicit user members,
//! local external mappings, prepared-genome neighborhoods, and deterministic
//! random samples into portable analysis operands.

use super::*;
use crate::gene_groups::LoadedGeneGroupRecord;
use gentle_protocol::{
    GENE_SET_CO_REGULATED_CACHE_SCHEMA, GENE_SET_DIRECT_LIST_CACHE_SCHEMA,
    GENE_SET_ONTOLOGY_ASSIGNMENT_CACHE_SCHEMA, GeneGroupMember,
    GeneSetCoRegulatedProducerMetadata, GeneSetProducerFilter, GeneSetProducerKind,
    GeneSetProducerProvenance, GeneSetProducerQueryMetadata,
};
use serde_json::Value;
use std::{collections::BTreeSet, path::Path};

#[derive(Clone, Debug)]
struct IndexedGeneRow {
    record: GenomeGeneRecord,
    dedup_key: String,
    symbol: String,
}

#[derive(Clone, Debug, Default)]
struct GeneSetDirectListCacheMetadata {
    provider_id: Option<String>,
    provider_label: Option<String>,
    provider_version: Option<String>,
    cache_id: Option<String>,
    cache_version: Option<String>,
    cache_digest: Option<String>,
    organism: Option<String>,
    taxon_id: Option<String>,
    symbol_namespace: Option<String>,
    review_status: Option<GeneSetResolutionReviewStatus>,
    filters: Vec<GeneSetProducerFilter>,
}

impl GeneSetDirectListCacheMetadata {
    fn merge_missing_from(&mut self, other: &Self) {
        if self.provider_id.is_none() {
            self.provider_id = other.provider_id.clone();
        }
        if self.provider_label.is_none() {
            self.provider_label = other.provider_label.clone();
        }
        if self.provider_version.is_none() {
            self.provider_version = other.provider_version.clone();
        }
        if self.cache_id.is_none() {
            self.cache_id = other.cache_id.clone();
        }
        if self.cache_version.is_none() {
            self.cache_version = other.cache_version.clone();
        }
        if self.cache_digest.is_none() {
            self.cache_digest = other.cache_digest.clone();
        }
        if self.organism.is_none() {
            self.organism = other.organism.clone();
        }
        if self.taxon_id.is_none() {
            self.taxon_id = other.taxon_id.clone();
        }
        if self.symbol_namespace.is_none() {
            self.symbol_namespace = other.symbol_namespace.clone();
        }
        if self.review_status.is_none() {
            self.review_status = other.review_status;
        }
        for filter in &other.filters {
            if !self.filters.contains(filter) {
                self.filters.push(filter.clone());
            }
        }
    }

    fn apply_overrides(
        &mut self,
        provider_id: Option<&str>,
        provider_label: Option<&str>,
        provider_version: Option<&str>,
        cache_id: Option<&str>,
        cache_version: Option<&str>,
        cache_digest: Option<&str>,
        organism: Option<&str>,
        taxon_id: Option<&str>,
        symbol_namespace: Option<&str>,
        review_status: Option<GeneSetResolutionReviewStatus>,
        filters: &[GeneSetProducerFilter],
    ) {
        if let Some(value) = Self::clean_option(provider_id) {
            self.provider_id = Some(value);
        }
        if let Some(value) = Self::clean_option(provider_label) {
            self.provider_label = Some(value);
        }
        if let Some(value) = Self::clean_option(provider_version) {
            self.provider_version = Some(value);
        }
        if let Some(value) = Self::clean_option(cache_id) {
            self.cache_id = Some(value);
        }
        if let Some(value) = Self::clean_option(cache_version) {
            self.cache_version = Some(value);
        }
        if let Some(value) = Self::clean_option(cache_digest) {
            self.cache_digest = Some(value);
        }
        if let Some(value) = Self::clean_option(organism) {
            self.organism = Some(value);
        }
        if let Some(value) = Self::clean_option(taxon_id) {
            self.taxon_id = Some(value);
        }
        if let Some(value) = Self::clean_option(symbol_namespace) {
            self.symbol_namespace = Some(value);
        }
        if let Some(value) = review_status {
            self.review_status = Some(value);
        }
        for filter in filters {
            if !self.filters.contains(filter) {
                self.filters.push(filter.clone());
            }
        }
    }

    fn clean_option(raw: Option<&str>) -> Option<String> {
        raw.map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
    }
}

#[derive(Clone, Debug, Default)]
struct GeneSetDirectListSelection {
    list_id: Option<String>,
    list_label: Option<String>,
    members: Vec<String>,
    metadata: GeneSetDirectListCacheMetadata,
    warnings: Vec<String>,
}

#[derive(Clone, Debug, Default)]
struct GeneSetOntologyAssignmentSelection {
    term_id: String,
    term_label: Option<String>,
    members: Vec<String>,
    metadata: GeneSetDirectListCacheMetadata,
    warnings: Vec<String>,
}

#[derive(Clone, Debug, Default)]
struct GeneSetCoRegulatedCohortSelection {
    members: Vec<String>,
    metadata: GeneSetDirectListCacheMetadata,
    dataset_ids: Vec<String>,
    contrast_labels: Vec<String>,
    condition_labels: Vec<String>,
    normalization_method: Option<String>,
    scoring_method: Option<String>,
    warnings: Vec<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GeneSetCoRegulatedThresholdMetric {
    Score,
    AbsoluteScore,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GeneSetCoRegulatedThresholdOperator {
    GreaterEqual,
    GreaterThan,
    LessEqual,
    LessThan,
}

#[derive(Clone, Copy, Debug)]
struct GeneSetCoRegulatedThreshold {
    metric: GeneSetCoRegulatedThresholdMetric,
    operator: GeneSetCoRegulatedThresholdOperator,
    value: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GeneSetCoRegulatedDirectionRule {
    Both,
    Positive,
    Negative,
}

impl GentleEngine {
    fn gene_set_normalize_lookup(raw: &str) -> String {
        crate::gene_groups::normalize_gene_group_lookup(raw)
    }

    fn gene_set_member_symbol(record: &GenomeGeneRecord) -> String {
        record
            .gene_name
            .as_deref()
            .or(record.gene_id.as_deref())
            .unwrap_or("unnamed_gene")
            .to_string()
    }

    fn gene_set_dedup_key(symbol: &str, gene_id: Option<&str>) -> String {
        gene_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| format!("gene_id:{}", value.to_ascii_lowercase()))
            .unwrap_or_else(|| format!("symbol:{}", Self::gene_set_normalize_lookup(symbol)))
    }

    fn gene_set_indexed_row(record: GenomeGeneRecord) -> IndexedGeneRow {
        let symbol = Self::gene_set_member_symbol(&record);
        let dedup_key = Self::gene_set_dedup_key(&symbol, record.gene_id.as_deref());
        IndexedGeneRow {
            record,
            dedup_key,
            symbol,
        }
    }

    fn gene_set_sort_indexed_genes(rows: &mut [IndexedGeneRow]) {
        rows.sort_by(|left, right| {
            left.record
                .chromosome
                .cmp(&right.record.chromosome)
                .then(left.record.start_1based.cmp(&right.record.start_1based))
                .then(left.record.gene_id.cmp(&right.record.gene_id))
                .then(left.symbol.cmp(&right.symbol))
        });
    }

    fn gene_set_sort_members(rows: &mut [GeneSetResolvedMember]) {
        rows.sort_by(
            |left, right| match (left.start_1based, right.start_1based) {
                (Some(left_start), Some(right_start)) => left
                    .chromosome
                    .cmp(&right.chromosome)
                    .then(left_start.cmp(&right_start))
                    .then(left.gene_id.cmp(&right.gene_id))
                    .then(left.symbol.cmp(&right.symbol)),
                (Some(_), None) => Ordering::Less,
                (None, Some(_)) => Ordering::Greater,
                (None, None) => left.symbol.cmp(&right.symbol),
            },
        );
    }

    fn gene_set_record_matches_query(record: &GenomeGeneRecord, query: &str) -> bool {
        let normalized = Self::gene_set_normalize_lookup(query);
        if normalized.is_empty() {
            return false;
        }
        record
            .gene_id
            .as_deref()
            .map(Self::gene_set_normalize_lookup)
            .is_some_and(|value| value == normalized)
            || record
                .gene_name
                .as_deref()
                .map(Self::gene_set_normalize_lookup)
                .is_some_and(|value| value == normalized)
    }

    fn gene_set_load_genome_genes(
        genome_id: Option<&str>,
        genome_catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<(Option<String>, Option<String>, Vec<IndexedGeneRow>), EngineError> {
        let Some(genome_id) = genome_id.map(str::trim).filter(|value| !value.is_empty()) else {
            return Ok((None, None, vec![]));
        };
        let effective_catalog_path =
            genome_catalog_path.unwrap_or(crate::genomes::default_catalog_discovery_token(false));
        let (catalog, catalog_label) =
            Self::open_reference_genome_catalog(Some(effective_catalog_path))?;
        let genes = catalog
            .list_gene_regions(genome_id, cache_dir)
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not load gene index for genome '{}': {}",
                    genome_id, e
                ),
                cause_chain: vec![],
            })?
            .into_iter()
            .map(Self::gene_set_indexed_row)
            .collect::<Vec<_>>();
        Ok((Some(genome_id.to_string()), Some(catalog_label), genes))
    }

    fn gene_set_resolve_genome_row<'a>(
        query: &str,
        genes: &'a [IndexedGeneRow],
    ) -> Vec<&'a IndexedGeneRow> {
        genes
            .iter()
            .filter(|row| Self::gene_set_record_matches_query(&row.record, query))
            .collect()
    }

    fn gene_set_member_from_indexed_gene(
        row: &IndexedGeneRow,
        member_status: Option<String>,
        confidence: Option<String>,
        contributing_group_ids: Vec<String>,
        provenance: Vec<GeneSetProvenanceRow>,
    ) -> GeneSetResolvedMember {
        GeneSetResolvedMember {
            dedup_key: row.dedup_key.clone(),
            symbol: row.symbol.clone(),
            gene_id: row.record.gene_id.clone(),
            aliases: vec![],
            chromosome: Some(row.record.chromosome.clone()),
            start_1based: Some(row.record.start_1based),
            end_1based: Some(row.record.end_1based),
            strand: row.record.strand.map(|value| value.to_string()),
            biotype: row.record.biotype.clone(),
            member_status,
            confidence,
            contributing_group_ids,
            provenance,
        }
    }

    fn gene_set_member_from_catalog_member(
        member: &GeneGroupMember,
        source_group_id: Option<&str>,
        source_path: Option<&str>,
        genes: &[IndexedGeneRow],
        warnings: &mut Vec<String>,
        unresolved: &mut Vec<GeneSetUnresolvedMember>,
        source_kind: &str,
    ) -> Option<GeneSetResolvedMember> {
        let symbol = member.symbol.trim();
        if symbol.is_empty() {
            return None;
        }
        match member.status.as_deref().map(str::trim).unwrap_or("") {
            "" | "included" => {}
            "draft" => warnings.push(format!(
                "Gene-set member '{}' is draft; including with warning",
                symbol
            )),
            "excluded" => {
                warnings.push(format!(
                    "Skipping excluded gene-set member '{}' from source '{}'",
                    symbol,
                    source_group_id.unwrap_or(source_kind)
                ));
                return None;
            }
            other => warnings.push(format!(
                "Unrecognized member status '{}' for '{}', treating as included",
                other, symbol
            )),
        }
        let provenance = vec![GeneSetProvenanceRow {
            source_kind: source_kind.to_string(),
            source_id: source_group_id.unwrap_or(source_kind).to_string(),
            source_label: None,
            source_path: source_path.map(str::to_string),
            note: member
                .provenance
                .clone()
                .or_else(|| member.evidence_note.clone()),
        }];
        let contributing_group_ids = source_group_id
            .map(|value| vec![value.to_string()])
            .unwrap_or_default();
        if genes.is_empty() {
            let dedup_key = Self::gene_set_dedup_key(symbol, member.gene_id.as_deref());
            return Some(GeneSetResolvedMember {
                dedup_key,
                symbol: symbol.to_string(),
                gene_id: member.gene_id.clone(),
                aliases: member.aliases.clone(),
                member_status: member.status.clone(),
                confidence: member.confidence.clone(),
                contributing_group_ids,
                provenance,
                ..GeneSetResolvedMember::default()
            });
        }
        let matches = Self::gene_set_resolve_genome_row(symbol, genes);
        match matches.as_slice() {
            [row] => {
                let mut resolved = Self::gene_set_member_from_indexed_gene(
                    row,
                    member.status.clone(),
                    member.confidence.clone(),
                    contributing_group_ids,
                    provenance,
                );
                resolved.aliases = member.aliases.clone();
                Some(resolved)
            }
            [] => {
                unresolved.push(GeneSetUnresolvedMember {
                    query: symbol.to_string(),
                    reason: "member did not resolve against the prepared genome gene index"
                        .to_string(),
                    source_kind: source_kind.to_string(),
                    source_id: source_group_id.map(str::to_string),
                });
                None
            }
            rows => {
                unresolved.push(GeneSetUnresolvedMember {
                    query: symbol.to_string(),
                    reason: format!(
                        "member resolved to {} loci; v1 does not guess among ambiguous loci",
                        rows.len()
                    ),
                    source_kind: source_kind.to_string(),
                    source_id: source_group_id.map(str::to_string),
                });
                None
            }
        }
    }

    fn gene_set_group_allowed(
        group: &LoadedGeneGroupRecord,
        allow_draft: bool,
        allow_deprecated: bool,
        warnings: &mut Vec<String>,
    ) -> Result<(), EngineError> {
        let status = group.record.curation_status.trim();
        match status {
            "" | "reviewed" | "included" | "curated" => Ok(()),
            "draft" if allow_draft => {
                warnings.push(format!(
                    "Using draft gene group '{}' because allow_draft=true",
                    group.record.id
                ));
                Ok(())
            }
            "draft" => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Gene group '{}' is draft; pass --allow-draft to use it as a gene set",
                    group.record.id
                ),
                cause_chain: vec![],
            }),
            "deprecated" if allow_deprecated => {
                warnings.push(format!(
                    "Using deprecated gene group '{}' because allow_deprecated=true",
                    group.record.id
                ));
                Ok(())
            }
            "deprecated" => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Gene group '{}' is deprecated; pass --allow-deprecated to use it as a gene set",
                    group.record.id
                ),
                cause_chain: vec![],
            }),
            other => {
                warnings.push(format!(
                    "Unrecognized curation_status '{}' for gene group '{}'; treating as usable",
                    other, group.record.id
                ));
                Ok(())
            }
        }
    }

    fn gene_set_insert_member(
        members: &mut Vec<GeneSetResolvedMember>,
        mut member: GeneSetResolvedMember,
    ) {
        if let Some(existing) = members
            .iter_mut()
            .find(|existing| existing.dedup_key == member.dedup_key)
        {
            for provenance in member.provenance.drain(..) {
                if !existing.provenance.contains(&provenance) {
                    existing.provenance.push(provenance);
                }
            }
            for group_id in member.contributing_group_ids {
                if !existing.contributing_group_ids.contains(&group_id) {
                    existing.contributing_group_ids.push(group_id);
                }
            }
            for alias in member.aliases {
                if !existing.aliases.contains(&alias) {
                    existing.aliases.push(alias);
                }
            }
        } else {
            members.push(member);
        }
    }

    fn gene_set_resolve_group_members(
        groups: &[LoadedGeneGroupRecord],
        genes: &[IndexedGeneRow],
        warnings: &mut Vec<String>,
        unresolved: &mut Vec<GeneSetUnresolvedMember>,
    ) -> Vec<GeneSetResolvedMember> {
        let mut members = vec![];
        for group in groups {
            for member in &group.record.members {
                if let Some(resolved) = Self::gene_set_member_from_catalog_member(
                    member,
                    Some(group.record.id.as_str()),
                    Some(group.source_path.as_str()),
                    genes,
                    warnings,
                    unresolved,
                    "catalog_group",
                ) {
                    Self::gene_set_insert_member(&mut members, resolved);
                }
            }
        }
        members
    }

    fn gene_set_resolve_explicit_members(
        requested: &[String],
        genes: &[IndexedGeneRow],
        unresolved: &mut Vec<GeneSetUnresolvedMember>,
    ) -> Vec<GeneSetResolvedMember> {
        let mut members = vec![];
        for raw in requested {
            let query = raw.trim();
            if query.is_empty() {
                continue;
            }
            let provenance = vec![GeneSetProvenanceRow {
                source_kind: "explicit_members".to_string(),
                source_id: query.to_string(),
                source_label: None,
                source_path: None,
                note: Some("User-supplied gene-set member".to_string()),
            }];
            if genes.is_empty() {
                Self::gene_set_insert_member(
                    &mut members,
                    GeneSetResolvedMember {
                        dedup_key: Self::gene_set_dedup_key(query, None),
                        symbol: query.to_string(),
                        provenance,
                        ..GeneSetResolvedMember::default()
                    },
                );
                continue;
            }
            let matches = Self::gene_set_resolve_genome_row(query, genes);
            match matches.as_slice() {
                [row] => Self::gene_set_insert_member(
                    &mut members,
                    Self::gene_set_member_from_indexed_gene(row, None, None, vec![], provenance),
                ),
                [] => unresolved.push(GeneSetUnresolvedMember {
                    query: query.to_string(),
                    reason: "explicit member did not resolve against the prepared genome gene index"
                        .to_string(),
                    source_kind: "explicit_members".to_string(),
                    source_id: None,
                }),
                rows => unresolved.push(GeneSetUnresolvedMember {
                    query: query.to_string(),
                    reason: format!(
                        "explicit member resolved to {} loci; v1 does not guess among ambiguous loci",
                        rows.len()
                    ),
                    source_kind: "explicit_members".to_string(),
                    source_id: None,
                }),
            }
        }
        members
    }

    fn stable_hash64(raw: &str) -> u64 {
        let mut hash = 0xcbf29ce484222325_u64;
        for byte in raw.as_bytes() {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(0x100000001b3);
        }
        hash
    }

    fn splitmix64_next(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }

    fn gene_set_resolve_random(
        genome_id: &str,
        count: usize,
        random_seed: u64,
        exclude_members: &[String],
        genes: &[IndexedGeneRow],
        gene_index_source: &str,
        warnings: &mut Vec<String>,
    ) -> (Vec<GeneSetResolvedMember>, GeneSetRandomProvenance) {
        let mut exclude_keys = BTreeSet::<String>::new();
        for token in exclude_members {
            let query = token.trim();
            if query.is_empty() {
                continue;
            }
            for row in Self::gene_set_resolve_genome_row(query, genes) {
                exclude_keys.insert(row.dedup_key.clone());
            }
            exclude_keys.insert(Self::gene_set_dedup_key(query, None));
        }
        let mut universe = genes
            .iter()
            .filter(|row| !exclude_keys.contains(&row.dedup_key))
            .cloned()
            .collect::<Vec<_>>();
        Self::gene_set_sort_indexed_genes(&mut universe);
        if universe.len() < count {
            warnings.push(format!(
                "Random gene-set universe has only {} non-foreground gene(s), fewer than requested {}",
                universe.len(),
                count
            ));
        }
        let mut state = random_seed ^ Self::stable_hash64(gene_index_source);
        for idx in (1..universe.len()).rev() {
            let swap_idx = (Self::splitmix64_next(&mut state) as usize) % (idx + 1);
            universe.swap(idx, swap_idx);
        }
        let mut chosen = universe.into_iter().take(count).collect::<Vec<_>>();
        Self::gene_set_sort_indexed_genes(&mut chosen);
        let members = chosen
            .iter()
            .map(|row| {
                Self::gene_set_member_from_indexed_gene(
                    row,
                    None,
                    None,
                    vec![],
                    vec![GeneSetProvenanceRow {
                        source_kind: "random".to_string(),
                        source_id: random_seed.to_string(),
                        source_label: None,
                        source_path: None,
                        note: Some(format!(
                            "Deterministic random sample from {}",
                            gene_index_source
                        )),
                    }],
                )
            })
            .collect::<Vec<_>>();
        (
            members,
            GeneSetRandomProvenance {
                genome_id: genome_id.to_string(),
                genome_build: Some(genome_id.to_string()),
                gene_index_source: gene_index_source.to_string(),
                random_seed,
                universe_size: genes.len(),
                foreground_exclusion_count: exclude_keys.len(),
            },
        )
    }

    fn gene_set_resolve_neighbors(
        anchor: &str,
        flank_gene_count: Option<usize>,
        flank_bp: Option<usize>,
        exclude_anchor: bool,
        genes: &[IndexedGeneRow],
        warnings: &mut Vec<String>,
        unresolved: &mut Vec<GeneSetUnresolvedMember>,
    ) -> Result<Vec<GeneSetResolvedMember>, EngineError> {
        let matches = Self::gene_set_resolve_genome_row(anchor, genes);
        let anchor_row = match matches.as_slice() {
            [row] => *row,
            [] => {
                unresolved.push(GeneSetUnresolvedMember {
                    query: anchor.to_string(),
                    reason:
                        "neighbor anchor did not resolve against the prepared genome gene index"
                            .to_string(),
                    source_kind: "genomic_neighbors".to_string(),
                    source_id: None,
                });
                return Ok(vec![]);
            }
            rows => {
                let loci = rows
                    .iter()
                    .map(|row| {
                        format!(
                            "{}:{}-{} ({})",
                            row.record.chromosome,
                            row.record.start_1based,
                            row.record.end_1based,
                            row.symbol
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(", ");
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Neighbor anchor '{}' resolved to {} loci; v1 requires a unique anchor: {}",
                        anchor,
                        rows.len(),
                        loci
                    ),
                    cause_chain: vec![],
                });
            }
        };
        let mut same_chrom = genes
            .iter()
            .filter(|row| row.record.chromosome == anchor_row.record.chromosome)
            .cloned()
            .collect::<Vec<_>>();
        same_chrom.sort_by(|left, right| {
            left.record
                .start_1based
                .cmp(&right.record.start_1based)
                .then(left.record.gene_id.cmp(&right.record.gene_id))
                .then(left.symbol.cmp(&right.symbol))
        });
        let Some(anchor_idx) = same_chrom
            .iter()
            .position(|row| row.dedup_key == anchor_row.dedup_key)
        else {
            return Ok(vec![]);
        };
        let selected = if let Some(window_bp) = flank_bp {
            let anchor_start = anchor_row.record.start_1based;
            same_chrom
                .into_iter()
                .filter(|row| {
                    let delta = row.record.start_1based.abs_diff(anchor_start);
                    delta <= window_bp && (!exclude_anchor || row.dedup_key != anchor_row.dedup_key)
                })
                .collect::<Vec<_>>()
        } else {
            let flank = flank_gene_count.unwrap_or(1);
            let start_idx = anchor_idx.saturating_sub(flank);
            let end_idx = (anchor_idx + flank + 1).min(same_chrom.len());
            if anchor_idx < flank {
                warnings.push(format!(
                    "Neighbor anchor '{}' is near the start of chromosome {}; only {} upstream gene(s) available",
                    anchor,
                    anchor_row.record.chromosome,
                    anchor_idx
                ));
            }
            let downstream_available = same_chrom.len().saturating_sub(anchor_idx + 1);
            if downstream_available < flank {
                warnings.push(format!(
                    "Neighbor anchor '{}' is near the end of chromosome {}; only {} downstream gene(s) available",
                    anchor,
                    anchor_row.record.chromosome,
                    downstream_available
                ));
            }
            same_chrom[start_idx..end_idx]
                .iter()
                .filter(|row| !exclude_anchor || row.dedup_key != anchor_row.dedup_key)
                .cloned()
                .collect::<Vec<_>>()
        };
        Ok(selected
            .iter()
            .map(|row| {
                Self::gene_set_member_from_indexed_gene(
                    row,
                    None,
                    None,
                    vec![],
                    vec![GeneSetProvenanceRow {
                        source_kind: "genomic_neighbors".to_string(),
                        source_id: anchor.to_string(),
                        source_label: None,
                        source_path: None,
                        note: Some("Same-chromosome prepared-genome neighborhood".to_string()),
                    }],
                )
            })
            .collect())
    }

    fn gene_set_direct_list_string(value: &Value, keys: &[&str]) -> Option<String> {
        keys.iter()
            .find_map(|key| {
                value.get(*key).and_then(|raw| match raw {
                    Value::String(text) => Some(text.trim().to_string()),
                    Value::Number(number) => Some(number.to_string()),
                    _ => None,
                })
            })
            .filter(|value| !value.is_empty())
    }

    fn gene_set_direct_list_review_status(raw: &str) -> Option<GeneSetResolutionReviewStatus> {
        match raw.trim().to_ascii_lowercase().replace('-', "_").as_str() {
            "unreviewed" | "" => Some(GeneSetResolutionReviewStatus::Unreviewed),
            "reviewed" => Some(GeneSetResolutionReviewStatus::Reviewed),
            "included" => Some(GeneSetResolutionReviewStatus::Included),
            "draft" => Some(GeneSetResolutionReviewStatus::Draft),
            "deprecated" => Some(GeneSetResolutionReviewStatus::Deprecated),
            _ => None,
        }
    }

    fn gene_set_direct_list_metadata_from_value(
        value: &Value,
        warnings: &mut Vec<String>,
    ) -> GeneSetDirectListCacheMetadata {
        let review_status = Self::gene_set_direct_list_string(
            value,
            &["review_status", "curation_status", "status"],
        )
        .and_then(|raw| {
            let parsed = Self::gene_set_direct_list_review_status(&raw);
            if parsed.is_none() {
                warnings.push(format!(
                    "Unrecognized direct-list review status '{}'; defaulting to unreviewed",
                    raw
                ));
            }
            parsed
        });
        GeneSetDirectListCacheMetadata {
            provider_id: Self::gene_set_direct_list_string(
                value,
                &["provider_id", "provider", "source_id"],
            ),
            provider_label: Self::gene_set_direct_list_string(
                value,
                &["provider_label", "provider_name", "source_label"],
            ),
            provider_version: Self::gene_set_direct_list_string(
                value,
                &["provider_version", "source_version"],
            ),
            cache_id: Self::gene_set_direct_list_string(value, &["cache_id"]),
            cache_version: Self::gene_set_direct_list_string(
                value,
                &["cache_version", "cache_format_version"],
            ),
            cache_digest: Self::gene_set_direct_list_string(
                value,
                &["cache_digest", "digest", "sha256"],
            ),
            organism: Self::gene_set_direct_list_string(value, &["organism", "species"]),
            taxon_id: Self::gene_set_direct_list_string(value, &["taxon_id", "tax_id"]),
            symbol_namespace: Self::gene_set_direct_list_string(
                value,
                &["symbol_namespace", "namespace"],
            ),
            review_status,
            filters: Self::gene_set_direct_list_filters_from_value(value, warnings),
        }
    }

    fn gene_set_direct_list_filters_from_value(
        value: &Value,
        warnings: &mut Vec<String>,
    ) -> Vec<GeneSetProducerFilter> {
        let Some(filters) = value.get("filters").and_then(Value::as_array) else {
            return vec![];
        };
        let mut parsed = vec![];
        for filter in filters {
            match filter {
                Value::String(raw) => {
                    if let Some(filter) = Self::gene_set_direct_list_filter_from_assignment(raw) {
                        parsed.push(filter);
                    } else {
                        warnings.push(format!(
                            "Skipping malformed direct-list filter '{}'; expected FIELD=VALUE",
                            raw
                        ));
                    }
                }
                Value::Object(_) => {
                    let field = Self::gene_set_direct_list_string(filter, &["field"]);
                    let value = Self::gene_set_direct_list_string(filter, &["value"]);
                    if let (Some(field), Some(value)) = (field, value) {
                        parsed.push(GeneSetProducerFilter {
                            field,
                            operator: Self::gene_set_direct_list_string(filter, &["operator"])
                                .unwrap_or_else(|| "equals".to_string()),
                            value,
                        });
                    } else {
                        warnings.push(
                            "Skipping malformed direct-list filter object; expected field and value"
                                .to_string(),
                        );
                    }
                }
                _ => warnings
                    .push("Skipping non-string/non-object direct-list filter entry".to_string()),
            }
        }
        parsed
    }

    fn gene_set_direct_list_filter_from_assignment(raw: &str) -> Option<GeneSetProducerFilter> {
        let (field, value) = raw.split_once('=')?;
        let field = field.trim();
        let value = value.trim();
        if field.is_empty() || value.is_empty() {
            return None;
        }
        Some(GeneSetProducerFilter {
            field: field.to_string(),
            operator: "equals".to_string(),
            value: value.to_string(),
        })
    }

    fn gene_set_direct_list_member_from_value(
        value: &Value,
        warnings: &mut Vec<String>,
    ) -> Option<String> {
        match value {
            Value::String(text) => {
                let trimmed = text.trim();
                (!trimmed.is_empty()).then(|| trimmed.to_string())
            }
            Value::Object(_) => Self::gene_set_direct_list_string(
                value,
                &[
                    "symbol",
                    "gene_symbol",
                    "member",
                    "gene",
                    "query",
                    "gene_id",
                    "id",
                ],
            ),
            _ => {
                warnings.push("Skipping non-string/non-object direct-list member".to_string());
                None
            }
        }
    }

    fn gene_set_direct_list_members_from_value(
        value: &Value,
        warnings: &mut Vec<String>,
    ) -> Vec<String> {
        value
            .as_array()
            .map(|members| {
                members
                    .iter()
                    .filter_map(|member| {
                        Self::gene_set_direct_list_member_from_value(member, warnings)
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default()
    }

    fn gene_set_direct_list_query_matches(value: &Value, query: &str) -> bool {
        let normalized_query = Self::gene_set_normalize_lookup(query);
        ["id", "label", "name", "query", "query_id"]
            .iter()
            .filter_map(|key| Self::gene_set_direct_list_string(value, &[*key]))
            .any(|candidate| Self::gene_set_normalize_lookup(&candidate) == normalized_query)
    }

    fn gene_set_direct_list_select_json(
        path: &str,
        text: &str,
        query: Option<&str>,
    ) -> Result<GeneSetDirectListSelection, EngineError> {
        let value: Value = serde_json::from_str(text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse direct gene-list cache '{}': {}", path, e),
            cause_chain: vec![],
        })?;
        let mut warnings = vec![];
        let mut root_metadata =
            Self::gene_set_direct_list_metadata_from_value(&value, &mut warnings);
        if value
            .get("schema")
            .and_then(Value::as_str)
            .is_some_and(|schema| schema == GENE_SET_DIRECT_LIST_CACHE_SCHEMA)
            && root_metadata.cache_version.is_none()
        {
            root_metadata.cache_version = Some("1".to_string());
        }

        if let Some(lists) = value
            .get("lists")
            .or_else(|| value.get("sets"))
            .and_then(Value::as_array)
        {
            let selected = match query {
                Some(query) => lists
                    .iter()
                    .find(|entry| Self::gene_set_direct_list_query_matches(entry, query)),
                None if lists.len() == 1 => lists.first(),
                None => {
                    let ids = lists
                        .iter()
                        .filter_map(|entry| {
                            Self::gene_set_direct_list_string(entry, &["id", "label", "name"])
                        })
                        .collect::<Vec<_>>()
                        .join(", ");
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Direct gene-list cache '{}' contains multiple lists; pass --query (available: {})",
                            path,
                            if ids.is_empty() {
                                "unlabeled".to_string()
                            } else {
                                ids
                            }
                        ),
                        cause_chain: vec![],
                    });
                }
            };
            let Some(selected) = selected else {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Direct gene-list cache '{}' does not contain list '{}'",
                        path,
                        query.unwrap_or("")
                    ),
                    cause_chain: vec![],
                });
            };
            let mut metadata =
                Self::gene_set_direct_list_metadata_from_value(selected, &mut warnings);
            metadata.merge_missing_from(&root_metadata);
            let members = selected
                .get("members")
                .map(|members| {
                    Self::gene_set_direct_list_members_from_value(members, &mut warnings)
                })
                .unwrap_or_default();
            return Ok(GeneSetDirectListSelection {
                list_id: Self::gene_set_direct_list_string(selected, &["id", "query_id", "name"])
                    .or_else(|| query.map(str::to_string)),
                list_label: Self::gene_set_direct_list_string(selected, &["label", "name"]),
                members,
                metadata,
                warnings,
            });
        }

        if let Some(members_value) = value.get("members") {
            if let Some(query) = query
                && !Self::gene_set_direct_list_query_matches(&value, query)
            {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Direct gene-list cache '{}' top-level members do not match query '{}'",
                        path, query
                    ),
                    cause_chain: vec![],
                });
            }
            let members =
                Self::gene_set_direct_list_members_from_value(members_value, &mut warnings);
            return Ok(GeneSetDirectListSelection {
                list_id: Self::gene_set_direct_list_string(&value, &["id", "query_id", "name"])
                    .or_else(|| query.map(str::to_string)),
                list_label: Self::gene_set_direct_list_string(&value, &["label", "name"]),
                members,
                metadata: root_metadata,
                warnings,
            });
        }

        Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Direct gene-list cache '{}' must contain top-level members[] or lists[]",
                path
            ),
            cause_chain: vec![],
        })
    }

    fn gene_set_direct_list_select_tsv(
        path: &str,
        text: &str,
        query: Option<&str>,
    ) -> Result<GeneSetDirectListSelection, EngineError> {
        let delimiter = if path.to_ascii_lowercase().ends_with(".csv") {
            b','
        } else {
            b'\t'
        };
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(delimiter)
            .flexible(true)
            .from_reader(text.as_bytes());
        let headers = reader
            .headers()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read direct gene-list cache header '{}': {}",
                    path, e
                ),
                cause_chain: vec![],
            })?
            .iter()
            .map(|header| header.trim().to_ascii_lowercase())
            .collect::<Vec<_>>();
        let find_col = |names: &[&str]| -> Option<usize> {
            names
                .iter()
                .find_map(|name| headers.iter().position(|header| header == *name))
        };
        let list_col = find_col(&["list_id", "set_id", "query_id", "group_id"]);
        let list_label_col = find_col(&["list_label", "set_label", "label"]);
        let symbol_col = find_col(&["symbol", "gene_symbol", "member", "gene"]);
        let gene_id_col = find_col(&["gene_id", "id"]);
        let provider_id_col = find_col(&["provider_id", "provider"]);
        let provider_label_col = find_col(&["provider_label", "provider_name"]);
        let provider_version_col = find_col(&["provider_version"]);
        let cache_id_col = find_col(&["cache_id"]);
        let cache_version_col = find_col(&["cache_version", "cache_format_version"]);
        let cache_digest_col = find_col(&["cache_digest", "digest", "sha256"]);
        let organism_col = find_col(&["organism", "species"]);
        let taxon_id_col = find_col(&["taxon_id", "tax_id"]);
        let namespace_col = find_col(&["symbol_namespace", "namespace"]);
        let mut selection = GeneSetDirectListSelection::default();
        let mut list_ids = BTreeSet::<String>::new();

        let get = |record: &csv::StringRecord, col: Option<usize>| -> Option<String> {
            col.and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string)
        };

        for row in reader.records() {
            let record = row.map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read direct gene-list cache row '{}': {}",
                    path, e
                ),
                cause_chain: vec![],
            })?;
            let row_list_id = get(&record, list_col);
            if let Some(list_id) = row_list_id.as_deref() {
                list_ids.insert(list_id.to_string());
            }
            if let Some(query) = query {
                let Some(row_list_id) = row_list_id.as_deref() else {
                    continue;
                };
                if Self::gene_set_normalize_lookup(row_list_id)
                    != Self::gene_set_normalize_lookup(query)
                {
                    continue;
                }
            }
            if let Some(value) = get(&record, list_col)
                && selection.list_id.is_none()
            {
                selection.list_id = Some(value);
            }
            if let Some(value) = get(&record, list_label_col)
                && selection.list_label.is_none()
            {
                selection.list_label = Some(value);
            }
            let member = get(&record, symbol_col).or_else(|| get(&record, gene_id_col));
            if let Some(member) = member {
                selection.members.push(member);
            } else {
                selection
                    .warnings
                    .push("Skipping direct-list TSV row without symbol or gene_id".to_string());
            }
            if selection.metadata.provider_id.is_none() {
                selection.metadata.provider_id = get(&record, provider_id_col);
            }
            if selection.metadata.provider_label.is_none() {
                selection.metadata.provider_label = get(&record, provider_label_col);
            }
            if selection.metadata.provider_version.is_none() {
                selection.metadata.provider_version = get(&record, provider_version_col);
            }
            if selection.metadata.cache_id.is_none() {
                selection.metadata.cache_id = get(&record, cache_id_col);
            }
            if selection.metadata.cache_version.is_none() {
                selection.metadata.cache_version = get(&record, cache_version_col);
            }
            if selection.metadata.cache_digest.is_none() {
                selection.metadata.cache_digest = get(&record, cache_digest_col);
            }
            if selection.metadata.organism.is_none() {
                selection.metadata.organism = get(&record, organism_col);
            }
            if selection.metadata.taxon_id.is_none() {
                selection.metadata.taxon_id = get(&record, taxon_id_col);
            }
            if selection.metadata.symbol_namespace.is_none() {
                selection.metadata.symbol_namespace = get(&record, namespace_col);
            }
        }

        if query.is_none() && list_ids.len() > 1 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Direct gene-list TSV cache '{}' contains multiple list_id values; pass --query (available: {})",
                    path,
                    list_ids.into_iter().collect::<Vec<_>>().join(", ")
                ),
                cause_chain: vec![],
            });
        }
        if selection.list_id.is_none() {
            selection.list_id = query.map(str::to_string);
        }
        Ok(selection)
    }

    fn gene_set_direct_list_select_cache(
        cache_path: &str,
        query: Option<&str>,
    ) -> Result<GeneSetDirectListSelection, EngineError> {
        let text = std::fs::read_to_string(cache_path).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Could not read direct gene-list cache '{}': {}",
                cache_path, e
            ),
            cause_chain: vec![],
        })?;
        let trimmed = text.trim_start();
        if trimmed.starts_with('{') || trimmed.starts_with('[') {
            Self::gene_set_direct_list_select_json(cache_path, &text, query)
        } else {
            Self::gene_set_direct_list_select_tsv(cache_path, &text, query)
        }
    }

    fn gene_set_direct_list_validate_cache_version(
        cache_version: Option<&str>,
    ) -> Result<(), EngineError> {
        let Some(version) = cache_version
            .map(str::trim)
            .filter(|value| !value.is_empty())
        else {
            return Ok(());
        };
        let major = version
            .split(['.', '-'])
            .next()
            .and_then(|value| value.parse::<usize>().ok());
        if major.is_some_and(|value| value > 1) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported direct gene-list cache_version '{}'; this engine supports major version 1",
                    version
                ),
                cause_chain: vec![],
            });
        }
        Ok(())
    }

    fn gene_set_direct_list_fallback_cache_id(cache_path: &str) -> String {
        Path::new(cache_path)
            .file_stem()
            .and_then(|stem| stem.to_str())
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("direct_gene_list_cache")
            .to_string()
    }

    fn gene_set_ontology_assignment_validate_cache_version(
        cache_version: Option<&str>,
    ) -> Result<(), EngineError> {
        let Some(version) = cache_version
            .map(str::trim)
            .filter(|value| !value.is_empty())
        else {
            return Ok(());
        };
        let major = version
            .split(['.', '-'])
            .next()
            .and_then(|value| value.parse::<usize>().ok());
        if major.is_some_and(|value| value > 1) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported ontology assignment cache_version '{}'; this engine supports major version 1",
                    version
                ),
                cause_chain: vec![],
            });
        }
        Ok(())
    }

    fn gene_set_ontology_assignment_fallback_cache_id(cache_path: &str) -> String {
        Path::new(cache_path)
            .file_stem()
            .and_then(|stem| stem.to_str())
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("ontology_assignment_cache")
            .to_string()
    }

    fn gene_set_ontology_split_term(raw: &str) -> (Option<String>, String) {
        let trimmed = raw.trim();
        if let Some((namespace, id)) = trimmed.split_once(':') {
            let namespace = namespace.trim();
            let id = id.trim();
            (
                (!namespace.is_empty()).then(|| namespace.to_ascii_uppercase()),
                if id.is_empty() {
                    trimmed.to_string()
                } else {
                    id.to_string()
                },
            )
        } else {
            (None, trimmed.to_string())
        }
    }

    fn gene_set_ontology_term_matches_candidates(
        term: &str,
        ontology_namespace: Option<&str>,
        candidate_namespace: Option<&str>,
        candidate_terms: &[String],
    ) -> bool {
        let (term_namespace, term_local_id) = Self::gene_set_ontology_split_term(term);
        let requested_namespace = ontology_namespace
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_ascii_uppercase)
            .or(term_namespace);
        let requested_full = requested_namespace
            .as_deref()
            .map(|namespace| format!("{}:{}", namespace, term_local_id));
        let normalized_term = Self::gene_set_normalize_lookup(term);
        let normalized_local_id = Self::gene_set_normalize_lookup(&term_local_id);
        let normalized_full = requested_full
            .as_deref()
            .map(Self::gene_set_normalize_lookup);
        let candidate_namespace = candidate_namespace
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_ascii_uppercase);

        candidate_terms
            .iter()
            .map(|candidate| candidate.trim())
            .filter(|candidate| !candidate.is_empty())
            .any(|candidate| {
                let normalized_candidate = Self::gene_set_normalize_lookup(candidate);
                if normalized_candidate == normalized_term
                    || normalized_full
                        .as_deref()
                        .is_some_and(|value| normalized_candidate == value)
                {
                    return true;
                }
                let (candidate_term_namespace, candidate_local_id) =
                    Self::gene_set_ontology_split_term(candidate);
                let effective_candidate_namespace =
                    candidate_term_namespace.or_else(|| candidate_namespace.clone());
                if Self::gene_set_normalize_lookup(&candidate_local_id) != normalized_local_id {
                    return false;
                }
                match requested_namespace.as_deref() {
                    Some(namespace) => effective_candidate_namespace
                        .as_deref()
                        .is_some_and(|candidate_namespace| candidate_namespace == namespace),
                    None => true,
                }
            })
    }

    fn gene_set_ontology_term_matches_value(
        value: &Value,
        term: &str,
        ontology_namespace: Option<&str>,
    ) -> bool {
        let namespace = Self::gene_set_direct_list_string(
            value,
            &[
                "ontology_namespace",
                "term_namespace",
                "namespace",
                "source_namespace",
            ],
        );
        let candidates = [
            "term",
            "term_id",
            "go_id",
            "id",
            "ontology_id",
            "accession",
            "query_id",
        ]
        .iter()
        .filter_map(|key| Self::gene_set_direct_list_string(value, &[*key]))
        .collect::<Vec<_>>();
        Self::gene_set_ontology_term_matches_candidates(
            term,
            ontology_namespace,
            namespace.as_deref(),
            &candidates,
        )
    }

    fn gene_set_ontology_filter_keys(field: &str) -> Vec<String> {
        let normalized = field.trim().to_ascii_lowercase().replace('-', "_");
        match normalized.as_str() {
            "evidence_code" | "evidence" => vec![
                "evidence_code".to_string(),
                "evidence_codes".to_string(),
                "evidence".to_string(),
                "go_evidence_code".to_string(),
            ],
            "assigned_by" | "provider" => vec![
                "assigned_by".to_string(),
                "provider".to_string(),
                "provider_id".to_string(),
            ],
            other => vec![other.to_string()],
        }
    }

    fn gene_set_ontology_value_strings(value: &Value, keys: &[String]) -> Vec<String> {
        let mut out = vec![];
        for key in keys {
            if let Some(raw) = value.get(key) {
                match raw {
                    Value::String(text) => {
                        let trimmed = text.trim();
                        if !trimmed.is_empty() {
                            out.push(trimmed.to_string());
                        }
                    }
                    Value::Number(number) => out.push(number.to_string()),
                    Value::Array(values) => {
                        for entry in values {
                            match entry {
                                Value::String(text) => {
                                    let trimmed = text.trim();
                                    if !trimmed.is_empty() {
                                        out.push(trimmed.to_string());
                                    }
                                }
                                Value::Number(number) => out.push(number.to_string()),
                                _ => {}
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
        out
    }

    fn gene_set_ontology_value_matches_filters(
        value: &Value,
        filters: &[GeneSetProducerFilter],
    ) -> bool {
        filters.iter().all(|filter| {
            let operator = filter.operator.trim().to_ascii_lowercase();
            if operator != "equals" && operator != "=" {
                return false;
            }
            let keys = Self::gene_set_ontology_filter_keys(&filter.field);
            Self::gene_set_ontology_value_strings(value, &keys)
                .iter()
                .any(|candidate| {
                    Self::gene_set_normalize_lookup(candidate)
                        == Self::gene_set_normalize_lookup(&filter.value)
                })
        })
    }

    fn gene_set_ontology_assignment_select_json(
        path: &str,
        text: &str,
        term: &str,
        ontology_namespace: Option<&str>,
        filters: &[GeneSetProducerFilter],
    ) -> Result<GeneSetOntologyAssignmentSelection, EngineError> {
        let value: Value = serde_json::from_str(text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Could not parse ontology assignment cache '{}': {}",
                path, e
            ),
            cause_chain: vec![],
        })?;
        let mut warnings = vec![];
        let mut root_metadata =
            Self::gene_set_direct_list_metadata_from_value(&value, &mut warnings);
        if value
            .get("schema")
            .and_then(Value::as_str)
            .is_some_and(|schema| schema == GENE_SET_ONTOLOGY_ASSIGNMENT_CACHE_SCHEMA)
            && root_metadata.cache_version.is_none()
        {
            root_metadata.cache_version = Some("1".to_string());
        }
        let entries = value
            .get("assignments")
            .or_else(|| value.get("terms"))
            .or_else(|| value.get("ontology_terms"))
            .or_else(|| value.get("entries"))
            .or_else(|| value.get("rows"))
            .and_then(Value::as_array)
            .map(|rows| rows.iter().collect::<Vec<_>>())
            .unwrap_or_else(|| vec![&value]);
        let mut selection = GeneSetOntologyAssignmentSelection {
            term_id: term.to_string(),
            metadata: root_metadata.clone(),
            warnings,
            ..GeneSetOntologyAssignmentSelection::default()
        };
        let mut matched_term = false;

        for entry in entries {
            if !Self::gene_set_ontology_term_matches_value(entry, term, ontology_namespace) {
                continue;
            }
            matched_term = true;
            if !Self::gene_set_ontology_value_matches_filters(entry, filters) {
                continue;
            }
            if selection.term_id == term {
                if let Some(term_id) = Self::gene_set_direct_list_string(
                    entry,
                    &["term", "term_id", "go_id", "id", "ontology_id", "accession"],
                ) {
                    selection.term_id = term_id;
                }
            }
            if selection.term_label.is_none() {
                selection.term_label = Self::gene_set_direct_list_string(
                    entry,
                    &["term_label", "label", "name", "term_name"],
                );
            }
            let mut metadata = Self::gene_set_direct_list_metadata_from_value(entry, &mut selection.warnings);
            metadata.merge_missing_from(&root_metadata);
            selection.metadata.merge_missing_from(&metadata);

            if let Some(members) = entry.get("members") {
                selection.members.extend(Self::gene_set_direct_list_members_from_value(
                    members,
                    &mut selection.warnings,
                ));
            } else if let Some(member) =
                Self::gene_set_direct_list_member_from_value(entry, &mut selection.warnings)
            {
                selection.members.push(member);
            }
        }

        if !matched_term {
            selection.warnings.push(format!(
                "Ontology assignment cache '{}' has no rows for term '{}'",
                path, term
            ));
        } else if selection.members.is_empty() {
            selection.warnings.push(format!(
                "Ontology assignment cache '{}' has no members for term '{}' after filters",
                path, term
            ));
        }
        Ok(selection)
    }

    fn gene_set_ontology_record_values(
        headers: &[String],
        record: &csv::StringRecord,
        field: &str,
    ) -> Vec<String> {
        Self::gene_set_ontology_filter_keys(field)
            .iter()
            .flat_map(|key| {
                headers
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, header)| {
                        (header == key)
                            .then(|| record.get(idx))
                            .flatten()
                            .map(str::trim)
                            .filter(|value| !value.is_empty())
                            .map(str::to_string)
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    fn gene_set_ontology_record_matches_filters(
        headers: &[String],
        record: &csv::StringRecord,
        filters: &[GeneSetProducerFilter],
    ) -> bool {
        filters.iter().all(|filter| {
            let operator = filter.operator.trim().to_ascii_lowercase();
            if operator != "equals" && operator != "=" {
                return false;
            }
            Self::gene_set_ontology_record_values(headers, record, &filter.field)
                .iter()
                .any(|candidate| {
                    Self::gene_set_normalize_lookup(candidate)
                        == Self::gene_set_normalize_lookup(&filter.value)
                })
        })
    }

    fn gene_set_ontology_assignment_select_tsv(
        path: &str,
        text: &str,
        term: &str,
        ontology_namespace: Option<&str>,
        filters: &[GeneSetProducerFilter],
    ) -> Result<GeneSetOntologyAssignmentSelection, EngineError> {
        let delimiter = if path.to_ascii_lowercase().ends_with(".csv") {
            b','
        } else {
            b'\t'
        };
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(delimiter)
            .flexible(true)
            .from_reader(text.as_bytes());
        let headers = reader
            .headers()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read ontology assignment cache header '{}': {}",
                    path, e
                ),
                cause_chain: vec![],
            })?
            .iter()
            .map(|header| header.trim().to_ascii_lowercase().replace('-', "_"))
            .collect::<Vec<_>>();
        let find_col = |names: &[&str]| -> Option<usize> {
            names
                .iter()
                .find_map(|name| headers.iter().position(|header| header == *name))
        };
        let namespace_col = find_col(&["ontology_namespace", "term_namespace", "source_namespace"]);
        let term_col = find_col(&["term", "term_id", "go_id", "ontology_id", "accession"]);
        let id_col = find_col(&["id", "query_id"]);
        let term_label_col = find_col(&["term_label", "term_name", "label", "name"]);
        let symbol_col = find_col(&["symbol", "gene_symbol", "member", "gene"]);
        let gene_id_col = find_col(&["gene_id"]);
        let provider_id_col = find_col(&["provider_id", "provider"]);
        let provider_label_col = find_col(&["provider_label", "provider_name"]);
        let provider_version_col = find_col(&["provider_version"]);
        let cache_id_col = find_col(&["cache_id"]);
        let cache_version_col = find_col(&["cache_version", "cache_format_version"]);
        let cache_digest_col = find_col(&["cache_digest", "digest", "sha256"]);
        let organism_col = find_col(&["organism", "species"]);
        let taxon_id_col = find_col(&["taxon_id", "tax_id"]);
        let symbol_namespace_col = find_col(&["symbol_namespace", "namespace"]);
        let get = |record: &csv::StringRecord, col: Option<usize>| -> Option<String> {
            col.and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string)
        };
        let mut selection = GeneSetOntologyAssignmentSelection {
            term_id: term.to_string(),
            ..GeneSetOntologyAssignmentSelection::default()
        };
        let mut matched_term = false;

        for row in reader.records() {
            let record = row.map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read ontology assignment cache row '{}': {}",
                    path, e
                ),
                cause_chain: vec![],
            })?;
            let namespace = get(&record, namespace_col);
            let candidates = [get(&record, term_col), get(&record, id_col)]
                .into_iter()
                .flatten()
                .collect::<Vec<_>>();
            if !Self::gene_set_ontology_term_matches_candidates(
                term,
                ontology_namespace,
                namespace.as_deref(),
                &candidates,
            ) {
                continue;
            }
            matched_term = true;
            if !Self::gene_set_ontology_record_matches_filters(&headers, &record, filters) {
                continue;
            }
            if selection.term_id == term {
                if let Some(value) = get(&record, term_col).or_else(|| get(&record, id_col)) {
                    selection.term_id = value;
                }
            }
            if selection.term_label.is_none() {
                selection.term_label = get(&record, term_label_col);
            }
            if let Some(member) = get(&record, symbol_col).or_else(|| get(&record, gene_id_col)) {
                selection.members.push(member);
            } else {
                selection.warnings.push(
                    "Skipping ontology assignment TSV row without symbol or gene_id".to_string(),
                );
            }
            if selection.metadata.provider_id.is_none() {
                selection.metadata.provider_id = get(&record, provider_id_col);
            }
            if selection.metadata.provider_label.is_none() {
                selection.metadata.provider_label = get(&record, provider_label_col);
            }
            if selection.metadata.provider_version.is_none() {
                selection.metadata.provider_version = get(&record, provider_version_col);
            }
            if selection.metadata.cache_id.is_none() {
                selection.metadata.cache_id = get(&record, cache_id_col);
            }
            if selection.metadata.cache_version.is_none() {
                selection.metadata.cache_version = get(&record, cache_version_col);
            }
            if selection.metadata.cache_digest.is_none() {
                selection.metadata.cache_digest = get(&record, cache_digest_col);
            }
            if selection.metadata.organism.is_none() {
                selection.metadata.organism = get(&record, organism_col);
            }
            if selection.metadata.taxon_id.is_none() {
                selection.metadata.taxon_id = get(&record, taxon_id_col);
            }
            if selection.metadata.symbol_namespace.is_none() {
                selection.metadata.symbol_namespace = get(&record, symbol_namespace_col);
            }
        }

        if !matched_term {
            selection.warnings.push(format!(
                "Ontology assignment cache '{}' has no rows for term '{}'",
                path, term
            ));
        } else if selection.members.is_empty() {
            selection.warnings.push(format!(
                "Ontology assignment cache '{}' has no members for term '{}' after filters",
                path, term
            ));
        }
        Ok(selection)
    }

    fn gene_set_ontology_assignment_select_cache(
        cache_path: &str,
        term: &str,
        ontology_namespace: Option<&str>,
        filters: &[GeneSetProducerFilter],
    ) -> Result<GeneSetOntologyAssignmentSelection, EngineError> {
        let text = std::fs::read_to_string(cache_path).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Could not read ontology assignment cache '{}': {}",
                cache_path, e
            ),
            cause_chain: vec![],
        })?;
        let trimmed = text.trim_start();
        if trimmed.starts_with('{') || trimmed.starts_with('[') {
            Self::gene_set_ontology_assignment_select_json(
                cache_path,
                &text,
                term,
                ontology_namespace,
                filters,
            )
        } else {
            Self::gene_set_ontology_assignment_select_tsv(
                cache_path,
                &text,
                term,
                ontology_namespace,
                filters,
            )
        }
    }

    fn gene_set_co_regulated_validate_cache_version(
        cache_version: Option<&str>,
    ) -> Result<(), EngineError> {
        let Some(version) = cache_version
            .map(str::trim)
            .filter(|value| !value.is_empty())
        else {
            return Ok(());
        };
        let major = version
            .split(['.', '-'])
            .next()
            .and_then(|value| value.parse::<usize>().ok());
        if major.is_some_and(|value| value > 1) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported co-regulated cohort cache_version '{}'; this engine supports major version 1",
                    version
                ),
                cause_chain: vec![],
            });
        }
        Ok(())
    }

    fn gene_set_co_regulated_fallback_cache_id(cache_path: &str) -> String {
        Path::new(cache_path)
            .file_stem()
            .and_then(|stem| stem.to_str())
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("co_regulated_cohort_cache")
            .to_string()
    }

    fn gene_set_co_regulated_string_list(value: &Value, keys: &[&str]) -> Vec<String> {
        let mut out = vec![];
        for key in keys {
            if let Some(raw) = value.get(*key) {
                match raw {
                    Value::String(text) => {
                        for item in text.split(',') {
                            let trimmed = item.trim();
                            if !trimmed.is_empty() {
                                out.push(trimmed.to_string());
                            }
                        }
                    }
                    Value::Number(number) => out.push(number.to_string()),
                    Value::Array(values) => {
                        for entry in values {
                            match entry {
                                Value::String(text) => {
                                    let trimmed = text.trim();
                                    if !trimmed.is_empty() {
                                        out.push(trimmed.to_string());
                                    }
                                }
                                Value::Number(number) => out.push(number.to_string()),
                                _ => {}
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
        out.sort();
        out.dedup();
        out
    }

    fn gene_set_co_regulated_merge_strings(target: &mut Vec<String>, source: &[String]) {
        for value in source {
            if !target.contains(value) {
                target.push(value.clone());
            }
        }
    }

    fn gene_set_co_regulated_parse_threshold(
        raw: &str,
    ) -> Result<GeneSetCoRegulatedThreshold, EngineError> {
        let compact = raw.split_whitespace().collect::<String>();
        for (operator_token, operator) in [
            (
                ">=",
                GeneSetCoRegulatedThresholdOperator::GreaterEqual,
            ),
            ("<=", GeneSetCoRegulatedThresholdOperator::LessEqual),
            (">", GeneSetCoRegulatedThresholdOperator::GreaterThan),
            ("<", GeneSetCoRegulatedThresholdOperator::LessThan),
        ] {
            if let Some((metric_raw, value_raw)) = compact.split_once(operator_token) {
                let metric = match metric_raw.to_ascii_lowercase().as_str() {
                    "score" | "value" => GeneSetCoRegulatedThresholdMetric::Score,
                    "abs" | "absolute" | "abs_score" | "absolute_score" => {
                        GeneSetCoRegulatedThresholdMetric::AbsoluteScore
                    }
                    _ => {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Unsupported co-regulated threshold metric '{}'; use score or abs",
                                metric_raw
                            ),
                            cause_chain: vec![],
                        });
                    }
                };
                let value = value_raw.parse::<f64>().map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Could not parse co-regulated threshold value '{}': {}",
                        value_raw, e
                    ),
                    cause_chain: vec![],
                })?;
                if !value.is_finite() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Co-regulated threshold value must be finite".to_string(),
                        cause_chain: vec![],
                    });
                }
                return Ok(GeneSetCoRegulatedThreshold {
                    metric,
                    operator,
                    value,
                });
            }
        }
        Err(EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Unsupported co-regulated threshold rule '{}'; expected score>=N, score>N, score<=N, score<N, abs>=N, or abs>N",
                raw
            ),
            cause_chain: vec![],
        })
    }

    fn gene_set_co_regulated_parse_direction(
        raw: &str,
    ) -> Result<GeneSetCoRegulatedDirectionRule, EngineError> {
        match raw.trim().to_ascii_lowercase().replace('-', "_").as_str() {
            "both" | "any" | "either" => Ok(GeneSetCoRegulatedDirectionRule::Both),
            "positive" | "up" | "upregulated" | "up_regulated" | "increase" => {
                Ok(GeneSetCoRegulatedDirectionRule::Positive)
            }
            "negative" | "down" | "downregulated" | "down_regulated" | "decrease" => {
                Ok(GeneSetCoRegulatedDirectionRule::Negative)
            }
            other => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported co-regulated direction rule '{}'; expected both, positive, or negative",
                    other
                ),
                cause_chain: vec![],
            }),
        }
    }

    fn gene_set_co_regulated_threshold_passes(
        score: f64,
        threshold: GeneSetCoRegulatedThreshold,
    ) -> bool {
        let value = match threshold.metric {
            GeneSetCoRegulatedThresholdMetric::Score => score,
            GeneSetCoRegulatedThresholdMetric::AbsoluteScore => score.abs(),
        };
        match threshold.operator {
            GeneSetCoRegulatedThresholdOperator::GreaterEqual => value >= threshold.value,
            GeneSetCoRegulatedThresholdOperator::GreaterThan => value > threshold.value,
            GeneSetCoRegulatedThresholdOperator::LessEqual => value <= threshold.value,
            GeneSetCoRegulatedThresholdOperator::LessThan => value < threshold.value,
        }
    }

    fn gene_set_co_regulated_direction_passes(
        score: f64,
        direction: GeneSetCoRegulatedDirectionRule,
    ) -> bool {
        match direction {
            GeneSetCoRegulatedDirectionRule::Both => true,
            GeneSetCoRegulatedDirectionRule::Positive => score > 0.0,
            GeneSetCoRegulatedDirectionRule::Negative => score < 0.0,
        }
    }

    fn gene_set_co_regulated_score_from_value(
        value: &Value,
        scoring_method: Option<&str>,
    ) -> Option<f64> {
        let mut keys = vec![];
        if let Some(method) = scoring_method.map(str::trim).filter(|value| !value.is_empty()) {
            keys.push(method.to_ascii_lowercase().replace('-', "_"));
        }
        keys.extend(
            ["score", "value", "logfc", "log2fc", "statistic", "effect_size"]
                .iter()
                .map(|value| value.to_string()),
        );
        for key in keys {
            if let Some(raw) = value.get(&key) {
                match raw {
                    Value::Number(number) => {
                        if let Some(value) = number.as_f64() {
                            return Some(value);
                        }
                    }
                    Value::String(text) => {
                        if let Ok(value) = text.trim().parse::<f64>() {
                            return Some(value);
                        }
                    }
                    _ => {}
                }
            }
        }
        None
    }

    fn gene_set_co_regulated_requested_matches_value(
        value: &Value,
        root: &Value,
        keys: &[&str],
        requested: &[String],
    ) -> bool {
        if requested.is_empty() {
            return true;
        }
        let mut candidates = Self::gene_set_co_regulated_string_list(value, keys);
        let root_candidates = Self::gene_set_co_regulated_string_list(root, keys);
        Self::gene_set_co_regulated_merge_strings(&mut candidates, &root_candidates);
        requested.iter().any(|requested| {
            candidates.iter().any(|candidate| {
                Self::gene_set_normalize_lookup(candidate)
                    == Self::gene_set_normalize_lookup(requested)
            })
        })
    }

    fn gene_set_co_regulated_selection_from_json(
        path: &str,
        text: &str,
        dataset_ids: &[String],
        contrast_labels: &[String],
        condition_labels: &[String],
        scoring_method: &str,
        threshold: GeneSetCoRegulatedThreshold,
        direction: GeneSetCoRegulatedDirectionRule,
        filters: &[GeneSetProducerFilter],
    ) -> Result<GeneSetCoRegulatedCohortSelection, EngineError> {
        let value: Value = serde_json::from_str(text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Could not parse co-regulated cohort cache '{}': {}",
                path, e
            ),
            cause_chain: vec![],
        })?;
        let mut warnings = vec![];
        let mut root_metadata =
            Self::gene_set_direct_list_metadata_from_value(&value, &mut warnings);
        if value
            .get("schema")
            .and_then(Value::as_str)
            .is_some_and(|schema| schema == GENE_SET_CO_REGULATED_CACHE_SCHEMA)
            && root_metadata.cache_version.is_none()
        {
            root_metadata.cache_version = Some("1".to_string());
        }
        let entries = value
            .get("rows")
            .or_else(|| value.get("records"))
            .or_else(|| value.get("members"))
            .or_else(|| value.get("entries"))
            .and_then(Value::as_array)
            .map(|rows| rows.iter().collect::<Vec<_>>())
            .unwrap_or_else(|| vec![&value]);
        let mut selection = GeneSetCoRegulatedCohortSelection {
            metadata: root_metadata.clone(),
            dataset_ids: Self::gene_set_co_regulated_string_list(
                &value,
                &["dataset_id", "dataset_ids", "dataset"],
            ),
            contrast_labels: Self::gene_set_co_regulated_string_list(
                &value,
                &["contrast", "contrast_label", "contrast_labels"],
            ),
            condition_labels: Self::gene_set_co_regulated_string_list(
                &value,
                &["condition", "condition_label", "condition_labels"],
            ),
            normalization_method: Self::gene_set_direct_list_string(
                &value,
                &["normalization_method", "normalization"],
            ),
            scoring_method: Self::gene_set_direct_list_string(
                &value,
                &["scoring_method", "score_kind", "score"],
            ),
            warnings,
            ..GeneSetCoRegulatedCohortSelection::default()
        };
        Self::gene_set_co_regulated_merge_strings(&mut selection.dataset_ids, dataset_ids);
        Self::gene_set_co_regulated_merge_strings(&mut selection.contrast_labels, contrast_labels);
        Self::gene_set_co_regulated_merge_strings(&mut selection.condition_labels, condition_labels);
        let mut matched_rows = 0usize;

        for entry in entries {
            if !Self::gene_set_co_regulated_requested_matches_value(
                entry,
                &value,
                &["dataset_id", "dataset_ids", "dataset"],
                dataset_ids,
            ) || !Self::gene_set_co_regulated_requested_matches_value(
                entry,
                &value,
                &["contrast", "contrast_label", "contrast_labels"],
                contrast_labels,
            ) || !Self::gene_set_co_regulated_requested_matches_value(
                entry,
                &value,
                &["condition", "condition_label", "condition_labels"],
                condition_labels,
            ) || !Self::gene_set_ontology_value_matches_filters(entry, filters)
            {
                continue;
            }
            matched_rows += 1;
            let Some(score) =
                Self::gene_set_co_regulated_score_from_value(entry, Some(scoring_method))
            else {
                selection
                    .warnings
                    .push("Skipping co-regulated JSON row without numeric score".to_string());
                continue;
            };
            if !Self::gene_set_co_regulated_threshold_passes(score, threshold)
                || !Self::gene_set_co_regulated_direction_passes(score, direction)
            {
                continue;
            }
            if let Some(member) =
                Self::gene_set_direct_list_member_from_value(entry, &mut selection.warnings)
            {
                selection.members.push(member);
            } else {
                selection
                    .warnings
                    .push("Skipping co-regulated JSON row without symbol or gene_id".to_string());
            }
            let mut metadata = Self::gene_set_direct_list_metadata_from_value(
                entry,
                &mut selection.warnings,
            );
            metadata.merge_missing_from(&root_metadata);
            selection.metadata.merge_missing_from(&metadata);
            Self::gene_set_co_regulated_merge_strings(
                &mut selection.dataset_ids,
                &Self::gene_set_co_regulated_string_list(
                    entry,
                    &["dataset_id", "dataset_ids", "dataset"],
                ),
            );
            Self::gene_set_co_regulated_merge_strings(
                &mut selection.contrast_labels,
                &Self::gene_set_co_regulated_string_list(
                    entry,
                    &["contrast", "contrast_label", "contrast_labels"],
                ),
            );
            Self::gene_set_co_regulated_merge_strings(
                &mut selection.condition_labels,
                &Self::gene_set_co_regulated_string_list(
                    entry,
                    &["condition", "condition_label", "condition_labels"],
                ),
            );
        }
        if matched_rows == 0 {
            selection.warnings.push(format!(
                "Co-regulated cohort cache '{}' had no rows matching the requested dataset/contrast/condition filters",
                path
            ));
        } else if selection.members.is_empty() {
            selection.warnings.push(format!(
                "Co-regulated cohort cache '{}' selected no members after threshold and direction filters",
                path
            ));
        }
        Ok(selection)
    }

    fn gene_set_co_regulated_selection_from_tsv(
        path: &str,
        text: &str,
        dataset_ids: &[String],
        contrast_labels: &[String],
        condition_labels: &[String],
        scoring_method: &str,
        threshold: GeneSetCoRegulatedThreshold,
        direction: GeneSetCoRegulatedDirectionRule,
        filters: &[GeneSetProducerFilter],
    ) -> Result<GeneSetCoRegulatedCohortSelection, EngineError> {
        let delimiter = if path.to_ascii_lowercase().ends_with(".csv") {
            b','
        } else {
            b'\t'
        };
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(delimiter)
            .flexible(true)
            .from_reader(text.as_bytes());
        let headers = reader
            .headers()
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read co-regulated cohort cache header '{}': {}",
                    path, e
                ),
                cause_chain: vec![],
            })?
            .iter()
            .map(|header| header.trim().to_ascii_lowercase().replace('-', "_"))
            .collect::<Vec<_>>();
        let find_col = |names: &[&str]| -> Option<usize> {
            names
                .iter()
                .find_map(|name| headers.iter().position(|header| header == *name))
        };
        let dataset_col = find_col(&["dataset_id", "dataset"]);
        let contrast_col = find_col(&["contrast", "contrast_label"]);
        let condition_col = find_col(&["condition", "condition_label"]);
        let symbol_col = find_col(&["symbol", "gene_symbol", "member", "gene"]);
        let gene_id_col = find_col(&["gene_id"]);
        let score_col = headers
            .iter()
            .position(|header| header == &scoring_method.to_ascii_lowercase().replace('-', "_"))
            .or_else(|| find_col(&["score", "value", "logfc", "log2fc", "statistic", "effect_size"]));
        let provider_id_col = find_col(&["provider_id", "provider"]);
        let provider_label_col = find_col(&["provider_label", "provider_name"]);
        let provider_version_col = find_col(&["provider_version"]);
        let cache_id_col = find_col(&["cache_id"]);
        let cache_version_col = find_col(&["cache_version", "cache_format_version"]);
        let cache_digest_col = find_col(&["cache_digest", "digest", "sha256"]);
        let organism_col = find_col(&["organism", "species"]);
        let taxon_id_col = find_col(&["taxon_id", "tax_id"]);
        let namespace_col = find_col(&["symbol_namespace", "namespace"]);
        let normalization_col = find_col(&["normalization_method", "normalization"]);
        let scoring_col = find_col(&["scoring_method", "score_kind"]);
        let get = |record: &csv::StringRecord, col: Option<usize>| -> Option<String> {
            col.and_then(|idx| record.get(idx))
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string)
        };
        let requested_matches = |record: &csv::StringRecord, col: Option<usize>, requested: &[String]| {
            if requested.is_empty() {
                return true;
            }
            let Some(candidate) = get(record, col) else {
                return false;
            };
            requested.iter().any(|requested| {
                Self::gene_set_normalize_lookup(&candidate)
                    == Self::gene_set_normalize_lookup(requested)
            })
        };
        let mut selection = GeneSetCoRegulatedCohortSelection::default();
        Self::gene_set_co_regulated_merge_strings(&mut selection.dataset_ids, dataset_ids);
        Self::gene_set_co_regulated_merge_strings(&mut selection.contrast_labels, contrast_labels);
        Self::gene_set_co_regulated_merge_strings(&mut selection.condition_labels, condition_labels);
        let mut matched_rows = 0usize;

        for row in reader.records() {
            let record = row.map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not read co-regulated cohort cache row '{}': {}",
                    path, e
                ),
                cause_chain: vec![],
            })?;
            if !requested_matches(&record, dataset_col, dataset_ids)
                || !requested_matches(&record, contrast_col, contrast_labels)
                || !requested_matches(&record, condition_col, condition_labels)
                || !Self::gene_set_ontology_record_matches_filters(&headers, &record, filters)
            {
                continue;
            }
            matched_rows += 1;
            let Some(score) = get(&record, score_col).and_then(|raw| raw.parse::<f64>().ok())
            else {
                selection
                    .warnings
                    .push("Skipping co-regulated TSV row without numeric score".to_string());
                continue;
            };
            if !score.is_finite()
                || !Self::gene_set_co_regulated_threshold_passes(score, threshold)
                || !Self::gene_set_co_regulated_direction_passes(score, direction)
            {
                continue;
            }
            if let Some(member) = get(&record, symbol_col).or_else(|| get(&record, gene_id_col)) {
                selection.members.push(member);
            } else {
                selection
                    .warnings
                    .push("Skipping co-regulated TSV row without symbol or gene_id".to_string());
            }
            if selection.metadata.provider_id.is_none() {
                selection.metadata.provider_id = get(&record, provider_id_col);
            }
            if selection.metadata.provider_label.is_none() {
                selection.metadata.provider_label = get(&record, provider_label_col);
            }
            if selection.metadata.provider_version.is_none() {
                selection.metadata.provider_version = get(&record, provider_version_col);
            }
            if selection.metadata.cache_id.is_none() {
                selection.metadata.cache_id = get(&record, cache_id_col);
            }
            if selection.metadata.cache_version.is_none() {
                selection.metadata.cache_version = get(&record, cache_version_col);
            }
            if selection.metadata.cache_digest.is_none() {
                selection.metadata.cache_digest = get(&record, cache_digest_col);
            }
            if selection.metadata.organism.is_none() {
                selection.metadata.organism = get(&record, organism_col);
            }
            if selection.metadata.taxon_id.is_none() {
                selection.metadata.taxon_id = get(&record, taxon_id_col);
            }
            if selection.metadata.symbol_namespace.is_none() {
                selection.metadata.symbol_namespace = get(&record, namespace_col);
            }
            if selection.normalization_method.is_none() {
                selection.normalization_method = get(&record, normalization_col);
            }
            if selection.scoring_method.is_none() {
                selection.scoring_method = get(&record, scoring_col);
            }
            if let Some(value) = get(&record, dataset_col) {
                Self::gene_set_co_regulated_merge_strings(&mut selection.dataset_ids, &[value]);
            }
            if let Some(value) = get(&record, contrast_col) {
                Self::gene_set_co_regulated_merge_strings(&mut selection.contrast_labels, &[value]);
            }
            if let Some(value) = get(&record, condition_col) {
                Self::gene_set_co_regulated_merge_strings(&mut selection.condition_labels, &[value]);
            }
        }
        if matched_rows == 0 {
            selection.warnings.push(format!(
                "Co-regulated cohort cache '{}' had no rows matching the requested dataset/contrast/condition filters",
                path
            ));
        } else if selection.members.is_empty() {
            selection.warnings.push(format!(
                "Co-regulated cohort cache '{}' selected no members after threshold and direction filters",
                path
            ));
        }
        Ok(selection)
    }

    fn gene_set_co_regulated_select_cache(
        cache_path: &str,
        dataset_ids: &[String],
        contrast_labels: &[String],
        condition_labels: &[String],
        scoring_method: &str,
        threshold_rule: &str,
        sign_direction_rule: &str,
        filters: &[GeneSetProducerFilter],
    ) -> Result<GeneSetCoRegulatedCohortSelection, EngineError> {
        let threshold = Self::gene_set_co_regulated_parse_threshold(threshold_rule)?;
        let direction = Self::gene_set_co_regulated_parse_direction(sign_direction_rule)?;
        let text = std::fs::read_to_string(cache_path).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "Could not read co-regulated cohort cache '{}': {}",
                cache_path, e
            ),
            cause_chain: vec![],
        })?;
        let trimmed = text.trim_start();
        if trimmed.starts_with('{') || trimmed.starts_with('[') {
            Self::gene_set_co_regulated_selection_from_json(
                cache_path,
                &text,
                dataset_ids,
                contrast_labels,
                condition_labels,
                scoring_method,
                threshold,
                direction,
                filters,
            )
        } else {
            Self::gene_set_co_regulated_selection_from_tsv(
                cache_path,
                &text,
                dataset_ids,
                contrast_labels,
                condition_labels,
                scoring_method,
                threshold,
                direction,
                filters,
            )
        }
    }

    pub(crate) fn produce_gene_set_direct_list(
        &self,
        cache_path: &str,
        query: Option<&str>,
        genome_id: Option<&str>,
        gene_group_catalog_path: Option<&str>,
        genome_catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        provider_id: Option<&str>,
        provider_label: Option<&str>,
        provider_version: Option<&str>,
        cache_id: Option<&str>,
        cache_version: Option<&str>,
        cache_digest: Option<&str>,
        organism: Option<&str>,
        taxon_id: Option<&str>,
        symbol_namespace: Option<&str>,
        review_status: Option<GeneSetResolutionReviewStatus>,
        filters: &[GeneSetProducerFilter],
    ) -> Result<GeneSetResolutionReport, EngineError> {
        let mut selection = Self::gene_set_direct_list_select_cache(cache_path, query)?;
        selection.metadata.apply_overrides(
            provider_id,
            provider_label,
            provider_version,
            cache_id,
            cache_version,
            cache_digest,
            organism,
            taxon_id,
            symbol_namespace,
            review_status,
            filters,
        );
        if selection.metadata.cache_id.is_none() {
            selection.metadata.cache_id =
                Some(Self::gene_set_direct_list_fallback_cache_id(cache_path));
        }
        Self::gene_set_direct_list_validate_cache_version(
            selection.metadata.cache_version.as_deref(),
        )?;
        let provider_id = selection
            .metadata
            .provider_id
            .clone()
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: "Direct gene-list producer requires provider_id in cache metadata or --provider-id"
                    .to_string(),
                cause_chain: vec![],
            })?;
        if selection.metadata.provider_version.is_none()
            && selection.metadata.cache_version.is_none()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Direct gene-list producer requires provider_version or cache_version metadata"
                        .to_string(),
                cause_chain: vec![],
            });
        }
        if selection.metadata.organism.is_none()
            && selection.metadata.taxon_id.is_none()
            && selection.metadata.symbol_namespace.is_none()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Direct gene-list producer requires organism, taxon_id, or symbol_namespace metadata"
                        .to_string(),
                cause_chain: vec![],
            });
        }
        if selection.members.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Direct gene-list cache '{}' selected no candidate members",
                    cache_path
                ),
                cause_chain: vec![],
            });
        }

        let mut report = self.resolve_gene_set(
            GeneSetRequest::ExplicitMembers {
                members: selection.members.clone(),
            },
            genome_id,
            gene_group_catalog_path,
            genome_catalog_path,
            cache_dir,
            false,
            false,
        )?;
        report.review_status = selection
            .metadata
            .review_status
            .unwrap_or(GeneSetResolutionReviewStatus::Unreviewed);
        report.organism = selection.metadata.organism.clone();
        report.taxon_id = selection.metadata.taxon_id.clone();
        report.symbol_namespace = selection.metadata.symbol_namespace.clone();
        let cache_id = selection
            .metadata
            .cache_id
            .clone()
            .unwrap_or_else(|| Self::gene_set_direct_list_fallback_cache_id(cache_path));
        report.producer = Some(GeneSetProducerProvenance {
            producer_kind: GeneSetProducerKind::DirectGeneList,
            provider_id,
            provider_label: selection.metadata.provider_label.clone(),
            provider_version: selection.metadata.provider_version.clone(),
            cache_id: Some(cache_id.clone()),
            cache_path: Some(cache_path.to_string()),
            cache_version: selection.metadata.cache_version.clone(),
            cache_digest: selection.metadata.cache_digest.clone(),
            import_op_id: None,
            imported_at_unix_ms: None,
        });
        report.query_metadata = Some(GeneSetProducerQueryMetadata {
            query_kind: "direct_gene_list".to_string(),
            query_id: selection
                .list_id
                .clone()
                .or_else(|| query.map(str::to_string)),
            query_label: selection.list_label.clone(),
            organism: selection.metadata.organism.clone(),
            taxon_id: selection.metadata.taxon_id.clone(),
            symbol_namespace: selection.metadata.symbol_namespace.clone(),
            filters: selection.metadata.filters.clone(),
        });
        let provenance = GeneSetProvenanceRow {
            source_kind: "direct_gene_list".to_string(),
            source_id: selection
                .list_id
                .clone()
                .unwrap_or_else(|| cache_id.clone()),
            source_label: selection.list_label.clone(),
            source_path: Some(cache_path.to_string()),
            note: Some(format!(
                "Direct gene-list producer selected {} candidate member(s)",
                selection.members.len()
            )),
        };
        if !report.provenance.contains(&provenance) {
            report.provenance.push(provenance.clone());
        }
        for member in &mut report.resolved_members {
            if !member.provenance.contains(&provenance) {
                member.provenance.push(provenance.clone());
            }
        }
        for unresolved in &mut report.unresolved_members {
            if unresolved.source_kind == "explicit_members" {
                unresolved.source_kind = "direct_gene_list".to_string();
                unresolved.source_id = selection.list_id.clone().or_else(|| Some(cache_id.clone()));
            }
        }
        report.warnings.extend(selection.warnings);
        Ok(report)
    }

    pub(crate) fn produce_gene_set_ontology_assignment(
        &self,
        cache_path: &str,
        term: &str,
        ontology_namespace: Option<&str>,
        genome_id: Option<&str>,
        gene_group_catalog_path: Option<&str>,
        genome_catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        provider_id: Option<&str>,
        provider_label: Option<&str>,
        provider_version: Option<&str>,
        cache_id: Option<&str>,
        cache_version: Option<&str>,
        cache_digest: Option<&str>,
        organism: Option<&str>,
        taxon_id: Option<&str>,
        symbol_namespace: Option<&str>,
        review_status: Option<GeneSetResolutionReviewStatus>,
        filters: &[GeneSetProducerFilter],
    ) -> Result<GeneSetResolutionReport, EngineError> {
        let mut selection = Self::gene_set_ontology_assignment_select_cache(
            cache_path,
            term,
            ontology_namespace,
            filters,
        )?;
        selection.metadata.apply_overrides(
            provider_id,
            provider_label,
            provider_version,
            cache_id,
            cache_version,
            cache_digest,
            organism,
            taxon_id,
            symbol_namespace,
            review_status,
            filters,
        );
        if selection.metadata.cache_id.is_none() {
            selection.metadata.cache_id = Some(
                Self::gene_set_ontology_assignment_fallback_cache_id(cache_path),
            );
        }
        Self::gene_set_ontology_assignment_validate_cache_version(
            selection.metadata.cache_version.as_deref(),
        )?;
        let provider_id = selection
            .metadata
            .provider_id
            .clone()
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Ontology assignment producer requires provider_id in cache metadata or --provider-id"
                        .to_string(),
                cause_chain: vec![],
            })?;
        if selection.metadata.provider_version.is_none()
            && selection.metadata.cache_version.is_none()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Ontology assignment producer requires provider_version or cache_version metadata"
                        .to_string(),
                cause_chain: vec![],
            });
        }
        if selection.metadata.organism.is_none()
            && selection.metadata.taxon_id.is_none()
            && selection.metadata.symbol_namespace.is_none()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Ontology assignment producer requires organism, taxon_id, or symbol_namespace metadata"
                        .to_string(),
                cause_chain: vec![],
            });
        }

        let mut report = self.resolve_gene_set(
            GeneSetRequest::ExplicitMembers {
                members: selection.members.clone(),
            },
            genome_id,
            gene_group_catalog_path,
            genome_catalog_path,
            cache_dir,
            false,
            false,
        )?;
        report.review_status = selection
            .metadata
            .review_status
            .unwrap_or(GeneSetResolutionReviewStatus::Unreviewed);
        report.organism = selection.metadata.organism.clone();
        report.taxon_id = selection.metadata.taxon_id.clone();
        report.symbol_namespace = selection.metadata.symbol_namespace.clone();
        let cache_id = selection
            .metadata
            .cache_id
            .clone()
            .unwrap_or_else(|| Self::gene_set_ontology_assignment_fallback_cache_id(cache_path));
        let term_source_id = if selection.term_id == term {
            ontology_namespace
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .filter(|_| !term.contains(':'))
                .map(|namespace| format!("{}:{}", namespace.to_ascii_uppercase(), term))
                .unwrap_or_else(|| term.to_string())
        } else {
            selection.term_id.clone()
        };
        report.producer = Some(GeneSetProducerProvenance {
            producer_kind: GeneSetProducerKind::OntologyAssignment,
            provider_id,
            provider_label: selection.metadata.provider_label.clone(),
            provider_version: selection.metadata.provider_version.clone(),
            cache_id: Some(cache_id.clone()),
            cache_path: Some(cache_path.to_string()),
            cache_version: selection.metadata.cache_version.clone(),
            cache_digest: selection.metadata.cache_digest.clone(),
            import_op_id: None,
            imported_at_unix_ms: None,
        });
        report.query_metadata = Some(GeneSetProducerQueryMetadata {
            query_kind: "ontology_assignment".to_string(),
            query_id: Some(term_source_id.clone()),
            query_label: selection.term_label.clone(),
            organism: selection.metadata.organism.clone(),
            taxon_id: selection.metadata.taxon_id.clone(),
            symbol_namespace: selection.metadata.symbol_namespace.clone(),
            filters: selection.metadata.filters.clone(),
        });
        let provenance = GeneSetProvenanceRow {
            source_kind: "ontology_assignment".to_string(),
            source_id: term_source_id.clone(),
            source_label: selection.term_label.clone(),
            source_path: Some(cache_path.to_string()),
            note: Some(format!(
                "Ontology assignment producer selected {} candidate member(s)",
                selection.members.len()
            )),
        };
        if !report.provenance.contains(&provenance) {
            report.provenance.push(provenance.clone());
        }
        for member in &mut report.resolved_members {
            if !member.provenance.contains(&provenance) {
                member.provenance.push(provenance.clone());
            }
        }
        for unresolved in &mut report.unresolved_members {
            if unresolved.source_kind == "explicit_members" {
                unresolved.source_kind = "ontology_assignment".to_string();
                unresolved.source_id = Some(term_source_id.clone());
            }
        }
        if selection.members.is_empty() {
            report.unresolved_members.push(GeneSetUnresolvedMember {
                query: term_source_id.clone(),
                reason: "no ontology assignment rows matched the requested term and filters"
                    .to_string(),
                source_kind: "ontology_assignment".to_string(),
                source_id: Some(term_source_id),
            });
        }
        report.resolved_member_count = report.resolved_members.len();
        report.unresolved_member_count = report.unresolved_members.len();
        report.warnings.extend(selection.warnings);
        Ok(report)
    }

    pub(crate) fn produce_gene_set_co_regulated_cohort(
        &self,
        cache_path: &str,
        dataset_ids: &[String],
        contrast_labels: &[String],
        condition_labels: &[String],
        normalization_method: Option<&str>,
        scoring_method: &str,
        threshold_rule: &str,
        sign_direction_rule: &str,
        relationship: GeneSetCohortRelationship,
        genome_id: Option<&str>,
        gene_group_catalog_path: Option<&str>,
        genome_catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        provider_id: Option<&str>,
        provider_label: Option<&str>,
        provider_version: Option<&str>,
        cache_id: Option<&str>,
        cache_version: Option<&str>,
        cache_digest: Option<&str>,
        organism: Option<&str>,
        taxon_id: Option<&str>,
        symbol_namespace: Option<&str>,
        review_status: Option<GeneSetResolutionReviewStatus>,
        filters: &[GeneSetProducerFilter],
    ) -> Result<GeneSetResolutionReport, EngineError> {
        if dataset_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Co-regulated cohort producer requires at least one dataset id"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        if contrast_labels.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Co-regulated cohort producer requires at least one contrast label"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        if scoring_method.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Co-regulated cohort producer requires a scoring method".to_string(),
                cause_chain: vec![],
            });
        }
        let mut selection = Self::gene_set_co_regulated_select_cache(
            cache_path,
            dataset_ids,
            contrast_labels,
            condition_labels,
            scoring_method,
            threshold_rule,
            sign_direction_rule,
            filters,
        )?;
        selection.metadata.apply_overrides(
            provider_id,
            provider_label,
            provider_version,
            cache_id,
            cache_version,
            cache_digest,
            organism,
            taxon_id,
            symbol_namespace,
            review_status,
            filters,
        );
        if selection.metadata.cache_id.is_none() {
            selection.metadata.cache_id =
                Some(Self::gene_set_co_regulated_fallback_cache_id(cache_path));
        }
        Self::gene_set_co_regulated_validate_cache_version(
            selection.metadata.cache_version.as_deref(),
        )?;
        let provider_id = selection
            .metadata
            .provider_id
            .clone()
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Co-regulated cohort producer requires provider_id in cache metadata or --provider-id"
                        .to_string(),
                cause_chain: vec![],
            })?;
        if selection.metadata.provider_version.is_none()
            && selection.metadata.cache_version.is_none()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Co-regulated cohort producer requires provider_version or cache_version metadata"
                        .to_string(),
                cause_chain: vec![],
            });
        }
        if selection.metadata.organism.is_none()
            && selection.metadata.taxon_id.is_none()
            && selection.metadata.symbol_namespace.is_none()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Co-regulated cohort producer requires organism, taxon_id, or symbol_namespace metadata"
                        .to_string(),
                cause_chain: vec![],
            });
        }

        let mut report = self.resolve_gene_set(
            GeneSetRequest::ExplicitMembers {
                members: selection.members.clone(),
            },
            genome_id,
            gene_group_catalog_path,
            genome_catalog_path,
            cache_dir,
            false,
            false,
        )?;
        report.review_status = selection
            .metadata
            .review_status
            .unwrap_or(GeneSetResolutionReviewStatus::Unreviewed);
        report.organism = selection.metadata.organism.clone();
        report.taxon_id = selection.metadata.taxon_id.clone();
        report.symbol_namespace = selection.metadata.symbol_namespace.clone();
        let cache_id = selection
            .metadata
            .cache_id
            .clone()
            .unwrap_or_else(|| Self::gene_set_co_regulated_fallback_cache_id(cache_path));
        let query_id = format!(
            "{}:{}",
            selection.dataset_ids.join(","),
            selection.contrast_labels.join(",")
        );
        report.producer = Some(GeneSetProducerProvenance {
            producer_kind: GeneSetProducerKind::CoRegulatedCohort,
            provider_id,
            provider_label: selection.metadata.provider_label.clone(),
            provider_version: selection.metadata.provider_version.clone(),
            cache_id: Some(cache_id),
            cache_path: Some(cache_path.to_string()),
            cache_version: selection.metadata.cache_version.clone(),
            cache_digest: selection.metadata.cache_digest.clone(),
            import_op_id: None,
            imported_at_unix_ms: None,
        });
        report.query_metadata = Some(GeneSetProducerQueryMetadata {
            query_kind: "co_regulated_cohort".to_string(),
            query_id: Some(query_id.clone()),
            query_label: Some(format!(
                "co-regulated cohort {}",
                selection.contrast_labels.join(",")
            )),
            organism: selection.metadata.organism.clone(),
            taxon_id: selection.metadata.taxon_id.clone(),
            symbol_namespace: selection.metadata.symbol_namespace.clone(),
            filters: selection.metadata.filters.clone(),
        });
        let normalization_method = normalization_method
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .or(selection.normalization_method.clone())
            .unwrap_or_else(|| "unspecified".to_string());
        let scoring_method = selection
            .scoring_method
            .clone()
            .unwrap_or_else(|| scoring_method.to_string());
        report.co_regulated_metadata = Some(GeneSetCoRegulatedProducerMetadata {
            dataset_ids: selection.dataset_ids.clone(),
            contrast_labels: selection.contrast_labels.clone(),
            condition_labels: selection.condition_labels.clone(),
            normalization_method,
            scoring_method,
            threshold_rule: threshold_rule.to_string(),
            sign_direction_rule: sign_direction_rule.to_string(),
            relationship,
            ..GeneSetCoRegulatedProducerMetadata::default()
        });
        let provenance = GeneSetProvenanceRow {
            source_kind: "co_regulated_cohort".to_string(),
            source_id: query_id,
            source_label: Some("co-regulated cohort producer".to_string()),
            source_path: Some(cache_path.to_string()),
            note: Some(format!(
                "Co-regulated cohort producer selected {} candidate member(s); this retrieval result does not prove regulation",
                selection.members.len()
            )),
        };
        if !report.provenance.contains(&provenance) {
            report.provenance.push(provenance.clone());
        }
        for member in &mut report.resolved_members {
            if !member.provenance.contains(&provenance) {
                member.provenance.push(provenance.clone());
            }
        }
        for unresolved in &mut report.unresolved_members {
            if unresolved.source_kind == "explicit_members" {
                unresolved.source_kind = "co_regulated_cohort".to_string();
                unresolved.source_id = report
                    .query_metadata
                    .as_ref()
                    .and_then(|metadata| metadata.query_id.clone());
            }
        }
        report.warnings.extend(selection.warnings);
        Ok(report)
    }

    pub(crate) fn resolve_gene_set(
        &self,
        source: GeneSetRequest,
        genome_id: Option<&str>,
        gene_group_catalog_path: Option<&str>,
        genome_catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        allow_draft: bool,
        allow_deprecated: bool,
    ) -> Result<GeneSetResolutionReport, EngineError> {
        let (resolved_genome_id, genome_catalog_label, mut genes) =
            Self::gene_set_load_genome_genes(genome_id, genome_catalog_path, cache_dir)?;
        Self::gene_set_sort_indexed_genes(&mut genes);
        let mut warnings = vec![];
        let mut unresolved_members = vec![];
        let mut contributing_group_ids = vec![];
        let mut provenance = vec![];
        let mut random = None;
        let mut gene_group_catalog_label = None;
        let mut requested_member_count = 0usize;

        let mut resolved_members = match &source {
            GeneSetRequest::CatalogGroup { query } => {
                let catalog = crate::gene_groups::load_gene_group_catalog(gene_group_catalog_path)
                    .map_err(|message| EngineError {
                        code: ErrorCode::InvalidInput,
                        message,
                        cause_chain: vec![],
                    })?;
                gene_group_catalog_label = Some(catalog.catalog_label().to_string());
                let groups = catalog.resolve_exact(query);
                if groups.is_empty() {
                    unresolved_members.push(GeneSetUnresolvedMember {
                        query: query.clone(),
                        reason: "no gene-group catalog entry matched the requested token"
                            .to_string(),
                        source_kind: "catalog_group".to_string(),
                        source_id: None,
                    });
                    vec![]
                } else if groups.len() > 1 {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Gene-set catalog query '{}' matched {} groups; use a stable group id",
                            query,
                            groups.len()
                        ),
                        cause_chain: vec![],
                    });
                } else {
                    let group = groups[0].clone();
                    Self::gene_set_group_allowed(
                        &group,
                        allow_draft,
                        allow_deprecated,
                        &mut warnings,
                    )?;
                    requested_member_count = group.record.members.len();
                    contributing_group_ids.push(group.record.id.clone());
                    provenance.push(GeneSetProvenanceRow {
                        source_kind: "catalog_group".to_string(),
                        source_id: group.record.id.clone(),
                        source_label: Some(group.record.label.clone()),
                        source_path: Some(group.source_path.clone()),
                        note: group.record.provenance.clone(),
                    });
                    Self::gene_set_resolve_group_members(
                        &[group],
                        &genes,
                        &mut warnings,
                        &mut unresolved_members,
                    )
                }
            }
            GeneSetRequest::ExplicitMembers { members } => {
                requested_member_count = members
                    .iter()
                    .filter(|value| !value.trim().is_empty())
                    .count();
                Self::gene_set_resolve_explicit_members(members, &genes, &mut unresolved_members)
            }
            GeneSetRequest::ExternalMapping { namespace, id } => {
                let catalog = crate::gene_groups::load_gene_group_catalog(gene_group_catalog_path)
                    .map_err(|message| EngineError {
                        code: ErrorCode::InvalidInput,
                        message,
                        cause_chain: vec![],
                    })?;
                gene_group_catalog_label = Some(catalog.catalog_label().to_string());
                let mut groups = catalog.groups_by_external_mapping(namespace, id);
                if groups.is_empty() {
                    warnings.push(format!(
                        "No local gene groups map to external mapping {}:{}",
                        namespace, id
                    ));
                    unresolved_members.push(GeneSetUnresolvedMember {
                        query: format!("{}:{}", namespace, id),
                        reason: "no local gene group maps to the requested external mapping"
                            .to_string(),
                        source_kind: "external_mapping".to_string(),
                        source_id: Some(format!("{}:{}", namespace, id)),
                    });
                    vec![]
                } else {
                    for group in &groups {
                        Self::gene_set_group_allowed(
                            group,
                            allow_draft,
                            allow_deprecated,
                            &mut warnings,
                        )?;
                    }
                    groups.sort_by(|a, b| a.record.id.cmp(&b.record.id));
                    for group in &groups {
                        requested_member_count += group.record.members.len();
                        contributing_group_ids.push(group.record.id.clone());
                        provenance.push(GeneSetProvenanceRow {
                            source_kind: "external_mapping".to_string(),
                            source_id: format!("{}:{}", namespace, id),
                            source_label: Some(group.record.label.clone()),
                            source_path: Some(group.source_path.clone()),
                            note: Some(format!("Contributing group {}", group.record.id)),
                        });
                    }
                    Self::gene_set_resolve_group_members(
                        &groups,
                        &genes,
                        &mut warnings,
                        &mut unresolved_members,
                    )
                }
            }
            GeneSetRequest::GenomicNeighbors {
                anchor,
                flank_gene_count,
                flank_bp,
                exclude_anchor,
            } => {
                if resolved_genome_id.is_none() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Genomic-neighbor gene sets require --genome GENOME_ID"
                            .to_string(),
                        cause_chain: vec![],
                    });
                }
                let members = Self::gene_set_resolve_neighbors(
                    anchor,
                    *flank_gene_count,
                    *flank_bp,
                    *exclude_anchor,
                    &genes,
                    &mut warnings,
                    &mut unresolved_members,
                )?;
                requested_member_count = members.len();
                members
            }
            GeneSetRequest::Random {
                count,
                random_seed,
                exclude_members,
            } => {
                let Some(genome_id) = resolved_genome_id.as_deref() else {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "Random gene sets require --genome GENOME_ID".to_string(),
                        cause_chain: vec![],
                    });
                };
                requested_member_count = *count;
                let gene_index_source = format!(
                    "{}|{}|{}",
                    genome_id,
                    genome_catalog_label.as_deref().unwrap_or("unknown_catalog"),
                    cache_dir.unwrap_or("")
                );
                let (members, random_provenance) = Self::gene_set_resolve_random(
                    genome_id,
                    *count,
                    *random_seed,
                    exclude_members,
                    &genes,
                    &gene_index_source,
                    &mut warnings,
                );
                random = Some(random_provenance);
                members
            }
        };

        Self::gene_set_sort_members(&mut resolved_members);
        Ok(GeneSetResolutionReport {
            schema: GENE_SET_RESOLUTION_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            request: source,
            review_status: GeneSetResolutionReviewStatus::Unreviewed,
            genome_id: resolved_genome_id,
            organism: None,
            taxon_id: None,
            symbol_namespace: None,
            gene_group_catalog_label,
            genome_catalog_label,
            producer: None,
            query_metadata: None,
            co_regulated_metadata: None,
            contributing_group_ids,
            random,
            requested_member_count,
            resolved_member_count: resolved_members.len(),
            unresolved_member_count: unresolved_members.len(),
            resolved_members,
            unresolved_members,
            warnings,
            provenance,
        })
    }

    pub(crate) fn build_gene_set_promoter_cohort(
        &self,
        genome_id: &str,
        resolution: GeneSetResolutionReport,
        relationship: GeneSetCohortRelationship,
        upstream_bp: usize,
        downstream_bp: usize,
        genome_catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<GeneSetPromoterCohortReport, EngineError> {
        let effective_catalog_path =
            genome_catalog_path.unwrap_or(crate::genomes::default_catalog_discovery_token(false));
        let (catalog, _) = Self::open_reference_genome_catalog(Some(effective_catalog_path))?;
        let effective_relationship = if relationship == GeneSetCohortRelationship::Unspecified {
            resolution
                .co_regulated_metadata
                .as_ref()
                .map(|metadata| metadata.relationship)
                .filter(|relationship| *relationship != GeneSetCohortRelationship::Unspecified)
                .unwrap_or(relationship)
        } else {
            relationship
        };
        let mut warnings = resolution.warnings.clone();
        let mut unresolved_members = resolution.unresolved_members.clone();
        let mut windows = vec![];
        for member in &resolution.resolved_members {
            let gene_query = member
                .gene_id
                .as_deref()
                .filter(|value| !value.trim().is_empty())
                .unwrap_or(&member.symbol);
            let resolved = match Self::resolve_genome_promoter_slice_request(
                &catalog,
                genome_id,
                gene_query,
                None,
                None,
                upstream_bp,
                downstream_bp,
                cache_dir,
            ) {
                Ok(resolved) => resolved,
                Err(err) => {
                    unresolved_members.push(GeneSetUnresolvedMember {
                        query: gene_query.to_string(),
                        reason: format!("promoter window could not be resolved: {}", err.message),
                        source_kind: "promoter_cohort".to_string(),
                        source_id: Some(member.dedup_key.clone()),
                    });
                    continue;
                }
            };
            warnings.extend(resolved.warnings.clone());
            let strand = resolved
                .selected_gene
                .strand
                .or(resolved.selected_transcript.strand)
                .unwrap_or('+');
            let promoter_length_bp = resolved
                .extract_end_1based
                .saturating_sub(resolved.extract_start_1based)
                .saturating_add(1);
            let display_label = resolved
                .selected_gene
                .gene_name
                .clone()
                .or_else(|| resolved.selected_gene.gene_id.clone())
                .unwrap_or_else(|| member.symbol.clone());
            windows.push(GeneSetPromoterWindow {
                member_dedup_key: member.dedup_key.clone(),
                symbol: member.symbol.clone(),
                gene_id: member
                    .gene_id
                    .clone()
                    .or(resolved.selected_gene.gene_id.clone()),
                gene_query: resolved.query.clone(),
                occurrence: resolved.occurrence,
                transcript_id_requested: None,
                transcript_id: resolved.selected_transcript.transcript_id.clone(),
                display_label,
                chromosome: resolved.selected_transcript.chromosome.clone(),
                strand: strand.to_string(),
                promoter_start_1based: resolved.extract_start_1based,
                promoter_end_1based: resolved.extract_end_1based,
                promoter_length_bp,
                tss_1based: resolved.tss_1based,
                sequence_orientation: "transcription_aligned".to_string(),
                used_fuzzy_gene_match: resolved.used_fuzzy_gene_match,
            });
        }
        Ok(GeneSetPromoterCohortReport {
            schema: GENE_SET_PROMOTER_COHORT_SCHEMA.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            genome_id: genome_id.to_string(),
            upstream_bp,
            downstream_bp,
            relationship: effective_relationship,
            requested_member_count: resolution.resolved_members.len(),
            returned_window_count: windows.len(),
            gene_set_resolution: resolution,
            windows,
            relationship_flags: vec![],
            unresolved_members,
            warnings,
        })
    }

    fn gene_set_artifact_id_from_report_parts(
        prefix: &str,
        op_id: Option<&str>,
        run_id: Option<&str>,
        generated_at_unix_ms: u128,
    ) -> String {
        op_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .or_else(|| {
                run_id
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(|run_id| format!("{prefix}:{run_id}:{generated_at_unix_ms}"))
            })
            .unwrap_or_else(|| format!("{prefix}:{generated_at_unix_ms}"))
    }

    pub(crate) fn gene_set_resolution_artifact_id(report: &GeneSetResolutionReport) -> String {
        Self::gene_set_artifact_id_from_report_parts(
            "resolution",
            report.op_id.as_deref(),
            report.run_id.as_deref(),
            report.generated_at_unix_ms,
        )
    }

    pub(crate) fn gene_set_promoter_cohort_artifact_id(
        report: &GeneSetPromoterCohortReport,
    ) -> String {
        Self::gene_set_artifact_id_from_report_parts(
            "promoter_cohort",
            report.op_id.as_deref(),
            report.run_id.as_deref(),
            report.generated_at_unix_ms,
        )
    }

    pub(crate) fn gene_set_cutrun_support_artifact_id(
        report: &GeneSetCutRunRegulatorySupportReport,
    ) -> String {
        Self::gene_set_artifact_id_from_report_parts(
            "cutrun_support",
            report.op_id.as_deref(),
            report.run_id.as_deref(),
            report.generated_at_unix_ms,
        )
    }

    fn read_gene_set_artifact_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> GeneSetArtifactStore {
        let mut store = value
            .cloned()
            .and_then(|value| serde_json::from_value::<GeneSetArtifactStore>(value).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = GENE_SET_ARTIFACTS_SCHEMA.to_string();
        }
        store
    }

    fn read_gene_set_artifact_store(&self) -> GeneSetArtifactStore {
        Self::read_gene_set_artifact_store_from_metadata(
            self.state.metadata.get(GENE_SET_ARTIFACTS_METADATA_KEY),
        )
    }

    fn write_gene_set_artifact_store(
        &mut self,
        mut store: GeneSetArtifactStore,
    ) -> Result<(), EngineError> {
        if store.resolutions.is_empty()
            && store.promoter_cohorts.is_empty()
            && store.cutrun_support_reports.is_empty()
        {
            self.state.metadata.remove(GENE_SET_ARTIFACTS_METADATA_KEY);
            return Ok(());
        }
        store.schema = GENE_SET_ARTIFACTS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize gene-set artifact metadata: {e}"),
            cause_chain: vec![],
        })?;
        self.state
            .metadata
            .insert(GENE_SET_ARTIFACTS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub(crate) fn upsert_gene_set_resolution_artifact(
        &mut self,
        report: GeneSetResolutionReport,
    ) -> Result<(), EngineError> {
        let mut store = self.read_gene_set_artifact_store();
        store
            .resolutions
            .insert(Self::gene_set_resolution_artifact_id(&report), report);
        self.write_gene_set_artifact_store(store)
    }

    pub(crate) fn upsert_gene_set_promoter_cohort_artifact(
        &mut self,
        report: GeneSetPromoterCohortReport,
    ) -> Result<(), EngineError> {
        let mut store = self.read_gene_set_artifact_store();
        store
            .promoter_cohorts
            .insert(Self::gene_set_promoter_cohort_artifact_id(&report), report);
        self.write_gene_set_artifact_store(store)
    }

    pub(crate) fn upsert_gene_set_cutrun_support_artifact(
        &mut self,
        report: GeneSetCutRunRegulatorySupportReport,
    ) -> Result<(), EngineError> {
        let mut store = self.read_gene_set_artifact_store();
        store.cutrun_support_reports.insert(
            Self::gene_set_cutrun_support_artifact_id(&report),
            report,
        );
        self.write_gene_set_artifact_store(store)
    }

    pub(crate) fn gene_set_resolution_artifacts_from_state(
        state: &ProjectState,
    ) -> Vec<GeneSetResolutionReport> {
        Self::read_gene_set_artifact_store_from_metadata(
            state.metadata.get(GENE_SET_ARTIFACTS_METADATA_KEY),
        )
        .resolutions
        .into_values()
        .collect()
    }

    pub(crate) fn gene_set_promoter_cohort_artifacts_from_state(
        state: &ProjectState,
    ) -> Vec<GeneSetPromoterCohortReport> {
        Self::read_gene_set_artifact_store_from_metadata(
            state.metadata.get(GENE_SET_ARTIFACTS_METADATA_KEY),
        )
        .promoter_cohorts
        .into_values()
        .collect()
    }

    pub(crate) fn gene_set_cutrun_support_artifacts_from_state(
        state: &ProjectState,
    ) -> Vec<GeneSetCutRunRegulatorySupportReport> {
        Self::read_gene_set_artifact_store_from_metadata(
            state.metadata.get(GENE_SET_ARTIFACTS_METADATA_KEY),
        )
        .cutrun_support_reports
        .into_values()
        .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gene_groups::LoadedGeneGroupRecord;
    use gentle_protocol::{GeneGroupMember, GeneGroupRecord};

    fn gene(symbol: &str, gene_id: &str, chromosome: &str, start_1based: usize) -> IndexedGeneRow {
        GentleEngine::gene_set_indexed_row(GenomeGeneRecord {
            chromosome: chromosome.to_string(),
            start_1based,
            end_1based: start_1based + 49,
            strand: Some('+'),
            gene_id: Some(gene_id.to_string()),
            gene_name: Some(symbol.to_string()),
            biotype: Some("protein_coding".to_string()),
        })
    }

    fn member(symbol: &str, status: Option<&str>) -> GeneGroupMember {
        GeneGroupMember {
            symbol: symbol.to_string(),
            gene_id: Some(format!("{symbol}_ID")),
            aliases: vec![format!("{symbol}_alias")],
            status: status.map(str::to_string),
            confidence: Some("reviewed".to_string()),
            provenance: Some("synthetic test member".to_string()),
            ..GeneGroupMember::default()
        }
    }

    fn loaded_group(id: &str, curation_status: &str) -> LoadedGeneGroupRecord {
        LoadedGeneGroupRecord {
            record: GeneGroupRecord {
                id: id.to_string(),
                label: id.to_string(),
                curation_status: curation_status.to_string(),
                ..GeneGroupRecord::default()
            },
            source_scope: "test".to_string(),
            source_path: "synthetic://gene-groups".to_string(),
        }
    }

    fn member_from_status(
        status: Option<&str>,
    ) -> (
        Option<GeneSetResolvedMember>,
        Vec<String>,
        Vec<GeneSetUnresolvedMember>,
    ) {
        let mut warnings = vec![];
        let mut unresolved = vec![];
        let resolved = GentleEngine::gene_set_member_from_catalog_member(
            &member("TP53", status),
            Some("test_group"),
            Some("synthetic://gene-groups"),
            &[],
            &mut warnings,
            &mut unresolved,
            "catalog_group",
        );
        (resolved, warnings, unresolved)
    }

    fn symbols(rows: &[GeneSetResolvedMember]) -> Vec<String> {
        rows.iter().map(|row| row.symbol.clone()).collect()
    }

    #[test]
    fn gene_set_member_status_gating_includes_warns_or_skips() {
        for status in [None, Some("included")] {
            let (resolved, warnings, unresolved) = member_from_status(status);
            assert!(resolved.is_some());
            assert!(warnings.is_empty());
            assert!(unresolved.is_empty());
        }

        let (resolved, warnings, unresolved) = member_from_status(Some("draft"));
        assert!(resolved.is_some());
        assert!(unresolved.is_empty());
        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("draft") && warning.contains("including"))
        );

        let (resolved, warnings, unresolved) = member_from_status(Some("excluded"));
        assert!(resolved.is_none());
        assert!(unresolved.is_empty());
        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("Skipping excluded gene-set member"))
        );

        let (resolved, warnings, unresolved) = member_from_status(Some("curator_maybe"));
        assert!(resolved.is_some());
        assert!(unresolved.is_empty());
        assert!(warnings.iter().any(|warning| {
            warning.contains("Unrecognized member status 'curator_maybe'")
                && warning.contains("treating as included")
        }));
    }

    #[test]
    fn gene_set_group_curation_gating_rejects_or_warns_by_flag() {
        for status in ["", "reviewed", "included", "curated"] {
            let mut warnings = vec![];
            GentleEngine::gene_set_group_allowed(
                &loaded_group("ok_group", status),
                false,
                false,
                &mut warnings,
            )
            .expect("accepted curation status");
            assert!(warnings.is_empty());
        }

        let mut warnings = vec![];
        let err = GentleEngine::gene_set_group_allowed(
            &loaded_group("draft_group", "draft"),
            false,
            false,
            &mut warnings,
        )
        .expect_err("draft rejected without allow flag");
        assert!(err.message.contains("--allow-draft"));
        assert!(warnings.is_empty());

        GentleEngine::gene_set_group_allowed(
            &loaded_group("draft_group", "draft"),
            true,
            false,
            &mut warnings,
        )
        .expect("draft accepted with allow flag");
        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("Using draft gene group 'draft_group'"))
        );

        warnings.clear();
        let err = GentleEngine::gene_set_group_allowed(
            &loaded_group("old_group", "deprecated"),
            false,
            false,
            &mut warnings,
        )
        .expect_err("deprecated rejected without allow flag");
        assert!(err.message.contains("--allow-deprecated"));
        assert!(warnings.is_empty());

        GentleEngine::gene_set_group_allowed(
            &loaded_group("old_group", "deprecated"),
            false,
            true,
            &mut warnings,
        )
        .expect("deprecated accepted with allow flag");
        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("Using deprecated gene group 'old_group'"))
        );

        warnings.clear();
        GentleEngine::gene_set_group_allowed(
            &loaded_group("unknown_group", "locally_reviewed"),
            false,
            false,
            &mut warnings,
        )
        .expect("unknown curation status is usable");
        assert!(warnings.iter().any(|warning| {
            warning.contains("Unrecognized curation_status 'locally_reviewed'")
                && warning.contains("treating as usable")
        }));
    }

    #[test]
    fn gene_set_random_is_deterministic_source_pinned_and_foreground_filtered() {
        let genes = (0..10)
            .map(|idx| {
                gene(
                    &format!("Gene{idx}"),
                    &format!("G{idx}"),
                    "chr1",
                    100 + idx * 100,
                )
            })
            .collect::<Vec<_>>();

        let mut first_warnings = vec![];
        let (first, first_random) = GentleEngine::gene_set_resolve_random(
            "ToyGenome",
            4,
            17,
            &[],
            &genes,
            "ToyGenome|catalog_a|cache",
            &mut first_warnings,
        );
        let mut second_warnings = vec![];
        let (second, second_random) = GentleEngine::gene_set_resolve_random(
            "ToyGenome",
            4,
            17,
            &[],
            &genes,
            "ToyGenome|catalog_a|cache",
            &mut second_warnings,
        );
        assert_eq!(symbols(&first), symbols(&second));
        assert_eq!(first_random.random_seed, second_random.random_seed);
        assert_eq!(
            first_random.gene_index_source,
            second_random.gene_index_source
        );
        assert!(first_warnings.is_empty());
        assert!(second_warnings.is_empty());

        let mut different_source_warnings = vec![];
        let (different_source, _) = GentleEngine::gene_set_resolve_random(
            "ToyGenome",
            4,
            17,
            &[],
            &genes,
            "ToyGenome|catalog_b|cache",
            &mut different_source_warnings,
        );
        assert_ne!(symbols(&first), symbols(&different_source));

        let mut excluded_warnings = vec![];
        let excluded = vec!["Gene1".to_string(), "G2".to_string()];
        let (filtered, filtered_random) = GentleEngine::gene_set_resolve_random(
            "ToyGenome",
            3,
            23,
            &excluded,
            &genes,
            "ToyGenome|catalog_a|cache",
            &mut excluded_warnings,
        );
        assert_eq!(filtered.len(), 3);
        assert!(filtered.iter().all(
            |row| row.gene_id.as_deref() != Some("G1") && row.gene_id.as_deref() != Some("G2")
        ));
        assert_eq!(filtered_random.universe_size, genes.len());
        assert!(filtered_random.foreground_exclusion_count >= 2);
        assert!(excluded_warnings.is_empty());

        let tiny_genes = genes.iter().take(2).cloned().collect::<Vec<_>>();
        let mut too_large_warnings = vec![];
        let (too_large, _) = GentleEngine::gene_set_resolve_random(
            "ToyGenome",
            5,
            42,
            &[],
            &tiny_genes,
            "ToyGenome|catalog_a|cache",
            &mut too_large_warnings,
        );
        assert_eq!(too_large.len(), tiny_genes.len());
        assert!(
            too_large_warnings
                .iter()
                .any(|warning| warning.contains("fewer than requested 5"))
        );
    }

    #[test]
    fn gene_set_neighbors_are_same_chromosome_and_report_edge_cases() {
        let genes = vec![
            gene("A", "GA", "chr1", 100),
            gene("B", "GB", "chr1", 180),
            gene("C", "GC", "chr1", 250),
            gene("D", "GD", "chr1", 330),
            gene("E", "GE", "chr1", 500),
            gene("Chr2Near", "GCHR2", "chr2", 260),
        ];

        let mut warnings = vec![];
        let mut unresolved = vec![];
        let rows = GentleEngine::gene_set_resolve_neighbors(
            "C",
            Some(2),
            None,
            false,
            &genes,
            &mut warnings,
            &mut unresolved,
        )
        .expect("resolve count neighbors");
        assert_eq!(symbols(&rows), vec!["A", "B", "C", "D", "E"]);
        assert!(!symbols(&rows).contains(&"Chr2Near".to_string()));
        assert!(warnings.is_empty());
        assert!(unresolved.is_empty());

        let rows = GentleEngine::gene_set_resolve_neighbors(
            "C",
            Some(1),
            Some(80),
            true,
            &genes,
            &mut vec![],
            &mut vec![],
        )
        .expect("resolve bp neighbors");
        assert_eq!(symbols(&rows), vec!["B", "D"]);

        let rows = GentleEngine::gene_set_resolve_neighbors(
            "A",
            Some(2),
            None,
            false,
            &genes,
            &mut warnings,
            &mut unresolved,
        )
        .expect("resolve boundary neighbors");
        assert_eq!(symbols(&rows), vec!["A", "B", "C"]);
        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("near the start of chromosome chr1"))
        );

        warnings.clear();
        unresolved.clear();
        let rows = GentleEngine::gene_set_resolve_neighbors(
            "missing",
            Some(1),
            None,
            false,
            &genes,
            &mut warnings,
            &mut unresolved,
        )
        .expect("missing anchor is unresolved, not fatal");
        assert!(rows.is_empty());
        assert_eq!(unresolved.len(), 1);
        assert!(
            unresolved[0]
                .reason
                .contains("neighbor anchor did not resolve")
        );

        let mut duplicated = genes.clone();
        duplicated.push(gene("C", "GC2", "chr2", 700));
        let err = GentleEngine::gene_set_resolve_neighbors(
            "C",
            Some(1),
            None,
            false,
            &duplicated,
            &mut vec![],
            &mut vec![],
        )
        .expect_err("ambiguous anchor is fatal");
        assert!(err.message.contains("resolved to 2 loci"));
        assert!(err.message.contains("chr1:250-299"));
        assert!(err.message.contains("chr2:700-749"));
    }

    #[test]
    fn gene_set_insert_member_collapses_identity_and_preserves_provenance() {
        let mut members = vec![];
        GentleEngine::gene_set_insert_member(
            &mut members,
            GeneSetResolvedMember {
                dedup_key: GentleEngine::gene_set_dedup_key("Tp53", Some("ENSG00000141510")),
                symbol: "Tp53".to_string(),
                gene_id: Some("ENSG00000141510".to_string()),
                aliases: vec!["p53".to_string()],
                chromosome: Some("chr17".to_string()),
                start_1based: Some(7_000),
                contributing_group_ids: vec!["group_a".to_string()],
                provenance: vec![GeneSetProvenanceRow {
                    source_kind: "catalog_group".to_string(),
                    source_id: "group_a".to_string(),
                    ..GeneSetProvenanceRow::default()
                }],
                ..GeneSetResolvedMember::default()
            },
        );
        GentleEngine::gene_set_insert_member(
            &mut members,
            GeneSetResolvedMember {
                dedup_key: GentleEngine::gene_set_dedup_key("TP53", Some("ENSG00000141510")),
                symbol: "TP53".to_string(),
                gene_id: Some("ENSG00000141510".to_string()),
                aliases: vec!["TRP53".to_string(), "p53".to_string()],
                chromosome: Some("chr17".to_string()),
                start_1based: Some(7_000),
                contributing_group_ids: vec!["group_b".to_string()],
                provenance: vec![GeneSetProvenanceRow {
                    source_kind: "external_mapping".to_string(),
                    source_id: "GO:0000123".to_string(),
                    ..GeneSetProvenanceRow::default()
                }],
                ..GeneSetResolvedMember::default()
            },
        );

        assert_eq!(members.len(), 1);
        let merged = &members[0];
        assert_eq!(merged.symbol, "Tp53");
        assert_eq!(merged.contributing_group_ids, vec!["group_a", "group_b"]);
        assert_eq!(merged.provenance.len(), 2);
        assert_eq!(merged.aliases, vec!["p53", "TRP53"]);

        members.push(GeneSetResolvedMember {
            dedup_key: "symbol:no coordinates".to_string(),
            symbol: "NoCoordinates".to_string(),
            ..GeneSetResolvedMember::default()
        });
        members.push(GeneSetResolvedMember {
            dedup_key: "gene_id:gata1".to_string(),
            symbol: "Gata1".to_string(),
            gene_id: Some("GATA1".to_string()),
            chromosome: Some("chrX".to_string()),
            start_1based: Some(1_000),
            ..GeneSetResolvedMember::default()
        });
        GentleEngine::gene_set_sort_members(&mut members);
        assert_eq!(symbols(&members), vec!["Tp53", "Gata1", "NoCoordinates"]);
    }
}
