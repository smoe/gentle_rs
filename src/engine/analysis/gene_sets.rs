//! Gene-set resolution and promoter-cohort construction.
//!
//! This module turns curated gene-group records, explicit user members,
//! local external mappings, prepared-genome neighborhoods, and deterministic
//! random samples into portable analysis operands.

use super::*;
use crate::gene_groups::LoadedGeneGroupRecord;
use gentle_protocol::GeneGroupMember;

#[derive(Clone, Debug)]
struct IndexedGeneRow {
    record: GenomeGeneRecord,
    dedup_key: String,
    symbol: String,
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
            genome_id: resolved_genome_id,
            gene_group_catalog_label,
            genome_catalog_label,
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
            relationship,
            requested_member_count: resolution.resolved_members.len(),
            returned_window_count: windows.len(),
            gene_set_resolution: resolution,
            windows,
            relationship_flags: vec![],
            unresolved_members,
            warnings,
        })
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
            .map(|idx| gene(&format!("Gene{idx}"), &format!("G{idx}"), "chr1", 100 + idx * 100))
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
        assert_eq!(first_random.gene_index_source, second_random.gene_index_source);
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
        assert!(filtered
            .iter()
            .all(|row| row.gene_id.as_deref() != Some("G1") && row.gene_id.as_deref() != Some("G2")));
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
        assert!(too_large_warnings
            .iter()
            .any(|warning| warning.contains("fewer than requested 5")));
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
        assert!(warnings
            .iter()
            .any(|warning| warning.contains("near the start of chromosome chr1")));

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
        assert!(unresolved[0].reason.contains("neighbor anchor did not resolve"));

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
