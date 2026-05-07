//! Catalog-extensible gene-group knowledge layer.
//!
//! Gene groups are deterministic local catalog records that may map to public
//! ontology/resource namespaces such as GO, but are not limited by them.

use gentle_protocol::{
    GENE_GROUP_CATALOG_SCHEMA, GENE_GROUP_DOCTOR_REPORT_SCHEMA, GENE_GROUP_DRAFT_REPORT_SCHEMA,
    GENE_GROUP_LIST_REPORT_SCHEMA, GENE_GROUP_RESOLVE_REPORT_SCHEMA, GENE_GROUP_SHOW_REPORT_SCHEMA,
    GeneGroupCatalog, GeneGroupCatalogSourceReport, GeneGroupDoctorReport, GeneGroupDraftReport,
    GeneGroupExternalMapping, GeneGroupExternalResource, GeneGroupListEntry, GeneGroupListReport,
    GeneGroupMember, GeneGroupRecord, GeneGroupResolveReport, GeneGroupShowReport,
};
use sha1::{Digest, Sha1};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

pub const DEFAULT_GENE_GROUP_CATALOG_PATH: &str = "assets/gene_groups.json";
pub const DEFAULT_GENE_GROUP_DISCOVERY_TOKEN: &str = "gentle://catalog/gene-groups/default";

#[derive(Clone, Debug)]
struct GeneGroupCatalogSourceCandidate {
    scope: &'static str,
    path: PathBuf,
}

#[derive(Clone, Debug)]
pub struct LoadedGeneGroupRecord {
    pub record: GeneGroupRecord,
    pub source_scope: String,
    pub source_path: String,
}

#[derive(Clone, Debug)]
pub struct GeneGroupCatalogIndex {
    catalog_label: String,
    groups: Vec<LoadedGeneGroupRecord>,
    external_resources: Vec<GeneGroupExternalResource>,
    warnings: Vec<String>,
}

#[derive(Clone, Debug, Default)]
pub struct GeneGroupDraftOptions {
    pub description: String,
    pub id: Option<String>,
    pub label: Option<String>,
    pub short_description: Option<String>,
    pub organism: Option<String>,
    pub taxon_id: Option<String>,
    pub symbol_namespace: Option<String>,
    pub aliases: Vec<String>,
    pub tags: Vec<String>,
    pub usages: Vec<String>,
    pub members: Vec<String>,
    pub candidate_members: Vec<String>,
    pub unresolved_candidates: Vec<String>,
    pub go_mappings: Vec<String>,
    pub provenance: Option<String>,
    pub agent_provider: Option<String>,
    pub agent_model: Option<String>,
    pub agent_generated_at_utc: Option<String>,
    pub output_path: Option<String>,
}

/// Stable human-facing label for default gene-group discovery.
pub fn default_gene_group_catalog_discovery_label() -> &'static str {
    "default gene-group catalog discovery"
}

fn normalize_lookup(raw: &str) -> String {
    raw.chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() {
                ch.to_ascii_lowercase()
            } else {
                ' '
            }
        })
        .collect::<String>()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

fn slugify_identifier(raw: &str) -> String {
    let mut out = String::new();
    let mut previous_underscore = false;
    for ch in raw.chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch.to_ascii_lowercase());
            previous_underscore = false;
        } else if !previous_underscore && !out.is_empty() {
            out.push('_');
            previous_underscore = true;
        }
    }
    while out.ends_with('_') {
        out.pop();
    }
    if out.is_empty() {
        "draft_gene_group".to_string()
    } else {
        out
    }
}

fn first_sentence(raw: &str) -> String {
    raw.split(['.', '\n'])
        .map(str::trim)
        .find(|value| !value.is_empty())
        .unwrap_or(raw.trim())
        .to_string()
}

fn derive_label_from_description(description: &str) -> String {
    let sentence = first_sentence(description);
    let words = sentence
        .split_whitespace()
        .take(10)
        .collect::<Vec<_>>()
        .join(" ");
    let label = words
        .trim_matches(|ch: char| !ch.is_ascii_alphanumeric())
        .trim();
    if label.is_empty() {
        "Draft gene group".to_string()
    } else {
        let mut chars = label.chars();
        match chars.next() {
            Some(first) => format!("{}{}", first.to_ascii_uppercase(), chars.as_str()),
            None => "Draft gene group".to_string(),
        }
    }
}

fn derive_short_description(description: &str) -> String {
    let sentence = first_sentence(description);
    if sentence.chars().count() <= 180 {
        sentence
    } else {
        let mut out = sentence.chars().take(177).collect::<String>();
        out.push_str("...");
        out
    }
}

fn normalized_unique(values: Vec<String>) -> Vec<String> {
    let mut seen = BTreeSet::<String>::new();
    let mut out = vec![];
    for value in values {
        let trimmed = value.trim();
        if trimmed.is_empty() {
            continue;
        }
        let key = normalize_lookup(trimmed);
        if seen.insert(key) {
            out.push(trimmed.to_string());
        }
    }
    out
}

fn parse_draft_candidate(raw: &str) -> Option<(String, Option<String>)> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let (symbol, evidence_note) = trimmed
        .split_once('=')
        .map(|(symbol, note)| (symbol.trim(), Some(note.trim())))
        .unwrap_or((trimmed, None));
    if symbol.is_empty() {
        return None;
    }
    let note = evidence_note
        .filter(|value| !value.is_empty())
        .map(str::to_string);
    Some((symbol.to_string(), note))
}

fn gene_group_catalog_discovery_candidates() -> Vec<GeneGroupCatalogSourceCandidate> {
    let mut candidates = vec![];
    let built_in_root = crate::genomes::configured_builtin_asset_root();
    candidates.push(GeneGroupCatalogSourceCandidate {
        scope: "built-in",
        path: built_in_root.join(DEFAULT_GENE_GROUP_CATALOG_PATH),
    });
    candidates.push(GeneGroupCatalogSourceCandidate {
        scope: "built-in",
        path: built_in_root.join("assets").join("gene_groups.d"),
    });

    let system_root = crate::genomes::configured_system_config_root();
    candidates.push(GeneGroupCatalogSourceCandidate {
        scope: "system",
        path: system_root.join("catalogs").join("gene_groups.json"),
    });
    candidates.push(GeneGroupCatalogSourceCandidate {
        scope: "system",
        path: system_root.join("catalogs").join("gene_groups.d"),
    });

    if let Some(user_root) = crate::genomes::configured_user_config_root() {
        candidates.push(GeneGroupCatalogSourceCandidate {
            scope: "user",
            path: user_root.join("catalogs").join("gene_groups.json"),
        });
        candidates.push(GeneGroupCatalogSourceCandidate {
            scope: "user",
            path: user_root.join("catalogs").join("gene_groups.d"),
        });
    }

    if let Some(project_root) = crate::genomes::configured_project_root() {
        candidates.push(GeneGroupCatalogSourceCandidate {
            scope: "project",
            path: project_root
                .join(".gentle")
                .join("catalogs")
                .join("gene_groups.json"),
        });
        candidates.push(GeneGroupCatalogSourceCandidate {
            scope: "project",
            path: project_root
                .join(".gentle")
                .join("catalogs")
                .join("gene_groups.d"),
        });
    }

    dedup_gene_group_catalog_sources(candidates)
}

fn dedup_gene_group_catalog_sources(
    candidates: Vec<GeneGroupCatalogSourceCandidate>,
) -> Vec<GeneGroupCatalogSourceCandidate> {
    let mut seen = BTreeSet::<String>::new();
    let mut out = vec![];
    for candidate in candidates {
        let key = candidate.path.to_string_lossy().to_string();
        if seen.insert(key) {
            out.push(candidate);
        }
    }
    out
}

fn discovered_gene_group_catalog_sources() -> Vec<GeneGroupCatalogSourceCandidate> {
    gene_group_catalog_discovery_candidates()
        .into_iter()
        .filter(|candidate| fs::metadata(&candidate.path).is_ok())
        .collect()
}

impl GeneGroupCatalogIndex {
    pub fn from_explicit_path(path: &str) -> Result<Self, String> {
        Self::from_sources(
            &[GeneGroupCatalogSourceCandidate {
                scope: "explicit",
                path: PathBuf::from(path),
            }],
            path.to_string(),
        )
    }

    pub fn from_default_discovery() -> Result<Self, String> {
        let sources = discovered_gene_group_catalog_sources();
        if sources.is_empty() {
            return Err(format!(
                "No gene-group catalog sources were found for {}",
                default_gene_group_catalog_discovery_label()
            ));
        }
        Self::from_sources(
            &sources,
            default_gene_group_catalog_discovery_label().to_string(),
        )
    }

    fn from_sources(
        sources: &[GeneGroupCatalogSourceCandidate],
        catalog_label: String,
    ) -> Result<Self, String> {
        let mut groups = vec![];
        let mut seen_group_ids = BTreeMap::<String, String>::new();
        let mut external_resources = BTreeMap::<String, GeneGroupExternalResource>::new();
        let mut warnings = vec![];
        for source in sources {
            Self::merge_catalog_source(
                source,
                &mut groups,
                &mut seen_group_ids,
                &mut external_resources,
                &mut warnings,
            )?;
        }
        groups.sort_by(|a, b| a.record.id.cmp(&b.record.id));
        Ok(Self {
            catalog_label,
            groups,
            external_resources: external_resources.into_values().collect(),
            warnings,
        })
    }

    fn merge_catalog_source(
        source: &GeneGroupCatalogSourceCandidate,
        groups: &mut Vec<LoadedGeneGroupRecord>,
        seen_group_ids: &mut BTreeMap<String, String>,
        external_resources: &mut BTreeMap<String, GeneGroupExternalResource>,
        warnings: &mut Vec<String>,
    ) -> Result<(), String> {
        let metadata = fs::metadata(&source.path).map_err(|e| {
            format!(
                "Could not read gene-group catalog '{}': {e}",
                source.path.display()
            )
        })?;
        if metadata.is_dir() {
            let mut json_files = json_files_in_dir(&source.path)?;
            json_files.sort();
            if json_files.is_empty() {
                return Err(format!(
                    "Gene-group catalog directory '{}' does not contain any .json files",
                    source.path.display()
                ));
            }
            for file_path in json_files {
                Self::merge_catalog_file(
                    source.scope,
                    &file_path,
                    groups,
                    seen_group_ids,
                    external_resources,
                    warnings,
                )?;
            }
            return Ok(());
        }
        Self::merge_catalog_file(
            source.scope,
            &source.path,
            groups,
            seen_group_ids,
            external_resources,
            warnings,
        )
    }

    fn merge_catalog_file(
        scope: &'static str,
        path: &Path,
        groups: &mut Vec<LoadedGeneGroupRecord>,
        seen_group_ids: &mut BTreeMap<String, String>,
        external_resources: &mut BTreeMap<String, GeneGroupExternalResource>,
        warnings: &mut Vec<String>,
    ) -> Result<(), String> {
        let text = fs::read_to_string(path).map_err(|e| {
            format!(
                "Could not read gene-group catalog '{}': {e}",
                path.display()
            )
        })?;
        let parsed = parse_gene_group_catalog(&text, path)?;
        for resource in parsed.external_resources {
            let resource_id = resource.id.trim().to_string();
            if resource_id.is_empty() {
                warnings.push(format!(
                    "Ignoring external resource without id in '{}'",
                    path.display()
                ));
                continue;
            }
            if let Some(previous) = external_resources.get(&resource_id) {
                if previous.namespace != resource.namespace {
                    return Err(format!(
                        "Duplicate external resource '{}' has namespaces '{}' and '{}' while loading '{}'",
                        resource_id,
                        previous.namespace,
                        resource.namespace,
                        path.display()
                    ));
                }
                continue;
            }
            external_resources.insert(resource_id, resource);
        }
        for group in parsed.groups {
            validate_group_minimum(&group, path)?;
            if let Some(previous) = seen_group_ids.get(&group.id) {
                return Err(format!(
                    "Duplicate gene-group id '{}' found in '{}' and '{}'",
                    group.id,
                    previous,
                    path.display()
                ));
            }
            seen_group_ids.insert(group.id.clone(), path.display().to_string());
            groups.push(LoadedGeneGroupRecord {
                record: group,
                source_scope: scope.to_string(),
                source_path: path.display().to_string(),
            });
        }
        Ok(())
    }

    pub fn list(&self, filter: Option<&str>) -> Vec<GeneGroupListEntry> {
        let mut rows = self
            .groups
            .iter()
            .filter(|row| group_matches_filter(row, filter))
            .map(loaded_group_to_list_entry)
            .collect::<Vec<_>>();
        rows.sort_by(|a, b| a.id.cmp(&b.id));
        rows
    }

    pub fn resolve_exact(&self, query: &str) -> Vec<LoadedGeneGroupRecord> {
        let normalized = normalize_lookup(query);
        if normalized.is_empty() {
            return vec![];
        }
        let mut rows = self
            .groups
            .iter()
            .filter(|row| group_lookup_keys(&row.record).contains(&normalized))
            .cloned()
            .collect::<Vec<_>>();
        rows.sort_by(|a, b| a.record.id.cmp(&b.record.id));
        rows
    }

    pub fn tf_query_groups(&self) -> Vec<LoadedGeneGroupRecord> {
        self.groups
            .iter()
            .filter(|row| {
                row.record
                    .usages
                    .iter()
                    .any(|usage| normalize_lookup(usage) == "tf query")
            })
            .cloned()
            .collect()
    }
}

fn parse_gene_group_catalog(text: &str, path: &Path) -> Result<GeneGroupCatalog, String> {
    let parsed: GeneGroupCatalog = serde_json::from_str(text).map_err(|e| {
        format!(
            "Could not parse gene-group catalog '{}': {e}",
            path.display()
        )
    })?;
    if !parsed.schema.is_empty() && parsed.schema != GENE_GROUP_CATALOG_SCHEMA {
        return Err(format!(
            "Gene-group catalog '{}' has unsupported schema '{}' (expected '{}')",
            path.display(),
            parsed.schema,
            GENE_GROUP_CATALOG_SCHEMA
        ));
    }
    Ok(parsed)
}

fn validate_group_minimum(group: &GeneGroupRecord, path: &Path) -> Result<(), String> {
    if group.id.trim().is_empty() {
        return Err(format!(
            "Gene-group catalog '{}' contains a group without id",
            path.display()
        ));
    }
    if group.label.trim().is_empty() {
        return Err(format!(
            "Gene-group '{}' in '{}' has an empty label",
            group.id,
            path.display()
        ));
    }
    if group
        .members
        .iter()
        .any(|member| member.symbol.trim().is_empty())
    {
        return Err(format!(
            "Gene-group '{}' in '{}' contains a member without symbol",
            group.id,
            path.display()
        ));
    }
    Ok(())
}

fn json_files_in_dir(path: &Path) -> Result<Vec<PathBuf>, String> {
    let mut files = vec![];
    for row in fs::read_dir(path).map_err(|e| {
        format!(
            "Could not read gene-group catalog directory '{}': {e}",
            path.display()
        )
    })? {
        let row = row.map_err(|e| {
            format!(
                "Could not inspect gene-group catalog directory '{}': {e}",
                path.display()
            )
        })?;
        if !row
            .file_type()
            .map_err(|e| format!("Could not inspect '{}': {e}", row.path().display()))?
            .is_file()
        {
            continue;
        }
        let file_path = row.path();
        if file_path
            .extension()
            .and_then(|value| value.to_str())
            .map(|value| value.eq_ignore_ascii_case("json"))
            .unwrap_or(false)
        {
            files.push(file_path);
        }
    }
    Ok(files)
}

fn group_lookup_keys(group: &GeneGroupRecord) -> BTreeSet<String> {
    let mut keys = BTreeSet::new();
    keys.insert(normalize_lookup(&group.id));
    keys.insert(normalize_lookup(&group.label));
    for alias in &group.aliases {
        keys.insert(normalize_lookup(alias));
    }
    keys.into_iter().filter(|key| !key.is_empty()).collect()
}

fn group_search_terms(group: &GeneGroupRecord) -> Vec<String> {
    let mut terms = vec![
        group.id.clone(),
        group.label.clone(),
        group.short_description.clone(),
        group.long_definition.clone(),
        group.organism.clone().unwrap_or_default(),
        group.taxon_id.clone().unwrap_or_default(),
        group.symbol_namespace.clone().unwrap_or_default(),
        group.curation_status.clone(),
    ];
    terms.extend(group.aliases.clone());
    terms.extend(group.usages.clone());
    terms.extend(group.tags.clone());
    terms.extend(group.members.iter().flat_map(|member| {
        let mut values = vec![member.symbol.clone()];
        values.extend(member.aliases.clone());
        values
    }));
    terms.extend(group.external_mappings.iter().flat_map(|mapping| {
        vec![
            mapping.namespace.clone(),
            mapping.id.clone(),
            mapping.label.clone().unwrap_or_default(),
            mapping.note.clone().unwrap_or_default(),
        ]
    }));
    terms
}

fn group_matches_filter(row: &LoadedGeneGroupRecord, filter: Option<&str>) -> bool {
    let Some(filter) = filter.map(str::trim).filter(|value| !value.is_empty()) else {
        return true;
    };
    let needle = filter.to_ascii_lowercase();
    group_search_terms(&row.record)
        .iter()
        .any(|term| term.to_ascii_lowercase().contains(&needle))
}

fn loaded_group_to_list_entry(row: &LoadedGeneGroupRecord) -> GeneGroupListEntry {
    GeneGroupListEntry {
        id: row.record.id.clone(),
        label: row.record.label.clone(),
        aliases: row.record.aliases.clone(),
        short_description: row.record.short_description.clone(),
        organism: row.record.organism.clone(),
        taxon_id: row.record.taxon_id.clone(),
        symbol_namespace: row.record.symbol_namespace.clone(),
        usages: row.record.usages.clone(),
        tags: row.record.tags.clone(),
        member_count: row.record.members.len(),
        curation_status: row.record.curation_status.clone(),
        external_mappings: row.record.external_mappings.clone(),
        source_scope: row.source_scope.clone(),
        source_path: row.source_path.clone(),
    }
}

pub fn load_gene_group_catalog(
    catalog_path: Option<&str>,
) -> Result<GeneGroupCatalogIndex, String> {
    if let Some(path) = catalog_path
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        if path == DEFAULT_GENE_GROUP_DISCOVERY_TOKEN {
            return GeneGroupCatalogIndex::from_default_discovery();
        }
        GeneGroupCatalogIndex::from_explicit_path(path)
    } else {
        GeneGroupCatalogIndex::from_default_discovery()
    }
}

pub fn list_gene_groups(
    catalog_path: Option<&str>,
    filter: Option<&str>,
) -> Result<GeneGroupListReport, String> {
    let catalog = load_gene_group_catalog(catalog_path)?;
    let groups = catalog.list(filter);
    Ok(GeneGroupListReport {
        schema: GENE_GROUP_LIST_REPORT_SCHEMA.to_string(),
        catalog_label: catalog.catalog_label,
        filter: filter.map(str::to_string),
        returned_group_count: groups.len(),
        groups,
        external_resources: catalog.external_resources,
        warnings: catalog.warnings,
    })
}

pub fn show_gene_group(
    group_id: &str,
    catalog_path: Option<&str>,
) -> Result<GeneGroupShowReport, String> {
    let catalog = load_gene_group_catalog(catalog_path)?;
    let rows = catalog.resolve_exact(group_id);
    let row = rows.first().ok_or_else(|| {
        format!(
            "Gene group '{}' was not found in {}",
            group_id, catalog.catalog_label
        )
    })?;
    Ok(GeneGroupShowReport {
        schema: GENE_GROUP_SHOW_REPORT_SCHEMA.to_string(),
        catalog_label: catalog.catalog_label,
        group: row.record.clone(),
        source_scope: row.source_scope.clone(),
        source_path: row.source_path.clone(),
        external_resources: catalog.external_resources,
        warnings: catalog.warnings,
    })
}

pub fn resolve_gene_group(
    query: &str,
    catalog_path: Option<&str>,
) -> Result<GeneGroupResolveReport, String> {
    let catalog = load_gene_group_catalog(catalog_path)?;
    let rows = catalog
        .resolve_exact(query)
        .iter()
        .map(loaded_group_to_list_entry)
        .collect::<Vec<_>>();
    Ok(GeneGroupResolveReport {
        schema: GENE_GROUP_RESOLVE_REPORT_SCHEMA.to_string(),
        catalog_label: catalog.catalog_label,
        query: query.to_string(),
        normalized_query: normalize_lookup(query),
        matched_group_count: rows.len(),
        groups: rows,
        warnings: catalog.warnings,
    })
}

pub fn draft_gene_group(options: GeneGroupDraftOptions) -> Result<GeneGroupDraftReport, String> {
    let description = options.description.trim();
    if description.is_empty() {
        return Err("gene-groups draft requires a non-empty --description".to_string());
    }
    let label = options
        .label
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .unwrap_or_else(|| derive_label_from_description(description));
    let id = options
        .id
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(slugify_identifier)
        .unwrap_or_else(|| slugify_identifier(&label));
    let short_description = options
        .short_description
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .unwrap_or_else(|| derive_short_description(description));
    let aliases = normalized_unique(options.aliases);
    let mut tags = normalized_unique(options.tags);
    if tags.is_empty() {
        tags.push("draft".to_string());
    }
    let mut usages = normalized_unique(options.usages);
    if usages.is_empty() {
        usages.push("gene_group".to_string());
    }
    let provenance = options.provenance.clone().unwrap_or_else(|| {
        "AI/user-assisted draft from gene-groups draft; requires review before promotion."
            .to_string()
    });
    let user_members = normalized_unique(options.members)
        .into_iter()
        .map(|symbol| GeneGroupMember {
            symbol,
            evidence_note: Some(
                "Draft candidate membership; review evidence before promoting this group."
                    .to_string(),
            ),
            confidence: Some("draft_candidate".to_string()),
            status: Some("draft".to_string()),
            provenance: Some(provenance.clone()),
            ..GeneGroupMember::default()
        })
        .collect::<Vec<_>>();
    let agent_provenance = {
        let mut parts = vec!["agent-assisted candidate draft".to_string()];
        if let Some(provider) = options
            .agent_provider
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            parts.push(format!("provider={provider}"));
        }
        if let Some(model) = options
            .agent_model
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            parts.push(format!("model={model}"));
        }
        if let Some(generated_at) = options
            .agent_generated_at_utc
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            parts.push(format!("generated_at_utc={generated_at}"));
        }
        parts.push("requires review before promotion".to_string());
        parts.join("; ")
    };
    let mut candidate_seen = BTreeSet::<String>::new();
    let candidate_members = options
        .candidate_members
        .iter()
        .filter_map(|raw| parse_draft_candidate(raw))
        .filter_map(|(symbol, evidence_note)| {
            let key = normalize_lookup(&symbol);
            if !candidate_seen.insert(key) {
                return None;
            }
            Some(GeneGroupMember {
                symbol,
                evidence_note: Some(evidence_note.unwrap_or_else(|| {
                    "Agent-suggested candidate membership; review evidence before promoting this group."
                        .to_string()
                })),
                confidence: Some("agent_suggested".to_string()),
                status: Some("draft".to_string()),
                provenance: Some(agent_provenance.clone()),
                ..GeneGroupMember::default()
            })
        })
        .collect::<Vec<_>>();
    let mut member_seen = BTreeSet::<String>::new();
    let mut members = vec![];
    for member in user_members.iter().chain(candidate_members.iter()) {
        if member_seen.insert(normalize_lookup(&member.symbol)) {
            members.push(member.clone());
        }
    }
    let unresolved_candidates = normalized_unique(options.unresolved_candidates);
    let external_mappings = normalized_unique(options.go_mappings)
        .into_iter()
        .map(|go_id| {
            if !is_valid_go_id(&go_id) {
                return Err(format!(
                    "Invalid --go mapping '{go_id}' for gene-groups draft (expected GO:0000000)"
                ));
            }
            Ok(GeneGroupExternalMapping {
                namespace: "GO".to_string(),
                id: go_id,
                relationship: Some("draft_external_anchor".to_string()),
                note: Some(
                    "Draft GO mapping supplied during gene-group drafting; review before promotion."
                        .to_string(),
                ),
                ..GeneGroupExternalMapping::default()
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let mut external_resources = vec![];
    if !external_mappings.is_empty() {
        external_resources.push(GeneGroupExternalResource {
            id: "gene_ontology".to_string(),
            label: "Gene Ontology".to_string(),
            namespace: "GO".to_string(),
            role: "external_ontology_mapping_namespace".to_string(),
            homepage_url: Some("https://geneontology.org/".to_string()),
            notes: vec![
                "Draft fragment declares GO as an external mapping namespace only.".to_string(),
            ],
        });
    }
    let group = GeneGroupRecord {
        id,
        label,
        aliases,
        short_description,
        long_definition: description.to_string(),
        organism: options.organism.filter(|value| !value.trim().is_empty()),
        taxon_id: options.taxon_id.filter(|value| !value.trim().is_empty()),
        symbol_namespace: options
            .symbol_namespace
            .filter(|value| !value.trim().is_empty()),
        usages,
        tags,
        members,
        external_mappings,
        curation_status: "draft".to_string(),
        source_kind: Some("ai_assisted_draft".to_string()),
        provenance: Some(provenance),
        notes: vec![
            "Generated as a review-gated draft; do not treat as curated until promoted by a user or project policy.".to_string(),
        ],
        ..GeneGroupRecord::default()
    };
    let catalog_fragment = GeneGroupCatalog {
        schema: GENE_GROUP_CATALOG_SCHEMA.to_string(),
        external_resources,
        groups: vec![group.clone()],
    };
    let input_description_sha1 = format!("{:x}", Sha1::digest(description.as_bytes()));
    let mut report = GeneGroupDraftReport {
        schema: GENE_GROUP_DRAFT_REPORT_SCHEMA.to_string(),
        generation_method: "deterministic_review_gated_draft".to_string(),
        review_required: true,
        input_description_sha1,
        agent_provider: options
            .agent_provider
            .filter(|value| !value.trim().is_empty()),
        agent_model: options.agent_model.filter(|value| !value.trim().is_empty()),
        agent_generated_at_utc: options
            .agent_generated_at_utc
            .filter(|value| !value.trim().is_empty()),
        user_member_count: user_members.len(),
        candidate_member_count: candidate_members.len(),
        candidate_members,
        unresolved_candidates,
        output_path: options.output_path.clone(),
        group,
        catalog_fragment,
        warnings: vec![],
    };
    if report.catalog_fragment.groups[0].members.is_empty() {
        report.warnings.push(
            "Draft has no candidate members yet; add members before using it for analysis."
                .to_string(),
        );
    }
    if !report.unresolved_candidates.is_empty() {
        report.warnings.push(format!(
            "Draft records {} unresolved candidate(s); review/normalize before promotion.",
            report.unresolved_candidates.len()
        ));
    }
    if let Some(path) = options
        .output_path
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        let parent = Path::new(path)
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        fs::create_dir_all(parent).map_err(|e| {
            format!("Could not create gene-group draft output directory for '{path}': {e}")
        })?;
        let mut text = serde_json::to_string_pretty(&report.catalog_fragment)
            .map_err(|e| format!("Could not serialize gene-group draft fragment: {e}"))?;
        text.push('\n');
        fs::write(path, text)
            .map_err(|e| format!("Could not write gene-group draft fragment '{path}': {e}"))?;
    }
    Ok(report)
}

pub fn resolve_tf_query_gene_group(query: &str) -> Option<LoadedGeneGroupRecord> {
    let catalog = load_gene_group_catalog(None).ok()?;
    let normalized = normalize_lookup(query);
    catalog
        .tf_query_groups()
        .into_iter()
        .find(|row| group_lookup_keys(&row.record).contains(&normalized))
}

pub fn list_tf_query_gene_groups() -> Vec<LoadedGeneGroupRecord> {
    load_gene_group_catalog(None)
        .map(|catalog| catalog.tf_query_groups())
        .unwrap_or_default()
}

pub fn doctor_gene_group_catalog(catalog_path: Option<&str>) -> GeneGroupDoctorReport {
    let explicit = catalog_path
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .filter(|value| *value != DEFAULT_GENE_GROUP_DISCOVERY_TOKEN);
    let catalog_label = explicit
        .map(str::to_string)
        .unwrap_or_else(|| default_gene_group_catalog_discovery_label().to_string());
    let candidates = if let Some(path) = explicit {
        vec![GeneGroupCatalogSourceCandidate {
            scope: "explicit",
            path: PathBuf::from(path),
        }]
    } else {
        gene_group_catalog_discovery_candidates()
    };
    let mut report = GeneGroupDoctorReport {
        schema: GENE_GROUP_DOCTOR_REPORT_SCHEMA.to_string(),
        catalog_label,
        source_count: candidates.len(),
        ..GeneGroupDoctorReport::default()
    };
    let mut seen_group_ids = BTreeMap::<String, String>::new();
    let mut seen_aliases = BTreeMap::<String, String>::new();
    let mut resource_namespaces = BTreeSet::<String>::new();
    let mut seen_resource_ids = BTreeMap::<String, String>::new();

    for candidate in candidates {
        inspect_doctor_source(
            &candidate,
            explicit.is_some(),
            &mut report,
            &mut seen_group_ids,
            &mut seen_aliases,
            &mut resource_namespaces,
            &mut seen_resource_ids,
        );
    }

    for source in &report.sources {
        report.warning_count += source.warnings.len();
        report.error_count += source.errors.len();
    }
    report.warning_count += report.warnings.len();
    report.error_count += report.errors.len();
    report.parsed_source_count = report
        .sources
        .iter()
        .filter(|source| source.status == "ok")
        .count();
    report
}

fn inspect_doctor_source(
    candidate: &GeneGroupCatalogSourceCandidate,
    explicit: bool,
    report: &mut GeneGroupDoctorReport,
    seen_group_ids: &mut BTreeMap<String, String>,
    seen_aliases: &mut BTreeMap<String, String>,
    resource_namespaces: &mut BTreeSet<String>,
    seen_resource_ids: &mut BTreeMap<String, String>,
) {
    let metadata = match fs::metadata(&candidate.path) {
        Ok(metadata) => metadata,
        Err(e) => {
            let mut source = GeneGroupCatalogSourceReport {
                scope: candidate.scope.to_string(),
                path: candidate.path.display().to_string(),
                status: "missing".to_string(),
                ..GeneGroupCatalogSourceReport::default()
            };
            if explicit {
                source.errors.push(format!(
                    "Explicit gene-group catalog '{}' is missing or unreadable: {e}",
                    candidate.path.display()
                ));
            }
            report.sources.push(source);
            return;
        }
    };
    if metadata.is_dir() {
        match json_files_in_dir(&candidate.path) {
            Ok(files) if files.is_empty() => {
                report.sources.push(GeneGroupCatalogSourceReport {
                    scope: candidate.scope.to_string(),
                    path: candidate.path.display().to_string(),
                    status: "empty_directory".to_string(),
                    warnings: vec!["Directory contains no .json catalog fragments".to_string()],
                    ..GeneGroupCatalogSourceReport::default()
                });
            }
            Ok(mut files) => {
                files.sort();
                for file_path in files {
                    inspect_doctor_file(
                        candidate.scope,
                        &file_path,
                        report,
                        seen_group_ids,
                        seen_aliases,
                        resource_namespaces,
                        seen_resource_ids,
                    );
                }
            }
            Err(e) => {
                report.sources.push(GeneGroupCatalogSourceReport {
                    scope: candidate.scope.to_string(),
                    path: candidate.path.display().to_string(),
                    status: "error".to_string(),
                    errors: vec![e],
                    ..GeneGroupCatalogSourceReport::default()
                });
            }
        }
        return;
    }
    inspect_doctor_file(
        candidate.scope,
        &candidate.path,
        report,
        seen_group_ids,
        seen_aliases,
        resource_namespaces,
        seen_resource_ids,
    );
}

fn inspect_doctor_file(
    scope: &'static str,
    path: &Path,
    report: &mut GeneGroupDoctorReport,
    seen_group_ids: &mut BTreeMap<String, String>,
    seen_aliases: &mut BTreeMap<String, String>,
    resource_namespaces: &mut BTreeSet<String>,
    seen_resource_ids: &mut BTreeMap<String, String>,
) {
    let mut source = GeneGroupCatalogSourceReport {
        scope: scope.to_string(),
        path: path.display().to_string(),
        status: "ok".to_string(),
        ..GeneGroupCatalogSourceReport::default()
    };
    let text = match fs::read_to_string(path) {
        Ok(text) => text,
        Err(e) => {
            source.status = "error".to_string();
            source
                .errors
                .push(format!("Could not read gene-group catalog: {e}"));
            report.sources.push(source);
            return;
        }
    };
    source.sha1 = Some(format!("{:x}", Sha1::digest(text.as_bytes())));
    let parsed = match parse_gene_group_catalog(&text, path) {
        Ok(parsed) => parsed,
        Err(e) => {
            source.status = "error".to_string();
            source.errors.push(e);
            report.sources.push(source);
            return;
        }
    };
    source.group_count = parsed.groups.len();
    source.external_resource_count = parsed.external_resources.len();
    report.group_count += parsed.groups.len();
    report.external_resource_count += parsed.external_resources.len();

    for resource in parsed.external_resources {
        if resource.id.trim().is_empty() {
            source
                .errors
                .push("External resource entry has empty id".to_string());
            continue;
        }
        if resource.namespace.trim().is_empty() {
            source.errors.push(format!(
                "External resource '{}' has empty namespace",
                resource.id
            ));
        }
        if let Some(previous) = seen_resource_ids.insert(resource.id.clone(), source.path.clone()) {
            source.errors.push(format!(
                "External resource '{}' is duplicated in '{}' and '{}'",
                resource.id, previous, source.path
            ));
        }
        resource_namespaces.insert(resource.namespace.clone());
    }

    let source_path = source.path.clone();
    for group in parsed.groups {
        validate_doctor_group(
            &group,
            &source_path,
            &mut source,
            seen_group_ids,
            seen_aliases,
            resource_namespaces,
        );
    }
    if !source.errors.is_empty() {
        source.status = "error".to_string();
    } else if !source.warnings.is_empty() {
        source.status = "warning".to_string();
    }
    report.sources.push(source);
}

fn validate_doctor_group(
    group: &GeneGroupRecord,
    source_path: &str,
    source: &mut GeneGroupCatalogSourceReport,
    seen_group_ids: &mut BTreeMap<String, String>,
    seen_aliases: &mut BTreeMap<String, String>,
    resource_namespaces: &BTreeSet<String>,
) {
    if group.id.trim().is_empty() {
        source.errors.push("Group has empty id".to_string());
        return;
    }
    if let Some(previous) = seen_group_ids.insert(group.id.clone(), source_path.to_string()) {
        source.errors.push(format!(
            "Group id '{}' is duplicated in '{}' and '{}'",
            group.id, previous, source_path
        ));
    }
    if group.label.trim().is_empty() {
        source
            .errors
            .push(format!("Group '{}' has empty label", group.id));
    }
    if group.members.is_empty() {
        source
            .warnings
            .push(format!("Group '{}' has no members", group.id));
    }
    let mut member_symbols = BTreeSet::<String>::new();
    for member in &group.members {
        let symbol = member.symbol.trim();
        if symbol.is_empty() {
            source.errors.push(format!(
                "Group '{}' has a member with empty symbol",
                group.id
            ));
        } else if !member_symbols.insert(normalize_lookup(symbol)) {
            source.warnings.push(format!(
                "Group '{}' repeats member symbol '{}'",
                group.id, symbol
            ));
        }
    }
    for key in group_lookup_keys(group) {
        if let Some(previous) = seen_aliases.insert(key.clone(), group.id.clone()) {
            if previous != group.id {
                source.errors.push(format!(
                    "Alias/key '{}' collides between groups '{}' and '{}'",
                    key, previous, group.id
                ));
            }
        }
    }
    if group.curation_status.trim().is_empty() {
        source
            .warnings
            .push(format!("Group '{}' has empty curation_status", group.id));
    }
    for mapping in &group.external_mappings {
        if mapping.namespace.trim().is_empty() || mapping.id.trim().is_empty() {
            source.errors.push(format!(
                "Group '{}' has an external mapping with empty namespace or id",
                group.id
            ));
        }
        if mapping.namespace == "GO" && !is_valid_go_id(&mapping.id) {
            source.errors.push(format!(
                "Group '{}' has malformed GO mapping id '{}'",
                group.id, mapping.id
            ));
        }
        if !mapping.namespace.trim().is_empty()
            && !resource_namespaces.contains(mapping.namespace.as_str())
        {
            source.warnings.push(format!(
                "Group '{}' maps to namespace '{}' but no external_resource declared that namespace in earlier/effective sources",
                group.id, mapping.namespace
            ));
        }
    }
}

fn is_valid_go_id(raw: &str) -> bool {
    let Some(rest) = raw.strip_prefix("GO:") else {
        return false;
    };
    rest.len() == 7 && rest.chars().all(|ch| ch.is_ascii_digit())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn builtin_gene_groups_resolve_yamanaka_and_report_go_resource() {
        let report = resolve_gene_group("Yamanaka factors", None).expect("resolve yamanaka");
        assert_eq!(report.schema, GENE_GROUP_RESOLVE_REPORT_SCHEMA);
        assert_eq!(report.matched_group_count, 1);
        assert_eq!(report.groups[0].id, "yamanaka_factors");

        let splicing = resolve_gene_group("regulation of alternative splicing", None)
            .expect("resolve splicing");
        assert_eq!(splicing.matched_group_count, 1);
        assert_eq!(splicing.groups[0].id, "regulation_of_alternative_splicing");
        assert!(
            splicing.groups[0]
                .external_mappings
                .iter()
                .any(|mapping| mapping.namespace == "GO" && mapping.id == "GO:0000381")
        );

        let list = list_gene_groups(None, Some("ontology")).expect("list ontology");
        assert!(
            list.external_resources
                .iter()
                .any(|resource| resource.namespace == "GO")
        );
    }

    #[test]
    fn draft_gene_group_writes_review_gated_fragment() {
        let dir = tempdir().expect("tempdir");
        let output = dir.path().join("draft.json");
        let report = draft_gene_group(GeneGroupDraftOptions {
            description: "Genes regulating alternative splice-site selection in a project-specific pancreatic cancer context.".to_string(),
            members: vec!["RBFOX2".to_string(), "PTBP1".to_string()],
            candidate_members: vec![
                "QKI=agent-suggested STAR-family splicing regulator".to_string()
            ],
            unresolved_candidates: vec!["ambiguous PTB family".to_string()],
            go_mappings: vec!["GO:0000381".to_string()],
            agent_provider: Some("Codex".to_string()),
            agent_model: Some("test-model".to_string()),
            output_path: Some(output.to_string_lossy().to_string()),
            ..GeneGroupDraftOptions::default()
        })
        .expect("draft gene group");
        assert_eq!(report.schema, GENE_GROUP_DRAFT_REPORT_SCHEMA);
        assert!(report.review_required);
        assert_eq!(report.user_member_count, 2);
        assert_eq!(report.candidate_member_count, 1);
        assert_eq!(report.candidate_members[0].symbol, "QKI");
        assert_eq!(report.unresolved_candidates, vec!["ambiguous PTB family"]);
        assert_eq!(report.catalog_fragment.groups.len(), 1);
        assert_eq!(report.catalog_fragment.groups[0].curation_status, "draft");
        assert_eq!(report.catalog_fragment.groups[0].members.len(), 3);
        let written = fs::read_to_string(output).expect("read draft");
        assert!(written.contains("\"schema\": \"gentle.gene_group_catalog.v1\""));
        assert!(written.contains("GO:0000381"));
    }

    #[test]
    fn doctor_reports_alias_collision_and_bad_go_mapping() {
        let dir = tempdir().expect("tempdir");
        let path = dir.path().join("groups.json");
        fs::write(
            &path,
            r#"{
              "schema":"gentle.gene_group_catalog.v1",
              "external_resources":[{"id":"gene_ontology","label":"Gene Ontology","namespace":"GO","role":"external_ontology_mapping_namespace"}],
              "groups":[
                {"id":"a","label":"Alpha","aliases":["shared"],"members":[{"symbol":"A"}],"external_mappings":[{"namespace":"GO","id":"bad"}]},
                {"id":"b","label":"Beta","aliases":["shared"],"members":[{"symbol":"B"}]}
              ]
            }"#,
        )
        .expect("write catalog");
        let report = doctor_gene_group_catalog(Some(path.to_string_lossy().as_ref()));
        assert!(report.error_count >= 2, "{report:#?}");
        assert!(
            report
                .sources
                .iter()
                .flat_map(|source| source.errors.iter())
                .any(|error| error.contains("collides"))
        );
        assert!(
            report
                .sources
                .iter()
                .flat_map(|source| source.errors.iter())
                .any(|error| error.contains("malformed GO"))
        );
    }
}
