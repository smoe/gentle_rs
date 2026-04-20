//! CUT&RUN catalog, prepare, status, and projection helpers.
//!
//! This module owns the V1 CUT&RUN processed-evidence workflow:
//! - catalog discovery/loading
//! - prepared processed-asset materialization under `data/cutrun`
//! - status/manifest inspection
//! - projection of prepared BED/BigWig evidence onto genome-anchored regions
//!
//! The implementation intentionally reuses the existing genome BED/BigWig
//! import paths so CUT&RUN evidence follows the same anchored-track semantics
//! as other projected signal data.

use super::*;
use crate::genomes::{BUILTIN_ASSET_ROOT_ENV, PROJECT_ROOT_ENV, SYSTEM_CONFIG_ROOT_ENV};
use reqwest::blocking::Client;
use sha1::{Digest, Sha1};
use std::fs;

const CUTRUN_MANIFEST_FILE_NAME: &str = "manifest.json";
const CUTRUN_CACHE_DIR_ENV: &str = "GENTLE_CUTRUN_CACHE_DIR";

#[derive(Clone, Debug)]
struct CutRunCatalogSourceCandidate {
    scope: &'static str,
    path: PathBuf,
}

#[derive(Clone, Debug)]
struct CutRunCatalog {
    entries: HashMap<String, CutRunCatalogEntry>,
    entry_base_dirs: HashMap<String, PathBuf>,
    catalog_origin_label: String,
    requested_catalog_path: Option<String>,
}

impl CutRunCatalog {
    fn from_json_file(path: &str) -> Result<Self, String> {
        let trimmed = path.trim();
        if trimmed.is_empty() {
            return Err("CUT&RUN catalog path must not be empty".to_string());
        }
        let source = PathBuf::from(trimmed);
        let requested_catalog_path = Some(trimmed.to_string());
        if source.is_dir() {
            Self::from_sources(
                &[CutRunCatalogSourceCandidate {
                    scope: "explicit",
                    path: source,
                }],
                trimmed.to_string(),
                requested_catalog_path,
            )
        } else {
            Self::from_sources(
                &[CutRunCatalogSourceCandidate {
                    scope: "explicit",
                    path: source,
                }],
                trimmed.to_string(),
                requested_catalog_path,
            )
        }
    }

    fn from_default_discovery() -> Result<Self, String> {
        let sources = discovered_cutrun_catalog_sources()?;
        let origin_label = cutrun_discovery_origin_label(&sources);
        Self::from_sources(&sources, origin_label, None)
    }

    fn from_sources(
        sources: &[CutRunCatalogSourceCandidate],
        origin_label: String,
        requested_catalog_path: Option<String>,
    ) -> Result<Self, String> {
        if sources.is_empty() {
            return Err("No CUT&RUN catalog sources were provided".to_string());
        }
        let mut entries: HashMap<String, CutRunCatalogEntry> = HashMap::new();
        let mut entry_base_dirs: HashMap<String, PathBuf> = HashMap::new();
        let mut source_files: HashMap<String, PathBuf> = HashMap::new();
        for source in sources {
            Self::merge_catalog_source(
                &source.path,
                &source.path.display().to_string(),
                &mut entries,
                &mut entry_base_dirs,
                &mut source_files,
            )?;
        }
        Ok(Self {
            entries,
            entry_base_dirs,
            catalog_origin_label: origin_label,
            requested_catalog_path,
        })
    }

    fn merge_catalog_source(
        path: &Path,
        display_path: &str,
        entries: &mut HashMap<String, CutRunCatalogEntry>,
        entry_base_dirs: &mut HashMap<String, PathBuf>,
        source_files: &mut HashMap<String, PathBuf>,
    ) -> Result<(), String> {
        let metadata = fs::metadata(path)
            .map_err(|e| format!("Could not read CUT&RUN catalog '{display_path}': {e}"))?;
        if metadata.is_dir() {
            Self::merge_json_directory(path, display_path, entries, entry_base_dirs, source_files)
        } else {
            Self::merge_single_json_file(path, display_path, entries, entry_base_dirs, source_files)
        }
    }

    fn merge_single_json_file(
        path: &Path,
        display_path: &str,
        entries: &mut HashMap<String, CutRunCatalogEntry>,
        entry_base_dirs: &mut HashMap<String, PathBuf>,
        source_files: &mut HashMap<String, PathBuf>,
    ) -> Result<(), String> {
        let text = fs::read_to_string(path)
            .map_err(|e| format!("Could not read CUT&RUN catalog '{display_path}': {e}"))?;
        let mut parsed: HashMap<String, CutRunCatalogEntry> = serde_json::from_str(&text)
            .map_err(|e| format!("Could not parse CUT&RUN catalog '{display_path}': {e}"))?;
        let base = path
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        let mut keys: Vec<String> = parsed.keys().cloned().collect();
        keys.sort_unstable();
        for key in keys {
            if let Some(previous) = source_files.get(&key) {
                return Err(format!(
                    "Duplicate CUT&RUN dataset '{}' found in '{}' and '{}'",
                    key,
                    previous.display(),
                    path.display()
                ));
            }
            let entry = parsed.remove(&key).ok_or_else(|| {
                format!(
                    "CUT&RUN dataset '{}' disappeared while loading '{}'",
                    key, display_path
                )
            })?;
            source_files.insert(key.clone(), path.to_path_buf());
            entry_base_dirs.insert(key.clone(), base.clone());
            entries.insert(key, entry);
        }
        Ok(())
    }

    fn merge_json_directory(
        path: &Path,
        display_path: &str,
        entries: &mut HashMap<String, CutRunCatalogEntry>,
        entry_base_dirs: &mut HashMap<String, PathBuf>,
        source_files: &mut HashMap<String, PathBuf>,
    ) -> Result<(), String> {
        let read_dir = fs::read_dir(path).map_err(|e| {
            format!("Could not read CUT&RUN catalog directory '{display_path}': {e}")
        })?;
        let mut json_files: Vec<PathBuf> = vec![];
        for row in read_dir {
            let row = row.map_err(|e| {
                format!("Could not iterate CUT&RUN catalog directory '{display_path}': {e}")
            })?;
            let file_type = row.file_type().map_err(|e| {
                format!(
                    "Could not inspect CUT&RUN catalog directory entry '{}' in '{}': {e}",
                    row.path().display(),
                    display_path
                )
            })?;
            if !file_type.is_file() {
                continue;
            }
            let file_path = row.path();
            if file_path
                .extension()
                .and_then(|value| value.to_str())
                .map(|value| value.eq_ignore_ascii_case("json"))
                .unwrap_or(false)
            {
                json_files.push(file_path);
            }
        }
        json_files.sort();
        if json_files.is_empty() {
            return Err(format!(
                "CUT&RUN catalog directory '{}' does not contain any .json files",
                display_path
            ));
        }
        for json_file in json_files {
            Self::merge_single_json_file(
                &json_file,
                &json_file.display().to_string(),
                entries,
                entry_base_dirs,
                source_files,
            )?;
        }
        Ok(())
    }

    fn list_entries(&self, filter: Option<&str>) -> Vec<CutRunCatalogListEntry> {
        let mut dataset_ids: Vec<String> = self.entries.keys().cloned().collect();
        dataset_ids.sort_unstable();
        dataset_ids
            .into_iter()
            .filter_map(|dataset_id| {
                let entry = self.entries.get(&dataset_id)?;
                if !Self::entry_matches_filter(&dataset_id, entry, filter) {
                    return None;
                }
                Some(CutRunCatalogListEntry {
                    dataset_id,
                    description: entry.description.clone(),
                    summary: entry.summary.clone(),
                    aliases: dedup_sorted_strings(entry.aliases.clone()),
                    tags: dedup_sorted_strings(entry.tags.clone()),
                    search_terms: dedup_sorted_strings(entry.search_terms.clone()),
                    species: entry.species.clone(),
                    assembly_label: entry.assembly_label.clone(),
                    target_factor: entry.target_factor.clone(),
                    sample_label: entry.sample_label.clone(),
                    tissue_or_cell_type: entry.tissue_or_cell_type.clone(),
                    condition: entry.condition.clone(),
                    replicate: entry.replicate.clone(),
                    assay_kind: entry.assay_kind.clone(),
                    supported_reference_genome_ids: dedup_sorted_strings(
                        entry.supported_reference_genome_ids.clone(),
                    ),
                    provider: entry.provider.clone(),
                    source_accession: entry.source_accession.clone(),
                    reference_url: entry.reference_url.clone(),
                    has_peaks_asset: cutrun_entry_source(
                        entry.peaks_local.as_deref(),
                        entry.peaks_remote.as_deref(),
                    )
                    .is_some(),
                    has_signal_asset: cutrun_entry_source(
                        entry.signal_local.as_deref(),
                        entry.signal_remote.as_deref(),
                    )
                    .is_some(),
                    has_raw_reads: cutrun_entry_source(
                        entry.reads_r1_local.as_deref(),
                        entry.reads_r1_remote.as_deref(),
                    )
                    .is_some()
                        || cutrun_entry_source(
                            entry.reads_r2_local.as_deref(),
                            entry.reads_r2_remote.as_deref(),
                        )
                        .is_some(),
                    read_layout: entry.read_layout,
                })
            })
            .collect()
    }

    fn entry_matches_filter(
        dataset_id: &str,
        entry: &CutRunCatalogEntry,
        filter: Option<&str>,
    ) -> bool {
        let Some(filter) = filter.map(str::trim).filter(|value| !value.is_empty()) else {
            return true;
        };
        let needle = filter.to_ascii_lowercase();
        Self::entry_search_terms(dataset_id, entry)
            .into_iter()
            .map(|value| value.to_ascii_lowercase())
            .any(|value| value.contains(&needle))
    }

    fn entry_search_terms(dataset_id: &str, entry: &CutRunCatalogEntry) -> Vec<String> {
        let mut values = vec![dataset_id.to_string()];
        values.extend(entry.aliases.clone());
        values.extend(entry.tags.clone());
        values.extend(entry.search_terms.clone());
        for value in [
            entry.description.as_ref(),
            entry.summary.as_ref(),
            entry.species.as_ref(),
            entry.assembly_label.as_ref(),
            entry.target_factor.as_ref(),
            entry.sample_label.as_ref(),
            entry.tissue_or_cell_type.as_ref(),
            entry.condition.as_ref(),
            entry.replicate.as_ref(),
            entry.assay_kind.as_ref(),
            entry.provider.as_ref(),
            entry.source_accession.as_ref(),
            entry.reference_url.as_ref(),
        ] {
            if let Some(value) = value {
                values.push(value.clone());
            }
        }
        values.extend(entry.supported_reference_genome_ids.clone());
        values
    }

    fn resolve_entry(
        &self,
        dataset_query: &str,
    ) -> Result<(String, &CutRunCatalogEntry, PathBuf), String> {
        let needle = dataset_query.trim().to_ascii_lowercase();
        if needle.is_empty() {
            return Err("CUT&RUN dataset id must not be empty".to_string());
        }
        for (dataset_id, entry) in &self.entries {
            if dataset_id.eq_ignore_ascii_case(&needle)
                || entry
                    .aliases
                    .iter()
                    .any(|alias| alias.eq_ignore_ascii_case(&needle))
            {
                let base_dir = self
                    .entry_base_dirs
                    .get(dataset_id)
                    .cloned()
                    .unwrap_or_else(|| PathBuf::from("."));
                return Ok((dataset_id.clone(), entry, base_dir));
            }
        }
        Err(format!(
            "CUT&RUN dataset '{}' was not found in {}",
            dataset_query, self.catalog_origin_label
        ))
    }

    fn cache_dir_for_entry(
        &self,
        entry: &CutRunCatalogEntry,
        cache_dir_override: Option<&str>,
    ) -> PathBuf {
        cache_dir_override
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(PathBuf::from)
            .or_else(configured_cutrun_cache_dir)
            .or_else(|| {
                entry.cache_dir.as_deref().map(str::trim).and_then(|value| {
                    if value.is_empty() {
                        None
                    } else {
                        Some(PathBuf::from(value))
                    }
                })
            })
            .unwrap_or_else(|| PathBuf::from(DEFAULT_CUTRUN_CACHE_DIR))
    }
}

fn configured_builtin_asset_root() -> PathBuf {
    std::env::var_os(BUILTIN_ASSET_ROOT_ENV)
        .filter(|value| !value.is_empty())
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")))
}

fn configured_system_config_root() -> PathBuf {
    std::env::var_os(SYSTEM_CONFIG_ROOT_ENV)
        .filter(|value| !value.is_empty())
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/etc/gentle"))
}

fn configured_user_config_root() -> Option<PathBuf> {
    if let Some(root) = std::env::var_os("XDG_CONFIG_HOME").filter(|value| !value.is_empty()) {
        return Some(PathBuf::from(root).join("gentle"));
    }
    std::env::var_os("HOME")
        .filter(|value| !value.is_empty())
        .map(|home| PathBuf::from(home).join(".config").join("gentle"))
}

fn discover_project_root_from_cwd() -> Option<PathBuf> {
    let cwd = std::env::current_dir().ok()?;
    let mut git_root: Option<PathBuf> = None;
    let mut cargo_root: Option<PathBuf> = None;
    for ancestor in cwd.ancestors() {
        if ancestor.join(".gentle").exists() {
            return Some(ancestor.to_path_buf());
        }
        if git_root.is_none() && ancestor.join(".git").exists() {
            git_root = Some(ancestor.to_path_buf());
        }
        if cargo_root.is_none() && ancestor.join("Cargo.toml").is_file() {
            cargo_root = Some(ancestor.to_path_buf());
        }
    }
    git_root.or(cargo_root)
}

fn configured_project_root() -> Option<PathBuf> {
    if let Some(root) = std::env::var_os(PROJECT_ROOT_ENV).filter(|value| !value.is_empty()) {
        return Some(PathBuf::from(root));
    }
    discover_project_root_from_cwd()
}

fn cutrun_catalog_discovery_candidates() -> Vec<CutRunCatalogSourceCandidate> {
    let mut candidates = vec![];
    let built_in_root = configured_builtin_asset_root();
    candidates.push(CutRunCatalogSourceCandidate {
        scope: "built-in",
        path: built_in_root.join(DEFAULT_CUTRUN_CATALOG_PATH),
    });
    candidates.push(CutRunCatalogSourceCandidate {
        scope: "built-in",
        path: built_in_root.join("assets").join("cutrun.d"),
    });

    let system_root = configured_system_config_root();
    candidates.push(CutRunCatalogSourceCandidate {
        scope: "system",
        path: system_root.join("catalogs").join("cutrun.json"),
    });
    candidates.push(CutRunCatalogSourceCandidate {
        scope: "system",
        path: system_root.join("catalogs").join("cutrun.d"),
    });

    if let Some(user_root) = configured_user_config_root() {
        candidates.push(CutRunCatalogSourceCandidate {
            scope: "user",
            path: user_root.join("catalogs").join("cutrun.json"),
        });
        candidates.push(CutRunCatalogSourceCandidate {
            scope: "user",
            path: user_root.join("catalogs").join("cutrun.d"),
        });
    }

    if let Some(project_root) = configured_project_root() {
        candidates.push(CutRunCatalogSourceCandidate {
            scope: "project",
            path: project_root
                .join(".gentle")
                .join("catalogs")
                .join("cutrun.json"),
        });
        candidates.push(CutRunCatalogSourceCandidate {
            scope: "project",
            path: project_root
                .join(".gentle")
                .join("catalogs")
                .join("cutrun.d"),
        });
    }

    candidates
}

fn discovered_cutrun_catalog_sources() -> Result<Vec<CutRunCatalogSourceCandidate>, String> {
    let mut existing = vec![];
    let mut searched = vec![];
    let mut seen_paths = BTreeSet::new();
    for candidate in cutrun_catalog_discovery_candidates() {
        let key = candidate.path.to_string_lossy().to_string();
        if !seen_paths.insert(key.clone()) {
            continue;
        }
        searched.push(format!("{}={}", candidate.scope, candidate.path.display()));
        if candidate.path.exists() {
            existing.push(candidate);
        }
    }
    if existing.is_empty() {
        return Err(format!(
            "Could not find any CUT&RUN catalog sources for default CUT&RUN catalog discovery. Searched: {}",
            searched.join(", ")
        ));
    }
    Ok(existing)
}

fn cutrun_discovery_origin_label(sources: &[CutRunCatalogSourceCandidate]) -> String {
    if sources.len() == 1 {
        return sources[0].path.display().to_string();
    }
    let scopes = sources
        .iter()
        .map(|source| format!("{}={}", source.scope, source.path.display()))
        .collect::<Vec<_>>();
    format!("default CUT&RUN catalog discovery [{}]", scopes.join(", "))
}

fn cutrun_entry_source(local: Option<&str>, remote: Option<&str>) -> Option<String> {
    local
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string())
        .or_else(|| {
            remote
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string())
        })
}

fn dedup_sorted_strings(values: Vec<String>) -> Vec<String> {
    let mut set = BTreeSet::new();
    for value in values {
        let trimmed = value.trim();
        if !trimmed.is_empty() {
            set.insert(trimmed.to_string());
        }
    }
    set.into_iter().collect()
}

fn configured_cutrun_cache_dir() -> Option<PathBuf> {
    std::env::var_os(CUTRUN_CACHE_DIR_ENV)
        .filter(|value| !value.is_empty())
        .map(PathBuf::from)
}

fn sanitize_for_path(value: &str) -> String {
    let mut out = String::new();
    let mut last_was_sep = false;
    for ch in value.chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch.to_ascii_lowercase());
            last_was_sep = false;
        } else if !last_was_sep {
            out.push('_');
            last_was_sep = true;
        }
    }
    let trimmed = out.trim_matches('_').to_string();
    if trimmed.is_empty() {
        "cutrun_dataset".to_string()
    } else {
        trimmed
    }
}

fn source_file_name(source: &str, fallback: &str) -> String {
    let without_query = source.split('?').next().unwrap_or(source);
    let localish = without_query
        .strip_prefix("file://")
        .unwrap_or(without_query);
    Path::new(localish)
        .file_name()
        .and_then(|value| value.to_str())
        .filter(|value| !value.trim().is_empty())
        .map(|value| value.to_string())
        .unwrap_or_else(|| fallback.to_string())
}

fn is_http_source(source: &str) -> bool {
    let lower = source.trim().to_ascii_lowercase();
    lower.starts_with("http://") || lower.starts_with("https://")
}

fn resolved_local_source_path(base_dir: &Path, source: &str) -> PathBuf {
    let trimmed = source.trim();
    if let Some(rest) = trimmed.strip_prefix("file://") {
        PathBuf::from(rest)
    } else {
        let path = PathBuf::from(trimmed);
        if path.is_absolute() {
            path
        } else {
            base_dir.join(path)
        }
    }
}

fn file_size_bytes(path: &Path) -> Result<u64, EngineError> {
    fs::metadata(path)
        .map(|metadata| metadata.len())
        .map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not inspect '{}': {e}", path.display()),
        })
}

fn checksum_sha1(path: &Path) -> Result<String, EngineError> {
    let mut file = File::open(path).map_err(|e| EngineError {
        code: ErrorCode::Io,
        message: format!(
            "Could not open '{}' to compute checksum: {e}",
            path.display()
        ),
    })?;
    let mut hasher = Sha1::new();
    let mut buffer = [0u8; 16 * 1024];
    loop {
        let read = file.read(&mut buffer).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not read '{}' while computing checksum: {e}",
                path.display()
            ),
        })?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(format!("sha1:{:x}", hasher.finalize()))
}

impl GentleEngine {
    fn open_cutrun_catalog(catalog_path: Option<&str>) -> Result<CutRunCatalog, EngineError> {
        let trimmed = catalog_path
            .map(str::trim)
            .filter(|value| !value.is_empty());
        let catalog = match trimmed {
            Some(path) => CutRunCatalog::from_json_file(path),
            None => CutRunCatalog::from_default_discovery(),
        }
        .map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: e,
        })?;
        Ok(catalog)
    }

    pub fn list_cutrun_datasets(
        filter: Option<&str>,
        catalog_path: Option<&str>,
    ) -> Result<CutRunDatasetListReport, EngineError> {
        let catalog = Self::open_cutrun_catalog(catalog_path)?;
        let datasets = catalog.list_entries(filter);
        Ok(CutRunDatasetListReport {
            schema: CUTRUN_DATASET_LIST_SCHEMA.to_string(),
            catalog_origin_label: catalog.catalog_origin_label.clone(),
            requested_catalog_path: catalog.requested_catalog_path.clone(),
            filter: filter
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            dataset_count: datasets.len(),
            datasets,
        })
    }

    pub fn show_cutrun_dataset_status(
        &self,
        dataset_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<CutRunDatasetStatus, EngineError> {
        let catalog = Self::open_cutrun_catalog(catalog_path)?;
        self.cutrun_dataset_status_from_catalog(&catalog, dataset_id, cache_dir)
    }

    pub fn prepare_cutrun_dataset(
        &self,
        dataset_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<CutRunDatasetStatus, EngineError> {
        let catalog = Self::open_cutrun_catalog(catalog_path)?;
        let (resolved_dataset_id, entry, entry_base_dir) =
            catalog.resolve_entry(dataset_id).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: e,
            })?;
        let effective_cache_dir = catalog.cache_dir_for_entry(entry, cache_dir);
        let install_dir = effective_cache_dir.join(sanitize_for_path(&resolved_dataset_id));
        fs::create_dir_all(&install_dir).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create CUT&RUN install directory '{}': {e}",
                install_dir.display()
            ),
        })?;

        let peaks_source =
            cutrun_entry_source(entry.peaks_local.as_deref(), entry.peaks_remote.as_deref());
        let signal_source = cutrun_entry_source(
            entry.signal_local.as_deref(),
            entry.signal_remote.as_deref(),
        );
        if peaks_source.is_none() && signal_source.is_none() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "CUT&RUN dataset '{}' does not declare any processed peaks/signal assets",
                    resolved_dataset_id
                ),
            });
        }

        let peaks = peaks_source
            .as_deref()
            .map(|source| {
                self.materialize_cutrun_asset(source, &entry_base_dir, &install_dir, "peaks.bed")
            })
            .transpose()?;
        let signal = signal_source
            .as_deref()
            .map(|source| {
                self.materialize_cutrun_asset(
                    source,
                    &entry_base_dir,
                    &install_dir,
                    "signal.bigwig",
                )
            })
            .transpose()?;

        let manifest = CutRunPreparedManifest {
            schema: CUTRUN_PREPARED_MANIFEST_SCHEMA.to_string(),
            dataset_id: resolved_dataset_id.clone(),
            prepared_at_unix_ms: Self::now_unix_ms(),
            catalog_origin_label: catalog.catalog_origin_label.clone(),
            install_dir: install_dir.display().to_string(),
            peaks,
            signal,
        };
        self.write_cutrun_prepared_manifest(
            &install_dir.join(CUTRUN_MANIFEST_FILE_NAME),
            &manifest,
        )?;
        self.cutrun_dataset_status_from_catalog(&catalog, &resolved_dataset_id, cache_dir)
    }

    pub fn project_cutrun_dataset(
        &mut self,
        seq_id: &str,
        dataset_id: &str,
        include_peaks: bool,
        include_signal: bool,
        clear_existing: bool,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<CutRunDatasetProjectionReport, EngineError> {
        if !include_peaks && !include_signal {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "ProjectCutRunDataset requires at least one of include_peaks/include_signal"
                        .to_string(),
            });
        }
        let catalog = Self::open_cutrun_catalog(catalog_path)?;
        let (resolved_dataset_id, entry, _) =
            catalog.resolve_entry(dataset_id).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: e,
            })?;
        let status =
            self.cutrun_dataset_status_from_catalog(&catalog, &resolved_dataset_id, cache_dir)?;
        let anchor = self.latest_genome_anchor_for_seq(seq_id)?;
        if !entry.supported_reference_genome_ids.is_empty()
            && !entry
                .supported_reference_genome_ids
                .iter()
                .any(|supported| supported.eq_ignore_ascii_case(&anchor.genome_id))
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "CUT&RUN dataset '{}' supports [{}], but sequence '{}' is anchored to '{}'",
                    resolved_dataset_id,
                    entry.supported_reference_genome_ids.join(", "),
                    seq_id,
                    anchor.genome_id
                ),
            });
        }

        let peak_path = status.peaks.local_path.clone();
        let signal_path = status.signal.local_path.clone();
        let mut report = CutRunDatasetProjectionReport {
            schema: CUTRUN_DATASET_PROJECTION_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            dataset_id: resolved_dataset_id.clone(),
            include_peaks,
            include_signal,
            clear_existing,
            projected_peak_features: 0,
            projected_signal_features: 0,
            peak_track_name: include_peaks.then(|| format!("{resolved_dataset_id} peaks")),
            signal_track_name: include_signal.then(|| format!("{resolved_dataset_id} signal")),
            warnings: vec![],
        };

        let mut cleared_once = false;
        let mut projected_any = false;
        {
            let dna = self
                .state
                .sequences
                .get_mut(seq_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!("Sequence '{}' was not found", seq_id),
                })?;

            if include_peaks {
                if let Some(path) = peak_path.as_deref() {
                    let import_report = Self::import_genome_bed_track(
                        dna,
                        &anchor,
                        path,
                        report.peak_track_name.as_deref(),
                        None,
                        None,
                        clear_existing && !cleared_once,
                        None,
                    )?;
                    report.projected_peak_features = import_report.imported_features;
                    report.warnings.extend(import_report.warnings);
                    cleared_once |= clear_existing;
                    projected_any |=
                        import_report.imported_features > 0 || import_report.parsed_records > 0;
                } else {
                    report.warnings.push(format!(
                        "CUT&RUN peaks asset for '{}' is not prepared",
                        resolved_dataset_id
                    ));
                }
            }
            if include_signal {
                if let Some(path) = signal_path.as_deref() {
                    let import_report = Self::import_genome_bigwig_track(
                        dna,
                        &anchor,
                        path,
                        report.signal_track_name.as_deref(),
                        None,
                        None,
                        clear_existing && !cleared_once,
                        None,
                    )?;
                    report.projected_signal_features = import_report.imported_features;
                    report.warnings.extend(import_report.warnings);
                    projected_any |=
                        import_report.imported_features > 0 || import_report.parsed_records > 0;
                } else {
                    report.warnings.push(format!(
                        "CUT&RUN signal asset for '{}' is not prepared",
                        resolved_dataset_id
                    ));
                }
            }
        }

        if !projected_any {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "No prepared CUT&RUN peaks/signal assets were available to project for '{}'",
                    resolved_dataset_id
                ),
            });
        }
        Ok(report)
    }

    fn cutrun_dataset_status_from_catalog(
        &self,
        catalog: &CutRunCatalog,
        dataset_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<CutRunDatasetStatus, EngineError> {
        let (resolved_dataset_id, entry, _) =
            catalog.resolve_entry(dataset_id).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: e,
            })?;
        let effective_cache_dir = catalog.cache_dir_for_entry(entry, cache_dir);
        let install_dir = effective_cache_dir.join(sanitize_for_path(&resolved_dataset_id));
        let manifest_path = install_dir.join(CUTRUN_MANIFEST_FILE_NAME);
        let manifest = if manifest_path.exists() {
            Some(self.read_cutrun_prepared_manifest(&manifest_path)?)
        } else {
            None
        };
        let peaks_source =
            cutrun_entry_source(entry.peaks_local.as_deref(), entry.peaks_remote.as_deref());
        let signal_source = cutrun_entry_source(
            entry.signal_local.as_deref(),
            entry.signal_remote.as_deref(),
        );
        let peaks_status = Self::cutrun_asset_status_from_manifest(
            "peaks",
            peaks_source,
            manifest.as_ref().and_then(|m| m.peaks.as_ref()),
        );
        let signal_status = Self::cutrun_asset_status_from_manifest(
            "signal",
            signal_source,
            manifest.as_ref().and_then(|m| m.signal.as_ref()),
        );
        let configured_assets =
            usize::from(peaks_status.configured) + usize::from(signal_status.configured);
        let prepared_assets =
            usize::from(peaks_status.prepared) + usize::from(signal_status.prepared);
        let prepared = configured_assets > 0 && configured_assets == prepared_assets;
        Ok(CutRunDatasetStatus {
            schema: CUTRUN_DATASET_STATUS_SCHEMA.to_string(),
            dataset_id: resolved_dataset_id,
            catalog_origin_label: catalog.catalog_origin_label.clone(),
            requested_catalog_path: catalog.requested_catalog_path.clone(),
            effective_cache_dir: effective_cache_dir.display().to_string(),
            install_dir: install_dir.display().to_string(),
            prepared,
            description: entry.description.clone(),
            summary: entry.summary.clone(),
            species: entry.species.clone(),
            assembly_label: entry.assembly_label.clone(),
            target_factor: entry.target_factor.clone(),
            sample_label: entry.sample_label.clone(),
            tissue_or_cell_type: entry.tissue_or_cell_type.clone(),
            condition: entry.condition.clone(),
            replicate: entry.replicate.clone(),
            assay_kind: entry.assay_kind.clone(),
            supported_reference_genome_ids: dedup_sorted_strings(
                entry.supported_reference_genome_ids.clone(),
            ),
            provider: entry.provider.clone(),
            source_accession: entry.source_accession.clone(),
            reference_url: entry.reference_url.clone(),
            read_layout: entry.read_layout,
            peaks: peaks_status,
            signal: signal_status,
            manifest_path: manifest_path
                .exists()
                .then(|| manifest_path.display().to_string()),
            manifest,
            warnings: vec![],
        })
    }

    fn cutrun_asset_status_from_manifest(
        _asset_kind: &str,
        source: Option<String>,
        manifest: Option<&CutRunPreparedAssetManifest>,
    ) -> CutRunPreparedAssetStatus {
        if let Some(manifest) = manifest {
            let local_path = PathBuf::from(&manifest.local_path);
            return CutRunPreparedAssetStatus {
                configured: source.is_some(),
                source,
                prepared: local_path.is_file(),
                local_path: Some(manifest.local_path.clone()),
                file_name: Some(manifest.file_name.clone()),
                file_size_bytes: Some(manifest.file_size_bytes),
                checksum_sha1: Some(manifest.checksum_sha1.clone()),
            };
        }
        CutRunPreparedAssetStatus {
            configured: source.is_some(),
            source,
            prepared: false,
            local_path: None,
            file_name: None,
            file_size_bytes: None,
            checksum_sha1: None,
        }
    }

    fn write_cutrun_prepared_manifest(
        &self,
        path: &Path,
        manifest: &CutRunPreparedManifest,
    ) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(manifest).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize CUT&RUN prepared manifest '{}': {e}",
                path.display()
            ),
        })?;
        fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write CUT&RUN prepared manifest '{}': {e}",
                path.display()
            ),
        })
    }

    fn read_cutrun_prepared_manifest(
        &self,
        path: &Path,
    ) -> Result<CutRunPreparedManifest, EngineError> {
        let text = fs::read_to_string(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not read CUT&RUN prepared manifest '{}': {e}",
                path.display()
            ),
        })?;
        let manifest =
            serde_json::from_str::<CutRunPreparedManifest>(&text).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not parse CUT&RUN prepared manifest '{}': {e}",
                    path.display()
                ),
            })?;
        if manifest.schema != CUTRUN_PREPARED_MANIFEST_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported CUT&RUN prepared manifest schema '{}' in '{}'",
                    manifest.schema,
                    path.display()
                ),
            });
        }
        Ok(manifest)
    }

    fn materialize_cutrun_asset(
        &self,
        source: &str,
        entry_base_dir: &Path,
        install_dir: &Path,
        fallback_name: &str,
    ) -> Result<CutRunPreparedAssetManifest, EngineError> {
        let file_name = source_file_name(source, fallback_name);
        let destination = install_dir.join(&file_name);
        if is_http_source(source) {
            let client = Client::builder().build().map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not build HTTP client for CUT&RUN prepare: {e}"),
            })?;
            let mut response = client.get(source).send().map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not download CUT&RUN asset '{}': {e}", source),
            })?;
            if !response.status().is_success() {
                return Err(EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not download CUT&RUN asset '{}': HTTP {}",
                        source,
                        response.status()
                    ),
                });
            }
            let mut output = File::create(&destination).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not create CUT&RUN asset destination '{}': {e}",
                    destination.display()
                ),
            })?;
            std::io::copy(&mut response, &mut output).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not write downloaded CUT&RUN asset '{}' to '{}': {e}",
                    source,
                    destination.display()
                ),
            })?;
        } else {
            let resolved_source = resolved_local_source_path(entry_base_dir, source);
            if !resolved_source.is_file() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "CUT&RUN asset source '{}' does not exist as a file",
                        resolved_source.display()
                    ),
                });
            }
            fs::copy(&resolved_source, &destination).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not copy CUT&RUN asset '{}' to '{}': {e}",
                    resolved_source.display(),
                    destination.display()
                ),
            })?;
        }
        Ok(CutRunPreparedAssetManifest {
            source: source.to_string(),
            local_path: destination.display().to_string(),
            file_name,
            file_size_bytes: file_size_bytes(&destination)?,
            checksum_sha1: checksum_sha1(&destination)?,
        })
    }
}
