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
use crate::genomes::{
    BUILTIN_ASSET_ROOT_ENV, GenomeChromosomeRecord, PROJECT_ROOT_ENV, SYSTEM_CONFIG_ROOT_ENV,
};
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

#[derive(Clone, Debug)]
struct ParsedCutRunInputRecord {
    record_index: usize,
    header_id: String,
    normalized_read_id: String,
    sequence: Vec<u8>,
}

#[derive(Clone, Debug)]
struct CutRunInputUnit {
    normalized_read_id: String,
    r1: Option<ParsedCutRunInputRecord>,
    r2: Option<ParsedCutRunInputRecord>,
}

#[derive(Clone, Debug)]
struct CutRunReferenceWindow {
    genome_id: String,
    chromosome: String,
    window_start_1based: usize,
    window_end_1based: usize,
    roi_local_start_1based: usize,
    roi_local_end_1based: usize,
    orientation: CutRunReadOrientation,
    sequence: Vec<u8>,
    warnings: Vec<String>,
}

impl GentleEngine {
    pub(super) fn normalize_cutrun_report_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "CUT&RUN report_id cannot be empty".to_string(),
            });
        }
        let mut out = String::with_capacity(trimmed.len());
        for ch in trimmed.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                out.push(ch);
            } else {
                out.push('_');
            }
        }
        if out.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "CUT&RUN report_id must contain at least one ASCII letter, digit, '-', '_' or '.'"
                    .to_string(),
            });
        }
        Ok(out)
    }

    pub(super) fn read_cutrun_read_report_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> CutRunReadReportStore {
        let mut store = value
            .cloned()
            .and_then(|raw| serde_json::from_value::<CutRunReadReportStore>(raw).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = CUTRUN_READ_REPORTS_SCHEMA.to_string();
        }
        store
    }

    pub(super) fn read_cutrun_read_report_store(&self) -> CutRunReadReportStore {
        Self::read_cutrun_read_report_store_from_metadata(
            self.state.metadata.get(CUTRUN_READ_REPORTS_METADATA_KEY),
        )
    }

    pub(super) fn write_cutrun_read_report_store(
        &mut self,
        mut store: CutRunReadReportStore,
    ) -> Result<(), EngineError> {
        if store.reports.is_empty() {
            self.state.metadata.remove(CUTRUN_READ_REPORTS_METADATA_KEY);
            return Ok(());
        }
        store.schema = CUTRUN_READ_REPORTS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize CUT&RUN read-report metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(CUTRUN_READ_REPORTS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub(super) fn upsert_cutrun_read_report(
        &mut self,
        report: CutRunReadReport,
    ) -> Result<(), EngineError> {
        let mut store = self.read_cutrun_read_report_store();
        store.reports.insert(report.report_id.clone(), report);
        self.write_cutrun_read_report_store(store)
    }

    fn cutrun_read_report_summary(report: &CutRunReadReport) -> CutRunReadReportSummary {
        CutRunReadReportSummary {
            report_id: report.report_id.clone(),
            seq_id: report.seq_id.clone(),
            generated_at_unix_ms: report.generated_at_unix_ms,
            input_format: report.input_format,
            read_layout: report.read_layout,
            roi_flank_bp: report.roi_flank_bp,
            reference_window_length_bp: report.reference_window_length_bp,
            total_units: report.total_units,
            mapped_units: report.mapped_units,
            fragment_count: report.fragment_count,
            concordant_pair_count: report.concordant_pair_count,
            orphan_unit_count: report
                .orphan_r1_count
                .saturating_add(report.orphan_r2_count),
            mean_read_length_bp: report.mean_read_length_bp,
        }
    }

    pub fn list_cutrun_read_reports(&self, seq_id: Option<&str>) -> Vec<CutRunReadReportSummary> {
        let seq_filter = seq_id.map(str::trim).filter(|value| !value.is_empty());
        let mut rows = self
            .read_cutrun_read_report_store()
            .reports
            .values()
            .filter(|report| {
                seq_filter
                    .map(|needle| report.seq_id == needle)
                    .unwrap_or(true)
            })
            .map(Self::cutrun_read_report_summary)
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            right
                .generated_at_unix_ms
                .cmp(&left.generated_at_unix_ms)
                .then_with(|| left.report_id.cmp(&right.report_id))
        });
        rows
    }

    pub fn get_cutrun_read_report(&self, report_id: &str) -> Result<CutRunReadReport, EngineError> {
        let report_id = Self::normalize_cutrun_report_id(report_id)?;
        self.read_cutrun_read_report_store()
            .reports
            .get(report_id.as_str())
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("CUT&RUN read report '{}' not found", report_id),
            })
    }

    fn normalize_cutrun_input_sequence(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .map(|base| match Self::normalize_nucleotide_base(*base) {
                b'A' | b'C' | b'G' | b'T' | b'N' => Self::normalize_nucleotide_base(*base),
                _ => b'N',
            })
            .collect()
    }

    fn normalize_cutrun_header_id(raw: &str, record_index: usize) -> String {
        raw.trim()
            .split_ascii_whitespace()
            .next()
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .unwrap_or_else(|| format!("record_{}", record_index + 1))
    }

    fn normalize_cutrun_read_id(header_id: &str) -> String {
        let first = header_id
            .trim()
            .split_ascii_whitespace()
            .next()
            .unwrap_or_default()
            .trim_start_matches('@')
            .trim_start_matches('>');
        if let Some(stripped) = first
            .strip_suffix("/1")
            .or_else(|| first.strip_suffix("/2"))
        {
            if !stripped.is_empty() {
                return stripped.to_string();
            }
        }
        if first.is_empty() {
            "unnamed_read".to_string()
        } else {
            first.to_string()
        }
    }

    fn parse_cutrun_fastq_records(path: &str) -> Result<Vec<ParsedCutRunInputRecord>, EngineError> {
        let file = File::open(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not open CUT&RUN FASTQ input '{}': {e}", path),
        })?;
        let lower = path.to_ascii_lowercase();
        let mut reader: Box<dyn BufRead> = if lower.ends_with(".gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(BufReader::new(file))))
        } else {
            Box::new(BufReader::new(file))
        };
        let mut out = Vec::<ParsedCutRunInputRecord>::new();
        let mut header = String::new();
        let mut sequence = String::new();
        let mut plus = String::new();
        let mut qualities = String::new();
        loop {
            header.clear();
            if reader.read_line(&mut header).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read CUT&RUN FASTQ input '{}': {e}", path),
            })? == 0
            {
                break;
            }
            sequence.clear();
            plus.clear();
            qualities.clear();
            if reader.read_line(&mut sequence).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not read CUT&RUN FASTQ sequence line '{}': {e}", path),
            })? == 0
                || reader.read_line(&mut plus).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not read CUT&RUN FASTQ separator line '{}': {e}",
                        path
                    ),
                })? == 0
                || reader.read_line(&mut qualities).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not read CUT&RUN FASTQ quality line '{}': {e}", path),
                })? == 0
            {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("CUT&RUN FASTQ input '{}' ended mid-record", path),
                });
            }
            if !header.starts_with('@') {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "CUT&RUN FASTQ input '{}' has a record without '@' header",
                        path
                    ),
                });
            }
            if !plus.starts_with('+') {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "CUT&RUN FASTQ input '{}' has a record without '+' separator",
                        path
                    ),
                });
            }
            let record_index = out.len();
            let header_text = header[1..].trim();
            let header_id = Self::normalize_cutrun_header_id(header_text, record_index);
            out.push(ParsedCutRunInputRecord {
                record_index,
                header_id: header_id.clone(),
                normalized_read_id: Self::normalize_cutrun_read_id(&header_id),
                sequence: Self::normalize_cutrun_input_sequence(sequence.trim().as_bytes()),
            });
        }
        if out.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("No FASTQ records found in '{}'", path),
            });
        }
        Ok(out)
    }

    fn parse_cutrun_input_records(
        path: &str,
        input_format: CutRunInputFormat,
    ) -> Result<Vec<ParsedCutRunInputRecord>, EngineError> {
        match input_format {
            CutRunInputFormat::Fasta => {
                let mut out = Vec::<ParsedCutRunInputRecord>::new();
                Self::visit_fasta_records_with_offsets(path, &mut |record, _progress| {
                    out.push(ParsedCutRunInputRecord {
                        record_index: record.record_index,
                        header_id: record.header_id.clone(),
                        normalized_read_id: Self::normalize_cutrun_read_id(&record.header_id),
                        sequence: Self::normalize_cutrun_input_sequence(&record.sequence),
                    });
                    Ok(())
                })?;
                Ok(out)
            }
            CutRunInputFormat::Fastq => Self::parse_cutrun_fastq_records(path),
        }
    }

    fn build_cutrun_input_units(
        input_r1: Vec<ParsedCutRunInputRecord>,
        input_r2: Vec<ParsedCutRunInputRecord>,
        read_layout: CutRunReadLayout,
    ) -> Vec<CutRunInputUnit> {
        match read_layout {
            CutRunReadLayout::SingleEnd => input_r1
                .into_iter()
                .map(|r1| CutRunInputUnit {
                    normalized_read_id: r1.normalized_read_id.clone(),
                    r1: Some(r1),
                    r2: None,
                })
                .collect(),
            CutRunReadLayout::PairedEnd => {
                let mut by_r1 = BTreeMap::<String, Vec<ParsedCutRunInputRecord>>::new();
                let mut by_r2 = BTreeMap::<String, Vec<ParsedCutRunInputRecord>>::new();
                for record in input_r1 {
                    by_r1
                        .entry(record.normalized_read_id.clone())
                        .or_default()
                        .push(record);
                }
                for record in input_r2 {
                    by_r2
                        .entry(record.normalized_read_id.clone())
                        .or_default()
                        .push(record);
                }
                let all_ids = by_r1
                    .keys()
                    .chain(by_r2.keys())
                    .cloned()
                    .collect::<BTreeSet<_>>();
                let mut units = Vec::<CutRunInputUnit>::new();
                for normalized_read_id in all_ids {
                    let left = by_r1.remove(&normalized_read_id).unwrap_or_default();
                    let right = by_r2.remove(&normalized_read_id).unwrap_or_default();
                    let pair_count = left.len().max(right.len());
                    for pair_index in 0..pair_count {
                        units.push(CutRunInputUnit {
                            normalized_read_id: normalized_read_id.clone(),
                            r1: left.get(pair_index).cloned(),
                            r2: right.get(pair_index).cloned(),
                        });
                    }
                }
                units
            }
        }
    }

    fn cutrun_match_fraction(left: &[u8], right: &[u8]) -> f64 {
        if left.is_empty() || left.len() != right.len() {
            return 0.0;
        }
        let matches = left
            .iter()
            .zip(right.iter())
            .filter(|(a, b)| a == b)
            .count();
        matches as f64 / left.len() as f64
    }

    fn cutrun_find_chromosome_length(
        catalog: &GenomeCatalog,
        genome_id: &str,
        chromosome: &str,
        cache_dir: Option<&str>,
    ) -> Result<usize, EngineError> {
        let records = catalog
            .list_chromosome_lengths(genome_id, cache_dir)
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not list chromosome lengths for genome '{}' while preparing CUT&RUN ROI window: {}",
                    genome_id, e
                ),
            })?;
        let normalized = Self::normalize_chromosome_alias(chromosome);
        records
            .into_iter()
            .find(|row: &GenomeChromosomeRecord| {
                row.chromosome.eq_ignore_ascii_case(chromosome)
                    || Self::normalize_chromosome_alias(&row.chromosome) == normalized
            })
            .map(|row| row.length_bp)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Chromosome '{}' is not available in prepared genome '{}'",
                    chromosome, genome_id
                ),
            })
    }

    fn build_cutrun_reference_window(
        &self,
        seq_id: &str,
        roi_flank_bp: usize,
    ) -> Result<CutRunReferenceWindow, EngineError> {
        let anchor = self.latest_genome_anchor_for_seq(seq_id)?;
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        let roi_sequence = Self::normalize_cutrun_input_sequence(dna.forward_bytes());
        let (catalog, _) = Self::open_reference_genome_catalog(anchor.catalog_path.as_deref())?;
        let chromosome_length = Self::cutrun_find_chromosome_length(
            &catalog,
            &anchor.genome_id,
            &anchor.chromosome,
            anchor.cache_dir.as_deref(),
        )?;
        let window_start_1based = anchor.start_1based.saturating_sub(roi_flank_bp).max(1);
        let window_end_1based = anchor
            .end_1based
            .saturating_add(roi_flank_bp)
            .min(chromosome_length);
        let forward_window = catalog
            .get_sequence_region_with_cache(
                &anchor.genome_id,
                &anchor.chromosome,
                window_start_1based,
                window_end_1based,
                anchor.cache_dir.as_deref(),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not load CUT&RUN ROI window for '{}': {}", seq_id, e),
            })?;
        let forward_window_bytes = Self::normalize_cutrun_input_sequence(forward_window.as_bytes());
        let anchor_offset_start = anchor.start_1based.saturating_sub(window_start_1based);
        let anchor_offset_end = anchor_offset_start.saturating_add(roi_sequence.len());
        if anchor_offset_end > forward_window_bytes.len() {
            return Err(EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "CUT&RUN ROI window for '{}' does not cover the anchored sequence span",
                    seq_id
                ),
            });
        }
        let forward_anchor = &forward_window_bytes[anchor_offset_start..anchor_offset_end];
        let reverse_anchor = Self::reverse_complement_bytes(forward_anchor);
        let forward_match = Self::cutrun_match_fraction(&roi_sequence, forward_anchor);
        let reverse_match = Self::cutrun_match_fraction(&roi_sequence, &reverse_anchor);
        let reverse_orientation = reverse_match > forward_match;
        let window_len = forward_window_bytes.len();
        let mut warnings = Vec::<String>::new();
        let (sequence, orientation, roi_local_start_1based, roi_local_end_1based) =
            if reverse_orientation {
                if reverse_match < 0.95 {
                    warnings.push(format!(
                        "CUT&RUN ROI window for '{}' matched the reverse-complement anchor weakly (identity={:.3})",
                        seq_id, reverse_match
                    ));
                }
                (
                    Self::reverse_complement_bytes(&forward_window_bytes),
                    CutRunReadOrientation::Reverse,
                    window_len
                        .saturating_sub(anchor_offset_end)
                        .saturating_add(1),
                    window_len.saturating_sub(anchor_offset_start),
                )
            } else {
                if forward_match < 0.95 {
                    warnings.push(format!(
                        "CUT&RUN ROI window for '{}' matched the forward anchor weakly (identity={:.3})",
                        seq_id, forward_match
                    ));
                }
                (
                    forward_window_bytes,
                    CutRunReadOrientation::Forward,
                    anchor_offset_start.saturating_add(1),
                    anchor_offset_end,
                )
            };
        Ok(CutRunReferenceWindow {
            genome_id: anchor.genome_id,
            chromosome: anchor.chromosome,
            window_start_1based,
            window_end_1based,
            roi_local_start_1based,
            roi_local_end_1based,
            orientation,
            sequence,
            warnings,
        })
    }

    fn cutrun_local_pos_to_genomic(
        window: &CutRunReferenceWindow,
        local_pos_1based: usize,
    ) -> usize {
        match window.orientation {
            CutRunReadOrientation::Forward => window
                .window_start_1based
                .saturating_add(local_pos_1based.saturating_sub(1)),
            CutRunReadOrientation::Reverse => window
                .window_end_1based
                .saturating_sub(local_pos_1based.saturating_sub(1)),
        }
    }

    fn cutrun_local_span_to_genomic(
        window: &CutRunReferenceWindow,
        local_start_1based: usize,
        local_end_1based: usize,
    ) -> (usize, usize) {
        let a = Self::cutrun_local_pos_to_genomic(window, local_start_1based);
        let b = Self::cutrun_local_pos_to_genomic(window, local_end_1based);
        (a.min(b), a.max(b))
    }

    fn build_cutrun_seed_index(reference: &[u8], kmer_len: usize) -> HashMap<Vec<u8>, Vec<usize>> {
        if kmer_len == 0 || reference.len() < kmer_len {
            return HashMap::new();
        }
        let mut index = HashMap::<Vec<u8>, Vec<usize>>::new();
        for start in 0..=reference.len() - kmer_len {
            index
                .entry(reference[start..start + kmer_len].to_vec())
                .or_default()
                .push(start);
        }
        index
    }

    fn cutrun_candidate_starts(
        read: &[u8],
        seed_filter: &CutRunSeedFilterConfig,
        seed_index: &HashMap<Vec<u8>, Vec<usize>>,
        reference_len: usize,
    ) -> Vec<(usize, usize)> {
        if read.is_empty() || read.len() > reference_len || seed_filter.kmer_len == 0 {
            return vec![];
        }
        if read.len() < seed_filter.kmer_len {
            return vec![];
        }
        let mut counts = HashMap::<usize, usize>::new();
        for offset in 0..=read.len() - seed_filter.kmer_len {
            let key = read[offset..offset + seed_filter.kmer_len].to_vec();
            let Some(reference_positions) = seed_index.get(&key) else {
                continue;
            };
            for reference_position in reference_positions {
                if *reference_position < offset {
                    continue;
                }
                let start = reference_position - offset;
                if start + read.len() > reference_len {
                    continue;
                }
                *counts.entry(start).or_insert(0) += 1;
            }
        }
        let mut out = counts
            .into_iter()
            .filter(|(_, seed_matches)| *seed_matches >= seed_filter.min_seed_matches.max(1))
            .collect::<Vec<_>>();
        out.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
        out
    }

    fn place_cutrun_read_against_window(
        read: &[u8],
        window: &CutRunReferenceWindow,
        seed_filter: &CutRunSeedFilterConfig,
        align_config: &CutRunAlignConfig,
    ) -> Option<CutRunReadPlacement> {
        if read.is_empty() || read.len() > window.sequence.len() {
            return None;
        }
        let normalized = Self::normalize_cutrun_input_sequence(read);
        let reverse = Self::reverse_complement_bytes(&normalized);
        let seed_index = Self::build_cutrun_seed_index(&window.sequence, seed_filter.kmer_len);
        let candidates = [
            (CutRunReadOrientation::Forward, normalized),
            (CutRunReadOrientation::Reverse, reverse),
        ];
        let mut best: Option<CutRunReadPlacement> = None;
        for (orientation, oriented_read) in candidates {
            let mut candidate_starts = Self::cutrun_candidate_starts(
                &oriented_read,
                seed_filter,
                &seed_index,
                window.sequence.len(),
            );
            if candidate_starts.is_empty() {
                candidate_starts = (0..=window.sequence.len() - oriented_read.len())
                    .map(|start| (start, 0usize))
                    .collect();
            }
            for (start_0based, seed_matches) in candidate_starts {
                let end_0based_exclusive = start_0based + oriented_read.len();
                let target = &window.sequence[start_0based..end_0based_exclusive];
                let mismatches = oriented_read
                    .iter()
                    .zip(target.iter())
                    .filter(|(a, b)| a != b)
                    .count();
                if mismatches > align_config.max_mismatches {
                    continue;
                }
                let identity_fraction = 1.0 - mismatches as f64 / oriented_read.len().max(1) as f64;
                if identity_fraction + f64::EPSILON < align_config.min_identity_fraction {
                    continue;
                }
                let local_start_1based = start_0based + 1;
                let local_end_1based = end_0based_exclusive;
                let (genomic_start_1based, genomic_end_1based) = Self::cutrun_local_span_to_genomic(
                    window,
                    local_start_1based,
                    local_end_1based,
                );
                let placement = CutRunReadPlacement {
                    orientation,
                    local_start_1based,
                    local_end_1based,
                    genomic_start_1based,
                    genomic_end_1based,
                    mismatches,
                    identity_fraction,
                    seed_matches,
                };
                let replace = best.as_ref().map_or(true, |current| {
                    placement
                        .identity_fraction
                        .partial_cmp(&current.identity_fraction)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .is_gt()
                        || (placement.identity_fraction - current.identity_fraction).abs()
                            < f64::EPSILON
                            && (placement.mismatches < current.mismatches
                                || (placement.mismatches == current.mismatches
                                    && (placement.seed_matches > current.seed_matches
                                        || (placement.seed_matches == current.seed_matches
                                            && (placement.local_start_1based
                                                < current.local_start_1based
                                                || (placement.local_start_1based
                                                    == current.local_start_1based
                                                    && placement.orientation.as_str()
                                                        < current.orientation.as_str()))))))
                });
                if replace {
                    best = Some(placement);
                }
            }
        }
        best
    }

    fn cutrun_fragment_from_row(row: &CutRunReadUnitRow) -> Option<CutRunFragmentSpan> {
        let local_start_1based = row.fragment_local_start_1based?;
        let local_end_1based = row.fragment_local_end_1based?;
        let genomic_start_1based = row.fragment_genomic_start_1based?;
        let genomic_end_1based = row.fragment_genomic_end_1based?;
        Some(CutRunFragmentSpan {
            fragment_id: format!("fragment_{}", row.unit_index + 1),
            normalized_read_id: row.normalized_read_id.clone(),
            status: row.status,
            local_start_1based,
            local_end_1based,
            genomic_start_1based,
            genomic_end_1based,
            length_bp: local_end_1based
                .saturating_sub(local_start_1based)
                .saturating_add(1),
            left_cut_site_local_1based: local_start_1based,
            right_cut_site_local_1based: local_end_1based,
            duplicate_count: row.duplicate_count.max(1),
        })
    }

    fn summarize_cutrun_support_clusters(
        coverage: &[u32],
        cut_site_counts: &[u32],
        fragments: &[CutRunFragmentSpan],
        window: &CutRunReferenceWindow,
    ) -> Vec<CutRunSupportCluster> {
        let mut clusters = Vec::<CutRunSupportCluster>::new();
        let mut idx = 0usize;
        while idx < coverage.len() {
            if coverage[idx] == 0 && cut_site_counts.get(idx).copied().unwrap_or(0) == 0 {
                idx += 1;
                continue;
            }
            let start = idx;
            let mut peak_coverage = 0u32;
            let mut total_cut_sites = 0u32;
            while idx < coverage.len()
                && (coverage[idx] > 0 || cut_site_counts.get(idx).copied().unwrap_or(0) > 0)
            {
                peak_coverage = peak_coverage.max(coverage[idx]);
                total_cut_sites =
                    total_cut_sites.saturating_add(cut_site_counts.get(idx).copied().unwrap_or(0));
                idx += 1;
            }
            let end = idx.saturating_sub(1);
            let local_start_1based = start + 1;
            let local_end_1based = end + 1;
            let (genomic_start_1based, genomic_end_1based) =
                Self::cutrun_local_span_to_genomic(window, local_start_1based, local_end_1based);
            let fragment_count = fragments
                .iter()
                .filter(|fragment| {
                    fragment.local_start_1based <= local_end_1based
                        && fragment.local_end_1based >= local_start_1based
                })
                .count();
            clusters.push(CutRunSupportCluster {
                cluster_index: clusters.len(),
                local_start_1based,
                local_end_1based,
                genomic_start_1based,
                genomic_end_1based,
                peak_coverage,
                total_cut_sites,
                fragment_count,
            });
        }
        clusters
    }

    fn build_cutrun_report_from_units(
        report_id: &str,
        seq_id: &str,
        input_r1_path: &str,
        input_r2_path: Option<String>,
        input_format: CutRunInputFormat,
        read_layout: CutRunReadLayout,
        roi_flank_bp: usize,
        deduplicate_fragments: bool,
        seed_filter: CutRunSeedFilterConfig,
        align_config: CutRunAlignConfig,
        window: &CutRunReferenceWindow,
        units: &[CutRunReadUnitRow],
        mean_read_length_bp: f64,
        mean_mapped_read_length_bp: f64,
        base_warnings: &[String],
    ) -> CutRunReadReport {
        let mut fragments = units
            .iter()
            .filter_map(Self::cutrun_fragment_from_row)
            .collect::<Vec<_>>();
        let mut unit_rows = units.to_vec();
        if deduplicate_fragments {
            let mut counts = BTreeMap::<(usize, usize, CutRunReadUnitStatus), usize>::new();
            for fragment in &fragments {
                *counts
                    .entry((
                        fragment.local_start_1based,
                        fragment.local_end_1based,
                        fragment.status,
                    ))
                    .or_insert(0) += 1;
            }
            let mut seen = BTreeSet::<(usize, usize, CutRunReadUnitStatus)>::new();
            fragments.retain_mut(|fragment| {
                let key = (
                    fragment.local_start_1based,
                    fragment.local_end_1based,
                    fragment.status,
                );
                fragment.duplicate_count = counts.get(&key).copied().unwrap_or(1);
                seen.insert(key)
            });
            let mut first_seen = BTreeSet::<(usize, usize, CutRunReadUnitStatus)>::new();
            for row in &mut unit_rows {
                let Some(start) = row.fragment_local_start_1based else {
                    continue;
                };
                let Some(end) = row.fragment_local_end_1based else {
                    continue;
                };
                let key = (start, end, row.status);
                row.duplicate_count = counts.get(&key).copied().unwrap_or(1);
                let first = first_seen.insert(key);
                row.deduplicated = !first && row.duplicate_count > 1;
            }
        } else {
            for row in &mut unit_rows {
                row.duplicate_count = row.duplicate_count.max(1);
            }
        }
        let mut coverage = vec![0u32; window.sequence.len()];
        let mut cut_site_counts = vec![0u32; window.sequence.len()];
        for fragment in &fragments {
            let start = fragment.local_start_1based.saturating_sub(1);
            let end = fragment.local_end_1based.min(coverage.len());
            for position in start..end {
                coverage[position] = coverage[position].saturating_add(1);
            }
            if fragment.left_cut_site_local_1based > 0
                && fragment.left_cut_site_local_1based <= cut_site_counts.len()
            {
                cut_site_counts[fragment.left_cut_site_local_1based - 1] =
                    cut_site_counts[fragment.left_cut_site_local_1based - 1].saturating_add(1);
            }
            if fragment.right_cut_site_local_1based > 0
                && fragment.right_cut_site_local_1based <= cut_site_counts.len()
            {
                cut_site_counts[fragment.right_cut_site_local_1based - 1] =
                    cut_site_counts[fragment.right_cut_site_local_1based - 1].saturating_add(1);
            }
        }
        let total_units = unit_rows.len();
        let mapped_units = unit_rows
            .iter()
            .filter(|row| row.fragment_local_start_1based.is_some())
            .count();
        let concordant_pair_count = unit_rows
            .iter()
            .filter(|row| matches!(row.status, CutRunReadUnitStatus::ConcordantPair))
            .count();
        let orphan_r1_count = unit_rows
            .iter()
            .filter(|row| matches!(row.status, CutRunReadUnitStatus::OrphanR1))
            .count();
        let orphan_r2_count = unit_rows
            .iter()
            .filter(|row| matches!(row.status, CutRunReadUnitStatus::OrphanR2))
            .count();
        let unmatched_pair_count = unit_rows
            .iter()
            .filter(|row| {
                matches!(
                    row.status,
                    CutRunReadUnitStatus::R1OnlyMapped
                        | CutRunReadUnitStatus::R2OnlyMapped
                        | CutRunReadUnitStatus::DiscordantPair
                        | CutRunReadUnitStatus::PairUnmapped
                )
            })
            .count();
        let support_clusters = Self::summarize_cutrun_support_clusters(
            &coverage,
            &cut_site_counts,
            &fragments,
            window,
        );
        CutRunReadReport {
            schema: CUTRUN_READ_REPORT_SCHEMA.to_string(),
            report_id: report_id.to_string(),
            seq_id: seq_id.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            input_r1_path: input_r1_path.to_string(),
            input_r2_path,
            input_format,
            read_layout,
            roi_flank_bp,
            deduplicate_fragments,
            seed_filter,
            align_config,
            genome_id: window.genome_id.clone(),
            chromosome: window.chromosome.clone(),
            reference_window_start_1based: window.window_start_1based,
            reference_window_end_1based: window.window_end_1based,
            reference_window_orientation: window.orientation.as_str().to_string(),
            roi_local_start_1based: window.roi_local_start_1based,
            roi_local_end_1based: window.roi_local_end_1based,
            reference_window_length_bp: window.sequence.len(),
            total_units,
            mapped_units,
            fragment_count: fragments.len(),
            concordant_pair_count,
            orphan_r1_count,
            orphan_r2_count,
            unmatched_pair_count,
            mean_read_length_bp,
            mean_mapped_read_length_bp,
            coverage,
            cut_site_counts,
            units: unit_rows,
            fragments,
            support_clusters,
            warnings: base_warnings.to_vec(),
        }
    }

    fn write_cutrun_report_snapshot(
        report: &CutRunReadReport,
        path: &str,
    ) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize CUT&RUN report snapshot '{}': {e}",
                path
            ),
        })?;
        fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write CUT&RUN report snapshot '{}': {e}", path),
        })
    }

    pub fn interpret_cutrun_reads(
        &mut self,
        seq_id: &str,
        input_r1_path: &str,
        input_r2_path: Option<&str>,
        input_format: CutRunInputFormat,
        read_layout: CutRunReadLayout,
        roi_flank_bp: usize,
        seed_filter: &CutRunSeedFilterConfig,
        align_config: &CutRunAlignConfig,
        deduplicate_fragments: bool,
        report_id: Option<&str>,
        checkpoint_path: Option<&str>,
        checkpoint_every_reads: usize,
    ) -> Result<CutRunReadReport, EngineError> {
        match read_layout {
            CutRunReadLayout::SingleEnd => {
                if input_r2_path.is_some() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "CUT&RUN single-end interpretation does not accept input_r2_path"
                            .to_string(),
                    });
                }
            }
            CutRunReadLayout::PairedEnd => {
                if input_r2_path
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .is_none()
                {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "CUT&RUN paired-end interpretation requires input_r2_path"
                            .to_string(),
                    });
                }
            }
        }
        if seed_filter.kmer_len == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "CUT&RUN seed_filter.kmer_len must be >= 1".to_string(),
            });
        }
        if !(0.0..=1.0).contains(&align_config.min_identity_fraction) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "CUT&RUN min_identity_fraction must be within 0.0..=1.0".to_string(),
            });
        }
        let report_id = match report_id {
            Some(raw) => Self::normalize_cutrun_report_id(raw)?,
            None => Self::normalize_cutrun_report_id(&format!(
                "cutrun_reads_{}_{}",
                seq_id,
                Self::now_unix_ms()
            ))?,
        };
        let checkpoint_every_reads = checkpoint_every_reads.max(1);
        let mut warnings = Vec::<String>::new();
        if let Some(path) = checkpoint_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            warnings.push(format!(
                "checkpoint snapshots enabled every {} unit(s) at '{}'",
                checkpoint_every_reads, path
            ));
        }
        let window = self.build_cutrun_reference_window(seq_id, roi_flank_bp)?;
        warnings.extend(window.warnings.clone());
        let input_r1_records = Self::parse_cutrun_input_records(input_r1_path, input_format)?;
        let input_r2_records = match input_r2_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            Some(path) => Self::parse_cutrun_input_records(path, input_format)?,
            None => vec![],
        };
        let total_record_count = input_r1_records
            .len()
            .saturating_add(input_r2_records.len())
            .max(1) as f64;
        let total_read_bases = input_r1_records
            .iter()
            .chain(input_r2_records.iter())
            .map(|record| record.sequence.len())
            .sum::<usize>();
        let mean_read_length_bp = total_read_bases as f64 / total_record_count;
        let units = Self::build_cutrun_input_units(input_r1_records, input_r2_records, read_layout);
        let mut mapped_read_bases = 0usize;
        let mut mapped_read_count = 0usize;
        let mut rows = Vec::<CutRunReadUnitRow>::new();
        for unit in units {
            let r1_placement = unit.r1.as_ref().and_then(|record| {
                Self::place_cutrun_read_against_window(
                    &record.sequence,
                    &window,
                    seed_filter,
                    align_config,
                )
            });
            let r2_placement = unit.r2.as_ref().and_then(|record| {
                Self::place_cutrun_read_against_window(
                    &record.sequence,
                    &window,
                    seed_filter,
                    align_config,
                )
            });
            if let Some(record) = unit.r1.as_ref().filter(|_| r1_placement.is_some()) {
                mapped_read_bases = mapped_read_bases.saturating_add(record.sequence.len());
                mapped_read_count = mapped_read_count.saturating_add(1);
            }
            if let Some(record) = unit.r2.as_ref().filter(|_| r2_placement.is_some()) {
                mapped_read_bases = mapped_read_bases.saturating_add(record.sequence.len());
                mapped_read_count = mapped_read_count.saturating_add(1);
            }
            let mut row = CutRunReadUnitRow {
                unit_index: rows.len(),
                normalized_read_id: unit.normalized_read_id.clone(),
                status: CutRunReadUnitStatus::PairUnmapped,
                r1_record_index: unit.r1.as_ref().map(|record| record.record_index),
                r1_header_id: unit.r1.as_ref().map(|record| record.header_id.clone()),
                r1_read_length_bp: unit.r1.as_ref().map(|record| record.sequence.len()),
                r1_placement,
                r2_record_index: unit.r2.as_ref().map(|record| record.record_index),
                r2_header_id: unit.r2.as_ref().map(|record| record.header_id.clone()),
                r2_read_length_bp: unit.r2.as_ref().map(|record| record.sequence.len()),
                r2_placement,
                fragment_local_start_1based: None,
                fragment_local_end_1based: None,
                fragment_genomic_start_1based: None,
                fragment_genomic_end_1based: None,
                deduplicated: false,
                duplicate_count: 1,
                notes: vec![],
            };
            match read_layout {
                CutRunReadLayout::SingleEnd => {
                    if let Some(placement) = row.r1_placement.clone() {
                        row.status = CutRunReadUnitStatus::SingleMapped;
                        row.fragment_local_start_1based = Some(placement.local_start_1based);
                        row.fragment_local_end_1based = Some(placement.local_end_1based);
                        row.fragment_genomic_start_1based = Some(placement.genomic_start_1based);
                        row.fragment_genomic_end_1based = Some(placement.genomic_end_1based);
                    } else {
                        row.status = CutRunReadUnitStatus::SingleUnmapped;
                    }
                }
                CutRunReadLayout::PairedEnd => {
                    match (
                        row.r1_placement.clone(),
                        row.r2_placement.clone(),
                        row.r1_record_index.is_some(),
                        row.r2_record_index.is_some(),
                    ) {
                        (Some(left), Some(right), true, true) => {
                            let fragment_start =
                                left.local_start_1based.min(right.local_start_1based);
                            let fragment_end = left.local_end_1based.max(right.local_end_1based);
                            let fragment_span_bp = fragment_end
                                .saturating_sub(fragment_start)
                                .saturating_add(1);
                            let opposite_orientation = left.orientation != right.orientation;
                            if opposite_orientation
                                && fragment_span_bp <= align_config.max_fragment_span_bp
                            {
                                row.status = CutRunReadUnitStatus::ConcordantPair;
                                row.fragment_local_start_1based = Some(fragment_start);
                                row.fragment_local_end_1based = Some(fragment_end);
                                let (genomic_start_1based, genomic_end_1based) =
                                    Self::cutrun_local_span_to_genomic(
                                        &window,
                                        fragment_start,
                                        fragment_end,
                                    );
                                row.fragment_genomic_start_1based = Some(genomic_start_1based);
                                row.fragment_genomic_end_1based = Some(genomic_end_1based);
                            } else {
                                row.status = CutRunReadUnitStatus::DiscordantPair;
                                row.notes.push(format!(
                                    "pair placements were not concordant (opposite_orientation={}, fragment_span_bp={})",
                                    opposite_orientation, fragment_span_bp
                                ));
                            }
                        }
                        (Some(placement), None, true, true) => {
                            row.status = CutRunReadUnitStatus::R1OnlyMapped;
                            row.fragment_local_start_1based = Some(placement.local_start_1based);
                            row.fragment_local_end_1based = Some(placement.local_end_1based);
                            row.fragment_genomic_start_1based =
                                Some(placement.genomic_start_1based);
                            row.fragment_genomic_end_1based = Some(placement.genomic_end_1based);
                        }
                        (None, Some(placement), true, true) => {
                            row.status = CutRunReadUnitStatus::R2OnlyMapped;
                            row.fragment_local_start_1based = Some(placement.local_start_1based);
                            row.fragment_local_end_1based = Some(placement.local_end_1based);
                            row.fragment_genomic_start_1based =
                                Some(placement.genomic_start_1based);
                            row.fragment_genomic_end_1based = Some(placement.genomic_end_1based);
                        }
                        (None, None, true, true) => {
                            row.status = CutRunReadUnitStatus::PairUnmapped;
                        }
                        (Some(placement), _, true, false) => {
                            row.status = CutRunReadUnitStatus::OrphanR1;
                            row.fragment_local_start_1based = Some(placement.local_start_1based);
                            row.fragment_local_end_1based = Some(placement.local_end_1based);
                            row.fragment_genomic_start_1based =
                                Some(placement.genomic_start_1based);
                            row.fragment_genomic_end_1based = Some(placement.genomic_end_1based);
                        }
                        (None, _, true, false) => {
                            row.status = CutRunReadUnitStatus::OrphanR1;
                            row.notes
                                .push("orphan R1 did not map inside the ROI window".to_string());
                        }
                        (_, Some(placement), false, true) => {
                            row.status = CutRunReadUnitStatus::OrphanR2;
                            row.fragment_local_start_1based = Some(placement.local_start_1based);
                            row.fragment_local_end_1based = Some(placement.local_end_1based);
                            row.fragment_genomic_start_1based =
                                Some(placement.genomic_start_1based);
                            row.fragment_genomic_end_1based = Some(placement.genomic_end_1based);
                        }
                        (_, None, false, true) => {
                            row.status = CutRunReadUnitStatus::OrphanR2;
                            row.notes
                                .push("orphan R2 did not map inside the ROI window".to_string());
                        }
                        _ => {}
                    }
                }
            }
            rows.push(row);
            if let Some(path) = checkpoint_path
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                if rows.len() % checkpoint_every_reads == 0 {
                    let snapshot = Self::build_cutrun_report_from_units(
                        &report_id,
                        seq_id,
                        input_r1_path,
                        input_r2_path.map(|value| value.to_string()),
                        input_format,
                        read_layout,
                        roi_flank_bp,
                        deduplicate_fragments,
                        seed_filter.clone(),
                        align_config.clone(),
                        &window,
                        &rows,
                        mean_read_length_bp,
                        if mapped_read_count == 0 {
                            0.0
                        } else {
                            mapped_read_bases as f64 / mapped_read_count as f64
                        },
                        &warnings,
                    );
                    Self::write_cutrun_report_snapshot(&snapshot, path)?;
                }
            }
        }
        let mean_mapped_read_length_bp = if mapped_read_count == 0 {
            0.0
        } else {
            mapped_read_bases as f64 / mapped_read_count as f64
        };
        let report = Self::build_cutrun_report_from_units(
            &report_id,
            seq_id,
            input_r1_path,
            input_r2_path.map(|value| value.to_string()),
            input_format,
            read_layout,
            roi_flank_bp,
            deduplicate_fragments,
            seed_filter.clone(),
            align_config.clone(),
            &window,
            &rows,
            mean_read_length_bp,
            mean_mapped_read_length_bp,
            &warnings,
        );
        if let Some(path) = checkpoint_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            Self::write_cutrun_report_snapshot(&report, path)?;
        }
        self.upsert_cutrun_read_report(report.clone())?;
        Ok(report)
    }

    pub fn export_cutrun_read_coverage(
        &self,
        report_id: &str,
        path: &str,
        kind: CutRunCoverageKind,
    ) -> Result<CutRunReadCoverageExport, EngineError> {
        let report = self.get_cutrun_read_report(report_id)?;
        let mut lines = Vec::<String>::new();
        match kind {
            CutRunCoverageKind::Coverage => {
                lines.push("local_pos_1based\tgenomic_pos_1based\tcoverage".to_string());
                for (idx, value) in report.coverage.iter().copied().enumerate() {
                    let local_pos_1based = idx + 1;
                    let genomic_pos_1based = if report.reference_window_orientation == "reverse" {
                        report
                            .reference_window_end_1based
                            .saturating_sub(local_pos_1based.saturating_sub(1))
                    } else {
                        report
                            .reference_window_start_1based
                            .saturating_add(local_pos_1based.saturating_sub(1))
                    };
                    lines.push(format!("{local_pos_1based}\t{genomic_pos_1based}\t{value}"));
                }
            }
            CutRunCoverageKind::CutSites => {
                lines.push("local_pos_1based\tgenomic_pos_1based\tcut_sites".to_string());
                for (idx, value) in report.cut_site_counts.iter().copied().enumerate() {
                    let local_pos_1based = idx + 1;
                    let genomic_pos_1based = if report.reference_window_orientation == "reverse" {
                        report
                            .reference_window_end_1based
                            .saturating_sub(local_pos_1based.saturating_sub(1))
                    } else {
                        report
                            .reference_window_start_1based
                            .saturating_add(local_pos_1based.saturating_sub(1))
                    };
                    lines.push(format!("{local_pos_1based}\t{genomic_pos_1based}\t{value}"));
                }
            }
            CutRunCoverageKind::Fragments => {
                lines.push("fragment_id\tnormalized_read_id\tstatus\tlocal_start_1based\tlocal_end_1based\tgenomic_start_1based\tgenomic_end_1based\tlength_bp\tduplicate_count".to_string());
                for fragment in &report.fragments {
                    lines.push(format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        fragment.fragment_id,
                        fragment.normalized_read_id,
                        fragment.status.as_str(),
                        fragment.local_start_1based,
                        fragment.local_end_1based,
                        fragment.genomic_start_1based,
                        fragment.genomic_end_1based,
                        fragment.length_bp,
                        fragment.duplicate_count
                    ));
                }
            }
        }
        fs::write(path, lines.join("\n")).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write CUT&RUN coverage export to '{}': {e}", path),
        })?;
        Ok(CutRunReadCoverageExport {
            schema: CUTRUN_READ_COVERAGE_EXPORT_SCHEMA.to_string(),
            path: path.to_string(),
            report_id: report.report_id,
            kind,
            row_count: lines.len().saturating_sub(1),
        })
    }
}
